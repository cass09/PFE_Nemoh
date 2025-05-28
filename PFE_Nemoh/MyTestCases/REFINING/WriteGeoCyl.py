import numpy as np

def write_nemoh_cylinder_mesh(radius, height, n_axial=4, n_circ=8, n_radial=3, nom="Cylinder", SYM=False):
    """
    Génère un mesh.dat compatible Nemoh pour un cylindre avec fond (z = -height) maillé.
    SYM : si True, ne génère que la moitié du cylindre (x >= 0)
    n_radial : nombre de divisions radiales pour la base
    """
    nodes = []
    
    # === PARAMÈTRES ANGULAIRES ===
    max_angle = np.pi if SYM else 2 * np.pi
    step_theta = max_angle / n_circ
    n_theta = n_circ + 1  # pour que la dernière ligne latérale ferme la maille

    # === NOEUDS SURFACE LATÉRALE ===
    for i in range(n_axial + 1):
        z = -i * height / n_axial
        for j in range(n_theta):
            theta = j * step_theta
            x = radius * np.cos(theta)
            y = radius * np.sin(theta)
            nodes.append((x, y, z))
    
    lateral_node_count = len(nodes)

    # === NOEUDS FACE INFÉRIEURE (z = -height) ===
    for i in range(n_radial + 1):  # n_radial + 1 cercles concentriques
        r = radius * i / n_radial
        for j in range(n_theta):
            theta = j * step_theta
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            z = -height
            nodes.append((x, y, z))

    bottom_node_offset = lateral_node_count
    total_node_count = len(nodes)

    elements = []

    # === PANNEAUX LATÉRAUX ===
    for i in range(n_axial):
        for j in range(n_circ):
            n1 = i * n_theta + j + 1
            n2 = n1 + 1
            n3 = n2 + n_theta
            n4 = n1 + n_theta
            elements.append((n1, n2, n3, n4))

    # === PANNEAUX FACE INFÉRIEURE ===
    for i in range(n_radial):
        for j in range(n_circ):
            # indices des 4 coins du rectangle entre deux cercles
            r1 = bottom_node_offset + i * n_theta
            r2 = bottom_node_offset + (i + 1) * n_theta

            n1 = r1 + j + 1
            n2 = r1 + (j + 1) % n_theta + 1
            n3 = r2 + (j + 1) % n_theta + 1
            n4 = r2 + j + 1

            elements.append((n1, n2, n3, n4))

    n_nodes = len(nodes)
    n_elements = len(elements)
    suffix = "SYM" if SYM else ""
    filename = f"{nom}R{radius}H{height}{suffix}"

    # === ÉCRITURE DU FICHIER ===
    with open(filename, "w") as f:
        f.write(f"{n_nodes}\n{n_elements}\n")
        for x, y, z in nodes:
            f.write(f"{x:.6f} {y:.6f} {z:.6f}\n")
        for el in elements:
            f.write(" ".join(map(str, el)) + "\n")

    print(f"✅ Fichier '{filename}' créé avec {n_nodes} nœuds et {n_elements} éléments.")

# 🔧 Exemple d'utilisation
write_nemoh_cylinder_mesh(radius=3, height=6, n_axial=6, n_circ=15, n_radial=3, nom="Cylinder", SYM=True)
