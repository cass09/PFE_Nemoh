import numpy as np

def write_nemoh_cylinder_mesh(GEO, radius, height, n_axial=4, n_circ=8, n_radial=3, nom="Cylinder", SYM=False, lid=False):
    """
    Génère un mesh.dat compatible Nemoh pour un cylindre avec fond (z = -height) et éventuellement un couvercle (z = 0).
    SYM : si True, ne génère que la moitié du cylindre (y >= 0)
    lid : si True, ajoute un couvercle à z = 0
    """
    nodes = []
    
    # === PARAMÈTRES ANGULAIRES ===
    max_angle = np.pi if SYM else 2 * np.pi
    step_theta = max_angle / n_circ
    n_theta = n_circ + 1  # pour fermer la maille

    # === NOEUDS SURFACE LATÉRALE ===
    for i in range(n_axial + 1):
        z = -i * height / n_axial
        for j in range(n_theta):
            theta = j * step_theta
            x = radius * np.cos(theta)
            y = radius * np.sin(theta)
            nodes.append((x, y, z))

    lateral_node_count = len(nodes)

    # === NOEUDS BASE INFÉRIEURE ===
    for i in range(n_radial + 1):
        r = radius * i / n_radial
        for j in range(n_theta):
            theta = j * step_theta
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            z = -height
            nodes.append((x, y, z))

    bottom_node_offset = lateral_node_count
    bottom_node_count = (n_radial + 1) * n_theta

    # === NOEUDS LID SUPÉRIEUR (OPTIONNEL) ===
    if lid:
        for i in range(n_radial + 1):
            r = radius * i / n_radial
            for j in range(n_theta):
                theta = j * step_theta
                x = r * np.cos(theta)
                y = r * np.sin(theta)
                z = 0
                nodes.append((x, y, z))
        lid_node_offset = bottom_node_offset + bottom_node_count

    elements = []

    # === PANNEAUX LATÉRAUX ===
    for i in range(n_axial):
        for j in range(n_circ):
            n1 = i * n_theta + j + 1
            n2 = n1 + 1
            n3 = n2 + n_theta
            n4 = n1 + n_theta
            elements.append((n4, n3, n2, n1))  # Inversion pour normales

    # === PANNEAUX BASE INFÉRIEURE ===
    for i in range(n_radial):
        for j in range(n_circ):
            r1 = bottom_node_offset + i * n_theta
            r2 = bottom_node_offset + (i + 1) * n_theta

            n1 = r1 + j + 1
            n2 = r1 + (j + 1) % n_theta + 1
            n3 = r2 + (j + 1) % n_theta + 1
            n4 = r2 + j + 1

            elements.append((n1, n2, n3, n4))

    # === PANNEAUX LID SUPÉRIEUR (OPTIONNEL) ===
    if lid:
        for i in range(n_radial):
            for j in range(n_circ):
                r1 = lid_node_offset + i * n_theta
                r2 = lid_node_offset + (i + 1) * n_theta

                n1 = r1 + j + 1
                n2 = r1 + (j + 1) % n_theta + 1
                n3 = r2 + (j + 1) % n_theta + 1
                n4 = r2 + j + 1

                elements.append((n4, n3, n2, n1))  # Inversion de l'ordre => normale vers l'extérieur

    n_nodes = len(nodes)
    n_elements = len(elements)
    suffix = "SYM" if SYM else ""
    suffix += "_lid" if lid else ""

    # === ÉCRITURE DU FICHIER ===
    if GEO : 
        filename = f"{nom}R{radius}H{height}{suffix}"
        with open(filename, "w") as f:
            f.write(f"{n_nodes}\n{n_elements}\n")
            for x, y, z in nodes:
                f.write(f"{x:.6f} {y:.6f} {z:.6f}\n")
            for el in elements:
                f.write(" ".join(map(str, el)) + "\n")
    else :
        filename = f"{nom}R{radius}H{height}{suffix}.dat"
        i_sym = 1 if SYM else 0
        with open(filename, "w") as f:
            f.write(f"2     {i_sym}\n")
            ID = 1
            for x, y, z in nodes:
                f.write(f"{ID} {x:.6f} {y:.6f} {z:.6f}\n")
                ID += 1
            f.write(f"0     0      0      0\n")
            for el in elements:
                f.write(" ".join(map(str, el)) + "\n")
            f.write(f"0     0      0      0\n")

        with open(f"{nom}R{radius}H{height}{suffix}_DATA.txt", "a") as f:
            f.write(f"{filename}     Nr={n_radial}     Ntheta={n_circ}   Nz={n_axial}\n")
            f.write(f" Nnoeuds={n_nodes}    Npanels={n_elements}\n\n")

    print(f"✅ Fichier '{filename}' créé avec {n_nodes} nœuds et {n_elements} éléments.")

# 🔧 Exemple d'utilisation avec lid :
write_nemoh_cylinder_mesh(GEO=False, radius=3, height=3, n_axial=14, n_circ=38, n_radial=10, nom="MNcyl", SYM=True, lid=False)
