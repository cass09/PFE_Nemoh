import numpy as np

def write_nemoh_parallelepiped_mesh(Lx=3, Ly=2, Lz=1, ndiv=3, nom="barge", SYM=False):
    """
    Génère un pavé rectangulaire (Lx × Ly × Lz) maillé en rectangles sur chaque face,
    et l'écrit au format mesh.dat de Nemoh.
    
    Lx, Ly, Lz : dimensions du pavé en mètres
    ndiv       : nombre de divisions par arête
    SYM        : si True, ne conserve que la moitié où x >= 0
    """
    nodes = []
    node_coords = {}
    step_x = Lx / ndiv
    step_y = Ly / ndiv
    step_z = -Lz / ndiv
    x_offset = -Lx / 2
    y_offset = -Ly / 2

    def add_face_nodes(fixed_val, step1, step2, coord_order):
        face_nodes = []
        for i in range(ndiv + 1):
            for j in range(ndiv + 1):
                coords = [0, 0, 0]
                coords[coord_order[0]] = i * step1
                coords[coord_order[1]] = j * step2
                coords[coord_order[2]] = fixed_val
                coords[0] += x_offset
                coords[1] += y_offset
                face_nodes.append(tuple(coords))
        return face_nodes

    faces = []
    faces.append(add_face_nodes(0, step_x, step_y, (0, 1, 2)))       # z = 0
    faces.append(add_face_nodes(-Lz, step_x, step_y, (0, 1, 2)))      # z = -Lz
    faces.append(add_face_nodes(0, step_y, step_z, (1, 2, 0)))       # x = 0
    faces.append(add_face_nodes(Lx, step_y, step_z, (1, 2, 0)))      # x = Lx
    faces.append(add_face_nodes(0, step_x, step_z, (0, 2, 1)))       # y = 0
    faces.append(add_face_nodes(Ly, step_x, step_z, (0, 2, 1)))      # y = Ly

    all_nodes = sum(faces, [])
    node_map = {}
    filtered_nodes = []
    for i, (x, y, z) in enumerate(all_nodes):
        if not SYM or y >= 0:
            idx = len(filtered_nodes) + 1  # Numérotation 1-based
            node_map[i] = idx
            filtered_nodes.append((x, y, z))

    elements = []
    offset = 0
    face_size = (ndiv + 1) ** 2

    for face in faces:
        for i in range(ndiv):
            for j in range(ndiv):
                n1 = offset + i * (ndiv + 1) + j
                n2 = offset + i * (ndiv + 1) + (j + 1)
                n3 = offset + (i + 1) * (ndiv + 1) + (j + 1)
                n4 = offset + (i + 1) * (ndiv + 1) + j
                nodes_idx = [n1, n2, n3, n4]
                if not SYM or all(all_nodes[k][1] >= 0 for k in nodes_idx):
                    try:
                        elements.append(tuple(node_map[k] for k in nodes_idx))
                    except KeyError:
                        continue  # Si un noeud est absent car x < 0
        offset += face_size

    n_nodes = len(filtered_nodes)
    n_elements = len(elements)
    suffix = "SYM" if SYM else ""
    filename = f"{nom}X{Lx}Y{Ly}H{Lz}{suffix}"

    with open(filename, "w") as f:
        f.write(f"{n_nodes}\n{n_elements}\n")
        for x, y, z in filtered_nodes:
            f.write(f"{x:.6f} {y:.6f} {z:.6f}\n")
        for el in elements:
            f.write(" ".join(map(str, el)) + "\n")

    print(f"✅ '{filename}' créé avec {n_nodes} noeuds et {n_elements} éléments.")

# 🔧 Exemple d'utilisation :
write_nemoh_parallelepiped_mesh(Lx=6, Ly=6, Lz=6, ndiv=6, nom="barge2", SYM=True)

