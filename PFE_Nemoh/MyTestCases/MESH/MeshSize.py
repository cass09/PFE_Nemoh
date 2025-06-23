import numpy as np

def read_meshfile(filename):
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]  # Supprime les lignes vides

    i = 0
    if not lines[i].startswith('2'):
        raise ValueError("Le fichier ne commence pas par '2 ...' comme attendu.")

    i += 1  # Skip la ligne '2 1' ou '2 0'

    # Lire les nœuds
    nodes = []
    while i < len(lines):
        parts = lines[i].split()
        if parts[0] == '0':
            break
        if len(parts) < 4:
            raise ValueError(f"Ligne invalide dans la section des nœuds : {lines[i]}")
        nodes.append([float(x) for x in parts[1:4]])
        i += 1
    nodes = np.array(nodes)
    i += 1  # ligne '0. 0. 0. 0.'

    # Lire les connectivités
    panels = []
    while i < len(lines):
        parts = lines[i].split()
        if parts[0] == '0':
            break
        if len(parts) < 3:
            raise ValueError(f"Ligne invalide dans la section des connectivités : {lines[i]}")
        panels.append([int(x) - 1 for x in parts])
        i += 1
    panels = np.array(panels)

    return nodes, panels


def max_mesh_size(nodes, panels):
    max_len = 0
    for panel in panels:
        for i in range(len(panel)):
            for j in range(i + 1, len(panel)):
                p1 = nodes[panel[i]]
                p2 = nodes[panel[j]]
                dist = np.linalg.norm(p1 - p2)
                max_len = max(max_len, dist)
    return max_len

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python script.py meshfile.dat")
        sys.exit(1)

    meshfile = sys.argv[1]
    nodes, panels = read_meshfile(meshfile)
    print(len(nodes), len(panels))
    max_size = max_mesh_size(nodes, panels)
    print(f"Max mesh size: {max_size:.5f} m")
