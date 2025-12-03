import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def read_nemoh_rectangular_dat(filename):
    with open(filename, 'r') as f:
        header = f.readline().strip().split()
        if header[0] != '2':
            raise ValueError("Ce fichier ne semble pas être au format géométrie .dat de Nemoh")
        symmetric = bool(int(header[1]))

        # Lecture des nœuds
        vertices = []
        while True:
            line = f.readline()
            if line.strip().startswith('0'):
                break
            parts = line.strip().split()
            x, y, z = map(float, parts[1:4])
            vertices.append([x, y, z])
        vertices = np.array(vertices)

        # Lecture des rectangles
        quads = []
        while True:
            line = f.readline()
            if line.strip().startswith('0'):
                break
            indices = list(map(int, line.strip().split()))
            if len(indices) != 4:
                raise ValueError("Ce fichier contient des connectivités qui ne sont pas des quadrilatères.")
            quads.append([i - 1 for i in indices])  # passage en indices 0-based
        quads = np.array(quads)

        # Traitement de la symétrie
        # if symmetric:
        #     mirrored_vertices = vertices.copy()
        #     mirrored_vertices[:, 1] *= -1

        #     all_vertices = np.vstack((vertices, mirrored_vertices))
        #     offset = len(vertices)

        #     # Inverser l’ordre des sommets pour garder les normales correctes
        #     mirrored_quads = quads[:, [0, 3, 2, 1]] + offset
        #     all_quads = np.vstack((quads, mirrored_quads))
        # else:
        all_vertices = vertices
        all_quads = quads

    return all_vertices, all_quads

def compute_normals_and_centers(vertices, quads):
    centers = []
    normals = []

    for quad in quads:
        p0, p1, p2, p3 = vertices[quad]

        # Centre (barycentre)
        center = (p0 + p1 + p2 + p3) / 4
        centers.append(center)

        # Normale : on prend deux diagonales
        v1 = p2 - p0
        v2 = p3 - p1
        n = np.cross(v1, v2)
        norm = np.linalg.norm(n)
        if norm != 0:
            n = n / norm  # normalisation
        normals.append(n)

    return np.array(centers), np.array(normals)
    
def plot_quadrilateral_mesh(vertices, quads):
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    mesh = Poly3DCollection(vertices[quads], alpha=0.7)
    mesh.set_edgecolor('k')
    mesh.set_facecolor('lightblue')
    ax.add_collection3d(mesh)

    # ✅ Calcul et affichage des normales
    centers, normals = compute_normals_and_centers(vertices, quads)
    ax.quiver(centers[:, 0], centers[:, 1], centers[:, 2],
              normals[:, 0], normals[:, 1], normals[:, 2],
              length=0.2, color='r', normalize=True)

    # Axes
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Échelle auto
    max_range = (vertices.max(axis=0) - vertices.min(axis=0)).max() / 2
    mid = vertices.mean(axis=0)
    ax.set_xlim(mid[0] - max_range, mid[0] + max_range)
    ax.set_ylim(mid[1] - max_range, mid[1] + max_range)
    ax.set_zlim(mid[2] - max_range, mid[2] + max_range)

    plt.tight_layout()
    plt.show()

# Exemple d'utilisation :
vertices, quads = read_nemoh_rectangular_dat("PontoonSYMm1.dat")
plot_quadrilateral_mesh(vertices, quads)
