import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import sys
def read_nemoh_rectangular_dat(filename):
    with open(filename, 'r') as f:
        header = f.readline().strip().split()
        if header[0] != '2':
            raise ValueError("Not right format. Expected Nemoh meshfile .dat")
        symmetric = bool(int(header[1]))
        vertices = []
        while True:
            line = f.readline()
            if line.strip().startswith('0'):
                break
            parts = line.strip().split()
            x, y, z = map(float, parts[1:4])
            vertices.append([x, y, z])
        vertices = np.array(vertices)
        quads = []
        while True:
            line = f.readline()
            if line.strip().startswith('0'):
                break
            indices = list(map(int, line.strip().split()))
            if len(indices) != 4:
                raise ValueError("Connectivités non quadrilatérales.")
            quads.append([i - 1 for i in indices])
        quads = np.array(quads)

        return vertices, quads


def compute_normals_and_centers(vertices, quads):
    centers = []
    normals = []

    for quad in quads:
        p0, p1, p2, p3 = vertices[quad]
        is_triangle = np.allclose(p2, p3)

        if is_triangle:
            center = (p0 + p1 + p2) / 3
            centers.append(center)
            v1 = p1 - p0
            v2 = p2 - p0
            n = np.cross(v1, v2)
        else:
            center = (p0 + p1 + p2 + p3) / 4
            centers.append(center)
            n1 = np.cross(p1 - p0, p2 - p0)
            n2 = np.cross(p2 - p0, p3 - p0)
            n = n1 + n2
        norm = np.linalg.norm(n)
        if norm != 0:
            n = n / norm
        normals.append(n)
    return np.array(centers), np.array(normals)


def classify_normals(vertices, centers, normals):
    interior_point = vertices.mean(axis=0)
    classification = []
    for center, normal in zip(centers, normals):
        v = interior_point - center
        # Si dot(normal, v) < 0 → outward
        is_outward = np.dot(normal, v) < 0
        classification.append("outward" if is_outward else "inward")

    return classification


def plot_quadrilateral_mesh(vertices, quads):
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    mesh = Poly3DCollection(vertices[quads], alpha=0.7)
    mesh.set_edgecolor('k')
    mesh.set_facecolor('lightblue')
    ax.add_collection3d(mesh)
    centers, normals = compute_normals_and_centers(vertices, quads)
    classification = classify_normals(vertices, centers, normals)
    colors = ['g' if c == "outward" else 'r' for c in classification]

    ax.quiver(centers[:, 0], centers[:, 1], centers[:, 2],
              normals[:, 0], normals[:, 1], normals[:, 2],
              length=1, color=colors, normalize=True)
    green_patch = plt.Line2D([0], [0], marker='o', color='w', label='Normal outward',
                             markerfacecolor='g', markersize=10)
    red_patch = plt.Line2D([0], [0], marker='o', color='w', label='Normal inward',
                           markerfacecolor='r', markersize=10)
    ax.legend(handles=[green_patch, red_patch])

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    max_range = (vertices.max(axis=0) - vertices.min(axis=0)).max() / 2
    mid = vertices.mean(axis=0)
    ax.set_xlim(mid[0] - max_range, mid[0] + max_range)
    ax.set_ylim(mid[1] - max_range, mid[1] + max_range)
    ax.set_zlim(mid[2] - max_range, mid[2] + max_range)
    plt.tight_layout()
    plt.show()

    for i, c in enumerate(classification):
    	if c=="inward":
    	    print(f"Panel {i:4d} : normal {c}")
    return classification


if len(sys.argv) != 2:
    print("Use: python MeshPlot3D_normals.py meshfile.dat")
    sys.exit(1)

input_file = sys.argv[1]
vertices, quads = read_nemoh_rectangular_dat(input_file)
classification = plot_quadrilateral_mesh(vertices, quads)

