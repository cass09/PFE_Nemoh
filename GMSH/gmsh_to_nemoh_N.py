#!/usr/bin/env python3
"""
To convert a meshfile Gmsh (.msh) into a meshfile Nemoh (.dat)
"""
import sys
import gmsh
import numpy as np
import os

def compute_panel_normal(p0, p1, p2, p3):
    if np.allclose(p2, p3):  # quad degenerated for NEMOH
        return np.cross(p1 - p0, p2 - p0)
    else:  # quad
        return np.cross(p2 - p0, p3 - p1)

def force_outward_orientation(nodes, panels):
    global_center = nodes.mean(axis=0)
    corrected_panels = []

    for p in panels:
        i1, i2, i3, i4 = [k - 1 for k in p]  
        p0, p1, p2, p3 = nodes[i1], nodes[i2], nodes[i3], nodes[i4]
        n = compute_panel_normal(p0, p1, p2, p3)

        if np.allclose(p2, p3):  # triangle
            c = (p0 + p1 + p2) / 3
        else:  # quad
            c = (p0 + p1 + p2 + p3) / 4
        inward = global_center - c

        # if inward normal → inverse panel
        if np.dot(n, inward) > 0:
            corrected_panels.append([p[0], p[3], p[2], p[1]])  # flip
        else:
            corrected_panels.append(p)
    return corrected_panels

# --------------------------------------------------------------------
#  main program
# --------------------------------------------------------------------

if len(sys.argv) != 2:
    print("Use: python gmsh_to_nemoh_N.py meshfile.msh")
    sys.exit(1)

input_file = sys.argv[1]
base, _ = os.path.splitext(input_file)
output_file = base + ".dat"

gmsh.initialize()
gmsh.open(input_file)
model = gmsh.model

# ---- NODES ----
node_tags, node_coords, _ = model.mesh.getNodes()
nodes = node_coords.reshape(-1, 3)
tag_to_idx = {tag: i+1 for i, tag in enumerate(node_tags)}  # indices Nemoh = 1..N

# ---- ELEMENTS ----
elem_types, elem_tags, elem_nodes = model.mesh.getElements()

triangles = []
quads = []

for etype, tags, nodes_list in zip(elem_types, elem_tags, elem_nodes):
    if etype == 2:  # triangle
        tri = np.array(nodes_list).reshape(-1, 3)
        triangles.append(tri)
    elif etype == 3:  # quad
        quad = np.array(nodes_list).reshape(-1, 4)
        quads.append(quad)

all_panels = []

# Triangles to quads degenerated (i3 = i4)
for tri in triangles:
    for t in tri:
        i1, i2, i3 = t
        all_panels.append([
            tag_to_idx[i1],
            tag_to_idx[i2],
            tag_to_idx[i3],
            tag_to_idx[i3]
        ])

# quads
for quad in quads:
    for q in quad:
        all_panels.append([tag_to_idx[q[0]], tag_to_idx[q[1]],
                           tag_to_idx[q[2]], tag_to_idx[q[3]]])

# ---- Outwards normals ----
all_panels = force_outward_orientation(nodes, all_panels)

# ---- NEMOH mesh file ----
with open(output_file, "w") as f:
    f.write("2 1\n")  # 1 if SYM

    # points
    for idx, (x, y, z) in enumerate(nodes, start=1):
        f.write(f"{idx} {x:.6f} {y:.6f} {z:.6f}\n")

    f.write("0 0. 0. 0.\n")

    # panels
    for p in all_panels:
        f.write(" ".join(map(str, p)) + "\n")

    f.write("0 0 0 0\n")

print(f"Nemoh Mesh with OUTWARDS normals : {output_file}")

gmsh.finalize()

