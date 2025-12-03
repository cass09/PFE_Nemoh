import gmsh
import sys

# =======================
# Dimensions 
# =======================
L = 25      # X
B = 200      # Y
D = 7       # Z

deltaX=3.0
num="3"

# =======================
# Initialisation Gmsh
# =======================
gmsh.initialize(sys.argv)
gmsh.option.setNumber("General.Terminal", 1)
model = gmsh.model
occ = model.occ
model.add("Pontoon")

# =======================
# Points 
# =======================
p000 = occ.addPoint(-L/2,   0,   0)
p100 = occ.addPoint(L/2,   0,   0)
p010 = occ.addPoint(-L/2,   B/2,   0)
p110 = occ.addPoint(L/2,   B/2,   0)
p001 = occ.addPoint(-L/2,   0,   -D)
p101 = occ.addPoint(L/2,   0,   -D)
p011 = occ.addPoint(-L/2,   B/2,   -D)
p111 = occ.addPoint(L/2,   B/2,   -D)

# =======================
# Arêtes
# =======================
l1 = occ.addLine(p000, p100)
l2 = occ.addLine(p100, p110)
l3 = occ.addLine(p110, p010)
l4 = occ.addLine(p010, p000)

l5 = occ.addLine(p001, p101)
l6 = occ.addLine(p101, p111)
l7 = occ.addLine(p111, p011)
l8 = occ.addLine(p011, p001)

l9  = occ.addLine(p000, p001)
l10 = occ.addLine(p100, p101)
l11 = occ.addLine(p110, p111)
l12 = occ.addLine(p010, p011)

# =======================
# Faces 
# =======================
# top
# loop_bottom = occ.addCurveLoop([l1, l2, -l3, -l4])
# f_bottom = occ.addPlaneSurface([loop_bottom])

# bottom
loop_top = occ.addCurveLoop([l5, l6, -l7, -l8])
f_top = occ.addPlaneSurface([loop_top])

# side x=-L/2
loop_side1 = occ.addCurveLoop([l4, l12, -l8, -l9])
f_side1 = occ.addPlaneSurface([loop_side1])

# side x=L/2
loop_side2 = occ.addCurveLoop([l10, l6, -l11, -l2])
f_side2 = occ.addPlaneSurface([loop_side2])

# side y=0
# loop_side3 = occ.addCurveLoop([l1, l10, -l5, -l9])
# f_side3 = occ.addPlaneSurface([loop_side3])

# side y=B/2
loop_side4 = occ.addCurveLoop([l12, l7, -l11, -l3])
f_side4 = occ.addPlaneSurface([loop_side4])

# =======================
# Mesh
# =======================
occ.synchronize()
all_pts = gmsh.model.getEntities(0)
gmsh.model.mesh.setSize(all_pts, deltaX)    

model.mesh.generate(2)  # 2D mesh 

# =======================
# Save
# =======================
gmsh.write(f"PontoonSYMm{num}.msh")
gmsh.fltk.run()
gmsh.finalize()

