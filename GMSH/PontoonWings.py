import gmsh
import sys
import numpy as np

# =======================
# Dimensions 
# =======================
L = 200      # Y
Dd = 7       # Z
Bd=30	     # X bottom lenght
Br=2.5	     # wings lenght
Dt=0.5	     # wings height
Bu = Bd-2*Br      # X top lenght

num="3"
deltaX=3.0
# =======================
# Initialisation Gmsh
# =======================
gmsh.initialize(sys.argv)
gmsh.option.setNumber("General.Terminal", 1)
model = gmsh.model
occ = model.occ
model.add("PontoonWings")

# =======================
# Points 
# =======================
# top body 
p000 = occ.addPoint(-Bu/2, 0,  0)
p100 = occ.addPoint(-Bu/2,  L/2,  0)
p010 = occ.addPoint(Bu/2,  0,  0)
p110 = occ.addPoint(Bu/2,   L/2,  0)

# bottom body
p001 = occ.addPoint(-Bd/2,  0, -Dd)
p101 = occ.addPoint(-Bd/2,   L/2, -Dd)
p011 = occ.addPoint(Bd/2,   0, -Dd)
p111 = occ.addPoint(Bd/2,    L/2, -Dd)

# middle
p002 = occ.addPoint(-Bd/2,   0,  -(Dd-Dt))
p102 = occ.addPoint(-Bd/2,    L/2,  -(Dd-Dt))
p012 = occ.addPoint(Bd/2,    0,  -(Dd-Dt))
p112 = occ.addPoint(Bd/2,     L/2,  -(Dd-Dt))

# quare corners
p003 = occ.addPoint(-Bu/2,   0,  -(Dd-Dt))
p103 = occ.addPoint(-Bu/2,    L/2,  -(Dd-Dt))
p013 = occ.addPoint(Bu/2,    0,  -(Dd-Dt))
p113 = occ.addPoint(Bu/2,     L/2,  -(Dd-Dt))

# =======================
# Arêtes
# =======================
##dessus
l1 = occ.addLine(p000, p100)
l2 = occ.addLine(p100, p110)
l3 = occ.addLine(p110, p010)
l4 = occ.addLine(p010, p000)

##dessous
l5 = occ.addLine(p001, p101)
l6 = occ.addLine(p101, p111)
l7 = occ.addLine(p111, p011)
l8 = occ.addLine(p011, p001)

## aretes verticales
l9  = occ.addLine(p002, p001)
l10 = occ.addLine(p102, p101)
l11 = occ.addLine(p112, p111)
l12 = occ.addLine(p012, p011)

## aretes horizontales
l13 = occ.addLine(p102, p002)
l14 = occ.addLine(p112, p012)
l15 = occ.addLine(p103, p003)
l16 = occ.addLine(p113, p013)

## aretes verticales
l17  = occ.addLine(p000, p003)
l18 = occ.addLine(p100, p103)
l19 = occ.addLine(p110, p113)
l20 = occ.addLine(p010, p013)

## aretes horizontales
l21 = occ.addLine(p002, p003)
l22 = occ.addLine(p102, p103)
l23 = occ.addLine(p113, p112)
l24 = occ.addLine(p013, p012)

# =======================
# Faces 
# =======================
# top
# loop_top = occ.addCurveLoop([l1, l2, -l3, -l4])
# f_top = occ.addPlaneSurface([loop_top])

# bottom
loop_bottom = occ.addCurveLoop([l5, l6, -l7, -l8])
f_bottom = occ.addPlaneSurface([loop_bottom])

# side y<0
# loop_side1 = occ.addCurveLoop([l4, arc3, l12, -l8, -l9, -arc1])
# f_side1 = occ.addPlaneSurface([loop_side1])

# side y>0
loop_side2 = occ.addCurveLoop([l18, l22, l10, l6, -l11, -l23, -l19, -l2])
f_side2 = occ.addPlaneSurface([loop_side2])

# side x<0
loop_side3 = occ.addCurveLoop([l13, l10, -l5, -l9])
f_side3 = occ.addPlaneSurface([loop_side3])
loop_starev1 = occ.addCurveLoop([l1, l18, -l15, -l17])
f_starev1 = occ.addPlaneSurface([loop_starev1])
loop_stareh1 = occ.addCurveLoop([l15, l22, -l13, -l21])
f_stareh1 = occ.addPlaneSurface([loop_stareh1])

# side x>0
loop_side4 = occ.addCurveLoop([l12, l7, -l11, -l14])
f_side4 = occ.addPlaneSurface([loop_side4])
loop_starev2 = occ.addCurveLoop([l3, l19, -l16, -l20])
f_starev2 = occ.addPlaneSurface([loop_starev2])
loop_stareh2 = occ.addCurveLoop([l16, l23, -l14, -l24])
f_stareh2 = occ.addPlaneSurface([loop_stareh2])

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
gmsh.write(f"PontoonWingsSYMm{num}.msh")
gmsh.fltk.run()
gmsh.finalize()

