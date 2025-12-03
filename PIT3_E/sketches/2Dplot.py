import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Rectangle

# Paramètres
a = 3                      # Rayon du cylindre ou demi-côté du carré
Rc = 3.25                     # Rayon du cylindre extérieur
# Distance = a + Rc          # Distance entre les centres
Distance=2*a
Nb=2
# Coordonnées des centres (config Nb9_C)
coord = np.zeros((Nb, 2))
coord[0, :] = [0, 0]
coord[1, :] = [Distance, 0]
# coord[2, :] = [Distance, -Distance]
# coord[3, :] = [0, -Distance]
# coord[4, :] = [-Distance, -Distance]
# coord[5, :] = [-Distance, 0]
# coord[6, :] = [-Distance, Distance]
# coord[7, :] = [0, Distance]
# coord[8, :] = [Distance, Distance]

# Choix du type de corps : "cylinder" ou "square"
body_type = "cylinder"  # ou "square"

# Création de la figure
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal')

# Couleurs
body_color = "black"
outer_color = "blue"

# Tracé des corps
for (x, y) in coord:
    # Outer cylinder
    outer = Circle((x, y), Rc, edgecolor=outer_color, facecolor='none', linestyle='--', linewidth=1.5)
    ax.add_patch(outer)

    # Corps principal
    if body_type == "cylinder":
        body = Circle((x, y), a, color=body_color, alpha=0.6)
    elif body_type == "square":
        body = Rectangle((x - a, y - a), 2 * a, 2 * a, color=body_color, alpha=0.6)
    else:
        raise ValueError("Type de corps inconnu : choisir 'cylinder' ou 'square'")
    
    ax.add_patch(body)

# Ajustement de la vue
# margin = Distance + Rc + 2
margin1=2*a
margin2=Distance+2*a
ax.set_xlim(-margin1, margin2)
ax.set_ylim(-margin1, margin1)
# ax.set_title(f"Configuration Nb9_C ({body_type})")
ax.set_xlabel("X")
ax.set_ylabel("Y")
plt.grid(True)
plt.savefig(f"Config_Nb2_d2a_{body_type}")
plt.show()

