import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

# =======================
# Paramètres schématiques
# =======================
a = 1.0        # Taille du corps (rayon si cylindre, demi-côté si carré) en UNITS arbitraires
Rc = 1.2       # Rayon du cylindre extérieur (par rapport à a, arbitraire aussi)
d = 5.0        # Distance arbitraire entre les corps (pas à l’échelle)

Nb = 2         # Nombre de corps
body_type = "cylinder"  # "cylinder" ou "square"

# =======================
# Coordonnées des centres
# =======================
coord = []
coord.append((0, 0))     # premier corps à l'origine
coord.append((d, 0))     # deuxième corps à une distance arbitraire d sur X

# =======================
# Création de la figure
# =======================
fig, ax = plt.subplots(figsize=(8, 4))
ax.set_aspect("equal")

body_color = "black"
outer_color = "blue"

# =======================
# Tracé des corps
# =======================
for (x, y) in coord:
    # Outer "cylinder"
    outer = Circle((x, y), Rc, edgecolor=outer_color, facecolor='none',
                   linestyle='--', linewidth=1.5)
    ax.add_patch(outer)

    # Corps principal
    if body_type == "cylinder":
        body = Circle((x, y), a, color=body_color, alpha=0.6)
    elif body_type == "square":
        body = Rectangle((x - a, y - a), 2 * a, 2 * a, color=body_color, alpha=0.6)
    else:
        raise ValueError("Type de corps inconnu : choisir 'cylinder' ou 'square'")
    
    ax.add_patch(body)

# =======================
# Ajustement de la vue
# =======================
ax.set_xlim(-2*a, d + 2*a)
ax.set_ylim(-2*a, 2*a)
# Supprimer les axes et les graduations
ax.axis("off")  # supprime les axes, graduations et labels

ax.set_xlabel("X ")
ax.set_ylabel("Y")
# ax.set_title(f"Configuration avec distance d = {d} (schématique, non à l'échelle)")
ax.grid(True)

plt.savefig(f"Config_schematic_{body_type}.png", dpi=150)
plt.show()

