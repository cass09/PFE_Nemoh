import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def set_axes_equal(ax):
    '''Fixe les échelles des axes X, Y, Z pour qu'elles soient égales.'''
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

def cylindre_maillage(R=1.0, Ntheta=30, Nz=10, Depth=5.0, show_body='cylinder'):
    """
    show_body : 'cylinder', 'cube' ou None
    """

    # Maillage du cylindre extérieur
    theta = np.linspace(0, 2*np.pi, Ntheta, endpoint=False)
    i_vals = np.arange(Nz)
    z = -Depth * (1 - np.cos((np.pi / 2) * i_vals / (Nz - 1))) if Nz > 1 else np.array([0.0])

    X, Y, Z = [], [], []
    for zi in z:
        for t in theta:
            x = R * np.cos(t)
            y = R * np.sin(t)
            X.append(x)
            Y.append(y)
            Z.append(zi)

    # Affichage 3D
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X, Y, Z, color='blue', s=5, label="Outer Cylinder")

    # ----------------------
    # Ajout du body au centre
    # ----------------------
    if show_body == 'cylinder':
        # Cylindre immergé : rayon 3 m, hauteur 6 m (draft)
        Rb = 3.0
        Hb = 6.0
        z_body = np.linspace(-Hb, 0, 20)
        theta_body = np.linspace(0, 2*np.pi, 30)
        Zb, Tb = np.meshgrid(z_body, theta_body)
        Xb = Rb * np.cos(Tb)
        Yb = Rb * np.sin(Tb)

        ax.plot_surface(Xb, Yb, Zb, color='red', alpha=0.6, label="Body (cylinder)")

    elif show_body == 'barge':
        # Cube immergé : côté 6 m, centré sur z = -3 m (pour être entre z=0 et z=-6)
        L = 6.0
        cube_origin = -L / 2
        # Sommets du cube
        x = [cube_origin, cube_origin + L]
        y = [cube_origin, cube_origin + L]
        z = [-L, 0.0]

        # Faces
        faces = [
            # Bottom
            [(x[0], y[0], z[0]), (x[1], y[0], z[0]), (x[1], y[1], z[0]), (x[0], y[1], z[0])],
            # Top
            [(x[0], y[0], z[1]), (x[1], y[0], z[1]), (x[1], y[1], z[1]), (x[0], y[1], z[1])],
            # Sides
            [(x[0], y[0], z[0]), (x[0], y[0], z[1]), (x[1], y[0], z[1]), (x[1], y[0], z[0])],
            [(x[1], y[0], z[0]), (x[1], y[0], z[1]), (x[1], y[1], z[1]), (x[1], y[1], z[0])],
            [(x[1], y[1], z[0]), (x[1], y[1], z[1]), (x[0], y[1], z[1]), (x[0], y[1], z[0])],
            [(x[0], y[1], z[0]), (x[0], y[1], z[1]), (x[0], y[0], z[1]), (x[0], y[0], z[0])]
        ]

        cube = Poly3DCollection(faces, facecolors='red', linewidths=0.5, edgecolors='black', alpha=0.6)
        ax.add_collection3d(cube)

    # ----------------------
    # Style
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Outer Cylinder with immersed body : {show_body}')
    ax.set_box_aspect([1, 1, Depth / R])
    ax.set_xlim(-R, R)
    ax.set_ylim(-R, R)
    ax.set_zlim(-Depth, 0)
    # set_axes_equal(ax)  # <- Ajout ici pour échelle 1:1:1
    # --------------------
    # Forcer le cube de bounding box pour échelle 1:1:1 réelle
    # --------------------
    max_range = max(R, Depth, 6.0)  # max entre cylindre et body

    # Points d’un cube invisible pour forcer le moteur de rendu
    for x in [-max_range, max_range]:
        for y in [-max_range, max_range]:
            for z in [-max_range, 0]:
                ax.scatter([x], [y], [z], alpha=0)  # point invisible

    ax.set_box_aspect([1, 1, 2])  # maintenant ça prend bien !
    ax.view_init(elev=10, azim=10)  # vue du dessus
    plt.tight_layout()
    plt.savefig(f"OuterCylinder_withBody_{show_body}.pdf")
    plt.show()

# Exemples :
cylindre_maillage(R=3.5, Ntheta=25, Nz=12, Depth=12.0, show_body='cylinder')  # pour le cylindre immergé
# cylindre_maillage(R=4.48, Ntheta=32, Nz=12, Depth=12.0, show_body='cube')     # pour le cube immergé

