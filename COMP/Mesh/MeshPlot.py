import numpy as np
import matplotlib.pyplot as plt

param_d=1
Nb=2
def load_and_plot_mesh(filename, color, label):
    coords = np.loadtxt(filename, skiprows=1, max_rows=264)
    x, y = coords[:, 1], coords[:, 2]
    plt.scatter(x, y, s=20, c=color, label=label)
    # Pour relier les points dans l'ordre, si tu veux
    # plt.plot(x, y, c=color, alpha=0.5)

plt.figure(figsize=(8,8))

# Charger et afficher premier maillage
load_and_plot_mesh("barge.dat", color='blue', label='Body 1')

# Charger et afficher deuxième maillage
load_and_plot_mesh(f"barge_{param_d}.dat", color='red', label='Body 2')

theta = np.linspace(0, 2 * np.pi, 200)
Rcyl=6.36
x_cyl = Rcyl * np.cos(theta)
y_cyl = Rcyl * np.sin(theta)
plt.plot(x_cyl, y_cyl, color="black", linestyle='--', linewidth=1, label=f"Cylinder")
x_cyl = (2*Rcyl+param_d) + Rcyl * np.cos(theta)
y_cyl = Rcyl * np.sin(theta)
plt.plot(x_cyl, y_cyl, color="black", linestyle='--', linewidth=1)
            

plt.xlabel('x')
plt.ylabel('y')
plt.title(f'Barge - Nb={Nb} - d={param_d}m')
plt.axis('equal')
plt.grid(True)
plt.legend()
# plt.show()
plt.savefig(f"Mesh_Nb{Nb}_X_d{param_d}.png")
plt.close()