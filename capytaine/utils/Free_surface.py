import matplotlib.pyplot as plt
import numpy as np
import os
import capytaine as cpt
from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
import cmath

def PLOT_FS(FSE, num, total, grid, folder):
    plt.pcolormesh(grid[0], grid[1], np.real(FSE))
    plt.xlabel("x")
    plt.ylabel("y")
    os.makedirs(folder, exist_ok=True)
    if total:
        filename = os.path.join(folder, f"C_FS_tot_numW{num}.png")
    else:
        filename = os.path.join(folder, f"C_FS_pb{num}.png")
    plt.savefig(filename)
    plt.close()  # bon réflexe pour libérer la figure
    
def SAVE_FS_TXT(FSE, num, total, grid, folder):
    X_grid, Y_grid = grid  # grid[0] = X, grid[1] = Y

    os.makedirs(folder, exist_ok=True)
    filename = f"C_FS_{'tot_numW' if total else 'pb'}{num}.txt"
    file = os.path.join(folder, filename)
    with open(file, 'w') as f:
        for i in range(X_grid.shape[0]):
            for j in range(X_grid.shape[1]):
                x = X_grid[i, j]
                y = Y_grid[i, j]
                Re = np.real(FSE[i, j])
                Im = np.imag(FSE[i, j])
                Abs = np.abs(FSE[i, j])
                Ph = cmath.phase(FSE[i, j])
                f.write(f"{x:.6e} {y:.6e} {Abs:.6e} {Ph:.6e} {Re:.6e} {Im:.6e}\n")
                             