import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_cylinder(ax, config, Nb, Rcyl, param_d):
    # === Tracer le cylindre ===
    theta = np.linspace(0, 2 * np.pi, 200)
    if (config=="X") :
        for i in range(Nb) : # à génélariser !!!!
            x_cyl = i*(2*Rcyl+param_d) + Rcyl * np.cos(theta)
            y_cyl = Rcyl * np.sin(theta)
            ax.plot(x_cyl, y_cyl, color="white", linestyle='--', linewidth=1)
    if (config=="Y") :
        for i in range(-(Nb-2), Nb-1, 1) : # à génélariser !!!!
            y_cyl = i*(2*Rcyl+param_d) + Rcyl * np.sin(theta)
            x_cyl = Rcyl * np.cos(theta)
            ax.plot(x_cyl, y_cyl, color="white", linestyle='--', linewidth=1, label=f"cylinder {i+1}")
    if (config=="T") :
        y_cyl = (2*Rcyl+param_d)/2 + Rcyl * np.sin(theta)
        x_cyl = Rcyl * np.cos(theta)
        ax.plot(x_cyl, y_cyl, color="white", linestyle='--', linewidth=1, label=f"cylinder {i+1}")
        y_cyl = -(2*Rcyl+param_d)/2 + Rcyl * np.sin(theta)
        x_cyl = Rcyl * np.cos(theta)
        ax.plot(x_cyl, y_cyl, color="white", linestyle='--', linewidth=1, label=f"cylinder {i+1}")
        x_cyl = (2*Rcyl+param_d)*np.sqrt(3)/2 + Rcyl * np.sin(theta)
        y_cyl = Rcyl * np.cos(theta)
        ax.plot(x_cyl, y_cyl, color="white", linestyle='--', linewidth=1, label=f"cylinder {i+1}")
    if (config=="Tinv") :
        y_cyl = Rcyl * np.sin(theta)
        x_cyl = Rcyl * np.cos(theta)
        ax.plot(x_cyl, y_cyl, color="white", linestyle='--', linewidth=1, label=f"cylinder {i+1}")
        y_cyl = -(2*Rcyl+param_d)/2 + Rcyl * np.sin(theta)
        x_cyl = (2*Rcyl+param_d)*np.sqrt(3)/2 + Rcyl * np.cos(theta)
        ax.plot(x_cyl, y_cyl, color="white", linestyle='--', linewidth=1, label=f"cylinder {i+1}")
        x_cyl = (2*Rcyl+param_d)*np.sqrt(3)/2 + Rcyl * np.sin(theta)
        y_cyl = (2*Rcyl+param_d)/2 + Rcyl * np.cos(theta)
        ax.plot(x_cyl, y_cyl, color="white", linestyle='--', linewidth=1, label=f"cylinder {i+1}")
    if (config=="C") :
        y_cyl = -(2*Rcyl+param_d)/2+Rcyl * np.sin(theta)
        x_cyl = -(2*Rcyl+param_d)/2+Rcyl * np.cos(theta)
        ax.plot(x_cyl, y_cyl, color="white", linestyle='--', linewidth=1, label=f"cylinder {i+1}")
        y_cyl = -(2*Rcyl+param_d)/2 + Rcyl * np.sin(theta)
        x_cyl = (2*Rcyl+param_d)/2 + Rcyl * np.cos(theta)
        ax.plot(x_cyl, y_cyl, color="white", linestyle='--', linewidth=1, label=f"cylinder {i+1}")
        y_cyl = (2*Rcyl+param_d)/2 + Rcyl * np.sin(theta)
        x_cyl = -(2*Rcyl+param_d)/2 + Rcyl * np.cos(theta)
        ax.plot(x_cyl, y_cyl, color="white", linestyle='--', linewidth=1, label=f"cylinder {i+1}")
        y_cyl = (2*Rcyl+param_d)/2+Rcyl * np.sin(theta)
        x_cyl = (2*Rcyl+param_d)/2+Rcyl * np.cos(theta)
        ax.plot(x_cyl, y_cyl, color="white", linestyle='--', linewidth=1, label=f"cylinder {i+1}")
    
    # === Charger les données du contour de la surface ===
    fichier_contour = "geo_barge.dat" 
    data_contour = np.loadtxt(fichier_contour)

    # Extraire les coordonnées X et Y du contour
    X_contour = data_contour[:, 1]  # Coordonnée X
    Y_contour = data_contour[:, 2]  # Coordonnée Y

    if (config=="X" or config=="Y") :
        ax.plot(X_contour, Y_contour, color='r', linestyle='-', lw=1, label="barge")

    return 


chemin="/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/"
PIT_N3_file=f"{chemin}PIT3_E/wec_inputs/PIT3_barge_w2pi"

param_d=16
limite=20

Nb=2
config="X"
Nw=1
min_w=0.3
max_w=2
dof=6


CYL=False
data_CYL="Module"
FREE_S=True
FREE_S_K=True
data_FS=data_CYL
 # Remplace par ton fichier des coordonnées du contour
Rcyl = 6.36  # Rayon du cylindre IT

Np=Nw*(dof*Nb+1)
print("Number of Problem :", Np)

omega = np.linspace(min_w, max_w, Nw)  # vecteur des fréquences
omega_per_pb = np.repeat(omega, Np // Nw) 
# num : 1 à 6 dof body 1 ; 7 à 12 dof body 2 ; 13 probleme diff beta=0
plot_1pb=False
type_p=1
num_w=1

while param_d<limite : 
    print("Distance = ", param_d)
    # BEM_file=f"{chemin}PFE_Nemoh/MyTestCases/BEM_barge_FS_d1/BEM_FS_Nb{Nb}_{config}_d{param_d}"
    BEM_file=f"{chemin}PFE_Nemoh/MyTestCases/BEM_kochin_Nb{Nb}_{config}_d{param_d}"
    titre=f"Barge - Nb={Nb}, Nw={Nw}, dof={dof}, Ndir=1, d={param_d}"
    centres = [(0, 0), (param_d, 0)]   # liste des centres (xg, yg)
    if FREE_S : 
        if plot_1pb : 
            range_plot = [type_p, Np + 1, int(Np/Nw)]
        else :
            range_plot = [(num_w-1)*int(Np/Nw)+1, (num_w)*int(Np/Nw)+1, 1]
    else :
        range_plot = [1, 12, 1]
    PROBLEMS = [f"{i:05d}" for i in range(range_plot[0], range_plot[1], range_plot[2])]
    for i, num_pb in enumerate(PROBLEMS) : 
        if (int(num_pb) - 1) % (Np // Nw) == 0:
            print("Problem n°", num_pb, "-> Diffraction")
            CAS="Diffraction"
        else:
            print("Problem n°", num_pb, "-> Radiation")
            CAS="Radiation"
        fichier_FS = f"{BEM_file}/results/freesurface.{num_pb}.dat"
        fichier_FS_K = f"{BEM_file}/results/WaveField_{num_pb}.dat"
        filename_CYL = f"{PIT_N3_file}/results/cylsurface.{num_pb}.dat"

        if CYL : 
            # Lire tout le fichier
            with open(filename_CYL, "r") as f:
                lines = f.readlines()

            # Trouver début des données et du maillage
            data_start = 0
            for i, line in enumerate(lines):
                if "ZONE" in line:
                    data_start = i + 1
                    break

            # Lire les points (N = 899)
            n_points = 899
            points_data = np.loadtxt(lines[data_start:data_start + n_points])
            x, y, z = points_data[:, 0], points_data[:, 1], points_data[:, 2]
            module, phase = points_data[:, 3], points_data[:, 4]
            potentiel = module * np.exp(1j * phase)

            # Créer les points
            points = np.column_stack((x, y, z))

            # Lire les connectivités (E = 868 quadrilatères)
            n_elems = 868
            conn_lines = lines[data_start + n_points:data_start + n_points + n_elems]
            conn = [list(map(int, l.strip().split())) for l in conn_lines]
            conn = np.array(conn) - 1  # Passage de Fortran (1-based) à Python (0-based)

            # Convertir en format PyVista
            cells = np.hstack([[4, *quad] for quad in conn])  # 4 pour quadrilatère
            celltypes = np.full(n_elems, pv.CellType.QUAD, dtype=np.uint8)

            # Créer la structure maillée
            grid = pv.UnstructuredGrid(cells, celltypes, points)

            # Ajouter les champs scalaires
            grid["Module"] = np.abs(potentiel)
            grid["Phase"] = np.angle(potentiel)
            grid["Reel"] = potentiel.real
            grid["Imag"] = potentiel.imag

            # Affichage
            plotter = pv.Plotter()
            plotter.add_mesh(grid, scalars=f"{data_CYL}", cmap="plasma", show_edges=False)
            plotter.add_axes()
            plotter.show(screenshot=f"VIS_CYL_potential_{data_CYL}_Nb{Nb}_{config}_pb{num_pb}.png")

        if FREE_S : 
            # === Charger les données ===
            with open(fichier_FS, "r") as f:
                lignes = f.readlines()

            points = []
            connect = []
            lecture_conn = False

            for ligne in lignes:
                if ligne.startswith('VARIABLES') or ligne.startswith('ZONE'):
                    continue
                vals = ligne.strip().split()
                if not lecture_conn and len(vals) == 7:
                    points.append([float(v) for v in vals])
                elif len(vals) == 4:
                    connect.append([int(v) - 1 for v in vals])
                    lecture_conn = True

            points = np.array(points)
            quads = np.array(connect)

            # === Transformer quadrilatères en triangles ===
            # Chaque quad : [n1, n2, n3, n4] → triangles : [n1, n2, n3] et [n1, n3, n4]
            triangles = []
            for q in quads:
                triangles.append([q[0], q[1], q[2]])
                triangles.append([q[0], q[2], q[3]])
            triangles = np.array(triangles)

            # === Coordonnées + données ===
            X = points[:, 0]
            Y = points[:, 1]
            if data_FS=="Module" : 
                Z = points[:, 3]  # abs(eta)
            elif data_FS=="Phase" :
                Z = points[:, 4]
            elif data_FS=="Reel" :
                Z = points[:, 5]
            elif data_FS=="Imag" :
                Z = points[:, 6]
            
            masque = np.ones_like(Z, dtype=bool)
            for (xg, yg) in centres:
                masque &= ((X - xg)**2 + (Y - yg)**2 > Rcyl**2)
            Z_out = Z[masque]

            # === Affichage ===
            fig, ax = plt.subplots(figsize=(8, 6))
            tri = Triangulation(X, Y, triangles)
            tcf = ax.tricontourf(tri, Z, levels=100, cmap="viridis", vmax=max(Z))
            if (config=="Y") :
                fig.colorbar(tcf, ax=ax, label="|η| (m)")
                # plt.legend()
            else : 
                # Créer un axe pour la colorbar au-dessus
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("bottom", size="10%", pad=0.5)
                # Ajouter la colorbar horizontale
                cbar = fig.colorbar(tcf, cax=cax, orientation='horizontal')
                cbar.set_label("|η| (m)")
                # Déplacer les ticks et le label en haut
                cax.xaxis.set_ticks_position('bottom')
                cax.xaxis.set_label_position('bottom')

            plot_cylinder(ax, config, Nb, Rcyl, param_d)
            
            ax.set_aspect("equal")
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            ax.set_title(f"Free surface Potentiel ({data_FS}) - {CAS} Problem \n for w = {omega_per_pb[int(num_pb)-1]:.{4}g}rad/s")
            # plt.legend()
            plt.tight_layout()
            plt.savefig(f"VIS_FS_{data_FS}_Nb{Nb}_{config}_d{param_d}_pb{num_pb}.png", dpi=300)
            plt.close()

        if FREE_S_K : 
            # === Charger les données ===
            with open(fichier_FS_K, "r") as f:
                lignes = f.readlines()

            points_K = []
            connect = []
            lecture_conn = False

            for ligne in lignes:
                if ligne.startswith('VARIABLES') or ligne.startswith('ZONE'):
                    continue
                vals = ligne.strip().split()
                if not lecture_conn and len(vals) == 6:
                    points_K.append([float(v) for v in vals])
                elif len(vals) == 4:
                    connect.append([int(v) - 1 for v in vals])
                    lecture_conn = True

            points_K = np.array(points_K)
            quads = np.array(connect)

            # === Transformer quadrilatères en triangles ===
            # Chaque quad : [n1, n2, n3, n4] → triangles : [n1, n2, n3] et [n1, n3, n4]
            triangles = []
            for q in quads:
                triangles.append([q[0], q[1], q[2]])
                triangles.append([q[0], q[2], q[3]])
            triangles = np.array(triangles)

            # === Coordonnées + données ===
            X = points_K[:, 0]
            Y = points_K[:, 1]
            if data_FS=="Module" : 
                Z_K = points_K[:, 2]  # abs(eta)
            elif data_FS=="Phase" :
                Z_K = points_K[:, 3]
            elif data_FS=="Reel" :
                Z_K = points_K[:, 4]
            elif data_FS=="Imag" :
                Z_K = points_K[:, 5]

            masque = np.ones_like(Z_K, dtype=bool)
            for (xg, yg) in centres:
                masque &= ((X - xg)**2 + (Y - yg)**2 > Rcyl**2)
            Z_Kout = Z_K[masque]
            # === Affichage ===
            fig, ax = plt.subplots(figsize=(8, 6))
            tri = Triangulation(X, Y, triangles)
            tcf = ax.tricontourf(tri, Z_K, levels=100, cmap="viridis", vmax=max(Z_K))
            if (config=="Y") :
                fig.colorbar(tcf, ax=ax, label="|η| (m)")
                # plt.legend()
            else : 
                # Créer un axe pour la colorbar au-dessus
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("bottom", size="10%", pad=0.5)
                # Ajouter la colorbar horizontale
                cbar = fig.colorbar(tcf, cax=cax, orientation='horizontal')
                cbar.set_label("|η| (m)")
                # Déplacer les ticks et le label en haut
                cax.xaxis.set_ticks_position('bottom')
                cax.xaxis.set_label_position('bottom')

            plot_cylinder(ax, config, Nb, Rcyl, param_d)
            
            ax.set_aspect("equal")
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            ax.set_title(f"Free surface Potentiel ({data_FS}) - {CAS} Problem \n for w = {omega_per_pb[int(num_pb)-1]:.{4}g}rad/s")
            # plt.legend()
            plt.tight_layout()
            plt.savefig(f"VIS_FS_K_{data_FS}_Nb{Nb}_{config}_d{param_d}_pb{num_pb}.png", dpi=300)
            plt.close()

    param_d=param_d*2