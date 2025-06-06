import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from mpl_toolkits.axes_grid1 import make_axes_locatable

chemin="/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/"
PIT_N3_file=f"{chemin}PIT3_E/wec_inputs/PIT3_kochin"

param_d=16
limite=20

Nb=2
config="X"
Nw=1
min_w=0.3
max_w=2
dof=6

FREE_S=True
COMP=False
data_FS="Module"
print(data_FS)
fichier_contour = "geo_barge.dat"  # Remplace par ton fichier des coordonnées du contour
Rcyl = 6.36  # Rayon du cylindre IT

Np=Nw*(dof*Nb+1)
print("Number of Problem :", Np)

ETA=0

if FREE_S : 
    while param_d<limite : 
        print("Distance = ", param_d)
        BEM_file=f"{chemin}PFE_Nemoh/MyTestCases/BEM_kochin_Nb{Nb}_{config}_d{param_d}"
        titre=f"Barge - Nb={Nb}, Nw={Nw}, dof={dof}, Ndir=1, d={param_d}"


        PROBLEMS = [f"{i:05d}" for i in range(1, Np+1)]

        for i, num_pb in enumerate(PROBLEMS) : 
            if (int(num_pb) - 1) % (Np // Nw) == 0:
                print("Problem n°", num_pb, "-> Diffraction")
                CAS="Diffraction"
            else:
                print("Problem n°", num_pb, "-> Radiation")
                CAS="Radiation"
            fichier_FS = f"{BEM_file}/results/freesurface.{num_pb}.dat"
            
                # === Charger les données ===
            with open(fichier_FS, "r") as f:
                lignes = f.readlines()

            points = []
            connect = []
            lecture_conn = False

            for ligne in lignes:
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

            ETA=ETA+Z

        # === Charger les données du contour de la surface ===
        data_contour = np.loadtxt(fichier_contour)

        # Extraire les coordonnées X et Y du contour
        X_contour = data_contour[:, 1]  # Coordonnée X
        Y_contour = data_contour[:, 2]  # Coordonnée Y

        # === Affichage ===
        fig, ax = plt.subplots(figsize=(8, 6))
        tri = Triangulation(X, Y, triangles)
        tcf = ax.tricontourf(tri, ETA, levels=100, cmap="viridis")
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
            
            # Tracer le contour de la surface
        if (config=="X" or config=="Y") :
            ax.plot(X_contour, Y_contour, color='r', linestyle='-', lw=1, label="barge")

        ax.set_aspect("equal")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_title(f"Free surface Potentiel ({data_FS})")
        # plt.legend()
        plt.tight_layout()
        plt.savefig(f"VIS_FS_{data_FS}_Nb{Nb}_{config}_d{param_d}_Np{Np}.png", dpi=300)
        plt.close()


        fichier_Ks = f"{BEM_file}/results/WaveField_S.tec"
        fichier_Ki = f"{BEM_file}/results/WaveField_I.tec"
        fichier_K = f"{BEM_file}/results/WaveField.tec"
        
        Files=[fichier_K, fichier_Ki, fichier_Ks]
        nom=["eta", "etaI", "etaS"]

        for ind, file in enumerate(Files) :
            # === Charger les données ===
            points = []
            connect = []
            lecture_conn = False
            with open(file, "r") as f:
                # lignes = f.readlines()

                for ligne in f:
                    ligne = ligne.strip()
                    if not ligne or ligne.startswith("VARIABLES") or ligne.startswith("ZONE"):
                        continue
                    vals = ligne.split()

                    if not lecture_conn and len(vals) == 6:
                        try:
                            points.append([float(v) for v in vals])
                        except ValueError:
                            continue
                    elif len(vals) == 4:
                        try:
                            connect.append([int(v) - 1 for v in vals])
                            lecture_conn = True
                        except ValueError:
                            continue
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
                Z = points[:, 2]  # abs(eta)
            elif data_FS=="Phase" :
                Z = points[:, 3]
            elif data_FS=="Reel" :
                Z = points[:, 4]
            elif data_FS=="Imag" :
                Z = points[:, 5]


            # === Charger les données du contour de la surface ===
            data_contour = np.loadtxt(fichier_contour)

            # Extraire les coordonnées X et Y du contour
            X_contour = data_contour[:, 1]  # Coordonnée X
            Y_contour = data_contour[:, 2]  # Coordonnée Y

            # === Affichage ===
            fig, ax = plt.subplots(figsize=(8, 6))
            tri = Triangulation(X, Y, triangles)
            tcf = ax.tricontourf(tri, Z, levels=100, cmap="viridis", vmax=max(ETA))
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

            # === Tracer le cylindre ===
            theta = np.linspace(0, 2 * np.pi, 200)
            if (config=="X") :
                for i in range(Nb) : # à génélariser !!!!
                    x_cyl = i*(2*Rcyl+param_d) + Rcyl * np.cos(theta)
                    y_cyl = Rcyl * np.sin(theta)
                    ax.plot(x_cyl, y_cyl, color="white", linestyle='--', linewidth=1)
            
                # Tracer le contour de la surface
            if (config=="X" or config=="Y") :
                ax.plot(X_contour, Y_contour, color='r', linestyle='-', lw=1, label="barge")

            ax.set_aspect("equal")
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            ax.set_title(f"Free surface Potentiel ({data_FS}) - {nom[ind]}")
            # plt.legend()
            plt.tight_layout()
            plt.savefig(f"VIS_{nom[ind]}_{data_FS}_Nb{Nb}_{config}_d{param_d}_Np{Np}.png", dpi=300)
            plt.close()


            if COMP : 
                # === Affichage ===
                fig, ax = plt.subplots(figsize=(8, 6))
                tri = Triangulation(X, Y, triangles)
                tcf = ax.tricontourf(tri, np.abs(ETA-Z)/ETA, levels=100, cmap="viridis")
                if (config=="Y") :
                    fig.colorbar(tcf, ax=ax, label="|ηF-ηK|/ηF")
                    # plt.legend()
                else : 
                    # Créer un axe pour la colorbar au-dessus
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("bottom", size="10%", pad=0.5)
                    # Ajouter la colorbar horizontale
                    cbar = fig.colorbar(tcf, cax=cax, orientation='horizontal')
                    cbar.set_label("|ηF-ηK|/ηF")
                    # Déplacer les ticks et le label en haut
                    cax.xaxis.set_ticks_position('bottom')
                    cax.xaxis.set_label_position('bottom')

                # === Tracer le cylindre ===
                theta = np.linspace(0, 2 * np.pi, 200)
                if (config=="X") :
                    for i in range(Nb) : # à génélariser !!!!
                        x_cyl = i*(2*Rcyl+param_d) + Rcyl * np.cos(theta)
                        y_cyl = Rcyl * np.sin(theta)
                        ax.plot(x_cyl, y_cyl, color="white", linestyle='--', linewidth=1)
                
                    # Tracer le contour de la surface
                if (config=="X" or config=="Y") :
                    ax.plot(X_contour, Y_contour, color='r', linestyle='-', lw=1, label="barge")

                ax.set_aspect("equal")
                ax.set_xlabel("X")
                ax.set_ylabel("Y")
                ax.set_title(f"Free surface Potentiel ({data_FS}) - {nom[ind]}")
                # plt.legend()
                plt.tight_layout()
                plt.savefig(f"VIS_COMP_{nom[ind]}_{data_FS}_Nb{Nb}_{config}_d{param_d}_Np{Np}.png", dpi=300)
                plt.close()


        param_d=param_d*2