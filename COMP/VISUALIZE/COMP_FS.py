import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from mpl_toolkits.axes_grid1 import make_axes_locatable


def READ_BEM(file) : 
    with open(file, "r") as f:
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
    return points, triangles

def PLOT(X, Y, Z, Zmax, triangles, NOM, legende, config, Nb, Rcyl, titre, titreFILE) :

    # === Charger les données du contour de la surface ===
    data_contour = np.loadtxt(fichier_contour)
    # Extraire les coordonnées X et Y du contour
    X_contour = data_contour[:, 1]  # Coordonnée X
    Y_contour = data_contour[:, 2]  # Coordonnée Y

    # === Affichage ===
    fig, ax = plt.subplots(figsize=(8, 6))
    if NOM == "PIT":
        # Scatter plot pour les points (car pas de connectivités)
        # sc = ax.scatter(X, Y, c=Z, cmap="viridis", s=15, marker='o')
        tri = Triangulation(X, Y)
    else:
        # Triangulation classique quand on a les connectivités
        tri = Triangulation(X, Y, triangles)
    tcf = ax.tricontourf(tri, Z, levels=100, cmap="plasma", vmax=Zmax)

    if config == "Y":
        fig.colorbar(tcf, ax=ax, label="|η| (m)")
    else:
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="10%", pad=0.5)
        cbar = fig.colorbar(tcf, cax=cax, orientation='horizontal')
        cbar.set_label(f"{legende}")
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
        for i in range(Nb) : # à génélariser !!!!
            y_cyl = i*(2*Rcyl+param_d) + Rcyl * np.sin(theta)
            x_cyl = Rcyl * np.cos(theta)
            ax.plot(x_cyl, y_cyl, color="white", linestyle='--', linewidth=1, label=f"cylinder {i+1}")
    
        # Tracer le contour de la surface
    if (config=="X" or config=="Y") :
        ax.plot(X_contour, Y_contour, color='r', linestyle='-', lw=1, label="barge")

    ax.set_aspect("equal")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title(titre)
    # plt.legend()
    plt.tight_layout()
    plt.savefig(titreFILE, dpi=300)
    plt.close()

chemin="/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/"

param_d=16
limite=20

Nb=2
config="X"
Nw=3
min_w=0.3
max_w=2
dof=6

data_FS="Reel"
fichier_contour = "geo_barge.dat"  # Remplace par ton fichier des coordonnées du contour
Rcyl = 6.36  # Rayon du cylindre IT

Np=Nw*(dof*Nb+1)
print("Number of Problem :", Np)

omega = np.linspace(min_w, max_w, Nw)  # vecteur des fréquences
omega_per_pb = np.repeat(omega, Np // Nw) 
# num : 1 à 6 dof body 1 ; 7 à 12 dof body 2 ; 13 probleme diff beta=0
plot_1pb=True
type_p=1
num_w=2
Ne=0
if Ne>0 :
    E=f"E{Ne}_"
else : 
    E=""
while param_d<limite : 
    print("Distance = ", param_d)
    print("--> ", data_FS)
    BEM_file=f"{chemin}PFE_Nemoh/MyTestCases/BEM_FS_Nb{Nb}_{config}_d{param_d}/"
    PIT_file=f"{chemin}PIT3_E/wec_inputs/PIT3_barge_FS/"
    titre=f"Barge - Nb={Nb}, Nw={Nw}, dof={dof}, Ndir=1, d={param_d}"
    Z_tot=0
    Zpit_tot=0
    if plot_1pb : 
        range_plot = [type_p, Np + 1, int(Np/Nw)]
        print("Same problem")
    else :
        range_plot = [(num_w-1)*int(Np/Nw)+1, (num_w)*int(Np/Nw)+1, 1]
        print("Same frequency")
   
    PROBLEMS = [f"{i:05d}" for i in range(range_plot[0], range_plot[1], range_plot[2])]
    for i, num_pb in enumerate(PROBLEMS) : 
        if (int(num_pb) - 1) % (Np // Nw) == 0:
            print("Problem n°", num_pb, "-> Diffraction")
            CAS="Diffraction"
        else:
            print("Problem n°", num_pb, "-> Radiation")
            CAS="Radiation"
        fichier_BEM = f"{BEM_file}/results/freesurface.{num_pb}.dat"
        fichier_PIT = f"{PIT_file}/Motion/Global_{E}FS_d{param_d}.00_pb{num_pb}.dat"
        fichier_PIT2 = f"{PIT_file}/Motion/Global_{E}ETA_d{param_d}.00_pb{num_w:05d}.dat"

        
        points, triangles = READ_BEM(fichier_BEM)
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

        pointsPIT= np.loadtxt(fichier_PIT)
        X = pointsPIT[:, 0]
        Y = pointsPIT[:, 1]
        if data_FS=="Module" : 
            Zpit = pointsPIT[:, 2]  # abs(eta)
        elif data_FS=="Phase" :
            Zpit = pointsPIT[:, 3]
        elif data_FS=="Reel" :
            Zpit = pointsPIT[:, 4]
        elif data_FS=="Imag" :
            Zpit = pointsPIT[:, 5]

        mask_total = np.ones_like(Zpit, dtype=bool)
        marge=0.01*Rcyl
        for i in range(Nb):
            if config == "X":
                cx = i * (2 * Rcyl + param_d)
                cy = 0
            elif config == "Y":
                cx = 0
                cy = i * (2 * Rcyl + param_d)
            else:
                raise ValueError("config must be 'X' or 'Y'")

            dist = np.sqrt((X - cx)**2 + (Y - cy)**2)
            mask_total &= (dist > Rcyl-marge)  # Garder les points hors cylindre (distance strictement > Rcyl)
        # Appliquer le masque : points dans le cylindre deviennent 0, les autres restent Zpit
        Zpit = np.where(mask_total, Zpit, 0.0)
        diff=np.abs(np.where(mask_total, Z, 0.0)-Zpit)
        print("Zmax :", max(np.where(mask_total, Z, 0.0)), max(Zpit), max(diff))
        print(max(np.where(mask_total, (Zpit/Z), 0.0)), min(np.where(mask_total, Zpit/Z, 100.0)), max(np.where(mask_total, np.abs(1-Zpit/Z), 0.0)))
        Z = np.where(mask_total, Z, 0.0)
        if plot_1pb : 
            NOM="BEM"
            legende=f"η (m) - {NOM}"
            Zmax=max(np.max(Z), np.max(Zpit))
            titre = f"Free surface Potentiel ({data_FS}) - {CAS} Problem \n for w = {omega_per_pb[int(num_pb)-1]:.{4}g}rad/s"
            titreFILE = f"VIS_FS_{data_FS}_BEM_Nb{Nb}_{config}_d{param_d}_pb{num_pb}.png"
            PLOT(X, Y, Z, np.max(Z), triangles, NOM, legende, config, Nb, Rcyl, titre, titreFILE)

            NOM="PIT"
            legende=f"η (m) - {NOM}"
            titre = f"Free surface Potentiel ({data_FS}) - {CAS} Problem \n for w = {omega_per_pb[int(num_pb)-1]:.{4}g}rad/s"
            titreFILE = f"VIS_FS_{E}{data_FS}_PIT_Nb{Nb}_{config}_d{param_d}_pb{num_pb}.png"
            PLOT(X, Y, Zpit, np.max(Zpit), None, NOM, legende, config, Nb, Rcyl, titre, titreFILE)

            legende="|η_BEM - η_PIT|/A  (%) "
            titre = f"Relative error of free surface Potentiel ({data_FS}) \n  {CAS} Problem - for w = {omega_per_pb[int(num_pb)-1]:.{4}g}rad/s"
            titreFILE = f"VIS_Error_FS_{E}{data_FS}_Nb{Nb}_{config}_d{param_d}_pb{num_pb}.png"
            PLOT(X, Y, diff*100, max(diff)*100, triangles, "PIT", legende, config, Nb, Rcyl, titre, titreFILE)

        if not plot_1pb : 
            Z_tot+=Z
            Zpit_tot+=Zpit

    if not plot_1pb :
        NOM="BEM"
        legende=f"η (m) - {NOM}"
        titre = f"Free surface Potentiel ({data_FS}) \n for w = {omega[num_w-1]:.{4}g}rad/s"
        titreFILE = f"VIS_FS_{data_FS}_BEM_Nb{Nb}_{config}_d{param_d}_numW{num_w}.png"
        PLOT(X, Y, Z_tot, max(Z_tot), triangles, NOM, legende, config, Nb, Rcyl, titre, titreFILE)
        NOM="PIT"
        legende=f"η (m) - {NOM}"
        titre = f"Free surface Potentiel ({data_FS})  \n for w = {omega[num_w-1]:.{4}g}rad/s"
        titreFILE = f"VIS_FS_{data_FS}_PIT_Nb{Nb}_{config}_d{param_d}_numW{num_w}.png"
        PLOT(X, Y, Zpit_tot, max(Zpit_tot), None, NOM, legende, config, Nb, Rcyl, titre, titreFILE)

        diff=np.abs(np.where(mask_total, Z_tot, 0.0)-Zpit_tot)
        legende="|η_BEM - η_PIT| (m) "
        titre = f"Absolute error of free surface Potentiel ({data_FS})  \n for w = {omega[num_w-1]:.{4}g}rad/s"
        titreFILE = f"VIS_Error_FS_{data_FS}_Nb{Nb}_{config}_d{param_d}_numW{num_w}.png"
        PLOT(X, Y, diff, max(diff), triangles, "PIT", legende, config, Nb, Rcyl, titre, titreFILE)
        
        NOM="PIT"
        legende=f"η (m) - {NOM}"
        pointsPIT= np.loadtxt(fichier_PIT2, skiprows=1)
        X = pointsPIT[:, 0]
        Y = pointsPIT[:, 1]
        if data_FS=="Module" : 
            Zpit2 = pointsPIT[:, 2]  # abs(eta)
        elif data_FS=="Phase" :
            Zpit2 = pointsPIT[:, 3]
        elif data_FS=="Reel" :
            Zpit2 = pointsPIT[:, 4]
        elif data_FS=="Imag" :
            Zpit2 = pointsPIT[:, 5]
        Zpit2 = np.where(mask_total, Zpit2, 0.0)
        titre = f"Free surface Potentiel ({data_FS})  \n for w = {omega[num_w-1]:.{4}g}rad/s"
        titreFILE = f"VIS_ETA_{data_FS}_PIT_Nb{Nb}_{config}_d{param_d}_numW{num_w}.png"
        PLOT(X, Y, Zpit, max(Zpit), None, NOM, legende, config, Nb, Rcyl, titre, titreFILE)

    param_d=param_d*2
