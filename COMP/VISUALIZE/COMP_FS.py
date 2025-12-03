import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import Normalize
from matplotlib.ticker import FormatStrFormatter
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

    # === Affichage ===
    fig, ax = plt.subplots(figsize=(8, 6))
    if NOM == "PIT" or triangles is None:
        # Scatter plot pour les points (car pas de connectivités)
        # sc = ax.scatter(X, Y, c=Z, cmap="viridis", s=15, marker='o')
        tri = Triangulation(X, Y)
    else:
        # Triangulation classique quand on a les connectivités
        tri = Triangulation(X, Y, triangles)
    # tcf = ax.tricontourf(tri, Z, levels=100, cmap="plasma", vmin=0, vmax=Zmax)
    norm = Normalize(vmin=0, vmax=Zmax)
    levels = np.linspace(0, Zmax, 100)
    tcf = ax.tricontourf(tri, Z, levels=levels, cmap="plasma", norm=norm)
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
        cbar.ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))


    # === Tracer le cylindre ===
    theta = np.linspace(0, 2 * np.pi, 200)
    if (config=="X" or config=="C") :
    	# for dist in [[0, 0], [param_d,0]]:
    	for dist in [[-param_d/2,-param_d/2], [param_d/2,-param_d/2], [param_d/2,param_d/2], [-param_d/2,param_d/2]]:
            x_cyl = dist[0] + Rcyl * np.cos(theta)
            y_cyl = dist[1] + Rcyl * np.sin(theta)
        # for i in range(Nb) : # à génélariser !!!!
            # x_cyl = i*(2*Rcyl+param_d) + Rcyl * np.cos(theta)
            # y_cyl = Rcyl * np.sin(theta)
            ax.plot(x_cyl, y_cyl, color="white", linestyle='--', linewidth=1, label="outer cylinder")
    if (config=="Y") :
        for i in range(Nb) : # à génélariser !!!!
            y_cyl = i*(2*Rcyl+param_d) + Rcyl * np.sin(theta)
            x_cyl = Rcyl * np.cos(theta)
            ax.plot(x_cyl, y_cyl, color="white", linestyle='--', linewidth=1, label=f"cylinder {i+1}")
    
        # Tracer le contour de la surface
    if (config=="X" or config=="C") :
        if mesh=="barge":
            # === Charger les données du contour de la surface ===
            data_contour = np.loadtxt(fichier_contour)
            # Extraire les coordonnées X et Y du contour
            X_contour = data_contour[:, 1]  # Coordonnée X
            Y_contour = data_contour[:, 2]  # Coordonnée Y
            ax.plot(X_contour, Y_contour, color='r', linestyle='-', lw=1, label="barge")
        elif mesh=="CylR3" :
            a=3
            # for dist in [[0, 0], [param_d,0]]:
            for dist in [[-param_d/2,-param_d/2], [param_d/2,-param_d/2], [param_d/2,param_d/2], [-param_d/2,param_d/2]]:
                x_cyl = dist[0] + a * np.cos(theta)
                y_cyl = dist[1] + a * np.sin(theta)
                ax.plot(x_cyl, y_cyl, color="black", linestyle='-', linewidth=1, label="body")
            # x_cyl2 = param_d + a * np.cos(theta)
            # ax.plot(x_cyl2, y_cyl, color="black", linestyle='-', linewidth=1, label="body")
        elif mesh=="BargeX6" :
            a=3
            x_cyl = [-a, a, a, -a, -a]
            y_cyl = [a, a, -a, -a, a]
            ax.plot(x_cyl, y_cyl, color="black", linestyle='-', linewidth=1, label="body")
    ax.set_aspect("equal")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title(titre)
    # plt.legend()
    plt.tight_layout()
    plt.savefig(titreFILE, dpi=300)
    plt.close()

chemin="/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/"

param_d=12
limite=param_d+0.1

Nb=4
config="C"
Nw=1
min_w=3.20515
max_w=3.20515
dof=1
mesh="CylR3"
data_FS="Module"
fichier_contour = "geo_barge.dat"  # Remplace par ton fichier des coordonnées du contour
a=3
if mesh=="BargeX6":
    Rcyl = 4.48  # Rayon du cylindre IT
else: 
    Rcyl=3.25
g=9.81
Nbeta=1
Np=Nw*(dof*Nb+Nbeta)
print("Number of Problem :", Np)
Nemoh=True
# num : 1 à 6 dof body 1 ; 7 à 12 dof body 2 ; 13 probleme diff beta=0
plot_1pb=False
type_p=1
num_w=1
Lambda=5
Ne=6
layout="Nb4_C_"
print("Ne=", Ne)
if Ne>0 :
    E=f"E{Ne}_"
else : 
    E=""
if mesh=="MNcylR3":
    if Lambda==3:
        min_w=2.616997
        max_w=2.616997
    elif Lambda==10:
        min_w=1.433383
        max_w=1.433383
    elif Lambda==30:
        min_w=0.815111
        max_w=0.515111
else : 
    if Lambda==1:
        min_w=4.53277
        max_w=4.53277
    elif Lambda==2:
        min_w=3.20515
        max_w=3.20515
    elif Lambda==5:
        min_w=2.02665
        max_w=2.02665
    elif Lambda==10:
        min_w=1.42401
        max_w=1.42401

ErrorFILE=f"Error_N{Ne}.txt"

omega = np.linspace(min_w, max_w, Nw)  # vecteur des fréquences
omega_per_pb = np.repeat(omega, Np // Nw) 
while param_d<=limite : 
    print("Distance = ", param_d)
    print("--> ", data_FS)
    if Nemoh:
        BEM_file=f"{chemin}PFE_Nemoh/MyTestCases/FS/{mesh}_surge/BEM_{mesh}_L{Lambda}a_surge_{layout}d{param_d}"
        PIT_file=f"{chemin}PIT3_E/FS_inputs/PIT3_{mesh}_surge_L{Lambda}a"
    else :
        BEM_file=f"{chemin}capytaine/MyTestCases/C_BEM_{mesh}_FS_L{Lambda}a"
        PIT_file=f"{chemin}capytaine/MyTestCases/CAP_{mesh}_FS_L{Lambda}a"

    # titre=f"Cylinder - Nb={Nb}, Nw={Nw}, dof={dof}, Ndir=1, d={param_d}"
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
        if Nemoh:
            fichier_BEM = f"{BEM_file}/results/freesurface.{num_pb}.dat"
            fichier_PIT = f"{PIT_file}/Motion/Global_{E}FS_{layout}d{param_d}.00_pb{num_pb}.dat"
            fichier_PIT2 = f"{PIT_file}/Motion/Global_{E}ETA_{layout}d{param_d}.00_pb{num_w:05d}.dat"
        else: 
            fichier_BEM = f"{BEM_file}/Motion/C_FS_pb{int(num_pb)}.txt"
            fichier_BEM2 = f"{BEM_file}/Motion/C_FS_tot_numW{int(num_w)}.txt"
            fichier_PIT = f"{PIT_file}/Motion/CapytaineIT_{E}FS_d{param_d}.00_pb{num_pb}.dat"
            fichier_PIT2 = f"{PIT_file}/Motion/CapytaineIT_{E}ETA_d{param_d}.00_pb{num_w:05d}.dat"
        
        print(f"w={omega_per_pb[int(num_pb)-1]:.{4}g}")
        if Nemoh:
            points, triangles = READ_BEM(fichier_BEM)
        else:
            points= np.loadtxt(fichier_BEM)
            triangles=None
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
        # for i in range(Nb):
        for dist in [[-param_d/2,-param_d/2], [param_d/2,-param_d/2], [param_d/2,param_d/2], [-param_d/2,param_d/2]]:
            cx = dist[0]
            cy = dist[1]
            dist = np.sqrt((X - cx)**2 + (Y - cy)**2)
            mask_total &= (dist > Rcyl-marge)  # Garder les points hors cylindre (distance strictement > Rcyl)
        # Appliquer le masque : points dans le cylindre deviennent 0, les autres restent Zpit
        Zpit = np.where(mask_total, Zpit, 0.0)
        diff=np.abs(np.where(mask_total, Z, 0.0)-Zpit)
        error = np.where(mask_total, np.abs((Z - Zpit) / Z), np.nan)
        print("Zmax :", max(np.where(mask_total, Z, 0.0)), max(Zpit), max(diff))
        # print(max(np.where(mask_total, (Zpit/Z), 0.0)), min(np.where(mask_total, Zpit/Z, 100.0)), max(np.where(mask_total, np.abs(1-Zpit/Z), 0.0)))
        print("mean diff (%)", np.mean(diff)*100)
        print("mean error (%)", np.nanmean(error)*100)
        with open(ErrorFILE, "a") as file:
            file.write(f"{num_pb} - mean error (%) : {np.nanmean(error)*100} \n")
        error = np.nan_to_num(error, nan=0.0)
        Z = np.where(mask_total, Z, 0.0)
        if plot_1pb : 
            NOM="BEM"
            legende=f"η (m) - {NOM}"
            Zmax = max(np.max(np.array(Z)), np.max(np.array(Zpit)))
            titre = f"Free surface Potentiel ({data_FS}) - {CAS} Problem \n for w = {omega_per_pb[int(num_pb)-1]:.{4}g}rad/s"
            titreFILE = f"VIS_FS_{data_FS}_L{Lambda}_BEM_Nb{Nb}_{config}_d{param_d}_pb{num_pb}.png"
            PLOT(X, Y, Z, Zmax, triangles, NOM, legende, config, Nb, Rcyl, titre, titreFILE)

            NOM="PIT"
            legende=f"η (m) - {NOM}"
            titre = f"Free surface Potentiel ({data_FS}) - {CAS} Problem \n for w = {omega_per_pb[int(num_pb)-1]:.{4}g}rad/s"
            titreFILE = f"VIS_FS_{E}{data_FS}_L{Lambda}_PIT_Nb{Nb}_{config}_d{param_d}_pb{num_pb}.png"
            PLOT(X, Y, Zpit, Zmax, None, NOM, legende, config, Nb, Rcyl, titre, titreFILE)

            legende="|η_BEM - η_PIT|/A  (%) "
            titre = f"Relative error of free surface Potentiel ({data_FS}) \n  {CAS} Problem - for w = {omega_per_pb[int(num_pb)-1]:.{4}g}rad/s"
            titreFILE = f"VIS_Diff_FS_{E}{data_FS}_L{Lambda}_Nb{Nb}_{config}_d{param_d}_pb{num_pb}.png"
            # PLOT(X, Y, diff*100, max(diff)*100, triangles, "PIT", legende, config, Nb, Rcyl, titre, titreFILE)
            
            legende="|η_BEM - η_PIT|/η_BEM (%)"
            titre = f"Relative error of free surface Potentiel ({data_FS}) \n  {CAS} Problem - for w = {omega_per_pb[int(num_pb)-1]:.{4}g}rad/s"
            titreFILE = f"VIS_Error_FS_{E}{data_FS}_L{Lambda}_Nb{Nb}_{config}_d{param_d}_pb{num_pb}.png"
            PLOT(X, Y, error*100, max(error)*100, triangles, "PIT", legende, config, Nb, Rcyl, titre, titreFILE)
            
            if CAS=="Diffraction":
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
                NOM="BEM"
                legende=f"η (m) - {NOM}"
                Zmax2 = max(np.max(np.array(Z+Zpit2)), np.max(np.array(Zpit+Zpit2)))
                titre = f"Free surface Potentiel ({data_FS}) - {CAS} Problem \n for w = {omega_per_pb[int(num_pb)-1]:.{4}g}rad/s"
                titreFILE = f"VIS_FS_{data_FS}_L{Lambda}_BEM_Nb{Nb}_{config}_d{param_d}_pb{num_pb}_withI.png"
                PLOT(X, Y, Z+Zpit2, Zmax2, triangles, NOM, legende, config, Nb, Rcyl, titre, titreFILE)

                NOM="PIT"
                legende=f"η (m) - {NOM}"
                titre = f"Free surface Potentiel ({data_FS}) - {CAS} Problem \n for w = {omega_per_pb[int(num_pb)-1]:.{4}g}rad/s"
                titreFILE = f"VIS_FS_{E}{data_FS}_L{Lambda}_PIT_Nb{Nb}_{config}_d{param_d}_pb{num_pb}_withI.png"
                PLOT(X, Y, Zpit+Zpit2, Zmax2, None, NOM, legende, config, Nb, Rcyl, titre, titreFILE)
                
                error2 = np.where(mask_total, np.abs((Z - Zpit) / (Z+Zpit2)), np.nan)
                print("mean error (%) with I", np.nanmean(error2)*100)
                error2 = np.nan_to_num(error2, nan=0.0)
                legende="|η_BEM - η_PIT|/η_BEM (%)"
                titre = f"Relative error of free surface Potentiel ({data_FS}) \n  {CAS} Problem - for w = {omega_per_pb[int(num_pb)-1]:.{4}g}rad/s"
                titreFILE = f"VIS_Error_FS_{E}{data_FS}_L{Lambda}_Nb{Nb}_{config}_d{param_d}_pb{num_pb}_withI.png"
                PLOT(X, Y, error2*100, max(error2)*100, triangles, "PIT", legende, config, Nb, Rcyl, titre, titreFILE)
            
        if not plot_1pb : 
            Z_tot+=Z
            Zpit_tot+=Zpit

    if not plot_1pb :
        print("--------- total free surface -----")
        maxZ_tot=max(Z_tot)
        maxZpit_tot=max(Zpit_tot)
        NOM="BEM"
        legende=f"η (m) - {NOM}"
        titre = f"Free surface Potentiel ({data_FS}) \n for w = {omega[num_w-1]:.{4}g}rad/s"
        titreFILE = f"FS_{data_FS}_BEM_Nb{Nb}_{config}_d{param_d}_numW{num_w}.png"
        # PLOT(X, Y, Z_tot, max(Z_tot), triangles, NOM, legende, config, Nb, Rcyl, titre, titreFILE)
        NOM="PIT"
        legende=f"η (m) - {NOM}"
        titre = f"Free surface Potentiel ({data_FS})  \n for w = {omega[num_w-1]:.{4}g}rad/s"
        titreFILE = f"FS_{data_FS}_{E}PIT_Nb{Nb}_{config}_d{param_d}_numW{num_w}.png"
        # PLOT(X, Y, Zpit_tot, max(Zpit_tot), None, NOM, legende, config, Nb, Rcyl, titre, titreFILE)

        diff=np.abs(np.where(mask_total, Z_tot, 0.0)-Zpit_tot)
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio_nan = np.divide(np.abs(Z_tot - Zpit_tot), Z_tot, out=np.full_like(Z_tot, np.nan), where=(Z_tot != 0))
            error_nan = np.where(mask_total, ratio_nan, np.nan)
            mean_error = np.nanmean(error_nan)
        error = np.nan_to_num(error_nan, nan=0.0)  # remplace NaN par 0
        legende="|η_BEM - η_PIT|/A (%)"
        titre = f"Absolute error of free surface Potentiel ({data_FS})  \n for w = {omega[num_w-1]:.{4}g}rad/s"
        titreFILE = f"FS_Diff_{data_FS}_{E}Nb{Nb}_{config}_d{param_d}_numW{num_w}.png"
        # PLOT(X, Y, diff*100, max(diff)*100, triangles, "PIT", legende, config, Nb, Rcyl, titre, titreFILE)
        print("mean diff (%)", np.mean(diff)*100)
        legende="|η_BEM - η_PIT|/η_BEM (%)"
        titre = f"Relative error of free surface Potentiel ({data_FS})  \n for w = {omega[num_w-1]:.{4}g}rad/s"
        titreFILE = f"FS_Error_{data_FS}_{E}Nb{Nb}_{config}_d{param_d}_numW{num_w}.png"
        # PLOT(X, Y, error*100, max(error)*100, triangles, "PIT", legende, config, Nb, Rcyl, titre, titreFILE)
        print("mean error (%)", mean_error*100)
        with open(ErrorFILE, "a") as file:
            file.write(f"tot - mean error (%) {mean_error*100} \n")
        print("--------- total free surface with incident wave -----")
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
        Z_tot=Z_tot+Zpit2
        Zpit_tot=Zpit_tot+Zpit2
        maxZ_tot=max(Z_tot)
        maxZpit_tot=max(Zpit_tot)
        NOM="BEM"
        legende=f"η (m) - {NOM}"
        titre = f"Free surface Potentiel ({data_FS}) \n for w = {omega[num_w-1]:.{4}g}rad/s"
        titreFILE = f"FS_TOT_{data_FS}_BEM_Nb{Nb}_{config}_d{param_d}_numW{num_w}.png"
        PLOT(X, Y, Z_tot, max(Z_tot), triangles, NOM, legende, config, Nb, Rcyl, titre, titreFILE)
        NOM="PIT"
        legende=f"η (m) - {NOM}"
        titre = f"Free surface Potentiel ({data_FS})  \n for w = {omega[num_w-1]:.{4}g}rad/s"
        titreFILE = f"FS_TOT_{data_FS}_{E}PIT_Nb{Nb}_{config}_d{param_d}_numW{num_w}.png"
        PLOT(X, Y, Zpit_tot, max(Zpit_tot), None, NOM, legende, config, Nb, Rcyl, titre, titreFILE)

        diff=np.abs(np.where(mask_total, Z_tot, 0.0)-Zpit_tot)
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio_nan = np.divide(np.abs(Z_tot - Zpit_tot), Z_tot, out=np.full_like(Z_tot, np.nan), where=(Z_tot != 0))
            error_nan = np.where(mask_total, ratio_nan, np.nan)
            mean_error = np.nanmean(error_nan)
        error = np.nan_to_num(error_nan, nan=0.0)  # remplace NaN par 0
        legende="|η_BEM - η_PIT|/A (%)"
        titre = f"Absolute error of free surface Potentiel ({data_FS})  \n for w = {omega[num_w-1]:.{4}g}rad/s"
        titreFILE = f"FS_TOT_Diff_{data_FS}_{E}Nb{Nb}_{config}_d{param_d}_numW{num_w}.png"
        # PLOT(X, Y, diff*100, max(diff)*100, triangles, "PIT", legende, config, Nb, Rcyl, titre, titreFILE)
        print("mean diff (%)", np.mean(diff)*100)
        legende="|η_BEM - η_PIT|/η_BEM (%)"
        titre = f"Relative error of free surface Potentiel ({data_FS})  \n for w = {omega[num_w-1]:.{4}g}rad/s"
        titreFILE = f"FS_TOT_Error_{data_FS}_{E}Nb{Nb}_{config}_d{param_d}_numW{num_w}.png"
        PLOT(X, Y, error*100, max(error)*100, triangles, "PIT", legende, config, Nb, Rcyl, titre, titreFILE)
        print("mean error (%)", mean_error*100)
        with open(ErrorFILE, "a") as file:
            file.write(f"tot with I - mean error (%) {mean_error*100} \n")
        legende="η_I"
        titre = f"Free surface Incident Potentiel ({data_FS})  \n for w = {omega[num_w-1]:.{4}g}rad/s"
        titreFILE = f"ETA_I_{data_FS}_Nb{Nb}_{config}_d{param_d}_numW{num_w}.png"
        # PLOT(X, Y, Zpit, max(Zpit), None, NOM, legende, config, Nb, Rcyl, titre, titreFILE)

    param_d=param_d+1
