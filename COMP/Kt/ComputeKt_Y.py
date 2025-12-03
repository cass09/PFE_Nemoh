import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import Normalize
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

def READ_BEM(file) : 
    with open(file, "r") as f:
        lignes = f.readlines()

    points, connect = [], []
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

    # Quadrilatères → triangles
    triangles = []
    for q in quads:
        triangles.append([q[0], q[1], q[2]])
        triangles.append([q[0], q[2], q[3]])
    return points, np.array(triangles)

def plot_eta_vs_x(X, Y, Z1_tot, Z2_tot, y_values, titre, titreFile, dimX, meshes):
    plt.figure(figsize=(8, 6))
    colors=['b', 'r', 'g','m']
    for i, y_val in enumerate(y_values):
        mask = (np.isclose(Y, y_val))
        x_vals = X[mask]
        z1_vals = np.abs(Z1_tot[mask])
        z2_vals = np.abs(Z2_tot[mask])
        sort_idx = np.argsort(x_vals)
        plt.plot(x_vals[sort_idx], z1_vals[sort_idx], color=colors[i], marker='o', markersize=3, linestyle='-', label=f"{meshes[0]} - y = {y_val:.2f}")
        plt.plot(x_vals[sort_idx], z2_vals[sort_idx], color=colors[i], marker='x', linestyle=':', label=f"{meshes[1]} - y = {y_val:.2f}")

    plt.axvline(x=dimX, color='black', linestyle='-', linewidth=1)
    plt.axvline(x=-dimX, color='black', linestyle='-', linewidth=1)
    plt.xlabel("x")
    plt.ylabel("|η|")
    plt.title(f"{titre}")
    plt.grid(True)
    plt.legend(loc='upper left')
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f"{titreFile}.png", dpi=300)
    plt.close()
    
def PLOT_2D(X, Y, Z, triangles, titre, titreFILE, bodyDIM) :
    fig, ax = plt.subplots(figsize=(8, 6))
    tri = Triangulation(X, Y, triangles)
    tcf = ax.tricontourf(tri, Z, levels=200, cmap="plasma")
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="10%", pad=0.5)
    cbar = fig.colorbar(tcf, cax=cax, orientation='horizontal')
    cbar.set_label(f"Free Surface Amplitude")
    cax.xaxis.set_ticks_position('bottom')
    cax.xaxis.set_label_position('bottom')
    cbar.ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    lx=bodyDIM[0]
    ly=bodyDIM[1]
    x_cyl = [-lx, lx, lx, -lx, -lx]
    y_cyl = [ly, ly, -ly, -ly, ly]
    ax.plot(x_cyl, y_cyl, color="black", linestyle='-', linewidth=1, label="FB")
    ax.set_aspect("equal")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title(titre)
    # plt.legend()
    plt.tight_layout()
    plt.savefig(titreFILE, dpi=300)
    plt.close()

def plot_mean(means_eta, omega, nom, title, titleFILE, meshes):
    # === TRACÉ FINAL : <|η|> = f(w) ===
    if nom == "w":
        unit = "[rad/s]"
        value = omega
        sort_idx = None
    else:
        unit = "[s]"
        T = 2 * np.pi / omega
        sort_idx = np.argsort(T)
        value = T[sort_idx]
    styles = ['o-', 'x:']  # premier champ, deuxième champ
    colors = ['b', 'r', 'g','m']
    plt.figure(figsize=(9, 6))
    for i, y_val in enumerate(means_eta[0].keys()): 
        for j, mean_dict in enumerate(means_eta):
            mean_list = mean_dict[y_val]
            mean_sorted = np.array(mean_list)[sort_idx] if sort_idx is not None else mean_list
            plt.plot(value, mean_sorted, styles[j], color=colors[i], label=f"{meshes[j]} - y={y_val:.1f}")

    plt.xlabel(f"{nom} {unit}")
    plt.ylabel("<|η|> (on y)")
    plt.title(f"{title}")
    plt.grid(True)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(f"{titleFILE}.png", dpi=300)
    plt.close()
    
            
def reduce_y_positions(y_positions, Ny):
    y_unique = np.sort(np.unique(y_positions))
    if Ny < 1 or Ny % 2 == 0:
        raise ValueError("Ny doit être impair et >=1")
    if Ny==1:
        y_final=[0.0]
    else :
        idxs = np.linspace(0, len(y_unique)-1, Ny, dtype=int)
        y_selected = y_unique[idxs].tolist()
        pos_selected = [y for y in y_selected if y > 0]
        y_final = [-y for y in reversed(pos_selected)] + ([0.0] if 0 in y_selected else []) + pos_selected
    
    print("y_final:", y_final)
    return np.array(sorted(y_final))

def recup_RAO(file, Nw, vect_dof):
    fileRAO = f"{file}/Motion/RAO.dat"
    RAOvalues = np.loadtxt(fileRAO, skiprows=2, max_rows=Nw)
    all_RAO = {
        1: RAOvalues[:, 1] * np.exp(1j * np.deg2rad(RAOvalues[:, 7])),   # surge
        2: RAOvalues[:, 2] * np.exp(1j * np.deg2rad(RAOvalues[:, 8])),   # 
        3: RAOvalues[:, 3] * np.exp(1j * np.deg2rad(RAOvalues[:, 9])),   # heave
        4: np.deg2rad(RAOvalues[:, 4]) * np.exp(1j * np.deg2rad(RAOvalues[:, 10])),  # 
        5: np.deg2rad(RAOvalues[:, 5]) * np.exp(1j * np.deg2rad(RAOvalues[:, 11])),  # pitch
        6: np.deg2rad(RAOvalues[:, 6]) * np.exp(1j * np.deg2rad(RAOvalues[:, 12]))  # pitch
    }
    RAO_complex = np.column_stack([all_RAO[d] for d in vect_dof])

    return RAO_complex
    # surge = RAOvalues[:, 1] * np.exp(1j * np.deg2rad(RAOvalues[:, 7]))
    # heave = RAOvalues[:, 3] * np.exp(1j * np.deg2rad(RAOvalues[:, 9]))
    # pitch =  np.deg2rad(RAOvalues[:, 5]) * np.exp(1j * np.deg2rad(RAOvalues[:, 11]))
    # RAO_complex = np.column_stack((surge, heave, pitch))

def recup_wave_number(file1, file2):
    fileWaveNumber1 = f"{file1}/WaveNumber.dat"
    fileWaveNumber2 = f"{file2}/WaveNumber.dat"
    data1 = np.loadtxt(fileWaveNumber1)
    data2 = np.loadtxt(fileWaveNumber2)
    if not np.array_equal(data1, data2):
        print("error : different data")
    omega = data1[:,0]
    number_k = data1[:,1]
    return omega, number_k

def recup_free_surface(file, num_pb):
    fichier_BEM = f"{file}/results/freesurface.{num_pb}.dat"
    points, triangles = READ_BEM(fichier_BEM)
    X = points[:,0]
    Y = points[:,1]
    Zreel = points[:,5]
    Zimag = points[:,6]
    return X, Y, (Zreel + 1j*Zimag), triangles
       
# === PARAMÈTRES ===
chemin="/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/Matlab/BEM_matlab/myProjets/"

mesh="Pontoon"
mesh2=f"{mesh}Wings"
meshes_name = [mesh, mesh2]
BEM_file1=f"{chemin}BEM_{mesh}_T"
BEM_file2=f"{chemin}BEM_{mesh2}_T"

Plot_1w = True
Plot_D = False
One_DOF=False
RAO=True
num_w=[1,3, 5]  
vect_dof=[1, 3, 5]
Ny_reduc = 1
freqType="T"
Nbeta=1
beta=0

folderPlot = f"{chemin}/FS_2Dplots/"
folderCOMP = f"{chemin}/COMP_FS/"

if not os.path.exists(folderPlot):
    os.makedirs(folderPlot)
    os.makedirs(folderCOMP)
omega, number_k = recup_wave_number(BEM_file1, BEM_file2)
Nw = len(omega)
Np=Nw*(6+Nbeta)

RAO_complex1 = recup_RAO(BEM_file1, Nw, vect_dof)
RAO_complex2 = recup_RAO(BEM_file2, Nw, vect_dof)
print("Ndof=", len(RAO_complex1[0,:]))

if mesh=='FB':
    bodyDIM=[22, 60]
elif mesh=='FB1':
    bodyDIM=[10, 100]
elif mesh=='Pontoon' or mesh=='PontoonWings' :
    bodyDIM=[12.5, 100]
if RAO :
    Scale="ScaleRAO_"
else:
    Scale=""
    
TypeResults = ["Total", "DR", "Diffraction"]
# TypeResults = ["Total", "Rpitch", "Diffraction"]
Titres = ["Incident + Diffracted + Radiated (surge, heave, pitch)", "Diffracted + Radiated (surge, heave, pitch)", "Incident + Diffracted"] 
# Titres = ["Incident + Diffracted + Radiated (surge, heave, pitch)", "Radiated (pitch)", "Diffracted"] 

fields = ["Z1_tot", "Z1_D", "Z1_R", "Z2_tot", "Z2_D", "Z2_R"]
mean_before = {}   # structure : mean_before[y][field] = [...]
mean_after  = {}
y_positions = None 

# === BOUCLE SUR LES FRÉQUENCES ===
print("-> DOF", vect_dof)
for iw in range(Nw):
    k = number_k[iw]
    w = omega[iw]
    print(f"--> Calcul fréquence {iw+1}/{Nw} : w = {w:.3f}, k = {k:.4f}, T={2*np.pi/w:.2f}")

    # bornes des fichiers freesurface
    start = iw * int(Np/Nw) + 1
    end   = (iw+1) * int(Np/Nw) + 1
    PROBLEMS = [f"{i:05d}" for i in range(start, end)]

    Z1_tot = 0
    Z1_D = 0
    Z1_R = 0
    Z2_tot = 0
    Z2_D = 0
    Z2_R = 0
    for num_pb in PROBLEMS:
        X, Y, Z1, triangles = recup_free_surface(BEM_file1, num_pb)
        X, Y, Z2, triangles = recup_free_surface(BEM_file2, num_pb)
        
        type_pb = (int(num_pb) - 1) % (Np // Nw)
        if  type_pb == 0:
            if Plot_D and iw+1 in num_w :
                print("Problem n°", num_pb, "-> Diffraction")
            Z1_D = Z1
            Z2_D = Z2
            i_dof=0
        elif type_pb in vect_dof:
            if One_DOF and iw+1 in num_w :
                print("Problem n°", num_pb, "-> Radiation - DOF", vect_dof)
            if RAO :
                RAO_dof_1 = RAO_complex1[iw,i_dof]
                RAO_dof_2 = RAO_complex2[iw,i_dof]
                print("RAO", i_dof, RAO_dof_1, RAO_dof_2)
                Z1_R += Z1*RAO_dof_1
                Z2_R += Z2*RAO_dof_2
                i_dof+=1
            else :
                Z1_R += Z1
                Z2_R += Z2

    # Ajout de l'onde incidente
    Z_I = np.exp(1j * k * (X*np.cos(beta) + Y*np.sin(beta)))
    Z1_tot = Z_I +Z1_D + Z1_R
    Z1_D += Z_I
    Z1_R += Z1_D
    Z2_tot = Z_I +Z2_D + Z2_R
    Z2_D += Z_I
    Z2_R += Z2_D

    # Déterminer les x distincts la première fois
    if y_positions is None:
        y_pos = np.unique(Y)
        y_pos_in = y_pos[(y_pos > -bodyDIM[1]) & (y_pos < bodyDIM[1])]
        y_positions = reduce_y_positions(y_pos_in, Ny_reduc)
        
    
    # Si c’est la fréquence choisie → tracer η(y)
    if iw+1 in num_w and Plot_1w :
        print("plot")
        plot_eta_vs_x(X, Y, Z1_tot, Z2_tot, y_positions, f"Free Surface : {Titres[0]} - T = {2*np.pi/w:.1f} s", f"{folderCOMP}/ComparisonX_{TypeResults[0]}_{Scale}Ny{Ny_reduc}_T{2*np.pi/w:.1f}", bodyDIM[0], meshes_name)
        mask_out = (X < -bodyDIM[0]) | (X > bodyDIM[0]) | (Y < -bodyDIM[1]) | (Y > bodyDIM[1])
        Zout1 = np.where(mask_out, Z1_tot, 0.0)
        PLOT_2D(X, Y, np.abs(Zout1), triangles, f"{mesh} - Free Surface : {Titres[0]} - T={2*np.pi/w:.1f}s", f"{folderPlot}/FS_{Scale}2d_{mesh}_T{2*np.pi/w:.1f}.png", bodyDIM)
        Zout2 = np.where(mask_out, Z2_tot, 0.0)
        PLOT_2D(X, Y, np.abs(Zout2), triangles, f"{mesh2} - Free Surface : {Titres[0]} - T={2*np.pi/w:.1f}s", f"{folderPlot}/FS_{Scale}2d_{mesh2}_T{2*np.pi/w:.1f}.png", bodyDIM)
        if Plot_D :
            plot_eta_vs_x(X, Y, Z1_D, Z2_D, y_positions, f"Free Surface : {Titres[2]} - T = {2*np.pi/w:.1f} s", f"{folderCOMP}/ComparisonX_{TypeResults[2]}_Ny{Ny_reduc}_T{2*np.pi/w:.1f}", bodyDIM[0], meshes_name)
            ZoutD1 = np.where(mask_out, Z1_D, 0.0)
            PLOT_2D(X, Y, np.abs(ZoutD1), triangles, f"{mesh} - Free Surface : {Titres[2]} - T={2*np.pi/w:.1f}s", f"{folderPlot}/FS_{TypeResults[2]}_2d_{mesh}_T{2*np.pi/w:.1f}.png", bodyDIM)
            ZoutD2 = np.where(mask_out, Z2_D, 0.0)
            PLOT_2D(X, Y, np.abs(ZoutD2), triangles, f"{mesh2} - Free Surface : {Titres[2]} - T={2*np.pi/w:.1f}s", f"{folderPlot}/FS_{TypeResults[2]}_2d_{mesh2}_T{2*np.pi/w:.1f}.png", bodyDIM)
        if One_DOF : 
            plot_eta_vs_x(X, Y, Z1_R, Z2_R, y_positions, f"Free Surface : {Titres[1]} - T = {2*np.pi/w:.1f} s", f"{folderCOMP}/ComparisonX_{TypeResults[1]}_{Scale}Ny{Ny_reduc}_T{2*np.pi/w:.1f}", bodyDIM[0], meshes_name)
            ZoutR1 = np.where(mask_out, Z1_R, 0.0)
            PLOT_2D(X, Y, np.abs(ZoutR1), triangles, f"{mesh} - Free Surface : {Titres[1]} - T={2*np.pi/w:.1f}s", f"{folderPlot}/FS_{Scale}{TypeResults[1]}_2d_{mesh}_T{2*np.pi/w:.1f}.png", bodyDIM)
            ZoutR2 = np.where(mask_out, Z2_R, 0.0)
            PLOT_2D(X, Y, np.abs(ZoutR2), triangles, f"{mesh2} - Free Surface : {Titres[1]} - T={2*np.pi/w:.1f}s", f"{folderPlot}/FS_{Scale}{TypeResults[1]}_2d_{mesh2}_T{2*np.pi/w:.1f}.png", bodyDIM)

    FS_types = { "Z1_tot": Z1_tot, "Z1_D": Z1_D, "Z1_R": Z1_R,
        "Z2_tot": Z2_tot, "Z2_D": Z2_D, "Z2_R": Z2_R}
    for y_val in y_positions:
        if y_val not in mean_before:
            mean_before[y_val] = {name: [] for name in fields}
            mean_after[y_val]  = {name: [] for name in fields}

        mask_before = np.isclose(Y, y_val) & (X < -bodyDIM[0])
        mask_after  = np.isclose(Y, y_val) & (X >  bodyDIM[0])

        for name, value in FS_types.items():
            mean_before[y_val][name].append(np.mean(np.abs(value[mask_before])))
            mean_after[y_val][name].append(np.mean(np.abs(value[mask_after])))

if not Plot_1w :
    meansB_tot=[{y:mean_before[y]["Z1_tot"] for y in mean_before}, {y:mean_before[y]["Z2_tot"] for y in mean_before} ]
    plot_mean(meansB_tot,omega,freqType, f"Before body - Mean Free Surface : {Titres[0]}", f"{folderCOMP}/ComparisonMean_Before_{TypeResults[0]}_{Scale}Ny{Ny_reduc}", meshes_name)
    meansA_tot=[{y:mean_after[y]["Z1_tot"] for y in mean_after}, {y:mean_after[y]["Z2_tot"] for y in mean_after} ]
    plot_mean(meansA_tot,omega,freqType, f"After body - Mean Free Surface : {Titres[0]}", f"{folderCOMP}/ComparisonMean_After_{TypeResults[0]}_{Scale}Ny{Ny_reduc}", meshes_name)
    if Plot_D :
        meansB_D=[{y:mean_before[y]["Z1_D"] for y in mean_before}, {y:mean_before[y]["Z2_D"] for y in mean_before} ]
        plot_mean(meansB_D,omega,freqType, f"Before body - Mean Free Surface : {Titres[2]}", f"{folderCOMP}/ComparisonMean_Before_{TypeResults[2]}_Ny{Ny_reduc}", meshes_name)
        meansA_D=[{y:mean_after[y]["Z1_D"] for y in mean_after}, {y:mean_after[y]["Z2_D"] for y in mean_after} ]
        plot_mean(meansA_D,omega,freqType, f"After body - Mean Free Surface : {Titres[2]}", f"{folderCOMP}/ComparisonMean_After_{TypeResults[2]}_Ny{Ny_reduc}", meshes_name)
    if One_DOF :
        meansB_R=[{y:mean_before[y]["Z1_R"] for y in mean_before}, {y:mean_before[y]["Z2_R"] for y in mean_before} ]
        plot_mean(meansB_R,omega,freqType, f"Before body - Mean Free Surface : {Titres[1]}", f"{folderCOMP}/ComparisonMean_Before_{TypeResults[1]}_{Scale}Ny{Ny_reduc}", meshes_name)
        meansA_R=[{y:mean_after[y]["Z1_R"] for y in mean_after}, {y:mean_after[y]["Z2_R"] for y in mean_after} ]
        plot_mean(meansA_R,omega,freqType, f"After body - Mean Free Surface : {Titres[1]}", f"{folderCOMP}/ComparisonMean_After_{TypeResults[1]}_{Scale}Ny{Ny_reduc}", meshes_name)
    


