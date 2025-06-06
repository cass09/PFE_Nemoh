import numpy as np
import cmath
import os
import h5py
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def PlotData(WECArr, farm, distance, directionMB, vect_dof, PlotFolder):

    w = 2*np.pi/WECArr.period
    Nb=farm['N_bodies']
    if len(directionMB)==1 :
        beta_value=f"beta={directionMB[0]}°"
    else :
        beta_value=""
    if farm['Configuration'] :
        type=farm["Type"]
    else :
        type=f"Nb={Nb}"
    if farm['Ne_modes']>0 : 
        test= f"Ne={farm['Ne_modes']}"
    else : 
        test="Ne=0"
    
    titre=f"{type} - Nw={len(w)} - {beta_value} - {test} - distance={distance}m"
    
    for i, dof in enumerate(vect_dof):
        num_dof=[]
        indice=[]
        for j in range(Nb):
            indice.append([dof+j*6,dof+j*6])
            num_dof.append(dof+j*6)
            for k in range(j+1, Nb):
                indice.append([dof+j*6,dof+k*6])
        Trace_RAD(w, WECArr.Madd, "Added_Mass", indice, titre, distance, PlotFolder)
        Trace_RAD(w, WECArr.Crad, "Damping", indice, titre, distance, PlotFolder)
        Trace_Fex(w, WECArr.Fex, "Fex", num_dof, titre, distance, PlotFolder)
        if farm['RAO'] :
            Trace_Fex(w, WECArr.RAO, "RAO", num_dof, titre, distance, PlotFolder)
    return True

def Trace_RAD(freq, data, nom, indice, titre, d, PlotFolder):
   
    couleurs = ['b', 'g', 'r', 'gold', 'lime', 'c', 'm', 'y', 'k']  # palette de couleurs (réutilisée si plus de 7 courbes)
    markers = ['s', 'o', '<', '^','x', '+', '*']
    
    plt.figure(figsize=(8, 6))
    plt.title(titre)
    plt.xlabel("Frequency (rad/s)")
    if nom=="Damping":
        unit="(kg/s)"
    else :
        unit="(kg)"
    plt.ylabel(f"{nom} {unit}")
    plt.grid(True)

    label_indice=''
    for k, ij in enumerate(indice): 
        i=ij[0]-1
        j=ij[1]-1
        color = couleurs[k % len(couleurs)]
        marker = markers[k % len(markers)]
        plt.plot(freq, data[:,i,j], linestyle='--', marker=marker, color=color, markersize=5, label=f"{ij[0]}-{ij[1]}")
        label_indice=f"{label_indice}_{ij[0]}-{ij[1]}"
    plt.legend(loc='best', fontsize='small', frameon=True)
    filename = os.path.join(PlotFolder, f"ITM_{nom}_d{d}_M{label_indice}.png")
    plt.savefig(filename, dpi=300)
    plt.close()
    return

def Trace_Fex(freq, data, nom, num_dof, titre, d, PlotFolder):
    couleurs = ['b', 'g', 'r', 'gold', 'lime', 'b', 'c', 'm', 'y', 'k']  # palette de couleurs (réutilisée si plus de 7 courbes)
    markers = ['s', 'o', '<', '^','x', '+', '*']
    
    # ONE wave direction !!
    i=0
    plt.figure(figsize=(10, 8))
    plt.title(titre)
    plt.xlabel("Frequency (rad/s)")
    plt.ylabel(f"|{nom}|")
    nom_dof=""
    for k, dof in enumerate(num_dof) : 
        j=dof-1
        nom_dof=nom_dof+"_"+f"{dof}"
        color = couleurs[k % len(couleurs)]
        marker = markers[k % len(markers)]
        plt.plot(freq, np.abs(data[:, i, j]), linestyle='--', marker=marker, color=color, markersize=4, label=f"dof {j+1}")   
    plt.legend()
    plt.grid(True)
    filename = os.path.join(PlotFolder, f"ITM_{nom}_abs_d{d}_dof{nom_dof}.png")
    plt.savefig(filename, dpi=300)
    plt.close()

    plt.figure(figsize=(10, 8))
    plt.title(titre)
    plt.xlabel("Frequency (rad/s)")
    plt.ylabel(f"phase({nom})")
    for k, dof in enumerate(num_dof) : 
        j=dof-1
        color = couleurs[k % len(couleurs)]
        marker = markers[k % len(markers)]
        plt.plot(freq, np.angle(data[:, i, j]), linestyle='--', marker=marker, color=color, markersize=4, label=f"dof {j+1}")   
    plt.legend()
    plt.grid(True)
    filename = os.path.join(PlotFolder, f"ITM_{nom}_phase_d{d}_dof{nom_dof}.png")
    plt.savefig(filename, dpi=300)
    plt.close()
    return
