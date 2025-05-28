import numpy as np
import re  # Pour utiliser les expressions régulières
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import itertools



def plot_coef(file, nom, titre, mesh, Nw, Ne, indices, AR) : 
    data = np.loadtxt(file, delimiter=None, skiprows=0) 
    colors = itertools.cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    omegas = data[:, 0]
    Nbeta=11
    ka=True
    if ka :
        abscisse=np.loadtxt(f"Data_progressive_waves_{mesh}.dat", skiprows=1)
        wave_numbers=abscisse[:,1]
    
    plt.figure(figsize=(8, 6))
    for k, ij in enumerate(indices): 
        coef=np.zeros(Nw)
        w=np.zeros(Nw)
        for i in range(Nw):
            # ! signe pour la convention -> comparaison avec résultats WAMIT - McNatt
            coef[i]=data[ij[0]+i*Nbeta,ij[1]]
            w[i]=omegas[i*Nbeta]
        mode1=ij[0] - int((Nbeta-1)/2)
        if AR :
            mode2=ij[1]
        else :
            mode2=ij[1] - int((Nbeta-1)/2)
        if ka:
            w=wave_numbers*5
        if nom=="Imag(D)":
            coef=-coef
        plt.plot(w, coef, marker='o', color=next(colors), linestyle='-', label=f'({mode1},{mode2})')
    
    if ka :
        plt.xlabel("ka (with a=5m)")
    else :
        plt.xlabel("frequencies (rad/s)")
    plt.ylabel(f"{nom}")
    if AR :
        plt.title(f"{nom} dof n°{ij[1]} \n ({titre})")
        title=f"Graphe_{mesh}_{nom}_{ij[1]}.png"
    else:
        plt.title(f"{nom} \n ({titre})")
        title=f"Graphe_{mesh}_{nom}.png"
    plt.legend()
    plt.grid()        
    plt.savefig(title)
    plt.close()  # Ferme la figure pour éviter trop d'ouvertures


chemin="/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/PIT3_E/wec_inputs/"

PLOT_D=True
PLOT_AR=True
mesh="Cylinder"
Nw=10
indicesD=[[5, 5], [6, 6], [4, 4], [7, 7], [3,3], [8,8], [9,9]]
# indicesD=[[6, 6], [7, 7], [8,8], [9,9], [10, 10]]
DOF=[1, 2]
Ne=0
type=""
titre=f"{mesh} - Nw={Nw}, Nbeta=11, dof=2, Ndir=1"
PIT= f"{chemin}PIT3_McNatt_{mesh}/results/OneBodyProblem{type}"

matrix_D_abs=f"{PIT}_DiffractionMatrix_abs.dat"
matrix_D_ph=f"{PIT}_DiffractionMatrix_ph.dat"
matrix_D_Re=f"{PIT}_DiffractionMatrix_Re.dat"
matrix_D_Imag=f"{PIT}_DiffractionMatrix_Imag.dat"

coefs_AR_abs=f"{PIT}_RadiationCoefficients_abs.dat"
coefs_AR_ph=f"{PIT}_RadiationCoefficients_ph.dat"
coefs_AR_Re=f"{PIT}_RadiationCoefficients_Re.dat"
coefs_AR_Imag=f"{PIT}_RadiationCoefficients_Imag.dat"

Nbeta=11
# IndMode=np.linspace(1, Nbeta, Nbeta) 
# Mode=IndMode -1 - (Nbeta-1)/2
# print(IndMode)
# print(Mode)
# Mode=[[0, 0], [1, 1]]
# indices=Mode+(Nbeta-1)/2+1
# print(indices)
if PLOT_D:
    plot_coef(matrix_D_Re, "Re(D)", titre, mesh, Nw, Ne, indicesD, False)
    plot_coef(matrix_D_Imag, "Imag(D)", titre, mesh, Nw, Ne, indicesD, False)

if PLOT_AR:
    for c, i in enumerate(DOF) : 
        indicesAR=[[5, i], [6, i], [4, i], [8, i], [2, i]]
        plot_coef(coefs_AR_Re, "Re(aR)", titre, mesh, Nw, Ne, indicesAR, True)
        plot_coef(coefs_AR_Imag, "Imag(aR)", titre, mesh, Nw, Ne, indicesAR, True)
