import numpy as np
import re  # Pour utiliser les expressions régulières
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import itertools



def plot_coef(file, fileS, fileC, nom, titre, mesh, Nw, Nbeta, indices, AR) : 
    data = np.loadtxt(file, delimiter=None, skiprows=0) 
    dataS = np.loadtxt(fileS, delimiter=None, skiprows=0) 
    # dataC = np.loadtxt(fileC, delimiter=None, skiprows=0) 
    colors = itertools.cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    omegas = data[:, 0]
    ka=False
    mat=data[:,1:]
    matS=dataS[:,1:]
    # matC=dataC[:,1:]
    # print(matS[:Nbeta+1,:])
    if ka :
        abscisse=np.loadtxt(f"Data_progressive_waves_{mesh}.dat", skiprows=1)
        wave_numbers=abscisse[:,1]
    
    plt.figure(figsize=(8, 6))
    for k, ij in enumerate(indices): 
        coef=np.zeros(Nw)
        coefS=np.zeros(Nw)
        # coefC=np.zeros(Nw)
        w=np.zeros(Nw)
        M=int((Nbeta-1)/2)
        mode1=ij[0] - M
        if AR :
            mode2=ij[1]+1
        else :
            mode2=ij[1] - M
        print(nom)
        print("w[i]", "coef[i]", "coefS[i]", "coefC[i]")
        for i in range(Nw):
            # ! signe pour la convention -> comparaison avec résultats WAMIT - McNatt
            coef[i]=mat[ij[0]+i*Nbeta,ij[1]]
            coefS[i]=matS[ij[0]+i*Nbeta,ij[1]]
            # coefC[i]=matC[ij[0]+i*Nbeta,ij[1]]
            	#coef[i]=data[-ij[0]+2*M+i*Nbeta,-ij[1]+2*M]
            w[i]=omegas[i*Nbeta]
            
            # print(w[i], coef[i], coefS[i], coefC[i])
        #if nom=="D_Imag":
            #coef=-coef
        if ka:
            w=wave_numbers*3
        
        plt.plot(w, coef, marker='o', color="black", linestyle='-', label=f'Outer Cylinder Method with Nemoh ({mode1},{mode2})')    
        plt.plot(w, coefS, marker='x', color="b", linestyle=':', label=f'Source Terms Method with Nemoh ({mode1},{mode2})')    
        # plt.plot(w, coefC, marker='+', color="m", linestyle=':', label=f'Source Terms Method with Capytaine({mode1},{mode2})')    
        if ka :
            plt.xlabel("ka (with a=3m)")
        else :
            plt.xlabel("frequencies (rad/s)")
        plt.ylabel(f"{nom}")
        if AR :
            if nom=="aR_Imag":
               data="Radiation Coefficient - Imaginary part"
            elif nom=="aR_Re":
               data="Radiation Coefficient - Real part"
            elif nom=="aR_abs":
               data="Radiation Coefficient - Absolute value"
            elif nom=="aR_ph":
               data="Radiation Coefficient - Phase"
            elif nom=="G_Imag":
               data="Force Transfers Matrix - Imaginary part"
            elif nom=="G_Re":
               data="Force Transfers Matrix - Real part"
            plt.title(f"{data} - DOF n°{ij[1]+1} - {titre}")
            title=f"Graphe_m0_{mesh}_{nom}_{mode1}_{ij[1]+1}.png"
        else:
            if nom=="D_Imag":
               data="Diffraction Matrix - Imaginary part"
            elif nom=="D_Re":
               data="Diffraction Matrix - Real part"
            plt.title(f"{data} - {titre}")
            title=f"Graphe_m0_{mesh}_{nom}_{mode1}.png"
        plt.legend()
        plt.grid()        
        plt.savefig(title)
        plt.close()  # Ferme la figure pour éviter trop d'ouvertures

def plot_coefS(file, fileS, fileC, nom, titre, mesh, Nw, Nbeta, indices, AR) : 
    data = np.loadtxt(file, delimiter=None, skiprows=0) 
    dataS = np.loadtxt(fileS, delimiter=None, skiprows=0) 
    dataC = np.loadtxt(fileC, delimiter=None, skiprows=0) 
    colors = itertools.cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    markers=["s", "o", "+", "o", "+", "o", "+"]
    markers2=["s", "^", "x", "^", "x", "^", "x"]
    omegas = data[:, 0]
    ka=False
    mat=data[:,1:]
    matS=dataS[:,1:]
    matC=dataC[:,1:]
    rho=1000
    g=9.81
    a=3
    # AD=rho*g*np.pi*2*a**2
    AD=1
    # print(matS[:Nbeta+1,:])
    if ka :
        abscisse=np.loadtxt(f"Data_progressive_waves_{mesh}.dat", skiprows=1)
        wave_numbers=abscisse[:,1]
    
    plt.figure(figsize=(8, 6))
    for k, ij in enumerate(indices): 
        coef=np.zeros(Nw)
        coefS=np.zeros(Nw)
        coefC=np.zeros(Nw)
        w=np.zeros(Nw)
        M=int((Nbeta-1)/2)
        mode1=ij[0] - M
        if AR :
            mode2=ij[1]+1
        else :
            mode2=ij[1] - M
        print(nom)
        print("w[i]", "coef[i]", "coefS[i]", "coefC[i]")
        for i in range(Nw):
            # ! signe pour la convention -> comparaison avec résultats WAMIT - McNatt
            coef[i]=mat[ij[0]+i*Nbeta,ij[1]]/AD
            coefS[i]=matS[ij[0]+i*Nbeta,ij[1]]/AD
            coefC[i]=matC[ij[0]+i*Nbeta,ij[1]]/AD
            	#coef[i]=data[-ij[0]+2*M+i*Nbeta,-ij[1]+2*M]
            w[i]=omegas[i*Nbeta]
            # if AR : 
            #     coefS[i]=coefS[i]*np.sqrt(w[i])
            # if not AR :
            #     coefS[i]=coefS[i]*np.sqrt(w[i])
            print(w[i], coef[i], coefS[i], coefC[i])
        #if nom=="D_Imag":
            #coef=-coef
        if ka:
            w=wave_numbers*a
        color_mode=next(colors)
        plt.plot(w, coef, marker='o', markersize=2, color=color_mode, linestyle='-', label=f'OCM with Nemoh ({mode1},{mode2})')    
        plt.plot(w, coefS, marker='^', markersize=4, color=color_mode, linestyle=':', label=f'STM with Nemoh ({mode1},{mode2})')    
        plt.plot(w, coefC, marker='+', markersize=8, color=color_mode, linestyle=':', label=f'STM with Capytaine({mode1},{mode2})')    
    if ka :
        plt.xlabel("ka (with a=3m)")
    else :
        plt.xlabel("frequencies (rad/s)")
    if AD>1:
        plt.ylabel(f"{nom}/(rho*g*pi*2a^2)")
    else : 
        plt.ylabel(f"{nom}")
    if AR :
    	if nom=="aR_Imag":
    	    data="Radiation Coefficient - Imaginary part"
    	elif nom=="aR_Re":
    	    data="Radiation Coefficient - Real part"
    	elif nom=="G_Imag":
    	    data="Force Transfers Matrix - Imaginary part"
    	elif nom=="G_Re":
    	    data="Force Transfers Matrix - Real part"
    	plt.title(f"{data} - DOF n°{ij[1]+1} - {titre}")
    	title=f"Graphe_{mesh}_{nom}_DOF{ij[1]+1}.png"
    else:
    	if nom=="D_Imag":
    	    data="Diffraction Matrix - Imaginary part"
    	elif nom=="D_Re":
    	    data="Diffraction Matrix - Real part"
    	plt.title(f"{data} - {titre}")
    	title=f"Graphe_{mesh}_{nom}.png"
    plt.legend()
    plt.grid()        
    plt.savefig(title)
    plt.close()  # Ferme la figure pour éviter trop d'ouvertures
    
def plot_coefS_mesh(file, file2, nom, titre, mesh, Nw, Nbeta, indices, AR) : 
    data = np.loadtxt(file, delimiter=None, skiprows=0) 
    data2 = np.loadtxt(file2, delimiter=None, skiprows=0) 
    colors = itertools.cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    markers=["s", "o", "+", "o", "+", "o", "+"]
    markers2=["s", "^", "x", "^", "x", "^", "x"]
    omegas = data[:, 0]
    ka=True
    mat=data[:,1:]
    mat2=data2[:,1:]
    if ka :
        abscisse=np.loadtxt(f"Data_progressive_waves_{mesh[0]}.dat", skiprows=1)
        wave_numbers=abscisse[:,1]
        abscisse2=np.loadtxt(f"Data_progressive_waves_{mesh[1]}.dat", skiprows=1)
        wave_numbers2=abscisse2[:,1]
    
    plt.figure(figsize=(8, 6))
    for k, ij in enumerate(indices): 
        coef=np.zeros(Nw)
        coef2=np.zeros(Nw)
        w=np.zeros(Nw)
        M=int((Nbeta-1)/2)
        mode1=ij[0] - M
        if AR :
            mode2=ij[1]+1
        else :
            mode2=ij[1] - M
        print(nom)
        print("w[i]", "coef[i]", "coefS[i]", "coefC[i]")
        for i in range(Nw):
            # ! signe pour la convention -> comparaison avec résultats WAMIT - McNatt
            coef[i]=mat[ij[0]+i*Nbeta,ij[1]]
            coef2[i]=mat2[ij[0]+i*Nbeta,ij[1]]
            w[i]=omegas[i*Nbeta]
            # print(w[i], coef[i], coefS[i], coefC[i])
        if ka:
            a=3
            w=wave_numbers*a
            w2=wave_numbers2*a
        
        plt.plot(w, coef, marker=markers[k], color=next(colors), linestyle='--', label=f'{mesh[0]} ({mode1},{mode2})')    
        plt.plot(w2, coef2, marker=markers2[k], color=next(colors), linestyle=':', label=f'{mesh[1]} ({mode1},{mode2})')     
    if ka :
        plt.xlabel("ka (with a=3m)")
    else :
        plt.xlabel("frequencies (rad/s)")
    plt.ylabel(f"{nom}")
    if AR :
        plt.title(f"Outer Cylinder Method - {nom} - dof n°{ij[1]+1} \n ({titre})")
        title=f"Graphe_{mesh}_{nom}_{ij[1]+1}.png"
    else:
        plt.title(f"Outer Cylinder Method - {nom} \n ({titre})")
        title=f"Graphe_{mesh}_{nom}.png"
    plt.legend()
    plt.grid()        
    plt.savefig(title)
    plt.close()  # Ferme la figure pour éviter trop d'ouvertures
    
    
def plot_G_fex(file_Re, file_Im, fileS_Re, fileS_Im, fileC_Re, fileC_Im, titre, mesh, Nw, Nbeta, indices) : 
    data_R = np.loadtxt(file_Re, delimiter=None, skiprows=0) 
    data_I = np.loadtxt(file_Im, delimiter=None, skiprows=0) 
    dataS_R = np.loadtxt(fileS_Re, delimiter=None, skiprows=0) 
    dataS_I = np.loadtxt(fileS_Im, delimiter=None, skiprows=0) 
    dataC_R = np.loadtxt(fileC_Re, delimiter=None, skiprows=0) 
    dataC_I = np.loadtxt(fileC_Im, delimiter=None, skiprows=0) 
    colors = itertools.cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    omegas = dataS_R[:, 0]
    mat=data_R[:,1:]+1j*data_I[:,1:]
    matS=dataS_R[:,1:]+1j*dataS_I[:,1:]
    matC=dataC_R[:,1:]+1j*dataC_I[:,1:]
    # print(matS[:Nbeta+1,:])
    for k, ij in enumerate(indices): 
        coef=np.zeros(Nw)
        coefS=np.zeros(Nw)
        coefC=np.zeros(Nw)
        w=np.zeros(Nw)
        M=int((Nbeta-1)/2)
        mode1=ij[0] - M
        mode2=ij[1]+1
        
        print("w[i]", "coefS[i]", "coefC[i]")
        for i in range(Nw):
            coef[i]=mat[ij[0]+i*Nbeta,ij[1]]
            coefS[i]=matS[ij[0]+i*Nbeta,ij[1]]
            coefC[i]=matC[ij[0]+i*Nbeta,ij[1]]
            w[i]=omegas[i*Nbeta]
            # if AR : 
            #     coef[i]=coef[i]/w[i]
            print(w[i], np.abs(coefS[i]), np.abs(coefC[i]))
        plt.figure(figsize=(8, 6))
        plt.plot(w, np.abs(coef), marker='o', color=next(colors), linestyle='--', label=f'Outer Cylinder Method  ({mode1},{mode2})')    
        plt.plot(w, np.abs(coefS), marker='^', color=next(colors), linestyle=':', label=f'Source Terms Method with Nemoh ({mode1},{mode2})')    
        plt.plot(w, np.abs(coefC), marker='+', color=next(colors), linestyle=':', label=f'Source Terms Method with Capytaine({mode1},{mode2})')    
        plt.xlabel("frequencies (rad/s)")
        plt.ylabel(f"G amplitude")
        plt.title(f"G abs dof n°{ij[1]+1} \n ({titre})")
        title=f"Graphe_{mesh}_G_abs_{mode1}_{ij[1]+1}.png"
        plt.legend()
        plt.grid()        
        plt.savefig(title)
        plt.close()  # Ferme la figure pour éviter trop d'ouvertures

        plt.figure(figsize=(8, 6))
        plt.plot(w, np.angle(coef), marker='o', color=next(colors), linestyle='--', label=f'Outer Cylinder Method  ({mode1},{mode2})')    
        plt.plot(w, np.angle(coefS), marker='^', color=next(colors), linestyle=':', label=f'Source Terms Method with Nemoh ({mode1},{mode2})')    
        plt.plot(w, np.angle(coefC), marker='+', color=next(colors), linestyle=':', label=f'Source Terms Method with Capytaine({mode1},{mode2})')    
        plt.xlabel("frequencies (rad/s)")
        plt.ylabel(f"G phase")
        plt.title(f"G phase dof n°{ij[1]+1} \n ({titre})")
        title=f"Graphe_{mesh}_G_ph_{mode1}_{ij[1]+1}.png"
        plt.legend()
        plt.grid()        
        plt.savefig(title)
        plt.close()  # Ferme la figure pour éviter trop d'ouvertures


chemin="/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/PIT3_E/wec_inputs/"
cheminC="/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/capytaine/MyTestCases/"

PLOT_D=False
PLOT_AR=True
PLOT_G=False
PLOT_G_fex=False
COMP_mesh=False
mesh="CylR3"
Nw=20
depth=""
# indicesD=[[6, 6], [5, 5], [7, 7], [4, 4], [8,8], [3,3], [9,9]]
# indicesD=[[5, 5], [4, 4], [6, 6]]
indicesD=[[4, 4], [5,5], [6,6]]
DOF=[3]
Nbeta=9

Ne=0
type=""
titre=f"Cylinder"
# PIT= f"{chemin}PIT3_{mesh}_h10/results/OneBodyProblem{type}"
PIT= f"{chemin}PIT3_{mesh}{depth}/results/OneBodyProblem{type}"
PIT_OC= f"{cheminC}CAP_{mesh}_OC/resultsBEM/OneBodyProblem{type}"
PIT_S= f"{chemin}PIT3_{mesh}{depth}_source/results/S_OneBodyProblem{type}"
PIT_C= f"{cheminC}CAP_{mesh}{depth}_source/resultsBEM/S_OneBodyProblem{type}"

matrix_D_abs=f"{PIT}_DiffractionMatrix_abs.dat"
matrix_D_ph=f"{PIT}_DiffractionMatrix_ph.dat"
matrix_D_Re=f"{PIT}_DiffractionMatrix_Re.dat"
matrix_D_Imag=f"{PIT}_DiffractionMatrix_Imag.dat"

matrix_G_Re=f"{PIT}_ForceTransferMatrix_Re.dat"
matrix_G_Imag=f"{PIT}_ForceTransferMatrix_Imag.dat"

coefs_AR_abs=f"{PIT}_RadiationCoefficients_abs.dat"
coefs_AR_ph=f"{PIT}_RadiationCoefficients_ph.dat"
coefs_AR_Re=f"{PIT}_RadiationCoefficients_Re.dat"
coefs_AR_Imag=f"{PIT}_RadiationCoefficients_Imag.dat"


matrix_D_abs_S=f"{PIT_S}_DiffractionMatrix_abs.dat"
matrix_D_ph_S=f"{PIT_S}_DiffractionMatrix_ph.dat"
matrix_D_Re_S=f"{PIT_S}_DiffractionMatrix_Re.dat"
matrix_D_Imag_S=f"{PIT_S}_DiffractionMatrix_Imag.dat"

matrix_G_Re_S=f"{PIT_S}_ForceTransferMatrix_Re.dat"
matrix_G_Imag_S=f"{PIT_S}_ForceTransferMatrix_Imag.dat"

coefs_AR_abs_S=f"{PIT_S}_RadiationCoefficients_abs.dat"
coefs_AR_ph_S=f"{PIT_S}_RadiationCoefficients_ph.dat"
coefs_AR_Re_S=f"{PIT_S}_RadiationCoefficients_Re.dat"
coefs_AR_Imag_S=f"{PIT_S}_RadiationCoefficients_Imag.dat"

matrix_D_abs_C=f"{PIT_C}_DiffractionMatrix_abs.dat"
matrix_D_ph_C=f"{PIT_C}_DiffractionMatrix_ph.dat"
matrix_D_Re_C=f"{PIT_C}_DiffractionMatrix_Re.dat"
matrix_D_Imag_C=f"{PIT_C}_DiffractionMatrix_Imag.dat"

matrix_G_Re_C=f"{PIT_C}_ForceTransferMatrix_Re.dat"
matrix_G_Imag_C=f"{PIT_C}_ForceTransferMatrix_Imag.dat"

coefs_AR_abs_C=f"{PIT_C}_RadiationCoefficients_abs.dat"
coefs_AR_ph_C=f"{PIT_C}_RadiationCoefficients_ph.dat"
coefs_AR_Re_C=f"{PIT_C}_RadiationCoefficients_Re.dat"
coefs_AR_Imag_C=f"{PIT_C}_RadiationCoefficients_Imag.dat"


matrix_D_abs_OC=f"{PIT_OC}_DiffractionMatrix_abs.dat"
matrix_D_ph_OC=f"{PIT_OC}_DiffractionMatrix_ph.dat"
matrix_D_Re_OC=f"{PIT_OC}_DiffractionMatrix_Re.dat"
matrix_D_Imag_OC=f"{PIT_OC}_DiffractionMatrix_Imag.dat"

matrix_G_Re_OC=f"{PIT_OC}_ForceTransferMatrix_Re.dat"
matrix_G_Imag_OC=f"{PIT_OC}_ForceTransferMatrix_Imag.dat"

coefs_AR_abs_OC=f"{PIT_OC}_RadiationCoefficients_abs.dat"
coefs_AR_ph_OC=f"{PIT_OC}_RadiationCoefficients_ph.dat"
coefs_AR_Re_OC=f"{PIT_OC}_RadiationCoefficients_Re.dat"
coefs_AR_Imag_OC=f"{PIT_OC}_RadiationCoefficients_Imag.dat"
# IndMode=np.linspace(1, Nbeta, Nbeta) 
# Mode=IndMode -1 - (Nbeta-1)/2
# print(IndMode)
# print(Mode)
# Mode=[[0, 0], [1, 1]]
# indices=Mode+(Nbeta-1)/2+1
# print(indices)
if PLOT_D:
    plot_coef(matrix_D_Re, matrix_D_Re_S, matrix_D_Re_C, "D_Re", titre, mesh, Nw, Nbeta, indicesD, False)
    plot_coef(matrix_D_Imag, matrix_D_Imag_S, matrix_D_Imag_C, "D_Imag", titre, mesh, Nw, Nbeta, indicesD, False)

for c, i in enumerate(DOF) : 
    indicesAR=[[7, i-1],[8, i-1], [0, i-1]]
    # indicesAR=[[4, i-1]]
    if PLOT_AR:
        plot_coef(coefs_AR_abs, coefs_AR_abs_S, coefs_AR_abs_C, "aR_abs", titre, mesh, Nw, Nbeta, indicesAR, True)
        plot_coef(coefs_AR_ph, coefs_AR_ph_S, coefs_AR_ph_C, "aR_ph", titre, mesh, Nw, Nbeta, indicesAR, True)
        plot_coef(coefs_AR_Re, coefs_AR_Re_S, coefs_AR_Re_C, "aR_Re", titre, mesh, Nw, Nbeta, indicesAR, True)
        plot_coef(coefs_AR_Imag, coefs_AR_Imag_S, coefs_AR_Imag_C, "aR_Imag", titre, mesh, Nw, Nbeta, indicesAR, True)
    if PLOT_G:
        plot_coef(matrix_G_Re, matrix_G_Re_S, matrix_G_Re_C, "G_Re", titre, mesh, Nw, Nbeta, indicesAR, True)
        plot_coef(matrix_G_Imag, matrix_G_Imag_S, matrix_G_Imag_C, "G_Imag", titre, mesh, Nw, Nbeta, indicesAR, True)
    if PLOT_G_fex:
        plot_G_fex(matrix_G_Re, matrix_G_Imag, matrix_G_Re_S, matrix_G_Imag_S, matrix_G_Re_C, matrix_G_Imag_C, titre, mesh, Nw, Nbeta, indicesAR)
        
        
if COMP_mesh :
    mesh2="cylinderR3_lid"
    meshes=[mesh, mesh2]
    PIT2= f"{chemin}PIT3_{mesh2}/results/OneBodyProblem{type}"
    matrix_D_Re2=f"{PIT2}_DiffractionMatrix_Re.dat"
    matrix_D_Imag2=f"{PIT2}_DiffractionMatrix_Imag.dat"

    matrix_G_Re2=f"{PIT2}_ForceTransferMatrix_Re.dat"
    matrix_G_Imag2=f"{PIT2}_ForceTransferMatrix_Imag.dat"
    coefs_AR_Re2=f"{PIT2}_RadiationCoefficients_Re.dat"
    coefs_AR_Imag2=f"{PIT2}_RadiationCoefficients_Imag.dat"
    
    plot_coefS_mesh(matrix_D_Re, matrix_D_Re2, "D_Re", titre, meshes, Nw, Nbeta, indicesD, False)
    plot_coefS_mesh(matrix_D_Imag, matrix_D_Imag2, "D_Imag", titre, meshes, Nw, Nbeta, indicesD, False)

    for c, i in enumerate(DOF) : 
        indicesAR=[[6, i-1],[7, i-1], [5, i-1]]
        # indicesAR=[[6, i-1]]
        plot_coefS_mesh(coefs_AR_Re, coefs_AR_Re2, "aR_Re", titre, meshes, Nw, Nbeta, indicesAR, True)
        plot_coefS_mesh(coefs_AR_Imag, coefs_AR_Imag2, "aR_Imag", titre, meshes, Nw, Nbeta, indicesAR, True)
        plot_coefS_mesh(matrix_G_Re, matrix_G_Re2, "G_Re", titre, meshes, Nw, Nbeta, indicesAR, True)
        plot_coefS_mesh(matrix_G_Imag, matrix_G_Imag2, "G_Imag", titre, meshes, Nw, Nbeta, indicesAR, True)
