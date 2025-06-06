import numpy as np
import re  # Pour utiliser les expressions régulières
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import itertools


def CompareResultsE(file, Nw, Ne) :
    dof=6
    Nbeta=11
    data= np.loadtxt(file, skiprows=0, max_rows=Nw*Nbeta*(Ne+1))
    freq = data[:, 0]        
    Matrix = data[:,1:dof+1]
    relative_errorL2=np.zeros((Nw, Ne))
    relative_errorINF=np.zeros((Nw, Ne))
    w=np.zeros(Nw)
    for i in range(Nw) :
        w[i]=freq[i*Nbeta*(Ne+1)+2]
        Matw=Matrix[i*(Ne+1)*Nbeta:(i+1)*Nbeta*(Ne+1), :]
        Pwave=Matw[0:Nbeta, 0:Nbeta]
        for j in range(Ne) :
            Ewave=Matw[Nbeta*(j+1):Nbeta*(j+2), :]
            errorL2 = np.sqrt(np.sum((np.array(Pwave) - np.array(Ewave))**2))
            erreurs = np.abs(np.array(Pwave) - np.array(Ewave))
            errorINF = np.max(erreurs)
            
            normL2_BEM = np.sqrt(np.sum(np.array(Pwave)**2))
            relative_errorL2[i, j] = errorL2 / normL2_BEM if normL2_BEM != 0 else np.nan

            max_BEM = np.max(np.abs(np.array(Pwave)))
            relative_errorINF[i, j] = errorINF / max_BEM if max_BEM != 0 else np.nan
  
    return w, relative_errorL2, relative_errorINF

def Ratio_E(file, type, Nw, Ne) : 
    dof=6
    Nbeta=11
    data= np.loadtxt(file, skiprows=0, max_rows=Nw*Nbeta*(Ne+1))
    freq = data[:, 0]     
    if type :    
        Matrix = data[:,1:Nbeta*(Ne+1)]
    else : 
        Matrix = data[:,1:dof+1]
    Ratio=np.zeros(Nw)
    Ratio_tot=np.zeros(Nw)
    w=np.zeros(Nw)
    for i in range(Nw) :
        w[i]=freq[i*Nbeta*(Ne+1)+2]
        Matw=Matrix[i*(Ne+1)*Nbeta:(i+1)*Nbeta*(Ne+1), :]
        MatP=Matw[:Nbeta, :]
        MatE=Matw[Nbeta:,:]
        # Norme de Frobenius
        norm_ev = np.linalg.norm(MatE)
        norm_p = np.linalg.norm(MatP)
        norm_tot = np.linalg.norm(Matrix)
        Ratio[i]=norm_ev/norm_p
        Ratio_tot[i]=norm_ev/norm_tot
    return w, Ratio, Ratio_tot

def plot_Ratio(file, nom, titre, Ne) : 
    data = np.loadtxt(file, delimiter=None, skiprows=1) 

    w = data[:, 0]
    # 1 abs Ratio / progressive, 2 abs Ratio / tot
    # 3 phase Ratio / progressive, 4 phase Ratio / tot
    p_abs_tot=data[:,2]
    plt.figure(figsize=(8, 6))
    plt.plot(w, p_abs_tot, marker='o', linestyle='-')
    plt.xlabel("frequencies (rad/s)")
    plt.ylabel(f"Ratio")
    plt.title(f"Norm ratio of {nom} for Ne={Ne} \n ({titre})")
    # plt.legend()
    plt.grid()
    plt.savefig(f"Ratio_E{Ne}_{nom}.png")
    plt.close()  # Ferme la figure pour éviter trop d'ouvertures

def plot_error(file, nom, titre, Ne) : 
    data = np.loadtxt(file, delimiter=None, skiprows=1) 
    colors = itertools.cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    w = data[:, 0]
    # (error INF - error L2 ) x Ne =2
    plt.figure(figsize=(8, 6))
    for i in range(Ne) :
        error=data[:,i+1]
        plt.plot(w, error, marker='o', color=next(colors), linestyle='-', label=f'mode E {i+1}')
    plt.xlabel("frequencies (rad/s)")
    plt.ylabel(f"Relative Errors")
    plt.title(f"Relative Errors of {nom} for Ne={Ne} \n ({titre})")
    plt.legend()
    plt.grid()
    plt.savefig(f"EcartsRelatifs_E{Ne}_{nom}.png")
    plt.close()  # Ferme la figure pour éviter trop d'ouvertures


chemin="/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/PIT3_E/wec_inputs/"
Calcul_Ratio = True
Error=False
PLOT=True
Nw=20
Ne=6
type=f"E{Ne}"
titre=f"Cylinder - Nw={Nw}, Nbeta=11, dof=6, Ndir=1"
PIT= f"{chemin}PIT3_Cylinder_dof6/results/OneBodyProblem_{type}"

matrix_D_abs=f"{PIT}_DiffractionMatrix_abs.dat"
matrix_D_ph=f"{PIT}_DiffractionMatrix_ph.dat"
matrix_G_abs=f"{PIT}_ForceTransferMatrix_abs.dat"
matrix_G_ph=f"{PIT}_ForceTransferMatrix_ph.dat"
coefs_AR_abs=f"{PIT}_RadiationCoefficients_abs.dat"
coefs_AR_ph=f"{PIT}_RadiationCoefficients_ph.dat"
coefs_AS_abs=f"{PIT}_ScatteringCoefficients_abs.dat"
coefs_AS_ph=f"{PIT}_ScatteringCoefficients_ph.dat"

if Calcul_Ratio : 
    print("Calcul Ratio")
    fileD=f"Ratio_{type}_D.dat"
    with open(fileD, "w") as f:        
        w, Ratio, Ratio_tot = Ratio_E(matrix_D_abs, True, Nw, Ne)
        w, Ratio_ph, Ratio_tot_ph = Ratio_E(matrix_D_ph, True, Nw, Ne)
        f.write(f"matrix D for Ne={Ne} - freq  - abs (Ratio - Ratio total) - phase  (Ratio - Ratio total) \n")
        for i in range(Nw) :
            f.write(f"{w[i]:.4f}    {Ratio[i]:16.8f}        {Ratio_tot[i]:16.8f}           {Ratio_ph[i]:16.8f}        {Ratio_tot_ph[i]:16.8f} ")  
            f.write("\n") 
    fileG=f"Ratio_{type}_G.dat"
    with open(fileG, "w") as f:        
        w, Ratio, Ratio_tot = Ratio_E(matrix_G_abs, False, Nw, Ne)
        w, Ratio_ph, Ratio_tot_ph = Ratio_E(matrix_G_ph, False, Nw, Ne)
        f.write(f"matrix G for Ne={Ne} - freq  - abs (Ratio - Ratio total) - phase  (Ratio - Ratio total) \n")
        for i in range(Nw) :
            f.write(f"{w[i]:.4f}    {Ratio[i]:16.8f}        {Ratio_tot[i]:16.8f}           {Ratio_ph[i]:16.8f}        {Ratio_tot_ph[i]:16.8f} ")  
            f.write("\n") 
    fileAR=f"Ratio_{type}_AR.dat"
    with open(fileAR, "w") as f:        
        w, Ratio, Ratio_tot = Ratio_E(coefs_AR_abs, False, Nw, Ne)
        w, Ratio_ph, Ratio_tot_ph = Ratio_E(coefs_AR_ph, False, Nw, Ne)
        f.write(f"coefs AR for Ne={Ne} - freq  - abs (Ratio - Ratio total) - phase  (Ratio - Ratio total) \n")
        for i in range(Nw) :
            f.write(f"{w[i]:.4f}    {Ratio[i]:16.8f}        {Ratio_tot[i]:16.8f}           {Ratio_ph[i]:16.8f}        {Ratio_tot_ph[i]:16.8f} ")  
            f.write("\n") 
    fileAS=f"Ratio_{type}_AS.dat"
    with open(fileAS, "w") as f:        
        w, Ratio, Ratio_tot = Ratio_E(coefs_AS_abs, False, Nw, Ne)
        w, Ratio_ph, Ratio_tot_ph = Ratio_E(coefs_AS_ph, False, Nw, Ne)
        f.write(f"coefs AS for Ne={Ne} - freq  - abs (Ratio - Ratio total) - phase  (Ratio - Ratio total) \n")
        for i in range(Nw) :
            f.write(f"{w[i]:.4f}    {Ratio[i]:16.8f}        {Ratio_tot[i]:16.8f}           {Ratio_ph[i]:16.8f}        {Ratio_tot_ph[i]:16.8f} ")  
            f.write("\n") 
    if PLOT :
        plot_Ratio(fileD, "matrix_D", titre, Ne)
        plot_Ratio(fileG, "matrix_G", titre, Ne)
        plot_Ratio(fileAR, "coefs_AR", titre, Ne)
        plot_Ratio(fileAS, "coefs_AS", titre, Ne)

if Error :                
    print("Calcul Relative Error")
    fileGabs=f"EcartsRelatifs_{type}_G_abs.dat"
    with open(fileGabs, "w") as f:        
        w, L_inf_error, L2_error = CompareResultsE(matrix_G_abs, Nw, Ne)
        f.write(f"matrix G - freq  - (error INF - error L2 ) x Ne={Ne} \n")
        for i in range(Nw) :
            f.write(f"{w[i]:.4f}")
            for j in range(Ne) :
                f.write(f"    {L_inf_error[i, j]:16.8f}     {L2_error[i, j]:16.8f}      ")        
            f.write("\n")
    fileARabs=f"EcartsRelatifs_{type}_AR_abs.dat"
    with open(fileARabs, "w") as f:        
        w, L_inf_error, L2_error = CompareResultsE(coefs_AR_abs, Nw, Ne)
        f.write(f"coef AR - freq  - (error INF - error L2 ) x Ne ={Ne}\n")
        for i in range(Nw) :
            f.write(f"{w[i]:.4f}")
            for j in range(Ne) :
                f.write(f"    {L_inf_error[i, j]:16.8f}     {L2_error[i, j]:16.8f}      ")        
            f.write("\n")
    fileASabs=f"EcartsRelatifs_{type}_AS_abs.dat"
    with open(fileASabs, "w") as f:        
        w, L_inf_error, L2_error = CompareResultsE(coefs_AS_abs, Nw, Ne)
        f.write(f"coef AS - freq  - (error INF - error L2 ) x Ne ={Ne}\n")
        for i in range(Nw) :
            f.write(f"{w[i]:.4f}")
            for j in range(Ne) :
                f.write(f"    {L_inf_error[i, j]:16.8f}     {L2_error[i, j]:16.8f}      ")        
            f.write("\n")
    fileGph=f"EcartsRelatifs_{type}_G_ph.dat"
    with open(fileGph, "w") as f:        
        w, L_inf_error, L2_error = CompareResultsE(matrix_G_ph, Nw, Ne)
        f.write(f"matrix G - freq  - (error INF - error L2 ) x Ne={Ne} \n")
        for i in range(Nw) :
            f.write(f"{w[i]:.4f}")
            for j in range(Ne) :
                f.write(f"    {L_inf_error[i, j]:16.8f}     {L2_error[i, j]:16.8f}      ")        
            f.write("\n")    
    fileARph=f"EcartsRelatifs_{type}_AR_ph.dat"
    with open(fileARph, "w") as f:        
        w, L_inf_error, L2_error = CompareResultsE(coefs_AR_ph, Nw, Ne)
        f.write(f"coef AR - freq  - (error INF - error L2 ) x Ne={Ne} \n")
        for i in range(Nw) :
            f.write(f"{w[i]:.4f}")
            for j in range(Ne) :
                f.write(f"    {L_inf_error[i, j]:16.8f}     {L2_error[i, j]:16.8f}      ")        
            f.write("\n")
    fileASph=f"EcartsRelatifs_{type}_AS_ph.dat"
    with open(fileASph, "w") as f:        
        w, L_inf_error, L2_error = CompareResultsE(coefs_AS_ph, Nw, Ne)
        f.write(f"coef AS - freq  - (error INF - error L2 ) x Ne={Ne} \n")
        for i in range(Nw) :
            f.write(f"{w[i]:.4f}")
            for j in range(Ne) :
                f.write(f"    {L_inf_error[i, j]:16.8f}     {L2_error[i, j]:16.8f}      ")        
            f.write("\n")
    if PLOT :
        plot_error(fileGabs, "G_abs", titre, Ne)
        plot_error(fileGph, "G_ph", titre, Ne)
        plot_error(fileARabs, "AR_abs", titre, Ne)
        plot_error(fileARph, "AR_ph", titre, Ne)
        plot_error(fileASabs, "AS_abs", titre, Ne)
        plot_error(fileASph, "AS_ph", titre, Ne)
