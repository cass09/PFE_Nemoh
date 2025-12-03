import numpy as np
import re  # Pour utiliser les expressions régulières
import matplotlib.pyplot as plt

def read_RAD_PIT(fileRAD):
    frequencies = []  
    matrices = []     
    current_matrix = [] 
    previous_period = None  # Pour suivre la fréquence précédente

    with open(fileRAD, 'r') as file:
        for line in file:
            line = line.strip()  # Enlever les espaces et les retours à la ligne
            if not line or line.startswith('#'):
                continue
            if 'Period' in line:
                    continue 
            if 'Frequency' in line:
                    continue
            # Convertir la ligne en une liste de valeurs flottantes
            values = list(map(float, line.split()))
            
            period = values[0]
            matrix_row = values[1:] 

            if previous_period is None or period == previous_period:
                # Si c'est la première ligne ou la même fréquence, ajouter la ligne à la matrice courante
                current_matrix.append(matrix_row)
            else:
                # Si la fréquence change, sauvegarder la matrice précédente et recommencer pour la nouvelle fréquence
                matrices.append(current_matrix)
                frequencies.append(previous_period)
                current_matrix = [matrix_row]  # Commencer une nouvelle matrice pour la nouvelle fréquence
            previous_period = period  # Mettre à jour la fréquence précédente

        if current_matrix:
            matrices.append(current_matrix)
            frequencies.append(previous_period)
    return frequencies, matrices

def read_Fex_PIT(fileFex):
    beta_values = []  
    matrices = []     
    frequencies = []  
    current_matrix = []
    previous_beta = None

    with open(fileFex, 'r') as file:
        for line in file:
            line = line.strip()  
            if not line or line.startswith('#'):
                continue
            if 'direction' in line:
                    continue 
            # Convertir la ligne en une liste de valeurs flottantes
            values = list(map(float, line.split()))
            beta = values[0]
            period = values[1]
            matrix_row = values[2:]  

            if previous_beta is None or beta == previous_beta:
                # Si c'est la première ligne ou la même fréquence, ajouter la ligne à la matrice courante
                if period not in frequencies:
                    frequencies.append(period)
                current_matrix.append(matrix_row)
            else:
                # Si la fréquence change, sauvegarder la matrice précédente et recommencer pour la nouvelle fréquence
                matrices.append(current_matrix)
                beta_values.append(previous_beta)
                current_matrix = [matrix_row]  
            previous_beta = beta  # Mettre à jour la fréquence précédente
        # Ajouter la dernière matrice lue (si elle existe)
        if current_matrix:
            matrices.append(current_matrix)
            beta_values.append(previous_beta)
    return beta_values, frequencies, matrices

def Trace_RAD(files, nom, indice, titre, distance, NOMS, mesh):
    
    couleurs = ['b', 'C1', 'g', 'r', 'gold', 'C4', 'c', 'm', 'y', 'k']  
    plt.figure(figsize=(8, 6))
    if nom=="Added_Mass" :
        plt.title(f"{titre} - Added Mass A_{indice[0]}_{indice[1]} ")
    else :
        plt.title(f"{titre} - Damping B_{indice[0]}_{indice[1]} ")
    
    plt.xlabel("Frequency (rad/s)")
    if nom=="Added_Mass" :
        unit="Added Mass (kg)"
    else :
        unit="Damping (kg/s)"
    plt.ylabel(f"{unit} {indice[0]}-{indice[1]}")
    plt.grid(True)

    # Labels des coefficients
    labels = ["M_{11}", "M_{12}", "M_{21}", "M_{22}"]
    styles = ['o-', 's-', 'o', 's']  

    # Tracé des courbes pour BEM
    # for idx, (i, j) in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)]):
    # for idx, (i, j) in enumerate([(0, 0)]):
        # Accède aux éléments dans chaque matrice de Coef_BEM pour chaque fréquence
    i=indice[0]-1
    j=indice[1]-1

    for k, file in enumerate(files) : 
        frequencies, Coef = read_RAD_PIT(file)
    # freq_N2, Coef_N2 = read_RAD_BEM(fileN2)
        plt.plot(frequencies, [mat[i][j] for mat in Coef], 'o-', color=couleurs[k], linewidth=0.8, markersize=4, label=f"{NOMS[k]}")
    # plt.plot(freq_N2, [mat[i][j] for mat in Coef_N2], 's-', label=f"{NOMS[0]}_{indice[0]}{indice[1]}")

    # Légende et affichage
    plt.legend(loc='best', fontsize='small', frameon=True)
    # plt.savefig(f"IT_{mesh}_{nom}_{distance}_M{indice[0]}{indice[1]}.png", dpi=300)
    plt.savefig(f"IT_{mesh}_{nom}_{distance}_M{indice[0]}{indice[1]}.pdf")
    # plt.show()
    return


def Trace_Fex(files, files_ph, num_dof, num_beta, titre, distance, NOMS, mesh):
    # beta_BEM, freq_BEM, Coef_BEM_abs, Coef_BEM_ph = read_Fex_BEM(fileBEM)
    # beta_N2, freq_N2, Coef_N2_abs, Coef_N2_ph = read_Fex_BEM(fileN2)
    i=num_beta-1
    # Tracé des matrices absolues
    couleurs = ['b', 'C1', 'g', 'r', 'gold', 'C4', 'c', 'm', 'y', 'k'] 
    plt.figure(figsize=(10, 8))
    plt.title(f"{titre} - Excitation Force Amplitude |F|_{num_dof}")
    plt.xlabel("Frequency (rad/s)")
    plt.ylabel(f"|F| (N) - DOF {num_dof}")
    
    
    j=num_dof-1

    for k, file in enumerate(files) : 
        beta, freq, Coef_abs = read_Fex_PIT(file)
        # print(k, Coef_abs[i][:, j], Coef_ph[i][:, j])
        plt.plot(freq, [row[j] for row in Coef_abs[i]], 'o-', color=couleurs[k], label=f"{NOMS[k]}")  # Exemple pour la première colonne (à adapter si nécessaire)

    plt.legend()
    plt.grid(True)
    plt.savefig(f"IT_{mesh}_Fex_abs_{distance}_b{num_dof}.pdf")

    # Tracé des matrices absolues
    i=num_beta-1
    plt.figure(figsize=(10, 8))
    plt.title(f"{titre} - Excitation Force Phase ph(F)_{num_dof}")
    plt.xlabel("Frequency (rad/s)")
    plt.ylabel(f"phase(F) (N) - DOF {num_dof}")

    
    j=num_dof-1
    for k, file in enumerate(files_ph) : 
        beta, freq, Coef_ph = read_Fex_PIT(file)
        plt.plot(freq, [row[j] for row in Coef_ph[i]], 'o-', color=couleurs[k], label=f"{NOMS[k]}")  # Exemple pour la première colonne (à adapter si nécessaire)

    plt.legend()
    plt.grid(True)
    plt.savefig(f"IT_{mesh}_Fex_phase_{distance}_b{num_dof}.pdf")
    # plt.show()
    return


def Calcul_GOF(REF, files, Nom_data, dof_vect, NOMS, mesh):
    min_CoF = np.zeros((len(dof_vect), len(files)))
    mean_CoF = np.zeros((len(dof_vect), len(files)))
    use_nrmse=True
    for k in range(len(files)):
        frequencies, COEF = read_RAD_PIT(files[k])
        frequencies, COEF_REF = read_RAD_PIT(REF)
        for j, dof in enumerate(dof_vect):
            diffs_squared = []
            ref_values = []
            gof_local = []
            for i in range(len(frequencies)):
                if frequencies[i] < 2:
                    value = COEF[i][dof - 1][dof - 1]
                    value_ref = COEF_REF[i][dof - 1][dof - 1]

                    diff_sq = np.abs(value - value_ref)**2
                    ref_abs = np.abs(value_ref)

                    diffs_squared.append(diff_sq)
                    ref_values.append(ref_abs)

                    if not use_nrmse:
                        gof_i = 100 * (1 - (np.sqrt(diff_sq) / ref_abs))
                        gof_local.append(gof_i)

            if use_nrmse:
                mse = np.mean(diffs_squared)
                rmse = np.sqrt(mse)
                mean_ref = np.mean(ref_values)
                nrmse = rmse / mean_ref
                gof = (1 - nrmse) * 100
                mean_CoF[j, k] = gof
            else:
                gof = np.mean(gof_local)

                min_gof = np.min(gof_local)
                mean_CoF[j, k] = gof
                min_CoF[j, k] = min_gof
            if gof < 99:
                print(f"dof={dof} - mesh={NOMS[k]} - mean GoF={gof:.2f}% ")
            if gof < 0:
                print(f" !!!!!! dof={dof} - mesh={NOMS[k]} - mean GoF={gof:.2f}% ")
                print(f"{len(diffs_squared)} - diffs_squared={diffs_squared}")
    Trace_GOF(mean_CoF, "mean", Nom_data, NOMS, dof_vect, len(files), mesh, True)
    # Trace_GOF(min_CoF, "min", Nom_data, NOMS, dof_vect, len(files), mesh)


def Calcul_GOF_Fex(REF, files, Nom_data, dof_vect, NOMS, mesh):
    min_CoF=np.zeros((len(dof_vect), len(files)))
    mean_CoF=np.zeros((len(dof_vect), len(files)))
    use_nrmse = True  # ⬅️ Change à False pour utiliser la méthode locale

    for k in range(len(files)) : 
        beta, frequencies, Coef_abs = read_Fex_PIT(files[k])
        beta, frequencies, Coef_abs_REF = read_Fex_PIT(REF)
        for j, dof in enumerate(dof_vect):
            abs_diffs_squared = []
            abs_refs = []
            gof_local_abs = []
            for i in range(len(frequencies)):
                value = Coef_abs[0][i][dof - 1]
                value_ref = Coef_abs_REF[0][i][dof - 1]
                diff_sq_abs = np.abs(value - value_ref)**2
                ref_abs = np.abs(value_ref)

                if frequencies[i] < 2:
                    abs_diffs_squared.append(diff_sq_abs)
                    abs_refs.append(ref_abs)
                    if not use_nrmse:
                        gof_i_abs = 100 * (1 - (np.sqrt(diff_sq_abs) / ref_abs))
                        gof_local_abs.append(gof_i_abs)

            if use_nrmse:
                rmse_abs = np.sqrt(np.mean(abs_diffs_squared))
                mean_abs_ref = np.mean(abs_refs)
                nrmse_abs = rmse_abs / mean_abs_ref
                gof_abs = (1 - nrmse_abs) * 100
            else:
                gof_abs = np.mean(gof_local_abs)

            mean_CoF[j, k] = gof_abs

            if gof_abs < 99:
                print(f"dof={dof} - mesh={NOMS[k]} - ABS mean GoF={gof_abs:.2f}%")
            
    Trace_GOF(mean_CoF, "mean", Nom_data, NOMS, dof_vect, len(files), mesh, False)
    # Trace_GOF(min_CoF, "min", "|Fex|", NOMS, dof_vect, len(files), mesh)
   
def Trace_GOF(GoF, type, data, NOMS, dof_vect, Nfiles, mesh, Coefs_croises) :
    # with open(f"GoodnessOfFit_{mesh}_{data}_{type}.dat", "w") as f:
    #     # En-tête avec les noms de fichiers
    #     f.write(f"# Goodness of Fit - {type} value ({data}) \n")
    #     f.write("DOF\\Mesh \n"+"".join([f"\t{nom}" for nom in NOMS]) + "\n")

    #     for j, dof in enumerate(dof_vect):
    #         f.write(f"{dof}")
    #         for k in range(Nfiles):
    #             f.write(f"\t{GoF[j, k]:.4f}")
    #         f.write("\n")
    
    # Création du graphique
    for j, dof in enumerate(dof_vect):
        plt.figure(figsize=(8, 5))
        bars = plt.bar(NOMS[:-1], GoF[j,:-1], color='skyblue', edgecolor='black')

        # Ajouter les valeurs au-dessus des barres
        for bar in bars:
            yval = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2, yval + 0.5, f'{yval:.2f}%', ha='center', va='bottom')

        # Mise en forme
        plt.ylim(0, 105)
        plt.xlabel('Mesh')
        plt.ylabel(f'{type} GoF (%)')
        if Coefs_croises :
            if data=="Added_Mass" :
                plt.title(f"Goodness of Fit - {mesh} - Added Mass A_{dof}_{dof} (reference value {NOMS[-1]})")
            else :
                plt.title(f"Goodness of Fit - {mesh} - Damping B_{dof}_{dof} (reference value {NOMS[-1]})")
        else :
            if data=="|Fex|" :
                plt.title(f"Goodness of Fit - {mesh} - Excitaion Forces Amplitude |F|_{dof} (reference value {NOMS[-1]})")
            else :
                plt.title(f"Goodness of Fit - {mesh} - Excitaion Forces Phase ph(F)_{dof} (reference value {NOMS[-1]})")
        plt.grid(axis='y', linestyle='--', alpha=0.7)

        plt.tight_layout()
        # plt.show()
        if Coefs_croises :
            plt.savefig(f"GOF_{mesh}_{data}_{type}_{dof}_{dof}.pdf")
        else :
            plt.savefig(f"GOF_{mesh}_{data}_{type}_dof{dof}.pdf")
        plt.close()


def Trace_CPU(Time, NOMS, mesh, unit) :
    # Création du graphe à bâtons
    plt.figure(figsize=(8, 5))
    bars = plt.bar(NOMS, Time, color='coral', edgecolor='black')

    # Ajout des valeurs au-dessus des barres
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval + 1, f'{yval:.2f}', ha='center', va='bottom')

    # Mise en forme
    plt.xlabel('Mesh')
    plt.ylabel(f'Computation Time ({unit})')
    plt.title(f'Mesh Convergence - {mesh} - Computation Time')
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()
    # plt.show()
    plt.savefig(f"CPU_{mesh}_{unit}.pdf")
    plt.close()


chemin= "/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/"

distance="Nb1"
geo = "FB1"
GRAPHS = True
GoF=True
CPU=False
if geo=="cylinderR3" :
    layout="Nb2_X_d1"
    folder="capytaine/MyTestCases/"
    mesh="cylinder"    
    file=[f"{folder}C_BEM_cylinderR3", f"{folder}C_BEM_cylinderR3_SYM"]
    Nm=2
    NOMS=["CAP Cylinder", "CAP Cylinder SYM"]
    CPU_data=[3550.07861, 4950.71436, 5967.88965, ]
elif geo=="cylinderR3SYM_lid" :
    layout="Nb1_d0"
    folder="capytaine/MyTestCases/Mesh_Convergence/C_BEM_cylinderR3_SYM_lid_m"
    mesh="cylinder"    
    file=[f"{folder}1", f"{folder}2", f"{folder}3", f"{folder}4", f"{folder}5", f"{folder}6", f"{folder}7", f"{folder}8", "capytaine/MyTestCases/Mesh_Convergence/C_BEM_cylinderR3_lid_m5"]
    # file=[f"{folder}1", f"{folder}3", f"{folder}5", f"{folder}8", "capytaine/MyTestCases/Mesh_Convergence/C_BEM_cylinderR3_lid_m5"]
    Nm=8
    NOMS=["Nm=456", "Nm=638", "Nm=918", "Nm=1080", "Nm=1386", "Nm=1540", "Nm=1824", "Nm=2376"]
    # NOMS=["Nm=456", "Nm=918", "Nm=1386", "Nm=2376", "NO SYM"]
    CPU_data=[910.220459, 1727.05688, ]
elif geo=="cylinderR3_Nbeta" :
    layout="Nb2_X_d1"
    folder="capytaine/MyTestCases/CAP_cylinderR3_Nbeta"
    mesh="cylinder"    
    # file=[f"{folder}3", f"{folder}7", f"{folder}11", f"{folder}13", f"{folder}15", f"{folder}17", f"{folder}19"]
    file=[f"{folder}11", f"{folder}13", f"{folder}15", f"{folder}17", f"{folder}19"]
    Nm=3
    # NOMS=["Nbeta=3", "Nbeta=7", "Nbeta=11", "Nbeta=13", "Nbeta=15", "Nbeta=17", "Nbeta=19"]
    NOMS=["Nbeta=11", "Nbeta=13", "Nbeta=15"]
    CPU_data=[910.220459, 1727.05688, ]
elif geo=="CylR3" :
    layout="Nb1_d0"
    folder="capytaine/MyTestCases/Mesh_CylR3/C_BEM_CylR3_m"
    mesh="cylinder"  
    file=[f"{folder}1", f"{folder}2", f"{folder}3", f"{folder}4", f"{folder}5", f"{folder}6", f"{folder}7", f"{folder}8"]
    Nm=8
    NOMS=["Np=361", "Np=506", "Np=729", "Np=870", "Np=1089", "Np=1225", "Np=1444", "Np=1892"]
    CPU_data=[375.908936, 674.773438, 1341.02722, 1862.98352, 2914.78955, 3654.30127, 5382.15234, 9549.78027]
elif geo=="BargeX6" :
    layout="Nb1_d0"
    folder="capytaine/MyTestCases/Mesh_BargeX6/C_BEM_BargeX6_m"
    mesh="barge"  
    file=[f"{folder}1", f"{folder}2", f"{folder}3", f"{folder}4", f"{folder}5", f"{folder}6", f"{folder}7", f"{folder}8"]
    Nm=8
    NOMS=["Np=403", "Np=540", "Np=810", "Np=1000", "Np=1210", "Np=1525", "Np=1782", "Np=2059"]
    CPU_data=[460.933411, 748.617188, 1637.47546, 2293.85913, 4028.12866, 5731.54297, 8317.33984, 9756.90918]
elif geo=="FB1" :
    layout="Nb1_d0"
    folder="capytaine/MyTestCases/C_BEM_FB1_m"
    mesh="FB1"  
    file=[f"{folder}1", f"{folder}2", f"{folder}3", f"{folder}4", f"{folder}5", f"{folder}6", f"{folder}7", f"{folder}8"]
    Nm=7
    NOMS=["Np=148", "Np=592", "Np=925", "Np=1332", "Np=1813", "Np=2368", "Np=3700"]
    # NOMS=["Np=148", "Np=592", "Np=925", "Np=1332", "Np=1813", "Np=2368", "Np=3700", "Np=5328"]
    CPU_data=[460.933411, 748.617188, 1637.47546, 2293.85913, 4028.12866, 5731.54297, 8317.33984, 9756.90918]

# titre=f"Nw=20, Nbeta=11, dof=6, Ndir=1, {distance}"
titre=f"Mesh Convergence - {mesh}"




CMfiles=[]
CAfiles=[]
Fefiles=[]
Fepfiles=[]
for i in range(Nm) :
    CMfiles.append(f"{chemin}{file[i]}/resultsBEM/{layout}/Capytaine_Madd.dat")
    CAfiles.append(f"{chemin}{file[i]}/resultsBEM/{layout}/Capytaine_Crad.dat")
    Fefiles.append(f"{chemin}{file[i]}/resultsBEM/{layout}/Capytaine_Fe_abs.dat")
    Fepfiles.append(f"{chemin}{file[i]}/resultsBEM/{layout}/Capytaine_Fe_phase.dat")
    # CMfiles.append(f"{chemin}{file[i]}/resultsIT/CapytaineIT_S_Madd_{layout}.00.dat")
    # CAfiles.append(f"{chemin}{file[i]}/resultsIT/CapytaineIT_S_Crad_{layout}.00.dat")
    # Fefiles.append(f"{chemin}{file[i]}/resultsIT/CapytaineIT_S_Fe_abs_{layout}.00.dat")
    # Fepfiles.append(f"{chemin}{file[i]}/resultsIT/CapytaineIT_S_Fe_phase_{layout}.00.dat")


if GRAPHS :
    for ind, dof in enumerate([1, 3, 5]) : 
        Trace_RAD(CAfiles, "Damping", [dof,dof], titre, distance, NOMS, mesh)
        Trace_RAD(CMfiles, "Added_Mass", [dof,dof], titre, distance, NOMS, mesh)
        Trace_Fex(Fefiles, Fepfiles, dof, 1, titre, distance, NOMS, mesh)

if GoF:
    dof=[1, 3, 5]
    print("------------GoF-----------------")
    print("Added Mass")
    Calcul_GOF(CMfiles[-1], CMfiles, "Added_Mass", dof, NOMS, mesh)
    print("Damping")
    Calcul_GOF(CAfiles[-1], CAfiles, "Damping", dof, NOMS, mesh)
    print("Fex")
    Calcul_GOF_Fex(Fefiles[-1], Fefiles, "|Fex|", dof, NOMS, mesh)
    Calcul_GOF_Fex(Fepfiles[-1], Fepfiles, "phase(Fex)", dof, NOMS, mesh)

if CPU:
    # Trace_CPU(CPU_data, NOMS, mesh, "s")
    CPU_data=np.array(CPU_data)/60
    Trace_CPU(CPU_data, NOMS, mesh, "min")

