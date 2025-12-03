import numpy as np
import re  # Pour utiliser les expressions régulières
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def read_RAD_Nemoh(fileRAD):
    frequencies = []  
    matrices = []     
    current_matrix = []  

    with open(fileRAD, 'r') as file:
        for line in file:
            line = line.strip()  # Enlever les espaces et les retours à la ligne
            if not line or line.startswith('#'):
                continue 

            try:
                if 'Nb' in line:
                    continue  
                frequency = float(line)
                if current_matrix:  # Si une matrice est en cours, l'ajouter aux matrices
                    matrices.append(current_matrix)
                    current_matrix = []  # Réinitialiser la matrice pour la prochaine fréquence
                frequencies.append(frequency)  
            except ValueError:
                # Si la ligne ne peut pas être convertie en nombre, on suppose que c'est une ligne de la matrice
                matrix_row = list(map(float, line.split()))  # Convertir la ligne en une liste de flottants
                current_matrix.append(matrix_row)  # Ajouter la ligne à la matrice en cours

        # Ajouter la dernière matrice (si elle existe) à la fin du fichier
        if current_matrix:
            matrices.append(current_matrix)
        
        # Case ONE body !!!
        ONE_BODY = False
        if ONE_BODY:
            matrices = [[[value]] for value in frequencies[1::2]]  # Indices impairs
            frequencies = frequencies[::2]   # Indices pairs
    return frequencies, matrices

def read_RAD_Capytaine(fileRAD):
    frequencies = []  
    matrices = []     
    current_matrix = [] 
    previous_period = None  # Pour suivre la fréquence précédente

    with open(fileRAD, 'r') as file:
        for line in file:
            line = line.strip()  # Enlever les espaces et les retours à la ligne
            if not line or line.startswith('#'):
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

def CompareResultsRAD(fileNemoh, fileCapytaine, nom, comment, param_d, type, Nb, tol, precision) :
    freq_Nemoh, Coef_Nemoh = read_RAD_Nemoh(fileNemoh)
    freq_Capytaine, Coef_Capytaine = read_RAD_Capytaine(fileCapytaine)
    if np.all(np.isclose(freq_Nemoh, freq_Capytaine, atol=1e-4)):
        if np.shape(Coef_Nemoh) != np.shape(Coef_Capytaine):
            print("Nemoh", np.shape(Coef_Nemoh), "Capytaine", np.shape(Coef_Capytaine))
            raise ValueError("tailles différentes")
        else :
            errorL2 = np.sqrt(np.sum((np.array(Coef_Nemoh) - np.array(Coef_Capytaine))**2))
            errorINF = np.max(np.abs(np.array(Coef_Nemoh) - np.array(Coef_Capytaine)))

            normL2_Nemoh = np.sqrt(np.sum(np.array(Coef_Nemoh)**2))
            relative_errorL2 = errorL2 / normL2_Nemoh if normL2_Nemoh != 0 else np.nan

            max_Nemoh = np.max(np.abs(np.array(Coef_Nemoh)))
            relative_errorINF = errorINF / max_Nemoh if max_Nemoh != 0 else np.nan

            All_errors=True
            if All_errors:
                dof=6*Nb
                n_error=0
                for i in range(dof):
                    for j in range(dof):
                        Nemoh = []
                        Capytaine = []                       

                        Nemoh = [0 if np.abs(sublist[i][j]) < precision else sublist[i][j] for sublist in Coef_Nemoh]
                        Capytaine = [0 if np.abs(sublist[i][j]) < precision else sublist[i][j] for sublist in Coef_Capytaine]
                        # Remplacer les valeurs < precision par 0
                        if np.max(np.abs(np.array(Nemoh))) != 0 and np.max(np.abs(np.array(Capytaine)))!=0:
                            error=np.max(np.abs(np.array(Nemoh)-np.array(Capytaine)))/np.max(np.abs(np.array(Nemoh)))
                        else :
                            error=0
                        with open(f"EcartsRelatifs_{comment}{nom}_{type}.dat", "a") as f:
                            if i==0 and j==0:
                                f.write(f"\n ")
                            if j==0:
                                f.write(f" \n {i}     ")
                            f.write(f"{error:16.8f}     ")
                        if error<=tol:
                            n_error=n_error+1
                print("d=", param_d, "--------> n_error", n_error, "/", dof*dof, "<", tol*100, "%")             
    else :
        print("fréquences différentes")
        print('Nemoh', freq_Nemoh)
        print('Capytaine', freq_Capytaine)
        relative_errorL2=[]
        relative_errorINF=[]
    
    return relative_errorL2, relative_errorINF

def read_Fex_Nemoh(fileFex):
    beta_values = []  
    matrices = []     
    matrices_abs = []
    matrices_phase = []
    frequencies = [] 

    with open(fileFex, 'r') as f:
        beta = None  # Initialiser une variable pour beta
        matrix = []   
        for line in f:
            line = line.strip() 
            if 'Zone' in line:
                # Chercher la valeur de beta
                match = re.search(r'beta\s*=\s*([-\d\.eE]+)', line)
                if match:
                    # Si un beta est trouvé, le stocker
                    if beta is not None and matrix:
                        # Si une matrice est déjà en cours pour le précédent beta, on l'ajoute
                        beta_values.append(beta) 
                        matrices.append(np.array(matrix))
                        matrix_abs = []
                        matrix_phase = []
                        for row in matrix:
                            abs_values = row[::2]  # Valeurs absolues (indices 0, 2, 4, ...)
                            phase_values = row[1::2]  # Phases (indices 1, 3, 5, ...)
                            matrix_abs.append(abs_values)
                            matrix_phase.append(phase_values)
                        matrices_abs.append(np.array(matrix_abs))
                        matrices_phase.append(np.array(matrix_phase))
                        matrix = [] 
                    beta = float(match.group(1))  # Récupérer la nouvelle valeur de beta
            # Si la ligne n'est pas vide et qu'on est dans une matrice
            elif line and beta is not None:
                # Convertir la ligne en une liste de valeurs flottantes
                values = list(map(float, line.split()))
                if values:
                    matrix.append(values[1:])  
                    if values[0] not in frequencies:
                        frequencies.append(values[0])  
        # Ajouter la dernière matrice et beta restant à la fin du fichier
        if beta is not None and matrix:
            beta_values.append(beta)
            matrices.append(np.array(matrix)) 
        # Séparer les colonnes de la matrice en valeurs absolues et phases
            matrix_abs = []
            matrix_phase = []
            for row in matrix:
                abs_values = row[::2]  # Valeurs absolues (indices 0, 2, 4, ...)
                phase_values = row[1::2]  # Phases (indices 1, 3, 5, ...)
                matrix_abs.append(abs_values)
                matrix_phase.append(phase_values)
        # Ajouter les résultats dans les listes correspondantes
            matrices_abs.append(np.array(matrix_abs))
            matrices_phase.append(np.array(matrix_phase))   
    return beta_values, frequencies, matrices_abs, matrices_phase

def read_Fex_Capytaine(fileFex):
    beta_values_set = set()      # Pour collecter tous les β distincts
    matrices = []                # Liste de matrices Fex (une par fréquence)
    omega_values = []           # Liste des fréquences ω
    current_matrix = []
    previous_omega = None

    with open(fileFex, 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith('#') or 'direction' in line:
                continue

            values = list(map(float, line.split()))
            beta = values[0]
            omega = values[1]
            matrix_row = values[2:]

            beta_values_set.add(beta)

            if previous_omega is None or omega == previous_omega:
                current_matrix.append(matrix_row)
            else:
                # Changement d’omega => stocker la matrice précédente
                matrices.append(current_matrix)
                omega_values.append(previous_omega)
                current_matrix = [matrix_row]

            previous_omega = omega

        # Ajouter la dernière matrice
        if current_matrix:
            matrices.append(current_matrix)
            omega_values.append(previous_omega)

    beta_values = sorted(beta_values_set)
    return beta_values, omega_values, matrices

def CompareResultsFex(fileNemoh, fileCapytaine, fileCapytaine_phase, comment, param_d, type, Nb, tol, precision) :
    beta_Nemoh, freq_Nemoh, Coef_Nemoh_abs, Coef_Nemoh_ph = read_Fex_Nemoh(fileNemoh)
    beta_Capytaine, freq_Capytaine, Coef_Capytaine = read_Fex_Capytaine(fileCapytaine)
    beta_Capytaine, freq_Capytaine, Phase_Capytaine = read_Fex_Capytaine(fileCapytaine_phase)
    # print("Nemoh", Coef_Nemoh_ph)
    # print("Capytaine", Phase_Capytaine)
    Coef_Capytaine=np.transpose(Coef_Capytaine, (1, 0, 2))
    Phase_Capytaine=np.transpose(Phase_Capytaine, (1, 0, 2))
    if len(freq_Nemoh)==len(freq_Capytaine) and len(beta_Nemoh)==len(beta_Capytaine):
        if np.shape(Coef_Nemoh_abs) != np.shape(Coef_Capytaine):
            print(np.shape(Coef_Nemoh_abs), np.shape(Coef_Capytaine))
            raise ValueError("tailles différentes")
        if np.shape(Coef_Nemoh_ph) != np.shape(Phase_Capytaine):
            print(np.shape(Coef_Nemoh_ph), np.shape(Phase_Capytaine))
            raise ValueError("tailles différentes")
        else :
            errorL2 = np.sqrt(np.sum((np.array(Coef_Nemoh_abs) - np.array(Coef_Capytaine))**2))
            # errorINF = np.max(np.abs(np.array(Coef_Nemoh_abs) - np.array(Coef_Capytaine)))
            erreurs = np.abs(np.array(Coef_Nemoh_abs) - np.array(Coef_Capytaine))
            errorINF = np.max(erreurs)
            indice_max = np.argmax(erreurs)
            # print("indice_max =", indice_max)
            # print("Fex_abs_Nemoh_max =", np.array(Coef_Nemoh_abs)[indice_max])
            # print("Fex_abs_Capytaine_max =", np.array(Coef_Capytaine)[indice_max])

            normL2_Nemoh = np.sqrt(np.sum(np.array(Coef_Nemoh_abs)**2))
            relative_errorL2_abs = errorL2 / normL2_Nemoh if normL2_Nemoh != 0 else np.nan

            max_Nemoh = np.max(np.abs(np.array(Coef_Nemoh_abs)))
            relative_errorINF_abs = errorINF / max_Nemoh if max_Nemoh != 0 else np.nan

            # PHASE 
            errorL2_p = np.sqrt(np.sum((np.array(Coef_Nemoh_ph) - np.array(Phase_Capytaine))**2))
            # errorINF_p = np.max(np.abs(np.array(Coef_Nemoh_ph) - np.array(Phase_Capytaine)))

            erreurs_p = np.abs(np.array(Coef_Nemoh_ph) - np.array(Phase_Capytaine))
            errorINF_p = np.max(erreurs_p)
            # indice_max_p = np.argmax(erreurs_p)
            # print("indice_max =", indice_max_p, Coef_Nemoh_ph[0][6][1], Phase_Capytaine[0][6][1])

            normL2_Nemoh_p = np.sqrt(np.sum(np.array(Coef_Nemoh_ph)**2))
            relative_errorL2_phase = errorL2_p / normL2_Nemoh_p if normL2_Nemoh_p != 0 else np.nan

            max_Nemoh_p = np.max(np.abs(np.array(Coef_Nemoh_ph)))
            relative_errorINF_phase = errorINF_p / max_Nemoh_p if max_Nemoh_p != 0 else np.nan

            # Calcul erreur pour chaque dof 
            All_errors=True
            if All_errors:
                dof=6*Nb
                n_err_abs=0
                n_err_ph=0
                for i in range(dof):
                    B_abs = [N[i] for N in Coef_Nemoh_abs[0]]
                    P_abs = [N[i] for N in Coef_Capytaine[0]]
                    B_ph = [N[i] for N in Coef_Nemoh_ph[0]]
                    P_ph = [N[i] for N in Phase_Capytaine[0]]
                    # Remplacer les valeurs < precision par 0
                    B_abs = [0 if np.abs(x) < precision else x for x in B_abs]
                    P_abs = [0 if np.abs(x) < precision else x for x in P_abs]
                    B_ph = [0 if np.abs(b) < precision else ph for b, ph in zip(B_abs, B_ph)]
                    P_ph = [0 if np.abs(b) < precision else ph for b, ph in zip(P_abs, P_ph)]
                    if np.max(np.abs(np.array(B_abs))) != 0 and np.max(np.abs(np.array(P_abs))) != 0 : 
                        err_abs=np.max(np.abs(np.array(B_abs)-np.array(P_abs)))/np.max(np.abs(np.array(B_abs))) 
                    else : 
                        err_abs=0
                    if np.max(np.abs(np.array(B_ph))) != 0 and np.max(np.abs(np.array(P_ph))) != 0 : 
                        err_ph=np.max(np.abs(np.array(B_ph)-np.array(P_ph)))/np.max(np.abs(np.array(B_ph))) 
                    else : 
                        err_ph=0
                    
                    with open(f"EcartsRelatifs_{comment}Fex_abs_{type}.dat", "a") as f:
                        if i==0:
                            f.write("\n")
                            f.write(f"{param_d}     ")
                        f.write(f"{err_abs:16.8f}     ")
                    if err_abs<=tol:
                        n_err_abs=n_err_abs+1
                        
                    with open(f"EcartsRelatifs_{comment}Fex_ph_{type}.dat", "a") as f:
                        if i==0:
                            f.write("\n")
                            f.write(f"{param_d}     ")
                        f.write(f"{err_ph:16.8f}      ")
                    if err_ph<=tol:
                        n_err_ph=n_err_ph+1
                print("d=", param_d,"--------> n_error (abs, ph)", n_err_abs, n_err_ph, "/", dof, "<", tol*100, "%")
    else :
        print("fréquences différentes")
        print('Nemoh', freq_Nemoh)
        print('Capytaine', freq_Capytaine)
        relative_errorL2_abs=[]
        relative_errorL2_phase=[]
        relative_errorINF_abs=[]
        relative_errorINF_phase=[]
    
    return relative_errorL2_abs, relative_errorINF_abs, relative_errorL2_phase, relative_errorINF_phase


def Trace_RAD(filePIT, fileCapytaine, filePIT_REF, file, nom, indice, titre, comp):
    freq_PIT, Coef_PIT = read_RAD_Nemoh(filePIT)
    freq_Capytaine, Coef_Capytaine = read_RAD_Capytaine(fileCapytaine)
    REF=False
    if filePIT_REF != filePIT :
        freq_PIT_REF, Coef_PIT_REF = read_RAD_Nemoh(filePIT_REF)
        REF=True
    couleurs = ['b', 'g', 'r', 'gold', 'lime', 'c', 'm', 'y', 'k']  # palette de couleurs (réutilisée si plus de 7 courbes)
    markers = ['s', 'o', '<', '^','x', '+', '*']
    
    spectre=freq_PIT
    axeX="Frequency (rad/s)"
    if nom=="Damping":
        unit="(kg/s)"
    else :
        unit="(kg)"
    rho=1000
    g=9.81
    a=3
    # AD=(rho*g*np.pi*a**2*(2*a))/(g/a)**(0.5)
    AD=1
    plt.figure(figsize=(8, 6))
    plt.title(titre)
    plt.xlabel(f"{axeX}")
    plt.ylabel(f"{nom} {unit}")
    # plt.ylabel(f"{nom} normalized (-)")
    plt.grid(True)

    for k, ij in enumerate(indice): 
        i=ij[0]-1
        j=ij[1]-1
        color = couleurs[k % len(couleurs)]
        marker = markers[k % len(markers)]
        color2 = couleurs[k+len(indice) % len(couleurs)]
        marker2 = markers[k+len(indice) % len(markers)]
        if REF :
            color3 = couleurs[k+len(indice)+1 % len(couleurs)]
            marker3 = markers[k+len(indice)+2 % len(markers)]
            plt.plot(freq_PIT_REF, [mat[i][j]/AD for mat in Coef_PIT_REF], linestyle='-', marker=marker3, color=color3, markersize=5, label=f"Outer Cylinder Method {ij[0]}_{ij[1]}")
   
        plt.plot(spectre, [mat[i][j]/AD for c, mat in enumerate(Coef_PIT)], linestyle=':', marker=marker, color=color, markersize=5, label=f"{comp} with Nemoh {ij[0]}_{ij[1]}")
        plt.plot(spectre, [mat[i][j]/AD for c, mat in enumerate(Coef_Capytaine)], linestyle=':', marker=marker2, color=color2, markersize=4, label=f"{comp} with Capytaine {ij[0]}_{ij[1]}")
        plt.legend(loc='best', fontsize='small', frameon=True)
    plt.savefig(f"{file}_{nom}_M{indice}.png", dpi=300)
    plt.close()
    return

def Trace_Fex(filePIT, fileCapytaine, fileCapytaine_phase, filePIT_REF, file, num_dof, beta, titre, comp):
    beta_PIT, freq_PIT, Coef_PIT_abs, Coef_PIT_ph = read_Fex_Nemoh(filePIT)
    beta_Capytaine, freq_Capytaine, Coef_Capytaine = read_Fex_Capytaine(fileCapytaine)
    beta_Capytaine, freq_Capytaine, Phase_Capytaine = read_Fex_Capytaine(fileCapytaine_phase)
    REF=False
    if filePIT_REF != filePIT :
        beta_PIT_REF, freq_PIT_REF, Coef_PIT_REF, Phase_PIT_REF = read_Fex_Nemoh(filePIT_REF)
        REF=True
    couleurs = ['b', 'm', 'r', 'gold', 'lime', 'g', 'c', 'm', 'y', 'k']  # palette de couleurs (réutilisée si plus de 7 courbes)
    markers = ['s', 'o', '<', '^','x', '+', '*']
    if beta : # Several wave directions
        plt.figure(figsize=(10, 8))
        plt.title(titre + f", DOF {num_dof[0]}")
        plt.xlabel("Frequency (rad/s)")
        plt.ylabel("|Fex| (N)")
        beta_Capytaine=np.array(beta_Capytaine) /np.pi*180
        j=num_dof[0]-1
        for k, beta in enumerate(beta_Capytaine) :
            color = couleurs[k % len(couleurs)]
            marker = markers[k % len(markers)]
            plt.plot(freq_PIT, Coef_PIT_abs[k][:, j], linestyle='-', marker=marker, color=color, markersize=3, label=f"PIT beta= {beta_PIT[k]}°")   
            plt.plot(freq_Capytaine, [row[j] for row in Coef_Capytaine[k]], linestyle=':', marker=marker, color=color, markersize=5, label=f"Capytaine beta= {beta_Capytaine[k]:.1f}°")   
        plt.legend()
        plt.grid(True)
        plt.savefig(f"{file}_BETAS{k+1}_Fex_abs_B{num_dof}.png", dpi=300)
        plt.close()

        plt.figure(figsize=(10, 8))
        plt.title(titre+ f", DOF {num_dof[0]}")
        plt.xlabel("Frequency (rad/s)")
        plt.ylabel("phase(Fex)")
        for k, beta in enumerate(beta_PIT) :
            color = couleurs[k % len(couleurs)]
            marker = markers[k % len(markers)]
            plt.plot(freq_PIT, Coef_PIT_ph[k][:, j], linestyle='-', marker=marker, color=color, markersize=3, label=f"PIT beta= {beta_PIT[k]}°")   
            plt.plot(freq_Capytaine, [row[j] for row in Phase_Capytaine[k]], linestyle=':', marker=marker, color=color, markersize=5, label=f"Capytaine beta= {beta_Capytaine[k]:.1f}°")   
        plt.legend()
        plt.grid(True)
        plt.savefig(f"{file}_BETAS{k+1}_Fex_phase_B{num_dof}.png", dpi=300)
        plt.close()

    else : # ONE wave direction !!
        i=0
        plt.figure(figsize=(10, 8))
        plt.title(titre + f", beta={beta_PIT[i]:.2f}°")
        plt.xlabel("Frequency (rad/s)")
        plt.ylabel("|Fex| (N)")
        nom_dof=""
        for k, dof in enumerate(num_dof) : 
            j=dof-1
            nom_dof=nom_dof+"_"+f"{dof}"
            color = couleurs[k % len(couleurs)]
            marker = markers[k % len(markers)]
            color2 = couleurs[k+len(num_dof) % len(couleurs)]
            marker2 = markers[k+len(num_dof) % len(markers)]
            if REF :
                color3 = couleurs[k+len(num_dof)+1 % len(couleurs)]
                marker3 = markers[k+len(num_dof)+2 % len(markers)]
                plt.plot(freq_PIT_REF, Coef_PIT_REF[i][:,j], linestyle='-', marker=marker3, color=color3, markersize=4, label=f"Outer Cylinder Method dof {j+1}")   
            plt.plot(freq_PIT, Coef_PIT_abs[i][:, j], linestyle=':', marker=marker, color=color, markersize=5, label=f"{comp} with Nemoh dof {j+1}")   
            plt.plot(freq_Capytaine, [row[i][j] for row in Coef_Capytaine], linestyle=':', marker=marker2, color=color2, markersize=4, label=f"{comp} with Capytaine dof {j+1}")   
            plt.legend()
        plt.grid(True)
        plt.savefig(f"{file}_Fex_abs_B{nom_dof}.png", dpi=300)
        plt.close()

        plt.figure(figsize=(10, 8))
        plt.title(titre + f", beta={beta_PIT[i]:.2f}°")
        plt.xlabel("Frequency (rad/s)")
        plt.ylabel("phase(Fex) (N)")
        i=0
        nom_dof=""
        for k, dof in enumerate(num_dof) : 
            j=dof-1
            nom_dof=nom_dof+"_"+f"{dof}"
            color = couleurs[k % len(couleurs)]
            marker = markers[k % len(markers)]
            color2 = couleurs[k+len(num_dof) % len(couleurs)]
            marker2 = markers[k+len(num_dof) % len(markers)]
            if REF :
                color3 = couleurs[k+len(num_dof)+1 % len(couleurs)]
                marker3 = markers[k+len(num_dof)+2 % len(markers)]
                plt.plot(freq_PIT_REF, Phase_PIT_REF[i][:,j], linestyle='-', marker=marker3, color=color3, markersize=4, label=f"Outer Cylinder Method dof {j+1}")   
        
            plt.plot(freq_PIT, Coef_PIT_ph[i][:, j], linestyle=':', marker=marker, color=color, markersize=5, label=f"{comp} with Nemoh dof {j+1}")   
            plt.plot(freq_Capytaine, [row[i][j] for row in Phase_Capytaine], linestyle=':', marker=marker2, color=color2, markersize=4, label=f"{comp} with Capytaine dof {j+1}")   
            plt.legend()
        plt.grid(True)
        plt.savefig(f"{file}_Fex_phase_B{nom_dof}.png", dpi=300)
        plt.close()
    return


#################################################################
#################################################################
######################### PARAMETERS ############################
#################################################################
#################################################################

N3 = False

GRAPHS_beta=False

if GRAPHS_beta or N3 : 
    GRAPHS_N3_F = False
    GRAPHS_N3_C=False
else : 
    GRAPHS_N3_F = True
    GRAPHS_N3_C = True

distance=1
mesh="BargeX6"
LID=""
config_a=False
Nb=1
layout=f"Nb{Nb}_X"
Nw=20
Nbeta=11
depth="_h12"
test=""
beta_value=""
Capytaine_type=f""
chemin="/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/"
type="BEM"
if Nb==1 : 
    # titre=f"{mesh} - Nb={Nb}, Nw={Nw}, Ndof=6, Nbeta={Nbeta}"
    titre=f"BEM Comparison - Barge - Nb={Nb}"
    if type=="S":
        comp="Source Terms Method"
        chemin_Capytaine = f"{chemin}capytaine/MyTestCases/CAP_{mesh}{depth}/resultsBEM/"
        chemin_PIT=f"{chemin}PIT3_E/wec_inputs/PIT3_{mesh}{depth}_source/results/"
        chemin_PIT_REF=f"{chemin}PIT3_E/wec_inputs/PIT3_{mesh}{depth}/results/"
    elif type=='BEM':
        comp="BEM"
        chemin_Capytaine = f"{chemin}capytaine/MyTestCases/C_BEM_{mesh}_Nb1/resultsBEM/Nb1_d0/"
        chemin_PIT=f"{chemin}PFE_Nemoh/MyTestCases/BEM_{mesh}_Nb1/results/"
        chemin_PIT_REF=chemin_PIT
    elif type=="OC":
        comp="Outer Cylinder Method"
        chemin_Capytaine = f"{chemin}capytaine/MyTestCases/CAP_{mesh}_OC/resultsBEM/"
        chemin_PIT=f"{chemin}PIT3_E/wec_inputs/PIT3_{mesh}{depth}/results/"
        chemin_PIT_REF=chemin_PIT
else : 
    titre=f"{mesh} - Nb={Nb}, Nw={Nw}, Ndof=6, Ndir=1, d={distance}m"
    comp="BEM"
    chemin_Capytaine = f"{chemin}capytaine/MyTestCases/C_BEM_{mesh}/resultsBEM/{layout}_d{distance}/"
    chemin_PIT=f"{chemin}PFE_Nemoh/MyTestCases/{mesh}/BEM_{mesh}_{layout}/BEM_{mesh}_{layout}_d{distance}/results/"
    # chemin_PIT=f"{chemin}PFE_Nemoh/MyTestCases/BEM_{mesh}_SYM_{layout}_d{distance}/results/"
    chemin_PIT_REF=chemin_PIT
    # chemin_PIT_REF=f"{chemin}PFE_Nemoh/MyTestCases/BEM_{mesh}_{layout}_d{distance}/results/"
   
CM_PIT = f"{chemin_PIT}CM.dat"
CM_Capytaine_N3 = f"{chemin_Capytaine}Capytaine_Madd.dat"
CM_PIT_REF = f"{chemin_PIT_REF}CM.dat"

CA_PIT = f"{chemin_PIT}CA.dat"
CA_Capytaine_N3 = f"{chemin_Capytaine}Capytaine_Crad.dat"
CA_PIT_REF = f"{chemin_PIT_REF}CA.dat"

Fex_PIT = f"{chemin_PIT}ExcitationForce.tec"  
Fex_Capytaine_N3 = f"{chemin_Capytaine}Capytaine_Fe_abs.dat"
Fex_PIT_REF = f"{chemin_PIT_REF}ExcitationForce.tec"
Fex_phase_Capytaine_N3 = f"{chemin_Capytaine}Capytaine_Fe_phase.dat"
if N3 : 
    param_d=0
    with open(f"EcartsRelatifs_Nb{Nb}.dat", "a") as f:
        L_inf_error, L2_error = CompareResultsRAD(CM_PIT, CM_Capytaine_N3, "CM", test, param_d, f"Nb{Nb}_{layout}", Nb, 0.01, 1) # tol sur erreur, precision sur valeur
        f.write(f"{L_inf_error:16.8f}     {L2_error:16.8f}      ")
        L_inf_error, L2_error = CompareResultsRAD(CA_PIT, CA_Capytaine_N3, "CA", test, param_d, f"Nb{Nb}_{layout}", Nb, 0.01, 1)
        f.write(f"{L_inf_error:16.8f}     {L2_error:16.8f}      ")
        L_inf_error_abs, L2_error_abs, L_inf_error_ph, L2_error_ph = CompareResultsFex(Fex_PIT, Fex_Capytaine_N3, Fex_phase_Capytaine_N3, test, param_d, f"Nb{Nb}_{layout}", Nb, 0.01, 1)
        f.write(f"{L_inf_error_abs:16.8f}     {L2_error_abs:16.8f}      ")
        f.write(f"{L_inf_error_ph:16.8f}      {L2_error_ph:16.8f}\n")

if GRAPHS_N3_F :
    for ind, i in enumerate([[1], [3], [5]]) :
        print("dof", i)
        Trace_Fex(Fex_PIT, Fex_Capytaine_N3, Fex_phase_Capytaine_N3, Fex_PIT_REF, f"BEM_{mesh}_Nb{Nb}", i, False, titre, comp)
if GRAPHS_N3_C : 
    # indices=[(5,5), (11,11), (5,11), (3,9), (3,5), (3,11), (5,9), (9,11)]
    # indices=[[(3,9)], [(1, 1),(3,3)], [(5,5)], [(5,11)], [(1, 7)]]
    indices=[[(1, 1)], [(3,3)], [(5,5)], [(1,5)]]
    for ij in indices : 
        print("indice", ij)
        Trace_RAD(CA_PIT, CA_Capytaine_N3, CA_PIT_REF, f"BEM_{mesh}_Nb{Nb}", "Damping", ij, titre, comp)
        Trace_RAD(CM_PIT, CM_Capytaine_N3, CM_PIT_REF, f"BEM_{mesh}_Nb{Nb}", "Added_Mass", ij, titre, comp)
if GRAPHS_beta : 
    for ind, i in enumerate([[3],[9], [5],[11], [1],[7]])  :
        print("dof", i)
        Trace_Fex(Fex_PIT, Fex_Capytaine_N3, Fex_phase_Capytaine_N3, Fex_PIT_REF, f"BEM_{mesh}_Nb{Nb}", i, True, titre, comp)
            

