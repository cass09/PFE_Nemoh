import numpy as np
import re  # Pour utiliser les expressions régulières
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def read_RAD_BEM(fileRAD):
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

def CompareResultsRAD(fileBEM, filePIT, nom, comment, param_d, type, Nb, tol, precision) :
    freq_BEM, Coef_BEM = read_RAD_BEM(fileBEM)
    freq_PIT, Coef_PIT = read_RAD_PIT(filePIT)
    if np.all(np.isclose(freq_BEM, freq_PIT, atol=1e-4)):
        if np.shape(Coef_BEM) != np.shape(Coef_PIT):
            print("BEM", np.shape(Coef_BEM), "PIT", np.shape(Coef_PIT))
            raise ValueError("tailles différentes")
        else :
            errorL2 = np.sqrt(np.sum((np.array(Coef_BEM) - np.array(Coef_PIT))**2))
            errorINF = np.max(np.abs(np.array(Coef_BEM) - np.array(Coef_PIT)))

            normL2_BEM = np.sqrt(np.sum(np.array(Coef_BEM)**2))
            relative_errorL2 = errorL2 / normL2_BEM if normL2_BEM != 0 else np.nan

            max_BEM = np.max(np.abs(np.array(Coef_BEM)))
            relative_errorINF = errorINF / max_BEM if max_BEM != 0 else np.nan

            All_errors=True
            if All_errors:
                dof=6*Nb
                n_error=0
                for i in range(dof):
                    for j in range(dof):
                        BEM = []
                        PIT = []                       

                        BEM = [0 if np.abs(sublist[i][j]) < precision else sublist[i][j] for sublist in Coef_BEM]
                        PIT = [0 if np.abs(sublist[i][j]) < precision else sublist[i][j] for sublist in Coef_PIT]
                        # Remplacer les valeurs < precision par 0
                        if np.max(np.abs(np.array(BEM))) != 0 and np.max(np.abs(np.array(PIT)))!=0:
                            error=np.max(np.abs(np.array(BEM)-np.array(PIT)))/np.max(np.abs(np.array(BEM)))
                        else :
                            error=0
                        with open(f"EcartsRelatifs_{comment}{nom}_{type}.dat", "a") as f:
                            if i==0 and j==0:
                                f.write(f"\n ")
                            if j==0:
                                f.write(f" \n {param_d}     ")
                            f.write(f"{error:16.8f}     ")
                        if error<=tol:
                            n_error=n_error+1
                print("d=", param_d, "--------> n_error", n_error, "/", dof*dof, "<", tol*100, "%")             
    else :
        print("fréquences différentes")
        print('BEM', freq_BEM)
        print('PIT', freq_PIT)
        relative_errorL2=[]
        relative_errorINF=[]
    
    return relative_errorL2, relative_errorINF

def read_Fex_BEM(fileFex):
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
            if 'Period' in line:
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

def CompareResultsFex(fileBEM, filePIT, filePIT_phase, comment, param_d, type, Nb, tol, precision) :
    beta_BEM, freq_BEM, Coef_BEM_abs, Coef_BEM_ph = read_Fex_BEM(fileBEM)
    beta_PIT, freq_PIT, Coef_PIT = read_Fex_PIT(filePIT)
    beta_PIT, freq_PIT, Phase_PIT = read_Fex_PIT(filePIT_phase)
    # print("BEM", Coef_BEM_ph)
    # print("PIT", Phase_PIT)

    if len(freq_BEM)==len(freq_PIT) and len(beta_BEM)==len(beta_PIT):
        if np.shape(Coef_BEM_abs) != np.shape(Coef_PIT):
            print(np.shape(Coef_BEM_abs), np.shape(Coef_PIT))
            raise ValueError("tailles différentes")
        if np.shape(Coef_BEM_ph) != np.shape(Phase_PIT):
            print(np.shape(Coef_BEM_ph), np.shape(Phase_PIT))
            raise ValueError("tailles différentes")
        else :
            errorL2 = np.sqrt(np.sum((np.array(Coef_BEM_abs) - np.array(Coef_PIT))**2))
            # errorINF = np.max(np.abs(np.array(Coef_BEM_abs) - np.array(Coef_PIT)))
            erreurs = np.abs(np.array(Coef_BEM_abs) - np.array(Coef_PIT))
            errorINF = np.max(erreurs)
            indice_max = np.argmax(erreurs)
            # print("indice_max =", indice_max)
            # print("Fex_abs_BEM_max =", np.array(Coef_BEM_abs)[indice_max])
            # print("Fex_abs_PIT_max =", np.array(Coef_PIT)[indice_max])

            normL2_BEM = np.sqrt(np.sum(np.array(Coef_BEM_abs)**2))
            relative_errorL2_abs = errorL2 / normL2_BEM if normL2_BEM != 0 else np.nan

            max_BEM = np.max(np.abs(np.array(Coef_BEM_abs)))
            relative_errorINF_abs = errorINF / max_BEM if max_BEM != 0 else np.nan

            # PHASE 
            errorL2_p = np.sqrt(np.sum((np.array(Coef_BEM_ph) - np.array(Phase_PIT))**2))
            # errorINF_p = np.max(np.abs(np.array(Coef_BEM_ph) - np.array(Phase_PIT)))

            erreurs_p = np.abs(np.array(Coef_BEM_ph) - np.array(Phase_PIT))
            errorINF_p = np.max(erreurs_p)
            # indice_max_p = np.argmax(erreurs_p)
            # print("indice_max =", indice_max_p, Coef_BEM_ph[0][6][1], Phase_PIT[0][6][1])

            normL2_BEM_p = np.sqrt(np.sum(np.array(Coef_BEM_ph)**2))
            relative_errorL2_phase = errorL2_p / normL2_BEM_p if normL2_BEM_p != 0 else np.nan

            max_BEM_p = np.max(np.abs(np.array(Coef_BEM_ph)))
            relative_errorINF_phase = errorINF_p / max_BEM_p if max_BEM_p != 0 else np.nan

            # Calcul erreur pour chaque dof 
            All_errors=True
            if All_errors:
                dof=6*Nb
                n_err_abs=0
                n_err_ph=0
                for i in range(dof):
                    B_abs = [N[i] for N in Coef_BEM_abs[0]]
                    P_abs = [N[i] for N in Coef_PIT[0]]
                    B_ph = [N[i] for N in Coef_BEM_ph[0]]
                    P_ph = [N[i] for N in Phase_PIT[0]]
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
        print('BEM', freq_BEM)
        print('PIT', freq_PIT)
        relative_errorL2_abs=[]
        relative_errorL2_phase=[]
        relative_errorINF_abs=[]
        relative_errorINF_phase=[]
    
    return relative_errorL2_abs, relative_errorINF_abs, relative_errorL2_phase, relative_errorINF_phase


def Trace_RAD(fileBEM, filePIT, filePIT_REF, file, nom, indice, titre, d):
    freq_BEM, Coef_BEM = read_RAD_BEM(fileBEM)
    freq_PIT, Coef_PIT = read_RAD_PIT(filePIT)
    REF=False
    if filePIT_REF != filePIT :
        freq_PIT_REF, Coef_PIT_REF = read_RAD_PIT(filePIT_REF)
        REF=True
    couleurs = ['b', 'g', 'r', 'gold', 'lime', 'c', 'm', 'y', 'k']  # palette de couleurs (réutilisée si plus de 7 courbes)
    markers = ['s', 'o', '<', '^','x', '+', '*']
    AD=False
    if AD :
        data=np.loadtxt("Data_waves.dat", delimiter=None, skiprows=1) 
        a=5
        rho=1000
        h=2*a
        norm=rho*np.pi*a**2*h/3
        if nom=="Damping":
            norm=norm*np.array(freq_BEM)
            unit="/(ρπa²hw)"
        else :
            unit="/(ρπa²h)"
        bis="AD_"
        axeX="ka (a=5m)"
        if np.isscalar(norm):
            norm_array = np.full_like(freq_BEM, norm, dtype=float)
        else:
            norm_array = np.asarray(norm, dtype=float)
    else :
        axeX="Frequency (rad/s)"
        bis=""
        if nom=="Damping":
            unit="(kg/s)"
        else :
            unit="(kg)"

    plt.figure(figsize=(8, 6))
    plt.title(titre)
    plt.xlabel(f"{axeX}")
    plt.ylabel(f"{nom} {unit}")
    plt.grid(True)

    for k, ij in enumerate(indice): 
        i=ij[0]-1
        j=ij[1]-1
        color = couleurs[k % len(couleurs)]
        marker = markers[k % len(markers)]
        color2 = couleurs[k+len(indice) % len(couleurs)]
        marker2 = markers[k+len(indice) % len(markers)]
        if AD :
            plt.plot(a*data[:,1], [mat[i][j]/ norm_array[c] for c, mat in enumerate(Coef_BEM)], linestyle='-', marker=marker, color=color, markersize=5, label=f"BEM {ij[0]}_{ij[1]}")
            plt.plot(a*data[:,1], [mat[i][j]/ norm_array[c] for c, mat in enumerate(Coef_PIT)], linestyle=':', marker=marker2, color=color2, markersize=5, label=f"ITM {ij[0]}_{ij[1]}")
        else :
            plt.plot(freq_BEM, [mat[i][j] for mat in Coef_BEM], linestyle='-', marker="o", color=color, markersize=3, label=f"BEM {ij[0]}_{ij[1]}")
            plt.plot(freq_PIT, [mat[i][j] for mat in Coef_PIT], linestyle=':', marker="x", color=color, markersize=6, label=f"STM {ij[0]}_{ij[1]}")
        
        if REF :
            color3 = couleurs[k+len(indice)+1 % len(couleurs)]
            marker3 = markers[k+len(indice)+2 % len(markers)]
            plt.plot(freq_PIT_REF, [mat[i][j] for mat in Coef_PIT_REF], linestyle=':', marker=marker3, color='black', markersize=3, label=f"OCM {ij[0]}_{ij[1]}")
    plt.legend(loc='best', fontsize='small', frameon=True)
    plt.savefig(f"{file}_{bis}{nom}_{d}_M{indice}.pdf")
    plt.close()
    return

def Trace_Fex(fileBEM, filePIT, filePIT_phase, filePIT_REF, filePIT_phase_REF, file, num_dof, beta, titre, d):
    beta_BEM, freq_BEM, Coef_BEM_abs, Coef_BEM_ph = read_Fex_BEM(fileBEM)
    beta_PIT, freq_PIT, Coef_PIT = read_Fex_PIT(filePIT)
    beta_PIT, freq_PIT, Phase_PIT = read_Fex_PIT(filePIT_phase)
    REF=False
    if filePIT_REF != filePIT :
        beta_PIT_REF, freq_PIT_REF, Coef_PIT_REF = read_Fex_PIT(filePIT_REF)
        beta_PIT_REF, freq_PIT_REF, Phase_PIT_REF = read_Fex_PIT(filePIT_phase_REF)
        REF=True
    couleurs = ['b', 'r', 'm','gold', 'lime', 'g', 'c', 'm', 'y', 'k']  # palette de couleurs (réutilisée si plus de 7 courbes)
    markers = ['s', 'o', '<', '^','x', '+', '*']
    if beta : # Several wave directions
        plt.figure(figsize=(10, 8))
        plt.title(titre + f", DOF {num_dof[0]}")
        plt.xlabel("Frequency (rad/s)")
        plt.ylabel("|Fex| (N)")
        beta_PIT=np.array(beta_PIT) /np.pi*180
        j=num_dof[0]-1
        for k, beta in enumerate(beta_PIT) :
            color = couleurs[k % len(couleurs)]
            marker = markers[k % len(markers)]
            plt.plot(freq_BEM, Coef_BEM_abs[k][:, j], linestyle='-', marker=marker, color=color, markersize=3, label=f"BEM beta= {beta_BEM[k]}°")   
            plt.plot(freq_PIT, [row[j] for row in Coef_PIT[k]], linestyle=':', marker=marker, color=color, markersize=5, label=f"PIT beta= {beta_PIT[k]:.1f}°")   
        plt.legend()
        plt.grid(True)
        plt.savefig(f"{file}_BETAS{k+1}_Fex_abs_{d}_B{num_dof}.pdf")
        plt.close()

        plt.figure(figsize=(10, 8))
        plt.title(titre+ f", DOF {num_dof[0]}")
        plt.xlabel("Frequency (rad/s)")
        plt.ylabel("phase(Fex)")
        for k, beta in enumerate(beta_BEM) :
            color = couleurs[k % len(couleurs)]
            marker = markers[k % len(markers)]
            plt.plot(freq_BEM, Coef_BEM_ph[k][:, j], linestyle='-', marker=marker, color=color, markersize=3, label=f"BEM beta= {beta_BEM[k]}°")   
            plt.plot(freq_PIT, [row[j] for row in Phase_PIT[k]], linestyle=':', marker=marker, color=color, markersize=5, label=f"PIT beta= {beta_PIT[k]:.1f}°")   
        plt.legend()
        plt.grid(True)
        plt.savefig(f"{file}_BETAS{k+1}_Fex_phase_{d}_B{num_dof}.png", dpi=300)
        plt.close()

    else : # ONE wave direction !!
        i=0
        plt.figure(figsize=(10, 8))
        plt.title(titre + f", beta={beta_BEM[i]:.2f}°")
        plt.xlabel("Frequency (rad/s)")
        plt.ylabel("|Fex| (N)")
        nom_dof=""
        for k, dof in enumerate(num_dof) : 
            j=dof-1
            nom_dof=nom_dof+"_"+f"{dof}"
            color = couleurs[k % len(couleurs)]
            marker = markers[k % len(markers)]
            plt.plot(freq_BEM, Coef_BEM_abs[i][:, j], linestyle='-', marker="o", color=color, markersize=3, label=f"BEM dof {j+1}")   
            plt.plot(freq_PIT, [row[j] for row in Coef_PIT[i]], linestyle='--', marker="x", color=color, markersize=6, label=f"STM dof {j+1}")   
            if REF :
                color3 = couleurs[k+len(num_dof)+1 % len(couleurs)]
                marker3 = markers[k+len(num_dof)+2 % len(markers)]
                plt.plot(freq_PIT_REF, [row[j] for row in Coef_PIT_REF[i]], linestyle=':', marker=marker3, color='black', markersize=5, label=f"OCM dof {j+1}")   
        plt.legend()
        plt.grid(True)
        plt.savefig(f"{file}_Fex_abs_{d}_B{nom_dof}.pdf")
        plt.close()

        plt.figure(figsize=(10, 8))
        plt.title(titre + f", beta={beta_BEM[i]:.2f}°")
        plt.xlabel("Frequency (rad/s)")
        plt.ylabel("phase(Fex) (N)")
        i=0
        nom_dof=""
        for k, dof in enumerate(num_dof) : 
            j=dof-1
            nom_dof=nom_dof+"_"+f"{dof}"
            color = couleurs[k % len(couleurs)]
            marker = markers[k % len(markers)]
            plt.plot(freq_BEM, Coef_BEM_ph[i][:, j], linestyle='-', marker="o", color=color, markersize=3, label=f"BEM dof {j+1}")   
            plt.plot(freq_PIT, [row[j] for row in Phase_PIT[i]], linestyle='--', marker="x", color=color, markersize=6, label=f"STM dof {j+1}")   
            if REF :
                color3 = couleurs[k+len(num_dof)+1 % len(couleurs)]
                marker3 = markers[k+len(num_dof)+2 % len(markers)]
                plt.plot(freq_PIT_REF, [row[j] for row in Phase_PIT_REF[i]], linestyle=':', marker=marker3, color='black', markersize=5, label=f"OCM dof {j+1}")   
        plt.legend()
        plt.grid(True)
        plt.savefig(f"{file}_Fex_phase_{d}_B{nom_dof}.pdf")
        plt.close()
    return

def plot_RAD_E(fileBEM, filePIT_REF, PIT_files, nom, indice, distance, titre, Ne):
    freq_BEM, Coef_BEM = read_RAD_BEM(fileBEM)
    freq_PIT_REF, Coef_PIT_REF = read_RAD_PIT(filePIT_REF)
    couleurs = ['g', 'gold', 'm', 'c', 'lime', 'y', 'k']  # palette de couleurs (réutilisée si plus de 7 courbes)
    markers = ['s', '^','x','<', '+', '*']
    if nom=="Damping":
        unit="(kg/s)"
    else :
        unit="(kg)"
    for k, ij in enumerate(indice): 
        plt.figure(figsize=(8, 6))
        plt.title(titre + f", Coefficient {ij[0]}-{ij[1]}")
        plt.xlabel("Frequency (rad/s)")
        plt.ylabel(f"{nom} {unit}")
        plt.grid(True)
        i=ij[0]-1
        j=ij[1]-1
        marker = markers[k % len(markers)]
        marker2 = markers[k+1 % len(markers)]
        plt.plot(freq_BEM, [mat[i][j] for mat in Coef_BEM], linestyle='-', marker=marker, color='b', markersize=4, label=f"BEM")
        plt.plot(freq_PIT_REF, [mat[i][j] for mat in Coef_PIT_REF], linestyle='--', marker=marker2, color='r', markersize=4, label=f"PIT Ne=0")
        for idx, file in enumerate(PIT_files) :
            color = couleurs[k+idx % len(couleurs)]
            freq_PIT, Coef_PIT = read_RAD_PIT(file)
            plt.plot(freq_PIT, [mat[i][j] for mat in Coef_PIT], linestyle=':', marker="o", color=color, markersize=2, label=f"PIT - Ne={Ne[idx]}")
        plt.legend(loc='best', fontsize='small', frameon=True)
        plt.savefig(f"N3_E_{nom}_d{distance}_M{ij}.png", dpi=300)
        plt.close()
    return

def plot_Fex_E(fileBEM, filePIT_REF, filePIT_phase_REF, filesPIT, filesPIT_phase, num_dof, distance, titre, Ne):
    beta_BEM, freq_BEM, Coef_BEM_abs, Coef_BEM_ph = read_Fex_BEM(fileBEM)
    beta_PIT_REF, freq_PIT_REF, Coef_PIT_REF = read_Fex_PIT(filePIT_REF)
    beta_PIT_REF, freq_PIT_REF, Phase_PIT_REF = read_Fex_PIT(filePIT_phase_REF)
    couleurs = ['g', 'gold', 'm', 'c', 'lime', 'y', 'k']  # palette de couleurs (réutilisée si plus de 7 courbes)
    markers = ['s', '<', '^','x', '+', '*']
    i=0
    
    for k, dof in enumerate(num_dof) : 
        plt.figure(figsize=(10, 8))
        plt.title(titre + f", beta={beta_BEM[i]:.2f}°, DOF {dof}")
        plt.xlabel("Frequency (rad/s)")
        plt.ylabel("|Fex| (N)")
        j=dof-1
        marker = markers[k % len(markers)]
        marker2 = markers[k+1 % len(markers)]
        plt.plot(freq_BEM, Coef_BEM_abs[i][:, j], linestyle='-', marker=marker, color='b', markersize=4, label=f"BEM ")   
        plt.plot(freq_PIT_REF, [row[j] for row in Coef_PIT_REF[i]], linestyle='--', marker=marker2, color='red', markersize=4, label=f"PIT - Ne=0")   
        for idx, file in enumerate(filesPIT) :
            color = couleurs[k+idx % len(couleurs)]
            beta_PIT, freq_PIT, Coef_PIT = read_Fex_PIT(file)
            plt.plot(freq_PIT, [row[j] for row in Coef_PIT[i]], linestyle=':', marker='o', color=color, markersize=2, label=f"PIT - Ne={Ne[idx]}")   
        plt.legend()
        plt.grid(True)
        plt.savefig(f"N3_E_Fex_abs_d{distance}_B{dof}.png", dpi=300)
        plt.close()

    
    i=0
    for k, dof in enumerate(num_dof) : 
        plt.figure(figsize=(10, 8))
        plt.title(titre + f", beta={beta_BEM[i]:.2f}°, DOF {dof}")
        plt.xlabel("Frequency (rad/s)")
        plt.ylabel("phase(Fex) (N)")
        j=dof-1
        marker = markers[k % len(markers)]
        marker2 = markers[k+1 % len(markers)]
        plt.plot(freq_BEM, Coef_BEM_ph[i][:, j], linestyle='-', marker=marker, color='b', markersize=4, label=f"BEM")   
        plt.plot(freq_PIT_REF, [row[j] for row in Phase_PIT_REF[i]], linestyle='--', marker=marker2, color='red', markersize=4, label=f"PIT -Ne=0")   
        for idx, file in enumerate(filesPIT_phase) :
            color = couleurs[k+idx % len(couleurs)]
            beta_PIT, freq_PIT, Phase_PIT = read_Fex_PIT(file)
            plt.plot(freq_PIT, [row[j] for row in Phase_PIT[i]], linestyle=':', marker='o', color=color, markersize=2, label=f"PIT - Ne={Ne[idx]}")   
        plt.legend()
        plt.grid(True)
        plt.savefig(f"N3_E_Fex_phase_d{distance}_B{dof}.png", dpi=300)
        plt.close()
    return

def Calcul_GOF(BEMfiles, PITfiles, Nom_data, indices, NOMS, mesh):
    min_CoF = np.zeros((len(indices), len(PITfiles)))
    mean_CoF = np.zeros((len(indices), len(PITfiles)))
    use_nrmse=True
    for k in range(len(BEMfiles)):
        frequencies, COEF = read_RAD_PIT(PITfiles[k])
        frequencies, COEF_REF = read_RAD_BEM(BEMfiles[k])
        for j, ij in enumerate(indices):
            diffs_squared = []
            ref_values = []
            gof_local = []
            for i in range(len(frequencies)):
                if frequencies[i]<4:
                    value = COEF[i][ij[0] - 1][ij[1] - 1]
                    value_ref = COEF_REF[i][ij[0] - 1][ij[1] - 1]

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
            if gof < 90:
                print(f"mesh={NOMS[k]} - ij={ij} - mean GoF={gof:.2f}% ")
            if gof < 0:
                print(f" !!!!!! mesh={NOMS[k]} - ij={ij} - mean GoF={gof:.2f}%")
                print(f"{len(diffs_squared)} - diffs_squared={diffs_squared}")
    
    # Création du graphique
    for j, ij in enumerate(indices):
        plt.figure(figsize=(8, 5))
        bars = plt.bar(NOMS[:], mean_CoF[j,:], color='skyblue', edgecolor='black')

        # Ajouter les valeurs au-dessus des barres
        for bar in bars:
            yval = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2, yval + 0.5, f'{yval:.2f}%', ha='center', va='bottom')

        # Mise en forme
        plt.ylim(0, 105)
        plt.xlabel('d/a')
        plt.ylabel(f'mean GoF (%)')
        if Nom_data=="Damping":
            plt.title(f'Goodness of Fit - Two Cylinders - {Nom_data} B_{ij[0]}_{ij[1]}')
        else : 
            plt.title(f'Goodness of Fit - Two Cylinders - Added Mass A_{ij[0]}_{ij[1]}')
        plt.grid(axis='y', linestyle='--', alpha=0.7)

        plt.tight_layout()
        # plt.show()
        plt.savefig(f"GOF_{mesh}_{Nom_data}_{ij[0]}_{ij[1]}.pdf")
        plt.close()

def Calcul_GOF_Fex(BEMfiles, PITfilesA, PITfilesP, dof_vect, NOMS, mesh):
    min_CoF=np.zeros((len(dof_vect), len(PITfilesA)))
    mean_CoF=np.zeros((len(dof_vect), len(PITfilesA)))
    min_CoF_ph=np.zeros((len(dof_vect), len(PITfilesP)))
    mean_CoF_ph=np.zeros((len(dof_vect), len(PITfilesP)))
    for k in range(len(BEMfiles)) : 
        beta, frequencies, Coef_abs_REF, Coef_ph_REF = read_Fex_BEM(BEMfiles[k])
        beta, frequencies, Coef_abs = read_Fex_PIT(PITfilesA[k])
        beta, frequencies, Coef_ph = read_Fex_PIT(PITfilesP[k])
        use_nrmse = True  # ⬅️ Change à False pour utiliser la méthode locale

        for j, dof in enumerate(dof_vect):
            abs_diffs_squared = []
            abs_refs = []
            ph_diffs_squared = []
            ph_refs = []

            gof_local_abs = []
            gof_local_ph = []

            for i in range(len(frequencies)):
                # --- Amplitude ---
                value = Coef_abs[0][i][dof - 1]
                value_ref = Coef_abs_REF[0][i, dof - 1]
                diff_sq_abs = np.abs(value - value_ref)**2
                ref_abs = np.abs(value_ref)

                if frequencies[i]<4:
                    abs_diffs_squared.append(diff_sq_abs)
                    abs_refs.append(ref_abs)
                    if not use_nrmse:
                        gof_i_abs = 100 * (1 - (np.sqrt(diff_sq_abs) / ref_abs))
                        gof_local_abs.append(gof_i_abs)

                # --- Phase ---
                value_ph = Coef_ph[0][i][dof - 1]
                value_ref_ph = Coef_ph_REF[0][i, dof - 1]
                diff_sq_ph = np.abs(value_ph - value_ref_ph)**2
                ref_ph = np.abs(value_ref_ph)

                if frequencies[i]<4:
                    ph_diffs_squared.append(diff_sq_ph)
                    ph_refs.append(ref_ph)
                    if not use_nrmse:
                        gof_i_ph = 100 * (1 - (np.sqrt(diff_sq_ph) / ref_ph))
                        gof_local_ph.append(gof_i_ph)

            # --- Amplitude GoF ---
            if use_nrmse:
                rmse_abs = np.sqrt(np.mean(abs_diffs_squared))
                mean_abs_ref = np.mean(abs_refs)
                nrmse_abs = rmse_abs / mean_abs_ref
                gof_abs = (1 - nrmse_abs) * 100
            else:
                gof_abs = np.mean(gof_local_abs)

            mean_CoF[j, k] = gof_abs

            # --- Phase GoF ---
            if use_nrmse:
                rmse_ph = np.sqrt(np.mean(ph_diffs_squared))
                mean_ph_ref = np.mean(ph_refs)
                nrmse_ph = rmse_ph / mean_ph_ref
                gof_ph = (1 - nrmse_ph) * 100
            else:
                gof_ph = np.mean(gof_local_ph)

            mean_CoF_ph[j, k] = gof_ph

            # Optionnel : avertissement si GoF < 90%
            if gof_abs < 90:
                print(f"dof={dof} - mesh={NOMS[k]} - ABS mean GoF={gof_abs:.2f}%")
            if gof_ph < 90:
                print(f"dof={dof} - mesh={NOMS[k]} - PHASE mean GoF={gof_ph:.2f}%")
    Trace_GOF(mean_CoF, False, NOMS, dof_vect, mesh)
    Trace_GOF(mean_CoF_ph, True, NOMS, dof_vect, mesh)
    
def Trace_GOF(GoF, Phase, NOMS, dof_vect, mesh) :
    for j, dof in enumerate(dof_vect):
        plt.figure(figsize=(8, 5))
        bars = plt.bar(NOMS[:], GoF[j,:], color='skyblue', edgecolor='black')

        # Ajouter les valeurs au-dessus des barres
        for bar in bars:
            yval = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2, yval + 0.5, f'{yval:.2f}%', ha='center', va='bottom')

        # Mise en forme
        plt.ylim(0, 105)
        plt.xlabel('d/a')
        plt.ylabel(f'mean GoF (%)')
        if Phase:
            data="Fph"
            plt.title(f'Goodness of Fit - Two Cylinders - Excitation Forces Phase ph(F)_{dof}')
        else : 
            data="Fabs"
            plt.title(f'Goodness of Fit - Two Cylinders - Excitation Forces Amplitude |F|_{dof}')
        plt.grid(axis='y', linestyle='--', alpha=0.7)

        plt.tight_layout()
        # plt.show()
        plt.savefig(f"GOF_{mesh}_{data}_dof{dof}.pdf")
        plt.close()

#################################################################
#################################################################
######################### PARAMETERS ############################
#################################################################
#################################################################

param_d=4

N3 = False
COMP_E = False
GRAPHS_E=False
ERROR_E=False

GRAPHS_beta=False

GOF=False
if N3 or GOF: 
    param_d=2
    limite=16
    GRAPHS_N3_F = False
    GRAPHS_N3_C=False
elif GRAPHS_beta or COMP_E : 
    limite=param_d
    GRAPHS_N3_F = False
    GRAPHS_N3_C=False
else : 
    limite=param_d
    GRAPHS_N3_F = True
    GRAPHS_N3_C = True


mesh="CylR3"
LID=""
config_a=True
Nb=2
layout="X"
Nw=20
Ndir=1
test="S_"
# beta_value="beta_180.0_"
beta_value=""
PIT_type=f"Nb{Nb}_{layout}_"
# PIT_type=""
PIT_N3_file=f"PIT3_{mesh}{LID}_source"
# PIT_N3_file=f"CAP_test"
PIT_N3_file_REF=f"PIT3_{mesh}{LID}"
chemin="/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/"
chemin_PIT_N3 = f"{chemin}PIT3_E/wec_inputs/{PIT_N3_file}/resultsIT/"
chemin_PIT_N3_REF = f"{chemin}PIT3_E/wec_inputs/{PIT_N3_file_REF}/resultsIT/"
# BEM_file=f"BEM_{mesh}_BETAS"
# chemin_BEM=f"BEM_{mesh}/BEM_{mesh}_Nb{Nb}_{layout}_"
# chemin_BEM=f"{mesh}{LID}/BEM_{mesh}_Nb{Nb}_{layout}/BEM_{mesh}_Nb{Nb}_{layout}_"
chemin_BEM=f"{mesh}{LID}/BEM_{mesh}_Nb{Nb}_{layout}_"
# titre=f"{mesh} - Nb={Nb} - {layout}, Nw={Nw}, Nbeta=13, Ndof=6, Ndir={Ndir}"
titre = "Four Cylinders"
if GRAPHS_beta : 
    Ndir=5
    beta_value=f"multi_beta_{Ndir}_"
    PIT_type=""


BEM_A=[]
BEM_B=[]
BEM_F=[]
PIT_A=[]
PIT_B=[]
PIT_Fa=[]
PIT_Fp=[]
NOMS=[]
while param_d<=limite : 
    if config_a : 
        distance=f"da{param_d}"  
        titre=f"{titre} - d/a={param_d}"
    else :
        distance=f"d{param_d}"  
        titre=f"{titre}, d=a+Rc={param_d}m"
    BEM_file=f"{chemin_BEM}da{param_d}"
    CM_BEM = f"{chemin}PFE_Nemoh/MyTestCases/{BEM_file}/results/CM.dat"
    CM_PIT_N3 = f"{chemin_PIT_N3}Global_{test}Madd_{PIT_type}{beta_value}{distance}.00.dat"
    CM_PIT_N3_REF = f"{chemin_PIT_N3_REF}Global_Madd_{PIT_type}{beta_value}{distance}.00.dat"

    CA_BEM = f"{chemin}PFE_Nemoh/MyTestCases/{BEM_file}/results/CA.dat"
    CA_PIT_N3 = f"{chemin_PIT_N3}Global_{test}Crad_{PIT_type}{beta_value}{distance}.00.dat"
    CA_PIT_N3_REF = f"{chemin_PIT_N3_REF}Global_Crad_{PIT_type}{beta_value}{distance}.00.dat"

    Fex_BEM = f"{chemin}PFE_Nemoh/MyTestCases/{BEM_file}/results/ExcitationForce.tec"  
    Fex_PIT_N3 = f"{chemin_PIT_N3}Global_{test}Fe_abs_{PIT_type}{beta_value}{distance}.00.dat"
    Fex_PIT_N3_REF = f"{chemin_PIT_N3_REF}Global_Fe_abs_{PIT_type}{beta_value}{distance}.00.dat"
    Fex_phase_PIT_N3 = f"{chemin_PIT_N3}Global_{test}Fe_phase_{PIT_type}{beta_value}{distance}.00.dat"
    Fex_phase_PIT_N3_REF = f"{chemin_PIT_N3_REF}Global_Fe_phase_{PIT_type}{beta_value}{distance}.00.dat"
    if GOF : 
        BEM_A.append(CM_BEM)
        BEM_B.append(CA_BEM)
        BEM_F.append(Fex_BEM)
        PIT_A.append(CM_PIT_N3)
        PIT_B.append(CA_PIT_N3)
        PIT_Fa.append(Fex_PIT_N3)
        PIT_Fp.append(Fex_phase_PIT_N3)
        NOMS.append(f"d/a={param_d}")
    if N3 : 
        print("\n ------------Distance=", param_d, "-----------------")
        with open(f"EcartsRelatifs_{test}Nb{Nb}_{layout}.dat", "a") as f:
            f.write(f"{param_d}     ")
            L_inf_error, L2_error = CompareResultsRAD(CM_BEM, CM_PIT_N3, "CM", test, param_d, f"Nb{Nb}_{layout}", Nb, 0.05, 10) # tol sur erreur, precision sur valeur
            f.write(f"{L_inf_error:16.8f}     {L2_error:16.8f}      ")
            L_inf_error, L2_error = CompareResultsRAD(CA_BEM, CA_PIT_N3, "CA", test, param_d, f"Nb{Nb}_{layout}", Nb, 0.05, 10)
            f.write(f"{L_inf_error:16.8f}     {L2_error:16.8f}      ")
            L_inf_error_abs, L2_error_abs, L_inf_error_ph, L2_error_ph = CompareResultsFex(Fex_BEM, Fex_PIT_N3, Fex_phase_PIT_N3, test, param_d, f"Nb{Nb}_{layout}", Nb, 0.01, 10)
            f.write(f"{L_inf_error_abs:16.8f}     {L2_error_abs:16.8f}      ")
            f.write(f"{L_inf_error_ph:16.8f}      {L2_error_ph:16.8f}\n")

    if GRAPHS_N3_F :
        for ind, i in enumerate([[1, 7], [3, 9], [5,11]]) :
            print("dof", i)
            Trace_Fex(Fex_BEM, Fex_PIT_N3, Fex_phase_PIT_N3, Fex_PIT_N3_REF, Fex_phase_PIT_N3_REF, f"N_{mesh}_{test}Nb{Nb}_{layout}", i, False, titre, distance)
    if GRAPHS_N3_C : 
        # indices=[(5,5), (11,11), (5,11), (3,9), (3,5), (3,11), (5,9), (9,11)]
        # indices=[[(1, 7)], [(3, 9), (3, 15), (3, 21)], [(5,11), (5, 17), (5, 23)]]
        indices=[[(1, 1)], [(1, 7)], [(3, 3)], [(3,9)], [(5, 5)], [(5,11)]]
        for ij in indices : 
            print("indice", ij)
            Trace_RAD(CA_BEM, CA_PIT_N3, CA_PIT_N3_REF, f"N_{mesh}_{test}Nb{Nb}_{layout}", "Damping", ij, titre, distance)
            Trace_RAD(CM_BEM, CM_PIT_N3, CM_PIT_N3_REF, f"N_{mesh}_{test}Nb{Nb}_{layout}", "Added_Mass", ij, titre, distance)
    if GRAPHS_beta : 
        for ind, i in enumerate([[3],[9], [5],[11], [1],[7]])  :
            print("dof", i)
            Trace_Fex(Fex_BEM, Fex_PIT_N3, Fex_phase_PIT_N3, Fex_PIT_N3_REF, Fex_phase_PIT_N3_REF, f"N_{mesh}_{test}Nb{Nb}_{layout}", i, True, titre, distance)
            
    param_d=param_d*2

if GOF : 
    dof=[1, 3, 5, 7, 9, 11]
    indices=[[1,7],[3,9],[5,11]]
    print(f"------------GoF -> {NOMS}")
    print("Added Mass")
    Calcul_GOF(BEM_A, PIT_A, "Added_Mass", indices, NOMS, mesh)
    print("Damping")
    Calcul_GOF(BEM_B, PIT_B, "Damping", indices, NOMS, mesh)
    print("Fex")
    Calcul_GOF_Fex(BEM_F, PIT_Fa, PIT_Fp, dof, NOMS, mesh)
    
#################################################################
#################################################################
################### EVANESCENT PROBLEM ########################
#################################################################
#################################################################


if COMP_E :   
    if GRAPHS_E : 
        print("\n ------------ PLOT -----------------")
        # dist_graphs= [[1], [8], [16], [64]]
        vect_d=[1]
        Ne=[1,6]
        DOF=[3]
        # indice=[(3,3), (3,9), (1, 7), (5, 11)]
        indice=[(3,9)]
        # for j, vect_d in enumerate(dist_graphs): 
        print("\n Ne=", Ne)
        print("\n DOF=", DOF )
        print("\n indice=", indice )

        for c, distance in enumerate(vect_d):
            print("\n  Distance=", distance)
            BEM_file=f"BEM_{mesh}_Nb{Nb}_{layout}/BEM_{mesh}_Nb{Nb}_{layout}_d{distance}"
            titre=f"{mesh} - Nb={Nb} {layout}, Nw=20, Nbeta=11, Ndof=6, Ndir=1, d={distance}m"

            CM_BEM = f"{chemin}PFE_Nemoh/MyTestCases/{BEM_file}/results/CM.dat"
            CM_PIT_N3_REF = f"{chemin_PIT_N3}Global_Madd_{PIT_type}{beta_value}d{distance}.00.dat"
            CA_BEM = f"{chemin}PFE_Nemoh/MyTestCases/{BEM_file}/results/CA.dat"
            CA_PIT_N3_REF = f"{chemin_PIT_N3}Global_Crad_{PIT_type}{beta_value}d{distance}.00.dat"

            Fex_BEM = f"{chemin}PFE_Nemoh/MyTestCases/{BEM_file}/results/ExcitationForce.tec"  
            Fex_PIT_N3_REF = f"{chemin_PIT_N3}Global_Fe_abs_{PIT_type}{beta_value}d{distance}.00.dat"
            Fex_phase_PIT_N3_REF = f"{chemin_PIT_N3}Global_Fe_phase_{PIT_type}{beta_value}d{distance}.00.dat"

            PIT_files_CM=[]
            PIT_files_CA=[]
            PIT_files_Fex=[]
            PIT_files_Fex_ph=[]
            
            for i, E in enumerate(Ne) : 
                CM_PIT_N3 = f"{chemin_PIT_N3}Global_E{E}_Madd_{PIT_type}{beta_value}d{distance}.00.dat"
                CA_PIT_N3 = f"{chemin_PIT_N3}Global_E{E}_Crad_{PIT_type}{beta_value}d{distance}.00.dat"
                Fex_PIT_N3 = f"{chemin_PIT_N3}Global_E{E}_Fe_abs_{PIT_type}{beta_value}d{distance}.00.dat"
                Fex_phase_PIT_N3 = f"{chemin_PIT_N3}Global_E{E}_Fe_phase_{PIT_type}{beta_value}d{distance}.00.dat"
                PIT_files_CM.append(CM_PIT_N3) 
                PIT_files_CA.append(CA_PIT_N3) 
                PIT_files_Fex.append(Fex_PIT_N3) 
                PIT_files_Fex_ph.append(Fex_phase_PIT_N3) 
            plot_RAD_E(CM_BEM, CM_PIT_N3_REF, PIT_files_CM, "Added_Mass", indice, distance, titre, Ne)
            plot_RAD_E(CA_BEM, CA_PIT_N3_REF, PIT_files_CA, "Damping", indice, distance, titre, Ne)
            plot_Fex_E(Fex_BEM, Fex_PIT_N3_REF, Fex_phase_PIT_N3_REF, PIT_files_Fex, PIT_files_Fex_ph, DOF, distance, titre, Ne)
    
    if ERROR_E :
        print("\n ------------ ERROR -----------------")
        vect_d=[]
        for k in range(0, 5):
            vect_d.append(2**k) 
        Ne=[1, 6, 12]
        DOF=[1, 3, 5]
        indice=[(3,9), (1, 7)]
        # for j, vect_d in enumerate(dist_graphs): 
        print("\n Ne=", Ne)
        print("\n DOF=", DOF )
        print("\n indice=", indice )

        for c, distance in enumerate(vect_d):
            print("\n  Distance=", distance)
            BEM_file=f"BEM_{mesh}_Nb{Nb}_{layout}/BEM_{mesh}_Nb{Nb}_{layout}_d{distance}"

            CM_BEM = f"{chemin}PFE_Nemoh/MyTestCases/{BEM_file}/results/CM.dat"
            CM_PIT_N3_REF = f"{chemin_PIT_N3}Global_Madd_{PIT_type}{beta_value}d{distance}.00.dat"
            CA_BEM = f"{chemin}PFE_Nemoh/MyTestCases/{BEM_file}/results/CA.dat"
            CA_PIT_N3_REF = f"{chemin_PIT_N3}Global_Crad_{PIT_type}{beta_value}d{distance}.00.dat"

            Fex_BEM = f"{chemin}PFE_Nemoh/MyTestCases/{BEM_file}/results/ExcitationForce.tec"  
            Fex_PIT_N3_REF = f"{chemin_PIT_N3}Global_Fe_abs_{PIT_type}{beta_value}d{distance}.00.dat"
            Fex_phase_PIT_N3_REF = f"{chemin_PIT_N3}Global_Fe_phase_{PIT_type}{beta_value}d{distance}.00.dat"
            
            for i, E in enumerate(Ne) : 
                print("\n  Ne=", E)
                CM_PIT_N3 = f"{chemin_PIT_N3}Global_E{E}_Madd_{PIT_type}{beta_value}d{distance}.00.dat"
                CA_PIT_N3 = f"{chemin_PIT_N3}Global_E{E}_Crad_{PIT_type}{beta_value}d{distance}.00.dat"
                Fex_PIT_N3 = f"{chemin_PIT_N3}Global_E{E}_Fe_abs_N{PIT_type}{beta_value}d{distance}.00.dat"
                Fex_phase_PIT_N3 = f"{chemin_PIT_N3}Global_E{E}_Fe_phase_{PIT_type}{beta_value}d{distance}.00.dat"
                
                testE=f"E{E}_"
                
                L_inf_error, L2_error = CompareResultsRAD(CM_BEM, CM_PIT_N3, "CM", testE, distance, f"Nb{Nb}_{layout}", Nb, 0.05, 10) # tol sur erreur, precision sur valeur
                L_inf_error, L2_error = CompareResultsRAD(CA_BEM, CA_PIT_N3, "CA", testE, distance, f"Nb{Nb}_{layout}", Nb, 0.05, 10)
                L_inf_error_abs, L2_error_abs, L_inf_error_ph, L2_error_ph = CompareResultsFex(Fex_BEM, Fex_PIT_N3, Fex_phase_PIT_N3, testE, distance, f"Nb{Nb}_{layout}", Nb, 0.01, 10)
            print("\n  Ne=0")
            L_inf_error, L2_error = CompareResultsRAD(CM_BEM, CM_PIT_N3_REF, "CM", "", distance, f"Nb{Nb}_{layout}", Nb, 0.05, 10) # tol sur erreur, precision sur valeur
            L_inf_error, L2_error = CompareResultsRAD(CA_BEM, CA_PIT_N3_REF, "CA", "", distance, f"Nb{Nb}_{layout}", Nb, 0.05, 10)
            L_inf_error_abs, L2_error_abs, L_inf_error_ph, L2_error_ph = CompareResultsFex(Fex_BEM, Fex_PIT_N3_REF, Fex_phase_PIT_N3_REF, "", distance, f"Nb{Nb}_{layout}", Nb, 0.01, 10)
            
            # Error_RAD_E(CM_BEM, CM_PIT_N3_REF, PIT_files_CM, "Added_Mass", indice, distance, titre, Ne)
            # Error_RAD_E(CA_BEM, CA_PIT_N3_REF, PIT_files_CA, "Damping", indice, distance, titre, Ne)
            # Error_Fex_E(Fex_BEM, Fex_PIT_N3_REF, Fex_phase_PIT_N3_REF, PIT_files_Fex, PIT_files_Fex_ph, DOF, distance, titre, Ne)
    
