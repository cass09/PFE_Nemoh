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

def CompareResultsRAD(fileCAP, filePIT, nom, comment, param_d, type, Nb, tol, precision) :
    freq_CAP, Coef_CAP = read_RAD_PIT(fileCAP)
    freq_PIT, Coef_PIT = read_RAD_PIT(filePIT)
    if np.all(np.isclose(freq_CAP, freq_PIT, atol=1e-4)):
        if np.shape(Coef_CAP) != np.shape(Coef_PIT):
            print("CAP", np.shape(Coef_CAP), "PIT", np.shape(Coef_PIT))
            raise ValueError("tailles différentes")
        else :
            errorL2 = np.sqrt(np.sum((np.array(Coef_CAP) - np.array(Coef_PIT))**2))
            errorINF = np.max(np.abs(np.array(Coef_CAP) - np.array(Coef_PIT)))

            normL2_CAP = np.sqrt(np.sum(np.array(Coef_CAP)**2))
            relative_errorL2 = errorL2 / normL2_CAP if normL2_CAP != 0 else np.nan

            max_CAP = np.max(np.abs(np.array(Coef_CAP)))
            relative_errorINF = errorINF / max_CAP if max_CAP != 0 else np.nan

            All_errors=True
            if All_errors:
                dof=6*Nb
                n_error=0
                for i in range(dof):
                    for j in range(dof):
                        CAP = []
                        PIT = []                       

                        CAP = [0 if np.abs(sublist[i][j]) < precision else sublist[i][j] for sublist in Coef_CAP]
                        PIT = [0 if np.abs(sublist[i][j]) < precision else sublist[i][j] for sublist in Coef_PIT]
                        # Remplacer les valeurs < precision par 0
                        if np.max(np.abs(np.array(CAP))) != 0 and np.max(np.abs(np.array(PIT)))!=0:
                            error=np.max(np.abs(np.array(CAP)-np.array(PIT)))/np.max(np.abs(np.array(CAP)))
                        else :
                            error=0
                        with open(f"EcartsRelatifs_{nom}_{type}.dat", "a") as f:
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
        print('CAP', freq_CAP)
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

def CompareResultsFex(fileBEM, fileBEM_phase, filePIT, filePIT_phase, nom, comment, param_d, type, Nb, tol, precision) :
    beta_BEM, freq_BEM, Coef_BEM_abs = read_Fex_PIT(fileBEM)
    beta_BEM, freq_BEM, Coef_BEM_ph = read_Fex_PIT(fileBEM_phase)
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
                    
                    with open(f"EcartsRelatifs_{nom}_Fex_abs_{type}.dat", "a") as f:
                        if i==0:
                            f.write("\n")
                            f.write(f"{param_d}     ")
                        f.write(f"{err_abs:16.8f}     ")
                    if err_abs<=tol:
                        n_err_abs=n_err_abs+1
                        
                    with open(f"EcartsRelatifs_{nom}_Fex_ph_{type}.dat", "a") as f:
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


def Trace_RAD(fileBEM, fileCAP, filePIT, filePIT_REF, file, nom, indice, titre, d):
    freq_BEM, Coef_BEM = read_RAD_PIT(fileBEM)
    freq_CAP, Coef_CAP = read_RAD_PIT(fileCAP)
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
            norm=norm*np.array(freq_CAP)
            unit="/(ρπa²hw)"
        else :
            unit="/(ρπa²h)"
        bis="AD_"
        axeX="ka (a=5m)"
        if np.isscalar(norm):
            norm_array = np.full_like(freq_CAP, norm, dtype=float)
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
    plt.xlabel(f"{axeX}")
    plt.ylabel(f"{nom} {unit}")
    plt.grid(True)

    for k, ij in enumerate(indice): 
        i=ij[0]-1
        j=ij[1]-1
        plt.title(titre + f", Coefficent {ij[0]}_{ij[1]}")
        color = couleurs[k % len(couleurs)]
        marker = markers[k % len(markers)]
        color2 = couleurs[k+len(indice) % len(couleurs)]
        marker2 = markers[k+len(indice) % len(markers)]
        plt.plot(freq_BEM, [mat[i][j] for mat in Coef_BEM], linestyle='-', marker=marker2, color='black', markersize=3, label=f"BEM with Capytaine")
        if REF :
            color3 = couleurs[k+len(indice)+1 % len(couleurs)]
            marker3 = markers[k+len(indice)+2 % len(markers)]
            plt.plot(freq_PIT_REF, [mat[i][j] for mat in Coef_PIT_REF], linestyle='--', marker=marker3, color='b', markersize=6, label=f"IT Outer Cylinder Method")
    
        if AD :
            plt.plot(a*data[:,1], [mat[i][j]/ norm_array[c] for c, mat in enumerate(Coef_PIT)], linestyle=':', marker=marker2, color=color2, markersize=5, label=f"PIT {ij[0]}_{ij[1]}")
            plt.plot(a*data[:,1], [mat[i][j]/ norm_array[c] for c, mat in enumerate(Coef_CAP)], linestyle=':', marker=marker, color=color, markersize=5, label=f"CAP {ij[0]}_{ij[1]}")
        else :
            plt.plot(freq_PIT, [mat[i][j] for mat in Coef_PIT], linestyle=':', marker=marker2, color='r', markersize=5, label=f"IT Source Terms Method with Nemoh")
            plt.plot(freq_CAP, [mat[i][j] for mat in Coef_CAP], linestyle=':', marker=marker, color='m', markersize=3, label=f"IT Source Terms Method with Capytaine")
        plt.legend(loc='best', fontsize='small', frameon=True)
    plt.savefig(f"{file}_{bis}{nom}_{d}_M{indice}.png", dpi=300)
    plt.close()
    return

def Trace_Fex(fileBEM, fileBEM_phase, fileCAP, fileCAP_phase, filePIT, filePIT_phase, filePIT_REF, filePIT_phase_REF, file, num_dof, beta, titre, d):
    beta_BEM, freq_BEM, Coef_BEM_abs = read_Fex_PIT(fileBEM)
    beta_BEM, freq_BEM, Coef_BEM_ph = read_Fex_PIT(fileBEM_phase)
    beta_CAP, freq_CAP, Coef_CAP_abs = read_Fex_PIT(fileCAP)
    beta_CAP, freq_CAP, Coef_CAP_ph = read_Fex_PIT(fileCAP_phase)
    beta_PIT, freq_PIT, Coef_PIT = read_Fex_PIT(filePIT)
    beta_PIT, freq_PIT, Phase_PIT = read_Fex_PIT(filePIT_phase)
    REF=False
    if filePIT_REF != filePIT :
        beta_PIT_REF, freq_PIT_REF, Coef_PIT_REF = read_Fex_PIT(filePIT_REF)
        beta_PIT_REF, freq_PIT_REF, Phase_PIT_REF = read_Fex_PIT(filePIT_phase_REF)
        REF=True
    couleurs = ['r', 'b', 'm', 'gold', 'lime', 'g', 'c', 'm', 'y', 'k']  # palette de couleurs (réutilisée si plus de 7 courbes)
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
            plt.plot(freq_CAP, Coef_CAP_abs[k][:, j], linestyle='-', marker=marker, color=color, markersize=3, label=f"CAP beta= {beta_CAP[k]}°")   
            plt.plot(freq_PIT, [row[j] for row in Coef_PIT[k]], linestyle=':', marker=marker, color=color, markersize=5, label=f"PIT beta= {beta_PIT[k]:.1f}°")   
        plt.legend()
        plt.grid(True)
        plt.savefig(f"{file}_BETAS{k+1}_Fex_abs_{d}_B{num_dof}.png", dpi=300)
        plt.close()

        plt.figure(figsize=(10, 8))
        plt.title(titre+ f", DOF {num_dof[0]}")
        plt.xlabel("Frequency (rad/s)")
        plt.ylabel("phase(Fex)")
        for k, beta in enumerate(beta_CAP) :
            color = couleurs[k % len(couleurs)]
            marker = markers[k % len(markers)]
            plt.plot(freq_CAP, Coef_CAP_ph[k][:, j], linestyle='-', marker=marker, color=color, markersize=3, label=f"CAP beta= {beta_CAP[k]}°")   
            plt.plot(freq_PIT, [row[j] for row in Phase_PIT[k]], linestyle=':', marker=marker, color=color, markersize=5, label=f"PIT beta= {beta_PIT[k]:.1f}°")   
        plt.legend()
        plt.grid(True)
        plt.savefig(f"{file}_BETAS{k+1}_Fex_phase_{d}_B{num_dof}.png", dpi=300)
        plt.close()

    else : # ONE wave direction !!
        i=0
        plt.figure(figsize=(10, 8))
        plt.xlabel("Frequency (rad/s)")
        plt.ylabel("|Fex| (N)")
        nom_dof=""
        for k, dof in enumerate(num_dof) : 
            j=dof-1
            plt.title(titre + f", beta={beta_CAP[i]:.2f}°" + f", DOF {j+1}")
            nom_dof=nom_dof+"_"+f"{dof}"
            color = couleurs[k % len(couleurs)]
            marker = markers[k % len(markers)]
            color2 = couleurs[k+len(num_dof) % len(couleurs)]
            marker2 = markers[k+len(num_dof) % len(markers)]
            plt.plot(freq_BEM, [row[j] for row in Coef_BEM_abs[i]], linestyle='-', marker=marker2, color='black', markersize=3, label=f"BEM with Capytaine")   
            if REF :
                color3 = couleurs[k+len(num_dof)+1 % len(couleurs)]
                marker3 = markers[k+len(num_dof)+2 % len(markers)]
                plt.plot(freq_PIT_REF, [row[j] for row in Coef_PIT_REF[i]], linestyle='--', marker=marker3, color='b', markersize=6, label=f"IT Outer Cylinder Method")   
            plt.plot(freq_PIT, [row[j] for row in Coef_PIT[i]], linestyle=':', marker=marker2, color='r', markersize=5, label=f"IT Source Terms Method with Nemoh")   
            plt.plot(freq_CAP, [row[j] for row in Coef_CAP_abs[i]], linestyle=':', marker=marker, color='m', markersize=3, label=f"IT Source Terms Method with Capytaine")   
            plt.legend()
        plt.grid(True)
        plt.savefig(f"{file}_Fex_abs_{d}_B{nom_dof}.png", dpi=300)
        plt.close()

        plt.figure(figsize=(10, 8))
        plt.xlabel("Frequency (rad/s)")
        plt.ylabel("phase(Fex) (N)")
        i=0
        nom_dof=""
        for k, dof in enumerate(num_dof) : 
            j=dof-1
            plt.title(titre + f", beta={beta_CAP[i]:.2f}°, DOF {j+1}")
            nom_dof=nom_dof+"_"+f"{dof}"
            color = couleurs[k % len(couleurs)]
            marker = markers[k % len(markers)]
            color2 = couleurs[k+len(num_dof) % len(couleurs)]
            marker2 = markers[k+len(num_dof) % len(markers)]
            plt.plot(freq_BEM, [row[j] for row in Coef_BEM_ph[i]], linestyle='-', marker=marker2, color='black', markersize=3, label=f"BEM with Capytaine")   
            if REF :
                color3 = couleurs[k+len(num_dof)+1 % len(couleurs)]
                marker3 = markers[k+len(num_dof)+2 % len(markers)]
                plt.plot(freq_PIT_REF, [row[j] for row in Phase_PIT_REF[i]], linestyle='--', marker=marker3, color='b', markersize=6, label=f"IT Outer Cylinder Method")   
            plt.plot(freq_PIT, [row[j] for row in Phase_PIT[i]], linestyle=':', marker=marker2, color='r', markersize=5, label=f"IT Source Terms Method with Nemoh")   
            plt.plot(freq_CAP, [row[j] for row in Coef_CAP_ph[i]], linestyle=':', marker=marker, color='m', markersize=3, label=f"IT Source Terms Method with Capytaine")   
            plt.legend()
        plt.grid(True)
        plt.savefig(f"{file}_Fex_phase_{d}_B{nom_dof}.png", dpi=300)
        plt.close()
    return

def plot_RAD_E(fileBEM, filePIT_REF, PIT_files, nom, indice, distance, titre, Ne):
    freq_BEM, Coef_BEM = read_RAD_PIT(fileBEM)
    freq_PIT_REF, Coef_PIT_REF = read_RAD_PIT(filePIT_REF)
    couleurs = ['g', 'gold', 'c', 'lime', 'b', 'm','y', 'k']  # palette de couleurs (réutilisée si plus de 7 courbes)
    markers = ['s','<','x', '+', '*']
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
        plt.plot(freq_BEM, [mat[i][j] for mat in Coef_BEM], linestyle='-', marker='o', color='black', markersize=3, label=f"BEM CAP")
        plt.plot(freq_PIT_REF, [mat[i][j] for mat in Coef_PIT_REF], linestyle='--', marker='^', color='r', markersize=4, label=f"IT CAP L=0")
        for idx, file in enumerate(PIT_files) :
            color = couleurs[idx % len(couleurs)]
            marker = markers[idx % len(markers)]
            freq_PIT, Coef_PIT = read_RAD_PIT(file)
            plt.plot(freq_PIT, [mat[i][j] for mat in Coef_PIT], linestyle=':', marker=marker, color=color, markersize=6, label=f"IT CAP - L={Ne[idx]}")
        plt.legend(loc='best', fontsize='small', frameon=True)
        plt.savefig(f"CAP_E_{nom}_d{distance}_M{ij}.png", dpi=300)
        plt.close()
    return

def plot_Fex_E(fileBEM, fileBEM_phase, filePIT_REF, filePIT_phase_REF, filesPIT, filesPIT_phase, num_dof, distance, titre, Ne):
    beta_BEM, freq_BEM, Coef_BEM_abs = read_Fex_PIT(fileBEM)
    beta_BEM, freq_BEM, Coef_BEM_ph = read_Fex_PIT(fileBEM_phase)
    beta_PIT_REF, freq_PIT_REF, Coef_PIT_REF = read_Fex_PIT(filePIT_REF)
    beta_PIT_REF, freq_PIT_REF, Phase_PIT_REF = read_Fex_PIT(filePIT_phase_REF)
    couleurs = ['g', 'gold', 'c', 'b', 'm', 'lime', 'y', 'k']  # palette de couleurs (réutilisée si plus de 7 courbes)
    markers = ['s', '<', 'x', '+', '*']
    i=0
    
    for k, dof in enumerate(num_dof) : 
        plt.figure(figsize=(10, 8))
        plt.title(titre + f", beta={beta_BEM[i]:.2f}°, DOF {dof}")
        plt.xlabel("Frequency (rad/s)")
        plt.ylabel("|Fex| (N)")
        j=dof-1
        plt.plot(freq_BEM, [row[j] for row in Coef_BEM_abs[i]], linestyle='-', marker='o', color='black', markersize=3, label=f"BEM CAP")   
        plt.plot(freq_PIT_REF, [row[j] for row in Coef_PIT_REF[i]], linestyle='--', marker='^', color='r', markersize=4, label=f"IT CAP - L=0")   
        for idx, file in enumerate(filesPIT) :
            marker = markers[idx % len(markers)]
            color = couleurs[idx % len(couleurs)]
            beta_PIT, freq_PIT, Coef_PIT = read_Fex_PIT(file)
            plt.plot(freq_PIT, [row[j] for row in Coef_PIT[i]], linestyle=':', marker=marker, color=color, markersize=6, label=f"IT CAP - L={Ne[idx]}")   
        plt.legend()
        plt.grid(True)
        plt.savefig(f"CAP_E_Fex_abs_d{distance}_B{dof}.png", dpi=300)
        plt.close()

    
    i=0
    for k, dof in enumerate(num_dof) : 
        plt.figure(figsize=(10, 8))
        plt.title(titre + f", beta={beta_BEM[i]:.2f}°, DOF {dof}")
        plt.xlabel("Frequency (rad/s)")
        plt.ylabel("phase(Fex) (N)")
        j=dof-1
        plt.plot(freq_BEM, [row[j] for row in Coef_BEM_ph[i]], linestyle='-', marker='o', color='black', markersize=3, label=f"BEM CAP")   
        plt.plot(freq_PIT_REF, [row[j] for row in Phase_PIT_REF[i]], linestyle='--', marker='^', color='r', markersize=4, label=f"IT CAP -L=0")   
        for idx, file in enumerate(filesPIT_phase) :
            marker = markers[idx % len(markers)]
            color = couleurs[idx % len(couleurs)]
            beta_PIT, freq_PIT, Phase_PIT = read_Fex_PIT(file)
            plt.plot(freq_PIT, [row[j] for row in Phase_PIT[i]], linestyle=':', marker=marker, color=color, markersize=6, label=f"IT CAP - L={Ne[idx]}")   
        plt.legend()
        plt.grid(True)
        plt.savefig(f"CAP_E_Fex_phase_d{distance}_B{dof}.png", dpi=300)
        plt.close()
    return

def CompareResultsRAD_E(fileBEM, filePIT, nom, param_d, type, Nb, Ne, tol, precision) :
    freq_BEM, Coef_BEM = read_RAD_PIT(fileBEM)
    freq_PIT, Coef_PIT = read_RAD_PIT(filePIT)
    if np.all(np.isclose(freq_BEM, freq_PIT, atol=1e-4)):
        if np.shape(Coef_BEM) != np.shape(Coef_PIT):
            print("BEM", np.shape(Coef_BEM), "PIT", np.shape(Coef_PIT))
            raise ValueError("tailles différentes")
        else :
            dof=6*Nb
            n_error=0
            for i in range(dof):
                for j in range(dof):
                    BEM = []
                    PIT = []                       
                    BEM = [0 if np.abs(sublist[i][j]) < precision else sublist[i][j] for sublist in Coef_BEM]
                    PIT = [0 if np.abs(sublist[i][j]) < precision else sublist[i][j] for sublist in Coef_PIT]
                    if np.max(np.abs(np.array(BEM))) != 0 and np.max(np.abs(np.array(PIT)))!=0:
                        error=np.max(np.abs(np.array(BEM)-np.array(PIT)))/np.max(np.abs(np.array(BEM)))
                    else :
                        error=0
                    with open(f"EcartsRelatifs_Ne{Ne}_{nom}_{type}.dat", "a") as f:
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
    
    return "OK"
def CompareResultsFex_E(fileBEM, fileBEM_phase, filePIT, filePIT_phase, param_d, type, Nb, Ne, tol, precision) :
    beta_BEM, freq_BEM, Coef_BEM_abs = read_Fex_PIT(fileBEM)
    beta_BEM, freq_BEM, Coef_BEM_ph = read_Fex_PIT(fileBEM_phase)
    beta_PIT, freq_PIT, Coef_PIT = read_Fex_PIT(filePIT)
    beta_PIT, freq_PIT, Phase_PIT = read_Fex_PIT(filePIT_phase)
   
    if len(freq_BEM)==len(freq_PIT) and len(beta_BEM)==len(beta_PIT):
        if np.shape(Coef_BEM_abs) != np.shape(Coef_PIT):
            print(np.shape(Coef_BEM_abs), np.shape(Coef_PIT))
            raise ValueError("tailles différentes")
        if np.shape(Coef_BEM_ph) != np.shape(Phase_PIT):
            print(np.shape(Coef_BEM_ph), np.shape(Phase_PIT))
            raise ValueError("tailles différentes")
        else :
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
                    
                with open(f"EcartsRelatifs_Ne{Ne}_Fex_abs_{type}.dat", "a") as f:
                    if i==0:
                        f.write("\n")
                        f.write(f"{param_d}     ")
                    f.write(f"{err_abs:16.8f}     ")
                if err_abs<=tol:
                    n_err_abs=n_err_abs+1
                    
                with open(f"EcartsRelatifs_Ne{Ne}_Fex_ph_{type}.dat", "a") as f:
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
    
    return "OK"

#################################################################
#################################################################
######################### PARAMETERS ############################
#################################################################
#################################################################

param_d=2

N3 = False
COMP_E = True
GRAPHS_E=True
ERROR_E=False

GRAPHS_beta=False

if N3 : 
    param_d=1
    limite=600
    GRAPHS_N3_F=False
    GRAPHS_N3_C=False
elif GRAPHS_beta or COMP_E : 
    limite=param_d
    GRAPHS_N3_F = False
    GRAPHS_N3_C=False
else : 
    limite=param_d
    GRAPHS_N3_F = True
    GRAPHS_N3_C = False


mesh="CylR3"
LID="_h12"
config_a=True
Nb=2
layout="X"
Nw=20
Ndir=1
test="S_"
# beta_value="beta_180.0_"
beta_value=""
PIT_type=f"Nb{Nb}_{layout}_"
PIT_N3_file=f"PIT3_{mesh}{LID}_source"
CAP_file=f"CAP_{mesh}{LID}"
PIT_N3_file_REF=f"PIT3_{mesh}{LID}"
chemin="/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/"
chemin_PIT_N3 = f"{chemin}PIT3_E/wec_inputs/{PIT_N3_file}/resultsIT/"
chemin_PIT_N3_REF = f"{chemin}PIT3_E/wec_inputs/{PIT_N3_file_REF}/resultsIT/"
chemin_CAP=f"{chemin}capytaine/MyTestCases/{CAP_file}/resultsIT/"
chemin_CAP_BEM=f"{chemin}capytaine/MyTestCases/C_BEM_{mesh}{LID}/resultsBEM/Nb{Nb}_{layout}"
titre=f"{mesh} - Nb={Nb} - {layout}, Nw={Nw}, Nbeta=13, Ndof=6, Ndir={Ndir}"

if GRAPHS_beta : 
    Ndir=5
    beta_value=f"multi_beta_{Ndir}_"
    PIT_type=""

while param_d<=limite : 
    if config_a : 
        distance=f"da{param_d}"  
        titre=f"{titre}, d={param_d*3}m"
    else :
        distance=f"d{param_d}"  
        titre=f"{titre}, d={param_d}m"
    CM_BEM = f"{chemin_CAP_BEM}_d{param_d}/Capytaine_Madd.dat"
    CM_CAP = f"{chemin_CAP}CapytaineIT_{test}Madd_{PIT_type}{beta_value}{distance}.00.dat"
    CM_PIT_N3 = f"{chemin_PIT_N3}Global_{test}Madd_{PIT_type}{beta_value}{distance}.00.dat"
    CM_PIT_N3_REF = f"{chemin_PIT_N3_REF}Global_Madd_{PIT_type}{beta_value}{distance}.00.dat"
    
    CA_BEM = f"{chemin_CAP_BEM}_d{param_d}/Capytaine_Crad.dat"
    CA_CAP = f"{chemin_CAP}CapytaineIT_{test}Crad_{PIT_type}{beta_value}{distance}.00.dat"
    CA_PIT_N3 = f"{chemin_PIT_N3}Global_{test}Crad_{PIT_type}{beta_value}{distance}.00.dat"
    CA_PIT_N3_REF = f"{chemin_PIT_N3_REF}Global_Crad_{PIT_type}{beta_value}{distance}.00.dat"

    Fex_BEM = f"{chemin_CAP_BEM}_d{param_d}/Capytaine_Fe_abs.dat"  
    Fex_phase_BEM = f"{chemin_CAP_BEM}_d{param_d}/Capytaine_Fe_phase.dat"  
            
    Fex_CAP =f"{chemin_CAP}CapytaineIT_{test}Fe_abs_{PIT_type}{beta_value}{distance}.00.dat"
    Fex_PIT_N3 = f"{chemin_PIT_N3}Global_{test}Fe_abs_{PIT_type}{beta_value}{distance}.00.dat"
    Fex_PIT_N3_REF = f"{chemin_PIT_N3_REF}Global_Fe_abs_{PIT_type}{beta_value}{distance}.00.dat"
    Fex_phase_CAP = f"{chemin_CAP}CapytaineIT_{test}Fe_phase_{PIT_type}{beta_value}{distance}.00.dat"
    Fex_phase_PIT_N3 = f"{chemin_PIT_N3}Global_{test}Fe_phase_{PIT_type}{beta_value}{distance}.00.dat"
    Fex_phase_PIT_N3_REF = f"{chemin_PIT_N3_REF}Global_Fe_phase_{PIT_type}{beta_value}{distance}.00.dat"

    if N3 : 
        print("\n ------------Distance=", param_d, "-----------------")
        list_files_CM=[CM_CAP, CM_PIT_N3, CM_PIT_N3_REF]
        list_files_CA=[CA_CAP, CA_PIT_N3, CA_PIT_N3_REF]
        list_files_F_abs=[Fex_CAP, Fex_PIT_N3, Fex_PIT_N3_REF]
        list_files_F_ph=[Fex_phase_CAP, Fex_phase_PIT_N3, Fex_phase_PIT_N3_REF]
        list_noms=["CAP_S", "Nemoh_S", "Nemoh_OC"]
        for ind in range(3):
            print(list_noms[ind])
            with open(f"EcartsRelatifs_{list_noms[ind]}_Nb{Nb}_{layout}.dat", "a") as f:
                f.write(f"{param_d}     ")
                L_inf_error, L2_error = CompareResultsRAD(CM_BEM, list_files_CM[ind], f"{list_noms[ind]}_CM", test, param_d, f"Nb{Nb}_{layout}", Nb, 0.02, 0) # tol sur erreur, precision sur valeur
                f.write(f"{L_inf_error:16.8f}     {L2_error:16.8f}      ")
                L_inf_error, L2_error = CompareResultsRAD(CA_BEM, list_files_CA[ind], f"{list_noms[ind]}_CA", test, param_d, f"Nb{Nb}_{layout}", Nb, 0.02, 0)
                f.write(f"{L_inf_error:16.8f}     {L2_error:16.8f}      ")
                L_inf_error_abs, L2_error_abs, L_inf_error_ph, L2_error_ph = CompareResultsFex(Fex_BEM, Fex_phase_BEM, list_files_F_abs[ind], list_files_F_ph[ind], list_noms[ind], test, param_d, f"Nb{Nb}_{layout}", Nb, 0.02, 0)
                f.write(f"{L_inf_error_abs:16.8f}     {L2_error_abs:16.8f}      ")
                f.write(f"{L_inf_error_ph:16.8f}      {L2_error_ph:16.8f}\n")

    if GRAPHS_N3_F :
        for ind, i in enumerate([[1], [3], [9]]) :
            print("dof", i)
            Trace_Fex(Fex_BEM, Fex_phase_BEM, Fex_CAP, Fex_phase_CAP, Fex_PIT_N3, Fex_phase_PIT_N3, Fex_PIT_N3_REF, Fex_phase_PIT_N3_REF, f"IT_{test}Nb{Nb}_{layout}", i, False, titre, distance)
    if GRAPHS_N3_C : 
        indices=[[(1,1)], [(3,3)], [(1, 7)], [(3, 9)], [(5,11)]]
        for ij in indices : 
            print("indice", ij)
            Trace_RAD(CA_BEM, CA_CAP, CA_PIT_N3, CA_PIT_N3_REF, f"IT_{test}Nb{Nb}_{layout}", "Damping", ij, titre, distance)
            Trace_RAD(CM_BEM, CM_CAP, CM_PIT_N3, CM_PIT_N3_REF, f"IT_{test}Nb{Nb}_{layout}", "Added_Mass", ij, titre, distance)
    if GRAPHS_beta : 
        for ind, i in enumerate([[3],[9], [5],[11], [1],[7]])  :
            print("dof", i)
            Trace_Fex(Fex_CAP, Fex_PIT_N3, Fex_phase_PIT_N3, Fex_PIT_N3_REF, Fex_phase_PIT_N3_REF, f"IT_{test}Nb{Nb}_{layout}", i, True, titre, distance)
            
    param_d=param_d*2


#################################################################
#################################################################
################### EVANESCENT PROBLEM ########################
#################################################################
#################################################################


if COMP_E :   
    if GRAPHS_E : 
        print("\n ------------ PLOT -----------------")
        vect_d=[4]
    if ERROR_E : 
        print("\n ------------ ERROR -----------------")
        vect_d=[]
        for k in range(0, 10):
            vect_d.append(2**k) 
    Ne=[1,6]
    DOF=[1, 3, 5]
    indice=[(1, 1), (3,9), (1,7), (5,11)]
    print("\n Ne=", Ne)
    print("\n DOF=", DOF )
    print("\n indice=", indice )
    chemin_PIT_N3=f"{chemin}capytaine/MyTestCases/CAP_{mesh}_h12/resultsIT/"
    chemin_CAP_BEM = f"{chemin}capytaine/MyTestCases/C_BEM_{mesh}_h12/resultsBEM/Nb{Nb}_{layout}"

    for c, distance in enumerate(vect_d):
        print("\n  Distance=", distance)
        titre=f"{mesh} - Nb={Nb} {layout}, Nw=20, Nbeta=13, Ndof=6, Ndir=1, d={distance*3}m"

        CM_BEM = f"{chemin_CAP_BEM}_da{distance}/Capytaine_Madd.dat"
        CM_PIT_N3_REF = f"{chemin_PIT_N3}CapytaineIT_S_Madd_{PIT_type}{beta_value}da{distance}.00.dat"
        CA_BEM = f"{chemin_CAP_BEM}_da{distance}/Capytaine_Crad.dat"
        CA_PIT_N3_REF = f"{chemin_PIT_N3}CapytaineIT_S_Crad_{PIT_type}{beta_value}da{distance}.00.dat"

        Fex_BEM = f"{chemin_CAP_BEM}_da{distance}/Capytaine_Fe_abs.dat"  
        Fex_phase_BEM = f"{chemin_CAP_BEM}_da{distance}/Capytaine_Fe_phase.dat"  
        Fex_PIT_N3_REF = f"{chemin_PIT_N3}CapytaineIT_S_Fe_abs_{PIT_type}{beta_value}da{distance}.00.dat"
        Fex_phase_PIT_N3_REF = f"{chemin_PIT_N3}CapytaineIT_S_Fe_phase_{PIT_type}{beta_value}da{distance}.00.dat"

        PIT_files_CM=[]
        PIT_files_CA=[]
        PIT_files_Fex=[]
        PIT_files_Fex_ph=[]
            
        for i, E in enumerate(Ne) : 
            CM_PIT_N3 = f"{chemin_PIT_N3}CapytaineIT_S_E{E}_Madd_{PIT_type}{beta_value}da{distance}.00.dat"
            CA_PIT_N3 = f"{chemin_PIT_N3}CapytaineIT_S_E{E}_Crad_{PIT_type}{beta_value}da{distance}.00.dat"
            Fex_PIT_N3 = f"{chemin_PIT_N3}CapytaineIT_S_E{E}_Fe_abs_{PIT_type}{beta_value}da{distance}.00.dat"
            Fex_phase_PIT_N3 = f"{chemin_PIT_N3}CapytaineIT_S_E{E}_Fe_phase_{PIT_type}{beta_value}da{distance}.00.dat"
            if ERROR_E : 
                CompareResultsRAD_E(CM_BEM, CM_PIT_N3, "CM", distance, f"Nb{Nb}_{layout}", Nb, E, 0.02, 0) # tol sur erreur, precision sur valeur
                CompareResultsRAD_E(CA_BEM, CA_PIT_N3, "CA", distance, f"Nb{Nb}_{layout}", Nb, E, 0.02, 0)
                CompareResultsFex_E(Fex_BEM, Fex_phase_BEM, Fex_PIT_N3, Fex_phase_PIT_N3, distance, f"Nb{Nb}_{layout}", Nb, E, 0.01, 0)
          
            PIT_files_CM.append(CM_PIT_N3) 
            PIT_files_CA.append(CA_PIT_N3) 
            PIT_files_Fex.append(Fex_PIT_N3) 
            PIT_files_Fex_ph.append(Fex_phase_PIT_N3) 
        if GRAPHS_E : 
            plot_RAD_E(CM_BEM, CM_PIT_N3_REF, PIT_files_CM, "Added_Mass", indice, distance, titre, Ne)
            plot_RAD_E(CA_BEM, CA_PIT_N3_REF, PIT_files_CA, "Damping", indice, distance, titre, Ne)
            plot_Fex_E(Fex_BEM, Fex_phase_BEM, Fex_PIT_N3_REF, Fex_phase_PIT_N3_REF, PIT_files_Fex, PIT_files_Fex_ph, DOF, distance, titre, Ne)
    
