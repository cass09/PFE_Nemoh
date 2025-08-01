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


def Trace_RAD(fileBEM, filePIT, filePIT_REF, file, nom, indice, titre, distances, i_w):
    Coef_BEM=[]
    Coef_PIT=[]
    Coef_PIT_REF=[]
    for ind in range(len(fileBEM)) : 
        freq_BEM, Coef_BEM_w = read_RAD_PIT(fileBEM[ind])
        freq_PIT, Coef_PIT_w = read_RAD_PIT(filePIT[ind])
        Coef_BEM.append(Coef_BEM_w[i_w])
        Coef_PIT.append(Coef_PIT_w[i_w])
        REF=False
        if filePIT_REF != filePIT :
            freq_PIT_REF, Coef_PIT_REF_w = read_RAD_PIT(filePIT_REF[ind])
            REF=True
            Coef_PIT_REF.append(Coef_PIT_REF_w[i_w])
    freq=freq_BEM[i_w]
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
        axeX="Distance (m)"
        bis=""
        if nom=="Damping":
            unit="(kg/s)"
        else :
            unit="(kg)"

    plt.figure(figsize=(8, 6))
    plt.title(titre +f", w={freq} rad/s")
    plt.xlabel(f"{axeX}")
    plt.ylabel(f"{nom} {unit}")
    plt.grid(True)
    plt.xscale("log", base=2) 

    for k, ij in enumerate(indice): 
        i=ij[0]-1
        j=ij[1]-1
        color = couleurs[k % len(couleurs)]
        marker = markers[k % len(markers)]
        color2 = couleurs[k+len(indice) % len(couleurs)]
        marker2 = markers[k+len(indice) % len(markers)]
        if AD :
            plt.plot(a*data[:,1], [mat[i][j]/ norm_array[c] for c, mat in enumerate(Coef_BEM)], linestyle='-', marker=marker, color=color, markersize=5, label=f"BEM {ij[0]}_{ij[1]}")
            plt.plot(a*data[:,1], [mat[i][j]/ norm_array[c] for c, mat in enumerate(Coef_PIT)], linestyle=':', marker=marker2, color=color2, markersize=5, label=f"PIT {ij[0]}_{ij[1]}")
        else :
            plt.plot(distances, [mat[i][j] for mat in Coef_BEM], linestyle='-', marker=marker, color=color, markersize=5, label=f"CAP BEM {ij[0]}_{ij[1]}")
            plt.plot(distances, [mat[i][j] for mat in Coef_PIT], linestyle=':', marker=marker2, color=color2, markersize=5, label=f"CAP IT {ij[0]}_{ij[1]}")
        
        if REF :
            marker3 = markers[k+len(indice)+2 % len(markers)]
            plt.plot(distances, [mat[i][j] for mat in Coef_PIT_REF], linestyle=':', marker=marker3, color='black', markersize=6, label=f"Nemoh IT Outer Cylinder {ij[0]}_{ij[1]}")
    plt.legend(loc='best', fontsize='small', frameon=True)
    plt.savefig(f"{file}_{bis}{nom}_w{i_w+1}_M{indice}.png", dpi=300)
    plt.close()
    return

def Trace_Fex(fileBEM, fileBEM_phase, filePIT, filePIT_phase, filePIT_REF, filePIT_phase_REF, file, num_dof, beta, titre, distances, i_w):
    Coef_BEM_abs=[]
    Coef_BEM_ph=[]
    Coef_PIT=[]
    Coef_PIT_REF=[]
    Phase_PIT=[]
    Phase_PIT_REF=[]
    for ind in range(len(fileBEM)) :  
        beta_BEM, freq_BEM, Coef_BEM_abs_w =read_Fex_PIT(fileBEM[ind])
        beta_BEM, freq_BEM, Coef_BEM_ph_w =read_Fex_PIT(fileBEM_phase[ind])
        beta_PIT, freq_PIT, Coef_PIT_w = read_Fex_PIT(filePIT[ind])
        beta_PIT, freq_PIT, Phase_PIT_w = read_Fex_PIT(filePIT_phase[ind])
        REF=False
        Coef_BEM_abs_w = np.array(Coef_BEM_abs_w)
        Coef_BEM_ph_w = np.array(Coef_BEM_ph_w)
        Coef_BEM_abs.append(Coef_BEM_abs_w[:,i_w,:])
        Coef_BEM_ph.append(Coef_BEM_ph_w[:,i_w,:])
        Coef_PIT.append(np.array([param[i_w] for param in Coef_PIT_w]))
        Phase_PIT.append(np.array([param[i_w] for param in Phase_PIT_w]))
        if filePIT_REF != filePIT :
            beta_PIT_REF, freq_PIT_REF, Coef_PIT_REF_w = read_Fex_PIT(filePIT_REF[ind])
            beta_PIT_REF, freq_PIT_REF, Phase_PIT_REF_w = read_Fex_PIT(filePIT_phase_REF[ind])
            REF=True
            Coef_PIT_REF.append(np.array([param[i_w] for param in Coef_PIT_REF_w]))
            Phase_PIT_REF.append(np.array([param[i_w] for param in Phase_PIT_REF_w]))
    freq=freq_BEM[i_w]
    Coef_BEM_abs = np.array(Coef_BEM_abs)
    Coef_BEM_ph = np.array(Coef_BEM_ph)
    Coef_PIT = np.array(Coef_PIT)
    Phase_PIT = np.array(Phase_PIT)
    Coef_PIT_REF = np.array(Coef_PIT_REF)
    Phase_PIT_REF = np.array(Phase_PIT_REF)
    couleurs = ['b', 'm', 'r', 'gold', 'lime', 'g', 'c', 'm', 'y', 'k']  # palette de couleurs (réutilisée si plus de 7 courbes)
    markers = ['s', 'o', '<', '^','x', '+', '*']
    if beta : # Several wave directions
        plt.figure(figsize=(10, 8))
        plt.title(titre + f", DOF {num_dof[0]}, w={freq} rad/s")
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
        plt.savefig(f"{file}_BETAS{k+1}_Fex_abs_w{i_w}_B{num_dof}.png", dpi=300)
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
        plt.savefig(f"{file}_BETAS{k+1}_Fex_phase_w{i_w}_B{num_dof}.png", dpi=300)
        plt.close()

    else : # ONE wave direction !!
        i=0
        plt.figure(figsize=(10, 8))
        plt.title(titre + f", beta={beta_BEM[i]:.2f}°, w={freq} rad/s")
        plt.xlabel("Frequency (rad/s)")
        plt.ylabel("|Fex| (N)")
        plt.xscale("log", base=2) 
        nom_dof=""
        for k, dof in enumerate(num_dof) : 
            j=dof-1
            nom_dof=nom_dof+"_"+f"{dof}"
            color = couleurs[k % len(couleurs)]
            marker = markers[k % len(markers)]
            color2 = couleurs[k+len(num_dof) % len(couleurs)]
            marker2 = markers[k+len(num_dof) % len(markers)]
            plt.plot(distances, Coef_BEM_abs[:,i, j], linestyle='-', marker=marker, color=color, markersize=5, label=f"CAP BEM dof {j+1}")   
            plt.plot(distances, Coef_PIT[:,i, j], linestyle=':', marker=marker2, color=color2, markersize=4, label=f" CAP IT dof {j+1}")   
            if REF :
                color3 = couleurs[k+len(num_dof)+1 % len(couleurs)]
                marker3 = markers[k+len(num_dof)+2 % len(markers)]
                plt.plot(distances, Coef_PIT_REF[:,i, j], linestyle=':', marker=marker3, color='black', markersize=6, label=f"Nemoh PIT Outer Cylinder dof {j+1}")   
        plt.legend()
        plt.grid(True)
        plt.savefig(f"{file}_Fex_abs_w{i_w+1}_B{nom_dof}.png", dpi=300)
        plt.close()

        plt.figure(figsize=(10, 8))
        plt.title(titre + f", beta={beta_BEM[i]:.2f}°")
        plt.xlabel("Frequency (rad/s)")
        plt.ylabel("phase(Fex) (N)")
        plt.xscale("log", base=2) 
        i=0
        nom_dof=""
        for k, dof in enumerate(num_dof) : 
            j=dof-1
            nom_dof=nom_dof+"_"+f"{dof}"
            color = couleurs[k % len(couleurs)]
            marker = markers[k % len(markers)]
            color2 = couleurs[k+len(num_dof) % len(couleurs)]
            marker2 = markers[k+len(num_dof) % len(markers)]
            plt.plot(distances, Coef_BEM_ph[:,i, j], linestyle='-', marker=marker, color=color, markersize=5, label=f"CAP BEM dof {j+1}")   
            plt.plot(distances, Phase_PIT[:,i, j], linestyle=':', marker=marker2, color=color2, markersize=4, label=f"CAP IT dof {j+1}")   
            if REF :
                color3 = couleurs[k+len(num_dof)+1 % len(couleurs)]
                marker3 = markers[k+len(num_dof)+2 % len(markers)]
                plt.plot(distances, Phase_PIT_REF[:,i, j], linestyle=':', marker=marker3, color='black', markersize=6, label=f"Nemoh PIT Outer Cylinder dof {j+1}")   
        plt.legend()
        plt.grid(True)
        plt.savefig(f"{file}_Fex_phase_w{i_w+1}_B{nom_dof}.png", dpi=300)
        plt.close()
    return
def plot_RAD_E(fileBEM, filePIT, filePIT_E, file, nom, indice, titre, distances, Ne, i_w):
    Coef_BEM=[]
    Coef_PIT=[]
    Coef_PIT_REF=[[] for _ in range(len(filePIT_E))]
    for ind in range(len(fileBEM)) : 
        freq_BEM, Coef_BEM_w = read_RAD_PIT(fileBEM[ind])
        freq_PIT, Coef_PIT_w = read_RAD_PIT(filePIT[ind])
        Coef_BEM.append(Coef_BEM_w[i_w])
        Coef_PIT.append(Coef_PIT_w[i_w])
        for e in range(len(filePIT_E)) :
            freq_PIT_REF, Coef_PIT_REF_w = read_RAD_PIT(filePIT_E[e][ind])
            Coef_PIT_REF[e].append(Coef_PIT_REF_w[i_w])
    freq=freq_BEM[i_w]
    couleurs = ['g', 'r', 'gold', 'lime', 'c', 'm', 'y', 'k']  # palette de couleurs (réutilisée si plus de 7 courbes)
    markers = ['<', '^','x', '+', '*']
  
    if nom=="Damping":
        unit="(kg/s)"
    else :
        unit="(kg)"

    for k, ij in enumerate(indice):
        plt.figure(figsize=(8, 6))
        plt.xlabel("d/a")
        plt.ylabel(f"{nom} {unit}")
        plt.grid(True) 
        i=ij[0]-1
        j=ij[1]-1
        # plt.xscale("log", base=2) 
        plt.title(titre +f", w={freq} rad/s, Coefficient {ij[0]}_{ij[1]}")
        plt.plot(distances, [mat[i][j] for mat in Coef_BEM], linestyle='--', marker='s', color='b', markersize=5, label=f"CAP BEM")
        plt.plot(distances, [mat[i][j] for mat in Coef_PIT], linestyle='--', marker='o', color='black', markersize=5, label=f"CAP IT L=0")
        for e in range(len(filePIT_E)) :
            marker = markers[e % len(markers)]
            color = couleurs[e % len(couleurs)]
            plt.plot(distances, [mat[i][j] for mat in Coef_PIT_REF[e]], linestyle=':', marker=marker, color=color, markersize=4, label=f"CAP IT L={Ne[e]}")
        plt.legend(loc='best', fontsize='small', frameon=True)
        plt.savefig(f"{file}_{nom}_w{i_w+1}_M{ij}.png", dpi=300)
        plt.close()
    return
def plot_Fex_E(fileBEM, filePIT, filePIT_E, file, num_dof, titre, distances, Ne, i_w):
    Coef_BEM=[]
    Coef_PIT=[]
    Coef_PIT_E=[[] for _ in range(len(filePIT_E))]
    for ind in range(len(fileBEM)) :  
        beta_BEM, freq_BEM, Coef_BEM_abs_w =read_Fex_PIT(fileBEM[ind])
        beta_PIT, freq_PIT, Coef_PIT_w = read_Fex_PIT(filePIT[ind])
        Coef_BEM_abs_w = np.array(Coef_BEM_abs_w)
        Coef_BEM.append(Coef_BEM_abs_w[:,i_w,:])
        Coef_PIT.append(np.array([param[i_w] for param in Coef_PIT_w]))
        for e in range(len(filePIT_E)) :
            beta_PIT_REF, freq_PIT_REF, Coef_PIT_REF_w = read_Fex_PIT(filePIT_E[e][ind])
            Coef_PIT_E[e].append(np.array([param[i_w] for param in Coef_PIT_REF_w]))
    freq=freq_BEM[i_w]
    Coef_BEM = np.array(Coef_BEM)
    Coef_PIT = np.array(Coef_PIT)
    Coef_PIT_E = np.array(Coef_PIT_E)
    couleurs = ['r', 'gold', 'lime', 'g', 'c', 'm', 'y', 'k']  # palette de couleurs (réutilisée si plus de 7 courbes)
    markers = ['o', '<', '^','x', '+', '*']
    
    i=0
    for k, dof in enumerate(num_dof) : 
        j=dof-1
        plt.figure(figsize=(10, 8))
        plt.title(titre + f", beta={beta_BEM[i]:.2f}°, w={freq} rad/s, dof {j+1}")
        plt.xlabel("d/a")
        plt.ylabel("|Fex| (N)")
        # plt.xscale("log", base=2) 
        plt.plot(distances, Coef_BEM[:,i, j], linestyle='--', marker='s', color='b', markersize=5, label=f"CAP BEM")   
        plt.plot(distances, Coef_PIT[:,i, j], linestyle='--', marker='o', color='black', markersize=5, label=f" CAP IT L=0")   
        for e in range(len(filePIT_E)) :
            marker = markers[e % len(markers)]
            color = couleurs[e % len(couleurs)]
            plt.plot(distances, Coef_PIT_E[e][:,i, j], linestyle=':', marker=marker, color=color, markersize=4, label=f"CAP PIT L={Ne[e]}")   
        plt.legend()
        plt.grid(True)
        plt.savefig(f"{file}_w{i_w+1}_B{dof}.png", dpi=300)
        plt.close()    
    return

def plot_RAD_old(fileBEM, filePIT_REF, PIT_files, nom, indice, distance, titre, Ne):
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

def plot_Fex_old(fileBEM, fileBEM_phase, filePIT_REF, filePIT_phase_REF, filesPIT, filesPIT_phase, num_dof, distance, titre, Ne):
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
limite=10
COMP_E = True

GRAPHS_RAD=True
GRAPHS_Fex=True

mesh="CylR3"
config_a=True
Nb=2
layout="X"
Nw=20
Ndir=1
test="S_"
# beta_value="beta_180.0_"
beta_value=""
PIT_type=f"Nb{Nb}_{layout}_"
PIT_N3_file=f"PIT3_{mesh}_source"
PIT_N3_file_REF=f"PIT3_{mesh}"
chemin="/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/"
chemin_PIT_N3 = f"{chemin}PIT3_E/wec_inputs/{PIT_N3_file}/resultsIT/"
chemin_PIT_N3_REF = f"{chemin}PIT3_E/wec_inputs/{PIT_N3_file_REF}/resultsIT/"
# BEM_file=f"BEM_{mesh}_BETAS"
chemin_CAP=f"{chemin}capytaine/MyTestCases/CAP_{mesh}_h12/resultsIT/"
chemin_CAP_BEM=f"{chemin}capytaine/MyTestCases/C_BEM_{mesh}_h12/resultsBEM/Nb{Nb}_{layout}"
titre=f"{mesh} - Nb={Nb} - {layout}, Nw={Nw}, Nbeta=13, Ndof=6, Ndir={Ndir}"


CM_BEM=[]
CM_PIT_N3=[]
CM_PIT_N3_REF=[]
CA_BEM=[]
CA_PIT_N3=[]
CA_PIT_N3_REF=[]
Fex_BEM=[]
Fex_BEM_phase=[]
Fex_PIT_N3=[]
Fex_PIT_N3_REF=[]
Fex_phase_PIT_N3=[]
Fex_phase_PIT_N3_REF=[]
vect_distance=[]
while param_d<=limite : 
    if config_a : 
        distance=f"da{param_d}"  
    else :
        distance=f"d{param_d}"  
    BEM_file=f"{chemin_CAP_BEM}_da{param_d}/"
    CM_BEM_d = f"{BEM_file}Capytaine_Madd.dat"
    CM_PIT_N3_d = f"{chemin_CAP}CapytaineIT_{test}Madd_{PIT_type}{beta_value}{distance}.00.dat"
    CM_PIT_N3_REF_d = f"{chemin_PIT_N3_REF}Global_Madd_{PIT_type}{beta_value}{distance}.00.dat"

    CA_BEM_d = f"{BEM_file}Capytaine_Crad.dat"
    CA_PIT_N3_d = f"{chemin_CAP}CapytaineIT_{test}Crad_{PIT_type}{beta_value}{distance}.00.dat"
    CA_PIT_N3_REF_d = f"{chemin_PIT_N3_REF}Global_Crad_{PIT_type}{beta_value}{distance}.00.dat"

    Fex_BEM_d = f"{BEM_file}Capytaine_Fe_abs.dat"  
    Fex_BEM_phase_d = f"{BEM_file}Capytaine_Fe_phase.dat"  
    Fex_PIT_N3_d = f"{chemin_CAP}CapytaineIT_{test}Fe_abs_{PIT_type}{beta_value}{distance}.00.dat"
    Fex_PIT_N3_REF_d = f"{chemin_PIT_N3_REF}Global_Fe_abs_{PIT_type}{beta_value}{distance}.00.dat"
    Fex_phase_PIT_N3_d = f"{chemin_CAP}CapytaineIT_{test}Fe_phase_{PIT_type}{beta_value}{distance}.00.dat"
    Fex_phase_PIT_N3_REF_d = f"{chemin_PIT_N3_REF}Global_Fe_phase_{PIT_type}{beta_value}{distance}.00.dat"

    CM_BEM.append(CM_BEM_d) 
    CM_PIT_N3.append(CM_PIT_N3_d) 
    CM_PIT_N3_REF.append(CM_PIT_N3_REF_d) 
    CA_BEM.append(CA_BEM_d) 
    CA_PIT_N3.append(CA_PIT_N3_d) 
    CA_PIT_N3_REF.append(CA_PIT_N3_REF_d) 
    Fex_BEM.append(Fex_BEM_d) 
    Fex_BEM_phase.append(Fex_BEM_phase_d) 
    Fex_PIT_N3.append(Fex_PIT_N3_d) 
    Fex_PIT_N3_REF.append(Fex_PIT_N3_REF_d) 
    Fex_phase_PIT_N3.append(Fex_phase_PIT_N3_d) 
    Fex_phase_PIT_N3_REF.append(Fex_phase_PIT_N3_REF_d) 
    vect_distance.append(param_d) 
    param_d=param_d+1

indices_w=[3]
if not COMP_E : 
    for i_w in indices_w : 
        print("i_w : ", i_w)
        if GRAPHS_Fex :
            for ind, i in enumerate([[1], [9]]) :
                print("dof", i)
                Trace_Fex(Fex_BEM, Fex_BEM_phase, Fex_PIT_N3, Fex_phase_PIT_N3, Fex_PIT_N3_REF, Fex_phase_PIT_N3_REF, f"C_D_{test}Nb{Nb}_{layout}", i, False, titre, vect_distance, i_w)
        if GRAPHS_RAD : 
            indices=[[(1, 7)], [(3,9)], [(5,11)]]
            for ij in indices : 
                print("indice", ij)
                Trace_RAD(CA_BEM, CA_PIT_N3, CA_PIT_N3_REF, f"C_D_{test}Nb{Nb}_{layout}", "Damping", ij, titre, vect_distance, i_w)
                Trace_RAD(CM_BEM, CM_PIT_N3, CM_PIT_N3_REF, f"C_D_{test}Nb{Nb}_{layout}", "Added_Mass", ij, titre, vect_distance, i_w)


#################################################################
#################################################################
################### EVANESCENT PROBLEM ########################
#################################################################
#################################################################
Ne=[1, 6]
DOF=[1, 3, 5]
indice=[(1, 1), (3,9), (1,7), (5,11)]
CA_IT_files=[]
CM_IT_files=[]
Fex_IT_files=[]
Fex_ph_IT_files=[]
if COMP_E :   
    print("\n ------------ PLOT -----------------")
    print("\n Ne=", Ne)
    for i, E in enumerate(Ne) : 
        PIT_files_CM=[]
        PIT_files_CA=[]
        PIT_files_Fex=[]
        PIT_files_Fex_ph=[]
        for c, distance in enumerate(vect_distance):
            CM_PIT_E = f"{chemin_CAP}CapytaineIT_S_E{E}_Madd_{PIT_type}{beta_value}da{distance}.00.dat"
            CA_PIT_E = f"{chemin_CAP}CapytaineIT_S_E{E}_Crad_{PIT_type}{beta_value}da{distance}.00.dat"
            Fex_PIT_E = f"{chemin_CAP}CapytaineIT_S_E{E}_Fe_abs_{PIT_type}{beta_value}da{distance}.00.dat"
            Fex_phase_PIT_E = f"{chemin_CAP}CapytaineIT_S_E{E}_Fe_phase_{PIT_type}{beta_value}da{distance}.00.dat"
            PIT_files_CM.append(CM_PIT_E) 
            PIT_files_CA.append(CA_PIT_E) 
            PIT_files_Fex.append(Fex_PIT_E) 
            PIT_files_Fex_ph.append(Fex_phase_PIT_E) 
        CA_IT_files.append(PIT_files_CM)
        CM_IT_files.append(PIT_files_CA)
        Fex_IT_files.append(PIT_files_Fex)
        Fex_ph_IT_files.append(PIT_files_Fex_ph)

    indices_w=[3]
    for i_w in indices_w : 
        print("i_w : ", i_w)
        if GRAPHS_RAD : 
            print("\n indice=", indice )
            plot_RAD_E(CM_BEM, CA_PIT_N3, CA_IT_files, f"C_E_D_{test}Nb{Nb}_{layout}", "Added_Mass", indice, titre, vect_distance, Ne, i_w)
            plot_RAD_E(CA_BEM, CM_PIT_N3, CM_IT_files, f"C_E_D_{test}Nb{Nb}_{layout}", "Damping", indice, titre, vect_distance, Ne, i_w)
        if GRAPHS_Fex :
            print("\n DOF=", DOF )
            plot_Fex_E(Fex_BEM, Fex_PIT_N3, Fex_IT_files, f"C_E_D_{test}Nb{Nb}_{layout}_Fex_abs", DOF, titre, vect_distance, Ne, i_w)
            plot_Fex_E(Fex_BEM_phase, Fex_phase_PIT_N3,Fex_ph_IT_files, f"C_E_D_{test}Nb{Nb}_{layout}_Fex_ph", DOF, titre, vect_distance, Ne, i_w)
