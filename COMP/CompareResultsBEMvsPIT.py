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


def Trace_RAD(fileBEM, filePIT, file, nom, indice, titre, d):
    freq_BEM, Coef_BEM = read_RAD_BEM(fileBEM)
    freq_PIT, Coef_PIT = read_RAD_PIT(filePIT)

    plt.figure(figsize=(8, 6))
    plt.title(titre)
    plt.xlabel("Frequency")
    plt.ylabel(nom)
    plt.grid(True)
    labels = ["M_{11}", "M_{12}", "M_{21}", "M_{22}"]
    styles = ['o-', 's-', 'o', 's']  

    i=indice[0]-1
    j=indice[1]-1
    plt.plot(freq_BEM, [mat[i][j] for mat in Coef_BEM], 'o-', label=f"BEM_{indice[0]}_{indice[1]}")
    plt.plot(freq_PIT, [mat[i][j] for mat in Coef_PIT], 's--', label=f"PIT_{file}_{indice[0]}_{indice[1]}")
    plt.legend(loc='best', fontsize='small', frameon=True)
    plt.savefig(f"{file}_{nom}_{d}_M{indice[0]}_{indice[1]}.png", dpi=300)
    plt.close()
    return

def Trace_Fex(fileBEM, filePIT, filePIT_phase, file, num_dof, beta, titre, d):
    beta_BEM, freq_BEM, Coef_BEM_abs, Coef_BEM_ph = read_Fex_BEM(fileBEM)
    beta_PIT, freq_PIT, Coef_PIT = read_Fex_PIT(filePIT)
    beta_PIT, freq_PIT, Phase_PIT = read_Fex_PIT(filePIT_phase)

    if beta : # Several wave directions
        plt.figure(figsize=(10, 8))
        plt.title(titre)
        plt.xlabel("Wave Direction")
        plt.ylabel("|Fex|")
        i=0
        j=num_dof-1
        plt.plot(beta_BEM, [array[i, j] for array in Coef_BEM_abs], 'o-', label=f"BEM_doj{j}_w = {freq_BEM[i]:.2f}")  
        plt.plot(beta_PIT, [Coef_PIT[i][0][j] for i in range(19)], 's--', label=f"PIT_dof{j}_w = {freq_PIT[i]:.2f}")  
        plt.legend()
        plt.grid(True)
        plt.savefig(f"{file}_Fex_abs_{d}_B{num_dof}.png", dpi=300)
        plt.close()

        plt.figure(figsize=(10, 8))
        plt.title(titre)
        plt.xlabel("Wave Direction")
        plt.ylabel("phase(Fex)")
        i=0
        j=num_dof-1
        plt.plot(beta_BEM, [array[i, j] for array in Coef_BEM_ph], 'o-', label=f"BEM_dof{j}_w = {freq_BEM[i]:.2f}")  
        plt.plot(beta_PIT, [Phase_PIT[i][0][j] for i in range(19)], 's-', label=f"PIT_dof{j}_w = {freq_PIT[i]:.2f}")   

        plt.legend()
        plt.grid(True)
        plt.savefig(f"{file}_Fex_phase_{d}_B{num_dof}.png", dpi=300)
        plt.close()

    else : # ONE wave direction !!
        plt.figure(figsize=(10, 8))
        plt.title(titre)
        plt.xlabel("Frequency")
        plt.ylabel("|Fex|")
        i=0
        j=num_dof-1
        plt.plot(freq_BEM, Coef_BEM_abs[i][:, j], 'o-', label=f"BEM_dof{j+1}_beta = {beta_BEM[i]:.2f}°")   
        plt.plot(freq_PIT, [row[j] for row in Coef_PIT[i]], 's--', label=f"PIT_dof{j+1}_beta = {beta_PIT[i]:.2f}°")   
        plt.legend()
        plt.grid(True)
        plt.savefig(f"{file}_Fex_abs_{d}_B{num_dof}.png", dpi=300)
        plt.close()

        plt.figure(figsize=(10, 8))
        plt.title(titre)
        plt.xlabel("Frequency")
        plt.ylabel("phase(Fex)")
        i=0
        j=num_dof-1
        plt.plot(freq_BEM, Coef_BEM_ph[i][:, j], 'o-', label=f"BEM_dof{j+1}_beta = {beta_BEM[i]:.2f}°")   
        plt.plot(freq_PIT, [row[j] for row in Phase_PIT[i]], 's--', label=f"PIT_dof{j+1}_beta = {beta_PIT[i]:.2f}°")   
        plt.legend()
        plt.grid(True)
        plt.savefig(f"{file}_Fex_phase_{d}_B{num_dof}.png", dpi=300)
        plt.close()
    return


param_d=64

N3 = True
GRAPHS_beta=False

if N3 : 
    param_d=1
    limite=600
    GRAPHS_N3 = False
else : 
    limite=param_d
    GRAPHS_N3 = True

while param_d<=limite : 
    print("\n ------------Distance=", param_d, "-----------------")
    distance=f"d{param_d}"
    # BEM_file=f"BEM_barge_dof2_R6.36/BEM_barge_{distance}"
    # BEM_file=f"BEM_Nb1_dof6"
    test=""
    PIT_N3_file="PIT3_barge_dof6"
    beta_value=""
    Nb=2
    config="X"
    PIT_type=f"Nb{Nb}_{config}"
    # BEM_file=f"BEM_cylinderC/BEM_cylinderC_Nb{Nb}_{config}_{distance}"
    BEM_file=f"BEM_barge_dof6_R6.36/BEM_barge_{distance}"

    titre=f"Barge - {PIT_type}, Nw=20, Nbeta=11, dof=6, Ndir=1, {distance}"


    chemin="/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/"

    CM_BEM = f"{chemin}PFE_Nemoh/MyTestCases/{BEM_file}/results/CM.dat"
    CM_PIT_N3 = f"{chemin}PIT3_E/wec_inputs/{PIT_N3_file}/results/Global_{test}Madd_{PIT_type}_{beta_value}d{param_d}.00.dat"

    CA_BEM = f"{chemin}PFE_Nemoh/MyTestCases/{BEM_file}/results/CA.dat"
    CA_PIT_N3 = f"{chemin}PIT3_E/wec_inputs/{PIT_N3_file}/results/Global_{test}Crad_{PIT_type}_{beta_value}d{param_d}.00.dat"

    Fex_BEM = f"{chemin}PFE_Nemoh/MyTestCases/{BEM_file}/results/ExcitationForce.tec"  
    Fex_PIT_N3 = f"{chemin}PIT3_E/wec_inputs/{PIT_N3_file}/results/Global_{test}Fe_abs_{PIT_type}_{beta_value}d{param_d}.00.dat"
    Fex_phase_PIT_N3 = f"{chemin}PIT3_E/wec_inputs/{PIT_N3_file}/results/Global_{test}Fe_phase_{PIT_type}_{beta_value}d{param_d}.00.dat"

    if N3 : 
        # print("Norme infini - Norme L2 ")
        with open(f"EcartsRelatifs_{test}{PIT_type}.dat", "a") as f:
            f.write(f"{param_d}     ")

            print("BEM vs PIT_N3")
            print("Comparaison sur Added Mass")
            L_inf_error, L2_error = CompareResultsRAD(CM_BEM, CM_PIT_N3, "CM", test, param_d, PIT_type, Nb, 0.05, 10) # tol sur erreur, precision sur valeur
            # print(L_inf_error, L2_error)
            f.write(f"{L_inf_error:16.8f}     {L2_error:16.8f}      ")
            # print(" ")

            print("Comparaison sur Damping")
            L_inf_error, L2_error = CompareResultsRAD(CA_BEM, CA_PIT_N3, "CA", test, param_d, PIT_type, Nb, 0.05, 10)
            # print(L_inf_error, L2_error)
            f.write(f"{L_inf_error:16.8f}     {L2_error:16.8f}      ")
            # print(" ")   

            print("Comparaison sur Fex")
            L_inf_error_abs, L2_error_abs, L_inf_error_ph, L2_error_ph = CompareResultsFex(Fex_BEM, Fex_PIT_N3, Fex_phase_PIT_N3, test, param_d, PIT_type, Nb, 0.01, 10)
            # print("Ecart relatif entre les données pour |Fex| :")
            # print(L_inf_error_abs, L2_error_abs)
            f.write(f"{L_inf_error_abs:16.8f}     {L2_error_abs:16.8f}      ")
            # print("Ecart relatif entre les données pour phase(Fex) :")
            # print(L_inf_error_ph, L2_error_ph)
            f.write(f"{L_inf_error_ph:16.8f}      {L2_error_ph:16.8f}\n")

    if GRAPHS_N3 :
        # for i in range(2, 18, 6):
        for ind, i in enumerate([5, 11]) :
            print("dof", i)
            # Trace_RAD(CA_BEM, CA_PIT_N3, f"N3_{test}{PIT_type}", "Damping", [i+1,i+1], titre, distance)
            # Trace_RAD(CM_BEM, CM_PIT_N3, f"N3_{test}{PIT_type}", "Added_Mass", [i+1,i+1], titre, distance)
            Trace_Fex(Fex_BEM, Fex_PIT_N3, Fex_phase_PIT_N3, f"N3_{test}{PIT_type}", i, False, titre, distance)
            
        indices=[(5,5), (11,11), (5,11), (3,9), (3,5), (3,11), (5,9), (9,11)]
        # indices=[(3,9), (9,3), (3, 11), (3, 7)]
        for j, i in indices : 
            print("indice", j, i)
            Trace_RAD(CA_BEM, CA_PIT_N3, f"N3_{test}{PIT_type}", "Damping", [j,i], titre, distance)
            Trace_RAD(CM_BEM, CM_PIT_N3, f"N3_{test}{PIT_type}", "Added_Mass", [j,i], titre, distance)
    
    if GRAPHS_beta : 
        for i in range(0, 12, 2):
            print("dof", i+1)
            # Trace_RAD(CA_BEM, CA_PIT_N3, "N3", "Damping", [i+1,i+1], titre, distance)
            # Trace_RAD(CM_BEM, CM_PIT_N3, "N3", "Added_Mass", [i+1,i+1], titre, distance)
            Trace_Fex(Fex_BEM, Fex_PIT_N3, Fex_phase_PIT_N3, f"N3_{PIT_type}", i+1, True, titre, distance)
            
    
    param_d=param_d*2

