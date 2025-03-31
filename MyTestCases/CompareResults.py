import numpy as np
import re  # Pour utiliser les expressions régulières
import matplotlib.pyplot as plt


def read_RAD_BEM(fileRAD):
    frequencies = []  # Liste pour stocker les fréquences
    matrices = []     # Liste pour stocker les matrices
    current_matrix = []  # Temporaire pour stocker la matrice courante

    with open(fileRAD, 'r') as file:
        for line in file:
            line = line.strip()  # Enlever les espaces et les retours à la ligne

            # Si la ligne est vide ou commence par un caractère de commentaire (par exemple, '#')
            if not line or line.startswith('#'):
                continue  # Ignorer cette ligne

            try:
                if 'Nb' in line:
                    continue  # Ignorer cette ligne, elle ne contient pas de nombre
                # Essayer de convertir la ligne en un nombre (float) pour une fréquence
                frequency = float(line)
                if current_matrix:  # Si une matrice est en cours, l'ajouter aux matrices
                    matrices.append(current_matrix)
                    current_matrix = []  # Réinitialiser la matrice pour la prochaine fréquence
                frequencies.append(frequency)  # Ajouter la fréquence à la liste
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
        # print('BEM - Nw = ', len(frequencies))
        # print('BEM - Ndof x Nforce = ', len(matrices))
    return frequencies, matrices

def read_RAD_PIT(fileRAD):
    frequencies = []  # Liste pour stocker les fréquences
    matrices = []     # Liste pour stocker les matrices
    current_matrix = []  # Temporaire pour stocker la matrice courante

    previous_period = None  # Pour suivre la fréquence précédente

    with open(fileRAD, 'r') as file:
        for line in file:
            line = line.strip()  # Enlever les espaces et les retours à la ligne

            # Ignorer les lignes vides ou les commentaires
            if not line or line.startswith('#'):
                continue
            if 'Period' in line:
                    continue  # Ignorer cette ligne, elle ne contient pas de nombre

            # Convertir la ligne en une liste de valeurs flottantes
            values = list(map(float, line.split()))
            
            # La première valeur est la fréquence
            period = values[0]
            matrix_row = values[1:]  # Les autres valeurs forment une ligne de la matrice

            if previous_period is None or period == previous_period:
                # Si c'est la première ligne ou la même fréquence, ajouter la ligne à la matrice courante
                current_matrix.append(matrix_row)
            else:
                # Si la fréquence change, sauvegarder la matrice précédente et recommencer pour la nouvelle fréquence
                matrices.append(current_matrix)
                frequencies.append(previous_period)
                current_matrix = [matrix_row]  # Commencer une nouvelle matrice pour la nouvelle fréquence

            previous_period = period  # Mettre à jour la fréquence précédente

        # Ajouter la dernière matrice lue (si elle existe)
        if current_matrix:
            matrices.append(current_matrix)
            frequencies.append(previous_period)
    # print('PIT - Nw = ', len(frequencies))
    # print('PIT - Ndof x Nforce = ', len(matrices))
    return frequencies, matrices

# # Exemple d'utilisation
# fileRAD = 'test.dat'  # Remplace par le nom de ton fichier
# frequencies, matrices = read_RAD_PIT(fileRAD)

# # Afficher les résultats
# print("Fréquences lues :")
# print(frequencies)
# print("Matrices lues :")
# print(matrices)
# print("freq 1 : ", frequencies[0], matrices[0])

def CompareResultsRAD(fileBEM, filePIT) :
    freq_BEM, Coef_BEM = read_RAD_BEM(fileBEM)
    freq_PIT, Coef_PIT = read_RAD_PIT(filePIT)
    # print('BEM', Coef_BEM)
    # print('PIT', Coef_PIT)
    if np.all(np.isclose(freq_BEM, freq_PIT, atol=1e-4)):
        print(" OK même fréquences")
        if np.shape(Coef_BEM) != np.shape(Coef_PIT):
            print(np.shape(Coef_BEM), np.shape(Coef_PIT))
            raise ValueError("tailles différentes")
        # Calcul de l'erreur L2 (norme 2) entre les deux matrices
        else :
            errorL2 = np.sqrt(np.sum((np.array(Coef_BEM) - np.array(Coef_PIT))**2))
            errorINF = np.max(np.abs(np.array(Coef_BEM) - np.array(Coef_PIT)))

            normL2_BEM = np.sqrt(np.sum(np.array(Coef_BEM)**2))
            relative_errorL2 = errorL2 / normL2_BEM if normL2_BEM != 0 else np.nan

            max_BEM = np.max(np.abs(np.array(Coef_BEM)))
            relative_errorINF = errorINF / max_BEM if max_BEM != 0 else np.nan
    else :
        print("fréquences différentes")
        print('BEM', freq_BEM)
        print('PIT', freq_PIT)
        relative_errorL2=[]
        relative_errorINF=[]
    
    return relative_errorL2, relative_errorINF

def read_Fex_BEM(fileFex):
    beta_values = []  # Liste pour stocker les valeurs de beta
    matrices = []     # Liste pour stocker les matrices pour chaque beta
    matrices_abs = []
    matrices_phase = []
    frequencies = []  # Liste pour stocker les fréquences associées

    with open(fileFex, 'r') as f:
        beta = None  # Initialiser une variable pour beta
        matrix = []  # Matrice vide pour stocker les valeurs associées à chaque beta
        
        for line in f:
            line = line.strip()
            
            # Si la ligne contient un 'Zone', chercher le beta
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
                        matrix = []  # Réinitialiser la matrice pour le prochain beta
                    
                    beta = float(match.group(1))  # Récupérer la nouvelle valeur de beta

            # Si la ligne n'est pas vide et qu'on est dans une matrice
            elif line and beta is not None:
                # Convertir la ligne en une liste de valeurs flottantes
                values = list(map(float, line.split()))
                if values:
                    matrix.append(values[1:])  # Les autres colonnes sont les valeurs de la matrice
                    if values[0] not in frequencies:
                        frequencies.append(values[0])  # Ajouter la fréquence à la liste si elle n'y est pas

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
        # print(matrices[0][0], matrices_abs[0][0], matrices_phase[0][0])
    # Retourner les valeurs de beta, les fréquences et les matrices
    # print('BEM - Nbeta = ', len(beta_values))
    # print('BEM - Nforce = ', len(matrices_abs))
    return beta_values, frequencies, matrices_abs, matrices_phase

def read_Fex_PIT(fileFex):
    beta_values = []  # Liste pour stocker les valeurs de beta
    matrices = []     # Liste pour stocker les matrices pour chaque beta
    frequencies = []  # Liste pour stocker les fréquences associées
    current_matrix = []
    previous_beta = None

    with open(fileFex, 'r') as file:
        for line in file:
            line = line.strip()  # Enlever les espaces et les retours à la ligne

            # Ignorer les lignes vides ou les commentaires
            if not line or line.startswith('#'):
                continue
            if 'Period' in line:
                    continue  # Ignorer cette ligne, elle ne contient pas de nombre

            # Convertir la ligne en une liste de valeurs flottantes
            values = list(map(float, line.split()))

            # La première valeur est la fréquence
            beta = values[0]
            period = values[1]
            matrix_row = values[2:]  # Les autres valeurs forment une ligne de la matrice

            if previous_beta is None or beta == previous_beta:
                # Si c'est la première ligne ou la même fréquence, ajouter la ligne à la matrice courante
                if period not in frequencies:
                    frequencies.append(period)
                current_matrix.append(matrix_row)
            else:
                # Si la fréquence change, sauvegarder la matrice précédente et recommencer pour la nouvelle fréquence
                matrices.append(current_matrix)
                beta_values.append(previous_beta)
                current_matrix = [matrix_row]  # Commencer une nouvelle matrice pour la nouvelle fréquence
                
            previous_beta = beta  # Mettre à jour la fréquence précédente

        # Ajouter la dernière matrice lue (si elle existe)
        if current_matrix:
            matrices.append(current_matrix)
            beta_values.append(previous_beta)
    # print('PIT - Nbeta = ', len(beta_values))
    # print('PIT - Nforce = ', len(matrices))
    return beta_values, frequencies, matrices

def CompareResultsFex(fileBEM, filePIT) :
    beta_BEM, freq_BEM, Coef_BEM_abs, Coef_BEM_ph = read_Fex_BEM(fileBEM)
    beta_PIT, freq_PIT, Coef_PIT = read_Fex_PIT(filePIT)

    # if np.all(np.isclose(freq_BEM, freq_PIT, atol=1e-2)) and np.all(np.isclose(beta_BEM, beta_PIT, atol=1e-2)):
    if len(freq_BEM)==len(freq_PIT) and len(beta_BEM)==len(beta_PIT):
        print("même fréquences")
        if np.shape(Coef_BEM_abs) != np.shape(Coef_PIT):
            print(np.shape(Coef_BEM_abs), np.shape(Coef_PIT))
            raise ValueError("tailles différentes")
        # Calcul de l'erreur L2 (norme 2) entre les deux matrices
        else :
            errorL2 = np.sqrt(np.sum((np.array(Coef_BEM_abs) - np.array(Coef_PIT))**2))
            errorINF = np.max(np.abs(np.array(Coef_BEM_abs) - np.array(Coef_PIT)))

            normL2_BEM = np.sqrt(np.sum(np.array(Coef_BEM_abs)**2))
            relative_errorL2 = errorL2 / normL2_BEM if normL2_BEM != 0 else np.nan

            max_BEM = np.max(np.abs(np.array(Coef_BEM_abs)))
            relative_errorINF = errorINF / max_BEM if max_BEM != 0 else np.nan
    else :
        print("fréquences différentes")
        print('BEM', freq_BEM)
        print('PIT', freq_PIT)
        relative_errorL2=[]
        relative_errorINF=[]
    
    return relative_errorL2, relative_errorINF


CM_BEM = "BEM_Nb2_dof1_d1/results/CM.dat"
CM_PIT = "/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/DTOcean_raw/wec_inputs/PIT_Nb2_dof1/results/Global_Madd_d1.00.dat"

CA_BEM = "BEM_Nb2_dof1_d1/results/CA.dat"
CA_PIT = "/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/DTOcean_raw/wec_inputs/PIT_Nb2_dof1/results/Global_Crad_d1.00.dat"

Fex_BEM = "BEM_Nb2_dof1_d1/results/ExcitationForce.tec"  
Fex_PIT = "/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/DTOcean_raw/wec_inputs/PIT_Nb2_dof1/results/Global_Fe_d1.00.dat"

CM = True
CA = False
Fex = False

if CM : 
    print("Comparaison sur Added Mass")
    L_inf_error, L2_error = CompareResultsRAD(CM_BEM, CM_PIT)

    print("Ecart relatif entre les données pour Added Mass :")
    print("Norme infini - Norme L2 ")
    print(L_inf_error, L2_error)
    print(" ")

if CA : 
    print("Comparaison sur Damping")
    L_inf_error, L2_error = CompareResultsRAD(CA_BEM, CA_PIT)

    print("Ecart relatif entre les données pour Damping :")
    print("Norme infini - Norme L2 ")
    print(L_inf_error, L2_error)
    print(" ")

if Fex : 
    print("Comparaison sur Fex")
    L_inf_error, L2_error = CompareResultsFex(Fex_BEM, Fex_PIT)

    print("Ecart relatif entre les données pour Fex :")
    print("Norme infini - Norme L2 ")
    print(L_inf_error, L2_error)


def Trace_RAD(fileBEM, filePIT, nom, param):
    freq_BEM, Coef_BEM = read_RAD_BEM(fileBEM)
    freq_PIT, Coef_PIT = read_RAD_PIT(filePIT)

    plt.figure(figsize=(8, 6))
    plt.title(param)
    plt.xlabel("Frequency")
    plt.ylabel(nom)
    plt.grid(True)

    # # Tracé des courbes
    # plt.plot(freq_BEM, Coef_BEM[:][0][0], marker='o', linestyle='-', label="BEM - M_{11}")
    # plt.plot(freq_BEM, Coef_BEM[:][0][1], marker='o', linestyle='-', label="BEM - M_{12}")
    # plt.scatter(freq_BEM, Coef_BEM[:][1][0], label="BEM - M_{21}")
    # plt.scatter(freq_BEM, Coef_BEM[:][1][1], label="BEM - M_{22}")

    # plt.plot(freq_PIT, Coef_PIT[:][0][0], marker='o', linestyle='-', label="PIT - M_{11}")
    # plt.plot(freq_PIT, Coef_PIT[:][0][1], marker='o', linestyle='-', label="PIT - M_{12}")
    # plt.scatter(freq_PIT, Coef_PIT[:][1][0], label="PIT - M_{21}")
    # plt.scatter(freq_PIT, Coef_PIT[:][1][1], label="PIT - M_{22}")
    # Labels des coefficients
    labels = ["M_{11}", "M_{12}", "M_{21}", "M_{22}"]
    styles = ['o-', 's-', 'o', 's']  

    # Tracé des courbes pour BEM
    # for idx, (i, j) in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)]):
    for idx, (i, j) in enumerate([(0, 0), (0, 1)]):
        # Accède aux éléments dans chaque matrice de Coef_BEM pour chaque fréquence
        plt.plot(freq_BEM, [mat[i][j] for mat in Coef_BEM], styles[idx], label=f"BEM - {labels[idx]}")

    # Tracé des points pour PIT
    # for idx, (i, j) in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)]):
    for idx, (i, j) in enumerate([(0, 0), (0, 1)]):
        # Accède aux éléments dans chaque matrice de Coef_PIT pour chaque fréquence
        plt.plot(freq_PIT, [mat[i][j] for mat in Coef_PIT], styles[idx], label=f"PIT - {labels[idx]}")

    # Légende et affichage
    plt.legend(loc='best', fontsize='small', frameon=True)
    plt.savefig(f"{nom}.png", dpi=300)
    # plt.show()
    return

def Trace_Fex(fileBEM, filePIT, Nombre, param):
    beta_BEM, freq_BEM, Coef_BEM_abs, Coef_BEM_ph = read_Fex_BEM(fileBEM)
    beta_PIT, freq_PIT, Coef_PIT = read_Fex_PIT(filePIT)

    # Tracé des matrices absolues
    plt.figure(figsize=(10, 8))
    plt.title(param)
    plt.xlabel("Frequency")
    plt.ylabel("|Fex|")
    
    # On suppose que les matrices sont une liste de matrices pour chaque valeur de beta
    # for i, beta in enumerate(beta_BEM):
    #     plt.plot(freq_BEM, Coef_BEM_abs[i][:, 0], 'o-', label=f"BEM_beta = {beta:.2f}°")  # Exemple pour la première colonne (à adapter si nécessaire)
    # for i, beta in enumerate(beta_PIT):
    #     plt.plot(freq_BEM, [row[0] for row in Coef_PIT[i]], 'o-', label=f"PIT_beta = {beta:.2f}°")  # Exemple pour la première colonne (à adapter si nécessaire)
    for i in range(Nombre[0],Nombre[1]) :
        plt.plot(freq_BEM, Coef_BEM_abs[i][:, 0], 'o-', label=f"BEM_beta = {beta_BEM[i]:.2f}°")  # Exemple pour la première colonne (à adapter si nécessaire)
        plt.plot(freq_BEM, [row[0] for row in Coef_PIT[i]], 's-', label=f"PIT_beta = {beta_PIT[i]:.2f}°")  # Exemple pour la première colonne (à adapter si nécessaire)

    plt.legend()
    plt.grid(True)
    plt.savefig(f"Fex_abs_N{Nombre[0]}-{Nombre[1]}.png", dpi=300)
    # plt.show()
    return

GRAPHS = False
if GRAPHS :
    Trace_RAD(CA_BEM, CA_PIT, "Damping_c2", "Nw=20, Nbeta=11, dof=1, Ndir=5, d=1")
    Trace_RAD(CM_BEM, CM_PIT, "Added_Mass_c2", "Nw=20, Nbeta=11, dof=1, Ndir=5, d=1")
    Trace_Fex(Fex_BEM, Fex_PIT, [3, 5], "Nw=20, Nbeta=11, dof=1, Ndir=5, d=1")


