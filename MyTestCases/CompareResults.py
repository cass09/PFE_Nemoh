import numpy as np
import re  # Pour utiliser les expressions régulières

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
        ONE_BODY = True
        if ONE_BODY:
            matrices = [[[value]] for value in frequencies[1::2]]  # Indices impairs
            frequencies = frequencies[::2]   # Indices pairs
        print('BEM - Nw = ', len(frequencies))
        print('BEM - Ndof x Nforce = ', len(matrices))
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
    print('PIT - Nw = ', len(frequencies))
    print('PIT - Ndof x Nforce = ', len(matrices))
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
    else :
        print("fréquences différentes")
        print('BEM', freq_BEM)
        print('PIT', freq_PIT)
        errorL2=[]
        errorINF=[]
    
    return errorL2, errorINF

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
                        beta_values.append(beta*180/np.pi) 
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
            beta_values.append(beta*180/np.pi)
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
    print('BEM - Nbeta = ', len(beta_values))
    print('BEM - Nforce = ', len(matrices_abs))
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
    print('PIT - Nbeta = ', len(beta_values))
    print('PIT - Nforce = ', len(matrices))
    return beta_values, frequencies, matrices

def CompareResultsFex(fileBEM, filePIT) :
    beta_BEM, freq_BEM, Coef_BEM_abs, Coef_BEM_ph = read_Fex_BEM(Fex_BEM)
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
    else :
        print("fréquences différentes")
        print('BEM', freq_BEM)
        print('PIT', freq_PIT)
        errorL2=[]
        errorINF=[]
    
    return errorL2, errorINF


CM_BEM = "BEM_Nb1_dof1/results/CM.dat"
CM_PIT = "PythonResults/BARGE/results/Global_Madd_d0.00.dat"

CA_BEM = "BEM_Nb1_dof1/results/CA.dat"
CA_PIT = "PythonResults/BARGE/results/Global_Crad_d0.00.dat"

Fex_BEM = "BEM_Nb1_dof1/results/DiffractionForce.tec"  
Fex_PIT = "PythonResults/BARGE/results/Global_Fe_d0.00.dat"

CM = False
CA = True
Fex = True

if CM : 
    print("Comparaison sur Added Mass")
    L_inf_error, L2_error = CompareResultsRAD(CM_BEM, CM_PIT)

    print("Erreur entre les données pour Added Mass :")
    print("Norme infini - Norme L2 ")
    print(L_inf_error, L2_error)
    print(" ")

if CA : 
    print("Comparaison sur Damping")
    L_inf_error, L2_error = CompareResultsRAD(CA_BEM, CA_PIT)

    print("Erreur entre les données pour Damping :")
    print("Norme infini - Norme L2 ")
    print(L_inf_error, L2_error)
    print(" ")

# Exemple d'utilisation de la fonction
if Fex : 
    print("Comparaison sur Fex")
    L_inf_error, L2_error = CompareResultsFex(Fex_BEM, Fex_PIT)

    print("Erreur entre les données pour Fex :")
    print("Norme infini - Norme L2 ")
    print(L_inf_error, L2_error)

test = False
if test : 
    # Exemple de valeurs pour Coef_BEM et Coef_PIT
    Coef_BEM = [0.3000,  0.266654E+04, 0.3895,  0.641170E+04, 0.4789,  0.134547E+05, 0.5684,  0.260890E+05, 0.6579,  0.482161E+05, 0.7474,  0.863001E+05, 0.8368,  0.150076E+06, 0.9263,  0.251168E+06, 1.0158,  0.395789E+06, 1.1053,  0.570975E+06, 1.1947,  0.739102E+06, 1.2842,  0.858127E+06, 1.3737,  0.909674E+06, 1.4632,  0.902934E+06, 1.5526,  0.858128E+06, 1.6421,  0.793974E+06, 1.7316,  0.722598E+06, 1.8211,  0.650887E+06, 1.9105,  0.581523E+06, 2.0000]
    Coef_PIT = [0.3, 2.658450e+03, 0.3895,   6.385620e+03, 0.4789 , 1.338850e+04  ,0.5684 ,  2.593260e+04  ,0.6579 , 4.779290e+04  ,0.7474,  8.540450e+04  ,0.8368,   1.482600e+05  ,0.9263  , 2.474380e+05  ,1.01578,  3.903800e+05  ,1.1053 ,  5.625200e+05  ,1.1947,   7.293760e+05  ,1.2842  , 8.471140e+05  ,1.3737 ,  8.986340e+05  ,1.4632,   8.914700e+05  ,1.5526  , 8.460750e+05  ,1.6421 ,  7.823980e+05  ,1.7316 ,  7.092290e+05  ,1.8211  , 6.367760e+05  ,1.9105 ,  5.681760e+05 , 2.0 ]
    # Coef_BEM = [1.2e3, 4.5e2, 6.7e1, 8.9e-1, 1.1e-2]
    # Coef_PIT = [0.12e4, 4.4e2, 6.6e1, 9.0e-1, 1.0e-2]
    # Calcul de l'erreur L2 (norme euclidienne)
    errorL2 = np.sqrt(np.sum((np.array(Coef_BEM) - np.array(Coef_PIT))**2))

    # Calcul de l'erreur infinie (valeur absolue maximale)
    errors = np.abs(np.array(Coef_BEM) - np.array(Coef_PIT))
    errorINF = np.max(errors)

    # Trouver l'indice de l'erreur maximale
    max_error_index = np.argmax(errors)

    # Affichage des résultats en notation scientifique
    print(f"Erreur L2: {errorL2:.3e}")  # Format scientifique avec 3 décimales
    print(f"Erreur INF: {errorINF:.3e}")  # Format scientifique avec 3 décimales

    # Afficher les valeurs pour lesquelles l'erreur est maximale
    print(f"Valeur de Coef_BEM à l'indice max: {Coef_BEM[max_error_index]}")
    print(f"Valeur de Coef_PIT à l'indice max: {Coef_PIT[max_error_index]}")
    print(f"Erreur maximale: {errors[max_error_index]:.3e}")