import numpy as np
import re  # Pour utiliser les expressions régulières
import matplotlib.pyplot as plt


distance="Nb1"
CM = False
CA = False
Fex = False
folder="ConvergenceMesh"
mesh = "barge"
N1_file=f"{folder}/BEM_{mesh}_Nb1_MeshOG"
N2_file=f"{folder}/BEM_{mesh}_Nb1_Mesh2"
N3_file=f"{folder}/BEM_{mesh}_Nb1_Mesh3"
N4_file=f"{folder}/BEM_{mesh}_Nb1_Mesh4"
N5_file=f"{folder}/BEM_{mesh}_Nb1_Mesh5"
# NOMS=["Np=250", "Np=404", "Np=584", "Np=1000", "Np=1308", "Np=1626"]
NOMS=["Np=300", "Np=624", "Np=1128", "Np=1764", "Np=2448"]
file=[N1_file, N2_file, N3_file, N4_file, N5_file]
Nmesh=5
GRAPHS = False
CoF=True
# titre=f"Nw=20, Nbeta=11, dof=6, Ndir=1, {distance}"
titre=f"Convergence study on {mesh} mesh"
COMP1=N1_file
COMP2=N2_file

chemin= "/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/PFE_Nemoh/MyTestCases/"
CM_BEM = f"{chemin}{COMP1}/results/CM.dat"
CM_N2 = f"{chemin}{COMP2}/results/CM.dat"

CA_BEM = f"{chemin}{COMP1}/results/CA.dat"
CA_N2 = f"{chemin}{COMP2}/results/CA.dat"

Fex_BEM = f"{chemin}{COMP1}/results/ExcitationForce.tec"  
Fex_N2 = f"{chemin}{COMP2}/results/ExcitationForce.tec"

CMfiles=[]
CAfiles=[]
Fefiles=[]
for i in range(Nmesh) :
    CMfiles.append(f"{chemin}{file[i]}/results/CM.dat")
    CAfiles.append(f"{chemin}{file[i]}/results/CA.dat")
    Fefiles.append(f"{chemin}{file[i]}/results/ExcitationForce.tec")

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


def CompareResultsRAD(fileBEM, fileN2) :
    freq_BEM, Coef_BEM = read_RAD_BEM(fileBEM)
    freq_N2, Coef_N2 = read_RAD_BEM(fileN2)
    # print('BEM', Coef_BEM)
    # print('PIT', Coef_N2)
    if np.all(np.isclose(freq_BEM, freq_N2, atol=1e-4)):
        print(" OK même fréquences")
        if np.shape(Coef_BEM) != np.shape(Coef_N2):
            print(np.shape(Coef_BEM), np.shape(Coef_N2))
            raise ValueError("tailles différentes")
        # Calcul de l'erreur L2 (norme 2) entre les deux matrices
        else :
            errorL2 = np.sqrt(np.sum((np.array(Coef_BEM) - np.array(Coef_N2))**2))
            errorINF = np.max(np.abs(np.array(Coef_BEM) - np.array(Coef_N2)))

            normL2_BEM = np.sqrt(np.sum(np.array(Coef_BEM)**2))
            relative_errorL2 = errorL2 / normL2_BEM if normL2_BEM != 0 else np.nan

            max_BEM = np.max(np.abs(np.array(Coef_BEM)))
            relative_errorINF = errorINF / max_BEM if max_BEM != 0 else np.nan
    else :
        print("fréquences différentes")
        print('BEM', freq_BEM)
        print('PIT', freq_N2)
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

def CompareResultsFex(fileBEM, fileN2) :
    beta_BEM, freq_BEM, Coef_BEM_abs, Coef_BEM_ph = read_Fex_BEM(fileBEM)
    beta_N2, freq_N2, Coef_N2_abs, Coef_N2_ph = read_Fex_BEM(fileN2)

    # if np.all(np.isclose(freq_BEM, freq_N2, atol=1e-2)) and np.all(np.isclose(beta_BEM, beta_N2, atol=1e-2)):
    if len(freq_BEM)==len(freq_N2) and len(beta_BEM)==len(beta_N2):
        print("même fréquences")
        if np.shape(Coef_BEM_abs) != np.shape(Coef_N2_abs):
            print(np.shape(Coef_BEM_abs), np.shape(Coef_N2_abs))
            raise ValueError("tailles différentes")
        if np.shape(Coef_BEM_ph) != np.shape(Coef_N2_ph):
            print(np.shape(Coef_BEM_ph), np.shape(Coef_N2_ph))
            raise ValueError("tailles différentes")
        # Calcul de l'erreur L2 (norme 2) entre les deux matrices
        else :
            errorL2 = np.sqrt(np.sum((np.array(Coef_BEM_abs) - np.array(Coef_N2_abs))**2))
            errorINF = np.max(np.abs(np.array(Coef_BEM_abs) - np.array(Coef_N2_abs)))

            normL2_BEM = np.sqrt(np.sum(np.array(Coef_BEM_abs)**2))
            relative_errorL2_abs = errorL2 / normL2_BEM if normL2_BEM != 0 else np.nan

            max_BEM = np.max(np.abs(np.array(Coef_BEM_abs)))
            relative_errorINF_abs = errorINF / max_BEM if max_BEM != 0 else np.nan

            errorL2_p = np.sqrt(np.sum((np.array(Coef_BEM_ph) - np.array(Coef_N2_ph))**2))
            errorINF_p = np.max(np.abs(np.array(Coef_BEM_ph) - np.array(Coef_N2_ph)))

            normL2_BEM_p = np.sqrt(np.sum(np.array(Coef_BEM_ph)**2))
            relative_errorL2_phase = errorL2_p / normL2_BEM_p if normL2_BEM_p != 0 else np.nan

            max_BEM_p = np.max(np.abs(np.array(Coef_BEM_ph)))
            relative_errorINF_phase = errorINF_p / max_BEM_p if max_BEM_p != 0 else np.nan
    else :
        print("fréquences différentes")
        print('BEM', freq_BEM)
        print('PIT', freq_N2)
        relative_errorL2_abs=[]
        relative_errorL2_phase=[]
        relative_errorINF_abs=[]
        relative_errorINF_phase=[]
    
    return relative_errorL2_abs, relative_errorINF_abs, relative_errorL2_phase, relative_errorINF_phase


def Trace_RAD(files, nom, indice, titre, distance, NOMS, mesh):
    
    couleurs = ['b', 'C1', 'g', 'r', 'gold', 'C4', 'c', 'm', 'y', 'k']  
    plt.figure(figsize=(8, 6))
    plt.title(f"{titre}")
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
        frequencies, Coef = read_RAD_BEM(file)
    # freq_N2, Coef_N2 = read_RAD_BEM(fileN2)
        plt.plot(frequencies, [mat[i][j] for mat in Coef], 'o-', color=couleurs[k], linewidth=0.8, markersize=4, label=f"{NOMS[k]}")
    # plt.plot(freq_N2, [mat[i][j] for mat in Coef_N2], 's-', label=f"{NOMS[0]}_{indice[0]}{indice[1]}")

    # Légende et affichage
    plt.legend(loc='best', fontsize='small', frameon=True)
    plt.savefig(f"N_{mesh}_{nom}_{distance}_M{indice[0]}{indice[1]}.png", dpi=300)
    # plt.show()
    return

def Trace_Fex(files, num_dof, num_beta, titre, distance, NOMS, mesh):
    # beta_BEM, freq_BEM, Coef_BEM_abs, Coef_BEM_ph = read_Fex_BEM(fileBEM)
    # beta_N2, freq_N2, Coef_N2_abs, Coef_N2_ph = read_Fex_BEM(fileN2)
    i=num_beta-1
    # Tracé des matrices absolues
    couleurs = ['b', 'C1', 'g', 'r', 'gold', 'C4', 'c', 'm', 'y', 'k'] 
    plt.figure(figsize=(10, 8))
    plt.title(f"{titre}")
    plt.xlabel("Frequency (rad/s)")
    plt.ylabel(f"|Fex| (N) - DOF {num_dof}")
    
    
    j=num_dof-1

    for k, file in enumerate(files) : 
        beta, freq, Coef_abs, Coef_ph = read_Fex_BEM(file)
        # print(k, Coef_abs[i][:, j], Coef_ph[i][:, j])
        plt.plot(freq, Coef_abs[i][:, j], 'o-', color=couleurs[k], label=f"{NOMS[k]}")  # Exemple pour la première colonne (à adapter si nécessaire)

    plt.legend()
    plt.grid(True)
    plt.savefig(f"N_{mesh}_Fex_abs_{distance}_b{num_dof}.png", dpi=300)

    # Tracé des matrices absolues
    i=num_beta-1
    plt.figure(figsize=(10, 8))
    plt.title(f"{titre}")
    plt.xlabel("Frequency (rad/s)")
    plt.ylabel(f"phase(Fex) (N) - DOF {num_dof}")

    
    j=num_dof-1
    for k, file in enumerate(files) : 
        beta, freq, Coef_abs, Coef_ph = read_Fex_BEM(file)
        plt.plot(freq, Coef_ph[i][:, j], 'o-', color=couleurs[k], label=f"{NOMS[k]}")  # Exemple pour la première colonne (à adapter si nécessaire)

    plt.legend()
    plt.grid(True)
    plt.savefig(f"N_{mesh}_Fex_phase_{distance}_b{num_dof}.png", dpi=300)
    # plt.show()
    return

def Calcul_COF(REF, files, Nom_data, dof_vect, NOMS, mesh):
    max_CoF=np.zeros((len(dof_vect), len(files)))
    for k in range(len(files)) : 
        frequencies, COEF = read_RAD_BEM(files[k])
        frequencies, COEF_REF = read_RAD_BEM(REF)
        for j, dof in enumerate(dof_vect) :
            CoF_value = []
            for i in range(len(frequencies)) :
                with np.errstate(divide='ignore', invalid='ignore'):
                    ratio = np.where(COEF_REF[i][dof-1][dof - 1] != 0, COEF[i][dof-1][dof - 1] / COEF_REF[i][dof-1][dof - 1], 0.0)
                CoF_value.append(ratio)
            max_CoF[j, k] = np.max(CoF_value)  # ou np.mean(CoF_value) si tu préfères la moyenne
        Trace_COF(max_CoF, Nom_data, NOMS, dof_vect, len(files), mesh)

def Calcul_COF_Fex(REF, files, dof_vect, NOMS, mesh):
    max_CoF=np.zeros((len(dof_vect), len(files)))
    max_CoF_ph=np.zeros((len(dof_vect), len(files)))
    for k in range(len(files)) : 
        beta, frequencies, Coef_abs, Coef_ph = read_Fex_BEM(files[k])
        beta, frequencies, Coef_abs_REF, Coef_ph_REF = read_Fex_BEM(REF)
        for j, dof in enumerate(dof_vect) :
            CoF_value = []
            CoF_value_ph = []
            for i in range(len(frequencies)) :
                # CoF_value[i]=COEF[i,dof-1]/COEF_REF[i,dof-1]
                # Pour éviter une division par zéro
                with np.errstate(divide='ignore', invalid='ignore'):
                    ratio = np.where(Coef_abs_REF[0][i, dof - 1] != 0, Coef_abs[0][i,dof - 1] / Coef_abs_REF[0][i,dof - 1], 0.0)
                    ratio_ph = np.where(Coef_ph_REF[0][i, dof - 1] != 0, Coef_ph[0][i,dof - 1] / Coef_ph_REF[0][i,dof - 1], 0.0)
                CoF_value.append(ratio)
                CoF_value_ph.append(ratio_ph)
            max_CoF[j, k] = np.max(CoF_value)  # ou np.mean(CoF_value) si tu préfères la moyenne
            max_CoF_ph[j, k] = np.max(CoF_value_ph)  # ou np.mean(CoF_value) si tu préfères la moyenne
        Trace_COF(max_CoF, "|Fex|", NOMS, dof_vect, len(files), mesh)
        Trace_COF(max_CoF_ph, "phase(Fex)", NOMS, dof_vect, len(files), mesh)

def Trace_COF(CoF, data, NOMS, dof_vect, Nfiles, mesh) :
    with open(f"CoefficientOfFidelity_{mesh}_{data}.dat", "w") as f:
        # En-tête avec les noms de fichiers
        f.write(f"# Coefficient of Fidelity ({data})\n")
        f.write("DOF\\Mesh \n"+"".join([f"\t{nom}" for nom in NOMS]) + "\n")

        for j, dof in enumerate(dof_vect):
            f.write(f"{dof}")
            for k in range(Nfiles):
                f.write(f"\t{CoF[j, k]:.4f}")
            f.write("\n")

if CM : 
    print("Comparaison sur Added Mass")
    L_inf_error, L2_error = CompareResultsRAD(CM_BEM, CM_N2)

    print("Ecart relatif entre les données pour Added Mass :")
    print("Norme infini - Norme L2 ")
    print(L_inf_error, L2_error)
    print(" ")

if CA : 
    print("Comparaison sur Damping")
    L_inf_error, L2_error = CompareResultsRAD(CA_BEM, CA_N2)

    print("Ecart relatif entre les données pour Damping :")
    print("Norme infini - Norme L2 ")
    print(L_inf_error, L2_error)
    print(" ")

if Fex : 
    print("Comparaison sur Fex")
    L_inf_error_abs, L2_error_abs, L_inf_error_ph, L2_error_ph = CompareResultsFex(Fex_BEM, Fex_N2)

    print("Ecart relatif entre les données pour |Fex| :")
    print("Norme infini - Norme L2 ")
    print(L_inf_error_abs, L2_error_abs)
    print("Ecart relatif entre les données pour phase(Fex) :")
    print("Norme infini - Norme L2 ")
    print(L_inf_error_ph, L2_error_ph)


if GRAPHS :
    for ind, dof in enumerate([1, 3, 5]) : 
        Trace_RAD(CAfiles, "Damping", [dof,dof], titre, distance, NOMS, mesh)
        Trace_RAD(CMfiles, "Added_Mass", [dof,dof], titre, distance, NOMS, mesh)
        Trace_Fex(Fefiles, dof, 1, titre, distance, NOMS, mesh)

if CoF:
    dof=[1, 3, 5]
    Calcul_COF(CAfiles[-1], CAfiles, "Added_Mass", dof, NOMS, mesh)
    Calcul_COF(CMfiles[-1], CMfiles, "Damping", dof, NOMS, mesh)
    Calcul_COF_Fex(Fefiles[-1], Fefiles, dof, NOMS, mesh)
