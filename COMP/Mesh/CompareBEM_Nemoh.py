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
    # plt.xlim([0, 2])
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
        # Supprimer l'élément d'indice 15 (16ᵉ élément)
        # frequencies = [f for i, f in enumerate(frequencies) if i != 15]
        # Coef = [mat for i, mat in enumerate(Coef) if i != 15]

        plt.plot(frequencies, [mat[i][j] for mat in Coef], 'o-', color=couleurs[k], linewidth=0.8, markersize=4, label=f"{NOMS[k]}")
    # plt.plot(freq_N2, [mat[i][j] for mat in Coef_N2], 's-', label=f"{NOMS[0]}_{indice[0]}{indice[1]}")

    # Légende et affichage
    plt.legend(loc='best', fontsize='small', frameon=True)
    plt.savefig(f"N_{mesh}_{nom}_{distance}_M{indice[0]}{indice[1]}.pdf")
    # plt.show()
    return

def Trace_Fex(files, num_dof, num_beta, titre, distance, NOMS, mesh):
    i=num_beta-1
    couleurs = ['b', 'C1', 'g', 'r', 'gold', 'C4', 'c', 'm', 'y', 'k'] 
    plt.figure(figsize=(10, 8))
    plt.title(f"{titre}")
    plt.xlabel("Frequency (rad/s)")
    plt.ylabel(f"|Fex| (N) - DOF {num_dof}")
    j=num_dof-1
    for k, file in enumerate(files) : 
        beta, freq, Coef_abs, Coef_ph = read_Fex_BEM(file)
        # Supprimer l'indice 15
        # freq = [f for i, f in enumerate(freq) if i != 15]
        # Coef_abs = [np.delete(mat, 15, axis=0) for mat in Coef_abs]
        # Convertir la matrice i en array NumPy
        mat_i = np.array(Coef_abs[i])
        plt.plot(freq, mat_i[:, j], 'o-', color=couleurs[k], label=f"{NOMS[k]}")
    plt.legend()
    plt.grid(True)
    plt.savefig(f"N_{mesh}_Fex_abs_{distance}_b{num_dof}.pdf")
    plt.close()
    i=num_beta-1
    plt.figure(figsize=(10, 8))
    plt.title(f"{titre}")
    plt.xlabel("Frequency (rad/s)")
    plt.ylabel(f"phase(Fex) (N) - DOF {num_dof}")    
    j=num_dof-1
    for k, file in enumerate(files) : 
        beta, freq, Coef_abs, Coef_ph = read_Fex_BEM(file)
        # freq = [f for i, f in enumerate(freq) if i != 15]
        # Coef_ph = [np.delete(mat, 15, axis=0) for mat in Coef_ph]
        mat_i = np.array(Coef_ph[i])
        plt.plot(freq, mat_i[:, j], 'o-', color=couleurs[k], label=f"{NOMS[k]}")
    plt.legend()
    plt.grid(True)
    plt.savefig(f"N_{mesh}_Fex_phase_{distance}_b{num_dof}.pdf")
    plt.close()
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

# def Calcul_GOF(REF, files, Nom_data, dof_vect, NOMS, mesh):
#     min_CoF=np.zeros((len(dof_vect), len(files)))
#     mean_CoF=np.zeros((len(dof_vect), len(files)))
#     for k in range(len(files)) : 
#         frequencies, COEF = read_RAD_BEM(files[k])
#         frequencies, COEF_REF = read_RAD_BEM(REF)
#         for j, dof in enumerate(dof_vect) :
#             CoF_value = []
#             for i in range(len(frequencies)) :
#                 value=COEF[i][dof-1][dof - 1]
#                 value_ref=COEF_REF[i][dof-1][dof - 1]
#                 ratio = np.abs(1-np.abs((value-value_ref)/value_ref))
#                 if ratio <0.9:
#                     print("dof=",dof, " - mesh=", NOMS[k]," - w=", frequencies[i], "GoF<90%", ratio)
#                 CoF_value.append(ratio)
#                 mean_CoF[j,k]=mean_CoF[j,k]+ratio
#             min_CoF[j, k] = np.min(CoF_value)*100  # ou np.mean(CoF_value) si tu préfères la moyenne
#     mean_CoF=mean_CoF/len(frequencies)*100
#     Trace_GOF(mean_CoF, "mean", Nom_data, NOMS, dof_vect, len(files), mesh)
#     Trace_GOF(min_CoF, "min", Nom_data, NOMS, dof_vect, len(files), mesh)

def Calcul_GOF(REF, files, Nom_data, dof_vect, NOMS, mesh):
    min_CoF = np.zeros((len(dof_vect), len(files)))
    mean_CoF = np.zeros((len(dof_vect), len(files)))
    use_nrmse=True
    for k in range(len(files)):
        frequencies, COEF = read_RAD_BEM(files[k])
        frequencies, COEF_REF = read_RAD_BEM(REF)
        for j, dof in enumerate(dof_vect):
            diffs_squared = []
            ref_values = []
            gof_local = []
            for i in range(len(frequencies)):
                if frequencies[i]<4:
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
            if gof < 96:
                print(f"mesh={NOMS[k]} - dof={dof} - mean GoF={gof:.2f}% ")
            if gof < 0:
                print(f" !!!!!! mesh={NOMS[k]} - dof={dof} - mean GoF={gof:.2f}%")
                print(f"{len(diffs_squared)} - diffs_squared={diffs_squared}")
    Trace_GOF(mean_CoF, "mean", Nom_data, NOMS, dof_vect, len(files), mesh)
    # Trace_GOF(min_CoF, "min", Nom_data, NOMS, dof_vect, len(files), mesh)

def Calcul_GOF_Fex(REF, files, dof_vect, NOMS, mesh):
    min_CoF=np.zeros((len(dof_vect), len(files)))
    mean_CoF=np.zeros((len(dof_vect), len(files)))
    min_CoF_ph=np.zeros((len(dof_vect), len(files)))
    mean_CoF_ph=np.zeros((len(dof_vect), len(files)))
    for k in range(len(files)) : 
        beta, frequencies, Coef_abs, Coef_ph = read_Fex_BEM(files[k])
        beta, frequencies, Coef_abs_REF, Coef_ph_REF = read_Fex_BEM(REF)

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
                value = Coef_abs[0][i, dof - 1]
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
                value_ph = Coef_ph[0][i, dof - 1]
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
            if gof_abs < 96:
                print(f"dof={dof} - mesh={NOMS[k]} - ABS mean GoF={gof_abs:.2f}%")
            if gof_ph < 96:
                print(f"dof={dof} - mesh={NOMS[k]} - PHASE mean GoF={gof_ph:.2f}%")

    Trace_GOF(mean_CoF, "mean", "abs(Fex)", NOMS, dof_vect, len(files), mesh)
    # Trace_GOF(min_CoF, "min", "|Fex|", NOMS, dof_vect, len(files), mesh)
    Trace_GOF(mean_CoF_ph, "mean", "phase(Fex)", NOMS, dof_vect, len(files), mesh)
    # Trace_GOF(min_CoF_ph, "min", "phase(Fex)", NOMS, dof_vect, len(files), mesh)

def Trace_GOF(GoF, type, data, NOMS, dof_vect, Nfiles, mesh) :
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
        plt.title(f'Goodness of Fit of {data} for DOF {dof} (reference value {NOMS[-1]})')
        plt.grid(axis='y', linestyle='--', alpha=0.7)

        plt.tight_layout()
        # plt.show()
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
    plt.title(f'Computation Time ({unit})')
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()
    # plt.show()
    plt.savefig(f"CPU_{mesh}_{unit}.pdf")
    plt.close()

def Trace_DX(DX, lambda_min, NOMS, mesh) :
    # Création du graphe à bâtons
    data=lambda_min/np.array(DX)
    plt.figure(figsize=(8, 5))
    bars = plt.bar(NOMS, data, color='coral', edgecolor='black')

    # Ajout des valeurs au-dessus des barres
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval + 0.05, f'{yval:.2f}', ha='center', va='bottom')

    # Mise en forme
    plt.xlabel('Mesh')
    plt.ylabel(f'Ratio (-)')
    plt.title(f'Ratio between the min length wave and the mesh size')
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()
    # plt.show()
    plt.savefig(f"Ratio_DX_{mesh}.png", dpi=300)
    plt.close()

def Trace_RAO(files, Nw, Phase, DOF, NOMS, Period, titre, prefix="N_RAO"):
    
    # Tracer l'amplitude
    plt.figure(figsize=(10, 8))
    couleurs = ['b', 'C1', 'g', 'r', 'gold', 'C4', 'c', 'm', 'y', 'k'] 
    styles = ['o', 's']  
    lignes = ['-', '--', ':']
    for k, file in enumerate(files): 
        data = np.loadtxt(file, skiprows=2, max_rows=Nw)
        if Period : 
            freq = 2*np.pi/data[:, 0]
        else :
            freq = data[:, 0]
        plt.plot(freq, data[:, DOF],  '-o', color=couleurs[k], label=f"{NOMS[k]}")
    
    plt.title(f"{titre}")
    if Period : 
        plt.xlabel("Period (s)")
    else :
        plt.xlabel("Frequency (rad/s)")
    plt.ylabel(f"|RAO| - DOF {DOF}")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{prefix}_amplitude_{DOF}.png", dpi=300)
    plt.close()
    
    # Tracer la phase si demandé
    if Phase:
        plt.figure(figsize=(10, 8))
        for k, file in enumerate(files): 
            data = np.loadtxt(file, skiprows=2, max_rows=Nw)
            if Period : 
                freq = 2*np.pi/data[:, 0]
            else :
                freq = data[:, 0]
            plt.plot(freq, data[:, DOF+6], '-o', color=couleurs[k], label=f"{NOMS[k]}")
        plt.title(f"{titre}")
        if Period :
            plt.xlabel("Period (s)")
        else :
            plt.xlabel("Frequency (rad/s)")
        plt.ylabel(f"Phase(RAO) - DOF {DOF}")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{prefix}_phase_{DOF}.png", dpi=300)
        plt.close()

def Calcul_GOF_RAO(REF, files, dof_vect, NOMS, mesh):
    min_CoF=np.zeros((len(dof_vect), len(files)))
    mean_CoF=np.zeros((len(dof_vect), len(files)))
    min_CoF_ph=np.zeros((len(dof_vect), len(files)))
    mean_CoF_ph=np.zeros((len(dof_vect), len(files)))
    for k in range(len(files)) : 
        dataREF = np.loadtxt(REF, skiprows=2, max_rows=20)
        frequencies=dataREF[:, 0]
        data = np.loadtxt(files[k], skiprows=2, max_rows=20)
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
                value = data[i, dof]
                value_ref = dataREF[i, dof]
                diff_sq_abs = np.abs(value - value_ref)**2
                ref_abs = np.abs(value_ref)

                if frequencies[i]<4:
                    abs_diffs_squared.append(diff_sq_abs)
                    abs_refs.append(ref_abs)
                    if not use_nrmse:
                        gof_i_abs = 100 * (1 - (np.sqrt(diff_sq_abs) / ref_abs))
                        gof_local_abs.append(gof_i_abs)

                # --- Phase ---
                value_ph = data[i, dof +6]
                value_ref_ph = dataREF[i, dof +6]
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
            if gof_abs < 95:
                print(f"dof={dof} - mesh={NOMS[k]} - ABS mean GoF={gof_abs:.2f}%")
            if gof_ph < 95:
                print(f"dof={dof} - mesh={NOMS[k]} - PHASE mean GoF={gof_ph:.2f}%")

    Trace_GOF(mean_CoF, "mean", "abs(RAO)", NOMS, dof_vect, len(files), mesh)
    Trace_GOF(mean_CoF_ph, "mean", "phase(RAO)", NOMS, dof_vect, len(files), mesh)

distance="Nb1"
CM = False
CA = False
Fex = False
geo = "PontoonWings"
GRAPHS = True
CoF=False
GoF=True
CPU=True
DX=True
N1_file=""
N2_file=""
if geo=="FB" :
    lambda_min=2*np.pi #k=1, a=3, h=30, w=pi
    folder="CV_FB/BEM_FB_m"
    mesh="FB"    
    file=[f"{folder}0",f"{folder}1", f"{folder}2", f"{folder}3", f"{folder}4",  f"{folder}5"]
    Nmesh=6
    NOMS=["Np=112", "Np=448", "Np=1008", "Np=1792", "Np=2800", "Np=4032"]
    CPU_data=[65.0066681, 467.856537, 2116.52271, 6783.87207, 17088.7832, 36038.5820]
    DX_data=[]  
elif geo=="FB0" :
    folder="CV_FB0/BEM_FB0_m"
    mesh="FB0"    
    file=[f"{folder}1", f"{folder}2", f"{folder}3", f"{folder}4",  f"{folder}5", f"{folder}6", f"{folder}7",  f"{folder}8"]
    Nmesh=8
    NOMS=["Np=88", "Np=162", "Np=304", "Np=430", "Np=648", "Np=826", "Np=1120", "Np=1350"]
    CPU_data=[]
    DX_data=[] 
elif geo=="FB1" :
    folder="CV_FB1/BEM_FB1_m"
    mesh="FB1"    
    file=[f"{folder}1", f"{folder}2", f"{folder}3", f"{folder}4",  f"{folder}5"]
    Nmesh=5
    NOMS=["Np=42", "Np=168", "Np=378", "Np=1050", "Np=1728"]
    CPU_data=[]
    DX_data=[]
elif geo=="Pontoon" :
    lambda_min=13
    folder="Pontoon/BEM_Pontoon_m"
    mesh="Pontoon"    
    file=[f"{folder}1", f"{folder}2", f"{folder}3", f"{folder}4",  f"{folder}5",  f"{folder}6",  f"{folder}7"]
    Nmesh=6
    NOMS=["Np=160", "Np=446", "Np=1200","Np=1674", "Np=2404", "Np=3118"]
    CPU_data=[30.6429844 , 269.450134, 1718.43469, 3649.03540, 8507.04785, 15021.0869]
    DX_data=[10, 5, 3, 2.5, 2, 1.8]  
elif geo=="PontoonWings" :
    lambda_min=13
    folder="PontoonWings/BEM_PontoonWings_m"
    mesh="PontoonWings"    
    file=[f"{folder}1", f"{folder}2", f"{folder}3", f"{folder}4",  f"{folder}5",  f"{folder}6"]
    Nmesh=6
    NOMS=["Np=258", "Np=729", "Np=1739", "Np=2200", "Np=3324", "Np=4017"]
    CPU_data=[186.005844 , 977.346497,  4603.67725 , 6965.42529, 15921.1934, 20957.5566]
    DX_data=[10, 5, 3, 2.5, 2, 1.8]  

# titre=f"Nw=20, Nbeta=11, dof=6, Ndir=1, {distance}"
titre=f"Convergence study on {mesh} mesh"
COMP1=N1_file
COMP2=N2_file
Nw=20
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
RAOfiles=[]
for i in range(Nmesh) :
    CMfiles.append(f"{chemin}{file[i]}/results/CM.dat")
    CAfiles.append(f"{chemin}{file[i]}/results/CA.dat")
    Fefiles.append(f"{chemin}{file[i]}/results/ExcitationForce.tec")
    RAOfiles.append(f"{chemin}{file[i]}/Motion/RAO.dat")


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
        Trace_RAO(RAOfiles, Nw, True, dof, NOMS, True, titre)
if CoF:
    dof=[1, 3, 5]
    Calcul_COF(CAfiles[-1], CAfiles, "Damping", dof, NOMS, mesh)
    Calcul_COF(CMfiles[-1], CMfiles, "Added_Mass", dof, NOMS, mesh)
    Calcul_COF_Fex(Fefiles[-1], Fefiles, dof, NOMS, mesh)

if GoF:
    dof=[1, 3, 5]
    print("------------GoF-----------------")
    print("Added Mass")
    Calcul_GOF(CMfiles[-1], CMfiles, "Added_Mass", dof, NOMS, mesh)
    print("Damping")
    Calcul_GOF(CAfiles[-1], CAfiles, "Damping", dof, NOMS, mesh)
    print("Fex")
    Calcul_GOF_Fex(Fefiles[-1], Fefiles, dof, NOMS, mesh)
    print("RAO")
    Calcul_GOF_RAO(RAOfiles[-1], RAOfiles, dof, NOMS, mesh)

if CPU:
    Trace_CPU(CPU_data, NOMS, mesh, "s")
    CPU_data=np.array(CPU_data)/60
    Trace_CPU(CPU_data, NOMS, mesh, "min")

if DX:
    Trace_DX(DX_data, lambda_min, NOMS, mesh)
    
