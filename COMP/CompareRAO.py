import numpy as np
import re  # Pour utiliser les expressions régulières
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def read_RAO(file):
    beta_values = []  
    matrices = []     
    matrices_abs = []
    matrices_phase = []
    frequencies = [] 

    with open(file, 'r') as f:
        beta = None  # Initialiser une variable pour beta
        matrix = []   
        for line in f:
            line = line.strip() 
            if 'ZONE' in line:
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
                            abs_values = row[1:7]    # colonnes 2 à 7 (indices 1 à 6)
                            phase_values = row[7:13] # colonnes 8 à 13 (indices 7 à 12)
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
                abs_values = row[1:7]  # Valeurs absolues 
                phase_values = row[7:13]  # Phases 
                matrix_abs.append(abs_values)
                matrix_phase.append(phase_values)
        # Ajouter les résultats dans les listes correspondantes
            matrices_abs.append(np.array(matrix_abs))
            matrices_phase.append(np.array(matrix_phase))   
    # return beta_values, frequencies, matrices_abs, matrices_phase
    # Filtrer les matrices en fonction de la première valeur de beta_values
    first_beta_value = beta_values[0] if beta_values else None
    if first_beta_value is not None:
        # Trouver les indices des matrices correspondant à la première valeur de beta
        indices_to_keep = [i for i, beta in enumerate(beta_values) if beta == first_beta_value]
        
        filtered_matrices = [matrices[i] for i in indices_to_keep]
        filtered_matrices_abs = [matrices_abs[i] for i in indices_to_keep]
        filtered_matrices_phase = [matrices_phase[i] for i in indices_to_keep]

        return first_beta_value, frequencies, filtered_matrices_abs, filtered_matrices_phase
    else:
        return None, [], [], []

def CompareResultsRAO(fileBEM, filePIT, param_d, type, Nw, Nb, tol, precision) :
    # beta_BEM, freq_BEM, Coef_BEM_abs, Coef_BEM_ph = read_RAO(fileBEM)
    # beta_PIT, freq_PIT, Coef_PIT_abs, Coef_PIT_ph = read_RAO(filePIT)
    dof=6
    dataBEM = np.loadtxt(fileBEM, skiprows=2, max_rows=Nw)
    freq_BEM = dataBEM[:, 0]
    Coef_BEM_abs = dataBEM[:,1:dof*Nb+1]
    Coef_BEM_ph = dataBEM[:,dof*Nb+1:2*dof*Nb+1]
    dataPIT = np.loadtxt(filePIT, skiprows=2, max_rows=Nw)
    freq_PIT = dataPIT[:, 0]
    Coef_PIT_abs = dataPIT[:,1:dof*Nb+1]
    Coef_PIT_ph = dataPIT[:,dof*Nb+1:2*dof*Nb+1]
        
    # print("BEM", beta_BEM, freq_BEM)
    # print("PIT", beta_PIT, freq_PIT)

    if len(freq_BEM)==len(freq_PIT) :
        if np.shape(Coef_BEM_abs) != np.shape(Coef_PIT_abs):
            print(np.shape(Coef_BEM_abs), np.shape(Coef_PIT_abs))
            raise ValueError("tailles différentes")
        if np.shape(Coef_BEM_ph) != np.shape(Coef_PIT_ph):
            print(np.shape(Coef_BEM_ph), np.shape(Coef_PIT_ph))
            raise ValueError("tailles différentes")
        else :
            errorL2 = np.sqrt(np.sum((np.array(Coef_BEM_abs) - np.array(Coef_PIT_abs))**2))
            # errorINF = np.max(np.abs(np.array(Coef_BEM_abs) - np.array(Coef_PIT)))
            erreurs = np.abs(np.array(Coef_BEM_abs) - np.array(Coef_PIT_abs))
            errorINF = np.max(erreurs)
           
            normL2_BEM = np.sqrt(np.sum(np.array(Coef_BEM_abs)**2))
            relative_errorL2_abs = errorL2 / normL2_BEM if normL2_BEM != 0 else np.nan

            max_BEM = np.max(np.abs(np.array(Coef_BEM_abs)))
            relative_errorINF_abs = errorINF / max_BEM if max_BEM != 0 else np.nan

            # PHASE 
            errorL2_p = np.sqrt(np.sum((np.array(Coef_BEM_ph) - np.array(Coef_PIT_ph))**2))
            # errorINF_p = np.max(np.abs(np.array(Coef_BEM_ph) - np.array(Phase_PIT)))
            erreurs_p = np.abs(np.array(Coef_BEM_ph) - np.array(Coef_PIT_ph))
            errorINF_p = np.max(erreurs_p)
           
            normL2_BEM_p = np.sqrt(np.sum(np.array(Coef_BEM_ph)**2))
            relative_errorL2_phase = errorL2_p / normL2_BEM_p if normL2_BEM_p != 0 else np.nan

            max_BEM_p = np.max(np.abs(np.array(Coef_BEM_ph)))
            relative_errorINF_phase = errorINF_p / max_BEM_p if max_BEM_p != 0 else np.nan

            # Calcul erreur pour chaque dof 
            All_errors=True
            if All_errors:
                dof=dof*Nb
                n_err_abs=0
                n_err_ph=0
                for i in range(dof):
                    B_abs = Coef_BEM_abs[:,i]
                    P_abs = Coef_PIT_abs[:,i]
                    B_ph = Coef_BEM_ph[:,i]
                    P_ph = Coef_PIT_ph[:,i]
                    # Remplacer les valeurs < precision par 0
                    # B_abs = [0 if np.abs(x) < precision else x for x in B_abs]
                    # P_abs = [0 if np.abs(x) < precision else x for x in P_abs]
                    # B_ph = [0 if np.abs(b) < precision else ph for b, ph in zip(B_abs, B_ph)]
                    # P_ph = [0 if np.abs(b) < precision else ph for b, ph in zip(P_abs, P_ph)]
                    if np.max(np.abs(np.array(B_abs))) != 0 and np.max(np.abs(np.array(P_abs))) != 0 : 
                        err_abs=np.max(np.abs(np.array(B_abs)-np.array(P_abs)))/np.max(np.abs(np.array(B_abs))) 
                    else : 
                        err_abs=0
                    if np.max(np.abs(np.array(B_ph))) != 0 and np.max(np.abs(np.array(P_ph))) != 0 : 
                        err_ph=np.max(np.abs(np.array(B_ph)-np.array(P_ph)))/np.max(np.abs(np.array(B_ph))) 
                    else : 
                        err_ph=0
                    
                    with open(f"EcartsRelatifs_RAO_abs_{type}.dat", "a") as f:
                        if i==0:
                            f.write("\n")
                            f.write(f"{param_d}     ")
                        f.write(f"{err_abs:16.8f}     ")
                    if err_abs<=tol:
                        n_err_abs=n_err_abs+1
                        
                    with open(f"EcartsRelatifs_RAO_ph_{type}.dat", "a") as f:
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



def plot_multiple_RAO(file_refs, noms, distances, Nb, Nw, titre, Phase, max_X):
    """
    file_refs : liste de chemins vers les fichiers REF (RAO)
    noms      : liste de noms pour les légendes et fichiers (ex: ['REF', 'Simulation1', 'Simulation2'])
    distances : liste des identifiants ou distances pour les titres/sauvegardes
    """
    couleurs = ['g', 'r', 'c', 'm', 'y', 'k']  # palette de couleurs (réutilisée si plus de 7 courbes)

    # Initialiser les figures pour les 12 courbes
    dof=6
    nb_courbe=dof
    if Phase : 
        nb_courbe=2*nb_courbe
    figs = [plt.figure(i) for i in range(1, nb_courbe+1)]
    axes = [fig.add_subplot(1,1,1) for fig in figs]
    
    for idx, fileREF in enumerate(file_refs):
        data = np.loadtxt(fileREF, skiprows=2, max_rows=Nw)
        freq = data[:, 0]
        

        # Tracer les RAO amplitudes (1 à 6)
        for i in range(dof):
            if noms[idx]=="REF":
                axes[i].plot(freq, data[:, i+1], f'b-o',markersize=2, label=f'{noms[idx]} d=0')
            else : 
                for j in range(Nb) : 
                    color = couleurs[j % len(couleurs)]
                    axes[i].plot(freq, data[:, (j*6)+i+1], f'{color}-o',markersize=2, label=f'{noms[idx]} - body{j+1} for d={distances[idx]}')

        # Tracer les RAO phases (7 à 12)
        if Phase : 
            for i in range(dof, 2*dof):
                if noms[idx]=="REF":
                    axes[i].plot(freq, data[:, i+1], f'b-o',markersize=2, label=f'{noms[idx]} d=0')
                else : 
                    for k in range(Nb) : 
                        color = couleurs[k % len(couleurs)]
                        j=i+6
                        axes[i].plot(freq, data[:, (k*6)+j+1], f'{color}-o',markersize=2, label=f'{noms[idx]} - body{k+1} d={distances[idx]}')
                    # axes[i].plot(freq, data[:, 6+j+1], f'{colorbis}-o',markersize=4, label=f'body2 d={distances[idx]}')

    # Mise en forme des figures et sauvegarde
    for i in range(dof):
        axes[i].set_xlim(0, max_X)
        # axes[i].set_ylim(-0.1, 0.1)
        axes[i].set_title(f"RAO |Amplitude| DOF {i+1} \n {titre}")
        axes[i].set_xlabel("Fréquence (rad/s)")
        axes[i].set_ylabel(f"|RAO| DOF {i+1}")
        axes[i].grid(True)
        axes[i].legend()
        figs[i].tight_layout()
        figs[i].savefig(f"N3_RAO_ref_{noms[1]}_amplitude_DOF{i+1}_d{distances}.png", dpi=300)

    if Phase : 
        for i in range(dof, 2*dof):
            axes[i].set_xlim(0, max_X)
            axes[i].set_title(f"RAO Phase DOF {i-5}")
            axes[i].set_xlabel("Fréquence (rad/s)")
            axes[i].set_ylabel(f"Phase(RAO) DOF {i-5} \n {titre}")
            axes[i].grid(True)
            axes[i].legend()
            figs[i].tight_layout()
            figs[i].savefig(f"N3_RAO_ref_{noms[1]}_phase_DOF{i-5}_d{distances}.png", dpi=300)

    # Fermer les figures (optionnel si on ne veut pas les afficher)
    for fig in figs:
        plt.close(fig)




def plot_RAO_multiple_d(fileREF, filesBEM, filesPIT, vect_DOF, distances, Nb, Nw, titre, Phase, max_X):
    """
    files : liste de chemins vers les fichiers REF (RAO)
    distances : liste des identifiants ou distances pour les titres/sauvegardes
    """
    couleurs = ['b', 'g', 'r', 'c', 'm', 'y', 'k']  # palette de couleurs (réutilisée si plus de 7 courbes)

    # Initialiser les figures pour les 12 courbes
    nb_figures = len(vect_DOF) * Nb
    if Phase : 
        nb_figures=2*nb_figures
    figs = [plt.figure(i) for i in range(1, nb_figures+1)]
    axes = [fig.add_subplot(1,1,1) for fig in figs]
    nb_dof=6

    dataREF = np.loadtxt(fileREF, skiprows=2, max_rows=Nw)
    freqREF = dataREF[:, 0]
    indice_graph=0
    for i, dof in enumerate(vect_DOF):
        for j in range(Nb) : 
            axes[indice_graph].plot(freqREF, dataREF[:, dof], linestyle='--', color='black', label=f'REF d=0')

            for idx, fileRAO in enumerate(filesBEM):
                data = np.loadtxt(fileRAO, skiprows=2, max_rows=Nw)
                freq = data[:, 0]
                color = couleurs[idx % len(couleurs)]
                axes[indice_graph].plot(freq, data[:, (j*6)+dof], f'{color}-o',markersize=2, label=f'BEM - body{j+1} for d={distances[idx]}')
            for idx, filePIT in enumerate(filesPIT):
                dataPIT = np.loadtxt(filePIT, skiprows=2, max_rows=Nw)
                freqIT = dataPIT[:, 0]
                color = couleurs[idx+2 % len(couleurs)]
                axes[indice_graph].plot(freqIT, dataPIT[:, (j*6)+dof], f'{color}--+',markersize=4, label=f'PIT - body{j+1} for d={distances[idx]}')
            axes[indice_graph].set_xlim(0, max_X)
            axes[indice_graph].set_title(f"RAO |Amplitude| - body{j+1} - DOF {dof+j*6} \n {titre}")
            axes[indice_graph].set_xlabel("Fréquence (rad/s)")
            axes[indice_graph].set_ylabel(f"|RAO| DOF {dof}")
            axes[indice_graph].grid(True)
            axes[indice_graph].legend()
            figs[indice_graph].tight_layout()
            figs[indice_graph].savefig(f"N3_RAO_amplitude_body{j+1}_DOF{dof+j*6}_d{distances}.png", dpi=300)
            indice_graph=indice_graph+1

        # Tracer les RAO phases (7 à 12)
    if Phase : 
        offset = nb_dof*Nb
        for i, dof in enumerate(vect_DOF):
            for j in range(Nb) : 
                k=dof-1+offset + 6*j
                axes[indice_graph].plot(freqREF, dataREF[:, dof+6], linestyle='--', color='black', label=f'REF d=0')

                for idx, fileRAO in enumerate(filesBEM):
                    data = np.loadtxt(fileRAO, skiprows=2, max_rows=Nw)
                    freq = data[:, 0]
                    color = couleurs[idx % len(couleurs)]
                    axes[indice_graph].plot(freq, data[:, k+1], f'{color}-o',markersize=2, label=f'BEM - body{j+1} d={distances[idx]}')
                for idx, filePIT in enumerate(filesPIT):
                    dataPIT = np.loadtxt(filePIT, skiprows=2, max_rows=Nw)
                    freqIT = dataPIT[:, 0]
                    color = couleurs[idx+2 % len(couleurs)]
                    axes[indice_graph].plot(freqIT, dataPIT[:, k+1], f'{color}-+',markersize=4, label=f'PIT - body{j+1} for d={distances[idx]}')
            
                # print("indice",i, dof, j, k, idx, indice_graph)
                axes[indice_graph].set_xlim(0, max_X)
                axes[indice_graph].set_title(f"RAO Phase - body{j+1} - DOF {dof+j*6} \n {titre}")
                axes[indice_graph].set_xlabel("Fréquence (rad/s)")
                axes[indice_graph].set_ylabel(f"Phase(RAO) DOF {dof} ")
                axes[indice_graph].grid(True)
                axes[indice_graph].legend()
                figs[indice_graph].tight_layout()
                figs[indice_graph].savefig(f"N3_RAO_phase_body{j+1}_DOF{dof+j*6}_d{distances}.png", dpi=300)
                indice_graph=indice_graph+1
    # Fermer les figures (optionnel si on ne veut pas les afficher)
    for fig in figs:
        plt.close(fig)


chemin="/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/"
N3 = False
PLOT_multi_d=True
PLOT_REF=False
Plot_phase=False
Nb=2
config="X"
Nw=60
max_X=2*np.pi
test="De_"
beta="beta_0.0_"
PIT_type=f"Nb{Nb}_{config}"
titre=f"Barge - {PIT_type}, Nw={Nw}, Nbeta=11, dof=6, Ndir=1"
# BEM= f"{chemin}PFE_Nemoh/MyTestCases/BEM_{PIT_type}"
BEM= f"{chemin}PFE_Nemoh/MyTestCases/BEM_barge_dof6_R6.36_w2pi"
PIT= f"{chemin}PIT3_E/wec_inputs/PIT3_barge_w2pi"
RAO_REF = f"{chemin}COMP/Motion_1body/RAO.dat"

if PLOT_REF : 
    d_graphs= [[1], [16], [64]]
    
    
    print("\n ------------ PLOT with REF -----------------")

    for j, vect_d in enumerate(d_graphs): 
        files = [RAO_REF]
        filesPIT = [RAO_REF]
        names = ['REF']
        namesPIT = ['REF']
        distances = [0]

        param_d=vect_d
        print("\n ------------Distances=", param_d, "-----------------")
        for i, d in enumerate(param_d) : 
            # RAO_BEM = f"{BEM}/BEM_barge_Nb{Nb}{config}_d{d}/Motion/RAO.dat" 
            RAO_BEM = f"{BEM}/BEM_barge_d{d}/Motion/RAO.dat" 
            RAO_PIT= f"{PIT}/Motion/Global_{test}RAO_{PIT_type}_d{d}.00.dat"
            files.append(RAO_BEM) 
            filesPIT.append(RAO_PIT) 
            distances.append(d)
            names.append('BEM') 
            namesPIT.append('PIT') 
        plot_multiple_RAO(files, names, distances, Nb, Nw, titre, Plot_phase, max_X)
        plot_multiple_RAO(filesPIT, namesPIT, distances, Nb, Nw, titre, Plot_phase, max_X)

if PLOT_multi_d : 
    print("\n ------------ PLOT without REF -----------------")

    dist_graphs= [[1], [16], [64]]
    DOF=[1, 3, 5]
    for j, vect_d in enumerate(dist_graphs): 
        BEM_files=[]
        PIT_files=[]
        print("\n ---DOF=", DOF, "Distances=", vect_d, "--------")

        for i, d in enumerate(vect_d) : 
            # RAO_BEM = f"{BEM}/BEM_barge_Nb{Nb}{config}_d{d}/Motion/RAO.dat"
            RAO_BEM = f"{BEM}/BEM_barge_d{d}/Motion/RAO.dat"
            RAO_PIT= f"{PIT}/Motion/Global_{test}RAO_{PIT_type}_{beta}d{d}.00.dat"
            BEM_files.append(RAO_BEM) 
            PIT_files.append(RAO_PIT) 
        plot_RAO_multiple_d(RAO_REF, BEM_files, PIT_files, DOF, vect_d, Nb, Nw, titre, Plot_phase, max_X)





if N3 : 
    param_d=8
    limite=70
    print("Comparaison sur RAO")
    while param_d<=limite : 
        print("\n ------------Distance=", param_d, "-----------------")
        distance=f"d{param_d}"
        # BEM_file=f"BEM_cylinder_{PIT_type}/BEM_cylinder_Nb{Nb}_{config}_{distance}"

        titre=f"{titre}, {distance}"

        # RAO_BEM = f"{BEM}/BEM_barge_Nb{Nb}{config}_d{param_d}/Motion/RAO.dat"
        RAO_BEM = f"{BEM}/BEM_barge_d{param_d}/Motion/RAO.dat"
        RAO_PIT= f"{PIT}/Motion/Global_{test}RAO_{PIT_type}_beta_0.0_d{param_d}.00.dat"

        with open(f"EcartsRelatifs_RAO_{PIT_type}.dat", "a") as f:
            f.write(f"{param_d}     ")
            
            L_inf_error_abs, L2_error_abs, L_inf_error_ph, L2_error_ph = CompareResultsRAO(RAO_BEM, RAO_PIT, param_d, PIT_type, Nw, Nb, 0.1, 10)
            # print("Ecart relatif entre les données pour |Fex| :")
            print(L_inf_error_abs, L2_error_abs)
            f.write(f"{L_inf_error_abs:16.8f}     {L2_error_abs:16.8f}      ")
            # print("Ecart relatif entre les données pour phase(Fex) :")
            # print(L_inf_error_ph, L2_error_ph)
            f.write(f"{L_inf_error_ph:16.8f}      {L2_error_ph:16.8f}\n")
        
        param_d=param_d*2
