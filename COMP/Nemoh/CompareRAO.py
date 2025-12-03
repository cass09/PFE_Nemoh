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
    Coef_PIT_abs[:,3:6]=Coef_PIT_abs[:,3:6]*180/np.pi
    Coef_PIT_ph[:,:]=Coef_PIT_ph[:,:]*180/np.pi
    Coef_PIT_abs[:,9:12]=Coef_PIT_abs[:,9:12]*180/np.pi
    Coef_PIT_abs[:,16:18]=Coef_PIT_abs[:,16:18]*180/np.pi
        
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
                    axes[i].plot(freq, data[:, (j*6)+i+1], f'{color}-o',markersize=2, label=f'{noms[idx]} - body{j+1} for d={distances[idx]}m')

        # Tracer les RAO phases (7 à 12)
        if Phase : 
            for i in range(dof, 2*dof):
                if noms[idx]=="REF":
                    axes[i].plot(freq, data[:, i+1], f'b-o',markersize=2, label=f'{noms[idx]} d=0')
                else : 
                    for k in range(Nb) : 
                        color = couleurs[k % len(couleurs)]
                        j=i+6
                        axes[i].plot(freq, data[:, (k*6)+j+1], f'{color}-o',markersize=2, label=f'{noms[idx]} - body{k+1} d={distances[idx]}m')
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




def plot_RAO_multiple_d(fileREF, filesBEM, filesPIT, filesPIT_REF, vect_DOF, distances, Nb, Nw, titre, Phase, max_X, test, MultiBody):
    """
    files : liste de chemins vers les fichiers REF (RAO)
    distances : liste des identifiants ou distances pour les titres/sauvegardes
    """
    couleurs = ['b', 'r', 'g', 'gold', 'm', 'darkcyan', 'orange', 'y', 'k']  # palette de couleurs (réutilisée si plus de 7 courbes)

    # Initialiser les figures pour les 12 courbes
    nb_figures = len(vect_DOF) * Nb
    if Phase : 
        nb_figures=2*nb_figures
    figs = [plt.figure(i) for i in range(1, nb_figures+1)]
    axes = [fig.add_subplot(1,1,1) for fig in figs]
    nb_dof=6

    # dataREF = np.loadtxt(fileREF, skiprows=2, max_rows=Nw)
    # freqREF = dataREF[:, 0]
    indice_graph=0
    for i, dof in enumerate(vect_DOF):
        # axes[indice_graph].plot(freqREF, dataREF[:, dof], linestyle='--', color='black', label=f'REF')
        for j in range(Nb) : 
            if MultiBody :
                legende=f" - DOF {dof+j*6}"
            else :
                legende=""
            for idx, fileRAO in enumerate(filesBEM):
                data = np.loadtxt(fileRAO, skiprows=2, max_rows=Nw)
                freq = data[:, 0]
                color = couleurs[idx+j % len(couleurs)]
                axes[indice_graph].plot(freq, data[:, (j*6)+dof], linestyle='-', color=color,marker='o', markersize=3, label=f'BEM {legende}')
            for idx, filePIT in enumerate(filesPIT):
                dataPIT = np.loadtxt(filePIT, skiprows=2, max_rows=Nw)
                freqIT = dataPIT[:, 0]
                # color = couleurs[idx+Nb+j % len(couleurs)]
                RAOvalue=dataPIT[:, (j*6)+dof]
                if dof>3 :
                    RAOvalue=RAOvalue*180/np.pi
                axes[indice_graph].plot(freqIT, RAOvalue, linestyle='--', color=color,marker='x',markersize=6, label=f'ITM {legende}')
                if filesPIT_REF[idx] != filePIT :
                    dataPIT_REF = np.loadtxt(filesPIT_REF[idx], skiprows=2, max_rows=Nw)
                    freqIT = dataPIT_REF[:, 0]
                    color2 = couleurs[idx+Nb+j+2 % len(couleurs)]
                    RAOvalueREF=dataPIT_REF[:, (j*6)+dof]
                    if dof>3 :
                        RAOvalueREF=RAOvalueREF*180/np.pi
                    axes[indice_graph].plot(freqIT, RAOvalueREF, linestyle=':', color=color,marker='<',markersize=4, label=f'ITM with Capytaine {legende}')
                
            # axes[indice_graph].set_xlim(0, max_X)
            if MultiBody :
                num_dof=f" DOF {dof} and {dof+j*6} "
            else :
                num_dof=f" DOF {dof+j*6}"
            # axes[indice_graph].set_title(f"RAO |Amplitude| - {num_dof} - d={distances[0]}m - {test}\n {titre}")
            axes[indice_graph].set_title(f"{titre}")
            axes[indice_graph].set_xlabel("Frequency (rad/s)")
            axes[indice_graph].set_ylabel(f"|RAO|")
            axes[indice_graph].grid(True)
            axes[indice_graph].legend()
            figs[indice_graph].tight_layout()
            if not MultiBody :
                figs[indice_graph].savefig(f"N3_{test}RAO_amplitude_DOF{dof+j*6}_d{distances}.png", dpi=300)
                indice_graph=indice_graph+1
        if MultiBody :
            figs[indice_graph].savefig(f"N3_{test}RAO_amplitude_DOFs{dof}_{dof+j*6}_d{distances[0]}.pdf")
            indice_graph=indice_graph+1

        # Tracer les RAO phases (7 à 12)
    if Phase : 
        offset = nb_dof*Nb
        for i, dof in enumerate(vect_DOF):
            # axes[indice_graph].plot(freqREF, dataREF[:, dof+6], linestyle='--', color='black', label=f'REF')
            for j in range(Nb) : 
                if MultiBody :
                    legende=f" - DOF {dof+j*6}"
                else :
                    legende=""
                k=dof-1+offset + 6*j
                for idx, fileRAO in enumerate(filesBEM):
                    data = np.loadtxt(fileRAO, skiprows=2, max_rows=Nw)
                    freq = data[:, 0]
                    color = couleurs[idx+j % len(couleurs)]
                    axes[indice_graph].plot(freq, data[:, k+1], linestyle='-', color=color,marker='o',markersize=3, label=f'BEM {legende}')
                for idx, filePIT in enumerate(filesPIT):
                    dataPIT = np.loadtxt(filePIT, skiprows=2, max_rows=Nw)
                    freqIT = dataPIT[:, 0]
                    # color = couleurs[idx+Nb+j % len(couleurs)]
                    axes[indice_graph].plot(freqIT, dataPIT[:, k+1]*180/np.pi,linestyle='--', color=color, marker='x', markersize=6, label=f'ITM {legende}')
                    if filesPIT_REF[idx] != filePIT :
                        dataPIT_REF = np.loadtxt(filesPIT_REF[idx], skiprows=2, max_rows=Nw)
                        freqIT = dataPIT_REF[:, 0]
                        color2 = couleurs[idx+Nb+j+2 % len(couleurs)]
                        axes[indice_graph].plot(freqIT, dataPIT_REF[:, k+1]*180/np.pi, linestyle=':', color=color,marker='<',markersize=4, label=f'ITM with Capytaine {legende}')
                
                # print("indice",i, dof, j, k, idx, indice_graph)
                # axes[indice_graph].set_xlim(0, max_X)
                if MultiBody :
                    num_dof=f" DOF {dof} and {dof+j*6} "
                else :
                    num_dof=f" DOF {dof+j*6}"
                axes[indice_graph].set_title(f"{titre}")
                axes[indice_graph].set_xlabel("Frequency (rad/s)")
                axes[indice_graph].set_ylabel(f"Phase(RAO)")
                axes[indice_graph].grid(True)
                axes[indice_graph].legend()
                figs[indice_graph].tight_layout()
                if not MultiBody :
                    figs[indice_graph].savefig(f"N3_{test}RAO_phase_DOF{dof+j*6}_d{distances}.png", dpi=300)
                    indice_graph=indice_graph+1
            if MultiBody :
                figs[indice_graph].savefig(f"N3_{test}RAO_phase_DOFs{dof}_{dof+j*6}_d{distances[0]}.pdf")
                indice_graph=indice_graph+1
    # Fermer les figures (optionnel si on ne veut pas les afficher)
    for fig in figs:
        plt.close(fig)

def plot_RAO_multiple_E(fileREF, fileBEM, filesPIT, filePIT_REF, vect_DOF, distances, Nb, Nw, titre, Phase, max_X, Ne, MultiBody):
    """
    files : liste de chemins vers les fichiers REF (RAO)
    distances : liste des identifiants ou distances pour les titres/sauvegardes
    """
    couleurs = ['b', 'g', 'gold', 'm', 'darkcyan', 'orange', 'y', 'k']  # palette de couleurs (réutilisée si plus de 7 courbes)

    # Initialiser les figures pour les 12 courbes
    nb_figures = len(vect_DOF) * Nb
    if Phase : 
        nb_figures=2*nb_figures
    figs = [plt.figure(i) for i in range(1, nb_figures+1)]
    axes = [fig.add_subplot(1,1,1) for fig in figs]
    nb_dof=6

    dataREF = np.loadtxt(fileREF, skiprows=2, max_rows=Nw)
    freqREF = dataREF[:, 0]
    dataREF_PIT = np.loadtxt(filePIT_REF, skiprows=2, max_rows=Nw)
    freqREF_PIT = dataREF_PIT[:, 0]
    # for idx, fileRAO in enumerate(filesBEM):
    dataBEM = np.loadtxt(fileBEM, skiprows=2, max_rows=Nw)
    freqBEM = dataBEM[:, 0]    
    indice_graph=0
    for i, dof in enumerate(vect_DOF):
        axes[indice_graph].plot(freqREF, dataREF[:, dof], linestyle='--', color='black', label=f'REF')
        for j in range(Nb) : 
            if MultiBody :
                legende=f" - body{j+1}"
            else :
                legende=""
            axes[indice_graph].plot(freqBEM, dataBEM[:, (j*6)+dof], linestyle='--', color='b',marker='o', markersize=4, label=f'BEM {legende}')
            RAO_PIT_REF=dataREF_PIT[:, (j*6)+dof]
            if dof>3 :
                RAO_PIT_REF=RAO_PIT_REF*180/np.pi
            axes[indice_graph].plot(freqREF_PIT, RAO_PIT_REF, linestyle='-', color='red', linewidth=0.8, label=f'PIT Ne=0')
            for idx, filePIT in enumerate(filesPIT):
                # print(idx, "Ne=", Ne[idx])
                dataPIT = np.loadtxt(filePIT, skiprows=2, max_rows=Nw)
                freqIT = dataPIT[:, 0]
                color = couleurs[idx+Nb+j % len(couleurs)]
                RAOvalue=dataPIT[:, (j*6)+dof]
                if dof>3 :
                    RAOvalue=RAOvalue*180/np.pi
                axes[indice_graph].plot(freqIT, RAOvalue, linestyle=':', color=color,marker='+',markersize=8, label=f'PIT - Ne={Ne[idx]}{legende}')
            axes[indice_graph].set_xlim(0, max_X)
            if MultiBody :
                num_dof=f" DOF {dof} and {dof+j*6} "
            else :
                num_dof=f" DOF {dof+j*6}"
            axes[indice_graph].set_title(f"RAO |Amplitude| - {num_dof} - d={distances}m - Ne={Ne}\n {titre}")
            axes[indice_graph].set_xlabel("Frequency (rad/s)")
            axes[indice_graph].set_ylabel(f"|RAO| DOF {dof}")
            axes[indice_graph].grid(True)
            axes[indice_graph].legend()
            figs[indice_graph].tight_layout()
            if not MultiBody :
                figs[indice_graph].savefig(f"N3_E_RAO_amplitude_DOF{dof+j*6}_d{distances}.png", dpi=300)
                indice_graph=indice_graph+1
        if MultiBody :
            figs[indice_graph].savefig(f"N3_E_RAO_amplitude_DOFs{dof}_{dof+j*6}_d{distances}.png", dpi=300)
            indice_graph=indice_graph+1

        # Tracer les RAO phases (7 à 12)
    if Phase : 
        offset = nb_dof*Nb
        for i, dof in enumerate(vect_DOF):
            axes[indice_graph].plot(freqREF, dataREF[:, dof+6], linestyle='--', color='black', label=f'REF')
            for j in range(Nb) : 
                if MultiBody :
                    legende=f" - body{j+1}"
                else :
                    legende=""
                k=dof-1+offset + 6*j
                axes[indice_graph].plot(freqBEM, dataBEM[:, k+1], linestyle='--', color='b',marker='o',markersize=2, label=f'BEM {legende}')
                axes[indice_graph].plot(freqREF_PIT, dataREF_PIT[:, k+1]*180/np.pi, linestyle='-', color='red', linewidth=0.8, label=f'PIT only P waves')
                for idx, filePIT in enumerate(filesPIT):
                    dataPIT = np.loadtxt(filePIT, skiprows=2, max_rows=Nw)
                    freqIT = dataPIT[:, 0]
                    color = couleurs[idx+Nb+j % len(couleurs)]
                    axes[indice_graph].plot(freqIT, dataPIT[:, k+1]*180/np.pi,linestyle=':', color=color,marker='+',markersize=6, label=f'PIT - Ne={Ne[idx]}{legende}')
                # print("indice",i, dof, j, k, idx, indice_graph)
                axes[indice_graph].set_xlim(0, max_X)
                if MultiBody :
                    num_dof=f" DOF {dof} and {dof+j*6} "
                else :
                    num_dof=f" DOF {dof+j*6}"
                axes[indice_graph].set_title(f"RAO Phase - {num_dof} - d={distances}m - Ne={Ne}\n {titre}")
                axes[indice_graph].set_xlabel("Fréquence (rad/s)")
                axes[indice_graph].set_ylabel(f"Phase(RAO) DOF {dof} ")
                axes[indice_graph].grid(True)
                axes[indice_graph].legend()
                figs[indice_graph].tight_layout()
                if not MultiBody :
                    figs[indice_graph].savefig(f"N3_E_RAO_phase_DOF{dof+j*6}_d{distances}.png", dpi=300)
                    indice_graph=indice_graph+1
        if MultiBody :
            figs[indice_graph].savefig(f"N3_E_RAO_phase_DOFs{dof}_{dof+j*6}_d{distances}.png", dpi=300)
            indice_graph=indice_graph+1
    # Fermer les figures (optionnel si on ne veut pas les afficher)
    for fig in figs:
        plt.close(fig)


def Calcul_GOF(BEMfiles, PITfiles, dof_vect, NOMS, mesh):
    Nw=20
    min_CoF=np.zeros((len(dof_vect), len(BEMfiles)))
    mean_CoF=np.zeros((len(dof_vect), len(BEMfiles)))
    min_CoF_ph=np.zeros((len(dof_vect), len(BEMfiles)))
    mean_CoF_ph=np.zeros((len(dof_vect), len(BEMfiles)))
    for k in range(len(BEMfiles)) : 
        use_nrmse = True  # ⬅️ Change à False pour utiliser la méthode locale
        data_BEM=np.loadtxt(BEMfiles[k], skiprows=2, max_rows=Nw)
        data_PIT=np.loadtxt(PITfiles[k], skiprows=2, max_rows=Nw)
        frequencies = data_PIT[:, 0]
                
        for j, dof in enumerate(dof_vect):
            abs_diffs_squared = []
            abs_refs = []
            ph_diffs_squared = []
            ph_refs = []

            gof_local_abs = []
            gof_local_ph = []
            Coef_abs_REF=data_BEM[:, dof]
            Coef_ph_REF=data_BEM[:, 12+dof]
            Coef_abs=data_PIT[:, dof]
            Coef_ph=data_PIT[:, 12+dof]*180/np.pi
            if dof in [4, 5, 6, 10, 11, 12] :
                Coef_abs=Coef_abs*180/np.pi

            for i in range(len(frequencies)):
                # --- Amplitude ---
                value = Coef_abs[i]
                value_ref = Coef_abs_REF[i]
                diff_sq_abs = np.abs(value - value_ref)**2
                ref_abs = np.abs(value_ref)

                if frequencies[i]<4:
                    abs_diffs_squared.append(diff_sq_abs)
                    abs_refs.append(ref_abs)
                    if not use_nrmse:
                        gof_i_abs = 100 * (1 - (np.sqrt(diff_sq_abs) / ref_abs))
                        gof_local_abs.append(gof_i_abs)

                # --- Phase ---
                value_ph = Coef_ph[i]
                value_ref_ph = Coef_ph_REF[i]
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
            data="RAOph"
            plt.title(f'Goodness of Fit - Two Cylinders - Response Amplitude Operator Phase ph(RAO)_{dof}')
        else : 
            data="RAOabs"
            plt.title(f'Goodness of Fit - Two Cylinders - Response Amplitude Operator Magnitude |RAO|_{dof}')
        plt.grid(axis='y', linestyle='--', alpha=0.7)

        plt.tight_layout()
        # plt.show()
        plt.savefig(f"GOF_{mesh}_{data}_dof{dof}.pdf")
        plt.close()

chemin="/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/"
N3 = False
PLOT_multi_d=True
PLOT_multi_E=False
GOF=False
# PLOT_REF=False
Plot_phase=True
MultiBody=True
Nb=4
config="C"
Nw=20
max_X=3.2
test=""
beta=""
mesh="CylR3"
PIT_type=f"Nb{Nb}_{config}_"
# titre=f"{mesh} - {PIT_type}, Nw={Nw}, Nbeta=11, dof=6, Ndir=1"
titreG = "Four Cylinders"
BEM= f"{chemin}/PFE_Nemoh/MyTestCases/{mesh}/BEM_{mesh}_Nb{Nb}_{config}_"
PIT= f"{chemin}PIT3_E/wec_inputs/PIT3_{mesh}_h12_bis/"
CAP= f"{chemin}capytaine/MyTestCases/CAP_{mesh}_h12/"
RAO_REF = f"{chemin}COMP/Motion_1body/RAO_{mesh}.dat"

if GOF : 
    BEMfiles=[]
    PITfiles=[]
    NOMS=[]
    distances=[2, 4, 8, 16]
    for i, d in enumerate(distances) :
        RAO_BEM = f"{BEM}da{d}/Motion/RAO.dat"
        RAO_PIT= f"{PIT}Motion/Global_{test}RAO_{PIT_type}_{beta}da{d}.00.dat"
        BEMfiles.append(RAO_BEM)
        PITfiles.append(RAO_PIT)
        NOMS.append(f"d/a={d}")
    dof=[1, 3, 5, 7, 9, 11]
    print(f"------------GoF -> {NOMS}")
    Calcul_GOF(BEMfiles, PITfiles, dof, NOMS, mesh)
    
if PLOT_multi_d : 
    print("\n ------------ PLOT -----------------")

    # dist_graphs= [[1], [8], [16], [64]]
    dist_graphs= [[6.25]]
    DOF=[1, 3, 5]
    for j, vect_d in enumerate(dist_graphs): 
        BEM_files=[]
        PIT_files=[]
        PIT_files_CAP=[]
        print("\n ---DOF=", DOF, "Distances=", vect_d, "--------")
	
        for i, d in enumerate(vect_d) : 
            titre=f"{titreG} - d=a+Rc=6.25m"
            RAO_BEM = f"{BEM}d6.25/Motion/RAO.dat"
            # RAO_BEM = f"{BEM}/BEM_{mesh}_d{d}/Motion/RAO.dat"
            # PIT_type=""
            RAO_PIT= f"{PIT}Motion/Global_{test}RAO_{PIT_type}{beta}d{d}.dat"
            RAO_PIT_CAP=RAO_PIT
            # RAO_PIT_CAP= f"{CAP}Motion/CapytaineIT_RAO_{PIT_type}_{beta}da{d}.00.dat"
            BEM_files.append(RAO_BEM) 
            PIT_files.append(RAO_PIT) 
            PIT_files_CAP.append(RAO_PIT_CAP) 
        plot_RAO_multiple_d(RAO_REF, BEM_files, PIT_files, PIT_files_CAP, DOF, vect_d, Nb, Nw, titre, Plot_phase, max_X, test, MultiBody)

if PLOT_multi_E : 
    print("\n ------------ PLOT -----------------")

    # dist_graphs= [[1], [8], [16], [64]]
    vect_d=[1, 16]
    Ne=[1, 6]
    DOF=[1, 3, 5]
    # for j, vect_d in enumerate(dist_graphs): 
    print("\n ---DOF=", DOF, "- Ne=", Ne)

    for c, distance in enumerate(vect_d):
        print("\n  Distance=", distance)
        PIT_files=[]
        RAO_BEM = f"{BEM}/BEM_{mesh}_Nb{Nb}_{config}_d{distance}/Motion/RAO.dat"
        # RAO_BEM = f"{BEM}/BEM_{mesh}_d{d}/Motion/RAO.dat"
        RAO_PIT_REF= f"{PIT}/Motion/Global_RAO_{PIT_type}_{beta}d{distance}.00.dat"
        for i, E in enumerate(Ne) : 
            RAO_PIT= f"{PIT}/Motion/Global_E{E}_RAO_{PIT_type}_{beta}d{distance}.00.dat"
            PIT_files.append(RAO_PIT) 
        plot_RAO_multiple_E(RAO_REF, RAO_BEM, PIT_files, RAO_PIT_REF, DOF, distance, Nb, Nw, titre, Plot_phase, max_X, Ne, MultiBody)




if N3 : 
    param_d=1
    limite=400
    print("Comparaison sur RAO")
    while param_d<=limite : 
        print("\n ------------Distance=", param_d, "-----------------")
        distance=f"d{param_d}"
        # BEM_file=f"BEM_cylinder_{PIT_type}/BEM_cylinder_Nb{Nb}_{config}_{distance}"

        titre=f"{titre}, d={param_d}m"

        RAO_BEM = f"{BEM}/BEM_cylinder_Nb{Nb}_{config}_d{param_d}/Motion/RAO.dat"
        # RAO_BEM = f"{BEM}/BEM_barge_d{param_d}/Motion/RAO.dat"
        RAO_PIT= f"{PIT}/Motion/Global_{test}RAO_{PIT_type}_{beta}d{param_d}.00.dat"

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


# if PLOT_REF : 
#     d_graphs= [[1], [16], [64]]
    
    
#     print("\n ------------ PLOT with REF -----------------")

#     for j, vect_d in enumerate(d_graphs): 
#         files = [RAO_REF]
#         filesPIT = [RAO_REF]
#         names = ['REF']
#         namesPIT = ['REF']
#         distances = [0]

#         param_d=vect_d
#         print("\n ------------Distances=", param_d, "-----------------")
#         for i, d in enumerate(param_d) : 
#             # RAO_BEM = f"{BEM}/BEM_barge_Nb{Nb}{config}_d{d}/Motion/RAO.dat" 
#             RAO_BEM = f"{BEM}/BEM_barge_d{d}/Motion/RAO.dat" 
#             RAO_PIT= f"{PIT}/Motion/Global_{test}RAO_{PIT_type}_d{d}.00.dat"
#             files.append(RAO_BEM) 
#             filesPIT.append(RAO_PIT) 
#             distances.append(d)
#             names.append('BEM') 
#             namesPIT.append('PIT') 
#         plot_multiple_RAO(files, names, distances, Nb, Nw, titre, Plot_phase, max_X)
#         plot_multiple_RAO(filesPIT, namesPIT, distances, Nb, Nw, titre, Plot_phase, max_X)
