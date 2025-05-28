import numpy as np
import matplotlib.pyplot as plt
import itertools


def plot_global_error(config, titre, test):
    # config="Nb2_Y"
    filename = f"EcartsRelatifs_{test}{config}.dat" 
    data = np.loadtxt(filename, delimiter=None, skiprows=0) 

    d = data[:, 0]
    Madd, Crad, Fex_abs, Fex_ph = data[:, 1], data[:, 3], data[:, 5], data[:, 7]

    # Liste des noms et des données
    y_values = [Madd, Crad, Fex_abs, Fex_ph]
    labels = ["Added_Mass", "Damping", "Fex_abs", "Fex_phase"]
    # Générer et sauvegarder chaque graphique
    for i in range(4):
        plt.figure(figsize=(8, 6))
        plt.plot(d, y_values[i], marker='o', linestyle='-', label=labels[i])
        plt.xscale("log", base=2) 
        plt.xlabel("Distance")
        plt.ylabel(labels[i])
        plt.title(f"Global relative errors of {labels[i]} depending on the distance \n ({titre})")
        plt.legend()
        plt.grid()

        # Sauvegarder en PDF
        plt.savefig(f"Errors_{config}_{labels[i]}.png")
        plt.close()  # Ferme la figure pour éviter trop d'ouvertures


def plot_error_Fex(config, titre) : 
    # config="Nb2_Y"
    filename = f"EcartsRelatifs_Fex_abs_{config}.dat" 
    filename2 = f"EcartsRelatifs_Fex_ph_{config}.dat" 
    data = np.loadtxt(filename, delimiter=None, skiprows=1) 
    data2 = np.loadtxt(filename2, delimiter=None, skiprows=1) 

    d = data[:, 0]
    for i in range(0, 12, 2):
        plt.figure(figsize=(8, 6))
        plt.plot(d, data[:,i+1], marker='o', linestyle='-')
        plt.xscale("log", base=2) 
        plt.xlabel("Distance")
        plt.ylabel(f"|Fex|_DOF_{i+1}")
        plt.title(f"Relative errors of |Fex| for DOF_{i+1} depending on the distance \n ({titre})")
        # plt.legend()
        plt.grid()
        # print(i+1, data[:,i+1])
        # Sauvegarder en PDF
        plt.savefig(f"Errors_{config}_Fe_abs_{i+1}_.png")
        plt.close()  # Ferme la figure pour éviter trop d'ouvertures

        plt.figure(figsize=(8, 6))
        plt.plot(d, data2[:,i+1], marker='o', linestyle='-')
        plt.xscale("log", base=2) 
        plt.xlabel("Distance")
        plt.ylabel(f"phase_DOF_{i+1}")
        plt.title(f"Relative errors of Fex phase for DOF_{i+1} depending on the distance \n({titre})")
        # plt.legend()
        plt.grid()

        # Sauvegarder en PDF
        plt.savefig(f"Errors_{config}_Fe_ph_{i+1}_.png")
        plt.close()  # Ferme la figure pour éviter trop d'ouvertures

# Fonction pour charger les matrices depuis un fichier
def load_matrices_from_file(file_path):
    matrices = []
    distances = []
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        matrix = []
        dist = None
        
        for line in lines:
            if line.strip() == "":
                # Une ligne vide sépare les matrices
                if matrix:
                    matrices.append(np.array(matrix))
                    distances.append(dist)
                    matrix = []  # Réinitialiser la matrice pour la prochaine
                continue
            else:
                values = list(map(float, line.split()))
                if dist is None:
                    dist = values[0]  # La distance est la première valeur de la première colonne
                matrix.append(values[1:])  # Ajouter les autres valeurs sans la première colonne (distance)
        
        # Ajouter la dernière matrice et distance si elles existent
        if matrix:
            matrices.append(np.array(matrix))
            distances.append(dist)
    
    return matrices, distances

# Fonction pour tracer un coefficient (i, j) en fonction de la distance
def plot_coefficient_vs_distance(matrices, i, j, nom, config):
    coef_values = []
    dist_values = []
    
    for matrix, dist in zip(matrices, distances):
        coef_values.append(matrix[i-1, j-1])  # Récupérer le coefficient M(i,j)
        # dist_values.append(dist)  # Ajouter la distance associée
    for k in range(9):
        dist_values.append(2**k)  # Ajouter la distance associée
    # Tracer le graphique
    plt.plot(dist_values, coef_values, marker='o', label=f'M({i},{j})')
    plt.xscale("log", base=2) 
    plt.xlabel('Distance')
    plt.ylabel(f'Coef({i},{j})')
    plt.title(f'{nom}_{config} ')
    plt.grid(True)
    plt.legend()
    plt.savefig(f"Errors_{nom}_{config}_{i}_{j}_.png")
    plt.close()  # Ferme la figure pour éviter trop d'ouvertures

def plot_many_coefficient_vs_distance(matrices, indices_coef, nom, config):
    indice=""
    linestyle_cycle = itertools.cycle(['-', '--', '-.', ':'])
    colors = itertools.cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    markers = itertools.cycle(['o', 's', 'D', 'v', '*', '^'])
    for nb_coef in indices_coef : 
        print(nb_coef)
        i=nb_coef[0]
        j=nb_coef[1]
        coef_values = []
        dist_values = []
        for matrix, dist in zip(matrices, distances):
            coef_values.append(matrix[i-1, j-1])  # Récupérer le coefficient M(i,j)
        for k in range(0, 9):
            dist_values.append(2**k)  # Ajouter la distance associée
        # Tracer le graphique
        plt.plot(dist_values, coef_values, color=next(colors),
             marker=next(markers),
             linestyle=next(linestyle_cycle),
             linewidth=2,
             markersize=6,
             alpha=0.9, label=f'M({i},{j})')
        indice=f"{indice}_{i}-{j}"
    plt.xscale("log", base=2) 
    plt.xlabel('Distance')
    plt.ylabel(f'{nom}')
    plt.title(f'{nom}_{config} ')
    plt.grid(True)
    plt.legend()
    plt.savefig(f"Errors_{config}_{nom}{indice}.png")
    plt.close()  # Ferme la figure pour éviter trop d'ouvertures

def plot_many_error_Fex(config, indices_coef, titre, test) : 
    # config="Nb2_Y"
    filename = f"EcartsRelatifs_{test}Fex_abs_{config}.dat" 
    filename2 = f"EcartsRelatifs_{test}Fex_ph_{config}.dat" 
    data = np.loadtxt(filename, delimiter=None, skiprows=1) 
    data2 = np.loadtxt(filename2, delimiter=None, skiprows=1) 
    indice=""
    linestyle_cycle = itertools.cycle(['-', '--', '-.', ':'])
    vivid_colors = [
    "#e6194b",  # rouge vif
    "#3cb44b",  # vert vif
    "#0082c8",  # bleu vif
    "#f58231",  # orange
    "#911eb4",  # violet
    "#000000",  # noir
]
    # colors = itertools.cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    colors = itertools.cycle(vivid_colors)
    markers = itertools.cycle(['o', 's', 'D', 'v', '*', '^'])
    d = data[:, 0]
    # for i in range(0, 12, 2):
    plt.figure(figsize=(8, 6))
    for i in indices_coef : 
        plt.plot(d, data[:,i],color=next(colors),
             marker=next(markers),
             linestyle=next(linestyle_cycle),
             linewidth=2,
             markersize=6,
             alpha=0.9, label=f'DOF_{i}')
        indice=f"{indice}_{i}"
    plt.xscale("log", base=2) 
    plt.xlabel("Distance")
    plt.ylabel(f"|Fex|")
    plt.title(f"Relative errors of |Fex| depending on the distance \n ({titre})")
    plt.legend()
    plt.grid()
    # print(i+1, data[:,i+1])
    # Sauvegarder en PDF
    plt.savefig(f"Errors_{config}_Fex_abs{indice}.png")
    plt.close()  # Ferme la figure pour éviter trop d'ouvertures

    plt.figure(figsize=(8, 6))
    for i in indices_coef :
        plt.plot(d, data2[:,i], color=next(colors),
             marker=next(markers),
             linestyle=next(linestyle_cycle),
             linewidth=2,
             markersize=6,
             alpha=0.9, label=f'DOF_{i}')
        plt.xscale("log", base=2) 
    plt.xlabel("Distance")
    plt.ylabel(f"phase_Fex")
    # plt.ylim(0, 0.3)
    plt.title(f"Relative errors of Fex phase depending on the distance \n({titre})")
    plt.legend()
    plt.grid()

    # Sauvegarder en PDF
    plt.savefig(f"Errors_{config}_Fex_phase{indice}.png")
    plt.close()  # Ferme la figure pour éviter trop d'ouvertures

def plot_many_error_RAO(config, indices_coef, titre, test) : 
    # config="Nb2_Y"
    filename = f"EcartsRelatifs_{test}RAO_abs_{config}.dat" 
    filename2 = f"EcartsRelatifs_{test}RAO_ph_{config}.dat" 
    data = np.loadtxt(filename, delimiter=None, skiprows=1) 
    data2 = np.loadtxt(filename2, delimiter=None, skiprows=1) 
    indice=""
    linestyle_cycle = itertools.cycle(['-', '--', '-.', ':'])
    vivid_colors = [
    "#e6194b",  # rouge vif
    "#3cb44b",  # vert vif
    "#0082c8",  # bleu vif
    "#f58231",  # orange
    "#911eb4",  # violet
    "#000000",  # noir
]
    # colors = itertools.cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    colors = itertools.cycle(vivid_colors)
    markers = itertools.cycle(['o', 's', 'D', 'v', '*', '^'])
    d = data[:, 0]
    # for i in range(0, 12, 2):
    plt.figure(figsize=(8, 6))
    for i in indices_coef : 
        plt.plot(d, data[:,i],color=next(colors),
             marker=next(markers),
             linestyle=next(linestyle_cycle),
             linewidth=2,
             markersize=6,
             alpha=0.9, label=f'DOF_{i}')
        indice=f"{indice}_{i}"
    plt.xscale("log", base=2) 
    plt.xlabel("Distance")
    plt.ylabel(f"|RAO|")
    plt.title(f"Relative errors of |RAO| depending on the distance \n ({titre})")
    plt.legend()
    plt.grid()
    # print(i+1, data[:,i+1])
    # Sauvegarder en PDF
    plt.savefig(f"Errors_{config}_RAO_abs{indice}.png")
    plt.close()  # Ferme la figure pour éviter trop d'ouvertures

    plt.figure(figsize=(8, 6))
    for i in indices_coef :
        plt.plot(d, data2[:,i], color=next(colors),
             marker=next(markers),
             linestyle=next(linestyle_cycle),
             linewidth=2,
             markersize=6,
             alpha=0.9, label=f'DOF_{i}')
        plt.xscale("log", base=2) 
    plt.xlabel("Distance")
    plt.ylabel(f"phase_RAO")
    # plt.ylim(0, 0.3)
    plt.title(f"Relative errors of RAO phase depending on the distance \n({titre})")
    plt.legend()
    plt.grid()

    # Sauvegarder en PDF
    plt.savefig(f"Errors_{config}_RAO_phase{indice}.png")
    plt.close()  # Ferme la figure pour éviter trop d'ouvertures


config="Nb2_X"
test=""
All_error=False
Fex_error=True
C_error=True
RAO_error=True
# i=3
titre=f"Cylinder - {config}, Nw=20, Nbeta=11, Ndof=6, Ndir=1"


if All_error:
    plot_global_error(config, titre, test)

for j, i in enumerate([1, 3, 5]) : 
    if Fex_error:
        # plot_many_error_Fex(config, [i, i+6, i+12], titre)
        plot_many_error_Fex(config, [i, i+6], titre, test)
    if RAO_error:
        plot_many_error_RAO(config, [i, i+6], titre, test)

    file_path = f'EcartsRelatifs_{test}CA_{config}.dat'  # Remplacez par le chemin de votre fichier
    file_path2 = f'EcartsRelatifs_{test}CM_{config}.dat'  # Remplacez par le chemin de votre fichier
    if C_error:
        matrices, distances = load_matrices_from_file(file_path)
        matrices2, distances = load_matrices_from_file(file_path2)
        # Tracer le coefficient choisi en fonction de la distance
        # plot_coefficient_vs_distance(matrices, i+1, j+1, "CA", config)
        # plot_coefficient_vs_distance(matrices2, i+1, j+1, "CM", config)
        coefs=[[i, i], [i+6, i+6]]
        # coefs=[[i, i], [i+2, i+2], [i+6, i+6], [i+8, i+8]]
        plot_many_coefficient_vs_distance(matrices, coefs, "Damping", config)
        plot_many_coefficient_vs_distance(matrices2, coefs, "Added_Mass", config)
        # coefs=[[i, i+6], [i+6, i+12], [i, i+12]]
        coefs=[[i, i+6]]
        # coefs=[[5, 11], [5, 9], [9, 11]]
        plot_many_coefficient_vs_distance(matrices, coefs, "Damping", config)
        plot_many_coefficient_vs_distance(matrices2, coefs, "Added_Mass", config)
