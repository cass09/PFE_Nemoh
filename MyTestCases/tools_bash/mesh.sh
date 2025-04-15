#!/bin/bash
# Charger les variables depuis config.sh

echo "Mesh"


# Calculation of Nnodes and Npanels
for ((i=1; i<=BEM_Nb; i++)); do
    mesh_fileI="mesh_file$i"
    mesh_file="${!mesh_fileI}"  
    cp "MESH/$mesh_file" "$Dossier_Project"


    # Utiliser awk pour détecter les lignes contenant uniquement des zéros
    # Afficher le numéro de la ligne (NR) et la ligne complète ($0)
    # index=$(awk '{if (($1 < 1e-6 && $1 > -1e-6) && ($2 < 1e-6 && $2 > -1e-6) && ($3 < 1e-6 && $3 > -1e-6) && ($4 < 1e-6 && $4 > -1e-6)) print NR}' "$Dossier_Project/$mesh_file")
    # Récupérer les indices dans un tableau
    mapfile -t indices < <(awk '{if (($1 < 1e-6 && $1 > -1e-6) && ($2 < 1e-6 && $2 > -1e-6) && ($3 < 1e-6 && $3 > -1e-6) && ($4 < 1e-6 && $4 > -1e-6)) print NR}' "$Dossier_Project/$mesh_file")

    # Vérification du tableau
    echo "Indices : ${indices[@]}"
    # Récupérer les indices
    index1=${indices[0]}
    index2=${indices[1]}

    # Calculer Nnoeuds et Npanels
    Nnoeuds=$((index1 - 2))
    Npanels=$((index2 - index1 - 1))

    echo " $Nnoeuds, $Npanels"
    # if [ "$mesh_file" = "Cylinder.dat" ]; then
    #     Nnoeuds=540
    #     Npanels=300
    # elif [ "$mesh_file" = "barge.dat" ]; then
    #     Nnoeuds=264
    #     Npanels=250
    # else 
    #     Nnoeuds=0
    #     Npanels=0
    #     found=0

    #     while IFS= read -r line; do
    #         # Suppression des espaces en début et fin (sans modifier les espaces au milieu)
    #         trimmed_line=$(echo "$line" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

    #         # Vérifie si la ligne est "0          0.00          0.00          0.00" avec espaces multiples
    #         # if [[ "$trimmed_line" =~ ^0[[:space:]]+0\.00[[:space:]]+0\.00[[:space:]]+0\.00$ ]]; then
    #         if [[ "$trimmed_line" =~ ^[[:space:]]*0\  ]]; then           
    #             if [[ $found -eq 0 ]]; then
    #                 found=1
    #                 ((Nnoeuds--))
    #                 # echo "Première section : $Nnoeuds lignes"
    #             fi
    #         elif [[ $found -eq 0 ]]; then
    #             ((Nnoeuds++))
    #         else
    #             ((Npanels++))
    #         fi
    #     done < "MESH/$mesh_file"
    # fi



    # total_lines=$(wc -l < "MESH/$mesh_file")
    # # Calculer la somme de Nnoeuds + Npanels + 3
    # expected_lines=$(($Nnoeuds + $Npanels + 3))
    # echo "$expected_lines"
    # # Vérification que le nombre total de lignes dans le fichier est égal à Nnoeuds + Npanels + 3
    # if [[ $total_lines -ne $expected_lines ]]; then
    #     echo "Erreur : Problème dans le calcul du nombre de noeuds et de panels du maillage !" >&2
    #     exit 1 
    # fi
done
