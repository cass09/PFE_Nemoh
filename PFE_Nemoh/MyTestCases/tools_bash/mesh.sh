#!/bin/bash

echo "--------------Mesh-----------------"
export LC_NUMERIC=C  # Important pour éviter les virgules comme séparateurs décimaux

# Calculation of Nnodes and Npanels
for ((i=1; i<=BEM_Nb; i++)); do
    # echo "Ibody = " $i
    mesh_fileI="mesh_file$i"
    mesh_file="${!mesh_fileI}"  
    cp "MESH/$mesh_file" "$Dossier_Project"

    mapfile -t indices < <(awk '{if (($1 < 1e-6 && $1 > -1e-6) && ($2 < 1e-6 && $2 > -1e-6) && ($3 < 1e-6 && $3 > -1e-6) && ($4 < 1e-6 && $4 > -1e-6)) print NR}' "$Dossier_Project/$mesh_file")
    index1=${indices[0]}
    index2=${indices[1]}

    # Calculer Nnoeuds et Npanels
    Nnoeuds=$((index1 - 2))
    Npanels=$((index2 - index1 - 1))

    # echo " $Nnoeuds, $Npanels"

    # Récupération dynamique du tableau de translation
    translate_var="translate$i"
    eval "translate=(\"\${$translate_var[@]}\")"

    tx=${translate[0]}
    ty=${translate[1]}
    tz=${translate[2]}
    # echo "Translation : tx=$tx, ty=$ty, tz=$tz"
    if [ "$BEM_Nb" -gt 1 ]; then
        translated_file="${mesh_file%.*}_$i.${mesh_file##*.}"
    else
        translated_file=$mesh_file
    fi

    if (( $(echo "$tx == 0.0 && $ty == 0.0 && $tz == 0.0" | bc -l) )); then     
        # echo "Translation nulle, copie simple du fichier."
        if [ "$translated_file" != "$mesh_file" ]; then
            cp "$Dossier_Project/$mesh_file" "$Dossier_Project/$translated_file"
        fi
        
        eval mesh_file$i=\"$(basename "$translated_file")\"
        continue
    fi
    awk -v tx="$tx" -v ty="$ty" -v tz="$tz" -v n_nodes="$Nnoeuds" -v index1="$index1" -v index2="$index2" '
    NR == 1 {
        print $1, $2
        next
    }
    NR > 1 && NR <= (1 + n_nodes) {
        # Translating node coordinates
        printf "%d %.6f %.6f %.6f\n", $1, $2 + tx, $3 + ty, $4 + tz
        next
    }
    NR == index1 {
        print "0 0 0 0"
        next
    }
    NR > index1 && NR <= index2 {
        print $1, $2, $3, $4
    }
    ' "$Dossier_Project/$mesh_file" > "$Dossier_Project/$translated_file"

    cp "$Dossier_Project/$mesh_file" "$translated_file"
    eval mesh_file$i=\"$(basename "$translated_file")\"
    rm "$Dossier_Project/$mesh_file"

    # Center of gravity
    cdg_var="CdG$i"
    eval "cdg=(\"\${$cdg_var[@]}\")"
    CdGx=$(echo "${cdg[0]} + $tx" | bc -l)
    CdGy=$(echo "${cdg[1]} + $ty" | bc -l)
    CdGz=$(echo "${cdg[2]} + $tz" | bc -l)
    eval $cdg_var="($CdGx $CdGy $CdGz)"
    # echo $tx $ty $tz
done


rm -f barge_*











































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

