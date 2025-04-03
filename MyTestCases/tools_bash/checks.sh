#!/bin/bash
# Charger les variables depuis config.sh


echo "CHECKS"
# Checks
DOF_ones=$(echo "${DOF[@]}" | tr ' ' '\n' | grep -c '^1$')
if [[ $DOF_ones -ne $N_DOF ]]; then
    echo "Error : The number of '1' in DOF ($DOF_ones) does not match N_DOF ($N_DOF) !" >&2
    exit 1 
fi
FORCES_ones=$(echo "${FORCES[@]}" | tr ' ' '\n' | grep -c '^1$')
if [[ $FORCES_ones -ne $N_Forces ]]; then
    echo "Error : The number of '1' in FORCES ($FORCES_ones) does not match N_Forces ($N_Forces) !" >&2
    exit 1 
fi

for ((i=1; i<=IT_Nb; i++)); do
    if ! declare -p "Bcoord$i" &>/dev/null; then
        echo "Error : Bcoord$i does not exist !" >&2
        exit 1 
    fi
done

for ((i=1; i<=BEM_Nb; i++)); do
    if ! declare -p "CdG$i" &>/dev/null; then
        echo "Error : CdG$i does not exist !" >&2
        exit 1 
    fi
    if ! declare -p "mesh_file$i" &>/dev/null; then
        echo "Error : mesh_file$i does not exist !" >&2
        exit 1  
    fi
done

if [ -z "$w_min" ] || (( $(echo "$w_min == 0" | bc -l) )); then
    echo "Error : w_min can not be 0 !"
    exit 1  
fi


calculate_distance() {
    local x1=$1
    local y1=$2
    local x2=$3
    local y2=$4
    result=$(echo "sqrt(($x2 - $x1)^2 + ($y2 - $y1)^2)" | bc -l)
    result_rounded=$(echo "scale=1; $result / 1" | bc)
    echo $result_rounded
}
epsilon=0.001
if [ "$run_IT" -eq 1 ]; then
    for ((i=1; i<=IT_Nb; i++)); do
        for ((j=i+1; j<=IT_Nb; j++)); do
            eval "coord1=(\"\${Bcoord$i[@]}\")"
            eval "coord2=(\"\${Bcoord$j[@]}\")"
            x1=${coord1[0]}
            y1=${coord1[1]}
            x2=${coord2[0]}
            y2=${coord2[1]}

            distance=$(calculate_distance $x1 $y1 $x2 $y2)
            # echo "Comparaison distance <= 2 * cylR : $distance <= 2 * $cylR-$epsilon"
            if (( $(echo "$distance <= 2 * $cylR-$epsilon" | bc -l) )); then
                echo "Error : The bodies ($i et $j) are too close. Distance = $distance < (2 * $cylR)."
                exit 1  
            fi
        done
    done
fi
