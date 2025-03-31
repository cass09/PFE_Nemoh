#!/bin/bash
# Charger les variables depuis config.sh
source Nemoh_project.sh
Origine=(0.0 0.0 0.0)

# TO DO : 
# - pour Nemoh.cal : checks (bodies no overlap)
# - calcul Nnoeuds et Npanels très lent
# problem memoire quand dof<6

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

if [ ! -d "$Dossier_Project" ]; then
    mkdir "$Dossier_Project"
else 
    if [ -z "$run_BEM" ] || [ "$run_BEM" -eq 0 ]; then
        echo "Error : the project does not exist, BEM must be run !"
        exit 1  
    fi 
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


cp -r "RUN"/* "$Dossier_Project"  
# Calculation of Nnodes and Npanels
for ((i=1; i<=BEM_Nb; i++)); do
    mesh_fileI="mesh_file$i"
    mesh_file="${!mesh_fileI}"  
    cp "MESH/$mesh_file" "$Dossier_Project"

    if [ "$mesh_file" = "Cylinder.dat" ]; then
        Nnoeuds=540
        Npanels=300
    elif [ "$mesh_file" = "barge.dat" ]; then
        Nnoeuds=264
        Npanels=250
    else 
        Nnoeuds=0
        Npanels=0
        found=0

        while IFS= read -r line; do
            # Suppression des espaces en début et fin (sans modifier les espaces au milieu)
            trimmed_line=$(echo "$line" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

            # Vérifie si la ligne est "0          0.00          0.00          0.00" avec espaces multiples
            # if [[ "$trimmed_line" =~ ^0[[:space:]]+0\.00[[:space:]]+0\.00[[:space:]]+0\.00$ ]]; then
            if [[ "$trimmed_line" =~ ^[[:space:]]*0\  ]]; then           
                if [[ $found -eq 0 ]]; then
                    found=1
                    ((Nnoeuds--))
                    # echo "Première section : $Nnoeuds lignes"
                fi
            elif [[ $found -eq 0 ]]; then
                ((Nnoeuds++))
            else
                ((Npanels++))
            fi
        done < "MESH/$mesh_file"
    fi
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

# Nemoh.cal creation
NemohFile="$Dossier_Project/Nemoh.cal"

echo "--- Environment ------------------------------------------------------------------------------------------------------------------"> $NemohFile
echo "$RHO					! RHO 		! KG/M**3 	! Fluid specific volume">> $NemohFile
echo "$G					! G			! M/S**2	! Gravity">> $NemohFile
echo "$DEPTH						! DEPTH			! M		! Water depth">> $NemohFile
echo "$XEFF	$YEFF					! XEFF YEFF		! M		! Wave measurement point">> $NemohFile

echo "--- Description of floating bodies -----------------------------------------------------------------------------------------------">> $NemohFile
echo "$BEM_Nb						! Number of bodies">> $NemohFile
for ((i=1; i<=$BEM_Nb; i++)); do
    CdGi="CdG$i"
    eval "CdG=(\"\${$CdGi[@]}\")"
    mesh_fileI="mesh_file$i"
    eval "mesh_file=\$$mesh_fileI"
    echo "--- Body $i -----------------------------------------------------------------------------------------------------------------------">> $NemohFile
    echo "$mesh_file			! Name of mesh file">> $NemohFile
    echo "$Nnoeuds	$Npanels					! Number of points and number of panels 	">> $NemohFile
    echo "$N_DOF					! Number of degrees of freedom">> $NemohFile
    if [ ${DOF[0]} -eq 1 ]; then
        echo "1 1. 0.	0. ${Origine[@]}		! Surge">> $NemohFile
    fi
    if [ ${DOF[1]} -eq 1 ]; then
        echo "1 0. 1.	0. ${Origine[@]}		! Sway">> $NemohFile
    fi
    if [ ${DOF[2]} -eq 1 ]; then
        echo "1 0. 0. 1. ${Origine[@]}		! Heave">> $NemohFile
    fi
    if [ ${DOF[3]} -eq 1 ]; then
        echo "2 1. 0. 0. ${CdG[@]}	! Roll about CdG">> $NemohFile
    fi
    if [ ${DOF[4]} -eq 1 ]; then
        echo "2 0. 1. 0. ${CdG[@]}	! Pitch about CdG">> $NemohFile
    fi
    if [ ${DOF[5]} -eq 1 ]; then
        echo "2 0. 0. 1. ${CdG[@]}	! Yaw about CdG">> $NemohFile
    fi
    echo "$N_Forces						! Number of resulting generalised forces">> $NemohFile
    if [ ${FORCES[0]} -eq 1 ]; then
        echo "1 1. 0.	0. ${Origine[@]}		! Force in x direction">> $NemohFile
    fi 
    if [ ${FORCES[1]} -eq 1 ]; then
        echo "1 0. 1.	0. ${Origine[@]}		! Force in y direction">> $NemohFile
    fi 
    if [ ${FORCES[2]} -eq 1 ]; then
        echo "1 0. 0. 1. ${Origine[@]}		! Force in z direction">> $NemohFile
    fi 
    if [ ${FORCES[3]} -eq 1 ]; then
        echo "2 1. 0. 0. ${CdG[@]}	! Moment force in x direction about CdG">> $NemohFile
    fi 
    if [ ${FORCES[4]} -eq 1 ]; then
        echo "2 0. 1. 0. ${CdG[@]}	! Moment force in y direction about CdG">> $NemohFile
    fi 
    if [ ${FORCES[5]} -eq 1 ]; then
        echo "2 0. 0. 1. ${CdG[@]} ! Moment force in z direction about CdG">> $NemohFile
    fi
    echo "0						! Number of lines of additional information ">> $NemohFile
done

echo "--- Load cases to be solved -------------------------------------------------------------------------------------------------------">> $NemohFile
echo "$w_type $Nw $w_min $w_max		! Freq type 1,2,3=[rad/s,Hz,s], Number of wave frequencies/periods, Min, and Max">> $NemohFile
echo "$BEM_Nbeta $BEM_BetaMin $BEM_BetaMax				! Number of wave directions, Min and Max (degrees)">> $NemohFile

echo "--- Post processing ---------------------------------------------------------------------------------------------------------------">> $NemohFile
echo "$IRF	$dt	$tf				! IRF calculation (0 for no calculation), time step and duration">> $NemohFile
echo "$show_pressure						! Show pressure">> $NemohFile
echo "$KochinN $KochinMin $KochinMax			! Kochin function 		! Number of directions of calculation (0 for no calculations), Min and Max (degrees)">> $NemohFile
echo "$FS_Nx	$FS_Ny	$FS_Lx	$FS_Ly	! Free surface elevation 	! Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction"	>> $NemohFile
echo "$RAO						! Response Amplitude Operator (RAO), 0 no calculation, 1 calculated -> Inertia.cal">> $NemohFile
echo "$w_type_output						! output freq type, 1,2,3=[rad/s,Hz,s]">> $NemohFile

echo "---Interaction Theory---">> $NemohFile
echo "$run_IT			! run IT">> $NemohFile
if [ $run_IT -eq 1 ]; then
    echo "$run_BEM			! run BEM">> $NemohFile
    echo "$cylR $cylNtheta $cylNz		! Interaction Theory Cylindrical Envelop">> $NemohFile
    echo "$IT_Nb			! Nb bodies">> $NemohFile
    for ((i=1; i<=BEM_Nb+1; i++)); do
        Bcoord_var="Bcoord$i"
        eval "coords=(\"\${$Bcoord_var[@]}\")"
        echo "${coords[@]}			! Coord body $i" >> "$NemohFile"
    done
    echo "$IT_Nbeta   $IT_BetaMin   $IT_BetaMax		! wave directions Nb, min, max">> $NemohFile
fi
echo "---QTF---">> $NemohFile
echo "$QTF         				! QTF flag, 1 is calculated ">> $NemohFile

# input_solver.txt creation
SolverFile="$Dossier_Project/input_solver.txt"
echo "$Gauss_N				! Gauss quadrature (GQ) surface integration, N^2 GQ Nodes, specify N(1,4)" > $SolverFile
echo "$eps_zmin			! eps_zmin for determine minimum z of flow and source points of panel, zmin=eps_zmin*body_diameter" >> $SolverFile
echo "$Solver_type 				! 0 GAUSS ELIM.; 1 LU DECOMP.: 2 GMRES	!Linear system solver">> $SolverFile
echo "$Restart $Tol $MaxIter  	! Restart parameter, Relative Tolerance, max iter -> additional input for GMRES">> $SolverFile


# Makefile
if ! grep -q "run_$Dossier_Project:" Makefile; then
    echo ".PHONY: run_$Dossier_Project clean_$Dossier_Project" >> Makefile
    echo "run_$Dossier_Project:" >> Makefile
    echo "	\$(MAKE) -C \$(testdir)/$Dossier_Project/ run" >> Makefile
fi

if ! grep -q "clean_$Dossier_Project:" Makefile; then
    echo "clean_$Dossier_Project:" >> Makefile
    echo "	\$(MAKE) -C \$(testdir)/$Dossier_Project/ clean" >> Makefile
fi

# CMakeLists.txt
if ! grep -q "add_subdirectory($Dossier_Project)" CMakeLists.txt; then
    echo "add_subdirectory($Dossier_Project)" >> CMakeLists.txt
fi