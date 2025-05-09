#!/bin/bash
# Charger les variables depuis config.sh

echo "-------------Writing input files-----------------"
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
echo "$cylR $cylNtheta $cylNz		! Interaction Theory Cylindrical Envelop">> $NemohFile
echo "$run_IT			! run IT">> $NemohFile
if [ $run_IT -eq 1 ]; then
    echo "$run_BEM			! run BEM">> $NemohFile
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
