#!/bin/bash
echo "------------- CONFIG -----------------"

Dossier_Project="${Nom_projet}_Nb${BEM_Nb}_${config_type}_d${distance}"

lenght_0=30
lenght_A=$(echo "2*(15*$BEM_Nb + ($BEM_Nb-1)*$distance)" | bc -l)
lenght_B=$(echo "2*(15*$BEM_Nb + $distance)" | bc -l)

echo "Nb = " $BEM_Nb "Type : " $config_type
if [[ $BEM_Nb -ge 2 ]]; then  
    x_t1=0
    y_t1=0
    if [[ "$config_type" == "X" ]]; then 
        x_t2=$(echo "2 * $Rcyl + $distance" | bc -l)
        y_t2=$y_t1
        lenghtX=$lenght_A
        lenghtY=$lenght_0
    elif [[ "$config_type" == "Y" ]]; then 
        y_t2=$(echo "2 * $Rcyl + $distance" | bc -l)
        x_t2=$x_t1
        lenghtY=$lenght_A
        lenghtX=$lenght_0
    fi
fi
if [[ $BEM_Nb -ge 3 ]]; then  
    if [[ "$config_type" == "X" ]]; then 
        x_t3=$(echo "-1*$x_t2" | bc -l)
        y_t3=$y_t2
    elif [[ "$config_type" == "Y" ]]; then 
        y_t3=$(echo "-1*$y_t2" | bc -l)
        x_t3=$x_t2
    elif [[ "$config_type" == "T" ]]; then 
        x_t1=$(echo "scale=8;(2 * $Rcyl + $distance)*sqrt(3)/2" | bc -l)
        y_t1=0
        y_t2=$(echo "(2 * $Rcyl + $distance)/2" | bc -l)
        x_t2=0
        y_t3=$(echo "-1*$y_t2"  | bc -l)
        x_t3=0
        lenghtX=$lenght_B
        lenghtY=$lenght_B
        
    elif [[ "$config_type" == "Tinv" ]]; then 
        x_t2=$(echo "scale=8;(2 * $Rcyl + $distance)*sqrt(3)/2" | bc -l)
        y_t2=$(echo "(2 * $Rcyl + $distance)/2" | bc -l)
        x_t3=$x_t2
        y_t3=$(echo "-1*$y_t2" | bc -l) 
        lenghtX=$lenght_B
        lenghtY=$lenght_B
    fi
fi
if [[ $BEM_Nb -ge 4 ]]; then 
    if [[ "$config_type" == "C" ]]; then 
        x_t1=$(echo "-(2 * $Rcyl + $distance)/2" | bc -l)
        y_t1=$x_t1
        x_t2=$x_t1
        y_t2=$(echo "-1*$x_t2"  | bc -l)
        x_t3=$y_t2
        y_t3=$x_t3
        x_t4=$x_t3
        y_t4=$y_t1
        lenghtX=$lenght_B
        lenghtY=$lenght_B
    fi 
    
fi 
if [[ $FS -eq 0 ]]; then 
    lenghtX=0
    lenghtY=0
fi 
translate1=($x_t1 $y_t1 0.0)
translate2=($x_t2 $y_t2 0.0)
translate3=($x_t3 $y_t3 0.0)
translate4=($x_t4 $y_t4 0.0)

# echo "scale FS" $lenghtX $lenghtY
for ((i=1; i<=BEM_Nb; i++)); do
    eval mesh_file$i=\$mesh_file1
    eval "CdG$i=(\"\${CdG1[@]}\")"
done 

# echo "CdG1=(${CdG1[@]})"
# echo "CdG2=(${CdG2[@]})"
# echo "CdG3=(${CdG3[@]})"
# echo "CdG4=(${CdG4[@]})"
