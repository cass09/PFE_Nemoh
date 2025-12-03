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
        x_t1=$(echo "-($distance)/2" | bc -l)
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
if [[ $BEM_Nb -ge 9 ]]; then 
    if [[ "$config_type" == "C" ]]; then 
        x_t1=0
        y_t1=0
        x_t2=$(echo "($Rcyl + $Ra)" | bc -l)
        y_t2=0
        x_t3=$x_t2
        y_t3=$(echo "-1*$x_t2"  | bc -l)
        x_t4=0
        y_t4=$y_t3
        x_t5=$y_t3
        y_t5=$y_t3
        x_t6=$y_t3
        y_t6=0
        x_t7=$y_t3
        y_t7=$x_t2
        x_t8=0
        y_t8=$x_t2
        x_t9=$x_t2
        y_t9=$x_t2
        lenghtX=0
        lenghtY=0
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
translate5=($x_t5 $y_t5 0.0)
translate6=($x_t6 $y_t6 0.0)
translate7=($x_t7 $y_t7 0.0)
translate8=($x_t8 $y_t8 0.0)
translate9=($x_t9 $y_t9 0.0)
for i in {1..9}; do
    eval "current_translate=(\"\${translate$i[@]}\")"
    echo "translate$i=${current_translate[@]}"
done
# echo "scale FS" $lenghtX $lenghtY
for ((i=1; i<=BEM_Nb; i++)); do
    eval mesh_file$i=\$mesh_file1
    eval "CdG$i=(\"\${CdG1[@]}\")"
done 

# echo "CdG1=(${CdG1[@]})"
# echo "CdG2=(${CdG2[@]})"
# echo "CdG3=(${CdG3[@]})"
# echo "CdG4=(${CdG4[@]})"
