#!/bin/bash
# Charger les variables depuis config.sh
source Nemoh_project.sh
Origine=(0.0 0.0 0.0)

# TO DO : 
# - pour Nemoh.cal : checks (bodies no overlap, R> dim body)
# - if IT_run - > CdG=Origine   ?= Coord1

source tools_bash/checks.sh
if [ $? -ne 0 ]; then
    exit 1
fi

if [ ! -d "$Dossier_Project" ]; then
    mkdir "$Dossier_Project"
else 
    if [ -z "$run_BEM" ] || [ "$run_BEM" -eq 0 ]; then
        echo "Error : the project does not exist, BEM must be run !"
        exit 1  
    fi 
    read -p "The folder '$Dossier_Project' already exist. Overwrite ? (o/n) : " answer
    if [[ "$answer" != "o" && "$answer" != "O" ]]; then
        echo "Cancelled"
        exit 1
    fi
fi

source tools_bash/mesh.sh
source tools_bash/input_files.sh

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

cp -r "RUN"/* "$Dossier_Project"  

if ! grep -q "run_$Dossier_Project " Allrun.sh; then
    echo "make run_$Dossier_Project " >> Allrun.sh
fi
echo "-> run_$Dossier_Project"