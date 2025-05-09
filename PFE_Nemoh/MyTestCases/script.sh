#!/bin/bash
if [ -z "$1" ]; then
    echo "Aucun nom de dossier fourni. Utilisation du nom de dossier par défaut."
    nom_dossier="test_script"  # Nom de dossier par défaut
else
    nom_dossier="$1"  # Utiliser l'argument passé dans le terminal
fi
# echo "Mise en place"
# Définir les chemins des dossiers source et destination

# cp -r "/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/BEM_Nemoh_matlab/$nom_dossier" "../MyTestCases"

cp -r "RUN"/* "$nom_dossier"  

# Makefile
if ! grep -q "run_$nom_dossier:" Makefile; then
    echo ".PHONY: run_$nom_dossier clean_$nom_dossier" >> Makefile
    echo "run_$nom_dossier:" >> Makefile
    echo "	\$(MAKE) -C \$(testdir)/$nom_dossier/ run" >> Makefile
fi

if ! grep -q "clean_$nom_dossier:" Makefile; then
    echo "clean_$nom_dossier:" >> Makefile
    echo "	\$(MAKE) -C \$(testdir)/$nom_dossier/ clean" >> Makefile
fi

# CMakeLists.txt
if ! grep -q "add_subdirectory($nom_dossier)" CMakeLists.txt; then
    echo "add_subdirectory($nom_dossier)" >> CMakeLists.txt
fi

if ! grep -q "run_$nom_dossier" Allrun.sh; then
    echo "make run_$nom_dossier" >> Allrun.sh
fi

echo "-> run_$nom_dossier in Allrun.sh"