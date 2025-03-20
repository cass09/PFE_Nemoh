#!/bin/bash
if [ -z "$1" ]; then
    echo "Aucun nom de dossier fourni. Utilisation du nom de dossier par défaut."
    nom_dossier="test_script"  # Nom de dossier par défaut
else
    nom_dossier="$1"  # Utiliser l'argument passé dans le terminal
fi
# echo "Mise en place"
# Définir les chemins des dossiers source et destination

cp -r "/home/cassandra/Documents/PFE_MOREnergy/Nemoh_myVersion/BEM_Nemoh_matlab/$nom_dossier" "../MyTestCases"

cp -r "RUN"/* "$nom_dossier"  
# echo "Le dossier créé est : $nom_dossier" >> toto
echo ".PHONY: run_$nom_dossier clean_$nom_dossier" >> Makefile
echo " run_$nom_dossier:" >> Makefile
echo "	\$(MAKE) -C \$(testdir)/$nom_dossier/ run" >> Makefile
echo "clean_$nom_dossier:" >> Makefile
echo "	\$(MAKE) -C \$(testdir)/$nom_dossier/ clean" >> Makefile
echo "add_subdirectory($nom_dossier)" >> CMakeLists.txt

echo "Terminé ! -> run_$nom_dossier"