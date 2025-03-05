#!/bin/bash

cd Calcul

extrait_coupe Prepare.data SONDE_VISC_TURB
extrait_coupe Prepare.data SONDE_EPS
extrait_coupe Prepare.data SONDE_K
extrait_coupe Prepare.data SONDE_VIT

python ../propertiesGeometry.py
python ../ligneTableau.py 
python ../courbes_reichardt.py Moyennes_spatiales_vitesse_rho_mu
python ../dernierTemps.py Moyennes_spatiales_vitesse_rho_mu
python ../dernierTemps.py Moyennes_spatiales_nut

#recuperation des valeurs de v,k,eps a injecter en entree du vrai calcul 
fichier_V=`ls -art | grep -i "pb_VITESSE_PERIO" | tail -1`
fichier_K=`ls -art | grep -i "pb_K_EPS_PERIO"   | tail -1`

#preparation du repertoire pour le calcul final sur le vrai maillage
cp $fichier_V FICHIER_V_PREPARE
cp $fichier_K FICHIER_K_PREPARE
mv u_tau.dat u_tau.dat.perio
mv ligneTableau.dat ligneTableau.dat.perio


awk '{print $3;}' Prepare.TU > NEW
head -7 NEW > NEW2
tail -1 NEW2 > NEW
mv NEW temps_calcul
rm NEW2		