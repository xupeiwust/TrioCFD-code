#!/bin/bash

cd Calcul

extrait_coupe mixedBC.data SONDE_VIT
extrait_coupe mixedBC.data SONDE_VIT2
extrait_coupe mixedBC.data SONDE_VIT3
extrait_coupe mixedBC.data SONDE_VIT4
extrait_coupe mixedBC.data SONDE_VISC_TURB1
extrait_coupe mixedBC.data SONDE_VISC_TURB2


python ../propertiesGeometry.py
python ../postInOut.py

