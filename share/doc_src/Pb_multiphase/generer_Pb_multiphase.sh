#!/bin/bash

# Generation documentation des modeles format note CEA
#DOC="TrioCFD_Modeles_clsDM2S"

# Generation documentation des modeles format document utilisateur
DOC="TrioCFD_Pb_multiphase"

# Generation of the pdf report
pdflatex ${DOC}.tex
bibtex ${DOC}
pdflatex ${DOC}.tex
pdflatex ${DOC}.tex

# Installation of the final pdf report
cp ${DOC}.pdf ../../doc

# Cleaning
for ext in aux bbl blg idx log lot toc bcf run.xml; do
    rm "${DOC}.${ext}" 2>/dev/null
done
rm "${DOC}-blx.bib" 2>/dev/null

exit 0
