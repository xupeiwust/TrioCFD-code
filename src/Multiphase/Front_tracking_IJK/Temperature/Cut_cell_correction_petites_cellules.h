/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
/////////////////////////////////////////////////////////////////////////////
//
// File      : Cut_cell_correction_petites_cellules.h
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Cut_cell_correction_petites_cellules_included
#define Cut_cell_correction_petites_cellules_included

enum class CORRECTION_PETITES_CELLULES : int
{
  DIRECTION_PRIVILEGIEE,   // Certaines directions sont privilegiees (le sens oppose)
  DIRECTION_PRIVILEGIEE_2, // Certaines directions sont privilegiees (le sens positif)
  CORRECTION_SYMETRIQUE,   // Pas de directions privilegiees
  DIRECTION_PRIVILEGIEE_AVEC_LIMITATION,   // Certaines directions sont privilegiees (le sens oppose)
  DIRECTION_PRIVILEGIEE_AVEC_LIMITATION_2, // Certaines directions sont privilegiees (le sens positif)
  CORRECTION_SYMETRIQUE_AVEC_LIMITATION,   // Pas de directions privilegiees
};

class Cut_cell_correction_petites_cellules
{
public:
  static void modification_flux_petites_cellules(CORRECTION_PETITES_CELLULES correction_petites_cellules, double quantite_totale, double flux[6]);
  static void limitation_flux_avec_flux_max(CORRECTION_PETITES_CELLULES correction_petites_cellules, double quantite_totale, double somme_flux, double flux_max[6], double flux[6]);
  static double calcul_somme_flux(const double flux[6]);
};


#endif /* Cut_cell_correction_petites_cellules_included */
