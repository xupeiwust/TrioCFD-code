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
// File      : Cut_cell_correction_petites_cellules.cpp
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#include <Cut_cell_correction_petites_cellules.h>
#include <IJK_Thermal_cut_cell.h>
#include <Process.h>

void Cut_cell_correction_petites_cellules::modification_flux_petites_cellules(CORRECTION_PETITES_CELLULES correction_petites_cellules, double quantite_totale, double flux[6])
{
  if (correction_petites_cellules == CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE || correction_petites_cellules == CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE_AVEC_LIMITATION)
    {
      for (int num_face = 0; num_face < 6; num_face++)
        {
          flux[num_face] = -flux[num_face];
        }
    }

  if (correction_petites_cellules == CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE || correction_petites_cellules == CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE_2 || correction_petites_cellules == CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE_AVEC_LIMITATION || correction_petites_cellules == CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE_AVEC_LIMITATION_2)
    {
      int direction_positive = quantite_totale > 0;
      double flux_max_positif = 0.;
      double flux_max_negatif = 0.;
      for (int num_face = 0; num_face < 6; num_face++)
        {
          flux_max_positif = std::max(flux_max_positif, flux[num_face]);
          flux_max_negatif = std::max(flux_max_negatif, -flux[num_face]);
        }
      if ((flux_max_negatif == 0.) && (flux_max_positif == 0.))
        {
          //Cerr << "Warning: for cell " << n << " no change is made due to lack of suitable fluxes; for information temperature_remplissage=" << temperature_remplissage << " and temperature_centre=" << temperature_centre << finl;
          return;
        }

      direction_positive = (flux_max_negatif < flux_max_positif/10.) ? 1 : direction_positive;
      direction_positive = (flux_max_positif < flux_max_negatif/10.) ? 0 : direction_positive;

      for (int num_face = 0; num_face < 6; num_face++)
        {
          flux[num_face] = direction_positive ? std::max(0., flux[num_face]) : std::max(0., -flux[num_face]);
        }
    }
  else if (correction_petites_cellules == CORRECTION_PETITES_CELLULES::CORRECTION_SYMETRIQUE || correction_petites_cellules == CORRECTION_PETITES_CELLULES::CORRECTION_SYMETRIQUE_AVEC_LIMITATION)
    {
      for (int num_face = 0; num_face < 6; num_face++)
        {
          flux[num_face] = std::abs(flux[num_face]);
        }
    }
  else
    {
      Cerr << "Valeur inconnue de correction_petites_cellules" << finl;
      Process::exit();
    }
}

void Cut_cell_correction_petites_cellules::limitation_flux_avec_flux_max(CORRECTION_PETITES_CELLULES correction_petites_cellules, double quantite_totale, double somme_flux, double flux_max[6], double flux[6])
{
  // Je pense que tous les flux et flux_max sont positifs
  assert((flux[0] >= 0 && flux[1] >= 0 && flux[2] >= 0 && flux[3] >= 0 && flux[4] >= 0 && flux[5] >= 0));
  assert((flux_max[0] >= 0 && flux_max[1] >= 0 && flux_max[2] >= 0 && flux_max[3] >= 0 && flux_max[4] >= 0 && flux_max[5] >= 0));

  // Ne fait rien si l'option n'inclue pas une limitation ou egalement
  // si le changmeent total a effectuer est nul
  if (quantite_totale != 0 && (correction_petites_cellules == CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE_AVEC_LIMITATION || correction_petites_cellules == CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE_AVEC_LIMITATION_2 || correction_petites_cellules == CORRECTION_PETITES_CELLULES::CORRECTION_SYMETRIQUE_AVEC_LIMITATION))
    {
      int number_of_unlimited_cells = 0.; // Compte le nombre de flux sans limite
      double excess_fluxes = 0.; // Ce qu'il faut repartir a la suite de la limitation
      double flux_room[6] = {0}; // Espace disponible pour la repartition les flux (0 = pas de limites)
      double total_room = 0.;

      // Les flux reels sont de la forme (flux[num_face]/somme_flux)*quantite_totale
      // On corrige flux_max pour qu'il soit comparable aux flux
      for (int num_face = 0; num_face < 6; num_face++)
        {
          flux_max[num_face] = flux_max[num_face]*std::abs(somme_flux/quantite_totale);
        }

      for (int num_face = 0; num_face < 6; num_face++)
        {
          if ((flux_max[num_face] != 0) && (flux[num_face] > flux_max[num_face]))
            {
              excess_fluxes += flux[num_face] - flux_max[num_face];
              flux[num_face] = flux_max[num_face];
            }
        }

      if (excess_fluxes != 0)
        {
          for (int num_face = 0; num_face < 6; num_face++)
            {
              if (flux_max[num_face] != 0)
                {
                  flux_room[num_face] = flux_max[num_face] - flux[num_face];
                  total_room += flux_room[num_face];
                }
              else if (flux_max[num_face] == 0 && flux[num_face] != 0)
                {
                  number_of_unlimited_cells += 1;
                }
            }

          // Si le compte n'est pas suffisant, les flux sans limites sont mis a contribution
          double left_over_fluxes = std::max(0., excess_fluxes - total_room);
          for (int num_face = 0; num_face < 6; num_face++)
            {
              if (flux_max[num_face] == 0 && flux[num_face] != 0)
                {
                  flux_room[num_face] = left_over_fluxes/number_of_unlimited_cells;
                  total_room += flux_room[num_face];
                }
            }

          // Correction des flux : Aucune correction si pas d'espace disponible, meme si l'exces subsiste
          for (int num_face = 0; num_face < 6; num_face++)
            {
              if (total_room != 0.)
                {
                  flux[num_face] += (flux_room[num_face]/total_room)*std::min(excess_fluxes, total_room);
                }
            }
        }

      for (int num_face = 0; num_face < 6; num_face++)
        {
          assert((flux_max[num_face] == 0.) || (flux[num_face] <= flux_max[num_face]*(1 + 1e-9)));
        }
    }
}

double Cut_cell_correction_petites_cellules::calcul_somme_flux(const double flux[6])
{
  double somme_flux = 0.;
  for (int num_face = 0; num_face < 6; num_face++)
    {
      somme_flux += flux[num_face];
    }
  return somme_flux;
}

