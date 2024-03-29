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
// File      : Cut_cell_convection_auxiliaire.cpp
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#include <Cut_cell_convection_auxiliaire.h>
#include <Cut_cell_FT_Disc.h>
#include <IJK_Thermal_base.h>
#include <Process.h>
#include <stat_counters.h>
#include <IJK_Navier_Stokes_tools.h>

void Cut_cell_convection_auxiliaire::add_convection_dying_cells(const Cut_field_vector& cut_field_velocity, Cut_field_scalar& cut_field_temperature)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();

  int statut_diphasique = static_cast<int>(cut_cell_disc.STATUT_DIPHASIQUE::MOURRANT);
  int index_min = cut_cell_disc.get_statut_diphasique_value_index(statut_diphasique);
  int index_max = cut_cell_disc.get_statut_diphasique_value_index(statut_diphasique+1);
  for (int index = index_min; index < index_max; index++)
    {
      int n = cut_cell_disc.get_n_from_statut_diphasique_index(index);

      Int3 ijk = cut_cell_disc.get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      if (!cut_cell_disc.within_ghost(i, j, k, 1, 1))
        continue;

      double old_indicatrice = cut_cell_disc.get_interfaces().I(i,j,k);
      double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);
      assert(cut_cell_disc.get_interfaces().est_pure(next_indicatrice));
      int phase = 1 - (int)next_indicatrice; // phase de la cellule mourrante

      double temperature_centre = (phase == 0) ? cut_field_temperature.diph_v_(n) : cut_field_temperature.diph_l_(n);

      double flux[6] = {0};
      for (int num_face = 0; num_face < 6; num_face++)
        {
          int dir = num_face%3;
          int decalage = num_face/3;
          int sign = decalage*2 -1;

          int di = decalage*(dir == 0);
          int dj = decalage*(dir == 1);
          int dk = decalage*(dir == 2);
          int di_decale = sign*(dir == 0);
          int dj_decale = sign*(dir == 1);
          int dk_decale = sign*(dir == 2);

          double old_indicatrice_decale = cut_cell_disc.get_interfaces().I(i+di_decale,j+dj_decale,k+dk_decale);
          double next_indicatrice_decale = cut_cell_disc.get_interfaces().In(i+di_decale,j+dj_decale,k+dk_decale);
          bool decale_also_dying = (cut_cell_disc.get_interfaces().devient_pure(old_indicatrice_decale, next_indicatrice_decale) && ((int)(1 - next_indicatrice_decale) == phase));
          bool decale_also_nascent = (cut_cell_disc.get_interfaces().devient_diphasique(old_indicatrice_decale, next_indicatrice_decale) && ((int)(1 - old_indicatrice_decale) == phase));
          bool decale_smaller = (phase == 0) ? cut_cell_disc.int_indicatrice(1 - old_indicatrice_decale) < cut_cell_disc.int_indicatrice(1 - old_indicatrice) : cut_cell_disc.int_indicatrice(old_indicatrice_decale) < cut_cell_disc.int_indicatrice(old_indicatrice);
          if (decale_also_nascent || (decale_also_dying && decale_smaller))
            {
              flux[num_face] = 0;
            }
          else
            {
              double temperature_decale = (phase == ((int)next_indicatrice_decale))*cut_field_temperature.pure_(i+di_decale,j+dj_decale,k+dk_decale);
              int n_decale = cut_cell_disc.get_n(i+di_decale, j+dj_decale, k+dk_decale);
              if (n_decale >= 0)
                {
                  temperature_decale = (phase == 0) ? cut_field_temperature.diph_v_(n_decale) : cut_field_temperature.diph_l_(n_decale);
                }

              int n_face = cut_cell_disc.get_n_face(num_face, n, i, j, k);
              if (n_face >= 0)
                {
                  double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir) : cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir);
                  double velocity = (phase == 0) ? cut_field_velocity.diph_v_(n_face, dir) : cut_field_velocity.diph_l_(n_face, dir);
                  double temperature = (sign*velocity < 0) ? temperature_decale : temperature_centre;
                  flux[num_face] = sign*surface_efficace*temperature*velocity;
                }
              else
                {
                  double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().I(i+di,j+dj,k+dk) : cut_cell_disc.get_interfaces().I(i+di,j+dj,k+dk);
                  double velocity = cut_field_velocity.pure_[dir](i+di,j+dj,k+dk);
                  double temperature = (sign*velocity < 0) ? temperature_decale : temperature_centre;
                  flux[num_face] = sign*surface_efficace*temperature*velocity;
                }
            }
        }

      double quantite_totale = (phase == 0) ? (1 - old_indicatrice)*cut_field_temperature.diph_v_(n) : old_indicatrice*cut_field_temperature.diph_l_(n);

      int direction_positive = quantite_totale > 0;
      double flux_max_positif = 0.;
      double flux_max_negatif = 0.;
      for (int num_face = 0; num_face < 6; num_face++)
        {
          flux_max_positif = std::max(flux_max_positif, flux[num_face]);
          flux_max_negatif = std::max(flux_max_negatif, -flux[num_face]);
        }
      assert((flux_max_negatif != 0.) || (flux_max_positif != 0.));
      direction_positive = (flux_max_negatif < flux_max_positif/10.) ? 1 : direction_positive;
      direction_positive = (flux_max_positif < flux_max_negatif/10.) ? 0 : direction_positive;

      double somme_flux = 0;
      for (int num_face = 0; num_face < 6; num_face++)
        {
          flux[num_face] = direction_positive ? std::max(0., flux[num_face]) : std::min(0., flux[num_face]);
          somme_flux += flux[num_face];
        }

      assert((flux[0] != 0) || (flux[1] != 0) || (flux[2] != 0) || (flux[3] != 0) || (flux[4] != 0) || (flux[5] != 0));

      for (int num_face = 0; num_face < 6; num_face++)
        {
          if (flux[num_face] != 0.)
            {
              int dir = num_face%3;
              int decalage = num_face/3;
              int sign = decalage*2 -1;

              int di_decale = sign*(dir == 0);
              int dj_decale = sign*(dir == 1);
              int dk_decale = sign*(dir == 2);

              int n_decale = cut_cell_disc.get_n(i+di_decale, j+dj_decale, k+dk_decale);
              double old_indicatrice_decale = cut_cell_disc.get_interfaces().I(i+di_decale,j+dj_decale,k+dk_decale);
              if (n_decale >= 0)
                {
                  if (phase == 0)
                    {
                      cut_field_temperature.diph_v_(n_decale) += (flux[num_face]/somme_flux)*quantite_totale/(1 - old_indicatrice_decale);
                      cut_field_temperature.diph_v_(n) -= (flux[num_face]/somme_flux)*quantite_totale/(1 - cut_cell_disc.get_interfaces().I(i,j,k));
                    }
                  else
                    {
                      cut_field_temperature.diph_l_(n_decale) += (flux[num_face]/somme_flux)*quantite_totale/old_indicatrice_decale;
                      cut_field_temperature.diph_l_(n) -= (flux[num_face]/somme_flux)*quantite_totale/cut_cell_disc.get_interfaces().I(i,j,k);
                    }
                }
              else
                {
                  assert((int)old_indicatrice_decale == 1 - (int)(1 - old_indicatrice_decale));
                  assert(phase == (int)old_indicatrice_decale);
                  cut_field_temperature.pure_(i+di_decale,j+dj_decale,k+dk_decale) += (flux[num_face]/somme_flux)*quantite_totale;
                  if (phase == 0)
                    {
                      cut_field_temperature.diph_v_(n) -= (flux[num_face]/somme_flux)*quantite_totale/(1 - cut_cell_disc.get_interfaces().I(i,j,k));
                    }
                  else
                    {
                      cut_field_temperature.diph_l_(n) -= (flux[num_face]/somme_flux)*quantite_totale/cut_cell_disc.get_interfaces().I(i,j,k);
                    }
                }
            }
        }
    }
}

void Cut_cell_convection_auxiliaire::add_convection_small_nascent_cells(const Cut_field_vector& cut_field_velocity, Cut_field_scalar& cut_field_temperature)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();

  int statut_diphasique_naissant = static_cast<int>(cut_cell_disc.STATUT_DIPHASIQUE::NAISSANT);
  int statut_diphasique_petit = static_cast<int>(cut_cell_disc.STATUT_DIPHASIQUE::DESEQUILIBRE_FINAL);
  assert(statut_diphasique_petit == statut_diphasique_naissant + 1);
  int index_min = cut_cell_disc.get_statut_diphasique_value_index(statut_diphasique_naissant);
  int index_max = cut_cell_disc.get_statut_diphasique_value_index(statut_diphasique_petit+1);
  for (int index = index_min; index < index_max; index++)
    {
      int n = cut_cell_disc.get_n_from_statut_diphasique_index(index);

      Int3 ijk = cut_cell_disc.get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      if (!cut_cell_disc.within_ghost(i, j, k, 1, 1))
        continue;

      double old_indicatrice = cut_cell_disc.get_interfaces().I(i,j,k);
      double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);
      int est_naissant = cut_cell_disc.get_interfaces().est_pure(old_indicatrice);
      int phase = est_naissant ? 1 - (int)old_indicatrice : ((cut_cell_disc.get_interfaces().below_small_threshold(next_indicatrice)) ? 1 : 0); // phase de la cellule petite ou naissante
      assert(est_naissant || cut_cell_disc.get_interfaces().next_below_small_threshold_for_phase(phase, cut_cell_disc.get_interfaces().I(i,j,k), next_indicatrice));

      if (est_naissant)
        {
          assert(((phase == 0) ? cut_field_temperature.diph_v_(n) : cut_field_temperature.diph_l_(n)) == 0.); // La cellule n'est normalement pas deja remplie
        }
      else
        {
          assert(((phase == 0) ? cut_field_temperature.diph_v_(n) : cut_field_temperature.diph_l_(n)) != 0.); // La cellule est normalement deja remplie
        }

      double flux[6] = {0};
      for (int num_face = 0; num_face < 6; num_face++)
        {
          int dir = num_face%3;
          int decalage = num_face/3;
          int sign = decalage*2 -1;

          int di = decalage*(dir == 0);
          int dj = decalage*(dir == 1);
          int dk = decalage*(dir == 2);
          int di_decale = sign*(dir == 0);
          int dj_decale = sign*(dir == 1);
          int dk_decale = sign*(dir == 2);

          double old_indicatrice_decale = cut_cell_disc.get_interfaces().I(i+di_decale,j+dj_decale,k+dk_decale);
          double next_indicatrice_decale = cut_cell_disc.get_interfaces().In(i+di_decale,j+dj_decale,k+dk_decale);
          bool decale_also_dying = (cut_cell_disc.get_interfaces().devient_pure(old_indicatrice_decale, next_indicatrice_decale) && ((int)(1 - next_indicatrice_decale) == phase));
          bool decale_nascent = (cut_cell_disc.get_interfaces().devient_diphasique(old_indicatrice_decale, next_indicatrice_decale) && ((int)(1 - old_indicatrice_decale) == phase));
          bool decale_small = cut_cell_disc.get_interfaces().next_below_small_threshold_for_phase(phase, old_indicatrice_decale, next_indicatrice_decale);
          bool decale_smaller = (phase == 0) ? cut_cell_disc.int_indicatrice(1 - next_indicatrice_decale) < cut_cell_disc.int_indicatrice(1 - next_indicatrice) : cut_cell_disc.int_indicatrice(next_indicatrice_decale) < cut_cell_disc.int_indicatrice(next_indicatrice);
          if (decale_also_dying || (est_naissant && decale_nascent && decale_smaller) || ((!est_naissant) && decale_nascent) || ((!est_naissant) && decale_small && decale_smaller))
            {
              flux[num_face] = 0;
            }
          else
            {
              double temperature_decale = (phase == ((int)next_indicatrice_decale))*cut_field_temperature.pure_(i+di_decale,j+dj_decale,k+dk_decale);
              int n_decale = cut_cell_disc.get_n(i+di_decale, j+dj_decale, k+dk_decale);
              if (n_decale >= 0)
                {
                  temperature_decale = (phase == 0) ? cut_field_temperature.diph_v_(n_decale) : cut_field_temperature.diph_l_(n_decale);
                }

              int n_face = cut_cell_disc.get_n_face(num_face, n, i, j, k);
              if (n_face >= 0)
                {
                  double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir) : cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir);
                  double velocity = (phase == 0) ? cut_field_velocity.diph_v_(n_face, dir) : cut_field_velocity.diph_l_(n_face, dir);
                  double temperature = temperature_decale;
                  flux[num_face] = sign*surface_efficace*temperature*velocity;
                }
              else
                {
                  double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().In(i+di,j+dj,k+dk) : cut_cell_disc.get_interfaces().In(i+di,j+dj,k+dk);
                  double velocity = cut_field_velocity.pure_[dir](i+di,j+dj,k+dk);
                  double temperature = temperature_decale;
                  flux[num_face] = sign*surface_efficace*temperature*velocity;
                }
            }
        }

      double temperature_centre = est_naissant ? 0. : ((phase == 0) ? cut_field_temperature.diph_v_(n) : cut_field_temperature.diph_l_(n));
      double temperature_remplissage = temperature_remplissage_(n);
      double quantite_totale = (phase == 0) ? (1 - next_indicatrice)*(temperature_remplissage - temperature_centre) : next_indicatrice*(temperature_remplissage - temperature_centre);

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
          continue;
        }

      direction_positive = (flux_max_negatif < flux_max_positif/10.) ? 1 : direction_positive;
      direction_positive = (flux_max_positif < flux_max_negatif/10.) ? 0 : direction_positive;

      double somme_flux = 0;
      for (int num_face = 0; num_face < 6; num_face++)
        {
          flux[num_face] = direction_positive ? std::max(0., flux[num_face]) : std::min(0., flux[num_face]);
          somme_flux += flux[num_face];
        }

      assert((flux[0] != 0) || (flux[1] != 0) || (flux[2] != 0) || (flux[3] != 0) || (flux[4] != 0) || (flux[5] != 0));

      for (int num_face = 0; num_face < 6; num_face++)
        {
          if (flux[num_face] != 0.)
            {
              int dir = num_face%3;
              int decalage = num_face/3;
              int sign = decalage*2 -1;

              int di_decale = sign*(dir == 0);
              int dj_decale = sign*(dir == 1);
              int dk_decale = sign*(dir == 2);

              int n_decale = cut_cell_disc.get_n(i+di_decale, j+dj_decale, k+dk_decale);
              double next_indicatrice_decale = cut_cell_disc.get_interfaces().In(i+di_decale,j+dj_decale,k+dk_decale);
              if (n_decale >= 0)
                {
                  if (phase == 0)
                    {
                      cut_field_temperature.diph_v_(n_decale) -= (flux[num_face]/somme_flux)*quantite_totale/(1 - next_indicatrice_decale);
                      cut_field_temperature.diph_v_(n) += (flux[num_face]/somme_flux)*quantite_totale/(1 - cut_cell_disc.get_interfaces().In(i,j,k));
                    }
                  else
                    {
                      cut_field_temperature.diph_l_(n_decale) -= (flux[num_face]/somme_flux)*quantite_totale/next_indicatrice_decale;
                      cut_field_temperature.diph_l_(n) += (flux[num_face]/somme_flux)*quantite_totale/cut_cell_disc.get_interfaces().In(i,j,k);
                    }
                }
              else
                {
                  assert((int)next_indicatrice_decale == 1 - (int)(1 - next_indicatrice_decale));
                  assert(phase == (int)next_indicatrice_decale);
                  cut_field_temperature.pure_(i+di_decale,j+dj_decale,k+dk_decale) -= (flux[num_face]/somme_flux)*quantite_totale;
                  if (phase == 0)
                    {
                      cut_field_temperature.diph_v_(n) += (flux[num_face]/somme_flux)*quantite_totale/(1 - cut_cell_disc.get_interfaces().In(i,j,k));
                    }
                  else
                    {
                      cut_field_temperature.diph_l_(n) += (flux[num_face]/somme_flux)*quantite_totale/cut_cell_disc.get_interfaces().In(i,j,k);
                    }
                }
            }
        }
    }
}

void Cut_cell_convection_auxiliaire::calcule_temperature_remplissage_ponderation_voisin(const Cut_field_vector& cut_field_velocity, Cut_field_scalar& cut_field_temperature)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();

  int statut_diphasique_naissant = static_cast<int>(cut_cell_disc.STATUT_DIPHASIQUE::NAISSANT);
  int statut_diphasique_petit = static_cast<int>(cut_cell_disc.STATUT_DIPHASIQUE::DESEQUILIBRE_FINAL);
  assert(statut_diphasique_petit == statut_diphasique_naissant + 1);
  int index_min = cut_cell_disc.get_statut_diphasique_value_index(statut_diphasique_naissant);
  int index_max = cut_cell_disc.get_statut_diphasique_value_index(statut_diphasique_petit+1);
  for (int index = index_min; index < index_max; index++)
    {
      int n = cut_cell_disc.get_n_from_statut_diphasique_index(index);

      Int3 ijk = cut_cell_disc.get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      if (!cut_cell_disc.within_ghost(i, j, k, 1, 1))
        continue;

      double old_indicatrice = cut_cell_disc.get_interfaces().I(i,j,k);
      double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);
      int est_naissant = cut_cell_disc.get_interfaces().est_pure(old_indicatrice);
      int phase = est_naissant ? 1 - (int)old_indicatrice : ((cut_cell_disc.get_interfaces().below_small_threshold(next_indicatrice)) ? 1 : 0); // phase de la cellule petite ou naissante
      assert(est_naissant || cut_cell_disc.get_interfaces().next_below_small_threshold_for_phase(phase, cut_cell_disc.get_interfaces().I(i,j,k), next_indicatrice));

      if (est_naissant)
        {
          assert(((phase == 0) ? cut_field_temperature.diph_v_(n) : cut_field_temperature.diph_l_(n)) == 0.); // La cellule n'est normalement pas deja remplie
        }
      else
        {
          assert(((phase == 0) ? cut_field_temperature.diph_v_(n) : cut_field_temperature.diph_l_(n)) != 0.); // La cellule est normalement deja remplie
        }

      double temp[6] = {0};
      double flux[6] = {0};
      for (int num_face = 0; num_face < 6; num_face++)
        {
          int dir = num_face%3;
          int decalage = num_face/3;
          int sign = decalage*2 -1;

          int di = decalage*(dir == 0);
          int dj = decalage*(dir == 1);
          int dk = decalage*(dir == 2);
          int di_decale = sign*(dir == 0);
          int dj_decale = sign*(dir == 1);
          int dk_decale = sign*(dir == 2);

          double next_indicatrice_decale = cut_cell_disc.get_interfaces().In(i+di_decale,j+dj_decale,k+dk_decale);

          double temperature_decale = (phase == ((int)next_indicatrice_decale))*cut_field_temperature.pure_(i+di_decale,j+dj_decale,k+dk_decale);
          int n_decale = cut_cell_disc.get_n(i+di_decale, j+dj_decale, k+dk_decale);
          if (n_decale >= 0)
            {
              temperature_decale = (phase == 0) ? cut_field_temperature.diph_v_(n_decale) : cut_field_temperature.diph_l_(n_decale);
            }

          int n_face = cut_cell_disc.get_n_face(num_face, n, i, j, k);
          if (n_face >= 0)
            {
              double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir) : cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir);
              double velocity = (phase == 0) ? cut_field_velocity.diph_v_(n_face, dir) : cut_field_velocity.diph_l_(n_face, dir);
              double temperature = temperature_decale;
              flux[num_face] = sign*surface_efficace*temperature*velocity;
              temp[num_face] = temperature;
            }
          else
            {
              double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().In(i+di,j+dj,k+dk) : cut_cell_disc.get_interfaces().In(i+di,j+dj,k+dk);
              double velocity = cut_field_velocity.pure_[dir](i+di,j+dj,k+dk);
              double temperature = temperature_decale;
              flux[num_face] = sign*surface_efficace*temperature*velocity;
              temp[num_face] = temperature;
            }
        }

      double temperature_remplissage = (temp[0]*std::abs(flux[0]) + temp[1]*std::abs(flux[1]) + temp[2]*std::abs(flux[2]) + temp[3]*std::abs(flux[3]) + temp[4]*std::abs(flux[4]) + temp[5]*std::abs(flux[5]))/((temp[0]*std::abs(flux[0]) != 0.)*std::abs(flux[0]) + (temp[1]*std::abs(flux[1]) != 0.)*std::abs(flux[1]) + (temp[2]*std::abs(flux[2]) != 0.)*std::abs(flux[2]) + (temp[3]*std::abs(flux[3]) != 0.)*std::abs(flux[3]) + (temp[4]*std::abs(flux[4]) != 0.)*std::abs(flux[4]) + (temp[5]*std::abs(flux[5]) != 0.)*std::abs(flux[5]));
      temperature_remplissage_(n) = temperature_remplissage;
    }
}

void Cut_cell_convection_auxiliaire::calcule_temperature_remplissage_semi_lagrangien(double timestep, double lambda_liquid, double lambda_vapour, const IJK_Field_double& flux_interface_ns, const Cut_field_vector& cut_field_velocity, Cut_field_scalar& cut_field_temperature)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();
  const IJK_Field_double& surface_interface_old = cut_cell_disc.get_interfaces().get_surface_interface_old();

  double dx = cut_cell_disc.get_splitting().get_grid_geometry().get_constant_delta(0);
  double dy = cut_cell_disc.get_splitting().get_grid_geometry().get_constant_delta(1);
  double dz = cut_cell_disc.get_splitting().get_grid_geometry().get_constant_delta(2);

  int statut_diphasique_naissant = static_cast<int>(cut_cell_disc.STATUT_DIPHASIQUE::NAISSANT);
  int statut_diphasique_petit = static_cast<int>(cut_cell_disc.STATUT_DIPHASIQUE::DESEQUILIBRE_FINAL);
  assert(statut_diphasique_petit == statut_diphasique_naissant + 1);
  int index_min = cut_cell_disc.get_statut_diphasique_value_index(statut_diphasique_naissant);
  int index_max = cut_cell_disc.get_statut_diphasique_value_index(statut_diphasique_petit+1);
  for (int index = index_min; index < index_max; index++)
    {
      int n = cut_cell_disc.get_n_from_statut_diphasique_index(index);

      Int3 ijk = cut_cell_disc.get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);
      int est_naissant = cut_cell_disc.get_interfaces().est_pure(cut_cell_disc.get_interfaces().I(i,j,k));
      int phase = est_naissant ? 1 - (int)cut_cell_disc.get_interfaces().I(i,j,k) : ((cut_cell_disc.get_interfaces().below_small_threshold(next_indicatrice)) ? 1 : 0); // phase de la cellule petite ou naissante

      double normal_x = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,0);
      double normal_y = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,1);
      double normal_z = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,2);

      double velocity_x = cut_cell_disc.get_interfaces().get_vitesse_deplacement_interface()(n,0);
      double velocity_y = cut_cell_disc.get_interfaces().get_vitesse_deplacement_interface()(n,1);
      double velocity_z = cut_cell_disc.get_interfaces().get_vitesse_deplacement_interface()(n,2);

      double next_bary_x = dx*cut_cell_disc.get_interfaces().get_barycentre_next(0, phase, i,j,k, cut_cell_disc.get_interfaces().I(i,j,k), next_indicatrice);
      double next_bary_y = dy*cut_cell_disc.get_interfaces().get_barycentre_next(1, phase, i,j,k, cut_cell_disc.get_interfaces().I(i,j,k), next_indicatrice);
      double next_bary_z = dz*cut_cell_disc.get_interfaces().get_barycentre_next(2, phase, i,j,k, cut_cell_disc.get_interfaces().I(i,j,k), next_indicatrice);

      double next_bary_deplace_x = next_bary_x - velocity_x*timestep;
      double next_bary_deplace_y = next_bary_y - velocity_y*timestep;
      double next_bary_deplace_z = next_bary_z - velocity_z*timestep;

      int i_old = i + static_cast <int>(std::floor(next_bary_deplace_x/dx));
      int j_old = j + static_cast <int>(std::floor(next_bary_deplace_y/dy));
      int k_old = k + static_cast <int>(std::floor(next_bary_deplace_z/dz));
      assert(i-i_old >= -1 && i-i_old <= 1);
      assert(j-j_old >= -1 && j-j_old <= 1);
      assert(k-k_old >= -1 && k-k_old <= 1);

      int n_old = cut_cell_disc.get_n(i_old,j_old,k_old);
      double old_indicatrice = cut_cell_disc.get_interfaces().I(i_old,j_old,k_old);

      double old_bary_x = dx*cut_cell_disc.get_interfaces().get_barycentre_old(0, phase, i_old,j_old,k_old, old_indicatrice, cut_cell_disc.get_interfaces().In(i_old,j_old,k_old));
      double old_bary_y = dy*cut_cell_disc.get_interfaces().get_barycentre_old(1, phase, i_old,j_old,k_old, old_indicatrice, cut_cell_disc.get_interfaces().In(i_old,j_old,k_old));
      double old_bary_z = dz*cut_cell_disc.get_interfaces().get_barycentre_old(2, phase, i_old,j_old,k_old, old_indicatrice, cut_cell_disc.get_interfaces().In(i_old,j_old,k_old));

      assert(n_old >= 0);
      double temperature_old = (phase == 0) ? cut_field_temperature.diph_v_(n_old) : cut_field_temperature.diph_l_(n_old);
      double temperature_next = (phase == 0) ? cut_field_temperature.diph_v_(n) : cut_field_temperature.diph_l_(n);
      if (temperature_next == 0.)
        {
          assert(temperature_old != 0);
        }
      assert(temperature_old != 0);

      double lambda = (phase == 0) ? lambda_liquid : lambda_vapour;
      double dTdn = (flux_interface_ns(i,j,k) == 0.) ? 0. : flux_interface_ns(i,j,k)/(lambda * surface_interface_old(i,j,k));
      assert((flux_interface_ns(i,j,k) == 0.) || (surface_interface_old(i,j,k) != 0));
      temperature_remplissage_(n) = temperature_old + dTdn * ((next_bary_deplace_x - old_bary_x)*normal_x + (next_bary_deplace_y - old_bary_y)*normal_y + (next_bary_deplace_z - old_bary_z)*normal_z);
    }
}

