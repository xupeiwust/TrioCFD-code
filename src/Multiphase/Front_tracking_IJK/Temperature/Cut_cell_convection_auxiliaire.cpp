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
#include <IJK_FT_cut_cell.h>
#include <stat_counters.h>
#include <IJK_Navier_Stokes_tools.h>
#include <IJK_Navier_Stokes_tools_cut_cell.h>
#include <Param.h>

Implemente_instanciable_sans_constructeur(Cut_cell_convection_auxiliaire, "Cut_cell_convection_auxiliaire", Objet_U) ;

Cut_cell_convection_auxiliaire::Cut_cell_convection_auxiliaire()
{
  methode_temperature_remplissage_ = METHODE_TEMPERATURE_REMPLISSAGE::NON_INITIALISE;
  convection_petites_cellules_ = CORRECTION_PETITES_CELLULES::CORRECTION_SYMETRIQUE;
}

Sortie& Cut_cell_convection_auxiliaire::printOn(Sortie& os) const
{
  Objet_U::printOn(os);
  return os;
}

Entree& Cut_cell_convection_auxiliaire::readOn(Entree& is)
{
  Param param(que_suis_je());

  param.ajouter("methode_temperature_remplissage", (int*)&methode_temperature_remplissage_);
  param.dictionnaire("non_initialise",(int)METHODE_TEMPERATURE_REMPLISSAGE::NON_INITIALISE);
  param.dictionnaire("ponderation_voisin", (int)METHODE_TEMPERATURE_REMPLISSAGE::PONDERATION_VOISIN);
  param.dictionnaire("ponderation_directionnelle_voisin", (int)METHODE_TEMPERATURE_REMPLISSAGE::PONDERATION_DIRECTIONNELLE_VOISIN);
  param.dictionnaire("semi_lagrangien", (int)METHODE_TEMPERATURE_REMPLISSAGE::SEMI_LAGRANGIEN);
  param.dictionnaire("semi_lagrangien_interpolate", (int)METHODE_TEMPERATURE_REMPLISSAGE::SEMI_LAGRANGIEN_INTERPOLATE);

  param.ajouter("convection_petites_cellules", (int*)&convection_petites_cellules_);
  param.dictionnaire("correction_directe", (int)CORRECTION_PETITES_CELLULES::CORRECTION_DIRECTE);
  param.dictionnaire("direction_privilegiee", (int)CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE);
  param.dictionnaire("direction_privilegiee_2", (int)CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE_2);
  param.dictionnaire("correction_symetrique", (int)CORRECTION_PETITES_CELLULES::CORRECTION_SYMETRIQUE);
  param.dictionnaire("direction_privilegiee_avec_limitation", (int)CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE_AVEC_LIMITATION);
  param.dictionnaire("direction_privilegiee_avec_limitation_2", (int)CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE_AVEC_LIMITATION_2);
  param.dictionnaire("correction_symetrique_avec_limitation", (int)CORRECTION_PETITES_CELLULES::CORRECTION_SYMETRIQUE_AVEC_LIMITATION);

  param.lire_avec_accolades(is);

  return is;
}

void Cut_cell_convection_auxiliaire::compute_flux_dying_cells(const Cut_field_vector& cut_field_total_velocity, Cut_field_scalar& cut_field_temperature)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();

  const IJK_Grid_Geometry& geom = cut_cell_disc.get_splitting().get_grid_geometry();
  const double delta_x = geom.get_constant_delta(0);
  const double delta_y = geom.get_constant_delta(1);
  const double delta_z = geom.get_constant_delta(2);
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(0));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(1));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(2));

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

          double f_dir = select(dir, delta_y*delta_z, delta_x*delta_z, delta_x*delta_y);

          double old_indicatrice_decale = cut_cell_disc.get_interfaces().I(i+di_decale,j+dj_decale,k+dk_decale);
          double next_indicatrice_decale = cut_cell_disc.get_interfaces().In(i+di_decale,j+dj_decale,k+dk_decale);
          bool decale_also_dying = (cut_cell_disc.get_interfaces().devient_pure(old_indicatrice_decale, next_indicatrice_decale) && ((int)(1 - next_indicatrice_decale) == phase));
          bool decale_also_nascent = (cut_cell_disc.get_interfaces().devient_diphasique(old_indicatrice_decale, next_indicatrice_decale) && ((int)(1 - old_indicatrice_decale) == phase));
          bool decale_smaller = (phase == 0) ? (1 - old_indicatrice_decale) < (1 - old_indicatrice) : (old_indicatrice_decale) < (old_indicatrice);
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
                  if (surface_efficace > 0)
                    {
                      double total_velocity = (phase == 0) ? cut_field_total_velocity.diph_v_(n_face, dir) : cut_field_total_velocity.diph_l_(n_face, dir);
                      double temperature = (sign*total_velocity < 0) ? temperature_decale : temperature_centre;
                      flux[num_face] = -sign*f_dir*surface_efficace*temperature*total_velocity;
                    }
                }
              else
                {
                  double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().I(i+di,j+dj,k+dk) : cut_cell_disc.get_interfaces().I(i+di,j+dj,k+dk);
                  if (surface_efficace > 0)
                    {
                      double total_velocity = cut_field_total_velocity.pure_[dir](i+di,j+dj,k+dk);
                      double temperature = (sign*total_velocity < 0) ? temperature_decale : temperature_centre;
                      flux[num_face] = -sign*f_dir*surface_efficace*temperature*total_velocity;
                    }
                }
            }

          flux_naive_(n,num_face) = flux[num_face];
        }

      // Correction du flux si aucun flux non-nul
      if ((flux[0] == 0) && (flux[1] == 0) && (flux[2] == 0) && (flux[3] == 0) && (flux[4] == 0) && (flux[5] == 0))
        {
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

              double f_dir = select(dir, delta_y*delta_z, delta_x*delta_z, delta_x*delta_y);

              double old_indicatrice_decale = cut_cell_disc.get_interfaces().I(i+di_decale,j+dj_decale,k+dk_decale);
              double next_indicatrice_decale = cut_cell_disc.get_interfaces().In(i+di_decale,j+dj_decale,k+dk_decale);
              bool decale_also_dying = (cut_cell_disc.get_interfaces().devient_pure(old_indicatrice_decale, next_indicatrice_decale) && ((int)(1 - next_indicatrice_decale) == phase));
              bool decale_also_nascent = (cut_cell_disc.get_interfaces().devient_diphasique(old_indicatrice_decale, next_indicatrice_decale) && ((int)(1 - old_indicatrice_decale) == phase));
              bool decale_smaller = (phase == 0) ? (1 - old_indicatrice_decale) < (1 - old_indicatrice) : (old_indicatrice_decale) < (old_indicatrice);
              if (decale_also_nascent || (decale_also_dying && decale_smaller))
                {
                  flux[num_face] = 0;
                }
              else
                {
                  int n_face = cut_cell_disc.get_n_face(num_face, n, i, j, k);
                  if (n_face >= 0)
                    {
                      double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir) : cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir);
                      if (surface_efficace > 0)
                        {
                          flux[num_face] = f_dir*surface_efficace;
                        }
                    }
                  else
                    {
                      double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().I(i+di,j+dj,k+dk) : cut_cell_disc.get_interfaces().I(i+di,j+dj,k+dk);
                      if (surface_efficace > 0)
                        {
                          flux[num_face] = f_dir*surface_efficace;
                        }
                    }
                }

              flux_naive_(n,num_face) = flux[num_face];
            }
        }
    }
}

void Cut_cell_convection_auxiliaire::compute_flux_small_nascent_cells(const Cut_field_vector& cut_field_total_velocity, Cut_field_scalar& cut_field_temperature)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();

  const IJK_Grid_Geometry& geom = cut_cell_disc.get_splitting().get_grid_geometry();
  const double delta_x = geom.get_constant_delta(0);
  const double delta_y = geom.get_constant_delta(1);
  const double delta_z = geom.get_constant_delta(2);
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(0));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(1));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(2));

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
          assert(std::abs((phase == 0) ? cut_field_temperature.diph_v_(n) : cut_field_temperature.diph_l_(n)) < 1e-16); // La cellule n'est normalement pas deja remplie
        }
      else
        {
          assert(std::abs((phase == 0) ? cut_field_temperature.diph_v_(n) : cut_field_temperature.diph_l_(n)) > 1e-16); // La cellule est normalement deja remplie
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

          double f_dir = select(dir, delta_y*delta_z, delta_x*delta_z, delta_x*delta_y);

          double old_indicatrice_decale = cut_cell_disc.get_interfaces().I(i+di_decale,j+dj_decale,k+dk_decale);
          double next_indicatrice_decale = cut_cell_disc.get_interfaces().In(i+di_decale,j+dj_decale,k+dk_decale);
          bool decale_also_dying = (cut_cell_disc.get_interfaces().devient_pure(old_indicatrice_decale, next_indicatrice_decale) && ((int)(1 - next_indicatrice_decale) == phase));
          bool decale_nascent = (cut_cell_disc.get_interfaces().devient_diphasique(old_indicatrice_decale, next_indicatrice_decale) && ((int)(1 - old_indicatrice_decale) == phase));
          bool decale_small = cut_cell_disc.get_interfaces().next_below_small_threshold_for_phase(phase, old_indicatrice_decale, next_indicatrice_decale);
          bool decale_smaller = (phase == 0) ? (1 - next_indicatrice_decale) < (1 - next_indicatrice) : (next_indicatrice_decale) < (next_indicatrice);
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
                  if (surface_efficace > 0)
                    {
                      double total_velocity = (phase == 0) ? cut_field_total_velocity.diph_v_(n_face, dir) : cut_field_total_velocity.diph_l_(n_face, dir);
                      double temperature = temperature_decale;
                      flux[num_face] = -sign*f_dir*surface_efficace*temperature*total_velocity;
                    }
                }
              else
                {
                  double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().In(i+di,j+dj,k+dk) : cut_cell_disc.get_interfaces().In(i+di,j+dj,k+dk);
                  if (surface_efficace > 0)
                    {
                      double total_velocity = cut_field_total_velocity.pure_[dir](i+di,j+dj,k+dk);
                      double temperature = temperature_decale;
                      flux[num_face] = -sign*f_dir*surface_efficace*temperature*total_velocity;
                    }
                }
            }

          flux_naive_(n,num_face) = flux[num_face];
        }

      // Correction du flux si aucun flux non-nul
      if ((flux[0] == 0) && (flux[1] == 0) && (flux[2] == 0) && (flux[3] == 0) && (flux[4] == 0) && (flux[5] == 0))
        {
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

              double f_dir = select(dir, delta_y*delta_z, delta_x*delta_z, delta_x*delta_y);

              double old_indicatrice_decale = cut_cell_disc.get_interfaces().I(i+di_decale,j+dj_decale,k+dk_decale);
              double next_indicatrice_decale = cut_cell_disc.get_interfaces().In(i+di_decale,j+dj_decale,k+dk_decale);
              bool decale_also_dying = (cut_cell_disc.get_interfaces().devient_pure(old_indicatrice_decale, next_indicatrice_decale) && ((int)(1 - next_indicatrice_decale) == phase));
              bool decale_nascent = (cut_cell_disc.get_interfaces().devient_diphasique(old_indicatrice_decale, next_indicatrice_decale) && ((int)(1 - old_indicatrice_decale) == phase));
              bool decale_small = cut_cell_disc.get_interfaces().next_below_small_threshold_for_phase(phase, old_indicatrice_decale, next_indicatrice_decale);
              bool decale_smaller = (phase == 0) ? (1 - next_indicatrice_decale) < (1 - next_indicatrice) : (next_indicatrice_decale) < (next_indicatrice);
              if (decale_also_dying || (est_naissant && decale_nascent && decale_smaller) || ((!est_naissant) && decale_nascent) || ((!est_naissant) && decale_small && decale_smaller))
                {
                  flux[num_face] = 0;
                }
              else
                {
                  int n_face = cut_cell_disc.get_n_face(num_face, n, i, j, k);
                  if (n_face >= 0)
                    {
                      double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir) : cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir);
                      if (surface_efficace > 0)
                        {
                          flux[num_face] = f_dir*surface_efficace;
                        }
                    }
                  else
                    {
                      double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().In(i+di,j+dj,k+dk) : cut_cell_disc.get_interfaces().In(i+di,j+dj,k+dk);
                      if (surface_efficace > 0)
                        {
                          flux[num_face] = f_dir*surface_efficace;
                        }
                    }
                }

              flux_naive_(n,num_face) = flux[num_face];
            }
        }
    }
}


void Cut_cell_convection_auxiliaire::add_convection_dying_cells(const Cut_field_vector& cut_field_total_velocity, Cut_field_scalar& cut_field_temperature)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();

  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(0));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(1));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(2));

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

      if (std::abs(temperature_centre) < 1e-24)
        {
          continue;
        }

      double quantite_totale = (phase == 0) ? -(1 - old_indicatrice)*cut_field_temperature.diph_v_(n) : -old_indicatrice*cut_field_temperature.diph_l_(n);


      double flux[6] = {flux_naive_(n,0), flux_naive_(n,1), flux_naive_(n,2), flux_naive_(n,3), flux_naive_(n,4), flux_naive_(n,5)};

      Cut_cell_correction_petites_cellules::modification_flux_petites_cellules(convection_petites_cellules_, quantite_totale, flux);
      double somme_flux = Cut_cell_correction_petites_cellules::calcul_somme_flux(flux);
      assert(somme_flux != 0.);


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
                      cut_field_temperature.diph_v_(n_decale) -= (flux[num_face]/somme_flux)*quantite_totale/(1 - old_indicatrice_decale);
                      cut_field_temperature.diph_v_(n) += (flux[num_face]/somme_flux)*quantite_totale/(1 - cut_cell_disc.get_interfaces().I(i,j,k));
                    }
                  else
                    {
                      cut_field_temperature.diph_l_(n_decale) -= (flux[num_face]/somme_flux)*quantite_totale/old_indicatrice_decale;
                      cut_field_temperature.diph_l_(n) += (flux[num_face]/somme_flux)*quantite_totale/cut_cell_disc.get_interfaces().I(i,j,k);
                    }
                }
              else
                {
                  assert((int)old_indicatrice_decale == 1 - (int)(1 - old_indicatrice_decale));
                  assert(phase == (int)old_indicatrice_decale);
                  cut_field_temperature.pure_(i+di_decale,j+dj_decale,k+dk_decale) -= (flux[num_face]/somme_flux)*quantite_totale;
                  if (phase == 0)
                    {
                      cut_field_temperature.diph_v_(n) += (flux[num_face]/somme_flux)*quantite_totale/(1 - cut_cell_disc.get_interfaces().I(i,j,k));
                    }
                  else
                    {
                      cut_field_temperature.diph_l_(n) += (flux[num_face]/somme_flux)*quantite_totale/cut_cell_disc.get_interfaces().I(i,j,k);
                    }
                }
            }
        }
    }
}

void Cut_cell_convection_auxiliaire::add_convection_small_nascent_cells(const Cut_field_vector& cut_field_total_velocity, Cut_field_scalar& cut_field_temperature)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();

  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(0));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(1));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(2));

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


      double temperature_centre = ((phase == 0) ? cut_field_temperature.diph_v_(n) : cut_field_temperature.diph_l_(n));
      double temperature_remplissage = temperature_remplissage_(n);
      double quantite_totale = (phase == 0) ? (1 - next_indicatrice)*(temperature_remplissage - temperature_centre) : next_indicatrice*(temperature_remplissage - temperature_centre);

      // The dummy value DMINFLOAT+1 indicates that no filling should be made.
      if (temperature_remplissage == DMINFLOAT+1)
        {
          if (phase == 0)
            {
              assert((cut_field_total_velocity.diph_v_(n, 0) == 0) && (cut_field_total_velocity.diph_v_(n, 1) == 0) && (cut_field_total_velocity.diph_v_(n, 2) == 0));
            }
          else
            {
              assert((cut_field_total_velocity.diph_l_(n, 0) == 0) && (cut_field_total_velocity.diph_l_(n, 1) == 0) && (cut_field_total_velocity.diph_l_(n, 2) == 0));
            }
          continue;
        }

      if ((std::abs(temperature_centre - temperature_remplissage)/temperature_remplissage) < 1e-24)
        {
          continue;
        }


      double flux[6] = {flux_naive_(n,0), flux_naive_(n,1), flux_naive_(n,2), flux_naive_(n,3), flux_naive_(n,4), flux_naive_(n,5)};

      // Calcul du flux maximum pour la limitation des flux
      double flux_max[6] = {0};
      double total_velocity[6] = {0};
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
          bool decale_smaller = (phase == 0) ? (1 - next_indicatrice_decale) < (1 - next_indicatrice) : (next_indicatrice_decale) < (next_indicatrice);
          if (decale_also_dying || (est_naissant && decale_nascent && decale_smaller) || ((!est_naissant) && decale_nascent) || ((!est_naissant) && decale_small && decale_smaller))
            {
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
                  if (surface_efficace > 0)
                    {
                      double next_volume_decale = (phase == 0) ? 1 - next_indicatrice_decale : next_indicatrice_decale;
                      double next_volume = (phase == 0) ? 1 - next_indicatrice : next_indicatrice;
                      flux_max[num_face] = (temperature_decale - temperature_centre)/(1/next_volume_decale + 1/next_volume);
                      flux_max[num_face] = flux[num_face] >= 0 ? std::max(0., flux_max[num_face]) : std::max(0., -flux_max[num_face]);
                    }
                  total_velocity[num_face] = (phase == 0) ? cut_field_total_velocity.diph_v_(n_face, dir) : cut_field_total_velocity.diph_l_(n_face, dir);
                }
              else
                {
                  double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().In(i+di,j+dj,k+dk) : cut_cell_disc.get_interfaces().In(i+di,j+dj,k+dk);
                  if (surface_efficace > 0)
                    {
                      double next_volume_decale = (phase == 0) ? 1 - next_indicatrice_decale : next_indicatrice_decale;
                      double next_volume = (phase == 0) ? 1 - next_indicatrice : next_indicatrice;
                      flux_max[num_face] = (temperature_decale - temperature_centre)/(1/next_volume_decale + 1/next_volume);
                      flux_max[num_face] = flux[num_face] >= 0 ? std::max(0., flux_max[num_face]) : std::max(0., -flux_max[num_face]);
                    }
                  total_velocity[num_face] = cut_field_total_velocity.pure_[dir](i+di,j+dj,k+dk);
                }
            }
        }

      if (((total_velocity[0] == 0) && (total_velocity[1] == 0) && (total_velocity[2] == 0) && (total_velocity[3] == 0) && (total_velocity[4] == 0) && (total_velocity[5] == 0)))
        {
          continue;
        }


      Cut_cell_correction_petites_cellules::modification_flux_petites_cellules(convection_petites_cellules_, quantite_totale, flux);
      double somme_flux = Cut_cell_correction_petites_cellules::calcul_somme_flux(flux);
      assert(somme_flux != 0.);

      if (somme_flux == 0.)
        {
          continue;
        }

      Cut_cell_correction_petites_cellules::limitation_flux_avec_flux_max(convection_petites_cellules_, quantite_totale, somme_flux, flux_max, flux);

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

void Cut_cell_convection_auxiliaire::calcule_temperature_remplissage_ponderation_voisin(bool est_directionnel, const Cut_field_vector& cut_field_total_velocity, const Cut_field_scalar& cut_field_temperature)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();

  const IJK_Grid_Geometry& geom = cut_cell_disc.get_splitting().get_grid_geometry();
  const double delta_x = geom.get_constant_delta(0);
  const double delta_y = geom.get_constant_delta(1);
  const double delta_z = geom.get_constant_delta(2);
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(0));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(1));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(2));

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


      double temp[6] = {0};
      double flux[6] = {0};
      double total_velocity[6] = {0};
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

          double f_dir = select(dir, delta_y*delta_z, delta_x*delta_z, delta_x*delta_y);

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
              total_velocity[num_face] = (phase == 0) ? cut_field_total_velocity.diph_v_(n_face, dir) : cut_field_total_velocity.diph_l_(n_face, dir);
              double temperature = temperature_decale;
              flux[num_face] = -sign*f_dir*surface_efficace*temperature*total_velocity[num_face];
              if (est_directionnel)
                {
                  flux[num_face] = (flux[num_face] < 0) ? 0. : flux[num_face];
                }
              temp[num_face] = temperature;
            }
          else
            {
              double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().In(i+di,j+dj,k+dk) : cut_cell_disc.get_interfaces().In(i+di,j+dj,k+dk);
              total_velocity[num_face] = cut_field_total_velocity.pure_[dir](i+di,j+dj,k+dk);
              double temperature = temperature_decale;
              flux[num_face] = -sign*f_dir*surface_efficace*temperature*total_velocity[num_face];
              if (est_directionnel)
                {
                  flux[num_face] = (flux[num_face] < 0) ? 0. : flux[num_face];
                }
              temp[num_face] = temperature;
            }
        }

      // Correction du flux si aucun flux non-nul
      if ((flux[0] == 0) && (flux[1] == 0) && (flux[2] == 0) && (flux[3] == 0) && (flux[4] == 0) && (flux[5] == 0))
        {
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

              double f_dir = select(dir, delta_y*delta_z, delta_x*delta_z, delta_x*delta_y);

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
                  double temperature = temperature_decale;
                  flux[num_face] = f_dir*surface_efficace*temperature;
                  temp[num_face] = temperature;
                }
              else
                {
                  double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().In(i+di,j+dj,k+dk) : cut_cell_disc.get_interfaces().In(i+di,j+dj,k+dk);
                  double temperature = temperature_decale;
                  flux[num_face] = f_dir*surface_efficace*temperature;
                  temp[num_face] = temperature;
                }
            }
        }

      if (((total_velocity[0] == 0) && (total_velocity[1] == 0) && (total_velocity[2] == 0) && (total_velocity[3] == 0) && (total_velocity[4] == 0) && (total_velocity[5] == 0)))
        {
          temperature_remplissage_(n) = DMINFLOAT+1; // Dummy value to indicate that no filling should be made.
        }
      else
        {
          double somme_T_flux = (temp[0]*std::abs(flux[0]) != 0.)*std::abs(flux[0]) + (temp[1]*std::abs(flux[1]) != 0.)*std::abs(flux[1]) + (temp[2]*std::abs(flux[2]) != 0.)*std::abs(flux[2]) + (temp[3]*std::abs(flux[3]) != 0.)*std::abs(flux[3]) + (temp[4]*std::abs(flux[4]) != 0.)*std::abs(flux[4]) + (temp[5]*std::abs(flux[5]) != 0.)*std::abs(flux[5]);
          assert(somme_T_flux != 0.);
          double temperature_remplissage = (temp[0]*std::abs(flux[0]) + temp[1]*std::abs(flux[1]) + temp[2]*std::abs(flux[2]) + temp[3]*std::abs(flux[3]) + temp[4]*std::abs(flux[4]) + temp[5]*std::abs(flux[5]))/somme_T_flux;
          temperature_remplissage_(n) = temperature_remplissage;
        }
    }
}

void Cut_cell_convection_auxiliaire::calcule_temperature_remplissage_semi_lagrangien(double timestep, double lambda_liquid, double lambda_vapour, const IJK_Field_double& flux_interface_ns, const Cut_field_scalar& cut_field_temperature)
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

      double next_bary_x = dx*cut_cell_disc.get_interfaces().get_barycentre(true, 0, phase, i,j,k, cut_cell_disc.get_interfaces().I(i,j,k), next_indicatrice);
      double next_bary_y = dy*cut_cell_disc.get_interfaces().get_barycentre(true, 1, phase, i,j,k, cut_cell_disc.get_interfaces().I(i,j,k), next_indicatrice);
      double next_bary_z = dz*cut_cell_disc.get_interfaces().get_barycentre(true, 2, phase, i,j,k, cut_cell_disc.get_interfaces().I(i,j,k), next_indicatrice);

      double next_bary_deplace_x = next_bary_x - velocity_x*timestep;
      double next_bary_deplace_y = next_bary_y - velocity_y*timestep;
      double next_bary_deplace_z = next_bary_z - velocity_z*timestep;

      int i_old = i + static_cast <int>(std::floor(next_bary_deplace_x/dx));
      int j_old = j + static_cast <int>(std::floor(next_bary_deplace_y/dy));
      int k_old = k + static_cast <int>(std::floor(next_bary_deplace_z/dz));
      assert(i-i_old >= -1 && i-i_old <= 1);
      assert(j-j_old >= -1 && j-j_old <= 1);
      assert(k-k_old >= -1 && k-k_old <= 1);

      // En raison du seuil sur l'indicatrice, il est possible que la cellule_old corresponde a la cellule initiale
      // bien que celle-ci soit consideree comme naissante.
      // Dans ce cas, on utilise un voisin ; qui pourra fournir une temperature de reference.
      if (est_naissant && (i_old == i) && (j_old == j) && (k_old == k))
        {
          double ecart_x = std::abs(next_bary_deplace_x/dx - std::round(next_bary_deplace_x/dx));
          double ecart_y = std::abs(next_bary_deplace_y/dy - std::round(next_bary_deplace_y/dy));
          double ecart_z = std::abs(next_bary_deplace_z/dz - std::round(next_bary_deplace_z/dz));
          int direction_a_modifier = (ecart_x <= ecart_y && ecart_x <= ecart_z) ? 0 : ((ecart_y <= ecart_x && ecart_y <= ecart_z) ? 1 : 2);
          double next_bary_deplace_normalise = select(direction_a_modifier, next_bary_deplace_x/dx, next_bary_deplace_y/dy, next_bary_deplace_z/dz);
          bool ajout_positif = std::floor(next_bary_deplace_normalise) > .5;

          if (direction_a_modifier == 0)
            {
              i_old += (2*ajout_positif - 1);
            }
          else if (direction_a_modifier == 1)
            {
              j_old += (2*ajout_positif - 1);
            }
          else if (direction_a_modifier == 2)
            {
              k_old += (2*ajout_positif - 1);
            }
        }

      int n_old = cut_cell_disc.get_n(i_old,j_old,k_old);
      double old_indicatrice = cut_cell_disc.get_interfaces().I(i_old,j_old,k_old);

      double old_bary_x = dx*((i_old - i) + cut_cell_disc.get_interfaces().get_barycentre(false, 0, phase, i_old,j_old,k_old, old_indicatrice, cut_cell_disc.get_interfaces().In(i_old,j_old,k_old)));
      double old_bary_y = dy*((j_old - j) + cut_cell_disc.get_interfaces().get_barycentre(false, 1, phase, i_old,j_old,k_old, old_indicatrice, cut_cell_disc.get_interfaces().In(i_old,j_old,k_old)));
      double old_bary_z = dz*((k_old - k) + cut_cell_disc.get_interfaces().get_barycentre(false, 2, phase, i_old,j_old,k_old, old_indicatrice, cut_cell_disc.get_interfaces().In(i_old,j_old,k_old)));

      if (n_old < 0)
        {
          Cerr << "Dans Cut_cell_convection_auxiliaire::calcule_temperature_remplissage_semi_lagrangien, on a n_old < 0." << finl;
          Cerr << "Ce cas est normalement pas possible mais pourrait peut-etre survenir a cause de la correction avec 'est_naissant && (i_old == i) && (j_old == j) && (k_old == k)'" << finl;
          Cerr << "Verification : est_naissant=" << est_naissant << " i=" << i << " j=" << j << " k=" << k << " i_old=" << i_old << " j_old=" << j_old << " k_old=" << k_old << finl;
          //Process::exit();
        }

      double temperature_old = (n_old < 0) ? cut_field_temperature.pure_(i_old,j_old,k_old) : ((phase == 0) ? cut_field_temperature.diph_v_(n_old) : cut_field_temperature.diph_l_(n_old));
      double temperature_next = (phase == 0) ? cut_field_temperature.diph_v_(n) : cut_field_temperature.diph_l_(n);
      if (temperature_next == 0.)
        {
          assert(temperature_old != 0);
        }
      assert(temperature_old != 0);

      assert((flux_interface_ns(i,j,k) == 0.) || (surface_interface_old(i,j,k) != 0.)); // La surface ne devrait pas etre nulle si il y a un flux
      double lambda = (phase == 0) ? lambda_vapour : lambda_liquid;
      double dTdn = (flux_interface_ns(i,j,k) == 0.) ? 0. : flux_interface_ns(i,j,k)/(lambda * surface_interface_old(i,j,k));
      assert((flux_interface_ns(i,j,k) == 0.) || (surface_interface_old(i,j,k) != 0));
      temperature_remplissage_(n) = temperature_old + dTdn * ((next_bary_deplace_x - old_bary_x)*normal_x + (next_bary_deplace_y - old_bary_y)*normal_y + (next_bary_deplace_z - old_bary_z)*normal_z);
    }
}

void Cut_cell_convection_auxiliaire::calcule_temperature_remplissage_semi_lagrangien_interpolate(double timestep, const ArrOfDouble& interfacial_temperature, const IJK_Field_double& temperature_ft, const Cut_field_scalar& cut_field_temperature)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();

  int count_status_not_found = 0;
  int count_status_below_10 = 0;
  int count_status_below_100 = 0;
  int count_status_below_1000 = 0;
  int count_status_below_10000 = 0;
  int count_status_below_100000 = 0;
  int count_status_below_1000000 = 0;
  int count_status_above_1000000 = 0;

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

      double velocity_x = cut_cell_disc.get_interfaces().get_vitesse_deplacement_interface()(n,0);
      double velocity_y = cut_cell_disc.get_interfaces().get_vitesse_deplacement_interface()(n,1);
      double velocity_z = cut_cell_disc.get_interfaces().get_vitesse_deplacement_interface()(n,2);

      double next_bary_x = dx*cut_cell_disc.get_interfaces().get_barycentre(true, 0, phase, i,j,k, cut_cell_disc.get_interfaces().I(i,j,k), next_indicatrice);
      double next_bary_y = dy*cut_cell_disc.get_interfaces().get_barycentre(true, 1, phase, i,j,k, cut_cell_disc.get_interfaces().I(i,j,k), next_indicatrice);
      double next_bary_z = dz*cut_cell_disc.get_interfaces().get_barycentre(true, 2, phase, i,j,k, cut_cell_disc.get_interfaces().I(i,j,k), next_indicatrice);

      double next_bary_deplace_x = next_bary_x - velocity_x*timestep;
      double next_bary_deplace_y = next_bary_y - velocity_y*timestep;
      double next_bary_deplace_z = next_bary_z - velocity_z*timestep;

      double coordinates_x = next_bary_deplace_x + (i + cut_cell_disc.get_splitting().get_offset_local(DIRECTION_I))*dx + cut_cell_disc.get_splitting().get_grid_geometry().get_origin(DIRECTION_I);
      double coordinates_y = next_bary_deplace_y + (j + cut_cell_disc.get_splitting().get_offset_local(DIRECTION_J))*dy + cut_cell_disc.get_splitting().get_grid_geometry().get_origin(DIRECTION_J);
      double coordinates_z = next_bary_deplace_z + (k + cut_cell_disc.get_splitting().get_offset_local(DIRECTION_K))*dz + cut_cell_disc.get_splitting().get_grid_geometry().get_origin(DIRECTION_K);
      double coordinates[3] = {coordinates_x, coordinates_y, coordinates_z};

      int status = -2;
      const int tolerate_not_within_tetrahedron = 1;
      double temperature_interpolate = ijk_interpolate_cut_cell_using_interface(false, phase, temperature_ft, cut_field_temperature, interfacial_temperature, coordinates, tolerate_not_within_tetrahedron, status);
      temperature_remplissage_(n) = temperature_interpolate;

      if (status == -1)
        {
          Cerr << "DEBUG status=-1 T_remplissage=" << temperature_remplissage_(n) << finl;
        }

      assert(status > -2);
      if (status == -1)
        {
          count_status_not_found++;
        }
      else if (status < 10)
        {
          count_status_below_10++;
        }
      else if (status < 100)
        {
          count_status_below_100++;
        }
      else if (status < 1000)
        {
          count_status_below_1000++;
        }
      else if (status < 10000)
        {
          count_status_below_10000++;
        }
      else if (status < 100000)
        {
          count_status_below_100000++;
        }
      else if (status < 1000000)
        {
          count_status_below_1000000++;
        }
      else if (status >= 1000000)
        {
          count_status_above_1000000++;
        }
      else
        {
          Cerr << "Un tel statut n'est pas possible." << finl;
          Process::exit();
        }
    }

  Cerr << "Bilan des statuts pour Cut_cell_convection_auxiliaire::calcule_temperature_remplissage_semi_lagrangien_interpolate : " << finl;
  Cerr << "    -1      " << count_status_not_found  << finl;
  Cerr << "   <10      " << count_status_below_10   << finl;
  Cerr << "  <100      " << count_status_below_100  << finl;
  Cerr << " <1000      " << count_status_below_1000 << finl;
  Cerr << " <10000     " << count_status_below_10000 << finl;
  Cerr << " <100000    " << count_status_below_100000 << finl;
  Cerr << " <1000000   " << count_status_below_1000000 << finl;
  Cerr << ">=1000000   " << count_status_above_1000000 << finl;
}

void Cut_cell_convection_auxiliaire::calcule_temperature_face_depuis_centre(double lambda_liquid, double lambda_vapour, const IJK_Field_double& flux_interface_ns, const Cut_field_scalar& cut_field_temperature, FixedVector<FixedVector<IJK_Field_double, 3>, 2>& temperature_face)
{
  const bool next_time = false;

  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();
  const IJK_Field_double& surface_interface = next_time ? cut_cell_disc.get_interfaces().get_surface_interface_next() : cut_cell_disc.get_interfaces().get_surface_interface_old();

  const FixedVector<IJK_Field_double, 3>& old_indicatrice_surfacique_face  = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_face_old();
  const FixedVector<IJK_Field_double, 3>& next_indicatrice_surfacique_face = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_face_next();

  double dx = cut_cell_disc.get_splitting().get_grid_geometry().get_constant_delta(0);
  double dy = cut_cell_disc.get_splitting().get_grid_geometry().get_constant_delta(1);
  double dz = cut_cell_disc.get_splitting().get_grid_geometry().get_constant_delta(2);

  temperature_face[0][0].data()=0;
  temperature_face[0][1].data()=0;
  temperature_face[0][2].data()=0;
  temperature_face[1][0].data()=0;
  temperature_face[1][1].data()=0;
  temperature_face[1][2].data()=0;

  for (int n = 0; n < cut_cell_disc.get_n_tot(); n++)
    {
      Int3 ijk = cut_cell_disc.get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      for (int dir = 0; dir < 3; dir++)
        {
          for (int phase = 0; phase < 2; phase++)
            {
              double old_indicatrice  = cut_cell_disc.get_interfaces().I(i,j,k);
              double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);

              double normal_x = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,0);
              double normal_y = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,1);
              double normal_z = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,2);

              double centre_bary_x = dx*cut_cell_disc.get_interfaces().get_barycentre(next_time, 0, phase, i,j,k, old_indicatrice, next_indicatrice);
              double centre_bary_y = dy*cut_cell_disc.get_interfaces().get_barycentre(next_time, 1, phase, i,j,k, old_indicatrice, next_indicatrice);
              double centre_bary_z = dz*cut_cell_disc.get_interfaces().get_barycentre(next_time, 2, phase, i,j,k, old_indicatrice, next_indicatrice);

              double old_indicatrice_surfacique  = old_indicatrice_surfacique_face[dir](i,j,k);
              double next_indicatrice_surfacique = next_indicatrice_surfacique_face[dir](i,j,k);

              double face_bary_x = dx*cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 0, phase, i,j,k, old_indicatrice_surfacique, next_indicatrice_surfacique);
              double face_bary_y = dy*cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 1, phase, i,j,k, old_indicatrice_surfacique, next_indicatrice_surfacique);
              double face_bary_z = dz*cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 2, phase, i,j,k, old_indicatrice_surfacique, next_indicatrice_surfacique);

              double temperature_centre = (phase == 0) ? cut_field_temperature.diph_v_(n) : cut_field_temperature.diph_l_(n);

              assert((flux_interface_ns(i,j,k) == 0.) || (surface_interface(i,j,k) != 0.)); // La surface ne devrait pas etre nulle si il y a un flux
              double lambda = (phase == 0) ? lambda_vapour : lambda_liquid;
              double dTdn = (flux_interface_ns(i,j,k) == 0.) ? 0. : flux_interface_ns(i,j,k)/(lambda * surface_interface(i,j,k));
              assert((flux_interface_ns(i,j,k) == 0.) || (surface_interface(i,j,k) != 0));

              const DoubleTabFT_cut_cell_vector3& indicatrice_surfacique = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face();
              double indicatrice_surface = ((phase == 0) ? 1 - indicatrice_surfacique(n,dir) : indicatrice_surfacique(n,dir));

              temperature_face[phase][dir](i,j,k) = (indicatrice_surface == 0.) ? 0. : (temperature_centre + dTdn * ((face_bary_x - centre_bary_x)*normal_x + (face_bary_y - centre_bary_y)*normal_y + (face_bary_z - centre_bary_z)*normal_z));
            }
        }
    }
}

void Cut_cell_convection_auxiliaire::calcule_temperature_face_depuis_facette(double lambda_liquid, double lambda_vapour, const ArrOfDouble& interfacial_temperature, const ArrOfDouble& interfacial_phin_ai, const Cut_field_scalar& cut_field_temperature, REF(IJK_FT_cut_cell)& ref_ijk_ft, FixedVector<FixedVector<IJK_Field_double, 3>, 2>& temperature_face_ft, FixedVector<FixedVector<IJK_Field_double, 3>, 2>& temperature_face_ns)
{
  const bool next_time = false;

  // Calcul de T_face sur le maillage FT
  {
    const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();
    const Maillage_FT_IJK& mesh = cut_cell_disc.get_interfaces().maillage_ft_ijk();
    const IntTab& facettes = next_time ? mesh.facettes() : mesh.facettes_old();
    const DoubleTab& sommets = next_time ? mesh.sommets() : mesh.sommets_old();
    const DoubleTab& normale_facettes = next_time ? mesh.get_update_normale_facettes() : mesh.get_normale_facettes_old();
    const ArrOfDouble& surface_facettes = next_time ? mesh.get_update_surface_facettes() : mesh.get_surface_facettes_old();
    const Intersections_Elem_Facettes& intersec = next_time ? mesh.intersections_elem_facettes() : mesh.intersections_elem_facettes_old();
    const IJK_Splitting& s = temperature_face_ft.get_splitting();

    double dx = s.get_grid_geometry().get_constant_delta(DIRECTION_I);
    double dy = s.get_grid_geometry().get_constant_delta(DIRECTION_J);
    double dz = s.get_grid_geometry().get_constant_delta(DIRECTION_K);

    const IJK_Grid_Geometry& geom = s.get_grid_geometry();
    double origin_x = geom.get_origin(DIRECTION_I);
    double origin_y = geom.get_origin(DIRECTION_J);
    double origin_z = geom.get_origin(DIRECTION_K);

    int offset_x = s.get_offset_local(DIRECTION_I);
    int offset_y = s.get_offset_local(DIRECTION_J);
    int offset_z = s.get_offset_local(DIRECTION_K);

    const FixedVector<IJK_Field_double, 3>& old_indicatrice_surfacique_face  = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_face_old();
    const FixedVector<IJK_Field_double, 3>& next_indicatrice_surfacique_face = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_face_next();

    const int ni_ft = temperature_face_ft[0][0].ni();
    const int nj_ft = temperature_face_ft[0][0].nj();
    const int nk_ft = temperature_face_ft[0][0].nk();

    // Initialisation
    {
      for (int k_ft = 0; k_ft < nk_ft; k_ft++)
        {
          for (int j_ft = 0; j_ft < nj_ft; j_ft++)
            {
              for (int i_ft = 0; i_ft < ni_ft; i_ft++)
                {
                  for (int dir = 0; dir < 3; dir++)
                    {
                      for (int phase = 0; phase < 2; phase++)
                        {
                          temperature_face_ft[phase][dir](i_ft, j_ft, k_ft) = 0.;
                        }
                    }
                }
            }
        }
    }

    // Calcul pour les faces coupees par l'interface
    {
      const ArrOfInt& index_elem = intersec.index_elem();
      for (int k_ft = 0; k_ft < nk_ft; k_ft++)
        {
          for (int j_ft = 0; j_ft < nj_ft; j_ft++)
            {
              for (int i_ft = 0; i_ft < ni_ft; i_ft++)
                {
                  if (next_time ? (!cut_cell_disc.get_interfaces().est_pure(cut_cell_disc.get_interfaces().In_ft(i_ft,j_ft,k_ft))) : (!cut_cell_disc.get_interfaces().est_pure(cut_cell_disc.get_interfaces().I_ft(i_ft,j_ft,k_ft))))
                    {
                      double x_centre = (i_ft + offset_x + .5)*dx + origin_x;
                      double y_centre = (j_ft + offset_y + .5)*dy + origin_y;
                      double z_centre = (k_ft + offset_z + .5)*dz + origin_z;
                      int i_ns = cut_cell_disc.get_i_selon_dir(0, x_centre);
                      int j_ns = cut_cell_disc.get_i_selon_dir(1, y_centre);
                      int k_ns = cut_cell_disc.get_i_selon_dir(2, z_centre);

                      int n = cut_cell_disc.get_n(i_ns, j_ns, k_ns);

                      assert(mesh.ref_splitting().valeur() == s);
                      const int num_elem = s.convert_ijk_cell_to_packed(i_ft,j_ft,k_ft);
                      int index = index_elem[num_elem];

                      int facette_la_plus_proche[2][3] = {{-1, -1, -1}, {-1, -1, -1}};
                      double min_distance_facette[2][3] = {{DMAXFLOAT, DMAXFLOAT, DMAXFLOAT}, {DMAXFLOAT, DMAXFLOAT, DMAXFLOAT}};

                      // Boucle sur les facettes qui traversent cet element
                      while (index >= 0)
                        {
                          const Intersections_Elem_Facettes_Data& data = intersec.data_intersection(index);
                          const int fa7 = data.numero_facette_;

                          double coord_facettes[3] = {0};
                          for (int dir = 0; dir < 3; dir++)
                            for (int som = 0; som < 3; som++)
                              coord_facettes[dir] += sommets(facettes(fa7, som), dir);
                          for (int dir = 0; dir < 3; dir++)
                            coord_facettes[dir] /= 3.;

                          for (int dir = 0; dir < 3; dir++)
                            {
                              for (int phase = 0; phase < 2; phase++)
                                {
                                  double old_indicatrice_surfacique  = old_indicatrice_surfacique_face[dir](i_ns,j_ns,k_ns);
                                  double next_indicatrice_surfacique = next_indicatrice_surfacique_face[dir](i_ns,j_ns,k_ns);
                                  double face_bary_x = origin_x + dx*(i_ft + offset_x + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 0, phase, i_ns,j_ns,k_ns, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                  double face_bary_y = origin_y + dy*(j_ft + offset_y + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 1, phase, i_ns,j_ns,k_ns, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                  double face_bary_z = origin_z + dz*(k_ft + offset_z + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 2, phase, i_ns,j_ns,k_ns, old_indicatrice_surfacique, next_indicatrice_surfacique)));

                                  double distance_facette = sqrt((face_bary_x - coord_facettes[0])*(face_bary_x - coord_facettes[0]) + (face_bary_y - coord_facettes[1])*(face_bary_y - coord_facettes[1]) + (face_bary_z - coord_facettes[2])*(face_bary_z - coord_facettes[2]));
                                  // Ce code est un debut d'implementation. Le calcul de la distance a la facette n'est pas termine (distance au centre uniquement).
                                  Cerr << "distance_facette : calcul pas implementee." << finl;
                                  Process::exit();
                                  facette_la_plus_proche[phase][dir] = distance_facette < min_distance_facette[phase][dir] ? fa7 : facette_la_plus_proche[phase][dir];
                                  min_distance_facette[phase][dir] = distance_facette < min_distance_facette[phase][dir] ? distance_facette : min_distance_facette[phase][dir];
                                }
                            }

                          index = data.index_facette_suivante_;
                        };

                      for (int dir = 0; dir < 3; dir++)
                        {
                          for (int phase = 0; phase < 2; phase++)
                            {
                              double coord_facettes[3] = {0};
                              for (int dir2 = 0; dir2 < 3; dir2++)
                                for (int som = 0; som < 3; som++)
                                  coord_facettes[dir2] += sommets(facettes(facette_la_plus_proche[phase][dir2], som), dir2);
                              for (int dir2 = 0; dir2 < 3; dir2++)
                                coord_facettes[dir2] /= 3.;


                              double old_indicatrice_surfacique  = old_indicatrice_surfacique_face[dir](i_ns,j_ns,k_ns);
                              double next_indicatrice_surfacique = next_indicatrice_surfacique_face[dir](i_ns,j_ns,k_ns);
                              double face_bary_x = origin_x + dx*(i_ft + offset_x + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 0, phase, i_ns,j_ns,k_ns, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                              double face_bary_y = origin_y + dy*(j_ft + offset_y + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 1, phase, i_ns,j_ns,k_ns, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                              double face_bary_z = origin_z + dz*(k_ft + offset_z + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 2, phase, i_ns,j_ns,k_ns, old_indicatrice_surfacique, next_indicatrice_surfacique)));

                              const DoubleTabFT_cut_cell_vector3& indicatrice_surfacique = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face();
                              double indicatrice_surface = ((phase == 0) ? 1 - indicatrice_surfacique(n,dir) : indicatrice_surfacique(n,dir));

                              double lambda = (phase == 0) ? lambda_vapour : lambda_liquid;
                              double dTdn = interfacial_phin_ai(facette_la_plus_proche[phase][dir])/(lambda * surface_facettes(facette_la_plus_proche[phase][dir]));

                              double normal_x = normale_facettes(facette_la_plus_proche[phase][dir], 0);
                              double normal_y = normale_facettes(facette_la_plus_proche[phase][dir], 1);
                              double normal_z = normale_facettes(facette_la_plus_proche[phase][dir], 2);
                              double norm_normal = sqrt(normal_x*normal_x + normal_y*normal_y + normal_z*normal_z);
                              normal_x /= norm_normal;
                              normal_y /= norm_normal;
                              normal_z /= norm_normal;
                              temperature_face_ft[phase][dir](i_ft,j_ft,k_ft) = (indicatrice_surface == 0.) ? 0. : (interfacial_temperature(facette_la_plus_proche[phase][dir])/surface_facettes(facette_la_plus_proche[phase][dir]) + dTdn * ((face_bary_x - coord_facettes[0])*normal_x + (face_bary_y - coord_facettes[1])*normal_y + (face_bary_z - coord_facettes[2])*normal_z));
                            }
                        }
                    }
                }
            }
        }
    }
  }

  for (int dir = 0; dir < 3; dir++)
    {
      for (int phase = 0; phase < 2; phase++)
        {
          temperature_face_ft[phase][dir].echange_espace_virtuel(temperature_face_ft[phase][dir].ghost());

          // Calcul du flux interface sur le domaine NS :
          ref_ijk_ft->redistrib_from_ft_elem().redistribute(temperature_face_ft[phase][dir], temperature_face_ns[phase][dir]);
          temperature_face_ns[phase][dir].echange_espace_virtuel(temperature_face_ns[phase][dir].ghost());
        }
    }
}

void Cut_cell_convection_auxiliaire::calcule_temperature_face_depuis_facette_interpolate(CUT_CELL_CONV_FACE_INTERPOLATION face_interp, double timestep, const ArrOfDouble& interfacial_temperature, const IJK_Field_double& temperature_ft, const Cut_field_scalar& cut_field_temperature, REF(IJK_FT_cut_cell)& ref_ijk_ft, FixedVector<FixedVector<IJK_Field_double, 3>, 2>& temperature_face_ft, FixedVector<FixedVector<IJK_Field_double, 3>, 2>& temperature_face_ns, const Cut_field_vector& cut_field_total_velocity)
{
  const bool next_time = false;

  int count_status_not_found = 0;
  int count_status_below_10 = 0;
  int count_status_below_100 = 0;
  int count_status_below_1000 = 0;
  int count_status_below_10000 = 0;
  int count_status_below_100000 = 0;
  int count_status_below_1000000 = 0;
  int count_status_above_1000000 = 0;

  // Calcul de T_face sur le maillage FT
  {
    const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();
    const IJK_Splitting& s = temperature_face_ft.get_splitting();

    double dx = s.get_grid_geometry().get_constant_delta(DIRECTION_I);
    double dy = s.get_grid_geometry().get_constant_delta(DIRECTION_J);
    double dz = s.get_grid_geometry().get_constant_delta(DIRECTION_K);

    const IJK_Grid_Geometry& geom = s.get_grid_geometry();
    double origin_x = geom.get_origin(DIRECTION_I);
    double origin_y = geom.get_origin(DIRECTION_J);
    double origin_z = geom.get_origin(DIRECTION_K);

    int offset_x = s.get_offset_local(DIRECTION_I);
    int offset_y = s.get_offset_local(DIRECTION_J);
    int offset_z = s.get_offset_local(DIRECTION_K);

    const FixedVector<IJK_Field_double, 3>& old_indicatrice_surfacique_face  = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_face_old();
    const FixedVector<IJK_Field_double, 3>& next_indicatrice_surfacique_face = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_face_next();

    const int ni_ft = temperature_face_ft[0][0].ni();
    const int nj_ft = temperature_face_ft[0][0].nj();
    const int nk_ft = temperature_face_ft[0][0].nk();

    // Initialisation
    {
      for (int k_ft = 0; k_ft < nk_ft; k_ft++)
        {
          for (int j_ft = 0; j_ft < nj_ft; j_ft++)
            {
              for (int i_ft = 0; i_ft < ni_ft; i_ft++)
                {
                  for (int dir = 0; dir < 3; dir++)
                    {
                      for (int phase = 0; phase < 2; phase++)
                        {
                          temperature_face_ft[phase][dir](i_ft, j_ft, k_ft) = 0.;
                        }
                    }
                }
            }
        }
    }

    // Calcul pour les faces coupees par l'interface
    {
      for (int k_ft = 0; k_ft < nk_ft; k_ft++)
        {
          for (int j_ft = 0; j_ft < nj_ft; j_ft++)
            {
              for (int i_ft = 0; i_ft < ni_ft; i_ft++)
                {
                  if (next_time ? (!cut_cell_disc.get_interfaces().est_pure(cut_cell_disc.get_interfaces().In_ft(i_ft,j_ft,k_ft))) : (!cut_cell_disc.get_interfaces().est_pure(cut_cell_disc.get_interfaces().I_ft(i_ft,j_ft,k_ft))))
                    {
                      double x_centre = (i_ft + offset_x + .5)*dx + origin_x;
                      double y_centre = (j_ft + offset_y + .5)*dy + origin_y;
                      double z_centre = (k_ft + offset_z + .5)*dz + origin_z;
                      int i_ns = cut_cell_disc.get_i_selon_dir(0, x_centre);
                      int j_ns = cut_cell_disc.get_i_selon_dir(1, y_centre);
                      int k_ns = cut_cell_disc.get_i_selon_dir(2, z_centre);

                      if (!cut_cell_disc.within_ghost(i_ns, j_ns, k_ns, 1, 1))
                        continue;

                      int n_centre = cut_cell_disc.get_n(i_ns, j_ns, k_ns);

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

                          int n = cut_cell_disc.get_n(i_ns+di, j_ns+dj, k_ns+dk);
                          int n_decale = cut_cell_disc.get_n(i_ns+di_decale, j_ns+dj_decale, k_ns+dk_decale);

                          double normal_x = 0.;
                          double normal_y = 0.;
                          double normal_z = 0.;
                          if (n_centre >= 0 && n_decale >= 0)
                            {
                              double normal_x_centre = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_centre,0);
                              double normal_y_centre = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_centre,1);
                              double normal_z_centre = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_centre,2);
                              double normal_x_decale = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_decale,0);
                              double normal_y_decale = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_decale,1);
                              double normal_z_decale = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_decale,2);
                              normal_x = .5*(normal_x_centre + normal_x_decale);
                              normal_y = .5*(normal_y_centre + normal_y_decale);
                              normal_z = .5*(normal_z_centre + normal_z_decale);
                            }
                          else if (n_centre >= 0)
                            {
                              double normal_x_centre = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_centre,0);
                              double normal_y_centre = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_centre,1);
                              double normal_z_centre = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_centre,2);
                              normal_x = .5*(normal_x_centre);
                              normal_y = .5*(normal_y_centre);
                              normal_z = .5*(normal_z_centre);
                            }
                          else if (n_decale >= 0)
                            {
                              double normal_x_decale = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_decale,0);
                              double normal_y_decale = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_decale,1);
                              double normal_z_decale = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_decale,2);
                              normal_x = .5*(normal_x_decale);
                              normal_y = .5*(normal_y_decale);
                              normal_z = .5*(normal_z_decale);
                            }
                          else
                            {
                              Cerr << "Cut_cell_convection_auxiliaire::calcule_temperature_face_depuis_facette_interpolate: Les deux n sont pures ?" << finl;
                              Process::exit();
                            }
                          double norm_normal = sqrt(normal_x*normal_x + normal_y*normal_y + normal_z*normal_z);
                          normal_x /= norm_normal;
                          normal_y /= norm_normal;
                          normal_z /= norm_normal;


                          for (int phase = 0; phase < 2; phase++)
                            {
                              //int n = cut_cell_disc.get_n(i_ns+di, j_ns+dj, k_ns+dk);
                              //const DoubleTabFT_cut_cell_vector3& indicatrice_surfacique = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face();
                              //double indicatrice_surface = ((phase == 0) ? 1 - indicatrice_surfacique(n,dir) : indicatrice_surfacique(n,dir));
                              double old_indicatrice_surfacique  = old_indicatrice_surfacique_face[dir](i_ns+di,j_ns+dj,k_ns+dk);
                              double next_indicatrice_surfacique = next_indicatrice_surfacique_face[dir](i_ns+di,j_ns+dj,k_ns+dk);
                              const double indicatrice_surfacique = next_time ? next_indicatrice_surfacique : old_indicatrice_surfacique;
                              double indicatrice_surface = ((phase == 0) ? 1 - indicatrice_surfacique : indicatrice_surfacique);
                              if (indicatrice_surface == 0.)
                                {
                                  continue; // Pas utile de determiner une temperature
                                }
                              else
                                {
                                  double coordinates[3] = {};

                                  if (face_interp == CUT_CELL_CONV_FACE_INTERPOLATION::INTERPOLATE)
                                    {
                                      double face_bary_x = origin_x + dx*(i_ft + di + offset_x + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 0, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_y = origin_y + dy*(j_ft + dj + offset_y + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 1, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_z = origin_z + dz*(k_ft + dk + offset_z + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 2, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));

                                      coordinates[0] = face_bary_x;
                                      coordinates[1] = face_bary_y;
                                      coordinates[2] = face_bary_z;
                                    }
                                  else if (face_interp == CUT_CELL_CONV_FACE_INTERPOLATION::INTERPOLATEAMONT0)
                                    {
                                      double v_x = (n < 0) ? cut_field_total_velocity.pure_[0](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,0) : cut_field_total_velocity.diph_l_(n,0));
                                      double v_y = (n < 0) ? cut_field_total_velocity.pure_[1](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,1) : cut_field_total_velocity.diph_l_(n,1));
                                      double v_z = (n < 0) ? cut_field_total_velocity.pure_[2](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,2) : cut_field_total_velocity.diph_l_(n,2));
                                      double v_dir = select(dir, v_x, v_y, v_z);

                                      int i_ft_amont = i_ft + di - (dir == 0)*(v_dir>0);
                                      int j_ft_amont = j_ft + dj - (dir == 1)*(v_dir>0);
                                      int k_ft_amont = k_ft + dk - (dir == 2)*(v_dir>0);

                                      int i_ns_amont = i_ns + di - (dir == 0)*(v_dir>0);
                                      int j_ns_amont = j_ns + dj - (dir == 1)*(v_dir>0);
                                      int k_ns_amont = k_ns + dk - (dir == 2)*(v_dir>0);

                                      double old_indicatrice_amont  = cut_cell_disc.get_interfaces().I(i_ns_amont,j_ns_amont,k_ns_amont);
                                      double next_indicatrice_amont = cut_cell_disc.get_interfaces().In(i_ns_amont,j_ns_amont,k_ns_amont);

                                      double amont_bary_x = origin_x + dx*(i_ft_amont + offset_x + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 0, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_y = origin_y + dy*(j_ft_amont + offset_y + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 1, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_z = origin_z + dz*(k_ft_amont + offset_z + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 2, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));

                                      coordinates[0] = amont_bary_x;
                                      coordinates[1] = amont_bary_y;
                                      coordinates[2] = amont_bary_z;
                                    }
                                  else if (face_interp == CUT_CELL_CONV_FACE_INTERPOLATION::INTERPOLATEAMONT1)
                                    {
                                      double v_x = (n < 0) ? cut_field_total_velocity.pure_[0](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,0) : cut_field_total_velocity.diph_l_(n,0));
                                      double v_y = (n < 0) ? cut_field_total_velocity.pure_[1](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,1) : cut_field_total_velocity.diph_l_(n,1));
                                      double v_z = (n < 0) ? cut_field_total_velocity.pure_[2](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,2) : cut_field_total_velocity.diph_l_(n,2));
                                      double v_dir = select(dir, v_x, v_y, v_z);

                                      int i_ft_amont = i_ft + di - (dir == 0)*(v_dir>0);
                                      int j_ft_amont = j_ft + dj - (dir == 1)*(v_dir>0);
                                      int k_ft_amont = k_ft + dk - (dir == 2)*(v_dir>0);

                                      int i_ns_amont = i_ns + di - (dir == 0)*(v_dir>0);
                                      int j_ns_amont = j_ns + dj - (dir == 1)*(v_dir>0);
                                      int k_ns_amont = k_ns + dk - (dir == 2)*(v_dir>0);

                                      double old_indicatrice_amont  = cut_cell_disc.get_interfaces().I(i_ns_amont,j_ns_amont,k_ns_amont);
                                      double next_indicatrice_amont = cut_cell_disc.get_interfaces().In(i_ns_amont,j_ns_amont,k_ns_amont);

                                      double amont_bary_x = origin_x + dx*(i_ft_amont + offset_x + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 0, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_y = origin_y + dy*(j_ft_amont + offset_y + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 1, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_z = origin_z + dz*(k_ft_amont + offset_z + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 2, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));

                                      double face_bary_x = origin_x + dx*(i_ft + di + offset_x + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 0, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_y = origin_y + dy*(j_ft + dj + offset_y + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 1, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_z = origin_z + dz*(k_ft + dk + offset_z + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 2, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));

                                      coordinates[0] = (dir==0) ? amont_bary_x : face_bary_x;
                                      coordinates[1] = (dir==1) ? amont_bary_y : face_bary_y;
                                      coordinates[2] = (dir==2) ? amont_bary_z : face_bary_z;
                                    }
                                  else if (face_interp == CUT_CELL_CONV_FACE_INTERPOLATION::INTERPOLATEAMONT2)
                                    {
                                      double v_x = (n < 0) ? cut_field_total_velocity.pure_[0](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,0) : cut_field_total_velocity.diph_l_(n,0));
                                      double v_y = (n < 0) ? cut_field_total_velocity.pure_[1](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,1) : cut_field_total_velocity.diph_l_(n,1));
                                      double v_z = (n < 0) ? cut_field_total_velocity.pure_[2](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,2) : cut_field_total_velocity.diph_l_(n,2));

                                      double face_bary_x = origin_x + dx*(i_ft + di + offset_x + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 0, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_y = origin_y + dy*(j_ft + dj + offset_y + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 1, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_z = origin_z + dz*(k_ft + dk + offset_z + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 2, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));

                                      coordinates[0] = face_bary_x - .5*v_x*timestep;
                                      coordinates[1] = face_bary_y - .5*v_y*timestep;
                                      coordinates[2] = face_bary_z - .5*v_z*timestep;
                                    }
                                  else if (face_interp == CUT_CELL_CONV_FACE_INTERPOLATION::INTERPOLATENORMALE0)
                                    {
                                      double v_x = (n < 0) ? cut_field_total_velocity.pure_[0](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,0) : cut_field_total_velocity.diph_l_(n,0));
                                      double v_y = (n < 0) ? cut_field_total_velocity.pure_[1](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,1) : cut_field_total_velocity.diph_l_(n,1));
                                      double v_z = (n < 0) ? cut_field_total_velocity.pure_[2](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,2) : cut_field_total_velocity.diph_l_(n,2));
                                      double v_dir = select(dir, v_x, v_y, v_z);

                                      int i_ft_amont = i_ft + di - (dir == 0)*(v_dir>0);
                                      int j_ft_amont = j_ft + dj - (dir == 1)*(v_dir>0);
                                      int k_ft_amont = k_ft + dk - (dir == 2)*(v_dir>0);

                                      int i_ns_amont = i_ns + di - (dir == 0)*(v_dir>0);
                                      int j_ns_amont = j_ns + dj - (dir == 1)*(v_dir>0);
                                      int k_ns_amont = k_ns + dk - (dir == 2)*(v_dir>0);

                                      double old_indicatrice_amont  = cut_cell_disc.get_interfaces().I(i_ns_amont,j_ns_amont,k_ns_amont);
                                      double next_indicatrice_amont = cut_cell_disc.get_interfaces().In(i_ns_amont,j_ns_amont,k_ns_amont);

                                      double amont_bary_x = origin_x + dx*(i_ft_amont + offset_x + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 0, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_y = origin_y + dy*(j_ft_amont + offset_y + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 1, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_z = origin_z + dz*(k_ft_amont + offset_z + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 2, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));

                                      double face_bary_x = origin_x + dx*(i_ft + di + offset_x + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 0, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_y = origin_y + dy*(j_ft + dj + offset_y + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 1, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_z = origin_z + dz*(k_ft + dk + offset_z + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 2, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));

                                      double normal_dir = select(dir, normal_x, normal_y, normal_z);
                                      bool case_perpendicular = (std::abs(normal_dir) > 0.72);

                                      coordinates[0] = case_perpendicular ? amont_bary_x : face_bary_x;
                                      coordinates[1] = case_perpendicular ? amont_bary_y : face_bary_y;
                                      coordinates[2] = case_perpendicular ? amont_bary_z : face_bary_z;
                                    }
                                  else if (face_interp == CUT_CELL_CONV_FACE_INTERPOLATION::INTERPOLATENORMALE1)
                                    {
                                      double v_x = (n < 0) ? cut_field_total_velocity.pure_[0](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,0) : cut_field_total_velocity.diph_l_(n,0));
                                      double v_y = (n < 0) ? cut_field_total_velocity.pure_[1](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,1) : cut_field_total_velocity.diph_l_(n,1));
                                      double v_z = (n < 0) ? cut_field_total_velocity.pure_[2](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,2) : cut_field_total_velocity.diph_l_(n,2));
                                      double v_dir = select(dir, v_x, v_y, v_z);

                                      int i_ft_amont = i_ft + di - (dir == 0)*(v_dir>0);
                                      int j_ft_amont = j_ft + dj - (dir == 1)*(v_dir>0);
                                      int k_ft_amont = k_ft + dk - (dir == 2)*(v_dir>0);

                                      int i_ns_amont = i_ns + di - (dir == 0)*(v_dir>0);
                                      int j_ns_amont = j_ns + dj - (dir == 1)*(v_dir>0);
                                      int k_ns_amont = k_ns + dk - (dir == 2)*(v_dir>0);

                                      double old_indicatrice_amont  = cut_cell_disc.get_interfaces().I(i_ns_amont,j_ns_amont,k_ns_amont);
                                      double next_indicatrice_amont = cut_cell_disc.get_interfaces().In(i_ns_amont,j_ns_amont,k_ns_amont);

                                      double amont_bary_x = origin_x + dx*(i_ft_amont + offset_x + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 0, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_y = origin_y + dy*(j_ft_amont + offset_y + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 1, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_z = origin_z + dz*(k_ft_amont + offset_z + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 2, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));

                                      double face_bary_x = origin_x + dx*(i_ft + di + offset_x + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 0, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_y = origin_y + dy*(j_ft + dj + offset_y + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 1, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_z = origin_z + dz*(k_ft + dk + offset_z + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 2, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));

                                      double normal_dir = select(dir, normal_x, normal_y, normal_z);
                                      bool case_perpendicular = (std::abs(normal_dir) > 0.72);

                                      coordinates[0] = case_perpendicular ? face_bary_x : amont_bary_x;
                                      coordinates[1] = case_perpendicular ? face_bary_y : amont_bary_y;
                                      coordinates[2] = case_perpendicular ? face_bary_z : amont_bary_z;
                                    }
                                  else if (face_interp == CUT_CELL_CONV_FACE_INTERPOLATION::INTERPOLATEDISTANCE0)
                                    {
                                      double v_x = (n < 0) ? cut_field_total_velocity.pure_[0](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,0) : cut_field_total_velocity.diph_l_(n,0));
                                      double v_y = (n < 0) ? cut_field_total_velocity.pure_[1](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,1) : cut_field_total_velocity.diph_l_(n,1));
                                      double v_z = (n < 0) ? cut_field_total_velocity.pure_[2](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,2) : cut_field_total_velocity.diph_l_(n,2));
                                      double v_dir = select(dir, v_x, v_y, v_z);

                                      int i_ft_amont = i_ft + di - (dir == 0)*(v_dir>0);
                                      int j_ft_amont = j_ft + dj - (dir == 1)*(v_dir>0);
                                      int k_ft_amont = k_ft + dk - (dir == 2)*(v_dir>0);

                                      int i_ns_amont = i_ns + di - (dir == 0)*(v_dir>0);
                                      int j_ns_amont = j_ns + dj - (dir == 1)*(v_dir>0);
                                      int k_ns_amont = k_ns + dk - (dir == 2)*(v_dir>0);

                                      double old_indicatrice_amont  = cut_cell_disc.get_interfaces().I(i_ns_amont,j_ns_amont,k_ns_amont);
                                      double next_indicatrice_amont = cut_cell_disc.get_interfaces().In(i_ns_amont,j_ns_amont,k_ns_amont);

                                      double amont_bary_x = origin_x + dx*(i_ft_amont + offset_x + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 0, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_y = origin_y + dy*(j_ft_amont + offset_y + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 1, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_z = origin_z + dz*(k_ft_amont + offset_z + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 2, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));

                                      double face_bary_x = origin_x + dx*(i_ft + di + offset_x + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 0, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_y = origin_y + dy*(j_ft + dj + offset_y + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 1, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_z = origin_z + dz*(k_ft + dk + offset_z + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 2, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));

                                      double dist_x = (amont_bary_x - face_bary_x);
                                      double dist_y = (amont_bary_y - face_bary_y);
                                      double dist_z = (amont_bary_z - face_bary_z);
                                      double dist_face_amont = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
                                      bool close = (dist_face_amont < 0.1);

                                      coordinates[0] = close ? amont_bary_x : face_bary_x;
                                      coordinates[1] = close ? amont_bary_y : face_bary_y;
                                      coordinates[2] = close ? amont_bary_z : face_bary_z;
                                    }
                                  else if (face_interp == CUT_CELL_CONV_FACE_INTERPOLATION::INTERPOLATEDISTANCE1)
                                    {
                                      double v_x = (n < 0) ? cut_field_total_velocity.pure_[0](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,0) : cut_field_total_velocity.diph_l_(n,0));
                                      double v_y = (n < 0) ? cut_field_total_velocity.pure_[1](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,1) : cut_field_total_velocity.diph_l_(n,1));
                                      double v_z = (n < 0) ? cut_field_total_velocity.pure_[2](i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity.diph_v_(n,2) : cut_field_total_velocity.diph_l_(n,2));
                                      double v_dir = select(dir, v_x, v_y, v_z);

                                      int i_ft_amont = i_ft + di - (dir == 0)*(v_dir>0);
                                      int j_ft_amont = j_ft + dj - (dir == 1)*(v_dir>0);
                                      int k_ft_amont = k_ft + dk - (dir == 2)*(v_dir>0);

                                      int i_ns_amont = i_ns + di - (dir == 0)*(v_dir>0);
                                      int j_ns_amont = j_ns + dj - (dir == 1)*(v_dir>0);
                                      int k_ns_amont = k_ns + dk - (dir == 2)*(v_dir>0);

                                      double old_indicatrice_amont  = cut_cell_disc.get_interfaces().I(i_ns_amont,j_ns_amont,k_ns_amont);
                                      double next_indicatrice_amont = cut_cell_disc.get_interfaces().In(i_ns_amont,j_ns_amont,k_ns_amont);

                                      double amont_bary_x = origin_x + dx*(i_ft_amont + offset_x + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 0, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_y = origin_y + dy*(j_ft_amont + offset_y + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 1, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_z = origin_z + dz*(k_ft_amont + offset_z + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 2, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));

                                      double face_bary_x = origin_x + dx*(i_ft + di + offset_x + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 0, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_y = origin_y + dy*(j_ft + dj + offset_y + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 1, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_z = origin_z + dz*(k_ft + dk + offset_z + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 2, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));

                                      double dist_x = (amont_bary_x - face_bary_x);
                                      double dist_y = (amont_bary_y - face_bary_y);
                                      double dist_z = (amont_bary_z - face_bary_z);
                                      double dist_face_amont = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
                                      bool close = (dist_face_amont < 0.1);

                                      coordinates[0] = close ? face_bary_x : amont_bary_x;
                                      coordinates[1] = close ? face_bary_y : amont_bary_y;
                                      coordinates[2] = close ? face_bary_z : amont_bary_z;
                                    }
                                  else
                                    {
                                      Cerr << "Cut_cell_convection_auxiliaire::calcule_temperature_face_depuis_facette_interpolate: face_interp inconnue." << finl;
                                      Process::exit();
                                    }

                                  int status = -2;
                                  const int tolerate_not_within_tetrahedron = 2;
                                  double temperature_interpolate = ijk_interpolate_cut_cell_using_interface(false, phase, temperature_ft, cut_field_temperature, interfacial_temperature, coordinates, tolerate_not_within_tetrahedron, status);

                                  if (status == -1)
                                    {
                                      temperature_interpolate = 0.;
                                    }

                                  // Note : Quand le statut est -1, la valeur de la temperature est souvent absolument pas pertinente
                                  temperature_face_ft[phase][dir](i_ft+di,j_ft+dj,k_ft+dk) = (status == -1 || indicatrice_surface == 0.) ? 0. : temperature_interpolate;

                                  assert(status > -2);
                                  if (status == -1)
                                    {
                                      count_status_not_found++;
                                    }
                                  else if (status < 10)
                                    {
                                      count_status_below_10++;
                                    }
                                  else if (status < 100)
                                    {
                                      count_status_below_100++;
                                    }
                                  else if (status < 1000)
                                    {
                                      count_status_below_1000++;
                                    }
                                  else if (status < 10000)
                                    {
                                      count_status_below_10000++;
                                    }
                                  else if (status < 100000)
                                    {
                                      count_status_below_100000++;
                                    }
                                  else if (status < 1000000)
                                    {
                                      count_status_below_1000000++;
                                    }
                                  else if (status >= 1000000)
                                    {
                                      count_status_above_1000000++;
                                    }
                                  else
                                    {
                                      Cerr << "Un tel statut n'est pas possible." << finl;
                                      Process::exit();
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
  }

  for (int dir = 0; dir < 3; dir++)
    {
      for (int phase = 0; phase < 2; phase++)
        {
          temperature_face_ft[phase][dir].echange_espace_virtuel(temperature_face_ft[phase][dir].ghost());

          // Calcul du flux interface sur le domaine NS :
          ref_ijk_ft->redistrib_from_ft_elem().redistribute(temperature_face_ft[phase][dir], temperature_face_ns[phase][dir]);
          temperature_face_ns[phase][dir].echange_espace_virtuel(temperature_face_ns[phase][dir].ghost());
        }
    }

  Cerr << "Bilan des statuts pour Cut_cell_convection_auxiliaire::calcule_temperature_face_depuis_facette_interpolate : " << finl;
  Cerr << "    -1      " << count_status_not_found  << finl;
  Cerr << "   <10      " << count_status_below_10   << finl;
  Cerr << "  <100      " << count_status_below_100  << finl;
  Cerr << " <1000      " << count_status_below_1000 << finl;
  Cerr << " <10000     " << count_status_below_10000 << finl;
  Cerr << " <100000    " << count_status_below_100000 << finl;
  Cerr << " <1000000   " << count_status_below_1000000 << finl;
  Cerr << ">=1000000   " << count_status_above_1000000 << finl;
}
