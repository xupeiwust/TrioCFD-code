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
// File      : Cut_cell_schema_auxiliaire.cpp
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#include <Cut_cell_schema_auxiliaire.h>
#include <Cut_cell_FT_Disc.h>
#include <IJK_Thermal_base.h>
#include <Process.h>
#include <Probleme_FTD_IJK_cut_cell.h>
#include <stat_counters.h>
#include <IJK_Navier_Stokes_tools.h>
#include <IJK_Navier_Stokes_tools_cut_cell.h>
#include <Param.h>

Implemente_base_sans_constructeur(Cut_cell_schema_auxiliaire, "Cut_cell_schema_auxiliaire", Objet_U) ;

Cut_cell_schema_auxiliaire::Cut_cell_schema_auxiliaire()
{
  tolerate_not_within_tetrahedron_ = -1;
}

Sortie& Cut_cell_schema_auxiliaire::printOn(Sortie& os) const
{
  return os;
}

Entree& Cut_cell_schema_auxiliaire::readOn(Entree& is)
{
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades(is);
  return is;
}

void Cut_cell_schema_auxiliaire::set_param(Param& param)
{
  param.ajouter("methode_valeur_remplissage", (int*)&methode_valeur_remplissage_);
  param.dictionnaire("non_initialise",(int)METHODE_TEMPERATURE_REMPLISSAGE::NON_INITIALISE);
  param.dictionnaire("copie_directe", (int)METHODE_TEMPERATURE_REMPLISSAGE::COPIE_DIRECTE);
  param.dictionnaire("ponderation_voisin", (int)METHODE_TEMPERATURE_REMPLISSAGE::PONDERATION_VOISIN);
  param.dictionnaire("ponderation_directionnelle_voisin", (int)METHODE_TEMPERATURE_REMPLISSAGE::PONDERATION_DIRECTIONNELLE_VOISIN);
  param.dictionnaire("semi_lagrangien", (int)METHODE_TEMPERATURE_REMPLISSAGE::SEMI_LAGRANGIEN);
  param.dictionnaire("semi_lagrangien_interpolate", (int)METHODE_TEMPERATURE_REMPLISSAGE::SEMI_LAGRANGIEN_INTERPOLATE);

  param.ajouter("correction_petites_cellules", (int*)&correction_petites_cellules_);
  param.dictionnaire("correction_directe", (int)CORRECTION_PETITES_CELLULES::CORRECTION_DIRECTE);
  param.dictionnaire("direction_privilegiee", (int)CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE);
  param.dictionnaire("direction_privilegiee_2", (int)CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE_2);
  param.dictionnaire("correction_symetrique", (int)CORRECTION_PETITES_CELLULES::CORRECTION_SYMETRIQUE);
  param.dictionnaire("direction_privilegiee_avec_limitation", (int)CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE_AVEC_LIMITATION);
  param.dictionnaire("direction_privilegiee_avec_limitation_2", (int)CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE_AVEC_LIMITATION_2);
  param.dictionnaire("correction_symetrique_avec_limitation", (int)CORRECTION_PETITES_CELLULES::CORRECTION_SYMETRIQUE_AVEC_LIMITATION);

  param.ajouter("tolerate_not_within_tetrahedron", (int*)&tolerate_not_within_tetrahedron_);
}

void Cut_cell_schema_auxiliaire::initialise(Cut_cell_FT_Disc& cut_cell_disc)
{
  temperature_remplissage_.associer_ephemere(cut_cell_disc);
  flux_naive_.associer_ephemere(cut_cell_disc);
}

void Cut_cell_schema_auxiliaire::calcule_valeur_remplissage(double timestep, double lambda_liquid, double lambda_vapour, const IJK_Field_double& flux_interface_ns, const ArrOfDouble& interfacial_temperature, const IJK_Field_double& temperature_ft, const Cut_field_vector3_double& cut_field_total_velocity, const Cut_field_double& cut_field_temperature)
{
  if (methode_valeur_remplissage_ == METHODE_TEMPERATURE_REMPLISSAGE::COPIE_DIRECTE)
    {
      calcule_valeur_remplissage_copie_directe(cut_field_temperature, temperature_remplissage_);
    }
  else if (methode_valeur_remplissage_ == METHODE_TEMPERATURE_REMPLISSAGE::PONDERATION_VOISIN)
    {
      calcule_valeur_remplissage_ponderation_voisin(false, cut_field_total_velocity, cut_field_temperature, temperature_remplissage_);
    }
  else if (methode_valeur_remplissage_ == METHODE_TEMPERATURE_REMPLISSAGE::PONDERATION_DIRECTIONNELLE_VOISIN)
    {
      calcule_valeur_remplissage_ponderation_voisin(true, cut_field_total_velocity, cut_field_temperature, temperature_remplissage_);
    }
  else if (methode_valeur_remplissage_ == METHODE_TEMPERATURE_REMPLISSAGE::SEMI_LAGRANGIEN)
    {
      calcule_valeur_remplissage_semi_lagrangien(timestep, lambda_liquid, lambda_vapour, flux_interface_ns, cut_field_temperature, temperature_remplissage_);
    }
  else if (methode_valeur_remplissage_ == METHODE_TEMPERATURE_REMPLISSAGE::SEMI_LAGRANGIEN_INTERPOLATE)
    {
      calcule_valeur_remplissage_semi_lagrangien_interpolate(timestep, interfacial_temperature, temperature_ft, cut_field_temperature, temperature_remplissage_);
    }
  else
    {
      Cerr << "Methode non reconnue pour le calcul de la valeur de remplissage." << finl;
      Process::exit();
    }
}

void Cut_cell_schema_auxiliaire::compute_flux_dying_cells(const Cut_field_vector3_double& cut_field_total_velocity, const Cut_field_double& cut_field_temperature)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();

  const IJK_Grid_Geometry& geom = cut_cell_disc.get_splitting().get_grid_geometry();
  const double delta_x = geom.get_constant_delta(0);
  const double delta_y = geom.get_constant_delta(1);
  const double delta_z = geom.get_constant_delta(2);
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(0));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(1));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(2));

  int statut_diphasique = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::MOURRANT);
  int index_min = cut_cell_disc.get_statut_diphasique_value_index(statut_diphasique);
  int index_max = cut_cell_disc.get_statut_diphasique_value_index(statut_diphasique+1);
  for (int index = index_min; index < index_max; index++)
    {
      int n = cut_cell_disc.get_n_from_statut_diphasique_index(index);

      Int3 ijk = cut_cell_disc.get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      if (!cut_cell_disc.get_splitting().within_ghost(i, j, k, 1, 1))
        continue;

      double old_indicatrice = cut_cell_disc.get_interfaces().I(i,j,k);
      double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);
      assert(cut_cell_disc.get_interfaces().est_pure(next_indicatrice));
      int phase = 1 - IJK_Interfaces::convert_indicatrice_to_phase(next_indicatrice); // phase de la cellule mourrante



      double flux[6] = {0};
      for (int num_face = 0; num_face < 6; num_face++)
        {
          int dir = num_face%3;
          int decalage = num_face/3;
          int sign = decalage*2 -1;

          int di_decale = sign*(dir == 0);
          int dj_decale = sign*(dir == 1);
          int dk_decale = sign*(dir == 2);

          double f_dir = select_dir(dir, delta_y*delta_z, delta_x*delta_z, delta_x*delta_y);

          double old_indicatrice_decale = cut_cell_disc.get_interfaces().I(i+di_decale,j+dj_decale,k+dk_decale);
          bool decale_also_dying = (cut_cell_disc.get_interfaces().phase_mourrante(phase, i+di_decale,j+dj_decale,k+dk_decale));
          bool decale_also_nascent = (cut_cell_disc.get_interfaces().phase_naissante(phase, i+di_decale,j+dj_decale,k+dk_decale));
          bool decale_smaller = (phase == 0) ? (1 - old_indicatrice_decale) < (1 - old_indicatrice) : (old_indicatrice_decale) < (old_indicatrice);
          if (decale_also_nascent || (decale_also_dying && decale_smaller))
            {
              flux[num_face] = 0;
            }
          else
            {
              flux[num_face] = f_dir*dying_cells_flux(num_face, phase, n, cut_field_total_velocity, cut_field_temperature);
            }

          flux_naive_(n,num_face) = flux[num_face];
        }

      // Correction du flux si aucun flux non-nul
      if ((std::abs(flux[0]) + std::abs(flux[1]) + std::abs(flux[2]) + std::abs(flux[3]) + std::abs(flux[4]) + std::abs(flux[5])) == 0)
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

              double f_dir = select_dir(dir, delta_y*delta_z, delta_x*delta_z, delta_x*delta_y);

              double old_indicatrice_decale = cut_cell_disc.get_interfaces().I(i+di_decale,j+dj_decale,k+dk_decale);
              bool decale_also_dying = (cut_cell_disc.get_interfaces().phase_mourrante(phase, i+di_decale,j+dj_decale,k+dk_decale));
              bool decale_also_nascent = (cut_cell_disc.get_interfaces().phase_naissante(phase, i+di_decale,j+dj_decale,k+dk_decale));
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
                      assert((surface_efficace == 0) || (surface_efficace == 1));
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

void Cut_cell_schema_auxiliaire::compute_flux_small_nascent_cells(const Cut_field_vector3_double& cut_field_total_velocity, Cut_field_double& cut_field_temperature)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();

  const IJK_Grid_Geometry& geom = cut_cell_disc.get_splitting().get_grid_geometry();
  const double delta_x = geom.get_constant_delta(DIRECTION_I);
  const double delta_y = geom.get_constant_delta(DIRECTION_J);
  const double delta_z = geom.get_constant_delta(DIRECTION_K);
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(0));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(1));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(2));

  int statut_diphasique_naissant = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::NAISSANT);
  int statut_diphasique_petit = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::DESEQUILIBRE_FINAL);
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

      if (!cut_cell_disc.get_splitting().within_ghost(i, j, k, 1, 1))
        continue;

      double old_indicatrice = cut_cell_disc.get_interfaces().I(i,j,k);
      double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);
      int est_naissant = cut_cell_disc.get_interfaces().est_pure(old_indicatrice);
      int phase = est_naissant ? 1 - IJK_Interfaces::convert_indicatrice_to_phase(old_indicatrice) : ((cut_cell_disc.get_interfaces().below_small_threshold(next_indicatrice)) ? 1 : 0); // phase de la cellule petite ou naissante
      assert(est_naissant || cut_cell_disc.get_interfaces().next_below_small_threshold_for_phase(phase, cut_cell_disc.get_interfaces().I(i,j,k), next_indicatrice));



      double flux[6] = {0};
      for (int num_face = 0; num_face < 6; num_face++)
        {
          int dir = num_face%3;
          int decalage = num_face/3;
          int sign = decalage*2 -1;

          int di_decale = sign*(dir == 0);
          int dj_decale = sign*(dir == 1);
          int dk_decale = sign*(dir == 2);

          double f_dir = select_dir(dir, delta_y*delta_z, delta_x*delta_z, delta_x*delta_y);

          double old_indicatrice_decale = cut_cell_disc.get_interfaces().I(i+di_decale,j+dj_decale,k+dk_decale);
          double next_indicatrice_decale = cut_cell_disc.get_interfaces().In(i+di_decale,j+dj_decale,k+dk_decale);
          bool decale_also_dying = (cut_cell_disc.get_interfaces().phase_mourrante(phase, i+di_decale,j+dj_decale,k+dk_decale));
          bool decale_nascent = (cut_cell_disc.get_interfaces().phase_naissante(phase, i+di_decale,j+dj_decale,k+dk_decale));
          bool decale_small = cut_cell_disc.get_interfaces().next_below_small_threshold_for_phase(phase, old_indicatrice_decale, next_indicatrice_decale);
          bool decale_smaller = (phase == 0) ? (1 - next_indicatrice_decale) < (1 - next_indicatrice) : (next_indicatrice_decale) < (next_indicatrice);
          if (decale_also_dying || (est_naissant && decale_nascent && decale_smaller) || ((!est_naissant) && decale_nascent) || ((!est_naissant) && decale_small && decale_smaller))
            {
              flux[num_face] = 0;
            }
          else
            {
              flux[num_face] = f_dir*small_nascent_cells_flux(num_face, phase, n, cut_field_total_velocity, cut_field_temperature);
            }

          flux_naive_(n,num_face) = flux[num_face];
        }

      // Correction du flux si aucun flux non-nul
      if ((std::abs(flux[0]) + std::abs(flux[1]) + std::abs(flux[2]) + std::abs(flux[3]) + std::abs(flux[4]) + std::abs(flux[5])) == 0)
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

              double f_dir = select_dir(dir, delta_y*delta_z, delta_x*delta_z, delta_x*delta_y);

              double old_indicatrice_decale = cut_cell_disc.get_interfaces().I(i+di_decale,j+dj_decale,k+dk_decale);
              double next_indicatrice_decale = cut_cell_disc.get_interfaces().In(i+di_decale,j+dj_decale,k+dk_decale);
              bool decale_also_dying = (cut_cell_disc.get_interfaces().phase_mourrante(phase, i+di_decale,j+dj_decale,k+dk_decale));
              bool decale_nascent = (cut_cell_disc.get_interfaces().phase_naissante(phase, i+di_decale,j+dj_decale,k+dk_decale));
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
                      assert((surface_efficace == 0) || (surface_efficace == 1));
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

void Cut_cell_schema_auxiliaire::calcule_valeur_remplissage_copie_directe(const Cut_field_double& cut_field_temperature, DoubleTabFT_cut_cell_scalar& valeur_remplissage)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();

  int statut_diphasique_naissant = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::NAISSANT);
  int statut_diphasique_petit = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::DESEQUILIBRE_FINAL);
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

      if (!cut_cell_disc.get_splitting().within_ghost(i, j, k, 1, 1))
        continue;

      double old_indicatrice = cut_cell_disc.get_interfaces().I(i,j,k);
      double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);
      int est_naissant = cut_cell_disc.get_interfaces().est_pure(old_indicatrice);
      int phase = est_naissant ? 1 - IJK_Interfaces::convert_indicatrice_to_phase(old_indicatrice) : ((cut_cell_disc.get_interfaces().below_small_threshold(next_indicatrice)) ? 1 : 0); // phase de la cellule petite ou naissante

      valeur_remplissage(n) = (phase == 0) ? cut_field_temperature.diph_v_(n) : cut_field_temperature.diph_l_(n);
    }
}

void Cut_cell_schema_auxiliaire::calcule_valeur_remplissage_ponderation_voisin(bool est_directionnel, const Cut_field_vector3_double& cut_field_total_velocity, const Cut_field_double& cut_field_temperature, DoubleTabFT_cut_cell_scalar& valeur_remplissage)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();

  const IJK_Grid_Geometry& geom = cut_cell_disc.get_splitting().get_grid_geometry();
  const double delta_x = geom.get_constant_delta(0);
  const double delta_y = geom.get_constant_delta(1);
  const double delta_z = geom.get_constant_delta(2);
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(0));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(1));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(2));

  int statut_diphasique_naissant = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::NAISSANT);
  int statut_diphasique_petit = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::DESEQUILIBRE_FINAL);
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

      if (!cut_cell_disc.get_splitting().within_ghost(i, j, k, 1, 1))
        continue;

      double old_indicatrice = cut_cell_disc.get_interfaces().I(i,j,k);
      double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);
      int est_naissant = cut_cell_disc.get_interfaces().est_pure(old_indicatrice);
      int phase = est_naissant ? 1 - IJK_Interfaces::convert_indicatrice_to_phase(old_indicatrice) : ((cut_cell_disc.get_interfaces().below_small_threshold(next_indicatrice)) ? 1 : 0); // phase de la cellule petite ou naissante
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

          double f_dir = select_dir(dir, delta_y*delta_z, delta_x*delta_z, delta_x*delta_y);

          double next_indicatrice_decale = cut_cell_disc.get_interfaces().In(i+di_decale,j+dj_decale,k+dk_decale);

          double temperature_decale = -1e37;
          int n_decale = cut_cell_disc.get_n(i+di_decale, j+dj_decale, k+dk_decale);
          if (n_decale >= 0)
            {
              temperature_decale = (phase == 0) ? cut_field_temperature.diph_v_(n_decale) : cut_field_temperature.diph_l_(n_decale);
            }
          else
            {
              temperature_decale = (phase == (IJK_Interfaces::convert_indicatrice_to_phase(next_indicatrice_decale)))*cut_field_temperature.pure_(i+di_decale,j+dj_decale,k+dk_decale);
            }

          total_velocity[num_face] = cut_field_total_velocity[dir].from_ijk_and_phase(i+di,j+dj,k+dk, phase);

          int n_face = cut_cell_disc.get_n_face(num_face, n, i, j, k);
          if (n_face >= 0)
            {
              double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir) : cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir);
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
              assert((surface_efficace == 0) || (surface_efficace == 1));
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
      if ((std::abs(flux[0]) + std::abs(flux[1]) + std::abs(flux[2]) + std::abs(flux[3]) + std::abs(flux[4]) + std::abs(flux[5])) == 0)
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

              double f_dir = select_dir(dir, delta_y*delta_z, delta_x*delta_z, delta_x*delta_y);

              double next_indicatrice_decale = cut_cell_disc.get_interfaces().In(i+di_decale,j+dj_decale,k+dk_decale);

              double temperature_decale = -1e37;
              int n_decale = cut_cell_disc.get_n(i+di_decale, j+dj_decale, k+dk_decale);
              if (n_decale >= 0)
                {
                  temperature_decale = (phase == 0) ? cut_field_temperature.diph_v_(n_decale) : cut_field_temperature.diph_l_(n_decale);
                }
              else
                {
                  temperature_decale = (phase == (IJK_Interfaces::convert_indicatrice_to_phase(next_indicatrice_decale)))*cut_field_temperature.pure_(i+di_decale,j+dj_decale,k+dk_decale);
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
                  assert((surface_efficace == 0) || (surface_efficace == 1));
                  double temperature = temperature_decale;
                  flux[num_face] = f_dir*surface_efficace*temperature;
                  temp[num_face] = temperature;
                }
            }
        }

      if ((std::abs(total_velocity[0]) + std::abs(total_velocity[1]) + std::abs(total_velocity[2]) + std::abs(total_velocity[3]) + std::abs(total_velocity[4]) + std::abs(total_velocity[5])) == 0)
        {
          // All velocities are zero
          valeur_remplissage(n) = DMINFLOAT+1; // Dummy value to indicate that no filling should be made.
        }
      if ((std::abs(temp[0]) + std::abs(temp[1]) + std::abs(temp[2]) + std::abs(temp[3]) + std::abs(temp[4]) + std::abs(temp[5])) < 1e-37)
        {
          // It looks like the input field is constant zero. In that case, the filling value should be zero as well.
          valeur_remplissage(n) = 0;
        }
      else
        {
          double somme_T_flux = (temp[0]*std::abs(flux[0]) != 0.)*std::abs(flux[0]) + (temp[1]*std::abs(flux[1]) != 0.)*std::abs(flux[1]) + (temp[2]*std::abs(flux[2]) != 0.)*std::abs(flux[2]) + (temp[3]*std::abs(flux[3]) != 0.)*std::abs(flux[3]) + (temp[4]*std::abs(flux[4]) != 0.)*std::abs(flux[4]) + (temp[5]*std::abs(flux[5]) != 0.)*std::abs(flux[5]);
          assert(somme_T_flux != 0.);
          double val_remplissage = (temp[0]*std::abs(flux[0]) + temp[1]*std::abs(flux[1]) + temp[2]*std::abs(flux[2]) + temp[3]*std::abs(flux[3]) + temp[4]*std::abs(flux[4]) + temp[5]*std::abs(flux[5]))/somme_T_flux;
          valeur_remplissage(n) = val_remplissage;
        }
    }
}

void Cut_cell_schema_auxiliaire::calcule_valeur_remplissage_semi_lagrangien(double timestep, double lambda_liquid, double lambda_vapour, const IJK_Field_double& flux_interface_ns, const Cut_field_double& cut_field_temperature, DoubleTabFT_cut_cell_scalar& valeur_remplissage)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();
  const IJK_Field_double& surface_interface_old = cut_cell_disc.get_interfaces().get_surface_interface_old();

  double dx = cut_cell_disc.get_splitting().get_grid_geometry().get_constant_delta(0);
  double dy = cut_cell_disc.get_splitting().get_grid_geometry().get_constant_delta(1);
  double dz = cut_cell_disc.get_splitting().get_grid_geometry().get_constant_delta(2);

  int statut_diphasique_naissant = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::NAISSANT);
  int statut_diphasique_petit = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::DESEQUILIBRE_FINAL);
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

      if (!cut_cell_disc.get_splitting().within_ghost(i, j, k, 1, 1))
        continue;

      double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);
      int est_naissant = cut_cell_disc.get_interfaces().est_pure(cut_cell_disc.get_interfaces().I(i,j,k));
      int phase = est_naissant ? 1 - IJK_Interfaces::convert_indicatrice_to_phase(cut_cell_disc.get_interfaces().I(i,j,k)) : ((cut_cell_disc.get_interfaces().below_small_threshold(next_indicatrice)) ? 1 : 0); // phase de la cellule petite ou naissante

      double normal_x = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,0);
      double normal_y = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,1);
      double normal_z = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,2);

      double velocity_x = cut_cell_disc.get_interfaces().get_vitesse_deplacement_interface()(n,0);
      double velocity_y = cut_cell_disc.get_interfaces().get_vitesse_deplacement_interface()(n,1);
      double velocity_z = cut_cell_disc.get_interfaces().get_vitesse_deplacement_interface()(n,2);

      double next_bary_x = dx*cut_cell_disc.get_interfaces().get_barycentre(true, 0, phase, i,j,k);
      double next_bary_y = dy*cut_cell_disc.get_interfaces().get_barycentre(true, 1, phase, i,j,k);
      double next_bary_z = dz*cut_cell_disc.get_interfaces().get_barycentre(true, 2, phase, i,j,k);

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
          double next_bary_deplace_normalise = select_dir(direction_a_modifier, next_bary_deplace_x/dx, next_bary_deplace_y/dy, next_bary_deplace_z/dz);
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

      double old_bary_x = dx*((i_old - i) + cut_cell_disc.get_interfaces().get_barycentre(false, 0, phase, i_old,j_old,k_old));
      double old_bary_y = dy*((j_old - j) + cut_cell_disc.get_interfaces().get_barycentre(false, 1, phase, i_old,j_old,k_old));
      double old_bary_z = dz*((k_old - k) + cut_cell_disc.get_interfaces().get_barycentre(false, 2, phase, i_old,j_old,k_old));

      if (n_old < 0)
        {
          Cerr << "Dans Cut_cell_schema_auxiliaire::calcule_valeur_remplissage_semi_lagrangien, on a n_old < 0." << finl;
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
      valeur_remplissage(n) = temperature_old + dTdn * ((next_bary_deplace_x - old_bary_x)*normal_x + (next_bary_deplace_y - old_bary_y)*normal_y + (next_bary_deplace_z - old_bary_z)*normal_z);
    }
}

void Cut_cell_schema_auxiliaire::calcule_valeur_remplissage_semi_lagrangien_interpolate(double timestep, const ArrOfDouble& interfacial_temperature, const IJK_Field_double& temperature_ft, const Cut_field_double& cut_field_temperature, DoubleTabFT_cut_cell_scalar& valeur_remplissage)
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

  int statut_diphasique_naissant = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::NAISSANT);
  int statut_diphasique_petit = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::DESEQUILIBRE_FINAL);
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

      if (!cut_cell_disc.get_splitting().within_ghost(i, j, k, 1, 1))
        continue;

      double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);
      int est_naissant = cut_cell_disc.get_interfaces().est_pure(cut_cell_disc.get_interfaces().I(i,j,k));
      int phase = est_naissant ? 1 - IJK_Interfaces::convert_indicatrice_to_phase(cut_cell_disc.get_interfaces().I(i,j,k)) : ((cut_cell_disc.get_interfaces().below_small_threshold(next_indicatrice)) ? 1 : 0); // phase de la cellule petite ou naissante

      double velocity_x = cut_cell_disc.get_interfaces().get_vitesse_deplacement_interface()(n,0);
      double velocity_y = cut_cell_disc.get_interfaces().get_vitesse_deplacement_interface()(n,1);
      double velocity_z = cut_cell_disc.get_interfaces().get_vitesse_deplacement_interface()(n,2);

      double next_bary_x = dx*cut_cell_disc.get_interfaces().get_barycentre(true, 0, phase, i,j,k);
      double next_bary_y = dy*cut_cell_disc.get_interfaces().get_barycentre(true, 1, phase, i,j,k);
      double next_bary_z = dz*cut_cell_disc.get_interfaces().get_barycentre(true, 2, phase, i,j,k);

      double next_bary_deplace_x = next_bary_x - velocity_x*timestep;
      double next_bary_deplace_y = next_bary_y - velocity_y*timestep;
      double next_bary_deplace_z = next_bary_z - velocity_z*timestep;

      double coordinates_x = next_bary_deplace_x + cut_cell_disc.get_splitting().get_coord_of_dof_along_dir(DIRECTION_I, i, IJK_Splitting::NODES);
      double coordinates_y = next_bary_deplace_y + cut_cell_disc.get_splitting().get_coord_of_dof_along_dir(DIRECTION_J, j, IJK_Splitting::NODES);
      double coordinates_z = next_bary_deplace_z + cut_cell_disc.get_splitting().get_coord_of_dof_along_dir(DIRECTION_K, k, IJK_Splitting::NODES);
      Vecteur3 coordinates(coordinates_x, coordinates_y, coordinates_z);

      int status = -2;
      const int tolerate_not_within_tetrahedron = std::max(1, tolerate_not_within_tetrahedron_);
      double temperature_interpolate = ijk_interpolate_cut_cell_using_interface(false, phase, temperature_ft, cut_field_temperature, interfacial_temperature, coordinates, tolerate_not_within_tetrahedron, status);
      valeur_remplissage(n) = temperature_interpolate;

      if (status == -1)
        {
          Cerr << "DEBUG status=-1 T_remplissage=" << valeur_remplissage(n) << finl;
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

  Cerr << "Bilan des statuts pour Cut_cell_schema_auxiliaire::calcule_valeur_remplissage_semi_lagrangien_interpolate : " << finl;
  Cerr << "    -1      " << count_status_not_found  << finl;
  Cerr << "   <10      " << count_status_below_10   << finl;
  Cerr << "  <100      " << count_status_below_100  << finl;
  Cerr << " <1000      " << count_status_below_1000 << finl;
  Cerr << " <10000     " << count_status_below_10000 << finl;
  Cerr << " <100000    " << count_status_below_100000 << finl;
  Cerr << " <1000000   " << count_status_below_1000000 << finl;
  Cerr << ">=1000000   " << count_status_above_1000000 << finl;
}

void Cut_cell_schema_auxiliaire::add_dying_cells(const Cut_field_vector3_double& cut_field_total_velocity, Cut_field_double& cut_field_temperature, bool write_flux, Cut_field_vector3_double& cut_field_current_fluxes)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();

  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(0));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(1));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(2));

  int statut_diphasique = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::MOURRANT);
  int index_min = cut_cell_disc.get_statut_diphasique_value_index(statut_diphasique);
  int index_max = cut_cell_disc.get_statut_diphasique_value_index(statut_diphasique+1);
  for (int index = index_min; index < index_max; index++)
    {
      int n = cut_cell_disc.get_n_from_statut_diphasique_index(index);

      Int3 ijk = cut_cell_disc.get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      if (!cut_cell_disc.get_splitting().within_ghost(i, j, k, 1, 1))
        continue;

      double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);
      assert(cut_cell_disc.get_interfaces().est_pure(next_indicatrice));
      int phase = 1 - IJK_Interfaces::convert_indicatrice_to_phase(next_indicatrice); // phase de la cellule mourrante

      // Cette variable next_nonzero_indicatrice est utilisee plutot que old_indicatrice,
      // par coherence avec l'indicatrice des cellules voisines, qui peut-etre est next mais
      // peut-etre est egalement old si elle est egalement mourrante.
      double next_nonzero_indicatrice = cut_cell_disc.get_interfaces().In_nonzero(phase,i,j,k);
      assert(next_nonzero_indicatrice == ((phase == 0) ? 1 - cut_cell_disc.get_interfaces().I(i,j,k) : cut_cell_disc.get_interfaces().I(i,j,k)));

      double temperature_centre = (phase == 0) ? cut_field_temperature.diph_v_(n) : cut_field_temperature.diph_l_(n);

      if (std::abs(temperature_centre) < 1e-24)
        {
          continue;
        }

      // L'indicatrice precedente est utilisee. En effet, pour ne pas perdre d'information la temperature des cellules
      // mortes reste associee a l'indicatrice precedente dans tous les cas.
      double quantite_totale = (phase == 0) ? -next_nonzero_indicatrice*cut_field_temperature.diph_v_(n) : -next_nonzero_indicatrice*cut_field_temperature.diph_l_(n);


      double flux[6] = {flux_naive_(n,0), flux_naive_(n,1), flux_naive_(n,2), flux_naive_(n,3), flux_naive_(n,4), flux_naive_(n,5)};

      Cut_cell_correction_petites_cellules::modification_flux_petites_cellules(correction_petites_cellules_, quantite_totale, flux);
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
              double next_nonzero_indicatrice_decale = cut_cell_disc.get_interfaces().In_nonzero(phase,i+di_decale,j+dj_decale,k+dk_decale);
              if (n_decale >= 0)
                {
                  // Si la cellule voisine est egalement morte, on utilise l'indicatrice precedente pour calculer
                  // la variation de la temperature, sinon, on utilise l'indicatrice du pas de temps actuel/suivant.
                  if (phase == 0)
                    {
                      cut_field_temperature.diph_v_(n_decale) -= (flux[num_face]/somme_flux)*quantite_totale/next_nonzero_indicatrice_decale;
                      cut_field_temperature.diph_v_(n) += (flux[num_face]/somme_flux)*quantite_totale/next_nonzero_indicatrice;
                    }
                  else
                    {
                      cut_field_temperature.diph_l_(n_decale) -= (flux[num_face]/somme_flux)*quantite_totale/next_nonzero_indicatrice_decale;
                      cut_field_temperature.diph_l_(n) += (flux[num_face]/somme_flux)*quantite_totale/next_nonzero_indicatrice;
                    }
                }
              else
                {
                  cut_field_temperature.pure_(i+di_decale,j+dj_decale,k+dk_decale) -= (flux[num_face]/somme_flux)*quantite_totale;
                  if (phase == 0)
                    {
                      cut_field_temperature.diph_v_(n) += (flux[num_face]/somme_flux)*quantite_totale/next_nonzero_indicatrice;
                    }
                  else
                    {
                      cut_field_temperature.diph_l_(n) += (flux[num_face]/somme_flux)*quantite_totale/next_nonzero_indicatrice;
                    }
                }


              if (write_flux)
                {
                  if (num_face < 3)
                    {
                      if (phase == 0)
                        {
                          cut_field_current_fluxes[dir].diph_v_(n) += (flux[num_face]/somme_flux)*quantite_totale;
                        }
                      else
                        {
                          cut_field_current_fluxes[dir].diph_l_(n) += (flux[num_face]/somme_flux)*quantite_totale;
                        }
                    }
                  else
                    {
                      if (n_decale >= 0)
                        {
                          if (phase == 0)
                            {
                              cut_field_current_fluxes[dir].diph_v_(n_decale) -= (flux[num_face]/somme_flux)*quantite_totale;
                            }
                          else
                            {
                              cut_field_current_fluxes[dir].diph_l_(n_decale) -= (flux[num_face]/somme_flux)*quantite_totale;
                            }
                        }
                      else
                        {
                          cut_field_current_fluxes[dir].pure_(i+di_decale,j+dj_decale,k+dk_decale) -= (flux[num_face]/somme_flux)*quantite_totale;
                        }
                    }
                }
            }
        }
    }

  if (write_flux)
    {
      cut_field_current_fluxes[0].echange_espace_virtuel(1);
      cut_field_current_fluxes[1].echange_espace_virtuel(1);
      cut_field_current_fluxes[2].echange_espace_virtuel(1);
    }
}

void Cut_cell_schema_auxiliaire::add_small_nascent_cells(const Cut_field_vector3_double& cut_field_total_velocity, Cut_field_double& cut_field_temperature, bool write_flux, Cut_field_vector3_double& cut_field_current_fluxes)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();

  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(0));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(1));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(2));

  int statut_diphasique_naissant = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::NAISSANT);
  int statut_diphasique_petit = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::DESEQUILIBRE_FINAL);
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

      if (!cut_cell_disc.get_splitting().within_ghost(i, j, k, 1, 1))
        continue;

      double old_indicatrice = cut_cell_disc.get_interfaces().I(i,j,k);
      double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);
      int est_naissant = cut_cell_disc.get_interfaces().est_pure(old_indicatrice);
      int phase = est_naissant ? 1 - IJK_Interfaces::convert_indicatrice_to_phase(old_indicatrice) : ((cut_cell_disc.get_interfaces().below_small_threshold(next_indicatrice)) ? 1 : 0); // phase de la cellule petite ou naissante
      assert(est_naissant || cut_cell_disc.get_interfaces().next_below_small_threshold_for_phase(phase, cut_cell_disc.get_interfaces().I(i,j,k), next_indicatrice));


      double temperature_centre = (phase == 0) ? cut_field_temperature.diph_v_(n) : cut_field_temperature.diph_l_(n);
      double val_remplissage = temperature_remplissage_(n);
      double quantite_totale = (phase == 0) ? (1 - next_indicatrice)*(val_remplissage - temperature_centre) : next_indicatrice*(val_remplissage - temperature_centre);

      // The dummy value DMINFLOAT+1 indicates that no filling should be made.
      if (val_remplissage == DMINFLOAT+1)
        {
          if (phase == 0)
            {
              assert((cut_field_total_velocity[0].from_ijk_and_phase(i,j,k, 0) == 0) && (cut_field_total_velocity[1].from_ijk_and_phase(i,j,k, 0) == 0) && (cut_field_total_velocity[2].from_ijk_and_phase(i,j,k, 0) == 0));
            }
          else
            {
              assert((cut_field_total_velocity[0].from_ijk_and_phase(i,j,k, 1) == 0) && (cut_field_total_velocity[1].from_ijk_and_phase(i,j,k, 1) == 0) && (cut_field_total_velocity[2].from_ijk_and_phase(i,j,k, 1) == 0));
            }
          continue;
        }

      if (std::abs(temperature_centre - val_remplissage) < 1e-24)
        {
          continue;
        }

      // Note: std::abs(val_remplissage) is compared to std::numeric_limits<double>::min(), as a naive std::abs(val_remplissage) > 0 does not prevent arithmetic errors upon division.
      if ((std::abs(val_remplissage) > std::numeric_limits<double>::min()) && (std::abs(temperature_centre - val_remplissage)/val_remplissage) < 1e-24)
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
          bool decale_also_dying = (cut_cell_disc.get_interfaces().phase_mourrante(phase, i+di_decale,j+dj_decale,k+dk_decale));
          bool decale_nascent = (cut_cell_disc.get_interfaces().phase_naissante(phase, i+di_decale,j+dj_decale,k+dk_decale));
          bool decale_small = cut_cell_disc.get_interfaces().next_below_small_threshold_for_phase(phase, old_indicatrice_decale, next_indicatrice_decale);
          bool decale_smaller = (phase == 0) ? (1 - next_indicatrice_decale) < (1 - next_indicatrice) : (next_indicatrice_decale) < (next_indicatrice);
          if (decale_also_dying || (est_naissant && decale_nascent && decale_smaller) || ((!est_naissant) && decale_nascent) || ((!est_naissant) && decale_small && decale_smaller))
            {
              flux_max[num_face] = 0;
            }
          else
            {
              double temperature_decale = -1e37;
              int n_decale = cut_cell_disc.get_n(i+di_decale, j+dj_decale, k+dk_decale);
              if (n_decale >= 0)
                {
                  temperature_decale = (phase == 0) ? cut_field_temperature.diph_v_(n_decale) : cut_field_temperature.diph_l_(n_decale);
                }
              else
                {
                  temperature_decale = (phase == (IJK_Interfaces::convert_indicatrice_to_phase(next_indicatrice_decale)))*cut_field_temperature.pure_(i+di_decale,j+dj_decale,k+dk_decale);
                }

              total_velocity[num_face] = cut_field_total_velocity[dir].from_ijk_and_phase(i+di,j+dj,k+dk, phase);

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
                }
              else
                {
                  double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().In(i+di,j+dj,k+dk) : cut_cell_disc.get_interfaces().In(i+di,j+dj,k+dk);
                  assert((surface_efficace == 0) || (surface_efficace == 1));
                  if (surface_efficace > 0)
                    {
                      double next_volume_decale = (phase == 0) ? 1 - next_indicatrice_decale : next_indicatrice_decale;
                      double next_volume = (phase == 0) ? 1 - next_indicatrice : next_indicatrice;
                      flux_max[num_face] = (temperature_decale - temperature_centre)/(1/next_volume_decale + 1/next_volume);
                      flux_max[num_face] = flux[num_face] >= 0 ? std::max(0., flux_max[num_face]) : std::max(0., -flux_max[num_face]);
                    }
                }
            }
        }

      if (no_static_update_)
        {
          if ((std::abs(total_velocity[0]) + std::abs(total_velocity[1]) + std::abs(total_velocity[2]) + std::abs(total_velocity[3]) + std::abs(total_velocity[4]) + std::abs(total_velocity[5])) == 0)
            {
              continue;
            }
        }

      Cut_cell_correction_petites_cellules::modification_flux_petites_cellules(correction_petites_cellules_, quantite_totale, flux);
      double somme_flux = Cut_cell_correction_petites_cellules::calcul_somme_flux(flux);
      assert(somme_flux != 0.);
      if (somme_flux == 0.)
        {
          continue;
        }

      Cut_cell_correction_petites_cellules::limitation_flux_avec_flux_max(correction_petites_cellules_, quantite_totale, somme_flux, flux_max, flux);

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
                  assert(phase == IJK_Interfaces::convert_indicatrice_to_phase(next_indicatrice_decale));
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


              if (write_flux)
                {
                  if (num_face < 3)
                    {
                      if (phase == 0)
                        {
                          cut_field_current_fluxes[dir].diph_v_(n) += (flux[num_face]/somme_flux)*quantite_totale;
                          if (n_decale < 0)
                            {
                              assert(phase == IJK_Interfaces::convert_indicatrice_to_phase(next_indicatrice_decale));
                              cut_field_current_fluxes[dir].pure_(i,j,k) += (flux[num_face]/somme_flux)*quantite_totale;
                            }
                        }
                      else
                        {
                          cut_field_current_fluxes[dir].diph_l_(n) += (flux[num_face]/somme_flux)*quantite_totale;
                          if (n_decale < 0)
                            {
                              assert(phase == IJK_Interfaces::convert_indicatrice_to_phase(next_indicatrice_decale));
                              cut_field_current_fluxes[dir].pure_(i,j,k) += (flux[num_face]/somme_flux)*quantite_totale;
                            }
                        }
                    }
                  else
                    {
                      if (n_decale >= 0)
                        {
                          if (phase == 0)
                            {
                              cut_field_current_fluxes[dir].diph_v_(n_decale) -= (flux[num_face]/somme_flux)*quantite_totale;
                            }
                          else
                            {
                              cut_field_current_fluxes[dir].diph_l_(n_decale) -= (flux[num_face]/somme_flux)*quantite_totale;
                            }
                        }
                      else
                        {
                          cut_field_current_fluxes[dir].pure_(i+di_decale,j+dj_decale,k+dk_decale) -= (flux[num_face]/somme_flux)*quantite_totale;
                        }
                    }
                }
            }
        }
    }

  if (write_flux)
    {
      cut_field_current_fluxes[0].echange_espace_virtuel(1);
      cut_field_current_fluxes[1].echange_espace_virtuel(1);
      cut_field_current_fluxes[2].echange_espace_virtuel(1);
    }
}
