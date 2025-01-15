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
#include <IJK_Field_vector.h>
#include <Cut_cell_FT_Disc.h>
#include <IJK_Thermal_base.h>
#include <Process.h>
#include <IJK_FT_cut_cell.h>
#include <stat_counters.h>
#include <IJK_Navier_Stokes_tools.h>
#include <IJK_Navier_Stokes_tools_cut_cell.h>
#include <Param.h>

Implemente_instanciable_sans_constructeur(Cut_cell_convection_auxiliaire, "Cut_cell_convection_auxiliaire", Cut_cell_schema_auxiliaire) ;

Cut_cell_convection_auxiliaire::Cut_cell_convection_auxiliaire()
{
  methode_temperature_remplissage_ = METHODE_TEMPERATURE_REMPLISSAGE::NON_INITIALISE;
  correction_petites_cellules_ = CORRECTION_PETITES_CELLULES::CORRECTION_SYMETRIQUE;

  no_static_update_ = true;
}

Sortie& Cut_cell_convection_auxiliaire::printOn(Sortie& os) const
{
  Objet_U::printOn(os);
  return os;
}

Entree& Cut_cell_convection_auxiliaire::readOn(Entree& is)
{
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades(is);
  return is;
}

void Cut_cell_convection_auxiliaire::set_param(Param& param)
{
  Cut_cell_schema_auxiliaire::set_param(param);
}

double Cut_cell_convection_auxiliaire::dying_cells_flux(int num_face, int phase, int n, const Cut_field_vector3_double& cut_field_total_velocity, const Cut_field_double& cut_field_temperature)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();

  Int3 ijk = cut_cell_disc.get_ijk(n);
  int i = ijk[0];
  int j = ijk[1];
  int k = ijk[2];

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

  double temperature_centre = (phase == 0) ? cut_field_temperature.diph_v_(n) : cut_field_temperature.diph_l_(n);

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
          double total_velocity = (phase == 0) ? cut_field_total_velocity[dir].diph_v_(n_face) : cut_field_total_velocity[dir].diph_l_(n_face);
          double temperature = (sign*total_velocity < 0) ? temperature_decale : temperature_centre;
          return -sign*surface_efficace*temperature*total_velocity;
        }
      else
        {
          return 0.;
        }
    }
  else
    {
      double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().I(i+di,j+dj,k+dk) : cut_cell_disc.get_interfaces().I(i+di,j+dj,k+dk);
      assert((surface_efficace == 0) || (surface_efficace == 1));
      if (surface_efficace > 0)
        {
          double total_velocity = cut_field_total_velocity[dir].pure_(i+di,j+dj,k+dk);
          double temperature = (sign*total_velocity < 0) ? temperature_decale : temperature_centre;
          return -sign*surface_efficace*temperature*total_velocity;
        }
      else
        {
          return 0.;
        }
    }
}

double Cut_cell_convection_auxiliaire::small_nascent_cells_flux(int num_face, int phase, int n, const Cut_field_vector3_double& cut_field_total_velocity, const Cut_field_double& cut_field_temperature)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();

  Int3 ijk = cut_cell_disc.get_ijk(n);
  int i = ijk[0];
  int j = ijk[1];
  int k = ijk[2];

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
      if (surface_efficace > 0)
        {
          double total_velocity = (phase == 0) ? cut_field_total_velocity[dir].diph_v_(n_face) : cut_field_total_velocity[dir].diph_l_(n_face);
          double temperature = temperature_decale;
          return -sign*surface_efficace*temperature*total_velocity;
        }
      else
        {
          return 0.;
        }
    }
  else
    {
      double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().In(i+di,j+dj,k+dk) : cut_cell_disc.get_interfaces().In(i+di,j+dj,k+dk);
      assert((surface_efficace == 0) || (surface_efficace == 1));
      if (surface_efficace > 0)
        {
          double total_velocity = cut_field_total_velocity[dir].pure_(i+di,j+dj,k+dk);
          double temperature = temperature_decale;
          return -sign*surface_efficace*temperature*total_velocity;
        }
      else
        {
          return 0.;
        }
    }
}


