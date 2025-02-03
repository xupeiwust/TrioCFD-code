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

#include <Cut_cell_diffusion_auxiliaire.h>
#include <Cut_cell_FT_Disc.h>
#include <IJK_Thermal_base.h>
#include <Probleme_FTD_IJK_cut_cell.h>
#include <Process.h>
#include <stat_counters.h>
#include <IJK_Navier_Stokes_tools.h>
#include <IJK_Navier_Stokes_tools_cut_cell.h>
#include <Param.h>

Implemente_instanciable_sans_constructeur(Cut_cell_diffusion_auxiliaire, "Cut_cell_diffusion_auxiliaire", Cut_cell_schema_auxiliaire) ;

Cut_cell_diffusion_auxiliaire::Cut_cell_diffusion_auxiliaire()
{
  methode_valeur_remplissage_ = METHODE_TEMPERATURE_REMPLISSAGE::COPIE_DIRECTE;

  deactivate_correction_petites_cellules_diffusion_ = 0;
  correction_petites_cellules_ = CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE_AVEC_LIMITATION_2;

  no_static_update_ = false;
}

Sortie& Cut_cell_diffusion_auxiliaire::printOn(Sortie& os) const
{
  Objet_U::printOn(os);
  return os;
}

Entree& Cut_cell_diffusion_auxiliaire::readOn(Entree& is)
{
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades(is);
  return is;
}

void Cut_cell_diffusion_auxiliaire::set_param(Param& param)
{
  Cut_cell_schema_auxiliaire::set_param(param);

  param.ajouter_flag("deactivate_correction_petites_cellules_diffusion", &deactivate_correction_petites_cellules_diffusion_);
}

double Cut_cell_diffusion_auxiliaire::dying_cells_flux(int num_face, int phase, int n, const Cut_field_vector3_double& cut_field_total_velocity, const Cut_field_double& cut_field_temperature)
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

  const DoubleTabFT_cut_cell_scalar& flux_interface_efficace = select_flux_interface(phase);

  double normal_x = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,0);
  double normal_y = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,1);
  double normal_z = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,2);
  int sign_flux_interf = (flux_interface_efficace(n) == 0.) ? 0 : (1 - 2*phase)*(2*(flux_interface_efficace(n) > 0) - 1);

  double normal_to_face = sign*select_dir(dir, normal_x, normal_y, normal_z);

  int n_face = cut_cell_disc.get_n_face(num_face, n, i, j, k);
  if (n_face >= 0)
    {
      double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir) : cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir);
      if (surface_efficace > 0)
        {
          return -sign*surface_efficace*normal_to_face*sign_flux_interf;
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
          return -sign*surface_efficace*normal_to_face*sign_flux_interf;
        }
      else
        {
          return 0.;
        }
    }
}

double Cut_cell_diffusion_auxiliaire::small_nascent_cells_flux(int num_face, int phase, int n, const Cut_field_vector3_double& cut_field_total_velocity, const Cut_field_double& cut_field)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field.get_cut_cell_disc();

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

  const DoubleTabFT_cut_cell_scalar& flux_interface_efficace = select_flux_interface(phase);

  double normal_x = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,0);
  double normal_y = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,1);
  double normal_z = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,2);
  int sign_flux_interf = (flux_interface_efficace(n) == 0.) ? 0 : (1 - 2*phase)*(2*(flux_interface_efficace(n) > 0) - 1);

  double normal_to_face = sign*select_dir(dir, normal_x, normal_y, normal_z);

  int n_face = cut_cell_disc.get_n_face(num_face, n, i, j, k);
  if (n_face >= 0)
    {
      double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir) : cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir);
      if (surface_efficace > 0)
        {
          return -sign*surface_efficace*normal_to_face*sign_flux_interf;
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
          return -sign*surface_efficace*normal_to_face*sign_flux_interf;
        }
      else
        {
          return 0.;
        }
    }
}

void Cut_cell_diffusion_auxiliaire::associer(DoubleTabFT_cut_cell_scalar& flux_interface_efficace)
{
  flux_interface_efficace_ptr_ = &flux_interface_efficace;
}

const DoubleTabFT_cut_cell_scalar& Cut_cell_diffusion_auxiliaire::select_flux_interface(int phase)
{
  if (flux_interface_efficace_ptr_ == nullptr)
    {
      Cerr << "Invalid pointer flux_interface_efficace_ in Cut_cell_diffusion_auxiliaire::select_flux_interface." << finl;
      Process::exit();
    }
  return *flux_interface_efficace_ptr_;
}
