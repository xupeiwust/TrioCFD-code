/****************************************************************************
* Copyright (c) 2023, CEA
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
// File      : IJK_Thermal_cut_cell.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Thermal_cut_cell.h>
#include <IJK_FT_cut_cell.h>
#include <IJK_FT.h>
#include <DebogIJK.h>
#include <IJK_Navier_Stokes_tools.h>
#include <IJK_Navier_Stokes_tools_cut_cell.h>
#include <Cut_cell_convection_auxiliaire.h>
#include <Cut_cell_diffusion_auxiliaire.h>
#include <OpConvQuickIJKScalar_cut_cell.h>

Implemente_instanciable_sans_constructeur( IJK_Thermal_cut_cell, "IJK_Thermal_cut_cell", IJK_Thermal_base ) ;

IJK_Thermal_cut_cell::IJK_Thermal_cut_cell()
{
  single_phase_=0;
  conserv_energy_global_=0; // Note : doit etre zero sinon la rustine est appliquee

  cut_cell_conv_scheme_.scheme = CUT_CELL_SCHEMA_CONVECTION::QUICK_OU_CENTRE2_STENCIL;


  temperature_ = std::make_shared<Cut_field_double>();
  div_coeff_grad_T_volume_ = std::make_shared<Cut_field_double>();
  d_temperature_ = std::make_shared<Cut_field_double>();

  temperature_post_dying_ = std::make_shared<Cut_field_double>();
  temperature_post_regular_ = std::make_shared<Cut_field_double>();
  temperature_post_convection_ = std::make_shared<Cut_field_double>();
  temperature_post_diff_regular_ = std::make_shared<Cut_field_double>();
}

Sortie& IJK_Thermal_cut_cell::printOn( Sortie& os ) const
{
  Nom front_space = "    ";
  Nom end_space = " ";
  Nom escape = "\n";

  IJK_Thermal_base::printOn( os );

  os<< "    type_T_source " << type_T_source_ << "\n";

  if (lambda_variable_)
    os<< "    lambda_variable \n";
  if (conserv_energy_global_)
    os<< "    conserv_energy_global \n";


  os<< "   convection_auxiliaire " << convective_correction_;

  os<< "   diffusion_auxiliaire " << diffusive_correction_;

  os<< "  \n}";
  return os;
}

Entree& IJK_Thermal_cut_cell::readOn( Entree& is )
{
  IJK_Thermal_base::readOn( is );
  return is;
}

void IJK_Thermal_cut_cell::set_param( Param& param )
{
  IJK_Thermal_base::set_param(param);
  param.ajouter("type_temperature_convection_form", &type_temperature_convection_form_);
  param.dictionnaire("non conservative",0);
  param.dictionnaire("conservative",1);
  param.ajouter("conserv_energy_global", &conserv_energy_global_);
  param.ajouter_flag("postraiter_champs_intermediaires", &postraiter_champs_intermediaires_);

  param.ajouter_flag("deactivate_diffusion_interface", &deactivate_diffusion_interface_);

  param.ajouter_flag("runge_kutta_fluxes_convection", &runge_kutta_fluxes_convection_);
  param.ajouter_flag("runge_kutta_fluxes_diffusion", &runge_kutta_fluxes_diffusion_);
  param.ajouter_flag("runge_kutta_fluxes_pas_de_correction_convection", &runge_kutta_fluxes_pas_de_correction_convection_);
  param.ajouter_flag("runge_kutta_fluxes_pas_de_correction_diffusion", &runge_kutta_fluxes_pas_de_correction_diffusion_);

  param.ajouter("runge_kutta_restriction_leniency_convection", &runge_kutta_restriction_leniency_convection_);
  param.ajouter("runge_kutta_restriction_leniency_diffusion", &runge_kutta_restriction_leniency_diffusion_);

  param.ajouter("convection_auxiliaire", &convective_correction_, Param::REQUIRED);
  param.ajouter("diffusion_auxiliaire", &diffusive_correction_, Param::REQUIRED);

  param.ajouter("cut_cell_schema_convection", (int*)&cut_cell_conv_scheme_.scheme);
  param.dictionnaire("quick_ou_centre2_stencil", (int)CUT_CELL_SCHEMA_CONVECTION::QUICK_OU_CENTRE2_STENCIL);
  param.dictionnaire("quick_ou_centre2_perpendicular_distance", (int)CUT_CELL_SCHEMA_CONVECTION::QUICK_OU_CENTRE2_PERPENDICULAR_DISTANCE);
  param.dictionnaire("quick_ou_lineaire2_stencil", (int)CUT_CELL_SCHEMA_CONVECTION::QUICK_OU_LINEAIRE2_STENCIL);
  param.dictionnaire("quick_ou_lineaire2_perpendicular_distance", (int)CUT_CELL_SCHEMA_CONVECTION::QUICK_OU_LINEAIRE2_PERPENDICULAR_DISTANCE);
  param.dictionnaire("quick_ou_amont_stencil", (int)CUT_CELL_SCHEMA_CONVECTION::QUICK_OU_AMONT_STENCIL);
  param.dictionnaire("quick_ou_amont_perpendicular_distance", (int)CUT_CELL_SCHEMA_CONVECTION::QUICK_OU_AMONT_PERPENDICULAR_DISTANCE);
  param.dictionnaire("centre2", (int)CUT_CELL_SCHEMA_CONVECTION::CENTRE2);
  param.dictionnaire("lineaire2", (int)CUT_CELL_SCHEMA_CONVECTION::LINEAIRE2);
  param.dictionnaire("amont", (int)CUT_CELL_SCHEMA_CONVECTION::AMONT);


  param.ajouter("verbosite", &verbosite_);
}

int IJK_Thermal_cut_cell::initialize(const IJK_Splitting& splitting, const int idx)
{
  Cout << que_suis_je() << "::initialize()" << finl;

  // Cut-cell variables
  ref_ijk_ft_cut_cell_ = ref_cast(IJK_FT_cut_cell, ref_ijk_ft_.valeur());

  Cut_field_double& cut_field_temperature                = static_cast<Cut_field_double&>(*temperature_);
  cut_field_temperature.allocate_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), splitting, IJK_Splitting::ELEM, ghost_cells_); // Overrides the allocate in IJK_Thermal_base::initialize

  Cut_field_vector3_double& cut_field_current_fluxes_conv = static_cast<Cut_field_vector3_double&>(current_fluxes_conv_);
  Cut_field_vector3_double& cut_field_current_fluxes_diff = static_cast<Cut_field_vector3_double&>(current_fluxes_diff_);
  Cut_field_vector3_double& cut_field_RK3_F_fluxes_conv   = static_cast<Cut_field_vector3_double&>(RK3_F_fluxes_conv_);
  Cut_field_vector3_double& cut_field_RK3_F_fluxes_diff   = static_cast<Cut_field_vector3_double&>(RK3_F_fluxes_diff_);
  if (runge_kutta_fluxes_convection_)
    {
      allocate_velocity_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), cut_field_current_fluxes_conv, splitting, 2);
      allocate_velocity_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), cut_field_RK3_F_fluxes_conv, splitting, 2);
    }

  if (runge_kutta_fluxes_diffusion_)
    {
      allocate_velocity_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), cut_field_current_fluxes_diff, splitting, 2);
      allocate_velocity_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), cut_field_RK3_F_fluxes_diff, splitting, 2);
    }

  cellule_rk_restreint_conv_l_.allocate(splitting, IJK_Splitting::ELEM, ghost_cells_);
  cellule_rk_restreint_conv_v_.allocate(splitting, IJK_Splitting::ELEM, ghost_cells_);
  cellule_rk_restreint_diff_l_.allocate(splitting, IJK_Splitting::ELEM, ghost_cells_);
  cellule_rk_restreint_diff_v_.allocate(splitting, IJK_Splitting::ELEM, ghost_cells_);

  if (postraiter_champs_intermediaires_)
    {
      Cut_field_double& cut_field_temperature_post_dying        = static_cast<Cut_field_double&>(*temperature_post_dying_);
      Cut_field_double& cut_field_temperature_post_regular      = static_cast<Cut_field_double&>(*temperature_post_regular_);
      Cut_field_double& cut_field_temperature_post_convection   = static_cast<Cut_field_double&>(*temperature_post_convection_);
      Cut_field_double& cut_field_temperature_post_diff_regular = static_cast<Cut_field_double&>(*temperature_post_diff_regular_);

      cut_field_temperature_post_dying.allocate_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), splitting, IJK_Splitting::ELEM, 2);
      cut_field_temperature_post_regular.allocate_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), splitting, IJK_Splitting::ELEM, 2);
      cut_field_temperature_post_convection.allocate_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), splitting, IJK_Splitting::ELEM, 2);
      cut_field_temperature_post_diff_regular.allocate_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), splitting, IJK_Splitting::ELEM, 2);
    }

  cut_cell_flux_diffusion_[0].associer_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());
  cut_cell_flux_diffusion_[1].associer_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());
  cut_cell_flux_diffusion_[2].associer_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());
  cut_cell_flux_convection_[0].associer_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());
  cut_cell_flux_convection_[1].associer_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());
  cut_cell_flux_convection_[2].associer_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());


  int nalloc = IJK_Thermal_base::initialize(splitting, idx);
  lambda_.allocate(splitting, IJK_Splitting::ELEM, 1);
  nalloc += 2;

  temperature_ft_.allocate(ref_ijk_ft_cut_cell_->get_splitting_ft(), IJK_Splitting::ELEM, 4);
  nalloc += 1;

  treatment_count_.allocate(splitting, IJK_Splitting::ELEM, 2);
  nalloc += 1;

  for (int next_time = 0; next_time < 2; next_time++)
    {
      diffusive_correction_.flux_interface_ft_[next_time].allocate(ref_ijk_ft_cut_cell_->get_splitting_ft(), IJK_Splitting::ELEM, 2);
      diffusive_correction_.flux_interface_ns_[next_time].allocate(ref_ijk_ft_cut_cell_->get_splitting_ns(), IJK_Splitting::ELEM, 2);
      diffusive_correction_.flux_interface_ft_[next_time].data() = 0.;
      diffusive_correction_.flux_interface_ns_[next_time].data() = 0.;
    }
  diffusive_correction_.flux_interface_efficace_.associer_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());

  convective_correction_.initialise(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());
  diffusive_correction_.initialise(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());

  Cut_field_double& cut_field_div_coeff_grad_T_volume = static_cast<Cut_field_double&>(*div_coeff_grad_T_volume_);
  cut_field_div_coeff_grad_T_volume.allocate_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), splitting, IJK_Splitting::ELEM, 2); // Overrides the allocate in IJK_Thermal_base::initialize

  Cut_field_double& cut_field_d_temperature           = static_cast<Cut_field_double&>(*d_temperature_);
  cut_field_d_temperature.allocate_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), splitting, IJK_Splitting::ELEM, 2); // Overrides the allocate in IJK_Thermal_base::initialize


  temperature_diffusion_op_.typer("OpDiffIJKScalar_cut_cell_double");
  temperature_diffusion_op_.initialize(splitting);
  temperature_diffusion_op_.set_uniform_lambda_liquid(lambda_liquid_);
  temperature_diffusion_op_.set_uniform_lambda_vapour(lambda_vapour_);
  temperature_diffusion_op_.set_lambda(lambda_);

  temperature_convection_op_.typer("OpConvQuickIJKScalar_cut_cell_double");
  temperature_convection_op_.initialize(splitting);

  if (fichier_reprise_temperature_ == "??")   // si on ne fait pas une reprise on initialise V
    {
      if (expression_T_init_ != "??")
        {
          compute_temperature_init();
        }
      else
        {
          Cerr << "Please provide initial conditions for temperature, either by an expression or a field for restart. "
               << "You should consider using either fichier_reprise_temperature or expression_T_init keywords. " << finl;
          Process::exit();
        }
    }

  // Already allocated if rho_cp_post
  if (!rho_cp_post_)
    {
      rho_cp_.allocate(splitting, IJK_Splitting::ELEM, 2);
      nalloc += 1;
    }

  // Compute initial energy :
  if (conserv_energy_global_)
    {
      Cerr << "conserv_energy_global_ is used." << finl;
      Process::exit();
    }
  if (type_temperature_convection_form_==1)
    {
      rho_cp_T_.allocate(splitting, IJK_Splitting::ELEM, 2);
      div_rho_cp_T_.allocate(splitting, IJK_Splitting::ELEM, 0);
      nalloc += 2;
    }
  Cout << "End of " << que_suis_je() << "::initialize()" << finl;
  return nalloc;
}

void IJK_Thermal_cut_cell::update_thermal_properties()
{
  //IJK_Thermal_base::update_thermal_properties();
  const IJK_Field_double& indic = ref_ijk_ft_cut_cell_->itfce().I();
  // Nombre de mailles du domaine NS :
  const int nx = indic.ni();
  const int ny = indic.nj();
  const int nz = indic.nk();
  const double rho_l = ref_ijk_ft_cut_cell_->get_rho_l();
  const double rho_v = ref_ijk_ft_cut_cell_->get_rho_v();
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        {
          double chi_l = indic(i,j,k);
          rho_cp_(i,j,k) = rho_l * cp_liquid_ * chi_l + rho_v * cp_vapour_ * (1.- chi_l);
          lambda_(i,j,k) = lambda_liquid_  * chi_l + (1.- chi_l) * lambda_vapour_ ;
        }

  lambda_.echange_espace_virtuel(lambda_.ghost());
  rho_cp_.echange_espace_virtuel(rho_cp_.ghost());
}

void IJK_Thermal_cut_cell::compute_temperature_init()
{
  Cut_field_double& cut_field_temperature = static_cast<Cut_field_double&>(*temperature_);

  if (cut_field_temperature.get_cut_cell_disc().get_n_tot() == 0)
    {
      IJK_Thermal_base::compute_temperature_init();
    }
  else
    {
      Cout << "Temperature initialization from expression \nTini = " << expression_T_init_ << finl;
      cut_field_temperature.set_field_data(expression_T_init_, ref_ijk_ft_->itfce().I(), 0.);
    }
}

void IJK_Thermal_cut_cell::recompute_temperature_init()
{
  Cut_field_double& cut_field_temperature = static_cast<Cut_field_double&>(*temperature_);

  if (cut_field_temperature.get_cut_cell_disc().get_n_tot() == 0)
    {
      IJK_Thermal_base::recompute_temperature_init();
    }
  else
    {
      Cout << "Temperature initialization from expression \nTini = " << expression_T_init_ << finl;
      cut_field_temperature.set_field_data(expression_T_init_, ref_ijk_ft_->itfce().In(), 0.);
    }
}

struct CutCell_GlobalInfo
{
  double overall;
  double overall_l;
  double overall_v;
  double pure;
  double diph_l;
  double diph_v;
  double diph_small;
  double diph_regular;
  double diph_nascent;
  double diph_dying;
};

// Cette fonction sert ajouter la pente liee a la correction a la pente liee au schema principal
// Le facteur vol_over_dt_surface existe car le flux calcule dans les routines de correction n'a
// pas la bonne unite.
void IJK_Thermal_cut_cell::add_flux_times_vol_over_dt_surface(double fractional_timestep, const Cut_field_vector3_double& cut_field_current_fluxes, Cut_field_vector3_double& cut_field_RK3_F_fluxes)
{
  for (int dir = 0; dir < 3; dir++)
    {
      const int ni = cut_field_current_fluxes[dir].ni();
      const int nj = cut_field_current_fluxes[dir].nj();
      const int nk = cut_field_current_fluxes[dir].nk();
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  double vol_over_dt = vol_/fractional_timestep;

                  cut_field_RK3_F_fluxes[dir].pure_(i,j,k) += cut_field_current_fluxes[dir].pure_(i,j,k)*vol_over_dt;
                }
            }
        }

      const Cut_cell_FT_Disc& cut_cell_disc = cut_field_current_fluxes[dir].get_cut_cell_disc();
      for (int n = 0; n < cut_cell_disc.get_n_tot(); n++)
        {
          double vol_over_dt = vol_/fractional_timestep;

          const DoubleTabFT_cut_cell_vector3& indicatrice_surfacique = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face();
          double indicatrice_surface = indicatrice_surfacique(n,dir);

          cut_field_RK3_F_fluxes[dir].diph_l_(n) += (indicatrice_surface == 0.)       ? 0. : cut_field_current_fluxes[dir].diph_l_(n)*vol_over_dt/indicatrice_surface;
          cut_field_RK3_F_fluxes[dir].diph_v_(n) += ((1 - indicatrice_surface) == 0.) ? 0. : cut_field_current_fluxes[dir].diph_v_(n)*vol_over_dt/(1 - indicatrice_surface);
        }
    }
}

void IJK_Thermal_cut_cell::set_rk_restreint(int rk_step, int rk_restriction_leniency, const Cut_cell_FT_Disc& cut_cell_disc, IJK_Field_int& cellule_rk_restreint_v, IJK_Field_int& cellule_rk_restreint_l)
{
  if (rk_step == 0)
    {
      cellule_rk_restreint_l.data()=0;
      cellule_rk_restreint_v.data()=0;
    }

  if (rk_restriction_leniency == 10)
    {
      cellule_rk_restreint_l.data()=1;
      cellule_rk_restreint_v.data()=1;
    }

  if (rk_restriction_leniency == 11)
    {
      for (int n = 0; n < cut_cell_disc.get_n_tot(); n++)
        {
          Int3 ijk = cut_cell_disc.get_ijk(n);
          int i = ijk[0];
          int j = ijk[1];
          int k = ijk[2];

          cellule_rk_restreint_v(i,j,k) = 1;
          cellule_rk_restreint_l(i,j,k) = 1;
        }
    }

  // Boucle sur les cellules qui disparaissent lors de ce sous pas de temps
  {
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

        if ((rk_restriction_leniency == 4) || (rk_restriction_leniency == 5))
          continue;

        if ((rk_restriction_leniency == 1) || (rk_restriction_leniency == 3))
          {
            double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);
            int phase = 1 - (int)next_indicatrice; // phase de la cellule mourrante

            if (phase == 0)
              {
                cellule_rk_restreint_v(i,j,k) = 1;
              }
            else
              {
                cellule_rk_restreint_l(i,j,k) = 1;
              }
          }
        else if ((rk_restriction_leniency == 0) || (rk_restriction_leniency == 2))
          {
            cellule_rk_restreint_v(i,j,k) = 1;
            cellule_rk_restreint_l(i,j,k) = 1;
          }
      }
  }

  // Boucle sur les cellules qui apparaissent ou bien qui sont petites lors de ce sous pas de temps
  {
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

        if (rk_step == 0 && ((rk_restriction_leniency == 5) || (rk_restriction_leniency == 4) || (rk_restriction_leniency == 3) || (rk_restriction_leniency == 2)))
          continue;

        if ((rk_restriction_leniency == 1) || (rk_restriction_leniency == 3) || (rk_restriction_leniency == 5))
          {
            double old_indicatrice = cut_cell_disc.get_interfaces().I(i,j,k);
            double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);
            int est_naissant = cut_cell_disc.get_interfaces().est_pure(old_indicatrice);
            int phase = est_naissant ? 1 - (int)old_indicatrice : ((cut_cell_disc.get_interfaces().below_small_threshold(next_indicatrice)) ? 1 : 0); // phase de la cellule petite ou naissante

            if (phase == 0)
              {
                cellule_rk_restreint_v(i,j,k) = 1;
              }
            else
              {
                cellule_rk_restreint_l(i,j,k) = 1;
              }
          }
        else if ((rk_restriction_leniency == 0) || (rk_restriction_leniency == 2) || (rk_restriction_leniency == 4))
          {
            cellule_rk_restreint_v(i,j,k) = 1;
            cellule_rk_restreint_l(i,j,k) = 1;
          }
      }
  }
}

void IJK_Thermal_cut_cell::perform_thermal_step(double total_timestep, int flag_rk, int rk_step)
{
  ref_cast(OpDiffIJKScalar_cut_cell_double, temperature_diffusion_op_.valeur()).initialise_cut_cell(false, cut_cell_flux_diffusion_, treatment_count_, new_treatment_);
  ref_cast(OpConvQuickIJKScalar_cut_cell_double, temperature_convection_op_.valeur()).initialise_cut_cell(cut_cell_conv_scheme_, false, cut_cell_flux_convection_, treatment_count_, new_treatment_);

  double fractional_timestep = (flag_rk) ? compute_fractionnal_timestep_rk3(total_timestep, rk_step) : total_timestep;

  if (debug_)
    Cerr << "Thermal Euler time-step" << finl;

  Cut_field_double& cut_field_temperature                   = static_cast<Cut_field_double&>(*temperature_);
  Cut_field_double& cut_field_temperature_post_regular      = static_cast<Cut_field_double&>(*temperature_post_regular_);
  Cut_field_double& cut_field_temperature_post_dying        = static_cast<Cut_field_double&>(*temperature_post_dying_);
  Cut_field_double& cut_field_temperature_post_convection   = static_cast<Cut_field_double&>(*temperature_post_convection_);
  Cut_field_double& cut_field_temperature_post_diff_regular = static_cast<Cut_field_double&>(*temperature_post_diff_regular_);
  Cut_field_double& cut_field_d_temperature                 = static_cast<Cut_field_double&>(*d_temperature_);

  Cut_field_vector3_double& cut_field_current_fluxes_conv = static_cast<Cut_field_vector3_double&>(current_fluxes_conv_);
  Cut_field_vector3_double& cut_field_current_fluxes_diff = static_cast<Cut_field_vector3_double&>(current_fluxes_diff_);
  Cut_field_vector3_double& cut_field_RK3_F_fluxes_conv = static_cast<Cut_field_vector3_double&>(RK3_F_fluxes_conv_);
  Cut_field_vector3_double& cut_field_RK3_F_fluxes_diff = static_cast<Cut_field_vector3_double&>(RK3_F_fluxes_diff_);

  CutCell_GlobalInfo ene_ini;
  CutCell_GlobalInfo d_ene_Conv;
  CutCell_GlobalInfo ene_postconv_switch;
  CutCell_GlobalInfo ene_postconv_dying;
  CutCell_GlobalInfo ene_postconv_small_nascent;
  CutCell_GlobalInfo d_ene_Diffu;
  CutCell_GlobalInfo ene_postdiff_regular;
  CutCell_GlobalInfo ene_postdiff_dying;
  CutCell_GlobalInfo ene_postdiff_small;

  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();


  if (runge_kutta_fluxes_convection_)
    {
      cut_field_current_fluxes_conv[0].set_to_uniform_value(0);
      cut_field_current_fluxes_conv[1].set_to_uniform_value(0);
      cut_field_current_fluxes_conv[2].set_to_uniform_value(0);
      if (rk_step == 0)
        {
          cut_field_RK3_F_fluxes_conv[0].set_to_uniform_value(0);
          cut_field_RK3_F_fluxes_conv[1].set_to_uniform_value(0);
          cut_field_RK3_F_fluxes_conv[2].set_to_uniform_value(0);
        }
      cut_field_current_fluxes_conv[0].echange_espace_virtuel(1);
      cut_field_current_fluxes_conv[1].echange_espace_virtuel(1);
      cut_field_current_fluxes_conv[2].echange_espace_virtuel(1);
      cut_field_RK3_F_fluxes_conv[0].echange_espace_virtuel(1);
      cut_field_RK3_F_fluxes_conv[1].echange_espace_virtuel(1);
      cut_field_RK3_F_fluxes_conv[2].echange_espace_virtuel(1);
    }

  if (runge_kutta_fluxes_diffusion_)
    {
      cut_field_current_fluxes_diff[0].set_to_uniform_value(0);
      cut_field_current_fluxes_diff[1].set_to_uniform_value(0);
      cut_field_current_fluxes_diff[2].set_to_uniform_value(0);
      if (rk_step == 0)
        {
          cut_field_RK3_F_fluxes_diff[0].set_to_uniform_value(0);
          cut_field_RK3_F_fluxes_diff[1].set_to_uniform_value(0);
          cut_field_RK3_F_fluxes_diff[2].set_to_uniform_value(0);
        }
      cut_field_current_fluxes_diff[0].echange_espace_virtuel(1);
      cut_field_current_fluxes_diff[1].echange_espace_virtuel(1);
      cut_field_current_fluxes_diff[2].echange_espace_virtuel(1);
      cut_field_RK3_F_fluxes_diff[0].echange_espace_virtuel(1);
      cut_field_RK3_F_fluxes_diff[1].echange_espace_virtuel(1);
      cut_field_RK3_F_fluxes_diff[2].echange_espace_virtuel(1);
    }




  update_thermal_properties();

  const Cut_field_vector3_double& cut_field_total_velocity = ref_ijk_ft_cut_cell_->get_cut_field_total_velocity();
  //calculer_dT_cut_cell(cut_field_total_velocity); // Note : ne fait rien

  const double current_time = ref_ijk_ft_cut_cell_->get_current_time();

  // Toujours necessaire, ou presque
  diffusive_correction_.calculer_flux_interface_old(lambda_liquid_, lambda_vapour_, coord_facettes_, interfacial_temperature_, interfacial_phin_ai_, cut_field_temperature, ref_ijk_ft_cut_cell_, *temperature_, temperature_ft_);

  if (!conv_temperature_negligible_)
    {
      convective_correction_.calcule_temperature_remplissage(fractional_timestep, lambda_liquid_, lambda_vapour_, diffusive_correction_.flux_interface_ns_[0], interfacial_temperature_.centre, temperature_ft_, cut_field_total_velocity, cut_field_temperature);
    }


  ene_ini = compute_global_energy_cut_cell(cut_field_temperature, 0);

  if (verbosite_ >= 7)
    print_Tmin_Tmax_cut_cell(cut_field_temperature, 0, current_time, "Temperature Min/Max au debut du pas de temps");

  if (!diff_temperature_negligible_)
    {
      diffusive_correction_.compute_flux_dying_cells(cut_field_total_velocity, cut_field_temperature);
    }

  if (!conv_temperature_negligible_)
    {
      convective_correction_.compute_flux_dying_cells(cut_field_total_velocity, cut_field_temperature);
      convective_correction_.compute_flux_small_nascent_cells(cut_field_total_velocity, cut_field_temperature);
    }

  if (runge_kutta_fluxes_convection_)
    {
      set_rk_restreint(rk_step, runge_kutta_restriction_leniency_convection_, cut_cell_disc, cellule_rk_restreint_conv_v_, cellule_rk_restreint_conv_l_);
      ref_cast(OpConvQuickIJKScalar_cut_cell_double, temperature_convection_op_.valeur()).set_runge_kutta(rk_step, total_timestep, current_fluxes_conv_, cut_field_RK3_F_fluxes_conv, cellule_rk_restreint_conv_v_, cellule_rk_restreint_conv_l_);
    }

  compute_temperature_convection_cut_cell(cut_field_total_velocity, cut_field_d_temperature);
  d_ene_Conv = compute_d_global_energy_cut_cell(cut_field_d_temperature, 0);

  euler_explicit_update_cut_cell_transport(fractional_timestep, cut_field_d_temperature, cut_field_temperature);
  cut_field_temperature.echange_espace_virtuel(temperature_->ghost());

  if (postraiter_champs_intermediaires_)
    {
      cut_field_temperature_post_regular.copy_from(cut_field_temperature);
    }

  ene_postconv_switch = compute_global_energy_cut_cell(cut_field_temperature, 1);

  if (verbosite_ >= 7)
    print_Tmin_Tmax_cut_cell(cut_field_temperature, 1, current_time, "Temperature Min/Max apres la convection principale");

  if (!conv_temperature_negligible_)
    {
      if (runge_kutta_fluxes_convection_ && (!runge_kutta_fluxes_pas_de_correction_convection_))
        {
          cut_field_current_fluxes_conv[0].set_to_uniform_value(0);
          cut_field_current_fluxes_conv[1].set_to_uniform_value(0);
          cut_field_current_fluxes_conv[2].set_to_uniform_value(0);
        }

      convective_correction_.add_dying_cells(cut_field_total_velocity, cut_field_temperature, runge_kutta_fluxes_convection_, cut_field_current_fluxes_conv);
      cut_field_temperature.echange_espace_virtuel(temperature_->ghost());
    }

  if (postraiter_champs_intermediaires_)
    {
      cut_field_temperature_post_dying.copy_from(cut_field_temperature);
    }
  ene_postconv_dying = compute_global_energy_cut_cell(cut_field_temperature, 1);

  if (verbosite_ >= 7)
    print_Tmin_Tmax_cut_cell(cut_field_temperature, 1, current_time, "Temperature Min/Max apres la convection des cellules mourrantes");

  if (!conv_temperature_negligible_)
    {
      convective_correction_.add_small_nascent_cells(cut_field_total_velocity, cut_field_temperature, runge_kutta_fluxes_convection_, cut_field_current_fluxes_conv);
      cut_field_temperature.echange_espace_virtuel(temperature_->ghost());
    }

  if (postraiter_champs_intermediaires_)
    {
      cut_field_temperature_post_convection.copy_from(cut_field_temperature);
    }

  if (!diff_temperature_negligible_ && (!diffusive_correction_.deactivate_correction_petites_cellules_diffusion_))
    {
      diffusive_correction_.calcule_temperature_remplissage(fractional_timestep, lambda_liquid_, lambda_vapour_, diffusive_correction_.flux_interface_ns_[0], interfacial_temperature_.centre, temperature_ft_, cut_field_total_velocity, cut_field_temperature);
    }
  ene_postconv_small_nascent = compute_global_energy_cut_cell(cut_field_temperature, 1);

  if (verbosite_ >= 7)
    print_Tmin_Tmax_cut_cell(cut_field_temperature, 1, current_time, "Temperature Min/Max apres la convection des cellules naissantes");

  if (runge_kutta_fluxes_convection_ && (!runge_kutta_fluxes_pas_de_correction_convection_))
    {
      add_flux_times_vol_over_dt_surface(fractional_timestep, cut_field_current_fluxes_conv, cut_field_RK3_F_fluxes_conv);
    }


  if (!diff_temperature_negligible_ && (!diffusive_correction_.deactivate_correction_petites_cellules_diffusion_))
    {
      diffusive_correction_.compute_flux_small_nascent_cells(cut_field_total_velocity, cut_field_temperature);
    }

  if (!deactivate_diffusion_interface_)
    {
      diffusive_correction_.calculer_flux_interface_next(lambda_liquid_, lambda_vapour_, coord_facettes_, interfacial_temperature_, interfacial_phin_ai_, cut_field_temperature, ref_ijk_ft_cut_cell_, *temperature_, temperature_ft_);
      diffusive_correction_.calculer_flux_interface_efficace();
    }
  else
    {
      //Cerr << "Warning: The diffusion at the interface is not activated." << finl;
    }

  if (runge_kutta_fluxes_diffusion_)
    {
      set_rk_restreint(rk_step, runge_kutta_restriction_leniency_diffusion_, cut_cell_disc, cellule_rk_restreint_diff_v_, cellule_rk_restreint_diff_l_);
      ref_cast(OpDiffIJKScalar_cut_cell_double, temperature_diffusion_op_.valeur()).set_runge_kutta(rk_step, total_timestep, current_fluxes_diff_, cut_field_RK3_F_fluxes_diff, cellule_rk_restreint_diff_v_, cellule_rk_restreint_diff_l_);
    }

  add_temperature_diffusion();
  d_ene_Diffu = compute_d_global_energy_cut_cell(cut_field_d_temperature, 1);

  euler_explicit_update_cut_cell_notransport(fractional_timestep, true, cut_field_d_temperature, cut_field_temperature);

  cut_field_temperature.echange_espace_virtuel(temperature_->ghost());

  if (postraiter_champs_intermediaires_)
    {
      cut_field_temperature_post_diff_regular.copy_from(cut_field_temperature);
    }
  ene_postdiff_regular = compute_global_energy_cut_cell(cut_field_temperature, 1);

  if (verbosite_ >= 7)
    print_Tmin_Tmax_cut_cell(cut_field_temperature, 1, current_time, "Temperature Min/Max apres la diffusion principale");

  if (!diff_temperature_negligible_)
    {
      if (runge_kutta_fluxes_diffusion_ && (!runge_kutta_fluxes_pas_de_correction_diffusion_))
        {
          cut_field_current_fluxes_diff[0].set_to_uniform_value(0);
          cut_field_current_fluxes_diff[1].set_to_uniform_value(0);
          cut_field_current_fluxes_diff[2].set_to_uniform_value(0);
        }

      diffusive_correction_.add_dying_cells(cut_field_total_velocity, cut_field_temperature, runge_kutta_fluxes_diffusion_, cut_field_current_fluxes_diff);
      cut_field_temperature.echange_espace_virtuel(temperature_->ghost());
    }

  ene_postdiff_dying = compute_global_energy_cut_cell(cut_field_temperature, 1);

  if (verbosite_ >= 7)
    print_Tmin_Tmax_cut_cell(cut_field_temperature, 1, current_time, "Temperature Min/Max apres la diffusion des cellules mourrantes");

  if (!diff_temperature_negligible_ && (!diffusive_correction_.deactivate_correction_petites_cellules_diffusion_))
    {
      diffusive_correction_.add_small_nascent_cells(cut_field_total_velocity, cut_field_temperature, runge_kutta_fluxes_diffusion_, cut_field_current_fluxes_diff);
      cut_field_temperature.echange_espace_virtuel(temperature_->ghost());
    }

  ene_postdiff_small = compute_global_energy_cut_cell(cut_field_temperature, 1);

  Cerr.precision(16);
  if (verbosite_ >= 5)
    {
      Cerr << "dT_budget, overall time: " << current_time << " dT_conv: " << d_ene_Conv.overall << " dT_diff: " << d_ene_Diffu.overall << finl;
    }
  if (verbosite_ >= 6)
    {
      Cerr << "dT_budget, overall_l time: " << current_time << " dT_conv: " << d_ene_Conv.overall_l << " dT_diff: " << d_ene_Diffu.overall_l << finl;
      Cerr << "dT_budget, overall_v time: " << current_time << " dT_conv: " << d_ene_Conv.overall_v << " dT_diff: " << d_ene_Diffu.overall_v << finl;
      Cerr << "dT_budget, pure.......... time: " << current_time << " dT_conv: " << d_ene_Conv.pure         << " dT_diff: " << d_ene_Diffu.pure         << finl;
      Cerr << "dT_budget, diph_l........ time: " << current_time << " dT_conv: " << d_ene_Conv.diph_l       << " dT_diff: " << d_ene_Diffu.diph_l       << finl;
      Cerr << "dT_budget, diph_v........ time: " << current_time << " dT_conv: " << d_ene_Conv.diph_v       << " dT_diff: " << d_ene_Diffu.diph_v       << finl;
      Cerr << "dT_budget, diph_small.... time: " << current_time << " dT_conv: " << d_ene_Conv.diph_small   << " dT_diff: " << d_ene_Diffu.diph_small   << finl;
      Cerr << "dT_budget, diph_regular.. time: " << current_time << " dT_conv: " << d_ene_Conv.diph_regular << " dT_diff: " << d_ene_Diffu.diph_regular << finl;
      Cerr << "dT_budget, diph_nascent.. time: " << current_time << " dT_conv: " << d_ene_Conv.diph_nascent << " dT_diff: " << d_ene_Diffu.diph_nascent << finl;
      Cerr << "dT_budget, diph_dying.... time: " << current_time << " dT_conv: " << d_ene_Conv.diph_dying   << " dT_diff: " << d_ene_Diffu.diph_dying   << finl;
    }

  if (verbosite_ >= 1)
    {
      Cerr << "T_Budget, overall time: " << current_time << " initial: " << ene_ini.overall << " post_conv_dying: " << ene_postconv_dying.overall << " post_conv_switch: " << ene_postconv_switch.overall << " post_conv_small_nascent: " << ene_postconv_small_nascent.overall << " post_diff_regular: " << ene_postdiff_regular.overall << " post_diff_dying: " << ene_postdiff_dying.overall << " post_diff_small: " << ene_postdiff_small.overall << finl;
    }
  if (verbosite_ >= 3)
    {
      Cerr << "T_Budget, overall_l time: " << current_time << " initial: " << ene_ini.overall_l << " post_conv_switch: " << ene_postconv_switch.overall_l << " post_conv_dying: " << ene_postconv_dying.overall_l << " post_conv_small_nascent: " << ene_postconv_small_nascent.overall_l << " post_diff_regular: " << ene_postdiff_regular.overall_l << " post_diff_dying: " << ene_postdiff_dying.overall_l << " post_diff_small: " << ene_postdiff_small.overall_l << finl;
      Cerr << "T_Budget, overall_v time: " << current_time << " initial: " << ene_ini.overall_v << " post_conv_switch: " << ene_postconv_switch.overall_v << " post_conv_dying: " << ene_postconv_dying.overall_v << " post_conv_small_nascent: " << ene_postconv_small_nascent.overall_v << " post_diff_regular: " << ene_postdiff_regular.overall_v << " post_diff_dying: " << ene_postdiff_dying.overall_v << " post_diff_small: " << ene_postdiff_small.overall_v << finl;
      Cerr << "T_Budget, pure.......... time: " << current_time << " initial: " << ene_ini.pure         << " post_conv_switch: " << ene_postconv_switch.pure         << " post_conv_dying: " << ene_postconv_dying.pure         << " post_conv_small_nascent: " << ene_postconv_small_nascent.pure         << " post_diff_regular: " << ene_postdiff_regular.pure << " post_diff_dying: " << ene_postdiff_dying.pure         << " post_diff_small: " << ene_postdiff_small.pure         << finl;
      Cerr << "T_Budget, diph_l........ time: " << current_time << " initial: " << ene_ini.diph_l       << " post_conv_switch: " << ene_postconv_switch.diph_l       << " post_conv_dying: " << ene_postconv_dying.diph_l       << " post_conv_small_nascent: " << ene_postconv_small_nascent.diph_l       << " post_diff_regular: " << ene_postdiff_regular.diph_l << " post_diff_dying: " << ene_postdiff_dying.diph_l       << " post_diff_small: " << ene_postdiff_small.diph_l       << finl;
      Cerr << "T_Budget, diph_v........ time: " << current_time << " initial: " << ene_ini.diph_v       << " post_conv_switch: " << ene_postconv_switch.diph_v       << " post_conv_dying: " << ene_postconv_dying.diph_v       << " post_conv_small_nascent: " << ene_postconv_small_nascent.diph_v       << " post_diff_regular: " << ene_postdiff_regular.diph_v << " post_diff_dying: " << ene_postdiff_dying.diph_v       << " post_diff_small: " << ene_postdiff_small.diph_v       << finl;
      Cerr << "T_Budget, diph_small.... time: " << current_time << " initial: " << ene_ini.diph_small   << " post_conv_switch: " << ene_postconv_switch.diph_small   << " post_conv_dying: " << ene_postconv_dying.diph_small   << " post_conv_small_nascent: " << ene_postconv_small_nascent.diph_small   << " post_diff_regular: " << ene_postdiff_regular.diph_small << " post_diff_dying: " << ene_postdiff_dying.diph_small   << " post_diff_small: " << ene_postdiff_small.diph_small   << finl;
      Cerr << "T_Budget, diph_regular.. time: " << current_time << " initial: " << ene_ini.diph_regular << " post_conv_switch: " << ene_postconv_switch.diph_regular << " post_conv_dying: " << ene_postconv_dying.diph_regular << " post_conv_small_nascent: " << ene_postconv_small_nascent.diph_regular << " post_diff_regular: " << ene_postdiff_regular.diph_regular << " post_diff_dying: " << ene_postdiff_dying.diph_regular << " post_diff_small: " << ene_postdiff_small.diph_regular << finl;
      Cerr << "T_Budget, diph_nascent.. time: " << current_time << " initial: " << ene_ini.diph_nascent << " post_conv_switch: " << ene_postconv_switch.diph_nascent << " post_conv_dying: " << ene_postconv_dying.diph_nascent << " post_conv_small_nascent: " << ene_postconv_small_nascent.diph_nascent << " post_diff_regular: " << ene_postdiff_regular.diph_nascent << " post_diff_dying: " << ene_postdiff_dying.diph_nascent << " post_diff_small: " << ene_postdiff_small.diph_nascent << finl;
      Cerr << "T_Budget, diph_dying.... time: " << current_time << " initial: " << ene_ini.diph_dying   << " post_conv_switch: " << ene_postconv_switch.diph_dying   << " post_conv_dying: " << ene_postconv_dying.diph_dying   << " post_conv_small_nascent: " << ene_postconv_small_nascent.diph_dying   << " post_diff_regular: " << ene_postdiff_regular.diph_dying << " post_diff_dying: " << ene_postdiff_dying.diph_dying   << " post_diff_small: " << ene_postdiff_small.diph_dying   << finl;
    }

  print_Tmin_Tmax_cut_cell(cut_field_temperature, 1, current_time, "");

  if (runge_kutta_fluxes_diffusion_ && (!runge_kutta_fluxes_pas_de_correction_diffusion_))
    {
      add_flux_times_vol_over_dt_surface(fractional_timestep, cut_field_current_fluxes_diff, cut_field_RK3_F_fluxes_diff);
    }



}

void IJK_Thermal_cut_cell::euler_time_step(const double timestep)
{
  perform_thermal_step(timestep, 0, 0);
}

void IJK_Thermal_cut_cell::rk3_sub_step(const int rk_step, const double total_timestep,
                                        const double time)
{
  Cut_field_double& cut_field_temperature                = static_cast<Cut_field_double&>(*temperature_);

  if (debug_)
    Cerr << "Thermal Runge-Kutta3 time-step" << finl;

  perform_thermal_step(total_timestep, 1, rk_step);

  if (verbosite_ >= 7)
    print_Tmin_Tmax_cut_cell(cut_field_temperature, 1, time, "Temperature Min/Max apres etape Runge-Kutta (= etat final)");
}

void IJK_Thermal_cut_cell::sauvegarder_temperature(Nom& lata_name, int idx, const int& stop)
{
  fichier_reprise_temperature_ = lata_name;
  timestep_reprise_temperature_ = 1;
  dumplata_scalar_cut_cell(true, lata_name, Nom("TEMPERATURE_") + Nom(idx) , temperature_, 0 /*we store a 0 */);
  if (stop)
    latastep_reprise_ = latastep_reprise_ini_ + ref_ijk_ft_->get_tstep() + 1;
}


void IJK_Thermal_cut_cell::compute_diffusion_increment()
{
  const Cut_field_double& cut_field_div_coeff_grad_T_volume = static_cast<const Cut_field_double&>(*div_coeff_grad_T_volume_);
  Cut_field_double& cut_field_d_temperature                 = static_cast<Cut_field_double&>(*d_temperature_);

  // Update d_temperature
  const int ni = cut_field_d_temperature.ni();
  const int nj = cut_field_d_temperature.nj();
  const int nk = cut_field_d_temperature.nk();
  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              double rhocpV = rho_cp_(i,j,k) * vol_;

              const double ope = cut_field_div_coeff_grad_T_volume.pure_(i,j,k);
              const double resu = ope/rhocpV;
              cut_field_d_temperature.pure_(i,j,k) = resu ;
            }
        }
    }

  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_div_coeff_grad_T_volume.get_cut_cell_disc();
  const double rho_l = ref_ijk_ft_cut_cell_->get_rho_l();
  const double rho_v = ref_ijk_ft_cut_cell_->get_rho_v();
  for (int n = 0; n < cut_cell_disc.get_n_loc(); n++)
    {
      // On pre-divise par le volume cartesien vol_ et non par par le vrai
      // volume car on travaille en volume*temperature pour le theoreme de Reynolds.
      double rhocpV_l = rho_l * cp_liquid_ * vol_;
      double rhocpV_v = rho_v * cp_vapour_ * vol_;

      const double ope_l = cut_field_div_coeff_grad_T_volume.diph_l_(n);
      const double resu_l = ope_l/rhocpV_l;
      cut_field_d_temperature.diph_l_(n) = resu_l;

      const double ope_v = cut_field_div_coeff_grad_T_volume.diph_v_(n);
      const double resu_v = ope_v/rhocpV_v;
      cut_field_d_temperature.diph_v_(n) = resu_v;
    }
}

void cut_cell_reinit_streamObj(std::ostringstream& streamObj, const double& param)
{
  streamObj.str("");
  streamObj.clear();
  streamObj << (double) param;
}


void IJK_Thermal_cut_cell::calculer_dT_cut_cell(const Cut_field_vector3_double& cut_field_total_velocity)
{
  // Note : Cette fonction regroupe tous les elements de calculer_dT qui ne sont pas utilisees en cut cell

  /*
   * Clean_subproblems !
   */
  if (debug_)
    Cerr << "Clean thermal subproblems" << finl;
  clean_thermal_subproblems();

  // Correct the vapour and mixed cells values
  if (debug_)
    Cerr << "Store temperature before extrapolation" << finl;
  store_temperature_before_extrapolation();

  /*
   * Correct the temperature field using either the ghost-fluid
   * approach or the laminar sub-resolution approach (and zero values for debug)
   */
  if (debug_)
    Cerr << "Start the Ghost-fluid (GFM) approach" << finl;
  if (debug_)
    Cerr << "Br0 (GFM approach)" << finl;
  compute_eulerian_grad_T_interface();
  if (debug_)
    Cerr << "Br1 (GFM approach)" << finl;
  propagate_eulerian_grad_T_interface();
  if (debug_)
    Cerr << "Br2 (GFM approach)" << finl;
  compute_eulerian_temperature_ghost();
  if (debug_)
    Cerr << "Br3 (GFM approach)" << finl;
  compute_eulerian_bounding_box_fill_compo();
  if (debug_)
    Cerr << "End the Ghost-fluid (GFM) approach" << finl;

  /*
   * Compute gradients and hessian of the temperature after the ghost fluid extension
   */
  if (debug_)
    Cerr << "Compute temperature derivatives" << finl;
  compute_temperature_gradient_elem();
  compute_temperature_hessian_diag_elem();
  compute_temperature_hessian_cross_elem();

  /*
   * Compute sub-problems (For Subresolution Child classes only !)
   */
  if (debug_)
    Cerr << "Compute thermal subproblems" << finl;
  compute_thermal_subproblems();

  /*
   * Interpolate a value for the QUICK SCHEME (first call)
   */
  if (debug_)
    Cerr << "Compute temperature mixed cell centres" << finl;
  compute_temperature_cell_centres(0);

  /*
   * Convective and Diffusive fluxes
   */
  if (debug_)
    Cerr << "Compute thermal convective and diffusive fluxes from subproblems" << finl;
  compute_convective_diffusive_fluxes_face_centre();

  if (debug_)
    Cerr << "Prepare ij fluxes" << finl;
  if (!conv_temperature_negligible_ || !diff_temperature_negligible_)
    prepare_ij_fluxes_k_layers();


  double nb_diam_upstream_velocity = ref_ijk_ft_cut_cell_->get_nb_diam_upstream();
  if (nb_diam_upstream_ == 0.)
    nb_diam_upstream_ = nb_diam_upstream_velocity;
  if (upstream_temperature_ > -1e20 && ref_ijk_ft_cut_cell_->get_vitesse_upstream() > -1e20)
    force_upstream_temperature(*temperature_, upstream_temperature_,
                               ref_ijk_ft_cut_cell_->get_interface(), nb_diam_upstream_,
                               ref_ijk_ft_cut_cell_->get_upstream_dir(), ref_ijk_ft_cut_cell_->get_direction_gravite(),
                               ref_ijk_ft_cut_cell_->get_upstream_stencil());

  return;
}

void IJK_Thermal_cut_cell::compute_temperature_convection_cut_cell(const Cut_field_vector3_double& cut_field_total_velocity, Cut_field_double& cut_field_d_temperature)
{
  const Cut_field_double& cut_field_temperature = static_cast<const Cut_field_double&>(*temperature_);

  static Stat_Counter_Id cnt_conv_temp = statistiques().new_counter(1, "FT convection rho");
  statistiques().begin_count(cnt_conv_temp);
  if (conv_temperature_negligible_)
    {
      cut_field_d_temperature.set_to_uniform_value(0);
      u_T_convective_volume_.data() = 0;
    }
  else
    {
      temperature_convection_op_.calculer(cut_field_temperature,
                                          cut_field_total_velocity[0],
                                          cut_field_total_velocity[1],
                                          cut_field_total_velocity[2],
                                          cut_field_d_temperature);
      const int ni = cut_field_d_temperature.ni();
      const int nj = cut_field_d_temperature.nj();
      const int nk = cut_field_d_temperature.nk();
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  u_T_convective_volume_(i,j,k) = cut_field_d_temperature.pure_(i,j,k);
                  const double resu = cut_field_d_temperature.pure_(i,j,k) / vol_;
                  cut_field_d_temperature.pure_(i,j,k) = resu;
                  //if (liste_post_instantanes_.contient_("U_T_CONVECTIVE"))
                  //  {
                  //    u_T_convective_(i,j,k) = resu;
                  //  }
                }
            }
        }

      const Cut_cell_FT_Disc& cut_cell_disc = cut_field_d_temperature.get_cut_cell_disc();
      for (int n = 0; n < cut_cell_disc.get_n_loc(); n++)
        {
          // On pre-divise par le volume cartesien vol_ et non par par le vrai
          // volume car on travaille en volume*temperature pour le theoreme de Reynolds.
          double V_l = vol_;
          double V_v = vol_;

          const double resu_l = (V_l == 0) ? 0. : cut_field_d_temperature.diph_l_(n)/V_l;
          cut_field_d_temperature.diph_l_(n) = resu_l;

          const double resu_v = (V_v == 0) ? 0. : cut_field_d_temperature.diph_v_(n)/V_v;
          cut_field_d_temperature.diph_v_(n) = resu_v;
        }
    }
  statistiques().end_count(cnt_conv_temp);
  DebogIJK::verifier("op_conv(rho)", cut_field_d_temperature);
  return;
}


void IJK_Thermal_cut_cell::add_temperature_diffusion()
{
  const Cut_field_double& cut_field_temperature            = static_cast<const Cut_field_double&>(*temperature_);
  Cut_field_double& cut_field_div_coeff_grad_T_volume      = static_cast<Cut_field_double&>(*div_coeff_grad_T_volume_);
  Cut_field_double& cut_field_d_temperature                = static_cast<Cut_field_double&>(*d_temperature_);


  lambda_.echange_espace_virtuel(lambda_.ghost());
  DebogIJK::verifier("lambda", lambda_);

  temperature_diffusion_op_.set_lambda(lambda_);

  if (boundary_conditions_.get_bctype_k_min() == Boundary_Conditions_Thermique::Paroi_Temperature_imposee)
    {
      const double T_paroi_impose_kmin = boundary_conditions_.get_temperature_kmin();
      double lambda_de_t_paroi_kmin = lambda_liquid_;
      // calculer ici div(lambda*grad(T))*volume)
      calculer_flux_thermique_bord(*temperature_, lambda_de_t_paroi_kmin,
                                   T_paroi_impose_kmin, boundary_flux_kmin_, 0 /* boundary kmin */);
    }
  else if (boundary_conditions_.get_bctype_k_min() == Boundary_Conditions_Thermique::Paroi_Flux_impose)
    {
      const double flux_paroi_impose_kmin = boundary_conditions_.get_flux_kmin();
      imposer_flux_thermique_bord(*temperature_,
                                  flux_paroi_impose_kmin, boundary_flux_kmin_, 0 /* boundary kmin */);
    }
  else
    {
      // Perio... not needed. The flux on the "boundary" (in that case , it's not truely a boundary) will be computed as inside...
    }
  if (boundary_conditions_.get_bctype_k_max() == Boundary_Conditions_Thermique::Paroi_Temperature_imposee)
    {
      const double T_paroi_impose_kmax = boundary_conditions_.get_temperature_kmax();
      double lambda_de_t_paroi_kmax = lambda_liquid_;
      //TODO: calculer ici div(lambda*grad(T))*volume)
      calculer_flux_thermique_bord(*temperature_, lambda_de_t_paroi_kmax,
                                   T_paroi_impose_kmax, boundary_flux_kmax_, 1 /* boundary kmax */);
    }
  else if (boundary_conditions_.get_bctype_k_max() == Boundary_Conditions_Thermique::Paroi_Flux_impose)
    {
      const double flux_paroi_impose_kmax = boundary_conditions_.get_flux_kmax();
      imposer_flux_thermique_bord(*temperature_,
                                  flux_paroi_impose_kmax, boundary_flux_kmax_, 1 /* boundary kmax */);
      // Cerr << "not coded yet" << finl;
      // Process::exit();
    }
  else
    {
      // Perio... not needed. The flux on the "boundary" (in that case , it's not truely a boundary) will be computed as inside...
    }
  temperature_->echange_espace_virtuel(temperature_->ghost());
  DebogIJK::verifier("temp", *temperature_);

  if (diff_temperature_negligible_)
    {
      cut_field_d_temperature.set_to_uniform_value(0);
    }
  else
    {
      // Performance counters:
      static Stat_Counter_Id cnt_diff_temp = statistiques().new_counter(1, "FT diffusion temperature");
      statistiques().begin_count(cnt_diff_temp);
      /*
       * Correct the diffusive fluxes here or in the operator ?
       */
      temperature_diffusion_op_.calculer(cut_field_temperature,
                                         cut_field_div_coeff_grad_T_volume,
                                         boundary_flux_kmin_,
                                         boundary_flux_kmax_);
      diffusive_correction_.ajout_flux_interface_a_divergence_simple(cut_field_div_coeff_grad_T_volume);
      compute_diffusion_increment();
      statistiques().end_count(cnt_diff_temp);
      DebogIJK::verifier("div_coeff_grad_T_volume_", *div_coeff_grad_T_volume_);
    }
}

CutCell_GlobalInfo IJK_Thermal_cut_cell::compute_d_global_energy_cut_cell(Cut_field_double& cut_field_d_temperature, bool next)
{
  double global_energy_overall = 0.;
  double global_energy_overall_l = 0.;
  double global_energy_overall_v = 0.;
  double global_energy_pure = 0.;
  double global_energy_diph_l = 0.;
  double global_energy_diph_v = 0.;
  double global_energy_diph_small = 0.;
  double global_energy_diph_regular = 0.;
  double global_energy_diph_nascent = 0.;
  double global_energy_diph_dying = 0.;
  double count_overall = 0.;
  double count_overall_l = 0.;
  double count_overall_v = 0.;
  double count_pure = 0.;
  double count_diph_l = 0.;
  double count_diph_v = 0.;
  double count_diph_small = 0.;
  double count_diph_regular = 0.;
  double count_diph_nascent = 0.;
  double count_diph_dying = 0.;
  const IJK_Field_double& indic_old = ref_ijk_ft_cut_cell_->itfce().I();
  const IJK_Field_double& indic_next = ref_ijk_ft_cut_cell_->itfce().In();
  const int nx = cut_field_d_temperature.ni();
  const int ny = cut_field_d_temperature.nj();
  const int nz = cut_field_d_temperature.nk();
  const IJK_Grid_Geometry& geom = indic_next.get_splitting().get_grid_geometry();
  assert(geom.get_constant_delta(DIRECTION_K) >0); // To be sure we're on a regular mesh
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_d_temperature.get_cut_cell_disc();
  for (int k=0; k < nz ; k++)
    {
      for (int j=0; j< ny; j++)
        {
          for (int i=0; i < nx; i++)
            {
              int n = cut_cell_disc.get_n(i,j,k);
              if (n < 0)
                {
                  global_energy_overall += cut_field_d_temperature.pure_(i,j,k);
                  count_overall += 1;

                  global_energy_overall_l += cut_field_d_temperature.pure_(i,j,k);
                  count_overall_l += 1;

                  global_energy_overall_v += cut_field_d_temperature.pure_(i,j,k);
                  count_overall_v += 1;

                  global_energy_pure += cut_field_d_temperature.pure_(i,j,k);
                  count_pure += 1;
                }
              else
                {
                  double chi_T_l = cut_field_d_temperature.diph_l_(n);
                  double chi_T_v = cut_field_d_temperature.diph_v_(n);

                  global_energy_overall += chi_T_l;
                  global_energy_overall += chi_T_v;
                  count_overall += 1;

                  global_energy_overall_l += chi_T_l;
                  count_overall_l += 1;

                  global_energy_overall_v += chi_T_v;
                  count_overall_v += 1;

                  global_energy_diph_l += chi_T_l;
                  count_diph_l += 1;

                  global_energy_diph_v += chi_T_v;
                  count_diph_v += 1;

                  if (ref_ijk_ft_cut_cell_->itfce().devient_pure(i,j,k))
                    {
                      int phase_dying = (int)(1 - indic_next(i,j,k));
                      if (phase_dying == 1)
                        {
                          global_energy_diph_dying += chi_T_l;
                          count_diph_dying += 1;
                        }
                      else
                        {
                          global_energy_diph_dying += chi_T_v;
                          count_diph_dying += 1;
                        }
                    }
                  else if (ref_ijk_ft_cut_cell_->itfce().devient_diphasique(i,j,k))
                    {
                      int phase_nascent = (int)(1 - indic_old(i,j,k));
                      if (phase_nascent == 1)
                        {
                          global_energy_diph_nascent += chi_T_l;
                          count_diph_nascent += 1;
                        }
                      else
                        {
                          global_energy_diph_nascent += chi_T_v;
                          count_diph_nascent += 1;
                        }
                    }
                  else if (ref_ijk_ft_cut_cell_->itfce().next_below_small_threshold_for_phase(1, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      global_energy_diph_small += chi_T_l;
                      count_diph_small += 1;

                      global_energy_diph_regular += chi_T_v;
                      count_diph_regular += 1;
                    }
                  else if (ref_ijk_ft_cut_cell_->itfce().next_below_small_threshold_for_phase(0, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      global_energy_diph_small += chi_T_v;
                      count_diph_small += 1;

                      global_energy_diph_regular += chi_T_l;
                      count_diph_regular += 1;
                    }
                  else
                    {
                      global_energy_diph_regular += chi_T_l;
                      global_energy_diph_regular += chi_T_v;
                      count_diph_regular += 1;
                    }

                }
            }
        }
    }
  count_overall = mp_sum(count_overall);
  count_overall_l = mp_sum(count_overall_l);
  count_overall_v = mp_sum(count_overall_v);
  count_pure = mp_sum(count_pure);
  count_diph_l = mp_sum(count_diph_l);
  count_diph_v = mp_sum(count_diph_v);
  count_diph_small = mp_sum(count_diph_small);
  count_diph_regular = mp_sum(count_diph_regular);
  count_diph_nascent = mp_sum(count_diph_nascent);
  count_diph_dying = mp_sum(count_diph_dying);
  assert(count_overall == cut_field_d_temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_I)
         *cut_field_d_temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_J)
         *cut_field_d_temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K));
  const double vol_cell = geom.get_constant_delta(DIRECTION_I)*geom.get_constant_delta(DIRECTION_J)*geom.get_constant_delta(DIRECTION_K);
  global_energy_overall      = vol_cell * mp_sum(global_energy_overall);
  global_energy_overall_l    = vol_cell * mp_sum(global_energy_overall_l);
  global_energy_overall_v    = vol_cell * mp_sum(global_energy_overall_v);
  global_energy_pure         = vol_cell * mp_sum(global_energy_pure);
  global_energy_diph_l       = vol_cell * mp_sum(global_energy_diph_l);
  global_energy_diph_v       = vol_cell * mp_sum(global_energy_diph_v);
  global_energy_diph_small   = vol_cell * mp_sum(global_energy_diph_small);
  global_energy_diph_regular = vol_cell * mp_sum(global_energy_diph_regular);
  global_energy_diph_nascent = vol_cell * mp_sum(global_energy_diph_nascent);
  global_energy_diph_dying   = vol_cell * mp_sum(global_energy_diph_dying);
  return {global_energy_overall, global_energy_overall_l, global_energy_overall_v, global_energy_pure, global_energy_diph_l, global_energy_diph_v, global_energy_diph_small, global_energy_diph_regular, global_energy_diph_nascent, global_energy_diph_dying};
}

CutCell_GlobalInfo IJK_Thermal_cut_cell::compute_global_energy_cut_cell(Cut_field_double& cut_field_temperature, bool next)
{
  double global_energy_overall = 0.;
  double global_energy_overall_l = 0.;
  double global_energy_overall_v = 0.;
  double global_energy_pure = 0.;
  double global_energy_diph_l = 0.;
  double global_energy_diph_v = 0.;
  double global_energy_diph_small = 0.;
  double global_energy_diph_regular = 0.;
  double global_energy_diph_nascent = 0.;
  double global_energy_diph_dying = 0.;
  double count_overall = 0.;
  double count_overall_l = 0.;
  double count_overall_v = 0.;
  double count_pure = 0.;
  double count_diph_l = 0.;
  double count_diph_v = 0.;
  double count_diph_small = 0.;
  double count_diph_regular = 0.;
  double count_diph_nascent = 0.;
  double count_diph_dying = 0.;
  const IJK_Field_double& indic_old = ref_ijk_ft_cut_cell_->itfce().I();
  const IJK_Field_double& indic_next = ref_ijk_ft_cut_cell_->itfce().In();
  const double rhocpl = get_rhocp_l();
  const double rhocpv = get_rhocp_v();
  const int nx = cut_field_temperature.ni();
  const int ny = cut_field_temperature.nj();
  const int nz = cut_field_temperature.nk();
  const IJK_Grid_Geometry& geom = indic_next.get_splitting().get_grid_geometry();
  assert(geom.get_constant_delta(DIRECTION_K) >0); // To be sure we're on a regular mesh
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();
  for (int k=0; k < nz ; k++)
    {
      for (int j=0; j< ny; j++)
        {
          for (int i=0; i < nx; i++)
            {
              double chi_l = next ? ref_ijk_ft_cut_cell_->itfce().In_nonzero(1,i,j,k) : ref_ijk_ft_cut_cell_->itfce().I_nonzero(1,i,j,k);
              double chi_v = next ? ref_ijk_ft_cut_cell_->itfce().In_nonzero(0,i,j,k) : ref_ijk_ft_cut_cell_->itfce().I_nonzero(0,i,j,k);
              int n = cut_cell_disc.get_n(i,j,k);
              if (n < 0)
                {
                  global_energy_overall += (chi_l * rhocpl + chi_v * rhocpv) * cut_field_temperature.pure_(i,j,k);
                  count_overall += 1;

                  global_energy_overall_l += chi_l * rhocpl * cut_field_temperature.pure_(i,j,k);
                  count_overall_l += chi_l;

                  global_energy_overall_v += chi_v * rhocpv * cut_field_temperature.pure_(i,j,k);
                  count_overall_v += chi_v;

                  global_energy_pure += (chi_l * rhocpl + chi_v * rhocpv) * cut_field_temperature.pure_(i,j,k);
                  count_pure += 1;
                }
              else
                {
                  double chi_T_l = chi_l * cut_field_temperature.diph_l_(n);
                  double chi_T_v = chi_v * cut_field_temperature.diph_v_(n);

                  global_energy_overall += rhocpl * chi_T_l;
                  global_energy_overall += rhocpv * chi_T_v;
                  count_overall += 1;

                  global_energy_overall_l += rhocpl * chi_T_l;
                  count_overall_l += 1;

                  global_energy_overall_v += rhocpv * chi_T_v;
                  count_overall_v += 1;

                  global_energy_diph_l += rhocpl * chi_T_l;
                  count_diph_l += chi_l;

                  global_energy_diph_v += rhocpv * chi_T_v;
                  count_diph_v += chi_v;

                  if (ref_ijk_ft_cut_cell_->itfce().devient_pure(i,j,k))
                    {
                      int phase_dying = (int)(1 - indic_next(i,j,k));
                      if (phase_dying == 1)
                        {
                          global_energy_diph_dying += rhocpl * chi_T_l;
                          count_diph_dying += chi_l;
                        }
                      else
                        {
                          global_energy_diph_dying += rhocpv * chi_T_v;
                          count_diph_dying += chi_v;
                        }
                    }
                  else if (ref_ijk_ft_cut_cell_->itfce().devient_diphasique(i,j,k))
                    {
                      int phase_nascent = (int)(1 - indic_old(i,j,k));
                      if (phase_nascent == 1)
                        {
                          global_energy_diph_nascent += rhocpl * chi_T_l;
                          count_diph_nascent += chi_l;
                        }
                      else
                        {
                          global_energy_diph_nascent += rhocpv * chi_T_v;
                          count_diph_nascent += chi_v;
                        }
                    }
                  else if (ref_ijk_ft_cut_cell_->itfce().next_below_small_threshold_for_phase(1, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      global_energy_diph_small += rhocpl * chi_T_l;
                      count_diph_small += chi_l;

                      global_energy_diph_regular += rhocpv * chi_T_v;
                      count_diph_regular += chi_v;
                    }
                  else if (ref_ijk_ft_cut_cell_->itfce().next_below_small_threshold_for_phase(0, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      global_energy_diph_small += rhocpv * chi_T_v;
                      count_diph_small += chi_v;

                      global_energy_diph_regular += rhocpl * chi_T_l;
                      count_diph_regular += chi_l;
                    }
                  else
                    {
                      global_energy_diph_regular += rhocpl * chi_T_l;
                      global_energy_diph_regular += rhocpv * chi_T_v;
                      count_diph_regular += 1;
                    }

                }
            }
        }
    }
  count_overall = mp_sum(count_overall);
  count_overall_l = mp_sum(count_overall_l);
  count_overall_v = mp_sum(count_overall_v);
  count_pure = mp_sum(count_pure);
  count_diph_l = mp_sum(count_diph_l);
  count_diph_v = mp_sum(count_diph_v);
  count_diph_small = mp_sum(count_diph_small);
  count_diph_regular = mp_sum(count_diph_regular);
  count_diph_nascent = mp_sum(count_diph_nascent);
  count_diph_dying = mp_sum(count_diph_dying);
  assert(count_overall == cut_field_temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_I)
         *cut_field_temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_J)
         *cut_field_temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K));
  const double vol_cell = geom.get_constant_delta(DIRECTION_I)*geom.get_constant_delta(DIRECTION_J)*geom.get_constant_delta(DIRECTION_K);
  global_energy_overall      = vol_cell * mp_sum(global_energy_overall);
  global_energy_overall_l    = vol_cell * mp_sum(global_energy_overall_l);
  global_energy_overall_v    = vol_cell * mp_sum(global_energy_overall_v);
  global_energy_pure         = vol_cell * mp_sum(global_energy_pure);
  global_energy_diph_l       = vol_cell * mp_sum(global_energy_diph_l);
  global_energy_diph_v       = vol_cell * mp_sum(global_energy_diph_v);
  global_energy_diph_small   = vol_cell * mp_sum(global_energy_diph_small);
  global_energy_diph_regular = vol_cell * mp_sum(global_energy_diph_regular);
  global_energy_diph_nascent = vol_cell * mp_sum(global_energy_diph_nascent);
  global_energy_diph_dying   = vol_cell * mp_sum(global_energy_diph_dying);
  return {global_energy_overall, global_energy_overall_l, global_energy_overall_v, global_energy_pure, global_energy_diph_l, global_energy_diph_v, global_energy_diph_small, global_energy_diph_regular, global_energy_diph_nascent, global_energy_diph_dying};
}

Nom get_T_location(double T, const Cut_field_double& cut_field_temperature)
{
  const int nx = cut_field_temperature.ni();
  const int ny = cut_field_temperature.nj();
  const int nz = cut_field_temperature.nk();
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();
  for (int k=0; k < nz ; k++)
    {
      for (int j=0; j< ny; j++)
        {
          for (int i=0; i < nx; i++)
            {
              int n = cut_cell_disc.get_n(i,j,k);
              if (n < 0)
                {
                  if (T == cut_field_temperature.pure_(i,j,k))
                    {
                      return Nom("ijk_pure_") + Nom(i) + Nom("_") + Nom(j) + Nom("_") + Nom(k);
                    }
                }
              else
                {
                  if (T == cut_field_temperature.diph_l_(n))
                    {
                      return Nom("ijk_l_") + Nom(i) + Nom("_") + Nom(j) + Nom("_") + Nom(k);
                    }
                  if (T == cut_field_temperature.diph_v_(n))
                    {
                      return Nom("ijk_v_") + Nom(i) + Nom("_") + Nom(j) + Nom("_") + Nom(k);
                    }
                }
            }
        }
    }
  return Nom("_not_here_");
}

void IJK_Thermal_cut_cell::print_Tmin_Tmax_cut_cell(const Cut_field_double& cut_field_temperature, bool next, double current_time, const std::string& heading)
{
  CutCell_GlobalInfo Tmin = compute_Tmin_cut_cell(cut_field_temperature, next);
  CutCell_GlobalInfo Tmax = compute_Tmax_cut_cell(cut_field_temperature, next);
  if (verbosite_ >= 7)
    {
      Cerr << heading << finl;
    }
  if (verbosite_ >= 2)
    {
      Cerr << "T_MinMax, overall time: " << current_time << " Tmin: " << Tmin.overall << " Tmax: " << Tmax.overall << finl;
    }
  if (verbosite_ >= 4)
    {
      Cerr << "T_MinMax, overall_l time: " << current_time << " Tmin: " << Tmin.overall_l << " Tmax: " << Tmax.overall_l << finl;
      Cerr << "T_MinMax, overall_v time: " << current_time << " Tmin: " << Tmin.overall_v << " Tmax: " << Tmax.overall_v << finl;
      Cerr << "T_MinMax, pure.......... time: " << current_time << " Tmin: " << Tmin.pure         << " Tmax: " << Tmax.pure         << finl;
      Cerr << "T_MinMax, diph_l........ time: " << current_time << " Tmin: " << Tmin.diph_l       << " Tmax: " << Tmax.diph_l       << finl;
      Cerr << "T_MinMax, diph_v........ time: " << current_time << " Tmin: " << Tmin.diph_v       << " Tmax: " << Tmax.diph_v       << finl;
      Cerr << "T_MinMax, diph_small.... time: " << current_time << " Tmin: " << Tmin.diph_small   << " Tmax: " << Tmax.diph_small   << finl;
      Cerr << "T_MinMax, diph_regular.. time: " << current_time << " Tmin: " << Tmin.diph_regular << " Tmax: " << Tmax.diph_regular << finl;
      Cerr << "T_MinMax, diph_nascent.. time: " << current_time << " Tmin: " << Tmin.diph_nascent << " Tmax: " << Tmax.diph_nascent << finl;
      Cerr << "T_MinMax, diph_dying.... time: " << current_time << " Tmin: " << Tmin.diph_dying   << " Tmax: " << Tmax.diph_dying   << finl;
    }
  if (verbosite_ >= 10)
    {
      Cerr << "Location T_MinMax, overall_l time: " << current_time << " Tmin: " << get_T_location(Tmin.overall_l, cut_field_temperature) << " Tmax: " << get_T_location(Tmax.overall_l, cut_field_temperature) << finl;
      Cerr << "Location T_MinMax, overall_v time: " << current_time << " Tmin: " << get_T_location(Tmin.overall_v, cut_field_temperature) << " Tmax: " << get_T_location(Tmax.overall_v, cut_field_temperature) << finl;
      Cerr << "Location T_MinMax, pure.......... time: " << current_time << " Tmin: " << get_T_location(Tmin.pure, cut_field_temperature)         << " Tmax: " << get_T_location(Tmax.pure, cut_field_temperature)         << finl;
      Cerr << "Location T_MinMax, diph_l........ time: " << current_time << " Tmin: " << get_T_location(Tmin.diph_l, cut_field_temperature)       << " Tmax: " << get_T_location(Tmax.diph_l, cut_field_temperature)       << finl;
      Cerr << "Location T_MinMax, diph_v........ time: " << current_time << " Tmin: " << get_T_location(Tmin.diph_v, cut_field_temperature)       << " Tmax: " << get_T_location(Tmax.diph_v, cut_field_temperature)       << finl;
      Cerr << "Location T_MinMax, diph_small.... time: " << current_time << " Tmin: " << get_T_location(Tmin.diph_small, cut_field_temperature)   << " Tmax: " << get_T_location(Tmax.diph_small, cut_field_temperature)   << finl;
      Cerr << "Location T_MinMax, diph_regular.. time: " << current_time << " Tmin: " << get_T_location(Tmin.diph_regular, cut_field_temperature) << " Tmax: " << get_T_location(Tmax.diph_regular, cut_field_temperature) << finl;
      Cerr << "Location T_MinMax, diph_nascent.. time: " << current_time << " Tmin: " << get_T_location(Tmin.diph_nascent, cut_field_temperature) << " Tmax: " << get_T_location(Tmax.diph_nascent, cut_field_temperature) << finl;
      Cerr << "Location T_MinMax, diph_dying.... time: " << current_time << " Tmin: " << get_T_location(Tmin.diph_dying, cut_field_temperature)   << " Tmax: " << get_T_location(Tmax.diph_dying, cut_field_temperature)   << finl;
    }
}

CutCell_GlobalInfo IJK_Thermal_cut_cell::compute_Tmin_cut_cell(const Cut_field_double& cut_field_temperature, bool next)
{
  double Tmin_overall = 1.e20;
  double Tmin_overall_l = 1.e20;
  double Tmin_overall_v = 1.e20;
  double Tmin_pure = 1.e20;
  double Tmin_diph_l = 1.e20;
  double Tmin_diph_v = 1.e20;
  double Tmin_diph_small = 1.e20;
  double Tmin_diph_regular = 1.e20;
  double Tmin_diph_nascent = 1.e20;
  double Tmin_diph_dying = 1.e20;
  const IJK_Field_double& indic_old = ref_ijk_ft_cut_cell_->itfce().I();
  const IJK_Field_double& indic_next = ref_ijk_ft_cut_cell_->itfce().In();
  const IJK_Field_double& indic = next ? indic_next : indic_old;
  const int nx = cut_field_temperature.ni();
  const int ny = cut_field_temperature.nj();
  const int nz = cut_field_temperature.nk();
  assert(indic.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_K) >0); // To be sure we're on a regular mesh
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();
  for (int k=0; k < nz ; k++)
    {
      for (int j=0; j< ny; j++)
        {
          for (int i=0; i < nx; i++)
            {
              double chi_l = indic(i,j,k);
              int n = cut_cell_disc.get_n(i,j,k);
              if (n < 0)
                {
                  Tmin_overall = std::min(Tmin_overall, cut_field_temperature.pure_(i,j,k));
                  Tmin_overall_l = (chi_l == 0) ? Tmin_overall_l : std::min(Tmin_overall_l, cut_field_temperature.pure_(i,j,k));
                  Tmin_overall_v = (chi_l == 1) ? Tmin_overall_v : std::min(Tmin_overall_v, cut_field_temperature.pure_(i,j,k));
                  Tmin_pure = std::min(Tmin_pure, cut_field_temperature.pure_(i,j,k));
                }
              else
                {
                  // Excluding the value of the phase in dying cells
                  bool exclude_l = (next && (ref_ijk_ft_cut_cell_->itfce().phase_mourrante(1, i,j,k))) || ((!next) && (ref_ijk_ft_cut_cell_->itfce().phase_naissante(1, i,j,k)));
                  bool exclude_v = (next && (ref_ijk_ft_cut_cell_->itfce().phase_mourrante(0, i,j,k))) || ((!next) && (ref_ijk_ft_cut_cell_->itfce().phase_naissante(0, i,j,k)));

                  Tmin_overall = exclude_l ? Tmin_overall : std::min(Tmin_overall, cut_field_temperature.diph_l_(n));
                  Tmin_overall = exclude_v ? Tmin_overall : std::min(Tmin_overall, cut_field_temperature.diph_v_(n));

                  Tmin_overall_l = exclude_l ? Tmin_overall_l : std::min(Tmin_overall_l, cut_field_temperature.diph_l_(n));
                  Tmin_overall_v = exclude_v ? Tmin_overall_v : std::min(Tmin_overall_v, cut_field_temperature.diph_v_(n));

                  Tmin_diph_l = exclude_l ? Tmin_diph_l : std::min(Tmin_diph_l, cut_field_temperature.diph_l_(n));
                  Tmin_diph_v = exclude_v ? Tmin_diph_v : std::min(Tmin_diph_v, cut_field_temperature.diph_v_(n));

                  if (ref_ijk_ft_cut_cell_->itfce().next_below_small_threshold_for_phase(1, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      Tmin_diph_small = exclude_l ? Tmin_diph_small : std::min(Tmin_diph_small, cut_field_temperature.diph_l_(n));
                      Tmin_diph_regular = exclude_v ? Tmin_diph_regular : std::min(Tmin_diph_regular, cut_field_temperature.diph_v_(n));
                    }
                  else if (ref_ijk_ft_cut_cell_->itfce().next_below_small_threshold_for_phase(0, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      Tmin_diph_small = exclude_v ? Tmin_diph_small : std::min(Tmin_diph_small, cut_field_temperature.diph_v_(n));
                      Tmin_diph_regular = exclude_l ? Tmin_diph_regular : std::min(Tmin_diph_regular, cut_field_temperature.diph_l_(n));
                    }
                  else
                    {
                      Tmin_diph_regular = exclude_l ? Tmin_diph_regular : std::min(Tmin_diph_regular, cut_field_temperature.diph_l_(n));
                      Tmin_diph_regular = exclude_v ? Tmin_diph_regular : std::min(Tmin_diph_regular, cut_field_temperature.diph_v_(n));
                    }

                  if (ref_ijk_ft_cut_cell_->itfce().devient_pure(i,j,k))
                    {
                      int phase_dying = (int)(1 - indic_next(i,j,k));
                      if (phase_dying == 1)
                        {
                          Tmin_diph_dying = std::min(Tmin_diph_dying, cut_field_temperature.diph_l_(n));
                        }
                      else
                        {
                          Tmin_diph_dying = std::min(Tmin_diph_dying, cut_field_temperature.diph_v_(n));
                        }
                    }

                  if (ref_ijk_ft_cut_cell_->itfce().devient_diphasique(i,j,k))
                    {
                      int phase_nascent = (int)(1 - indic_old(i,j,k));
                      if (phase_nascent == 1)
                        {
                          Tmin_diph_nascent = std::min(Tmin_diph_nascent, cut_field_temperature.diph_l_(n));
                        }
                      else
                        {
                          Tmin_diph_nascent = std::min(Tmin_diph_nascent, cut_field_temperature.diph_v_(n));
                        }
                    }
                }
            }
        }
    }
  Tmin_overall      = Process::mp_min(Tmin_overall);
  Tmin_overall_l      = Process::mp_min(Tmin_overall_l);
  Tmin_overall_v      = Process::mp_min(Tmin_overall_v);
  Tmin_pure         = Process::mp_min(Tmin_pure);
  Tmin_diph_l       = Process::mp_min(Tmin_diph_l);
  Tmin_diph_v       = Process::mp_min(Tmin_diph_v);
  Tmin_diph_small   = Process::mp_min(Tmin_diph_small);
  Tmin_diph_regular   = Process::mp_min(Tmin_diph_regular);
  Tmin_diph_nascent = Process::mp_min(Tmin_diph_nascent);
  Tmin_diph_dying   = Process::mp_min(Tmin_diph_dying);
  return {Tmin_overall, Tmin_overall_l, Tmin_overall_v, Tmin_pure, Tmin_diph_l, Tmin_diph_v, Tmin_diph_small, Tmin_diph_regular, Tmin_diph_nascent, Tmin_diph_dying};
}

CutCell_GlobalInfo IJK_Thermal_cut_cell::compute_Tmax_cut_cell(const Cut_field_double& cut_field_temperature, bool next)
{
  double Tmax_overall = -1.e20;
  double Tmax_overall_l = -1.e20;
  double Tmax_overall_v = -1.e20;
  double Tmax_pure = -1.e20;
  double Tmax_diph_l = -1.e20;
  double Tmax_diph_v = -1.e20;
  double Tmax_diph_small = -1.e20;
  double Tmax_diph_regular = -1.e20;
  double Tmax_diph_nascent = -1.e20;
  double Tmax_diph_dying = -1.e20;
  const IJK_Field_double& indic_old = ref_ijk_ft_cut_cell_->itfce().I();
  const IJK_Field_double& indic_next = ref_ijk_ft_cut_cell_->itfce().In();
  const IJK_Field_double& indic = next ? indic_next : indic_old;
  const int nx = cut_field_temperature.ni();
  const int ny = cut_field_temperature.nj();
  const int nz = cut_field_temperature.nk();
  // To be sure we're on a regular mesh
  assert(indic.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_K) >0);
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();
  for (int k=0; k < nz ; k++)
    {
      for (int j=0; j< ny; j++)
        {
          for (int i=0; i < nx; i++)
            {
              double chi_l = indic(i,j,k);
              int n = cut_cell_disc.get_n(i,j,k);
              if (n < 0)
                {
                  Tmax_overall = std::max(Tmax_overall, cut_field_temperature.pure_(i,j,k));
                  Tmax_overall_l = (chi_l == 0) ? Tmax_overall_l : std::max(Tmax_overall_l, cut_field_temperature.pure_(i,j,k));
                  Tmax_overall_v = (chi_l == 1) ? Tmax_overall_v : std::max(Tmax_overall_v, cut_field_temperature.pure_(i,j,k));
                  Tmax_pure = std::max(Tmax_pure, cut_field_temperature.pure_(i,j,k));
                }
              else
                {
                  // Excluding the value of the phase in dying cells
                  bool exclude_l = (next && (ref_ijk_ft_cut_cell_->itfce().phase_mourrante(1, i,j,k))) || ((!next) && (ref_ijk_ft_cut_cell_->itfce().phase_naissante(1, i,j,k)));
                  bool exclude_v = (next && (ref_ijk_ft_cut_cell_->itfce().phase_mourrante(0, i,j,k))) || ((!next) && (ref_ijk_ft_cut_cell_->itfce().phase_naissante(0, i,j,k)));

                  Tmax_overall = exclude_l ? Tmax_overall : std::max(Tmax_overall, cut_field_temperature.diph_l_(n));
                  Tmax_overall = exclude_v ? Tmax_overall : std::max(Tmax_overall, cut_field_temperature.diph_v_(n));

                  Tmax_overall_l = exclude_l ? Tmax_overall_l : std::max(Tmax_overall_l, cut_field_temperature.diph_l_(n));
                  Tmax_overall_v = exclude_v ? Tmax_overall_v : std::max(Tmax_overall_v, cut_field_temperature.diph_v_(n));

                  Tmax_diph_l = exclude_l ? Tmax_diph_l : std::max(Tmax_diph_l, cut_field_temperature.diph_l_(n));
                  Tmax_diph_v = exclude_v ? Tmax_diph_v : std::max(Tmax_diph_v, cut_field_temperature.diph_v_(n));

                  if (ref_ijk_ft_cut_cell_->itfce().next_below_small_threshold_for_phase(1, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      Tmax_diph_small = exclude_l ? Tmax_diph_small : std::max(Tmax_diph_small, cut_field_temperature.diph_l_(n));
                      Tmax_diph_regular = exclude_v ? Tmax_diph_regular : std::max(Tmax_diph_regular, cut_field_temperature.diph_v_(n));
                    }
                  else if (ref_ijk_ft_cut_cell_->itfce().next_below_small_threshold_for_phase(0, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      Tmax_diph_small = exclude_v ? Tmax_diph_small : std::max(Tmax_diph_small, cut_field_temperature.diph_v_(n));
                      Tmax_diph_regular = exclude_l ? Tmax_diph_regular : std::max(Tmax_diph_regular, cut_field_temperature.diph_l_(n));
                    }
                  else
                    {
                      Tmax_diph_regular = exclude_l ? Tmax_diph_regular : std::max(Tmax_diph_regular, cut_field_temperature.diph_l_(n));
                      Tmax_diph_regular = exclude_v ? Tmax_diph_regular : std::max(Tmax_diph_regular, cut_field_temperature.diph_v_(n));
                    }

                  if (ref_ijk_ft_cut_cell_->itfce().devient_pure(i,j,k))
                    {
                      int phase_dying = (int)(1 - indic_next(i,j,k));
                      if (phase_dying == 1)
                        {
                          Tmax_diph_dying = std::max(Tmax_diph_dying, cut_field_temperature.diph_l_(n));
                        }
                      else
                        {
                          Tmax_diph_dying = std::max(Tmax_diph_dying, cut_field_temperature.diph_v_(n));
                        }
                    }

                  if (ref_ijk_ft_cut_cell_->itfce().devient_diphasique(i,j,k))
                    {
                      int phase_nascent = (int)(1 - indic_old(i,j,k));
                      if (phase_nascent == 1)
                        {
                          Tmax_diph_nascent = std::max(Tmax_diph_nascent, cut_field_temperature.diph_l_(n));
                        }
                      else
                        {
                          Tmax_diph_nascent = std::max(Tmax_diph_nascent, cut_field_temperature.diph_v_(n));
                        }
                    }
                }
            }
        }
    }
  Tmax_overall      = Process::mp_max(Tmax_overall);
  Tmax_overall_l      = Process::mp_max(Tmax_overall_l);
  Tmax_overall_v      = Process::mp_max(Tmax_overall_v);
  Tmax_pure         = Process::mp_max(Tmax_pure);
  Tmax_diph_l       = Process::mp_max(Tmax_diph_l);
  Tmax_diph_v       = Process::mp_max(Tmax_diph_v);
  Tmax_diph_small   = Process::mp_max(Tmax_diph_small);
  Tmax_diph_regular   = Process::mp_max(Tmax_diph_regular);
  Tmax_diph_nascent = Process::mp_max(Tmax_diph_nascent);
  Tmax_diph_dying   = Process::mp_max(Tmax_diph_dying);
  return {Tmax_overall, Tmax_overall_l, Tmax_overall_v, Tmax_pure, Tmax_diph_l, Tmax_diph_v, Tmax_diph_small, Tmax_diph_regular, Tmax_diph_nascent, Tmax_diph_dying};
}

void IJK_Thermal_cut_cell::lire_temperature(const IJK_Splitting& splitting, int idx)
{
  Cout << "Reading initial temperature field T" << rang_ << " from file " << fichier_reprise_temperature_ << " timestep= " << timestep_reprise_temperature_ << finl;
  const Nom& geom_name = splitting.get_grid_geometry().le_nom();
  lire_dans_lata_cut_cell(true, fichier_reprise_temperature_, timestep_reprise_temperature_, geom_name, Nom("TEMPERATURE_") + Nom(idx),
                          temperature_); // fonction qui lit un champ a partir d'un lata .
  temperature_->echange_espace_virtuel(temperature_->ghost()); // It is essential to fill the EV because the first call to convection needs them.
}

void IJK_Thermal_cut_cell::compute_interfacial_temperature2(ArrOfDouble& interfacial_temperature,
                                                            ArrOfDouble& flux_normal_interp)
{
  const int nb_facettes = ref_ijk_ft_->itfce().maillage_ft_ijk().nb_facettes();
  interfacial_temperature.resize_array(nb_facettes);
  flux_normal_interp.resize_array(nb_facettes);

  // Peut-etre que des facettes virtuelles ont ete ajoutees, mais je pense que le nombre de facettes reeles n'a pas change.
  // Mise-a-jour des tableaux pour cette eventualite :
  interfacial_temperature_.centre.resize(nb_facettes);
  interfacial_phin_ai_.resize(nb_facettes);
  ref_ijk_ft_->itfce().maillage_ft_ijk().desc_facettes().echange_espace_virtuel(interfacial_temperature_.centre);
  ref_ijk_ft_->itfce().maillage_ft_ijk().desc_facettes().echange_espace_virtuel(interfacial_phin_ai_);

  int dimension_temp = interfacial_temperature_.centre.dimension(0);
  assert(interfacial_phin_ai_.dimension(0) == dimension_temp);
  if ((dimension_temp == 0) || (deactivate_diffusion_interface_))
    {
      // Fields not already computed; write zero
      for (int fa7 = 0; fa7 < nb_facettes; fa7++)
        {
          interfacial_temperature(fa7) = 0.;
          flux_normal_interp(fa7) = 0.;
        }
    }
  else if (dimension_temp == nb_facettes)
    {
      for (int fa7 = 0; fa7 < nb_facettes; fa7++)
        {
          interfacial_temperature(fa7) = interfacial_temperature_.centre(fa7);
          flux_normal_interp(fa7) = interfacial_phin_ai_(fa7);
        }
    }
  else
    {
      Cerr << "IJK_Thermal_cut_cell::compute_interfacial_temperature2: Unexpected discrepancy in the number of facets." << finl;
      Process::exit();
    }
}

