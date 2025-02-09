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

#include <IJK_Thermal_cut_cell.h>
#include <Probleme_FTD_IJK_cut_cell.h>
#include <Cut_cell_tools.h>
#include <Probleme_FTD_IJK.h>
#include <DebogIJK.h>
#include <IJK_Navier_Stokes_tools.h>
#include <IJK_Navier_Stokes_tools_cut_cell.h>
#include <Cut_cell_tools.h>
#include <Cut_cell_convection_auxiliaire.h>
#include <Cut_cell_diffusion_auxiliaire.h>
#include <OpConvQuickIJKScalar_cut_cell.h>
#include <Cut_cell_diffusion_flux_interface.h>

Implemente_instanciable_sans_constructeur( IJK_Thermal_cut_cell, "IJK_Thermal_cut_cell", IJK_Thermal_base ) ;

IJK_Thermal_cut_cell::IJK_Thermal_cut_cell()
{
  single_phase_=0;
  conserv_energy_global_=0; // Note : doit etre zero sinon la rustine est appliquee

  cut_cell_conv_scheme_.scheme = CUT_CELL_SCHEMA_CONVECTION::QUICK_OU_CENTRE2_STENCIL;


  temperature_ = std::make_shared<Cut_field_double>();
  div_coeff_grad_T_volume_ = std::make_shared<Cut_field_double>();
  d_temperature_ = std::make_shared<Cut_field_double>();
}

Sortie& IJK_Thermal_cut_cell::printOn( Sortie& os ) const
{

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

  param.ajouter("methode_flux_interface", (int*)&methode_flux_interface_);
  param.dictionnaire("non_initialise", (int)METHODE_FLUX_INTERFACE::NON_INITIALISE);
  param.dictionnaire("interp_pure", (int)METHODE_FLUX_INTERFACE::INTERP_PURE);
  param.dictionnaire("interp_pure_no_jump", (int)METHODE_FLUX_INTERFACE::INTERP_PURE_NO_JUMP);
  param.dictionnaire("interp_cut_cell", (int)METHODE_FLUX_INTERFACE::INTERP_CUT_CELL);
  param.dictionnaire("interp_cut_cell_no_jump", (int)METHODE_FLUX_INTERFACE::INTERP_CUT_CELL_NO_JUMP);

  param.ajouter("verbosite", &verbosite_);
}

int IJK_Thermal_cut_cell::initialize(const Domaine_IJK& splitting, const int idx)
{
  Cout << que_suis_je() << "::initialize()" << finl;

  // Cut-cell variables
  ref_ijk_ft_cut_cell_ = ref_cast(Probleme_FTD_IJK_cut_cell, ref_ijk_ft_.valeur());

  Cut_field_double& cut_field_temperature                = static_cast<Cut_field_double&>(*temperature_);
  cut_field_temperature.allocate_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), splitting, Domaine_IJK::ELEM, ghost_cells_); // Overrides the allocate in IJK_Thermal_base::initialize

  if (runge_kutta_fluxes_convection_)
    {
      allocate_velocity_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), current_fluxes_conv_, splitting, 2);
      allocate_velocity_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), RK3_F_fluxes_conv_, splitting, 2);
    }

  if (runge_kutta_fluxes_diffusion_)
    {
      allocate_velocity_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), current_fluxes_diff_, splitting, 2);
      allocate_velocity_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), RK3_F_fluxes_diff_, splitting, 2);
    }

  cellule_rk_restreint_conv_.allocate_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), splitting, Domaine_IJK::ELEM, ghost_cells_);
  cellule_rk_restreint_diff_.allocate_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), splitting, Domaine_IJK::ELEM, ghost_cells_);

  if (postraiter_champs_intermediaires_)
    {
      temperature_post_dying_.allocate_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), splitting, Domaine_IJK::ELEM, 2);
      temperature_post_regular_.allocate_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), splitting, Domaine_IJK::ELEM, 2);
      temperature_post_convection_.allocate_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), splitting, Domaine_IJK::ELEM, 2);
      temperature_post_diff_regular_.allocate_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), splitting, Domaine_IJK::ELEM, 2);
    }

  cut_cell_flux_diffusion_[0].associer_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());
  cut_cell_flux_diffusion_[1].associer_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());
  cut_cell_flux_diffusion_[2].associer_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());
  cut_cell_flux_convection_[0].associer_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());
  cut_cell_flux_convection_[1].associer_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());
  cut_cell_flux_convection_[2].associer_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());


  int nalloc = IJK_Thermal_base::initialize(splitting, idx);
  lambda_.allocate(splitting, Domaine_IJK::ELEM, 1);
  nalloc += 2;

  temperature_ft_.allocate(ref_ijk_ft_cut_cell_->get_domaine_ft(), Domaine_IJK::ELEM, 4);
  nalloc += 1;

  for (int next_time = 0; next_time < 2; next_time++)
    {
      flux_interface_ft_scalar_old_.allocate(ref_ijk_ft_cut_cell_->get_domaine_ft(), Domaine_IJK::ELEM, 2);
      flux_interface_ns_scalar_old_.allocate(ref_ijk_ft_cut_cell_->get_domaine(), Domaine_IJK::ELEM, 2);
      flux_interface_ft_scalar_next_.allocate(ref_ijk_ft_cut_cell_->get_domaine_ft(), Domaine_IJK::ELEM, 2);
      flux_interface_ns_scalar_next_.allocate(ref_ijk_ft_cut_cell_->get_domaine(), Domaine_IJK::ELEM, 2);
      flux_interface_ft_scalar_old_.data() = 0.;
      flux_interface_ns_scalar_old_.data() = 0.;
      flux_interface_ft_scalar_next_.data() = 0.;
      flux_interface_ns_scalar_next_.data() = 0.;
    }
  flux_interface_efficace_scalar_.associer_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());

  convective_correction_.initialise(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());
  diffusive_correction_.initialise(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());
  diffusive_correction_.associer(flux_interface_efficace_scalar_);

  Cut_field_double& cut_field_div_coeff_grad_T_volume = static_cast<Cut_field_double&>(*div_coeff_grad_T_volume_);
  cut_field_div_coeff_grad_T_volume.allocate_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), splitting, Domaine_IJK::ELEM, 2);
  nalloc += 1;

  Cut_field_double& cut_field_d_temperature           = static_cast<Cut_field_double&>(*d_temperature_);
  cut_field_d_temperature.allocate_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc(), splitting, Domaine_IJK::ELEM, 2); // Overrides the allocate in IJK_Thermal_base::initialize


  if (temperature_diffusion_op_.non_nul())
    {
      temperature_diffusion_op_.typer("OpDiffIJKScalar_cut_cell_double");
      temperature_diffusion_op_.initialize(splitting);
      temperature_diffusion_op_->set_uniform_lambda_liquid(lambda_liquid_);
      temperature_diffusion_op_->set_uniform_lambda_vapour(lambda_vapour_);
      temperature_diffusion_op_->set_lambda(lambda_);
    }

  if (temperature_convection_op_.non_nul())
    {
      temperature_convection_op_.typer("OpConvQuickIJKScalar_cut_cell_double");
      temperature_convection_op_.initialize(splitting);
    }

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
      rho_cp_.allocate(splitting, Domaine_IJK::ELEM, 2);
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
      rho_cp_T_.allocate(splitting, Domaine_IJK::ELEM, 2);
      div_rho_cp_T_.allocate(splitting, Domaine_IJK::ELEM, 0);
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
  const double rho_l = ref_ijk_ft_cut_cell_->milieu_ijk().get_rho_liquid();
  const double rho_v = ref_ijk_ft_cut_cell_->milieu_ijk().get_rho_vapour();
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

void IJK_Thermal_cut_cell::perform_thermal_step(double total_timestep, int flag_rk, int rk_step)
{
  if (temperature_diffusion_op_.non_nul())
    ref_cast(OpDiffIJKScalar_cut_cell_double, temperature_diffusion_op_.valeur()).initialise_cut_cell(false, cut_cell_flux_diffusion_, ref_ijk_ft_cut_cell_->treatment_count(), ref_ijk_ft_cut_cell_->new_treatment());
  if (temperature_convection_op_.non_nul())
    ref_cast(OpConvQuickIJKScalar_cut_cell_double, temperature_convection_op_.valeur()).initialise_cut_cell(cut_cell_conv_scheme_, false, cut_cell_flux_convection_, ref_ijk_ft_cut_cell_->treatment_count(), ref_ijk_ft_cut_cell_->new_treatment());

  double fractional_timestep = (flag_rk) ? compute_fractionnal_timestep_rk3(total_timestep, rk_step) : total_timestep;

  if (debug_)
    Cerr << "Thermal Euler time-step" << finl;

  Cut_field_double &cut_field_temperature = static_cast<Cut_field_double&>(*temperature_);
  Cut_field_double &cut_field_d_temperature = static_cast<Cut_field_double&>(*d_temperature_);

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
      current_fluxes_conv_.set_to_uniform_value(0);
      if (rk_step == 0)
        {
          RK3_F_fluxes_conv_.set_to_uniform_value(0);
        }
      current_fluxes_conv_.echange_espace_virtuel(1);
      RK3_F_fluxes_conv_.echange_espace_virtuel(1);
    }

  if (runge_kutta_fluxes_diffusion_)
    {
      current_fluxes_diff_.set_to_uniform_value(0);
      if (rk_step == 0)
        {
          RK3_F_fluxes_diff_.set_to_uniform_value(0);
        }
      current_fluxes_diff_.echange_espace_virtuel(1);
      RK3_F_fluxes_diff_.echange_espace_virtuel(1);
    }




  update_thermal_properties();

  const Cut_field_vector3_double& cut_field_total_velocity = ref_ijk_ft_cut_cell_->get_cut_field_velocity();

  const double current_time = ref_ijk_ft_cut_cell_->schema_temps_ijk().get_current_time();

  // Toujours necessaire, ou presque
  calculer_flux_interface_old(methode_flux_interface_, lambda_liquid_, lambda_vapour_, interfacial_temperature_, interfacial_phin_ai_, cut_field_temperature, ref_ijk_ft_cut_cell_->get_cut_cell_facettes_interpolation(), flux_interface_ft_scalar_old_);

  // Calcul du flux interface sur le domaine NS :
  ref_ijk_ft_cut_cell_->eq_ns().redistrib_from_ft_elem().redistribute(flux_interface_ft_scalar_old_, flux_interface_ns_scalar_old_);
  flux_interface_ns_scalar_old_.echange_espace_virtuel(flux_interface_ns_scalar_old_.ghost());

  if (!conv_temperature_negligible_)
    {
      convective_correction_.calcule_valeur_remplissage(fractional_timestep, lambda_liquid_, lambda_vapour_, flux_interface_ns_scalar_old_, interfacial_temperature_, temperature_ft_, cut_field_total_velocity, cut_field_temperature);
    }


  ene_ini = cut_field_temperature.compute_global_energy_cut_cell(0, get_rhocp_l(), get_rhocp_v());

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
      set_rk_restreint(rk_step, runge_kutta_restriction_leniency_convection_, cut_cell_disc, cellule_rk_restreint_conv_);
      ref_cast(OpConvQuickIJKScalar_cut_cell_double, temperature_convection_op_.valeur()).set_runge_kutta(rk_step, total_timestep, current_fluxes_conv_, RK3_F_fluxes_conv_, cellule_rk_restreint_conv_);
    }

  compute_temperature_convection_cut_cell(cut_field_total_velocity, cut_field_d_temperature);
  d_ene_Conv = cut_field_d_temperature.compute_d_global_energy_cut_cell(0);

  euler_explicit_update_cut_cell_transport(fractional_timestep, cut_field_d_temperature, cut_field_temperature);
  cut_field_temperature.echange_espace_virtuel(temperature_->ghost());

  if (postraiter_champs_intermediaires_)
    {
      temperature_post_regular_.copy_from(cut_field_temperature);
    }

  ene_postconv_switch = cut_field_temperature.compute_global_energy_cut_cell(1, get_rhocp_l(), get_rhocp_v());

  if (verbosite_ >= 7)
    print_Tmin_Tmax_cut_cell(cut_field_temperature, 1, current_time, "Temperature Min/Max apres la convection principale");

  if (!conv_temperature_negligible_)
    {
      if (runge_kutta_fluxes_convection_ && (!runge_kutta_fluxes_pas_de_correction_convection_))
        {
          current_fluxes_conv_.set_to_uniform_value(0);
        }

      convective_correction_.add_dying_cells(cut_field_total_velocity, cut_field_temperature, runge_kutta_fluxes_convection_, current_fluxes_conv_);
      cut_field_temperature.echange_espace_virtuel(temperature_->ghost());
    }

  if (postraiter_champs_intermediaires_)
    {
      temperature_post_dying_.copy_from(cut_field_temperature);
    }
  ene_postconv_dying = cut_field_temperature.compute_global_energy_cut_cell(1, get_rhocp_l(), get_rhocp_v());

  if (verbosite_ >= 7)
    print_Tmin_Tmax_cut_cell(cut_field_temperature, 1, current_time, "Temperature Min/Max apres la convection des cellules mourrantes");

  if (!conv_temperature_negligible_)
    {
      convective_correction_.add_small_nascent_cells(cut_field_total_velocity, cut_field_temperature, runge_kutta_fluxes_convection_, current_fluxes_conv_);
      cut_field_temperature.echange_diph_vers_pure_cellules_finalement_pures();
      cut_field_temperature.remplir_tableau_pure_cellules_diphasiques(true /* next time */);
      cut_field_temperature.echange_espace_virtuel(temperature_->ghost());
    }

  if (postraiter_champs_intermediaires_)
    {
      temperature_post_convection_.copy_from(cut_field_temperature);
    }

  if (!diff_temperature_negligible_ && (!diffusive_correction_.deactivate_correction_petites_cellules_diffusion_))
    {
      diffusive_correction_.calcule_valeur_remplissage(fractional_timestep, lambda_liquid_, lambda_vapour_, flux_interface_ns_scalar_old_, interfacial_temperature_, temperature_ft_, cut_field_total_velocity, cut_field_temperature);
    }
  ene_postconv_small_nascent = cut_field_temperature.compute_global_energy_cut_cell(1, get_rhocp_l(), get_rhocp_v());

  if (verbosite_ >= 7)
    print_Tmin_Tmax_cut_cell(cut_field_temperature, 1, current_time, "Temperature Min/Max apres la convection des cellules naissantes");

  if (runge_kutta_fluxes_convection_ && (!runge_kutta_fluxes_pas_de_correction_convection_))
    {
      add_flux_times_vol_over_dt_surface(fractional_timestep, current_fluxes_conv_, RK3_F_fluxes_conv_);
    }


  if (!diff_temperature_negligible_ && (!diffusive_correction_.deactivate_correction_petites_cellules_diffusion_))
    {
      diffusive_correction_.compute_flux_small_nascent_cells(cut_field_total_velocity, cut_field_temperature);
    }

  if (!deactivate_diffusion_interface_)
    {
      calculer_flux_interface_next(methode_flux_interface_, lambda_liquid_, lambda_vapour_, interfacial_temperature_, interfacial_phin_ai_, cut_field_temperature, ref_ijk_ft_cut_cell_->get_cut_cell_facettes_interpolation(), flux_interface_ft_scalar_next_);

      // Calcul du flux interface sur le domaine NS :
      ref_ijk_ft_->eq_ns().redistrib_from_ft_elem().redistribute(flux_interface_ft_scalar_next_, flux_interface_ns_scalar_next_);
      flux_interface_ns_scalar_next_.echange_espace_virtuel(flux_interface_ns_scalar_next_.ghost());

      calculer_flux_interface_efficace(flux_interface_ns_scalar_old_, flux_interface_ns_scalar_next_, flux_interface_efficace_scalar_);
    }
  else
    {
      //Cerr << "Warning: The diffusion at the interface is not activated." << finl;
    }

  if (runge_kutta_fluxes_diffusion_)
    {
      set_rk_restreint(rk_step, runge_kutta_restriction_leniency_diffusion_, cut_cell_disc, cellule_rk_restreint_diff_);
      ref_cast(OpDiffIJKScalar_cut_cell_double, temperature_diffusion_op_.valeur()).set_runge_kutta(rk_step, total_timestep, current_fluxes_diff_, RK3_F_fluxes_diff_, cellule_rk_restreint_diff_);
    }

  add_temperature_diffusion();
  d_ene_Diffu = cut_field_d_temperature.compute_d_global_energy_cut_cell(1);

  euler_explicit_update_cut_cell_notransport(fractional_timestep, true, cut_field_d_temperature, cut_field_temperature);

  cut_field_temperature.echange_espace_virtuel(temperature_->ghost());

  if (postraiter_champs_intermediaires_)
    {
      temperature_post_diff_regular_.copy_from(cut_field_temperature);
    }
  ene_postdiff_regular = cut_field_temperature.compute_global_energy_cut_cell(1, get_rhocp_l(), get_rhocp_v());

  if (verbosite_ >= 7)
    print_Tmin_Tmax_cut_cell(cut_field_temperature, 1, current_time, "Temperature Min/Max apres la diffusion principale");

  if (!diff_temperature_negligible_)
    {
      if (runge_kutta_fluxes_diffusion_ && (!runge_kutta_fluxes_pas_de_correction_diffusion_))
        {
          current_fluxes_diff_.set_to_uniform_value(0);
        }

      diffusive_correction_.add_dying_cells(cut_field_total_velocity, cut_field_temperature, runge_kutta_fluxes_diffusion_, current_fluxes_diff_);
      cut_field_temperature.echange_espace_virtuel(temperature_->ghost());
    }

  ene_postdiff_dying = cut_field_temperature.compute_global_energy_cut_cell(1, get_rhocp_l(), get_rhocp_v());

  if (verbosite_ >= 7)
    print_Tmin_Tmax_cut_cell(cut_field_temperature, 1, current_time, "Temperature Min/Max apres la diffusion des cellules mourrantes");

  if (!diff_temperature_negligible_ && (!diffusive_correction_.deactivate_correction_petites_cellules_diffusion_))
    {
      diffusive_correction_.add_small_nascent_cells(cut_field_total_velocity, cut_field_temperature, runge_kutta_fluxes_diffusion_, current_fluxes_diff_);
      cut_field_temperature.echange_diph_vers_pure_cellules_finalement_pures();
      cut_field_temperature.remplir_tableau_pure_cellules_diphasiques(true /* next time */);
      cut_field_temperature.echange_espace_virtuel(temperature_->ghost());
    }

  ene_postdiff_small = cut_field_temperature.compute_global_energy_cut_cell(1, get_rhocp_l(), get_rhocp_v());

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
      add_flux_times_vol_over_dt_surface(fractional_timestep, current_fluxes_diff_, RK3_F_fluxes_diff_);
    }



}

void IJK_Thermal_cut_cell::euler_time_step(const double timestep)
{
  perform_thermal_step(timestep, 0, 0);
}

void IJK_Thermal_cut_cell::rk3_sub_step(const int rk_step, const double total_timestep,
                                        const double time)
{
  Cut_field_double& cut_field_temperature = static_cast<Cut_field_double&>(*temperature_);

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
    latastep_reprise_ = latastep_reprise_ini_ + ref_ijk_ft_->schema_temps_ijk().get_tstep() + 1;
}


void cut_cell_reinit_streamObj(std::ostringstream& streamObj, const double& param)
{
  streamObj.str("");
  streamObj.clear();
  streamObj << (double) param;
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
      temperature_convection_op_->calculer(cut_field_temperature,
                                           cut_field_total_velocity[0],
                                           cut_field_total_velocity[1],
                                           cut_field_total_velocity[2],
                                           cut_field_d_temperature);
      cut_field_d_temperature.divide_by_scalar(vol_, vol_);
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

  temperature_diffusion_op_->set_lambda(lambda_);

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
      temperature_diffusion_op_->calculer(cut_field_temperature,
                                          cut_field_d_temperature, // On utilise directement cut_field_d_temperature, sans passer par div_coeff_grad_T_volume_
                                          boundary_flux_kmin_,
                                          boundary_flux_kmax_);

      ajout_flux_interface_a_divergence(flux_interface_efficace_scalar_, cut_field_d_temperature);

      // On copie cut_field_d_temperature dans cut_field_div_coeff_grad_T_volume.
      // Cela ne sert a rien a part conserver le comportement historique de ces champs.
      cut_field_div_coeff_grad_T_volume.copy_from(cut_field_d_temperature);

      const double rho_l = ref_ijk_ft_cut_cell_->milieu_ijk().get_rho_liquid();
      const double rho_v = ref_ijk_ft_cut_cell_->milieu_ijk().get_rho_vapour();
      cut_field_d_temperature.divide_by_scalar(rho_l*cp_liquid_*vol_, rho_v*cp_vapour_*vol_);

      statistiques().end_count(cnt_diff_temp);
      DebogIJK::verifier("div_coeff_grad_T_volume_", *div_coeff_grad_T_volume_);
    }
}

void IJK_Thermal_cut_cell::print_Tmin_Tmax_cut_cell(const Cut_field_double& cut_field_temperature, bool next, double current_time, const std::string& heading)
{
  CutCell_GlobalInfo Tmin = cut_field_temperature.compute_min_cut_cell(next);
  CutCell_GlobalInfo Tmax = cut_field_temperature.compute_max_cut_cell(next);
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
      Cerr << "Location T_MinMax, overall_l time: " << current_time << " Tmin: " << cut_field_temperature.get_value_location(Tmin.overall_l) << " Tmax: " << cut_field_temperature.get_value_location(Tmax.overall_l) << finl;
      Cerr << "Location T_MinMax, overall_v time: " << current_time << " Tmin: " << cut_field_temperature.get_value_location(Tmin.overall_v) << " Tmax: " << cut_field_temperature.get_value_location(Tmax.overall_v) << finl;
      Cerr << "Location T_MinMax, pure.......... time: " << current_time << " Tmin: " << cut_field_temperature.get_value_location(Tmin.pure)         << " Tmax: " << cut_field_temperature.get_value_location(Tmax.pure)         << finl;
      Cerr << "Location T_MinMax, diph_l........ time: " << current_time << " Tmin: " << cut_field_temperature.get_value_location(Tmin.diph_l)       << " Tmax: " << cut_field_temperature.get_value_location(Tmax.diph_l)       << finl;
      Cerr << "Location T_MinMax, diph_v........ time: " << current_time << " Tmin: " << cut_field_temperature.get_value_location(Tmin.diph_v)       << " Tmax: " << cut_field_temperature.get_value_location(Tmax.diph_v)       << finl;
      Cerr << "Location T_MinMax, diph_small.... time: " << current_time << " Tmin: " << cut_field_temperature.get_value_location(Tmin.diph_small)   << " Tmax: " << cut_field_temperature.get_value_location(Tmax.diph_small)   << finl;
      Cerr << "Location T_MinMax, diph_regular.. time: " << current_time << " Tmin: " << cut_field_temperature.get_value_location(Tmin.diph_regular) << " Tmax: " << cut_field_temperature.get_value_location(Tmax.diph_regular) << finl;
      Cerr << "Location T_MinMax, diph_nascent.. time: " << current_time << " Tmin: " << cut_field_temperature.get_value_location(Tmin.diph_nascent) << " Tmax: " << cut_field_temperature.get_value_location(Tmax.diph_nascent) << finl;
      Cerr << "Location T_MinMax, diph_dying.... time: " << current_time << " Tmin: " << cut_field_temperature.get_value_location(Tmin.diph_dying)   << " Tmax: " << cut_field_temperature.get_value_location(Tmax.diph_dying)   << finl;
    }
}

void IJK_Thermal_cut_cell::lire_temperature(const Domaine_IJK& geom, int idx)
{
  Cout << "Reading initial temperature field T" << rang_ << " from file " << fichier_reprise_temperature_ << " timestep= " << timestep_reprise_temperature_ << finl;
  const Nom& geom_name = geom.le_nom();
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
  interfacial_temperature_.resize(nb_facettes);
  interfacial_phin_ai_.resize(nb_facettes);
  ref_ijk_ft_->itfce().maillage_ft_ijk().desc_facettes().echange_espace_virtuel(interfacial_temperature_);
  ref_ijk_ft_->itfce().maillage_ft_ijk().desc_facettes().echange_espace_virtuel(interfacial_phin_ai_);

  int dimension_temp = interfacial_temperature_.dimension(0);
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
          interfacial_temperature(fa7) = interfacial_temperature_(fa7);
          flux_normal_interp(fa7) = interfacial_phin_ai_(fa7);
        }
    }
  else
    {
      Cerr << "IJK_Thermal_cut_cell::compute_interfacial_temperature2: Unexpected discrepancy in the number of facets." << finl;
      Process::exit();
    }
}

