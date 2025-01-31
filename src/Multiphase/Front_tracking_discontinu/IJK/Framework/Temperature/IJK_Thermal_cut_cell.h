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

#ifndef IJK_Thermal_cut_cell_included
#define IJK_Thermal_cut_cell_included

#include <IJK_Thermal_base.h>
#include <IJK_Field_vector.h>
#include <IJK_Field.h>
#include <Boundary_Conditions_Thermique.h>
#include <IJK_Splitting.h>
#include <IJK_Field.h>
#include <Parser.h>
#include <IJK_Lata_writer.h>
#include <OpConvQuickIJKScalar.h>
#include <OpConvCentre2IJKScalar.h>
#include <Ouvrir_fichier.h>
#include <TRUST_Ref.h>
#include <Operateur_IJK_elem_diff_base.h>
#include <OpConvAmontIJK.h>
#include <OpConvDiscQuickIJKScalar.h>
#include <OpConvCentre4IJK.h>
#include <Cut_cell_correction_petites_cellules.h>
#include <Cut_cell_convection_auxiliaire.h>
#include <Cut_cell_diffusion_auxiliaire.h>


class IJK_Thermal_cut_cell : public IJK_Thermal_base
{

  Declare_instanciable( IJK_Thermal_cut_cell ) ;

public :

  int initialize(const IJK_Splitting& splitting, const int idx) override;
  double compute_timestep(const double timestep,
                          const double rho_l, const double rho_v,
                          const double dxmin);
  void update_thermal_properties() override;
  void set_param( Param& param ) override;
  void compute_temperature_init() override;
  void recompute_temperature_init() override;

  void perform_thermal_step(double total_timestep, int flag_rk, int rk_step);
  void euler_time_step(const double timestep) override;
  void rk3_sub_step(const int rk_step, const double total_timestep, const double time) override;
  void sauvegarder_temperature(Nom& lata_name, int idx, const int& stop=0) override;

  double compute_global_energy() override
  {
    Cerr << "I want to make sure the function compute_global_energy is not used." << finl;
    Process::exit();
    return 0.;
  }
  void print_Tmin_Tmax_cut_cell(const Cut_field_double& cut_field_temperature, bool next, double current_time, const std::string& heading);
  void calculer_flux_interface();

  void copie_pure_vers_diph_sans_interpolation() override
  {
    Cut_field_double& cut_field_temperature = static_cast<Cut_field_double&>(*temperature_);
    cut_field_temperature.copie_pure_vers_diph_sans_interpolation();

    if (runge_kutta_fluxes_convection_)
      {
        RK3_F_fluxes_conv_[0].copie_pure_vers_diph_sans_interpolation();
        RK3_F_fluxes_conv_[1].copie_pure_vers_diph_sans_interpolation();
        RK3_F_fluxes_conv_[2].copie_pure_vers_diph_sans_interpolation();
      }

    if (runge_kutta_fluxes_diffusion_)
      {
        RK3_F_fluxes_diff_[0].copie_pure_vers_diph_sans_interpolation();
        RK3_F_fluxes_diff_[1].copie_pure_vers_diph_sans_interpolation();
        RK3_F_fluxes_diff_[2].copie_pure_vers_diph_sans_interpolation();
      }

    cellule_rk_restreint_conv_.copie_pure_vers_diph_sans_interpolation();
    cellule_rk_restreint_diff_.copie_pure_vers_diph_sans_interpolation();
  }
  void echange_pure_vers_diph_cellules_initialement_pures() override
  {
    Cut_field_double& cut_field_temperature = static_cast<Cut_field_double&>(*temperature_);
    cut_field_temperature.echange_pure_vers_diph_cellules_initialement_pures();

    if (runge_kutta_fluxes_convection_)
      {
        RK3_F_fluxes_conv_[0].echange_pure_vers_diph_cellules_initialement_pures();
        RK3_F_fluxes_conv_[1].echange_pure_vers_diph_cellules_initialement_pures();
        RK3_F_fluxes_conv_[2].echange_pure_vers_diph_cellules_initialement_pures();
      }

    if (runge_kutta_fluxes_diffusion_)
      {
        RK3_F_fluxes_diff_[0].echange_pure_vers_diph_cellules_initialement_pures();
        RK3_F_fluxes_diff_[1].echange_pure_vers_diph_cellules_initialement_pures();
        RK3_F_fluxes_diff_[2].echange_pure_vers_diph_cellules_initialement_pures();
      }

    cellule_rk_restreint_conv_.echange_pure_vers_diph_cellules_initialement_pures();
    cellule_rk_restreint_diff_.echange_pure_vers_diph_cellules_initialement_pures();
  }
  void echange_diph_vers_pure_cellules_finalement_pures() override
  {
    Cut_field_double& cut_field_temperature = static_cast<Cut_field_double&>(*temperature_);
    cut_field_temperature.echange_diph_vers_pure_cellules_finalement_pures();

    if (runge_kutta_fluxes_convection_)
      {
        RK3_F_fluxes_conv_[0].echange_diph_vers_pure_cellules_finalement_pures();
        RK3_F_fluxes_conv_[1].echange_diph_vers_pure_cellules_finalement_pures();
        RK3_F_fluxes_conv_[2].echange_diph_vers_pure_cellules_finalement_pures();
      }

    if (runge_kutta_fluxes_diffusion_)
      {
        RK3_F_fluxes_diff_[0].echange_diph_vers_pure_cellules_finalement_pures();
        RK3_F_fluxes_diff_[1].echange_diph_vers_pure_cellules_finalement_pures();
        RK3_F_fluxes_diff_[2].echange_diph_vers_pure_cellules_finalement_pures();
      }

    cellule_rk_restreint_conv_.echange_diph_vers_pure_cellules_finalement_pures();
    cellule_rk_restreint_diff_.echange_diph_vers_pure_cellules_finalement_pures();
  }
  void vide_phase_invalide_cellules_diphasiques() override
  {
    Cut_field_double& cut_field_temperature = static_cast<Cut_field_double&>(*temperature_);
    cut_field_temperature.vide_phase_invalide_cellules_diphasiques();

    if (runge_kutta_fluxes_convection_)
      {
        RK3_F_fluxes_conv_[0].vide_phase_invalide_cellules_diphasiques();
        RK3_F_fluxes_conv_[1].vide_phase_invalide_cellules_diphasiques();
        RK3_F_fluxes_conv_[2].vide_phase_invalide_cellules_diphasiques();
      }

    if (runge_kutta_fluxes_diffusion_)
      {
        RK3_F_fluxes_diff_[0].vide_phase_invalide_cellules_diphasiques();
        RK3_F_fluxes_diff_[1].vide_phase_invalide_cellules_diphasiques();
        RK3_F_fluxes_diff_[2].vide_phase_invalide_cellules_diphasiques();
      }

    cellule_rk_restreint_conv_.vide_phase_invalide_cellules_diphasiques();
    cellule_rk_restreint_diff_.vide_phase_invalide_cellules_diphasiques();
  }
  void remplir_tableau_pure_cellules_diphasiques(bool next_time) override
  {
    Cut_field_double& cut_field_temperature = static_cast<Cut_field_double&>(*temperature_);
    cut_field_temperature.remplir_tableau_pure_cellules_diphasiques(next_time);

    if (runge_kutta_fluxes_convection_)
      {
        RK3_F_fluxes_conv_[0].remplir_tableau_pure_cellules_diphasiques(next_time);
        RK3_F_fluxes_conv_[1].remplir_tableau_pure_cellules_diphasiques(next_time);
        RK3_F_fluxes_conv_[2].remplir_tableau_pure_cellules_diphasiques(next_time);
      }

    if (runge_kutta_fluxes_diffusion_)
      {
        RK3_F_fluxes_diff_[0].remplir_tableau_pure_cellules_diphasiques(next_time);
        RK3_F_fluxes_diff_[1].remplir_tableau_pure_cellules_diphasiques(next_time);
        RK3_F_fluxes_diff_[2].remplir_tableau_pure_cellules_diphasiques(next_time);
      }

    cellule_rk_restreint_conv_.remplir_tableau_pure_cellules_diphasiques_max(next_time);
    cellule_rk_restreint_diff_.remplir_tableau_pure_cellules_diphasiques_max(next_time);
  }

protected :
  friend class IJK_FT_Post;

  void lire_temperature(const IJK_Splitting& splitting, int idx) override;

  void compute_interfacial_temperature2(ArrOfDouble& interfacial_temperature, ArrOfDouble& flux_normal_interp) override;

  void compute_temperature_convection_cut_cell(const Cut_field_vector3_double& cut_field_total_velocity, Cut_field_double& cut_field_d_temperature);
  void add_temperature_diffusion() override;
  void compute_diffusion_increment() override { ; }; // Unused function in the cut-cell case
  void correct_temperature_for_eulerian_fluxes() override { ; };


  int type_temperature_convection_form_ = 0;

  OBS_PTR(Probleme_FTD_IJK_cut_cell) ref_ijk_ft_cut_cell_;

  IJK_Field_double div_rho_cp_T_;
  IJK_Field_double lambda_;

  FixedVector<Cut_cell_double, 3> cut_cell_flux_diffusion_;
  FixedVector<Cut_cell_double, 3> cut_cell_flux_convection_;

  Cut_cell_conv_scheme cut_cell_conv_scheme_;

  METHODE_FLUX_INTERFACE methode_flux_interface_ = METHODE_FLUX_INTERFACE::INTERP_CUT_CELL;
  IJK_Field_double flux_interface_ns_scalar_old_;
  IJK_Field_double flux_interface_ft_scalar_old_;
  IJK_Field_double flux_interface_ns_scalar_next_;
  IJK_Field_double flux_interface_ft_scalar_next_;
  DoubleTabFT_cut_cell_scalar flux_interface_efficace_scalar_;

  DoubleTabFT interfacial_temperature_;
  DoubleTabFT interfacial_phin_ai_;

  Cut_field_vector3_double current_fluxes_conv_;
  Cut_field_vector3_double current_fluxes_diff_;
  Cut_field_vector3_double RK3_F_fluxes_conv_;
  Cut_field_vector3_double RK3_F_fluxes_diff_;
  int runge_kutta_fluxes_convection_ = 0;
  int runge_kutta_fluxes_diffusion_ = 0;
  int runge_kutta_fluxes_pas_de_correction_convection_ = 0;
  int runge_kutta_fluxes_pas_de_correction_diffusion_ = 0;
  Cut_field_int cellule_rk_restreint_conv_;
  Cut_field_int cellule_rk_restreint_diff_;

  // Temporary fields, to inspect each step of the time advance
  int postraiter_champs_intermediaires_ = 0;
  Cut_field_double temperature_post_dying_;
  Cut_field_double temperature_post_regular_;
  Cut_field_double temperature_post_convection_;
  Cut_field_double temperature_post_diff_regular_;

  Cut_cell_convection_auxiliaire convective_correction_;
  Cut_cell_diffusion_auxiliaire diffusive_correction_;

  int deactivate_diffusion_interface_ = 0;

  int runge_kutta_restriction_leniency_convection_ = 0;
  int runge_kutta_restriction_leniency_diffusion_ = 0;


  int verbosite_ = 2;
};

#endif /* IJK_Thermal_cut_cell_included */
