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
// File      : IJK_Thermal_cut_cell.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_Thermal_cut_cell_included
#define IJK_Thermal_cut_cell_included

#include <IJK_Thermal_base.h>
#include <IJK_Field.h>
#include <Boundary_Conditions_Thermique.h>
#include <IJK_Splitting.h>
#include <IJK_Field.h>
#include <Parser.h>
#include <IJK_Lata_writer.h>
#include <OpConvQuickIJKScalar.h>
#include <OpConvCentre2IJKScalar.h>
#include <Ouvrir_fichier.h>
#include <Corrige_flux_FT.h>
#include <TRUST_Ref.h>
#include <Operateur_IJK_elem_diff_base.h>
#include <OpConvAmontIJK.h>
#include <OpConvDiscQuickIJKScalar.h>
#include <OpConvCentre4IJK.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_Thermal_cut_cell
//
// <Description of class IJK_Thermal_cut_cell>
//
/////////////////////////////////////////////////////////////////////////////

struct CutCell_GlobalInfo;

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

  void calculer_flux_interface();

  void euler_time_step(const double timestep) override;
  void rk3_sub_step(const int rk_step, const double total_timestep, const double time) override;

  CutCell_GlobalInfo compute_global_energy_cut_cell(Cut_field_scalar& cut_field_temperature, bool next);
  double compute_global_energy() override
  {
    Cerr << "I want to make sure the function compute_global_energy is not used." << finl;
    Process::exit();
    return 0.;
  }
  CutCell_GlobalInfo compute_Tmin_cut_cell(Cut_field_scalar& cut_field_temperature, bool next);
  CutCell_GlobalInfo compute_Tmax_cut_cell(Cut_field_scalar& cut_field_temperature, bool next);

  void remplir_cellules_diphasiques() override
  {
    cut_field_temperature_.remplir_cellules_diphasiques();
  }
  void remplir_cellules_devenant_diphasiques() override
  {
    cut_field_temperature_.remplir_cellules_devenant_diphasiques();
  }
  void remplir_cellules_maintenant_pures() override
  {
    cut_field_temperature_.remplir_cellules_maintenant_pures();
  }
  void transfert_diphasique_vers_pures() override
  {
    cut_field_temperature_.transfert_diphasique_vers_pures();
  }

protected :
  friend class IJK_FT_Post;

  void correct_any_temperature_fields_for_eulerian_fluxes(IJK_Field_double& temperature);
  void compare_temperature_fields(const IJK_Field_double& temperature,
                                  const IJK_Field_double& temperature_ana,
                                  IJK_Field_double& error_temperature_ana,
                                  IJK_Field_double& error_temperature_ana_rel);
  void evaluate_total_liquid_absolute_parameter(const IJK_Field_double& field,
                                                double& total_parameter);
  void evaluate_total_liquid_parameter_squared(const IJK_Field_double& field,
                                               double& total_parameter);
  void correct_any_temperature_field_for_visu(IJK_Field_double& temperature);

  Nom compute_quasi_static_spherical_diffusion_expression(const double& time_scope, const int index_bubble, const int index_bubble_real);
  Nom generate_expression_temperature_ini(const double& time_scope, const double x, const double y, const double z);
  void calculer_ecart_T_ana() override { ; };

  int computed_centred_bubble_start_;
  double single_centred_bubble_radius_ini_;
  int spherical_diffusion_;
  int allow_temperature_correction_for_visu_;
  int override_vapour_mixed_values_; // For debug purposes
  int source_terms_type_;
  Motcles source_terms_type_dict_;
  int source_terms_correction_;
  double delta_T_subcooled_overheated_=-1.;

  double error_temperature_ana_total_;
  double error_temperature_ana_squared_total_;
  double error_temperature_ana_rel_total_;

  void calculer_dT_cut_cell(const Cut_field_vector& cut_field_velocity);
  void compute_temperature_convection_cut_cell(const Cut_field_vector& cut_field_velocity);
  void add_temperature_diffusion() override;
  void compute_diffusion_increment() override;
  /* correct_temperature_for_eulerian_fluxes() May be clearly overridden later */
  void correct_temperature_for_eulerian_fluxes() override { ; };
  double get_rho_cp_ijk(int i, int j, int k) const;
  double get_rho_cp_u_ijk(const IJK_Field_double& vx, int i, int j, int k) const override;
  double compute_rho_cp_u_mean(const IJK_Field_double& vx) override;
  double get_div_lambda_ijk(int i, int j, int k) const override;
  double compute_temperature_dimensionless_theta_mean(const IJK_Field_double& vx) override;

  void add_convection_dying_cells(const Cut_field_vector& cut_field_velocity);
  void add_convection_nascent_cells(const Cut_field_vector& cut_field_velocity);

  //Rustine
  double E0_;//volumique
  IJK_Field_double T_rust_;
  void compute_T_rust(const FixedVector<IJK_Field_double, 3>& velocity);

  int type_temperature_convection_form_;

  REF(IJK_FT_cut_cell) ref_ijk_ft_cut_cell_;

  IJK_Field_double div_rho_cp_T_;
  IJK_Field_double lambda_;

  IJK_Field_double RK3_F_temperature_diff_;

  Cut_field_scalar cut_field_temperature_;
  Cut_field_scalar cut_field_RK3_F_temperature_diff_;
  Cut_field_scalar cut_field_RK3_F_temperature_conv_;
  Cut_cell_vector cut_cell_flux_diffusion_;
  Cut_cell_vector cut_cell_flux_convection_;
  Cut_field_scalar cut_field_div_coeff_grad_T_volume_;
  Cut_field_scalar cut_field_d_temperature_;

  IJK_Field_double flux_interface_ns_;
  IJK_Field_double flux_interface_ft_;
  DoubleTabFT_cut_cell_scalar flux_interface_;

  // Cut cell options
  int activate_diffusion_interface_;
};

#endif /* IJK_Thermal_cut_cell_included */
