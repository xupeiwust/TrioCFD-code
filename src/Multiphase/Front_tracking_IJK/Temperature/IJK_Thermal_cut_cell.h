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

  void add_flux_times_vol_over_dt_surface(double fractional_timestep, const Cut_field_vector3_double& cut_field_current_fluxes, Cut_field_vector3_double& cut_field_RK3_F_fluxes);
  void set_rk_restreint(int rk_step, int rk_restriction_leniency, const Cut_cell_FT_Disc& cut_cell_disc, IJK_Field_int& cellule_rk_restreint_v, IJK_Field_int& cellule_rk_restreint_l);
  void perform_thermal_step(double total_timestep, int flag_rk, int rk_step);
  void euler_time_step(const double timestep) override;
  void rk3_sub_step(const int rk_step, const double total_timestep, const double time) override;
  void sauvegarder_temperature(Nom& lata_name, int idx, const int& stop=0) override;

  CutCell_GlobalInfo compute_global_energy_cut_cell(Cut_field_double& cut_field_temperature, bool next);
  CutCell_GlobalInfo compute_d_global_energy_cut_cell(Cut_field_double& cut_field_d_temperature, bool next);
  double compute_global_energy() override
  {
    Cerr << "I want to make sure the function compute_global_energy is not used." << finl;
    Process::exit();
    return 0.;
  }
  void print_Tmin_Tmax_cut_cell(const Cut_field_double& cut_field_temperature, bool next, double current_time, const std::string& heading);
  CutCell_GlobalInfo compute_Tmin_cut_cell(const Cut_field_double& cut_field_temperature, bool next);
  CutCell_GlobalInfo compute_Tmax_cut_cell(const Cut_field_double& cut_field_temperature, bool next);
  void calculer_flux_interface();

  void remplir_cellules_diphasiques() override
  {
    Cut_field_double& cut_field_temperature = static_cast<Cut_field_double&>(*temperature_);
    cut_field_temperature.remplir_cellules_diphasiques();

    if (runge_kutta_fluxes_convection_)
      {
        Cut_field_vector3_double& cut_field_RK3_F_fluxes_conv = static_cast<Cut_field_vector3_double&>(RK3_F_fluxes_conv_);
        cut_field_RK3_F_fluxes_conv[0].remplir_cellules_diphasiques();
        cut_field_RK3_F_fluxes_conv[1].remplir_cellules_diphasiques();
        cut_field_RK3_F_fluxes_conv[2].remplir_cellules_diphasiques();
      }

    if (runge_kutta_fluxes_diffusion_)
      {
        Cut_field_vector3_double& cut_field_RK3_F_fluxes_diff = static_cast<Cut_field_vector3_double&>(RK3_F_fluxes_diff_);
        cut_field_RK3_F_fluxes_diff[0].remplir_cellules_diphasiques();
        cut_field_RK3_F_fluxes_diff[1].remplir_cellules_diphasiques();
        cut_field_RK3_F_fluxes_diff[2].remplir_cellules_diphasiques();
      }
  }
  void remplir_cellules_devenant_diphasiques() override
  {
    Cut_field_double& cut_field_temperature = static_cast<Cut_field_double&>(*temperature_);
    cut_field_temperature.remplir_cellules_devenant_diphasiques();

    if (runge_kutta_fluxes_convection_)
      {
        Cut_field_vector3_double& cut_field_RK3_F_fluxes_conv = static_cast<Cut_field_vector3_double&>(RK3_F_fluxes_conv_);
        cut_field_RK3_F_fluxes_conv[0].remplir_cellules_devenant_diphasiques();
        cut_field_RK3_F_fluxes_conv[1].remplir_cellules_devenant_diphasiques();
        cut_field_RK3_F_fluxes_conv[2].remplir_cellules_devenant_diphasiques();
      }

    if (runge_kutta_fluxes_diffusion_)
      {
        Cut_field_vector3_double& cut_field_RK3_F_fluxes_diff = static_cast<Cut_field_vector3_double&>(RK3_F_fluxes_diff_);
        cut_field_RK3_F_fluxes_diff[0].remplir_cellules_devenant_diphasiques();
        cut_field_RK3_F_fluxes_diff[1].remplir_cellules_devenant_diphasiques();
        cut_field_RK3_F_fluxes_diff[2].remplir_cellules_devenant_diphasiques();
      }
  }
  void remplir_cellules_maintenant_pures() override
  {
    Cut_field_double& cut_field_temperature = static_cast<Cut_field_double&>(*temperature_);
    cut_field_temperature.remplir_cellules_maintenant_pures();

    if (runge_kutta_fluxes_convection_)
      {
        Cut_field_vector3_double& cut_field_RK3_F_fluxes_conv = static_cast<Cut_field_vector3_double&>(RK3_F_fluxes_conv_);
        cut_field_RK3_F_fluxes_conv[0].remplir_cellules_maintenant_pures();
        cut_field_RK3_F_fluxes_conv[1].remplir_cellules_maintenant_pures();
        cut_field_RK3_F_fluxes_conv[2].remplir_cellules_maintenant_pures();
      }

    if (runge_kutta_fluxes_diffusion_)
      {
        Cut_field_vector3_double& cut_field_RK3_F_fluxes_diff = static_cast<Cut_field_vector3_double&>(RK3_F_fluxes_diff_);
        cut_field_RK3_F_fluxes_diff[0].remplir_cellules_maintenant_pures();
        cut_field_RK3_F_fluxes_diff[1].remplir_cellules_maintenant_pures();
        cut_field_RK3_F_fluxes_diff[2].remplir_cellules_maintenant_pures();
      }
  }
  void transfert_diphasique_vers_pures() override
  {
    Cut_field_double& cut_field_temperature = static_cast<Cut_field_double&>(*temperature_);
    cut_field_temperature.transfert_diphasique_vers_pures();

    if (runge_kutta_fluxes_convection_)
      {
        Cut_field_vector3_double& cut_field_RK3_F_fluxes_conv = static_cast<Cut_field_vector3_double&>(RK3_F_fluxes_conv_);
        cut_field_RK3_F_fluxes_conv[0].transfert_diphasique_vers_pures();
        cut_field_RK3_F_fluxes_conv[1].transfert_diphasique_vers_pures();
        cut_field_RK3_F_fluxes_conv[2].transfert_diphasique_vers_pures();
      }

    if (runge_kutta_fluxes_diffusion_)
      {
        Cut_field_vector3_double& cut_field_RK3_F_fluxes_diff = static_cast<Cut_field_vector3_double&>(RK3_F_fluxes_diff_);
        cut_field_RK3_F_fluxes_diff[0].transfert_diphasique_vers_pures();
        cut_field_RK3_F_fluxes_diff[1].transfert_diphasique_vers_pures();
        cut_field_RK3_F_fluxes_diff[2].transfert_diphasique_vers_pures();
      }
  }

protected :
  friend class IJK_FT_Post;

  void lire_temperature(const IJK_Splitting& splitting, int idx) override;

  void compute_interfacial_temperature2(ArrOfDouble& interfacial_temperature, ArrOfDouble& flux_normal_interp) override;

  void calculer_dT_cut_cell(const Cut_field_vector3_double& cut_field_total_velocity);
  void compute_temperature_convection_cut_cell(const Cut_field_vector3_double& cut_field_total_velocity, Cut_field_double& cut_field_d_temperature);
  void add_temperature_diffusion() override;
  void compute_diffusion_increment() override;
  void correct_temperature_for_eulerian_fluxes() override { ; };
  double get_temperature_interfaciale_moyenne() const { return temperature_interfaciale_moyenne_ ; }

  //Rustine
  double E0_ = 0;//volumique
  IJK_Field_double T_rust_;
  void compute_T_rust(const IJK_Field_vector3_double& velocity);

  int type_temperature_convection_form_ = 0;

  OBS_PTR(IJK_FT_cut_cell) ref_ijk_ft_cut_cell_;
  IJK_Splitting::Localisation localisation_ = IJK_Splitting::ELEM;

  IJK_Field_double div_rho_cp_T_;
  IJK_Field_double lambda_;

  FixedVector<Cut_cell_double, 3> cut_cell_flux_diffusion_;
  FixedVector<Cut_cell_double, 3> cut_cell_flux_convection_;

  Cut_cell_conv_scheme cut_cell_conv_scheme_;

  double flux_interfacial_moyen_;
  double temperature_interfaciale_moyenne_;
  Facettes_data coord_facettes_;
  Facettes_data interfacial_temperature_;
  DoubleTabFT interfacial_phin_ai_;

  IJK_Field_vector3_double current_fluxes_conv_;
  IJK_Field_vector3_double current_fluxes_diff_;
  IJK_Field_vector3_double RK3_F_fluxes_conv_;
  IJK_Field_vector3_double RK3_F_fluxes_diff_;
  int runge_kutta_fluxes_convection_ = 0;
  int runge_kutta_fluxes_diffusion_ = 0;
  int runge_kutta_fluxes_pas_de_correction_convection_ = 0;
  int runge_kutta_fluxes_pas_de_correction_diffusion_ = 0;
  IJK_Field_int cellule_rk_restreint_conv_l_;
  IJK_Field_int cellule_rk_restreint_conv_v_;
  IJK_Field_int cellule_rk_restreint_diff_l_;
  IJK_Field_int cellule_rk_restreint_diff_v_;

  // Temporary fields, to inspect each step of the time advance
  int postraiter_champs_intermediaires_ = 0;
  std::shared_ptr<IJK_Field_double> temperature_post_dying_;
  std::shared_ptr<IJK_Field_double> temperature_post_regular_;
  std::shared_ptr<IJK_Field_double> temperature_post_convection_;
  std::shared_ptr<IJK_Field_double> temperature_post_diff_regular_;

  Cut_cell_convection_auxiliaire convective_correction_;
  Cut_cell_diffusion_auxiliaire diffusive_correction_;

  int deactivate_diffusion_interface_ = 0;

  int runge_kutta_restriction_leniency_convection_ = 0;
  int runge_kutta_restriction_leniency_diffusion_ = 0;

  // Champ IJK_Field notant les cellules parcouru lors d'un traitement,
  // c'est-a-dire pour eviter de recalculer plusieurs fois les memes cases lors du calculs des flux.
  IJK_Field_int treatment_count_;

  // Compteur du dernier traitement effectue dans treatment_count_
  int new_treatment_ = 0;

  int verbosite_ = 2;
};

#endif /* IJK_Thermal_cut_cell_included */
