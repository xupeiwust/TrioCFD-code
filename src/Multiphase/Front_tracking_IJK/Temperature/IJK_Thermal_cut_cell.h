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

  void euler_time_step(const double timestep) override;
  void rk3_sub_step(const int rk_step, const double total_timestep, const double time) override;

  CutCell_GlobalInfo compute_global_energy_cut_cell(Cut_field_scalar& cut_field_temperature, bool next);
  CutCell_GlobalInfo compute_d_global_energy_cut_cell(Cut_field_scalar& cut_field_d_temperature, bool next);
  double compute_global_energy() override
  {
    Cerr << "I want to make sure the function compute_global_energy is not used." << finl;
    Process::exit();
    return 0.;
  }
  CutCell_GlobalInfo compute_Tmin_cut_cell(Cut_field_scalar& cut_field_temperature, bool next);
  CutCell_GlobalInfo compute_Tmax_cut_cell(Cut_field_scalar& cut_field_temperature, bool next);
  void calculer_flux_interface();

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

  void compute_interfacial_temperature2(ArrOfDouble& interfacial_temperature, ArrOfDouble& flux_normal_interp) override;

  void calculer_dT_cut_cell(const Cut_field_vector& cut_field_total_velocity);
  void compute_temperature_convection_cut_cell(const Cut_field_vector& cut_field_total_velocity);
  void add_temperature_diffusion() override;
  void compute_diffusion_increment() override;
  void correct_temperature_for_eulerian_fluxes() override { ; };
  double get_temperature_interfaciale_moyenne() const { return temperature_interfaciale_moyenne_ ; }

  //Rustine
  double E0_ = 0;//volumique
  IJK_Field_double T_rust_;
  void compute_T_rust(const FixedVector<IJK_Field_double, 3>& velocity);

  int type_temperature_convection_form_ = 0;

  REF(IJK_FT_cut_cell) ref_ijk_ft_cut_cell_;

  IJK_Field_double div_rho_cp_T_;
  IJK_Field_double lambda_;

  IJK_Field_double div_coeff_grad_T_volume_temp_;

  Cut_field_scalar cut_field_temperature_;
  Cut_field_scalar cut_field_RK3_F_temperature_;
  Cut_cell_vector cut_cell_flux_diffusion_;
  Cut_cell_vector cut_cell_flux_convection_;
  Cut_field_scalar cut_field_div_coeff_grad_T_volume_;
  Cut_field_scalar cut_field_div_coeff_grad_T_volume_temp_;
  Cut_field_scalar cut_field_d_temperature_;

  Cut_cell_conv_scheme cut_cell_conv_scheme_;
  FixedVector<FixedVector<IJK_Field_double, 3>, 2> temperature_face_;
  FixedVector<FixedVector<IJK_Field_double, 3>, 2> temperature_face_ft_;

  double flux_interfacial_moyen_;
  double temperature_interfaciale_moyenne_;
  Facettes_data coord_facettes_;
  Facettes_data interfacial_temperature_;
  DoubleTabFT interfacial_phin_ai_;

  IJK_Field_double temperature_debut_sous_pas_;
  Cut_field_scalar cut_field_temperature_debut_sous_pas_;

  IJK_Field_int cellule_rk_restreint_;

  // Temporary fields, to inspect each step of the time advance
  int postraiter_champs_intermediaires_ = 0;
  IJK_Field_double temperature_post_dying_;
  IJK_Field_double temperature_post_regular_;
  IJK_Field_double temperature_post_convection_;
  IJK_Field_double temperature_post_diff_regular_;
  Cut_field_scalar cut_field_temperature_post_dying_;
  Cut_field_scalar cut_field_temperature_post_regular_;
  Cut_field_scalar cut_field_temperature_post_convection_;
  Cut_field_scalar cut_field_temperature_post_diff_regular_;

  Cut_cell_convection_auxiliaire convective_correction_;
  Cut_cell_diffusion_auxiliaire diffusive_correction_;

  int deactivate_diffusion_interface_ = 0;
  ETALEMENT_DIFFUSION etalement_diffusion_ = ETALEMENT_DIFFUSION::AUCUN_ETALEMENT;

  // Champ IJK_Field notant les cellules parcouru lors d'un traitement,
  // c'est-a-dire pour eviter de recalculer plusieurs fois les memes cases lors du calculs des flux.
  IJK_Field_int treatment_count_;

  // Compteur du dernier traitement effectue dans treatment_count_
  int new_treatment_;

  int verbosite_ = 2;
};

#endif /* IJK_Thermal_cut_cell_included */
