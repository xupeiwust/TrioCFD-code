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
// File      : IJK_Thermal_base.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_Thermal_base_included
#define IJK_Thermal_base_included

#include <IJK_Field.h>
#include <IJK_Field_vector.h>
#include <Objet_U.h>
#include <Boundary_Conditions_Thermique.h>
#include <Domaine_IJK.h>
#include <Parser.h>
#include <IJK_Lata_writer.h>
#include <Operateur_IJK_elem_conv.h>
#include <Operateur_IJK_elem_diff.h>
#include <OpGradCentre2IJKScalar.h>
#include <OpHessCentre2IJKScalar.h>
#include <OpGradQuickIJKScalar.h>
#include <Ouvrir_fichier.h>
#include <Corrige_flux_FT_base.h>
#include <TRUST_Ref.h>
#include <IJK_FT_Post.h>
#include <IJK_Ghost_Fluid_Fields.h>
#include <IJK_One_Dimensional_Subproblems_Interfaces_Fields.h>
#include <Fluide_Diphasique_IJK.h>


class Probleme_FTD_IJK_base;
class Switch_FT_double;
class IJK_Interfaces;

class IJK_Thermal_base : public Objet_U
{
  Declare_base( IJK_Thermal_base ) ;
  friend class IJK_One_Dimensional_Subproblems;
  friend class IJK_One_Dimensional_Subproblem;
public:
  /*
   * Initialisation
   */
  const Milieu_base& milieu() const;
  Milieu_base& milieu();
  void associer_milieu_base(const Milieu_base&);
  inline Fluide_Diphasique_IJK& milieu_ijk() { return ref_cast(Fluide_Diphasique_IJK, milieu()); }
  inline const Fluide_Diphasique_IJK& milieu_ijk() const { return ref_cast(Fluide_Diphasique_IJK, milieu()); }

  virtual void set_param(Param& param);
  virtual int initialize(const Domaine_IJK& splitting, const int idx);
  virtual int initialize_switch(const Domaine_IJK& splitting, const int idx);
  virtual void update_thermal_properties();
  double compute_timestep(const double timestep,
                          const double dxmin);
  void set_fichier_reprise(const char *lataname);
  const Nom& get_fichier_reprise() const { return fichier_reprise_temperature_; }
  void associer(const Probleme_FTD_IJK_base& ijk_ft);
  void associer_post(const IJK_FT_Post& ijk_ft_post);
  void associer_switch(const Switch_FT_double& ijk_ft_switch);
  void associer_interface_intersections(const Intersection_Interface_ijk_cell& intersection_ijk_cell,
                                        const Intersection_Interface_ijk_face& intersection_ijk_face);
  void associer_ghost_fluid_fields(const IJK_Ghost_Fluid_Fields& ghost_fluid_fields);
  void retrieve_ghost_fluid_params(int& compute_distance, int& compute_curvature, int& n_iter_distance, int& avoid_gfm_parallel_calls);
  void get_boundary_fluxes(IJK_Field_local_double& boundary_flux_kmin, IJK_Field_local_double& boundary_flux_kmax);
  virtual void euler_time_step(const double timestep);
  virtual void rk3_sub_step(const int rk_step,
                            const double total_timestep,
                            const double time);
  virtual void sauvegarder_temperature(Nom& lata_name, int idx, const int& stop=0);

  double compute_global_energy(const IJK_Field_double& temperature);
  virtual double compute_global_energy()
  {
    return compute_global_energy(*temperature_); // changes the attribute global_energy [J/m3]
  }
  int calculer_k_pour_bord(const IJK_Field_double& temperature, const bool bord_kmax);
  int calculer_flux_thermique_bord(const IJK_Field_double& temperature,
                                   const double lambda_de_t_paroi,
                                   const double T_paroi_impose,
                                   IJK_Field_local_double& flux_bord,
                                   const bool bord_kmax);
  int imposer_flux_thermique_bord(const IJK_Field_double& temperature,
                                  const double flux_paroi_impose,
                                  IJK_Field_local_double& flux_bord,
                                  const bool bord_kmax);
  virtual int get_first_step_thermals_post() { return 0; };
  void set_latastep_reprise(const int latastep)
  {
    latastep_reprise_ = latastep;
  }
  const int& get_latastep_reprise() const
  {
    return latastep_reprise_;
  }
  const int& get_latastep_reprise_ini() const
  {
    return latastep_reprise_;
  }
  /*
   * Getters and setters
   */
  double get_rhocp_l() const;
  double get_rhocp_v() const;
  const int& get_rank() const { return rang_; };
  const std::shared_ptr<IJK_Field_double>& get_temperature() const
  {
    return temperature_ ;
  }
  const IJK_Field_double& get_temperature_before_extrapolation() const
  {
    return temperature_before_extrapolation_ ;
  }
  IJK_Field_double& get_temperature_ft()
  {
    return temperature_ft_ ;
  }
  const IJK_Field_vector3_double& get_grad_T() const
  {
    return grad_T_ ;
  }
  IJK_Field_double& set_temperature()
  {
    return *temperature_ ;
  }
  const IJK_Field_double& get_temperature_ana() const
  {
    return temperature_ana_ ;
  }
  const IJK_Field_double& get_ecart_t_ana() const
  {
    return ecart_t_ana_ ;
  }
  const IJK_Field_double& get_ecart_t_ana_rel() const
  {
    return ecart_t_ana_rel_;
  }
  const IJK_Field_double& get_div_lambda_grad_T() const
  {
    return div_coeff_grad_T_ ;
  }
  const std::shared_ptr<IJK_Field_double>& get_div_lambda_grad_T_volume() const
  {
    return div_coeff_grad_T_volume_ ;
  }
  const IJK_Field_double& get_u_T_convective() const
  {
    return u_T_convective_;
  }
  const IJK_Field_double& get_u_T_convective_volume() const
  {
    return u_T_convective_volume_;
  }
  const IJK_Field_double& get_eulerian_distance_ft() const
  {
    return ghost_fluid_fields_->get_eulerian_distance_ft();
    // return eulerian_distance_ft_;
  }
  const IJK_Field_double& get_eulerian_curvature_ft() const
  {
    return ghost_fluid_fields_->get_eulerian_curvature_ft();
    // return eulerian_curvature_ft_ ;
  }
  const IJK_Field_double& get_interfacial_area_ft() const
  {
    return ghost_fluid_fields_->get_eulerian_interfacial_area_ft();
    // return eulerian_interfacial_area_ft_;
  }
  const IJK_Field_double& get_grad_T_interface_ft() const
  {
    return eulerian_grad_T_interface_ft_;
  }
  const IJK_Field_double& get_eulerian_compo_connex_ft() const;
  const IJK_Field_double& get_eulerian_compo_connex_ghost_ft() const;
  const IJK_Field_double& get_eulerian_compo_connex_from_interface_ft() const;
  const IJK_Field_double& get_eulerian_compo_connex_from_interface_ghost_ft() const;

  const IJK_Field_double& get_eulerian_compo_connex_ns() const;
  const IJK_Field_double& get_eulerian_compo_connex_ghost_ns() const;
  const IJK_Field_double& get_eulerian_compo_connex_from_interface_ns() const;
  const IJK_Field_double& get_eulerian_compo_connex_from_interface_ghost_ns() const;
  const IJK_Field_int& get_eulerian_compo_connex_int_from_interface_ns() const;
  const IJK_Field_int& get_eulerian_compo_connex_int_from_interface_ghost_ns() const;

  const IJK_Field_double& get_eulerian_distance_ns() const
  {
    return ghost_fluid_fields_->get_eulerian_distance_ns();
    // return eulerian_distance_ns_;
  }
  const IJK_Field_double& get_eulerian_curvature_ns() const
  {
    return ghost_fluid_fields_->get_eulerian_curvature_ns();
    // return eulerian_curvature_ns_ ;
  }
  const IJK_Field_double& get_interfacial_area_ns() const
  {
    return ghost_fluid_fields_->get_eulerian_interfacial_area_ns();
    // return eulerian_interfacial_area_ns_;
  }
  const IJK_Field_double& get_grad_T_interface_ns() const
  {
    return eulerian_grad_T_interface_ns_;
  }
  const IJK_Field_double& get_eulerian_rising_velocities() const
  {
    return *eulerian_rising_velocities_;
  }
  const IJK_Field_double& get_temperature_adim_bulles() const
  {
    return temperature_adim_bulles_;
  }
  const IJK_Field_vector3_double& get_gradient_temperature() const
  {
    return grad_T_ ;
  }
  const IJK_Field_vector3_double& get_gradient_temperature_elem() const
  {
    return grad_T_elem_ ;
  }
  const IJK_Field_vector3_double& get_gradient_temperature_elem_smooth() const
  {
    if (smooth_grad_T_elem_)
      return grad_T_elem_smooth_;
    else
      return dummy_double_vect_;
  }
  const IJK_Field_vector3_double& get_tangential_gradient_temperature_elem_smooth() const
  {
    if (smooth_grad_T_elem_)
      return grad_T_elem_tangential_;
    else
      return dummy_double_vect_;
  }
  const IJK_Field_double& get_temperature_elem_smooth() const
  {
    if (smooth_grad_T_elem_)
      return temperature_gaussian_filtered_;
    else
      return dummy_double_field_;
  }
  const IJK_Field_vector3_double& get_normal_vector_ns() const
  {
    return ghost_fluid_fields_->get_eulerian_normal_vectors_ns();
    // return eulerian_normal_vectors_ns_;
  }
  const IJK_Field_vector3_double& get_normal_vector_ns_normed() const
  {
    return ghost_fluid_fields_->get_eulerian_normal_vectors_ns_normed();
    // return eulerian_normal_vectors_ns_normed_;
  }
  const IJK_Field_vector3_double& get_normal_vector_ft() const
  {
    return ghost_fluid_fields_->get_eulerian_normal_vectors_ft();
    // return eulerian_normal_vectors_ft_;
  }
  const IJK_Field_vector3_double& get_hessian_diag_temperature_elem() const
  {
    return hess_diag_T_elem_ ;
  }
  const IJK_Field_vector3_double& get_hessian_cross_temperature_elem() const
  {
    return hess_cross_T_elem_ ;
  }
  const IJK_Field_vector3_double& get_bary() const
  {
    return ghost_fluid_fields_->get_eulerian_facets_barycentre_ft();
    // return eulerian_facets_barycentre_ft_;
  }
  const int& get_ghost_fluid_flag() const
  {
    return ghost_fluid_;
  };
  const int& get_ghost_cells() const
  {
    return ghost_cells_;
  };
  const int& get_debug() const
  {
    return debug_;
  };
  virtual const IJK_Field_double& get_temperature_cell_neighbours_debug() const
  {
    return dummy_double_field_; //dummy
  }
  virtual const IJK_Field_double& get_temperature_cell_neighbours() const
  {
    return dummy_double_field_; //dummy
  }
  virtual const IJK_Field_int& get_cell_neighbours_corrected() const
  {
    return dummy_int_field_; //dummy
  }
  virtual const IJK_Field_double& get_neighbours_temperature_colinearity_weighting() const
  {
    return dummy_double_field_; //dummy
  }
  virtual const IJK_Field_double& get_debug_lrs_cells() const
  {
    return dummy_double_field_; //dummy
  };
  virtual int get_disable_post_processing_probes_out_files() const
  {
    return 1;
  };
  virtual const IJK_Field_vector3_double& get_cell_faces_corrected_diffusive() const
  {
    return dummy_double_vect_; //dummy
  }
  virtual const IJK_Field_vector3_double& get_cell_faces_corrected_convective() const
  {
    return dummy_double_vect_; //dummy
  }
  virtual const IJK_Field_vector3_int& get_cell_faces_corrected_bool() const
  {
    return dummy_int_vect_;
  }
  virtual const IJK_Field_vector3_int& get_cell_faces_neighbours_corrected_diag_bool() const
  {
    return dummy_int_vect_;
  }
  virtual const IJK_Field_vector3_int& get_cell_faces_neighbours_corrected_all_bool() const
  {
    return dummy_int_vect_;
  }
  virtual const IJK_Field_vector3_int& get_cell_faces_neighbours_corrected_min_max_bool() const
  {
    return dummy_int_vect_;
  }
  virtual const IJK_Field_vector3_double& get_cell_faces_neighbours_corrected_velocity_temperature() const
  {
    return dummy_double_vect_;
  }
  virtual const IJK_Field_vector3_double& get_cell_faces_neighbours_corrected_convective() const
  {
    return dummy_double_vect_;
  }
  virtual const IJK_Field_vector3_double& get_cell_faces_neighbours_corrected_convective_frame_of_ref() const
  {
    return dummy_double_vect_;
  }
  virtual const IJK_Field_vector3_double& get_cell_faces_neighbours_corrected_diffusive() const
  {
    return dummy_double_vect_;
  }
  virtual const IJK_Field_vector3_double& get_neighbours_faces_weighting_colinearity() const
  {
    return dummy_double_vect_;
  }
  virtual const IJK_Field_int& get_cell_neighbours_corrected_trimmed() const
  {
    return dummy_int_field_;
  }
  virtual const IJK_Field_double& get_probe_collision_debug_field() const
  {
    return dummy_double_field_; //dummy
  }
  const IJK_Field_vector3_double& get_rho_cp_u_T_convective_fluxes() const
  {
    if (store_flux_operators_for_energy_balance_)
      return rho_cp_u_T_convective_raw_;
    else
      return dummy_double_vect_;
  }
  const IJK_Field_vector3_double& get_div_coeff_grad_T_diffusive_fluxes() const
  {
    if (store_flux_operators_for_energy_balance_)
      return div_coeff_grad_T_raw_;
    else
      return dummy_double_vect_;
  }
  virtual const IJK_Field_vector3_double& get_interfacial_heat_flux_dispatched() const
  {
    return dummy_double_vect_;
  }
  virtual const IJK_Field_vector3_double& get_interfacial_heat_flux_contrib() const
  {
    return dummy_double_vect_;
  }
  virtual const IJK_Field_vector3_double& get_interfacial_heat_flux_current() const
  {
    return dummy_double_vect_;
  }

  virtual double get_modified_time();
  void get_rising_velocities_parameters(int& compute_rising_velocities,
                                        int& fill_rising_velocities,
                                        int& use_bubbles_velocities_from_interface,
                                        int& use_bubbles_velocities_from_barycentres);

  virtual double get_rho_cp_u_ijk(const IJK_Field_double& vx, int i, int j, int k) const;
  virtual double get_div_lambda_ijk(int i, int j, int k) const { return 0; };
  virtual double compute_temperature_dimensionless_theta_mean(const IJK_Field_double& vx);

  const char * get_fichier_sauvegarde() const
  {
    return fichier_reprise_temperature_;
  }
  void set_fichier_sauvegarde(const char *lataname)
  {
    fichier_reprise_temperature_ = lataname;
  }
  virtual void set_field_T_ana();

  void euler_rustine_step(const double timestep, const double dE);
  void rk3_rustine_sub_step(const int rk_step, const double total_timestep,
                            const double fractionnal_timestep, const double time, const double dE);
  virtual int& get_conserv_energy_global() { return conserv_energy_global_; };
  const double& get_E0() const { return E0_; };

  void compute_dT_rustine(const double dE);
  void compute_T_rust(const IJK_Field_vector3_double& velocity);

  virtual void calculer_ecart_T_ana();
  virtual void compute_interfacial_temperature2(
    ArrOfDouble& interfacial_temperature,
    ArrOfDouble& flux_normal_interp); //const ;
#if 0
  void ecrire_reprise_thermique(SFichier& fichier);
#endif
  virtual void compute_ghost_cell_numbers_for_subproblems(const Domaine_IJK& splitting, int ghost_init) { ghost_cells_ = ghost_init; };

  void compute_eulerian_distance();
  void compute_eulerian_curvature();
  void compute_eulerian_curvature_from_interface();

  virtual void update_intersections() {  }
  virtual void clean_ijk_intersections() {  }
  virtual void post_process_thermal_wake_slices(const Nom& local_quantities_thermal_slices_time_index_folder) { ; };
  virtual void post_process_thermal_downstream_lines(const Nom& local_quantities_thermal_lines_time_index_folder) { ; };
  virtual void set_thermal_subresolution_outputs(const Nom& interfacial_quantities_thermal_probes,
                                                 const Nom& shell_quantities_thermal_probes,
                                                 const Nom& overall_bubbles_quantities,
                                                 const Nom& local_quantities_thermal_probes_time_index_folder) { }
  virtual void compute_temperature_init();
  virtual void recompute_temperature_init();
  virtual int set_subproblems_interfaces_fields(const Domaine_IJK& splitting) { return 0; };
  void copy_previous_interface_state();
  int post_process_quantities_from_subresolution(const Motcles& liste_post_instantanes,
                                                 const char *lata_name,
                                                 const int latastep);

  virtual void copie_pure_vers_diph_sans_interpolation();
  virtual void echange_pure_vers_diph_cellules_initialement_pures();
  virtual void echange_diph_vers_pure_cellules_finalement_pures();
  virtual void vide_phase_invalide_cellules_diphasiques();
  virtual void remplir_tableau_pure_cellules_diphasiques(bool next_time);

  virtual void compare_fluxes_thermal_subproblems() { }

  void post_process_std_thermal_field(const Motcles& liste_post_instantanes, const char *lata_name, const int latastep, const double current_time, const int idx, const Motcles& tested_names, const Nom& name_field, const Motcle& lata_suffix, const IJK_Field_double& field, std::ostringstream& oss, int& counter, const int& first_thermal_rank = 0);
  int posttraiter_champs_instantanes_thermal(const Motcles& liste_post_instantanes, const char *lata_name, const int latastep, const double current_time, const int idx);

  inline Nom& get_thermal_problem_type() { return thermal_problem_type_; }
  inline int& get_thermal_rank() { return thermal_rank_; }
  inline const Motcles& get_thermal_words() const { return thermal_words_; }
  inline const Motcles& get_thermal_suffix() const { return lata_suffix_; }

  void posttraiter_tous_champs_thermal(Motcles& liste, const int idx) const;
  void ecrire_statistiques_bulles(int reset, const Nom& nom_cas, const double current_time, const ArrOfDouble& surface, const int idx);

  int posttraiter_champs_instantanes_thermal_interface(const Motcles& liste_post_instantanes, const char *lata_name, const int latastep, const double current_time, const int idx);
  int posttraiter_champs_instantanes_thermal_interface_ref(const Motcles& liste_post_instantanes, const char *lata_name, const int latastep, const double current_time, const int idx);
  void thermal_subresolution_outputs(const Nom& interfacial_quantities_thermal_probes, const Nom& shell_quantities_thermal_probes, const Nom& overall_bubbles_quantities, const Nom& local_quantities_thermal_probes_time_index_folder, const Nom& local_quantities_thermal_slices_time_index_folder, const Nom& local_quantities_thermal_lines_time_index_folder);

  static void typer_lire_thermal_equation(OWN_PTR(IJK_Thermal_base)&, Entree&);

protected:
  OBS_PTR(Milieu_base) le_fluide_;
  int thermal_rank_ = 0;
  Nom thermal_problem_type_ = "subresolution";
  Motcles thermal_words_, lata_suffix_;
  enum THERMAL_TYPE {SUBRES, MSUBRES, ONEFLUID, ONEFLUIDE, CUTCELL};

  void compute_cell_volume();
  void compute_min_cell_delta();
  void compute_cell_diagonal(const Domaine_IJK& splitting);

  virtual void lire_temperature(const Domaine_IJK& splitting, int idx);

  void calculer_dT(const IJK_Field_vector3_double& velocity);
  virtual void post_process_after_temperature_increment();

  void compute_temperature_convective_fluxes(const IJK_Field_vector3_double& velocity);
  void compute_temperature_convection(const IJK_Field_vector3_double& velocity);

  void compute_boundary_conditions_thermal();
  void compute_temperature_diffusive_fluxes();
  virtual void add_temperature_diffusion();
  virtual void compute_diffusion_increment()=0;

  virtual void correct_temperature_for_eulerian_fluxes()=0;
  virtual void store_temperature_before_extrapolation() { }
  virtual void correct_temperature_increment_for_interface_leaving_cell() { }
  void enforce_zero_value_eulerian_distance();
  void enforce_zero_value_eulerian_curvature();
  void enforce_max_value_eulerian_curvature();

  void compute_eulerian_grad_T_interface(const int on_splitting_ns=0);
  void propagate_eulerian_grad_T_interface();

  void compute_eulerian_temperature_ghost(const int on_splitting_ns=0);

  void compute_eulerian_bounding_box_fill_compo();
  void compute_rising_velocities();
  //  void enforce_zero_value_eulerian_field(IJK_Field_double& eulerian_field);
  //  void enforce_max_value_eulerian_field(IJK_Field_double& eulerian_field);
  //  void enforce_min_value_eulerian_field(IJK_Field_double& eulerian_field);
  void compute_temperature_gradient_elem();
  void compute_temperature_hessian_diag_elem();
  void compute_temperature_hessian_cross_elem();
  virtual void correct_temperature_for_visu() { }
  virtual void correct_operators_for_visu() { }
  virtual void clip_temperature_values() { }
  virtual void clip_min_temperature_values() {  }
  virtual void clip_max_temperature_values() { }
  virtual void compute_thermal_subproblems() { }
  virtual void compute_convective_diffusive_fluxes_face_centre() { }
  virtual void compute_convective_fluxes_face_centre() { }
  virtual void compute_diffusive_fluxes_face_centre() { }
  virtual void prepare_ij_fluxes_k_layers() { }
  virtual void compute_temperature_cell_centres(const int first_corr) { }
  virtual void set_zero_temperature_increment() { }
  virtual void clean_thermal_subproblems() { }

  void calculer_gradient_temperature(const IJK_Field_double& temperature,
                                     IJK_Field_vector3_double& grad_T);
  void calculer_energies(double& E_liq_pure,
                         double& E_liq,
                         double& E_vap_pure,
                         double& E_vap,
                         double& E_mixt,
                         double& E_tot);

  void source_callback();
  void calculer_temperature_physique_T(const IJK_Field_double&  vx, const double dTm);
  void calculer_temperature_adim_bulles();
  void add_temperature_source();
  void calculer_Nusselt(const IJK_Field_double& vx);
  void calculer_temperature_adimensionnelle_theta(const IJK_Field_double&  vx, const double qw);
  void calculer_source_temperature_ana();
  virtual double compute_rho_cp_u_mean(const IJK_Field_double& vx);
  double compute_variable_wall_temperature(const int kmin, const int kmax);

  void force_upstream_temperature(IJK_Field_double& temperature, double T_imposed,
                                  const IJK_Interfaces& interfaces, double nb_diam, int upstream_dir,
                                  int gravity_dir, int upstream_stencil);

  virtual void enforce_periodic_temperature_boundary_value() {  }

  int debug_ = 0;
  int latastep_reprise_ = 0;
  int latastep_reprise_ini_ = 0;
  /*
   * Patch to conserve energy
   */
  double E0_ = 0.; //volumique
  IJK_Field_double T_rust_;
  IJK_Field_double d_T_rustine_; // Temperature increment to conserve the energy.
  IJK_Field_double RK3_F_rustine_; // Temporary storage for substeps in the RK3 algorithm for the rustine calculation.

  OBS_PTR(Probleme_FTD_IJK_base) ref_ijk_ft_;
  OBS_PTR(IJK_FT_Post) ref_ijk_ft_post_;
  OBS_PTR(Switch_FT_double) ref_ijk_ft_switch_;
  OBS_PTR(Intersection_Interface_ijk_cell) ref_intersection_ijk_cell_;
  OBS_PTR(Intersection_Interface_ijk_face) ref_intersection_ijk_face_;
  OWN_PTR(Corrige_flux_FT_base) corrige_flux_;
  const IJK_Field_double& get_IJK_field(const Nom& nom) const;
  int rang_ = 0; // default value, used as an index for the list of thermal sub-problems

  /*
   * Physical parameters and inputs
   */
  double dt_fo_ = 1.e20;
  double fo_ = 1.; // Fourier number
  double cp_liquid_ = -123., cp_vapour_ = -123., lambda_liquid_ = -123., lambda_vapour_ = -123.;
  int single_phase_ = 1.;
  double uniform_lambda_ = 0.;
  double uniform_alpha_ = 0.;
  double prandtl_number_ = 0.;
  /*
   * Initialisation (B.Cs, expression)
   */
  Boundary_Conditions_Thermique boundary_conditions_;
  IJK_Field_local_double boundary_flux_kmin_;
  IJK_Field_local_double boundary_flux_kmax_;
  Nom expression_T_init_ = "??";
  double upstream_temperature_ = -1.1e20;
  double nb_diam_upstream_ = 0;
  int side_temperature_ = 0;
  int stencil_side_ = 0;

  /*
   * Settings to resume a calculation
   */
  Nom fichier_reprise_temperature_ = "??";
  int timestep_reprise_temperature_ = 1;
  int rank_reprise_temperature_ = -1;

  /*
   * Source of temperature, wall heating,
   * integral per boundary (Tryggvason)
   */
  /*
   * Ceci est une initialisation des derivees des temperatures moyenne de chaque phase
   * Il n'est peut-etre pas pertinent de les mettre ici
   */
  Nom expression_source_temperature_;
  Nom type_T_source_ = "??";
  int lambda_variable_ = 0; // terme source variable
  int wall_flux_ = 0;
  IJK_Field_double source_temperature_;
  IJK_Field_double source_temperature_v_;
  IJK_Field_double source_temperature_l_;
  IJK_Field_double d_source_Tl_;
  IJK_Field_double d_source_Tv_;
  double dTl_ = 0.;
  double dTv_ = 1.;
  double Tl_ = 0.;
  double Tv_ = 1.;
  double Tv0_ = 1.; // Serait-ce plutot Tref (une temperature de reference pour reconstruire le champ dim??)
  double kl_ = -100000000000000.;
  double kv_ = -200000000000000.;
  double T0v_ = 1.;
  double T0l_ = 0.;
  IJK_Field_double source_temperature_ana_;
  IJK_Field_double ecart_source_t_ana_;

  /*
   * Dimensionless temperature
   */
  IJK_Field_double temperature_physique_T_;
  IJK_Field_double temperature_adimensionnelle_theta_;
  IJK_Field_double temperature_adim_bulles_;

  /*
   * Storage for operators & time scheme
   */
  int diff_temperature_negligible_ = 0;
  int conv_temperature_negligible_ = 0;

  /*
   * type_temperature_convection_op_:
   * 1 : Quick
   * 2 : Centre2
   */
  Operateur_IJK_elem_conv temperature_convection_op_;
  Operateur_IJK_elem_diff temperature_diffusion_op_;
  IJK_Field_vector3_double div_coeff_grad_T_raw_;
  std::shared_ptr<IJK_Field_double> div_coeff_grad_T_volume_;
  IJK_Field_double div_coeff_grad_T_;
  IJK_Field_vector3_double rho_cp_u_T_convective_raw_;
  IJK_Field_double u_T_convective_volume_;
  IJK_Field_double u_T_convective_;
  OpGradFluxQuickIJKScalar_double temperature_grad_flux_op_quick_;
  OpGradCentre2IJKScalar_double temperature_grad_op_centre_;
  OpHessCentre2IJKScalar_double temperature_hess_op_centre_;
  OpHessFluxCentre2IJKScalar_double temperature_hess_flux_op_centre_;


  /*
   * Fields
   */
  double vol_ = 0.;
  double min_delta_xyz_ = 0.;
  double cell_diagonal_ = 0.;
  int ghost_cells_ = 4;
  IJK_Field_double rho_cp_;
  IJK_Field_double rho_cp_T_;
  std::shared_ptr<IJK_Field_double> temperature_;
  IJK_Field_double temperature_for_ini_per_bubble_;
  IJK_Field_double temperature_before_extrapolation_;
  std::shared_ptr<IJK_Field_double> d_temperature_; // Temperature increment.
  std::shared_ptr<IJK_Field_double> RK3_F_temperature_; // Temporary storage for substeps in the RK3 algorithm.
  IJK_Field_vector3_double storage_; // Temporary storage for fluxes calculation.
  int calculate_local_energy_ = 0;
  int conserv_energy_global_ = 0;

  /*
   * Fields FT
   * TODO: Clean FT_fields and redistribute curvature, interfacial_area if necessary
   */
  IJK_Field_double temperature_ft_;

  /*
   * Post-processing
   */
  Nom expression_T_ana_ = "??";
  Motcles liste_post_instantanes_; // liste des champs instantanes a postraiter
  IJK_Field_double temperature_ana_, ecart_t_ana_, ecart_t_ana_rel_;
  IJK_Field_vector3_double grad_T_;
  int calulate_grad_T_ = 0;
  int rho_cp_post_ = 0;

  /*
   * For Ghost fluid method & Subresolution or Post-processing
   */
  int ghost_fluid_ = 0;
  int n_iter_distance_ = 3;
  int gfm_recompute_field_ini_ = 1;
  int gfm_zero_neighbour_value_mean_ = 0;
  int gfm_vapour_mixed_only_ = 1;
  int gfm_vapour_liquid_vapour_ = 0;
  int gfm_smooth_factor_ = 20;
  int avoid_gfm_parallel_calls_ = 0;

  int compute_distance_ = 0;
  int compute_curvature_ = 0;
  int compute_grad_T_interface_ = 0;
  /*
   * TODO: Move fields and avoid redundancies in IJK_Interfaces
   * Clean FT_fields
   */
  const DoubleTab * bounding_box_ = nullptr;
  const DoubleTab * min_max_larger_box_ = nullptr;
  const IJK_Ghost_Fluid_Fields * ghost_fluid_fields_ = nullptr;
  const IJK_Field_double * eulerian_distance_ft_ = nullptr;
  const IJK_Field_double * eulerian_distance_ns_ = nullptr;
  const IJK_Field_vector3_double * eulerian_normal_vectors_ft_ = nullptr;
  const IJK_Field_vector3_double * eulerian_facets_barycentre_ft_ = nullptr;
  const IJK_Field_vector3_double * eulerian_normal_vectors_ns_ = nullptr;
  const IJK_Field_vector3_double * eulerian_normal_vectors_ns_normed_ = nullptr;
  const IJK_Field_vector3_double * eulerian_facets_barycentre_ns_ = nullptr;
  const IJK_Field_double * eulerian_curvature_ft_ = nullptr;
  const IJK_Field_double * eulerian_curvature_ns_ = nullptr;
  const IJK_Field_double * eulerian_interfacial_area_ft_ = nullptr;
  const IJK_Field_double * eulerian_interfacial_area_ns_ = nullptr;
  IJK_Field_double eulerian_grad_T_interface_ft_;
  IJK_Field_double eulerian_grad_T_interface_ns_;
  int compute_grad_T_elem_ = 0;
  IJK_Field_vector3_double grad_T_elem_;
  int smooth_grad_T_elem_ = 0;
  IJK_Field_vector3_double grad_T_elem_smooth_;
  /*
   * hess(T) = grad(grad(T))
   * Only 6 coefficients in
   * cartesian coordinate system
   * | * * * |
   * | - * * |
   * | - - * |
   *   FixedVector<IJK_Field_double, 6> grad_grad_T_elem_;
   */
  IJK_Field_vector3_double hess_diag_T_elem_;
  IJK_Field_vector3_double hess_cross_T_elem_;
  IJK_Field_vector3_double facets_barycentre;
  IJK_Field_double d_temperature_uncorrected_;
  IJK_Field_double div_coeff_grad_T_volume_uncorrected_;
  Operateur_IJK_elem_conv temperature_convection_op_uncorrected_;
  Operateur_IJK_elem_diff temperature_diffusion_op_uncorrected_;
  int compute_hess_T_elem_ = 0;
  int compute_hess_diag_T_elem_ = 0;
  int compute_hess_cross_T_elem_ = 0;

  int mixed_cells_number_ = 0;
  void compute_mixed_cells_number(const IJK_Field_double& indicator);
  int compute_eulerian_compo_ = 1;

//  IJK_Field_double eulerian_compo_connex_ft_;
//  IJK_Field_double eulerian_compo_connex_ns_;
//  IJK_Field_double eulerian_compo_connex_ghost_ft_;
//  IJK_Field_double eulerian_compo_connex_ghost_ns_;
  int spherical_approx_ = 1;
  int spherical_exact_ = 0;
  const IJK_Field_double * eulerian_compo_connex_ft_ = nullptr;
  const IJK_Field_double * eulerian_compo_connex_ns_ = nullptr;
  const IJK_Field_double * eulerian_compo_connex_ghost_ft_ = nullptr;
  const IJK_Field_double * eulerian_compo_connex_ghost_ns_ = nullptr;
  const IJK_Field_double * eulerian_compo_connex_from_interface_ft_ = nullptr;
  const IJK_Field_double * eulerian_compo_connex_from_interface_ns_ = nullptr;
  const IJK_Field_double * eulerian_compo_connex_from_interface_ghost_ft_ = nullptr;
  const IJK_Field_double * eulerian_compo_connex_from_interface_ghost_ns_ = nullptr;
  const IJK_Field_int * eulerian_compo_connex_from_interface_int_ns_ = nullptr;
  const IJK_Field_int * eulerian_compo_connex_from_interface_ghost_int_ns_ = nullptr;

  int compute_rising_velocities_ = 0;
  int fill_rising_velocities_ = 0;
  int use_bubbles_velocities_from_interface_ = 0;
  int use_bubbles_velocities_from_barycentres_ = 0;

  const Vecteur3 * liquid_velocity_ = nullptr;
  const Vecteur3 * rising_velocity_overall_ = nullptr;
  const ArrOfDouble * rising_velocities_ = nullptr;
  const ArrOfDouble * rising_velocities_from_barycentres_ = nullptr;
  const DoubleTab * rising_vectors_ = nullptr;
  const DoubleTab * rising_vectors_from_barycentres_ = nullptr;
  const IJK_Field_double * eulerian_rising_velocities_ = nullptr;
  const ArrOfDouble * bubbles_volume_ = nullptr;
  const DoubleTab * bubbles_barycentre_ = nullptr;
  const DoubleTab * bubbles_barycentres_old_ = nullptr;
  const DoubleTab * bubbles_barycentres_new_ = nullptr;

  IJK_Field_vector3_int dummy_int_vect_;
  IJK_Field_vector3_double dummy_double_vect_;
  IJK_Field_int dummy_int_field_;
  IJK_Field_double dummy_double_field_;

  int store_flux_operators_for_energy_balance_ = 0;
  int disable_relative_velocity_energy_balance_ = 0;

  IJK_One_Dimensional_Subproblems_Interfaces_Fields thermal_local_subproblems_interfaces_fields_;
  IJK_Field_double temperature_gaussian_filtered_;
  IJK_Field_double tmp_smoothing_field_;
  IJK_Field_vector3_double grad_T_elem_tangential_;
  // IJK_Field_vector3_double grad_T_elem_gaussian_filtered_;
  int smoothing_numbers_ = 1;
  int smoothing_remove_normal_compo_ = 0;
  int smoothing_use_unique_phase_ = 0;
  double direct_smoothing_factors_[7] = {1.,1.,1.,1.,1.,1.,2.};
  double gaussian_smoothing_factors_[3][3][3] = {{{1,2,1},
      {2,4,2},
      {1,2,1}
    },
    { {1,4,1},
      {4,8,4},
      {2,4,2}
    },
    { {1,2,1},
      {2,4,2},
      {1,2,1}
    }
  };
  double sharpen_smoothing_factors_[3][3][3] = {{{0,0,0},
      {0,-1,0},
      {0,0,0}
    },
    { {0,-1,0},
      {-1,8,-1},
      {0,-1,0}
    },
    { {0,0,0},
      {0,-1,0},
      {0,0,0}
    }
  };
};

#endif /* IJK_Thermal_base_included */
