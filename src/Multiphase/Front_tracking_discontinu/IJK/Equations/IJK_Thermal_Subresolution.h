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

#ifndef IJK_Thermal_Subresolution_included
#define IJK_Thermal_Subresolution_included

#include <IJK_Thermal_base.h>
#include <IJK_Field_vector.h>
#include <IJK_Field.h>
#include <Boundary_Conditions_Thermique.h>
#include <Domaine_IJK.h>
#include <Parser.h>
#include <IJK_Lata_writer.h>
#include <OpConvQuickIJKScalar.h>
#include <OpConvCentre2IJKScalar.h>
#include <Ouvrir_fichier.h>
#include <Corrige_flux_FT_base.h>
#include <TRUST_Ref.h>
#include <Operateur_IJK_elem_diff_base.h>
#include <OpConvAmontIJK.h>
#include <OpConvDiscQuickIJKScalar.h>
#include <OpConvCentre4IJK.h>
#include <IJK_One_Dimensional_Subproblems.h>
#include <IJK_Finite_Difference_One_Dimensional_Matrix_Assembler.h>
#include <IJK_SolveSys_FD_thermal.h>
#include <MD_Vector_tools.h>


class IJK_Thermal_Subresolution : public IJK_Thermal_base
{

  Declare_instanciable( IJK_Thermal_Subresolution ) ;
  friend class IJK_One_Dimensional_Subproblems;
  friend class IJK_One_Dimensional_Subproblem;

public :

  void initialize(const Domaine_IJK& splitting, const int idx) override;
  // void sauvegarder_temperature(Nom& lata_name, int idx) override;
  void update_thermal_properties() override;
  void post_process_after_temperature_increment() override;
  void set_param(Param& param) override;
  void compute_ghost_cell_numbers_for_subproblems(const Domaine_IJK& splitting, int ghost_init) override;

  double get_probes_length();
  // Entree& read_fd_solver(Entree& is);
  // void read_fd_solver(const Motcle& mot, Entree& is);
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  void set_thermal_subresolution_outputs(const Nom& interfacial_quantities_thermal_probes,
                                         const Nom& shell_quantities_thermal_probes,
                                         const Nom& overall_bubbles_quantities,
                                         const Nom& local_quantities_thermal_probes_time_index_folder) override;
  void post_process_thermal_downstream_lines(const Nom& local_quantities_thermal_lines_time_index_folder) override;
  void initialise_thermal_dowstreamlines_tabs(std::vector<std::vector<FixedVector<ArrOfInt,3>>>& parameters,
                                              const int& nb_thermal_circles,
                                              const int& nb_thermal_lines);
  void initialise_thermal_dowstreamlines_tabs(std::vector<std::vector<FixedVector<ArrOfDouble,2>>>& parameters,
                                              const int& nb_thermal_circles,
                                              const int& nb_thermal_lines);
  void initialise_thermal_dowstreamlines_tabs(std::vector<std::vector<ArrOfInt>>& parameters,
                                              const int& nb_thermal_circles,
                                              const int& nb_thermal_lines);
  void initialise_thermal_dowstreamlines_tabs(std::vector<std::vector<ArrOfDouble>>& parameters,
                                              const int& nb_thermal_circles,
                                              const int& nb_thermal_lines);
  void initialise_thermal_line_points(const int& line_dir,
                                      ArrOfDouble& linear_coord,
                                      FixedVector<ArrOfDouble,3>& coordinates_line,
                                      double& diameter);
  void find_thermal_line_points_in_procs(std::vector<std::vector<ArrOfInt>>& parameters);
  void find_cocentric_line_coordinates(const int& nb_thermal_circles,
                                       const int& nb_thermal_lines,
                                       const double& diameter_approx,
                                       std::vector<std::vector<FixedVector<ArrOfDouble,2>>>& coordinates_sides);
  void find_points_on_proc(std::vector<std::vector<ArrOfInt>>& is_point_on_proc,
                           std::vector<std::vector<FixedVector<ArrOfInt,3>>>& ijk_indices,
                           const FixedVector<ArrOfDouble,3>& coordinates_line,
                           const std::vector<std::vector<FixedVector<ArrOfDouble,2>>>& coordinates_sides,
                           const int& line_dir);
  void min_max_ldir(const int& dir,
                    const Domaine_IJK& geom,
                    double& min_dir, double& max_dir);
  void interpolate_fields_on_downstream_line(const int& dir,
                                             const int& nb_thermal_circles,
                                             const int& index_circle,
                                             const int& index_line,
                                             const std::vector<std::vector<ArrOfInt>>& is_point_on_proc,
                                             const FixedVector<ArrOfDouble,3>& coordinates_line,
                                             const std::vector<std::vector<FixedVector<ArrOfDouble,2>>>& coordinates_sides,
                                             const IJK_Field_double& field,
                                             const IJK_Field_vector3_double& field_gradient,
                                             const IJK_Field_vector3_double& velocity,
                                             DoubleVect& values,
                                             const int field_type);
  void interpolate_temperature_on_downstream_line(const int& dir,
                                                  const int& nb_thermal_circles,
                                                  const int& index_circle,
                                                  const int& index_line,
                                                  const std::vector<std::vector<ArrOfInt>>& is_point_on_proc,
                                                  const FixedVector<ArrOfDouble,3>& coordinates_line,
                                                  const std::vector<std::vector<FixedVector<ArrOfDouble,2>>>& coordinates_sides,
                                                  DoubleVect& values);
  void interpolate_velocity_on_downstream_line(const int& dir,
                                               const int& nb_thermal_circles,
                                               const int& index_circle,
                                               const int& index_line,
                                               const std::vector<std::vector<ArrOfInt>>& is_point_on_proc,
                                               const FixedVector<ArrOfDouble,3>& coordinates_line,
                                               const std::vector<std::vector<FixedVector<ArrOfDouble,2>>>& coordinates_sides,
                                               DoubleVect& values);
  void interpolate_convective_term_on_downstream_line(const int& dir,
                                                      const int& nb_thermal_circles,
                                                      const int& index_circle,
                                                      const int& index_line,
                                                      const std::vector<std::vector<ArrOfInt>>& is_point_on_proc,
                                                      const FixedVector<ArrOfDouble,3>& coordinates_line,
                                                      const std::vector<std::vector<FixedVector<ArrOfDouble,2>>>& coordinates_sides,
                                                      DoubleVect& values);
  void interpolate_diffusive_term_on_downstream_line(const int& dir,
                                                     const int& nb_thermal_circles,
                                                     const int& index_circle,
                                                     const int& index_line,
                                                     const std::vector<std::vector<ArrOfInt>>& is_point_on_proc,
                                                     const FixedVector<ArrOfDouble,3>& coordinates_line,
                                                     const std::vector<std::vector<FixedVector<ArrOfDouble,2>>>& coordinates_sides,
                                                     DoubleVect& values);
  void interpolate_temperature_increment_on_downstream_line(const int& dir,
                                                            const int& nb_thermal_circles,
                                                            const int& index_circle,
                                                            const int& index_line,
                                                            const std::vector<std::vector<ArrOfInt>>& is_point_on_proc,
                                                            const FixedVector<ArrOfDouble,3>& coordinates_line,
                                                            const std::vector<std::vector<FixedVector<ArrOfDouble,2>>>& coordinates_sides,
                                                            DoubleVect& values);
  void post_processed_fields_on_downstream_line(const Nom& local_quantities_thermal_lines_time_index_folder,
                                                const int& line_dir,
                                                const int& linear_circle_line_index,
                                                const int& index_circle,
                                                const int& index_line,
                                                const double& diameter_approx,
                                                std::vector<std::vector<ArrOfInt>>& is_point_on_proc,
                                                const std::vector<std::vector<FixedVector<ArrOfInt,3>>>& indices_ijk,
                                                const ArrOfDouble& linear_coord,
                                                const FixedVector<ArrOfDouble,3>& coordinates_line,
                                                const std::vector<std::vector<FixedVector<ArrOfDouble,2>>>& coordinates_sides,
                                                const DoubleVect& temperature_line,
                                                const DoubleVect& velocity_line,
                                                const DoubleVect& convective_term_line,
                                                const DoubleVect& diffusive_term_line,
                                                const DoubleVect& temperature_incr_line);
  void post_process_thermal_wake_slices(const Nom& local_quantities_thermal_slices_time_index_folder) override;
  void post_process_thermal_wake_slice(const int& slice,
                                       const double& nb_diam_slice,
                                       const Nom& local_quantities_thermal_slices_time_index_folder);
  double post_process_thermal_wake_slice_index_dir(int& index_dir_local,
                                                   int& index_dir_global,
                                                   int& n_cross_section_1,
                                                   int& n_cross_section_2,
                                                   int& dir,
                                                   const double& nb_diam,
                                                   int upstream_dir,
                                                   int gravity_dir,
                                                   double& diameter);
  void complete_field_thermal_wake_slice_ij_values(int& index_dir_local,
                                                   const int& dir,
                                                   const double& slice_pos,
                                                   FixedVector<IntTab, 2>& ij_indices,
                                                   FixedVector<DoubleTab, 3>& ij_coords,
                                                   const IJK_Field_double& field,
                                                   const IJK_Field_vector3_double& field_gradient,
                                                   const IJK_Field_vector3_double& velocity,
                                                   DoubleTab& values,
                                                   const int field_type,
                                                   const int slice_to_nearest_plane,
                                                   const int compute_indices = 0);
  void complete_field_thermal_wake_slice_ij_indices_coords(const int& slice,
                                                           int& index_dir_local,
                                                           const int& dir,
                                                           const double& slice_pos,
                                                           FixedVector<IntTab, 2>& ij_indices,
                                                           FixedVector<DoubleTab, 3>& ij_coords,
                                                           DoubleTab& values,
                                                           const Nom& local_quantities_thermal_slices_time_index_folder);
  void complete_field_thermal_wake_slice_ij_temperature(const int& slice,
                                                        int& index_dir_local,
                                                        const int& dir,
                                                        const double& slice_pos,
                                                        FixedVector<IntTab, 2>& ij_indices,
                                                        FixedVector<DoubleTab, 3>& ij_coords,
                                                        DoubleTab& values,
                                                        const Nom& local_quantities_thermal_slices_time_index_folder);
  void complete_field_thermal_wake_slice_ij_velocity(const int& slice,
                                                     int& index_dir_local,
                                                     const int& dir,
                                                     const double& slice_pos,
                                                     FixedVector<IntTab, 2>& ij_indices,
                                                     FixedVector<DoubleTab, 3>& ij_coords,
                                                     DoubleTab& velocity_values,
                                                     const Nom& local_quantities_thermal_slices_time_index_folder);
  void complete_field_thermal_wake_slice_ij_convection(const int& slice,
                                                       int& index_dir_local,
                                                       const int& dir,
                                                       const double& slice_pos,
                                                       FixedVector<IntTab, 2>& ij_indices,
                                                       FixedVector<DoubleTab, 3>& ij_coords,
                                                       DoubleTab& values,
                                                       DoubleTab& velocity_values,
                                                       const Nom& local_quantities_thermal_slices_time_index_folder);
  void complete_field_thermal_wake_slice_ij_diffusion(const int& slice,
                                                      int& index_dir_local,
                                                      const int& dir,
                                                      const double& slice_pos,
                                                      FixedVector<IntTab, 2>& ij_indices,
                                                      FixedVector<DoubleTab, 3>& ij_coords,
                                                      DoubleTab& values,
                                                      const Nom& local_quantities_thermal_slices_time_index_folder);
  void complete_field_thermal_wake_slice_ij_temperature_incr(const int& slice,
                                                             int& index_dir_local,
                                                             const int& dir,
                                                             const double& slice_pos,
                                                             FixedVector<IntTab, 2>& ij_indices,
                                                             FixedVector<DoubleTab, 3>& ij_coords,
                                                             DoubleTab& values,
                                                             const Nom& local_quantities_thermal_slices_time_index_folder);
  void post_processed_field_thermal_wake_slice_ij(const int& slice,
                                                  const Nom& local_quantities_thermal_slices_time_index_folder,
                                                  const double& diameter_approx,
                                                  const double& nb_diam_slice,
                                                  const int& n_cross_section_1,
                                                  const int& n_cross_section_2,
                                                  const FixedVector<IntTab, 2> ij_indices,
                                                  const FixedVector<DoubleTab, 3>& ij_coords,
                                                  const DoubleTab& temperature_slice,
                                                  const DoubleTab& velocity_slice,
                                                  const DoubleTab& convection_slice,
                                                  const DoubleTab& diffusion_slice,
                                                  const DoubleTab& temperature_incr_slice);

  const IJK_Field_double& get_debug_lrs_cells() const override
  {
    return debug_LRS_cells_;
  }
  const IJK_Field_double& get_temperature_cell_neighbours_debug() const override
  {
    if (find_temperature_cell_neighbours_ && debug_)
      return temperature_cell_neighbours_debug_;
    else
      return dummy_double_field_;
  }
  const IJK_Field_double& get_temperature_cell_neighbours() const override
  {
    if (find_temperature_cell_neighbours_)
      return temperature_cell_neighbours_;
    else
      return dummy_double_field_;
  }
  const IJK_Field_int& get_cell_neighbours_corrected() const override
  {
    if (find_temperature_cell_neighbours_)
      return neighbours_temperature_to_correct_;
    else
      return dummy_int_field_;
  }
  const IJK_Field_double& get_neighbours_temperature_colinearity_weighting() const override
  {
    if (find_temperature_cell_neighbours_)
      return neighbours_temperature_colinearity_weighting_;
    else
      return dummy_double_field_;
  }
  const IJK_Field_vector3_double& get_cell_faces_corrected_diffusive() const override
  {
    if((convective_flux_correction_ || diffusive_flux_correction_) && store_cell_faces_corrected_)
      return cell_faces_corrected_diffusive_;
    else
      return dummy_double_vect_;
  }
  const IJK_Field_vector3_double& get_cell_faces_corrected_convective() const override
  {
    if((convective_flux_correction_ || diffusive_flux_correction_) && store_cell_faces_corrected_)
      return cell_faces_corrected_convective_;
    else
      return dummy_double_vect_;
  }
  const IJK_Field_vector3_int& get_cell_faces_corrected_bool() const override
  {
    if ((convective_flux_correction_ || diffusive_flux_correction_) && store_cell_faces_corrected_)
      return cell_faces_corrected_bool_;
    else
      return dummy_int_vect_;
  }
  const IJK_Field_vector3_int& get_cell_faces_neighbours_corrected_diag_bool() const override
  {
    if (find_cell_neighbours_for_fluxes_spherical_correction_)
      return cell_faces_neighbours_corrected_diag_bool_;
    else
      return dummy_int_vect_;
  }
  const IJK_Field_vector3_int& get_cell_faces_neighbours_corrected_all_bool() const override
  {
    if (find_reachable_fluxes_)
      return cell_faces_neighbours_corrected_all_bool_;
    else
      return dummy_int_vect_;
  }
  const IJK_Field_vector3_int& get_cell_faces_neighbours_corrected_min_max_bool() const override
  {
    if (find_reachable_fluxes_)
      return cell_faces_neighbours_corrected_min_max_bool_;
    else
      return dummy_int_vect_;
  }
  const IJK_Field_vector3_double& get_cell_faces_neighbours_corrected_convective() const override
  {
    if (use_reachable_fluxes_ && convective_flux_correction_)
      return cell_faces_neighbours_corrected_convective_;
    else
      return dummy_double_vect_;
  }
  const IJK_Field_vector3_double& get_cell_faces_neighbours_corrected_convective_frame_of_ref() const override
  {
    if (use_reachable_fluxes_ && convective_flux_correction_ && store_flux_operators_for_energy_balance_)
      return cell_faces_neighbours_corrected_convective_frame_of_reference_;
    else
      return dummy_double_vect_;
  }
  const IJK_Field_vector3_double& get_cell_faces_neighbours_corrected_velocity_temperature() const override
  {
    if (use_reachable_fluxes_ && convective_flux_correction_ && store_flux_operators_for_energy_balance_)
      return cell_faces_neighbours_corrected_velocity_temperature_;
    else
      return dummy_double_vect_;
  }
  const IJK_Field_vector3_double& get_cell_faces_neighbours_corrected_diffusive() const override
  {
    if (use_reachable_fluxes_ && diffusive_flux_correction_)
      return cell_faces_neighbours_corrected_diffusive_;
    else
      return dummy_double_vect_;
  }
  const IJK_Field_vector3_double& get_neighbours_faces_weighting_colinearity() const override
  {
    if (use_reachable_fluxes_ && neighbours_colinearity_weighting_)
      return neighbours_faces_weighting_colinearity_;
    else
      return dummy_double_vect_;
  }
  const IJK_Field_int& get_cell_neighbours_corrected_trimmed() const override
  {
    if (find_reachable_fluxes_)
      return neighbours_temperature_to_correct_trimmed_;
    else
      return dummy_int_field_;
  }
  const IJK_Field_double& get_probe_collision_debug_field() const override
  {
    if (debug_probe_collision_)
      return probe_collision_debug_field_;
    else
      return dummy_double_field_;
  }
  const IJK_Field_vector3_double& get_interfacial_heat_flux_dispatched() const override
  {
    if (fluxes_correction_conservations_)
      return interfacial_heat_flux_dispatched_;
    else
      return dummy_double_vect_;
  }
  const IJK_Field_vector3_double& get_interfacial_heat_flux_contrib() const override
  {
    if (fluxes_correction_conservations_)
      return interfacial_heat_flux_contrib_;
    else
      return dummy_double_vect_;
  }

  const IJK_Field_vector3_double& get_interfacial_heat_flux_current() const override
  {
    if (fluxes_correction_conservations_)
      return interfacial_heat_flux_current_;
    else
      return dummy_double_vect_;
  }

  int get_disable_post_processing_probes_out_files() const override
  {
    return disable_post_processing_probes_out_files_;
  }
  int get_first_step_thermals_post() override { return first_step_thermals_post_; }
  void set_subproblems_interfaces_fields(const Domaine_IJK& splitting) override;
protected :
  void reset_subresolution_distributed_vectors();
  void compute_thermal_subproblems() override;
  void compute_diffusion_increment() override;
  void correct_temperature_increment_for_interface_leaving_cell() override;
  void correct_any_temperature_fields_for_eulerian_fluxes(IJK_Field_double& temperature);
  void correct_temperature_for_eulerian_fluxes() override;
  void store_temperature_before_extrapolation() override;
  void compare_temperature_fields(const IJK_Field_double& temperature,
                                  const IJK_Field_double& temperature_ana,
                                  IJK_Field_double& error_temperature_ana,
                                  IJK_Field_double& error_temperature_ana_rel);
  void evaluate_total_liquid_absolute_parameter(const IJK_Field_double& field,
                                                double& total_parameter);
  void evaluate_total_liquid_parameter_squared(const IJK_Field_double& field,
                                               double& total_parameter);
  void correct_any_temperature_field_for_visu(IJK_Field_double& temperature);
  void correct_temperature_for_visu() override;
  void clip_min_temperature_values() override;
  void clip_max_temperature_values() override;
  void compute_mean_liquid_temperature();
  void compute_overall_probes_parameters();

  void pre_initialise_thermal_subproblems_any_matrices();
  void pre_initialise_thermal_subproblems_matrices();

  void interpolate_indicator_on_probes();
  void clear_sort_problems_colliding_bubbles();
  void interpolate_project_velocities_on_probes();
  void reajust_probes_length_collisions();
  void reajust_probes_length();
  void compute_radial_subresolution_convection_diffusion_operators();
  void compute_local_substep();
  void prepare_temporal_schemes();
  void compute_source_terms_impose_subresolution_boundary_conditions();
  void compute_add_subresolution_source_terms();
  void compute_subresolution_temporal_explicit_implicit_matrices();
  void approximate_temperature_increment_material_derivative();
  void compute_radial_first_second_order_operators(Matrice& radial_first_order_operator_raw,
                                                   Matrice& radial_second_order_operator_raw,
                                                   Matrice& radial_first_order_operator,
                                                   Matrice& radial_second_order_operator);
  void compute_first_order_operator_raw(Matrice& radial_first_order_operator);
  void compute_first_order_operator(Matrice& radial_first_order_operator, double dr);
  void compute_second_order_operator(Matrice& radial_second_order_operator, double dr);
  void compute_second_order_operator_raw(Matrice& radial_second_order_operator);
  void initialise_identity_matrices(Matrice& identity_matrix,
                                    Matrice& identity_matrix_subproblems);
  void initialise_identity_matrices_sparse(Matrice& identity_matrix,
                                           Matrice& identity_matrix_subproblems);
  void initialise_radial_convection_operator(Matrice& radial_first_order_operator,
                                             Matrice& radial_convection_matrix);
  void initialise_radial_convection_operator_sparse(Matrice& radial_first_order_operator,
                                                    Matrice& radial_convection_matrix);
  void initialise_radial_diffusion_operator(Matrice& radial_second_order_operator,
                                            Matrice& radial_diffusion_matrix);
  void initialise_radial_diffusion_operator_sparse(Matrice& radial_second_order_operator,
                                                   Matrice& radial_diffusion_matrix);
  // int copy_local_unknwowns_rhs();
  void convert_into_sparse_matrix();
  void compute_md_vector();
  void retrieve_temperature_solution();
  void store_previous_temperature_indicator_velocities();
  void check_wrong_values_rhs();
  void initialise_thermal_subproblems_list();
  void initialise_thermal_subproblems();
  void detect_probe_collision();
  void solve_thermal_subproblems();
  void prepare_thermal_flux_correction();
  void compute_min_max_reachable_fluxes();
  void complete_convective_flux_frame_of_reference();
  void update_intersections() override;
  void compute_convective_diffusive_fluxes_face_centre() override;
  void compute_convective_fluxes_face_centre() override;
  void compute_diffusive_fluxes_face_centre() override;

  void complete_thermal_fluxes_face_centre(const int& fluxes_correction_conservations);

  void compute_temperature_cell_centres(const int first_corr) override;
  void compute_temperature_cell_centres_first_correction();
  void compute_temperature_cell_centres_second_correction();
  void replace_temperature_cell_centres_neighbours(const int& use_neighbours_temperature_to_correct_trimmed);
  void prepare_ij_fluxes_k_layers() override;
  void set_zero_temperature_increment() override;
  void clean_thermal_subproblems() override;
  void clean_ijk_intersections() override;
  void clean_add_thermal_subproblems();
  void enforce_periodic_temperature_boundary_value() override;
  void correct_operators_for_visu() override;

  double get_modified_time() override;
  void compute_temperature_init() override;
  void recompute_temperature_init() override;
  void set_field_temperature_per_bubble(const int index_bubble);
  Nom compute_quasi_static_spherical_diffusion_expression(const double& time_scope, const int index_bubble, const int index_bubble_real);
  Nom generate_expression_temperature_ini(const double& time_scope, const double x, const double y, const double z);
  void approx_erf_inverse(const double& x, double& res);
  void set_field_T_ana() override;
  void calculer_ecart_T_ana() override { ; };
  double compute_spherical_steady_dirichlet_left_right_value(const double& r);
  double compute_spherical_steady_dirichlet_left_right_derivative_value(const double& r, const double& temperature_prev);
  double compute_spherical_steady_dirichlet_left_right_integral(const double& temperature_end_prev);
  double find_time_dichotomy_integral(const double& temperature_integral, double& temperature_end_prev);
  void compute_Nusselt_spherical_diffusion();
  double get_time_inflection_derivative(double& temperature_end_min);
  double find_time_dichotomy_derivative(const double& temperature_derivative, double& temperature_limit_left, double& temperature_limit_right);

  /* compute_rho_cp_u_mean() May be clearly overridden later */
  double compute_rho_cp_u_mean(const IJK_Field_double& vx) override { return IJK_Thermal_base::compute_rho_cp_u_mean(vx); };

  void compare_fluxes_thermal_subproblems() override;

  int enable_probe_collision_detection_ = 0;
  int enable_resize_probe_collision_ = 0;
  int debug_probe_collision_ = 0;
  IJK_Field_double probe_collision_debug_field_;
  int reference_gfm_on_probes_ = 0;
  int compute_normal_derivatives_on_reference_probes_ = 0;

  int disable_probe_weak_gradient_ = 0;
  int disable_probe_weak_gradient_gfm_ = 0;

  int reconstruct_previous_probe_field_ = 0;
  int implicit_solver_from_previous_probe_field_ = 0;

  int disable_spherical_diffusion_start_ = 0;
  int single_centred_bubble_ = 1;
  int computed_centred_bubble_start_ = 1;
  double single_centred_bubble_radius_ini_ = 1.e-3;
  double probes_end_value_start_ = -1;
  double probes_end_value_coeff_ = 0.05;
  int temperature_ini_type_ = 1;
  double modified_time_init_ = 0.;
  int spherical_diffusion_ = 1;
  double nusselt_spherical_diffusion_ = 2.;
  double nusselt_spherical_diffusion_liquid_ = 2.;
  double heat_flux_spherical_ = 0.;
  enum temperature_ini_dict { local_criteria, integral_criteria, derivative_criteria, time_criteria };
  double mean_liquid_temperature_ = -1;
  double time_ini_user_ = 0.;

  int disable_mixed_cells_increment_ = 0;
  int enable_mixed_cells_increment_ = 1;
  int allow_temperature_correction_for_visu_ = 0;
  int disable_subresolution_ = 0;
  int diffusive_flux_correction_ = 0;
  int convective_flux_correction_ = 0;
  int fluxes_correction_conservations_ = 0;
  int conserve_max_interfacial_fluxes_ = 0;
  int fluxes_corrections_weighting_ = 0;
  int impose_fo_flux_correction_ = 1;
  int disable_fo_flux_correction_ = 0;
  int subproblem_temperature_extension_ = 0; // ghost fluid extension based on the interfacial gradient computed with the subproblem

  int override_vapour_mixed_values_ = 0; // For debug purposes

  IJK_One_Dimensional_Subproblems thermal_local_subproblems_;
  int points_per_thermal_subproblem_ = 32;
  double coeff_distance_diagonal_ = 2.;
  double probe_length_ = 0.;
  double dr_ = 0.;
  DoubleVect radial_coordinates_;
  Matrice identity_matrix_explicit_implicit_;
  Matrice radial_first_order_operator_raw_;
  Matrice radial_second_order_operator_raw_;
  Matrice radial_first_order_operator_;
  Matrice radial_second_order_operator_;
  Matrice identity_matrix_subproblems_;
  Matrice radial_diffusion_matrix_;
  Matrice radial_convection_matrix_;

  Matrice radial_diffusion_matrix_test_;

  FixedVector<ArrOfInt, 6> first_indices_sparse_matrix_;

  int initialise_sparse_indices_ = 0;
  /*
   * Thermal subproblems are regrouped in a single linear system AX=b
   * on each processor !
   */
  IJK_Finite_Difference_One_Dimensional_Matrix_Assembler finite_difference_assembler_;
  Matrice thermal_subproblems_matrix_assembly_;
  DoubleVect thermal_subproblems_temperature_solution_ini_;
  DoubleVect thermal_subproblems_rhs_assembly_;
  DoubleVect thermal_subproblems_temperature_solution_;

  /*
   * TODO: Cast the matrice with Matrice Morse directly (not Matrice Bloc)
   */
  Matrice * thermal_subproblems_matrix_assembly_for_solver_ref_ = nullptr;
  Matrice thermal_subproblems_matrix_assembly_for_solver_;
  Matrice thermal_subproblems_matrix_assembly_for_solver_reduced_;

  // SolveurSys one_dimensional_advection_diffusion_thermal_solver_;
  IJK_SolveSys_FD_thermal one_dimensional_advection_diffusion_thermal_solver_;
  IJK_SolveSys_FD_thermal one_dimensional_advection_diffusion_thermal_solver_implicit_;
  MD_Vector md_;
  Motcles fd_solvers_;
  Motcles fd_solvers_jdd_;
  int fd_solver_rank_ = 0;
  Nom fd_solver_type_;
  int discrete_integral_ = 0;
  int quadtree_levels_ = 1;

  int compute_tangential_variables_ = 0;

  int boundary_condition_interface_ = -1;
  Motcles boundary_condition_interface_dict_;
  double interfacial_boundary_condition_value_ = 0.;
  int impose_boundary_condition_interface_from_simulation_ = 0;
  int boundary_condition_end_ = 0;
  Motcles boundary_condition_end_dict_;
  double end_boundary_condition_value_ = -1.;
  int impose_user_boundary_condition_end_value_ = 0;
  int source_terms_type_ = 2;
  Motcles source_terms_type_dict_;
  int source_terms_correction_ = 0;
  int advected_frame_of_reference_ = 0;
  int neglect_frame_of_reference_radial_advection_ = 0;
  int approximate_temperature_increment_ = 0;

  /*
   * Some tries to do explicit temporal variations at the beginning
   */
  bool is_first_time_step_ = false;
  int first_time_step_temporal_ = 0;
  int first_time_step_explicit_ = 0;
  int first_time_step_implicit_ = 0;
  int local_diffusion_fourier_priority_ = 0;
  int nb_iter_explicit_local_ = 0;
  double local_fourier_ = 1.;
  double local_cfl_ = 1.;
  double delta_T_subcooled_overheated_=-1.;
  double convection_diffusion_time_scale_factor_ = 1.;
  double local_fourier_time_step_probe_length_ = 0.;
  double local_cfl_time_step_probe_length_ = 0.;
  double local_dt_cfl_ = 0.;
  double local_dt_cfl_min_delta_xyz_ = 0.;
  double local_dt_cfl_counter_ = 0.;
  double local_dt_fourier_counter_ = 0.;

  double error_temperature_ana_total_ = 0.;
  double error_temperature_ana_squared_total_ = 0.;
  double error_temperature_ana_rel_total_ = 0.;

  /*
   * Some tries to make the probe length varies at the beginning of the simulation
   */
  int first_time_step_varying_probes_ = 0;
  int probe_variations_enabled_ = 0;
  int probe_variations_priority_ = 0;
  int disable_interpolation_in_mixed_cells_ = 0;
  int keep_temperature_extrapolated_from_LRS_ = 0;
  int max_u_radial_=0;

  IJK_Field_double debug_LRS_cells_;
  int distance_cell_faces_from_lrs_ = 1;
  int  disable_distance_cell_faces_from_lrs_ = 0;

  int pre_initialise_thermal_subproblems_list_ = 0;
  double pre_factor_subproblems_number_ = 3.;
  int remove_append_subproblems_ = 0;
  int use_sparse_matrix_ = 0;
  int global_probes_characteristics_ = 1;

  int correct_temperature_cell_neighbours_first_iter_ = 0;
  int correct_first_iter_deactivate_cell_neighbours_ = 0;
  int find_temperature_cell_neighbours_ = 0;
  int use_temperature_cell_neighbours_ = 0;
  int correct_neighbours_using_probe_length_ = 0;
  int neighbours_corrected_rank_ = 1;
  int neighbours_weighting_ = 0;
  int neighbours_colinearity_weighting_ = 0;
  int neighbours_distance_weighting_ = 0;
  int neighbours_colinearity_distance_weighting_ = 0;
  int smooth_temperature_field_ = 0;
  int readjust_probe_length_from_vertices_ = 0;
  IJK_Field_double temperature_cell_neighbours_;
  IJK_Field_double temperature_cell_neighbours_debug_;
  IJK_Field_int neighbours_temperature_to_correct_;
  IJK_Field_double neighbours_temperature_colinearity_weighting_;

  int keep_max_flux_correction_ = 0;

  int clip_temperature_values_ = 0;
  int enforce_periodic_boundary_value_ = 0;
  int stencil_periodic_boundary_value_ = 2;

  int disable_post_processing_probes_out_files_ = 0;

  /*
   * Pure cells corrected for visualisation
   */
  IJK_Field_vector3_double cell_faces_corrected_diffusive_;
  IJK_Field_vector3_double cell_faces_corrected_convective_;
  IJK_Field_vector3_int cell_faces_corrected_bool_;
  int store_cell_faces_corrected_ = 0;

  /*
   * Neighbouring faces in the diagonal
   */
  IJK_Field_vector3_int cell_faces_neighbours_corrected_diag_bool_;
  int find_cell_neighbours_for_fluxes_spherical_correction_ = 0;
  int use_cell_neighbours_for_fluxes_spherical_correction_ = 0;

  /*
   * All reachable faces to correct
   */
  int find_reachable_fluxes_ = 0;
  int use_reachable_fluxes_ = 0;
  int keep_first_reachable_fluxes_ = 0;
  IJK_Field_vector3_int cell_faces_neighbours_corrected_all_bool_;
  IJK_Field_vector3_double cell_faces_neighbours_corrected_velocity_temperature_;
  IJK_Field_vector3_double cell_faces_neighbours_corrected_convective_frame_of_reference_;
  IJK_Field_vector3_double cell_faces_neighbours_corrected_convective_;
  IJK_Field_vector3_double cell_faces_neighbours_corrected_diffusive_;
  IJK_Field_vector3_double neighbours_faces_weighting_colinearity_;
  IJK_Field_vector3_int cell_faces_neighbours_corrected_min_max_bool_;
  IJK_Field_int neighbours_temperature_to_correct_trimmed_;
  int neighbours_last_faces_weighting_ = 0;
  int neighbours_last_faces_colinearity_weighting_ = 0;
  int neighbours_last_faces_colinearity_face_weighting_ = 0;
  int neighbours_last_faces_distance_weighting_ = 0;
  int neighbours_last_faces_distance_colinearity_weighting_ = 0;
  int neighbours_last_faces_distance_colinearity_face_weighting_ = 0;

  int post_process_all_probes_ = 0;
  int nb_theta_post_pro_ = 10;
  int nb_phi_post_pro_ = 4;
  int nb_probes_post_pro_ = 40;

  int interp_eulerian_ = 0;
  int first_step_thermals_post_ = 1;
  int disable_first_step_thermals_post_ = 0;

  int copy_fluxes_on_every_procs_ = 1;
  int copy_temperature_on_every_procs_ = 1;

  int post_process_thermal_slices_ = 0;
  int thermal_slices_regions_ = 0;
  double nb_diam_slice_ = -1;
  int nb_slices_ = 1;
  int upstream_dir_slice_ = -1;
  int disable_slice_to_nearest_plane_ = 0;

  int post_process_thermal_lines_ = 0;
  int upstream_dir_line_ = -1;
  int nb_thermal_lines_ = 1;
  int nb_thermal_concentric_circles_ = 1;
  int nb_thermal_line_points_ = 100;
  double nb_diam_thermal_line_length_ = -1;

  int use_corrected_velocity_convection_ = 0;
  int use_velocity_cartesian_grid_ = 0;
  int compute_radial_displacement_ = 0;
  int use_normal_gradient_for_flux_corr_ = 0;

  IJK_Field_vector3_double interfacial_heat_flux_dispatched_;
  FixedVector<ArrOfInt, 3> ijk_indices_flux_out_;
  FixedVector<ArrOfDouble, 3> thermal_flux_out_;

  IJK_Field_int zero_liquid_neighbours_;
  IJK_Field_vector3_double interfacial_heat_flux_contrib_;
  IJK_Field_vector3_double interfacial_heat_flux_current_;
  FixedVector<ArrOfInt, 4> ijk_indices_flux_contrib_;
  ArrOfDouble thermal_flux_out_contrib_;
};

#endif /* IJK_Thermal_Subresolution_included */
