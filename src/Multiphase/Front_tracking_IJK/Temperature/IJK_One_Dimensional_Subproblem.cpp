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
// File      : IJK_One_Dimensional_Subproblem.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_One_Dimensional_Subproblem.h>
#include <IJK_Navier_Stokes_tools.h>
#include <Ouvrir_fichier.h>
#include <IJK_FT.h>
#include <IJK_Thermal_base.h>
#include <IJK_Thermal_Subresolution.h>
#include <IJK_One_Dimensional_Subproblems.h>


Implemente_instanciable_sans_constructeur( IJK_One_Dimensional_Subproblem, "IJK_One_Dimensional_Subproblem", Objet_U ) ;

IJK_One_Dimensional_Subproblem::IJK_One_Dimensional_Subproblem()
{

}

IJK_One_Dimensional_Subproblem::IJK_One_Dimensional_Subproblem(const IJK_FT_double& ijk_ft) : IJK_One_Dimensional_Subproblem()
{
  associer(ijk_ft);
}

Sortie& IJK_One_Dimensional_Subproblem::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  return os;
}

Entree& IJK_One_Dimensional_Subproblem::readOn( Entree& is )
{
  Objet_U::readOn( is );
  return is;
}

static int compute_periodic_index(const int index, const int n)
{
  return (n + index % n) % n;
}

void IJK_One_Dimensional_Subproblem::get_ijk_indices(int& i, int& j, int& k) const
{
  i = (int) index_i_;
  j = (int) index_j_;
  k = (int) index_k_;
}

void IJK_One_Dimensional_Subproblem::reinit_variable(DoubleVect& vect)
{
  vect.resize(*points_per_thermal_subproblem_);
  vect *= 0.;
}

/*
 * TODO: Remplacer cette methode pour eviter de fournir tous les attributs
 */
void IJK_One_Dimensional_Subproblem::associate_sub_problem_to_inputs(IJK_Thermal_Subresolution& ref_thermal_subresolution,
                                                                     IJK_One_Dimensional_Subproblems& ref_one_dimensional_subproblems,
                                                                     int i, int j, int k,
                                                                     int init,
                                                                     int sub_problem_index,
                                                                     double global_time_step,
                                                                     double current_time,
                                                                     int compo_connex,
                                                                     double distance,
                                                                     double curvature,
                                                                     double interfacial_area,
                                                                     ArrOfDouble facet_barycentre,
                                                                     ArrOfDouble normal_vector,
                                                                     double bubble_rising_velocity,
                                                                     ArrOfDouble bubble_rising_vector,
                                                                     ArrOfDouble bubble_barycentre,
                                                                     const double& indicator,
                                                                     const IJK_Interfaces& interfaces,
                                                                     const FixedVector<IJK_Field_double, 3>& velocity,
                                                                     const FixedVector<IJK_Field_double, 3>& velocity_ft,
                                                                     const IJK_Field_double& pressure)
{
  /*
   * Should not change over iterations
   */
  init_ = init;
  sub_problem_index_ = sub_problem_index;

  if (init_)
    {
      associate_thermal_subproblem_parameters(ref_thermal_subresolution.reference_gfm_on_probes_,
                                              ref_thermal_subresolution.debug_,
                                              ref_thermal_subresolution.n_iter_distance_,
                                              ref_thermal_subresolution.delta_T_subcooled_overheated_,
                                              ref_thermal_subresolution.pre_initialise_thermal_subproblems_list_,
                                              ref_thermal_subresolution.use_sparse_matrix_,
                                              ref_thermal_subresolution.compute_normal_derivatives_on_reference_probes_,
                                              ref_thermal_subresolution.latastep_reprise_ini_);
      associate_thermal_subproblem_sparse_matrix(ref_thermal_subresolution.first_indices_sparse_matrix_);
      associate_eulerian_fields_references(interfaces,
                                           ref_thermal_subresolution.eulerian_distance_ns_,
                                           ref_thermal_subresolution.eulerian_curvature_ns_,
                                           ref_thermal_subresolution.eulerian_interfacial_area_ns_,
                                           ref_thermal_subresolution.eulerian_normal_vectors_ns_,
                                           ref_thermal_subresolution.eulerian_facets_barycentre_ns_,
                                           ref_thermal_subresolution.temperature_,
                                           ref_thermal_subresolution.temperature_ft_,
                                           ref_thermal_subresolution.temperature_before_extrapolation_,
                                           velocity,
                                           velocity_ft,
                                           pressure,
                                           ref_thermal_subresolution.grad_T_elem_,
                                           ref_thermal_subresolution.grad_T_elem_smooth_,
                                           ref_thermal_subresolution.hess_diag_T_elem_,
                                           ref_thermal_subresolution.hess_cross_T_elem_,
                                           ref_thermal_subresolution.eulerian_grad_T_interface_ns_,
                                           ref_thermal_subresolution.probe_collision_debug_field_,
                                           ref_thermal_subresolution.zero_liquid_neighbours_,
                                           ref_thermal_subresolution.smooth_grad_T_elem_);
      associate_probe_parameters(ref_thermal_subresolution.points_per_thermal_subproblem_,
                                 ref_thermal_subresolution.cp_liquid_,
                                 ref_thermal_subresolution.uniform_alpha_,
                                 ref_thermal_subresolution.uniform_lambda_,
                                 ref_thermal_subresolution.prandtl_number_,
                                 ref_thermal_subresolution.coeff_distance_diagonal_,
                                 ref_thermal_subresolution.cell_diagonal_,
                                 ref_thermal_subresolution.dr_,
                                 ref_thermal_subresolution.radial_coordinates_);
      associate_bubble_parameters(ref_one_dimensional_subproblems.total_surface_per_bubble_,
                                  ref_one_dimensional_subproblems.radius_from_surfaces_per_bubble_,
                                  ref_one_dimensional_subproblems.radius_from_volumes_per_bubble_,
                                  ref_one_dimensional_subproblems.delta_temperature_,
                                  ref_thermal_subresolution.mean_liquid_temperature_,
                                  ref_thermal_subresolution.bubbles_volume_,
                                  ref_thermal_subresolution.rising_vectors_);
      associate_global_subproblems_parameters(ref_thermal_subresolution.reconstruct_previous_probe_field_,
                                              ref_thermal_subresolution.implicit_solver_from_previous_probe_field_,
                                              ref_one_dimensional_subproblems.subproblem_to_ijk_indices_previous_,
                                              ref_one_dimensional_subproblems.temperature_probes_previous_,
                                              ref_one_dimensional_subproblems.indicator_probes_previous_,
                                              ref_one_dimensional_subproblems.velocities_probes_previous_,
                                              ref_one_dimensional_subproblems.normal_vector_compo_probes_previous_);
      associate_finite_difference_operators(ref_thermal_subresolution.radial_first_order_operator_raw_,
                                            ref_thermal_subresolution.radial_second_order_operator_raw_,
                                            ref_thermal_subresolution.radial_first_order_operator_,
                                            ref_thermal_subresolution.radial_second_order_operator_,
                                            ref_thermal_subresolution.identity_matrix_explicit_implicit_,
                                            ref_thermal_subresolution.identity_matrix_subproblems_,
                                            ref_thermal_subresolution.radial_diffusion_matrix_,
                                            ref_thermal_subresolution.radial_convection_matrix_);
      associate_finite_difference_solver_solution(ref_thermal_subresolution.finite_difference_assembler_,
                                                  ref_thermal_subresolution.thermal_subproblems_matrix_assembly_,
                                                  ref_thermal_subresolution.thermal_subproblems_rhs_assembly_,
                                                  ref_thermal_subresolution.thermal_subproblems_temperature_solution_,
                                                  ref_thermal_subresolution.thermal_subproblems_temperature_solution_ini_);
      associate_source_terms_parameters(ref_thermal_subresolution.source_terms_type_,
                                        ref_thermal_subresolution.source_terms_correction_,
                                        ref_thermal_subresolution.source_terms_correction_,
                                        ref_thermal_subresolution.advected_frame_of_reference_,
                                        ref_thermal_subresolution.neglect_frame_of_reference_radial_advection_,
                                        ref_thermal_subresolution.compute_tangential_variables_);
      associate_flux_correction_parameters((ref_thermal_subresolution.convective_flux_correction_
                                            || ref_thermal_subresolution.diffusive_flux_correction_),
                                           ref_thermal_subresolution.distance_cell_faces_from_lrs_,
                                           ref_thermal_subresolution.interp_eulerian_,
                                           ref_thermal_subresolution.use_corrected_velocity_convection_,
                                           ref_thermal_subresolution.use_velocity_cartesian_grid_,
                                           ref_thermal_subresolution.compute_radial_displacement_,
                                           ref_thermal_subresolution.fluxes_correction_conservations_,
                                           ref_thermal_subresolution.fluxes_corrections_weighting_,
                                           ref_thermal_subresolution.use_normal_gradient_for_flux_corr_);
      associate_varying_probes_params(ref_thermal_subresolution.readjust_probe_length_from_vertices_,
                                      ref_thermal_subresolution.first_time_step_varying_probes_,
                                      ref_thermal_subresolution.probe_variations_priority_,
                                      ref_thermal_subresolution.disable_interpolation_in_mixed_cells_);
      associate_flags_neighbours_correction(ref_thermal_subresolution.find_temperature_cell_neighbours_,
                                            !ref_thermal_subresolution.correct_neighbours_using_probe_length_,
                                            ref_thermal_subresolution.neighbours_corrected_rank_,
                                            ref_thermal_subresolution.neighbours_colinearity_weighting_,
                                            ref_thermal_subresolution.neighbours_distance_weighting_,
                                            ref_thermal_subresolution.neighbours_colinearity_distance_weighting_,
                                            ref_thermal_subresolution.neighbours_last_faces_colinearity_weighting_,
                                            ref_thermal_subresolution.neighbours_last_faces_colinearity_face_weighting_,
                                            ref_thermal_subresolution.neighbours_last_faces_distance_weighting_,
                                            ref_thermal_subresolution.neighbours_last_faces_distance_colinearity_weighting_,
                                            ref_thermal_subresolution.neighbours_last_faces_distance_colinearity_face_weighting_,
                                            ref_thermal_subresolution.find_reachable_fluxes_,
                                            ref_thermal_subresolution.find_cell_neighbours_for_fluxes_spherical_correction_);
      associate_tweaked_parameters(ref_thermal_subresolution.disable_probe_weak_gradient_,
                                   ref_thermal_subresolution.disable_probe_weak_gradient_gfm_);
      associate_collisions_parameters(ref_thermal_subresolution.enable_probe_collision_detection_,
                                      ref_thermal_subresolution.enable_resize_probe_collision_,
                                      ref_thermal_subresolution.debug_probe_collision_);

//      // TODO: If the calculation of the distance is changed in intersection_ijk, it will be useless...
      associate_sub_problem_temporal_params(ref_thermal_subresolution.is_first_time_step_,
                                            ref_thermal_subresolution.first_time_step_temporal_,
                                            ref_thermal_subresolution.first_time_step_explicit_,
                                            ref_thermal_subresolution.local_fourier_,
                                            ref_thermal_subresolution.local_cfl_,
                                            ref_thermal_subresolution.min_delta_xyz_,
                                            ref_thermal_subresolution.max_u_radial_);
      disable_relative_velocity_energy_balance_ = ref_thermal_subresolution.disable_relative_velocity_energy_balance_;
    }

  /*
   * Should be reinitialised at each time step.
   */
  clear_vectors();
  reset_counters();
  reset_flags();
  set_global_index(0);
  reset_post_processing_theta_phi_scope();
  associate_temporal_parameters(global_time_step, current_time);
  associate_cell_ijk(i, j, k);
  associate_eulerian_field_values(compo_connex, indicator);
  associate_interface_related_parameters(distance, curvature, interfacial_area, facet_barycentre, normal_vector);
  associate_rising_velocity(bubble_rising_velocity, bubble_rising_vector, bubble_barycentre);
  initialise_thermal_probe();
  if (!global_probes_characteristics_)
    (*first_time_step_temporal_) = 0;
  recompute_finite_difference_matrices();
}

void IJK_One_Dimensional_Subproblem::clear_vectors()
{
  pure_neighbours_to_correct_.clear();
  pure_neighbours_corrected_distance_.clear();
  pure_neighbours_corrected_colinearity_.clear();
  pure_neighbours_last_faces_to_correct_.clear();
  pure_neighbours_last_faces_corrected_distance_.clear();
  pure_neighbours_last_faces_corrected_colinearity_.clear();
  for (int l=0; l<6; l++)
    {
      convective_flux_op_lrs_[l] = 0.;
      diffusive_flux_op_lrs_[l] = 0.;
      corrective_flux_current_[l] = 0.;
      corrective_flux_to_neighbours_[l] = 0.;
      corrective_flux_from_neighbours_[l] = 0.;
      temperature_interp_conv_flux_[l] = 0.;
    }
  sum_convective_flux_op_lrs_ = 0.;
  sum_convective_flux_op_leaving_lrs_ = 0.;
  sum_convective_flux_op_entering_lrs_ = 0.;
  sum_diffusive_flux_op_lrs_ = 0.;
  sum_diffusive_flux_op_leaving_lrs_ = 0.;
  sum_diffusive_flux_op_entering_lrs_ = 0.;
}

void IJK_One_Dimensional_Subproblem::reset_counters()
{
  velocities_calculation_counter_ = 0;
}

void IJK_One_Dimensional_Subproblem::reset_flags()
{
  has_computed_cell_centre_distance_ = false;
  has_computed_cell_faces_distance_ = false;
  has_computed_liquid_neighbours_ = false;
  has_computed_lrs_flux_frame_of_ref_terms_ = false;
  disable_probe_because_collision_ = 0;
  disable_probe_weak_gradient_local_= 0;
}

void IJK_One_Dimensional_Subproblem::associate_thermal_subproblem_parameters(const int& reference_gfm_on_probes,
                                                                             const int& debug,
                                                                             const int& n_iter_distance,
                                                                             const double& delta_T_subcooled_overheated,
                                                                             const int& pre_initialise_thermal_subproblems_list,
                                                                             const int& use_sparse_matrix,
                                                                             const int& compute_normal_derivative_on_reference_probes,
                                                                             const int& latastep_reprise)
{
  reference_gfm_on_probes_ = reference_gfm_on_probes;
  debug_ = debug;
  n_iter_distance_ = n_iter_distance;
  delta_T_subcooled_overheated_ = delta_T_subcooled_overheated;
  pre_initialise_thermal_subproblems_list_ = pre_initialise_thermal_subproblems_list;
  use_sparse_matrix_ = use_sparse_matrix;
  compute_normal_derivative_on_reference_probes_ = compute_normal_derivative_on_reference_probes;
  latastep_reprise_ = &latastep_reprise;
}

void IJK_One_Dimensional_Subproblem::associate_thermal_subproblem_sparse_matrix(FixedVector<ArrOfInt,6>& first_indices_sparse_matrix)
{
  first_indices_sparse_matrix_ = &first_indices_sparse_matrix;
}

void IJK_One_Dimensional_Subproblem::associate_eulerian_fields_references(const IJK_Interfaces& interfaces,
                                                                          const IJK_Field_double * eulerian_distance,
                                                                          const IJK_Field_double * eulerian_curvature,
                                                                          const IJK_Field_double * eulerian_interfacial_area,
                                                                          const FixedVector<IJK_Field_double, 3> * eulerian_normal_vect,
                                                                          const FixedVector<IJK_Field_double, 3> * eulerian_facets_barycentre,
                                                                          const IJK_Field_double& temperature,
                                                                          const IJK_Field_double& temperature_ft,
                                                                          const IJK_Field_double& temperature_before_extrapolation,
                                                                          const FixedVector<IJK_Field_double, 3>& velocity,
                                                                          const FixedVector<IJK_Field_double, 3>& velocity_ft,
                                                                          const IJK_Field_double& pressure,
                                                                          const FixedVector<IJK_Field_double, 3>& grad_T_elem,
                                                                          const FixedVector<IJK_Field_double, 3>& grad_T_elem_smooth,
                                                                          const FixedVector<IJK_Field_double, 3>& hess_diag_T_elem,
                                                                          const FixedVector<IJK_Field_double, 3>& hess_cross_T_elem,
                                                                          const IJK_Field_double& eulerian_grad_T_interface_ns,
                                                                          IJK_Field_double& probe_collision_debug_field,
                                                                          IJK_Field_int& zero_liquid_neighbours,
                                                                          const int& smooth_grad_T_elem)
{
  interfaces_ = &interfaces;
  eulerian_distance_ = eulerian_distance;
  eulerian_curvature_ = eulerian_curvature;
  eulerian_interfacial_area_ = eulerian_interfacial_area;
  eulerian_normal_vect_ = eulerian_normal_vect;
  eulerian_facets_barycentre_ = eulerian_facets_barycentre;
  temperature_ = &temperature ;
  temperature_ft_ = &temperature_ft ;
  temperature_before_extrapolation_ = &temperature_before_extrapolation;
  velocity_ = &velocity ;
  velocity_ft_ = &velocity_ft;
  pressure_ = &pressure;
  grad_T_elem_ = &grad_T_elem;
  grad_T_elem_smooth_ = &grad_T_elem_smooth;
  smooth_grad_T_elem_ = smooth_grad_T_elem;
  if (smooth_grad_T_elem_)
    grad_T_elem_solver_ = &grad_T_elem_smooth;
  else
    grad_T_elem_solver_ = &grad_T_elem;
  hess_diag_T_elem_ = &hess_diag_T_elem;
  hess_cross_T_elem_ = &hess_cross_T_elem;
  eulerian_grad_T_interface_ns_ = &eulerian_grad_T_interface_ns;
  probe_collision_debug_field_ = &probe_collision_debug_field;
  zero_liquid_neighbours_ = &zero_liquid_neighbours;
}

void IJK_One_Dimensional_Subproblem::associate_tweaked_parameters(const int& disable_probe_weak_gradient,
                                                                  const int& disable_probe_weak_gradient_gfm)
{
  disable_probe_weak_gradient_ = disable_probe_weak_gradient;
  disable_probe_weak_gradient_gfm_ = disable_probe_weak_gradient_gfm;
  if (disable_probe_weak_gradient_gfm_)
    disable_probe_weak_gradient_ = 1;
}

void IJK_One_Dimensional_Subproblem::associate_collisions_parameters(const int& enable_probe_collision_detection,
                                                                     const int& enable_resize_probe_collision,
                                                                     const int& debug_probe_collision)
{
  enable_probe_collision_detection_ = enable_probe_collision_detection;
  enable_resize_probe_collision_ = enable_resize_probe_collision;
  debug_probe_collision_ = debug_probe_collision;
}

void IJK_One_Dimensional_Subproblem::associate_sub_problem_temporal_params(const bool& is_first_time_step,
                                                                           int& first_time_step_temporal,
                                                                           const int& first_time_step_explicit,
                                                                           const double& local_fourier,
                                                                           const double& local_cfl,
                                                                           const double& min_delta_xyz,
                                                                           int max_u_radial)
{
  is_first_time_step_ = is_first_time_step;
  first_time_step_temporal_ = &first_time_step_temporal;
  first_time_step_explicit_ = first_time_step_explicit;
  local_fourier_ = local_fourier;
  local_cfl_ = local_cfl;
  min_delta_xyz_ = min_delta_xyz;
  max_u_radial_ = max_u_radial;
  max_u_cartesian_ = !max_u_radial_;
}

void IJK_One_Dimensional_Subproblem::associate_varying_probes_params(const int& readjust_probe_length_from_vertices,
                                                                     const int& first_time_step_varying_probes,
                                                                     const int& probe_variations_priority,
                                                                     const int& disable_interpolation_in_mixed_cells)
{
  readjust_probe_length_from_vertices_ = readjust_probe_length_from_vertices;
  first_time_step_varying_probes_ = first_time_step_varying_probes;
  probe_variations_priority_ = probe_variations_priority;
  disable_interpolation_in_mixed_cells_ = disable_interpolation_in_mixed_cells;
}

void IJK_One_Dimensional_Subproblem::associate_flux_correction_parameters(const int& correct_fluxes,
                                                                          const int& distance_cell_faces_from_lrs,
                                                                          const int& interp_eulerian,
                                                                          const int& use_corrected_velocity_convection,
                                                                          const int& use_velocity_cartesian_grid,
                                                                          const int& compute_radial_displacement,
                                                                          const int& fluxes_correction_conservations,
                                                                          const int& fluxes_corrections_weighting,
                                                                          const int& use_normal_gradient_for_flux_corr)
{
  correct_fluxes_ = correct_fluxes;
  distance_cell_faces_from_lrs_ = distance_cell_faces_from_lrs;
  interp_eulerian_ = interp_eulerian;
  use_corrected_velocity_convection_ = use_corrected_velocity_convection;
  use_velocity_cartesian_grid_ = use_velocity_cartesian_grid;
  compute_radial_displacement_ = compute_radial_displacement;
  fluxes_correction_conservations_ = fluxes_correction_conservations;
  fluxes_corrections_weighting_ = fluxes_corrections_weighting;
  use_normal_gradient_for_flux_corr_ = use_normal_gradient_for_flux_corr;
}

void IJK_One_Dimensional_Subproblem::associate_source_terms_parameters(const int& source_terms_type,
                                                                       const int& correct_tangential_temperature_gradient,
                                                                       const int& correct_tangential_temperature_hessian,
                                                                       const int& advected_frame_of_reference,
                                                                       const int& neglect_frame_of_reference_radial_advection,
                                                                       const int& compute_tangential_variables)
{
  source_terms_type_ = source_terms_type;
  correct_tangential_temperature_gradient_ = correct_tangential_temperature_gradient;
  correct_tangential_temperature_hessian_ = correct_tangential_temperature_hessian;
  advected_frame_of_reference_ = advected_frame_of_reference;
  neglect_frame_of_reference_radial_advection_ = neglect_frame_of_reference_radial_advection;
  pure_thermal_diffusion_ = (source_terms_type_ == linear_diffusion
                             || source_terms_type_ == spherical_diffusion
                             || source_terms_type_ == spherical_diffusion_approx);
  if (pure_thermal_diffusion_)
    compute_tangential_variables_ = compute_tangential_variables;
}

void 	IJK_One_Dimensional_Subproblem::associate_finite_difference_solver_solution(IJK_Finite_Difference_One_Dimensional_Matrix_Assembler& finite_difference_assembler,
                                                                                  Matrice& thermal_subproblems_matrix_assembly,
                                                                                  DoubleVect& thermal_subproblems_rhs_assembly,
                                                                                  DoubleVect& thermal_subproblems_temperature_solution,
                                                                                  DoubleVect& thermal_subproblems_temperature_solution_ini)
{
  finite_difference_assembler_ = &finite_difference_assembler;
  thermal_subproblems_matrix_assembly_ = &thermal_subproblems_matrix_assembly;
  thermal_subproblems_rhs_assembly_ = &thermal_subproblems_rhs_assembly;
  thermal_subproblems_temperature_solution_ = &thermal_subproblems_temperature_solution;
  thermal_subproblems_temperature_solution_ini_ = &thermal_subproblems_temperature_solution_ini;
}

void IJK_One_Dimensional_Subproblem::associate_temporal_parameters(const double& global_time_step, const double& current_time)
{
  global_time_step_ = global_time_step;
  current_time_ = current_time;
}

void IJK_One_Dimensional_Subproblem::associate_flags_neighbours_correction(const int& correct_temperature_cell_neighbours,
                                                                           const int& correct_neighbours_rank,
                                                                           const int& neighbours_corrected_rank,
                                                                           const int& neighbours_colinearity_weighting,
                                                                           const int& neighbours_distance_weighting,
                                                                           const int& neighbours_colinearity_distance_weighting,
                                                                           const int& neighbours_last_faces_colinearity_weighting,
                                                                           const int& neighbours_last_faces_colinearity_face_weighting,
                                                                           const int& neighbours_last_faces_distance_weighting,
                                                                           const int& neighbours_last_faces_distance_colinearity_weighting,
                                                                           const int& neighbours_last_faces_distance_colinearity_face_weighting,
                                                                           const int& compute_reachable_fluxes,
                                                                           const int& find_cell_neighbours_for_fluxes_spherical_correction)
{
  /*
   * Keep it under IJK_One_Dimensional_Subproblem::associate_flux_correction_parameters()
   */
  correct_temperature_cell_neighbours_ = correct_temperature_cell_neighbours;
  correct_neighbours_rank_ = correct_neighbours_rank;
  neighbours_corrected_rank_ = neighbours_corrected_rank;
  neighbours_colinearity_weighting_ = neighbours_colinearity_weighting;
  neighbours_distance_weighting_ = neighbours_distance_weighting;
  neighbours_colinearity_distance_weighting_ = neighbours_colinearity_distance_weighting;
  neighbours_weighting_ = (neighbours_colinearity_weighting_
                           || neighbours_distance_weighting_
                           || neighbours_colinearity_distance_weighting_);
  neighbours_last_faces_colinearity_weighting_ = neighbours_last_faces_colinearity_weighting;
  neighbours_last_faces_colinearity_face_weighting_ = neighbours_last_faces_colinearity_face_weighting;
  neighbours_last_faces_distance_weighting_ = neighbours_last_faces_distance_weighting;
  neighbours_last_faces_distance_colinearity_weighting_ = neighbours_last_faces_distance_colinearity_weighting;
  neighbours_last_faces_distance_colinearity_face_weighting_ = neighbours_last_faces_distance_colinearity_face_weighting;
  neighbours_last_faces_weighting_ = (neighbours_last_faces_colinearity_weighting_
                                      || neighbours_last_faces_colinearity_face_weighting_
                                      || neighbours_last_faces_distance_weighting_
                                      || neighbours_last_faces_distance_colinearity_weighting_
                                      || neighbours_last_faces_distance_colinearity_face_weighting_);
  compute_reachable_fluxes_ = compute_reachable_fluxes;
  find_cell_neighbours_for_fluxes_spherical_correction_ = find_cell_neighbours_for_fluxes_spherical_correction;
  correct_temperature_cell_neighbours_ = (correct_temperature_cell_neighbours_ && distance_cell_faces_from_lrs_);
}

void IJK_One_Dimensional_Subproblem::associate_probe_parameters(const int& points_per_thermal_subproblem,
                                                                const double& cp_liquid,
                                                                const double& alpha,
                                                                const double& lambda,
                                                                const double& prandtl_number,
                                                                const double& coeff_distance_diagonal,
                                                                const double& cell_diagonal,
                                                                const double& dr_base,
                                                                const DoubleVect& radial_coordinates)
{
  /*
   * If the probe characteristics are the same (don't copy attributes)
   */
  points_per_thermal_subproblem_base_ = &points_per_thermal_subproblem;
  coeff_distance_diagonal_ = &coeff_distance_diagonal;
  cp_liquid_ = &cp_liquid;
  alpha_ = &alpha;
  lambda_ = &lambda;
  prandtl_number_ = &prandtl_number;
  cell_diagonal_ = &cell_diagonal;
  dr_base_ = &dr_base;
  radial_coordinates_base_ = &radial_coordinates;

  /*
   * If each probe differs (create attributes !)
   */
  if (global_probes_characteristics_)
    points_per_thermal_subproblem_ = points_per_thermal_subproblem_base_;
  else
    points_per_thermal_subproblem_ = increase_number_of_points(); //copy if modified later
}

void IJK_One_Dimensional_Subproblem::associate_bubble_parameters(const ArrOfDouble& bubbles_surface,
                                                                 const ArrOfDouble& radius_from_surfaces_per_bubble,
                                                                 const ArrOfDouble& radius_from_volumes_per_bubble,
                                                                 const double& delta_temperature,
                                                                 const double& mean_liquid_temperature,
                                                                 const ArrOfDouble * bubbles_volume,
                                                                 const DoubleTab * rising_vectors)
{
  bubbles_surface_ = &bubbles_surface;
  radius_from_surfaces_per_bubble_ = &radius_from_surfaces_per_bubble;
  radius_from_volumes_per_bubble_ = &radius_from_volumes_per_bubble;
  delta_temperature_ = &delta_temperature;
  mean_liquid_temperature_ = &mean_liquid_temperature;
  bubbles_volume_ = bubbles_volume;
  bubbles_rising_vectors_per_bubble_ = rising_vectors;
}

void  IJK_One_Dimensional_Subproblem::associate_global_subproblems_parameters(const int& reconstruct_previous_probe_field,
                                                                              const int& implicit_solver_from_previous_probe_field,
                                                                              const std::map<int, std::map<int, std::map<int, int>>>& subproblem_to_ijk_indices_previous,
                                                                              const std::vector<DoubleVect>& temperature_probe_previous,
                                                                              const std::vector<double>& indicator_probes_previous,
                                                                              const std::vector<Vecteur3>& velocities_probes_previous,
                                                                              const std::vector<Vecteur3>& normal_vector_compo_probes_previous)

{
  reconstruct_previous_probe_field_ = reconstruct_previous_probe_field;
  implicit_solver_from_previous_probe_field_ = implicit_solver_from_previous_probe_field;
  subproblem_to_ijk_indices_previous_ = &subproblem_to_ijk_indices_previous;
  temperature_probes_previous_ = &temperature_probe_previous;
  indicator_probes_previous_ = &indicator_probes_previous;
  velocities_probes_previous_ = &velocities_probes_previous;
  normal_vector_compo_probes_previous_ = &normal_vector_compo_probes_previous;
}

void IJK_One_Dimensional_Subproblem::associate_finite_difference_operators(const Matrice& radial_first_order_operator_raw,
                                                                           const Matrice& radial_second_order_operator_raw,
                                                                           const Matrice& radial_first_order_operator,
                                                                           const Matrice& radial_second_order_operator,
                                                                           const Matrice& identity_matrix_explicit_implicit_raw,
                                                                           Matrice& identity_matrix_subproblems,
                                                                           Matrice& radial_diffusion_matrix,
                                                                           Matrice& radial_convection_matrix)
{
  radial_first_order_operator_raw_base_ = &radial_first_order_operator_raw;
  radial_second_order_operator_raw_base_ = &radial_second_order_operator_raw;
  radial_first_order_operator_base_ = &radial_first_order_operator;
  radial_second_order_operator_base_ = &radial_second_order_operator;
  identity_matrix_explicit_implicit_base_ = &identity_matrix_explicit_implicit_raw;

  radial_first_order_operator_ = radial_first_order_operator_base_;
  radial_second_order_operator_ = radial_second_order_operator_base_;
  identity_matrix_explicit_implicit_ = identity_matrix_explicit_implicit_base_;

  identity_matrix_subproblems_ = &identity_matrix_subproblems;
  radial_diffusion_matrix_base_ = &radial_diffusion_matrix;
  radial_convection_matrix_base_ = &radial_convection_matrix;
}

void IJK_One_Dimensional_Subproblem::initialise_thermal_probe()
{
  if (debug_)
    Cerr << "Compute interface basis vectors" << finl;
  compute_interface_basis_vectors();

  if (debug_)
    Cerr << "Compute pure spherical basis vectors" << finl;
  compute_pure_spherical_basis_vectors();

  if (debug_)
    Cerr << "Compute probe parameters" << finl;

  probe_length_ = (*coeff_distance_diagonal_) * (*cell_diagonal_);

  /*
   *  Curvature is negative for a convex bubble
   *  but R should be positive in that case
   *  FIXME: What happen with highly deformable bubbles (concave interface portions) ?
   */
  if (fabs(curvature_) > DMINFLOAT)
    {
      osculating_radius_  = 2 / curvature_;
      if (curvature_ >= DMINFLOAT)
        {
          const double max_length = osculating_radius_ - probe_length_;
          if (curvature_ >= DMINFLOAT && max_length < 0)
            osculating_radius_ = 1.e30;
        }
      else
        osculating_radius_ = fabs(osculating_radius_);
    }

  if (debug_)
    Cerr << "Compute local discretisation" << finl;
  compute_local_discretisation();

  if (!reference_gfm_on_probes_)
    if (distance_cell_faces_from_lrs_)
      {
        compute_distance_cell_centre();
        if (debug_)
          Cerr << "Compute cell and faces distance to the interface" << finl;
        if (correct_fluxes_ || correct_temperature_cell_neighbours_
            || find_cell_neighbours_for_fluxes_spherical_correction_
            || compute_reachable_fluxes_)
          {
            compute_distance_faces_centres();
            if (!readjust_probe_length_from_vertices_ || !enable_resize_probe_collision_)
              {
                if (debug_)
                  Cerr << "Compute distance cell neighbours" << finl;
                if (correct_temperature_cell_neighbours_ || find_cell_neighbours_for_fluxes_spherical_correction_)
                  compute_distance_cell_centres_neighbours();
                if (debug_)
                  Cerr << "Compute distance faces neighbours" << finl;
                if (compute_reachable_fluxes_)
                  compute_distance_last_cell_faces_neighbours();
              }
          }
      }

  if (fluxes_correction_conservations_)
    {
      compute_pure_liquid_neighbours();
      locate_pure_mixed_neighbours_without_pure_liquid_faces();
    }

  surface_ = (*eulerian_interfacial_area_)(index_i_, index_j_, index_k_);
  // Resize won't change the values if size is not changing...
  reinit_variable(rhs_assembly_);
  // rhs_assembly_.resize(*points_per_thermal_subproblem_);
  // rhs_assembly_ *= 0.;
}

void IJK_One_Dimensional_Subproblem::compute_interface_basis_vectors()
{
  /*
   * TODO: Associate a basis to each subproblem
   * Use Rodrigues' rotation formula to determine ephi ?
   * Needs an axis of (rotation gravity_dir x relative_vectors)
   * and an angle (gravity_dir dot relative_vectors) / (norm(gravity_dir)*norm(relative_vectors))
   * ephi is determined in the gravity_align rising direction
   * 		 | gravity_dir
   * 		 |
   *   *****
   * ***   ***
   * **     **
   * ***   ***
   *   *****
   *     |
   *     |
   */

  facet_barycentre_relative_ = facet_barycentre_ - bubble_barycentre_;
  if (debug_)
    {
      Cerr << "bubble_barycentre_"<< bubble_barycentre_[0] << " ; " << bubble_barycentre_[1] << " ; " << bubble_barycentre_[2] << finl;
      Cerr << "facet_barycentre_"<< facet_barycentre_[0] << " ; " << facet_barycentre_[1] << " ; " << facet_barycentre_[2] << finl;
      Cerr << "facet_barycentre_relative_"<< facet_barycentre_relative_[0] << " ; " << facet_barycentre_relative_[1] << " ; " << facet_barycentre_relative_[2] << finl;
    }
  Vecteur3 facet_barycentre_relative_normed = facet_barycentre_relative_;
  const double facet_barycentre_relative_norm = facet_barycentre_relative_normed.length();
  facet_barycentre_relative_normed *= (1 / facet_barycentre_relative_norm);
  Vecteur3 normal_contrib;
  const double normal_vector_compo_norm = normal_vector_compo_.length();
  normal_vector_compo_ *= (1 / normal_vector_compo_norm);

  if (debug_)
    Cerr << "Normal vector norm:" << normal_vector_compo_norm << finl;
  /*
   * First method with tangential direction of maximum velocity variations
   */
  DoubleTab facet_barycentre(1, 3);
  interfacial_velocity_compo_ = 0.;
  for (int dir=0; dir<3; dir++)
    facet_barycentre(0, dir) = facet_barycentre_[dir];
  for (int dir=0; dir<3; dir++)
    {
      DoubleVect interfacial_velocity_component(1);
      ijk_interpolate_skip_unknown_points((*velocity_)[dir], facet_barycentre, interfacial_velocity_component, INVALID_INTERP);
      interfacial_velocity_compo_[dir] = interfacial_velocity_component[0];
    }
  if (interfacial_velocity_compo_.length() < INVALID_VELOCITY)
    {
      normal_contrib = normal_vector_compo_;
      normal_contrib *= Vecteur3::produit_scalaire(facet_barycentre_relative_normed, normal_vector_compo_);
      first_tangential_vector_compo_ = facet_barycentre_relative_normed - normal_contrib;
    }
  else
    {
      // Should I remove the rising velocity ?
      interfacial_velocity_compo_ = interfacial_velocity_compo_ - bubble_rising_velocity_compo_;
      normal_contrib = normal_vector_compo_;
      normal_contrib *= Vecteur3::produit_scalaire(interfacial_velocity_compo_, normal_vector_compo_);
      interfacial_tangential_velocity_compo_ = interfacial_velocity_compo_ - normal_contrib;
      first_tangential_vector_compo_ = interfacial_tangential_velocity_compo_;
    }
  const double norm_first_tangential_vector = first_tangential_vector_compo_.length();
  first_tangential_vector_compo_ *= (1 / norm_first_tangential_vector);
  Vecteur3::produit_vectoriel(normal_vector_compo_, first_tangential_vector_compo_, second_tangential_vector_compo_);
  const double norm_second_tangential_vector = second_tangential_vector_compo_.length();
  second_tangential_vector_compo_ *= (1 / norm_second_tangential_vector);

  /*
   * Second method with rising velocity
   */
  Vecteur3::produit_vectoriel(bubble_rising_vector_, facet_barycentre_relative_, azymuthal_vector_compo_raw_);

  azymuthal_vector_compo_ = azymuthal_vector_compo_raw_;
  const double norm_azymuthal_vector_compo_raw_ = azymuthal_vector_compo_raw_.length();
//  const int sign_vector = signbit(Vecteur3::produit_scalaire(bubble_rising_vector_, normal_vector_compo_));
//  if (sign_vector)
//    azymuthal_vector_compo_ *= -1;
  azymuthal_vector_compo_ *= (1 / norm_azymuthal_vector_compo_raw_);

  normal_contrib = normal_vector_compo_;
  normal_contrib *=	Vecteur3::produit_scalaire(azymuthal_vector_compo_, normal_vector_compo_);
  azymuthal_vector_compo_ = azymuthal_vector_compo_ - normal_contrib;
  Vecteur3::produit_vectoriel(azymuthal_vector_compo_, normal_vector_compo_, first_tangential_vector_compo_from_rising_dir_);
  const double norm_first_tangential_vector_from_rising_dir = first_tangential_vector_compo_from_rising_dir_.length();
  first_tangential_vector_compo_from_rising_dir_ *= (1 / norm_first_tangential_vector_from_rising_dir);


  if (tangential_from_rising_vel_)
    {
      first_tangential_vector_compo_solver_ = &first_tangential_vector_compo_from_rising_dir_;
      second_tangential_vector_compo_solver_ = &azymuthal_vector_compo_;
      first_tangential_velocity_not_corrected_ = &first_tangential_velocity_from_rising_dir_;
      second_tangential_velocity_not_corrected_ = &azymuthal_velocity_;
      first_tangential_velocity_solver_ = &first_tangential_velocity_from_rising_dir_corrected_;
      second_tangential_velocity_solver_ = &azymuthal_velocity_corrected_;
      tangential_temperature_gradient_first_solver_ = &tangential_temperature_gradient_first_from_rising_dir_;
      tangential_temperature_gradient_second_solver_ = &azymuthal_temperature_gradient_;
    }
  else
    {
      // By default
      first_tangential_vector_compo_solver_= &first_tangential_vector_compo_;
      second_tangential_vector_compo_solver_ = &second_tangential_vector_compo_;
      first_tangential_velocity_not_corrected_ = &first_tangential_velocity_;
      second_tangential_velocity_not_corrected_ = &second_tangential_velocity_;
      first_tangential_velocity_solver_ = &first_tangential_velocity_corrected_;
      second_tangential_velocity_solver_ = &second_tangential_velocity_corrected_;
      tangential_temperature_gradient_first_solver_ = &tangential_temperature_gradient_first_;
      tangential_temperature_gradient_second_solver_ = &tangential_temperature_gradient_second_;
    }
}

void IJK_One_Dimensional_Subproblem::compute_pure_spherical_basis_vectors()
{
  /*
   * FIXME: It is align with gravity z but it should be modified to be align with the gravity dir ?
   */
  if (debug_)
    Cerr << "r_sph_ calculation"  << finl;
  r_sph_ = sqrt(facet_barycentre_relative_[0] * facet_barycentre_relative_[0]
                + facet_barycentre_relative_[1] * facet_barycentre_relative_[1]
                + facet_barycentre_relative_[2] * facet_barycentre_relative_[2]);
  if (debug_)
    {
      Cerr << "r_sph_ = " << r_sph_ << finl;
      Cerr << "theta_sph_ calculation"  << finl;
    }

//  theta_sph_ = atan(sqrt(facet_barycentre_relative_[0] * facet_barycentre_relative_[0]
//                         + facet_barycentre_relative_[1] * facet_barycentre_relative_[1])/ facet_barycentre_relative_[2]);
  theta_sph_ = atan2(sqrt(facet_barycentre_relative_[0] * facet_barycentre_relative_[0]
                          + facet_barycentre_relative_[1] * facet_barycentre_relative_[1]), facet_barycentre_relative_[2]);
  const double atan_theta_incr_ini = M_PI / 2;
  const double atan_incr_factor = -1;
  theta_sph_ = (theta_sph_ - atan_theta_incr_ini) * atan_incr_factor;

  if (debug_probe_collision_)
    if (theta_sph_ < 0)
      disable_probe_because_collision_ = 1;

  if (debug_)
    {
      Cerr << "theta_sph_ = " << theta_sph_ << finl;
      Cerr << "phi_sph_ calculation"  << finl;
    }
  phi_sph_ = atan2(facet_barycentre_relative_[1], facet_barycentre_relative_[0]);

  if (debug_)
    {
      Cerr << "phi_sph_ = " << phi_sph_ << finl;
      Cerr << "er_sph_ calculation"  << finl;
    }
  for (int dir=0; dir<3; dir++)
    er_sph_[dir] = facet_barycentre_relative_[dir] / r_sph_;

  const double length = sqrt(facet_barycentre_relative_[0] * facet_barycentre_relative_[0]
                             + facet_barycentre_relative_[1] * facet_barycentre_relative_[1]);

  if (debug_)
    {
      Cerr << "er_sph_ = " << er_sph_[0] << finl;
      Cerr << "etheta_sph_ calculation"  << finl;
    }
  for (int dir=0; dir<2; dir++)
    etheta_sph_[dir] = facet_barycentre_relative_[dir] * facet_barycentre_relative_[2] / (r_sph_ * length);
  etheta_sph_[2] = - facet_barycentre_relative_[2] * length / r_sph_;

  ephi_sph_ = {0., 0., 0.};
  ephi_sph_[0] = - facet_barycentre_relative_[1];
  ephi_sph_[1] = facet_barycentre_relative_[0];
}

void IJK_One_Dimensional_Subproblem::compute_local_discretisation()
{
  int i;
  if (global_probes_characteristics_)
    {
      if (!probe_variations_enabled_)
        {
          if (!velocities_calculation_counter_)
            {
              radial_coordinates_ = radial_coordinates_base_;
              dr_ = *dr_base_;
            }
        }
      else
        {
          dr_ = probe_length_ / (*points_per_thermal_subproblem_ - 1);
          radial_coordinates_modified_.resize(*points_per_thermal_subproblem_);
          for (i=0; i < *points_per_thermal_subproblem_; i++)
            radial_coordinates_modified_(i) = i * dr_;
          radial_coordinates_ = &radial_coordinates_modified_;
        }
    }
  else
    {
      /*
       * coeff_distance_diagonal_ as well as
       * points_per_thermal_subproblem_ could be adapted
       */
      if (!probe_variations_enabled_)
        radial_coordinates_modified_.resize(*points_per_thermal_subproblem_);
      dr_ = probe_length_ / (*points_per_thermal_subproblem_ - 1);
      for (i=0; i < *points_per_thermal_subproblem_; i++)
        radial_coordinates_modified_(i) = i * dr_;
      radial_coordinates_ = &radial_coordinates_modified_;
    }
  /*
   * Following attributes differ anyway !
   */
  if (!velocities_calculation_counter_ || probe_variations_enabled_)
    {
      dr_inv_ = 1 / dr_;
      osculating_radial_coordinates_ = (*radial_coordinates_);
      osculating_radial_coordinates_ += osculating_radius_;
      if (!probe_variations_enabled_)
        {
          radial_coordinates_cartesian_compo_.resize(*points_per_thermal_subproblem_, 3);
          coordinates_cartesian_compo_.resize(*points_per_thermal_subproblem_, 3);
          osculating_radial_coordinates_cartesian_compo_.resize(*points_per_thermal_subproblem_, 3);
          osculating_radial_coordinates_inv_.resize(*points_per_thermal_subproblem_);
        }
      for (i=0; i < *points_per_thermal_subproblem_; i++)
        {
          osculating_radial_coordinates_inv_[i] = 1 / osculating_radial_coordinates_[i];
          for (int dir=0; dir<3; dir++)
            {
              radial_coordinates_cartesian_compo_(i, dir) = (*radial_coordinates_)(i) * normal_vector_compo_[dir];
              osculating_radial_coordinates_cartesian_compo_(i, dir) = osculating_radial_coordinates_(i) * normal_vector_compo_[dir];
              coordinates_cartesian_compo_(i, dir) = radial_coordinates_cartesian_compo_(i, dir) + facet_barycentre_[dir];
            }
        }
    }
}

void IJK_One_Dimensional_Subproblem::compute_local_time_step()
{
  if (*first_time_step_temporal_)
    {
      double max_u_inv = 1.e20;
      if (max_u_ > INVALID_VELOCITY_CFL)
        max_u_inv = 1 / max_u_;
      local_fourier_time_step_probe_length_ = probe_length_ * probe_length_ * local_fourier_ / (*alpha_) * 0.125; // factor 1/8 in 3D ?
      local_cfl_time_step_probe_length_ = probe_length_ * max_u_inv  * local_cfl_ * 0.5; // factor 1/2 in 3D ?
      local_dt_cfl_min_delta_xyz_ = min_delta_xyz_ * max_u_inv  * local_cfl_ * 0.5;
      if (first_time_step_explicit_)
        {
          local_dt_fo_ = dr_ * dr_ * local_fourier_ / (*alpha_);

          local_dt_cfl_ = dr_ * max_u_inv  * local_cfl_;
          local_time_step_ = std::min(local_dt_cfl_, local_dt_fo_);
          local_time_step_ = std::min(local_time_step_, global_time_step_);
          if (is_first_time_step_)
            nb_iter_explicit_ = (int) (global_time_step_ / local_time_step_);
          else
            nb_iter_explicit_ = (int) ((current_time_ + global_time_step_) / local_time_step_);
          nb_iter_explicit_++;
          if (is_first_time_step_)
            local_time_step_round_ = global_time_step_ / (double) nb_iter_explicit_;
          else
            local_time_step_round_ = (current_time_ + global_time_step_) / (double) nb_iter_explicit_;
          local_time_step_round_ /= (double) (*points_per_thermal_subproblem_);
          nb_iter_explicit_ *= (*points_per_thermal_subproblem_);
        }
      else
        {
          local_time_step_ = global_time_step_;
          local_time_step_round_ = local_time_step_;
        }
    }
}

const int * IJK_One_Dimensional_Subproblem::increase_number_of_points()
{
  increased_point_numbers_ = *points_per_thermal_subproblem_base_;
  return &increased_point_numbers_;
}

void IJK_One_Dimensional_Subproblem::compute_identity_matrix_local(Matrice& identity_matrix_explicit_implicit)
{
  int check_nb_elem;
  check_nb_elem = (*finite_difference_assembler_).build(identity_matrix_explicit_implicit, *points_per_thermal_subproblem_, -1);
  if (debug_)
    Cerr << "Check_nb_elem" << check_nb_elem << finl;
}

void IJK_One_Dimensional_Subproblem::compute_first_order_operator_local(Matrice& radial_first_order_operator)
{
  int check_nb_elem;
  check_nb_elem = (*finite_difference_assembler_).build(radial_first_order_operator, *points_per_thermal_subproblem_, 0);
  if (debug_)
    Cerr << "Check_nb_elem" << check_nb_elem << finl;
  radial_first_order_operator_local_ *= dr_inv_;
}

void IJK_One_Dimensional_Subproblem::compute_second_order_operator_local(Matrice& radial_second_order_operator)
{
  int check_nb_elem;
  check_nb_elem = (*finite_difference_assembler_).build(radial_second_order_operator, *points_per_thermal_subproblem_, 0);
  if (debug_)
    Cerr << "Check_nb_elem" << check_nb_elem << finl;
  const double dr_squared_inv = 1 / pow(dr_, 2);
  radial_second_order_operator_local_ *= dr_squared_inv;
}

void IJK_One_Dimensional_Subproblem::recompute_finite_difference_matrices()
{
  if (!global_probes_characteristics_)
    {
      compute_identity_matrix_local(identity_matrix_explicit_implicit_local_);
      compute_first_order_operator_local(radial_first_order_operator_local_);
      compute_second_order_operator_local(radial_second_order_operator_local_);
      identity_matrix_explicit_implicit_local_ = (*identity_matrix_explicit_implicit_base_);
      identity_matrix_explicit_implicit_ = &identity_matrix_explicit_implicit_local_;
      radial_first_order_operator_ = &radial_first_order_operator_local_;
      radial_second_order_operator_ = &radial_second_order_operator_local_;
    }
}

void IJK_One_Dimensional_Subproblem::compute_first_order_operator_local_varying_probe_length(const Matrice * radial_first_order_operator)
{
  radial_first_order_operator_local_ = (*radial_first_order_operator);
  radial_first_order_operator_local_ *= dr_inv_;
}

void IJK_One_Dimensional_Subproblem::compute_second_order_operator_local_varying_probe_length(const Matrice * radial_second_order_operator)
{
  radial_second_order_operator_local_ = (*radial_second_order_operator);
  const double dr_squared_inv = 1 / pow(dr_, 2);
  radial_second_order_operator_local_ *= dr_squared_inv;
}

void IJK_One_Dimensional_Subproblem::recompute_finite_difference_matrices_varying_probe_length()
{
  compute_first_order_operator_local_varying_probe_length(radial_first_order_operator_raw_base_);
  compute_second_order_operator_local_varying_probe_length(radial_second_order_operator_raw_base_);
  radial_first_order_operator_= &radial_first_order_operator_local_;
  radial_second_order_operator_= &radial_second_order_operator_local_;
}

void IJK_One_Dimensional_Subproblem::interpolate_indicator_on_probes()
{
  if (enable_probe_collision_detection_)
    {
      if (readjust_probe_length_from_vertices_)
        compute_modified_probe_length_vertex_condition();
      indicator_interp_.resize(*points_per_thermal_subproblem_);
      const IJK_Field_double& indicator = interfaces_->I();
      ijk_interpolate_skip_unknown_points(indicator, coordinates_cartesian_compo_, indicator_interp_, INVALID_INTERP);
      if (enable_resize_probe_collision_)
        {
          for (int i=(*points_per_thermal_subproblem_)-1; i>=0; i--)
            {
              const double indic_last = find_cell_related_indicator_on_probes(i);
              if (indic_last > LIQUID_INDICATOR_TEST)
                {
                  resize_probe_collision_index_ = i;
                  if (i != (*points_per_thermal_subproblem_)-1)
                    resize_probe_collision_ = 1;
                  return;
                }
              disable_probe_because_collision_ = (i == 0);
            }
        }
      else
        {
          if (!disable_find_cell_centre_probe_tip_)
            {

              const int last_index = (*points_per_thermal_subproblem_) - 1;
              const double indic_last = find_cell_related_indicator_on_probes(last_index);
              if (indic_last < LIQUID_INDICATOR_TEST)
                {
                  disable_probe_because_collision_ = 1;
                  return;
                }
            }
          else
            {
              for (int i=(*points_per_thermal_subproblem_)-1; i>=0; i--)
                {
                  const double indicator_val = indicator_interp_(i);
                  if (indicator_val < LIQUID_INDICATOR_TEST)
                    {
                      disable_probe_because_collision_ = 1;
                      return;
                    }
                }
            }
        }
    }
}

double IJK_One_Dimensional_Subproblem::find_cell_related_indicator_on_probes(const int& last_index)
{
  const IJK_Field_double& indicator = interfaces_->I();
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();

  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  Vecteur3 xyz_cart_end = {coordinates_cartesian_compo_(last_index, 0),
                           coordinates_cartesian_compo_(last_index, 1),
                           coordinates_cartesian_compo_(last_index, 2)
                          };
  Vecteur3 centre_elem = ref_ijk_ft_->get_splitting_ns().get_coords_of_dof(index_i_, index_j_, index_k_, IJK_Splitting::ELEM);
  Vecteur3 displacement_centre_probe = centre_elem;
  displacement_centre_probe *= (-1);
  displacement_centre_probe += xyz_cart_end;
  Vecteur3 displacement_factor = {displacement_centre_probe[0] / dx,
                                  displacement_centre_probe[1] / dy,
                                  displacement_centre_probe[2] / dz
                                 };
  const int offset_x = (int) displacement_factor[0];
  const int offset_y = (int) displacement_factor[1];
  const int offset_z = (int) displacement_factor[2];
  const int offset_x_elem = ((abs(displacement_factor[0] - offset_x) >= 0.5) ? 1 : 0);
  const int offset_y_elem = ((abs(displacement_factor[1] - offset_y) >= 0.5) ? 1 : 0);
  const int offset_z_elem = ((abs(displacement_factor[2] - offset_z) >= 0.5) ? 1 : 0);
  const int real_offset_x = offset_x + (signbit(offset_x) ? - offset_x_elem : offset_x_elem);
  const int real_offset_y = offset_y + (signbit(offset_y) ? - offset_y_elem : offset_y_elem);
  const int real_offset_z = offset_z + (signbit(offset_z) ? - offset_z_elem : offset_z_elem);
  const double indic_last = indicator(index_i_ + real_offset_x,
                                      index_j_ + real_offset_y,
                                      index_k_ + real_offset_z);
  return indic_last;
}

void IJK_One_Dimensional_Subproblem::interpolate_project_velocities_on_probes()
{
  if (!velocities_calculation_counter_ || probe_variations_enabled_)
    {
      radial_velocity_.resize(*points_per_thermal_subproblem_);
      first_tangential_velocity_.resize(*points_per_thermal_subproblem_);
      second_tangential_velocity_.resize(*points_per_thermal_subproblem_);
      azymuthal_velocity_.resize(*points_per_thermal_subproblem_);
      first_tangential_velocity_from_rising_dir_.resize(*points_per_thermal_subproblem_);

      radial_velocity_corrected_.resize(*points_per_thermal_subproblem_);
      first_tangential_velocity_corrected_.resize(*points_per_thermal_subproblem_);
      second_tangential_velocity_corrected_.resize(*points_per_thermal_subproblem_);
      azymuthal_velocity_corrected_.resize(*points_per_thermal_subproblem_);
      first_tangential_velocity_from_rising_dir_corrected_.resize(*points_per_thermal_subproblem_);

      pressure_interp_.resize(*points_per_thermal_subproblem_);

      x_velocity_.resize(*points_per_thermal_subproblem_);
      y_velocity_.resize(*points_per_thermal_subproblem_);
      z_velocity_.resize(*points_per_thermal_subproblem_);

      x_velocity_corrected_.resize(*points_per_thermal_subproblem_);
      y_velocity_corrected_.resize(*points_per_thermal_subproblem_);
      z_velocity_corrected_.resize(*points_per_thermal_subproblem_);

      interpolate_pressure_on_probes();
      interpolate_cartesian_velocities_on_probes();
      interpolate_velocity_at_cell_centre();
      compute_velocity_magnitude();
      project_velocities_on_probes();
      velocities_calculation_counter_++;
    }
}

void IJK_One_Dimensional_Subproblem::reajust_probe_length()
{
  if (readjust_probe_length_from_vertices_)
    compute_modified_probe_length_condition(1);
  else if (first_time_step_varying_probes_)
    compute_modified_probe_length_condition(2);
  else
    Cerr << "This strategy for readjusting the probe length does not exist" << finl;

}

void IJK_One_Dimensional_Subproblem::compute_modified_probe_length_condition(const int probe_length_condition)
{
  switch(probe_length_condition)
    {
    case 0:
      compute_modified_probe_length_collision();
      break;
    case 1:
      compute_modified_probe_length_vertex_condition();
      break;
    case 2:
      compute_modified_probe_length_temporal_condition();
      break;
    default:
      break;
    }
}

void IJK_One_Dimensional_Subproblem::compute_modified_probe_length_collision()
{
  modified_probe_length_from_collision_ = (*radial_coordinates_)(resize_probe_collision_index_);
}

void IJK_One_Dimensional_Subproblem::compute_modified_probe_length_vertex_condition()
{
  compute_distance_faces_centres();
  bool has_liquid_neighbours = 1;
  for (int i=0; i<6; i++)
    has_liquid_neighbours = has_liquid_neighbours && pure_liquid_neighbours_[i];
  // const double max_distance_pure_face_centre = compute_max_distance_pure_face_centre();
  // const double max_distance_face_centre_vertex = std::max(max_distance_pure_face_centre, max_distance_pure_vertex_centre);
  int lmax, mmax;
  const double max_distance_pure_vertex_centre = compute_max_distance_pure_face_vertices(lmax, mmax);
  //  Vecteur3 normal_contrib = normal_vector_compo_;
  //  normal_contrib *= vertices_centres_distance_[lmax][mmax];
  //  Vecteur3 tangential_distance = vertices_tangential_distance_vector_[lmax][mmax];
  //  tangential_distance *= vertices_centres_tangential_distance_[lmax][mmax];
  //  Vecteur3 facet_to_vertex = facet_to_vertex;
  //  facet_to_vertex += tangential_distance;

  modified_probe_length_from_vertices_ = vertices_centres_distance_[lmax][mmax];
  modified_probe_length_from_vertices_ += (*cell_diagonal_) / 2;
  probe_variations_enabled_ = 1;
  if (debug_)
    Cerr << "Maximum vertex distance:" << max_distance_pure_vertex_centre << finl;
}

void IJK_One_Dimensional_Subproblem::compute_modified_probe_length_temporal_condition()
{
  const double current_time = ref_ijk_ft_->get_current_time();
  cfl_probe_length_ = (max_u_ * current_time) / local_cfl_; // Add 3D constants ?
  fourier_probe_length_ = sqrt(((*alpha_ * (current_time + global_time_step_))) / local_fourier_); // Add 3D constants ?
  max_cfl_fourier_probe_length_ = std::max(cfl_probe_length_, fourier_probe_length_);
  /*
   * TODO: ADD constraint on the temperature because of the displacement of the interface !!!!
   */
  cell_temperature_ = (*temperature_before_extrapolation_)(index_i_, index_j_, index_k_);
  if (max_cfl_fourier_probe_length_ < probe_length_)
    {
      if (debug_)
        Cerr << "Probe length should be modified" << finl;
      if (!correct_fluxes_)
        {
          if (indicator_ < 0.5)
            {
              compute_distance_cell_centre();
              assert(cell_centre_distance_ >= 0.);
              if (max_cfl_fourier_probe_length_ < cell_centre_distance_)
                short_probe_condition_ = 1;
              else
                short_probe_condition_ = 0;
            }
        }
      else
        {
          compute_distance_faces_centres();
          bool has_liquid_neighbours = 1;
          for (int i=0; i<6; i++)
            has_liquid_neighbours = has_liquid_neighbours && pure_liquid_neighbours_[i];
          const double max_distance_pure_face_centre = compute_max_distance_pure_face_centre();
          const double max_distance_pure_vertex_centre = compute_max_distance_pure_face_vertices();
          const double max_distance_face_centre_vertex = std::max(max_distance_pure_face_centre, max_distance_pure_vertex_centre);
          if (max_cfl_fourier_probe_length_ < max_distance_face_centre_vertex || !has_liquid_neighbours)
            {
              short_probe_condition_ = 1;
              if (cell_temperature_ != delta_T_subcooled_overheated_)
                temperature_probe_condition_ = 1;
              else
                temperature_probe_condition_ = 0;
            }
          else
            short_probe_condition_= 0;
        }
      probe_variations_enabled_ = 1;
    }
  else
    probe_variations_enabled_ = 0;
}

/*
 * TODO: Use ijk_intersections_interface instead !
 * and avoid redundancy...
 * Compare Calculation with mean_over_compo()
 * not weighted by the surface
 * Separate geometric attributes and physical attributes -> in the future
 */

void IJK_One_Dimensional_Subproblem::compute_distance_cell_centre()
{
  if (!has_computed_cell_centre_distance_)
    {
      Vecteur3 centre = ref_ijk_ft_->get_splitting_ns().get_coords_of_dof(index_i_, index_j_, index_k_, IJK_Splitting::ELEM);

      Vecteur3 facet_to_cell_centre = facet_barycentre_;
      facet_to_cell_centre *= -1;
      facet_to_cell_centre += centre;

      Vecteur3 facet_to_osculating_cell_centre = normal_vector_compo_;
      facet_to_osculating_cell_centre *= - (- osculating_radius_);
      facet_to_osculating_cell_centre += facet_to_cell_centre;

      cell_centre_distance_ = Vecteur3::produit_scalaire(facet_to_cell_centre, normal_vector_compo_);
      cell_centre_radius_difference_ = facet_to_cell_centre.length() - cell_centre_distance_;
      cell_centre_osculating_radius_difference_ = (facet_to_osculating_cell_centre.length()
                                                   - osculating_radius_ - cell_centre_distance_);

      Vecteur3 normal_contrib = normal_vector_compo_;
      normal_contrib *= cell_centre_distance_;
      Vecteur3 tangential_displacement = normal_contrib;
      tangential_displacement *= (-1);
      tangential_displacement += facet_to_cell_centre;
      cell_centre_tangential_distance_ = tangential_displacement.length();
      tangential_distance_vector_ = tangential_displacement;
      if (cell_centre_tangential_distance_ > 1e-16)
        tangential_distance_vector_ *= (1 / cell_centre_tangential_distance_);
      has_computed_cell_centre_distance_ = true;
    }
  else if (debug_)
    Cerr << "Cell centre distances have already been computed" << finl;
}

void IJK_One_Dimensional_Subproblem::compute_distance_faces_centres()
{
  if (!has_computed_cell_faces_distance_)
    {
      Vecteur3 bary_face {0., 0., .0};
      Vecteur3 vector_relative {0., 0., 0.};
      Vecteur3 bary_vertex {0., 0., 0.};
      const int neighbours_i[6] = NEIGHBOURS_I;
      const int neighbours_j[6] = NEIGHBOURS_J;
      const int neighbours_k[6] = NEIGHBOURS_K;
      const int neighbours_faces_i[6] = NEIGHBOURS_FACES_I;
      const int neighbours_faces_j[6] = NEIGHBOURS_FACES_J;
      const int neighbours_faces_k[6] = NEIGHBOURS_FACES_K;
      const int face_dir[6] = FACES_DIR;
      int m;
      for (int l=0; l<6; l++)
        {
          const int ii = neighbours_i[l];
          const int jj = neighbours_j[l];
          const int kk = neighbours_k[l];
          const double indic_neighbour = ref_ijk_ft_->itfce().I()(index_i_+ii, index_j_+jj, index_k_+kk);
          if (fabs(indic_neighbour) > LIQUID_INDICATOR_TEST)
            {
              const int ii_f = neighbours_faces_i[l];
              const int jj_f = neighbours_faces_j[l];
              const int kk_f = neighbours_faces_k[l];
              pure_vapour_neighbours_[l] = 0;
              pure_liquid_neighbours_[l] = 1;
              if (ii)
                bary_face = ref_ijk_ft_->get_splitting_ns().get_coords_of_dof(index_i_+ii_f, index_j_+jj_f, index_k_+kk_f, IJK_Splitting::FACES_I);
              if (jj)
                bary_face = ref_ijk_ft_->get_splitting_ns().get_coords_of_dof(index_i_+ii_f, index_j_+jj_f, index_k_+kk_f, IJK_Splitting::FACES_J);
              if (kk)
                bary_face = ref_ijk_ft_->get_splitting_ns().get_coords_of_dof(index_i_+ii_f, index_j_+jj_f, index_k_+kk_f, IJK_Splitting::FACES_K);

              // Normal distance
              vector_relative = facet_barycentre_;
              vector_relative *= (-1);
              vector_relative += bary_face;
              {
                const double distance_face_centre = Vecteur3::produit_scalaire(vector_relative, normal_vector_compo_);
                face_centres_radius_difference_[l] = vector_relative.length() - distance_face_centre;
                face_centres_distance_[l] = distance_face_centre;
                // Tangential distance
                Vecteur3 normal_contrib = normal_vector_compo_;
                normal_contrib *= distance_face_centre;
                Vecteur3 tangential_displacement = normal_contrib;
                tangential_displacement *= (-1);
                tangential_displacement += vector_relative;
                face_centres_tangential_distance_[l] = tangential_displacement.length();
                if (face_centres_tangential_distance_[l] > 1e-16)
                  tangential_displacement *= (1 / face_centres_tangential_distance_[l]);
                face_tangential_distance_vector_[l] = tangential_displacement;
              }
              // Distance to vertex
              for (m=0; m<4; m++)
                {
                  double distance_vertex_centre = 0.;
                  double tangential_distance_vertex_centre = 0.;
                  Vecteur3 tangential_distance_vector_vertex_centre = {0., 0., 0.};
                  bary_vertex = vector_relative;
                  compute_vertex_position(m,
                                          face_dir[l],
                                          bary_vertex,
                                          distance_vertex_centre,
                                          tangential_distance_vertex_centre,
                                          tangential_distance_vector_vertex_centre);
                  vertices_centres_distance_[l][m] = distance_vertex_centre;
                  vertices_centres_tangential_distance_[l][m] = tangential_distance_vertex_centre;
                  vertices_tangential_distance_vector_[l][m] = tangential_distance_vector_vertex_centre;
                }
            }
          else
            {
              if (fabs(indic_neighbour) < VAPOUR_INDICATOR_TEST)
                pure_vapour_neighbours_[l] = 1;
              else
                pure_vapour_neighbours_[l] = 0;
              pure_liquid_neighbours_[l] = 0;
              face_centres_distance_[l] = 0.;
              face_centres_tangential_distance_[l] = 0.;
              face_tangential_distance_vector_[l] = {0., 0., 0.};
              for (m=0; m<4; m++)
                {
                  vertices_centres_distance_[l][m] = 0.;
                  vertices_centres_tangential_distance_[l][m] = 0.;
                  vertices_tangential_distance_vector_[l][m] = {0., 0., 0.};
                }
            }
        }
      has_computed_liquid_neighbours_ = true;
      has_computed_cell_faces_distance_ = true;
    }
  else if (debug_)
    Cerr << "Cell face distances have already been computed" << finl;
}

double IJK_One_Dimensional_Subproblem::compute_min_distance_pure_face_centre()
{
  double min_face_centre_distance = 1.e20;
  for (int l=0; l<6; l++)
    if (pure_liquid_neighbours_[l])
      if(face_centres_distance_[l] > 0)
        min_face_centre_distance = std::min(min_face_centre_distance, face_centres_distance_[l]);
  return min_face_centre_distance;
}

double IJK_One_Dimensional_Subproblem::compute_min_distance_pure_face_vertices()
{
  double min_face_vertex_distance = 1.e20;
  int m;
  for (int l=0; l<6; l++)
    if (pure_liquid_neighbours_[l])
      for (m=0; m<4; m++)
        if(vertices_centres_distance_[l][m] > 0)
          min_face_vertex_distance = std::min(min_face_vertex_distance, vertices_centres_distance_[l][m]);
  return min_face_vertex_distance;
}

double IJK_One_Dimensional_Subproblem::compute_max_distance_pure_face_centre()
{
  double max_face_centre_distance = 0.;
  for (int l=0; l<6; l++)
    if (pure_liquid_neighbours_[l])
      if(face_centres_distance_[l] > 0)
        max_face_centre_distance = std::max(max_face_centre_distance, face_centres_distance_[l]);
  return max_face_centre_distance;
}

double IJK_One_Dimensional_Subproblem::compute_max_distance_pure_face_vertices()
{
  double max_face_vertex_distance = 0.;
  int m;
  for (int l=0; l<6; l++)
    if (pure_liquid_neighbours_[l])
      for (m=0; m<4; m++)
        if(vertices_centres_distance_[l][m] > 0)
          max_face_vertex_distance = std::max(max_face_vertex_distance, vertices_centres_distance_[l][m]);
  return max_face_vertex_distance;
}

double IJK_One_Dimensional_Subproblem::compute_max_distance_pure_face_vertices(int& lmax, int& mmax)
{
  double max_face_vertex_distance = 0.;
  int m;
  lmax = 0;
  mmax = 0;
  for (int l=0; l<6; l++)
    if (pure_liquid_neighbours_[l])
      for (m=0; m<4; m++)
        if(vertices_centres_distance_[l][m] > 0)
          {
            max_face_vertex_distance = std::max(max_face_vertex_distance, vertices_centres_distance_[l][m]);
            if (vertices_centres_distance_[l][m] >= max_face_vertex_distance)
              {
                lmax = l;
                mmax = m;
              }
          }
  return max_face_vertex_distance;
}

void IJK_One_Dimensional_Subproblem::compute_vertex_position(const int& vertex_number,
                                                             const int& face_dir,
                                                             Vecteur3& bary_vertex,
                                                             double& distance_vertex_centre,
                                                             double& tangential_distance_vertex_centre,
                                                             Vecteur3& tangential_distance_vector_vertex_centre)
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  const double neighbours_first_dir[4] = NEIGHBOURS_FIRST_DIR;
  const double neighbours_second_dir[4] = NEIGHBOURS_SECOND_DIR;
  Vecteur3 point_coords {0., 0., 0.};
  double dl1;
  double dl2;
  switch(face_dir)
    {
    case 0:
      dl1 = dy / 2.;
      dl2 = dz / 2.;
      point_coords[1] = dl1 * neighbours_first_dir[vertex_number];
      point_coords[2] = dl2 * neighbours_second_dir[vertex_number];
      break;
    case 1:
      dl1 = dx / 2.;
      dl2 = dz / 2.;
      point_coords[0] = dl1 * neighbours_first_dir[vertex_number];
      point_coords[2] = dl2 * neighbours_second_dir[vertex_number];
      break;
    case 2:
      dl1 = dx / 2.;
      dl2 = dy / 2.;
      point_coords[0] = dl1 * neighbours_first_dir[vertex_number];
      point_coords[1] = dl2 * neighbours_second_dir[vertex_number];
      break;
    default:
      dl1 = dx / 2.;
      dl2 = dy / 2.;
      point_coords[0] = dl1 * neighbours_first_dir[vertex_number];
      point_coords[1] = dl2 * neighbours_second_dir[vertex_number];
      break;
    }
  bary_vertex += point_coords;
  distance_vertex_centre = Vecteur3::produit_scalaire(bary_vertex, normal_vector_compo_);
  Vecteur3 tangential_distance_vector = normal_vector_compo_;
  tangential_distance_vector *= distance_vertex_centre;
  tangential_distance_vector *= (-1);
  tangential_distance_vector += bary_vertex;
  tangential_distance_vertex_centre = tangential_distance_vector.length();
  if (tangential_distance_vertex_centre > 1e-16)
    tangential_distance_vector *= (1 / tangential_distance_vertex_centre);
  tangential_distance_vector_vertex_centre = tangential_distance_vector;
}

void IJK_One_Dimensional_Subproblem::compute_distance_cell_centres_neighbours()
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  int l, m, n;

  int dxyz_increment_max = get_dxyz_increment_max();
  /*
   * 8-1 values for one neighbour in each dir... (Too much, enhance later)
   * Positive OR Negative dir depending on the normal vector
   */
  pure_neighbours_to_correct_.resize(dxyz_increment_max + 1);
  pure_neighbours_corrected_distance_.resize(dxyz_increment_max + 1);
  if (neighbours_weighting_)
    pure_neighbours_corrected_colinearity_.resize(dxyz_increment_max + 1);
  for (l=dxyz_increment_max; l>=0; l--)
    {
      pure_neighbours_to_correct_[l].resize(dxyz_increment_max + 1);
      pure_neighbours_corrected_distance_[l].resize(dxyz_increment_max + 1);
      if (neighbours_weighting_)
        pure_neighbours_corrected_colinearity_[l].resize(dxyz_increment_max + 1);
      for (m=dxyz_increment_max; m>=0; m--)
        {
          pure_neighbours_to_correct_[l][m].resize(dxyz_increment_max + 1);
          pure_neighbours_corrected_distance_[l][m].resize(dxyz_increment_max + 1);
          if (neighbours_weighting_)
            pure_neighbours_corrected_colinearity_[l][m].resize(dxyz_increment_max + 1);
          for (n=dxyz_increment_max; n>=0; n--)
            {
              pure_neighbours_to_correct_[l][m][n] = false;
              pure_neighbours_corrected_distance_[l][m][n] = 0.;
              if (neighbours_weighting_)
                pure_neighbours_corrected_colinearity_[l][m][n] = 0.;
            }
        }
    }

  double remaining_distance_diag = probe_length_ - cell_centre_distance_;
  Vecteur3 remaining_distance_diag_projected = normal_vector_compo_;
  remaining_distance_diag_projected *= remaining_distance_diag;
  for (int i=0; i<3; i++)
    pure_neighbours_corrected_sign_[i] = signbit(normal_vector_compo_[i]);
  int dx_increment = (int) abs(remaining_distance_diag_projected[0] / dx);
  int dy_increment = (int) abs(remaining_distance_diag_projected[1] / dy);
  int dz_increment = (int) abs(remaining_distance_diag_projected[2] / dz);
  if (correct_neighbours_rank_)
    {
      dx_increment = std::min(dx_increment, dxyz_increment_max);
      dy_increment = std::min(dy_increment, dxyz_increment_max);
      dz_increment = std::min(dz_increment, dxyz_increment_max);
    }
  dxyz_increment_bool_ = (dx_increment || dy_increment || dz_increment);
  for (l=dx_increment; l>=0; l--)
    for (m=dy_increment; m>=0; m--)
      for (n=dz_increment; n>=0; n--)
        if (l!=0 || m!=0 || n!=0)
          {
            const int l_dir = (pure_neighbours_corrected_sign_[0]) ? l * (-1) : l;
            const int m_dir = (pure_neighbours_corrected_sign_[1]) ? m * (-1) : m;
            const int n_dir = (pure_neighbours_corrected_sign_[2]) ? n * (-1) : n;
            const double indic_neighbour = ref_ijk_ft_->itfce().I()(index_i_ + l_dir, index_j_ + m_dir, index_k_ + n_dir);
            if (indic_neighbour > LIQUID_INDICATOR_TEST)
              {
                pure_neighbours_to_correct_[l][m][n] = true;
                const double dx_contrib = l_dir * dx;
                const double dy_contrib = m_dir * dy;
                const double dz_contrib = n_dir * dz;
                Vecteur3 distance_contrib = {dx_contrib, dy_contrib, dz_contrib};
                pure_neighbours_corrected_distance_[l][m][n] = cell_centre_distance_
                                                               + Vecteur3::produit_scalaire(normal_vector_compo_, distance_contrib);
                if (neighbours_weighting_)
                  {
                    const double colinearity = compute_cell_weighting(dx_contrib, dy_contrib, dz_contrib);
                    pure_neighbours_corrected_colinearity_[l][m][n] = colinearity;
                  }
              }
          }
}

double IJK_One_Dimensional_Subproblem::compute_cell_weighting(const double& dx_contrib,
                                                              const double& dy_contrib,
                                                              const double& dz_contrib)
{
  if (neighbours_colinearity_weighting_)
    return compute_colinearity(dx_contrib, dy_contrib, dz_contrib);
  if (neighbours_distance_weighting_)
    return compute_distance_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  if (neighbours_colinearity_distance_weighting_)
    return compute_colinearity(dx_contrib, dy_contrib, dz_contrib) * compute_distance_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  return 1;
}

void IJK_One_Dimensional_Subproblem::compute_distance_last_cell_faces_neighbours()
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  const double dx_over_two = dx / 2.;
  const double dy_over_two = dy / 2.;
  const double dz_over_two = dz / 2.;
  int l, m, n;
  int l_cell, m_cell, n_cell;

  int dxyz_increment_max = get_dxyz_increment_max();
  int dxyz_over_two_increment_max = get_dxyz_over_two_increment_max();
  const int first_increment[3] = {dxyz_over_two_increment_max + 1, dxyz_increment_max, dxyz_increment_max};
  const int second_increment[3] = {dxyz_increment_max, dxyz_over_two_increment_max + 1, dxyz_increment_max};
  const int third_increment[3] = {dxyz_increment_max, dxyz_increment_max, dxyz_over_two_increment_max + 1};
  //  dxyz_over_two_increment_max *= 2;
  //  if (!dxyz_over_two_increment_max%2)
  //    dxyz_over_two_increment_max -= 1;

  /*
   * 8-1 values for one neighbour in each dir... (Too much, enhance later)
   * Positive OR Negative dir depending on the normal vector
   */
  pure_neighbours_last_faces_to_correct_.resize(3);
  pure_neighbours_last_faces_corrected_distance_.resize(3);
  //	if (neighbours_last_faces_weighting_) ?
  pure_neighbours_last_faces_corrected_colinearity_.resize(3);
  for (int c=0; c<3; c++)
    {
      const int first_incr = first_increment[c];
      const int second_incr = second_increment[c];
      const int third_incr = third_increment[c];
      pure_neighbours_last_faces_to_correct_[c].resize(first_incr + 1);
      pure_neighbours_last_faces_corrected_distance_[c].resize(first_incr + 1);
      //	if (neighbours_last_faces_weighting_) ?
      pure_neighbours_last_faces_corrected_colinearity_[c].resize(first_incr + 1);
      for (l=first_incr; l>=0; l--)
        {
          pure_neighbours_last_faces_to_correct_[c][l].resize(second_incr + 1);
          pure_neighbours_last_faces_corrected_distance_[c][l].resize(second_incr + 1);
          //	if (neighbours_last_faces_weighting_) ?
          pure_neighbours_last_faces_corrected_colinearity_[c][l].resize(second_incr + 1);
          for (m=second_incr; m>=0; m--)
            {
              pure_neighbours_last_faces_to_correct_[c][l][m].resize(third_incr + 1);
              pure_neighbours_last_faces_corrected_distance_[c][l][m].resize(third_incr + 1);
              //	if (neighbours_last_faces_weighting_) ?
              pure_neighbours_last_faces_corrected_colinearity_[c][l][m].resize(third_incr + 1);
              for (n=third_incr; n>=0; n--)
                {
                  pure_neighbours_last_faces_to_correct_[c][l][m][n] = false;
                  pure_neighbours_last_faces_corrected_distance_[c][l][m][n] = 0.;
                  //	if (neighbours_last_faces_weighting_) ?
                  pure_neighbours_last_faces_corrected_colinearity_[c][l][m][n] = 0.;
                }
            }
        }
    }

  double remaining_distance_diag = probe_length_ - cell_centre_distance_;
  Vecteur3 remaining_distance_diag_projected = normal_vector_compo_;
  remaining_distance_diag_projected *= remaining_distance_diag;
  for (int i=0; i<3; i++)
    pure_neighbours_corrected_sign_[i] = signbit(normal_vector_compo_[i]);
  int dx_over_two_increment = (int) abs(remaining_distance_diag_projected[0] / dx_over_two);
  int dy_over_two_increment = (int) abs(remaining_distance_diag_projected[1] / dy_over_two);
  int dz_over_two_increment = (int) abs(remaining_distance_diag_projected[2] / dz_over_two);
  int dx_increment = (int) (dx_over_two_increment / 2);
  int dy_increment = (int) (dy_over_two_increment / 2);
  int dz_increment = (int) (dz_over_two_increment / 2);
  dx_over_two_increment -= dx_increment;
  dy_over_two_increment -= dy_increment;
  dz_over_two_increment -= dz_increment;
  if (correct_neighbours_rank_)
    {
      dx_over_two_increment = std::min(dx_over_two_increment, dxyz_over_two_increment_max);
      dy_over_two_increment = std::min(dy_over_two_increment, dxyz_over_two_increment_max);
      dz_over_two_increment = std::min(dz_over_two_increment, dxyz_over_two_increment_max);
      dx_increment = std::min(dx_increment, dxyz_increment_max);
      dy_increment = std::min(dy_increment, dxyz_increment_max);
      dz_increment = std::min(dz_increment, dxyz_increment_max);
    }
  dxyz_over_two_increment_bool_ = (dx_over_two_increment >0 || dy_over_two_increment>0 || dz_over_two_increment >0 );
  if (dx_over_two_increment>0)
    dx_over_two_increment--;
  if (dy_over_two_increment>0)
    dy_over_two_increment--;
  if (dz_over_two_increment>0)
    dz_over_two_increment--;

  /*
   * TODO: Should we look to cells faces that are not in the normal vector direction ?
   */
  for (l=dx_over_two_increment + 1; l>=0; l--)
    for (m_cell=dy_increment; m_cell>=0; m_cell--)
      for (n_cell=dz_increment; n_cell>=0; n_cell--)
        {
          const int l_dir = (pure_neighbours_corrected_sign_[0]) ? l * (-1) + 1 : l;
          const int m_dir = (pure_neighbours_corrected_sign_[1]) ? m_cell * (-1) : m_cell;
          const int n_dir = (pure_neighbours_corrected_sign_[2]) ? n_cell * (-1) : n_cell;
          const int l_dir_elem = (pure_neighbours_corrected_sign_[0]) ? (l + 1) * (-1) + 1 : l;
          const double indic_neighbour = ref_ijk_ft_->itfce().I()(index_i_ + l_dir_elem, index_j_ + m_dir, index_k_ + n_dir);
          if (indic_neighbour > LIQUID_INDICATOR_TEST)
            {
              pure_neighbours_last_faces_to_correct_[0][l][m_cell][n_cell] = true;
              const double lmn_zero = (l > 0) ? 1. : 0.;
              const double contrib_factor = (pure_neighbours_corrected_sign_[0]) ? (lmn_zero * (2 * abs(l_dir) + 1) - (1. - lmn_zero)) * (-1):
                                            lmn_zero * (2 * (l_dir - 1) + 1) - (1. - lmn_zero);
              const double dx_contrib = contrib_factor * dx_over_two;
              const double dy_contrib = m_dir * dy;
              const double dz_contrib = n_dir * dz;
              Vecteur3 distance_contrib = {dx_contrib, dy_contrib, dz_contrib};
              pure_neighbours_last_faces_corrected_distance_[0][l][m_cell][n_cell] = cell_centre_distance_
                                                                                     + Vecteur3::produit_scalaire(normal_vector_compo_, distance_contrib);
              if (pure_neighbours_last_faces_corrected_distance_[0][l][m_cell][n_cell] < 0)
                pure_neighbours_last_faces_to_correct_[0][l][m_cell][n_cell] = false;
              if (neighbours_last_faces_weighting_)
                {
                  const double colinearity = compute_cell_faces_weighting(dx_contrib, dy_contrib, dz_contrib, 0);
                  pure_neighbours_last_faces_corrected_colinearity_[0][l][m_cell][n_cell] = colinearity;
                }
            }
        }

  for (l_cell=dx_increment; l_cell>=0; l_cell--)
    for (m=dy_over_two_increment + 1; m>=0; m--)
      for (n_cell=dz_increment; n_cell>=0; n_cell--)
        {
          const int l_dir = (pure_neighbours_corrected_sign_[0]) ? l_cell * (-1) : l_cell;
          const int m_dir = (pure_neighbours_corrected_sign_[1]) ? m * (-1) + 1: m;
          const int n_dir = (pure_neighbours_corrected_sign_[2]) ? n_cell * (-1) : n_cell;
          const int m_dir_elem = (pure_neighbours_corrected_sign_[1]) ? (m + 1) * (-1) + 1 : m;
          const double indic_neighbour = ref_ijk_ft_->itfce().I()(index_i_ + l_dir, index_j_ + m_dir_elem, index_k_ + n_dir);
          if (indic_neighbour > LIQUID_INDICATOR_TEST)
            {
              pure_neighbours_last_faces_to_correct_[1][l_cell][m][n_cell] = true;
              const double lmn_zero = (m > 0) ? 1. : 0.;
              const double contrib_factor = (pure_neighbours_corrected_sign_[1]) ? (lmn_zero * (2 * abs(m_dir) + 1) - (1. - lmn_zero)) * (-1):
                                            lmn_zero * (2 * (m_dir - 1) + 1) - (1. - lmn_zero);
              const double dx_contrib = l_dir * dx;
              const double dy_contrib = contrib_factor * dy_over_two;
              const double dz_contrib = n_dir * dz;
              Vecteur3 distance_contrib = {dx_contrib, dy_contrib, dz_contrib};
              pure_neighbours_last_faces_corrected_distance_[1][l_cell][m][n_cell] = cell_centre_distance_
                                                                                     + Vecteur3::produit_scalaire(normal_vector_compo_, distance_contrib);
              if (pure_neighbours_last_faces_corrected_distance_[1][l_cell][m][n_cell] < 0)
                pure_neighbours_last_faces_to_correct_[1][l_cell][m][n_cell] = false;
              if (neighbours_last_faces_weighting_)
                {
                  const double colinearity = compute_cell_faces_weighting(dx_contrib, dy_contrib, dz_contrib, 1);
                  pure_neighbours_last_faces_corrected_colinearity_[1][l_cell][m][n_cell] = colinearity;
                }
            }
        }

  for (l_cell=dx_increment; l_cell>=0; l_cell--)
    for (m_cell=dy_increment; m_cell>=0; m_cell--)
      for (n=dz_over_two_increment + 1; n>=0; n--)
        {
          const int l_dir = (pure_neighbours_corrected_sign_[0]) ? l_cell * (-1) : l_cell;
          const int m_dir = (pure_neighbours_corrected_sign_[1]) ? m_cell * (-1) : m_cell;
          const int n_dir = (pure_neighbours_corrected_sign_[2]) ? n * (-1) + 1: n;
          const int n_dir_elem = (pure_neighbours_corrected_sign_[2]) ? (n + 1) * (-1) + 1 : n;
          const double indic_neighbour = ref_ijk_ft_->itfce().I()(index_i_ + l_dir, index_j_ + m_dir, index_k_ + n_dir_elem);
          if (indic_neighbour > LIQUID_INDICATOR_TEST)
            {
              pure_neighbours_last_faces_to_correct_[2][l_cell][m_cell][n] = true;
              const double lmn_zero = (n > 0) ? 1. : 0.;
              const double contrib_factor = (pure_neighbours_corrected_sign_[2]) ? (lmn_zero * (2 * abs(n_dir) + 1) - (1. - lmn_zero)) * (-1):
                                            lmn_zero * (2 * (n_dir - 1) + 1) - (1. - lmn_zero);
              const double dx_contrib = l_dir * dx;
              const double dy_contrib = m_dir * dy;
              const double dz_contrib = contrib_factor * dz_over_two;
              Vecteur3 distance_contrib = {dx_contrib, dy_contrib, dz_contrib};
              pure_neighbours_last_faces_corrected_distance_[2][l_cell][m_cell][n] = cell_centre_distance_ +
                                                                                     Vecteur3::produit_scalaire(normal_vector_compo_, distance_contrib);
              if (pure_neighbours_last_faces_corrected_distance_[2][l_cell][m_cell][n] < 0)
                pure_neighbours_last_faces_to_correct_[2][l_cell][m_cell][n] = false;
              if (neighbours_last_faces_weighting_)
                {
                  const double colinearity = compute_cell_faces_weighting(dx_contrib, dy_contrib, dz_contrib, 2);
                  pure_neighbours_last_faces_corrected_colinearity_[2][l_cell][m_cell][n] = colinearity;
                }
            }
        }
  //  if (neighbours_last_faces_colinearity_weighting_)
  //    {
  //      Vecteur3 relative_vector = normal_vector_compo_;
  //      relative_vector *= cell_centre_distance_;
  //      relative_vector[0] += ((l + 1) * normal_vector_compo_[0] * dx_over_two);
  //      relative_vector[1] += ((m + 1) * normal_vector_compo_[1] * dy_over_two);
  //      relative_vector[2] += ((n + 1) * normal_vector_compo_[2] * dz_over_two);
  //      const double relative_vector_norm = relative_vector.length();
  //      relative_vector *= (1 / relative_vector_norm);
  //      // pure_neighbours_corrected_colinearity_[l][m][n] = Vecteur3::produit_scalaire(normal_vector_compo_, relative_vector);
  //    }
}

double IJK_One_Dimensional_Subproblem::compute_cell_faces_weighting(const double& dx_contrib,
                                                                    const double& dy_contrib,
                                                                    const double& dz_contrib,
                                                                    const int& dir)
{
  if (neighbours_last_faces_colinearity_weighting_)
    return compute_colinearity(dx_contrib, dy_contrib, dz_contrib);
  if (neighbours_last_faces_colinearity_face_weighting_)
    return compute_colinearity_cell_faces(dx_contrib, dy_contrib, dz_contrib, dir);
  if (neighbours_last_faces_distance_weighting_)
    return compute_distance_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  if (neighbours_last_faces_distance_colinearity_weighting_)
    return compute_colinearity(dx_contrib, dy_contrib, dz_contrib) * compute_distance_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  if (neighbours_last_faces_distance_colinearity_face_weighting_)
    return compute_colinearity_cell_faces(dx_contrib, dy_contrib, dz_contrib, dir) * compute_distance_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  return 1;
}

Vecteur3 IJK_One_Dimensional_Subproblem::compute_relative_vector_cell_faces(const double& dx_contrib,
                                                                            const double& dy_contrib,
                                                                            const double& dz_contrib)
{
  Vecteur3 relative_vector = normal_vector_compo_;
  relative_vector *= cell_centre_distance_;
  Vecteur3 tangential_relative_vector = tangential_distance_vector_;
  tangential_relative_vector *= cell_centre_tangential_distance_;
  relative_vector += tangential_relative_vector;
  Vecteur3 dxyz_contrib = {dx_contrib, dy_contrib, dz_contrib};
  relative_vector += dxyz_contrib;
  return relative_vector;
}

double IJK_One_Dimensional_Subproblem::compute_colinearity(const double& dx_contrib, const double& dy_contrib, const double& dz_contrib)
{
  Vecteur3 relative_vector = compute_relative_vector_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  const double relative_vector_norm = relative_vector.length();
  relative_vector *= (1 / relative_vector_norm);
  const double colinearity = Vecteur3::produit_scalaire(normal_vector_compo_, relative_vector);
  return abs(colinearity);
}

double IJK_One_Dimensional_Subproblem::compute_colinearity_cell_faces(const double& dx_contrib,
                                                                      const double& dy_contrib,
                                                                      const double& dz_contrib,
                                                                      const int& dir)
{
  Vecteur3 relative_vector = compute_relative_vector_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  const double relative_vector_norm = relative_vector.length();
  relative_vector *= (1 / relative_vector_norm);
  return abs(relative_vector[dir]);
}

double IJK_One_Dimensional_Subproblem::compute_distance_cell_faces(const double& dx_contrib,
                                                                   const double& dy_contrib,
                                                                   const double& dz_contrib)
{
  Vecteur3 relative_vector = compute_relative_vector_cell_faces(dx_contrib, dy_contrib, dz_contrib);
  const double distance = Vecteur3::produit_scalaire(tangential_distance_vector_, relative_vector);
  return abs(1 / (distance + 1e-16));
}

int IJK_One_Dimensional_Subproblem::get_dxyz_increment_max()
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);

  int dxyz_increment_max;
  if (!correct_neighbours_rank_)
    {
      const int dx_increment_max = (int) ((dx + probe_length_) / dx);
      const int dy_increment_max = (int) ((dy + probe_length_) / dy);
      const int dz_increment_max = (int) ((dz + probe_length_) / dz);
      dxyz_increment_max = std::max(std::max(dx_increment_max, dy_increment_max), dz_increment_max);
    }
  else
    {
      dxyz_increment_max = neighbours_corrected_rank_;
    }
  return dxyz_increment_max;
}

int IJK_One_Dimensional_Subproblem::get_dxyz_over_two_increment_max()
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);

  int dxyz_over_two_increment_max;
  if (!correct_neighbours_rank_)
    {
      const int dx_increment_max = (int) ((probe_length_ + dx) / (dx / 2.));
      const int dy_increment_max = (int) ((probe_length_ + dy) / (dy / 2.));
      const int dz_increment_max = (int) ((probe_length_ + dz) / (dz / 2.));
      dxyz_over_two_increment_max = std::max(std::max(dx_increment_max, dy_increment_max), dz_increment_max);
    }
  else
    {
      dxyz_over_two_increment_max = neighbours_face_corrected_rank_;
    }
  return dxyz_over_two_increment_max;
}

void IJK_One_Dimensional_Subproblem::compute_modified_probe_length(const int& probe_variations_enabled)
{
  if (probe_variations_enabled && probe_variations_enabled_)
    {
      if (enable_resize_probe_collision_ && readjust_probe_length_from_vertices_)
        {
          double min_probe_length_ = std::min(modified_probe_length_from_collision_, modified_probe_length_from_vertices_);
          if (min_probe_length_ == modified_probe_length_from_vertices_)
            probe_length_ = min_probe_length_;
          else
            {
              disable_probe_because_collision_ = 1;
              return;
            }
        }
      else if (enable_resize_probe_collision_)
        probe_length_ = modified_probe_length_from_collision_;
      else if (readjust_probe_length_from_vertices_)
        probe_length_ = modified_probe_length_from_vertices_;
      else if (first_time_step_varying_probes_)
        probe_length_ = max_cfl_fourier_probe_length_;
      else
        probe_length_ = (*coeff_distance_diagonal_) * (*cell_diagonal_);

      compute_local_discretisation();
      recompute_finite_difference_matrices_varying_probe_length();

      if (readjust_probe_length_from_vertices_ || enable_resize_probe_collision_)
        {
          clear_vectors();
          if (debug_)
            Cerr << "Compute distance cell neighbours" << finl;
          if (correct_temperature_cell_neighbours_ || find_cell_neighbours_for_fluxes_spherical_correction_)
            compute_distance_cell_centres_neighbours();
          if (debug_)
            Cerr << "Compute distance faces neighbours" << finl;
          if (compute_reachable_fluxes_)
            compute_distance_last_cell_faces_neighbours();
        }
    }
  else
    {
      probe_variations_enabled_ = 0;
      first_time_step_varying_probes_ = 0;
    }
}


void IJK_One_Dimensional_Subproblem::interpolate_velocity_at_cell_centre()
{
  Vecteur3 bary_cell = ref_ijk_ft_->get_splitting_ns().get_coords_of_dof(index_i_, index_j_, index_k_, IJK_Splitting::ELEM);
  xyz_velocity_cell_ = {0., 0., 0.};
  double x_velocity_cell, y_velocity_cell, z_velocity_cell;
  interpolate_quantities_at_point((*velocity_)[0], bary_cell, x_velocity_cell);
  interpolate_quantities_at_point((*velocity_)[1], bary_cell, y_velocity_cell);
  interpolate_quantities_at_point((*velocity_)[2], bary_cell, z_velocity_cell);
  xyz_velocity_cell_ = {x_velocity_cell, y_velocity_cell, z_velocity_cell};
}

void IJK_One_Dimensional_Subproblem::interpolate_quantities_at_point(const IJK_Field_double& eulerian_field, const Vecteur3& compo_xyz, double& field_value)
{
  DoubleVect field_interp(1);
  DoubleTab coordinates_point = get_single_point_coordinates(compo_xyz);
  ijk_interpolate_skip_unknown_points(eulerian_field, coordinates_point, field_interp, INVALID_INTERP);
  field_value = field_interp(0);
}

DoubleTab IJK_One_Dimensional_Subproblem::get_single_point_coordinates(const Vecteur3& compo_xyz)
{
  DoubleTab coordinates_point;
  coordinates_point.resize(1,3);
  for (int c=0; c<3; c++)
    coordinates_point(0,c) = compo_xyz[c];
  return coordinates_point;
}

void IJK_One_Dimensional_Subproblem::interpolate_pressure_on_probes()
{
  ijk_interpolate_skip_unknown_points((*pressure_), coordinates_cartesian_compo_, pressure_interp_, INVALID_INTERP);
}

void IJK_One_Dimensional_Subproblem::interpolate_cartesian_velocities_on_probes()
{
  DoubleVect * vel_compo = nullptr;
  for (int dir = 0; dir < 3; dir++)
    {
      switch (dir)
        {
        case 0:
          vel_compo = &x_velocity_;
          break;
        case 1:
          vel_compo = &y_velocity_;
          break;
        case 2:
          vel_compo = &z_velocity_;
          break;
        default:
          vel_compo = &x_velocity_;
          break;
        }
      ijk_interpolate_skip_unknown_points((*velocity_)[dir], coordinates_cartesian_compo_, *vel_compo, INVALID_INTERP);
    }
}

void IJK_One_Dimensional_Subproblem::compute_velocity_magnitude()
{
  velocity_magnitude_ = x_velocity_;
  velocity_magnitude_ *= x_velocity_;
  DoubleVect velocity_magnitude_y = y_velocity_;
  velocity_magnitude_y *= y_velocity_;
  velocity_magnitude_ += velocity_magnitude_y;
  DoubleVect velocity_magnitude_z = z_velocity_;
  velocity_magnitude_z *= z_velocity_;
  velocity_magnitude_ += velocity_magnitude_z;
  for(int i=0; i<(*points_per_thermal_subproblem_); i++)
    {
      const double velocity_magnitude_sqrt = sqrt(velocity_magnitude_[i]);
      velocity_magnitude_[i] = velocity_magnitude_sqrt;
    }
}

void IJK_One_Dimensional_Subproblem::project_velocities_on_probes()
{
  project_cartesian_onto_basis_vector(x_velocity_, y_velocity_, z_velocity_, normal_vector_compo_, radial_velocity_);
  project_cartesian_onto_basis_vector(x_velocity_, y_velocity_, z_velocity_, first_tangential_vector_compo_, first_tangential_velocity_);
  project_cartesian_onto_basis_vector(x_velocity_, y_velocity_, z_velocity_, second_tangential_vector_compo_, second_tangential_velocity_);
  project_cartesian_onto_basis_vector(x_velocity_, y_velocity_, z_velocity_, azymuthal_vector_compo_, azymuthal_velocity_);
  project_cartesian_onto_basis_vector(x_velocity_, y_velocity_, z_velocity_, first_tangential_vector_compo_from_rising_dir_,
                                      first_tangential_velocity_from_rising_dir_);
  //  project_cartesian_onto_basis_vector(x_velocity_, y_velocity_, z_velocity_, *first_tangential_vector_compo_solver_, first_tangential_velocity_);
  //	project_cartesian_onto_basis_vector(x_velocity_, y_velocity_, z_velocity_, *second_tangential_vector_compo_solver_, second_tangential_velocity_);

  correct_velocities();

  max_u_ = INVALID_VELOCITY_CFL * 0.9;
  for (int i=0; i<radial_velocity_corrected_.size(); i++)
    if (max_u_cartesian_)
      max_u_ = std::max(max_u_, fabs(velocity_magnitude_[i]));
    else
      max_u_ = std::max(max_u_, fabs(radial_velocity_corrected_[i]));
  compute_local_time_step();

  reinit_variable(x_velocity_corrected_);
  reinit_variable(y_velocity_corrected_);
  reinit_variable(z_velocity_corrected_);
  project_basis_vector_onto_cartesian_dir(0, first_tangential_velocity_, second_tangential_velocity_, radial_velocity_corrected_,
                                          *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                          x_velocity_corrected_);
  project_basis_vector_onto_cartesian_dir(1, first_tangential_velocity_, second_tangential_velocity_, radial_velocity_corrected_,
                                          *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                          y_velocity_corrected_);
  project_basis_vector_onto_cartesian_dir(2, first_tangential_velocity_, second_tangential_velocity_, radial_velocity_corrected_,
                                          *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                          z_velocity_corrected_);
}

void IJK_One_Dimensional_Subproblem::correct_velocities()
{
  correct_velocity(radial_velocity_, radial_velocity_advected_frame_);
  correct_velocity(first_tangential_velocity_, first_tangential_velocity_advected_frame_);
  correct_velocity(second_tangential_velocity_, second_tangential_velocity_advected_frame_);
  correct_velocity(first_tangential_velocity_from_rising_dir_, first_tangential_velocity_from_rising_dir_advected_frame_);
  correct_velocity(azymuthal_velocity_, azymuthal_velocity_advected_frame_);

  correct_velocity_rise(radial_velocity_, normal_vector_compo_, radial_velocity_static_frame_);
  correct_velocity_rise(first_tangential_velocity_, first_tangential_vector_compo_, first_tangential_velocity_static_frame_);
  correct_velocity_rise(second_tangential_velocity_, second_tangential_vector_compo_, second_tangential_velocity_static_frame_);
  correct_velocity_rise(first_tangential_velocity_from_rising_dir_, first_tangential_vector_compo_from_rising_dir_, first_tangential_velocity_from_rising_dir_static_frame_);
  correct_velocity_rise(azymuthal_velocity_, azymuthal_vector_compo_, azymuthal_velocity_static_frame_);

  if (advected_frame_of_reference_)
    {
      radial_velocity_corrected_ = radial_velocity_advected_frame_;
      first_tangential_velocity_corrected_ = first_tangential_velocity_advected_frame_;
      second_tangential_velocity_corrected_ = second_tangential_velocity_advected_frame_;
      first_tangential_velocity_from_rising_dir_corrected_ = first_tangential_velocity_from_rising_dir_advected_frame_;
      azymuthal_velocity_corrected_ = azymuthal_velocity_advected_frame_;
    }
  else
    {
      if (neglect_frame_of_reference_radial_advection_)
        radial_velocity_corrected_ = radial_velocity_static_frame_;
      else
        radial_velocity_corrected_ = radial_velocity_advected_frame_;
      first_tangential_velocity_corrected_ = first_tangential_velocity_static_frame_;
      second_tangential_velocity_corrected_ = second_tangential_velocity_static_frame_;
      first_tangential_velocity_from_rising_dir_corrected_ = first_tangential_velocity_from_rising_dir_static_frame_;
      azymuthal_velocity_corrected_ = azymuthal_velocity_static_frame_;
    }

  if (compute_radial_displacement_)
    {
      radial_displacement_over_time_step_ = 0.5 * global_time_step_ * radial_velocity_[0]; // or * 0.5 ? Like crank nicholson ?
      if (distance_cell_faces_from_lrs_)
        {
          cell_centre_distance_corrected_ = cell_centre_distance_ - radial_displacement_over_time_step_;
          cell_centre_distance_corrected_ = cell_centre_distance_corrected_ < 0 ? cell_centre_distance_: cell_centre_distance_corrected_;
          if (correct_fluxes_ || correct_temperature_cell_neighbours_
              || find_cell_neighbours_for_fluxes_spherical_correction_
              || compute_reachable_fluxes_)
            {
              face_centres_distance_corrected_ = face_centres_distance_;
              for (int l=0; l<6; l++)
                {
                  face_centres_distance_corrected_[l] += (- radial_displacement_over_time_step_);
                  face_centres_distance_corrected_[l] = face_centres_distance_corrected_[l] < 0 ? face_centres_distance_[l] : face_centres_distance_corrected_[l];
                }
            }
        }
    }
}

void IJK_One_Dimensional_Subproblem::correct_velocity(const DoubleVect& velocity, DoubleVect& velocity_corrected)
{
  velocity_corrected = velocity;
  for (int i=0; i<velocity_corrected.size(); i++)
    velocity_corrected[i] -= velocity[0];
}

void IJK_One_Dimensional_Subproblem::correct_velocity_rise(const DoubleVect& velocity, const Vecteur3& basis, DoubleVect& velocity_corrected)
{
  DoubleVect bubble_rising_velocity_projection(1);
  DoubleVect bubble_rising_velocity_compo_x(1);
  DoubleVect bubble_rising_velocity_compo_y(1);
  DoubleVect bubble_rising_velocity_compo_z(1);
  bubble_rising_velocity_compo_x[0] = bubble_rising_velocity_compo_[0];
  bubble_rising_velocity_compo_y[0] = bubble_rising_velocity_compo_[1];
  bubble_rising_velocity_compo_z[0] = bubble_rising_velocity_compo_[2];
  project_cartesian_onto_basis_vector(bubble_rising_velocity_compo_x,bubble_rising_velocity_compo_y, bubble_rising_velocity_compo_z,
                                      basis, bubble_rising_velocity_projection);
  velocity_corrected = velocity;
  for (int i=0; i<velocity_corrected.size(); i++)
    velocity_corrected[i] -= bubble_rising_velocity_projection[0];
}

void IJK_One_Dimensional_Subproblem::correct_radial_velocity_probe()
{
  correct_velocity(radial_velocity_, radial_velocity_corrected_);
}

void IJK_One_Dimensional_Subproblem::project_cartesian_onto_basis_vector(const DoubleVect& compo_x,
                                                                         const DoubleVect& compo_y,
                                                                         const DoubleVect& compo_z,
                                                                         const Vecteur3& basis,
                                                                         DoubleVect& projection)
{
  const int size_x = compo_x.size();
  for (int i=0; i<size_x; i++)
    projection[i] = compo_x[i] * basis[0] + compo_y[i] * basis[1] + compo_z[i] * basis[2];
}

void IJK_One_Dimensional_Subproblem::project_basis_vector_onto_cartesian_dir(const int& dir,
                                                                             const DoubleVect& compo_u,
                                                                             const DoubleVect& compo_v,
                                                                             const DoubleVect& compo_w,
                                                                             const Vecteur3& basis_u,
                                                                             const Vecteur3& basis_v,
                                                                             const Vecteur3& basis_w,
                                                                             DoubleVect& projection)
{
  const int size_u = compo_u.size();
  const int size_v = compo_v.size();
  const int size_w = compo_w.size();
  for (int i=0; i<projection.size(); i++)
    {
      if (i< size_u)
        projection[i] += (compo_u[i] * basis_u[dir]);
      if (i< size_v)
        projection[i] += (compo_v[i] * basis_v[dir]);
      if (i< size_w)
        projection[i] += (compo_w[i] * basis_w[dir]);
    }
}

void IJK_One_Dimensional_Subproblem::compute_integral_quantity(DoubleVect& quantity, double& integrated_quantity)
{
  double integrated_quantity_tmp = 0.;
  for (int i=0; i<quantity.size()-1; i++)
    {
      integrated_quantity_tmp += 0.5*(quantity[i] + quantity[i+1]) * dr_;
    }
  integrated_quantity = integrated_quantity_tmp;
}

void IJK_One_Dimensional_Subproblem::compute_integral_quantity_on_probe(DoubleVect& quantity, double& integrated_quantity)
{
  compute_integral_quantity(quantity, integrated_quantity);
  integrated_quantity /= (probe_length_);
}

void IJK_One_Dimensional_Subproblem::compute_energy_from_temperature_interp()
{
  compute_integral_quantity_on_probe(temperature_interp_, energy_temperature_interp_);
}

void IJK_One_Dimensional_Subproblem::retrieve_previous_temperature_on_probe()
{
  if (!reconstruct_previous_probe_field_)
    temperature_previous_ = temperature_interp_;
  else
    {
      if (ref_ijk_ft_->get_tstep() == 0)
        {
          temperature_previous_.resize((*points_per_thermal_subproblem_));
          const double radius_ini = osculating_radius_;
          for (int i=0; i<(*points_per_thermal_subproblem_); i++)
            temperature_previous_[i] = delta_T_subcooled_overheated_ - delta_T_subcooled_overheated_
                                       * (radius_ini / osculating_radial_coordinates_[i])
                                       * (1 - erf( (*radial_coordinates_)[i] / (2 * sqrt((*alpha_) * ref_ijk_ft_->get_current_time())) ) );
        }
      else
        {
          const int count_index_ijk = is_in_map_index_ijk((*subproblem_to_ijk_indices_previous_), index_i_, index_j_, index_k_);
          if (count_index_ijk)
            {
              const int previous_rank = (*subproblem_to_ijk_indices_previous_).at(index_i_).at(index_j_).at(index_k_);
              temperature_previous_ = (*temperature_probes_previous_)[previous_rank];
            }
          else
            {
              const int neighbours_i[6] = NEIGHBOURS_I;
              const int neighbours_j[6] = NEIGHBOURS_J;
              const int neighbours_k[6] = NEIGHBOURS_K;
              int counter_mixed_neighbours = 0;
              //              double indicator_mixed_neighbours = 0;
              //              double colinearity_mixed_neighbours = 0;
              //              double velocity_mixed_neighbours = 0;
              double averaging_weight = 0.;
              DoubleVect temperature_previous, temperature_previous_options;
              temperature_previous.resize((*points_per_thermal_subproblem_));
              temperature_previous_options.resize((*points_per_thermal_subproblem_));
              for (int l=0; l<6; l++)
                {
                  const int ii = neighbours_i[l];
                  const int jj = neighbours_j[l];
                  const int kk = neighbours_k[l];
                  int index_i_neighbour_prev = index_i_ + neighbours_i[l];
                  int index_j_neighbour_prev = index_j_ + neighbours_j[l];
                  int index_k_neighbour_prev = index_k_ + neighbours_k[l];
                  const int count_index_ijk_neighbour = is_in_map_index_ijk((*subproblem_to_ijk_indices_previous_), index_i_neighbour_prev, index_j_neighbour_prev, index_k_neighbour_prev);
                  if (count_index_ijk_neighbour)
                    {
                      const int previous_rank = (*subproblem_to_ijk_indices_previous_).at(index_i_neighbour_prev).at(index_j_neighbour_prev).at(index_k_neighbour_prev);
                      const double indicator_prev = (*indicator_probes_previous_)[previous_rank];
                      const double indicator_prev_vapour = 1 - indicator_prev;
                      double best_indicator_prev = 0;
                      best_indicator_prev = (indicator_ > 0.5) ? indicator_prev_vapour : indicator_prev;
                      double velocity = 0.;
                      Vecteur3 velocities_xyz = (*velocities_probes_previous_)[previous_rank];
                      Vecteur3 normal_vector_compo_previous = (*normal_vector_compo_probes_previous_)[previous_rank];
                      const double colinearity = Vecteur3::produit_scalaire(normal_vector_compo_, normal_vector_compo_previous);
                      // int incr = 0;
                      if (ii)
                        {
                          velocity = velocities_xyz[0];
                          // incr = ii;
                        }
                      if (jj)
                        {
                          velocity = velocities_xyz[1];
                          // incr = jj;
                        }
                      if (kk)
                        {
                          velocity = velocities_xyz[2];
                          // incr = kk;
                        }
                      const double velocity_eval = abs(velocity) > 1e-10 ? abs(velocity) : 1.;
                      // if (signbit(velocity) != signbit(incr) || (abs(velocity) <= 1e-10))
                      {
                        retrieve_previous_temperature_on_probe_type(0, previous_rank,
                                                                    best_indicator_prev,
                                                                    colinearity,
                                                                    velocity_eval,
                                                                    temperature_previous,
                                                                    temperature_previous_options,
                                                                    averaging_weight);
                        counter_mixed_neighbours++;
                      }
                    }
                }
              if (counter_mixed_neighbours)
                {
                  if (averaging_weight < 1e-12)
                    {
                      if (debug_)
                        Cerr << "temperature_previous_mixed_neighbour: " << counter_mixed_neighbours << finl;
                      temperature_previous_ = temperature_previous;
                      temperature_previous_ *= (1 / counter_mixed_neighbours);
                    }
                  else
                    {
                      if (debug_)
                        Cerr << "temperature_previous_averaging_weight: " << averaging_weight << finl;
                      temperature_previous_ = temperature_previous_options;
                      temperature_previous_ *= (1 / averaging_weight);
                    }
                }
              else
                temperature_previous_ = temperature_interp_;
            }
        }
    }
}

void IJK_One_Dimensional_Subproblem::retrieve_previous_temperature_on_probe_type(const int computation_type,
                                                                                 const int& previous_rank,
                                                                                 const double& best_indicator_prev,
                                                                                 const double& colinearity,
                                                                                 const double& velocity_eval,
                                                                                 DoubleVect& temperature_previous,
                                                                                 DoubleVect& temperature_previous_options,
                                                                                 double& averaging_weight)
{
  DoubleVect temperature_previous_tmp = (*temperature_probes_previous_)[previous_rank];
  temperature_previous +=  temperature_previous_tmp;
  // temperature_previous_tmp *= best_indicator_prev;
  // temperature_previous_tmp *= velocity_eval;
  temperature_previous_tmp *= colinearity;
  temperature_previous_options += temperature_previous_tmp;
  //	averaging_weight += best_indicator_prev;
  //	averaging_weight += velocity_eval;
  averaging_weight += colinearity;
}

int IJK_One_Dimensional_Subproblem::is_in_map_index_ijk(const std::map<int, std::map<int, std::map<int, int>>>& subproblem_to_ijk_indices,
                                                        const int& index_i,
                                                        const int& index_j,
                                                        const int& index_k)
{
  const int count_index_i = (int) subproblem_to_ijk_indices.count(index_i);
  int count_index_j = 0;
  int count_index_k = 0;
  if (count_index_i)
    {
      count_index_j = (int) subproblem_to_ijk_indices.at(index_i).count(index_j);
      if (count_index_j)
        count_index_k = (int) subproblem_to_ijk_indices.at(index_i).at(index_j).count(index_k);
    }
  return (count_index_i && count_index_j && count_index_k);
}

void IJK_One_Dimensional_Subproblem::interpolate_temperature_on_probe()
{
  temperature_interp_.resize(*points_per_thermal_subproblem_);
  if (first_time_step_varying_probes_)
    ijk_interpolate_skip_unknown_points(*temperature_before_extrapolation_, coordinates_cartesian_compo_, temperature_interp_, INVALID_INTERP);
  else
    ijk_interpolate_skip_unknown_points(*temperature_, coordinates_cartesian_compo_, temperature_interp_, INVALID_INTERP);
  temperature_cell_ = (*temperature_)(index_i_, index_j_, index_k_);
  compute_energy_from_temperature_interp();
  retrieve_previous_temperature_on_probe();
}

void IJK_One_Dimensional_Subproblem::interpolate_temperature_gradient_on_probe()
{
  for (int dir = 0; dir < 3; dir++)
    {
      grad_T_elem_interp_[dir].resize(*points_per_thermal_subproblem_);
      ijk_interpolate_skip_unknown_points((*grad_T_elem_solver_)[dir], coordinates_cartesian_compo_, grad_T_elem_interp_[dir], INVALID_INTERP);
    }
}

void IJK_One_Dimensional_Subproblem::project_temperature_gradient_on_probes()
{
  normal_temperature_gradient_interp_.resize(*points_per_thermal_subproblem_);
  project_cartesian_onto_basis_vector(grad_T_elem_interp_[0], grad_T_elem_interp_[1], grad_T_elem_interp_[2],
                                      normal_vector_compo_, normal_temperature_gradient_interp_);

  tangential_temperature_gradient_first_.resize(*points_per_thermal_subproblem_);
  tangential_temperature_gradient_second_.resize(*points_per_thermal_subproblem_);
  project_cartesian_onto_basis_vector(grad_T_elem_interp_[0], grad_T_elem_interp_[1], grad_T_elem_interp_[2],
                                      first_tangential_vector_compo_, tangential_temperature_gradient_first_);
  project_cartesian_onto_basis_vector(grad_T_elem_interp_[0], grad_T_elem_interp_[1], grad_T_elem_interp_[2],
                                      second_tangential_vector_compo_, tangential_temperature_gradient_second_);

  azymuthal_temperature_gradient_.resize(*points_per_thermal_subproblem_);
  tangential_temperature_gradient_first_from_rising_dir_.resize(*points_per_thermal_subproblem_);
  project_cartesian_onto_basis_vector(grad_T_elem_interp_[0], grad_T_elem_interp_[1], grad_T_elem_interp_[2],
                                      azymuthal_vector_compo_, azymuthal_temperature_gradient_);
  project_cartesian_onto_basis_vector(grad_T_elem_interp_[0], grad_T_elem_interp_[1], grad_T_elem_interp_[2],
                                      first_tangential_vector_compo_from_rising_dir_, tangential_temperature_gradient_first_from_rising_dir_);
}

void IJK_One_Dimensional_Subproblem::interpolate_temperature_hessian_on_probe()
{
  for (int dir = 0; dir < 3; dir++)
    {
      hess_diag_T_elem_interp_[dir].resize(*points_per_thermal_subproblem_);
      hess_cross_T_elem_interp_[dir].resize(*points_per_thermal_subproblem_);
      hess_diag_T_elem_spherical_[dir].resize(*points_per_thermal_subproblem_);
      hess_cross_T_elem_spherical_[dir].resize(*points_per_thermal_subproblem_);
      hess_diag_T_elem_spherical_from_rising_[dir].resize(*points_per_thermal_subproblem_);
      hess_cross_T_elem_spherical_from_rising_[dir].resize(*points_per_thermal_subproblem_);
      ijk_interpolate_skip_unknown_points((*hess_diag_T_elem_)[dir], coordinates_cartesian_compo_,
                                          hess_diag_T_elem_interp_[dir], INVALID_INTERP);
      ijk_interpolate_skip_unknown_points((*hess_cross_T_elem_)[dir], coordinates_cartesian_compo_,
                                          hess_cross_T_elem_interp_[dir], INVALID_INTERP);
    }
  temperature_diffusion_hessian_cartesian_trace_ = hess_diag_T_elem_interp_[0];
  temperature_diffusion_hessian_cartesian_trace_ += hess_diag_T_elem_interp_[1];
  temperature_diffusion_hessian_cartesian_trace_ += hess_diag_T_elem_interp_[2];
}

void IJK_One_Dimensional_Subproblem::project_temperature_hessian_on_probes()
{
  compute_projection_matrix_cartesian_to_local_spherical();
  /*
   * A' = P^tAP
   */
  Matrice33 temperature_hessian_cartesian;
  Matrice33 temperature_hessian_spherical;
  Matrice33 temperature_hessian_spherical_from_rising;
  for (int i=0; i<*points_per_thermal_subproblem_; i++)
    {
      temperature_hessian_cartesian = Matrice33(hess_diag_T_elem_interp_[0](i), hess_cross_T_elem_interp_[2](i), hess_cross_T_elem_interp_[1](i),
                                                hess_cross_T_elem_interp_[2](i), hess_diag_T_elem_interp_[1](i), hess_cross_T_elem_interp_[0](i),
                                                hess_cross_T_elem_interp_[1](i), hess_cross_T_elem_interp_[0](i), hess_diag_T_elem_interp_[2](i));

      project_matrix_on_basis(projection_matrix_, inverse_projection_matrix_, temperature_hessian_cartesian, temperature_hessian_spherical);
      hess_diag_T_elem_spherical_[0](i) = temperature_hessian_spherical(0, 0);
      hess_diag_T_elem_spherical_[1](i) = temperature_hessian_spherical(1, 1);
      hess_diag_T_elem_spherical_[2](i) = temperature_hessian_spherical(2, 2);
      hess_cross_T_elem_spherical_[0](i) = temperature_hessian_spherical(1, 2);
      hess_cross_T_elem_spherical_[1](i) = temperature_hessian_spherical(0, 2);
      hess_cross_T_elem_spherical_[2](i) = temperature_hessian_spherical(0, 1);

      project_matrix_on_basis(projection_matrix_from_rising_, inverse_projection_matrix_from_rising_, temperature_hessian_cartesian, temperature_hessian_spherical_from_rising);
      hess_diag_T_elem_spherical_from_rising_[0](i) = temperature_hessian_spherical_from_rising(0, 0);
      hess_diag_T_elem_spherical_from_rising_[1](i) = temperature_hessian_spherical_from_rising(1, 1);
      hess_diag_T_elem_spherical_from_rising_[2](i) = temperature_hessian_spherical_from_rising(2, 2);
      hess_cross_T_elem_spherical_from_rising_[0](i) = temperature_hessian_spherical_from_rising(1, 2);
      hess_cross_T_elem_spherical_from_rising_[1](i) = temperature_hessian_spherical_from_rising(0, 2);
      hess_cross_T_elem_spherical_from_rising_[2](i) = temperature_hessian_spherical_from_rising(0, 1);
    }
}

void IJK_One_Dimensional_Subproblem::retrieve_temperature_diffusion_spherical_on_probes()
{
  temperature_diffusion_hessian_trace_ = hess_diag_T_elem_spherical_[0];
  temperature_diffusion_hessian_trace_ += hess_diag_T_elem_spherical_[1];
  temperature_diffusion_hessian_trace_ += hess_diag_T_elem_spherical_[2];
  radial_temperature_diffusion_ = normal_temperature_gradient_interp_;
  radial_temperature_diffusion_ *= osculating_radial_coordinates_inv_;
  radial_temperature_diffusion_ *= 2;
  radial_temperature_diffusion_ += hess_diag_T_elem_spherical_[0];
  tangential_temperature_diffusion_ = radial_temperature_diffusion_;
  tangential_temperature_diffusion_ *= -1;
  tangential_temperature_diffusion_ += temperature_diffusion_hessian_trace_;
}

void IJK_One_Dimensional_Subproblem::compute_projection_matrix_cartesian_to_local_spherical()
{
  /*
   * X = PX' <-> X'=P^tX
   * Easier to calculate the transpose/inverse P^t
   * Change of basis matrix send the coordinate to the basis of the LHS
   * X=P(beta->beta')X'
   * A' = P^tAP
   */
  inverse_projection_matrix_ = Matrice33(normal_vector_compo_[0], normal_vector_compo_[1], normal_vector_compo_[2],
                                         first_tangential_vector_compo_[0], first_tangential_vector_compo_[1], first_tangential_vector_compo_[2],
                                         second_tangential_vector_compo_[0], second_tangential_vector_compo_[1], second_tangential_vector_compo_[2]);
  Matrice33::inverse(inverse_projection_matrix_, projection_matrix_, 0);

  inverse_projection_matrix_from_rising_ = Matrice33(normal_vector_compo_[0], normal_vector_compo_[1], normal_vector_compo_[2],
                                                     first_tangential_vector_compo_from_rising_dir_[0], first_tangential_vector_compo_from_rising_dir_[1], first_tangential_vector_compo_from_rising_dir_[2],
                                                     azymuthal_vector_compo_[0], azymuthal_vector_compo_[1], azymuthal_vector_compo_[2]);
  Matrice33::inverse(inverse_projection_matrix_from_rising_, projection_matrix_from_rising_, 0);
}

void IJK_One_Dimensional_Subproblem::project_matrix_on_basis(const Matrice33& projection_matrix, const Matrice33& inverse_projection_matrix, const Matrice33& matrix, Matrice33& projected_matrix)
{
  Matrice33 matrix_ap_tmp = matrix;
  // AP
  Matrice33::produit_matriciel(matrix, projection_matrix, matrix_ap_tmp);
  // P^tAP
  Matrice33::produit_matriciel(inverse_projection_matrix, matrix_ap_tmp, projected_matrix);
}

void IJK_One_Dimensional_Subproblem::initialise_radial_convection_operator_local()
{
  if (!global_probes_characteristics_ ||
      first_time_step_varying_probes_ ||
      readjust_probe_length_from_vertices_ ||
      pre_initialise_thermal_subproblems_list_)
    (*finite_difference_assembler_).reinitialise_any_matrix_subproblem(radial_convection_matrix_base_,
                                                                       radial_first_order_operator_,
                                                                       sub_problem_index_,
                                                                       use_sparse_matrix_,
                                                                       first_indices_sparse_matrix_,
                                                                       operators_reinitialisation_,
                                                                       global_probes_characteristics_);
  operators_reinitialisation_ = 0;
}

void IJK_One_Dimensional_Subproblem::initialise_radial_diffusion_operator_local()
{
  if (!global_probes_characteristics_ ||
      first_time_step_varying_probes_ ||
      readjust_probe_length_from_vertices_ ||
      pre_initialise_thermal_subproblems_list_)
    (*finite_difference_assembler_).reinitialise_any_matrix_subproblem(radial_diffusion_matrix_base_,
                                                                       radial_second_order_operator_,
                                                                       sub_problem_index_,
                                                                       use_sparse_matrix_,
                                                                       first_indices_sparse_matrix_,
                                                                       operators_reinitialisation_,
                                                                       global_probes_characteristics_);
  operators_reinitialisation_ = 0;
}

void IJK_One_Dimensional_Subproblem::initialise_identity_operator_local()
{
  if ((!global_probes_characteristics_ || pre_initialise_thermal_subproblems_list_) && (*first_time_step_temporal_))
    (*finite_difference_assembler_).reinitialise_any_matrix_subproblem(identity_matrix_subproblems_,
                                                                       identity_matrix_explicit_implicit_,
                                                                       sub_problem_index_,
                                                                       use_sparse_matrix_,
                                                                       first_indices_sparse_matrix_,
                                                                       operators_reinitialisation_,
                                                                       global_probes_characteristics_);
  operators_reinitialisation_ = 0;
}

void IJK_One_Dimensional_Subproblem::compute_radial_convection_diffusion_operators()
{
  const double alpha_inv = - 1 / *alpha_;
  DoubleVect osculating_radial_coefficient = osculating_radial_coordinates_inv_;
  osculating_radial_coefficient *= 2;
  reinit_variable(radial_convection_prefactor_);
  // radial_convection_prefactor_.resize(*points_per_thermal_subproblem_);
  // interpolate_project_velocities_on_probes();
  if (source_terms_type_ != linear_diffusion
      && source_terms_type_ != spherical_diffusion
      && source_terms_type_ != spherical_diffusion_approx)
    {
      if (correct_radial_velocity_)
        radial_convection_prefactor_ = radial_velocity_corrected_;
      else
        radial_convection_prefactor_ = radial_velocity_;
      radial_convection_prefactor_ *= alpha_inv;
    }
  else
    Cerr << "Diffusion case : don't compute the radial convection pre-factor" << finl;
  if (source_terms_type_ != linear_diffusion)
    {
      if (source_terms_type_ == spherical_diffusion_approx)
        radial_convection_prefactor_ +=	2 / osculating_radius_;
      else
        radial_convection_prefactor_ +=	osculating_radial_coefficient;
    }
  const int boundary_conditions = 0;
  operators_reinitialisation_ = 1;
  initialise_radial_convection_operator_local();
  initialise_radial_diffusion_operator_local();
  initialise_identity_operator_local();
  (*finite_difference_assembler_).scale_matrix_subproblem_by_vector(radial_convection_matrix_base_,
                                                                    radial_convection_prefactor_,
                                                                    sub_problem_index_,
                                                                    boundary_conditions,
                                                                    use_sparse_matrix_,
                                                                    first_indices_sparse_matrix_);
}

void IJK_One_Dimensional_Subproblem::prepare_temporal_schemes()
{
  if (implicit_solver_from_previous_probe_field_)
    {
      if (sub_problem_index_ == 0)
        {
          (*radial_convection_matrix_base_) *= (- (*alpha_) * global_time_step_);
          (*radial_diffusion_matrix_base_) *= (- (*alpha_) * global_time_step_);
        }
    }

  if (*first_time_step_temporal_)
    {
      if (first_time_step_explicit_)
        (*finite_difference_assembler_).apply_euler_time_step(radial_convection_matrix_base_,
                                                              radial_diffusion_matrix_base_,
                                                              sub_problem_index_,
                                                              local_time_step_overall_,
                                                              (*alpha_));
      else
        (*finite_difference_assembler_).correct_sign_temporal_schemes_subproblems(radial_convection_matrix_base_,
                                                                                  radial_diffusion_matrix_base_,
                                                                                  sub_problem_index_,
                                                                                  local_time_step_overall_,
                                                                                  (*alpha_));
    }
}

void IJK_One_Dimensional_Subproblem::prepare_boundary_conditions(DoubleVect * thermal_subproblems_rhs_assembly,
                                                                 DoubleVect * thermal_subproblems_temperature_solution_ini,
                                                                 int& boundary_condition_interface,
                                                                 const double& interfacial_boundary_condition_value,
                                                                 const int& impose_boundary_condition_interface_from_simulation,
                                                                 int& boundary_condition_end,
                                                                 const double& end_boundary_condition_value,
                                                                 const int& impose_user_boundary_condition_end_value)
{
  interpolate_temperature_on_probe();
  interpolate_temperature_gradient_on_probe();
  project_temperature_gradient_on_probes();

  for (int i=0; i<temperature_interp_.size(); i++)
    if (fabs(temperature_interp_[i]) > INVALID_INTERP_TEST)
      Cerr << "Error in the temperature_interpolation" << temperature_interp_[i] << finl;

  /*
   * Written for Dirichlet only
   * TODO: Write for other B.Cs or vapour-liquid coupled system
   */
  if (!impose_boundary_condition_interface_from_simulation)
    interfacial_boundary_condition_value_ = interfacial_boundary_condition_value;
  else
    {
      switch (boundary_condition_interface)
        {
        case default_bc:
          interfacial_boundary_condition_value_ = temperature_interp_[0];
          break;
        case dirichlet:
          interfacial_boundary_condition_value_ = temperature_interp_[0];
          break;
        case neumann:
          interfacial_boundary_condition_value_ = normal_temperature_gradient_interp_[0];
          break;
        case flux_jump:
          interfacial_boundary_condition_value_ = normal_temperature_gradient_interp_[0];
          break;
        default:
          break;
        }
    }


  if (impose_user_boundary_condition_end_value)
    end_boundary_condition_value_ = end_boundary_condition_value;
  else
    {
      switch (boundary_condition_end)
        {
        case default_bc:
          end_boundary_condition_value_ = temperature_interp_[temperature_interp_.size() -1];
          break;
        case dirichlet:
          end_boundary_condition_value_ = temperature_interp_[temperature_interp_.size() -1];
          break;
        case neumann:
          end_boundary_condition_value_ = normal_temperature_gradient_interp_[temperature_interp_.size() -1];
          break;
        case flux_jump:
          end_boundary_condition_value_ = normal_temperature_gradient_interp_[temperature_interp_.size() -1];
          break;
        default:
          break;
        }
    }


  const int thermal_subproblems_rhs_size = (*thermal_subproblems_rhs_assembly_).size();
  if ((use_sparse_matrix_ || pre_initialise_thermal_subproblems_list_) && global_probes_characteristics_)
    start_index_ = (int) (sub_problem_index_ * (*points_per_thermal_subproblem_));
  else
    {
      start_index_ = thermal_subproblems_rhs_size;
      (*thermal_subproblems_rhs_assembly_).resize(thermal_subproblems_rhs_size + *points_per_thermal_subproblem_);
    }
  end_index_ = start_index_ + (*points_per_thermal_subproblem_);

  if (implicit_solver_from_previous_probe_field_)
    {
      /*
       * Some tests !
       */
      boundary_condition_end = implicit;
      end_boundary_condition_value_ = 0.;

      boundary_condition_end = neumann;
      end_boundary_condition_value_ = normal_temperature_gradient_interp_[(*points_per_thermal_subproblem_) - 1];
      normal_temperature_gradient_previous_.resize(temperature_previous_.size());
      (*finite_difference_assembler_).compute_operator(radial_first_order_operator_, temperature_previous_, normal_temperature_gradient_previous_);
      // end_boundary_condition_value_ = normal_temperature_gradient_previous_[(*points_per_thermal_subproblem_) - 1];
    }

  /*
   * Some tests...
   */
  if (*first_time_step_temporal_)
    {
      if (first_time_step_explicit_)
        if (!(pre_initialise_thermal_subproblems_list_ && global_probes_characteristics_))
          (*thermal_subproblems_temperature_solution_ini_).resize(thermal_subproblems_rhs_size + *points_per_thermal_subproblem_);

      temperature_ini_temporal_schemes_.resize(*points_per_thermal_subproblem_);
      if (boundary_condition_interface == dirichlet || boundary_condition_interface == default_bc)
        temperature_ini_temporal_schemes_[0] = interfacial_boundary_condition_value;
      else
        temperature_ini_temporal_schemes_[0] = temperature_interp_[0];

      double end_field_value_temporal=0.;
      if (boundary_condition_end == dirichlet || boundary_condition_end == default_bc)
        if (impose_user_boundary_condition_end_value)
          end_field_value_temporal = end_boundary_condition_value;
        else
          end_field_value_temporal = delta_T_subcooled_overheated_;
      else
        end_field_value_temporal = temperature_interp_[(*points_per_thermal_subproblem_) - 1];
      for (int i=1; i<(*points_per_thermal_subproblem_); i++)
        temperature_ini_temporal_schemes_[i] = end_field_value_temporal;
      // temperature_ini_temporal_schemes_[(*points_per_thermal_subproblem_) - 1] = delta_T_subcooled_overheated_;
    }
}

void IJK_One_Dimensional_Subproblem::compute_source_terms_impose_boundary_conditions(const int& boundary_condition_interface,
                                                                                     const double& interfacial_boundary_condition_value,
                                                                                     const int& impose_boundary_condition_interface_from_simulation,
                                                                                     const int& boundary_condition_end,
                                                                                     const double& end_boundary_condition_value,
                                                                                     const int& impose_user_boundary_condition_end_value)
{
  int boundary_condition_interface_local = boundary_condition_interface;
  int boundary_condition_end_local = boundary_condition_end;

  prepare_boundary_conditions(thermal_subproblems_rhs_assembly_,
                              thermal_subproblems_temperature_solution_ini_,
                              boundary_condition_interface_local,
                              interfacial_boundary_condition_value,
                              impose_boundary_condition_interface_from_simulation,
                              boundary_condition_end_local,
                              end_boundary_condition_value,
                              impose_user_boundary_condition_end_value);

  compute_source_terms();

  (*finite_difference_assembler_).impose_boundary_conditions_subproblem(thermal_subproblems_matrix_assembly_,
                                                                        thermal_subproblems_rhs_assembly_,
                                                                        rhs_assembly_,
                                                                        boundary_condition_interface_local,
                                                                        interfacial_boundary_condition_value_,
                                                                        boundary_condition_end_local,
                                                                        end_boundary_condition_value_,
                                                                        sub_problem_index_,
                                                                        dr_inv_,
                                                                        (*first_time_step_temporal_),
                                                                        first_time_step_explicit_,
                                                                        temperature_ini_temporal_schemes_,
                                                                        start_index_,
                                                                        first_indices_sparse_matrix_,
                                                                        use_sparse_matrix_);

  add_source_terms(boundary_condition_interface_local, boundary_condition_end_local);

  add_source_terms_temporal_tests(boundary_condition_interface_local, boundary_condition_end_local);

}

void IJK_One_Dimensional_Subproblem::compute_source_terms()
{
  if (!pure_thermal_diffusion_ || !avoid_post_processing_all_terms_ || compute_tangential_variables_)
    {
      //      if (source_terms_type_ == tangential_conv_2D_tangential_diffusion_3D ||
      //          source_terms_type_ == tangential_conv_3D_tangentual_diffusion_3D ||
      //          !avoid_post_processing_all_terms_)
      //        {
      interpolate_temperature_hessian_on_probe();
      project_temperature_hessian_on_probes();
      retrieve_temperature_diffusion_spherical_on_probes();
      //        }
    }
  if (compute_tangential_variables_)
    {
      /*
       * Compute at least the tangential convection
       */
      compute_tangential_convection_source_terms_first();
      switch (source_terms_type_)
        {
        case tangential_conv_2D:
          tangential_convection_source_terms_ = tangential_convection_source_terms_first_;
          source_terms_ = tangential_convection_source_terms_;
          break;
        case tangential_conv_3D:
          compute_tangential_convection_source_terms_second();
          tangential_convection_source_terms_ = tangential_convection_source_terms_first_;
          tangential_convection_source_terms_ += tangential_convection_source_terms_second_;
          source_terms_ = tangential_convection_source_terms_;
          break;
        case tangential_conv_2D_tangential_diffusion_3D:
          tangential_convection_source_terms_ = tangential_convection_source_terms_first_;
          source_terms_ = tangential_convection_source_terms_;
          compute_tangential_diffusion_source_terms();
          // Be careful to the sign
          source_terms_ += tangential_diffusion_source_terms_;
          break;
        case tangential_conv_3D_tangentual_diffusion_3D:
          compute_tangential_convection_source_terms_second();
          tangential_convection_source_terms_ = tangential_convection_source_terms_first_;
          tangential_convection_source_terms_ += tangential_convection_source_terms_second_;
          source_terms_ = tangential_convection_source_terms_;

          compute_tangential_diffusion_source_terms();
          // Be careful to the sign

          source_terms_ += tangential_diffusion_source_terms_;
          break;
        default:
          tangential_convection_source_terms_ = tangential_convection_source_terms_first_;
          source_terms_ = tangential_convection_source_terms_;
          break;
        }

      if (implicit_solver_from_previous_probe_field_)
        {
          source_terms_ *= (-1);
          source_terms_ *= (global_time_step_ * (*alpha_));
          source_terms_ += temperature_previous_;
          // rhs_assembly_ = source_terms_;
          // rhs_assembly_ += temperature_previous_;
        }

      if (*first_time_step_temporal_)
        {
          source_terms_ *= (-1);
          source_terms_ *= (local_time_step_overall_ * (*alpha_));
          if (first_time_step_explicit_)
            {
              rhs_assembly_ = source_terms_;
              rhs_assembly_[0] = 0.;
            }
          else
            source_terms_ += temperature_ini_temporal_schemes_;
        }
    }
}

void IJK_One_Dimensional_Subproblem::compute_tangential_convection_source_terms_first()
{
  tangential_convection_source_terms_first_ = (*tangential_temperature_gradient_first_solver_);
  if (correct_tangential_temperature_gradient_)
    correct_tangential_temperature_gradient(tangential_convection_source_terms_first_);
  tangential_convection_source_terms_first_ *= (*first_tangential_velocity_solver_);
  const double alpha_inv = 1 / (*alpha_);
  tangential_convection_source_terms_first_ *= alpha_inv;
}

void IJK_One_Dimensional_Subproblem::compute_tangential_convection_source_terms_second()
{
  tangential_convection_source_terms_second_ = (*tangential_temperature_gradient_second_solver_);
  if (correct_tangential_temperature_gradient_)
    correct_tangential_temperature_gradient(tangential_convection_source_terms_second_);
  tangential_convection_source_terms_second_ *= (*second_tangential_velocity_solver_);
  const double alpha_inv = 1 / (*alpha_);
  tangential_convection_source_terms_second_ *= alpha_inv;
}

void IJK_One_Dimensional_Subproblem::compute_tangential_diffusion_source_terms()
{
  tangential_diffusion_source_terms_ = tangential_temperature_diffusion_;
  tangential_diffusion_source_terms_ *= (-1);
  if (correct_tangential_temperature_hessian_)
    correct_tangential_temperature_hessian(tangential_diffusion_source_terms_);
}

void IJK_One_Dimensional_Subproblem::add_source_terms(const int& boundary_condition_interface, const int& boundary_condition_end)
{
  if (!(*first_time_step_temporal_))
    if (!pure_thermal_diffusion_)
      (*finite_difference_assembler_).add_source_terms(thermal_subproblems_rhs_assembly_,
                                                       rhs_assembly_,
                                                       source_terms_,
                                                       start_index_,
                                                       boundary_condition_interface,
                                                       boundary_condition_end);
}

void IJK_One_Dimensional_Subproblem::add_source_terms_temporal_tests(const int& boundary_condition_interface, const int& boundary_condition_end)
{
  if ((*first_time_step_temporal_))
    {
      if (first_time_step_explicit_)
        (*finite_difference_assembler_).add_source_terms(thermal_subproblems_temperature_solution_ini_,
                                                         rhs_assembly_,
                                                         temperature_ini_temporal_schemes_,
                                                         start_index_,
                                                         boundary_condition_interface,
                                                         boundary_condition_end);
      else
        (*finite_difference_assembler_).add_source_terms(thermal_subproblems_rhs_assembly_,
                                                         rhs_assembly_,
                                                         source_terms_,
                                                         start_index_,
                                                         boundary_condition_interface,
                                                         boundary_condition_end);
    }
}

void IJK_One_Dimensional_Subproblem::approximate_temperature_increment_material_derivative()
{
  approximate_partial_temperature_time_increment();
  approximate_temperature_material_derivatives();
}

void IJK_One_Dimensional_Subproblem::approximate_partial_temperature_time_increment()
{
  temperature_time_increment_.resize(*points_per_thermal_subproblem_);
  if (!avoid_post_processing_all_terms_)
    {
      // const double dt_iter = 0.;
      const double alpha_liq = *alpha_;
      for (int i=0; i<*points_per_thermal_subproblem_; i++)
        {
          temperature_time_increment_[i] = -(x_velocity_(i) * grad_T_elem_interp_[0](i) +
                                             y_velocity_(i) * grad_T_elem_interp_[1](i) +
                                             z_velocity_(i) * grad_T_elem_interp_[2](i)) +
                                           alpha_liq *
                                           (hess_diag_T_elem_interp_[0](i) +
                                            hess_diag_T_elem_interp_[1](i) +
                                            hess_diag_T_elem_interp_[2](i));
        }
    }
}

void IJK_One_Dimensional_Subproblem::approximate_temperature_material_derivatives()
{
  /*
   * Two frame reference calculations
   * Advected or not by the velocity tangentially
   * -> 4 cases
   */
  if (!avoid_post_processing_all_terms_)
    {
      approximate_temperature_material_derivatives(normal_vector_compo_,
                                                   first_tangential_vector_compo_,
                                                   second_tangential_vector_compo_,
                                                   radial_velocity_advected_frame_,
                                                   first_tangential_velocity_advected_frame_,
                                                   second_tangential_velocity_advected_frame_,
                                                   temperature_time_increment_,
                                                   convective_term_advected_frame_,
                                                   material_derivative_advected_frame_);
      approximate_temperature_material_derivatives(normal_vector_compo_,
                                                   first_tangential_vector_compo_,
                                                   second_tangential_vector_compo_,
                                                   radial_velocity_static_frame_,
                                                   first_tangential_velocity_static_frame_,
                                                   second_tangential_velocity_static_frame_,
                                                   temperature_time_increment_,
                                                   convective_term_static_frame_,
                                                   material_derivative_static_frame_);
      approximate_temperature_material_derivatives(normal_vector_compo_,
                                                   first_tangential_vector_compo_from_rising_dir_,
                                                   azymuthal_vector_compo_,
                                                   radial_velocity_advected_frame_,
                                                   first_tangential_velocity_advected_frame_,
                                                   second_tangential_velocity_advected_frame_,
                                                   temperature_time_increment_,
                                                   convective_term_advected_frame_rising_,
                                                   material_derivative_advected_frame_rising_);
      approximate_temperature_material_derivatives(normal_vector_compo_,
                                                   first_tangential_vector_compo_from_rising_dir_,
                                                   azymuthal_vector_compo_,
                                                   radial_velocity_static_frame_,
                                                   first_tangential_velocity_static_frame_,
                                                   second_tangential_velocity_static_frame_,
                                                   temperature_time_increment_,
                                                   convective_term_static_frame_rising_,
                                                   material_derivative_static_frame_rising_);
    }
}

void IJK_One_Dimensional_Subproblem::approximate_temperature_material_derivatives(const Vecteur3& normal_vector_compo,
                                                                                  const Vecteur3& first_tangential_vector_compo,
                                                                                  const Vecteur3& second_tangential_vector_compo,
                                                                                  const DoubleVect& radial_velocity_frame,
                                                                                  const DoubleVect& first_tangential_velocity_frame,
                                                                                  const DoubleVect& second_tangential_velocity_frame,
                                                                                  const DoubleVect& temperature_time_increment,
                                                                                  DoubleVect& convective_term,
                                                                                  DoubleVect& material_derivative)
{
  material_derivative.resize(*points_per_thermal_subproblem_);
  convective_term.resize(*points_per_thermal_subproblem_);
  for (int i=0; i<*points_per_thermal_subproblem_; i++)
    {
      Vecteur3 cartesian_velocity = {x_velocity_(i), y_velocity_(i), z_velocity_(i)};

      const double radial_velocity = radial_velocity_frame[i];
      const double first_tangential_velocity = first_tangential_velocity_frame[i];
      const double second_tangential_velocity = second_tangential_velocity_frame[i];
      Vecteur3 temperature_gradient = {grad_T_elem_interp_[0](i), grad_T_elem_interp_[1](i), grad_T_elem_interp_[2](i)};
      Vecteur3 radial_convection_frame = normal_vector_compo;
      radial_convection_frame *= -radial_velocity;
      Vecteur3 first_tangent_convection_frame = first_tangential_vector_compo;
      first_tangent_convection_frame *= -first_tangential_velocity;
      Vecteur3 second_tangent_convection_frame = second_tangential_vector_compo;
      second_tangent_convection_frame *= -second_tangential_velocity;
      radial_convection_frame += cartesian_velocity;
      first_tangent_convection_frame += cartesian_velocity;
      second_tangent_convection_frame += cartesian_velocity;
      convective_term[i] = Vecteur3::produit_scalaire(radial_convection_frame, temperature_gradient) +
                           Vecteur3::produit_scalaire(first_tangent_convection_frame, temperature_gradient) +
                           Vecteur3::produit_scalaire(second_tangent_convection_frame, temperature_gradient);
      material_derivative[i] = temperature_time_increment[i] + convective_term[i];
    }
}

void IJK_One_Dimensional_Subproblem::correct_tangential_temperature_gradient(DoubleVect& tangential_convection_source_terms)
{
  DoubleVect tangential_convection_source_terms_tmp = tangential_convection_source_terms;
  for (int i=0; i<tangential_convection_source_terms.size(); i++)
    tangential_convection_source_terms[i] -= tangential_convection_source_terms_tmp[0];
}

void IJK_One_Dimensional_Subproblem::correct_tangential_temperature_hessian(DoubleVect& tangential_diffusion_source_terms)
{
  DoubleVect tangential_diffusion_source_terms_tmp = tangential_diffusion_source_terms;
  for (int i=0; i<tangential_diffusion_source_terms_tmp.size(); i++)
    tangential_diffusion_source_terms_tmp[i] -= tangential_diffusion_source_terms_tmp[0];
}

void IJK_One_Dimensional_Subproblem::retrieve_radial_quantities()
{
  if (!reference_gfm_on_probes_)
    {
      retrieve_temperature_solution();
      compute_local_temperature_gradient_solution();
      compute_local_velocity_gradient();
      compute_local_pressure_gradient();
    }
  else
    {
      retrieve_variables_solution_gfm_on_probes();
    }
  compute_integral_quantities_solution();
  if (debug_probe_collision_ && disable_probe_because_collision_)
    (*probe_collision_debug_field_)(index_i_, index_j_, index_k_) = 1;
  is_updated_ = true;
}

void IJK_One_Dimensional_Subproblem::complete_tangential_source_terms_for_post_processings()
{
  if (compute_tangential_variables_)
    {
      switch (source_terms_type_)
        {
        case tangential_conv_2D:
          compute_tangential_convection_source_terms_second();
          compute_tangential_diffusion_source_terms();
          break;
        case tangential_conv_3D:
          compute_tangential_diffusion_source_terms();
          break;
        case tangential_conv_2D_tangential_diffusion_3D:
          compute_tangential_convection_source_terms_second();
          break;
        case tangential_conv_3D_tangentual_diffusion_3D:
          break;
        default:
          compute_tangential_convection_source_terms_first();
          compute_tangential_convection_source_terms_second();
          compute_tangential_diffusion_source_terms();
          break;
        }
    }
  else
    {
      tangential_convection_source_terms_first_.resize(*points_per_thermal_subproblem_);
      tangential_convection_source_terms_second_.resize(*points_per_thermal_subproblem_);
      tangential_diffusion_source_terms_.resize(*points_per_thermal_subproblem_);
    }
}

void IJK_One_Dimensional_Subproblem::compute_integral_quantities_solution()
{
  compute_integral_quantity_on_probe(temperature_solution_, energy_temperature_solution_);
  compute_integral_quantity_on_probe(normal_temperature_gradient_solution_, normal_temperature_gradient_solution_numerical_integral_);
  normal_temperature_gradient_solution_integral_exact_ = (temperature_solution_[(*points_per_thermal_subproblem_)-1] -
                                                          temperature_solution_[0]) / (probe_length_);
  compute_integral_quantity_on_probe(normal_temperature_double_derivative_solution_, normal_temperature_double_derivative_solution_numerical_integral_);
  normal_temperature_double_derivative_solution_integral_exact_ = (normal_temperature_gradient_solution_[(*points_per_thermal_subproblem_)-1] -
                                                                   normal_temperature_gradient_solution_[0]) / (probe_length_);
  compute_integral_quantity_on_probe(radial_scale_factor_solution_, radial_scale_factor_solution_integral_);
  compute_integral_quantity_on_probe(radial_convection_solution_, radial_convection_solution_integral_);

  if (compute_tangential_variables_)
    {
      compute_integral_quantity_on_probe(tangential_convection_source_terms_first_, tangential_convection_first_integral_);
      compute_integral_quantity_on_probe(tangential_convection_source_terms_second_, tangential_convection_second_integral_);
      compute_integral_quantity_on_probe(tangential_diffusion_source_terms_, tangential_diffusion_integral_);
      tangential_diffusion_integral_ *= (- 1);
      DoubleVect total_source_terms_local = tangential_convection_source_terms_first_;
      total_source_terms_local += tangential_convection_source_terms_second_;
      total_source_terms_local += tangential_diffusion_source_terms_;
      compute_integral_quantity_on_probe(total_source_terms_local, tangential_source_terms_integral_);
    }
  energy_increment_times_dt =  (normal_temperature_double_derivative_solution_integral_exact_ + radial_scale_factor_solution_integral_) -
                               (radial_convection_solution_integral_ + tangential_source_terms_integral_);
  // time_increment_from_energy_increment_ = (energy_temperature_solution_ - energy_temperature_interp_) / energy_increment_times_dt;
}

void IJK_One_Dimensional_Subproblem::retrieve_variables_solution_gfm_on_probes()
{
  /*
   * Put zeros in DoubleVect or the temperature interpolated
   * for easier post-processing with other simulations
   */
  initialise_empty_variables_for_post_processing();
  interpolate_temperature_on_probe();
  interpolate_temperature_gradient_on_probe();
  project_temperature_gradient_on_probes();
  interpolate_temperature_hessian_on_probe();
  project_temperature_hessian_on_probes();
  retrieve_temperature_diffusion_spherical_on_probes();
  copy_interpolations_on_solution_variables_for_post_processing();
  if (compute_normal_derivative_on_reference_probes_)
    {
      compute_local_temperature_gradient_solution();
      compute_local_velocity_gradient();
      compute_local_pressure_gradient();
    }
}

void IJK_One_Dimensional_Subproblem::copy_interpolations_on_solution_variables_for_post_processing()
{
  const double flux_coeff = ((*lambda_) * surface_);
  temperature_solution_ = temperature_interp_;
  normal_temperature_gradient_solution_ = normal_temperature_gradient_interp_;
  temperature_x_gradient_solution_ = grad_T_elem_interp_[0];
  temperature_y_gradient_solution_ = grad_T_elem_interp_[1];
  temperature_z_gradient_solution_ = grad_T_elem_interp_[2];
  thermal_flux_ = normal_temperature_gradient_solution_;
  const double grad_T_elem_gfm = (*eulerian_grad_T_interface_ns_)(index_i_, index_j_, index_k_);
  thermal_flux_[0] = grad_T_elem_gfm;
  thermal_flux_*= flux_coeff;
  thermal_flux_interp_gfm_ = thermal_flux_;
  radial_temperature_diffusion_solution_ = radial_temperature_diffusion_;
  const double sign_temp = signbit(*delta_temperature_) ? -1 : 1;
  thermal_flux_abs_ = thermal_flux_total_ * sign_temp;

  thermal_flux_gfm_ = grad_T_elem_gfm * flux_coeff;
  thermal_flux_raw_ = normal_temperature_gradient_interp_[0] * flux_coeff;
  thermal_flux_lrs_ = normal_temperature_gradient_solution_[0] * flux_coeff;

  const double sign_flux = signbit(thermal_flux_raw_) ? -1. : 1.;
  thermal_flux_max_raw_ = sign_flux * std::max(abs(thermal_flux_raw_), abs(thermal_flux_lrs_));
  thermal_flux_max_gfm_ = sign_flux * std::max(abs(thermal_flux_gfm_), abs(thermal_flux_lrs_));
  thermal_flux_max_ = sign_flux * std::max(std::max(abs(thermal_flux_gfm_), abs(thermal_flux_raw_)), abs(thermal_flux_lrs_));
}

void IJK_One_Dimensional_Subproblem::retrieve_temperature_solution()
{
  temperature_solution_.resize(rhs_assembly_.size());
  for (int i=start_index_; i<end_index_; i++)
    temperature_solution_[i - start_index_] = (*thermal_subproblems_temperature_solution_)[i];
}

void IJK_One_Dimensional_Subproblem::compute_local_temperature_gradient_solution()
{
  //	for (int dir=0; dir<3; dir++)
  //		temperature_gradient_solution_[dir].resize(temperature_solution_.size());
  normal_temperature_gradient_solution_.resize(temperature_solution_.size());
  (*finite_difference_assembler_).compute_operator(radial_first_order_operator_, temperature_solution_, normal_temperature_gradient_solution_);

  normal_temperature_double_derivative_solution_.resize(temperature_solution_.size());
  (*finite_difference_assembler_).compute_operator(radial_second_order_operator_, temperature_solution_, normal_temperature_double_derivative_solution_);

  reinit_variable(temperature_x_gradient_solution_);
  reinit_variable(temperature_y_gradient_solution_);
  reinit_variable(temperature_z_gradient_solution_);
  if (((source_terms_type_ == linear_diffusion
        || source_terms_type_ == spherical_diffusion
        || source_terms_type_ == spherical_diffusion_approx)
       && !compute_tangential_variables_) || use_normal_gradient_for_flux_corr_)
    {
      DoubleVect dummy_tangential_deriv;
      dummy_tangential_deriv.resize(*points_per_thermal_subproblem_);
      temperature_x_gradient_solution_.resize(normal_temperature_gradient_solution_.size());
      project_basis_vector_onto_cartesian_dir(0, dummy_tangential_deriv, dummy_tangential_deriv, normal_temperature_gradient_solution_,
                                              *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                              temperature_x_gradient_solution_);
      temperature_y_gradient_solution_.resize(normal_temperature_gradient_solution_.size());
      project_basis_vector_onto_cartesian_dir(1, dummy_tangential_deriv, dummy_tangential_deriv, normal_temperature_gradient_solution_,
                                              *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                              temperature_y_gradient_solution_);
      temperature_z_gradient_solution_.resize(normal_temperature_gradient_solution_.size());
      project_basis_vector_onto_cartesian_dir(2, dummy_tangential_deriv, dummy_tangential_deriv, normal_temperature_gradient_solution_,
                                              *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                              temperature_z_gradient_solution_);
    }
  else
    {
      temperature_x_gradient_solution_.resize(normal_temperature_gradient_solution_.size());
      project_basis_vector_onto_cartesian_dir(0, tangential_temperature_gradient_first_, tangential_temperature_gradient_second_, normal_temperature_gradient_solution_,
                                              *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                              temperature_x_gradient_solution_);
      temperature_y_gradient_solution_.resize(normal_temperature_gradient_solution_.size());
      project_basis_vector_onto_cartesian_dir(1, tangential_temperature_gradient_first_, tangential_temperature_gradient_second_, normal_temperature_gradient_solution_,
                                              *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                              temperature_y_gradient_solution_);
      temperature_z_gradient_solution_.resize(normal_temperature_gradient_solution_.size());
      project_basis_vector_onto_cartesian_dir(2, tangential_temperature_gradient_first_, tangential_temperature_gradient_second_, normal_temperature_gradient_solution_,
                                              *first_tangential_vector_compo_solver_, *second_tangential_vector_compo_solver_, normal_vector_compo_,
                                              temperature_z_gradient_solution_);
    }


  disable_probe_weak_gradient_local_ = 0;
  if (disable_probe_weak_gradient_)
    {
      const double interfacial_gradient_solution = abs(normal_temperature_gradient_solution_[0]);
      const double interfacial_gradient_interp = abs(normal_temperature_gradient_interp_[0]);
      const double interfacial_gradient_gfm = abs((*eulerian_grad_T_interface_ns_)(index_i_, index_j_, index_k_));
      if (disable_probe_weak_gradient_gfm_)
        {
          if (interfacial_gradient_solution < interfacial_gradient_gfm)
            disable_probe_weak_gradient_local_ = 1;
        }
      else
        {
          if (interfacial_gradient_solution < interfacial_gradient_interp)
            disable_probe_weak_gradient_local_ = 1;
        }
    }

  if ((source_terms_type_ == linear_diffusion
       || source_terms_type_ == spherical_diffusion
       || source_terms_type_ == spherical_diffusion_approx)
      && avoid_post_processing_all_terms_)
    initialise_empty_variables_for_post_processing();

  compute_radial_convection_scale_factor_solution();
  compute_radial_temperature_diffusion_solution();
  complete_tangential_source_terms_for_post_processings();

  if (disable_probe_weak_gradient_local_)
    thermal_flux_ = normal_temperature_gradient_interp_;
  else
    thermal_flux_ = normal_temperature_gradient_solution_;

  const double grad_T_elem_gfm = (*eulerian_grad_T_interface_ns_)(index_i_, index_j_, index_k_);
  if (reference_gfm_on_probes_)
    {
      normal_temperature_gradient_solution_[0] = grad_T_elem_gfm;
      thermal_flux_[0] = grad_T_elem_gfm;
    }
  thermal_flux_interp_gfm_ = normal_temperature_gradient_interp_;
  thermal_flux_interp_gfm_[0] = grad_T_elem_gfm;

  const double flux_coeff = ((*lambda_) * surface_);
  thermal_flux_*= flux_coeff;
  thermal_flux_interp_gfm_ *= flux_coeff;
  thermal_flux_total_ = thermal_flux_[0];
  const double sign_temp = signbit(*delta_temperature_) ? -1 : 1;
  thermal_flux_abs_ = thermal_flux_total_ * sign_temp;

  thermal_flux_gfm_ = grad_T_elem_gfm * flux_coeff;
  thermal_flux_raw_ = normal_temperature_gradient_interp_[0] * flux_coeff;
  thermal_flux_lrs_ = normal_temperature_gradient_solution_[0] * flux_coeff;

  const double sign_flux = signbit(thermal_flux_raw_) ? -1. : 1.;
  thermal_flux_max_raw_ = sign_flux * std::max(abs(thermal_flux_raw_), abs(thermal_flux_lrs_));
  thermal_flux_max_gfm_ = sign_flux * std::max(abs(thermal_flux_gfm_), abs(thermal_flux_lrs_));
  thermal_flux_max_ = sign_flux * std::max(std::max(abs(thermal_flux_gfm_), abs(thermal_flux_raw_)), abs(thermal_flux_lrs_));
}

void IJK_One_Dimensional_Subproblem::compute_radial_convection_scale_factor_solution()
{
  radial_scale_factor_interp_ = osculating_radial_coordinates_inv_;
  radial_scale_factor_interp_ *= 2;
  radial_scale_factor_interp_ *= normal_temperature_gradient_interp_;
  radial_scale_factor_solution_ = osculating_radial_coordinates_inv_;
  radial_scale_factor_solution_ *= 2;
  radial_scale_factor_solution_ *= normal_temperature_gradient_solution_;
  radial_convection_interp_ = normal_temperature_gradient_interp_;
  radial_convection_solution_ = normal_temperature_gradient_solution_;
  radial_convection_interp_ *= radial_velocity_;
  radial_convection_solution_ *= radial_velocity_;
  const double alpha_inv = 1 / (*alpha_);
  radial_convection_interp_ *= alpha_inv;
  radial_convection_solution_ *= alpha_inv;
}

void IJK_One_Dimensional_Subproblem::compute_radial_temperature_diffusion_solution()
{
  radial_temperature_diffusion_solution_ = normal_temperature_gradient_solution_;
  radial_temperature_diffusion_solution_ *= osculating_radial_coordinates_inv_;
  radial_temperature_diffusion_solution_ *= 2;
  radial_temperature_diffusion_solution_ += normal_temperature_double_derivative_solution_;
}

void IJK_One_Dimensional_Subproblem::initialise_empty_variables_for_post_processing()
{
  normal_temperature_gradient_interp_.resize(*points_per_thermal_subproblem_);
  normal_temperature_double_derivative_solution_.resize(*points_per_thermal_subproblem_);

  tangential_temperature_gradient_first_.resize(*points_per_thermal_subproblem_);
  tangential_temperature_gradient_second_.resize(*points_per_thermal_subproblem_);
  tangential_temperature_gradient_first_from_rising_dir_.resize(*points_per_thermal_subproblem_);
  azymuthal_temperature_gradient_.resize(*points_per_thermal_subproblem_);

  for (int dir = 0; dir < 3; dir++)
    {
      hess_diag_T_elem_interp_[dir].resize(*points_per_thermal_subproblem_);
      hess_cross_T_elem_interp_[dir].resize(*points_per_thermal_subproblem_);
      hess_diag_T_elem_spherical_[dir].resize(*points_per_thermal_subproblem_);
      hess_cross_T_elem_spherical_[dir].resize(*points_per_thermal_subproblem_);
      hess_diag_T_elem_spherical_from_rising_[dir].resize(*points_per_thermal_subproblem_);
      hess_cross_T_elem_spherical_from_rising_[dir].resize(*points_per_thermal_subproblem_);
    }

  temperature_diffusion_hessian_trace_.resize(*points_per_thermal_subproblem_);
  temperature_diffusion_hessian_cartesian_trace_.resize(*points_per_thermal_subproblem_);
  radial_temperature_diffusion_.resize(*points_per_thermal_subproblem_);
  tangential_temperature_diffusion_.resize(*points_per_thermal_subproblem_);
  radial_temperature_diffusion_solution_.resize(*points_per_thermal_subproblem_);

  radial_scale_factor_interp_.resize(*points_per_thermal_subproblem_);
  radial_scale_factor_solution_.resize(*points_per_thermal_subproblem_);

  radial_convection_interp_.resize(*points_per_thermal_subproblem_);
  radial_convection_solution_.resize(*points_per_thermal_subproblem_);

  tangential_convection_source_terms_first_.resize(*points_per_thermal_subproblem_);
  tangential_convection_source_terms_second_.resize(*points_per_thermal_subproblem_);
  tangential_diffusion_source_terms_.resize(*points_per_thermal_subproblem_);

  source_terms_.resize(*points_per_thermal_subproblem_);
  pressure_interp_.resize(*points_per_thermal_subproblem_);

  normal_velocity_normal_gradient_.resize(*points_per_thermal_subproblem_);
  first_tangential_velocity_normal_gradient_.resize(*points_per_thermal_subproblem_);
  second_tangential_velocity_normal_gradient_.resize(*points_per_thermal_subproblem_);
  first_tangential_velocity_normal_gradient_from_rising_dir_.resize(*points_per_thermal_subproblem_);
  azymuthal_velocity_normal_gradient_.resize(*points_per_thermal_subproblem_);

  shear_stress_.resize(*points_per_thermal_subproblem_);
  shear_stress_from_rising_dir_.resize(*points_per_thermal_subproblem_);
}

double IJK_One_Dimensional_Subproblem::get_interfacial_gradient_corrected() const
{
  return normal_temperature_gradient_solution_[0];
}

double IJK_One_Dimensional_Subproblem::get_interfacial_double_derivative_corrected() const
{
  return normal_temperature_double_derivative_solution_[0];
}

void IJK_One_Dimensional_Subproblem::compute_local_velocity_gradient()
{
  normal_velocity_normal_gradient_.resize(*points_per_thermal_subproblem_);
  first_tangential_velocity_normal_gradient_.resize(*points_per_thermal_subproblem_);
  second_tangential_velocity_normal_gradient_.resize(*points_per_thermal_subproblem_);
  azymuthal_velocity_normal_gradient_.resize(*points_per_thermal_subproblem_);
  first_tangential_velocity_normal_gradient_from_rising_dir_.resize(*points_per_thermal_subproblem_);

  (*finite_difference_assembler_).compute_operator(radial_first_order_operator_,
                                                   radial_velocity_,
                                                   normal_velocity_normal_gradient_);
  (*finite_difference_assembler_).compute_operator(radial_first_order_operator_,
                                                   first_tangential_velocity_,
                                                   first_tangential_velocity_normal_gradient_);
  (*finite_difference_assembler_).compute_operator(radial_first_order_operator_,
                                                   second_tangential_velocity_,
                                                   second_tangential_velocity_normal_gradient_);
  (*finite_difference_assembler_).compute_operator(radial_first_order_operator_,
                                                   first_tangential_velocity_from_rising_dir_,
                                                   first_tangential_velocity_normal_gradient_from_rising_dir_);
  (*finite_difference_assembler_).compute_operator(radial_first_order_operator_,
                                                   azymuthal_velocity_,
                                                   azymuthal_velocity_normal_gradient_);
  compute_local_shear_stress();
}

void IJK_One_Dimensional_Subproblem::compute_local_shear_stress()
{
  shear_stress_.resize(*points_per_thermal_subproblem_);
  shear_stress_from_rising_dir_.resize(*points_per_thermal_subproblem_);
  for (int i=0; i<(*points_per_thermal_subproblem_); i++)
    {
      Vecteur3 local_shear_stress = first_tangential_vector_compo_;
      Vecteur3 local_shear_stress_second = second_tangential_vector_compo_;
      local_shear_stress *= first_tangential_velocity_normal_gradient_(i);
      local_shear_stress_second *= second_tangential_velocity_normal_gradient_(i);
      local_shear_stress += local_shear_stress_second;
      shear_stress_(i) = local_shear_stress.length();

      local_shear_stress = first_tangential_vector_compo_from_rising_dir_;
      local_shear_stress_second = azymuthal_vector_compo_;
      local_shear_stress *= first_tangential_velocity_normal_gradient_from_rising_dir_(i);
      local_shear_stress_second *= azymuthal_velocity_normal_gradient_(i);
      local_shear_stress += local_shear_stress_second;
      shear_stress_from_rising_dir_(i) = local_shear_stress.length();
    }
  const double mu_liquid = ref_ijk_ft_->get_mu_liquid();
  shear_stress_ *= mu_liquid;
  shear_stress_from_rising_dir_ *= mu_liquid;
  velocity_shear_stress_ = shear_stress_[0];
  velocity_shear_force_ = velocity_shear_stress_ * surface_;
}

void IJK_One_Dimensional_Subproblem::compute_local_pressure_gradient()
{
  pressure_normal_gradient_.resize(*points_per_thermal_subproblem_);
  (*finite_difference_assembler_).compute_operator(radial_first_order_operator_,
                                                   pressure_interp_,
                                                   pressure_normal_gradient_);
  pressure_gradient_ = pressure_normal_gradient_[0];
}

double IJK_One_Dimensional_Subproblem::get_normal_velocity_normal_gradient() const
{
  return normal_velocity_normal_gradient_[0];
}

double IJK_One_Dimensional_Subproblem::get_tangential_velocity_normal_gradient() const
{
  return first_tangential_velocity_normal_gradient_[0];
}

double IJK_One_Dimensional_Subproblem::get_second_tangential_velocity_normal_gradient() const
{
  return second_tangential_velocity_normal_gradient_[0];
}

double IJK_One_Dimensional_Subproblem::get_azymuthal_velocity_normal_gradient() const
{
  return azymuthal_velocity_normal_gradient_[0];
}

double IJK_One_Dimensional_Subproblem::get_field_profile_at_point(const double& dist,
                                                                  const DoubleVect& field,
                                                                  const DoubleVect& field_weak_gradient,
                                                                  const IJK_Field_double& eulerian_field,
                                                                  const int temp_bool,
                                                                  const int weak_gradient_variable,
                                                                  const int interp_eulerian) const
{
  double field_value = INVALID_TEMPERATURE;
  const DoubleVect * field_ref = &field;
  if (disable_probe_weak_gradient_local_ && weak_gradient_variable)
    {
      if (temp_bool)
        return temperature_cell_;
      field_ref = &field_weak_gradient;
    }
  if (debug_ && temp_bool)
    {
      Cerr << "Thermal subproblem index: " << sub_problem_index_ << finl;
      Cerr << "Radial ini: " << (*radial_coordinates_)[0] << finl;
      Cerr << "Radial coordinate end: " << (*radial_coordinates_)[*points_per_thermal_subproblem_-1] << finl;
      Cerr << "Field ini: " << field[0] << finl;
      Cerr << "Field end: " << field[*points_per_thermal_subproblem_-1] << finl;
      Cerr << "Field interp ini: " << temperature_interp_[0] << finl;
      Cerr << "Field interp end: " << temperature_interp_[*points_per_thermal_subproblem_-1] << finl;
      Cerr << "End_boundary_condition_value: " << end_boundary_condition_value_ << finl;
      Cerr << "Curvature: " << curvature_ << finl;
      Cerr << "Osculating radius: " << osculating_radius_ << finl;
    }
  if (dist >= (*radial_coordinates_)[0] && dist <= (*radial_coordinates_)[*points_per_thermal_subproblem_-1])
    {
      /*
       * Dummy dichotomy and linear interpolation along the probe
       */
      int left_interval = 0;
      int right_interval = *points_per_thermal_subproblem_-1;
      find_interval(dist, left_interval, right_interval);
      const double field_interp = ((*field_ref)[right_interval] - (*field_ref)[left_interval])
                                  / ((*radial_coordinates_)[right_interval] - (*radial_coordinates_)[left_interval]) *
                                  (dist-(*radial_coordinates_)[left_interval]) + (*field_ref)[left_interval];
      field_value = field_interp;
    }
  else if (dist < (*radial_coordinates_)[0])
    {
      if (temp_bool)
        {
          const double interfacial_temperature_gradient_solution = get_interfacial_gradient_corrected();
          const double interfacial_temperature_double_derivative_solution = get_interfacial_double_derivative_corrected();
          switch(order_approx_temperature_ext_)
            {
            case 1:
              field_value = interfacial_temperature_gradient_solution * dist;
              break;
            case 2:
              field_value = (interfacial_temperature_gradient_solution * dist
                             + 0.5 * interfacial_temperature_double_derivative_solution * (dist * dist));
              break;
            default:
              field_value = interfacial_temperature_gradient_solution * dist;
            }
        }
      else
        field_value = field[0];
    }
  else
    {
      if(interp_eulerian)
        {
          DoubleTab coordinates_point;
          DoubleVect field_interp(1);
          coordinates_point.resize(1,3);
          Vecteur3 compo_xyz = normal_vector_compo_;
          compo_xyz *= dist;
          compo_xyz += facet_barycentre_;
          for (int c=0; c<3; c++)
            coordinates_point(0,c) = compo_xyz[c];
          ijk_interpolate_skip_unknown_points(eulerian_field, coordinates_point, field_interp, INVALID_INTERP);
          field_value = field_interp(0);
        }
      else
        field_value = field[*points_per_thermal_subproblem_ - 1];
    }
  return field_value;
}

double IJK_One_Dimensional_Subproblem::get_field_profile_at_point(const double& dist, const DoubleVect& field, const int temp_bool) const
{
  double field_value = INVALID_TEMPERATURE;// temperature_solution_;
  if (debug_ && ((dist < 0 && indicator_>0.5) || (dist > (*radial_coordinates_)[*points_per_thermal_subproblem_-1] && indicator_>0.5)))
    {
      Cerr << "Probe length: " << probe_length_ << finl;
      Cerr << "Distance d: " << dist << finl;
      Cerr << "Indicator I: " << indicator_ << finl;
    }
  if (debug_ && temp_bool)
    {
      Cerr << "Thermal subproblem index: " << sub_problem_index_ << finl;
      Cerr << "Radial 	 ini: " << (*radial_coordinates_)[0] << finl;
      Cerr << "Radial coordinate end: " << (*radial_coordinates_)[*points_per_thermal_subproblem_-1] << finl;
      Cerr << "Field ini: " << field[0] << finl;
      Cerr << "Field end: " << field[*points_per_thermal_subproblem_-1] << finl;
      Cerr << "Field interp ini: " << temperature_interp_[0] << finl;
      Cerr << "Field interp end: " << temperature_interp_[*points_per_thermal_subproblem_-1] << finl;
      Cerr << "End_boundary_condition_value: " << end_boundary_condition_value_ << finl;
      Cerr << "Curvature: " << curvature_ << finl;
      Cerr << "Osculating radius: " << osculating_radius_ << finl;
    }
  if (dist >= (*radial_coordinates_)[0] && dist <= (*radial_coordinates_)[*points_per_thermal_subproblem_-1])
    {
      /*
       * Dummy dichotomy and linear interpolation along the probe
       */
      int left_interval = 0;
      int right_interval = *points_per_thermal_subproblem_-1;
      find_interval(dist, left_interval, right_interval);
      const double field_interp = (field[right_interval] - field[left_interval])
                                  / ((*radial_coordinates_)[right_interval] - (*radial_coordinates_)[left_interval]) *
                                  (dist-(*radial_coordinates_)[left_interval]) + field[left_interval];
      field_value = field_interp;
    }
  else if (dist < (*radial_coordinates_)[0])
    {
      if (temp_bool)
        {
          if (short_probe_condition_)
            field_value = 0.;
          else
            {
              const double interfacial_temperature_gradient_solution = get_interfacial_gradient_corrected();
              const double interfacial_temperature_double_derivative_solution = get_interfacial_double_derivative_corrected();
              switch(order_approx_temperature_ext_)
                {
                case 1:
                  field_value = interfacial_temperature_gradient_solution * dist;
                  break;
                case 2:
                  field_value = (interfacial_temperature_gradient_solution * dist
                                 + 0.5 * interfacial_temperature_double_derivative_solution * (dist * dist));
                  break;
                default:
                  field_value = interfacial_temperature_gradient_solution * dist;
                }
            }
        }
      else
        field_value = 0.;
    }
  else
    {
      if (temp_bool)
        {
          if (short_probe_condition_)
            {
              if (temperature_probe_condition_)
                field_value = cell_temperature_;
              else
                field_value = delta_T_subcooled_overheated_;
            }
          else
            field_value = field[(*points_per_thermal_subproblem_) - 1];
        }
      else
        field_value = field[(*points_per_thermal_subproblem_) - 1];
    }
  return field_value;
}

double IJK_One_Dimensional_Subproblem::get_temperature_profile_at_point(const double& dist) const
{
  return get_field_profile_at_point(dist, temperature_solution_, temperature_interp_, *temperature_, 1, 1, interp_eulerian_);
  // return get_field_profile_at_point(dist, temperature_solution_, 1);
}

double IJK_One_Dimensional_Subproblem::get_velocity_component_at_point(const double& dist,
                                                                       const int& dir,
                                                                       const int& index_i,
                                                                       const int& index_j,
                                                                       const int& index_k) const
{
  double velocity = 0;
  if (use_velocity_cartesian_grid_)
    velocity = get_velocity_cartesian_grid_value(dist,
                                                 dir,
                                                 signbit(normal_vector_compo_[dir]),
                                                 index_i,
                                                 index_j,
                                                 index_k);
  else
    {
      switch(dir)
        {
        case 0:
          velocity = get_field_profile_at_point(dist, x_velocity_, x_velocity_, (*velocity_)[0] , 0, 0, interp_eulerian_);
          break;
        case 1:
          velocity = get_field_profile_at_point(dist, y_velocity_, y_velocity_, (*velocity_)[1] , 0, 0, interp_eulerian_);
          break;
        case 2:
          velocity = get_field_profile_at_point(dist, z_velocity_, z_velocity_, (*velocity_)[2] , 0, 0, interp_eulerian_);
          break;
        default:
          velocity = get_field_profile_at_point(dist, x_velocity_, x_velocity_, (*velocity_)[0] , 0, 0, interp_eulerian_);
          break;
        }
    }
  return velocity;
}

double IJK_One_Dimensional_Subproblem::get_velocity_cartesian_grid_value(const double& dist,
                                                                         const int& dir,
                                                                         const int& sign_dir,
                                                                         const int& index_i,
                                                                         const int& index_j,
                                                                         const int& index_k) const
{
  /*
   * Access directly the value ?
   */
  double velocity = 0.;
  const int test_index_i = index_i != INVALID_INDEX;
  const int test_index_j = index_j != INVALID_INDEX;
  const int test_index_k = index_k != INVALID_INDEX;
  int i = test_index_i ? index_i: index_i_;
  int j = test_index_j ? index_j: index_j_;
  int k = test_index_k ? index_k: index_k_;
  // const double delta = (dist - cell_centre_distance_) * normal_vector_compo_[dir];
  //  double incr_dir;
  //  const IJK_Grid_Geometry& geom= ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  //  double ddir = geom.get_constant_delta(dir);
  //  incr_dir = (int) delta / ddir;
  if (!(test_index_i && test_index_j && test_index_k))
    {
      switch(dir)
        {
        case 0:
          if (!sign_dir)
            i += 1;
          break;
        case 1:
          if (!sign_dir)
            j += 1;
          break;
        case 2:
          if (!sign_dir)
            k += 1;
          break;
        default:
          if (!sign_dir)
            i += 1;
          break;
        }
    }
  velocity = (*velocity_)[dir](i,j,k);

  /*
   * Use inteporlations
   */
  //  DoubleTab coordinates_point;
  //	DoubleVect field_interp(1);
  //	coordinates_point.resize(1,3);
  //  Vecteur3 bary_cell = ref_ijk_ft_->get_splitting_ns().get_coords_of_dof(index_i_, index_j_, index_k_, IJK_Splitting::ELEM);
  //  const IJK_Grid_Geometry& geom= ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  //  double ddir = geom.get_constant_delta(dir) / 2.;
  //	Vecteur3 compo_xyz = bary_cell;
  //	ddir = sign_dir ? - ddir: ddir;
  //	compo_xyz[dir] += bary_cell[dir] + ddir;
  //	for (int c=0; c<3; c++)
  //		coordinates_point(0,c) = compo_xyz[c];
  //	ijk_interpolate_skip_unknown_points((*velocity_)[dir], coordinates_point, field_interp, INVALID_INTERP);
  //	velocity = field_interp(0);

  return velocity;
}

double IJK_One_Dimensional_Subproblem::get_temperature_gradient_profile_at_point(const double& dist, const int& dir) const
{
  double temperature_gradient = 0;
  switch(dir)
    {
    case 0:
      temperature_gradient = get_field_profile_at_point(dist,
                                                        temperature_x_gradient_solution_,
                                                        grad_T_elem_interp_[0],
                                                        (*grad_T_elem_)[0],
                                                        0, 1,
                                                        interp_eulerian_);
      break;
    case 1:
      temperature_gradient = get_field_profile_at_point(dist,
                                                        temperature_y_gradient_solution_,
                                                        grad_T_elem_interp_[1],
                                                        (*grad_T_elem_)[1],
                                                        0, 1,
                                                        interp_eulerian_);
      break;
    case 2:
      temperature_gradient = get_field_profile_at_point(dist,
                                                        temperature_z_gradient_solution_,
                                                        grad_T_elem_interp_[2],
                                                        (*grad_T_elem_)[2],
                                                        0, 1,
                                                        interp_eulerian_);
      break;
    default:
      temperature_gradient = get_field_profile_at_point(dist,
                                                        temperature_x_gradient_solution_,
                                                        grad_T_elem_interp_[0],
                                                        (*grad_T_elem_)[0],
                                                        0, 1,
                                                        interp_eulerian_);
      break;
    }
  return temperature_gradient;
}

double IJK_One_Dimensional_Subproblem::get_temperature_times_velocity_profile_at_point(const double& dist,
                                                                                       const int& dir,
                                                                                       const int& l,
                                                                                       const int& index_i,
                                                                                       const int& index_j,
                                                                                       const int& index_k,
                                                                                       const int& temperature)
{
  double temperature_interp = get_field_profile_at_point(dist, temperature_solution_, temperature_interp_, *temperature_,
                                                         1, 1, interp_eulerian_);
  if (temperature)
    return temperature_interp;
  double velocity_interp = get_velocity_component_at_point(dist, dir, index_i, index_j, index_k);

  if (use_corrected_velocity_convection_)
    {
      const double radial_corr_ = radial_velocity_[0] - radial_velocity_corrected_[0];
      const double first_tangential_corr = (*first_tangential_velocity_not_corrected_)[0] - (*first_tangential_velocity_solver_)[0];
      const double second_tangential_corr = (*second_tangential_velocity_not_corrected_)[0] - (*second_tangential_velocity_solver_)[0];
      velocity_interp -= (radial_corr_ * normal_vector_compo_[dir]);
      velocity_interp -= (first_tangential_corr * (*first_tangential_vector_compo_solver_)[dir]);
      velocity_interp -= (second_tangential_corr * (*second_tangential_vector_compo_solver_)[dir]);
    }
  const double temperature_interp_conv = temperature_interp_conv_flux_[dir];
  // temperature_interp_conv_flux_[dir] = (temperature_interp_conv + temperature_interp);
  if ( l!=-1 )
    temperature_interp_conv_flux_[l] = (temperature_interp_conv + temperature_interp);
  return temperature_interp * velocity_interp;
}

double IJK_One_Dimensional_Subproblem::get_temperature_gradient_times_conductivity_profile_at_point(const double& dist, const int& dir) const
{
  double diffusive_flux = 0;
  diffusive_flux = get_temperature_gradient_profile_at_point(dist, dir);
  diffusive_flux *= (*lambda_);
  return diffusive_flux;
}

DoubleVect IJK_One_Dimensional_Subproblem::get_field_discrete_integral_velocity_weighting_at_point(const double& dist,
                                                                                                   const int& levels,
                                                                                                   const int& dir,
                                                                                                   const DoubleVect& field,
                                                                                                   const DoubleVect& field_weak_gradient,
                                                                                                   const IJK_Field_double& eulerian_field,
                                                                                                   const int temp_bool,
                                                                                                   const int weak_gradient_variable,
                                                                                                   const int vel,
                                                                                                   const int& l)
{
  const int nb_values = (int) pow(4., (double) levels);
  DoubleVect discrete_values(nb_values);
  double surface = get_discrete_surface_at_level(dir, levels);
  double value;
  double velocity;
  int value_counter = 0;
  if (levels==0)
    {
      velocity = get_velocity_weighting(dist, dir, vel);
      value = get_field_profile_at_point(dist, field, field_weak_gradient, eulerian_field, temp_bool, weak_gradient_variable, interp_eulerian_);
      discrete_values(0) = value * surface * velocity;
    }
  else
    {
      double dl1_ini = 0.;
      double dl2_ini = 0.;
      Vecteur3 point_coords_ini = {0., 0., 0.};
      get_field_discrete_value_recursive(1, levels + 1, dir, dist, vel, surface,
                                         field, field_weak_gradient, eulerian_field,
                                         temp_bool, weak_gradient_variable,
                                         dl1_ini, dl2_ini,
                                         point_coords_ini, discrete_values, value_counter);
    }
  if (l != -1 && !disable_relative_velocity_energy_balance_)
    {
      for (int c=0; c<discrete_values.size(); c++)
        temperature_interp_conv_flux_[l] += discrete_values[c];
    }
  return discrete_values;
}

void IJK_One_Dimensional_Subproblem::get_field_discrete_value_recursive(const int& ilevel, const int& max_level,
                                                                        const int& dir, const double& dist,
                                                                        const int& vel,
                                                                        const double& surface,
                                                                        const DoubleVect& field,
                                                                        const DoubleVect& field_weak_gradient,
                                                                        const IJK_Field_double& eulerian_field,
                                                                        const int temp_bool,
                                                                        const int weak_gradient_variable,
                                                                        const double dl1_parent,
                                                                        const double dl2_parent,
                                                                        Vecteur3& point_coords_parent,
                                                                        DoubleVect& discrete_values,
                                                                        int& value_counter) const
{
  if (ilevel != max_level)
    {
      const double neighbours_first_dir[4] = NEIGHBOURS_FIRST_DIR;
      const double neighbours_second_dir[4] = NEIGHBOURS_SECOND_DIR;
      for(int i=ilevel; i<max_level; i++)
        {
          if (i==ilevel)
            {
              for(int l=0; l<4; l++)
                {
                  const double first_dir = neighbours_first_dir[l];
                  const double second_dir = neighbours_second_dir[l];
                  double dl1;
                  double dl2;
                  Vecteur3 point_coords = {0., 0., 0.};
                  get_discrete_two_dimensional_spacing(dir, ilevel, first_dir, second_dir, dl1, dl2, point_coords);
                  dl1 += dl1_parent;
                  dl2 += dl2_parent;
                  point_coords += point_coords_parent;
                  get_field_discrete_value_recursive(i+1, max_level, dir, dist, vel, surface,
                                                     field, field_weak_gradient, eulerian_field,
                                                     temp_bool, weak_gradient_variable,
                                                     dl1, dl2, point_coords, discrete_values, value_counter);
                }
            }
        }
    }
  else
    {
      const double dist_increment = Vecteur3::produit_scalaire(point_coords_parent, normal_vector_compo_);
      const double dist_value = dist + dist_increment;
      // const double value = get_field_profile_at_point(dist_value, field, 0);
      const double value = get_field_profile_at_point(dist_value, field, field_weak_gradient, eulerian_field,
                                                      temp_bool, weak_gradient_variable, interp_eulerian_);
      const double velocity = get_velocity_weighting(dist, dir, vel);
      discrete_values(value_counter) = value * surface * velocity;
      value_counter++;
    }
}

double IJK_One_Dimensional_Subproblem::get_velocity_weighting(const double& dist, const int& dir, const int vel) const
{
  if (vel)
    return get_velocity_component_at_point(dist, dir);
  else
    return 1.;
}

DoubleVect IJK_One_Dimensional_Subproblem::get_field_discrete_integral_at_point(const double& dist,
                                                                                const int& levels,
                                                                                const int& dir,
                                                                                const DoubleVect& field,
                                                                                const DoubleVect& field_weak_gradient,
                                                                                const IJK_Field_double& eulerian_field,
                                                                                const int temp_bool,
                                                                                const int weak_gradient_variable)
{
  return get_field_discrete_integral_velocity_weighting_at_point(dist, levels, dir,
                                                                 field, field_weak_gradient,
                                                                 eulerian_field,
                                                                 temp_bool, weak_gradient_variable, 0);
}

DoubleVect IJK_One_Dimensional_Subproblem::get_field_times_velocity_discrete_integral_at_point(const double& dist,
                                                                                               const int& levels,
                                                                                               const int& dir,
                                                                                               const DoubleVect& field,
                                                                                               const DoubleVect& field_weak_gradient,
                                                                                               const IJK_Field_double& eulerian_field,
                                                                                               const int& l)
{
  return get_field_discrete_integral_velocity_weighting_at_point(dist, levels, dir,
                                                                 field, field_weak_gradient, eulerian_field,
                                                                 1, 1, 1, l);
}

DoubleVect IJK_One_Dimensional_Subproblem::get_temperature_profile_discrete_integral_at_point(const double& dist,
                                                                                              const int& levels,
                                                                                              const int& dir)
{
  return get_field_discrete_integral_at_point(dist, levels, dir, temperature_solution_, temperature_interp_, *temperature_, 1, 1);
}

DoubleVect IJK_One_Dimensional_Subproblem::get_temperature_times_velocity_profile_discrete_integral_at_point(const double& dist,
                                                                                                             const int& levels,
                                                                                                             const int& dir,
                                                                                                             const int& l)
{
  return get_field_times_velocity_discrete_integral_at_point(dist, levels, dir, temperature_solution_, temperature_interp_, *temperature_, l);
}

DoubleVect IJK_One_Dimensional_Subproblem::get_temperature_gradient_profile_discrete_integral_at_point(const double& dist,
                                                                                                       const int& levels,
                                                                                                       const int& dir)
{
  DoubleVect temperature_gradient;
  switch(dir)
    {
    case 0:
      temperature_gradient = get_field_discrete_integral_at_point(dist, levels, dir,
                                                                  temperature_x_gradient_solution_, grad_T_elem_interp_[0],
                                                                  (*grad_T_elem_)[0], 0, 1);
      break;
    case 1:
      temperature_gradient = get_field_discrete_integral_at_point(dist, levels, dir,
                                                                  temperature_y_gradient_solution_, grad_T_elem_interp_[1],
                                                                  (*grad_T_elem_)[1], 0, 1);
      break;
    case 2:
      temperature_gradient = get_field_discrete_integral_at_point(dist, levels, dir,
                                                                  temperature_z_gradient_solution_, grad_T_elem_interp_[2],
                                                                  (*grad_T_elem_)[2], 0, 1);
      break;
    default:
      temperature_gradient = get_field_discrete_integral_at_point(dist, levels, dir,
                                                                  temperature_x_gradient_solution_, grad_T_elem_interp_[0],
                                                                  (*grad_T_elem_)[0], 0, 1);
      break;
    }
  return temperature_gradient;
}

DoubleVect IJK_One_Dimensional_Subproblem::get_temperature_gradient_times_conductivity_profile_discrete_integral_at_point(const double& dist, const int& levels, const int& dir)
{
  DoubleVect diffusive_flux = get_temperature_gradient_profile_discrete_integral_at_point(dist, levels, dir);
  diffusive_flux *= (*lambda_);
  return diffusive_flux;
}

void IJK_One_Dimensional_Subproblem::find_interval(const double& dist, int& left_interval, int& right_interval) const
{
  int mid_interval = left_interval + (right_interval - left_interval) / 2;
  while ((right_interval - left_interval) != 1)
    {
      if (dist > (*radial_coordinates_)[mid_interval])
        left_interval = mid_interval;
      else
        right_interval = mid_interval;
      mid_interval = left_interval + (right_interval - left_interval) / 2;
    }
}

void IJK_One_Dimensional_Subproblem::get_discrete_two_dimensional_spacing(const int& dir, const int& levels,
                                                                          const double& first_dir, const double& second_dir,
                                                                          double& dl1, double& dl2,
                                                                          Vecteur3& point_coords) const
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  switch(dir)
    {
    case 0:
      dl1 = dy;
      dl2 = dz;
      point_coords[1] = dl1 * first_dir;
      point_coords[2] = dl2 * second_dir;
      break;
    case 1:
      dl1 = dx;
      dl2 = dz;
      point_coords[0] = dl1 * first_dir;
      point_coords[2] = dl2 * second_dir;
      break;
    case 2:
      dl1 = dx;
      dl2 = dy;
      point_coords[0] = dl1 * first_dir;
      point_coords[1] = dl2 * second_dir;
      break;
    default:
      dl1 = dx;
      dl2 = dy;
      point_coords[0] = dl1 * first_dir;
      point_coords[1] = dl2 * second_dir;
      break;
    }
  dl1 /= pow(2., (double) levels + 1.);
  dl2 /= pow(2., (double) levels + 1.);
  dl1 *= first_dir;
  dl2 *= second_dir;
  point_coords *= (1 / pow(2., (double) levels + 1.));
}

double IJK_One_Dimensional_Subproblem::get_discrete_surface_at_level(const int& dir, const int& level) const
{
  const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  double surface = 0.;
  switch(dir)
    {
    case 0:
      surface = (dy*dz);
      break;
    case 1:
      surface = (dx*dz);
      break;
    case 2:
      surface = (dx*dy);
      break;
    default:
      surface = (dx*dy);
      break;
    }
  surface /= pow(pow(2., (double) level), 2.);
  return surface;
}

void IJK_One_Dimensional_Subproblem::compute_temperature_integral_subproblem_probe()
{
  temperature_integral_ = compute_temperature_integral_subproblem(probe_length_);
}

double IJK_One_Dimensional_Subproblem::compute_temperature_integral_subproblem(const double& distance)
{
  double max_distance = distance;
  if (distance <= 0 || distance > probe_length_)
    max_distance = probe_length_;
  ArrOfDouble discrete_int_eval(*points_per_thermal_subproblem_ - 1);
  double integral_eval = 0.;
  const double radial_incr = max_distance / (*points_per_thermal_subproblem_ - 1);
  for (int i=0; i<(*points_per_thermal_subproblem_) - 1; i++)
    {
      discrete_int_eval(i) = get_temperature_profile_at_point(radial_incr * i);
      discrete_int_eval(i) += get_temperature_profile_at_point(radial_incr * (i + 1));
      discrete_int_eval(i) *= (radial_incr / 2.);
      integral_eval += discrete_int_eval(i);
    }
  integral_eval *= (1 / (radial_incr + 1e-20));
  return integral_eval;
}

void IJK_One_Dimensional_Subproblem::compute_bubble_related_quantities()
{
  nusselt_number_ = normal_temperature_gradient_solution_;
  nusselt_number_liquid_temperature_ = normal_temperature_gradient_solution_;
  nusselt_number_ *=  ((2* (*radius_from_volumes_per_bubble_)(compo_connex_)) / (*delta_temperature_));
  nusselt_number_liquid_temperature_ *=  ((2 * (*radius_from_volumes_per_bubble_)(compo_connex_)) / (*mean_liquid_temperature_));
  const double atan_theta_incr_ini = M_PI / 2;
  const double atan_incr_factor = -1;
  const double theta = (theta_sph_ * atan_incr_factor) + atan_theta_incr_ini;
  DoubleVect radius_squared = osculating_radial_coordinates_;
  radius_squared *= radius_squared;
  radius_squared *= (double) std::sin(theta);
  nusselt_number_integrand_ = nusselt_number_;
  nusselt_number_integrand_ *= radius_squared;
  nusselt_number_liquid_temperature_integrand_ = nusselt_number_liquid_temperature_;
  nusselt_number_liquid_temperature_integrand_ *= radius_squared;
}


void IJK_One_Dimensional_Subproblem::thermal_subresolution_outputs(SFichier& fic, const int rank, const Nom& local_quantities_thermal_probes_time_index_folder)
{
  post_process_interfacial_quantities(fic, rank);
  post_process_radial_quantities(rank, local_quantities_thermal_probes_time_index_folder);
}

void IJK_One_Dimensional_Subproblem::thermal_subresolution_outputs_parallel(const int rank, const Nom& local_quantities_thermal_probes_time_index_folder)
{
  post_process_radial_quantities(rank, local_quantities_thermal_probes_time_index_folder);
}

void IJK_One_Dimensional_Subproblem::retrieve_interfacial_quantities(const int rank,
                                                                     const int& itr,
                                                                     std::vector<std::string> key_results_int,
                                                                     std::vector<std::string> key_results_double,
                                                                     std::map<std::string, ArrOfInt>& results_probes_int,
                                                                     std::map<std::string, ArrOfDouble>& results_probes_double,
                                                                     const int& coord)
{
  const double last_time = ref_ijk_ft_->get_current_time() - ref_ijk_ft_->get_timestep();
  const int last_time_index = ref_ijk_ft_->get_tstep() + (*latastep_reprise_);
  std::vector<int> results_int =
  {
    last_time_index, rank, index_post_processing_, global_subproblem_index_, sub_problem_index_
  };

  std::vector<double> results_double =
  {
    last_time,
    normal_vector_compo_[0], normal_vector_compo_[1], normal_vector_compo_[2],
    first_tangential_vector_compo_[0], first_tangential_vector_compo_[1], first_tangential_vector_compo_[2],
    second_tangential_vector_compo_[0], second_tangential_vector_compo_[1], second_tangential_vector_compo_[2],
    first_tangential_vector_compo_from_rising_dir_[0], first_tangential_vector_compo_from_rising_dir_[1], first_tangential_vector_compo_from_rising_dir_[2],
    azymuthal_vector_compo_[0], azymuthal_vector_compo_[1], azymuthal_vector_compo_[2],
    r_sph_, theta_sph_, phi_sph_,
    temperature_interp_[coord], temperature_solution_[coord], temperature_previous_[coord],
    normal_temperature_gradient_interp_[coord], normal_temperature_gradient_solution_[coord],
    (*eulerian_grad_T_interface_ns_)(index_i_, index_j_, index_k_),
    hess_diag_T_elem_spherical_[0][coord], normal_temperature_double_derivative_solution_[coord],
    tangential_temperature_gradient_first_[coord],
    tangential_temperature_gradient_second_[coord],
    tangential_temperature_gradient_first_from_rising_dir_[coord],
    azymuthal_temperature_gradient_[coord],
    temperature_diffusion_hessian_cartesian_trace_[coord],
    temperature_diffusion_hessian_trace_[coord],
    radial_temperature_diffusion_[coord],
    radial_temperature_diffusion_solution_[coord],
    tangential_temperature_diffusion_[coord],
    radial_scale_factor_interp_[coord], radial_scale_factor_solution_[coord],
    radial_convection_interp_[coord], 	radial_convection_solution_[coord],
    tangential_convection_source_terms_first_[coord], tangential_convection_source_terms_second_[coord],
    surface_, thermal_flux_[coord],
    thermal_flux_gfm_, thermal_flux_raw_,
    thermal_flux_lrs_, thermal_flux_max_,
    (*lambda_), (*alpha_), (*prandtl_number_),
    nusselt_number_[coord], nusselt_number_liquid_temperature_[coord],
    nusselt_number_integrand_[coord], nusselt_number_liquid_temperature_integrand_[coord],
    velocity_shear_force_, velocity_shear_stress_,
    pressure_interp_[coord],
    x_velocity_[coord], y_velocity_[coord], z_velocity_[coord],
    radial_velocity_[coord], radial_velocity_corrected_[coord],
    radial_velocity_static_frame_[coord], radial_velocity_advected_frame_[coord],
    first_tangential_velocity_[coord], first_tangential_velocity_corrected_[coord],
    first_tangential_velocity_static_frame_[coord], first_tangential_velocity_advected_frame_[coord],
    second_tangential_velocity_[coord], second_tangential_velocity_corrected_[coord],
    second_tangential_velocity_static_frame_[coord], second_tangential_velocity_advected_frame_[coord],
    first_tangential_velocity_from_rising_dir_[coord], first_tangential_velocity_from_rising_dir_corrected_[coord],
    first_tangential_velocity_from_rising_dir_static_frame_[coord], first_tangential_velocity_from_rising_dir_advected_frame_[coord],
    azymuthal_velocity_[coord], azymuthal_velocity_corrected_[coord],
    azymuthal_velocity_static_frame_[coord], azymuthal_velocity_advected_frame_[coord],
    normal_velocity_normal_gradient_[coord],
    first_tangential_velocity_normal_gradient_[coord],
    second_tangential_velocity_normal_gradient_[coord],
    first_tangential_velocity_normal_gradient_from_rising_dir_[coord],
    azymuthal_velocity_normal_gradient_[coord],
    (*bubbles_surface_)(compo_connex_), (*bubbles_volume_)(compo_connex_),
    (*radius_from_surfaces_per_bubble_)(compo_connex_), (*radius_from_volumes_per_bubble_)(compo_connex_),
    (*delta_temperature_), (*mean_liquid_temperature_),
    (*bubbles_rising_vectors_per_bubble_)(compo_connex_, 0),
    (*bubbles_rising_vectors_per_bubble_)(compo_connex_, 1),
    (*bubbles_rising_vectors_per_bubble_)(compo_connex_, 2)
  };
  int i;
  assert(key_results_int.size() == results_int.size());
  int size_int = (int) key_results_int.size();
  for (i=0; i<size_int; i++)
    results_probes_int[key_results_int[i]](itr) = results_int[i];
  assert(key_results_double.size() == results_double.size());
  int size_double = (int) key_results_double.size();
  for (i=0; i<size_double; i++)
    results_probes_double[key_results_double[i]](itr) = results_double[i];

}

void IJK_One_Dimensional_Subproblem::post_process_interfacial_quantities(SFichier& fic, const int rank, const int& coord) //SFichier& fic)
{
  // if (Process::je_suis_maitre())
  {
    if (is_updated_)
      {
        const double last_time = ref_ijk_ft_->get_current_time() - ref_ijk_ft_->get_timestep();
        const int last_time_index = ref_ijk_ft_->get_tstep() + (*latastep_reprise_);
        fic << last_time_index << " ";
        fic << rank << " " << index_post_processing_ << " " << global_subproblem_index_ << " " << sub_problem_index_ << " ";
        fic << last_time << " ";
        fic << normal_vector_compo_[0] << " " << normal_vector_compo_[1] << " " << normal_vector_compo_[2] << " ";
        fic << first_tangential_vector_compo_[0] << " " << first_tangential_vector_compo_[1] << " " << first_tangential_vector_compo_[2] << " ";
        fic << second_tangential_vector_compo_[0] << " " << second_tangential_vector_compo_[1] << " " << second_tangential_vector_compo_[2] << " ";
        fic << first_tangential_vector_compo_from_rising_dir_[0] << " " << first_tangential_vector_compo_from_rising_dir_[1] << " " << first_tangential_vector_compo_from_rising_dir_[2] << " ";
        fic << azymuthal_vector_compo_[0] << " " << azymuthal_vector_compo_[1] << " " << azymuthal_vector_compo_[2] << " ";
        fic << r_sph_ << " " << theta_sph_ << " " << phi_sph_ << " ";
        fic << temperature_interp_[coord] << " " << temperature_solution_[coord] << " " << temperature_previous_[coord] << " ";
        fic << normal_temperature_gradient_interp_[coord] << " " << normal_temperature_gradient_solution_[coord] << " ";
        fic << (*eulerian_grad_T_interface_ns_)(index_i_, index_j_, index_k_) << " ";
        fic << hess_diag_T_elem_spherical_[0][coord] << " " << normal_temperature_double_derivative_solution_[coord] << " ";
        fic << tangential_temperature_gradient_first_[coord] << " ";
        fic << tangential_temperature_gradient_second_[coord] << " ";
        fic << tangential_temperature_gradient_first_from_rising_dir_[coord] << " ";
        fic << azymuthal_temperature_gradient_[coord] << " ";
        fic << temperature_diffusion_hessian_cartesian_trace_[coord] << " ";
        fic << temperature_diffusion_hessian_trace_[coord] << " ";
        fic << radial_temperature_diffusion_[coord] << " ";
        fic << radial_temperature_diffusion_solution_[coord] << " ";
        fic << tangential_temperature_diffusion_[coord] << " ";
        fic << radial_scale_factor_interp_[coord] << " " << radial_scale_factor_solution_[coord] << " ";
        fic << radial_convection_interp_[coord] << " " << radial_convection_solution_[coord] << " ";
        fic << tangential_convection_source_terms_first_[coord] << " " << tangential_convection_source_terms_second_[coord] << " ";
        fic << surface_ << " " << thermal_flux_[coord] << " ";
        fic << thermal_flux_gfm_ << " " << thermal_flux_raw_ << " ";
        fic << thermal_flux_lrs_ << " " << thermal_flux_max_ << " ";
        fic << *lambda_ << " " << *alpha_ << " " << *prandtl_number_ << " ";
        fic << nusselt_number_[coord] << " " << nusselt_number_liquid_temperature_[coord] << " ";
        fic << nusselt_number_integrand_[coord] << " " << nusselt_number_liquid_temperature_integrand_[coord] << " ";
        fic << velocity_shear_force_ << " " << velocity_shear_stress_ << " ";
        fic << pressure_interp_[coord] << " ";
        fic << x_velocity_[coord] << " " << y_velocity_[coord] << " " << z_velocity_[coord] << " ";
        fic << radial_velocity_[coord] << " " << radial_velocity_corrected_[coord] << " ";
        fic << radial_velocity_static_frame_[coord] << " " << radial_velocity_advected_frame_[coord] << " ";
        fic << first_tangential_velocity_[coord] << " " << first_tangential_velocity_corrected_[coord] << " ";
        fic << first_tangential_velocity_static_frame_[coord] << " " << first_tangential_velocity_advected_frame_[coord] << " ";
        fic << second_tangential_velocity_[coord] << " " << second_tangential_velocity_corrected_[coord] << " ";
        fic << second_tangential_velocity_static_frame_[coord] << " " << second_tangential_velocity_advected_frame_[coord] << " ";
        fic << first_tangential_velocity_from_rising_dir_[coord] << " " << first_tangential_velocity_from_rising_dir_corrected_[coord] << " ";
        fic << first_tangential_velocity_from_rising_dir_static_frame_[coord] << " " << first_tangential_velocity_from_rising_dir_advected_frame_[coord] << " ";
        fic << azymuthal_velocity_[coord] << " " << azymuthal_velocity_corrected_[coord] << " ";
        fic << azymuthal_velocity_static_frame_[coord] << " " << azymuthal_velocity_advected_frame_[coord] << " ";
        fic << normal_velocity_normal_gradient_[coord] << " ";
        fic << first_tangential_velocity_normal_gradient_[coord] << " ";
        fic << second_tangential_velocity_normal_gradient_[coord] << " ";
        fic << first_tangential_velocity_normal_gradient_from_rising_dir_[coord] << " ";
        fic << azymuthal_velocity_normal_gradient_[coord] << " ";
        fic << (*bubbles_surface_)(compo_connex_) << " " << (*bubbles_volume_)(compo_connex_) << " ";
        fic << (*radius_from_surfaces_per_bubble_)(compo_connex_) << " " << (*radius_from_volumes_per_bubble_)(compo_connex_) << " ";
        fic << (*delta_temperature_) << " " << (*mean_liquid_temperature_) << " ";
        fic << (*bubbles_rising_vectors_per_bubble_)(compo_connex_, 0) << " ";
        fic << (*bubbles_rising_vectors_per_bubble_)(compo_connex_, 1) << " ";
        fic << (*bubbles_rising_vectors_per_bubble_)(compo_connex_, 2) << " ";
        fic << finl;
      }
  }
}

void IJK_One_Dimensional_Subproblem::post_process_radial_quantities(const int rank, const Nom& local_quantities_thermal_probes_time_index_folder)
{
  // if (Process::je_suis_maitre())
  {
    if (is_updated_ && is_post_processed_local_)
      {
        const int reset = 1;
        const int max_digit = 8;
        const int last_time_index = (*latastep_reprise_) + ref_ijk_ft_->get_tstep();
        const int nb_digit_index_post_pro = index_post_processing_ < 1 ? 1 : (int) (log10(index_post_processing_) + 1);
        const int nb_digit_index_global = global_subproblem_index_ < 1 ? 1 : (int) (log10(global_subproblem_index_) + 1);
        const int nb_digit_tstep = last_time_index < 1 ? 1 : (int) (log10(last_time_index) + 1);
        const int max_digit_rank = 3;
        const int nb_digit_rank = rank < 1 ? 1 : (int) (log10(rank) + 1);
        Nom probe_name = Nom("_thermal_rank_") + Nom(std::string(max_digit_rank - nb_digit_rank, '0')) + Nom(rank) + Nom("_thermal_subproblem_") + Nom(std::string(max_digit - nb_digit_index_post_pro, '0'))
                         + Nom(index_post_processing_)
                         + Nom("_global_") + Nom(std::string(max_digit - nb_digit_index_global, '0')) + Nom(global_subproblem_index_)
                         + Nom("_radial_quantities_time_index_") +
                         + Nom(std::string(max_digit - nb_digit_tstep, '0')) + Nom(last_time_index) + Nom(".out");
        Nom probe_header = Nom("tstep\tthermal_rank\tpost_pro_index\tglobal_subproblem\tlocal_subproblem\ttime"
                               "\tnx\tny\tnz"
                               "\tt1x\tt1y\tt1z"
                               "\tt2x\tt2y\tt2z"
                               "\ts1x\ts1y\ts1z"
                               "\ts2x\ts2y\ts2z"
                               "\tr_sph\ttheta_sph\tphi_sph"
                               "\tradial_coord"
                               "\ttemperature_interp\ttemperature_sol\ttemperature_prev"
                               "\ttemperature_gradient\ttemperature_gradient_sol"
                               "\ttemperature_double_deriv\ttemperature_double_deriv_sol"
                               "\ttemperature_gradient_tangential\ttemperature_gradient_tangential2"
                               "\ttemperature_gradient_tangential_rise\ttemperature_gradient_azymuthal"
                               "\ttemperature_diffusion_hessian_cartesian_trace"
                               "\ttemperature_diffusion_hessian_trace"
                               "\tradial_temperature_diffusion"
                               "\tradial_temperature_diffusion_sol"
                               "\ttangential_temperature_diffusion"
                               "\tradial_scale_factor_interp\tradial_scale_factor_sol"
                               "\tradial_convection_interp\tradial_convection_sol"
                               "\ttangential_convection_first\ttangential_convection_second"
                               "\tsurface\tthermal_flux"
                               "\tlambda_liq\talpha_liq\tprandtl_liq"
                               "\tnusselt_number\tnusselt_number_liquid_temperature"
                               "\tnusselt_number_integrand\tnusselt_number_liquid_temperature_integrand"
                               "\tshear\tforce"
                               "\tpressure"
                               "\tu_x\tu_y\tu_z"
                               "\tu_r\tu_r_corr\tu_r_static\tu_r_advected"
                               "\tu_theta\tu_theta_corr\tu_theta_static\tu_theta_advected"
                               "\tu_theta2\tu_theta2_corr\tu_theta2_static\tu_theta2_advected"
                               "\tu_theta_rise\tu_theta_rise_corr\tu_theta_rise_static\tu_theta_rise_advected"
                               "\tu_phi\tu_phi_corr\tu_phi_static\tu_phi_advected"
                               "\tdu_r_dr\tdu_theta_dr\tdu_theta2_dr\tdu_theta_rise_dr\tdu_phi_dr"
                               "\ttotal_surface\ttotal_volume\tradius_from_surface\tradius_from_volume"
                               "\tdelta_temperature\tmean_liquid_temperature"
                               "\trising_dir_x\trising_dir_y\trising_dir_z");
        SFichier fic = Open_file_folder(local_quantities_thermal_probes_time_index_folder, probe_name, probe_header, reset);
        const double last_time = ref_ijk_ft_->get_current_time() - ref_ijk_ft_->get_timestep();
        for (int i=0; i<(*points_per_thermal_subproblem_); i++)
          {
            fic << last_time_index << " ";
            fic << rank << " " << index_post_processing_ << " " << global_subproblem_index_ << " " << sub_problem_index_ << " ";
            fic << last_time << " ";
            fic << normal_vector_compo_[0] << " " << normal_vector_compo_[1] << " " << normal_vector_compo_[2] << " ";
            fic << first_tangential_vector_compo_[0] << " " << first_tangential_vector_compo_[1] << " " << first_tangential_vector_compo_[2] << " ";
            fic << second_tangential_vector_compo_[0] << " " << second_tangential_vector_compo_[1] << " " << second_tangential_vector_compo_[2] << " ";
            fic << first_tangential_vector_compo_from_rising_dir_[0] << " " << first_tangential_vector_compo_from_rising_dir_[1] << " " << first_tangential_vector_compo_from_rising_dir_[2] << " ";
            fic << azymuthal_vector_compo_[0] << " " << azymuthal_vector_compo_[1] << " " << azymuthal_vector_compo_[2] << " ";
            fic << r_sph_ << " " << theta_sph_ << " " << phi_sph_ << " ";
            fic << (*radial_coordinates_)[i] << " ";
            fic << temperature_interp_[i] << " " << temperature_solution_[i] << " " << temperature_previous_[i] << " ";
            fic << normal_temperature_gradient_interp_[i] << " " << normal_temperature_gradient_solution_[i] << " ";
            fic << hess_diag_T_elem_spherical_[0][i] << " " <<  normal_temperature_double_derivative_solution_[i] << " ";
            fic << tangential_temperature_gradient_first_[i] << " ";
            fic << tangential_temperature_gradient_second_[i] << " ";
            fic << tangential_temperature_gradient_first_from_rising_dir_[i] << " ";
            fic << azymuthal_temperature_gradient_[i] << " ";
            fic << temperature_diffusion_hessian_cartesian_trace_[i] << " ";
            fic << temperature_diffusion_hessian_trace_[i] << " ";
            fic << radial_temperature_diffusion_[i] << " ";
            fic << radial_temperature_diffusion_solution_[i] << " ";
            fic << tangential_temperature_diffusion_[i] << " ";
            fic << radial_scale_factor_interp_[i] << " " << radial_scale_factor_solution_[i] << " ";
            fic << radial_convection_interp_[i] << " " << radial_convection_solution_[i] << " ";
            fic << tangential_convection_source_terms_first_[i] << " " << tangential_convection_source_terms_second_[i] << " ";
            fic << surface_ << " " << thermal_flux_[i] << " ";
            // fic <<  << " " <<  << " ";
            // fic <<  << " " <<  << " ";
            fic << *lambda_ << " " << *alpha_ << " " << *prandtl_number_ << " ";
            fic << nusselt_number_[i] << " " << nusselt_number_liquid_temperature_[i] << " ";
            fic << nusselt_number_integrand_[i] << " " << nusselt_number_liquid_temperature_integrand_[i] << " ";
            fic << shear_stress_[i] << " " << (shear_stress_[i] * surface_) << " ";
            fic << pressure_interp_[i] << " ";
            fic << x_velocity_[i] << " " << y_velocity_[i] << " " << z_velocity_[i] << " ";
            fic << radial_velocity_[i] << " " << radial_velocity_corrected_[i] << " ";
            fic << radial_velocity_static_frame_[i] << " " << radial_velocity_advected_frame_[i] << " ";
            fic << first_tangential_velocity_[i] << " " << first_tangential_velocity_corrected_[i] << " ";
            fic << first_tangential_velocity_static_frame_[i] << " " << first_tangential_velocity_advected_frame_[i] << " ";
            fic << second_tangential_velocity_[i] << " " << second_tangential_velocity_corrected_[i] << " ";
            fic << second_tangential_velocity_static_frame_[i] << " " << second_tangential_velocity_advected_frame_[i] << " ";
            fic << first_tangential_velocity_from_rising_dir_[i] << " " << first_tangential_velocity_from_rising_dir_corrected_[i] << " ";
            fic << first_tangential_velocity_from_rising_dir_static_frame_[i] << " " << first_tangential_velocity_from_rising_dir_advected_frame_[i] << " ";
            fic << azymuthal_velocity_[i] << " " << azymuthal_velocity_corrected_[i] << " ";
            fic << azymuthal_velocity_static_frame_[i] << " " << azymuthal_velocity_advected_frame_[i] << " ";
            fic << normal_velocity_normal_gradient_[i] << " ";
            fic << first_tangential_velocity_normal_gradient_[i] << " ";
            fic << second_tangential_velocity_normal_gradient_[i] << " ";
            fic << first_tangential_velocity_normal_gradient_from_rising_dir_[i] << " ";
            fic << azymuthal_velocity_normal_gradient_[i] << " ";
            fic << (*bubbles_surface_)(compo_connex_) << " " << (*bubbles_volume_)(compo_connex_) << " ";
            fic << (*radius_from_surfaces_per_bubble_)(compo_connex_) << " " << (*radius_from_volumes_per_bubble_)(compo_connex_) << " ";
            fic << (*delta_temperature_) << " " << (*mean_liquid_temperature_) << " ";
            fic << (*bubbles_rising_vectors_per_bubble_)(compo_connex_, 0) << " ";
            fic << (*bubbles_rising_vectors_per_bubble_)(compo_connex_, 1) << " ";
            fic << (*bubbles_rising_vectors_per_bubble_)(compo_connex_, 2) << " ";
            fic << finl;
          }
        fic.close();
      }
  }
}

double IJK_One_Dimensional_Subproblem::get_min_temperature() const
{
  double min_temperature_value=1e20;
  for (int i=0; i<temperature_solution_.size(); i++)
    min_temperature_value = std::min(min_temperature_value, temperature_solution_[i]);
  return min_temperature_value;
}

double IJK_One_Dimensional_Subproblem::get_max_temperature() const
{
  double max_temperature_value=-1e20;
  for (int i=0; i<temperature_solution_.size(); i++)
    max_temperature_value = std::max(max_temperature_value, temperature_solution_[i]);
  return max_temperature_value;
}

double IJK_One_Dimensional_Subproblem::get_min_temperature_domain_ends() const
{
  double min_temperature_value = std::min(temperature_solution_[0], temperature_solution_[temperature_solution_.size()-1]);
  return min_temperature_value;
}

double IJK_One_Dimensional_Subproblem::get_max_temperature_domain_ends() const
{
  double max_temperature_value = std::max(temperature_solution_[0], temperature_solution_[temperature_solution_.size()-1]);
  return max_temperature_value;
}

void IJK_One_Dimensional_Subproblem::complete_frame_of_reference_lrs_fluxes_eval()
{
  if (!disable_relative_velocity_energy_balance_ && !has_computed_lrs_flux_frame_of_ref_terms_)
    {
      const IJK_Grid_Geometry& geom = ref_ijk_ft_->get_splitting_ns().get_grid_geometry();
      const int face_dir[6] = FACES_DIR;
      const int flux_out[6] = FLUXES_OUT;
      const double rho_cp = ref_ijk_ft_->get_rho_l() * (*cp_liquid_);
      FixedVector<double, 6> convective_term_frame_of_ref = temperature_interp_conv_flux_;
      for (int l=0; l<6; l++)
        convective_term_frame_of_ref[l] *= bubble_rising_velocity_compo_[face_dir[l]];
      convective_term_frame_of_ref *= (-1) * rho_cp;
      double surf_face;
      for (int l=0; l<6; l++)
        {
          surf_face = 1.;
          if (pure_liquid_neighbours_[l])
            {
              for (int c = 0; c < 3; c++)
                if (c!=face_dir[l])
                  surf_face *= geom.get_constant_delta(c);
              double flux_val = convective_term_frame_of_ref[face_dir[l]] * flux_out[l] * surf_face;
              const double sign_temp = signbit(*delta_temperature_) ? -1. : 1.;
              flux_val *= sign_temp;
              sum_convective_flux_op_lrs_ += flux_val;
              if (sign_temp * flux_val >= 0)
                sum_convective_flux_op_leaving_lrs_ += flux_val;
              else
                sum_convective_flux_op_entering_lrs_ += flux_val;
            }
        }
      has_computed_lrs_flux_frame_of_ref_terms_ = true;
    }
}

//const double& IJK_One_Dimensional_Subproblem::get_sum_convective_diffusive_flux_op_value_lrs(const int flux_type)

void IJK_One_Dimensional_Subproblem::set_pure_flux_corrected(const double& flux_face, const int& l, const int flux_type)
{
  /*
   * Positive contributions for flux outward
   */
  const double sign_temp = signbit(*delta_temperature_) ? -1. : 1.;
  if (flux_type==0)
    {
      const double rho_cp = ref_ijk_ft_->get_rho_l() * (*cp_liquid_);
      const double rho_cp_flux = rho_cp * flux_face * sign_temp;
      convective_flux_op_lrs_[l] = rho_cp_flux;
      sum_convective_flux_op_lrs_ += rho_cp_flux;
      if (rho_cp_flux >= 0)
        sum_convective_flux_op_leaving_lrs_ += rho_cp_flux;
      else
        sum_convective_flux_op_entering_lrs_ += rho_cp_flux;
    }
  else
    {
      const double flux_face_sign = sign_temp * flux_face;
      diffusive_flux_op_lrs_[l] = flux_face_sign;
      sum_diffusive_flux_op_lrs_ += flux_face_sign;
      if (flux_face_sign >= 0)
        sum_diffusive_flux_op_leaving_lrs_ += flux_face_sign;
      else
        sum_diffusive_flux_op_entering_lrs_ += flux_face_sign;
    }
}

void IJK_One_Dimensional_Subproblem::compute_error_flux_interface()
{
  const double sign_temp = signbit(*delta_temperature_) ? -1. : 1.;

  double total_flux_error = 0.;
  total_flux_error = (- thermal_flux_abs_);
  total_flux_error += (- sum_convective_flux_op_lrs_);
  total_flux_error += sum_diffusive_flux_op_lrs_;

  double weight_tot = 0.;
  const int face_dir[6] = FACES_DIR;
  const int flux_out[6] = FLUXES_OUT;

  const int neighbours_i[6] = NEIGHBOURS_I;
  const int neighbours_j[6] = NEIGHBOURS_J;
  const int neighbours_k[6] = NEIGHBOURS_K;

  int counter_assert = 0;
  double weight = 0.;
  std::vector<int> mixed_neighbours;
  for (int l=0; l<6; l++)
    if (!pure_liquid_neighbours_[l] && !pure_vapour_neighbours_[l])
      {
        double normal_compo = normal_vector_compo_[face_dir[l]];
        if (signbit(flux_out[l]) == signbit(normal_compo))
          {
            const int ii = neighbours_i[l];
            const int jj = neighbours_j[l];
            const int kk = neighbours_k[l];
            const int isolated_mixed_neighbours = (*zero_liquid_neighbours_)(index_i_ + ii, index_j_ + jj, index_k_ + kk);
            if (!isolated_mixed_neighbours)
              {
                compute_weighting_coefficient(l, weight, fluxes_corrections_weighting_);
                weight_tot += abs(weight);
                mixed_neighbours.push_back(l);
                counter_assert++;
              }
            else
              {
                if (debug_)
                  Cerr << "The neighbour is isolated" << finl;
              }
          }
      }

  for (int l=0; l<6; l++)
    if (pure_liquid_neighbours_[l])
      {
        double normal_compo = normal_vector_compo_[face_dir[l]];
        if (signbit(flux_out[l]) == signbit(normal_compo))
          {
            compute_weighting_coefficient(l, weight, fluxes_corrections_weighting_);
            corrective_flux_current_[l] = (total_flux_error * abs(weight) * flux_out[l] * sign_temp);
            weight_tot += abs(weight);
            counter_assert++;
          }
      }

  if (debug_)
    Cerr << "counter_assert: " << counter_assert << finl;
  assert(counter_assert <= 3);
  if (counter_assert)
    {
      corrective_flux_current_ *= (1 / weight_tot);
      for (int m=0; m<(int) mixed_neighbours.size(); m++)
        {
          const int mixed_neighbour = mixed_neighbours[m];
          // const double normal_compo = abs(normal_vector_compo_[face_dir[mixed_neighbour]]);
          compute_weighting_coefficient(mixed_neighbour, weight, fluxes_corrections_weighting_);
          corrective_flux_to_neighbours_[mixed_neighbour] = (total_flux_error * abs(weight)) / weight_tot;
          // * flux_out[mixed_neighbour]
        }
    }
  else
    {
      Cerr << "Some fluxes contributions are not relocated !" << finl;
    }

  for (int l=0; l<6; l++)
    diffusive_flux_op_lrs_[l] = 0.;
  sum_diffusive_flux_op_lrs_ = 0.;
  sum_diffusive_flux_op_leaving_lrs_ = 0.;
  sum_diffusive_flux_op_entering_lrs_ = 0.;
}

void IJK_One_Dimensional_Subproblem::compute_weighting_coefficient(const int& l, double& weight, const int& weight_type)
{
  weight = 0.;
  const int face_dir[6] = FACES_DIR;
  const int dir = face_dir[l];
  if (weight_type == 0)
    weight = normal_vector_compo_[dir];
  else
    {
      const int flux_out[6] = FLUXES_OUT;
      const double kinematic_viscosity = ref_ijk_ft_->get_mu_liquid() / ref_ijk_ft_->get_rho_l();
      const double velocity = (*first_tangential_velocity_solver_)[0] * (*first_tangential_vector_compo_solver_)[dir]
                              + (*second_tangential_velocity_solver_)[0] * (*second_tangential_vector_compo_solver_)[dir];
      Vecteur3 first_vel = (*first_tangential_vector_compo_solver_);
      first_vel *= (*first_tangential_velocity_solver_)[0];
      Vecteur3 second_vel = (*second_tangential_vector_compo_solver_);
      second_vel *= (*second_tangential_velocity_solver_)[0];
      Vecteur3 vel_vect = first_vel;
      vel_vect += second_vel;
      const double vel_vect_norm = vel_vect.length();
      if (vel_vect_norm > 1e-12)
        vel_vect *= (1 / vel_vect_norm);
      const int vel_sign = signbit(velocity);
      const int vel_effect = (vel_sign == signbit(flux_out[l])) ? 1. : 0.;
      if (weight_type == 1)
        weight = kinematic_viscosity * vel_effect * vel_vect[dir];
      else
        weight = (*alpha_) * normal_vector_compo_[dir] + kinematic_viscosity * vel_effect * vel_vect[dir];
    }
}

void IJK_One_Dimensional_Subproblem::compare_flux_interface(std::vector<double>& radial_flux_error)
{
  const int flux_out[6] = FLUXES_OUT;
  sum_convective_diffusive_flux_op_lrs_ = sum_diffusive_flux_op_lrs_ - sum_convective_flux_op_lrs_;
  radial_flux_error_lrs_ = sum_convective_diffusive_flux_op_lrs_ - thermal_flux_abs_;
  double weight_tot = 0.;
  for (int l=0; l<3; l++)
    weight_tot += abs(normal_vector_compo_[l]);
  const int face_dir[6] = FACES_DIR;
  std::vector<double> thermal_flux_neighbour;
  for (int l=0; l<6; l++)
    {
      // FIXME: There are more than 3 faces sometimes !!!!!
      if (pure_liquid_neighbours_[l])
        {
          const double weight = abs(normal_vector_compo_[face_dir[l]]);
          //          const double flux_sum = convective_flux_op_lrs_[l] + diffusive_flux_op_lrs_[l];
          //          const double weight = abs(flux_sum);
          //          weight_tot += weight;
          const double radial_flux_contrib = (weight * radial_flux_error_lrs_ ) * flux_out[l];
          thermal_flux_neighbour.push_back(- thermal_flux_dir_[face_dir[l]] * flux_out[l]);
          radial_flux_error.push_back(radial_flux_contrib);
        }
    }
  for (int l=0; l<(int) radial_flux_error.size(); l++)
    {
      radial_flux_error[l] /= weight_tot;
      radial_flux_error[l] += thermal_flux_neighbour[l];
    }
  // sum_convective_diffusive_flux_op_lrs_ = thermal_flux_[0];
  // TODO: TMP for post-processing
  //  sum_convective_flux_op_lrs_ = 0.;
  //  sum_diffusive_flux_op_lrs_ = sum_convective_diffusive_flux_op_lrs_;
}

void IJK_One_Dimensional_Subproblem::dispatch_interfacial_heat_flux_correction(FixedVector<IJK_Field_double,3>& interfacial_heat_flux_dispatched,
                                                                               FixedVector<ArrOfInt, 4>& ijk_indices_out,
                                                                               ArrOfDouble& thermal_flux_out,
                                                                               FixedVector<IJK_Field_double,3>& interfacial_heat_flux_current)
{
  if (!has_computed_liquid_neighbours_)
    compute_pure_liquid_neighbours();

  const int ni = ref_ijk_ft_->itfce().I().ni();
  const int nj = ref_ijk_ft_->itfce().I().nj();
  const int nk = ref_ijk_ft_->itfce().I().nk();

  const IJK_Splitting& splitting_ns = ref_ijk_ft_->itfce().I().get_splitting();
  const int offset_i = splitting_ns.get_offset_local(0);
  const int offset_j = splitting_ns.get_offset_local(1);
  const int offset_k = splitting_ns.get_offset_local(2);

  const IJK_Grid_Geometry& geometry = splitting_ns.get_grid_geometry();
  const int ni_tot = geometry.get_nb_elem_tot(0);
  const int nj_tot = geometry.get_nb_elem_tot(1);
  const int nk_tot = geometry.get_nb_elem_tot(2);

  const int neighbours_i[6] = NEIGHBOURS_I;
  const int neighbours_j[6] = NEIGHBOURS_J;
  const int neighbours_k[6] = NEIGHBOURS_K;

  const int face_dir[6] = FACES_DIR;
  // const int flux_out[6] = FLUXES_OUT;

  int index_i_neighbour_global, index_j_neighbour_global, index_k_neighbour_global;
  int index_i_procs, index_j_procs, index_k_procs;
  for (int l=0; l<6; l++)
    {
      interfacial_heat_flux_current[face_dir[l]](index_i_,index_j_,index_k_) += corrective_flux_current_[l];
      const double flux_corr = corrective_flux_to_neighbours_[l];
      if (abs(flux_corr) > 1e-16)
        {
          const int ii = neighbours_i[l];
          const int jj = neighbours_j[l];
          const int kk = neighbours_k[l];
          const int i = index_i_ + ii;
          const int j = index_j_ + jj;
          const int k = index_k_ + kk;
          index_i_neighbour_global = compute_periodic_index((i + offset_i), ni_tot);
          index_j_neighbour_global = compute_periodic_index((j + offset_j), nj_tot);
          index_k_neighbour_global = compute_periodic_index((k + offset_k), nk_tot);
          index_i_procs = compute_periodic_index(i, ni);
          index_j_procs = compute_periodic_index(j, nj);
          index_k_procs = compute_periodic_index(k, nk);
          if (index_i_procs == i
              && index_j_procs == j
              && index_k_procs == k)
            {
              interfacial_heat_flux_dispatched[face_dir[l]](i,j,k) += flux_corr;
            }
          else
            {
              ijk_indices_out[0].append_array(index_i_neighbour_global);
              ijk_indices_out[1].append_array(index_j_neighbour_global);
              ijk_indices_out[2].append_array(index_k_neighbour_global);
              ijk_indices_out[3].append_array(face_dir[l]);
              thermal_flux_out.append_array(flux_corr);
            }
        }
    }
}

void IJK_One_Dimensional_Subproblem::dispatch_interfacial_heat_flux(FixedVector<IJK_Field_double,3>& interfacial_heat_flux_dispatched,
                                                                    FixedVector<ArrOfInt, 3>& ijk_indices_out,
                                                                    FixedVector<ArrOfDouble, 3>& thermal_flux_out)
{
  if (!has_computed_liquid_neighbours_)
    compute_pure_liquid_neighbours();

  const int ni = ref_ijk_ft_->itfce().I().ni();
  const int nj = ref_ijk_ft_->itfce().I().nj();
  const int nk = ref_ijk_ft_->itfce().I().nk();

  const int neighbours_i[6] = NEIGHBOURS_I;
  const int neighbours_j[6] = NEIGHBOURS_J;
  const int neighbours_k[6] = NEIGHBOURS_K;
  bool is_all_mix = true;
  double weight_tot = 0.;
  for (int l=0; l<3; l++)
    weight_tot += abs(normal_vector_compo_[l]);
  const int face_dir[6] = FACES_DIR;
  const int flux_out[6] = FLUXES_OUT;
  std::vector<int> mixed_neighbours;
  for (int l=0; l<6; l++)
    {
      if (pure_liquid_neighbours_[l])
        is_all_mix = false;
      else if (!pure_vapour_neighbours_[l])
        {
          if (signbit(flux_out[l]) == signbit(normal_vector_compo_[face_dir[l]]))
            mixed_neighbours.push_back(l);
          // weight_tot += normal_vector_compo_[face_dir[l]];
        }
    }
  if (is_all_mix)
    for (int l=0; l<(int) mixed_neighbours.size(); l++)
      for (int m=0; m<(int) mixed_neighbours.size(); m++)
        if (m!=l)
          {
            const int mixed_neighbour = mixed_neighbours[l];
            const int neighbour_dir = mixed_neighbours[m];
            const int ii = neighbours_i[mixed_neighbour];
            const int jj = neighbours_j[mixed_neighbour];
            const int kk = neighbours_k[mixed_neighbour];
            const int i = index_i_ + ii;
            const int j = index_j_ + jj;
            const int k = index_k_ + kk;
            const double weight = abs(normal_vector_compo_[face_dir[mixed_neighbour]]) * abs(normal_vector_compo_[face_dir[neighbour_dir]]);
            const double heat_flux_dispatch = thermal_flux_total_ * weight / (weight_tot * weight_tot);
            if ((i==ni || i==-1) || (j==nj || j==-1) || (k==nk || k==-1))
              {
                ijk_indices_out[0].append_array(i);
                ijk_indices_out[1].append_array(j);
                ijk_indices_out[2].append_array(k);
                thermal_flux_out[face_dir[neighbour_dir]].append_array(heat_flux_dispatch);
                for (int n=0; n<(int) mixed_neighbours.size(); n++)
                  if (n != m)
                    {
                      const int neighbour_dir_tmp = mixed_neighbours[n];
                      thermal_flux_out[face_dir[neighbour_dir_tmp]].append_array(0.);
                    }
              }
            else
              interfacial_heat_flux_dispatched[face_dir[neighbour_dir]](i, j, k) += heat_flux_dispatch;
          }
}

void IJK_One_Dimensional_Subproblem::add_interfacial_heat_flux_neighbours_correction(FixedVector<IJK_Field_double,3>& interfacial_heat_flux_dispatched,
                                                                                     FixedVector<IJK_Field_double,3>& interfacial_heat_flux_current)
{
  const int isolated_mixed_cell = (*zero_liquid_neighbours_)(index_i_, index_j_, index_k_);
  if (!isolated_mixed_cell)
    {
      const double sign_temp = signbit(*delta_temperature_) ? -1. : 1.;
      const int flux_out[6] = FLUXES_OUT;

      Vecteur3 flux_error = {0.,0.,0.};
      for (int c=0; c<3; c++)
        flux_error[c] = interfacial_heat_flux_dispatched[c](index_i_, index_j_, index_k_);

      const int face_dir[6] = FACES_DIR;
      double weight_dir_tot;
      double weight = 0.;
      std::vector<int> pure_faces;
      std::vector<double> weight_dir;
      for (int c=0; c<3; c++)
        {
          weight_dir_tot = 0.;
          pure_faces.clear();
          weight_dir.clear();
          for (int l=0; l<6; l++)
            {
              double normal_compo = normal_vector_compo_[face_dir[l]];
              compute_weighting_coefficient(l, weight, fluxes_corrections_weighting_);
              if (signbit(flux_out[l]) == signbit(normal_compo))
                if (pure_liquid_neighbours_[l] && face_dir[l] != c)
                  {
                    // normal_compo = abs(normal_compo);
                    weight = abs(weight);
                    pure_faces.push_back(l);
                    weight_dir.push_back(weight);
                    weight_dir_tot += weight;
                  }
            }
          if (weight_dir.size() == 0)
            for (int l=0; l<6; l++)
              {
                double normal_compo = normal_vector_compo_[face_dir[l]];
                compute_weighting_coefficient(l, weight, fluxes_corrections_weighting_);
                if (signbit(flux_out[l]) == signbit(normal_compo))
                  if (pure_liquid_neighbours_[l])
                    {
                      // normal_compo = abs(normal_compo);
                      weight = abs(weight);
                      pure_faces.push_back(l);
                      weight_dir.push_back(weight);
                      weight_dir_tot += weight;
                    }
              }
          if (weight_dir.size() == 0)
            if (debug_)
              for (int l=0; l<6; l++)
                {
                  Cerr << "Neighbour face index l: " << l << finl;
                  Cerr << "Pure liquid neighbour:" << (int) pure_liquid_neighbours_[l] << finl;
                  Cerr << "Normal vector compo:" << normal_vector_compo_[face_dir[l]] << finl;
                  Cerr << "Flux value:" << flux_error[c] << finl;
                }
          assert(weight_dir.size() != 0);
          for (int m=0; m<(int) pure_faces.size(); m++)
            {
              const int pure_face = pure_faces[m];
              corrective_flux_from_neighbours_[pure_face] += ((sign_temp * flux_out[pure_face] * flux_error[c])
                                                              * (weight_dir[m] / weight_dir_tot));
            }
        }
      for (int l=0; l<6; l++)
        interfacial_heat_flux_current[face_dir[l]](index_i_,index_j_,index_k_) += corrective_flux_from_neighbours_[l];
    }
}

void IJK_One_Dimensional_Subproblem::add_interfacial_heat_flux_neighbours(FixedVector<IJK_Field_double,3>& interfacial_heat_flux_dispatched)
{
  for (int c=0; c<3; c++)
    thermal_flux_dir_[c] = interfacial_heat_flux_dispatched[c](index_i_, index_j_, index_k_);
}

void IJK_One_Dimensional_Subproblem::compute_pure_liquid_neighbours()
{
  if (!has_computed_liquid_neighbours_)
    {
      const int neighbours_i[6] = NEIGHBOURS_I;
      const int neighbours_j[6] = NEIGHBOURS_J;
      const int neighbours_k[6] = NEIGHBOURS_K;
      for (int l=0; l<6; l++)
        {
          const int ii = neighbours_i[l];
          const int jj = neighbours_j[l];
          const int kk = neighbours_k[l];
          const double indic_neighbour = ref_ijk_ft_->itfce().I()(index_i_+ii, index_j_+jj, index_k_+kk);
          if (fabs(indic_neighbour) > LIQUID_INDICATOR_TEST)
            {
              pure_liquid_neighbours_[l] = 1;
              pure_vapour_neighbours_[l] = 0;
            }
          else
            {
              pure_liquid_neighbours_[l] = 0;
              if (fabs(indic_neighbour) < VAPOUR_INDICATOR_TEST)
                pure_vapour_neighbours_[l] = 1;
              else
                pure_vapour_neighbours_[l] = 0;
            }
        }
      has_computed_liquid_neighbours_ = true;
    }
}

void IJK_One_Dimensional_Subproblem::locate_pure_mixed_neighbours_without_pure_liquid_faces()
{
  int liquid_faces = 0;
  for (int l=0; l<6; l++)
    if (pure_liquid_neighbours_[l])
      liquid_faces++;
  if (!liquid_faces)
    (*zero_liquid_neighbours_)(index_i_, index_j_, index_k_) = 1;
}

void IJK_One_Dimensional_Subproblem::compare_fluxes_thermal_subproblems(const FixedVector<IJK_Field_double, 3>& convective_diffusive_fluxes_raw,
                                                                        const int flux_type,
                                                                        const int inv_sign)
{
  FixedVector<double, 6>* convective_diffusive_flux_op_value = nullptr;
  FixedVector<double, 6>* convective_diffusive_flux_op_value_vap = nullptr;
  FixedVector<double, 6>* convective_diffusive_flux_op_value_mixed = nullptr;
  FixedVector<double, 6>* convective_diffusive_flux_op_value_normal_contrib = nullptr;
  FixedVector<double, 6>* convective_diffusive_flux_op_value_leaving = nullptr;
  FixedVector<double, 6>* convective_diffusive_flux_op_value_entering = nullptr;

  double * sum_convective_diffusive_flux_op_value = nullptr;
  double * sum_convective_diffusive_flux_op_value_vap = nullptr;
  double * sum_convective_diffusive_flux_op_value_mixed = nullptr;
  double * sum_convective_diffusive_flux_op_value_normal_contrib = nullptr;
  double * sum_convective_diffusive_flux_op_value_leaving = nullptr;
  double * sum_convective_diffusive_flux_op_value_entering = nullptr;

  switch(flux_type)
    {
    case 0:
      convective_diffusive_flux_op_value = &convective_flux_op_value_;
      convective_diffusive_flux_op_value_vap = &convective_flux_op_value_vap_;
      convective_diffusive_flux_op_value_mixed = &convective_flux_op_value_mixed_;
      convective_diffusive_flux_op_value_normal_contrib = &convective_flux_op_value_normal_contrib_;
      convective_diffusive_flux_op_value_leaving = &convective_flux_op_leaving_value_;
      convective_diffusive_flux_op_value_entering = &convective_flux_op_entering_value_;

      sum_convective_diffusive_flux_op_value = &sum_convective_flux_op_value_;
      sum_convective_diffusive_flux_op_value_vap = &sum_convective_flux_op_value_vap_;
      sum_convective_diffusive_flux_op_value_mixed = &sum_convective_flux_op_value_mixed_;
      sum_convective_diffusive_flux_op_value_normal_contrib = &sum_convective_flux_op_value_normal_contrib_;
      sum_convective_diffusive_flux_op_value_leaving = &sum_convective_flux_op_leaving_value_;
      sum_convective_diffusive_flux_op_value_entering = &sum_convective_flux_op_entering_value_;
      break;
    case 1:
      convective_diffusive_flux_op_value = &diffusive_flux_op_value_;
      convective_diffusive_flux_op_value_vap = &diffusive_flux_op_value_vap_;
      convective_diffusive_flux_op_value_mixed = &diffusive_flux_op_value_mixed_;
      convective_diffusive_flux_op_value_normal_contrib = &diffusive_flux_op_value_normal_contrib_;
      convective_diffusive_flux_op_value_leaving = &diffusive_flux_op_leaving_value_;
      convective_diffusive_flux_op_value_entering = &diffusive_flux_op_entering_value_;

      sum_convective_diffusive_flux_op_value = &sum_diffusive_flux_op_value_;
      sum_convective_diffusive_flux_op_value_vap = &sum_diffusive_flux_op_value_vap_;
      sum_convective_diffusive_flux_op_value_mixed = &sum_diffusive_flux_op_value_mixed_;
      sum_convective_diffusive_flux_op_value_normal_contrib = &sum_diffusive_flux_op_value_normal_contrib_;
      sum_convective_diffusive_flux_op_value_leaving = &sum_diffusive_flux_op_leaving_value_;
      sum_convective_diffusive_flux_op_value_entering = &sum_diffusive_flux_op_entering_value_;
      break;
    }

  (*sum_convective_diffusive_flux_op_value) = 0.;
  (*sum_convective_diffusive_flux_op_value_vap) = 0.;
  (*sum_convective_diffusive_flux_op_value_mixed) = 0.;
  (*sum_convective_diffusive_flux_op_value_normal_contrib) = 0.;
  (*sum_convective_diffusive_flux_op_value_leaving) = 0.;
  (*sum_convective_diffusive_flux_op_value_entering) = 0.;

  if (!has_computed_liquid_neighbours_)
    compute_pure_liquid_neighbours();

  const int flux_out[6] = FLUXES_OUT;
  // const int neighbours_ijk_sign[6] = NEIGHBOURS_SIGN;
  const int neighbours_faces_i[6] = NEIGHBOURS_FACES_I;
  const int neighbours_faces_j[6] = NEIGHBOURS_FACES_J;
  const int neighbours_faces_k[6] = NEIGHBOURS_FACES_K;
  const int face_dir[6] = FACES_DIR;
  for (int l=0; l<6; l++)
    {
      (*convective_diffusive_flux_op_value)[l] = 0.;
      (*convective_diffusive_flux_op_value_vap)[l] = 0.;
      (*convective_diffusive_flux_op_value_mixed)[l] = 0.;
      (*convective_diffusive_flux_op_value_normal_contrib)[l] = 0.;
      (*convective_diffusive_flux_op_value_leaving)[l] = 0.;
      (*convective_diffusive_flux_op_value_entering)[l] = 0.;

      const int ii_f = neighbours_faces_i[l];
      const int jj_f = neighbours_faces_j[l];
      const int kk_f = neighbours_faces_k[l];

      double flux_val = convective_diffusive_fluxes_raw[face_dir[l]](index_i_ + ii_f,
                                                                     index_j_ + jj_f,
                                                                     index_k_ + kk_f);
      if (inv_sign)
        flux_val = -flux_val;
      flux_val *= flux_out[l];

      // flux_val = neighbours_ijk_sign[l] ? -flux_val: flux_val;

      // Keep normals out
      // const double sign_conv = flux_type ? 1.: -1.;
      // flux_val *= sign_conv;

      const double sign_temp = signbit(*delta_temperature_) ? -1 : 1;
      flux_val *= sign_temp;

      /*
       * TODO: Count only positive contributions !
       */
      if (pure_liquid_neighbours_[l])
        {
          // const bool sign_check_bool = (signbit(normal_vector_compo_[face_dir[l]]) == signbit(flux_out[l]));
          // const int sign_check = sign_check_bool ? 1: -1;
          // flux_val *= sign_check;

          // if (!sign_check_bool && flux_type && debug_)
          //   {
          //     Cerr << "The diffusive flux is negative" << finl;
          //      Cerr << "normal_vector_compo: " << normal_vector_compo_[0] << "; ";
          //      Cerr << normal_vector_compo_[1] << "; " << normal_vector_compo_[2] << "; ";
          //      Cerr << "l: " << l << finl;
          //      Cerr << "flux_val: " << flux_val << finl;
          //  }

          (*convective_diffusive_flux_op_value)[l] = flux_val;
          (*sum_convective_diffusive_flux_op_value) += flux_val;

          if (flux_val >= 0)
            {
              (*convective_diffusive_flux_op_value_leaving)[l] = flux_val;
              (*sum_convective_diffusive_flux_op_value_leaving) += flux_val;
            }
          else
            {
              /*
               * Some cases where there are more than 4 faces to correct !!!
               * | |_| |
               * |/| |\|
               */
              (*convective_diffusive_flux_op_value_entering)[l] = flux_val;
              (*sum_convective_diffusive_flux_op_value_entering) += flux_val;
            }

          (*convective_diffusive_flux_op_value_normal_contrib)[l] = flux_val * normal_vector_compo_[face_dir[l]] * flux_out[l];
          (*sum_convective_diffusive_flux_op_value_normal_contrib) += flux_val * normal_vector_compo_[face_dir[l]] * flux_out[l];
        }
      else if (pure_vapour_neighbours_[l])
        {
          (*convective_diffusive_flux_op_value_vap)[l] = flux_val;
          (*sum_convective_diffusive_flux_op_value_vap) += flux_val;
        }
      else
        {
          (*convective_diffusive_flux_op_value_mixed)[l] = flux_val;
          (*sum_convective_diffusive_flux_op_value_mixed) += flux_val;
        }
    }
}
