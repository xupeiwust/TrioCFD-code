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
// File      : IJK_Ghost_Fluid_tools.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Field.h>
#include <IJK_Field_vector.h>
#include <IJK_Interfaces.h>
#include <Maillage_FT_IJK.h>

#define INVALID_TEST -1.e30
//Be coherent with LOCAL EPS of Intersection Interface IJK
#define LIQUID_INDICATOR_TEST 1.-1.e-12
#define VAPOUR_INDICATOR_TEST 1.e-12
#define NEIGHBOURS_I {-1, 1, 0, 0, 0, 0}
#define NEIGHBOURS_J {0, 0, -1, 1, 0, 0}
#define NEIGHBOURS_K {0, 0, 0, 0, -1, 1}
#define select_dir(a,x,y,z) ((a==0)?(x):((a==1)?(y):(z)))
#define selectxy(a,x,y) ((a==0)?(x):(y))

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_Ghost_Fluid_tools
//
// <Description of class IJK_Ghost_Fluid_tools>
//
/////////////////////////////////////////////////////////////////////////////

void compute_eulerian_normal_distance_facet_barycentre_field(const IJK_Interfaces& interface,
                                                             IJK_Field_double& distance,
                                                             IJK_Field_vector3_double& normal_vect,
                                                             IJK_Field_vector3_double& facets_barycentre,
                                                             IJK_Field_vector3_double& tmp_old_vector_val,
                                                             IJK_Field_vector3_double& tmp_new_vector_val,
                                                             IJK_Field_double& tmp_old_val,
                                                             IJK_Field_double& tmp_new_val,
                                                             IJK_Field_int& tmp_interf_cells,
                                                             IJK_Field_int& tmp_propagated_cells,
                                                             FixedVector<ArrOfInt,3>& interf_cells_indices,
                                                             FixedVector<ArrOfInt,3>& gfm_first_cells_indices_,
                                                             FixedVector<ArrOfInt,3>& propagated_cells_indices,
                                                             const int& n_iter,
                                                             const int& avoid_gfm_parallel_calls=0);

void compute_eulerian_curvature_field_from_distance_field(const IJK_Field_double& distance,
                                                          IJK_Field_double& curvature,
                                                          const IJK_Field_local_double& boundary_flux_kmin,
                                                          const IJK_Field_local_double& boundary_flux_kmax);

void compute_eulerian_curvature_field_from_normal_vector_field(const IJK_Field_vector3_double& normal_vect,
                                                               IJK_Field_double& curvature);

void compute_eulerian_curvature_field_from_interface(const IJK_Field_vector3_double& normal_vect,
                                                     const IJK_Interfaces& interfaces,
                                                     IJK_Field_double& interfacial_area,
                                                     IJK_Field_double& curvature,
                                                     IJK_Field_double& tmp_old_val,
                                                     IJK_Field_double& tmp_new_val,
                                                     const int& n_iter,
                                                     int igroup);

void compute_eulerian_normal_temperature_gradient_interface(const IJK_Field_double& distance,
                                                            const IJK_Field_double& indicator,
                                                            const IJK_Field_double& interfacial_area,
                                                            const IJK_Field_double& curvature,
                                                            const	IJK_Field_double& temperature,
                                                            IJK_Field_double& grad_T_interface,
                                                            const int& spherical_approx,
                                                            const double& temperature_interf=0);

void propagate_eulerian_normal_temperature_gradient_interface(const IJK_Interfaces& interfaces,
                                                              const IJK_Field_double& distance,
                                                              IJK_Field_double& grad_T_interface,
                                                              const int& stencil_width,
                                                              const int& recompute_field_ini=1,
                                                              const int& zero_neighbour_value_mean=0,
                                                              const int& vapour_mixed_only=1,
                                                              const int& smooth_factor=10);

void compute_eulerian_extended_temperature(const IJK_Field_double& indicator,
                                           const IJK_Field_double& distance,
                                           const IJK_Field_double& curvature,
                                           IJK_Field_double& grad_T_interface,
                                           IJK_Field_double& temperature,
                                           const int& spherical_approx,
                                           const double& temperature_interf=0);

void smooth_vector_field(IJK_Field_vector3_double& vector_field,
                         IJK_Field_vector3_double& vector_field_init,
                         const IJK_Field_vector3_double * eulerian_normal_vectors_ns_normed,
                         const IJK_Interfaces& interfaces,
                         const double (&direct_smoothing_factors) [7],
                         const double (&gaussian_smoothing_factors) [3][3][3],
                         const int& smooth_numbers=1,
                         const int& remove_normal_compo=0,
                         const int& direct_neighbours=0,
                         const int& use_field_init=1,
                         const int& use_unique_phase=1);

void smooth_eulerian_field(IJK_Field_double& field,
                           IJK_Field_double& field_init,
                           const int& c,
                           IJK_Field_vector3_double& vector_field_init,
                           const IJK_Field_vector3_double * eulerian_normal_vectors_ns_normed,
                           const IJK_Interfaces& interfaces,
                           const double (&direct_smoothing_factors) [7],
                           const double (&gaussian_smoothing_factors) [3][3][3],
                           const int& smooth_numbers=1,
                           const int& remove_normal_compo=0,
                           const int& direct_neighbours=0,
                           const int& use_field_init=1,
                           const int& use_unique_phase=1);

void fill_tangential_gradient(const IJK_Field_vector3_double& vector_field,
                              const IJK_Field_vector3_double * eulerian_normal_vectors_ns_normed,
                              IJK_Field_vector3_double& tangential_vector_field);

void fill_tangential_gradient_compo(const IJK_Field_vector3_double& vector_field,
                                    const IJK_Field_vector3_double * eulerian_normal_vectors_ns_normed,
                                    IJK_Field_double& tangential_field,
                                    const int& dir);

