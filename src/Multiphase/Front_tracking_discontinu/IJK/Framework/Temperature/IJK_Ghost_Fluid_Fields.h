/****************************************************************************
* Copyright (c) 2024, CEA
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

#ifndef IJK_Ghost_Fluid_Fields_included
#define IJK_Ghost_Fluid_Fields_included

#include <Objet_U.h>
#include <IJK_Field_vector.h>
#include <IJK_Field.h>

class Probleme_FTD_IJK_base;

class IJK_Ghost_Fluid_Fields : public Objet_U
{

  Declare_instanciable( IJK_Ghost_Fluid_Fields ) ;

public :

  void associer(const Probleme_FTD_IJK_base& ijk_ft);
  void initialize(int& nalloc, const Domaine_IJK& splitting);

  void compute_eulerian_distance();
  void enforce_zero_value_eulerian_distance();

  void enforce_distance_curvature_values_for_post_processings();

  void compute_eulerian_curvature_from_interface();
  void compute_eulerian_curvature();
  void enforce_zero_value_eulerian_curvature();
  void enforce_max_value_eulerian_curvature();

  void reset_field_computations()
  {
    has_computed_distance_ = false;
    has_computed_curvature_ = false;
  }
  void retrieve_ghost_fluid_params(const int& compute_distance,
                                   const int& compute_curvature,
                                   const int& n_iter_distance,
                                   const int& avoid_gfm_parallel_calls,
                                   const IJK_Field_local_double& boundary_flux_kmin,
                                   const IJK_Field_local_double& boundary_flux_kmax)
  {
    compute_distance_ = compute_distance;
    compute_curvature_ = compute_curvature;
    n_iter_distance_ = n_iter_distance;
    avoid_gfm_parallel_calls_ = avoid_gfm_parallel_calls;
    use_n_iter_distance_ = 0;
    nb_cells_gfm_parallel_calls_ = 2;
    boundary_flux_kmin_ = boundary_flux_kmin;
    boundary_flux_kmax_ = boundary_flux_kmax;
  }

  const IJK_Field_double& get_eulerian_distance_ft() const
  {
    return eulerian_distance_ft_;
  }
  const IJK_Field_double& get_eulerian_distance_ns() const
  {
    return eulerian_distance_ns_;
  }

  const IJK_Field_vector3_double& get_eulerian_normal_vectors_ft() const
  {
    return eulerian_normal_vectors_ft_;
  }
  const IJK_Field_vector3_double& get_eulerian_facets_barycentre_ft() const
  {
    return eulerian_facets_barycentre_ft_;
  }
  const IJK_Field_vector3_double& get_eulerian_normal_vectors_ns() const
  {
    return eulerian_normal_vectors_ns_;
  }
  const IJK_Field_vector3_double& get_eulerian_normal_vectors_ns_normed() const
  {
    return eulerian_normal_vectors_ns_normed_;
  }
  const IJK_Field_vector3_double& get_eulerian_facets_barycentre_ns() const
  {
    return eulerian_facets_barycentre_ns_;
  }

  const IJK_Field_double& get_eulerian_curvature_ft() const
  {
    return eulerian_curvature_ft_;
  }
  const IJK_Field_double& get_eulerian_curvature_ns() const
  {
    return eulerian_curvature_ns_;
  }
  const IJK_Field_double& get_eulerian_interfacial_area_ft() const
  {
    return eulerian_interfacial_area_ft_;
  }
  const IJK_Field_double& get_eulerian_interfacial_area_ns() const
  {
    return eulerian_interfacial_area_ns_;
  }

  const IJK_Field_double& get_eulerian_rising_velocities() const
  {
    return eulerian_rising_velocities_;
  }

  const IJK_Field_vector3_int& get_dummy_int_vect() const
  {
    return dummy_int_vect_;
  }
  const IJK_Field_vector3_double& get_dummy_double_vect() const
  {
    return dummy_double_vect_;
  }
  const IJK_Field_int& get_dummy_int_field() const
  {
    return dummy_int_field_;
  }
  const IJK_Field_double& get_dummy_double_field() const
  {
    return dummy_double_field_;
  }

protected :
  IJK_Field_double eulerian_distance_ft_;
  IJK_Field_double eulerian_distance_ns_;

  IJK_Field_vector3_double eulerian_normal_vectors_ft_;
  IJK_Field_vector3_double eulerian_facets_barycentre_ft_;
  IJK_Field_vector3_double eulerian_normal_vectors_ns_;
  IJK_Field_vector3_double eulerian_normal_vectors_ns_normed_;
  IJK_Field_vector3_double eulerian_facets_barycentre_ns_;

  IJK_Field_vector3_double tmp_old_vector_val_;
  IJK_Field_vector3_double tmp_new_vector_val_;

  IJK_Field_int tmp_interf_cells_;
  IJK_Field_int tmp_propagated_cells_;
  FixedVector<ArrOfInt,3> interf_cells_indices_;
  FixedVector<ArrOfInt,3> gfm_first_cells_indices_;
  FixedVector<ArrOfInt,3> propagated_cells_indices_;

  IJK_Field_double tmp_old_dist_val_;
  IJK_Field_double tmp_new_dist_val_;
  IJK_Field_double tmp_old_curv_val_;
  IJK_Field_double tmp_new_curv_val_;

  IJK_Field_double eulerian_curvature_ft_;
  IJK_Field_double eulerian_curvature_ns_;
  IJK_Field_double eulerian_interfacial_area_ft_;
  IJK_Field_double eulerian_interfacial_area_ns_;

  IJK_Field_double eulerian_rising_velocities_;

  IJK_Field_vector3_int dummy_int_vect_;
  IJK_Field_vector3_double dummy_double_vect_;
  IJK_Field_int dummy_int_field_;
  IJK_Field_double dummy_double_field_;

  OBS_PTR(Probleme_FTD_IJK_base) ref_ijk_ft_;
  int compute_distance_ = 1;
  int compute_curvature_ = 1;
  int n_iter_distance_ = 8;
  int avoid_gfm_parallel_calls_ = 0;
  int use_n_iter_distance_ = 0;
  int nb_cells_gfm_parallel_calls_ = 2;
  IJK_Field_local_double boundary_flux_kmin_;
  IJK_Field_local_double boundary_flux_kmax_;

  bool has_computed_distance_ = false;
  bool has_computed_curvature_ = false;
};

#endif /* IJK_Ghost_Fluid_Fields_included */
