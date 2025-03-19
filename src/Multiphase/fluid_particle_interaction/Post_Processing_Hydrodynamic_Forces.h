/****************************************************************************
* Copyright (c) 2025, CEA
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

#ifndef POST_PROCESSING_HYDRODYNAMIC_FORCES_H
#define POST_PROCESSING_HYDRODYNAMIC_FORCES_H

#include <Convection_Diffusion_Temperature_FT_Disc.h>
#include <Maillage_FT_Disc.h>
#include <Equation_base.h>
#include <Matrice_Dense.h>
#include <TRUST_Ref.h>
#include <Objet_U.h>

class Transport_Interfaces_FT_Disc;
class Navier_Stokes_FT_Disc;

class Post_Processing_Hydrodynamic_Forces: public Objet_U
{
  Declare_instanciable_sans_constructeur(Post_Processing_Hydrodynamic_Forces);

public:
  Post_Processing_Hydrodynamic_Forces();

  void set_param(Param& p);
  int lire_motcle_non_standard(const Motcle&, Entree&) override;

  void associate_transport_equation(Transport_Interfaces_FT_Disc& ptr_eq_transport)
  { ptr_eq_transport_=ptr_eq_transport; }
  void associate_ns_equation(const Navier_Stokes_FT_Disc& eqn)
  { ptr_eq_ns_=eqn; }
  void associate_temp_equation(const Convection_Diffusion_Temperature_FT_Disc& eqn)
  { ptr_eq_temp_ = eqn; }

  void drop_the_flag() { flag_force_computation_=0; }
  void raise_the_flag() { flag_force_computation_=1; }
  void drop_the_flag_heat_transfer() { flag_heat_transfer_computation_=0; }
  void raise_the_flag_heat_transfer() { flag_heat_transfer_computation_=1; }
  virtual void compute_hydrodynamic_forces();
  void compute_heat_transfer();

  const DoubleVect& get_pressure_fa7() const { return pressure_fa7_; }
  const DoubleVect& get_total_surface_interf() const { return total_surface_interf_; }
  const DoubleVect& get_heat_transfer() const { return total_heat_transfer_; }
  const DoubleTab& get_pressure_force() const { return total_pressure_force_; }
  const DoubleTab& get_friction_force() const { return total_friction_force_; }
  const DoubleTab& get_pressure_force_fa7() const { return pressure_force_fa7_; }
  const DoubleTab& get_friction_force_fa7() const { return friction_force_fa7_; }
  const DoubleTab& get_heat_transfer_fa7() const { return heat_transfer_fa7_; }
  const DoubleTab& get_U_P1() const { return U_P1_; }
  const DoubleTab& get_U_P2() const { return U_P2_; }
  const DoubleTab& get_U_P2_moy() const { return U_P2_moy_; }
  const DoubleTab& get_prop_fa7_ok_P2() const { return proportion_fa7_ok_UP2_; }
  const DoubleTab& get_prop_P2_fluid() const { return prop_P2_fluid_compo_; }
  const DoubleTab& get_sigma_xx_fa7() const { return sigma_xx_fa7_; }
  const DoubleTab& get_sigma_xy_fa7() const { return sigma_xy_fa7_; }
  const DoubleTab& get_sigma_xz_fa7() const { return sigma_xz_fa7_; }
  const DoubleTab& get_sigma_yx_fa7() const { return sigma_yx_fa7_; }
  const DoubleTab& get_sigma_yy_fa7() const { return sigma_yy_fa7_; }
  const DoubleTab& get_sigma_yz_fa7() const { return sigma_yz_fa7_; }
  const DoubleTab& get_sigma_zx_fa7() const { return sigma_zx_fa7_; }
  const DoubleTab& get_sigma_zy_fa7() const { return sigma_zy_fa7_; }
  const DoubleTab& get_sigma_zz_fa7() const { return sigma_zz_fa7_; }
  const DoubleTab& get_dUdx_P1() { return dUdx_P1_; }
  const DoubleTab& get_dUdy_P1() { return dUdy_P1_; }
  const DoubleTab& get_dUdz_P1() { return dUdz_P1_; }
  const DoubleTab& get_dVdx_P1() { return dVdx_P1_; }
  const DoubleTab& get_dVdy_P1() { return dVdy_P1_; }
  const DoubleTab& get_dVdz_P1() { return dVdz_P1_; }
  const DoubleTab& get_dWdx_P1() { return dWdx_P1_; }
  const DoubleTab& get_dWdy_P1() { return dWdy_P1_; }
  const DoubleTab& get_dWdz_P1() { return dWdz_P1_; }
  const DoubleTab& get_dUdx_P2() { return dUdx_P2_; }
  const DoubleTab& get_dUdy_P2() { return dUdy_P2_; }
  const DoubleTab& get_dUdz_P2() { return dUdz_P2_; }
  const DoubleTab& get_dVdx_P2() { return dVdx_P2_; }
  const DoubleTab& get_dVdy_P2() { return dVdy_P2_; }
  const DoubleTab& get_dVdz_P2() { return dVdz_P2_; }
  const DoubleTab& get_dWdx_P2() { return dWdx_P2_; }
  const DoubleTab& get_dWdy_P2() { return dWdy_P2_; }
  const DoubleTab& get_dWdz_P2() { return dWdz_P2_; }

  int get_is_compute_forces() const { return is_compute_forces_; }
  int get_is_compute_forces_Stokes_th() const { return is_compute_stokes_theoretical_forces_; }
  int get_is_compute_heat_transfer() const { return is_compute_heat_transfer_; }
  int get_is_post_process_pressure_force_fa7() const { return is_post_process_pressure_force_fa7_; }
  int get_is_post_process_friction_force_fa7() const { return is_post_process_friction_force_fa7_; }
  int get_is_post_process_stress_tensor_fa7() const { return is_post_process_stress_tensor_fa7_; }
  int get_is_post_process_pressure_fa7() const { return is_post_process_pressure_fa7_; }

  const double& get_interpolation_distance_pressure_P1() const { return interpolation_distance_pressure_P1_; }
  const double& get_interpolation_distance_pressure_P2() const { return interpolation_distance_pressure_P2_; }
  const double& get_interpolation_distance_gradU_P1() const { return interpolation_distance_gradU_P1_; }
  const double& get_interpolation_distance_gradU_P2() const { return interpolation_distance_gradU_P2_; }

  virtual void resize_and_init_tables(int nb_compo_tot);

  enum class Method_friction_force_computation
  {
    TRILINEAR_LINEAR_COMPLET_TENSOR,
    TRILINEAR_LINEAR_PROJECTED_TENSOR
  };
  enum class Method_pressure_force_computation { TRILINEAR_LINEAR };
  enum class Location_stress_tensor { FACES_NORMALE_X, ELEMENTS };
  const Method_pressure_force_computation& get_method_pressure_force_computation()
  const { return method_pressure_force_computation_; }
  const Method_friction_force_computation& get_method_friction_force_computation()
  const { return method_friction_force_computation_; }
  const Location_stress_tensor& get_location_stress_tensor()
  const { return location_stress_tensor_; }

protected:
  OBS_PTR(Transport_Interfaces_FT_Disc) ptr_eq_transport_;
  OBS_PTR(Navier_Stokes_FT_Disc) ptr_eq_ns_;
  OBS_PTR(Convection_Diffusion_Temperature_FT_Disc) ptr_eq_temp_;

  int elem_faces_for_interp(int num_face,int i) const;
  int face_voisins_for_interp(int num_face,int i) const;

  int trilinear_interpolation_face(const DoubleTab& valeurs_champ,
                                   DoubleTab& coord,
                                   DoubleTab& resu);
  int trilinear_interpolation_elem(const DoubleTab& valeurs_champ,
                                   DoubleTab& coord,
                                   DoubleTab& resu);
  int trilinear_interpolation_elem(const DoubleTab& valeurs_champ,
                                   DoubleTab& coord,
                                   DoubleTab& resu,
                                   const int is_P2,
                                   const Convection_Diffusion_Temperature_FT_Disc::
                                   Thermal_correction_discretization_method discr);
  int trilinear_interpolation_gradU_face(const DoubleTab& valeurs_champ,
                                         DoubleTab& coord,
                                         DoubleTab& resu);
  int trilinear_interpolation_gradU_elem_P1(const DoubleTab& valeurs_champ,
                                            DoubleTab& coord,
                                            DoubleTab& resu);
  int trilinear_interpolation_gradU_elem(const DoubleTab& valeurs_champ,
                                         DoubleTab& coord,
                                         DoubleTab& resu);
  int trilinear_interpolation_face_sommets(const DoubleTab& valeurs_champ,
                                           DoubleTab& coord,
                                           DoubleTab& resu);

  double compute_viscosity_edges_sphere(int face1, int face2, int particle_id);
  double find_neighboring_elements(DoubleVect& coord_elem_interp,
                                   IntVect& neighboring_elements,
                                   const int sauv_list_P1=0,
                                   const int num_fa7=-1);

  void find_neighboring_faces (DoubleVect& coord_elem_interp,
                               IntVect& neighboring_faces,
                               int orientation) const;
  void find_neighboring_faces_xyz (DoubleVect& coord_elem_interp, IntTab& neighboring_faces) const;

  virtual void resize_sigma(int nb_fa7);
  virtual void resize_gradU_P1(int nb_fa7);
  virtual void resize_gradU_P2(int nb_fa7);
  virtual void resize_data_fa7(int nb_fa7);
  void resize_coord_neighbor_fluid_fa7(int nb_fa7);
  void fill_absurd_value_coord_neighbor_fluid_fa7(int fa7);
  void fill_sigma(int fa7, Matrice_Dense stress_tensor);
  void fill_gradU_P1(int fa7, DoubleTab gradU_P1);
  void fill_gradU_P2(int fa7, DoubleTab gradU_P2);

  Method_pressure_force_computation method_pressure_force_computation_=
    Method_pressure_force_computation::TRILINEAR_LINEAR;
  Method_friction_force_computation method_friction_force_computation_=
    Method_friction_force_computation::TRILINEAR_LINEAR_PROJECTED_TENSOR;
  Location_stress_tensor location_stress_tensor_=Location_stress_tensor::ELEMENTS;

  virtual void compute_pressure_force_trilinear_linear(int nb_fa7,
                                                       const Maillage_FT_Disc& mesh,
                                                       Convection_Diffusion_Temperature_FT_Disc::
                                                       Thermal_correction_discretization_method
                                                       thermal_correction_discretization_method,
                                                       const IntVect& compo_connexes_fa7,
                                                       const ArrOfDouble& fa7_surface,
                                                       const DoubleTab& tab_fa7_normal);

  virtual void compute_friction_force_complet_tensor(int nb_fa7,
                                                     const Maillage_FT_Disc& mesh,
                                                     const IntVect& compo_connexes_fa7,
                                                     const ArrOfDouble& fa7_surface,
                                                     const DoubleTab& tab_fa7_normal);


  virtual void compute_friction_force_projected_tensor(int nb_fa7,
                                                       const Maillage_FT_Disc& mesh,
                                                       const IntVect& compo_connexes_fa7,
                                                       const ArrOfDouble& fa7_surface,
                                                       const DoubleTab& tab_fa7_normal);

  virtual void compute_neighbors_coordinates_fluid_fa7(const int nb_fa7,
                                                       const int is_discr_elem_diph,
                                                       const DoubleTab& gravity_center_fa7,
                                                       const Maillage_FT_Disc& mesh,
                                                       const DoubleTab& tab_fa7_normal,
                                                       const IntTab& particles_eulerian_id_number);

  void compute_U_P2_moy(const int nb_particles_tot);
  void compute_proportion_fa7_ok_and_is_fluid_P2(const int nb_particles_tot);

  int is_compute_stokes_theoretical_forces_=0;
  int is_post_process_pressure_force_fa7_=0;
  int is_post_process_friction_force_fa7_=0;
  int is_post_process_stress_tensor_fa7_=0;
  int is_post_process_pressure_fa7_=0;
  int is_compute_heat_transfer_=0;
  int is_compute_forces_=0;

  int flag_force_computation_=0;
  int flag_heat_transfer_computation_=0;

  double interpolation_distance_temperature_P1_=0;
  double interpolation_distance_temperature_P2_=0;
  double interpolation_distance_pressure_P1_=0;
  double interpolation_distance_pressure_P2_=0;
  double interpolation_distance_gradU_P1_=0;
  double interpolation_distance_gradU_P2_=0;

  DoubleVect phase_indicator_function_P1_;
  DoubleVect phase_indicator_function_P2_;
  DoubleVect total_surface_interf_;
  DoubleVect total_heat_transfer_;
  DoubleVect pressure_fa7_;

  mutable IntTab list_elem_P1_all_;
  mutable IntTab list_elem_P1_;
  IntTab Nb_fa7_tot_par_compo_;
  IntTab list_elem_diph_;

  DoubleTab proportion_fa7_ok_UP2_;
  DoubleTab total_pressure_force_;
  DoubleTab total_friction_force_;
  DoubleTab prop_P2_fluid_compo_;
  DoubleTab pressure_force_fa7_;
  DoubleTab friction_force_fa7_;
  DoubleTab heat_transfer_fa7_;
  DoubleTab U_P2_moy_;
  DoubleTab U_P1_;
  DoubleTab U_P2_;

  DoubleTab sigma_xx_fa7_;
  DoubleTab sigma_xy_fa7_;
  DoubleTab sigma_xz_fa7_;
  DoubleTab sigma_yx_fa7_;
  DoubleTab sigma_yy_fa7_;
  DoubleTab sigma_yz_fa7_;
  DoubleTab sigma_zx_fa7_;
  DoubleTab sigma_zy_fa7_;
  DoubleTab sigma_zz_fa7_;

  DoubleTab dUdx_P1_;
  DoubleTab dUdy_P1_;
  DoubleTab dUdz_P1_;
  DoubleTab dVdx_P1_;
  DoubleTab dVdy_P1_;
  DoubleTab dVdz_P1_;
  DoubleTab dWdx_P1_;
  DoubleTab dWdy_P1_;
  DoubleTab dWdz_P1_;

  DoubleTab dUdx_P2_;
  DoubleTab dUdy_P2_;
  DoubleTab dUdz_P2_;
  DoubleTab dVdx_P2_;
  DoubleTab dVdy_P2_;
  DoubleTab dVdz_P2_;
  DoubleTab dWdx_P2_;
  DoubleTab dWdy_P2_;
  DoubleTab dWdz_P2_;

  DoubleTab coord_neighbor_fluid_fa7_pressure_1_;
  DoubleTab coord_neighbor_fluid_fa7_pressure_2_;
  DoubleTab coord_neighbor_fluid_fa7_gradU_1_;
  DoubleTab coord_neighbor_fluid_fa7_gradU_2_;
  DoubleTab coord_neighbor_fluid_fa7_temp_1_;
  DoubleTab coord_neighbor_fluid_fa7_temp_2_;
};

#endif
