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

#ifndef POST_PROCESSING_HYDRODYNAMIC_FORCES_STOKES
#define POST_PROCESSING_HYDRODYNAMIC_FORCES_STOKES

#include <Post_Processing_Hydrodynamic_Forces.h>

class Post_Processing_Hydrodynamic_Forces_Stokes: public Post_Processing_Hydrodynamic_Forces
{
  Declare_instanciable_sans_constructeur(Post_Processing_Hydrodynamic_Forces_Stokes);

public:
  Post_Processing_Hydrodynamic_Forces_Stokes();
  void compute_hydrodynamic_forces() override;
  void set_param(const Post_Processing_Hydrodynamic_Forces& post_process_hydro_forces_);

  const DoubleTab& get_pressure_force_Stokes_th() const { return total_pressure_force_Stokes_th_; }
  const DoubleTab& get_friction_force_Stokes_th() const { return total_friction_force_Stokes_th_; }
  const DoubleTab& get_pressure_force_Stokes_th_fa7() const { return pressure_force_Stokes_th_fa7_; }
  const DoubleTab& get_friction_force_Stokes_th_fa7() const { return friction_force_Stokes_th_fa7_; }

  void resize_and_init_tables(int nb_particles_tot) override;
  const DoubleTab& get_sigma_xx_fa7_Stokes_th() const { return sigma_xx_fa7_Stokes_th_; }
  const DoubleTab& get_sigma_xy_fa7_Stokes_th() const { return sigma_xy_fa7_Stokes_th_; }
  const DoubleTab& get_sigma_xz_fa7_Stokes_th() const { return sigma_xz_fa7_Stokes_th_; }
  const DoubleTab& get_sigma_yy_fa7_Stokes_th() const { return sigma_yy_fa7_Stokes_th_; }
  const DoubleTab& get_sigma_yz_fa7_Stokes_th() const { return sigma_yz_fa7_Stokes_th_; }
  const DoubleTab& get_sigma_zz_fa7_Stokes_th() const { return sigma_zz_fa7_Stokes_th_; }
  const DoubleTab& get_dUdx_P1_Stokes_th() { return dUdx_P1_Stokes_th_; }
  const DoubleTab& get_dUdz_P1_Stokes_th() { return dUdz_P1_Stokes_th_; }
  const DoubleTab& get_dVdz_P1_Stokes_th() { return dVdz_P1_Stokes_th_; }
  const DoubleTab& get_dWdx_P1_Stokes_th() { return dWdx_P1_Stokes_th_; }
  const DoubleTab& get_dWdy_P1_Stokes_th() { return dWdy_P1_Stokes_th_; }
  const DoubleTab& get_dWdz_P1_Stokes_th() { return dWdz_P1_Stokes_th_; }
  const DoubleTab& get_dUdx_P2_Stokes_th() { return dUdx_P2_Stokes_th_; }
  const DoubleTab& get_dUdz_P2_Stokes_th() { return dUdz_P2_Stokes_th_; }
  const DoubleTab& get_dVdz_P2_Stokes_th() { return dVdz_P2_Stokes_th_; }
  const DoubleTab& get_dWdx_P2_Stokes_th() { return dWdx_P2_Stokes_th_; }
  const DoubleTab& get_dWdy_P2_Stokes_th() { return dWdy_P2_Stokes_th_; }
  const DoubleTab& get_dWdz_P2_Stokes_th() { return dWdz_P2_Stokes_th_; }
  const DoubleTab& get_U_P1_Stokes_th() const { return U_P1_Stokes_th_; }
  const DoubleTab& get_U_P2_Stokes_th() const { return U_P2_Stokes_th_; }
  const DoubleVect& get_pressure_fa7_Stokes_th() const { return pressure_fa7_Stokes_th_; };

protected:
  void resize_data_fa7(int nb_fa7) override;
  void resize_sigma(int nb_fa7) override;
  void resize_gradU_P1(int nb_fa7) override;
  void resize_gradU_P2(int nb_fa7) override;
  void fill_Stokes_velocity_field();
  void fill_Stokes_pressure_field();
  void compute_vinf_Stokes();
  void compute_UP1_UP2_Stokes(int fa7,double vinf_Stokes,
                              double particle_radius, double phi_mu);
  double compute_pressure_interf(double x, double y, double z);
  double compute_Stokes_Ux_fluid(double x, double y, double z,
                                 double vinf_Stokes, double particle_radius, double phi_mu);
  double compute_Stokes_Uy_fluid(double x, double y, double z,
                                 double vinf_Stokes, double particle_radius, double phi_mu);
  double compute_Stokes_Uz_fluid(double x, double y, double z,
                                 double vinf_Stokes, double particle_radius, double phi_mu);

  const double& get_vinf_Stokes() const { return vinf_Stokes_; }

  void compute_neighbors_coordinates_fluid_fa7(const int nb_fa7,
                                               const int is_discr_elem_diph,
                                               const DoubleTab& gravity_center_fa7,
                                               const Maillage_FT_Disc& mesh,
                                               const DoubleTab& tab_fa7_normal,
                                               const IntTab& particles_eulerian_id_number) override;

  void compute_pressure_force_trilinear_linear(int nb_fa7,
                                               const Maillage_FT_Disc& mesh,
                                               Convection_Diffusion_Temperature_FT_Disc::
                                               Thermal_correction_discretization_method
                                               dummy_value,
                                               const IntVect& compo_connexes_fa7,
                                               const ArrOfDouble& fa7_surface,
                                               const DoubleTab& tab_fa7_normal) override;

  void compute_friction_force_complet_tensor(int nb_fa7,
                                             const Maillage_FT_Disc& mesh,
                                             const IntVect& compo_connexes_fa7,
                                             const ArrOfDouble& fa7_surface,
                                             const DoubleTab& tab_fa7_normal) override;

  void compute_friction_force_projected_tensor(int nb_fa7,
                                               const Maillage_FT_Disc& mesh,
                                               const IntVect& compo_connexes_fa7,
                                               const ArrOfDouble& fa7_surface,
                                               const DoubleTab& tab_fa7_normal) override;

  double compute_dUdx_Stokes_th(int fa7, double x_fa7, double y_fa7, double z_fa7,
                                double r_fa7, double phi_mu, double particle_radius);
  double compute_dUdz_Stokes_th(int fa7, double x_fa7, double y_fa7, double z_fa7,
                                double r_fa7, double phi_mu, double particle_radius);
  double compute_dVdz_Stokes_th(int fa7, double x_fa7, double y_fa7, double z_fa7,
                                double r_fa7, double phi_mu, double particle_radius);
  double compute_dWdx_Stokes_th(int fa7, double x_fa7, double y_fa7, double z_fa7,
                                double r_fa7, double phi_mu, double particle_radius);
  double compute_dWdy_Stokes_th(int fa7, double x_fa7, double y_fa7, double z_fa7,
                                double r_fa7, double phi_mu, double particle_radius);
  double compute_dWdz_Stokes_th(int fa7, double x_fa7, double y_fa7, double z_fa7,
                                double r_fa7, double phi_mu, double particle_radius);

  void fill_gradU_P1_Stokes_th(int fa7, double phi_mu, double particle_radius);
  void fill_gradU_P2_Stokes_th(int fa7, double phi_mu, double particle_radius);
  void fill_sigma_Stokes_th(int fa7);



private:
  double vinf_Stokes_;

  DoubleTab total_pressure_force_Stokes_th_;
  DoubleTab total_friction_force_Stokes_th_;
  DoubleVect pressure_fa7_Stokes_th_;

  DoubleTab pressure_force_Stokes_th_fa7_;
  DoubleTab friction_force_Stokes_th_fa7_;

  DoubleTab sigma_xx_fa7_Stokes_th_;
  DoubleTab sigma_xy_fa7_Stokes_th_;
  DoubleTab sigma_xz_fa7_Stokes_th_;
  DoubleTab sigma_yy_fa7_Stokes_th_;
  DoubleTab sigma_yz_fa7_Stokes_th_;
  DoubleTab sigma_zz_fa7_Stokes_th_;

  DoubleTab dUdx_P1_Stokes_th_;
  DoubleTab dUdz_P1_Stokes_th_;
  DoubleTab dVdz_P1_Stokes_th_;
  DoubleTab dWdx_P1_Stokes_th_;
  DoubleTab dWdy_P1_Stokes_th_;
  DoubleTab dWdz_P1_Stokes_th_;
  DoubleTab dUdx_P2_Stokes_th_;
  DoubleTab dUdz_P2_Stokes_th_;
  DoubleTab dVdz_P2_Stokes_th_;
  DoubleTab dWdx_P2_Stokes_th_;
  DoubleTab dWdy_P2_Stokes_th_;
  DoubleTab dWdz_P2_Stokes_th_;

  DoubleTab U_P1_Stokes_th_;
  DoubleTab U_P2_Stokes_th_;
};

#endif
