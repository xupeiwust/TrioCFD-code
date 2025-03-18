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

#include <Post_Processing_Hydrodynamic_Forces_Stokes.h>

#include <Convection_Diffusion_Temperature_FT_Disc.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Navier_Stokes_FT_Disc.h>
#include <Fluide_Incompressible.h>
#include <Connex_components_FT.h>
#include <Solid_Particle_base.h>
#include <Fluide_Diphasique.h>
#include <Maillage_FT_Disc.h>
#include <Matrice_Dense.h>
#include <Domaine_VDF.h>
#include <TRUSTTab.h>

Implemente_instanciable_sans_constructeur(Post_Processing_Hydrodynamic_Forces_Stokes,
                                          "Post_Processing_Hydrodynamic_Forces_Stokes",
                                          Post_Processing_Hydrodynamic_Forces);

Post_Processing_Hydrodynamic_Forces_Stokes::Post_Processing_Hydrodynamic_Forces_Stokes()
{
}

Entree& Post_Processing_Hydrodynamic_Forces_Stokes::readOn(Entree& is)
{
  Cerr << "Error: Post_Processing_Hydrodynamic_Forces::printOn is not implemented." << finl;
  Process::exit();
  return is;
}

Sortie& Post_Processing_Hydrodynamic_Forces_Stokes::printOn(Sortie& os) const
{
  Cerr << "Error: Post_Processing_Hydrodynamic_Forces::printOn is not implemented." << finl;
  Process::exit();
  return os;
}

void Post_Processing_Hydrodynamic_Forces_Stokes::set_param(const Post_Processing_Hydrodynamic_Forces& post_process_hydro_forces)
{
  is_compute_forces_=post_process_hydro_forces.get_is_compute_forces();
  is_compute_stokes_theoretical_forces_=post_process_hydro_forces.get_is_compute_forces_Stokes_th();
  is_post_process_pressure_fa7_=post_process_hydro_forces.get_is_post_process_pressure_fa7();
  is_post_process_pressure_force_fa7_=post_process_hydro_forces.get_is_post_process_pressure_force_fa7();
  is_post_process_friction_force_fa7_=post_process_hydro_forces.get_is_post_process_friction_force_fa7();
  is_post_process_stress_tensor_fa7_=post_process_hydro_forces.get_is_post_process_stress_tensor_fa7();
  interpolation_distance_pressure_P1_=post_process_hydro_forces.get_interpolation_distance_pressure_P1();
  interpolation_distance_pressure_P2_=post_process_hydro_forces.get_interpolation_distance_pressure_P2();
  interpolation_distance_gradU_P1_=post_process_hydro_forces.get_interpolation_distance_gradU_P1();
  interpolation_distance_gradU_P2_=post_process_hydro_forces.get_interpolation_distance_gradU_P2();
  method_pressure_force_computation_=post_process_hydro_forces.get_method_pressure_force_computation();
  location_stress_tensor_=post_process_hydro_forces.get_location_stress_tensor();
}

void Post_Processing_Hydrodynamic_Forces_Stokes::compute_hydrodynamic_forces()
{
  if (!flag_force_computation_)
    {
      /***************************************************************
      *	  	  	  	  	  	RECUPERATION DES GRANDEURS		         *
      ***************************************************************/
      Cerr << "Post_Processing_Hydrodynamic_Forces_Stokes::compute_hydrodynamic_forces" << finl;
      Transport_Interfaces_FT_Disc& eq_transport = ptr_eq_transport_.valeur();
      const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
      const Maillage_FT_Disc& mesh = eq_transport.maillage_interface_pour_post();
      const int nb_fa7 = mesh.nb_facettes();

      IntVect compo_connexes_fa7(nb_fa7); // Init a zero
      int n = search_connex_components_local_FT(mesh, compo_connexes_fa7);
      int nb_particles_tot=compute_global_connex_components_FT(mesh, compo_connexes_fa7, n);

      resize_and_init_tables(nb_particles_tot);
      compute_vinf_Stokes();
      fill_Stokes_velocity_field();
      fill_Stokes_pressure_field();

      if (nb_fa7>0)
        {
          resize_data_fa7(nb_fa7);
          resize_coord_neighbor_fluid_fa7(nb_fa7);

          const DoubleTab& gravity_center_fa7=mesh.get_gravity_center_fa7();
          const ArrOfDouble& fa7_surface = mesh.get_update_surface_facettes();
          const DoubleTab& tab_fa7_normal = mesh.get_update_normale_facettes();

          const IntTab& particles_eulerian_id_number=eq_ns.get_particles_eulerian_id_number();
          compute_neighbors_coordinates_fluid_fa7(nb_fa7,
                                                  0,
                                                  gravity_center_fa7,
                                                  mesh,
                                                  tab_fa7_normal,
                                                  particles_eulerian_id_number);

          // ----------------------------- Pressure force  -----------------------------
          compute_pressure_force_trilinear_linear(nb_fa7,mesh,Convection_Diffusion_Temperature_FT_Disc::
                                                  Thermal_correction_discretization_method::P1,
                                                  compo_connexes_fa7, fa7_surface, tab_fa7_normal);
          // ----------------------------- Friction force -----------------------------
          if (method_friction_force_computation_==Method_friction_force_computation::
              TRILINEAR_LINEAR_COMPLET_TENSOR)
            {
              compute_friction_force_complet_tensor(nb_fa7, mesh, compo_connexes_fa7,
                                                    fa7_surface, tab_fa7_normal);
            }
          else if (method_friction_force_computation_==Method_friction_force_computation::
                   TRILINEAR_LINEAR_PROJECTED_TENSOR)
            {
              compute_friction_force_projected_tensor(nb_fa7, mesh, compo_connexes_fa7,
                                                      fa7_surface, tab_fa7_normal);
            }
        }
      mp_sum_for_each_item(total_pressure_force_Stokes_th_);
      mp_sum_for_each_item(total_friction_force_Stokes_th_);
      mp_sum_for_each_item(total_pressure_force_);
      mp_sum_for_each_item(total_friction_force_);

      raise_the_flag();
    }
}

void Post_Processing_Hydrodynamic_Forces_Stokes::resize_and_init_tables(int nb_particles_tot)
{
  total_pressure_force_Stokes_th_.resize(nb_particles_tot,dimension);
  total_friction_force_Stokes_th_.resize(nb_particles_tot,dimension);
  total_pressure_force_.resize(nb_particles_tot,dimension);
  total_friction_force_.resize(nb_particles_tot,dimension);
  total_pressure_force_Stokes_th_=0;
  total_friction_force_Stokes_th_=0;
  total_pressure_force_=0;
  total_friction_force_=0;
}

void Post_Processing_Hydrodynamic_Forces_Stokes::compute_vinf_Stokes()
{
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  const Fluide_Diphasique& two_phase_fluid = eq_ns.fluide_diphasique();
  const Solid_Particle_base& solid_particle=ref_cast(Solid_Particle_base,
                                                     two_phase_fluid.fluide_phase(0));
  const double particle_radius = solid_particle.get_equivalent_radius();
  const int id_fluid_phase=two_phase_fluid.get_id_fluid_phase();
  const double mu_f =two_phase_fluid.fluide_phase(id_fluid_phase).
                     viscosite_dynamique().valeurs()(0, 0);
  const double mu_p = two_phase_fluid.fluide_phase(1-id_fluid_phase).
                      masse_volumique().valeurs()(0, 0);
  const double rho_f =two_phase_fluid.fluide_phase(id_fluid_phase).
                      viscosite_dynamique().valeurs()(0, 0);
  const double rho_p = two_phase_fluid.fluide_phase(1-id_fluid_phase).
                       masse_volumique().valeurs()(0, 0);
  const double phi_mu=mu_p/mu_f;
  const DoubleTab& gravite = eq_ns.milieu().gravite().valeurs();
  DoubleVect vect_g(dimension);
  for (int i=0; i<dimension; i++)
    vect_g(i)=gravite(0,i);
  const double norme_g = sqrt(local_carre_norme_vect(vect_g));

  vinf_Stokes_=2./3.*(pow(particle_radius,2)*norme_g/mu_f)*((1+phi_mu)/(2.+3.*phi_mu)*(rho_f-rho_p));
}

void Post_Processing_Hydrodynamic_Forces_Stokes::fill_Stokes_velocity_field()
{
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  Navier_Stokes_FT_Disc& eq_ns_non_const = ptr_eq_ns_.valeur();

  const Domaine_VDF& domain_vdf = ref_cast(Domaine_VDF, eq_ns.domaine_dis());
  const DoubleTab& xv = domain_vdf.xv();
  const Fluide_Diphasique& two_phase_fluid = eq_ns.fluide_diphasique();
  const Solid_Particle_base& solid_particle=ref_cast(Solid_Particle_base,
                                                     two_phase_fluid.fluide_phase(0));
  const double particle_radius = solid_particle.get_equivalent_radius();
  const int id_fluid_phase=two_phase_fluid.get_id_fluid_phase();
  const double mu_f =two_phase_fluid.fluide_phase(id_fluid_phase).
                     viscosite_dynamique().valeurs()(0, 0);
  const double mu_p = two_phase_fluid.fluide_phase(1-id_fluid_phase).
                      masse_volumique().valeurs()(0, 0);
  const double phi_mu=mu_p/mu_f;
  const double& vinf_Stokes=get_vinf_Stokes();

  DoubleTab& vitesse_Stokes_th=eq_ns_non_const.get_set_velocity_field_Stokes_th();

  vitesse_Stokes_th=0;

  for (int num_face=0; num_face<vitesse_Stokes_th.dimension_tot(0); num_face++)
    {
      double x=xv(num_face,0);
      double y=xv(num_face,1);
      double z=xv(num_face,2);
      if (sqrt(x*x+y*y+z*z)>=particle_radius)
        {
          if (domain_vdf.orientation(num_face)==0)
            vitesse_Stokes_th(num_face)=compute_Stokes_Ux_fluid(x, y, z, vinf_Stokes,
                                                                particle_radius, phi_mu);
          else if (domain_vdf.orientation(num_face)==1)
            vitesse_Stokes_th(num_face)=compute_Stokes_Uy_fluid(x, y, z, vinf_Stokes,
                                                                particle_radius, phi_mu);
          else
            vitesse_Stokes_th(num_face)=compute_Stokes_Uz_fluid(x, y, z, vinf_Stokes,
                                                                particle_radius, phi_mu);
        }
      else
        {
          if (domain_vdf.orientation(num_face)==0)
            {
              vitesse_Stokes_th(num_face)=vinf_Stokes/(2*pow(particle_radius,2))*
                                          (1/(1+phi_mu))*x*z;
            }
          else if (domain_vdf.orientation(num_face)==1)
            {
              vitesse_Stokes_th(num_face)=vinf_Stokes/(2*pow(particle_radius,2))*
                                          (1/(1+phi_mu))*y*z;
            }
          else
            {
              vitesse_Stokes_th(num_face)=vinf_Stokes/(2*(1+phi_mu))*
                                          (1+z*z/(pow(particle_radius,2))-2*(x*x+y*y+z*z)/
                                           (pow(particle_radius,2)));
            }
        }
    }
  vitesse_Stokes_th.echange_espace_virtuel();
}

void Post_Processing_Hydrodynamic_Forces_Stokes::fill_Stokes_pressure_field()
{
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  Navier_Stokes_FT_Disc& eq_ns_non_const = ptr_eq_ns_.valeur();

  const Domaine_VDF& domain_vdf = ref_cast(Domaine_VDF, eq_ns.domaine_dis());
  const DoubleTab& xp = domain_vdf.xp();
  const Fluide_Diphasique& two_phase_fluid = eq_ns.fluide_diphasique();
  const Solid_Particle_base& solid_particle=ref_cast(Solid_Particle_base,
                                                     two_phase_fluid.fluide_phase(0));
  const double particle_radius = solid_particle.get_equivalent_radius();
  const int id_fluid_phase=two_phase_fluid.get_id_fluid_phase();
  const double mu_f =two_phase_fluid.fluide_phase(id_fluid_phase).
                     viscosite_dynamique().valeurs()(0, 0);
  const double mu_p = two_phase_fluid.fluide_phase(1-id_fluid_phase).
                      masse_volumique().valeurs()(0, 0);
  const double phi_mu=mu_p/mu_f;
  const double& vinf_Stokes=get_vinf_Stokes();
  DoubleTab& pressure_Stokes_th=eq_ns_non_const.get_set_pressure_field_Stokes_th();

  for (int num_elem=0; num_elem<pressure_Stokes_th.dimension_tot(0); num_elem++)
    {
      double x= xp(num_elem,0);
      double y= xp(num_elem,1);
      double z= xp(num_elem,2);
      if (sqrt(x*x+y*y+z*z)>=particle_radius)
        {
          pressure_Stokes_th(num_elem)=mu_f*vinf_Stokes*((2+3*phi_mu)/(1+phi_mu))*
                                       1/2*particle_radius*z/(pow(x*x+y*y+z*z,1.5));
        }
      else
        {
          pressure_Stokes_th(num_elem)=mu_p*vinf_Stokes*5/((1+phi_mu)*(pow(particle_radius,2)))*z;
        }
    }

  pressure_Stokes_th.echange_espace_virtuel();
}

void Post_Processing_Hydrodynamic_Forces_Stokes::resize_data_fa7(int nb_fa7)
{
  if (is_post_process_pressure_fa7_)
    {
      pressure_fa7_.resize(nb_fa7);
      pressure_fa7_=3e15;
      pressure_fa7_Stokes_th_.resize(nb_fa7);
      pressure_fa7_Stokes_th_=3e15;
    }

  pressure_force_Stokes_th_fa7_.resize(nb_fa7,dimension);
  pressure_force_Stokes_th_fa7_=3e15;

  friction_force_Stokes_th_fa7_.resize(nb_fa7,dimension);
  friction_force_Stokes_th_fa7_=3e15;

  pressure_force_fa7_.resize(nb_fa7,dimension);
  pressure_force_fa7_=3e15;

  friction_force_fa7_.resize(nb_fa7,dimension);
  friction_force_fa7_=3e15;

  if (is_post_process_stress_tensor_fa7_)
    {
      resize_sigma(nb_fa7);
      resize_gradU_P1(nb_fa7);
      resize_gradU_P2(nb_fa7);
    }

  U_P1_Stokes_th_.resize(nb_fa7,dimension);
  U_P2_Stokes_th_.resize(nb_fa7,dimension);
  U_P1_.resize(nb_fa7,dimension);
  U_P2_.resize(nb_fa7,dimension);
}

void Post_Processing_Hydrodynamic_Forces_Stokes::resize_sigma(int nb_fa7)
{
  sigma_xx_fa7_.resize(nb_fa7);
  sigma_xy_fa7_.resize(nb_fa7);
  sigma_xz_fa7_.resize(nb_fa7);
  sigma_yx_fa7_.resize(nb_fa7);
  sigma_yy_fa7_.resize(nb_fa7);
  sigma_yz_fa7_.resize(nb_fa7);
  sigma_zx_fa7_.resize(nb_fa7);
  sigma_zy_fa7_.resize(nb_fa7);
  sigma_zz_fa7_.resize(nb_fa7);

  sigma_xx_fa7_Stokes_th_.resize(nb_fa7);
  sigma_xy_fa7_Stokes_th_.resize(nb_fa7);
  sigma_xz_fa7_Stokes_th_.resize(nb_fa7);
  sigma_yy_fa7_Stokes_th_.resize(nb_fa7);
  sigma_yz_fa7_Stokes_th_.resize(nb_fa7);
  sigma_zz_fa7_Stokes_th_.resize(nb_fa7);
}

void Post_Processing_Hydrodynamic_Forces_Stokes::resize_gradU_P1(int nb_fa7)
{
  dUdx_P1_.resize(nb_fa7);
  dUdy_P1_.resize(nb_fa7);
  dUdz_P1_.resize(nb_fa7);
  dVdx_P1_.resize(nb_fa7);
  dVdy_P1_.resize(nb_fa7);
  dVdz_P1_.resize(nb_fa7);
  dWdx_P1_.resize(nb_fa7);
  dWdy_P1_.resize(nb_fa7);
  dWdz_P1_.resize(nb_fa7);

  dUdx_P1_Stokes_th_.resize(nb_fa7);
  dUdz_P1_Stokes_th_.resize(nb_fa7);
  dVdz_P1_Stokes_th_.resize(nb_fa7);
  dWdx_P1_Stokes_th_.resize(nb_fa7);
  dWdy_P1_Stokes_th_.resize(nb_fa7);
  dWdz_P1_Stokes_th_.resize(nb_fa7);
}

void Post_Processing_Hydrodynamic_Forces_Stokes::resize_gradU_P2(int nb_fa7)
{
  dUdx_P2_.resize(nb_fa7);
  dUdy_P2_.resize(nb_fa7);
  dUdz_P2_.resize(nb_fa7);
  dVdx_P2_.resize(nb_fa7);
  dVdy_P2_.resize(nb_fa7);
  dVdz_P2_.resize(nb_fa7);
  dWdx_P2_.resize(nb_fa7);
  dWdy_P2_.resize(nb_fa7);
  dWdz_P2_.resize(nb_fa7);

  dUdx_P2_Stokes_th_.resize(nb_fa7);
  dUdz_P2_Stokes_th_.resize(nb_fa7);
  dVdz_P2_Stokes_th_.resize(nb_fa7);
  dWdx_P2_Stokes_th_.resize(nb_fa7);
  dWdy_P2_Stokes_th_.resize(nb_fa7);
  dWdz_P2_Stokes_th_.resize(nb_fa7);
}

double Post_Processing_Hydrodynamic_Forces_Stokes::compute_pressure_interf(double x, double y, double z)
{
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  const Fluide_Diphasique& two_phase_fluid = eq_ns.fluide_diphasique();
  const Solid_Particle_base& solid_particle=ref_cast(Solid_Particle_base,
                                                     two_phase_fluid.fluide_phase(0));
  const double& vinf_Stokes=get_vinf_Stokes();
  const double particle_radius = solid_particle.get_equivalent_radius();
  const int id_fluid_phase=two_phase_fluid.get_id_fluid_phase();
  const double mu_f =two_phase_fluid.fluide_phase(id_fluid_phase).
                     viscosite_dynamique().valeurs()(0, 0);
  const double mu_p = two_phase_fluid.fluide_phase(1-id_fluid_phase).
                      masse_volumique().valeurs()(0, 0);
  const double phi_mu=mu_p/mu_f;
  return(mu_f*vinf_Stokes*((2+3*phi_mu)/(1+phi_mu))*1/2*particle_radius*z/(pow(x*x+y*y+z*z,1.5)));
}

double Post_Processing_Hydrodynamic_Forces_Stokes::compute_Stokes_Ux_fluid(double x, double y, double z,
                                                                           double vinf_Stokes, double particle_radius, double phi_mu)
{
  return (vinf_Stokes*z*x*((2+3*phi_mu)/(1+phi_mu)*(particle_radius)/(4*(pow(x*x+y*y+z*z,1.5)))-
                           3*phi_mu/(1+phi_mu)*pow(particle_radius,3)/(4*pow(x*x+y*y+z*z,2.5))));
}

double Post_Processing_Hydrodynamic_Forces_Stokes::compute_Stokes_Uy_fluid(double x, double y, double z,
                                                                           double vinf_Stokes, double particle_radius, double phi_mu)
{
  return (vinf_Stokes*z*y*((2+3*phi_mu)/(1+phi_mu)*(particle_radius)/(4*(pow(x*x+y*y+z*z,1.5)))-
                           3*phi_mu/(1+phi_mu)*pow(particle_radius,3)/(4*pow(x*x+y*y+z*z,2.5))));
}

double Post_Processing_Hydrodynamic_Forces_Stokes::compute_Stokes_Uz_fluid(double x, double y, double z,
                                                                           double vinf_Stokes, double particle_radius, double phi_mu)
{
  return (-vinf_Stokes*(+(z*z)*(1/(x*x+y*y+z*z)-(2+3*phi_mu)/(1+phi_mu)*
                                (particle_radius)/(2*pow(x*x+y*y+z*z,1.5))+phi_mu/(1+phi_mu)*
                                pow(particle_radius,3)/(2*pow(x*x+y*y+z*z,2.5)))+(1-z*z/(x*x+y*y+z*z))*
                        (1-(2+3*phi_mu)/(1+phi_mu)*particle_radius/(4*sqrt(x*x+y*y+z*z))-
                         phi_mu/(1+phi_mu)*pow(particle_radius,3)/(4*pow(x*x+y*y+z*z,1.5)))));
}

void Post_Processing_Hydrodynamic_Forces_Stokes::compute_UP1_UP2_Stokes(int fa7,double vinf_Stokes,
                                                                        double particle_radius, double phi_mu)
{
  double x=coord_neighbor_fluid_fa7_pressure_1_(fa7,0);
  double y=coord_neighbor_fluid_fa7_pressure_1_(fa7,1);
  double z=coord_neighbor_fluid_fa7_pressure_1_(fa7,2);

  U_P1_Stokes_th_(fa7,0) = compute_Stokes_Ux_fluid(x, y, z, vinf_Stokes, particle_radius, phi_mu);
  U_P1_Stokes_th_(fa7,1) = compute_Stokes_Uy_fluid(x, y, z, vinf_Stokes, particle_radius, phi_mu);
  U_P1_Stokes_th_(fa7,2) = compute_Stokes_Uz_fluid(x, y, z, vinf_Stokes, particle_radius, phi_mu);

  x=coord_neighbor_fluid_fa7_pressure_2_(fa7,0);
  y=coord_neighbor_fluid_fa7_pressure_2_(fa7,1);
  z=coord_neighbor_fluid_fa7_pressure_2_(fa7,2);

  U_P2_Stokes_th_(fa7,0) = compute_Stokes_Ux_fluid(x, y, z, vinf_Stokes, particle_radius, phi_mu);
  U_P2_Stokes_th_(fa7,1) = compute_Stokes_Uy_fluid(x, y, z, vinf_Stokes, particle_radius, phi_mu);
  U_P2_Stokes_th_(fa7,2) = compute_Stokes_Uz_fluid(x, y, z, vinf_Stokes, particle_radius, phi_mu);
}

void Post_Processing_Hydrodynamic_Forces_Stokes::compute_neighbors_coordinates_fluid_fa7(const int nb_fa7,
                                                                                         const int is_discr_elem_diph,
                                                                                         const DoubleTab& gravity_center_fa7,
                                                                                         const Maillage_FT_Disc& mesh,
                                                                                         const DoubleTab& tab_fa7_normal,
                                                                                         const IntTab& particles_eulerian_id_number)
{
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  const Domaine_VDF& domain_vdf = ref_cast(Domaine_VDF, eq_ns.domaine_dis());
  const Domaine& domain = domain_vdf.domaine();
  const double& vinf_Stokes=get_vinf_Stokes();
  const Fluide_Diphasique& two_phase_fluid = eq_ns.fluide_diphasique();
  const int id_fluid_phase=two_phase_fluid.get_id_fluid_phase();
  const Solid_Particle_base& solid_particle=ref_cast(Solid_Particle_base,
                                                     two_phase_fluid.fluide_phase(1-id_fluid_phase));
  const double particle_radius = solid_particle.get_equivalent_radius(); // On ne se sert de cette fonction que lorsqu'il y a une seule particule dans le domaine
  const double mu_f =two_phase_fluid.fluide_phase(id_fluid_phase).
                     viscosite_dynamique().valeurs()(0, 0);
  const double mu_p = two_phase_fluid.fluide_phase(1-id_fluid_phase).
                      masse_volumique().valeurs()(0, 0);
  const double phi_mu=mu_p/mu_f;

  for (int fa7 =0 ; fa7<nb_fa7 ; fa7++)
    {
      if (!mesh.facette_virtuelle(fa7))
        {
          DoubleVect normale_fa7(dimension);
          int elem_diph=domain.chercher_elements(gravity_center_fa7(fa7,0),
                                                 gravity_center_fa7(fa7,1),
                                                 gravity_center_fa7(fa7,2));
          DoubleVect delta_i(dimension);

          // On calcule les epaisseurs des mailles euleriennes  dans lesquelles se trouvent les facettes
          // Si on y a acces, on prend l'epaisseur a l'exterieur de la particule
          // Sinon, on prend l'epaisseur dans la particule
          // Cela revient simplement a choisir la maille juxtaposee a la maille diphasique
          for (int dim=0; dim<dimension; dim++)
            {
              int elem_haut=face_voisins_for_interp(elem_faces_for_interp(elem_diph, dim+dimension),1);
              int elem_bas=face_voisins_for_interp(elem_faces_for_interp(elem_diph, dim),0);
              if (tab_fa7_normal(fa7,dim)>0)
                delta_i(dim) =  (elem_haut>=0) ?
                                fabs(domain_vdf.dist_elem(elem_diph,elem_haut, dim)) :
                                fabs(domain_vdf.dist_elem(elem_diph,elem_bas, dim));
              else
                delta_i(dim) =  (elem_bas>=0) ?
                                fabs(domain_vdf.dist_elem(elem_diph,elem_bas, dim)) :
                                fabs(domain_vdf.dist_elem(elem_diph,elem_haut, dim));
            }

          double epsilon=0;
          for (int dim=0; dim<dimension; dim++)
            {
              epsilon+= fabs(delta_i(dim)*fabs(tab_fa7_normal(fa7,dim))); // la distance d'interpolation varie en fonction du raffinement du maillage
            }

          for (int dim=0; dim<dimension; dim++)
            {
              normale_fa7(dim)=tab_fa7_normal(fa7,dim);
              coord_neighbor_fluid_fa7_pressure_1_(fa7,dim)=gravity_center_fa7(fa7,dim)+
                                                            interpolation_distance_pressure_P1_*epsilon*normale_fa7(dim);
              coord_neighbor_fluid_fa7_pressure_2_(fa7,dim)=gravity_center_fa7(fa7,dim)+
                                                            interpolation_distance_pressure_P2_*epsilon*normale_fa7(dim);
              coord_neighbor_fluid_fa7_gradU_1_(fa7,dim)=gravity_center_fa7(fa7,dim)+
                                                         interpolation_distance_gradU_P1_*epsilon*normale_fa7(dim);
              coord_neighbor_fluid_fa7_gradU_2_(fa7,dim)=gravity_center_fa7(fa7,dim)+
                                                         interpolation_distance_gradU_P2_*epsilon*normale_fa7(dim);
            }

          compute_UP1_UP2_Stokes(fa7, vinf_Stokes, particle_radius, phi_mu);
        }
    }
}

void Post_Processing_Hydrodynamic_Forces_Stokes::compute_pressure_force_trilinear_linear(int nb_fa7,
                                                                                         const Maillage_FT_Disc& mesh,
                                                                                         Convection_Diffusion_Temperature_FT_Disc::
                                                                                         Thermal_correction_discretization_method
                                                                                         dummy_value,
                                                                                         const IntVect& compo_connexes_fa7,
                                                                                         const ArrOfDouble& fa7_surface,
                                                                                         const DoubleTab& tab_fa7_normal)
{
  DoubleTab pressure_P1(nb_fa7);
  DoubleTab pressure_P2(nb_fa7);
  Navier_Stokes_FT_Disc& eq_ns_non_const = ptr_eq_ns_.valeur();
  DoubleTab& pressure_Stokes_th=eq_ns_non_const.get_set_pressure_field_Stokes_th();

  if (trilinear_interpolation_elem(pressure_Stokes_th,
                                   coord_neighbor_fluid_fa7_pressure_1_,pressure_P1) && trilinear_interpolation_elem(
        pressure_Stokes_th,
        coord_neighbor_fluid_fa7_pressure_2_, pressure_P2)) // soit on est capable d'interpoler en P1 et P1 et on calcule la force de p
    {
      const DoubleTab& gravity_center_fa7=mesh.get_gravity_center_fa7();
      const double& vinf_Stokes=get_vinf_Stokes();
      DoubleTab extrapolated_pressure(nb_fa7);
      const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
      const Fluide_Diphasique& two_phase_fluid = eq_ns.fluide_diphasique();
      const int id_fluid_phase=two_phase_fluid.get_id_fluid_phase();
      const Solid_Particle_base& solid_particle=ref_cast(Solid_Particle_base,
                                                         two_phase_fluid.fluide_phase(1-id_fluid_phase));
      const double particle_radius = solid_particle.get_equivalent_radius(); // On ne se sert de cette fonction que lorsqu'il y a une seule particule dans le domaine
      const double mu_f =two_phase_fluid.fluide_phase(id_fluid_phase).
                         viscosite_dynamique().valeurs()(0, 0);
      const double mu_p = two_phase_fluid.fluide_phase(1-id_fluid_phase).
                          masse_volumique().valeurs()(0, 0);
      const double phi_mu=mu_p/mu_f;

      for (int fa7=0; fa7<nb_fa7; fa7++)
        {
          if (!mesh.facette_virtuelle(fa7))
            {
              extrapolated_pressure(fa7) = pressure_P2(fa7)-interpolation_distance_pressure_P2_*
                                           (pressure_P2(fa7)-pressure_P1(fa7))/
                                           (interpolation_distance_pressure_P2_-
                                            interpolation_distance_pressure_P1_); //3*pressure_P2(compo,fa7)-2*pressure_P1(compo,fa7); // Si on n'est pas capable d'interpoler en P1 ET P2, alors on ne calcule pas la force de pressure
              pressure_fa7_(fa7)=extrapolated_pressure(fa7);
              pressure_fa7_Stokes_th_(fa7)=compute_pressure_interf(gravity_center_fa7(fa7,0),
                                                                   gravity_center_fa7(fa7,1),gravity_center_fa7(fa7,2));
            }
          else
            {
              pressure_fa7_(fa7)=1e15;
              extrapolated_pressure(fa7)=-1e15;
            }
        }

      for (int fa7=0; fa7<nb_fa7; fa7++)
        {
          int compo=compo_connexes_fa7(fa7);
          double coeff=-extrapolated_pressure(fa7)*fa7_surface(fa7);
          DoubleVect pressure_force_fa7(dimension);
          for (int dim=0; dim<dimension; dim++)
            pressure_force_fa7(dim)=coeff*tab_fa7_normal(fa7,dim);
          if (!mesh.facette_virtuelle(fa7))
            {

              for (int dim=0; dim<dimension; dim++)
                pressure_force_fa7_(fa7,dim)=pressure_force_fa7(dim);

              pressure_force_Stokes_th_fa7_(fa7,0)=-(mu_f*vinf_Stokes*(2.+3.*phi_mu)/
                                                     (1.+phi_mu)*1./(2.*particle_radius)*
                                                     gravity_center_fa7(fa7,2)/particle_radius)*
                                                   tab_fa7_normal(fa7,0)*fa7_surface(fa7);
              pressure_force_Stokes_th_fa7_(fa7,1)=-(mu_f*vinf_Stokes*(2.+3.*phi_mu)/(1.+phi_mu)*
                                                     1./(2.*particle_radius)*
                                                     gravity_center_fa7(fa7,2)/particle_radius)*
                                                   tab_fa7_normal(fa7,1)*fa7_surface(fa7);
              pressure_force_Stokes_th_fa7_(fa7,2)=-(mu_f*vinf_Stokes*(2.+3.*phi_mu)/(1.+phi_mu)*
                                                     1./(2.*particle_radius)*
                                                     gravity_center_fa7(fa7,2)/particle_radius)*
                                                   tab_fa7_normal(fa7,2)*fa7_surface(fa7);

              for (int dim=0; dim<dimension; dim++)
                {
                  total_pressure_force_(compo,dim)+=pressure_force_fa7(dim);
                  total_pressure_force_Stokes_th_(compo,dim)+=pressure_force_Stokes_th_fa7_(fa7,dim);
                }
            }
        }
    }
}

void Post_Processing_Hydrodynamic_Forces_Stokes::compute_friction_force_complet_tensor(int nb_fa7,
                                                                                       const Maillage_FT_Disc& mesh,
                                                                                       const IntVect& compo_connexes_fa7,
                                                                                       const ArrOfDouble& fa7_surface,
                                                                                       const DoubleTab& tab_fa7_normal)
{
  int interp_gradU_P1_ok=0;
  int interp_gradU_P2_ok=0;
  DoubleTab grad_U_P1(nb_fa7, dimension, dimension);
  DoubleTab grad_U_P2(nb_fa7, dimension, dimension);
  grad_U_P1=-1e15;
  grad_U_P2=-1e30;
  const double& vinf_Stokes=get_vinf_Stokes();
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  const Fluide_Diphasique& two_phase_fluid = eq_ns.fluide_diphasique();
  const int id_fluid_phase=two_phase_fluid.get_id_fluid_phase();
  const Solid_Particle_base& solid_particle=ref_cast(Solid_Particle_base,
                                                     two_phase_fluid.fluide_phase(1-id_fluid_phase));
  const double particle_radius = solid_particle.get_equivalent_radius(); // On ne se sert de cette fonction que lorsqu'il y a une seule particule dans le domaine
  const double mu_f =two_phase_fluid.fluide_phase(id_fluid_phase).
                     viscosite_dynamique().valeurs()(0, 0);
  const double mu_p = two_phase_fluid.fluide_phase(1-id_fluid_phase).
                      masse_volumique().valeurs()(0, 0);
  const double phi_mu=mu_p/mu_f;
  const DoubleTab& gravity_center_fa7=mesh.get_gravity_center_fa7();
  Navier_Stokes_FT_Disc& eq_ns_non_const = ptr_eq_ns_.valeur();
  DoubleTab& velocity_Stokes_th=eq_ns_non_const.get_set_velocity_field_Stokes_th();
  if (location_stress_tensor_== Location_stress_tensor::FACES_NORMALE_X)
    {
      interp_gradU_P1_ok=trilinear_interpolation_gradU_face(velocity_Stokes_th,
                                                            coord_neighbor_fluid_fa7_gradU_1_,
                                                            grad_U_P1);
      interp_gradU_P2_ok=trilinear_interpolation_gradU_face(velocity_Stokes_th,
                                                            coord_neighbor_fluid_fa7_gradU_2_,
                                                            grad_U_P2);
    }
  else if (location_stress_tensor_== Location_stress_tensor::ELEMENTS)
    {
      interp_gradU_P1_ok=trilinear_interpolation_gradU_elem(velocity_Stokes_th,
                                                            coord_neighbor_fluid_fa7_gradU_1_,
                                                            grad_U_P1);
      interp_gradU_P2_ok=trilinear_interpolation_gradU_elem(velocity_Stokes_th,
                                                            coord_neighbor_fluid_fa7_gradU_2_,
                                                            grad_U_P2);
    }
  if ( interp_gradU_P1_ok==1 && interp_gradU_P2_ok==1 )
    {
      DoubleTab grad_U_extrapole(nb_fa7, dimension, dimension);
      grad_U_extrapole=1e20;

      for (int fa7=0; fa7<nb_fa7; fa7++)
        {
          if (!mesh.facette_virtuelle(fa7))
            {
              fill_gradU_P1(fa7, grad_U_P1);
              fill_gradU_P2(fa7, grad_U_P2);
              fill_gradU_P1_Stokes_th(fa7,phi_mu,particle_radius);
              fill_gradU_P2_Stokes_th(fa7,phi_mu,particle_radius);
            }
          for (int i=0; i<dimension; i++)
            {
              for (int j=0; j<dimension; j++)
                {
                  grad_U_extrapole(fa7,i,j) = grad_U_P2(fa7,i,j)-
                                              interpolation_distance_gradU_P2_*(grad_U_P2(fa7,i,j)-
                                                                                grad_U_P1(fa7,i,j))/(interpolation_distance_gradU_P2_-
                                                                                                     interpolation_distance_gradU_P1_);
                }
            }
        }

      for (int fa7=0; fa7<nb_fa7; fa7++)
        {
          int compo=compo_connexes_fa7(fa7);
          Matrice_Dense stress_tensor(dimension,dimension);
          for (int i=0; i<dimension; i++)
            {
              for (int j=0; j<dimension; j++)
                {
                  stress_tensor(i,j) = (grad_U_extrapole(fa7,i,j) +
                                        grad_U_extrapole(fa7,j,i));
                }
            }

          if (is_post_process_stress_tensor_fa7_)
            {
              fill_sigma(fa7,stress_tensor);
              fill_sigma_Stokes_th(fa7);
            }

          DoubleTab la_normale_fa7_x_surface(dimension);

          for (int dim=0; dim<dimension; dim++)
            la_normale_fa7_x_surface(dim) = fa7_surface(fa7)*tab_fa7_normal(fa7,dim);

          DoubleVect friction_force_fa7=stress_tensor*la_normale_fa7_x_surface;

          if (!mesh.facette_virtuelle(fa7))
            {
              for (int dim=0; dim<dimension; dim++)
                friction_force_fa7_(fa7,dim)=friction_force_fa7(dim);

              friction_force_Stokes_th_fa7_(fa7,0)=-mu_f*vinf_Stokes*gravity_center_fa7(fa7,0)*
                                                   gravity_center_fa7(fa7,2)/(pow(particle_radius,3))*(9.*phi_mu+8.)/
                                                   (1.+phi_mu)*fa7_surface(fa7);
              friction_force_Stokes_th_fa7_(fa7,1)=-mu_f*vinf_Stokes*
                                                   gravity_center_fa7(fa7,1)*gravity_center_fa7(fa7,2)/(pow(particle_radius,3))*
                                                   (9.*phi_mu+8.)/(1.+phi_mu)*fa7_surface(fa7);
              friction_force_Stokes_th_fa7_(fa7,2)=mu_f*vinf_Stokes*1./(2.*
                                                                        particle_radius*(1.+phi_mu))*(-3.*phi_mu+pow(gravity_center_fa7(fa7,2)/
                                                                                                                     particle_radius,2)*(3.*phi_mu-4.))*fa7_surface(fa7);

              for (int dim=0; dim<dimension; dim++)
                {
                  total_friction_force_(compo,dim)+=friction_force_fa7_(dim);
                  total_friction_force_Stokes_th_(compo,dim)+=
                    friction_force_Stokes_th_fa7_(fa7,dim);
                }
            }
        }
    }
}

double Post_Processing_Hydrodynamic_Forces_Stokes::compute_dUdx_Stokes_th(int fa7, double x_fa7, double y_fa7, double z_fa7,
                                                                          double r_fa7, double phi_mu, double particle_radius)
{
  double vinf_Stokes=get_vinf_Stokes();
  double dUdx=(vinf_Stokes*(pow(x_fa7,2)*z_fa7/pow(r_fa7,3)*
                            ( -(2+3*phi_mu)/(1+phi_mu)*3*particle_radius/(pow(2*r_fa7,2)) +
                              phi_mu/(1+phi_mu)*9*pow(particle_radius,3)/(4*pow(r_fa7,4)) )  +
                            pow(x_fa7,2)/(pow(r_fa7,2)-
                                          pow(z_fa7,2))*z_fa7/r_fa7*
                            ( (1-pow(z_fa7/r_fa7,2))*( (2+3*phi_mu)/
                                                       (1+phi_mu)*particle_radius/
                                                       (4*pow(r_fa7,2)) + phi_mu/(1+phi_mu)*
                                                       3*pow(particle_radius,3)/(4*pow(r_fa7,4)) )
                              + pow(z_fa7/r_fa7,2)*( (2+3*phi_mu)/(1+phi_mu)*
                                                     particle_radius/(4*pow(r_fa7,2))-
                                                     phi_mu/(1+phi_mu)*3*pow(particle_radius,3)/
                                                     (4*pow(r_fa7,4)) ) )
                            + pow(y_fa7,2)*z_fa7/(r_fa7*(pow(r_fa7,2)
                                                         -pow(z_fa7,2)))*((2+3*phi_mu)/(1+phi_mu)*
                                                                          particle_radius/(pow(2*r_fa7,2))-
                                                                          phi_mu/(1+phi_mu)*3*
                                                                          pow(particle_radius,3)/(4*pow(r_fa7,4)))));
  return dUdx;
}
double Post_Processing_Hydrodynamic_Forces_Stokes::compute_dUdz_Stokes_th(int fa7, double x_fa7, double y_fa7, double z_fa7,
                                                                          double r_fa7, double phi_mu, double particle_radius)
{
  double vinf_Stokes=get_vinf_Stokes();
  double dUdz=vinf_Stokes*(x_fa7/r_fa7*(pow(z_fa7/r_fa7,2)*
                                        (-(2+3*phi_mu)/(1+phi_mu)*particle_radius/(2*pow(r_fa7,2))+phi_mu/
                                         (1+phi_mu)*3*pow(particle_radius,3)/(2*pow(r_fa7,4)))
                                        +(1-pow(z_fa7/r_fa7,2))*
                                        ((2+3*phi_mu)/(1+phi_mu)*particle_radius/
                                         (4*pow(r_fa7,2))-phi_mu/(1+phi_mu)*3*
                                         pow(particle_radius,3)/(4*pow(r_fa7,4))))+
                           pow(z_fa7,2)*x_fa7/pow(r_fa7,3)*
                           phi_mu/(1+phi_mu)*3*pow(particle_radius,3)/
                           (2*pow(r_fa7,4)));
  return dUdz;
}

double Post_Processing_Hydrodynamic_Forces_Stokes::compute_dVdz_Stokes_th(int fa7, double x_fa7, double y_fa7, double z_fa7,
                                                                          double r_fa7, double phi_mu, double particle_radius)
{
  double vinf_Stokes=get_vinf_Stokes();
  double dVdz=vinf_Stokes*(y_fa7/r_fa7*(pow(z_fa7/r_fa7,2)*
                                        (-(2+3*phi_mu)/(1+phi_mu)*particle_radius/(2*pow(r_fa7,2))+phi_mu/
                                         (1+phi_mu)*3*pow(particle_radius,3)/(2*pow(r_fa7,4)))
                                        +(1-pow(z_fa7/r_fa7,2))*((2+3*phi_mu)/(1+phi_mu)*
                                                                 particle_radius/(4*pow(r_fa7,2))-phi_mu/(1+phi_mu)*
                                                                 3*pow(particle_radius,3)/(4*pow(r_fa7,4))))+
                           pow(z_fa7,2)*y_fa7/pow(r_fa7,3)*
                           phi_mu/(1+phi_mu)*3*pow(particle_radius,3)/(2*pow(r_fa7,4)));

  return dVdz;
}

double Post_Processing_Hydrodynamic_Forces_Stokes::compute_dWdx_Stokes_th(int fa7, double x_fa7, double y_fa7, double z_fa7,
                                                                          double r_fa7, double phi_mu, double particle_radius)
{
  double vinf_Stokes=get_vinf_Stokes();
  double dWdx=vinf_Stokes*(pow(z_fa7,2)*x_fa7/pow(r_fa7,3)*
                           (-(2+3*phi_mu)/(1+phi_mu)*3*particle_radius/(pow(2*r_fa7,2))+
                            phi_mu/(1+phi_mu)*9*pow(particle_radius,3)/(4*pow(r_fa7,4)))-
                           x_fa7/r_fa7*((1-pow(z_fa7/r_fa7,2))*((2+3*phi_mu)/
                                                                (1+phi_mu)*particle_radius/pow(2*r_fa7,2)+
                                                                phi_mu/(1+phi_mu)*3*pow(particle_radius,3)/
                                                                (4*pow(r_fa7,4)))+pow(z_fa7/r_fa7,2)*
                                        ((2+3*phi_mu)/(1+phi_mu)*particle_radius/pow(2*r_fa7,2)-
                                         phi_mu/(1+phi_mu)*3*pow(particle_radius,3)/
                                         (4*pow(r_fa7,4)))));

  return dWdx;
}

double Post_Processing_Hydrodynamic_Forces_Stokes::compute_dWdy_Stokes_th(int fa7, double x_fa7, double y_fa7, double z_fa7,
                                                                          double r_fa7, double phi_mu, double particle_radius)
{
  double vinf_Stokes=get_vinf_Stokes();
  double dWdy=vinf_Stokes*(pow(z_fa7,2)*y_fa7/pow(r_fa7,3)*
                           (-(2+3*phi_mu)/(1+phi_mu)*3*particle_radius/(pow(2*r_fa7,2))+
                            phi_mu/(1+phi_mu)*9*pow(particle_radius,3)/(4*pow(r_fa7,4)))-
                           y_fa7/r_fa7*((1-pow(z_fa7/r_fa7,2))*((2+3*phi_mu)/
                                                                (1+phi_mu)*particle_radius/pow(2*r_fa7,2)+phi_mu/
                                                                (1+phi_mu)*3*pow(particle_radius,3)/(4*pow(r_fa7,4)))
                                        +pow(z_fa7/r_fa7,2)*((2+3*phi_mu)/(1+phi_mu)*
                                                             particle_radius/pow(2*r_fa7,2)-phi_mu/(1+phi_mu)*
                                                             3*pow(particle_radius,3)/(4*pow(r_fa7,4)))));

  return dWdy;
}

double Post_Processing_Hydrodynamic_Forces_Stokes::compute_dWdz_Stokes_th(int fa7, double x_fa7, double y_fa7, double z_fa7,
                                                                          double r_fa7, double phi_mu, double particle_radius)
{
  double vinf_Stokes=get_vinf_Stokes();
  double dWdz=vinf_Stokes*(z_fa7/r_fa7*(pow(z_fa7/r_fa7,2)*
                                        (-(2+3*phi_mu)/(1+phi_mu)*particle_radius/(2*pow(r_fa7,2))+
                                         phi_mu/(1+phi_mu)*3*pow(particle_radius,3)/(2*pow(r_fa7,4)))+
                                        (1-pow(z_fa7/r_fa7,2))*((2+3*phi_mu)/(1+phi_mu)*
                                                                particle_radius/(pow(2*r_fa7,2))-phi_mu/(1+phi_mu)*
                                                                3*pow(particle_radius,3)/(4*pow(r_fa7,4))))-
                           z_fa7/r_fa7*(1-pow(z_fa7/r_fa7,2))*
                           phi_mu/(1+phi_mu)*3*pow(particle_radius,3)/(2*pow(r_fa7,4)));

  return dWdz;
}

void Post_Processing_Hydrodynamic_Forces_Stokes::fill_gradU_P1_Stokes_th(int fa7,
                                                                         double phi_mu, double particle_radius)
{
  double x_fa7_P1=coord_neighbor_fluid_fa7_gradU_1_(fa7,0);
  double y_fa7_P1=coord_neighbor_fluid_fa7_gradU_1_(fa7,1);
  double z_fa7_P1=coord_neighbor_fluid_fa7_gradU_1_(fa7,2);
  double r_fa7_P1=sqrt(x_fa7_P1*x_fa7_P1+y_fa7_P1*y_fa7_P1+z_fa7_P1*z_fa7_P1);

  dUdx_P1_Stokes_th_(fa7)=compute_dUdx_Stokes_th(fa7, x_fa7_P1, y_fa7_P1, z_fa7_P1,
                                                 r_fa7_P1, phi_mu, particle_radius);
  dUdz_P1_Stokes_th_(fa7)=compute_dUdz_Stokes_th(fa7, x_fa7_P1, y_fa7_P1, z_fa7_P1,
                                                 r_fa7_P1, phi_mu, particle_radius);
  dVdz_P1_Stokes_th_(fa7)=compute_dVdz_Stokes_th(fa7, x_fa7_P1, y_fa7_P1, z_fa7_P1,
                                                 r_fa7_P1, phi_mu, particle_radius);
  dWdx_P1_Stokes_th_(fa7)=compute_dWdx_Stokes_th(fa7, x_fa7_P1, y_fa7_P1, z_fa7_P1,
                                                 r_fa7_P1, phi_mu, particle_radius);
  dWdy_P1_Stokes_th_(fa7)=compute_dWdy_Stokes_th(fa7, x_fa7_P1, y_fa7_P1, z_fa7_P1,
                                                 r_fa7_P1, phi_mu, particle_radius);
  dWdz_P1_Stokes_th_(fa7)=compute_dWdz_Stokes_th(fa7, x_fa7_P1, y_fa7_P1, z_fa7_P1,
                                                 r_fa7_P1, phi_mu, particle_radius);
}

void Post_Processing_Hydrodynamic_Forces_Stokes::fill_gradU_P2_Stokes_th(int fa7,
                                                                         double phi_mu, double particle_radius)
{
  double x_fa7_P2=coord_neighbor_fluid_fa7_gradU_2_(fa7,0);
  double y_fa7_P2=coord_neighbor_fluid_fa7_gradU_2_(fa7,1);
  double z_fa7_P2=coord_neighbor_fluid_fa7_gradU_2_(fa7,2);
  double r_fa7_P2=sqrt(x_fa7_P2*x_fa7_P2+y_fa7_P2*y_fa7_P2+z_fa7_P2*z_fa7_P2);

  dUdx_P2_Stokes_th_(fa7)=compute_dUdx_Stokes_th(fa7, x_fa7_P2, y_fa7_P2, z_fa7_P2,
                                                 r_fa7_P2, phi_mu, particle_radius);
  dUdz_P2_Stokes_th_(fa7)=compute_dUdz_Stokes_th(fa7, x_fa7_P2, y_fa7_P2, z_fa7_P2,
                                                 r_fa7_P2, phi_mu, particle_radius);
  dVdz_P2_Stokes_th_(fa7)=compute_dVdz_Stokes_th(fa7, x_fa7_P2, y_fa7_P2, z_fa7_P2,
                                                 r_fa7_P2, phi_mu, particle_radius);
  dWdx_P2_Stokes_th_(fa7)=compute_dWdx_Stokes_th(fa7, x_fa7_P2, y_fa7_P2, z_fa7_P2,
                                                 r_fa7_P2, phi_mu, particle_radius);
  dWdy_P2_Stokes_th_(fa7)=compute_dWdy_Stokes_th(fa7, x_fa7_P2, y_fa7_P2, z_fa7_P2,
                                                 r_fa7_P2, phi_mu, particle_radius);
  dWdz_P2_Stokes_th_(fa7)=compute_dWdz_Stokes_th(fa7, x_fa7_P2, y_fa7_P2, z_fa7_P2,
                                                 r_fa7_P2, phi_mu, particle_radius);
}

void Post_Processing_Hydrodynamic_Forces_Stokes::fill_sigma_Stokes_th(int fa7)
{
  const Transport_Interfaces_FT_Disc& eq_transport = ptr_eq_transport_.valeur();
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  const Maillage_FT_Disc& mesh = eq_transport.maillage_interface_pour_post();
  const Fluide_Diphasique& two_phase_fluid = eq_ns.fluide_diphasique();
  const int id_fluid_phase=two_phase_fluid.get_id_fluid_phase();
  const Solid_Particle_base& solid_particle=ref_cast(Solid_Particle_base,
                                                     two_phase_fluid.fluide_phase(1-id_fluid_phase));
  const double particle_radius = solid_particle.get_equivalent_radius(); // On ne se sert de cette fonction que lorsqu'il y a une seule particule dans le domaine
  const double mu_f =two_phase_fluid.fluide_phase(id_fluid_phase).
                     viscosite_dynamique().valeurs()(0, 0);
  const double mu_p = two_phase_fluid.fluide_phase(1-id_fluid_phase).
                      masse_volumique().valeurs()(0, 0);
  const double phi_mu=mu_p/mu_f;
  const double& vinf_Stokes=get_vinf_Stokes();

  const DoubleTab& gravity_center_fa7=mesh.get_gravity_center_fa7();
  double X_p=gravity_center_fa7(fa7,0);
  double Y_p=gravity_center_fa7(fa7,1);
  double Z_p=gravity_center_fa7(fa7,2);

  sigma_xx_fa7_Stokes_th_(fa7)=(mu_f*vinf_Stokes/(particle_radius))*
                               (Z_p/(particle_radius*(1.+phi_mu)))*(pow(X_p/particle_radius,2)*
                                                                    (Z_p*Z_p/(pow(particle_radius,2)*(1.-
                                                                                                      pow(Z_p/particle_radius,2)))+3.*phi_mu-2.)+
                                                                    (pow(Y_p/particle_radius,2)*(1.-
                                                                                                 pow(Z_p/particle_radius,2))));
  sigma_xy_fa7_Stokes_th_(fa7)=(mu_f*vinf_Stokes/(particle_radius))*
                               (Z_p/(particle_radius*(1.+phi_mu)))*X_p*Y_p/(particle_radius*
                                                                            particle_radius)*(3.*phi_mu-3.);
  sigma_xz_fa7_Stokes_th_(fa7)=(mu_f*vinf_Stokes/(particle_radius))*
                               X_p/(particle_radius*(1.+phi_mu))*(-3./2.*phi_mu*(1.-
                                                                                 pow(Z_p/particle_radius,2))+3.*pow(Z_p/particle_radius,2)*
                                                                  (phi_mu/2.-1.));
  sigma_yy_fa7_Stokes_th_(fa7)=(mu_f*vinf_Stokes/(particle_radius))*
                               (Z_p/(particle_radius*(1.+phi_mu)))*(Y_p*Y_p/pow(particle_radius,2)*
                                                                    (3.*phi_mu-2.+(Z_p*Z_p/pow(particle_radius,2))/
                                                                     (1.-Z_p*Z_p/pow(particle_radius,2)))+
                                                                    X_p*X_p/(pow(particle_radius,2)*(1.-Z_p*Z_p/
                                                                                                     pow(particle_radius,2))));
  sigma_yz_fa7_Stokes_th_(fa7)=(mu_f*vinf_Stokes/(particle_radius))*
                               Y_p/(particle_radius*(1.+phi_mu))*(-3./2.*phi_mu*(1.-Z_p*
                                                                                 Z_p/pow(particle_radius,2))+3.*Z_p*Z_p/
                                                                  pow(particle_radius,2)*(phi_mu/2.-1.));
  sigma_zz_fa7_Stokes_th_(fa7)=(mu_f*vinf_Stokes/(particle_radius))*
                               (Z_p/(particle_radius*(1.+phi_mu)))*(-3.*phi_mu/2.*
                                                                    (1.-Z_p*Z_p/pow(particle_radius,2))+1.-3.*
                                                                    Z_p*Z_p/pow(particle_radius,2));
}

void Post_Processing_Hydrodynamic_Forces_Stokes::compute_friction_force_projected_tensor(int nb_fa7,
                                                                                         const Maillage_FT_Disc& mesh,
                                                                                         const IntVect& compo_connexes_fa7,
                                                                                         const ArrOfDouble& fa7_surface,
                                                                                         const DoubleTab& tab_fa7_normal)
{
  const Transport_Interfaces_FT_Disc& eq_transport = ptr_eq_transport_.valeur();
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  const Fluide_Diphasique& two_phase_fluid = eq_ns.fluide_diphasique();
  const Domaine_VDF& domain_vdf = ref_cast(Domaine_VDF, eq_ns.domaine_dis());
  const Domaine& domain = domain_vdf.domaine();
  const int id_fluid_phase=two_phase_fluid.get_id_fluid_phase();
  const Solid_Particle_base& solid_particle=ref_cast(Solid_Particle_base,
                                                     two_phase_fluid.fluide_phase(1-id_fluid_phase));
  const double particle_radius = solid_particle.get_equivalent_radius(); // On ne se sert de cette fonction que lorsqu'il y a une seule particule dans le domaine
  const double mu_f =two_phase_fluid.fluide_phase(id_fluid_phase).
                     viscosite_dynamique().valeurs()(0, 0);
  const double mu_p = two_phase_fluid.fluide_phase(1-id_fluid_phase).
                      masse_volumique().valeurs()(0, 0);
  const double phi_mu=mu_p/mu_f;
  const double& vinf_Stokes=get_vinf_Stokes();

  const DoubleTab& gravity_center_fa7=mesh.get_gravity_center_fa7();
  int interp_U_P1_ok=0;
  int interp_U_P2_ok=0;
  U_P1_.resize(nb_fa7, dimension);
  U_P2_.resize(nb_fa7, dimension);

  DoubleTab U_P1_spherique(nb_fa7, dimension);
  DoubleTab U_P2_spherique(nb_fa7, dimension);
  DoubleTab U_cg_spherique(nb_fa7, dimension);
  DoubleTab Urr(nb_fa7);
  DoubleTab Uthetar(nb_fa7);
  DoubleTab Uphir(nb_fa7);

  U_P1_=-1e15;
  U_P2_=-1e30;
  U_P1_spherique=-1e15;
  U_P2_spherique=-1e30;
  U_cg_spherique=-1e20;
  Urr=1e8;
  Uthetar=1e12;
  Uphir=1e15;
  const DoubleTab& positions_compo=eq_transport.get_particles_position();
  double theta=0;
  double phi=0;
  double distance_au_cg=0;
  Navier_Stokes_FT_Disc& eq_ns_non_const = ptr_eq_ns_.valeur();
  DoubleTab& velocity_Stokes_th=eq_ns_non_const.get_set_velocity_field_Stokes_th();
  if (location_stress_tensor_== Location_stress_tensor::ELEMENTS)
    {
      // 1. On calcule (interpolation trilineaire) la vitesse en P1 et P2 en coord cartesiennes
      interp_U_P1_ok=trilinear_interpolation_face(velocity_Stokes_th,
                                                  coord_neighbor_fluid_fa7_gradU_1_, U_P1_);
      interp_U_P2_ok=trilinear_interpolation_face(velocity_Stokes_th,
                                                  coord_neighbor_fluid_fa7_gradU_2_, U_P2_);
    }
  if ( interp_U_P1_ok && interp_U_P2_ok )
    {
      // 2. On passe en coordonnees spheriques
      for (int fa7=0; fa7<nb_fa7; fa7++)
        {
          int compo=compo_connexes_fa7(fa7);
          if (!mesh.facette_virtuelle(fa7))
            {
              DoubleVect distance_cg_vect(dimension);
              for (int i=0; i<dimension; i++)
                {
                  distance_cg_vect(i)=coord_neighbor_fluid_fa7_gradU_1_(fa7,i)-
                                      positions_compo(compo,i);
                }
              distance_au_cg=sqrt(local_carre_norme_vect(distance_cg_vect));

              if (fabs((coord_neighbor_fluid_fa7_gradU_1_(fa7,2)-positions_compo(compo,2))/
                       distance_au_cg)<=1)
                {
                  theta=acos((coord_neighbor_fluid_fa7_gradU_1_(fa7,2)-
                              positions_compo(compo,2))/distance_au_cg);
                }
              else if ((coord_neighbor_fluid_fa7_gradU_1_(fa7,2)-positions_compo(compo,2))/
                       distance_au_cg>1)
                {
                  theta=0;
                }
              else if ((coord_neighbor_fluid_fa7_gradU_1_(fa7,2)-positions_compo(compo,2))/
                       distance_au_cg<-1)
                {
                  theta=M_PI;
                }

              if ((coord_neighbor_fluid_fa7_gradU_1_(fa7,0)-positions_compo(compo,0))>0 &&
                  (coord_neighbor_fluid_fa7_gradU_1_(fa7,1)-positions_compo(compo,1))>=0)
                {
                  phi=atan((coord_neighbor_fluid_fa7_gradU_1_(fa7,1)-
                            positions_compo(compo,1))/(coord_neighbor_fluid_fa7_gradU_1_(fa7,0)
                                                       -positions_compo(compo,0)));
                }
              else if ((coord_neighbor_fluid_fa7_gradU_1_(fa7,0)-positions_compo(compo,0))>0
                       && (coord_neighbor_fluid_fa7_gradU_1_(fa7,1)-positions_compo(compo,1))<0)
                {
                  phi=atan((coord_neighbor_fluid_fa7_gradU_1_(fa7,1)-
                            positions_compo(compo,1))/(coord_neighbor_fluid_fa7_gradU_1_(fa7,0)
                                                       -positions_compo(compo,0)))+2*M_PI;
                }
              else if ((coord_neighbor_fluid_fa7_gradU_1_(fa7,0)-positions_compo(compo,0))<0)
                {
                  phi=atan((coord_neighbor_fluid_fa7_gradU_1_(fa7,1)-
                            positions_compo(compo,1))/(coord_neighbor_fluid_fa7_gradU_1_(fa7,0)
                                                       -positions_compo(compo,0)))+M_PI;
                }
              else if ((coord_neighbor_fluid_fa7_gradU_1_(fa7,0)-positions_compo(compo,0))==0
                       && (coord_neighbor_fluid_fa7_gradU_1_(fa7,1)-positions_compo(compo,1))>0)
                {
                  phi=M_PI/2.;
                }
              else if ((coord_neighbor_fluid_fa7_gradU_1_(fa7,0)-positions_compo(compo,0))==0
                       && (coord_neighbor_fluid_fa7_gradU_1_(fa7,1)-positions_compo(compo,1))<0)
                {
                  phi=3.*M_PI/2.;
                }

              U_P1_spherique(fa7,0)=sin(theta)*cos(phi)*U_P1_(fa7,0)+sin(theta)*sin(phi)*
                                    U_P1_(fa7,1)+cos(theta)*U_P1_(fa7,2);
              U_P1_spherique(fa7,1)=cos(theta)*cos(phi)*U_P1_(fa7,0)+cos(theta)*sin(phi)*
                                    U_P1_(fa7,1)-sin(theta)*U_P1_(fa7,2);
              U_P1_spherique(fa7,2)=-sin(phi)*U_P1_(fa7,0)+cos(phi)*U_P1_(fa7,1);

              U_P2_spherique(fa7,0)=sin(theta)*cos(phi)*U_P2_(fa7,0)+sin(theta)*sin(phi)*
                                    U_P2_(fa7,1)+cos(theta)*U_P2_(fa7,2);
              U_P2_spherique(fa7,1)=cos(theta)*cos(phi)*U_P2_(fa7,0)+cos(theta)*sin(phi)*
                                    U_P2_(fa7,1)-sin(theta)*U_P2_(fa7,2);
              U_P2_spherique(fa7,2)=-sin(phi)*U_P2_(fa7,0)+cos(phi)*U_P2_(fa7,1);

              U_cg_spherique(fa7,0)=0;//sin(theta)*cos(phi)*vitesses_compo(compo,0)+sin(theta)*sin(phi)*vitesses_compo(compo,1)+cos(theta)*vitesses_compo(compo,2);
              U_cg_spherique(fa7,1)=0;//cos(theta)*cos(phi)*vitesses_compo(compo,0)+cos(theta)*sin(phi)*vitesses_compo(compo,1)-sin(theta)*vitesses_compo(compo,2);
              U_cg_spherique(fa7,2)=0;//-sin(phi)*vitesses_compo(compo,0)+cos(phi)*vitesses_compo(compo,1);

              // On recalcule delta --> epsilon
              int elem_diph=domain.chercher_elements(gravity_center_fa7(fa7,0),
                                                     gravity_center_fa7(fa7,1),
                                                     gravity_center_fa7(fa7,2));
              DoubleVect delta_i(dimension);
              delta_i(0) = fabs(domain_vdf.dist_elem(elem_diph, domain_vdf.face_voisins(
                                                       domain_vdf.elem_faces(elem_diph, 0+dimension),1), 0));
              delta_i(1) = fabs(domain_vdf.dist_elem(elem_diph, domain_vdf.face_voisins(
                                                       domain_vdf.elem_faces(elem_diph, 1+dimension),1), 1));

              if (tab_fa7_normal(fa7,2)>0)
                {
                  delta_i(2) = fabs(domain_vdf.dist_elem(elem_diph, domain_vdf.face_voisins(
                                                           domain_vdf.elem_faces(elem_diph, 2+dimension),1), 2));
                }
              else
                delta_i(2) = fabs(domain_vdf.dist_elem(elem_diph, domain_vdf.face_voisins(
                                                         domain_vdf.elem_faces(elem_diph, 2),0), 2));
              double epsilon=0;
              for (int dim=0; dim<dimension; dim++)
                epsilon+= fabs(delta_i(dim)*fabs(tab_fa7_normal(fa7,dim))); // la distance d'interpolation varie en fonction du raffinement du maillage

              // On calcule les composantes de la force de frottements en coord spherique (apres simplifications) : ff=mu*(2*Urr, Uthetar, Uphir)
              Urr(fa7)=(-U_P2_spherique(fa7,0)+4.*U_P1_spherique(fa7,0)-
                        3.*U_cg_spherique(fa7,0))/(2.*epsilon);
              Uthetar(fa7)=(-U_P2_spherique(fa7,1)+4.*U_P1_spherique(fa7,1)-
                            3.*U_cg_spherique(fa7,1))/(2.*epsilon);
              Uphir(fa7)=(-U_P2_spherique(fa7,2)+4.*U_P1_spherique(fa7,2)-
                          3.*U_cg_spherique(fa7,2))/(2.*epsilon);

              friction_force_fa7_(fa7,0)=mu_f*fa7_surface(fa7)*(2.*sin(theta)*
                                                                cos(phi)*Urr(fa7)+cos(theta)*cos(phi)*Uthetar(fa7)-sin(phi)*Uphir(fa7));
              friction_force_fa7_(fa7,1)=mu_f*fa7_surface(fa7)*(2.*sin(theta)*
                                                                sin(phi)*Urr(fa7)+cos(theta)*sin(phi)*Uthetar(fa7)+cos(phi)*Uphir(fa7));
              friction_force_fa7_(fa7,2)=mu_f*fa7_surface(fa7)*(2.*cos(theta)*
                                                                Urr(fa7)-sin(theta)*Uthetar(fa7));

              friction_force_Stokes_th_fa7_(fa7,0)=-mu_f*vinf_Stokes*gravity_center_fa7(fa7,0)*
                                                   gravity_center_fa7(fa7,2)/(pow(particle_radius,3))*(9.*phi_mu+8.)/
                                                   (1.+phi_mu)*fa7_surface(fa7);
              friction_force_Stokes_th_fa7_(fa7,1)=-mu_f*vinf_Stokes*gravity_center_fa7(fa7,1)*
                                                   gravity_center_fa7(fa7,2)/(pow(particle_radius,3))*(9.*phi_mu+8.)/
                                                   (1.+phi_mu)*fa7_surface(fa7);
              friction_force_Stokes_th_fa7_(fa7,2)=mu_f*vinf_Stokes*1./(2.*particle_radius*(1.+phi_mu))*
                                                   (-3.*phi_mu+pow(gravity_center_fa7(fa7,2)/particle_radius,2)*(3.*phi_mu-4.))*fa7_surface(fa7);

              for (int dim=0; dim<dimension; dim++)
                {
                  total_friction_force_(compo,dim)+=friction_force_fa7_(fa7,dim);
                  total_friction_force_Stokes_th_(compo,dim)+=
                    friction_force_Stokes_th_fa7_(fa7,dim);
                }
            }
        }
    }
}
