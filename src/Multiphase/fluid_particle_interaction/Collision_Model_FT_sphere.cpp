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
*****************************************************************************/

#include <Collision_Model_FT_sphere.h>
#include <Probleme_FT_Disc_gen.h>
#include <Solid_Particle_sphere.h>


Implemente_instanciable_sans_constructeur(Collision_Model_FT_sphere,"Collision_Model_FT_sphere",Collision_Model_FT_base);

Collision_Model_FT_sphere::Collision_Model_FT_sphere()
{
}

Entree& Collision_Model_FT_sphere::readOn (Entree& is)
{
  Collision_Model_FT_base::readOn(is);
  return is;
}

Sortie& Collision_Model_FT_sphere::printOn(Sortie& os) const
{
  Cerr << "Error::printOn is not implemented." << finl;
  Process::exit();
  return os;
}


void Collision_Model_FT_sphere::compute_lagrangian_contact_forces(const Fluide_Diphasique& two_phase_fluid,
                                                                  const DoubleTab& particles_position,
                                                                  const DoubleTab& particles_velocity,
                                                                  const double& deltat_simu)
{
  const int& id_fluid_phase= two_phase_fluid.get_id_fluid_phase();
  const int& id_solid_phase=1-id_fluid_phase;
  const auto& solid_particle=ref_cast(Solid_Particle_sphere,two_phase_fluid.fluide_phase(id_solid_phase));
  const auto& incompressible_fluid=ref_cast(Fluide_Incompressible,
                                            two_phase_fluid.fluide_phase(id_fluid_phase));
  const double& solid_density = solid_particle.masse_volumique().valeurs()(0, 0);
  const double& fluid_density = incompressible_fluid.masse_volumique().valeurs()(0, 0);
  const double& fluid_viscosity  = fluid_density
                                   * incompressible_fluid.viscosite_cinematique().valeurs()(0, 0);
  const double& radius_sphere=solid_particle.get_radius();
  const double& volume_sphere=solid_particle.get_volume();
  const double& e_dry=solid_particle.get_e_dry();
  const double min_threshold=1e-10;
  DoubleTab dX(dimension), dU(dimension);
  lagrangian_contact_forces_=0;
  for (int ind_particle_i = 0; ind_particle_i < nb_real_particles_; ind_particle_i++)
    {
      int particle_i=get_particle_i(ind_particle_i);
      int nb_particles_j=get_nb_particles_j(ind_particle_i);
      int ind_start_part_j=get_ind_start_particles_j(ind_particle_i);
      for (int ind_particle_j =ind_start_part_j; ind_particle_j < nb_particles_j; ind_particle_j++)
        {
          dX = 0;
          dU = 0;
          int particle_j=get_particle_j(ind_particle_i,ind_particle_j);
          int is_particle_particle_collision = particle_j < nb_particles_tot_;
          compute_dX_dU(dX, dU, particle_i, particle_j, particles_position,
                        particles_velocity, is_particle_particle_collision);
          double dist_gravity_center = sqrt(local_carre_norme_vect(dX));
          double dist_between_particles = dist_gravity_center - 2 * radius_sphere;

          F_now_(particle_i, particle_j) = 0;
          if (dist_between_particles <= 0) // contact
            {
              DoubleTab norm(dimension);
              for (int d = 0; d < dimension; d++)
                norm(d) = dX(d) / dist_gravity_center;

              double dX_scal_dU = local_prodscal(dX,dU);
              DoubleTab dUn(dimension);
              for (int d = 0; d < dimension; d++)
                dUn(d) = (dX_scal_dU / dist_gravity_center) * norm(d);

              const double impact_velocity = sqrt(local_carre_norme_vect(dUn));

              F_now_(particle_i, particle_j) = 1;
              int is_start_of_collision = F_now_(particle_i, particle_j) >
                                          F_old_(particle_i, particle_j); // We need to know
              // if this is the first time step of the collision to compute the impact velocity

              DoubleTab next_dX(dimension);
              for (int d = 0; d < dimension; d++)
                next_dX(d) = dX(d) + deltat_simu * dU(d);
              const double next_dist_gravity_center = sqrt(local_carre_norme_vect(next_dX));
              const double next_dist_between_particles = next_dist_gravity_center - 2 * radius_sphere;
              const double effective_radius = is_particle_particle_collision ? radius_sphere/2 :
                                              radius_sphere;
              const double impact_Stokes = solid_density * 2 * effective_radius * impact_velocity /
                                           (9 * fluid_viscosity);
              if (is_start_of_collision)
                e_eff_(particle_i,particle_j)=e_dry*compute_ewet_legendre(impact_Stokes);
              DoubleTab force_contact=compute_contact_force(
                                        next_dist_between_particles,
                                        norm,
                                        dUn,
                                        particle_i,
                                        particle_j,
                                        dX_scal_dU<=0,
                                        is_particle_particle_collision);

              for (int d = 0; d < dimension; d++)
                {
                  lagrangian_contact_forces_(particle_i, d) += fabs(force_contact(d)) <=
                                                               min_threshold ? 0 : force_contact(d) / volume_sphere;
                  if (!is_particle_particle_collision)
                    continue; // wall collision, no force to apply on the wall
                  lagrangian_contact_forces_(particle_j, d) -= fabs(force_contact(d)) <=
                                                               min_threshold ? 0 :  force_contact(d) / volume_sphere;
                }
            }
          F_old_(particle_i, particle_j) = F_now_(particle_i, particle_j);
        }
    }

  if (detection_method_==Detection_method::LC_VERLET)
    {
      mp_sum_for_each_item(lagrangian_contact_forces_);
      mp_max_for_each_item(F_old_);
      mp_max_for_each_item(F_now_);
      mp_max_for_each_item(e_eff_);
    }
}

void Collision_Model_FT_sphere::discretize_contact_forces_eulerian_field(
  const DoubleTab& volumic_phase_indicator_function,
  const Domaine_VF& domain_vf,
  const IntTab& particles_eulerian_id_number,
  DoubleTab& contact_force_source_term)
{
  const DoubleVect& interlaced_volumes=domain_vf.volumes_entrelaces();
  const int nb_faces=interlaced_volumes.size_array();
  const IntVect& orientation = domain_vf.orientation();
  const IntTab& face_voisins=domain_vf.face_voisins();
  for (int face=0; face<nb_faces; face++)
    {
      const int left_elem=face_voisins(face,0);
      const int right_elem=face_voisins(face,1);
      int id_left = left_elem != -1 ? particles_eulerian_id_number(left_elem) : -1;
      int id_right = right_elem != -1 ? particles_eulerian_id_number(right_elem) : -1;
      const int id_number=std::max(id_left,id_right);
      if (id_number!=-1)
        {
          const int ori=orientation(face);
          contact_force_source_term(face)=(1-volumic_phase_indicator_function(face))
                                          *interlaced_volumes(face)*lagrangian_contact_forces_(id_number,ori);
        }
    }
}

