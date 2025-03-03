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

#include <EFichier.h>
#include <EcritureLectureSpecial.h>
#include <Probleme_FT_Disc_gen.h>

#include <type_traits>
#include <Collision_Model_FT_base.h>

Implemente_base_sans_constructeur(Collision_Model_FT_base,"Collision_Model_FT_base",Objet_U);

Collision_Model_FT_base::Collision_Model_FT_base()
{
  fictive_wall_coordinates_.resize(2*dimension);
}
Entree& Collision_Model_FT_base::readOn(Entree& is)
{
  Param p(que_suis_je());
  set_param(p);
  p.lire_avec_accolades_depuis(is);
  return is;
}

Sortie& Collision_Model_FT_base::printOn(Sortie& os) const
{
  Cerr << "Error::printOn is not implemented." << finl;
  Process::exit();
  return os;
}

void Collision_Model_FT_base::set_param(Param& p)
{
  p.ajouter_non_std("collision_model", (this),Param::REQUIRED);
  p.ajouter_non_std("detection_method", (this),Param::REQUIRED);
  p.ajouter("collision_duration", &collision_duration_, Param::REQUIRED); // XD_ADD_P duration of the collision in seconds;
  p.ajouter("activate_collision_before_impact", &is_collision_activated_before_impact_, Param::REQUIRED);
  p.ajouter("activation_distance_percentage_diameter", &activation_distance_percentage_diameter_, Param::REQUIRED);
  p.ajouter("force_on_two_phase_elem", &is_force_on_two_phase_elem_, Param::REQUIRED);
}

int Collision_Model_FT_base::lire_motcle_non_standard(const Motcle& word, Entree& is)
{
  if (word=="collision_model")
    {
      Motcles words;
      words.add("hybrid_esi");
      words.add("breugem");
      Motcle secondword;
      is >> secondword;
      Cerr << "Reading collision_model attributes: " << secondword << finl;
      const int r = words.search(secondword);
      switch(r)
        {
        case 0:
          collision_model_ = Collision_model::HYBRID_ESI;
          break;
        case 1:
          collision_model_ = Collision_model::BREUGEM;
          break;
        default:
          Cerr << "Error " << words << "was expected whereas " << secondword <<
               " has been found."<< finl;
          Process::exit();
        }
      return 1;
    }
  else if (word=="detection_method")
    {

      Motcle openbrace ="{";
      Motcle closedbrace="}";
      Motcle secondword;
      is >> secondword;
      if (secondword==openbrace)
        {
          is >> secondword;

          Motcles words;
          words.add("check_all");
          words.add("verlet");
          words.add("lc_verlet");
          const int r = words.search(secondword);
          switch(r)
            {
            case 0:
              detection_method_ = Detection_method::CHECK_ALL;
              break;
            case 1:
              detection_method_ = Detection_method::VERLET;
              break;
            case 2:
              detection_method_ = Detection_method::LC_VERLET;
              break;
            default:
              Cerr << "Error " << words << "was expected whereas "
                   << secondword << " has been found."<< finl;
              Process::exit();
            }

          is >> secondword;
          while (secondword != closedbrace)
            {
              Motcles otherwords;
              otherwords.add("detection_thickness_Verlet");
              otherwords.add("nb_pas_dt_max_Verlet");
              int rang2 = otherwords.search(secondword);
              switch(rang2)
                {
                case 0:
                  is >> detection_thickness_Verlet_;
                  break;
                case 1:
                  is >> nb_dt_max_Verlet_;
                  break;
                default:
                  Cerr << "Collision_Model_FT_base::lire_motcle_non_standard\n"
                       << " options of collision_detection are:\n"
                       << otherwords;
                  Process::exit();
                }
              is >> secondword;
            }
        }
      return 1;
    }
  else
    {
      Cerr << word << " is not a keyword understood by " << que_suis_je() <<
           " in lire_motcle_non_standard"<< finl;
      Process::exit();
    }
  return -1;
}

void Collision_Model_FT_base::reset()
{
  const int nb_boundaries=2*dimension;
  F_old_.resize(nb_particles_tot_,nb_particles_tot_+nb_boundaries);
  F_now_.resize(nb_particles_tot_,nb_particles_tot_+nb_boundaries);
  e_eff_.resize(nb_particles_tot_,nb_particles_tot_+nb_boundaries);
  Cerr << "WARNING: Collision_Model_FT_base::reset of F_old_, F_now_, "
       "lagrangian_contact_forces_ and e_eff_." << finl;
}

int Collision_Model_FT_base::reprendre(Entree& is)
{
  Nom readword;
  const int format_xyz = EcritureLectureSpecial::is_lecture_special();

  is >> readword;
  if (readword != que_suis_je())
    {
      Cerr << "Error in Collision_Model_FT_base::reprendre\n";
      Cerr << "We was expecting " << que_suis_je();
      Cerr << "\n We found " << readword << finl;
      Process::exit();
    }
  is >> nb_particles_tot_;
  e_eff_.resize(nb_particles_tot_,nb_particles_tot_+2*dimension);
  e_eff_.resize(nb_particles_tot_,nb_particles_tot_+2*dimension);

  if (format_xyz)
    {
      for (int i=0; i<F_old_.dimension(0); i++)
        for (int j=0; j<F_old_.dimension(1); j++)
          is>>F_old_(i,j);
      for (int i=0; i<e_eff_.dimension(0); i++)
        for (int j=0; j<e_eff_.dimension(1); j++)
          is>>e_eff_(i,j);
    }
  else
    {
      is>>F_old_;
      is>>e_eff_;
    }
  return 1;
}

int Collision_Model_FT_base::sauvegarder(Sortie& os) const
{
  int special, afaire;
  const int format_xyz = EcritureLectureSpecial::is_ecriture_special(special, afaire);
  int bytes = 0;
  if (format_xyz)
    {
      if (Process::je_suis_maitre())
        {
          os << Nom(que_suis_je()) << finl;
          os << nb_particles_tot_;
          for (int i=0; i<F_old_.dimension(0); i++)
            for (int j=0; j<F_old_.dimension(1); j++)
              os<<F_old_(i,j);
          for (int i=0; i<e_eff_.dimension(0); i++)
            for (int j=0; j<e_eff_.dimension(1); j++)
              os<<e_eff_(i,j);
        }
      return 0;
    }
  else
    {
      os << que_suis_je() << finl;
      os << nb_particles_tot_;
      os << F_old_;
      bytes += 8 * F_old_.size_array();
      os << e_eff_;
      bytes += 8 * e_eff_.size_array();
      return bytes;
    }
  return 0;
}

void Collision_Model_FT_base::resize_geometric_parameters()
{
  domain_dimensions_.resize(dimension);
  nb_nodes_.resize(dimension);
  origin_.resize(dimension);
}

void Collision_Model_FT_base::set_spring_properties(const Solid_Particle_base& solid_particle)
{
  if (collision_duration_>0)
    {
      const double mass_sphere = solid_particle.get_mass();
      const double e_dry = solid_particle.get_e_dry();
      const double mass_eff_part_part = mass_sphere/2;
      const double mass_eff_wall_part = mass_sphere;

      stiffness_breugem_part_part_ = compute_stiffness_breugem(mass_eff_part_part,e_dry);
      stiffness_breugem_wall_part_ = compute_stiffness_breugem(mass_eff_wall_part,e_dry);

      if (collision_model_ == Collision_model::BREUGEM)
        {
          damper_breugem_part_part_ = compute_damper_breugem(mass_eff_part_part,e_dry);
          damper_breugem_wall_part_ = compute_damper_breugem(mass_eff_wall_part,e_dry);
        }
    }
}

void Collision_Model_FT_base::associate_transport_equation(const Equation_base& equation)
{
  const Transport_Interfaces_FT_Disc& eq = ref_cast(Transport_Interfaces_FT_Disc,equation);
  refequation_transport_ = eq;
}

void Collision_Model_FT_base::compute_fictive_wall_coordinates(const double& radius)
{
  DoubleVect offset_values(2*dimension);

  switch(is_collision_activated_before_impact_)
    {
    case 0:
      offset_values=0;
      break;
    case 1:
      if (activation_distance_>0) offset_values=activation_distance_;
      break;
    default:
      Cerr << "Collision_Model_FT_base::compute_fictive_wall_coordinates error"  <<finl;
      Process::exit();
      break;
    }

  // an offset is computed for the walls to activate the collision process before the impact
  fictive_wall_coordinates_(0) = origin_(0) - radius + offset_values(0);
  fictive_wall_coordinates_(1) = origin_(1) - radius + offset_values(1);
  fictive_wall_coordinates_(2) = origin_(2) - radius + offset_values(2);
  fictive_wall_coordinates_(3) = origin_(0) + domain_dimensions_(0) + radius - offset_values(3);
  fictive_wall_coordinates_(4) = origin_(1) + domain_dimensions_(1) + radius - offset_values(4);
  fictive_wall_coordinates_(5) = origin_(2) + domain_dimensions_(2) + radius - offset_values(5);
}

/*! @brief Recover the geometric parameters of the domain:
 * number of nodes in each direction
 * origin of the domain
 * dimensions of the domain
 */
void Collision_Model_FT_base::set_geometric_parameters(const Domaine_VDF& domaine_vdf)
{
  const Domaine& domain = domaine_vdf.domaine();
  const DoubleTab BB=domain.getBoundingBox();
  const Bords& bords=domain.faces_bord();
  DoubleVect NiNj(dimension); // NiNj=(NyNz NxNz NxNy )
  resize_geometric_parameters();

  // 1. Number of nodes per direction
  NiNj=0;
  for (int i=0; i<bords.nb_bords(); i++)
    {
      int nb_boundary_faces = mp_sum(ref_cast(Frontiere,bords(i)).nb_faces());
      int nb_boundary_faces_local=ref_cast(Frontiere,bords(i)).nb_faces();
      if (nb_boundary_faces_local>0)
        {
          int face1=ref_cast(Frontiere,bords(i)).num_premiere_face();
          int orientation_face1=domaine_vdf.orientation(face1);
          NiNj(orientation_face1)=nb_boundary_faces;
        }
    }
  long long NxNy= static_cast<int>(mp_max(NiNj(2))) ;
  long long NxNz= static_cast<int>(mp_max(NiNj(1)));
  long long NyNz= static_cast<int>(mp_max(NiNj(0)));

  int Nx,Ny,Nz;

  Ny= NxNz>0 ? static_cast<int>(sqrt(NxNy*NyNz/NxNz)) : 0; // nb elem in the y-direction
  Nz= NxNy>0 ? static_cast<int>((NxNz*Ny/NxNy)) : 0; // nb elem in the z-direction, WARNING: operations order matter because Nz is an integer
  Nx= Ny>0 ? static_cast<int>(NxNy/Ny) : 0; // nb elem in the x-direction

  nb_nodes_(0)=Nx++; // nb nodes in the x-direction
  nb_nodes_(1)=Ny++; // nb nodes in the y-direction
  nb_nodes_(2)=Nz++; // nb nodes in the z-direction

  // 2. Origin and Domain_dimensions
  DoubleVect Origin(dimension);
  DoubleVect Domain_dimensions(dimension);
  Origin=0.;
  Domain_dimensions=0.;

  for (int j=0; j<dimension; j++)
    {
      double min_ = mp_min(BB(j,0));
      double max_ = mp_max(BB(j,1));

      Origin(j)=min_;
      Domain_dimensions(j)=max_-min_;
    }

  set_origin(Origin);
  set_domain_dimensions(Domain_dimensions);

  Cerr << "Origin " << Origin << finl;
  Cerr << "Domain length" << Domain_dimensions << finl;
  Cerr << "Nx Ny Nz " << nb_nodes_ << finl;
}

/*! @brief Check if two particles have the same ID
 * Very important function to stop computation if
 * two particles coalesce.
 */
int Collision_Model_FT_base::check_for_duplicates(ArrOfInt& vector)
{
  int flag =0;
  ArrOfInt copy_vector(vector);
  const int size = copy_vector.size_array();
  copy_vector.ordonne_array();
  for (int i = 0; i < size-1; i++)
    {
      if (copy_vector(i)==copy_vector(i+1))
        {
          flag = 1;
          Cerr << copy_vector(i) << " is duplicate !!" << finl ;
        }
    }
  return flag;
}


/* @brief we send the list of real particles to proc from lower zones
 */
ArrOfInt send_receive_list_id_particles(ArrOfInt& list_id_real_particles_to_send,
                                        const ArrOfInt& list_lower_zone,
                                        const Schema_Comm_FT& schema_comm)
{
  ArrOfInt list_id_particle_recv;
  list_id_particle_recv.resize_array(0);
  int nb_elem_recv=0;
  const int nb_real_particles_to_send=list_id_real_particles_to_send.size_array();

  schema_comm.begin_comm();

  for (int ind_pe_dest=0; ind_pe_dest<list_lower_zone.size_array(); ind_pe_dest++)
    {
      int PE_dest=list_lower_zone(ind_pe_dest);
      for(int ind_particle=0; ind_particle<nb_real_particles_to_send; ind_particle++)
        {
          const int id_particle=list_id_real_particles_to_send(ind_particle);
          assert(PE_dest!=Process::me());
          schema_comm.send_buffer(PE_dest) << id_particle;
        }
    }

  schema_comm.echange_taille_et_messages();

  const ArrOfInt& recv_pe_list = schema_comm.get_recv_pe_list();
  const int nb_recv_pe = recv_pe_list.size_array();
  for (int i=0; i<nb_recv_pe; i++)
    {
      const int pe_source = recv_pe_list[i];
      Entree& buffer = schema_comm.recv_buffer(pe_source);
      while(1)
        {
          int id_particle_recv=-1;
          buffer >> id_particle_recv;
          if (buffer.eof())
            break;
          if (id_particle_recv<0)
            Process::exit();
          nb_elem_recv++;
          list_id_particle_recv.append_array(id_particle_recv);
        }
    }
  schema_comm.end_comm();

  return list_id_particle_recv;
}


void Collision_Model_FT_base::research_collision_pairs_Verlet(const Navier_Stokes_FT_Disc&
                                                              eq_ns,const Transport_Interfaces_FT_Disc& eq_transport)
{
  ArrOfInt list_particles_to_check_LC;
  const ArrOfInt& gravity_center_elem=eq_transport.get_gravity_center_elem();
  const Domaine_VF& domain_vf=ref_cast(Domaine_VF,eq_ns.domaine_dis());
  const int& nb_elem=domain_vf.nb_elem();
  const Maillage_FT_Disc& mesh=eq_transport.maillage_interface();
  const Schema_Comm_FT& schema_comm=mesh.get_schema_comm_FT();
  const DoubleTab& particles_velocity=eq_transport.get_particles_velocity();
  const DoubleTab& particles_position=eq_transport.get_particles_position();
  const Fluide_Diphasique& two_phase_elem=eq_ns.fluide_diphasique();
  const int id_solid=1-two_phase_elem.get_id_fluid_phase();
  const auto& solid_particle=ref_cast(Solid_Particle_base,two_phase_elem.fluide_phase(id_solid));
  const double& radius=solid_particle.get_equivalent_radius();
  const Schema_Temps_base& time_scheme=eq_ns.schema_temps();
  const double& time_step=time_scheme.pas_de_temps();

  Cerr << "Update of Verlet tables - " ;
  // Step 1 : update of Linked Cells (for more info on the method, see X. Fang et al.,
  // J. Sound Vib., 2007).
  // Here, each linked cell match a CPU domain. We seed the id of real particles to LC
  // from lower zones (to receive LC from upper zones)
  if (detection_method_==Detection_method::LC_VERLET)
    {
      nb_real_particles_=0;
      list_real_particles_.resize_array(0);
      for (int particle=0; particle<nb_particles_tot_; particle++)
        {
          if (gravity_center_elem(particle)>-1 && gravity_center_elem(particle)<nb_elem)
            {
              list_real_particles_.append_array(particle);
              nb_real_particles_++;
            }
        }
      list_real_particles_.resize_array(nb_real_particles_);

      //Process::barrier(); // we wait that all procs have updated their list to do the exchange
      ArrOfInt list_virtual_particles;
      const ArrOfInt& list_lower_zone=get_list_lower_zone();
      list_virtual_particles=send_receive_list_id_particles(list_real_particles_,
                                                            list_lower_zone,
                                                            schema_comm);
      if (nb_real_particles_>0)
        {
          Cerr << "list_real_particles_ " << list_real_particles_<< finl;
          Cerr << "list_virtual_particles " << list_virtual_particles<< finl;
          list_particles_to_check_LC.resize(list_real_particles_.size_array()+
                                            list_virtual_particles.size_array());
          for (int k=0; k<list_real_particles_.size_array()
               +list_virtual_particles.size_array(); k++)
            {
              if (k<list_real_particles_.size_array()) list_particles_to_check_LC(k)=
                  list_real_particles_(k);
              else list_particles_to_check_LC(k)=
                  list_virtual_particles(k-list_real_particles_.size_array());
            }
          Cerr << "list_particles_to_check_LC " << list_particles_to_check_LC<< finl; // do not sort the list
        }
    }
  // ETAPE B : update of Verlet tables
  // we compute the distance between particles and we save the pairs that are within
  // the distance detection_thickness_Verlet_. I usually define it as 30% of the particle
  // diameter in the data file.
  // For Verlet table without LC, we test all particle pairs
  // For Verlet method with LC, for all real particles, we test all particles
  // pairs from the upper zone (proc domain + upper procs domain)
  double max_vi_glob=0;
  double max_vi=0.;
  if (nb_real_particles_>0)
    {
      Verlet_tables_=0;
      Verlet_tables_.dimensionner(nb_real_particles_);

      compute_Verlet_tables(particles_position,
                            particles_velocity,
                            max_vi,
                            radius,
                            list_particles_to_check_LC);
    }
  max_vi_glob=mp_max(max_vi); // with the linked cell method, we compute the max.
  // With only Verlet method,
  // there is no need to, but max_vi=mp_max(max_vi) for each proc.
  //  Computing the next Verlet table update time step based on the maximum velocity
  // of a particle in the domain
  int deltat_compute_min = (max_vi_glob>0) ?
                           static_cast<int>(floor((detection_thickness_Verlet_/(2*max_vi_glob))/
                                                  time_step)) :
                           nb_dt_max_Verlet_;

  nb_dt_compute_Verlet_=std::min(deltat_compute_min,nb_dt_max_Verlet_);
  Cerr << "Next update in " <<nb_dt_compute_Verlet_ << " time steps." << finl;
  nb_dt_Verlet_=0;
  // Cerr to check if Verlet tables are correctly filled
  if (nb_particles_tot_<20)
    {
      for (int ind_part_i=0; ind_part_i <nb_real_particles_; ind_part_i++)
        {
          Cerr << "part " << ind_part_i << " : ";
          Cerr << "Verlet_tables_[part_i].size() " << Verlet_tables_[ind_part_i].size() << " -- ";
          for (int ind_part_j=0; ind_part_j<Verlet_tables_[ind_part_i].size(); ind_part_j++)
            Cerr << " " << Verlet_tables_[ind_part_i][ind_part_j];
          Cerr << finl;
        }
    }
}

int Collision_Model_FT_base::get_last_id(const ArrOfInt& list_particles_to_check_LC)
{
  if (detection_method_==Detection_method::CHECK_ALL)
    Process::exit("Collision_Model_FT_base::get_last_id Should not be in here!");
  else if (detection_method_==Detection_method::VERLET ||
           detection_method_==Detection_method::LC_VERLET)
    return(list_particles_to_check_LC.size_array());
  return 0;
}


int Collision_Model_FT_base::get_id(const ArrOfInt& list_particle, const int ind_id_particle)
{
  if (detection_method_==Detection_method::CHECK_ALL)
    Process::exit("Collision_Model_FT_base::get_id Should not be in here!");
  else if (detection_method_==Detection_method::VERLET)
    return(ind_id_particle);
  else if (detection_method_==Detection_method::LC_VERLET)
    return(list_particle(ind_id_particle));
  else
    Process::exit("Collision_Model_FT_base::get_id unkwnown detection_method_!");
  return 0;
}

int Collision_Model_FT_base::get_particle_j(const int ind_particle_i, const int ind_particle_j)
{
  if (detection_method_==Detection_method::CHECK_ALL)
    return(ind_particle_j);
  else if (detection_method_==Detection_method::VERLET ||
           detection_method_==Detection_method::LC_VERLET)
    return (Verlet_tables_[ind_particle_i][ind_particle_j]);
  return 0;
}

int Collision_Model_FT_base::get_nb_particles_j(const int ind_particle_i) const
{
  if (detection_method_==Detection_method::CHECK_ALL)
    return (nb_particles_tot_+2*dimension);
  else if (detection_method_==Detection_method::VERLET ||
           detection_method_==Detection_method::LC_VERLET)
    return (Verlet_tables_[ind_particle_i].size());
  return 0;
}

int Collision_Model_FT_base::get_particle_i(const int ind_particle_i)
{
  return(detection_method_==Detection_method::LC_VERLET ?
         list_real_particles_(ind_particle_i) : ind_particle_i);
}

int Collision_Model_FT_base::get_ind_start_particles_j(const int ind_particle_i) const
{
  return(detection_method_==Detection_method::CHECK_ALL ?
         (ind_particle_i+1) : 0);
}


void Collision_Model_FT_base::compute_Verlet_tables(const DoubleTab& particles_position,
                                                    const DoubleTab& particles_velocity,
                                                    double& max_vi,
                                                    const double& radius,
                                                    const ArrOfInt& list_particles_to_check_LC)
{
  int last_id = get_last_id(list_particles_to_check_LC);

  for (int ind_id_particle_i=0; ind_id_particle_i< nb_real_particles_; ind_id_particle_i++)
    {
      int id_particle_i=get_id(list_particles_to_check_LC,ind_id_particle_i);
      const double vi=fabs(sqrt(
                             pow(particles_velocity(id_particle_i,0),2) +
                             pow(particles_velocity(id_particle_i,1),2) +
                             pow(particles_velocity(id_particle_i,2),2)));
      if (vi>max_vi)
        max_vi=vi;
      // particle-particle distance
      for (int ind_id_particle_j=ind_id_particle_i+1; ind_id_particle_j< last_id; ind_id_particle_j++)
        {
          int id_particle_j= get_id(list_particles_to_check_LC, ind_id_particle_j);
          double dij=sqrt(
                       pow(particles_position(id_particle_j,0)-particles_position(id_particle_i,0),2) +
                       pow(particles_position(id_particle_j,1)-particles_position(id_particle_i,1),2) +
                       pow(particles_position(id_particle_j,2)-particles_position(id_particle_i,2),2));

          if ((dij-2*radius)<=detection_thickness_Verlet_)
            Verlet_tables_[ind_id_particle_i].add(id_particle_j);

        }

      // wall-particle distance
      for (int ind_wall=0; ind_wall<2*dimension; ind_wall++)
        {
          const int ori = ind_wall < dimension ? ind_wall : ind_wall - dimension;
          const double dij=fabs(fictive_wall_coordinates_(ind_wall)-
                                particles_position(id_particle_i,ori));

          if ((dij-2*radius)<=detection_thickness_Verlet_)
            Verlet_tables_[ind_id_particle_i].add(nb_particles_tot_+ind_wall);
        }
    }
}


void Collision_Model_FT_base::compute_dX_dU(DoubleTab& dX,
                                            DoubleTab& dU,
                                            const int& particle,
                                            const int& neighbor,
                                            const DoubleTab& particles_position,
                                            const DoubleTab& particles_velocity,
                                            const bool is_particle_particle_collision)
{
  if (is_particle_particle_collision)
    {
      for (int d = 0; d < dimension; d++)
        {
          dX(d) = particles_position(particle, d) - particles_position(neighbor, d);
          dU(d) = particles_velocity(particle, d) - particles_velocity(neighbor, d);
        }
    }
  else
    {
      int ind_wall = neighbor - nb_particles_tot_;
      int ori = ind_wall < dimension ? ind_wall : ind_wall - dimension;
      dX(ori) = particles_position(particle, ori) - fictive_wall_coordinates_(ind_wall);
      for (int d = 0; d < dimension; d++)
        dU(d) = particles_velocity(particle, d);
    }
}

DoubleTab Collision_Model_FT_base::compute_contact_force(
  const double& next_dist_int,
  const DoubleTab& norm,
  const DoubleTab& dUn,
  const int& particle_i,
  const int& particle_j,
  const int& is_compression_step,
  const double& is_collision_part_part)
{
  DoubleTab force_contact(dimension);

  if (collision_model_==Collision_model::HYBRID_ESI) // see Hamidi et al., IJMF, 2023.
    {
      const double e_eff_particle = is_compression_step ? 1 : e_eff_(particle_i, particle_j);
      const double stiffness = is_collision_part_part ? stiffness_breugem_part_part_:
                               stiffness_breugem_wall_part_;
      for (int d = 0; d < dimension; d++)
        force_contact(d)= -pow(e_eff_particle,2) * stiffness * next_dist_int * norm(d);
    }
  else if (collision_model_==Collision_model::BREUGEM) // See. W-P. Breugem, 2010.
    {
      const double stiffness = is_collision_part_part ? stiffness_breugem_part_part_:
                               stiffness_breugem_wall_part_;
      const double damper = is_collision_part_part ? damper_breugem_part_part_:
                            damper_breugem_wall_part_;
      for (int d = 0; d < dimension; d++)
        force_contact(d)= -stiffness*next_dist_int*norm(d) -damper*dUn(d);
    }
  else
    Process::exit("Collision_Model_FT_base::compute_contact_force unknown collision_model_");

  return force_contact;
}

void Collision_Model_FT_base::discretize_contact_forces_eulerian_field(
  const DoubleTab& volumic_phase_indicator_function,
  const Domaine_VF& domain_vf,
  const IntTab& particles_eulerian_id_number,
  DoubleTab& contact_force_source_term)
{
  Process::exit("Collision_Model_FT_base::discretize_contact_forces_eulerian_field "
                "not coded for particles of arbitrary shape");
}

bool Collision_Model_FT_base::is_Verlet_activated()
{
  if (detection_method_ == Detection_method::CHECK_ALL)
    return false;
  return true;
}

bool Collision_Model_FT_base::is_LC_activated()
{
  if (detection_method_ == Detection_method::LC_VERLET)
    return true;
  return false;
}

void Collision_Model_FT_base::set_LC_zones(const Domaine_VF& domain_vf, const Schema_Comm_FT& schema_com)
{
  int nsup=0;
  int ninf=0;

  int nb_joints=domain_vf.nb_joints();

  double x0 = domain_vf.xp(0,0);
  double y0 = domain_vf.xp(0,1);
  double z0 = domain_vf.xp(0,2);
  double epsilon=1e-10;
  const int nb_elems=domain_vf.nb_elem();
  const int nb_elems_tot=domain_vf.nb_elem_tot();

  Cerr << "number or reals elements " << nb_elems << finl;
  Cerr << "total number of elements " << nb_elems_tot << finl;
  Cerr << "(x0 y0 z0) = (" << x0 << " " << y0 << " " << z0 << ")" << finl;
  Cerr << "joints number " << nb_joints << finl;

  // Firstly, every proc send coordinates of its first elem to other procs
  DoubleTabFT list_coord_recv(0,dimension);
  IntTabFT list_pe_recv(0);
  Cerr << "list_coord_recv.dimensions " << list_coord_recv.dimension(0) << " "
       << list_coord_recv.dimension(1) << finl;
  int nb_elem_recv=0;
  schema_com.begin_comm();
  for (int ind_pe_dest=0; ind_pe_dest<nb_joints; ind_pe_dest++)
    {
      const Joint& joint_temp = domain_vf.joint(ind_pe_dest);
      const int pe_dest = joint_temp.PEvoisin();
      assert(pe_dest!=Process::me());
      schema_com.send_buffer(pe_dest)  << x0 << y0 << z0;
    }
  schema_com.echange_taille_et_messages();
  const ArrOfInt& recv_pe_list = schema_com.get_recv_pe_list();
  const int nb_recv_pe = recv_pe_list.size_array();
  for (int i=0; i<nb_recv_pe; i++)
    {
      const int pe_source = recv_pe_list[i];
      Entree& buffer = schema_com.recv_buffer(pe_source);
      while(1)
        {
          double x0_rcv=-1, y0_rcv=-1, z0_rcv=-1;
          buffer >> x0_rcv >> y0_rcv >> z0_rcv;
          if (buffer.eof())
            break;
          nb_elem_recv++;
          list_coord_recv.append_line(x0_rcv,y0_rcv,z0_rcv);
          list_pe_recv.append_line(pe_source);
        }
    }
  schema_com.end_comm();
  list_coord_recv.resize(nb_elem_recv,dimension);
  list_pe_recv.resize(nb_elem_recv);

  if (nb_joints != nb_elem_recv) Process::exit("EB : Collision_Model_FT_base::"
                                                 "set_LC_zones -- erreur identification du 1er element"
                                                 " des procs de la zone de joint.");
  for(int ind_pe=0; ind_pe<nb_elem_recv; ind_pe++)
    {
      const int pe_voisin = list_pe_recv(ind_pe);
      double x0_joint =  list_coord_recv(ind_pe,0);
      double y0_joint =  list_coord_recv(ind_pe,1);
      double z0_joint =  list_coord_recv(ind_pe,2);

      if ( y0_joint>y0 || (fabs(y0_joint-y0)<epsilon && x0_joint>x0) ||
           ( fabs(x0_joint-x0) < epsilon && fabs(y0_joint-y0)<epsilon && (z0_joint>z0)) )
        list_upper_zone_.append_array(pe_voisin), nsup++;
      if ( y0_joint<=y0 && (fabs(y0_joint-y0)>=epsilon || x0_joint<=x0)
           && ( fabs(x0_joint-x0) >= epsilon || fabs(y0_joint-y0)>=epsilon
                || (z0_joint<=z0)) )
        list_lower_zone_.append_array(pe_voisin), ninf++;

    }

  list_upper_zone_.resize_array(nsup);
  list_lower_zone_.resize_array(ninf);

  Cerr << "list_upper_zone_ " << list_upper_zone_ << finl;
  Cerr << "list_lower_zone_ " << list_lower_zone_ << finl;
  Cerr << "max number of upper zone  = " << mp_max(nsup) << finl;
  Cerr << "max number of lower zone  = " << mp_max(ninf) << finl;
}


