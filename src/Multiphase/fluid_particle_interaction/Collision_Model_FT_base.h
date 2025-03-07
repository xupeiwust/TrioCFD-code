
#ifndef COLLISION_MODEL_FT_BASE_included
#define COLLISION_MODEL_FT_BASE_included


#include <TRUSTTabFT_forward.h>
#include <Fluide_Diphasique.h>
#include <Schema_Comm_FT.h>
#include <Champ_Don_base.h>
#include <Matrice_Morse.h>
#include <Domaine_VDF.h>
#include <TRUSTTabFT.h>
#include <TRUSTLists.h>
#include <TRUST_Ref.h>
#include <type_traits>


class Param;
class Maillage_FT_Disc;
class Transport_Interfaces_FT_Disc;
class Navier_Stokes_FT_Disc;


/*! @brief : class Collision_Model_FT
 *
 *  Description: This class enables to compute solid-solid
 *  interactions for fpi module under the framework of
 *  soft-sphere collision model. Under this framework,
 *  multiple collisions can occurs at the same time (ie a
 *  particle can collide with 2 or more particles). The
 *  collision is spread out on multiple time steps. A slight
 *  overlap (less than the mesh grid size) occurs during the
 *  process.
 */
class Collision_Model_FT_base: public Objet_U
{
  Declare_base_sans_constructeur(Collision_Model_FT_base);

public:

  Collision_Model_FT_base();
  // override functions
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  int reprendre(Entree& is) override;
  int sauvegarder(Sortie& os) const override;
  int preparer_calcul(const Domaine_VDF& domain_vdf,
                      const int nb_particles_tot,
                      const Navier_Stokes_FT_Disc& ns,
                      const Transport_Interfaces_FT_Disc& eq_transport,
                      const Schema_Comm_FT& schema_comm_FT);
  void set_param(Param& p);
  void reset(); // Tables must have the right dimension to be correctly read during the restart
  void resize_geometric_parameters();
  void resize_lagrangian_contact_force()
  {
    lagrangian_contact_forces_.resize(nb_particles_tot_,dimension);
  }
  void resize_particles_collision_number() {particles_collision_number_.resize(nb_particles_tot_);}
  void associate_transport_equation(const Equation_base& equation);
  int check_for_duplicates(ArrOfInt& vector);
  void compute_fictive_wall_coordinates(const double& radius);

  virtual void compute_lagrangian_contact_forces(const Fluide_Diphasique& two_phase_fluid,
                                                 const DoubleTab& particles_position,
                                                 const DoubleTab& particles_velocity,
                                                 const double& deltat_simu)=0;

  virtual void discretize_contact_forces_eulerian_field(const DoubleTab& volumic_phase_indicator_function,
                                                        const Domaine_VF& domain_vf,
                                                        const IntTab& particles_eulerian_id_number,
                                                        DoubleTab& contact_force_source_term)=0;

  void research_collision_pairs_Verlet(const Navier_Stokes_FT_Disc& eq_ns,
                                       const Transport_Interfaces_FT_Disc& eq_transport);

  void compute_Verlet_tables(const DoubleTab& particles_position,
                             const DoubleTab& particles_velocity,
                             double& max_vi,
                             const double& radius,
                             const ArrOfInt& list_particles_to_check_LC);

  DoubleTab compute_contact_force(
    const double& next_dist_int,
    const DoubleTab& norm,
    const DoubleTab& dUn,
    const int& particle_i,
    const int& particle_j,
    const int& is_compression_step,
    const double& is_collision_part_part);

  double compute_ewet_legendre(const double& St) {return exp(-35 / (St + 1e-6));} // See: D. Legendre et al, Chem. Eng. Sci., (2006).

  double compute_stiffness_breugem(const double& mass_eff, const double& e_dry) // See. W-P. Breugem, 2010.
  {return (mass_eff * (pow(M_PI,2)+ pow(log(e_dry), 2))) / pow(collision_duration_, 2);}

  double compute_damper_breugem(const double& mass_eff, const double& e_dry) // See. W-P. Breugem, 2010.
  {return -2*(mass_eff * log(e_dry)) / (collision_duration_);}

  void add_collision(const int part_i, const int part_j, const bool is_part_part_collision);
  bool is_Verlet_activated();
  bool is_LC_activated();

  // setters
  void set_nb_particles_tot(int nb_particles_tot) { nb_particles_tot_=nb_particles_tot; }
  void set_nb_real_particles(int nb_real_particles) { nb_real_particles_=nb_real_particles; }

  void set_activation_distance(const double& diameter)
  {
    activation_distance_=
      diameter*activation_distance_percentage_diameter_/100;
  }
  void set_spring_properties(const Solid_Particle_base& solid_particle);
  void set_domain_dimensions(DoubleVect& Longueurs) { domain_dimensions_=Longueurs; }
  void set_origin(DoubleVect& Origin) { origin_=Origin; }
  void set_geometric_parameters(const Domaine_VDF& domaine_vdf);
  void set_LC_zones(const Domaine_VF& domaine_vf, const Schema_Comm_FT& schema_com);

  // getters
  const int& get_nb_dt_Verlet() const { return nb_dt_Verlet_; }
  const int& get_nb_dt_compute_Verlet() const { return nb_dt_compute_Verlet_; }
  const int& get_nb_dt_max_Verlet() const { return nb_dt_max_Verlet_; }
  const int& get_is_force_on_two_phase_elem() const { return is_force_on_two_phase_elem_; }
  const int& get_collision_number() const { return collision_number_; }
  const int& get_nb_real_particles() const { return nb_real_particles_; }


  int get_last_id(const ArrOfInt& list_particles_to_check_LC);
  int get_id(const ArrOfInt& list_particle, const int ind_id_particle);

  const double& get_duration_collision() const { return collision_duration_; }
  const double& get_delta_n() const { return activation_distance_percentage_diameter_; }
  const double& get_s_Verlet() const  { return detection_thickness_Verlet_; }
  const DoubleVect& get_wall_coordinates() const { return fictive_wall_coordinates_; }
  const DoubleVect& get_origin() const { return origin_; }
  const DoubleVect& get_domain_dimensions() const { return domain_dimensions_; }
  const DoubleTab& get_lagrangian_contact_forces() const { return lagrangian_contact_forces_; }
  const DoubleTab& get_particles_collision_number() const { return particles_collision_number_; }
  const ArrOfInt& get_list_upper_zone() const { return list_upper_zone_; }
  const ArrOfInt& get_list_lower_zone() const { return list_lower_zone_; }
  const IntLists& get_Verlet_table() const { return Verlet_tables_; }
  // getters to fill tab and vectors
  DoubleVect& get_set_collisions_detected() { return collision_detected_; }
  DoubleTab& get_set_e_eff() { return e_eff_; }
  DoubleTab& get_set_F_old() { return F_old_; }
  DoubleTab& get_set_F_now() { return F_now_; }
  DoubleTab& get_set_lagrangian_contact_forces() { return lagrangian_contact_forces_; }
  int& get_set_nb_dt_Verlet() { return nb_dt_Verlet_; }


protected:
  int is_collision_activated_before_impact_ = 1; // to activate, or not, the collision process before the impact
  int is_force_on_two_phase_elem_ = 0; // eulerian discretization of the collision forces on the total particles volume (including two-phase cells)
  int nb_particles_tot_ = 0; // number of particles in the whole domain
  int nb_dt_Verlet_ = 0; // number of time steps since last compute of Verlet tables
  int nb_dt_compute_Verlet_ = 0; // number of time steps to achieve before recomputing Verlet tables
  int nb_dt_max_Verlet_ = 0; // maximum number of time step before recomputing Verlet tables
  int collision_number_ = 0;
  int no_collision_model_=0;
  double activation_distance_percentage_diameter_ = 0; // distance to the wall as a percentage of the particle diameter
  double activation_distance_ = 0;
  double collision_duration_ = 0.; // duration of the collision in seconds
  double detection_thickness_Verlet_ = 0.;

  DoubleTab particles_collision_number_;
  DoubleTab e_eff_; // effective restitution coefficient
  DoubleTab F_old_;
  DoubleTab F_now_;
  DoubleTab lagrangian_contact_forces_;
  DoubleVect collision_detected_;
  ArrOfInt list_upper_zone_;
  ArrOfInt list_lower_zone_;

  OBS_PTR(Transport_Interfaces_FT_Disc) refequation_transport_;

  OBS_PTR(Domaine) ref_domaine;

  void compute_dX_dU(DoubleTab& dX, DoubleTab& dU, const int& particle,\
                     const int& neighbor, const DoubleTab& particles_position, const\
                     DoubleTab& particles_velocity, const bool is_particle_particle_collision );


  int get_nb_particles_j(const int ind_particle_i) const;
  int get_ind_start_particles_j(const int ind_particle_i) const;
  int get_particle_i(const int ind_particle_i);
  int get_particle_j(const int ind_particle_i, const int ind_particle_j);


  IntVect nb_nodes_;
  DoubleVect origin_;
  DoubleVect domain_dimensions_;
  DoubleVect fictive_wall_coordinates_; // an offset is applied to compute collision forces before the impact
  IntLists Verlet_tables_;
  ArrOfInt list_real_particles_;
  int nb_real_particles_; // = nb_particles_tot_ for CHECK_ALL and VERLET

  double stiffness_breugem_part_part_ = 0;
  double stiffness_breugem_wall_part_ = 0;
  double damper_breugem_part_part_ = 0;
  double damper_breugem_wall_part_ = 0;

  enum class Collision_model { HYBRID_ESI, BREUGEM };
  Collision_model collision_model_ = Collision_model::HYBRID_ESI;

  enum class Detection_method { CHECK_ALL, VERLET, LC_VERLET};
  Detection_method detection_method_ = Detection_method::CHECK_ALL;
};

#endif
