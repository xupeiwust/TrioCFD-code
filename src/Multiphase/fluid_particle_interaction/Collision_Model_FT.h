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

#ifndef Modele_Collision_FT_included
#define Modele_Collision_FT_included

#include <TRUSTTabFT_forward.h>
#include <TRUSTTabFT.h>
#include <Fluide_Diphasique.h>
#include <Domaine_VDF.h>

class Param;
class Maillage_FT_Disc;
class Transport_Interfaces_FT_Disc;

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

class Collision_Model_FT : public Objet_U
{
  Declare_instanciable_sans_constructeur(Collision_Model_FT);

public:

  Collision_Model_FT();

  // override functions
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  int reprendre(Entree& is) override;
  int sauvegarder(Sortie& os) const override;

  void reset(); // Tables must have the right dimension to be correctly read during the restart
  void resize_geometric_parameters();
  void associate_transport_equation(const Equation_base& equation);
  int check_for_duplicates(ArrOfInt& vector);
  void compute_contact_force(DoubleTab& force_contact, int& isFirstStepOfCollision, double& dist_int, double& next_dist_int, DoubleTab& norm, DoubleTab& dUn, double& masse_eff, int& compo, int& voisin, double& Stb, double& ed, double& vitesseRelNorm, double& dt, double& prod_scal);
  void compute_fictive_wall_coordinates(const double& radius);
  double compute_ewet_legendre(const double& St) {return exp(-35 / (St + 1e-6));} // See: D. Legendre et al, Chem. Eng. Sci., (2006).
  double compute_stiffness_breugem(const double& mass_eff, const double& e_dry) {return (mass_eff * (pow(M_PI,2)+ pow(log(e_dry), 2))) / pow(collision_duration_, 2);} // See. W-P. Breugem, 2010.
  double compute_damper_breugem(const double& mass_eff, const double& e_dry) {return -2*(mass_eff * log(e_dry)) / (collision_duration_);}  // See. W-P. Breugem, 2010.

  // setters
  void set_param(Param& p);
  void set_nb_compo_tot(int nb_compo_tot) {nb_compo_tot_=nb_compo_tot;}
  void detection_thickness_Verlet(double s_Verlet) {detection_thickness_Verlet_=s_Verlet;}
  void set_domain_dimensions(DoubleVect& Longueurs) {domain_dimensions_=Longueurs;}
  void set_origin(DoubleVect& Origin) {origin_=Origin;}
  void set_geometric_parameters(Domaine_VDF& domaine_vdf);

  // getters
  const int& get_nb_dt_Verlet() const { return nb_dt_Verlet_; }
  const int& get_dt_compute_Verlet() const { return dt_compute_Verlet_; }
  const int& get_nb_pas_dt_max_Verlet() const { return nb_pas_dt_max_Verlet_; }
  const int& get_is_force_on_two_phase_elem() const { return is_force_on_two_phase_elem_; }
  const int& get_is_detection_Verlet() const { return is_detection_Verlet_;}
  const int& get_is_LC_activated() const { return is_linked_cell_activated_; }
  const double& get_duration_collision() const { return collision_duration_; }
  const double& get_delta_n() const { return activation_distance_percentage_diameter_; }
  const double& get_s_Verlet() const  { return detection_thickness_Verlet_; }
  const DoubleVect& get_wall_coordinates() const { return fictive_wall_coordinates_; }
  const ArrOfIntFT& get_list_upper_zone() const { return list_upper_zone_; }
  const ArrOfIntFT& get_list_lower_zone() const { return list_lower_zone_; }

  // getters to fill tab and vectors
  DoubleVect& get_collisions_detected() { return collision_detected_; }
  DoubleTab& get_stiffness() { return stiffness_; }
  DoubleTab& get_e_eff() { return e_eff_; }
  DoubleTab& get_F_old() { return F_old_; }
  DoubleTab& get_F_now() { return F_now_; }
  DoubleTab& get_solid_forces() { return solid_forces_; }

protected:
  int is_collision_activated_before_impact_=1; // to activate, or not, the collision process before the impact
  int is_force_on_two_phase_elem_=0; // eulerian discretization of the collision forces on the total particles volume (including two-phase cells)
  int nb_compo_tot_=0; // number of particles in the whole domain
  int nb_dt_Verlet_=0;
  int dt_compute_Verlet_=0;
  int nb_pas_dt_max_Verlet_=0; // maximum number of time step before recomputing Verlet tables
  int is_detection_Verlet_=0; // to activate, or not, the detection with Verlet method
  int is_linked_cell_activated_=0;

  double activation_distance_percentage_diameter_=0; // distance to the wall as a percentage of the particle diameter
  double collision_duration_=0.; // duration of the collision in seconds
  double detection_thickness_Verlet_=0.;

  DoubleTab stiffness_; // stiffness of the spring
  DoubleTab e_eff_; // effective restitution coefficient
  DoubleTab F_old_;
  DoubleTab F_now_;
  DoubleTab solid_forces_;
  DoubleVect fictive_wall_coordinates_; // an offset is applied to compute collision forces before the impact
  DoubleVect collision_detected_;
  ArrOfIntFT list_upper_zone_;
  ArrOfIntFT list_lower_zone_;

  OBS_PTR(Transport_Interfaces_FT_Disc) refequation_transport_;

private:
  OBS_PTR(Domaine) ref_domaine;

  Nom filename_data_fpi_="data_fpi.sauv";

  IntVect nb_nodes_;
  DoubleVect origin_;
  DoubleVect domain_dimensions_;

  enum Collision_model { HYBRID_ESI, BREUGEM };
  Collision_model collision_model_=HYBRID_ESI;
};

#endif

