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

#ifndef Navier_Stokes_FT_Disc_included
#define Navier_Stokes_FT_Disc_included

#include <Convection_Diffusion_Temperature_FT_Disc.h>
#include <Navier_Stokes_FT_Disc_interne.h>
#include <Navier_Stokes_Turbulent.h>
#include <type_traits>
#include <TRUST_Ref.h>
#include <Collision_Model_FT_base.h>

class Probleme_FT_Disc_gen;
class Maillage_FT_Disc;
class Fluide_Diphasique;

class Navier_Stokes_FT_Disc: public Navier_Stokes_Turbulent
{
  Declare_instanciable(Navier_Stokes_FT_Disc);
public:
  friend class Post_Processing_Hydrodynamic_Forces;

  void set_param(Param& titi) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  const Milieu_base& milieu() const override;
  Milieu_base& milieu() override;
  void associer_pb_base(const Probleme_base& probleme) override;
  void discretiser() override;
  std::vector<YAML_data> data_a_sauvegarder() const override;
  int sauvegarder(Sortie&) const override;
  int reprendre(Entree&) override;
  int preparer_calcul() override;
  void preparer_pas_de_temps();
  void mettre_a_jour(double temps) override;
  void calculer_la_pression_en_pa() override;
  DoubleTab& derivee_en_temps_inco(DoubleTab& vpoint) override;
  void projeter() override;
  virtual const Champ_base& calculer_div_normale_interface();
  void correct_at_exit_bad_gradient(DoubleTab& u0) const;
  void calculer_delta_u_interface(Champ_base& u0, int phase_pilote, int ordre);
  const Champ_Don_base& diffusivite_pour_transport() const override;

  virtual const Champ_base* get_delta_vitesse_interface() const;
  virtual const Fluide_Diphasique& fluide_diphasique() const;

  void compute_boussinesq_additional_gravity(const Convection_Diffusion_Temperature_FT_Disc& eq, const Fluide_Diphasique& fluide_diphasique, const IntTab& face_voisins,
                                             const DoubleVect& volumes_entrelaces, const IntVect& orientation, const DoubleTab& indicatrice, const ArrOfDouble& g, DoubleTab& gravite_face) const;

  int is_terme_gravite_rhog() const;
  const Champ_Fonc_base& champ_rho_faces() const;

  virtual void calculer_dI_dt(DoubleVect& dI_dt); // const;
  const int& get_is_penalized() const;
  const int& get_new_mass_source() const;
  const DoubleTab& get_interfacial_area() const;
  DoubleTab& get_set_interfacial_area();  // Open access  in write-mode..
  const DoubleTab& get_mpoint() const;
  DoubleTab& get_set_mpoint(); // Open access to mpoint in write-mode...
  //void corriger_mpoint(); // Apply correction based on TCL model

  const SolveurSys& get_solveur_pression() const;

  const bool& get_is_solid_particle() const {return is_solid_particle_;}
  const IntTab& get_particles_eulerian_id_number() const { return particles_eulerian_id_number_; }
  // both following methods should be protected but Transport_Interfaces_FT_Disc is not a friend of NS_FT_Disc...
  void compute_particles_eulerian_id_number(const OWN_PTR(Collision_Model_FT_base)& collision_model_ptr);
  void swap_particles_eulerian_id_number(const ArrOfInt& gravity_center_elem);

  DoubleTab& get_set_velocity_field_Stokes_th() { return velocity_field_Stokes_th_->valeurs(); }
  DoubleTab& get_set_pressure_field_Stokes_th() { return pressure_field_Stokes_th_->valeurs(); }

protected:

  // Methode surchargee de Navier_Stokes_std :
  void discretiser_assembleur_pression() override;
  void associer_milieu_base(const Milieu_base& fluide) override;

  // Nouvelles methodes
  virtual const Probleme_FT_Disc_gen& probleme_ft() const;
  virtual Probleme_FT_Disc_gen& probleme_ft();
  virtual void calculer_champ_forces_superficielles(const Maillage_FT_Disc& maillage, const Champ_base& gradient_indicatrice, Champ_base& potentiel_elements, Champ_base& potentiel_faces,
                                                    Champ_base& champ);
  virtual void calculer_gradient_indicatrice(const Champ_base& indicatrice, const DoubleTab& distance_interface_sommets, Champ_base& gradient_i);

  // for fpi module
  void set_is_solid_particle(const bool& is_solid_particle) { is_solid_particle_=is_solid_particle; }
  void set_id_fluid_phase(const int& id_fluid_phase) { id_fluid_phase_=id_fluid_phase; }

  void compute_eulerian_field_contact_forces(const Maillage_FT_Disc& mesh, const Champ_base& field_volumic_phase_indicator_function);

  void eulerian_discretization_contact_forces(const DoubleTab& volumic_phase_indicator_function, const DoubleTab& interlaced_volumes, const DoubleTab& eu);
  OBS_PTR(Probleme_FT_Disc_gen) probleme_ft_;

  // Masse volumique calculee aux elements
  OWN_PTR(Champ_Fonc_base)  champ_rho_elem_;
  // Masse volumique calculee pour les volumes de controle de la vitesse
  // (pour division   v = (rho.v) / rho et pour matrice de pression)
  OWN_PTR(Champ_Fonc_base)  champ_rho_faces_;
  // Viscosite dynamique (calcul dans preparer_pas_de_temps)
  // champ du type requis pour l'operateur diffusion.
  OWN_PTR(Champ_Don_base) champ_mu_;
  // Viscosite cinematique pour le calcul du pas de temps de diffusion
  OWN_PTR(Champ_Don_base) champ_nu_;

  // for fpi_module
  bool is_solid_particle_=false;  // pointer to Fluide_Diphasique::is_solid_particle_
  int id_fluid_phase_=1; // number (0 or 1) of the Fluid_Incompressible phase
  IntTab particles_eulerian_id_number_;
  OWN_PTR(Champ_Fonc_base)  particles_eulerian_id_number_post_; // for post-processing only
  OWN_PTR(Champ_Fonc_base) velocity_field_Stokes_th_;
  OWN_PTR(Champ_Fonc_base) pressure_field_Stokes_th_;

private:
  const Navier_Stokes_FT_Disc_interne& variables_internes() const;
  Navier_Stokes_FT_Disc_interne& variables_internes();

  // Ne pas utiliser ce pointeur : utiliser variables_internes() a la place !
  Navier_Stokes_FT_Disc_interne variables_internes_;

  double minx = -123., maxx = -123., pente = -123.;
  int is_repulsion = 0;

};

#endif /* Navier_Stokes_FT_Disc_included */
