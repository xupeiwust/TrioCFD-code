/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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

#ifndef Probleme_FTD_IJK_base_included
#define Probleme_FTD_IJK_base_included

#include <IJK_Field_vector.h>
#include <IJK_Field.h>
#include <Domaine_IJK.h>
#include <Operateur_IJK_faces_diff.h>
#include <Operateur_IJK_faces_conv.h>
#include <Multigrille_Adrien.h>
#include <Interprete.h>
#include <Linear_algebra_tools.h>
#include <Boundary_Conditions.h>
#include <IJK_Interfaces.h>
#include <Intersection_Interface_ijk.h>
#include <Redistribute_Field.h>
#include <Parser.h>
#include <IJK_FT_Post.h>
#include <IJK_Thermique.h>
#include <Cut_cell_FT_Disc.h>
#include <init_forcage_THI.h>
#include <corrections_qdm.h>
#include <Force_sp.h>
#include <TRUST_List.h>
#include <IJK_Energie.h>
#include <IJK_Thermals.h>
#include <TRUST_Ref.h>
#include <Objet_U.h>
#include <Cut_cell_surface_efficace.h>
#include <Probleme_FT_Disc_gen.h>
#include <Fluide_Diphasique_IJK.h>

class Domaine_IJK;

class Probleme_FTD_IJK_base : public Probleme_FT_Disc_gen
{
  Declare_base_sans_constructeur(Probleme_FTD_IJK_base) ;
  Probleme_FTD_IJK_base();
  Probleme_FTD_IJK_base(const Probleme_FTD_IJK_base& x);

public :
  // We take too much advantage of it ...:
  friend class IJK_Thermique;
  friend class IJK_Thermique_cut_cell;
  friend class Statistiques_dns_ijk_FT;

  /*
   * ===============================================================
   */

  /*
     * ===============================================================
     */
  void preparer_calcul() override { }

  void typer_lire_milieu(Entree& is) override;
  void lire_solved_equations(Entree& is) override;

  void completer() override { }
  void mettre_a_jour(double temps) override { }
  virtual bool updateGivenFields() override { return false; }

  /*
   * Milieu_IJK
   */
  inline Fluide_Diphasique_IJK& milieu_ijk() { return ref_cast(Fluide_Diphasique_IJK, le_milieu_[0].valeur()); }
  inline const Fluide_Diphasique_IJK& milieu_ijk() const { return ref_cast(Fluide_Diphasique_IJK, le_milieu_[0].valeur()); }


  //////////////////////////////////////////////////
  //                                              //
  // Implementation de l'interface de Probleme_U  //
  //                                              //
  //////////////////////////////////////////////////

  // interface Problem
  void initialize() override {  }
  void terminate() override {  }

  // interface UnsteadyProblem
  double presentTime() const override { return 0.; }
  double computeTimeStep(bool& stop) const override { return 0.0; }
  bool initTimeStep(double dt) override { return true; }
  bool solveTimeStep() override { return true; }
  void validateTimeStep() override {  }
  void setStationary(bool flag) override {  }
  void abortTimeStep() override {  }
  void resetTime(double time) override { }

  // interface IterativeUnsteadyProblem
  bool iterateTimeStep(bool& converged) override { return iterateTimeStep_impl(*this, converged); }



  /*
     * ===============================================================
     *//*
* ===============================================================
*/


  /*
   * This
   */
  const Probleme_base& probleme(const Domaine_IJK& dom) const
  {
    if (dom == domaine_ft_)
      return refprobleme_ft_disc_.valeur();
    Cerr << "Unrecognized domain provided" << finl;
    Process::exit();
    return refprobleme_ns_.valeur();
  }
  const Probleme_FTD_IJK_base& operator=(const Probleme_FTD_IJK_base& a) { throw ; }
  int initialise();

  void ecrire_donnees(const IJK_Field_vector3_double& f3compo, SFichier& le_fichier, const int compo, bool binary) const;
  void dumpxyz_vector(const IJK_Field_vector3_double& f3compo, const char * filename, bool binary) const;
  void sauvegarder_probleme(const char *fichier_sauvegarde,
                            const int& stop); //  const;

  void reprendre_probleme(const char *fichier_reprise);



  void discretiser(Discretisation_base& dis) override;


  /*
   *
   * Domaine_ijk
   */


  int associer_(Objet_U&) override;
  const Domaine_IJK& get_domaine_ft() const { return domaine_ft_ ; }
  const Domaine_IJK& get_domaine() const { return domaine_ijk_.valeur(); }

  void redistribute_to_splitting_ft_elem(const IJK_Field_double& input_field,
                                         IJK_Field_double& output_field);

  void redistribute_from_splitting_ft_elem(const IJK_Field_double& input_field,
                                           IJK_Field_double& output_field);

  const Domaine_IJK& domaine_ijk() const { return domaine_ijk_.valeur(); }
  Domaine_IJK& domaine_ijk() { return domaine_ijk_.valeur(); }

  OBS_PTR(Domaine_IJK) domaine_ijk_;
  Domaine_IJK domaine_ft_;

  // Classe outil pour passer entre domain_ijk_ et domain_ft_
  // Une instance par direction des faces:
  FixedVector<Redistribute_Field, 3> redistribute_to_splitting_ft_faces_;
  FixedVector<Redistribute_Field, 3> redistribute_from_splitting_ft_faces_;
  Redistribute_Field redistribute_from_splitting_ft_elem_;
  Redistribute_Field redistribute_from_splitting_ft_elem_ghostz_;
  Redistribute_Field redistribute_from_splitting_ft_elem_ghostz_min_;
  Redistribute_Field redistribute_from_splitting_ft_elem_ghostz_max_;
  /*
   * TIME SCHEME
   *
   *    */
  // Methodes d'acces :
  enum TimeScheme { EULER_EXPLICITE, RK3_FT };
  TimeScheme get_time_scheme() const;
  const double& get_current_time() const
  {
    return current_time_ ;
  }
  const double& get_timestep() const
  {
    return timestep_ ;
  }
  const int& get_nb_timesteps() const
  {
    return nb_timesteps_ ;
  }

  int get_tstep() const
  {
    return tstep_;
  }
  const double& get_dt_cfl() const
  {
    return dt_cfl_;
  }

  const double& get_dt_fo() const
  {
    return dt_fo_;
  }
  const double& get_dt_oh() const
  {
    return dt_oh_;
  }
  const double& get_dt_cfl_liq() const
  {
    return dt_cfl_liq_;
  }
  const double& get_dt_cfl_vap_() const
  {
    return dt_cfl_vap_;
  }
  const double& get_dt_fo_liq() const
  {
    return dt_fo_liq_;
  }
  const double& get_dt_fo_vap_() const
  {
    return dt_fo_vap_;
  }
  virtual void euler_time_step(ArrOfDouble& var_volume_par_bulle) = 0;
  virtual void rk3_sub_step(const int rk_step, const double total_timestep, const double fractionnal_timestep, const double time) = 0;
  // MODIF Aymeric : c'est plus pratique de mettre cette methode publique,
  // c'est une methode generique cont.
  void euler_explicit_update(const IJK_Field_double& dv, IJK_Field_double& v,
                             const int k_layer) const;


  void write_check_etapes_et_termes(const int rk_step);

  double find_timestep(const double max_timestep, const double cfl, const double fo, const double oh);

  double dt_cfl_ = 1.e20;
  double dt_fo_ = 1.e20;
  double dt_oh_ = 1.e20;
  double dt_fo_liq_ = 1.e20;
  double dt_fo_vap_ = 1.e20;
  double dt_cfl_liq_ = 1.e20;
  double dt_cfl_vap_ = 1.e20;


  int enable_dt_oh_ideal_length_factor_ = 0;

  int first_step_interface_smoothing_ = 0;
  int time_scheme_ = EULER_EXPLICITE;
  double store_RK3_source_acc_ = 0.;
  double store_RK3_fac_sv_ = 1.;
  double modified_time_ini_ = 0.;
  double current_time_ = 0.;
  double current_time_at_rk3_step_ = 0;
  int tstep_ = 0; // The iteration number
  int tstep_sauv_ = 0;
  int nb_timesteps_ = 0;
  double max_simu_time_ = 1e6;
  int tstep_init_ = 0;
  int use_tstep_init_ = 0;
  double timestep_ = 0.;
  double timestep_facsec_ = 1.;
  double cfl_ = 1.;
  double fo_ = 1.;
  double oh_ = 1.;



  OBS_PTR(Schema_Temps_base) schema_temps_;

  /*
   * NS
   */
  const IJK_Field_vector3_double& get_velocity() const
  {
    return velocity_;
  }
  const IJK_Field_vector3_double& get_velocity_ft() const
  {
    return velocity_ft_;
  }
  const IJK_Field_double& get_pressure() const
  {
    return pressure_;
  }
  const IJK_Field_double& get_pressure_ghost_cells() const
  {
    return pressure_ghost_cells_;
  }

  const double& get_vitesse_upstream() const
  {
    return vitesse_upstream_;
  }


  const double& get_nb_diam_upstream() const
  {
    return nb_diam_upstream_;
  }
  const int& get_upstream_dir() const
  {
    return upstream_dir_;
  }
  const int& get_upstream_stencil() const
  {
    return upstream_stencil_;
  }


  const int& get_disable_diffusion_qdm() const
  {
    return disable_diffusion_qdm_;
  }
  const int& get_disable_convection_qdm() const
  {
    return disable_convection_qdm_;
  }
  const int& get_compute_rising_velocities() const
  {
    return compute_rising_velocities_;
  }
  const int& get_fill_rising_velocities() const
  {
    return fill_rising_velocities_;
  }


  const IJK_Field_double& get_IJK_field(const Nom& nom) const;

  double calculer_true_moyenne_de_phase_vap(const IJK_Field_double& vx);
  double calculer_moyenne_de_phase_vap(const IJK_Field_double& vx);
  double calculer_true_moyenne_de_phase_liq(const IJK_Field_double& vx);
  double calculer_moyenne_de_phase_liq(const IJK_Field_double& vx);
  void set_time_for_corrections();
  void compute_and_add_qdm_corrections();
  void compute_and_add_qdm_corrections_monophasic();
  void compute_var_volume_par_bulle(ArrOfDouble& var_volume_par_bulle);

  void write_qdm_corrections_information();
  int disable_diphasique() const {return disable_diphasique_;}
  void update_rho_v();
  void update_v_ghost_from_rho_v();
  void update_pressure_phase();
  void terme_source_gravite(IJK_Field_double& dv, int k_index, int dir) const;
  void calculer_dv(const double timestep, const double time, const int rk_step);
  static void force_entry_velocity(IJK_Field_double& vx,
                                   IJK_Field_double& vy,
                                   IJK_Field_double& vz,
                                   double v_imposed,
                                   const int& dir,
                                   const int& compo,
                                   const int& stencil);
  static void force_upstream_velocity(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                      double v_imposed,const IJK_Interfaces& interfaces, double nb_diam,
                                      int upstream_dir, int gravity_dir, int upstream_stencil);
  static void force_upstream_velocity_shear_perio(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                                  double v_imposed,
                                                  const IJK_Interfaces& interfaces,
                                                  double nb_diam, Boundary_Conditions& bc, double nb_diam_ortho_shear_perio,double Ux0,double Uy0,double Uz0,
                                                  int epaisseur_maille);


  void calculer_terme_source_acceleration(IJK_Field_double& vx, const double time, const double timestep,
                                          const int rk_step);
  // Correcteur PID
  void calculer_terme_asservissement(double& ax, double& ay, double& az);
  void calculer_vitesse_gauche(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz, double& vx_moy, double& vy_moy, double& vz_moy);
  void calculer_vitesse_droite(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz, double& vx_moy, double& vy_moy, double& vz_moy);

  void compute_correction_for_momentum_balance(const int rk_step);
  void compute_add_external_forces(const int dir);
  void fill_variable_source_and_potential_phi(const double time);




  Vecteur3 calculer_inv_rho_grad_p_moyen(const IJK_Field_double& inv_rho, const IJK_Field_double& pression);
  Vecteur3 calculer_grad_p_moyen(const IJK_Field_double& pression);
  Vecteur3 calculer_grad_p_over_rho_moyen(const IJK_Field_double& pression);
  IJK_Field_vector3_double terme_convection_mass_solver_;
  IJK_Field_vector3_double terme_diffusion_mass_solver_;
  IJK_Field_vector3_double rho_u_euler_av_prediction_champ_;
  IJK_Field_vector3_double rho_du_euler_ap_prediction_champ_;
  IJK_Field_vector3_double rho_u_euler_ap_projection_champ_;
  IJK_Field_vector3_double rho_du_euler_ap_projection_champ_;
  IJK_Field_vector3_double rho_u_euler_av_rho_mu_ind_champ_;
  IJK_Field_vector3_double rho_u_euler_ap_rho_mu_ind_champ_;
  IJK_Field_vector3_double terme_diffusion_local_;
  IJK_Field_vector3_double terme_pression_local_;
  IJK_Field_vector3_double terme_pression_in_ustar_local_;
  IJK_Field_vector3_double d_v_diff_et_conv_;
  Vecteur3 rho_u_euler_av_prediction_ = {0.,0.,0.};
  Vecteur3 rho_du_euler_ap_prediction_ = {0.,0.,0.};
  Vecteur3 rho_u_euler_ap_projection_ = {0.,0.,0.};
  Vecteur3 rho_du_euler_ap_projection_ = {0.,0.,0.};
  Vecteur3 rho_u_euler_av_rho_mu_ind_ = {0.,0.,0.};
  Vecteur3 rho_u_euler_ap_rho_mu_ind_ = {0.,0.,0.};
  Vecteur3 u_euler_ap_rho_mu_ind_ = {0.,0.,0.};
  Vecteur3 terme_diffusion_ = {0.,0.,0.};
  Vecteur3 terme_convection_ = {0.,0.,0.};
  Vecteur3 terme_pression_ = {0.,0.,0.};
  Vecteur3 terme_pression_bis_ = {0.,0.,0.};
  Vecteur3 terme_pression_ter_ = {0.,0.,0.};
  Vecteur3 terme_interfaces_;
  Vecteur3 terme_pression_in_ustar_ = {0.,0.,0.};
  Vecteur3 terme_moyen_convection_mass_solver_ = {0.,0.,0.};
  Vecteur3 terme_moyen_diffusion_mass_solver_ = {0.,0.,0.};
  double pression_ap_proj_ = 0.;
  double vap_velocity_tmoy_ = 0.;
  double reprise_vap_velocity_tmoy_ = 0.;
  double liq_velocity_tmoy_ = 0.;
  double reprise_liq_velocity_tmoy_ = 0.;

  int compute_rising_velocities_ = 0;
  int fill_rising_velocities_ = 0;
  corrections_qdm qdm_corrections_;



  int disable_solveur_poisson_ = 0;
  int resolution_fluctuations_ = 0;
  int projection_initiale_demandee_ = 0;
  int disable_diffusion_qdm_ = 0;
  int disable_convection_qdm_ = 0;
  int disable_source_interf_ = 0;
  int frozen_velocity_ = 0;
  int velocity_reset_ = 0;

  // Correcteur PID
  double Kp_ = 0.;
  double Kd_ = 0.;
  double Ki_ = 0.;
  double int_x_ = 0.;
  double int_y_ = 0.;
  double int_z_ = 0.;
  int epaisseur_maille_ = 8;

  int vitesse_entree_dir_ = DIRECTION_I;
  int vitesse_entree_compo_to_force_ = -1;
  double vitesse_entree_ = -1.1e20;
  int stencil_vitesse_entree_ = 3;
  double vitesse_upstream_ = -1.1e20;
  double vitesse_upstream_reprise_ = -1.1e20;
  double velocity_bubble_new_ = 0.;
  double velocity_bubble_old_ = -1.1e20;
  double velocity_bubble_integral_err_ = 0.;
  double velocity_bubble_scope_ = 0.;
  double upstream_velocity_bubble_factor_ = 1.;
  double upstream_velocity_bubble_factor_deriv_ = 0.;
  double upstream_velocity_bubble_factor_integral_ = 0.;
  Nom expression_vitesse_upstream_ = "??";
  // Molecular diffusivity (see diffusion operator)
  IJK_Field_double molecular_mu_;
  // right hand side for pressure solver
  IJK_Field_double pressure_rhs_;
  // Operators and pressure solver

  /* Velocity diffusion operator:
   * simple_arithmetic : div(mu grad(u))
   * full_arithmetic    : Tenseur des contraintes complet : div[mu (grad(u)+grad^T(u))]
   *                      mu : moyenne arithmetique
   * full_adaptative    : Tenseur des contraintes complet : div[mu (grad(u)+grad^T(u))]
   *     mu : switch from arithmetic to geometric mean depending on the direction (Not available yet)
   */
  int use_harmonic_viscosity_ = 0;
  Operateur_IJK_faces_diff velocity_diffusion_op_;
  enum velocity_diffusion_options_ { simple_arithmetic, full_arithmetic, full_adaptative};

  /*
   * Velocity convection operator
   * non_conservative_simple : rho div(u u)
   * non_conservative_rhou   : div(rho u u) - u div(rho u)
   * conservative            : div(rho u u)
   */
  Operateur_IJK_faces_conv velocity_convection_op_;
  enum velocity_convection_options_ { non_conservative_simple, non_conservative_rhou, conservative};

  Multigrille_Adrien poisson_solver_;
  // Simulation parameters

  // Pressure field
  IJK_Field_double pressure_;
  IJK_Field_double pressure_ghost_cells_;

  int upstream_dir_ = -1; // static
  int upstream_stencil_ = 3;
  int upstream_velocity_measured_ = 0;
  double nb_diam_upstream_ = 0.;
  double nb_diam_ortho_shear_perio_= -1.1e20;

  int harmonic_nu_in_diff_operator_ = 0;
  int harmonic_nu_in_calc_with_indicatrice_ = 0;




public:


  /*
   * Thermique
   */

  void compute_add_THI_force(const IJK_Field_vector3_double& vitesse,
                             const int time_iteration,
                             const double dt,
                             const double current_time,
                             const Domaine_IJK& my_splitting
                            );
  void compute_add_THI_force_sur_d_velocity(const IJK_Field_vector3_double& vitesse,
                                            const int time_iteration,
                                            const double dt,
                                            const double current_time,
                                            const Domaine_IJK& my_splitting,
                                            const int facteur
                                           );
  // Dealing with thermal aspects:
  LIST(IJK_Thermique) thermique_;
  LIST(IJK_Energie) energie_;
  IJK_Thermals thermals_;
  int thermal_probes_ghost_cells_ = 2;
  init_forcage_THI forcage_;

  /*
   * Milieu
   */
  const IJK_Field_double& get_rho_field() const
  {
    return rho_field_;
  }


  double get_rho_field_ijk(int i, int j, int k) const
  {
    return rho_field_(i,j,k);
  }
  int get_disable_diphasique() const
  {
    return disable_diphasique_;
  }


  double rho_moyen_ = 0.;
  //








  /*
   * Post
   */


  double t_debut_statistiques() const
  {
    return post_.t_debut_statistiques();
  }
  int get_reprise() const
  {
    return reprise_;
  }
  const IJK_FT_Post& get_post() const {return post_;}



  /*
   * FT
   */



  const Maillage_FT_IJK& get_maillage_ft_ijk() const
  {
    return interfaces_.maillage_ft_ijk();
  }
  const Remaillage_FT_IJK& get_remaillage_ft_ijk() const
  {
    return interfaces_.remaillage_ft_ijk();
  }

  const IJK_Interfaces& get_interface() const
  {
    return interfaces_;
  }

  virtual void update_indicator_field();
  void update_pre_remeshing_indicator_field();
  void update_post_remeshing_indicator_field();
  virtual void update_twice_indicator_field();

  void update_old_intersections();

  virtual Cut_cell_FT_Disc* get_cut_cell_disc()
  {
    Cerr << "No cut fields are found." << finl;
    Process::exit();
    return nullptr;
  }
  const IJK_Interfaces& itfce() const {return interfaces_;}
  IJK_Interfaces& get_set_interface() {return interfaces_;}


  void parcourir_maillage();
  void calculer_rho_mu_indicatrice(const bool parcourir = true);
  void maj_indicatrice_rho_mu(const bool parcourir = true);
  Redistribute_Field& redistrib_to_ft_elem() {return redistribute_to_splitting_ft_elem_;}
  Redistribute_Field& redistrib_from_ft_elem() {return redistribute_from_splitting_ft_elem_;}
  // FixedVector<Redistribute_Field, 3>& redistrib_from_ft_face() const {return redistribute_from_splitting_ft_faces_;}
  void get_redistribute_from_splitting_ft_faces(
    const IJK_Field_vector3_double& faces_ft,
    IJK_Field_vector3_double& faces_ns
  )
  {
    for (int dir = 0; dir < 3; dir++)
      {
        redistribute_from_splitting_ft_faces_[dir].redistribute(
          faces_ft[dir],
          faces_ns[dir]);
      }
  }


  virtual void deplacer_interfaces(const double timestep,
                                   const int rk_step,
                                   ArrOfDouble& var_volume_par_bulle,
                                   const int first_step_interface_smoothing);
  virtual void deplacer_interfaces_rk3(const double timestep, const int rk_step, ArrOfDouble& var_volume_par_bulle);
  void calculer_gradient_indicatrice_et_repul_ns(const IJK_Field_double& indic);


  void transfer_ft_to_ns();
  // att
  Vecteur3 terme_interfaces_bf_mass_solver_ = {0.,0.,0.};
  Vecteur3 terme_interfaces_bf_mass_solver_bis_ = {0.,0.,0.};
  Vecteur3 terme_interfaces_af_mass_solver_ = {0.,0.,0.};
  Vecteur3 terme_interfaces_conv_diff_mass_solver_ = {0.,0.,0.};
  int use_bubbles_velocities_from_interface_ = 0;
  int use_bubbles_velocities_from_barycentres_ = 0;
  TYPE_SURFACE_EFFICACE_FACE type_surface_efficace_face_ = TYPE_SURFACE_EFFICACE_FACE::NON_INITIALISE;
  TYPE_SURFACE_EFFICACE_INTERFACE type_surface_efficace_interface_ = TYPE_SURFACE_EFFICACE_INTERFACE::NON_INITIALISE;
  int deactivate_remeshing_velocity_ = 0;

  DoubleTab vitesses_translation_bulles_; // Vecteur de translation rigide pour chaque bulle
  DoubleTab mean_bubble_rotation_vector_; // Vecteur de rotation rigide pour chaque bulle
  DoubleTab centre_gravite_bulles_;       // Position du centre de gravite pour chaque bulle (associee a la rotation)
  int correction_semi_locale_volume_bulle_ = 0;


  IJK_Interfaces interfaces_;
  IJK_Field_vector3_double velocity_ft_;

  Redistribute_Field redistribute_to_splitting_ft_elem_;





  /*
   * Tools
   */
  void copy_field_values(IJK_Field_double& field, const IJK_Field_double& field_to_copy);

  IJK_Field_double scalar_product(const IJK_Field_vector3_double& V1, const IJK_Field_vector3_double& V2);
  IJK_Field_vector3_double scalar_times_vector(const IJK_Field_double& Sca, const IJK_Field_vector3_double& Vec);
  IJK_Field_double scalar_fields_product(const IJK_Field_double& S1, const IJK_Field_double& S2, int dir);





  /*
   * ===============================================================
   */



  Nom expression_derivee_acceleration_ = "0"; // par defaut pas de terme d'acceleration
  Parser parser_derivee_acceleration_;
  Noms expression_variable_source_; // on attend trois expressions
  Nom expression_potential_phi_ = "??"; // source variable formulee en gradient
  Vecteur3 store_rhov_moy_;
  Vecteur3 integrated_residu_ = {0.,0.,0.};
  // terme source qdm pour pousser le fluide dans le canal (en m/s/s)
  double terme_source_acceleration_ = 0.; // par defaut, zero
  int compute_force_init_ = 0;

  // Vecteurs de taille 3 a lire dans le jeu de donnees :
  ArrOfDouble terme_source_correction_; // Valeur de la force de correction moyenne a appliquer
  ArrOfInt correction_force_; // 3 flags d'activation de la correction ou non

  // Storage for rhov for the evaluation of the acceleration source term :
  IJK_Field_vector3_double rho_v_;
  IJK_Field_vector3_double psi_velocity_; // for storage.
  //ab-forcage-control-ecoulement-fin

  IJK_Field_vector3_double variable_source_;
  Nom expression_derivee_facteur_variable_source_ = "0";
  Parser parser_derivee_facteur_variable_source_;
  double facteur_variable_source_ = 1.; // ArrOfDouble? vecteur de taille 3 a lire dans le jeu de donnees

  IJK_Field_double potential_phi_;

  /*
   * Post-processing helper class:
   */
  friend class IJK_FT_Post;
  IJK_FT_Post post_;

  int check_divergence_ = 0;
  int rk_step_ = -1; // default value

  Nom check_stop_file_; // Nom du fichier stop

  //ab-sauv/repr-deb

  Nom fichier_post_ = "??"; // Nom du fichier post
  int dt_sauvegarde_ = 2000000000;
  int sauvegarder_xyz_ = 0; // drapeau 0 ou 1
  Nom nom_sauvegarde_;
  Nom nom_reprise_;
  int reprise_ = 0;// flag pour indiquer si on fait une reprise
  // Le jeu de donnees doit fournir soit des fichiers de reprise: ..
  Nom fichier_reprise_vitesse_ = "??"; // par defaut, invalide
  int timestep_reprise_vitesse_ = 1;
  // ... soit des expressions f(x,y,z)
  Noms expression_vitesse_initiale_; // on attend trois expressions
  Nom expression_pression_initiale_ = "??"; // useless, unless post-pro OR pressure_increment.
  //ab-sauv/repr-fin

  // Le probleme ft disc qui porte le maillage vdf pour les algorithmes front-tracking
  OBS_PTR(Probleme_base) refprobleme_ft_disc_;
  // Creation d'un probleme sur le domaine d'origine pour les sondes et pour faire leur VDF...
  OBS_PTR(Probleme_base) refprobleme_ns_;

  ArrOfDouble_with_ghost delta_z_local_;

  Boundary_Conditions boundary_conditions_;

  // Inconnues du probleme (a sauvegarder et a reprendre)
  // Velocity field:
  IJK_Field_vector3_double velocity_;

  // Masse volumique:
  IJK_Field_double rho_field_;
  IJK_Field_double inv_rho_field_;

  // Pour les cas a bulles fixes
  // valeurs par default des parametres de bulles fixes
  double coef_immobilisation_ = 0.;
  double coef_ammortissement_ = 0.;
  double coef_mean_force_ = 0.;
  double coef_force_time_n_ = 0.;
  double coef_rayon_force_rappel_ = 0.;
  double p_seuil_max_ = 10000000;
  double p_seuil_min_ = -10000000;
  IJK_Field_vector3_double force_rappel_;
  IJK_Field_vector3_double force_rappel_ft_;

  double vol_bulle_monodisperse_ = -1; // Pour imposer le volume des bulles
  double diam_bulle_monodisperse_ = -1; // Pour imposer le volume des bulles
  ArrOfDouble vol_bulles_;   // Le volume impose individuellement a chaque bulle.
  double coeff_evol_volume_ = 0.;

  // Field only needed for the option type_velocity_convection_form_== Nom("non_conservative_rhou")
  IJK_Field_double div_rhou_;

  // Temporary storage for the derivative
  IJK_Field_vector3_double d_velocity_;
  // Temporary storage for the fractional timestep in rk3 :
  IJK_Field_vector3_double RK3_F_velocity_;

  // To work in pressure increment and potentially help the solver for high density ratios
  IJK_Field_double d_pressure_;
  // Only if pressure increment formulation AND RK3 time_scheme.
  IJK_Field_double RK3_F_pressure_;

  // Celui la est discretise sur le maillage etendu:
  IJK_Field_vector3_double terme_source_interfaces_ft_;
  IJK_Field_vector3_double terme_repulsion_interfaces_ft_;
  IJK_Field_vector3_double terme_abs_repulsion_interfaces_ft_;

  // Celui la est discretise sur le maillage ns:
  IJK_Field_vector3_double terme_source_interfaces_ns_; // [ dt (rho u )/dt = N/m^3 ]
  IJK_Field_vector3_double backup_terme_source_interfaces_ns_; // [ dt (rho u )/dt = N/m^3 ]
  IJK_Field_vector3_double backup_terme_source_interfaces_ft_; // [ dt (rho u )/dt = N/m^3 ]
  IJK_Field_vector3_double terme_repulsion_interfaces_ns_;
  IJK_Field_vector3_double terme_abs_repulsion_interfaces_ns_;

  // Pour l'option alternative du calcul de la diffusion :
  IJK_Field_double unit_;
  IJK_Field_vector3_double zero_field_ft_;
  IJK_Field_vector3_double laplacien_velocity_;


  // Booleen pour savoir si on a mis a jour l'indicatrice avec indicatrice next
  // bool indicatrice_is_indicatrice_next_;



  IJK_Field_double kappa_ft_;
  IJK_Field_double kappa_ns_;

  IJK_Field_double I_ns_;







  // Pour la premiere projection, on initialise la pression au champ de pression a l'equilibre diphasique :
  int improved_initial_pressure_guess_ = 0;
  // travail en increment de pression pour aider le solveur :
  int include_pressure_gradient_in_ustar_ = 0;
  // Discretisation du champ (1/rho) et utilisation dans le calcul de rho*v*v, dans le mass_solver
  // et dans pressure_projection_with_inv_rho :
  int use_inv_rho_for_mass_solver_and_calculer_rho_v_ = 0;
  int use_inv_rho_in_poisson_solver_ = 0;
  int use_inv_rho_ = 0;

  int correction_bilan_qdm_ = 0;

  int refuse_patch_conservation_QdM_RK3_source_interf_ = 0; // Par defaut, on utilise le patch!
  // GAB, qdm
  int test_etapes_et_bilan_ = 0;
  //
  // GAB, champ de reprise + champ initial
  int add_initial_field_ = 0;
  //
  int diffusion_alternative_ = 0;
  int suppression_rejetons_ = 0;  // By defaults, break-ups are not fixed on restart. (no deletion of smaller fractions)
  // Supprime l'appel a quelques fonctions qui n'ont pas de sens en monophasique :
  // comme par exemple : deplacer_interfaces,
  // calculer_rho_mu_indicatrice, ecrire_statistiques_bulles
  int disable_diphasique_ = 0;



};

#endif /* Probleme_FTD_IJK_base_included */
