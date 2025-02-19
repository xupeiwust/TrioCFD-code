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

#ifndef Navier_Stokes_FTD_IJK_included
#define Navier_Stokes_FTD_IJK_included

#include <Cut_cell_surface_efficace.h>
#include <Operateur_IJK_faces_diff.h>
#include <Operateur_IJK_faces_conv.h>
#include <Fluide_Diphasique_IJK.h>
#include <Schema_Temps_IJK_base.h>
#include <Boundary_Conditions.h>
#include <Multigrille_Adrien.h>
#include <init_forcage_THI.h>
#include <IJK_Field_vector.h>
#include <corrections_qdm.h>
#include <IJK_Interfaces.h>
#include <Equation_base.h>
#include <IJK_Field.h>
#include <TRUST_Ref.h>


class Probleme_FTD_IJK_base;
class Fluide_base;

class Navier_Stokes_FTD_IJK: public Equation_base
{
  Declare_instanciable_sans_constructeur(Navier_Stokes_FTD_IJK);
public:

  friend class Postprocessing_IJK;
  friend class Statistiques_dns_ijk_FT;

  Navier_Stokes_FTD_IJK();

  void set_param(Param& titi) override;
  void set_param_reprise_pb(Param& );

  void completer() override;
  void associer_pb_base(const Probleme_base&) override;
  void discretiser() override { }
  int preparer_calcul() override;
  void projeter();
  const Milieu_base& milieu() const override;
  Milieu_base& milieu() override;
  void associer_milieu_base(const Milieu_base& ) override;

  int nombre_d_operateurs() const override { return -123; }
  const Operateur& operateur(int) const override { throw; }
  Operateur& operateur(int) override { throw; }
  const Champ_Inc_base& inconnue() const override { throw; }
  Champ_Inc_base& inconnue() override { throw; }

  Probleme_FTD_IJK_base& probleme_ijk();
  const Probleme_FTD_IJK_base& probleme_ijk() const;

  const IJK_Field_double& get_IJK_field(const Nom& nom) const;
  bool has_IJK_field(const Nom& nom) const;
  void initialise_ijk_fields();
  void initialise_ns_fields();
  void complete_initialise_ijk_fields();

  const Boundary_Conditions& get_boundary_conditions() const { return boundary_conditions_; }
  const IJK_Field_double& get_pressure() const { return pressure_; }
  const IJK_Field_double& get_pressure_ghost_cells() const { return pressure_ghost_cells_; }

  double get_vitesse_upstream() const { return vitesse_upstream_; }
  double get_nb_diam_upstream() const { return nb_diam_upstream_; }

  int get_upstream_dir() const { return upstream_dir_; }
  int get_upstream_stencil() const { return upstream_stencil_; }

  const IJK_Field_double& get_rho_field() const { return rho_field_; }
  double get_rho_field_ijk(int i, int j, int k) const { return rho_field_(i,j,k); }

  const IJK_Field_vector3_double& get_velocity() const { return velocity_; }
  IJK_Field_vector3_double& get_velocity()  { return velocity_; }

  int get_disable_diffusion_qdm() const { return disable_diffusion_qdm_; }
  int get_disable_convection_qdm() const { return disable_convection_qdm_; }

  int get_compute_rising_velocities() const { return compute_rising_velocities_; }
  int get_fill_rising_velocities() const { return fill_rising_velocities_; }
  int get_use_bubbles_velocities_from_interface() const { return use_bubbles_velocities_from_interface_; }
  int get_use_bubbles_velocities_from_barycentres() const { return use_bubbles_velocities_from_barycentres_; }
  int get_upstream_velocity_measured() const { return upstream_velocity_measured_; }
  int& get_compute_rising_velocities() { return compute_rising_velocities_; }
  int& get_fill_rising_velocities() { return fill_rising_velocities_; }
  int& get_use_bubbles_velocities_from_interface() { return use_bubbles_velocities_from_interface_; }
  int& get_use_bubbles_velocities_from_barycentres() { return use_bubbles_velocities_from_barycentres_; }
  int& get_upstream_velocity_measured() { return upstream_velocity_measured_; }

  const IJK_Field_vector3_double& get_velocity_ft() const { return velocity_ft_; }

  void associer_interfaces(const IJK_Interfaces& inter) { interfaces_ = inter; }
//   void associer_domaine_ft(const Domaine_IJK& dom) { domaine_ft_ = dom; }

  inline Fluide_Diphasique_IJK& milieu_ijk() { return ref_cast(Fluide_Diphasique_IJK, milieu()); }
  inline const Fluide_Diphasique_IJK& milieu_ijk() const { return ref_cast(Fluide_Diphasique_IJK, milieu()); }
  inline Schema_Temps_IJK_base& schema_temps_ijk() { return ref_cast(Schema_Temps_IJK_base, le_schema_en_temps.valeur()); }
  inline const Schema_Temps_IJK_base& schema_temps_ijk() const { return ref_cast(Schema_Temps_IJK_base, le_schema_en_temps.valeur()); }

  void redistribute_to_splitting_ft_elem(const IJK_Field_double& input_field, IJK_Field_double& output_field);
  void redistribute_from_splitting_ft_elem(const IJK_Field_double& input_field, IJK_Field_double& output_field);

  Redistribute_Field& redistrib_to_ft_elem() { return redistribute_to_splitting_ft_elem_; }
  Redistribute_Field& redistrib_from_ft_elem() { return redistribute_from_splitting_ft_elem_; }
  // FixedVector<Redistribute_Field, 3>& redistrib_from_ft_face() const {return redistribute_from_splitting_ft_faces_;}
  void get_redistribute_from_splitting_ft_faces(const IJK_Field_vector3_double& faces_ft, IJK_Field_vector3_double& faces_ns)
  {
    for (int dir = 0; dir < 3; dir++)
      redistribute_from_splitting_ft_faces_[dir].redistribute(faces_ft[dir], faces_ns[dir]);
  }

  void maj_indicatrice_rho_mu(const bool parcourir = true);

  void update_v_ghost_from_rho_v();
  void update_rho_v();

  void transfer_ft_to_ns();

  // Correcteur PID
  void calculer_terme_asservissement(double& ax, double& ay, double& az);
  void calculer_vitesse_gauche(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz, double& vx_moy, double& vy_moy, double& vz_moy);
  void calculer_vitesse_droite(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz, double& vx_moy, double& vy_moy, double& vz_moy);

  void calculer_terme_source_acceleration(IJK_Field_double& vx, const double time, const double timestep, const int rk_step);

  void calculer_terme_source_acceleration(const double time, const double timestep, const int rk_step, const int );

  void compute_correction_for_momentum_balance(const int rk_step);

  void calculer_dv(const double timestep, const double time, const int rk_step);

  void compute_add_THI_force(const IJK_Field_vector3_double& vitesse, const int time_iteration, const double dt, const double current_time, const Domaine_IJK& my_splitting);

  void compute_add_THI_force_sur_d_velocity(const IJK_Field_vector3_double& vitesse, const int time_iteration, const double dt, const double current_time, const Domaine_IJK& my_splitting,
                                            const int facteur);

  double calculer_moyenne_de_phase_liq(const IJK_Field_double& vx);

  void compute_and_add_qdm_corrections();

  void fill_variable_source_and_potential_phi(const double time);

  void write_check_etapes_et_termes(const int rk_step);

  void compute_add_external_forces(const int dir);

  double calculer_true_moyenne_de_phase_vap(const IJK_Field_double& vx);
  double calculer_moyenne_de_phase_vap(const IJK_Field_double& vx);
  double calculer_true_moyenne_de_phase_liq(const IJK_Field_double& vx);

  void terme_source_gravite(IJK_Field_double& dv, int k_index, int dir) const;

  void set_time_for_corrections();
  void compute_and_add_qdm_corrections_monophasic();
  void compute_var_volume_par_bulle(ArrOfDouble& var_volume_par_bulle);
  void write_qdm_corrections_information();

  Vecteur3 calculer_inv_rho_grad_p_moyen(const IJK_Field_double& inv_rho, const IJK_Field_double& pression);
  Vecteur3 calculer_grad_p_moyen(const IJK_Field_double& pression);
  Vecteur3 calculer_grad_p_over_rho_moyen(const IJK_Field_double& pression);

  void euler_explicit_update(const IJK_Field_double& dv, IJK_Field_double& v, const int k_layer) const;

  void forcage_control_ecoulement();
  const IJK_Field_double& get_molecular_mu() const { return molecular_mu_; }
  int get_improved_initial_pressure_guess() const { return improved_initial_pressure_guess_; }
  int get_suppression_rejetons() const { return suppression_rejetons_; }

  void rk3_sub_step(const int rk_step, const double total_timestep, const double fractionnal_timestep, const double time);
  void euler_time_step(ArrOfDouble& var_volume_par_bulle);
  void update_v_or_rhov(bool with_p = false);

  void corriger_qdm();

  // TODO FIXME DANS DOMAINE
  void build_redistribute_extended_splitting_ft();

  void test_etapes_et_bilan_rho_u_euler(bool apres);

  void deplacer_interfaces(const double timestep, const int rk_step, ArrOfDouble& var_volume_par_bulle, const int first_step_interface_smoothing);
  void deplacer_interfaces_rk3(const double timestep, const int rk_step, ArrOfDouble& var_volume_par_bulle);
  void calculer_vitesse_ft();
  void update_indicatrice_variables_monofluides();

  void sauvegarder_equation(const Nom&, SFichier& ) const;

  void set_fichier_reprise_vitesse(const Nom& prefix) { fichier_reprise_vitesse_ = prefix + fichier_reprise_vitesse_; }

  void create_forced_dilation();

  bool get_flag_variable_source() const { return flag_variable_source_; }

protected:
  OBS_PTR(IJK_Interfaces) interfaces_;

  void initialise_velocity_using_expression(const Noms& expression_vitesse_initiale);
  void initialise_velocity_from_file(const Nom& fichier_reprise_vitesse);

  // TODO FIXME ELie : move to milieu
  // Masse volumique:
  double rho_moyen_ = 0.;
  IJK_Field_double rho_field_;
  IJK_Field_double inv_rho_field_;
  // Molecular diffusivity (see diffusion operator)
  IJK_Field_double molecular_mu_;


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

  // TODO FIXME pas ici ...
  // Classe outil pour passer entre domain_ijk_ et domain_ft_
  // Une instance par direction des faces:
  FixedVector<Redistribute_Field, 3> redistribute_to_splitting_ft_faces_;
  FixedVector<Redistribute_Field, 3> redistribute_from_splitting_ft_faces_;
  Redistribute_Field redistribute_from_splitting_ft_elem_;
  Redistribute_Field redistribute_from_splitting_ft_elem_ghostz_;
  Redistribute_Field redistribute_from_splitting_ft_elem_ghostz_min_;
  Redistribute_Field redistribute_from_splitting_ft_elem_ghostz_max_;
  Redistribute_Field redistribute_to_splitting_ft_elem_;

  OBS_PTR(Milieu_base) le_fluide_;
  Multigrille_Adrien poisson_solver_;

  Noms expression_vitesse_initiale_; // on attend trois expressions
  Nom expression_pression_initiale_; // useless, unless post-pro OR pressure_increment.
  Nom expression_vitesse_upstream_;

  int upstream_dir_ = -1;
  int upstream_stencil_ = 3;
  int upstream_velocity_measured_ = 0;
  int harmonic_nu_in_diff_operator_ = 0;
  int harmonic_nu_in_calc_with_indicatrice_ = 0;
  int vitesse_entree_dir_ = DIRECTION_I;
  int vitesse_entree_compo_to_force_ = -1;
  int stencil_vitesse_entree_ = 3;
  int test_etapes_et_bilan_ = 0;
  int add_initial_field_ = 0;
  int diffusion_alternative_ = 0;
  int suppression_rejetons_ = 0;  // By defaults, break-ups are not fixed on restart. (no deletion of smaller fractions)

  // Discretisation du champ (1/rho) et utilisation dans le calcul de rho*v*v, dans le mass_solver
  // et dans pressure_projection_with_inv_rho :
  int use_inv_rho_for_mass_solver_and_calculer_rho_v_ = 0;
  int use_inv_rho_in_poisson_solver_ = 0;
  int use_inv_rho_ = 0;

  // Pour la premiere projection, on initialise la pression au champ de pression a l'equilibre diphasique :
  int improved_initial_pressure_guess_ = 0;
  // travail en increment de pression pour aider le solveur :
  int include_pressure_gradient_in_ustar_ = 0;
  int correction_bilan_qdm_ = 0;
  int refuse_patch_conservation_QdM_RK3_source_interf_ = 0;

  int disable_solveur_poisson_ = 0;
  int resolution_fluctuations_ = 0;
  int projection_initiale_demandee_ = 0;
  int disable_diffusion_qdm_ = 0;
  int disable_convection_qdm_ = 0;
  int disable_source_interf_ = 0;
  int frozen_velocity_ = 0;
  int velocity_reset_ = 0;

  double nb_diam_upstream_ = 0.;
  double nb_diam_ortho_shear_perio_= -1.1e20;
  double vitesse_entree_ = -1.1e20;
  double vitesse_upstream_ = -1.1e20;
  double vitesse_upstream_reprise_ = -1.1e20;
  double velocity_bubble_new_ = 0.;
  double velocity_bubble_old_ = -1.1e20;
  double velocity_bubble_integral_err_ = 0.;
  double velocity_bubble_scope_ = 0.;
  double upstream_velocity_bubble_factor_ = 1.;
  double upstream_velocity_bubble_factor_deriv_ = 0.;
  double upstream_velocity_bubble_factor_integral_ = 0.;


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

  // right hand side for pressure solver
  IJK_Field_double pressure_rhs_;
  IJK_Field_double pressure_;
  IJK_Field_double pressure_ghost_cells_;
  IJK_Field_vector3_double velocity_;

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

  // Field only needed for the option type_velocity_convection_form_== Nom("non_conservative_rhou")
  IJK_Field_double div_rhou_;

  // Temporary storage for the derivative
  IJK_Field_vector3_double d_velocity_;
  // Temporary storage for the fractional timestep in rk3 :
  IJK_Field_vector3_double RK3_F_velocity_;


  double pression_ap_proj_ = 0.;
  double vap_velocity_tmoy_ = 0.;
  double liq_velocity_tmoy_ = 0.;

  int compute_rising_velocities_ = 0;
  int fill_rising_velocities_ = 0;
  corrections_qdm qdm_corrections_;


  double vol_bulle_monodisperse_ = -1; // Pour imposer le volume des bulles
  double diam_bulle_monodisperse_ = -1; // Pour imposer le volume des bulles
  ArrOfDouble vol_bulles_;   // Le volume impose individuellement a chaque bulle.
  double coeff_evol_volume_ = 0.;


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

  IJK_Field_vector3_double velocity_ft_;

  // Correcteur PID
  double Kp_ = 0.;
  double Kd_ = 0.;
  double Ki_ = 0.;
  double int_x_ = 0.;
  double int_y_ = 0.;
  double int_z_ = 0.;
  int epaisseur_maille_ = 8;

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

  Boundary_Conditions boundary_conditions_;

  // Le jeu de donnees doit fournir soit des fichiers de reprise: ..
  Nom fichier_reprise_vitesse_ = "??"; // par defaut, invalide
  int timestep_reprise_vitesse_ = 1;

  init_forcage_THI forcage_;

  Champs_compris_T<IJK_Field_double> champs_compris_;

  bool flag_variable_source_ = false;
};

#endif /* Navier_Stokes_FTD_IJK_included */
