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

#ifndef Postprocessing_IJK_included
#define Postprocessing_IJK_included

#include <IJK_Field.h>
#include <IJK_Field_vector.h>
#include <Motcle.h>
#include <TRUST_Vector.h>
#include <Statistiques_dns_ijk_FT.h>
#include <Statistiques_dns_ijk_monophasique.h>
#include <IJK_Interfaces.h>
#include <Multigrille_Adrien.h>
#include <Postraitement_ft_lata.h>
#include <Champs_compris_IJK_interface.h>
#include <Champs_compris_IJK.h>

class Probleme_FTD_IJK_base;
class Navier_Stokes_FTD_IJK;
class Domaine_IJK;

class IJK_Thermals;

/**
 * Post-processing stuff of Probleme_FTD_IJK_base.
 */
class Postprocessing_IJK: public Postraitement_ft_lata, public Champs_compris_IJK_interface
{
  Declare_instanciable(Postprocessing_IJK);

  friend class Statistiques_dns_ijk_FT;

public:
  // Name / Localisation (elem, face, ...) / Nature (scalare, vector) / Needs interpolation
  using FieldInfo_t = Champs_compris_IJK_interface::FieldInfo_t;

  static std::vector<FieldInfo_t>& Get_champs_postraitables() { return champs_postraitables_; }

  void set_param(Param& param) override;
  void lire_entete_bloc_interface(Entree& is) override;
  int lire_champs_a_postraiter(Entree& is, bool expect_acco) override;
  void register_interface_field(const Motcle& nom_champ, const Motcle& loc) override;

  void init() override;
  void completer() override { /* Does nothing */  }
  void postraiter(int forcer) override;
  int postraiter_champs() override;
  void prepare_lata_and_stats(); // merge with init() ?
  int write_extra_mesh() override;

  void resetTime(double t, std::string dirname) override { /* not impl. */ throw; }

  void associer_probleme(const Probleme_FTD_IJK_base& );

  void associer_domaines(Domaine_IJK& dom_ijk, Domaine_IJK& dom_ft);
  void init_integrated_and_ana(bool reprise);
  void fill_indic(bool reprise=0);
  void initialise_stats(Domaine_IJK& splitting, ArrOfDouble& vol_bulles, const double vol_bulle_monodisperse);
  void init_indicatrice_non_perturbe();

  void posttraiter_champs_instantanes(const char * lata_name, double time, int time_iteration);
  void posttraiter_statistiques_plans(double time);
  void ecrire_statistiques_bulles(int reset, const Nom& nom_cas, const double current_time) const;
  void ecrire_statistiques_cisaillement(int reset, const Nom& nom_cas, const double current_time) const;
  void ecrire_statistiques_rmf(int reset, const Nom& nom_cas, const double current_time) const;
  void update_stat_ft(const double dt);
  void get_update_lambda2();
  void get_update_lambda2_and_rot_and_Q();
  void activate_cut_cell() { cut_cell_activated_ = 1; };

  IJK_Field_double& rebuilt_indic() { return rebuilt_indic_;  }
  IJK_Field_vector3_double& coords()  { return coords_;  }
  IJK_Field_double& integrated_timescale() { return integrated_timescale_; }
  bool postraiter_sous_pas_de_temps() const { return postraiter_sous_pas_de_temps_; }

  int post_par_paires() const { return post_par_paires_; }
  double t_debut_statistiques() const { return t_debut_statistiques_; }

  inline int sondes_demande() { return sondes_demande_; }

  bool is_post_required(const Motcle& nom) const;

  // Interface Champs_compris_IJK_interface
  bool has_champ(const Motcle& nom) const override { return champs_compris_.has_champ(nom);  }
  bool has_champ_vectoriel(const Motcle& nom) const override { return champs_compris_.has_champ_vectoriel(nom); }
  const IJK_Field_vector3_double& get_IJK_field_vector(const Motcle& nom) override;
  const IJK_Field_double& get_IJK_field(const Motcle& nom) override;

  static void Fill_postprocessable_fields(std::vector<FieldInfo_t>& chps);
  void get_noms_champs_postraitables(Noms& noms,Option opt=NONE) const;

  const int& get_IJK_flag(const Nom& nom) const;

  inline IJK_Field_vector3_double& get_grad_I_ns() { return grad_I_ns_; }

  void sauvegarder_post(const Nom& lata_name);
  void sauvegarder_post_maitre(const Nom& lata_name, SFichier& fichier) const;
  void reprendre_post(Param& param);

  void fill_op_conv();
  void fill_surface_force(IJK_Field_vector3_double& the_field_you_know);
  void fill_surface_force_bis(const char * lata_name, double time, int time_iteration);
  IJK_Field_vector3_double get_rho_Ssigma();

  void calculer_gradient_indicatrice_et_pression(const IJK_Field_double& indic);

  void improved_initial_pressure_guess(bool imp);

  void posttraiter_tous_champs_thermique(Motcles& liste,  const int idx) const;
  void posttraiter_tous_champs_energie(Motcles& liste,  const int idx) const;
  void posttraiter_tous_champs_thermal(Motcles& liste, const int idx) const;

  Motcles get_liste_post_instantanes() const { return liste_post_instantanes_; }

  void alloc_fields();
  void alloc_velocity_and_co();

  void compute_extended_pressures();

  bool is_stats_bulles_activated() const;
  bool is_stats_plans_activated() const;
  bool is_stats_cisaillement_activated() const;
  bool is_stats_rmf_activated() const;


  double get_max_timestep_for_post(double current_time) const;

protected:
  static std::vector<FieldInfo_t> champs_postraitables_;  ///< list of fields that can be potentially postprocessed
  /** Index in 'champs_postraitables_' of each of the requested field for post-processing
   * and flag indicating if interpolation will be needed:
   */
  std::vector<std::pair<int,bool>> field_post_idx_;

  std::vector<Motcle> list_post_required_;

  Champs_compris_IJK champs_compris_;  ///< the actual fields registered and managed by the post-processing part (=all the extra fields, not the main unknowns)

  // Storage of all the extra fields created for post processing:
  std::map<Motcle, IJK_Field_double> scalar_post_fields_;
  std::map<Motcle, IJK_Field_vector3_double> vect_post_fields_;

  void compute_phase_pressures_based_on_poisson(const int phase);
  Statistiques_dns_ijk_FT statistiques_FT_;

  // Post-traitement selon un nombre de pas de temps
  // The main postprocessing freq (nb_pas_dt_post) is in the base class
  int nb_pas_dt_post_thermals_probes_ = -1;
  int nb_pas_dt_post_stats_bulles_ = -1; // intervalle de posttraitement des donnees par bulles
  int nb_pas_dt_post_stats_plans_ = -1; // intervalle de posttraitement des donnees par plan (pour les statistiques de canal)
  int nb_pas_dt_post_stats_cisaillement_ = -1; // intervalle de posttraitement des données liés au cisaillement
  int nb_pas_dt_post_stats_rmf_ = -1; // intervalle de posttraitement des données liés au au rmf

  // Post-traitement selon un intervale de temps (en secondes)
  double time_interval_post_ = -1.0;
  double time_interval_post_thermals_probes_ = -1.0;
  double time_interval_post_stats_bulles_ = -1.0;
  double time_interval_post_stats_plans_ = -1.0;
  double time_interval_post_stats_cisaillement_ = -1.0;
  double time_interval_post_stats_rmf_ = -1.0;

  Motcles liste_post_instantanes_; // liste des champs instantanes a postraiter
  // Pour numeroter les fichiers .lata il faut compter combien on en a ecrit:
  int compteur_post_instantanes_ = 0;
  int postraiter_sous_pas_de_temps_ = 0; // drapeau 0 ou 1
  // Pour reconstruire au post-traitement la grandeur du/dt, on peut choisir de relever u^{dt_post} et u^{dt_post+1} :
  int post_par_paires_ = 0; // drapeau 0 ou 1

  // Pour des fiches de validation, on post-traite le champ analytique attendu dans le lata pour calcul de l'erreur:
  Noms expression_vitesse_analytique_; // on attend trois expressions
  Nom expression_pression_analytique_ = "??"; // par defaut, invalide, on attend une expression
  Noms expression_dvitesse_analytique_; // on attend trois expressions

  // Pour check_stats (and_grads)
  int check_stats_ = 0;
  Noms expression_gradP_analytique_; // on attend trois expressions
  Noms expression_gradU_analytique_; // on attend trois expressions
  Noms expression_gradV_analytique_; // on attend trois expressions
  Noms expression_gradW_analytique_; // on attend trois expressions
  // Second gradient (laplacian_P or Hessienne
  Noms expression_grad2P_analytique_; // on attend 6 expressions (car tens sym)
  // And for the 3 components of velocity :
  Noms expression_grad2U_analytique_; // on attend 6 expressions (car tens sym)
  Noms expression_grad2V_analytique_; // on attend 6 expressions (car tens sym)
  Noms expression_grad2W_analytique_; // on attend 6 expressions (car tens sym)

  // -------------------------------------------------
  // Statistiques temporelles
  // Pour les groupes de bulles : on cree un vecteur d'objet stat de taille dimensionnee a 0 puis nb_groups():
  VECT(Statistiques_dns_ijk_FT) groups_statistiques_FT_;
  //Statistiques_dns_ijk_FT statistiques_FT_;
  //
  // 2020.03.12. CHOIX : Meme en disable_diphasique, on fait appel a la classe fille stats FT.
  // La classe de stats monophasique n'est plus maintenue. Suppression du membre.
  // Statistiques_dns_ijk_monophasique statistiques_;
  double t_debut_statistiques_ = -1.0;
  // -------------------------------------------------

  // Pour les cas a bulles fixes
  IJK_Field_vector3_double integrated_velocity_;
  IJK_Field_double integrated_pressure_;
  IJK_Field_double indicatrice_non_perturbe_;
  IJK_Field_double integrated_timescale_;

  bool reset_reprise_integrated_ = false;
  // Pour la reprise bulles fixes, parametres de lecture de champ de condition initiale pour variables de post-tt:
  Nom fichier_reprise_integrated_velocity_ = "??";
  Nom fichier_reprise_integrated_pressure_ = "??";
  Nom fichier_reprise_indicatrice_non_perturbe_ = "??";
  Nom fichier_reprise_integrated_timescale_ = "??";

  // Temporary storage for the coords (for postprocessing) :
  IJK_Field_vector3_double coords_;
  IJK_Field_vector3_double velocity_ana_;
  IJK_Field_vector3_double ecart_ana_;
  IJK_Field_vector3_double op_conv_;
  IJK_Field_vector3_double cell_op_conv_;
  IJK_Field_vector3_double cell_rho_Ssigma_;

  IJK_Field_vector3_double d_velocity_ana_;
  IJK_Field_double ecart_p_ana_;

  // Celui la est discretise sur le maillage etendu:
  IJK_Field_vector3_double grad_I_ft_;

  // Pour postraitement :
  IJK_Field_double rebuilt_indic_;
  int extended_pressure_computed_ = 0;
  IJK_Field_double pressure_ft_;
  IJK_Field_double extended_pl_ft_;
  IJK_Field_double extended_pv_ft_;
  IJK_Field_double extended_pl_;
  IJK_Field_double extended_pv_;
  // Pour le calcul des stats  :
  IJK_Field_double kappa_ai_ft_;
  IJK_Field_vector3_double normale_cell_ft_;
  IJK_Field_double ai_ns_;
  IJK_Field_double kappa_ai_ns_;
  IJK_Field_vector3_double normale_cell_ns_;
  // For lambda and curl
  IJK_Field_double dudx_;
  IJK_Field_double dvdy_;
  IJK_Field_double dwdx_;
  IJK_Field_double dudz_;
  IJK_Field_double dvdz_;
  IJK_Field_double dwdz_;
//  IJK_Field_double lambda2_;
//  IJK_Field_double critere_Q_;
//  IJK_Field_vector3_double rot_;
  IJK_Field_vector3_double grad_I_ns_;
  IJK_Field_vector3_double grad_P_;
//  IJK_Field_double num_compo_ft_;

  // Pour la verification des stats :
  // Le gradient de pression aux faces :
  //  IJK_Field_vector3_double gradP_;
  // Les gradients des compo de vitesses aux elems : (sont finalement stockes dans statistiques_FT_ si besoin)
  //IJK_Field_vector3_double dUd_;
  //IJK_Field_vector3_double dVd_;
  //IJK_Field_vector3_double dWd_;
  // Et leurs solutions analytiques :
  IJK_Field_vector3_double ana_gradP_;
  IJK_Field_vector3_double ana_dUd_;
  IJK_Field_vector3_double ana_dVd_;
  IJK_Field_vector3_double ana_dWd_;
  // Pour les deriv secondes :
  IJK_Field_vector3_double ana_grad2Pi_; // Partie diagonale de la jacobienne
  IJK_Field_vector3_double ana_grad2Pc_; // contient les deriv croisees
  IJK_Field_vector3_double ana_grad2Ui_; // Partie diagonale de la jacobienne
  IJK_Field_vector3_double ana_grad2Uc_; // contient les deriv croisees
  IJK_Field_vector3_double ana_grad2Vi_; // Partie diagonale de la jacobienne
  IJK_Field_vector3_double ana_grad2Vc_; // contient les deriv croisees
  IJK_Field_vector3_double ana_grad2Wi_; // Partie diagonale de la jacobienne
  IJK_Field_vector3_double ana_grad2Wc_; // contient les deriv croisees

  IJK_Field_double IFd_source_spectraleX_;
  IJK_Field_double AOD_source_spectraleX_;
  IJK_Field_double source_spectraleY_;
  IJK_Field_double source_spectraleZ_;
  // Pour post-traitement :
  IJK_Field_double dudy_, dvdx_, dwdy_;
  IJK_Field_vector3_double cell_velocity_;
  IJK_Field_vector3_double cell_source_spectrale_;
  IJK_Field_vector3_double cell_bk_tsi_ns_;
  //  IJK_Field_vector3_double cell_source_interface_totale_;   // non-const because some echange_espace_virtuel()
  IJK_Field_vector3_double cell_grad_p_;
  IJK_Field_vector3_double cell_source_interface_; // toujours en _ns_
  IJK_Field_vector3_double cell_backup_source_interface_; // toujours en _ns_
  IJK_Field_vector3_double cell_repulsion_interface_; // toujours en _ns_


  int sondes_demande_ = 0;

  /*
   * References to various members of Probleme_FTD_IJK_base that are heavily used in the post:
   */
  OBS_PTR( Probleme_FTD_IJK_base) ref_ijk_ft_;

  OBS_PTR( IJK_Interfaces) interfaces_;
  OBS_PTR(IJK_Field_double) pressure_;                   // non-const because some echange_espace_virtuel()
  OBS_PTR(IJK_Field_vector3_double) velocity_;   // non-const because some echange_espace_virtuel()
  OBS_PTR(IJK_Field_vector3_double) source_spectrale_;   // non-const because some echange_espace_virtuel()
  OBS_PTR(IJK_Field_vector3_double) bk_tsi_ns_;
  IJK_Field_vector3_double source_interface_ft_;   // non-const because some echange_espace_virtuel()
  IJK_Field_vector3_double source_interface_ns_;   // non-const because some echange_espace_virtuel()
  IJK_Field_vector3_double repulsion_interface_ns_;   // non-const because some echange_espace_virtuel()
  OBS_PTR(IJK_Field_vector3_double) d_velocity_;

  OBS_PTR(Domaine_IJK) domaine_ijk_;
  OBS_PTR(Domaine_IJK) domaine_ft_;
  OBS_PTR(IJK_Thermals) thermals_;
  int first_step_thermals_post_=0;

  Multigrille_Adrien poisson_solver_post_;

  // Pour le post-traitement des champs cut-cell
  int cut_cell_activated_ = 0;

  void register_one_field(const Motcle& fld_nam, const Motcle& loc);

private:
  IJK_Field_vector3_double post_projected_field_; ///< Temporary storage space used when invoking 'interpolate_to_center'

  void postraiter_thermals(bool stop);

  void postraiter_stats(bool stop);
};

#endif /* Postprocessing_IJK_included */
