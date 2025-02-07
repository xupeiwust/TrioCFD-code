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

#ifndef IJK_FT_Post_included
#define IJK_FT_Post_included

#include <IJK_Field.h>
#include <IJK_Field_vector.h>
#include <Motcle.h>
#include <TRUST_Vector.h>
#include <Statistiques_dns_ijk_FT.h>
#include <Statistiques_dns_ijk_monophasique.h>
#include <Sondes_IJK.h>
#include <IJK_Interfaces.h>
#include <Multigrille_Adrien.h>

class Probleme_FTD_IJK_base;
class Navier_Stokes_FTD_IJK;
class Domaine_IJK;

/*
 * TODO: Demander à Aymeric l'interet (obsolete ??)
 */
class IJK_Thermique;
class IJK_Energie;
class IJK_Thermals;

/**
 * All the post-processing stuff of Probleme_FTD_IJK_base delegated into this helper class:
 */
class IJK_FT_Post
{
  friend class Statistiques_dns_ijk_FT;

public:
  IJK_FT_Post();
  void complete_interpreter(Param& param, Entree& e);
  void associer_probleme(const Probleme_FTD_IJK_base& );

  void associer_domaines(Domaine_IJK& dom_ijk, Domaine_IJK& dom_ft);
  int initialise(int reprise);
  void complete(int reprise);
  int initialise_stats(Domaine_IJK& splitting, ArrOfDouble& vol_bulles, const double vol_bulle_monodisperse);
  void init_indicatrice_non_perturbe();

  void posttraiter_champs_instantanes(const char * lata_name, double time, int time_iteration);
  void posttraiter_statistiques_plans(double time);
  void ecrire_statistiques_bulles(int reset, const Nom& nom_cas, const DoubleTab& gravite, const double current_time) const;
  void ecrire_statistiques_cisaillement(int reset, const Nom& nom_cas, const double current_time) const;
  void ecrire_statistiques_rmf(int reset, const Nom& nom_cas, const double current_time) const;
  void update_stat_ft(const double dt);
  void get_update_lambda2();
  void get_update_lambda2_and_rot_and_curl();
  void activate_cut_cell() { cut_cell_activated_ = 1; };

  IJK_Field_double& rebuilt_indic()
  {
    return rebuilt_indic_;
  }
  IJK_Field_double& potentiel()
  {
    return potentiel_;
  }
  IJK_Field_vector3_double& coords()
  {
    return coords_;
  }
  IJK_Field_double& integrated_timescale()
  {
    return integrated_timescale_;
  }
  bool postraiter_sous_pas_de_temps() const
  {
    return postraiter_sous_pas_de_temps_;
  }
  double get_timestep_simu_post(double current_time, double max_simu_time) const
  {
    // Note : the (1+1e-12) safety factor ensures that the simulation reaches the target.
    // Otherwise, the simulation time might fall just below the target due to numerical errors, not triggering the desired post.
    double max_simu_timestep = (max_simu_time - current_time)*(1+1e-12);
    double max_post_timestep                 = ((std::floor(current_time/time_interval_post_) + 1)*time_interval_post_ - current_time)*(1+1e-12);
    double max_post_thermals_probes_timestep = ((std::floor(current_time/time_interval_post_thermals_probes_) + 1)*time_interval_post_thermals_probes_ - current_time)*(1+1e-12);
    double max_post_stats_bulles_timestep    = ((std::floor(current_time/time_interval_post_stats_bulles_) + 1)*time_interval_post_stats_bulles_ - current_time)*(1+1e-12);
    double max_post_stats_plans_timestep     = ((std::floor(current_time/time_interval_post_stats_plans_) + 1)*time_interval_post_stats_plans_ - current_time)*(1+1e-12);
    double max_post_stats_cisaillement_timestep     = ((std::floor(current_time/time_interval_post_stats_cisaillement_) + 1)*time_interval_post_stats_cisaillement_ - current_time)*(1+1e-12);
    double max_post_stats_rmf_timestep     = ((std::floor(current_time/time_interval_post_stats_rmf_) + 1)*time_interval_post_stats_rmf_ - current_time)*(1+1e-12);
    if (max_post_timestep == 0)
      {
        max_post_timestep = max_simu_timestep;
      }
    if (max_post_thermals_probes_timestep == 0)
      {
        max_post_thermals_probes_timestep = max_simu_timestep;
      }
    if (max_post_stats_bulles_timestep == 0)
      {
        max_post_stats_bulles_timestep = max_simu_timestep;
      }
    if (max_post_stats_plans_timestep == 0)
      {
        max_post_stats_plans_timestep = max_simu_timestep;
      }
    if (max_post_stats_cisaillement_timestep == 0)
      {
        max_post_stats_cisaillement_timestep = max_simu_timestep;
      }
    if (max_post_stats_rmf_timestep == 0)
      {
        max_post_stats_rmf_timestep = max_simu_timestep;
      }

    return std::min(max_simu_timestep, std::min(max_post_timestep, std::min(max_post_thermals_probes_timestep, std::min(max_post_stats_plans_timestep, std::min(max_post_stats_bulles_timestep, std::min(max_post_stats_cisaillement_timestep, max_post_stats_rmf_timestep))))));
  }
  int dt_post() const
  {
    return dt_post_;
  }
  int post_par_paires() const
  {
    return post_par_paires_;
  }
  double t_debut_statistiques() const
  {
    return t_debut_statistiques_;
  }

  inline int sondes_demande()
  {
    return sondes_demande_;
  }
  const IJK_Field_double& get_IJK_field(const Nom& nom) const;
  const int& get_IJK_flag(const Nom& nom) const;

  const IJK_Field_vector3_double& get_IJK_vector_field(const Nom& nom) const;

  inline IJK_Field_vector3_double& get_grad_I_ns()
  {
    return grad_I_ns_;
  }

  void sauvegarder_post(const Nom& lata_name);
  void sauvegarder_post_maitre(const Nom& lata_name, SFichier& fichier) const;
  void reprendre_post(Param& param);

  void fill_op_conv();
  void fill_surface_force(IJK_Field_vector3_double& the_field_you_know);//const Nom lata_name, double instant, int iteration);
  void fill_surface_force_bis(const char * lata_name, double time, int time_iteration);
  IJK_Field_vector3_double get_rho_Ssigma();

  void calculer_gradient_indicatrice_et_pression(const IJK_Field_double& indic);

  // Part of the run() method in Probleme_FTD_IJK_base:
  int alloc_fields();
  int alloc_velocity_and_co(bool flag_variable_source);
  void completer_sondes();
  void postraiter_sondes();
  void improved_initial_pressure_guess(bool imp);
  void postraiter_ci(const Nom& lata_name, const double current_time);
  void postraiter_fin(bool stop, int tstep, const int& tstep_init, double current_time, double timestep, const Nom& lata_name,
                      const DoubleTab& gravite, const Nom& nom_cas);
  //void ijk_interpolate_implementation_bis(const IJK_Field_double& field, const DoubleTab& coordinates, ArrOfDouble& result,
  //                                        int skip_unknown_points, double value_for_bad_points,const IJK_Field_double& indic);
  //void  ijk_interpolate_skip_unknown_points_bis(const IJK_Field_double& field, const DoubleTab& coordinates, ArrOfDouble& result,
  //                                              const double value_for_bad_points,const IJK_Field_double& indic);
  void compute_extended_pressures(const Maillage_FT_IJK& mesh);
  //IJK_Field_double& extended_p);
  /*
   * TODO:
   */
  void posttraiter_tous_champs_thermique(Motcles& liste,  const int idx) const;
  void posttraiter_tous_champs_energie(Motcles& liste,  const int idx) const;
  void posttraiter_tous_champs_thermal(Motcles& liste, const int idx) const;

  /*
   * TODO:
   */
  int posttraiter_champs_instantanes_thermique(const Motcles& liste_post_instantanes,
                                               const char *lata_name,
                                               const int lata_step, const double current_time,
                                               IJK_Thermique& itr, const int idx);
  int posttraiter_champs_instantanes_energie(const Motcles& liste_post_instantanes,
                                             const char *lata_name,
                                             const int lata_step, const double current_time,
                                             IJK_Energie& itr,  const int idx);
  /*
   * TODO:
   */
  int posttraiter_champs_instantanes_thermique_interfaciaux(const Motcles& liste_post_instantanes,
                                                            const char *lata_name,
                                                            const int lata_step, const double current_time,
                                                            IJK_Thermique& ,  const int idx);
  int posttraiter_champs_instantanes_energie_interfaciaux(const Motcles& liste_post_instantanes,
                                                          const char *lata_name,
                                                          const int lata_step, const double current_time,
                                                          IJK_Energie& ,  const int idx);

//  void calculer_gradient_temperature(const IJK_Field_double& temperature, IJK_Field_vector3_double& grad_T);

  Motcles get_liste_post_instantanes() const
  {
    return liste_post_instantanes_;
  }
protected:
  void compute_phase_pressures_based_on_poisson(const int phase);
  Statistiques_dns_ijk_FT statistiques_FT_;

  // Post-traitement selon un nombre de pas de temps
  int dt_post_ = 100;
  int dt_post_thermals_probes_ = 100;
  int dt_post_stats_bulles_ = 1; // intervalle de posttraitement des donnees par bulles
  int dt_post_stats_plans_ = 1; // intervalle de posttraitement des donnees par plan (pour les statistiques de canal)
  int dt_post_stats_cisaillement_ = 100; // intervalle de posttraitement des données liés au cisaillement
  int dt_post_stats_rmf_ = 100; // intervalle de posttraitement des données liés au au rmf

  // Post-traitement selon un intervale de temps (en secondes)
  double time_interval_post_ = DMAXFLOAT;
  double time_interval_post_thermals_probes_ = DMAXFLOAT;
  double time_interval_post_stats_bulles_ = DMAXFLOAT;
  double time_interval_post_stats_plans_ = DMAXFLOAT;
  double time_interval_post_stats_cisaillement_ = DMAXFLOAT;
  double time_interval_post_stats_rmf_ = DMAXFLOAT;

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
  double t_debut_statistiques_ = 1.e20;
  // -------------------------------------------------

  // Pour les cas a bulles fixes
  IJK_Field_vector3_double integrated_velocity_;
  IJK_Field_double integrated_pressure_;
  IJK_Field_double indicatrice_non_perturbe_;
  IJK_Field_double integrated_timescale_;
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
  IJK_Field_vector3_double rho_Ssigma_;
  IJK_Field_vector3_double cell_rho_Ssigma_;

  IJK_Field_vector3_double d_velocity_ana_;
  IJK_Field_double pressure_ana_,ecart_p_ana_;

  // Celui la est discretise sur le maillage etendu:
  IJK_Field_vector3_double grad_I_ft_;

  // Pour postraitement :
  IJK_Field_double rebuilt_indic_;
  IJK_Field_double potentiel_;
  IJK_Field_double ai_ft_;
  int extended_pressure_computed_ = 0;
  IJK_Field_double pressure_ft_;
  IJK_Field_double extended_pl_ft_;
  IJK_Field_double extended_pv_ft_;
  IJK_Field_double extended_pl_;
  IJK_Field_double extended_pv_;
  //For the liquid pressure gradient
  // FixedVector<IJK_Field_double 3> dP_ft_;
  // FixedVector<IJK_Field_double 3> dP_;
  // Pour le calcul des stats  :
  IJK_Field_double kappa_ai_ft_;
  IJK_Field_vector3_double normale_cell_ft_;
  IJK_Field_double ai_ns_;
  IJK_Field_double kappa_ai_ns_;
  IJK_Field_vector3_double normale_cell_ns_;
  // The following three fields are needed too for the gradient extension
// /IJK_Field_double dudy_;
  //IJK_Field_double dvdx_;//
  //IJK_Field_double dwdy_;
  // For lambda and curl
  IJK_Field_double dudx_;
  IJK_Field_double dvdy_;
  IJK_Field_double dwdx_;
  IJK_Field_double dudz_;
  IJK_Field_double dvdz_;
  IJK_Field_double dwdz_;
  IJK_Field_double critere_Q_;
  IJK_Field_vector3_double rot_;
  IJK_Field_vector3_double grad_I_ns_;
  IJK_Field_vector3_double grad_P_;
  //  IJK_Field_vector3_double grad_U_ns_;
  //  IJK_Field_vector3_double grad_V_ns_;
  //  IJK_Field_vector3_double grad_W_ns_;
  IJK_Field_double num_compo_ft_;

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

  // GAB
  IJK_Field_double IFd_source_spectraleX_;
  IJK_Field_double AOD_source_spectraleX_;
  IJK_Field_double source_spectraleY_;
  IJK_Field_double source_spectraleZ_;
  // Pour post-traitement :
  IJK_Field_double lambda2_, dudy_, dvdx_, dwdy_;
  IJK_Field_vector3_double cell_velocity_;
  IJK_Field_vector3_double cell_source_spectrale_;
  IJK_Field_vector3_double cell_bk_tsi_ns_;
  //  IJK_Field_vector3_double cell_source_interface_totale_;   // non-const because some echange_espace_virtuel()
  IJK_Field_vector3_double cell_grad_p_;
  IJK_Field_vector3_double cell_source_interface_; // toujours en _ns_
  IJK_Field_vector3_double cell_backup_source_interface_; // toujours en _ns_
  IJK_Field_vector3_double cell_repulsion_interface_; // toujours en _ns_


  int sondes_demande_ = 0;
  Sondes_IJK les_sondes_;  // Sondes a traiter

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
  OBS_PTR(LIST(IJK_Thermique)) thermique_;
  OBS_PTR(LIST(IJK_Energie)) energie_;
  OBS_PTR(IJK_Thermals) thermals_;
  int first_step_thermals_post_=0;

  /* IJK_Field_double temperature_ana_, ecart_t_ana_;
    Nom expression_T_ana_;
    IJK_Field_double source_temperature_ana_, ecart_source_t_ana_; */
  // IJK_Field_vector3_double grad_T_;

  Multigrille_Adrien poisson_solver_post_;

  // Pour le post-traitement des champs cut-cell
  int cut_cell_activated_ = 0;
};


#endif /* IJK_FT_Post_included */
