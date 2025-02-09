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
#include <Schema_Euler_explicite_IJK.h>
#include <Schema_RK3_IJK.h>
#include <Navier_Stokes_FTD_IJK.h>

class Domaine_IJK;

class Probleme_FTD_IJK_base : public Probleme_FT_Disc_gen
{
  Declare_base_sans_constructeur(Probleme_FTD_IJK_base) ;
public :
  // We take too much advantage of it ...:
  friend class IJK_Thermique;
  friend class IJK_Thermique_cut_cell;
  friend class Statistiques_dns_ijk_FT;
  friend class IJK_FT_Post;

  Probleme_FTD_IJK_base();
  Probleme_FTD_IJK_base(const Probleme_FTD_IJK_base& x);
  int associer_(Objet_U&) override;

  virtual void set_param(Param& param);

  const Domaine_IJK& get_domaine_ft() const { return domaine_ft_ ; }
  Domaine_IJK& get_domaine_ft() { return domaine_ft_ ; } // non-constpour NS :: // TODO FIXME // associer comme avant !

  const Domaine_IJK& get_domaine() const { return domaine_ijk_.valeur(); }
  const Domaine_IJK& domaine_ijk() const { return domaine_ijk_.valeur(); }
  Domaine_IJK& domaine_ijk() { return domaine_ijk_.valeur(); }

  void discretiser(Discretisation_base& dis) override;
  void preparer_calcul() override;

  void typer_lire_milieu(Entree& is) override;
  void lire_solved_equations(Entree& is) override;

  void completer() override;
  void mettre_a_jour(double temps) override { }
  virtual bool updateGivenFields() override { return false; }

  inline Fluide_Diphasique_IJK& milieu_ijk() { return ref_cast(Fluide_Diphasique_IJK, le_milieu_.front().valeur()); }
  inline const Fluide_Diphasique_IJK& milieu_ijk() const { return ref_cast(Fluide_Diphasique_IJK, le_milieu_.front().valeur()); }
  const Schema_Temps_IJK_base& schema_temps_ijk() const { return ref_cast(Schema_Temps_IJK_base, le_schema_en_temps_.valeur()); }
  Schema_Temps_IJK_base& schema_temps_ijk() { return ref_cast(Schema_Temps_IJK_base, le_schema_en_temps_.valeur()); }

  // interface Problem
  bool run() override;
  void initialize() override;
  void terminate() override;
  int postraiter(int force=1) override;
  void sauver() const override;

  // interface UnsteadyProblem
  double presentTime() const override { return 0.; }
  double computeTimeStep(bool& stop) const override;
  bool initTimeStep(double dt) override { return true; }
  bool solveTimeStep() override;
  void validateTimeStep() override;
  void setStationary(bool flag) override {  }
  void abortTimeStep() override {  }
  void resetTime(double time) override { }

  // interface IterativeUnsteadyProblem
  bool iterateTimeStep(bool& converged) override { return iterateTimeStep_impl(*this, converged); }

  Entree& lire_equations(Entree& is, Motcle& dernier_mot) override { return is; }

  IJK_Field_int& treatment_count() { return treatment_count_; }
  const IJK_Field_int& treatment_count() const { return treatment_count_; }
  int new_treatment() const { return new_treatment_; }
  int& new_treatment() { return new_treatment_; }

  const Nom& nom_sauvegarde() const { return nom_sauvegarde_; }

  const IJK_Field_double& get_IJK_field(const Nom& nom) const;

  int initialise_ijk_fields();
  int initialise_interfaces();

  void sauvegarder_probleme(const char *fichier_sauvegarde, const int& stop); //  const;
  void reprendre_probleme(const char *fichier_reprise);

  void update_thermal_properties();

  const LIST(IJK_Thermique)& get_list_thermique() const { return thermique_; }
  LIST(IJK_Thermique)& get_list_thermique() { return thermique_; }

  const LIST(IJK_Energie)& get_list_energie() const { return energie_; }
  LIST(IJK_Energie)& get_list_energie() { return energie_; }

  const IJK_Thermals& get_ijk_thermals() const { return thermals_; }
  IJK_Thermals& get_ijk_thermals() { return thermals_; }

  int get_thermal_probes_ghost_cells() const { return thermal_probes_ghost_cells_; }

  double t_debut_statistiques() const { return post_.t_debut_statistiques(); }
  int get_reprise() const { return reprise_; }
  const IJK_FT_Post& get_post() const { return post_; }
  IJK_FT_Post& get_post() { return post_; }

  const Maillage_FT_IJK& get_maillage_ft_ijk() const { return interfaces_.maillage_ft_ijk(); }
  const Remaillage_FT_IJK& get_remaillage_ft_ijk() const { return interfaces_.remaillage_ft_ijk(); }

  const IJK_Interfaces& get_interface() const { return interfaces_; }
  const IJK_Interfaces& itfce() const { return interfaces_; }
  IJK_Interfaces& get_set_interface() { return interfaces_; }

  ArrOfDouble_with_ghost& get_delta_z_local() { return delta_z_local_; }
  const ArrOfDouble_with_ghost& get_delta_z_local() const { return delta_z_local_; }

  virtual void deplacer_interfaces(const double timestep, const int rk_step, ArrOfDouble& var_volume_par_bulle, const int first_step_interface_smoothing);
  virtual void deplacer_interfaces_rk3(const double timestep, const int rk_step, ArrOfDouble& var_volume_par_bulle);
  virtual void update_indicator_field();
  virtual void update_twice_indicator_field();
  void update_pre_remeshing_indicator_field();
  void update_post_remeshing_indicator_field();
  void update_old_intersections();
  void parcourir_maillage();
  void calculer_rho_mu_indicatrice(const bool parcourir = true);

  const Nom& get_lata_name() const { return lata_name_; }

  const Probleme_base& probleme(const Domaine_IJK& dom) const
  {
    if (dom == domaine_ft_)
      return refprobleme_ft_disc_.valeur();
    Process::exit("Unrecognized domain provided");
    throw;
  }
  const Probleme_FTD_IJK_base& operator=(const Probleme_FTD_IJK_base& a) { throw ; }

  virtual Cut_cell_FT_Disc* get_cut_cell_disc()
  {
    Process::exit("No cut fields are found.");
    throw;
  }

  const Navier_Stokes_FTD_IJK& eq_ns() const { return ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur()); }
  Navier_Stokes_FTD_IJK& eq_ns() { return ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur()); }

protected:
  ArrOfDouble_with_ghost delta_z_local_;
  IJK_Interfaces interfaces_;
  OBS_PTR(Domaine_IJK) domaine_ijk_;
  Domaine_IJK domaine_ft_;
  Nom lata_name_;
  bool stop_ = false;

  IJK_FT_Post post_;
  Nom fichier_post_, nom_sauvegarde_, nom_reprise_;
  int sauvegarder_xyz_ = 0; // drapeau 0 ou 1
  int reprise_ = 0;// flag pour indiquer si on fait une reprise

  // Le probleme ft disc qui porte le maillage vdf pour les algorithmes front-tracking
  OBS_PTR(Probleme_base) refprobleme_ft_disc_;
  // Creation d'un probleme sur le domaine d'origine pour les sondes et pour faire leur VDF...
  OBS_PTR(Probleme_base) refprobleme_ns_;

  void euler_time_step(ArrOfDouble& var_volume_par_bulle);
  void rk3_sub_step(const int rk_step, const double total_timestep, const double fractionnal_timestep, const double time);
  virtual void create_forced_dilation() { }
  void solveTimeStep_Euler(DoubleTrav&);
  void solveTimeStep_RK3(DoubleTrav&);
  void build_vdf_domaine();

  // Champ IJK_Field notant les cellules parcouru lors d'un traitement,
  // c'est-a-dire pour eviter de recalculer plusieurs fois les memes cases lors du calculs des flux.
  // Le champ est public pour faciliter l'utilisation dans IJK_Thermal_cut_cell
  IJK_Field_int treatment_count_;

  // Compteur du dernier traitement effectue dans treatment_count_
  int new_treatment_ = 0;

  /*
   * XXX WILL BE REMOVED
   */
  LIST(IJK_Thermique) thermique_;
  LIST(IJK_Energie) energie_;

  /*
   * PRIORITE
   */
  IJK_Thermals thermals_;
  int thermal_probes_ghost_cells_ = 2;
};

#endif /* Probleme_FTD_IJK_base_included */
