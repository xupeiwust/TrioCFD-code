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

#include <Fluide_Diphasique_IJK.h>
#include <Schema_Temps_IJK_base.h>
#include <Navier_Stokes_FTD_IJK.h>
#include <Probleme_FT_Disc_gen.h>
#include <IJK_Interfaces.h>
#include <IJK_Thermals.h>
#include <Domaine_IJK.h>
#include <Postprocessing_IJK.h>
#include <Champs_compris.h>

class Domaine_IJK;

class Probleme_FTD_IJK_base : public Probleme_FT_Disc_gen
{
  Declare_base(Probleme_FTD_IJK_base) ;
public :
  // We take too much advantage of it ...:
  friend class IJK_Thermique_cut_cell;
  friend class Statistiques_dns_ijk_FT;
  friend class Postprocessing_IJK;

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
//  int postraiter(int force=1) override;
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

  void initialise_ijk_fields();
  void initialise_interfaces();

  void sauvegarder_probleme(const char *fichier_sauvegarde, const int& stop); //  const;
  void reprendre_probleme(const char *fichier_reprise);

  void update_thermal_properties();

  int get_thermal_probes_ghost_cells() const { return thermal_probes_ghost_cells_; }

  double t_debut_statistiques() const { return get_post().t_debut_statistiques(); }
  bool get_reprise() const { return reprise_; }

  void get_noms_champs_postraitables(Noms& noms,Option opt) const override;
  const Postprocessing_IJK& get_post() const { return ref_cast(Postprocessing_IJK,  les_postraitements_.front().valeur()); }
  Postprocessing_IJK& get_post() { return ref_cast(Postprocessing_IJK,  les_postraitements_.front().valeur()); }

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

  virtual Cut_cell_FT_Disc* get_cut_cell_disc()
  {
    Process::exit("No cut fields are found.");
    throw;
  }

  bool has_ns() const { return has_ns_; }
  bool has_interface() const { return has_interface_; }
  bool has_thermals() const { return has_thermals_; }

  const Navier_Stokes_FTD_IJK& eq_ns() const
  {
    assert (has_ns_);
    if (!has_ns_) { throw; };
    return ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur());
  }
  Navier_Stokes_FTD_IJK& eq_ns()
  {
    assert (has_ns_);
    if (!has_ns_) { throw; };
    return ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur());
  }

  const IJK_Interfaces& get_interface() const
  {
    if (!has_interface_) { return interface_to_remove_later_when_clean_; }; // TODO (teo.boutin) (just in case we forget, and I'm not the one who did that)
    return ref_cast(IJK_Interfaces, equations_[1].valeur());
  }
  IJK_Interfaces& get_interface()
  {
    if (!has_interface_) {  return interface_to_remove_later_when_clean_; };
    return ref_cast(IJK_Interfaces, equations_[1].valeur());
  }

  const Maillage_FT_IJK& get_maillage_ft_ijk() const
  {
    assert (has_interface_);
    if (!has_interface_) { throw; };
    return get_interface().maillage_ft_ijk();
  }
  const Remaillage_FT_IJK& get_remaillage_ft_ijk() const
  {
    assert (has_interface_);
    if (!has_interface_) { throw; };
    return get_interface().remaillage_ft_ijk();
  }

  const IJK_Thermals& get_ijk_thermals() const
  {
    assert (has_thermals_);
    if (has_interface_ && has_thermals_)
      return ref_cast(IJK_Thermals, equations_[2].valeur());
    else if (!has_interface_ && has_thermals_)
      return ref_cast(IJK_Thermals, equations_[1].valeur());
    else
      throw;
  }

  IJK_Thermals& get_ijk_thermals()
  {
    assert (has_thermals_);
    if (has_interface_ && has_thermals_)
      return ref_cast(IJK_Thermals, equations_[2].valeur());
    else if (!has_interface_ && has_thermals_)
      return ref_cast(IJK_Thermals, equations_[1].valeur());
    else
      throw;
  }


protected:
  IJK_Interfaces interface_to_remove_later_when_clean_; // TODO FIXME
  bool has_interface_ = false, has_ns_ = false, has_thermals_ = false;
  ArrOfDouble_with_ghost delta_z_local_;
  OBS_PTR(Domaine_IJK) domaine_ijk_;
  Domaine_IJK domaine_ft_;
  Nom lata_name_;
  bool stop_ = false;

  // Name / Localisation (elem, face, ...) / Nature (scalare, vector) / Needs interpolation
  using FieldInfo_t = std::tuple<Motcle, Entity, Nature_du_champ, bool>;
  std::vector<FieldInfo_t> champs_postraitables_; // list of fields that may be postprocessed

  Nom fichier_post_;  // TODO a virer une fois le clean du post fini
  Nom nom_sauvegarde_, nom_reprise_;
  bool sauvegarder_xyz_ = false; // drapeau 0 ou 1
  bool reprise_ = false;// flag pour indiquer si on fait une reprise

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

  Champs_compris_T<IJK_Field_double> champs_compris_;

  // Champ IJK_Field notant les cellules parcouru lors d'un traitement,
  // c'est-a-dire pour eviter de recalculer plusieurs fois les memes cases lors du calculs des flux.
  // Le champ est public pour faciliter l'utilisation dans IJK_Thermal_cut_cell
  IJK_Field_int treatment_count_;

  // Compteur du dernier traitement effectue dans treatment_count_
  int new_treatment_ = 0;
  int thermal_probes_ghost_cells_ = 2;

  void fill_post_fields();
};

#endif /* Probleme_FTD_IJK_base_included */
