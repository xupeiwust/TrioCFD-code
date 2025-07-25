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

#ifndef IJK_Interfaces_included
#define IJK_Interfaces_included

#include <Connectivite_frontieres.h>
#include <IJK_Field_vector.h>
#include <IJK_Field.h> // est-ce que j'en ai vraiment besoin ?
#include <Linear_algebra_tools_impl.h>
#include <Maillage_FT_IJK.h>
#include <Equation_base.h>
#include <Parcours_interface.h>
#include <Remaillage_FT_IJK.h>
#include <SFichier.h>
#include <Vecteur3.h>
#include <Linear_algebra_tools_impl.h>
#include <SurfaceVapeurIJKComputation.h>
#include <ComputeValParCompoInCell.h>
#include <TRUST_Ref.h>
#include <Intersection_Interface_ijk.h>
#include <IJK_Composantes_Connex.h>
#include <TRUSTTabFT_cut_cell.h>
#include <Cut_field.h>
#include <Cut_cell_surface_efficace.h>
#include <Domaine_dis_base.h>
#include <Option_IJK.h>
#include <Champs_compris_IJK_interface.h>
#include <Champs_compris_IJK.h>

class Probleme_FTD_IJK_base;
class Switch_FT_double;

#define VERIF_INDIC 0

/*! @brief : class IJK_Interfaces
 *
 *  Cette classe rassemble tous les algorithmes de gestion des interfaces pour le ijk
 *  (le maillage, les algo de remaillage, sauvegarde, reprise, etc)
 *
 */
class IJK_Interfaces : public Equation_base, public Champs_compris_IJK_interface
{
  Declare_instanciable(IJK_Interfaces);
  friend class IJK_Composantes_Connex;
public :
  /*
   * Surcharge de l'eq base
   */
  void set_param(Param& titi) override { /* Do nothing */ }
  void set_param_reprise_pb(Param& );

  void completer() override { /* Do nothing */ }
  void associer_pb_base(const Probleme_base&) override;
  void discretiser() override { /* Do nothing */ }
  int preparer_calcul() override { return 1; }
  const Milieu_base& milieu() const override;
  Milieu_base& milieu() override;
  void associer_milieu_base(const Milieu_base& ) override { /* Do nothing */ }
  int nombre_d_operateurs() const override { return -123; }
  const Operateur& operateur(int) const override { throw; }
  Operateur& operateur(int) override { throw; }
  const Champ_Inc_base& inconnue() const override { throw; }
  Champ_Inc_base& inconnue() override { throw; }

  Probleme_FTD_IJK_base& probleme_ijk();
  const Probleme_FTD_IJK_base& probleme_ijk() const;

  // Interface Champs_compris_IJK_interface
  bool has_champ(const Motcle& nom) const override  { return champs_compris_.has_champ(nom); }
  bool has_champ(const Motcle& nom, OBS_PTR(Champ_base)& ref_champ) const override { /* not used */ throw; }
  bool has_champ_vectoriel(const Motcle& nom) const override { return champs_compris_.has_champ_vectoriel(nom); }
  const IJK_Field_double& get_IJK_field(const Motcle& nom) override;
  const IJK_Field_vector3_double& get_IJK_field_vector(const Motcle& nom) override;
  static void Fill_postprocessable_fields(std::vector<FieldInfo_t>& chps);
  void get_noms_champs_postraitables(Noms& noms,Option opt=NONE) const override;

  void register_fields();

  void dumplata_ft_mesh(const char *filename, const char *meshname, int step) const;

  void initialize(const Domaine_IJK& splitting_FT,
                  const Domaine_IJK& splitting_NS,
                  const Domaine_dis_base& domaine_dis,
                  const int thermal_probes_ghost_cells=0,
                  const bool compute_vint=true,
                  const bool is_switch=false);
  void associer_switch(const Switch_FT_double& ijk_ft_switch);
  void posttraiter_tous_champs(Motcles& liste) const;
  int posttraiter_champs_instantanes(const Motcles& liste_post_instantanes,
                                     const char *lata_name,
                                     const int lata_step) const;
  void sauvegarder_interfaces(const char *lata_name, const Nom& interf_name="??"); // const;
  void calculer_color(ArrOfInt& color) const;
  void postraiter_colors(Sortie& os, const double current_time) const;

  void calculer_var_volume_remaillage(double timestep,
                                      const DoubleTab& vitesses_translation_bulles,
                                      const DoubleTab& mean_bubble_rotation_vector,
                                      const DoubleTab& centre_gravite,
                                      ArrOfDouble& var_volume);
  void calculer_vecteurs_de_deplacement_rigide(DoubleTab& vitesses_translation_bulles,
                                               DoubleTab& mean_bubble_rotation_vector,
                                               DoubleTab& centre_gravite,
                                               const int first_step_interface_smoothing = 0);
  void transporter_maillage_deformation(const int correction_semi_locale_volume_bulle,
                                        const DoubleTab& vitesses_translation_bulles,
                                        const DoubleTab& mean_bubble_rotation_vector,
                                        const DoubleTab& centre_gravite,
                                        const double dt_tot,
                                        ArrOfDouble& dvol,
                                        const int rk_step,
                                        const int first_step_interface_smoothing = 0);
  void transporter_maillage_remaillage(int correction_semi_locale_volume_bulle,
                                       const DoubleTab& vitesses_translation_bulles,
                                       const DoubleTab& mean_bubble_rotation_vector,
                                       const DoubleTab& centre_gravite,
                                       double dt_tot,
                                       ArrOfDouble& dvol,
                                       const int rk_step,
                                       const double temps);
  void transporter_maillage_rigide(const double dt_tot,
                                   const DoubleTab& vitesses_translation_bulles,
                                   const DoubleTab& mean_bubble_rotation_vector,
                                   const DoubleTab& centre_gravite,
                                   const int rk_step,
                                   const int first_step_interface_smoothing = 0);
  void calculer_vitesse_de_deformation(int compo,
                                       const DoubleTab& bounding_box_bulles,
                                       const Cut_field_vector3_double& cut_field_velocity,
                                       const DoubleTab& vitesses_translation_bulles,
                                       const DoubleTab& mean_bubble_rotation_vector,
                                       const DoubleTab& positions_bulles);

  void calculer_bounding_box_bulles(DoubleTab& bounding_box, int option_shear = 0) const;
  void preparer_duplicata_bulles(const DoubleTab& bounding_box_of_bubbles,
                                 const DoubleTab& bounding_box_offsetp,
                                 const DoubleTab& bounding_box_offsetm,
                                 const DoubleTab& authorized_bounding_box,
                                 ArrOfInt& masque_duplicata_pour_compo_reel);

  void preparer_duplicata_bulles_masque_6bit(const DoubleTab& bounding_box,
                                             const DoubleTab& authorized_bounding_box,
                                             ArrOfInt& masque_duplicata_pour_compo);
  void dupliquer_bulle_perio(ArrOfInt& masque_duplicata_pour_compo);
  void creer_duplicata_bulles();
  void supprimer_duplicata_bulles();
  void supprimer_certaines_bulles_reelles();

  void deplacer_bulle_perio(const ArrOfInt& masque_deplacement_par_compo);
  void transferer_bulle_perio(); // Les bulles trop proche du bord sont
  // deplacees a l'oppose, vers l'autre bord perio.

  void update_indicatrice_variables_monofluides();

  void compute_vinterp();
  // methode pour bulles fixes
  void compute_external_forces_(IJK_Field_vector3_double& rappel_ft,
                                IJK_Field_vector3_double& rappel,
                                const IJK_Field_vector3_double& vitesse,
                                const IJK_Field_double& indic_ns,
                                const IJK_Field_double& indic_ft,
                                const double coef_immo,
                                const int tstep,
                                const double current_time,
                                const double coef_ammortissement,
                                const double coef_rayon_force_rappel,
                                double compteur,
                                double coef_mean_force,
                                double coef_force_time_n);
  void compute_external_forces_parser(IJK_Field_vector3_double& rappel,
                                      const IJK_Field_double& indic_ns,
                                      const DoubleTab& individual_forces,
                                      const ArrOfDouble& volume_reel,
                                      const DoubleTab& position,
                                      const double coef_rayon_force_rappel);
  void compute_external_forces_color_function(IJK_Field_vector3_double& rappel,
                                              const IJK_Field_double& indic_ns,
                                              const IJK_Field_double& indic_ft,
                                              DoubleTab& individual_forces,
                                              const ArrOfDouble& volume_reel,
                                              const DoubleTab& position);

  void compute_indicatrice_non_perturbe(IJK_Field_double& indic_np,
                                        const IJK_Field_double& indic,
                                        const ArrOfDouble& volume_reel,
                                        const DoubleTab& position) const;

  int lire_motcle_non_standard(const Motcle& un_mot, Entree& is) override;

  void activate_cut_cell();
  void imprime_bilan_indicatrice();

  void set_fichier_reprise_interface(const Nom& prefix)
  {
    set_reprise(1);
    fichier_reprise_interface_ = prefix + fichier_reprise_interface_;
  }

  void calcul_surface_efficace_face(TYPE_SURFACE_EFFICACE_FACE type_surface_efficace_face, double timestep, const Cut_field_vector3_double& total_velocity);
  void calcul_surface_efficace_interface(TYPE_SURFACE_EFFICACE_INTERFACE type_surface_efficace_interface, double timestep, const Cut_field_vector3_double& velocity);
  void calcul_surface_efficace_face_initial(TYPE_SURFACE_EFFICACE_FACE type_surface_efficace_face);
  void calcul_surface_efficace_interface_initial(TYPE_SURFACE_EFFICACE_INTERFACE type_surface_efficace_interface);

  // fin de methode pour bulles fixes
  const Domaine_dis_base& get_domaine_dis() const { return refdomaine_dis_.valeur(); }
  int get_dt_impression_bilan_indicatrice() const { return dt_impression_bilan_indicatrice_; }
  int get_nb_bulles_reelles() const { return nb_bulles_reelles_; }
  int get_flag_positions_reference() const { return flag_positions_reference_; }
  int is_frozen() const { return frozen_; }
  void freeze() { frozen_ = 1; }

  int get_nb_bulles_ghost(const int print = 0) const
  {
    if (print)
      Cerr << "IJK_Interfaces::get_nb_bulles_ghost : On a " << nb_bulles_ghost_ << " bulles ghost." << finl;
    return nb_bulles_ghost_;
  }

  int get_forcing_method() const { return parser_; }
  int get_recompute_indicator() const { return recompute_indicator_; }
  void set_recompute_indicator(int i) { recompute_indicator_ = i; }
  int follow_colors() const { return follow_colors_; }
  const ArrOfInt& get_colors() const { return through_yminus_; }
  int ghost_compo_converter(const int i) const { return ghost_compo_converter_[i]; }

// Methodes d'acces aux parametres de la repulsion a la paroi :
  double portee_wall_repulsion() const { return portee_wall_repulsion_; }
  double delta_p_wall_max_repulsion() const { return delta_p_wall_max_repulsion_; }

  void set_reprise(const int i)
  {
    reprise_ = i;
    return;
  }

  void set_fichier_sauvegarde(const char *lataname) { fichier_sauvegarde_interface_ = lataname; }
  void set_fichier_reprise(const char *lataname) { fichier_reprise_interface_ = lataname; }
  void set_seuil_indicatrice_petite(double seuil_indicatrice_petite)
  {
    seuil_indicatrice_petite_ = seuil_indicatrice_petite;

    if (seuil_indicatrice_petite_ <= seuil_indicatrice_negligeable_)
      {
        Cerr << "Erreur IJK_Interfaces: Incoherence parametres jeu de donnees, seuil_indicatrice_petite_ <= seuil_indicatrice_negligeable_." << finl;
        Process::exit();
      }
  };

  Nom get_fichier_reprise() { return fichier_reprise_interface_; }

  // methode non-const. Pour appel depuis l'exterieur (comme IJK_switch)
  void lire_maillage_ft_dans_lata()
  {
    maillage_ft_ijk_.lire_maillage_ft_dans_lata(fichier_reprise_interface_,
                                                timestep_reprise_interface_,
                                                lata_interfaces_meshname_);
    return;
  };

  // methodes non-const.
  void parcourir_maillage()
  {
    maillage_ft_ijk_.parcourir_maillage();
    if (!maillage_ft_ijk_.Surfactant_facettes().get_disable_surfactant())
      {
        maillage_ft_ijk_.update_gradient_laplacien_Surfactant();
        maillage_ft_ijk_.update_sigma_and_interfacial_source_term_sommet(ref_domaine_, false, use_tryggvason_interfacial_source_);
      }
    return;
  };

  void RK3_G_store_vi_resize(int n, int n2)
  {
    RK3_G_store_vi_.resize(n, n2);
    return;
  };

  void RK3_G_store_vi_echange_esp_vect()
  {
    maillage_ft_ijk_.desc_sommets().echange_espace_virtuel(RK3_G_store_vi_);
    return;
  };

  // methode non-const. Remplit le tableau des positions de reference
  // vers lesquelles la force de rappel attire les bulles):
  // (a partir de la position actuelle des bulles)
  // pour calculs a bulles fixes
  void set_positions_reference()
  {
    ArrOfDouble vol;
    calculer_volume_bulles(vol, positions_reference_);
    flag_positions_reference_ = 1;
    return;
  };

  const Maillage_FT_IJK& maillage_ft_ijk() const { return maillage_ft_ijk_; }
  const Maillage_FT_IJK& old_maillage_ft_ijk() const { return old_maillage_ft_ijk_; }
  const Remaillage_FT_IJK& remaillage_ft_ijk() const { return remaillage_ft_ijk_; }
  const DoubleTab& RK3_G_store_vi() const { return  RK3_G_store_vi_; }
  const IntVect& get_num_compo() const { return num_compo_; }
  const ArrOfInt& get_compo_to_group() const { return compo_to_group_; }

  // Methode qui parcourt les facettes contenues dans la cellule num_elem
  // pour trouver la phase de sa cellule voisine
  // Donnees d'entree:
  //   - num_elem   : indice de l'element traverse par l'interface
  //   - direction  : precise dans quelle direction se trouve le voisin dont on
  //   cherche la phase
  //   - face_plus  : precise dans quel sens se trouve le voisin dont on cherche
  //   la phase
  //                  (+1 pour le voisin d'indice plus eleve, -1 pour l'autre )
  int compute_cell_phase_with_interface_normal(int num_elem, int direction, int face_plus);

  void calculer_kappa_ft(IJK_Field_double& kappa_ft);

  void calculer_normales_et_aires_interfaciales(IJK_Field_double& ai,
                                                IJK_Field_double& kappa_ai,
                                                IJK_Field_vector3_double& normale_cell,
                                                const int igroup) const;

  int compute_list_compo_connex_in_element(const Maillage_FT_IJK& mesh,
                                           const int elem,
                                           ArrOfInt& liste_composantes_connexes_dans_element) const;
  const int& nb_groups() const { return nb_groups_ ; }
  const ArrOfDoubleFT& get_distance_autres_interfaces() const { return distance_autres_interfaces_; }
  int get_ghost_number_from_compo(const int compo) const;

  void calculer_surface_bulles(ArrOfDouble& surfaces) const;
  void compute_surface_average_per_bubble(const ArrOfDouble& surfaces, const ArrOfDouble& in, ArrOfDouble& out) const;
  void read_bubbles_barycentres_old_new(const Nom& interf_name);
  bool read_bubbles_barycentres_vel(const Nom& interf_name,
                                    FixedVector<ArrOfDouble,3>& bubbles_rising_dir,
                                    FixedVector<ArrOfDouble,3>& bubbles_rising_vel,
                                    ArrOfDouble& bubbles_rising_vel_mag);
  bool read_bubbles_barycentres(const Nom& interf_name, const Nom& suffix, FixedVector<ArrOfDouble,3>& bubbles_bary);
  void store_bubbles_barycentres(const Nom& interf_name);
  void compute_bubbles_volume_and_barycentres(ArrOfDouble& volumes,
                                              DoubleTab& barycentres,
                                              const int& store_values);
  void calculer_volume_bulles(ArrOfDouble& volumes,
                              DoubleTab& centre_gravite) const;

  void calculer_aspect_ratio(ArrOfDouble& aspect_ratio) const;
  void calculer_surfactant(ArrOfDouble& surfactant,ArrOfDouble& surfactant_min,ArrOfDouble& surfactant_max) const;
  void calculer_poussee_bulles(const DoubleTab& gravite, DoubleTab& poussee) const;
  void calculer_aire_interfaciale(IJK_Field_double& ai) const;
  void calculer_aire_interfaciale_for_compo(IJK_Field_double& ai, const int compo) const;
  double calculer_aire_interfaciale_for_compo(const int compo, const int i_ref, const int j_ref, const int k_ref) const;

  void calculer_normale_et_aire_interfaciale(IJK_Field_double& ai,
                                             IJK_Field_double& kappa_ai,
                                             IJK_Field_vector3_double& normale_cell) const;

  void compute_drapeaux_vapeur_v2(const IntVect& vecteur_composantes,
                                  ArrOfInt& drapeau_liquide) const;
  void compute_drapeaux_vapeur_v3(const Maillage_FT_IJK& mesh,
                                  const Domaine_IJK& split,
                                  const IntVect& vecteur_composantes,
                                  ArrOfInt& drapeau_vapeur) const;
  void compute_drapeaux_vapeur_v4(const IntVect& vecteur_composantes,
                                  ArrOfInt& drapeau_vapeur) const;

  void convert_to_IntVect(const ArrOfInt& in, IntVect& out) const;

  void ajouter_terme_source_interfaces(
    IJK_Field_vector3_double& vpoint,
    IJK_Field_vector3_double& vrepul,
    IJK_Field_vector3_double& vabsrepul
  ) const;

  void remailler_interface(const double temps,
                           Maillage_FT_IJK& maillage,
                           ArrOfDouble& var_volume,
                           Remaillage_FT_IJK& algo_remaillage_local);

  int is_terme_gravite_rhog() const;
  void detecter_et_supprimer_rejeton(bool duplicatas_etaient_presents);
  void update_surface_normale() const;

  // Permet de recuperer un mcu qui est une vue de maillage_bulles_ft_ijk. Il
  // n'y a pas de copie memoire, seulement des passages de case memoire.
  // splitting.get_origin, get_constant_delta(DIRECTION_X),
  // ...
  static void get_maillage_MED_from_IJK_FT(MEDCouplingUMesh *maillage_bulles_mcu,
                                           const Maillage_FT_IJK& maillage_bulles_ft_ijk);

  const IJK_Field_double& get_surface_interface_old_ft() const { return surface_interface_ft_[old()]; }
  const IJK_Field_double& get_surface_interface_old() const { return surface_interface_ns_[old()]; }
  const IJK_Field_double& get_surface_interface_next_ft() const { return surface_interface_ft_[next()]; }
  const IJK_Field_double& get_surface_interface_next() const { return surface_interface_ns_[next()]; }

  const IJK_Field_vector3_double& get_barycentre_phase1_old_ft() const { return barycentre_phase1_ft_[old()]; }
  const IJK_Field_vector3_double& get_barycentre_phase1_old() const { return barycentre_phase1_ns_[old()]; }
  const IJK_Field_vector3_double& get_barycentre_phase1_next_ft() const { return barycentre_phase1_ft_[next()]; }
  const IJK_Field_vector3_double& get_barycentre_phase1_next() const { return barycentre_phase1_ns_[next()]; }

  double get_barycentre(bool next_time, int bary_compo, int phase, int i, int j, int k) const;
  double get_barycentre_face(bool next_time, int face_dir, int bary_compo, int phase, int i, int j, int k) const;

  // Getter des surfaces par face
  const IJK_Field_vector3_double& get_surface_vapeur_par_face_ft() const { return surface_vapeur_par_face_[old()]; }
  const IJK_Field_vector3_double& get_surface_vapeur_par_face() const { return surface_vapeur_par_face_ns_[old()]; }
  const IJK_Field_vector3_double& get_indicatrice_surfacique_face_ft() const { return indicatrice_surfacique_face_ft_[old()]; }
  const IJK_Field_vector3_double& get_indicatrice_surfacique_face_old() const { return indicatrice_surfacique_face_ns_[old()]; }
  const IJK_Field_vector3_double& get_indicatrice_surfacique_face_next() const { return indicatrice_surfacique_face_ns_[next()]; }
  const DoubleTabFT_cut_cell_vector3& get_indicatrice_surfacique_efficace_face() const { return indicatrice_surfacique_efficace_face_; }
  const DoubleTabFT_cut_cell_scalar& get_surface_efficace_interface() const { return surface_efficace_interface_; }
  const DoubleTabFT_cut_cell_vector3& get_vitesse_deplacement_interface() const { return vitesse_deplacement_interface_; }
  const DoubleTabFT_cut_cell_vector3& get_normale_deplacement_interface() const { return normale_deplacement_interface_; }
  const DoubleTabFT_cut_cell_vector3& get_coord_deplacement_interface() const { return coord_deplacement_interface_; }

  // Getter des surfaces par face
  // void get_surface_vapeur_par_face_ns(IJK_Field_vector3_double &surfs) const ;
  // Getter des barycentres par face
  const FixedVector<IJK_Field_vector3_double, 3>& get_barycentre_vapeur_par_face_ft() const { return barycentre_vapeur_par_face_[old()]; }
  const FixedVector<IJK_Field_vector3_double, 3>& get_barycentre_vapeur_par_face() const { return barycentre_vapeur_par_face_ns_[old()]; }
  const FixedVector<FixedVector<IJK_Field_double, 2>, 3>& get_barycentre_phase1_face_ft() const { return barycentre_phase1_face_ft_[old()]; }
  const FixedVector<FixedVector<IJK_Field_double, 2>, 3>& get_barycentre_phase1_face() const { return barycentre_phase1_face_ns_[old()]; }
  static inline double opposing_barycentre(double initial_barycentre, double initial_area)
  {
    double weighted_barycentre = initial_barycentre*initial_area;
    double opposing_area = 1 - initial_area;
    double opposing_barycentre = ((opposing_area == 0.) || (opposing_area == 1.)) ? 1./2. : (1./2. - weighted_barycentre)/opposing_area;
    //if (opposing_barycentre == .5 && (opposing_area > 0. && opposing_area < 1.))
    //  {
    //    assert(false);
    //  }
    return opposing_barycentre;
  }

  int get_nb_face_mouillees() const { return n_faces_mouilles_[old()]; }

  // TODO: utiliser le allocate de allocate_velociter dans
  // IJK_Navier_Stockes_tools.cpp utiliser le pslitting de NS, pas le FT

  //////////////////////////////////////////////////////////////////////
  // API indicatrice et variables du maillage par cellule et par face //
  //////////////////////////////////////////////////////////////////////

  static double mean_over_compo(
    const FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& field_for_compo,
    const IJK_Field_int& nb_compo_traversante,
    const int i,
    const int j,
    const int k
  )
  {
    double res = 0.;
    for (int compo = 0; compo <= nb_compo_traversante(i, j, k); compo++)
      {
        res += field_for_compo[compo](i, j, k);
      }
    if (nb_compo_traversante(i,j,k) > 0)
      return res / nb_compo_traversante(i, j, k);
    else
      return 0.;
  }

  static void mean_over_compo(
    const FixedVector<IJK_Field_double, 3 * max_authorized_nb_of_components_>& field_for_compo,
    const IJK_Field_int& nb_compo_traversante,
    IJK_Field_vector3_double& mean_par_compo_field
  )
  {
    const int ni = nb_compo_traversante.ni();
    const int nj = nb_compo_traversante.nj();
    const int nk = nb_compo_traversante.nk();
    for (int dir=0; dir < 3; dir ++)
      for (int i=0; i < ni; i++)
        for (int j=0; j < nj; j++)
          for (int k=0; k < nk; k++)
            {
              double res = 0.;
              for (int compo = 0; compo <= nb_compo_traversante(i, j, k); compo++)
                {
                  int idx = 3 * compo + dir;
                  const double last_val = field_for_compo[idx](i, j, k);
                  res += last_val;
                }
              if (nb_compo_traversante(i,j,k) > 0)
                res = res / nb_compo_traversante(i, j, k);
              else
                res = 0.;
              mean_par_compo_field[dir](i,j,k) = res;
            }
  }

  static double mean_over_compo(
    const FixedVector<IJK_Field_double, 3 * max_authorized_nb_of_components_>& field_for_compo,
    const IJK_Field_int& nb_compo_traversante,
    const int dir,
    const int i,
    const int j,
    const int k
  )
  {
    double res = 0.;
    for (int compo = 0; compo <= nb_compo_traversante(i, j, k); compo++)
      {
        int idx = 3 * compo + dir;
        const double last_val = field_for_compo[idx](i, j, k);
        res += last_val;
      }
    if (nb_compo_traversante(i,j,k) > 0)
      return res / nb_compo_traversante(i, j, k);
    else
      return 0.;
  }

  int old() const { return 1 - old_en_premier_; }
  int next() const { return old_en_premier_; }

  // TODO: retirer l'acces a FT
  const FixedVector<IJK_Field_double, max_authorized_nb_of_components_ * 3>& get_bary_par_compo_itfc_in_cell_ft() const { return bary_par_compo_[old()]; }
  // TODO: retirer l'acces a FT
  const FixedVector<IJK_Field_double, max_authorized_nb_of_components_ * 3>& get_norm_par_compo_itfc_in_cell_ft() const { return normale_par_compo_[old()]; }

  // TODO: retirer l'acces a FT
  const IJK_Field_double& I_ft() const { return indicatrice_ft_[old()]; }
  const double& I_ft(const int i, const int j, const int k) const { return indicatrice_ft_[old()](i, j, k); }
  const IJK_Field_double& In_ft() const { return indicatrice_ft_[next()]; }
  const double& In_ft(const int i, const int j, const int k) const { return indicatrice_ft_[next()](i, j, k); }

  const IJK_Field_vector3_double& BoI() const { return bary_of_interf_ns_[old()]; }
  const IJK_Field_vector3_double& BoIn() const { return bary_of_interf_ns_[next()]; }

  const IJK_Field_double& I() const { return indicatrice_ns_[old()]; }
  const IJK_Field_double& In() const { return indicatrice_ns_[next()]; }
  inline double I(const int i, const int j, const int k) const { return indicatrice_ns_[old()](i, j, k); }
  inline double In(const int i, const int j, const int k) const { return indicatrice_ns_[next()](i, j, k); }

  static int convert_indicatrice_to_phase(double indicatrice)
  {
    assert(est_pure(indicatrice));
    int phase = (int)indicatrice;
    return phase;
  }

  // Indicatrice non-zero : Cette indicatrice est utilisee pour calculer l'energie vol*indicatrice*T d'une cellule coupee en Probleme_FTD_IJK_cut_cell.
  // Pour ne pas perdre d'information, le volume de l'autre temps est utilise si le volume est zero pour le temps consideree.
  inline double I_nonzero(const int phase, const int i, const int j, const int k) const
  {
    double current_indic = (phase == 0) ? 1 - I(i, j, k) : I(i, j, k);
    if (current_indic == 0)
      {
        double other_indic = (phase == 0) ? 1 - In(i, j, k) : In(i, j, k);
        return (other_indic == 0) ? 1. : other_indic;
      }
    else
      {
        return current_indic;
      }
  }
  inline double In_nonzero(const int phase, const int i, const int j, const int k) const
  {
    double current_indic = (phase == 0) ? 1 - In(i, j, k) : In(i, j, k);
    if (current_indic == 0)
      {
        double other_indic = (phase == 0) ? 1 - I(i, j, k) : I(i, j, k);
        return (other_indic == 0) ? 1. : other_indic;
      }
    else
      {
        return current_indic;
      }
  }

  static inline int est_pure(double indicatrice) { return ((indicatrice == 0.) || (indicatrice == 1.)); }
  static inline int devient_pure(double old_indicatrice, double next_indicatrice) { return ((!est_pure(old_indicatrice)) && (est_pure(next_indicatrice))); }
  static inline int devient_diphasique(double old_indicatrice, double next_indicatrice) { return ((est_pure(old_indicatrice)) && (!est_pure(next_indicatrice))); }
  static inline int phase_mourrante(int phase, double old_indicatrice, double next_indicatrice) { return devient_pure(old_indicatrice, next_indicatrice) && (IJK_Interfaces::convert_indicatrice_to_phase(1 - next_indicatrice) == phase); }
  static inline int phase_naissante(int phase, double old_indicatrice, double next_indicatrice) { return devient_diphasique(old_indicatrice, next_indicatrice) && (IJK_Interfaces::convert_indicatrice_to_phase(1 - old_indicatrice) == phase); }
  inline double devient_pure(const int i, const int j, const int k) const { return devient_pure(indicatrice_ns_[old()](i, j, k), indicatrice_ns_[next()](i, j, k)); }
  inline double devient_diphasique(const int i, const int j, const int k) const { return devient_diphasique(indicatrice_ns_[old()](i, j, k), indicatrice_ns_[next()](i, j, k)); }
  inline double phase_mourrante(const int phase, const int i, const int j, const int k) const { return phase_mourrante(phase, indicatrice_ns_[old()](i, j, k), indicatrice_ns_[next()](i, j, k)); }
  inline double phase_naissante(const int phase, const int i, const int j, const int k) const { return phase_naissante(phase, indicatrice_ns_[old()](i, j, k), indicatrice_ns_[next()](i, j, k)); }

  inline int est_reguliere(double old_indicatrice, double next_indicatrice) const { return ((old_indicatrice >= seuil_indicatrice_petite_) && (old_indicatrice <= 1-seuil_indicatrice_petite_) && (next_indicatrice >= seuil_indicatrice_petite_) && (next_indicatrice <= 1-seuil_indicatrice_petite_)); }

  inline int est_desequilibre(double indicatrice) const { return (((indicatrice < seuil_indicatrice_petite_) || (indicatrice > 1-seuil_indicatrice_petite_)) && (!est_pure(indicatrice))); }
  inline int a_desequilibre_final(double old_indicatrice, double next_indicatrice) const { return est_desequilibre(next_indicatrice) && (!devient_diphasique(old_indicatrice, next_indicatrice)); }
  inline int a_desequilibre_initial_uniquement(double old_indicatrice, double next_indicatrice) const { return est_desequilibre(old_indicatrice) && (!a_desequilibre_final(old_indicatrice, next_indicatrice)) && (!devient_pure(old_indicatrice, next_indicatrice)); }

  inline int below_small_threshold(double indicatrice) const { return ((indicatrice < seuil_indicatrice_petite_) && (!est_pure(indicatrice))); }
  inline int next_below_small_threshold(double old_indicatrice, double next_indicatrice) const { return below_small_threshold(next_indicatrice) && (!devient_diphasique(old_indicatrice, next_indicatrice)); }
  inline int only_old_below_small_threhshold(double old_indicatrice, double next_indicatrice) const { return below_small_threshold(old_indicatrice) && (!next_below_small_threshold(old_indicatrice, next_indicatrice)) && (!devient_pure(old_indicatrice, next_indicatrice)); }
  inline int next_below_small_threshold_for_phase(int phase, double old_indicatrice, double next_indicatrice) const { return (phase ==0) ? next_below_small_threshold(1 - old_indicatrice, 1 - next_indicatrice) : next_below_small_threshold(old_indicatrice, next_indicatrice); }
  inline int only_old_below_small_threshold_for_phase(int phase, double old_indicatrice, double next_indicatrice) const { return (phase ==0) ? only_old_below_small_threhshold(1 - old_indicatrice, 1 - next_indicatrice) : only_old_below_small_threhshold(old_indicatrice, next_indicatrice); }

  const double& SI(const int compo, const int i, const int j, const int k) const
  {
    return surface_par_compo_[old()][compo](i, j, k);
  }
  double SI(const int i, const int j, const int k) const
  {
    double res = mean_over_compo(surface_par_compo_[old()], nb_compo_traversante_[old()], i, j, k);
    return res;
  }

  const IJK_Field_vector<double, max_authorized_nb_of_groups_>& groups_indicatrice_ft() const { return groups_indicatrice_ft_[old()]; }
  const IJK_Field_vector<double, max_authorized_nb_of_groups_>& groups_indicatrice_ns() const { return groups_indicatrice_ns_[old()]; }
  const IJK_Field_vector<double, max_authorized_nb_of_groups_>& groups_indicatrice_n_ft() const { return groups_indicatrice_ft_[next()]; }
  const IJK_Field_vector<double, max_authorized_nb_of_groups_>& groups_indicatrice_n_ns() const { return groups_indicatrice_ns_[next()]; }
  // const double& SIn(const int compo, const int i, const int j, const int k)
  // const {
  //   return surface_par_compo_next_[compo](i,j,k);
  // }
  // double SIn(const int i, const int j, const int k) const {
  //   double res = mean_over_compo(surface_par_compo_next_, i,j,k);
  //   return res;
  // }

  // Les composantes (dir) du vecteur normal moyen pour chaque compo comnnexe
  // (c) de bulle dans le maille (compo = max_nb_of_compo_ * c + dir)
  const double& nI(const int compo, const int i, const int j, const int k) const
  {
    return normal_of_interf_ns_[old()][compo](i, j, k);
  }
  Vecteur3 nI(const int i, const int j, const int k) const
  {
    const Vecteur3 res(normal_of_interf_ns_[old()][0](i, j, k),
                       normal_of_interf_ns_[old()][1](i, j, k),
                       normal_of_interf_ns_[old()][2](i, j, k));
    return res;
  }
  const double& nIn(const int compo, const int i, const int j, const int k) const
  {
    return normal_of_interf_ns_[next()][compo](i, j, k);
  }
  Vecteur3 nIn(const int i, const int j, const int k) const
  {
    const Vecteur3 res(normal_of_interf_ns_[next()][0](i, j, k),
                       normal_of_interf_ns_[next()][1](i, j, k),
                       normal_of_interf_ns_[next()][2](i, j, k));
    return res;
  }

  // Les composantes (dir) du barycentre de l'interface pour chaque compo
  // comnnexe (c) de bulle dans le maille (compo = max_nb_of_compo_ * c + dir)
  const double& xI(const int compo, const int i, const int j, const int k) const
  {
    return bary_of_interf_ns_[old()][compo](i, j, k);
  }
  Vecteur3 xIn(const int i, const int j, const int k) const
  {
    const Vecteur3 res(bary_of_interf_ns_[next()][0](i, j, k),
                       bary_of_interf_ns_[next()][1](i, j, k),
                       bary_of_interf_ns_[next()][2](i, j, k));
    return res;
  }
  Vecteur3 xI(const int i, const int j, const int k) const
  {
    const Vecteur3 res(bary_of_interf_ns_[old()][0](i, j, k),
                       bary_of_interf_ns_[old()][1](i, j, k),
                       bary_of_interf_ns_[old()][2](i, j, k));
    return res;
  }

  // Surface de la vapeur sur la face dans la direction compo de la cellule ijk
  const double& Sf(const int compo, const int i, const int j, const int k) const
  {
    return surface_vapeur_par_face_ns_[old()][compo](i, j, k);
  }
  // Surface de la vapeur sur la face dans la direction compo de la cellule ijk
  // au temps n+1
  const double& Sfn(const int compo, const int i, const int j, const int k) const
  {
    return surface_vapeur_par_face_ns_[next()][compo](i, j, k);
  }
  // Surface de la vapeur sur la face dans la direction compo de la cellule ijk
  // moyennee entre n et n+1.
  // TODO: faire un tableau surface_vapeur_par_face_moy et le calculer en
  // faisant la moyenne sur plusieurs petits pas.
  double Sfm(const int compo, const int i, const int j, const int k) const
  {
    return (surface_vapeur_par_face_ns_[old()][compo](i, j, k) + surface_vapeur_par_face_ns_[next()][compo](i, j, k)) * 0.5;
  }

  void update_old_intersections(); // Copie de l'interface et des intersections associees au pas de temps precedent
  void switch_indicatrice_next_old();
  void calculer_indicatrice_next(
    const DoubleTab& gravite,
    const double delta_rho,
    const double sigma,
    const double time,
    const int itstep,
    const bool parcourir = true
  );
  void calculer_indicatrice_intermediaire(
    IJK_Field_double& indicatrice_intermediaire_ft_,
    IJK_Field_double& indicatrice_intermediaire_ns_,
    IJK_Field_vector3_double& indicatrice_surfacique_intermediaire_face_ft_,
    IJK_Field_vector3_double& indicatrice_surfacique_intermediaire_face_ns_,
    const bool parcourir = true
  );
  void calculer_indicatrice_avant_remaillage(const bool parcourir = true)
  {
    calculer_indicatrice_intermediaire(
      indicatrice_avant_remaillage_ft_,
      indicatrice_avant_remaillage_ns_,
      indicatrice_surfacique_avant_remaillage_face_ft_,
      indicatrice_surfacique_avant_remaillage_face_ns_,
      parcourir
    );
  }
  void calculer_indicatrice_apres_remaillage(const bool parcourir = true)
  {
    calculer_indicatrice_intermediaire(
      indicatrice_apres_remaillage_ft_,
      indicatrice_apres_remaillage_ns_,
      indicatrice_surfacique_apres_remaillage_face_ft_,
      indicatrice_surfacique_apres_remaillage_face_ns_,
      parcourir
    );
  }
  void set_compute_surfaces_mouillees() { surface_vapeur_par_face_computation_.set_compute_surfaces_mouillees(); }

  const int& nb_compo_traversantes(const int i, const int j, const int k) const
  {
    return nb_compo_traversante_[next()](i,j,k);
  }

  const Intersection_Interface_ijk_cell& get_intersection_ijk_cell() const { return intersection_ijk_cell_; }
  const Intersection_Interface_ijk_face& get_intersection_ijk_face() const { return intersection_ijk_face_; }
  Intersection_Interface_ijk_cell& get_set_intersection_ijk_cell() { return intersection_ijk_cell_; }
  Intersection_Interface_ijk_face& get_set_intersection_ijk_face() { return intersection_ijk_face_; }
  const IJK_Composantes_Connex& get_ijk_compo_connex() const { return ijk_compo_connex_; }

  void compute_compo_connex_from_bounding_box()
  {
    if (Option_IJK::DISABLE_DIPHASIQUE)
      return;
    ijk_compo_connex_.compute_bounding_box_fill_compo_connex();
  }

  void compute_compo_connex_from_interface()
  {
    if (Option_IJK::DISABLE_DIPHASIQUE)
      return;
    ijk_compo_connex_.compute_compo_connex_from_interface();
  }

  void initialise_ijk_compo_connex_bubbles_params()
  {
    if (Option_IJK::DISABLE_DIPHASIQUE)
      return;
    ijk_compo_connex_.initialise_bubbles_params();
  }

  void allocate_ijk_compo_connex_fields(const Domaine_IJK& splitting, const int& allocate_compo_fields)
  {
    if (Option_IJK::DISABLE_DIPHASIQUE)
      return;
    ijk_compo_connex_.allocate_fields(splitting, allocate_compo_fields);
  }

  void associate_rising_velocities_parameters(const Domaine_IJK& splitting,
                                              const int& compute_rising_velocities,
                                              const int& fill_rising_velocities,
                                              const int& use_bubbles_velocities_from_interface,
                                              const int& use_bubbles_velocities_from_barycentres)
  {
    if (Option_IJK::DISABLE_DIPHASIQUE)
      return;
    ijk_compo_connex_.associate_rising_velocities_parameters(splitting,
                                                             compute_rising_velocities,
                                                             fill_rising_velocities,
                                                             use_bubbles_velocities_from_interface,
                                                             use_bubbles_velocities_from_barycentres);
  }

  void compute_rising_velocities_from_compo()
  {
    if (Option_IJK::DISABLE_DIPHASIQUE)
      return;
    ijk_compo_connex_.compute_rising_velocities();
  }

  const DoubleTab& get_bubble_barycentres_old_new(const int& get_new) const
  {
    if (get_new)
      return bubbles_bary_new_;
    else
      return bubbles_bary_old_;
  }

  const DoubleTab& get_bubble_velocities_from_interface() const { return bubbles_velocities_; }
  const DoubleTab& get_bubble_velocities_from_barycentres() const { return bubbles_velocities_bary_; }
  const DoubleTab& get_bubble_rising_vectors_from_barycentres() const { return bubbles_rising_vectors_bary_; }
  const ArrOfDouble& get_bubbles_velocities_magnitude_from_barycentres() const { return bubbles_velocities_bary_; }

  void reset_flags_and_counters();

protected:
  // Met a jour les valeurs de surface_vapeur_par_face_ et barycentre_vapeur_par_face_
  SurfaceVapeurIJKComputation surface_vapeur_par_face_computation_;

  ComputeValParCompoInCell val_par_compo_in_cell_computation_;

  void verif_indic() ;

  void calculer_phi_repuls_sommet(
    ArrOfDouble& potentiels_sommets,
    ArrOfDouble& repulsions_sommets,
    const DoubleTab& gravite,
    const double delta_rho,
    const double sigma,
    const double time,
    const int itstep
  );

  void calculer_phi_repuls_par_compo(
    FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& surf_par_compo,
    FixedVector<FixedVector<IJK_Field_double, max_authorized_nb_of_components_>, 3>& source_interf_par_compo,
    FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& phi_par_compo,
    FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& repuls_par_compo,
    const DoubleTab& gravite,
    const double delta_rho,
    const double sigma,
    const double time,
    const int itstep
  );

  void calculer_indicatrice(IJK_Field_double& indic);
  void calculer_indicatrice_optim(IJK_Field_double& indic);
  void calculer_indicatrices(IJK_Field_vector3_double& indic);
  void calculer_indicatrices_optim(IJK_Field_vector3_double& indic);

  // Methode qui parcourt tous les elements de indic et met a jour uniquement
  // ceux qui etaient traverses par l'interface a l'iteration precedente et qui
  // ne le sont plus a l'iteration courante ces elements ont ete marques au
  // prealable
  int update_indicatrice(IJK_Field_double& indic);

  void calculer_surface_interface(IJK_Field_double& surf_interface, IJK_Field_double& indic);
  void calculer_barycentre(IJK_Field_vector3_double& baric, IJK_Field_double& indic);
  void calculer_indicatrice_surfacique_face(IJK_Field_vector3_double& indic_surfacique_face, IJK_Field_double& indic, IJK_Field_vector3_double& norme);
  void calculer_indicatrice_surfacique_barycentre_face(IJK_Field_vector3_double& indic_surfacique_face, FixedVector<FixedVector<IJK_Field_double, 2>, 3>& baric_face, IJK_Field_double& indic, IJK_Field_vector3_double& norme);


  // Cette methode appelle la methode statique get_maillage_MED_from_IJK_FT sur
  // ses propres membres. Elle met donc a jour le maillage maillage_bulles_med_.

  // TODO: utiliser le allocate de allocate_velocity dans IJK_Navier_Stokes_tools.cpp utiliser le pslitting de NS, pas le FT

  void calculer_vmoy_translation_composantes_connexes(const Maillage_FT_IJK& maillage,
                                                      const ArrOfDouble& surface_facette,
                                                      const ArrOfDouble& surface_par_bulle,
                                                      const ArrOfInt& compo_connexes_facettes,
                                                      const int nbulles_reelles,
                                                      const int nbulles_ghost,
                                                      const DoubleTab& vitesse_sommets,
                                                      DoubleTab& vitesses_translation_sommets) const;
  void calculer_vmoy_rotation_composantes_connexes(const Maillage_FT_IJK& maillage,
                                                   const ArrOfDouble& surface_facette,
                                                   const ArrOfDouble& surface_par_bulle,
                                                   const ArrOfInt& compo_connexes_facettes,
                                                   const int nbulles_reelles,
                                                   const int nbulles_ghost,
                                                   const DoubleTab& centre_gravite,
                                                   const DoubleTab& vitesse_sommets,
                                                   const DoubleTab& vitesse_translation_sommets,
                                                   DoubleTab& mean_bubble_rotation_vector) const;
  void recursive_calcul_distance_chez_voisin(DoubleTab& vinterp_tmp,
                                             int dir,
                                             const Maillage_FT_IJK& mesh,
                                             DoubleTab& coord_sommets,
                                             ArrOfInt& compo_sommet,
                                             ArrOfDouble& distance,
                                             DoubleTab& v_closer,
                                             double distmax);
  void calculer_distance_autres_compo_connexe2(ArrOfDouble& distance,
                                               DoubleTab& v_closer);
  void calculer_distance_autres_compo_connexe_octree(const DoubleTab& sommets_a_tester,
                                                     const ArrOfInt& compo_connexe_sommets,
                                                     const DoubleTab& vinterp_tmp,
                                                     const Maillage_FT_IJK& mesh,
                                                     ArrOfDouble& distance,
                                                     DoubleTab& v_closer,
                                                     const double distmax);

  void calculer_distance_autres_compo_connexe_ijk(const DoubleTab& sommets_a_tester,
                                                  const ArrOfInt& compo_connexe_sommets,
                                                  const DoubleTab& vinterp_tmp,
                                                  const Maillage_FT_IJK& mesh,
                                                  ArrOfDouble& distance,
                                                  DoubleTab& v_closer,
                                                  const double distmax);

// reference vers le splitting_ft_ pour les interfaces :
  OBS_PTR(Domaine_IJK) ref_domaine_;
  OBS_PTR(Domaine_dis_base) refdomaine_dis_;
  OBS_PTR(Probleme_FTD_IJK_base) ref_ijk_ft_;
  OBS_PTR(Switch_FT_double) ref_ijk_ft_switch_;

  // Ou lire le maillage initial, dans la methode initialize():
  // On peut reprendre un fichier lata ou sauv.lata :

  IntVect num_compo_;
  int nb_compo_in_num_compo_ = 0;
  // variales pour calcul a bulles fixes
  DoubleTab mean_force_;
  DoubleTab force_time_n_;
  DoubleTab positions_reference_;
  int flag_positions_reference_ = 0; // Pas de position de reference imposee

  Nom fichier_reprise_interface_;
  int timestep_reprise_interface_ = 1;
  Nom lata_interfaces_meshname_ = "INTERFACES";

  // Pour ecrire dans le fichier sauv :
  Nom fichier_sauvegarde_interface_;
  int timestep_sauvegarde_interface_ = 1;

  // Activation du suivi des couleurs des bulles
  int follow_colors_ = 0;

  // Activer la repulsion aux parois :
  int active_repulsion_paroi_ = 0;
  int use_tryggvason_interfacial_source_ = 0;
  // La repulsion paroi est desactive par defaut,
  // meme si l'inter-bulles l'est

  // Pour calculer le terme source comme grad(potentiel*I) au lieu de
  // potentiel_face*gradI
  // Modification de l'evaluation du potentiel :
  int correction_gradient_potentiel_ = 0;

  int compute_distance_autres_interfaces_ = 0;
  //  recompute_indicator_ = 1; // doit-on calculer l'indicatrice avec une methode
  //  de debug (1) ou optimisee (0) ?

  // Nombres de bulles reeles :
  int nb_bulles_reelles_ = 0;
  int nb_bulles_ghost_ = 0;
  int nb_bulles_ghost_before_ = 0;
  int recompute_indicator_ = 1;
  int parser_ = 0;
  // doit-on calculer le forcage avec une methode de parser (1,
  // lente) ou optimisee basee sur le num_compo_ (0) ?
  //             Dans le cas ou parser_ est a zero, il faut que
  //             recompute_indicator_ soit a 1, car c'est cette methode qui
  //             rempli num_compo_

  // Stockage du maillage:
  Maillage_FT_IJK maillage_ft_ijk_;
  Maillage_FT_IJK old_maillage_ft_ijk_;

  // Tableau intermediaire pour le deplacement des marqueurs en RK3 :
  DoubleTab RK3_G_store_vi_;
  DoubleTab vinterp_;

  int disable_rigid_translation_ = 0; // Desactive la partie translation du mouvement rigide
  int disable_rigid_rotation_ = 1;    // Desactive la partie rotation du mouvement rigide

  ArrOfDouble var_volume_deformation_;  // Variation de volume observee sur chaque sommet lors de la deformation de la bulle
  ArrOfDouble var_volume_remaillage_;   // Variation de volume cible pour l'operation de remaillage
  ArrOfDouble var_volume_correction_globale_;  // Variation de volume cible pour la correction globale de volume

  IJK_Field_vector3_double deformation_velocity_;

  // Algorithmes de parcours de l'interface (intersections Eulerien/Lagrangien)
  Parcours_interface parcours_;
  // Le parcours a besoin d'une connectivite des faces de bords du maillage
  // eulerien
  Connectivite_frontieres connectivite_frontieres_;
  // Algorithme de remaillage local
  Remaillage_FT_IJK remaillage_ft_ijk_;

  // Donnees sur la maillage NS : sa bounding_box et la periodicite :
  DoubleTab bounding_box_NS_domain_;
  bool perio_NS_[3];

  int avoid_duplicata_ = 0;
  double factor_length_duplicata_ = 1.;

  // Domaine autorise pour les bulles :
  // c'est le geom du splitting_FT reduit de ncells_forbidden_
  // dans toutes les direction ou le domaine NS est perio.
  // ncells_forbidden_ est le nombre de mailles au bord du domaine etendu ou on
  // interdit a des bulles d'entrer (bulles detruites et remplacees par leur
  // duplicata de l'autre cote)
  int ncells_forbidden_ = 3; // Valeur recommandee par defaut.
  // Suppression des bulles sur le pourtour du domaine lors de la sauvegarde
  // finale.
  int ncells_deleted_ = -1; // Valeur recommandee par defaut. On ne veut pas supprimer de bulles.
  int frozen_ = 0; // flag to disable the interfaces motion. By default, we want the motion of the interfaces.
  DoubleTab bounding_box_forbidden_criteria_;
  DoubleTab bounding_box_delete_criteria_;

  // Domaine a l'interieur duquel la duplication d'une bulle est inutile
  DoubleTab bounding_box_duplicate_criteria_;

  // Distance max en metres a laquelle agit la force de repulsion entre les
  // bulles
  double portee_force_repulsion_ = 1.e-8;
  // delta de pression maxi cree par la force de repulsion
  // (pour l'instant lineaire, valeur max quand la distance est nulle)
  double delta_p_max_repulsion_ = 0.; // desactive par defaut

  // Si souhaite, une valeur differente pour les parois :
  double portee_wall_repulsion_ = 1.e-8;
  double delta_p_wall_max_repulsion_ = 0.; // desactive par defaut
  int no_octree_method_ = 0;    // to use the IJK-discretization to search for closest faces of vertices instead of the octree method (disabled by default)

  ArrOfDoubleFT distance_autres_interfaces_;

  // Pour chaque bulle, a-t-elle traverse la frontiere perio y?
  ArrOfInt through_yminus_; // 0 : On ne sait pas
  // 1 : entree par yminus
  // -1 : entree par yplus

  ArrOfInt ghost_compo_converter_;

  int reprise_ = 0; // Flag indiquant si on fait une reprise

  enum Terme_Gravite { GRAVITE_RHO_G, GRAVITE_GRAD_I };
  // Terme_Gravite terme_gravite_;
  int terme_gravite_ = GRAVITE_GRAD_I; // Par defaut terme gravite ft sans courants parasites

  int nb_groups_ = 1;           // Nombre de groupes/classes de bulles. Par defaut toutes les bulles sont dans le meme group.
  ArrOfInt compo_to_group_; // Tableau de conversion: numero_de_groupe =
  // compo_to_group_[icompo]

  /////////////////////////////////////////////////////////////////
  // Tableau des valeurs aux faces mouillees (calcules avec med) //
  /////////////////////////////////////////////////////////////////

  bool compute_surf_mouillees_ = false; // active seulement dans le cas
  // ou il y a des champs thermique ou d energie.
  // attention, ca desactive seulement le calcul, pas l'allocation.
  FixedVector<int,2> n_faces_mouilles_;

  // Surfaces vapeur des faces du maillage IJK
  FixedVector<IJK_Field_vector3_double, 2> surface_vapeur_par_face_;
  FixedVector<IJK_Field_vector3_double, 2> surface_vapeur_par_face_ns_;

  // Normale de l'interface par maille ijk sur domaine NS
  FixedVector<IJK_Field_vector3_double, 2> normal_of_interf_;
  FixedVector<IJK_Field_vector3_double, 2> normal_of_interf_ns_;
  // Barycentre de l'interface par maille ijk sur domaine NS
  FixedVector<IJK_Field_vector3_double, 2> bary_of_interf_;
  FixedVector<IJK_Field_vector3_double, 2> bary_of_interf_ns_;

  // Barycentres vapeur des faces du maillage IJK
  FixedVector<FixedVector<IJK_Field_vector3_double, 3>, 2> barycentre_vapeur_par_face_;
  FixedVector<FixedVector<IJK_Field_vector3_double, 3>, 2> barycentre_vapeur_par_face_ns_;

  /////////////////////////////////////
  // indicatrice et var moy par cell //
  /////////////////////////////////////

  int n_cell_diph_ = 0;
  bool old_en_premier_ = true;

  FixedVector<IJK_Field_double, 2> indicatrice_ns_;
  FixedVector<IJK_Field_double, 2> indicatrice_ft_;

  // Indicatrice apres la deformation de l'interface mais avant le remaillage/lissage de l'interface
  IJK_Field_double indicatrice_avant_remaillage_ns_;
  IJK_Field_double indicatrice_avant_remaillage_ft_;

  // Indicatrice apres le remaillage/lissage de l'interface mais avant le deplacement rigide
  IJK_Field_double indicatrice_apres_remaillage_ns_;
  IJK_Field_double indicatrice_apres_remaillage_ft_;

  // Indicatrice cible apres le remaillage. Le remaillage vise a atteindre cette indicatrice.
  IJK_Field_double delta_volume_theorique_bilan_ns_;

  FixedVector<IJK_Field_double, 2> surface_interface_ns_;
  FixedVector<IJK_Field_double, 2> surface_interface_ft_;

  FixedVector<IJK_Field_vector3_double, 2> barycentre_phase1_ns_;
  FixedVector<IJK_Field_vector3_double, 2> barycentre_phase1_ft_;


  IJK_Field_double field_repulsion_;

  // Indicatrice surfacique aux faces du maillage cartesien,
  // indiquant la fraction de la surface de chaque face associee a la phase.
  // Note : similaire a surface_vapeur_par_face_, mais sans medcoupling.
  FixedVector<IJK_Field_vector3_double, 2> indicatrice_surfacique_face_ns_;
  FixedVector<IJK_Field_vector3_double, 2> indicatrice_surfacique_face_ft_;

  // Indicatrice surfacique apres la deformation de l'interface mais avant le remaillage/lissage de l'interface
  IJK_Field_vector3_double indicatrice_surfacique_avant_remaillage_face_ns_;
  IJK_Field_vector3_double indicatrice_surfacique_avant_remaillage_face_ft_;

  // Indicatrice surfacique apres le remaillage/lissage de l'interface mais avant le deplacement rigide (pas utilisee)
  IJK_Field_vector3_double indicatrice_surfacique_apres_remaillage_face_ns_;
  IJK_Field_vector3_double indicatrice_surfacique_apres_remaillage_face_ft_;

  FixedVector<FixedVector<FixedVector<IJK_Field_double, 2>, 3>, 2> barycentre_phase1_face_ns_;
  FixedVector<FixedVector<FixedVector<IJK_Field_double, 2>, 3>, 2> barycentre_phase1_face_ft_;

  // On prevoie un tableau assez grand pour contenir tous les groupes.
  FixedVector<IJK_Field_vector<double, max_authorized_nb_of_groups_>, 2> groups_indicatrice_ft_;
  FixedVector<IJK_Field_vector<double, max_authorized_nb_of_groups_>, 2> groups_indicatrice_ns_;

#if VERIF_INDIC
  // pour verifier le calcul optimise de l'indicatrice
  IJK_Field_double indicatrice_ft_test_;
  IJK_Field_vector<double, max_authorized_nb_of_groups_> groups_indicatrice_ft_test_;
#endif

  // Vecteur des composantes normale dans chaque cellule par composante connexe
  FixedVector<IJK_Field_int, 2> nb_compo_traversante_;
  FixedVector<FixedVector<IJK_Field_int, max_authorized_nb_of_components_>, 2> compos_traversantes_;
  FixedVector<FixedVector<IJK_Field_double, max_authorized_nb_of_components_>, 2> indicatrice_par_compo_;
  FixedVector<FixedVector<IJK_Field_double, max_authorized_nb_of_components_>, 2> surface_par_compo_;
  FixedVector<FixedVector<IJK_Field_double, max_authorized_nb_of_components_>, 2> courbure_par_compo_;
  FixedVector<FixedVector<IJK_Field_double, max_authorized_nb_of_components_>, 2> phi_par_compo_;
  FixedVector<FixedVector<IJK_Field_double, max_authorized_nb_of_components_>, 2> surf_par_compo_;
  FixedVector<FixedVector<FixedVector<IJK_Field_double, max_authorized_nb_of_components_>, 3>, 2> source_interf_par_compo_;
  FixedVector<FixedVector<IJK_Field_double, max_authorized_nb_of_components_>, 2> gradx_sigma_par_compo_;
  FixedVector<FixedVector<IJK_Field_double, max_authorized_nb_of_components_>, 2> grady_sigma_par_compo_;
  FixedVector<FixedVector<IJK_Field_double, max_authorized_nb_of_components_>, 2> gradz_sigma_par_compo_;
  FixedVector<FixedVector<IJK_Field_double, max_authorized_nb_of_components_>, 2> repuls_par_compo_;
  FixedVector<FixedVector<IJK_Field_double, 3 * max_authorized_nb_of_components_>, 2> normale_par_compo_;
  FixedVector<FixedVector<IJK_Field_double, 3 * max_authorized_nb_of_components_>, 2> bary_par_compo_;

  Intersection_Interface_ijk_cell intersection_ijk_cell_;
  Intersection_Interface_ijk_face intersection_ijk_face_;

  IJK_Composantes_Connex ijk_compo_connex_;
  DoubleTab bubbles_velocities_;
  DoubleTab bubbles_velocities_bary_;
  DoubleTab bubbles_rising_vectors_bary_;
  ArrOfDouble bubbles_velocities_bary_magnitude_;
  DoubleTab bubbles_bary_old_;
  DoubleTab bubbles_bary_new_;
  int read_barycentres_velocity_ = 0;
  int use_barycentres_velocity_ = 0;
  bool has_computed_bubble_barycentres_ = false;
  bool has_readen_barycentres_prev_ = false;

  int dt_impression_bilan_indicatrice_ = -1;
  int verbosite_surface_efficace_face_ = 1;
  int verbosite_surface_efficace_interface_ = 1;
  double seuil_indicatrice_petite_ = -1;
  double seuil_indicatrice_negligeable_ = 1e-6;

  // Pour le calcul des champs cut-cell
  int cut_cell_activated_ = 0;
  DoubleTabFT_cut_cell_vector3 indicatrice_surfacique_efficace_face_;
  DoubleTabFT_cut_cell_vector3 indicatrice_surfacique_efficace_face_initial_;
  DoubleTabFT_cut_cell_vector6 indicatrice_surfacique_efficace_face_correction_;
  DoubleTabFT_cut_cell_scalar indicatrice_surfacique_efficace_face_absolute_error_;
  IJK_Field_vector3_double indicatrice_surfacique_efficace_deformation_face_; // ne peut etre sur la structure diphasique, car cree et utilisee avant
  DoubleTabFT_cut_cell_scalar surface_efficace_interface_;
  DoubleTabFT_cut_cell_scalar surface_efficace_interface_initial_;
  DoubleTabFT_cut_cell_vector3 coord_deplacement_interface_;
  DoubleTabFT_cut_cell_vector3 vitesse_deplacement_interface_;
  DoubleTabFT_cut_cell_vector3 normale_deplacement_interface_;

  std::map<Motcle, IJK_Field_double> scalar_post_fields_;

private:

  Champs_compris_IJK champs_compris_;

};

#endif /* IJK_Interfaces_included */
