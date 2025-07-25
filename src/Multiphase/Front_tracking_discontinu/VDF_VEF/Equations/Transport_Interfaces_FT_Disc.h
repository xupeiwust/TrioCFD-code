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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Transport_Interfaces_FT_Disc.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/38
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Transport_Interfaces_FT_Disc_included
#define Transport_Interfaces_FT_Disc_included

#include <Equation_base.h>
#include <Transport_Interfaces_base.h>
#include <Postraitement_base.h>

#include <Remaillage_FT.h>
#include <Parcours_interface.h>
#include <Marching_Cubes.h>
#include <Connectivite_frontieres.h>
#include <Topologie_Maillage_FT.h>
#include <Algorithmes_Transport_FT_Disc.h>

#include <Navier_Stokes_FT_Disc.h>
#include <Proprietes_part_vol.h>
#include <TRUSTTabs_forward.h>
#include <TRUSTTabFT_forward.h>
#include <TRUST_Ref.h>

#include <Collision_Model_FT_base.h>
#include <Post_Processing_Hydrodynamic_Forces.h>
#include <Post_Processing_Hydrodynamic_Forces_Stokes.h>

#include <map>

class Probleme_base;
class Milieu_base;
class Navier_Stokes_FT_Disc;
class Loi_horaire;

class Transport_Interfaces_FT_Disc_interne;

class Transport_Interfaces_FT_Disc : public Transport_Interfaces_base
{
  Declare_instanciable_sans_constructeur(Transport_Interfaces_FT_Disc);
public:

  Transport_Interfaces_FT_Disc();
  friend class Post_Processing_Hydrodynamic_Forces;
  //
  void set_param(Param& titi) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  // Methodes virtuelles pures de Equation_base
  //
  int            nombre_d_operateurs() const override; // Zero, y'en a pas.
  const Operateur& operateur(int i) const override;    // Erreur
  Operateur&        operateur(int i) override;         // Erreur
  const Champ_Inc_base& inconnue() const override;         // C'est l'indicatrice
  Champ_Inc_base&        inconnue() override;

  // Methodes surchargees de Equation_base
  //
  void                associer_milieu_base(const Milieu_base& milieu) override;

  void                 associer_equation_ns(const Navier_Stokes_FT_Disc& ns);

  Milieu_base&        milieu() override;       // Erreur
  const Milieu_base& milieu() const override;  // Erreur
  void    associer_pb_base(const Probleme_base& probleme) override;
  void    discretiser() override;
  Entree& lire_cond_init(Entree& is) override;
  int  verif_Cl() const override;
  double  calculer_pas_de_temps() const override;
  DoubleTab& derivee_en_temps_inco(DoubleTab& derivee) override;
  void assembler( Matrice_Morse& mat_morse, const DoubleTab& present, DoubleTab& secmem) override ;

  void    mettre_a_jour(double temps) override;
  std::vector<YAML_data> data_a_sauvegarder() const override;
  int  sauvegarder(Sortie& ) const override;
  int  reprendre(Entree&) override;
  int impr(Sortie& os) const override;
  void update_critere_statio();

  // To save front in a different file
  void init_save_file() override;
  void close_save_file() override;

  //
  // Nouvelles methodes
  //
  virtual void                            lire_maillage_ft_cao(Entree& is);
  int                          preparer_calcul() override;
  virtual void                            preparer_pas_de_temps();
  const Maillage_FT_Disc&                 maillage_interface() const;
  const Champ_base&               get_update_indicatrice() override;
  virtual const Champ_base&               get_indicatrice_faces();
  virtual const Champ_base&               get_compute_indicatrice_faces();
  virtual const Parcours_interface&       parcours_interface() const;
  virtual const Marching_Cubes&           marching_cubes() const;
  virtual const Algorithmes_Transport_FT_Disc& algorithmes_transport() const;
  virtual const Connectivite_frontieres& connectivite_frontieres() const;
  Remaillage_FT&                          remaillage_interface();
  const Remaillage_FT&                    remaillage_interface() const;
  const Topologie_Maillage_FT&            topologie_interface() const;
  virtual double calculer_integrale_indicatrice(const DoubleVect& indicatrice, double& v_ph0) const;

  const Proprietes_part_vol&           proprietes_particules() const;
  const Maillage_FT_Disc&              maillage_inject() const;
  const Proprietes_part_vol&           proprietes_inject() const;

  void nettoyer_proprietes_particules(const ArrOfInt& som_utilises);

  virtual void calculer_vitesse_transport_interpolee(const Champ_base& champ_vitesse,
                                                     const Maillage_FT_Disc& m,
                                                     DoubleTab& vitesse_noeuds,
                                                     int nv_calc) const
  {
    calculer_vitesse_transport_interpolee(champ_vitesse, m, vitesse_noeuds, nv_calc, 1);
  };

  virtual void calculer_vitesse_transport_interpolee(const Champ_base& champ_vitesse,
                                                     const Maillage_FT_Disc&,
                                                     DoubleTab& vitesse_noeuds,
                                                     int nv_calc,
                                                     int standard) const;
  void calculer_scalaire_interpole(const Champ_base& ch_scal,
                                   const Maillage_FT_Disc&,
                                   DoubleTab& ch_scal_noeuds,
                                   int nv_calc) const;

  virtual void remailler_interface();

  //methodes utilisees pour le post-traitement
  virtual int get_champ_post_FT(const Motcle& champ, Postraitement_base::Localisation loc, DoubleTab *dtab = 0) const;
  virtual int get_champ_post_FT(const Motcle& champ, Postraitement_base::Localisation loc, IntTab    *itab = 0) const;
  virtual const Maillage_FT_Disc& maillage_interface_pour_post() const;
  virtual const int& get_n_iterations_distance() const;
  int get_mesh_tag() const override
  {
    return maillage_interface_pour_post().get_mesh_tag();
  };


  //Methode d acces au probleme
  const Probleme_base& get_probleme_base() const;

  //Modifie vpoint (pour N_S) pour imposer au fluide la vitesse de l interface
  void modifier_vpoint_pour_imposer_vit(const DoubleTab& inco_val,DoubleTab& vpoint0,
                                        DoubleTab& vpoint,const DoubleTab& rho_faces,
                                        DoubleTab& terme_source,const double temps, const double dt,
                                        const int is_explicite,const double eta) override;

  //Methode outil utilisee par modifier_vpoint_pour_imposer_vit(...)
  void calcul_source(const DoubleTab& inco_val,
                     const DoubleTab& vpoint,
                     const DoubleTab& rho_faces,
                     DoubleTab& source_val,
                     const DoubleTab& vit_imposee,
                     const DoubleTab& indicatrice_faces,
                     const int is_QC,
                     const double dt,
                     const int is_explicite,
                     const double eta);
  void modifie_source(DoubleTab& so_modif,const DoubleTab& so_val,const DoubleTab& rho_faces,
                      const int n,const int m, const int is_QC,
                      const DoubleVect& vol_entrelaces,const Solveur_Masse_base& solv_masse);

  void calcul_effort_fluide_interface(const DoubleTab& vpoint,const DoubleTab& rho_faces,
                                      DoubleTab& source_val,const int is_explicite,const double eta);

  void impr_effort_fluide_interface( DoubleTab& source_val, DoubleTab& pressure_part, DoubleTab& friction_part, DoubleTab& diff_part ) ;

  //Calcul la vitesse imposee a l interface a partir de expression_vitesse_imposee
  virtual void calcul_vitesse(DoubleTab& vitesse_imp, const DoubleTab& champ_vitesse,
                              const DoubleTab& vpoint, const double temps, const double dt);
  virtual void get_expression_vitesse_imposee(DoubleTab& vitesse_imp);
  //Effectue l integration d un ensemble de points (sans notion de facettes)
  void integrer_ensemble_lagrange(const double temps) override;

  virtual void interpoler_vitesse_face(const DoubleTab& distance_interface,
                                       const int phase, const int stencil_width,
                                       DoubleTab& champ, DoubleTab& gradient,
                                       const double t, const double dt ) ;

  void interpoler_simple_vitesse_face(const DoubleTab& distance_interface,
                                      const int phase, const int stencil_width,
                                      DoubleTab& champ, DoubleTab& gradient,
                                      const double t, const double dt ) ;

  virtual void calcul_nb_traverse(  const DoubleTab& xe, const double dx,
                                    const int dim, const int ori,
                                    Maillage_FT_Disc& maillage, int elem,
                                    int& traverse ) ;
  virtual void calcul_OutElemFa7( Maillage_FT_Disc& maillage,
                                  const DoubleTab& indicatrice,
                                  const int nb_elem,
                                  int& nb_fa7_accepted,
                                  IntList& OutElem, IntList& OutFa7 ) ;
  virtual void PPP_face_interface( Maillage_FT_Disc& maillage, const DoubleTab& indicatrice,
                                   const DoubleTab& indicatrice_face, DoubleTab& Vertex ) ;

  virtual void PPP_face_interface_voisin( const DoubleTab& indicatrice, const DoubleTab& indicatrice_face,
                                          DoubleTab& Vertex, DoubleTab& PPP ) ;
  virtual void PPP_face_voisin( const DoubleTab& indicatrice, const DoubleTab& indicatrice_face, DoubleTab& PPP ) ;

  virtual void calcul_maxfa7( Maillage_FT_Disc& maillage, const DoubleTab& indicatrice,
                              const int nb_elem, int& max_fa7, const int exec_planfa7existan ) ;
  virtual void RenumFa7( DoubleTab& Barycentre, DoubleTab& Tab110,DoubleTab& Tab111,
                         DoubleTab& Tab112, IntTab& Tab12, IntTab& CptFacette,
                         const int nb_facettes, const int nb_facettes_dim ) ;
  virtual void StockageFa7( Maillage_FT_Disc& maillage, IntTab& CptFacette, DoubleTab& Tab100,
                            DoubleTab& Tab101,DoubleTab& Tab102, DoubleTab& Tab103,
                            DoubleTab& Tab110,DoubleTab& Tab111,DoubleTab& Tab112,
                            IntTab& Tab12, DoubleTab& Barycentre, const DoubleTab& indicatrice,
                            IntList& OutElem, ArrOfBit& fa7, const int exec_planfa7existant) ;
  virtual void StockageFa7( Maillage_FT_Disc& maillage, DoubleTab& Tab100, DoubleTab& Tab101,
                            DoubleTab& Tab102, DoubleTab& Tab103, DoubleTab& Tab110,
                            DoubleTab& Tab111, DoubleTab& Tab112, IntTab& Tab12,
                            DoubleTab& Barycentre, IntList& OutElem, IntTab& TabOutFa7, ArrOfBit& fa7 ) ;
  virtual void BaryFa7( Maillage_FT_Disc& maillage, const int i_facette, DoubleTab& Barycentre ) ;
  virtual void plan_facette_existant( Maillage_FT_Disc& maillage,
                                      DoubleList A, DoubleList B, DoubleList C,
                                      DoubleList D, const int i_facette,
                                      int& test_liste ) ;
  virtual void calcul_eq_plan_facette(Maillage_FT_Disc& maillage, const int i_facette,
                                      double& a, double& b, double& c, double& d);
  virtual void calcul_eq_plan_facette(const int i_facette,
                                      double& a, double& b, double& c, double& d);

  virtual void calcul_tolerance_projete_monophasique( const int i_face, const int ori, const int voisin0,
                                                      const int voisin1, const DoubleTab& indicatrice_face,
                                                      const DoubleTab& indicatrice, double& tol ) ;

  virtual void calcul_tolerance_projete_diphasique( const int i_face, const int ori, const int voisin0,
                                                    const int voisin1, const DoubleTab& indicatrice, double& tol ) ;

  void verifprojete(const int monophasique,const double Lref, double d, const DoubleTab& x,
                    const DoubleTab& V, DoubleTab& coord_projete, int& cpt ) ;


  virtual void uzawa(const double d,const DoubleTab& matrice, const DoubleTab& x,
                     const DoubleTab& secmem, DoubleTab& solution) const ;

  virtual void projete_point_face_fluide( int& nb_proj_modif, const int dim_fa7,
                                          const DoubleTab& indicatrice_face, const DoubleTab& indicatrice,
                                          const DoubleTab& dist_face, const double t, const double dt,
                                          DoubleTab& Tab100, DoubleTab& Tab101,DoubleTab& Tab102,
                                          DoubleTab& Tab103, IntTab& Tab12, IntTab& CptFacette,
                                          DoubleTab& v_imp, DoubleTab& Vertex,
                                          Parser& parser_x, Parser& parser_y,Parser& parser_z );
  virtual void projete_point_face_interface(   int& nb_proj_modif,const int dim_fa7,
                                               const DoubleTab& indicatrice_face,
                                               const DoubleTab& indicatrice,
                                               const DoubleTab& dist_face, const double t,
                                               const double dt, DoubleTab& Tab100,
                                               DoubleTab& Tab101,DoubleTab& Tab102,
                                               DoubleTab& Tab103, IntTab& Tab12,
                                               IntTab& CptFacette, DoubleTab& v_imp, DoubleTab& Vertex,
                                               Parser& parser_x, Parser& parser_y,Parser& parser_z) ;

  virtual void transporter_sans_changement_topologie(DoubleTab& vitesse,
                                                     const double coeff,const double temps);

  virtual int calculer_composantes_connexes_pour_suppression(IntVect& num_compo);
  virtual double suppression_interfaces(const IntVect& num_compo, const ArrOfInt& flags_compo_a_supprimer);

  const int& get_vimp_regul() const;

  virtual const Champ_base& get_update_distance_interface() const;
  virtual const Champ_base& get_update_distance_interface_faces() const;
  virtual const Champ_base& get_update_normale_interface() const;
  // renvoie DoubleTab parce qu'il n'existe pas de champ aux sommets en VDF ! Zut...
  virtual const DoubleTab&   get_update_distance_interface_sommets() const;
  void ramasse_miettes(const Maillage_FT_Disc& maillage,
                       DoubleVect& flux,
                       DoubleVect& valeurs);
  void nettoyer_maillage()
  {
    maillage_interface().nettoyer_maillage();
  };

  // On utilise des OWN_PTR() pour ne pas avoir a inclure la definition
  // de ces classes (pour reduire les dependances).
  static void transfert_conservatif_eulerien_vers_lagrangien_sommets(const Maillage_FT_Disc& maillage,
                                                                     const DoubleVect& valeurs_euler,
                                                                     ArrOfDouble& valeurs_lagrange);

  // for fpi module
  // getters/setters
  Collision_Model_FT_base& get_set_collision_model() { return collision_model_.valeur(); }

  // getters
  const int& get_nb_particles_tot() const { return nb_particles_tot_; }
  const Collision_Model_FT_base& get_collision_model() const { return collision_model_.valeur(); }
  const OWN_PTR(Collision_Model_FT_base)& get_ptr_collision_model() const { return collision_model_; }
  const DoubleTab& get_particles_position() const { return particles_position_collision_; }
  const DoubleTab& get_particles_velocity() const { return particles_velocity_collision_; }
  const ArrOfInt& get_gravity_center_elem() const { return gravity_center_elem_; }
  Post_Processing_Hydrodynamic_Forces& get_post_process_hydro_forces()
  { return post_process_hydro_forces_; }
  void associate_temp_equation_post_processing(OBS_PTR(Convection_Diffusion_Temperature_FT_Disc) ref_eq_temp)
  { post_process_hydro_forces_.associate_temp_equation(ref_eq_temp); }
  const bool& get_is_solid_particle() const { return is_solid_particle_; }

protected:

  virtual void calculer_vmoy_composantes_connexes(const Maillage_FT_Disc& maillage,
                                                  const ArrOfInt& compo_connexes_facettes,
                                                  const int nb_compo_tot,
                                                  const DoubleTab& vitesse_sommets,
                                                  DoubleTab& vitesses,
                                                  DoubleTab& positions) const;

  void ajouter_contribution_saut_vitesse(DoubleTab& deplacement) const;
  virtual void deplacer_maillage_ft_v_fluide(const double temps);

  virtual void calculer_distance_interface(const Maillage_FT_Disc& maillage,
                                           DoubleTab& distance_elements,
                                           DoubleTab& normale_elements,
                                           const int n_iter) const;

  virtual void calculer_distance_interface_sommets(const DoubleTab& dist_elem,
                                                   const DoubleTab& normale_elem,
                                                   DoubleTab&        dist_som) const;


  virtual void calculer_vitesse_repere_local(const Maillage_FT_Disc& maillage,
                                             DoubleTab& deplacement,
                                             DoubleTab& Positions,
                                             DoubleTab& Vitesses) const;
  virtual void test_suppression_interfaces_sous_domaine();


  virtual void calculer_distance_interface_faces(const DoubleTab& dist_elem,
                                                 const DoubleTab& normale_elem,
                                                 DoubleTab&        dist_faces) const;
  // Nouvelle methodes
  // Methode outil utilisee par modifier_vpoint_pour_imposer_vit(...)
  // Calcul l'indicatrice sur chaque face
  void calcul_indicatrice_faces(const DoubleTab& indicatrice,
                                const IntTab& face_voisins);

  // for fpi module
  // getters
  const DoubleTab& get_mean_particles_volumic_velocity() const
  { return mean_particles_volumic_velocity_; }
  const DoubleTab& get_mean_particles_volumic_squared_velocity() const
  { return mean_particles_volumic_squared_velocity_; }
  const DoubleTab& get_rms_particles_volumic_velocity() const
  { return rms_particles_volumic_velocity_; }
  const DoubleTab& get_particles_purely_solid_mesh_volume() const
  { return particles_purely_solid_mesh_volume_; }
  //setters
  void set_is_solid_particle(const bool is_solid_particle) { is_solid_particle_=is_solid_particle; }
  void init_particles_position_velocity();
  void swap_particles_lagrangian_position_velocity();
  void compute_particles_rms();

  void add_fields_to_post_FT(Motcles& fields) const;
  void fill_ftab_vertices_curvature(DoubleTab *ftab,const DoubleTab& dummytab) const;
  void fill_ftab_velocity(DoubleTab *ftab,const DoubleTab& dummytab) const;
  void fill_ftab_local_reference_frame_velocity(DoubleTab *ftab,const DoubleTab& dummytab) const;
  void fill_ftab_normal_unit(DoubleTab *ftab,const DoubleTab& dummytab) const;
  void fill_ftab_pressure(DoubleTab *ftab,const DoubleTab& dummytab) const;
  void fill_ftab_pressure_force(DoubleTab *ftab,const DoubleTab& dummytab) const;
  void fill_ftab_friction_force(DoubleTab *ftab,const DoubleTab& dummytab) const;
  void fill_ftab_Stokes_pressure_interp(DoubleTab* ftab, const DoubleTab& values) const;
  void fill_ftab_Stokes_pressure_th(DoubleTab* ftab, const DoubleTab& values) const;
  void fill_ftab_Stokes(DoubleTab* ftab, const DoubleTab& values) const;

  OBS_PTR(Probleme_base) probleme_base_;
  OBS_PTR(Navier_Stokes_FT_Disc) equation_ns_;
  // L'inconnue du probleme
  OWN_PTR(Champ_Inc_base) indicatrice_;
  OWN_PTR(Champ_Inc_base) indicatrice_faces_;
  // Utiliser ces accesseurs :
  Maillage_FT_Disc& maillage_interface();
  // Utiliser ces accesseurs :
  Marching_Cubes& marching_cubes();
  // Utiliser ces accesseurs :
  Topologie_Maillage_FT& topologie_interface();

  Proprietes_part_vol&         proprietes_particules();
  Maillage_FT_Disc&            maillage_inject();
  Proprietes_part_vol&         proprietes_inject();

  DoubleTab& tableaux_positions();
  IntVect& vecteur_elements();
  DoubleTab&    deplacement_som();

  Nom suppression_interfaces_sous_domaine_;

  OWN_PTR(Champ_Fonc_base)  vitesse_imp_interp_;

  // for fpi module
  bool is_solid_particle_=false; // pointer to NS_FT_Disc::is_solid_particle_
  int compute_particles_rms_=0;
  OWN_PTR(Collision_Model_FT_base) collision_model_;
  mutable DoubleTab particles_position_collision_; // for contact forces computation
  mutable DoubleTab particles_velocity_collision_; // for contact forces computation
  mutable ArrOfInt gravity_center_elem_; // number of the element which contains the gravity center of the particles
  DoubleTab mean_particles_volumic_velocity_;
  DoubleTab mean_particles_volumic_squared_velocity_;
  DoubleTab rms_particles_volumic_velocity_;
  DoubleTab particles_purely_solid_mesh_volume_;

  mutable Post_Processing_Hydrodynamic_Forces post_process_hydro_forces_;
  mutable Post_Processing_Hydrodynamic_Forces_Stokes post_process_hydro_forces_Stokes_;




private:
  // Variables internes a la methode de transport
  Transport_Interfaces_FT_Disc_interne *variables_internes_;

  double temps_debut_;

  OBS_PTR(Milieu_base) ref_milieu_;

  int interpolation_repere_local_;
  ArrOfDouble force_;
  ArrOfDouble moment_;

  void compute_nb_particles_tot();
  // for map_element_post_FT ...
  void fill_ftab_scalar(DoubleTab *ftab, const ArrOfDouble& values) const;
  void fill_ftab_scalar(DoubleTab *ftab, const DoubleVect& values) const;
  void fill_ftab_scalar(DoubleTab *ftab, const DoubleTab& values) const;
  void fill_ftab_vector(DoubleTab *ftab, const DoubleTab& values) const;

  int nb_particles_tot_=0;

  struct map_element_post_FT
  {
    using func_type=  void (Transport_Interfaces_FT_Disc::*)(DoubleTab*,const DoubleTab&) const;
    map_element_post_FT() {};
    map_element_post_FT(const Motcle& location, func_type function, DoubleTab* ptr, const DoubleTab&  values):
      location_(location),
      function_(function),
      ptr_(ptr),
      values_(values)
    {};
    Motcle location_;
    func_type function_;
    DoubleTab* ptr_;
    DoubleTab  values_;
  };

  using my_map=std::map<Motcle, map_element_post_FT>;

  void fill_map_post_FT(my_map& map_post, DoubleTab *ftab) const;

};

class Transport_Interfaces_FT_Disc_interne : public Objet_U
{
  Declare_instanciable_sans_constructeur_ni_destructeur(Transport_Interfaces_FT_Disc_interne);
public:
  Transport_Interfaces_FT_Disc_interne() :
    indicatrice_cache_tag(-1),
    iterations_correction_volume(0),
    VOFlike_correction_volume(0),
    nb_lissage_correction_volume(0),
    nb_iterations_correction_volume(3),
    volume_impose_phase_1(-1.),
    n_iterations_distance(3),
    n_iterations_interpolation_ibc(5),
    distance_normale_cache_tag(-1),
    distance_sommets_cache_tag(-1),
    distance_faces_cache_tag(-1),
    methode_transport(INDEFINI),
    methode_interpolation_v(VALEUR_A_ELEM),
    injection_interfaces_last_time_(0.),
    interpolation_champ_face(BASE),
    vf_explicite(0),
    is_extra_diphasique_solide(0),
    is_extra_solide(0),
    is_distance_projete_face(0),
    nomb_fa7_accepted(3),
    seuil_uzawa(1.e-8),
    nb_iter_uzawa(30),
    vimp_regul(1),
    type_indic_faces_(STANDARD),
    modified_indic_faces_position(0.),
    modified_indic_faces_thickness(1.),
    type_vitesse_imposee(UNIFORME),
    type_distance_calculee(DIST_INITIALE),
    type_projete_calcule(PROJETE_INITIAL),
    expression_vitesse_imposee(Objet_U::dimension)


  {};
  ~Transport_Interfaces_FT_Disc_interne() override
  {};
private:
  Transport_Interfaces_FT_Disc_interne(const Transport_Interfaces_FT_Disc_interne& a): Objet_U(a)
  {}; // Interdit
  const Transport_Interfaces_FT_Disc_interne& operator=(const Transport_Interfaces_FT_Disc_interne& a)
  {
    return *this;
  }; // Interdit

public:
  std::vector<YAML_data> data_a_sauvegarder() const;
  int sauvegarder(Sortie& os) const override;
  int reprendre(Entree& is) override;
  void init_save_file();
  void close_save_file();
  Nom maillage_interface_xyz_filename(int restart) const;

  // Les membres suivantes sont sauvegardes et repris:
  OWN_PTR(Champ_Inc_base)        indicatrice_cache;     // L'indicatrice calculee par get_update_indicatrice
  int           indicatrice_cache_tag; // Le tag du maillage correspondant
  Maillage_FT_Disc maillage_interface;          // Objet qui peut se reduire a un ensemble de sommets
  // quand il represente les positions de particules
  Remaillage_FT    remaillage_interface_;
  // Fin des membres sauvegardes / repris

  // For now, the front can not be saved with the PDI format
  // (as the mesh changes with every time step, which means the arrays do not have constant dimensions during the simulation
  // and this configuration is not supported yet with the current implementation)
  // so the front will be saved in a different file
  OWN_PTR(Sortie_Fichier_base) fic_front_sauv_;

  Nom pb_name_;
  inline void set_pb_name(const Nom& name) { pb_name_ = name; }
  Nom restart_fname_;
  inline void set_restart_fname(const Nom& name) { restart_fname_ = name; }
  Nom checkpoint_fname_;
  inline void set_checkpoint_fname(const Nom& name) { checkpoint_fname_ = name; }
  Nom ident_;
  inline void set_ident(const Nom& name) { ident_ = name; }
  double time_;
  inline void set_time(double t) { time_ = t; }

  Proprietes_part_vol proprietes_particules_; //Proprietes physiques de particules

  Maillage_FT_Disc maillage_inject_;              //Ensemble de particules a injecter periodiquement
  Proprietes_part_vol proprietes_inject_;     //Proprietes physiques des particules injectees

  // FORMER Keyword (obsolete):
  // Si iterations_correction_volume == 0, le maillage est transporte par le fluide
  // a l'aide du champ de vitesse L2 interpole (pas de conservation du volume).
  // Si iterations_correction_volume > 0, on calcule une correction de volume aux
  // sommets de l'interface et on l'etale par autant d'iterations d'un lisseur
  // Voir Transport_Interfaces_FT_Disc::mettre_a_jour
  int iterations_correction_volume;

  // NEW Keywords/parameters for volume preserving correction in agreement with phase change :
  int VOFlike_correction_volume;
  // Si VOFlike_correction_volume == 0, le maillage est transporte par le fluide
  // a l'aide du champ de vitesse L2 interpole (pas de conservation du volume).
  // Si VOFlike_correction_volume > 0, on calcule une correction de volume aux
  // sommets de l'interface et on l'etale par nb_lissage_correction_volume iterations d'un lisseur,
  // Voir Transport_Interfaces_FT_Disc::mettre_a_jour
  // pour eviter l'apparition de pic aux interfaces (ie, on lisse legement la correction de volume
  // (uniquement s'il y en a une))
  int nb_lissage_correction_volume;
  // La correction est iterative car on ne corrige pas exactement du volume demande en deplacant les
  // noeuds sequentiellement. On fait nb_iterations_correction_volume ou jusqu'a ce que l'erreur
  // soit inferieure au seuil de correction de volume de Remaillage_FT
  int nb_iterations_correction_volume;

  // ADDITIONAL GLOBAL mass conservation with-out phase change:
  // Rustine introduite pour corriger les pertes de masse:
  // Si cette valeur est positive, on deplace toute l'interface d'une certaine distance
  // pour que le volume de la phase 1 reste toujours egal a la valeur prescrite
  double volume_impose_phase_1;
  // Si non nul, le calcul de l'integrale pour le volume impose porte sur cette sous-domaine
  Nom    nom_domaine_volume_impose_;

  OWN_PTR(Champ_Inc_base) vitesse_filtree;
  DoubleTab doubletab_pos;
  DoubleTab doubletab_vitesses;
  IntVect   intvect_elements;

  OWN_PTR(Champ_Fonc_base)  distance_interface; // Distance a l'interface (aux elements)
  OWN_PTR(Champ_Fonc_base)  normale_interface;  // Une normale etalee
  OWN_PTR(Champ_Fonc_base)  surface_interface;  // GB : La surface d'interface dans chaque element
  OWN_PTR(Champ_Fonc_base)  tmp_flux;           // Tableau temporaire pour le ramasse-miettes
  OWN_PTR(Champ_Fonc_base)  distance_interface_faces; // CF : Distance a l'interface (aux faces)
  DoubleTab  distance_interface_sommets; // Distance a l'interface (aux sommets)
  OWN_PTR(Champ_Fonc_base)  distance_interface_faces_corrigee; // CI : Distance a l'interface corrigee (aux faces)
  OWN_PTR(Champ_Fonc_base)  distance_interface_faces_difference; // CI : Distance a l'interface corrigee - Distance a l'interface calculee (aux faces)
  OWN_PTR(Champ_Fonc_base)  index_element; // CI : indexation des elements
  OWN_PTR(Champ_Fonc_base)  nelem_par_direction; // CI : nombre d'elements par direction
  // Note de B.M. : zut, y'a pas de champ aux sommets en VDF... donc DoubleTab et
  // du coup on ne peut pas postraiter facilement. C'est trop con...
  int     n_iterations_distance;
  int     n_iterations_interpolation_ibc;
  int     distance_normale_cache_tag;
  int     distance_sommets_cache_tag;
  int     distance_faces_cache_tag;
  // Le maillage postraite n'est pas forcement le maillage
  Maillage_FT_Disc maillage_pour_post;
  DoubleTabFT   deplacement_sommets;

  enum Methode_transport { INDEFINI, VITESSE_IMPOSEE, LOI_HORAIRE, VITESSE_INTERPOLEE };
  Methode_transport      methode_transport;
  OBS_PTR(Navier_Stokes_std) refequation_vitesse_transport;

  enum Methode_interpolation_v { VALEUR_A_ELEM, VDF_LINEAIRE, MEAN_VOLUMIC_VELOCITY };
  Methode_interpolation_v methode_interpolation_v;


  // Injecteur d'interfaces au cours du temps
  ArrOfDouble injection_interfaces_temps_;
  ArrOfInt    injection_interfaces_phase_;
  Noms        injection_interfaces_expressions_;
  // temps physique reel auquel a eu lieu la derniere injection (ce n'est pas le temps demande)
  double      injection_interfaces_last_time_;

  enum Interpolation_champ_face { BASE, LINEAIRE };
  Interpolation_champ_face interpolation_champ_face;

  int vf_explicite, is_extra_diphasique_solide, is_extra_solide, is_distance_projete_face;
  int nomb_fa7_accepted ;
  double seuil_uzawa ;
  int nb_iter_uzawa ;
  int vimp_regul ;

  enum Type_indic_faces { STANDARD, MODIFIEE, AI_BASED };
  Type_indic_faces type_indic_faces_;
  double modified_indic_faces_position ;
  double modified_indic_faces_thickness ;

  enum Type_vitesse_imposee { UNIFORME, ANALYTIQUE };
  Type_vitesse_imposee type_vitesse_imposee;

  enum Type_distance_calculee { DIST_INITIALE, DIST_MODIFIEE };
  Type_distance_calculee type_distance_calculee;

  enum Type_projete_calcule { PROJETE_INITIAL, PROJETE_MODIFIE, PROJETE_SIMPLIFIE };
  Type_projete_calcule type_projete_calcule;

  Noms                      expression_vitesse_imposee;
  // Reference a une loi horaire eventuelle
  OBS_PTR(Loi_horaire)         loi_horaire_;

  // Integration de la vitesse a partir du point de depart (x,y,z)
  //  pendant un temps dt.
  void integrer_vitesse_imposee(double temps, double dt, double& x, double& y, double& z) const;

  // Les objets-algorithmes :
  Parcours_interface      parcours_interface_;
  Marching_Cubes          marching_cubes_;
  Connectivite_frontieres connectivite_frontieres_;
  Topologie_Maillage_FT   topologie_interface_;
  // Cet objet est type en fonction de la discretisation:
  OWN_PTR(Algorithmes_Transport_FT_Disc) algorithmes_transport_;

};
#endif
