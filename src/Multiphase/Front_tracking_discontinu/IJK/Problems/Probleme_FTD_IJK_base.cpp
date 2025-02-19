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

#include <Navier_Stokes_FTD_IJK_tools.h>
#include <Schema_Euler_explicite_IJK.h>
#include <IJK_Navier_Stokes_tools.h>
#include <EcritureLectureSpecial.h>
#include <Probleme_FTD_IJK_tools.h>
#include <Probleme_FTD_IJK_base.h>
#include <Navier_Stokes_FTD_IJK.h>
#include <IJK_discretisation.h>
#include <LecFicDiffuse_JDD.h>
#include <IJK_Field_vector.h>
#include <init_forcage_THI.h>
#include <IJK_Lata_writer.h>
#include <Interprete_bloc.h>
#include <MaillerParallel.h>
#include <corrections_qdm.h>
#include <communications.h>
#include <Ouvrir_fichier.h>
#include <Schema_RK3_IJK.h>
#include <Init_spectral.h>
#include <Probleme_base.h>
#include <Domaine_IJK.h>
#include <Option_IJK.h>
#include <Domaine_VF.h>
#include <Type_info.h>
#include <Sonde_IJK.h>
#include <EFichier.h>
#include <SFichier.h>
#include <Force_sp.h>
#include <Force_ph.h>
#include <Vecteur3.h>
#include <EChaine.h>
#include <SChaine.h>
#include <Param.h>

Implemente_base(Probleme_FTD_IJK_base, "Probleme_FTD_IJK_base", Probleme_FT_Disc_gen);
// XD Probleme_FTD_IJK_base interprete Probleme_FTD_IJK_base 1 not_set

Sortie& Probleme_FTD_IJK_base::printOn(Sortie& os) const
{
  return os << " " << que_suis_je() << finl;
}

Entree& Probleme_FTD_IJK_base::readOn(Entree& is)
{
  if (Objet_U::dimension != 3)
    Process::exit("Probleme_FTD_IJK_base currently available in 3D !\n");

  Cerr << "Reading of the problem " << le_nom() << finl;

  Motcle motlu;
  is >> motlu;
  if (motlu != "{")
    {
      Cerr << "We expected { to start to read the problem" << finl;
      Process::exit();
    }

  /* 1 : solved_equations + milieu : NEW SYNTAX */
  lire_solved_equations(is);
  Cerr << "Probleme_FTD_IJK_base::readOn => We expect to read " << (int)equations_.size() << " equations."  << finl;

  typer_lire_milieu(is);

  /* 2 : On lit les equations */
  for (auto& itr : equations_)
    {
      Cerr << "Probleme_FTD_IJK_base::readOn => Will read equation of type : " << itr->que_suis_je() << " ..." << finl;
      is >> motlu;
      is >> getset_equation_by_name(motlu);
      itr->associer_milieu_base(milieu_ijk());
    }

  /* 3 : Tous les params encore non ranges proprement - a changer plus tard */
  Param param(que_suis_je());
  nom_sauvegarde_ = nom_du_cas() + ".sauv";
  set_param(param);
  param.lire_avec_accolades(is);

  /* 4 : Les postraitements */
  is >> motlu;  // Read next word
  // Si le postraitement comprend le mot, on en lit un autre...
  while (les_postraitements_.lire_postraitements(is, motlu, *this))
    is >> motlu;

  completer();

  if (motlu != "}")
    {
      Cerr << "We expected } to start to read the problem" << finl;
      Process::exit();
    }
  return is;
}

void Probleme_FTD_IJK_base::set_param(Param& param)
{
  param.ajouter("nom_sauvegarde", &nom_sauvegarde_); // XD_ADD_P chaine Definition of filename to save the calculation
  param.ajouter_flag("sauvegarder_xyz", &sauvegarder_xyz_); // XD_ADD_P rien save in xyz format
  param.ajouter("nom_reprise", &nom_reprise_); // XD_ADD_P chaine Enable restart from filename given
}

int Probleme_FTD_IJK_base::associer_(Objet_U& obj)
{
  if (sub_type(Domaine_IJK, obj))
    {
      domaine_ijk_ = ref_cast(Domaine_IJK, obj);
      return 1;
    }
  else
    return Probleme_FT_Disc_gen::associer_(obj);
}

void Probleme_FTD_IJK_base::lire_solved_equations(Entree& is)
{
  /* Step 1 : special FT : read the list of equations to solve ... */
  // Here are all possible equations !!!
  Noms noms_eq, noms_eq_maj;
  Type_info::les_sous_types(Nom("Equation_base"), noms_eq);
  for (auto &itr : noms_eq) noms_eq_maj.add(Motcle(itr)); //ha ha ha

  Motcle read_mc;
  Nom nom_eq;
  is >> read_mc;

  if (read_mc != "SOLVED_EQUATIONS")
    {
      Cerr << "Error in Probleme_FTD_IJK_base::lire_solved_equations !!! We expected reading the SOLVED_EQUATIONS bloc instead of " << read_mc << " !!!" << finl;
      Cerr << "Fix your data file !!!" << finl;
      Process::exit();
    }

  is >> read_mc;
  if (read_mc != "{")
    {
      Cerr << "Error in Probleme_FTD_IJK_base::lire_solved_equations !!! We expected { instead of " << read_mc << " !!!" << finl;
      Cerr << "Fix your data file !!!" << finl;
      Process::exit();
    }

  std::vector<Nom> eq_types, eq_name;

  for (is >> read_mc; read_mc != "}"; is >> read_mc)
    {
      const bool non_accepted_eqs = (!read_mc.contient("_FTD_IJK") && (read_mc != "IJK_INTERFACES") && (read_mc != "IJK_THERMALS") );

      if (noms_eq_maj.rang(read_mc) == -1 || non_accepted_eqs)
        {
          Cerr << "Error in Probleme_FTD_IJK_base::lire_solved_equations !!! The equation " << read_mc << " could not be used with a problem of type " << que_suis_je() << " !!!" << finl;
          Cerr << "You can only use the following equations :" << finl;
          for (auto &itr : noms_eq_maj)
            if (itr.contient("_FTD_IJK") || itr == "IJK_INTERFACES" || itr == "IJK_THERMALS")
              Cerr << "  - " << itr << finl;
          Process::exit();
        }

      eq_types.push_back(read_mc);
      is >> nom_eq;
      eq_name.push_back(nom_eq);
    }

  if (eq_types.size() != eq_name.size())
    {
      Cerr << "Error in Probleme_FTD_IJK_base::lire_solved_equations !!! The number of strings read in the bloc SOLVED_EQUATIONS is not correct !!!" << finl;
      Cerr << "Fix your data file !!!" << finl;
      Process::exit();
    }

  /* Step 2 : add equations to the list ... */
  /* Add Navier_Stokes_FTD_IJK at first */
  for (int i = 0; i < static_cast<int>(eq_types.size()); i++)
    if (eq_types[i] == "NAVIER_STOKES_FTD_IJK")
      {
        has_ns_ = true;
        add_FT_equation(eq_name[i], eq_types[i]);
      }

  /* Add Transport_Interfaces at second */
  for (int i = 0; i < static_cast<int>(eq_types.size()); i++)
    if (eq_types[i] == "IJK_INTERFACES")
      {
        has_interface_ = true;
        add_FT_equation(eq_name[i], eq_types[i]);
      }

  /* Add thermals */
  for (int i = 0; i < static_cast<int>(eq_types.size()); i++)
    if (eq_types[i] == "IJK_THERMALS")
      {
        has_thermals_ = true;
        add_FT_equation(eq_name[i], eq_types[i]);
      }
}

void Probleme_FTD_IJK_base::typer_lire_milieu(Entree& is)
{
  const int nb_milieu = 1;
  le_milieu_.resize(nb_milieu);

  for (int i = 0; i < nb_milieu; i++)
    is >> le_milieu_[i]; // On commence par la lecture du milieu
}

void Probleme_FTD_IJK_base::completer()
{
  // Preparer le fichier de postraitement
  lata_name_ = nom_du_cas();
  if (fichier_post_ != "??")
    lata_name_ = fichier_post_;

  lata_name_ += Nom(".lata");

  // FIXME
  Navier_Stokes_FTD_IJK& ns = ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur());
  if (has_interface_)
    ns.associer_interfaces(get_interface());
//  ns.associer_domaine_ft(domaine_ft_);

  milieu_ijk().calculate_direction_gravite(); // pour rotation

  build_vdf_domaine();

  for (auto& itr : equations_)
    itr->completer();

  les_postraitements_.completer();

  get_post().associer_probleme(*this);

  if (nom_reprise_ != "??")
    reprendre_probleme(nom_reprise_);

  schema_temps_ijk().completer();

  // Register domains for post processing object
  get_post().associer_domaines(domaine_ijk_.valeur(), domaine_ft_);
}

const IJK_Field_double& Probleme_FTD_IJK_base::get_IJK_field(const Nom& nom) const
{
// TODO : FIXME : faut boucler plus tard sur les equations IJK
  const Navier_Stokes_FTD_IJK& ns = ref_cast(Navier_Stokes_FTD_IJK, equations_[0].valeur());

  if (ns.has_IJK_field(nom))
    return ns.get_IJK_field(nom);

  if (nom== "INDICATRICE")
    return get_interface().I_ft();


  if (has_thermals_ && get_ijk_thermals().has_IJK_field(nom))
    return get_ijk_thermals().get_IJK_field(nom);

  return get_post().get_IJK_field(nom);
}

void Probleme_FTD_IJK_base::sauvegarder_probleme(const char *fichier_sauvegarde, const int& stop)  //  const
{
  statistiques().begin_count(sauvegarde_counter_);

  const double current_time = schema_temps_ijk().get_current_time();

  // TODO : FIXME : faut boucler plus tard sur les equations IJK
  const Navier_Stokes_FTD_IJK& ns = ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur());
  const auto& velocity = ns.get_velocity();

  Nom lata_name(fichier_sauvegarde);
  Nom interf_name = lata_name + ".interfaces";
  lata_name += ".lata";
  dumplata_header(lata_name, velocity[0] /* on passe un champ pour ecrire la geometrie */);
  dumplata_newtime(lata_name, current_time);
  dumplata_vector(lata_name, "VELOCITY", velocity[0], velocity[1], velocity[2], 0);

  get_post().sauvegarder_post(lata_name);

  if (sauvegarder_xyz_)
    {
      Nom xyz_name(fichier_sauvegarde);
      xyz_name += ".xyz";
      // Nom xyz_name_ascii = xyz_name + "_ascii";
      dumpxyz_vector(*this, velocity, xyz_name, true);
      // dumpxyz_vector(velocity_, xyz_name_ascii, false);
    }
  if (!Option_IJK::DISABLE_DIPHASIQUE)
    get_interface().sauvegarder_interfaces(lata_name, interf_name);

  if (has_thermals_)
    get_ijk_thermals().sauvegarder_temperature(lata_name, stop);

  // curseur = thermique_; // RAZ : Remise au depart du curseur. GB -> Anida : Ne marche pas sur une liste vide? Je dois grader le curseur_bis ensuite.
  SFichier fichier;
  if (Process::je_suis_maitre())
    {
      fichier.ouvrir(fichier_sauvegarde);
      Cerr << "T= " << current_time << " Checkpointing dans le fichier " << fichier_sauvegarde << finl;
      fichier.precision(17);
      fichier.setf(std::ios_base::scientific);

//      Param param(que_suis_je());
//      param.ajouter("tinit", &current_time_);
//      param.ajouter("terme_acceleration_init", &terme_source_acceleration_);
//      // param.ajouter("force_init", &force_init_);
//      param.ajouter("fichier_reprise_vitesse", &fichier_reprise_vitesse_);
//      param.ajouter("timestep_reprise_vitesse", &timestep_reprise_vitesse_);
//      //  param.ajouter("timestep_reprise_rho", &timestep_reprise_rho_);
//      param.ajouter("interfaces", & interfaces_);
//      param.ajouter("statistiques", &statistiques_);
//      param.ajouter("statistiques_FT", &statistiques_FT_);
//      param.ajouter("groups_statistiques_FT", &groups_statistiques_FT_);
//      // S'il y a plusieurs groups, on s'occupe des objets stats pour chaque group:
//      // (en ecrivant directement le vecteur d'objets)
//      param.print(fichier);
      fichier << "{\n";
      if (schema_temps_ijk().get_use_tstep_init())
        fichier << " tstep_init " << (schema_temps_ijk().get_tstep() + schema_temps_ijk().get_tstep_init() + 1) << "\n";

      ns.sauvegarder_equation(lata_name, fichier);

      /*
       * Thermals problems
       */
      //thermals_.sauvegarder_thermals(fichier);

      get_post().sauvegarder_post_maitre(lata_name, fichier);
      fichier << "}\n" ;
      Cerr << "T= " << current_time << " Checkpointing dans le fichier l.1168 " << fichier_sauvegarde << finl;
    }
  statistiques().end_count(sauvegarde_counter_);

}

void Probleme_FTD_IJK_base::reprendre_probleme(const char *fichier_reprise)
{
  // Lecture par tous les processeurs, on retire les commentaires etc...
  LecFicDiffuse_JDD fichier(fichier_reprise);

  // TODO : FIXME : faut boucler plus tard sur les equations IJK
  Navier_Stokes_FTD_IJK& ns = ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur());
  IJK_Interfaces& interf = get_interface();

  Param param(que_suis_je());
  ns.set_param_reprise_pb(param);
  schema_temps_ijk().set_param_reprise_pb(param);
  interf.set_param_reprise_pb(param);

  get_post().reprendre_post(param);

  param.lire_avec_accolades(fichier);

  // Appeler ensuite initialize() pour lire les fichiers lata etc...
  Cerr << "Reprise des donnees a t=" << schema_temps_ijk().get_current_time() << "\n" << finl;

  IJK_Shear_Periodic_helpler::shear_x_time_ = ns.get_boundary_conditions().get_dU_perio()*
                                              (schema_temps_ijk().get_current_time() + ns.get_boundary_conditions().get_t0_shear());

  reprise_ = 1;

  Nom prefix = dirname(fichier_reprise);

  interf.set_fichier_reprise_interface(prefix);

  if (has_thermals_)
    get_ijk_thermals().set_fichier_reprise(prefix + get_ijk_thermals().get_fichier_reprise());

  ns.set_fichier_reprise_vitesse(prefix);
}

// C'est ici aussi qu'on alloue les champs de temperature.
void Probleme_FTD_IJK_base::initialise_ijk_fields()
{
  Cerr << que_suis_je() << "::initialise_ijk_fields()" << finl;

  get_post().initialise(reprise_);
// TODO : FIXME : faut boucler plus tard sur les equations IJK
  Navier_Stokes_FTD_IJK& eq_ns = ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur());
  IJK_Interfaces& interf = get_interface();
  eq_ns.initialise_ijk_fields();

  // L'indicatrice non-perturbee est remplie (si besoin, cad si post-traitement) par le post.complete()
  get_post().fill_indic(reprise_);

  /*
   * Thermal problems
   */
  interf.initialise_ijk_compo_connex_bubbles_params();

  int ghost_fluid_flag = 0;
  if (has_thermals_)
    {
      IJK_Thermals& thermals = get_ijk_thermals();
      thermals.initialize(domaine_ijk_.valeur());
      thermals.get_rising_velocities_parameters(eq_ns.get_compute_rising_velocities(), eq_ns.get_fill_rising_velocities(), eq_ns.get_use_bubbles_velocities_from_interface(),
                                                eq_ns.get_use_bubbles_velocities_from_barycentres());

      ghost_fluid_flag = thermals.ghost_fluid_flag();
    }


  interf.allocate_ijk_compo_connex_fields(domaine_ijk_.valeur(), ghost_fluid_flag || eq_ns.get_upstream_velocity_measured());
  interf.associate_rising_velocities_parameters(domaine_ijk_.valeur(), eq_ns.get_compute_rising_velocities() || eq_ns.get_upstream_velocity_measured(),
                                                eq_ns.get_fill_rising_velocities(), eq_ns.get_use_bubbles_velocities_from_interface(), eq_ns.get_use_bubbles_velocities_from_barycentres());


  eq_ns.complete_initialise_ijk_fields();
}

// Deplacement des interfaces par le champ de vitesse :
// 1. Calcule vitesse_ft (etendue) a partir du champ de vitesse.
// 2. Supprime les duplicatas.
// 3. Transporte le maillage avec velocity_ft_.
// 4. Transfere les bulles a travers la periodicite si besoin.
// 5. Cree les nouveaux duplicatas.
//
// ATTENTION : rho_mu_indicatrice ne sont pas mis a jours.
//
// Mettre rk_step = -1 si schema temps different de rk3.
void Probleme_FTD_IJK_base::deplacer_interfaces(const double timestep, const int rk_step, ArrOfDouble& var_volume_par_bulle, const int first_step_interface_smoothing)
{
  if (Option_IJK::DISABLE_DIPHASIQUE || get_interface().is_frozen())
    return;

  static Stat_Counter_Id deplacement_interf_counter_ = statistiques().new_counter(1, "Deplacement de l'interface");
  statistiques().begin_count(deplacement_interf_counter_);

  // FIXME
  Navier_Stokes_FTD_IJK& ns = ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur());

  //  Calculer vitesse_ft (etendue) a partir du champ de vitesse.
  ns.calculer_vitesse_ft();

  if (has_thermals_)
    {
      IJK_Thermals& thermals = get_ijk_thermals();
      /*
       * Calculation of intersections on interface at time (n)
       */
      Cerr << "Compute Eulerian distance and curvature fields" << finl;
      thermals.compute_eulerian_distance_curvature();
      Cerr << "Clean IJK intersections" << finl;
      thermals.clean_ijk_intersections();
      Cerr << "Copy interface state for post-processing on surface" << finl;
      thermals.copy_previous_interface_state();
      // thermals_.update_intersections(); // no need as IJK_intersections call interfaces_nI interfaces_xI
    }

  get_interface().update_indicatrice_variables_monofluides();

  ns.deplacer_interfaces(timestep, rk_step, var_volume_par_bulle, first_step_interface_smoothing);

  // On supprime les fragments de bulles.
  //interfaces_.detecter_et_supprimer_rejeton(false);

  // On met a jour l'indicatrice du pas de temps d'apres.
  // On met aussi a jour le surf et bary des faces mouillees,
  // les valeurs moyennes en ijk, les val moy en ijkf, etc

  update_indicator_field();

  // mise a jour de l'indicatrice pour les variables monofluides
  ns.update_indicatrice_variables_monofluides();

  statistiques().end_count(deplacement_interf_counter_);
}

// Nouvelle version ou le transport se fait avec les ghost...
void Probleme_FTD_IJK_base::deplacer_interfaces_rk3(const double timestep, const int rk_step,
                                                    ArrOfDouble& var_volume_par_bulle)
{
  if (Option_IJK::DISABLE_DIPHASIQUE || get_interface().is_frozen())
    return;

  // FIXME
  ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur()).deplacer_interfaces_rk3(timestep, rk_step, var_volume_par_bulle);
}

//  Parcourir_maillage cree des noeuds et facettes virtuelles.
//  Pour maintenir un tableau a jour entre avant et apres,
//  il suffit de resizer le tableau a la sortie de la methode
//  (forcement plus grand qu'avant) et de faire un echange_espace_virtuel.
void Probleme_FTD_IJK_base::parcourir_maillage()
{
  IJK_Interfaces& interf = get_interface();
  //const int nbsom_before = interfaces_.maillage_ft_ijk().nb_sommets();
  interf.parcourir_maillage();

  const int nbsom = interf.maillage_ft_ijk().nb_sommets();
  // const int size_store = interfaces_.RK3_G_store_vi().dimension(0);
  // if (!((nbsom >= nbsom_before) &&
  //       ((nbsom_before == size_store) || (0 == size_store) )))
  //   {
  //     Cerr << "Une des tailles de tableau n'est pas bonne... "
  //          << " size_store = " << size_store
  //          << " nbsom_before = " << nbsom_before
  //          << " nbsom = " << nbsom
  //          << finl;
  //     Process::exit();
  //   }
  interf.RK3_G_store_vi_resize(nbsom, 3);
  interf.RK3_G_store_vi_echange_esp_vect();
}

void Probleme_FTD_IJK_base::update_indicator_field()
{
  IJK_Interfaces& interf = get_interface();
  const double delta_rho = milieu_ijk().get_delta_rho();
  interf.switch_indicatrice_next_old();
  interf.calculer_indicatrice_next(get_post().potentiel(),
                                   milieu_ijk().gravite().valeurs(),
                                   delta_rho,
                                   milieu_ijk().sigma(),
                                   schema_temps_ijk().get_current_time(), schema_temps_ijk().get_tstep()
                                  );
}

void Probleme_FTD_IJK_base::update_pre_remeshing_indicator_field()
{
  get_interface().calculer_indicatrice_avant_remaillage();
}

void Probleme_FTD_IJK_base::update_post_remeshing_indicator_field()
{
  get_interface().calculer_indicatrice_apres_remaillage();
}

void Probleme_FTD_IJK_base::update_twice_indicator_field()
{
  for(int i=0; i<2; i++)
    {
      update_indicator_field();
      update_old_intersections();
    }
}

void Probleme_FTD_IJK_base::update_old_intersections()
{
  get_interface().update_old_intersections();
}

void Probleme_FTD_IJK_base::initialise_interfaces()
{
  Cerr << que_suis_je() << "::initialise_interfaces()" << finl;

  const Domaine_dis_base& domaine_dis_ft = refprobleme_ft_disc_->domaine_dis();

  get_interface().initialize(domaine_ft_, domaine_ijk_.valeur(), domaine_dis_ft, thermal_probes_ghost_cells_);

  // On la met a jour 2 fois, une fois next et une fois old
  if (!Option_IJK::DISABLE_DIPHASIQUE)
    update_twice_indicator_field();
}

/*! Nothing to be done here for now ...
 */
void Probleme_FTD_IJK_base::discretiser(Discretisation_base& dis)
{
  Cerr << "Discretization (IJK) of the domain associated with the problem '" << le_nom() << "' (nothing to do ...)" << finl;
  if (!sub_type(IJK_discretisation, dis))
    Process::exit("Error!! IJK problem must be associated with an IJK discretisation!!");
}

void Probleme_FTD_IJK_base::initialize()
{
  Cerr << "Probleme_FTD_IJK_base::initialize()" << finl;

  // FIXME
  Navier_Stokes_FTD_IJK& ns = ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur());
  ns.initialise_ns_fields();

  if (que_suis_je() == "Probleme_FTD_IJK_cut_cell") // TODO should do virtual method here ?
    treatment_count_.allocate(domaine_ijk_.valeur(), Domaine_IJK::ELEM, 2);

  initialise_interfaces();
  initialise_ijk_fields();

  Cerr << " Allocating " << IJK_Field_double::alloc_counter() << " IJK_FT_double objects." << finl;

// Les champs ont etes alloues.
// On peut completer les sondes car les ijk_field.get_domaine() sont a present remplis.
  get_post().completer_sondes();
  get_post().improved_initial_pressure_guess(ns.get_improved_initial_pressure_guess());
}

/*
 * TODO: Change this block with OWN_PTR CLASS IJK_Thermal
 */
void Probleme_FTD_IJK_base::update_thermal_properties()
{
  if (has_thermals_)
    get_ijk_thermals().update_thermal_properties();
}

void Probleme_FTD_IJK_base::preparer_calcul()
{
  // FIXME
  IJK_Interfaces& interf = get_interface();
  Navier_Stokes_FTD_IJK& ns = ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur());
  ns.preparer_calcul();

  if ((!Option_IJK::DISABLE_DIPHASIQUE) && (get_post().get_liste_post_instantanes().contient_("VI") || get_post().get_liste_post_instantanes().contient_("TOUS")))
    interf.compute_vinterp();

  get_post().postraiter_ci(lata_name_, schema_temps_ijk().get_current_time());

  get_post().compute_extended_pressures(interf.maillage_ft_ijk());

  if (has_thermals_)
    schema_temps_ijk().set_modified_time_ini(  get_ijk_thermals().get_modified_time() );

  if (!reprise_ && schema_temps_ijk().get_current_time() == 0.)
    schema_temps_ijk().set_current_time(schema_temps_ijk().get_modified_time_ini());

  if (!schema_temps_ijk().get_first_step_interface_smoothing())
    {
      Cout << "BF posttraiter_champs_instantanes " << schema_temps_ijk().get_current_time() << " " << schema_temps_ijk().get_tstep() << finl;
      get_post().posttraiter_champs_instantanes(lata_name_, schema_temps_ijk().get_current_time(), schema_temps_ijk().get_tstep());
      if (has_thermals_)
        get_ijk_thermals().thermal_subresolution_outputs(); // for thermal counters
      Cout << "AF posttraiter_champs_instantanes" << finl;
    }

  // GB 2019.01.01 Why immobilisation? if (!Option_IJK::DISABLE_DIPHASIQUE && coef_immobilisation_==0.)
  if ((!Option_IJK::DISABLE_DIPHASIQUE) && ns.get_suppression_rejetons())
    interf.detecter_et_supprimer_rejeton(true);

  if (reprise_)
    {
      // On ecrit a la suite du fichier. Cela suppose qu'il est bien a jour.
      // L'instant initial a deja ete ecrit a la fin du calcul precedent donc on ne le reecrit pas.
    }
  else
    {
      // On creer de nouveaux fichiers :
      Cout << "BF ecrire_statistiques_bulles" << finl;
      get_post().ecrire_statistiques_bulles(1 /* reset files */, nom_du_cas(), milieu_ijk().gravite().valeurs(), schema_temps_ijk().get_current_time());
      Cout << "AF ecrire_statistiques_bulles" << finl;
    }

  // Ecrire la valeur initiale dans les sondes :
  // Ecriture de la valeur initiale seulement hors reprise
  if (!reprise_)
    get_post().postraiter_sondes();

  ns.forcage_control_ecoulement();
}

void Probleme_FTD_IJK_base::euler_time_step(ArrOfDouble& var_volume_par_bulle)
{
  // FIXME
  Navier_Stokes_FTD_IJK& ns = ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur());

  static Stat_Counter_Id euler_rk3_counter_ = statistiques().new_counter(2, "Mise a jour de la vitesse");
  statistiques().begin_count(euler_rk3_counter_);
  if (has_thermals_)
    {
      ns.update_v_or_rhov();
      get_ijk_thermals().euler_time_step(schema_temps_ijk().get_timestep());
    }


  ns.euler_time_step(var_volume_par_bulle);
  statistiques().end_count(euler_rk3_counter_);
}

// Perform one sub-step of rk3 for FT algorithm, called 3 times per time step. rk_step = 0, 1 or 2
// total_timestep = not the fractionnal timestep !
void Probleme_FTD_IJK_base::rk3_sub_step(const int rk_step, const double total_timestep, const double fractionnal_timestep, const double time)
{
  assert(rk_step >= 0 && rk_step < 3);
  static Stat_Counter_Id euler_rk3_counter_ = statistiques().new_counter(2, "Mise a jour de la vitesse");
  statistiques().begin_count(euler_rk3_counter_);

  if (has_thermals_)
    get_ijk_thermals().rk3_sub_step(rk_step, total_timestep, time);

  // FIXME
  ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur()).rk3_sub_step(rk_step, total_timestep, fractionnal_timestep, time);

  statistiques().end_count(euler_rk3_counter_);
}

void Probleme_FTD_IJK_base::sauver() const
{
  Probleme_FTD_IJK_base& pb_non_cst = const_cast<Probleme_FTD_IJK_base&>(*this);
  Schema_Temps_IJK_base& sh_non_cst = const_cast<Schema_Temps_IJK_base&>(schema_temps_ijk()); // FIXME
  IJK_Interfaces& interf = const_cast<IJK_Interfaces&>(get_interface());

  sh_non_cst.set_tstep_sauv(schema_temps_ijk().get_tstep() + schema_temps_ijk().get_tstep_init());
  const int dt_sauvegarde = schema_temps_ijk().get_dt_sauvegarde();

  if (schema_temps_ijk().get_tstep_sauv() % dt_sauvegarde == dt_sauvegarde - 1 || stop_)
    {
      // Choix : On supprime les duplicatas pour la sauvegarde.
      // On pourrait tres bien tout garder. ca serait plus leger en CPU, plus lourd en espace disque.
      if (!Option_IJK::DISABLE_DIPHASIQUE)
        interf.supprimer_duplicata_bulles();

      pb_non_cst.sauvegarder_probleme(nom_sauvegarde_, stop_);
      if (!Option_IJK::DISABLE_DIPHASIQUE)
        {
          // On les recree :
          interf.creer_duplicata_bulles();

          // Be on the safe side, on met a jour :
          //   A la suppression des duplicatas, on avait fait mesh.supprimer_facettes qui remet le maillage
          //   a l'etat MINIMAL. Pour les post-tt sur l'interface (eg ai_ft_), il faut que le statut du maillage
          //   soit >= PARCOURU. C'est fait au debut de maj_indicatrice_rho_mu dans
          //   IJK_Interfaces::calculer_indicatrice.
          const double delta_rho = milieu_ijk().get_delta_rho();
          interf.calculer_indicatrice_next(pb_non_cst.get_post().potentiel(), milieu_ijk().gravite().valeurs(), delta_rho, milieu_ijk().sigma(),
                                           schema_temps_ijk().get_current_time(), schema_temps_ijk().get_tstep());
        }
    }
}

int Probleme_FTD_IJK_base::postraiter(int force)
{
  get_post().postraiter_fin(stop_, schema_temps_ijk().get_tstep(), schema_temps_ijk().get_tstep_init(),
                            schema_temps_ijk().get_current_time(), schema_temps_ijk().get_timestep(),
                            lata_name_, milieu_ijk().gravite().valeurs(), nom_du_cas());

  return 1;
}

double Probleme_FTD_IJK_base::computeTimeStep(bool& stop) const
{
  return schema_temps_ijk().computeTimeStep(stop);
}

void Probleme_FTD_IJK_base::validateTimeStep()
{
  // TODO: on pourrait mutualiser tous les parcourir maillages dans IJK_Interface au moment du transport de l'interface
  // et le supprimer de IJK_FT
  // interfaces_.parcourir_maillage();
  if ((!Option_IJK::DISABLE_DIPHASIQUE) && (get_post().get_liste_post_instantanes().contient_("VI")))
    get_interface().compute_vinterp();
}

void Probleme_FTD_IJK_base::terminate()
{
  if (Process::je_suis_maitre())
    {
      SFichier master_file;
      master_file.ouvrir(lata_name_, ios::app);
      master_file << "FIN" << finl;
      master_file.close();
    }
}

bool Probleme_FTD_IJK_base::solveTimeStep()
{
  // Tableau permettant de calculer la variation de volume au cours du pas de temps :
  // Si on veut le mettre en optionel, il faut faire attention a faire vivre la taille de ce tableau avec les
  // creations et destructions de ghosts :
  IJK_Interfaces& interf = get_interface();
  const int nbulles_tot = interf.get_nb_bulles_reelles() + interf.get_nb_bulles_ghost(1/*print=1*/);
  DoubleTrav var_volume_par_bulle(nbulles_tot);
  var_volume_par_bulle = 0.; // Je ne suis pas sur que ce soit un bon choix. Si on ne le remet pas a zero a chaque dt, on corrigera la petite erreur qui pouvait rester d'avant...

  // FIXME
  Navier_Stokes_FTD_IJK& ns = ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur());
  ns.compute_var_volume_par_bulle(var_volume_par_bulle);

  // Au cas ou on soit dans un cas ou des duplicatas sont necessaires mais n'ont pas ete crees, on les cree :
  if (!interf.get_nb_bulles_ghost() && !Option_IJK::DISABLE_DIPHASIQUE)
    interf.creer_duplicata_bulles();

  if ( sub_type(Schema_Euler_explicite_IJK, schema_temps_ijk()) )
    solveTimeStep_Euler(var_volume_par_bulle);
  else if ( sub_type(Schema_RK3_IJK, schema_temps_ijk()) )
    solveTimeStep_RK3(var_volume_par_bulle);
  else
    {
      Cerr << "Erreur dans Probleme_FTD_IJK_base::solveTimeStep() : time_scheme " << schema_temps_ijk().que_suis_je() << " inconnu !" << finl;
      Process::exit();
    }

  ns.corriger_qdm();

  //
  // !!!  TODO - the below should not be done here !!! solveTimeStep() should not update current time:
  //
  schema_temps_ijk().set_current_time(schema_temps_ijk().get_current_time() + schema_temps_ijk().get_timestep()); // update tn

  // stock dans le spliting le decallage periodique total avec condition de shear (current_time_) et celui du pas de temps (timestep_)
  IJK_Shear_Periodic_helpler::shear_x_time_ = ns.get_boundary_conditions().get_dU_perio() *
                                              (schema_temps_ijk().get_current_time() + ns.get_boundary_conditions().get_t0_shear());

  create_forced_dilation(); // rien si Probleme_FTD_IJK_cut_cell

  if (schema_temps_ijk().get_current_time() >= get_post().t_debut_statistiques())
    {
      ns.update_v_or_rhov(true /* with pressure */);

      get_post().update_stat_ft(schema_temps_ijk().get_timestep());
      if (!Option_IJK::DISABLE_DIPHASIQUE)
        get_post().compute_extended_pressures(interf.maillage_ft_ijk());
    }

  return true;
}

void Probleme_FTD_IJK_base::solveTimeStep_Euler(DoubleTrav& var_volume_par_bulle)
{
  // FIXME
  Navier_Stokes_FTD_IJK& ns = ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur());
  IJK_Interfaces& interf = get_interface();

  // Deplacement des interfaces par le champ de vitesse de l'instant n :
  if (!Option_IJK::DISABLE_DIPHASIQUE)
    {
      int counter_first_iter = 1;
      int& first_step_interface_smoothing = schema_temps_ijk().get_first_step_interface_smoothing(); // attention ref

      do
        {
          first_step_interface_smoothing = first_step_interface_smoothing && counter_first_iter;
          deplacer_interfaces(schema_temps_ijk().get_timestep(), -1 /* le numero du sous pas de temps est -1 si on n'est pas en rk3 */, var_volume_par_bulle, first_step_interface_smoothing);
          counter_first_iter--;
          if (first_step_interface_smoothing)
            {
              if (has_thermals_)
                get_ijk_thermals().set_temperature_ini();

              get_post().posttraiter_champs_instantanes(lata_name_, schema_temps_ijk().get_current_time(), schema_temps_ijk().get_tstep());
              ns.compute_var_volume_par_bulle(var_volume_par_bulle);

              if (has_thermals_)
                get_ijk_thermals().set_post_pro_first_call();
            }
        }
      while (first_step_interface_smoothing);

      parcourir_maillage();
    }

  // Mise a jour de la vitesse (utilise les positions des marqueurs, rho, mu et indic a l'instant n)
  // Retourne une vitesse mise a jour et projetee a div nulle
  euler_time_step(var_volume_par_bulle);

  ns.calculer_terme_source_acceleration(schema_temps_ijk().get_current_time() + schema_temps_ijk().get_timestep(),
                                        schema_temps_ijk().get_timestep(),
                                        -1 /* Euler */,
                                        milieu_ijk().get_direction_gravite() /* direction */);

  // Deplacement des interfaces par le champ de vitesse :
  // met a jour la position des marqueurs, la vitesse_ft, et gere les duplicatas.
  // Ne met pas a jour rho_mu_indicatrice

  if (!Option_IJK::DISABLE_DIPHASIQUE) // && !marker_advection_first_)
    {
      // Les sous-pas de temps sont termines. Il n'est plus necessaire de gerer le tableau
      // RK3_G_store_vi_. On peut donc transferer les bulles et re-creer les duplicatas :
      interf.supprimer_duplicata_bulles();
      interf.transferer_bulle_perio();
      // On supprime les fragments de bulles.
      //interfaces_.detecter_et_supprimer_rejeton(false);
      interf.creer_duplicata_bulles();

      ns.test_etapes_et_bilan_rho_u_euler(false /* avant */);
      ns.maj_indicatrice_rho_mu();

      if (has_thermals_)
        get_ijk_thermals().euler_rustine_step(schema_temps_ijk().get_timestep());

      ns.test_etapes_et_bilan_rho_u_euler(true /* apres */);
    }
  else
    ns.test_etapes_et_bilan_rho_u_euler(true /* apres */);
}

void Probleme_FTD_IJK_base::solveTimeStep_RK3(DoubleTrav& var_volume_par_bulle)
{
  Schema_RK3_IJK& rk3 = ref_cast(Schema_RK3_IJK, schema_temps_ijk());
  double current_time_at_rk3_step = rk3.get_current_time();
  // GAB, qdm : passe en attribut de classe car utilise au moment de l'ecriture de mon out
  rk3.get_current_time_at_rk3_step() = rk3.get_current_time();
  // Evaluation de la variation de volume accumule au cours des sous pas de temps.
  // On la laisse croitre pendant les sous dt 0 et 1 puis on la corrige a la fin du 2eme :

  int& rk_step = rk3.get_rk_step();
  const double timestep = rk3.get_timestep(), current_time = rk3.get_current_time();

  // FIXME
  Navier_Stokes_FTD_IJK& ns = ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur());
  IJK_Interfaces& interf = get_interface();

  for (rk_step = 0; rk_step < 3; rk_step++)
    {
      const double fractionnal_timestep = compute_fractionnal_timestep_rk3(timestep /* total*/, rk_step);

      // Mise a jour des positions des marqueurs.
      // Deplacement des interfaces par le champ de vitesse au sous pas de temps k :
      if (!Option_IJK::DISABLE_DIPHASIQUE)
        {
          deplacer_interfaces_rk3(timestep /* total */, rk_step, var_volume_par_bulle);
          parcourir_maillage();
        }
      // Mise a jour de la temperature et de la vitesse :
      rk3_sub_step(rk_step, timestep, fractionnal_timestep, current_time_at_rk3_step);

      ns.test_etapes_et_bilan_rho_u_euler(false /* avant */);

      // Mise a jour rho, mu et l'indicatrice a partir de la nouvelle position de l'interface :
      // (sauf au dernier sous pas de temps pour lequel c'est fait a la fin du pas de temps)
      // TODO: verifier qu'on doit bien le faire aussi au dernier sous pas de temps : rk_step != 2 &&
      // TODO aym: verifier ce bloc, qui applique les sous pas de temps RK3 de la rustine a la temperature
      if (rk_step != 2 && !Option_IJK::DISABLE_DIPHASIQUE)
        {
          // Attention, il faut que les duplicatas soient present pour faire maj_indicatrice_rho_mu :
          ns.maj_indicatrice_rho_mu();

          if (has_thermals_)
            get_ijk_thermals().rk3_rustine_sub_step(rk_step, timestep, fractionnal_timestep, current_time_at_rk3_step);
        }

      ns.calculer_terme_source_acceleration(current_time_at_rk3_step, timestep /*total*/, rk_step, milieu_ijk().get_direction_gravite()  /* direction */);

      current_time_at_rk3_step += fractionnal_timestep;
      // GAB, qdm : passe en attribut de classe car utilise au moment de l'ecriture de mon out
      rk3.get_current_time_at_rk3_step() += fractionnal_timestep;

      // On ne postraite pas le sous-dt 2 car c'est fait plus bas si on post-traite le pas de temps :
      if (get_post().postraiter_sous_pas_de_temps()
          && ((rk3.get_tstep() % get_post().nb_pas_dt_post() == get_post().nb_pas_dt_post() - 1)
              || (std::floor((current_time - timestep) / get_post().get_timestep_simu_post(current_time, rk3.get_max_simu_time()))
                  < std::floor(current_time / get_post().get_timestep_simu_post(current_time, rk3.get_max_simu_time())))) && (rk_step != 2))
        {
          get_post().posttraiter_champs_instantanes(lata_name_, current_time_at_rk3_step, rk3.get_tstep());
        }
    }
  if (!Option_IJK::DISABLE_DIPHASIQUE)
    {
      // Les sous-pas de temps sont termines. Il n'est plus necessaire de gerer le tableau
      // RK3_G_store_vi_. On peut donc transferer les bulles et re-creer les duplicatas :
      interf.supprimer_duplicata_bulles();
      interf.transferer_bulle_perio();
      // On supprime les fragments de bulles.
      //interfaces_.detecter_et_supprimer_rejeton(false);
      interf.creer_duplicata_bulles();

      // Mise a jour rho, mu et l'indicatrice a partir de la nouvelle position de l'interface :
      ns.maj_indicatrice_rho_mu();

      if (has_thermals_)
        get_ijk_thermals().update_thermal_properties();
    }

  ns.test_etapes_et_bilan_rho_u_euler(true /* apres */);
}

bool Probleme_FTD_IJK_base::run()
{
  Cerr << "Probleme_FTD_IJK_base::run()" << finl;

  preparer_calcul();

  schema_temps_ijk().set_max_timestep(schema_temps_ijk().get_timestep());

  statistiques().end_count(initialisation_calcul_counter_);

  if (!disable_TU)
    {
      if (GET_COMM_DETAILS)
        statistiques().print_communciation_tracking_details("Statistiques d'initialisation du calcul", 0);

      statistiques().dump("Statistiques d'initialisation du calcul", 0);
      print_statistics_analyse("Statistiques d'initialisation du calcul", 0);
    }
  statistiques().reset_counters();
  statistiques().begin_count(temps_total_execution_counter_);

  bool ok = true;
  int& tstep = schema_temps_ijk().get_tstep();

  for (tstep = 0; tstep < schema_temps_ijk().get_nb_timesteps() && !stop_; tstep++)
    {
      statistiques().begin_count(timestep_counter_);

      schema_temps_ijk().set_timestep() = computeTimeStep(stop_);

      /* Contrary to what is done in the classical Probleme_base::run() method, here we do not stop directly
       * we still go through the end of the loop:
       */
      //      if (stop_) /* stop file detected ? */
      //        break;

      // Prepare the next time step
      if (!initTimeStep(schema_temps_ijk().get_timestep()))
        return false;

      // Solve the next time step
      ok = solveTimeStep();

      // Should we stop the computation ? TODO This should be all done in computeTimeStep above ... not ICoCo compliant at the moment ...
      schema_temps_ijk().check_stop_criteria(stop_);

      sauver();

      if (!ok)   // The resolution failed, try with a new time interval.
        {
          abortTimeStep();
          schema_temps_ijk().set_timestep() = computeTimeStep(stop_);
        }
      else // The resolution was successful, validate and go to the next time step.
        validateTimeStep();

      postraiter(stop_);

      statistiques().end_count(timestep_counter_);

      if (JUMP_3_FIRST_STEPS && tstep < 3)
        {
          //demarrage des compteurs CPU
          if (tstep == 2)
            statistiques().set_three_first_steps_elapsed(true);
        }
      else
        statistiques().compute_avg_min_max_var_per_step(tstep);
    }

  if (!disable_TU)
    {
      if (GET_COMM_DETAILS)
        statistiques().print_communciation_tracking_details("Statistiques de resolution du probleme", 1);

      statistiques().dump("Statistiques de resolution du probleme", 1);
      print_statistics_analyse("Statistiques de resolution du probleme", 1);
    }

  statistiques().reset_counters();
  statistiques().begin_count(temps_total_execution_counter_);

  return ok;
}

void Probleme_FTD_IJK_base::build_vdf_domaine()
{
  Cerr << "Construction du domaine VDF NS pour les sondes..." << finl;
  refprobleme_ns_ = creer_domaine_vdf(domaine_ijk_.valeur(), "DOM_NS_VDF");

  Cerr << "Construction du domaine VDF..." << finl;
  build_extended_splitting(domaine_ijk_.valeur(), domaine_ft_, domaine_ijk_->ft_extension());
  refprobleme_ft_disc_ = creer_domaine_vdf(domaine_ft_, "DOM_VDF");

  // FIXME
  ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur()).build_redistribute_extended_splitting_ft();
}
