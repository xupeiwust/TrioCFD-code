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

#include <Check_espace_virtuel.h>
#include <Postprocessing_IJK.h>
#include <Schema_Euler_explicite_IJK.h>
#include <IJK_Navier_Stokes_tools.h>
#include <Probleme_FTD_IJK_base.h>
#include <IJK_Thermal_cut_cell.h>
#include <IJK_Field_vector.h>
#include <IJK_Lata_writer.h>
#include <Schema_RK3_IJK.h>
#include <Domaine_IJK.h>
#include <IJK_Interfaces.h>
#include <Option_IJK.h>

Implemente_instanciable_sans_constructeur(Postprocessing_IJK, "Postprocessing_IJK", Postraitement_ft_lata);

// list of fields that may be postprocessed
std::vector<Postprocessing_IJK::FieldInfo_t> Postprocessing_IJK::champs_postraitables_ = {};


namespace
{
// Could not find it elsewhere but surely must already exist?
inline Entity str_to_entity(const Motcle& loc)
{
  if (loc == "ELEM") return Entity::ELEMENT;
  if (loc == "SOM") return Entity::NODE;
  if (loc == "FACES") return Entity::FACE;

  Cerr << "Invalid localisation for field postprocessing : '" << loc << "' !!" << finl;
  Process::exit();
  return Entity::ELEMENT; // for compilers
}

inline Nom entity_to_str(const Entity& e)
{
  if (e == Entity::ELEMENT) return "ELEM";
  if (e == Entity::NODE) return "SOM";
  if (e == Entity::FACE) return "FACES";
  return "??";
}

inline int extract_component(const Nom& fld_name)
{
  const std::vector<std::string> compos = {"_X", "_Y", "_Z"};
  std::string s = fld_name.getString();
  std::string end = s.substr(s.size()-2, 2);
  auto it = std::find(compos.begin(), compos.end(), end);
  if(it == compos.end())
    {
      Cerr << "ERROR: field name '" << fld_name << "' does not end with a component name (_X, _Y or _Z)!! Should not happend?" << finl;
      Process::exit();
    }
  return (int)std::distance(compos.begin(), it);
}

}

Postprocessing_IJK::Postprocessing_IJK()
{
  groups_statistiques_FT_.dimensionner(0);
}

Sortie& Postprocessing_IJK::printOn(Sortie& os) const { return os; }

Entree& Postprocessing_IJK::readOn(Entree& is)
{
  return Postraitement_base::readOn(is);
}

void Postprocessing_IJK::associer_probleme(const Probleme_FTD_IJK_base& ijk_ft)
{
  ref_ijk_ft_ = ijk_ft;
  statistiques_FT_.associer_probleme(ijk_ft);
  interfaces_ = ref_ijk_ft_->get_interface();
  pressure_ = ref_ijk_ft_->eq_ns().pressure_;
  velocity_ = ref_ijk_ft_->eq_ns().velocity_;
  source_spectrale_ = ref_ijk_ft_->eq_ns().forcage_.get_force_ph2();
  bk_tsi_ns_ = ref_ijk_ft_->eq_ns().backup_terme_source_interfaces_ns_;
  d_velocity_ = ref_ijk_ft_->eq_ns().d_velocity_;

  if (ref_ijk_ft_->has_thermals())
    thermals_ = ref_ijk_ft_->get_ijk_thermals();
}

void Postprocessing_IJK::set_param(Param& param)
{
  Postraitement_ft_lata::set_param(param);

  param.ajouter("nb_pas_dt_post_thermals_probes", &nb_pas_dt_post_thermals_probes_);
  param.ajouter("nb_pas_dt_post_stats_bulles", &nb_pas_dt_post_stats_bulles_);
  param.ajouter("nb_pas_dt_post_stats_plans", &nb_pas_dt_post_stats_plans_);
  param.ajouter("nb_pas_dt_post_stats_cisaillement", &nb_pas_dt_post_stats_cisaillement_);
  param.ajouter("nb_pas_dt_post_stats_rmf", &nb_pas_dt_post_stats_rmf_);

  param.ajouter("time_interval_post_", &time_interval_post_);
  param.ajouter("time_interval_post_thermals_probes", &time_interval_post_thermals_probes_);
  param.ajouter("time_interval_post_stats_bulles", &time_interval_post_stats_bulles_);
  param.ajouter("time_interval_post_stats_plans", &time_interval_post_stats_plans_);
  param.ajouter("time_interval_post_stats_cisaillement", &time_interval_post_stats_cisaillement_);
  param.ajouter("time_interval_post_stats_rmf", &time_interval_post_stats_rmf_);

  param.ajouter("champs_a_postraiter", &liste_post_instantanes_);

  param.ajouter_flag("check_stats", &check_stats_);
  param.ajouter_flag("postraiter_sous_pas_de_temps", &postraiter_sous_pas_de_temps_);
  // Pour reconstruire au post-traitement la grandeur du/dt, on peut choisir de relever u^{dt_post} et u^{dt_post+1} :
  param.ajouter_flag("post_par_paires", &post_par_paires_);

  expression_vitesse_analytique_.dimensionner_force(3);
  param.ajouter("expression_vx_ana", &expression_vitesse_analytique_[0]);
  param.ajouter("expression_vy_ana", &expression_vitesse_analytique_[1]);
  param.ajouter("expression_vz_ana", &expression_vitesse_analytique_[2]);

  expression_dvitesse_analytique_.dimensionner_force(3);
  param.ajouter("expression_dvx_ana", &expression_dvitesse_analytique_[0]);
  param.ajouter("expression_dvy_ana", &expression_dvitesse_analytique_[1]);
  param.ajouter("expression_dvz_ana", &expression_dvitesse_analytique_[2]);
  param.ajouter("expression_p_ana", &expression_pression_analytique_);

  expression_gradP_analytique_.dimensionner_force(3);
  param.ajouter("expression_dPdx_ana", &expression_gradP_analytique_[0]);
  param.ajouter("expression_dPdy_ana", &expression_gradP_analytique_[1]);
  param.ajouter("expression_dPdz_ana", &expression_gradP_analytique_[2]);

  expression_gradU_analytique_.dimensionner_force(3);
  param.ajouter("expression_dUdx_ana", &expression_gradU_analytique_[0]);
  param.ajouter("expression_dUdy_ana", &expression_gradU_analytique_[1]);
  param.ajouter("expression_dUdz_ana", &expression_gradU_analytique_[2]);
  expression_gradV_analytique_.dimensionner_force(3);
  param.ajouter("expression_dVdx_ana", &expression_gradV_analytique_[0]);
  param.ajouter("expression_dVdy_ana", &expression_gradV_analytique_[1]);
  param.ajouter("expression_dVdz_ana", &expression_gradV_analytique_[2]);
  expression_gradW_analytique_.dimensionner_force(3);
  param.ajouter("expression_dWdx_ana", &expression_gradW_analytique_[0]);
  param.ajouter("expression_dWdy_ana", &expression_gradW_analytique_[1]);
  param.ajouter("expression_dWdz_ana", &expression_gradW_analytique_[2]);

  // Pour les seconds gradients :
  expression_grad2P_analytique_.dimensionner_force(6);
  param.ajouter("expression_ddPdxdx_ana", &expression_grad2P_analytique_[0]);
  param.ajouter("expression_ddPdydy_ana", &expression_grad2P_analytique_[1]);
  param.ajouter("expression_ddPdzdz_ana", &expression_grad2P_analytique_[2]);
  param.ajouter("expression_ddPdxdy_ana", &expression_grad2P_analytique_[3]);
  param.ajouter("expression_ddPdxdz_ana", &expression_grad2P_analytique_[4]);
  param.ajouter("expression_ddPdydz_ana", &expression_grad2P_analytique_[5]);
  // And for velocities :
  expression_grad2U_analytique_.dimensionner_force(6);
  param.ajouter("expression_ddUdxdx_ana", &expression_grad2U_analytique_[0]);
  param.ajouter("expression_ddUdydy_ana", &expression_grad2U_analytique_[1]);
  param.ajouter("expression_ddUdzdz_ana", &expression_grad2U_analytique_[2]);
  param.ajouter("expression_ddUdxdy_ana", &expression_grad2U_analytique_[3]);
  param.ajouter("expression_ddUdxdz_ana", &expression_grad2U_analytique_[4]);
  param.ajouter("expression_ddUdydz_ana", &expression_grad2U_analytique_[5]);

  expression_grad2V_analytique_.dimensionner_force(6);
  param.ajouter("expression_ddVdxdx_ana", &expression_grad2V_analytique_[0]);
  param.ajouter("expression_ddVdydy_ana", &expression_grad2V_analytique_[1]);
  param.ajouter("expression_ddVdzdz_ana", &expression_grad2V_analytique_[2]);
  param.ajouter("expression_ddVdxdy_ana", &expression_grad2V_analytique_[3]);
  param.ajouter("expression_ddVdxdz_ana", &expression_grad2V_analytique_[4]);
  param.ajouter("expression_ddVdydz_ana", &expression_grad2V_analytique_[5]);

  expression_grad2W_analytique_.dimensionner_force(6);
  param.ajouter("expression_ddWdxdx_ana", &expression_grad2W_analytique_[0]);
  param.ajouter("expression_ddWdydy_ana", &expression_grad2W_analytique_[1]);
  param.ajouter("expression_ddWdzdz_ana", &expression_grad2W_analytique_[2]);
  param.ajouter("expression_ddWdxdy_ana", &expression_grad2W_analytique_[3]);
  param.ajouter("expression_ddWdxdz_ana", &expression_grad2W_analytique_[4]);
  param.ajouter("expression_ddWdydz_ana", &expression_grad2W_analytique_[5]);

  param.ajouter("t_debut_statistiques", &t_debut_statistiques_);

  // mots clés pour la reprise de stats
  param.ajouter_flag("reset_reprise_integrated", &reset_reprise_integrated_);
}

void Postprocessing_IJK::lire_entete_bloc_interface(Entree& is)
{
  Motcle motlu;
  is >> motlu;

  if (Process::je_suis_maitre())
    Cerr << "Post-processing for the interface of IJK_interfaces object named: '" << motlu << "'" << finl;

  // Check valid interface equation name:
  Motcle interf_nam;
  for (int i=0; i < mon_probleme->nombre_d_equations(); i++)
    {
      const Equation_base& eb = mon_probleme->equation(i);
      if (motlu == eb.le_nom() && eb.que_suis_je() == "IJK_Interfaces")
        interf_nam = eb.le_nom();
    }

  if (motlu != interf_nam)
    {
      Cerr << "ERROR: the requested interface equation '" << motlu << "' is not available for post-processing!!!" << finl;
      Process::exit();
    }

  is >> motlu;
  if (motlu != "{")
    {
      Cerr << "ERROR: Postprocessing_IJK::lire_entete_bloc_interface()\n";
      Cerr << " { was expected after the keyword interfaces\n";
      Cerr << " We found '" << motlu << "'" << finl;
      Process::exit();
    }
}

void Postprocessing_IJK::register_one_field(const Motcle& fld_nam, const Motcle& reqloc_s)
{
  Entity reqloc = str_to_entity(reqloc_s);
  noms_champs_a_post_.add(fld_nam);

  Motcle true_field_name = fld_nam;

  // IJK_Thermals gives a list of possible prefixes. The field names are appended with _index when added to champs_compris by individual thermal equations
  // Because of that, we have to get the prefix from the field name asked
  const std::string& str = fld_nam.getString();
  std::string key ("_");
  std::size_t found = str.rfind(key);
  if (found!=std::string::npos)
    {
      Nom var = Nom(str.substr(0,found));
      found++;
      try
        {
          // this may raise exception if field name did not finish with integer. Then we go to catch block below
          std::stoi(str.substr(found));
          // Cerr << "field was suffixed with an integer: " << fld_nam<<finl;
          // if field was ending with integer, check that the prefix is known by ijk thermals
          std::vector<FieldInfo_t> prefix_thermals;
          IJK_Thermals::Fill_postprocessable_fields(prefix_thermals);
          for (const auto& f : prefix_thermals)
            {
              Nom prefix = std::get<0>(f);
              if (prefix == var)
                {
                  // if we found a matching prefix, we know that the field location must be defined by the prefix
                  true_field_name=prefix;
                  Cerr << "   For IJK_thermals: postpro field name '" << fld_nam << "' matched with prefix " << true_field_name << " defined in IJK_thermals::Fill_postprocessable_fields " << finl;

                  // then we must add an entry to champs_postraitables_ so that the field can be found
                  std::vector<FieldInfo_t> c = {{ fld_nam, std::get<1>(f), std::get<2>(f), std::get<3>(f)}};
                  champs_postraitables_.insert(champs_postraitables_.end(), c.begin(), c.end());
                }
            }
        }
      catch ( const std::invalid_argument& )
        {
          // in that case, everything is fine, field name should be treated as normal (not from thermals)
          // Cerr << "field was not suffixed with an integer: " << fld_nam<<finl;
        }
    }

  // Lookup the field in champs_postraitables_
  int idx = 0;
  for (const auto& f : champs_postraitables_)
    {
      Motcle nom2 = std::get<0>(f);
      Entity loc2 = std::get<1>(f);

      // Prepare entries in storage map (some of them might not be used for example when postprocessing
      // an unknown directly):
      // Teo Boutin: moved here because sometime a post may depend on another.
      // Then some fields in the map may not be prepared
      // For example: if CURL is asked in post but not CRITERE_Q, we try in alloc_fields to allocate scalar_post_fields_.at("CRITERE_Q") which was not pre initialized.
      // now we pre fill for every post processable field. Probably there is a better solution
      if (!scalar_post_fields_.count(nom2))
        scalar_post_fields_[nom2] = IJK_Field_double();
      if(!vect_post_fields_.count(nom2))
        vect_post_fields_[nom2] = IJK_Field_vector3_double();

      bool want_interp = (reqloc == Entity::ELEMENT && loc2 == Entity::FACE);
      if(nom2 == fld_nam
          && (reqloc == loc2 || want_interp)) // allow FACE to ELEMENT interpolation
        {
          // If already there, error:
          FieldIndex_t k = {idx, want_interp};
          if (std::find(field_post_idx_.begin(), field_post_idx_.end(), k) != field_post_idx_.end())
            {
              Cerr << "ERROR: field '" << fld_nam << "' at localisation '"<< reqloc_s <<"' duplicated in list of fields to be postprocessed!!" << finl;
              Process::exit();
            }
          field_post_idx_.push_back(k);
          list_post_required_.push_back(fld_nam);
          break;
        }

      idx++;
    }

  if (idx == (int)champs_postraitables_.size())
    {
      Cerr << "ERROR: field '" << fld_nam << "' at localisation '"<< reqloc_s <<"' is not available for postprocessing!!" << finl;
      Process::exit();
    }
}

/** Override. Called from base class.
 */
void Postprocessing_IJK::register_interface_field(const Motcle& nom_champ, const Motcle& loc_lu)
{
  Entity e = str_to_entity(loc_lu);

  // Check the field is indeed on the interface:
  std::vector<FieldInfo_t> flds;
  IJK_Interfaces::Fill_postprocessable_fields(flds);
  bool ok=false;
  for (const auto& f : flds)
    {
      Motcle nom2 = get<0>(f);
      Entity loc2 = get<1>(f);
      bool isOnInterface = get<3>(f); // must check that it is truly on interface, as IJK_Interfaces also defines eulerian fields
      if(isOnInterface && nom2 == nom_champ && e == loc2)
        {
          ok = true;
          break;
        }
    }
  if(!ok)
    {
      Cerr << "ERROR: on the interface, field '" << nom_champ << "' at localisation '"<< loc_lu <<"' is not available for postprocessing!!" << finl;
      Process::exit();
    }
  // Everything ok, we register the field for post:
  register_one_field(nom_champ, loc_lu);
  interface_post_required_ = true;
}

/** Override to have a simpler logic than base class. We really want to retrieve names + location.
 */
int Postprocessing_IJK::lire_champs_a_postraiter(Entree& is, bool expect_acco)
{
  Motcle accolade_fermee("}"), motlu;

  is >> motlu;

  while (motlu != accolade_fermee)
    {
      Motcle loc;
      is >> loc;
      register_one_field(motlu, loc);
      is >> motlu;
    }

  return 1;
}

/** Initialise lata file and various other stuff
 */
void Postprocessing_IJK::init()
{
  if(!nom_fich_.finit_par(".lata"))
    nom_fich_ = nom_fich_ + Nom(".lata");

  // Post_processing field allocations:
  alloc_fields();
  alloc_velocity_and_co();

  // Integrated field initialisation
  init_integrated_and_ana(ref_ijk_ft_->get_reprise());

  completer_sondes();

  prepare_lata_and_stats();

  compute_extended_pressures();
}

/**
 * Write the master lata file and prepare statistics and other stuff
 */
void Postprocessing_IJK::prepare_lata_and_stats()
{
  const Nom& lata_name = nom_fich_;
  const double current_time = ref_ijk_ft_->schema_temps_ijk().get_current_time();
  dumplata_header(lata_name);
  dumplata_add_geometry(lata_name, velocity_.valeur()[0]);
  dumplata_add_geometry(lata_name, ref_ijk_ft_->eq_ns().velocity_ft_[0]);

  // Calcul des moyennes spatiales sur la condition initiale:
  if (is_stats_plans_activated() && current_time >= t_debut_statistiques_)
    {
      // FA AT 16/07/2013 pensent que necessaire pour le calcul des derivees dans statistiques_.update_stat_k(...)
      // Je ne sais pas si c'est utile, mais j'assure...
      velocity_.valeur()[0].echange_espace_virtuel(2 /*, IJK_Field_ST::EXCHANGE_GET_AT_RIGHT_I*/);
      velocity_.valeur()[1].echange_espace_virtuel(2 /*, IJK_Field_ST::EXCHANGE_GET_AT_RIGHT_J*/);
      velocity_.valeur()[2].echange_espace_virtuel(2 /*, IJK_Field_ST::EXCHANGE_GET_AT_RIGHT_K*/);
      pressure_->echange_espace_virtuel(1);

      // C'est update_stat_ft qui gere s'il y a plusieurs groupes
      // pour faire la vraie indicatrice + les groupes
      update_stat_ft(0.);
    }
}

void Postprocessing_IJK::postraiter(int forcer)
{
  // Take care of fields and probes - new mode
  Postraitement_ft_lata::postraiter(forcer);

//  // Take care of fields - OLD MODE:
//  {
//    Schema_Temps_IJK_base& sch = ref_ijk_ft_->schema_temps_ijk();
//    Cout << "BF posttraiter_champs_instantanes " << sch.get_current_time() << " " << sch.get_tstep() << finl;
//    posttraiter_champs_instantanes(nom_fich_, sch.get_current_time(), sch.get_tstep());
//    if (ref_ijk_ft_->has_thermals())
//      ref_ijk_ft_->get_ijk_thermals().thermal_subresolution_outputs(); // for thermal counters
//    Cout << "AF posttraiter_champs_instantanes" << finl;
//  }

  // Thermal part
  postraiter_thermals(forcer);

  // Statistics (bubbles, etc.)
  postraiter_stats(forcer);

  // Probes update -
  // TODO : not necessary? done in base class?
  Schema_Temps_IJK_base& sch = ref_ijk_ft_->schema_temps_ijk();
  double current_time = sch.get_current_time(),
         timestep = sch.get_timestep();
  les_sondes_.mettre_a_jour(current_time, timestep);
}

/** Override. Write the interface mesh if present, and the integer field 'COMPO_CONNEXE' on it.
 */
int Postprocessing_IJK::write_extra_mesh()
{
  if(!interface_post_required_ || !ref_ijk_ft_->has_interface())
    return 1;
  Schema_Temps_IJK_base& sch = ref_ijk_ft_->schema_temps_ijk();
  int latastep = sch.get_tstep();
  const IJK_Interfaces& interf = ref_ijk_ft_->get_interface();
  interf.dumplata_ft_mesh(nom_fich_.getChar(), "INTERFACES", latastep);

  // Writing systematically COMPO_CONNEXE, the only integer field:
  const ArrOfInt& comp_c = ref_ijk_ft_->get_interface().maillage_ft_ijk().compo_connexe_facettes();
  dumplata_ft_field(nom_fich_, "INTERFACES", "COMPO_CONNEXE", "ELEM", comp_c, latastep);

  return 1;
}

static inline bool check_loc_compat(Entity requested_loc, Domaine_IJK::Localisation loc_ijk, bool want_interp)
{
  using LocIJK = Domaine_IJK::Localisation;
  if (requested_loc == Entity::ELEMENT && loc_ijk == LocIJK::ELEM)
    return true;
  if (requested_loc == Entity::NODE && loc_ijk == LocIJK::NODES)
    return true;
  if ((loc_ijk == LocIJK::FACES_I || loc_ijk == LocIJK::FACES_J || loc_ijk == LocIJK::FACES_K))
    {
      if(requested_loc == Entity::FACE) return true;
      if(requested_loc == Entity::ELEMENT)  // natural localisation is face, user wants interpolation to ELEM - OK
        {
          assert(want_interp);
          return true;
        }
    }
  return false;
}

/*! Override from 'Postraitement' since the logic is simpler here
 */
int Postprocessing_IJK::postraiter_champs()
{
  int latastep = ref_ijk_ft_->schema_temps_ijk().get_tstep();

  // Write out a new time in the lata:
  dumplata_newtime(nom_fich_, ref_ijk_ft_->schema_temps_ijk().get_current_time());

  // Write INTERFACES if there:
  write_extra_mesh();

  // Dump requested fields:
  for (const FieldIndex_t& v : field_post_idx_)
    {
      int idx = get<0>(v);
      bool want_interpolation = get<1>(v);
      const FieldInfo_t& fld_info = champs_postraitables_[idx];
      Motcle fld_nam = get<0>(fld_info);
      Entity fld_loc = get<1>(fld_info);  // the natural localisation of field
      Nature_du_champ nat = get<2>(fld_info);
      bool is_on_interf = get<3>(fld_info);

      if (!ref_ijk_ft_->has_champ(fld_nam) && !ref_ijk_ft_->has_champ_vectoriel(fld_nam) )
        {
          Cerr << finl << "ERROR Field '" << fld_nam << "' is not available for postprocessing!" << finl << finl;
          Cerr << "(Fields currently available are : " << finl;
          Noms ns;
          ref_ijk_ft_->get_noms_champs_postraitables(ns, Option::DESCRIPTION);
          for (const auto& n: ns) Cerr << n << " ";
          Cerr << ")" << finl;
          Process::exit();
        }

      if (nat == Nature_du_champ::scalaire)
        {
          const IJK_Field_double* fld = &(ref_ijk_ft_->get_IJK_field(fld_nam));
          assert(fld->nature_du_champ() == Nature_du_champ::scalaire);

          if (is_on_interf)     // for now, loc validity is not checked for interf fields ...
            dumplata_ft_field(nom_fich_, "INTERFACES", fld_nam, entity_to_str(fld_loc), fld->data(), latastep);
          else
            {
              if (!check_loc_compat(fld_loc, fld->get_localisation(), want_interpolation))
                {
                  Cerr << "ERROR Field '" << fld_nam << "' is available for postprocessing, but NOT at localisation '" << entity_to_str(fld_loc) << "'!" << finl;
                  Process::exit();
                }
              if(want_interpolation)
                {
                  // First time we interpolate - needs to allocate storage:
                  if (post_projected_field_.get_ptr(0) == nullptr)
                    allocate_cell_vector(post_projected_field_, domaine_ijk_, 0);
                  // Identify the correct component to use for temp storage (we could have had separate temporary members for
                  // component interpolation, but I don't think this would have made the below any better/shorter ...):
                  int compo_num = extract_component(fld_nam);
                  interpolate_to_center_compo(post_projected_field_[compo_num], *fld);
                  post_projected_field_[compo_num].dumplata_scalar(nom_fich_, latastep);
                }
              else
                fld->dumplata_scalar(nom_fich_, latastep);
            }
        }
      else if (nat == Nature_du_champ::vectoriel)
        {
          if(is_on_interf)
            {
              Cerr << "ERROR: post-processing of vectorial field on the interface not implemented (here '" << fld_nam << "') !" << finl;
              Process::exit();
            }

          const IJK_Field_vector3_double& fld = ref_ijk_ft_->get_IJK_field_vector(fld_nam);
          const Entity loc = fld.localisation();
          assert(fld.nature_du_champ() == Nature_du_champ::vectoriel);
          if (loc == Entity::ELEMENT && fld_loc == Entity::FACE)
            {
              Cerr << "ERROR Field '" << fld_nam << "' is available for postprocessing on '" << entity_to_str(loc) << "' and can not be projected to '" << entity_to_str(fld_loc) << "'" << finl;
              Process::exit();
            }
          if (want_interpolation)   // needs interpolation
            {
              assert(loc == Entity::FACE && fld_loc == Entity::FACE);
              // First time we interpolate - needs to allocate storage:
              if (post_projected_field_.get_ptr(0) == nullptr)
                allocate_cell_vector(post_projected_field_, domaine_ijk_, 0);
              // TODO: pour ABN from GB : problem avec le const ! -> faire un ref_cast_non_const tout moche?
              // for (int dir2 = 0; dir2 < 3; dir2++)
              //   fld[dir2].echange_espace_virtuel(fld[dir2].ghost());
              /* Exit in error if the virtual spaces of the distributed array are not up to date */
              // for (int dir2 = 0; dir2 < 3; dir2++) assert(check_espace_virtuel_vect(fld[dir2]));
              interpolate_to_center(post_projected_field_, fld);
              dumplata_cellvector(nom_fich_,Nom("CELL_") + fld_nam, post_projected_field_, latastep);
            }
          if(fld_loc != loc) // Incompatible localisations, not interpolation possible
            {
              Cerr << "ERROR Field '" << fld_nam << "' is available for postprocessing, but NOT at localisation '" << entity_to_str(fld_loc) << "' !" << finl;
              Process::exit();
            }

          switch(loc)
            {
            case Entity::NODE:
              Process::exit("ERROR: post-processing vectorial field on NODES not implemented!");
              break;
            case Entity::ELEMENT:
              dumplata_cellvector(nom_fich_, fld_nam, fld, latastep);
              break;
            case Entity::FACE:
              dumplata_vector(nom_fich_, fld_nam, fld[0], fld[1], fld[2], latastep);
              break;
            default:
              Process::exit("ERROR: Unexpected error!");
            }
        } // if(nat == vectoriel)
      else // nat ==  Nature_du_champ::multi_scalaire
        Process::exit("ERROR: no multiscalar support!!");
    }
  return 1;
}

void Postprocessing_IJK::associer_domaines(Domaine_IJK& dom_ijk, Domaine_IJK& dom_ft)
{
  domaine_ijk_ = dom_ijk;
  domaine_ft_ = dom_ft;
  statistiques_FT_.associer_domaine(dom_ijk);
}

void Postprocessing_IJK::init_integrated_and_ana(bool reprise)
{
  Navier_Stokes_FTD_IJK& ns = ref_ijk_ft_->eq_ns();


  // En reprise, il se peut que le champ ne soit pas dans la liste des posts, mais qu'on l'ait quand meme.
  // Dans ce cas, on choisi de le lire, remplir le field et le re-sauvegarder a la fin (on n'en a rien fait de plus entre temps...)
  if (
    (( ns.coef_immobilisation_ > 1e-16) && is_stats_plans_activated())
    || is_post_required("INDICATRICE_PERTURBE")
    || ((reprise) && ((ns.fichier_reprise_vitesse_ != "??"))) // TODO teo boutin: read from fichier_reprise_interface maybe ??
  )
    {
      indicatrice_non_perturbe_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 0, "INDICATRICE_PERTURBE");
      champs_compris_.ajoute_champ(indicatrice_non_perturbe_);
      fill_indic(ref_ijk_ft_->get_reprise());
    }

  // pour relire les champs de temps integres:
  if (is_post_required("INTEGRATED_TIMESCALE"))
    {
      integrated_timescale_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 0, "INTEGRATED_TIMESCALE");
      champs_compris_.ajoute_champ(integrated_timescale_);
      if ((reprise) && (!reset_reprise_integrated_))
        {
          if (ns.fichier_reprise_vitesse_ == "??")
            {
              Cerr << "fichier_reprise_vitesse_ should be specified in the restart file in Navier_Stokes_FTD_IJK object" << endl;
              Process::exit();
            }
          const int timestep_reprise_integrated_timescale_ = 1;
          const Nom& geom_name = integrated_timescale_.get_domaine().le_nom();

          if (!lata_has_field(ns.fichier_reprise_vitesse_, timestep_reprise_integrated_timescale_, geom_name, "INTEGRATED_TIMESCALE"))
            {
              Cerr << "fichier_reprise_vitesse_ " <<  ns.fichier_reprise_vitesse_ << " does not contain field INTEGRATED_TIMESCALE to restart statistics. You may specify the flag reset_reprise_integrated in postprocessing to reset statistics" << endl;
              Process::exit();
            }
          lire_dans_lata(ns.fichier_reprise_vitesse_, timestep_reprise_integrated_timescale_, geom_name, "INTEGRATED_TIMESCALE", integrated_timescale_); // fonction qui lit un champ a partir d'un lata .
        }
      else
        {
          integrated_timescale_.data() = 0.;
          // Question GB pour Antoine : c'est une precaution pour pas qu'il vaille 0 au debut?
          // Mais du coup, on le compte deux fois...
          update_integral_indicatrice(interfaces_->In(), 1. /* Should be the integration timestep */, integrated_timescale_);
        }
    }

  // Pour relire les champs de vitesse et pression integres :
  if ((( ns.coef_immobilisation_ > 1e-16) && is_stats_plans_activated()) || is_post_required("INTEGRATED_VELOCITY"))
    {
      allocate_velocity(integrated_velocity_, domaine_ijk_, 2, "INTEGRATED_VELOCITY");
      champs_compris_.ajoute_champ_vectoriel(integrated_velocity_);
    }
  if (is_post_required("INTEGRATED_VELOCITY"))
    {
      if ((reprise) && (!reset_reprise_integrated_))
        {
          if (ns.fichier_reprise_vitesse_ == "??")
            {
              Cerr << "fichier_reprise_vitesse_ should be specified in the restart file in Navier_Stokes_FTD_IJK object" << endl;
              Process::exit();
            }
          const int timestep_reprise_integrated_velocity_ = 1;
          const Nom& geom_name = velocity_.valeur()[0].get_domaine().le_nom();

          if (!lata_has_field(ns.fichier_reprise_vitesse_, timestep_reprise_integrated_velocity_, geom_name, "INTEGRATED_VELOCITY"))
            {
              Cerr << "fichier_reprise_vitesse_ " <<  ns.fichier_reprise_vitesse_ << " does not contain field INTEGRATED_VELOCITY to restart statistics. You may specify the flag reset_reprise_integrated in postprocessing to reset statistics" << endl;
              Process::exit();
            }

          cout << "Lecture vitesse integree initiale dans fichier " << ns.fichier_reprise_vitesse_ << " timestep= " << timestep_reprise_integrated_velocity_ << endl;
          lire_dans_lata(ns.fichier_reprise_vitesse_, timestep_reprise_integrated_velocity_, geom_name, "INTEGRATED_VELOCITY", integrated_velocity_[0], integrated_velocity_[1],
                         integrated_velocity_[2]); // fonction qui lit un champ a partir d'un lata .
        }
      else
        {
          for (int i = 0; i < 3; i++)
            integrated_velocity_[i].data() = 0.;
          velocity_.valeur()[0].echange_espace_virtuel(velocity_.valeur()[0].ghost());
          velocity_.valeur()[1].echange_espace_virtuel(velocity_.valeur()[1].ghost());
          velocity_.valeur()[2].echange_espace_virtuel(velocity_.valeur()[2].ghost());

          update_integral_velocity(velocity_, integrated_velocity_, interfaces_->In(), integrated_timescale_);

        }
    }

  if ((( ns.coef_immobilisation_ > 1e-16) && is_stats_plans_activated()) || is_post_required("INTEGRATED_PRESSURE"))
    {
      integrated_pressure_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 0, "INTEGRATED_PRESSURE");
      champs_compris_.ajoute_champ(integrated_pressure_);
    }

  if (is_post_required("INTEGRATED_PRESSURE"))
    {
      if ((reprise) && (!reset_reprise_integrated_))
        {
          if (ns.fichier_reprise_vitesse_ == "??")
            {
              Cerr << "fichier_reprise_vitesse_ should be specified in the restart file in Navier_Stokes_FTD_IJK object" << endl;
              Process::exit();
            }
          const int timestep_reprise_integrated_pressure_ = 1;
          const Nom& geom_name = pressure_->get_domaine().le_nom();

          if (!lata_has_field(ns.fichier_reprise_vitesse_, timestep_reprise_integrated_pressure_, geom_name, "INTEGRATED_PRESSURE"))
            {
              Cerr << "fichier_reprise_vitesse_ " <<  ns.fichier_reprise_vitesse_ << " does not contain field INTEGRATED_PRESSURE to restart statistics. You may specify the flag reset_reprise_integrated in postprocessing to reset statistics" << endl;
              Process::exit();
            }
          cout << "Lecture pression integree initiale dans fichier " << ns.fichier_reprise_vitesse_ << " timestep= " << timestep_reprise_integrated_pressure_ << endl;
          lire_dans_lata(ns.fichier_reprise_vitesse_, timestep_reprise_integrated_pressure_, geom_name, "INTEGRATED_PRESSURE", integrated_pressure_); // fonction qui lit un champ a partir d'un lata .

        }
      else
        {
          integrated_pressure_.data() = 0.;

          // Le champ de pression initial ne vaut-il pas forcemment 0?
          update_integral_pressure(pressure_, integrated_pressure_, interfaces_->In(), integrated_timescale_);

        }
    }

  // Pour le post-traitement de lambda2  -  ALREADY DONE IN alloc_fields() !!!
//  if (liste_post_instantanes_.contient_("LAMBDA2"))
//    {
////      lambda2_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 0);
//      dudy_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 0);
//      dvdx_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 0);
//      dwdy_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 0);
//    }

  // Pour le check_stats_ :
  if (check_stats_)
    {
      Cout << "Initialisation champs analytiques (derivee P)" << "\ndPdx = " << expression_gradP_analytique_[0] << "\ndPdy = " << expression_gradP_analytique_[1] << "\ndPdz = "
           << expression_gradP_analytique_[2] << finl;

      Cout << "Initialisation champs analytiques (derivee U)" << "\ndUdx = " << expression_gradU_analytique_[0] << "\ndUdy = " << expression_gradU_analytique_[1] << "\ndUdz = "
           << expression_gradU_analytique_[2] << finl;

      Cout << "Initialisation champs analytiques (derivee V)" << "\ndVdx = " << expression_gradV_analytique_[0] << "\ndVdy = " << expression_gradV_analytique_[1] << "\ndVdz = "
           << expression_gradV_analytique_[2] << finl;

      Cout << "Initialisation champs analytiques (derivee W)" << "\ndWdx = " << expression_gradW_analytique_[0] << "\ndWdy = " << expression_gradW_analytique_[1] << "\ndWdz = "
           << expression_gradW_analytique_[2] << finl;

      Cout << "Initialisation champs analytiques (derivees secondes P) " << "\nddPdxdx = " << expression_grad2P_analytique_[0] << "\nddPdydy = " << expression_grad2P_analytique_[1] << "\nddPdzdz = "
           << expression_grad2P_analytique_[2];
      Cout << "\nddPdxdy = " << expression_grad2P_analytique_[3] << "\nddPdxdz = " << expression_grad2P_analytique_[4] << "\nddPdydz = " << expression_grad2P_analytique_[5] << finl;

      Cout << "Initialisation champs analytiques (derivees secondes U) " << "\nddUdxdx = " << expression_grad2U_analytique_[0] << "\nddUdydy = " << expression_grad2U_analytique_[1] << "\nddUdzdz = "
           << expression_grad2U_analytique_[2];
      Cout << "\nddUdxdy = " << expression_grad2U_analytique_[3] << "\nddUdxdz = " << expression_grad2U_analytique_[4] << "\nddUdydz = " << expression_grad2U_analytique_[5] << finl;

      Cout << "Initialisation champs analytiques (derivees secondes V) " << "\nddVdxdx = " << expression_grad2V_analytique_[0] << "\nddVdydy = " << expression_grad2V_analytique_[1] << "\nddVdzdz = "
           << expression_grad2V_analytique_[2];
      Cout << "\nddVdxdy = " << expression_grad2V_analytique_[3] << "\nddVdxdz = " << expression_grad2V_analytique_[4] << "\nddVdydz = " << expression_grad2V_analytique_[5] << finl;

      Cout << "Initialisation champs analytiques (derivees secondes W) " << "\nddWdxdx = " << expression_grad2W_analytique_[0] << "\nddWdydy = " << expression_grad2W_analytique_[1] << "\nddWdzdz = "
           << expression_grad2W_analytique_[2];
      Cout << "\nddWdxdy = " << expression_grad2W_analytique_[3] << "\nddWdxdz = " << expression_grad2W_analytique_[4] << "\nddWdydz = " << expression_grad2W_analytique_[5] << finl;

      // Le remplissage des set_field est fait par l'objet Statistiques au moment de l'appel a get_IJK_...
    }
}

void Postprocessing_IJK::fill_indic(bool reprise)
{
  // Meme if que pour l'allocation.
  // On ne fait le calcul/remplissage du champ que dans un deuxieme temps car on
  // n'avait pas les interfaces avant (lors de l'init)
  if (((ref_ijk_ft_->eq_ns().coef_immobilisation_ > 1e-16) && is_stats_plans_activated()) || is_post_required("INDICATRICE_PERTURBE")
      || (reprise && (fichier_reprise_indicatrice_non_perturbe_ != "??")))
    {
      init_indicatrice_non_perturbe();
    }
}

void Postprocessing_IJK::initialise_stats(Domaine_IJK& splitting, ArrOfDouble& vol_bulles, const double vol_bulle_monodisperse)
{
  cout << "Initialisation des statistiques. T_debut_statistiques=" << t_debut_statistiques_ << endl;
  statistiques_FT_.initialize(ref_ijk_ft_, splitting, check_stats_);
  // Si on utilise un seul groupe et qu'on impose un volume unique a toutes les bulles,
  if (vol_bulle_monodisperse >= 0.)
    {
      // on redimensionne le tableau a nb bulles reelles'
      vol_bulles.resize_array(interfaces_->get_nb_bulles_reelles());
      vol_bulles = vol_bulle_monodisperse;
    }
  // S'il n'y a pas qu'un group, on s'occupe des objets stats pour chaque group:
  const int nb_groups = interfaces_->nb_groups();
  if (nb_groups > 1)
    {
      groups_statistiques_FT_.dimensionner(nb_groups);
      for (int igroup = 0; igroup < nb_groups; igroup++)
        groups_statistiques_FT_[igroup].initialize(ref_ijk_ft_, splitting, check_stats_);
    }
}

void Postprocessing_IJK::init_indicatrice_non_perturbe()
{
  // Est-il deja rempli et stocke?
  // Si on n'est pas en reprise de calcul, le fichier "fichier_reprise_indicatrice_non_perturbe_" est forcement a "??"
  if ((fichier_reprise_indicatrice_non_perturbe_ != "??") && (!reset_reprise_integrated_))
    {
      const int timestep_reprise_indicatrice_non_perturbe = 1; // 1 ou 0 est le premier? attention au get_db ou latadb...
      cout << "Lecture indicatrice non perturbee dans fichier " << fichier_reprise_indicatrice_non_perturbe_ << " timestep= " << timestep_reprise_indicatrice_non_perturbe << endl;
      const Nom& geom_name = indicatrice_non_perturbe_.get_domaine().le_nom();
      lire_dans_lata(fichier_reprise_indicatrice_non_perturbe_, timestep_reprise_indicatrice_non_perturbe, geom_name, "INDICATRICE_PERTURBE", indicatrice_non_perturbe_); // fonction qui lit un champ a partir d'un lata .
    }
  else if (
    (( ref_ijk_ft_->eq_ns().coef_immobilisation_ > 1e-16) && is_stats_plans_activated())
    || is_post_required("INDICATRICE_PERTURBE")
  )
    {

      // Sinon, on le calcule une fois pour toute (cas bulles fixe = le champ ne varie pas en temps...)
      ArrOfDouble volume_reel;
      DoubleTab position;
      interfaces_->calculer_volume_bulles(volume_reel, position);
      interfaces_->compute_indicatrice_non_perturbe(indicatrice_non_perturbe_, ref_ijk_ft_->get_interface().I(), volume_reel, position);
      supprimer_chevauchement(indicatrice_non_perturbe_);
    }
}

void Postprocessing_IJK::posttraiter_champs_instantanes(const char *lata_name, double current_time, int time_iteration)
{
  throw; // THIS METHOD SHOULD NOT BE CALLED ANYMORE
}

// 2020.03.12. CHOIX : Meme en disable_diphasique, on fait appel a la classe fille stats FT
void Postprocessing_IJK::posttraiter_statistiques_plans(double current_time)
{
  statistiques().begin_count(postraitement_counter_);

  if (Process::je_suis_maitre())
    {
      Nom n("");
      // post-traitement des stats diphasiques :
      // (si calcul monophasique, on devrait avoir chi=1)
      for (int flag_valeur_instantanee = 0; flag_valeur_instantanee < 2; flag_valeur_instantanee++)
        {
          if (Option_IJK::DISABLE_DIPHASIQUE)
            n = "monophasique_";
          else
            n = "diphasique_";

          if (flag_valeur_instantanee == 0)
            n += "statistiques_";
          else
            n += "moyenne_spatiale_";

          n += Nom(current_time);
          if ((flag_valeur_instantanee) || (statistiques_FT_.t_integration() > 0.))
            {
              SFichier f(n + Nom(".txt"));
              f.setf(ios::scientific);      // precision pour allez chercher les 4 ordres
              f.precision(15);
              statistiques_FT_.postraiter(f, flag_valeur_instantanee /* flag pour ecrire la moyenne instantanee ou la moyenne */);
              statistiques_FT_.postraiter_thermique(current_time); /* moyenne instantanee et temporelle */

              // S'il n'y a pas qu'un group, on posttraite les objets stats pour chaque group:
              if ((!Option_IJK::DISABLE_DIPHASIQUE) && (interfaces_->nb_groups() > 1))
                {
                  for (int igroup = 0; igroup < interfaces_->nb_groups(); igroup++)
                    {
                      SFichier figroup(n + Nom("_grp") + Nom(igroup) + Nom(".txt"));
                      figroup.setf(ios::scientific);
                      figroup.precision(15);
                      groups_statistiques_FT_[igroup].postraiter(figroup, flag_valeur_instantanee /* flag pour ecrire la moyenne instantanee ou la moyenne */);
                      if (flag_valeur_instantanee == 1)
                        groups_statistiques_FT_[igroup].postraiter_thermique(current_time); /* moyenne instantanee et temporelle */
                    }
                }
            }
        }
      statistiques_FT_.postraiter_thermique(current_time); /* moyenne instantanee et temporelle */
    }
  statistiques().end_count(postraitement_counter_);

}

// Le nom du fichier est base sur le nom du cas...
// Si reset!=0, on efface le fichier avant d'ecrire, sinon on ajoute...
void Postprocessing_IJK::ecrire_statistiques_bulles(int reset, const Nom& nom_cas, const double current_time) const
{
  if (Option_IJK::DISABLE_DIPHASIQUE)
    return;

  statistiques().begin_count(postraitement_counter_);

  const DoubleTab& gravite = ref_ijk_ft_->milieu_ijk().gravite().valeurs();
  ArrOfDouble volume;
  DoubleTab position;
  ArrOfDouble surface;
  ArrOfDouble aspect_ratio;
  ArrOfDouble surfactant;
  ArrOfDouble surfactant_min;
  ArrOfDouble surfactant_max;
  const int nbulles = interfaces_->get_nb_bulles_reelles();
  DoubleTab hauteurs_bulles(nbulles, 3);
  DoubleTab bounding_box;
  interfaces_->calculer_bounding_box_bulles(bounding_box);
  for (int ib = 0; ib < nbulles; ib++)
    for (int dir = 0; dir < 3; dir++)
      hauteurs_bulles(ib, dir) = bounding_box(ib, dir, 1) - bounding_box(ib, dir, 0);

  // La methode calcule a present les surfaces meme pour les bulles ghost.
  // Pour les enlever, il suffit simplement de reduire la taille du tableau :
  interfaces_->calculer_surface_bulles(surface);
  surface.resize_array(nbulles);

  // La methode calcule a present les volumes meme pour les bulles ghost.
  // Pour les enlever, il suffit simplement de reduire la taille du tableau :
  interfaces_->calculer_volume_bulles(volume, position);
  interfaces_->calculer_aspect_ratio(aspect_ratio);
  interfaces_->calculer_surfactant(surfactant, surfactant_min, surfactant_max);
  volume.resize_array(nbulles);
  position.resize(nbulles, 3);

  DoubleTab poussee;
  interfaces_->calculer_poussee_bulles(gravite, poussee);

  if(ref_ijk_ft_->has_thermals())
    const_cast<Postprocessing_IJK*>(this)->thermals_->ecrire_statistiques_bulles(reset, nom_cas, current_time, surface);

  if (Process::je_suis_maitre())
    {
      char s[1000];
      const char *nomcas = nom_cas;
      SFichier fic;
      const int n = position.dimension(0);
      IOS_OPEN_MODE mode = (reset) ? ios::out : ios::app;

      auto write_func_tab = [&](const char * fnam, const DoubleTab& tab, int col)
      {
        snprintf(s, 1000, fnam, nomcas);
        fic.ouvrir(s, mode);
        snprintf(s, 1000, "%.16e ", current_time);
        fic << s;
        for (int i = 0; i < n; i++)
          {
            snprintf(s, 1000, "%.16e ", tab(i, col));
            fic << s;
          }
        fic << endl;
        fic.close();
      };

      auto write_func_arr = [&](const char * fnam, const ArrOfDouble& tab)
      {
        snprintf(s, 1000, fnam, nomcas);
        fic.ouvrir(s, mode);
        snprintf(s, 1000, "%.16e ", current_time);
        fic << s;
        for (int i = 0; i < n; i++)
          {
            snprintf(s, 1000, "%.16e ", tab[i]);
            fic << s;
          }
        fic << endl;
        fic.close();
      };

      write_func_tab("%s_bulles_pousseex.out", poussee, 0);
      write_func_tab("%s_bulles_hx.out", hauteurs_bulles,0);
      write_func_tab("%s_bulles_hy.out", hauteurs_bulles,1);
      write_func_tab("%s_bulles_hz.out", hauteurs_bulles,2);

      write_func_tab("%s_bulles_centre_x.out", position, 0);
      write_func_tab("%s_bulles_centre_y.out", position, 0);
      write_func_tab("%s_bulles_centre_z.out", position, 0);

      write_func_arr("%s_bulles_surface.out", surface);
      write_func_arr("%s_bulles_volume.out", volume);
      write_func_arr("%s_bulles_aspect_ratio.out", aspect_ratio);
      write_func_arr("%s_bulles_volume.out", volume);

      if (!interfaces_->maillage_ft_ijk().Surfactant_facettes().get_disable_surfactant())
        {
          write_func_arr("%s_bulles_surfactant.out", surfactant);
          write_func_arr("%s_bulles_surfactant_min.out", surfactant_min);
          write_func_arr("%s_bulles_surfactant_max.out", surfactant_max);
        }
      if (interfaces_->follow_colors())
        {
          const ArrOfInt& colors = interfaces_->get_colors();
          snprintf(s, 1000, "%s_bulles_colors.out", nomcas);
          fic.ouvrir(s, mode);
          snprintf(s, 1000, "%.16e ", current_time);
          fic << s;
          for (int i = 0; i < n; i++)
            {
              snprintf(s, 1000, "%d ", (True_int) colors[i]);
              fic << s;
            }
          fic << endl;
          fic.close();
        }
    }
  statistiques().end_count(postraitement_counter_);

}

void Postprocessing_IJK::ecrire_statistiques_cisaillement(int reset, const Nom& nom_cas, const double current_time) const
{
  if (Option_IJK::DISABLE_DIPHASIQUE)
    return;

  statistiques().begin_count(postraitement_counter_);

  double v_x_droite;
  double v_y_droite;
  double v_z_droite;

  double v_x_gauche;
  double v_y_gauche;
  double v_z_gauche;

  Navier_Stokes_FTD_IJK& ns = const_cast<Postprocessing_IJK*>(this)->ref_ijk_ft_->eq_ns();
  ns.calculer_vitesse_gauche(velocity_.valeur()[0],velocity_.valeur()[1],velocity_.valeur()[2],v_x_gauche,v_y_gauche,v_z_gauche);
  ns.calculer_vitesse_droite(velocity_.valeur()[0],velocity_.valeur()[1],velocity_.valeur()[2],v_x_droite,v_y_droite,v_z_droite);

  if (Process::je_suis_maitre())
    {
      char s[1000];
      const char *nomcas = nom_cas;
      SFichier fic;
      IOS_OPEN_MODE mode = (reset) ? ios::out : ios::app;

      snprintf(s, 1000, "%s_cisaillement.out", nomcas);
      // Cerr << "Ecriture des donnees par bulles: fichier " << s << endl;
      fic.ouvrir(s, mode);

      for (const auto v: {current_time, v_x_droite, v_y_droite, v_z_droite, v_x_gauche, v_y_gauche, v_z_gauche})
      {
        snprintf(s, 1000, "%.16e ", v);
        fic << s;
      }

      fic << endl;
      fic.close();
    }
  statistiques().end_count(postraitement_counter_);

}

void Postprocessing_IJK::ecrire_statistiques_rmf(int reset, const Nom& nom_cas, const double current_time) const
{
  if (Option_IJK::DISABLE_DIPHASIQUE)
    return;

  statistiques().begin_count(postraitement_counter_);

  double ax_PID, ay_PID, az_PID;

  const_cast<Postprocessing_IJK*>(this)->ref_ijk_ft_->eq_ns().calculer_terme_asservissement(ax_PID,ay_PID,az_PID);

  if (Process::je_suis_maitre())
    {
      char s[1000];
      const char *nomcas = nom_cas;
      SFichier fic;
      IOS_OPEN_MODE mode = (reset) ? ios::out : ios::app;

      snprintf(s, 1000, "%s_rmf.out", nomcas);
      // Cerr << "Ecriture des donnees par bulles: fichier " << s << endl;
      fic.ouvrir(s, mode);
      for (const auto v: {current_time, ax_PID, ay_PID, az_PID})
      {
        snprintf(s, 1000, "%.16e ", v);
        fic << s;
      }
      fic << endl;
      fic.close();
    }
  statistiques().end_count(postraitement_counter_);
}

/** Methode qui met a jour l'indicatrice, les termes de repulsion
 * ainsi que les termes interfaciaux : ai, kappa*ai, n(aux cellules)
 *
 * Par definition, mettre igroup a -1 pour inclure toutes les bulles
 * Dans ce cas, la methode met a jour l'ev de l'indicatrice au lieu de celui de interfaces_.groups_indicatrice_n_ns()[igroup]
 *
 * Attention: de nombreux tableaux sont modifies par cette methode en sortie.
 * Ils peuvent etre des tableaux de travail. Si on veut qu'il soient correctent
 * pour la suite, il faut faire l'appel avec les champs globaux (incluant tous
 * les groupes a la fin). Sinon, les champs en ai, normale ou grad_I ne contiendront qu'un groupe.
 */
void Postprocessing_IJK::update_stat_ft(const double dt)
{
  Navier_Stokes_FTD_IJK& ns = ref_ijk_ft_->eq_ns();
  //ArrOfDouble volume;
  //DoubleTab position;
  //interfaces_.calculer_volume_bulles(volume, position);
  static Stat_Counter_Id updtstat_counter_ = statistiques().new_counter(2, "update statistiques");
  statistiques().begin_count(updtstat_counter_);
  if (Option_IJK::DISABLE_DIPHASIQUE)
    {
      // Calcul du champ grad_P_ est fait dans update_stat
      statistiques_FT_.update_stat(ref_ijk_ft_, dt);
      return;
    }
  int nb_groups = interfaces_->nb_groups();
  // Boucle debute a -1 pour faire l'indicatrice globale.
  // S'il n'y a pas de groupes de bulles (monophasique ou monodisperse), on passe exactement une fois dans la boucle
  if (nb_groups == 1)
    nb_groups = 0; // Quand il n'y a qu'un groupe, on ne posttraite pas les choses pour ce groupe unique puisque c'est identique au cas global
  for (int igroup = -1; igroup < nb_groups; igroup++)
    {
      if (is_post_required("AIRE_INTERF"))
        {
          IJK_Field_double& ai_ft = scalar_post_fields_.at("AIRE_INTERF");
          interfaces_->calculer_normales_et_aires_interfaciales(ai_ft, kappa_ai_ft_, normale_cell_ft_, igroup);
          // Puis les redistribue sur le ns :
          ns.redistribute_from_splitting_ft_elem_.redistribute(ai_ft, ai_ns_);
          ns.redistribute_from_splitting_ft_elem_.redistribute(kappa_ai_ft_, kappa_ai_ns_);
          ns.redistribute_from_splitting_ft_elem_.redistribute(normale_cell_ft_, normale_cell_ns_);
        }
      // (pas besoin d'echange EV car ils n'ont pas de ghost).

      // Calcul du gradient de l'indicatrice et de vitesse :
      if (igroup == -1)
        {
          // interfaces_.In().echange_espace_virtuel(1);
          // Calcul des champs grad_I_ns_
          calculer_gradient_indicatrice(interfaces_->In());
          ns.transfer_ft_to_ns(); // pour remplir : terme_repulsion_interfaces_ft_ et terme_abs_repulsion_interfaces_ft_
          // Calcul des champs grad_I_ns_, terme_repulsion_interfaces_ns_, terme_abs_repulsion_interfaces_ns_
          // a partir de interfaces_.In(), et terme_*_ft_
          statistiques_FT_.update_stat(ref_ijk_ft_, dt);
        }
      else
        {
          // interfaces_.groups_indicatrice_n_ns()[igroup].echange_espace_virtuel(1);
          // Calcul des champs grad_I_ns_, terme_repulsion_interfaces_ns_, terme_abs_repulsion_interfaces_ns_
          // a partir de interfaces_.In(), et terme_*_ft_
          calculer_gradient_indicatrice(interfaces_->groups_indicatrice_n_ns()[igroup]);
          ns.transfer_ft_to_ns();
          groups_statistiques_FT_[igroup].update_stat(ref_ijk_ft_, dt);
        }
    }
  statistiques().end_count(updtstat_counter_);
}

// Calcul du lambda2 a partir du gradient.
// A optimiser simplement en mutualisant avec la methode d'update_stats.
// Et en ne faisant le calcul que si besoin, cad si les champs de gradient ne sont pas a jour...
void Postprocessing_IJK::update_gradU_lambda2(const bool need_lambda2)
{
  IJK_Field_double& lambda2 = scalar_post_fields_.at("LAMBDA2");
  // TODO : Clean theses : dudx_  and vectorise with what's in statistiques ... gradU[2]
  if (need_lambda2)
    compute_and_store_gradU_cell(velocity_.valeur()[0], velocity_.valeur()[1], velocity_.valeur()[2],
                                 /* Et les champs en sortie */
                                 dudx_, dvdy_, dwdx_, dudz_, dvdz_, dwdz_, 1 /* yes compute_all */,
                                 dudy_, dvdx_, dwdy_, lambda2);
  statistiques_FT_.compute_and_store_gradU_cell(velocity_.valeur()[0], velocity_.valeur()[1], velocity_.valeur()[2]);
}

void Postprocessing_IJK::get_update_lambda2_and_rot_and_Q()
{
  IJK_Field_vector3_double& rot = vect_post_fields_.at("CURL");
  IJK_Field_double& critere_Q = scalar_post_fields_.at("CRITERE_Q");
  //IJK_Field_double& lambda2 = scalar_post_fields_.at("LAMBDA2");

  const IJK_Field_vector3_double& gradU=statistiques_FT_.get_IJK_field_vector("dUd");
  const IJK_Field_vector3_double& gradV=statistiques_FT_.get_IJK_field_vector("dVd");
  const IJK_Field_vector3_double& gradW=statistiques_FT_.get_IJK_field_vector("dWd");
  const IJK_Field_double& dudx = gradU[0];
  const IJK_Field_double& dudy = gradU[1];
  const IJK_Field_double& dudz = gradU[2];
  const IJK_Field_double& dvdx = gradV[0];
  const IJK_Field_double& dvdy = gradV[1];
  const IJK_Field_double& dvdz = gradV[2];
  const IJK_Field_double& dwdx = gradW[0];
  const IJK_Field_double& dwdy = gradW[1];
  const IJK_Field_double& dwdz = gradW[2];
  update_gradU_lambda2();
  // Nombre local de mailles en K
  const int kmax = critere_Q.nk();
  const int imax = critere_Q.ni();
  const int jmax = critere_Q.nj();
  for (int k = 0; k < kmax; k++)
    for (int j = 0; j < jmax; j++)
      for (int i = 0; i < imax; i++)
        {
          rot[0](i, j, k) = dwdy(i, j, k) - dvdz(i, j, k);
          rot[1](i, j, k) = dudz(i, j, k) - dwdx(i, j, k);
          rot[2](i, j, k) = dvdx(i, j, k) - dudy(i, j, k);
          // Calcul du critere Q selon (Jeong & Hussain 1995)
          critere_Q(i, j, k) = -0.5
                               * (dudx(i, j, k) * dudx(i, j, k)
                                  + 2. * dudy(i, j, k) * dvdx(i, j, k)
                                  + 2. * dudz(i, j, k) * dwdx(i, j, k)
                                  + dvdy(i, j, k) * dvdy(i, j, k)
                                  + 2. * dvdz(i, j, k) * dwdy(i, j, k)
                                  + dwdz(i, j, k) * dwdz(i, j, k));
        }
}

/** Was the field of name 'nom' requested for postprocessing?
 */
bool Postprocessing_IJK::is_post_required(const Motcle& nom) const
{
  return (std::find(list_post_required_.begin(), list_post_required_.end(), nom) != list_post_required_.end());
}

// Pour qu'un champ vectoriel puisse etre interpole a l'element, il faut qu'il ait au moins 1 ghost.
void Postprocessing_IJK::Fill_postprocessable_fields(std::vector<FieldInfo_t>& chps)
{
  std::vector<FieldInfo_t> c =
  {
    // Name     /     Localisation (elem, face, ...) /    Nature (scalare, vector)   /  Located on interface?

    { "FORCE_PH", Entity::FACE, Nature_du_champ::vectoriel, false },
    { "CURL", Entity::ELEMENT, Nature_du_champ::vectoriel, false },
    { "CRITERE_Q", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "NUM_COMPO", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "FORCE_PH", Entity::FACE, Nature_du_champ::vectoriel, false },
    { "COORDS", Entity::FACE, Nature_du_champ::vectoriel, false },
    { "LAMBDA2", Entity::ELEMENT, Nature_du_champ::scalaire, false },

    { "VELOCITY_ANA", Entity::FACE, Nature_du_champ::vectoriel, false },
    { "ECART_ANA", Entity::FACE, Nature_du_champ::vectoriel, false },
    { "PRESSURE_ANA", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "ECART_P_ANA", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "D_VELOCITY_ANA", Entity::FACE, Nature_du_champ::vectoriel, false },

    { "D_VELOCITY", Entity::FACE, Nature_du_champ::vectoriel, false },
    { "OP_CONV", Entity::FACE, Nature_du_champ::vectoriel, false },
    { "D_PRESSURE", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "MU", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "RHO", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "INDICATRICE_NS", Entity::ELEMENT, Nature_du_champ::scalaire, false },

    // A faire sauter avec le travail de William:
    { "INDICATRICE_FT", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "GRAD_INDICATRICE_FT", Entity::FACE, Nature_du_champ::vectoriel, false },
    { "REBUILT_INDICATRICE_FT", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "GROUPS_FT", Entity::ELEMENT, Nature_du_champ::vectoriel, false },

    // Dans NS :
    // { "SOURCE_QDM_INTERF", Entity::FACE, Nature_du_champ::vectoriel, false },
    // { "SHIELD_REPULSION", Entity::FACE, Nature_du_champ::vectoriel, false },
    { "RHO_SOURCE_QDM_INTERF", Entity::FACE, Nature_du_champ::vectoriel, false },
    { "RHO_SOURCE_QDM_INTERF", Entity::ELEMENT, Nature_du_champ::vectoriel, false },
    { "BK_SOURCE_QDM_INTERF", Entity::FACE, Nature_du_champ::vectoriel, false },

    { "AIRE_INTERF", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "COURBURE_AIRE_INTERF", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "NORMALE_EULER", Entity::ELEMENT, Nature_du_champ::vectoriel, false },
    { "PRESSURE_LIQ", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "PRESSURE_VAP", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "GROUPS", Entity::ELEMENT, Nature_du_champ::vectoriel, false },
    { "SURFACE_VAPEUR_PAR_FACE", Entity::FACE, Nature_du_champ::vectoriel, false },
    { "BARYCENTRE_VAPEUR_PAR_FACE", Entity::FACE, Nature_du_champ::vectoriel, false },

    { "VARIABLE_SOURCE", Entity::FACE, Nature_du_champ::vectoriel, false },
    { "INTEGRATED_VELOCITY", Entity::FACE, Nature_du_champ::vectoriel, false },
    { "INTEGRATED_PRESSURE", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "INDICATRICE_PERTURBE", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "INTEGRATED_TIMESCALE", Entity::ELEMENT, Nature_du_champ::scalaire, false }
  };

  chps.insert(chps.end(), c.begin(), c.end());
}

void Postprocessing_IJK::get_noms_champs_postraitables(Noms& noms,Option opt) const
{
  statistiques_FT_.get_noms_champs_postraitables(noms, opt);
  for (const auto& n : champs_compris_.liste_noms_compris())
    noms.add(n);
  for (const auto& n : champs_compris_.liste_noms_compris_vectoriel())
    noms.add(n);
}
bool Postprocessing_IJK::has_champ(const Motcle& nom) const
{
  if (statistiques_FT_.has_champ(nom))
    {
      return true;
    }
  return champs_compris_.has_champ(nom);
}
bool Postprocessing_IJK::has_champ_vectoriel(const Motcle& nom) const
{
  if (statistiques_FT_.has_champ_vectoriel(nom))
    {
      return true;
    }
  return champs_compris_.has_champ_vectoriel(nom);
}
/** Retrieve requested field for postprocessing, potentially updating it.
 */
const IJK_Field_double& Postprocessing_IJK::get_IJK_field(const Motcle& nom)
{

  if (statistiques_FT_.has_champ(nom))
    {
      return statistiques_FT_.get_IJK_field(nom);
    }

  if (!has_champ(nom))
    {
      Cerr << "ERROR in Postprocessing_IJK::get_IJK_field : " << finl;
      Cerr << "Requested field '" << nom << "' is not recognized by Postprocessing_IJK::get_IJK_field()." << finl;
      throw;
    }

  double current_time = ref_ijk_ft_->schema_temps_ijk().get_current_time();

  // TODO pas optimal :
  if (nom == "LAMBDA2")
    update_gradU_lambda2();
  if (nom == "CRITERE_Q")
    get_update_lambda2_and_rot_and_Q();

  if (nom == "NUM_COMPO")
    {
      IJK_Field_double& num_compo_ft = scalar_post_fields_.at("NUM_COMPO");
      const int ni = num_compo_ft.ni(), nj = num_compo_ft.nj(), nk = num_compo_ft.nk();
      const IntVect& num_compo = interfaces_->get_num_compo();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const int num_elem = domaine_ft_->convert_ijk_cell_to_packed(i, j, k);
              num_compo_ft(i, j, k) = num_compo[num_elem];
            }
    }

  if (nom == "AIRE_INTERF")
    {
      IJK_Field_double& ai_ft = scalar_post_fields_.at("AIRE_INTERF");
      interfaces_->calculer_aire_interfaciale(ai_ft);
    }
  if (nom == "PRESSURE_ANA")
    {
      set_field_data(scalar_post_fields_.at("PRESSURE_ANA"), expression_pression_analytique_, current_time);
    }
  if (nom == "ECART_P_ANA")
    {
      auto& pressure_ana=scalar_post_fields_.at("PRESSURE_ANA");
      double ct = current_time;
      if ( sub_type(Schema_Euler_explicite_IJK, ref_ijk_ft_->schema_temps_ijk()) )
        {
          ct -= ref_ijk_ft_->schema_temps_ijk().get_timestep();
        }
      else if ( sub_type(Schema_RK3_IJK, ref_ijk_ft_->schema_temps_ijk()) )
        {
          Schema_RK3_IJK& rk3 = ref_cast(Schema_RK3_IJK, ref_ijk_ft_->schema_temps_ijk());
          Cerr << "rkstep " << rk3.get_rk_step() << finl;
          int rk_step_before = rk3.get_rk_step();
          if ((rk_step_before == 0) || (rk_step_before == 3))
            rk_step_before = 2;
          else if (rk_step_before == 1)
            rk_step_before = 0;
          else
            /* ici, c'est rk_step_before=2 */
            rk_step_before = 1;
          Cerr << "rkstep_before " << rk_step_before << finl;
          const double intermediate_dt = compute_fractionnal_timestep_rk3(rk3.get_timestep(), rk_step_before);
          ct -= intermediate_dt;
        }
      else
        {
          Cerr << "To do for other time scheme" << endl;
        }
      Cerr << "GB: ERROR P FIELD " << ct;
      double err = 0.;
      set_field_data(pressure_ana, expression_pression_analytique_, ct);
      const int ni = pressure_->ni();
      const int nj = pressure_->nj();
      const int nk = pressure_->nk();
      const trustIdType ntot = Process::mp_sum(ni * nj * nk);
      // La pression est definie a une constante pres:
      const double cst_press = pressure_ana(0, 0, 0) - pressure_.valeur()(0, 0, 0);
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double val = pressure_ana(i, j, k) - pressure_.valeur()(i, j, k) - cst_press;
              ecart_p_ana_(i, j, k) = val;
              err += val * val;
            }
      err = Process::mp_sum(err);
      err = sqrt(err / static_cast<double>(ntot));
      Cerr << " " << err;
      if (!Process::je_suis_maitre())
        {
          Process::Journal() << "Postprocessing_IJK::posttraiter_champs_instantanes : OWN_PTR(Champ_base) ECART_P_ANA sur ce proc (ni,nj,nk,ntot):" << " " << ni << " " << nj << " " << nk << " " << ntot << finl;
        }
      ecart_p_ana_.echange_espace_virtuel(ecart_p_ana_.ghost());
      Cerr << finl;
    }

  return champs_compris_.get_champ(nom);

  //  if (nom.debute_par("dU"))
  //    {
  //      const IJK_Field_vector3_double& gradU = statistiques_FT_.get_IJK_field_vector("gradU");
  //      if (nom == "DUDX")
  //        return gradU[0];
  //      if (nom == "DUDY")
  //        return gradU[1];
  //      if (nom == "DUDZ")
  //        return gradU[2];
  //    }
  //  if (nom.debute_par("dV"))
  //    {
  //      const IJK_Field_vector3_double& gradV = statistiques_FT_.get_IJK_field_vector("gradV");
  //      if (nom == "DVDX")
  //        return gradV[0];
  //      if (nom == "DVDY")
  //        return gradV[1];
  //      if (nom == "DVDZ")
  //        return gradV[2];
  //    }
  //  if (nom.debute_par("dW"))
  //    {
  //      const IJK_Field_vector3_double& gradW = statistiques_FT_.get_IJK_field_vector("gradW");
  //      if (nom == "DWDX")
  //        return gradW[0];
  //      if (nom == "DWDY")
  //        return gradW[1];
  //      if (nom == "DWDZ")
  //        return gradW[2];
  //    }
  //
  //  // GAB, sondes THI
  //  if (nom.debute_par("FORCE_PH"))
  //    {
  //      // A priori inutile : tester si sonde ok avec champ_a_postrer sans FORCE_PH
  //      // Reponse GB : le vrai test a faire c'est si le field force_ph existe, ie s'il y a un forcage
  //      if (!liste_post_instantanes_.contient_("FORCE_PH"))
  //        {
  //          Cerr << "A probe is attempting to access a field FORCE_PH while it has not been computed in the post-processed fields" << endl;
  //          Process::exit();
  //        }
  ////      IJK_Field_vector3_double& source_spectrale = ref_ijk_ft_.forcage_.get_force_ph2();
  //      if (nom == "FORCE_PH_X")
  //        return source_spectrale_.valeur()[0];
  //      if (nom == "FORCE_PH_Y")
  //        return source_spectrale_.valeur()[1];
  //      if (nom == "FORCE_PH_Z")
  //        return source_spectrale_.valeur()[2];
  //    }
  //  //
  //  // if (Option_IJK::DISABLE_DIPHASIQUE)
  //  {
  //    if (nom == "VELOCITY_ANA_X")
  //      return velocity_ana_[0];
  //    if (nom == "VELOCITY_ANA_Y")
  //      return velocity_ana_[1];
  //    if (nom == "VELOCITY_ANA_Z")
  //      return velocity_ana_[2];
  //    if (nom == "ECART_ANA_X")
  //      return ecart_ana_[0];
  //    if (nom == "ECART_ANA_Y")
  //      return ecart_ana_[1];
  //    if (nom == "ECART_ANA_Z")
  //      return ecart_ana_[2];
  //  }
  //
  //  const int idx_wanted = convert_suffix_to_int(nom);
  //  const Nom field_name = nom.getPrefix(Nom("_") + Nom(idx_wanted));
  //  Cerr << "In get_IJK_field by name : " << nom << " read as : (" << field_name << " ; " << idx_wanted << ")" << finl;
  //
  //// Remplir la liste de tous les possibles :
  //  Motcles liste_champs_thermiques_possibles;
  //  posttraiter_tous_champs_thermique(liste_champs_thermiques_possibles, 0);
  //  int rang = liste_champs_thermiques_possibles.rang(field_name);
  //  if (rang == -1)
  //    {
  //      Cerr << field_name << " not found as possible for field name. Should be in the list: " << liste_champs_thermiques_possibles << finl;
  //      Process::exit();
  //    }
  //
  //  const Motcle& mot = liste_champs_thermiques_possibles[rang];
  //  if ((mot == field_name) && (idx_wanted >= 0))
  //    {
  //      Cerr << "found as planned " << endl;
  //    }
  //  else
  //    {
  //      Cerr << "Some issue with the name provided for the sonde. Unrecognised." << finl;
  //      Process::exit();
  //    }
  //
  //  Cerr << "Erreur dans Postprocessing_IJK::get_IJK_field : " << endl;
  //  Cerr << "Champ demande : " << nom << endl;
  //  Cerr << "Liste des champs possibles pour la thermique : " << liste_champs_thermiques_possibles << finl;
  //  Process::exit();
  //  throw;


}

const IJK_Field_vector3_double& Postprocessing_IJK::get_IJK_field_vector(const Motcle& nom)
{
  if (nom == "CURL")
    get_update_lambda2_and_rot_and_Q();

  if (statistiques_FT_.has_champ_vectoriel(nom))
    {
      return statistiques_FT_.get_IJK_field_vector(nom);
    }
  if (!has_champ_vectoriel(nom))
    {
      Cerr << "ERROR in Postprocessing_IJK::get_IJK_field_vector : " << finl;
      Cerr << "Requested field '" << nom << "' is not recognized by Postprocessing_IJK::get_IJK_field_vector()." << finl;
      throw;
    }

  double current_time = ref_ijk_ft_->schema_temps_ijk().get_current_time();

  if (nom == "VELOCITY_ANA")
    {
      for (int i = 0; i < 3; i++)
        set_field_data(velocity_ana_[i], expression_vitesse_analytique_[i], current_time);
    }
  if (nom == "ECART_ANA")
    {
      Cerr << "GB: ERROR FIELD " << current_time;
      for (int dir = 0; dir < 3; dir++)
        {
          double err = 0.;
          set_field_data(velocity_ana_[dir], expression_vitesse_analytique_[dir], current_time);
          const int ni = ref_ijk_ft_->eq_ns().velocity_[dir].ni();
          const int nj = ref_ijk_ft_->eq_ns().velocity_[dir].nj();
          const int nk = ref_ijk_ft_->eq_ns().velocity_[dir].nk();
          const double ntot = Process::mp_sum_as_double(ni * nj * nk);
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  const double val = velocity_ana_[dir](i, j, k) - ref_ijk_ft_->eq_ns().velocity_[dir](i, j, k);
                  ecart_ana_[dir](i, j, k) = val;
                  err += val * val;
                }
          err = Process::mp_sum(err);
          err = sqrt(err / ntot);
          Cerr << " " << err;
          if (!Process::je_suis_maitre())
            {
              Process::Journal() << "Postprocessing_IJK::posttraiter_champs_instantanes : OWN_PTR(Champ_base) ECART_ANA sur ce proc (ni,nj,nk,ntot):" << " " << ni << " " << nj << " " << nk << " " << ntot << finl;
            }
        }
      Cerr << finl;
    }

  if (nom == "INTEGRATED_VELOCITY")
    {
      update_integral_velocity(velocity_, integrated_velocity_, interfaces_->In(), integrated_timescale_);
    }
  if (nom == "INTEGRATED_PRESSURE")
    {
      update_integral_pressure(pressure_, integrated_pressure_, interfaces_->In(), integrated_timescale_);
    }
  if (nom == "INTEGRATED_TIMESCALE")
    {
      update_integral_indicatrice(interfaces_->In(), 1. /* Should be the integration timestep */, integrated_timescale_);
    }
  if (nom == "INDICATRICE_PERTURBE")
    {
      //Faut-il faire un update_...( indicatrice_non_perturbe_) avant?
    }



  return champs_compris_.get_champ_vectoriel(nom);
}

const int& Postprocessing_IJK::get_IJK_flag(const Nom& nom) const
{
  Cerr << "Erreur dans Postprocessing_IJK::get_IJK_variable : " << "Variable demandee : " << nom << " Liste des variables possibles : " << finl;
  Process::exit();
  throw;
}

void Postprocessing_IJK::sauvegarder_post(const Nom& lata_name)
{
  if (is_post_required("INTEGRATED_VELOCITY"))
    dumplata_vector(lata_name, "INTEGRATED_VELOCITY", integrated_velocity_[0], integrated_velocity_[1], integrated_velocity_[2], 0);

  if (is_post_required("INTEGRATED_PRESSURE"))
    dumplata_scalar(lata_name, "INTEGRATED_PRESSURE", integrated_pressure_, 0);

  if (is_post_required("INDICATRICE_PERTURBE"))
    dumplata_scalar(lata_name, "INDICATRICE_PERTURBE", indicatrice_non_perturbe_, 0);

  if (is_post_required("INTEGRATED_TIMESCALE"))
    dumplata_scalar(lata_name, "INTEGRATED_TIMESCALE", integrated_timescale_, 0);
}

void Postprocessing_IJK::sauvegarder_post_maitre(const Nom& lata_name, SFichier& fichier) const
{

  if (statistiques_FT_.t_integration() > 0.)
    {
      Cerr << "All bubbles : " << endl;
      fichier << " statistiques_FT " << statistiques_FT_;
      // S'il y a plusieurs groups, on s'occupe des objets stats pour chaque group:
      // (en ecrivant directement le vecteur d'objets)
      if (interfaces_->nb_groups() > 1)
        {
          Cerr << "Group by group :" << endl;
          fichier << " groups_statistiques_FT " << groups_statistiques_FT_;
        }
    }
}

void Postprocessing_IJK::reprendre_post(Param& param)
{
  param.ajouter("statistiques_FT", &statistiques_FT_);
  param.ajouter("groups_statistiques_FT", &groups_statistiques_FT_);
}


void Postprocessing_IJK::fill_op_conv()
{
  if (liste_post_instantanes_.contient_("OP_CONV"))
    for (int i = 0; i < 3; i++)
      op_conv_[i].data() = d_velocity_.valeur()[i].data();

  if (liste_post_instantanes_.contient_("CELL_OP_CONV"))
    interpolate_to_center(cell_op_conv_,d_velocity_);
}

void Postprocessing_IJK::fill_surface_force(IJK_Field_vector3_double& the_field_you_know)
{
  double volume = 1.;
  for (int i = 0; i < 3; i++)
    volume *= domaine_ijk_->get_constant_delta(i);

  if (is_post_required("RHO_SOURCE_QDM_INTERF"))
    {
      IJK_Field_vector3_double& rho_Ssigma = vect_post_fields_.at("RHO_SOURCE_QDM_INTERF");
      for (int dir = 0; dir < 3; dir++)
        {
          IJK_Field_double& source = ref_ijk_ft_->eq_ns().terme_source_interfaces_ns_[dir];
          for (int k = 0; k < source.nk(); k++)
            for (int j = 0; j < source.nj(); j++)
              for (int i = 0; i < source.ni(); i++)
                rho_Ssigma[dir](i,j,k) = source(i,j,k)/volume;
        }
    }

  // TODO : probablement inutile nouvelle syntaxe
  if (liste_post_instantanes_.contient_("CELL_RHO_SOURCE_QDM_INTERF"))
    {
      interpolate_to_center(cell_rho_Ssigma_,ref_ijk_ft_->eq_ns().terme_source_interfaces_ns_);
      for (int dir = 0; dir < 3; dir++)
        {
          IJK_Field_double& source = cell_rho_Ssigma_[dir];
          for (int k = 0; k < source.nk(); k++)
            for (int j = 0; j < source.nj(); j++)
              for (int i = 0; i < source.ni(); i++)
                cell_rho_Ssigma_[dir](i,j,k) = source(i,j,k)/volume;
        }
    }
}

// Calcul du gradient de l'indicatrice et de pression :
//   Attention, il faut que la pression et l'indicatrice soient a jour
//   dans leur espaces virtuels avant d'appeler cette methode
// Methode qui calcule des champs grad_I_ns_,
// a partir de pressure_ et indicatrice_ns__
void Postprocessing_IJK::calculer_gradient_indicatrice(const IJK_Field_double& indic)
{
  // Remise a zero :
  for (int dir = 0; dir < 3; dir++)
    {
      grad_I_ns_[dir].data() = 0.;
    }

  // From IJK_Navier_Stokes_Tools.cpp
  // interfaces_.In().echange_espace_virtuel(1);
  add_gradient_times_constant(indic, 1. /*Constante multiplicative*/, grad_I_ns_[DIRECTION_I], grad_I_ns_[DIRECTION_J], grad_I_ns_[DIRECTION_K]);

  for (int dir = 0; dir < 3; dir++)
    {
      grad_I_ns_[dir].echange_espace_virtuel(1);
    }
}

void Postprocessing_IJK::alloc_fields()
{
  rebuilt_indic_.allocate(domaine_ft_, Domaine_IJK::ELEM, 0);
  if (is_post_required("AIRE_INTERF") && (!Option_IJK::DISABLE_DIPHASIQUE || is_stats_plans_activated()))
    {
      IJK_Field_double& ai_ft = scalar_post_fields_.at("AIRE_INTERF");
      ai_ft.allocate(domaine_ft_, Domaine_IJK::ELEM, 0, "AIRE_INTERF");
      champs_compris_.ajoute_champ(ai_ft);
    }
  // Pour les stats, on calcule kappa*ai :
  if (!Option_IJK::DISABLE_DIPHASIQUE && is_stats_plans_activated())
    {
      kappa_ai_ft_.allocate(domaine_ft_, Domaine_IJK::ELEM, 0);
      kappa_ai_ns_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 0);
      ai_ns_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 0);
      allocate_cell_vector(normale_cell_ns_, domaine_ijk_, 0);
    }

  // For the pressure field extension:
  if (!Option_IJK::DISABLE_DIPHASIQUE
      && ((is_post_required("PRESSURE_LIQ")) || (is_post_required("PRESSURE_VAP")) || is_stats_plans_activated()))
    {
      extended_pressure_computed_ = 1;
      pressure_ft_.allocate(domaine_ft_, Domaine_IJK::ELEM, 5);
      extended_pl_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 0, "PRESSURE_LIQ");
      champs_compris_.ajoute_champ(extended_pl_);
      extended_pv_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 0, "PRESSURE_VAP");
      champs_compris_.ajoute_champ(extended_pv_);
      extended_pl_ft_.allocate(domaine_ft_, Domaine_IJK::ELEM, 0);
      extended_pv_ft_.allocate(domaine_ft_, Domaine_IJK::ELEM, 0);
    }
  else
    extended_pressure_computed_ = 0;

  // Allocation du champ de normale aux cellules :
  if (!Option_IJK::DISABLE_DIPHASIQUE && ((is_post_required("NORMALE_INTERF")) || (is_post_required("PRESSURE_LIQ")) // Je ne suis pas sur que ce soit necessaire. Seulement si on l'utilise dans le calcul de p_ext
                                          || (is_post_required("PRESSURE_VAP")) // Je ne suis pas sur que ce soit necessaire.
                                          || is_stats_plans_activated()))
    {
      allocate_cell_vector(normale_cell_ft_, domaine_ft_, 0);
    }


  if (is_post_required("NUM_COMPO"))
    {
      IJK_Field_double& num_c = scalar_post_fields_.at("NUM_COMPO");
      num_c.allocate(domaine_ft_, Domaine_IJK::ELEM, 0, "NUM_COMPO");
      champs_compris_.ajoute_champ(num_c);
    }

  if (is_post_required("VELOCITY_ANA") or is_post_required("ECART_ANA") )
    {
      allocate_velocity(velocity_ana_, domaine_ijk_, 1, "VELOCITY_ANA");
      champs_compris_.ajoute_champ_vectoriel(velocity_ana_);
    }

  if (is_post_required("ECART_ANA"))
    {
      allocate_velocity(ecart_ana_, domaine_ijk_, 1, "ECART_ANA");
      champs_compris_.ajoute_champ_vectoriel(ecart_ana_);
    }

  if (is_post_required("PRESSURE_ANA") or is_post_required("ECART_P_ANA") )
    {
      scalar_post_fields_.at("PRESSURE_ANA").allocate(domaine_ijk_, Domaine_IJK::ELEM, 1, "PRESSURE_ANA");
      champs_compris_.ajoute_champ(scalar_post_fields_.at("PRESSURE_ANA"));
    }

  if (is_post_required("ECART_P_ANA"))
    {
      ecart_p_ana_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 1, "ECART_P_ANA");
      champs_compris_.ajoute_champ(ecart_p_ana_);
    }



  // Allocation des champs derivee de vitesse :
  if (is_post_required("LAMBDA2")|| is_post_required("CRITERE_Q") || is_post_required("CURL"))
    {
      dudx_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 1);
      dudy_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 1);
      dudz_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 1);
      dvdx_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 1);
      dvdy_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 1);
      dvdz_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 1);
      dwdx_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 1);
      dwdy_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 1);
      dwdz_.allocate(domaine_ijk_, Domaine_IJK::ELEM, 1);
      if (is_post_required("LAMBDA2"))
        {
          scalar_post_fields_.at("LAMBDA2").allocate(domaine_ijk_, Domaine_IJK::ELEM, 0, "LAMBDA2");
          champs_compris_.ajoute_champ(scalar_post_fields_.at("LAMBDA2"));
        }
      if (is_post_required("CRITERE_Q") || is_post_required("CURL"))
        {
          scalar_post_fields_.at("CRITERE_Q").allocate(domaine_ijk_, Domaine_IJK::ELEM, 0, "CRITERE_Q");
          champs_compris_.ajoute_champ(scalar_post_fields_.at("CRITERE_Q"));
          // Le rotationnel, aux elems aussi :
          allocate_cell_vector(vect_post_fields_.at("CURL"), domaine_ijk_, 0, "CURL");
          champs_compris_.ajoute_champ_vectoriel(vect_post_fields_.at("CURL"));
        }
    }

  if (interfaces_->nb_groups() > 1)
    {
      // On alloue un tableau assez grand pour contenir tous les groupes.
      if (interfaces_->nb_groups() > 3)
        {
          Cerr << "More than 3 groups are planned, but the allocated fields has only 3 components" << endl;
          Process::exit();
        }
      // TODO AYM: allocate fait dans IJK_Interfaces a l init
      // allocate_cell_vector(interfaces_.groups_indicatrice_n_ns(),splitting_, 1); // Besoin d'un ghost pour le calcul du grad
      // allocate_cell_vector(interfaces_.groups_indicatrice_n_ft(),splitting_ft_, 1); // peut-etre qu'un ghost=0 suffit et fonctionne pour le redistribute.
    }

  // Pour verification des stats :
  if (check_stats_)
    {
      // Le gradient de pression aux faces : (il est deja alloue par defaut)
      //    allocate_velocity(gradP_, splitting_, 1);
      // Le gradient de vitesse aux elems : (il est stocke, seulement si besoin, dans l'objet statistiques_FT_)
      //    allocate_cell_vector(dUd_, splitting_, 1);
      //    allocate_cell_vector(dVd_, splitting_, 1);
      //    allocate_cell_vector(dWd_, splitting_, 1);
      // Et leurs solutions analytiques sans ghost :
      allocate_velocity(ana_gradP_, domaine_ijk_, 1);
      allocate_cell_vector(ana_dUd_, domaine_ijk_, 0);
      allocate_cell_vector(ana_dVd_, domaine_ijk_, 0);
      allocate_cell_vector(ana_dWd_, domaine_ijk_, 0);
      // Pour les deriv secondes :
      allocate_cell_vector(ana_grad2Pi_, domaine_ijk_, 0);
      allocate_cell_vector(ana_grad2Pc_, domaine_ijk_, 0);
      // And for velocities components :
      allocate_cell_vector(ana_grad2Ui_, domaine_ijk_, 0);
      allocate_cell_vector(ana_grad2Uc_, domaine_ijk_, 0);
      allocate_cell_vector(ana_grad2Vi_, domaine_ijk_, 0);
      allocate_cell_vector(ana_grad2Vc_, domaine_ijk_, 0);
      allocate_cell_vector(ana_grad2Wi_, domaine_ijk_, 0);
      allocate_cell_vector(ana_grad2Wc_, domaine_ijk_, 0);
    }

  if (liste_post_instantanes_.contient_("CELL_FORCE_PH")||liste_post_instantanes_.contient_("TOUS"))
    allocate_cell_vector(cell_source_spectrale_, domaine_ijk_, 0);
  if (liste_post_instantanes_.contient_("CELL_GRAD_P"))
    allocate_cell_vector(cell_grad_p_, domaine_ijk_, 0);
  if (liste_post_instantanes_.contient_("CELL_SOURCE_QDM_INTERF")||liste_post_instantanes_.contient_("TOUS"))
    allocate_cell_vector(cell_source_interface_,domaine_ijk_, 0);
  if (liste_post_instantanes_.contient_("CELL_SHIELD_REPULSION")||liste_post_instantanes_.contient_("TOUS"))
    allocate_cell_vector(cell_repulsion_interface_,domaine_ijk_, 0);
  if (liste_post_instantanes_.contient_("CELL_RHO_SOURCE_QDM_INTERF")||liste_post_instantanes_.contient_("TOUS"))
    {
      allocate_cell_vector(cell_bk_tsi_ns_,domaine_ijk_, 1);
      allocate_cell_vector(cell_rho_Ssigma_,domaine_ijk_, 1);
    }
}

void Postprocessing_IJK::alloc_velocity_and_co()
{
  bool flag_variable_source = ref_ijk_ft_->eq_ns().get_flag_variable_source();
  // Le mot cle TOUS n'a pas encore ete compris comme tel.
  if ((liste_post_instantanes_.contient_("GRAD_INDICATRICE_FT")) || (liste_post_instantanes_.contient_("TOUS")))
    allocate_velocity(grad_I_ft_, domaine_ft_, 2);

  if (liste_post_instantanes_.contient_("D_VELOCITY_ANA"))
    allocate_velocity(d_velocity_ana_, domaine_ijk_, 1);
  if (liste_post_instantanes_.contient_("OP_CONV"))
    {
      allocate_velocity(op_conv_, domaine_ijk_, ref_ijk_ft_->eq_ns().d_velocity_[0].ghost()); // Il y a 1 ghost chez d_velocity_
      //                                          On veut qqch d'aligne pour copier les data() l'un dans l'autre
    }
  if (liste_post_instantanes_.contient_("CELL_OP_CONV"))
    {
      allocate_cell_vector(cell_op_conv_, domaine_ijk_, ref_ijk_ft_->eq_ns().d_velocity_[0].ghost()); // Il y a 1 ghost chez d_velocity_
      //                                          On veut qqch d'aligne pour copier les data() l'un dans l'autre
    }

  if (is_post_required("RHO_SOURCE_QDM_INTERF"))
    {
      IJK_Field_vector3_double& rho_Ssigma = vect_post_fields_.at("RHO_SOURCE_QDM_INTERF");
      allocate_velocity(rho_Ssigma, domaine_ijk_, 1, "RHO_SOURCE_QDM_INTERF");
      champs_compris_.ajoute_champ_vectoriel(rho_Ssigma);
      // TODO : probablement inutile maintenant:
      allocate_cell_vector(cell_rho_Ssigma_, domaine_ijk_, 0);
    }

  // Pour le calcul des statistiques diphasiques :
  // (si le t_debut_stat a ete initialise... Sinon, on ne va pas les calculer au cours de ce calcul)
  if (is_stats_plans_activated())
    {
      allocate_velocity(grad_I_ns_, domaine_ijk_, 1);
    }
  else if (flag_variable_source)
    allocate_velocity(grad_I_ns_, domaine_ijk_, 1);
}

void Postprocessing_IJK::improved_initial_pressure_guess(bool imp)
{
  if ((imp) || (liste_post_instantanes_.contient_("COORDS")))
    {
      Noms noms_coords; // on attend trois expressions
      noms_coords.dimensionner_force(3);
      noms_coords[0] = "X";
      noms_coords[1] = "Y";
      noms_coords[2] = "Z";
      allocate_velocity(coords_, domaine_ijk_, 1);
      for (int i = 0; i < 3; i++)
        {
          // Cette methode parcours ni(), nj() et nk() et donc pas les ghost...
          set_field_data(coords_[i], noms_coords[i]);
        }
    }
}

void Postprocessing_IJK::postraiter_thermals(bool stop)
{
  if (ref_ijk_ft_->has_thermals())
    thermals_->set_first_step_thermals_post(first_step_thermals_post_);

  // TODO review this:

//   if (ref_ijk_ft_->has_thermals())
//     if (stop || first_step_thermals_post_
//         || (nb_pas_dt_post_thermals_probes_ >= 0 && tstep_sauv % nb_pas_dt_post_thermals_probes_ == nb_pas_dt_post_thermals_probes_ - 1)
//         || (std::floor((current_time-timestep)/time_interval_post_thermals_probes_) < std::floor(current_time/time_interval_post_thermals_probes_)))
//       {
//         Cout << "tstep : " << tstep << finl;
//         thermals_->thermal_subresolution_outputs(nb_pas_dt_post_thermals_probes_);
//       }
}

bool Postprocessing_IJK::is_stats_bulles_activated() const
{
  return nb_pas_dt_post_stats_bulles_ > 0 || time_interval_post_stats_bulles_ >= 0.;
}

// Warning: this one looking at t_debut_statistiques_ too
bool Postprocessing_IJK::is_stats_plans_activated() const
{
  return (nb_pas_dt_post_stats_plans_ > 0 || time_interval_post_stats_plans_ >= 0.) && t_debut_statistiques_ >= 0.0;
}

bool Postprocessing_IJK::is_stats_cisaillement_activated() const
{
  return nb_pas_dt_post_stats_cisaillement_ > 0 || time_interval_post_stats_cisaillement_ >= 0.;
}

bool Postprocessing_IJK::is_stats_rmf_activated() const
{
  return nb_pas_dt_post_stats_rmf_ > 0 || time_interval_post_stats_rmf_ >= 0.;
}

void Postprocessing_IJK::postraiter_stats(bool stop)
{
  Schema_Temps_IJK_base& sch = ref_ijk_ft_->schema_temps_ijk();
  int tstep = sch.get_tstep(),
      tstep_init = sch.get_tstep_init();
  double current_time = sch.get_current_time(),
         timestep = sch.get_timestep();
  const Nom& nom_cas = nom_du_cas();

  const int tstep_sauv = tstep + tstep_init;

  // Helper function, indicating for a given nb_pas_dt_XXX if post for this XXX category should be triggered:
  auto should_post = [&](int nb_pas_dt, double time_interv)
  {
    if (stop
        || (nb_pas_dt >= 0 && tstep_sauv % nb_pas_dt == nb_pas_dt - 1)
        || (std::floor((current_time-timestep)/time_interv) < std::floor(current_time/time_interv)))
      return true;
    return false;
  };

  if (is_stats_bulles_activated() && should_post(nb_pas_dt_post_stats_bulles_, time_interval_post_stats_bulles_))
    ecrire_statistiques_bulles(0, nom_cas, current_time);
  if (is_stats_plans_activated() && should_post(nb_pas_dt_post_stats_plans_, time_interval_post_stats_plans_))
    posttraiter_statistiques_plans(current_time);
  if (is_stats_cisaillement_activated() && should_post(nb_pas_dt_post_stats_cisaillement_, time_interval_post_stats_cisaillement_))
    ecrire_statistiques_cisaillement(0, nom_cas, current_time);
  if (is_stats_rmf_activated() && should_post(nb_pas_dt_post_stats_rmf_, time_interval_post_stats_rmf_))
    ecrire_statistiques_rmf(0, nom_cas, current_time);
}

// NEW INTERPOLATING FUNCTION TO AVOID EULERIAN POINTS LYING IN THE BUBBLES OR ON THE INTERFACE
static void ijk_interpolate_implementation_bis(const IJK_Field_double& field, const DoubleTab& coordinates, ArrOfDouble& result, int skip_unknown_points, double value_for_bad_points,
                                               const IJK_Field_double& indic)
{
  const int ghost = field.ghost();
  const int ni = field.ni();
  const int nj = field.nj();
  const int nk = field.nk();
  const Domaine_IJK& geom = field.get_domaine();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  const Domaine_IJK::Localisation loc = field.get_localisation();
  double origin_x = geom.get_origin(DIRECTION_I) + ((loc == Domaine_IJK::FACES_J || loc == Domaine_IJK::FACES_K || loc == Domaine_IJK::ELEM) ? (dx * 0.5) : 0.);
  double origin_y = geom.get_origin(DIRECTION_J) + ((loc == Domaine_IJK::FACES_K || loc == Domaine_IJK::FACES_I || loc == Domaine_IJK::ELEM) ? (dy * 0.5) : 0.);
  double origin_z = geom.get_origin(DIRECTION_K) + ((loc == Domaine_IJK::FACES_I || loc == Domaine_IJK::FACES_J || loc == Domaine_IJK::ELEM) ? (dz * 0.5) : 0.);
  const int nb_coords = coordinates.dimension(0);
  const int gh = indic.ghost();
  result.resize_array(nb_coords);
  for (int idx = 0; idx < nb_coords; idx++)
    {
      const double x = coordinates(idx, 0);  // coordinate of the point where the interpolation is needed
      const double y = coordinates(idx, 1);
      const double z = coordinates(idx, 2);
      const double x2 = (x - origin_x) / dx;
      const double y2 = (y - origin_y) / dy;
      const double z2 = (z - origin_z) / dz;
      const int index_i = (int) (floor(x2)) - geom.get_offset_local(DIRECTION_I);  //index of the closest cell center is this defined only in the single processor?
      const int index_j = (int) (floor(y2)) - geom.get_offset_local(DIRECTION_J);
      const int index_k = (int) (floor(z2)) - geom.get_offset_local(DIRECTION_K);
      const int index_kp1 = index_k + 1;
      const double xfact = x2 - floor(x2);  // non-dimension distance inside the cell where the interpolation is requested
      const double yfact = y2 - floor(y2);
      const double zfact = z2 - floor(z2);

      // Compute the phase indicator function in all the eulerian points involved in the interpolation
      double c_1 = 0.;
      double c_2 = 0.;
      double c_3 = 0.;
      double c_4 = 0.;
      double c_5 = 0.;
      double c_6 = 0.;
      double c_7 = 0.;
      double c_8 = 0.;
      //The following loop is used to take into account the non-periodic condition along the z-axis
      // if the index of the point is outside the wall the cell is considered as vapor, therefore putting an invalid value
      if (!geom.get_periodic_flag(DIRECTION_K))
        // In this case the two domains (FT and NS) are equal in the z direction
        {
          const int kmin = geom.get_offset_local(DIRECTION_K);
          const int nktot = geom.get_nb_items_global(loc, DIRECTION_K);
          //Serial calculation
          if (nktot == nk)
            {
              if ((kmin == 0) || (kmin + nk == nktot)) //The conditions on walls may be imposed at the same time since there is no misunderstanding on the indices
                {
                  // The phase indicator function is correctly evaluated only inside the domain of the single processor
                  // otherwise its value remains equal to 0, leading to ad invalid value for the interpolation
                  if (index_k >= 0 && index_k < nk && index_kp1 >= 0 && index_kp1 < nk && index_i >= 0 && index_i < ni && index_j >= 0 && index_j < nj)
                    {
                      c_1 = indic(index_i, index_j, index_k);
                      if (index_i != ni - 1)
                        c_2 = indic(index_i + 1, index_j, index_k);
                      if (index_j != nj - 1)
                        c_3 = indic(index_i, index_j + 1, index_k);
                      if (index_i != ni - 1 && index_j != nj - 1)
                        c_4 = indic(index_i + 1, index_j + 1, index_k);

                      c_5 = indic(index_i, index_j, index_kp1);
                      if (index_i != ni - 1)
                        c_6 = indic(index_i + 1, index_j, index_kp1);
                      if (index_j != nj - 1)
                        c_7 = indic(index_i, index_j + 1, index_kp1);
                      if (index_i != ni - 1 && index_j != nj - 1)
                        c_8 = indic(index_i + 1, index_j + 1, index_kp1);
                    }
                }
            }
          //Parallel calculation
          //NB The following procedure does not fit for one processor
          //since the first condition assigns a value to the point outside the right wall which should be invalid
          //THE LOOP SHOULD BE IMPROVED
          //The two walls are addressed separately since they are elaborated by two different processors.
          else
            {
              if (kmin == 0) //on the left wall
                {
                  if (index_k >= 0 && index_kp1 >= 0 && index_i >= -gh && index_i < ni + gh && index_j >= -gh && index_j < nj + gh && index_kp1 < nk && index_k < nk)
                    {
                      c_1 = indic(index_i, index_j, index_k);
                      if (index_i != ni + gh - 1)
                        c_2 = indic(index_i + 1, index_j, index_k);
                      if (index_j != nj + gh - 1)
                        c_3 = indic(index_i, index_j + 1, index_k);
                      if (index_i != ni + gh - 1 && index_j != nj + gh - 1)
                        c_4 = indic(index_i + 1, index_j + 1, index_k);

                      c_5 = indic(index_i, index_j, index_kp1);
                      if (index_i != ni + gh - 1)
                        c_6 = indic(index_i + 1, index_j, index_kp1);
                      if (index_j != nj + gh - 1)
                        c_7 = indic(index_i, index_j + 1, index_kp1);
                      if (index_i != ni + gh - 1 && index_j != nj + gh - 1)
                        c_8 = indic(index_i + 1, index_j + 1, index_kp1);
                    }
                }
              else if (kmin + nk == nktot) // on the right wall
                {
                  if (index_k < nk && index_kp1 < nk && index_i >= -gh && index_i < ni + gh && index_j >= -gh && index_j < nj + gh && index_kp1 >= 0 && index_k >= 0)
                    {
                      //Process::Journal() << "Proc: " << Process::me() << " i,j,k: "
                      //                   << index_i << " " << index_j << " " << index_k << endl;
                      c_1 = indic(index_i, index_j, index_k);
                      if (index_i != ni + gh - 1)
                        c_2 = indic(index_i + 1, index_j, index_k);
                      if (index_j != nj + gh - 1)
                        c_3 = indic(index_i, index_j + 1, index_k);
                      if (index_i != ni + gh - 1 && index_j != nj + gh - 1)
                        c_4 = indic(index_i + 1, index_j + 1, index_k);

                      c_5 = indic(index_i, index_j, index_kp1);
                      if (index_i != ni + gh - 1)
                        c_6 = indic(index_i + 1, index_j, index_kp1);
                      if (index_j != nj + gh - 1)
                        c_7 = indic(index_i, index_j + 1, index_kp1);
                      if (index_i != ni + gh - 1 && index_j != nj + gh - 1)
                        c_8 = indic(index_i + 1, index_j + 1, index_kp1);
                    }
                }
              else  // points in processor not treating walls
                {
                  if (index_k < nk + gh && index_kp1 < nk + gh && index_i >= -gh && index_i < ni + gh && index_j >= -gh && index_j < nj + gh && index_kp1 >= -gh && index_k >= -gh)
                    {
                      c_1 = indic(index_i, index_j, index_k);
                      if (index_i < ni + gh - 1)
                        c_2 = indic(index_i + 1, index_j, index_k);
                      if (index_j < nj + gh - 1)
                        c_3 = indic(index_i, index_j + 1, index_k);
                      if (index_i < ni + gh - 1 && index_j < nj + gh - 1)
                        c_4 = indic(index_i + 1, index_j + 1, index_k);

                      c_5 = indic(index_i, index_j, index_kp1);
                      if (index_i < ni + gh - 1)
                        c_6 = indic(index_i + 1, index_j, index_kp1);
                      if (index_j < nj + gh - 1)
                        c_7 = indic(index_i, index_j + 1, index_kp1);
                      if (index_i < ni + gh - 1 && index_j < nj + gh - 1)
                        c_8 = indic(index_i + 1, index_j + 1, index_kp1);
                    }
                }
            }
        }

      else // ghost cells may be interrogated when the direction is periodic
        {
          if (index_k < nk + gh && index_kp1 < nk + gh && index_i >= -gh && index_i < ni + gh && index_j >= -gh && index_j < nj + gh && index_kp1 >= -gh && index_k >= -gh)
            {
              c_1 = indic(index_i, index_j, index_k);
              if (index_i < ni + gh - 1)
                c_2 = indic(index_i + 1, index_j, index_k);
              if (index_j < nj + gh - 1)
                c_3 = indic(index_i, index_j + 1, index_k);
              if (index_i < ni + gh - 1 && index_j < nj + gh - 1)
                c_4 = indic(index_i + 1, index_j + 1, index_k);

              c_5 = indic(index_i, index_j, index_kp1);
              if (index_i < ni + gh - 1)
                c_6 = indic(index_i + 1, index_j, index_kp1);
              if (index_j < nj + gh - 1)
                c_7 = indic(index_i, index_j + 1, index_kp1);
              if (index_i < ni + gh - 1 && index_j < nj + gh - 1)
                c_8 = indic(index_i + 1, index_j + 1, index_kp1);
            }
        }

      // Where at least of the points used in the interpolation is outside the liquid phase (at the interface or in the vapor phase)
      // the invalid value is assigned by the interpolation function
      bool ok_2 = (c_1 + c_2 + c_3 + c_4 + c_5 + c_6 + c_7 + c_8 >= 7.99);

// is point in the domain ? (ghost cells ok...)
      bool ok = (index_i >= -ghost && index_i < ni + ghost - 1) && (index_j >= -ghost && index_j < nj + ghost - 1) && (index_k >= -ghost && index_k < nk + ghost - 1);
      if (!ok or !ok_2)
        {
          if (skip_unknown_points)
            {
              result[idx] = value_for_bad_points;
              continue; // go to next point
            }
          else
            {
              // Error!
              Cerr << "Error in ijk_interpolate_implementation: request interpolation of point " << x << " " << y << " " << z << " which is outside of the domain on processor " << Process::me()
                   << finl;
              Process::exit();
            }
        }

      double r = (((1. - xfact) * field(index_i, index_j, index_k) + xfact * field(index_i + 1, index_j, index_k)) * (1. - yfact)
                  + ((1. - xfact) * field(index_i, index_j + 1, index_k) + xfact * field(index_i + 1, index_j + 1, index_k)) * (yfact)) * (1. - zfact)
                 + (((1. - xfact) * field(index_i, index_j, index_kp1) + xfact * field(index_i + 1, index_j, index_kp1)) * (1. - yfact)
                    + ((1. - xfact) * field(index_i, index_j + 1, index_kp1) + xfact * field(index_i + 1, index_j + 1, index_kp1)) * (yfact)) * (zfact);
      result[idx] = r;
    }
}
void ijk_interpolate_skip_unknown_points_bis(const IJK_Field_double& field, const DoubleTab& coordinates, ArrOfDouble& result, const double value_for_bad_points, const IJK_Field_double& indic)
{
  ijk_interpolate_implementation_bis(field, coordinates, result, 1 /* yes:skip unknown points */, value_for_bad_points, indic);
}

// Here begins the computation of the extended pressure fields
// ALL the calculations are performed on the extended domain (FT) and the final field is the redistributed on the physical one (Navier-Stokes)
void Postprocessing_IJK::compute_extended_pressures()
{
  if (!extended_pressure_computed_)
    return; // Leave the function if the extended fields are not necessary...

  Navier_Stokes_FTD_IJK& ns = ref_ijk_ft_->eq_ns();
  statistiques().begin_count(postraitement_counter_);
  // The following calculation is defined on the extended domain ft
  const Domaine_IJK& split_ft = domaine_ft_;
  const int ni = split_ft.get_nb_elem_local(DIRECTION_I);
  const int nj = split_ft.get_nb_elem_local(DIRECTION_J);
  const int nk = split_ft.get_nb_elem_local(DIRECTION_K);
  const double dx = split_ft.get_constant_delta(DIRECTION_I);
  const double dy = split_ft.get_constant_delta(DIRECTION_J);
  const double dz = split_ft.get_constant_delta(DIRECTION_K);
  int nbsom = 0;
  // ArrOfInt liste_composantes_connexes_dans_element;
  DoubleTab positions_liq(2 * nbsom, 3); // Table of coordinates where interpolation needs to be computed
  DoubleTab positions_vap(2 * nbsom, 3);
  IntTab crossed_cells(nbsom, 3); // Table to store i,j,k of cells crossed by the interface.

  //pressure field has to be extended from ns to ft
  ns.redistribute_to_splitting_ft_elem_.redistribute(pressure_, pressure_ft_);
  pressure_ft_.echange_espace_virtuel(pressure_ft_.ghost()); // 5 ghost cells needed to avoid invalid points in the vapor phase
  //dP_ft_.echange_espace_virtuel(gradP_ft_.ghost());

  //initialisation

  // The other points of the domain ( the one outside the crossed cells are still equal to original value of pressure)
  // This generates the discontinuities seen in Visit
  extended_pl_ft_ = pressure_ft_;
  extended_pv_ft_ = pressure_ft_;
  //for (int dir = 0; dir < 3; dir++)
  //  {
  // dP_ft_[dir].data() = 0.;
  //}

  int errcount_pext = 0;
  // i,j,k are the indices if the cells in the extended domain, for each processor
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          if ((interfaces_->In_ft()(i, j, k) > 1.e-6) && (1. - interfaces_->In_ft()(i, j, k) > 1.e-6))
            {
              Vecteur3 bary_facettes_dans_elem;
              Vecteur3 normale;
              double norm = 1.;
              double dist = 0.;
              for (int c = 0; c < 3; c++)
                {
                  normale[c] = interfaces_->get_norm_par_compo_itfc_in_cell_ft()[c](i, j, k);
                  bary_facettes_dans_elem[c] = interfaces_->get_bary_par_compo_itfc_in_cell_ft()[c](i, j, k);
                }
              norm = sqrt(normale[0] * normale[0] + normale[1] * normale[1] + normale[2] * normale[2]);
              //if (norm<0.95)
              //    Process::Journal() << "[WARNING-NORM] " << "Indicatrice[" << i <<"," << j << "," << k << "] = " << interfaces_.In_ft()(i,j,k)
              //                       << " norm= " << norm << finl;
              //  }
              if (norm < 1.e-8)
                {
                  // Process::Journal() << " nb_compo_traversantes " << nb_compo_traversantes << finl;
                  Process::Journal() << "Indicatrice[" << i << "," << j << "," << k << "] = " << interfaces_->In_ft()(i, j, k) << finl;
                  Process::Journal() << "[WARNING-Extended-pressure] on Proc. " << Process::me() << "Floating Point Exception is barely avoided (" << " normale " << normale[0] << " " << normale[1]
                                     << " " << normale[2] << " )" << finl;
                  Process::Journal() << " But we have no distance to extrapolate the pressure" << finl;
                  dist = 1.52 * sqrt(dx * dx + dy * dy + dz * dz) / 3.;
                  if (interfaces_->In_ft()(i, j, k) * (1 - interfaces_->In_ft()(i, j, k)) > 1.e-6)
                    {
                      Process::Journal() << "[WARNING-Extended-pressure] " << "Indicatrice[" << i << "," << j << "," << k << "] = " << interfaces_->In_ft()(i, j, k) << finl;
                      if (interfaces_->In_ft()(i, j, k) > 0.99)
                        {
                          Process::Journal() << "[WARNING-Extended-pressure] " << "Pressure_ft_ will be kept as an extension for p_liq pressure_[" << i << "," << j << "," << k << "] = "
                                             << pressure_ft_(i, j, k) << finl;
                          extended_pv_ft_(i, j, k) = 1.e20;
                        }
                      else
                        {
                          Process::Journal() << "[WARNING-Extended-pressure] " << "Pressure_ft_ will be kept as an extension for p_vap pressure_[" << i << "," << j << "," << k << "] = "
                                             << pressure_ft_(i, j, k) << finl;
                          extended_pl_ft_(i, j, k) = 1.e20;
                        }
                      continue; // This cell is not added to crossed cells.
                    }
                  else
                    {
                      // We still need a unit normal
                      for (int dir = 0; dir < 3; dir++)
                        if (normale[dir] != 0)
                          normale[dir] = (normale[dir] > 0.) ? 1.0 : -1.0;

                      norm = sqrt(normale[0] * normale[0] + normale[1] * normale[1] + normale[2] * normale[2]);
                      if (std::fabs(norm) < 1.e-10)
                        {
                          Process::Journal() << "[WARNING-Extended-pressure] ||normal|| < 1.e-10" << finl;
                        }
                    }
                }

              if (std::fabs(norm) < 1.e-10)
                {
                  Process::Journal() << "[WARNING-Extended-pressure] Even with precaution, the normal truely is zero in compute_extended_pressures()... " << finl;
                  Cerr << "We ignore the extrapolation and keep the local value... (hopefully rare enough)" << finl;
                  dist = 0.;
                  errcount_pext++;
                }
              else
                {
                  dist = 1.52 * (std::fabs(dx * normale[0]) + std::fabs(dy * normale[1]) + std::fabs(dz * normale[2])) / norm;
                  // 2020.04.15 : GB correction for non-isotropic meshes and closer extrapolation points:
                  // The previous version was looking very far away in the direction where the mesh is fine (dz in channel flows).
                  // Then, a lot of ghost would be required.
                  // This new version still goes to the same value of dist when nx=1 or when nz=1, but is sin between
                  // double eps = 1.e-20;
                  // dist = 1.52*std::min(min(dx/(std::fabs(normale[0])+eps),dy/(std::fabs(normale[1])+eps)),dz/(std::fabs(normale[2])+eps));
                }

              nbsom++;
              crossed_cells.resize(nbsom, 3, RESIZE_OPTIONS::COPY_INIT);
              positions_liq.resize(2 * nbsom, 3, RESIZE_OPTIONS::COPY_INIT);
              positions_vap.resize(2 * nbsom, 3, RESIZE_OPTIONS::COPY_INIT);

              crossed_cells(nbsom - 1, 0) = i;
              crossed_cells(nbsom - 1, 1) = j;
              crossed_cells(nbsom - 1, 2) = k;

              for (int dir = 0; dir < 3; dir++)
                {
                  // Four image points are calculated, two on each side of the interface
                  //liquid phase
                  positions_liq(2 * nbsom - 2, dir) = bary_facettes_dans_elem[dir] + dist * normale[dir]; // 1st point to be done...
                  positions_liq(2 * nbsom - 1, dir) = bary_facettes_dans_elem[dir] + 2 * dist * normale[dir]; // 2nd point to be done...
                  //vapor_phase
                  positions_vap(2 * nbsom - 2, dir) = bary_facettes_dans_elem[dir] - dist * normale[dir]; // 1st point to be done...
                  positions_vap(2 * nbsom - 1, dir) = bary_facettes_dans_elem[dir] - 2 * dist * normale[dir]; // 2nd point to be done...
                }
            }
        }



  errcount_pext = static_cast<int>(Process::mp_sum(errcount_pext));
  if (Process::je_suis_maitre() && errcount_pext)
    Cerr << "[WARNING-Extended-pressure] Error Count = " << errcount_pext << endl;

  // Interpolation on the image points
  // All the quantities are evaluated on the extended domain, both the pressure field, both the image points coordinates

  ArrOfDouble p_interp_liq(2 * nbsom);
  ArrOfDouble p_interp_vap(2 * nbsom);
  ijk_interpolate_skip_unknown_points(pressure_ft_, positions_vap, p_interp_vap, 1.e5 /*value for unknown points*/);
  ijk_interpolate_skip_unknown_points_bis(pressure_ft_, positions_liq, p_interp_liq, 1.e5 /* value for unknown points */, interfaces_->In_ft());

  // Extrapolation in the eulerian cells crossed by the interface
  int inval_pl_count = 0.;
  int inval_pv_count = 0.;
  for (int icell = 0; icell < nbsom; icell++)
    {
      const int i = crossed_cells(icell, 0);
      const int j = crossed_cells(icell, 1);
      const int k = crossed_cells(icell, 2);
      // const int nb_compo_traversantes = interfaces_.compute_list_compo_connex_in_element(mesh, elem, liste_composantes_connexes_dans_element);
      const int nb_compo_traversantes = interfaces_->nb_compo_traversantes(i, j, k);
      if (nb_compo_traversantes != 1)
        {
          extended_pv_ft_(i, j, k) = 1.e20;
          inval_pv_count++;
        }
      else
        {
          extended_pv_ft_(i, j, k) = 2 * p_interp_vap[2 * icell] - 1 * p_interp_vap[2 * icell + 1];
        }

      //for(int dir=0; dir<3; dir++)
      if (p_interp_liq[2 * icell + 1] == 1.e5)
        {
          if (p_interp_liq[2 * icell] == 1.e5)   //May changing the order affect the result??
            {
              extended_pl_ft_(i, j, k) = 1.e20;
              inval_pl_count++;
              //dP_ft_[dir](i,j,k) = 1.e20
            }
          else
            {
              extended_pl_ft_(i, j, k) = p_interp_liq[2 * icell];
              //dP_ft_[dir](i,j,k) = 1.e20
            }
        }
      else
        {
          if (p_interp_liq[2 * icell] != 1.e5)
            {
              extended_pl_ft_(i, j, k) = 2 * p_interp_liq[2 * icell] - 1 * p_interp_liq[2 * icell + 1];
              //dP_ft_[dir](i,j,k) = (p_interp_liq(2*icell) - 1*p_interp_liq(2*icell+1))*normale[dir]/dist;
            }
          else
            {
              extended_pl_ft_(i, j, k) = p_interp_liq[2 * icell + 1];
              //dP_ft_[dir](i,j,k) = 1.e20
            }
        }
    }

  inval_pl_count = static_cast<int>(Process::mp_sum(inval_pl_count));
  inval_pv_count = static_cast<int>(Process::mp_sum(inval_pv_count));
  if (Process::je_suis_maitre() && inval_pl_count)
    Cerr << "[WARNING-Extended-pressure] Invalid p_l cells Count = " << inval_pl_count << endl;
  if (Process::je_suis_maitre() && inval_pv_count)
    Cerr << "[WARNING-Extended-pressure] Invalid p_v cells Count = " << inval_pv_count << endl;

  // The previous evaluated extended pressure has to be recomputed on the real NS domain
  ns.redistribute_from_splitting_ft_elem_.redistribute(extended_pl_ft_, extended_pl_);
  ns.redistribute_from_splitting_ft_elem_.redistribute(extended_pv_ft_, extended_pv_);
  //ref_ijk_ft_.redistribute_from_splitting_ft_elem_.redistribute(dP_ft, dP_);

  extended_pl_.echange_espace_virtuel(extended_pl_.ghost());
  extended_pv_.echange_espace_virtuel(extended_pv_.ghost());
  statistiques().end_count(postraitement_counter_);
}

// Methode appelee lorsqu'on a mis "TOUS" dans la liste des champs a postraiter.
// Elle ajoute a la liste tous les noms de champs postraitables par IJK_Interfaces

void Postprocessing_IJK::posttraiter_tous_champs_thermique(Motcles& liste, const int idx) const
{
  liste.add("TEMPERATURE");
  liste.add("CP");
  liste.add("LAMBDA");
  liste.add("TEMPERATURE_ANA");
  liste.add("ECART_T_ANA");
  liste.add("DIV_LAMBDA_GRAD_T_VOLUME");
  liste.add("SOURCE_TEMPERATURE");
  liste.add("TEMPERATURE_ADIM_BULLES");
  liste.add("TEMPERATURE_PHYSIQUE_T");
  liste.add("TEMPERATURE_ADIMENSIONNELLE_THETA");
  liste.add("SOURCE_TEMPERATURE_ANA");
  liste.add("ECART_SOURCE_TEMPERATURE_ANA");
  liste.add("GRAD_T");
  liste.add("GRAD_T0");
  liste.add("GRAD_T1");
  liste.add("GRAD_T2");
  liste.add("T_RUST");
  liste.add("DIV_RHO_CP_T_V");

}

// Methode appelee lorsqu'on a mis "TOUS" dans la liste des champs a postraiter.
// Elle ajoute a la liste tous les noms de champs postraitables par IJK_Interfaces
void Postprocessing_IJK::posttraiter_tous_champs_energie(Motcles& liste, const int idx) const
{
  liste.add("TEMPERATURE");
  liste.add("CP");
  liste.add("LAMBDA");
  liste.add("TEMPERATURE_ANA");
  liste.add("ECART_T_ANA");
  liste.add("DIV_LAMBDA_GRAD_T_VOLUME");
  liste.add("SOURCE_TEMPERATURE");
  liste.add("TEMPERATURE_ADIM_BULLES");
  liste.add("TEMPERATURE_PHYSIQUE_T");
  liste.add("TEMPERATURE_ADIMENSIONNELLE_THETA");
  liste.add("SOURCE_TEMPERATURE_ANA");
  liste.add("ECART_SOURCE_TEMPERATURE_ANA");
  liste.add("GRAD_T");
  liste.add("GRAD_T0");
  liste.add("GRAD_T1");
  liste.add("GRAD_T2");
  liste.add("T_RUST");
  liste.add("DIV_RHO_CP_T_V");
}

// local utilitary functions for the method Postprocessing_IJK::get_max_timestep_for_post
namespace
{

// returns remainder of division between real numbers
double modulo(double dividend, double divisor)
{
  return ((std::floor(dividend/divisor) + 1)*divisor - dividend);
}

// Note : the (1+1e-12) safety factor ensures that the simulation reaches the target.
// Otherwise, the simulation time might fall just below the target due to numerical errors, not triggering the desired post.
// If you happen to miss an interval, you may adjust the factor to a higher value
void add_max_timestep(std::vector<double>& timesteps, double current_time, double interval)
{
  if (interval>0 && (modulo(current_time, interval) != 0))
    timesteps.push_back(modulo(current_time, interval)*(1+1e-12));

}

}

/// @brief Compute the max possible timestep to use during the next iteration
/// in order to not skip a time interval for postpro
/// @param current_time
/// @return max time step acceptable in order to land on a time for requirede post
double Postprocessing_IJK::get_max_timestep_for_post(double current_time) const
{

  std::vector<double> timesteps;
  timesteps.push_back(DMAXFLOAT);
  add_max_timestep(timesteps, current_time, time_interval_post_);
  add_max_timestep(timesteps, current_time, time_interval_post_thermals_probes_);
  add_max_timestep(timesteps, current_time, time_interval_post_stats_bulles_);
  add_max_timestep(timesteps, current_time, time_interval_post_stats_plans_);
  add_max_timestep(timesteps, current_time, time_interval_post_stats_cisaillement_);
  add_max_timestep(timesteps, current_time, time_interval_post_stats_rmf_);

  double min=*std::min_element(timesteps.begin(),timesteps.end());
  assert(min>0);
  return min;

}
