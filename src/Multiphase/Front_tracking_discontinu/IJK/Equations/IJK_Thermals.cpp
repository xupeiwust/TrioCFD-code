/****************************************************************************
* Copyright (c) 2023, CEA
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

#include <IJK_Thermals.h>
#include <Probleme_FTD_IJK.h>
#include <IJK_switch_FT.h>
#include <IJK_switch.h>
#include <Option_IJK.h>

Implemente_instanciable( IJK_Thermals, "IJK_Thermals", Equation_base) ;

Sortie& IJK_Thermals::printOn( Sortie& os ) const
{
  return os;
}

Entree& IJK_Thermals::readOn( Entree& is )
{
  if (liste_thermique_.size())
    liste_thermique_.vide();

  Nom accouverte = "{", accfermee = "}", virgule = ",";
  Motcle nom;
  is >> nom;

  if (nom != "{")
    Process::exit("Error while reading the IJK_Thermals list. One expected an opened bracket '{' to start !!! \n");

  while (1)
    {
      liste_thermique_.add(OWN_PTR(IJK_Thermal_base)());
      IJK_Thermal_base::typer_lire_thermal_equation(liste_thermique_.dernier(), is);
      is >> nom;
      if (nom == "}") return is;
      if (nom != ",")
        {
          Cerr << "Error while reading the IJK_Thermals list. One expected a comma ',' or a closing bracket '}' and not " << nom << finl;
          Process::exit("Check your data file ! \n");
        }
    }

  return is;
}

void IJK_Thermals::set_fichier_reprise(const char *lataname)
{
  if (!liste_thermique_.est_vide())
    for (auto& itr : liste_thermique_)
      itr->set_fichier_reprise(lataname);
}

const Nom& IJK_Thermals::get_fichier_reprise()
{
  assert(!liste_thermique_.est_vide());
  return liste_thermique_[0]->get_fichier_reprise();
}

const Milieu_base& IJK_Thermals::milieu() const
{
  if (le_fluide_.est_nul())
    {
      Cerr << "You forgot to associate a fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return le_fluide_.valeur();
}

Milieu_base& IJK_Thermals::milieu()
{
  if (le_fluide_.est_nul())
    {
      Cerr << "You forgot to associate a fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return le_fluide_.valeur();
}

void IJK_Thermals::verifie_milieu()
{
  const Fluide_Diphasique_IJK& mil = milieu_ijk();
  if (!mil.fluide_phase(0).has_capacite_calorifique() || !mil.fluide_phase(0).has_conductivite() ||
      !mil.fluide_phase(1).has_capacite_calorifique() || !mil.fluide_phase(1).has_conductivite())
    {
      Cerr << "Error in IJK_Thermals::verifie_milieu() !! " << finl;
      Cerr << "You forgot to define the cp and/or fields in your Fluide_Diphasique_IJK medium ! Update your data file." << finl;
      Process::exit();
    }

  const int nb_cp0 = mil.fluide_phase(0).capacite_calorifique().valeurs().size(),
            nb_lamb0 = mil.fluide_phase(0).conductivite().valeurs().size(),
            nb_cp1 = mil.fluide_phase(1).capacite_calorifique().valeurs().size(),
            nb_lamb1 = mil.fluide_phase(1).conductivite().valeurs().size();

  if (nb_cp0 != liste_thermique_.size() || nb_lamb0 != nb_cp0 || nb_cp1 != nb_cp0 || nb_lamb1 != nb_cp0)
    {
      Cerr << "Error in IJK_Thermals::verifie_milieu() !! " << finl;
      Cerr << "You define " << liste_thermique_.size() << " thermal equations in your list which is not coherent with the medium :"<< finl;
      Cerr << "       - Size cp vapour field = " << nb_cp0 << finl;
      Cerr << "       - Size cp liquid field = " << nb_cp1 << finl;
      Cerr << "       - Size lambda vapour field = " << nb_lamb0 << finl;
      Cerr << "       - Size lambda liquid field = " << nb_lamb1 << finl;
      Process::exit();
    }
}

void IJK_Thermals::associer_milieu_base(const Milieu_base& un_milieu)
{
  if (sub_type(Fluide_Diphasique_IJK, un_milieu))
    {
      const Milieu_base& un_fluide = ref_cast(Milieu_base, un_milieu);
      le_fluide_ = un_fluide;
    }
  else
    {
      Cerr << "Error of fluid type for the method IJK_Thermals::associer_milieu_base" << finl;
      Process::exit();
    }

  verifie_milieu();

  for (auto& itr : liste_thermique_)
    itr->associer_milieu_base(un_milieu);
}

void IJK_Thermals::associer_pb_base(const Probleme_base& pb)
{
  if (!sub_type(Probleme_FTD_IJK_base, pb))
    {
      Cerr << "Error for the method IJK_Interfaces::associer_pb_base\n";
      Cerr << " IJK_Thermals equation must be associated to\n";
      Cerr << " a Probleme_FTD_IJK_base problem type\n";
      Process::exit();
    }
  mon_probleme = pb;
  if (nom_ == "??")
    {
      nom_ = pb.le_nom();
      nom_ += que_suis_je();
    }

  ref_ijk_ft_ = ref_cast(Probleme_FTD_IJK_base, pb);
}

void IJK_Thermals::completer()
{
  associer_post(ref_ijk_ft_->get_post());
  associer_interface_intersections(ref_ijk_ft_->get_interface().get_intersection_ijk_cell(), ref_ijk_ft_->get_interface().get_intersection_ijk_face());
  for (auto& itr : liste_thermique_)
    {
      itr->associer(ref_ijk_ft_.valeur());
      itr->associer_ghost_fluid_fields(ghost_fluid_fields_);
    }
  ghost_fluid_fields_.associer(ref_ijk_ft_.valeur());
  if (!liste_thermique_.est_vide())
    retrieve_ghost_fluid_params();
}

void IJK_Thermals::retrieve_ghost_fluid_params()
{
  int compute_distance = 1;
  int compute_curvature = 1;
  int n_iter_distance = 3;
  int avoid_gfm_parallel_calls = 0;
  IJK_Field_local_double boundary_flux_kmin;
  IJK_Field_local_double boundary_flux_kmax;
  assert(!liste_thermique_.est_vide());
  liste_thermique_[0]->get_boundary_fluxes(boundary_flux_kmin, boundary_flux_kmax);
  for (auto& itr : liste_thermique_)
    itr->retrieve_ghost_fluid_params(compute_distance,
                                     compute_curvature,
                                     n_iter_distance,
                                     avoid_gfm_parallel_calls);
  ghost_fluid_fields_.retrieve_ghost_fluid_params(compute_distance,
                                                  compute_curvature,
                                                  n_iter_distance,
                                                  avoid_gfm_parallel_calls,
                                                  boundary_flux_kmin,
                                                  boundary_flux_kmax);
}

void IJK_Thermals::associer_post(const Postprocessing_IJK& ijk_ft_post)
{
  ref_ijk_ft_post_ = ijk_ft_post;
  for (auto& itr : liste_thermique_)
    itr->associer_post(ijk_ft_post);
}

void IJK_Thermals::associer_switch(const Switch_FT_double& ijk_ft_switch)
{
  ref_ijk_ft_switch_ = ijk_ft_switch;
  for (auto& itr : liste_thermique_)
    itr->associer_switch(ref_ijk_ft_switch_);
}

bool IJK_Thermals::has_IJK_field(const Nom& nom) const
{
  if (nom.contient("TEMPERATURE") || nom.contient("ECART_T") || nom== "INTERFACE_PHIN")
    return true;
  else
    return false;
}


const IJK_Field_double& IJK_Thermals::get_IJK_field(const Nom& nom) const
{
  for (int i = 0; i < (int) liste_thermique_.size(); i++)
    {
      if (nom== Nom("TEMPERATURE_")+Nom(i))
        return *((liste_thermique_.operator[](i))->get_temperature());

      if (nom== Nom("TEMPERATURE_ANA_")+Nom(i))
        return ((liste_thermique_.operator[](i))->get_temperature_ana());

      if (nom== Nom("ECART_T_ANA_")+Nom(i))
        return ((liste_thermique_.operator[](i))->get_ecart_t_ana());

      if (nom== Nom("INTERFACE_PHIN_")+Nom(i))
        return ((liste_thermique_.operator[](i))->get_grad_T_interface_ns());

      if (nom== Nom("TEMPERATURE_ADIMENSIONNELLE_THETA_")+Nom(i))
        return ((liste_thermique_.operator[](i))->get_temperature_adim_theta());
    }


  Cerr << "Erreur dans IJK_Thermals::get_IJK_field : " << finl;
  Cerr << "Le champ demande " << nom << " n'est pas connu par  IJK_Thermals::get_IJK_field." << finl;
  throw;
}

void IJK_Thermals::associer_interface_intersections(const Intersection_Interface_ijk_cell& intersection_ijk_cell,
                                                    const Intersection_Interface_ijk_face& intersection_ijk_face)
{
  ref_intersection_ijk_cell_ = intersection_ijk_cell;
  ref_intersection_ijk_face_ = intersection_ijk_face;
  for (auto& itr : liste_thermique_)
    itr->associer_interface_intersections(intersection_ijk_cell, intersection_ijk_face);
}

double IJK_Thermals::get_modified_time()
{
  double modified_time;
  if (liste_thermique_.est_vide())
    modified_time = ref_ijk_ft_->schema_temps_ijk().get_current_time();
  else
    modified_time = 0.;
  for (auto& itr : liste_thermique_)
    modified_time = std::max(modified_time, itr->get_modified_time());
  return modified_time;
}

void IJK_Thermals::get_rising_velocities_parameters(int& compute_rising_velocities,
                                                    int& fill_rising_velocities,
                                                    int& use_bubbles_velocities_from_interface,
                                                    int& use_bubbles_velocities_from_barycentres)
{
  if (!liste_thermique_.est_vide())
    for (auto& itr : liste_thermique_)
      itr->get_rising_velocities_parameters(compute_rising_velocities,
                                            fill_rising_velocities,
                                            use_bubbles_velocities_from_interface,
                                            use_bubbles_velocities_from_barycentres);
}

void IJK_Thermals::sauvegarder_temperature(Nom& lata_name,
                                           const int& stop)
{
  int idth = 0;
  for (auto& itr : liste_thermique_)
    {
      itr->sauvegarder_temperature(lata_name, idth, stop);
      idth++;
    }
}

void IJK_Thermals::sauvegarder_thermals(SFichier& fichier)
{
  int flag_list_not_empty_th = 0;
  if (liste_thermique_.size() > 0)
    {
      fichier << " thermals {\n" ;
      flag_list_not_empty_th = 1;
    }
  for(auto itr = liste_thermique_.begin(); itr != liste_thermique_.end(); )
    {
      fichier << *itr ;
      ++itr;
      if (itr != liste_thermique_.end())
        fichier << ", \n" ;
      else
        fichier << "\n" ;
    }
  if (flag_list_not_empty_th)
    fichier << " } \n" ;
}

void IJK_Thermals::compute_timestep(double& dt_thermals, const double dxmin)
{
  for (auto& itr : liste_thermique_)
    {
      const double dt_th = itr->compute_timestep(dt_thermals, dxmin);
      // We take the most restrictive of all thermal problems and use it for all:
      dt_thermals = std::min(dt_thermals, dt_th);
    }
}

void IJK_Thermals::initialize(const Domaine_IJK& splitting)
{
  if (!liste_thermique_.est_vide())
    {
      ghost_fluid_fields_.initialize(splitting);
      int idth =0;
      Nom thermal_outputs_rank_base = Nom("thermal_outputs_rank_");
      const int max_digit = 3;
      for (auto& itr : liste_thermique_)
        {
          itr->initialize(splitting, idth);
          if (!Option_IJK::DISABLE_DIPHASIQUE)
            itr->update_thermal_properties();
          const int max_rank_digit = idth < 1 ? 1 : (int) (log10(idth) + 1);
          thermal_rank_folder_.add(thermal_outputs_rank_base
                                   + Nom(std::string(max_digit - max_rank_digit, '0')) + Nom(idth));
          idth++;
        }
      if (!Option_IJK::DISABLE_DIPHASIQUE)
        {
          overall_bubbles_quantities_folder_ = Nom("overall_bubbles_quantities");
          interfacial_quantities_thermal_probes_folder_ = Nom("interfacial_quantities_thermal_probes");
          shell_quantities_thermal_probes_folder_ = Nom("shell_quantities_thermal_probes");
          local_quantities_thermal_probes_folder_ = Nom("local_quantities_thermal_probes");
          local_quantities_thermal_probes_time_index_folder_ = Nom("local_quantities_thermal_probes_time_index_");
          local_quantities_thermal_slices_folder_ = Nom("local_quantities_thermal_slices");
          local_quantities_thermal_slices_time_index_folder_ = Nom("local_quantities_thermal_slices_time_index_");
          local_quantities_thermal_lines_folder_ = Nom("local_quantities_thermal_lines");
          local_quantities_thermal_lines_time_index_folder_ = Nom("local_quantities_thermal_lines_time_index_");
        }
      for (auto& itr : liste_thermique_)
        {
          lata_step_reprise_.push_back(itr->get_latastep_reprise());
          lata_step_reprise_ini_.push_back(itr->get_latastep_reprise_ini());
        }
    }
}

void IJK_Thermals::recompute_temperature_init()
{
  for (auto& itr : liste_thermique_)
    itr->recompute_temperature_init();
}

int IJK_Thermals::size_thermal_problem(Nom thermal_problem)
{
  int size=0;
  for (auto& itr : liste_thermique_)
    {
      if (thermal_problem == itr->get_thermal_problem_type())
        size++;
    }
  return size;
}

void IJK_Thermals::update_thermal_properties()
{
  for (auto& itr : liste_thermique_)
    itr->update_thermal_properties();
}

void IJK_Thermals::euler_time_step(const double timestep)
{
  for (auto& itr : liste_thermique_)
    itr->euler_time_step(timestep);
  ghost_fluid_fields_.enforce_distance_curvature_values_for_post_processings();
}

void IJK_Thermals::euler_rustine_step(const double timestep)
{
  for (auto& itr : liste_thermique_)
    if (itr->get_thermal_problem_type() == Nom("onefluid"))
      {
        itr->update_thermal_properties();
        if (itr->get_conserv_energy_global())
          {
            const double dE = itr->get_E0() - itr->compute_global_energy();
            itr->euler_rustine_step(timestep, dE);
          }
      }
}

void IJK_Thermals::rk3_sub_step(const int rk_step, const double total_timestep, const double time)
{
  for (auto& itr : liste_thermique_)
    {
      int thermal_rank = itr->get_thermal_rank();
      switch (thermal_rank)
        {
        case 0:
          Cerr << "RK3 Time scheme is not implemented yet with" << itr->get_thermal_words()[thermal_rank] << finl;
          break;
        case 1:
          Cerr << "RK3 Time scheme is not implemented yet with" << itr->get_thermal_words()[thermal_rank] << finl;
          break;
        case 2:
          itr->rk3_sub_step(rk_step, total_timestep, time);
          Cerr << "RK3 Time scheme is implemented with" << itr->get_thermal_words()[thermal_rank] << finl;
          break;
        case 3:
          Cerr << "RK3 Time scheme is not implemented  with" << itr->get_thermal_words()[thermal_rank] << finl;
          break;
        case 4:
          itr->rk3_sub_step(rk_step, total_timestep, time);
          Cerr << "RK3 Time scheme is implemented with" << itr->get_thermal_words()[thermal_rank] << finl;
          break;
        default:
          Process::exit();
        }
    }
}

void IJK_Thermals::rk3_rustine_sub_step(const int rk_step, const double total_timestep,
                                        const double fractionnal_timestep, const double time)
{
  for (auto& itr : liste_thermique_)
    if (itr->get_thermal_problem_type() == Nom("onefluid") )
      {
        itr->update_thermal_properties();
        if (itr->get_conserv_energy_global())
          {
            const double dE = itr->get_E0() - itr->compute_global_energy();
            itr->rk3_rustine_sub_step(rk_step, total_timestep, fractionnal_timestep, time, dE);
          }
      }
}

void IJK_Thermals::ecrire_statistiques_bulles(int reset, const Nom& nom_cas, const double current_time, const ArrOfDouble& surface)
{
  int idx_th = 0;
  for (auto &itr : liste_thermique_)
    {
      itr->ecrire_statistiques_bulles(reset, nom_cas, current_time, surface, idx_th);
      ++idx_th;
    }
}

void IJK_Thermals::posttraiter_tous_champs_thermal(Motcles& liste_post_instantanes_)
{
  int idx_th = 0;
  for (auto &itr : liste_thermique_)
    {
      itr->posttraiter_tous_champs_thermal(liste_post_instantanes_, idx_th);
      ++idx_th;
    }
}

void IJK_Thermals::posttraiter_champs_instantanes_thermal(const Motcles& liste_post_instantanes,
                                                          const char *lata_name,
                                                          const int latastep,
                                                          const double current_time,
                                                          int& n)
{
  Cerr << "Post-process Eulerian fields related to the temperature resolution" << finl;
  int idx_th = 0;
  for (auto &itr : liste_thermique_)
    {
      int nb = itr->posttraiter_champs_instantanes_thermal(liste_post_instantanes,
                                                           lata_name,
                                                           latastep,
                                                           current_time,
                                                           idx_th);
      // Interfacial thermal fields :
      if (!Option_IJK::DISABLE_DIPHASIQUE)
        nb += itr->posttraiter_champs_instantanes_thermal_interface(liste_post_instantanes,
                                                                    lata_name,
                                                                    latastep,
                                                                    current_time,
                                                                    idx_th);
      if (idx_th == 0)
        n -= nb; // On compte comme "un" tous les CHAMPS_N (ou N est la longueur de la liste)
      ++idx_th;
    }
}

void IJK_Thermals::init_switch_thermals(const Domaine_IJK& splitting)
{
  int idx =0;
  for (auto& itr : liste_thermique_)
    {
      Cout << "Reading the old temperature field from " << Nom(itr->get_fichier_sauvegarde())
           << " to fill the liste_thermique_ field."<< finl;
      itr->initialize_switch(splitting, idx);
      idx++;
    }
}

void IJK_Thermals::prepare_thermals(const char *lata_name)
{
  for (auto& itr : liste_thermique_)
    itr->set_fichier_sauvegarde(lata_name);
}

void IJK_Thermals::ecrire_fichier_reprise(SFichier& fichier, const char *lata_name)
{
  Cerr << "  potentially saving temperature fields... " << finl;
  int flag_list_not_empty = 0;
  if ((int) liste_thermique_.size() > 0)
    {
      fichier << " thermals {\n" ;
      flag_list_not_empty = 1;
    }
  int idx =0;
  for (auto itr = liste_thermique_.begin(); itr != liste_thermique_.end(); )
    {
      (*itr)->set_fichier_sauvegarde(lata_name);
      fichier << *itr ;
      ++itr;
      if (itr != liste_thermique_.end() )
        fichier << ", \n" ;
      else
        fichier << "\n" ;
      Cerr << "  end of temperature field #" << idx << "... " << finl;
      ++idx;
    }
  if (flag_list_not_empty)
    fichier << " } \n" ;
}

int IJK_Thermals::ghost_fluid_flag()
{
  int ghost_fluid = 0;
  for (auto& itr : liste_thermique_)
    {
      ghost_fluid = itr->get_ghost_fluid_flag();
      if (ghost_fluid)
        return ghost_fluid;
    }
  return ghost_fluid;
}

void IJK_Thermals::compute_ghost_cell_numbers_for_subproblems(const Domaine_IJK& splitting, int ghost_init)
{
  for (auto& itr : liste_thermique_)
    itr->compute_ghost_cell_numbers_for_subproblems(splitting, ghost_init);
}

int IJK_Thermals::get_probes_ghost_cells(int ghost_init)
{
  int ghost_cells = ghost_init;
  for (auto& itr : liste_thermique_)
    {
      const int itr_ghost_cells = itr->get_ghost_cells();
      if (itr_ghost_cells > ghost_cells)
        ghost_cells = itr_ghost_cells;
    }
  return ghost_cells;
}

void IJK_Thermals::update_intersections()
{
  for (auto& itr : liste_thermique_)
    itr->update_intersections();
}

void IJK_Thermals::clean_ijk_intersections()
{
  for (auto& itr : liste_thermique_)
    itr->clean_ijk_intersections();
}

void IJK_Thermals::compute_eulerian_distance()
{
  assert(!liste_thermique_.est_vide());
  ghost_fluid_fields_.compute_eulerian_distance();
}

void IJK_Thermals::compute_eulerian_curvature()
{
  assert(!liste_thermique_.est_vide());
  ghost_fluid_fields_.compute_eulerian_curvature();
}

void IJK_Thermals::compute_eulerian_curvature_from_interface()
{
  assert(!liste_thermique_.est_vide());
  ghost_fluid_fields_.compute_eulerian_curvature_from_interface();
}

void IJK_Thermals::compute_eulerian_distance_curvature()
{
  if (!liste_thermique_.est_vide())
    if (ghost_fluid_flag())
      {
        compute_eulerian_distance();
        compute_eulerian_curvature_from_interface();
      }
}

int IJK_Thermals::get_disable_post_processing_probes_out_files() const
{
  int disable_post_processing_probes_out_files = 1;
  for (auto& itr : liste_thermique_)
    disable_post_processing_probes_out_files = (disable_post_processing_probes_out_files && itr->get_disable_post_processing_probes_out_files());
  return disable_post_processing_probes_out_files;
}

void IJK_Thermals::set_latastep_reprise(const bool stop)
{
  if (stop)
    for (auto& itr : liste_thermique_)
      itr->set_latastep_reprise(ref_ijk_ft_->schema_temps_ijk().get_tstep() + 1);
}

void IJK_Thermals::thermal_subresolution_outputs(const int& dt_post_thermals_probes)
{
  const int disable_post_processing_probes_out_files = get_disable_post_processing_probes_out_files();
  if (!disable_post_processing_probes_out_files && post_pro_first_call_)
    {
      if(!ini_folder_out_files_)
        {
          create_folders_for_probes();
          ini_folder_out_files_ = 1;
        }
      int rank = 0;
      for (auto& itr : liste_thermique_)
        {
          const int last_time = ref_ijk_ft_->schema_temps_ijk().get_tstep() + lata_step_reprise_ini_[rank];
          const int max_digit_time = 8;
          const int nb_digit_tstep = last_time < 1 ? 1 : (int) (log10(last_time) + 1);
          Nom prefix_local_quantities = thermal_rank_folder_[rank] + "/";
          Nom suffix_local_quantities = Nom(std::string(max_digit_time - nb_digit_tstep, '0')) + Nom(last_time);
          Nom local_quantities_thermal_probes_time_index_folder = prefix_local_quantities
                                                                  + local_quantities_thermal_probes_folder_ + "/"
                                                                  + local_quantities_thermal_probes_time_index_folder_
                                                                  + suffix_local_quantities;
          Nom overall_bubbles_quantities = thermal_rank_folder_[rank] + "/" + overall_bubbles_quantities_folder_;
          Nom interfacial_quantities_thermal_probes = thermal_rank_folder_[rank] + "/" + interfacial_quantities_thermal_probes_folder_;
          Nom shell_quantities_thermal_probes = thermal_rank_folder_[rank] + "/" + shell_quantities_thermal_probes_folder_;
          Nom local_quantities_thermal_slices_time_index_folder = prefix_local_quantities
                                                                  + local_quantities_thermal_slices_folder_ + "/"
                                                                  + local_quantities_thermal_slices_time_index_folder_
                                                                  + suffix_local_quantities;
          Nom local_quantities_thermal_lines_time_index_folder = prefix_local_quantities
                                                                 + local_quantities_thermal_lines_folder_ + "/"
                                                                 + local_quantities_thermal_lines_time_index_folder_
                                                                 + suffix_local_quantities;

          create_folders(local_quantities_thermal_probes_time_index_folder);
          create_folders(local_quantities_thermal_slices_time_index_folder);
          create_folders(local_quantities_thermal_lines_time_index_folder);

          itr->thermal_subresolution_outputs(interfacial_quantities_thermal_probes,
                                             shell_quantities_thermal_probes,
                                             overall_bubbles_quantities,
                                             local_quantities_thermal_probes_time_index_folder,
                                             local_quantities_thermal_slices_time_index_folder,
                                             local_quantities_thermal_lines_time_index_folder);
          // .sauv written before the post-processing on probes
          int latastep_reprise = lata_step_reprise_ini_[rank] + ref_ijk_ft_->schema_temps_ijk().get_tstep() + 2;
          const int nb_dt_max = ref_ijk_ft_->schema_temps_ijk().get_nb_timesteps();
          if ((ref_ijk_ft_->schema_temps_ijk().get_tstep() + dt_post_thermals_probes) >= nb_dt_max)
            latastep_reprise = nb_dt_max + 1;
          itr->set_latastep_reprise(latastep_reprise);
          rank++;
        }
    }
  post_pro_first_call_++;
}

void IJK_Thermals::create_folders_for_probes()
{
  for (int idth = 0; idth < liste_thermique_.size(); idth++)
    {
      create_folders(thermal_rank_folder_[idth]);
      create_folders(thermal_rank_folder_[idth] + "/" + overall_bubbles_quantities_folder_);
      create_folders(thermal_rank_folder_[idth] + "/" + interfacial_quantities_thermal_probes_folder_);
      create_folders(thermal_rank_folder_[idth] + "/" + shell_quantities_thermal_probes_folder_);
      create_folders(thermal_rank_folder_[idth] + "/" + local_quantities_thermal_probes_folder_);
      create_folders(thermal_rank_folder_[idth] + "/" + local_quantities_thermal_slices_folder_);
      create_folders(thermal_rank_folder_[idth] + "/" + local_quantities_thermal_lines_folder_);
    }
}

void IJK_Thermals::create_folders(Nom folder_name_base)
{
  if (Process::je_suis_maitre())
    {
      Nom spacing = " ";
      Nom folder_name = "mkdir -p";
      folder_name = folder_name + spacing + folder_name_base.getString(); // donothing";
      Cerr << "Shell command executed: " << folder_name << finl;
      int error = system(folder_name);
      assert(!error);
      if (error)
        Process::exit();
    }
}

void IJK_Thermals::set_first_step_thermals_post(int& first_step_thermals_post)
{
  first_step_thermals_post = 0;
  for (int idth = 0; idth < liste_thermique_.size(); idth++)
    first_step_thermals_post = (first_step_thermals_post || liste_thermique_[idth]->get_first_step_thermals_post());
}

void IJK_Thermals::set_temperature_ini()
{
  for (auto& itr : liste_thermique_)
    itr->compute_temperature_init();
}

void IJK_Thermals::recompute_interface_smoothing()
{
  set_temperature_ini();
  set_post_pro_first_call();
}

void IJK_Thermals::compute_new_thermal_field(Switch_FT_double& switch_double_ft,
                                             const Domaine_IJK& new_mesh,
                                             const Nom& lata_name,
                                             DoubleTab& coeff_i,
                                             IntTab Indice_i,
                                             DoubleTab& coeff_j,
                                             IntTab Indice_j,
                                             DoubleTab& coeff_k,
                                             IntTab Indice_k)
{
  IJK_Field_double new_thermal_field;
  if (liste_thermique_.size() > 0)
    {
      switch_double_ft.calculer_coords_elem();
      switch_double_ft.calculer_coeff(coeff_i,Indice_i,coeff_j,Indice_j,coeff_k,Indice_k);
      new_thermal_field.allocate(new_mesh /* it is in fact a splitting */, Domaine_IJK::ELEM, 0);
    }
  int idth = 0;
  for (auto& itr : liste_thermique_)
    {
      switch_double_ft.switch_scalar_field(*itr->get_temperature(),
                                           new_thermal_field,
                                           coeff_i, Indice_i,
                                           coeff_j ,Indice_j,
                                           coeff_k ,Indice_k);

      Cout << "Writing " << Nom("TEMPERATURE_") + Nom(idth) << " into " << lata_name << finl;
      dumplata_scalar(lata_name, Nom("TEMPERATURE_") + Nom(idth), new_thermal_field, 0 /*we store a 0 */);
      ++idth;
    }
}

void IJK_Thermals::copy_previous_interface_state()
{
  for (auto& itr : liste_thermique_)
    itr->copy_previous_interface_state();
}

void IJK_Thermals::copie_pure_vers_diph_sans_interpolation()
{
  for (auto& itr : liste_thermique_)
    itr->copie_pure_vers_diph_sans_interpolation();
}

void IJK_Thermals::echange_pure_vers_diph_cellules_initialement_pures()
{
  for (auto& itr : liste_thermique_)
    itr->echange_pure_vers_diph_cellules_initialement_pures();
}

void IJK_Thermals::echange_diph_vers_pure_cellules_finalement_pures()
{
  for (auto& itr : liste_thermique_)
    itr->echange_diph_vers_pure_cellules_finalement_pures();
}

void IJK_Thermals::vide_phase_invalide_cellules_diphasiques()
{
  for (auto& itr : liste_thermique_)
    itr->vide_phase_invalide_cellules_diphasiques();
}

void IJK_Thermals::remplir_tableau_pure_cellules_diphasiques(bool next_time)
{
  for (auto& itr : liste_thermique_)
    itr->remplir_tableau_pure_cellules_diphasiques(next_time);
}
