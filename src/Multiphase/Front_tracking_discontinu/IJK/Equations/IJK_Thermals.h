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
/////////////////////////////////////////////////////////////////////////////
//
// File      : IJK_Thermals.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_Thermals_included
#define IJK_Thermals_included

#include <IJK_Thermal_base.h>
#include <TRUST_List.h>
#include <System.h>


class Probleme_FTD_IJK_base;
class Switch_FT_double;

class IJK_Thermals : public Equation_base
{

  Declare_instanciable( IJK_Thermals ) ;

public :
  /*
   * Surcharge de l'eq base
   */
  void set_param(Param& titi) override { /* Do nothing */ }
  void completer() override;
  void associer_pb_base(const Probleme_base&) override;
  void discretiser() override { /* Do nothing */ }
  int preparer_calcul() override { return 1; }
  const Milieu_base& milieu() const override;
  Milieu_base& milieu() override;
  void associer_milieu_base(const Milieu_base& ) override ;
  int nombre_d_operateurs() const override { return -123; }
  const Operateur& operateur(int) const override { throw; }
  Operateur& operateur(int) override { throw; }
  const Champ_Inc_base& inconnue() const override { throw; }
  Champ_Inc_base& inconnue() override { throw; }

  void verifie_milieu();

  inline Fluide_Diphasique_IJK& milieu_ijk() { return ref_cast(Fluide_Diphasique_IJK, milieu()); }
  inline const Fluide_Diphasique_IJK& milieu_ijk() const { return ref_cast(Fluide_Diphasique_IJK, milieu()); }

  void set_fichier_reprise(const char *lataname);
  const Nom& get_fichier_reprise();
  void associer_post(const Postprocessing_IJK& ijk_ft_post);
  void associer_switch(const Switch_FT_double& ijk_ft_switch);
  bool has_IJK_field(const Nom& nom) const;
  const IJK_Field_double& get_IJK_field(const Nom& nom) const;
  void associer_interface_intersections(const Intersection_Interface_ijk_cell& intersection_ijk_cell_,
                                        const Intersection_Interface_ijk_face& intersection_ijk_face_);
  void retrieve_ghost_fluid_params();
  void sauvegarder_temperature(Nom& lata_name, const int& stop);
  void sauvegarder_thermals(SFichier& fichier);
  void compute_timestep(double& dt_thermals, const double dxmin);
  void initialize(const Domaine_IJK& splitting, int& nalloc);
  void recompute_temperature_init();
  void copie_pure_vers_diph_sans_interpolation();
  void echange_pure_vers_diph_cellules_initialement_pures();
  void echange_diph_vers_pure_cellules_finalement_pures();
  void vide_phase_invalide_cellules_diphasiques();
  void remplir_tableau_pure_cellules_diphasiques(bool next_time);
  int size_thermal_problem(Nom thermal_problem);
  void update_thermal_properties();
  void euler_time_step(const double timestep);
  void euler_rustine_step(const double timestep);
  void rk3_sub_step(const int rk_step, const double total_timestep, const double time);
  void rk3_rustine_sub_step(const int rk_step, const double total_timestep,
                            const double fractionnal_timestep, const double time);
  void ecrire_statistiques_bulles(int reset, const Nom& nom_cas, const double current_time, const ArrOfDouble& surface);
  void posttraiter_tous_champs_thermal(Motcles& liste_post_instantanes_);
  void posttraiter_champs_instantanes_thermal(const Motcles& liste_post_instantanes,
                                              const char *lata_name,
                                              const int latastep,
                                              const double current_time,
                                              int& n);
  int init_switch_thermals(const Domaine_IJK& splitting);
  void prepare_thermals(const char *lataname);
  int ghost_fluid_flag();
  void ecrire_fichier_reprise(SFichier& fichier, const char *lata_name);
  void compute_ghost_cell_numbers_for_subproblems(const Domaine_IJK& splitting, int ghost_init);
  int get_probes_ghost_cells(int ghost_init);

  void update_intersections();
  void clean_ijk_intersections();

  void compute_eulerian_distance();
  void compute_eulerian_curvature();
  void compute_eulerian_curvature_from_interface();
  void compute_eulerian_distance_curvature();

  void set_latastep_reprise(const bool stop);
  void thermal_subresolution_outputs(const int& dt_post_thermals_probes=0);
  int get_disable_post_processing_probes_out_files() const;
  double get_modified_time();
  void get_rising_velocities_parameters(int& compute_rising_velocities,
                                        int& fill_rising_velocities,
                                        int& use_bubbles_velocities_from_interface,
                                        int& use_bubbles_velocities_from_barycentres);
  void create_folders_for_probes();
  void create_folders(Nom folder_name_base);
  void set_first_step_thermals_post(int& first_step_thermals_post);
  void set_post_pro_first_call() { post_pro_first_call_ = 1; }
  void set_temperature_ini();
  void recompute_interface_smoothing();
  void compute_new_thermal_field(Switch_FT_double& switch_double_ft,
                                 const Domaine_IJK& new_mesh,
                                 const Nom& lata_name,
                                 DoubleTab& coeff_i,
                                 IntTab Indice_i,
                                 DoubleTab& coeff_j,
                                 IntTab Indice_j,
                                 DoubleTab& coeff_k,
                                 IntTab Indice_k);
  void copy_previous_interface_state();

  /*
   * Pour simplifier la vie
   */
  const LIST(OWN_PTR(IJK_Thermal_base))& get_liste_eqs() const { return liste_thermique_; }
  LIST(OWN_PTR(IJK_Thermal_base))& get_liste_eqs() { return liste_thermique_; }
  int size() const { return (int)liste_thermique_.size(); }
  int est_vide() const { return liste_thermique_.est_vide(); }

protected :
  // All Done here !!
  LIST(OWN_PTR(IJK_Thermal_base)) liste_thermique_;
  OBS_PTR(Milieu_base) le_fluide_;
  OBS_PTR(Probleme_FTD_IJK_base) ref_ijk_ft_;
  OBS_PTR(Postprocessing_IJK) ref_ijk_ft_post_;
  OBS_PTR(Switch_FT_double) ref_ijk_ft_switch_;
  OBS_PTR(Intersection_Interface_ijk_cell) ref_intersection_ijk_cell_;
  OBS_PTR(Intersection_Interface_ijk_face) ref_intersection_ijk_face_;

  IJK_Ghost_Fluid_Fields ghost_fluid_fields_;

  int post_pro_first_call_ = 0;

  System make_dir_for_out_files_;
  LIST(Nom) thermal_rank_folder_;
  Nom overall_bubbles_quantities_folder_;
  Nom interfacial_quantities_thermal_probes_folder_;
  Nom shell_quantities_thermal_probes_folder_;
  Nom local_quantities_thermal_probes_folder_;
  Nom local_quantities_thermal_probes_time_index_folder_;
  Nom local_quantities_thermal_slices_folder_;
  Nom local_quantities_thermal_slices_time_index_folder_;
  Nom local_quantities_thermal_lines_folder_;
  Nom local_quantities_thermal_lines_time_index_folder_;
  int ini_folder_out_files_ = 0;

  std::vector<int> lata_step_reprise_ini_, lata_step_reprise_;
};

#endif /* IJK_Thermals_included */
