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

#include <Schema_Euler_explicite_IJK.h>
#include <IJK_Ghost_Fluid_tools.h>
#include <Corrige_flux_FT_base.h>
#include <IJK_Thermal_base.h>
#include <IJK_Field_vector.h>
#include <Probleme_FTD_IJK.h>
#include <IJK_Bubble_tools.h>
#include <Cut_cell_tools.h>
#include <Schema_RK3_IJK.h>
#include <IJK_switch_FT.h>
#include <Postprocessing_IJK.h>
#include <Option_IJK.h>
#include <DebogIJK.h>
#include <Param.h>

Implemente_base_sans_constructeur( IJK_Thermal_base, "IJK_Thermal_base", Objet_U ) ;

/********************************************
 * Methods inherited from Objet_U
 ********************************************/

IJK_Thermal_base::IJK_Thermal_base()
{
  temperature_ = std::make_shared<IJK_Field_double>();
  RK3_F_temperature_ = std::make_shared<IJK_Field_double>();
  div_coeff_grad_T_volume_ = std::make_shared<IJK_Field_double>();
  d_temperature_ = std::make_shared<IJK_Field_double>();

  thermal_words_ = Motcles(5);
  {
    thermal_words_[0] = "subresolution";
    thermal_words_[1] = "multiplesubresolutions";
    thermal_words_[2] = "onefluid";
    thermal_words_[3] = "onefluidenergy";
    thermal_words_[4] = "cut_cell";
  }
  lata_suffix_ = Motcles(5);
  {
    lata_suffix_[0] = "SUBRES_";
    lata_suffix_[1] = "MSUBRES_";
    lata_suffix_[2] = "OF_";
    lata_suffix_[3] = "OFE_";
    lata_suffix_[4] = "CUT_CELL_";
  }
}

void IJK_Thermal_base::typer_lire_thermal_equation(OWN_PTR(IJK_Thermal_base)& eq, Entree& is)
{
  Cerr << "Read and Cast OWN_PTR(IJK_Thermal_base) => type : ";
  Nom prefix = "IJK_Thermal_";

  Motcles thermal_words = Motcles(5);
  {
    thermal_words[0] = "subresolution";
    thermal_words[1] = "multiplesubresolutions";
    thermal_words[2] = "onefluid";
    thermal_words[3] = "onefluidenergy";
    thermal_words[4] = "cut_cell";
  }
  Motcles lata_suffix = Motcles(5);
  {
    lata_suffix[0] = "SUBRES_";
    lata_suffix[1] = "MSUBRES_";
    lata_suffix[2] = "OF_";
    lata_suffix[3] = "OFE_";
    lata_suffix[4] = "CUT_CELL_";
  }

  Motcle word;
  is >> word;
  if (word.debute_par(prefix))
    word = Motcle((word.getString()).substr(prefix.getString().length()));

  Nom type = "";
  const int thermal_rank = thermal_words.search(word);
  type += prefix;

  switch(thermal_rank)
    {
    case 0:
      {
        type += "Subresolution";
        break;
      }
    case 1:
      {
        type += "Multiple_Subresolutions";
        break;
      }
    case 2:
      {
        type += "Onefluid";
        break;
      }
    case 3:
      {
        type += "Onefluid_Energy";
        break;
      }
    case 4:
      {
        type += "Cut_cell";
        break;
      }
    default:
      {
        Cerr << "ERROR : Thermal problems that are already implemented are:" << finl;
        Cerr << thermal_words << finl;
        abort();
      }
    }
  eq.typer(type);
  Cerr << type << " ... " << finl;
  eq->get_thermal_rank() = thermal_rank;
  eq->get_thermal_problem_type() = thermal_words[thermal_rank];
  is >> eq.valeur(); // Call the readOn
}

// X_D thermique listobj thermique -1 thermique_bloc 1 to add energy equation resolution if needed
// X_D thermique_bloc interprete nul 1 not_set
Entree& IJK_Thermal_base::readOn( Entree& is )
{
  /*
   * Parse the datafile
   */
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades(is);
  Cout << "IJK_Thermal_base::readOn : Parameters summary. " << finl;
  printOn(Cout);
  return is;
}

void IJK_Thermal_base::set_fichier_reprise(const char *lataname)
{
  fichier_reprise_temperature_ = lataname;
}

void IJK_Thermal_base::set_param(Param& param)
{
  param.ajouter("fo", &fo_); // X_D_ADD_P floattant not_set
  param.ajouter("expression_T_init", &expression_T_init_); // X_D_ADD_P chaine Expression of initial temperature (parser of x,y,z)
  param.ajouter("boundary_conditions", &boundary_conditions_, Param::REQUIRED); // X_D_ADD_P bloc_lecture boundary conditions
  param.ajouter("type_T_source", &type_T_source_); // X_D_ADD_P chaine(into=["dabiri","patch_dabiri","unweighted_dabiri"]) source term
  param.ajouter("expression_source_temperature", &expression_source_temperature_); // X_D_ADD_P chaine source terms
  param.ajouter_flag("lambda_variable", &lambda_variable_);
  param.ajouter_flag("wall_flux", &wall_flux_); // X_D_ADD_P rien not_set
  param.ajouter("kl_source", &kl_);
  param.ajouter("kv_source", &kv_);
  param.ajouter("T0l_source", &T0l_);
  param.ajouter("T0v_source", &T0v_);
  param.ajouter("fichier_reprise_temperature", &fichier_reprise_temperature_);
  param.ajouter("timestep_reprise_temperature", &timestep_reprise_temperature_);
  param.ajouter("rank_reprise_temperature", &rank_reprise_temperature_);
  param.ajouter("latastep_reprise", &latastep_reprise_ini_);
  param.ajouter_flag("conv_temperature_negligible", &conv_temperature_negligible_); // X_D_ADD_P rien neglect temperature convection
  param.ajouter_flag("diff_temperature_negligible", &diff_temperature_negligible_); // X_D_ADD_P rien neglect temperature diffusion
  param.ajouter("temperature_diffusion_op", &temperature_diffusion_op_);
  param.ajouter("temperature_convection_op", &temperature_convection_op_);
  param.ajouter("expression_T_ana", &expression_T_ana_); // X_D_ADD_P chaine Analytical expression T=f(x,y,z,t) for post-processing only
  param.ajouter("calculate_local_energy", &calculate_local_energy_);
  param.ajouter("upstream_temperature", &upstream_temperature_);
  param.ajouter("nb_diam_upstream", &nb_diam_upstream_);
  param.ajouter("side_temperature", &side_temperature_);
  param.ajouter("stencil_side", &stencil_side_);
  param.ajouter("n_iter_distance", &n_iter_distance_);
  param.ajouter_flag("ghost_fluid", &ghost_fluid_);
  param.ajouter_flag("spherical_exact", &spherical_exact_);
  param.ajouter_flag("debug", &debug_);
  param.ajouter_flag("compute_eulerian_compo", &compute_eulerian_compo_);

  param.ajouter_flag("store_flux_operators_for_energy_balance", &store_flux_operators_for_energy_balance_);
  param.ajouter_flag("disable_relative_velocity_energy_balance", &disable_relative_velocity_energy_balance_);

  param.ajouter_flag("use_bubbles_velocities_from_interface", &use_bubbles_velocities_from_interface_);
  param.ajouter_flag("use_bubbles_velocities_from_barycentres", &use_bubbles_velocities_from_barycentres_);

  param.ajouter_flag("smooth_grad_T_elem", &smooth_grad_T_elem_);
  param.ajouter("smoothing_numbers", &smoothing_numbers_);
  param.ajouter_flag("smoothing_remove_normal_compo", &smoothing_remove_normal_compo_);
  param.ajouter_flag("smoothing_use_unique_phase", &smoothing_use_unique_phase_);

  param.ajouter_flag("gfm_vapour_liquid_vapour", &gfm_vapour_liquid_vapour_);
  param.ajouter("gfm_smooth_factor", &gfm_smooth_factor_);
  param.ajouter_flag("avoid_gfm_parallel_calls", &avoid_gfm_parallel_calls_);
  //  param.ajouter_flag("gfm_recompute_field_ini", &gfm_recompute_field_ini_);
  //  param.ajouter_flag("gfm_zero_neighbour_value_mean", &gfm_zero_neighbour_value_mean_);

}

/********************************************
 * Public methods
 ********************************************/

const Milieu_base& IJK_Thermal_base::milieu() const
{
  if (le_fluide_.est_nul())
    {
      Cerr << "You forgot to associate a fluid to the problem named " << ref_ijk_ft_->le_nom() << finl;
      Process::exit();
    }
  return le_fluide_.valeur();
}

Milieu_base& IJK_Thermal_base::milieu()
{
  if (le_fluide_.est_nul())
    {
      Cerr << "You forgot to associate a fluid to the equation named " << ref_ijk_ft_->le_nom() << finl;
      Process::exit();
    }
  return le_fluide_.valeur();
}

void IJK_Thermal_base::associer_milieu_base(const Milieu_base& un_milieu)
{
  if (sub_type(Fluide_Diphasique_IJK, un_milieu))
    {
      const Milieu_base& un_fluide = ref_cast(Milieu_base, un_milieu);
      le_fluide_ = un_fluide;
    }
  else
    {
      Cerr << "Error of fluid type for the method IJK_Thermal_base::associer_milieu_base" << finl;
      Process::exit();
    }

  const Fluide_Diphasique_IJK& mil = milieu_ijk() ;
  cp_liquid_ = mil.get_cp_liquid(rang_);
  cp_vapour_ = mil.get_cp_vapour(rang_);
  lambda_liquid_ =mil.get_lambda_liquid(rang_);;
  lambda_vapour_ = mil.get_lambda_vapour(rang_);;
}

void IJK_Thermal_base::post_process_std_thermal_field(const Motcles& liste_post_instantanes,
                                                      const char *lata_name,
                                                      const int latastep,
                                                      const double current_time,
                                                      const int idx,
                                                      const Motcles& tested_names,
                                                      const Nom& name_field,
                                                      const Motcle& lata_suffix,
                                                      const IJK_Field_double& field,
                                                      std::ostringstream& oss,
                                                      int& counter,
                                                      const int& first_thermal_rank)
{
  oss << name_field << "_" << lata_suffix << idx;
  Nom name_field_tmp(oss.str().c_str());
  const int nb_names = tested_names.size();
  bool check_names = false;
  for (int i=0; i<nb_names; i++)
    check_names = check_names || (liste_post_instantanes.contient_(tested_names[i]));
  check_names = check_names || liste_post_instantanes.contient_(name_field_tmp);
  if (check_names)
    {
      counter++, dumplata_scalar(lata_name, name_field_tmp, field, latastep);
    }
  oss.str("");
}

int IJK_Thermal_base::initialize_switch(const Domaine_IJK& splitting, const int idx)
{
  int nalloc = 0;
  temperature_->allocate(splitting, Domaine_IJK::ELEM, 1);
  nalloc += 1;
  if (fichier_reprise_temperature_ == "??") // si on ne fait pas une reprise on initialise V
    {
      Cerr << "Please provide initial conditions for temperature, either by an expression or a field for restart."
           << "You should consider using either fichier_reprise_temperature or expression_T_init keywords. " << finl;
      Process::exit();
    }
  else
    {
      lire_temperature(splitting);
    }
  return nalloc;
}

int IJK_Thermal_base::initialize(const Domaine_IJK& splitting, const int idx)
{
  //  Cout << que_suis_je() << "::initialize()" << finl;
  rang_ = idx;
  int nalloc = 0;

  latastep_reprise_ = latastep_reprise_ini_;

  /*
   * Diffusion operator:
   * If temperature_diffusion_op_ is not written in the .data
   * the operator is initialised with uniform_lambda & centre2
   * in Operateur_IJK_elem_diff
   */
  if (needs_op_unform_)
    temperature_diffusion_op_.typer_diffusion_op("uniform");
  else
    temperature_diffusion_op_.typer_diffusion_op("standard");

  /*
   * Convection operator
   * If temperature_convection_op_ is not written in the .data
   * the operator is initialised with quick
   * in Operateur_IJK_elem_conv.h
   */

  /*
   * Initialise the operators
   */
  if (!conv_temperature_negligible_)
    temperature_convection_op_.initialize(splitting);
  if (!diff_temperature_negligible_)
    temperature_diffusion_op_.initialize(splitting);

  /*
   * Corrige Flux FT
   */
  corrige_flux_.typer("Corrige_flux_FT_temperature_conv");

  /*
   * Fields
   */
  temperature_->allocate(splitting, Domaine_IJK::ELEM, ghost_cells_);
  temperature_for_ini_per_bubble_.allocate(splitting, Domaine_IJK::ELEM, 1);

  d_temperature_->allocate(splitting, Domaine_IJK::ELEM, 2);
  nalloc += 3;
  compute_cell_volume();
  compute_min_cell_delta();

  /*
   * GFM
   */
  gfm_vapour_mixed_only_ = !gfm_vapour_liquid_vapour_;

  // if (!diff_temperature_negligible_)
  {
    div_coeff_grad_T_volume_->allocate(splitting, Domaine_IJK::ELEM, 2);
    nalloc += 1;
    div_coeff_grad_T_volume_->data() = 0.;
  }
  if (liste_post_instantanes_.contient_("DIV_LAMBDA_GRAD_T"))
    {
      div_coeff_grad_T_.allocate(splitting, Domaine_IJK::ELEM, 0);
      nalloc += 1;
      div_coeff_grad_T_.data() = 0.;
    }
  if (store_flux_operators_for_energy_balance_)
    {
      temperature_hess_flux_op_centre_.initialize(splitting);
      allocate_cell_vector(div_coeff_grad_T_raw_, splitting, 1);
      nalloc += 3;
      for (int c=0; c<3; c++)
        div_coeff_grad_T_raw_[c].data() = 0.;
    }
  // if (!conv_temperature_negligible_)
  {
    u_T_convective_volume_.allocate(splitting, Domaine_IJK::ELEM, 0);
    nalloc += 1;
    u_T_convective_volume_.data() = 0.;
  }
  if (liste_post_instantanes_.contient_("U_T_CONVECTIVE"))
    {
      u_T_convective_.allocate(splitting, Domaine_IJK::ELEM, 0);
      nalloc += 1;
      u_T_convective_.data() = 0.;
    }
  if (store_flux_operators_for_energy_balance_)
    {
      temperature_grad_flux_op_quick_.initialize(splitting);
      allocate_cell_vector(rho_cp_u_T_convective_raw_, splitting, 1);
      nalloc += 3;
      for (int c=0; c<3; c++)
        rho_cp_u_T_convective_raw_[c].data() = 0.;
    }

  rho_cp_post_ = (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("RHO_CP"));
  if (rho_cp_post_)
    {
      rho_cp_.allocate(splitting, Domaine_IJK::ELEM, 2);
      nalloc += 1;
    }
  if (calculate_local_energy_)
    {
      rho_cp_T_.allocate(splitting, Domaine_IJK::ELEM, 2);
      nalloc += 1;
    }

  /*
   * Storage for temperature gradient post-processing or method
   */
  if ((ref_ijk_ft_.non_nul()) && (!Option_IJK::DISABLE_DIPHASIQUE))
    {
      Cout << "Allocating fields temperature_ft_ and storage" << finl;
      allocate_cell_vector(storage_, ref_ijk_ft_->get_domaine_ft(), 1);
      nalloc += 3;
    }

  /*
   * Dimensionless temperature field (thermostat)
   */
  if ((wall_flux_) || liste_post_instantanes_.contient_("SOURCE_TEMPERATURE")
      || liste_post_instantanes_.contient_("TEMPERATURE_PHYSIQUE_T")
      || liste_post_instantanes_.contient_("TEMPERATURE_ADIMENSIONNELLE_THETA")
      || (type_T_source_ != "??"))
    {
      Cout << "Allocating field for the thermal source term & co. " << finl;
      source_temperature_.allocate(splitting, Domaine_IJK::ELEM, 1);
      source_temperature_v_.allocate(splitting, Domaine_IJK::ELEM, 1);
      source_temperature_l_.allocate(splitting, Domaine_IJK::ELEM, 1);
      d_source_Tv_.allocate(splitting, Domaine_IJK::ELEM, 1);
      d_source_Tl_.allocate(splitting, Domaine_IJK::ELEM, 1);
      temperature_physique_T_.allocate(splitting, Domaine_IJK::ELEM, 2);
      temperature_adimensionnelle_theta_.allocate(splitting, Domaine_IJK::ELEM, 2);
      nalloc += 7;
      // par defaut s'il n'y a pas de source renseignee, on utilise la source de Dabiri/Kawamura
      // cela veut dire que dans le cas des SWARMS il faut imperativement renseigner le nom de
      // la source
      if (type_T_source_ == "??")
        {
          Cerr << "Attention on demande des post-traitement sans avoir renseigner type_T_source" << finl;
          throw "Erreur post et type_T_source";
        }
    }
  if (liste_post_instantanes_.contient_("TEMPERATURE_ADIM_BULLES"))
    {
      temperature_adim_bulles_.allocate(splitting, Domaine_IJK::ELEM, 2);
      nalloc += 1;
    }

  /*
   * RK3 sub-steps
   * Check that pointer is not null:
   */
  if (ref_ijk_ft_.non_nul() && sub_type(Schema_RK3_IJK, ref_ijk_ft_->schema_temps_ijk()))
    {
      RK3_F_temperature_->allocate(splitting, Domaine_IJK::ELEM, 0);
      nalloc +=1;
    }

  if (liste_post_instantanes_.size() && (liste_post_instantanes_.contient_("TEMPERATURE_ANA")
                                         || liste_post_instantanes_.contient_("ECART_T_ANA") || liste_post_instantanes_.contient_("ECART_T_ANA_REL")))
    {
      temperature_ana_.allocate(splitting, Domaine_IJK::ELEM, 1);
      ecart_t_ana_.allocate(splitting, Domaine_IJK::ELEM, 1);
      ecart_t_ana_rel_.allocate(splitting, Domaine_IJK::ELEM, 1);
      nalloc +=3;
      temperature_ana_.data() = 0.;
      ecart_t_ana_.data() = 0.;
      ecart_t_ana_rel_.data() = 0.;
      temperature_ana_.echange_espace_virtuel(temperature_ana_.ghost());
      ecart_t_ana_.echange_espace_virtuel(ecart_t_ana_.ghost());
      ecart_t_ana_rel_.echange_espace_virtuel(ecart_t_ana_rel_.ghost());
    }

  /*
   * Resume the calculation
   */
  if (fichier_reprise_temperature_ == "??")   // si on ne fait pas une reprise on initialise V
    {
      if (expression_T_init_ != "??")
        {
          IJK_Thermal_base::compute_temperature_init();
        }
      else
        {
          Cerr << "Please provide initial conditions for temperature, either by an expression or a field for restart. "
               << "You should consider using either fichier_reprise_temperature or expression_T_init keywords. " << finl;
          Process::exit();
        }
    }
  else
    {
      lire_temperature(splitting);
    }

  /*
   * List of post-processed data
   */
  Cerr << " Initializing thermal fields dependant on the post-pro list : " ;
  if (liste_post_instantanes_.size())
    Cerr << liste_post_instantanes_;
  else
    Cerr << "empty";
  Cerr << finl;

  if ((liste_post_instantanes_.contient_("SOURCE_TEMPERATURE_ANA")) || (liste_post_instantanes_.contient_("ECART_SOURCE_TEMPERATURE_ANA")) )
    {
      source_temperature_ana_.allocate(splitting, Domaine_IJK::ELEM, 1);
      ecart_source_t_ana_.allocate(splitting, Domaine_IJK::ELEM, 1);
      nalloc += 2;
    }

  // TODO: Check with Aymeric
  calulate_grad_T_ = (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("GRAD_T"))
                     || (ref_ijk_ft_.non_nul() && ref_ijk_ft_->t_debut_statistiques() <  1.e10 );
  if (calulate_grad_T_)
    {
      allocate_velocity(grad_T_, splitting, 1);
      nalloc += 3;
    }

  compute_grad_T_interface_ = ghost_fluid_ || compute_grad_T_interface_ || (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("GRAD_T_INTERFACE"));
  compute_curvature_ = compute_curvature_ || compute_grad_T_interface_ || (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("CURVATURE"));
  compute_distance_ = compute_distance_ || compute_curvature_ || (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("DISTANCE"));

  if (compute_distance_)
    {
      eulerian_distance_ft_ = &(ghost_fluid_fields_->get_eulerian_distance_ft());
      eulerian_distance_ns_ = &(ghost_fluid_fields_->get_eulerian_distance_ns());
      eulerian_normal_vectors_ft_ = &(ghost_fluid_fields_->get_eulerian_normal_vectors_ft());
      eulerian_facets_barycentre_ft_ = &(ghost_fluid_fields_->get_eulerian_facets_barycentre_ft());
      eulerian_normal_vectors_ns_ = &(ghost_fluid_fields_->get_eulerian_normal_vectors_ns());
      eulerian_normal_vectors_ns_normed_ = &(ghost_fluid_fields_->get_eulerian_normal_vectors_ns_normed());
      eulerian_facets_barycentre_ns_ = &(ghost_fluid_fields_->get_eulerian_facets_barycentre_ns());
    }
  if (compute_curvature_)
    {
      eulerian_curvature_ft_ = &(ghost_fluid_fields_->get_eulerian_curvature_ft());
      eulerian_curvature_ns_ = &(ghost_fluid_fields_->get_eulerian_curvature_ns());
      eulerian_interfacial_area_ft_ = &(ghost_fluid_fields_->get_eulerian_interfacial_area_ft());
      eulerian_interfacial_area_ns_ = &(ghost_fluid_fields_->get_eulerian_interfacial_area_ns());
    }
  if (compute_grad_T_interface_)
    {
      // 1 ghost cell for eulerian_grad_T_interface_ and temperature_ft_ to access its neighbour
      eulerian_grad_T_interface_ft_.allocate(ref_ijk_ft_->get_domaine_ft(), Domaine_IJK::ELEM, 1);
      nalloc += 1;
      eulerian_grad_T_interface_ft_.echange_espace_virtuel(eulerian_grad_T_interface_ft_.ghost());
      //
      eulerian_grad_T_interface_ns_.allocate(splitting, Domaine_IJK::ELEM, 1);
      nalloc += 1;
      eulerian_grad_T_interface_ns_.echange_espace_virtuel(eulerian_grad_T_interface_ns_.ghost());
      //
    }

  temperature_ft_.allocate(ref_ijk_ft_->get_domaine_ft(), Domaine_IJK::ELEM, ghost_cells_);
  nalloc += 1;
  temperature_ft_.echange_espace_virtuel(eulerian_grad_T_interface_ft_.ghost());

  compute_eulerian_compo_ = compute_eulerian_compo_ || (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("EULERIAN_COMPO"))
                            || (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("EULERIAN_COMPO_NS"));
  if (compute_eulerian_compo_)
    {
      eulerian_compo_connex_ft_ = &(ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_eulerian_compo_connex_ft());
      eulerian_compo_connex_ns_ = &(ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_eulerian_compo_connex());
      eulerian_compo_connex_ghost_ft_ = &(ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_eulerian_compo_connex_ghost_ft());
      eulerian_compo_connex_ghost_ns_ = &(ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_eulerian_compo_connex_ghost());
      eulerian_compo_connex_from_interface_ft_ = &(ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_eulerian_compo_connex_from_interface_ft());
      eulerian_compo_connex_from_interface_ns_ = &(ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_eulerian_compo_connex_from_interface_ns());
      eulerian_compo_connex_from_interface_ghost_ft_ = &(ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_eulerian_compo_connex_from_interface_ghost_ft());
      eulerian_compo_connex_from_interface_ghost_ns_ = &(ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_eulerian_compo_connex_from_interface_ghost_ns());
      eulerian_compo_connex_from_interface_int_ns_ = &(ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_eulerian_compo_connex_int_from_interface_ns());
      eulerian_compo_connex_from_interface_ghost_int_ns_ = &(ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_eulerian_compo_connex_int_from_interface_ghost_ns());
    }

  compute_rising_velocities_ = compute_rising_velocities_ ||
                               (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("RISING_VELOCITIES"));
  fill_rising_velocities_ = compute_rising_velocities_ && (fill_rising_velocities_ ||
                                                           (liste_post_instantanes_.size() && liste_post_instantanes_.contient_("RISING_VELOCITIES")));

  if (compute_rising_velocities_)
    {
      rising_velocities_ = &(ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_rising_velocities());
      rising_velocities_from_barycentres_ = &(ref_ijk_ft_->get_interface().get_bubbles_velocities_magnitude_from_barycentres());
      rising_vectors_ = &(ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_rising_vectors());
      rising_vectors_from_barycentres_ = &(ref_ijk_ft_->get_interface().get_bubble_rising_vectors_from_barycentres());
      liquid_velocity_ = &(ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_liquid_velocity());
      rising_velocity_overall_ = &(ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_rising_velocity_overall());
      bubbles_barycentre_ = &(ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_bubbles_barycentre());
      bubbles_barycentres_old_ = &(ref_ijk_ft_->get_interface().get_bubble_barycentres_old_new(0));
      bubbles_barycentres_new_ = &(ref_ijk_ft_->get_interface().get_bubble_barycentres_old_new(1));
    }
  if (fill_rising_velocities_)
    eulerian_rising_velocities_ = &(ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_eulerian_rising_velocities());

  compute_hess_T_elem_ = compute_hess_T_elem_ || liste_post_instantanes_.contient_("HESS_T_ELEM");
  compute_hess_diag_T_elem_ = compute_hess_T_elem_ || compute_hess_diag_T_elem_ || liste_post_instantanes_.contient_("HESS_DIAG_T_ELEM")
                              || liste_post_instantanes_.contient_("HESS_XX_T_ELEM") || liste_post_instantanes_.contient_("HESS_YY_T_ELEM")
                              || liste_post_instantanes_.contient_("HESS_ZZ_T_ELEM");
  compute_hess_cross_T_elem_ = compute_hess_T_elem_ || compute_hess_cross_T_elem_ || liste_post_instantanes_.contient_("HESS_CROSS_T_ELEM")
                               || liste_post_instantanes_.contient_("HESS_XY_T_ELEM") || liste_post_instantanes_.contient_("HESS_XZ_T_ELEM")
                               || liste_post_instantanes_.contient_("HESS_YX_T_ELEM") || liste_post_instantanes_.contient_("HESS_YZ_T_ELEM")
                               || liste_post_instantanes_.contient_("HESS_ZX_T_ELEM") || liste_post_instantanes_.contient_("HESS_ZY_T_ELEM")
                               || liste_post_instantanes_.contient_("HESS_XZ_ZX_T_ELEM") || liste_post_instantanes_.contient_("HESS_ZX_XZ_T_ELEM")
                               || liste_post_instantanes_.contient_("HESS_YZ_ZY_T_ELEM") || liste_post_instantanes_.contient_("HESS_ZY_YZ_T_ELEM")
                               || liste_post_instantanes_.contient_("HESS_XY_YX_T_ELEM") || 	liste_post_instantanes_.contient_("HESS_YX_XY_T_ELEM");

  compute_hess_T_elem_ = compute_hess_diag_T_elem_ && compute_hess_cross_T_elem_;

  compute_grad_T_elem_ = compute_hess_cross_T_elem_ || compute_grad_T_elem_ || liste_post_instantanes_.contient_("GRAD_T_ELEM")
                         || liste_post_instantanes_.contient_("GRAD_T_DIR_X_ELEM") || liste_post_instantanes_.contient_("GRAD_T_DIR_Y_ELEM")
                         || liste_post_instantanes_.contient_("GRAD_T_DIR_Z_ELEM");
  smooth_grad_T_elem_ = smooth_grad_T_elem_ && compute_grad_T_elem_;
  if (compute_grad_T_elem_)
    {
      allocate_cell_vector(grad_T_elem_, splitting, ghost_cells_); // 1 or 0 ?
      nalloc += 3;
      grad_T_elem_.echange_espace_virtuel();
      temperature_grad_op_centre_.initialize(splitting);
    }
  if (smooth_grad_T_elem_)
    {
      temperature_gaussian_filtered_.allocate(splitting, Domaine_IJK::ELEM, ghost_cells_);
      nalloc += 1;
      temperature_gaussian_filtered_.echange_espace_virtuel(temperature_gaussian_filtered_.ghost());
      tmp_smoothing_field_.allocate(splitting, Domaine_IJK::ELEM, ghost_cells_);
      nalloc += 1;
      tmp_smoothing_field_.echange_espace_virtuel(tmp_smoothing_field_.ghost());
      allocate_cell_vector(grad_T_elem_smooth_, splitting, ghost_cells_); // 1 or 0 ?
      nalloc += 3;
      grad_T_elem_smooth_.echange_espace_virtuel();
      allocate_cell_vector(grad_T_elem_tangential_, splitting, ghost_cells_); // 1 or 0 ?
      nalloc += 3;
      grad_T_elem_tangential_.echange_espace_virtuel();

    }

  if (compute_hess_diag_T_elem_)
    {
      allocate_cell_vector(hess_diag_T_elem_, splitting, ghost_cells_);  // 1 or 0 ?
      nalloc += 3;
      hess_diag_T_elem_.echange_espace_virtuel();
      temperature_hess_op_centre_.initialize(splitting);
    }

  if (compute_hess_cross_T_elem_)
    {
      allocate_cell_vector(hess_cross_T_elem_, splitting, ghost_cells_);  // 1 or 0 ?
      nalloc += 3;
      hess_cross_T_elem_.echange_espace_virtuel();
      /*
       * TODO: Cross derivatives (adapt the diffusion operator ?)
       * Pb with Finite Volume Op, can not really relate on the fluxes
       * to derive the cross derivatives...
       */
    }
  spherical_approx_ = !spherical_exact_;

  /*
   * FIXME: Temporary need to rewrite IJK_Thermal.cpp posttraitements
   */
  allocate_cell_vector(dummy_int_vect_, splitting, 0); // 1 or 0 ?
  allocate_cell_vector(dummy_double_vect_, splitting, 0); // 1 or 0 ?
  dummy_int_field_.allocate(splitting, Domaine_IJK::ELEM, 0);
  dummy_double_field_.allocate(splitting, Domaine_IJK::ELEM, 0);
  nalloc += 8;
  for (int c=0; c<3; c++)
    {
      dummy_int_vect_[c].data() = 0;
      dummy_int_vect_[c].data() = 0.;
    }
  dummy_int_field_.data() = 0;
  dummy_double_field_.data() = 0;

  // ref_ijk_ft_->redistrib_from_ft_elem().redistribute(eulerian_grad_T_interface_, eulerian_grad_T_interface_);
  // Cout << "End of " << que_suis_je() << "::initialize()" << finl;
  return nalloc;
}

void IJK_Thermal_base::compute_temperature_init()
{
  Cout << "Temperature initialization from expression \nTini = " << expression_T_init_ << finl;
  set_field_data(*temperature_, expression_T_init_, ref_ijk_ft_->get_interface().I(), 0.);
}

void IJK_Thermal_base::recompute_temperature_init()
{
  compute_temperature_init();
  //  Cout << "Temperature initialization from expression \nTini = " << expression_T_init_ << finl;
  //  set_field_data(*temperature_, expression_T_init_, ref_ijk_ft_->itfce().In(), 0.);
}

void IJK_Thermal_base::copie_pure_vers_diph_sans_interpolation()
{
  Cerr << "copie_pure_vers_diph_sans_interpolation est seulement possible dans le cas IJK_Thermal_cut_cell." << finl;
  Process::exit();
}

void IJK_Thermal_base::echange_pure_vers_diph_cellules_initialement_pures()
{
  Cerr << "echange_pure_vers_diph_cellules_initialement_pures est seulement possible dans le cas IJK_Thermal_cut_cell." << finl;
  Process::exit();
}

void IJK_Thermal_base::echange_diph_vers_pure_cellules_finalement_pures()
{
  Cerr << "echange_diph_vers_pure_cellules_finalement_pures est seulement possible dans le cas IJK_Thermal_cut_cell." << finl;
  Process::exit();
}

void IJK_Thermal_base::vide_phase_invalide_cellules_diphasiques()
{
  Cerr << "vide_phase_invalide_cellules_diphasiques est seulement possible dans le cas IJK_Thermal_cut_cell." << finl;
  Process::exit();
}

void IJK_Thermal_base::remplir_tableau_pure_cellules_diphasiques(bool next_time)
{
  Cerr << "Remplir_cellules_diphasiques est seulement possible dans le cas IJK_Thermal_cut_cell." << finl;
  Process::exit();
}
Sortie& IJK_Thermal_base::printOn( Sortie& os ) const
{
  Nom front_space = "    ";
  Nom end_space = " ";
  Nom escape = "\n";
  Objet_U::printOn( os );
  os << "  {" << escape;

  os << escape;
  os << front_space << "# BASE PARAMS #" << escape;
  os << escape;

  os << "    boundary_conditions {"  << escape;
  /*
   * Boundary conditions (Periodicity or wall)
   */
  Nom bctype_kmin, bctype_kmax, bckmin, bckmax, valeur_kmin, valeur_kmax;
  if( boundary_conditions_.get_bctype_k_max()==boundary_conditions_.Paroi_Temperature_imposee)
    {
      bctype_kmax="Paroi_Temperature_imposee";
      bckmax = "temperature_imposee_kmax";
      valeur_kmax = boundary_conditions_.get_temperature_kmax();
    }
  else if( boundary_conditions_.get_bctype_k_max()==boundary_conditions_.Paroi_Flux_impose)
    {
      bctype_kmax="Paroi_Flux_impose";
      bckmax = "flux_impose_kmax";
      valeur_kmax = boundary_conditions_.get_flux_kmax();
    }
  else if( boundary_conditions_.get_bctype_k_max()==boundary_conditions_.Perio)
    {
      bctype_kmax="Perio";
      bctype_kmin="Perio";
      bckmax = " ";
      bckmin = " ";
      valeur_kmax = " ";
    }
  if( boundary_conditions_.get_bctype_k_min()==boundary_conditions_.Paroi_Temperature_imposee)
    {
      bctype_kmin="Paroi_Temperature_imposee";
      bckmin = "temperature_imposee_kmin";
      valeur_kmin = boundary_conditions_.get_temperature_kmin();
    }
  else if( boundary_conditions_.get_bctype_k_min()==boundary_conditions_.Paroi_Flux_impose)
    {
      bctype_kmin="Paroi_Flux_impose";
      bckmin = "flux_impose_kmin";
      valeur_kmin = boundary_conditions_.get_flux_kmin();
    }
  else if( boundary_conditions_.get_bctype_k_min()==boundary_conditions_.Perio)
    {
      bctype_kmin="Perio";
      bctype_kmin="Perio";
      bckmin = " ";
      bckmin = " ";
      valeur_kmin = " ";
    }
  os<< "      bctype_kmin" << " " << bctype_kmin << " \n";
  os<< "      bctype_kmax" << " " << bctype_kmax << " \n";
  os<< "      " << bckmin << " " << valeur_kmin << " \n";
  os<< "      " << bckmax << " " << valeur_kmax << " \n";
  os<< "    } \n" ;

  /*
   * Source term
   */
  if (type_T_source_!="??")
    os<< "    type_T_source " << type_T_source_ << "\n";
  if (type_T_source_=="SWARM")
    {
      os<< "      kl_source " <<  kl_ << "\n";
      os<< "      kv_source " <<  kv_ << "\n";
      os<< "      T0l_source " <<  T0l_ << "\n";
      os<< "      T0v_source " <<  T0v_ << "\n";
    }

  if( wall_flux_)
    os << front_space << "wall_flux" << escape;

  /*
   * Resume calculation
   */
  os << front_space << "fichier_reprise_temperature" << end_space << basename(fichier_reprise_temperature_)  << escape;
  os << front_space << "timestep_reprise_temperature" << end_space << timestep_reprise_temperature_ << escape;
  os << front_space << "rank_reprise_temperature" << end_space << rank_reprise_temperature_ << escape;
  os << front_space << "latastep_reprise" << end_space << latastep_reprise_ << escape;

  /*
   * Analytical expression of temperature at t_initial
   */
  if ( expression_T_ana_!="??")
    os << front_space << "expression_T_ana" <<  end_space << expression_T_ana_ << escape;


  os << front_space << "upstream_temperature" << end_space << upstream_temperature_ << escape;
  os << front_space << "nb_diam_upstream" << end_space << nb_diam_upstream_ << escape;
  os << front_space << "side_temperature" << end_space << side_temperature_ << escape;
  os << front_space << "stencil_side" << end_space << stencil_side_ << escape;
  os << front_space << "n_iter_distance" << end_space << n_iter_distance_ << escape;


  os << front_space << "temperature_diffusion_op" << end_space << temperature_diffusion_op_ << escape;
  os << front_space << "temperature_convection_op" << end_space << temperature_convection_op_ << escape;

  os << front_space << "gfm_smooth_factor" << end_space << gfm_smooth_factor_ << escape;

  os << front_space << "smoothing_numbers" << end_space << smoothing_numbers_ << escape;

  /*
   * Neglect an operator
   */

  os << escape;
  os << front_space << "# BASE FLAGS #" << escape;
  os << escape;

  if ( conv_temperature_negligible_)
    os << front_space << "conv_temperature_negligible" << escape;
  if ( diff_temperature_negligible_)
    os << front_space << "diff_temperature_negligible" << escape;
  if (ghost_fluid_)
    os << front_space << "ghost_fluid" <<  escape;
  if (spherical_exact_)
    os << front_space << "spherical_exact" <<  escape;
  if (debug_)
    os << front_space << "debug" <<  escape;
  if (calculate_local_energy_)
    os << front_space << "calculate_local_energy" <<  escape;

  if (store_flux_operators_for_energy_balance_)
    os << front_space << "store_flux_operators_for_energy_balance" <<  escape;
  if (disable_relative_velocity_energy_balance_)
    os << front_space << "disable_relative_velocity_energy_balance" <<  escape;

  if (use_bubbles_velocities_from_interface_)
    os << front_space << "use_bubbles_velocities_from_eulerian" <<  escape;
  if (use_bubbles_velocities_from_barycentres_)
    os << front_space << "use_bubbles_velocities_from_barycentres" <<  escape;

  if (smooth_grad_T_elem_)
    os << front_space << "smooth_grad_T_elem" <<  escape;
  if (smoothing_remove_normal_compo_)
    os << front_space << "smoothing_remove_normal_compo" <<  escape;
  if (smoothing_use_unique_phase_)
    os << front_space << "smoothing_use_unique_phase" <<  escape;
  if (compute_eulerian_compo_)
    os << front_space << "compute_eulerian_compo" <<  escape;

  if (avoid_gfm_parallel_calls_)
    os << front_space << "avoid_gfm_parallel_calls" <<  escape;

  return os;
}

double IJK_Thermal_base::get_modified_time()
{
  return ref_ijk_ft_->schema_temps_ijk().get_current_time();
}

void IJK_Thermal_base::get_rising_velocities_parameters(int& compute_rising_velocities,
                                                        int& fill_rising_velocities,
                                                        int& use_bubbles_velocities_from_interface,
                                                        int& use_bubbles_velocities_from_barycentres)
{
  compute_rising_velocities = compute_rising_velocities_ || compute_rising_velocities;
  fill_rising_velocities = fill_rising_velocities_ || fill_rising_velocities;
  use_bubbles_velocities_from_interface = use_bubbles_velocities_from_interface_ || use_bubbles_velocities_from_interface;
  use_bubbles_velocities_from_barycentres = use_bubbles_velocities_from_barycentres_ || use_bubbles_velocities_from_barycentres;
}

void IJK_Thermal_base::compute_cell_volume()
{
  const Domaine_IJK& geom = d_temperature_->get_domaine();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  vol_ = dx*dy*dz;
}

void IJK_Thermal_base::compute_min_cell_delta()
{
  const Domaine_IJK& geom = d_temperature_->get_domaine();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  min_delta_xyz_ = std::min(std::min(dx,dy),dz);
}

void IJK_Thermal_base::compute_cell_diagonal(const Domaine_IJK& geom)
{
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  cell_diagonal_ = sqrt(dx*dx + dy*dy + dz*dz);
}

const IJK_Field_double& IJK_Thermal_base::get_eulerian_compo_connex_ft() const
{
  if (!ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_compute_from_bounding_box())
    return dummy_double_field_;
  return *eulerian_compo_connex_ft_;
}
const IJK_Field_double& IJK_Thermal_base::get_eulerian_compo_connex_ghost_ft() const
{
  if (!ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_compute_from_bounding_box())
    return dummy_double_field_;
  return *eulerian_compo_connex_ghost_ft_;
}
const IJK_Field_double& IJK_Thermal_base::get_eulerian_compo_connex_from_interface_ft() const
{
  if (!ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_compute_compo_fields())
    return dummy_double_field_;
  return *eulerian_compo_connex_from_interface_ft_;
}
const IJK_Field_double& IJK_Thermal_base::get_eulerian_compo_connex_from_interface_ghost_ft() const
{
  if (!ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_compute_compo_fields())
    return dummy_double_field_;
  return *eulerian_compo_connex_from_interface_ghost_ft_;
}
const IJK_Field_double& IJK_Thermal_base::get_eulerian_compo_connex_ns() const
{
  if (!ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_compute_from_bounding_box())
    return dummy_double_field_;
  return *eulerian_compo_connex_ns_;
}
const IJK_Field_double& IJK_Thermal_base::get_eulerian_compo_connex_ghost_ns() const
{
  if (!ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_compute_from_bounding_box())
    return dummy_double_field_;
  return *eulerian_compo_connex_ghost_ns_;
}
const IJK_Field_double& IJK_Thermal_base::get_eulerian_compo_connex_from_interface_ns() const
{
  if (!ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_compute_compo_fields())
    return dummy_double_field_;
  return *eulerian_compo_connex_from_interface_ns_;
}
const IJK_Field_double& IJK_Thermal_base::get_eulerian_compo_connex_from_interface_ghost_ns() const
{
  if (!ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_compute_compo_fields())
    return dummy_double_field_;
  return *eulerian_compo_connex_from_interface_ghost_ns_;
}
const IJK_Field_int& IJK_Thermal_base::get_eulerian_compo_connex_int_from_interface_ns() const
{
  if (!ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_compute_compo_fields())
    return dummy_int_field_;
  return *eulerian_compo_connex_from_interface_int_ns_;
}
const IJK_Field_int& IJK_Thermal_base::get_eulerian_compo_connex_int_from_interface_ghost_ns() const
{
  if (!ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_compute_compo_fields())
    return dummy_int_field_;
  return *eulerian_compo_connex_from_interface_ghost_int_ns_;
}

void IJK_Thermal_base::update_thermal_properties()
{
  const IJK_Field_double& temperature = *temperature_;
  // Nombre de mailles du domaine NS :
  const int nx = temperature.ni();
  const int ny = temperature.nj();
  const int nz = temperature.nk();
  const double rhocpl = get_rhocp_l();
  const double ene_ini = compute_global_energy();
  if (Option_IJK::DISABLE_DIPHASIQUE)
    {
      if (rho_cp_post_) rho_cp_.data()= rhocpl;
      if (calculate_local_energy_)
        for (int k=0; k < nz ; k++)
          for (int j=0; j< ny; j++)
            for (int i=0; i < nx; i++)
              rho_cp_T_(i,j,k) = rhocpl*temperature(i,j,k);
    }
  else
    {
      const double rhocpv = get_rhocp_v();
      const IJK_Field_double& indic = ref_ijk_ft_->get_interface().I();
      for (int k=0; k < nz ; k++)
        for (int j=0; j< ny; j++)
          for (int i=0; i < nx; i++)
            {
              double chi_l = indic(i,j,k);
              if (rho_cp_post_)
                rho_cp_(i,j,k) = rhocpl*chi_l + rhocpv*(1-chi_l);
              if (calculate_local_energy_)
                rho_cp_T_(i,j,k) = (rhocpl*chi_l + rhocpv*(1-chi_l))*temperature(i,j,k);
            }
    }
  if (rho_cp_post_)
    rho_cp_.echange_espace_virtuel(rho_cp_.ghost());
  if (calculate_local_energy_)
    rho_cp_T_.echange_espace_virtuel(rho_cp_T_.ghost());

  // Semble un endroit approprie pour calculer la variation d'energie due au transport de l'interface:
  const double ene_post = compute_global_energy();
  Cerr << "[Energy-Budget-T"<<rang_<<"-2-TransportIndic] time t=" << ref_ijk_ft_->schema_temps_ijk().get_current_time()
       << " " << ene_ini
       << " " << ene_post
       << " delta=" << ene_post-ene_ini << " [W.m-3]." << finl;
}

// Methode de calcul du pas de temps max base sur Fo pour l'equation de thermique.
// CFL value is not computed as it is the same as for the velocity equation.
// The calculation should be stable if Fo <= 1.0 (thanks to the 0.5 in the formula below).
double IJK_Thermal_base::compute_timestep(const double timestep,
                                          const double dxmin)
{
  double alpha_max;
  double rho_l = ref_ijk_ft_->milieu_ijk().get_rho_liquid();
  double rho_v= ref_ijk_ft_->milieu_ijk().get_rho_vapour();
  if (Option_IJK::DISABLE_DIPHASIQUE)
    alpha_max = lambda_liquid_ / (rho_l * cp_liquid_);
  else
    alpha_max = std::max(lambda_liquid_ / (rho_l * cp_liquid_), lambda_vapour_ / (rho_v * cp_vapour_));
  dt_fo_ = dxmin * dxmin / (alpha_max + 1.e-20) * fo_ * 0.125; // 1/6 ou 1/8 ?
  if (diff_temperature_negligible_) dt_fo_ = 1.e20;
  return dt_fo_;
}

void IJK_Thermal_base::associer(const Probleme_FTD_IJK_base& ijk_ft)
{
  ref_ijk_ft_ = ijk_ft;
  liste_post_instantanes_ = ijk_ft.get_post().get_liste_post_instantanes();
  thermal_local_subproblems_interfaces_fields_.associer(ijk_ft);
}

void IJK_Thermal_base::associer_post(const Postprocessing_IJK& ijk_ft_post)
{
  ref_ijk_ft_post_ = ijk_ft_post;
}

void IJK_Thermal_base::associer_switch(const Switch_FT_double& ijk_ft_switch)
{
  ref_ijk_ft_switch_ = ijk_ft_switch;
}

void IJK_Thermal_base::associer_interface_intersections(const Intersection_Interface_ijk_cell& intersection_ijk_cell,
                                                        const Intersection_Interface_ijk_face& intersection_ijk_face)
{
  ref_intersection_ijk_cell_ = intersection_ijk_cell;
  ref_intersection_ijk_face_ = intersection_ijk_face;
}

void IJK_Thermal_base::associer_ghost_fluid_fields(const IJK_Ghost_Fluid_Fields& ghost_fluid_fields)
{
  ghost_fluid_fields_ = &ghost_fluid_fields;
}

void IJK_Thermal_base::retrieve_ghost_fluid_params(int& compute_distance,
                                                   int& compute_curvature,
                                                   int& n_iter_distance,
                                                   int& avoid_gfm_parallel_calls)
{
  compute_distance = compute_distance || compute_distance_;
  compute_curvature = compute_curvature || compute_curvature_;
  n_iter_distance = std::max(n_iter_distance, n_iter_distance_);
  avoid_gfm_parallel_calls = avoid_gfm_parallel_calls || avoid_gfm_parallel_calls_;
}

void IJK_Thermal_base::get_boundary_fluxes(IJK_Field_local_double& boundary_flux_kmin,
                                           IJK_Field_local_double& boundary_flux_kmax)
{
  boundary_flux_kmin = boundary_flux_kmin_;
  boundary_flux_kmax = boundary_flux_kmax_;
}

void IJK_Thermal_base::euler_time_step(const double timestep)
{
  static Stat_Counter_Id cnt_euler_thermal = statistiques().new_counter(1, "Euler Time Step - Temperature");
  static Stat_Counter_Id cnt_euler_thermal_post = statistiques().new_counter(2, "Euler Time Step - Temperature post");
  statistiques().begin_count(cnt_euler_thermal);

  if (debug_)
    Cerr << "Thermal Euler time-step" << finl;
  calculer_dT(ref_ijk_ft_->eq_ns().get_velocity());
  // Update the temperature :
  const int kmax = temperature_->nk();
  const double ene_ini = compute_global_energy();
  if (debug_)
    Cerr << "Apply temperature increment d_temperature" << finl;
  for (int k = 0; k < kmax; k++)
    ref_ijk_ft_->eq_ns().euler_explicit_update(*d_temperature_, *temperature_, k);

  /*
   * Erase the temperature increment (second call)
   */
  statistiques().begin_count(cnt_euler_thermal_post);
  if (debug_)
    Cerr << "Temperature after increment" << finl;
  post_process_after_temperature_increment();
  statistiques().end_count(cnt_euler_thermal_post);

  temperature_->echange_espace_virtuel(temperature_->ghost());
  const double ene_post = compute_global_energy();
  Cerr << "[Energy-Budget-T" << rang_ << "]"
       << " time t=" << ref_ijk_ft_->schema_temps_ijk().get_current_time()
       << " " << ene_ini
       << " " << ene_post << " [W.m-3]." << finl;
  source_callback();

  statistiques().end_count(cnt_euler_thermal);
}

void IJK_Thermal_base::rk3_sub_step(const int rk_step, const double total_timestep,
                                    const double time)
{
  calculer_dT(ref_ijk_ft_->eq_ns().get_velocity());
  // Update the temperature :
  const int kmax = temperature_->nk();
  const double ene_ini = compute_global_energy();
  for (int k = 0; k < kmax; k++)
    {
      runge_kutta3_update(*d_temperature_, *RK3_F_temperature_, *temperature_, rk_step, k, total_timestep);
    }
  temperature_->echange_espace_virtuel(temperature_->ghost());
  const double ene_post = compute_global_energy();
  Cerr << "[Energy-Budget-T"<<rang_<<"] time t=" << ref_ijk_ft_->schema_temps_ijk().get_current_time()
       << " " << ene_ini
       << " " << ene_post << " [W.m-3]. [step"<< rk_step << "]" << finl;
  source_callback();
}

void IJK_Thermal_base::sauvegarder_temperature(Nom& lata_name, int idx, const int& stop)
{
  fichier_reprise_temperature_ = lata_name;
  timestep_reprise_temperature_ = 1;
  rank_reprise_temperature_ = 0;
  dumplata_scalar(lata_name, Nom("TEMPERATURE_") + Nom(idx) , *temperature_, rank_reprise_temperature_);
  if (stop)
    latastep_reprise_ = latastep_reprise_ini_ + ref_ijk_ft_->schema_temps_ijk().get_tstep() + 1;
}

/********************************************
 * Protected methods
 ********************************************/

void IJK_Thermal_base::lire_temperature(const Domaine_IJK& splitting)
{
  if (rank_reprise_temperature_ ==-1) rank_reprise_temperature_ = rang_;
  Cout << "Reading initial temperature field T" << rang_ << " from file " << fichier_reprise_temperature_ << " timestep= " << timestep_reprise_temperature_
       << " from rank " << rank_reprise_temperature_ << finl;
  const Nom& geom_name = splitting.le_nom();
  lire_dans_lata(fichier_reprise_temperature_, timestep_reprise_temperature_, geom_name, Nom("TEMPERATURE_") + Nom(rank_reprise_temperature_),
                 *temperature_); // fonction qui lit un champ a partir d'un lata .
  temperature_->echange_espace_virtuel(temperature_->ghost()); // It is essential to fill the EV because the first call to convection needs them.
}

// Mettre rk_step = -1 si schema temps different de rk3.
void IJK_Thermal_base::calculer_dT(const IJK_Field_vector3_double& velocity)
{
  static Stat_Counter_Id cnt_gfm_temperature = statistiques().new_counter(2, "GFM - Temperature");
  static Stat_Counter_Id cnt_temperature_grad_hess = statistiques().new_counter(2, "Gradient and Hessian Temperature calculation");
  static Stat_Counter_Id cnt_lrs_temperature = statistiques().new_counter(2, "Solve the LRS thermal problems");
  static Stat_Counter_Id cnt_lrs_compute_flux = statistiques().new_counter(2, "Compute flux correction from LRS");
  static Stat_Counter_Id cnt_lrs_prepare_flux = statistiques().new_counter(2, "Prepare flux correction from LRS");
  static Stat_Counter_Id cnt_cell_temperature_first = statistiques().new_counter(2, "Temperature Cell centre - First call");
  static Stat_Counter_Id cnt_flux_balance = statistiques().new_counter(2, "Thermal flux balance");
  static Stat_Counter_Id cnt_upstream_temperature = statistiques().new_counter(2, "Upstream temperature");
  // static Stat_Counter_Id cnt_cell_temperature_second = statistiques().new_counter(1, "Temperature Cell centre - Second call");

  const double current_time = ref_ijk_ft_->schema_temps_ijk().get_current_time();
  const double ene_ini = compute_global_energy(*d_temperature_);

  /*
   * Clean_subproblems !
   */
  if (debug_)
    Cerr << "Clean thermal subproblems" << finl;
  clean_thermal_subproblems();

  // Correct the vapour and mixed cells values
  if (debug_)
    Cerr << "Store temperature before extrapolation" << finl;
  store_temperature_before_extrapolation();
  correct_temperature_for_eulerian_fluxes();

  /*
   * Correct the temperature field using either the ghost-fluid
   * approach or the laminar sub-resolution approach (and zero values for debug)
   */
  statistiques().begin_count(cnt_gfm_temperature);
  if (debug_)
    Cerr << "Start the Ghost-fluid (GFM) approach" << finl;
  if (debug_)
    Cerr << "Br0 (GFM approach)" << finl;
  compute_eulerian_grad_T_interface();
  if (debug_)
    Cerr << "Br1 (GFM approach)" << finl;
  propagate_eulerian_grad_T_interface();
  if (debug_)
    Cerr << "Br2 (GFM approach)" << finl;
  compute_eulerian_temperature_ghost();
  if (debug_)
    Cerr << "Br3 (GFM approach)" << finl;
  compute_eulerian_bounding_box_fill_compo();
  if (debug_)
    Cerr << "End the Ghost-fluid (GFM) approach" << finl;
  statistiques().end_count(cnt_gfm_temperature);

  /*
   * Compute gradients and hessian of the temperature after the ghost fluid extension
   */
  statistiques().begin_count(cnt_temperature_grad_hess);
  if (debug_)
    Cerr << "Compute temperature derivatives" << finl;
  compute_temperature_gradient_elem();
  compute_temperature_hessian_diag_elem();
  compute_temperature_hessian_cross_elem();
  statistiques().end_count(cnt_temperature_grad_hess);

  /*
   * Compute sub-problems (For Subresolution Child classes only !)
   */
  statistiques().begin_count(cnt_lrs_temperature);
  if (debug_)
    Cerr << "Compute thermal subproblems" << finl;
  compute_thermal_subproblems();
  statistiques().end_count(cnt_lrs_temperature);

  /*
   * Interpolate a value for the QUICK SCHEME (first call)
   */

  statistiques().begin_count(cnt_cell_temperature_first);
  if (debug_)
    Cerr << "Compute temperature mixed cell centres" << finl;
  compute_temperature_cell_centres(0);
  statistiques().end_count(cnt_cell_temperature_first);

  /*
   * Convective and Diffusive fluxes
   */
  statistiques().begin_count(cnt_lrs_compute_flux);
  if (debug_)
    Cerr << "Compute thermal convective and diffusive fluxes from subproblems" << finl;
  compute_convective_diffusive_fluxes_face_centre();
  statistiques().end_count(cnt_lrs_compute_flux);

  statistiques().begin_count(cnt_lrs_prepare_flux);
  if (debug_)
    Cerr << "Prepare ij fluxes" << finl;
  if (!conv_temperature_negligible_ || !diff_temperature_negligible_)
    prepare_ij_fluxes_k_layers();
  statistiques().end_count(cnt_lrs_prepare_flux);


  statistiques().begin_count(cnt_flux_balance);
  if (debug_)
    Cerr << "Compute fluxes balance" << finl;
  compute_temperature_diffusive_fluxes();
  compute_temperature_convective_fluxes(velocity);
  compare_fluxes_thermal_subproblems();
  statistiques().end_count(cnt_flux_balance);


  statistiques().begin_count(cnt_upstream_temperature);

  Navier_Stokes_FTD_IJK& ns = ref_ijk_ft_->eq_ns();
  if (debug_)
    Cerr << "Apply Upstream temperature" << finl;
  double nb_diam_upstream_velocity = ns.get_nb_diam_upstream();
  if (nb_diam_upstream_ == 0.)
    nb_diam_upstream_ = nb_diam_upstream_velocity;
  if (upstream_temperature_ > -1e20 && ns.get_vitesse_upstream() > -1e20)
    force_upstream_temperature(*temperature_, upstream_temperature_,
                               ref_ijk_ft_->get_interface(), nb_diam_upstream_,
                               ns.get_upstream_dir(), ref_ijk_ft_->milieu_ijk().get_direction_gravite(),
                               ns.get_upstream_stencil());
  statistiques().end_count(cnt_upstream_temperature);

  if (debug_)
    Cerr << "Convection of temperature" << finl;
  compute_temperature_convection(velocity);
  const double ene_postConv = compute_global_energy(*d_temperature_);

  if (debug_)
    Cerr << "Diffusion of temperature" << finl;
  add_temperature_diffusion();
  const double ene_postDiffu = compute_global_energy(*d_temperature_);

  if (debug_)
    Cerr << "Temperature source term" << finl;
  add_temperature_source();
  const double ene_postSource = compute_global_energy(*d_temperature_);

  /*
   * In case of the subresolution or not
   */
  set_zero_temperature_increment();
  // calculer_gradient_temperature(*temperature_, grad_T_); Routine Aymeric gradient sur faces

  Cerr << "[Energy-Budget-T"<<rang_<<"-1-TimeResolution] time t=" << current_time
       << " " << ene_ini
       << " " << ene_postConv
       << " " << ene_postDiffu
       << " " << ene_postSource
       << " delta=" << ene_postSource-ene_ini << " [W.m-3]." << finl;

  const IJK_Field_double& temperature = *temperature_;
  double Tmax = -1.e20;
  double Tmin = 1.e20;
  const int ni = temperature.ni();
  const int nj = temperature.nj();
  const int nk = temperature.nk();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          Tmax = std::max(Tmax, temperature(i,j,k));
          Tmin = std::min(Tmin, temperature(i,j,k));
        }
  Tmax = Process::mp_max(Tmax);
  Tmin = Process::mp_min(Tmin);
  Cerr <<"[Temperature-MinMax-" << rang_ <<"] t/Tmin/Tmax " << current_time << " "
       << Tmin << " " << Tmax
       << finl;
  return;
}

void IJK_Thermal_base::post_process_after_temperature_increment()
{
  if (debug_)
    Cerr << "Post-Processing routines for visu" << finl;
  Cerr << "Override temperature cell centres for visu" << finl;
  compute_temperature_cell_centres(1);
  if (debug_)
    Cerr << "Enforce periodic boundary conditions values" << finl;
  enforce_periodic_temperature_boundary_value();
  if (debug_)
    Cerr << "Clip temperature values" << finl;
  clip_min_temperature_values();
  clip_max_temperature_values();
  if (debug_)
    Cerr << "Correct temperature for visu" << finl;
  correct_temperature_for_visu();
  if (debug_)
    Cerr << "set analytical temperature field" << finl;
  set_field_T_ana();
  if (debug_)
    Cerr << "Correct operators in mixed cells for visu" << finl;
  correct_operators_for_visu();
}

void IJK_Thermal_base::compute_eulerian_grad_T_interface(const int on_splitting_ns)
{
  if (compute_grad_T_interface_)
    {
      temperature_ft_.data() = 0.;
      temperature_ft_.echange_espace_virtuel(temperature_ft_.ghost());
      temperature_->echange_espace_virtuel(temperature_->ghost());
      ref_ijk_ft_->eq_ns().redistribute_to_splitting_ft_elem(*temperature_, temperature_ft_);
      temperature_ft_.echange_espace_virtuel(temperature_ft_.ghost());

      compute_eulerian_normal_temperature_gradient_interface(*eulerian_distance_ft_,
                                                             ref_ijk_ft_->get_interface().I_ft(),
                                                             *eulerian_interfacial_area_ft_,
                                                             *eulerian_curvature_ft_,
                                                             temperature_ft_,
                                                             eulerian_grad_T_interface_ft_,
                                                             spherical_approx_);
      if (on_splitting_ns)
        {
          compute_eulerian_normal_temperature_gradient_interface(*eulerian_distance_ns_,
                                                                 ref_ijk_ft_->get_interface().I(),
                                                                 *eulerian_interfacial_area_ns_,
                                                                 *eulerian_curvature_ns_,
                                                                 *temperature_,
                                                                 eulerian_grad_T_interface_ns_,
                                                                 spherical_approx_);
        }
      else
        ref_ijk_ft_->eq_ns().redistribute_from_splitting_ft_elem(eulerian_grad_T_interface_ft_, eulerian_grad_T_interface_ns_);

    }
  else
    Cerr << "Don't compute the grad_T_interface field" << finl;
}

void IJK_Thermal_base::propagate_eulerian_grad_T_interface()
{
  if (compute_grad_T_interface_)
    {

      eulerian_grad_T_interface_ft_.echange_espace_virtuel(eulerian_grad_T_interface_ft_.ghost());
      propagate_eulerian_normal_temperature_gradient_interface(ref_ijk_ft_->get_interface(),
                                                               *eulerian_distance_ft_,
                                                               eulerian_grad_T_interface_ft_,
                                                               n_iter_distance_,
                                                               gfm_recompute_field_ini_,
                                                               gfm_zero_neighbour_value_mean_,
                                                               gfm_vapour_mixed_only_,
                                                               gfm_smooth_factor_);
      eulerian_grad_T_interface_ft_.echange_espace_virtuel(eulerian_grad_T_interface_ft_.ghost());
      ref_ijk_ft_->eq_ns().redistribute_from_splitting_ft_elem(eulerian_grad_T_interface_ft_, eulerian_grad_T_interface_ns_);
      eulerian_grad_T_interface_ns_.echange_espace_virtuel(eulerian_grad_T_interface_ns_.ghost());
    }
  else
    Cerr << "Don't compute the grad_T_interface field" << finl;
}

void IJK_Thermal_base::compute_eulerian_temperature_ghost(const int on_splitting_ns)
{
  if (ghost_fluid_)
    {

      temperature_ft_.echange_espace_virtuel(temperature_ft_.ghost());
      compute_eulerian_extended_temperature(ref_ijk_ft_->get_interface().I_ft(),
                                            *eulerian_distance_ft_,
                                            *eulerian_curvature_ft_,
                                            eulerian_grad_T_interface_ft_,
                                            temperature_ft_,
                                            spherical_approx_);
      if (on_splitting_ns)
        {
          eulerian_grad_T_interface_ns_.echange_espace_virtuel(eulerian_grad_T_interface_ns_.ghost());
          temperature_->echange_espace_virtuel(temperature_->ghost());
          compute_eulerian_extended_temperature(ref_ijk_ft_->get_interface().I(),
                                                *eulerian_distance_ns_,
                                                *eulerian_curvature_ns_,
                                                eulerian_grad_T_interface_ns_,
                                                *temperature_,
                                                spherical_approx_);
        }
      else
        {
          ref_ijk_ft_->eq_ns().redistribute_from_splitting_ft_elem(eulerian_grad_T_interface_ft_, eulerian_grad_T_interface_ns_);
          ref_ijk_ft_->eq_ns().redistribute_from_splitting_ft_elem(temperature_ft_, *temperature_);
        }
      temperature_->echange_espace_virtuel(temperature_->ghost());
    }
  else
    Cerr << "Don't compute the ghost temperature field" << finl;
}

void IJK_Thermal_base::compute_eulerian_bounding_box_fill_compo()
{
  if (compute_eulerian_compo_)
    {
      bubbles_barycentre_ = &(ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_bubbles_barycentre());
      bubbles_volume_ = &(ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_bubbles_volume());
      bounding_box_= &(ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_bounding_box());
      min_max_larger_box_ = &(ref_ijk_ft_->get_interface().get_ijk_compo_connex().get_min_max_larger_box());
    }
  else
    Cerr << "Don't compute the eulerian bubbles' components (composantes connexes)" << finl;
}

void IJK_Thermal_base::compute_temperature_gradient_elem()
{
  if (compute_grad_T_elem_)
    {
      temperature_->echange_espace_virtuel(temperature_->ghost());
      for (int dir=0; dir<3; dir++)
        grad_T_elem_[dir].data()=0;
      grad_T_elem_.echange_espace_virtuel();
      temperature_grad_op_centre_.calculer_grad(*temperature_, grad_T_elem_);
      grad_T_elem_.echange_espace_virtuel();

      if (smooth_grad_T_elem_)
        {
          // temperature_gaussian_filtered_.data() = temperature_->data();
          smooth_eulerian_field(temperature_gaussian_filtered_,
                                *temperature_,
                                0,
                                grad_T_elem_,
                                eulerian_normal_vectors_ns_normed_,
                                ref_ijk_ft_->get_interface(),
                                direct_smoothing_factors_,
                                gaussian_smoothing_factors_,
                                smoothing_numbers_, 0, 1, 0);
          temperature_gaussian_filtered_.echange_espace_virtuel(temperature_gaussian_filtered_.ghost());

          // grad_T_elem_smooth_ = grad_T_elem_;
          // grad_T_elem_smooth_.echange_espace_virtuel();
          smooth_vector_field(grad_T_elem_smooth_,
                              grad_T_elem_,
                              eulerian_normal_vectors_ns_normed_,
                              ref_ijk_ft_->get_interface(),
                              direct_smoothing_factors_,
                              gaussian_smoothing_factors_,
                              smoothing_numbers_,
                              smoothing_remove_normal_compo_, 0, 1, smoothing_use_unique_phase_);
          grad_T_elem_smooth_.echange_espace_virtuel();

          fill_tangential_gradient(grad_T_elem_smooth_,
                                   eulerian_normal_vectors_ns_normed_,
                                   grad_T_elem_tangential_);
          grad_T_elem_tangential_.echange_espace_virtuel();
        }
    }
  else
    Cerr << "The temperature gradient at the cell centres is not computed" << finl;
}

void IJK_Thermal_base::compute_temperature_hessian_diag_elem()
{
  if (compute_hess_diag_T_elem_)
    {
      temperature_->echange_espace_virtuel(temperature_->ghost());
      for (int dir=0; dir<3; dir++)
        hess_diag_T_elem_[dir].data()=0;
      hess_diag_T_elem_.echange_espace_virtuel();
      temperature_hess_op_centre_.calculer_hess(*temperature_, hess_diag_T_elem_,
                                                boundary_flux_kmin_, boundary_flux_kmax_);
      hess_diag_T_elem_.echange_espace_virtuel();
    }
  else
    Cerr << "The temperature gradient at the cell centres is not computed" << finl;
}

void IJK_Thermal_base::compute_temperature_hessian_cross_elem()
{
  if (compute_hess_cross_T_elem_)
    {
      temperature_->echange_espace_virtuel(temperature_->ghost());
      for (int dir=0; dir<3; dir++)
        hess_cross_T_elem_[dir].data()=0;
      hess_cross_T_elem_.echange_espace_virtuel();
      temperature_grad_op_centre_.calculer_grad_z(grad_T_elem_[1], hess_cross_T_elem_[0]);
      temperature_grad_op_centre_.calculer_grad_z(grad_T_elem_[0], hess_cross_T_elem_[1]);
      temperature_grad_op_centre_.calculer_grad_y(grad_T_elem_[0], hess_cross_T_elem_[2]);
      hess_cross_T_elem_.echange_espace_virtuel();
    }
  else
    Cerr << "The temperature gradient at the cell centres is not computed" << finl;
}

void IJK_Thermal_base::compute_temperature_convective_fluxes(const IJK_Field_vector3_double& velocity)
{
  static Stat_Counter_Id cnt_convective_flux_balance = statistiques().new_counter(2, "Thermal flux balance - Convection");
  static Stat_Counter_Id cnt_convective_flux_balance_op = statistiques().new_counter(3, "Thermal flux balance - Convection operator");
  static Stat_Counter_Id cnt_convective_flux_balance_factor = statistiques().new_counter(3, "Thermal flux balance - Convection prefactor");

  statistiques().begin_count(cnt_convective_flux_balance);

  if (store_flux_operators_for_energy_balance_)
    {
      for (int c=0; c<3; c++)
        rho_cp_u_T_convective_raw_[c].data() = 0.;
      if (!disable_relative_velocity_energy_balance_)
        {
          const Vecteur3 bubbles_velocity = (*rising_velocity_overall_);
          temperature_grad_flux_op_quick_.set_velocity_frame_of_reference(bubbles_velocity);
        }

      statistiques().begin_count(cnt_convective_flux_balance_op);
      temperature_grad_flux_op_quick_.calculer_grad_flux(*temperature_,
                                                         velocity[0],
                                                         velocity[1],
                                                         velocity[2],
                                                         rho_cp_u_T_convective_raw_);
      statistiques().end_count(cnt_convective_flux_balance_op);

      statistiques().begin_count(cnt_convective_flux_balance_factor);
      const double rhocp = get_rhocp_l();
      const int ni = d_temperature_->ni();
      const int nj = d_temperature_->nj();
      const int nk = d_temperature_->nk();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            if (store_flux_operators_for_energy_balance_)
              for (int c=0; c<3; c++)
                {
                  const double convective_flux = rho_cp_u_T_convective_raw_[c](i,j,k);
                  rho_cp_u_T_convective_raw_[c](i,j,k) = convective_flux * rhocp;
                }
      statistiques().end_count(cnt_convective_flux_balance_factor);
      rho_cp_u_T_convective_raw_.echange_espace_virtuel();
    }

  statistiques().end_count(cnt_convective_flux_balance);
}

// Convect temperature field by the velocity.
// The output is stored in *d_temperature_ (it is a volume integral over the CV)
void IJK_Thermal_base::compute_temperature_convection(const IJK_Field_vector3_double& velocity)
{
  static Stat_Counter_Id cnt_conv_temp = statistiques().new_counter(2, "FT convection temperature");
  static Stat_Counter_Id cnt_conv_temp_op = statistiques().new_counter(3, "FT convection temperature operator");
  static Stat_Counter_Id cnt_conv_temp_factor = statistiques().new_counter(3, "FT convection temperature factor");

  statistiques().begin_count(cnt_conv_temp);

  IJK_Field_double& d_temperature = *d_temperature_;

  if (conv_temperature_negligible_)
    {
      d_temperature.data()=0;
      u_T_convective_volume_.data() = 0;
    }
  else
    {
      statistiques().begin_count(cnt_conv_temp_op);
      temperature_convection_op_->calculer(*temperature_, velocity[0], velocity[1], velocity[2], d_temperature);
      statistiques().end_count(cnt_conv_temp_op);

      statistiques().begin_count(cnt_conv_temp_factor);
      const int post_pro_u_T_convective = liste_post_instantanes_.contient_("U_T_CONVECTIVE");
      const int ni = d_temperature.ni();
      const int nj = d_temperature.nj();
      const int nk = d_temperature.nk();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              u_T_convective_volume_(i,j,k) = d_temperature(i,j,k);
              const double resu = d_temperature(i,j,k) / vol_;
              d_temperature(i,j,k) = resu;

              if (post_pro_u_T_convective)
                u_T_convective_(i,j,k) = resu;
            }
      statistiques().end_count(cnt_conv_temp_factor);
    }

  statistiques().end_count(cnt_conv_temp);
  DebogIJK::verifier("op_conv(rho)", d_temperature);
  return;
}

void IJK_Thermal_base::compute_boundary_conditions_thermal()
{
  if (boundary_conditions_.get_bctype_k_min() == Boundary_Conditions_Thermique::Paroi_Temperature_imposee)
    {
      const double T_paroi_impose_kmin = boundary_conditions_.get_temperature_kmin();
      double lambda_de_t_paroi_kmin = lambda_liquid_;
      // calculer ici div(lambda*grad(T))*volume)
      calculer_flux_thermique_bord(*temperature_, lambda_de_t_paroi_kmin,
                                   T_paroi_impose_kmin, boundary_flux_kmin_, 0 /* boundary kmin */);
    }
  else if (boundary_conditions_.get_bctype_k_min() == Boundary_Conditions_Thermique::Paroi_Flux_impose)
    {
      const double flux_paroi_impose_kmin = boundary_conditions_.get_flux_kmin();
      imposer_flux_thermique_bord(*temperature_,
                                  flux_paroi_impose_kmin, boundary_flux_kmin_, 0 /* boundary kmin */);
    }
  else
    {
      // Perio... not needed. The flux on the "boundary" (in that case , it's not truely a boundary) will be computed as inside...
    }
  if (boundary_conditions_.get_bctype_k_max() == Boundary_Conditions_Thermique::Paroi_Temperature_imposee)
    {
      const double T_paroi_impose_kmax = boundary_conditions_.get_temperature_kmax();
      double lambda_de_t_paroi_kmax = lambda_liquid_;
      //TODO: calculer ici div(lambda*grad(T))*volume)
      calculer_flux_thermique_bord(*temperature_, lambda_de_t_paroi_kmax,
                                   T_paroi_impose_kmax, boundary_flux_kmax_, 1 /* boundary kmax */);
    }
  else if (boundary_conditions_.get_bctype_k_max() == Boundary_Conditions_Thermique::Paroi_Flux_impose)
    {
      const double flux_paroi_impose_kmax = boundary_conditions_.get_flux_kmax();
      imposer_flux_thermique_bord(*temperature_,
                                  flux_paroi_impose_kmax, boundary_flux_kmax_, 1 /* boundary kmax */);
      // Cerr << "not coded yet" << finl;
      // Process::exit();
    }
  else
    {
      // Perio... not needed. The flux on the "boundary" (in that case , it's not truely a boundary) will be computed as inside...
    }
  temperature_->echange_espace_virtuel(temperature_->ghost());
  DebogIJK::verifier("temp", *temperature_);
}

void IJK_Thermal_base::compute_temperature_diffusive_fluxes()
{
  static Stat_Counter_Id cnt_diffusive_flux_balance = statistiques().new_counter(2, "Thermal flux balance - Diffusion");
  statistiques().begin_count(cnt_diffusive_flux_balance);

  if (store_flux_operators_for_energy_balance_)
    {
      for (int dir=0; dir<3; dir++)
        div_coeff_grad_T_raw_[dir].data()=0;
      temperature_hess_flux_op_centre_.calculer_hess_flux(*temperature_,
                                                          div_coeff_grad_T_raw_,
                                                          boundary_flux_kmin_,
                                                          boundary_flux_kmax_);
      div_coeff_grad_T_raw_.echange_espace_virtuel();
    }

  statistiques().end_count(cnt_diffusive_flux_balance);
}

void IJK_Thermal_base::add_temperature_diffusion()
{
  static Stat_Counter_Id cnt_diff_temp = statistiques().new_counter(2, "FT diffusion temperature");
  static Stat_Counter_Id cnt_diff_temp_op = statistiques().new_counter(3, "FT diffusion temperature operator");
  static Stat_Counter_Id cnt_diff_temp_factor = statistiques().new_counter(3, "FT diffusion temperature factor");
  statistiques().begin_count(cnt_diff_temp);

  compute_boundary_conditions_thermal();
  if (!diff_temperature_negligible_)
    {
      /*
       * Correct the diffusive fluxes here or in the operator ?
       */
      statistiques().begin_count(cnt_diff_temp_op);
      temperature_diffusion_op_->calculer(*temperature_,
                                          *div_coeff_grad_T_volume_,
                                          boundary_flux_kmin_,
                                          boundary_flux_kmax_);
      statistiques().end_count(cnt_diff_temp_op);

      statistiques().begin_count(cnt_diff_temp_factor);
      compute_diffusion_increment();
      statistiques().end_count(cnt_diff_temp_factor);

      DebogIJK::verifier("div_coeff_grad_T_volume_", *div_coeff_grad_T_volume_);
    }

  statistiques().end_count(cnt_diff_temp);
}

//////////////////////////////////////////
//void IJK_Thermal_base::add_temperature_source() { ; }
double IJK_Thermal_base::compute_rho_cp_u_mean(const IJK_Field_double& vx)
{
  /*
   * By default use only the liquid phase (same for subresolution)
   * Overridden in Onefluid and others
   */
  const double rho_cp = ref_ijk_ft_->milieu_ijk().get_rho_liquid() * cp_liquid_;
  return calculer_rho_cp_u_moyen(vx, vx, vx, rho_cp, 2);
}

double IJK_Thermal_base::get_rho_cp_u_ijk(const IJK_Field_double& vx, int i, int j, int k) const
{
  return ref_ijk_ft_->milieu_ijk().get_rho_liquid() * cp_liquid_ * vx(i,j,k);
}

void IJK_Thermal_base::add_temperature_source()
{
  static Stat_Counter_Id cnt_source_temp = statistiques().new_counter(2, "FT source temperature");
  statistiques().begin_count(cnt_source_temp);
  // Dans le cas ou les flux entrants et sortants sont identiques :
  // DONE: changer cette condition non adaptee
  if (type_T_source_!="??")
    {
      const IJK_Field_double& temperature = *temperature_;
      IJK_Field_double& d_temperature     = *d_temperature_;

      // int gravity_dir = ref_ijk_ft_->milieu_ijk().get_direction_gravite();
      const int gravity_dir=0;
      /*
       * Modifications un jour peut-être ?!
       * Adaptation paroi dans une autre direction ?
       */
      const int wall_normal_dir = DIRECTION_K;
      const IJK_Field_double& vx = ref_ijk_ft_->eq_ns().get_velocity()[gravity_dir];
      double rho_cp_u_moy = compute_rho_cp_u_mean(vx);
      const Domaine_IJK& geom = temperature.get_domaine();
      const double dl = geom.get_constant_delta(gravity_dir);
      const double lwall = geom.get_domain_length(wall_normal_dir) ;
      const double h = lwall / 2.;
      const int nk = d_temperature_->nk();
      const int ni = d_temperature_->ni();
      const int nj = d_temperature_->nj();
      // TODO: faire une methode calculer_rho_cp
      //debut if source = ponderee par la vitesse
      if (type_T_source_=="dabiri")
        {
          const double wall_flux = boundary_conditions_.get_flux_kmax();
          // wall_flux in W.m^{-2}
          const double qw = wall_flux;
          // TODO: Ask Aymeric -> 2 wall with same B.Cs ?
          const double dTm = 2 * qw / rho_cp_u_moy;
          /*
           * Local correction using global dTm
           * TODO: faux, il manque un facteur 1/h pour homogeneite
           * MG: Pas sûr...
           */
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
//                  const double rho = ref_ijk_ft_->rho_field_(i,j,k);
//                  const double cp = cp_(i,j,k);
//                  const double u = (vx(i,j,k) +vx(i+1,j,k))/2;
//                  const double rho_cp_u = rho*cp*u;
                  /*
                   * TODO: Ask Aymeric : WTF lambda_variable_; at the wall ???
                   */
//                  if(lambda_variable_)
//                    {
//											const double div_lambda = compute_lambda_variations(dl);
//                      const double div_lambda = (lambda_(i+1,j,k)-lambda_(i-1,j,k))/(2*dl);
//                      source_temperature_(i,j,k) = (rho_cp_u - div_lambda) * dTm;
//                    }
//                  else
//                    {
//                      source_temperature_(i,j,k) = rho_cp_u * dTm;
//                    }
                  const double rho_cp_u = get_rho_cp_u_ijk(vx, i, j, k);
                  const double div_lambda = get_div_lambda_ijk(i,j,k) / (2*dl);
                  source_temperature_(i,j,k) = (rho_cp_u - div_lambda) * dTm;
                  // TODO: Each plate generates heat in a volume of half the plates distance ??
                  const double Sc = qw / h ;  // qw / lambda * h;
                  const double source = (source_temperature_(i,j,k) - Sc) * vol_; // -Sc ) * volume;
                  // TODO: faux, que vient faire Sc en plus de source ?
                  d_temperature(i,j,k) += source;
                }
          //
          calculer_temperature_physique_T(vx, dTm);
          calculer_temperature_adimensionnelle_theta(vx, wall_flux);
          calculer_Nusselt(vx);
          return;
        }
      // TODO: remplacer dabiri par patch_dabiri apres verif
      else if (type_T_source_=="patch_dabiri")
        {
          Cerr << "Type de source : patch_dabiri" << finl;
          const double wall_flux = boundary_conditions_.get_flux_kmax();
          const double qw = wall_flux;
          const double dTm = -2*qw/(2*h*rho_cp_u_moy) ;
          if (std::fabs(rho_cp_u_moy)<DMINFLOAT)
            {
              Cerr << "You cannot use this source " << type_T_source_ << " without flowrate!" << finl;
              Process::exit();
            }

          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  const double u = (vx(i,j,k) +vx(i+1,j,k))/2;
                  source_temperature_(i,j,k) = (u - get_div_lambda_ijk(i,j,k) / (2*dl)) * dTm;
                  const double source = source_temperature_(i,j,k);
                  d_temperature(i,j,k) += source;
                }
          //
          calculer_temperature_physique_T(vx, dTm);
          calculer_temperature_adimensionnelle_theta(vx, wall_flux);
          calculer_Nusselt(vx);
          return;
        }
      else if (type_T_source_=="unweighted_dabiri")
        {
          // Sans la ponderation de S par u (S=rhocpu.../<rhocpu>), cela permet d'avoir une source uniforme!
          // Mais ce n'est plus le meme changement de variable, c'est quoi alors?
          // TODO: Que faire de rho_cp en diphasique? Moyen ou local?
          Cerr << "Type de source : unweighted_dabiri" << finl;
          const double wall_flux = boundary_conditions_.get_flux_kmax();
          const double qw = wall_flux;
          const double liquid_fraction = calculer_v_moyen(ref_ijk_ft_->get_interface().I());
          const double rho_cp_l = ref_ijk_ft_->milieu_ijk().get_rho_liquid() * cp_liquid_;
          const double rho_cp_v = ref_ijk_ft_->milieu_ijk().get_rho_vapour() * cp_vapour_;
          const double rhocp_moy =  rho_cp_l*liquid_fraction + rho_cp_v*(1-liquid_fraction);
          const double dTm = -2*qw/(2*h*rhocp_moy) ;
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  if(lambda_variable_)
                    {
                      Cerr << "Veut-on vraiment calculer une partie en lambda spatialement variable??" <<finl;
                      Cerr << "Exit at IJK_Thermal_base::add_temperature_source" << finl;
                      Process::exit();
                    }
                  else
                    {
                      source_temperature_(i,j,k) = dTm;
                    }
                  const double source = source_temperature_(i,j,k);
                  d_temperature(i,j,k) += source;
                }
          //
          calculer_temperature_physique_T(vx, dTm);
          calculer_temperature_adimensionnelle_theta(vx, wall_flux);
          calculer_Nusselt(vx);
          return;
        }
      // Debut source = SWARM
      else if (type_T_source_=="SWARM")
        {
          // DONE: idem
          double Tv=0;
          double Tl=0;
          double Vv = 0;
          double Vl = 0;

          // calcul de Tv et Tl
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  const double chi = ref_ijk_ft_->get_interface().I(i, j, k);
                  const double T = temperature(i, j, k);
                  Tv += T*(1.-chi);
                  Vv += (1.-chi);
                  Tl += T*chi;
                  Vl += chi;
                }
          Tv = Process::mp_sum(Tv);
          Tl = Process::mp_sum(Tl);
          Vv = Process::mp_sum(Vv);
          Vl = Process::mp_sum(Vl);
          Tv /= Vv;
          Tl /= Vl;
          Cerr << "AY-test_source : " <<  Tv << finl;

          // Calcul de dT a partir de l'expression suivante : dT = k(T - Tm)/(rho*cp)
          const double kv = kv_;
          const double kl = kl_;
          const double T0v = T0v_;
          const double T0l = T0l_;
          const double rho_cp_l = ref_ijk_ft_->milieu_ijk().get_rho_liquid() * cp_liquid_;
          const double rho_cp_v = ref_ijk_ft_->milieu_ijk().get_rho_vapour() * cp_vapour_;
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  d_source_Tv_(i,j,k) = kv/rho_cp_v * (Tv - T0v);
                  d_source_Tl_(i,j,k) = kl/rho_cp_l * (Tl - T0l);
                }
          Cerr << "AY-test_source1 : " <<  Tv << finl;
          // TODO: Remplacer euler_explicit_update par l'utilisation de timestep_ et utiliser d_source_Tv_ comme une constante
          for (int k = 0; k < nk; k++)
            {
              ref_ijk_ft_->eq_ns().euler_explicit_update(d_source_Tv_, source_temperature_v_, k);
              ref_ijk_ft_->eq_ns().euler_explicit_update(d_source_Tl_, source_temperature_l_, k);
            }
          Cerr << "AY-test_source2 : " <<  Tv << finl;
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  // chi vaut 1. dans le liquide et 0 dans les bulles
                  const double chi = ref_ijk_ft_->get_interface().I(i, j, k);
                  source_temperature_(i,j,k) = (1.-chi)*source_temperature_v_(i,j,k) + chi*source_temperature_l_(i,j,k);
                }
          Cerr << "AY-test_source3 : " <<  Tv << finl;
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  d_temperature(i,j,k) += source_temperature_(i,j,k)*vol_;
                }
          Cerr << "AY-test_source 111 : " <<  source_temperature_(1,1,1) << finl;
          Cerr << "AY-test_source vol : " <<  vol_ << finl;
          Cerr << "source_temp " << " " << d_temperature(1,1,1) << finl;
          const double current_time = ref_ijk_ft_->schema_temps_ijk().get_current_time();

          // GB : Ma comprehension est que ce ne sont pas des champs, mais un scalaire unique
          const double Sl = source_temperature_l_(0,0,0);
          const double Sv = source_temperature_v_(0,0,0);
          const double dSl = d_source_Tl_(0,0,0);
          const double dSv = d_source_Tv_(0,0,0);
          Cerr <<"[ThermalInfo-" << rang_ <<"] t/Tl/Tv/Sl/Sv/dSldt/dSvdt " << current_time << " "
               << Tl << " " << Tv << " "
               << Sl << " " << Sv << " "
               << dSl << " " << dSv
               <<finl;
          /*
           * TODO: M.G c'est quoi ??
           * Si on utilise ca c'est pour remplir dans tous les cas d'utilisation d'une source le champs temperature_physique_T
           * calculer_temperature_physique_T_dummy();
           */
          temperature_physique_T_.data() = 0.;
          statistiques().end_count(cnt_source_temp);
          return;
        }
    }
  // Dans ce cas la ce ne sont pas des flux thermiques identiques
  else
    {
      Cerr << "no_source_for_temperature" << finl;
      statistiques().end_count(cnt_source_temp);
      return;
    }
}

void IJK_Thermal_base::source_callback()
{
  if (debug_)
    Cerr << "Source terms post-processing routines" << finl;
  if (liste_post_instantanes_.contient_("TEMPERATURE_ADIM_BULLES"))
    calculer_temperature_adim_bulles();
}

/*
 * Aymeric: Renommer pour expliciter qu'il s'agit de la transformation
 * inverse de Kawamura avec le gradient de temperature moyenne
 */
void IJK_Thermal_base::calculer_temperature_physique_T(const IJK_Field_double&  vx, const double dTm)
{
  if (wall_flux_)
    {
      const IJK_Field_double& temperature = *temperature_;

      const Domaine_IJK& geom = temperature.get_domaine();
      double dx =geom.get_constant_delta(DIRECTION_I);
      double origin_x = geom.get_origin(DIRECTION_I) + (dx * 0.5) ;
      const int offset_i = geom.get_offset_local(DIRECTION_I);

      const int nk = d_temperature_->nk();
      const int ni = d_temperature_->ni();
      const int nj = d_temperature_->nj();

      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double x = (i+ offset_i ) *dx + origin_x;    //MR A GB: ne doit-on pas soustraire par  -dx*0.5???
              temperature_physique_T_(i,j,k) = (x*dTm)-temperature(i,j,k);
            }
      temperature_physique_T_.echange_espace_virtuel(temperature_physique_T_.ghost());
      DebogIJK::verifier("temperature_physique_T", temperature_physique_T_);
      return;
    }
  else
    {
      Cerr << "No source for the temperature field" << finl;
      return;
    }
}

void IJK_Thermal_base::calculer_temperature_adim_bulles()
{
  const IJK_Field_double& temperature = *temperature_;

  const int nk = temperature.nk();
  const int ni = temperature.ni();
  const int nj = temperature.nj();

  // Calcul de moy1 = moy(chi*T+)/moy(chi) et moy2 = moy((1-chi)*T+) / moy(1-chi)
  double Tl = 0.;
  double Tv = 0.;
  double Vl = 0.;
  double Vv = 0.;
  const IJK_Field_double& chi = ref_ijk_ft_->get_interface().I(); // rappel : chi vaut 1. dans le liquide et 0 dans la vapeur

  // assuming uniform mesh : vol_cell=cste.
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          Tl += temperature(i,j,k)*chi(i, j, k);
          Tv += temperature(i,j,k)*(1.-chi(i, j, k));
          Vl += chi(i, j, k);
          Vv += (1.-chi(i, j, k));
        }
  Tl = Process::mp_sum(Tl);
  Tv = Process::mp_sum(Tv);
  Vv = Process::mp_sum(Vv);
  Vl = Process::mp_sum(Vl);
  Tl /= Vl;
  Tv /= Vv;
  // Calcul de Tl et Tv :
  // const double Tl = Tv0_ / (1 - Tl) * (1 - Tl / (1 + Tv)) / (1 + Tl * Tv / ((1 - Tl) * (1 + Tv)));
  // const double Tv = Tv0_ / (1 + Tv) * (1 + Tv / (1 - Tl)) / (1 + Tv * Tl / ((1 + Tv) * (1 - Tl)));
  const int ntot = temperature.get_domaine().get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_I)
                   *temperature.get_domaine().get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_J)
                   *temperature.get_domaine().get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_K);
  const double time = ref_ijk_ft_->schema_temps_ijk().get_current_time();
  Cerr << "Tl_test : time= "<< time << " alpha_l= " << Vl/ntot<< "  TI=" <<Tl*Vl/ntot<< " Tl="<< Tl << finl;
  Cerr << "Tv_test : time= "<< time << " alpha_v= " << Vv/ntot<< " TIv=" <<Tv*Vv/ntot<< " Tv="<< Tv << finl;

  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        temperature_adim_bulles_(i,j,k) = (temperature(i,j,k) - Tv0_)/(Tl - Tv);

  //  temperature_adim_bulles_.echange_espace_virtuel(temperature_adim_bulles_.ghost());
  //  TODO: ce qui suit ne devrait pas etre la, mais je le met ici temporairement avant de trouver une meilleure solution

  double E_tot = 0.;
  double E_liq_pure = 0., E_liq = 0;
  double E_vap_pure = 0., E_vap = 0;
  double E_mixt = 0.;
  calculer_energies(E_liq_pure, E_liq,
                    E_vap_pure, E_vap,
                    E_mixt, E_tot);

  /*
   * TODO: voir si on ne doit pas faire mieux, mais a priori les variations de Tl et Tv
   * sont lentes par rapport au reste donc ea devrait aller.
   * DONE: il y a manifestement un pb ici, car on ne peut pas avoir acces e Tv(n+1) encore,
   * donc il faut stocker Tv(n-1)
   */

  // Impression dans le fichier temperature_bulles.out
  if (Process::je_suis_maitre())
    {
      int reset = (!ref_ijk_ft_->get_reprise()) && (ref_ijk_ft_->schema_temps_ijk().get_tstep()==0);
      SFichier fic=Ouvrir_fichier(Nom("_source_temperature_")+Nom(rang_)+Nom("_bulles.out"),
                                  "tstep\ttime\tTl\tTv\tEtot\tElpu\tEl\tEvpu\tEv\tEm",
                                  reset);
      // la derivee_acceleration n'est connue que sur le maitre
      fic<< ref_ijk_ft_->schema_temps_ijk().get_tstep()<<" "<< ref_ijk_ft_->schema_temps_ijk().get_current_time() <<" "<< Tl << " " << Tv << " " << E_tot;
      fic<< " " << E_liq_pure <<  " " << E_liq;
      fic<< " " << E_vap_pure <<  " " << E_vap;
      fic<< " " << E_mixt;
      fic<<finl;
      fic.close();
    }
}

void IJK_Thermal_base::calculer_energies(double& E_liq_pure,
                                         double& E_liq,
                                         double& E_vap_pure,
                                         double& E_vap,
                                         double& E_mixt,
                                         double& E_tot)
{
  const IJK_Field_double& temperature = *temperature_;

  const int nk = temperature.nk();
  const int ni = temperature.ni();
  const int nj = temperature.nj();
  const IJK_Field_double& chi = ref_ijk_ft_->get_interface().I(); // rappel : chi vaut 1. dans le liquide et 0 dans la vapeur
  const double rhocpl = ref_ijk_ft_->milieu_ijk().get_rho_liquid()*cp_liquid_;
  const double rhocpv = ref_ijk_ft_->milieu_ijk().get_rho_vapour()*cp_vapour_;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double indic=chi(i, j, k);
          const double  rhocpm = indic*rhocpl+(1.-indic)*rhocpv;
          const double T = temperature(i,j,k);
          E_tot += T * rhocpm;
          E_liq += indic * rhocpl * T;
          E_vap += (1.-indic) * rhocpv * T;
          if (std::fabs(indic)<1.e-8)
            {
              // vap pure
              E_vap_pure += rhocpv * T;
            }
          else if (std::fabs(1.-indic)<1.e-8)
            {
              // liq pure
              E_liq_pure += rhocpl * T;
            }
          else
            {
              // mixte :
              E_mixt += rhocpm * T;
            }
        }

  const int ntot = (temperature.get_domaine().get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_I)
                    * temperature.get_domaine().get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_J)
                    * temperature.get_domaine().get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_K));
  E_vap_pure = Process::mp_sum(E_vap_pure)/ntot;
  E_liq_pure = Process::mp_sum(E_liq_pure)/ntot;
  E_vap = Process::mp_sum(E_vap)/ntot;
  E_liq = Process::mp_sum(E_liq)/ntot;
  E_tot = Process::mp_sum(E_tot)/ntot;
  E_mixt = Process::mp_sum(E_mixt)/ntot;
}

void IJK_Thermal_base::calculer_source_temperature_ana()
{
  if (liste_post_instantanes_.contient_("ECART_SOURCE_TEMPERATURE_ANA"))
    {
      if (!liste_post_instantanes_.contient_("SOURCE_TEMPERATURE_ANA"))
        set_field_data(source_temperature_ana_, expression_source_temperature_, ref_ijk_ft_->eq_ns().get_velocity()[0], ref_ijk_ft_->schema_temps_ijk().get_current_time());
      // do some work

      double ct = ref_ijk_ft_->schema_temps_ijk().get_current_time();
      Cerr << "MR: ERROR SOURCE T FIELD " << ct;
      double err = 0.;
      //set_field_data(source_temperature_ana_, curseur->expression_source_T_ana_, ct);
      const int ni = source_temperature_.ni();
      const int nj = source_temperature_.nj();
      const int nk = source_temperature_.nk();
      const trustIdType ntot=Process::mp_sum(ni*nj*nk);
      // La temperature est definie a une constante pres:
      // const double cst_temp = temperature_ana_(0,0,0) - curseur->temperature_(0,0,0);
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double val = source_temperature_ana_(i,j,k) - source_temperature_(i,j,k); //- cst_temp;
              ecart_source_t_ana_(i,j,k) = val;
              err += val*val;
            }
      err=Process::mp_sum(err);
      err=sqrt(err/static_cast<double>(ntot));
      Cerr << " " << err ;
      if (!Process::je_suis_maitre())
        {
          Process::Journal() << "IJK_FT_Post::posttraiter_champs_instantanes : OWN_PTR(Champ_base) ECART_SOURCE_TEMPERATURE_ANA sur ce proc (ni,nj,nk,ntot):"
                             << " " << ni << " " << nj << " " << nk << " " << ntot << finl;
        }
      // ecart_t_ana_.echange_espace_virtuel(ecart_t_ana_.ghost());
      Cerr << finl ;
      // n++,dumplata_scalar(lata_name,"ECART_SOURCE_TEMPERATURE_ANA", ecart_source_t_ana_, latastep);
    }
}

double IJK_Thermal_base::compute_variable_wall_temperature(const int kmin, const int kmax)
{
  return calculer_variable_wall(*temperature_, *temperature_, *temperature_, cp_liquid_ * ref_ijk_ft_->milieu_ijk().get_rho_liquid(), kmin, kmax, 2);
}

void IJK_Thermal_base::calculer_temperature_adimensionnelle_theta(const IJK_Field_double& vx, const double wall_flux)
{
  /*
   * TODO : M.G -> Don't understand anything ask Aymeric
   */
  if(wall_flux_)
    {
      const IJK_Field_double& temperature = *temperature_;

      const int wall_normal_dir = DIRECTION_K;
      const Domaine_IJK& geom = temperature.get_domaine();
      const double lwall = geom.get_domain_length(wall_normal_dir);
      const double h = lwall/2.;
      const double q_w = wall_flux;
      const int kmin = temperature.get_domaine().get_offset_local(wall_normal_dir);
      const int kmax = geom.get_nb_elem_tot(wall_normal_dir);
      const int nk = temperature.nk();
      const int ni = temperature.ni();
      const int nj = temperature.nj();
      double T_wall = compute_variable_wall_temperature(kmin, kmax);
      /*   if (Process::je_suis_maitre())
           {
           T_wall = calculer_variable_wall(temperature, cp_, ref_ijk_ft_->rho_field_, kmin, kmax);
           Cerr << "calcul de T_wall sur maitre" << finl;
           }
           envoyer_broadcast(T_wall, 0); */
      /*
       * if(kmin+ nk == kmax) //|| (kmin ==0)
       {
       rank = Process::me();
       T_wall = calculer_variable_wall(temperature, cp_, ref_ijk_ft_->rho_field_, kmin, kmax);
       envoyer_broadcast(T_wall, rank);
       }
       */
      for (int k = 0; k < nk; k++)
        {
          //  const double T_mean = compute_spatial_mean(vx, temperature, cp_, ref_ijk_ft_->rho_field_, kmin, nktot, k);
          for (int j = 0; j < nj; j++)
            for (int i = 0; i < ni; i++)
              {
                double theta = T_wall - temperature(i,j,k);
                // double theta = T_wall -T_mean;
                // temperature_adimensionnelle_theta_(i,j,k) = theta/theta_tau;
                const double lambda_l = lambda_liquid_;
                temperature_adimensionnelle_theta_(i,j,k) = theta * lambda_l / q_w / h;
                //    temperature_adimensionnelle_theta_(i,j,k) = -T_mean*lambda_l/q_w/h;
              }

        }
      return;
    }
  else
    {
      Cerr << "no_source_for_temperature" << finl;
      return;
    }
}

double IJK_Thermal_base::compute_temperature_dimensionless_theta_mean(const IJK_Field_double& vx)
{
  return calculer_temperature_adimensionnelle_theta_moy(vx, temperature_adimensionnelle_theta_, vx, vx, ref_ijk_ft_->milieu_ijk().get_rho_liquid() * cp_liquid_, 2);
}

void IJK_Thermal_base::calculer_Nusselt(const IJK_Field_double& vx)
{
  const double theta_adim_moy = compute_temperature_dimensionless_theta_mean(vx);
  double Nu = 0.;
  if (std::fabs(theta_adim_moy)>1.e-10)
    Nu = 2./theta_adim_moy;
  const double rho_cp_u_moy = compute_rho_cp_u_mean(vx);

  // Impression dans le fichier source_temperature.out
  if (Process::je_suis_maitre())
    {
      int reset = (!ref_ijk_ft_->get_reprise()) && (ref_ijk_ft_->schema_temps_ijk().get_tstep() == 0);
      const Nom name = Nom("_temperature_") + Nom(rang_) + Nom(".out");
      SFichier fic = Ouvrir_fichier(name,
                                    "tstep\ttime\ttheta_adim_moy\tNu\trho_cp_u",
                                    reset);
      fic << ref_ijk_ft_->schema_temps_ijk().get_tstep() << " " << ref_ijk_ft_->schema_temps_ijk().get_current_time() << " " << theta_adim_moy << " " << Nu << " " << rho_cp_u_moy  << finl;
      fic.close();
    }
}

void IJK_Thermal_base::set_field_T_ana()
{
  if (expression_T_ana_ != "??")
    {
      Cerr << "Setting analytical temperature "<< rang_ <<" field to "<< expression_T_ana_ << finl;
      set_field_data(temperature_ana_, expression_T_ana_, ref_ijk_ft_->schema_temps_ijk().get_current_time());
    }
}

void IJK_Thermal_base::calculer_ecart_T_ana()
{
  if (liste_post_instantanes_.contient_("ECART_T_ANA"))
    {
      const IJK_Field_double& temperature = *temperature_;

      if (!liste_post_instantanes_.contient_("TEMPERATURE_ANA"))
        {
          set_field_data(temperature_ana_, expression_T_ana_, ref_ijk_ft_->schema_temps_ijk().get_current_time());
        }
      // do some work

      double ct = ref_ijk_ft_->schema_temps_ijk().get_current_time();
      Cerr << "GB: ERROR T FIELD " << ct;
      double err = 0.;
      set_field_data(temperature_ana_, expression_T_ana_, ct);
      const int ni = temperature.ni();
      const int nj = temperature.nj();
      const int nk = temperature.nk();
      const trustIdType ntot=Process::mp_sum(ni*nj*nk);
      // La temperature est definie a une constante pres:
      // const double cst_temp = temperature_ana_(0,0,0) - curseur->temperature_(0,0,0);
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double val =   temperature_ana_(i,j,k) - temperature(i,j,k); //- cst_temp;
              ecart_t_ana_(i,j,k) = val;
              err += val*val;
            }
      err=Process::mp_sum(err);
      err=sqrt(err/static_cast<double>(ntot));
      Cerr << " " << err ;
      if (!Process::je_suis_maitre())
        {
          Process::Journal() << "IJK_FT_Post::posttraiter_champs_instantanes : OWN_PTR(Champ_base) ECART_T_ANA sur ce proc (ni,nj,nk,ntot):"
                             << " " << ni << " " << nj << " " << nk << " " << ntot << finl;
        }
      ecart_t_ana_.echange_espace_virtuel(ecart_t_ana_.ghost());
      Cerr << finl ;
      //  n++,dumplata_scalar(lata_name,"ECART_T_ANA", ecart_t_ana_, latastep);
    }
}

void IJK_Thermal_base::calculer_gradient_temperature(const IJK_Field_double& temperature, IJK_Field_vector3_double& grad_T)
{
  if ((liste_post_instantanes_.contient_("GRAD_T") || (calulate_grad_T_)))
    {
      /*
       * Re-initialisation of the gradient vector
       */
      for (int dir = 0; dir < 3; dir++)
        grad_T[dir].data() = 0.;

      //  add_gradient_temperature(temperature, 1. /*constant*/,  grad_T[0], grad_T[1], grad_T[2], boundary_conditions_, lambda_);
      for (int dir = 0; dir < 3; dir++)
        grad_T[dir].echange_espace_virtuel(1);
    }
}

// Results are intensive (ie prop to area)
// Method fills storage_ so it changes the class
// Les interfaces connaissent le splitting_ft_ donc la correspondance doit etre appliquee au splitting ft pour convertir :
// convert_packed_to_ijk_cell.
// Donc il faut un champ de T etendu...

double IJK_Thermal_base::compute_global_energy(const IJK_Field_double& temperature)
{
  double global_energy = 0.;
  const IJK_Field_double& indic = ref_ijk_ft_->get_interface().I();
  const double rhocpl = get_rhocp_l();
  const double rhocpv = get_rhocp_v();
  const int nx = temperature.ni();
  const int ny = temperature.nj();
  const int nz = temperature.nk();
  // To be sure we're on a regular mesh
  assert(indic.get_domaine().get_constant_delta(DIRECTION_K) > 0);
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        {
          double chi_l = indic(i,j,k);
          global_energy += (rhocpl * chi_l + (1.- chi_l) * rhocpv) * temperature(i,j,k);
        }
  const int ntot = temperature.get_domaine().get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_I)
                   * temperature.get_domaine().get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_J)
                   * temperature.get_domaine().get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_K);
  global_energy = mp_sum(global_energy)/(double)(ntot);
  return global_energy;
}

/*
 * Getters and setters
 */

double IJK_Thermal_base::get_rhocp_l() const
{
  return ref_ijk_ft_->milieu_ijk().get_rho_liquid()*ref_ijk_ft_->milieu_ijk().get_cp_liquid(0);
}

double IJK_Thermal_base::get_rhocp_v() const
{
  return ref_ijk_ft_->milieu_ijk().get_rho_vapour() * ref_ijk_ft_->milieu_ijk().get_cp_vapour(0);
}

/*
 * Methods that do not belong to the class
 */

// From DNS_QC; Vectorize code later?
int IJK_Thermal_base::calculer_k_pour_bord(const IJK_Field_double& temperature, const bool bord_kmax)
{
  const int kmin = temperature.get_domaine().get_offset_local(DIRECTION_K);
  const int nktot = temperature.get_domaine().get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_K);
  int k;
  // calcul l'indice k de la couche de mailles voisine du bord. Si je n'ai pas de bord, on met k = -1
  if (!bord_kmax)
    {
      // on veut le bord "k_global = 0"
      if (kmin == 0)
        {
          // ce bord est chez moi... et il est en k=0
          k = 0;
        }
      else
        {
          // ce bord n'est pas chez moi
          k = -1;
        }
    }
  else
    {
      // on veut le bord kmax
      if (kmin + temperature.nk() == nktot)
        {
          // ce bord est chez moi... et il est en k= truc...
          k = temperature.nk() - 1;
        }
      else
        {
          k = -1;
        }
    }
  return k;
}

// From DNS_QC; Vectorize code later?
// valeur de retour: indice local du plan de temperature voisin utilise,
//  -1 si on n'a pas le bord sur ce processeur
// Calcule l'integrale sur chaque face du bord demande du flux de chaleur a travers la face
// positif si le flux va vers les k positifs.
int IJK_Thermal_base::calculer_flux_thermique_bord(const IJK_Field_double& temperature,
                                                   const double lambda_de_t_paroi,
                                                   const double T_paroi_impose,
                                                   IJK_Field_local_double& flux_bord,
                                                   const bool bord_kmax)
{
  const int kmin = temperature.get_domaine().get_offset_local(DIRECTION_K);
  int k = calculer_k_pour_bord(temperature, bord_kmax);
  if (k == -1)
    return k;

  // redimensionne flux_bord avec ni * nj:
  const int ni = temperature.ni(); // nombre d'element local sur ce processeur
  const int nj = temperature.nj();
  flux_bord.allocate(ni, nj, 1 /* 1 seule couche de valeurs en k */, 0 /* pas d'elements fantomes */);

  const Domaine_IJK& geometry = temperature.get_domaine();
  const double delta_k = geometry.get_delta(DIRECTION_K)[k + kmin]; // k+kmin est l'indice global de la maille locale k
  double facteur = 2.0 / delta_k * geometry.get_constant_delta(DIRECTION_I) * geometry.get_constant_delta(DIRECTION_J);
  if (bord_kmax)
    facteur *= -1.; // he he... je vous laisse reflechir a ca :)
  // nan c'est pas simpa: la convention dans l'operateur de diffusion est
  // d/dt = flux(i,j) - flux(i+1,j) + ... + flux(i,j,k) - flux(i,j,k+1)
  // donc si la paroi inferieure (k=0) est plus froide que le fluide, il faut que le flux stocke soit negatif.
  // et si la paroi inferieure (k=kmax) est plus froide que le fluide, il faut que le flux stocke soit positif.
  for (int j = 0; j < nj; j++)
    {
      for (int i = 0; i < ni; i++)
        {
          // Temperature de la maille voisine
          const double t = temperature(i,j,k);
          // le flux est positif s'il va vers les k croissants
          flux_bord(i,j,0) = (T_paroi_impose - t) * lambda_de_t_paroi * facteur;
        }
    }
  return k;
}
int IJK_Thermal_base::imposer_flux_thermique_bord(const IJK_Field_double& temperature,
                                                  const double flux_paroi_impose,
                                                  IJK_Field_local_double& flux_bord,
                                                  const bool bord_kmax)
{
  int k = calculer_k_pour_bord(temperature, bord_kmax);
  if (k == -1)
    return k;

  // redimensionne flux_bord avec ni * nj:
  const int ni = temperature.ni(); // nombre d'element local sur ce processeur
  const int nj = temperature.nj();
  flux_bord.allocate(ni, nj, 1 /* 1 seule couche de valeurs en k */, 0 /* pas d'elements fantomes */);
  // MR je multiplie le flux par la surface dxdy
  const Domaine_IJK& geometry = temperature.get_domaine();

  double facteur = 1.* geometry.get_constant_delta(DIRECTION_I) * geometry.get_constant_delta(DIRECTION_J);
  if (bord_kmax)
    facteur *= -1.; // he he... je vous laisse reflechir a ca :)
  // nan c'est pas sympa: la convention dans l'operateur de diffusion est
  // d/dt = flux(i,j) - flux(i+1,j) + ... + flux(i,j,k) - flux(i,j,k+1)
  // donc si la paroi inferieure (k=0) est plus froide que le fluide, il faut que le flux stocke soit negatif.
  // et si la paroi superieure (k=kmax) est plus froide que le fluide, il faut que le flux stocke soit positif.
  for (int j = 0; j < nj; j++)
    {
      for (int i = 0; i < ni; i++)
        {
          // le flux est positif s'il va vers les k croissants
          flux_bord(i,j,0) = flux_paroi_impose * facteur;
        }
    }
  return k;
}

void IJK_Thermal_base::euler_rustine_step(const double timestep, const double dE)
{
  compute_dT_rustine(dE);
  // Update the temperature :
  const int kmax = temperature_->nk();
  const double ene_ini = compute_global_energy();
  for (int k = 0; k < kmax; k++)
    ref_ijk_ft_->eq_ns().euler_explicit_update(d_T_rustine_, *temperature_, k);
  temperature_->echange_espace_virtuel(temperature_->ghost());
  const double ene_post = compute_global_energy();
  Cerr << "[Energy-Budget-T"<<rang_<<" euler rustine] time t=" << ref_ijk_ft_->schema_temps_ijk().get_current_time()
       << " " << ene_ini
       << " " << ene_post << " [W.m-3]."
       << " dE "<< dE
       << finl;
  source_callback();
}

void IJK_Thermal_base::compute_dT_rustine(const double dE)
{
  const int ni = T_rust_.ni();
  const int nj = T_rust_.nj();
  const int nk = T_rust_.nk();
  const double rho_l = ref_ijk_ft_->milieu_ijk().get_rho_liquid();
  const double rho_v = ref_ijk_ft_->milieu_ijk().get_rho_vapour();
  const IJK_Field_double& indic = ref_ijk_ft_->get_interface().I();
  double int_rhocpTrust = 0;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          int_rhocpTrust +=  (rho_l*cp_liquid_*indic(i,j,k) + rho_v*cp_vapour_*(1.-indic(i,j,k)))*T_rust_(i,j,k);
        }
  const int ntot = T_rust_.get_domaine().get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_I)
                   *T_rust_.get_domaine().get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_J)
                   *T_rust_.get_domaine().get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_K);
  int_rhocpTrust = mp_sum(int_rhocpTrust)/(double)(ntot);
  Cerr << "Le coeff de manque d'energie dE/int_rhocpTrust vaut : " << dE/int_rhocpTrust << finl;
  if (int_rhocpTrust)
    {
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              d_T_rustine_(i,j,k) = dE/ int_rhocpTrust * T_rust_(i,j,k);
            }
    }
}

void IJK_Thermal_base::rk3_rustine_sub_step(const int rk_step, const double total_timestep,
                                            const double fractionnal_timestep, const double time, const double dE)
{
  compute_dT_rustine(dE);
  // Update the temperature :
  const int kmax = temperature_->nk();
  const double ene_ini = compute_global_energy();
  for (int k = 0; k < kmax; k++)
    {
      runge_kutta3_update(d_T_rustine_, RK3_F_rustine_, *temperature_, rk_step, k, total_timestep);
    }
  temperature_->echange_espace_virtuel(temperature_->ghost());
  const double ene_post = compute_global_energy();
  Cerr << "[Energy-Budget-T"<<rang_<<"RK3 rustine step "<<rk_step<<"] time t=" << ref_ijk_ft_->schema_temps_ijk().get_current_time()
       << " " << ene_ini
       << " " << ene_post << " [W.m-3]. [step"<< rk_step << "]" << finl;
  source_callback();
}

void IJK_Thermal_base::compute_interfacial_temperature2(ArrOfDouble& interfacial_temperature,
                                                        ArrOfDouble& flux_normal_interp)
{
  const Domaine_IJK& geom = ref_ijk_ft_->get_domaine();
  const double dist = 1.52 * std::pow(std::pow(geom.get_constant_delta(0), 2.) +
                                      std::pow(geom.get_constant_delta(1), 2.) +
                                      std::pow(geom.get_constant_delta(2), 2.),
                                      0.5);
  const Maillage_FT_IJK& maillage = ref_ijk_ft_->get_interface().maillage_ft_ijk();
  const DoubleTab& normale_facettes = maillage.get_update_normale_facettes();
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  const int nb_facettes = maillage.nb_facettes();
  const IntTab& facettes = maillage.facettes();
  const DoubleTab& sommets = maillage.sommets();
  DoubleTab coord_facettes;
  coord_facettes.resize(nb_facettes, 3);
  coord_facettes = 0.;
  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    for (int dir = 0; dir < 3; dir++)
      for (int som = 0; som < 3; som++)
        coord_facettes(fa7, dir) += sommets(facettes(fa7, som), dir);

  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    for (int dir = 0; dir < 3; dir++)
      coord_facettes(fa7, dir) /= 3.;

  DoubleTab coo_liqu, coo_vap;
  ArrOfDouble temp_liqu, temp_vap;
  corrige_flux_->calcul_temperature_flux_interface(temperature_ft_,
                                                   lambda_liquid_,
                                                   lambda_vapour_,
                                                   dist,
                                                   coord_facettes,
                                                   normale_facettes,
                                                   interfacial_temperature,
                                                   flux_normal_interp,
                                                   temp_liqu,
                                                   temp_vap,
                                                   coo_liqu,
                                                   coo_vap);
  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    {
      flux_normal_interp(fa7) *= surface_facettes(fa7);
      interfacial_temperature(fa7) *= surface_facettes(fa7);
    }
}

void IJK_Thermal_base::force_upstream_temperature(IJK_Field_double& temperature,
                                                  double T_imposed,
                                                  const IJK_Interfaces& interfaces,
                                                  double nb_diam,
                                                  int upstream_dir,
                                                  int gravity_dir,
                                                  int upstream_stencil)
{
  int dir = 0;
  if (upstream_dir == -1)
    {
      dir = gravity_dir;
      if (dir == -1)
        dir=0;
    }
  const Domaine_IJK& geom = temperature.get_domaine();

  bool perio =  geom.get_periodic_flag(dir);

  assert(interfaces.get_nb_bulles_reelles() == 1);
  DoubleTab bounding_box;
  // interfaces.calculer_bounding_box_bulles(bounding_box);
  bounding_box = interfaces.get_ijk_compo_connex().get_bounding_box();
  // Calcule la hauteur en x de la permiere bulle et la position de son cdg :
  const double Dbdir = bounding_box(0, dir, 1) - bounding_box(0, dir, 0);
  const double dirb  = ( bounding_box(0, dir, 1) + bounding_box(0, dir, 0) ) / 2.;
  const double ldir = geom.get_domain_length(dir) ;
  if (nb_diam == 0.)
    nb_diam = (ldir/Dbdir) / 2;
  double dirobj = dirb + nb_diam*Dbdir;

  // L'origine est sur un noeud. Donc que la premiere face en I est sur get_origin(DIRECTION_I)
  const double ddir = geom.get_constant_delta(dir);
  const double origin_dir = geom.get_origin(dir) ;
  const int offset_dir = geom.get_offset_local(dir);

  // FIXME: If nb_diam is too large it will iterate a lot
  if (perio)
    {
      while (dirobj<origin_dir)
        dirobj += ldir;
      while (dirobj>origin_dir+ldir)
        dirobj -= ldir;
    }

  // On devrait avoir xobj dans le domaine, sinon, on a choisi nb_diam trop grand :
  assert( ((dirobj>=origin_dir) && (dirobj <= origin_dir+ldir) ));

  const double x2 = (dirobj-origin_dir)/ ddir;
  int index_dir = (int)(floor(x2)) - offset_dir; // C'est l'index local, donc potentiellement negatif...
  const int& ndir = select_dir(dir, temperature.ni(), temperature.nj(), temperature.nk());

  // Cerr << "index_dir " << index_dir << finl;
  if ((index_dir >=0) && (index_dir < ndir))
    {
      // On est sur le bon proc...
      if (index_dir+upstream_stencil >= ndir)
        // On ne veut pas s'embeter sur 2 procs...
        index_dir = ndir-upstream_stencil;
    }
  else
    return;

  const double imposed = T_imposed;
  const int& imin = select_dir(dir, index_dir, 0, 0);
  const int& jmin = select_dir(dir, 0, index_dir, 0);
  const int& kmin = select_dir(dir, 0, 0, index_dir);
  const int& imax = select_dir(dir, imin + upstream_stencil, temperature.nj(), temperature.nk());
  const int& jmax = select_dir(dir, temperature.ni(), jmin + upstream_stencil, temperature.nk());
  const int& kmax = select_dir(dir, temperature.ni(), temperature.nj(), kmin + upstream_stencil);
  for (int k = kmin; k < kmax; k++)
    for (int j = jmin; j < jmax; j++)
      for (int i = imin; i < imax; i++)
        temperature(i,j,k) = imposed;
}


void IJK_Thermal_base::copy_previous_interface_state()
{
  thermal_local_subproblems_interfaces_fields_.copy_previous_interface_state();
}

int IJK_Thermal_base::post_process_quantities_from_subresolution(const Motcles& liste_post_instantanes,
                                                                 const char *lata_name,
                                                                 const int latastep)
{
  return thermal_local_subproblems_interfaces_fields_.posttraiter_champs_instantanes(liste_post_instantanes,
                                                                                     lata_name,
                                                                                     latastep);
}


void IJK_Thermal_base::posttraiter_tous_champs_thermal(Motcles& liste, const int idx) const
{
  /*
   * thermal_words_[0] = "subresolution";
   * thermal_words_[1] = "multiplesubresolutions";
   * thermal_words_[2] = "onefluid";
   * thermal_words_[3] = "onefluidenergy";
   * thermal_words_[4] = "cut_cell";
   */
  liste.add("TEMPERATURE");
  liste.add("TEMPERATURE_ANA");
  liste.add("ECART_T_ANA");
  liste.add("DIV_LAMBDA_GRAD_T_VOLUME");
  liste.add("U_T_CONVECTIVE_VOLUME");
  liste.add("GRAD_T");
  //
  liste.add("DISTANCE");
  liste.add("CURVATURE");
  if (thermal_rank_ == SUBRES || thermal_rank_ == MSUBRES)
    {
      /*
       * TODO: ADD SOME PARTICULAR FIELDS OR DO SWITCH CASE(thermal_rank)
       */
    }
  else
    {
      /*
       * TODO: CHECK IF GRAD_T0 MUST STILL BE POST-PROCESSED
       */
      liste.add("CP");
      liste.add("LAMBDA");
      //
      liste.add("SOURCE_TEMPERATURE");
      liste.add("TEMPERATURE_ADIM_BULLES");
      liste.add("TEMPERATURE_PHYSIQUE_T");
      liste.add("TEMPERATURE_ADIMENSIONNELLE_THETA");
      liste.add("SOURCE_TEMPERATURE_ANA");
      liste.add("ECART_SOURCE_TEMPERATURE_ANA");
      //
      liste.add("GRAD_T0");
      liste.add("GRAD_T1");
      liste.add("GRAD_T2");
      //
      liste.add("T_RUST");
      liste.add("DIV_RHO_CP_T_V");
    }
}

void IJK_Thermal_base::ecrire_statistiques_bulles(int reset, const Nom& nom_cas, const double current_time, const ArrOfDouble& surface, const int idx)
{
  ArrOfDouble interfacial_temperature;
  ArrOfDouble interfacial_phin_ai;
  // To transfer the field to FT splitting (because interfaces are there...) !!! NEEDED for compute_interfacial_temperature
  IJK_Field_double& temperature_ft = get_temperature_ft();
  ref_ijk_ft_->eq_ns().redistribute_to_splitting_ft_elem(*get_temperature(), temperature_ft);
  temperature_ft.echange_espace_virtuel(temperature_ft.ghost());
  //compute_interfacial_temperature(interfacial_temperature, interfacial_phin_ai, get_storage());
  compute_interfacial_temperature2(interfacial_temperature, interfacial_phin_ai);

  // Compute Bubble mean :
  ArrOfDouble Ti_per_bubble;
  ArrOfDouble phin_per_bubble;
  ref_ijk_ft_->get_interface().compute_surface_average_per_bubble(surface, interfacial_temperature, Ti_per_bubble);
  ref_ijk_ft_->get_interface().compute_surface_average_per_bubble(surface, interfacial_phin_ai, phin_per_bubble);
  if (Process::je_suis_maitre())
    {
      char s[1000];
      const char *nomcas = nom_cas;
      SFichier fic;
      const int n = Ti_per_bubble.size_array();
      IOS_OPEN_MODE mode = (reset) ? ios::out : ios::app;

#if INT_is_64_ == 1
      snprintf(s, 1000, "%s_bulles_Ti_%ld.out", nomcas, idx);
#else
      snprintf(s, 1000, "%s_bulles_Ti_%d.out", nomcas, idx);
#endif
      // Cerr << "Ecriture des donnees par bulles: fichier " << s << finl;
      fic.ouvrir(s, mode);
      snprintf(s, 1000, "%.16e ", current_time);
      fic << s;
      for (int i = 0; i < n; i++)
        {
          snprintf(s, 1000, "%.16e ", Ti_per_bubble[i]);
          fic << s;
        }
      fic << finl;
      fic.close();

      // Cerr << "Ecriture des donnees par bulles: fichier " << s << finl;
#if INT_is_64_ == 1
      snprintf(s, 1000, "%s_bulles_phin_%ld.out", nomcas, idx);
#else
      snprintf(s, 1000, "%s_bulles_phin_%d.out", nomcas, idx);
#endif
      fic.ouvrir(s, mode);
      snprintf(s, 1000, "%.16e ", current_time);
      fic << s;
      for (int i = 0; i < n; i++)
        {
          snprintf(s, 1000, "%.16e ", phin_per_bubble[i]);
          fic << s;
        }
      fic << finl;
      fic.close();

      Cerr << "Fin de l'ecriture des stats par bulles pour la temperature " << idx << finl;
    }
}



int IJK_Thermal_base::posttraiter_champs_instantanes_thermal_interface(const Motcles& liste_post_instantanes, const char *lata_name, const int latastep, const double current_time, const int idx)
{
  int n = 0; // nombre de champs postraites
  if (thermal_rank_ == SUBRES || thermal_rank_ == MSUBRES)
    {
      n = post_process_quantities_from_subresolution(liste_post_instantanes,
                                                     lata_name,
                                                     latastep);
      /*
       * TODO: COMPUTE INTERFACIAL GRADIENT
       */
    }
  else
    {
      n = posttraiter_champs_instantanes_thermal_interface_ref(liste_post_instantanes, lata_name, latastep, current_time, idx);
    }
  return n;
}

int IJK_Thermal_base::posttraiter_champs_instantanes_thermal_interface_ref(const Motcles& liste_post_instantanes, const char *lata_name, const int latastep, const double current_time, const int idx)
{
  Cerr << liste_post_instantanes << finl;
  int n = 0; // nombre de champs postraites
  Motcle lata_suffix = lata_suffix_[thermal_rank_];

  std::ostringstream oss;
  oss << "INTERFACE_TEMPERATURE_" << lata_suffix << idx;
  Nom nom_temp(oss.str().c_str());

  std::ostringstream oss2;
  oss2 << "INTERFACE_PHIN_" << lata_suffix << idx;
  Nom nom_phin(oss2.str().c_str());

  if ((liste_post_instantanes.contient_("INTERFACE_TEMPERATURE")) || (liste_post_instantanes.contient_("INTERFACE_PHIN")) || (liste_post_instantanes.contient_(nom_temp))
      || (liste_post_instantanes.contient_(nom_phin)))
    {
      //  Computing interfacial temperature at fa7 centre :
      const Maillage_FT_IJK& mesh = ref_ijk_ft_->get_maillage_ft_ijk(); // ref_ijk_ft_post_->interfaces_.maillage_ft_ijk();
      const ArrOfDouble& surface_facettes = mesh.get_update_surface_facettes();
      const int nb_facettes = mesh.nb_facettes();
      ArrOfDouble interfacial_temperature;
      ArrOfDouble interfacial_phin;
      // To transfer the field to FT splitting (because interfaces are there...) !!! NEEDED for compute_interfacial_temperature
      IJK_Field_double& temperature_ft = get_temperature_ft();
      ref_ijk_ft_->eq_ns().redistribute_to_splitting_ft_elem(*get_temperature(), temperature_ft);
      temperature_ft.echange_espace_virtuel(temperature_ft.ghost());
      // results are prop to the area :
      //itr.compute_interfacial_temperature(interfacial_temperature, interfacial_phin, itr.get_storage());
      compute_interfacial_temperature2(interfacial_temperature, interfacial_phin);
      for (int fa7 = 0; fa7 < nb_facettes; fa7++)
        {
          const double sf = surface_facettes[fa7];
          interfacial_temperature[fa7] /= sf;
          interfacial_phin[fa7] /= sf;
        }
      if ((liste_post_instantanes.contient_("INTERFACE_TEMPERATURE")) || (liste_post_instantanes.contient_(nom_temp)))
        n++, dumplata_ft_field(lata_name, "INTERFACES", nom_temp, "ELEM", interfacial_temperature, latastep);
      if ((liste_post_instantanes.contient_("INTERFACE_PHIN")) || (liste_post_instantanes.contient_(nom_phin)))
        n++, dumplata_ft_field(lata_name, "INTERFACES", nom_phin, "ELEM", interfacial_phin, latastep);
    }
  oss.str("");
  return n;
}

void IJK_Thermal_base::thermal_subresolution_outputs(const Nom& interfacial_quantities_thermal_probes, const Nom& shell_quantities_thermal_probes, const Nom& overall_bubbles_quantities, const Nom& local_quantities_thermal_probes_time_index_folder, const Nom& local_quantities_thermal_slices_time_index_folder, const Nom& local_quantities_thermal_lines_time_index_folder)
{
  if (thermal_rank_ == SUBRES || thermal_rank_ == MSUBRES)
    {
      set_thermal_subresolution_outputs(interfacial_quantities_thermal_probes, shell_quantities_thermal_probes, overall_bubbles_quantities, local_quantities_thermal_probes_time_index_folder);
      post_process_thermal_wake_slices(local_quantities_thermal_slices_time_index_folder);
      post_process_thermal_downstream_lines(local_quantities_thermal_lines_time_index_folder);
    }
}


int IJK_Thermal_base::posttraiter_champs_instantanes_thermal(const Motcles& liste_post_instantanes, const char *lata_name, const int latastep, const double current_time, const int idx)
{
  const int& rank = get_rank();
  Cerr << liste_post_instantanes << finl;
  int n = 0; // number of post-processed field
  std::ostringstream oss;
  Motcle lata_suffix = lata_suffix_[thermal_rank_];

  int cut_cell_activated = (thermal_rank_ == CUTCELL);

  /*
   * TEMPERATURE
   */
//  {
//    Motcles tested_names(1);
//    tested_names[0] = "TEMPERATURE";
//    post_process_std_thermal_field(liste_post_instantanes, lata_name, latastep, current_time, idx,
//                                   tested_names, "TEMPERATURE", lata_suffix, *get_temperature(), oss, n);
//  }
  /*
    for (auto &itr : liste_post_instantanes) {
  	  oss << nom << lata_suffix << idx;
  	  const IJK_Field_double& my_field = probleme().get_IJK_field(nom_complet);
        n++, dumplata_scalar(lata_name, nom_complet, my_field, latastep);
    }
  */
  oss << "TEMPERATURE_" << lata_suffix << idx;
  Nom nom_temp(oss.str().c_str());
  if ((liste_post_instantanes.contient_("TEMPERATURE")) || (liste_post_instantanes.contient_(nom_temp)))
    {
      n++, dumplata_scalar_cut_cell(cut_cell_activated, lata_name, nom_temp, get_temperature(), latastep);
    }

  oss.str("");

  oss << "SOURCE_TEMPERATURE_" << lata_suffix << idx;
  Nom nom_src_temp(oss.str().c_str());
  if ((liste_post_instantanes.contient_("SOURCE_TEMPERATURE")) || (liste_post_instantanes.contient_(nom_src_temp)))
    {
      n++, dumplata_scalar(lata_name, nom_src_temp, source_temperature_, latastep);
    }

  oss.str("");
  oss << "TEMPERATURE_ADIMENSIONNELLE_THETA_" << lata_suffix << idx;
  Nom nom_tempa(oss.str().c_str());
  if ((liste_post_instantanes.contient_("TEMPERATURE_ADIMENSIONNELLE_THETA")) || (liste_post_instantanes.contient_(nom_tempa)))
    {
      n++, dumplata_scalar(lata_name, nom_tempa, temperature_adimensionnelle_theta_, latastep);
    }
  oss.str("");

  oss << "TEMPERATURE_SMOOTH_" << lata_suffix << idx;
  Nom nom_temp_smooth(oss.str().c_str());
  if ((liste_post_instantanes.contient_("TEMPERATURE_SMOOTH")) || (liste_post_instantanes.contient_(nom_temp_smooth)))
    {
      n++, dumplata_scalar(lata_name, nom_temp_smooth, get_temperature_elem_smooth(), latastep);
    }
  oss.str("");

  /*
   * TEMPERATURE_BEFORE_EXTRAP
   */
  oss << "TEMPERATURE_BEFORE_EXTRAP_" << lata_suffix << idx;
  Nom nom_temp_before_extrap(oss.str().c_str());
  if ((liste_post_instantanes.contient_("TEMPERATURE_BEFORE_EXTRAP")) || (liste_post_instantanes.contient_(nom_temp_before_extrap)))
    {
      n++, dumplata_scalar(lata_name, nom_temp_before_extrap, get_temperature_before_extrapolation(), latastep);
    }
  oss.str("");

  /*
   * TEMPERATURE_CELL_NEIGHBOURS
   */
  oss << "TEMPERATURE_CELL_NEIGHBOURS_" << lata_suffix << idx;
  Nom nom_temp_cell_neighbours(oss.str().c_str());
  if ((liste_post_instantanes.contient_("TEMPERATURE_CELL_NEIGHBOURS")) || (liste_post_instantanes.contient_(nom_temp_cell_neighbours)))
    {
      n++, dumplata_scalar(lata_name, nom_temp_cell_neighbours, get_temperature_cell_neighbours(), latastep);
    }
  oss.str("");

  /*
   * TEMPERATURE_CELL_NEIGHBOURS_DEBUG
   */
  oss << "TEMPERATURE_CELL_NEIGHBOURS_DEBUG_" << lata_suffix << idx;
  Nom nom_temp_cell_neighbours_debug(oss.str().c_str());
  if ((liste_post_instantanes.contient_("TEMPERATURE_CELL_NEIGHBOURS_DEBUG")) || (liste_post_instantanes.contient_(nom_temp_cell_neighbours_debug)))
    {
      n++, dumplata_scalar(lata_name, nom_temp_cell_neighbours_debug, get_temperature_cell_neighbours_debug(), latastep);
    }
  oss.str("");

  /*
   * CELL_NEIGHBOURS_CORRECTED
   */
  oss << "CELL_NEIGHBOURS_CORRECTED_" << lata_suffix << idx;
  Nom nom_temp_cell_neighbours_corrected(oss.str().c_str());
  if ((liste_post_instantanes.contient_("CELL_NEIGHBOURS_CORRECTED")) || (liste_post_instantanes.contient_(nom_temp_cell_neighbours_corrected)))
    {
      n++, dumplata_scalar(lata_name, nom_temp_cell_neighbours_corrected, get_cell_neighbours_corrected(), latastep);
    }
  oss.str("");

  /*
   * CELL_NEIGHBOURS_COLINEARITY
   */
  oss << "CELL_NEIGHBOURS_COLINEARITY_" << lata_suffix << idx;
  Nom nom_temp_cell_neighbours_colinearity(oss.str().c_str());
  if ((liste_post_instantanes.contient_("CELL_NEIGHBOURS_COLINEARITY")) || (liste_post_instantanes.contient_(nom_temp_cell_neighbours_colinearity)))
    {
      n++, dumplata_scalar(lata_name, nom_temp_cell_neighbours_colinearity, get_neighbours_temperature_colinearity_weighting(), latastep);
    }
  oss.str("");

  /*
   * TEMPERATURE_ANA
   */
  oss << "TEMPERATURE_ANA_" << lata_suffix << idx;
  Nom nom_ana(oss.str().c_str());
  if ((liste_post_instantanes.contient_("TEMPERATURE_ANA")) || (liste_post_instantanes.contient_(nom_ana)) || (liste_post_instantanes.contient_("ECART_T_ANA")))
    {
      //set_field_data(itr.temperature_ana_, itr.expression_T_ana_, current_time);
      set_field_T_ana();
      n++, dumplata_scalar(lata_name, nom_ana, get_temperature_ana(), latastep);
    }
  oss.str("");

  /*
   * ECART_T_ANA
   */
  oss << "ECART_T_ANA_" << lata_suffix << idx;
  Nom nom_ecart_ana(oss.str().c_str());
  if ((liste_post_instantanes.contient_("ECART_T_ANA") || liste_post_instantanes.contient_(nom_ecart_ana)))
    {
      calculer_ecart_T_ana();
      n++, dumplata_scalar(lata_name, nom_ecart_ana, get_ecart_t_ana(), latastep);
    }
  oss.str("");

  /*
   * ECART_T_ANA_REL
   */
  oss << "ECART_T_ANA_REL_" << lata_suffix << idx;
  Nom nom_ecart_ana_rel(oss.str().c_str());
  if ((liste_post_instantanes.contient_("ECART_T_ANA_REL") || liste_post_instantanes.contient_(nom_ecart_ana_rel)))
    {
      // calculer_ecart_T_ana();
      n++, dumplata_scalar(lata_name, nom_ecart_ana_rel, get_ecart_t_ana_rel(), latastep);
    }
  oss.str("");

  /*
   * DIV_LAMBDA_GRAD_T_VOLUME
   */
  oss << "DIV_LAMBDA_GRAD_T_VOLUME_" << lata_suffix << idx;
  Nom nom_div_lambda_grad_T_volume(oss.str().c_str());
  if ((liste_post_instantanes.contient_("DIV_LAMBDA_GRAD_T_VOLUME") || liste_post_instantanes.contient_(nom_div_lambda_grad_T_volume)))
    {
      n++, dumplata_scalar_cut_cell(cut_cell_activated, lata_name, nom_div_lambda_grad_T_volume, get_div_lambda_grad_T_volume(), latastep);
    }
  oss.str("");

  /*
   * DIV_LAMBDA_GRAD_T
   */
  oss << "DIV_LAMBDA_GRAD_T_" << lata_suffix << idx;
  Nom nom_div_lambda_grad_T(oss.str().c_str());
  if ((liste_post_instantanes.contient_("DIV_LAMBDA_GRAD_T") || liste_post_instantanes.contient_(nom_div_lambda_grad_T)))
    {
      n++, dumplata_scalar(lata_name, nom_div_lambda_grad_T, get_div_lambda_grad_T(), latastep);
    }
  oss.str("");

  /*
   * U_T_CONVECTIVE_VOLUME
   */
  oss << "U_T_CONVECTIVE_VOLUME_" << lata_suffix << idx;
  Nom nom_u_T_convective_volume(oss.str().c_str());
  if ((liste_post_instantanes.contient_("U_T_CONVECTIVE_VOLUME") || liste_post_instantanes.contient_(nom_u_T_convective_volume)))
    {
      n++, dumplata_scalar(lata_name, nom_u_T_convective_volume, get_u_T_convective_volume(), latastep);
    }
  oss.str("");

  /*
   * U_T_CONVECTIVE
   */
  oss << "U_T_CONVECTIVE_" << lata_suffix << idx;
  Nom nom_u_T_convective(oss.str().c_str());
  if ((liste_post_instantanes.contient_("U_T_CONVECTIVE") || liste_post_instantanes.contient_(nom_u_T_convective)))
    {
      n++, dumplata_scalar(lata_name, nom_u_T_convective, get_u_T_convective(), latastep);
    }
  oss.str("");

  /*
   * GRAD_T
   */
  oss << "GRAD_T_" << lata_suffix << idx;
  Nom nom_grad(oss.str().c_str());
  if (liste_post_instantanes.contient_("GRAD_T") || (liste_post_instantanes.contient_(nom_grad)))
    {
      const IJK_Field_vector3_double& grad_T = get_grad_T();
      n++, dumplata_vector(lata_name, nom_grad, grad_T[0], grad_T[1], grad_T[2], latastep);
    }
  oss.str("");

  /*
   * GRAD_T_INTERFACE
   */
  oss << "GRAD_T_INTERFACE_" << lata_suffix << idx;
  Nom nom_grad_T_interface(oss.str().c_str());
  if ((liste_post_instantanes.contient_("GRAD_T_INTERFACE") || liste_post_instantanes.contient_(nom_grad_T_interface)))
    {
      n++, dumplata_scalar(lata_name, nom_grad_T_interface, get_grad_T_interface_ns(), latastep);
    }
  oss.str("");

  /*
   * GRAD_T_INTERFACE_FT
   */
  oss << "GRAD_T_INTERFACE_FT_" << lata_suffix << idx;
  Nom nom_grad_T_interface_ft(oss.str().c_str());
  if ((liste_post_instantanes.contient_("GRAD_T_INTERFACE_FT") || liste_post_instantanes.contient_(nom_grad_T_interface_ft)))
    {
      n++, dumplata_scalar(lata_name, nom_grad_T_interface_ft, get_grad_T_interface_ft(), latastep);
    }
  oss.str("");

  /*
   * TEMPERATURE_FT
   */
  oss << "TEMPERATURE_FT_" << lata_suffix << idx;
  Nom nom_temperature_ft(oss.str().c_str());
  if ((liste_post_instantanes.contient_("TEMPERATURE_FT") || liste_post_instantanes.contient_(nom_temperature_ft)))
    {
      n++, dumplata_scalar(lata_name, nom_temperature_ft, get_temperature_ft(), latastep);
    }
  oss.str("");

  /*
   * GRAD_T_DIR_X_ELEM
   */
  oss << "GRAD_T_DIR_X_ELEM_" << lata_suffix << idx;
  Nom nom_grad_T_dir_x(oss.str().c_str());
  if (liste_post_instantanes.contient_("GRAD_T_DIR_X_ELEM") || liste_post_instantanes.contient_("GRAD_T_ELEM") || liste_post_instantanes.contient_(nom_grad_T_dir_x))
    {
      n++, dumplata_scalar(lata_name, nom_grad_T_dir_x, get_gradient_temperature_elem()[0], latastep);
    }
  oss.str("");

  /*
   * GRAD_T_DIR_Y_ELEM
   */
  oss << "GRAD_T_DIR_Y_ELEM_" << lata_suffix << idx;
  Nom nom_grad_T_dir_y(oss.str().c_str());
  if (liste_post_instantanes.contient_("GRAD_T_DIR_Y_ELEM_") || liste_post_instantanes.contient_("GRAD_T_ELEM") || liste_post_instantanes.contient_(nom_grad_T_dir_y))
    {
      n++, dumplata_scalar(lata_name, nom_grad_T_dir_y, get_gradient_temperature_elem()[1], latastep);
    }
  oss.str("");

  /*
   * GRAD_T_DIR_Z_ELEM
   */
  oss << "GRAD_T_DIR_Z_ELEM_" << lata_suffix << idx;
  Nom nom_grad_T_dir_z(oss.str().c_str());
  if (liste_post_instantanes.contient_("GRAD_T_DIR_Z_ELEM_") || liste_post_instantanes.contient_("GRAD_T_ELEM") || liste_post_instantanes.contient_(nom_grad_T_dir_z))
    {
      n++, dumplata_scalar(lata_name, nom_grad_T_dir_z, get_gradient_temperature_elem()[2], latastep);
    }
  oss.str("");

  /*
   * GRAD_T_DIR_X_ELEM_SMOOTH
   */
  oss << "GRAD_T_DIR_X_ELEM_SMOOTH_" << lata_suffix << idx;
  Nom nom_grad_T_dir_x_smooth(oss.str().c_str());
  if (liste_post_instantanes.contient_("GRAD_T_DIR_X_ELEM_SMOOTH") || liste_post_instantanes.contient_("GRAD_T_ELEM_SMOOTH")
      || liste_post_instantanes.contient_(nom_grad_T_dir_x_smooth))
    {
      n++, dumplata_scalar(lata_name, nom_grad_T_dir_x_smooth, get_gradient_temperature_elem_smooth()[0], latastep);
    }
  oss.str("");

  /*
   * GRAD_T_DIR_Y_ELEM_SMOOTH
   */
  oss << "GRAD_T_DIR_Y_ELEM_SMOOTH_" << lata_suffix << idx;
  Nom nom_grad_T_dir_y_smooth(oss.str().c_str());
  if (liste_post_instantanes.contient_("GRAD_T_DIR_Y_ELEM_SMOOTH") || liste_post_instantanes.contient_("GRAD_T_ELEM_SMOOTH")
      || liste_post_instantanes.contient_(nom_grad_T_dir_y_smooth))
    {
      n++, dumplata_scalar(lata_name, nom_grad_T_dir_y_smooth, get_gradient_temperature_elem_smooth()[1], latastep);
    }
  oss.str("");

  /*
   * GRAD_T_DIR_Z_ELEM_SMOOTH
   */
  oss << "GRAD_T_DIR_Z_ELEM_SMOOTH_" << lata_suffix << idx;
  Nom nom_grad_T_dir_z_smooth(oss.str().c_str());
  if (liste_post_instantanes.contient_("GRAD_T_DIR_Z_ELEM_SMOOTH") || liste_post_instantanes.contient_("GRAD_T_ELEM_SMOOTH")
      || liste_post_instantanes.contient_(nom_grad_T_dir_z_smooth))
    {
      n++, dumplata_scalar(lata_name, nom_grad_T_dir_z_smooth, get_gradient_temperature_elem_smooth()[2], latastep);
    }
  oss.str("");

  /*
     * GRAD_T_DIR_X_ELEM_SMOOTH
     */
  oss << "GRAD_T_DIR_X_ELEM_TAN_SMOOTH_" << lata_suffix << idx;
  Nom nom_tangential_grad_T_dir_x_smooth(oss.str().c_str());
  if (liste_post_instantanes.contient_("GRAD_T_DIR_X_ELEM_TAN_SMOOTH") || liste_post_instantanes.contient_("GRAD_T_ELEM_TAN_SMOOTH")
      || liste_post_instantanes.contient_(nom_tangential_grad_T_dir_x_smooth))
    {
      n++, dumplata_scalar(lata_name, nom_tangential_grad_T_dir_x_smooth, get_tangential_gradient_temperature_elem_smooth()[0], latastep);
    }
  oss.str("");

  /*
   * GRAD_T_DIR_Y_ELEM_SMOOTH
   */
  oss << "GRAD_T_DIR_Y_ELEM_TAN_SMOOTH_" << lata_suffix << idx;
  Nom nom_tangential_grad_T_dir_y_smooth(oss.str().c_str());
  if (liste_post_instantanes.contient_("GRAD_T_DIR_Y_ELEM_TAN_SMOOTH") || liste_post_instantanes.contient_("GRAD_T_ELEM_TAN_SMOOTH")
      || liste_post_instantanes.contient_(nom_tangential_grad_T_dir_y_smooth))
    {
      n++, dumplata_scalar(lata_name, nom_tangential_grad_T_dir_y_smooth, get_tangential_gradient_temperature_elem_smooth()[1], latastep);
    }
  oss.str("");

  /*
   * GRAD_T_DIR_Z_ELEM_SMOOTH
   */
  oss << "GRAD_T_DIR_Z_ELEM_TAN_SMOOTH_" << lata_suffix << idx;
  Nom nom_tangential_grad_T_dir_z_smooth(oss.str().c_str());
  if (liste_post_instantanes.contient_("GRAD_T_DIR_Z_ELEM_TAN_SMOOTH") || liste_post_instantanes.contient_("GRAD_T_ELEM_TAN_SMOOTH")
      || liste_post_instantanes.contient_(nom_tangential_grad_T_dir_z_smooth))
    {
      n++, dumplata_scalar(lata_name, nom_tangential_grad_T_dir_z_smooth, get_tangential_gradient_temperature_elem_smooth()[2], latastep);
    }
  oss.str("");



  /*
   * HESS_T_DIR_XX_ELEM
   */
  oss << "HESS_T_DIR_XX_ELEM_" << lata_suffix << idx;
  Nom nom_hess_T_dir_xx(oss.str().c_str());
  if (liste_post_instantanes.contient_("HESS_T_DIR_XX_ELEM") || liste_post_instantanes.contient_("HESS_DIAG_T_ELEM") || liste_post_instantanes.contient_("HESS_T_ELEM")
      || liste_post_instantanes.contient_(nom_hess_T_dir_xx))
    {
      n++, dumplata_scalar(lata_name, nom_hess_T_dir_xx, get_hessian_diag_temperature_elem()[0], latastep);
    }
  oss.str("");

  /*
   * HESS_T_DIR_YY_ELEM
   */
  oss << "HESS_T_DIR_YY_ELEM_" << lata_suffix << idx;
  Nom nom_hess_T_dir_yy(oss.str().c_str());
  if (liste_post_instantanes.contient_("HESS_T_DIR_YY_ELEM") || liste_post_instantanes.contient_("HESS_DIAG_T_ELEM") || liste_post_instantanes.contient_("HESS_T_ELEM")
      || liste_post_instantanes.contient_(nom_hess_T_dir_yy))
    {
      n++, dumplata_scalar(lata_name, nom_hess_T_dir_yy, get_hessian_diag_temperature_elem()[1], latastep);
    }
  oss.str("");

  /*
   * HESS_T_DIR_ZZ_ELEM
   */
  oss << "HESS_T_DIR_ZZ_ELEM_" << lata_suffix << idx;
  Nom nom_hess_T_dir_zz(oss.str().c_str());
  if (liste_post_instantanes.contient_("HESS_T_DIR_ZZ_ELEM") || liste_post_instantanes.contient_("HESS_DIAG_T_ELEM") || liste_post_instantanes.contient_("HESS_T_ELEM")
      || liste_post_instantanes.contient_(nom_hess_T_dir_zz))
    {
      n++, dumplata_scalar(lata_name, nom_hess_T_dir_zz, get_hessian_diag_temperature_elem()[2], latastep);
    }
  oss.str("");

  /*
   * HESS_T_DIR_XY_YX_ELEM
   */
  oss << "HESS_T_DIR_XY_YX_ELEM_" << lata_suffix << idx;
  Nom nom_hess_T_dir_xy_yx(oss.str().c_str());
  if (liste_post_instantanes.contient_("HESS_T_DIR_XY_YX_ELEM") || liste_post_instantanes.contient_("HESS_CROSS_T_ELEM") || liste_post_instantanes.contient_("HESS_T_ELEM")
      || liste_post_instantanes.contient_(nom_hess_T_dir_xy_yx))
    {
      n++, dumplata_scalar(lata_name, nom_hess_T_dir_xy_yx, get_hessian_cross_temperature_elem()[2], latastep);
    }
  oss.str("");

  /*
   * HESS_T_DIR_XZ_ZX_ELEM
   */
  oss << "HESS_T_DIR_XZ_ZX_ELEM_" << lata_suffix << idx;
  Nom nom_hess_T_dir_xz_zx(oss.str().c_str());
  if (liste_post_instantanes.contient_("HESS_T_DIR_XZ_ZX_ELEM") || liste_post_instantanes.contient_("HESS_CROSS_T_ELEM") || liste_post_instantanes.contient_("HESS_T_ELEM")
      || liste_post_instantanes.contient_(nom_hess_T_dir_xz_zx))
    {
      n++, dumplata_scalar(lata_name, nom_hess_T_dir_xz_zx, get_hessian_cross_temperature_elem()[1], latastep);
    }
  oss.str("");

  /*
   * HESS_T_DIR_YZ_ELEM
   */
  oss << "HESS_T_DIR_YZ_ZY_ELEM_" << lata_suffix << idx;
  Nom nom_hess_T_dir_yz_zy(oss.str().c_str());
  if (liste_post_instantanes.contient_("HESS_T_DIR_YZ_ZY_ELEM") || liste_post_instantanes.contient_("HESS_CROSS_T_ELEM") || liste_post_instantanes.contient_("HESS_T_ELEM")
      || liste_post_instantanes.contient_(nom_hess_T_dir_yz_zy))
    {
      n++, dumplata_scalar(lata_name, nom_hess_T_dir_yz_zy, get_hessian_cross_temperature_elem()[0], latastep);
    }
  oss.str("");

  /*
   * RHO_CP_U_T_CONVECTIVE_FLUXES_X
   */
  oss << "RHO_CP_U_T_CONVECTIVE_FLUXES_X_" << lata_suffix << idx;
  Nom nom_rho_cp_u_T_convective_x(oss.str().c_str());
  if (liste_post_instantanes.contient_("RHO_CP_U_T_CONVECTIVE_FLUXES_X") || liste_post_instantanes.contient_(nom_rho_cp_u_T_convective_x))
    {
      n++, dumplata_scalar(lata_name, nom_rho_cp_u_T_convective_x, get_rho_cp_u_T_convective_fluxes()[0], latastep);
    }
  oss.str("");

  /*
   * RHO_CP_U_T_CONVECTIVE_FLUXES_Y
   */
  oss << "RHO_CP_U_T_CONVECTIVE_FLUXES_Y_" << lata_suffix << idx;
  Nom nom_rho_cp_u_T_convective_y(oss.str().c_str());
  if (liste_post_instantanes.contient_("RHO_CP_U_T_CONVECTIVE_FLUXES_Y") || liste_post_instantanes.contient_(nom_rho_cp_u_T_convective_y))
    {
      n++, dumplata_scalar(lata_name, nom_rho_cp_u_T_convective_y, get_rho_cp_u_T_convective_fluxes()[1], latastep);
    }
  oss.str("");

  /*
   * RHO_CP_U_T_CONVECTIVE_FLUXES_Z
   */
  oss << "RHO_CP_U_T_CONVECTIVE_FLUXES_Z_" << lata_suffix << idx;
  Nom nom_rho_cp_u_T_convective_z(oss.str().c_str());
  if (liste_post_instantanes.contient_("RHO_CP_U_T_CONVECTIVE_FLUXES_Z") || liste_post_instantanes.contient_(nom_rho_cp_u_T_convective_z))
    {
      n++, dumplata_scalar(lata_name, nom_rho_cp_u_T_convective_z, get_rho_cp_u_T_convective_fluxes()[2], latastep);
    }
  oss.str("");

  /*
   * LAMBDA_GRAD_T_DIFFUSIVE_FLUXES_X
   */
  oss << "LAMBDA_GRAD_T_DIFFUSIVE_FLUXES_X_" << lata_suffix << idx;
  Nom nom_coeff_grad_T_diffusive_x(oss.str().c_str());
  if (liste_post_instantanes.contient_("LAMBDA_GRAD_T_DIFFUSIVE_FLUXES_X") || liste_post_instantanes.contient_(nom_coeff_grad_T_diffusive_x))
    {
      n++, dumplata_scalar(lata_name, nom_coeff_grad_T_diffusive_x, get_div_coeff_grad_T_diffusive_fluxes()[0], latastep);
    }
  oss.str("");

  /*
   * LAMBDA_GRAD_T_DIFFUSIVE_FLUXES_Y
   */
  oss << "LAMBDA_GRAD_T_DIFFUSIVE_FLUXES_Y_" << lata_suffix << idx;
  Nom nom_coeff_grad_T_diffusive_y(oss.str().c_str());
  if (liste_post_instantanes.contient_("LAMBDA_GRAD_T_DIFFUSIVE_FLUXES_Y") || liste_post_instantanes.contient_(nom_coeff_grad_T_diffusive_y))
    {
      n++, dumplata_scalar(lata_name, nom_coeff_grad_T_diffusive_y, get_div_coeff_grad_T_diffusive_fluxes()[1], latastep);
    }
  oss.str("");

  /*
   * LAMBDA_GRAD_T_DIFFUSIVE_FLUXES_Z
   */
  oss << "LAMBDA_GRAD_T_DIFFUSIVE_FLUXES_Z_" << lata_suffix << idx;
  Nom nom_coeff_grad_T_diffusive_z(oss.str().c_str());
  if (liste_post_instantanes.contient_("LAMBDA_GRAD_T_DIFFUSIVE_FLUXES_Z") || liste_post_instantanes.contient_(nom_coeff_grad_T_diffusive_z))
    {
      n++, dumplata_scalar(lata_name, nom_coeff_grad_T_diffusive_z, get_div_coeff_grad_T_diffusive_fluxes()[2], latastep);
    }
  oss.str("");

  if (rank == 0)
    {
      /*
       * DISTANCE
       */
      // oss << "DISTANCE_" << lata_suffix << idx;
      oss << "DISTANCE";
      Nom nom_distance(oss.str().c_str());
      if ((liste_post_instantanes.contient_("DISTANCE") || liste_post_instantanes.contient_(nom_distance)))
        {
          n++, dumplata_scalar(lata_name, nom_distance, get_eulerian_distance_ns(), latastep);
        }
      oss.str("");

      /*
       * NORMAL_VECTOR_X
       */
      // oss << "NORMAL_VECTOR_X_" << lata_suffix << idx;
      oss << "NORMAL_VECTOR_X";
      Nom nom_normal_vector_x(oss.str().c_str());
      if (liste_post_instantanes.contient_("NORMAL_VECTOR_X") || liste_post_instantanes.contient_(nom_normal_vector_x))
        {
          n++, dumplata_scalar(lata_name, nom_normal_vector_x, get_normal_vector_ns()[0], latastep);
        }
      oss.str("");

      /*
       * NORMAL_VECTOR_Y
       */
      // oss << "NORMAL_VECTOR_Y_" << lata_suffix << idx;
      oss << "NORMAL_VECTOR_Y";
      Nom nom_normal_vector_y(oss.str().c_str());
      if (liste_post_instantanes.contient_("NORMAL_VECTOR_Y") || liste_post_instantanes.contient_(nom_normal_vector_y))
        {
          n++, dumplata_scalar(lata_name, nom_normal_vector_y, get_normal_vector_ns()[1], latastep);
        }
      oss.str("");

      /*
       * NORMAL_VECTOR_Z
       */
      // oss << "NORMAL_VECTOR_Z_" << lata_suffix << idx;
      oss << "NORMAL_VECTOR_Z";
      Nom nom_normal_vector_z(oss.str().c_str());
      if (liste_post_instantanes.contient_("NORMAL_VECTOR_Z") || liste_post_instantanes.contient_(nom_normal_vector_z))
        {
          n++, dumplata_scalar(lata_name, nom_normal_vector_z, get_normal_vector_ns()[2], latastep);
        }
      oss.str("");

      /*
       * NORMAL_VECTOR_X_NORMED
       */
      // oss << "NORMAL_VECTOR_X_NORMED_" << lata_suffix << idx;
      oss << "NORMAL_VECTOR_X_NORMED";
      Nom nom_normal_vector_x_normed(oss.str().c_str());
      if (liste_post_instantanes.contient_("NORMAL_VECTOR_X_NORMED") || liste_post_instantanes.contient_(nom_normal_vector_x_normed))
        {
          n++, dumplata_scalar(lata_name, nom_normal_vector_x_normed, get_normal_vector_ns_normed()[0], latastep);
        }
      oss.str("");

      /*
       * NORMAL_VECTOR_Y_NORMED
       */
      // oss << "NORMAL_VECTOR_Y_NORMED_" << lata_suffix << idx;
      oss << "NORMAL_VECTOR_Y_NORMED";
      Nom nom_normal_vector_y_normed(oss.str().c_str());
      if (liste_post_instantanes.contient_("NORMAL_VECTOR_Y_NORMED") || liste_post_instantanes.contient_(nom_normal_vector_y_normed))
        {
          n++, dumplata_scalar(lata_name, nom_normal_vector_y_normed, get_normal_vector_ns_normed()[1], latastep);
        }
      oss.str("");

      /*
       * NORMAL_VECTOR_Z_NORMED
       */
      // oss << "NORMAL_VECTOR_Z_NORMED_" << lata_suffix << idx;
      oss << "NORMAL_VECTOR_Z_NORMED";
      Nom nom_normal_vector_z_normed(oss.str().c_str());
      if (liste_post_instantanes.contient_("NORMAL_VECTOR_Z_NORMED") || liste_post_instantanes.contient_(nom_normal_vector_z_normed))
        {
          n++, dumplata_scalar(lata_name, nom_normal_vector_z_normed, get_normal_vector_ns_normed()[2], latastep);
        }
      oss.str("");

      /*
       * CURVATURE
       */
      // oss << "CURVATURE_" << lata_suffix << idx;
      oss << "CURVATURE";
      Nom nom_curvature(oss.str().c_str());
      if ((liste_post_instantanes.contient_("CURVATURE") || liste_post_instantanes.contient_(nom_curvature)))
        {
          n++, dumplata_scalar(lata_name, nom_curvature, get_eulerian_curvature_ns(), latastep);
        }
      oss.str("");

      /*
       * INTERFACIAL AREA
       */
      // oss << "INTERFACIAL_AREA_" << lata_suffix << idx;
      oss << "INTERFACIAL_AREA";
      Nom nom_interfacial_area(oss.str().c_str());
      if ((liste_post_instantanes.contient_("INTERFACIAL_AREA") || liste_post_instantanes.contient_(nom_interfacial_area)))
        {
          n++, dumplata_scalar(lata_name, nom_interfacial_area, get_interfacial_area_ns(), latastep);
        }
      oss.str("");

      /*
       * DISTANCE_FT
       */
      // oss << "DISTANCE_FT_" << lata_suffix << idx;
      oss << "DISTANCE_FT";
      Nom nom_distance_ft(oss.str().c_str());
      if ((liste_post_instantanes.contient_("DISTANCE_FT") || liste_post_instantanes.contient_(nom_distance_ft)))
        {
          n++, dumplata_scalar(lata_name, nom_distance_ft, get_eulerian_distance_ft(), latastep);
        }
      oss.str("");

      /*
       * NORMAL_VECTOR_X_FT
       */
      // oss << "NORMAL_VECTOR_X_FT_" << lata_suffix << idx;
      oss << "NORMAL_VECTOR_X_FT";
      Nom nom_normal_vector_x_ft(oss.str().c_str());
      if (liste_post_instantanes.contient_("NORMAL_VECTOR_X_FT") || liste_post_instantanes.contient_("NORMAL_VECTOR_X_FT") || liste_post_instantanes.contient_(nom_normal_vector_x_ft))
        {
          n++, dumplata_scalar(lata_name, nom_normal_vector_x_ft, get_normal_vector_ft()[0], latastep);
        }
      oss.str("");

      /*
       * NORMAL_VECTOR_Y_FT
       */
      // oss << "NORMAL_VECTOR_Y_FT_" << lata_suffix << idx;
      oss << "NORMAL_VECTOR_Y_FT";
      Nom nom_normal_vector_y_ft(oss.str().c_str());
      if (liste_post_instantanes.contient_("NORMAL_VECTOR_Y_FT") || liste_post_instantanes.contient_("NORMAL_VECTOR_Y_FT") || liste_post_instantanes.contient_(nom_normal_vector_y_ft))
        {
          n++, dumplata_scalar(lata_name, nom_normal_vector_y_ft, get_normal_vector_ft()[1], latastep);
        }
      oss.str("");

      /*
       * NORMAL_VECTOR_Z_FT
       */
      // oss << "NORMAL_VECTOR_Z_FT_" << lata_suffix << idx;
      oss << "NORMAL_VECTOR_Z_FT";
      Nom nom_normal_vector_z_ft(oss.str().c_str());
      if (liste_post_instantanes.contient_("NORMAL_VECTOR_Z_FT") || liste_post_instantanes.contient_("NORMAL_VECTOR_Z_FT") || liste_post_instantanes.contient_(nom_normal_vector_z_ft))
        {
          n++, dumplata_scalar(lata_name, nom_normal_vector_z_ft, get_normal_vector_ft()[2], latastep);
        }
      oss.str("");

      /*
       * CURVATURE_FT
       */
      // oss << "CURVATURE_FT_" << lata_suffix << idx;
      oss << "CURVATURE_FT";
      Nom nom_curvature_ft(oss.str().c_str());
      if ((liste_post_instantanes.contient_("CURVATURE_FT") || liste_post_instantanes.contient_(nom_curvature_ft)))
        {
          n++, dumplata_scalar(lata_name, nom_curvature_ft, get_eulerian_curvature_ft(), latastep);
        }
      oss.str("");

      /*
       * INTERFACIAL_AREA_FT
       */
      // oss << "INTERFACIAL_AREA_FT_" << lata_suffix << idx;
      oss << "INTERFACIAL_AREA_FT";
      Nom nom_interfacial_area_ft(oss.str().c_str());
      if ((liste_post_instantanes.contient_("INTERFACIAL_AREA_FT") || liste_post_instantanes.contient_(nom_interfacial_area_ft)))
        {
          n++, dumplata_scalar(lata_name, nom_interfacial_area_ft, get_interfacial_area_ft(), latastep);
        }
      oss.str("");

      /*
       * INDICATOR_FT
       */
      // oss << "INDICATOR_FT_" << lata_suffix << idx;
      oss << "INDICATOR_FT";
      Nom nom_indicator_ft(oss.str().c_str());
      if ((liste_post_instantanes.contient_("INDICATOR_FT") || liste_post_instantanes.contient_(nom_indicator_ft)))
        {
          n++, dumplata_scalar(lata_name, nom_indicator_ft, ref_ijk_ft_->get_interface().I_ft(), latastep);
        }
      oss.str("");

      /*
       * BARY_X
       */
      // oss << "BARY_X_" << lata_suffix << idx;
      oss << "BARY_X";
      Nom nom_bary_x(oss.str().c_str());
      if (liste_post_instantanes.contient_("BARY") || liste_post_instantanes.contient_("BARY_X") || liste_post_instantanes.contient_(nom_bary_x))
        {
          n++, dumplata_scalar(lata_name, nom_bary_x, get_bary()[0], latastep);
        }
      oss.str("");

      /*
       * BARY_Y
       */
      // oss << "BARY_Y_" << lata_suffix << idx;
      oss << "BARY_Y";
      Nom nom_bary_y(oss.str().c_str());
      if (liste_post_instantanes.contient_("BARY") || liste_post_instantanes.contient_("BARY_Y") || liste_post_instantanes.contient_(nom_bary_y))
        {
          n++, dumplata_scalar(lata_name, nom_bary_y, get_bary()[1], latastep);
        }
      oss.str("");

      /*
       * BARY_Z
       */
      // oss << "BARY_Z_" << lata_suffix << idx;
      oss << "BARY_Z";
      Nom nom_bary_z(oss.str().c_str());
      if (liste_post_instantanes.contient_("BARY") || liste_post_instantanes.contient_("BARY_Z") || liste_post_instantanes.contient_(nom_bary_z))

        {
          n++, dumplata_scalar(lata_name, nom_bary_z, get_bary()[2], latastep);
        }
      oss.str("");

      /*
       * EULERIAN_COMPO_FT
       */
      // oss << "EULERIAN_COMPO_FT_" << lata_suffix << idx;
      oss << "EULERIAN_COMPO_FT";
      Nom nom_eulerian_compo_ft(oss.str().c_str());
      if (liste_post_instantanes.contient_("EULERIAN_COMPO_FT") || liste_post_instantanes.contient_(nom_eulerian_compo_ft))

        {
          n++, dumplata_scalar(lata_name, nom_eulerian_compo_ft, get_eulerian_compo_connex_ft(), latastep);
        }
      oss.str("");

      /*
       * EULERIAN_COMPO
       */
      // oss << "EULERIAN_COMPO_" << lata_suffix << idx;
      oss << "EULERIAN_COMPO";
      Nom nom_eulerian_compo(oss.str().c_str());
      if (liste_post_instantanes.contient_("EULERIAN_COMPO") || liste_post_instantanes.contient_(nom_eulerian_compo))

        {
          n++, dumplata_scalar(lata_name, nom_eulerian_compo, get_eulerian_compo_connex_ns(), latastep);
        }
      oss.str("");

      /*
       * EULERIAN_COMPO_GHOST_FT
       */
      // oss << "EULERIAN_COMPO_GHOST_FT_" << lata_suffix << idx;
      oss << "EULERIAN_COMPO_GHOST_FT";
      Nom nom_eulerian_compo_ghost_ft(oss.str().c_str());
      if (liste_post_instantanes.contient_("EULERIAN_COMPO_GHOST_FT") || liste_post_instantanes.contient_(nom_eulerian_compo_ghost_ft))

        {
          n++, dumplata_scalar(lata_name, nom_eulerian_compo_ghost_ft, get_eulerian_compo_connex_ghost_ft(), latastep);
        }
      oss.str("");

      /*
       * EULERIAN_COMPO_GHOST
       */
      // oss << "EULERIAN_COMPO_GHOST_" << lata_suffix << idx;
      oss << "EULERIAN_COMPO_GHOST";
      Nom nom_eulerian_compo_ghost(oss.str().c_str());
      if (liste_post_instantanes.contient_("EULERIAN_COMPO_GHOST") || liste_post_instantanes.contient_(nom_eulerian_compo_ghost))

        {
          n++, dumplata_scalar(lata_name, nom_eulerian_compo_ghost, get_eulerian_compo_connex_ghost_ns(), latastep);
        }
      oss.str("");

      /*
       * EULERIAN_COMPO_INTERFACE_FT
       */
      // oss << "EULERIAN_COMPO_INTERFACE_FT_" << lata_suffix << idx;
      oss << "EULERIAN_COMPO_INTERFACE_FT";
      Nom nom_eulerian_compo_from_interface_ft(oss.str().c_str());
      if (liste_post_instantanes.contient_("EULERIAN_COMPO_INTERFACE_FT") || liste_post_instantanes.contient_(nom_eulerian_compo_from_interface_ft))

        {
          n++, dumplata_scalar(lata_name, nom_eulerian_compo_from_interface_ft, get_eulerian_compo_connex_from_interface_ft(), latastep);
        }
      oss.str("");

      /*
       * EULERIAN_COMPO_INTERFACE_GHOST_FT
       */
      // oss << "EULERIAN_COMPO_INTERFACE_GHOST_FT_" << lata_suffix << idx;
      oss << "EULERIAN_COMPO_INTERFACE_GHOST_FT";
      Nom nom_eulerian_compo_from_interface_ghost_ft(oss.str().c_str());
      if (liste_post_instantanes.contient_("EULERIAN_COMPO_INTERFACE_GHOST_FT") || liste_post_instantanes.contient_(nom_eulerian_compo_from_interface_ghost_ft))

        {
          n++, dumplata_scalar(lata_name, nom_eulerian_compo_from_interface_ghost_ft, get_eulerian_compo_connex_from_interface_ghost_ft(), latastep);
        }
      oss.str("");

      /*
       * EULERIAN_COMPO_INTERFACE_NS
       */
      // oss << "EULERIAN_COMPO_INTERFACE_NS_" << lata_suffix << idx;
      oss << "EULERIAN_COMPO_INTERFACE_NS";
      Nom nom_eulerian_compo_from_interface_ns(oss.str().c_str());
      if (liste_post_instantanes.contient_("EULERIAN_COMPO_INTERFACE_NS") || liste_post_instantanes.contient_(nom_eulerian_compo_from_interface_ns))

        {
          n++, dumplata_scalar(lata_name, nom_eulerian_compo_from_interface_ns, get_eulerian_compo_connex_from_interface_ns(), latastep);
        }
      oss.str("");

      /*
       * EULERIAN_COMPO_INTERFACE_GHOST_NS
       */
      // oss << "EULERIAN_COMPO_INTERFACE_GHOST_NS_" << lata_suffix << idx;
      oss << "EULERIAN_COMPO_INTERFACE_GHOST_NS";
      Nom nom_eulerian_compo_from_interface_ghost_ns(oss.str().c_str());
      if (liste_post_instantanes.contient_("EULERIAN_COMPO_INTERFACE_GHOST_NS") || liste_post_instantanes.contient_(nom_eulerian_compo_from_interface_ghost_ns))

        {
          n++, dumplata_scalar(lata_name, nom_eulerian_compo_from_interface_ghost_ns, get_eulerian_compo_connex_from_interface_ghost_ns(), latastep);
        }
      oss.str("");

      /*
       * EULERIAN_COMPO_INTERFACE_INT
       */
      // oss << "EULERIAN_COMPO_INTERFACE_INT_" << lata_suffix << idx;
      oss << "EULERIAN_COMPO_INTERFACE_INT";
      Nom nom_eulerian_compo_from_interface_int(oss.str().c_str());
      if (liste_post_instantanes.contient_("EULERIAN_COMPO_INTERFACE_INT") || liste_post_instantanes.contient_(nom_eulerian_compo_from_interface_int))

        {
          n++, dumplata_scalar(lata_name, nom_eulerian_compo_from_interface_int, get_eulerian_compo_connex_int_from_interface_ns(), latastep);
        }
      oss.str("");

      /*
       * EULERIAN_COMPO_INTERFACE_GHOST_INT
       */
      // oss << "EULERIAN_COMPO_INTERFACE_GHOST_INT_" << lata_suffix << idx;
      oss << "EULERIAN_COMPO_INTERFACE_GHOST_INT";
      Nom nom_eulerian_compo_from_interface_ghost_int(oss.str().c_str());
      if (liste_post_instantanes.contient_("EULERIAN_COMPO_INTERFACE_GHOST_INT") || liste_post_instantanes.contient_(nom_eulerian_compo_from_interface_ghost_int))

        {
          n++, dumplata_scalar(lata_name, nom_eulerian_compo_from_interface_ghost_int, get_eulerian_compo_connex_int_from_interface_ghost_ns(), latastep);
        }
      oss.str("");

      /*
       * RISING_VELOCITIES
       */
      // oss << "RISING_VELOCITIES_" << lata_suffix << idx;
      oss << "RISING_VELOCITIES";
      Nom nom_eulerian_rising_velocities(oss.str().c_str());
      if (liste_post_instantanes.contient_("RISING_VELOCITIES") || liste_post_instantanes.contient_(nom_eulerian_rising_velocities))

        {
          n++, dumplata_scalar(lata_name, nom_eulerian_rising_velocities, get_eulerian_rising_velocities(), latastep);
        }
      oss.str("");
    }


  /*
   * DEBUG_LRS_CELLS
   */
  if (get_debug())
    {
      oss << "DEBUG_LRS_CELLS_" << lata_suffix << idx;
      Nom nom_debug_lrs_cells(oss.str().c_str());
      if (liste_post_instantanes.contient_("DEBUG_LRS_CELLS") || liste_post_instantanes.contient_(nom_debug_lrs_cells))
        {
          n++, dumplata_scalar(lata_name, nom_debug_lrs_cells, get_debug_lrs_cells(), latastep);
        }
      oss.str("");
    }

  /*
   * CELL_FACES_CORRECTED_CONVECTIVE_X
   */
  oss << "CELL_FACES_CORRECTED_CONVECTIVE_X_" << lata_suffix << idx;
  Nom nom_cell_faces_corrected_convective_x(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_CORRECTED_CONVECTIVE_X") || liste_post_instantanes.contient_(nom_cell_faces_corrected_convective_x))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_corrected_convective_x, get_cell_faces_corrected_convective()[0], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_CORRECTED_CONVECTIVE_Y
   */
  oss << "CELL_FACES_CORRECTED_CONVECTIVE_Y_" << lata_suffix << idx;
  Nom nom_cell_faces_corrected_convective_y(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_CORRECTED_CONVECTIVE_Y") || liste_post_instantanes.contient_(nom_cell_faces_corrected_convective_y))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_corrected_convective_y, get_cell_faces_corrected_convective()[1], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_CORRECTED_CONVECTIVE_Z
   */
  oss << "CELL_FACES_CORRECTED_CONVECTIVE_Z_" << lata_suffix << idx;
  Nom nom_cell_faces_corrected_convective_z(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_CORRECTED_CONVECTIVE_Z") || liste_post_instantanes.contient_(nom_cell_faces_corrected_convective_z))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_corrected_convective_z, get_cell_faces_corrected_convective()[2], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_CORRECTED_DIFFUSIVE_X
   */
  oss << "CELL_FACES_CORRECTED_DIFFUSIVE_X_" << lata_suffix << idx;
  Nom nom_cell_faces_corrected_diffusive_x(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_CORRECTED_DIFFUSIVE_X") || liste_post_instantanes.contient_(nom_cell_faces_corrected_diffusive_x))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_corrected_diffusive_x, get_cell_faces_corrected_diffusive()[0], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_CORRECTED_DIFFUSIVE_Y
   */
  oss << "CELL_FACES_CORRECTED_DIFFUSIVE_Y_" << lata_suffix << idx;
  Nom nom_cell_faces_corrected_diffusive_y(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_CORRECTED_DIFFUSIVE_Y") || liste_post_instantanes.contient_(nom_cell_faces_corrected_diffusive_y))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_corrected_diffusive_y, get_cell_faces_corrected_diffusive()[1], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_CORRECTED_DIFFUSIVE_Z
   */
  oss << "CELL_FACES_CORRECTED_DIFFUSIVE_Z_" << lata_suffix << idx;
  Nom nom_cell_faces_corrected_diffusive_z(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_CORRECTED_DIFFUSIVE_Z") || liste_post_instantanes.contient_(nom_cell_faces_corrected_diffusive_z))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_corrected_diffusive_z, get_cell_faces_corrected_diffusive()[2], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_CORRECTED_BOOL_X
   */
  oss << "CELL_FACES_CORRECTED_BOOL_X_" << lata_suffix << idx;
  Nom nom_cell_faces_corrected_bool_x(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_CORRECTED_BOOL_X") || liste_post_instantanes.contient_(nom_cell_faces_corrected_bool_x))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_corrected_bool_x, get_cell_faces_corrected_bool()[0], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_CORRECTED_BOOL_Y
   */
  oss << "CELL_FACES_CORRECTED_BOOL_Y_" << lata_suffix << idx;
  Nom nom_cell_faces_corrected_bool_y(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_CORRECTED_BOOL_Y") || liste_post_instantanes.contient_(nom_cell_faces_corrected_bool_y))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_corrected_bool_y, get_cell_faces_corrected_bool()[1], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_CORRECTED_BOOL_Z
   */
  oss << "CELL_FACES_CORRECTED_BOOL_Z_" << lata_suffix << idx;
  Nom nom_cell_faces_corrected_bool_z(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_CORRECTED_BOOL_Z") || liste_post_instantanes.contient_(nom_cell_faces_corrected_bool_z))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_corrected_bool_z, get_cell_faces_corrected_bool()[2], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_DIAG_BOOL_X
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_DIAG_BOOL_X_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_diag_bool_x(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_DIAG_BOOL_X") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_diag_bool_x))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_diag_bool_x, get_cell_faces_neighbours_corrected_diag_bool()[0], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_DIAG_BOOL_Y
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_DIAG_BOOL_Y_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_diag_bool_y(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_DIAG_BOOL_Y") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_diag_bool_y))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_diag_bool_y, get_cell_faces_neighbours_corrected_diag_bool()[1], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_DIAG_BOOL_Z
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_DIAG_BOOL_Z_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_diag_bool_z(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_DIAG_BOOL_Z") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_diag_bool_z))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_diag_bool_z, get_cell_faces_neighbours_corrected_diag_bool()[2], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_ALL_BOOL_X
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_ALL_BOOL_X_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_all_bool_x(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_ALL_BOOL_X") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_all_bool_x))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_all_bool_x, get_cell_faces_neighbours_corrected_all_bool()[0], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_ALL_BOOL_Y
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_ALL_BOOL_Y_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_all_bool_y(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_ALL_BOOL_Y") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_all_bool_y))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_all_bool_y, get_cell_faces_neighbours_corrected_all_bool()[1], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_ALL_BOOL_Z
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_ALL_BOOL_Z_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_all_bool_z(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_ALL_BOOL_Z") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_all_bool_z))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_all_bool_z, get_cell_faces_neighbours_corrected_all_bool()[2], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_MIN_MAX_BOOL_X
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_MIN_MAX_BOOL_X_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_min_max_bool_x(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_MIN_MAX_BOOL_X") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_min_max_bool_x))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_min_max_bool_x, get_cell_faces_neighbours_corrected_min_max_bool()[0], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_MIN_MAX_BOOL_Y
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_MIN_MAX_BOOL_Y_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_min_max_bool_y(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_MIN_MAX_BOOL_Y") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_min_max_bool_y))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_min_max_bool_y, get_cell_faces_neighbours_corrected_min_max_bool()[1], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_MIN_MAX_BOOL_Z
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_MIN_MAX_BOOL_Z_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_min_max_bool_z(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_MIN_MAX_BOOL_Z") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_min_max_bool_z))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_min_max_bool_z, get_cell_faces_neighbours_corrected_min_max_bool()[2], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_U_T_X
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_U_T_X_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_u_T_x(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_U_T_X") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_u_T_x))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_u_T_x, get_cell_faces_neighbours_corrected_velocity_temperature()[0], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_U_T_Y
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_U_T_Y_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_u_T_y(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_U_T_Y") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_u_T_y))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_u_T_y, get_cell_faces_neighbours_corrected_velocity_temperature()[1], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_U_T_Z
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_U_T_Z_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_u_T_z(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_U_T_Z") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_u_T_z))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_u_T_z, get_cell_faces_neighbours_corrected_velocity_temperature()[2], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_CONVECTIVE_FRAME_REF_X
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_CONVECTIVE_FRAME_REF_X_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_convective_frame_ref_x(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_CONVECTIVE_FRAME_REF_X") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_convective_frame_ref_x))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_convective_frame_ref_x, get_cell_faces_neighbours_corrected_convective_frame_of_ref()[0], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_CONVECTIVE_FRAME_REF_Y
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_CONVECTIVE_FRAME_REF_Y_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_convective_frame_ref_y(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_CONVECTIVE_FRAME_REF_Y") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_convective_frame_ref_y))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_convective_frame_ref_y, get_cell_faces_neighbours_corrected_convective_frame_of_ref()[1], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_CONVECTIVE_FRAME_REF_Z
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_CONVECTIVE_FRAME_REF_Z_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_convective_frame_ref_z(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_CONVECTIVE_FRAME_REF_Z") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_convective_frame_ref_z))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_convective_frame_ref_z, get_cell_faces_neighbours_corrected_convective_frame_of_ref()[2], latastep);
    }
  oss.str("");


  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_CONVECTIVE_X
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_CONVECTIVE_X_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_convective_x(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_CONVECTIVE_X") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_convective_x))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_convective_x, get_cell_faces_neighbours_corrected_convective()[0], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_CONVECTIVE_Y
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_CONVECTIVE_Y_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_convective_y(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_CONVECTIVE_Y") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_convective_y))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_convective_y, get_cell_faces_neighbours_corrected_convective()[1], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_CONVECTIVE_Z
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_CONVECTIVE_Z_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_convective_z(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_CONVECTIVE_Z") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_convective_z))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_convective_z, get_cell_faces_neighbours_corrected_convective()[2], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_DIFFUSIVE_X
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_DIFFUSIVE_X_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_diffusive_x(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_DIFFUSIVE_X") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_diffusive_x))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_diffusive_x, get_cell_faces_neighbours_corrected_diffusive()[0], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_DIFFUSIVE_Y
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_DIFFUSIVE_Y_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_diffusive_y(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_DIFFUSIVE_Y") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_diffusive_y))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_diffusive_y, get_cell_faces_neighbours_corrected_diffusive()[1], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_DIFFUSIVE_Z
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_DIFFUSIVE_Z_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_diffusive_z(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_DIFFUSIVE_Z") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_diffusive_z))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_diffusive_z, get_cell_faces_neighbours_corrected_diffusive()[2], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_COLINEARITY_X
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_COLINEARITY_X_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_colinearity_x(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_COLINEARITY_X") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_colinearity_x))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_colinearity_x, get_neighbours_faces_weighting_colinearity()[0], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_COLINEARITY_Y
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_COLINEARITY_Y_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_colinearity_y(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_COLINEARITY_Y") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_colinearity_y))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_colinearity_y, get_neighbours_faces_weighting_colinearity()[1], latastep);
    }
  oss.str("");

  /*
   * CELL_FACES_NEIGHBOURS_CORRECTED_COLINEARITY_Z
   */
  oss << "CELL_FACES_NEIGHBOURS_CORRECTED_COLINEARITY_Z_" << lata_suffix << idx;
  Nom nom_cell_faces_neighbours_corrected_colinearity_z(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_FACES_NEIGHBOURS_CORRECTED_COLINEARITY_Z") || liste_post_instantanes.contient_(nom_cell_faces_neighbours_corrected_colinearity_z))
    {
      n++, dumplata_scalar(lata_name, nom_cell_faces_neighbours_corrected_colinearity_z, get_neighbours_faces_weighting_colinearity()[2], latastep);
    }
  oss.str("");

  /*
   * CELL_NEIGHBOURS_CORRECTED_TRIMMED
   */
  oss << "CELL_NEIGHBOURS_CORRECTED_TRIMMED_" << lata_suffix << idx;
  Nom nom_cell_neighbours_corrected_trimmed(oss.str().c_str());
  if (liste_post_instantanes.contient_("CELL_NEIGHBOURS_CORRECTED_TRIMMED") || liste_post_instantanes.contient_(nom_cell_neighbours_corrected_trimmed))
    {
      n++, dumplata_scalar(lata_name, nom_cell_neighbours_corrected_trimmed, get_cell_neighbours_corrected_trimmed(), latastep);
    }
  oss.str("");

  if (rank == 0)
    {
      /*
       * PROBE_COLLISION_DEBUG
       */
      oss << "PROBE_COLLISION_DEBUG";
      Nom nom_probe_collision_debug(oss.str().c_str());
      if (liste_post_instantanes.contient_("PROBE_COLLISION_DEBUG") || liste_post_instantanes.contient_(nom_probe_collision_debug))
        {
          n++, dumplata_scalar(lata_name, nom_probe_collision_debug, get_probe_collision_debug_field(), latastep);
        }
      oss.str("");
    }

  /*
   * INTERFACIAL_FLUX_DISPATCH_X
   */
  oss << "INTERFACIAL_FLUX_DISPATCH_X_" << lata_suffix << idx;
  Nom nom_interfacial_flux_dispatch_x(oss.str().c_str());
  if (liste_post_instantanes.contient_("INTERFACIAL_FLUX_DISPATCH_X") || liste_post_instantanes.contient_(nom_interfacial_flux_dispatch_x))
    {
      n++, dumplata_scalar(lata_name, nom_interfacial_flux_dispatch_x, get_interfacial_heat_flux_dispatched()[0], latastep);
    }
  oss.str("");

  /*
   * INTERFACIAL_FLUX_DISPATCH_Y
   */
  oss << "INTERFACIAL_FLUX_DISPATCH_Y_" << lata_suffix << idx;
  Nom nom_interfacial_flux_dispatch_y(oss.str().c_str());
  if (liste_post_instantanes.contient_("INTERFACIAL_FLUX_DISPATCH_Y") || liste_post_instantanes.contient_(nom_interfacial_flux_dispatch_y))
    {
      n++, dumplata_scalar(lata_name, nom_interfacial_flux_dispatch_y, get_interfacial_heat_flux_dispatched()[1], latastep);
    }
  oss.str("");

  /*
   * INTERFACIAL_FLUX_DISPATCH_Z
   */
  oss << "INTERFACIAL_FLUX_DISPATCH_Z_" << lata_suffix << idx;
  Nom nom_interfacial_flux_dispatch_z(oss.str().c_str());
  if (liste_post_instantanes.contient_("INTERFACIAL_FLUX_DISPATCH_Z") || liste_post_instantanes.contient_(nom_interfacial_flux_dispatch_z))
    {
      n++, dumplata_scalar(lata_name, nom_interfacial_flux_dispatch_z, get_interfacial_heat_flux_dispatched()[2], latastep);
    }
  oss.str("");


  /*
   * INTERFACIAL_FLUX_CONTRIB_X
   */
  oss << "INTERFACIAL_FLUX_CONTRIB_X_" << lata_suffix << idx;
  Nom nom_interfacial_flux_contrib_x(oss.str().c_str());
  if (liste_post_instantanes.contient_("INTERFACIAL_FLUX_CONTRIB_X") || liste_post_instantanes.contient_(nom_interfacial_flux_contrib_x))
    {
      n++, dumplata_scalar(lata_name, nom_interfacial_flux_contrib_x, get_interfacial_heat_flux_contrib()[0], latastep);
    }
  oss.str("");

  /*
   * INTERFACIAL_FLUX_CONTRIB_Y
   */
  oss << "INTERFACIAL_FLUX_CONTRIB_Y_" << lata_suffix << idx;
  Nom nom_interfacial_flux_contrib_y(oss.str().c_str());
  if (liste_post_instantanes.contient_("INTERFACIAL_FLUX_CONTRIB_Y") || liste_post_instantanes.contient_(nom_interfacial_flux_contrib_y))
    {
      n++, dumplata_scalar(lata_name, nom_interfacial_flux_contrib_y, get_interfacial_heat_flux_contrib()[1], latastep);
    }
  oss.str("");

  /*
   * INTERFACIAL_FLUX_CONTRIB_Z
   */
  oss << "INTERFACIAL_FLUX_CONTRIB_Z_" << lata_suffix << idx;
  Nom nom_interfacial_flux_contrib_z(oss.str().c_str());
  if (liste_post_instantanes.contient_("INTERFACIAL_FLUX_CONTRIB_Z") || liste_post_instantanes.contient_(nom_interfacial_flux_contrib_z))
    {
      n++, dumplata_scalar(lata_name, nom_interfacial_flux_contrib_z, get_interfacial_heat_flux_contrib()[2], latastep);
    }
  oss.str("");

  /*
   * INTERFACIAL_FLUX_CURRENT_X
   */
  oss << "INTERFACIAL_FLUX_CURRENT_X_" << lata_suffix << idx;
  Nom nom_interfacial_flux_current_x(oss.str().c_str());
  if (liste_post_instantanes.contient_("INTERFACIAL_FLUX_CURRENT_X") || liste_post_instantanes.contient_(nom_interfacial_flux_current_x))
    {
      n++, dumplata_scalar(lata_name, nom_interfacial_flux_current_x, get_interfacial_heat_flux_current()[0], latastep);
    }
  oss.str("");

  /*
   * INTERFACIAL_FLUX_CURRENT_Y
   */
  oss << "INTERFACIAL_FLUX_CURRENT_Y_" << lata_suffix << idx;
  Nom nom_interfacial_flux_current_y(oss.str().c_str());
  if (liste_post_instantanes.contient_("INTERFACIAL_FLUX_CURRENT_Y") || liste_post_instantanes.contient_(nom_interfacial_flux_current_y))
    {
      n++, dumplata_scalar(lata_name, nom_interfacial_flux_current_y, get_interfacial_heat_flux_current()[1], latastep);
    }
  oss.str("");

  /*
   * INTERFACIAL_FLUX_CURRENT_Z
   */
  oss << "INTERFACIAL_FLUX_CURRENT_Z_" << lata_suffix << idx;
  Nom nom_interfacial_flux_current_z(oss.str().c_str());
  if (liste_post_instantanes.contient_("INTERFACIAL_FLUX_CURRENT_Z") || liste_post_instantanes.contient_(nom_interfacial_flux_current_z))
    {
      n++, dumplata_scalar(lata_name, nom_interfacial_flux_current_z, get_interfacial_heat_flux_current()[2], latastep);
    }
  oss.str("");

  return n;
}
