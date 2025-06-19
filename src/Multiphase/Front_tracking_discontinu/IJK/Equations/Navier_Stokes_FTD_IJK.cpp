/****************************************************************************
* Copyright (c) 2024, CEA
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
#include <Probleme_FTD_IJK_tools.h>
#include <Navier_Stokes_FTD_IJK.h>
#include <Probleme_FTD_IJK_base.h>
#include <Fluide_Diphasique_IJK.h>
#include <Schema_RK3_IJK.h>
#include <Option_IJK.h>
#include <Param.h>

#define COMPLEMENT_ANTI_DEVIATION_RESIDU
// #define VARIABLE_DZ
// #define PROJECTION_DE_LINCREMENT_DV

Implemente_instanciable_sans_constructeur(Navier_Stokes_FTD_IJK, "Navier_Stokes_FTD_IJK", Equation_base);
// XD Navier_Stokes_FTD_IJK eqn_base Navier_Stokes_FTD_IJK -1 Navier-Stokes equations.

Navier_Stokes_FTD_IJK::Navier_Stokes_FTD_IJK()
{
  terme_source_correction_.resize_array(3); // Initialement a zero, puis sera calcule a chaque iter.
  terme_source_correction_ = 0.;
  correction_force_.resize_array(3); // Par defaut, les flags d'activations sont a zero (ie inactif).
  correction_force_ = 0;

  vol_bulles_.resize_array(0); // Initialement a zero, puis sera calcule a chaque iter.
  vol_bulles_ = 0.;

  expression_variable_source_.dimensionner_force(3);
  expression_vitesse_initiale_.dimensionner_force(3);
}

Sortie& Navier_Stokes_FTD_IJK::printOn(Sortie& os) const
{
  return os<< que_suis_je() << " " << le_nom();
}

Entree& Navier_Stokes_FTD_IJK::readOn(Entree& is)
{
  Cerr<<"Reading of data for a "<<que_suis_je()<<" equation"<<finl;
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Navier_Stokes_FTD_IJK::set_param(Param& param)
{
  /*
   * BASE
   */
  param.ajouter("multigrid_solver", &poisson_solver_, Param::REQUIRED); // XD_ADD_P multigrid_solver not_set
  param.ajouter("vitesse_entree_dir", &vitesse_entree_dir_);
  param.ajouter("vitesse_entree_compo_to_force", &vitesse_entree_compo_to_force_);
  param.ajouter("stencil_vitesse_entree", &stencil_vitesse_entree_);
  param.ajouter("vitesse_entree", &vitesse_entree_); // XD_ADD_P floattant Velocity to prescribe at inlet
  param.ajouter("expression_vx_init", &expression_vitesse_initiale_[0]); // XD_ADD_P chaine initial field for x-velocity component (parser of x,y,z)
  param.ajouter("expression_vy_init", &expression_vitesse_initiale_[1]); // XD_ADD_P chaine initial field for y-velocity component (parser of x,y,z)
  param.ajouter("expression_vz_init", &expression_vitesse_initiale_[2]); // XD_ADD_P chaine initial field for z-velocity component (parser of x,y,z)
  param.ajouter("expression_p_init", &expression_pression_initiale_); // XD_ADD_P chaine initial pressure field (optional)
  param.ajouter("velocity_diffusion_op", &velocity_diffusion_op_);
  param.ajouter("velocity_convection_op", &velocity_convection_op_); // XD_ADD_P chaine Type of velocity convection scheme

  // ATTENTION les fichiers reprises sont des fichiers .lata ou sauv.lata
  // On peut reprendre uniquement la vitesse ou uniquement rho dans un fichier de post:
  param.ajouter("fichier_reprise_vitesse", &fichier_reprise_vitesse_); // XD_ADD_P chaine not_set
  param.ajouter("timestep_reprise_vitesse", &timestep_reprise_vitesse_); // XD_ADD_P chaine not_set

  param.ajouter("boundary_conditions", &boundary_conditions_, Param::REQUIRED); // XD_ADD_P bloc_lecture BC

  param.ajouter_flag("projection_initiale", &projection_initiale_demandee_);
  param.ajouter_flag("disable_solveur_poisson", &disable_solveur_poisson_); // XD_ADD_P rien Disable pressure poisson solver
  param.ajouter_flag("disable_diffusion_qdm", &disable_diffusion_qdm_); // XD_ADD_P rien Disable diffusion operator in momentum
  param.ajouter_flag("disable_convection_qdm", &disable_convection_qdm_); // XD_ADD_P rien Disable convection operator in momentum

  // XXX Equation non resolu : renomer
  param.ajouter_flag("frozen_velocity", &frozen_velocity_); // XD_ADD_P chaine not_set

  param.ajouter_flag("velocity_reset", &velocity_reset_); // XD_ADD_P chaine not_set
  param.ajouter_flag("resolution_fluctuations", &resolution_fluctuations_); // XD_ADD_P rien Disable pressure poisson solver
  param.ajouter_flag("use_harmonic_viscosity", &use_harmonic_viscosity_);
  param.ajouter_flag("harmonic_nu_in_diff_operator", &harmonic_nu_in_diff_operator_); // XD_ADD_P rien Disable pressure poisson solver
  param.ajouter_flag("use_inv_rho_for_mass_solver_and_calculer_rho_v", &use_inv_rho_for_mass_solver_and_calculer_rho_v_); // XD_ADD_P chaine not_set
  param.ajouter_flag("use_inv_rho_in_poisson_solver", &use_inv_rho_in_poisson_solver_); // XD_ADD_P flag not_set
  param.ajouter_flag("diffusion_alternative", &diffusion_alternative_); // XD_ADD_P chaine not_set
  param.ajouter_flag("test_etapes_et_bilan", &test_etapes_et_bilan_); // XD_ADD_P chaine not_set
  param.ajouter_flag("ajout_init_a_reprise", &add_initial_field_); // XD_ADD_P chaine not_set
  param.ajouter_flag("improved_initial_pressure_guess", &improved_initial_pressure_guess_); // XD_ADD_P chaine not_set
  param.ajouter_flag("include_pressure_gradient_in_ustar", &include_pressure_gradient_in_ustar_); // XD_ADD_P chaine not_set


  /*
   * FT
   */
  param.ajouter("upstream_dir", &upstream_dir_); // XD_ADD_P entier Direction to prescribe the velocity
  param.ajouter("vitesse_upstream", &vitesse_upstream_); // XD_ADD_P floattant Velocity to prescribe at 'nb_diam_upstream_' before bubble 0.
  param.ajouter("expression_vitesse_upstream", &expression_vitesse_upstream_); // XD_ADD_P chaine Analytical expression to set the upstream velocity
  param.ajouter("upstream_velocity_bubble_factor", &upstream_velocity_bubble_factor_);
  param.ajouter("upstream_velocity_bubble_factor_deriv", &upstream_velocity_bubble_factor_deriv_);
  param.ajouter("upstream_velocity_bubble_factor_integral", &upstream_velocity_bubble_factor_integral_);
  param.ajouter("velocity_bubble_scope", &velocity_bubble_scope_);
  param.ajouter("upstream_stencil", &upstream_stencil_); // XD_ADD_P int Width on which the velocity is set
  param.ajouter("nb_diam_upstream", &nb_diam_upstream_); // XD_ADD_P floattant Number of bubble diameters upstream of bubble 0 to prescribe the velocity.
  param.ajouter("nb_diam_ortho_shear_perio", &nb_diam_ortho_shear_perio_); // XD_ADD_P chaine not_set
  param.ajouter("vol_bulle_monodisperse", &vol_bulle_monodisperse_); // XD_ADD_P chaine not_set
  param.ajouter("diam_bulle_monodisperse", &diam_bulle_monodisperse_); // XD_ADD_P chaine not_set
  param.ajouter("coeff_evol_volume", &coeff_evol_volume_); // XD_ADD_P chaine not_set
  param.ajouter("vol_bulles", &vol_bulles_); // XD_ADD_P chaine not_set
  param.ajouter("reprise_vap_velocity_tmoy", &vap_velocity_tmoy_); // XD_ADD_P chaine not_set
  param.ajouter("reprise_liq_velocity_tmoy", &liq_velocity_tmoy_); // XD_ADD_P chaine not_set
  param.ajouter_flag("disable_source_interf", &disable_source_interf_); // XD_ADD_P rien Disable computation of the interfacial source term
  param.ajouter_flag("upstream_velocity_measured", &upstream_velocity_measured_);
  param.ajouter_flag("harmonic_nu_in_calc_with_indicatrice", &harmonic_nu_in_calc_with_indicatrice_); // XD_ADD_P rien Disable pressure poisson solver
  param.ajouter_flag("refuse_patch_conservation_QdM_RK3_source_interf", &refuse_patch_conservation_QdM_RK3_source_interf_); // XD_ADD_P rien experimental Keyword, not for use
  param.ajouter_flag("suppression_rejetons", &suppression_rejetons_); // XD_ADD_P chaine not_set


  /*
   * Sources
   */
  param.ajouter("p_seuil_max", &p_seuil_max_); // XD_ADD_P floattant not_set, default 10000000
  param.ajouter("p_seuil_min", &p_seuil_min_); // XD_ADD_P floattant not_set, default -10000000
  param.ajouter("coef_ammortissement", &coef_ammortissement_); // XD_ADD_P floattant not_set
  param.ajouter("coef_immobilisation", &coef_immobilisation_); // XD_ADD_P floattant not_set
  param.ajouter("expression_derivee_force", &expression_derivee_acceleration_); // XD_ADD_P chaine expression of the time-derivative of the X-component of a source-term (see terme_force_ini for the initial value). terme_force_ini : initial value of the X-component of the source term (see expression_derivee_force  for time evolution)
  param.ajouter("terme_force_init", &terme_source_acceleration_); // XD_ADD_P chaine not_set
  param.ajouter("correction_force", &correction_force_); // XD_ADD_P chaine not_set
  param.ajouter_flag("compute_force_init", &compute_force_init_); // XD_ADD_P chaine not_set

  param.ajouter("expression_variable_source_x", &expression_variable_source_[0]); // XD_ADD_P chaine not_set
  param.ajouter("expression_variable_source_y", &expression_variable_source_[1]); // XD_ADD_P chaine not_set
  param.ajouter("expression_variable_source_z", &expression_variable_source_[2]); // XD_ADD_P chaine not_set
  param.ajouter("facteur_variable_source_init", &facteur_variable_source_); // XD_ADD_P chaine not_set
  param.ajouter("expression_derivee_facteur_variable_source", &expression_derivee_facteur_variable_source_); // XD_ADD_P chaine not_set

  param.ajouter("expression_potential_phi", &expression_potential_phi_); // XD_ADD_P chaine parser to define phi and make a momentum source Nabla phi.

  param.ajouter("forcage", &forcage_);  // XD_ADD_P chaine not_set
  param.ajouter("corrections_qdm", &qdm_corrections_); // XD_ADD_P chaine not_set


  // Correcteur PID
  param.ajouter("Kp", &Kp_);
  param.ajouter("Kd", &Kd_);
  param.ajouter("Ki", &Ki_);
  param.ajouter("epaisseur_maille", &epaisseur_maille_);



  /*
   * Sources FT
   */
  param.ajouter("coef_mean_force", &coef_mean_force_); // XD_ADD_P floattant not_set
  param.ajouter("coef_force_time_n", &coef_force_time_n_); // XD_ADD_P floattant not_set
  param.ajouter("coef_rayon_force_rappel", &coef_rayon_force_rappel_); // XD_ADD_P floattant not_set
  param.ajouter_flag("correction_semi_locale_volume_bulle", &correction_semi_locale_volume_bulle_);
}

void Navier_Stokes_FTD_IJK::set_param_reprise_pb(Param& param)
{
  /*
   * BASE
   */
  param.ajouter("terme_acceleration_init", &terme_source_acceleration_);
  param.ajouter("fichier_reprise_vitesse", &fichier_reprise_vitesse_);
  param.ajouter("timestep_reprise_vitesse", &timestep_reprise_vitesse_);
  param.ajouter("forcage", &forcage_);
  param.ajouter("corrections_qdm", &qdm_corrections_);

  /*
   * FT
   */
  param.ajouter("vitesse_upstream_reprise", &vitesse_upstream_reprise_);
  param.ajouter("velocity_bubble_old", &velocity_bubble_old_);
  param.ajouter("reprise_vap_velocity_tmoy", &vap_velocity_tmoy_);
  param.ajouter("reprise_liq_velocity_tmoy", &liq_velocity_tmoy_);
}

void Navier_Stokes_FTD_IJK::associer_pb_base(const Probleme_base& pb)
{
  if (!sub_type(Probleme_FTD_IJK_base, pb))
    {
      Cerr << "Error for the method Navier_Stokes_FTD_IJK::associer_pb_base\n";
      Cerr << " Navier_Stokes_FTD_IJK equation must be associated to\n";
      Cerr << " a Probleme_FTD_IJK_base problem type\n";
      Process::exit();
    }
  mon_probleme = pb;
  if (nom_ == "??")
    {
      nom_ = pb.le_nom();
      nom_ += que_suis_je();
    }
}

const Milieu_base& Navier_Stokes_FTD_IJK::milieu() const
{
  if (le_fluide_.est_nul())
    {
      Cerr << "You forgot to associate a fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return le_fluide_.valeur();
}

Milieu_base& Navier_Stokes_FTD_IJK::milieu()
{
  if (le_fluide_.est_nul())
    {
      Cerr << "You forgot to associate a fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return le_fluide_.valeur();
}

void Navier_Stokes_FTD_IJK::associer_milieu_base(const Milieu_base& un_milieu)
{
  if (sub_type(Fluide_Diphasique_IJK, un_milieu))
    {
      const Milieu_base& un_fluide = ref_cast(Milieu_base, un_milieu);
      le_fluide_ = un_fluide;
    }
  else
    {
      Cerr << "Error of fluid type for the method Navier_Stokes_FTD_IJK::associer_milieu_base" << finl;
      Process::exit();
    }
}

Probleme_FTD_IJK_base& Navier_Stokes_FTD_IJK::probleme_ijk()
{
  return ref_cast(Probleme_FTD_IJK_base, mon_probleme.valeur());
}

const Probleme_FTD_IJK_base& Navier_Stokes_FTD_IJK::probleme_ijk() const
{
  return ref_cast(Probleme_FTD_IJK_base, mon_probleme.valeur());
}

void Navier_Stokes_FTD_IJK::Fill_postprocessable_fields(std::vector<FieldInfo_t>& chps)
{
  std::vector<FieldInfo_t> c =
  {
    // Name     /     Localisation (elem, face, ...) /    Nature (scalare, vector)   /    Located on interface?

    { "VELOCITY", Entity::FACE, Nature_du_champ::vectoriel, false },
    { "VITESSE", Entity::FACE, Nature_du_champ::vectoriel, false },
    { "VELOCITY_FT", Entity::FACE, Nature_du_champ::vectoriel, false },
    { "PRESSURE", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "PRESSION", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "PRESSURE_RHS", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "PRESSION_RHS", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "SOURCE_QDM_INTERF", Entity::FACE, Nature_du_champ::vectoriel, false },
    { "SHIELD_REPULSTION", Entity::FACE, Nature_du_champ::vectoriel, false },
    { "SHIELD_REPULSTION_ABS", Entity::FACE, Nature_du_champ::vectoriel, false },
    { "RHO", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "DENSITY", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "MU", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "VISCOSITY", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "EXTERNAL_FORCE", Entity::FACE, Nature_du_champ::vectoriel, false },
  };
  chps.insert(chps.end(), c.begin(), c.end());
}

void Navier_Stokes_FTD_IJK::get_noms_champs_postraitables(Noms& noms,Option opt) const
{
  for (const auto& n : champs_compris_.liste_noms_compris())
    noms.add(n);
  for (const auto& n : champs_compris_.liste_noms_compris_vectoriel())
    noms.add(n);
}

const IJK_Field_vector3_double& Navier_Stokes_FTD_IJK::get_IJK_field_vector(const Motcle& nom)
{

  if (nom=="EXTERNAL_FORCE")
    {
      if ( coef_immobilisation_ > 1e-16)
        {
          if (!probleme_ijk().get_interface().get_forcing_method())
            for (int dir = 0; dir < 3; dir++)
              redistribute_from_splitting_ft_faces_[dir].redistribute( force_rappel_ft_[dir], force_rappel_[dir]);
        }
      else
        Cerr << "Posttraitement demande pour EXTERNAL_FORCE but ignored because coef_immobilisation_ <= 1e-16" << endl;
    }

  if (has_champ_vectoriel(nom))
    return champs_compris_.get_champ_vectoriel(nom);

  Cerr << "ERROR in Navier_Stokes_FTD_IJK::get_IJK_field_vector : " << finl;
  Cerr << "Requested field '" << nom << "' is not recognized by Navier_Stokes_FTD_IJK::get_IJK_field_vector()." << finl;
  throw;
}

const IJK_Field_double& Navier_Stokes_FTD_IJK::get_IJK_field(const Motcle& nom)
{
  if (has_champ(nom))
    return champs_compris_.get_champ(nom);

  Cerr << "ERROR in Navier_Stokes_FTD_IJK::get_IJK_field : " << finl;
  Cerr << "Requested field '" << nom << "' is not recognized by Navier_Stokes_FTD_IJK::get_IJK_field()." << finl;
  throw;
}

void Navier_Stokes_FTD_IJK::completer()
{
  if (frozen_velocity_)
    {
      disable_solveur_poisson_ = 1; // automatically force the suppression of the poisson solver
      if (!Option_IJK::DISABLE_DIPHASIQUE)
        {
          interfaces_->freeze();  // Stop the interfacial displacement.
          Cout << "The option frozen_velocity automatically freeze the interface motion " << "by activating the flag interfaces_.frozen_" << finl;
        }
    }

  if ((expression_potential_phi_ != "??") && ((expression_variable_source_[0] != "??") || (expression_variable_source_[1] != "??") || (expression_variable_source_[2] != "??")))
    {
      Cerr << "expression_potential_phi and expression_variable_source are used together" << "nabla(phi) will be added to the expression given for the variable source" << finl;
      //Process::exit();
    }

  if ((include_pressure_gradient_in_ustar_) && (expression_pression_initiale_ == "??"))
    {
      Cerr << "When using pressure increment in u^star, expression_p_init " << expression_pression_initiale_ << "becomes a required parameter. " << finl;
      Process::exit();
    }

  // Si on utilise un seul groupe et qu'on impose un volume unique a toutes les bulles,
  if (vol_bulle_monodisperse_ >= 0.)
    {
      if (vol_bulles_.size_array() != 0)
        {
          Cerr << "Attention, conflit entre les options : vol_bulle_monodisperse_ et vol_bulles." << "Merci de choisir" << finl;
          Process::exit();
        }
    }

  // On utilise inv_rho pour l'un ou l'autre... Il faut donc le calculer :
  use_inv_rho_ = use_inv_rho_for_mass_solver_and_calculer_rho_v_ + use_inv_rho_in_poisson_solver_;

  // Avec cette option, on travaille avec nu :
  if (diffusion_alternative_)
    {
      Cerr << "Option diffusion_alternative activee : le champ mu contient nu (la viscosite dynamique)" << finl;
      Cerr << "TODO FIXME " << finl;
      Process::exit();
    }

  if (vol_bulle_monodisperse_ != -1 || diam_bulle_monodisperse_ != -1)
    {
      if (vol_bulle_monodisperse_ != -1)
        diam_bulle_monodisperse_ = pow(6. * vol_bulle_monodisperse_ / (M_PI), 1. / 3.);
      else
        vol_bulle_monodisperse_ = M_PI * pow(diam_bulle_monodisperse_, 3) / 6.;

      probleme_ijk().get_domaine_ft().set_extension_from_bulle_param(vol_bulle_monodisperse_, diam_bulle_monodisperse_);
    }

  // Preparation de l'expression derivee de l'acceleration
  std::string tmpstring(expression_derivee_acceleration_);
  parser_derivee_acceleration_.setString(tmpstring);
  parser_derivee_acceleration_.setNbVar(9);
  parser_derivee_acceleration_.addVar("rappel_moyen");
  parser_derivee_acceleration_.addVar("force");
  parser_derivee_acceleration_.addVar("v_moyen");
  parser_derivee_acceleration_.addVar("ur");
  parser_derivee_acceleration_.addVar("ul");
  parser_derivee_acceleration_.addVar("uv");
  parser_derivee_acceleration_.addVar("T");
  parser_derivee_acceleration_.addVar("rhov_moyen");
  parser_derivee_acceleration_.addVar("tauw");
  parser_derivee_acceleration_.parseString();

  std::string tmpstring2(expression_derivee_facteur_variable_source_);
  parser_derivee_facteur_variable_source_.setString(tmpstring2);
  parser_derivee_facteur_variable_source_.setNbVar(6);
  parser_derivee_facteur_variable_source_.addVar("rappel_moyen");
  parser_derivee_facteur_variable_source_.addVar("facteur_sv");
  parser_derivee_facteur_variable_source_.addVar("v_moyen");
  parser_derivee_facteur_variable_source_.addVar("T");
  parser_derivee_facteur_variable_source_.addVar("rhov_moyen");
  parser_derivee_facteur_variable_source_.addVar("tauw");
  parser_derivee_facteur_variable_source_.parseString();

}

void Navier_Stokes_FTD_IJK::initialise_velocity_from_file(const Nom& fichier_reprise_vitesse)
{
  Cout << "Lecture vitesse initiale dans fichier " << fichier_reprise_vitesse << " timestep= " << timestep_reprise_vitesse_ << finl;
  const Nom& geom_name = velocity_[0].get_domaine().le_nom();
  lire_dans_lata(fichier_reprise_vitesse, timestep_reprise_vitesse_, geom_name, "VELOCITY", velocity_[0], velocity_[1], velocity_[2]); // fonction qui lit un champ a partir d'un lata .

  if (add_initial_field_)
    {
      compose_field_data(velocity_[0], expression_vitesse_initiale_[0]);
      compose_field_data(velocity_[1], expression_vitesse_initiale_[1]);
      compose_field_data(velocity_[2], expression_vitesse_initiale_[2]);
//          velocity_[0].data() += expression_vitesse_initiale_;
//          velocity_[1].data() += expression_vitesse_initiale_;
//          velocity_[2].data() += expression_vitesse_initiale_;
    }

  velocity_[0].echange_espace_virtuel(2);
  velocity_[1].echange_espace_virtuel(2);
  velocity_[2].echange_espace_virtuel(2);
#ifdef CONVERT_AT_READING_FROM_NURESAFE_TO_ADIM_TRYGGVASON_FOR_LIQUID_VELOCITY
  const double coef = 14.353432757182377;
  for (int dir=0; dir< 3; dir++)
    velocity_[dir].data() *= coef;
#endif
}

void Navier_Stokes_FTD_IJK::initialise_velocity_using_expression(const Noms& expression_vitesse_initiale)
{
  if (expression_vitesse_initiale_.size() != 3)
    {
      Cerr << "Erreur dans l'initialisation: la vitesse initiale doit etre fournie avec trois expressions" << finl;
      Process::exit();
    }
  else
    {
      Cout << "Initialisation vitesse \nvx = " << expression_vitesse_initiale[0] << "\nvy = " << expression_vitesse_initiale[1] << "\nvz = " << expression_vitesse_initiale[2] << finl;
      for (int i = 0; i < 3; i++)
        {
          // Cette methode parcours ni(), nj() et nk() et donc pas les ghost...
          set_field_data(velocity_[i], expression_vitesse_initiale[i]);
        }

      velocity_[0].echange_espace_virtuel(2);
      velocity_[1].echange_espace_virtuel(2);
      velocity_[2].echange_espace_virtuel(2);
    }
}

// Maj indicatrice rho mu met indicatrice a indicatrice next
// et maj rho et mu en fonction de la nouvelle indicatrice
void Navier_Stokes_FTD_IJK::maj_indicatrice_rho_mu(const bool parcourir)
{
  // En monophasique, les champs sont a jours donc on zap :
  if (Option_IJK::DISABLE_DIPHASIQUE)
    return;

  static Stat_Counter_Id calculer_rho_mu_indicatrice_counter_ = statistiques().new_counter(2, "Calcul Rho Mu Indicatrice");
  statistiques().begin_count(calculer_rho_mu_indicatrice_counter_);

  const double rho_l = milieu_ijk().get_rho_liquid(), rho_v = milieu_ijk().get_rho_vapour(), mu_l = milieu_ijk().get_mu_liquid(), mu_v = milieu_ijk().get_mu_vapour();

  // En diphasique sans bulle (pour cas tests)
  if (interfaces_->get_nb_bulles_reelles() == 0)
    {
      rho_field_.data() = rho_l;
      rho_field_.echange_espace_virtuel(rho_field_.ghost());
      if (use_inv_rho_)
        {
          inv_rho_field_.data() = 1. / rho_l;
          inv_rho_field_.echange_espace_virtuel(inv_rho_field_.ghost());
        }
      molecular_mu_.data() = mu_l;
      molecular_mu_.echange_espace_virtuel(molecular_mu_.ghost());

      if (parcourir)
        interfaces_->parcourir_maillage();
      return;
    }

  if (parcourir)
    probleme_ijk().parcourir_maillage();

  // Nombre de mailles du domaine NS :
  const int nx = interfaces_->I().ni();
  const int ny = interfaces_->I().nj();
  const int nz = interfaces_->I().nk();
  for (int k = 0; k < nz; k++)
    for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++)
        {
          double chi_l = interfaces_->I(i, j, k);
          rho_field_(i, j, k) = rho_l * chi_l + (1. - chi_l) * rho_v;
          if (harmonic_nu_in_calc_with_indicatrice_ == 1 and chi_l != 0. and chi_l != 1.)
            {
              molecular_mu_(i, j, k) = 1. / (mu_l / chi_l + mu_v / (1. - chi_l));
            }
          else
            {
              molecular_mu_(i, j, k) = mu_l * chi_l + (1. - chi_l) * mu_v;
            }
        }

  if (use_inv_rho_)
    {
      for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
          for (int i = 0; i < nx; i++)
            {
              double chi_l = interfaces_->I(i, j, k);
              inv_rho_field_(i, j, k) = 1. / rho_l * chi_l + (1. - chi_l) * 1. / rho_v;
            }
    }

//  if (use_harmonic_viscosity_)
//    for (int k=0; k < nz; k++)
//      for (int j=0; j < ny; j++)
//        for (int i=0; i < nx; i++)
//          {
//            double chi_l = interfaces_.I(i,j,k);
//            rho_field_(i,j,k)    = rho_liquide_ * chi_l + (1.- chi_l) * rho_vapeur_;
//            molecular_mu_(i,j,k) = (mu_liquide_ * mu_vapeur_) / (chi_l * mu_vapeur_ + (1.- chi_l) * mu_liquide_);
//          }
//  else
//    for (int k=0; k < nz; k++)
//      for (int j=0; j < ny; j++)
//        for (int i=0; i < nx; i++)
//          {
//            double chi_l = interfaces_.I(i,j,k);
//            rho_field_(i,j,k)    = rho_liquide_ * chi_l + (1.- chi_l) * rho_vapeur_;
//            molecular_mu_(i,j,k) = mu_liquide_ * chi_l + (1.- chi_l) * mu_vapeur_ ;
//          }

  //Mise a jour des espaces virtuels des champs :
  rho_field_.echange_espace_virtuel(rho_field_.ghost());
  if (use_inv_rho_)
    inv_rho_field_.echange_espace_virtuel(inv_rho_field_.ghost());
  molecular_mu_.echange_espace_virtuel(molecular_mu_.ghost());

  statistiques().end_count(calculer_rho_mu_indicatrice_counter_);
}

void Navier_Stokes_FTD_IJK::initialise_ns_fields()
{
  Cerr << "Navier_Stokes_FTD_IJK::initialise_ns_fields()" << finl;
  const Probleme_FTD_IJK_base& pb_ijk = probleme_ijk();
  const Domaine_IJK& dom_ijk = pb_ijk.domaine_ijk();

  if (IJK_Shear_Periodic_helpler::defilement_ == 1)
    {
      if (dom_ijk.get_nb_elem_local(2) < 4 )
        {
          Cerr << "Nb of cells / proc in z-direction must be >=4 for shear periodic run" << finl;
          Cerr << "Nb of cells / proc in z-direction must be >=8 for shear periodic run if only one proc on z-direction" << finl;
          Cerr << "Or find an other way to stock indic_ghost_zmin and zmax than IJK_Field" << finl;
          Process::exit();
        }
      if (dom_ijk.get_offset_local(0)!=0.)
        {
          Cerr << " Shear_periodic conditions works only without splitting in i-direction " << finl;
          Cerr << "if splitting in i-direction --> get_neighbour_processor has to be changed" << finl;
          Process::exit();
        }
    }


  allocate_velocity(d_velocity_, dom_ijk, 1, "D_VELOCITY");
  champs_compris_.ajoute_champ_vectoriel(d_velocity_);

  if (test_etapes_et_bilan_)
    {
      auto fields = { &rho_u_euler_av_prediction_champ_, &rho_u_euler_av_rho_mu_ind_champ_, &rho_du_euler_ap_prediction_champ_,
                      &rho_u_euler_ap_projection_champ_, &rho_du_euler_ap_projection_champ_, &rho_u_euler_ap_rho_mu_ind_champ_,
                      &terme_diffusion_local_, &terme_pression_local_, &terme_pression_in_ustar_local_, &d_v_diff_et_conv_, &terme_convection_mass_solver_, &terme_diffusion_mass_solver_
                    };

      for (auto field : fields)
        allocate_velocity(*field, dom_ijk, 1);
    }

  pressure_ghost_cells_.allocate(dom_ijk, Domaine_IJK::ELEM, pb_ijk.get_thermal_probes_ghost_cells());
  pressure_ghost_cells_.data() = 0.;
  pressure_ghost_cells_.echange_espace_virtuel(pressure_ghost_cells_.ghost());

  const double mu_l = milieu_ijk().get_mu_liquid(),
               rho_l = milieu_ijk().get_rho_liquid(),
               mu_v = milieu_ijk().get_mu_vapour(),
               rho_v = milieu_ijk().get_rho_vapour();

  // if interp_monofluide == 2 --> reconstruction uniquement sur rho, mu. Pas sur P !
  if (!Option_IJK::DISABLE_DIPHASIQUE && boundary_conditions_.get_correction_interp_monofluide() == 1)
    {
      pressure_.allocate(dom_ijk, Domaine_IJK::ELEM, 3);
      pressure_.allocate_shear_BC(1, rho_v, rho_l, use_inv_rho_in_poisson_solver_);
    }
  else
    pressure_.allocate(dom_ijk, Domaine_IJK::ELEM, 3);
  pressure_.nommer("PRESSURE");
  pressure_.add_synonymous("PRESSION");
  champs_compris_.ajoute_champ(pressure_);

  if (include_pressure_gradient_in_ustar_)
    {
      d_pressure_.allocate(dom_ijk, Domaine_IJK::ELEM, 1);
      if ( sub_type(Schema_RK3_IJK, schema_temps_ijk()) )
        RK3_F_pressure_.allocate(dom_ijk, Domaine_IJK::ELEM, 1);
    }

  // On utilise aussi rhov pour le bilan de forces et pour d'autres formes de convection...
  allocate_velocity(rho_v_, dom_ijk, 2);

  pressure_rhs_.allocate(dom_ijk, Domaine_IJK::ELEM, 1, "PRESSURE_RHS");
  pressure_rhs_.add_synonymous("PRESSION_RHS");
  champs_compris_.ajoute_champ(pressure_rhs_);


  I_ns_.allocate(dom_ijk, Domaine_IJK::ELEM, 2);
  kappa_ns_.allocate(dom_ijk, Domaine_IJK::ELEM, 2);

  if (!Option_IJK::DISABLE_DIPHASIQUE && (boundary_conditions_.get_correction_interp_monofluide() == 1 || boundary_conditions_.get_correction_interp_monofluide() == 2))
    {
      molecular_mu_.allocate(dom_ijk, Domaine_IJK::ELEM, 2, 0, 1, "VISCOSITY");
      molecular_mu_.allocate_shear_BC(2, mu_v, mu_l);
      rho_field_.allocate(dom_ijk, Domaine_IJK::ELEM, 2, 0, 1, "DENSITY");
      rho_field_.allocate_shear_BC( 2, rho_v, rho_l);
      IJK_Shear_Periodic_helpler::rho_vap_ref_for_poisson_ = rho_v;
      IJK_Shear_Periodic_helpler::rho_liq_ref_for_poisson_ = rho_l;
      if (use_inv_rho_)
        {
          inv_rho_field_.allocate(dom_ijk, Domaine_IJK::ELEM, 2, 0, 1);
          inv_rho_field_.allocate_shear_BC( 2, 1. / rho_v, 1. / rho_l);
          IJK_Shear_Periodic_helpler::rho_vap_ref_for_poisson_ = 1. / rho_v;
          IJK_Shear_Periodic_helpler::rho_liq_ref_for_poisson_ = 1. / rho_l;
        }
    }
  else
    {
      IJK_Shear_Periodic_helpler::rho_vap_ref_for_poisson_ = rho_v;
      IJK_Shear_Periodic_helpler::rho_liq_ref_for_poisson_ = rho_l;
      molecular_mu_.allocate(dom_ijk, Domaine_IJK::ELEM, 2, "VISCOSITY");
      rho_field_.allocate(dom_ijk, Domaine_IJK::ELEM, 2, "DENSITY");
      if (use_inv_rho_)
        {
          inv_rho_field_.allocate(dom_ijk, Domaine_IJK::ELEM, 2);
          IJK_Shear_Periodic_helpler::rho_vap_ref_for_poisson_ = 1. / rho_v;
          IJK_Shear_Periodic_helpler::rho_liq_ref_for_poisson_ = 1. / rho_l;
        }
    }
  rho_field_.nommer("DENSITY");
  rho_field_.add_synonymous("RHO");
  champs_compris_.ajoute_champ(rho_field_);
  molecular_mu_.nommer("VISCOSITY");
  molecular_mu_.add_synonymous("MU");
  champs_compris_.ajoute_champ(molecular_mu_);

  if (schema_temps_ijk().get_first_step_interface_smoothing())
    {
      allocate_velocity(zero_field_ft_, pb_ijk.get_domaine_ft(), pb_ijk.get_thermal_probes_ghost_cells());
      for (int dir = 0; dir < 3; dir++)
        zero_field_ft_[dir].data() = 0.;
      zero_field_ft_.echange_espace_virtuel();
    }

  if (diffusion_alternative_)
    {
      allocate_velocity(laplacien_velocity_, dom_ijk, 1);
      unit_.allocate(dom_ijk, Domaine_IJK::ELEM, 2);
      unit_.data() = 1.;
      unit_.echange_espace_virtuel(unit_.ghost());
    }

  if (velocity_convection_op_.get_convection_op_option_rank() == non_conservative_rhou)
    div_rhou_.allocate(dom_ijk, Domaine_IJK::ELEM, 1);
  allocate_velocity(psi_velocity_, dom_ijk, 2);


  // Allocation du terme source variable spatialement:
  if ((expression_variable_source_[0] != "??") || (expression_variable_source_[1] != "??") || (expression_variable_source_[2] != "??") || (expression_potential_phi_ != "??"))
    {
      allocate_velocity(variable_source_, dom_ijk, 1, "VARIABLE_SOURCE");
      flag_variable_source_ = true;
      potential_phi_.allocate(dom_ijk, Domaine_IJK::ELEM, 1);
      for (int dir = 0; dir < 3; dir++)
        variable_source_[dir].data() = 0.;
      potential_phi_.data() = 0.;
      champs_compris_.ajoute_champ_vectoriel(variable_source_);
    }


  // GB : Je ne sais pas si on a besoin d'un ghost... Je crois que oui. Lequel?
  // Si la a vitesse ft doit transporter les sommets virtuels des facettes reelles,
  // alors il faut un domaine ghost de la taille de la longueur maximale des arretes.
  // allocate_velocity(velocity_ft_, domaine_ft_, 0);
  /*
   * FIXME: Allocate based on the thermal subproblems
   * as the thermal probes necessitates several ghost cells to interpolate velocity !
   * Check the difference between elem and faces ? and for interpolation of the velocity ?
   */
  constexpr int ft_ghost_cells = 4;
  allocate_velocity(velocity_ft_, pb_ijk.get_domaine_ft(), ft_ghost_cells); // named at the end

  kappa_ft_.allocate(pb_ijk.get_domaine_ft(), Domaine_IJK::ELEM, 2, "KAPPA_FT");
  champs_compris_.ajoute_champ(kappa_ft_);

  if (!Option_IJK::DISABLE_DIPHASIQUE)
    {
      auto fields = { &backup_terme_source_interfaces_ft_, &terme_source_interfaces_ns_, &backup_terme_source_interfaces_ns_, &terme_repulsion_interfaces_ns_, &terme_abs_repulsion_interfaces_ns_ };

      for (auto field : fields)
        allocate_velocity(*field, dom_ijk, 1);

      //auto fields_ft = { &terme_source_interfaces_ft_, &terme_repulsion_interfaces_ft_, &terme_abs_repulsion_interfaces_ft_ };
      // for (auto field : fields_ft)
      //  allocate_velocity(*field, pb_ijk.get_domaine_ft(), 1);

      allocate_velocity(terme_source_interfaces_ft_, pb_ijk.get_domaine_ft(), 1, "SOURCE_QDM_INTERF");
      champs_compris_.ajoute_champ_vectoriel(terme_source_interfaces_ft_);
      allocate_velocity(terme_repulsion_interfaces_ft_, pb_ijk.get_domaine_ft(), 1, "SHIELD_REPULSION");
      champs_compris_.ajoute_champ_vectoriel(terme_repulsion_interfaces_ft_);
      allocate_velocity(terme_abs_repulsion_interfaces_ft_, pb_ijk.get_domaine_ft(), 1, "SHIELD_REPULSION_ABS");
      champs_compris_.ajoute_champ_vectoriel(terme_abs_repulsion_interfaces_ft_);
    }

  if ( sub_type(Schema_RK3_IJK, schema_temps_ijk()) )
    {
      Cerr << "Schema temps de type : RK3_FT" << finl;
      allocate_velocity(RK3_F_velocity_, dom_ijk, 1);
    }
  else
    Cerr << "Schema temps de type : Euler_Explicite" << finl;

  velocity_diffusion_op_.initialize(dom_ijk, harmonic_nu_in_diff_operator_);
  velocity_diffusion_op_->set_bc(boundary_conditions_);
  velocity_convection_op_.initialize(dom_ijk);

  // Economise la memoire si pas besoin
  if (!disable_solveur_poisson_)
    poisson_solver_.initialize(dom_ijk);

  // Register champs compris
  velocity_.nommer("VELOCITY");
  velocity_.add_synonymous("VITESSE");
  if (!Option_IJK::DISABLE_DIPHASIQUE)
    {
      velocity_ft_.nommer("VELOCITY_FT");
      velocity_ft_.add_synonymous("VITESSE_FT");
      champs_compris_.ajoute_champ_vectoriel(velocity_ft_);
    }
  champs_compris_.ajoute_champ_vectoriel(velocity_);



}

void Navier_Stokes_FTD_IJK::projeter()
{
  // PROJECTION INITIALE
  // Cette projection n'est pas utile en reprise.
  // Elle sert uniquement a rendre le champ de vitesse initial a divergence nulle lorsque son expression est analytique.
  if (!disable_solveur_poisson_)
    {
      if (improved_initial_pressure_guess_)
        {
          Cerr << "Improved initial pressure" << finl;
          maj_indicatrice_rho_mu();
          if (!Option_IJK::DISABLE_DIPHASIQUE)
            probleme_ijk().update_thermal_properties();

          // La pression n'est pas encore initialisee. elle est donc nulle.
          // Avec cette option, on essaye une initialisation basee sur le champ de pression diphasique
          // a l'equilibre, cad sans vitesse, ou a minima pour un champ a div(u)=0.

          if (!Option_IJK::DISABLE_DIPHASIQUE)
            {
              IJK_Field_vector3_double& coords = probleme_ijk().get_post().coords();

              // Calcul du potentiel.
              for (int dir = 0; dir < 3; dir++)
                {
                  terme_source_interfaces_ft_[dir].data() = 0.;
                  terme_repulsion_interfaces_ft_[dir].data() = 0.;
                  terme_abs_repulsion_interfaces_ft_[dir].data() = 0.;
                }
              const double delta_rho = milieu_ijk().get_delta_rho();
              interfaces_->ajouter_terme_source_interfaces(terme_source_interfaces_ft_, terme_repulsion_interfaces_ft_, terme_abs_repulsion_interfaces_ft_);

              assert(interfaces_->get_nb_bulles_reelles() == 1);
              DoubleTab bounding_box;
              interfaces_->calculer_bounding_box_bulles(bounding_box);
              // Calcul la hauteur en x de la permiere bulle :
              const double Dbx = bounding_box(0, 0, 1) - bounding_box(0, 0, 0);
              const double kappa = 2. / (Dbx / 2.);

              const int ni = pressure_.ni();
              const int nj = pressure_.nj();
              const int nk = pressure_.nk();

              const auto& gravite = milieu_ijk().gravite().valeurs();
              for (int k = 0; k < nk; k++)
                for (int j = 0; j < nj; j++)
                  for (int i = 0; i < ni; i++)
                    {
                      double phi = gravite(0,0) * coords[0](i, j, k) +  gravite(0,1) * coords[1](i, j, k) +  gravite(0,2) * coords[2](i, j, k);
                      double potentiel_elem = milieu_ijk().sigma() * kappa - delta_rho * phi;
                      // La pression est hydrostatique, cad : pressure_ = P - rho g z
                      pressure_(i, j, k) = potentiel_elem * interfaces_->I(i, j, k); // - rho_field_(i,j,k) * phi;
                    }

              // pressure gradient requires the "left" value in all directions:
              pressure_.echange_espace_virtuel(1 /*, IJK_Field_double::EXCHANGE_GET_AT_LEFT_IJK*/);

              // Mise a jour du champ de vitesse (avec dv = seulement le terme source)
              d_velocity_[0].data() = 0.;
              d_velocity_[1].data() = 0.;
              d_velocity_[2].data() = 0.;

              for (int dir = 0; dir < 3; dir++)
                redistribute_from_splitting_ft_faces_[dir].redistribute_add(terme_source_interfaces_ft_[dir], d_velocity_[dir]);

              for (int dir = 0; dir < 3; dir++)
                {
                  const int kmax = d_velocity_[dir].nk();
                  for (int k = 0; k < kmax; k++)
                    euler_explicit_update(d_velocity_[dir], velocity_[dir], k);
                }

              if (use_inv_rho_in_poisson_solver_)
                pressure_projection_with_inv_rho(inv_rho_field_, velocity_[0], velocity_[1], velocity_[2], pressure_, 1., pressure_rhs_, poisson_solver_);
              else
                pressure_projection_with_rho(rho_field_, velocity_[0], velocity_[1], velocity_[2], pressure_, 1., pressure_rhs_, poisson_solver_);
            }
          else
            pressure_projection(velocity_[0], velocity_[1], velocity_[2], pressure_, 1., pressure_rhs_, poisson_solver_);

          copy_field_values(pressure_ghost_cells_, pressure_);
        }
      else if (projection_initiale_demandee_)
        {
          Cerr << "*****************************************************************************\n"
               << "  Attention : projection du champ de vitesse initial sur div(u)=0\n"
               << "*****************************************************************************" << finl;

          if (!Option_IJK::DISABLE_DIPHASIQUE)
            {
              if (use_inv_rho_in_poisson_solver_)
                pressure_projection_with_inv_rho(inv_rho_field_, velocity_[0], velocity_[1], velocity_[2], pressure_, 1., pressure_rhs_, poisson_solver_);
              else
                pressure_projection_with_rho(rho_field_, velocity_[0], velocity_[1], velocity_[2], pressure_, 1., pressure_rhs_, poisson_solver_);
            }
          else
            pressure_projection(velocity_[0], velocity_[1], velocity_[2], pressure_, 1., pressure_rhs_, poisson_solver_);

          pressure_.data() = 0.;
          pressure_rhs_.data() = 0.;
        }
    }

  // Projection initiale sur div(u)=0, si demande: (attention, ne pas le faire en reprise)
  if (correction_semi_locale_volume_bulle_)
    {
      if (disable_solveur_poisson_)
        {
          Cerr << " Warning: Possible incoherence des mots-cles\n"
               << " ===========================================\n"
               << "\n"
               << "Avec correction_semi_locale_volume_bulle, la conservation du volume de la bulle repose entierement sur le fait que la divergence de la vitesse est numeriquement nulle en tout point.\n"
               << "Il est suspect d'utiliser cette option avec disable_solveur_poisson.\n" << finl;
        }
      if (!projection_initiale_demandee_)
        {
          Cerr << " Warning: Possible incoherence des mots-cles\n"
               << " ===========================================\n"
               << "\n"
               << "Avec correction_semi_locale_volume_bulle, la conservation du volume de la bulle repose entierement sur le fait que la divergence de la vitesse est numeriquement nulle en tout point.\n"
               << "Pour securite, il est recommande d'utiliser avec cette option une projection initiale du champ de vitesse, garantissant cette propriete [mot-cle : projection_initiale].\n" << finl;
        }
    }
  Cerr << "End of initial velocity projection" << finl;
}

int Navier_Stokes_FTD_IJK::preparer_calcul()
{
  projeter();

  if (!probleme_ijk().domaine_ijk().get_periodic_flag(DIRECTION_K)) /* Apply BC */
    force_zero_on_walls(velocity_[2]);

  const double mu_l = milieu_ijk().get_mu_liquid(), rho_l = milieu_ijk().get_rho_liquid(), rho_v = milieu_ijk().get_rho_vapour();

  // Si calcul monophasique, on initialise correctement rho, mu, I une fois pour toute :
  if (Option_IJK::DISABLE_DIPHASIQUE)
    {
      rho_field_.data() = rho_l;
      rho_moyen_ = rho_l;
      molecular_mu_.data() = mu_l;
      probleme_ijk().update_thermal_properties();
    }
  else
    {
      Cerr << "Cas normal diphasique Probleme_FTD_IJK_base::run()" << finl;
      probleme_ijk().update_thermal_properties();
      const double indic_moyen = calculer_v_moyen(interfaces_->I());
      rho_moyen_ = indic_moyen * rho_l + (1 - indic_moyen) * rho_v;
      if (probleme_ijk().get_post().is_post_required("EXTERNAL_FORCE"))
        for (int dir = 0; dir < 3; dir++)
          compute_add_external_forces(dir);
    }

  return 1;
}

void Navier_Stokes_FTD_IJK::forcage_control_ecoulement()
{
  update_rho_v(); // Peut-etre pas toujours necessaire selon la formulation pour la convection?
  for (int direction = 0; direction < 3; direction++)
    store_rhov_moy_[direction] = calculer_v_moyen(rho_v_[direction]);
}

void Navier_Stokes_FTD_IJK::initialise_ijk_fields()
{
  Probleme_FTD_IJK_base& pb_ijk = probleme_ijk();
  Cout << "forcage_.get_type_forcage() : " << forcage_.get_type_forcage() << finl;
  if (forcage_.get_type_forcage() > 0)
    {
      const Domaine_IJK& gbz_splitting = velocity_[0].get_domaine();
      const Domaine_IJK& my_geom = velocity_[0].get_domaine();

      const int my_ni = velocity_[0].ni();
      const int my_nj = velocity_[0].nj();
      const int my_nk = velocity_[0].nk();
      const int nproc_tot = Process::nproc();
      Cout << "BF compute_initial_chouippe" << finl;
      Cout << "ni : " << my_ni << " ,nj : " << my_nj << " ,nk : " << my_nk << finl;
      std::cout << "in initialise i_offset : " << gbz_splitting.get_offset_local(DIRECTION_I) << std::endl;
      std::cout << "Process::me()" << Process::me() << std::endl;
      forcage_.compute_initial_chouippe(nproc_tot, my_geom, my_ni, my_nj, my_nk, gbz_splitting, pb_ijk.nom_sauvegarde());




      // TODO (teo.boutin) move this by adding forcage_ to the tree of champs_compris
      champs_compris_.ajoute_champ_vectoriel(forcage_.get_force_ph2());


      statistiques().begin_count(m2_counter_);
      Cout << "AF compute_initial_chouippe" << finl;
    }

  if (fichier_reprise_vitesse_ == "??")   // si on ne fait pas une reprise on initialise V
    initialise_velocity_using_expression(expression_vitesse_initiale_);
  else
    initialise_velocity_from_file(fichier_reprise_vitesse_);

  // Pour le check_stats_ ou pour travailler en increment de pression, il faut connaitre la pression initiale :
  if (expression_pression_initiale_ != "??")
    {
      Cout << "Initialisation pression \nPini = " << expression_pression_initiale_ << finl;
      set_field_data(pressure_, expression_pression_initiale_);
      pressure_.echange_espace_virtuel(pressure_.ghost());
    }

  /*
   * Compute mean rho_g using the indicator function
   */
  if (compute_force_init_)
    {
      terme_source_acceleration_ = 0.;
      double indicator_sum = 0;
      double terme_source_acceleration_increment;
      const int ni = interfaces_->I().ni();
      const int nj = interfaces_->I().nj();
      const int nk = interfaces_->I().nk();
      const double rho_l = milieu_ijk().get_rho_liquid(), rho_v = milieu_ijk().get_rho_vapour();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double indic = interfaces_->I(i, j, k);
              indicator_sum += indic;
              terme_source_acceleration_increment = rho_l * indic + rho_v * (1 - indic);
              terme_source_acceleration_ += terme_source_acceleration_increment;
            }
      indicator_sum = Process::mp_sum(indicator_sum);
      terme_source_acceleration_ = Process::mp_sum(terme_source_acceleration_);
      const double gravite_norm = milieu_ijk().get_gravite_norm();
      terme_source_acceleration_ /= indicator_sum;
      terme_source_acceleration_ *= gravite_norm;
      Cout << "Calculation of force init due to gravity mean(rho_g): " << terme_source_acceleration_ << finl;
    }

  // statistiques...
  pb_ijk.get_post().initialise_stats(pb_ijk.domaine_ijk(), vol_bulles_, vol_bulle_monodisperse_);

  if (coef_immobilisation_ > 1e-16)
    {
      allocate_velocity(force_rappel_, pb_ijk.domaine_ijk(), 2, "EXTERNAL_FORCE");
      champs_compris_.ajoute_champ_vectoriel(force_rappel_);
      allocate_velocity(force_rappel_ft_, pb_ijk.get_domaine_ft(), 2);
      // A la reprise, c'est fait par le IJK_Interfaces::readOn
      if (interfaces_->get_flag_positions_reference() == 0) // (!reprise_)
        {
          Cerr << "Saving interfacial positions as references." << finl;
          interfaces_->set_positions_reference();
        }
    }

  // Maj des grandeurs shear perio
  if (IJK_Shear_Periodic_helpler::defilement_ == 1)
    {
      IJK_Shear_Periodic_helpler::shear_x_time_ = boundary_conditions_.get_dU_perio() * (schema_temps_ijk().get_current_time() + boundary_conditions_.get_t0_shear());
      redistribute_from_splitting_ft_elem_ghostz_min_.redistribute(interfaces_->I_ft(), I_ns_);
      redistribute_from_splitting_ft_elem_ghostz_max_.redistribute(interfaces_->I_ft(), I_ns_);
      rho_field_.get_shear_BC_helpler().set_indicatrice_ghost_zmin_(I_ns_, 0);
      rho_field_.get_shear_BC_helpler().set_indicatrice_ghost_zmax_(I_ns_, rho_field_.nk() - 4);
      molecular_mu_.get_shear_BC_helpler().set_indicatrice_ghost_zmin_(I_ns_, 0);
      molecular_mu_.get_shear_BC_helpler().set_indicatrice_ghost_zmax_(I_ns_, molecular_mu_.nk() - 4);
      if (use_inv_rho_)
        {
          inv_rho_field_.get_shear_BC_helpler().set_indicatrice_ghost_zmin_(I_ns_, 0);
          inv_rho_field_.get_shear_BC_helpler().set_indicatrice_ghost_zmax_(I_ns_, inv_rho_field_.nk() - 4);
        }
      if (boundary_conditions_.get_correction_interp_monofluide() == 1)
        {
          interfaces_->calculer_kappa_ft(kappa_ft_);
          redistribute_from_splitting_ft_elem_ghostz_min_.redistribute(kappa_ft_, kappa_ns_);
          redistribute_from_splitting_ft_elem_ghostz_max_.redistribute(kappa_ft_, kappa_ns_);
          pressure_.get_shear_BC_helpler().set_I_sig_kappa_zmin_(I_ns_, kappa_ns_, milieu_ijk().sigma(), 0);
          pressure_.get_shear_BC_helpler().set_I_sig_kappa_zmax_(I_ns_, kappa_ns_, milieu_ijk().sigma(), pressure_.nk() - 4);
        }
    }

  maj_indicatrice_rho_mu();
}

void Navier_Stokes_FTD_IJK::complete_initialise_ijk_fields()
{
  if (IJK_Shear_Periodic_helpler::defilement_ == 1)
    {
      IJK_Shear_Periodic_helpler::shear_x_time_ = boundary_conditions_.get_dU_perio() * (schema_temps_ijk().get_current_time() + boundary_conditions_.get_t0_shear());
      redistribute_from_splitting_ft_elem_ghostz_min_.redistribute(interfaces_->I_ft(), I_ns_);
      redistribute_from_splitting_ft_elem_ghostz_max_.redistribute(interfaces_->I_ft(), I_ns_);
      rho_field_.get_shear_BC_helpler().set_indicatrice_ghost_zmin_(I_ns_, 0);
      rho_field_.get_shear_BC_helpler().set_indicatrice_ghost_zmax_(I_ns_, rho_field_.nk() - 4);
      molecular_mu_.get_shear_BC_helpler().set_indicatrice_ghost_zmin_(I_ns_, 0);
      molecular_mu_.get_shear_BC_helpler().set_indicatrice_ghost_zmax_(I_ns_, molecular_mu_.nk() - 4);
      if (use_inv_rho_)
        {
          inv_rho_field_.get_shear_BC_helpler().set_indicatrice_ghost_zmin_(I_ns_, 0);
          inv_rho_field_.get_shear_BC_helpler().set_indicatrice_ghost_zmax_(I_ns_, inv_rho_field_.nk() - 4);
        }
      if (boundary_conditions_.get_correction_interp_monofluide() == 1)
        {
          interfaces_->calculer_kappa_ft(kappa_ft_);
          redistribute_from_splitting_ft_elem_ghostz_min_.redistribute(kappa_ft_, kappa_ns_);
          redistribute_from_splitting_ft_elem_ghostz_max_.redistribute(kappa_ft_, kappa_ns_);
          pressure_.get_shear_BC_helpler().set_I_sig_kappa_zmin_(I_ns_, kappa_ns_, milieu_ijk().sigma(), 0);
          pressure_.get_shear_BC_helpler().set_I_sig_kappa_zmax_(I_ns_, kappa_ns_, milieu_ijk().sigma(), pressure_.nk() - 4);
        }
    }
}

void Navier_Stokes_FTD_IJK::redistribute_to_splitting_ft_elem(const IJK_Field_double& input_field, IJK_Field_double& output_field)
{
  redistribute_to_splitting_ft_elem_.redistribute(input_field, output_field);
}

void Navier_Stokes_FTD_IJK::redistribute_from_splitting_ft_elem(const IJK_Field_double& input_field, IJK_Field_double& output_field)
{
  redistribute_from_splitting_ft_elem_.redistribute(input_field, output_field);
}

void Navier_Stokes_FTD_IJK::update_rho_v()
{
  if (use_inv_rho_for_mass_solver_and_calculer_rho_v_)
    {
      // rho_face = 2*(rho_gauche*rho_droite)/(rho_gauche+rho_droite)
      //          = 1./ (1/2 * (1/rho_g + 1/rho_d))
      // 1/rho_face est donc la moyenne geometrique de inv_rho.
      calculer_rho_harmonic_v(rho_field_, velocity_, rho_v_);
    }
  else
    calculer_rho_v(rho_field_, velocity_, rho_v_);

  rho_v_[0].echange_espace_virtuel(rho_v_[0].ghost());
  rho_v_[1].echange_espace_virtuel(rho_v_[1].ghost());
  rho_v_[2].echange_espace_virtuel(rho_v_[2].ghost());
}

void Navier_Stokes_FTD_IJK::calculer_terme_asservissement(double& ax, double& ay, double& az)
{
  // On trouve la vitesse moyenne de la phase vapeur pour la partie derivee du correcteur

  update_rho_v();
  double v_moyx = calculer_v_moyen(velocity_[0]);
  double v_moyy = calculer_v_moyen(velocity_[1]);
  double v_moyz = calculer_v_moyen(velocity_[2]);

  double rhov_moyx = calculer_v_moyen(rho_v_[DIRECTION_I]);
  double rhov_moyy = calculer_v_moyen(rho_v_[DIRECTION_J]);
  double rhov_moyz = calculer_v_moyen(rho_v_[DIRECTION_K]);

  const Domaine_IJK& geom = velocity_[milieu_ijk().get_direction_gravite()].get_domaine();
  double Lz = geom.get_domain_length(DIRECTION_K);
  double Lx = geom.get_domain_length(DIRECTION_I);
  double Ly = geom.get_domain_length(DIRECTION_J);

  double vol_dom = Lz * Lx * Ly;
  double alv = 0.;
  if (vol_bulle_monodisperse_ >= 0.)
    alv = interfaces_->get_nb_bulles_reelles() * vol_bulle_monodisperse_ / vol_dom;
  else
    alv = 1. - calculer_v_moyen(interfaces_->I());

  const double drho = milieu_ijk().get_delta_rho(), rho_l = milieu_ijk().get_rho_liquid();
  double facv = 0.;
  if (std::fabs(alv * drho) > DMINFLOAT)
    {
      facv = 1. / (alv * drho);
    }
  double uvx = facv * (rho_l * v_moyx - rhov_moyx);
  double uvy = facv * (rho_l * v_moyy - rhov_moyy);
  double uvz = facv * (rho_l * v_moyz - rhov_moyz);

  // On evalue la position de chaque bulles pour trouver le barycentre de la phase vapeur

  ArrOfDouble volumes;
  DoubleTab centre_gravite;

  const int nbulles_reelles = interfaces_->get_nb_bulles_reelles();
  const int nbulles_ghost = interfaces_->get_nb_bulles_ghost();
  const int nbulles_tot = nbulles_reelles + nbulles_ghost;

  volumes.resize_array(nbulles_tot);
  volumes = 0.;
  centre_gravite.resize(nbulles_tot, 3);
  centre_gravite = 0.;

  double centre_moyx = 0;
  double centre_moyy = 0;
  double centre_moyz = 0;

  interfaces_->calculer_volume_bulles(volumes, centre_gravite);

  for (int i = 0; i < nbulles_tot; i++)
    {
      centre_moyx += 1.0 * centre_gravite(i, 0) / nbulles_tot;
      centre_moyy += 1.0 * centre_gravite(i, 1) / nbulles_tot;
      centre_moyz += 1.0 * centre_gravite(i, 2) / nbulles_tot;
    }

  // On met a jour l'integrale du deplacement du barycentre

  const double timestep = schema_temps_ijk().get_timestep();
  int_x_ += (centre_moyx - Lx / 2) * timestep;
  int_y_ += (centre_moyy - Ly / 2) * timestep;
  int_z_ += (centre_moyz - Lz / 2) * timestep;

  ax = -Kp_ * (centre_moyx - Lx / 2) - Ki_ * int_x_ - Kd_ * uvx;
  ay = -Kp_ * (centre_moyy - Ly / 2) - Ki_ * int_y_ - Kd_ * uvy;
  az = -Kp_ * (centre_moyz - Lz / 2) - Ki_ * int_z_ - Kd_ * uvz;
}

void Navier_Stokes_FTD_IJK::update_v_ghost_from_rho_v()
{
  const Domaine_IJK& dom = probleme_ijk().domaine_ijk();

  for (int dir = 0; dir < 3; dir++)
    {
      const int imax = velocity_[dir].ni();
      const int jmax = velocity_[dir].nj();
      const int kmax = velocity_[dir].nk();
      const int ghost = velocity_[dir].ghost();
      const int last_global_k = dom.get_nb_items_global(Domaine_IJK::ELEM, 2);

      for (int j = 0; j < jmax; j++)
        for (int i = 0; i < imax; i++)
          {
            if (dom.get_offset_local(2) == 0)
              {
                for (int k = -ghost; k < 0; k++)
                  {
                    double rho = 0.;
                    double DU = 0.;
                    if (dir == 0)
                      {
                        rho = 0.5 * (rho_field_(i, j, k) + rho_field_(i - 1, j, k));
                        DU = boundary_conditions_.get_dU_perio();
                      }
                    else if (dir == 1)
                      rho = 0.5 * (rho_field_(i, j, k) + rho_field_(i, j - 1, k));
                    else if (dir == 2)
                      rho = 0.5 * (rho_field_(i, j, k) + rho_field_(i, j, k - 1));

                    velocity_[dir](i, j, k) = rho_v_[dir](i, j, k) / rho - DU;
                  }
              }
            if (dom.get_offset_local(2) + kmax == last_global_k)
              {
                for (int k = kmax; k < kmax + ghost; k++)
                  {
                    double rho = 0.;
                    double DU = 0.;
                    if (dir == 0)
                      {
                        rho = 0.5 * (rho_field_(i, j, k) + rho_field_(i - 1, j, k));
                        DU = boundary_conditions_.get_dU_perio();
                      }
                    else if (dir == 1)
                      rho = 0.5 * (rho_field_(i, j, k) + rho_field_(i, j - 1, k));
                    else if (dir == 2)
                      rho = 0.5 * (rho_field_(i, j, k) + rho_field_(i, j, k - 1));

                    velocity_[dir](i, j, k) = rho_v_[dir](i, j, k) / rho + DU;
                  }
              }
          }
    }
}

// Transfert du maillage ft vers ns de champs aux faces :
void Navier_Stokes_FTD_IJK::transfer_ft_to_ns()
{
  for (int dir = 0; dir < 3; dir++)
    {
      redistribute_from_splitting_ft_faces_[dir].redistribute(terme_repulsion_interfaces_ft_[dir], terme_repulsion_interfaces_ns_[dir]);
      redistribute_from_splitting_ft_faces_[dir].redistribute(terme_abs_repulsion_interfaces_ft_[dir], terme_abs_repulsion_interfaces_ns_[dir]);
    }
}

void Navier_Stokes_FTD_IJK::calculer_vitesse_droite(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz, double& vx_moy, double& vy_moy, double& vz_moy)
{
  /* Renvoie le vecteur vitesse moyen (spatial) en z = 0 */
  /* Ne fonctionne que pour des maillages uniformes */
  const Domaine_IJK& splitting = vx.get_domaine();
  const int ni = vx.ni();
  const int nj = vx.nj();
  const int nk = vx.nk();
  //const int nk = vx.nk();
  //double dz = splitting.get_constant_delta(DIRECTION_K);
  int z_index = splitting.get_local_slice_index(2);
  int z_index_max = splitting.get_nprocessor_per_direction(2) - 1;
  vx_moy = 0.;
  vy_moy = 0.;
  vz_moy = 0.;
  double alpha_l_moy = 0.;

  if (z_index == z_index_max)
    {
      for (int i = 0; i < ni; i++)
        {
          for (int j = 0; j < nj; j++)
            {
              vx_moy += vx(i, j, nk - 1) * interfaces_->I(i, j, nk - 1);
              vy_moy += vy(i, j, nk - 1) * interfaces_->I(i, j, nk - 1);
              vz_moy += vz(i, j, nk - 1) * interfaces_->I(i, j, nk - 1);
              alpha_l_moy += interfaces_->I(i, j, nk - 1);
            }
        }
    }
  vx_moy = Process::mp_sum(vx_moy);
  vy_moy = Process::mp_sum(vy_moy);
  vz_moy = Process::mp_sum(vz_moy);
  alpha_l_moy = Process::mp_sum(alpha_l_moy);

  vx_moy = vx_moy / alpha_l_moy;
  vy_moy = vy_moy / alpha_l_moy;
  vz_moy = vz_moy / alpha_l_moy;

  return;
}

void Navier_Stokes_FTD_IJK::calculer_vitesse_gauche(const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz, double& vx_moy, double& vy_moy, double& vz_moy)
{
  /* Renvoie le vecteur vitesse moyen (spatial) en z = 0 */
  /* Ne fonctionne que pour des maillages uniformes */
  const Domaine_IJK& splitting = vx.get_domaine();
  const int ni = vx.ni();
  const int nj = vx.nj();
  //const int nk = vx.nk();
  //double dz = splitting.get_constant_delta(DIRECTION_K);
  int z_index = splitting.get_local_slice_index(2);
  int z_index_min = 0;
  vx_moy = 0.;
  vy_moy = 0.;
  vz_moy = 0.;
  double alpha_l_moy = 0.;

  if (z_index == z_index_min)
    {
      for (int i = 0; i < ni; i++)
        {
          for (int j = 0; j < nj; j++)
            {
              vx_moy += vx(i, j, 0) * interfaces_->I(i, j, 0);
              vy_moy += vy(i, j, 0) * interfaces_->I(i, j, 0);
              vz_moy += vz(i, j, 0) * interfaces_->I(i, j, 0);
              alpha_l_moy += interfaces_->I(i, j, 0);
            }
        }
    }

  vx_moy = Process::mp_sum(vx_moy);
  vy_moy = Process::mp_sum(vy_moy);
  vz_moy = Process::mp_sum(vz_moy);
  alpha_l_moy = Process::mp_sum(alpha_l_moy);

  vx_moy = vx_moy / alpha_l_moy;
  vy_moy = vy_moy / alpha_l_moy;
  vz_moy = vz_moy / alpha_l_moy;
  return;
}

static double calculer_tau_wall(const IJK_Field_double& vx, const double mu_liquide)
{
  const int nj = vx.nj();
  const int ni = vx.ni();
  const Domaine_IJK& geom = vx.get_domaine();
  // Maillage uniforme, il suffit donc de diviser par le nombre total de mailles:
  // cast en double au cas ou on voudrait faire un maillage >2 milliards
  const double n_mailles_plan_xy = ((double) geom.get_nb_elem_tot(0)) * geom.get_nb_elem_tot(1);
  const int kmin = vx.get_domaine().get_offset_local(DIRECTION_K);
  const Domaine_IJK::Localisation loc = vx.get_localisation();
  const int nktot = vx.get_domaine().get_nb_items_global(loc, DIRECTION_K);
  double tauw = 0.;

#ifndef VARIABLE_DZ
  const double dz = geom.get_constant_delta(DIRECTION_K);
#else
  const ArrOfDouble& tab_dz=geom.get_delta(DIRECTION_K);
  //  const int nz = geom.get_nb_elem_tot(DIRECTION_K);
  double dz = -1.; // invalid value
  // On prend le dz sur la paroi basse :
  if (kmin == 0)
    dz = tab_dz[0];
#endif
  if (kmin == 0)
    {
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
          tauw += vx(i, j, 0);
    }
  if (kmin + vx.nk() == nktot)
    {
      const int k = vx.nk() - 1;
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
          tauw += vx(i, j, k);
    }

  // somme sur tous les processeurs.
  tauw = Process::mp_sum(tauw);
  // Il faut diviser par dz/2. et par 2*n_mailles_plan_xy
#ifdef VARIABLE_DZ
  // Pour definir un dz sur tous les procs, on recupere celui de la paroi basse
  // sur tous les procs, y compris ceux qui n'ont pas de mur, ou pas le mur gauche.
  dz = Process::mp_max(dz);
#endif
  tauw /= (dz * n_mailles_plan_xy);
  tauw *= mu_liquide;
  return tauw;
}

// Inspiree de la methode d'IJK_Navier_Stokes_Tools :
static void runge_kutta3_update_for_float(const double dx, double& store, double& v, const int step, double dt_tot)
{
  const double coeff_a[3] = { 0., -5. / 9., -153. / 128. };
  // Fk[0] = 1; Fk[i+1] = Fk[i] * a[i+1] + 1
  const double coeff_Fk[3] = { 1., 4. / 9., 15. / 32. };

  const double facteurF = coeff_a[step];
  const double intermediate_dt = compute_fractionnal_timestep_rk3(dt_tot, step);
  const double delta_t_divided_by_Fk = intermediate_dt / coeff_Fk[step];
  double x;
  switch(step)
    {
    case 0:
      x = dx;
      store = x;
      v += x * delta_t_divided_by_Fk;
      break;
    case 1:
      // general case, read and write F
      x = store * facteurF + dx;
      store = x;
      v += x * delta_t_divided_by_Fk;
      break;
    case 2:
      // do not write F
      x = store * facteurF + dx;
      v += x * delta_t_divided_by_Fk;
      break;
    default:
      Cerr << "Error in runge_kutta_update_for_float: wrong step" << finl;
      Process::exit();
    };
}

void Navier_Stokes_FTD_IJK::test_etapes_et_bilan_rho_u_euler(bool apres)
{
  if (test_etapes_et_bilan_)
    {
      if (apres)
        {
          calculer_rho_v(rho_field_, velocity_, rho_u_euler_ap_rho_mu_ind_champ_);
          for (int dir = 0; dir < 3; dir++)
            {
              rho_u_euler_ap_rho_mu_ind_[dir] = calculer_v_moyen(rho_u_euler_ap_rho_mu_ind_champ_[dir]);
              u_euler_ap_rho_mu_ind_[dir] = calculer_v_moyen(velocity_[dir]);
            }
        }
      else /* avant */
        {
          calculer_rho_v(rho_field_, velocity_, rho_u_euler_av_rho_mu_ind_champ_);
          for (int dir = 0; dir < 3; dir++)
            rho_u_euler_av_rho_mu_ind_[dir] = calculer_v_moyen(rho_u_euler_av_rho_mu_ind_champ_[dir]);
        }
    }
}

void Navier_Stokes_FTD_IJK::create_forced_dilation()
{
  if (coeff_evol_volume_ != 0.)
    {
      vol_bulle_monodisperse_ = vol_bulle_monodisperse_ * (1. + coeff_evol_volume_ * schema_temps_ijk().get_timestep());
      const int nb_reelles = interfaces_->get_nb_bulles_reelles();
      for (int ib = 0; ib < nb_reelles; ib++)
        vol_bulles_[ib] = vol_bulle_monodisperse_;
    }
}

void Navier_Stokes_FTD_IJK::calculer_terme_source_acceleration(const double time, const double timestep, const int rk_step, const int dir)
{
  calculer_terme_source_acceleration(velocity_[dir], time, timestep, rk_step);
}

void Navier_Stokes_FTD_IJK::calculer_terme_source_acceleration(IJK_Field_double& vx, const double time, const double timestep, const int rk_step)
{
  /*
   * Cette methode calcule la source de qdm a appliquer pour que les bulles soient
   * globalement fixes dans la direction de la gravite
   *  o Si le parametre source_qdm_gr vaut 1, la correction a appliquer est :
   *          vx = vx - moy^{xyz} (rho vx) - moy^{xyz} (rho vx_terminale) / moy^{xyz} (rho)
   *  o Si le parametre source_qdm_gr vaut 0, la source est determinee par :
   *          temre_force_init         --> temre_source_acceleration et par
   *          expression_derivee_force --> expression_derivee_acceleration
   * REMARQUE II : On peut envisager de faire une correction qui n'a pas besoin qu'on lui donne vx_terminale en
   *               entree. C'est ce qui a ete explore, mais qui n'a pas aboutit.
   *  */
  statistiques().begin_count(source_counter_);
  double new_time = time;
  double v_moy = calculer_v_moyen(vx);

  // S'il n'y a pas de derivee, la source est constante donc on peut sortir:
  if (expression_derivee_acceleration_ == Nom("0"))
    {
      //    terme_source_acceleration_ = 0.;
      statistiques().end_count(source_counter_);
      return;
    }

  update_rho_v();
  const int direction_gravite = milieu_ijk().get_direction_gravite();

  // GAB, rotation : pas sur de mon coup la
  // double rhov_moy = calculer_v_moyen(rho_v_[DIRECTION_I]);
  double rhov_moy = calculer_v_moyen(rho_v_[direction_gravite]);

  double moy_rappel = 0.;
  // GAB, rotation
  const Domaine_IJK& geom = velocity_[direction_gravite].get_domaine();
  double vol_dom = geom.get_domain_length(DIRECTION_I) * geom.get_domain_length(DIRECTION_J) * geom.get_domain_length(DIRECTION_K);
  if (coef_immobilisation_ > 1e-16)
    {
      // GAB, rotation, pas sur de mon coup
      // moy_rappel=calculer_v_moyen(force_rappel_[DIRECTION_I])*vol_dom;
      moy_rappel = calculer_v_moyen(force_rappel_[direction_gravite]) * vol_dom;
    }

  const double rho_l = milieu_ijk().get_rho_liquid(), rho_v = milieu_ijk().get_rho_vapour(), mu_l = milieu_ijk().get_mu_liquid();
  double tauw = calculer_tau_wall(vx, mu_l);
  double derivee_acceleration = 0.;
  double derivee_facteur_sv = 0.;

  double alv = 0.;
  if (probleme_ijk().has_interface())
    {
      if (vol_bulle_monodisperse_ >= 0.)
        alv = interfaces_->get_nb_bulles_reelles() * vol_bulle_monodisperse_ / vol_dom;
      else
        alv = 1. - calculer_v_moyen(interfaces_->I());
    }
  if (Process::je_suis_maitre())
    {
      //
      double drho = rho_l - rho_v;
      double facv = 0., facl = 1.;
      if (std::fabs(alv * drho) > DMINFLOAT)
        {
          facv = 1. / (alv * drho);
          facl = 1. / ((1. - alv) * drho);
        }
      double ul = facl * (rhov_moy - rho_v * v_moy);
      double uv = facv * (rho_l * v_moy - rhov_moy);

      // Mise a jour de l'acceleration
      parser_derivee_acceleration_.setVar("rappel_moyen", moy_rappel);
      parser_derivee_acceleration_.setVar("force", terme_source_acceleration_);
      parser_derivee_acceleration_.setVar("v_moyen", v_moy);
      parser_derivee_acceleration_.setVar("ur", uv - ul);
      parser_derivee_acceleration_.setVar("ul", ul);
      parser_derivee_acceleration_.setVar("uv", uv);
      parser_derivee_acceleration_.setVar("T", time);
      parser_derivee_acceleration_.setVar("rhov_moyen", rhov_moy);
      parser_derivee_acceleration_.setVar("tauw", tauw);
      // Pour utiliser rho_v il faudrait deplacer cette mise a jour a un endroit ou rho
      // est a jour en fonction de l'indicatrice
      // parser_derivee_acceleration_.setVar("rho_v_moyen", rho_v_moy);
      derivee_acceleration = parser_derivee_acceleration_.eval();

      // Mise a jour de la source variable
      if (expression_derivee_facteur_variable_source_ != Nom("0"))
        {
          parser_derivee_facteur_variable_source_.setVar("rappel_moyen", moy_rappel);
          parser_derivee_facteur_variable_source_.setVar("facteur_sv", facteur_variable_source_);
          parser_derivee_facteur_variable_source_.setVar("v_moyen", v_moy);
          parser_derivee_facteur_variable_source_.setVar("T", time);
          parser_derivee_facteur_variable_source_.setVar("rhov_moyen", rhov_moy);
          parser_derivee_facteur_variable_source_.setVar("tauw", tauw);
          derivee_facteur_sv = parser_derivee_facteur_variable_source_.eval();
        }      //

      if (qdm_corrections_.is_type_gb())
        {
          // calcul de terme_source_acceleration_ et de terme_source_acceleration_

          // ON NE VEUT PAS METTRE A JOUR TERME_SOURCE_ACCELERATION_ AVEC CETTE METHODE
          if (sub_type(Schema_Euler_explicite_IJK, schema_temps_ijk()))
            {
              terme_source_acceleration_ += derivee_acceleration * timestep;
              //terme_source_acceleration_ += 0;//derivee_acceleration * timestep;
              facteur_variable_source_ += derivee_facteur_sv * timestep;
              //facteur_variable_source_ += 0;//derivee_facteur_sv * timestep;
              new_time += timestep;
            }
          else if (sub_type(Schema_RK3_IJK, schema_temps_ijk()))
            {
              const double intermediate_dt = compute_fractionnal_timestep_rk3(timestep, rk_step);
              Schema_RK3_IJK& rk3 = ref_cast(Schema_RK3_IJK, schema_temps_ijk());
              runge_kutta3_update_for_float(derivee_acceleration, rk3.get_store_RK3_source_acc(), terme_source_acceleration_, rk_step, timestep);
              Cout << "terme_source_acceleration_" << terme_source_acceleration_ << finl;
              //terme_source_acceleration_ += 0;
              runge_kutta3_update_for_float(derivee_facteur_sv, rk3.get_store_RK3_fac_sv(), facteur_variable_source_, rk_step, timestep);
              Cout << "facteur_variable_source_" << facteur_variable_source_ << finl;
              //facteur_variable_source_ += 0;
              new_time += intermediate_dt;
            }
        }
      else if (sub_type(Schema_RK3_IJK, schema_temps_ijk()))
        {
          const double intermediate_dt = compute_fractionnal_timestep_rk3(timestep/*total */, rk_step);
          Schema_RK3_IJK& rk3 = ref_cast(Schema_RK3_IJK, schema_temps_ijk());
          runge_kutta3_update_for_float(derivee_acceleration, rk3.get_store_RK3_source_acc(), terme_source_acceleration_, rk_step, timestep/*total */);

          runge_kutta3_update_for_float(derivee_facteur_sv, rk3.get_store_RK3_fac_sv(), facteur_variable_source_, rk_step, timestep/*total */);
          new_time += intermediate_dt;
        }
    }
  envoyer_broadcast(terme_source_acceleration_, 0);

  // -----------------------------------------------------------
  // Force interface (:"force sigma") seloon x,y,z et u.Force_interface
  const Probleme_FTD_IJK_base& pb_ijk = probleme_ijk();
  double fs0(0), fs1(0), fs2(0), psn(0);
  if (!Option_IJK::DISABLE_DIPHASIQUE)
    {
      // FORCE INTERFACIALE : on veut un terme homogene a [rho.g]=[N.m^{-3}]
      // terme_source_interfaces_ns_       est homogene a [du/dt]=[m.s^{-2}]
      //   -> in calculer_dv : "~ velocity_ += force_interf * dt ~"
      //   --> force_interf a bien eu un mass_solver_with_rho plus haut
      fs0 = calculer_v_moyen(scalar_fields_product(pb_ijk, rho_field_, terme_source_interfaces_ns_[0], 0));
      fs1 = calculer_v_moyen(scalar_fields_product(pb_ijk, rho_field_, terme_source_interfaces_ns_[1], 1));
      fs2 = calculer_v_moyen(scalar_fields_product(pb_ijk, rho_field_, terme_source_interfaces_ns_[2], 2));
      psn = calculer_v_moyen(scalar_product(pb_ijk, velocity_, scalar_times_vector(pb_ijk, rho_field_, terme_source_interfaces_ns_)));
    }
  // energie cinetique (monophasique) et diphasique
  double uu(calculer_v_moyen(scalar_product(pb_ijk, velocity_, velocity_)));
  double uru(calculer_v_moyen(scalar_product(pb_ijk, velocity_, scalar_times_vector(pb_ijk, rho_field_, velocity_))));
  // Force exterieur (:"force thi") selon x,y,z et acceleration_thi.acceleration_thi, force_thi.foce_thi, u.Force_THI
  double ft0(0), ft1(0), ft2(0), atat(0), ftft(0), ptn(0);
  if (forcage_.get_type_forcage() > 0)
    {
      // FORCE IMPOSEE : on veut un terme homogene a [rho.g]=[N.m^{-3}]
      // forcage_.get_force_ph2()     est homogene a [du/dt]=[m.s^{-2}]
      //   -> in compute_add_THI_force_sur_d_velocity : "~ d_velocity += forcage_.get_force_ph2() ~"
      ft0 = calculer_v_moyen(scalar_fields_product(pb_ijk, rho_field_, forcage_.get_force_ph2()[0], 0));
      ft1 = calculer_v_moyen(scalar_fields_product(pb_ijk, rho_field_, forcage_.get_force_ph2()[1], 1));
      ft2 = calculer_v_moyen(scalar_fields_product(pb_ijk, rho_field_, forcage_.get_force_ph2()[2], 2));
      atat = calculer_v_moyen(scalar_product(pb_ijk, forcage_.get_force_ph2(), forcage_.get_force_ph2()));
      ftft = calculer_v_moyen(scalar_product(pb_ijk, scalar_times_vector(pb_ijk, rho_field_, forcage_.get_force_ph2()), scalar_times_vector(pb_ijk, rho_field_, forcage_.get_force_ph2())));
      ptn = calculer_v_moyen(scalar_product(pb_ijk, velocity_, scalar_times_vector(pb_ijk, rho_field_, forcage_.get_force_ph2())));
    }
  // -----------------------------------------------------------

  // Impression dans le fichier _acceleration.out
  if (Process::je_suis_maitre())
    {
      // GR : 07.01.22 : ce serai pas mal de mettre une condition if (tstep % dt_post_stats_acc_ == dt_post_stats_acc_ - 1 || stop)
      //      pour alleger le dossier OUT. Voir avec GB et AB.
      // double ff=0.;
      int reset = (!pb_ijk.get_reprise()) && (schema_temps_ijk().get_tstep() == 0);
      SFichier fic =
        Ouvrir_fichier("_acceleration.out",
                       "1.tstep\t2.time\t3.Vx\t4.rhoVx\t5.tauw\t6.da/dt\t7.NewT\t8.acceleration\t9.fac_var_source\t10.qdm_source\t11.vap_velocity_tmoy_\t12.liq_velocity_tmoy_\t13.qdm_patch_correction_[0]\t14.qdm_patch_correction_[1]\t15.qdm_patch_correction_[2]\t16.F_sigma_moyen[0]\t17.F_sigma_moyen[1]\t18.F_sigma_moyen[2]\t19.y.F_sigma\t20.u.u\t21.F_THI[0]\t22.F_THI[1]\t23.F_THI[2]\t24.A_THI.A_THI\t25.F_THI.F_THI\t26.u.F_THI\t27.u.rho.u",
                       reset);
      // la derivee_acceleration n'est connue que sur le maitre
      fic << schema_temps_ijk().get_tstep() << " " << time << " " << v_moy << " " << rhov_moy << " " << tauw;
      fic << " " << derivee_acceleration << " " << new_time << " " << terme_source_acceleration_;
      fic << " " << facteur_variable_source_;
      fic << " " << 0.; //qdm_source_;
      fic << " " << vap_velocity_tmoy_;
      fic << " " << liq_velocity_tmoy_;

      if (coef_immobilisation_ > 1e-16)
        fic << " " << moy_rappel;

      for (int dir = 0; dir < 3; dir++)
        fic << " " << 0.; //qdm_patch_correction_[dir];

      // Force interfaciale et puissance du travail des forces interfaciales
      // rho*terme_source_interfaces_ns_
      fic << " " << fs0; // F_sigma_moyen[0]
      fic << " " << fs1; // F_sigma_moyen[1]
      fic << " " << fs2; // F_sigma_moyen[2]
      // u.rho*terme_source_interfaces_ns_
      fic << " " << psn; // velocity.F_sigma

      // Energie cinetique (double)
      // u.u (qui est aussi accessible par les .txt)
      fic << " " << uu;

      // Force imposee et puissance du traveil de la force imposee
      // rho*F_THI
      fic << " " << ft0; // F_THI[0]
      fic << " " << ft1; // F_THI[1]
      fic << " " << ft2; // F_THI[2]
      fic << " " << atat; // A_THI.A_THI
      fic << " " << ftft; // F_THI.F_THI
      // u.rho*F_THI
      fic << " " << ptn; // velocity.F_THI
      // Energie cinetique en diphasiqeu (double)
      // u.rho.u (qui est aussi accessible par les .txt)
      fic << " " << uru;

      fic << finl;
      fic.close();
    }
  statistiques().end_count(source_counter_);
}

#ifdef COMPLEMENT_ANTI_DEVIATION_RESIDU
// left average - right average
static double calculer_wall_difference(const IJK_Field_double& vx)
{
  const int nj = vx.nj();
  const int ni = vx.ni();
  const Domaine_IJK& geom = vx.get_domaine();
  // Maillage uniforme, il suffit donc de diviser par le nombre total de mailles:
  // cast en double au cas ou on voudrait faire un maillage >2 milliards
  const double n_mailles_plan_xy = ((double) geom.get_nb_elem_tot(0)) * geom.get_nb_elem_tot(1);
  const int kmin = vx.get_domaine().get_offset_local(DIRECTION_K);
  const Domaine_IJK::Localisation loc = vx.get_localisation();
  const int nktot = vx.get_domaine().get_nb_items_global(loc, DIRECTION_K);

  double x = 0.;
  if (kmin == 0)
    {
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
          x += vx(i, j, 0);
    }
  if (kmin + vx.nk() == nktot)
    {
      const int k = vx.nk() - 1;
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
          x -= vx(i, j, k);
    }

  // somme sur tous les processeurs.
  x = Process::mp_sum(x);
  // Il faut diviser par n_mailles_plan_xy
  x /= n_mailles_plan_xy;
  return x;
}
#endif

// force_tot est homogene a une force volumique (comparable a la source S), comme tout ce qu'on evalue ici.
// Attention, cette methode ne fait pas que des evaluations, elle calcule aussi terme_source_correction_x_ ou _y_
void Navier_Stokes_FTD_IJK::compute_correction_for_momentum_balance(const int rk_step)
{
  // Toutes les interfaces etant fermees, on a donc la somme des forces
  // donnee par la poussee d'archimede : delta_rho * g * alpha
  const Probleme_FTD_IJK_base& pb_ijk = probleme_ijk();
  double vol_cell = 1., vol_NS = 1., vol_gaz = 0.;
  const Domaine_IJK& geom_NS = pb_ijk.get_domaine();
  for (int direction = 0; direction < 3; direction++)
    {
      vol_NS *= geom_NS.get_domain_length(direction);
      vol_cell *= geom_NS.get_constant_delta(direction);
    }

  ArrOfDouble volume;
  DoubleTab position;
  const int nbulles = interfaces_->get_nb_bulles_reelles();
  // La methode calcule a present les volumes meme pour les bulles ghost.
  // Pour les enlever, il suffit simplement de reduire la taille du tableau :
  interfaces_->calculer_volume_bulles(volume, position);
  volume.resize_array(nbulles);
  position.resize(nbulles, 3);
  for (int ib = 0; ib < nbulles; ib++)
    vol_gaz += volume[ib];

  update_rho_v();
  Vecteur3 force_tot, force_theo, rhov_moy, tauw;
  Vecteur3 acc, residu; // homogenes a d(rhou)/dt
#ifdef COMPLEMENT_ANTI_DEVIATION_RESIDU
  double moins_delta_Pwall_sur_h = 0.;
#endif
  tauw = 0.;
  // Flottaison opposee a la gravite :
  const double rho_l = milieu_ijk().get_rho_liquid(), rho_v = milieu_ijk().get_rho_vapour(), mu_l = milieu_ijk().get_mu_liquid();

  const double drho_alpha = -(rho_l - rho_v) * vol_gaz / vol_NS;
  const double un_sur_h = 2. / geom_NS.get_domain_length(DIRECTION_K);
  double fractional_dt;
  const double timestep = schema_temps_ijk().get_timestep();

  if (sub_type(Schema_RK3_IJK, schema_temps_ijk()))
    fractional_dt = compute_fractionnal_timestep_rk3(timestep /* total*/, rk_step);
  else
    // On suppose que c'est euler :
    fractional_dt = timestep;

  // Convertir en des contributions volumiques (homogene au terme source) :
  for (int direction = 0; direction < 3; direction++)
    {
      if (direction != DIRECTION_K)
        {
          // Calcul du frottement moyen sur un mur (selon X ou Y)
          tauw[direction] = calculer_tau_wall(velocity_[direction], mu_l);
#ifdef COMPLEMENT_ANTI_DEVIATION_RESIDU
        }
      else
        {
          // On n'est pas certain qu'il (le gradient normal de Uz moyenne sur un mur)
          // soit nul en instantane. Pour verifier, on le calcule:
          tauw[direction] = calculer_tau_wall(velocity_[direction], mu_l);
          // La fonction calculer_tau_wall suppose que le champ de vitesse est a dz/2 du mur.
          // Pour la vitesse uz, il est en fait situe a dz. Donc il faut le doubler:
          tauw[direction] *= 2.;
          // Il y a aussi potentiellement une partie due a la pression en instantanne:
          // (dans la direction normale aux parois seulement)
          // +/-?  (-Pw^+ + Pw^- ) / h :
          moins_delta_Pwall_sur_h = -calculer_wall_difference(pressure_) * un_sur_h;
#endif
        }
      force_tot[direction] = calculer_v_moyen(terme_source_interfaces_ns_[direction]) / vol_cell;
      if (interfaces_->is_terme_gravite_rhog())
        force_theo[direction] = 0.;
      else
        force_theo[direction] = drho_alpha * milieu_ijk().gravite().valeurs()(0, direction);

      rhov_moy[direction] = calculer_v_moyen(rho_v_[direction]);
      acc[direction] = (rhov_moy[direction] - store_rhov_moy_[direction]) / fractional_dt;
      store_rhov_moy_[direction] = rhov_moy[direction];
      // Une partie de la force en x ou en z et la force en y sont des erreurs de discretisation,
      // que l'on peut corriger globalement :
      terme_source_correction_[direction] = force_theo[direction] - force_tot[direction];
      // Le frottement moyen doit etre divise par h pour etre rendu "volumique"
      // Il est oriente vers z- donc signe "-"
      tauw[direction] *= -un_sur_h;
    }

  const int tstep = schema_temps_ijk().get_tstep();

  if (tstep == 0)
    {
      // L'acceleration n'est pas valide car  store_rhov_moy_^0 = rhov_moy^0 .. donc on ferait (0-0)/0...
      // Donc le residu est invalide.
      // Donc on ne l'ajoute pas a l'integrated_residu_, qui vaut 0 a ce moment la...
      residu = 0.;
      integrated_residu_ = 0.;
    }
  else
    {
      for (int dir = 0; dir < 3; dir++)
        {
          residu[dir] = acc[dir] - (force_tot[dir] + tauw[dir] + correction_force_[dir] * terme_source_correction_[dir]);

          if (dir == milieu_ijk().get_direction_gravite())
            residu[dir] -= terme_source_acceleration_;

          if (interfaces_->is_terme_gravite_rhog())
            residu[dir] = -(rho_l - (rho_l - rho_v) * vol_gaz / vol_NS) * milieu_ijk().gravite().valeurs()(0, dir);

          integrated_residu_[dir] += residu[dir] * fractional_dt;
        }
    }

  // Impression dans un fichier dedie :
  if (Process::je_suis_maitre())
    {
      int reset = (!pb_ijk.get_reprise()) && (tstep == 0);
      SFichier fic = Ouvrir_fichier("_bilan_qdm.out", "tstep\ttime\tFs_theo\tFtot\trhov\ttauw\tS\tacceleration\tres\tcumul_res\tCorrFs\nConvection\tDiffusion\tPression\n# Forces have 3 components.",
                                    reset, 20/*prec*/);

      fic << tstep << " " << schema_temps_ijk().get_current_time() << " " << force_theo[0] << " " << force_theo[1] << " " << force_theo[2] << " ";
      fic << force_tot[0] << " " << force_tot[1] << " ";
      fic << force_tot[2] << " " << rhov_moy[0] << " " << rhov_moy[1] << " " << rhov_moy[2] << " " << tauw[0] << " ";
      fic << tauw[1] << " " << tauw[2] << " ";
      switch(milieu_ijk().get_direction_gravite())
        {
        case 0:
          fic << terme_source_acceleration_ << " " << "0. " << "0. ";
          break;
        case 1:
          fic << "0. " << terme_source_acceleration_ << " " << "0. ";
          break;
        case 2:
          fic << "0. " << "0. " << terme_source_acceleration_ << " ";
          break;
        default:
          fic << terme_source_acceleration_ << " " << "0. " << "0. ";
        }
      fic << acc[0] << " " << acc[1] << " " << acc[2] << " " << residu[0] << " " << residu[1] << " ";
      fic << residu[2] << " " << integrated_residu_[0] << " " << integrated_residu_[1] << " " << integrated_residu_[2] << " ";
      fic << terme_source_correction_[0] << " " << terme_source_correction_[1] << " " << terme_source_correction_[2] << " ";
      // GAB, qdm
      // fic << terme_interfaces[0] << " "
      // << terme_interfaces[1] << " "
      // << terme_interfaces[2] << " ";
      fic << terme_convection_[0] << " " << terme_convection_[1] << " " << terme_convection_[2] << " ";
      fic << terme_diffusion_[0] << " " << terme_diffusion_[1] << " " << terme_diffusion_[2] << " ";
      fic << terme_pression_[0] << " " << terme_pression_[1] << " " << terme_pression_[2] << " ";
      //
#ifdef COMPLEMENT_ANTI_DEVIATION_RESIDU
      fic << moins_delta_Pwall_sur_h << " ";
#endif
      fic << finl;
      fic.close();
    }
}

// Mettre rk_step = -1 si schema temps different de rk3.
// /!\ rk_step = 0 signifie qu'on est en RK3, mais n'indique pas a quelle ss pas de temps de RK3 on se place !!!
void Navier_Stokes_FTD_IJK::calculer_dv(const double timestep, const double time, const int rk_step)
{
  // GAB : initialisation. On initialise que pour rk_step<=0, mais on pourrai initialiser a chaque ss pdt
  // Avec le post-traitement que je fais de mes termes, je dois les re-initialiser a chaque ss pdt,
  //   si je ne les re-initialise pas alors je n'aurai pas leur valeur pour chaque avancement (apres je pourrai me passer de ce niveau de detail)
  //   ==> Ce qui est certain c'est que je ne dois pas re initialiser ICI mes rho_u_... puisqu'ils sont
  //       evalues EN DEHORS des boucles de RK3.

  Probleme_FTD_IJK_base& pb_ijk = probleme_ijk();
  if (rk_step <= 0)
    {
      terme_convection_ = 0.;
      terme_pression_bis_ = 0.;
      terme_pression_ter_ = 0.;
      terme_diffusion_ = 0.;
      terme_interfaces_ = 0.;
      terme_interfaces_bf_mass_solver_ = 0.;
    }
  // GAB
  // /!\ Valable que pour un domaine uniforme, j'ai vu des choses que je ne comprends pas la ou on defini volume aussi ...
  //     on est dans variable dz et on ne prends pas en compte k dans le calcul du volume...
  double volume_cell_uniforme = 1.;
  for (int i = 0; i < 3; i++)
    volume_cell_uniforme *= pb_ijk.domaine_ijk().get_constant_delta(i);
  if (velocity_reset_)
    for (int dir = 0; dir < 3; dir++)
      velocity_[dir].data() = 0.; //Velocity reset for test

  static Stat_Counter_Id calcul_dv_counter = statistiques().new_counter(2, "maj vitesse : calcul derivee vitesse");
  statistiques().begin_count(calcul_dv_counter);
  // Calcul d_velocity = convection

  const double rho_l = milieu_ijk().get_rho_liquid();

  if (!disable_convection_qdm_)
    {
      if (velocity_convection_op_.get_convection_op_option_rank() == non_conservative_simple)
        {
          velocity_convection_op_->calculer(velocity_[0], velocity_[1], velocity_[2], velocity_[0], velocity_[1], velocity_[2], d_velocity_[0], d_velocity_[1], d_velocity_[2]);
          // Multiplication par rho (on va rediviser a la fin)
          // (a partir de rho aux elements et dv aux faces)
          if (use_inv_rho_for_mass_solver_and_calculer_rho_v_)
            {
              // rho_face = 2*(rho_gauche*rho_droite)/(rho_gauche+rho_droite)
              //          = 1./ (1/2 * (1/rho_g + 1/rho_d))
              // 1/rho_face est donc la moyenne geometrique de inv_rho.
              calculer_rho_harmonic_v(rho_field_, d_velocity_, d_velocity_);
            }
          else
            calculer_rho_v(rho_field_, d_velocity_, d_velocity_);
        }
      else if (velocity_convection_op_.get_convection_op_option_rank() == non_conservative_rhou)
        {
          update_rho_v();

          rho_v_[0].echange_espace_virtuel(2);
          rho_v_[1].echange_espace_virtuel(2);
          rho_v_[2].echange_espace_virtuel(2);

          // Non optimise car la methode calculer_avec_u_div_rhou inexistante.
          // Alors on initialise a 0 puis on fait ajouter :
          d_velocity_[0].data() = 0.;
          d_velocity_[1].data() = 0.;
          d_velocity_[2].data() = 0.;
          if (velocity_convection_op_.get_convection_op() != Nom("Centre4"))
            {
              velocity_convection_op_->ajouter_avec_u_div_rhou(rho_v_[0], rho_v_[1], rho_v_[2], // rhov_
                                                               velocity_[0], velocity_[1], velocity_[2], d_velocity_[0], d_velocity_[1], d_velocity_[2], div_rhou_);
            }
        }
      else if (velocity_convection_op_.get_convection_op_option_rank() == conservative)
        {
          update_rho_v();

          rho_v_[0].echange_espace_virtuel(2);
          rho_v_[1].echange_espace_virtuel(2);
          rho_v_[2].echange_espace_virtuel(2);

          velocity_convection_op_->calculer(rho_v_[0], rho_v_[1], rho_v_[2], velocity_[0], velocity_[1], velocity_[2], d_velocity_[0], d_velocity_[1], d_velocity_[2]);
        }
      else
        {
          Cerr << "Unknown velocity convection type! " << finl;
          Process::exit();
        }
      pb_ijk.get_post().fill_op_conv();
      // GAB, qdm

      // a priori homogene a int_{volume_cellule} (d rho v / dt) pour le moment
      // on divise le terme_convection par volume_cellule
      // Il ne faudrai pas un mass solveur sur mon terme de convection... avant de le moyenner meme ?

      for (int dir = 0; dir < 3; dir++)
        {
          // if (rk_step==-1 || rk_step==0) // euler ou premier pdt de rk3
          // Pourquoi diviser par volume_cell_uniforme ?
          terme_convection_[dir] = calculer_v_moyen(d_velocity_[dir]) / volume_cell_uniforme;

          if (test_etapes_et_bilan_)
            {
              terme_convection_mass_solver_[dir] = d_velocity_[dir];
              if (!Option_IJK::DISABLE_DIPHASIQUE)
                {
                  for (int k = 0; k < d_velocity_[dir].nk(); k++)
                    mass_solver_with_rho(terme_convection_mass_solver_[dir], rho_field_, pb_ijk.get_delta_z_local(), k);
                }
              else
                {
                  for (int k = 0; k < d_velocity_[dir].nk(); k++)
                    for (int j = 0; j < d_velocity_[dir].nj(); j++)
                      for (int i = 0; i < d_velocity_[dir].ni(); i++)
                        terme_convection_mass_solver_[dir](i, j, k) = terme_convection_mass_solver_[dir](i, j, k) / (rho_l * volume_cell_uniforme);
                }

              terme_moyen_convection_mass_solver_[dir] = calculer_v_moyen(terme_convection_mass_solver_[dir]);
            }
        }
    }
  else
    {
      d_velocity_[0].data() = 0.;
      d_velocity_[1].data() = 0.;
      d_velocity_[2].data() = 0.;
    }

  // Calcul diffusion
  if ((!diffusion_alternative_) && (!disable_diffusion_qdm_))
    {
      velocity_diffusion_op_->set_nu(molecular_mu_);
      velocity_diffusion_op_->ajouter(velocity_[0], velocity_[1], velocity_[2], d_velocity_[0], d_velocity_[1], d_velocity_[2]);
      // GAB, qdm
      // a priori homogene a int_{volume_cellule} (d rho v / dt) pour le moment
      // mais on le divise par volume_cell_uniforme donc homogene a d rho v / dt maintenant
      // en diffusion "normale", on ajoute la diffusion plus tard
      for (int dir = 0; dir < 3; dir++)
        {
          // if (rk_step==-1 || rk_step==0) // revient a rk_step>0 // euler ou premier pdt de rk3
          terme_diffusion_[dir] = calculer_v_moyen(d_velocity_[dir]) / volume_cell_uniforme - terme_convection_[dir];
          if (test_etapes_et_bilan_)
            {
              terme_diffusion_mass_solver_[dir] = d_velocity_[dir];
              // terme_diffusion_mass_solver_ contient la diffusion et la convection
              if (!Option_IJK::DISABLE_DIPHASIQUE)
                {
                  for (int k = 0; k < terme_diffusion_mass_solver_[dir].nk(); k++)
                    {
                      mass_solver_with_rho(terme_diffusion_mass_solver_[dir], rho_field_, pb_ijk.get_delta_z_local(), k);

                      // On retranche la convection
                      for (int j = 0; j < terme_diffusion_mass_solver_[dir].nj(); j++)
                        for (int i = 0; i < terme_diffusion_mass_solver_[dir].ni(); i++)
                          terme_diffusion_mass_solver_[dir](i, j, k) = terme_diffusion_mass_solver_[dir](i, j, k) - terme_convection_mass_solver_[dir](i, j, k);
                    }
                }
              else
                {
                  // Dans le cas monophasique, rho_field_ vaut partout rho_liquide
                  for (int k = 0; k < terme_diffusion_mass_solver_[dir].nk(); k++)
                    for (int j = 0; j < terme_diffusion_mass_solver_[dir].nj(); j++)
                      for (int i = 0; i < terme_diffusion_mass_solver_[dir].ni(); i++)
                        {
                          terme_diffusion_mass_solver_[dir](i, j, k) = terme_diffusion_mass_solver_[dir](i, j, k) / (rho_l * volume_cell_uniforme);
                          terme_diffusion_mass_solver_[dir](i, j, k) = terme_diffusion_mass_solver_[dir](i, j, k) - terme_convection_mass_solver_[dir](i, j, k);
                        }
                }
              // Moyenne sur le volume du domaine de simulation
              terme_moyen_diffusion_mass_solver_[dir] = calculer_v_moyen(terme_diffusion_mass_solver_[dir]);
            }
        }
    }

  // GAB, qdm
  if (test_etapes_et_bilan_)
    for (int i = 0; i < 3; i++)
      d_v_diff_et_conv_[i] = d_velocity_[i];

  // Calcul et ajout du grad(P^{n})*volume_cell (_entrelace?) :
  if (include_pressure_gradient_in_ustar_)
    {
      // pressure gradient requires the "left" value in all directions:
      pressure_.echange_espace_virtuel(1 /*, IJK_Field_double::EXCHANGE_GET_AT_LEFT_IJK*/);
#ifndef VARIABLE_DZ
      double volume = 1.;
      for (int i = 0; i < 3; i++)
        volume *= pb_ijk.domaine_ijk().get_constant_delta(i);
#else
      Cerr << "This methods does not support variable DZ yet... Contact trust support. ";
      Process::exit();
      const double volume = get_channel_control_volume(dv, k, delta_z_local_);
#endif
      add_gradient_times_constant(pressure_, -volume /*constant*/, d_velocity_[0], d_velocity_[1], d_velocity_[2]);
      Cout << "BF add_gradient_times_constant" << finl;
      // GAB, qdm
      // En utilisation normale de la pression, elle n'est pas ajoutee ici
      if (test_etapes_et_bilan_)
        {
          add_gradient_times_constant(pressure_, -volume, terme_pression_in_ustar_local_[0], terme_pression_in_ustar_local_[1], terme_pression_in_ustar_local_[2]);
          for (int dir_press = 0; dir_press < 3; dir_press++)
            terme_pression_in_ustar_[dir_press] = calculer_v_moyen(terme_pression_in_ustar_local_[dir_press]);
          Cout << "AF add_gradient_times_constant" << finl;
        }
      //
    }

  const double sch_timestep = schema_temps_ijk().get_timestep();
  const int sch_tstep = schema_temps_ijk().get_tstep();

  // Calcul du terme source aux interfaces pour l'ajouter a dv :
  // ATTENZIONE : Questa è una bugia. I termini di interfaccia vengono aggiunti direttamente a velocity_!
  if (!Option_IJK::DISABLE_DIPHASIQUE)
    {
      for (int dir = 0; dir < 3; dir++)
        {
          terme_source_interfaces_ft_[dir].data() = 0.;
          terme_repulsion_interfaces_ft_[dir].data() = 0.;
          terme_abs_repulsion_interfaces_ft_[dir].data() = 0.;
        }
      interfaces_->ajouter_terme_source_interfaces(terme_source_interfaces_ft_, terme_repulsion_interfaces_ft_, terme_abs_repulsion_interfaces_ft_);

      // Avant le solveur de masse, il faut un terme homogene a \int_vol {rho v }
      if (!disable_source_interf_)
        {
          for (int dir = 0; dir < 3; dir++)
            {
              if ((rk_step == -1) || (refuse_patch_conservation_QdM_RK3_source_interf_))
                {
                  // Avec le schema Euler, ou si on refuse le patch de conservation de la QdM en RK3 diphasique :
                  redistribute_from_splitting_ft_faces_[dir].redistribute_add(terme_source_interfaces_ft_[dir], d_velocity_[dir]);
                }
              else
                {
                  // On n'ajoute pas le terme source a d_velocity_...
                  // On le prendra en compte directement dans velocity_.
                }
            }

          for (int dir = 0; dir < 3; dir++)
            {
              redistribute_from_splitting_ft_faces_[dir].redistribute(terme_source_interfaces_ft_[dir], terme_source_interfaces_ns_[dir]);
              // GAB, qdm : les forces d'interface sont directement ajoutees a velocity... aucun effet sur d_velocity normalement
              //            Donc terme_interfaces_bf_mass_solver_bis est nul. C'EST A VERIFIER !!
              if (test_etapes_et_bilan_)
                {
                  terme_interfaces_bf_mass_solver_[dir] = calculer_v_moyen(terme_source_interfaces_ns_[dir]);
                  terme_interfaces_bf_mass_solver_bis_[dir] = calculer_v_moyen(d_velocity_[dir]) / volume_cell_uniforme - terme_convection_[dir] - terme_diffusion_[dir];
                }
            }
          pb_ijk.get_post().fill_surface_force(terme_source_interfaces_ns_);
          // Computing force_tot (Attention, il faut le faire avant d'appliquer le solver mass a terme_source_interfaces_ns_) :
          compute_correction_for_momentum_balance(rk_step);
          for (int dir = 0; dir < 3; dir++)
            {
              // On est en RK3 et on utilise le patch de conservation de la QdM (comportement par defaut du RK3)
              // Utilisation directe du terme source interf pour l'ajouter a velocity_.
              // On ne le met plus dans d_velocity_ car ce n'est pas conservatif globalement... (test quand sigm et drho)
              if (pb_ijk.get_post().get_liste_post_instantanes().contient_("REPULSION_FT") || pb_ijk.get_post().get_liste_post_instantanes().contient_("CELL_REPULSION_FT"))
                {
                  redistribute_from_splitting_ft_faces_[dir].redistribute(terme_repulsion_interfaces_ft_[dir], terme_repulsion_interfaces_ns_[dir]);
                }
              const int kmax = terme_source_interfaces_ns_[dir].nk();
              for (int k = 0; k < kmax; k++)
                {
                  // division par le produit (volume * rho_face)
                  if (use_inv_rho_for_mass_solver_and_calculer_rho_v_)
                    {
                      Cerr << "Je ne sais pas si inv_rho_field_ est a jour ici. A Verifier avant de l'activer." << finl;
                      //calculer_rho_harmonic_v(rho_field_, d_velocity_, d_velocity_);
                      mass_solver_with_inv_rho(terme_source_interfaces_ns_[dir], inv_rho_field_, pb_ijk.get_delta_z_local(), k);
                    }
                  else
                    {
                      // Division de terme_source_interfaces_ns_ par rho_field et par volume cellule
                      //cout << " code Delta : "<<terme_repulsion_interfaces_ns_[0](6,9,10);
                      //cout << " code Delta : "<<backup_terme_source_interfaces_ns_[0](6,9,10);
                      //cout << " code Delta : "<<post_.get_rho_Ssigma()[2](6,11,11);
                      //cout << endl;
                      mass_solver_with_rho(terme_source_interfaces_ns_[dir], rho_field_, pb_ijk.get_delta_z_local(), k);
                      //cout << " code Nabla : "<<terme_repulsion_interfaces_ns_[0](6,9,10);
                      //cout << " code Nabla : "<<backup_terme_source_interfaces_ns_[0](6,9,10);
                      //cout << " code Nabla : "<<post_.get_rho_Ssigma()[2](6,11,11);
                      //cout << endl;
                      //if (k==10)
                      //  {
                      //    cout << "code 2211 " << rho_field_(5,9,10) << " & " << rho_field_(6,9,10) << " & " << rho_field_(7,9,10)<<endl;
                      //    cout << "code Z,611 " << post_.get_rho_Ssigma()[2](6,11,11);
                      //    cout << ", code Z,611 " << backup_terme_source_interfaces_ns_[2](6,11,11);
                      //    cout << ", code Z,611 " << terme_source_interfaces_ns_[2](6,11,11)<<endl;
                      //    cout << "code X,6910 " << post_.get_rho_Ssigma()[0](6,9,10);
                      //    cout << ", code X,6910 " << backup_terme_source_interfaces_ns_[0](6,9,10);
                      //    cout << ", code X,6910 " << terme_source_interfaces_ns_[0](6,9,10)<<endl;
                      //  }
                      if (pb_ijk.get_post().get_liste_post_instantanes().contient_("REPULSION_FT") || pb_ijk.get_post().get_liste_post_instantanes().contient_("CELL_REPULSION_FT"))
                        {
                          // Division de terme_repulsion_interfaces_ns_ par rho_field_ et par volume cellule
                          mass_solver_with_rho(terme_repulsion_interfaces_ns_[dir], rho_field_, pb_ijk.get_delta_z_local(), k);
                          //terme_repulsion_interfaces_ns_[0](6,9,10)=50;
                        }
                      // Egalite aux dimensions :
                      // [terme_source_interfaces_ns_]=[terme_repulsion_interfaces_ns_]=[du/dt / Vcell] = m/s^2/m^3
                    }
                  if ((!refuse_patch_conservation_QdM_RK3_source_interf_) && (rk_step >= 0))
                    {
                      // puis
                      // comme euler_explicit_update mais avec un pas de temps partiel :
                      const double delta_t = compute_fractionnal_timestep_rk3(sch_timestep /* total*/, rk_step);
                      const int imax = terme_source_interfaces_ns_[dir].ni();
                      const int jmax = terme_source_interfaces_ns_[dir].nj();
                      for (int j = 0; j < jmax; j++)
                        for (int i = 0; i < imax; i++)
                          {
                            double x = terme_source_interfaces_ns_[dir](i, j, k);
                            velocity_[dir](i, j, k) += x * delta_t;
                          }
                    }
                  // On est dans une boucle sur les directions la, c ok
                  // GAB, qdm  ATTENTION on ne va ici que si on est en rk3
                  if (test_etapes_et_bilan_)
                    terme_interfaces_af_mass_solver_[dir] = calculer_v_moyen(terme_source_interfaces_ns_[dir]);
                }
            }
        }
    }

  // On laisse l'ecriture de ce fichier de sortie A L'INTERIEUR de calculer_dv car on souhaite relever
  // la valeur des differents termes A CHAQUE sous-pas de temps du schema RK3. On peut en revanche,

  fill_variable_source_and_potential_phi(time);

  // Correcteur PID

  ArrOfDouble acc_rmf;
  acc_rmf.resize_array(3);
  acc_rmf = 0.;

  if (Kp_ != 0. || Kd_ != 0. || Ki_ != 0.)
    {
      double ax_PID;
      double ay_PID;
      double az_PID;

      calculer_terme_asservissement(ax_PID, ay_PID, az_PID);

      acc_rmf[0] = ax_PID;
      acc_rmf[1] = ay_PID;
      acc_rmf[2] = az_PID;
    }

  for (int dir = 0; dir < 3; dir++)
    {
      // GAB question : pour quoi on cree dv, qu'est ce qui empeche de travailler sur d_velocity ???
      //                -> Est ce que c'est uniquement parce qu'on va faire un certain nb d'operations dessus
      //                   et que manipuler dv est plus leger que de manipuler d_velocity ?
      IJK_Field_double& dv = d_velocity_[dir];
      const int kmax = d_velocity_[dir].nk();
      const int ni = dv.ni();
      const int nj = dv.nj();
      // terme source acceleration homogene a une force volumique (gradient de pression uniforme)
      // Si la correction est activee, on oppose la force_moy
      double force_volumique = correction_force_[dir] * terme_source_correction_[dir];
      // if (dir == DIRECTION_I)
      //   {
      //     force_volumique += terme_source_acceleration_;
      //   }
      // GAB, rotation
      if (dir == milieu_ijk().get_direction_gravite())
        force_volumique += terme_source_acceleration_;

#ifndef VARIABLE_DZ
      double volume = 1.;
      for (int i = 0; i < 3; i++)
        volume *= pb_ijk.domaine_ijk().get_constant_delta(i);

      // GAB, qdm
      // dans d_velocity_moyen on a la contrib de interfaces, forces ajoutees
      terme_interfaces_conv_diff_mass_solver_[dir] = calculer_v_moyen(d_velocity_[dir]);

      Cerr << "disable_diffusion_qdm_ : " << disable_diffusion_qdm_ << finl;
      Cerr << "diffusion_alternative_ : " << diffusion_alternative_ << finl;
      Cerr << "type_velocity_diffusion_form : " << velocity_diffusion_op_.get_diffusion_op_option() << finl;
      for (int k = 0; k < kmax; k++)
        {
          // #else
          //       for (int k = 0; k < kmax; k++)
          //         {
          //           const double volume = get_channel_control_volume(dv, k, delta_z_local_);
          //         }
          const double f = force_volumique * volume;
          if ((expression_variable_source_[0] != "??") || (expression_variable_source_[1] != "??") || (expression_variable_source_[2] != "??") || (expression_potential_phi_ != "??"))
            {
              for (int j = 0; j < nj; j++)
                for (int i = 0; i < ni; i++)
                  dv(i, j, k) += facteur_variable_source_ * variable_source_[dir](i, j, k) * volume + f;
            }
          else
            {
              for (int j = 0; j < nj; j++)
                for (int i = 0; i < ni; i++)
                  dv(i, j, k) += f;
            }

          if (use_inv_rho_for_mass_solver_and_calculer_rho_v_)
            mass_solver_with_inv_rho(d_velocity_[dir], inv_rho_field_, pb_ijk.get_delta_z_local(), k);
          else
            mass_solver_with_rho(d_velocity_[dir], rho_field_, pb_ijk.get_delta_z_local(), k);

          if (Kp_ != 0. || Kd_ != 0. || Ki_ != 0.)
            {
              for (int j = 0; j < nj; j++)
                for (int i = 0; i < ni; i++)
                  d_velocity_[dir](i, j, k) += acc_rmf[dir];
            }

          // Terme source gravitaire en version simple "rho_g":
          // GAB : question a Guillaume : l 3235 -> "IJK_Field_double& dv = d_velocity_[dir];" ,  MAIS C'EST APRES
          // LA LIGNE 2865 ou d_velocity subi un calculer_rho_v.
          //       donc dv est homogene a rho*v et donc a rho*g. Or 8 lignes plus bas on a dv(i,j,k) += g; !!!
          // il faut appliquer un mass_solver_with_rho a dv pour que l'operation soit homogene d'apres gr...
          // a moins que le mass_solver soit applique dans une des methodes...
          // mass_solver_with... est applique sur d_velocity_ juste au dessus de ce long commentaire. Donc d_velocity
          // est bien homogene a g, et ainsi (jespere) dv est bien homogene a g.
          if (pb_ijk.has_interface() && interfaces_->is_terme_gravite_rhog())
            {
              const double g = milieu_ijk().gravite().valeurs()(0, dir);
              for (int j = 0; j < nj; j++)
                for (int i = 0; i < ni; i++)
                  dv(i, j, k) += g;  // c'est rho*g qu'il faut pour GR, 18.08.2021
            }

          if ((diffusion_alternative_) && (!disable_diffusion_qdm_))
            {
              velocity_diffusion_op_->set_nu(unit_);
              velocity_diffusion_op_->ajouter(velocity_[0], velocity_[1], velocity_[2], laplacien_velocity_[0], laplacien_velocity_[1], laplacien_velocity_[2]);
              for (int dir2 = 0; dir2 < 3; dir2++)
                {
                  const int kmax2 = d_velocity_[dir2].nk();
                  const int imax = d_velocity_[dir2].ni();
                  const int jmax = d_velocity_[dir2].nj();
                  for (int k2 = 0; k2 < kmax2; k2++)
                    for (int j = 0; j < jmax; j++)
                      for (int i = 0; i < imax; i++)
                        {
                          // double laplacien_u = laplacien_velocity_[dir2](i,j,k2);
                          // GAB, c'est assez dangereux de nommer nu une viscosite dynamique...
                          // const double nu = molecular_mu_(i,j,k2) ; // On stocke nu dans mu dans ce cas.
                          // GAB, qdm
                          if (test_etapes_et_bilan_)
                            {
                              terme_diffusion_local_[dir2](i, j, k2) = molecular_mu_(i, j, k2) * laplacien_velocity_[dir2](i, j, k2);
                              // Cerr << "terme diffusion local" << terme_diffusion_local[dir2](i,j,k2) << finl;
                              d_velocity_[dir2](i, j, k2) += terme_diffusion_local_[dir2](i, j, k2);
                              // d_velocity_[dir2](i,j,k2) += nu * laplacien_u;
                            }
                        }
                  // GAB, qdm
                  d_velocity_[dir2].echange_espace_virtuel(d_velocity_[dir2].ghost());
                  //
                }
            }

        }
      // GAB, TODO bilan qdm
      // d_velocity_moyen_ap_mass_solver[dir] = calculer_v_moyen(d_velocity_[dir]);
      // dans d_velocity_conv_et_diff_moy on a que la contrib de convection et de diffusion
      // d_velocity_conv_et_diff_moy[dir] = calculer_v_moyen(d_velocity_conv_et_diff[dir])

#endif
      // For the addition of the external forces to d_velocity_
      compute_add_external_forces(dir);
    } // end of loop [dir].
  ///////////////////////////////////////////////////////
  // GAB, THI
  // MODIF GAB.. /!\ tstep : le numero d'iteration; timestep : le pas de temps, le dt
  // est-on au moment ou velocity_ contient la vitesse ou est-on au moment ou velocity_ contient la qdm ?
  if (forcage_.get_type_forcage() > 0)  // && (rk_step==-1 || rk_step==0))
    {
      if (rk_step == -1)
        compute_add_THI_force_sur_d_velocity(velocity_, sch_tstep, sch_timestep, time, d_velocity_.get_domaine(), forcage_.get_facteur_forcage());  //, rk_step);
      else
        {
          const double intermediate_dt = compute_fractionnal_timestep_rk3(sch_timestep, rk_step);
          compute_add_THI_force_sur_d_velocity(velocity_, sch_tstep, intermediate_dt, time, d_velocity_.get_domaine(), forcage_.get_facteur_forcage());  //, rk_step);
        }
    }
  // verifier si mon terme de thi est bon en integrale
  ///////////////////////////////////////////////////////
  Cout << "G bilan qdm " << finl;
  if (test_etapes_et_bilan_)
    write_check_etapes_et_termes(rk_step);

  // Il est important de s'assurer a la fin que la derivee de la vitesse soit a zero sur les parois:
  if (!pb_ijk.domaine_ijk().get_periodic_flag(DIRECTION_K))
    force_zero_on_walls(d_velocity_[2]);

  statistiques().end_count(calcul_dv_counter);
}

void Navier_Stokes_FTD_IJK::compute_add_external_forces(const int dir)
{
  if (coef_immobilisation_ > 1e-16)
    {
      // maxValue n'avait pas le mp_max
      double integration_time = max_ijk(probleme_ijk().get_post().integrated_timescale());
      integration_time=std::max(1.,integration_time); // if integrated_timescale is missing, the value of integration will be -1.e30;
      //                                            The trick is to set it to 1, as it is in fact numbers ot timesteps stored (see dirty code)
      interfaces_->compute_external_forces_(force_rappel_ft_, force_rappel_, velocity_,interfaces_->I(),interfaces_->I_ft(),
                                            coef_immobilisation_, schema_temps_ijk().get_tstep(), schema_temps_ijk().get_current_time(),
                                            coef_ammortissement_, coef_rayon_force_rappel_,
                                            integration_time, coef_mean_force_, coef_force_time_n_);

      if (interfaces_->get_forcing_method())
        {
          // Si Parser
          const int kmax = d_velocity_[dir].nk();
          const int ni = d_velocity_[dir].ni();
          const int nj = d_velocity_[dir].nj();
          for (int k = 0; k < kmax; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                d_velocity_[dir](i,j,k) += force_rappel_[dir](i,j,k);
        }
      else
        {
          // Si la force de rappel est calculee sans le parser (par la color function), elle est evaluee sur le FT.
          redistribute_from_splitting_ft_faces_[dir].redistribute_add(
            force_rappel_ft_[dir], d_velocity_[dir]);
          // Je pense qu'il faut la redistribuer pour le calcul de la force de rappel dans le terme expression_derivee_acceleration_
          redistribute_from_splitting_ft_faces_[dir].redistribute(
            force_rappel_ft_[dir], force_rappel_[dir]);
        }
    }
  return;
}


void Navier_Stokes_FTD_IJK::fill_variable_source_and_potential_phi(const double time)
{
  // Ajout du terme source (force acceleration)
  // Solveur masse pour chaque composante du bilan de QdM
  const IJK_Field_vector3_double& grad_I_ns = probleme_ijk().get_post().get_grad_I_ns();
  for (int dir = 0; dir < 3; dir++)
    {
      // Si on est en presence d'une source analytique variable spatialement:
      if (expression_variable_source_[dir] != "??" && probleme_ijk().has_interface())
        set_field_data(variable_source_[dir], expression_variable_source_[dir], interfaces_->I(), grad_I_ns[dir], time);
      else if (expression_variable_source_[dir] != "??") // sans interface
        set_field_data(variable_source_[dir], expression_variable_source_[dir], grad_I_ns[dir], time);
      else if (expression_potential_phi_ != "??")
        {
          // Pour Remettre a zero la source:
          set_field_data(variable_source_[dir], Nom("0."), grad_I_ns[dir], time);
        }
    }

  if (expression_potential_phi_ != "??")
    {
      set_field_data(potential_phi_, expression_potential_phi_, time);
      potential_phi_.echange_espace_virtuel(potential_phi_.ghost());
      add_gradient_times_constant(potential_phi_, 1.,variable_source_[0],variable_source_[1],variable_source_[2]);
    }
}

void Navier_Stokes_FTD_IJK::write_check_etapes_et_termes(int rk_step)
{
  // GR : 07.01.22 : mettre une condition if (tstep % dt_post_stats_check_ == dt_post_stats_check_ - 1 || stop)
  //                 serai vraiment une bonne chose pour alleger les _check_....out... voir avec AB et GB.
  if (Process::je_suis_maitre())
    {
      if (test_etapes_et_bilan_)
        {
          int reset_test = (probleme_ijk().get_reprise()) && ( schema_temps_ijk().get_tstep() == 0 );

          // GAB, qdm RK3 : accurate_current_time_ = t0(1+ 1/4); t0(1+ 1/4 + 5/12); t0(1+ 1/4 + 5/12 + 1/3)

          double accurate_current_time = 0.0;
          if ( sub_type(Schema_Euler_explicite_IJK, schema_temps_ijk()) )
            accurate_current_time = schema_temps_ijk().get_current_time();
          else if ( sub_type(Schema_Euler_explicite_IJK, schema_temps_ijk()) )
            accurate_current_time = ref_cast(Schema_RK3_IJK, schema_temps_ijk()).get_current_time_at_rk3_step();

          SFichier fic_test=Ouvrir_fichier("_check_etapes_et_termes.out",
                                           "tstep\tcurrent_time\trk_step"
                                           "\trho_u_euler_av_prediction\trho_du_euler_ap_prediction"
                                           "\trho_u_euler_ap_projection\trho_du_euler_ap_projection"
                                           "\trho_u_euler_av_rho_mu_ind\trho_u_euler_ap_rho_mu_ind"
                                           "\tu_euler_ap_rho_mu_ind"
                                           "\ttemre_interfaces\tterme_convection\tterme_diffusion"
                                           "\tterme_pression_bis\tterme_pression_ter"
                                           "\tterme_interfaces_bf_mass_solver_bis\tterme_interfaces_bf_mass_solver\tterme_interfaces_bf_mass_solver"
                                           "\tpression_ap_proj"
                                           "\tdrho_u"
                                           "\tterme_moyen_convection_mass_solver"
                                           "\tterme_moyen_diffusion_mass_solver"
                                           "\n# Forces have 3 components.",
                                           reset_test, 20/*prec*/);


          // GAB, qdm
          Cout << "check etapes et termes" << finl;
          fic_test << schema_temps_ijk().get_tstep() << " " << accurate_current_time << " ";
          fic_test << rk_step << " ";
          // Inspection a gros grain
          fic_test << rho_u_euler_av_prediction_[0] << " "
                   << rho_u_euler_av_prediction_[1] << " "
                   << rho_u_euler_av_prediction_[2] << " ";  // colone 5
          fic_test << rho_du_euler_ap_prediction_[0] << " "
                   << rho_du_euler_ap_prediction_[1] << " "
                   << rho_du_euler_ap_prediction_[2] << " ";
          fic_test << rho_u_euler_ap_projection_[0] << " "
                   << rho_u_euler_ap_projection_[1] << " "  // colonne 10
                   << rho_u_euler_ap_projection_[2] << " ";
          fic_test << rho_du_euler_ap_projection_[0] << " "
                   << rho_du_euler_ap_projection_[1] << " "
                   << rho_du_euler_ap_projection_[2] << " ";
          fic_test << rho_u_euler_av_rho_mu_ind_[0] << " "  // colonne 15
                   << rho_u_euler_av_rho_mu_ind_[1] << " "
                   << rho_u_euler_av_rho_mu_ind_[2] << " ";
          fic_test << rho_u_euler_ap_rho_mu_ind_[0] << " "
                   << rho_u_euler_ap_rho_mu_ind_[1] << " "
                   << rho_u_euler_ap_rho_mu_ind_[2] << " ";  // colonne 20
          fic_test << u_euler_ap_rho_mu_ind_[0] << " "
                   << u_euler_ap_rho_mu_ind_[1] << " "
                   << u_euler_ap_rho_mu_ind_[2] << " ";          // Dans l'etape de prediction, inspection terme a terme
          fic_test << terme_interfaces_[0] << " "
                   << terme_interfaces_[1] << " "  // colonne 25
                   << terme_interfaces_[2] << " ";
          fic_test << terme_convection_[0] << " "
                   << terme_convection_[1] << " "
                   << terme_convection_[2] << " ";
          fic_test << terme_diffusion_[0] << " "  // colonne 30
                   << terme_diffusion_[1] << " "
                   << terme_diffusion_[2] << " ";
          // Moyenne_spatiale{ grad(p) }
          fic_test << terme_pression_bis_[0] << " "
                   << terme_pression_bis_[1] << " "
                   << terme_pression_bis_[2] << " ";  // colonne 35
          // Moyenne_spatiale{ 1/rho grad(p) }
          fic_test << terme_pression_ter_[0] << " "
                   << terme_pression_ter_[1] << " "
                   << terme_pression_ter_[2] << " ";
          // Inspection de l'effet du mass solver
          fic_test << terme_interfaces_bf_mass_solver_bis_[0] << " "
                   << terme_interfaces_bf_mass_solver_bis_[1] << " "  // colonne 40
                   << terme_interfaces_bf_mass_solver_bis_[2] << " ";
          fic_test << terme_interfaces_bf_mass_solver_[0] << " "
                   << terme_interfaces_bf_mass_solver_[1] << " "
                   << terme_interfaces_bf_mass_solver_[2] << " ";
          fic_test << terme_interfaces_af_mass_solver_[0] << " "  // colonne 45 /!\ au 17.01.22 une coquille a ete corrigee : terme_interfaces_bf_mass_solver-> terme_interfaces_af_mass_solver
                   << terme_interfaces_af_mass_solver_[1] << " "
                   << terme_interfaces_af_mass_solver_[2] << " ";
          fic_test << pression_ap_proj_ << " ";
          // On pourrai se passer de cette sortie et la creer uniquement dans python
          // ==> ( r^{n+1} - r^{n} ) * u^{n+1}
          fic_test << rho_u_euler_av_rho_mu_ind_[0]-rho_u_euler_ap_rho_mu_ind_[0] << " "
                   << rho_u_euler_av_rho_mu_ind_[1]-rho_u_euler_ap_rho_mu_ind_[1] << " "  // colonne 50
                   << rho_u_euler_av_rho_mu_ind_[2]-rho_u_euler_ap_rho_mu_ind_[2] << " ";
          fic_test << terme_moyen_convection_mass_solver_[0] << " "
                   << terme_moyen_convection_mass_solver_[1] << " "
                   << terme_moyen_convection_mass_solver_[2] << " ";
          fic_test << terme_moyen_diffusion_mass_solver_[0] << " "  // colonne 55
                   << terme_moyen_diffusion_mass_solver_[1] << " "
                   << terme_moyen_diffusion_mass_solver_[2] << " ";
          fic_test << finl;
          fic_test.close();
          Cout << "check etapes et termes fini" << finl;
        }
    }
}

// -----------------------------------------------------------------------------------
//  FORCAGE EXTERIEUR, DEFINI DANS L'ESPACE SPECTRAL
void Navier_Stokes_FTD_IJK::compute_add_THI_force(const IJK_Field_vector3_double& vitesse, const int time_iteration, const double dt, //tstep, /!\ ce dt est faux, je ne sais pas pk mais en comparant sa valeur avec celle du dt_ev, je vois que c'est faux
                                                  const double current_time, const Domaine_IJK& my_splitting)
{
  statistiques().begin_count(m2_counter_);
  if (forcage_.get_forced_advection() == -1)
    {
      ArrOfDouble mean_u_liq;
      mean_u_liq.resize_array(3);
      for (int dir = 0; dir < 3; dir++)
        mean_u_liq[dir] = calculer_moyenne_de_phase_liq(vitesse[dir]);
      forcage_.update_advection_velocity(mean_u_liq);
    }
  if (forcage_.get_forced_advection() != 0)
    {
      Cout << "forced_advection" << forcage_.get_forced_advection() << finl;
      Cout << "BF : update_advection_length" << finl;
      forcage_.update_advection_length(dt);
      Cout << "AF : update_advection_length" << finl;
    }
  forcage_.compute_THI_force(time_iteration, dt, current_time, probleme_ijk().domaine_ijk());

  statistiques().end_count(m2_counter_);

  statistiques().begin_count(m3_counter_);

  const IJK_Field_vector3_double& force = forcage_.get_force_ph2();
  for (int dir = 0; dir < 3; dir++)
    {
      // d_velocity_ est deja decoupe sur les differents procs; donc ni, nj, nk =! nb elem.
      const int kmax = d_velocity_[dir].nk();
      const int ni = d_velocity_[dir].ni();
      const int nj = d_velocity_[dir].nj();

      for (int k = 0; k < kmax + 1; k++)
        for (int j = 0; j < nj + 1; j++)
          for (int i = 0; i < ni + 1; i++)
            {
              // GAB error : No continuous phase found ...
              // appeler mass_solver_with_rho si je veux div par rho*V (tant que je garde f localise aux faces)
              // double cell_mass = rho_field_(i,j,k)*my_geom.get_constant_delta(DIRECTION_I)*my_geom.get_constant_delta(DIRECTION_J)*my_geom.get_constant_delta(DIRECTION_K);
              const double inv_cell_mass = 1.; //my_geom.get_constant_delta(DIRECTION_I)*my_geom.get_constant_delta(DIRECTION_J)*my_geom.get_constant_delta(DIRECTION_K);
              //d_velocity_[dir](i,j,k) += force[dir](i,j,k)*inv_cell_mass;
              velocity_[dir](i, j, k) += force[dir](i, j, k) * inv_cell_mass * dt;
            }
    }
  statistiques().end_count(m3_counter_);
}

void Navier_Stokes_FTD_IJK::compute_add_THI_force_sur_d_velocity(const IJK_Field_vector3_double& vitesse, const int time_iteration, const double dt, //tstep,  /!\ ce dt est faux, je ne sais pas pk mais en comparant sa valeur avec celle du dt_ev, je vois que c'est faux
                                                                 const double current_time, const Domaine_IJK& my_splitting, const int facteur)
{
  statistiques().begin_count(m2_counter_);
  if (forcage_.get_forced_advection() == -1)
    {
      /* Advection du champ de force par mean{u_l}^l */
      ArrOfDouble mean_u_liq;
      mean_u_liq.resize_array(3);
      for (int dir = 0; dir < 3; dir++)
        {
          // (22.02.22) FAUX : mean_u_liq[dir] contient alors alpha_liq*vitesse_liq[dir]
          // mean_u_liq[dir] = calculer_moyenne_de_phase_liq(vitesse_liq[dir]);
          // (22.02.22) JUSTE : mean_u liq[dir] contient bien vitesse_liq[dir]
          mean_u_liq[dir] = calculer_true_moyenne_de_phase_liq(vitesse[dir]);
        }
      forcage_.update_advection_velocity(mean_u_liq);
    }
  /* Advection du champ de force par advection_velocity_, donnee du jdd */
  if (forcage_.get_forced_advection() != 0)
    forcage_.update_advection_length(dt);

  forcage_.compute_THI_force(time_iteration, dt, current_time, probleme_ijk().domaine_ijk());
  statistiques().end_count(m2_counter_);

  statistiques().begin_count(m3_counter_);

  const IJK_Field_vector3_double& force = forcage_.get_force_ph2();

  for (int dir = 0; dir < 3; dir++)
    {
      // d_velocity_ est deja decoupe sur les differents procs; donc ni, nj, nk =! nb elem.
      const int kmax = d_velocity_[dir].nk();
      const int ni = d_velocity_[dir].ni();
      const int nj = d_velocity_[dir].nj();

      // La meme force appliquee partout
      if (facteur == 0)
        {
          for (int k = 0; k < kmax; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  // appeler mass_solver_with_rho si je veux div par rho*V (tant que je garde f localise aux faces)
                  // double cell_mass = rho_field_(i,j,k)*my_geom.get_constant_delta(DIRECTION_I)*my_geom.get_constant_delta(DIRECTION_J)*my_geom.get_constant_delta(DIRECTION_K);
                  //d_velocity_[dir](i,j,k) += force[dir](i,j,k)*inv_cell_mass;
                  d_velocity_[dir](i, j, k) += force[dir](i, j, k);
                }
        }
      // Force active uniquement dans le liquide. Nulle dans le gaz
      else if (facteur == 1)
        {
          for (int k = 0; k < kmax; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  // appeler mass_solver_with_rho si je veux div par rho*V (tant que je garde f localise aux faces)
                  // double cell_mass = rho_field_(i,j,k)*my_geom.get_constant_delta(DIRECTION_I)*my_geom.get_constant_delta(DIRECTION_J)*my_geom.get_constant_delta(DIRECTION_K);
                  // indicatrice == 0 dans le gaz et 1 dans le liquide. ATTN : ~0.5 aux interfaces !!!
                  d_velocity_[dir](i, j, k) += force[dir](i, j, k) * interfaces_->I(i, j, k);
                }
        }
      // Force ponderee par la masse volumique de la phase. Attention, il faut rester homogene a une vitesse
      else if (facteur == 2)
        {
          // S'assurer qu'on a bien calcule rho_moyen_
          for (int k = 0; k < kmax; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  // appeler mass_solver_with_rho si je veux div par rho*V (tant que je garde f localise aux faces)
                  // double cell_mass = rho_field_(i,j,k)*my_geom.get_constant_delta(DIRECTION_I)*my_geom.get_constant_delta(DIRECTION_J)*my_geom.get_constant_delta(DIRECTION_K);
                  d_velocity_[dir](i, j, k) += force[dir](i, j, k) * rho_field_(i, j, k) / rho_moyen_;
                }
        }
      d_velocity_[dir].echange_espace_virtuel(d_velocity_[dir].ghost());
    }
  statistiques().end_count(m3_counter_);
  Cout << "end of from_spect_to_phys_opti2_advection" << finl;
}
// -----------------------------------------------------------------------------------

double Navier_Stokes_FTD_IJK::calculer_moyenne_de_phase_liq(const IJK_Field_double& vx)
{
  /* Au 04.11.21 : Renvoi alpha_liq * vx_liq
   * */
  const Domaine_IJK& geom = vx.get_domaine();
  const int ni = vx.ni();
  const int nj = vx.nj();
  const int nk = vx.nk();
  double v_moy = 0.;

#ifndef VARIABLE_DZ
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        v_moy += vx(i, j, k) * interfaces_->I(i, j, k);

  // somme sur tous les processeurs.
  v_moy = Process::mp_sum(v_moy);
  // Maillage uniforme, il suffit donc de diviser par le nombre total de mailles:
  // cast en double au cas ou on voudrait faire un maillage >2 milliards
  const double n_mailles_tot = ((double) geom.get_nb_elem_tot(0)) * geom.get_nb_elem_tot(1) * geom.get_nb_elem_tot(2);
  v_moy /= n_mailles_tot;
#else
  const int offset = splitting.get_offset_local(DIRECTION_K);
  const ArrOfDouble& tab_dz=geom.get_delta(DIRECTION_K);
  for (int k = 0; k < nk; k++)
    {
      const double dz = tab_dz[k+offset];
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
          v_moy += vx(i,j,k)*dz;
    }
  // somme sur tous les processeurs.
  v_moy = Process::mp_sum(v_moy);
  // Maillage uniforme, il suffit donc de diviser par le nombre total de mailles:
  // cast en double au cas ou on voudrait faire un maillage >2 milliards
  const double n_mailles_xy = ((double) geom.get_nb_elem_tot(0)) * geom.get_nb_elem_tot(1);
  v_moy /= (n_mailles_xy * geom.get_domain_length(DIRECTION_K) );
#endif
  return v_moy;
}

void Navier_Stokes_FTD_IJK::compute_and_add_qdm_corrections()
{
  /*
   * Corrections to comply with the momentum budget
   * u_corrected = u - u_correction
   * u_corrected and u_corrrection are computed in the qdm_corrections_ object.
   * */
  double alpha_l = calculer_v_moyen(interfaces_->I());
  double rho_moyen = calculer_v_moyen(rho_field_);
  qdm_corrections_.set_rho_moyen_alpha_l(rho_moyen, alpha_l);
  qdm_corrections_.set_rho_liquide(milieu_ijk().get_rho_liquid());
  for (int dir = 0; dir < 3; dir++)
    {
      IJK_Field_double& vel = velocity_[dir];
      double rho_vel_moyen = calculer_v_moyen(rho_v_[dir]);
      qdm_corrections_.set_rho_vel_moyen(dir, rho_vel_moyen);
      Cout << "qdm_corrections_.get_need_for_vitesse_relative(" << dir << ")" << qdm_corrections_.get_need_for_vitesse_relative(dir) << finl;
      if (qdm_corrections_.get_need_for_vitesse_relative(dir))
        {
          double vel_rel = calculer_true_moyenne_de_phase_liq(vel) - calculer_true_moyenne_de_phase_vap(vel);
          qdm_corrections_.set_vitesse_relative(dir, vel_rel);
        }
      // TODO : Demander de l'aide a Guillaume :
      /* Pour moyenne glissante, je fais appel a des ArrOfDouble. En utilisation sequentielle, j'en
       * suis satisfait disons. En utilisation parallele, je ne sais pas vraiment comment sont gerees
       * mes listes. Sont-t-elles dupliquees ? Decoupees ? En tout cas la simu plante pour une liste
       * plus de 10 doubles, sur un maillage a 40^3 mailles, une seule bulle... */
      if (qdm_corrections_.get_need_to_compute_correction_value_one_direction(dir))
        qdm_corrections_.compute_correction_value_one_direction(dir);
      for (int k = 0; k < vel.nk(); ++k)
        for (int j = 0; j < vel.nj(); ++j)
          for (int i = 0; i < vel.ni(); ++i)
            {
              qdm_corrections_.compute_correct_velocity_one_direction(dir, vel(i, j, k));
              velocity_[dir](i, j, k) = qdm_corrections_.get_correct_velocitiy_one_direction(dir);
            }
    }

  Cout << "AF : compute_and_add_qdm_corrections" << finl;
}

// -----------------------------------------------------------------------------------
//  CORRECTION DE QdM
double Navier_Stokes_FTD_IJK::calculer_true_moyenne_de_phase_liq(const IJK_Field_double& vx)
{
  /* Au 04.11.21 : Renvoi vx_liq */
  double alpha_liq_vx_liq = calculer_moyenne_de_phase_liq(vx);
  double alpha_liq = calculer_v_moyen(interfaces_->I()); //en utilisant calculer_moyenne_de_phase_liq(indicatrice_ns_) on somme des chi^2 ce qui est genant aux mailles diphasiques
  return alpha_liq_vx_liq / alpha_liq;
}

double Navier_Stokes_FTD_IJK::calculer_true_moyenne_de_phase_vap(const IJK_Field_double& vx)
{
  /* Au 04.11.21 : Renvoi vx_vap */
  double alpha_vap_vx_vap = calculer_moyenne_de_phase_vap(vx);
  double alpha_vap = 1 - calculer_v_moyen(interfaces_->I()); //en utilisant 1-calculer_moyenne_de_phase_liq(indicatrice_ns_) on somme des chi^2 ce qui est genant aux mailles diphasiques
  return alpha_vap_vx_vap / alpha_vap;
}

double Navier_Stokes_FTD_IJK::calculer_moyenne_de_phase_vap(const IJK_Field_double& vx)
{
  /* Au 04.11.21 : Renvoi alpha_vap * vx_vap
   * */
  const Domaine_IJK& geom = vx.get_domaine();
  const int ni = vx.ni();
  const int nj = vx.nj();
  const int nk = vx.nk();
  double v_moy = 0.;

#ifndef VARIABLE_DZ
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        v_moy += vx(i, j, k) * (1 - interfaces_->I(i, j, k));

  // somme sur tous les processeurs.
  v_moy = Process::mp_sum(v_moy);
  // Maillage uniforme, il suffit donc de diviser par le nombre total de mailles:
  // cast en double au cas ou on voudrait faire un maillage >2 milliards
  const double n_mailles_tot = ((double) geom.get_nb_elem_tot(0)) * geom.get_nb_elem_tot(1) * geom.get_nb_elem_tot(2);
  v_moy /= n_mailles_tot;
#else
  const int offset = splitting.get_offset_local(DIRECTION_K);
  const ArrOfDouble& tab_dz=geom.get_delta(DIRECTION_K);
  for (int k = 0; k < nk; k++)
    {
      const double dz = tab_dz[k+offset];
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
          v_moy += vx(i,j,k)*dz;
    }
  // somme sur tous les processeurs.
  v_moy = Process::mp_sum(v_moy);
  // Maillage uniforme, il suffit donc de diviser par le nombre total de mailles:
  // cast en double au cas ou on voudrait faire un maillage >2 milliards
  const double n_mailles_xy = ((double) geom.get_nb_elem_tot(0)) * geom.get_nb_elem_tot(1);
  v_moy /= (n_mailles_xy * geom.get_domain_length(DIRECTION_K) );
#endif
  return v_moy;
}

// Hard coded constant pressure gradient in i direction, add contribution in m/s*volume of control volume
void Navier_Stokes_FTD_IJK::terme_source_gravite(IJK_Field_double& dv, int k_index, int dir) const
{
  const double constant = milieu_ijk().gravite().valeurs()(0, dir);
  const int imax = dv.ni();
  const int jmax = dv.nj();
  for (int j = 0; j < jmax; j++)
    for (int i = 0; i < imax; i++)
      dv(i, j, k_index) += constant;
}

void Navier_Stokes_FTD_IJK::set_time_for_corrections()
{
  qdm_corrections_.set_time(schema_temps_ijk().get_current_time(), schema_temps_ijk().get_timestep(), schema_temps_ijk().get_tstep());
}

void Navier_Stokes_FTD_IJK::compute_and_add_qdm_corrections_monophasic()
{
  /* For monophasic corrections only */
  /*
   * Corrections to comply with the momentum budget
   * u_corrected = u - u_correction
   * u_corrected and u_corrrection are computed in the qdm_corrections_ object.
   * */
  double alpha_l = 1.; // calculer_v_moyen(interfaces_.I())
  double rho_moyen = milieu_ijk().get_rho_liquid(); // calculer_v_moyen(rho_field_)
  qdm_corrections_.set_rho_moyen_alpha_l(rho_moyen,alpha_l);
  qdm_corrections_.set_rho_liquide(milieu_ijk().get_rho_liquid());
  for (int dir=0; dir<3; dir++)
    {
      IJK_Field_double& vel = velocity_[dir];
      double rho_vel_moyen = calculer_v_moyen(rho_v_[dir]);
      qdm_corrections_.set_rho_vel_moyen(dir,rho_vel_moyen);
      Cout << "qdm_corrections_.get_need_for_vitesse_relative("<<dir<<")" << qdm_corrections_.get_need_for_vitesse_relative(dir)<< finl;
      if (qdm_corrections_.get_need_for_vitesse_relative(dir))
        {
          double vel_rel = calculer_true_moyenne_de_phase_liq(vel); // - calculer_true_moyenne_de_phase_vap(vel);
          qdm_corrections_.set_vitesse_relative(dir, vel_rel);
        }
      // TODO : Demander de l'aide a Guillaume :
      /* Pour moyenne glissante, je fais appel a des ArrOfDouble. En utilisation sequentielle, j'en
       * suis satisfait disons. En utilisation parallele, je ne sais pas vraiment comment sont gerees
       * mes listes. Sont-t-elles dupliquees ? Decoupees ? En tout cas la simu plante pour une liste
       * plus de 10 doubles, sur un maillage a 40^3 mailles, une seule bulle... */
      if (qdm_corrections_.get_need_to_compute_correction_value_one_direction(dir))
        qdm_corrections_.compute_correction_value_one_direction(dir);
      for (int k=0; k<vel.nk(); k++)
        for (int j=0; j<vel.nj(); j++)
          for (int i=0; i<vel.ni(); i++)
            {
              qdm_corrections_.compute_correct_velocity_one_direction(dir, vel(i,j,k));
              velocity_[dir](i,j,k) = qdm_corrections_.get_correct_velocitiy_one_direction(dir);
            }
    }
  Cout << "AF : compute_and_add_qdm_corrections_monophasic" << finl;
}

static int decoder_numero_bulle(const int code)
{
  const int num_bulle = code >>6;
  return num_bulle;
}

void Navier_Stokes_FTD_IJK::compute_var_volume_par_bulle(ArrOfDouble& var_volume_par_bulle)
{
  if (vol_bulles_.size_array() > 0.)
    {
      ArrOfDouble volume_reel;
      DoubleTab position;
      interfaces_->calculer_volume_bulles(volume_reel, position);
      const int nb_reelles = interfaces_->get_nb_bulles_reelles();
      for (int ib = 0; ib < nb_reelles; ib++)
        var_volume_par_bulle[ib] = volume_reel[ib] - vol_bulles_[ib];
      // Pour les ghost : on retrouve leur vrai numero pour savoir quel est leur volume...
      for (int i = 0; i < interfaces_->get_nb_bulles_ghost(0 /* no print*/); i++)
        {
          const int ighost = interfaces_->ghost_compo_converter(i);
          const int ibulle_reelle = decoder_numero_bulle(-ighost);
          // Cerr << " aaaa " << i << " " << ighost << " " << ibulle_reelle << finl;
          var_volume_par_bulle[nb_reelles + i] = volume_reel[nb_reelles+ i] - vol_bulles_[ibulle_reelle];
        }
    }
  else
    {
      // Initialisation de vol_bulles_ a la valeur initiale du volume des bulles, dans le cas ou il n'est pas specifie dans le jeu de donnees.
      // Note : uniquement dans le cas de la correction_semi_locale du volume, pour ne pas impacter les cas tests
      if (correction_semi_locale_volume_bulle_)
        {
          ArrOfDouble volume_reel;
          DoubleTab position;
          interfaces_->calculer_volume_bulles(volume_reel, position);
          const int nb_reelles = interfaces_->get_nb_bulles_reelles();
          vol_bulles_.resize_array(nb_reelles);
          for (int ib = 0; ib < nb_reelles; ib++)
            {
              var_volume_par_bulle[ib] = 0.;
              vol_bulles_[ib] = volume_reel[ib];
            }
          // Pour les ghost : on retrouve leur vrai numero pour savoir quel est leur volume...
          for (int i = 0; i < interfaces_->get_nb_bulles_ghost(0 /* no print*/); i++)
            {
              const int ighost = interfaces_->ghost_compo_converter(i);
              const int ibulle_reelle = decoder_numero_bulle(-ighost);
              // Cerr << " aaaa " << i << " " << ighost << " " << ibulle_reelle << finl;
              var_volume_par_bulle[nb_reelles + i] = 0.;
              vol_bulles_[ibulle_reelle] = volume_reel[nb_reelles+ i];
            }
        }
    }
}

void Navier_Stokes_FTD_IJK::write_qdm_corrections_information()
{
  // Impression dans le fichier qdm_correction.out
  if ( sub_type(Schema_RK3_IJK, schema_temps_ijk()) ) // && (rk3_sub_step!=0)
    Cout << "in write_qdm_corrections_information, rk_step = " << ref_cast(Schema_RK3_IJK, schema_temps_ijk()).get_rk_step() << finl;

  Vecteur3 qdm_cible = qdm_corrections_.get_correction_values();
  Vecteur3 velocity_correction = qdm_corrections_.get_velocity_corrections();
  Vecteur3 rho_vel;
  for (int dir=0; dir<3; ++dir) rho_vel[dir] =  calculer_v_moyen(scalar_fields_product(probleme_ijk(), rho_field_,velocity_[dir],dir));

  if (Process::je_suis_maitre())
    {
      int reset = (!probleme_ijk().get_reprise()) && (schema_temps_ijk().get_tstep()==0);
      SFichier fic=Ouvrir_fichier("_qdm_correction.out",
                                  "1.iteration\t2.time\t3.qdm_cible[0]\t4.qdm_cible[1]\t5.qdm_cible[2]\t6.velocity_correction[0]\t7.velocity_correction[1]\t8.velocity_correction[2]\t9.qdm[0]\t10.qdm[1]\t11.qdm[2]",
                                  reset);
      // temps
      fic << schema_temps_ijk().get_tstep() << " ";
      fic << schema_temps_ijk().get_current_time() << " ";
      // CIBLE CONSTANE : qdm_cible = al.rl.u_cible
      fic << qdm_cible[0] << " ";
      fic << qdm_cible[1] << " ";
      fic << qdm_cible[2] << " ";
      // velocity_correction = (<r.u> - qdm_cible) / <r>
      fic << velocity_correction[0] << " ";
      fic << velocity_correction[1] << " ";
      fic << velocity_correction[2] << " ";
      // <r.u>
      fic << rho_vel[0] << " ";
      fic << rho_vel[1] << " ";
      fic << rho_vel[2] << " ";
      fic<<finl;
      fic.close();
      //     << finl;
    }
}

// GAB, qdm : construction du terme de pression pour le bilan de qdm. Je trouve ea plus logique de bouger cette fonction dans
//            IJK_Navier_Stokes_tool, mais etant donne qu'elle se trouve dans IJK_kernel, je ne sais pas trop si c'est propre que j'y touche
//            voir avec Guillaume ce qu'il en dit.
//            inspire de : void IJK_FT_Post::calculer_gradient_indicatrice_et_pression(const IJK_Field_double& indic)
Vecteur3 Navier_Stokes_FTD_IJK::calculer_inv_rho_grad_p_moyen(const IJK_Field_double& rho, const IJK_Field_double& pression)
{
  IJK_Field_vector3_double champ;
  allocate_velocity(champ, probleme_ijk().domaine_ijk(), 1);
  Vecteur3 resu;

  // Remise a zero :
  for (int dir = 0; dir < 3; dir++)
    {
      champ[dir].data() = 0.;
      // je ne sais pas a quoi ca sert, mais il est appele avant add_gradient_times_constant_over_rho dans IJK_NS_tool
      champ[dir].echange_espace_virtuel(1);
    }
  add_gradient_times_constant_over_rho(pression, rho, 1. /*constant*/, champ[0], champ[1], champ[2]);
  for (int dir = 0; dir < 3; dir++)
    {
      // Je ne sais pas e quoi cette ligne sert. elle est dans calculer_gradient_indicatrice_et_pression, alors je copie
      champ[dir].echange_espace_virtuel(1);
    }
  // calculer_rho_v(inv_rho, champ, inv_rho_champ);
  for (int dir = 0; dir < 3; dir++)
    resu[dir] = calculer_v_moyen(champ[dir]);

  return resu;
}

Vecteur3 Navier_Stokes_FTD_IJK::calculer_grad_p_moyen(const IJK_Field_double& pression)
{
  IJK_Field_vector3_double champ;
  allocate_velocity(champ, probleme_ijk().domaine_ijk(), 1);
  Vecteur3 resu;

  // Remise a zero :
  for (int dir = 0; dir < 3; dir++)
    {
      champ[dir].data() = 0.;
      // je ne sais pas a quoi ca sert, mais il est appele avant add_gradient_times_constant_over_rho dans IJK_NS_tool
      champ[dir].echange_espace_virtuel(1);
    }
  add_gradient_times_constant(pression, 1. /*constant*/, champ[0], champ[1], champ[2]);
  for (int dir = 0; dir < 3; dir++)
    {
      // Je ne sais pas e quoi cette ligne sert. elle est dans calculer_gradient_indicatrice_et_pression, alors je copie
      champ[dir].echange_espace_virtuel(1);
    }
  for (int dir = 0; dir < 3; dir++)
    resu[dir] = calculer_v_moyen(champ[dir]);

  return resu;
}

Vecteur3 Navier_Stokes_FTD_IJK::calculer_grad_p_over_rho_moyen(const IJK_Field_double& pression)
{
  /*
   * Calcule Moyenne_spatiale{ 1/rho * grad(p) }
   * */
  IJK_Field_vector3_double champ;
  allocate_velocity(champ, probleme_ijk().domaine_ijk(), 1);
  Vecteur3 resu;

  // Remise a zero :
  for (int dir = 0; dir < 3; dir++)
    {
      champ[dir].data() = 0.;
      // je ne sais pas a quoi ca sert, mais il est appele avant add_gradient_times_constant_over_rho dans IJK_NS_tool
      champ[dir].echange_espace_virtuel(1);
    }
  add_gradient_times_constant_over_rho(pression, rho_field_, 1. /*constant*/, champ[0], champ[1], champ[2]);
  for (int dir = 0; dir < 3; dir++)
    {
      // Je ne sais pas e quoi cette ligne sert. elle est dans calculer_gradient_indicatrice_et_pression, alors je copie
      champ[dir].echange_espace_virtuel(1);
    }
  for (int dir = 0; dir < 3; dir++)
    resu[dir] = calculer_v_moyen(champ[dir]);

  return resu;
}

void Navier_Stokes_FTD_IJK::euler_explicit_update(const IJK_Field_double& dv, IJK_Field_double& v, const int k_layer) const
{
  const double delta_t = schema_temps_ijk().get_timestep();
  const int imax = v.ni();
  const int jmax = v.nj();
  for (int j = 0; j < jmax; j++)
    for (int i = 0; i < imax; i++)
      {
        double x = dv(i, j, k_layer);
        v(i, j, k_layer) += x * delta_t;
      }
}

void Navier_Stokes_FTD_IJK::update_v_or_rhov(bool with_p)
{
  if (boundary_conditions_.get_correction_conserv_qdm() == 2)
    {
      update_rho_v();
      rho_field_.echange_espace_virtuel(rho_field_.ghost());
      update_v_ghost_from_rho_v();
    }
  else
    {
      // Protection to make sure that even without the activation of the flag check_divergence_, the EV of velocity is correctly field.
      // This protection MAY be necessary if convection uses ghost velocity (but I'm not sure it actually does)
      velocity_[0].echange_espace_virtuel(2);
      velocity_[1].echange_espace_virtuel(2);
      velocity_[2].echange_espace_virtuel(2);
    }

  if (with_p)
    pressure_.echange_espace_virtuel(1);
}

void Navier_Stokes_FTD_IJK::rk3_sub_step(const int rk_step, const double total_timestep, const double fractionnal_timestep, const double time)
{
  if (!frozen_velocity_)
    {
      update_v_or_rhov();

      // GAB TODO : voir dans euler_explicite ce qu'on a dit qu'on ferai pour voir
      // si le calculer_dv s'est bien passe
      Cout << "rk3ss: rk_step " << rk_step << finl;
      calculer_dv(total_timestep, time, rk_step);
      //
#ifdef PROJECTION_DE_LINCREMENT_DV
      // ajout du gradient de pression a dv
      if (!disable_solveur_poisson_)
        {
          if (include_pressure_gradient_in_ustar_)
            pressure_projection_with_rho(rho_field_, d_velocity_[0], d_velocity_[1],  d_velocity_[2],
                                         d_pressure_,1. ,  pressure_rhs_, check_divergence_, poisson_solver_);
          else
            pressure_projection_with_rho(rho_field_, d_velocity_[0], d_velocity_[1],  d_velocity_[2],
                                         pressure_,1. ,  pressure_rhs_, check_divergence_, poisson_solver_);
        }
#else
#endif

      // Mise a jour du champ de vitesse (etape de projection GAB, 28/06/21 : c'est la prediction, pas la projection non ?)
      for (int dir = 0; dir < 3; dir++)
        {
          const int kmax = d_velocity_[dir].nk();
          for (int k = 0; k < kmax; k++)
            runge_kutta3_update(d_velocity_[dir], RK3_F_velocity_[dir], velocity_[dir], rk_step, k, total_timestep);
        }

#ifdef PROJECTION_DE_LINCREMENT_DV
      // Mise a jour du champ de pression
      if ((!disable_solveur_poisson_) && (include_pressure_gradient_in_ustar_))
        {
          const int kmax = pressure_.nk();
          for (int k = 0; k < kmax; k++)
            runge_kutta3_update(d_pressure_ /* increment */,
                                RK3_F_pressure_ /* intermediate storage */,
                                pressure_ /* variable to update */, rk_step, k, total_timestep);
        }
#else
#endif

      // Conditions en entree
      if (vitesse_entree_ > -1e20)
        force_entry_velocity(velocity_[0], velocity_[1], velocity_[2], vitesse_entree_, vitesse_entree_dir_, vitesse_entree_compo_to_force_, stencil_vitesse_entree_);

      // Forcage de la vitesse en amont de la bulle :
      if (vitesse_upstream_ > -1e20)
        {
          if (IJK_Shear_Periodic_helpler::defilement_ == 1)
            {

              double vx;
              double vy;
              double vz;

              calculer_vitesse_gauche(velocity_[0], velocity_[1], velocity_[2], vx, vy, vz);

              force_upstream_velocity_shear_perio(velocity_[0], velocity_[1], velocity_[2], vitesse_upstream_, interfaces_, nb_diam_upstream_, boundary_conditions_, nb_diam_ortho_shear_perio_, vx, vy,
                                                  vz, epaisseur_maille_);
            }
          else
            force_upstream_velocity(velocity_[0], velocity_[1], velocity_[2], vitesse_upstream_, interfaces_, nb_diam_upstream_, upstream_dir_, milieu_ijk().get_direction_gravite(),
                                    upstream_stencil_);
        }
    } // end of if ! frozen_velocity
  //static Stat_Counter_Id projection_counter_ = statistiques().new_counter(0, "projection");
#ifdef PROJECTION_DE_LINCREMENT_DV
  if (0)
#else
  if (!disable_solveur_poisson_)
#endif
    {
      //statistiques().begin_count(projection_counter_);
      if (include_pressure_gradient_in_ustar_)
        {
          Cerr << "L'option include_pressure_gradient_in_ustar n'est pas encore implementee en RK3." << finl;
          Process::exit();
        }

      if (include_pressure_gradient_in_ustar_)
        {
          Cerr << "Methode incremental pour le grad(P)" << finl;
          Cerr << " Option codee uniquement pour le sch_euler... Tester et implementer si besoiN. " << finl;
          Process::exit();
        }
      // OPTION A SELECTIONNE DANS LE CAS DUN SHEAR PERIO POUR EVITER LES PBMS D INTERPOLATION DE RHO AU NIVEAU DU BORD PERIO_Z
      if (use_inv_rho_in_poisson_solver_)
        pressure_projection_with_inv_rho(inv_rho_field_, velocity_[0], velocity_[1], velocity_[2], pressure_, fractionnal_timestep, pressure_rhs_, poisson_solver_);
      else
        {
#ifdef PROJECTION_DE_LINCREMENT_DV
          // On l'a fait avant pour etre sur qu'elle soit bien dans la derivee stockee...
#else
          pressure_projection_with_rho(rho_field_, velocity_[0], velocity_[1], velocity_[2], pressure_, fractionnal_timestep, pressure_rhs_, poisson_solver_);

          // GAB TODO : cest a peu pres ici qu'il faudra travailler pour recuperer le
          // terme de pression
#endif
          // GAB TODO : checker si le passage de rho_n a rho_n+1 est bon
          // chercher ca pour l'etape de deplacement de rho : maj_indicatrice_rho_mu
        }

      // GAB, qdm : on recupere ici le terme grad(p),
      terme_pression_bis_ = calculer_grad_p_moyen(pressure_);
      // GAB, qdm : on recupere ici le terme de pression (1/rho * grad(p))
      terme_pression_ter_ = calculer_grad_p_over_rho_moyen(pressure_);
      pression_ap_proj_ += calculer_v_moyen(pressure_);

      //statistiques().end_count(projection_counter_);
    }

  if (Process::je_suis_maitre())
    {
      Cout << "Timings diff=" << statistiques().last_time(diffusion_counter_) << " conv=" << statistiques().last_time(convection_counter_);
      Cout << " src=" << statistiques().last_time(source_counter_) << finl;
    }
}

void Navier_Stokes_FTD_IJK::euler_time_step(ArrOfDouble& var_volume_par_bulle)
{
  if (!frozen_velocity_)
    {
      update_v_or_rhov();

      // GAB, qdm
      if (test_etapes_et_bilan_)
        {
          calculer_rho_v(rho_field_, velocity_, rho_u_euler_av_prediction_champ_);
          for (int dir = 0; dir < 3; dir++)
            rho_u_euler_av_prediction_[dir] = calculer_v_moyen(rho_u_euler_av_prediction_champ_[dir]);
        }

      // GAB, remarque : calculer dv calcule dv, MAIS NE L'APPLIQUE PAS au champ de vitesse !!!
      //                 l'increment de vitesse est ajoute au champ de vitesse avec euler_explicit_update
      calculer_dv(schema_temps_ijk().get_timestep(), schema_temps_ijk().get_current_time(), -1 /*rk_step = -1 pour sch euler... */);
      // GAB, qdm calculer_dv ne fait que l'etape de prediction)
      if (test_etapes_et_bilan_)
        {
          calculer_rho_v(rho_field_, d_velocity_, rho_du_euler_ap_prediction_champ_);
          for (int dir = 0; dir < 3; dir++)
            rho_du_euler_ap_prediction_[dir] = calculer_v_moyen(rho_du_euler_ap_prediction_champ_[dir]);
        }
#ifdef PROJECTION_DE_LINCREMENT_DV
      // ajout du gradient de pression a dv
      if (!disable_solveur_poisson_)
        {
          if (!include_pressure_gradient_in_ustar_)
            {
              pressure_projection_with_rho(rho_field_, d_velocity_[0], d_velocity_[1],  d_velocity_[2],
                                           pressure_, 1.,  pressure_rhs_, check_divergence_, poisson_solver_);
              // GAB --> C'est plutot ici que l(on ajoute le terme_pression !!)
            }
          else
            {
              pressure_projection_with_rho(rho_field_, d_velocity_[0], d_velocity_[1],  d_velocity_[2],
                                           d_pressure_, 1.,  pressure_rhs_, check_divergence_, poisson_solver_);

              // Then update the pressure field :
              const int kmax = pressure_.nk();
              for (int k = 0; k < kmax; k++)
                euler_explicit_update(d_pressure_, pressure_, k);
            }
        }
#endif
      // Mise a jour du champ de vitesse (etape de projection et de prediction)
      for (int dir = 0; dir < 3; dir++)
        {
          const int kmax = d_velocity_[dir].nk();
          for (int k = 0; k < kmax; k++)
            {
              // GAB, question : d_velocity est issu de calculer_dv. Il manque pas l'appliquation de l'operateur de divergence avant d'appliquer d_velocity e velocity ?
              //                 il manque au moins la multiplication par les surfaces je pense -> NON, lis bien les etapes de convection et de diffusion.
              euler_explicit_update(d_velocity_[dir], velocity_[dir], k);
            }
        }

      // GAB, qdm : cree le rho_n * v_n+1
      if (test_etapes_et_bilan_)
        {
          calculer_rho_v(rho_field_, d_velocity_, rho_du_euler_ap_projection_champ_);
          for (int dir = 0; dir < 3; dir++)
            rho_du_euler_ap_projection_[dir] = calculer_v_moyen(rho_du_euler_ap_projection_champ_[dir]);
        }

      if (test_etapes_et_bilan_)
        {
          calculer_rho_v(rho_field_, velocity_, rho_u_euler_ap_projection_champ_);
          for (int dir = 0; dir < 3; dir++)
            rho_u_euler_ap_projection_[dir] = calculer_v_moyen(rho_u_euler_ap_projection_champ_[dir]);
        }

      // Conditions en entree
      if (vitesse_entree_ > -1e20)
        force_entry_velocity(velocity_[0], velocity_[1], velocity_[2], vitesse_entree_, vitesse_entree_dir_, vitesse_entree_compo_to_force_, stencil_vitesse_entree_);

      // Forcage de la vitesse en amont de la bulle :
      if (vitesse_upstream_ > -1e20)
        {
          if (!upstream_velocity_measured_)
            {
              if (expression_vitesse_upstream_ != "??")
                {
                  std::string expr(expression_vitesse_upstream_);
                  Parser parser;
                  parser.setString(expr);
                  parser.setNbVar((int) 1);
                  parser.addVar("t");
                  parser.parseString();
                  parser.setVar((int) 0, schema_temps_ijk().get_current_time() - schema_temps_ijk().get_modified_time_ini());
                  vitesse_upstream_ = parser.eval();
                }
            }
          else
            {
              int dir = 0;
              if (upstream_dir_ == -1)
                {
                  dir = milieu_ijk().get_direction_gravite();
                  if (dir == -1)
                    dir = 0;
                }
              const DoubleTab& rising_vector = interfaces_->get_ijk_compo_connex().get_rising_vectors();
              const double velocity_magnitude = interfaces_->get_ijk_compo_connex().get_rising_velocities()[0];
              const Vecteur3& velocity_vector = interfaces_->get_ijk_compo_connex().get_rising_velocity_overall();
              velocity_bubble_new_ = velocity_vector[dir]; //  * rising_vector[dir];
              if (schema_temps_ijk().get_tstep() == 0)
                {
                  if (velocity_bubble_old_ < -1e20)
                    velocity_bubble_old_ = 0.;
                  else
                    velocity_bubble_new_ = velocity_bubble_old_;
                  if (vitesse_upstream_reprise_ < -1e20)
                    vitesse_upstream_ = -velocity_bubble_scope_;
                  else
                    vitesse_upstream_ = vitesse_upstream_reprise_;
                }
              const double delta_velocity = velocity_bubble_scope_ + velocity_bubble_new_;
              const double ddelta_velocity = (velocity_bubble_new_ - velocity_bubble_old_) / schema_temps_ijk().get_timestep();
              if (schema_temps_ijk().get_tstep() % 100)
                velocity_bubble_integral_err_ = 0.;

              velocity_bubble_integral_err_ += delta_velocity * schema_temps_ijk().get_timestep();
              vitesse_upstream_ -= delta_velocity * upstream_velocity_bubble_factor_;
              vitesse_upstream_ -= ddelta_velocity * upstream_velocity_bubble_factor_deriv_;
              vitesse_upstream_ -= velocity_bubble_integral_err_ * upstream_velocity_bubble_factor_integral_;
              Cerr << "Velocity bubble (old): " << velocity_bubble_old_ << finl;
              velocity_bubble_old_ = velocity_bubble_new_;
              Cerr << "Velocity upstream: " << vitesse_upstream_ << finl;
              Cerr << "Velocity bubble (new): " << velocity_bubble_new_ << finl;
              Cerr << "Velocity magnitude: " << velocity_magnitude << finl;
              Cerr << "Velocity dir upstream: " << rising_vector(0, dir) << finl;
            }
          vitesse_upstream_reprise_ = vitesse_upstream_;
          Cerr << "Force upstream velocity" << finl;

          if (IJK_Shear_Periodic_helpler::defilement_ == 1)
            {
              double vx;
              double vy;
              double vz;

              calculer_vitesse_gauche(velocity_[0], velocity_[1], velocity_[2], vx, vy, vz);

              force_upstream_velocity_shear_perio(velocity_[0], velocity_[1], velocity_[2], vitesse_upstream_, interfaces_.valeur(), nb_diam_upstream_, boundary_conditions_,
                                                  nb_diam_ortho_shear_perio_, vx, vy, vz, epaisseur_maille_);
            }
          else
            force_upstream_velocity(velocity_[0], velocity_[1], velocity_[2], vitesse_upstream_, interfaces_.valeur(), nb_diam_upstream_, upstream_dir_, milieu_ijk().get_direction_gravite(),
                                    upstream_stencil_);

        }
    } // end of if ! frozen_velocity
  // static Stat_Counter_Id projection_counter_ = statistiques().new_counter(0, "projection");
#ifdef PROJECTION_DE_LINCREMENT_DV
  if (0)
#else
  if (!disable_solveur_poisson_)
#endif
    {
      //  statistiques().begin_count(projection_counter_);
      if (include_pressure_gradient_in_ustar_)
        {
          Cerr << "Methode incrementale pour le grad(P)" << finl;
          Cerr << "Initialisation du d_pressure_ conservee depuis le pas de temps precedent... " << finl;
          Cerr << "Ce n'est probablement pas optimal. QQ idees dans les sources si divergence . " << finl;
          // Que vaut d_pressure ?
          // Important car c'est l'initialisation du solveur...

          // 1. raz :
          // d_pressure_.data() = 0.; // raz...
          // 2. dp = - timestep_ * (potentiel_elem - delta_rho * phi) * u . grad(I)
          //                       (sigma_ * courbure)                  on a ustar a dispo, par u^n.
          // 3. dp = - timestep_ * u . grad(P^n)
          //                       ici, on a ustar dispo, plus u^n.
          //                       Si on veut tester avec u^n (c mieux je pense), il faut initialiser dp
          //                       juste avant l'euler_explicit_update
          if (use_inv_rho_in_poisson_solver_)
            {
              pressure_projection_with_inv_rho(inv_rho_field_, velocity_[0], velocity_[1], velocity_[2], d_pressure_, schema_temps_ijk().get_timestep(), pressure_rhs_, poisson_solver_);
            }
          else
            {
#ifdef PROJECTION_DE_LINCREMENT_DV
              // On l'a fait avant pour etre sur qu'elle soit bien dans la derivee stockee...
#else
              pressure_projection_with_rho(rho_field_, velocity_[0], velocity_[1], velocity_[2], d_pressure_, schema_temps_ijk().get_timestep(), pressure_rhs_, poisson_solver_);
#endif
            }

          // Mise a jour de la pression :
          for (int dir = 0; dir < 3; dir++)
            {
              const int kmax = pressure_.nk();
              for (int k = 0; k < kmax; k++)
                euler_explicit_update(d_pressure_, pressure_, k);
            }

          Cerr << " Un exit pour voir avec gdb... " << finl;
          Cerr << " Si ca fonctionne, faire le meme en RK3... " << finl;
          // Process::exit();
        }
      else
        {
          if (use_inv_rho_in_poisson_solver_)
            pressure_projection_with_inv_rho(inv_rho_field_, velocity_[0], velocity_[1], velocity_[2], pressure_, schema_temps_ijk().get_timestep(), pressure_rhs_, poisson_solver_);
          else
            {
#ifdef PROJECTION_DE_LINCREMENT_DV
#else
              pressure_projection_with_rho(rho_field_, velocity_[0], velocity_[1], velocity_[2], pressure_, schema_temps_ijk().get_timestep(), pressure_rhs_, poisson_solver_);
#endif
            }
        }
      if (test_etapes_et_bilan_)
        {
          // GAB, qdm : recuperons le temre de pression (1/rho * grad(p)) si on fait le bilan en u (ca a du sens meme?)
          //                                                     grap(p) si on fait le bilan de qdm
          // terme_pression_bis = calculer_inv_rho_grad_p_moyen(rho_field_, pressure_);
          terme_pression_bis_ = calculer_grad_p_moyen(pressure_);
          // GAB, qdm : recuperons le terme de pression (1/rho * grad(p))
          terme_pression_ter_ = calculer_grad_p_over_rho_moyen(pressure_);
          pression_ap_proj_ = calculer_v_moyen(pressure_);
        }
    }

  Cerr << "Copy pressure on extended field for probes" << finl;
  copy_field_values(pressure_ghost_cells_, pressure_);

  Cout << "Timings diff=" << statistiques().last_time(diffusion_counter_) << " conv=" << statistiques().last_time(convection_counter_);
  Cout << " src=" << statistiques().last_time(source_counter_) << finl;
}

void Navier_Stokes_FTD_IJK::corriger_qdm()
{
  // CORRECTION DE QUANTITE DE MOUVEMENT : Correction de QdM : permet de controler la QdM globale, dans chaque direction
  if (!(qdm_corrections_.is_type_none()))
    {
      set_time_for_corrections();
      if (Option_IJK::DISABLE_DIPHASIQUE)
        compute_and_add_qdm_corrections_monophasic();
      else
        compute_and_add_qdm_corrections();
    }
  else
    {
      Cout << "qdm_corrections_.is_type_none() : " << qdm_corrections_.is_type_none() << finl;
      Cout << "terme_source_acceleration_" << terme_source_acceleration_ << finl;
    }

  if (qdm_corrections_.write_me())
    write_qdm_corrections_information();
}

// TODO FIXME DANS DOMAINE
void Navier_Stokes_FTD_IJK::build_redistribute_extended_splitting_ft()
{
  const Domaine_IJK& dom_ijk = probleme_ijk().domaine_ijk();
  const Domaine_IJK& dom_ft = probleme_ijk().get_domaine_ft();

  for (int dir = 0; dir < 3; dir++)
    {
      VECT(IntTab) map(3);
      Domaine_IJK::Localisation loc = (dir==0) ? Domaine_IJK::FACES_I : (dir==1) ? Domaine_IJK::FACES_J : Domaine_IJK::FACES_K;
      const int n_ext = dom_ijk.ft_extension();

      for (int dir2 = 0; dir2 < 3; dir2++)
        {
          const int n = dom_ijk.get_nb_items_global(loc, dir2);
          if (n_ext == 0 || !dom_ijk.get_periodic_flag(dir2))
            {
              map[dir2].resize(1,3);
              map[dir2](0,0) = 0;// source index
              map[dir2](0,1) = 0;// dest index
              map[dir2](0,2) = n;// size
            }
          else
            {
              map[dir2].resize(3,3);
              // copy NS field to central domaine of extended field
              map[dir2](0,0) = 0;
              map[dir2](0,1) = n_ext;
              map[dir2](0,2) = n;
              // copy right part of NS field to left part of extended field
              map[dir2](1,0) = n - n_ext;
              map[dir2](1,1) = 0;
              map[dir2](1,2) = n_ext;
              // copy left part of NS field to right of extended field
              map[dir2](2,0) = 0;
              map[dir2](2,1) = n + n_ext;
              map[dir2](2,2) = n_ext;
            }
        }
      redistribute_to_splitting_ft_faces_[dir].initialize(dom_ijk, dom_ft, loc, map);

      for (int dir2 = 0; dir2 < 3; dir2++)
        {
          const int n = dom_ijk.get_nb_items_global(loc, dir2);
          if (n_ext == 0 || !dom_ijk.get_periodic_flag(dir2))
            {
              map[dir2].resize(1,3);
              map[dir2](0,0) = 0;
              map[dir2](0,1) = 0;
              map[dir2](0,2) = n;
            }
          else
            {
              map[dir2].resize(1,3);
              // When copying back from extended splitting, ignore extended data, take only central part
              map[dir2](0,0) = n_ext;  // source index
              map[dir2](0,1) = 0; // dest index
              map[dir2](0,2) = n; // size
            }
        }
      redistribute_from_splitting_ft_faces_[dir].initialize(dom_ft, dom_ijk, loc, map);
    }

  // Pour les elements:
  {
    VECT(IntTab) map(3);
    Domaine_IJK::Localisation loc = Domaine_IJK::ELEM;
    const int n_ext = dom_ijk.ft_extension();

    for (int dir2 = 0; dir2 < 3; dir2++)
      {
        const int n = dom_ijk.get_nb_items_global(loc, dir2);
        if (n_ext == 0 || !dom_ijk.get_periodic_flag(dir2))
          {
            map[dir2].resize(1,3);
            map[dir2](0,0) = 0;// source index
            map[dir2](0,1) = 0;// dest index
            map[dir2](0,2) = n;// size
          }
        else
          {
            map[dir2].resize(3,3);
            // copy NS field to central domaine of extended field
            map[dir2](0,0) = 0;
            map[dir2](0,1) = n_ext;
            map[dir2](0,2) = n;
            // copy right part of NS field to left part of extended field
            map[dir2](1,0) = n - n_ext;
            map[dir2](1,1) = 0;
            map[dir2](1,2) = n_ext;
            // copy left part of NS field to right of extended field
            map[dir2](2,0) = 0;
            map[dir2](2,1) = n + n_ext;
            map[dir2](2,2) = n_ext;
          }
      }
    redistribute_to_splitting_ft_elem_.initialize(dom_ijk, dom_ft, loc, map);


    for (int dir2 = 0; dir2 < 3; dir2++)
      {
        const int n = dom_ijk.get_nb_items_global(loc, dir2);
        if (n_ext == 0 || !dom_ijk.get_periodic_flag(dir2))
          {
            map[dir2].resize(1,3);
            map[dir2](0,0) = 0;
            map[dir2](0,1) = 0;
            map[dir2](0,2) = n;
          }
        else
          {
            map[dir2].resize(1,3);
            // When copying back from extended splitting, ignore extended data, take only central part
            map[dir2](0,0) = n_ext;  // source index
            map[dir2](0,1) = 0; // dest index
            map[dir2](0,2) = n; // size
          }
      }
    redistribute_from_splitting_ft_elem_.initialize(dom_ft, dom_ijk, loc, map);

    for (int dir2 = 0; dir2 < 3; dir2++)
      {
        const int ghost_a_redistribute = 2 ;
        const int n = dom_ijk.get_nb_items_global(loc, dir2);
        if(dir2==2)
          {
            // on ne redistribue les ghost que sur z pour le shear perio
            map[dir2].resize(1,3);
            // envoyer les rangees -2, -1, 0, 1 dans 0, 1, 2, 3
            map[dir2](0,0) = n_ext-ghost_a_redistribute;  // source index
            map[dir2](0,1) = 0; // dest index
            map[dir2](0,2) = ghost_a_redistribute*2; // size
          }
        else
          {
            map[dir2].resize(1,3);
            map[dir2](0,0) = n_ext;  // source index
            map[dir2](0,1) = 0; // dest index
            map[dir2](0,2) = n; // size
          }

      }
    redistribute_from_splitting_ft_elem_ghostz_min_.initialize(dom_ft,dom_ijk, loc, map);

    for (int dir2 = 0; dir2 < 3; dir2++)
      {
        const int ghost_a_redistribute = 2 ;
        const int n = dom_ijk.get_nb_items_global(loc, dir2);
        const int n_ft = dom_ft.get_nb_items_global(loc, dir2);
        if(dir2==2)
          {
            // on ne redistribue les ghost que sur z pour le shear perio
            map[dir2].resize(1,3);
            // envoyer les rangees n-2, n-1, n, n+1 dans 4, 5, 6, 7
            map[dir2](0,0) = n_ft - n_ext - ghost_a_redistribute;  // source index
            map[dir2](0,1) = n - 2*ghost_a_redistribute; // dest index
            map[dir2](0,2) = ghost_a_redistribute*2; // size
          }
        else
          {
            map[dir2].resize(1,3);
            map[dir2](0,0) = n_ext;  // source index
            map[dir2](0,1) = 0; // dest index
            map[dir2](0,2) = n; // size
          }

      }
    redistribute_from_splitting_ft_elem_ghostz_max_.initialize(dom_ft, dom_ijk, loc, map);
  }
}

// Calcule vitesse_ft (etendue) a partir du champ de vitesse.
void Navier_Stokes_FTD_IJK::calculer_vitesse_ft()
{
  Probleme_FTD_IJK_base& pb_ijk = probleme_ijk();
  const Domaine_IJK& dom_ijk = pb_ijk.domaine_ijk();

  for (int dir = 0; dir < 3; dir++)
    redistribute_to_splitting_ft_faces_[dir].redistribute(velocity_[dir], velocity_ft_[dir]);

  redistribute_to_splitting_ft_faces_[2].redistribute(velocity_[2], velocity_ft_[2]);
  redistribute_to_splitting_ft_faces_[1].redistribute(velocity_[1], velocity_ft_[1]);
  redistribute_to_splitting_ft_faces_[0].redistribute(velocity_[0], velocity_ft_[0]);

  if (IJK_Shear_Periodic_helpler::defilement_ == 1)
    {
      // after redistribute, velocity in ft domain must be shifted by the shear
      velocity_ft_[0].redistribute_with_shear_domain_ft(velocity_[0], boundary_conditions_.get_dU_perio(boundary_conditions_.get_resolution_u_prime_()), dom_ijk.ft_extension());
      velocity_ft_[1].redistribute_with_shear_domain_ft(velocity_[1], 0., dom_ijk.ft_extension());
      velocity_ft_[2].redistribute_with_shear_domain_ft(velocity_[2], 0., dom_ijk.ft_extension());
    }

  for (int dir = 0; dir < 3; dir++)
    velocity_ft_[dir].echange_espace_virtuel(velocity_ft_[dir].ghost());
}

// Nouvelle version ou le transport se fait avec les ghost...
void Navier_Stokes_FTD_IJK::deplacer_interfaces_rk3(const double timestep, const int rk_step, ArrOfDouble& var_volume_par_bulle)
{
  Probleme_FTD_IJK_base& pb_ijk = probleme_ijk();

  //  Calculer vitesse_ft (etendue) a partir du champ de vitesse.
  static Stat_Counter_Id deplacement_interf_counter_ = statistiques().new_counter(1, "Deplacement de l'interface");
  statistiques().begin_count(deplacement_interf_counter_);

  calculer_vitesse_ft();

  // On conserve les duplicatas que l'on transporte comme le reste.
  // Normalement, transporter_maillage gere aussi les duplicatas...
  if (correction_semi_locale_volume_bulle_)
    {
      interfaces_->calculer_vecteurs_de_deplacement_rigide(vitesses_translation_bulles_, mean_bubble_rotation_vector_, centre_gravite_bulles_);

      interfaces_->transporter_maillage_deformation(correction_semi_locale_volume_bulle_, vitesses_translation_bulles_, mean_bubble_rotation_vector_, centre_gravite_bulles_,
                                                    timestep/* total meme si RK3*/, var_volume_par_bulle, rk_step);

      pb_ijk.update_pre_remeshing_indicator_field();

      interfaces_->transporter_maillage_remaillage(correction_semi_locale_volume_bulle_, vitesses_translation_bulles_, mean_bubble_rotation_vector_, centre_gravite_bulles_, timestep,
                                                   var_volume_par_bulle, rk_step, schema_temps_ijk().get_current_time());

      pb_ijk.update_post_remeshing_indicator_field();

      interfaces_->transporter_maillage_rigide(timestep, vitesses_translation_bulles_, mean_bubble_rotation_vector_, centre_gravite_bulles_, rk_step);
    }
  else
    {
      interfaces_->calculer_vecteurs_de_deplacement_rigide(vitesses_translation_bulles_, mean_bubble_rotation_vector_, centre_gravite_bulles_);

      interfaces_->transporter_maillage_deformation(correction_semi_locale_volume_bulle_, vitesses_translation_bulles_, mean_bubble_rotation_vector_, centre_gravite_bulles_,
                                                    timestep/* total meme si RK3*/, var_volume_par_bulle, rk_step);

      interfaces_->transporter_maillage_remaillage(correction_semi_locale_volume_bulle_, vitesses_translation_bulles_, mean_bubble_rotation_vector_, centre_gravite_bulles_, timestep,
                                                   var_volume_par_bulle, rk_step, schema_temps_ijk().get_current_time());
    }

  statistiques().end_count(deplacement_interf_counter_);
  // On verra a la fin du pas de temps si certaines bulles reeles sont trop proche du bord
  // du domaine etendu. Pour l'instant, dans les sous dt, on ne les transferts pas.

  // On a conserve les duplicatas donc pas besoin de les re-creer...
  // On calcule l'indicatrice du prochain pas de temps (qui correspond aux interfaces qu'on vient de deplacer.
  pb_ijk.update_indicator_field();

  update_indicatrice_variables_monofluides();
}


void Navier_Stokes_FTD_IJK::update_indicatrice_variables_monofluides()
{
  if (IJK_Shear_Periodic_helpler::defilement_ == 1)
    {
      redistribute_from_splitting_ft_elem_ghostz_min_.redistribute(interfaces_->I_ft(), I_ns_);
      redistribute_from_splitting_ft_elem_ghostz_max_.redistribute(interfaces_->I_ft(), I_ns_);
      rho_field_.get_shear_BC_helpler().set_indicatrice_ghost_zmin_(I_ns_, 0);
      rho_field_.get_shear_BC_helpler().set_indicatrice_ghost_zmax_(I_ns_, rho_field_.nk() - 4);
      molecular_mu_.get_shear_BC_helpler().set_indicatrice_ghost_zmin_(I_ns_, 0);
      molecular_mu_.get_shear_BC_helpler().set_indicatrice_ghost_zmax_(I_ns_, molecular_mu_.nk() - 4);
      if (use_inv_rho_)
        {
          inv_rho_field_.get_shear_BC_helpler().set_indicatrice_ghost_zmin_(I_ns_, 0);
          inv_rho_field_.get_shear_BC_helpler().set_indicatrice_ghost_zmax_(I_ns_, inv_rho_field_.nk() - 4);
        }
      if (boundary_conditions_.get_correction_interp_monofluide() == 1)
        {
          interfaces_->calculer_kappa_ft(kappa_ft_);
          redistribute_from_splitting_ft_elem_ghostz_min_.redistribute(kappa_ft_, kappa_ns_);
          redistribute_from_splitting_ft_elem_ghostz_max_.redistribute(kappa_ft_, kappa_ns_);
          pressure_.get_shear_BC_helpler().set_I_sig_kappa_zmin_(I_ns_, kappa_ns_, milieu_ijk().sigma(), 0);
          pressure_.get_shear_BC_helpler().set_I_sig_kappa_zmax_(I_ns_, kappa_ns_, milieu_ijk().sigma(), pressure_.nk() - 4);
        }
    }
}

void Navier_Stokes_FTD_IJK::deplacer_interfaces(const double timestep, const int rk_step, ArrOfDouble& var_volume_par_bulle, const int first_step_interface_smoothing)
{
  /*
   * TODO: Advect with zero velocity at the beggining of the simulation to use the remeshing algo
   */
  Probleme_FTD_IJK_base& pb_ijk = probleme_ijk();

  if (correction_semi_locale_volume_bulle_)
    {
      interfaces_->calculer_vecteurs_de_deplacement_rigide(vitesses_translation_bulles_, mean_bubble_rotation_vector_, centre_gravite_bulles_, first_step_interface_smoothing);

      interfaces_->transporter_maillage_deformation(correction_semi_locale_volume_bulle_, vitesses_translation_bulles_, mean_bubble_rotation_vector_, centre_gravite_bulles_,
                                                    timestep/* total meme si RK3*/, var_volume_par_bulle, rk_step, first_step_interface_smoothing);

      pb_ijk.update_pre_remeshing_indicator_field();

      interfaces_->transporter_maillage_remaillage(correction_semi_locale_volume_bulle_, vitesses_translation_bulles_, mean_bubble_rotation_vector_, centre_gravite_bulles_, timestep,
                                                   var_volume_par_bulle, rk_step, schema_temps_ijk().get_current_time());

      pb_ijk.update_post_remeshing_indicator_field();

      interfaces_->transporter_maillage_rigide(timestep, vitesses_translation_bulles_, mean_bubble_rotation_vector_, centre_gravite_bulles_, rk_step, first_step_interface_smoothing);
    }
  else
    {
      // Note : Pour ne pas alterer les cas tests, on supprime
      // les duplicatas avant le transport dans ce cas. Si on ne supprime pas
      // les duplicatas, cela fonctionne aussi, mais le resultat des cas tests est un peu different.
      // Dans ce cadre, pour que l'indicatrice intermediaire puisse etre calculee,
      // il faut dupliquer les bulles aux frontieres periodiques ; puis supprimer
      // ces bulles pour realiser le remailage de l'interface ; puis les dupliquer a nouveau pour calculer l'indicatrice finale.
      interfaces_->supprimer_duplicata_bulles();

      interfaces_->calculer_vecteurs_de_deplacement_rigide(vitesses_translation_bulles_, mean_bubble_rotation_vector_, centre_gravite_bulles_, first_step_interface_smoothing);

      interfaces_->transporter_maillage_deformation(correction_semi_locale_volume_bulle_, vitesses_translation_bulles_, mean_bubble_rotation_vector_, centre_gravite_bulles_,
                                                    timestep/* total meme si RK3*/, var_volume_par_bulle, rk_step, first_step_interface_smoothing);

      interfaces_->transporter_maillage_remaillage(correction_semi_locale_volume_bulle_, vitesses_translation_bulles_, mean_bubble_rotation_vector_, centre_gravite_bulles_, timestep,
                                                   var_volume_par_bulle, rk_step, schema_temps_ijk().get_current_time());

      interfaces_->transferer_bulle_perio();
      interfaces_->creer_duplicata_bulles();
    }
}

void Navier_Stokes_FTD_IJK::sauvegarder_equation(const Nom& lataname, SFichier& fichier) const
{
  fichier << " tinit " << schema_temps_ijk().get_current_time() << "\n"
          << " terme_acceleration_init " << terme_source_acceleration_ << "\n"
          // GAB : qdm_source. Les valeurs des attributs utiles pour le calcul de source_qdm_gr sont
          //       ecrits dans la reprise. Ils sont ecrits avec des mots-clefs qui n'ont pas vocation a
          //       etre dans un jdd. Ces mots doivent se trouver uniquement dans des fichiers sauv.
          << " reprise_vap_velocity_tmoy " << vap_velocity_tmoy_ << "\n"
          << " reprise_liq_velocity_tmoy " << liq_velocity_tmoy_ << "\n"
          << " fichier_reprise_vitesse " << basename(lataname) << "\n";
  fichier << " timestep_reprise_vitesse 1" << "\n";
  if (probleme_ijk().has_interface())
    fichier     << " interfaces " << interfaces_.valeur()  << "\n";
  fichier << " forcage " << forcage_ << "\n"
          << " corrections_qdm " << qdm_corrections_;

  fichier << "velocity_bubble_old " << velocity_bubble_old_ << "\n";
  fichier << "vitesse_upstream_reprise " << vitesse_upstream_reprise_ << "\n";
}
