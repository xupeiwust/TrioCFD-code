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
// File      : IJK_Thermal_cut_cell.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Thermal_cut_cell.h>
#include <IJK_FT_cut_cell.h>
#include <IJK_FT.h>
#include <DebogIJK.h>
#include <IJK_Navier_Stokes_tools.h>
#include <Cut_cell_convection_auxiliaire.h>
#include <Cut_cell_diffusion_auxiliaire.h>

Implemente_instanciable_sans_constructeur( IJK_Thermal_cut_cell, "IJK_Thermal_cut_cell", IJK_Thermal_base ) ;

IJK_Thermal_cut_cell::IJK_Thermal_cut_cell() :
  cut_field_temperature_(temperature_),
  cut_field_RK3_F_temperature_diff_(RK3_F_temperature_diff_),
  cut_field_RK3_F_temperature_conv_(RK3_F_temperature_),
  cut_field_div_coeff_grad_T_volume_(div_coeff_grad_T_volume_),
  cut_field_div_coeff_grad_T_volume_temp_(div_coeff_grad_T_volume_temp_),
  cut_field_d_temperature_(d_temperature_),
  cut_field_temperature_post_dying_(temperature_post_dying_),
  cut_field_temperature_post_regular_(temperature_post_regular_),
  cut_field_temperature_post_small_(temperature_post_small_),
  cut_field_temperature_post_diff_regular_(temperature_post_diff_regular_)
{
  single_phase_=0;
  type_temperature_convection_form_ = 0;  // Default value: 0 : non conservative
  conserv_energy_global_=0; // Note : doit etre zero sinon la rustine est appliquee
  E0_=0;
  computed_centred_bubble_start_ = 1.;
  single_centred_bubble_radius_ini_ = 1.e-3;
  allow_temperature_correction_for_visu_=0;
  override_vapour_mixed_values_ = 0;

  error_temperature_ana_total_ = 0.;
  error_temperature_ana_squared_total_ = 0.;
  error_temperature_ana_rel_total_ = 0.;

  source_terms_type_=2;
  source_terms_type_dict_ = Motcles(7);
  {
    source_terms_type_dict_[0] = "linear_diffusion";
    source_terms_type_dict_[1] = "spherical_diffusion";
    source_terms_type_dict_[2] = "spherical_diffusion_approx";
    source_terms_type_dict_[3] = "tangential_conv_2D";
    source_terms_type_dict_[4] = "tangential_conv_3D";
    source_terms_type_dict_[5] = "tangential_conv_2D_tangential_diffusion_3D";
    source_terms_type_dict_[6] = "tangential_conv_3D_tangentual_diffusion_3D";
  }
  source_terms_correction_=0;
  spherical_diffusion_ = 1;

  postraiter_champs_intermediaires_ = 0;

  activate_diffusion_interface_ = 0;
  etalement_diffusion_ = ETALEMENT_DIFFUSION::AUCUN_ETALEMENT;
  methode_temperature_remplissage_ = METHODE_TEMPERATURE_REMPLISSAGE::NON_INITIALISE;
  methode_flux_interface_ = METHODE_FLUX_INTERFACE::NON_INITIALISE;
  scaled_distance_flux_interface_ = 1.52*sqrt(3);

  correction_petites_cellules_diffusion_ = 0;
  diffusion_petites_cellules_ = CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE;
  convection_petites_cellules_ = CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE;

  cut_cell_conv_scheme_ = CUT_CELL_CONV_SCHEME::QUICK_OU_CENTRE2_STENCIL;

  interfacial_temperature_.set_smart_resize(1);
  interfacial_phin_ai_.set_smart_resize(1);
}

Sortie& IJK_Thermal_cut_cell::printOn( Sortie& os ) const
{
  Nom front_space = "    ";
  Nom end_space = " ";
  Nom escape = "\n";

  IJK_Thermal_base::printOn( os );
  os<< "  {\n";

  os<< "    type_T_source " << type_T_source_ << "\n";

  if (lambda_variable_)
    os<< "    lambda_variable \n";
  if (conserv_energy_global_)
    os<< "    conserv_energy_global \n";

  os<< "  \n}";
  if (override_vapour_mixed_values_)
    os << front_space << "override_vapour_mixed_values" << escape;
  if (allow_temperature_correction_for_visu_)
    os << front_space << "allow_temperature_correction_for_visu" << escape;
  if (source_terms_type_!= -1)
    os << front_space << "source_terms_type" << end_space << source_terms_type_dict_[source_terms_type_] << escape;
  os << front_space << "delta_T_subcooled_overheated" << end_space << delta_T_subcooled_overheated_ << escape;
  return os;
}

Entree& IJK_Thermal_cut_cell::readOn( Entree& is )
{
  IJK_Thermal_base::readOn( is );
  if (ghost_fluid_)
    override_vapour_mixed_values_ = 1;
  return is;
}

void IJK_Thermal_cut_cell::set_param( Param& param )
{
  IJK_Thermal_base::set_param(param);
  param.ajouter("type_temperature_convection_form", &type_temperature_convection_form_);
  param.dictionnaire("non conservative",0);
  param.dictionnaire("conservative",1);
  param.ajouter("conserv_energy_global", &conserv_energy_global_);
  param.ajouter_flag("override_vapour_mixed_values", &override_vapour_mixed_values_);
  param.ajouter_flag("allow_temperature_correction_for_visu", &allow_temperature_correction_for_visu_);
  param.ajouter("source_terms_type", &source_terms_type_);
  param.dictionnaire("linear_diffusion", 0);
  param.dictionnaire("spherical_diffusion",1);
  param.dictionnaire("spherical_diffusion_approx",2);
  param.dictionnaire("tangential_conv_2D", 3);
  param.dictionnaire("tangential_conv_3D", 4);
  param.dictionnaire("tangential_conv_2D_tangential_diffusion_3D", 5);
  param.dictionnaire("tangential_conv_3D_tangentual_diffusion_3D", 6);
  param.ajouter_flag("source_terms_correction", &source_terms_correction_);

  param.ajouter("delta_T_subcooled_overheated", &delta_T_subcooled_overheated_);

  param.ajouter_flag("postraiter_champs_intermediaires", &postraiter_champs_intermediaires_);

  param.ajouter("methode_temperature_remplissage", (int*)&methode_temperature_remplissage_);
  param.dictionnaire("non_initialise",(int)METHODE_TEMPERATURE_REMPLISSAGE::NON_INITIALISE);
  param.dictionnaire("ponderation_voisin", (int)METHODE_TEMPERATURE_REMPLISSAGE::PONDERATION_VOISIN);
  param.dictionnaire("semi_lagrangien", (int)METHODE_TEMPERATURE_REMPLISSAGE::SEMI_LAGRANGIEN);

  param.ajouter_flag("activate_diffusion_interface", &activate_diffusion_interface_);
  param.ajouter("etalement_diffusion", (int*)&etalement_diffusion_);
  param.dictionnaire("aucun_etalement",(int)ETALEMENT_DIFFUSION::AUCUN_ETALEMENT);
  param.dictionnaire("flux_interface_uniquement", (int)ETALEMENT_DIFFUSION::FLUX_INTERFACE_UNIQUEMENT);
  param.dictionnaire("divergence_uniquement", (int)ETALEMENT_DIFFUSION::DIVERGENCE_UNIQUEMENT);
  param.dictionnaire("flux_interface_et_divergence", (int)ETALEMENT_DIFFUSION::FLUX_INTERFACE_ET_DIVERGENCE);
  param.ajouter("methode_flux_interface", (int*)&methode_flux_interface_);
  param.dictionnaire("non_initialise", (int)METHODE_FLUX_INTERFACE::NON_INITIALISE);
  param.dictionnaire("interp_pure", (int)METHODE_FLUX_INTERFACE::INTERP_PURE);
  param.dictionnaire("interp_cut_cell", (int)METHODE_FLUX_INTERFACE::INTERP_CUT_CELL);
  param.dictionnaire("local_cellule", (int)METHODE_FLUX_INTERFACE::LOCAL_CELLULE);
  param.ajouter("scaled_distance_flux_interface", &scaled_distance_flux_interface_);

  param.ajouter_flag("correction_petites_cellules_diffusion", &correction_petites_cellules_diffusion_);
  param.ajouter("diffusion_petites_cellules", (int*)&diffusion_petites_cellules_);
  param.dictionnaire("direction_privilegiee", (int)CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE);
  param.dictionnaire("direction_privilegiee_2", (int)CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE_2);
  param.dictionnaire("correction_symetrique", (int)CORRECTION_PETITES_CELLULES::CORRECTION_SYMETRIQUE);
  param.dictionnaire("direction_privilegiee_avec_limitation", (int)CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE_AVEC_LIMITATION);
  param.dictionnaire("direction_privilegiee_avec_limitation_2", (int)CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE_AVEC_LIMITATION_2);
  param.dictionnaire("correction_symetrique_avec_limitation", (int)CORRECTION_PETITES_CELLULES::CORRECTION_SYMETRIQUE_AVEC_LIMITATION);
  param.ajouter("convection_petites_cellules", (int*)&convection_petites_cellules_);
  param.dictionnaire("direction_privilegiee", (int)CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE);
  param.dictionnaire("direction_privilegiee_2", (int)CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE_2);
  param.dictionnaire("correction_symetrique", (int)CORRECTION_PETITES_CELLULES::CORRECTION_SYMETRIQUE);
  param.dictionnaire("direction_privilegiee_avec_limitation", (int)CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE_AVEC_LIMITATION);
  param.dictionnaire("direction_privilegiee_avec_limitation_2", (int)CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE_AVEC_LIMITATION_2);
  param.dictionnaire("correction_symetrique_avec_limitation", (int)CORRECTION_PETITES_CELLULES::CORRECTION_SYMETRIQUE_AVEC_LIMITATION);

  param.ajouter("cut_cell_conv_scheme", (int*)&cut_cell_conv_scheme_);
  param.dictionnaire("quick_ou_centre2_stencil", (int)CUT_CELL_CONV_SCHEME::QUICK_OU_CENTRE2_STENCIL);
  param.dictionnaire("quick_ou_centre2_perpendicular_distance", (int)CUT_CELL_CONV_SCHEME::QUICK_OU_CENTRE2_PERPENDICULAR_DISTANCE);
  param.dictionnaire("quick_ou_lineaire2_stencil", (int)CUT_CELL_CONV_SCHEME::QUICK_OU_LINEAIRE2_STENCIL);
  param.dictionnaire("quick_ou_lineaire2_perpendicular_distance", (int)CUT_CELL_CONV_SCHEME::QUICK_OU_LINEAIRE2_PERPENDICULAR_DISTANCE);
  param.dictionnaire("quick_ou_amont_stencil", (int)CUT_CELL_CONV_SCHEME::QUICK_OU_AMONT_STENCIL);
  param.dictionnaire("quick_ou_amont_perpendicular_distance", (int)CUT_CELL_CONV_SCHEME::QUICK_OU_AMONT_PERPENDICULAR_DISTANCE);
  param.dictionnaire("centre2", (int)CUT_CELL_CONV_SCHEME::CENTRE2);
  param.dictionnaire("lineaire2", (int)CUT_CELL_CONV_SCHEME::LINEAIRE2);
  param.dictionnaire("amont", (int)CUT_CELL_CONV_SCHEME::AMONT);
  param.dictionnaire("interp_facette_ou_quick_ou_amont_stencil", (int)CUT_CELL_CONV_SCHEME::INTERP_FACETTE_OU_QUICK_OU_AMONT_STENCIL);
  param.dictionnaire("interp_point_ou_quick_ou_amont_stencil", (int)CUT_CELL_CONV_SCHEME::INTERP_POINT_OU_QUICK_OU_AMONT_STENCIL);
}

int IJK_Thermal_cut_cell::initialize(const IJK_Splitting& splitting, const int idx)
{
  Cout << que_suis_je() << "::initialize()" << finl;

  // Cut-cell variables
  ref_ijk_ft_cut_cell_ = ref_cast(IJK_FT_cut_cell, ref_ijk_ft_.valeur());

  int nalloc = 0;

  cut_field_temperature_.associer_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());

  cut_field_RK3_F_temperature_diff_.associer_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());
  cut_field_RK3_F_temperature_conv_.associer_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());

  temperature_post_small_.allocate(splitting, IJK_Splitting::ELEM, 2);
  cut_field_temperature_post_small_.associer_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());

  if (postraiter_champs_intermediaires_)
    {
      temperature_post_dying_.allocate(splitting, IJK_Splitting::ELEM, 2);
      temperature_post_regular_.allocate(splitting, IJK_Splitting::ELEM, 2);
      temperature_post_diff_regular_.allocate(splitting, IJK_Splitting::ELEM, 2);
      cut_field_temperature_post_dying_.associer_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());
      cut_field_temperature_post_regular_.associer_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());
      cut_field_temperature_post_diff_regular_.associer_persistant(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());
    }

  cut_cell_flux_diffusion_.associer_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());
  cut_cell_flux_convection_.associer_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());

  for (int d = 0; d < 3; d++)
    {
      temperature_face_[0][d].allocate(splitting, IJK_Splitting::ELEM, 2);
      temperature_face_ft_[0][d].allocate(ref_ijk_ft_cut_cell_->get_splitting_ft(), IJK_Splitting::ELEM, 4);
      temperature_face_[1][d].allocate(splitting, IJK_Splitting::ELEM, 2);
      temperature_face_ft_[1][d].allocate(ref_ijk_ft_cut_cell_->get_splitting_ft(), IJK_Splitting::ELEM, 4);
    }
  nalloc += 6;

  temperature_ft_.allocate(ref_ijk_ft_cut_cell_->get_splitting_ft(), IJK_Splitting::ELEM, 2);
  nalloc += 1;

  convective_correction_.temperature_remplissage_.associer_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());

  diffusive_correction_.flux_interface_ft_.allocate(ref_ijk_ft_cut_cell_->get_splitting_ft(), IJK_Splitting::ELEM, 2);
  diffusive_correction_.flux_interface_ns_.allocate(ref_ijk_ft_cut_cell_->get_splitting_ns(), IJK_Splitting::ELEM, 2);
  diffusive_correction_.flux_interface_ft_.data() = 0.;
  diffusive_correction_.flux_interface_ns_.data() = 0.;
  diffusive_correction_.flux_interface_.associer_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());

  cut_field_div_coeff_grad_T_volume_.associer_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());
  cut_field_d_temperature_.associer_ephemere(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());

  if (etalement_diffusion_ == ETALEMENT_DIFFUSION::DIVERGENCE_UNIQUEMENT || etalement_diffusion_ == ETALEMENT_DIFFUSION::FLUX_INTERFACE_ET_DIVERGENCE)
    {
      div_coeff_grad_T_volume_temp_.allocate(splitting, IJK_Splitting::ELEM, 2);
      cut_field_div_coeff_grad_T_volume_temp_.associer_paresseux(*ref_ijk_ft_cut_cell_->get_cut_cell_disc());
      nalloc += 1;
    }

  nalloc = IJK_Thermal_base::initialize(splitting, idx);
  lambda_.allocate(splitting, IJK_Splitting::ELEM, 1);
  nalloc += 2;

  RK3_F_temperature_diff_.allocate(splitting, IJK_Splitting::ELEM, 0);
  nalloc += 1;

  temperature_diffusion_op_.typer("OpDiffIJKScalar_cut_cell_double");
  temperature_diffusion_op_.initialize(splitting);
  temperature_diffusion_op_.set_uniform_lambda_liquid(lambda_liquid_);
  temperature_diffusion_op_.set_uniform_lambda_vapour(lambda_vapour_);
  temperature_diffusion_op_.set_lambda(lambda_);

  temperature_convection_op_.typer("OpConvQuickIJKScalar_cut_cell_double");
  temperature_convection_op_.initialize(splitting);

  // Already allocated if rho_cp_post
  if (!rho_cp_post_)
    {
      rho_cp_.allocate(splitting, IJK_Splitting::ELEM, 2);
      nalloc += 1;
    }

  // Compute initial energy :
  if (conserv_energy_global_)
    {
      Cerr << "conserv_energy_global_ is used." << finl;
      Process::exit();
    }
  if (type_temperature_convection_form_==1)
    {
      rho_cp_T_.allocate(splitting, IJK_Splitting::ELEM, 2);
      div_rho_cp_T_.allocate(splitting, IJK_Splitting::ELEM, 0);
      nalloc += 2;
    }
  Cout << "End of " << que_suis_je() << "::initialize()" << finl;
  return nalloc;
}

void IJK_Thermal_cut_cell::update_thermal_properties()
{
  //IJK_Thermal_base::update_thermal_properties();
  const IJK_Field_double& indic = ref_ijk_ft_cut_cell_->itfce().I();
  // Nombre de mailles du domaine NS :
  const int nx = indic.ni();
  const int ny = indic.nj();
  const int nz = indic.nk();
  const double rho_l = ref_ijk_ft_cut_cell_->get_rho_l();
  const double rho_v = ref_ijk_ft_cut_cell_->get_rho_v();
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        {
          double chi_l = indic(i,j,k);
          rho_cp_(i,j,k) = rho_l * cp_liquid_ * chi_l + rho_v * cp_vapour_ * (1.- chi_l);
          lambda_(i,j,k) = lambda_liquid_  * chi_l + (1.- chi_l) * lambda_vapour_ ;
        }

  lambda_.echange_espace_virtuel(lambda_.ghost());
  rho_cp_.echange_espace_virtuel(rho_cp_.ghost());
}

Nom IJK_Thermal_cut_cell::compute_quasi_static_spherical_diffusion_expression(const double& time_scope,
                                                                              const int index_bubble,
                                                                              const int index_bubble_real)
{

  if (computed_centred_bubble_start_)
    {
      const DoubleTab& bubbles_centres = ref_ijk_ft_cut_cell_->itfce().get_ijk_compo_connex().get_bubbles_barycentre();
      double x,y,z;
      x = bubbles_centres(index_bubble,0);
      y = bubbles_centres(index_bubble,1);
      z = bubbles_centres(index_bubble,2);
      if (index_bubble != index_bubble_real)
        {
          const double x_real = bubbles_centres(index_bubble_real,0);
          const double y_real = bubbles_centres(index_bubble_real,1);
          const double z_real = bubbles_centres(index_bubble_real,2);
          const double lx = temperature_.get_splitting().get_grid_geometry().get_domain_length(0);
          const double ly = temperature_.get_splitting().get_grid_geometry().get_domain_length(0);
          const double lz = temperature_.get_splitting().get_grid_geometry().get_domain_length(0);
          x = (abs(x - x_real)< (lx / 2.)) ? x_real : ((x_real < (lx / 2.)) ? x_real + lx : x_real - lx);
          y = (abs(y - y_real)< (ly / 2.)) ? y_real : ((y_real < (ly / 2.)) ? y_real + ly : y_real - ly);
          z = (abs(z - z_real)< (lz / 2.)) ? z_real : ((z_real < (lz / 2.)) ? z_real + lz : z_real - lz);
        }
      return generate_expression_temperature_ini(time_scope, x, y, z);
    }
  return generate_expression_temperature_ini(time_scope, 0., 0., 0.);
}

void IJK_Thermal_cut_cell::compute_temperature_init()
{
  if (cut_field_temperature_.get_cut_cell_disc().get_n_tot() == 0)
    {
      IJK_Thermal_base::compute_temperature_init();
    }
  else
    {
      Cout << "Temperature initialization from expression \nTini = " << expression_T_init_ << finl;
      cut_field_temperature_.set_field_data(expression_T_init_, ref_ijk_ft_->itfce().I(), 0.);
    }
}

void IJK_Thermal_cut_cell::recompute_temperature_init()
{
  if (cut_field_temperature_.get_cut_cell_disc().get_n_tot() == 0)
    {
      IJK_Thermal_base::recompute_temperature_init();
    }
  else
    {
      Cout << "Temperature initialization from expression \nTini = " << expression_T_init_ << finl;
      cut_field_temperature_.set_field_data(expression_T_init_, ref_ijk_ft_->itfce().In(), 0.);
    }

  // Calcul du flux a l'interface pour le premier pas de temps.
  // Cela n'est pas fait a l'initialisation car il faut que les structures diphasiques soient creees.
  // A deplacer dans un endroit plus approprie.
  diffusive_correction_.calculer_flux_interface(methode_flux_interface_, scaled_distance_flux_interface_, lambda_liquid_, lambda_vapour_, interfacial_temperature_, interfacial_phin_ai_, cut_field_temperature_, ref_ijk_ft_cut_cell_, temperature_, temperature_ft_);
}

struct CutCell_GlobalInfo
{
  double overall;
  double overall_l;
  double overall_v;
  double pure;
  double diph_l;
  double diph_v;
  double diph_small;
  double diph_regular;
  double diph_nascent;
  double diph_dying;
};

void IJK_Thermal_cut_cell::euler_time_step(const double timestep)
{
  // Calcul du flux a l'interface au premier pas de temps.
  // Pour les autres pas de temps, la valeur a la fin du pas de temps precedent est utilisee.
  // Cela n'est pas fait a l'initialisation car il faut que les structures diphasiques soient creees.

  if (debug_)
    Cerr << "Thermal Euler time-step" << finl;

  update_thermal_properties();

  const Cut_field_vector& cut_field_velocity = ref_ijk_ft_cut_cell_->get_cut_field_velocity();
  calculer_dT_cut_cell(cut_field_velocity); // Note : ne fait rien

  const double current_time = ref_ijk_ft_cut_cell_->get_current_time();

  if (!conv_temperature_negligible_)
    {
      if (methode_temperature_remplissage_ == METHODE_TEMPERATURE_REMPLISSAGE::PONDERATION_VOISIN)
        {
          convective_correction_.calcule_temperature_remplissage_ponderation_voisin(cut_field_velocity, cut_field_temperature_);
        }
      else if (methode_temperature_remplissage_ == METHODE_TEMPERATURE_REMPLISSAGE::SEMI_LAGRANGIEN)
        {
          convective_correction_.calcule_temperature_remplissage_semi_lagrangien(timestep, lambda_liquid_, lambda_vapour_, diffusive_correction_.flux_interface_ns_, cut_field_velocity, cut_field_temperature_);
        }
      else
        {
          Cerr << "Methode non reconnue pour le calcul de la temperature de remplissage." << finl;
          Process::exit();
        }
    }

  if (cut_cell_conv_scheme_ == CUT_CELL_CONV_SCHEME::INTERP_POINT_OU_QUICK_OU_AMONT_STENCIL)
    {
      convective_correction_.calcule_temperature_face_depuis_centre(lambda_liquid_, lambda_vapour_, diffusive_correction_.flux_interface_ns_, cut_field_temperature_, temperature_face_);
    }
  else if (cut_cell_conv_scheme_ == CUT_CELL_CONV_SCHEME::INTERP_FACETTE_OU_QUICK_OU_AMONT_STENCIL)
    {
      diffusive_correction_.calculer_flux_interface(methode_flux_interface_, scaled_distance_flux_interface_, lambda_liquid_, lambda_vapour_, interfacial_temperature_, interfacial_phin_ai_, cut_field_temperature_, ref_ijk_ft_cut_cell_, temperature_, temperature_ft_);

      convective_correction_.calcule_temperature_face_depuis_facette(lambda_liquid_, lambda_vapour_, interfacial_temperature_, interfacial_phin_ai_, cut_field_temperature_, ref_ijk_ft_cut_cell_, temperature_face_ft_, temperature_face_);
    }
  else
    {
      // Do nothing
    }

  const CutCell_GlobalInfo ene_ini = compute_global_energy_cut_cell(cut_field_temperature_, 0);

  if (!conv_temperature_negligible_)
    {
      convective_correction_.add_convection_dying_cells(convection_petites_cellules_, cut_field_velocity, cut_field_temperature_);
      cut_field_temperature_.echange_espace_virtuel(temperature_.ghost());
    }

  if (postraiter_champs_intermediaires_)
    {
      cut_field_temperature_post_dying_.copy_from(cut_field_temperature_);
    }
  const CutCell_GlobalInfo ene_postconv_dying = compute_global_energy_cut_cell(cut_field_temperature_, 0);

  compute_temperature_convection_cut_cell(cut_field_velocity);
  const CutCell_GlobalInfo d_ene_Conv = compute_global_energy_cut_cell(cut_field_d_temperature_, 1);
  ref_ijk_ft_cut_cell_->euler_explicit_update_cut_cell_transport(cut_field_d_temperature_, cut_field_temperature_);
  cut_field_temperature_.echange_espace_virtuel(temperature_.ghost());

  if (postraiter_champs_intermediaires_)
    {
      cut_field_temperature_post_regular_.copy_from(cut_field_temperature_);
    }
  const CutCell_GlobalInfo ene_postconv_regular = compute_global_energy_cut_cell(cut_field_temperature_, 1);

  if (!conv_temperature_negligible_)
    {
      convective_correction_.add_convection_small_nascent_cells(convection_petites_cellules_, cut_field_velocity, cut_field_temperature_);
      cut_field_temperature_.echange_espace_virtuel(temperature_.ghost());
    }

  cut_field_temperature_post_small_.copy_from(cut_field_temperature_);
  const CutCell_GlobalInfo ene_postconv_small_nascent = compute_global_energy_cut_cell(cut_field_temperature_, 1);

  if (activate_diffusion_interface_)
    {
      diffusive_correction_.calculer_flux_interface(methode_flux_interface_, scaled_distance_flux_interface_, lambda_liquid_, lambda_vapour_, interfacial_temperature_, interfacial_phin_ai_, cut_field_temperature_, ref_ijk_ft_cut_cell_, temperature_, temperature_ft_);
    }
  else
    {
      Cerr << "Warning: The diffusion at the interface is not activated." << finl;
    }
  add_temperature_diffusion();
  const CutCell_GlobalInfo d_ene_Diffu = compute_global_energy_cut_cell(cut_field_d_temperature_, 1);
  ref_ijk_ft_cut_cell_->euler_explicit_update_cut_cell_notransport(cut_field_d_temperature_, cut_field_temperature_);
  cut_field_temperature_.echange_espace_virtuel(temperature_.ghost());

  if (postraiter_champs_intermediaires_)
    {
      cut_field_temperature_post_diff_regular_.copy_from(cut_field_temperature_);
    }
  const CutCell_GlobalInfo ene_postdiff_regular = compute_global_energy_cut_cell(cut_field_temperature_, 1);


  if (!diff_temperature_negligible_ && correction_petites_cellules_diffusion_)
    {
      diffusive_correction_.add_diffusion_small_cells(diffusion_petites_cellules_, cut_field_temperature_post_small_, cut_field_temperature_);
      cut_field_temperature_.echange_espace_virtuel(temperature_.ghost());
    }

  const CutCell_GlobalInfo ene_postdiff_small = compute_global_energy_cut_cell(cut_field_temperature_, 1);

  // Recalcul du flux a l'interface pour l'estimation de la temperature de remplissage au prochain pas de temps
  diffusive_correction_.calculer_flux_interface(methode_flux_interface_, scaled_distance_flux_interface_, lambda_liquid_, lambda_vapour_, interfacial_temperature_, interfacial_phin_ai_, cut_field_temperature_, ref_ijk_ft_cut_cell_, temperature_, temperature_ft_);

  Cerr.precision(16);
  Cerr << "dT_budget, overall time: " << current_time << " dT_conv: " << d_ene_Conv.overall << " dT_diff: " << d_ene_Diffu.overall << finl;
  Cerr << "dT_budget, overall_l time: " << current_time << " dT_conv: " << d_ene_Conv.overall_l << " dT_diff: " << d_ene_Diffu.overall_l << finl;
  Cerr << "dT_budget, overall_v time: " << current_time << " dT_conv: " << d_ene_Conv.overall_v << " dT_diff: " << d_ene_Diffu.overall_v << finl;
  Cerr << "dT_budget, pure.......... time: " << current_time << " dT_conv: " << d_ene_Conv.pure         << " dT_diff: " << d_ene_Diffu.pure         << finl;
  Cerr << "dT_budget, diph_l........ time: " << current_time << " dT_conv: " << d_ene_Conv.diph_l       << " dT_diff: " << d_ene_Diffu.diph_l       << finl;
  Cerr << "dT_budget, diph_v........ time: " << current_time << " dT_conv: " << d_ene_Conv.diph_v       << " dT_diff: " << d_ene_Diffu.diph_v       << finl;
  Cerr << "dT_budget, diph_small.... time: " << current_time << " dT_conv: " << d_ene_Conv.diph_small   << " dT_diff: " << d_ene_Diffu.diph_small   << finl;
  Cerr << "dT_budget, diph_regular.. time: " << current_time << " dT_conv: " << d_ene_Conv.diph_regular << " dT_diff: " << d_ene_Diffu.diph_regular << finl;
  Cerr << "dT_budget, diph_nascent.. time: " << current_time << " dT_conv: " << d_ene_Conv.diph_nascent << " dT_diff: " << d_ene_Diffu.diph_nascent << finl;
  Cerr << "dT_budget, diph_dying.... time: " << current_time << " dT_conv: " << d_ene_Conv.diph_dying   << " dT_diff: " << d_ene_Diffu.diph_dying   << finl;

  Cerr << "T_Budget, overall time: " << current_time << " initial: " << ene_ini.overall << " post_conv_dying: " << ene_postconv_dying.overall << " post_conv_regular: " << ene_postconv_regular.overall << " post_conv_small_nascent: " << ene_postconv_small_nascent.overall << " post_diff_regular: " << ene_postdiff_regular.overall << " post_diff_small: " << ene_postdiff_small.overall << finl;
  Cerr << "T_Budget, overall_l time: " << current_time << " initial: " << ene_ini.overall_l << " post_conv_dying: " << ene_postconv_dying.overall_l << " post_conv_regular: " << ene_postconv_regular.overall_l << " post_conv_small_nascent: " << ene_postconv_small_nascent.overall_l << " post_diff_regular: " << ene_postdiff_regular.overall_l << " post_diff_small: " << ene_postdiff_small.overall_l << finl;
  Cerr << "T_Budget, overall_v time: " << current_time << " initial: " << ene_ini.overall_v << " post_conv_dying: " << ene_postconv_dying.overall_v << " post_conv_regular: " << ene_postconv_regular.overall_v << " post_conv_small_nascent: " << ene_postconv_small_nascent.overall_v << " post_diff_regular: " << ene_postdiff_regular.overall_v << " post_diff_small: " << ene_postdiff_small.overall_v << finl;
  Cerr << "T_Budget, pure.......... time: " << current_time << " initial: " << ene_ini.pure         << " post_conv_dying: " << ene_postconv_dying.pure         << " post_conv_regular: " << ene_postconv_regular.pure         << " post_conv_small_nascent: " << ene_postconv_small_nascent.pure         << " post_diff_regular: " << ene_postdiff_regular.pure << " post_diff_small: " << ene_postdiff_small.pure         << finl;
  Cerr << "T_Budget, diph_l........ time: " << current_time << " initial: " << ene_ini.diph_l       << " post_conv_dying: " << ene_postconv_dying.diph_l       << " post_conv_regular: " << ene_postconv_regular.diph_l       << " post_conv_small_nascent: " << ene_postconv_small_nascent.diph_l       << " post_diff_regular: " << ene_postdiff_regular.diph_l << " post_diff_small: " << ene_postdiff_small.diph_l       << finl;
  Cerr << "T_Budget, diph_v........ time: " << current_time << " initial: " << ene_ini.diph_v       << " post_conv_dying: " << ene_postconv_dying.diph_v       << " post_conv_regular: " << ene_postconv_regular.diph_v       << " post_conv_small_nascent: " << ene_postconv_small_nascent.diph_v       << " post_diff_regular: " << ene_postdiff_regular.diph_v << " post_diff_small: " << ene_postdiff_small.diph_v       << finl;
  Cerr << "T_Budget, diph_small.... time: " << current_time << " initial: " << ene_ini.diph_small   << " post_conv_dying: " << ene_postconv_dying.diph_small   << " post_conv_regular: " << ene_postconv_regular.diph_small   << " post_conv_small_nascent: " << ene_postconv_small_nascent.diph_small   << " post_diff_regular: " << ene_postdiff_regular.diph_small << " post_diff_small: " << ene_postdiff_small.diph_small   << finl;
  Cerr << "T_Budget, diph_regular.. time: " << current_time << " initial: " << ene_ini.diph_regular << " post_conv_dying: " << ene_postconv_dying.diph_regular << " post_conv_regular: " << ene_postconv_regular.diph_regular << " post_conv_small_nascent: " << ene_postconv_small_nascent.diph_regular << " post_diff_regular: " << ene_postdiff_regular.diph_regular << " post_diff_small: " << ene_postdiff_small.diph_regular << finl;
  Cerr << "T_Budget, diph_nascent.. time: " << current_time << " initial: " << ene_ini.diph_nascent << " post_conv_dying: " << ene_postconv_dying.diph_nascent << " post_conv_regular: " << ene_postconv_regular.diph_nascent << " post_conv_small_nascent: " << ene_postconv_small_nascent.diph_nascent << " post_diff_regular: " << ene_postdiff_regular.diph_nascent << " post_diff_small: " << ene_postdiff_small.diph_nascent << finl;
  Cerr << "T_Budget, diph_dying.... time: " << current_time << " initial: " << ene_ini.diph_dying   << " post_conv_dying: " << ene_postconv_dying.diph_dying   << " post_conv_regular: " << ene_postconv_regular.diph_dying   << " post_conv_small_nascent: " << ene_postconv_small_nascent.diph_dying   << " post_diff_regular: " << ene_postdiff_regular.diph_dying << " post_diff_small: " << ene_postdiff_small.diph_dying   << finl;

  CutCell_GlobalInfo Tmin = compute_Tmin_cut_cell(cut_field_temperature_, 1);
  CutCell_GlobalInfo Tmax = compute_Tmax_cut_cell(cut_field_temperature_, 1);
  Cerr << "T_MinMax, overall time: " << current_time << " Tmin: " << Tmin.overall << " Tmax: " << Tmax.overall << finl;
  Cerr << "T_MinMax, overall_l time: " << current_time << " Tmin: " << Tmin.overall_l << " Tmax: " << Tmax.overall_l << finl;
  Cerr << "T_MinMax, overall_v time: " << current_time << " Tmin: " << Tmin.overall_v << " Tmax: " << Tmax.overall_v << finl;
  Cerr << "T_MinMax, pure.......... time: " << current_time << " Tmin: " << Tmin.pure         << " Tmax: " << Tmax.pure         << finl;
  Cerr << "T_MinMax, diph_l........ time: " << current_time << " Tmin: " << Tmin.diph_l       << " Tmax: " << Tmax.diph_l       << finl;
  Cerr << "T_MinMax, diph_v........ time: " << current_time << " Tmin: " << Tmin.diph_v       << " Tmax: " << Tmax.diph_v       << finl;
  Cerr << "T_MinMax, diph_small.... time: " << current_time << " Tmin: " << Tmin.diph_small   << " Tmax: " << Tmax.diph_small   << finl;
  Cerr << "T_MinMax, diph_regular.. time: " << current_time << " Tmin: " << Tmin.diph_regular << " Tmax: " << Tmax.diph_regular << finl;
  Cerr << "T_MinMax, diph_nascent.. time: " << current_time << " Tmin: " << Tmin.diph_nascent << " Tmax: " << Tmax.diph_nascent << finl;
  Cerr << "T_MinMax, diph_dying.... time: " << current_time << " Tmin: " << Tmin.diph_dying   << " Tmax: " << Tmax.diph_dying   << finl;
}

void IJK_Thermal_cut_cell::rk3_sub_step(const int rk_step, const double total_timestep,
                                        const double time)
{
  if (debug_)
    Cerr << "Thermal Runge-Kutta3 time-step" << finl;
  update_thermal_properties();

  const Cut_field_vector& cut_field_velocity = ref_ijk_ft_cut_cell_->get_cut_field_velocity();
  calculer_dT_cut_cell(cut_field_velocity); // Note : ne fait rien

  const double current_time = ref_ijk_ft_cut_cell_->get_current_time();

  if (!conv_temperature_negligible_)
    {
      if (methode_temperature_remplissage_ == METHODE_TEMPERATURE_REMPLISSAGE::PONDERATION_VOISIN)
        {
          convective_correction_.calcule_temperature_remplissage_ponderation_voisin(cut_field_velocity, cut_field_temperature_);
        }
      else if (methode_temperature_remplissage_ == METHODE_TEMPERATURE_REMPLISSAGE::SEMI_LAGRANGIEN)
        {
          double fractional_timestep = compute_fractionnal_timestep_rk3(total_timestep, rk_step);
          convective_correction_.calcule_temperature_remplissage_semi_lagrangien(fractional_timestep, lambda_liquid_, lambda_vapour_, diffusive_correction_.flux_interface_ns_, cut_field_velocity, cut_field_temperature_);
        }
      else
        {
          Cerr << "Methode non reconnue pour le calcul de la temperature de remplissage." << finl;
          Process::exit();
        }
    }

  if (cut_cell_conv_scheme_ == CUT_CELL_CONV_SCHEME::INTERP_POINT_OU_QUICK_OU_AMONT_STENCIL)
    {
      convective_correction_.calcule_temperature_face_depuis_centre(lambda_liquid_, lambda_vapour_, diffusive_correction_.flux_interface_ns_, cut_field_temperature_, temperature_face_);
    }
  else if (cut_cell_conv_scheme_ == CUT_CELL_CONV_SCHEME::INTERP_FACETTE_OU_QUICK_OU_AMONT_STENCIL)
    {
      diffusive_correction_.calculer_flux_interface(methode_flux_interface_, scaled_distance_flux_interface_, lambda_liquid_, lambda_vapour_, interfacial_temperature_, interfacial_phin_ai_, cut_field_temperature_, ref_ijk_ft_cut_cell_, temperature_, temperature_ft_);

      convective_correction_.calcule_temperature_face_depuis_facette(lambda_liquid_, lambda_vapour_, interfacial_temperature_, interfacial_phin_ai_, cut_field_temperature_, ref_ijk_ft_cut_cell_, temperature_face_ft_, temperature_face_);
    }
  else
    {
      // Do nothing
    }

  const CutCell_GlobalInfo ene_ini = compute_global_energy_cut_cell(cut_field_temperature_, 0);

  if (!conv_temperature_negligible_)
    {
      convective_correction_.add_convection_dying_cells(convection_petites_cellules_, cut_field_velocity, cut_field_temperature_);
      cut_field_temperature_.echange_espace_virtuel(temperature_.ghost());
    }

  if (postraiter_champs_intermediaires_)
    {
      cut_field_temperature_post_dying_.copy_from(cut_field_temperature_);
    }
  const CutCell_GlobalInfo ene_postconv_dying = compute_global_energy_cut_cell(cut_field_temperature_, 0);

  compute_temperature_convection_cut_cell(cut_field_velocity);
  const CutCell_GlobalInfo d_ene_Conv = compute_global_energy_cut_cell(cut_field_d_temperature_, 1);
  ref_ijk_ft_cut_cell_->runge_kutta3_update_cut_cell_transport(cut_field_d_temperature_, cut_field_RK3_F_temperature_conv_, cut_field_temperature_, rk_step, total_timestep);
  cut_field_temperature_.echange_espace_virtuel(temperature_.ghost());

  if (postraiter_champs_intermediaires_)
    {
      cut_field_temperature_post_regular_.copy_from(cut_field_temperature_);
    }
  const CutCell_GlobalInfo ene_postconv_regular = compute_global_energy_cut_cell(cut_field_temperature_, 1);

  if (!conv_temperature_negligible_)
    {
      convective_correction_.add_convection_small_nascent_cells(convection_petites_cellules_, cut_field_velocity, cut_field_temperature_);
      cut_field_temperature_.echange_espace_virtuel(temperature_.ghost());
    }

  cut_field_temperature_post_small_.copy_from(cut_field_temperature_);
  const CutCell_GlobalInfo ene_postconv_small_nascent = compute_global_energy_cut_cell(cut_field_temperature_, 1);

  if (activate_diffusion_interface_)
    {
      diffusive_correction_.calculer_flux_interface(methode_flux_interface_, scaled_distance_flux_interface_, lambda_liquid_, lambda_vapour_, interfacial_temperature_, interfacial_phin_ai_, cut_field_temperature_, ref_ijk_ft_cut_cell_, temperature_, temperature_ft_);
    }
  else
    {
      Cerr << "Warning: The diffusion at the interface is not activated." << finl;
    }
  add_temperature_diffusion();
  const CutCell_GlobalInfo d_ene_Diffu = compute_global_energy_cut_cell(cut_field_d_temperature_, 1);
  ref_ijk_ft_cut_cell_->runge_kutta3_update_cut_cell_notransport(cut_field_d_temperature_, cut_field_RK3_F_temperature_diff_, cut_field_temperature_, rk_step, total_timestep);
  cut_field_temperature_.echange_espace_virtuel(temperature_.ghost());

  if (postraiter_champs_intermediaires_)
    {
      cut_field_temperature_post_diff_regular_.copy_from(cut_field_temperature_);
    }
  const CutCell_GlobalInfo ene_postdiff_regular = compute_global_energy_cut_cell(cut_field_temperature_, 1);

  if (!diff_temperature_negligible_ && correction_petites_cellules_diffusion_)
    {
      diffusive_correction_.add_diffusion_small_cells(diffusion_petites_cellules_, cut_field_temperature_post_small_, cut_field_temperature_);
      cut_field_temperature_.echange_espace_virtuel(temperature_.ghost());
    }

  const CutCell_GlobalInfo ene_postdiff_small = compute_global_energy_cut_cell(cut_field_temperature_, 1);

  // Recalcul du flux a l'interface pour l'estimation de la temperature de remplissage au prochain sous pas de temps
  diffusive_correction_.calculer_flux_interface(methode_flux_interface_, scaled_distance_flux_interface_, lambda_liquid_, lambda_vapour_, interfacial_temperature_, interfacial_phin_ai_, cut_field_temperature_, ref_ijk_ft_cut_cell_, temperature_, temperature_ft_);

  Cerr.precision(16);
  Cerr << "dT_budget, overall time: " << current_time << " dT_conv: " << d_ene_Conv.overall << " dT_diff: " << d_ene_Diffu.overall << finl;
  Cerr << "dT_budget, overall_l time: " << current_time << " dT_conv: " << d_ene_Conv.overall_l << " dT_diff: " << d_ene_Diffu.overall_l << finl;
  Cerr << "dT_budget, overall_v time: " << current_time << " dT_conv: " << d_ene_Conv.overall_v << " dT_diff: " << d_ene_Diffu.overall_v << finl;
  Cerr << "dT_budget, pure.......... time: " << current_time << " dT_conv: " << d_ene_Conv.pure         << " dT_diff: " << d_ene_Diffu.pure         << finl;
  Cerr << "dT_budget, diph_l........ time: " << current_time << " dT_conv: " << d_ene_Conv.diph_l       << " dT_diff: " << d_ene_Diffu.diph_l       << finl;
  Cerr << "dT_budget, diph_v........ time: " << current_time << " dT_conv: " << d_ene_Conv.diph_v       << " dT_diff: " << d_ene_Diffu.diph_v       << finl;
  Cerr << "dT_budget, diph_small.... time: " << current_time << " dT_conv: " << d_ene_Conv.diph_small   << " dT_diff: " << d_ene_Diffu.diph_small   << finl;
  Cerr << "dT_budget, diph_regular.. time: " << current_time << " dT_conv: " << d_ene_Conv.diph_regular << " dT_diff: " << d_ene_Diffu.diph_regular << finl;
  Cerr << "dT_budget, diph_nascent.. time: " << current_time << " dT_conv: " << d_ene_Conv.diph_nascent << " dT_diff: " << d_ene_Diffu.diph_nascent << finl;
  Cerr << "dT_budget, diph_dying.... time: " << current_time << " dT_conv: " << d_ene_Conv.diph_dying   << " dT_diff: " << d_ene_Diffu.diph_dying   << finl;

  Cerr << "T_Budget, overall time: " << current_time << " initial: " << ene_ini.overall << " post_conv_dying: " << ene_postconv_dying.overall << " post_conv_regular: " << ene_postconv_regular.overall << " post_conv_small_nascent: " << ene_postconv_small_nascent.overall << " post_diff_regular: " << ene_postdiff_regular.overall << " post_diff_small: " << ene_postdiff_small.overall << finl;
  Cerr << "T_Budget, overall_l time: " << current_time << " initial: " << ene_ini.overall_l << " post_conv_dying: " << ene_postconv_dying.overall_l << " post_conv_regular: " << ene_postconv_regular.overall_l << " post_conv_small_nascent: " << ene_postconv_small_nascent.overall_l << " post_diff_regular: " << ene_postdiff_regular.overall_l << " post_diff_small: " << ene_postdiff_small.overall_l << finl;
  Cerr << "T_Budget, overall_v time: " << current_time << " initial: " << ene_ini.overall_v << " post_conv_dying: " << ene_postconv_dying.overall_v << " post_conv_regular: " << ene_postconv_regular.overall_v << " post_conv_small_nascent: " << ene_postconv_small_nascent.overall_v << " post_diff_regular: " << ene_postdiff_regular.overall_v << " post_diff_small: " << ene_postdiff_small.overall_v << finl;
  Cerr << "T_Budget, pure.......... time: " << current_time << " initial: " << ene_ini.pure         << " post_conv_dying: " << ene_postconv_dying.pure         << " post_conv_regular: " << ene_postconv_regular.pure         << " post_conv_small_nascent: " << ene_postconv_small_nascent.pure         << " post_diff_regular: " << ene_postdiff_regular.pure << " post_diff_small: " << ene_postdiff_small.pure         << finl;
  Cerr << "T_Budget, diph_l........ time: " << current_time << " initial: " << ene_ini.diph_l       << " post_conv_dying: " << ene_postconv_dying.diph_l       << " post_conv_regular: " << ene_postconv_regular.diph_l       << " post_conv_small_nascent: " << ene_postconv_small_nascent.diph_l       << " post_diff_regular: " << ene_postdiff_regular.diph_l << " post_diff_small: " << ene_postdiff_small.diph_l       << finl;
  Cerr << "T_Budget, diph_v........ time: " << current_time << " initial: " << ene_ini.diph_v       << " post_conv_dying: " << ene_postconv_dying.diph_v       << " post_conv_regular: " << ene_postconv_regular.diph_v       << " post_conv_small_nascent: " << ene_postconv_small_nascent.diph_v       << " post_diff_regular: " << ene_postdiff_regular.diph_v << " post_diff_small: " << ene_postdiff_small.diph_v       << finl;
  Cerr << "T_Budget, diph_small.... time: " << current_time << " initial: " << ene_ini.diph_small   << " post_conv_dying: " << ene_postconv_dying.diph_small   << " post_conv_regular: " << ene_postconv_regular.diph_small   << " post_conv_small_nascent: " << ene_postconv_small_nascent.diph_small   << " post_diff_regular: " << ene_postdiff_regular.diph_small << " post_diff_small: " << ene_postdiff_small.diph_small   << finl;
  Cerr << "T_Budget, diph_regular.. time: " << current_time << " initial: " << ene_ini.diph_regular << " post_conv_dying: " << ene_postconv_dying.diph_regular << " post_conv_regular: " << ene_postconv_regular.diph_regular << " post_conv_small_nascent: " << ene_postconv_small_nascent.diph_regular << " post_diff_regular: " << ene_postdiff_regular.diph_regular << " post_diff_small: " << ene_postdiff_small.diph_regular << finl;
  Cerr << "T_Budget, diph_nascent.. time: " << current_time << " initial: " << ene_ini.diph_nascent << " post_conv_dying: " << ene_postconv_dying.diph_nascent << " post_conv_regular: " << ene_postconv_regular.diph_nascent << " post_conv_small_nascent: " << ene_postconv_small_nascent.diph_nascent << " post_diff_regular: " << ene_postdiff_regular.diph_nascent << " post_diff_small: " << ene_postdiff_small.diph_nascent << finl;
  Cerr << "T_Budget, diph_dying.... time: " << current_time << " initial: " << ene_ini.diph_dying   << " post_conv_dying: " << ene_postconv_dying.diph_dying   << " post_conv_regular: " << ene_postconv_regular.diph_dying   << " post_conv_small_nascent: " << ene_postconv_small_nascent.diph_dying   << " post_diff_regular: " << ene_postdiff_regular.diph_dying << " post_diff_small: " << ene_postdiff_small.diph_dying   << finl;

  CutCell_GlobalInfo Tmin = compute_Tmin_cut_cell(cut_field_temperature_, 1);
  CutCell_GlobalInfo Tmax = compute_Tmax_cut_cell(cut_field_temperature_, 1);
  Cerr << "T_MinMax, overall time: " << current_time << " Tmin: " << Tmin.overall << " Tmax: " << Tmax.overall << finl;
  Cerr << "T_MinMax, overall_l time: " << current_time << " Tmin: " << Tmin.overall_l << " Tmax: " << Tmax.overall_l << finl;
  Cerr << "T_MinMax, overall_v time: " << current_time << " Tmin: " << Tmin.overall_v << " Tmax: " << Tmax.overall_v << finl;
  Cerr << "T_MinMax, pure.......... time: " << current_time << " Tmin: " << Tmin.pure         << " Tmax: " << Tmax.pure         << finl;
  Cerr << "T_MinMax, diph_l........ time: " << current_time << " Tmin: " << Tmin.diph_l       << " Tmax: " << Tmax.diph_l       << finl;
  Cerr << "T_MinMax, diph_v........ time: " << current_time << " Tmin: " << Tmin.diph_v       << " Tmax: " << Tmax.diph_v       << finl;
  Cerr << "T_MinMax, diph_small.... time: " << current_time << " Tmin: " << Tmin.diph_small   << " Tmax: " << Tmax.diph_small   << finl;
  Cerr << "T_MinMax, diph_regular.. time: " << current_time << " Tmin: " << Tmin.diph_regular << " Tmax: " << Tmax.diph_regular << finl;
  Cerr << "T_MinMax, diph_nascent.. time: " << current_time << " Tmin: " << Tmin.diph_nascent << " Tmax: " << Tmax.diph_nascent << finl;
  Cerr << "T_MinMax, diph_dying.... time: " << current_time << " Tmin: " << Tmin.diph_dying   << " Tmax: " << Tmax.diph_dying   << finl;
}

void IJK_Thermal_cut_cell::compute_diffusion_increment()
{
  // Update d_temperature
  const int ni = cut_field_d_temperature_.pure_.ni();
  const int nj = cut_field_d_temperature_.pure_.nj();
  const int nk = cut_field_d_temperature_.pure_.nk();
  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              double rhocpV = rho_cp_(i,j,k) * vol_;

              const double ope = cut_field_div_coeff_grad_T_volume_.pure_(i,j,k);
              const double resu = ope/rhocpV;
              cut_field_d_temperature_.pure_(i,j,k) = resu ;
            }
        }
    }

  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_div_coeff_grad_T_volume_.get_cut_cell_disc();
  const double rho_l = ref_ijk_ft_cut_cell_->get_rho_l();
  const double rho_v = ref_ijk_ft_cut_cell_->get_rho_v();
  for (int n = 0; n < cut_cell_disc.get_n_loc(); n++)
    {
      // On pre-divise par le volume cartesien vol_ et non par par le vrai
      // volume car on travaille en volume*temperature pour le theoreme de Reynolds.
      double rhocpV_l = rho_l * cp_liquid_ * vol_;
      double rhocpV_v = rho_v * cp_vapour_ * vol_;

      const double ope_l = cut_field_div_coeff_grad_T_volume_.diph_l_(n);
      const double resu_l = ope_l/rhocpV_l;
      cut_field_d_temperature_.diph_l_(n) = resu_l;

      const double ope_v = cut_field_div_coeff_grad_T_volume_.diph_v_(n);
      const double resu_v = ope_v/rhocpV_v;
      cut_field_d_temperature_.diph_v_(n) = resu_v;
    }
}

void cut_cell_reinit_streamObj(std::ostringstream& streamObj, const double& param)
{
  streamObj.str("");
  streamObj.clear();
  streamObj << (double) param;
}


Nom IJK_Thermal_cut_cell::generate_expression_temperature_ini(const double& time_scope, const double x, const double y, const double z)
{
  const double rho_l = ref_ijk_ft_cut_cell_->get_rho_l();
  const double alpha_liq = lambda_liquid_ / (rho_l * cp_liquid_);
  std::ostringstream streamObj;
  Nom expression_T = "(";
  streamObj << delta_T_subcooled_overheated_;
  expression_T += streamObj.str().c_str();
  expression_T += ")-(";
  expression_T += streamObj.str().c_str();
  expression_T += ")*(";
  cut_cell_reinit_streamObj(streamObj, single_centred_bubble_radius_ini_);
  expression_T += streamObj.str().c_str();
  Nom expression_tmp = "sqrt((x-(";
  cut_cell_reinit_streamObj(streamObj, (signbit(x) ? x-1e-16 : x+1e-16));
  expression_tmp += streamObj.str().c_str();
  expression_tmp += "))^2+(y-(";
  cut_cell_reinit_streamObj(streamObj, (signbit(y) ? y-1e-16 : y+1e-16));
  expression_tmp += streamObj.str().c_str();
  expression_tmp += "))^2+(z-(";
  cut_cell_reinit_streamObj(streamObj, (signbit(z) ? z-1e-16 : z+1e-16));
  expression_tmp += streamObj.str().c_str();
  expression_tmp += "))^2)";
  expression_T += "/";
  expression_T += expression_tmp;
  expression_T += "*(1-erf((";
  expression_T += expression_tmp;
  expression_T += "-";
  cut_cell_reinit_streamObj(streamObj, single_centred_bubble_radius_ini_);
  expression_T += streamObj.str().c_str();
  expression_T += ")/(2.*sqrt(";
  streamObj.str("");
  streamObj.clear();
  cut_cell_reinit_streamObj(streamObj, ((alpha_liq * time_scope) + 1e-16));
  expression_T += streamObj.str().c_str();
  expression_T += ")))))";
  return expression_T;
}

void IJK_Thermal_cut_cell::calculer_dT_cut_cell(const Cut_field_vector& cut_field_velocity)
{
  // Note : Cette fonction regroupe tous les elements de calculer_dT qui ne sont pas utilisees en cut cell

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

  /*
   * Compute gradients and hessian of the temperature after the ghost fluid extension
   */
  if (debug_)
    Cerr << "Compute temperature derivatives" << finl;
  compute_temperature_gradient_elem();
  compute_temperature_hessian_diag_elem();
  compute_temperature_hessian_cross_elem();

  /*
   * Compute sub-problems (For Subresolution Child classes only !)
   */
  if (debug_)
    Cerr << "Compute thermal subproblems" << finl;
  compute_thermal_subproblems();

  /*
   * Interpolate a value for the QUICK SCHEME (first call)
   */
  if (debug_)
    Cerr << "Compute temperature mixed cell centres" << finl;
  compute_temperature_cell_centres(0);

  /*
   * Convective and Diffusive fluxes
   */
  if (debug_)
    Cerr << "Compute thermal convective and diffusive fluxes from subproblems" << finl;
  compute_convective_diffusive_fluxes_face_centre();

  if (debug_)
    Cerr << "Prepare ij fluxes" << finl;
  if (!conv_temperature_negligible_ || !diff_temperature_negligible_)
    prepare_ij_fluxes_k_layers();


  double nb_diam_upstream_velocity = ref_ijk_ft_cut_cell_->get_nb_diam_upstream();
  if (nb_diam_upstream_ == 0.)
    nb_diam_upstream_ = nb_diam_upstream_velocity;
  if (upstream_temperature_ > -1e20 && ref_ijk_ft_cut_cell_->get_vitesse_upstream() > -1e20)
    force_upstream_temperature(temperature_, upstream_temperature_,
                               ref_ijk_ft_cut_cell_->get_interface(), nb_diam_upstream_,
                               ref_ijk_ft_cut_cell_->get_upstream_dir(), ref_ijk_ft_cut_cell_->get_direction_gravite(),
                               ref_ijk_ft_cut_cell_->get_upstream_stencil());

  return;
}

void IJK_Thermal_cut_cell::compute_temperature_convection_cut_cell(const Cut_field_vector& cut_field_velocity)
{
  static Stat_Counter_Id cnt_conv_temp = statistiques().new_counter(1, "FT convection rho");
  statistiques().begin_count(cnt_conv_temp);
  if (conv_temperature_negligible_)
    {
      cut_field_d_temperature_.pure_.data()=0;
      cut_field_d_temperature_.set_valeur_cellules_diphasiques(0);
      u_T_convective_volume_.data() = 0;
    }
  else
    {
      temperature_convection_op_.calculer_cut_cell(false, cut_cell_conv_scheme_, cut_field_temperature_, cut_field_velocity, temperature_face_, cut_cell_flux_convection_, cut_field_d_temperature_);
      const int ni = cut_field_d_temperature_.pure_.ni();
      const int nj = cut_field_d_temperature_.pure_.nj();
      const int nk = cut_field_d_temperature_.pure_.nk();
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  u_T_convective_volume_(i,j,k) = cut_field_d_temperature_.pure_(i,j,k);
                  const double resu = cut_field_d_temperature_.pure_(i,j,k) / vol_;
                  cut_field_d_temperature_.pure_(i,j,k) = resu;
                  //if (liste_post_instantanes_.contient_("U_T_CONVECTIVE"))
                  //  {
                  //    u_T_convective_(i,j,k) = resu;
                  //  }
                }
            }
        }

      const Cut_cell_FT_Disc& cut_cell_disc = cut_field_d_temperature_.get_cut_cell_disc();
      for (int n = 0; n < cut_cell_disc.get_n_loc(); n++)
        {
          // On pre-divise par le volume cartesien vol_ et non par par le vrai
          // volume car on travaille en volume*temperature pour le theoreme de Reynolds.
          double V_l = vol_;
          double V_v = vol_;

          const double resu_l = (V_l == 0) ? 0. : cut_field_d_temperature_.diph_l_(n)/V_l;
          cut_field_d_temperature_.diph_l_(n) = resu_l;

          const double resu_v = (V_v == 0) ? 0. : cut_field_d_temperature_.diph_v_(n)/V_v;
          cut_field_d_temperature_.diph_v_(n) = resu_v;
        }
    }
  statistiques().end_count(cnt_conv_temp);
  DebogIJK::verifier("op_conv(rho)", cut_field_d_temperature_.pure_);
  return;
}


void IJK_Thermal_cut_cell::add_temperature_diffusion()
{
  lambda_.echange_espace_virtuel(lambda_.ghost());
  DebogIJK::verifier("lambda", lambda_);

  temperature_diffusion_op_.set_lambda(lambda_);

  if (boundary_conditions_.get_bctype_k_min() == Boundary_Conditions_Thermique::Paroi_Temperature_imposee)
    {
      const double T_paroi_impose_kmin = boundary_conditions_.get_temperature_kmin();
      double lambda_de_t_paroi_kmin = lambda_liquid_;
      // calculer ici div(lambda*grad(T))*volume)
      calculer_flux_thermique_bord(temperature_, lambda_de_t_paroi_kmin,
                                   T_paroi_impose_kmin, boundary_flux_kmin_, 0 /* boundary kmin */);
    }
  else if (boundary_conditions_.get_bctype_k_min() == Boundary_Conditions_Thermique::Paroi_Flux_impose)
    {
      const double flux_paroi_impose_kmin = boundary_conditions_.get_flux_kmin();
      imposer_flux_thermique_bord(temperature_,
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
      calculer_flux_thermique_bord(temperature_, lambda_de_t_paroi_kmax,
                                   T_paroi_impose_kmax, boundary_flux_kmax_, 1 /* boundary kmax */);
    }
  else if (boundary_conditions_.get_bctype_k_max() == Boundary_Conditions_Thermique::Paroi_Flux_impose)
    {
      const double flux_paroi_impose_kmax = boundary_conditions_.get_flux_kmax();
      imposer_flux_thermique_bord(temperature_,
                                  flux_paroi_impose_kmax, boundary_flux_kmax_, 1 /* boundary kmax */);
      // Cerr << "not coded yet" << finl;
      // Process::exit();
    }
  else
    {
      // Perio... not needed. The flux on the "boundary" (in that case , it's not truely a boundary) will be computed as inside...
    }
  temperature_.echange_espace_virtuel(temperature_.ghost());
  DebogIJK::verifier("temp", temperature_);

  if (diff_temperature_negligible_)
    {
      cut_field_d_temperature_.pure_.data()=0;
      cut_field_d_temperature_.set_valeur_cellules_diphasiques(0);
    }
  else
    {
      // Performance counters:
      static Stat_Counter_Id cnt_diff_temp = statistiques().new_counter(1, "FT diffusion temperature");
      statistiques().begin_count(cnt_diff_temp);
      /*
       * Correct the diffusive fluxes here or in the operator ?
       */
      temperature_diffusion_op_.calculer_cut_cell(false,
                                                  cut_field_temperature_,
                                                  cut_cell_flux_diffusion_,
                                                  cut_field_div_coeff_grad_T_volume_,
                                                  boundary_flux_kmin_,
                                                  boundary_flux_kmax_);
      if (etalement_diffusion_ == ETALEMENT_DIFFUSION::FLUX_INTERFACE_UNIQUEMENT || etalement_diffusion_ == ETALEMENT_DIFFUSION::FLUX_INTERFACE_ET_DIVERGENCE)
        {
          diffusive_correction_.ajout_flux_interface_a_divergence_etale(cut_field_div_coeff_grad_T_volume_);
        }
      else
        {
          diffusive_correction_.ajout_flux_interface_a_divergence_simple(cut_field_div_coeff_grad_T_volume_);
        }
      if (etalement_diffusion_ == ETALEMENT_DIFFUSION::DIVERGENCE_UNIQUEMENT || etalement_diffusion_ == ETALEMENT_DIFFUSION::FLUX_INTERFACE_ET_DIVERGENCE)
        {
          diffusive_correction_.etalement_divergence_flux_diffusifs(cut_field_div_coeff_grad_T_volume_, cut_field_div_coeff_grad_T_volume_temp_);
        }
      compute_diffusion_increment();
      statistiques().end_count(cnt_diff_temp);
      DebogIJK::verifier("div_coeff_grad_T_volume_", div_coeff_grad_T_volume_);
    }
}

CutCell_GlobalInfo IJK_Thermal_cut_cell::compute_global_energy_cut_cell(Cut_field_scalar& cut_field_temperature, bool next)
{
  double global_energy_overall = 0.;
  double global_energy_overall_l = 0.;
  double global_energy_overall_v = 0.;
  double global_energy_pure = 0.;
  double global_energy_diph_l = 0.;
  double global_energy_diph_v = 0.;
  double global_energy_diph_small = 0.;
  double global_energy_diph_regular = 0.;
  double global_energy_diph_nascent = 0.;
  double global_energy_diph_dying = 0.;
  double count_overall = 0.;
  double count_overall_l = 0.;
  double count_overall_v = 0.;
  double count_pure = 0.;
  double count_diph_l = 0.;
  double count_diph_v = 0.;
  double count_diph_small = 0.;
  double count_diph_regular = 0.;
  double count_diph_nascent = 0.;
  double count_diph_dying = 0.;
  const IJK_Field_double& indic_old = ref_ijk_ft_cut_cell_->itfce().I();
  const IJK_Field_double& indic_next = ref_ijk_ft_cut_cell_->itfce().In();
  const IJK_Field_double& indic = next ? indic_next : indic_old;
  const double rhocpl = get_rhocp_l();
  const double rhocpv = get_rhocp_v();
  const int nx = cut_field_temperature.pure_.ni();
  const int ny = cut_field_temperature.pure_.nj();
  const int nz = cut_field_temperature.pure_.nk();
  // To be sure we're on a regular mesh
  assert(indic.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_K) >0);
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();
  for (int k=0; k < nz ; k++)
    {
      for (int j=0; j< ny; j++)
        {
          for (int i=0; i < nx; i++)
            {
              double chi_l = indic(i,j,k);
              double old_chi_l = indic_old(i,j,k);
              double next_chi_l = indic_next(i,j,k);
              int n = cut_cell_disc.get_n(i,j,k);
              if (n < 0)
                {
                  global_energy_overall += (chi_l * rhocpl + (1.- chi_l) * rhocpv) * cut_field_temperature.pure_(i,j,k);
                  count_overall += 1;

                  global_energy_overall_l += chi_l * rhocpl * cut_field_temperature.pure_(i,j,k);
                  count_overall_l += chi_l;

                  global_energy_overall_v += (1.- chi_l) * rhocpv * cut_field_temperature.pure_(i,j,k);
                  count_overall_v += (1.- chi_l);

                  global_energy_pure += (chi_l * rhocpl + (1.- chi_l) * rhocpv) * cut_field_temperature.pure_(i,j,k);
                  count_pure += 1;
                }
              else
                {
                  global_energy_overall += chi_l * rhocpl * cut_field_temperature.diph_l_(n);
                  global_energy_overall += (1.- chi_l) * rhocpv * cut_field_temperature.diph_v_(n);
                  count_overall += 1;

                  global_energy_overall_l += chi_l * rhocpl * cut_field_temperature.diph_l_(n);
                  count_overall_l += 1;

                  global_energy_overall_v += (1.- chi_l) * rhocpv * cut_field_temperature.diph_v_(n);
                  count_overall_v += 1;

                  global_energy_diph_l += chi_l * rhocpl * cut_field_temperature.diph_l_(n);
                  count_diph_l += chi_l;

                  global_energy_diph_v += (1.- chi_l) * rhocpv * cut_field_temperature.diph_v_(n);
                  count_diph_v += (1.- chi_l);

                  if (ref_ijk_ft_cut_cell_->itfce().next_below_small_threshold_for_phase(1, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      global_energy_diph_small += chi_l * rhocpl * cut_field_temperature.diph_l_(n);
                      count_diph_small += chi_l;

                      global_energy_diph_regular += (1.- chi_l) * rhocpv * cut_field_temperature.diph_v_(n);
                      count_diph_regular += (1.- chi_l);
                    }
                  else if (ref_ijk_ft_cut_cell_->itfce().next_below_small_threshold_for_phase(0, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      global_energy_diph_small += (1.- chi_l) * rhocpv * cut_field_temperature.diph_v_(n);
                      count_diph_small += (1.- chi_l);

                      global_energy_diph_regular += chi_l * rhocpl * cut_field_temperature.diph_l_(n);
                      count_diph_regular += chi_l;
                    }
                  else
                    {
                      global_energy_diph_regular += chi_l * rhocpl * cut_field_temperature.diph_l_(n);
                      global_energy_diph_regular += (1.- chi_l) * rhocpv * cut_field_temperature.diph_v_(n);
                      count_diph_regular += 1;
                    }

                  if (ref_ijk_ft_cut_cell_->itfce().devient_pure(indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      int phase_dying = (int)(1 - indic_next(i,j,k));
                      if (phase_dying == 1)
                        {
                          global_energy_diph_dying += old_chi_l * rhocpl * cut_field_temperature.diph_l_(n);
                          count_diph_dying += old_chi_l;
                        }
                      else
                        {
                          global_energy_diph_dying += (1.- old_chi_l) * rhocpv * cut_field_temperature.diph_v_(n);
                          count_diph_dying += (1.- old_chi_l);
                        }
                    }

                  if (ref_ijk_ft_cut_cell_->itfce().devient_diphasique(indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      int phase_nascent = (int)(1 - indic_old(i,j,k));
                      if (phase_nascent == 1)
                        {
                          global_energy_diph_nascent += next_chi_l * rhocpl * cut_field_temperature.diph_l_(n);
                          count_diph_nascent += next_chi_l;
                        }
                      else
                        {
                          global_energy_diph_nascent += (1.- next_chi_l) * rhocpv * cut_field_temperature.diph_v_(n);
                          count_diph_nascent += (1.- next_chi_l);
                        }
                    }
                }
            }
        }
    }
  count_overall = mp_sum(count_overall);
  count_overall_l = mp_sum(count_overall_l);
  count_overall_v = mp_sum(count_overall_v);
  count_pure = mp_sum(count_pure);
  count_diph_l = mp_sum(count_diph_l);
  count_diph_v = mp_sum(count_diph_v);
  count_diph_small = mp_sum(count_diph_small);
  count_diph_regular = mp_sum(count_diph_regular);
  count_diph_nascent = mp_sum(count_diph_nascent);
  count_diph_dying = mp_sum(count_diph_dying);
  assert(count_overall == cut_field_temperature.pure_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_I)
         *cut_field_temperature.pure_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_J)
         *cut_field_temperature.pure_.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K));
  global_energy_overall      = mp_sum(global_energy_overall);
  global_energy_overall_l    = mp_sum(global_energy_overall_l);
  global_energy_overall_v    = mp_sum(global_energy_overall_v);
  global_energy_pure         = mp_sum(global_energy_pure);
  global_energy_diph_l       = mp_sum(global_energy_diph_l);
  global_energy_diph_v       = mp_sum(global_energy_diph_v);
  global_energy_diph_small   = mp_sum(global_energy_diph_small);
  global_energy_diph_regular = mp_sum(global_energy_diph_regular);
  global_energy_diph_nascent = mp_sum(global_energy_diph_nascent);
  global_energy_diph_dying   = mp_sum(global_energy_diph_dying);
  return {global_energy_overall, global_energy_overall_l, global_energy_overall_v, global_energy_pure, global_energy_diph_l, global_energy_diph_v, global_energy_diph_small, global_energy_diph_regular, global_energy_diph_nascent, global_energy_diph_dying};
}

CutCell_GlobalInfo IJK_Thermal_cut_cell::compute_Tmin_cut_cell(Cut_field_scalar& cut_field_temperature, bool next)
{
  double Tmin_overall = 1.e20;
  double Tmin_overall_l = 1.e20;
  double Tmin_overall_v = 1.e20;
  double Tmin_pure = 1.e20;
  double Tmin_diph_l = 1.e20;
  double Tmin_diph_v = 1.e20;
  double Tmin_diph_small = 1.e20;
  double Tmin_diph_regular = 1.e20;
  double Tmin_diph_nascent = 1.e20;
  double Tmin_diph_dying = 1.e20;
  int n_T_min = -1;
  const IJK_Field_double& indic_old = ref_ijk_ft_cut_cell_->itfce().I();
  const IJK_Field_double& indic_next = ref_ijk_ft_cut_cell_->itfce().In();
  const IJK_Field_double& indic = next ? indic_next : indic_old;
  const int nx = cut_field_temperature.pure_.ni();
  const int ny = cut_field_temperature.pure_.nj();
  const int nz = cut_field_temperature.pure_.nk();
  // To be sure we're on a regular mesh
  assert(indic.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_K) >0);
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();
  for (int k=0; k < nz ; k++)
    {
      for (int j=0; j< ny; j++)
        {
          for (int i=0; i < nx; i++)
            {
              double chi_l = indic(i,j,k);
              int n = cut_cell_disc.get_n(i,j,k);
              if (n < 0)
                {
                  Tmin_overall = std::min(Tmin_overall, cut_field_temperature.pure_(i,j,k));
                  Tmin_overall_l = (chi_l == 0) ? Tmin_overall_l : std::min(Tmin_overall_l, cut_field_temperature.pure_(i,j,k));
                  Tmin_overall_v = (chi_l == 1) ? Tmin_overall_v : std::min(Tmin_overall_v, cut_field_temperature.pure_(i,j,k));
                  Tmin_pure = std::min(Tmin_pure, cut_field_temperature.pure_(i,j,k));
                  if (Tmin_overall == cut_field_temperature.pure_(i,j,k))
                    {
                      n_T_min = n;
                    }
                }
              else
                {
                  // Excluding the value of the phase in dying cells
                  bool exclude_l = (ref_ijk_ft_cut_cell_->itfce().devient_pure(indic_old(i,j,k), indic_next(i,j,k))) && ((int)(1 - indic_next(i,j,k)) == 1.);
                  bool exclude_v = (ref_ijk_ft_cut_cell_->itfce().devient_pure(indic_old(i,j,k), indic_next(i,j,k))) && ((int)(1 - indic_next(i,j,k)) == 0.);

                  Tmin_overall = exclude_l ? Tmin_overall : std::min(Tmin_overall, cut_field_temperature.diph_l_(n));
                  Tmin_overall = exclude_v ? Tmin_overall : std::min(Tmin_overall, cut_field_temperature.diph_v_(n));

                  Tmin_overall_l = exclude_l ? Tmin_overall_l : std::min(Tmin_overall_l, cut_field_temperature.diph_l_(n));
                  Tmin_overall_v = exclude_v ? Tmin_overall_v : std::min(Tmin_overall_v, cut_field_temperature.diph_v_(n));

                  if (Tmin_overall == cut_field_temperature.diph_l_(n))
                    {
                      n_T_min = n;
                    }
                  if (Tmin_overall == cut_field_temperature.diph_v_(n))
                    {
                      n_T_min = n;
                    }

                  Tmin_diph_l = exclude_l ? Tmin_diph_l : std::min(Tmin_diph_l, cut_field_temperature.diph_l_(n));
                  Tmin_diph_v = exclude_v ? Tmin_diph_v : std::min(Tmin_diph_v, cut_field_temperature.diph_v_(n));

                  if (ref_ijk_ft_cut_cell_->itfce().next_below_small_threshold_for_phase(1, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      Tmin_diph_small = exclude_l ? Tmin_diph_small : std::min(Tmin_diph_small, cut_field_temperature.diph_l_(n));
                      Tmin_diph_regular = exclude_v ? Tmin_diph_regular : std::min(Tmin_diph_regular, cut_field_temperature.diph_v_(n));
                    }
                  else if (ref_ijk_ft_cut_cell_->itfce().next_below_small_threshold_for_phase(0, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      Tmin_diph_small = exclude_v ? Tmin_diph_small : std::min(Tmin_diph_small, cut_field_temperature.diph_v_(n));
                      Tmin_diph_regular = exclude_l ? Tmin_diph_regular : std::min(Tmin_diph_regular, cut_field_temperature.diph_l_(n));
                    }
                  else
                    {
                      Tmin_diph_regular = exclude_l ? Tmin_diph_regular : std::min(Tmin_diph_regular, cut_field_temperature.diph_l_(n));
                      Tmin_diph_regular = exclude_v ? Tmin_diph_regular : std::min(Tmin_diph_regular, cut_field_temperature.diph_v_(n));
                    }

                  if (ref_ijk_ft_cut_cell_->itfce().devient_pure(indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      int phase_dying = (int)(1 - indic_next(i,j,k));
                      if (phase_dying == 1)
                        {
                          Tmin_diph_dying = std::min(Tmin_diph_dying, cut_field_temperature.diph_l_(n));
                        }
                      else
                        {
                          Tmin_diph_dying = std::min(Tmin_diph_dying, cut_field_temperature.diph_v_(n));
                        }
                    }

                  if (ref_ijk_ft_cut_cell_->itfce().devient_diphasique(indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      int phase_nascent = (int)(1 - indic_old(i,j,k));
                      if (phase_nascent == 1)
                        {
                          Tmin_diph_nascent = std::min(Tmin_diph_nascent, cut_field_temperature.diph_l_(n));
                        }
                      else
                        {
                          Tmin_diph_nascent = std::min(Tmin_diph_nascent, cut_field_temperature.diph_v_(n));
                        }
                    }
                }
            }
        }
    }
  Tmin_overall      = Process::mp_min(Tmin_overall);
  Tmin_overall_l      = Process::mp_min(Tmin_overall_l);
  Tmin_overall_v      = Process::mp_min(Tmin_overall_v);
  Tmin_pure         = Process::mp_min(Tmin_pure);
  Tmin_diph_l       = Process::mp_min(Tmin_diph_l);
  Tmin_diph_v       = Process::mp_min(Tmin_diph_v);
  Tmin_diph_small   = Process::mp_min(Tmin_diph_small);
  Tmin_diph_regular   = Process::mp_min(Tmin_diph_regular);
  Tmin_diph_nascent = Process::mp_min(Tmin_diph_nascent);
  Tmin_diph_dying   = Process::mp_min(Tmin_diph_dying);
  Cerr << " n Tmin " << n_T_min << finl;
  return {Tmin_overall, Tmin_overall_l, Tmin_overall_v, Tmin_pure, Tmin_diph_l, Tmin_diph_v, Tmin_diph_small, Tmin_diph_regular, Tmin_diph_nascent, Tmin_diph_dying};
}

CutCell_GlobalInfo IJK_Thermal_cut_cell::compute_Tmax_cut_cell(Cut_field_scalar& cut_field_temperature, bool next)
{
  double Tmax_overall = -1.e20;
  double Tmax_overall_l = -1.e20;
  double Tmax_overall_v = -1.e20;
  double Tmax_pure = -1.e20;
  double Tmax_diph_l = -1.e20;
  double Tmax_diph_v = -1.e20;
  double Tmax_diph_small = -1.e20;
  double Tmax_diph_regular = -1.e20;
  double Tmax_diph_nascent = -1.e20;
  double Tmax_diph_dying = -1.e20;
  const IJK_Field_double& indic_old = ref_ijk_ft_cut_cell_->itfce().I();
  const IJK_Field_double& indic_next = ref_ijk_ft_cut_cell_->itfce().In();
  const IJK_Field_double& indic = next ? indic_next : indic_old;
  const int nx = cut_field_temperature.pure_.ni();
  const int ny = cut_field_temperature.pure_.nj();
  const int nz = cut_field_temperature.pure_.nk();
  // To be sure we're on a regular mesh
  assert(indic.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_K) >0);
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();
  for (int k=0; k < nz ; k++)
    {
      for (int j=0; j< ny; j++)
        {
          for (int i=0; i < nx; i++)
            {
              double chi_l = indic(i,j,k);
              int n = cut_cell_disc.get_n(i,j,k);
              if (n < 0)
                {
                  Tmax_overall = std::max(Tmax_overall, cut_field_temperature.pure_(i,j,k));
                  Tmax_overall_l = (chi_l == 0) ? Tmax_overall_l : std::max(Tmax_overall_l, cut_field_temperature.pure_(i,j,k));
                  Tmax_overall_v = (chi_l == 1) ? Tmax_overall_v : std::max(Tmax_overall_v, cut_field_temperature.pure_(i,j,k));
                  Tmax_pure = std::max(Tmax_pure, cut_field_temperature.pure_(i,j,k));
                }
              else
                {
                  // Excluding the value of the phase in dying cells
                  bool exclude_l = (ref_ijk_ft_cut_cell_->itfce().devient_pure(indic_old(i,j,k), indic_next(i,j,k))) && ((int)(1 - indic_next(i,j,k)) == 1.);
                  bool exclude_v = (ref_ijk_ft_cut_cell_->itfce().devient_pure(indic_old(i,j,k), indic_next(i,j,k))) && ((int)(1 - indic_next(i,j,k)) == 0.);

                  Tmax_overall = exclude_l ? Tmax_overall : std::max(Tmax_overall, cut_field_temperature.diph_l_(n));
                  Tmax_overall = exclude_v ? Tmax_overall : std::max(Tmax_overall, cut_field_temperature.diph_v_(n));

                  Tmax_overall_l = exclude_l ? Tmax_overall_l : std::max(Tmax_overall_l, cut_field_temperature.diph_l_(n));
                  Tmax_overall_v = exclude_v ? Tmax_overall_v : std::max(Tmax_overall_v, cut_field_temperature.diph_v_(n));

                  Tmax_diph_l = exclude_l ? Tmax_diph_l : std::max(Tmax_diph_l, cut_field_temperature.diph_l_(n));
                  Tmax_diph_v = exclude_v ? Tmax_diph_v : std::max(Tmax_diph_v, cut_field_temperature.diph_v_(n));

                  if (ref_ijk_ft_cut_cell_->itfce().next_below_small_threshold_for_phase(1, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      Tmax_diph_small = exclude_l ? Tmax_diph_small : std::max(Tmax_diph_small, cut_field_temperature.diph_l_(n));
                      Tmax_diph_regular = exclude_v ? Tmax_diph_regular : std::max(Tmax_diph_regular, cut_field_temperature.diph_v_(n));
                    }
                  else if (ref_ijk_ft_cut_cell_->itfce().next_below_small_threshold_for_phase(0, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      Tmax_diph_small = exclude_v ? Tmax_diph_small : std::max(Tmax_diph_small, cut_field_temperature.diph_v_(n));
                      Tmax_diph_regular = exclude_l ? Tmax_diph_regular : std::max(Tmax_diph_regular, cut_field_temperature.diph_l_(n));
                    }
                  else
                    {
                      Tmax_diph_regular = exclude_l ? Tmax_diph_regular : std::max(Tmax_diph_regular, cut_field_temperature.diph_l_(n));
                      Tmax_diph_regular = exclude_v ? Tmax_diph_regular : std::max(Tmax_diph_regular, cut_field_temperature.diph_v_(n));
                    }

                  if (ref_ijk_ft_cut_cell_->itfce().devient_pure(indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      int phase_dying = (int)(1 - indic_next(i,j,k));
                      if (phase_dying == 1)
                        {
                          Tmax_diph_dying = std::max(Tmax_diph_dying, cut_field_temperature.diph_l_(n));
                        }
                      else
                        {
                          Tmax_diph_dying = std::max(Tmax_diph_dying, cut_field_temperature.diph_v_(n));
                        }
                    }

                  if (ref_ijk_ft_cut_cell_->itfce().devient_diphasique(indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      int phase_nascent = (int)(1 - indic_old(i,j,k));
                      if (phase_nascent == 1)
                        {
                          Tmax_diph_nascent = std::max(Tmax_diph_nascent, cut_field_temperature.diph_l_(n));
                        }
                      else
                        {
                          Tmax_diph_nascent = std::max(Tmax_diph_nascent, cut_field_temperature.diph_v_(n));
                        }
                    }
                }
            }
        }
    }
  Tmax_overall      = Process::mp_max(Tmax_overall);
  Tmax_overall_l      = Process::mp_max(Tmax_overall_l);
  Tmax_overall_v      = Process::mp_max(Tmax_overall_v);
  Tmax_pure         = Process::mp_max(Tmax_pure);
  Tmax_diph_l       = Process::mp_max(Tmax_diph_l);
  Tmax_diph_v       = Process::mp_max(Tmax_diph_v);
  Tmax_diph_small   = Process::mp_max(Tmax_diph_small);
  Tmax_diph_regular   = Process::mp_max(Tmax_diph_regular);
  Tmax_diph_nascent = Process::mp_max(Tmax_diph_nascent);
  Tmax_diph_dying   = Process::mp_max(Tmax_diph_dying);
  return {Tmax_overall, Tmax_overall_l, Tmax_overall_v, Tmax_pure, Tmax_diph_l, Tmax_diph_v, Tmax_diph_small, Tmax_diph_regular, Tmax_diph_nascent, Tmax_diph_dying};
}

double IJK_Thermal_cut_cell::compute_rho_cp_u_mean(const IJK_Field_double& vx)
{
  double rho_cp_u_mean = 0.;
  rho_cp_u_mean = calculer_rho_cp_u_moyen(vx, rho_cp_, vx, 0., 1);
  return rho_cp_u_mean;
}

double IJK_Thermal_cut_cell::get_rho_cp_ijk(int i, int j, int k) const
{
  double rho_cp = 0.;
  rho_cp = rho_cp_(i,j,k);
  return rho_cp;
}

double IJK_Thermal_cut_cell::get_rho_cp_u_ijk(const IJK_Field_double& vx, int i, int j, int k) const
{
  return get_rho_cp_ijk(i,j,k) * vx(i,j,k);
}

double IJK_Thermal_cut_cell::get_div_lambda_ijk(int i, int j, int k) const
{
  return (lambda_(i+1,j,k)-lambda_(i-1,j,k));
}

double IJK_Thermal_cut_cell::compute_temperature_dimensionless_theta_mean(const IJK_Field_double& vx)
{
  double theta_adim_moy = 0.;
  theta_adim_moy = calculer_temperature_adimensionnelle_theta_moy(vx, temperature_adimensionnelle_theta_, rho_cp_, vx, 0., 1);
  return theta_adim_moy;
}

void IJK_Thermal_cut_cell::correct_any_temperature_fields_for_eulerian_fluxes(IJK_Field_double& temperature)
{
  const int ni = temperature.ni();
  const int nj = temperature.nj();
  const int nk = temperature.nk();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double indic = ref_ijk_ft_cut_cell_->itfce().I(i,j,k);
          if (fabs(indic) < LIQUID_INDICATOR_TEST) // Mixed cells and pure vapour cells
            temperature(i,j,k) = 0.;
        }
  temperature.echange_espace_virtuel(temperature.ghost());
}


void IJK_Thermal_cut_cell::compare_temperature_fields(const IJK_Field_double& temperature,
                                                      const IJK_Field_double& temperature_ana,
                                                      IJK_Field_double& error_temperature_ana,
                                                      IJK_Field_double& error_temperature_ana_rel)
{
  const int ni = temperature.ni();
  const int nj = temperature.nj();
  const int nk = temperature.nk();
  error_temperature_ana.data() = 0.;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          error_temperature_ana(i,j,k) = (temperature(i,j,k) - temperature_ana(i,j,k));
          error_temperature_ana_rel(i,j,k) = error_temperature_ana(i,j,k) / (temperature_ana(i,j,k) + 1.e-16) * 100.;
        }
}

void IJK_Thermal_cut_cell::evaluate_total_liquid_absolute_parameter(const IJK_Field_double& field,
                                                                    double& total_parameter)
{
  const int ni = field.ni();
  const int nj = field.nj();
  const int nk = field.nk();
  total_parameter = 0;
  double liquid_volume = 0.;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double indic = ref_ijk_ft_cut_cell_->itfce().I(i,j,k);
          liquid_volume += (vol_ * indic);
          total_parameter += abs(field(i,j,k)) * (vol_ * indic);
        }
  total_parameter = mp_sum(total_parameter);
  liquid_volume = mp_sum(liquid_volume);
  total_parameter = total_parameter / liquid_volume;
}

void IJK_Thermal_cut_cell::evaluate_total_liquid_parameter_squared(const IJK_Field_double& field,
                                                                   double& total_parameter)
{
  const int ni = field.ni();
  const int nj = field.nj();
  const int nk = field.nk();
  total_parameter = 0;
  double liquid_volume = 0.;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double indic = ref_ijk_ft_cut_cell_->itfce().I(i,j,k);
          liquid_volume += (vol_ * indic);
          total_parameter += pow(field(i,j,k) * (vol_ * indic), 2);
        }
  total_parameter = mp_sum(total_parameter);
  liquid_volume = mp_sum(liquid_volume);
  total_parameter = total_parameter / pow(liquid_volume, 2);
}

void IJK_Thermal_cut_cell::correct_any_temperature_field_for_visu(IJK_Field_double& temperature)
{
  const int ni = temperature.ni();
  const int nj = temperature.nj();
  const int nk = temperature.nk();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          // const double temperature = temperature_(i,j,k);
          const double indic = ref_ijk_ft_cut_cell_->itfce().I(i,j,k);
          // if (temperature > 0)
          if (indic < VAPOUR_INDICATOR_TEST)
            temperature(i,j,k) = 0;
        }
  temperature.echange_espace_virtuel(temperature.ghost());
}
