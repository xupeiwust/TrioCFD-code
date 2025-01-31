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

#include <Probleme_FTD_IJK_base.h>
#include <IJK_Field_vector.h>
#include <Probleme_FTD_IJK_cut_cell.h>
#include <Cut_cell_FT_Disc.h>
#include <IJK_FT_Post.h>
#include <IJK_Navier_Stokes_tools.h>
#include <EFichier.h>
#include <IJK_Splitting.h>
#include <Cut_cell_tools.h>
#include <Cut_cell_diffusion_flux_interface.h>

Implemente_instanciable(Probleme_FTD_IJK_cut_cell, "Probleme_FTD_IJK_cut_cell", Probleme_FTD_IJK_base);

Sortie& Probleme_FTD_IJK_cut_cell::printOn(Sortie& os) const
{
  return os;
}

Entree& Probleme_FTD_IJK_cut_cell::readOn(Entree& is)
{
  Probleme_FTD_IJK_base::readOn(is);

  // Determination pour le seuil des petites cellules en cut-cell
  if ((seuil_indicatrice_petite_fixe_ == -1) && (seuil_indicatrice_petite_facsec_ == -1))
    {
      seuil_indicatrice_petite_facsec_ = 0.125; // Default value = facsec/8. -- Note: The value facsec/20. = 0.01 would often work but not always be stable
    }

  double seuil_indicatrice_petite;
  if (seuil_indicatrice_petite_facsec_ != -1)
    {
      seuil_indicatrice_petite = timestep_facsec_*seuil_indicatrice_petite_facsec_;
    }
  else
    {
      seuil_indicatrice_petite = seuil_indicatrice_petite_fixe_;
    }
  interfaces_.set_seuil_indicatrice_petite(seuil_indicatrice_petite);
  Cerr << "Le seuil pour l'indicatrice des petites cellules est : " << seuil_indicatrice_petite << finl;

  return is;
}


void Probleme_FTD_IJK_cut_cell::set_param(Param& param)
{
  Probleme_FTD_IJK_base::set_param(param);

  param.ajouter("seuil_indicatrice_petite_fixe", &seuil_indicatrice_petite_fixe_);
  param.ajouter("seuil_indicatrice_petite_facsec", &seuil_indicatrice_petite_facsec_);

  param.ajouter("type_surface_efficace_face", (int*)&type_surface_efficace_face_);
  param.dictionnaire("non_initialise",(int)TYPE_SURFACE_EFFICACE_FACE::NON_INITIALISE);
  param.dictionnaire("explicite",(int)TYPE_SURFACE_EFFICACE_FACE::EXPLICITE);
  param.dictionnaire("algebrique_simple",(int)TYPE_SURFACE_EFFICACE_FACE::ALGEBRIQUE_SIMPLE);
  param.dictionnaire("conservation_volume_iteratif", (int)TYPE_SURFACE_EFFICACE_FACE::CONSERVATION_VOLUME_ITERATIF);
  param.ajouter("type_surface_efficace_interface", (int*)&type_surface_efficace_interface_);
  param.dictionnaire("non_initialise",(int)TYPE_SURFACE_EFFICACE_INTERFACE::NON_INITIALISE);
  param.dictionnaire("explicite",(int)TYPE_SURFACE_EFFICACE_INTERFACE::EXPLICITE);
  param.dictionnaire("algebrique_simple",(int)TYPE_SURFACE_EFFICACE_INTERFACE::ALGEBRIQUE_SIMPLE);
  param.dictionnaire("conservation_volume", (int)TYPE_SURFACE_EFFICACE_INTERFACE::CONSERVATION_VOLUME);

  param.ajouter("facettes_interpolation", &cut_cell_facettes_interpolation_);
}

bool Probleme_FTD_IJK_cut_cell::run()
{
  // Activation des champs cut-cell de post_ et interfaces_ (obligatoirement avant l'initialisation)
  cut_cell_disc_.initialise(interfaces_, splitting_, IJK_Splitting::ELEM);
  post_.activate_cut_cell();
  interfaces_.activate_cut_cell();
  cut_cell_facettes_interpolation_.associer(interfaces_, cut_cell_disc_, splitting_ft_, interfaces_.maillage_ft_ijk(), interfaces_.old_maillage_ft_ijk());

  splitting_.get_local_mesh_delta(DIRECTION_K, 2 /* ghost cells */,
                                  delta_z_local_);
  Cerr << "Probleme_FTD_IJK_cut_cell::run()" << finl;
  int nalloc = 0;
  thermal_probes_ghost_cells_ = 4;
  thermals_.compute_ghost_cell_numbers_for_subproblems(splitting_, thermal_probes_ghost_cells_);
  thermal_probes_ghost_cells_ = thermals_.get_probes_ghost_cells(thermal_probes_ghost_cells_);

  Cut_field_vector3_double& cut_field_velocity = static_cast<Cut_field_vector3_double&>(velocity_);
  if (IJK_Shear_Periodic_helpler::defilement_ == 1)
    {
      allocate_velocity_persistant(cut_cell_disc_, cut_field_velocity, splitting_, 2, boundary_conditions_.get_dU_perio(boundary_conditions_.get_resolution_u_prime_()));
    }
  else
    {
      allocate_velocity_persistant(cut_cell_disc_, cut_field_velocity, splitting_, thermal_probes_ghost_cells_);
    }

  if (IJK_Shear_Periodic_helpler::defilement_ == 1)
    {
      if (splitting_.get_nb_elem_local(2) < 4 )
        {
          std::cout << "Nb of cells / proc in z-direction must be >=4 for shear periodic run" << std::endl;
          std::cout << "Nb of cells / proc in z-direction must be >=8 for shear periodic run if only one proc on z-direction" << std::endl;
          std::cout << "Or find an other way to stock indic_ghost_zmin and zmax than IJK_Field" << std::endl;
          Process::exit();
        }
      if (splitting_.get_offset_local(0)!=0.)
        {
          std::cout << " Shear_periodic conditions works only without splitting in i-direction " << std::endl;
          std::cout << "if splitting in i-direction --> get_neighbour_processor has to be changed" << std::endl;
          Process::exit();
        }
    }

  allocate_velocity(d_velocity_, splitting_, 1);
  nalloc += 6;
  // GAB, qdm
  if (test_etapes_et_bilan_)
    {
      allocate_velocity(rho_u_euler_av_prediction_champ_, splitting_, 1);
      allocate_velocity(rho_u_euler_av_rho_mu_ind_champ_, splitting_, 1);
      allocate_velocity(rho_du_euler_ap_prediction_champ_, splitting_, 1);
      allocate_velocity(rho_u_euler_ap_projection_champ_, splitting_, 1);
      allocate_velocity(rho_du_euler_ap_projection_champ_, splitting_, 1);
      allocate_velocity(rho_u_euler_ap_rho_mu_ind_champ_, splitting_, 1);
      allocate_velocity(terme_diffusion_local_, splitting_, 1);
      allocate_velocity(terme_pression_local_, splitting_, 1);
      allocate_velocity(terme_pression_in_ustar_local_, splitting_, 1);
      allocate_velocity(d_v_diff_et_conv_, splitting_, 1);
      allocate_velocity(terme_convection_mass_solver_, splitting_, 1);
      allocate_velocity(terme_diffusion_mass_solver_, splitting_, 1);
      nalloc += 36;
    }
  //
  pressure_ghost_cells_.allocate(splitting_, IJK_Splitting::ELEM, thermal_probes_ghost_cells_);
  pressure_ghost_cells_.data() = 0.;
  pressure_ghost_cells_.echange_espace_virtuel(pressure_ghost_cells_.ghost());
  nalloc += 1;

  // if interp_monofluide == 2 --> reconstruction uniquement sur rho, mu. Pas sur P !
  if (!disable_diphasique_ && boundary_conditions_.get_correction_interp_monofluide()==1)
    {
      pressure_.allocate(splitting_, IJK_Splitting::ELEM, 3, 0 ,1, false, 1, rho_vapeur_, rho_liquide_, use_inv_rho_in_poisson_solver_);
    }
  else
    pressure_.allocate(splitting_, IJK_Splitting::ELEM, 3);
  nalloc += 1;

  if (include_pressure_gradient_in_ustar_)
    {
      d_pressure_.allocate(splitting_, IJK_Splitting::ELEM, 1);
      if (get_time_scheme() == RK3_FT)
        {
          RK3_F_pressure_.allocate(splitting_, IJK_Splitting::ELEM, 1);
          nalloc += 1;
        }
      nalloc += 1;
    }

  // On utilise aussi rhov pour le bilan de forces et pour d'autres formes de convection...
  //  if (!(expression_derivee_acceleration_ == Nom("0")))
  allocate_velocity(rho_v_, splitting_, 2);
  nalloc += 3;

  pressure_rhs_.allocate(splitting_, IJK_Splitting::ELEM, 1);
  nalloc += 1;
  I_ns_.allocate(splitting_, IJK_Splitting::ELEM, 2);
  kappa_ns_.allocate(splitting_, IJK_Splitting::ELEM, 2);
  if (!disable_diphasique_ && (boundary_conditions_.get_correction_interp_monofluide()==1 || boundary_conditions_.get_correction_interp_monofluide()==2))
    {
      molecular_mu_.allocate(splitting_, IJK_Splitting::ELEM, 2, 0 ,1, false, 2, mu_vapeur_, mu_liquide_);
      rho_field_.allocate(splitting_, IJK_Splitting::ELEM, 2, 0 ,1, false, 2, rho_vapeur_, rho_liquide_);
      nalloc += 2;
      IJK_Shear_Periodic_helpler::rho_vap_ref_for_poisson_=rho_vapeur_;
      IJK_Shear_Periodic_helpler::rho_liq_ref_for_poisson_=rho_liquide_;
      if (use_inv_rho_)
        {
          inv_rho_field_.allocate(splitting_, IJK_Splitting::ELEM, 2, 0 ,1, false, 2, 1./rho_vapeur_, 1./rho_liquide_);
          IJK_Shear_Periodic_helpler::rho_vap_ref_for_poisson_=1./rho_vapeur_;
          IJK_Shear_Periodic_helpler::rho_liq_ref_for_poisson_=1./rho_liquide_;
          nalloc += 1;
        }
    }
  else
    {

      IJK_Shear_Periodic_helpler::rho_vap_ref_for_poisson_=rho_vapeur_;
      IJK_Shear_Periodic_helpler::rho_liq_ref_for_poisson_=rho_liquide_;
      molecular_mu_.allocate(splitting_, IJK_Splitting::ELEM, 2);
      rho_field_.allocate(splitting_, IJK_Splitting::ELEM, 2);
      nalloc += 2;
      if (use_inv_rho_)
        {
          inv_rho_field_.allocate(splitting_, IJK_Splitting::ELEM, 2);
          nalloc += 1;
          IJK_Shear_Periodic_helpler::rho_vap_ref_for_poisson_=1./rho_vapeur_;
          IJK_Shear_Periodic_helpler::rho_liq_ref_for_poisson_=1./rho_liquide_;
        }
    }

  //  rho_batard_.allocate(splitting_, IJK_Splitting::ELEM, 2);

  if (first_step_interface_smoothing_)
    {
      allocate_velocity(zero_field_ft_, splitting_ft_, thermal_probes_ghost_cells_);
      for (int dir = 0; dir < 3; dir++)
        zero_field_ft_[dir].data() = 0.;
      zero_field_ft_.echange_espace_virtuel();
      nalloc += 3;
    }

  if (diffusion_alternative_)
    {
      allocate_velocity(laplacien_velocity_, splitting_, 1);
      unit_.allocate(splitting_, IJK_Splitting::ELEM, 2);
      unit_.data() = 1.;
      unit_.echange_espace_virtuel(unit_.ghost());
      nalloc += 4;
    }

  if (velocity_convection_op_.get_convection_op_option_rank() == non_conservative_rhou)
    {
      div_rhou_.allocate(splitting_, IJK_Splitting::ELEM, 1);
      nalloc += 1;
    }
  allocate_velocity(psi_velocity_, splitting_, 2);
  nalloc += 3;

#ifdef SMOOTHING_RHO
  // Pour le smoothing :
  rho_field_ft_.allocate(splitting_ft_, IJK_Splitting::ELEM, 2);
  nalloc += 1;
#endif
  // champs pour post-traitement :
  // post_.alloc_fields();
  nalloc += post_.alloc_fields();

  // Allocation du terme source variable spatialement:
  int flag_variable_source = false;
  if ((expression_variable_source_[0] != "??")
      || (expression_variable_source_[1] != "??")
      || (expression_variable_source_[2] != "??")
      || (expression_potential_phi_ != "??"))
    {
      allocate_velocity(variable_source_, splitting_, 1);
      flag_variable_source = true;
      potential_phi_.allocate(splitting_, IJK_Splitting::ELEM, 1);
      nalloc += 4;
      for (int dir = 0; dir < 3; dir++)
        variable_source_[dir].data() = 0.;
      potential_phi_.data() = 0.;
    }

  // thermals_.initialize(splitting_, nalloc);

  // GB : Je ne sais pas si on a besoin d'un ghost... Je crois que oui. Lequel?
  // Si la a vitesse ft doit transporter les sommets virtuels des facettes reelles,
  // alors il faut un domaine ghost de la taille de la longueur maximale des arretes.
  // allocate_velocity(velocity_ft_, splitting_ft_, 0);
  /*
   * FIXME: Allocate based on the thermal subproblems
   * as the thermal probes necessitates several ghost cells to interpolate velocity !
   * Check the difference between elem and faces ? and for interpolation of the velocity ?
   */
  /*
   * Finally use the ns velocity field for the thermal sub-problems
   */
  // thermals_.compute_ghost_cell_numbers_for_subproblems(splitting_, ft_ghost_cells);
  // ft_ghost_cells = thermals_.get_probes_ghost_cells(ft_ghost_cells);
  int ft_ghost_cells = 4;
  allocate_velocity(velocity_ft_, splitting_ft_, ft_ghost_cells);
  //  allocate_velocity(velocity_ft_, splitting_ft_, 4);
  nalloc += 3;

  kappa_ft_.allocate(splitting_ft_, IJK_Splitting::ELEM, 2);

  if (!disable_diphasique_)
    {
      allocate_velocity(terme_source_interfaces_ft_, splitting_ft_, 2);
      allocate_velocity(backup_terme_source_interfaces_ft_, splitting_, 2);
      // Seulement pour le calcul du bilan de forces :
      allocate_velocity(terme_source_interfaces_ns_, splitting_, 1);
      allocate_velocity(backup_terme_source_interfaces_ns_, splitting_, 1);
      // Seulement pour le calcul des statistiques :
      allocate_velocity(terme_repulsion_interfaces_ns_, splitting_, 1);
      allocate_velocity(terme_repulsion_interfaces_ft_, splitting_ft_, 1);
      allocate_velocity(terme_abs_repulsion_interfaces_ns_, splitting_, 1);
      allocate_velocity(terme_abs_repulsion_interfaces_ft_, splitting_ft_, 1);
      nalloc += 18;
    }

  // FIXME: on a oublie pleins de choses la !
  // int nalloc = 24;
  nalloc += post_.alloc_velocity_and_co(flag_variable_source);
  if (get_time_scheme() == RK3_FT)
    {
      allocate_velocity(RK3_F_velocity_, splitting_, 1);
      nalloc += 3;
      Cout << "Schema temps de type : RK3_FT" << finl;
    }
  else
    Cout << "Schema temps de type : euler_explicite" << finl;

  velocity_diffusion_op_.initialize(splitting_, harmonic_nu_in_diff_operator_);
  velocity_diffusion_op_->set_bc(boundary_conditions_);
  velocity_convection_op_.initialize(splitting_);

  treatment_count_.allocate(splitting_, IJK_Splitting::ELEM, 2);
  nalloc += 1;


  // Economise la memoire si pas besoin
  if (!disable_solveur_poisson_)
    poisson_solver_.initialize(splitting_);


  nalloc += initialise_interfaces();

  // C'est ici aussi qu'on alloue les champs de temperature.
  nalloc += initialise();

//  rho_field_.echange_espace_virtuel(2);
//  recalculer_rho_de_chi(chi_, rho_field_, 2);
  Cerr << " Allocating " << nalloc << " arrays, approx total size= "
       << (double)(molecular_mu_.data().size_array() * (int)sizeof(double) * nalloc)
       * 9.537E-07 << " MB per core" << finl;

// Les champs ont etes alloues.
// On peut completer les sondes car les ijk_field.get_splitting() sont a present remplis.
  post_.completer_sondes();
  post_.improved_initial_pressure_guess(improved_initial_pressure_guess_);

// Cette projection n'est pas utile en reprise.
// Elle sert uniquement a rendre le champ de vitesse initial a divergence nulle
// lorsque son expression est analytique.

  if (!disable_solveur_poisson_)
    {
      if (improved_initial_pressure_guess_)
        {
          Cerr << "Improved initial pressure" << finl;
          maj_indicatrice_rho_mu();
          if (!disable_diphasique_)
            {
              /*
               * TODO: Change this block with OWN_PTR CLASS IJK_Thermal
               */
              for (auto& itr : thermique_)
                itr.update_thermal_properties();

              for (auto& itr : energie_)
                itr.update_thermal_properties();

              thermals_.update_thermal_properties();
            }
          // La pression n'est pas encore initialisee. elle est donc nulle.
          // Avec cette option, on essaye une initialisation basee sur le champ de pression diphasique
          // a l'equilibre, cad sans vitesse, ou a minima pour un champ a div(u)=0.

          if (!disable_diphasique_)
            {
              IJK_Field_vector3_double& coords = post_.coords();

              // Calcul du potentiel.
              for (int dir = 0; dir < 3; dir++)
                {
                  terme_source_interfaces_ft_[dir].data() = 0.;
                  terme_repulsion_interfaces_ft_[dir].data() = 0.;
                  terme_abs_repulsion_interfaces_ft_[dir].data() = 0.;
                }
              const double delta_rho = rho_liquide_ - rho_vapeur_;
              interfaces_.ajouter_terme_source_interfaces(
                terme_source_interfaces_ft_,
                terme_repulsion_interfaces_ft_,
                terme_abs_repulsion_interfaces_ft_
              );

              assert(interfaces_.get_nb_bulles_reelles() == 1);
              DoubleTab bounding_box;
              interfaces_.calculer_bounding_box_bulles(bounding_box);
              // Calcul la hauteur en x de la permiere bulle :
              const double Dbx = bounding_box(0, 0, 1)
                                 - bounding_box(0, 0, 0);
              const double kappa = 2. / (Dbx / 2.);

              const int ni = pressure_.ni();
              const int nj = pressure_.nj();
              const int nk = pressure_.nk();
              for (int k = 0; k < nk; k++)
                for (int j = 0; j < nj; j++)
                  for (int i = 0; i < ni; i++)
                    {
                      double phi = gravite_[0] * coords[0](i, j, k)
                                   + gravite_[1] * coords[1](i, j, k)
                                   + gravite_[2] * coords[2](i, j, k);
                      double potentiel_elem = sigma_ * kappa
                                              - delta_rho * phi;
                      // La pression est hydrostatique, cad : pressure_ = P - rho g z
                      pressure_(i, j, k) = potentiel_elem
                                           * interfaces_.I(i, j, k); // - rho_field_(i,j,k) * phi;
                    }

              // pressure gradient requires the "left" value in all directions:
              pressure_.echange_espace_virtuel(
                1 /*, IJK_Field_double::EXCHANGE_GET_AT_LEFT_IJK*/);

              // Mise a jour du champ de vitesse (avec dv = seulement le terme source)
              d_velocity_[0].data() = 0.;
              d_velocity_[1].data() = 0.;
              d_velocity_[2].data() = 0.;

              for (int dir = 0; dir < 3; dir++)
                redistribute_from_splitting_ft_faces_[dir].redistribute_add(
                  terme_source_interfaces_ft_[dir], d_velocity_[dir]);

              for (int dir = 0; dir < 3; dir++)
                {
                  const int kmax = d_velocity_[dir].nk();
                  for (int k = 0; k < kmax; k++)
                    {
                      euler_explicit_update(d_velocity_[dir], velocity_[dir],
                                            k);
                    }
                }

              if (use_inv_rho_in_poisson_solver_)
                {
                  pressure_projection_with_inv_rho(inv_rho_field_,
                                                   velocity_[0], velocity_[1], velocity_[2], pressure_,
                                                   1., pressure_rhs_, check_divergence_,
                                                   poisson_solver_);

                }
              else
                {



                  pressure_projection_with_rho(rho_field_, velocity_[0],
                                               velocity_[1], velocity_[2], pressure_, 1.,
                                               pressure_rhs_, check_divergence_, poisson_solver_);

                }

            }
          else
            {

              pressure_projection(velocity_[0], velocity_[1], velocity_[2],
                                  pressure_, 1., pressure_rhs_, check_divergence_,
                                  poisson_solver_);
            }
          copy_field_values(pressure_ghost_cells_, pressure_);
        }
    }

  const double max_timestep = timestep_;

// Si calcul monophasique, on initialise correctement rho, mu, I une fois pour toute :
  if (disable_diphasique_)
    {
      rho_field_.data() = rho_liquide_;
      rho_moyen_ = rho_liquide_;
      molecular_mu_.data() = mu_liquide_;

      // C'est deja fait dans l'initialize (aucune raison de ne pas le faire)
      // indicatrice_ns_.data() = 1.;
      // indicatrice_ns_next_.data() = 1.;

      /*
       * TODO: Change this block with OWN_PTR CLASS IJK_Thermal
       */
      for (auto& itr : thermique_)
        {
          // To fill in fields for cp (with cp_liq) and lambda (with lambda_liq)
          itr.update_thermal_properties();
        }
      for (auto& itr : energie_)
        {
          // To fill in fields for cp (with cp_liq) and lambda (with lambda_liq)
          itr.update_thermal_properties();
        }
      thermals_.update_thermal_properties();
    }
  else
    {
      Cerr << "Cas normal diphasique Probleme_FTD_IJK_cut_cell::run()" << finl;

      /*
       * TODO: Change this block with OWN_PTR CLASS IJK_Thermal
       */
      for (auto& itr : thermique_)
        itr.update_thermal_properties();

      for (auto& itr : energie_)
        itr.update_thermal_properties();

      thermals_.update_thermal_properties();

      const double indic_moyen = calculer_v_moyen(interfaces_.I());
      rho_moyen_ = indic_moyen*rho_liquide_ + (1-indic_moyen)*rho_vapeur_;
      if (post_.get_liste_post_instantanes().contient_("EXTERNAL_FORCE"))
        {
          for (int dir=0; dir<3; dir++)
            compute_add_external_forces(dir);
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
  if (!disable_solveur_poisson_)
    {
      if (projection_initiale_demandee_)
        {
          Cerr << "*****************************************************************************\n"
               << "  Attention : projection du champ de vitesse initial sur div(u)=0\n"
               << "*****************************************************************************" << finl;

          pressure_projection_with_rho(rho_field_, velocity_[0],
                                       velocity_[1], velocity_[2], pressure_, 1.,
                                       pressure_rhs_, check_divergence_, poisson_solver_);
          pressure_.data() = 0.;
          pressure_rhs_.data() = 0.;
        }
    }


  if ((!disable_diphasique_) && (post_.get_liste_post_instantanes().contient_("VI")
                                 || post_.get_liste_post_instantanes().contient_("TOUS")))
    interfaces_.compute_vinterp();

  // Preparer le fichier de postraitement et postraiter la condition initiale:
  Nom lata_name = nom_du_cas();
  if (fichier_post_ != "??")
    {
      lata_name = fichier_post_;
    }
  lata_name += Nom(".lata");
  post_.postraiter_ci(lata_name, current_time_);

//  if ( 0 && disable_diphasique_
//       && (liste_post_instantanes_.contient_("CURL")))
//    {
//      Cerr << " Dans un calcul monophasique, on ne calcule pas le rotationnel, "
//           << "donc ce n'est pas la peine de le demander dans les posts!" << finl;
//      Process::exit();
//    }

  post_.compute_extended_pressures(interfaces_.maillage_ft_ijk());
//post_.compute_phase_pressures_based_on_poisson(0);
//post_.compute_phase_pressures_based_on_poisson(1);

  modified_time_ini_ = thermals_.get_modified_time();
  if (!reprise_ && current_time_ == 0.)
    current_time_ = modified_time_ini_;

  if (!first_step_interface_smoothing_)
    {
      Cout << "BF posttraiter_champs_instantanes "
           << current_time_ << " " << tstep_ << finl;
      post_.posttraiter_champs_instantanes(lata_name, current_time_, tstep_);
      thermals_.thermal_subresolution_outputs(); // for thermal counters
      Cout << "AF posttraiter_champs_instantanes" << finl;
    }

// GB 2019.01.01 Why immobilisation? if (!disable_diphasique_ && coef_immobilisation_==0.)
  if ((!disable_diphasique_) && suppression_rejetons_)
    interfaces_.detecter_et_supprimer_rejeton(true);
  if (reprise_)
    {
      // On ecrit a la suite du fichier. Cela suppose qu'il est bien a jour.
      // L'instant initial a deja ete ecrit a la fin du calcul precedent donc on
      // ne le reecrit pas.
    }
  else
    {
      // On creer de nouveaux fichiers :
      Cout << "BF ecrire_statistiques_bulles" << finl;
      post_.ecrire_statistiques_bulles(1 /* reset files */, nom_du_cas(),
                                       gravite_, current_time_);
      Cout << "AF ecrire_statistiques_bulles" << finl;
    }

// Ecrire la valeur initiale dans les sondes :
// Ecriture de la valeur initiale seulement hors reprise
  if (!reprise_)
    post_.postraiter_sondes();

//ab-forcage-control-ecoulement-deb
  update_rho_v(); // Peut-etre pas toujours necessaire selon la formulation pour la convection?
  for (int direction = 0; direction < 3; direction++)
    store_rhov_moy_[direction] = calculer_v_moyen(rho_v_[direction]);
//ab-forcage-control-ecoulement-fin

  statistiques().end_count(initialisation_calcul_counter_);

  if (!disable_TU)
    {
      if(GET_COMM_DETAILS)
        statistiques().print_communciation_tracking_details("Statistiques d'initialisation du calcul", 0);

      statistiques().dump("Statistiques d'initialisation du calcul", 0);
      print_statistics_analyse("Statistiques d'initialisation du calcul", 0);
    }
  statistiques().reset_counters();
  statistiques().begin_count(temps_total_execution_counter_);

  int stop = 0;

// Variation de volume de chaque bulle integree au cours du pas de temps :
  ArrOfDouble var_volume_par_bulle;
  for (tstep_ = 0; tstep_ < nb_timesteps_ && stop == 0; tstep_++)
    {
      statistiques().begin_count(timestep_counter_);
      if (!splitting_.get_grid_geometry().get_periodic_flag(DIRECTION_K))
        {
          force_zero_on_walls(velocity_[2]);
        }

      if (timestep_facsec_ > 0.)
        {
          double max_post_simu_timestep = post_.get_timestep_simu_post(current_time_, max_simu_time_);
          timestep_ = find_timestep(std::min(max_timestep, max_post_simu_timestep), cfl_, fo_, oh_);
        }

      // Tableau permettant de calculer la variation de volume au cours du pas de temps :
      // Si on veut le mettre en optionel, il faut faire attention a faire vivre la taille de ce tableau avec les
      // creations et destructions de ghosts :
      const int nbulles_tot = interfaces_.get_nb_bulles_reelles()
                              + interfaces_.get_nb_bulles_ghost(1/*print=1*/);
      var_volume_par_bulle.resize_array(nbulles_tot);
      var_volume_par_bulle = 0.; // Je ne suis pas sur que ce soit un bon choix. Si on ne le remet pas a zero
      //                          a chaque dt, on corrigera la petite erreur qui pouvait rester d'avant...

      compute_var_volume_par_bulle(var_volume_par_bulle);

      // Au cas ou on soit dans un cas ou des duplicatas sont necessaires mais n'ont pas ete
      // crees, on les cree :
      if (!interfaces_.get_nb_bulles_ghost() && !disable_diphasique_)
        {
          interfaces_.creer_duplicata_bulles();
        }

      // Choix de l'avancement en temps :
      // euler_explicite ou RK3.
      if (get_time_scheme() == EULER_EXPLICITE)
        {
          // Deplacement des interfaces par le champ de vitesse de l'instant n :
          if (!disable_diphasique_)
            {
              // TODO: aym pour GAB, si tu veux gagner en memoire et virer le doublon n/np1 il faut
              // inserer une methode ici style "mettre_a_jour_valeur_interface_temps_n()"
              int counter_first_iter = 1;
              do
                {
                  first_step_interface_smoothing_ = first_step_interface_smoothing_ && counter_first_iter;
                  deplacer_interfaces(timestep_,
                                      -1 /* le numero du sous pas de temps est -1 si on n'est pas en rk3 */,
                                      var_volume_par_bulle,
                                      first_step_interface_smoothing_);
                  counter_first_iter--;
                  if(first_step_interface_smoothing_)
                    {
                      thermals_.set_temperature_ini();
                      Cout << "BF posttraiter_champs_instantanes " << current_time_ << " " << tstep_ << finl;
                      post_.posttraiter_champs_instantanes(lata_name, current_time_, tstep_);
                      Cout << "AF posttraiter_champs_instantanes" << finl;
                      compute_var_volume_par_bulle(var_volume_par_bulle);
                      thermals_.set_post_pro_first_call();
                    }
                }
              while (first_step_interface_smoothing_);
              parcourir_maillage();
            }
          // Mise a jour de la vitesse (utilise les positions des marqueurs, rho, mu et indic a l'instant n)
          // Retourne une vitesse mise a jour et projetee a div nulle
          euler_time_step(var_volume_par_bulle);

          // Calcul du terme source force acceleration :
          // GAB : question a Guillaume, on fait time + time_step ? 'est pas homogene non ?'
          // GAB, rotation
          // /!\ On  laisse ce calcul active meme pour source_qdm_gr_!=-1 pour toujours avoir un fichier acceleration.out rempli correctement
          calculer_terme_source_acceleration(velocity_[direction_gravite_],
                                             current_time_ + timestep_, timestep_, -1);

          // Deplacement des interfaces par le champ de vitesse :
          // met a jour la position des marqueurs, la vitesse_ft, et gere les duplicatas.
          // Ne met pas a jour rho_mu_indicatrice

          if (!disable_diphasique_) // && !marker_advection_first_)
            {
              // Les sous-pas de temps sont termines. Il n'est plus necessaire de gerer le tableau
              // RK3_G_store_vi_. On peut donc transferer les bulles et re-creer les duplicatas :
              interfaces_.supprimer_duplicata_bulles();
              interfaces_.transferer_bulle_perio();
              // On supprime les fragments de bulles.
              //interfaces_.detecter_et_supprimer_rejeton(false);
              interfaces_.creer_duplicata_bulles();

              // indicatrice (and rho, mu...) are updated from the new interface position.
              // GB 2019.01.01 It is important to keep that calculation, because without it, the interface status would be
              // set to "minimal" where it should be "parcouru".
              // GAB, qdm : rho_n v_n+1
              if (test_etapes_et_bilan_)
                {
                  calculer_rho_v(rho_field_,velocity_,rho_u_euler_av_rho_mu_ind_champ_);
                  for (int dir=0; dir<3; dir++)
                    rho_u_euler_av_rho_mu_ind_[dir] = calculer_v_moyen(rho_u_euler_av_rho_mu_ind_champ_[dir]);
                }
              maj_indicatrice_rho_mu();

              /*
               * TODO: Change this block with OWN_PTR CLASS IJK_Thermal
               */
              for (auto& itr : thermique_)
                {
                  itr.update_thermal_properties();
                  if (itr.conserv_energy_global_)
                    {
                      const double dE = itr.E0_ - itr.compute_global_energy();
                      itr.euler_rustine_step(timestep_, dE);
                    }
                }

              thermals_.euler_rustine_step(timestep_);

              // GAB, qdm rho_n+1 v_n+1 :
              if (test_etapes_et_bilan_)
                {
                  calculer_rho_v(rho_field_,velocity_,rho_u_euler_ap_rho_mu_ind_champ_);
                  for (int dir=0; dir<3; dir++)
                    {
                      rho_u_euler_ap_rho_mu_ind_[dir] = calculer_v_moyen(rho_u_euler_ap_rho_mu_ind_champ_[dir]);
                      u_euler_ap_rho_mu_ind_[dir] = calculer_v_moyen(velocity_[dir]);
                    }
                }
            }
          else
            {
              if (test_etapes_et_bilan_)
                {
                  calculer_rho_v(rho_field_,velocity_,rho_u_euler_ap_rho_mu_ind_champ_);
                  for (int dir=0; dir<3; dir++)
                    {
                      rho_u_euler_ap_rho_mu_ind_[dir] = calculer_v_moyen(rho_u_euler_ap_rho_mu_ind_champ_[dir]);
                      u_euler_ap_rho_mu_ind_[dir] = calculer_v_moyen(velocity_[dir]);
                    }
                }
            }
        }
      else if (get_time_scheme() == RK3_FT)
        {
          double current_time_at_rk3_step = current_time_;
          // GAB, qdm : passe en attribut de classe car utilise au moment de l'ecriture de mon out
          current_time_at_rk3_step_ = current_time_;
          // Evaluation de la variation de volume accumule au cours des sous pas de temps.
          // On la laisse croitre pendant les sous dt 0 et 1 puis on la corrige a la fin du 2eme :

          for (rk_step_ = 0; rk_step_ < 3; rk_step_++)
            {
              const double fractionnal_timestep =
                compute_fractionnal_timestep_rk3(timestep_ /* total*/,
                                                 rk_step_);

              // Mise a jour des positions des marqueurs.
              // Deplacement des interfaces par le champ de vitesse au sous pas de temps k :
              if (!disable_diphasique_)
                {
                  deplacer_interfaces_rk3(timestep_ /* total */, rk_step_,
                                          var_volume_par_bulle);
                  parcourir_maillage();
                }
              // Cerr << "RK3 : step " << rk_step << finl;
              // Mise a jour de la temperature et de la vitesse :
              rk3_sub_step(rk_step_, timestep_, fractionnal_timestep,
                           current_time_at_rk3_step);

              // GAB patch qdm : choix 1

              // GAB, qdm : rho_n v_n+1
              if (test_etapes_et_bilan_)
                {
                  calculer_rho_v(rho_field_,velocity_,rho_u_euler_av_rho_mu_ind_champ_);
                  for (int dir=0; dir<3; dir++)
                    rho_u_euler_av_rho_mu_ind_[dir] = calculer_v_moyen(rho_u_euler_av_rho_mu_ind_champ_[dir]);
                }

              // Mise a jour rho, mu et l'indicatrice a partir de la nouvelle position de l'interface :
              // (sauf au dernier sous pas de temps pour lequel c'est fait a la fin du pas de temps)
              // TODO: verifier qu'on doit bien le faire aussi au dernier sous pas de temps : rk_step != 2 &&
              // TODO aym: verifier ce bloc, qui applique les sous pas de temps RK3 de la rustine a la temperature
              if (rk_step_ != 2 && !disable_diphasique_)
                {
                  // Attention, il faut que les duplicatas soient present pour faire maj_indicatrice_rho_mu :
                  maj_indicatrice_rho_mu();
                  for (auto& itr : thermique_)
                    {
                      itr.update_thermal_properties();
                      if (itr.conserv_energy_global_)
                        {
                          const double dE = itr.E0_ - itr.compute_global_energy();
                          itr.rk3_rustine_sub_step(rk_step_, timestep_, fractionnal_timestep,
                                                   current_time_at_rk3_step, dE);
                        }
                    }

                  thermals_.rk3_rustine_sub_step(rk_step_, timestep_, fractionnal_timestep,
                                                 current_time_at_rk3_step);

                }
              // Calcul du terme source force acceleration :
              // GAB, rotation
              calculer_terme_source_acceleration(velocity_[direction_gravite_],
                                                 current_time_at_rk3_step, timestep_ /*total*/, rk_step_);


              current_time_at_rk3_step += fractionnal_timestep;
              // GAB, qdm : passe en attribut de classe car utilise au moment de l'ecriture de mon out
              current_time_at_rk3_step_ += fractionnal_timestep;
              // On ne postraite pas le sous-dt 2 car c'est fait plus bas si on post-traite le pas de temps :
              if (post_.postraiter_sous_pas_de_temps()
                  && ((tstep_ % post_.dt_post() == post_.dt_post() - 1)
                      || (std::floor((current_time_-timestep_)/post_.get_timestep_simu_post(current_time_, max_simu_time_)) < std::floor(current_time_/post_.get_timestep_simu_post(current_time_, max_simu_time_))))
                  && (rk_step_ != 2))
                {
                  cut_field_velocity[0].copie_pure_vers_diph_sans_interpolation();
                  cut_field_velocity[1].copie_pure_vers_diph_sans_interpolation();
                  cut_field_velocity[2].copie_pure_vers_diph_sans_interpolation();

                  post_.posttraiter_champs_instantanes(lata_name, current_time_at_rk3_step, tstep_);
                }
            }
          if (!disable_diphasique_)
            {
              // Les sous-pas de temps sont termines. Il n'est plus necessaire de gerer le tableau
              // RK3_G_store_vi_. On peut donc transferer les bulles et re-creer les duplicatas :
              interfaces_.supprimer_duplicata_bulles();
              interfaces_.transferer_bulle_perio();
              // On supprime les fragments de bulles.
              //interfaces_.detecter_et_supprimer_rejeton(false);
              interfaces_.creer_duplicata_bulles();

              // Mise a jour rho, mu et l'indicatrice a partir de la nouvelle position de l'interface :
              maj_indicatrice_rho_mu();

              for (auto& itr : thermique_)
                itr.update_thermal_properties();

              for (auto& itr : energie_)
                itr.update_thermal_properties();

              thermals_.update_thermal_properties();
            }
          // GAB, qdm rho_n+1 v_n+1 :
          if (test_etapes_et_bilan_)
            {
              calculer_rho_v(rho_field_,velocity_,rho_u_euler_ap_rho_mu_ind_champ_);
              for (int dir=0; dir<3; dir++)
                {
                  rho_u_euler_ap_rho_mu_ind_[dir] = calculer_v_moyen(rho_u_euler_ap_rho_mu_ind_champ_[dir]);
                  u_euler_ap_rho_mu_ind_[dir] = calculer_v_moyen(velocity_[dir]);
                }
            }
        }
      else
        {
          Cerr << "Erreur dans le run: time_scheme " << time_scheme_
               << " inconnu!" << finl;
          Process::exit();
        }
      // ------------------------------------------------------------------
      // CORRECTION DE QUANTITE DE MOUVEMENT
      // ------------------------------------------------------------------
      // Correction de QdM : permet de controler la QdM globale, dans chaque direction
      if (!(qdm_corrections_.is_type_none()) )
        {
          set_time_for_corrections();
          if (disable_diphasique_)
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
        {write_qdm_corrections_information();}
      else
        {;}
      // ------------------------------------------------------------------

      //ab-forcage-control-ecoulement-deb
      // Quel que soit le schema en temps, on corrige le bilan de qdm par le residu integre :
      // integrated_residu_ est homogene a rho*u.
      // Il faut donc appliquer le solveur masse a integrated_residu_. Pour cela, on a besoin d'un champ.
      //                              On prend psi_velocity_ qui est dispo.
      // Attention, en entree du solveur mass, il faut qqch homogene a rho*u*volume_cell...
      // On rempli donc psi_velocity avec vol * integrated_residu_

      // static Stat_Counter_Id bilanQdM_counter_ = statistiques().new_counter(2, "Bilan QdM & Corrections");
      // statistiques().begin_count(bilanQdM_counter_);
      // statistiques().end_count(bilanQdM_counter_);

      //ab-forcage-control-ecoulement-fin
      current_time_ += timestep_;
      // stock dans le spliting le decallage periodique total avec condition de shear (current_time_) et celui du pas de temps (timestep_)
      IJK_Shear_Periodic_helpler::shear_x_time_ = boundary_conditions_.get_dU_perio()*(current_time_ + boundary_conditions_.get_t0_shear());

      if (current_time_ >= post_.t_debut_statistiques())
        {
          if (boundary_conditions_.get_correction_conserv_qdm()==2)
            {
              update_rho_v();
              rho_field_.echange_espace_virtuel(rho_field_.ghost());
              update_v_ghost_from_rho_v();
            }
          else
            {
              // FA AT 16/07/2013 pensent que necessaire pour le calcul des derivees dans statistiques_.update_stat_k(...)
              // Je ne sais pas si c'est utile, mais j'assure...
              velocity_[0].echange_espace_virtuel(
                2 /*, IJK_Field_ST::EXCHANGE_GET_AT_RIGHT_I*/);
              velocity_[1].echange_espace_virtuel(
                2 /*, IJK_Field_ST::EXCHANGE_GET_AT_RIGHT_J*/);
              velocity_[2].echange_espace_virtuel(
                2 /*, IJK_Field_ST::EXCHANGE_GET_AT_RIGHT_K*/);
            }


          pressure_.echange_espace_virtuel(1);

          // C'est update_stat_ft qui gere s'il y a plusieurs groupes
          // pour faire la vraie indicatrice + les groupes
          post_.update_stat_ft(timestep_);
          if (!disable_diphasique_)
            {
              post_.compute_extended_pressures(interfaces_.maillage_ft_ijk());
              //post_.compute_phase_pressures_based_on_poisson(0);
              //post_.compute_phase_pressures_based_on_poisson(1);
            }
        }

      // Calcul du terme source d'acceleration deplacee dans les iterations du rk3 ou dans l'iteration d'euler.

      // verification du fichier stop
      stop = 0;
      if (check_stop_file_ != "??")
        {
          if (je_suis_maitre())
            {
              EFichier f;
              stop = f.ouvrir(check_stop_file_);
              if (stop)
                {
                  // file exists, check if it contains 1:
                  f >> stop;
                }
            }
          envoyer_broadcast(stop, 0);
        }
      if (tstep_ == nb_timesteps_ - 1)
        stop = 1;
      if (current_time_ >= max_simu_time_)
        stop = 1;

      tstep_sauv_ = tstep_ + tstep_init_;
      if (tstep_sauv_ % dt_sauvegarde_ == dt_sauvegarde_ - 1 || stop)
        {
          // Choix : On supprime les duplicatas pour la sauvegarde.
          // On pourrait tres bien tout garder. ca serait plus leger en CPU, plus lourd en espace disque.
          if (!disable_diphasique_)
            interfaces_.supprimer_duplicata_bulles();

          sauvegarder_probleme(nom_sauvegarde_, stop);
          if (!disable_diphasique_)
            {
              // On les recree :
              interfaces_.creer_duplicata_bulles();

              // Be on the safe side, on met a jour :
              //   A la suppression des duplicatas, on avait fait mesh.supprimer_facettes qui remet le maillage
              //   a l'etat MINIMAL. Pour les post-tt sur l'interface (eg ai_ft_), il faut que le statut du maillage
              //   soit >= PARCOURU. C'est fait au debut de maj_indicatrice_rho_mu dans
              //   IJK_Interfaces::calculer_indicatrice.
              const double delta_rho = rho_liquide_ - rho_vapeur_;
              interfaces_.calculer_indicatrice_next(
                post_.potentiel(),
                gravite_,
                delta_rho,
                sigma_,
                /*Pour post-traitement : post_.rebuilt_indic()
                */
#ifdef SMOOTHING_RHO
                /* Pour le smoothing : */
                rho_field_ft_,
                rho_vapeur_,
                smooth_density_,
#endif
                current_time_, tstep_
              );
            }
        }

      // TODO: on pourrait mutualiser tous les parcourir maillages dans IJK_Interface au moment du transport de l'interface
      // et le supprimer de Probleme_FTD_IJK_cut_cell
      // interfaces_.parcourir_maillage();
      if ((!disable_diphasique_) && (post_.get_liste_post_instantanes().contient_("VI")))
        interfaces_.compute_vinterp();

      cut_field_velocity[0].copie_pure_vers_diph_sans_interpolation();
      cut_field_velocity[1].copie_pure_vers_diph_sans_interpolation();
      cut_field_velocity[2].copie_pure_vers_diph_sans_interpolation();

      post_.postraiter_fin(stop, tstep_, tstep_init_, current_time_, timestep_, lata_name,
                           gravite_, nom_du_cas());
      statistiques().end_count(timestep_counter_);

      if(JUMP_3_FIRST_STEPS && tstep_ < 3)
        {
          //demarrage des compteurs CPU
          if(tstep_ == 2)
            statistiques().set_three_first_steps_elapsed(true);
        }
      else
        statistiques().compute_avg_min_max_var_per_step(tstep_);


    }
  if (Process::je_suis_maitre())
    {
      SFichier master_file;
      master_file.ouvrir(lata_name, ios::app);
      master_file << "FIN" << finl;
      master_file.close();
    }
// Pour forcer l'ecriture du dernier pas de temps dans la sonde (peut-etre deja ecrit...)
// Alan 2020/03/02 : effectivement, deja ecrit
// post_.postraiter_sondes();

  if (!disable_TU)
    {
      if(GET_COMM_DETAILS)
        statistiques().print_communciation_tracking_details("Statistiques de resolution du probleme", 1);

      statistiques().dump("Statistiques de resolution du probleme", 1);
      print_statistics_analyse("Statistiques de resolution du probleme", 1);
    }

  statistiques().reset_counters();
  statistiques().begin_count(temps_total_execution_counter_);

  return true;
}

void Probleme_FTD_IJK_cut_cell::update_indicator_field()
{
  // La suppression des cellules mortes est vraiment au tout dernier moment,
  // pour laisser la possibilite d'utiliser ces cellules lors des bilans
  // pour determiner l'indicatrice cible du remaillage.
  cut_cell_disc_.remove_dead_and_virtual_cells(interfaces_.In());

  Probleme_FTD_IJK_base::update_indicator_field();
}

void Probleme_FTD_IJK_cut_cell::update_twice_indicator_field()
{
  for(int i=0; i<2; i++)
    {
      update_indicator_field();
      update_old_intersections();
    }

  // Mise a jour des structures cut-cell
  cut_cell_disc_.update(interfaces_.I(), interfaces_.In());

  // Mise a jour des indices et coefficients des points d'interpolation a une certaine distance des facettes de l'interface
  cut_cell_perform_interpolation_facettes();

  // Calcul pour le temps old() egalement, de telle maniere a ce que les coefficients next() et old() sont initialises
  cut_cell_facettes_interpolation_.cut_cell_perform_interpolation_facettes_old(interfaces_.old());

  Cut_field_vector3_double& cut_field_velocity = static_cast<Cut_field_vector3_double&>(velocity_);
  cut_field_velocity[0].copie_pure_vers_diph_sans_interpolation();
  cut_field_velocity[1].copie_pure_vers_diph_sans_interpolation();
  cut_field_velocity[2].copie_pure_vers_diph_sans_interpolation();
}

void Probleme_FTD_IJK_cut_cell::deplacer_interfaces(const double timestep, const int rk_step,
                                                    ArrOfDouble& var_volume_par_bulle,
                                                    const int first_step_interface_smoothing)
{
  Cut_field_vector3_double& cut_field_velocity = static_cast<Cut_field_vector3_double&>(velocity_);

  thermals_.echange_diph_vers_pure_cellules_finalement_pures();
  thermals_.vide_phase_invalide_cellules_diphasiques();
  update_old_intersections(); // Pour conserver les donnees sur l'interface au temps t_{n} (en plus de t_{n+1})

  cut_field_velocity[0].copie_pure_vers_diph_sans_interpolation();
  cut_field_velocity[1].copie_pure_vers_diph_sans_interpolation();
  cut_field_velocity[2].copie_pure_vers_diph_sans_interpolation();

  Probleme_FTD_IJK_base::deplacer_interfaces(timestep, rk_step, var_volume_par_bulle, first_step_interface_smoothing);

  // Mise a jour des structures cut-cell
  cut_cell_disc_.update(interfaces_.I(), interfaces_.In());

  // Mise a jour des indices et coefficients des points d'interpolation a une certaine distance des facettes de l'interface
  cut_cell_perform_interpolation_facettes();

  thermals_.echange_pure_vers_diph_cellules_initialement_pures();
  cut_field_velocity[0].copie_pure_vers_diph_sans_interpolation();
  cut_field_velocity[1].copie_pure_vers_diph_sans_interpolation();
  cut_field_velocity[2].copie_pure_vers_diph_sans_interpolation();

  interfaces_.calcul_surface_efficace_face_initial(type_surface_efficace_face_);
  interfaces_.calcul_surface_efficace_interface_initial(type_surface_efficace_interface_);

  interfaces_.calcul_surface_efficace_face(type_surface_efficace_face_, timestep_, cut_field_velocity);
  interfaces_.calcul_surface_efficace_interface(type_surface_efficace_interface_, timestep_, cut_field_velocity);

  if (interfaces_.get_dt_impression_bilan_indicatrice() >= 0 && tstep_ % interfaces_.get_dt_impression_bilan_indicatrice() == interfaces_.get_dt_impression_bilan_indicatrice() - 1)
    {
      interfaces_.imprime_bilan_indicatrice();
    }
}

void Probleme_FTD_IJK_cut_cell::deplacer_interfaces_rk3(const double timestep, const int rk_step,
                                                        ArrOfDouble& var_volume_par_bulle)
{
  Cut_field_vector3_double& cut_field_velocity = static_cast<Cut_field_vector3_double&>(velocity_);

  thermals_.echange_diph_vers_pure_cellules_finalement_pures();
  thermals_.vide_phase_invalide_cellules_diphasiques();
  update_old_intersections(); // Pour conserver les donnees sur l'interface au temps t_{n} (en plus de t_{n+1})

  cut_field_velocity[0].copie_pure_vers_diph_sans_interpolation();
  cut_field_velocity[1].copie_pure_vers_diph_sans_interpolation();
  cut_field_velocity[2].copie_pure_vers_diph_sans_interpolation();

  Probleme_FTD_IJK_base::deplacer_interfaces_rk3(timestep, rk_step, var_volume_par_bulle);

  // Mise a jour des structures cut-cell
  cut_cell_disc_.update(interfaces_.I(), interfaces_.In());

  // Mise a jour des indices et coefficients des points d'interpolation a une certaine distance des facettes de l'interface
  cut_cell_perform_interpolation_facettes();

  thermals_.echange_pure_vers_diph_cellules_initialement_pures();
  cut_field_velocity[0].copie_pure_vers_diph_sans_interpolation();
  cut_field_velocity[1].copie_pure_vers_diph_sans_interpolation();
  cut_field_velocity[2].copie_pure_vers_diph_sans_interpolation();

  interfaces_.calcul_surface_efficace_face_initial(type_surface_efficace_face_);
  interfaces_.calcul_surface_efficace_interface_initial(type_surface_efficace_interface_);

  const double fractionnal_timestep = compute_fractionnal_timestep_rk3(timestep, rk_step_);

  interfaces_.calcul_surface_efficace_face(type_surface_efficace_face_, fractionnal_timestep, cut_field_velocity);
  interfaces_.calcul_surface_efficace_interface(type_surface_efficace_interface_, fractionnal_timestep, cut_field_velocity);

  if (interfaces_.get_dt_impression_bilan_indicatrice() >= 0 && tstep_ % interfaces_.get_dt_impression_bilan_indicatrice() == interfaces_.get_dt_impression_bilan_indicatrice() - 1)
    {
      interfaces_.imprime_bilan_indicatrice();
    }
}

void Probleme_FTD_IJK_cut_cell::euler_time_step(ArrOfDouble& var_volume_par_bulle)
{
  static Stat_Counter_Id euler_rk3_counter_ = statistiques().new_counter(2, "Mise a jour de la vitesse");
  statistiques().begin_count(euler_rk3_counter_);
  if (thermals_.size())
    {
      if (boundary_conditions_.get_correction_conserv_qdm()==2)
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


    }

  /*
   * TODO: Change this block with OWN_PTR CLASS IJK_Thermal
   */
  for (auto& itr : thermique_)
    itr.euler_time_step(timestep_);

  for (auto& itr : energie_)
    itr.euler_time_step(velocity_);

  thermals_.euler_time_step(timestep_);

  if (!frozen_velocity_)
    {
      if (boundary_conditions_.get_correction_conserv_qdm()==2)
        {
          update_rho_v();
          rho_field_.echange_espace_virtuel(rho_field_.ghost());
          update_v_ghost_from_rho_v();
        }
      else
        {
          velocity_[0].echange_espace_virtuel(2);
          velocity_[1].echange_espace_virtuel(2);
          velocity_[2].echange_espace_virtuel(2);
        }

      // GAB, qdm
      if (test_etapes_et_bilan_)
        {
          calculer_rho_v(rho_field_,velocity_,rho_u_euler_av_prediction_champ_);
          for (int dir = 0; dir<3; dir++)
            rho_u_euler_av_prediction_[dir] = calculer_v_moyen(rho_u_euler_av_prediction_champ_[dir]);
        }
      // GAB, remarque : calculer dv calcule dv, MAIS NE L'APPLIQUE PAS au champ de vitesse !!!
      //                 l'increment de vitesse est ajoute au champ de vitesse avec euler_explicit_update
      calculer_dv(timestep_, current_time_, -1 /*rk_step = -1 pour sch euler... */);
      // GAB, qdm calculer_dv ne fait que l'etape de prediction)
      if (test_etapes_et_bilan_)
        {
          calculer_rho_v(rho_field_,d_velocity_,rho_du_euler_ap_prediction_champ_);
          for (int dir = 0; dir<3; dir++)
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

#else
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
          calculer_rho_v(rho_field_,d_velocity_,rho_du_euler_ap_projection_champ_);
          for (int dir=0; dir<3; dir++)
            rho_du_euler_ap_projection_[dir] = calculer_v_moyen(rho_du_euler_ap_projection_champ_[dir]);
        }

      if (test_etapes_et_bilan_)
        {
          calculer_rho_v(rho_field_,velocity_,rho_u_euler_ap_projection_champ_);
          for (int dir=0; dir<3; dir++)
            rho_u_euler_ap_projection_[dir] = calculer_v_moyen(rho_u_euler_ap_projection_champ_[dir]);
        }

      // Conditions en entree
      if (vitesse_entree_ > -1e20)
        force_entry_velocity(velocity_[0],
                             velocity_[1],
                             velocity_[2],
                             vitesse_entree_,
                             vitesse_entree_dir_,
                             vitesse_entree_compo_to_force_,
                             stencil_vitesse_entree_);

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
                  parser.setVar((int) 0, (*this).current_time_ - modified_time_ini_);
                  vitesse_upstream_ = parser.eval();
                }
            }
          else
            {
              int dir = 0;
              if (upstream_dir_ == -1)
                {
                  dir = get_direction_gravite();
                  if (dir == -1)
                    dir=0;
                }
              const DoubleTab& rising_vector = interfaces_.get_ijk_compo_connex().get_rising_vectors();
              const double velocity_magnitude = interfaces_.get_ijk_compo_connex().get_rising_velocities()[0];
              const Vecteur3& velocity_vector = interfaces_.get_ijk_compo_connex().get_rising_velocity_overall();
              velocity_bubble_new_ = velocity_vector[dir]; //  * rising_vector[dir];
              if (tstep_ == 0)
                {
                  if (velocity_bubble_old_ < -1e20)
                    velocity_bubble_old_ = 0.;
                  else
                    velocity_bubble_new_ = velocity_bubble_old_;
                  if (vitesse_upstream_reprise_ < -1e20)
                    vitesse_upstream_ = - velocity_bubble_scope_;
                  else
                    vitesse_upstream_ = vitesse_upstream_reprise_;
                }
              const double delta_velocity = velocity_bubble_scope_ + velocity_bubble_new_;
              const double ddelta_velocity = (velocity_bubble_new_ - velocity_bubble_old_) / timestep_;
              if (tstep_ % 100)
                velocity_bubble_integral_err_ = 0.;
              velocity_bubble_integral_err_ += delta_velocity * timestep_;
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

              calculer_vitesse_gauche(velocity_[0],velocity_[1],velocity_[2],vx,vy,vz);

              force_upstream_velocity_shear_perio(velocity_[0], velocity_[1], velocity_[2],
                                                  vitesse_upstream_, interfaces_, nb_diam_upstream_, boundary_conditions_, nb_diam_ortho_shear_perio_,
                                                  vx,vy,vz,epaisseur_maille_);
            }
          else
            {
              force_upstream_velocity(velocity_[0], velocity_[1], velocity_[2],
                                      vitesse_upstream_, interfaces_, nb_diam_upstream_,
                                      upstream_dir_, get_direction_gravite(), upstream_stencil_);
            }

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

              pressure_projection_with_inv_rho(inv_rho_field_, velocity_[0], velocity_[1],  velocity_[2], d_pressure_, timestep_,
                                               pressure_rhs_, check_divergence_, poisson_solver_);

            }
          else
            {
#ifdef PROJECTION_DE_LINCREMENT_DV
              // On l'a fait avant pour etre sur qu'elle soit bien dans la derivee stockee...
#else

              pressure_projection_with_rho(rho_field_, velocity_[0], velocity_[1],  velocity_[2], d_pressure_, timestep_,
                                           pressure_rhs_, check_divergence_, poisson_solver_);
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
            {

              pressure_projection_with_inv_rho(inv_rho_field_, velocity_[0], velocity_[1],  velocity_[2], pressure_, timestep_,
                                               pressure_rhs_, check_divergence_, poisson_solver_);

            }
          else
            {
#ifdef PROJECTION_DE_LINCREMENT_DV
#else

              pressure_projection_with_rho(rho_field_, velocity_[0], velocity_[1],  velocity_[2], pressure_, timestep_,
                                           pressure_rhs_, check_divergence_, poisson_solver_);

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
      //statistiques().end_count(projection_counter_);
    }

  Cerr << "Copy pressure on extended field for probes" << finl;
  copy_field_values(pressure_ghost_cells_, pressure_);

  if (Process::je_suis_maitre())
    {
      Cout << "Timings diff=" << statistiques().last_time(diffusion_counter_)
           << " conv=" << statistiques().last_time(convection_counter_);
      Cout << " src=" << statistiques().last_time(source_counter_)
           << finl;
    }
  statistiques().end_count(euler_rk3_counter_);
}

// Perform one sub-step of rk3 for FT algorithm, called 3 times per time step.
// rk_step = 0, 1 or 2
// total_timestep = not the fractionnal timestep !
void Probleme_FTD_IJK_cut_cell::rk3_sub_step(const int rk_step, const double total_timestep,
                                             const double fractionnal_timestep, const double time )
{
  assert(rk_step>=0 && rk_step<3);
  static Stat_Counter_Id euler_rk3_counter_ = statistiques().new_counter(2, "Mise a jour de la vitesse");
  statistiques().begin_count(euler_rk3_counter_);

  /*
   * TODO: Change this block with OWN_PTR CLASS IJK_Thermal
   */
  for (auto& itr : thermique_)
    {
      itr.rk3_sub_step(rk_step, total_timestep, time);
    }
  for (auto&& itr = energie_.begin(); itr != energie_.begin(); ++itr)
    {
      // curseur->rk3_sub_step(rk_step, total_timestep, time);
      // ++curseur;
      Cerr << "Le schema RK3 n est pas implemente avec des champs d energie" << finl;
      Process::exit();
    }

  thermals_.rk3_sub_step(rk_step, total_timestep, time);

  if (!frozen_velocity_)
    {
      if (boundary_conditions_.get_correction_conserv_qdm()==2)
        {
          update_rho_v();
          rho_field_.echange_espace_virtuel(rho_field_.ghost());
          update_v_ghost_from_rho_v();
        }
      else
        {
          velocity_[0].echange_espace_virtuel(2);
          velocity_[1].echange_espace_virtuel(2);
          velocity_[2].echange_espace_virtuel(2);
        }
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
            {
              pressure_projection_with_rho(rho_field_, d_velocity_[0], d_velocity_[1],  d_velocity_[2],
                                           d_pressure_,1. ,  pressure_rhs_, check_divergence_, poisson_solver_);
            }
          else
            {
              pressure_projection_with_rho(rho_field_, d_velocity_[0], d_velocity_[1],  d_velocity_[2],
                                           pressure_,1. ,  pressure_rhs_, check_divergence_, poisson_solver_);
            }
        }
#else
#endif

      // Mise a jour du champ de vitesse (etape de projection GAB, 28/06/21 : c'est la prediction, pas la projection non ?)
      for (int dir = 0; dir < 3; dir++)
        {
          const int kmax = d_velocity_[dir].nk();
          for (int k = 0; k < kmax; k++)
            {
              runge_kutta3_update(d_velocity_[dir], RK3_F_velocity_[dir], velocity_[dir], rk_step, k, total_timestep);
              // GAB, correction qdm a posteriori
            }
        }

#ifdef PROJECTION_DE_LINCREMENT_DV
      // Mise a jour du champ de pression
      if ((!disable_solveur_poisson_) && (include_pressure_gradient_in_ustar_))
        {
          const int kmax = pressure_.nk();
          for (int k = 0; k < kmax; k++)
            {
              runge_kutta3_update(d_pressure_ /* increment */,
                                  RK3_F_pressure_ /* intermediate storage */,
                                  pressure_ /* variable to update */, rk_step, k, total_timestep);
            }
        }
#else
#endif

      // Conditions en entree
      if (vitesse_entree_ > -1e20)
        force_entry_velocity(velocity_[0],
                             velocity_[1],
                             velocity_[2],
                             vitesse_entree_,
                             vitesse_entree_dir_,
                             vitesse_entree_compo_to_force_,
                             stencil_vitesse_entree_);


      // Forcage de la vitesse en amont de la bulle :
      if (vitesse_upstream_ > -1e20)
        {
          if (IJK_Shear_Periodic_helpler::defilement_ == 1)
            {

              double vx;
              double vy;
              double vz;

              calculer_vitesse_gauche(velocity_[0],velocity_[1],velocity_[2],vx,vy,vz);

              force_upstream_velocity_shear_perio(velocity_[0], velocity_[1], velocity_[2],
                                                  vitesse_upstream_, interfaces_, nb_diam_upstream_, boundary_conditions_, nb_diam_ortho_shear_perio_,
                                                  vx,vy,vz,epaisseur_maille_);
            }
          else
            {
              force_upstream_velocity(velocity_[0], velocity_[1], velocity_[2],
                                      vitesse_upstream_, interfaces_, nb_diam_upstream_,
                                      upstream_dir_, get_direction_gravite(), upstream_stencil_);
            }
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
        {
          // Cerr << "Methode basee sur inv rho pour le grad(P) en RK3" << finl;
          // Cerr << " Option a tester si besoin. " << finl;

          pressure_projection_with_inv_rho(inv_rho_field_, velocity_[0], velocity_[1],  velocity_[2], pressure_,
                                           fractionnal_timestep,
                                           pressure_rhs_, check_divergence_, poisson_solver_);

        }
      else
        {
#ifdef PROJECTION_DE_LINCREMENT_DV
          // On l'a fait avant pour etre sur qu'elle soit bien dans la derivee stockee...
#else


          pressure_projection_with_rho(rho_field_, velocity_[0], velocity_[1],  velocity_[2], pressure_,
                                       fractionnal_timestep,
                                       pressure_rhs_, check_divergence_, poisson_solver_);

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
      Cout << "Timings diff=" << statistiques().last_time(diffusion_counter_)
           << " conv=" << statistiques().last_time(convection_counter_);
      Cout << " src=" << statistiques().last_time(source_counter_)
           << finl;
    }
  statistiques().end_count(euler_rk3_counter_);

}

