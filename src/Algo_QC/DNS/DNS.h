/****************************************************************************
 * Copyright (c) 2015 - 2016, CEA
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 *modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *this list of conditions and the following disclaimer in the documentation
 *and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *may be used to endorse or promote products derived from this software without
 *specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *****************************************************************************/
/////////////////////////////////////////////////////////////////////////////
//
// File      : DNS_QC_double.h
// Directory : $NEW_ALGO_QC_ROOT/src/DNS_QC_double
//
/////////////////////////////////////////////////////////////////////////////
#ifndef DNS_QC_double_H
#define DNS_QC_double_H
#include <Boundary_Conditions.h>
#include <IJK_Field_vector.h>
#include <IJK_Field.h>
#include <Domaine_IJK.h>
#include <Operateur_IJK_elem_diff_base.h>
#include <Operateur_IJK_faces_diff.h>
#include <Operateur_IJK_faces_conv.h>
// #include <OpDiffTurbIJK.h>
// #include <OpDiffTurbIJKScalar.h>
#include <OpConvCentre4IJK.h>
#include <OpConvCentre2IJK.h>
#include <OpConvQuickIJKScalar.h>
#include <OpConvQuickSharpIJK.h>
#include <OpConvAmontIJK.h>
#include <Multigrille_Adrien.h>
#include <Interprete.h>
#include <IJK_Lata_writer.h>
#include <Linear_algebra_tools.h>
#include <Parser.h>
#include <Boundary_Conditions.h>
#include <Statistiques_dns_qc_ijk.h>
#include <Filter_kernel.h>
#include <OpConvCentre2IJKScalar.h>

#define WITH_FFTW

#if defined(WITH_FFTW)
#include <Fourier_trans.h>
#include <Redistribute_Field.h>
#endif


class DNS_QC_double : public Interprete
{
  Declare_instanciable(DNS_QC_double);
public:
  Entree& interpreter(Entree&) override;
  void run();
  void rk3_time_step();
protected:
  void initialise();
  void calculer_moyennes_flux();
  void rk3_sub_step(const int rk_step, const double total_timestep);

  template<class T>
  void calculer_convection_vitesse(IJK_Field_vector3_double& rho_v,
                                   IJK_Field_vector3_double& velocity,
                                   const ArrOfDouble_with_ghost& delta_z,
                                   const double facteur_delta_x,
                                   const double facteur_delta_y,
                                   const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                   T& kernel,
                                   FixedVector<IJK_Field_local_double, 18>& tmp_b,
                                   FixedVector<IJK_Field_local_double, 18>& tmp_a,
                                   IJK_Field_vector3_double& d_velocity_tmp,
                                   IJK_Field_vector3_double& d_velocity,
                                   IJK_Field_double& u_div_rho_u);
  template<class T>
  void calculer_turbulent_diffusion_vitesse(IJK_Field_vector3_double& velocity,
                                            const IJK_Field_double& turbulent_mu_xx,
                                            const IJK_Field_double& turbulent_mu_xy,
                                            const IJK_Field_double& turbulent_mu_xz,
                                            const IJK_Field_double& turbulent_mu_yy,
                                            const IJK_Field_double& turbulent_mu_yz,
                                            const IJK_Field_double& turbulent_mu_zz,
                                            const ArrOfDouble_with_ghost& delta_z,
                                            const double facteur_delta_x,
                                            const double facteur_delta_y,
                                            const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                            T& kernel,
                                            FixedVector<IJK_Field_local_double, 18>& tmp_b,
                                            FixedVector<IJK_Field_local_double, 18>& tmp_a,
                                            IJK_Field_vector3_double& d_velocity_tmp,
                                            IJK_Field_vector3_double& d_velocity);
  template<class T>
  void calculer_structural_diffusion_vitesse(IJK_Field_vector3_double& velocity,
                                             const FixedVector<IJK_Field_double, 6>& structural_uu_tensor,
                                             const ArrOfDouble_with_ghost& delta_z,
                                             const double facteur_delta_x,
                                             const double facteur_delta_y,
                                             const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                             T& kernel,
                                             FixedVector<IJK_Field_local_double, 18>& tmp_b,
                                             FixedVector<IJK_Field_local_double, 18>& tmp_a,
                                             IJK_Field_vector3_double& d_velocity_tmp,
                                             IJK_Field_vector3_double& d_velocity);
  void calculer_diffusion_scalar(const IJK_Field_double& rho,
                                 const IJK_Field_double& turbulent_kappa_x,
                                 const IJK_Field_double& turbulent_kappa_y,
                                 const IJK_Field_double& turbulent_kappa_z,
                                 IJK_Field_double& d_rho,
                                 const IJK_Field_local_double& boundary_flux_kmin,
                                 const IJK_Field_local_double& boundary_flux_kmax);


  void compute_sum_entering_heat_flux(double& flux_sum);
  void calcul_p_thermo_et_bilan(const IJK_Field_double& rho,
                                IJK_Field_double& temperature,
                                const int turbulent_diffusivity,
                                const IJK_Field_double& lambda_turbulent,
                                const int flag_lambda_anisotropic,
                                const int structural_uscalar,
                                const IJK_Field_double& structural_uscalar_z,
                                const double P_th_initial,
                                double& P_th_final,
                                const double fractionnal_timestep,
                                double& d_Pth_divise_par_gammamoins1) const;
  // initialise lambda aux paroi, crash si appeller apres le dt Numero 0 ;
  void calculer_lambda_paroi_air(double& lambda_de_t_paroi_kmin_ , double& lambda_de_t_paroi_kmax_) const;

  void ecrire_fichier_sauv(const char *fichier_sauvegarde,
                           const char *lata_name);
  void sauvegarder_qc(const char *fichier_sauvegarde);
  void reprendre_qc(const char *fichier_reprise);
  void calculer_velocity_elem();
  void posttraiter_champs_instantanes(const char * lata_name, double time);
  void save_raw_data(const char *lata_name, double current_time);
  void save_stats(double current_time);

  void fixer_reference_pression(double& reference_pression);
  void translation_pression(const double& reference_pression);
  void compute_rho_bulk(double& rho_bulk_);
  void compute_t_bulk(double& t_bulk_);
  void compute_rho_t_bulk(double& rho_bulk_, double& t_bulk_);

  Domaine_IJK domaine_;
  // taille des mailles de ce processeur en z, avec des mailles fantomes:
  ArrOfDouble_with_ghost delta_z_local_;
  // Cond.lim diffusion qdm, par defaut paroi fixe:
  Boundary_Conditions boundary_conditions_;
  // -------------------------------------------------
  // Champs inconnus (variables principales des equations differentielles
  //  qui doivent etre remplies au debut du pas de temps)
  // Velocity field:
  IJK_Field_vector3_double velocity_;
  // masse volumique
  IJK_Field_double rho_;
  // pression thermodynamique
  double P_thermodynamique_;
  // terme source qdm pour pousser le fluide dans le canal (en m/s/s)
  double terme_source_acceleration_;
  double terme_source_acceleration_constant_;
  int mode_terme_source_impose_;
  // when are we ?
  double current_time_;
  // source term in spanwise direction
  double spanwise_acceleration_source_term_ = 0;
  double spanwise_mass_flux_aim_ = 0;
  double spanwise_mass_flux_now_ = 0;

  // -------------------------------------------------
  // Statistiques temporelles
  Statistiques_dns_qc_ijk statistiques_;
  double t_debut_statistiques_;

  // Statistiques spectrale de la turbulence
  Fourier_trans partie_fourier_;
  int dt_post_spectral_;
  int reprise_fourier_;
  Domaine_IJK post_splitting_;

  // Sauvegarde des lata par plan
  Domaine_IJK sauvegarde_splitting_;
  Nom sauvegarde_splitting_name_;
  FixedVector<Redistribute_Field, 3> redistribute_to_sauvegarde_splitting_faces_;
  Redistribute_Field redistribute_to_sauvegarde_splitting_elem_;
  IJK_Field_vector3_double velocity_sauvegarde_;
  IJK_Field_double temperature_sauvegarde_;
  IJK_Field_double molecular_lambda_sauvegarde_;
  IJK_Field_double molecular_mu_sauvegarde_;
  IJK_Field_double rho_sauvegarde_;
  IJK_Field_double pressure_sauvegarde_;
  IJK_Field_double velocity_elem_X_sauvegarde_;
  IJK_Field_double velocity_elem_Y_sauvegarde_;
  IJK_Field_double velocity_elem_Z_sauvegarde_;


  // -------------------------------------------------
  // Champs supplementaires (variables temporaires de calcul deduites des inconnues)
  IJK_Field_double velocity_elem_X_;
  // Vitesse aux elements
  IJK_Field_double velocity_elem_Y_;
  // Vitesse aux elements
  IJK_Field_double velocity_elem_Z_;
  // Vitesse aux elements
  // momentum:
  IJK_Field_vector3_double rho_v_;
  // Temporary storage for the derivative
  IJK_Field_vector3_double d_velocity_;
  // Temporary storage for the RungeKutta algorithm
  IJK_Field_vector3_double RK3_F_velocity_;
  // Pressure field
  IJK_Field_double pressure_;
  // viscosite dynamique: div(mu*grad(v)) => d/dt (rho*v)
  IJK_Field_double molecular_mu_;
  // right hand side for pressure solver
  IJK_Field_double pressure_rhs_;
  // derivee de la masse volumique par rapport au temps
  IJK_Field_double d_rho_;
  // temperature
  IJK_Field_double temperature_;
  // temporary storage for the rk3 algorithm
  IJK_Field_double RK3_F_rho_;
  // conductivite thermique
  IJK_Field_double molecular_lambda_;
  // resultat de l'operateur diffusion temperature:
  IJK_Field_double div_lambda_grad_T_volume_;
  // tableau temporaire utilise par l'operateur de convection pour y stocker u*div(rho*u)
  IJK_Field_double u_div_rho_u_;

  // pour rajouter le terme en divergence a la diffusion
  IJK_Field_double divergence_;

  // Flux thermique au bord (integrale par face de bord)
  IJK_Field_local_double boundary_flux_kmin_;
  IJK_Field_local_double boundary_flux_kmax_;

  // Operators and pressure solver
  // DD,2017-04-27: diffusion modifie en vue de l'ajout de modeles
  OpDiffIJK_double velocity_diffusion_op_simple_;
  OpDiffStdWithLaminarTransposeIJK_double velocity_diffusion_op_simple_with_transpose_;
  OpDiffStdWithLaminarTransposeAndDivergenceIJK_double velocity_diffusion_op_full_;
  Nom type_velocity_diffusion_;

  OpDiffTensorialZeroatwallIJK_double velocity_turbulent_diffusion_op_simple_;
  OpDiffStdWithLaminarTransposeTensorialZeroatwallIJK_double velocity_turbulent_diffusion_op_simple_with_transpose_;
  OpDiffStdWithLaminarTransposeAndDivergenceTensorialZeroatwallIJK_double velocity_turbulent_diffusion_op_full_;
  OpDiffTensorialAnisotropicZeroatwallIJK_double velocity_turbulent_diffusion_op_simple_anisotropic_;
  OpDiffStdWithLaminarTransposeTensorialAnisotropicZeroatwallIJK_double velocity_turbulent_diffusion_op_simple_with_transpose_anisotropic_;
  OpDiffStdWithLaminarTransposeAndDivergenceTensorialAnisotropicZeroatwallIJK_double velocity_turbulent_diffusion_op_full_anisotropic_;
  Nom type_velocity_turbulent_diffusion_;
  Nom turbulent_viscosity_model_;
  Nom turbulent_viscosity_dynamic_type_;
  double turbulent_viscosity_model_constant_;
  ArrOfDouble turbulent_viscosity_tensor_coefficients_;
  int flag_nu_tensorial_;
  int flag_nu_anisotropic_;
  IJK_Field_double turbulent_mu_;
  FixedVector<IJK_Field_double, 6> turbulent_mu_tensor_;
  int turbulent_viscosity_;
  int variation_cste_modele_fonctionnel_;
  // M.D 10/2021
  double smoothing_center_fr_;
  double smoothing_factor_fr_;
  double Re_tau_fr_;
  double Re_tau_ch_;
  double pond_fr_;
  double pond_ch_;
  double center_constant_;

  OpDiffVectorialIJKScalar_double operateur_diffusion_turbulent_scalar_;
  OpDiffVectorialAnisotropicIJKScalar_double operateur_diffusion_turbulent_scalar_anisotropic_;
  Nom type_scalar_turbulent_diffusion_;
  Nom turbulent_diffusivity_model_;
  Nom turbulent_diffusivity_dynamic_type_;
  double turbulent_diffusivity_model_constant_;
  ArrOfDouble turbulent_diffusivity_vector_coefficients_;
  int flag_kappa_vectorial_;
  int flag_kappa_anisotropic_;
  IJK_Field_double turbulent_kappa_;
  IJK_Field_vector3_double turbulent_kappa_vector_;
  int turbulent_diffusivity_;

  Nom filter_kernel_name_;
  Filter_kernel_base* kernel_;
  int flag_filtrage_convection_qdm_;
  int flag_filtrage_turbulent_diffusion_qdm_;
  int flag_filtrage_structural_diffusion_qdm_;
  int flag_convection_qdm_sans_rho_;
  int flag_convection_qdm_sans_divergence_;
  IJK_Field_vector3_double velocity_filtre_;
  IJK_Field_double rho_filtre_;
  IJK_Field_double rho_velocity_i_filtre_ ;
  IJK_Field_double rho_velocity_j_filtre_ ;
  IJK_Field_double rho_velocity_k_filtre_ ;
  IJK_Field_double temperature_filtre_;
  IJK_Field_double turbulent_mu_filtre_;
  FixedVector<IJK_Field_double, 6> turbulent_mu_filtre_tensor_;
  IJK_Field_double turbulent_kappa_filtre_;
  IJK_Field_vector3_double turbulent_kappa_filtre_vector_;
  FixedVector<IJK_Field_double, 6> structural_uu_filtre_tensor_; // Vector with 6 components, 0:xx 1:xy 2:xz 3:yy 4:yz 5:zz
  IJK_Field_vector3_double structural_uscalar_filtre_vector_;
  int flag_u_filtre_;
  int flag_rho_filtre_;
  int flag_temperature_filtre_;
  int flag_turbulent_mu_filtre_;
  int flag_turbulent_kappa_filtre_;
  int flag_structural_uu_filtre_;
  int flag_structural_uscalar_filtre_;

  OpDiffStructuralOnlyZeroatwallIJK_double velocity_turbulent_diffusion_op_structural_;
  Nom structural_uu_model_;
  Nom structural_uu_dynamic_type_;
  double structural_uu_model_constant_;
  ArrOfDouble structural_uu_tensor_coefficients_;
  FixedVector<IJK_Field_double, 6> structural_uu_tensor_; // Vector with 6 components, 0:xx 1:xy 2:xz 3:yy 4:yz 5:zz
  int flag_structural_uu_tmp_;
  FixedVector<IJK_Field_double, 6> structural_uu_tmp_tensor_; // Vector with 6 components, 0:xx 1:xy 2:xz 3:yy 4:yz 5:zz
  int structural_uu_;

  Nom structural_uscalar_model_;
  Nom structural_uscalar_dynamic_type_;
  double structural_uscalar_model_constant_;
  ArrOfDouble structural_uscalar_vector_coefficients_;
  IJK_Field_vector3_double structural_uscalar_vector_;
  int flag_structural_uscalar_tmp_;
  IJK_Field_vector3_double structural_uscalar_tmp_vector_;
  int structural_uscalar_;

  Nom large_eddy_simulation_formulation_;
  double facteur_delta_x_;
  double facteur_delta_y_;
  ArrOfDouble_with_ghost delta_z_local_pour_delta_;
  double facteur_delta_filtre_x_;
  double facteur_delta_filtre_y_;
  ArrOfDouble_with_ghost delta_z_local_pour_delta_filtre_;
  int formulation_favre_;
  int formulation_velocity_;
  ArrOfDouble_with_ghost constante_modele_;
  FixedVector<IJK_Field_local_double, 18> tmp_b_; // Temporary array used to compute the filter
  FixedVector<IJK_Field_local_double, 18> tmp_a_; // Temporary array used to compute the filter
  FixedVector<FixedVector<ArrOfDouble, 7>, 8> ml_; // Vector with 8 components, 0:l 1:m 2:h 3:mm 4:hh 5:ml 6:hl 7:mh. Each is a vector with 7 components, 0:xx 1:xy 2:xz 3:yy 4:yz 5:zz 6:sum
  int flag_d_velocity_tmp_;
  IJK_Field_vector3_double d_velocity_tmp_;

  OpConvCentre4IJK_double velocity_convection_op_;
  OpConvCentre2IJK_double velocity_convection_op_centre_2_;

  OpConvQuickSharpIJK_double velocity_convection_op_quicksharp_;
  OpConvAmontIJK_double velocity_convection_op_amont_;

  OpConvQuickIJKScalar_double rho_convection_op_;
  OpConvCentre2IJKScalar_double rho_convection_op_centre2_;
  OpConvCentre4IJK_double rho_convection_op_centre4_;
  OpConvAmontIJK_double rho_convection_op_amont_;
  OpDiffIJKScalar_double operateur_diffusion_temperature_;
  OpDiffStructuralOnlyIJKScalar_double operateur_diffusion_temperature_structural_;
  Multigrille_Adrien poisson_solver_;

  // Simulation parameters
  // ----------------------------------------------
  // Valeurs lues dans le jeu de donnees:
  int nb_timesteps_;

  double timestep_max_;
  double timestep_facsec_;
  double timestep_;
  double old_timestep_;
  int dt_post_;
  Motcles liste_post_instantanes_; // liste des champs instantanes a postraiter
  int dt_sauvegarde_;
  int dt_raw_data_;
  int dt_stats_;

  double dt_save_oscillating_cycle_raw_data_; // saves raw_data every physical dt
  double dt_save_cycle_;                      // incremented timestep for saving
  int postraiter_sous_pas_de_temps_; // drapeau 0 ou 1
  Nom nom_sauvegarde_;
  Nom nom_reprise_;
  // Pour numeroter les fichiers .lata il faut compter combien on en a ecrit:
  int compteur_post_instantanes_;
  int calcul_2d_; // drapeau 0 ou 1
  Nom check_stop_file_; // Nom du fichier stop
  int check_divergence_;
  int projection_initiale_demandee_;
  double puit_;

  // Le jeu de donnees doit fournir soit des fichiers de reprise: ..
  Nom fichier_reprise_rho_;
  int timestep_reprise_rho_;
  Nom fichier_reprise_vitesse_;
  int timestep_reprise_vitesse_;
  // ... soit des expressions f(x,y,z)
  Nom expression_temperature_initiale_;
  Noms expression_vitesse_initiale_; // on attend trois expressions

  Nom expression_derivee_acceleration_;
  Parser parser_derivee_acceleration_;
  // F.A nouvelle source
  double dump_factor_ ;
  double debit_actuel_ ;
  double debit_cible_ ;
  // Variables needed for forcing T bulk
  double aim_t_bulk_;
  double actual_t_bulk_;
  bool flag_t_bulk_forced_;
  double dump_factor_2_;
  double rho_bulk_;
  /* double sum_entering_flux_; */

  int convection_rho_centre2_;
  int convection_rho_amont_;
  int convection_rho_centre4_;

  int convection_velocity_amont_;
  int convection_velocity_centre2_;
  int convection_velocity_quicksharp_;



  // Nom expression_derivee_acceleration_;
  // Parser parser_derivee_acceleration_;
  // pour me simplifier la vie
  double Lx_tot_ ;
  double Ly_tot_ ;
  double Lz_tot_ ;

  double T_paroi_impose_kmin_;
  double T_paroi_impose_kmax_;
  double Cp_gaz_;
  double gamma_;

  // F.A 17/03/14 ajout diffs et convection negligeables
  int diff_qdm_negligeable_;
  int diff_temp_negligeable_;

  int conv_qdm_negligeable_;
  int conv_rho_negligeable_;

  // DD,2016-10-14: ajout disable_solveur_poisson_
  int disable_solveur_poisson_;

  // F.A 17/03/14 ajout possibiltie de dt_start
  double dt_start_;

  // ----------------------------------------------
  // Valeurs recalculees a partir des donnees:
  double rho_paroi_impose_kmin_;
  double rho_paroi_impose_kmax_;
  double lambda_de_t_paroi_kmin_;
  double lambda_de_t_paroi_kmax_;
  double constante_specifique_gaz_;
  double volume_total_domaine_;

  // DD,2016-28-01: statistiques a partir de fichiers lata uniquement
  Noms statlata_namelist_;
  int sauvegarde_post_instantanes_;
  int lecture_post_instantanes_;
  int lecture_post_instantanes_filtrer_u_;
  int lecture_post_instantanes_filtrer_rho_;
  int lecture_post_instantanes_filtrer_p_;
  int lecture_post_instantanes_filtrer_tous_;

  // DD,2016-06-15: changement des pas de temps de stabilite,
  // pour faciliter les tests, j'ajoute un parametre pour avoir l'ancienne version
  int old_dtstab_;

  // Adding expression for oscillating boundary condition
  // Nom expression_oscillating_boundary;
  double amplitude_oscillating_boundary_;
  double frequency_oscillating_boundary_;
  bool flag_oscillating_boundary;

  // Saving each delta_t time
  bool flag_save_each_delta_t_ = false;
};
#endif
