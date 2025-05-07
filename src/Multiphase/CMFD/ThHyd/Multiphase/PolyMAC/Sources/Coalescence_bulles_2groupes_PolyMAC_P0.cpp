/****************************************************************************
* Copyright (c) 2021, CEA
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
// File:        Coalescence_Yao_Morel.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Coalescence_bulles_2groupes_PolyMAC_P0.h>
#include <Pb_Multiphase.h>
#include <Milieu_composite.h>
#include <Op_Diff_Turbulent_PolyMAC_P0_Face.h>
#include <Viscosite_turbulente_base.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <Coalescence_bulles_2groupes_base.h>
#include <math.h>
#include <Champ_Elem_PolyMAC_P0.h>
#include <Champ_Face_base.h>
#include <Domaine_PolyMAC_P0.h>

Implemente_instanciable(Coalescence_bulles_2groupes_PolyMAC_P0, "Coalescence_bulles_2groupes_elem_PolyMAC_P0", Source_base);

Sortie& Coalescence_bulles_2groupes_PolyMAC_P0::printOn(Sortie& os) const
{
  return os;
}

Entree& Coalescence_bulles_2groupes_PolyMAC_P0::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("beta_k", &beta_k_);
  param.ajouter("Dh", &dh_);
  param.lire_avec_accolades_depuis(is);


  const Pb_Multiphase *pbm = sub_type(Pb_Multiphase, equation().probleme()) ? &ref_cast(Pb_Multiphase, equation().probleme()) : nullptr;

  if (!pbm || pbm->nb_phases() == 1) Process::exit(que_suis_je() + " : not needed for single-phase flow!");
  for (int n = 0; n < pbm->nb_phases(); n++)  //recherche de n_l, n_g : phase {liquide,gaz}_continu en priorite
    {
      if (pbm->nom_phase(n).debute_par("liquide") && (n_l < 0 || pbm->nom_phase(n).finit_par("continu")))  n_l = n;
      if (( pbm->nom_phase(n).finit_par("group1")))  n_g1 = n;
      if (( pbm->nom_phase(n).finit_par("group2")))  n_g2 = n;
    }
  if (n_l < 0) Process::exit(que_suis_je() + " : liquid phase not found!");
  if (n_g1 < 0) Process::exit(que_suis_je() + " : group 1 not found!");
  if (n_g2 < 0) Process::exit(que_suis_je() + " : group 2 not found!");

  if (pbm->has_correlation("Coalescence_bulles_2groupes")) correlation_ = pbm->get_correlation("Coalescence_bulles_2groupes"); //correlation fournie par le bloc correlation
  else Correlation_base::typer_lire_correlation(correlation_, *pbm, "Coalescence_bulles_2groupes", is); //sinon -> on la lit

  return is;
}

void Coalescence_bulles_2groupes_PolyMAC_P0::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Domaine_PolyMAC_P0& domaine = ref_cast(Domaine_PolyMAC_P0, equation().domaine_dis());
  const int ne = domaine.nb_elem(), ne_tot = domaine.nb_elem_tot(), N = equation().inconnue().valeurs().line_size();

  for (auto &&n_m : matrices)
    if (n_m.first == "alpha" || n_m.first == "k" || n_m.first == "tau" || n_m.first == "omega" || n_m.first == "interfacial_area")
      {
        Matrice_Morse& mat = *n_m.second, mat2;
        const DoubleTab& dep = equation().probleme().get_champ(n_m.first.c_str()).valeurs();
        int nc = dep.dimension_tot(0),
            M  = dep.line_size();
        IntTab sten(0, 2);
        if (n_m.first == "alpha")
          for (int e = 0; e < ne; e++)
            for (int n = 0; n < N; n++) sten.append_line(N * e + n, N * e + n);
        if (n_m.first == "interfacial_area" ) // N <= M
          for (int e = 0; e < ne; e++)
            for (int n = 0; n < N; n++) sten.append_line(N * e + n, M * e + n);
        if (n_m.first == "k" || n_m.first == "tau" || n_m.first == "omega") // N <= M
          for (int e = 0; e < ne; e++)
            for (int n = 0; n < N; n++)
              for (int k = 0; k < M; k++) sten.append_line(N * e + n, M * e + k);
        //tableau_trier_retirer_doublons(sten);
        Matrix_tools::allocate_morse_matrix(N * ne_tot, M * nc, sten, mat2);
        mat.nb_colonnes() ? mat += mat2 : mat = mat2;
      }
}

void Coalescence_bulles_2groupes_PolyMAC_P0::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Pb_Multiphase& pbm = ref_cast(Pb_Multiphase, equation().probleme());
  const Domaine_PolyMAC_P0& domaine = ref_cast(Domaine_PolyMAC_P0, equation().domaine_dis());
  const DoubleVect& pe = equation().milieu().porosite_elem(), &ve = domaine.volumes();


  const DoubleTab& inco = equation().inconnue().valeurs(),
                   &inco_p = equation().inconnue().passe(),
                    &d_b = equation().probleme().get_champ("diametre_bulles").passe(),
                     &alpha = pbm.equation_masse().inconnue().valeurs(),
                      &alpha_p = pbm.equation_masse().inconnue().passe(),
                       &press_p = ref_cast(QDM_Multiphase,pbm.equation_qdm()).pression().passe(),
                        &temp_p  = pbm.equation_energie().inconnue().passe(),
                         &rho_p   = equation().milieu().masse_volumique().passe(),
                          &nu_p = equation().probleme().get_champ("viscosite_cinematique").passe(),
                           *tab_k_p = equation().probleme().has_champ("k") ? &equation().probleme().get_champ("k").passe() : nullptr,
                            *tab_k = equation().probleme().has_champ("k") ? &equation().probleme().get_champ("k").valeurs() : nullptr,
                             *tau = equation().probleme().has_champ("tau") ? &equation().probleme().get_champ("tau").valeurs() : nullptr,
                              *omega = equation().probleme().has_champ("omega") ? &equation().probleme().get_champ("omega").valeurs() : nullptr ;

  const Milieu_composite& milc = ref_cast(Milieu_composite, equation().milieu());
  int N = pbm.nb_phases(), Nk = (tab_k_p) ? (*tab_k_p).line_size() : -1, Np = equation().probleme().get_champ("pression").valeurs().line_size();

  // Models use epsilon but with omega and tau we induce new variations of tau/omega and k
  std::string Type_diss = "other"; // omega, tau or other dissipation
  if (tau) Type_diss = "tau";
  else if (omega) Type_diss = "omega";

  DoubleTrav epsilon(alpha);
  const Op_Diff_Turbulent_PolyMAC_P0_Face& op_diff 	= ref_cast(Op_Diff_Turbulent_PolyMAC_P0_Face, equation().probleme().equation(0).operateur(0).l_op_base());
  const Viscosite_turbulente_base&   	visc_turb 		= ref_cast(Viscosite_turbulente_base, op_diff.correlation());
  visc_turb.eps(epsilon);  // Epsilon is in the past
  double limiter = visc_turb.limiteur();
  double dh = dh_;


  Matrice_Morse  *Ma = matrices.count("alpha") ? matrices.at("alpha") : nullptr,
                  *Mk = matrices.count("k") ? matrices.at("k") : nullptr,
                   *Mtau = matrices.count("tau") ? matrices.at("tau") : nullptr,
                    *Momega = matrices.count("omega") ? matrices.at("omega") : nullptr,
                     *Mai = matrices.count("interfacial_area") ? matrices.at("interfacial_area") : nullptr;

  int cR = (rho_p.dimension_tot(0) == 1), cM = (nu_p.dimension_tot(0) == 1), n, k, e, d, D = dimension;
  DoubleTrav a_l(N), a_l_p(N), p_l(N), T_l(N), rho_l(N), nu_l(N), sigma_l(N,N), dv(N, N), d_bulles(N), eps_l(Nk), k_l(Nk), coeff_RC(N, N), coeff_WE(N, N); //arguments pour coeff
  const Coalescence_bulles_2groupes_base& correlation_coal = ref_cast(Coalescence_bulles_2groupes_base, correlation_.valeur());

  // fill velocity at elem tab
  DoubleTab pvit_elem(0, N * D);
  domaine.domaine().creer_tableau_elements(pvit_elem);
  const Champ_Face_base& ch_vit = ref_cast(Champ_Face_base,ref_cast(Pb_Multiphase, equation().probleme()).equation_qdm().inconnue());
  ch_vit.get_elem_vector_field(pvit_elem);

  const double fac_sec = 1.e4 ; // numerical security
  const double alpha_sec = 0.509 ; // cf C.Bazin thesis
  const double alpha_min = 1.e-3 ; // to avoid numerical problems

  /* elements */
  for (e = 0; e < domaine.nb_elem(); e++)
    {
      // Get field values for correlations-------------------------------------------------------------------------------------------------------------------------------------------------------
      for (n = 0; n < N; n++)
        {
          a_l(n)   = alpha(e, n);
          a_l_p(n)   = alpha_p(e, n);
          p_l(n)   = press_p(e, n * (Np > 1));
          T_l(n)   =  temp_p(e, n);
          rho_l(n) =   rho_p(!cR * e, n);
          nu_l(n)  =    nu_p(!cM * e, n);
          for (k = 0; k < N; k++)
            if(milc.has_interface(n, k))
              {
                Interface_base& sat = milc.get_interface(n, k);
                sigma_l(n,k) = sat.sigma(temp_p(e,n), press_p(e,n * (Np > 1)));
              }
            else if (milc.has_saturation(n, k))
              {
                Saturation_base& z_sat = milc.get_saturation(n, k);

                DoubleTab& sig = z_sat.get_sigma_tab();
                sigma_l(n,k) = sig(e);

              }

          d_bulles(n) = d_b(e,n);

        }

      for (n = 0; n < Nk; n++)
        {
          eps_l(n) =epsilon(e, n) ;
          k_l(n)   = (tab_k_p)   ? (*tab_k_p)(e,n) : 0;
        }

      for (dv =0, d = 0; d < D; d++)
        for (n = 0; n < N; n++)
          for (k = 0 ; k<N ; k++) dv(n, k) += (pvit_elem(e, N * d + n) - ((n!=k) ? pvit_elem(e, N * d + k) : 0) ) * (pvit_elem(e, N * d + n) - ((n!=k) ? pvit_elem(e, N * d + k) : 0) ); // nv(n,n) = ||v(n)||, nv(n, k!=n) = ||v(n)-v(k)||
      for (n = 0; n < N; n++)
        for (k = 0 ; k<N ; k++) dv(n, k) = sqrt(dv(n, k)) ;

      // Get correlations-------------------------------------------------------------------------------------------------------------------------------------------------------
      correlation_coal.coefficient_RC(a_l, a_l_p, p_l, T_l, rho_l, nu_l, sigma_l, dh, dv, d_bulles, eps_l, k_l, n_l, n_g1, n_g2,  coeff_RC); // Semi-Explicit coeff for Random Collisions : alpha is implicit
      correlation_coal.coefficient_WE(a_l_p, p_l, T_l, rho_l, nu_l, sigma_l, dh, dv, d_bulles, eps_l, k_l, n_l, n_g1, n_g2, coeff_WE); // Semi-Explicit coeff for Wake entrainment : alpha is implicit

      // Get epsilon according to model
      double eps_valeurs = epsilon(e, n_l) ;
      if (Type_diss == "tau")        eps_valeurs = beta_k_ * ((*tab_k)(e, n_l)>1.e-8 ? (*tab_k)(e, n_l)*(*tab_k)(e, n_l)/ std::max((*tab_k)(e, n_l) * (*tau)(e, n_l), limiter * nu_p(e, n_l)) : 0 );
      else if (Type_diss == "omega") eps_valeurs = beta_k_ * ((*tab_k)(e, n_l)*(*omega)(e, n_l)) ;
      else eps_valeurs = epsilon(e, n_l);


      // Recuring functions for source terms-------------------------------------------------------------------------------------------------------------------------------------------------------
      // fac = prefactor, da = void derivative

      const double fac_RC1 = (alpha(e, n_g1)>alpha_min) ? pe(e) * ve(e) * coeff_RC(n_g1, n_l) : 0.;
      const double dfac_RC1da = (alpha(e, n_g1)>alpha_min) ? pe(e) * ve(e) * coeff_RC(n_l, n_g1) : 0. ;
      const double fac_WE1 = (alpha(e, n_g1)>alpha_min) ? pe(e) * ve(e) * coeff_WE(n_g1, n_l) : 0. ;
      const double dfac_WE1da = (alpha(e, n_g1)>alpha_min) ?  pe(e) * ve(e) * coeff_WE(n_l, n_g1) : 0.;
      const double fac_RC2 = (alpha(e, n_g2)>alpha_min) ?  pe(e) * ve(e) * coeff_RC(n_g2, n_l) : 0.;
      const double dfac_RC2da = (alpha(e, n_g2)>alpha_min) ? pe(e) * ve(e) * coeff_RC(n_l, n_g2) : 0. ;
      const double fac_WE2 = (alpha(e, n_g2)>alpha_min) ? pe(e) * ve(e) * coeff_WE(n_g2, n_l) : 0. ;
      const double fac_RC112 = (alpha(e, n_g1)>alpha_min) ? pe(e) * ve(e) * coeff_RC(n_g1, n_g1) : 0.;
      const double fac_RC122_1 =(alpha(e, n_g2)>alpha_min) ? pe(e) * ve(e) * coeff_RC(n_g1, n_g2) : 0.;
      const double dfac_RC122_1da = (alpha(e, n_g2)>alpha_min) ? pe(e) * ve(e) * coeff_RC(n_l, n_l) : 0. ;
      const double fac_RC122_2 = (alpha(e, n_g2)>alpha_min) ? pe(e) * ve(e) * coeff_RC(n_g2, n_g1) : 0. ;
      const double fac_WE112 = (alpha(e, n_g1)>alpha_min) ? pe(e) * ve(e) * coeff_WE(n_g1, n_g1) : 0.;
      const double fac_WE122_1 = (alpha(e, n_g2)>alpha_min) ? pe(e) * ve(e) * coeff_WE(n_g1, n_g2) : 0.;
      const double fac_WE122_2 = (alpha(e, n_g2)>alpha_min) ? pe(e) * ve(e) * coeff_WE(n_g2, n_g1) : 0.;


      const double ai1_p = std::max(inco_p(e, n_g1), 0.) ; // security inco negative
      const double ai1 = std::max(inco(e, n_g1), 0.) ; // security inco negative

      const double ai2_p = std::max(inco_p(e, n_g2), 0.) ; // security inco negative
      const double ai2 = std::max(inco(e, n_g2), 0.) ; // security inco negative
      const double eps_1_over3 = std::cbrt(eps_valeurs) ;

      const double alphag1_1_over3 = std::cbrt(std::min(alpha(e, n_g1), alpha_sec)) ;
      const double alphag1_p_1_over3 = std::cbrt(std::min(alpha_p(e, n_g1), alpha_sec)) ;

      const double alphag2_1_over3 = std::cbrt(alpha(e, n_g2)) ;
      const double alphag2p_1_over3 = std::cbrt(alpha_p(e, n_g2)) ;

      const double ai2_1_over3 = std::cbrt(ai2) ;
      const double ai2p_1_over3 = std::cbrt(ai2_p) ;

      const double ai1_5_over3 = std::cbrt(ai1) * std::cbrt(ai1) * std::cbrt(ai1) * std::cbrt(ai1) * std::cbrt(ai1) ;
      const double ai1p_5_over3 = std::cbrt(ai1_p) * std::cbrt(ai1_p) * std::cbrt(ai1_p) * std::cbrt(ai1_p) * std::cbrt(ai1_p) ;
      const double ai2_5_over3 = std::cbrt(ai2) * std::cbrt(ai2) * std::cbrt(ai2) * std::cbrt(ai2) * std::cbrt(ai2) ;

      // Fill the matrix--------------------------------------------------------------------------------------------------------------------------------------------------

      // RC (1)--------------------------------------------------------------------------------------------------------------------------

      secmem(e , n_g1) += (fac_RC1 > 0. ) ? fac_RC1 * std::min(alpha(e, n_g1), alpha_sec) *  ai1_5_over3 * eps_1_over3 : 0. ; // (alpha1, ai1, epsilon) implicit dependance

      // RC (2)---------------------------------------------------------------------------------------------------------------------------

      secmem(e , n_g2) += (fac_RC2 > 0. ) ? fac_RC2 * (alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3) /std::min(ai2_1_over3,fac_sec) * eps_1_over3 : 0. ; // (alpha2, ai2, epsilon) implicit dependance

      // RC (112)----------------------------------------------------------------------------------------------------------------------------

      secmem(e , n_g2) += (fac_RC112 > 0. ) ? fac_RC112 * eps_1_over3 * std::min(alpha_p(e, n_g1), alpha_sec) * ai1p_5_over3: 0. ; // (epsilon) implicit dependance

      // RC (122) for g1----------------------------------------------------------------------------------------------------------------------------

      secmem(e , n_g1) += (fac_RC122_1 > 0. ) ? fac_RC122_1 * eps_1_over3 * (alphag1_1_over3 * alphag1_1_over3) * (alphag2p_1_over3 * alphag2p_1_over3 * alphag2p_1_over3 * alphag2p_1_over3) * ai1 * (ai2p_1_over3 * ai2p_1_over3) : 0.; // (alpha1, ai1, epsilon) implicit dependance

      // RC (122) for g2----------------------------------------------------------------------------------------------------------------------------

      secmem(e , n_g2) += (fac_RC122_2 > 0. ) ? fac_RC122_2 * eps_1_over3 * (alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3) * (alphag2_1_over3)  * (ai2_5_over3) : 0.  ; // (alpha2, ai2, epsilon) implicit dependance

      // WE (1)----------------------------------------------------------------------------------------------------------------------------

      secmem(e , n_g1) += (fac_WE1 > 0. ) ? fac_WE1  * ai1 * ai1 : 0. ; // (ai1) implicit dependance

      // WE (2)----------------------------------------------------------------------------------------------------------------------------

      secmem(e , n_g2) += (fac_WE2 > 0. ) ? fac_WE2  * ai2 * ai2 / std::min(alpha(e, n_g2) ,fac_sec): 0.; // (alpha2, ai2) implicit dependance

      // WE (112)----------------------------------------------------------------------------------------------------------------------------

      secmem(e , n_g2) += (fac_WE112 > 0. ) ? fac_WE112  * ai1_p * ai1_p: 0. ;// no implicit dependance

      // WE (122) for g1----------------------------------------------------------------------------------------------------------------------------

      secmem(e , n_g1) += (fac_WE122_1 > 0. ) ? fac_WE122_1  * ai1 * ai2_p: 0. ; // (ai1) implicit dependance

      // WE (122) for g2----------------------------------------------------------------------------------------------------------------------------

      secmem(e , n_g2) += (fac_WE122_2 > 0. ) ? fac_WE122_2  * ai2 * std::min(alpha_p(e, n_g1), alpha_sec)  / std::min(alpha(e, n_g2),fac_sec): 0.; // (alpha2, ai2) implicit dependance


      if (Ma)//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        {
          // RC (1)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

          (*Ma)(N * e + n_g1 , N * e + n_g1) -= (fac_RC1 > 0. ) ? fac_RC1 *  ai1_5_over3 * eps_1_over3 + dfac_RC1da * std::min(alpha(e, n_g1), alpha_sec) *  ai1_5_over3 * eps_1_over3: 0.  ;

          // RC (2)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

          (*Ma)(N * e + n_g2 , N * e + n_g2) -= (fac_RC2 > 0. ) ? fac_RC2 * 7./3. * (alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3) / std::min(ai2_1_over3,fac_sec) * eps_1_over3 + dfac_RC2da * (alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3) /  std::min(ai2_1_over3,fac_sec) * eps_1_over3 : 0. ;

          // RC (122) for g1---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

          (*Ma)(N * e + n_g1 , N * e + n_g1) -= (fac_RC122_1 > 0. ) ? fac_RC122_1 * eps_1_over3 * 2./3. / std::min(alphag1_1_over3,fac_sec) * (ai2p_1_over3 * ai2p_1_over3 * ai2p_1_over3 * ai2p_1_over3) * ai1 * (ai2p_1_over3 * ai2p_1_over3) + dfac_RC122_1da * eps_1_over3 * (alphag1_1_over3 * alphag1_1_over3) * (alphag2p_1_over3 * alphag2p_1_over3 * alphag2p_1_over3 * alphag2p_1_over3) * ai1 * (ai2p_1_over3 * ai2p_1_over3) : 0. ;

          // RC (122) for g2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

          (*Ma)(N * e + n_g2 , N * e + n_g2) -= (fac_RC122_2 > 0. ) ? fac_RC122_2 * eps_1_over3 / std::min(alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3 ,fac_sec) * 1./3./ std::min( alphag2_1_over3 * alphag2_1_over3,fac_sec)  *  ai2_5_over3 : 0. ;//alpha

          // WE (1)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

          (*Ma)(N * e + n_g1 , N * e + n_g1) -=(dfac_WE1da > 0. ) ? dfac_WE1da  * ai1 * ai1 : 0. ;

          // WE (2)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

          (*Ma)(N * e + n_g2 , N * e + n_g2) -=(fac_WE2 > 0. ) ? -fac_WE2 * (ai2 * ai2) / std::min(alpha(e, n_g2) * alpha(e, n_g2), fac_sec) : 0. ;


          //WE (122) for g2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

          (*Ma)(N * e + n_g2 , N * e + n_g2) -=(fac_WE122_2 > 0. ) ? - fac_WE122_2  * ai2 * alpha_p(e, n_g1) / std::min(alpha(e, n_g2) * alpha(e, n_g2),fac_sec) : 0. ;

        }
      if (Mai)//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        {
          // RC (1)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

          (*Mai)(N * e + n_g1 , N * e + n_g1) -=(fac_RC1 > 0. ) ? fac_RC1 * std::min(alpha(e, n_g1), alpha_sec) * 5./3. * std::cbrt(ai1) * std::cbrt(ai1) * eps_1_over3 : 0. ;

          // RC (2)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

          (*Mai)(N * e + n_g2 , N * e + n_g2) -=(fac_RC2 > 0. ) ? - fac_RC2 * (alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3) * 1./3.  / std::min(ai2_1_over3 * ai2_1_over3 * ai2_1_over3 * ai2_1_over3,fac_sec) * eps_1_over3 : 0. ;

          // RC (122) for g1---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

          (*Mai)(N * e + n_g1 , N * e + n_g1) -=(fac_RC122_1 > 0. ) ? fac_RC122_1 * eps_1_over3 * (alphag1_1_over3 * alphag1_1_over3) * (alphag2p_1_over3 * alphag2p_1_over3 * alphag2p_1_over3 * alphag2p_1_over3) * (ai2p_1_over3 * ai2p_1_over3) : 0. ;

          // RC (122) for g2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

          (*Mai)(N * e + n_g2 , N * e + n_g2) -=(fac_RC122_2 > 0. ) ? fac_RC122_2 * eps_1_over3 * (alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3) * alphag2_1_over3  * 5./3. * (ai2_1_over3 * ai2_1_over3) : 0. ;

          // WE (1)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

          (*Mai)(N * e + n_g1 , N * e + n_g1) -=(fac_WE1 > 0. ) ? fac_WE1 * 2. * ai1 : 0. ;
          // WE (2)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

          (*Mai)(N * e + n_g2 , N * e + n_g2) -=(fac_WE2 > 0. ) ? fac_WE2 * 2. * ai2 / std::max(alpha(e, n_g2),1./fac_sec) : 0. ;

          // WE (122) for g1---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

          (*Mai)(N * e + n_g1 , N * e + n_g1) -=(fac_WE122_1 > 0. ) ? fac_WE122_1  * ai2_p : 0. ;

          // WE (122) for g2---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

          (*Mai)(N * e + n_g2 , N * e + n_g2) -=(fac_WE122_2 > 0. ) ? fac_WE122_2  * 2. * ai2 * alpha_p(e, n_g1) / std::max(alpha(e, n_g2),1./fac_sec)  : 0. ;
        }
      if (Type_diss == "tau")//-------------------------------------------------------------------------------------------------------------------------------------------------------
        {
          if ((*tab_k)(e, n_l) * (*tau)(e, n_l) > limiter * nu_p(e, n_l)) // derivee en k ; depend de l'activation ou non du limiteur
            {
              if (Mk)//--------------------------------------------------------------------------------------------------------------------------------------------------------------------
                {
                  const double deps = 1./3. * std::cbrt(beta_k_) / std::min(std::cbrt((*tab_k)(e, n_l)) * std::cbrt((*tab_k)(e, n_l)),fac_sec) / std::min(std::cbrt((*tau)(e, n_l)),fac_sec) ;

                  // RC (1)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                  (*Mk)(N * e + n_g1, Nk * e + n_l) -=(fac_RC1 > 0. ) ?  fac_RC1 * std::min(alpha(e, n_g1), alpha_sec) * ai1_5_over3 * 1./3. * deps : 0. ;

                  // RC (2)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                  (*Mk)(N * e + n_g2, Nk * e + n_l)   -=(fac_RC2 > 0. ) ? fac_RC2 * (alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3) /  std::min(ai2_1_over3,fac_sec) * deps : 0.;

                  // RC (112)--------------------------------------------------------------------------------------------------------------------------------------------------------------

                  (*Mk)(N * e + n_g2, Nk * e + n_l) -=(fac_RC112 > 0. ) ?  fac_RC112 * std::min(alpha_p(e, n_g1), alpha_sec) * ai1p_5_over3 * deps: 0. ;

                  // RC (122) for g1----------------------------------------------------------------------------------------------------------------------------------------------------------------

                  (*Mk)(N * e + n_g1, Nk * e + n_l) -=(fac_RC122_1 > 0. ) ? fac_RC122_1 * (alphag1_1_over3 * alphag1_1_over3) * (alphag2p_1_over3 * alphag2p_1_over3 * alphag2p_1_over3 * alphag2p_1_over3) * ai1 * (ai2p_1_over3 * ai2p_1_over3) * deps : 0.;

                  // RC (122) for g2-------------------------------------------------------------------------------------------------------------------------------------------------------------

                  (*Mk)(N * e + n_g2, Nk * e + n_l) -=(fac_RC122_2 > 0. ) ? fac_RC122_2 * (alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3) * alphag2_1_over3  * (ai2_5_over3) * deps : 0.;

                }
              if (Mtau)//--------------------------------------------------------------------------------------------------------------------------------------------------------------------
                {
                  const double deps = -1./3. * std::cbrt(beta_k_) * std::cbrt((*tab_k)(e, n_l)) / std::min(std::cbrt((*tau)(e, n_l)) * std::cbrt((*tau)(e, n_l)) * std::cbrt((*tau)(e, n_l)) * std::cbrt((*tau)(e, n_l)),fac_sec) ;

                  // RC (1)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                  (*Mtau)(N * e + n_g1, Nk * e + n_l)-=(fac_RC1 > 0. ) ? fac_RC1 * std::min(alpha(e, n_g1), alpha_sec) * ai1_5_over3 * deps : 0.;

                  // RC (2)-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                  (*Mtau)(N * e + n_g2, Nk * e + n_l)-=(fac_RC2 > 0. ) ? fac_RC2 * (alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3) / std::min(ai2_1_over3,fac_sec) * deps : 0.;

                  // RC (112)----------------------------------------------------------------------------------------------------------------------------------------------------------------------

                  (*Mtau)(N * e + n_g2, Nk * e + n_l)-=(fac_RC112 > 0. ) ? fac_RC112 * std::min(alpha_p(e, n_g1), alpha_sec) * ai1p_5_over3 * deps : 0.;

                  // RC (122) for g1------------------------------------------------------------------------------------------------------------------------------------------------------------------

                  (*Mtau)(N * e + n_g1, Nk * e + n_l) -=(fac_RC122_1 > 0. ) ? fac_RC122_1 * (alphag1_1_over3 * alphag1_1_over3) * (alphag2p_1_over3 * alphag2p_1_over3 * alphag2p_1_over3 * alphag2p_1_over3) * ai1 * (ai2p_1_over3 * ai2p_1_over3) * deps : 0.;

                  // RC (122) for g2------------------------------------------------------------------------------------------------------------------------------------------------------------------

                  (*Mtau)(N * e + n_g2, Nk * e + n_l) -=(fac_RC122_2 > 0. ) ? fac_RC122_2 * (alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3) * alphag2_1_over3  * (ai2_5_over3) * deps : 0.;

                }
            }
          else if (Mk)//---------------------------------------------------------------------------------------------------------------------------------------------------------------------
            {
              const double deps = 1./3. * std::cbrt(beta_k_) / std::min(std::cbrt((*tab_k)(e, n_l)),fac_sec) / std::min(std::cbrt(limiter * nu_p(e, n_l)),fac_sec) ;

              // RC (1)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

              (*Mk)(N * e + n_g1, Nk * e + n_l)   -=(fac_RC1 > 0. ) ? fac_RC1 * std::min(alpha(e, n_g1), alpha_sec) * ai1_5_over3 * deps : 0.;

              // RC (2)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

              (*Mk)(N * e + n_g2, Nk * e + n_l)   -= (fac_RC2 > 0. ) ? fac_RC2 * (alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3) /  std::min(ai2_1_over3,fac_sec) * deps : 0.;

              // RC (112)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

              (*Mk)(N * e + n_g2, Nk * e + n_l)   -=(fac_RC112 > 0. ) ? fac_RC112 * std::min(alpha_p(e, n_g1), alpha_sec) * ai1p_5_over3 * deps : 0.;

              // RC (122) for g1----------------------------------------------------------------------------------------------------------------------------------------------------------------

              (*Mk)(N * e + n_g1, Nk * e + n_l) -=(fac_RC122_1 > 0. ) ? fac_RC122_1 * (alphag1_1_over3 * alphag1_1_over3) * (alphag2p_1_over3 * alphag2p_1_over3 * alphag2p_1_over3 * alphag2p_1_over3) * ai1 * (ai2p_1_over3 * ai2p_1_over3) * deps : 0.;

              //RC (122) for g2--------------------------------------------------------------------------------------------------------------------------------------------------------------------

              (*Mk)(N * e + n_g2, Nk * e + n_l) -=(fac_RC122_2 > 0. ) ? fac_RC122_2 * (alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3) * alphag2_1_over3  * (ai2_5_over3) * deps : 0.;

            }
        }
      if (Type_diss == "omega")//-------------------------------------------------------------------------------------------------------------------------------------------------------
        {
          if (Momega)//-----------------------------------------------------------------------------------------------------------------------------------------------------------------
            {
              const double deps = 1./3. * std::cbrt(beta_k_) * std::cbrt((*tab_k)(e, n_l)) / std::min(std::cbrt((*omega)(e, n_l)) * std::cbrt((*omega)(e, n_l)),fac_sec) ;

              // RC (1)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

              (*Momega)(N * e + n_g1 , Nk * e + n_l) -=(fac_RC1 > 0. ) ? fac_RC1 * std::min(alpha(e, n_g1), alpha_sec) * ai1_5_over3 * deps : 0.;

              // RC (2)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

              (*Momega)(N * e + n_g2 , Nk * e + n_l) -=(fac_RC2 > 0. ) ? fac_RC2 * (alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3) / std::min(ai2_1_over3,fac_sec) * deps : 0.;

              // RC (112)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

              (*Momega)(N * e + n_g2 , Nk * e + n_l) -=(fac_RC112 > 0. ) ? fac_RC112 * std::min(alpha_p(e, n_g1), alpha_sec) * ai1p_5_over3 * deps : 0.;

              // RC (122) for g1----------------------------------------------------------------------------------------------------------------------------------------------------------------------

              (*Momega)(N * e + n_g1, Nk * e + n_l) -=(fac_RC122_1 > 0. ) ? fac_RC122_1 * (alphag1_1_over3 * alphag1_1_over3) * (alphag2p_1_over3 * alphag2p_1_over3 * alphag2p_1_over3 * alphag2p_1_over3) * ai1 * (ai2p_1_over3 * ai2p_1_over3) * deps : 0.;

              // RC (122) for g2----------------------------------------------------------------------------------------------------------------------------------------------------------------

              (*Momega)(N * e + n_g2, Nk * e + n_l) -=(fac_RC1 > 0. ) ? fac_RC122_2 * (alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3) * alphag2_1_over3  * (ai2_5_over3) * deps : 0.;

            }
          if (Mk)//---------------------------------------------------------------------------------------------------------------------------------------------------------------------
            {
              const double deps = 1./3. * std::cbrt(beta_k_) / std::min(std::cbrt((*tab_k)(e, n_l)) * std::cbrt((*tab_k)(e, n_l)) ,fac_sec) * std::cbrt((*omega)(e, n_l)) ;

              // RC (1)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

              (*Mk)(N * e + n_g1 , Nk * e + n_l) -=(fac_RC1 > 0. ) ? fac_RC1 * std::min(alpha(e, n_g1), alpha_sec) * ai1_5_over3 * deps : 0.;

              // RC (2)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

              (*Mk)(N * e + n_g2 , Nk * e + n_l) -=(fac_RC2 > 0. ) ? fac_RC2 * (alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3) /  std::min(ai2_1_over3,fac_sec) * deps : 0.;

              // RC (112)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

              (*Mk)(N * e + n_g2 , Nk * e + n_l) -=(fac_RC112 > 0. ) ? fac_RC112 * std::min(alpha_p(e, n_g1), alpha_sec) * ai1p_5_over3 * deps : 0.;

              // RC (122) for g1---------------------------------------------------------------------------------------------------------------------------------------------------------------------

              (*Mk)(N * e + n_g1, Nk * e + n_l) -=(fac_RC122_1 > 0. ) ? fac_RC122_1 * (alphag1_1_over3 * alphag1_1_over3) * (alphag2p_1_over3 * alphag2p_1_over3 * alphag2p_1_over3 * alphag2p_1_over3) * ai1 * (ai2p_1_over3 * ai2p_1_over3) * deps : 0.;

              // RC (122) for g2--------------------------------------------------------------------------------------------------------------------------------------------------------------------

              (*Mk)(N * e + n_g2, Nk * e + n_l) -=(fac_RC122_2 > 0. ) ? fac_RC122_2 * (alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3 * alphag1_p_1_over3) * alphag2_1_over3  * (ai2_5_over3) * deps : 0.;

            }
        }
    }
}
