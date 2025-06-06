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
// File:        Rupture_Yao_Morel.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Rupture_bulles_1groupe_PolyMAC_P0.h>
#include <Pb_Multiphase.h>
#include <Champ_Elem_PolyMAC_P0.h>
#include <Milieu_composite.h>
#include <Op_Diff_Turbulent_PolyMAC_P0_Face.h>
#include <Viscosite_turbulente_base.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <Rupture_bulles_1groupe_base.h>
#include <Domaine_PolyMAC_P0.h>
#include <Champ_Face_base.h>
#include <cmath>
#include <math.h>

Implemente_instanciable(Rupture_bulles_1groupe_PolyMAC_P0, "Rupture_bulles_1groupe_elem_PolyMAC_P0", Source_base);

Sortie& Rupture_bulles_1groupe_PolyMAC_P0::printOn(Sortie& os) const
{
  return os;
}

Entree& Rupture_bulles_1groupe_PolyMAC_P0::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("beta_k", &beta_k_);
  param.lire_avec_accolades_depuis(is);


  const Pb_Multiphase *pbm = sub_type(Pb_Multiphase, equation().probleme()) ? &ref_cast(Pb_Multiphase, equation().probleme()) : nullptr;

  if (!pbm || pbm->nb_phases() == 1)
    Process::exit(que_suis_je() + " : not needed for single-phase flow!");

  for (int n = 0; n < pbm->nb_phases(); n++) //recherche de n_l, n_g : phase {liquide,gaz}_continu en priorite
    if (pbm->nom_phase(n).debute_par("liquide")
        && (n_l < 0 || pbm->nom_phase(n).finit_par("continu")))
      n_l = n;

  if (n_l < 0)
    Process::exit(que_suis_je() + " : liquid phase not found!");

  if (pbm->has_correlation("Rupture_bulles_1groupe"))
    correlation_ = pbm->get_correlation("Rupture_bulles_1groupe"); //correlation fournie par le bloc correlation
  else
    Correlation_base::typer_lire_correlation(correlation_, *pbm, "Rupture_bulles_1groupe", is); //sinon -> on la lit

  return is;
}

void Rupture_bulles_1groupe_PolyMAC_P0::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Domaine_PolyMAC_P0& domaine = ref_cast(Domaine_PolyMAC_P0, equation().domaine_dis());
  const int ne = domaine.nb_elem();
  const int ne_tot = domaine.nb_elem_tot();
  const int N = equation().inconnue().valeurs().line_size();

  for (auto &&n_m : matrices)
    if (n_m.first == "alpha" || n_m.first == "k" || n_m.first == "tau" || n_m.first == "omega" || n_m.first == "interfacial_area")
      {
        Matrice_Morse& mat = *n_m.second, mat2;
        const DoubleTab& dep = equation().probleme().get_champ(n_m.first.c_str()).valeurs();
        const int nc = dep.dimension_tot(0);
        const int M  = dep.line_size();
        IntTab sten(0, 2);

        if (n_m.first == "alpha")
          for (int e = 0; e < ne; e++)
            for (int n = 0; n < N; n++)
              {
                sten.append_line(N * e + n, N * e + n_l);
                if (n != n_l)
                  sten.append_line(N * e + n, N * e + n);
              }
        if (n_m.first == "interfacial_area" ) // N <= M
          for (int e = 0; e < ne; e++)
            for (int n = 0; n < N; n++) sten.append_line(N * e + n, M * e + n);
        if (n_m.first == "k" || n_m.first == "tau" || n_m.first == "omega") // N <= M
          for (int e = 0; e < ne; e++)
            for (int n = 0; n < N; n++)
              sten.append_line(N * e + n, M * e +n_l);

        Matrix_tools::allocate_morse_matrix(N * ne_tot, M * nc, sten, mat2);
        mat.nb_colonnes() ? mat += mat2 : mat = mat2;
      }
}

void Rupture_bulles_1groupe_PolyMAC_P0::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Domaine_PolyMAC_P0& domaine = ref_cast(Domaine_PolyMAC_P0, equation().domaine_dis());
  const DoubleVect& pe = equation().milieu().porosite_elem(), &ve = domaine.volumes();
  const Pb_Multiphase& pbm = ref_cast(Pb_Multiphase, equation().probleme());
  const DoubleTab& inco = equation().inconnue().valeurs();
  const DoubleTab& d_bulles = equation().probleme().get_champ("diametre_bulles").passe();
  const DoubleTab& alpha = pbm.equation_masse().inconnue().valeurs();
  const DoubleTab& alpha_p = pbm.equation_masse().inconnue().passe();
  const DoubleTab& press_p = ref_cast(QDM_Multiphase, pbm.equation_qdm()).pression().passe();
  const DoubleTab& temp_p  = pbm.equation_energie().inconnue().passe();
  const DoubleTab& rho_p   = equation().milieu().masse_volumique().passe();
  const DoubleTab& nu_p = equation().probleme().get_champ("viscosite_cinematique").passe();
  const DoubleTab *tab_k_p = equation().probleme().has_champ("k") ? &equation().probleme().get_champ("k").passe() : nullptr;
  const DoubleTab *tab_k = equation().probleme().has_champ("k") ? &equation().probleme().get_champ("k").valeurs() : nullptr;
  const DoubleTab *tau = equation().probleme().has_champ("tau") ? &equation().probleme().get_champ("tau").valeurs() : nullptr;
  const DoubleTab *omega = equation().probleme().has_champ("omega") ? &equation().probleme().get_champ("omega").valeurs() : nullptr ;

  const Milieu_composite& milc = ref_cast(Milieu_composite, equation().milieu());
  const int N = pbm.nb_phases();
  const int Nk = (tab_k) ? (*tab_k).line_size() : -1;
  const int Np = equation().probleme().get_champ("pression").valeurs().line_size();

  // Models use epsilon but with omega and tau we induce new variations of tau/omega and k
  std::string Type_diss = "other"; // omega, tau or other dissipation
  if (tau)
    Type_diss = "tau";
  else if (omega)
    Type_diss = "omega";

  DoubleTrav epsilon(alpha);
  const Op_Diff_Turbulent_PolyMAC_P0_Face& op_diff = ref_cast(Op_Diff_Turbulent_PolyMAC_P0_Face, equation().probleme().equation(0).operateur(0).l_op_base());
  const Viscosite_turbulente_base& visc_turb = ref_cast(Viscosite_turbulente_base, op_diff.correlation());
  visc_turb.eps(epsilon); // Epsilon is in the past
  const double limiter = visc_turb.limiteur();
  const double dh = 0;

  Matrice_Morse *Ma = matrices.count("alpha") ? matrices.at("alpha") : nullptr;
  Matrice_Morse *Mk = matrices.count("k") ? matrices.at("k") : nullptr;
  Matrice_Morse *Mtau = matrices.count("tau") ? matrices.at("tau") : nullptr;
  Matrice_Morse *Momega = matrices.count("omega") ? matrices.at("omega") : nullptr;
  Matrice_Morse *Mai = matrices.count("interfacial_area") ? matrices.at("interfacial_area") : nullptr;

  const int cR = (rho_p.dimension_tot(0) == 1);
  const int cM = (nu_p.dimension_tot(0) == 1);
  const int D = dimension;

  DoubleTrav a_l(N), p_l(N), T_l(N), rho_l(N), nu_l(N), sigma_l(N,N), dv(N, N), d_bulles_l(N), eps_l(Nk), k_l(Nk), coeff(N, N); //arguments pour coeff
  const Rupture_bulles_1groupe_base& correlation_rupt = ref_cast(Rupture_bulles_1groupe_base, correlation_.valeur());

  // fill velocity at elem tab
  DoubleTab pvit_elem(0, N * D);
  domaine.domaine().creer_tableau_elements(pvit_elem);
  const Champ_Face_base& ch_vit = ref_cast(Champ_Face_base,ref_cast(Pb_Multiphase, equation().probleme()).equation_qdm().inconnue());
  ch_vit.get_elem_vector_field(pvit_elem);

  const double fac_sec = 1.e4 ; // numerical security
  const double alpha_min = 1.e-3 ; // to avoid numerical problems

  /* elements */
  for (int e = 0; e < domaine.nb_elem(); e++)
    {

      for (int n = 0; n < N; n++)
        {
          a_l(n)   = alpha(e, n); // to further implicit the source term
          p_l(n)   = press_p(e, n * (Np > 1));
          T_l(n)   =  temp_p(e, n);
          rho_l(n) =   rho_p(!cR * e, n);
          nu_l(n)  =    nu_p(!cM * e, n);
          for (int k = 0; k < N; k++)
            if(milc.has_interface(n, k))
              {
                Interface_base& sat = milc.get_interface(n, k);
                sigma_l(n, k) = sat.sigma(temp_p(e,n), press_p(e,n * (Np > 1)));
              }
          d_bulles_l(n) = d_bulles(e,n);
        }

      for (int n = 0; n < Nk; n++)
        {
          eps_l(n) =epsilon(e, n) ;
          k_l(n)   = (tab_k_p)   ? (*tab_k_p)(e,n) : 0;
        }

      dv = 0;
      for (int d = 0; d < D; d++)
        for (int n = 0; n < N; n++)
          for (int k = 0; k < N; k++)
            dv(n, k) += (pvit_elem(e, N * d + n) - ((n!=k) ? pvit_elem(e, N * d + k) : 0) ) * (pvit_elem(e, N * d + n) - ((n!=k) ? pvit_elem(e, N * d + k) : 0) ); // nv(n,n) = ||v(n)||, nv(n, k!=n) = ||v(n)-v(k)||

      for (int n = 0; n < N; n++)
        for (int k = 0; k < N ; k++)
          dv(n, k) = sqrt(dv(n, k)) ;

      // Get correlations-------------------------------------------------------------------------------------------------------------------------------------------------------
      correlation_rupt.coefficient(a_l, p_l, T_l, rho_l, nu_l, sigma_l, dh, dv, d_bulles_l, eps_l, k_l, coeff); // Semi-Explicit coeff : alpha is implicit

      for (int k = 0; k < N ; k++)
        {
          double eps_valeurs {0.0};

          if (Type_diss == "tau")
            eps_valeurs = beta_k_ * ((*tab_k)(e, n_l)>1.e-8 ? (*tab_k)(e, n_l)*(*tab_k)(e, n_l)/ std::max((*tab_k)(e, n_l) * (*tau)(e, n_l), limiter * nu_p(e, n_l)) : 0 );
          else if (Type_diss == "omega")
            eps_valeurs = beta_k_ * ((*tab_k)(e, n_l)*(*omega)(e, n_l)) ;
          else
            eps_valeurs = epsilon(e, n_l);

          // Recuring functions for source terms

          const double cbrt_6 = std::cbrt(6.);
          const double factor_tmp = M_PI /( 3. * cbrt_6 * cbrt_6 * cbrt_6 * cbrt_6 * cbrt_6 );
          const double fac = (alpha(e, k)>alpha_min) ? pe(e) * ve(e) * factor_tmp * coeff(k, n_l) : 0. ;
          const double dalpha_fac = (alpha(e, k)>alpha_min) ? pe(e) * ve(e) * factor_tmp * coeff(n_l, k) : 0.;
          const double ai = std::max(inco(e, k), 0.) ; //security inco negative
          const double eps_1_over3 = std::cbrt(eps_valeurs) ;
          const double ai_5_over3 = std::cbrt(ai) * std::cbrt(ai) * std::cbrt(ai) * std::cbrt(ai) * std::cbrt(ai) ;

          const double cbrt_alpha_ek = std::cbrt(alpha(e, k));
          const double one_overalpha_2_over3 = (alpha(e, k) > alpha_min) ? 1./std::min(cbrt_alpha_ek * cbrt_alpha_ek, fac_sec) :0. ;
          const double one_overalpha_5_over3 = (alpha(e, k) > alpha_min) ? 1./std::min(cbrt_alpha_ek * cbrt_alpha_ek * cbrt_alpha_ek * cbrt_alpha_ek * cbrt_alpha_ek, fac_sec) : 0.;

          // Fill the matrix-------------------------------------------------------------------------------------------------------------------------------------------------------------

          secmem(e , k) += (fac > 0. ) ? fac * one_overalpha_2_over3 * alpha_p(e, n_l) * ai_5_over3 * eps_1_over3 : 0. ;// (alpha, ai, epsilon) implicit dependance

          if (Ma)//-----------------------------------------------------------------------------------------------------------------------------------------------------------------
            (*Ma)(N * e + k , N * e + k) -= (fac > 0. ) ? fac * -2./3.* one_overalpha_5_over3 * alpha_p(e, n_l) * ai_5_over3 * eps_1_over3 + dalpha_fac * one_overalpha_2_over3 * alpha_p(e, n_l) * ai_5_over3 * eps_1_over3 : 0. ;

          if (Mai)//----------------------------------------------------------------------------------------------------------------------------------------------------------------------
            (*Mai)(N * e + k , N * e + k) -= (fac > 0. ) ? fac * one_overalpha_2_over3 * alpha_p(e, n_l) * 5./3. * std::cbrt(ai) * std::cbrt(ai) * eps_1_over3 : 0. ;

          // Turbulent dissipation
          if (Type_diss == "tau")
            {
              if ((*tab_k)(e, n_l) * (*tau)(e, n_l) > limiter * nu_p(e, n_l)) // derivative according to k due to epsilon
                {
                  if (Mk)
                    (*Mk)(N * e + k, Nk * e + n_l)   -= (fac > 0. ) ? fac * one_overalpha_2_over3 * alpha_p(e, n_l) * ai_5_over3 * 1./3. * std::cbrt(beta_k_)  / std::min(std::cbrt((*tab_k)(e, n_l)) * std::cbrt((*tab_k)(e, n_l)),fac_sec) / std::min(std::cbrt((*tau)(e, n_l)),fac_sec) : 0.;

                  if (Mtau)
                    (*Mtau)(N * e + k, Nk * e + n_l) -= (fac > 0. ) ? fac * one_overalpha_2_over3 * alpha_p(e, n_l) * ai_5_over3 *-1./3. * std::cbrt(beta_k_)  * std::cbrt((*tab_k)(e, n_l)) / std::min(std::cbrt((*tau)(e, n_l)) * std::cbrt((*tau)(e, n_l)) * std::cbrt((*tau)(e, n_l)) * std::cbrt((*tau)(e, n_l)),fac_sec) : 0. ;
                }

              else
                {
                  if (Mk)//--------------------------------------------------------------------------------------------------------------------------------------------------------
                    (*Mk)(N * e + k, Nk * e + n_l) -= (fac > 0. ) ? fac * one_overalpha_2_over3 * alpha_p(e, n_l) * ai_5_over3 * 1./3. * std::cbrt(beta_k_) / std::min(std::cbrt((*tab_k)(e, n_l)),fac_sec) / std::min(std::cbrt(limiter * nu_p(e, n_l)),fac_sec) : 0.;
                }
            }
          if (Type_diss == "omega")
            {
              if (Momega)
                (*Momega)(N * e + k , Nk * e + n_l) -= (fac > 0. ) ? fac * one_overalpha_2_over3 * alpha_p(e, n_l) * ai_5_over3 * 1./3. * std::cbrt(beta_k_) * std::cbrt((*tab_k)(e, n_l)) / std::min(std::cbrt((*omega)(e, n_l)) * std::cbrt((*omega)(e, n_l)),fac_sec) : 0.;

              if (Mk)
                (*Mk)(N * e + k , Nk * e + n_l) -= (fac > 0. ) ? fac * one_overalpha_2_over3 * alpha_p(e, n_l) * ai_5_over3 * 1./3. * std::cbrt(beta_k_) / std::min(std::cbrt((*tab_k)(e, n_l)) * std::cbrt((*tab_k)(e, n_l)) ,fac_sec) * std::cbrt((*omega)(e, n_l)) : 0.;
            }

        }
    }
}
