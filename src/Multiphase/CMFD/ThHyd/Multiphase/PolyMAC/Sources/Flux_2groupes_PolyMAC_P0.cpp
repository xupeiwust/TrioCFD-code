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
////////////////////////////////////////////////////////////////////////////////

///Joseph L. Bottini, Taiyang Zhang, Caleb S. Brooks, Validation of two-group interfacial area transport equation in boiling flow,
/// International Journal of Heat and Mass Transfer, 2024, https://doi.org/10.1016/j.ijheatmasstransfer.2024.125515

////////////////////////////////////////////////////////////////////////////////

#include <Op_Diff_Turbulent_PolyMAC_P0_Face.h>
#include <Domaine_Cl_PolyMAC.h>
#include <Neumann_paroi.h>
#include <Flux_2groupes_PolyMAC_P0.h>
#include <Viscosite_turbulente_base.h>
#include <Aire_interfaciale.h>
#include <Flux_2groupes_base.h>
#include <Milieu_composite.h>
#include <Pb_Multiphase.h>
#include <Array_tools.h>
#include <Matrix_tools.h>
#include <math.h>
#include <Changement_phase_base.h>
#include <Sources.h>
#include <Flux_parietal_base.h>
#include <Flux_interfacial_PolyMAC_P0P1NC.h>
#include <Champ_Elem_PolyMAC_P0.h>
#include <Champ_Face_PolyMAC_P0.h>
#include <Champ_Face_base.h>
#include <Domaine_PolyMAC_P0.h>

Implemente_instanciable(Flux_2groupes_PolyMAC_P0,"Flux_2groupes_elem_PolyMAC_P0", Source_base);

Sortie& Flux_2groupes_PolyMAC_P0::printOn(Sortie& os) const { return os; }

Entree& Flux_2groupes_PolyMAC_P0::readOn(Entree& is)
{

  Param param(que_suis_je());
  param.ajouter("Dh", &dh_);
  param.lire_avec_accolades_depuis(is);

  const Pb_Multiphase *pbm = sub_type(Pb_Multiphase, equation().probleme()) ? &ref_cast(Pb_Multiphase, equation().probleme()) : nullptr;

  for (int n = 0; n < pbm->nb_phases(); n++)  //recherche de n_l, n_g : phase {liquide,gaz}_continu en priorite
    {
      if (pbm->nom_phase(n).debute_par("liquide") && (n_l < 0 || pbm->nom_phase(n).finit_par("continu")))  n_l = n;
      if (( pbm->nom_phase(n).finit_par("group1")))  n_g1 = n;
      if (( pbm->nom_phase(n).finit_par("group2")))  n_g2 = n;
    }
  if (n_l < 0) Process::exit(que_suis_je() + " : liquid phase not found!");
  if (n_g1 < 0) Process::exit(que_suis_je() + " : group 1 not found!");
  if (n_g2 < 0) Process::exit(que_suis_je() + " : group 2 not found!");

  if (pbm->has_correlation("flux_2groupes")) correlation_ = pbm->get_correlation("flux_2groupes"); //correlation fournie par le bloc correlation
  else Correlation_base::typer_lire_correlation(correlation_, *pbm, "flux_2groupes", is); //sinon -> on la lit

  const Sources& les_sources_loc = pbm->equation(2).sources();
  for (int j = 0 ; j<les_sources_loc.size(); j++)
    {
      if sub_type(Flux_interfacial_PolyMAC_P0P1NC, les_sources_loc(j).valeur()) src_flux_interfacial_ = les_sources_loc(j).valeur();
    }

  return is;
}

void Flux_2groupes_PolyMAC_P0::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Domaine_PolyMAC_P0& domaine = ref_cast(Domaine_PolyMAC_P0, equation().domaine_dis());
  const int ne = domaine.nb_elem(), ne_tot = domaine.nb_elem_tot(), N = equation().inconnue().valeurs().line_size();

  for (auto &&n_m : matrices)
    if (n_m.first == "temperature" ||  n_m.first == "pression" || n_m.first == "alpha"  ||  n_m.first == "interfacial_area" )
      {
        Matrice_Morse& mat = *n_m.second, mat2;
        const DoubleTab& dep = equation().probleme().get_champ(n_m.first.c_str()).valeurs();
        int nc = dep.dimension_tot(0),
            M  = dep.line_size();
        IntTab sten(0, 2);
        if (n_m.first == "alpha")
          for (int e = 0; e < ne; e++)
            for (int n = 0; n < N; n++) sten.append_line(N * e + n, N * e + n);
        if (n_m.first == "interfacial_area" || n_m.first == "temperature") // N <= M
          for (int e = 0; e < ne; e++)
            for (int n = 0; n < N; n++) sten.append_line(N * e + n, M * e + n);
        if (n_m.first == "pression" )
          for (int e = 0; e < ne; e++)
            for (int n = 0, m = 0; n < N; n++, m+=(M>1)) sten.append_line(N * e + n, M * e + m);

        //tableau_trier_retirer_doublons(sten);
        Matrix_tools::allocate_morse_matrix(N * ne_tot, M * nc, sten, mat2);
        mat.nb_colonnes() ? mat += mat2 : mat = mat2;
      }

}


void Flux_2groupes_PolyMAC_P0::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Pb_Multiphase& pbm = ref_cast(Pb_Multiphase, equation().probleme());
  const Milieu_composite& milc = ref_cast(Milieu_composite, equation().milieu());
  const Domaine_PolyMAC_P0& domaine = ref_cast(Domaine_PolyMAC_P0, equation().domaine_dis());
  const DoubleVect& pe = milc.porosite_elem(), &ve = domaine.volumes();
  const Champ_base& ch_rho = milc.masse_volumique();
  const Champ_Elem_PolyMAC_P0& ch = ref_cast(Champ_Elem_PolyMAC_P0, pbm.equation_energie().inconnue());
  const Champ_Inc_base *pch_rho = sub_type(Champ_Inc_base, ch_rho) ? &ref_cast(Champ_Inc_base, ch_rho) : nullptr;
  const IntTab& el_f = domaine.elem_faces(), &fcl = ch.fcl(), &f_e = domaine.face_voisins();
  const Conds_lim& cls = ch.domaine_Cl_dis().les_conditions_limites();

  const DoubleTab& inco = equation().inconnue().valeurs(),
                   &alpha = pbm.equation_masse().inconnue().valeurs(),
                    &alpha_p = pbm.equation_masse().inconnue().passe(),
                     &ai = equation().probleme().get_champ("interfacial_area").passe(),
                      &press = ref_cast(QDM_Multiphase,pbm.equation_qdm()).pression().valeurs(),
                       &press_p = ref_cast(QDM_Multiphase,pbm.equation_qdm()).pression().passe(),
                        &temp_p  = pbm.equation_energie().inconnue().passe(),
                         &rho = equation().milieu().masse_volumique().valeurs(),
                          &rho_p = equation().milieu().masse_volumique().passe(),
                           &d_b_p = equation().probleme().get_champ("diametre_bulles").passe(),
                            *k_turb = (equation().probleme().has_champ("k")) ? &equation().probleme().get_champ("k").passe() : nullptr ,
                             &lambda = milc.conductivite().passe(),
                              &Cp = milc.capacite_calorifique().passe(),
                               &mu_p = milc.viscosite_dynamique().passe();


  Matrice_Morse *Mp = matrices.count("pression")    ? matrices.at("pression")    : nullptr,
                 *Mt = matrices.count("temperature") ? matrices.at("temperature") : nullptr,
                  *Ma = matrices.count("alpha") ? matrices.at("alpha") : nullptr,
                   *Mai = matrices.count("interfacial_area") ? matrices.at("interfacial_area") : nullptr;


  int N = inco.line_size(), Np = press.line_size(), Nk = (k_turb) ? (*k_turb).line_size() : -1, D = dimension, n, k, e,d;
  const int cR = (rho_p.dimension_tot(0) == 1), cM = (mu_p.dimension_tot(0) == 1), cL = (lambda.dimension_tot(0) == 1), cCp = (Cp.dimension_tot(0) == 1) ;

  const Flux_2groupes_base& correlation_flux = ref_cast(Flux_2groupes_base, correlation_.valeur());
  DoubleTrav nut(domaine.nb_elem_tot(), N);

  double pas_tps = equation().probleme().schema_temps().pas_de_temps();


  DoubleTrav epsilon(alpha);
  const Op_Diff_Turbulent_PolyMAC_P0_Face& op_diff = ref_cast(Op_Diff_Turbulent_PolyMAC_P0_Face, equation().probleme().equation(0).operateur(0).l_op_base());
  const Viscosite_turbulente_base&   visc_turb = ref_cast(Viscosite_turbulente_base, op_diff.correlation());
  visc_turb.eps(epsilon);

  Flux_2groupes_base::input_coeffs in_coeffs;
  Flux_2groupes_base::output_coeffs out_coeffs;
  Flux_2groupes_base::input_therms in_therms;
  Flux_2groupes_base::output_therms out_therms;

  in_coeffs.n_g1 = n_g1, in_therms.n_g1 = n_g1;
  in_coeffs.n_g2=n_g2, in_therms.n_g2=n_g2;
  in_coeffs.n_l=n_l, in_therms.n_l=n_l;

  in_coeffs.alpha.resize(N), in_coeffs.p.resize(N) , in_coeffs.rho.resize(N) , in_coeffs.mu.resize(N) , in_coeffs.d_bulles.resize(N) , in_coeffs.sigma.resize(N, N) , in_coeffs.epsilon.resize(Nk) , in_coeffs.k_turb.resize(Nk) , in_coeffs.nv.resize(N, N), in_coeffs.a_i.resize(N) ;
  in_therms.alpha.resize(N), in_therms.T.resize(N), in_therms.p.resize(N) , in_therms.d_bulles.resize(N), in_therms.lambda.resize(N), in_therms.mu.resize(N), in_therms.rho.resize(N), in_therms.Cp.resize(N), in_therms.sigma.resize(N, N);

  out_coeffs.gamma.resize(N, N), out_coeffs.da_gamma.resize(N, N),out_coeffs.dai_gamma.resize(N, N);
  out_therms.dT_G1.resize(N), out_therms.dT_G2.resize(N) , out_therms.da_G1.resize(N), out_therms.da_G2.resize(N) , out_therms.dT_etaph1.resize(N) , out_therms.dT_etaph2.resize(N) , out_therms.da_etaph1.resize(N), out_therms.da_etaph2.resize(N) ;



  // fill velocity at elem tab
  DoubleTab pvit_elem(0, N * D);
  domaine.domaine().creer_tableau_elements(pvit_elem);
  const Champ_Face_base& ch_vit = ref_cast(Champ_Face_base,ref_cast(Pb_Multiphase, equation().probleme()).equation_qdm().inconnue());
  ch_vit.get_elem_vector_field(pvit_elem);

  const double fac_sec = 1.e4; // numerical security factor


  for (e = 0; e < domaine.nb_elem(); e++)
    {
      // Get field values for correlations-------------------------------------------------------------------------------------------------------------------------------------------------------

      for (in_coeffs.nv =0, d = 0; d < D; d++)
        for ( n = 0; n < N; n++)
          for (k = 0 ; k<N ; k++) in_coeffs.nv(n, k) += (pvit_elem(e, N * d + n) - ((n!=k) ? pvit_elem(e, N * d + k) : 0)) * (pvit_elem(e, N * d + n) - ((n!=k) ? pvit_elem(e, N * d + k) : 0)); // nv(n,n) = ||v(n)||, nv(n, k!=n) = ||v(n)-v(k)||
      for (n = 0; n < N; n++)
        for ( k = 0 ; k<N ; k++) in_coeffs.nv(n, k) = sqrt(in_coeffs.nv(n, k)) ;
      /* arguments de coeff */
      for (n = 0; n < N; n++)
        {
          in_coeffs.alpha[n] = alpha_p(e, n), in_therms.alpha[n] = alpha_p(e, n) ;
          in_coeffs.p[n]     = press_p(e, n * (Np > 1)), in_therms.p[n]     = press_p(e, n * (Np > 1));
          in_therms.T[n]   =  temp_p(e, n);
          in_therms.lambda[n] = lambda(!cL * e, n) ;
          in_therms.Cp[n] = Cp(!cCp * e, n) ;
          in_coeffs.rho[n]   = rho_p(!cR * e, n), in_therms.rho[n]   = rho_p(!cR * e, n) ;
          in_coeffs.mu[n]   = mu_p(!cM * e, n), in_therms.mu[n]   = mu_p(!cM * e, n);
          in_coeffs.d_bulles[n] = d_b_p(e,n), in_therms.d_bulles[n] = d_b_p(e,n) ;
          in_coeffs.a_i[n] = std::max(ai(e,n),0.) ;
          for (k = 0; k < N; k++)
            {
              if(milc.has_interface(n, k))
                {
                  Interface_base& sat = milc.get_interface(n, k);
                  in_coeffs.sigma(n,k) = sat.sigma(temp_p(e,n), press_p(e,n * (Np > 1)));
                }
              else if (milc.has_saturation(n, k))
                {
                  Saturation_base& z_sat = milc.get_saturation(n, k);

                  DoubleTab& sig = z_sat.get_sigma_tab();
                  in_coeffs.sigma(n,k) = sig(e);
                  in_therms.sigma(n,k) = sig(e);

                }
            }
        }
      for (n = 0; n < Nk; n++)
        {
          in_coeffs.epsilon[n] = epsilon(e, n) ;
          in_coeffs.k_turb[n]   = (k_turb)   ? (*k_turb)(e,n) : 0;
        }

      // Initialisation-----------------------------------------------------------------------------------------------------

      double bilaninterg1 = 0. ;
      double bilaninterg2 = 0.  ;
      double dag1_bilaninterg1 =  0. ;
      double dag2_bilaninterg2 =   0. ;


      if (milc.has_saturation(n_l, n_g1)) // If eos exit we can compute thermal effects
        {

          Saturation_base& sat1 = milc.get_saturation(n_l, n_g1);
          Saturation_base& sat2 = milc.get_saturation(n_l, n_g2);

          in_therms.Lvap = sat1.Lvap(press_p(e,n_l * (Np > 1))) ;
          in_therms.hlsat = sat1.Hls(press_p(e,n_l * (Np > 1)));
          in_therms.Tsatg1= sat1.Tsat(   press_p(e,n_l * (Np > 1))) ;
          in_therms.Tsatg2 = sat2.Tsat(press_p(e,n_l * (Np > 1))) ;
          in_therms.dp_Tsat1= sat1.dP_Tsat(press_p(e,n_l * (Np > 1))) ;
          in_therms.dp_Tsat2 = sat2.dP_Tsat(press_p(e,n_l * (Np > 1))) ;
          in_therms.qp = 0. ;
          for (int j = 0; j < int( el_f.dimension(1) ); j++)
            {
              int fb = el_f(e, j);
              if (fb < 0) continue;

              int el = f_e(fb, 1);
              if (el < 0)
                {
                  if (fcl(fb,0) == 4)
                    {
                      in_therms.qp =  ref_cast(Neumann_paroi, cls[fcl(fb, 1)].valeur()).flux_impose(fcl(fb, 2), 0); // Get wall heat
                    }
                }
            }

          // Get correlations for thermal effects-----------------------------------------------------------------------------------------------------------------------------------------------------

          correlation_flux.therm(in_therms, out_therms);

          // Thermal exchange balance for each group-------------------------------------------------------------------------------------------------------------------------------------------------

          bilaninterg1 = out_therms.G1/rho(e, n_g1) - out_therms.etaph1  - alpha(e, n_g1) * ( 1. - rho_p(e, n_g1)/rho(e, n_g1))/pas_tps; // for group 1
          bilaninterg2 = out_therms.G2/rho(e, n_g1) - out_therms.etaph2 - alpha(e, n_g2) * ( 1. - rho_p(e, n_g2)/rho(e, n_g2))/pas_tps  ; // for group 2
          dag1_bilaninterg1 = (bilaninterg1> 0.) ? out_therms.da_G1(n_g1)/rho(e, n_g1) - out_therms.da_etaph1(n_g1)  - ( 1. - rho_p(e, n_g1)/rho(e, n_g1))/pas_tps : 0. ; // void fraction variations
          dag2_bilaninterg2 = (0.> bilaninterg2) ? out_therms.da_G2(n_g2)/rho(e, n_g1) - out_therms.da_etaph2(n_g2) - ( 1. - rho_p(e, n_g2)/rho(e, n_g2))/pas_tps : 0. ; // void fraction variations


        }



      const double vol = pe(e) * ve(e) ;

      in_coeffs.dh = dh_, in_coeffs.e = e, in_therms.dh = dh_  ;

      // Get correlations for non thermal effects-----------------------------------------------------------------------------------------------------------------------------------------------

      correlation_flux.coeffs(in_coeffs, out_coeffs);


// Mass balance equation-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      if (sub_type(Masse_Multiphase, equation()))
        {

          // Interphase transfer between group 1 and group 2 : G 1->2 = G 2->1
          // For group1

          const double gammag1 = - rho(e, n_g1) * (out_coeffs.gamma(n_g1,n_g1) + out_coeffs.inter3g1 * std::max(bilaninterg1,0.) + out_coeffs.inter3g2 * std::min(bilaninterg2,0.) ) ;

          // Physical/Numerical security because correlations are not asymptotically valid : 1 is ok, 0 is cancelled

          double limg1 =  1.;
          if (alpha(e, n_g1) < 1./fac_sec) // if not enough group 1
            {

              limg1 = (gammag1 < 0.) ? 0. : 1.;
            }

          // For group2

          const double gammag2 = -gammag1  ;

          // Physical/Numerical security because correlations are not asymptotically valid : 1 is ok, 0 is cancelled

          double limg2 = 1. ;

          if (alpha(e, n_g2) < 1./fac_sec) // if not enough group 2
            {

              limg2 = (gammag2 < 0.) ? 0. : 1.;
            }


          // Fill matrix with mass transfer---------------------------------------------------------------------------------------------------------------------------------------------------------

          // Non-thermal mass tranfers : only between groups----------------------------------------------------------------------------------------------------------------------------

          secmem(e, n_g1) +=  limg1 *  vol * gammag1  ;// (alpha1, ai1, T, P) implicit dependance

          secmem(e, n_g2) +=  limg2 *  vol * gammag2 ;// (alpha2, ai2, T, P) implicit dependance

          if (milc.has_saturation(n_l, n_g1))// If eos exit we can compute thermal effects----------------------------------------------------------------------------------------------
            {
              // Thermal mass transfers between gas and liquid

              secmem(e, n_g1) +=  limg1 *  vol * out_therms.G1  ;// (alpha1, T, P) implicit dependance

              secmem(e, n_g2) +=  limg2 *  vol * out_therms.G2  ;// (alpha2, T, P) implicit dependance

              secmem(e, n_l) += - limg1 *  vol * out_therms.G1 - limg2 *  vol * out_therms.G2 ;// (alpha_l, T, P) implicit dependance

            }

          if (Ma) // void fraction dependance-------------------------------------------------------------------------------------------------------------------------------------------------------
            {
              // Non-thermal mass tranfers : only between groups----------------------------------------------------------------------------------------------------------------------------

              (*Ma)(N * e + n_g1, N * e + n_g1) -= -limg1 * vol * rho(e, n_g1) * (  out_coeffs.da_gamma(n_g1,n_g1) + out_coeffs.inter3g1 * dag1_bilaninterg1 + out_coeffs.da_inter3g1 * std::max(bilaninterg1,0.) );

              (*Ma)(N * e + n_g2, N * e + n_g2) -=  limg2 * vol * rho(e, n_g1) * ( out_coeffs.da_gamma(n_g2,n_g2) + out_coeffs.inter3g2 * dag2_bilaninterg2 + out_coeffs.da_inter3g2 * std::max(bilaninterg2,0.)  );


              if (milc.has_saturation(n_l, n_g1)) // If eos exit we can compute thermal effects----------------------------------------------------------------------------------------------
                {

                  (*Ma)(N * e + n_g1, N * e + n_g1) -= limg1 *  vol * out_therms.da_G1(n_g1) ;

                  (*Ma)(N * e + n_g2, N * e + n_g2) -= limg2 *  vol * out_therms.da_G2(n_g2) ;

                  (*Ma)(N * e + n_l, N * e + n_l) -= - limg1 *  vol * out_therms.da_G1(n_l) - limg2 *  vol * out_therms.da_G2(n_l) ;

                }
            }

          if (Mai) // interfacial area dependance only for non thermal effects--------------------------------------------------------------------------------------------------------------------
            {
              // Non-thermal mass tranfers : only between groups----------------------------------------------------------------------------------------------------------------------------

              (*Mai)(N * e + n_g1, N * e + n_g1) -=  - limg1 * vol * rho(e, n_g1) * out_coeffs.dai_gamma(n_g1,n_g1);

              (*Mai)(N * e + n_g2, N * e + n_g2) -= limg2 * vol * rho(e, n_g1) * out_coeffs.dai_gamma(n_g2,n_g2);

            }

          // (T,P) dependance only for thermal effects---------------------------------------------------------------------------------------------------------------------------------------------

          if (milc.has_saturation(n_l, n_g1)) // If eos exit we can compute thermal effects----------------------------------------------------------------------------------------------
            {
              // dT : Temperature variations, dp : Pressure variations

              const double dT_bilaninterg1 = (bilaninterg1> 0.) ? out_therms.dT_G1(n_g1)/rho(e, n_g1)- out_therms.G1/rho(e, n_g1)/rho(e, n_g1) * pch_rho->derivees().at("temperature")(e, n_g1) - out_therms.dT_etaph1(n_g1)  - alpha(e, n_g1) * ( rho_p(e, n_g1)/rho(e, n_g1)/rho(e, n_g1))* pch_rho->derivees().at("temperature")(e, n_g1)/pas_tps : 0.;

              const double dp_bilaninterg1 = (bilaninterg1> 0.) ? out_therms.dp_G1/rho(e, n_g1)- out_therms.G1/rho(e, n_g1)/rho(e, n_g1) * pch_rho->derivees().at("pression")(e, n_g1) - out_therms.dp_etaph1  - alpha(e, n_g1) * ( rho_p(e, n_g1)/rho(e, n_g1)/rho(e, n_g1))* pch_rho->derivees().at("pression")(e, n_g1)/pas_tps : 0.;

              const double dT_bilaninterg2 = (bilaninterg2> 0.) ? out_therms.dT_G2(n_g2)/rho(e, n_g1) - out_therms.G2/rho(e, n_g2)/rho(e, n_g2) * pch_rho->derivees().at("temperature")(e, n_g2) - out_therms.dT_etaph2(n_g2) - alpha(e, n_g2) * ( rho_p(e, n_g2)/rho(e, n_g2)/rho(e, n_g2))*pch_rho->derivees().at("temperature")(e, n_g2)/pas_tps : 0. ;

              const double dp_bilaninterg2 = (bilaninterg2> 0.) ?  out_therms.dp_G2/rho(e, n_g1) - out_therms.G2/rho(e, n_g2)/rho(e, n_g2) * pch_rho->derivees().at("temperature")(e, n_g2) - out_therms.dp_etaph2 - alpha(e, n_g2) * ( rho_p(e, n_g2)/rho(e, n_g2)/rho(e, n_g2))*pch_rho->derivees().at("pression")(e, n_g2)/pas_tps : 0.;

              if (Mt)//--------------------------------------------------------------------------------------------------------------------------------------------------------
                {

                  // Intergoup thermal effects-------------------------------------------------------------------------------------------

                  (*Mt)(N * e + n_g1, N * e + n_g1) -=   -limg1 * vol * ( pch_rho->derivees().at("temperature")(e, n_g1) * (  out_coeffs.gamma(n_g1,n_g1) + out_coeffs.inter3g1 * std::max(bilaninterg1,0.) - out_coeffs.inter3g2 * std::min(bilaninterg2,0.) ) +  rho(e, n_g1) * (  out_coeffs.inter3g1 * dT_bilaninterg1 + out_coeffs.inter3g2 * dT_bilaninterg2 ));

                  (*Mt)(N * e + n_g2, N * e + n_g2) -=  limg2 * vol * ( pch_rho->derivees().at("temperature")(e, n_g1) * ( out_coeffs.gamma(n_g1,n_g1) + out_coeffs.inter3g1 * std::max(bilaninterg1,0.) + out_coeffs.inter3g2 * std::min(bilaninterg2,0.) ) + rho(e, n_g1) * (  out_coeffs.inter3g1 * dT_bilaninterg1 + out_coeffs.inter3g2 * dT_bilaninterg2 ) );

                  // Exchanges gas-liquid---------------------------------------------------------------------------------------------------

                  (*Mt)(N * e + n_g1, N * e + n_g1) -= limg1 *  vol * out_therms.dT_G1(n_g1) ;

                  (*Mt)(N * e + n_g2, N * e + n_g2) -= limg2 *  vol * out_therms.dT_G2(n_g2) ;

                  (*Mt)(N * e + n_l, N * e + n_l) -= - limg1 *  vol * out_therms.dT_G1(n_l) - limg2 *  vol * out_therms.dT_G2(n_l) ;

                }
              if (Mp)//--------------------------------------------------------------------------------------------------------------------------------------------------------------
                {

                  // Intergoup thermal effects-------------------------------------------------------------------------------------------

                  (*Mp)(N * e + n_g1, e) -= -limg1 * vol * ( pch_rho->derivees().at("pression")(e, n_g1) * (  out_coeffs.gamma(n_g1,n_g1) + out_coeffs.inter3g1 * std::max(bilaninterg1,0.) - out_coeffs.inter3g2 * std::min(bilaninterg2,0.) ) + rho(e, n_g1) * (   out_coeffs.inter3g1 * dp_bilaninterg1 + out_coeffs.inter3g2 * dp_bilaninterg2 ) );

                  (*Mp)(N * e + n_g2, e) -=  limg2 * vol * ( pch_rho->derivees().at("pression")(e, n_g1) * ( out_coeffs.gamma(n_g1,n_g1) + out_coeffs.inter3g1 * std::max(bilaninterg1,0.) + out_coeffs.inter3g2 * std::min(bilaninterg2,0.) ) +  rho(e, n_g1) * (  out_coeffs.inter3g1 * dp_bilaninterg1 + out_coeffs.inter3g2 * dp_bilaninterg2 ) );

                  // Exchanges gas-liquid---------------------------------------------------------------------------------------------------

                  (*Mp)(N * e + n_g1, e) -= limg1 *  vol * out_therms.dp_G1 ;

                  (*Mp)(N * e + n_g2, e) -= limg2 *  vol * out_therms.dp_G2 ;

                  (*Mp)(N * e + n_l, e) -= - limg1 *  vol * out_therms.dp_G1 - limg2 *  vol * out_therms.dp_G2 ;

                }
            }
        }

// IATE -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      else if (sub_type(Aire_interfaciale, equation()))
        {

          const double ai1 = std::max(inco(e, n_g1), 0.) ; // security inco negative
          const double ai2 = std::max(inco(e, n_g2), 0.) ; // security inco negative

          // Recurring term with numerical issue when void fraction is 0 : ai/alpha

          const double aisural_1 = (1./fac_sec < alpha(e, n_g1) ) ? ai1 / std::max(alpha(e, n_g1),1./fac_sec) : 0. ;
          const double dal_aisural_1 = (0. < aisural_1 ) ? - ai1 / std::max(alpha(e, n_g1) /alpha(e, n_g1),1./fac_sec) : 0. ;
          const double dai_aisural_1 = (0. < aisural_1 ) ? 1. / std::max(alpha(e, n_g1),1./fac_sec) : 0. ;
          const double aisural_2 = (1./fac_sec < alpha(e, n_g2) ) ? ai2 / std::max(alpha(e, n_g2),1./fac_sec) : 0. ;
          const double dal_aisural_2 = (0. < aisural_2 ) ? - ai2 / std::max(alpha(e, n_g2) /alpha(e, n_g2),1./fac_sec) : 0. ;
          const double dai_aisural_2 = (0. < aisural_2 ) ? 1. / std::max(alpha(e, n_g2), 1./fac_sec) : 0. ;

          // Fill matrix with transfers---------------------------------------------------------------------------------------------------------------------------------------------------------

          secmem(e, n_g1) +=   vol * ( std::max(bilaninterg1,0.) * aisural_1 * (2./3. - out_coeffs.inter2g1) - std::min(bilaninterg2,0.) * aisural_2 * out_coeffs.inter2g2 + aisural_1 * out_therms.etaph1 )  ;

          secmem(e, n_g2) +=   vol * (std::max(bilaninterg2,0.) * aisural_2 * (2./3. + out_coeffs.inter2g2) + std::min(bilaninterg1,0.) * aisural_1 * out_coeffs.inter2g1 + aisural_2 * out_therms.etaph2 )  ;


          if (Ma)//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
            {

              (*Ma)(N * e + n_g1 , N * e + n_g1) -=   vol * ( std::max(bilaninterg1,0.) * dal_aisural_1 * (2./3. - out_coeffs.inter2g1) + dag1_bilaninterg1 * aisural_1 * (2./3. - out_coeffs.inter2g1) - std::max(bilaninterg1,0.) * aisural_1 * out_coeffs.da_inter2g1 + dal_aisural_1 * out_therms.etaph1 + aisural_1 * out_therms.da_etaph1(n_g1) ) ;

              (*Ma)(N * e + n_g2 , N * e + n_g2) -=   vol * ( std::max(bilaninterg2,0.) * dal_aisural_2 * (2./3. - out_coeffs.inter2g2) +  dag2_bilaninterg2 * aisural_2 * (2./3. + out_coeffs.inter2g2) +  std::max(bilaninterg2,0.) * aisural_2 * out_coeffs.da_inter2g2 + dal_aisural_2 * out_therms.etaph2 + aisural_2 * out_therms.da_etaph2(n_g2) )  ;

            }

          if (Mai)//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
            {

              (*Mai)(N * e + n_g1 , N * e + n_g1) -= vol * ( std::max(bilaninterg1,0.) * dai_aisural_1 * (2./3. - out_coeffs.inter2g1) + dai_aisural_1 * out_therms.etaph1 )  ;

              (*Mai)(N * e + n_g2 , N * e + n_g2) -= vol * ( std::max(bilaninterg2,0.) * dai_aisural_2 * (2./3. + out_coeffs.inter2g2) + dai_aisural_2 * out_therms.etaph2 ) ;
            }

          // (T,P) dependance only for thermal effects-----------------------------------------------------------------------------------------------------------------------------------

          if (milc.has_saturation(n_l, n_g1))// If eos exit we can compute thermal effects
            {
              // dT : Temperature variations, dp : Pressure variations

              const double dT_bilaninterg1 = (bilaninterg1> 0.) ? out_therms.dT_G1(n_g1)/rho(e, n_g1)- out_therms.G1/rho(e, n_g1)/rho(e, n_g1) * pch_rho->derivees().at("temperature")(e, n_g1) - out_therms.dT_etaph1(n_g1)  - alpha(e, n_g1) * ( rho_p(e, n_g1)/rho(e, n_g1)/rho(e, n_g1))* pch_rho->derivees().at("temperature")(e, n_g1)/pas_tps : 0.;

              const double dp_bilaninterg1 = (bilaninterg1> 0.) ? out_therms.dp_G1/rho(e, n_g1)- out_therms.G1/rho(e, n_g1)/rho(e, n_g1) * pch_rho->derivees().at("pression")(e, n_g1) - out_therms.dp_etaph1  - alpha(e, n_g1) * ( rho_p(e, n_g1)/rho(e, n_g1)/rho(e, n_g1))* pch_rho->derivees().at("pression")(e, n_g1)/pas_tps : 0.;
              const double dT_bilaninterg2 = (bilaninterg2> 0.) ? out_therms.dT_G2(n_g2)/rho(e, n_g1) - out_therms.G2/rho(e, n_g2)/rho(e, n_g2) * pch_rho->derivees().at("temperature")(e, n_g2) - out_therms.dT_etaph2(n_g2) - alpha(e, n_g2) * ( rho_p(e, n_g2)/rho(e, n_g2)/rho(e, n_g2))*pch_rho->derivees().at("temperature")(e, n_g2)/pas_tps : 0. ;

              const double dp_bilaninterg2 = (bilaninterg2> 0.) ?  out_therms.dp_G2/rho(e, n_g1) - out_therms.G2/rho(e, n_g2)/rho(e, n_g2) * pch_rho->derivees().at("temperature")(e, n_g2) - out_therms.dp_etaph2 - alpha(e, n_g2) * ( rho_p(e, n_g2)/rho(e, n_g2)/rho(e, n_g2))*pch_rho->derivees().at("pression")(e, n_g2)/pas_tps : 0.;


              if (Mt)//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
                {

                  (*Mt)(N * e + n_g1 , N * e + n_g1) -=  vol * ( dT_bilaninterg1 * aisural_1 * (2./3. - out_coeffs.inter2g1) -  dT_bilaninterg2 * aisural_2 * out_coeffs.inter2g2 + aisural_1 * out_therms.dT_etaph1(n_g1) ) ;

                  (*Mt)(N * e + n_g2 , N * e + n_g2) -=  vol * ( dT_bilaninterg2 * aisural_2 * (2./3. + out_coeffs.inter2g2) +  dT_bilaninterg1 * aisural_1 * out_coeffs.inter2g1 + aisural_2 * out_therms.dT_etaph1(n_g2) ) ;

                }

              if (Mp)//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
                {

                  (*Mp)(N * e + n_g1 , e ) -= vol * ( dp_bilaninterg1 * aisural_1 * (2./3. - out_coeffs.inter2g1) - dp_bilaninterg2 * aisural_2 * out_coeffs.inter2g2 + aisural_1 * out_therms.dp_etaph1 ) ;

                  (*Mp)(N * e + n_g2 ,  e ) -= vol * ( dp_bilaninterg2 * aisural_2 * (2./3. + out_coeffs.inter2g2) +  dp_bilaninterg1 * aisural_1 * out_coeffs.inter2g1 + aisural_2 * out_therms.dp_etaph2 ) ;

                }
            }

        }

    }
}
