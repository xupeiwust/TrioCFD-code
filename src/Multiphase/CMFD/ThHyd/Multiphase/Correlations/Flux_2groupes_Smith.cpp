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
/////////////////////////////////////////////////////////////////////////////

/// T.R. Smith, J.P. Schlegel, T. Hibiki, M. Ishii, Mechanistic modeling of interfacial area transport in large diameter pipes,
/// International Journal of Multiphase Flow, Volume 47, 2012, https://doi.org/10.1016/j.ijmultiphaseflow.2012.06.009
/// and Joseph L. Bottini, Taiyang Zhang, Caleb S. Brooks, Validation of two-group interfacial area transport equation in boiling flow,
/// International Journal of Heat and Mass Transfer, 2024, https://doi.org/10.1016/j.ijheatmasstransfer.2024.125515

/////////////////////////////////////////////////////////////////////////////

#ifndef __APPLE__
#include <boost/math/special_functions/gamma.hpp>
#endif

#include <Flux_2groupes_Smith.h>
#include <Pb_Multiphase.h>
Implemente_instanciable(Flux_2groupes_Smith, "Flux_2groupes_Smith", Flux_2groupes_base);

Sortie& Flux_2groupes_Smith::printOn(Sortie& os) const
{
  return os;
}

Entree& Flux_2groupes_Smith::readOn(Entree& is)
{
#ifdef __APPLE__
  Process::exit("Flux_2groupes_Smith could not be used with Mac. Contact the support !");
#endif

  Param param(que_suis_je());
  param.ajouter("hPNVG", &hPNVG_);
  param.ajouter("Xi_h", &Xi_h_);
  param.ajouter("A_c", &A_c_);
  param.ajouter("R", &R_);
  param.ajouter("theta", &theta_rad);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Flux_2groupes_Smith::coeffs(const input_coeffs& in, output_coeffs& out) const
{
#ifndef __APPLE__
  // Initialisation--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  double eta = 0 ;
  double deta_dalpha1 = 0 ;
  double deta_dai1 = 0 ;
  double deta_dalpha2 = 0 ;
  double deta_dai2 = 0 ;

  // Numerical security----------------------------------------------------------------------------------------------------------------------------------------------------------------------

  const double fac_sec = 1.e4;
  int clipping = 0;

  // Number often used in source terms ------------------------------------------------------------------------------------------------------------------------------------------------------

  const double alpha_coal_max = 0.509 ;
  const double alpha_max_cbrt = std::cbrt(alpha_max) ;
  const double alpha_clipped = std::min(in.alpha(in.n_g1), alpha_coal_max) ;
  const double alpha_clipped_cbrt = std::cbrt(alpha_clipped) ;
  const double alphag1_5_over3 = std::cbrt(in.alpha(in.n_g1)) * std::cbrt(in.alpha(in.n_g1)) * std::cbrt(in.alpha(in.n_g1)) * std::cbrt(in.alpha(in.n_g1)) * std::cbrt(in.alpha(in.n_g1))  ;
  const double alphag2_1_over3 = std::cbrt(in.alpha(in.n_g2))  ;
  const double alphag2_4_over3 = alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3 * alphag2_1_over3  ;
  const double ai1_cbrt = std::cbrt(in.a_i(in.n_g1)) ;
  const double ai2_cbrt = std::cbrt(in.a_i(in.n_g2)) ;
  const double eps_1_over3 = std::cbrt(in.epsilon(in.n_l)) ;
  const double alphafonction_3_over2 = std::sqrt(1.-in.alpha(in.n_g1)) * std::sqrt(1.-in.alpha(in.n_g1)) * std::sqrt(1.-in.alpha(in.n_g1)) ;


  // Physical constants------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  const double D_crit = 4. * std::sqrt(in.sigma(in.n_l, in.n_g1) / g / (in.rho(in.n_l)-in.rho(in.n_g1)));
  const double D_etoile = D_crit/in.d_bulles(in.n_g1);
  const double D_etoile2 = D_crit/in.d_bulles(in.n_g2);
  const double khi_d =  b/3. *((3./2. * D_etoile) * std::pow(b*3./2. * D_etoile , p) * std::exp(-b * 3./2.* D_etoile)) / (boost::math::tgamma(p+1.,0.)-boost::math::tgamma(p+1.,b * 3./2. * D_etoile));
  const double C_D1 = std::max(2./3. * in.d_bulles(in.n_g1) * std::sqrt( g * std::abs(in.rho(in.n_l) - in.rho(in.n_g1) ) / in.sigma(in.n_g1, in.n_l) ) * ( (1.+ 17.67 * std::pow(1.-in.alpha(in.n_g1), 9./7.) )/ (18.67 * alphafonction_3_over2 ) ) * ( (1.+ 17.67 * std::pow(1.-in.alpha(in.n_g1), 9./7.) )/ (18.67 * alphafonction_3_over2 ) ) , 1./fac_sec);
  const double Ur1 = std::min(in.nv(in.n_l,in.n_g1),std::sqrt(4./3. * in.d_bulles(in.n_g1) / C_D1 * g *std::abs(in.rho(in.n_l) - in.rho(in.n_g1) ) / in.rho(in.n_l) * (in.alpha(in.n_l))));
  const double Ur2 = std::max(in.nv(in.n_l,in.n_g2),std::sqrt(1./2. * D_crit *  g *std::abs(in.rho(in.n_l) - in.rho(in.n_g2) ) / in.rho(in.n_l) ));
  const double We2 = 2. * in.rho(in.n_l) * std::cbrt(in.epsilon(in.n_l)*std::max(in.d_bulles(in.n_g2),D_crit)) * std::cbrt(in.epsilon(in.n_l)*std::max(in.d_bulles(in.n_g2),D_crit)) * std::max(in.d_bulles(in.n_g2),D_crit) / in.sigma(in.n_g2, in.n_l) ;

  // Fonctions in source terms--------------------------------------------------------------------------------------------------------------------------------------------------------------------

  const double lambda_RC1 = ( in.alpha(in.n_g1)>1./fac_sec ) ? std::exp( - C_RC0 * std::sqrt(in.d_bulles(in.n_g1) * in.rho(in.n_l) / in.sigma(in.n_g1, in.n_l)) * std::cbrt(in.d_bulles(in.n_g1) * in.epsilon(in.n_l) ) ) : 0. ;
  const double formfunc_alpha = ( in.alpha(in.n_g1)>1./fac_sec ) ?  (1. - std::exp(-C_RC1 * (alpha_max_cbrt * alpha_clipped_cbrt ) / ( alpha_max_cbrt - alpha_clipped_cbrt  ) ) ) : 0. ;
  const double dformfunc_alpha_dalpha = ((in.alpha(in.n_g1)>1./fac_sec) && (in.alpha(in.n_g1) < alpha_coal_max))? ( 1./ alpha_max_cbrt / 3. * (1.-std::exp(-C_RC0 * (alpha_max_cbrt * alpha_clipped_cbrt ) / ( alpha_max_cbrt - alpha_clipped_cbrt  ) )) / (alpha_clipped_cbrt * alpha_clipped_cbrt) / ( alpha_max_cbrt - alpha_clipped_cbrt) / ( alpha_max_cbrt - alpha_clipped_cbrt) + std::exp(-C_RC0 * (alpha_max_cbrt * alpha_clipped_cbrt ) / ( alpha_max_cbrt - alpha_clipped_cbrt  ) ) / (alpha_max_cbrt - alpha_clipped_cbrt ) * (C_RC0 / alpha_clipped_cbrt / 3. ) / (alpha_max_cbrt - alpha_clipped_cbrt) * (  1. / std::min(in.alpha(in.n_g1)*in.alpha(in.n_g1), alpha_coal_max * alpha_coal_max) + 1. / ( alpha_max_cbrt - alpha_clipped_cbrt  ) )) : 0.;
  const double Dc_etoile_16_over3 = std::cbrt(D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) * std::cbrt(D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) * std::cbrt(D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) * std::cbrt(D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) * std::cbrt(D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) * std::cbrt(D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) * std::cbrt(D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) * std::cbrt(D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) * std::cbrt(D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) * std::cbrt(D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) * std::cbrt(D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) * std::cbrt(D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) * std::cbrt(D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) * std::cbrt(D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) * std::cbrt(D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) * std::cbrt(D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) ;
  const double Dc_etoile_6 = (D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) * (D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) * (D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) * (D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) * (D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) * (D_crit/std::max(in.d_bulles(in.n_g2),D_crit)) ;

  //Eta RC (112) coefficient for secmem---------------------------------------------------------------------------------------------------------------------------------------------------------------

  eta +=  (in.alpha(in.n_g1)>1./fac_sec) ?  3.15 * C_RC1 * lambda_RC1 * formfunc_alpha /(alpha_max_cbrt * alpha_max_cbrt ) * (1.-2./3. * D_etoile) * eps_1_over3 * alpha_clipped * alpha_clipped * (ai1_cbrt * ai1_cbrt)   : 0;

  //dEta RC (112) coefficient for mat ----------------------------------------------------------------------------------------------------------------------------------------------------------------

  deta_dalpha1 += (in.alpha(in.n_g1)>1./fac_sec) ?  3.15 * C_RC1 * lambda_RC1 * formfunc_alpha /(alpha_max_cbrt * alpha_max_cbrt) * (1.-2./3. * D_etoile) * eps_1_over3 * 2.* alpha_clipped * (ai1_cbrt * ai1_cbrt) + 3.15 * C_RC1 * lambda_RC1 * dformfunc_alpha_dalpha /(alpha_max_cbrt * alpha_max_cbrt ) * (1.-2./3. * D_etoile) * eps_1_over3 * alpha_clipped * alpha_clipped * (ai1_cbrt * ai1_cbrt) : 0; // dEta/dalpha1

  deta_dai1 += (in.alpha(in.n_g1)>1./fac_sec) ?  3.15 * C_RC1 * lambda_RC1 * formfunc_alpha /(alpha_max_cbrt * alpha_max_cbrt ) * (1.-2./3. * D_etoile) * eps_1_over3 * alpha_clipped * alpha_clipped * 2./3. / std::min(ai1_cbrt,fac_sec) : 0; // dEta/dai1

  //Eta RC (122) coefficient for secmem--------------------------------------------------------------------------------------------------------------------------------------------------------------

  eta +=  (in.alpha(in.n_g1)>1./fac_sec) ?  1.44 * C_RC122 * lambda_RC1 * formfunc_alpha * eps_1_over3 * alphag1_5_over3 * alphag2_4_over3 * (ai2_cbrt * ai2_cbrt)   : 0;

  //dEta RC (122) coefficient for mat -----------------------------------------------------------------------------------------------------------------------------------------------------------------

  deta_dalpha1 += (in.alpha(in.n_g2)>1./fac_sec) ?  1.44 * C_RC122 * lambda_RC1 * formfunc_alpha * eps_1_over3 * 5./3. * (alpha_clipped_cbrt * alpha_clipped_cbrt) * (alpha_clipped * alpha_clipped_cbrt * alpha_clipped_cbrt * alpha_clipped_cbrt) * (ai2_cbrt * ai2_cbrt) + 1.44 * C_RC122 * lambda_RC1 * dformfunc_alpha_dalpha * eps_1_over3 * alphag1_5_over3 * (alpha_clipped_cbrt * alpha_clipped_cbrt * alpha_clipped_cbrt * alpha_clipped_cbrt) * (ai2_cbrt * ai2_cbrt) : 0; // dEta/dalpha1

  deta_dalpha2 +=  (in.alpha(in.n_g2)>1./fac_sec) ?  1.44 * C_RC122 * lambda_RC1 * formfunc_alpha * eps_1_over3 * alphag1_5_over3 * 4./3.* std::cbrt(in.alpha(in.n_g2)) * (ai2_cbrt * ai2_cbrt)   : 0;  // dEta/dalpha2

  deta_dai2 += (in.alpha(in.n_g2)>1./fac_sec) ?  1.44 * C_RC122 * lambda_RC1 * formfunc_alpha * eps_1_over3 * alphag1_5_over3 * alphag2_4_over3 * 2./3. / std::min(ai2_cbrt,fac_sec) : 0; // dEta/dai2


  //Eta WE coefficient for secmem --------------------------------------------------------------------------------------------------------------------------------------------------------------------

  eta += (in.alpha(in.n_g1)>1./fac_sec) ?  3.85 * C_WE1 * std::cbrt(C_D1) * Ur1 * in.alpha(in.n_g1) * in.a_i(in.n_g1) : 0;

  //dEta WE coefficient for mat--------------------------------------------------------------------------------------------------------------------------------------------------------

  deta_dalpha1 += (in.alpha(in.n_g1)>1./fac_sec) ?  3.85 * C_WE1 * std::cbrt(C_D1) * Ur1  * in.a_i(in.n_g1) : 0;

  deta_dai1 += (in.alpha(in.n_g1)>1./fac_sec) ?  3.85 * C_WE1 * std::cbrt(C_D1) * Ur1 * in.alpha(in.n_g1) : 0;


  //Eta TI coefficient for secmem-----------------------------------------------------------------------------------------------------------------------------------------------------------

  eta += (in.alpha(in.n_g2)>1./fac_sec) ?  - 11.65 * C_TI21 * std::exp(-We_cr2/We2) * std::sqrt(1.-std::min(We_cr1/We2 ,1.))*std::max(0.15 * Dc_etoile_16_over3 - 0.117 * Dc_etoile_6,0.) * eps_1_over3 * in.alpha(in.n_l) * alphag2_1_over3 * (ai2_cbrt * ai2_cbrt) : 0;

  //dEta TI coefficient for mat--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  deta_dalpha2 += (in.alpha(in.n_g2)>1./fac_sec) ?  - 11.65 * C_TI21 * std::exp(-We_cr2/We2) * std::sqrt(1.-std::min(We_cr1/We2 ,1.))*std::max(0.15 * Dc_etoile_16_over3 - 0.117 * Dc_etoile_6,0.) * eps_1_over3 * in.alpha(in.n_l) * 1./3. / (std::min(alphag2_1_over3 * alphag2_1_over3,fac_sec)) * (ai2_cbrt * ai2_cbrt) : 0; // dEta/dalpha2

  deta_dai2 += (in.alpha(in.n_g2)>1./fac_sec) ?  - 11.65 * C_TI21 * std::exp(-We_cr2/We2) * std::sqrt(1.-std::min(We_cr1/We2 ,1.))*std::max(0.15 * Dc_etoile_16_over3 - 0.117 * Dc_etoile_6,0.) * eps_1_over3 * in.alpha(in.n_l) * alphag2_1_over3 * 2./3. / (std::min(ai2_cbrt,fac_sec)) : 0; // dEta/dai2


  //Eta SO coefficient for secmem----------------------------------------------------------------------------------------------------------------------------------------------------------------

  eta += (in.alpha(in.n_g2)>1./fac_sec) ?  -2.33 * C_S0 * in.sigma(in.n_g2, in.n_l) / in.rho(in.n_g2) / Ur2 * std::max(1. - We_SO / We2 ,0.) * (in.a_i(in.n_g2)) * (in.a_i(in.n_g2)) / in.alpha(in.n_g2) : 0. ;

  //dEta SO coefficient for mat---------------------------------------------------------------------------------------------------------------------------------------------------------------

  deta_dalpha2 += (in.alpha(in.n_g2)>1./fac_sec) ?  2.33 * C_S0 * in.sigma(in.n_g2, in.n_l) / in.rho(in.n_g2) / Ur2 * std::max(1. - We_SO / We2 ,0.) * (in.a_i(in.n_g2) ) * (in.a_i(in.n_g2)) / (in.alpha(in.n_g2) * in.alpha(in.n_g2)) : 0; // dEta/dalpha2

  deta_dai2 += (in.alpha(in.n_g2)>1./fac_sec) ?  -2.33 * C_S0 * in.sigma(in.n_g2, in.n_l) / in.rho(in.n_g2) / Ur2 * std::max(1. - We_SO / We2 ,0.) * 2. * (in.a_i(in.n_g2)) / in.alpha(in.n_g2) : 0; //dEta/dai2


  // Physical/numerical security for mass exchange between phases --------------------------------------------------------------------------------------------------------------------------
  if (eta > 0.) // if transfer from g1 to g2
    {
      if (1./fac_sec > in.alpha(in.n_g1)) // if not enough g1
        {

          clipping = 1;
        }
    }
  if (0. > eta ) // if transfer from g2 to g1
    {
      if (1./fac_sec >in.alpha(in.n_g2)) // if not enough g2
        {

          clipping = 1;
        }
    }

  // Put in  out-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  out.gamma(in.n_g1, in.n_g1) = (clipping == 0) ? eta : 0. ; // G1
  out.gamma(in.n_g2, in.n_g2) =  (clipping == 0) ? -eta : 0. ; // G2 = -G1
  out.da_gamma(in.n_g1, in.n_g1) =  (clipping == 0) ? deta_dalpha1 : 0. ; // dG1/dalpha1
  out.da_gamma(in.n_g2, in.n_g2) = (clipping == 0) ? deta_dalpha2 : 0. ; // dG1/dalpha2
  out.dai_gamma(in.n_g1, in.n_g1) = (clipping == 0) ? deta_dai1 : 0. ; // dG1/dai1
  out.dai_gamma(in.n_g2, in.n_g2) = (clipping == 0) ? deta_dai2 : 0. ; // dG1/dai2


  out.inter2g1 = khi_d * D_etoile * D_etoile ;
  out.da_inter2g1 = 0.;

  out.inter2g2=  khi_d * D_etoile2 * D_etoile2 * ((in.alpha(in.n_g2)>1./fac_sec) ? in.alpha(in.n_g1)/in.alpha(in.n_g2) * (in.d_bulles(in.n_g2)/in.d_bulles(in.n_g1))*(in.d_bulles(in.n_g2)/in.d_bulles(in.n_g1))*(in.d_bulles(in.n_g2)/in.d_bulles(in.n_g1)) : 0.) ;
  out.da_inter2g2 =  khi_d * D_etoile2 * D_etoile2 * ((in.alpha(in.n_g2)>1./fac_sec) ? -in.alpha(in.n_g1)/in.alpha(in.n_g2)/in.alpha(in.n_g2) * (in.d_bulles(in.n_g2)/in.d_bulles(in.n_g1))*(in.d_bulles(in.n_g2)/in.d_bulles(in.n_g1))*(in.d_bulles(in.n_g2)/in.d_bulles(in.n_g1)) : 0.)  ;

  out.inter3g1 =  khi_d * D_etoile * D_etoile * D_etoile  ;
  out.da_inter3g1 = 0.;

  out.inter3g2 =  khi_d * D_etoile2 * D_etoile2 * D_etoile2 * ((in.alpha(in.n_g2)>1./fac_sec) ? in.alpha(in.n_g1)/in.alpha(in.n_g2) * (in.d_bulles(in.n_g2)/in.d_bulles(in.n_g1))*(in.d_bulles(in.n_g2)/in.d_bulles(in.n_g1))*(in.d_bulles(in.n_g2)/in.d_bulles(in.n_g1)) : 0.)  ;
  out.da_inter3g2 =  khi_d * D_etoile2 * D_etoile2 * D_etoile2 * ((in.alpha(in.n_g2)>1./fac_sec) ? -in.alpha(in.n_g1)/in.alpha(in.n_g2)/in.alpha(in.n_g2) * (in.d_bulles(in.n_g2)/in.d_bulles(in.n_g1))*(in.d_bulles(in.n_g2)/in.d_bulles(in.n_g1))*(in.d_bulles(in.n_g2)/in.d_bulles(in.n_g1)) : 0.) ;
#endif
}

void Flux_2groupes_Smith::therm(const input_therms& in, output_therms& out) const
{

  // Numerical security----------------------------------------------------------------------------------------------------------------------------------------------------------------------

  const double fac_sec = 1.e4;

  // For lisibility----------------------------------------------------------------------------------------------------------------------------------------------------------------------

  const double dsm1 = in.d_bulles(in.n_g1);
  const double dsm2 = in.d_bulles(in.n_g2);


  // Physical constants ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

  const double Prl = in.mu(in.n_l) * in.Cp(in.n_l) / in.lambda(in.n_l) ;
  const double Ja1 = std::max(in.rho(in.n_l) * in.Cp(in.n_l) * (in.Tsatg1-in.T(in.n_l))/in.rho(in.n_g1)/in.Lvap,0.);
  const double Ja2 = std::max(in.rho(in.n_l) * in.Cp(in.n_l)*(in.Tsatg2-in.T(in.n_l))/in.rho(in.n_g2)/in.Lvap,0.);
  const double dTl_Ja1 = ( 0.< Ja1) ? - in.rho(in.n_l) * in.Cp(in.n_l) / in.rho(in.n_g1) / in.Lvap : 0.  ;
  const double dp_Ja1 = ( 0.< Ja1) ? in.rho(in.n_l) * in.Cp(in.n_l) * (in.dp_Tsat1)/in.rho(in.n_g1)/in.Lvap : 0. ;
  const double dTl_Ja2 = ( 0.< Ja2) ? - in.rho(in.n_l) * in.Cp(in.n_l) / in.rho(in.n_g2) / in.Lvap : 0. ;
  const double dp_Ja2 = ( 0.< Ja2) ? in.rho(in.n_l) * in.Cp(in.n_l) * (in.dp_Tsat2)/in.rho(in.n_g2)/in.Lvap : 0. ;


// Initialisation----------------------------------------------------------------------------------------------------------------------------------------------------------------------

  double etaCo = 0.;
  double dak_etaCo = 0.;
  double dal_etaCo = 0.;
  double dTk_etaCo = 0. ;
  double dTl_etaCo =  0.;
  double dp_etaCo = 0.;

  double etaPc =  0.;
  double dak_etaPc = 0.;
  double dal_etaPc = 0.;
  double dTk_etaPc = 0. ;
  double dTl_etaPc = 0.;
  double dp_etaPc =  0.;

  double GWn = 0.;
  double dak_GWn = 0.;
  double dal_GWn = 0.;
  double dTk_GWn = 0. ;
  double dTl_GWn =  0.;
  double dp_GWn = 0.;

  double Reb = 0.; // used for both phases cant be const
  double Nuc = 0.; // used for both phases cant be const


  // Interfacial area monitored by Dsm for numerical safety -----------------------------------------------------------------------------------------------------------------------------------------------

  const double Xai = (in.alpha(in.n_g2)>1./fac_sec) ? (6. * in.alpha(in.n_g1) / dsm1) / ((6. * in.alpha(in.n_g1) / dsm1) + (6. * in.alpha(in.n_g2) / dsm2)) : 1. ;
  const double da1_Xai = (in.alpha(in.n_g2)>1./fac_sec) ? 36. / dsm1 / dsm2 * in.alpha(in.n_g2) / ((6. * in.alpha(in.n_g1) / dsm1) + (6. * in.alpha(in.n_g2) / dsm2)* (6. * in.alpha(in.n_g1) / dsm1) + (6. * in.alpha(in.n_g2) / dsm2)) : 0.;
  const double da2_1_Xai = (in.alpha(in.n_g2)>1./fac_sec) ? 36. / dsm2 / dsm1 * in.alpha(in.n_g1) / ((6. * in.alpha(in.n_g2) / dsm1) + (6. * in.alpha(in.n_g1) / dsm1)* (6. * in.alpha(in.n_g2) / dsm2) + (6. * in.alpha(in.n_g1) / dsm1)) : 0.;

// Group 1 transfer -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  Reb = std::sqrt(2.)*std::pow(( in.rho(in.n_l)- in.rho(in.n_g1) )*g * in.sigma(in.n_g1, in.n_l)/in.rho(in.n_l) / in.rho(in.n_l), 1./4.) * in.rho(in.n_l) * betac * dsm1 / in.mu(in.n_l);
  Nuc = 4./M_PI*std::sqrt(Reb)*std::sqrt(Prl);

  // Condensation -------------------------------------------------------------------------------------------------------------------------------------------------------------------

  const double fac_Co1 = betac*betac*betac/(1.-betac*betac)/9. ;
  etaCo = -fac_Co1* (6. * in.alpha(in.n_g1) / dsm1 ) *(1.-in.alpha(in.n_l)) / dsm1 * Nuc * Ja1 ;
  dak_etaCo = -fac_Co1* (6. / dsm1 ) *(1.-in.alpha(in.n_l)) / dsm1 * Nuc * Ja1 ;
  dal_etaCo = fac_Co1* (6. * in.alpha(in.n_g1) / dsm1 ) / dsm1 * Nuc * Ja1 ;
  dTk_etaCo = 0. ;
  dTl_etaCo =  -fac_Co1* (6. * in.alpha(in.n_g1) / dsm1 ) *(1.-in.alpha(in.n_l)) / dsm1 * Nuc * dTl_Ja1 ;
  dp_etaCo =  -fac_Co1* (6. * in.alpha(in.n_g1) / dsm1 ) *(1.-in.alpha(in.n_l)) / dsm1 * Nuc * dp_Ja1 ;

  // Phase change ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

  const double fac_Pc1 = (1.-pc);
  etaPc =  -fac_Pc1 * (6. * in.alpha(in.n_g1) / dsm1 ) /   dsm1 * (1.-in.alpha(in.n_l)) * Nuc * Ja1 ;
  dak_etaPc = -fac_Pc1 * (6. / dsm1 ) /   dsm1 * (1.-in.alpha(in.n_l)) * Nuc * Ja1 ;
  dal_etaPc =-fac_Pc1 * (6. * in.alpha(in.n_g1) / dsm1 ) /   dsm1 * Nuc * Ja1 ;
  dTk_etaPc = 0. ;
  dTl_etaPc = -fac_Pc1 * (6. * in.alpha(in.n_g1) / dsm1 ) /   dsm1 * (1.-in.alpha(in.n_l)) * Nuc * dTl_Ja1 ;
  dp_etaPc =  -fac_Pc1 * (6. * in.alpha(in.n_g1) / dsm1 ) /   dsm1 * (1.-in.alpha(in.n_l)) * Nuc * dp_Ja1 ;

  // Wall nucleation -------------------------------------------------------------------------------------------------------------------------------------------------------------------

  GWn = ( 0.< in.hl - hPNVG_) ?  Xai * Xi_h_* in.qp / A_c_ * (in.hl - hPNVG_)/( in.hlsat - hPNVG_ ) /  in.Lvap : 0. ;
  dak_GWn = ( 0.< in.hl - hPNVG_) ?  da1_Xai * Xi_h_* in.qp / A_c_ * (in.hl - hPNVG_)/( in.hlsat - hPNVG_ ) /  in.Lvap : 0. ;

  // Balance for G1 ------------------------------------------------------------------------------------------------------------------------------------------------------------------

  out.G1 = in.rho(in.n_g1) * (etaCo + etaPc ) + GWn;
  out.etaph1 = etaCo + GWn / in.rho(in.n_g1) ;

  out.dT_G1(in.n_g1) = in.rho(in.n_g1) * (dTk_etaPc + dTk_etaCo ) + dTk_GWn;
  out.dT_G1(in.n_l) = in.rho(in.n_g1) * (dTl_etaPc + dTl_etaCo ) + dTl_GWn ;
  out.dT_etaph1(in.n_g1) = dTk_etaCo + dTk_GWn / in.rho(in.n_g1) ;
  out.dT_etaph1(in.n_l) = dTl_etaCo + dTl_GWn / in.rho(in.n_g1) ;


  out.da_G1(in.n_g1) = in.rho(in.n_g1) * (dak_etaPc + dak_etaCo ) + dak_GWn ;
  out.da_G1(in.n_l) = in.rho(in.n_g1) * (dal_etaPc + dal_etaCo ) + dal_GWn ;
  out.da_etaph1(in.n_g1) = dak_etaCo + dak_GWn / in.rho(in.n_g1) ;
  out.da_etaph1(in.n_l) = dal_etaCo + dal_GWn / in.rho(in.n_g1) ;

  out.dp_G1 = in.rho(in.n_g1) * (dp_etaPc + dp_etaCo ) + dp_GWn ;
  out.dp_etaph1 = dp_etaCo + dp_GWn / in.rho(in.n_g1) ;

// Group 2 transfer ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  Reb = 0.35*std::pow((in.rho(in.n_l)-in.rho(in.n_g2))* g * in.dh / in.rho(in.n_l) , 1./4.) * in.rho(in.n_l) * dsm2 / in.mu(in.n_l);
  Nuc = 0.185 * std::pow(Reb , 0.7) * std::sqrt(Prl) ;

  // No condensation
  // Phase change ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

  etaPc = - (6. * in.alpha(in.n_g2) / dsm2 ) /dsm2 * Nuc * Ja2 * (1.-in.alpha(in.n_l)) ;
  dak_etaPc = - (6. / dsm2 ) / dsm2 * Nuc * Ja2 * (1.-in.alpha(in.n_l)) ;
  dal_etaPc = (6. * in.alpha(in.n_g2) / dsm2 ) / dsm2 * Nuc * Ja2 ;
  dTk_etaPc = 0. ;
  dTl_etaPc = - (6. * in.alpha(in.n_g2) / dsm2 ) / dsm2 * Nuc  * (1.-in.alpha(in.n_l)) * dTl_Ja2 ;
  dp_etaPc = - (6. * in.alpha(in.n_g2) / dsm2 ) /dsm2 * Nuc  * (1.-in.alpha(in.n_l)) * dp_Ja2 ;

  // Wall nucleation -------------------------------------------------------------------------------------------------------------------------------------------------------------------

  GWn = ( 0.< in.hl - hPNVG_) ?  (1.-Xai) * Xi_h_* in.qp / A_c_ * (in.hl - hPNVG_)/( in.hlsat - hPNVG_ ) /  in.Lvap : 0.;
  dak_GWn = ( 0.< in.hl - hPNVG_) ?  da2_1_Xai * Xi_h_* in.qp / A_c_ * (in.hl - hPNVG_)/( in.hlsat - hPNVG_ ) /  in.Lvap : 0. ;

  // Balance for G2 ------------------------------------------------------------------------------------------------------------------------------------------------------------------

  out.G2 = in.rho(in.n_g2) * ( etaPc ) + GWn;
  out.etaph2 = GWn / in.rho(in.n_g2) ;
  out.dT_G2(in.n_g2) = in.rho(in.n_g2) * (dTk_etaPc  ) + dTk_GWn;
  out.dT_G2(in.n_l) = in.rho(in.n_g2) * (dTl_etaPc  ) + dTl_GWn ;
  out.dT_etaph2(in.n_g2) =  dTk_GWn / in.rho(in.n_g2) ;
  out.dT_etaph2(in.n_l) =  dTl_GWn / in.rho(in.n_g2) ;

  out.da_G2(in.n_g2) = in.rho(in.n_g2) * (dak_etaPc  ) + dak_GWn ;
  out.da_G2(in.n_l) = in.rho(in.n_g2) * (dal_etaPc  ) + dal_GWn ;
  out.da_etaph2(in.n_g2) =   dak_GWn / in.rho(in.n_g2) ;
  out.da_etaph2(in.n_l) =   dal_GWn / in.rho(in.n_g2) ;

  out.dp_G2 = in.rho(in.n_g2) * (dp_etaPc ) + dp_GWn ;
  out.dp_etaph2 =  dp_GWn / in.rho(in.n_g2) ;


}
