/****************************************************************************
* Copyright (c) 2022, CEA
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

/////////////////////////////////////////////////////////////////////////////


#include <Coalescence_bulles_2groupes_Smith.h>
#include <Pb_Multiphase.h>


Implemente_instanciable(Coalescence_bulles_2groupes_Smith, "Coalescence_bulles_2groupes_Smith", Coalescence_bulles_2groupes_base);


Sortie& Coalescence_bulles_2groupes_Smith::printOn(Sortie& os) const
{
  return os;
}

Entree& Coalescence_bulles_2groupes_Smith::readOn(Entree& is)
{

  return is;
}

void Coalescence_bulles_2groupes_Smith::coefficient_RC(const DoubleTab& alpha, const DoubleTab& alpha_p, const DoubleTab& p, const DoubleTab& T,
                                                       const DoubleTab& rho, const DoubleTab& nu, const DoubleTab& sigma, const double Dh,
                                                       const DoubleTab& ndv, const DoubleTab& d_bulles,
                                                       const DoubleTab& eps, const DoubleTab& k_turb, const int n_l, const int n_g1, const int n_g2,
                                                       DoubleTab& coeff) const
{

  const double fac_sec = 1.e4 ; // numerical security
  const double lambda_RC1 = std::exp( - C_RC0 * std::sqrt(d_bulles(n_g1) * rho(n_l) / sigma(n_g1, n_l)) * std::cbrt(d_bulles(n_g1) * eps(n_l)) ) ;
  const double D_crit = 4. * std::sqrt(sigma(n_g2, n_l) / g / (rho(n_l)-rho(n_g2)));
  const double lambda_RC2 = std::exp( - C_RC0 * std::sqrt(std::max(d_bulles(n_g2),D_crit) * rho(n_l) / sigma(n_g2, n_l)) * std::cbrt(std::max(d_bulles(n_g2),D_crit) * eps(n_l) ) ) ;
  const double alpha_coal_max = 0.509 ;
  const double alpha_max_cbrt = std::cbrt(alpha_max) ;
  const double alpha_clipped_cbrt = std::cbrt(std::min(alpha(n_g1), alpha_coal_max)) ;
  const double alpha_p_clipped_cbrt = std::cbrt(std::min(alpha_p(n_g1), alpha_coal_max)) ;

  const double formfunc_alpha = ( alpha(n_g1)>1./fac_sec ) ? (1. - std::exp(-C_RC0 * (alpha_max_cbrt * alpha_clipped_cbrt ) / std::max( alpha_max_cbrt - alpha_clipped_cbrt ,1./fac_sec ) ) ) : 0. ;
  const double formfunc_alpha_p = ( alpha(n_g1)>1./fac_sec ) ? (1. - std::exp(-C_RC0 * (alpha_max_cbrt * alpha_p_clipped_cbrt ) / std::max( alpha_max_cbrt - alpha_p_clipped_cbrt  ,1./fac_sec) ) ) : 0. ;
  const double dformfunc_alpha_dalpha = ( alpha(n_g1)>1./fac_sec ) ? ( 1./ alpha_max_cbrt / 3. * (1.-std::exp(-C_RC0 * (alpha_max_cbrt * alpha_clipped_cbrt ) / ( alpha_max_cbrt - alpha_clipped_cbrt  ) )) / (alpha_clipped_cbrt * alpha_clipped_cbrt) / ( alpha_max_cbrt - alpha_clipped_cbrt) / ( alpha_max_cbrt - alpha_clipped_cbrt) + std::exp(-C_RC0 * (alpha_max_cbrt * alpha_clipped_cbrt ) / ( alpha_max_cbrt - alpha_clipped_cbrt  ) ) / (alpha_max_cbrt - alpha_clipped_cbrt ) * (C_RC0 / alpha_clipped_cbrt / 3. ) / (alpha_max_cbrt - alpha_clipped_cbrt) * (  1. / std::min(alpha(n_g1)*alpha(n_g1), alpha_coal_max) + 1. / ( alpha_max_cbrt - alpha_clipped_cbrt  ) )) : 0.;
  const double formfunc_alpha2 = std::max(1. - std::exp(-C_RC0 * (alpha(n_g2) * alpha(n_g2) ) ) ,0.) * std::max(1.- 0.37 * ( D_crit/std::max(d_bulles(n_g2),D_crit) * D_crit/std::max(d_bulles(n_g2),D_crit) * D_crit/std::max(d_bulles(n_g2),D_crit) ),0.) ;

  // RC (1) coefficient for secmem---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  coeff(n_g1, n_l) = ( alpha(n_g1)>1./fac_sec ) ?  - 0.17 * C_RC1 * lambda_RC1 * formfunc_alpha /(alpha_max_cbrt / (alpha_max_cbrt - alpha_clipped_cbrt ) ) : 0.  ;

  // dRC (1) coefficient for mat---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  if (alpha(n_g1) < alpha_coal_max)
    {
      coeff(n_l, n_g1) = ( alpha(n_g1)>1./fac_sec ) ? - 0.17 * C_RC1 * lambda_RC1 * dformfunc_alpha_dalpha : 0.; // dRC1/dalpha
      coeff(n_l, n_l) =  ( alpha(n_g1)>1./fac_sec ) ?  -1.14 * C_RC122 * ( lambda_RC2 * C_RC0 * ( alpha_max_cbrt * alpha_max_cbrt )  / 3. / ( alpha_clipped_cbrt * alpha_clipped_cbrt ) / ( (alpha_max_cbrt - alpha_clipped_cbrt) * (alpha_max_cbrt - alpha_clipped_cbrt)) * std::exp(-C_RC0 * (alpha_max_cbrt * alpha_clipped_cbrt ) / ( alpha_max_cbrt - alpha_clipped_cbrt  ) ) )    : 0. ; //dRC122 in 1 /dalpha1
    }
  else
    {
      coeff(n_l, n_g1) = 0. ;
      coeff(n_l, n_l) = 0. ;

    }

  // RC (2) coefficient for secmem---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  coeff(n_g2, n_l) = ( alpha(n_g2)>1./fac_sec ) ? - 95.7 * C_RC2 * lambda_RC2 * formfunc_alpha2 / Dh / Dh  : 0. ;

  // dRC (2) coefficient for mat---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  coeff(n_l, n_g2) = ( alpha(n_g2)>1./fac_sec ) ? - 95.7 * C_RC2 * lambda_RC2 * std::max(1.- 0.37 * (D_crit/std::max(d_bulles(n_g2),D_crit)) * (D_crit/std::max(d_bulles(n_g2),D_crit)) * (D_crit/std::max(d_bulles(n_g2),D_crit)) * (D_crit/std::max(d_bulles(n_g2),D_crit)),0.) / Dh / Dh * C_RC0 * std::exp(-C_RC0 * (alpha(n_g2) *alpha(n_g2) ) ) / 2. / std::sqrt(std::min(alpha(n_g2), 0.81)) : 0.  ; // dRC2/dlapha2

  // RC (112) coefficient for 1st group for secmem---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  coeff(n_g1, n_g1) = ( alpha(n_g1)>1./fac_sec ) ? 4.1 * C_RC1 * lambda_RC1 / (alpha_max_cbrt * alpha_max_cbrt) * formfunc_alpha_p * std::max( 1. - 2./3. * D_crit / d_bulles(n_g1),0.)  : 0.;

  // RC (122) coefficient for 1st group for secmem---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  coeff(n_g1, n_g2) = ( alpha(n_g1)>1./fac_sec  ) ? -1.14 * C_RC122 * formfunc_alpha * lambda_RC2 : 0.;

  // RC (122) coefficient for 2nd group for secmem---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  coeff(n_g2, n_g1) = ( alpha(n_g1)>1./fac_sec ) ? 1.80 * C_RC122 * lambda_RC2 * formfunc_alpha_p  : 0.;

}

void Coalescence_bulles_2groupes_Smith::coefficient_WE(const DoubleTab& alpha, const DoubleTab& p, const DoubleTab& T,
                                                       const DoubleTab& rho, const DoubleTab& nu, const DoubleTab& sigma, const double Dh,
                                                       const DoubleTab& ndv, const DoubleTab& d_bulles,
                                                       const DoubleTab& eps, const DoubleTab& k_turb, const int n_l, const int n_g1, const int n_g2,
                                                       DoubleTab& coeff) const
{

  const double fac_sec = 1.e4 ;// numerical security
  const double alphafonction_3_over2 = std::sqrt(1.-alpha(n_g1)) * std::sqrt(1.-alpha(n_g1)) * std::sqrt(1.-alpha(n_g1)) ;
  const double alphafonction_9_over7 = std::pow(1.-alpha(n_g1), 9./7.) ;
  const double alphafonction_5_over2 =  std::sqrt(1.-alpha(n_g1)) * std::sqrt(1.-alpha(n_g1)) * std::sqrt(1.-alpha(n_g1)) * std::sqrt(1.-alpha(n_g1)) * std::sqrt(1.-alpha(n_g1)) ;
  const double C_D1 = 2./3. * d_bulles(n_g1) * std::sqrt( g * std::abs(rho(n_l) - rho(n_g1) ) / sigma(n_g1, n_l) ) * ( (1.+ 17.67 * alphafonction_9_over7 )/ (18.67 * alphafonction_3_over2 ) ) * ( (1.+ 17.67 * alphafonction_9_over7 ) / (18.67 * alphafonction_3_over2 ) );
  const double Ur1 = std::min(ndv(n_l,n_g1),std::sqrt(4./3. * d_bulles(n_g1) / C_D1 * g *std::abs(rho(n_l) - rho(n_g1) ) / rho(n_l) * (alpha(n_l))));
  const double D_crit = 4. * std::sqrt(sigma(n_g2, n_l) / g / (rho(n_l)-rho(n_g2)));
  const double C_D2 = 8./3. * (1.-alpha(n_g2) ) * (1.-alpha(n_g2));
  const double Ur2 = std::min(ndv(n_l,n_g2),std::sqrt(4./3. * std::max(d_bulles(n_g2),D_crit) / C_D2 * g *std::abs(rho(n_l) - rho(n_g2) ) / rho(n_l) * (alpha(n_l))));
  const double formfunc_alpha = std::max(1. - std::exp(-0.7 * alpha(n_g2) ) ,0.) * std::max(1.- 0.1 * ( D_crit/std::max(d_bulles(n_g2),D_crit))* (D_crit/std::max(d_bulles(n_g2),D_crit)) , 0.);
  const double Uw12 = std::max(0.94 * Ur2 * std::cbrt(C_D2) + Ur1 - Ur2 ,0.) ;


  // WE (1) coefficient for secmem---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  coeff(n_g1, n_l) = ( alpha(n_g1)>1./fac_sec ) ?   - 0.17 * C_WE1 * std::cbrt(C_D1) * Ur1   : 0.; //WE1

  // dWE (1) coefficient for mat---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  coeff(n_l, n_g1) = ( alpha(n_g1)>1./fac_sec ) ? - 0.17 * C_WE1 * ( Ur1 * 0.097282 * (3. * (17.67 * alphafonction_9_over7 + 1. ) / 2. / alphafonction_5_over2  - 22.7186 * std::pow(1.-alpha(n_g1), 17./14.) ) / std::cbrt((17.67 * alphafonction_9_over7 + 1.) / alphafonction_3_over2 ) ) : 0.; //dWE1/alpha1

  // WE (2) coefficient for secmem---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  coeff(n_g2, n_l) = ( alpha(n_g2)>1./fac_sec ) ?  - 1.02 * C_WE2 * (0.94 * Ur2 * std::cbrt(C_D2)) * formfunc_alpha   : 0.; //WE2
  coeff(n_l, n_g2) = ( alpha(n_g2)>1./fac_sec ) ? - 1.02 * C_WE2 * (0.94 * Ur2 * std::cbrt(C_D2)) * formfunc_alpha  : 0.;

  // WE (112) coefficient for 1st group for secmem---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  coeff(n_g1, n_g1) = ( alpha(n_g1)>1./fac_sec ) ? 2.57 * C_WE112 * C_D1 * Ur1 * std::max(1. - 2./3. * D_crit / d_bulles(n_g1),0.) : 0.;

  // WE (122) coefficient for 1st group for secmem---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  coeff(n_g1, n_g2) = ( alpha(n_g1)>1./fac_sec ) ?  -0.33 * C_WE112 * Uw12 : 0.;

  // WE (122) coefficient for 2nd group for secmem---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  coeff(n_g2, n_g1) =( alpha(n_g1)>1./fac_sec ) ?   0.922 * C_WE112 * Uw12 : 0.; //WE122 in 2



}

