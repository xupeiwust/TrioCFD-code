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


#include <Rupture_bulles_2groupes_Smith.h>
#include <Pb_Multiphase.h>

Implemente_instanciable(Rupture_bulles_2groupes_Smith, "Rupture_bulles_2groupes_Smith", Rupture_bulles_2groupes_base);

Sortie& Rupture_bulles_2groupes_Smith::printOn(Sortie& os) const
{
  return os;
}

Entree& Rupture_bulles_2groupes_Smith::readOn(Entree& is)
{
  return is;
}

void Rupture_bulles_2groupes_Smith::coefficient_TI(const DoubleTab& alpha, const DoubleTab& p, const DoubleTab& T,
                                                   const DoubleTab& rho, const DoubleTab& nu, const DoubleTab& sigma, const double Dh,
                                                   const DoubleTab& ndv, const DoubleTab& d_bulles,
                                                   const DoubleTab& eps, const DoubleTab& k_turb, const int n_l, const int n_g1, const int n_g2,
                                                   DoubleTab& coeff) const
{

  const double fac_sec = 1.e4 ; // numerical security
  const double We1 = 2. * rho(n_l) * ( std::cbrt(eps(n_l)*d_bulles(n_g1)) * std::cbrt(eps(n_l)*d_bulles(n_g1)) )  * d_bulles(n_g1) / sigma(n_g1, n_l) ;
  const double D_crit = 4. * std::sqrt(sigma(n_g1, n_l) / g / (rho(n_l)-rho(n_g1)));
  const double We2 = 2. * rho(n_l) * ( std::cbrt(eps(n_l)*std::max(d_bulles(n_g2),D_crit)) * std::cbrt(eps(n_l)*std::max(d_bulles(n_g2),D_crit)) ) * std::max(d_bulles(n_g2),D_crit) / sigma(n_g2, n_l) ;
  const double Dcrit_overD2_13_over3 = std::cbrt(D_crit/std::max(d_bulles(n_g2),D_crit)) * std::cbrt(D_crit/std::max(d_bulles(n_g2),D_crit)) * std::cbrt(D_crit/std::max(d_bulles(n_g2),D_crit)) * std::cbrt(D_crit/std::max(d_bulles(n_g2),D_crit)) * std::cbrt(D_crit/std::max(d_bulles(n_g2),D_crit)) * std::cbrt(D_crit/std::max(d_bulles(n_g2),D_crit)) * std::cbrt(D_crit/std::max(d_bulles(n_g2),D_crit)) * std::cbrt(D_crit/std::max(d_bulles(n_g2),D_crit)) * std::cbrt(D_crit/std::max(d_bulles(n_g2),D_crit)) * std::cbrt(D_crit/std::max(d_bulles(n_g2),D_crit)) * std::cbrt(D_crit/std::max(d_bulles(n_g2),D_crit)) * std::cbrt(D_crit/std::max(d_bulles(n_g2),D_crit)) * std::cbrt(D_crit/std::max(d_bulles(n_g2),D_crit)) ;
  const double Dcrit_overD2_5 = (D_crit/std::max(d_bulles(n_g2),D_crit)) * (D_crit/std::max(d_bulles(n_g2),D_crit)) * (D_crit/std::max(d_bulles(n_g2),D_crit)) * (D_crit/std::max(d_bulles(n_g2),D_crit)) * (D_crit/std::max(d_bulles(n_g2),D_crit)) ;

  // TI (1) coefficient for secmem-------------------------------------------------------------------------------------------------------------------------------------------------------------------

  coeff(n_g1, n_l) = ( alpha(n_g1)>1./fac_sec ) ?  0.12 * C_TI1 * std::exp(-We_cr1/We1) * std::sqrt(1.-std::min(We_cr1/We1 ,1.)) : 0.;
  coeff(n_l, n_g1) =  coeff(n_g1, n_l);

  // TI (2) coefficient for secmem--------------------------------------------------------------------------------------------------------------------------------------------------------------------

  coeff(n_g2, n_l) = ( alpha(n_g2)>1./fac_sec ) ? 0.378 * C_TI2 * std::exp(-We_cr2/We2) * std::sqrt(1.-std::min(We_cr1/We2 ,1.))*std::max(1.-0.212 * Dcrit_overD2_13_over3,0.) : 0. ;
  coeff(n_l, n_g2) = coeff(n_g2, n_l);

  // TI (21) coefficient for secmem------------------------------------------------------------------------------------------------------------------------------------------------------------------

  coeff(n_g1, n_g2) = ( alpha(n_g2)>1./fac_sec ) ? 6.165 * C_TI21 * std::exp(-We_cr2/We2) * std::sqrt(1.-std::min(We_cr1/We2 ,1.))*std::max(0.212 * Dcrit_overD2_13_over3 - 0.167 * Dcrit_overD2_5,0.) : 0.;
  coeff(n_g2, n_g1) = coeff(n_g1, n_g2);

}
void Rupture_bulles_2groupes_Smith::coefficient_SI(const DoubleTab& alpha, const DoubleTab& p, const DoubleTab& T,
                                                   const DoubleTab& rho, const DoubleTab& nu, const DoubleTab& sigma, const double Dh,
                                                   const DoubleTab& ndv, const DoubleTab& d_bulles,
                                                   const DoubleTab& eps, const DoubleTab& k_turb, const int n_l, const int n_g1, const int n_g2,
                                                   DoubleTab& coeff) const
{

  const double La_square = sigma(n_g2, n_l) / g / (rho(n_l)-rho(n_g2));
  const double C_D2 = 8./3. * (1.-std::min(alpha(n_g2),0.81)) * (1.-std::min(alpha(n_g2),0.81));
  const double D_crit = 4. * std::sqrt(sigma(n_g2, n_l) / g / (rho(n_l)-rho(n_g1)));
  const double Ur2 = std::min(ndv(n_l,n_g2),std::sqrt(4./3. * std::max(d_bulles(n_g2),D_crit) / C_D2 * g *std::abs(rho(n_l) - rho(n_g2) ) / rho(n_l) * (alpha(n_l))));
  const double fac_sec = 1.e4 ; // numerical security

  //SI coefficient for secmem----------------------------------------------------------------------------------------------------------------------------------------------------------------------

  coeff(n_g2, n_l) = ( alpha(n_g2)>1./fac_sec ) ? 2.616e-4 * C_RC2 * std::cbrt(eps(n_l)) / Dh / Dh * std::pow(La_square , 1./6.) * std::max(1-std::exp(-C_RC0 * std::sqrt(alpha(n_g2))),0.) + 1.425e-7 * C_WE2 * (1.-std::exp(-0.7 * alpha(n_g2))) * 0.94 * Ur2 * std::cbrt(C_D2) / La_square : 0. ;
  coeff(n_l, n_g2) = coeff(n_g2, n_l);

}
void Rupture_bulles_2groupes_Smith::coefficient_SO(const DoubleTab& alpha, const DoubleTab& p, const DoubleTab& T,
                                                   const DoubleTab& rho, const DoubleTab& nu, const DoubleTab& sigma, const double Dh,
                                                   const DoubleTab& ndv, const DoubleTab& d_bulles,
                                                   const DoubleTab& eps, const DoubleTab& k_turb, const int n_l,  const int n_g1,  const int n_g2,
                                                   DoubleTab& coeff) const
{
  const double fac_sec = 1.e4 ; // numerical security
  const double D_crit = 4. * std::sqrt(sigma(n_g2, n_l) / g / (rho(n_l)-rho(n_g1)));
  const double We2 = 2. * rho(n_l) * (std::cbrt(eps(n_l)*std::max(d_bulles(n_g2),D_crit)) * std::cbrt(eps(n_l)*std::max(d_bulles(n_g2),D_crit)) ) * std::max(d_bulles(n_g2),D_crit) / sigma(n_g2, n_l) ;
  const double Ur2 = std::max(ndv(n_l,n_g2),std::sqrt(1./2. * D_crit *  g *std::abs(rho(n_l) - rho(n_g2) ) / rho(n_l) ));


  //SO1 coefficient for secmem-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

  coeff(n_g1, n_l) = ( alpha(n_g2)>1./fac_sec ) ?  8.0 * C_S0 * std::pow(rho(n_l) , 3./5.) * std::pow(ndv(n_l,n_g2),1./5.) * std::pow(sigma(n_g2, n_l) , 2./5.) / rho(n_g2) / std::pow(Dh,2./5.) /  std::pow(We_cr2,3./5.) * std::max(1. - (We_SO / We2 * We_SO / We2 * We_SO / We2 * We_SO / We2 ),0.) : 0.; // SO1
  coeff(n_l, n_g1) =  coeff(n_g1, n_l);

  //SO2 coefficient for secmem-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

  coeff(n_g2, n_l) = ( alpha(n_g2)>1./fac_sec ) ? -0.36 * C_S0 * sigma(n_g2, n_l) / rho(n_g2) / Ur2 * std::max(1. - We_SO / We2 ,0.) : 0.; // SO2
  coeff(n_l,n_g2) = coeff(n_g2, n_l)  ;


}
