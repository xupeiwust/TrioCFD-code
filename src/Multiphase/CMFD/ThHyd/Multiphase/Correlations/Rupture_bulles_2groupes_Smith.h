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

#ifndef Rupture_bulles_2groupes_Smith_included
#define Rupture_bulles_2groupes_Smith_included

#include <TRUSTTabs_forward.h>
#include <Rupture_bulles_2groupes_base.h>
#include <TRUSTTab.h>

/*! @brief
 *
 */

class Rupture_bulles_2groupes_Smith : public Rupture_bulles_2groupes_base
{
  Declare_instanciable(Rupture_bulles_2groupes_Smith);
public:
  void coefficient_TI(const DoubleTab& alpha, const DoubleTab& p, const DoubleTab& T,
                      const DoubleTab& rho, const DoubleTab& nu, const DoubleTab& sigma, const double Dh,
                      const DoubleTab& ndv, const DoubleTab& d_bulles,
                      const DoubleTab& eps, const DoubleTab& k_turb, const int n_l, const int n_g1, const int n_g2,
                      DoubleTab& coeff) const override;
  void coefficient_SI(const DoubleTab& alpha, const DoubleTab& p, const DoubleTab& T,
                      const DoubleTab& rho, const DoubleTab& nu, const DoubleTab& sigma, const double Dh,
                      const DoubleTab& ndv, const DoubleTab& d_bulles,
                      const DoubleTab& eps, const DoubleTab& k_turb, const int n_l, const int n_g1, const int n_g2,
                      DoubleTab& coeff) const override;

  void coefficient_SO(const DoubleTab& alpha, const DoubleTab& p, const DoubleTab& T,
                      const DoubleTab& rho, const DoubleTab& nu, const DoubleTab& sigma, const double Dh,
                      const DoubleTab& ndv, const DoubleTab& d_bulles,
                      const DoubleTab& eps, const DoubleTab& k_turb, const int n_l, const int n_g1, const int n_g2,
                      DoubleTab& coeff) const override;

private:
  const double C_TI1 = 0.1 ;
  const double C_TI21 = 0.02 ;
  const double We_cr1 = 6.5 ;
  const double C_TI2 = 0.02 ;
  const double We_cr2 = 7.0 ;
  const double C_RC2 = 0.005 ;
  const double C_WE2 = 0.005 ;
  const double C_RC0 = 3.0 ;
  const double C_S0 = 0.000038 ;
  const double We_SO = 4500. ;
  const double g= 9.81 ;

};

#endif
