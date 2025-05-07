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

#ifndef Flux_2groupes_Smith_included
#define Flux_2groupes_Smith_included

#include <Flux_2groupes_base.h>
#include <TRUSTTabs_forward.h>
#include <TRUSTTab.h>

/*! @brief Flux interfacial a coefficient constant par phase
 *
 */

class Flux_2groupes_Smith : public Flux_2groupes_base
{
  Declare_instanciable(Flux_2groupes_Smith);
public:
  void coeffs(const input_coeffs& input, output_coeffs& output) const override;
  void therm(const input_therms& input, output_therms& output) const override;
protected:
  const double g = 9.81;
  const double b = 0.597;
  const double p = 0.01;
  const double alpha_max = 0.62 ;
  const double C_RC1 = 0.005 ;
  const double C_RC0 = 3.0 ;
  const double C_RC122 = 0.005 ;
  const double C_WE1 = 0.002 ;
  const double C_TI21 = 0.02 ;
  const double C_S0 = 0.000038 ;
  const double We_SO = 4500. ;
  const double We_cr1 = 6.5 ;
  const double We_cr2 = 7.0 ;
  double Xi_h_ = 0.;
  double A_c_ = 1.;
  double hPNVG_ = 0.;
  double betac=0.4;
  double pc = 0.055;
  double R_ = 462. ;
  double theta_rad = M_PI/2. ;

};

#endif
