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

#ifndef Coalescence_bulles_1groupe_Yao_Morel_included
#define Coalescence_bulles_1groupe_Yao_Morel_included

#include <TRUSTTabs_forward.h>
#include <Coalescence_bulles_1groupe_base.h>
#include <TRUSTTab.h>

/*! @brief Model for bubble coalescence from Yao and Morel (2003)
 *
 * Model for bubble coalescence using the turbulence induced coalescence given in
 * Yao, W. and Morel, C. (2003) Volumetric interfacial area prediction in upward bubbly
 * two-phase ﬂow. International Journal of Heat and Mass Transfer 47 (2004) 307--328.
 * The coefficient formula is given at page 312, equation (22).
 */
class Coalescence_bulles_1groupe_Yao_Morel: public Coalescence_bulles_1groupe_base
{
  Declare_instanciable(Coalescence_bulles_1groupe_Yao_Morel);

public:
  // TODO change to input_t/output_t format
  void coefficient(const DoubleTab& alpha, const DoubleTab& p, const DoubleTab& T,
                   const DoubleTab& rho, const DoubleTab& nu, const DoubleTab& sigma, double Dh,
                   const DoubleTab& ndv, const DoubleTab& d_bulles,
                   const DoubleTab& eps, const DoubleTab& k_turb,
                   DoubleTab& coeff) const override;

private:
  int n_l = -1;

  static constexpr double Kc1 = 2.86;
  static constexpr double Kc2 = 1.922;
  static constexpr double Kc3 = 1.017;
  const double alpha_max_1_3 = std::cbrt(M_PI/6.);
  static constexpr double We_cr = 1.24;
  static constexpr double alpha_sec = 2./3.;


};

#endif
