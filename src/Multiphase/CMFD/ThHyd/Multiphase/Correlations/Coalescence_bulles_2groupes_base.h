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

#ifndef Coalescence_bulles_2groupes_base_included
#define Coalescence_bulles_2groupes_base_included

#include <TRUSTTabs_forward.h>
#include <Correlation_base.h>
#include <TRUSTTab.h>

/*! @brief
 *
 */

class Coalescence_bulles_2groupes_base : public Correlation_base
{
  Declare_base(Coalescence_bulles_2groupes_base);
public:
  virtual void coefficient_WE(const DoubleTab& alpha, const DoubleTab& p, const DoubleTab& T,
                              const DoubleTab& rho, const DoubleTab& nu, const DoubleTab& sigma, const double Dh,
                              const DoubleTab& ndv, const DoubleTab& d_bulles,
                              const DoubleTab& eps, const DoubleTab& k_turb, const int n_l, const int n_g1, const int n_g2,
                              DoubleTab& coeff) const  = 0;
  virtual void coefficient_RC(const DoubleTab& alpha, const DoubleTab& alpha_p, const DoubleTab& p, const DoubleTab& T,
                              const DoubleTab& rho, const DoubleTab& nu, const DoubleTab& sigma, double Dh,
                              const DoubleTab& ndv, const DoubleTab& d_bulles,
                              const DoubleTab& eps, const DoubleTab& k_turb, const int n_l, const int n_g1, const int n_g2,
                              DoubleTab& coeff) const  = 0;

};

#endif
