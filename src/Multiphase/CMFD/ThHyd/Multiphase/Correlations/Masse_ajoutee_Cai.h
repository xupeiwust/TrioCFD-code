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

#ifndef Masse_ajoutee_Cai_included
#define Masse_ajoutee_Cai_included

#include <Masse_ajoutee_base.h>
#include <TRUSTTabs_forward.h>

/*! @brief Masse ajoutee de la forme ma(k, l) = +/- beta * alpha_k * alpha_l * rho_m
 *
 *     avec beta un coefficient constant (0.5 par defaut) et rho_m la masse volumique du melange
 *
 *
 */
class Masse_ajoutee_Cai : public Masse_ajoutee_base
{
  Declare_instanciable(Masse_ajoutee_Cai);

public:
  void ajouter(const double *alpha, const double *rho, DoubleTab& a_r) const override;
  void ajouter_inj(const double *flux_alpha, const double *alpha, const double *rho, DoubleTab& f_a_r) const override;
  void coeff(const DoubleTab& alpha, const DoubleTab& rho, DoubleTab& coeff) const override;

protected:
  double beta = 0.5;
  int n_l = -1; //liquid phase
  double inj_ajoutee_liquide_ = 1.;
  double inj_ajoutee_gaz_ = 1.;
};

#endif
