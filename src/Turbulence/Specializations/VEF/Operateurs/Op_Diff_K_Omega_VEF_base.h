/****************************************************************************
* Copyright (c) 2023, CEA
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
/*! @brief class Op_Diff_K_Omega_VEF_Face Cette classe represente l'operateur de diffusion turbulente
 *
 *    La discretisation est VEF
 *   Les methodes pour l'implicite sont codees.
 *
 */
#ifndef Op_Diff_K_Omega_VEF_base_included
#define Op_Diff_K_Omega_VEF_base_included

#define PRDT_K_DEFAUT 1/0.6 // from Wilcox STD model, 1/Prandtl_K_ = sigma_k = 0.6
#define PRDT_OMEGA_DEFAUT 1/0.5 // from Wilcox STD model, 1/Prandtl_Omega_ = sigma_omega = 0.5

#include <Op_Dift_VEF_base.h>
#include <Op_VEF_Face.h>
#include <Modele_turbulence_hyd_K_Omega.h>


//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Op_Diff_K_Omega_VEF_base
//
//////////////////////////////////////////////////////////////////////////////

class Op_Diff_K_Omega_VEF_base : public Op_Dift_VEF_base
{

  Declare_instanciable_sans_constructeur(Op_Diff_K_Omega_VEF_base);

public:

  Op_Diff_K_Omega_VEF_base(double Prandt_K = PRDT_K_DEFAUT ,
                           double Prandt_Omega = PRDT_OMEGA_DEFAUT ) : Prdt_K(Prandt_K) , Prdt_Omega(Prandt_Omega)
  { }

  void completer() override;


protected:
  double Prdt_K;
  double Prdt_Omega;

  // For SST model
  static constexpr double SIGMA_K1 = 0.85;
  static constexpr double SIGMA_K2 = 1.0;
  static constexpr double SIGMA_OMEGA1 = 0.5;
  static constexpr double SIGMA_OMEGA2 = 0.856;

  OBS_PTR(Modele_turbulence_hyd_K_Omega) turbulence_model;
};

#endif
