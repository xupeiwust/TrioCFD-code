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

#ifndef Fluide_Diphasique_IJK_included
#define Fluide_Diphasique_IJK_included

#include <Fluide_Diphasique.h>

class Fluide_Diphasique_IJK: public Fluide_Diphasique
{
  Declare_instanciable(Fluide_Diphasique_IJK);
public:

  /*
   * ATTENTION : phase0_ => VAPEUR, phase1_ => LIQUIDE
   */
  double get_mu_liquid() const
  {
    return ref_cast(Fluide_Incompressible, phase1_.valeur()).viscosite_dynamique().valeurs()(0,0);
  }

  double get_mu_vapour() const
  {
    return ref_cast(Fluide_Incompressible, phase0_.valeur()).viscosite_dynamique().valeurs()(0,0);;
  }

  double get_delta_rho() const
  {
    /* rho_l - rho_v */
    return (ref_cast(Fluide_Incompressible, phase1_.valeur()).masse_volumique().valeurs()(0,0) -
            ref_cast(Fluide_Incompressible, phase0_.valeur()).masse_volumique().valeurs()(0,0));
  }

  double get_rho_liquid() const
  {
    return ref_cast(Fluide_Incompressible, phase1_.valeur()).masse_volumique().valeurs()(0,0);
  }

  double get_rho_vapour() const
  {
    return ref_cast(Fluide_Incompressible, phase0_.valeur()).masse_volumique().valeurs()(0,0);
  }

  double get_gravite_norm() const;

  int get_direction_gravite() const { return direction_gravite_; }
  void calculate_direction_gravite();

protected:
  int direction_gravite_ = 0;
};

#endif /* Fluide_Diphasique_IJK_included */
