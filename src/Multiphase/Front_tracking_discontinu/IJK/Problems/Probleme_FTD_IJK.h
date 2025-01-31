/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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

#ifndef Probleme_FTD_IJK_included
#define Probleme_FTD_IJK_included

#include <Probleme_FTD_IJK_base.h>

/*! @brief : class Probleme_FTD_IJK
 *
 *  <Description of class Probleme_FTD_IJK>
 *
 *
 *  La classe Probleme_FTD_IJK herite de la classe Probleme_FTD_IJK_base.
 *
 */
class Probleme_FTD_IJK : public Probleme_FTD_IJK_base
{
  friend class IJK_Thermique;
  friend class Statistiques_dns_ijk_FT;
  Declare_instanciable(Probleme_FTD_IJK) ;

public :

  bool run() override;
  void euler_time_step(ArrOfDouble& var_volume_par_bulle) override;
  void rk3_sub_step(const int rk_step, const double total_timestep, const double fractionnal_timestep, const double time) override;
};

#endif /* Probleme_FTD_IJK_included */
