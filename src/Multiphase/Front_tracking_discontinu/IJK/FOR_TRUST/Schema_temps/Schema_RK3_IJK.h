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

#ifndef Schema_RK3_IJK_included
#define Schema_RK3_IJK_included

#include <Schema_Temps_IJK_base.h>

class Schema_RK3_IJK: public Schema_Temps_IJK_base
{
  Declare_instanciable(Schema_RK3_IJK);
public :
  int faire_un_pas_de_temps_eqn_base(Equation_base&) override;

  double& get_store_RK3_source_acc() { return store_RK3_source_acc_; }
  double get_store_RK3_source_acc() const { return store_RK3_source_acc_; }

  double& get_store_RK3_fac_sv() { return store_RK3_fac_sv_; }
  double get_store_RK3_fac_sv() const { return store_RK3_fac_sv_; }

  double& get_current_time_at_rk3_step() { return current_time_at_rk3_step_; }
  double get_current_time_at_rk3_step() const { return current_time_at_rk3_step_; }

  int& get_rk_step() { return rk_step_; }
  int get_rk_step() const { return rk_step_; }

protected:
  double store_RK3_source_acc_ = 0., store_RK3_fac_sv_ = 1.;
  double current_time_at_rk3_step_ = 0.;
  int rk_step_ = -1; // default value
};

#endif /* Schema_RK3_IJK_included */
