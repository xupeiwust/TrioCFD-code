/****************************************************************************
* Copyright (c) 2021, CEA
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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Production_HZDR_PolyMAC_P0.h
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/PolyMAC_P0
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Production_HZDR_PolyMAC_P0_included
#define Production_HZDR_PolyMAC_P0_included

#include <Source_base.h>
#include <TRUSTTab.h>
#include <TRUST_Ref.h>


/*! @brief Classe Production_HZDR_PolyMAC_P0 : Cette classe implemente, dans PolyMAC_P0, la production de BIA dans l'equation de k.
 * L'agitation induite par les bulles (BIA) est prise en compte dans l'équation de k grace à un terme de production.
 * Ce modèle est donné par Roland Rzehak et Eckhard Krepper (2013) du HZDR (d'où le nom donné au modèle).
 *
 * @sa Operateur_PolyMAC_P0_base Operateur_base
 */
class Production_HZDR_PolyMAC_P0: public Source_base
{
  Declare_instanciable(Production_HZDR_PolyMAC_P0);
public :
  int has_interface_blocs() const override
  {
    return 1;
  }
  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = {}) const override;
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {}) const override;
  void check_multiphase_compatibility() const override {}; //of course

  void associer_domaines(const Domaine_dis_base& ,const Domaine_Cl_dis_base& ) override { };
  void associer_pb(const Probleme_base& ) override { };
  void mettre_a_jour(double temps) override { };

protected:
  double g_ = 9.81;
  double C_k = 1.;
  int n_l = -1 ;
};

#endif
