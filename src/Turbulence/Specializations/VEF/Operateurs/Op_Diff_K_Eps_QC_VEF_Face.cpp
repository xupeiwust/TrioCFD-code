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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Op_Diff_K_Eps_QC_VEF_Face.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Operateurs
//
//////////////////////////////////////////////////////////////////////////////

#include <Op_Diff_K_Eps_QC_VEF_Face.h>
#include <Modele_turbulence_hyd_K_Eps.h>
#include <Champ_P1NC.h>
#include <Periodique.h>
#include <Paroi_hyd_base_VEF.h>
#include <Fluide_Quasi_Compressible.h>
#include <TRUSTTrav.h>

Implemente_instanciable(Op_Diff_K_Eps_QC_VEF_Face,"Op_Diff_K_Eps_QC_VEF_P1NC",Op_Diff_K_Eps_VEF_Face);

////  printOn
//


Sortie& Op_Diff_K_Eps_QC_VEF_Face::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

//// readOn
//

Entree& Op_Diff_K_Eps_QC_VEF_Face::readOn(Entree& s )
{
  return s ;
}


DoubleTab& Op_Diff_K_Eps_QC_VEF_Face::ajouter(const DoubleTab& inconnue, DoubleTab& resu) const
{
  DoubleTrav KEps_divided_by_rho(inconnue);
  const Transport_K_Eps& eqn_transport = ref_cast(Transport_K_Eps,mon_equation.valeur());
  const Fluide_Quasi_Compressible& mil = ref_cast(Fluide_Quasi_Compressible,eqn_transport.milieu());
  const DoubleTab& rho=mil.masse_volumique().valeurs();
  int size = inconnue.dimension_tot(0);
  for (int i=0; i<size; i++)
    {
      KEps_divided_by_rho(i,0) = inconnue(i,0) / rho(i); // k/rho
      KEps_divided_by_rho(i,1) = inconnue(i,1) / rho(i); // eps/rho
    }
  return Op_Diff_K_Eps_VEF_Face::ajouter(KEps_divided_by_rho, resu);
}
