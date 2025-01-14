/****************************************************************************
* Copyright (c) 2015, CEA
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
// File:        Traitement_particulier_NS_canal_VDF.h
// Directory:   $TRIO_U_ROOT/src/VDF/Turbulence
// Version:     /main/16
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Traitement_particulier_NS_canal_VDF_included
#define Traitement_particulier_NS_canal_VDF_included

#include <Traitement_particulier_NS_canal.h>

class Navier_Stokes_Turbulent;

/*! @brief classe Traitement_particulier_NS_canal_VDF Cette classe permet de faire les traitements particuliers
 *
 *      pour le calcul d'un canal plan :
 *          * conservation du debit
 *          * calculs de moyennes
 *
 *
 * @sa Navier_Stokes_Turbulent, Traitement_particulier_base,, Traitement_particulier_VDF
 */
class Traitement_particulier_NS_canal_VDF : public Traitement_particulier_NS_canal
{
  Declare_instanciable(Traitement_particulier_NS_canal_VDF);

public :

  Entree& lire(const Motcle& , Entree& );
  Entree& lire(Entree& );

protected :
  void remplir_Y(DoubleVect&, DoubleVect&, int& ) const;
  void remplir_Tab_recap(IntTab&) const;
  void calculer_moyenne_spatiale_vitesse_rho_mu(DoubleTab&) const;
  void calculer_moyenne_spatiale_nut(DoubleTab&) const;
  void calculer_moyenne_spatiale_Temp(DoubleTab&) const;
  void calculer_Temp(DoubleTab&,const DoubleTab&) const;
};
#endif
