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
// File:        Transport_K_ou_Eps.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Incompressible/Equations/RANS
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Transport_K_ou_Eps_included
#define Transport_K_ou_Eps_included

#include <Transport_K_ou_Eps_base.h>
#include <Operateur_Conv.h>


class Motcle;

/*! @brief classe Transport_K_ou_Eps Cette classe represente l'equation de transport de l'energie cinetique
 *
 *     turbulente K et du taux de dissipation epsilon (eps) associee au modele
 *     de turbulence (k,eps).
 *     On traite en une seule equation le transport des deux
 *     grandeurs turbulentes. Il s'agit donc d'une equation vectorielle, dont
 *     le champ inconnue possede 2 composantes: K et epsilon.
 *     On instanciera un objet de cette classe uniquement si
 *     on utilise le modele (k,eps) pour traiter la turbulence en mode "bicephale"
 *
 * @sa Equation_base
 */
class Transport_K_ou_Eps : public Transport_K_ou_Eps_base
{

  Declare_instanciable(Transport_K_ou_Eps);

public:
  void set_param(Param& param) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  inline void associer_Champ_Inconnu(const Champ_Inc_base& );
  void associer_modele_turbulence(const Modele_turbulence_hyd_2_eq_base& ) override;
  int nombre_d_operateurs() const override;
  const Operateur& operateur(int) const override;
  Operateur& operateur(int) override;
  const Motcle& domaine_application() const override;
  DoubleTab& corriger_derivee_impl(DoubleTab& d) override;

protected :
  int with_nu_;

  OBS_PTR(Champ_Inc_base) inco_eqn_associee;
  OWN_PTR(Champ_Don_base) Champ_don_nul_;  // on y met 0 si on ne veut pas de nu

};


/*! @brief Associe un champ de vitesse (transportante)a l'equation.
 *
 * @param (Champ_Inc_base& vit) le champ de vitesse a associer a l'equation
 */
inline void Transport_K_ou_Eps::associer_Champ_Inconnu(const Champ_Inc_base& vit)
{
  inco_eqn_associee = vit;
}

#endif


