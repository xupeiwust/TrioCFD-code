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
// File:        Pb_Hydraulique_Melange_Binaire_Turbulent_QC.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Quasi_Compressible
// Version:     /main/15
//
//////////////////////////////////////////////////////////////////////////////

#include <Pb_Hydraulique_Melange_Binaire_Turbulent_QC.h>
#include <Verif_Cl.h>
#include <Les_mod_turb.h>
#include <Verif_Cl_Turb.h>
#include <Modifier_nut_pour_QC.h>

Implemente_instanciable(Pb_Hydraulique_Melange_Binaire_Turbulent_QC,"Pb_Hydraulique_Melange_Binaire_Turbulent_QC",Pb_QC_base);

// XD pb_hydraulique_melange_binaire_turbulent_qc Pb_base pb_hydraulique_melange_binaire_turbulent_qc -1 Resolution of turbulent binary mixture problem under low Mach number.
// XD attr navier_stokes_turbulent_qc navier_stokes_turbulent_qc navier_stokes_turbulent_qc 0 Navier-Stokes equations under low Mach number as well as the associated turbulence model equations.
// XD attr convection_diffusion_fraction_massique_mb_turbulent_qc convection_diffusion_fraction_massique_mb_turbulent_qc convection_diffusion_fraction_massique_mb_turbulent_qc 0 Species conservation equation under low Mach number as well as the associated turbulence model equations.

// Description:
//    Simple appel a: Probleme_base::printOn(Sortie&)
//    Ecrit le probleme sur un flot de sortie.
// Precondition:
// Parametre: Sortie& os
//    Signification: un flot de sortie
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Sortie&
//    Signification: le flot de sortie modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
Sortie& Pb_Hydraulique_Melange_Binaire_Turbulent_QC::printOn(Sortie& os) const
{
  return Probleme_base::printOn(os);
}

// Description:
//    Simple appel a: Probleme_base::readOn(Entree&)
//    Lit le probleme a partir d'un flot d'entree.
// Precondition:
// Parametre: Entree& is
//    Signification: un flot d'entree
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Entree&
//    Signification: le flot d'entree modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
Entree& Pb_Hydraulique_Melange_Binaire_Turbulent_QC::readOn(Entree& is)
{
  return Probleme_base::readOn(is);
}

// Description:
//    Renvoie le nombre d'equation,
//    Renvoie 2 car il y a 2 equations a un probleme de
//    Hydraulique Melange Binaire turbulent :
//        l'equation de Navier Stokes
//        l' equation de conv/diff fraction massique
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: int
//    Signification: le nombre d'equation
//    Contraintes: toujours 2 car il y a 2 equations au probleme
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
int Pb_Hydraulique_Melange_Binaire_Turbulent_QC::nombre_d_equations() const
{
  return 2;
}

// Description:
//    Renvoie l'equation d'hydraulique de type Navier_Stokes_Turbulent_QC si i=0
//    Renvoie l'equation de conv/diff fraction massique de type
//    Convection_Diffusion_fraction_massique_Turbulent_QC si i=1
//    (version const)
// Precondition:
// Parametre: int i
//    Signification: l'index de l'equation a renvoyer
//    Valeurs par defaut:
//    Contraintes: 0 <= i <= 1
//    Acces:
// Retour: Equation_base&
//    Signification: l'equation correspondante a l'index
//    Contraintes: reference constante
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
const Equation_base& Pb_Hydraulique_Melange_Binaire_Turbulent_QC::equation(int i) const
{
  if ( !( i==0 || i==1 ) )
    {
      Cerr << "\nError in Pb_Hydraulique_Melange_Binaire_Turbulent_QC::equation() : Wrong number of equation !" << finl;
      Process::exit();
    }
  if (i == 0)
    return eq_hydraulique;
  else
    return eq_frac_mass;
}

// Description:
//    Renvoie l'equation d'hydraulique de type Navier_Stokes_Turbulent_QC si i=0
//    Renvoie l'equation de conv/diff fraction massique de type
//    Convection_Diffusion_fraction_massique_Turbulent_QC si i=1
// Precondition:
// Parametre: int i
//    Signification: l'index de l'equation a renvoyer
//    Valeurs par defaut:
//    Contraintes: 0 <= i <= 1
//    Acces:
// Retour: Equation_base&
//    Signification: l'equation correspondante a l'index
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
Equation_base& Pb_Hydraulique_Melange_Binaire_Turbulent_QC::equation(int i)
{
  if ( !( i==0 || i==1 ) )
    {
      Cerr << "\nError in Pb_Hydraulique_Melange_Binaire_Turbulent_QC::equation() : Wrong number of equation !" << finl;
      Process::exit();
    }
  if (i == 0)
    return eq_hydraulique;
  else
    return eq_frac_mass;
}

// Description:
//    Teste la compatibilite des equations de la thermique
//    et de l'hydraulique. Le test se fait sur les conditions
//    aux limites discretisees de chaque equation.
//    Appel la fonction de librairie hors classe:
//      tester_compatibilite_hydr_fraction_massique(const Zone_Cl_dis&,const Zone_Cl_dis&)
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: int
//    Signification: renvoie toujours 1
//    Contraintes:
// Exception: modeles de turbulence de famille differente pour
//            l'hydraulique et la thermique
// Effets de bord:
// Postcondition:
int Pb_Hydraulique_Melange_Binaire_Turbulent_QC::verifier()
{
  const Zone_Cl_dis& zone_Cl_hydr = eq_hydraulique.zone_Cl_dis();
  const Zone_Cl_dis& zone_Cl_fm = eq_frac_mass.zone_Cl_dis();
  // Verification de la compatibilite des conditions aux limites:
  tester_compatibilite_hydr_fraction_massique(zone_Cl_hydr,zone_Cl_fm);

  if ( sub_type(Mod_turb_hyd_RANS, eq_hydraulique.get_modele(TURBULENCE).valeur() ) )
    {
      const Mod_turb_hyd_RANS& le_mod_RANS = ref_cast(Mod_turb_hyd_RANS, eq_hydraulique.get_modele(TURBULENCE).valeur());
      const Transport_K_Eps_base& eqn = ref_cast(Transport_K_Eps_base, le_mod_RANS.eqn_transp_K_Eps());
      const Zone_Cl_dis& zone_Cl_turb = eqn.zone_Cl_dis();
      tester_compatibilite_hydr_turb(zone_Cl_hydr, zone_Cl_turb);
    }

  // Verification de la compatibilite des modeles de turbulence:
  const Mod_turb_hyd& le_mod_turb_hyd = eq_hydraulique.modele_turbulence();
  const Modele_turbulence_scal_base& le_mod_turb_th =
    ref_cast(Modele_turbulence_scal_base,eq_frac_mass.get_modele(TURBULENCE).valeur());

  if  (sub_type(Modele_turbulence_hyd_K_Eps,le_mod_turb_hyd.valeur()))
    {
      if (!sub_type(Modele_turbulence_scal_Prandtl,le_mod_turb_th))
        {
          Cerr << "Les modeles de turbulence ne sont pas de la meme famille" << finl;
          Cerr << "pour l'hydraulique et la fraction massique" << finl;
          Process::exit();
        }
    }

  return 1;
}

int Pb_Hydraulique_Melange_Binaire_Turbulent_QC::expression_predefini(const Motcle& motlu, Nom& expression)
{
  if (motlu=="VISCOSITE_TURBULENTE")
    {
      expression  = "predefini { pb_champ ";
      expression += le_nom();
      expression += " viscosite_turbulente } ";
      return 1;
    }
  else
    return Pb_QC_base::expression_predefini(motlu,expression);
  return 0;
}

