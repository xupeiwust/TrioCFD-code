/****************************************************************************
* Copyright (c) 2020, CEA
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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Pb_Hydraulique_sensibility.cpp
// Directory : $Sensitivity_analysis/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Pb_Hydraulique_sensibility.h>
#include <Fluide_Incompressible.h>

Implemente_instanciable(Pb_Hydraulique_sensibility,"Pb_Hydraulique_sensibility",Pb_Fluide_base);
// XD Pb_Hydraulique_sensibility Pb_base Pb_Hydraulique_sensibility -1 Resolution of hydraulic sensibility problems
// XD  attr fluide_incompressible fluide_incompressible fluide_incompressible 0 The fluid medium associated with the problem.
// XD  attr Navier_Stokes_standard_sensibility Navier_Stokes_standard_sensibility Navier_Stokes_standard_sensibility 0 Navier-Stokes sensibility equations


/*! @brief Simple appel a: Pb_Fluide_base::printOn(Sortie&) Ecrit le probleme sur un flot de sortie.
 *
 * @param (Sortie& os) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Pb_Hydraulique_sensibility::printOn(Sortie& os) const
{
  return Pb_Fluide_base::printOn(os);
}


/*! @brief Simple appel a: Pb_Fluide_base::readOn(Entree&) Lit le probleme a partir d'un flot d'entree.
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Pb_Hydraulique_sensibility::readOn(Entree& is)
{
  return Pb_Fluide_base::readOn(is);
}

/*! @brief Renvoie le nombre d'equation, Renvoie 1 car il y a seulement 1 equation a un probleme
 *
 *     hydraulique standard: l'equation de Navier Stokes
 *
 * @return (int) le nombre d'equation
 */
int Pb_Hydraulique_sensibility::nombre_d_equations() const
{
  return 1;
}

/*! @brief Renvoie l'equation d'hydraulique de type Navier_Stokes_std si i=0, sort (exit) sinon.
 *
 *     (version const)
 *
 * @param (int i) l'index de l'equation a renvoyer
 * @return (Equation_base&) l'equation d'hydraulique de type Navier_Stokes_std
 */
const Equation_base& Pb_Hydraulique_sensibility::equation(int i) const
{
  if (i!=0)
    {
      Cerr << "\nError in Pb_Hydraulique::equation() : The problem has only one equation !" << finl;
      Process::exit();
    }
  return eq_hydraulique;
}

/*! @brief Renvoie l'equation d'hydraulique de type Navier_Stokes_std si i=0, sort (exit) sinon.
 *
 * @param (int i) l'index de l'equation a renvoyer
 * @return (Equation_base&) l'equation d'hydraulique de type Navier_Stokes_std
 */
Equation_base& Pb_Hydraulique_sensibility::equation(int i)
{
  if (i!=0)
    {
      Cerr << "\nPb_Hydraulique::equation() : The problem has only one equation !" << finl;
      Process::exit();
    }
  return eq_hydraulique;
}



/*! @brief Associe le milieu au probleme.
 *
 * Le milieu doit etre de type fluide incompressible.
 *
 * @param (Milieu_base& mil) le milieu physique a associer au probleme
 * @throws le milieu n'est pas du type Fluide_Incompressible
 */
void Pb_Hydraulique_sensibility::associer_milieu_base(const Milieu_base& mil)
{
  if (sub_type(Fluide_Incompressible,mil))
    eq_hydraulique.associer_milieu_base(mil);
  else
    {
      Cerr << "Un milieu de type " << mil.que_suis_je() << " ne peut etre associe a "
           << finl;
      Cerr << "un probleme de type Pb_Hydraulique_sensibility " << finl;
      Process::exit();
    }
}
