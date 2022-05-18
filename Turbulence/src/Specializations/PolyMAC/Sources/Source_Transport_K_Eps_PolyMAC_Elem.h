/****************************************************************************
* Copyright (c) 2019, CEA
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
// File:        Source_Transport_K_Eps_PolyMAC_Elem.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/PolyMAC/Sources
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Source_Transport_K_Eps_PolyMAC_Elem_included
#define Source_Transport_K_Eps_PolyMAC_Elem_included

#define C1_DEFAULT 1.44   // Valeurs par defaut des constantes qui interviennent
#define C2_DEFAULT 1.92   // dans le calcul des termes sources des equations
#define C3_DEFAULT 1.0    // de transport de K et Eps source: Chabard et N3S

#include <Source_base.h>
#include <Ref_Zone_PolyMAC.h>
#include <Ref_Champ_Don.h>
#include <Ref_Champ_Don_base.h>
#include <Ref_Convection_Diffusion_Temperature.h>
#include <Ref_Convection_Diffusion_Concentration.h>
#include <Ref_Equation_base.h>
#include <Ref_Transport_K_Eps.h>
#include <Calcul_Production_K_PolyMAC.h>

class Probleme_base;
class Champ_Don_base;
#include <TRUSTTabs_forward.h>
#include <TRUSTTabs_forward.h>
class Zone_dis;
class Zone_Cl_dis;
class Zone_Cl_PolyMAC;
class Champ_Face;

//////////////////////////////////////////////////////////////////////////////
//.DESCRIPTION class Source_Transport_K_Eps_PolyMAC_Elem
//
// Cette classe represente le terme source qui figure dans l'equation
// de transport du couple (k,eps) dans le cas ou les equations de Navier-Stokes
// ne sont pas couplees a la thermique ou a l'equation de convection-diffusion
// d'une concentration.
//
//////////////////////////////////////////////////////////////////////////////

class Source_Transport_K_Eps_PolyMAC_Elem : public Source_base,
  public Calcul_Production_K_PolyMAC
{

  Declare_instanciable_sans_constructeur(Source_Transport_K_Eps_PolyMAC_Elem);

public:

  inline Source_Transport_K_Eps_PolyMAC_Elem(double cte1 = C1_DEFAULT,
                                             double cte2 = C2_DEFAULT );
  DoubleTab& ajouter(DoubleTab& ) const override;
  DoubleTab& calculer(DoubleTab& ) const override;
  void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const override ;
  void mettre_a_jour(double temps) override ;

protected:

  double C1;
  double C2;
  REF(Zone_PolyMAC) la_zone_PolyMAC;
  REF(Equation_base) eq_hydraulique;
  REF(Transport_K_Eps)  mon_eq_transport_K_Eps;

  void associer_pb(const Probleme_base& pb) override;
  void associer_zones(const Zone_dis& ,const Zone_Cl_dis& ) override;
};

//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Source_Transport_K_Eps_anisotherme_PolyMAC_Elem
//
// Cette classe represente le terme source qui figure dans l'equation
// de transport du couple (k,eps) dans le cas ou les equations de Navier_Stokes
// sont couplees a l'equation de la thermique
// On suppose que le coefficient de variation de la masse volumique
// du fluide en fonction de ce scalaire est un coefficient uniforme.
//
//////////////////////////////////////////////////////////////////////////////

class Source_Transport_K_Eps_anisotherme_PolyMAC_Elem :
  public Source_Transport_K_Eps_PolyMAC_Elem
{

  Declare_instanciable_sans_constructeur(Source_Transport_K_Eps_anisotherme_PolyMAC_Elem);

public:

  inline Source_Transport_K_Eps_anisotherme_PolyMAC_Elem(double cte1 = C1_DEFAULT,
                                                         double cte2 = C2_DEFAULT,
                                                         double cte3 = C3_DEFAULT);
  void associer_pb(const Probleme_base& ) override;
  DoubleTab& ajouter(DoubleTab& ) const override;
  DoubleTab& calculer(DoubleTab& ) const override;

protected:

  double C3;
  REF(Convection_Diffusion_Temperature) eq_thermique;
  REF(Champ_Don) beta_t;
  REF(Champ_Don_base) gravite;

};

//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Source_Transport_K_Eps_concen_PolyMAC_Elem
//
// Cette classe represente le terme source qui figure dans l'equation
// de transport du couple (k,eps) dans le cas ou les equations de Navier_Stokes
// sont couplees a l'equation d'un transport d'un constituant
// On suppose que le coefficient de variation de la masse volumique
// du fluide en fonction de ce scalaire est un coefficient uniforme.
//
//////////////////////////////////////////////////////////////////////////////

class Source_Transport_K_Eps_concen_PolyMAC_Elem :
  public Source_Transport_K_Eps_PolyMAC_Elem
{

  Declare_instanciable_sans_constructeur(Source_Transport_K_Eps_concen_PolyMAC_Elem);

public:

  inline Source_Transport_K_Eps_concen_PolyMAC_Elem(double cte1 = C1_DEFAULT,
                                                    double cte2 = C2_DEFAULT,
                                                    double cte3 = C3_DEFAULT);
  void associer_pb(const Probleme_base& ) override;
  DoubleTab& ajouter(DoubleTab& ) const override;
  DoubleTab& calculer(DoubleTab& ) const override;

protected:
  double C3;
  REF(Convection_Diffusion_Concentration) eq_concentration;
  REF(Champ_Don) beta_c;
  REF(Champ_Don_base) gravite;

};

//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Source_Transport_K_Eps_aniso_therm_concen_PolyMAC_Elem
//
// Cette classe represente le terme source qui figure dans l'equation
// de transport du couple (k,eps) dans le cas ou les equations de
// Navier_Stokes sont couplees a l'equation de convection diffusion
// d'une concentration et a l'equation de la thermique
// Les champs beta_t et beta_c sont uniformes
//
//////////////////////////////////////////////////////////////////////////////

class Source_Transport_K_Eps_aniso_therm_concen_PolyMAC_Elem :
  public Source_Transport_K_Eps_PolyMAC_Elem
{

  Declare_instanciable_sans_constructeur(Source_Transport_K_Eps_aniso_therm_concen_PolyMAC_Elem);

public:

  inline Source_Transport_K_Eps_aniso_therm_concen_PolyMAC_Elem(double cte1 = C1_DEFAULT,
                                                                double cte2 = C2_DEFAULT,
                                                                double cte3 = C3_DEFAULT);
  void associer_pb(const Probleme_base& ) override;
  DoubleTab& ajouter(DoubleTab& ) const override;
  DoubleTab& calculer(DoubleTab& ) const override;

protected:

  double C3;
  REF(Convection_Diffusion_Temperature) eq_thermique;
  REF(Convection_Diffusion_Concentration) eq_concentration;
  REF(Champ_Don) beta_t;
  REF(Champ_Don) beta_c;
  REF(Champ_Don_base) gravite;

};


//////////////////////////////////////////////////////////////////////////////
//
//   Fonctions inline de la classe Source_Transport_K_Eps_PolyMAC_Elem
//
//////////////////////////////////////////////////////////////////////////////

inline Source_Transport_K_Eps_PolyMAC_Elem::
Source_Transport_K_Eps_PolyMAC_Elem(double cte1,double cte2)

  : C1(cte1), C2(cte2) {}


//////////////////////////////////////////////////////////////////////////////
//
//                     Fonctions inline de la classe
//
//             Source_Transport_K_Eps_anisotherme_PolyMAC_Elem
//
//////////////////////////////////////////////////////////////////////////////

inline Source_Transport_K_Eps_anisotherme_PolyMAC_Elem::
Source_Transport_K_Eps_anisotherme_PolyMAC_Elem(double cte1,double cte2,double cte3)

  : Source_Transport_K_Eps_PolyMAC_Elem(cte1,cte2) , C3(cte3) {}

//////////////////////////////////////////////////////////////////////////////
//
//                     Fonctions inline de la classe
//
//             Source_Transport_K_Eps_concen_PolyMAC_Elem
//
//////////////////////////////////////////////////////////////////////////////

inline Source_Transport_K_Eps_concen_PolyMAC_Elem::
Source_Transport_K_Eps_concen_PolyMAC_Elem(double cte1,double cte2,double cte3)

  : Source_Transport_K_Eps_PolyMAC_Elem(cte1,cte2) , C3(cte3) {}

//////////////////////////////////////////////////////////////////////////////
//
//                        Fonctions inline de la classe
//
//                Source_Transport_K_Eps_aniso_therm_concen_PolyMAC_Elem
//
//////////////////////////////////////////////////////////////////////////////

inline Source_Transport_K_Eps_aniso_therm_concen_PolyMAC_Elem::
Source_Transport_K_Eps_aniso_therm_concen_PolyMAC_Elem(double cte1,double cte2,double cte3)

  : Source_Transport_K_Eps_PolyMAC_Elem(cte1,cte2) , C3(cte3) {}

#endif