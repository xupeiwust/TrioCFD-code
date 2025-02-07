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
// File:        Op_Diff_K_Eps_VEF_Face.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Operateurs
//
//////////////////////////////////////////////////////////////////////////////

#include <Op_Diff_K_Eps_VEF_Face.h>
#include <Paroi_hyd_base_VEF.h>
#include <Debog.h>

Implemente_instanciable(Op_Diff_K_Eps_VEF_Face,"Op_Diff_K_Eps_VEF_P1NC",Op_Diff_K_Eps_VEF_base);

Sortie& Op_Diff_K_Eps_VEF_Face::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}
Entree& Op_Diff_K_Eps_VEF_Face::readOn(Entree& s )
{
  return s ;
}

DoubleTab& Op_Diff_K_Eps_VEF_Face::ajouter(const DoubleTab& inconnue_org, DoubleTab& resu) const
{
  remplir_nu(nu_); // On remplit le tableau nu car ajouter peut se faire avant le premier pas de temps

  const DoubleTab& nu_turb = diffusivite_turbulente().valeurs();
  const int nb_comp = resu.line_size();

  // On dimensionne et initialise le tableau des bilans de flux:
  flux_bords_.resize(le_dom_vef->nb_faces_bord(), nb_comp);
  flux_bords_ = 0.;

  int n_tot = nu_.dimension_tot(0); //TODO peut mieux faire
  DoubleTab nu_turb_m(n_tot, 2);
  for (int k=0; k<n_tot; k++)
    {
      nu_turb_m(k,0) = nu_turb(k)/Prdt[0];
      nu_turb_m(k,1) = nu_turb(k)/Prdt[1];
    }

  ajouter_bord_gen<Type_Champ::SCALAIRE, true>(inconnue_org, resu, flux_bords_, nu_, nu_turb_m);
  ajouter_interne_gen<Type_Champ::SCALAIRE, true>(inconnue_org, resu, flux_bords_, nu_, nu_turb_m);

  modifier_flux(*this);

  return resu;
}

void Op_Diff_K_Eps_VEF_Face::contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& matrice) const
{
  modifier_matrice_pour_periodique_avant_contribuer(matrice, equation());
  remplir_nu(nu_); // On remplit le tableau nu car l'assemblage d'une matrice avec ajouter_contribution peut se faire avant le premier pas de temps

  const DoubleTab& nu_turb = diffusivite_turbulente().valeurs();

  int n_tot = nu_.dimension_tot(0);
  DoubleTab nu_turb_m(n_tot, 2);
  for (int k=0; k<n_tot; k++)
    {
      nu_turb_m(k,0) = nu_turb(k)/Prdt[0];
      nu_turb_m(k,1) = nu_turb(k)/Prdt[1];
    }

  int marq = phi_psi_diffuse(equation());

  DoubleVect porosite_eventuelle(equation().milieu().porosite_face());
  if (!marq) porosite_eventuelle = 1;

  ajouter_contribution_bord_gen<Type_Champ::SCALAIRE, false, true>(inco, matrice, nu_, nu_turb_m, porosite_eventuelle);
  ajouter_contribution_interne_gen<Type_Champ::SCALAIRE, false, true>(inco, matrice, nu_, nu_turb_m, porosite_eventuelle);

  modifier_matrice_pour_periodique_apres_contribuer(matrice, equation());
}


void Op_Diff_K_Eps_VEF_Face::modifier_pour_Cl(Matrice_Morse& matrice, DoubleTab& secmem) const
{
  Op_Dift_VEF_base::modifier_pour_Cl(matrice, secmem);

  const Turbulence_paroi_base& mod=le_modele_turbulence->loi_paroi();
  const Paroi_hyd_base_VEF& paroi=ref_cast(Paroi_hyd_base_VEF,mod);
  const ArrOfInt& face_keps_imposee=paroi.face_keps_imposee();
  int size=secmem.dimension(0);
  const IntVect& tab1=matrice.get_tab1();
  DoubleVect& coeff = matrice.get_set_coeff();
  const DoubleTab& val=equation().inconnue().valeurs();
  const int nb_comp = equation().inconnue().valeurs().line_size();

  if (face_keps_imposee.size_array()>0) //TODO a reformuler (Kokkos ?)
    // en plus des dirichlets ????
    // on change la matrice et le resu sur toutes les lignes ou k_eps_ est imposee....
    for (int face=0; face<size; face++)
      {
        if (face_keps_imposee[face]!=-2)
          {
            for (int comp=0; comp<nb_comp; comp++)
              {
                // on doit remettre la ligne a l'identite et le secmem a l'inconnue
                int idiag = tab1[face*nb_comp+comp]-1;
                coeff[idiag]=1;
                // pour les voisins
                int nbvois = tab1[face*nb_comp+1+comp] - tab1[face*nb_comp+comp];
                for (int k=1; k < nbvois; k++)
                  {
                    coeff[idiag+k]=0;
                  }
                secmem(face,comp)=val(face,comp);
              }
          }
      }

}

