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
// File:        Source_Transport_K_Realisable_VEF_Face.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Sources
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_K_Eps_Realisable_Bicephale.h>
#include <Source_Transport_K_Realisable_VEF_Face.h>
#include <Schema_Temps_base.h>
#include <Champ_Uniforme.h>
#include <Fluide_base.h>
#include <Champ_P1NC.h>
#include <Domaine_VEF.h>
#include <Debog.h>

Implemente_instanciable(Source_Transport_K_Realisable_VEF_Face,"Source_Transport_K_Realisable_VEF_P1NC",Source_base);

Sortie& Source_Transport_K_Realisable_VEF_Face::printOn(Sortie& s ) const { return s << que_suis_je() ; }

Entree& Source_Transport_K_Realisable_VEF_Face::readOn(Entree& is ) { return Source_Transport_Realisable_VEF_Face_base::readOn_nothing(is,que_suis_je()); }

void Source_Transport_K_Realisable_VEF_Face::associer_pb(const Probleme_base& pb)
{
  Source_Transport_Realisable_VEF_Face_base::associer_pb(pb);
  eqn_k_Rea = ref_cast(Transport_K_ou_Eps_Realisable, equation());
  const Modele_turbulence_hyd_K_Eps_Realisable_Bicephale& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Realisable_Bicephale, eqn_eps_Rea->modele_turbulence());
  eqn_eps_Rea = ref_cast(Transport_K_ou_Eps_Realisable, mod_turb.eqn_transp_Eps());
}

DoubleTab& Source_Transport_K_Realisable_VEF_Face::ajouter(DoubleTab& resu) const
{
  Debog::verifier("Source_Transport_K_Realisable_VEF_Face::ajouter resu 0", resu);

  const DoubleTab& K_Rea = eqn_k_Rea->inconnue().valeurs(), &eps_Rea = eqn_eps_Rea->inconnue().valeurs();
  const Modele_turbulence_hyd_K_Eps_Realisable_Bicephale& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Realisable_Bicephale, eqn_k_Rea->modele_turbulence());
  const DoubleTab& visco_turb = mod_turb.viscosite_turbulente().valeurs(), &vit = eq_hydraulique->inconnue().valeurs();
  const DoubleVect& vol_ent = le_dom_VEF->volumes_entrelaces();

  DoubleTab vitesse_filtree(vit);
  ref_cast(Champ_P1NC,eq_hydraulique->inconnue()).filtrer_L2(vitesse_filtree);

  const int nb_faces = le_dom_VEF->nb_faces();
  DoubleTrav P(nb_faces);

  calculer_terme_production_K_BiK(le_dom_VEF.valeur(), le_dom_Cl_VEF.valeur(), P, K_Rea, eps_Rea, vitesse_filtree, visco_turb, _interpolation_viscosite_turbulente, _coefficient_limiteur);

  Debog::verifier("Source_Transport_K_Realisable_VEF_Face::ajouter P 0", P);

  for (int num_face = 0; num_face < nb_faces; num_face++)
    resu(num_face) += (P(num_face) - eps_Rea(num_face)) * vol_ent(num_face);

  return resu;
}

void Source_Transport_K_Realisable_VEF_Face::mettre_a_jour(double temps)
{
  Modele_turbulence_hyd_K_Eps_Realisable_Bicephale& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Realisable_Bicephale, eqn_k_Rea->modele_turbulence());
  Modele_Fonc_Realisable_base& mon_modele_fonc = mod_turb.associe_modele_fonction();
  const DoubleTab& visco_turb = mod_turb.viscosite_turbulente().valeurs();
  const DoubleTab& vit = eq_hydraulique->inconnue().valeurs();
  const double epsilon_minimum = eqn_k_Rea->modele_turbulence().get_EPS_MIN();
  const Champ_Don_base& ch_visco_cin = ref_cast(Fluide_base,eqn_k_Rea->milieu()).viscosite_cinematique();
  const DoubleTab& tab_visco = ch_visco_cin.valeurs();

  /*Paroi*/
  DoubleTab visco_tab(visco_turb.dimension_tot(0));
  assert(sub_type(Champ_Uniforme,ch_visco_cin));
  visco_tab = tab_visco(0, 0);
  const int idt = eq_hydraulique->schema_temps().nb_pas_dt();
  const DoubleTab& tab_paroi = mod_turb.loi_paroi().Cisaillement_paroi();

  const Domaine_Cl_dis_base& zcl_keps = eqn_k_Rea->domaine_Cl_dis();
  const Domaine_dis_base& domaine_dis_keps = eqn_k_Rea->domaine_dis();
  const DoubleTab& K_Rea = eqn_k_Rea->inconnue().valeurs(), & eps_Rea = eqn_eps_Rea->inconnue().valeurs();
  mon_modele_fonc.Contributions_Sources_Paroi_BiK(domaine_dis_keps, zcl_keps, vit, K_Rea, eps_Rea, epsilon_minimum, visco_tab, visco_turb, tab_paroi, idt);

  Calcul_Production_K_VEF::mettre_a_jour(temps);
}

void Source_Transport_K_Realisable_VEF_Face::fill_coeff_matrice(const int face, const DoubleVect& porosite_face, const DoubleVect& volumes_entrelaces, const double visco, Matrice_Morse& matrice) const
{
  const DoubleTab& K_Rea = eqn_k_Rea->inconnue().valeurs(), &eps_Rea = eqn_eps_Rea->inconnue().valeurs();
  const Modele_turbulence_hyd_K_Eps_Realisable_Bicephale& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Realisable_Bicephale, eqn_k_Rea->modele_turbulence());
  const double LeK_MIN = mod_turb.get_K_MIN(), LeEPS_MIN = mod_turb.get_EPS_MIN();

  if ((K_Rea(face) >= LeK_MIN) && (eps_Rea(face) >= LeEPS_MIN))
    matrice(face, face) += porosite_face(face) * volumes_entrelaces(face) * eps_Rea(face) / (K_Rea(face) + sqrt(visco * eps_Rea(face)));
}
