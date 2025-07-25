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
// File:        Source_Transport_K_Omega_VEF_Face.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Sources
//
//////////////////////////////////////////////////////////////////////////////

#include <Discretisation_tools.h>
#include <Operateur_Grad.h>
#include <TRUSTTabs_forward.h>
#include <Source_Transport_K_Omega_VEF_Face.h>
#include <Modele_turbulence_hyd_K_Omega.h>
#include <Transport_K_Omega.h>
#include <Navier_Stokes_std.h>
#include <Pb_Hydraulique_Turbulent.h>
#include <Milieu_base.h>
#include <VEF_discretisation.h>
#include <Domaine_VEF.h>
#include <Domaine_Cl_VEF.h>

Implemente_instanciable_sans_constructeur(Source_Transport_K_Omega_VEF_Face,
                                          "Source_Transport_K_Omega_VEF_P1NC",
                                          Source_Transport_K_Omega_VEF_Face_base);

Sortie& Source_Transport_K_Omega_VEF_Face::printOn(Sortie& s) const
{
  return s << que_suis_je();
}

Entree& Source_Transport_K_Omega_VEF_Face::readOn(Entree& is)
{
  Source_Transport_K_Omega_VEF_Face_base::verifier_pb_komega(mon_equation->probleme(),
                                                             que_suis_je());
  Source_Transport_K_Omega_VEF_Face_base::readOn(is);
  return (is);
}

void Source_Transport_K_Omega_VEF_Face::get_noms_champs_postraitables(Noms& nom,
                                                                      Option opt) const
{
  Source_Transport_K_Omega_VEF_Face_base::get_noms_champs_postraitables(nom, opt);
  Noms noms_compris = champs_compris_.liste_noms_compris();
  noms_compris.add("grad_k_grad_omega");

  if (opt == DESCRIPTION)
    Cerr << " Source_Transport_K_Omega_VEF_Face : " << noms_compris << finl;
  else
    nom.add(noms_compris);
}

void Source_Transport_K_Omega_VEF_Face::creer_champ(const Motcle& nom)
{
  Source_Transport_K_Omega_VEF_Face_base::creer_champ(nom);

  if (grad_k_omega_.est_nul())
    {
      const VEF_discretisation& disc = ref_cast(VEF_discretisation, equation().discretisation());
      Noms noms(1), unites(1);
      noms[0] = "grad_k_grad_omega";
      disc.discretiser_champ("champ_elem", equation().domaine_dis(), scalaire,
                             noms , unites, 1, 0, grad_k_omega_);
      champs_compris_.ajoute_champ(grad_k_omega_);
    }
}

void Source_Transport_K_Omega_VEF_Face::associer_pb(const Probleme_base& pb)
{
  Source_Transport_K_Omega_VEF_Face_base::associer_pb(pb);
  eqn_K_Omega = ref_cast(Transport_K_Omega, equation());
  turbulence_model = ref_cast(Modele_turbulence_hyd_K_Omega, eqn_K_Omega->modele_turbulence());
}

const DoubleTab& Source_Transport_K_Omega_VEF_Face::get_visc_turb() const
{
  return eqn_K_Omega->modele_turbulence().viscosite_turbulente().valeurs();
}

const DoubleTab& Source_Transport_K_Omega_VEF_Face::get_cisaillement_paroi() const
{
  return turbulence_model->loi_paroi().Cisaillement_paroi();
}

const DoubleTab& Source_Transport_K_Omega_VEF_Face::get_K_pour_production() const
{
  return eqn_K_Omega->inconnue().valeurs();
}

const Nom Source_Transport_K_Omega_VEF_Face::get_type_paroi() const
{
  return turbulence_model->loi_paroi().que_suis_je();
}

void Source_Transport_K_Omega_VEF_Face::compute_blending_F1(DoubleTab& gradKgradOmega) const
{
  const DoubleTab& K_Omega = eqn_K_Omega->inconnue().valeurs();
  const DoubleTab& kinematic_viscosity = get_visc_turb();
  const DoubleTab& distmin = le_dom_VEF->y_faces(); // Minimum distance to the edge

  DoubleTab& tabF1 = ref_cast_non_const(DoubleTab, turbulence_model->get_tabF1());
  DoubleTab& tabF2 = ref_cast_non_const(DoubleTab, turbulence_model->get_tabF2());

  DoubleTab visc_face(le_dom_VEF->nb_faces_tot()); // dimension_tot(0)
  // Discretisation_tools::cells_to_faces(le_dom_VEF.valeur(), kinematic_viscosity, visc_face);
  elem_to_face(le_dom_VEF.valeur(), kinematic_viscosity, visc_face);

  for (int face = 0; face < le_dom_VEF->nb_faces(); face++)
    {
      double const dmin = std::max(distmin(face), 1e-12);
      double const enerK = K_Omega(face, 0);
      double const omega = K_Omega(face, 1);

      double const tmp1 = sqrt(enerK)/(BETA_K*omega*dmin);
      double const tmp2 = 500.0*visc_face(face)/(omega*dmin*dmin);
      double const maxval = std::max(2*SIGMA_OMEGA2*gradKgradOmega(face)/omega, 1e-20);
      double const tmp3 = 4.0*SIGMA_OMEGA2*enerK/(maxval*dmin*dmin);

      double const arg1 = std::min(std::max(tmp1, tmp2), tmp3); // Common name of the variable
      tabF1(face) = std::tanh(arg1*arg1*arg1*arg1);

      double const arg2 = std::max(2.*tmp1, tmp2);
      tabF2(face) = std::tanh(arg2*arg2);
    }
}

double Source_Transport_K_Omega_VEF_Face::blender(double const val1, double const val2,
                                                  int const face) const
{
  const DoubleTab& F1 = turbulence_model->get_tabF1();
  return F1(face)*val1 + (1 - F1(face))*val2;
}


void Source_Transport_K_Omega_VEF_Face::compute_cross_diffusion(DoubleTab& gradKgradOmega) const
{
  // Same structure than Source_WC_Chaleur_VEF.cpp in TRUST

  // = First, store K and Omega in two separated Tab for the gradient computation, not ideal
  const DoubleTab& K_Omega = eqn_K_Omega->inconnue().valeurs();
  DoubleTab enerK; // field on faces
  DoubleTab omega; // field on faces
  enerK.resize(K_Omega.dimension_tot(0));
  omega.resize(K_Omega.dimension_tot(0));
  for (int num_face = 0; num_face < le_dom_VEF->nb_faces_tot(); ++num_face)
    {
      enerK(num_face) = K_Omega(num_face, 0);
      omega(num_face) = K_Omega(num_face, 1);
      gradKgradOmega(num_face) = 0;
    }

  // = Definition of the gradients
  // Get number of componants, for gradient resize
  // We use the velocity field to resize the gradient as the velocity is on faces
  const Navier_Stokes_std& eqHyd = ref_cast(Navier_Stokes_std, mon_equation->probleme().equation(0));
  const DoubleTab& velocity_field_face = eqHyd.vitesse().valeurs(); // Velocity on faces
  const int nbr_velocity_components = velocity_field_face.dimension(1); // dimension_tot ?

  DoubleTab gradK_elem; // field on elements
  DoubleTab gradOmega_elem; // field on elements
  gradK_elem.resize(le_dom_VEF->nb_elem_tot(), nbr_velocity_components);
  gradOmega_elem.resize(le_dom_VEF->nb_elem_tot(), nbr_velocity_components);

  // Compute the two gradients
  const Operateur_Grad& Op_Grad_komega = eqn_K_Omega->gradient_operator_komega();
  Op_Grad_komega.calculer(enerK, gradK_elem);
  Op_Grad_komega.calculer(omega, gradOmega_elem);

  DoubleTab& chmp_post = ref_cast_non_const(DoubleTab, grad_k_omega_->valeurs());
  for (int num_elem = 0; num_elem < le_dom_VEF->nb_elem_tot(); ++num_elem)
    for (int ncompo = 0; ncompo < nbr_velocity_components; ++ncompo)
      chmp_post(num_elem) += gradK_elem(num_elem, ncompo) * gradOmega_elem(num_elem, ncompo);

  // Interpolate from elem to face
  const Domaine_dis_base& domaine_dis = mon_equation->inconnue().domaine_dis_base();
  const Domaine_VF& domaine = ref_cast(Domaine_VF, domaine_dis);

  DoubleTab gradK_face;
  gradK_face.resize(velocity_field_face.dimension_tot(0), nbr_velocity_components);
  Discretisation_tools::cells_to_faces(domaine, gradK_elem, gradK_face);

  DoubleTab gradOmega_face;
  gradOmega_face.resize(velocity_field_face.dimension_tot(0), nbr_velocity_components);
  Discretisation_tools::cells_to_faces(domaine, gradOmega_elem, gradOmega_face);

  // Dot Product gradKgradOmega
  for (int num_face = 0; num_face < le_dom_VEF->nb_faces_tot(); ++num_face)
    for (int ncompo = 0; ncompo < nbr_velocity_components; ++ncompo)
      gradKgradOmega(num_face) +=
        gradK_face(num_face, ncompo) * gradOmega_face(num_face, ncompo);
}

void Source_Transport_K_Omega_VEF_Face::fill_resu(const DoubleVect& volumes_entrelaces,
                                                  const DoubleTrav& ProdK,
                                                  const DoubleTab& gradKgradOmega,
                                                  DoubleTab& resu) const
{
  const DoubleTab& K_Omega = eqn_K_Omega->inconnue().valeurs();
  const double LeK_MIN = eqn_K_Omega->modele_turbulence().get_K_MIN();

  for (int face = 0; face < le_dom_VEF->nb_faces(); face++)
    {
      const double tke = K_Omega(face, 0);
      const double omega = K_Omega(face, 1);

      resu(face, 0) += (ProdK(face) - BETA_K*tke*omega)*volumes_entrelaces(face);

      if (tke >= LeK_MIN)
        {
          const double cALPHA = turbulence_model->is_SST()
                                ? blender(GAMMA1, GAMMA2, face)
                                : ALPHA_OMEGA;
          const double cBETA = turbulence_model->is_SST()
                               ? blender(BETA1, BETA2, face)
                               : BETA_OMEGA;
          const double cSIGMA = turbulence_model->is_SST()
                                ? 2*(1 - turbulence_model->get_tabF1()(face)*SIGMA_OMEGA2)
                                : (gradKgradOmega(face) > 0)*1/8;

          double contrib {0};
          contrib += cALPHA*ProdK(face)*omega/tke; // production
          contrib += - cBETA*omega*omega; // dissipation
          contrib += cSIGMA/omega*gradKgradOmega(face); // cross diffusion
          contrib *= volumes_entrelaces(face);
          resu(face, 1) += contrib;
        }
    }
}

DoubleTab& Source_Transport_K_Omega_VEF_Face::ajouter(DoubleTab& resu) const
{
  return Source_Transport_K_Omega_VEF_Face_base::ajouter_komega(resu);
}

void Source_Transport_K_Omega_VEF_Face::contribuer_a_avec(const DoubleTab& a,
                                                          Matrice_Morse& matrice) const
{
  const double LeK_MIN = eqn_K_Omega->modele_turbulence().get_K_MIN();
  const DoubleTab& K_Omega = equation().inconnue().valeurs();
  const DoubleVect& porosite_face = eqn_K_Omega->milieu().porosite_face();
  const DoubleVect& volumes_entrelaces = le_dom_VEF->volumes_entrelaces();

  DoubleTrav production_TKE {le_dom_VEF->nb_faces()};

  const Domaine_Cl_VEF& domaine_Cl_VEF = ref_cast(Domaine_Cl_VEF,
                                                  eq_hydraulique->domaine_Cl_dis());
  const DoubleTab& visco_turb = get_visc_turb(); // voir les classes filles
  const DoubleTab& velocity = eq_hydraulique->inconnue().valeurs();
  const DoubleTab& TKE = get_K_pour_production(); // voir les classes filles
  calculer_terme_production_K(le_dom_VEF.valeur(), domaine_Cl_VEF, production_TKE,
                              TKE, velocity, visco_turb,
                              _interpolation_viscosite_turbulente,
                              _coefficient_limiteur);

  DoubleTrav gradKgradOmega(le_dom_VEF->nb_faces_tot());
  compute_cross_diffusion(gradKgradOmega);

  for (int face = 0; face < K_Omega.dimension(0); face++)
    if (K_Omega(face, 0) >= LeK_MIN)
      {
        const double tke = K_Omega(face, 0);
        const double omega = K_Omega(face, 1);

        const double volporo = porosite_face(face) * volumes_entrelaces(face);

        const double coef_k = BETA_K*omega*volporo; // K_Omega(face, 1)/K_Omega(face, 0)*volporo;
        matrice(face*2, face*2) += coef_k;

        const double cALPHA = turbulence_model->is_SST()
                              ? blender(GAMMA1, GAMMA2, face)
                              : ALPHA_OMEGA;

        const double cBETA = turbulence_model->is_SST()
                             ? blender(BETA1, BETA2, face)
                             : BETA_OMEGA;

        const double cSIGMA = turbulence_model->is_SST()
                              ? 2*(1 - turbulence_model->get_tabF1()(face)*SIGMA_OMEGA2)
                              : (gradKgradOmega(face) > 0)*1/8;

        const double coef_omega = (-cALPHA*production_TKE(face)/tke
                                   + cBETA*omega
                                   - cSIGMA/(omega*omega)*gradKgradOmega(face) ) * volporo;
        matrice(face*2 + 1, face*2 + 1) += coef_omega;
      }
}
