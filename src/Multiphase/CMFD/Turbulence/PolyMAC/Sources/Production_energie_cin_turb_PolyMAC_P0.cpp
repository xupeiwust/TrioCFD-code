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
// File:        Production_energie_cin_turb_PolyMAC_P0.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/PolyMAC_P0
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Production_energie_cin_turb_PolyMAC_P0.h>

#include <Op_Diff_Turbulent_PolyMAC_P0_Face.h>
#include <Echelle_temporelle_turbulente.h>
#include <Taux_dissipation_turbulent.h>
#include <Viscosite_turbulente_base.h>
#include <Domaine_PolyMAC_P0.h>
#include <Navier_Stokes_std.h>
#include <Pb_Multiphase.h>

Implemente_instanciable(Production_energie_cin_turb_PolyMAC_P0,"Production_energie_cin_turb_Elem_PolyMAC_P0", Source_Production_energie_cin_turb);

Sortie& Production_energie_cin_turb_PolyMAC_P0::printOn(Sortie& os) const {return Source_Production_energie_cin_turb::printOn(os);}
Entree& Production_energie_cin_turb_PolyMAC_P0::readOn(Entree& is) { return Source_Production_energie_cin_turb::readOn(is);}

void Production_energie_cin_turb_PolyMAC_P0::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Domaine_PolyMAC_P0&             domaine = ref_cast(Domaine_PolyMAC_P0, equation().domaine_dis());
  const Probleme_base&                       pb = ref_cast(Probleme_base, equation().probleme());
  const Navier_Stokes_std&               eq_qdm = ref_cast(Navier_Stokes_std, pb.equation(0));
  const Viscosite_turbulente_base&    visc_turb = ref_cast(Viscosite_turbulente_base, ref_cast(Op_Diff_Turbulent_PolyMAC_P0_Face, eq_qdm.operateur(0).l_op_base()).correlation());
  const DoubleVect& pe = equation().milieu().porosite_elem(), &ve = domaine.volumes();

  std::string Type_diss = ""; // omega or tau dissipation
  for (int i = 0 ; i < equation().probleme().nombre_d_equations() ; i++)
    {
      if      sub_type(Echelle_temporelle_turbulente, equation().probleme().equation(i)) Type_diss = "tau";
      else if sub_type(Taux_dissipation_turbulent, equation().probleme().equation(i)) Type_diss = "omega";
    }

  int Nph = pb.get_champ("vitesse").valeurs().dimension(1), nb_elem = domaine.nb_elem(), D = dimension, nf_tot = domaine.nb_faces_tot() ;
  int N = equation().inconnue().valeurs().line_size(),
      Np = equation().probleme().get_champ("pression").valeurs().line_size(),
      Na = sub_type(Pb_Multiphase,equation().probleme()) ? equation().probleme().get_champ("alpha").valeurs().line_size() : 1;
  int e, n, mp;

  double limiter_ = visc_turb.limiteur();
  double nut_l = -10000., fac;
  const bool is_multi = sub_type(Pb_Multiphase,equation().probleme());
  const DoubleTab& tab_rho = equation().probleme().get_champ("masse_volumique").passe(),
                   *palp = is_multi ? &equation().probleme().get_champ("alpha").passe() : nullptr,
                    *alp = is_multi ? &equation().probleme().get_champ("alpha").valeurs() : nullptr,
                     &nu =  equation().probleme().get_champ("viscosite_cinematique").passe(),
                      &k = equation().probleme().get_champ("k").valeurs(),
                       &tab_grad = pb.get_champ("gradient_vitesse").passe(),
                        *diss = equation().probleme().has_champ(Type_diss) ? &equation().probleme().get_champ(Type_diss).valeurs() : nullptr,
                         *pdiss = equation().probleme().has_champ(Type_diss) ? &equation().probleme().get_champ(Type_diss).passe() : nullptr;

  const int cnu = nu.dimension(0) == 1;
  const Champ_base&   ch_alpha_rho = is_multi ? ref_cast(Pb_Multiphase,equation().probleme()).equation_masse().champ_conserve() : equation().milieu().masse_volumique();

  if (Type_diss == "")
    {
      DoubleTrav nut(0, Nph);
      MD_Vector_tools::creer_tableau_distribue(eq_qdm.pression().valeurs().get_md_vector(), nut); //Necessary to compare size in eddy_viscosity()
      visc_turb.eddy_viscosity(nut);

      Matrice_Morse *mat = matrices.count("k") ? matrices.at("k") : nullptr;

      for( e = 0 ; e < nb_elem ; e++)
        for( n = 0; n<N ; n++)
          {
            double secmem_en = 0.;
            for (int d_U = 0; d_U < D; d_U++)
              for (int d_X = 0; d_X < D; d_X++)
                secmem_en += ( tab_grad(nf_tot + d_U + e * D , D * n + d_X) + tab_grad(nf_tot + d_X + e * D , D * n + d_U) ) * tab_grad(nf_tot + d_X + e * D , D * n + d_U) ;
            secmem_en *= pe(e) * ve(e) * (palp ? (*palp)(e, n) : 1.0) * tab_rho(e, n) * nut(e, n) ;

            secmem(e, n) += std::max(secmem_en, 0.);//  secmem(e, n) += fac_(e, n, 0) * std::max(secmem_en, 0.);

            if (mat) (*mat)(N * e + n, N * e + n) -= 0.;//fac_(e, n, 1) * std::max(secmem_en, 0.);
          }
    }

  else
    {
      for( e = 0 ; e < nb_elem ; e++)
        for(n = 0, mp = 0; n<N ; n++, mp += (Np > 1))
          {
            double grad_grad = 0.;
            for (int d_U = 0; d_U < D; d_U++)
              for (int d_X = 0; d_X < D; d_X++)
                grad_grad += ( tab_grad(nf_tot + d_U + e * D , D * n + d_X) + tab_grad(nf_tot + d_X + e * D , D * n + d_U) ) * tab_grad(nf_tot + d_X + e * D , D * n + d_U) ;

            fac = std::max(grad_grad, 0.) * pe(e) * ve(e) ;

            if      (Type_diss == "tau")   nut_l = k(e, n) * (*diss)(e, n) + limiter_ * nu(!cnu * e, n);
            else if (Type_diss == "omega") nut_l = k(e, n) / std::max((*pdiss)(e, n), omega_min_) + 0.*(2-(*diss)(e, n)/std::max((*pdiss)(e, n), omega_min_)) + limiter_ * nu(!cnu * e, n);
            else Process::exit(que_suis_je() + " : ajouter_blocs : probleme !!!") ;

            const double alp_en = (alp ? (*alp)(e, n) : 1.0);
            secmem(e, n) += fac * nut_l;

            if (is_multi)
              {
                const int Nt = equation().probleme().get_champ("temperature").valeurs().line_size();

                const tabs_t&      der_alpha_rho = ref_cast(Champ_Inc_base, ch_alpha_rho).derivees(); // dictionnaire des derivees
                for (auto &&i_m : matrices)
                  {
                    Matrice_Morse& mat = *i_m.second;
                    if (i_m.first == "alpha") 	    mat(N * e + n, Na * e + n) -= 0.*fac * nut_l * alp_en * (der_alpha_rho.count("alpha") ?       der_alpha_rho.at("alpha")(e, n) : 0 );			      // derivee par rapport au taux de vide
                    if (i_m.first == "alpha") 	    mat(N * e + n, Na * e + n) -= fac * nut_l ;			      // derivee par rapport au taux de vide
                    if (i_m.first == "temperature") mat(N * e + n, Nt * e + n) -= 0.*fac * nut_l * alp_en * (der_alpha_rho.count("temperature") ? der_alpha_rho.at("temperature")(e, n) : 0 );// derivee par rapport a la temperature
                    if (i_m.first == "pression")    mat(N * e + n, Np * e + mp)-= 0.*fac * nut_l * alp_en * (der_alpha_rho.count("pression") ?    der_alpha_rho.at("pression")(e, mp) : 0 );		  // derivee par rapport a la pression
                  }
              }

            if (Type_diss == "tau")
              for (auto &&i_m : matrices)
                {
                  Matrice_Morse& mat = *i_m.second;
                  if (i_m.first == "k")         mat(N * e + n,  N * e + n) -= fac * alp_en *(*diss)(e,n);
                  if (i_m.first == "tau")       mat(N * e + n,  N * e + n) -= fac * alp_en * k(e,n);
                }
            else if (Type_diss == "omega")
              for (auto &&i_m : matrices)
                {
                  Matrice_Morse& mat = *i_m.second;
                  if (i_m.first == "k")         mat(N * e + n,  N * e + n) -= fac * alp_en *      1./ std::max((*pdiss)(e, n), omega_min_) + 0. * std::max((*pdiss)(e, n), omega_min_)*(2-(*diss)(e, n)/std::max((*pdiss)(e, n), omega_min_)) ;
                  if (i_m.first == "omega")     mat(N * e + n,  N * e + n) -= fac * alp_en * -k(e,n)/(std::max((*pdiss)(e, n), omega_min_)*std::max((*pdiss)(e, n), omega_min_)) *0.;
                }
          }
    }
}
