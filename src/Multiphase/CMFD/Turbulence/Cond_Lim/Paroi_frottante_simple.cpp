/****************************************************************************
* Copyright (c) 2023, CEA
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

#include <Paroi_frottante_simple.h>

#include <Op_Dift_Multiphase_VDF_Face.h>
#include <Op_Diff_PolyMAC_P0_base.h>
#include <Loi_paroi_adaptative.h>
#include <Frontiere_dis_base.h>
#include <Navier_Stokes_std.h>
#include <Pb_Multiphase.h>
#include <Domaine_VF.h>
#include <Frontiere.h>
#include <Motcle.h>

#include <math.h>

Implemente_instanciable(Paroi_frottante_simple,"Paroi_frottante_simple", Frottement_global_impose);
// XD Paroi_frottante_simple condlim_base Paroi_frottante_simple 1 Adaptive wall-law boundary condition for velocity

Sortie& Paroi_frottante_simple::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Paroi_frottante_simple::readOn(Entree& s )
{
  if (app_domains.size() == 0)
    app_domains = { Motcle("turbulence") };
  le_champ_front.typer("Champ_front_vide");
  return s;
}

void Paroi_frottante_simple::completer()
{
  if (!ref_cast(Operateur_Diff_base, domaine_Cl_dis().equation().operateur(0).l_op_base()).is_turb())
    Process::exit(que_suis_je() + " : diffusion operator must be turbulent !");
}

void Paroi_frottante_simple::liste_faces_loi_paroi(IntTab& tab)
{
  int nf = la_frontiere_dis->frontiere().nb_faces(), f1 = la_frontiere_dis->frontiere().num_premiere_face();
  int N = tab.line_size();

  for (int f = 0 ; f < nf ; f++)
    for (int n = 0 ; n < N ; n++)
      tab(f + f1, n) |= 1;
}

int Paroi_frottante_simple::initialiser(double temps)
{
  const int nbp = sub_type(Pb_Multiphase, domaine_Cl_dis().equation().probleme()) ? ref_cast(Pb_Multiphase, domaine_Cl_dis().equation().probleme()).nb_phases() : 1;
  valeurs_coeff_.resize(0, nbp);
  la_frontiere_dis->frontiere().creer_tableau_faces(valeurs_coeff_);
  valeurs_coeff_ = 0 ;

  valeurs_coeff_grad_.resize(0, nbp);
  la_frontiere_dis->frontiere().creer_tableau_faces(valeurs_coeff_grad_);
  valeurs_coeff_grad_ = 0 ;

  correlation_loi_paroi_ = domaine_Cl_dis().equation().probleme().get_correlation("Loi_paroi");
  return 1;
}

void Paroi_frottante_simple::mettre_a_jour(double tps)
{
  if (mon_temps < tps)
    {
      me_calculer();
      mon_temps=tps;
    }
}

void Paroi_frottante_simple::me_calculer()
{
  Loi_paroi_adaptative& corr_loi_paroi = ref_cast(Loi_paroi_adaptative, correlation_loi_paroi_.valeur());
  const Domaine_VF& domaine = ref_cast(Domaine_VF, domaine_Cl_dis().equation().probleme().domaine_dis());

  const DoubleTab& u_tau = corr_loi_paroi.get_tab("u_tau"); // y_p est numerote selon les faces du domaine
  const DoubleTab& nu_visc  = ref_cast(Fluide_base, domaine_Cl_dis().equation().probleme().milieu()).viscosite_cinematique().valeurs(),
                   &mu_visc = ref_cast(Fluide_base, domaine_Cl_dis().equation().probleme().milieu()).viscosite_dynamique().valeurs(),
                    &vit   = domaine_Cl_dis().equation().probleme().get_champ("vitesse").passe(),
                     &rho = domaine_Cl_dis().equation().probleme().get_champ("masse_volumique").passe(),
                      *alp = sub_type(Pb_Multiphase, domaine_Cl_dis().equation().probleme()) ? &domaine_Cl_dis().equation().probleme().get_champ("alpha").passe() : nullptr;

  const int cnu = nu_visc.dimension(0) == 1, cmu = mu_visc.dimension(0) == 1, cr = rho.dimension(0) == 1;
  // On va chercher le mu turbulent de polymac et celui de vdf et on prend le bon dans la suite
  const DoubleTab* mu_poly = domaine.que_suis_je().debute_par("Domaine_PolyMAC") ? &ref_cast(Op_Diff_PolyMAC_P0_base, domaine_Cl_dis().equation().operateur(0).l_op_base()).nu() : nullptr,
                   *mu_vdf = domaine.que_suis_je().debute_par("Domaine_VDF") ? &ref_cast(Op_Dift_Multiphase_VDF_Face, domaine_Cl_dis().equation().operateur(0).l_op_base()).get_diffusivite_turbulente() : nullptr;
  assert((mu_poly) || (mu_vdf));

  int nf = la_frontiere_dis->frontiere().nb_faces(), nf_tot = domaine.nb_faces_tot(), f1 = la_frontiere_dis->frontiere().num_premiere_face();
  int N = domaine_Cl_dis().equation().inconnue().valeurs().line_size(), D = dimension;

  const DoubleTab& n_f = domaine.face_normales();
  const DoubleVect& fs = domaine.face_surfaces();
  const IntTab& f_e = domaine.face_voisins();

  DoubleTab pvit_elem(0, N * dimension);
  if (mu_vdf)
    {
      const Champ_Face_base& ch = ref_cast(Champ_Face_base, domaine_Cl_dis().equation().probleme().equation(0).inconnue());
      domaine.domaine().creer_tableau_elements(pvit_elem);
      ch.get_elem_vector_field(pvit_elem, true);
    }

  int n = 0; // The turbulent phase rubs.
  for (int f = 0; f < nf; f++)
    {
      const int f_domaine = f + f1; // number of the face in the domaine
      const int e = f_e(f_domaine,0) >= 0 ? f_e(f_domaine, 0) : f_e(f_domaine, 1);

      double u_orth = 0 ;
      DoubleTrav u_parallel(D);
      if (mu_vdf)
        {
          for (int d = 0; d <D ; d++)
            u_orth -= pvit_elem(e, N*d + n)*n_f(f_domaine,d)/fs(f_domaine); // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -
          for (int d = 0 ; d < D ; d++)
            u_parallel(d) = pvit_elem(e, N*d+n) - u_orth*(-n_f(f_domaine,d))/fs(f_domaine) ; // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -
        }
      else // mu_polymac
        {
          for (int d = 0; d < D ; d++)
            u_orth -= vit(nf_tot + e * D+d, n)*n_f(f_domaine,d)/fs(f_domaine); // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -
          for (int d = 0 ; d < D ; d++)
            u_parallel(d) = vit(nf_tot + e * D + d, n) - u_orth*(-n_f(f_domaine, d))/fs(f_domaine) ; // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -
        }
      const double norm_u_parallel = std::sqrt(domaine.dot(&u_parallel(0), &u_parallel(0)));

      double y_loc = f_e(f_domaine,0) >=0  ? domaine.dist_face_elem0(f_domaine,e) : domaine.dist_face_elem1(f_domaine,e);
      double y_plus_loc = y_loc * u_tau(f_domaine, n)/ nu_visc(!cnu * e, n) ;
      double fac_coeff_grad_ = fac_coeff_grad(y_plus_loc);
      double mu_tot_loc = (mu_poly) ? ((alp ? 1.0 : rho(!cr * e, n)) * (*mu_poly)(e,n)) : (mu_vdf) ? (*mu_vdf)(e,n) + mu_visc(!cmu * e,n) : -1;
      if (y_plus_loc>1)
        {
          valeurs_coeff_(f, n) = (alp ? (*alp)(e, n) * rho(!cr * e, n) : 1) * u_tau(f_domaine, n)*u_tau(f_domaine, n)/norm_u_parallel; // f_tau = - alpha_k rho_k u_tau**2 n_par, coeff = u_tau**2 /u_par
          valeurs_coeff_grad_(f, n) =  fac_coeff_grad_*(alp ? (*alp)(e, n) : 1) * std::min(1./y_loc, 1/mu_tot_loc * rho(!cr * e, n) * u_tau(f_domaine, n)*u_tau(f_domaine, n)/norm_u_parallel); // f_tau = - alpha_k rho_k u_tau**2 n_par, coeff = u_tau**2 /u_par
        }
      else
        {
          valeurs_coeff_(f, n) = (alp ? (*alp)(e, n) * rho(!cr * e, n) : 1) * nu_visc(!cnu * e, n)/y_loc; // viscous case : if u_tau is small
          valeurs_coeff_grad_(f, n) = fac_coeff_grad_*(alp ? (*alp)(e, n) : 1) * 1./y_loc ; // dirichlet for calculation of gradient
        }
    }

  // Treatment of other phases which are supposed not to be turbulent
  for (n = 1; n < N; n++)
    for (int f = 0; f < nf; f++)
      {
        valeurs_coeff_(f, n) = 0; // les phases non turbulentes sont non porteuses : pas de contact paroi => des symmetries
        valeurs_coeff_grad_(f, n) = 0; // les phases non turbulentes sont non porteuses : pas de contact paroi => des symmetries
      }

  valeurs_coeff_.echange_espace_virtuel();
  valeurs_coeff_grad_.echange_espace_virtuel();
}
