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
// File:        Injection_QDM_nulle_PolyMAC_P0.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Injection_QDM_nulle_VDF.h>
#include <Source_Flux_interfacial_base.h>
#include <Masse_ajoutee_base.h>
#include <Milieu_composite.h>
#include <Neumann_val_ext.h>
#include <Saturation_base.h>
#include <Champ_Face_VDF.h>
#include <Neumann_paroi.h>
#include <Pb_Multiphase.h>
#include <Champ_P0_VDF.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <Domaine_VF.h>
#include <Conds_lim.h>
#include <cfloat>
#include <math.h>

Implemente_instanciable(Injection_QDM_nulle_VDF, "Injection_QDM_nulle_VDF_Face", Source_injection_QDM_base);
// XD Injection_QDM_nulle source_base Injection_QDM_nulle 1 not_set


Sortie& Injection_QDM_nulle_VDF::printOn(Sortie& os) const {  return os; }

Entree& Injection_QDM_nulle_VDF::readOn(Entree& is) { return Source_injection_QDM_base::readOn(is);}

void Injection_QDM_nulle_VDF::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Champ_Face_VDF& ch = ref_cast(Champ_Face_VDF, equation().inconnue());
  const Champ_P0_VDF& cha= ref_cast(Champ_P0_VDF, equation().probleme().equation(1).inconnue()); // volume fraction
  const Domaine_VF&  domaine = ref_cast(Domaine_VF, equation().domaine_dis());
  const Conds_lim&      clsa = cha.domaine_Cl_dis().les_conditions_limites();
  const Milieu_composite& milc = ref_cast(Milieu_composite, equation().milieu());

  const IntTab&  fcl = ch.fcl(),
                 &fcla = cha.fcl(),
                  &f_e = domaine.face_voisins(),
                   &e_f = domaine.elem_faces();
  const DoubleVect& vf = domaine.volumes_entrelaces(),
                    &fs = domaine.face_surfaces();
  const DoubleTab& vf_dir = domaine.volumes_entrelaces_dir();

  const DoubleTab& vit = ch.valeurs(),
                   &rho   = equation().milieu().masse_volumique().passe(), // passe car qdm
                    &alpha = cha.passe();

  Matrice_Morse *mat = matrices.count(ch.le_nom().getString()) ? matrices.at(ch.le_nom().getString()) : nullptr; // Derivee locale/QDM

  const Pb_Multiphase *pbm = sub_type(Pb_Multiphase, equation().probleme()) ? &ref_cast(Pb_Multiphase, equation().probleme()) : nullptr;
  const Masse_ajoutee_base *corr = pbm && pbm->has_correlation("masse_ajoutee") ? &ref_cast(Masse_ajoutee_base, pbm->get_correlation("masse_ajoutee")) : nullptr;

  int N = vit.line_size(), nf_tot = domaine.nb_faces_tot(), nf = domaine.nb_faces(), ne_tot = domaine.nb_elem_tot();

  // Cas adiabatique : injection de bulles a la paroi (manip de Gabillet et al.)
  for (int f = 0 ; f< nf_tot ; f ++)
    if (fcla(f, 0) == 4) // Neumann_paroi
      {
        int e = f_e(f, 0)>-1 ? f_e(f, 0) : f_e(f, 1) ;
        if (e>=0)
          {
            DoubleTrav f_a_masse(N, N) ;
            DoubleTrav f_a(1, N);
            for (int n = 0 ; n<N ; n++)
              {
                f_a_masse(n, n) = ref_cast(Neumann_paroi, clsa[fcla(f, 1)].valeur()).flux_impose(fcla(f, 2), n) * rho(e, n) ;   // Pas de porosite ; unite : kg/m3 m/s
                f_a(0, n)       = ref_cast(Neumann_paroi, clsa[fcla(f, 1)].valeur()).flux_impose(fcla(f, 2), n) ;               // Pas de porosite ; unite : m/s
              }
            corr->ajouter_inj( &f_a(0,0) , &alpha(e, 0),   &rho(e, 0)    ,   f_a_masse   );

            for (int i=0 ; i < e_f.line_size() ; i++)
              {
                int f2 = e_f(e, i);
                if ( (f2 >= 0) && (f2<nf) && (fcl(f2, 0) < 2)  ) // Si pas face de bord
                  {
                    int c = ( e == f_e(f2, 0) ) ? 0 : 1 ;
                    for (int n = 0 ; n<N ; n++)
                      for (int m = 0 ; m<N ; m++)
                        {
                          secmem(f2, n) -= fs(f) * f_a_masse(n, m) * vf_dir(f2, c) / vf(f2)  * vit( f2, m) * beta_;
                          if (mat)
                            (*mat)( N * f2 + n , N * f2 + m ) += fs(f) * f_a_masse(n, m) * vf_dir(f2, c) / vf(f2) * beta_;
                        }
                  }
              }
          }
      }

  // Cas bouillant : on passe par qpi
  if ( (pbm) && pbm->has_correlation("flux_parietal") )
    {
      const DoubleTab& qpi = ref_cast(Source_Flux_interfacial_base, pbm->equation_energie().sources().dernier().valeur()).qpi(),
                       &press = ref_cast(QDM_Multiphase, pbm->equation_qdm()).pression().passe();
      int nb_max_sat = N*(N-1)/2;

      /* limiteur de changement de phase : pas mis dans la V0 de cette force */

      // On a besoin juste de l'enthalpie de changement d'etat
      DoubleTrav Lvap_tab(ne_tot,nb_max_sat);
      DoubleTrav f_a_masse(N, N) ;
      DoubleTrav f_a(1, N);

      for (int k = 0; k < N; k++)
        for (int l = k + 1; l < N; l++)
          if (milc.has_saturation(k, l))
            {
              Saturation_base& z_sat = milc.get_saturation(k, l);
              const int ind_trav = (k*(N-1)-(k-1)*(k)/2) + (l-k-1); // Et oui ! matrice triang sup !
              assert (press.line_size() == 1);

              z_sat.Lvap(press.get_span_tot(), Lvap_tab.get_span_tot(), nb_max_sat, ind_trav);
            }

      for (int f = 0 ; f< nf_tot ; f ++)
        if ( (f_e(f, 0)<= 0) || (f_e(f, 1)<= 0)) // Si face de bord
          if ( fs(f) > DBL_MIN ) // Si pas dans le milieu d'un bidim_axi
            {
              int e = f_e(f, 0)>-1 ? f_e(f, 0) : f_e(f, 1) ;
              f_a_masse = 0;
              f_a = 0;
              double G=0;

              for (int k = 0; k < N; k++)
                for (int l = k + 1; l < N; l++)
                  if (milc.has_saturation(k, l)) //flux phase k <-> phase l si saturation
                    {
                      const int i_sat = (k*(N-1)-(k-1)*(k)/2) + (l-k-1); // Et oui ! matrice triang sup !

                      G = qpi(e, k, l) / Lvap_tab(e, i_sat) ; // Flux de matiere de la phase k vers l ; attention : qpi est en W, donc G en kgs-1 : il n'est pas volumique !

                      f_a(0, k)       -= G/(fs(f)*rho(e, k)) ;
                      f_a_masse(k, k) -= G/fs(f) ;
                      f_a(0, l)       += G/(fs(f)*rho(e, l)) ;
                      f_a_masse(l, l) += G/fs(f) ;

                    }

              corr->ajouter_inj( &f_a(0,0) , &alpha(e, 0),   &rho(e, 0)    ,   f_a_masse   );

              for (int i=0 ; i < e_f.line_size() ; i++)
                {
                  int f2 = e_f(e, i);
                  if ( (f2 >= 0) && (f2<nf) && (fcl(f2, 0) < 2)  ) // Si pas face de bord
                    {
                      int c = ( e == f_e(f2, 0) ) ? 0 : 1 ;
                      for (int n = 0 ; n<N ; n++)
                        for (int m = 0 ; m<N ; m++)
                          {
                            secmem(f2, n) -= fs(f) * f_a_masse(n, m) * vf_dir(f2, c) / vf(f2)  * vit( f2, m) * beta_ ;
                            if (mat)
                              (*mat)( N * f2 + n , N * f2 + m ) += fs(f) * f_a_masse(n, m) * vf_dir(f2, c) / vf(f2) * beta_ ;
                          }
                    }
                }
            }
    }
}
