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
// File:        Dispersion_bulles_turbulente_Lavieville.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

// NeptuneCFD like force no publication

#include <Dispersion_bulles_turbulente_Lavieville.h>
#include <Pb_Multiphase.h>
#include <QDM_Multiphase.h>
#include <TRUSTTrav.h>
#include <Frottement_interfacial_base.h>
#include <Masse_ajoutee_base.h>
#include <math.h>

Implemente_instanciable(Dispersion_bulles_turbulente_Lavieville, "Dispersion_bulles_turbulente_Lavieville", Dispersion_bulles_base);

Sortie& Dispersion_bulles_turbulente_Lavieville::printOn(Sortie& os) const
{
  return os;
}

Entree& Dispersion_bulles_turbulente_Lavieville::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.lire_avec_accolades_depuis(is);

  const Pb_Multiphase *pbm = sub_type(Pb_Multiphase, pb_.valeur()) ? &ref_cast(Pb_Multiphase, pb_.valeur()) : nullptr;

  if (!pbm || pbm->nb_phases() == 1) Process::exit(que_suis_je() + " : not needed for single-phase flow!");
  for (int n = 0; n < pbm->nb_phases(); n++) //recherche de n_l, n_g : phase {liquide,gaz}_continu en priorite
    if (pbm->nom_phase(n).debute_par("liquide") && (n_l < 0 || pbm->nom_phase(n).finit_par("continu")))  n_l = n;
  if (n_l < 0) Process::exit(que_suis_je() + " : liquid phase not found!");

  if (pbm->has_correlation("frottement_interfacial")) correlation_drag_ = pbm->get_correlation("frottement_interfacial"); //correlation fournie par le bloc correlation
  else Correlation_base::typer_lire_correlation(correlation_drag_, *pbm, "frottement_interfacial", is);
  if (pbm->has_correlation("masse_ajoutee")) correlation_masse_ajoutee_ = pbm->get_correlation("masse_ajoutee"); //correlation fournie par le bloc correlation
  else Correlation_base::typer_lire_correlation(correlation_masse_ajoutee_, *pbm, "masse_ajoutee", is);

  return is;
}


void Dispersion_bulles_turbulente_Lavieville::coefficient(const input_t& in, output_t& out) const
{
  const Frottement_interfacial_base& corr = ref_cast(Frottement_interfacial_base, correlation_drag_.valeur());
  const Masse_ajoutee_base& corrCam = ref_cast(Masse_ajoutee_base, correlation_masse_ajoutee_.valeur());
  int N = out.Ctd.dimension(0);

  DoubleTab coeff_drag(N, N);
  DoubleTab coeff_Cam(N); // a zero ?

  corrCam.coeff(in.alpha, in.rho, coeff_Cam) ; //

  corr.coefficient_CD( in.alpha, in.p, in.T, in.rho, in.mu, in.sigma, in.dh, in.nv, in.d_bulles, coeff_drag);

  out.Ctd = 0;

  for (int k = 0; k < N; k++)
    if (k!=n_l)
      {

        double gamma = 0.5 * in.rho[n_l] / in.rho[k] ;
        double u_r = in.nv(k,n_l) ;
        double tau_t = 3. / 2. * in.nut[n_l] / in.k_turb[n_l] * std::pow (1. + 2.07 * u_r * u_r / in.k_turb[n_l], -0.5 ) ;
        double tau_f = ( 1. / (2. * gamma) + coeff_Cam(k) ) * 4./3. * in.d_bulles[k] / ( coeff_drag(k, n_l) * std::max(u_r, 1e-5) ) ;

        double b = (1. + coeff_Cam(k) )/(1. / (2.*gamma) + coeff_Cam(k) ) ;
        double eta_r = tau_t / tau_f ;



        // Calcul des coefficients de dispersion turbulente
        out.Ctd(k, n_l) = coeff_Cam(k) * ( b * (b-1.) / (1.+eta_r) + eta_r * ( (b + eta_r) / (1. + eta_r) + in.alpha[k] / std::max(1.-in.alpha[k],1e-5) ) ) + 1. / (2. * gamma) * ((b * b + eta_r) / (1.+eta_r) + eta_r * ( (b + eta_r) / (1. + eta_r) + in.alpha[k] / std::max(1.-in.alpha[k],1e-3))) - ( (b + eta_r) / (1. + eta_r) + in.alpha[k] / std::max(1.-in.alpha[k],1e-5) ) * in.rho[n_l] * 2./3. * in.k_turb[n_l];
      }
}
