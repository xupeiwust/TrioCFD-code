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
// File:        Flux_parietal_Hibiki.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

//Caleb S. Brooks, Takashi Hibiki, Wall nucleation modeling in subcooled boiling flow, International Journal of Heat and Mass Transfer, 2015, https://doi.org/10.1016/j.ijheatmasstransfer.2015.03.005

#include <Flux_parietal_Hibiki.h>
#include <Flux_parietal_adaptatif.h>
#include <Loi_paroi_adaptative.h>
#include <Correlation_base.h>
#include <Pb_Multiphase.h>
#include <Domaine_dis_base.h>
#include <Domaine_VF.h>
#include <TRUSTTrav.h>
#include <Milieu_composite.h>
#include <Saturation_base.h>

#include <math.h>

Implemente_instanciable(Flux_parietal_Hibiki, "Flux_parietal_Hibiki", Flux_parietal_base);

Sortie& Flux_parietal_Hibiki::printOn(Sortie& os) const { return Flux_parietal_base::printOn(os); }

Entree& Flux_parietal_Hibiki::readOn(Entree& is)
{
  const Pb_Multiphase& pbm = ref_cast(Pb_Multiphase, pb_.valeur());
  Correlation_base::typer_lire_correlation(correlation_monophasique_, pbm, "Flux_parietal", is);
  Cout << que_suis_je() << " : single-phase wall heat flux is " << correlation_monophasique_->que_suis_je() << finl;

  Param param(que_suis_je());
  param.ajouter("contact_angle_deg",&theta_);
  param.ajouter("molar_mass",&molar_mass_,Param::REQUIRED);
  param.ajouter("Qw",&Qw_,Param::REQUIRED);
  param.ajouter("G",&G_,Param::REQUIRED);
  param.lire_avec_accolades(is);

  if ( !sub_type(Milieu_composite, pb_->milieu())) Process::exit("Flux_parietal_Hibiki::readOn : the medium must be composite !");
  if (!pbm.nom_phase(0).debute_par("liquide")) Process::exit("Flux_parietal_Hibiki::readOn : the first phase must be liquid !");

  for (int n = 0; n < pbm.nb_phases(); n++)  //recherche de n_l, n_g : phase {liquide,gaz}_continu en priorite
    {
      if (pbm.nom_phase(n).debute_par("liquide") && (n_l < 0 || pbm.nom_phase(n).finit_par("continu")))  n_l = n;
      if (( pbm.nom_phase(n).finit_par("group1")))  n_g1 = n;
      if (( pbm.nom_phase(n).finit_par("group2")))  n_g2 = n;
    }
  if (n_l < 0) Process::exit(que_suis_je() + " : liquid phase not found!");
  if (n_g1 < 0) Process::exit(que_suis_je() + " : group 1 not found!");
  if (n_g2 < 0) Process::exit(que_suis_je() + " : group 2 not found!");

  return is;
}

void Flux_parietal_Hibiki::completer()
{
  correlation_monophasique_->completer();
}

void Flux_parietal_Hibiki::qp(const input_t& in, output_t& out) const
{
  // On met tout a 0 a tout hasard
  if (out.qpk)     (*out.qpk)    = 0.;
  if (out.da_qpk)  (*out.da_qpk) = 0.;
  if (out.dp_qpk)  (*out.dp_qpk) = 0.;
  if (out.dv_qpk)  (*out.dv_qpk) = 0.;
  if (out.dTf_qpk) (*out.dTf_qpk)= 0.;
  if (out.dTp_qpk) (*out.dTp_qpk)= 0.;
  if (out.qpi)     (*out.qpi)    = 0.;
  if (out.da_qpi)  (*out.da_qpi) = 0.;
  if (out.dp_qpi)  (*out.dp_qpi) = 0.;
  if (out.dv_qpi)  (*out.dv_qpi) = 0.;
  if (out.dTf_qpi) (*out.dTf_qpi)= 0.;
  if (out.dTp_qpi) (*out.dTp_qpi)= 0.;
  if (out.nonlinear) (*out.nonlinear) = 1;

  // On remplit le monophasique ; pas besoin du flux interfacial normalement
  ref_cast(Flux_parietal_base, correlation_monophasique_.valeur()).qp(in, out);

  // Ici la phase liquide est forcement la phase 0 car la correlation monophasique ne remplit que la phase 0
  const Milieu_composite& milc = ref_cast(Milieu_composite, pb_->milieu());

  if (milc.has_saturation(n_l, n_g1))
    {
      int ind_sat = n_g1<n_l ? ( n_g1 *(in.N-1)-( n_g1 -1)*( n_g1 )/2) + (n_l- n_g1 -1) :
                    (n_l*(in.N-1)-(n_l-1)*(n_l)/2) + ( n_g1 -n_l-1);

      double Delta_T_sup = in.Tp - in.Tsat[ind_sat]; // Wall superheat

      if (Delta_T_sup > 0) // Else : no wall superheat => no nucleation => single phase heat transfer only
        {

          double JaT = in.Cp[n_l] * std::max(in.Tp - in.T[n_l], 0.) / in.Lvap[ind_sat] ;
          double dTp_JaT = in.Tp - in.T[n_l] > 0. ? in.Cp[n_l] * std::max(in.Tp , 1.e-8) / in.Lvap[ind_sat] : 0. ;
          double dTl_JaT = in.Tp - in.T[n_l] > 0. ? - in.Cp[n_l] * in.T[n_l] / in.Lvap[ind_sat] : 0. ;
          double Jaw = in.Cp[n_l] * std::max(in.Tp - in.Tsat[ind_sat], 1.e-8) / in.Lvap[ind_sat] ;
          double dTp_Jaw = in.Tp - in.Tsat[ind_sat] > 0. ? in.Cp[n_l] * in.Tp / in.Lvap[ind_sat] : 0. ;
          double Prl = (in.mu[n_l]*in.Cp[n_l])/in.lambda[n_l];
          double Bo = Qw_ / G_ / in.Lvap[ind_sat] ;


          // Nucleation site density (Hibiki Ishii 2003)
          double N_sites, dTp_N_sites, dTl_N_sites, dTg_N_sites;
          N_sites     =     Hibiki_Ishii_Site_density(in.rho[n_g1], in.rho[n_l], in.T[n_g1], in.T[n_l], in.p, in.Tp, in.Lvap[ind_sat], in.Tsat[ind_sat], in.Sigma[ind_sat], theta_, molar_mass_);
          dTp_N_sites = dTp_Hibiki_Ishii_Site_density(in.rho[n_g1], in.rho[n_l], in.T[n_g1], in.T[n_l], in.p, in.Tp, in.Lvap[ind_sat], in.Tsat[ind_sat], in.Sigma[ind_sat], theta_, molar_mass_);
          dTl_N_sites = dTl_Hibiki_Ishii_Site_density(in.rho[n_g1], in.rho[n_l], in.T[n_g1], in.T[n_l], in.p, in.Tp, in.Lvap[ind_sat], in.Tsat[ind_sat], in.Sigma[ind_sat], theta_, molar_mass_);
          dTg_N_sites = dTg_Hibiki_Ishii_Site_density(in.rho[n_g1], in.rho[n_l], in.T[n_g1], in.T[n_l], in.p, in.Tp, in.Lvap[ind_sat], in.Tsat[ind_sat], in.Sigma[ind_sat], theta_, molar_mass_);

          // Departure diameter
          double D_d     = in.Tp - in.T[n_l] > 0.? 2.11e-3 * std::pow(JaT,-0.49) * std::pow(in.rho[n_g1]/in.rho[n_l],-0.78) * std::pow(Bo,0.44) * std::pow(Prl,1.72) : 1e-8 ;
          double dTp_D_d = in.Tp - in.T[n_l] > 0. ?  2.11e-3 * -0.49 * dTp_JaT * std::pow(JaT,-1.49) * std::pow(in.rho[n_g1]/in.rho[n_l],-0.78) * std::pow(Bo,0.44) * std::pow(Prl,1.72) :  0.;
          double dTl_D_d = in.Tp - in.T[n_l] > 0. ? 2.11e-3 * -0.49 * dTl_JaT *std::pow(JaT,-1.49) * std::pow(in.rho[n_g1]/in.rho[n_l],-0.78) * std::pow(Bo,0.44) * std::pow(Prl,1.72) : 0.;

          // Bubble departure frequency
          double f_dep = in.Tp - in.T[n_l] > 0.? D_d * D_d / (in.lambda[n_l]/(in.rho[n_l]*in.Cp[n_l])) * std::pow(Jaw,2.28) * std::pow(in.rho[n_g1]/in.rho[n_l],-0.93) * std::pow(JaT,-1.46) * std::pow(Prl,2.36) : 0. ;
          double dTp_f_dep = in.Tp - in.T[n_l] > 0. ?  2. * D_d * dTp_D_d   / (in.lambda[n_l]/(in.rho[n_l]*in.Cp[n_l])) * std::pow(Jaw,2.28) * std::pow(in.rho[n_g1]/in.rho[n_l],-0.93) * std::pow(JaT,-1.46) * std::pow(Prl,2.36)
                             + D_d * D_d / (in.lambda[n_l]/(in.rho[n_l]*in.Cp[n_l])) * 2.28 * dTp_Jaw * std::pow(Jaw,1.28) * std::pow(in.rho[n_g1]/in.rho[n_l],-0.93) * std::pow(JaT,-1.46) * std::pow(Prl,2.36)
                             +  D_d * D_d / (in.lambda[n_l]/(in.rho[n_l]*in.Cp[n_l])) * std::pow(Jaw,2.28) * std::pow(in.rho[n_g1]/in.rho[n_l],-0.93) * -1.46 * dTp_JaT * std::pow(JaT,-2.46) * std::pow(Prl,2.36) : 0. ;
          double dTl_f_dep = in.Tp - in.T[n_l] > 0. ? 2. * D_d * dTl_D_d  / (in.lambda[n_l]/(in.rho[n_l]*in.Cp[n_l])) * std::pow(Jaw,2.28) * std::pow(in.rho[n_g1]/in.rho[n_l],-0.93) * std::pow(JaT,-1.46) * std::pow(Prl,2.36)
                             +  D_d * D_d / (in.lambda[n_l]/(in.rho[n_l]*in.Cp[n_l])) * std::pow(Jaw,2.28) * std::pow(in.rho[n_g1]/in.rho[n_l],-0.93) * -1.46 *dTl_JaT * std::pow(JaT,-2.46) * std::pow(Prl,2.36) : 0. ;


          // Evaporation
          if (out.qpi)         (*out.qpi)(n_l, n_g1)      =  1./6.*M_PI * in.rho[n_g1] * in.Lvap[ind_sat] * std::pow(D_d,3.) * f_dep * N_sites;
          if (out.dTp_qpi) (*out.dTp_qpi)(n_l, n_g1)      =1./6.*M_PI * in.rho[n_g1] * in.Lvap[ind_sat] * (3.*dTp_D_d * std::pow(D_d,2.) * f_dep * N_sites
                                                                                                             + std::pow(D_d,3.) * dTp_f_dep *     N_sites
                                                                                                             + std::pow(D_d,3.) * f_dep * dTp_N_sites);
          if (out.dTf_qpi) (*out.dTf_qpi)(n_l, n_g1, n_l) = 1./6.*M_PI * in.rho[n_g1] * in.Lvap[ind_sat] * (3.*dTl_D_d * std::pow(D_d,2.) *     f_dep *     N_sites+  std::pow(D_d,3.) * dTl_f_dep * N_sites +  std::pow(D_d,3.) * f_dep * dTl_N_sites);
          if (out.dTf_qpi) (*out.dTf_qpi)(n_l, n_g1,   n_g1) = 1./6.*M_PI * in.rho[n_g1] * in.Lvap[ind_sat] * ( std::pow(D_d,3.) *     f_dep * dTg_N_sites);

          if (out.d_nuc) (*out.d_nuc)(n_g1) = D_d;

        }
    }
}

// See Hibiki, Ishii 2003
double Flux_parietal_Hibiki::Hibiki_Ishii_Site_density(double rho_v, double rho_l, double T_v, double T_l, double p, double Tp, double L_vap, double T_sat, double sigma, double theta, double molar_mass) const
{
  double N_bar = 4.72e5;
  double mu2 = 0.722*180/M_PI;
  mu2 *= mu2;

  double rho_p = std::log( (rho_l-rho_v)/rho_v);
  double f_rho_plus = -0.01064 + 0.48246*rho_p - 0.22712*rho_p*rho_p + 0.05468*rho_p*rho_p*rho_p;

  double lambda_prime = 2.5e-6;
  double R = 8.314462618 / molar_mass ;

  double inv_Rc = p * (-1 + std::exp(L_vap*std::max(T_v-T_sat, 0.)/( R * T_v * T_sat ))) / (2 * sigma * (1+rho_v/rho_l));

  return N_bar * (1 - std::exp(-theta*theta/(8.*mu2))) * (-1. + std::exp(f_rho_plus * lambda_prime * inv_Rc));
}

double Flux_parietal_Hibiki::dTg_Hibiki_Ishii_Site_density(double rho_v, double rho_l, double T_v, double T_l, double p, double Tp, double L_vap, double T_sat, double sigma, double theta, double molar_mass) const
{
  double N_bar = 4.72e5;
  double mu2 = 0.722*180/M_PI;
  mu2 *= mu2;

  double rho_p = std::log( (rho_l-rho_v)/rho_v);
  double f_rho_plus = -0.01064 + 0.48246*rho_p - 0.22712*rho_p*rho_p + 0.05468*rho_p*rho_p*rho_p;

  double lambda_prime = 2.5e-6;
  double R = 8.314462618 / molar_mass ;

  double inv_Rc = p * (-1 + std::exp(L_vap*std::max(T_v-T_sat, 0.)/( R * T_v * T_sat ))) / (2. * sigma * (1+rho_v/rho_l));
  double dTg_inv_Rc = 0.;
  if (T_v-T_sat> 0.)
    dTg_inv_Rc = p * L_vap / ( R * T_v*T_v ) * std::exp(L_vap*(T_v-T_sat)/(R*T_v*T_sat)) / (2. * sigma * (1+rho_v/rho_l));

  return N_bar * (1 - std::exp(-theta*theta/(8.*mu2))) * (0. + f_rho_plus * lambda_prime * dTg_inv_Rc* std::exp(f_rho_plus * lambda_prime * inv_Rc));
}

double Flux_parietal_Hibiki::dTl_Hibiki_Ishii_Site_density(double rho_v, double rho_l, double T_v, double T_l, double p, double Tp, double L_vap, double T_sat, double sigma, double theta, double molar_mass) const
{
  return 0.;
}


double Flux_parietal_Hibiki::dTp_Hibiki_Ishii_Site_density(double rho_v, double rho_l, double T_v, double T_l, double p, double Tp, double L_vap, double T_sat, double sigma, double theta, double molar_mass) const
{
  return 0.;
}
