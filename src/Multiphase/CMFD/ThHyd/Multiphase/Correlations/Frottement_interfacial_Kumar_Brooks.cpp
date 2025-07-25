/****************************************************************************
* Copyright (c) 2024, CEA
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


/// Vineet Kumar, Caleb S. Brooks, Inter-group mass transfer modeling in the two-group two-fluid model with interfacial area transport equation in condensing flow, International Journal of Heat and Mass Transfer, 2018, https://doi.org/10.1016/j.ijheatmasstransfer.2017.11.087


////////////////////////////////////////////////////////////////////////////




#include <Frottement_interfacial_Kumar_Brooks.h>
#include <Pb_Multiphase.h>
#include <Milieu_composite.h>
#include <math.h>

Implemente_instanciable(Frottement_interfacial_Kumar_Brooks, "Frottement_interfacial_Kumar_Brooks", Frottement_interfacial_base);

Sortie& Frottement_interfacial_Kumar_Brooks::printOn(Sortie& os) const
{
  return os;
}

Entree& Frottement_interfacial_Kumar_Brooks::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("Dh", &Dh_);
  param.ajouter("constante_gravitation", &g_);
  param.lire_avec_accolades_depuis(is);

  const Pb_Multiphase *pbm = sub_type(Pb_Multiphase, pb_.valeur()) ? &ref_cast(Pb_Multiphase, pb_.valeur()) : nullptr;

  if (!pbm || pbm->nb_phases() == 1) Process::exit(que_suis_je() + " : not needed for single-phase flow!");
  for (int n = 0; n < pbm->nb_phases(); n++) //recherche de n_l, n_g : phase {liquide,gaz}_continu en priorite
    if (pbm->nom_phase(n).debute_par("liquide") && (n_l < 0 || pbm->nom_phase(n).finit_par("continu")))  n_l = n;
  if (n_l < 0) Process::exit(que_suis_je() + " : liquid phase not found!");



  return is;
}

void Frottement_interfacial_Kumar_Brooks::completer()
{

}


void Frottement_interfacial_Kumar_Brooks::coefficient(const DoubleTab& alpha, const DoubleTab& p, const DoubleTab& T,
                                                      const DoubleTab& rho, const DoubleTab& mu, const DoubleTab& sigma, double Dh,
                                                      const DoubleTab& ndv, const DoubleTab& d_bulles, DoubleTab& coeff) const
{
  int N = ndv.dimension(0);

  coeff = 0;

  for (int k = 0; k < N; k++)
    if (k!=n_l)
      {
        int ind_trav = (k>n_l) ? (n_l*(N-1)-(n_l-1)*(n_l)/2) + (k-n_l-1) : (k*(N-1)-(k-1)*(k)/2) + (n_l-k-1);

        double Eo = g_ * std::abs(rho(n_l)-rho(k)) * d_bulles(k)*d_bulles(k)/sigma(ind_trav);
        double Cd = (std::sqrt(Eo)<4.)? 2./3.* std::sqrt(Eo) / std::sqrt(1.-alpha(k)) : 4./3. / 0.35 / 0.35 * std::sqrt(g_ * std::abs(rho(n_l)-rho(k))/rho(n_l)/Dh_) * d_bulles(k) /(1.-alpha(k));

        coeff(k, n_l, 1) = (alpha(n_l) < 1.e-6) ? 3./4.*Cd/d_bulles(k) * alpha(k) * rho(n_l) * alpha(n_l) * 1.e6
                           : 3./4.* Cd / d_bulles(k) * alpha(k) * rho(n_l);
        coeff(k, n_l, 0) = coeff(k, n_l, 1) * ndv(n_l,k);
        coeff(n_l, k, 0) = coeff(k, n_l, 0);
        coeff(n_l, k, 1) = coeff(k, n_l, 1);

      }
}


void Frottement_interfacial_Kumar_Brooks::coefficient_CD(const DoubleTab& alpha, const DoubleTab& p, const DoubleTab& T,
                                                         const DoubleTab& rho, const DoubleTab& mu, const DoubleTab& sigma, double Dh,
                                                         const DoubleTab& ndv, const DoubleTab& d_bulles, DoubleTab& coeff) const
{
  int N = ndv.dimension(0);

  coeff = 0;

  for (int k = 0; k < N; k++)
    if (k!=n_l)
      {
        int ind_trav = (k>n_l) ? (n_l*(N-1)-(n_l-1)*(n_l)/2) + (k-n_l-1) : (k*(N-1)-(k-1)*(k)/2) + (n_l-k-1);


        double Eo = g_ * std::abs(rho(n_l)-rho(k)) * d_bulles(k)*d_bulles(k)/sigma(ind_trav);
        double Cd = (std::sqrt(Eo)<4.)? 2./3.*std::sqrt(Eo) / std::sqrt(1.-alpha(k)) : 4./3. / 0.35 / 0.35 * std::sqrt(g_ * std::abs(rho(n_l)-rho(k))/rho(n_l)/Dh_) * d_bulles(k) / (1.-alpha(k));

        coeff(k, n_l) = (coeff(n_l, k) = Cd);
      }
}
