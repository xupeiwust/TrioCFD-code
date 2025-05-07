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

#ifndef Flux_2groupes_base_included
#define Flux_2groupes_base_included

#include <Correlation_base.h>
#include <TRUSTTabs.h>
#include <Saturation_base.h>

/*! @brief classe Flux_2groupes_base correlations de flux entre 2 groupes
 *
 *
 */
class Flux_2groupes_base : public Correlation_base
{
  Declare_base(Flux_2groupes_base);
public:
  /* parametres d'entree coeffs */
  struct input_coeffs
  {
    double dh ;            // diametre hyd
    DoubleTab alpha ;  // alpha[n] : taux de vide de la phase n
    DoubleTab p ;             // pression
    DoubleTab nv ;     // nv : norme de ||v_k - v_l||
    DoubleTab mu ;     // mu[n]         : viscosite dynamique de la phase n
    DoubleTab rho ;    // rho[n]        : masse volumique de la phase n
    DoubleTab k_turb ; // k_turb   : energie cinetique turbulente de la phase n
    DoubleTab epsilon ;    // nut[n]        : viscosite turbulente de bulles de la phase n
    DoubleTab sigma ;  //sigma[ind_trav]:tension superficielle sigma(ind_trav), ind_trav = (n*(N-1)-(n-1)*(n)/2) + (m-n-1)
    DoubleTab d_bulles ;//d_bulles[n]   : diametre de bulles de la phase n
    DoubleTab a_i ;
    int e;                // indice d'element
    int n_l ;
    int n_g1 ;
    int n_g2 ;
  };
  /* valeurs de sortie coeffs */
  struct output_coeffs
  {
    DoubleTab gamma ;
    DoubleTab da_gamma;
    DoubleTab dai_gamma;
    double inter2g1 ;
    double inter3g1 ;
    double inter2g2 ;
    double inter3g2 ;
    double da_inter2g1 ;
    double da_inter3g1 ;
    double da_inter2g2 ;
    double da_inter3g2 ;

  };
  /* parametres d'entree therms */
  struct input_therms
  {
    int n_l;
    int n_g1;
    int n_g2;
    double dh;
    DoubleTab alpha;
    DoubleTab T;
    DoubleTab p;
    DoubleTab d_bulles;
    DoubleTab lambda;
    DoubleTab mu;
    DoubleTab rho;
    DoubleTab Cp;
    DoubleTab sigma ;
    double Lvap;
    double qp ;
    double hl;
    double hlsat ;
    double Tsatg1 ;
    double Tsatg2 ;
    double dp_Tsat1 ;
    double dp_Tsat2 ;
  };
  /* valeurs de sortie therms */
  struct output_therms
  {
    DoubleTab dT_G1;
    DoubleTab dT_G2;
    DoubleTab da_G1;
    DoubleTab da_G2;
    double dp_G1;
    double dp_G2;
    double G1;
    double G2;
    DoubleTab dT_etaph1;
    DoubleTab dT_etaph2;
    DoubleTab da_etaph1;
    DoubleTab da_etaph2;
    double dp_etaph1;
    double dp_etaph2;
    double etaph1;
    double etaph2;

  };

  virtual void coeffs(const input_coeffs& input, output_coeffs& output) const = 0;
  virtual void therm(const input_therms& input, output_therms& output) const = 0;

};

#endif
