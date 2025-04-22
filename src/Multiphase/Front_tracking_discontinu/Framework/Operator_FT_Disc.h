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

#ifndef Operator_FT_Disc_included
#define Operator_FT_Disc_included

#include <TRUSTTabFT.h>
#include <Maillage_FT_Disc.h>

/*! @brief : class Operator_FT_Disc
 */

class Operator_FT_Disc
{
public:
  void Operator_Laplacian_FT_element(const ArrOfDouble& Phi_Facet,const Maillage_FT_Disc& FTmesh, ArrOfDouble& Laplacian_Phi_Facet,DoubleTab& Grad_Phi_Sommet);

  /*void Operator_Gradient_FT_sommets(const ArrOfDouble& Phi_Facet, const Maillage_FT_Disc& FTmesh,
                                    DoubleTab& Grad_Phi_Sommet, bool Normalised_with_Surface);*/

  void Operator_Gradient_FT_sommets(const ArrOfDouble& Phi_Facet, const Maillage_FT_Disc& FTmesh,
                                    DoubleTab& Grad_Phi_Sommet, bool Normalised_with_Surface=true);

  void Compute_interfaciale_source(const ArrOfDouble& sigma_Facet, const Maillage_FT_Disc& FTmesh,
                                   DoubleTab& df_sigma, bool Normalised_with_Surface, bool use_tryggvason_formulation, bool with_marangoni=false);
  void Operator_integral_bord_facette_phi_p_dl(const ArrOfDouble& Surface_sommet, const ArrOfDouble& Phi_sommet, const ArrOfDouble& Phi_Facet, const Maillage_FT_Disc& FTmesh,
                                               DoubleTab& int_phi_p_dl, DoubleTab& int_p_dl);
  void Sommets_to_Facettes(DoubleTab& Phi_Facet, const DoubleTab& Phi_Som, const Maillage_FT_Disc& FTmesh, bool Normalised_with_Surface);
  void Sommets_to_Facettes(ArrOfDouble& Phi_Facet, const ArrOfDouble& Phi_Som, const Maillage_FT_Disc& FTmesh, bool Normalised_with_Surface);
  void Facette_to_Sommets(ArrOfDouble& Surface_sommet, DoubleTab& Phi_Som, const DoubleTab& Phi_Facet, const Maillage_FT_Disc& FTmesh, bool Normalised_with_Surface);
  void Facette_to_Sommets(ArrOfDouble& Surface_sommet, ArrOfDouble& Phi_Som, const ArrOfDouble& Phi_Facet, const Maillage_FT_Disc& FTmesh, bool Normalised_with_Surface);
  //void Operator_Gradient_FT_sommets(const ArrOfDouble& Phi_Facet, DoubleTab& Grad_Phi_Sommet, const Maillage_FT_Disc& FTmesh);
  void produit_vectoriel(const ArrOfDouble& a, const ArrOfDouble& b, ArrOfDouble& resu);
  double norme(const ArrOfDouble& a);
  void unitarisation(ArrOfDouble& a);
  ArrOfDouble Phi_sommet_;
  ArrOfDouble Surface_sommet_;
  DoubleTab Kappa_n_, n_sommet_;
  const ArrOfDouble& get_Phi_sommet() const
  {
    return Phi_sommet_;
  };
  const ArrOfDouble& get_Surface_sommet() const
  {
    return Surface_sommet_;
  };
  const DoubleTab& get_n_sommet() const
  {
    return n_sommet_;
  };

  ArrOfDouble get_n_sommet(int dir) const
  {
    int nbsom = n_sommet_.dimension(0);
    ArrOfDouble  n_dir_sommet;
    n_dir_sommet.resize(nbsom);
    for (int som=0 ; som<nbsom ; som++)
      {
        n_dir_sommet(som)=n_sommet_(som, dir);
      }
    return n_dir_sommet;
  };

};

#endif /* Operator_FT_Disc_included */
