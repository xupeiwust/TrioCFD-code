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


#include <Operator_FT_Disc.h>
using namespace std;
void Operator_FT_Disc::Operator_Laplacian_FT_element(const ArrOfDouble& Phi_Facet,const Maillage_FT_Disc& FTmesh, ArrOfDouble& Laplacian_Phi_Facet,DoubleTab& Grad_Phi_Sommet)
{
  Operator_Gradient_FT_sommets(Phi_Facet, FTmesh, Grad_Phi_Sommet, true);
  const int dim = Objet_U::dimension;
  int nbfa7=FTmesh.nb_facettes();
  DoubleTab sommets=FTmesh.sommets();
  IntTab facettes=FTmesh.facettes();
  const DoubleTab& nfa7 = FTmesh.get_update_normale_facettes();
  const Desc_Structure_FT& desc_facettes = FTmesh.desc_facettes();
  Laplacian_Phi_Facet.resize(nbfa7);
  for (int fa7=0 ; fa7<nbfa7 ; fa7++)
    Laplacian_Phi_Facet[fa7]=0.;

  for (int fa7=0 ; fa7<nbfa7 ; fa7++)
    {
      if(! FTmesh.facette_virtuelle(fa7))
        {
          ArrOfDouble nfac(dim);
          nfac[0]=nfa7(fa7,0);
          nfac[1]=nfa7(fa7,1);
          nfac[2]=nfa7(fa7,2);

          // positions des sommets de la fa7
          DoubleTab x_sommets(dim,dim);
          for (int sommet_fa7=0 ; sommet_fa7<dim ; sommet_fa7++)
            {
              int indice_sommet = facettes(fa7,sommet_fa7);
              for (int dir=0 ; dir<dim ; dir++)
                x_sommets(sommet_fa7,dir)=sommets(indice_sommet,dir);
            }
          // positions du barycentre de la fa7
          ArrOfDouble x_g(dim);
          for (int dir=0 ; dir<dim ; dir++)
            x_g[dir] = (x_sommets(0, dir)+x_sommets(1, dir)+x_sommets(2, dir))/dim;

          if (dim==3)
            {
              for (int sommet_fa7=0 ; sommet_fa7<dim ; sommet_fa7++)
                {
                  ArrOfDouble p1(dim);
                  ArrOfDouble t1(dim);
                  ArrOfDouble x_midpoint(dim);
                  ArrOfDouble Grad_Phi_midpoint(dim);
                  for (int dir=0 ; dir<dim ; dir++)
                    {
                      t1[dir]=x_sommets(sommet_fa7,dir)-x_sommets(((sommet_fa7-1)%dim+dim)%dim,dir);
                      x_midpoint[dir]=(x_sommets(sommet_fa7,dir)+x_sommets(((sommet_fa7-1)%dim+dim)%dim,dir))/2.;
                    }

                  int indice_sommet = facettes(fa7,sommet_fa7);
                  int other_sommet = facettes(fa7,((sommet_fa7-1)%dim+dim)%dim);
                  for (int dir=0 ; dir<dim ; dir++)
                    Grad_Phi_midpoint(dir)=(Grad_Phi_Sommet(indice_sommet, dir)+Grad_Phi_Sommet(other_sommet, dir))/2.;

                  produit_vectoriel(nfac,t1,p1);

                  // il faut s'assurer que p1 soit dans la direction exterieur a l element de surface de la fa7.
                  // il faut donc que p1.(x_midpoint-x_g)>0
                  // sinon, il faut inverser p1
                  if (p1[0]*(x_midpoint[0]-x_g[0])+p1[1]*(x_midpoint[1]-x_g[1])+p1[2]*(x_midpoint[2]-x_g[2])<0)
                    {
                      for (int dir=0 ; dir<dim ; dir++)
                        p1[dir]=-p1[dir];
                    }
                  for (int dir=0 ; dir<dim ; dir++)
                    Laplacian_Phi_Facet[fa7]+=Grad_Phi_midpoint(dir)*p1[dir];
                }
            }
        }
    }
  desc_facettes.echange_espace_virtuel(Laplacian_Phi_Facet);
}

void Operator_FT_Disc::Compute_interfaciale_source(const ArrOfDouble& sigma_Facet, const Maillage_FT_Disc& FTmesh,
                                                   DoubleTab& df_sigma, bool Normalised_with_Surface, bool use_tryggvason_formulation, bool with_marangoni)
{
  int dim = Objet_U::dimension;
  int nbfa7=FTmesh.nb_facettes();
  int nbsom=FTmesh.nb_sommets();
  //DoubleTab sommets=FTmesh.sommets();
  // IntTab facettes=FTmesh.facettes();
  //const ArrOfDouble& Sfa7 = FTmesh.get_update_surface_facettes();
  const Desc_Structure_FT& desc_sommets = FTmesh.desc_sommets();

  /* initialisation des tableaux locaux */
  ArrOfDouble Unit_Facet, Unit_Somm, Surface_sommet, sigma_sommet ;
  DoubleTab df_sigma_bis ;
  df_sigma.resize(nbsom, dim);
  df_sigma_bis.resize(nbsom, dim);
  Unit_Somm.resize(nbsom);
  Unit_Facet.resize(nbfa7);
  sigma_sommet.resize(nbsom);
  Surface_sommet.resize(nbsom);
  for (int fa7=0 ; fa7<nbfa7 ; fa7++)
    {
      Unit_Facet(fa7)=1. ;
    }
  for (int som=0 ; som<nbsom ; som++)
    {
      Unit_Somm(som)=1.;
      sigma_sommet(som)=0.;
      for (int dir=0 ; dir<dim ; dir++)
        {
          df_sigma(som, dir)= 0. ;
          df_sigma_bis(som, dir)= 0. ;
        }
    }
  /* Interpolation aux sommets */
  Facette_to_Sommets(Surface_sommet, sigma_sommet, sigma_Facet, FTmesh, true);
  Operator_integral_bord_facette_phi_p_dl(Surface_sommet, sigma_sommet, sigma_Facet, FTmesh, df_sigma, df_sigma_bis);

  if (!with_marangoni and use_tryggvason_formulation)
    {
      // on ne veut pas du terme de marangoni
      // Cela revient à extraire sigma de l'integrale
      // on veut donc sigma * int_p_dl
      for (int som=0 ; som<nbsom ; som++)
        if(! FTmesh.sommet_virtuel(som))
          for (int dir=0 ; dir<3 ; dir++)
            df_sigma(som, dir) = df_sigma_bis(som, dir)*sigma_sommet[som];

    }
  else if(with_marangoni and !use_tryggvason_formulation)
    {
      // on ne veut QUE du terme de marangoni
      // Cela revient à soustraire les deux contribution int_sigma_p_dl - sigma * int_p_dl
      for (int som=0 ; som<nbsom ; som++)
        if(! FTmesh.sommet_virtuel(som))
          for (int dir=0 ; dir<3 ; dir++)
            df_sigma(som, dir)-=df_sigma_bis(som, dir)*sigma_sommet[som];

    }
  else if (!with_marangoni and !use_tryggvason_formulation)
    {
      // on ne veut rien du tout
      for (int som=0 ; som<nbsom ; som++)
        for (int dir=0 ; dir<dim ; dir++)
          df_sigma(som, dir)= 0. ;
    }
  // dans le cas with_marangoni and use_tryggvason_formulation, on veut df_sigma=int_sigma_p_dl
  // Il n'y a donc rien a changer

  desc_sommets.echange_espace_virtuel(df_sigma);
}

void Operator_FT_Disc::Operator_Gradient_FT_sommets(const ArrOfDouble& Phi_Facet, const Maillage_FT_Disc& FTmesh,
                                                    DoubleTab& Grad_Phi_Sommet, bool Normalised_with_Surface)
{


  int dim = Objet_U::dimension;
  int nbfa7=FTmesh.nb_facettes();
  int nbsom=FTmesh.nb_sommets();
  //DoubleTab sommets=FTmesh.sommets();
  //IntTab facettes=FTmesh.facettes();
  //const ArrOfDouble& Sfa7 = FTmesh.get_update_surface_facettes();
  const Desc_Structure_FT& desc_sommets = FTmesh.desc_sommets();

  /* initialisation des tableaux locaux */
  ArrOfDouble Unit_Facet, Unit_Somm, Phi_sommet, Surface_sommet ;
  DoubleTab int_phi_p_dl, int_p_dl ;
  Grad_Phi_Sommet.resize(nbsom, dim);
  int_phi_p_dl.resize(nbsom, dim);
  int_p_dl.resize(nbsom, dim);
  Unit_Somm.resize(nbsom);
  Unit_Facet.resize(nbfa7);
  Phi_sommet.resize(nbsom);
  Surface_sommet.resize(nbsom);
  for (int fa7=0 ; fa7<nbfa7 ; fa7++)
    {
      Unit_Facet(fa7)=1. ;
    }
  for (int som=0 ; som<nbsom ; som++)
    {
      Unit_Somm(som)=1.;
      for (int dir=0 ; dir<dim ; dir++)
        {
          Grad_Phi_Sommet(som, dir)= 0. ;
          int_p_dl(som, dir)=0.;
          int_phi_p_dl(som, dir)=0.;
        }
    }

  /* Interpolation aux sommets */
  Facette_to_Sommets(Surface_sommet, Phi_sommet, Phi_Facet, FTmesh, true);

  /* Calcul des integrales surfaciques */
  Operator_integral_bord_facette_phi_p_dl(Surface_sommet, Phi_sommet, Phi_Facet, FTmesh, int_phi_p_dl, int_p_dl);

  /* Assemblage du gradient */

  for (int som=0 ; som<nbsom ; som++)
    if(! FTmesh.sommet_virtuel(som))
      for (int dir=0 ; dir<dim ; dir++)
        Grad_Phi_Sommet(som, dir) = int_phi_p_dl(som, dir)  - int_p_dl(som, dir) *Phi_sommet[som];

  /* normalisation si souhaitee */
  if (!Normalised_with_Surface)
    for (int som=0 ; som<nbsom ; som++)
      if(! FTmesh.sommet_virtuel(som))
        for (int dir=0 ; dir<dim ; dir++)
          Grad_Phi_Sommet(som, dir)*= Surface_sommet[som];

  desc_sommets.echange_espace_virtuel(Grad_Phi_Sommet);
}

void Operator_FT_Disc::Operator_integral_bord_facette_phi_p_dl(const ArrOfDouble& Surface_sommet, const ArrOfDouble& Phi_sommet, const ArrOfDouble& Phi_Facet, const Maillage_FT_Disc& FTmesh,
                                                               DoubleTab& int_phi_p_dl, DoubleTab& int_p_dl)
{
  int dim = Objet_U::dimension;
  double R = 1.; // TODO GUILLAUME
  double pi = 3.1415; // TODO GUILLAUME
  int nbfa7=FTmesh.nb_facettes();
  int nbsom=FTmesh.nb_sommets();
  DoubleTab sommets=FTmesh.sommets();
  IntTab facettes=FTmesh.facettes();
  const DoubleTab& nfa7 = FTmesh.get_update_normale_facettes();
  const Desc_Structure_FT& desc_sommets = FTmesh.desc_sommets();
// Calcul du gradient aux sommets
  for (int fa7=0 ; fa7<nbfa7 ; fa7++)
    {
      if(! FTmesh.facette_virtuelle(fa7))
        {
          ArrOfDouble nfac(dim);
          for (int dir=0 ; dir<dim ; dir++)
            nfac[dir]=nfa7(fa7,dir);

          // positions du barycentre de la fa7
          ArrOfDouble x_g(dim);
          for (int dir=0 ; dir<dim ; dir++)
            {
              for (int sommet_fa7=0 ; sommet_fa7<dim ; sommet_fa7++)
                {
                  int indice_sommet = facettes(fa7,sommet_fa7);
                  x_g[dir] += sommets(indice_sommet, dir);
                }
              x_g[dir]/=dim;
            }

          // positions des 2 midpoints adjacents a chaques sommets de la fa7
          // Moyenne de Phi a ces midpoints
          DoubleTab x_midpoint1(dim,dim);
          DoubleTab x_midpoint2(dim,dim);
          ArrOfDouble Phi_midpoint1(dim);
          ArrOfDouble Phi_midpoint2(dim);

          for (int sommet_fa7=0 ; sommet_fa7<dim ; sommet_fa7++)
            {
              int indice_sommet = facettes(fa7,sommet_fa7);
              int indice_sommet_second = facettes(fa7,((sommet_fa7-1)%dim+dim)%dim);
              Phi_midpoint1[sommet_fa7] = (Phi_sommet[indice_sommet]+Phi_sommet[indice_sommet_second])/2.;
              for (int dir=0 ; dir<dim ; dir++)
                {
                  x_midpoint1(sommet_fa7, dir)=(sommets(indice_sommet, dir)+sommets(indice_sommet_second, dir))/2.;
                }
              if(dim==3)
                {
                  int indice_sommet_third= facettes(fa7,((sommet_fa7+1)%dim+dim)%dim);
                  Phi_midpoint2[sommet_fa7] = (Phi_sommet[indice_sommet]+Phi_sommet[indice_sommet_third])/2.;
                  for (int dir=0 ; dir<dim ; dir++)
                    {
                      x_midpoint2(sommet_fa7, dir)=(sommets(indice_sommet, dir)+sommets(indice_sommet_third, dir))/2.;
                    }
                }
            }
          // calcul de l integrale lineique de Phi.p sur le bord de la sous zone e du volument du controle du sommet.
          // La sous-zone e est definie par le croisement des mediane de l element triangulaire.
          // Le bord de la sous-zone est composee de 2 segments --> int = Phi1.p1.DS1 + Phi2.p2.DS2
          // Voir Muradoglu et Tryggvason 2014 pour plus de details.

          for (int sommet_fa7=0 ; sommet_fa7<dim ; sommet_fa7++)
            {
              int indice_sommet = facettes(fa7,sommet_fa7);
              ArrOfDouble p1(dim), p2(dim);
              ArrOfDouble t1(dim), t2(dim);
              double Phi1 = 0.;
              double Phi2 = 0.;
              Phi1 = (Phi_Facet(fa7)+Phi_midpoint1(sommet_fa7))/2.;
              if (dim==3)
                Phi2 = (Phi_Facet(fa7)+Phi_midpoint2(sommet_fa7))/2.;

              if (dim==3)
                {
                  for (int dir=0 ; dir<dim ; dir++)
                    {
                      t1[dir]=x_midpoint1(sommet_fa7,dir)-x_g[dir];
                      t2[dir]=x_midpoint2(sommet_fa7,dir)-x_g[dir];
                    }
                  produit_vectoriel(nfac,t1,p1);
                  produit_vectoriel(nfac,t2,p2);
                  // il faut s'assurer que p1 et p2 sont dans la direction exterieur a la surface de controle du sommet.
                  // il faut donc que p1.(xg-x_sommet)>0
                  // sinon, il faut inverser p1, idem pour p2
                  // Uniquement en 3D
                  if (p1[0]*(x_g[0]-sommets(indice_sommet,0))+p1[1]*(x_g[1]-sommets(indice_sommet,1))+p1[2]*(x_g[2]-sommets(indice_sommet,2))<0)
                    {
                      for (int dir=0 ; dir<dim ; dir++)
                        p1[dir]=-p1[dir];
                    }
                  if (p2[0]*(x_g[0]-sommets(indice_sommet,0))+p2[1]*(x_g[1]-sommets(indice_sommet,1))+p2[2]*(x_g[2]-sommets(indice_sommet,2))<0)
                    {
                      for (int dir=0 ; dir<dim ; dir++)
                        p2[dir]=-p2[dir];
                    }
                }
              else if (dim==2)
                {
                  for (int dir=0 ; dir<dim ; dir++)
                    p1[dir]=x_midpoint1(sommet_fa7,dir)-sommets(indice_sommet,dir);
                  unitarisation(p1);
                  for (int dir=0 ; dir<dim ; dir++)
                    p1[dir]*=2*pi*R;
                }

              for (int dir=0 ; dir<dim ; dir++)
                {
                  int_phi_p_dl(indice_sommet, dir)+=Phi1*p1[dir];
                  if (dim==3)
                    int_phi_p_dl(indice_sommet, dir)+=Phi2*p2[dir];
                  int_p_dl(indice_sommet, dir)+=p1[dir];
                  if (dim==3)
                    int_p_dl(indice_sommet, dir)+=p2[dir];
                }
            }
        }
    }

  desc_sommets.collecter_espace_virtuel(int_phi_p_dl, MD_Vector_tools::EV_SOMME);
  desc_sommets.collecter_espace_virtuel(int_p_dl, MD_Vector_tools::EV_SOMME);

  for (int som=0 ; som<nbsom ; som++)
    if(! FTmesh.sommet_virtuel(som))
      for (int dir=0 ; dir<dim ; dir++)
        {
          int_phi_p_dl(som, dir)/= Surface_sommet[som];
          int_p_dl(som, dir)/= Surface_sommet[som];
        }

  desc_sommets.echange_espace_virtuel(int_phi_p_dl);
  desc_sommets.echange_espace_virtuel(int_p_dl);
}


void Operator_FT_Disc::Facette_to_Sommets(ArrOfDouble& Surface_sommet, DoubleTab& Phi_Som, const DoubleTab& Phi_Facet, const Maillage_FT_Disc& FTmesh, bool Normalised_with_Surface)
{
  int dim = Objet_U::dimension;
  int nbfa7=FTmesh.nb_facettes();
  int nbsom=FTmesh.nb_sommets();
  ArrOfDouble Phi_Facet_dir ;
  ArrOfDouble Phi_Som_dir ;
  Phi_Facet_dir.resize(nbfa7);
  Phi_Som_dir.resize(nbsom);
  Phi_Som.resize(nbsom, dim);

  for (int dir=0 ; dir<dim ; dir++)
    {
      for (int fa7=0 ; fa7<nbfa7 ; fa7++)
        Phi_Facet_dir(fa7) = Phi_Facet(fa7, dir);

      Facette_to_Sommets(Surface_sommet, Phi_Som_dir, Phi_Facet_dir, FTmesh, Normalised_with_Surface);

      for (int som=0 ; som<nbsom ; som++)
        Phi_Som(som, dir) = Phi_Som_dir(som) ;
    }
}

void Operator_FT_Disc::Facette_to_Sommets(ArrOfDouble& Surface_sommet, ArrOfDouble& Phi_Som, const ArrOfDouble& Phi_Facet, const Maillage_FT_Disc& FTmesh, bool Normalised_with_Surface)
{
  int dim = Objet_U::dimension;
  int nbfa7=FTmesh.nb_facettes();
  int nbsom=FTmesh.nb_sommets();
  IntTab facettes=FTmesh.facettes();
  const ArrOfDouble& Sfa7 = FTmesh.get_update_surface_facettes();
  Phi_Som.resize(nbsom);
  Surface_sommet.resize(nbsom);

  for (int som=0 ; som<nbsom ; som++)
    {
      Phi_Som[som]=0.;
      Surface_sommet[som]=0.;
    }

  // interpolation des valeur de Phi et de la normale aux sommets en moyennant les contributions adjacentes
  for (int fa7=0 ; fa7<nbfa7 ; fa7++)
    {
      if(! FTmesh.facette_virtuelle(fa7))
        {
          for (int sommet_fa7=0 ; sommet_fa7<dim ; sommet_fa7++)
            {
              int indice_sommet = facettes(fa7,sommet_fa7);
              Phi_Som[indice_sommet]+=Phi_Facet[fa7]*Sfa7[fa7]/dim; // TODO Guillaume : il faut que Sfa7 en 2D contienne * 2 pi R ; sinon il faut le rajouter ?
              Surface_sommet[indice_sommet]+=Sfa7[fa7]/dim; // TODO Guillaume : il faut que Sfa7 en 2D contienne * 2 pi R ; sinon il faut le rajouter ?
            }
        }
    }

  // On a calcule la contribution de chaque facette reelle aux differents sommets.
  // Certaines contributions ont ete ajoutees a des sommets virtuels, il
  // faut recuperer ces contributions sur le sommet reel.
  const Desc_Structure_FT& desc_sommets = FTmesh.desc_sommets();
  desc_sommets.collecter_espace_virtuel(Phi_Som, MD_Vector_tools::EV_SOMME);
  desc_sommets.collecter_espace_virtuel(Surface_sommet, MD_Vector_tools::EV_SOMME);

  if (Normalised_with_Surface)
    for (int som=0 ; som<nbsom ; som++)
      {
        if(! FTmesh.sommet_virtuel(som))
          {
            Phi_Som[som]/= Surface_sommet[som];
          }
      }


  desc_sommets.echange_espace_virtuel(Phi_Som);
}

void Operator_FT_Disc::Sommets_to_Facettes(DoubleTab& Phi_Facet, const DoubleTab& Phi_Som, const Maillage_FT_Disc& FTmesh, bool Normalised_with_Surface)
{
  int dim = Objet_U::dimension;
  int nbfa7=FTmesh.nb_facettes();
  int nbsom=FTmesh.nb_sommets();
  ArrOfDouble Phi_Facet_dir ;
  ArrOfDouble Phi_Som_dir ;
  Phi_Facet_dir.resize(nbfa7);
  Phi_Som_dir.resize(nbsom);
  Phi_Facet.resize(nbfa7, dim);

  for (int dir=0 ; dir<dim ; dir++)
    {
      for (int som=0 ; som<nbsom ; som++)
        Phi_Som_dir(som) = Phi_Som(som, dir);

      Sommets_to_Facettes(Phi_Facet_dir, Phi_Som_dir, FTmesh, Normalised_with_Surface);

      for (int fa7=0 ; fa7<nbfa7 ; fa7++)
        Phi_Facet(fa7, dir)=Phi_Facet_dir(fa7);
    }
}

void Operator_FT_Disc::Sommets_to_Facettes(ArrOfDouble& Phi_Facet, const ArrOfDouble& Phi_Som, const Maillage_FT_Disc& FTmesh, bool Normalised_with_Surface)
{
  int dim = Objet_U::dimension;
  int nbfa7=FTmesh.nb_facettes();
  IntTab facettes=FTmesh.facettes();
  const ArrOfDouble& Sfa7 = FTmesh.get_update_surface_facettes();

  Phi_Facet.resize(nbfa7);
  for (int fa=0 ; fa<nbfa7 ; fa++)
    Phi_Facet[fa]=0.;

  // interpolation des valeur de Phi et de la normale aux sommets en moyennant les contributions adjacentes
  for (int fa7=0 ; fa7<nbfa7 ; fa7++)
    {
      if(! FTmesh.facette_virtuelle(fa7))
        {
          for (int sommet_fa7=0 ; sommet_fa7<dim ; sommet_fa7++)
            {
              int indice_sommet = facettes(fa7,sommet_fa7);
              Phi_Facet[fa7]+= Phi_Som[indice_sommet]/dim ;
            }
          if(!Normalised_with_Surface)
            Phi_Facet[fa7]*=Sfa7[fa7];
        }
    }
  const Desc_Structure_FT& desc_facettes = FTmesh.desc_facettes();
  desc_facettes.echange_espace_virtuel(Phi_Facet);
}


void Operator_FT_Disc::produit_vectoriel(const ArrOfDouble& a, const ArrOfDouble& b, ArrOfDouble& resu)
{
  if (Objet_U::dimension == 3)
    {
      resu[0] = a[1]*b[2] - a[2]*b[1];
      resu[1] = a[2]*b[0] - a[0]*b[2];
      resu[2] = a[0]*b[1] - a[1]*b[0];
    }
}

double Operator_FT_Disc::norme(const ArrOfDouble& a)
{
  double z =(Objet_U::dimension == 3) ? a[2]*a[2] : 0.;
  return std::sqrt(a[0]*a[0]+a[1]*a[1]+z);
}
void Operator_FT_Disc::unitarisation(ArrOfDouble& a)
{
  double normea=norme(a);
  a[0] /= normea ;
  a[1] /= normea ;
  if (Objet_U::dimension == 3)
    a[2] /= normea ;
}
