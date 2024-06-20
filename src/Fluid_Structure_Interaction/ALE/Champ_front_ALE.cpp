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
// File:        Champ_front_ALE.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/ALE/src
// Version:     /main/17
//
//////////////////////////////////////////////////////////////////////////////

#include <Champ_front_ALE.h>
#include <Domaine.h>
#include <Frontiere_dis_base.h>
#include <MD_Vector_tools.h>

Implemente_instanciable(Champ_front_ALE,"Champ_front_ALE",Ch_front_var_instationnaire_dep);
// XD Champ_front_ale front_field_base Champ_front_ale 0 Class to define a boundary condition on a moving boundary of a mesh (only for the Arbitrary Lagrangian-Eulerian framework ).
// XD attr val listchaine val 0 NL2 Example:  2 -y*0.01 x*0.01

/*Champ_front_ALE::Champ_front_ALE()
{
  const Frontiere& front=la_frontiere_dis->frontiere();
  const Domaine& domaine=front.domaine();
  const Domaine& domaine=domaine.domaine();
  vit_som_bord_ALE.resize(domaine.nb_som(),nb_comp());
  const MD_Vector& md = domaine.domaine().md_vector_sommets();
  MD_Vector_tools::creer_tableau_distribue(md, vit_som_bord_ALE);
}*/

/*! @brief Impression sur un flot de sortie au format: taille
 *
 *     valeur(0) ... valeur(i)  ... valeur(taille-1)
 *
 * @param (Sortie& os) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Champ_front_ALE::printOn(Sortie& os) const
{
  //   const DoubleTab& tab=valeurs();
  //   os << tab.size() << " ";
  //   for(int i=0; i<tab.size(); i++)
  //     os << tab(0,i);
  return os;
}




/*! @brief Lecture a partir d'un flot d'entree au format: nombre_de_composantes
 *
 *     fonction dependant de xyzt pour chaque composante
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 * @throws accolade ouvrante attendue
 * @throws mot clef inconnu a cet endroit
 * @throws accolade fermante attendue
 */
Entree& Champ_front_ALE::readOn(Entree& is)
{
  int dim;
  is >> dim;
  fixer_nb_comp(dim);

  fxyzt.dimensionner(dim);


  //        Cout << "dim = " << dim << finl;
  for (int i = 0; i<dim; i++)
    {
      Nom tmp;
      //Cout << "i = " << i << finl;
      is >> tmp;
      //Cout << "fonc = " << tmp << finl;
      Cerr << "Lecture et interpretation de la fonction " << tmp << finl;
      fxyzt[i].setNbVar(4);
      fxyzt[i].setString(tmp);
      fxyzt[i].addVar("x");
      fxyzt[i].addVar("y");
      fxyzt[i].addVar("z");
      fxyzt[i].addVar("t");
      fxyzt[i].parseString();
      Cerr << "Interpretation de la fonction " << tmp << " Ok" << finl;
      //        Cout << "end = " << tmp << finl;
    }
  return is;
}


/*! @brief Pas code !!
 *
 * @param (Champ_front_base& ch)
 * @return (Champ_front_base&)
 */
Champ_front_base& Champ_front_ALE::affecter_(const Champ_front_base& ch)
{
  return *this;
}

int Champ_front_ALE::initialiser(double temps, const Champ_Inc_base& inco)
{
  if (!Ch_front_var_instationnaire_dep::initialiser(temps,inco))
    return 0;

  mettre_a_jour(temps);
  return 1;
}

/*! @brief Mise a jour du champ front et remplie un tableau de dimension egale au nombre total de sommets du domaine.
 *
 * Ce tableau a des
 *      valeurs nulles aux sommets qui n'appartiennent pas au bord ALE
 *      traite.
 *
 * @param (double tps) le temps de mise a jour
 */
void Champ_front_ALE::mettre_a_jour(double temps)
{
  //Cerr << "Champ_front_ALE ::mettre_a_jour" << finl;

  const Frontiere& front=la_frontiere_dis->frontiere();
  int nb_faces=front.nb_faces();
  const Faces& faces=front.faces();
  int nbsf=faces.nb_som_faces();
  int i,j,k;
  DoubleTab& tab=valeurs_au_temps(temps);
  tab=0.;

  // appel systematique a remplir_vit_som_bord_ALE()
  remplir_vit_som_bord_ALE(temps);

  for( i=0; i<nb_faces; i++)
    {
      for( j=0; j<nb_comp(); j++)
        {
          for( k=0; k<nbsf; k++)
            {
              tab(i,j)+=vit_som_bord_ALE(faces.sommet(i,k),j);
            }
        }
    }
  tab/=nbsf;
  tab.echange_espace_virtuel();
}

void Champ_front_ALE::remplir_vit_som_bord_ALE(double tps)
{
  //Cerr<<"Champ_front_ALE::remplir_vit_som_bord_ALE"<<finl;
  const Frontiere& front=la_frontiere_dis->frontiere();
  int nb_faces=front.nb_faces();
  const Domaine& domaine=front.domaine();
  const Faces& faces=front.faces();
  double x,y,z;
  int nbsf=faces.nb_som_faces();
  int i,j,k;
  int nb_som_tot=domaine.nb_som_tot();
  vit_som_bord_ALE.resize(nb_som_tot,nb_comp());
  /*if (vit_som_bord_ALE.dimension(0) != domaine.nb_som())
    {
      vit_som_bord_ALE.resize(domaine.nb_som(),nb_comp());
      const MD_Vector& md = domaine.domaine().md_vector_sommets();
      MD_Vector_tools::creer_tableau_distribue(md, vit_som_bord_ALE);
    }*/

  vit_som_bord_ALE=0.;

  for( i=0; i<nb_faces; i++)
    {
      x=y=z=0;
      for( k=0; k<nbsf; k++)
        {
          x=domaine.coord(faces.sommet(i,k),0);
          if(dimension>1)
            y=domaine.coord(faces.sommet(i,k),1);
          if(dimension>2)
            z=domaine.coord(faces.sommet(i,k),2);
          for( j=0; j<nb_comp(); j++)
            {
              fxyzt[j].setVar("x",x);
              fxyzt[j].setVar("y",y);
              fxyzt[j].setVar("z",z);
              fxyzt[j].setVar("t",tps);
              vit_som_bord_ALE(faces.sommet(i,k),j)=fxyzt[j].eval();
              //cout << " x y  " << x << " " << y << " " << z << " " << vit_som_bord_ALE(faces.sommet(i,k),j) << endl;
            }
        }
    }
  //vit_som_bord_ALE.echange_espace_virtuel();
}
