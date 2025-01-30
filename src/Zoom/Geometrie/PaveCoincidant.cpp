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
// File:        PaveCoincidant.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Geometrie
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////

#include <PaveCoincidant.h>
#include <math.h>
#include <Motcle.h>
#include <Domaine.h>
#include <EChaine.h>
#include <Interprete.h>

Implemente_instanciable(PaveCoincidant,"PaveCoincidant",Pave);

/*! @brief Simple appel a: Domaine::printOn(Sortie&)
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& PaveCoincidant::printOn(Sortie& s ) const
{
  return Domaine::printOn(s) ;
}


/*! @brief Lit les specifications d'un pave a partir d'un flot d'entree.
 *
 *     Le format de lecture d'un pave dans le jeu de donnee est le suivant:
 *      PaveCoincidant nom_pave
 *      {
 *      Origine OX OY (OZ)
 *      Longueurs LX LY (LZ)
 *      Nombre_de_noeuds NX NY (NZ)
 *      Facteurs Fx Fy (Fz)
 *      (Symx)
 *      (Symy)
 *      (Symz)
 *      }
 *      {
 *      (Bord)  nom X = X0 Y0 <= Y <= Y1 Z0 <= Z <= Z1
 *      ...
 *      (Raccord)  local homogene nom X = X0 Y0 <= Y <= Y1 Z0 <= Z <= Z1
 *      ...
 *      (Internes)  nom X = X0 Y0 <= Y <= Y1 Z0 <= Z <= Z1
 *      ...
 *      (Joint)  nom X = X0 Y0 <= Y <= Y1 Z0 <= Z <= Z1 PE_voisin
 *      ...
 *      }
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 * @throws dimension d'espace necessaire pour mailler
 * @throws accolade ouvrante attendue
 * @throws Symy n'a de sens que pour une dimension >= 2
 * @throws Symz n'a de sens qu'en dimension 3
 * @throws Les facteurs de progression doivent etre positifs
 * @throws Il doit y avoir au moins deux mailles en x
 * @throws accolade oubvrante attendue avant lecture des bords
 * @throws mot cle non reconnu
 */
Entree& PaveCoincidant::readOn(Entree& is)
{
  Nom ledom,bh,bb,bd,bg,bder,bdev;

  is >> nom_;
  is >> ledom;
  Cerr << "Lecture du pave " << nom_ << " qui doit coincider avec le domaine " << ledom << finl;

  dom_qui_coincide=ref_cast(Domaine, Interprete::objet(ledom));

  Motcle motlu;
  int rang,nmx,nmy,nmz;
  is >> motlu;
  if (motlu!="{")
    {
      Cerr << "On attendait une { apres " << nom_ << finl;
      Process::exit();
    }

  Motcles les_mots(14);
  {
    les_mots[0]="xmin";
    les_mots[1]="xmax";
    les_mots[2]="ymin";
    les_mots[3]="ymax";
    les_mots[4]="zmin";
    les_mots[5]="zmax";
    les_mots[6]="Bord_Haut";
    les_mots[7]="Bord_Bas";
    les_mots[8]="Bord_gauche";
    les_mots[9]="Bord_droit";
    les_mots[10]="Bord_devant";
    les_mots[11]="Bord_derriere";
    les_mots[12]="Nombre_de_Noeuds";
    les_mots[13]="}";
  }
  while(motlu != "}")
    {
      is >> motlu;
      rang = les_mots.search(motlu);
      switch(rang)
        {
        case 0:
          is >> xmin;
          break;
        case 1:
          is >> xmax;
          break;
        case 2:
          is >> ymin;
          break;
        case 3:
          is >> ymax;
          break;
        case 4:
          is >> zmin;
          break;
        case 5:
          is >> zmax;
          break;
        case 6:
          is >> bh;
          break;
        case 7:
          is >> bb;
          break;
        case 8:
          is >> bg;
          break;
        case 9:
          is >> bd;
          break;
        case 10:
          is >> bdev;
          break;
        case 11:
          is >> bder;
          break;
        case 12:
          is >> nmx >> nmy;
          if (dimension ==3) is >> nmz;
          break;
        case 13: // Fin
          break;
        default:
          Cerr << motlu << "  n'est pas compris " << finl;
          Cerr << les_mots;
          Process::exit();
        }
    }

  determine_bornes();

  //char* str=nom;
  //strcat(str," { Origine " );

  //EChaine le_pave(str);
  Nom desc_pave;
  if (dimension == 2)
    {
      Nom nom_nx(nmx);
      Nom nom_ny(nmy);
      desc_pave = nom_;
      desc_pave+=" { Origine ";
      desc_pave+= (Nom)xmin;
      desc_pave+= " ";
      desc_pave+= (Nom)ymin;
      desc_pave+= " Nombre_de_Noeuds ";
      desc_pave+= nom_nx;
      desc_pave+= " ";
      desc_pave+= nom_ny;
      desc_pave+=" Longueurs ";
      desc_pave+= (Nom)(xmax - xmin);
      desc_pave+= " ";
      desc_pave+= (Nom)(ymax - ymin);
      desc_pave+= " } { bord " ;
      desc_pave+= (Nom)bg;
      desc_pave+=" X = ";
      desc_pave+= (Nom)xmin;
      desc_pave+=" ";
      desc_pave+= (Nom)ymin;
      desc_pave+="  <= Y <= ";
      desc_pave+= (Nom)ymax;
      desc_pave+="  bord ";
      desc_pave+= (Nom)bh;
      desc_pave+="   Y = ";
      desc_pave+= (Nom)ymax;
      desc_pave+=" ";
      desc_pave+= (Nom)xmin;
      desc_pave+=" <= X <= ";
      desc_pave+= (Nom)xmax;
      desc_pave+=" bord ";
      desc_pave+= (Nom)bb;
      desc_pave+="    Y = ";
      desc_pave+= (Nom)ymin;
      desc_pave+=" ";
      desc_pave+= (Nom)xmin;
      desc_pave+=" <= X <= ";
      desc_pave+= (Nom)xmax;
      desc_pave+=" bord ";
      desc_pave+= (Nom)bd;
      desc_pave+=" X = ";
      desc_pave+= (Nom)xmax;
      desc_pave+=" ";
      desc_pave+= (Nom)ymin;
      desc_pave+=" <= Y <= ";
      desc_pave+= (Nom)ymax;
      desc_pave+=" }";
    }
  else if (dimension == 3)
    {
      Nom nom_nx(nmx);
      Nom nom_ny(nmy);
      Nom nom_nz(nmz);
      desc_pave = nom_;
      desc_pave+= " { Origine " ;
      desc_pave+= (Nom)xmin ;
      desc_pave+= " " ;
      desc_pave+= (Nom)ymin ;
      desc_pave+= " " ;
      desc_pave+=(Nom)zmin;
      desc_pave+= " Nombre_de_Noeuds " ;
      desc_pave+=   nom_nx ;
      desc_pave+= " " ;
      desc_pave+= nom_ny ;
      desc_pave+= " " ;
      desc_pave+= nom_nz ;
      desc_pave+= " Longueurs " ;
      desc_pave+= (Nom)(xmax - xmin) ;
      desc_pave+= " " ;
      desc_pave+= (Nom)(ymax - ymin) ;
      desc_pave+= " " ;
      desc_pave+= (Nom)(zmax - zmin);
      desc_pave+= " } { bord " ;
      desc_pave+= (Nom)bg ;
      desc_pave+= " X = " ;
      desc_pave+= (Nom)xmin ;
      desc_pave+= " " ;
      desc_pave+= (Nom)ymin ;
      desc_pave+= "  <= Y <= " ;
      desc_pave+= (Nom)ymax ;
      desc_pave+= " " ;
      desc_pave+=(Nom)zmin;
      desc_pave+="  <= Z <= " ;
      desc_pave+= (Nom)zmax ;
      desc_pave+="  bord " ;
      desc_pave+= (Nom)bh ;
      desc_pave+= "   Z = " ;
      desc_pave+= (Nom)zmax  ;
      desc_pave+= " " ;
      desc_pave+= (Nom)xmin ;
      desc_pave+=" <= X <= " ;
      desc_pave+= (Nom)xmax ;
      desc_pave+= " " ;
      desc_pave+= (Nom)ymin ;
      desc_pave+=" <= Y <= " ;
      desc_pave+= (Nom)ymax ;

      //desc_pave= desc_pave;
      desc_pave+=" bord " ;
      desc_pave+= (Nom)bb ;
      desc_pave+= "    Z = " ;
      desc_pave+= (Nom)zmin ;
      desc_pave+= " " ;
      desc_pave+= (Nom)xmin  ;
      desc_pave+= " <= X <= " ;
      desc_pave+= (Nom)xmax ;
      desc_pave+= " " ;
      desc_pave+= (Nom)ymin  ;
      desc_pave+= " <= Y <= " ;
      desc_pave+= (Nom)ymax ;
      desc_pave+= " bord " ;
      desc_pave+= (Nom)bd ;
      desc_pave+=  " X = " ;
      desc_pave+= (Nom)xmax ;
      desc_pave+= " " ;
      desc_pave+= (Nom)ymin ;
      desc_pave+= " <= Y <= " ;
      desc_pave+= (Nom)ymax ;
      desc_pave+= " " ;
      desc_pave+=(Nom)zmin;
      desc_pave+= " <= Z <= " ;
      desc_pave+= (Nom)zmax ;
      desc_pave+= " bord " ;
      desc_pave+= (Nom)bdev ;
      desc_pave+= " Y = ";
      desc_pave+=(Nom)ymin;
      desc_pave+=" " ;
      desc_pave+= (Nom)xmin ;
      desc_pave+= "  <= X <= " ;
      desc_pave+= (Nom)xmax ;
      desc_pave+= " " ;
      desc_pave+= (Nom)zmin ;
      desc_pave+= "  <= Z <= " ;
      desc_pave+= (Nom)zmax ;
      desc_pave+=  " bord ";
      desc_pave+=(Nom)bder;
      desc_pave+=" Y = ";
      desc_pave+=(Nom)ymax;
      desc_pave+=" " ;
      desc_pave+= (Nom)xmin ;
      desc_pave+= "  <= X <= " ;
      desc_pave+= (Nom)xmax ;
      desc_pave+= " " ;
      desc_pave+=(Nom)zmin ;
      desc_pave+= "  <= Z <= " ;
      desc_pave+= (Nom)zmax ;
      desc_pave+=" }";

    }

  EChaine le_pave(desc_pave);

  //Cerr << " Chaine creee = " << nom+" { Origine " + xmin + " " + ymin + " Nombre_de_Noeuds 17 17 Longueurs " + (xmax - xmin) + " " + (ymax - ymin) + " } { bord " + bg + " X = " + xmin + " " + ymin + "  <= Y <= " + ymax + "  bord " + bh + "   Y = " + ymax  + " " + xmin +" <= X <= " + xmax + " bord " + bb + "    Y = " + ymin + " " + xmin  + " <= X <= " + xmax + " bord " + bd +  " X = " + xmax + " " + ymin + " <= Y <= " + ymax + " }" << finl;

  Pave::readOn(le_pave);

  return is;
}

void PaveCoincidant::determine_bornes()
{
  if(dom_qui_coincide->type_elem()->que_suis_je() == "Rectangle" || dom_qui_coincide->type_elem()->que_suis_je() == "Hexaedre" )
    {
      if (dimension == 2)
        dom_qui_coincide->typer("Rectangle");
      else
        dom_qui_coincide->typer("Hexaedre");

      IntTab& leselems=dom_qui_coincide->les_elems();
      int oldsz=leselems.dimension(0);
      const DoubleTab& les_coord = dom_qui_coincide->coord_sommets();


      double xm=0, xM=0, ym=0, yM=0, zm=0, zM=0; // bornes trouvees coincidantes avec le domaine souhaite
      double dxm, dxM, dym, dyM, dzm=0, dzM=0; // variables tempos pour determiner les bornes


      // init
      dxm = std::fabs(xmin-les_coord(leselems(0,0),0));
      dxM = std::fabs(xmin-les_coord(leselems(0,0),0));
      dym = std::fabs(xmin-les_coord(leselems(0,0),1));
      dyM = std::fabs(xmin-les_coord(leselems(0,0),1));
      if (dimension == 3)
        {
          dzm = std::fabs(xmin-les_coord(leselems(0,0),2));
          dzM = std::fabs(xmin-les_coord(leselems(0,0),2));
        }

      // Determination du nouveau nb d'elements
      for(int i=0; i< oldsz; i++)
        {
          DoubleTab xg(Objet_U::dimension);
          for(int j=0; j<dom_qui_coincide->nb_som_elem(); j++)
            for(int k=0; k<Objet_U::dimension; k++)
              xg(k)+=les_coord(leselems(i,j),k);
          xg/=dom_qui_coincide->nb_som_elem();

          if (xg(0) > xmin && xg(0) < xmax && xg(1) < ymax && xg(1) > ymin ) // si le centre de gravite de l element est dans le domaine a zoomer ...
            {
              if ( dxm > std::fabs(xmin-les_coord(leselems(i,0),0)) )
                {
                  dxm = std::fabs(xmin-les_coord(leselems(i,0),0));
                  xm = les_coord(leselems(i,0),0);
                }
              if ( dxM > std::fabs(xmax-les_coord(leselems(i,3),0)) )
                {
                  dxM = std::fabs(xmax-les_coord(leselems(i,3),0));
                  xM = les_coord(leselems(i,3),0);
                }
              if ( dym > std::fabs(ymin-les_coord(leselems(i,0),1)) )
                {
                  dym = std::fabs(ymin-les_coord(leselems(i,0),1));
                  ym = les_coord(leselems(i,0),1);
                }
              if ( dyM > std::fabs(ymax-les_coord(leselems(i,3),1)) )
                {
                  dyM = std::fabs(ymax-les_coord(leselems(i,3),1));
                  yM = les_coord(leselems(i,3),1);
                }
              if (dimension == 3)
                {
                  if ( dzm > std::fabs(zmin-les_coord(leselems(i,0),2)) )
                    {
                      dzm = std::fabs(zmin-les_coord(leselems(i,0),2));
                      zm = les_coord(leselems(i,0),2);
                    }
                  if ( dzM > std::fabs(zmax-les_coord(leselems(i,3),2)) )
                    {
                      dzM = std::fabs(zmax-les_coord(leselems(i,3),2));
                      zM = les_coord(leselems(i,3),2);
                    }
                }

            }
        }

      xmin = xm;
      xmax = xM;
      ymin = ym;
      ymax = yM;
      zmin = zm;
      zmax = zM;
    }
  else
    {
      Cerr << "On ne sait pas encore creer un pave coincidant avec un "
           << dom_qui_coincide->type_elem()->que_suis_je() <<"s"<<finl;
    }
}


