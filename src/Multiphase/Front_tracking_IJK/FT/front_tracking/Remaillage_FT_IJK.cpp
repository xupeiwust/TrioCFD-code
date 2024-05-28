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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Remaillage_FT_IJK.cpp
// Directory : $IJK_ROOT/src/FT/front_tracking
//
/////////////////////////////////////////////////////////////////////////////

#include <Remaillage_FT_IJK.h>
#include <Param.h>
#include <Maillage_FT_Disc.h>
#include <Comm_Group.h>
#include <stat_counters.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <variant>

using namespace std;

Implemente_instanciable_sans_constructeur(Remaillage_FT_IJK,"Remaillage_FT_IJK",Remaillage_FT) ;

Remaillage_FT_IJK::Remaillage_FT_IJK()
{
  nb_iter_barycentrage_ = 0; // Par defaut, pas de barycentrage
  nb_iter_bary_volume_seul_ = 0; // Pas de correction de volume
  nb_iter_remaillage_ = 1;
  facteur_longueur_ideale_ = -1; // Pas de redecoupage ou suppression d'arretes.
  //                                Sinon, la valeur recommandee en 3D est 1.5.
  critere_arete_ = 0.35;   // Je crois que c'est sa valeur par defaut en FTD. Je n'arrive pas a le retrouver.
  equilateral_ = 1;        // Par defaut, sur maillage etire, on fait des triangles equilateraux avec
  //                         un facteur_longueur qui vaut 1 si l'arrete a la meme taille que la diagonale.
  //                         Comportement different du FTD classique de trio (equilateral = 0) qui regarde
  //                         l'orientation de la facette pour ajuster sa taille a l'element eulerien.

}

Sortie& Remaillage_FT_IJK::printOn(Sortie& os) const
{
  //  Objet_U::printOn(os);
  os << "{\n"
     << "     pas_remaillage " << dt_remaillage_ << "\n"
     << "     nb_iter_barycentrage " << nb_iter_barycentrage_ << "\n"
     << "     relax_barycentrage " << relax_barycentrage_ << "\n";
  os << "     critere_arete " << critere_arete_ << "\n"
     << "     seuil_dvolume_residuel " << seuil_dvolume_residuel_ << "\n";
  os << "     nb_iter_correction_volume " << nb_iter_bary_volume_seul_ << "\n"
     << "     nb_iter_remaillage " << nb_iter_remaillage_ << "\n"
     << "     facteur_longueur_ideale " << facteur_longueur_ideale_ << "\n";
  os << "     equilateral " << equilateral_ << "\n"
     << "     lissage_courbure_coeff " << lissage_courbure_coeff_ << "\n"
     << "     lissage_courbure_iterations_systematique " << lissage_courbure_iterations_systematique_ << "\n"
     << "     lissage_courbure_iterations_si_remaillage " << lissage_courbure_iterations_si_remaillage_ << "\n";
  os << "   }\n" ;
  return os;
}
// XD remaillage_ft_ijk interprete nul 1 not_set
Entree& Remaillage_FT_IJK::readOn(Entree& is)
{
  Param p(que_suis_je());
  p.ajouter("pas_remaillage", &dt_remaillage_); // XD_ADD_P floattant not_set
  p.ajouter("nb_iter_barycentrage", &nb_iter_barycentrage_); // XD_ADD_P entier not_set
  p.ajouter("relax_barycentrage", &relax_barycentrage_); // XD_ADD_P floattant not_set
  p.ajouter("critere_arete", &critere_arete_); // XD_ADD_P floattant not_set
  p.ajouter("seuil_dvolume_residuel", &seuil_dvolume_residuel_); // XD_ADD_P floattant not_set
  p.ajouter("nb_iter_correction_volume",  &nb_iter_bary_volume_seul_); // XD_ADD_P entier not_set
  p.ajouter("nb_iter_remaillage", &nb_iter_remaillage_); // XD_ADD_P entier not_set
  p.ajouter("facteur_longueur_ideale", &facteur_longueur_ideale_); // XD_ADD_P floattant not_set
  p.ajouter("equilateral", &equilateral_); // XD_ADD_P entier not_set
  p.ajouter("lissage_courbure_coeff", &lissage_courbure_coeff_); // XD_ADD_P floattant not_set
  p.ajouter("lissage_courbure_iterations_systematique", &lissage_courbure_iterations_systematique_); // XD_ADD_P entier not_set
  p.ajouter("lissage_courbure_iterations_si_remaillage", &lissage_courbure_iterations_si_remaillage_); // XD_ADD_P entier not_set
  p.lire_avec_accolades_depuis(is);

  Cout << "Remaillage_FT_IJK::readOn : Les options lues sont : " << finl;
  p.print(Cout);

  if ((dt_remaillage_ > 0. ) && (facteur_longueur_ideale_<0.))
    {
      Cerr << "Erreur dans les parametres de Remaillage_FT_IJK."
           << " Il faut specifier la valeur de facteur_longueur_ideale si pas_remaillage>0."
           << " La valeur recommandee est : facteur_longueur_ideale 1.5 (en 3D)." << finl;
      Process::exit();
    }
  return is;
}




/* methode static duplique depuis Remaillage_FT --> a faire autrement ? */
static void SPA_choisir_sommets_remplacement(const Maillage_FT_IJK& maillage,
                                             const IntTab& tab_aretesMarquees,
                                             IntTab& sommets_remplacement)
{
  const int       nb_sommets       = maillage.nb_sommets();
  const IntTab&    facettes         = maillage.facettes();
  const ArrOfInt& sommet_PE_owner  = maillage.sommet_PE_owner();
  const ArrOfInt& sommet_num_owner = maillage.sommet_num_owner();

  sommets_remplacement.resize(nb_sommets, 2);
  // On initialise le tableau a "no_PE": on mettra ensuite
  // les processeurs d'accord entre eux avec un
  //   "collecter(Descripteur_FT::MIN_COLONNE1)"
  // Si tout le monde contient "no_PE", alors le noeud ne devra pas etre remplace.
  const int no_PE = Process::nproc();
  sommets_remplacement = no_PE;

  // Pour chaque sommet, y a-t-il un processeur qui veut le conserver ?
  // (le sommet sert a remplacer d'autres sommets)
  // 0 => non  sinon => oui
  ArrOfIntFT sommets_conserves(nb_sommets);
  sommets_conserves = 0;

  // Parcours de la liste des aretes a supprimer. Pour chacune, on choisit
  // un sommet a conserver et un sommet a supprimer (le plus petit en numerotation
  // globale pour que le choix soit identique pour les deux triangles adjacents).
  const int n = tab_aretesMarquees.dimension(0);
  const int nb_som_par_facette = facettes.dimension(1);
  int i;
  for (i = 0; i < n; i++)
    {
      const int fa7    = tab_aretesMarquees(i,0);
      const int isom   = tab_aretesMarquees(i,1);
      assert(isom >= 0 && isom < nb_som_par_facette);
      const int isom_s = (isom + 1 < nb_som_par_facette) ? isom + 1 : 0;
      const int som    = facettes(fa7, isom);
      const int som_s  = facettes(fa7, isom_s);
      const int bord   = maillage.sommet_ligne_contact(som);
      const int bord_s = maillage.sommet_ligne_contact(som_s);

      int somRempl = -1; // Le sommet de remplacement (num_owner)
      int somSupp  = -1; // Le sommet supprime (indice local)
      if (bord == bord_s)
        {
          //les deux sommets sont sur le bord ou les 2 sont internes
          //on compare alors leur numero global
          const int pe         = sommet_PE_owner[som];
          const int pe_s       = sommet_PE_owner[som_s];
          const int numOwner   = sommet_num_owner[som];
          const int numOwner_s = sommet_num_owner[som_s];
          if (FTd_compare_sommets_global(pe, numOwner, pe_s, numOwner_s) > 0)
            {
              //som_s>som : on garde som, on supprime som_s
              somRempl = som;
              somSupp = som_s;
            }
          else
            {
              //som_s<=som : on garde som_s, on supprime som
              somRempl = som_s;
              somSupp = som;
            }
        }
      else if (bord==1)
        {
          //som est de bord : on garde som, on supprimer som_s
          somRempl = som;
          somSupp = som_s;
        }
      else
        {
          //som_s est de bord : on garde som_s, on supprime som
          somRempl = som_s;
          somSupp = som;
        }
      // Remplacement du sommet si
      // * le sommet a supprimer ne l'a pas encore ete
      // * le sommet de remplacement n'a pas encore ete supprime
      // * le sommet a supprimer n'est pas utilise pour en remplacer un autre
      if (   sommets_remplacement(somRempl, 0) == no_PE
             && sommets_remplacement(somSupp, 0) == no_PE
             && sommets_conserves[somSupp] == 0)
        {
          const int pe  = sommet_PE_owner[somRempl];
          const int le_som = sommet_num_owner[somRempl];
          sommets_remplacement(somSupp, 0) = pe;
          sommets_remplacement(somSupp, 1) = le_som;
          sommets_conserves[somRempl] = 1;
        }
    }

  // Si un sommet est conserve sur un processeur, il ne doit etre supprime
  // sur aucun processeur. On fait la somme pour tous les processeurs
  // de "sommets_conserves" : si un sommet est conserve sur l'un, il devra
  // etre conserve sur tous.
  const Desc_Structure_FT& desc_sommets = maillage.desc_sommets();
  desc_sommets.collecter_espace_virtuel(sommets_conserves,
                                        MD_Vector_tools::EV_SOMME);
  desc_sommets.echange_espace_virtuel(sommets_conserves);

  // On met aussi tous les processeurs d'accord sur le sommet de remplacement.
  // Choix arbitraire de l'un des sommets de remplacement parmi ceux
  // proposes par les differents processeurs.
  desc_sommets.collecter_espace_virtuel(sommets_remplacement,
                                        MD_Vector_tools::EV_MINCOL1);
  desc_sommets.echange_espace_virtuel(sommets_remplacement);

  // Mise a jour de sommets_remplacement : -1 pour les sommets a conserver
  for (i = 0; i < nb_sommets; i++)
    {
      const int a_conserver = sommets_conserves[i];
      const int non_remplace = (sommets_remplacement(i,0) == no_PE);
      if (a_conserver || non_remplace)
        {
          sommets_remplacement(i,0) = -1;
          sommets_remplacement(i,1) = -1;
        }
    }
}

void Remaillage_FT_IJK::barycentrer_lisser_systematique_ijk(Maillage_FT_IJK& maillage,
                                                            ArrOfDouble& var_volume)
{
  static Stat_Counter_Id barycentre_lissage_sys_counter_ = statistiques().new_counter(2, "Remaillage interf: bary/lissage systematiques");
  FT_Field& Surfactant = maillage.Surfactant_facettes_non_const();
  ArrOfDouble surfactant_avant_remaillage ;
  if (!Surfactant.get_disable_surfactant())
    {
      var_volume.resize_array(maillage.nb_sommets());
      surfactant_avant_remaillage = Surfactant.check_conservation(maillage);
      maillage.set_barycentrage(true);
      //nettoyer_maillage(maillage);
      Surfactant.passer_variable_intensive(maillage);
    }
  statistiques().begin_count(barycentre_lissage_sys_counter_);
  regulariser_maillage(maillage,
                       var_volume,
                       relax_barycentrage_,
                       lissage_courbure_coeff_,
                       nb_iter_barycentrage_,
                       lissage_courbure_iterations_systematique_,
                       nb_iter_bary_volume_seul_,
                       seuil_dvolume_residuel_);
  supprimer_facettes_bord(maillage);
  if (!Surfactant.get_disable_surfactant())
    {
      Surfactant.passer_variable_extensive(maillage);
      maillage.set_barycentrage(false);
      ArrOfDouble surfactant_apres_remaillage = Surfactant.check_conservation(maillage);
      Surfactant.correction_conservation_globale(maillage,  surfactant_avant_remaillage,  surfactant_apres_remaillage);
    }
  statistiques().end_count(barycentre_lissage_sys_counter_);
}

void Remaillage_FT_IJK::barycentrer_lisser_apres_remaillage(Maillage_FT_IJK& maillage, ArrOfDouble& var_volume)
{
  static Stat_Counter_Id barycentre_lissage_apres_counter_ = statistiques().new_counter(2, "Remaillage local: bary/lissage apres remaillage");
  statistiques().begin_count(barycentre_lissage_apres_counter_);
  FT_Field& surfactant = maillage.Surfactant_facettes_non_const();
  if (!surfactant.get_disable_surfactant())
    {
      //ArrOfDouble surfactant_avant_remaillage = surfactant.check_conservation(maillage);
      maillage.set_barycentrage(true);
      //nettoyer_maillage(maillage);
      surfactant.passer_variable_intensive(maillage);
    }
  regulariser_maillage(maillage, var_volume,
                       relax_barycentrage_,
                       lissage_courbure_coeff_,
                       nb_iter_barycentrage_,
                       lissage_courbure_iterations_si_remaillage_,
                       nb_iter_bary_volume_seul_,
                       seuil_dvolume_residuel_);
  supprimer_facettes_bord(maillage);
  if (!surfactant.get_disable_surfactant())
    {
      surfactant.passer_variable_extensive(maillage);
      maillage.set_barycentrage(false);
      //ArrOfDouble surfactant_apres_remaillage = surfactant.check_conservation(maillage);
      //surfactant.correction_conservation_globale(maillage,  surfactant_avant_remaillage,  surfactant_apres_remaillage);
    }

  // Dans le doute, je laisse l'appel a nettoyer_maillage :
  nettoyer_maillage(maillage);
  statistiques().end_count(barycentre_lissage_apres_counter_);
}

// Surcharge de Remaillage_FT::diviser_grandes_aretes(Maillage_FT_Disc& maillage) const
// A la creation des facettes, il faut leur attribuer un numero de compo connexe dans compo_connex_facettes_.
int Remaillage_FT_IJK::diviser_grandes_aretes(Maillage_FT_IJK& maillage) const
{
  static Stat_Counter_Id sup_div_aretes_counter_ = statistiques().new_counter(2, "Remaillage local: suppressions / divisions aretes");
  FT_Field& surfactant = maillage.Surfactant_facettes_non_const();
  statistiques().begin_count(sup_div_aretes_counter_);
  static int compteur = 0;
  static int test_val = -1;
  if (!surfactant.get_disable_surfactant())
    {
      maillage.corriger_proprietaires_facettes();
      surfactant.echange_espace_virtuel(maillage);
    }

  Process::Journal()<<"Remaillage_FT_IJK::diviser_grandes_aretes "<<temps_<<"  nb_som="<<maillage.nb_sommets()<<finl;
  Process::Journal()<<" Compteur = " << compteur << finl;
  compteur++;
  if (compteur == test_val)
    {
      Process::Journal() << " STOP." << finl;
    }
  //  int res = 1;

  maillage.nettoyer_elements_virtuels();

  //tableaux de stockage
  IntTabFT tab_aretesMarquees;
  ArrOfIntFT tab_somD;
  DoubleTabFT tab_deplacement_somD;

  //on commence par marquer les grandes aretes
  marquer_aretes(maillage,
                 tab_aretesMarquees,
                 tab_somD,
                 tab_deplacement_somD,
                 1 /* marquage des aretes trop grandes */);
  // resultat =  tab_aretesMarquees : [ fa7 iarete pe som ]

  const int nb_aretes_divis = tab_aretesMarquees.dimension(0);


  int nb_facettes = maillage.nb_facettes();
  const int nb_facettes0 = nb_facettes;
  IntTab& facettes = maillage.facettes_;
  const int nb_som_par_facette = facettes.dimension(1);

  const ArrOfInt& sommet_num_owner = maillage.sommet_num_owner_;
  int fa7,iarete, isom,som,isom_s,som_s,isom_ss,som_ss, pe_somD,numOwner_somD,somD,somD_s,somD_ss;

  const int dimension3 = (dimension==3);
  const int nb_aretes_par_facette = (dimension3)?3:1;
  //tableau stockant le sommet servant a decouper l'arete (ou -1 si arete non divisee)
  IntTabFT tab_fa7Divis(nb_facettes,nb_aretes_par_facette);
  for (fa7=0 ; fa7<nb_facettes ; fa7++)
    {
      for (isom=0 ; isom<nb_aretes_par_facette ; isom++)
        {
          tab_fa7Divis(fa7,isom) = -1;
        }
    }

  //on va balayer les aretes a memoriser l'ensemble des aretes a scinder
  for (iarete=0 ; iarete<nb_aretes_divis ; iarete++)
    {
      fa7 = tab_aretesMarquees(iarete,0);
      isom = tab_aretesMarquees(iarete,1);
      isom_s = (isom+1)%nb_som_par_facette;
      //sommet qui va rester dans fa7
      som = facettes(fa7,isom);
      //sommet qui va aller dans une nouvelle facette
      som_s = facettes(fa7,isom_s);
      //sommet a inserer dans l'arete
      pe_somD = tab_aretesMarquees(iarete,2);
      numOwner_somD = tab_aretesMarquees(iarete,3);

      if (pe_somD!=me())
        {
          //je ne connais pas l'indice du sommet a inserer dans me()
          maillage.convertir_numero_distant_local(maillage.desc_sommets(),sommet_num_owner,numOwner_somD,pe_somD,somD);
          assert(somD >= 0);
        }
      else
        {
          //je suis le proprietaire de somD : je connais donc le bon indice du sommet a inserer
          somD = numOwner_somD;
        }

      tab_fa7Divis(fa7,isom) = somD;
    }

  // Specifique IJK :
  ArrOfInt& compo_connexe_facettes = maillage.compo_connexe_facettes_non_const();

  //on va ensuite balayer les facettes et les scinder
  //la configuration depend du nb d'aretes a scinder par facette
  int nb_areteScinder, isom0=-1,isom1=-1;


  if (!surfactant.get_disable_surfactant())
    {
      for (fa7=0 ; fa7<nb_facettes0 ; fa7++)
        {
          //on compte le nombre d'aretes a scinder
          nb_areteScinder = 0;
          for (isom=0 ; isom<nb_aretes_par_facette ; isom++)
            {
              if (tab_fa7Divis(fa7,isom)>=0)
                {
                  if (nb_areteScinder==0)
                    {
                      isom0 = isom;
                    }
                  else
                    {
                      isom1 = isom;
                    }
                  nb_areteScinder++;
                }
            }
          if (nb_areteScinder==1)
            {
              double concentration_surfactant_avant_decoupage = surfactant[fa7];
              //s'il n'y a qu'une arete a scinder
              somD = tab_fa7Divis(fa7,isom0);
              //on modifie la facettes d'origine
              isom_s = (isom0+1)%nb_som_par_facette;
              som_s = facettes(fa7,isom_s);
              facettes(fa7,isom_s) = somD;
              surfactant[fa7]=concentration_surfactant_avant_decoupage;
              //on cree une nouvelle facette
              if (nb_facettes>=facettes.dimension(0))
                {
                  //Redimensionnement
                  facettes.resize(nb_facettes+10,facettes.dimension(1));
                  compo_connexe_facettes.resize_array(nb_facettes+10);
                  surfactant.resize_array(nb_facettes+10);
                }
              isom_ss = (isom_s+1)%nb_som_par_facette;
              som_ss = facettes(fa7,isom_ss);
              compo_connexe_facettes[nb_facettes] = compo_connexe_facettes[fa7];
              surfactant[nb_facettes]=concentration_surfactant_avant_decoupage;
              if (dimension==2)
                {
                  facettes(nb_facettes,isom0) = somD;
                  facettes(nb_facettes,isom_s) = som_s;
                }
              else
                {
                  facettes(nb_facettes,0) = som_ss;
                  facettes(nb_facettes,1) = somD; //sommet insere
                  facettes(nb_facettes,2) = som_s;
                }
              nb_facettes++;
            }
          else if (nb_areteScinder==2)
            {
              //si 2 aretes sont a scinder
              //on cree deux nouvelles facettes
              double concentration_surfactant_avant_decoupage  = surfactant[fa7];
              if (nb_facettes+2>=facettes.dimension(0))
                {
                  //Redimensionnement
                  facettes.resize(nb_facettes+12,facettes.dimension(1));
                  compo_connexe_facettes.resize_array(nb_facettes+12);
                  surfactant.resize_array(nb_facettes+12);
                }
              //on positionne isom et isom_s tq elles soient les aretes a scinder
              if (isom1==(isom0+1)%nb_som_par_facette)
                {
                  isom = isom0;
                  isom_s = isom1;
                }
              else
                {
                  isom = isom1;
                  isom_s = isom0;
                }
              isom_ss = ((isom_s+1)%nb_som_par_facette);
              som = facettes(fa7,isom);
              som_s = facettes(fa7,isom_s);
              som_ss = facettes(fa7,isom_ss);
              somD = tab_fa7Divis(fa7,isom);
              somD_s = tab_fa7Divis(fa7,isom_s);
              somD_ss = tab_fa7Divis(fa7,isom_ss);
              //on modifie la facette existante
              facettes(fa7,0) = som;
              facettes(fa7,1) = somD;
              facettes(fa7,2) = som_ss;
              surfactant[fa7]=concentration_surfactant_avant_decoupage;
              //premiere nouvelle facette
              compo_connexe_facettes[nb_facettes] = compo_connexe_facettes[fa7];
              facettes(nb_facettes,0) = somD;
              facettes(nb_facettes,1) = som_s;
              facettes(nb_facettes,2) = somD_s;
              surfactant[nb_facettes]=concentration_surfactant_avant_decoupage;
              nb_facettes++;
              //seconde nouvelle facette
              compo_connexe_facettes[nb_facettes] = compo_connexe_facettes[fa7];
              facettes(nb_facettes,0) = somD;
              facettes(nb_facettes,1) = somD_s;
              facettes(nb_facettes,2) = som_ss;
              surfactant[nb_facettes]=concentration_surfactant_avant_decoupage;
              nb_facettes++;
            }
          else if (nb_areteScinder==3)
            {
              //si toutes les aretes sont a scinder
              //on cree trois nouvelles facettes
              double concentration_surfactant_avant_decoupage = surfactant[fa7];
              if (nb_facettes+3>=facettes.dimension(0))
                {
                  //Redimensionnement
                  facettes.resize(nb_facettes+13,facettes.dimension(1));
                  compo_connexe_facettes.resize_array(nb_facettes+13);
                  surfactant.resize_array(nb_facettes+13);
                }
              isom = 0;
              isom_s = 1;
              isom_ss = 2;
              som = facettes(fa7,isom);
              som_s = facettes(fa7,isom_s);
              som_ss = facettes(fa7,isom_ss);
              somD = tab_fa7Divis(fa7,isom);
              somD_s = tab_fa7Divis(fa7,isom_s);
              somD_ss = tab_fa7Divis(fa7,isom_ss);
              //on modifie la facette existante
              facettes(fa7,0) = som;
              facettes(fa7,1) = somD;
              facettes(fa7,2) = somD_ss;
              surfactant[fa7]=concentration_surfactant_avant_decoupage;
              //premiere nouvelle facette
              compo_connexe_facettes[nb_facettes] = compo_connexe_facettes[fa7];
              facettes(nb_facettes,0) = somD;
              facettes(nb_facettes,1) = som_s;
              facettes(nb_facettes,2) = somD_s;
              surfactant[nb_facettes]=concentration_surfactant_avant_decoupage;
              nb_facettes++;
              //seconde nouvelle facette
              compo_connexe_facettes[nb_facettes] = compo_connexe_facettes[fa7];
              facettes(nb_facettes,0) = somD_s;
              facettes(nb_facettes,1) = som_ss;
              facettes(nb_facettes,2) = somD_ss;
              surfactant[nb_facettes]=concentration_surfactant_avant_decoupage;
              nb_facettes++;
              //troisieme nouvelle facette
              compo_connexe_facettes[nb_facettes] = compo_connexe_facettes[fa7];
              facettes(nb_facettes,0) = somD;
              facettes(nb_facettes,1) = somD_s;
              facettes(nb_facettes,2) = somD_ss;
              surfactant[nb_facettes]=concentration_surfactant_avant_decoupage;
              nb_facettes++;
            }
        }
    }
  else
    {
      for (fa7=0 ; fa7<nb_facettes0 ; fa7++)
        {
          //on compte le nombre d'aretes a scinder
          nb_areteScinder = 0;
          for (isom=0 ; isom<nb_aretes_par_facette ; isom++)
            {
              if (tab_fa7Divis(fa7,isom)>=0)
                {
                  if (nb_areteScinder==0)
                    {
                      isom0 = isom;
                    }
                  else
                    {
                      isom1 = isom;
                    }
                  nb_areteScinder++;
                }
            }
          if (nb_areteScinder==1)
            {
              //s'il n'y a qu'une arete a scinder
              somD = tab_fa7Divis(fa7,isom0);
              //on modifie la facettes d'origine
              isom_s = (isom0+1)%nb_som_par_facette;
              som_s = facettes(fa7,isom_s);
              facettes(fa7,isom_s) = somD;
              //on cree une nouvelle facette
              if (nb_facettes>=facettes.dimension(0))
                {
                  //Redimensionnement
                  facettes.resize(nb_facettes+10,facettes.dimension(1));
                  compo_connexe_facettes.resize_array(nb_facettes+10);
                }
              isom_ss = (isom_s+1)%nb_som_par_facette;
              som_ss = facettes(fa7,isom_ss);
              compo_connexe_facettes[nb_facettes] = compo_connexe_facettes[fa7];
              if (dimension==2)
                {
                  facettes(nb_facettes,isom0) = somD;
                  facettes(nb_facettes,isom_s) = som_s;
                }
              else
                {
                  facettes(nb_facettes,0) = som_ss;
                  facettes(nb_facettes,1) = somD; //sommet insere
                  facettes(nb_facettes,2) = som_s;
                }
              nb_facettes++;
            }
          else if (nb_areteScinder==2)
            {
              //si 2 aretes sont a scinder
              //on cree deux nouvelles facettes
              if (nb_facettes+2>=facettes.dimension(0))
                {
                  //Redimensionnement
                  facettes.resize(nb_facettes+12,facettes.dimension(1));
                  compo_connexe_facettes.resize_array(nb_facettes+12);
                }
              //on positionne isom et isom_s tq elles soient les aretes a scinder
              if (isom1==(isom0+1)%nb_som_par_facette)
                {
                  isom = isom0;
                  isom_s = isom1;
                }
              else
                {
                  isom = isom1;
                  isom_s = isom0;
                }
              isom_ss = ((isom_s+1)%nb_som_par_facette);
              som = facettes(fa7,isom);
              som_s = facettes(fa7,isom_s);
              som_ss = facettes(fa7,isom_ss);
              somD = tab_fa7Divis(fa7,isom);
              somD_s = tab_fa7Divis(fa7,isom_s);
              somD_ss = tab_fa7Divis(fa7,isom_ss);
              //on modifie la facette existante
              facettes(fa7,0) = som;
              facettes(fa7,1) = somD;
              facettes(fa7,2) = som_ss;
              //premiere nouvelle facette
              compo_connexe_facettes[nb_facettes] = compo_connexe_facettes[fa7];
              facettes(nb_facettes,0) = somD;
              facettes(nb_facettes,1) = som_s;
              facettes(nb_facettes,2) = somD_s;
              nb_facettes++;
              //seconde nouvelle facette
              compo_connexe_facettes[nb_facettes] = compo_connexe_facettes[fa7];
              facettes(nb_facettes,0) = somD;
              facettes(nb_facettes,1) = somD_s;
              facettes(nb_facettes,2) = som_ss;
              nb_facettes++;
            }
          else if (nb_areteScinder==3)
            {
              //si toutes les aretes sont a scinder
              //on cree trois nouvelles facettes
              if (nb_facettes+3>=facettes.dimension(0))
                {
                  //Redimensionnement
                  facettes.resize(nb_facettes+13,facettes.dimension(1));
                  compo_connexe_facettes.resize_array(nb_facettes+13);
                }
              isom = 0;
              isom_s = 1;
              isom_ss = 2;
              som = facettes(fa7,isom);
              som_s = facettes(fa7,isom_s);
              som_ss = facettes(fa7,isom_ss);
              somD = tab_fa7Divis(fa7,isom);
              somD_s = tab_fa7Divis(fa7,isom_s);
              somD_ss = tab_fa7Divis(fa7,isom_ss);
              //on modifie la facette existante
              facettes(fa7,0) = som;
              facettes(fa7,1) = somD;
              facettes(fa7,2) = somD_ss;
              //premiere nouvelle facette
              compo_connexe_facettes[nb_facettes] = compo_connexe_facettes[fa7];
              facettes(nb_facettes,0) = somD;
              facettes(nb_facettes,1) = som_s;
              facettes(nb_facettes,2) = somD_s;
              nb_facettes++;
              //seconde nouvelle facette
              compo_connexe_facettes[nb_facettes] = compo_connexe_facettes[fa7];
              facettes(nb_facettes,0) = somD_s;
              facettes(nb_facettes,1) = som_ss;
              facettes(nb_facettes,2) = somD_ss;
              nb_facettes++;
              //troisieme nouvelle facette
              compo_connexe_facettes[nb_facettes] = compo_connexe_facettes[fa7];
              facettes(nb_facettes,0) = somD;
              facettes(nb_facettes,1) = somD_s;
              facettes(nb_facettes,2) = somD_ss;
              nb_facettes++;
            }
        }
    }
  //Redimensionnement
  facettes.resize(nb_facettes,facettes.dimension(1));
  compo_connexe_facettes.resize_array(nb_facettes);
  if (!surfactant.get_disable_surfactant())
    surfactant.resize_array(nb_facettes);
  maillage.desc_facettes_.calcul_schema_comm(nb_facettes);
  maillage.corriger_proprietaires_facettes();

  ArrOfIntFT liste_sommets_sortis;
  ArrOfIntFT numero_face_sortie;
  if (!surfactant.get_disable_surfactant())
    {
      surfactant.echange_espace_virtuel(maillage);
    }
  maillage.deplacer_sommets(tab_somD,tab_deplacement_somD,liste_sommets_sortis,numero_face_sortie);
  maillage.corriger_proprietaires_facettes();
  if (!surfactant.get_disable_surfactant())
    surfactant.echange_espace_virtuel(maillage);

  int nb_aretes_divis_tot = Process::mp_sum(nb_aretes_divis);
  Process::Journal()<<"FIN Remaillage_FT::diviser_grandes_aretes "<<temps_<<"  nb_som="<<maillage.nb_sommets()
                    <<"  nb_aretes_divisees_on_proc="<< nb_aretes_divis;
  Process::Journal()<< "  nb_aretes_divisees_tot="<< nb_aretes_divis_tot<<finl;
  maillage.maillage_modifie(Maillage_FT_Disc::MINIMAL);

  statistiques().end_count(sup_div_aretes_counter_);

  return nb_aretes_divis_tot;
  //  return res;
}


/* Surcharge de supprimer_petites_aretes de Remaillage_FT */
// pour gestion des Surfactants IJK
int Remaillage_FT_IJK::supprimer_petites_aretes(Maillage_FT_IJK& maillage,
                                                ArrOfDouble& varVolume) const
{
  FT_Field& surfactant = maillage.Surfactant_facettes_non_const();
  if (surfactant.get_disable_surfactant())
    {
      return Remaillage_FT::supprimer_petites_aretes(maillage,varVolume);
    }

  static const Stat_Counter_Id stat_counter = statistiques().new_counter(3, "Supprimer_petites_aretes");
  statistiques().begin_count(stat_counter);
  maillage.check_mesh();
  int nb_sommets_supprimes_tot = 0;
  int nb_sommets_supprimes = 0;
  do
    {
      // Pour chaque sommet a remplacer:
      // colonne0: indice local du sommet a remplacer
      // colonne1: indice local du sommet de remplacement
      IntTabFT remplacement_ilocal(0,2);

      {
        // ******************************************************
        // Recherche des aretes trop petites
        IntTabFT tab_aretesMarquees;
        {
          ArrOfIntFT tab_somD;             // Inutilise
          DoubleTabFT tab_deplacement_somD; // Inutilise
          marquer_aretes(maillage,
                         tab_aretesMarquees,
                         tab_somD,
                         tab_deplacement_somD,
                         -1 /* marquage des aretes trop petites */);
        }

        // ******************************************************
        // On supprime les aretes en remplacant une extremite par l'autre
        // (genere des triangles dont deux sommets sont confondus)
        // Determination des sommets a remplacer et des sommets de remplacement
        // en numerotation globale.
        IntTabFT sommets_remplacement;
        SPA_choisir_sommets_remplacement(maillage,
                                         tab_aretesMarquees,
                                         sommets_remplacement);

        // ******************************************************
        // Certains sommets de remplacement n'existent pas encore sur le processeur local.
        // Creation de la liste des sommets de remplacement virtuels dont on aura besoin
        ArrOfIntFT request_sommets_pe;
        ArrOfIntFT request_sommets_num;
        int i;
        const int nb_sommets = maillage.nb_sommets();
        const int moi = Process::me();
        assert(nb_sommets == sommets_remplacement.dimension(0));
        for (i = 0; i < nb_sommets; i++)
          {
            const int pe = sommets_remplacement(i,0);
            const int som= sommets_remplacement(i,1);
            if (pe >= 0)
              {
                int j;
                if (pe != moi)   // Le sommet de remplacement est virtuel
                  {
                    request_sommets_pe.append_array(pe);
                    request_sommets_num.append_array(som);
                    j = -1;
                  }
                else
                  {
                    j = som;
                  }
                remplacement_ilocal.append_line(i, j);
              }
          }
        // ******************************************************
        // Creation des sommets virtuels qui remplaceront les sommets supprimes
        // (certains existent peut-etre deja, ils ne seront pas recrees)

        maillage.creer_sommets_virtuels_numowner(request_sommets_pe,
                                                 request_sommets_num);
        // On a cree de nouveaux sommets virtuels. Mise a jour de varVolume:
        // Inutile de mettre a jour l'espace virtuel, l'echange sera fait
        // apres le calcul de dvolume
        {
          const int old_size = varVolume.size_array();
          const int new_size = maillage.nb_sommets();
          varVolume.resize_array(new_size);
          int ii;
          for (ii = old_size; ii < new_size; ii++)
            varVolume[ii] = 0.;
        }
        // ******************************************************
        // Calcul de l'indice local des sommets de remplacement
        ArrOfIntFT request_sommets_ilocal;

        maillage.convertir_numero_distant_local(maillage.desc_sommets(),
                                                maillage.sommet_num_owner(),
                                                request_sommets_num,
                                                request_sommets_pe,
                                                request_sommets_ilocal);
        int j = 0;
        const int n_rempl = remplacement_ilocal.dimension(0);
        for (i = 0; i < n_rempl; i++)
          {
            int som_new = remplacement_ilocal(i,1);
            if (som_new < 0)   // Le sommet de remplacement est virtuel ?
              {
                som_new = request_sommets_ilocal[j];
                remplacement_ilocal(i,1) = som_new;
                j++;
              }
          }
      }
      // ******************************************************
      // Determination de la variation de volume liee a la suppression des aretes
      // On calcule la variation de volume correspondant au deplacement des
      // noeuds remplaces vers les noeuds de remplacement.
      {
        const DoubleTab& sommets = maillage.sommets();
        // Nombre de sommets a remplacer
        const int n_rempl = remplacement_ilocal.dimension(0);
        DoubleTabFT position_finale(sommets); // Copie du tableau.
        int i;
        const int dim = Objet_U::dimension;
        for (i = 0; i < n_rempl; i++)
          {
            const int som_old = remplacement_ilocal(i,0);
            const int som_new = remplacement_ilocal(i,1);
            int j;
            for (j = 0; j < dim; j++)
              position_finale(som_old, j) = position_finale(som_new, j);
          }
        {
          ArrOfDoubleFT dvolume;
          calculer_variation_volume(maillage,
                                    position_finale,
                                    dvolume);

          // Avant d'affecter la variation de volume des noeuds remplaces aux noeuds
          // de remplacement: certains sommets de remplacement sont virtuels,
          // il faut donc "collecter" les contributions a la fin. Donc il faut
          // annuler la valeur de varVolume pour les sommets virtuels avant
          // de commencer.
          const int n = maillage.nb_sommets();
          for (i = 0; i < n; i++)
            {
              if (maillage.sommet_virtuel(i))
                varVolume[i] = 0.;
              else
                varVolume[i] += dvolume[i];
            }
        }
        for (i = 0; i < n_rempl; i++)
          {
            const int som_old = remplacement_ilocal(i,0);
            const int som_new = remplacement_ilocal(i,1);
            const double dv = varVolume[som_old];
            varVolume[som_new] += dv;
            varVolume[som_old] = 0;
          }
        const Desc_Structure_FT& desc = maillage.desc_sommets();
        desc.collecter_espace_virtuel(varVolume, MD_Vector_tools::EV_SOMME);
        desc.echange_espace_virtuel(varVolume);
      }

      // ******************************************************
      // Remplacement des sommets dans le tableau des facettes,
      // On rend invalides les facettes qui disparaissent (deux sommets confondus)


      {
        const int nb_som = maillage.nb_sommets();
        // Creation d'une table de remplacement de sommets:
        ArrOfIntFT table_old_new(nb_som);
        table_old_new = -1;
        const int nb_rempl = remplacement_ilocal.dimension(0);
        for (int i = 0; i < nb_rempl; i++)
          {
            const int som_old = remplacement_ilocal(i,0);
            const int som_new = remplacement_ilocal(i,1);
            table_old_new[som_old] = som_new;
          }
        IntTab& facettes = maillage.facettes_;
        const DoubleTab& sommets = maillage.sommets_;
        //const DoubleTab& sommets = maillage.sommets();
        const int nb_facettes = facettes.dimension(0);
        const int nb_som_par_facette = facettes.dimension(1);


        // doit mettre a jour sommet


        // int nb_proc = splitting->get_nprocessor_per_direction(0)*splitting->get_nprocessor_per_direction(1)*splitting->get_nprocessor_per_direction(2);

        //////////////////////// GESTION DU CHAMP DE SURFACTANT QUAND ON REMPLACE LES SOMMETS DES FACETTES //////////////////
        /* Avant de changer le sommet, on stocke les triangles (facettes) initiaux avant changement (i.e. les triangles dont l'un des trois sommets va changer)
         * Apres la modification du sommet, on stocke les triangles modifiés (i.e. les triangles dont l'un des trois sommets va changer)
         * On calcule ensuite toute les intersections entre les triangles initiaux et les triangles finaux
         * Ces intersections permettent de reconstruire le champ de surfactant apres transport des sommets de maniere conservative et sans diffusion numerique
         * Si plusieurs arretes successives sont supprimees, le nombre de triangles impliqué sur un seul sommet modifie peut être très grand
         * La parallèlisation est alors très compliquee. On peut avoir besoin de voisins non direct.
         * C'est pourquoi, en attendant une indexation déterministe des facettes et des sommets, on est obligés de partager les tableaux complets sur l'ensemble des procs
         */

        // on stocke les coord des sommets reels qui bougent

        // premier partage entre proc : coordonnees new_sommet sur les proc compo
        //surfactant.dimensionner_remaillage_FT_Field(maillage, table_old_new);
        //DoubleTab sommet_bouge = surfactant.get_sommet_bouge();
        // construction des tableau avant/apres deplacement
        const ArrOfDouble& Sfa7 = maillage.get_surface_facettes();
        double longueur_cara_fa7 = 0.;
        int n = 0 ;
        for (int fa = 0; fa < Sfa7.size_array(); fa++)
          {
            longueur_cara_fa7 += Sfa7(fa);
            n++;
          }
        longueur_cara_fa7= Process::mp_sum(longueur_cara_fa7);
        n= Process::mp_sum(n);
        if (n>0)
          {
            longueur_cara_fa7 /= n ;
            longueur_cara_fa7 = std::sqrt(longueur_cara_fa7);
          }

        surfactant.set_tolerance_point_identique(longueur_cara_fa7);
        surfactant.echange_espace_virtuel(maillage);

        OBS_PTR(IJK_Splitting) splitting = maillage.ref_splitting();
        DoubleTab liste_sommets_avant_deplacement(0, 3);
        DoubleTab liste_sommets_apres_deplacement(0, 3);
        ArrOfInt compo_connexe_sommets_deplace(0);
        ArrOfInt& compo_connexe_facettes = maillage.compo_connexe_facettes_non_const();
        surfactant.champ_sommet_from_facettes(compo_connexe_facettes, maillage);
        ArrOfInt compo_connexe_sommet = surfactant.get_compo_connexe_sommets();
        int index = 0 ;
        for (int i = 0; i < nb_rempl; i++)
          {
            const int som_old = remplacement_ilocal(i,0);
            const int som_new = remplacement_ilocal(i,1);
            if (som_new >= 0)
              {
                liste_sommets_avant_deplacement.resize(index+1,3);
                liste_sommets_apres_deplacement.resize(index+1,3);
                compo_connexe_sommets_deplace.resize(index+1);
                liste_sommets_avant_deplacement(index, 0) = sommets(som_old,0);
                liste_sommets_avant_deplacement(index, 1) = sommets(som_old,1);
                liste_sommets_avant_deplacement(index, 2) = sommets(som_old,2);
                liste_sommets_apres_deplacement(index, 0) = sommets(som_new,0);
                liste_sommets_apres_deplacement(index, 1) = sommets(som_new,1);
                liste_sommets_apres_deplacement(index, 2) = sommets(som_new,2);
                compo_connexe_sommets_deplace(index) = compo_connexe_sommet(som_old);
                //std::cout << " som avant = " << liste_sommets_avant_deplacement(index, 0) << liste_sommets_avant_deplacement(index, 1) << liste_sommets_avant_deplacement(index, 2);
                //std::cout << " som apres = " << liste_sommets_apres_deplacement(index, 0) << liste_sommets_apres_deplacement(index, 1) << liste_sommets_apres_deplacement(index, 2) << std::endl;
                index++;
              }
          }


        surfactant.completer_compo_connexe_partielle(maillage, splitting, liste_sommets_apres_deplacement, liste_sommets_avant_deplacement, compo_connexe_sommets_deplace);
        DoubleTab& facettes_sommets_full_compo = surfactant.get_facettes_sommets_full_compo_non_const();
        DoubleTab& liste_sommets_et_deplacements_full_compo = surfactant.get_liste_sommets_et_deplacements_non_const();
        ArrOfInt sorted_index = surfactant.get_sorted_index();

        int nb_facettes_compo_complete = facettes_sommets_full_compo.dimension(0);
        int nbsom_compo_complete = liste_sommets_et_deplacements_full_compo.dimension(0);


        /*Process::barrier();
        int nb_proc = splitting->get_nprocessor_per_direction(0)*splitting->get_nprocessor_per_direction(1)*splitting->get_nprocessor_per_direction(2);
        int ny = splitting->get_nprocessor_per_direction(1);
        int nz = splitting->get_nprocessor_per_direction(2);
        int x = splitting->get_local_slice_index(0);
        int y = splitting->get_local_slice_index(1);
        int z = splitting->get_local_slice_index(2);
        int proc_index = z + y * nz + x * ny * nz ;
        int index_global = 0;
        for (int proc = 0; proc < nb_proc; proc++)
          {
            index_global = Process::mp_max(index_global);
            Process::barrier();
            if (proc == proc_index)
              {*/



        for (int indice_sommet_desordre = 0; indice_sommet_desordre < nbsom_compo_complete; indice_sommet_desordre++)
          {
            int indice_sommet=sorted_index(indice_sommet_desordre);
            Vecteur3 pos(liste_sommets_et_deplacements_full_compo(indice_sommet,0),liste_sommets_et_deplacements_full_compo(indice_sommet,1),liste_sommets_et_deplacements_full_compo(indice_sommet,2));
            Vecteur3 pos_apres_dep(liste_sommets_et_deplacements_full_compo(indice_sommet,3),liste_sommets_et_deplacements_full_compo(indice_sommet,4),liste_sommets_et_deplacements_full_compo(indice_sommet,5));
            surfactant.reinit_remeshing_table();
            for (int fa = 0; fa < nb_facettes_compo_complete; fa++)
              {
                /*bool facette_virtuelle_locale = false;
                if(facettes_sommets_full_compo(index, 11)==1.)
                  {
                    facette_virtuelle_locale = true;
                  }*/

                for (int somfa7 = 0; somfa7 < nb_som_par_facette; somfa7++)
                  {
                    double sx = facettes_sommets_full_compo(fa, 3*somfa7+0);
                    double sy = facettes_sommets_full_compo(fa, 3*somfa7+1);
                    double sz = facettes_sommets_full_compo(fa, 3*somfa7+2);
                    Point3D p1 = {sx,sy,sz};
                    Point3D p2 = {pos[0],pos[1],pos[2]};
                    if (p1==p2)
                      {
                        sx = facettes_sommets_full_compo(fa, 0);
                        sy = facettes_sommets_full_compo(fa, 1);
                        sz = facettes_sommets_full_compo(fa, 2);
                        Point3D s0 = {sx,sy,sz};
                        sx = facettes_sommets_full_compo(fa, 3);
                        sy = facettes_sommets_full_compo(fa, 4);
                        sz = facettes_sommets_full_compo(fa, 5);
                        Point3D s1 = {sx,sy,sz};
                        sx = facettes_sommets_full_compo(fa, 6);
                        sy = facettes_sommets_full_compo(fa, 7);
                        sz = facettes_sommets_full_compo(fa, 8);
                        Point3D s2 = {sx,sy,sz};

                        if (!(s0 == s1 or s1 == s2 or s0 == s2))
                          {
                            // on sauvegarde le triangle avant deplacement
                            // on deplace le sommet du triangle
                            // on sauvegarde le triangle apres deplacement ssi il nest pas de surface nulle
                            //if(!facette_virtuelle_locale)
                            surfactant.sauvegarder_triangle(maillage, fa, 0);

                            Point3D n_apres_dep = surfactant.calculer_normale_apres_deplacement(fa, somfa7, pos_apres_dep);

                            // on bouge le sommet
                            for (int dir = 0; dir < 3; dir++)
                              facettes_sommets_full_compo(fa, 3*somfa7+dir)=pos_apres_dep[dir];
                            // on met a jour la normale au sommet
                            facettes_sommets_full_compo(fa, 12)=n_apres_dep.x;
                            facettes_sommets_full_compo(fa, 13)=n_apres_dep.y;
                            facettes_sommets_full_compo(fa, 14)=n_apres_dep.z;


                            const int indice_fa_locale = int(facettes_sommets_full_compo(fa, 10));
                            if(indice_fa_locale!=-1)
                              {
                                // alors la fa7 est reelle ou virtuelle locale sur le proc et il faut la modifier
                                const int i_sommet = facettes(indice_fa_locale,somfa7);
                                const int new_sommet = table_old_new[i_sommet];
                                facettes(indice_fa_locale,somfa7) = new_sommet;
                              }

                            sx = facettes_sommets_full_compo(fa, 0);
                            sy = facettes_sommets_full_compo(fa, 1);
                            sz = facettes_sommets_full_compo(fa, 2);
                            s0 = {sx,sy,sz};
                            sx = facettes_sommets_full_compo(fa, 3);
                            sy = facettes_sommets_full_compo(fa, 4);
                            sz = facettes_sommets_full_compo(fa, 5);
                            s1 = {sx,sy,sz};
                            sx = facettes_sommets_full_compo(fa, 6);
                            sy = facettes_sommets_full_compo(fa, 7);
                            sz = facettes_sommets_full_compo(fa, 8);
                            s2 = {sx,sy,sz};
                            if (!(s0 == s1 or s1 == s2 or s0 == s2))
                              {
                                //if(!facette_virtuelle_locale)
                                surfactant.sauvegarder_triangle(maillage, fa, 1);
                              }
                            else
                              {
                                // la fa7 doit devenir invalide
                                // on rempli avec des valeurs bidons
                                for (int somfa7bis = 0; somfa7bis < nb_som_par_facette; somfa7bis++)
                                  {
                                    for (int dir = 0; dir < 3; dir++)
                                      facettes_sommets_full_compo(fa, 3*somfa7bis+dir)=-123.;
                                  }
                                facettes_sommets_full_compo(fa, 9)=-123.;
                              }
                            break;
                          }
                      }
                  }
              }
            surfactant.remailler_FT_Field(maillage);
          }

        //}
        //}
        surfactant.update_FT_Field_local_from_full_compo(maillage);
        surfactant.echange_espace_virtuel(maillage);


        for (int i = 0; i < nb_facettes; i++)
          {
            // En dimension2, une facette est invalide si les deux sommets
            // sont identiques. En dimension 3, si deux sommets sont confondus,
            // il faut invalider la facette (sommet0 == sommet1)
            if (nb_som_par_facette == 3)
              {
                // Facette modifiee, on teste si elle disparait
                const int s0 = facettes(i,0);
                const int s1 = facettes(i,1);
                const int s2 = facettes(i,2);
                if (s0 == s2 || s1 == s2)
                  {
                    facettes(i,1) = s0;
                    facettes(i,2) = s0;
                  }
              }
          }
      }


      nb_sommets_supprimes = remplacement_ilocal.dimension(0);
      nb_sommets_supprimes = Process::mp_sum(nb_sommets_supprimes);
      nb_sommets_supprimes_tot += nb_sommets_supprimes;
      maillage.corriger_proprietaires_facettes();
      surfactant.echange_espace_virtuel(maillage);
    }
  while (nb_sommets_supprimes > 0);

  //maillage.corriger_proprietaires_facettes();
  //surfactant.echange_espace_virtuel(maillage);
  // corriger_proprietaires_facettes() cree eventuellement des sommets virtuels.
  // => Mise a jour de varVolume
  varVolume.resize_array(maillage.nb_sommets());
  maillage.desc_sommets().echange_espace_virtuel(varVolume);

  maillage.maillage_modifie(Maillage_FT_Disc::MINIMAL);
  statistiques().end_count(stat_counter);
  return nb_sommets_supprimes_tot;
}

void Remaillage_FT_IJK::remaillage_local_interface(double temps, Maillage_FT_IJK& maillage)
{
  static Stat_Counter_Id remaillage_loc_interf_counter_ = statistiques().new_counter(2, "Remaillage interf: remaillage local");
  FT_Field& surfactant = maillage.Surfactant_facettes_non_const();
  statistiques().begin_count(remaillage_loc_interf_counter_);
  temps_dernier_remaillage_ = temps_dernier_lissage_ = temps_ = temps;
  ArrOfDouble surfactant_avant_remaillage ;
  if (!surfactant.get_disable_surfactant())
    {
      surfactant_avant_remaillage = surfactant.check_conservation(maillage);
    }

  maillage.nettoyer_elements_virtuels();
  maillage.check_mesh();
  //boucle sur les remaillages
  int iter;
  ArrOfDoubleFT varVolume;

  for (iter = 0; iter < nb_iter_remaillage_; iter++)
    {
      if (!surfactant.get_disable_surfactant())
        {
          maillage.corriger_proprietaires_facettes();
          surfactant.echange_espace_virtuel(maillage);
          surfactant.passer_variable_intensive(maillage);
          //maillage.nettoyer_elements_virtuels();
          maillage.check_mesh();
        }

      const int nb_sommets = maillage.nb_sommets();
      varVolume.resize_array(nb_sommets);
      varVolume = 0.;
      variation_volume_ = 0.;
      // n = nombre de sommets supprimes
      int n = supprimer_petites_aretes(maillage,varVolume);
      if (!surfactant.get_disable_surfactant())
        {
          supprimer_doublons_facettes(maillage);
          //nettoyer_maillage(maillage);
          surfactant.passer_variable_extensive(maillage);
          maillage.corriger_proprietaires_facettes();
          surfactant.echange_espace_virtuel(maillage);
        }
      if (Comm_Group::check_enabled()) maillage.check_mesh();
      if (n > 0)
        {
          if (!surfactant.get_disable_surfactant())
            varVolume.resize_array(maillage.nb_sommets());

          supprimer_doublons_facettes(maillage);
          if (Comm_Group::check_enabled()) maillage.check_mesh();
          // On a supprime les petites aretes en deplacant un noeud sur
          //  un autre, la variation de volume engendree a ete mise dans varVolume:
          // On barycentre et on lisse, notamment pour recuperer cette variation.
          if(!surfactant.get_only_remaillage())
            barycentrer_lisser_apres_remaillage(maillage, varVolume);
          if (Comm_Group::check_enabled()) maillage.check_mesh();
        }

      int m = diviser_grandes_aretes(maillage);
      if (Comm_Group::check_enabled()) maillage.check_mesh();
      if (m > 0)
        {
          varVolume.resize_array(maillage.nb_sommets());
          varVolume = 0.;
          if(!surfactant.get_only_remaillage())
            barycentrer_lisser_apres_remaillage(maillage, varVolume);
          if (Comm_Group::check_enabled()) maillage.check_mesh();
        }
      if (Process::je_suis_maitre())
        Journal() << "remaillage_local_interface t= " << temps << " suppressions: " << n << " divisions: " << m << finl;
      if (!surfactant.get_disable_surfactant())
        nettoyer_maillage(maillage);
    }
  if (!surfactant.get_disable_surfactant())
    {
      ArrOfDouble surfactant_apres_remaillage = surfactant.check_conservation(maillage);
      surfactant.correction_conservation_globale(maillage, surfactant_avant_remaillage, surfactant_apres_remaillage);
    }
  nettoyer_maillage(maillage);
  statistiques().end_count(remaillage_loc_interf_counter_);
}

Vecteur3 Remaillage_FT_IJK::get_delta_euler(const Maillage_FT_IJK& maillage) const
{
  OBS_PTR(IJK_Splitting) s = maillage.ref_splitting();
  const IJK_Grid_Geometry& geom = s->get_grid_geometry();
  //  const IJK_Grid_Geometry & geom = s.get_grid_geometry();
  Vecteur3 delta(0., 0., 0.);
  const int dim = Objet_U::dimension;
  for (int k = 0; k < dim; k++)
    delta[k] = geom.get_constant_delta(k);
  // Le carre de la diagonale des elements :
  //  diago_elem2_ = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
  return delta;
}

double Remaillage_FT_IJK::calculer_longueurIdeale2_arete(const Maillage_FT_Disc& maillage,
                                                         int som0,
                                                         double x, double y, double z) const
{
  double lgrId2 = 0.;
  if (facteur_longueur_ideale_ > 0.)
    {
      const Maillage_FT_IJK& maillage_ijk = ref_cast(Maillage_FT_IJK, maillage);
      const Vecteur3 delta_xv = get_delta_euler(maillage_ijk);
      const int dim = Objet_U::dimension;
      int k;
      if (equilateral_)
        {
          // On calcul la diagonale de l'element eulerien. Puis longueur au carre.
          for (k = 0; k < dim; k++)
            {
              lgrId2 += delta_xv[k]*delta_xv[k];
            }
        }
      else
        {
          const DoubleTab& sommets = maillage.sommets();
          const Vecteur3 xyz(x, y, z);
          Vecteur3 v(0., 0., 0.);
          double norme2 = 0.;
          for (k = 0; k < dim; k++)
            {
              v[k] = xyz[k] - sommets(som0, k);
              norme2 += v[k] * v[k];
            }
          if (norme2 == 0)
            {
              v[0] = 1.;
              v[1] = 1.;
              v[2] = dim==3 ? 1. : 0.;
              norme2 = dim;
            }
          double f = 1. / sqrt(norme2);
          norme2 = 0.;
          for (k = 0; k < dim; k++)
            {
              v[k] *= f * delta_xv[k];
              norme2 += v[k] * v[k];
            }
          lgrId2 = norme2;

        }
      lgrId2 *= facteur_longueur_ideale_ * facteur_longueur_ideale_;
    }
  else
    {
      Cerr << "Erreur Remaillage_FT_IJK::calculer_longueurIdeale2_arete" << finl;
      exit();
    }
  return lgrId2;
}
