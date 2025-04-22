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

#ifndef Maillage_FT_IJK_included
#define Maillage_FT_IJK_included

#include <Objet_U.h>
#include <Maillage_FT_Disc.h>
#include <Linear_algebra_tools.h>
#include <Domaine_IJK.h>
#include <TRUSTTab.h>
#include <FT_Field.h>
#include <Operator_FT_Disc.h>
class Domaine_dis_base;

class Parcours_interface;

/*! @brief : class Maillage_FT_IJK
 *
 */
class Maillage_FT_IJK : public Maillage_FT_Disc
{

  Declare_instanciable(Maillage_FT_IJK) ;

public:
  FT_Field Surfactant_facettes_;
  Maillage_FT_IJK(const Maillage_FT_IJK&) = default;
  void initialize(const Domaine_IJK&, const Domaine_dis_base&, const Parcours_interface&, const bool use_tryggvason_interfacial_source=false);
  const ArrOfInt& compo_connexe_facettes() const
  {
    return compo_connexe_facettes_;
  };
  ArrOfInt& compo_connexe_facettes_non_const()
  {
    return compo_connexe_facettes_;
  };
  const FT_Field& Surfactant_facettes() const
  {
    return Surfactant_facettes_;
  };
  FT_Field& Surfactant_facettes_non_const()
  {
    return Surfactant_facettes_;
  };

  void update_gradient_laplacien_Surfactant()
  {
    Surfactant_facettes_.update_gradient_laplacien_FT(*this);
  };
  DoubleTab update_sigma_and_interfacial_source_term_sommet(const Domaine_IJK& splitting, bool compute_interfacial_source, bool use_tryggvason_formulation, const double sigma_const = -1.);
  void set_Surfactant_facettes(ArrOfDouble Surfactant_field);
  void set_Surfactant_facettes_sommets(ArrOfDouble Surfactant_field);
// Surcharge de Maillage_FT_Disc:
  void supprimer_facettes(const ArrOfInt& liste_facettes);


  void nettoyer_maillage() override;
  void sauv_facette_indexation_avant_transport();
  void corriger_proprietaires_facettes();
  void update_surfactant_apres_transport();
  void parcourir_maillage();
  void transporter(const DoubleTab& deplacement) override;
  void deplacer_sommets(const ArrOfInt& liste_sommets_initiale,
                        const DoubleTab& deplacement_initial,
                        ArrOfInt& liste_sommets_sortis,
                        ArrOfInt& numero_face_sortie, int skip_facettes=0) override;
// Surcharge de maillage_ft_disc pour conserver les composantes connexes.
  void recopie(const Maillage_FT_Disc& source_mesh, Statut_Maillage niveau_copie) override;

// Surcharge de Maillage_FT_Disc :
// surcharge: initialise les composantes connexes et appelle la methode ajouter_maillage_IJK.
  void ajouter_maillage(const Maillage_FT_Disc& maillage_tmp,int skip_facettes=0) override;
  void echanger_facettes(const ArrOfInt& liste_facettes,
                         const ArrOfInt& liste_elem_arrivee,
                         ArrOfInt& facettes_recues_numfacettes,
                         ArrOfInt& facettes_recues_numelement);
// ajout et gestion du tableau des compo_connexes a partir du tableau compo_connex du maillage source
  void ajouter_maillage_IJK(const Maillage_FT_IJK& added_mesh);

  void lire_maillage_ft_dans_lata(const char *filename_with_path, int tstep,
                                  const char *geometryname);

  void set_ijk_cell_index(int num_sommet, Int3 ijk)
  {
    const Domaine_IJK& s = ref_domaine_.valeur();
    sommet_elem_[num_sommet] = s.convert_ijk_cell_to_packed(ijk[0],ijk[1],ijk[2]);
  }
  void set_barycentrage(bool bary)
  {
    during_barycentrage_=bary;
  }
  Int3 get_ijk_cell_index(int num_sommet) const
  {
    const Domaine_IJK& s = ref_domaine_.valeur();
    int index = sommet_elem_[num_sommet];
    return s.convert_packed_to_ijk_cell(index);
  }

  void initialize_processor_neighbourhood();
// Surcharges de Maillage_FT_Disc :
  int check_sommets(int error_is_fatal = 1) const override;
  int check_mesh(int error_is_fatal = 1, int skip_facette_owner = 0, int skip_facettes = 0) const override;

  void calculer_compo_connexe_sommets(ArrOfIntFT& compo_connexe_sommets) const;

  void recopie_force_compo(const Maillage_FT_IJK& source_mesh, const int icompo);
// Methode qui modifie l'attribut compo_connexe_facettes_
  void set_composante_connexe(const int i_facette, const int icompo)
  {
    compo_connexe_facettes_[i_facette] = icompo;
  };

  double minimum_longueur_arrete() const;
  int nb_facettes_sans_duplicata() const;

  const Domaine_IJK& get_domaine() const { return ref_domaine_.valeur(); }

protected:
  // Surcharge de Maillage_FT_Disc :
  bool during_barycentrage_ = false;
  bool use_tryggvason_interfacial_source_=false;
  void   calculer_costheta_minmax(DoubleTab& costheta) const override;

  const Maillage_FT_IJK& operator=(const Maillage_FT_IJK&)
  {
    Cerr << "Maillage_FT_IJK& operator=" << finl;
    Process::exit();
    return *this;
  }    // Interdit !
  void creer_facettes_virtuelles(const ArrOfInt& liste_facettes,
                                 const ArrOfInt& liste_pe,
                                 const ArrOfInt& facettes_send_pe_list,
                                 const ArrOfInt& facettes_recv_pe_list) override;

  OBS_PTR(Domaine_IJK) ref_domaine_;

// Taille du domaine IJK sur chaque proc:
  int nbmailles_euler_i_ = 0;
  int nbmailles_euler_j_ = 0;
  int nbmailles_euler_k_ = 0;
// Tableaux des processeurs voisins :
  ArrOfInt voisinage_processeur_;
  ArrOfInt liste_processeurs_voisins_faces_;
  ArrOfInt liste_processeurs_voisins_aretes_;
  ArrOfInt liste_processeurs_voisins_coins_;

  // Pour chaque facette d'interface, numero de la composante connexe (numero de la bulle)
  ArrOfIntFT compo_connexe_facettes_;

  DoubleTab indexation_facettes_avant_transport_;

};
#endif /* Maillage_FT_IJK_included */
