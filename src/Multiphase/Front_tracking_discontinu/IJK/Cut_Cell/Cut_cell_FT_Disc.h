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

#ifndef Cut_cell_FT_Disc_included
#define Cut_cell_FT_Disc_included

#include <Descripteur_FT.h>
#include <Maillage_FT_Disc.h>
#include <IJK_Field.h>
#include <Domaine_IJK.h>
#include <ConstIJK_ptr.h>
#include <TRUSTTabFT_cut_cell.h>
#include <IJK_Field_simd_tools.h>
#include <IJK_Navier_Stokes_tools.h>

class IJK_FT_Post;
class Probleme_FTD_IJK_cut_cell;
class IJK_Interfaces;

/*! @brief : class Cut_cell_FT_Disc
 *
 *  Cette classe decrit un maillage cut-cell, c'est-a-dire permettant
 *  des tableaux a valeurs dans les mailles diphasiques uniquement.
 *
 */

class IJK_FT_Post;
class IJK_FT_cut_cell;
class IJK_Interfaces;

class Cut_cell_FT_Disc
{
public:
  Cut_cell_FT_Disc();

  void add_to_persistent_data(DoubleTabFT_cut_cell& field, int dimension);
  void add_to_transient_data(DoubleTabFT_cut_cell& field, int dimension);
  void add_to_lazy_data(DoubleTabFT_cut_cell& field, int dimension);
  void add_to_persistent_data(IntTabFT_cut_cell& field, int dimension);
  void add_to_transient_data(IntTabFT_cut_cell& field, int dimension);
  void add_to_lazy_data(IntTabFT_cut_cell& field, int dimension);

  void initialise(IJK_Interfaces& interfaces, Domaine_IJK& splitting, Domaine_IJK::Localisation loc);
  void initialise(IJK_Interfaces& interfaces, Domaine_IJK& splitting, Domaine_IJK::Localisation loc, const IJK_Field_double& old_indicatrice, const IJK_Field_double& next_indicatrice);

  int initialise_independent_index(const IJK_Field_double& old_indicatrice, const IJK_Field_double& next_indicatrice);
  void initialise_permutation();
  void initialise_processed();

  void set_coord();

  void update(const IJK_Field_double& old_indicatrice, const IJK_Field_double& next_indicatrice);

  int add_and_remove_local_elements(const IJK_Field_double& old_indicatrice, const IJK_Field_double& next_indicatrice);

  void resize_data(int size);
  void compute_virtual_independent_index();

  void initialise_schema_comm();
  int initialise_communications();
  void echange_espace_virtuel();
  void echange_espace_virtuel(MD_Vector_tools::Operations_echange op);

  void update_index_sorted_by_k();
  void update_index_sorted_by_indicatrice(const IJK_Field_double& old_indicatrice, const IJK_Field_double& next_indicatrice);
  void update_index_sorted_by_statut_diphasique(const IJK_Field_double& old_indicatrice, const IJK_Field_double& next_indicatrice);

  template<typename T>
  void fill_buffer_with_variable(const TRUSTTabFT<T>& array, int component = 0) const; // Attention : fonction const, mais modifie write_buffer_

  template<typename T>
  void fill_variable_with_buffer(TRUSTTabFT<T>& array, int component = 0) const;

  void remplir_indice_diphasique();
  void remove_dead_and_virtual_cells(const IJK_Field_double& next_indicatrice);

  bool verifier_coherence_coord_independent_index();
  bool verifier_taille_tableaux();

  void imprime_elements_diphasiques();
  void imprime_elements_distants();

  enum STATUT_DIPHASIQUE : int
  {
    REGULIER = 0,             // Les deux phases sont regulieres (hors des cas diphasiques suivants)
    MOURRANT = 1,             // Une des deux phases est mourrante, l'autre occupera toute la cellule au temps suivant
    MONOPHASIQUE = 2,         // La cellule est purement monophasique. Une des deux phases n'existe pas sur la cellule.
    NAISSANT = 3,             // Une des deux phases est naissante, l'autre occupe initialement toute la cellule
    DESEQUILIBRE_FINAL = 4,   // Une des deux phases forme une cellule petite au temps suivant, sans indication sur l'etat initial
    DESEQUILIBRE_INITIAL = 5, // Une des deux phases forme une cellule petite au temps initial uniquement
    count = 6
  };

  const IntTabFT_cut_cell& get_independent_index() const { return independent_index_; }
  int get_independent_index(int n) const { return independent_index_(n); }
  int get_ghost_size() const { return ghost_size_; }
  int get_n_loc() const { return n_loc_; }
  int get_n_tot() const { return n_tot_; }
  const IJK_Interfaces& get_interfaces() const;
  const Domaine_IJK& get_domaine() const;
  const Desc_Structure_FT& get_desc_structure() const;

  IJK_Field_double& get_write_buffer() const { return write_buffer_; }

  inline Int3 ijk_per_of_index(int i, int j, int k, int index) const;
  inline int next_index_ijk_per(int i, int j, int k, int index, int negative_ghost_size, int positive_ghost_size) const;
  inline Int3 get_ijk(int n) const;
  inline int get_n(int i, int j, int k) const;
  inline int get_n_face(int num_face, int n, int i, int j, int k) const;

  double indic_pure(const int i, const int j, const int k) const;

  inline bool within_ghost(int n, int negative_ghost_size, int positive_ghost_size) const;

  template <int _DIR_>
  inline bool within_ghost_(int n, int negative_ghost_size, int positive_ghost_size) const;

  inline int get_k_value_index(int k) const;
  inline int get_n_from_k_index(int index) const;
  inline int get_n_from_indicatrice_index(int index) const;
  inline int get_statut_diphasique_value_index(int statut_diphasique) const;
  inline int get_n_from_statut_diphasique_index(int index) const;

protected:
  template<typename T>
  static int find_value_unsorted(T value, const TRUSTTabFT<T>& array, int imin, int imax);

  static bool verifier_toujours_meme_signe_et_non_nul(const IntTabFT_cut_cell& array);
  static bool verifier_tableau_jamais_nul(const DoubleTabFT_cut_cell& array);
  static bool verifier_valide_permutation(const IntTabFT_cut_cell& array);
  static bool verifier_pas_de_doublons(const IntTabFT_cut_cell& array);

  template<typename T>
  static void apply_permutation(TRUSTTabFT<T>& array, const IntTabFT_cut_cell& permutation_indices, IntTabFT_cut_cell& processed);

  static bool verifier_si_ordonne(const IntTabFT_cut_cell& array, int imin, int imax);

  friend IJK_FT_Post;
  friend Probleme_FTD_IJK_cut_cell;

  // Champ IJK_Field utilise pour le post-traitement des champs cut-cell
  mutable IJK_Field_double write_buffer_;

  // Champ IJK_Field de l'indice dans la structure diphasique
  // Ce tableau semble pertinent pour acceder aux cellules diphasiques voisines.
  IJK_Field_int indice_diphasique_;

  Schema_Comm schema_comm_;
  Desc_Structure_FT desc_;

  OBS_PTR(IJK_Interfaces) ref_interfaces_;
  OBS_PTR(Domaine_IJK) ref_domaine_;
  Domaine_IJK::Localisation localisation_;
  int ghost_size_;

  // Compteur du nombre de mise a jour des champs diphasiques
  int processed_count_;

  // Les champs diphasiques sont stockes dans l'ordre :
  //    0 < i < n_loc_        elements locaux
  //   n_loc_ < i < n_tot_    elements virtuels
  int n_loc_;
  int n_tot_;

  // Tableaux pour permettre le parcours des cellules diphasiques pour une
  // hauteur k donnee (au sein du calcul des flux)
  IntTabFT_cut_cell_vector2 index_sorted_by_k_;
  IntTabFT k_value_index_;

  // Tableaux pour permettre le parcours des cellules diphasiques dans le
  // sens d'une indicatrice croissante
  DoubleTabFT_cut_cell_vector3 index_sorted_by_indicatrice_;

  // Tableaux pour permettre le parcours des cellules diphasiques selon le
  // statut_diphasique (pour les cellules naissantes ou mourrantes)
  IntTabFT_cut_cell_vector3 index_sorted_by_statut_diphasique_;
  IntTabFT statut_diphasique_value_index_;

  // independent_index_ est l'indice lineaire IJK de la cellule diphasique.
  // Pour les elements locaus, les indices sont determines directement
  // et les cellules sont ordonnees dans le sens des indices croissants.
  // Pour les elements virtuels, les indices sont recalcules a partir
  // du vecteur des coordonnees sans ordre particulier.
  IntTabFT_cut_cell_scalar independent_index_;

  // permutation_ est un tableau d'indices utilise pour le reordonnement des
  // cellules diphasiques dans le sens des indices lineaires croissants.
  IntTabFT_cut_cell_scalar permutation_;

  // processed_ indique pour chaque cellule l'iteration de l'ajout.
  // Egal a processed_count_ pour les elements nouvellement diphasiques.
  // Le signe de processed est modifie lors de la permutation des champs.
  // Il peut etre positif ou negatif mais tous les elements de processed_
  // doivent etre de meme signe
  IntTabFT_cut_cell_scalar processed_;

  // coord_ est utilise pour donner l'emplacement absolue des cellules lors
  // des communications entre processeurs.
  DoubleTabFT_cut_cell_vector3 coord_;

  // Tableau des champs diphasiques (persistant, avec reordonnement a chaque pas de temps)
  LIST(OBS_PTR(DoubleTabFT_cut_cell)) persistent_double_data_;
  LIST(OBS_PTR(IntTabFT_cut_cell)) persistent_int_data_;

  // Tableau des champs diphasiques (ephemeres, avec reinitilisation a chaque pas de temps)
  LIST(OBS_PTR(DoubleTabFT_cut_cell)) transient_double_data_;
  LIST(OBS_PTR(IntTabFT_cut_cell)) transient_int_data_;

  // Tableau des champs diphasiques (paresseux, sans modifications entre chaque pas de temps)
  LIST(OBS_PTR(DoubleTabFT_cut_cell)) lazy_double_data_;
  LIST(OBS_PTR(IntTabFT_cut_cell)) lazy_int_data_;
};

inline void Cut_cell_FT_Disc::add_to_persistent_data(DoubleTabFT_cut_cell& field, int dimension)
{
  field.resize(n_tot_, dimension);
  persistent_double_data_.add(field);
}

inline void Cut_cell_FT_Disc::add_to_transient_data(DoubleTabFT_cut_cell& field, int dimension)
{
  field.resize(n_tot_, dimension);
  transient_double_data_.add(field);
}

inline void Cut_cell_FT_Disc::add_to_lazy_data(DoubleTabFT_cut_cell& field, int dimension)
{
  field.resize(n_tot_, dimension);
  lazy_double_data_.add(field);
}

inline void Cut_cell_FT_Disc::add_to_persistent_data(IntTabFT_cut_cell& field, int dimension)
{
  field.resize(n_tot_, dimension);
  persistent_int_data_.add(field);
}

inline void Cut_cell_FT_Disc::add_to_transient_data(IntTabFT_cut_cell& field, int dimension)
{
  field.resize(n_tot_, dimension);
  transient_int_data_.add(field);
}

inline void Cut_cell_FT_Disc::add_to_lazy_data(IntTabFT_cut_cell& field, int dimension)
{
  field.resize(n_tot_, dimension);
  lazy_int_data_.add(field);
}

inline Int3 Cut_cell_FT_Disc::ijk_per_of_index(int i, int j, int k, int index) const
{
  assert(index >= 0);
  assert(index <= 26);
  if (index == 0)
    {
      return {i, j, k};
    }
  else
    {
      int per_z = index/9;
      int per_y = (index%9)/3;
      int per_x = (index%9)%3;

      assert((per_x == 0) || ref_domaine_->get_periodic_flag(0));
      assert((per_y == 0) || ref_domaine_->get_periodic_flag(1));
      assert((per_z == 0) || ref_domaine_->get_periodic_flag(2));

      int n_dir_x = ref_domaine_->get_nb_elem_local(0);
      int n_dir_y = ref_domaine_->get_nb_elem_local(1);
      int n_dir_z = ref_domaine_->get_nb_elem_local(2);
      assert((per_x == 0) || (n_dir_x == ref_domaine_->get_nb_elem_tot(0)));
      assert((per_y == 0) || (n_dir_y == ref_domaine_->get_nb_elem_tot(1)));
      assert((per_z == 0) || (n_dir_z == ref_domaine_->get_nb_elem_tot(2)));

      assert((per_x == 0) || ((i < ghost_size_) || (i >= n_dir_x - ghost_size_)));
      assert((per_y == 0) || ((j < ghost_size_) || (j >= n_dir_y - ghost_size_)));
      assert((per_z == 0) || ((k < ghost_size_) || (k >= n_dir_z - ghost_size_)));
      int i_per = (per_x == 0)*i + (per_x == 1)*(n_dir_x + i) + (per_x == 2)*(i - n_dir_x);
      int j_per = (per_y == 0)*j + (per_y == 1)*(n_dir_y + j) + (per_y == 2)*(j - n_dir_y);
      int k_per = (per_z == 0)*k + (per_z == 1)*(n_dir_z + k) + (per_z == 2)*(k - n_dir_z);
      return {i_per, j_per, k_per};

    }
}

inline int Cut_cell_FT_Disc::next_index_ijk_per(int i, int j, int k, int index, int negative_ghost_size, int positive_ghost_size) const
{
  // per =
  //   0: no change in the given direction
  //   1: periodic towards the negative values in the given direction
  //   2: periodic towards the positive values in the given direction
  //
  // index = per_z * 9 + per_y * 3 + per_x
  //
  if (index == 26)
    {
      return -1;
    }
  else
    {
      int max_dir = 3;
      int min_dir = ((index >= 16) ? 2 : ((index >= 10) ? 1 : 0));
      assert(index >= 0);
      assert(min_dir >= 0 && min_dir <=2);
      for (int dir_1 = min_dir; dir_1 < max_dir ; dir_1++)
        {
          if (ref_domaine_->get_periodic_flag(dir_1))
            {
              int n_dir_1 = ref_domaine_->get_nb_elem_local(dir_1);
              int n_dir_tot_1 = ref_domaine_->get_nb_elem_tot(dir_1);

              if (n_dir_1 == n_dir_tot_1)
                {
                  int i_dir_1 = select_dir(dir_1, i, j, k);

                  if (i_dir_1 < positive_ghost_size)
                    {
                      int per_x = (dir_1 == 0)*1;
                      int per_y = (dir_1 == 1)*1;
                      int per_z = (dir_1 == 2)*1;

                      int next_index = per_x + per_y*3 + per_z*9;
                      if (next_index > index)
                        {
                          return next_index;
                        }

                      int per_x_debut_boucle_2 = per_x;
                      int per_y_debut_boucle_2 = per_y;
                      int per_z_debut_boucle_2 = per_z;
                      for (int dir_2 = 0; dir_2 < dir_1 ; dir_2++)
                        {
                          if (ref_domaine_->get_periodic_flag(dir_2))
                            {
                              int n_dir_2 = ref_domaine_->get_nb_elem_local(dir_2);
                              int n_dir_tot_2 = ref_domaine_->get_nb_elem_tot(dir_2);

                              if (n_dir_2 == n_dir_tot_2)
                                {
                                  int i_dir_2 = select_dir(dir_2, i, j, k);

                                  if (i_dir_2 < positive_ghost_size)
                                    {
                                      per_x = per_x_debut_boucle_2 + (dir_2 == 0)*1;
                                      per_y = per_y_debut_boucle_2 + (dir_2 == 1)*1;
                                      per_z = per_z_debut_boucle_2 + (dir_2 == 2)*1;

                                      next_index = per_x + per_y*3 + per_z*9;
                                      if (next_index > index)
                                        {
                                          return next_index;
                                        }

                                      int per_x_debut_boucle_3 = per_x;
                                      int per_y_debut_boucle_3 = per_y;
                                      int per_z_debut_boucle_3 = per_z;
                                      for (int dir_3 = 0; dir_3 < dir_2 ; dir_3++)
                                        {
                                          if (ref_domaine_->get_periodic_flag(dir_3))
                                            {
                                              int n_dir_3 = ref_domaine_->get_nb_elem_local(dir_3);
                                              int n_dir_tot_3 = ref_domaine_->get_nb_elem_tot(dir_3);

                                              if (n_dir_3 == n_dir_tot_3)
                                                {
                                                  int i_dir_3 = select_dir(dir_3, i, j, k);

                                                  if (i_dir_3 < positive_ghost_size)
                                                    {
                                                      per_x = per_x_debut_boucle_3 + (dir_3 == 0)*1;
                                                      per_y = per_y_debut_boucle_3 + (dir_3 == 1)*1;
                                                      per_z = per_z_debut_boucle_3 + (dir_3 == 2)*1;

                                                      next_index = per_x + per_y*3 + per_z*9;
                                                      if (next_index > index)
                                                        {
                                                          return next_index;
                                                        }
                                                    }
                                                  if (i_dir_3 >= n_dir_3 - negative_ghost_size)
                                                    {
                                                      per_x = per_x_debut_boucle_3 + (dir_3 == 0)*2;
                                                      per_y = per_y_debut_boucle_3 + (dir_3 == 1)*2;
                                                      per_z = per_z_debut_boucle_3 + (dir_3 == 2)*2;

                                                      next_index = per_x + per_y*3 + per_z*9;
                                                      if (next_index > index)
                                                        {
                                                          return next_index;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                  if (i_dir_2 >= n_dir_2 - negative_ghost_size)
                                    {
                                      per_x = per_x_debut_boucle_2 + (dir_2 == 0)*2;
                                      per_y = per_y_debut_boucle_2 + (dir_2 == 1)*2;
                                      per_z = per_z_debut_boucle_2 + (dir_2 == 2)*2;

                                      next_index = per_x + per_y*3 + per_z*9;
                                      if (next_index > index)
                                        {
                                          return next_index;
                                        }

                                      int per_x_debut_boucle_3 = per_x;
                                      int per_y_debut_boucle_3 = per_y;
                                      int per_z_debut_boucle_3 = per_z;
                                      for (int dir_3 = 0; dir_3 < dir_2 ; dir_3++)
                                        {
                                          if (ref_domaine_->get_periodic_flag(dir_3))
                                            {
                                              int n_dir_3 = ref_domaine_->get_nb_elem_local(dir_3);
                                              int n_dir_tot_3 = ref_domaine_->get_nb_elem_tot(dir_3);

                                              if (n_dir_3 == n_dir_tot_3)
                                                {
                                                  int i_dir_3 = select_dir(dir_3, i, j, k);

                                                  if (i_dir_3 < positive_ghost_size)
                                                    {
                                                      per_x = per_x_debut_boucle_3 + (dir_3 == 0)*1;
                                                      per_y = per_y_debut_boucle_3 + (dir_3 == 1)*1;
                                                      per_z = per_z_debut_boucle_3 + (dir_3 == 2)*1;

                                                      next_index = per_x + per_y*3 + per_z*9;
                                                      if (next_index > index)
                                                        {
                                                          return next_index;
                                                        }
                                                    }
                                                  if (i_dir_3 >= n_dir_3 - negative_ghost_size)
                                                    {
                                                      per_x = per_x_debut_boucle_3 + (dir_3 == 0)*2;
                                                      per_y = per_y_debut_boucle_3 + (dir_3 == 1)*2;
                                                      per_z = per_z_debut_boucle_3 + (dir_3 == 2)*2;

                                                      next_index = per_x + per_y*3 + per_z*9;
                                                      if (next_index > index)
                                                        {
                                                          return next_index;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                  if (i_dir_1 >= n_dir_1 - negative_ghost_size)
                    {
                      int per_x = (dir_1 == 0)*2;
                      int per_y = (dir_1 == 1)*2;
                      int per_z = (dir_1 == 2)*2;

                      int next_index = per_x + per_y*3 + per_z*9;
                      if (next_index > index)
                        {
                          return next_index;
                        }

                      int per_x_debut_boucle_2 = per_x;
                      int per_y_debut_boucle_2 = per_y;
                      int per_z_debut_boucle_2 = per_z;
                      for (int dir_2 = 0; dir_2 < dir_1 ; dir_2++)
                        {
                          if (ref_domaine_->get_periodic_flag(dir_2))
                            {
                              int n_dir_2 = ref_domaine_->get_nb_elem_local(dir_2);
                              int n_dir_tot_2 = ref_domaine_->get_nb_elem_tot(dir_2);

                              if (n_dir_2 == n_dir_tot_2)
                                {
                                  int i_dir_2 = select_dir(dir_2, i, j, k);

                                  if (i_dir_2 < positive_ghost_size)
                                    {
                                      per_x = per_x_debut_boucle_2 + (dir_2 == 0)*1;
                                      per_y = per_y_debut_boucle_2 + (dir_2 == 1)*1;
                                      per_z = per_z_debut_boucle_2 + (dir_2 == 2)*1;

                                      next_index = per_x + per_y*3 + per_z*9;
                                      if (next_index > index)
                                        {
                                          return next_index;
                                        }

                                      int per_x_debut_boucle_3 = per_x;
                                      int per_y_debut_boucle_3 = per_y;
                                      int per_z_debut_boucle_3 = per_z;
                                      for (int dir_3 = 0; dir_3 < dir_2 ; dir_3++)
                                        {
                                          if (ref_domaine_->get_periodic_flag(dir_3))
                                            {
                                              int n_dir_3 = ref_domaine_->get_nb_elem_local(dir_3);
                                              int n_dir_tot_3 = ref_domaine_->get_nb_elem_tot(dir_3);

                                              if (n_dir_3 == n_dir_tot_3)
                                                {
                                                  int i_dir_3 = select_dir(dir_3, i, j, k);

                                                  if (i_dir_3 < positive_ghost_size)
                                                    {
                                                      per_x = per_x_debut_boucle_3 + (dir_3 == 0)*1;
                                                      per_y = per_y_debut_boucle_3 + (dir_3 == 1)*1;
                                                      per_z = per_z_debut_boucle_3 + (dir_3 == 2)*1;

                                                      next_index = per_x + per_y*3 + per_z*9;
                                                      if (next_index > index)
                                                        {
                                                          return next_index;
                                                        }
                                                    }
                                                  if (i_dir_3 >= n_dir_3 - negative_ghost_size)
                                                    {
                                                      per_x = per_x_debut_boucle_3 + (dir_3 == 0)*2;
                                                      per_y = per_y_debut_boucle_3 + (dir_3 == 1)*2;
                                                      per_z = per_z_debut_boucle_3 + (dir_3 == 2)*2;

                                                      next_index = per_x + per_y*3 + per_z*9;
                                                      if (next_index > index)
                                                        {
                                                          return next_index;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                  if (i_dir_2 >= n_dir_2 - negative_ghost_size)
                                    {
                                      per_x = per_x_debut_boucle_2 + (dir_2 == 0)*2;
                                      per_y = per_y_debut_boucle_2 + (dir_2 == 1)*2;
                                      per_z = per_z_debut_boucle_2 + (dir_2 == 2)*2;

                                      next_index = per_x + per_y*3 + per_z*9;
                                      if (next_index > index)
                                        {
                                          return next_index;
                                        }

                                      int per_x_debut_boucle_3 = per_x;
                                      int per_y_debut_boucle_3 = per_y;
                                      int per_z_debut_boucle_3 = per_z;
                                      for (int dir_3 = 0; dir_3 < dir_2 ; dir_3++)
                                        {
                                          if (ref_domaine_->get_periodic_flag(dir_3))
                                            {
                                              int n_dir_3 = ref_domaine_->get_nb_elem_local(dir_3);
                                              int n_dir_tot_3 = ref_domaine_->get_nb_elem_tot(dir_3);

                                              if (n_dir_3 == n_dir_tot_3)
                                                {
                                                  int i_dir_3 = select_dir(dir_3, i, j, k);

                                                  if (i_dir_3 < positive_ghost_size)
                                                    {
                                                      per_x = per_x_debut_boucle_3 + (dir_3 == 0)*1;
                                                      per_y = per_y_debut_boucle_3 + (dir_3 == 1)*1;
                                                      per_z = per_z_debut_boucle_3 + (dir_3 == 2)*1;

                                                      next_index = per_x + per_y*3 + per_z*9;
                                                      if (next_index > index)
                                                        {
                                                          return next_index;
                                                        }
                                                    }
                                                  if (i_dir_3 >= n_dir_3 - negative_ghost_size)
                                                    {
                                                      per_x = per_x_debut_boucle_3 + (dir_3 == 0)*2;
                                                      per_y = per_y_debut_boucle_3 + (dir_3 == 1)*2;
                                                      per_z = per_z_debut_boucle_3 + (dir_3 == 2)*2;

                                                      next_index = per_x + per_y*3 + per_z*9;
                                                      if (next_index > index)
                                                        {
                                                          return next_index;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
      return -1;
    }
}

inline Int3 Cut_cell_FT_Disc::get_ijk(int n) const
{
  int independent_index = independent_index_(n);
  Int3 ijk = ref_domaine_->get_ijk_from_independent_index(independent_index);

  return ijk;
}

inline int Cut_cell_FT_Disc::get_n(int i, int j, int k) const
{
  // Note: La condition sur n_tot_ est ajoute pour le cas ou indice_diphasique_ n'est pas encore alloue.
  int indice_diphasique = n_tot_ == 0 ? -1 : indice_diphasique_(i,j,k);
  return indice_diphasique;
}

inline int Cut_cell_FT_Disc::get_n_face(int num_face, int n, int i, int j, int k) const
{
  int decalage = num_face/3;    //  0 pour 0,1,2 et  1 pour 3,4,5
  if (!decalage)
    {
      return n;
    }
  else
    {
      int dir = num_face%3;
      int di = decalage*(dir == 0);
      int dj = decalage*(dir == 1);
      int dk = decalage*(dir == 2);

      int n_face = get_n(i+di,j+dj,k+dk);
      return n_face;
    }
}

inline bool Cut_cell_FT_Disc::within_ghost(int n, int negative_ghost_size, int positive_ghost_size) const
{
  Int3 ijk = get_ijk(n);

  int i = ijk[0];
  int j = ijk[1];
  int k = ijk[2];

  return ref_domaine_->within_ghost(i, j, k, negative_ghost_size, positive_ghost_size);
}

template <int _DIR_>
inline bool Cut_cell_FT_Disc::within_ghost_(int n, int negative_ghost_size, int positive_ghost_size) const
{
  Int3 ijk = get_ijk(n);

  int i = ijk[0];
  int j = ijk[1];
  int k = ijk[2];

  return ref_domaine_->within_ghost_<_DIR_>(i, j, k, negative_ghost_size, positive_ghost_size);
}

inline int Cut_cell_FT_Disc::get_k_value_index(int k) const
{
  assert(k <= ref_domaine_->get_nb_elem_local(2) + 2*ghost_size_);
  return k_value_index_(k + ghost_size_);
}

inline int Cut_cell_FT_Disc::get_n_from_k_index(int index) const
{
  assert(index >= 0);
  return index_sorted_by_k_(index,0);
}

inline int Cut_cell_FT_Disc::get_n_from_indicatrice_index(int index) const
{
  assert(index >= 0);

  union
  {
    double d;
    long long int i;
  } u;
  u.d = index_sorted_by_indicatrice_(index,0);

  // Retrieving the bits of the index, which were stored into a double variable.
  return (int)(u.i);
}

inline int Cut_cell_FT_Disc::get_statut_diphasique_value_index(int statut_diphasique) const
{
  assert(statut_diphasique <= static_cast<int>(STATUT_DIPHASIQUE::count));
  return statut_diphasique_value_index_(statut_diphasique);
}

inline int Cut_cell_FT_Disc::get_n_from_statut_diphasique_index(int index) const
{
  assert(index >= 0);
  return index_sorted_by_statut_diphasique_(index,0);
}

#endif /* Cut_cell_FT_Disc_included */
