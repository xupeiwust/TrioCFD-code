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
#include <IJK_FT_Post.h>

/*! @brief : class Cut_cell_FT_Disc
 *
 *  <Description of class Cut_cell_FT_Disc>
 *
 *
 *  Cette classe decrit un maillage cut-cell, c'est-a-dire permettant
 *  des tableaux a valeurs dans les mailles diphasiques uniquement.
 *
 */

struct IJK;

class Cut_cell_FT_Disc
{
public:
  Cut_cell_FT_Disc(IJK_Splitting& splitting);

  void initialise(const IJK_Field_double& indicatrice, const FixedVector<IJK_Field_double, 3>& coord, const IJK_Field_double& rho);

  int initialise_linear_index(const IJK_Field_double& indicatrice);
  void initialise_permutation();
  void initialise_processed();

  void set_coord(const FixedVector<IJK_Field_double, 3>& coord);
  void set_rho(const IJK_Field_double& rho);

  void update(const IJK_Field_double& indicatrice, const FixedVector<IJK_Field_double, 3>& coord, const IJK_Field_double& rho);

  int add_and_remove_local_elements(const IJK_Field_double& indicatrice);

  void resize_data(int size);
  void compute_virtual_linear_index();

  void initialise_schema_comm();
  int initialise_communications();
  void echange_espace_virtuel();

  template<typename T>
  void fill_buffer_with_variable(const TRUSTTabFT<T>& array, int component = 0);

  bool verifier_coherence_coord_linear_index();
  bool verifier_taille_tableaux();

  void imprime_elements_diphasiques();
  void imprime_elements_distants();

protected:
  static inline int get_i_selon_dir(const int direction, const double coord_dir, const int ghost_size, const IJK_Splitting& splitting);
  static inline IJK get_ijk_from_coord(const double coord_x, const double coord_y, const double coord_z, const int ghost_size, const IJK_Splitting& splitting, bool expect_local_value);
  static inline IJK get_ijk_from_linear_index(const int linear_index, const int ghost_size, const IJK_Splitting& splitting, bool expect_local_value);
  static inline int get_linear_index(const int i, const int j, const int k, const int ghost_size, const IJK_Splitting& splitting, bool expect_local_value);

  template<typename T>
  static int find_value(T value, const TRUSTTabFT<T>& array, int imin, int imax);

  template<typename T>
  static int find_value_unsorted(T value, const TRUSTTabFT<T>& array, int imin, int imax);

  static bool verifier_toujours_meme_signe_et_non_nul(const IntTabFT& array);
  static bool verifier_tableau_jamais_nul(const DoubleTabFT& array);
  static bool verifier_valide_permutation(const IntTabFT& array);

  template<typename T>
  static void apply_permutation(TRUSTTabFT<T>& array, const IntTabFT& permutation_indices, IntTabFT& processed);

  static bool verifier_si_ordonne(const IntTabFT& array, int imin, int imax);

  friend IJK_FT_Post;

  // Champ IJK_Field utilise pour le post-traitement des champs cut-cell
  IJK_Field_double write_buffer_;

  Schema_Comm_FT schema_comm_;
  Desc_Structure_FT desc_;

  IJK_Splitting& splitting_;
  int ghost_size_;

  // Delai pour la suppression des elements anciennement diphasiques
  int cell_removal_time_offset_;

  // Compteur du nombre de mise a jour des champs diphasiques
  int processed_count_;

  // Les champs diphasiques sont stockes dans l'ordre :
  //    0 < i < n_loc_        elements locaux
  //   n_loc_ < i < n_tot_    elements virtuels
  int n_loc_;
  int n_tot_;

  // linear_index_ est l'indice lineaire IJK de la cellule diphasique.
  // Pour les elements locaus, les indices sont determines directement
  // et les cellules sont ordonnees dans le sens des indices croissants.
  // Pour les elements virtuels, les indices sont recalcules a partir
  // du vecteur des coordonnees sans ordre particulier.
  IntTabFT linear_index_;

  // permutation_ est un tableau d'indices utilise pour le reordonnement des
  // cellules diphasiques dans le sens des indices lineaires croissants.
  IntTabFT permutation_;

  // processed_ indique pour chaque cellule l'iteration de la derniere
  // mise a jour (pour la suppression des anciens elements).
  // Egal a processed_count_ pour les elements actuellement diphasiques.
  // Le signe de processed est modifie lors de la permutation des champs.
  // Il peut etre positif ou negatif mais tous les elements de processed_
  // doivent etre de meme signe
  IntTabFT processed_;

  // Tableau des champs diphasiques
  DoubleTabFT *data_[2];
  const int number_of_fields_ = 2;

  // Champs diphasiques inclus dans data_
  DoubleTabFT coord_;
  DoubleTabFT rho_;
};

#endif /* Cut_cell_FT_Disc_included */
