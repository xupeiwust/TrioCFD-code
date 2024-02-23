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
#include <IJK_Interfaces.h>
#include <IJK_Splitting.h>
#include <Champ_diphasique.h>

/*! @brief : class Cut_cell_FT_Disc
 *
 *  <Description of class Cut_cell_FT_Disc>
 *
 *
 *  Cette classe decrit un maillage cut-cell, c'est-a-dire permettant
 *  des tableaux a valeurs dans les mailles diphasiques uniquement.
 *
 */

class IJK_FT_Post;
class IJK_FT_cut_cell;

class Cut_cell_FT_Disc
{
public:
  Cut_cell_FT_Disc(IJK_Interfaces& interfaces, IJK_Splitting& splitting);

  void add_to_data(DoubleTabFT_cut_cell& field, int dimension);
  void add_to_lazy_data(DoubleTabFT_cut_cell& field, int dimension);
  void add_to_int_data(IntTabFT_cut_cell& field, int dimension);
  void add_to_lazy_int_data(IntTabFT_cut_cell& field, int dimension);

  void initialise(const IJK_Field_double& indicatrice_old, const IJK_Field_double& indicatrice_next, const FixedVector<IJK_Field_double, 3>& coord);

  int initialise_linear_index(const IJK_Field_double& indicatrice_old, const IJK_Field_double& indicatrice_next);
  void initialise_permutation();
  void initialise_processed();

  void set_coord(const FixedVector<IJK_Field_double, 3>& coord);

  void update(const IJK_Field_double& indicatrice_old, const IJK_Field_double& indicatrice_next, const FixedVector<IJK_Field_double, 3>& coord);

  int add_and_remove_local_elements(const IJK_Field_double& indicatrice_old, const IJK_Field_double& indicatrice_next);

  void resize_data(int size);
  void compute_virtual_linear_index();

  void initialise_schema_comm();
  int initialise_communications();
  void echange_espace_virtuel();

  template<typename T>
  void fill_buffer_with_variable(const TRUSTTabFT<T>& array, int component = 0);
  void remplir_indice_diphasique();

  bool verifier_coherence_coord_linear_index();
  bool verifier_taille_tableaux();

  void imprime_elements_diphasiques();
  void imprime_elements_distants();

  const IntTabFT_cut_cell& get_linear_index() const { return linear_index_; }
  int get_linear_index(int n) const { return linear_index_(n); }
  int get_ghost_size() const { return ghost_size_; }
  int get_n_loc() const { return n_loc_; }
  int get_n_tot() const { return n_tot_; }
  const IJK_Interfaces& get_interfaces() const { return interfaces_; }
  const IJK_Splitting& get_splitting() const { return splitting_; }
  const Desc_Structure_FT& get_desc_structure() const { return desc_; }

  inline Int3 get_ijk(int n) const;
  inline int get_n(int i, int j, int k) const;
  inline int get_linear_index(int i, int j, int k) const;
  inline int get_n_face(int num_face, int n, int i, int j, int k) const;

  static inline int get_i_selon_dir(int direction, double coord_dir, int ghost_size, const IJK_Splitting& splitting);
  static inline Int3 get_ijk_from_coord(double coord_x, double coord_y, double coord_z, int ghost_size, const IJK_Splitting& splitting, bool expect_local_value);
  static inline Int3 get_ijk_from_linear_index(int linear_index, int ghost_size, const IJK_Splitting& splitting, bool expect_local_value);
  static inline int get_linear_index(int i, int j, int k, int ghost_size, const IJK_Splitting& splitting, bool expect_local_value);


protected:
  template<typename T>
  static int find_value_unsorted(T value, const TRUSTTabFT<T>& array, int imin, int imax);

  static bool verifier_toujours_meme_signe_et_non_nul(const IntTabFT_cut_cell& array);
  static bool verifier_tableau_jamais_nul(const DoubleTabFT_cut_cell& array);
  static bool verifier_valide_permutation(const IntTabFT_cut_cell& array);

  template<typename T>
  static void apply_permutation(TRUSTTabFT<T>& array, const IntTabFT_cut_cell& permutation_indices, IntTabFT_cut_cell& processed);

  static bool verifier_si_ordonne(const IntTabFT_cut_cell& array, int imin, int imax);

  friend IJK_FT_Post;
  friend IJK_FT_cut_cell;

  // Champ IJK_Field utilise pour le post-traitement des champs cut-cell
  IJK_Field_double write_buffer_;

  // Champ IJK_Field de l'indice dans la structure diphasique
  // Ce tableau semble pertinent pour acceder aux cellules diphasiques voisines.
  IJK_Field_int indice_diphasique_;

  Schema_Comm_FT schema_comm_;
  Desc_Structure_FT desc_;

  IJK_Interfaces& interfaces_;
  IJK_Splitting& splitting_;
  int ghost_size_;

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
  IntTabFT_cut_cell_scalar linear_index_;

  // permutation_ est un tableau d'indices utilise pour le reordonnement des
  // cellules diphasiques dans le sens des indices lineaires croissants.
  IntTabFT_cut_cell_scalar permutation_;

  // processed_ indique pour chaque cellule l'iteration de l'ajout.
  // Egal a processed_count_ pour les elements nouvellement diphasiques.
  // Le signe de processed est modifie lors de la permutation des champs.
  // Il peut etre positif ou negatif mais tous les elements de processed_
  // doivent etre de meme signe
  IntTabFT_cut_cell_scalar processed_;

  // Tableau des champs diphasiques (persistant, avec reordonnement a chaque pas de temps)
  LIST(REF(DoubleTabFT_cut_cell)) data_;
  LIST(REF(IntTabFT_cut_cell)) int_data_;

  // Champs diphasiques inclus dans data_
  DoubleTabFT_cut_cell_vector3 coord_;

  // Tableau des champs diphasiques (paresseux, sans modifications entre chaque pas de temps)
  LIST(REF(DoubleTabFT_cut_cell)) lazy_data_;
  LIST(REF(IntTabFT_cut_cell)) lazy_int_data_;
};

inline void Cut_cell_FT_Disc::add_to_data(DoubleTabFT_cut_cell& field, int dimension)
{
  field.set_smart_resize(1);
  field.resize(n_tot_, dimension);
  data_.add(field);
}

inline void Cut_cell_FT_Disc::add_to_lazy_data(DoubleTabFT_cut_cell& field, int dimension)
{
  field.set_smart_resize(1);
  field.resize(n_tot_, dimension);
  lazy_data_.add(field);
}

inline void Cut_cell_FT_Disc::add_to_int_data(IntTabFT_cut_cell& field, int dimension)
{
  field.set_smart_resize(1);
  field.resize(n_tot_, dimension);
  int_data_.add(field);
}

inline void Cut_cell_FT_Disc::add_to_lazy_int_data(IntTabFT_cut_cell& field, int dimension)
{
  field.set_smart_resize(1);
  field.resize(n_tot_, dimension);
  lazy_int_data_.add(field);
}

inline int Cut_cell_FT_Disc::get_i_selon_dir(int direction, double coord_dir, int ghost_size, const IJK_Splitting& splitting)
{
  const int offset_dir = splitting.get_offset_local(direction);
  double origin_dir = splitting.get_grid_geometry().get_origin(direction);

  const int n = splitting.get_nb_elem_local(direction);

  int index = -1;
  if (splitting.get_grid_geometry().is_uniform(direction))
    {
      const double d = splitting.get_grid_geometry().get_constant_delta(direction);
      index = (int)((coord_dir - origin_dir)/d) - offset_dir;

      // Si invalide, essaye le periodique
      if ((index < -ghost_size) || (index > n + ghost_size))
        {
          if (splitting.get_grid_geometry().get_periodic_flag(direction))
            {
              const double domain_length_dir = splitting.get_grid_geometry().get_domain_length(direction);
              index = (int)(((domain_length_dir - coord_dir) - origin_dir)/d) - offset_dir;
              if ((index < -ghost_size) || (index > n + ghost_size))
                {
                  Cerr << "Error: In get_ijk_from_coord(), invalid index along direction " << direction << " (which is periodic)." << finl;
                  assert(false);
                }
            }
          else
            {
              Cerr << "Error: In get_ijk_from_coord(), invalid index along direction " << direction << " (which is not periodic)." << finl;
              assert(false);
            }
        }
    }
  else
    {
      Cerr << "Error: In get_ijk_from_coord(), the case of a non-uniform mesh along direction " << direction << " is not implemented." << finl;
      assert(false);
    }

  return index;
}

inline Int3 Cut_cell_FT_Disc::get_ijk_from_coord(double coord_x, double coord_y, double coord_z, int ghost_size, const IJK_Splitting& splitting, bool expect_local_value)
{
  int i = get_i_selon_dir(0, coord_x, ghost_size, splitting);
  int j = get_i_selon_dir(1, coord_y, ghost_size, splitting);
  int k = get_i_selon_dir(2, coord_z, ghost_size, splitting);

  if (expect_local_value)
    {
      assert((i >= 0) && (i < splitting.get_nb_elem_local(0)));
      assert((j >= 0) && (j < splitting.get_nb_elem_local(1)));
      assert((k >= 0) && (k < splitting.get_nb_elem_local(2)));
    }
  else
    {
      assert((i >= -ghost_size) && (i < splitting.get_nb_elem_local(0) + ghost_size));
      assert((j >= -ghost_size) && (j < splitting.get_nb_elem_local(1) + ghost_size));
      assert((k >= -ghost_size) && (k < splitting.get_nb_elem_local(2) + ghost_size));
    }

  Int3 ijk = {i, j, k};
  return ijk;
}

inline Int3 Cut_cell_FT_Disc::get_ijk_from_linear_index(int linear_index, int ghost_size, const IJK_Splitting& splitting, bool expect_local_value)
{
  const int ni = splitting.get_nb_elem_local(0);
  const int nj = splitting.get_nb_elem_local(1);

  // Computation of the indices disregarding the offset (i goes from 0 to ni + 2*ghost_size)
  int k = (linear_index)/((ni + 2*ghost_size)*(nj + 2*ghost_size));
  int j = ((linear_index) - (ni + 2*ghost_size)*(nj + 2*ghost_size)*k)/(ni + 2*ghost_size);
  int i = (linear_index) - (ni + 2*ghost_size)*(nj + 2*ghost_size)*k - (ni + 2*ghost_size)*j;

  // Computation of the indices with the offset (i goes from -ghost_size_ to ni + ghost_size)
  k -= ghost_size;
  j -= ghost_size;
  i -= ghost_size;

  if (expect_local_value)
    {
      assert((i >= 0) && (i < splitting.get_nb_elem_local(0)));
      assert((j >= 0) && (j < splitting.get_nb_elem_local(1)));
      assert((k >= 0) && (k < splitting.get_nb_elem_local(2)));
    }
  else
    {
      assert((i >= -ghost_size) && (i < splitting.get_nb_elem_local(0) + ghost_size));
      assert((j >= -ghost_size) && (j < splitting.get_nb_elem_local(1) + ghost_size));
      assert((k >= -ghost_size) && (k < splitting.get_nb_elem_local(2) + ghost_size));
    }

  Int3 ijk = {i, j, k};
  return ijk;
}

inline int Cut_cell_FT_Disc::get_linear_index(int i, int j, int k, int ghost_size, const IJK_Splitting& splitting, bool expect_local_value)
{
  const int ni = splitting.get_nb_elem_local(0);
  const int nj = splitting.get_nb_elem_local(1);

  int offset = ghost_size + (ni + 2*ghost_size)*ghost_size + (ni + 2*ghost_size)*(nj + 2*ghost_size)*ghost_size;
  int linear_index = offset + i + (ni + 2*ghost_size)*j + (ni + 2*ghost_size)*(nj + 2*ghost_size)*k;

  if (expect_local_value)
    {
      assert((linear_index >= 0) && (linear_index <= offset + (ni - 1) + (ni + 2*ghost_size)*(nj - 1) + (ni + 2*ghost_size)*(nj + 2*ghost_size)*(splitting.get_nb_elem_local(2) - 1)));
    }
  else
    {
      assert((linear_index >= 0) && (linear_index <= offset + (ni - 1) + (ni + 2*ghost_size)*(nj - 1) + (ni + 2*ghost_size)*(nj + 2*ghost_size)*(splitting.get_nb_elem_local(2) + 2*ghost_size - 1)));
    }

  return linear_index;
}

inline Int3 Cut_cell_FT_Disc::get_ijk(int n) const
{
  int linear_index = linear_index_(n);
  Int3 ijk = get_ijk_from_linear_index(linear_index, ghost_size_, splitting_, false);

  return ijk;
}

inline int Cut_cell_FT_Disc::get_n(int i, int j, int k) const
{
  return indice_diphasique_(i,j,k);
}

inline int Cut_cell_FT_Disc::get_linear_index(int i, int j, int k) const
{
  return get_linear_index(i, j, k, ghost_size_, splitting_, false);
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

      int n_decale = get_n(i+di,j+dj,k+dk);
      return n_decale;
    }
}

#endif /* Cut_cell_FT_Disc_included */
