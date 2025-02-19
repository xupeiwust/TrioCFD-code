/****************************************************************************
* Copyright (c) 2019, CEA
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

#ifndef Cut_field_included
#define Cut_field_included

#include <Objet_U.h>
#include <IJK_Field.h>
#include <IJK_Field_vector.h>
#include <TRUSTTabFT_cut_cell.h>

class Cut_cell_FT_Disc;

struct CutCell_GlobalInfo
{
  double overall;
  double overall_l;
  double overall_v;
  double pure;
  double diph_l;
  double diph_v;
  double diph_small;
  double diph_regular;
  double diph_nascent;
  double diph_dying;
};

/*! @brief : class Cut_field_template
 *
 */
template<typename _TYPE_, typename _TYPE_ARRAY_>
class Cut_field_template : public IJK_Field_template<_TYPE_,_TYPE_ARRAY_>
{
public :
  Cut_field_template();

  TRUSTTabFT_cut_cell<_TYPE_> diph_l_;
  TRUSTTabFT_cut_cell<_TYPE_> diph_v_;

  void echange_espace_virtuel(int ghost);
  void copie_pure_vers_diph_sans_interpolation();

  void echange_pure_vers_diph_cellules_initialement_pures();

  bool check_agreement_diph_pure_cellules_initialement_pures() const;
  bool check_agreement_diph_pure_cellules_finalement_pures() const;
  void echange_diph_vers_pure_cellules_finalement_pures();

  void vide_phase_invalide_cellules_diphasiques();

  bool check_agreement_tableau_pure_cellules_diphasiques(bool next_time) const;
  void remplir_tableau_pure_cellules_diphasiques(bool next_time);
  void remplir_tableau_pure_cellules_diphasiques_max(bool next_time);

  void set_field_data(const Nom& parser_expression_of_x_y_z_and_t, const IJK_Field_template<_TYPE_,_TYPE_ARRAY_>& input_f, const double current_time);
  void set_to_uniform_value(_TYPE_ valeur);
  void set_to_sum(const Cut_field_template<_TYPE_,_TYPE_ARRAY_>& data_1, const Cut_field_template<_TYPE_,_TYPE_ARRAY_>& data_2);

  void associer_persistant(Cut_cell_FT_Disc& cut_cell_disc);
  void associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc);
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc);

  void copy_from(const Cut_field_template<_TYPE_,_TYPE_ARRAY_>& data);
  void add_from(const Cut_field_template<_TYPE_,_TYPE_ARRAY_>& data, _TYPE_ constant = 1);

  const Cut_cell_FT_Disc& get_cut_cell_disc() const { return cut_cell_disc_.valeur(); }

  void allocate_persistant(Cut_cell_FT_Disc& cut_cell_disc, const Domaine_IJK& splitting, Domaine_IJK::Localisation loc, int ghost_size, int additional_k_layers = 0, int nb_compo = 1, bool external_storage = false, int monofluide=0, double rov=0., double rol=0., int use_inv_rho_in_pressure_solver=0)
  {
    IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::allocate(splitting, loc, ghost_size, additional_k_layers, nb_compo, external_storage, monofluide, rov, rol, use_inv_rho_in_pressure_solver);
    associer_persistant(cut_cell_disc);
  }

  void allocate_ephemere(Cut_cell_FT_Disc& cut_cell_disc, const Domaine_IJK& splitting, Domaine_IJK::Localisation loc, int ghost_size, int additional_k_layers = 0, int nb_compo = 1, bool external_storage = false, int monofluide=0, double rov=0., double rol=0., int use_inv_rho_in_pressure_solver=0)
  {
    IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::allocate(splitting, loc, ghost_size, additional_k_layers, nb_compo, external_storage, monofluide, rov, rol, use_inv_rho_in_pressure_solver);
    associer_ephemere(cut_cell_disc);
  }

  void allocate_paresseux(Cut_cell_FT_Disc& cut_cell_disc, const Domaine_IJK& splitting, Domaine_IJK::Localisation loc, int ghost_size, int additional_k_layers = 0, int nb_compo = 1, bool external_storage = false, int monofluide=0, double rov=0., double rol=0., int use_inv_rho_in_pressure_solver=0)
  {
    IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::allocate(splitting, loc, ghost_size, additional_k_layers, nb_compo, external_storage, monofluide, rov, rol, use_inv_rho_in_pressure_solver);
    associer_paresseux(cut_cell_disc);
  }

  _TYPE_& from_signed_independent_index(int signed_independent_index);
  const _TYPE_& from_signed_independent_index(int signed_independent_index) const;

  _TYPE_& from_ijk_and_phase(int i, int j, int k, bool phase);
  const _TYPE_& from_ijk_and_phase(int i, int j, int k, bool phase) const;

  _TYPE_& from_n_num_face_and_phase(int n, int num_face, bool phase);
  const _TYPE_& from_n_num_face_and_phase(int n, int num_face, bool phase) const;

  CutCell_GlobalInfo compute_global_energy_cut_cell(bool next, double constant_l, double constant_v) const;
  CutCell_GlobalInfo compute_d_global_energy_cut_cell(bool next) const;
  CutCell_GlobalInfo compute_norm_cut_cell(bool next) const;
  CutCell_GlobalInfo compute_min_cut_cell(bool next) const;
  CutCell_GlobalInfo compute_max_cut_cell(bool next) const;
  Nom get_value_location(_TYPE_ T) const;

  _TYPE_& pure_(int i, int j, int k)
  {
    return IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::operator()(i,j,k);
  }

  const _TYPE_& pure_(int i, int j, int k) const
  {
    return IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::operator()(i,j,k);
  }

  _TYPE_& operator()(int i, int j, int k)
  {
    Cerr << "Disabling operator() for the derived class Cut_field_template<_TYPE_,_TYPE_ARRAY_> of IJK_Field_template<_TYPE_,_TYPE_ARRAY_>. Please use pure_() instead." << finl;
    Process::exit();
    return IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::operator()(i,j,k);
  }

  const _TYPE_& operator()(int i, int j, int k) const
  {
    Cerr << "Disabling operator() for the derived class Cut_field_template<_TYPE_,_TYPE_ARRAY_> of IJK_Field_template<_TYPE_,_TYPE_ARRAY_>. Please use pure_() instead." << finl;
    Process::exit();
    return IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::operator()(i,j,k);
  }

  _TYPE_& operator()(int i, int j, int k, int compo)
  {
    Cerr << "Disabling operator() for the derived class Cut_field_template<_TYPE_,_TYPE_ARRAY_> of IJK_Field_template<_TYPE_,_TYPE_ARRAY_>. Please use pure_() instead." << finl;
    Process::exit();
    return IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::operator()(i,j,k,compo);
  }

  const _TYPE_& operator()(int i, int j, int k, int compo) const
  {
    Cerr << "Disabling operator() for the derived class Cut_field_template<_TYPE_,_TYPE_ARRAY_> of IJK_Field_template<_TYPE_,_TYPE_ARRAY_>. Please use pure_() instead." << finl;
    Process::exit();
    return IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::operator()(i,j,k,compo);
  }

  void multiply_by_scalar(_TYPE_ scalar_l, _TYPE_ scalar_v);
  void divide_by_scalar(_TYPE_ scalar_l, _TYPE_ scalar_v);

protected :
  OBS_PTR(Cut_cell_FT_Disc) cut_cell_disc_;
};

using Cut_field_int = Cut_field_template<int,ArrOfInt>;
using Cut_field_double = Cut_field_template<double,ArrOfDouble>;

/*! @brief : class Cut_field_vector
 *
 */
template<class T, int N>
class Cut_field_vector : public IJK_Field_vector<T, N>
{
public :
  Cut_field_vector();

  void set_to_uniform_value(int valeur);
  void echange_espace_virtuel(int ghost);
  int ghost();

  void echange_pure_vers_diph_cellules_initialement_pures();
  void echange_diph_vers_pure_cellules_finalement_pures();
  void vide_phase_invalide_cellules_diphasiques();

  Cut_field_template<T,TRUSTArray<T>>& operator[](int i)
  {
    assert(i>=0 && i<N);
    IJK_Field_template<T,TRUSTArray<T>>& field = IJK_Field_vector<T, N>::operator[](i);
    return static_cast<Cut_field_template<T,TRUSTArray<T>>&>(field);
  }
  const Cut_field_template<T,TRUSTArray<T>>& operator[](int i) const
  {
    assert(i>=0 && i<N);
    const IJK_Field_template<T,TRUSTArray<T>>& field = IJK_Field_vector<T, N>::operator[](i);
    return static_cast<const Cut_field_template<T,TRUSTArray<T>>&>(field);
  }

protected :
};

using Cut_field_vector3_int = Cut_field_vector<int, 3>;
using Cut_field_vector3_double = Cut_field_vector<double, 3>;

#endif /* Cut_field_included */

