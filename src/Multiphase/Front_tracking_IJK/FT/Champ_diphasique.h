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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Champ_diphasique.h
// Directory : $IJK_ROOT/src/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Champ_diphasique_included
#define Champ_diphasique_included

#include <Objet_U.h>
#include <IJK_Field.h>
#include <FixedVector.h>
#include <TRUSTTabFT.h>

class Cut_cell_FT_Disc;

/*! @brief : class IntTabFT_cut_cell
 *
 *  <Description of class IntTabFT_cut_cell>
 *
 *
 *
 */
class IntTabFT_cut_cell : public IntTabFT
{
public :
  void associer_persistant(Cut_cell_FT_Disc& cut_cell_disc, int dimension);
  void associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc, int dimension);
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc, int dimension);
  const Cut_cell_FT_Disc& get_cut_cell_disc() const { return cut_cell_disc_.valeur(); }

  void sort_tot(int colum);
  void sort_tot(int column_1, int column_2);

  void echange_espace_virtuel() override;
  void echange_espace_virtuel(MD_Vector_tools::Operations_echange op);

protected :
  REF(Cut_cell_FT_Disc) cut_cell_disc_;
};

/*! @brief : class DoubleTabFT_cut_cell
 *
 *  <Description of class DoubleTabFT_cut_cell>
 *
 *
 *
 */
class DoubleTabFT_cut_cell : public DoubleTabFT
{
public :
  void associer_persistant(Cut_cell_FT_Disc& cut_cell_disc, int dimension);
  void associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc, int dimension);
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc, int dimension);
  const Cut_cell_FT_Disc& get_cut_cell_disc() const { return cut_cell_disc_.valeur(); }

  void sort_tot(int colum);
  void sort_tot(int column_1, int column_2);

  void echange_espace_virtuel() override;
  void echange_espace_virtuel(MD_Vector_tools::Operations_echange op);

protected :
  REF(Cut_cell_FT_Disc) cut_cell_disc_;
};

/*! @brief : class IntTabFT_cut_cell_scalar
 *
 *  <Description of class IntTabFT_cut_cell_scalar>
 *
 *
 *
 */
class IntTabFT_cut_cell_scalar : public IntTabFT_cut_cell
{
public :
  void associer_persistant(Cut_cell_FT_Disc& cut_cell_disc) { IntTabFT_cut_cell::associer_persistant(cut_cell_disc, 1); }
  void associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc) { IntTabFT_cut_cell::associer_ephemere(cut_cell_disc, 1); }
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc) { IntTabFT_cut_cell::associer_paresseux(cut_cell_disc, 1); }

protected :
};

/*! @brief : class IntTabFT_cut_cell_vector2
 *
 *  <Description of class IntTabFT_cut_cell_vector2>
 *
 *
 *
 */
class IntTabFT_cut_cell_vector2 : public IntTabFT_cut_cell
{
public :
  void associer_persistant(Cut_cell_FT_Disc& cut_cell_disc) { IntTabFT_cut_cell::associer_persistant(cut_cell_disc, 2); }
  void associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc) { IntTabFT_cut_cell::associer_ephemere(cut_cell_disc, 2); }
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc) { IntTabFT_cut_cell::associer_paresseux(cut_cell_disc, 2); }

protected :
};

/*! @brief : class IntTabFT_cut_cell_vector3
 *
 *  <Description of class IntTabFT_cut_cell_vector3>
 *
 *
 *
 */
class IntTabFT_cut_cell_vector3 : public IntTabFT_cut_cell
{
public :
  void associer_persistant(Cut_cell_FT_Disc& cut_cell_disc) { IntTabFT_cut_cell::associer_persistant(cut_cell_disc, 3); }
  void associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc) { IntTabFT_cut_cell::associer_ephemere(cut_cell_disc, 3); }
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc) { IntTabFT_cut_cell::associer_paresseux(cut_cell_disc, 3); }

protected :
};

/*! @brief : class IntTabFT_cut_cell_vector6
 *
 *  <Description of class IntTabFT_cut_cell_vector6>
 *
 *
 *
 */
class IntTabFT_cut_cell_vector6 : public IntTabFT_cut_cell
{
public :
  void associer_persistant(Cut_cell_FT_Disc& cut_cell_disc) { IntTabFT_cut_cell::associer_persistant(cut_cell_disc, 6); }
  void associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc) { IntTabFT_cut_cell::associer_ephemere(cut_cell_disc, 6); }
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc) { IntTabFT_cut_cell::associer_paresseux(cut_cell_disc, 6); }

protected :
};

/*! @brief : class DoubleTabFT_cut_cell_scalar
 *
 *  <Description of class DoubleTabFT_cut_cell_scalar>
 *
 *
 *
 */
class DoubleTabFT_cut_cell_scalar : public DoubleTabFT_cut_cell
{
public :
  void associer_persistant(Cut_cell_FT_Disc& cut_cell_disc) { DoubleTabFT_cut_cell::associer_persistant(cut_cell_disc, 1); }
  void associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc) { DoubleTabFT_cut_cell::associer_ephemere(cut_cell_disc, 1); }
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc) { DoubleTabFT_cut_cell::associer_paresseux(cut_cell_disc, 1); }

protected :
};

/*! @brief : class DoubleTabFT_cut_cell_vector2
 *
 *  <Description of class DoubleTabFT_cut_cell_vector2>
 *
 *
 *
 */
class DoubleTabFT_cut_cell_vector2 : public DoubleTabFT_cut_cell
{
public :
  void associer_persistant(Cut_cell_FT_Disc& cut_cell_disc) { DoubleTabFT_cut_cell::associer_persistant(cut_cell_disc, 2); }
  void associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc) { DoubleTabFT_cut_cell::associer_ephemere(cut_cell_disc, 2); }
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc) { DoubleTabFT_cut_cell::associer_paresseux(cut_cell_disc, 2); }

protected :
};


/*! @brief : class DoubleTabFT_cut_cell_vector3
 *
 *  <Description of class DoubleTabFT_cut_cell_vector3>
 *
 *
 *
 */
class DoubleTabFT_cut_cell_vector3 : public DoubleTabFT_cut_cell
{
public :
  void associer_persistant(Cut_cell_FT_Disc& cut_cell_disc) { DoubleTabFT_cut_cell::associer_persistant(cut_cell_disc, 3); }
  void associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc) { DoubleTabFT_cut_cell::associer_ephemere(cut_cell_disc, 3); }
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc) { DoubleTabFT_cut_cell::associer_paresseux(cut_cell_disc, 3); }

protected :
};

/*! @brief : class DoubleTabFT_cut_cell_vector6
 *
 *  <Description of class DoubleTabFT_cut_cell_vector6>
 *
 *
 *
 */
class DoubleTabFT_cut_cell_vector6 : public DoubleTabFT_cut_cell
{
public :
  void associer_persistant(Cut_cell_FT_Disc& cut_cell_disc) { DoubleTabFT_cut_cell::associer_persistant(cut_cell_disc, 6); }
  void associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc) { DoubleTabFT_cut_cell::associer_ephemere(cut_cell_disc, 6); }
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc) { DoubleTabFT_cut_cell::associer_paresseux(cut_cell_disc, 6); }

protected :
};

/*! @brief : class Cut_cell_data
 *
 *  <Description of class Cut_cell_data>
 *
 *
 *
 */
class Cut_cell_data
{
public :
  void associer_persistant(Cut_cell_FT_Disc& cut_cell_disc, int dimension);
  void associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc, int dimension);
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc, int dimension);
  const Cut_cell_FT_Disc& get_cut_cell_disc() const { return cut_cell_disc_.valeur(); }

  DoubleTabFT_cut_cell diph_l_;
  DoubleTabFT_cut_cell diph_v_;

  void echange_espace_virtuel();
  void echange_espace_virtuel(MD_Vector_tools::Operations_echange op);

protected :
  REF(Cut_cell_FT_Disc) cut_cell_disc_;
};

/*! @brief : class Cut_cell_scalar
 *
 *  <Description of class Cut_cell_scalar>
 *
 *
 *
 */
class Cut_cell_scalar : public Cut_cell_data
{
public :
  void associer_persistant(Cut_cell_FT_Disc& cut_cell_disc);
  void associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc);
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc);

protected :
};

/*! @brief : class Cut_field_scalar
 *
 *  <Description of class Cut_field_scalar>
 *
 *
 *
 */
class Cut_field_scalar : public IJK_Field_double
{
public :
  Cut_field_scalar();

  DoubleTabFT_cut_cell diph_l_;
  DoubleTabFT_cut_cell diph_v_;

  void echange_espace_virtuel(int ghost);
  void remplir_cellules_diphasiques();
  void remplir_cellules_devenant_diphasiques();
  void remplir_cellules_maintenant_pures();
  void transfert_diphasique_vers_pures();
  void set_field_data(const Nom& parser_expression_of_x_y_z_and_t, const IJK_Field_double& input_f, const double current_time);
  void set_to_uniform_value(double valeur);
  void set_to_sum(const Cut_field_scalar& data_1, const Cut_field_scalar& data_2);

  void associer_persistant(Cut_cell_FT_Disc& cut_cell_disc);
  void associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc);
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc);

  void copy_from(Cut_field_scalar& data);
  void add_from(Cut_field_scalar& data);

  const Cut_cell_FT_Disc& get_cut_cell_disc() const { return cut_cell_disc_.valeur(); }

  // :integration(Dorian) De maniere temporaire avant la modification du type des champs
  // vectoriels IJK, on ajoute des routines a la classe Cut_field_scalar,
  // qui permettent a l'objet de referer aux donnees d'un objet IJK_Field_double existant.
  // Ces fonctionnalites sont temporaires et devront etre supprimees.
  void set_ijk_field(IJK_Field_double& field)
  {
    data_.ref_array(field.data());
    //data_ = field.data();
    ni_ = field.ni();
    nj_ = field.nj();
    nk_ = field.nk();
    nb_compo_ = field.nb_compo();
    j_stride_ = field.j_stride();
    compo_stride_ = field.compo_stride();
    ghost_size_ = field.ghost();
    k_layer_shift_ = field.k_shift();
    additional_k_layers_ = field.k_shift_max();
    allocated_size_ = field.get_allocated_size();
    offset_ = (int)(((long int)field.k_layer(0) - (long int)field.data().addr())/8);
    splitting_ref_ = field.get_splitting();
    localisation_ = field.get_localisation();
  }

  double& pure_(int i, int j, int k)
  {
    return IJK_Field_double::operator()(i,j,k);
  }

  const double& pure_(int i, int j, int k) const
  {
    return IJK_Field_double::operator()(i,j,k);
  }

  double& operator()(int i, int j, int k)
  {
    Cerr << "Disabling operator() for the derived class Cut_field_scalar of IJK_Field_double." << finl;
    Process::exit();
    return IJK_Field_double::operator()(i,j,k);
  }

  const double& operator()(int i, int j, int k) const
  {
    Cerr << "Disabling operator() for the derived class Cut_field_scalar of IJK_Field_double." << finl;
    Process::exit();
    return IJK_Field_double::operator()(i,j,k);
  }

  double& operator()(int i, int j, int k, int compo)
  {
    Cerr << "Disabling operator() for the derived class Cut_field_scalar of IJK_Field_double." << finl;
    Process::exit();
    return IJK_Field_double::operator()(i,j,k,compo);
  }

  const double& operator()(int i, int j, int k, int compo) const
  {
    Cerr << "Disabling operator() for the derived class Cut_field_scalar of IJK_Field_double." << finl;
    Process::exit();
    return IJK_Field_double::operator()(i,j,k,compo);
  }

protected :
  REF(Cut_cell_FT_Disc) cut_cell_disc_;
};

#endif /* Champ_diphasique_included */
