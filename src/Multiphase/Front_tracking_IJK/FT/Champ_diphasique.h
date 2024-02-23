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
  void associer(Cut_cell_FT_Disc& cut_cell_disc, int dimension);
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc, int dimension);
  const Cut_cell_FT_Disc& get_cut_cell_disc() const { return cut_cell_disc_.valeur(); }

  void echange_espace_virtuel() override;

protected :
  REF(Cut_cell_FT_Disc) cut_cell_disc_;
  bool has_virtual_elements_;
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
  void associer(Cut_cell_FT_Disc& cut_cell_disc, int dimension);
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc, int dimension);
  const Cut_cell_FT_Disc& get_cut_cell_disc() const { return cut_cell_disc_.valeur(); }

  void echange_espace_virtuel() override;

protected :
  REF(Cut_cell_FT_Disc) cut_cell_disc_;
  bool has_virtual_elements_;
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
  void associer(Cut_cell_FT_Disc& cut_cell_disc) { IntTabFT_cut_cell::associer(cut_cell_disc, 1); }
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc) { IntTabFT_cut_cell::associer_paresseux(cut_cell_disc, 1); }

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
  void associer(Cut_cell_FT_Disc& cut_cell_disc) { IntTabFT_cut_cell::associer(cut_cell_disc, 3); }
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
  void associer(Cut_cell_FT_Disc& cut_cell_disc) { IntTabFT_cut_cell::associer(cut_cell_disc, 6); }
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
  void associer(Cut_cell_FT_Disc& cut_cell_disc) { DoubleTabFT_cut_cell::associer(cut_cell_disc, 1); }
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc) { DoubleTabFT_cut_cell::associer_paresseux(cut_cell_disc, 1); }

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
  void associer(Cut_cell_FT_Disc& cut_cell_disc) { DoubleTabFT_cut_cell::associer(cut_cell_disc, 3); }
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
  void associer(Cut_cell_FT_Disc& cut_cell_disc) { DoubleTabFT_cut_cell::associer(cut_cell_disc, 6); }
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
  Cut_cell_data() {}
  Cut_cell_data(Cut_cell_FT_Disc& cut_cell_disc, int dimension);

  void associer(Cut_cell_FT_Disc& cut_cell_disc, int dimension);
  const Cut_cell_FT_Disc& get_cut_cell_disc() { return cut_cell_disc_.valeur(); }

  DoubleTabFT_cut_cell diph_l_;
  DoubleTabFT_cut_cell diph_v_;

  void echange_espace_virtuel();

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
  Cut_cell_scalar() : Cut_cell_data() {}
  Cut_cell_scalar(Cut_cell_FT_Disc& cut_cell_disc);

  void associer(Cut_cell_FT_Disc& cut_cell_disc);

protected :
};

/*! @brief : class Cut_cell_vector
 *
 *  <Description of class Cut_cell_vector>
 *
 *
 *
 */
class Cut_cell_vector : public Cut_cell_data
{
public :
  Cut_cell_vector() : Cut_cell_data() {}
  Cut_cell_vector(Cut_cell_FT_Disc& cut_cell_disc);

  void associer(Cut_cell_FT_Disc& cut_cell_disc);

protected :
};

/*! @brief : class Cut_field_scalar
 *
 *  <Description of class Cut_field_scalar>
 *
 *
 *
 */
class Cut_field_scalar : public Cut_cell_scalar
{
public :
  Cut_field_scalar(IJK_Field_double& IJK_Field, Cut_cell_FT_Disc& cut_cell_disc);

  IJK_Field_double& pure_;

  void echange_espace_virtuel(int ghost, double Shear_DU=0.);
  void remplir_cellules_diphasiques();

protected :
};

/*! @brief : class Cut_field_vector
 *
 *  <Description of class Cut_field_vector>
 *
 *
 *
 */
class Cut_field_vector : public Cut_cell_vector
{
public :
  Cut_field_vector(FixedVector<IJK_Field_double, 3>& IJK_Field, Cut_cell_FT_Disc& cut_cell_disc);

  FixedVector<IJK_Field_double, 3>& pure_;

  void echange_espace_virtuel(int ghost, double Shear_DU=0.);
  void remplir_cellules_diphasiques();

protected :
};

#endif /* Champ_diphasique_included */
