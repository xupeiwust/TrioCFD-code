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
// File      : TRUSTTabFT_cut_cell.h
// Directory : $IJK_ROOT/src/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef TRUSTTabFT_cut_cell_included
#define TRUSTTabFT_cut_cell_included

#include <Objet_U.h>
#include <TRUSTTabFT.h>
#include <TRUST_Ref.h>
#include <IJK_Field_local_template.h>

class Cut_cell_FT_Disc;

/*! @brief : class TRUSTTabFT_cut_cell
 *
 */
template<typename _TYPE_>
class TRUSTTabFT_cut_cell : public TRUSTTabFT<_TYPE_>
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
  OBS_PTR(Cut_cell_FT_Disc) cut_cell_disc_;
};

/*! @brief : class TRUSTTabFT_cut_cell_scalar
 *
 */
template<typename _TYPE_>
class TRUSTTabFT_cut_cell_scalar : public TRUSTTabFT_cut_cell<_TYPE_>
{
public :
  void associer_persistant(Cut_cell_FT_Disc& cut_cell_disc) { TRUSTTabFT_cut_cell<_TYPE_>::associer_persistant(cut_cell_disc, 1); }
  void associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc) { TRUSTTabFT_cut_cell<_TYPE_>::associer_ephemere(cut_cell_disc, 1); }
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc) { TRUSTTabFT_cut_cell<_TYPE_>::associer_paresseux(cut_cell_disc, 1); }

protected :
};

/*! @brief : class TRUSTTabFT_cut_cell_vector2
 *
 */
template<typename _TYPE_>
class TRUSTTabFT_cut_cell_vector2 : public TRUSTTabFT_cut_cell<_TYPE_>
{
public :
  void associer_persistant(Cut_cell_FT_Disc& cut_cell_disc) { TRUSTTabFT_cut_cell<_TYPE_>::associer_persistant(cut_cell_disc, 2); }
  void associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc) { TRUSTTabFT_cut_cell<_TYPE_>::associer_ephemere(cut_cell_disc, 2); }
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc) { TRUSTTabFT_cut_cell<_TYPE_>::associer_paresseux(cut_cell_disc, 2); }

protected :
};


/*! @brief : class TRUSTTabFT_cut_cell_vector3
 *
 */
template<typename _TYPE_>
class TRUSTTabFT_cut_cell_vector3 : public TRUSTTabFT_cut_cell<_TYPE_>
{
public :
  void associer_persistant(Cut_cell_FT_Disc& cut_cell_disc) { TRUSTTabFT_cut_cell<_TYPE_>::associer_persistant(cut_cell_disc, 3); }
  void associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc) { TRUSTTabFT_cut_cell<_TYPE_>::associer_ephemere(cut_cell_disc, 3); }
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc) { TRUSTTabFT_cut_cell<_TYPE_>::associer_paresseux(cut_cell_disc, 3); }

protected :
};

/*! @brief : class TRUSTTabFT_cut_cell_vector6
 *
 */
template<typename _TYPE_>
class TRUSTTabFT_cut_cell_vector6 : public TRUSTTabFT_cut_cell<_TYPE_>
{
public :
  void associer_persistant(Cut_cell_FT_Disc& cut_cell_disc) { TRUSTTabFT_cut_cell<_TYPE_>::associer_persistant(cut_cell_disc, 6); }
  void associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc) { TRUSTTabFT_cut_cell<_TYPE_>::associer_ephemere(cut_cell_disc, 6); }
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc) { TRUSTTabFT_cut_cell<_TYPE_>::associer_paresseux(cut_cell_disc, 6); }

protected :
};

using IntTabFT_cut_cell = TRUSTTabFT_cut_cell<int>;
using DoubleTabFT_cut_cell = TRUSTTabFT_cut_cell<double>;

using IntTabFT_cut_cell_scalar = TRUSTTabFT_cut_cell_scalar<int>;
using DoubleTabFT_cut_cell_scalar = TRUSTTabFT_cut_cell_scalar<double>;

using IntTabFT_cut_cell_vector2 = TRUSTTabFT_cut_cell_vector2<int>;
using DoubleTabFT_cut_cell_vector2 = TRUSTTabFT_cut_cell_vector2<double>;

using IntTabFT_cut_cell_vector3 = TRUSTTabFT_cut_cell_vector3<int>;
using DoubleTabFT_cut_cell_vector3 = TRUSTTabFT_cut_cell_vector3<double>;

using IntTabFT_cut_cell_vector6 = TRUSTTabFT_cut_cell_vector6<int>;
using DoubleTabFT_cut_cell_vector6 = TRUSTTabFT_cut_cell_vector6<double>;

/*! @brief : class Cut_cell_data
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
  OBS_PTR(Cut_cell_FT_Disc) cut_cell_disc_;
};

/*! @brief : class Cut_cell_double
 *
 */
class Cut_cell_double : public Cut_cell_data
{
public :
  void associer_persistant(Cut_cell_FT_Disc& cut_cell_disc);
  void associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc);
  void associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc);

  bool verify_consistency_within_layer(int dir, int k_layer, const IJK_Field_local_double& flux);

protected :
};

#endif /* TRUSTTabFT_cut_cell_included */
