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
// File      : Champ_diphasique.cpp
// Directory : $IJK_ROOT/src/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Champ_diphasique_TPP_included
#define Champ_diphasique_TPP_included

#include <Champ_diphasique.h>
#include <Cut_cell_FT_Disc.h>
#include <IJK_Field.h>

void IntTabFT_cut_cell::associer(Cut_cell_FT_Disc& cut_cell_disc, int dimension)
{
  cut_cell_disc_ = cut_cell_disc;
  cut_cell_disc_->add_to_int_data(*this, dimension);
}

void IntTabFT_cut_cell::associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc, int dimension)
{
  cut_cell_disc_ = cut_cell_disc;
  cut_cell_disc_->add_to_lazy_int_data(*this, dimension);
}

void IntTabFT_cut_cell::echange_espace_virtuel()
{
  cut_cell_disc_->get_desc_structure().echange_espace_virtuel(*this);
}

void DoubleTabFT_cut_cell::associer(Cut_cell_FT_Disc& cut_cell_disc, int dimension)
{
  cut_cell_disc_ = cut_cell_disc;
  cut_cell_disc_->add_to_data(*this, dimension);
}

void DoubleTabFT_cut_cell::associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc, int dimension)
{
  cut_cell_disc_ = cut_cell_disc;
  cut_cell_disc_->add_to_lazy_data(*this, dimension);
}

void DoubleTabFT_cut_cell::echange_espace_virtuel()
{
  cut_cell_disc_->get_desc_structure().echange_espace_virtuel(*this);
}

Cut_cell_data::Cut_cell_data(Cut_cell_FT_Disc& cut_cell_disc, int dimension)
{
  associer(cut_cell_disc, dimension);
}

void Cut_cell_data::associer(Cut_cell_FT_Disc& cut_cell_disc, int dimension)
{
  cut_cell_disc_ = cut_cell_disc;
  diph_l_.associer(cut_cell_disc_, dimension);
  diph_v_.associer(cut_cell_disc_, dimension);
}

void Cut_cell_data::echange_espace_virtuel()
{
  diph_l_.echange_espace_virtuel();
  diph_v_.echange_espace_virtuel();
}

Cut_cell_scalar::Cut_cell_scalar(Cut_cell_FT_Disc& cut_cell_disc) :
  Cut_cell_data(cut_cell_disc, 1)
{
}

void Cut_cell_scalar::associer(Cut_cell_FT_Disc& cut_cell_disc)
{
  Cut_cell_data::associer(cut_cell_disc, 1);
}

Cut_cell_vector::Cut_cell_vector(Cut_cell_FT_Disc& cut_cell_disc) :
  Cut_cell_data(cut_cell_disc, 3)
{
}

void Cut_cell_vector::associer(Cut_cell_FT_Disc& cut_cell_disc)
{
  Cut_cell_data::associer(cut_cell_disc, 3);
}

Cut_field_scalar::Cut_field_scalar(IJK_Field_double& field, Cut_cell_FT_Disc& cut_cell_disc) :
  Cut_cell_scalar(cut_cell_disc),
  pure_(field)
{
}

void Cut_field_scalar::echange_espace_virtuel(int le_ghost, double jump_i)
{
  Cut_cell_scalar::echange_espace_virtuel();
  pure_.echange_espace_virtuel(le_ghost, jump_i);
}

void Cut_field_scalar::remplir_cellules_diphasiques()
{
  for (int n = 0; n < cut_cell_disc_->get_n_loc(); n++)
    {
      Int3 ijk = Cut_cell_FT_Disc::get_ijk_from_linear_index(cut_cell_disc_->get_linear_index(n), cut_cell_disc_->get_ghost_size(), pure_.get_splitting(), true);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      diph_l_(n) = pure_(i,j,k);
      diph_v_(n) = pure_(i,j,k);
    }
  diph_l_.echange_espace_virtuel();
  diph_v_.echange_espace_virtuel();
}

Cut_field_vector::Cut_field_vector(FixedVector<IJK_Field_double, 3>& field, Cut_cell_FT_Disc& cut_cell_disc) :
  Cut_cell_vector(cut_cell_disc),
  pure_(field)
{
}

void Cut_field_vector::echange_espace_virtuel(int le_ghost, double jump_i)
{
  Cut_cell_vector::echange_espace_virtuel();
  pure_[0].echange_espace_virtuel(le_ghost, jump_i);
  pure_[1].echange_espace_virtuel(le_ghost, jump_i);
  pure_[2].echange_espace_virtuel(le_ghost, jump_i);
}

void Cut_field_vector::remplir_cellules_diphasiques()
{
  for (int n = 0; n < cut_cell_disc_->get_n_loc(); n++)
    {
      Int3 ijk = Cut_cell_FT_Disc::get_ijk_from_linear_index(cut_cell_disc_->get_linear_index(n), cut_cell_disc_->get_ghost_size(), pure_.get_splitting(), true);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      for (int dir = 0; dir < 3; dir++)
        {
          diph_l_(n, dir) = pure_[dir](i,j,k);
          diph_v_(n, dir) = pure_[dir](i,j,k);
        }
    }
  diph_l_.echange_espace_virtuel();
  diph_v_.echange_espace_virtuel();
}
#endif
