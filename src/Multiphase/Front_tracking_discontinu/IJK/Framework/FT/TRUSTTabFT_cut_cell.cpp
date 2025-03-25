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

#include <Cut_cell_FT_Disc.h>
#include <TRUSTTabFT_cut_cell.h>
#include <IJK_Interfaces.h>
#include <IJK_Navier_Stokes_tools.h>

template <typename _TYPE_>
void compare_second_column(_TYPE_* ptr, int sz)
{
  using pair = std::array<_TYPE_, 2>;
  pair* tmp = reinterpret_cast<pair*>(ptr);
  std::sort(tmp, tmp+sz, [&](const pair& p1, const pair& p2)
  {
    return (p1[1] < p2[1]);
  });
}

template <typename _TYPE_>
void compare_second_then_third_column(_TYPE_* ptr, int sz)
{
  using triplet = std::array<_TYPE_, 3>;
  triplet* tmp = reinterpret_cast<triplet*>(ptr);
  std::sort(tmp, tmp+sz, [&](const triplet& p1, const triplet& p2)
  {
    if(p1[1] == p2[1])
      return (p1[2] < p2[2]);
    else
      return (p1[1] < p2[1]);
  });
}

template<typename _TYPE_>
void TRUSTTabFT_cut_cell<_TYPE_>::associer_persistant(Cut_cell_FT_Disc& cut_cell_disc, int dimension)
{
  cut_cell_disc_ = cut_cell_disc;
  cut_cell_disc_->add_to_persistent_data(*this, dimension);
}

template<typename _TYPE_>
void TRUSTTabFT_cut_cell<_TYPE_>::associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc, int dimension)
{
  cut_cell_disc_ = cut_cell_disc;
  cut_cell_disc_->add_to_transient_data(*this, dimension);
}

template<typename _TYPE_>
void TRUSTTabFT_cut_cell<_TYPE_>::associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc, int dimension)
{
  cut_cell_disc_ = cut_cell_disc;
  cut_cell_disc_->add_to_lazy_data(*this, dimension);
}

template<>
void TRUSTTabFT_cut_cell<int>::sort_tot(int column)
{
  if (dimension(1) == 2)
    {
      if (column == 1)
        {
          compare_second_column<int>(addr(), cut_cell_disc_->get_n_tot());
        }
      else
        {
          Cerr << "NotImplementedError: TRUSTTabFT_cut_cell<int>::sort_tot(int) with sorting other than the second column." << finl;
          Process::exit();
        }
    }
  else
    {
      Cerr << "NotImplementedError: TRUSTTabFT_cut_cell<int>::sort_tot(int) with other than 2 columns." << finl;
      Process::exit();
    }
}

template<>
void TRUSTTabFT_cut_cell<double>::sort_tot(int column)
{
  if (dimension(1) == 2)
    {
      if (column == 1)
        {
          compare_second_column<double>(addr(), cut_cell_disc_->get_n_tot());
        }
      else
        {
          Cerr << "NotImplementedError: TRUSTTabFT_cut_cell<double>::sort_tot(double) with sorting other than the second column." << finl;
          Process::exit();
        }
    }
  else
    {
      Cerr << "NotImplementedError: TRUSTTabFT_cut_cell<double>::sort_tot(double) with other than 2 columns." << finl;
      Process::exit();
    }
}

template<>
void TRUSTTabFT_cut_cell<int>::sort_tot(int column_1, int column_2)
{
  if (dimension(1) == 3)
    {
      if (column_1 == 1 && column_2 == 2)
        {
          compare_second_then_third_column<int>(addr(), cut_cell_disc_->get_n_tot());
        }
      else
        {
          Cerr << "NotImplementedError: TRUSTTabFT_cut_cell<int>::sort_tot(int,int) with sorting other than the second, then third column." << finl;
          Process::exit();
        }
    }
  else
    {
      Cerr << "NotImplementedError: TRUSTTabFT_cut_cell<int>::sort_tot(int,int) with other than 3 columns." << finl;
      Process::exit();
    }
}

template<>
void TRUSTTabFT_cut_cell<double>::sort_tot(int column_1, int column_2)
{
  if (dimension(1) == 3)
    {
      if (column_1 == 1 && column_2 == 2)
        {
          compare_second_then_third_column<double>(addr(), cut_cell_disc_->get_n_tot());
        }
      else
        {
          Cerr << "NotImplementedError: TRUSTTabFT_cut_cell<double>::sort_tot(double,double) with sorting other than the second, then third column." << finl;
          Process::exit();
        }
    }
  else
    {
      Cerr << "NotImplementedError: TRUSTTabFT_cut_cell<double>::sort_tot(double,double) with other than 3 columns." << finl;
      Process::exit();
    }
}

template<typename _TYPE_>
void TRUSTTabFT_cut_cell<_TYPE_>::echange_espace_virtuel()
{
  cut_cell_disc_->get_desc_structure().echange_espace_virtuel(*this);
}

template<typename _TYPE_>
void TRUSTTabFT_cut_cell<_TYPE_>::echange_espace_virtuel(MD_Vector_tools::Operations_echange op)
{
  cut_cell_disc_->get_desc_structure().echange_espace_virtuel(*this, op);
}


void Cut_cell_data::associer_persistant(Cut_cell_FT_Disc& cut_cell_disc, int dimension)
{
  cut_cell_disc_ = cut_cell_disc;
  diph_l_.associer_persistant(cut_cell_disc_, dimension);
  diph_v_.associer_persistant(cut_cell_disc_, dimension);
}

void Cut_cell_data::associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc, int dimension)
{
  cut_cell_disc_ = cut_cell_disc;
  diph_l_.associer_ephemere(cut_cell_disc_, dimension);
  diph_v_.associer_ephemere(cut_cell_disc_, dimension);
}

void Cut_cell_data::associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc, int dimension)
{
  cut_cell_disc_ = cut_cell_disc;
  diph_l_.associer_paresseux(cut_cell_disc_, dimension);
  diph_v_.associer_paresseux(cut_cell_disc_, dimension);
}

void Cut_cell_data::echange_espace_virtuel()
{
  diph_l_.echange_espace_virtuel();
  diph_v_.echange_espace_virtuel();
}

void Cut_cell_data::echange_espace_virtuel(MD_Vector_tools::Operations_echange op)
{
  diph_l_.echange_espace_virtuel(op);
  diph_v_.echange_espace_virtuel(op);
}

void Cut_cell_double::associer_persistant(Cut_cell_FT_Disc& cut_cell_disc)
{
  Cut_cell_data::associer_persistant(cut_cell_disc, 1);
}

void Cut_cell_double::associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc)
{
  Cut_cell_data::associer_ephemere(cut_cell_disc, 1);
}

void Cut_cell_double::associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc)
{
  Cut_cell_data::associer_paresseux(cut_cell_disc, 1);
}

// Verifie que le champ pure_, suppose etre localise sur une surface, est coherent avec
// le champ diph_v_ ou diph_l_ si la dite surface appartient a une cellule pure.
bool Cut_cell_double::verify_consistency_within_layer(int dir, int k_layer, const IJK_Field_local_double& flux)
{
  int ni = (dir == 0) ? flux.ni() : flux.ni() - 1;
  int nj = (dir == 1) ? flux.nj() : flux.nj() - 1;
  int di = (dir == 0)*(-1);
  int dj = (dir == 1)*(-1);
  int dk = (dir == 2)*(-1);
  for (int j = 0; j < nj; j++)
    {
      for (int i = 0; i < ni; i++)
        {
          int n = cut_cell_disc_->get_n(i, j, k_layer);
          if (n >= 0)
            {
              int n_decale = cut_cell_disc_->get_n(i+di, j+dj, k_layer+dk);
              if (n_decale < 0)
                {
                  // Surface between a pure and non-pure cells
                  // There should be a consistency
                  int phase = IJK_Interfaces::convert_indicatrice_to_phase(cut_cell_disc_->indic_pure(i+di, j+dj, k_layer+dk));

                  if (phase == 0)
                    {
                      if (flux(i,j,0) != diph_v_(n))
                        {
                          assert(0);
                          return 0;
                        }
                    }
                  else
                    {
                      if (flux(i,j,0) != diph_l_(n))
                        {
                          assert(0);
                          return 0;
                        }
                    }
                }
            }
        }
    }

  return 1;
}


// Explicit instantiations

template class TRUSTTabFT_cut_cell<int>;
template class TRUSTTabFT_cut_cell<double>;

