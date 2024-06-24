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
#include <IJK_Navier_Stokes_tools.h>
#include <IJK_Field.h>

void IntTabFT_cut_cell::associer_persistant(Cut_cell_FT_Disc& cut_cell_disc, int dimension)
{
  cut_cell_disc_ = cut_cell_disc;
  cut_cell_disc_->add_to_persistent_int_data(*this, dimension);
}

void IntTabFT_cut_cell::associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc, int dimension)
{
  cut_cell_disc_ = cut_cell_disc;
  cut_cell_disc_->add_to_transient_int_data(*this, dimension);
}

void IntTabFT_cut_cell::associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc, int dimension)
{
  cut_cell_disc_ = cut_cell_disc;
  cut_cell_disc_->add_to_lazy_int_data(*this, dimension);
}

struct struct_two_int_columns
{
  int index;
  int value;
};

struct struct_three_int_columns
{
  int index;
  int first_value;
  int second_value;
};

int compare_second_int_column(const void *a, const void *b)
{
  struct_two_int_columns *a1 = (struct_two_int_columns *)a;
  struct_two_int_columns *a2 = (struct_two_int_columns *)b;
  if ((*a1).value < (*a2).value)
    return -1;
  else if ((*a1).value > (*a2).value)
    return 1;
  else
    return 0;
}

int compare_second_then_third_int_column(const void *a, const void *b)
{
  struct_three_int_columns *a1 = (struct_three_int_columns *)a;
  struct_three_int_columns *a2 = (struct_three_int_columns *)b;
  if ((*a1).first_value < (*a2).first_value)
    {
      return -1;
    }
  else if ((*a1).first_value > (*a2).first_value)
    {
      return 1;
    }
  else
    {
      if ((*a1).second_value < (*a2).second_value)
        {
          return -1;
        }
      else if ((*a1).second_value > (*a2).second_value)
        {
          return 1;
        }
      else
        {
          return 0;
        }
    }
}

void IntTabFT_cut_cell::sort_tot(int column)
{
  if (dimension(1) == 2)
    {
      if (column == 1)
        {
          qsort(addr(), cut_cell_disc_->get_n_tot(), dimension(1)*sizeof(int), compare_second_int_column);
        }
      else
        {
          Cerr << "NotImplementedError: IntTabFT_cut_cell::sort_tot(int) with sorting other than the second column." << finl;
          Process::exit();
        }
    }
  else
    {
      Cerr << "NotImplementedError: IntTabFT_cut_cell::sort_tot(int) with other than 2 columns." << finl;
      Process::exit();
    }
}

void IntTabFT_cut_cell::sort_tot(int column_1, int column_2)
{
  if (dimension(1) == 3)
    {
      if (column_1 == 1 && column_2 == 2)
        {
          qsort(addr(), cut_cell_disc_->get_n_tot(), dimension(1)*sizeof(int), compare_second_then_third_int_column);
        }
      else
        {
          Cerr << "NotImplementedError: IntTabFT_cut_cell::sort_tot(int,int) with sorting other than the second, then third column." << finl;
          Process::exit();
        }
    }
  else
    {
      Cerr << "NotImplementedError: IntTabFT_cut_cell::sort_tot(int,int) with other than 3 columns." << finl;
      Process::exit();
    }
}

void IntTabFT_cut_cell::echange_espace_virtuel()
{
  cut_cell_disc_->get_desc_structure().echange_espace_virtuel(*this);
}

void IntTabFT_cut_cell::echange_espace_virtuel(MD_Vector_tools::Operations_echange op)
{
  cut_cell_disc_->get_desc_structure().echange_espace_virtuel(*this, op);
}

void DoubleTabFT_cut_cell::associer_persistant(Cut_cell_FT_Disc& cut_cell_disc, int dimension)
{
  cut_cell_disc_ = cut_cell_disc;
  cut_cell_disc_->add_to_persistent_double_data(*this, dimension);
}

void DoubleTabFT_cut_cell::associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc, int dimension)
{
  cut_cell_disc_ = cut_cell_disc;
  cut_cell_disc_->add_to_transient_double_data(*this, dimension);
}

void DoubleTabFT_cut_cell::associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc, int dimension)
{
  cut_cell_disc_ = cut_cell_disc;
  cut_cell_disc_->add_to_lazy_double_data(*this, dimension);
}

struct struct_two_double_columns
{
  double index;
  double value;
};

struct struct_three_double_columns
{
  double index;
  double first_value;
  double second_value;
};

int compare_second_double_column(const void *a, const void *b)
{
  struct_two_double_columns *a1 = (struct_two_double_columns *)a;
  struct_two_double_columns *a2 = (struct_two_double_columns *)b;
  if ((*a1).value < (*a2).value)
    return -1;
  else if ((*a1).value > (*a2).value)
    return 1;
  else
    return 0;
}

int compare_second_then_third_double_column(const void *a, const void *b)
{
  struct_three_double_columns *a1 = (struct_three_double_columns *)a;
  struct_three_double_columns *a2 = (struct_three_double_columns *)b;
  if ((*a1).first_value < (*a2).first_value)
    {
      return -1;
    }
  else if ((*a1).first_value > (*a2).first_value)
    {
      return 1;
    }
  else
    {
      if ((*a1).second_value < (*a2).second_value)
        {
          return -1;
        }
      else if ((*a1).second_value > (*a2).second_value)
        {
          return 1;
        }
      else
        {
          return 0;
        }
    }
}

void DoubleTabFT_cut_cell::sort_tot(int column)
{
  if (dimension(1) == 2)
    {
      if (column == 1)
        {
          qsort(addr(), cut_cell_disc_->get_n_tot(), dimension(1)*sizeof(double), compare_second_double_column);
        }
      else
        {
          Cerr << "NotImplementedError: DoubleTabFT_cut_cell::sort_tot(int) with sorting other than the second column." << finl;
          Process::exit();
        }
    }
  else
    {
      Cerr << "NotImplementedError: DoubleTabFT_cut_cell::sort_tot(int) with other than 2 columns." << finl;
      Process::exit();
    }
}

void DoubleTabFT_cut_cell::sort_tot(int column_1, int column_2)
{
  if (dimension(1) == 3)
    {
      if (column_1 == 1 && column_2 == 2)
        {
          qsort(addr(), cut_cell_disc_->get_n_tot(), dimension(1)*sizeof(double), compare_second_then_third_double_column);
        }
      else
        {
          Cerr << "NotImplementedError: DoubleTabFT_cut_cell::sort_tot(int,int) with sorting other than the second, then the third column." << finl;
          Process::exit();
        }
    }
  else
    {
      Cerr << "NotImplementedError: DoubleTabFT_cut_cell::sort_tot(int,int) with other than 3 columns." << finl;
      Process::exit();
    }
}

void DoubleTabFT_cut_cell::echange_espace_virtuel()
{
  cut_cell_disc_->get_desc_structure().echange_espace_virtuel(*this);
}

void DoubleTabFT_cut_cell::echange_espace_virtuel(MD_Vector_tools::Operations_echange op)
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

Cut_field_double::Cut_field_double()
{
}

void Cut_field_double::echange_espace_virtuel(int le_ghost)
{
  IJK_Field_double::echange_espace_virtuel(le_ghost);
  diph_l_.echange_espace_virtuel();
  diph_v_.echange_espace_virtuel();
}

void Cut_field_double::remplir_cellules_diphasiques()
{
  for (int n = 0; n < cut_cell_disc_->get_n_loc(); n++)
    {
      Int3 ijk = Cut_cell_FT_Disc::get_ijk_from_linear_index(cut_cell_disc_->get_linear_index(n), cut_cell_disc_->get_ghost_size(), IJK_Field_double::get_splitting(), true);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      diph_l_(n) = pure_(i,j,k);
      diph_v_(n) = pure_(i,j,k);
    }
  diph_l_.echange_espace_virtuel();
  diph_v_.echange_espace_virtuel();
}

void Cut_field_double::remplir_cellules_devenant_diphasiques()
{
  int statut_diphasique = static_cast<int>(cut_cell_disc_->STATUT_DIPHASIQUE::NAISSANT);
  int index_min = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique);
  int index_max = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique+1);
  for (int index = index_min; index < index_max; index++)
    {
      int n = cut_cell_disc_->get_n_from_statut_diphasique_index(index);

      Int3 ijk = Cut_cell_FT_Disc::get_ijk_from_linear_index(cut_cell_disc_->get_linear_index(n), cut_cell_disc_->get_ghost_size(), IJK_Field_double::get_splitting(), false);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      double old_indicatrice = cut_cell_disc_->get_interfaces().I(i,j,k);
      assert(cut_cell_disc_->get_interfaces().devient_diphasique(i,j,k));

      // On garde les donnees de l'ancienne phase pour la nouvelle cellule_diphasique
      int ancienne_phase = (int)old_indicatrice;
      if (ancienne_phase == 1)
        {
          diph_l_(n) = pure_(i,j,k);
        }
      else if (ancienne_phase == 0)
        {
          diph_v_(n) = pure_(i,j,k);
        }
    }
}

void Cut_field_double::remplir_cellules_maintenant_pures()
{
  int statut_diphasique = static_cast<int>(cut_cell_disc_->STATUT_DIPHASIQUE::MOURRANT);
  int index_min = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique);
  int index_max = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique+1);
  for (int index = index_min; index < index_max; index++)
    {
      int n = cut_cell_disc_->get_n_from_statut_diphasique_index(index);

      Int3 ijk = Cut_cell_FT_Disc::get_ijk_from_linear_index(cut_cell_disc_->get_linear_index(n), cut_cell_disc_->get_ghost_size(), IJK_Field_double::get_splitting(), false);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      double indicatrice = cut_cell_disc_->get_interfaces().In(i,j,k); // Note : In car on est avant l'inversion
      assert(cut_cell_disc_->get_interfaces().est_pure(indicatrice));
      // On garde les donnees de la cellule diphasique pour la nouvelle cellule_pure
      int phase_pure = (int)indicatrice;
      if (phase_pure == 1)
        {
          pure_(i,j,k) = diph_l_(n);
        }
      else if (phase_pure == 0)
        {
          pure_(i,j,k) = diph_v_(n);
        }
    }
}

void Cut_field_double::transfert_diphasique_vers_pures()
{
  for (int n = 0; n < cut_cell_disc_->get_n_loc(); n++)
    {
      Int3 ijk = Cut_cell_FT_Disc::get_ijk_from_linear_index(cut_cell_disc_->get_linear_index(n), cut_cell_disc_->get_ghost_size(), IJK_Field_double::get_splitting(), true);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      double indicatrice = cut_cell_disc_->get_interfaces().I(i,j,k);

      pure_(i,j,k) = indicatrice*diph_l_(n) + (1 - indicatrice)*diph_v_(n);
    }
  IJK_Field_double::echange_espace_virtuel(IJK_Field_double::ghost());
}

void Cut_field_double::set_field_data(const Nom& parser_expression_of_x_y_z_and_t, const IJK_Field_double& input_f, const double current_time)
{
  ArrOfDouble coord_i, coord_j, coord_k;
  build_local_coords(*this, coord_i, coord_j, coord_k);

  const int ni = IJK_Field_double::ni();
  const int nj = IJK_Field_double::nj();
  const int nk = IJK_Field_double::nk();

  std::string expr(parser_expression_of_x_y_z_and_t);
  Parser parser;
  parser.setString(expr);
  parser.setNbVar(5);
  parser.addVar("x");
  parser.addVar("y");
  parser.addVar("z");
  parser.addVar("t");
  parser.addVar("ff");
  parser.parseString();
  parser.setVar(3, current_time);
  for (int k = 0; k < nk; k++)
    {
      double z = coord_k[k];
      parser.setVar(2, z);
      for (int j = 0; j < nj; j++)
        {
          double y = coord_j[j];
          parser.setVar(1, y);
          for (int i = 0; i < ni; i++)
            {
              double x = coord_i[i];
              parser.setVar((int) 0, x);

              int n = cut_cell_disc_->get_n(i, j, k);
              if (n >= 0)
                {
                  double vol_l = cut_cell_disc_->get_interfaces().I()(i,j,k);

                  double dx = cut_cell_disc_->get_splitting().get_grid_geometry().get_constant_delta(0);
                  double dy = cut_cell_disc_->get_splitting().get_grid_geometry().get_constant_delta(1);
                  double dz = cut_cell_disc_->get_splitting().get_grid_geometry().get_constant_delta(2);

                  double bary_x = cut_cell_disc_->get_interfaces().get_barycentre_phase1_next()[0](i,j,k);
                  double bary_y = cut_cell_disc_->get_interfaces().get_barycentre_phase1_next()[1](i,j,k);
                  double bary_z = cut_cell_disc_->get_interfaces().get_barycentre_phase1_next()[2](i,j,k);

                  parser.setVar((int) 0, x + (bary_x - .5)*dx);
                  parser.setVar(1, y + (bary_y - .5)*dy);
                  parser.setVar(2, z + (bary_z - .5)*dz);
                  parser.setVar(4, 1.);
                  diph_l_(n) = (vol_l > 0)*parser.eval();

                  double opposing_bar_x = IJK_Interfaces::opposing_barycentre(bary_x, vol_l);
                  double opposing_bar_y = IJK_Interfaces::opposing_barycentre(bary_y, vol_l);
                  double opposing_bar_z = IJK_Interfaces::opposing_barycentre(bary_z, vol_l);

                  parser.setVar((int) 0, x + (opposing_bar_x - .5)*dx);
                  parser.setVar(1, y + (opposing_bar_y - .5)*dy);
                  parser.setVar(2, z + (opposing_bar_z - .5)*dz);
                  parser.setVar(4, 0.);
                  diph_v_(n) = (vol_l < 1)*parser.eval();

                  // Re-setting the Cartesian coordinates into the parser
                  parser.setVar((int) 0, x);
                  parser.setVar(1, y);
                  parser.setVar(2, z);
                }
              else
                {
                  parser.setVar(4, input_f(i, j, k));
                  pure_(i, j, k) = parser.eval();
                  assert((input_f(i, j, k) == 0.) || (input_f(i, j, k) == 1.));
                }
            }
        }
    }
  IJK_Field_double::echange_espace_virtuel(IJK_Field_double::ghost());
}

void Cut_field_double::set_to_uniform_value(double valeur)
{
  Cut_field_double::data() = valeur;
  for (int n = 0; n < cut_cell_disc_->get_n_tot(); n++)
    {
      diph_l_(n) = valeur;
      diph_v_(n) = valeur;
    }
}


void Cut_field_double::associer_persistant(Cut_cell_FT_Disc& cut_cell_disc)
{
  cut_cell_disc_ = cut_cell_disc;
  diph_l_.associer_persistant(cut_cell_disc_, 1);
  diph_v_.associer_persistant(cut_cell_disc_, 1);
}

void Cut_field_double::associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc)
{
  cut_cell_disc_ = cut_cell_disc;
  diph_l_.associer_ephemere(cut_cell_disc_, 1);
  diph_v_.associer_ephemere(cut_cell_disc_, 1);
}

void Cut_field_double::associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc)
{
  cut_cell_disc_ = cut_cell_disc;
  diph_l_.associer_paresseux(cut_cell_disc_, 1);
  diph_v_.associer_paresseux(cut_cell_disc_, 1);
}

void Cut_field_double::copy_from(Cut_field_double& data)
{
  const int ni = IJK_Field_double::ni();
  const int nj = IJK_Field_double::nj();
  const int nk = IJK_Field_double::nk();
  const int ghost = IJK_Field_double::ghost();
  assert(ni == data.ni());
  assert(nj == data.nj());
  assert(nk == data.nk());
  assert(data.ghost() >= ghost);
  for (int k = -ghost; k < nk+ghost; k++)
    {
      for (int j = -ghost; j < nj+ghost; j++)
        {
          for (int i = -ghost; i < ni+ghost; i++)
            {
              pure_(i,j,k) = data.pure_(i,j,k);
            }
        }
    }
  for (int n = 0; n < cut_cell_disc_->get_n_tot(); n++)
    {
      diph_l_(n) = data.diph_l_(n);
      diph_v_(n) = data.diph_v_(n);
    }
}

void Cut_field_double::add_from(Cut_field_double& data)
{
  const int ni = data.ni();
  const int nj = data.nj();
  const int nk = data.nk();
  const int ghost = data.ghost();
  assert(ni == IJK_Field_double::ni());
  assert(nj == IJK_Field_double::nj());
  assert(nk == IJK_Field_double::nk());
  assert(ghost == IJK_Field_double::ghost());
  for (int k = -ghost; k < nk+ghost; k++)
    {
      for (int j = -ghost; j < nj+ghost; j++)
        {
          for (int i = -ghost; i < ni+ghost; i++)
            {
              pure_(i,j,k) += data.pure_(i,j,k);
            }
        }
    }
  for (int n = 0; n < cut_cell_disc_->get_n_tot(); n++)
    {
      diph_l_(n) += data.diph_l_(n);
      diph_v_(n) += data.diph_v_(n);
    }
}

void Cut_field_double::set_to_sum(const Cut_field_double& data_1, const Cut_field_double& data_2)
{
  const int ni = data_1.ni();
  const int nj = data_1.nj();
  const int nk = data_1.nk();
  const int ghost = data_1.ghost();
  assert(ni == data_2.ni());
  assert(nj == data_2.nj());
  assert(nk == data_2.nk());
  assert(ghost == data_2.ghost());
  assert(ni == IJK_Field_double::ni());
  assert(nj == IJK_Field_double::nj());
  assert(nk == IJK_Field_double::nk());
  assert(ghost == IJK_Field_double::ghost());
  for (int k = -ghost; k < nk+ghost; k++)
    {
      for (int j = -ghost; j < nj+ghost; j++)
        {
          for (int i = -ghost; i < ni+ghost; i++)
            {
              pure_(i,j,k) = data_1.pure_(i,j,k) + data_2.pure_(i,j,k);
            }
        }
    }
  for (int n = 0; n < cut_cell_disc_->get_n_tot(); n++)
    {
      diph_l_(n) = data_1.diph_l_(n) + data_2.diph_l_(n);
      diph_v_(n) = data_1.diph_v_(n) + data_2.diph_v_(n);
    }
}

#endif
