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

void Cut_cell_scalar::associer_persistant(Cut_cell_FT_Disc& cut_cell_disc)
{
  Cut_cell_data::associer_persistant(cut_cell_disc, 1);
}

void Cut_cell_scalar::associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc)
{
  Cut_cell_data::associer_ephemere(cut_cell_disc, 1);
}

void Cut_cell_scalar::associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc)
{
  Cut_cell_data::associer_paresseux(cut_cell_disc, 1);
}

void Cut_cell_scalar::set_valeur_cellules_diphasiques(double valeur)
{
  for (int n = 0; n < cut_cell_disc_->get_n_tot(); n++)
    {
      diph_l_(n) = valeur;
      diph_v_(n) = valeur;
    }
}


void Cut_cell_vector::associer_persistant(Cut_cell_FT_Disc& cut_cell_disc)
{
  Cut_cell_data::associer_persistant(cut_cell_disc, 3);
}

void Cut_cell_vector::associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc)
{
  Cut_cell_data::associer_ephemere(cut_cell_disc, 3);
}

void Cut_cell_vector::associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc)
{
  Cut_cell_data::associer_paresseux(cut_cell_disc, 3);
}

void Cut_cell_vector::set_valeur_cellules_diphasiques(double valeur)
{
  for (int n = 0; n < cut_cell_disc_->get_n_tot(); n++)
    {
      for (int dir = 0; dir < 3; dir++)
        {
          diph_l_(n,dir) = valeur;
          diph_v_(n,dir) = valeur;
        }
    }
}

Cut_field_scalar::Cut_field_scalar(IJK_Field_double& field) :
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

void Cut_field_scalar::remplir_cellules_devenant_diphasiques()
{
  int statut_diphasique = static_cast<int>(cut_cell_disc_->STATUT_DIPHASIQUE::NAISSANT);
  int index_min = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique);
  int index_max = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique+1);
  for (int index = index_min; index < index_max; index++)
    {
      int n = cut_cell_disc_->get_n_from_statut_diphasique_index(index);

      Int3 ijk = Cut_cell_FT_Disc::get_ijk_from_linear_index(cut_cell_disc_->get_linear_index(n), cut_cell_disc_->get_ghost_size(), pure_.get_splitting(), false);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      double old_indicatrice = cut_cell_disc_->get_interfaces().I(i,j,k);
      assert(cut_cell_disc_->get_interfaces().devient_diphasique(old_indicatrice, cut_cell_disc_->get_interfaces().In(i,j,k)));

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

void Cut_field_scalar::remplir_cellules_maintenant_pures()
{
  int statut_diphasique = static_cast<int>(cut_cell_disc_->STATUT_DIPHASIQUE::MOURRANT);
  int index_min = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique);
  int index_max = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique+1);
  for (int index = index_min; index < index_max; index++)
    {
      int n = cut_cell_disc_->get_n_from_statut_diphasique_index(index);

      Int3 ijk = Cut_cell_FT_Disc::get_ijk_from_linear_index(cut_cell_disc_->get_linear_index(n), cut_cell_disc_->get_ghost_size(), pure_.get_splitting(), false);
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

void Cut_field_scalar::transfert_diphasique_vers_pures()
{
  for (int n = 0; n < cut_cell_disc_->get_n_loc(); n++)
    {
      Int3 ijk = Cut_cell_FT_Disc::get_ijk_from_linear_index(cut_cell_disc_->get_linear_index(n), cut_cell_disc_->get_ghost_size(), pure_.get_splitting(), true);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      double indicatrice = cut_cell_disc_->get_interfaces().I(i,j,k);

      pure_(i,j,k) = indicatrice*diph_l_(n) + (1 - indicatrice)*diph_v_(n);
    }
  pure_.echange_espace_virtuel(pure_.ghost());
}

void Cut_field_scalar::set_field_data(const Nom& parser_expression_of_x_y_z_and_t, const IJK_Field_double& input_f, const double current_time)
{
  ArrOfDouble coord_i, coord_j, coord_k;
  build_local_coords(pure_, coord_i, coord_j, coord_k);

  const int ni = pure_.ni();
  const int nj = pure_.nj();
  const int nk = pure_.nk();

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
  pure_.echange_espace_virtuel(pure_.ghost());
}

void Cut_field_scalar::copy_from(Cut_field_scalar& data)
{
  const int ni = pure_.ni();
  const int nj = pure_.nj();
  const int nk = pure_.nk();
  const int ghost = pure_.ghost();
  assert(ni == data.pure_.ni());
  assert(nj == data.pure_.nj());
  assert(nk == data.pure_.nk());
  assert(data.pure_.ghost() >= ghost);
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

void Cut_field_scalar::add_from(Cut_field_scalar& data)
{
  const int ni = data.pure_.ni();
  const int nj = data.pure_.nj();
  const int nk = data.pure_.nk();
  const int ghost = data.pure_.ghost();
  assert(ni == pure_.ni());
  assert(nj == pure_.nj());
  assert(nk == pure_.nk());
  assert(ghost == pure_.ghost());
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

Cut_field_vector::Cut_field_vector(FixedVector<IJK_Field_double, 3>& field) :
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

void Cut_field_vector::remplir_cellules_devenant_diphasiques()
{
  int statut_diphasique = static_cast<int>(cut_cell_disc_->STATUT_DIPHASIQUE::NAISSANT);
  int index_min = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique);
  int index_max = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique+1);
  for (int index = index_min; index < index_max; index++)
    {
      int n = cut_cell_disc_->get_n_from_statut_diphasique_index(index);

      Int3 ijk = Cut_cell_FT_Disc::get_ijk_from_linear_index(cut_cell_disc_->get_linear_index(n), cut_cell_disc_->get_ghost_size(), pure_.get_splitting(), false);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      double old_indicatrice = cut_cell_disc_->get_interfaces().I(i,j,k);
      assert(cut_cell_disc_->get_interfaces().devient_diphasique(old_indicatrice, cut_cell_disc_->get_interfaces().In(i,j,k)));
      // On garde les donnees de l'ancienne phase pour la nouvelle cellule_diphasique
      int ancienne_phase = (int)old_indicatrice;
      if (ancienne_phase == 1)
        {
          for (int dir = 0; dir < 3; dir++)
            {
              diph_l_(n, dir) = pure_[dir](i,j,k);
            }
        }
      else if (ancienne_phase == 0)
        {
          for (int dir = 0; dir < 3; dir++)
            {
              diph_v_(n, dir) = pure_[dir](i,j,k);
            }
        }
    }
}

void Cut_field_vector::remplir_cellules_maintenant_pures()
{
  int statut_diphasique = static_cast<int>(cut_cell_disc_->STATUT_DIPHASIQUE::MOURRANT);
  int index_min = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique);
  int index_max = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique+1);
  for (int index = index_min; index < index_max; index++)
    {
      int n = cut_cell_disc_->get_n_from_statut_diphasique_index(index);

      Int3 ijk = Cut_cell_FT_Disc::get_ijk_from_linear_index(cut_cell_disc_->get_linear_index(n), cut_cell_disc_->get_ghost_size(), pure_.get_splitting(), false);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      double indicatrice = cut_cell_disc_->get_interfaces().I(i,j,k);
      assert(cut_cell_disc_->get_interfaces().est_pure(indicatrice));
      // On garde les donnees de la cellule diphasique pour la nouvelle cellule_pure
      int phase_pure = (int)indicatrice;
      if (phase_pure == 1)
        {
          for (int dir = 0; dir < 3; dir++)
            {
              pure_[dir](i,j,k) = diph_l_(n, dir);
            }
        }
      else if (phase_pure == 0)
        {
          for (int dir = 0; dir < 3; dir++)
            {
              pure_[dir](i,j,k) = diph_v_(n, dir);
            }
        }
    }
}

void Cut_field_vector::transfert_diphasique_vers_pures()
{
  for (int n = 0; n < cut_cell_disc_->get_n_loc(); n++)
    {
      Int3 ijk = Cut_cell_FT_Disc::get_ijk_from_linear_index(cut_cell_disc_->get_linear_index(n), cut_cell_disc_->get_ghost_size(), pure_.get_splitting(), true);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      double indicatrice = cut_cell_disc_->get_interfaces().I(i,j,k);

      for (int dir = 0; dir < 3; dir++)
        {
          pure_[dir](i,j,k) = indicatrice*diph_l_(n, dir) + (1 - indicatrice)*diph_v_(n, dir);
        }
    }
  pure_[0].echange_espace_virtuel(pure_[0].ghost());
  pure_[1].echange_espace_virtuel(pure_[1].ghost());
  pure_[2].echange_espace_virtuel(pure_[2].ghost());
}

void Cut_field_vector::set_to_sum(const Cut_field_vector& data_1, const Cut_field_vector& data_2)
{
  for (int dir = 0; dir < 3; dir++)
    {
      const int ni = data_1.pure_[dir].ni();
      const int nj = data_1.pure_[dir].nj();
      const int nk = data_1.pure_[dir].nk();
      const int ghost = data_1.pure_[dir].ghost();
      assert(ni == data_2.pure_[dir].ni());
      assert(nj == data_2.pure_[dir].nj());
      assert(nk == data_2.pure_[dir].nk());
      assert(ghost == data_2.pure_[dir].ghost());
      assert(ni == pure_[dir].ni());
      assert(nj == pure_[dir].nj());
      assert(nk == pure_[dir].nk());
      assert(ghost == pure_[dir].ghost());
      for (int k = -ghost; k < nk+ghost; k++)
        {
          for (int j = -ghost; j < nj+ghost; j++)
            {
              for (int i = -ghost; i < ni+ghost; i++)
                {
                  pure_[dir](i,j,k) = data_1.pure_[dir](i,j,k) + data_2.pure_[dir](i,j,k);
                }
            }
        }
    }
  for (int n = 0; n < cut_cell_disc_->get_n_tot(); n++)
    {
      for (int dir = 0; dir < 3; dir++)
        {
          diph_l_(n, dir) = data_1.diph_l_(n, dir) + data_2.diph_l_(n, dir);
          diph_v_(n, dir) = data_1.diph_v_(n, dir) + data_2.diph_v_(n, dir);
        }
    }
}



#endif
