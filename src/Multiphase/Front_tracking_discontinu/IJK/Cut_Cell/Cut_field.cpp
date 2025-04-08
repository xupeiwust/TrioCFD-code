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

#include <Cut_field.h>
#include <IJK_Interfaces.h>
#include <Cut_cell_FT_Disc.h>
#include <IJK_Navier_Stokes_tools.h>
#include <IJK_Field.h>

template<typename _TYPE_, typename _TYPE_ARRAY_>
Cut_field_template<_TYPE_,_TYPE_ARRAY_>::Cut_field_template()
{
}


template<typename _TYPE_, typename _TYPE_ARRAY_>
void Cut_field_template<_TYPE_,_TYPE_ARRAY_>::echange_espace_virtuel(int le_ghost)
{
  IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::echange_espace_virtuel(le_ghost);
  diph_l_.echange_espace_virtuel();
  diph_v_.echange_espace_virtuel();
}

template<typename _TYPE_, typename _TYPE_ARRAY_>
void Cut_field_template<_TYPE_,_TYPE_ARRAY_>::copie_pure_vers_diph_sans_interpolation()
{
  for (int n = 0; n < cut_cell_disc_->get_n_loc(); n++)
    {
      Int3 ijk = cut_cell_disc_->get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      diph_l_(n) = pure_(i,j,k);
      diph_v_(n) = pure_(i,j,k);
    }
  diph_l_.echange_espace_virtuel();
  diph_v_.echange_espace_virtuel();
}

template<typename _TYPE_, typename _TYPE_ARRAY_>
void Cut_field_template<_TYPE_,_TYPE_ARRAY_>::echange_pure_vers_diph_cellules_initialement_pures()
{
  if (cut_cell_disc_.est_nul())
    {
      // Si le champ n'est pas initialise, c'est normal de ne pas avoir de cut_cell_disc et on ne fait rien
      assert((IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ni() == 0));
      assert((IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::nj() == 0));
      assert((IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::nk() == 0));
    }
  else
    {
      IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::echange_espace_virtuel(IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ghost());

      int statut_diphasique_monophasique = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::MONOPHASIQUE);
      int statut_diphasique_naissant = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::NAISSANT);
      assert(statut_diphasique_naissant == statut_diphasique_monophasique + 1);
      int index_min = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique_monophasique);
      int index_max = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique_naissant+1);
      for (int index = index_min; index < index_max; index++)
        {
          int n = cut_cell_disc_->get_n_from_statut_diphasique_index(index);

          Int3 ijk = cut_cell_disc_->get_ijk(n);
          int i = ijk[0];
          int j = ijk[1];
          int k = ijk[2];

          if (!cut_cell_disc_->get_domaine().within_ghost(i, j, k, IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ghost(), IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ghost()))
            continue;

          double old_indicatrice = cut_cell_disc_->get_interfaces().I(i,j,k);

          // On garde les donnees de l'ancienne phase pour la nouvelle cellule_diphasique
          int ancienne_phase = IJK_Interfaces::convert_indicatrice_to_phase(old_indicatrice);
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
}

template<typename _TYPE_, typename _TYPE_ARRAY_>
bool Cut_field_template<_TYPE_,_TYPE_ARRAY_>::check_agreement_diph_pure_cellules_initialement_pures() const
{
  int statut_diphasique_monophasique = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::MONOPHASIQUE);
  int statut_diphasique_naissant = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::NAISSANT);
  assert(statut_diphasique_naissant == statut_diphasique_monophasique + 1);
  int index_min = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique_monophasique);
  int index_max = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique_naissant+1);
  for (int index = index_min; index < index_max; index++)
    {
      int n = cut_cell_disc_->get_n_from_statut_diphasique_index(index);

      Int3 ijk = cut_cell_disc_->get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      if (!cut_cell_disc_->get_domaine().within_ghost(i, j, k, IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ghost(), IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ghost()))
        continue;

      double indicatrice = cut_cell_disc_->get_interfaces().I(i,j,k);
      assert(cut_cell_disc_->get_interfaces().est_pure(indicatrice));
      // On garde les donnees de la cellule diphasique pour pure
      int phase_pure = IJK_Interfaces::convert_indicatrice_to_phase(indicatrice);
      double diph_value = (phase_pure == 0) ? diph_v_(n) : diph_l_(n);
      if (diph_value != pure_(i,j,k))
        {
          assert(0);
          return 0;
        }
    }

  return 1;
}

template<typename _TYPE_, typename _TYPE_ARRAY_>
bool Cut_field_template<_TYPE_,_TYPE_ARRAY_>::check_agreement_diph_pure_cellules_finalement_pures() const
{
  int statut_diphasique_mourrant = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::MOURRANT);
  int statut_diphasique_monophasique = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::MONOPHASIQUE);
  assert(statut_diphasique_monophasique == statut_diphasique_mourrant + 1);
  int index_min = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique_mourrant);
  int index_max = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique_monophasique+1);
  for (int index = index_min; index < index_max; index++)
    {
      int n = cut_cell_disc_->get_n_from_statut_diphasique_index(index);

      Int3 ijk = cut_cell_disc_->get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      if (!cut_cell_disc_->get_domaine().within_ghost(i, j, k, IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ghost(), IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ghost()))
        continue;

      double indicatrice = cut_cell_disc_->get_interfaces().In(i,j,k); // Note : In car on est avant l'inversion
      assert(cut_cell_disc_->get_interfaces().est_pure(indicatrice));
      // On garde les donnees de la cellule diphasique pour pure
      int phase_pure = IJK_Interfaces::convert_indicatrice_to_phase(indicatrice);
      double diph_value = (phase_pure == 0) ? diph_v_(n) : diph_l_(n);
      if (diph_value != pure_(i,j,k))
        {
          assert(0);
          return 0;
        }
    }

  return 1;
}

template<typename _TYPE_, typename _TYPE_ARRAY_>
void Cut_field_template<_TYPE_,_TYPE_ARRAY_>::echange_diph_vers_pure_cellules_finalement_pures()
{
  if (cut_cell_disc_.est_nul())
    {
      // Si le champ n'est pas initialise, c'est normal de ne pas avoir de cut_cell_disc et on ne fait rien
      assert((IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ni() == 0));
      assert((IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::nj() == 0));
      assert((IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::nk() == 0));
    }
  else
    {
      int statut_diphasique_mourrant = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::MOURRANT);
      int statut_diphasique_monophasique = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::MONOPHASIQUE);
      assert(statut_diphasique_monophasique == statut_diphasique_mourrant + 1);
      int index_min = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique_mourrant);
      int index_max = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique_monophasique+1);
      for (int index = index_min; index < index_max; index++)
        {
          int n = cut_cell_disc_->get_n_from_statut_diphasique_index(index);

          Int3 ijk = cut_cell_disc_->get_ijk(n);
          int i = ijk[0];
          int j = ijk[1];
          int k = ijk[2];

          if (!cut_cell_disc_->get_domaine().within_ghost(i, j, k, IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ghost(), IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ghost()))
            continue;

          double indicatrice = cut_cell_disc_->get_interfaces().In(i,j,k);
          assert(cut_cell_disc_->get_interfaces().est_pure(indicatrice));
          // On garde les donnees de la cellule diphasique pour pure
          int phase_pure = IJK_Interfaces::convert_indicatrice_to_phase(indicatrice);
          _TYPE_ diph_value = (phase_pure == 0) ? diph_v_(n) : diph_l_(n);
          if (diph_value != pure_(i,j,k))
            {
              pure_(i,j,k) = diph_value;
            }
        }
    }
}

template<typename _TYPE_, typename _TYPE_ARRAY_>
void Cut_field_template<_TYPE_,_TYPE_ARRAY_>::vide_phase_invalide_cellules_diphasiques()
{
  if (cut_cell_disc_.est_nul())
    {
      // Si le champ n'est pas initialise, c'est normal de ne pas avoir de cut_cell_disc et on ne fait rien
      assert((IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ni() == 0));
      assert((IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::nj() == 0));
      assert((IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::nk() == 0));
    }
  else
    {
      int statut_diphasique_mourrant = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::MOURRANT);
      int statut_diphasique_monophasique = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::MONOPHASIQUE);
      assert(statut_diphasique_monophasique == statut_diphasique_mourrant + 1);
      int index_min = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique_mourrant);
      int index_max = cut_cell_disc_->get_statut_diphasique_value_index(statut_diphasique_monophasique+1);
      for (int index = index_min; index < index_max; index++)
        {
          int n = cut_cell_disc_->get_n_from_statut_diphasique_index(index);

          Int3 ijk = cut_cell_disc_->get_ijk(n);
          int i = ijk[0];
          int j = ijk[1];
          int k = ijk[2];

          if (!cut_cell_disc_->get_domaine().within_ghost(i, j, k, IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ghost(), IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ghost()))
            continue;

          double next_indicatrice = cut_cell_disc_->get_interfaces().In(i,j,k);

          int phase_valide = IJK_Interfaces::convert_indicatrice_to_phase(next_indicatrice);
          if (phase_valide == 1)
            {
              //assert(std::abs(diph_v_(n)) < 1e-12);
              diph_v_(n) = 0.;
            }
          else if (phase_valide == 0)
            {
              //assert(std::abs(diph_l_(n)) < 1e-12);
              diph_l_(n) = 0.;
            }
        }
    }
}

template<>
bool Cut_field_template<double,ArrOfDouble>::check_agreement_tableau_pure_cellules_diphasiques(bool next_time) const
{
  for (int n = 0; n < cut_cell_disc_->get_n_loc(); n++)
    {
      Int3 ijk = cut_cell_disc_->get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      double indicatrice = next_time ? cut_cell_disc_->get_interfaces().In(i,j,k) : cut_cell_disc_->get_interfaces().I(i,j,k);

      if (diph_l_(n) == diph_v_(n))
        {
          if (pure_(i,j,k) != diph_l_(n))
            {
              assert(0);
              return 0;
            }
        }
      else
        {
          if (pure_(i,j,k) != indicatrice*diph_l_(n) + (1 - indicatrice)*diph_v_(n))
            {
              assert(0);
              return 0;
            }
        }
    }

  return 1;
}

template<>
void Cut_field_template<double,ArrOfDouble>::remplir_tableau_pure_cellules_diphasiques(bool next_time)
{
  for (int n = 0; n < cut_cell_disc_->get_n_loc(); n++)
    {
      Int3 ijk = cut_cell_disc_->get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      double indicatrice = next_time ? cut_cell_disc_->get_interfaces().In(i,j,k) : cut_cell_disc_->get_interfaces().I(i,j,k);

      if (diph_l_(n) == diph_v_(n))
        {
          // The average based on the indicatrice does not guarantee the property pure_(i,j,k) = diph_l_(n) in the case diph_l_(n) == diph_v_(n).
          pure_(i,j,k) = diph_l_(n);
        }
      else
        {
          pure_(i,j,k) = indicatrice*diph_l_(n) + (1 - indicatrice)*diph_v_(n);
        }
    }
  IJK_Field_template<double,ArrOfDouble>::echange_espace_virtuel(IJK_Field_template<double,ArrOfDouble>::ghost());
}

template<>
void Cut_field_template<int,ArrOfInt>::remplir_tableau_pure_cellules_diphasiques_max(bool next_time)
{
  for (int n = 0; n < cut_cell_disc_->get_n_loc(); n++)
    {
      Int3 ijk = cut_cell_disc_->get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      if (diph_l_(n) == diph_v_(n))
        {
          // The average based on the indicatrice does not guarantee the property pure_(i,j,k) = diph_l_(n) in the case diph_l_(n) == diph_v_(n).
          pure_(i,j,k) = diph_l_(n);
        }
      else
        {
          pure_(i,j,k) = std::max(diph_l_(n), diph_v_(n));
        }
    }
  IJK_Field_template<int,ArrOfInt>::echange_espace_virtuel(IJK_Field_template<int,ArrOfInt>::ghost());
}


template<>
void Cut_field_template<double,ArrOfDouble>::set_field_data(const Nom& parser_expression_of_x_y_z_and_t, const IJK_Field_template<double,ArrOfDouble>& input_f, const double current_time)
{
  ArrOfDouble coord_i, coord_j, coord_k;
  build_local_coords(*this, coord_i, coord_j, coord_k);

  const int ni = IJK_Field_template<double,ArrOfDouble>::ni();
  const int nj = IJK_Field_template<double,ArrOfDouble>::nj();
  const int nk = IJK_Field_template<double,ArrOfDouble>::nk();

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
                  double vol_l = input_f(i, j, k);

                  double dx = cut_cell_disc_->get_domaine().get_constant_delta(DIRECTION_I);
                  double dy = cut_cell_disc_->get_domaine().get_constant_delta(DIRECTION_J);
                  double dz = cut_cell_disc_->get_domaine().get_constant_delta(DIRECTION_K);

                  assert(input_f.get_localisation() == localisation_);
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
  IJK_Field_template<double,ArrOfDouble>::echange_espace_virtuel(IJK_Field_template<double,ArrOfDouble>::ghost());
}

template<>
void Cut_field_template<double,ArrOfDouble>::dumplata_scalar(const char *filename, int step) const
{
  Cut_field_template<double,ArrOfDouble> *this_non_const = const_cast<Cut_field_template<double,ArrOfDouble>*>(this);
  this_non_const->echange_diph_vers_pure_cellules_finalement_pures();
  this_non_const->remplir_tableau_pure_cellules_diphasiques(true);

  IJK_Field_template<double,ArrOfDouble>::dumplata_scalar(filename, step);

  this_non_const->cut_cell_disc_->fill_buffer_with_variable(diph_l_);
  this_non_const->cut_cell_disc_->get_write_buffer().nommer(this->le_nom() + "_CUT_L");
  this_non_const->cut_cell_disc_->get_write_buffer().dumplata_scalar(filename, step);

  this_non_const->cut_cell_disc_->fill_buffer_with_variable(diph_v_);
  this_non_const->cut_cell_disc_->get_write_buffer().nommer(this->le_nom() + "_CUT_V");
  this_non_const->cut_cell_disc_->get_write_buffer().dumplata_scalar(filename, step);

}
template<>
void Cut_field_template<int,ArrOfInt>::dumplata_scalar(const char *filename, int step) const
{
  Cerr << "not implemented"<<finl;
  throw;

}

template<typename _TYPE_, typename _TYPE_ARRAY_>
void Cut_field_template<_TYPE_,_TYPE_ARRAY_>::set_to_uniform_value(_TYPE_ valeur)
{
  Cut_field_template<_TYPE_,_TYPE_ARRAY_>::data() = valeur;
  for (int n = 0; n < cut_cell_disc_->get_n_tot(); n++)
    {
      diph_l_(n) = valeur;
      diph_v_(n) = valeur;
    }
}


template<typename _TYPE_, typename _TYPE_ARRAY_>
void Cut_field_template<_TYPE_,_TYPE_ARRAY_>::associer_persistant(Cut_cell_FT_Disc& cut_cell_disc)
{
  cut_cell_disc_ = cut_cell_disc;
  diph_l_.associer_persistant(cut_cell_disc_, 1);
  diph_v_.associer_persistant(cut_cell_disc_, 1);
}

template<typename _TYPE_, typename _TYPE_ARRAY_>
void Cut_field_template<_TYPE_,_TYPE_ARRAY_>::associer_ephemere(Cut_cell_FT_Disc& cut_cell_disc)
{
  cut_cell_disc_ = cut_cell_disc;
  diph_l_.associer_ephemere(cut_cell_disc_, 1);
  diph_v_.associer_ephemere(cut_cell_disc_, 1);
}

template<typename _TYPE_, typename _TYPE_ARRAY_>
void Cut_field_template<_TYPE_,_TYPE_ARRAY_>::associer_paresseux(Cut_cell_FT_Disc& cut_cell_disc)
{
  cut_cell_disc_ = cut_cell_disc;
  diph_l_.associer_paresseux(cut_cell_disc_, 1);
  diph_v_.associer_paresseux(cut_cell_disc_, 1);
}

template<typename _TYPE_, typename _TYPE_ARRAY_>
void Cut_field_template<_TYPE_,_TYPE_ARRAY_>::copy_from(const Cut_field_template<_TYPE_,_TYPE_ARRAY_>& data)
{
  const int ni = IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ni();
  const int nj = IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::nj();
  const int nk = IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::nk();
  const int ghost = IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ghost();
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

template<typename _TYPE_, typename _TYPE_ARRAY_>
void Cut_field_template<_TYPE_,_TYPE_ARRAY_>::add_from(const Cut_field_template<_TYPE_,_TYPE_ARRAY_>& data, _TYPE_ constant)
{
  assert(constant == 1 || constant == -1);
  const int ni = data.ni();
  const int nj = data.nj();
  const int nk = data.nk();
  assert(((ni == IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ni())));
  assert(((nj == IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::nj())));
  assert(((nk == IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::nk())));

  const int ghost = IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ghost();
  if (ghost <= data.ghost())
    {
      for (int k = -ghost; k < nk+ghost; k++)
        {
          for (int j = -ghost; j < nj+ghost; j++)
            {
              for (int i = -ghost; i < ni+ghost; i++)
                {
                  pure_(i,j,k) += constant*data.pure_(i,j,k);
                }
            }
        }
    }
  else
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  pure_(i,j,k) += constant*data.pure_(i,j,k);
                }
            }
        }
      IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::echange_espace_virtuel(ghost);
    }

  for (int n = 0; n < cut_cell_disc_->get_n_tot(); n++)
    {
      diph_l_(n) += constant*data.diph_l_(n);
      diph_v_(n) += constant*data.diph_v_(n);
    }
}

template<typename _TYPE_, typename _TYPE_ARRAY_>
void Cut_field_template<_TYPE_,_TYPE_ARRAY_>::set_to_sum(const Cut_field_template<_TYPE_,_TYPE_ARRAY_>& data_1, const Cut_field_template<_TYPE_,_TYPE_ARRAY_>& data_2)
{
  const int ni = data_1.ni();
  const int nj = data_1.nj();
  const int nk = data_1.nk();
  assert(ni == data_2.ni());
  assert(nj == data_2.nj());
  assert(nk == data_2.nk());
  assert(((ni == IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ni())));
  assert(((nj == IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::nj())));
  assert(((nk == IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::nk())));

  const int ghost = IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ghost();
  if (ghost <= data_1.ghost() && ghost <= data_2.ghost())
    {
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
    }
  else
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  pure_(i,j,k) = data_1.pure_(i,j,k) + data_2.pure_(i,j,k);
                }
            }
        }
      IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::echange_espace_virtuel(ghost);
    }

  for (int n = 0; n < cut_cell_disc_->get_n_tot(); n++)
    {
      diph_l_(n) = data_1.diph_l_(n) + data_2.diph_l_(n);
      diph_v_(n) = data_1.diph_v_(n) + data_2.diph_v_(n);
    }
}

template<typename _TYPE_, typename _TYPE_ARRAY_>
_TYPE_& Cut_field_template<_TYPE_,_TYPE_ARRAY_>::from_signed_independent_index(int signed_independent_index)
{
  bool liquid = (signed_independent_index < 0);
  int independent_index = Cut_field_template<_TYPE_,_TYPE_ARRAY_>::domaine_ref_->get_independent_index_from_signed_independent_index(signed_independent_index);
  Int3 ijk = Cut_field_template<_TYPE_,_TYPE_ARRAY_>::domaine_ref_->get_ijk_from_independent_index(independent_index);
  int i = ijk[0];
  int j = ijk[1];
  int k = ijk[2];
  int n = cut_cell_disc_->get_n(i,j,k);
  assert((cut_cell_disc_->get_domaine() == Cut_field_template<_TYPE_,_TYPE_ARRAY_>::domaine_ref_));

  _TYPE_& value = (liquid && n >= 0) ? diph_l_(n) : ((n < 0) ? pure_(i,j,k) : diph_v_(n));
  return value;
}

template<typename _TYPE_, typename _TYPE_ARRAY_>
const _TYPE_& Cut_field_template<_TYPE_,_TYPE_ARRAY_>::from_signed_independent_index(int signed_independent_index) const
{
  bool liquid = (signed_independent_index < 0);
  int independent_index = Cut_field_template<_TYPE_,_TYPE_ARRAY_>::domaine_ref_->get_independent_index_from_signed_independent_index(signed_independent_index);
  Int3 ijk = Cut_field_template<_TYPE_,_TYPE_ARRAY_>::domaine_ref_->get_ijk_from_independent_index(independent_index);
  int i = ijk[0];
  int j = ijk[1];
  int k = ijk[2];
  int n = cut_cell_disc_->get_n(i,j,k);
  assert((cut_cell_disc_->get_domaine() == Cut_field_template<_TYPE_,_TYPE_ARRAY_>::domaine_ref_));

  const _TYPE_& value = (liquid && n >= 0) ? diph_l_(n) : ((n < 0) ? pure_(i,j,k) : diph_v_(n));
  return value;
}

template<typename _TYPE_, typename _TYPE_ARRAY_>
_TYPE_& Cut_field_template<_TYPE_,_TYPE_ARRAY_>::from_ijk_and_phase(int i, int j, int k, bool phase)
{
  assert((phase == 0) || (phase == 1));
  assert((cut_cell_disc_->get_domaine() == Cut_field_template<_TYPE_,_TYPE_ARRAY_>::domaine_ref_));

  int n = cut_cell_disc_->get_n(i, j, k);
  if (n < 0)
    {
      return pure_(i, j, k);
    }
  else
    {
      return (phase == 0) ? diph_v_(n) : diph_l_(n);
    }
}

template<typename _TYPE_, typename _TYPE_ARRAY_>
const _TYPE_& Cut_field_template<_TYPE_,_TYPE_ARRAY_>::from_ijk_and_phase(int i, int j, int k, bool phase) const
{
  assert((phase == 0) || (phase == 1));
  assert((cut_cell_disc_->get_domaine() == Cut_field_template<_TYPE_,_TYPE_ARRAY_>::domaine_ref_));

  int n = cut_cell_disc_->get_n(i, j, k);
  if (n < 0)
    {
      return pure_(i, j, k);
    }
  else
    {
      return (phase == 0) ? diph_v_(n) : diph_l_(n);
    }
}

template<typename _TYPE_, typename _TYPE_ARRAY_>
_TYPE_& Cut_field_template<_TYPE_,_TYPE_ARRAY_>::from_n_num_face_and_phase(int n, int num_face, bool phase)
{
  assert((phase == 0) || (phase == 1));
  assert((cut_cell_disc_->get_domaine() == Cut_field_template<_TYPE_,_TYPE_ARRAY_>::domaine_ref_));

  assert(n >= 0);
  return (phase == 0) ? diph_v_(n) : diph_l_(n);
}

template<typename _TYPE_, typename _TYPE_ARRAY_>
const _TYPE_& Cut_field_template<_TYPE_,_TYPE_ARRAY_>::from_n_num_face_and_phase(int n, int num_face, bool phase) const
{
  assert((phase == 0) || (phase == 1));
  assert((cut_cell_disc_->get_domaine() == Cut_field_template<_TYPE_,_TYPE_ARRAY_>::domaine_ref_));

  assert(n >= 0);
  return (phase == 0) ? diph_v_(n) : diph_l_(n);
}

template<>
CutCell_GlobalInfo Cut_field_template<double,ArrOfDouble>::compute_norm_cut_cell(bool next) const
{
  double norm_overall = 0.;
  double norm_overall_l = 0.;
  double norm_overall_v = 0.;
  double norm_pure = 0.;
  double norm_diph_l = 0.;
  double norm_diph_v = 0.;
  double norm_diph_small = 0.;
  double norm_diph_regular = 0.;
  double norm_diph_nascent = 0.;
  double norm_diph_dying = 0.;
  double count_overall = 0.;
  double count_overall_l = 0.;
  double count_overall_v = 0.;
  double count_pure = 0.;
  double count_diph_l = 0.;
  double count_diph_v = 0.;
  double count_diph_small = 0.;
  double count_diph_regular = 0.;
  double count_diph_nascent = 0.;
  double count_diph_dying = 0.;
  const IJK_Field_template<double,ArrOfDouble>& indic_old = cut_cell_disc_->get_interfaces().I();
  const IJK_Field_template<double,ArrOfDouble>& indic_next = cut_cell_disc_->get_interfaces().In();
  const int nx = IJK_Field_template<double,ArrOfDouble>::ni();
  const int ny = IJK_Field_template<double,ArrOfDouble>::nj();
  const int nz = IJK_Field_template<double,ArrOfDouble>::nk();
  for (int k=0; k < nz ; k++)
    {
      for (int j=0; j< ny; j++)
        {
          for (int i=0; i < nx; i++)
            {
              int n = cut_cell_disc_->get_n(i,j,k);
              if (n < 0)
                {
                  norm_overall += std::abs(pure_(i,j,k));
                  count_overall += 1;

                  norm_overall_l += (indic_old(i,j,k) == 0) ? 0. : std::abs(pure_(i,j,k));
                  count_overall_l += (indic_old(i,j,k) == 0) ? 0 : 1;

                  norm_overall_v += (indic_old(i,j,k) == 0) ? std::abs(pure_(i,j,k)) : 0.;
                  count_overall_v += (indic_old(i,j,k) == 0) ? 1 : 0;

                  norm_pure += std::abs(pure_(i,j,k));
                  count_pure += 1;
                }
              else if (IJK_Interfaces::est_pure(.5*(cut_cell_disc_->get_interfaces().I(i,j,k) + cut_cell_disc_->get_interfaces().In(i,j,k))))
                {
                  bool phase_invalide_l = (indic_old(i,j,k) == 0);
                  if (phase_invalide_l && diph_l_(n) != 0.)
                    {
                      Cerr << "compute_d_norm_cut_cell: There is a non-zero diph_l_(" << n << ") in an invalid cell (in the non-existant phase of a purely monophasic cell)." << finl;
                      Process::exit();
                    }

                  bool phase_invalide_v = (indic_old(i,j,k) == 1);
                  if (phase_invalide_v && diph_v_(n) != 0.)
                    {
                      Cerr << "compute_d_norm_cut_cell: There is a non-zero diph_v_(" << n << ") in an invalid cell (in the non-existant phase of a purely monophasic cell)." << finl;
                      Process::exit();
                    }

                  norm_overall += (indic_old(i,j,k) == 0) ? std::abs(diph_v_(n)) : std::abs(diph_l_(n));
                  count_overall += 1;

                  norm_overall_l += (indic_old(i,j,k) == 0) ? 0. : std::abs(diph_l_(n));
                  count_overall_l += (indic_old(i,j,k) == 0) ? 0 : 1;

                  norm_overall_v += (indic_old(i,j,k) == 0) ? std::abs(diph_v_(n)) : 0.;
                  count_overall_v += (indic_old(i,j,k) == 0) ? 1 : 0;

                  norm_pure += (indic_old(i,j,k) == 0) ? std::abs(diph_v_(n)) : std::abs(diph_l_(n));
                  count_pure += 1;
                }
              else
                {
                  double chi_T_l = std::abs(diph_l_(n));
                  double chi_T_v = std::abs(diph_v_(n));

                  norm_overall += chi_T_l;
                  norm_overall += chi_T_v;
                  count_overall += 1;

                  norm_overall_l += chi_T_l;
                  count_overall_l += 1;

                  norm_overall_v += chi_T_v;
                  count_overall_v += 1;

                  norm_diph_l += chi_T_l;
                  count_diph_l += 1;

                  norm_diph_v += chi_T_v;
                  count_diph_v += 1;

                  if (cut_cell_disc_->get_interfaces().devient_pure(i,j,k))
                    {
                      int phase_dying = IJK_Interfaces::convert_indicatrice_to_phase(1 - indic_next(i,j,k));
                      if (phase_dying == 1)
                        {
                          norm_diph_dying += chi_T_l;
                          count_diph_dying += 1;
                        }
                      else
                        {
                          norm_diph_dying += chi_T_v;
                          count_diph_dying += 1;
                        }
                    }
                  else if (cut_cell_disc_->get_interfaces().devient_diphasique(i,j,k))
                    {
                      int phase_nascent = IJK_Interfaces::convert_indicatrice_to_phase(1 - indic_old(i,j,k));
                      if (phase_nascent == 1)
                        {
                          norm_diph_nascent += chi_T_l;
                          count_diph_nascent += 1;
                        }
                      else
                        {
                          norm_diph_nascent += chi_T_v;
                          count_diph_nascent += 1;
                        }
                    }
                  else if (cut_cell_disc_->get_interfaces().next_below_small_threshold_for_phase(1, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      norm_diph_small += chi_T_l;
                      count_diph_small += 1;

                      norm_diph_regular += chi_T_v;
                      count_diph_regular += 1;
                    }
                  else if (cut_cell_disc_->get_interfaces().next_below_small_threshold_for_phase(0, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      norm_diph_small += chi_T_v;
                      count_diph_small += 1;

                      norm_diph_regular += chi_T_l;
                      count_diph_regular += 1;
                    }
                  else
                    {
                      norm_diph_regular += chi_T_l;
                      norm_diph_regular += chi_T_v;
                      count_diph_regular += 1;
                    }

                }
            }
        }
    }
  count_overall = Process::mp_sum(count_overall);
  count_overall_l = Process::mp_sum(count_overall_l);
  count_overall_v = Process::mp_sum(count_overall_v);
  count_pure = Process::mp_sum(count_pure);
  count_diph_l = Process::mp_sum(count_diph_l);
  count_diph_v = Process::mp_sum(count_diph_v);
  count_diph_small = Process::mp_sum(count_diph_small);
  count_diph_regular = Process::mp_sum(count_diph_regular);
  count_diph_nascent = Process::mp_sum(count_diph_nascent);
  count_diph_dying = Process::mp_sum(count_diph_dying);
  assert(count_overall == domaine_ref_->get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_I)
         *domaine_ref_->get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_J)
         *domaine_ref_->get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_K));
  norm_overall      = Process::mp_sum(norm_overall);
  norm_overall_l    = Process::mp_sum(norm_overall_l);
  norm_overall_v    = Process::mp_sum(norm_overall_v);
  norm_pure         = Process::mp_sum(norm_pure);
  norm_diph_l       = Process::mp_sum(norm_diph_l);
  norm_diph_v       = Process::mp_sum(norm_diph_v);
  norm_diph_small   = Process::mp_sum(norm_diph_small);
  norm_diph_regular = Process::mp_sum(norm_diph_regular);
  norm_diph_nascent = Process::mp_sum(norm_diph_nascent);
  norm_diph_dying   = Process::mp_sum(norm_diph_dying);

  if (count_overall      != 0)    norm_overall      /= count_overall;
  if (count_overall_l    != 0)    norm_overall_l    /= count_overall_l;
  if (count_overall_v    != 0)    norm_overall_v    /= count_overall_v;
  if (count_pure         != 0)    norm_pure         /= count_pure;
  if (count_diph_l       != 0)    norm_diph_l       /= count_diph_l;
  if (count_diph_v       != 0)    norm_diph_v       /= count_diph_v;
  if (count_diph_small   != 0)    norm_diph_small   /= count_diph_small;
  if (count_diph_regular != 0)    norm_diph_regular /= count_diph_regular;
  if (count_diph_nascent != 0)    norm_diph_nascent /= count_diph_nascent;
  if (count_diph_dying   != 0)    norm_diph_dying   /= count_diph_dying;

  return {norm_overall, norm_overall_l, norm_overall_v, norm_pure, norm_diph_l, norm_diph_v, norm_diph_small, norm_diph_regular, norm_diph_nascent, norm_diph_dying};
}

template<>
CutCell_GlobalInfo Cut_field_template<double,ArrOfDouble>::compute_d_global_energy_cut_cell(bool next) const
{
  double global_energy_overall = 0.;
  double global_energy_overall_l = 0.;
  double global_energy_overall_v = 0.;
  double global_energy_pure = 0.;
  double global_energy_diph_l = 0.;
  double global_energy_diph_v = 0.;
  double global_energy_diph_small = 0.;
  double global_energy_diph_regular = 0.;
  double global_energy_diph_nascent = 0.;
  double global_energy_diph_dying = 0.;
  double count_overall = 0.;
  double count_overall_l = 0.;
  double count_overall_v = 0.;
  double count_pure = 0.;
  double count_diph_l = 0.;
  double count_diph_v = 0.;
  double count_diph_small = 0.;
  double count_diph_regular = 0.;
  double count_diph_nascent = 0.;
  double count_diph_dying = 0.;
  const IJK_Field_template<double,ArrOfDouble>& indic_old = cut_cell_disc_->get_interfaces().I();
  const IJK_Field_template<double,ArrOfDouble>& indic_next = cut_cell_disc_->get_interfaces().In();
  const int nx = IJK_Field_template<double,ArrOfDouble>::ni();
  const int ny = IJK_Field_template<double,ArrOfDouble>::nj();
  const int nz = IJK_Field_template<double,ArrOfDouble>::nk();
  const Domaine_IJK& geom = indic_next.get_domaine();
  assert(geom.get_constant_delta(DIRECTION_K) >0); // To be sure we're on a regular mesh
  for (int k=0; k < nz ; k++)
    {
      for (int j=0; j< ny; j++)
        {
          for (int i=0; i < nx; i++)
            {
              int n = cut_cell_disc_->get_n(i,j,k);
              if (n < 0)
                {
                  global_energy_overall += pure_(i,j,k);
                  count_overall += 1;

                  global_energy_overall_l += (indic_old(i,j,k) == 0) ? 0. : pure_(i,j,k);
                  count_overall_l += (indic_old(i,j,k) == 0) ? 0 : 1;

                  global_energy_overall_v += (indic_old(i,j,k) == 0) ? pure_(i,j,k) : 0.;
                  count_overall_v += (indic_old(i,j,k) == 0) ? 1 : 0;

                  global_energy_pure += pure_(i,j,k);
                  count_pure += 1;
                }
              else if (IJK_Interfaces::est_pure(.5*(cut_cell_disc_->get_interfaces().I(i,j,k) + cut_cell_disc_->get_interfaces().In(i,j,k))))
                {
                  bool phase_invalide_l = (indic_old(i,j,k) == 0);
                  if (phase_invalide_l && diph_l_(n) != 0.)
                    {
                      Cerr << "compute_d_global_energy_cut_cell: There is a non-zero diph_l_(" << n << ") in an invalid cell (in the non-existant phase of a purely monophasic cell)." << finl;
                      Process::exit();
                    }

                  bool phase_invalide_v = (indic_old(i,j,k) == 1);
                  if (phase_invalide_v && diph_v_(n) != 0.)
                    {
                      Cerr << "compute_d_global_energy_cut_cell: There is a non-zero diph_v_(" << n << ") in an invalid cell (in the non-existant phase of a purely monophasic cell)." << finl;
                      Process::exit();
                    }

                  global_energy_overall += (indic_old(i,j,k) == 0) ? diph_v_(n) : diph_l_(n);
                  count_overall += 1;

                  global_energy_overall_l += (indic_old(i,j,k) == 0) ? 0. : diph_l_(n);
                  count_overall_l += (indic_old(i,j,k) == 0) ? 0 : 1;

                  global_energy_overall_v += (indic_old(i,j,k) == 0) ? diph_v_(n) : 0.;
                  count_overall_v += (indic_old(i,j,k) == 0) ? 1 : 0;

                  global_energy_pure += (indic_old(i,j,k) == 0) ? diph_v_(n) : diph_l_(n);
                  count_pure += 1;
                }
              else
                {
                  double chi_T_l = diph_l_(n);
                  double chi_T_v = diph_v_(n);

                  global_energy_overall += chi_T_l;
                  global_energy_overall += chi_T_v;
                  count_overall += 1;

                  global_energy_overall_l += chi_T_l;
                  count_overall_l += 1;

                  global_energy_overall_v += chi_T_v;
                  count_overall_v += 1;

                  global_energy_diph_l += chi_T_l;
                  count_diph_l += 1;

                  global_energy_diph_v += chi_T_v;
                  count_diph_v += 1;

                  if (cut_cell_disc_->get_interfaces().devient_pure(i,j,k))
                    {
                      int phase_dying = IJK_Interfaces::convert_indicatrice_to_phase(1 - indic_next(i,j,k));
                      if (phase_dying == 1)
                        {
                          global_energy_diph_dying += chi_T_l;
                          count_diph_dying += 1;
                        }
                      else
                        {
                          global_energy_diph_dying += chi_T_v;
                          count_diph_dying += 1;
                        }
                    }
                  else if (cut_cell_disc_->get_interfaces().devient_diphasique(i,j,k))
                    {
                      int phase_nascent = IJK_Interfaces::convert_indicatrice_to_phase(1 - indic_old(i,j,k));
                      if (phase_nascent == 1)
                        {
                          global_energy_diph_nascent += chi_T_l;
                          count_diph_nascent += 1;
                        }
                      else
                        {
                          global_energy_diph_nascent += chi_T_v;
                          count_diph_nascent += 1;
                        }
                    }
                  else if (cut_cell_disc_->get_interfaces().next_below_small_threshold_for_phase(1, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      global_energy_diph_small += chi_T_l;
                      count_diph_small += 1;

                      global_energy_diph_regular += chi_T_v;
                      count_diph_regular += 1;
                    }
                  else if (cut_cell_disc_->get_interfaces().next_below_small_threshold_for_phase(0, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      global_energy_diph_small += chi_T_v;
                      count_diph_small += 1;

                      global_energy_diph_regular += chi_T_l;
                      count_diph_regular += 1;
                    }
                  else
                    {
                      global_energy_diph_regular += chi_T_l;
                      global_energy_diph_regular += chi_T_v;
                      count_diph_regular += 1;
                    }

                }
            }
        }
    }
  count_overall = Process::mp_sum(count_overall);
  count_overall_l = Process::mp_sum(count_overall_l);
  count_overall_v = Process::mp_sum(count_overall_v);
  count_pure = Process::mp_sum(count_pure);
  count_diph_l = Process::mp_sum(count_diph_l);
  count_diph_v = Process::mp_sum(count_diph_v);
  count_diph_small = Process::mp_sum(count_diph_small);
  count_diph_regular = Process::mp_sum(count_diph_regular);
  count_diph_nascent = Process::mp_sum(count_diph_nascent);
  count_diph_dying = Process::mp_sum(count_diph_dying);
  assert(count_overall == domaine_ref_->get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_I)
         *domaine_ref_->get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_J)
         *domaine_ref_->get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_K));
  const double vol_cell = geom.get_constant_delta(DIRECTION_I)*geom.get_constant_delta(DIRECTION_J)*geom.get_constant_delta(DIRECTION_K);
  global_energy_overall      = vol_cell * mp_sum(global_energy_overall);
  global_energy_overall_l    = vol_cell * mp_sum(global_energy_overall_l);
  global_energy_overall_v    = vol_cell * mp_sum(global_energy_overall_v);
  global_energy_pure         = vol_cell * mp_sum(global_energy_pure);
  global_energy_diph_l       = vol_cell * mp_sum(global_energy_diph_l);
  global_energy_diph_v       = vol_cell * mp_sum(global_energy_diph_v);
  global_energy_diph_small   = vol_cell * mp_sum(global_energy_diph_small);
  global_energy_diph_regular = vol_cell * mp_sum(global_energy_diph_regular);
  global_energy_diph_nascent = vol_cell * mp_sum(global_energy_diph_nascent);
  global_energy_diph_dying   = vol_cell * mp_sum(global_energy_diph_dying);
  return {global_energy_overall, global_energy_overall_l, global_energy_overall_v, global_energy_pure, global_energy_diph_l, global_energy_diph_v, global_energy_diph_small, global_energy_diph_regular, global_energy_diph_nascent, global_energy_diph_dying};
}

template<>
CutCell_GlobalInfo Cut_field_template<double,ArrOfDouble>::compute_global_energy_cut_cell(bool next, double constant_l, double constant_v) const
{
  double global_energy_overall = 0.;
  double global_energy_overall_l = 0.;
  double global_energy_overall_v = 0.;
  double global_energy_pure = 0.;
  double global_energy_diph_l = 0.;
  double global_energy_diph_v = 0.;
  double global_energy_diph_small = 0.;
  double global_energy_diph_regular = 0.;
  double global_energy_diph_nascent = 0.;
  double global_energy_diph_dying = 0.;
  double count_overall = 0.;
  double count_overall_l = 0.;
  double count_overall_v = 0.;
  double count_pure = 0.;
  double count_diph_l = 0.;
  double count_diph_v = 0.;
  double count_diph_small = 0.;
  double count_diph_regular = 0.;
  double count_diph_nascent = 0.;
  double count_diph_dying = 0.;
  const IJK_Field_template<double,ArrOfDouble>& indic_old = cut_cell_disc_->get_interfaces().I();
  const IJK_Field_template<double,ArrOfDouble>& indic_next = cut_cell_disc_->get_interfaces().In();
  const int nx = IJK_Field_template<double,ArrOfDouble>::ni();
  const int ny = IJK_Field_template<double,ArrOfDouble>::nj();
  const int nz = IJK_Field_template<double,ArrOfDouble>::nk();
  const Domaine_IJK& geom = indic_next.get_domaine();
  assert(geom.get_constant_delta(DIRECTION_K) >0); // To be sure we're on a regular mesh
  for (int k=0; k < nz ; k++)
    {
      for (int j=0; j< ny; j++)
        {
          for (int i=0; i < nx; i++)
            {
              double chi_l = next ? cut_cell_disc_->get_interfaces().In(i,j,k) : cut_cell_disc_->get_interfaces().I(i,j,k);
              double chi_v = next ? 1-cut_cell_disc_->get_interfaces().In(i,j,k) : 1-cut_cell_disc_->get_interfaces().I(i,j,k);
              double chi_nonzero_l = next ? cut_cell_disc_->get_interfaces().In_nonzero(1,i,j,k) : cut_cell_disc_->get_interfaces().I_nonzero(1,i,j,k);
              double chi_nonzero_v = next ? cut_cell_disc_->get_interfaces().In_nonzero(0,i,j,k) : cut_cell_disc_->get_interfaces().I_nonzero(0,i,j,k);
              int n = cut_cell_disc_->get_n(i,j,k);
              if (n < 0)
                {
                  global_energy_overall += (chi_l * constant_l + chi_v * constant_v) * pure_(i,j,k);
                  count_overall += 1;

                  global_energy_overall_l += chi_l * constant_l * pure_(i,j,k);
                  count_overall_l += chi_l;

                  global_energy_overall_v += chi_v * constant_v * pure_(i,j,k);
                  count_overall_v += chi_v;

                  global_energy_pure += (chi_l * constant_l + chi_v * constant_v) * pure_(i,j,k);
                  count_pure += 1;
                }
              else if (IJK_Interfaces::est_pure(.5*(cut_cell_disc_->get_interfaces().I(i,j,k) + cut_cell_disc_->get_interfaces().In(i,j,k))))
                {
                  bool phase_invalide_l = (indic_old(i,j,k) == 0);
                  if (phase_invalide_l && diph_l_(n) != 0.)
                    {
                      Cerr << "compute_global_energy_cut_cell: There is a non-zero diph_l_(" << n << ") in an invalid cell (in the non-existant phase of a purely monophasic cell)." << finl;
                      Process::exit();
                    }

                  bool phase_invalide_v = (indic_old(i,j,k) == 1);
                  if (phase_invalide_v && diph_v_(n) != 0.)
                    {
                      Cerr << "compute_global_energy_cut_cell: There is a non-zero diph_v_(" << n << ") in an invalid cell (in the non-existant phase of a purely monophasic cell)." << finl;
                      Process::exit();
                    }

                  global_energy_overall += (chi_l * constant_l + chi_v * constant_v) * ((indic_old(i,j,k) == 0) ? diph_v_(n) : diph_l_(n));
                  count_overall += 1;

                  global_energy_overall_l += chi_l * constant_l * diph_l_(n);
                  count_overall_l += chi_l;

                  global_energy_overall_v += chi_v * constant_v * diph_v_(n);
                  count_overall_v += chi_v;

                  global_energy_pure += (chi_l * constant_l + chi_v * constant_v) * ((indic_old(i,j,k) == 0) ? diph_v_(n) : diph_l_(n));
                  count_pure += 1;
                }
              else
                {
                  double chi_T_l = chi_nonzero_l * diph_l_(n);
                  double chi_T_v = chi_nonzero_v * diph_v_(n);

                  global_energy_overall += constant_l * chi_T_l;
                  global_energy_overall += constant_v * chi_T_v;
                  count_overall += 1;

                  global_energy_overall_l += constant_l * chi_T_l;
                  count_overall_l += 1;

                  global_energy_overall_v += constant_v * chi_T_v;
                  count_overall_v += 1;

                  global_energy_diph_l += constant_l * chi_T_l;
                  count_diph_l += chi_nonzero_l;

                  global_energy_diph_v += constant_v * chi_T_v;
                  count_diph_v += chi_nonzero_v;

                  if (cut_cell_disc_->get_interfaces().devient_pure(i,j,k))
                    {
                      int phase_dying = IJK_Interfaces::convert_indicatrice_to_phase(1 - indic_next(i,j,k));
                      if (phase_dying == 1)
                        {
                          global_energy_diph_dying += constant_l * chi_T_l;
                          count_diph_dying += chi_nonzero_l;
                        }
                      else
                        {
                          global_energy_diph_dying += constant_v * chi_T_v;
                          count_diph_dying += chi_nonzero_v;
                        }
                    }
                  else if (cut_cell_disc_->get_interfaces().devient_diphasique(i,j,k))
                    {
                      int phase_nascent = IJK_Interfaces::convert_indicatrice_to_phase(1 - indic_old(i,j,k));
                      if (phase_nascent == 1)
                        {
                          global_energy_diph_nascent += constant_l * chi_T_l;
                          count_diph_nascent += chi_nonzero_l;
                        }
                      else
                        {
                          global_energy_diph_nascent += constant_v * chi_T_v;
                          count_diph_nascent += chi_nonzero_v;
                        }
                    }
                  else if (cut_cell_disc_->get_interfaces().next_below_small_threshold_for_phase(1, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      global_energy_diph_small += constant_l * chi_T_l;
                      count_diph_small += chi_nonzero_l;

                      global_energy_diph_regular += constant_v * chi_T_v;
                      count_diph_regular += chi_nonzero_v;
                    }
                  else if (cut_cell_disc_->get_interfaces().next_below_small_threshold_for_phase(0, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      global_energy_diph_small += constant_v * chi_T_v;
                      count_diph_small += chi_nonzero_v;

                      global_energy_diph_regular += constant_l * chi_T_l;
                      count_diph_regular += chi_nonzero_l;
                    }
                  else
                    {
                      global_energy_diph_regular += constant_l * chi_T_l;
                      global_energy_diph_regular += constant_v * chi_T_v;
                      count_diph_regular += 1;
                    }

                }
            }
        }
    }
  count_overall = Process::mp_sum(count_overall);
  count_overall_l = Process::mp_sum(count_overall_l);
  count_overall_v = Process::mp_sum(count_overall_v);
  count_pure = Process::mp_sum(count_pure);
  count_diph_l = Process::mp_sum(count_diph_l);
  count_diph_v = Process::mp_sum(count_diph_v);
  count_diph_small = Process::mp_sum(count_diph_small);
  count_diph_regular = Process::mp_sum(count_diph_regular);
  count_diph_nascent = Process::mp_sum(count_diph_nascent);
  count_diph_dying = Process::mp_sum(count_diph_dying);
  assert(count_overall == domaine_ref_->get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_I)
         *domaine_ref_->get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_J)
         *domaine_ref_->get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_K));
  const double vol_cell = geom.get_constant_delta(DIRECTION_I)*geom.get_constant_delta(DIRECTION_J)*geom.get_constant_delta(DIRECTION_K);
  global_energy_overall      = vol_cell * mp_sum(global_energy_overall);
  global_energy_overall_l    = vol_cell * mp_sum(global_energy_overall_l);
  global_energy_overall_v    = vol_cell * mp_sum(global_energy_overall_v);
  global_energy_pure         = vol_cell * mp_sum(global_energy_pure);
  global_energy_diph_l       = vol_cell * mp_sum(global_energy_diph_l);
  global_energy_diph_v       = vol_cell * mp_sum(global_energy_diph_v);
  global_energy_diph_small   = vol_cell * mp_sum(global_energy_diph_small);
  global_energy_diph_regular = vol_cell * mp_sum(global_energy_diph_regular);
  global_energy_diph_nascent = vol_cell * mp_sum(global_energy_diph_nascent);
  global_energy_diph_dying   = vol_cell * mp_sum(global_energy_diph_dying);
  return {global_energy_overall, global_energy_overall_l, global_energy_overall_v, global_energy_pure, global_energy_diph_l, global_energy_diph_v, global_energy_diph_small, global_energy_diph_regular, global_energy_diph_nascent, global_energy_diph_dying};
}

template<>
CutCell_GlobalInfo Cut_field_template<double,ArrOfDouble>::compute_min_cut_cell(bool next) const
{
  double Tmin_overall = 1.e20;
  double Tmin_overall_l = 1.e20;
  double Tmin_overall_v = 1.e20;
  double Tmin_pure = 1.e20;
  double Tmin_diph_l = 1.e20;
  double Tmin_diph_v = 1.e20;
  double Tmin_diph_small = 1.e20;
  double Tmin_diph_regular = 1.e20;
  double Tmin_diph_nascent = 1.e20;
  double Tmin_diph_dying = 1.e20;
  const IJK_Field_template<double,ArrOfDouble>& indic_old = cut_cell_disc_->get_interfaces().I();
  const IJK_Field_template<double,ArrOfDouble>& indic_next = cut_cell_disc_->get_interfaces().In();
  const IJK_Field_template<double,ArrOfDouble>& indic = next ? indic_next : indic_old;
  const int nx = IJK_Field_template<double,ArrOfDouble>::ni();
  const int ny = IJK_Field_template<double,ArrOfDouble>::nj();
  const int nz = IJK_Field_template<double,ArrOfDouble>::nk();
  assert(indic.get_domaine().get_constant_delta(DIRECTION_K) >0); // To be sure we're on a regular mesh
  for (int k=0; k < nz ; k++)
    {
      for (int j=0; j< ny; j++)
        {
          for (int i=0; i < nx; i++)
            {
              double chi_l = indic(i,j,k);
              int n = cut_cell_disc_->get_n(i,j,k);
              if (n < 0)
                {
                  Tmin_overall = std::min(Tmin_overall, pure_(i,j,k));
                  Tmin_overall_l = (chi_l == 0) ? Tmin_overall_l : std::min(Tmin_overall_l, pure_(i,j,k));
                  Tmin_overall_v = (chi_l == 1) ? Tmin_overall_v : std::min(Tmin_overall_v, pure_(i,j,k));
                  Tmin_pure = std::min(Tmin_pure, pure_(i,j,k));
                }
              else if (IJK_Interfaces::est_pure(.5*(cut_cell_disc_->get_interfaces().I(i,j,k) + cut_cell_disc_->get_interfaces().In(i,j,k))))
                {
                  Tmin_overall = std::min(Tmin_overall, (chi_l == 0) ? diph_v_(n) : diph_l_(n));
                  Tmin_overall_l = (chi_l == 0) ? Tmin_overall_l : std::min(Tmin_overall_l, diph_l_(n));
                  Tmin_overall_v = (chi_l == 1) ? Tmin_overall_v : std::min(Tmin_overall_v, diph_v_(n));
                  Tmin_pure = std::min(Tmin_pure, (chi_l == 0) ? diph_v_(n) : diph_l_(n));
                }
              else
                {
                  // Excluding the value of the phase in dying cells
                  bool exclude_l = (next && (cut_cell_disc_->get_interfaces().phase_mourrante(1, i,j,k))) || ((!next) && (cut_cell_disc_->get_interfaces().phase_naissante(1, i,j,k)));
                  bool exclude_v = (next && (cut_cell_disc_->get_interfaces().phase_mourrante(0, i,j,k))) || ((!next) && (cut_cell_disc_->get_interfaces().phase_naissante(0, i,j,k)));

                  Tmin_overall = exclude_l ? Tmin_overall : std::min(Tmin_overall, diph_l_(n));
                  Tmin_overall = exclude_v ? Tmin_overall : std::min(Tmin_overall, diph_v_(n));

                  Tmin_overall_l = exclude_l ? Tmin_overall_l : std::min(Tmin_overall_l, diph_l_(n));
                  Tmin_overall_v = exclude_v ? Tmin_overall_v : std::min(Tmin_overall_v, diph_v_(n));

                  Tmin_diph_l = exclude_l ? Tmin_diph_l : std::min(Tmin_diph_l, diph_l_(n));
                  Tmin_diph_v = exclude_v ? Tmin_diph_v : std::min(Tmin_diph_v, diph_v_(n));

                  if (cut_cell_disc_->get_interfaces().next_below_small_threshold_for_phase(1, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      Tmin_diph_small = exclude_l ? Tmin_diph_small : std::min(Tmin_diph_small, diph_l_(n));
                      Tmin_diph_regular = exclude_v ? Tmin_diph_regular : std::min(Tmin_diph_regular, diph_v_(n));
                    }
                  else if (cut_cell_disc_->get_interfaces().next_below_small_threshold_for_phase(0, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      Tmin_diph_small = exclude_v ? Tmin_diph_small : std::min(Tmin_diph_small, diph_v_(n));
                      Tmin_diph_regular = exclude_l ? Tmin_diph_regular : std::min(Tmin_diph_regular, diph_l_(n));
                    }
                  else
                    {
                      Tmin_diph_regular = exclude_l ? Tmin_diph_regular : std::min(Tmin_diph_regular, diph_l_(n));
                      Tmin_diph_regular = exclude_v ? Tmin_diph_regular : std::min(Tmin_diph_regular, diph_v_(n));
                    }

                  if (cut_cell_disc_->get_interfaces().devient_pure(i,j,k))
                    {
                      int phase_dying = IJK_Interfaces::convert_indicatrice_to_phase(1 - indic_next(i,j,k));
                      if (phase_dying == 1)
                        {
                          Tmin_diph_dying = std::min(Tmin_diph_dying, diph_l_(n));
                        }
                      else
                        {
                          Tmin_diph_dying = std::min(Tmin_diph_dying, diph_v_(n));
                        }
                    }

                  if (cut_cell_disc_->get_interfaces().devient_diphasique(i,j,k))
                    {
                      int phase_nascent = IJK_Interfaces::convert_indicatrice_to_phase(1 - indic_old(i,j,k));
                      if (phase_nascent == 1)
                        {
                          Tmin_diph_nascent = std::min(Tmin_diph_nascent, diph_l_(n));
                        }
                      else
                        {
                          Tmin_diph_nascent = std::min(Tmin_diph_nascent, diph_v_(n));
                        }
                    }
                }
            }
        }
    }
  Tmin_overall      = Process::mp_min(Tmin_overall);
  Tmin_overall_l      = Process::mp_min(Tmin_overall_l);
  Tmin_overall_v      = Process::mp_min(Tmin_overall_v);
  Tmin_pure         = Process::mp_min(Tmin_pure);
  Tmin_diph_l       = Process::mp_min(Tmin_diph_l);
  Tmin_diph_v       = Process::mp_min(Tmin_diph_v);
  Tmin_diph_small   = Process::mp_min(Tmin_diph_small);
  Tmin_diph_regular   = Process::mp_min(Tmin_diph_regular);
  Tmin_diph_nascent = Process::mp_min(Tmin_diph_nascent);
  Tmin_diph_dying   = Process::mp_min(Tmin_diph_dying);
  return {Tmin_overall, Tmin_overall_l, Tmin_overall_v, Tmin_pure, Tmin_diph_l, Tmin_diph_v, Tmin_diph_small, Tmin_diph_regular, Tmin_diph_nascent, Tmin_diph_dying};
}

template<>
CutCell_GlobalInfo Cut_field_template<double,ArrOfDouble>::compute_max_cut_cell(bool next) const
{
  double Tmax_overall = -1.e20;
  double Tmax_overall_l = -1.e20;
  double Tmax_overall_v = -1.e20;
  double Tmax_pure = -1.e20;
  double Tmax_diph_l = -1.e20;
  double Tmax_diph_v = -1.e20;
  double Tmax_diph_small = -1.e20;
  double Tmax_diph_regular = -1.e20;
  double Tmax_diph_nascent = -1.e20;
  double Tmax_diph_dying = -1.e20;
  const IJK_Field_template<double,ArrOfDouble>& indic_old = cut_cell_disc_->get_interfaces().I();
  const IJK_Field_template<double,ArrOfDouble>& indic_next = cut_cell_disc_->get_interfaces().In();
  const IJK_Field_template<double,ArrOfDouble>& indic = next ? indic_next : indic_old;
  const int nx = IJK_Field_template<double,ArrOfDouble>::ni();
  const int ny = IJK_Field_template<double,ArrOfDouble>::nj();
  const int nz = IJK_Field_template<double,ArrOfDouble>::nk();
  // To be sure we're on a regular mesh
  assert(indic.get_domaine().get_constant_delta(DIRECTION_K) >0);
  for (int k=0; k < nz ; k++)
    {
      for (int j=0; j< ny; j++)
        {
          for (int i=0; i < nx; i++)
            {
              double chi_l = indic(i,j,k);
              int n = cut_cell_disc_->get_n(i,j,k);
              if (n < 0)
                {
                  Tmax_overall = std::max(Tmax_overall, pure_(i,j,k));
                  Tmax_overall_l = (chi_l == 0) ? Tmax_overall_l : std::max(Tmax_overall_l, pure_(i,j,k));
                  Tmax_overall_v = (chi_l == 1) ? Tmax_overall_v : std::max(Tmax_overall_v, pure_(i,j,k));
                  Tmax_pure = std::max(Tmax_pure, pure_(i,j,k));
                }
              else if (IJK_Interfaces::est_pure(.5*(cut_cell_disc_->get_interfaces().I(i,j,k) + cut_cell_disc_->get_interfaces().In(i,j,k))))
                {
                  Tmax_overall = std::max(Tmax_overall, (chi_l == 0) ? diph_v_(n) : diph_l_(n));
                  Tmax_overall_l = (chi_l == 0) ? Tmax_overall_l : std::max(Tmax_overall_l, diph_l_(n));
                  Tmax_overall_v = (chi_l == 1) ? Tmax_overall_v : std::max(Tmax_overall_v, diph_v_(n));
                  Tmax_pure = std::max(Tmax_pure, (chi_l == 0) ? diph_v_(n) : diph_l_(n));
                }
              else
                {
                  // Excluding the value of the phase in dying cells
                  bool exclude_l = (next && (cut_cell_disc_->get_interfaces().phase_mourrante(1, i,j,k))) || ((!next) && (cut_cell_disc_->get_interfaces().phase_naissante(1, i,j,k)));
                  bool exclude_v = (next && (cut_cell_disc_->get_interfaces().phase_mourrante(0, i,j,k))) || ((!next) && (cut_cell_disc_->get_interfaces().phase_naissante(0, i,j,k)));

                  Tmax_overall = exclude_l ? Tmax_overall : std::max(Tmax_overall, diph_l_(n));
                  Tmax_overall = exclude_v ? Tmax_overall : std::max(Tmax_overall, diph_v_(n));

                  Tmax_overall_l = exclude_l ? Tmax_overall_l : std::max(Tmax_overall_l, diph_l_(n));
                  Tmax_overall_v = exclude_v ? Tmax_overall_v : std::max(Tmax_overall_v, diph_v_(n));

                  Tmax_diph_l = exclude_l ? Tmax_diph_l : std::max(Tmax_diph_l, diph_l_(n));
                  Tmax_diph_v = exclude_v ? Tmax_diph_v : std::max(Tmax_diph_v, diph_v_(n));

                  if (cut_cell_disc_->get_interfaces().next_below_small_threshold_for_phase(1, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      Tmax_diph_small = exclude_l ? Tmax_diph_small : std::max(Tmax_diph_small, diph_l_(n));
                      Tmax_diph_regular = exclude_v ? Tmax_diph_regular : std::max(Tmax_diph_regular, diph_v_(n));
                    }
                  else if (cut_cell_disc_->get_interfaces().next_below_small_threshold_for_phase(0, indic_old(i,j,k), indic_next(i,j,k)))
                    {
                      Tmax_diph_small = exclude_v ? Tmax_diph_small : std::max(Tmax_diph_small, diph_v_(n));
                      Tmax_diph_regular = exclude_l ? Tmax_diph_regular : std::max(Tmax_diph_regular, diph_l_(n));
                    }
                  else
                    {
                      Tmax_diph_regular = exclude_l ? Tmax_diph_regular : std::max(Tmax_diph_regular, diph_l_(n));
                      Tmax_diph_regular = exclude_v ? Tmax_diph_regular : std::max(Tmax_diph_regular, diph_v_(n));
                    }

                  if (cut_cell_disc_->get_interfaces().devient_pure(i,j,k))
                    {
                      int phase_dying = IJK_Interfaces::convert_indicatrice_to_phase(1 - indic_next(i,j,k));
                      if (phase_dying == 1)
                        {
                          Tmax_diph_dying = std::max(Tmax_diph_dying, diph_l_(n));
                        }
                      else
                        {
                          Tmax_diph_dying = std::max(Tmax_diph_dying, diph_v_(n));
                        }
                    }

                  if (cut_cell_disc_->get_interfaces().devient_diphasique(i,j,k))
                    {
                      int phase_nascent = IJK_Interfaces::convert_indicatrice_to_phase(1 - indic_old(i,j,k));
                      if (phase_nascent == 1)
                        {
                          Tmax_diph_nascent = std::max(Tmax_diph_nascent, diph_l_(n));
                        }
                      else
                        {
                          Tmax_diph_nascent = std::max(Tmax_diph_nascent, diph_v_(n));
                        }
                    }
                }
            }
        }
    }
  Tmax_overall      = Process::mp_max(Tmax_overall);
  Tmax_overall_l      = Process::mp_max(Tmax_overall_l);
  Tmax_overall_v      = Process::mp_max(Tmax_overall_v);
  Tmax_pure         = Process::mp_max(Tmax_pure);
  Tmax_diph_l       = Process::mp_max(Tmax_diph_l);
  Tmax_diph_v       = Process::mp_max(Tmax_diph_v);
  Tmax_diph_small   = Process::mp_max(Tmax_diph_small);
  Tmax_diph_regular   = Process::mp_max(Tmax_diph_regular);
  Tmax_diph_nascent = Process::mp_max(Tmax_diph_nascent);
  Tmax_diph_dying   = Process::mp_max(Tmax_diph_dying);
  return {Tmax_overall, Tmax_overall_l, Tmax_overall_v, Tmax_pure, Tmax_diph_l, Tmax_diph_v, Tmax_diph_small, Tmax_diph_regular, Tmax_diph_nascent, Tmax_diph_dying};
}

template<typename _TYPE_, typename _TYPE_ARRAY_>
Nom Cut_field_template<_TYPE_,_TYPE_ARRAY_>::get_value_location(_TYPE_ T) const
{
  const int nx = IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::ni();
  const int ny = IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::nj();
  const int nz = IJK_Field_template<_TYPE_,_TYPE_ARRAY_>::nk();
  const Cut_cell_FT_Disc& cut_cell_disc = get_cut_cell_disc();
  for (int k=0; k < nz ; k++)
    {
      for (int j=0; j< ny; j++)
        {
          for (int i=0; i < nx; i++)
            {
              int n = cut_cell_disc.get_n(i,j,k);
              if (n < 0)
                {
                  if (T == pure_(i,j,k))
                    {
                      return Nom("ijk_pure_") + Nom(i) + Nom("_") + Nom(j) + Nom("_") + Nom(k);
                    }
                }
              else
                {
                  if (T == diph_l_(n))
                    {
                      return Nom("ijk_l_") + Nom(i) + Nom("_") + Nom(j) + Nom("_") + Nom(k);
                    }
                  if (T == diph_v_(n))
                    {
                      return Nom("ijk_v_") + Nom(i) + Nom("_") + Nom(j) + Nom("_") + Nom(k);
                    }
                }
            }
        }
    }
  return Nom("_not_here_");
}

template<>
void Cut_field_template<double,ArrOfDouble>::multiply_by_scalar(double scalar_l, double scalar_v)
{
  const int ni = IJK_Field_template<double,ArrOfDouble>::ni();
  const int nj = IJK_Field_template<double,ArrOfDouble>::nj();
  const int nk = IJK_Field_template<double,ArrOfDouble>::nk();
  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              int n = cut_cell_disc_->get_n(i,j,k);
              if (n < 0)
                {
                  int phase = IJK_Interfaces::convert_indicatrice_to_phase(cut_cell_disc_->indic_pure(i,j,k));
                  pure_(i,j,k) *= (phase == 0) ? scalar_v : scalar_l;
                }
              else
                {
                  diph_l_(n) *= scalar_l;
                  diph_v_(n) *= scalar_v;
                }
            }
        }
    }
}

template<>
void Cut_field_template<double,ArrOfDouble>::divide_by_scalar(double scalar_l, double scalar_v)
{
  multiply_by_scalar(1./scalar_l, 1./scalar_v);
}

template<class T, int N>
Cut_field_vector<T, N>::Cut_field_vector()
{
}

template<class T, int N>
void Cut_field_vector<T, N>::set_to_uniform_value(int valeur)
{
  for (int i = 0; i < N; i++)
    {
      Cut_field_vector<T, N>::operator[](i).set_to_uniform_value(valeur);
    }
}

template<class T, int N>
void Cut_field_vector<T, N>::echange_espace_virtuel(int ghost)
{
  for (int i = 0; i < N; i++)
    {
      Cut_field_vector<T, N>::operator[](i).echange_espace_virtuel(ghost);
    }
}

template<class T, int N>
int Cut_field_vector<T, N>::ghost()
{
  for (int i = 0; i < N; i++)
    {
      assert((Cut_field_vector<T, N>::operator[](i).ghost() == Cut_field_vector<T, N>::operator[](0).ghost()));
    }
  return Cut_field_vector<T, N>::operator[](0).ghost();
}

template<class T, int N>
void Cut_field_vector<T, N>::echange_pure_vers_diph_cellules_initialement_pures()
{
  for (int i = 0; i < N; i++)
    {
      Cut_field_vector<T, N>::operator[](i).echange_pure_vers_diph_cellules_initialement_pures();
    }
}

template<class T, int N>
void Cut_field_vector<T, N>::echange_diph_vers_pure_cellules_finalement_pures()
{
  for (int i = 0; i < N; i++)
    {
      Cut_field_vector<T, N>::operator[](i).echange_diph_vers_pure_cellules_finalement_pures();
    }
}

template<class T, int N>
void Cut_field_vector<T, N>::vide_phase_invalide_cellules_diphasiques()
{
  for (int i = 0; i < N; i++)
    {
      Cut_field_vector<T, N>::operator[](i).vide_phase_invalide_cellules_diphasiques();
    }
}


// Explicit instantiations

template class Cut_field_template<int,ArrOfInt>;
template class Cut_field_template<double,ArrOfDouble>;

template class Cut_field_vector<int, 3>;
template class Cut_field_vector<double, 3>;


