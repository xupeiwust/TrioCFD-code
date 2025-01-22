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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Facettes_Interp_FT.cpp
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#include <Facettes_Interp_FT.h>
#include <Intersection_Interface_ijk.h>
#include <IJK_Navier_Stokes_tools_cut_cell.h>

Implemente_instanciable_sans_constructeur(Facettes_Interp_FT, "Facettes_Interp_FT", Objet_U) ;

Facettes_Interp_FT::Facettes_Interp_FT()
{
}

Sortie& Facettes_Interp_FT::printOn(Sortie& os) const
{
  return os;
}

Entree& Facettes_Interp_FT::readOn(Entree& is)
{
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades(is);
  return is;
}

void Facettes_Interp_FT::set_param(Param& param)
{
  param.ajouter("number_of_interpolation_points", &number_of_interpolation_points_);
  param.ajouter("scaled_distance_interpolation_1", &scaled_distance_interpolation_1_);
  param.ajouter("scaled_distance_interpolation_2", &scaled_distance_interpolation_2_);
}

void Facettes_Interp_FT::associer(const IJK_Interfaces& interfaces, const Cut_cell_FT_Disc& cut_cell_disc, const IJK_Splitting& splitting_ft, const Maillage_FT_IJK& maillage_ft_ijk, const Maillage_FT_IJK& old_maillage_ft_ijk)
{
  ref_interfaces_ = interfaces;
  ref_cut_cell_disc_ = cut_cell_disc;
  ref_splitting_ = splitting_ft;
  ref_maillage_ft_ijk_ = maillage_ft_ijk;
  ref_old_maillage_ft_ijk_ = old_maillage_ft_ijk;
}

void Facettes_Interp_FT::cut_cell_perform_interpolation_facettes_next(int old_en_premier /* next() */)
{
  // Synchronisation de old_en_premier_ avec le old_en_premier_ de IJK_Interfaces
  old_en_premier_ = old_en_premier;
  assert(next() == ref_interfaces_->next());

  cut_cell_perform_interpolation_facettes(true,
                                          ref_cut_cell_disc_,
                                          ref_splitting_->get_grid_geometry(),
                                          ref_maillage_ft_ijk_.valeur(),
                                          interpolation_coord_[next()],
                                          signed_independent_index_[next()],
                                          coefficient_[next()]);
}

void Facettes_Interp_FT::cut_cell_perform_interpolation_facettes_old(bool not_old_en_premier /* old() */)
{
  // Synchronisation de old_en_premier_ avec le old_en_premier_ de IJK_Interfaces
  old_en_premier_ = not not_old_en_premier;
  assert(old() == ref_interfaces_->old());

  cut_cell_perform_interpolation_facettes(false,
                                          ref_cut_cell_disc_,
                                          ref_splitting_->get_grid_geometry(),
                                          ref_old_maillage_ft_ijk_.valeur(),
                                          interpolation_coord_[old()],
                                          signed_independent_index_[old()],
                                          coefficient_[old()]);
}

void Facettes_Interp_FT::cut_cell_perform_interpolation_facettes(bool next_time,
                                                                 const Cut_cell_FT_Disc& cut_cell_disc,
                                                                 const IJK_Grid_Geometry& geom,
                                                                 const Maillage_FT_IJK& maillage,
                                                                 FixedVector<DoubleTabFT, 4>& interpolation_coord,
                                                                 FixedVector<IntTabFT, 4>& interpolation_signed_independent_index,
                                                                 FixedVector<DoubleTabFT, 4>& interpolation_coefficient)
{
  const double dist_1 = get_distance_interpolation_1();
  const double dist_2 = get_distance_interpolation_2();
  const DoubleTab& normale_facettes = maillage.get_update_normale_facettes();

  if (number_of_interpolation_points_ == 1)
    {
      calcul_coefficient_interpolation_facette_cut_cell(next_time,
                                                        dist_1,
                                                        cut_cell_disc,
                                                        maillage,
                                                        interpolation_coord,
                                                        normale_facettes,
                                                        interpolation_signed_independent_index,
                                                        interpolation_coefficient);
    }
  else if (number_of_interpolation_points_ == 2)
    {
      calcul_coefficient_interpolation_facette_cut_cell_second_order(next_time,
                                                                     dist_1,
                                                                     dist_2,
                                                                     cut_cell_disc,
                                                                     maillage,
                                                                     interpolation_coord,
                                                                     normale_facettes,
                                                                     interpolation_signed_independent_index,
                                                                     interpolation_coefficient);
    }
  else
    {
      Cerr << "Nombre de points d'interpolation non reconnu." << finl;
      Process::exit();
    }
}

void Facettes_Interp_FT::calcul_coefficient_interpolation_facette_cut_cell(
  bool next_time,
  const double dist,
  const Cut_cell_FT_Disc& cut_cell_disc,
  const Maillage_FT_IJK& maillage,
  FixedVector<DoubleTabFT, 4>& interpolation_coord,
  const DoubleTab& normal_on_interf,
  FixedVector<IntTabFT, 4>& cut_cell_facettes_interpolation_signed_independent_index,
  FixedVector<DoubleTabFT, 4>& cut_cell_facettes_interpolation_coefficient)
{
  const int nb_facettes = maillage.nb_facettes();
  interpolation_coord[1].resize(nb_facettes, 3);
  interpolation_coord[0].resize(nb_facettes, 3);
  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    {
      Vecteur3 coords_fa7 = maillage.coords_fa7(fa7);
      Vecteur3 normal(normal_on_interf, fa7);

      Vecteur3 interp_1 = Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(coords_fa7, normal, dist);
      interpolation_coord[1](fa7, 0) = interp_1[0];
      interpolation_coord[1](fa7, 1) = interp_1[1];
      interpolation_coord[1](fa7, 2) = interp_1[2];

      Vecteur3 interp_0 = Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(coords_fa7, normal, -dist);
      interpolation_coord[0](fa7, 0) = interp_0[0];
      interpolation_coord[0](fa7, 1) = interp_0[1];
      interpolation_coord[0](fa7, 2) = interp_0[2];
    }

  cut_cell_facettes_interpolation_signed_independent_index[1].resize(interpolation_coord[1].dimension(0), 4);
  cut_cell_facettes_interpolation_coefficient[1].resize(interpolation_coord[1].dimension(0), 4);
  ijk_interpolate_cut_cell_skip_unknown_points(next_time, 1, cut_cell_disc, interpolation_coord[1], cut_cell_facettes_interpolation_signed_independent_index[1], cut_cell_facettes_interpolation_coefficient[1], 1.e31);

  cut_cell_facettes_interpolation_signed_independent_index[0].resize(interpolation_coord[0].dimension(0), 4);
  cut_cell_facettes_interpolation_coefficient[0].resize(interpolation_coord[0].dimension(0), 4);
  ijk_interpolate_cut_cell_skip_unknown_points(next_time, 0, cut_cell_disc, interpolation_coord[0], cut_cell_facettes_interpolation_signed_independent_index[0], cut_cell_facettes_interpolation_coefficient[0], 1.e31);
}

void Facettes_Interp_FT::calcul_coefficient_interpolation_facette_cut_cell_second_order(
  bool next_time,
  const double dist_1,
  const double dist_2,
  const Cut_cell_FT_Disc& cut_cell_disc,
  const Maillage_FT_IJK& maillage,
  FixedVector<DoubleTabFT, 4>& interpolation_coord,
  const DoubleTab& normal_on_interf,
  FixedVector<IntTabFT, 4>& cut_cell_facettes_interpolation_signed_independent_index,
  FixedVector<DoubleTabFT, 4>& cut_cell_facettes_interpolation_coefficient)
{
  const int nb_facettes = maillage.nb_facettes();
  interpolation_coord[1].resize(nb_facettes, 3);
  interpolation_coord[0].resize(nb_facettes, 3);
  interpolation_coord[3].resize(nb_facettes, 3);
  interpolation_coord[2].resize(nb_facettes, 3);
  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    {
      Vecteur3 coords_fa7 = maillage.coords_fa7(fa7);
      Vecteur3 normal(normal_on_interf, fa7);

      Vecteur3 interp_1 = Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(coords_fa7, normal, dist_1);
      interpolation_coord[1](fa7, 0) = interp_1[0];
      interpolation_coord[1](fa7, 1) = interp_1[1];
      interpolation_coord[1](fa7, 2) = interp_1[2];

      Vecteur3 interp_0 = Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(coords_fa7, normal, -dist_1);
      interpolation_coord[0](fa7, 0) = interp_0[0];
      interpolation_coord[0](fa7, 1) = interp_0[1];
      interpolation_coord[0](fa7, 2) = interp_0[2];

      Vecteur3 interp_3 = Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(coords_fa7, normal, dist_2);
      interpolation_coord[3](fa7, 0) = interp_3[0];
      interpolation_coord[3](fa7, 1) = interp_3[1];
      interpolation_coord[3](fa7, 2) = interp_3[2];

      Vecteur3 interp_2 = Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(coords_fa7, normal, -dist_2);
      interpolation_coord[2](fa7, 0) = interp_2[0];
      interpolation_coord[2](fa7, 1) = interp_2[1];
      interpolation_coord[2](fa7, 2) = interp_2[2];
    }

  cut_cell_facettes_interpolation_signed_independent_index[1].resize(interpolation_coord[1].dimension(0), 4);
  cut_cell_facettes_interpolation_coefficient[1].resize(interpolation_coord[1].dimension(0), 4);
  ijk_interpolate_cut_cell_skip_unknown_points(next_time, 1, cut_cell_disc, interpolation_coord[1], cut_cell_facettes_interpolation_signed_independent_index[1], cut_cell_facettes_interpolation_coefficient[1], 1.e31);

  cut_cell_facettes_interpolation_signed_independent_index[0].resize(interpolation_coord[0].dimension(0), 4);
  cut_cell_facettes_interpolation_coefficient[0].resize(interpolation_coord[0].dimension(0), 4);
  ijk_interpolate_cut_cell_skip_unknown_points(next_time, 0, cut_cell_disc, interpolation_coord[0], cut_cell_facettes_interpolation_signed_independent_index[0], cut_cell_facettes_interpolation_coefficient[0], 1.e31);

  cut_cell_facettes_interpolation_signed_independent_index[3].resize(interpolation_coord[3].dimension(0), 4);
  cut_cell_facettes_interpolation_coefficient[3].resize(interpolation_coord[3].dimension(0), 4);
  ijk_interpolate_cut_cell_skip_unknown_points(next_time, 1, cut_cell_disc, interpolation_coord[3], cut_cell_facettes_interpolation_signed_independent_index[3], cut_cell_facettes_interpolation_coefficient[3], 1.e31);

  cut_cell_facettes_interpolation_signed_independent_index[2].resize(interpolation_coord[2].dimension(0), 4);
  cut_cell_facettes_interpolation_coefficient[2].resize(interpolation_coord[2].dimension(0), 4);
  ijk_interpolate_cut_cell_skip_unknown_points(next_time, 0, cut_cell_disc, interpolation_coord[2], cut_cell_facettes_interpolation_signed_independent_index[2], cut_cell_facettes_interpolation_coefficient[2], 1.e31);
}

double Facettes_Interp_FT::get_distance_interpolation_1() const
{
  const IJK_Grid_Geometry& geom = ref_splitting_->get_grid_geometry();
  const double dist_1 = scaled_distance_interpolation_1_/sqrt(3) * std::pow(std::pow(geom.get_constant_delta(DIRECTION_I), 2.) +
                                                                            std::pow(geom.get_constant_delta(DIRECTION_J), 2.) +
                                                                            std::pow(geom.get_constant_delta(DIRECTION_K), 2.),
                                                                            0.5);
  return dist_1;
}

double Facettes_Interp_FT::get_distance_interpolation_2() const
{
  const IJK_Grid_Geometry& geom = ref_splitting_->get_grid_geometry();
  const double dist_2 = scaled_distance_interpolation_2_/sqrt(3) * std::pow(std::pow(geom.get_constant_delta(DIRECTION_I), 2.) +
                                                                            std::pow(geom.get_constant_delta(DIRECTION_J), 2.) +
                                                                            std::pow(geom.get_constant_delta(DIRECTION_K), 2.),
                                                                            0.5);
  return dist_2;
}

