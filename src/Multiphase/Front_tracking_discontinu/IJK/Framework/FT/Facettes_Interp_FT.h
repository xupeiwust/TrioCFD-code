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
// File      : Facettes_Interp_FT.h
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Facettes_Interp_FT_included
#define Facettes_Interp_FT_included

#include <FixedVector.h>
#include <Maillage_FT_IJK.h>
#include <TRUSTTabFT.h>
#include <Cut_cell_FT_Disc.h>
#include <Objet_U.h>

class IJK_Interfaces;
class Cut_cell_FT_Disc;
class IJK_Splitting;

class Facettes_Interp_FT : public Objet_U
{
  Declare_instanciable(Facettes_Interp_FT);

public:
  void associer(const IJK_Interfaces& interfaces, const Cut_cell_FT_Disc& cut_cell_disc, const IJK_Splitting& splitting_ft, const Maillage_FT_IJK& maillage_ft_ijk, const Maillage_FT_IJK& old_maillage_ft_ijk);
  void set_param(Param& param);

  void cut_cell_perform_interpolation_facettes_next(int old_en_premier /* next() */);

  void cut_cell_perform_interpolation_facettes_old(bool not_old_en_premier /* old() */);

  const FixedVector<IntTabFT, 4>& get_signed_independent_index_next() const { return signed_independent_index_[next()]; }
  const FixedVector<IntTabFT, 4>& get_signed_independent_index_old() const { return signed_independent_index_[old()]; }
  const FixedVector<DoubleTabFT, 4>& get_coefficient_next() const { return coefficient_[next()]; }
  const FixedVector<DoubleTabFT, 4>& get_coefficient_old() const { return coefficient_[old()]; }

  int old() const { return 1 - old_en_premier_; }
  int next() const { return old_en_premier_; }

  int get_number_of_interpolation_points() const { return number_of_interpolation_points_; }
  double get_scaled_distance_interpolation_1() const { return scaled_distance_interpolation_1_; }
  double get_scaled_distance_interpolation_2() const { return scaled_distance_interpolation_2_; }

  const Maillage_FT_IJK& maillage_ft_ijk() const { return ref_maillage_ft_ijk_; }
  const Maillage_FT_IJK& old_maillage_ft_ijk() const { return ref_old_maillage_ft_ijk_; }

  double get_distance_interpolation_1() const;
  double get_distance_interpolation_2() const;

protected:

  void cut_cell_perform_interpolation_facettes(bool next_time,
                                               const Cut_cell_FT_Disc& cut_cell_disc,
                                               const IJK_Grid_Geometry& geom,
                                               const Maillage_FT_IJK& maillage,
                                               FixedVector<DoubleTabFT, 4>& interpolation_coord,
                                               FixedVector<IntTabFT, 4>& interpolation_signed_independent_index,
                                               FixedVector<DoubleTabFT, 4>& interpolation_coefficient);

  void calcul_coefficient_interpolation_facette_cut_cell(bool next_time,
                                                         const double dist,
                                                         const Cut_cell_FT_Disc& cut_cell_disc,
                                                         const Maillage_FT_IJK& maillage,
                                                         FixedVector<DoubleTabFT, 4>& interpolation_coord,
                                                         const DoubleTab& normal_on_interf,
                                                         FixedVector<IntTabFT, 4>& cut_cell_facettes_interpolation_signed_independent_index,
                                                         FixedVector<DoubleTabFT, 4>& cut_cell_facettes_interpolation_coefficient);

  void calcul_coefficient_interpolation_facette_cut_cell_second_order(bool next_time,
                                                                      const double dist_1,
                                                                      const double dist_2,
                                                                      const Cut_cell_FT_Disc& cut_cell_disc,
                                                                      const Maillage_FT_IJK& maillage,
                                                                      FixedVector<DoubleTabFT, 4>& interpolation_coord,
                                                                      const DoubleTab& normal_on_interf,
                                                                      FixedVector<IntTabFT, 4>& cut_cell_facettes_interpolation_signed_independent_index,
                                                                      FixedVector<DoubleTabFT, 4>& cut_cell_facettes_interpolation_coefficient);

  // Stockage des indices et coefficients des points d'interpolation a une certaine distance des facettes de l'interface
  // Chaque tableau FixedVector<___, 4> est compose de la facon suivante :
  //     [0] Premier point d'interpolation dans la phase 0
  //     [1] Premier point d'interpolation dans la phase 1
  //     [2] Second point d'interpolation dans la phase 0
  //     [3] Second point d'interpolation dans la phase 1
  // Les tableaux IntTabFT/DoubleTabFT ont la dimension (nombre_de_facettes, nombre_de_points_utilises_par_l_interpolation_pour_chaque_facette).
  FixedVector<FixedVector<IntTabFT, 4>, 2> signed_independent_index_;
  FixedVector<FixedVector<DoubleTabFT, 4>, 2> coefficient_;

  // Intermediate variables
  FixedVector<FixedVector<DoubleTabFT, 4>, 2> interpolation_coord_;

  bool old_en_premier_ = true; // Doit etre synchronise avec old_en_premier_ dans IJK_Interfaces
  OBS_PTR(IJK_Interfaces) ref_interfaces_;
  OBS_PTR(Cut_cell_FT_Disc) ref_cut_cell_disc_;
  OBS_PTR(IJK_Splitting) ref_splitting_;
  OBS_PTR(Maillage_FT_IJK) ref_maillage_ft_ijk_;
  OBS_PTR(Maillage_FT_IJK) ref_old_maillage_ft_ijk_;

  int number_of_interpolation_points_ = 1;
  double scaled_distance_interpolation_1_ = 1.0; // Distance a l'interface du premier point d'interpolation
  double scaled_distance_interpolation_2_ = 2.0; // Distance a l'interface du second point d'interpolation
};

#endif /* Facettes_Interp_FT_included */
