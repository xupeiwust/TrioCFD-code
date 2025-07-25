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

#ifndef OpConvQuickIJKScalar_cut_cell_TPP_included
#define OpConvQuickIJKScalar_cut_cell_TPP_included

#include <IJK_Field.h>
#include <Cut_cell_convection_auxiliaire.h>
#include <simd_tools.h>

template <DIRECTION _DIR_>
void OpConvQuickIJKScalar_cut_cell_double::compute_flux_(IJK_Field_local_double& resu, const int k_layer)
{
  if (_DIR_==DIRECTION::Z)
    {
      // The previous layer of curv and fram values might have been computed already
      if (stored_curv_fram_layer_z_ != k_layer-1)
        {
          compute_curv_fram(_DIR_, k_layer-1);
          shift_curv_fram(tmp_curv_fram_);
        }
    }
  compute_curv_fram(_DIR_, k_layer);

  ConstIJK_double_ptr velocity_dir(get_input_velocity(_DIR_), 0, 0, k_layer);
  ConstIJK_double_ptr input_field(*input_field_, 0, 0, k_layer);
  ConstIJK_double_ptr curv_values(tmp_curv_fram_, 0, 0, 1); /* if z direction, "left" will be in layer 0 */
  ConstIJK_double_ptr fram_values(tmp_curv_fram_, 0, 0, 3); /* if z direction, "left" is in layer 2 */
  IJK_double_ptr resu_ptr(resu, 0, 0, 0);
  const int nx = _DIR_==DIRECTION::X ? input_field_->ni() + 1 : input_field_->ni();
  const int ny = _DIR_==DIRECTION::Y ? input_field_->nj() + 1 : input_field_->nj();

  const double delta_xyz = _DIR_==DIRECTION::Z ? (channel_data_.get_delta_z()[k_layer-1] + channel_data_.get_delta_z()[k_layer]) * 0.5 : channel_data_.get_delta((int)_DIR_);
  const double surface = channel_data_.get_surface(k_layer, 1, (int)_DIR_);
  if (_DIR_==DIRECTION::Z)
    {
      // Are we on the wall ?
      const int global_k_layer = k_layer + channel_data_.offset_to_global_k_layer();
      // global index of the layer of flux of the wall
      //  (assume one walls at zmin and zmax)
      const int first_global_k_layer = 0;
      const int last_global_k_layer = channel_data_.nb_elem_k_tot();

      // GB 21/12/2020 : Similarly to velocity, we make the same adjustments.
      // if (global_k_layer == first_global_k_layer || global_k_layer == last_global_k_layer) {
      // ie (i) replace the former "==" by "<=" and ">=" to be identic to the condition in OpCentre4IJK
      // and (ii) add the condition on perio_k
      if (!perio_k_ && (global_k_layer <= first_global_k_layer || global_k_layer >= last_global_k_layer))
        {
          // We are on (or worse inside) the wall, zero flux
          putzero(resu);
          return;
        }
    }

  const double delta_xyz_squared_over_8 = delta_xyz * delta_xyz * 0.125;
  const int imax = nx;
  const int jmax = ny;
  const int vsize = Simd_double::size();
  const Simd_double zero = 0.;
  for (int j = 0; ; j++)
    {
      for (int i = 0; i < imax; i += vsize)
        {
          Simd_double velocity;
          velocity_dir.get_center(i, velocity);
          Simd_double T0, T1; // scalar value at left and at right of the computed flux
          Simd_double fram0, fram1;
          Simd_double curv0, curv1;
          input_field.get_left_center(_DIR_, i, T0, T1);
          fram_values.get_left_center(_DIR_, i, fram0, fram1);
          curv_values.get_left_center(_DIR_, i, curv0, curv1);
          Simd_double fram = max(fram0, fram1);
          Simd_double curv	   = select_double(velocity, zero, curv1, curv0);
          Simd_double T_amont = select_double(velocity, zero, T1 /* if velocity < 0 */, T0 /* if velocity > 0 */);
          Simd_double flux	   = (T0 + T1) * 0.5 - delta_xyz_squared_over_8 * curv;
          flux		   = ((1. - fram) * flux + fram * T_amont) * velocity * surface;
          resu_ptr.put_val(i, flux);
        }
      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      input_field.next_j();
      velocity_dir.next_j();
      fram_values.next_j();
      curv_values.next_j();
      resu_ptr.next_j();
    }

  if (_DIR_==DIRECTION::Z)
    {
      // store curv and fram for next layer of fluxes in z direction
      shift_curv_fram(tmp_curv_fram_);
      stored_curv_fram_layer_z_ = k_layer;
    }
}

//template <DIRECTION _DIR_>
//void OpConvQuickIJKScalar_cut_cell_double::compute_flux_(IJK_Field_local_double& resu, const int k_layer)
//{
//  const int nx = _DIR_==DIRECTION::X ? input_field_->ni() + 1 : input_field_->ni();
//  const int ny = _DIR_==DIRECTION::Y ? input_field_->nj() + 1 : input_field_->nj();
//
//  const int imax = nx;
//  const int jmax = ny;
//  for (int j = 0; j < jmax; j++)
//    {
//      for (int i = 0; i < imax; i++)
//        {
//          double flux = compute_flux_local_<_DIR_>(i,j,k_layer);
//          resu(i,j,0) = flux;
//        }
//    }
//}

template <DIRECTION _DIR_>
double OpConvQuickIJKScalar_cut_cell_double::compute_flux_local_(int i, int j, int k)
{
  const int dir_i = (_DIR_ == DIRECTION::X);
  const int dir_j = (_DIR_ == DIRECTION::Y);
  const int dir_k = (_DIR_ == DIRECTION::Z);

  const IJK_Field_local_double& velocity_dir = get_input_velocity(_DIR_);
  const IJK_Field_local_double& input_field = *input_field_;

  const double delta_xyz = _DIR_==DIRECTION::Z ? (channel_data_.get_delta_z()[k-1] + channel_data_.get_delta_z()[k]) * 0.5 : channel_data_.get_delta((int)_DIR_);
  const double surface = channel_data_.get_surface(k, 1, (int)_DIR_);

  if (flux_determined_by_wall_<_DIR_>(k))
    {
      return 0;
    }
  else
    {
      double velocity = velocity_dir(i,j,k);
      double input_left_left = input_field(i-2*dir_i,j-2*dir_j,k-2*dir_k);
      double input_left = input_field(i-dir_i,j-dir_j,k-dir_k);
      double input_centre = input_field(i,j,k);
      double input_right = input_field(i+dir_i,j+dir_j,k+dir_k);

      double flux = OpConvQuickIJKScalar_cut_cell_double::compute_flux_local_<_DIR_>(k, delta_xyz, surface, velocity, input_left_left, input_left, input_centre, input_right);
      return flux;
    }
}

template <DIRECTION _DIR_>
double OpConvQuickIJKScalar_cut_cell_double::compute_flux_local_(int k_layer, double delta_xyz, double surface, double velocity, double input_left_left, double input_left, double input_centre, double input_right)
{
  const int dir_k = (_DIR_ == DIRECTION::Z);

  const double delta_xyz_squared_over_8 = delta_xyz * delta_xyz * 0.125;

  Vecteur3 curv_fram_left = compute_curv_fram_local_<_DIR_>(k_layer-dir_k, input_left_left, input_left, input_centre);
  double curv_left = curv_fram_left[0];
  double fram_left = curv_fram_left[1];

  Vecteur3 curv_fram_centre = compute_curv_fram_local_<_DIR_>(k_layer, input_left, input_centre, input_right);
  double curv_centre = curv_fram_centre[0];
  double fram_centre = curv_fram_centre[1];

  double fram = std::max(fram_left, fram_centre);
  double curv = velocity < 0. ? curv_centre : curv_left;
  double T_amont = velocity < 0. ? input_centre /* if velocity < 0 */ : input_left /* if velocity > 0 */;
  double flux = (input_left + input_centre) * 0.5 - delta_xyz_squared_over_8 * curv;
  flux = ((1. - fram) * flux + fram * T_amont) * velocity * surface;
  return flux;
}

template <DIRECTION _DIR_>
double OpConvQuickIJKScalar_cut_cell_double::compute_flux_local_(double surface, double velocity, double input)
{
  double flux = input * velocity * surface;
  return flux;
}

template <DIRECTION _DIR_>
bool OpConvQuickIJKScalar_cut_cell_double::flux_determined_by_wall_(int k)
{
  if (_DIR_==DIRECTION::Z)
    {
      // Are we on the wall ?
      const int global_k_layer = k + channel_data_.offset_to_global_k_layer();
      // global index of the layer of flux of the wall
      //  (assume one walls at zmin and zmax)
      const int first_global_k_layer = 0;
      const int last_global_k_layer = channel_data_.nb_elem_k_tot();

      // GB 21/12/2020 : Similarly to velocity, we make the same adjustments.
      // if (global_k_layer == first_global_k_layer || global_k_layer == last_global_k_layer) {
      // ie (i) replace the former "==" by "<=" and ">=" to be identic to the condition in OpCentre4IJK
      // and (ii) add the condition on perio_k
      if (!perio_k_ && (global_k_layer <= first_global_k_layer || global_k_layer >= last_global_k_layer))
        {
          // We are on (or worse inside) the wall, zero flux
          return 1;
        }
    }
  return 0;
}

template <DIRECTION _DIR_>
Vecteur3 OpConvQuickIJKScalar_cut_cell_double::compute_curv_fram_local_(int k_layer, double input_left, double input_centre, double input_right)
{
  if (_DIR_ == DIRECTION::Z)
    {
      // Are we on the wall ?
      const int global_k_layer = k_layer + channel_data_.offset_to_global_k_layer();
      // global index of the layer of flux of the wall
      //  (assume one walls at zmin and zmax)
      const int first_global_k_layer = 0;
      const int last_global_k_layer = channel_data_.nb_elem_k_tot();

      if (!perio_k_ && (global_k_layer <= first_global_k_layer || global_k_layer >= last_global_k_layer))
        {
          double curv = 0.;
          double fram = 1.;
          return {curv, fram, 0.};
        }
    }

  double factor01 = -1.0;
  double factor12 = -1.0;
  if (_DIR_==DIRECTION::X)
    {
      // Compute invd1 and invd2 factors:
      const double inv_dx = 1./channel_data_.get_delta_x();
      factor01 = inv_dx * inv_dx;
      factor12 = factor01;
    }
  if (_DIR_==DIRECTION::Y)
    {
      // Compute invd1 and invd2 factors:
      const double inv_dy = 1. / channel_data_.get_delta_y();
      factor01 = inv_dy * inv_dy;
      factor12 = factor01;
    }
  if (_DIR_==DIRECTION::Z)
    {
      // Compute invd1 and invd2 factors:
      const double dz0 = channel_data_.get_delta_z()[k_layer-1];
      const double dz1 = channel_data_.get_delta_z()[k_layer];
      const double dz2 = channel_data_.get_delta_z()[k_layer+1];
      factor01 = 1. / (dz1 * (dz0 + dz1) * 0.5);
      factor12 = 1. / (dz1 * (dz1 + dz2) * 0.5);
    }

  double curv = (input_right - input_centre) * factor12 - (input_centre - input_left) * factor01;
  double smin = std::min(input_left, input_right);
  double smax = std::max(input_left, input_right);
  // Compared to original code (Eval_Quick_VDF_Elem.h), we first compute the 4th power,
  // then take the max (dabs is then useless)
  double dsabs = 0. < smax - smin ? smax - smin : smin - smax;
  double ds = dsabs < DMINFLOAT ? 1. : smax - smin;
  double sr = dsabs < DMINFLOAT ? 0. : ((input_centre - smin) / ds - 0.5) * 2.;
  sr *= sr;
  sr *= sr;
  sr = std::min(sr, 1.);
  double fram = sr;
  return {curv, fram, 0.};
}

template <DIRECTION _DIR_>
void OpConvQuickIJKScalar_cut_cell_double::correct_flux_(IJK_Field_local_double *const flux, int k_layer)
{
  const int dir = static_cast<int>(_DIR_);

  const Cut_field_double& velocity_dir = static_cast<const Cut_field_double&>(get_input_velocity(_DIR_));
  const Cut_field_double& input_cut_field = static_cast<const Cut_field_double&>(*input_field_);
  assert(&(*cut_cell_flux_)[0].get_cut_cell_disc() == &(*cut_cell_flux_)[1].get_cut_cell_disc());
  assert(&(*cut_cell_flux_)[0].get_cut_cell_disc() == &(*cut_cell_flux_)[2].get_cut_cell_disc());
  const Cut_cell_FT_Disc& cut_cell_disc = (*cut_cell_flux_)[0].get_cut_cell_disc();

  IJK_Field_int& treatment_count = *treatment_count_;
  int& new_treatment = *new_treatment_;
  new_treatment += 1;

  int backward_receptive_stencil = 1;
  int forward_receptive_stencil = 2;
  assert(backward_receptive_stencil <= cut_cell_disc.get_ghost_size());
  assert(forward_receptive_stencil <= cut_cell_disc.get_ghost_size());

  if (_DIR_ == DIRECTION::Z)
    {
      if (cut_cell_disc.get_domaine().get_periodic_flag(dir))
        {
          const int kmax = cut_cell_disc.get_interfaces().I().nk();
          int n_dir = cut_cell_disc.get_domaine().get_nb_elem_local(dir);
          int n_dir_tot = cut_cell_disc.get_domaine().get_nb_elem_tot(dir);

          // Le processeur contient deux fois les valeurs sur les bords
          if (n_dir == n_dir_tot)
            {
              if (k_layer == kmax)
                {
                  k_layer = 0;
                }
            }
        }
    }

  int min_k = (_DIR_ == DIRECTION::Z) ? k_layer-forward_receptive_stencil : k_layer;
  int max_k = (_DIR_ == DIRECTION::Z) ? k_layer+backward_receptive_stencil : k_layer;
  for (int k_c = min_k; k_c <= max_k; k_c++)
    {
      for (int index = cut_cell_disc.get_k_value_index(k_c); index < cut_cell_disc.get_k_value_index(k_c+1); index++)
        {
          int n = cut_cell_disc.get_n_from_k_index(index);

          Int3 ijk_no_per = cut_cell_disc.get_ijk(n);
          assert(k_c == ijk_no_per[2]);

          int min_decalage = (_DIR_ == DIRECTION::Z) ? (k_layer - k_c) : -backward_receptive_stencil;
          int max_decalage = (_DIR_ == DIRECTION::Z) ? (k_layer - k_c) : forward_receptive_stencil;
          for (int decalage = min_decalage; decalage <= max_decalage; decalage++)
            {
              int i = ijk_no_per[0] + (_DIR_ == DIRECTION::X)*decalage;
              int j = ijk_no_per[1] + (_DIR_ == DIRECTION::Y)*decalage;
              int k = ijk_no_per[2] + (_DIR_ == DIRECTION::Z)*decalage;
              assert(k_layer == k);

              if (!cut_cell_disc.get_domaine().within_ghost_<dir>(i, j, k, 0, 1))
                continue;

              if (treatment_count(i,j,k) == new_treatment)
                continue;

              {
                int index_ijk_per = 0;
                while (index_ijk_per >= 0)
                  {
                    Int3 ijk = cut_cell_disc.ijk_per_of_index(i, j, k, index_ijk_per);
                    index_ijk_per = cut_cell_disc.next_index_ijk_per(i, j, k, index_ijk_per, 0, 1);

                    treatment_count(ijk[0],ijk[1],ijk[2]) = new_treatment;
                  }
              }

              double indicatrice_centre = cut_cell_disc.get_interfaces().I(i,j,k);
              double next_indicatrice_centre = cut_cell_disc.get_interfaces().In(i,j,k);
              int n_centre = cut_cell_disc.get_n(i, j, k);

              if (flux_determined_by_wall_<_DIR_>(k))
                {
                  // Le flux est dans ce cas determine par la paroi
                  // Aucune correction n'est donc necessaire
                  assert(n_centre < 0); // Le cas d'une cellule diphasique avec flux de paroi n'est pas traite
                }
              else
                {
                  assert((n_centre >= 0) || ((indicatrice_centre == 0.) || (indicatrice_centre == 1.)));
                  bool centre_monophasique = IJK_Interfaces::est_pure(.5*(indicatrice_centre + next_indicatrice_centre));

                  int phase_min = (int)std::floor(.5*(indicatrice_centre + next_indicatrice_centre));
                  int phase_max = (int)std::ceil(.5*(indicatrice_centre + next_indicatrice_centre));
                  for (int phase = phase_min ; phase <= phase_max ; phase++)
                    {
                      const DoubleTabFT_cut_cell& diph_input = (phase == 0) ? input_cut_field.diph_v_ : input_cut_field.diph_l_;
                      const DoubleTabFT_cut_cell& diph_velocity = (phase == 0) ? velocity_dir.diph_v_ : velocity_dir.diph_l_;
                      DoubleTabFT_cut_cell& diph_flux = (phase == 0) ? (*cut_cell_flux_)[dir].diph_v_ : (*cut_cell_flux_)[dir].diph_l_;

                      const int dir_i = (_DIR_ == DIRECTION::X);
                      const int dir_j = (_DIR_ == DIRECTION::Y);
                      const int dir_k = (_DIR_ == DIRECTION::Z);

                      int n_left_left = cut_cell_disc.get_n(i-2*dir_i,j-2*dir_j,k-2*dir_k);
                      int n_left = cut_cell_disc.get_n(i-dir_i,j-dir_j,k-dir_k);
                      int n_right = cut_cell_disc.get_n(i+dir_i,j+dir_j,k+dir_k);

                      double indicatrice_left_left = cut_cell_disc.get_interfaces().I(i-2*dir_i,j-2*dir_j,k-2*dir_k);
                      double indicatrice_left = cut_cell_disc.get_interfaces().I(i-dir_i,j-dir_j,k-dir_k);
                      double indicatrice_right = cut_cell_disc.get_interfaces().I(i+dir_i,j+dir_j,k+dir_k);

                      double next_indicatrice_left_left = cut_cell_disc.get_interfaces().In(i-2*dir_i,j-2*dir_j,k-2*dir_k);
                      double next_indicatrice_left = cut_cell_disc.get_interfaces().In(i-dir_i,j-dir_j,k-dir_k);
                      double next_indicatrice_right = cut_cell_disc.get_interfaces().In(i+dir_i,j+dir_j,k+dir_k);

                      bool left_left_monophasique = IJK_Interfaces::est_pure(.5*(indicatrice_left_left + next_indicatrice_left_left));
                      bool left_monophasique = IJK_Interfaces::est_pure(.5*(indicatrice_left + next_indicatrice_left));
                      bool right_monophasique = IJK_Interfaces::est_pure(.5*(indicatrice_right + next_indicatrice_right));
                      bool phase_invalide_left_left = (left_left_monophasique && (phase != IJK_Interfaces::convert_indicatrice_to_phase(indicatrice_left_left)));
                      bool phase_invalide_left = (left_monophasique && (phase != IJK_Interfaces::convert_indicatrice_to_phase(indicatrice_left)));
                      bool phase_invalide_centre = (centre_monophasique && (phase != IJK_Interfaces::convert_indicatrice_to_phase(indicatrice_centre)));
                      bool phase_invalide_right = (right_monophasique && (phase != IJK_Interfaces::convert_indicatrice_to_phase(indicatrice_right)));

                      const double delta_xyz = _DIR_==DIRECTION::Z ? (channel_data_.get_delta_z()[k-1] + channel_data_.get_delta_z()[k]) * 0.5 : channel_data_.get_delta((int)_DIR_);

                      const DoubleTabFT_cut_cell_vector3& indicatrice_surfacique = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face();
                      double indicatrice_surface = (n_centre < 0) ? ((phase == 0) ? 1 - indicatrice_centre : indicatrice_centre) : ((phase == 0) ? 1 - indicatrice_surfacique(n_centre,dir) : indicatrice_surfacique(n_centre,dir));

                      const double surface = channel_data_.get_surface(k, 1, (int)_DIR_) * indicatrice_surface;

                      double velocity = (n_centre < 0) ? velocity_dir.pure_(i,j,k) : diph_velocity(n_centre);

                      int phase_left_left = (n_left_left < 0) ? IJK_Interfaces::convert_indicatrice_to_phase(indicatrice_left_left) : phase;
                      int phase_left = (n_left < 0) ? IJK_Interfaces::convert_indicatrice_to_phase(indicatrice_left) : phase;
                      int phase_right = (n_right < 0) ? IJK_Interfaces::convert_indicatrice_to_phase(indicatrice_right) : phase;
                      assert((phase_left == phase) || (indicatrice_surface == 0.));

                      double input_left_left = (n_left_left < 0) ? input_cut_field.pure_(i-2*dir_i,j-2*dir_j,k-2*dir_k) : diph_input(n_left_left);
                      double input_left = (n_left < 0) ? input_cut_field.pure_(i-dir_i,j-dir_j,k-dir_k) : diph_input(n_left);
                      double input_centre = (n_centre < 0) ? input_cut_field.pure_(i,j,k) : diph_input(n_centre);
                      double input_right = (n_right < 0) ? input_cut_field.pure_(i+dir_i,j+dir_j,k+dir_k) : diph_input(n_right);

                      double flux_2 = 0.;
                      if        ((cut_cell_conv_scheme_.scheme == CUT_CELL_SCHEMA_CONVECTION::CENTRE2)
                                 || (cut_cell_conv_scheme_.scheme == CUT_CELL_SCHEMA_CONVECTION::QUICK_OU_CENTRE2_STENCIL)
                                 || (cut_cell_conv_scheme_.scheme == CUT_CELL_SCHEMA_CONVECTION::QUICK_OU_CENTRE2_PERPENDICULAR_DISTANCE))
                        {
                          double input_milieu = (input_left + input_centre) * 0.5;
                          flux_2 = OpConvQuickIJKScalar_cut_cell_double::compute_flux_local_<_DIR_>(surface, velocity, input_milieu);
                        }
                      else if ((cut_cell_conv_scheme_.scheme == CUT_CELL_SCHEMA_CONVECTION::LINEAIRE2)
                               || (cut_cell_conv_scheme_.scheme == CUT_CELL_SCHEMA_CONVECTION::QUICK_OU_LINEAIRE2_STENCIL)
                               || (cut_cell_conv_scheme_.scheme == CUT_CELL_SCHEMA_CONVECTION::QUICK_OU_LINEAIRE2_PERPENDICULAR_DISTANCE))
                        {
                          double bar_dir_left = cut_cell_disc.get_interfaces().get_barycentre(true, dir, phase, i-dir_i,j-dir_j,k-dir_k);
                          assert((n_left >= 0) || (bar_dir_left == .5));

                          double bar_dir_centre = cut_cell_disc.get_interfaces().get_barycentre(true, dir, phase, i,j,k);
                          assert((n_centre >= 0) || (bar_dir_centre == .5));

                          // Note : suppose un maillage uniforme, splitting_.is_uniform(_DIR_)
                          double input_milieu = (input_left + input_centre) * 0.5;
                          double input_lineaire = (1 - bar_dir_left + bar_dir_centre) == 0. ? input_milieu : input_left + (1 - bar_dir_left)/(1 - bar_dir_left + bar_dir_centre)*(input_centre - input_left);
                          assert(std::abs(input_lineaire - input_left) <= std::abs(input_centre - input_left));
                          assert((input_lineaire - input_left)*(input_centre - input_left) >= 0);
                          flux_2 = OpConvQuickIJKScalar_cut_cell_double::compute_flux_local_<_DIR_>(surface, velocity, input_lineaire);
                        }
                      else if ((cut_cell_conv_scheme_.scheme == CUT_CELL_SCHEMA_CONVECTION::AMONT)
                               || (cut_cell_conv_scheme_.scheme == CUT_CELL_SCHEMA_CONVECTION::QUICK_OU_AMONT_STENCIL)
                               || (cut_cell_conv_scheme_.scheme == CUT_CELL_SCHEMA_CONVECTION::QUICK_OU_AMONT_PERPENDICULAR_DISTANCE))
                        {
                          double input_amont = velocity < 0. ? input_centre : input_left;
                          flux_2 = OpConvQuickIJKScalar_cut_cell_double::compute_flux_local_<_DIR_>(surface, velocity, input_amont);
                        }
                      else
                        {
                          Cerr << "OpConvQuickIJKScalar_cut_cell.tpp: cut_cell_conv_scheme_ inconnu pour flux_2." << finl;
                          Process::exit();
                        }


                      double flux_1l = OpConvQuickIJKScalar_cut_cell_double::compute_flux_local_<_DIR_>(surface, velocity, input_left);
                      double flux_1c = OpConvQuickIJKScalar_cut_cell_double::compute_flux_local_<_DIR_>(surface, velocity, input_centre);

                      double flux_4 = OpConvQuickIJKScalar_cut_cell_double::compute_flux_local_<_DIR_>(k, delta_xyz, surface, velocity, input_left_left, input_left, input_centre, input_right);
                      if      ((cut_cell_conv_scheme_.scheme == CUT_CELL_SCHEMA_CONVECTION::QUICK_OU_CENTRE2_STENCIL)
                               || (cut_cell_conv_scheme_.scheme == CUT_CELL_SCHEMA_CONVECTION::QUICK_OU_LINEAIRE2_STENCIL)
                               || (cut_cell_conv_scheme_.scheme == CUT_CELL_SCHEMA_CONVECTION::QUICK_OU_AMONT_STENCIL)
                               || (left_left_monophasique && left_monophasique && centre_monophasique && right_monophasique))
                        {
                          // Ne fait rien : on garde le flux_4 quick
                        }
                      else if ((cut_cell_conv_scheme_.scheme == CUT_CELL_SCHEMA_CONVECTION::QUICK_OU_CENTRE2_PERPENDICULAR_DISTANCE)
                               || (cut_cell_conv_scheme_.scheme == CUT_CELL_SCHEMA_CONVECTION::QUICK_OU_LINEAIRE2_PERPENDICULAR_DISTANCE)
                               || (cut_cell_conv_scheme_.scheme == CUT_CELL_SCHEMA_CONVECTION::QUICK_OU_AMONT_PERPENDICULAR_DISTANCE))
                        {
                          double bar_dir_perp1_left_left = cut_cell_disc.get_interfaces().get_barycentre(true, (dir+1)%3, phase, i-2*dir_i,j-2*dir_j,k-2*dir_k);
                          double bar_dir_perp2_left_left = cut_cell_disc.get_interfaces().get_barycentre(true, (dir+2)%3, phase, i-2*dir_i,j-2*dir_j,k-2*dir_k);
                          assert((n_left_left >= 0) || (bar_dir_perp1_left_left == .5));
                          assert((n_left_left >= 0) || (bar_dir_perp2_left_left == .5));

                          double bar_dir_perp1_left = cut_cell_disc.get_interfaces().get_barycentre(true, (dir+1)%3, phase, i-dir_i,j-dir_j,k-dir_k);
                          double bar_dir_perp2_left = cut_cell_disc.get_interfaces().get_barycentre(true, (dir+2)%3, phase, i-dir_i,j-dir_j,k-dir_k);
                          assert((n_left >= 0) || (bar_dir_perp1_left == .5));
                          assert((n_left >= 0) || (bar_dir_perp2_left == .5));

                          double bar_dir_perp1_centre = cut_cell_disc.get_interfaces().get_barycentre(true, (dir+1)%3, phase, i,j,k);
                          double bar_dir_perp2_centre = cut_cell_disc.get_interfaces().get_barycentre(true, (dir+2)%3, phase, i,j,k);
                          assert((n_centre >= 0) || (bar_dir_perp1_centre == .5));
                          assert((n_centre >= 0) || (bar_dir_perp2_centre == .5));

                          double bar_dir_perp1_right = cut_cell_disc.get_interfaces().get_barycentre(true, (dir+1)%3, phase, i+dir_i,j+dir_j,k+dir_k);
                          double bar_dir_perp2_right = cut_cell_disc.get_interfaces().get_barycentre(true, (dir+2)%3, phase, i+dir_i,j+dir_j,k+dir_k);
                          assert((n_right >= 0) || (bar_dir_perp1_right == .5));
                          assert((n_right >= 0) || (bar_dir_perp2_right == .5));

                          bool small_perpendicular_distance_left_left = ((std::abs(bar_dir_perp1_left_left - .5) < 0.1) && (std::abs(bar_dir_perp2_left_left - .5) < 0.1));
                          bool small_perpendicular_distance_left      = ((std::abs(bar_dir_perp1_left      - .5) < 0.1) && (std::abs(bar_dir_perp2_left      - .5) < 0.1));
                          bool small_perpendicular_distance_centre    = ((std::abs(bar_dir_perp1_centre    - .5) < 0.1) && (std::abs(bar_dir_perp2_centre    - .5) < 0.1));
                          bool small_perpendicular_distance_right     = ((std::abs(bar_dir_perp1_right     - .5) < 0.1) && (std::abs(bar_dir_perp2_right     - .5) < 0.1));

                          if (small_perpendicular_distance_left_left && small_perpendicular_distance_left && small_perpendicular_distance_centre && small_perpendicular_distance_right)
                            {
                              // Ne fait rien : on garde le flux_4 quick
                            }
                          else
                            {
                              flux_4 = flux_2;
                            }
                        }
                      else if ((cut_cell_conv_scheme_.scheme == CUT_CELL_SCHEMA_CONVECTION::CENTRE2)
                               || (cut_cell_conv_scheme_.scheme == CUT_CELL_SCHEMA_CONVECTION::LINEAIRE2)
                               || (cut_cell_conv_scheme_.scheme == CUT_CELL_SCHEMA_CONVECTION::AMONT))
                        {
                          flux_4 = flux_2;
                        }
                      else
                        {
                          Cerr << "OpConvQuickIJKScalar_cut_cell.tpp: cut_cell_conv_scheme_ inconnu pour flux_4." << finl;
                          Process::exit();
                        }

                      int phase_mourrante_left_left = cut_cell_disc.get_interfaces().phase_mourrante(phase, indicatrice_left_left, next_indicatrice_left_left);
                      int phase_mourrante_left = cut_cell_disc.get_interfaces().phase_mourrante(phase, indicatrice_left, next_indicatrice_left);
                      int phase_mourrante_centre = cut_cell_disc.get_interfaces().phase_mourrante(phase, indicatrice_centre, next_indicatrice_centre);
                      int phase_mourrante_right = cut_cell_disc.get_interfaces().phase_mourrante(phase, indicatrice_right, next_indicatrice_right);
                      int phase_naissante_left_left = cut_cell_disc.get_interfaces().phase_naissante(phase, indicatrice_left_left, next_indicatrice_left_left);
                      int phase_naissante_left = cut_cell_disc.get_interfaces().phase_naissante(phase, indicatrice_left, next_indicatrice_left);
                      int phase_naissante_centre = cut_cell_disc.get_interfaces().phase_naissante(phase, indicatrice_centre, next_indicatrice_centre);
                      int phase_naissante_right = cut_cell_disc.get_interfaces().phase_naissante(phase, indicatrice_right, next_indicatrice_right);
                      //int petit_left_left = cut_cell_disc.get_interfaces().next_below_small_threshold_for_phase(phase, indicatrice_left_left, next_indicatrice_left_left);
                      int petit_left = cut_cell_disc.get_interfaces().next_below_small_threshold_for_phase(phase, indicatrice_left, next_indicatrice_left);
                      int petit_centre = cut_cell_disc.get_interfaces().next_below_small_threshold_for_phase(phase, indicatrice_centre, next_indicatrice_centre);
                      //int petit_right = cut_cell_disc.get_interfaces().next_below_small_threshold_for_phase(phase, indicatrice_right, next_indicatrice_right);

                      double flux_value;
                      if (phase_naissante_centre)
                        {
                          assert(std::abs(input_centre) < 1e-16);
                          //assert(input_centre == 0);
                          flux_value = flux_1l;
                        }
                      else if (phase_naissante_left)
                        {
                          assert(std::abs(input_left) < 1e-16);
                          //assert(input_left == 0);
                          flux_value = flux_1c;
                        }
                      else if (ignore_small_cells_ && (petit_centre || petit_left || phase_mourrante_centre || phase_mourrante_left))
                        {
                          flux_value = 0.;
                        }
                      else if (phase_mourrante_left_left || phase_mourrante_right || phase_naissante_left_left || phase_naissante_right) //|| petit_left_left || petit_right)
                        {
                          assert(phase_invalide_left || (input_left != 0.));
                          assert(phase_invalide_centre || (input_centre != 0.));
                          flux_value = flux_2;
                        }
                      else if ((phase_left_left == phase && (!phase_invalide_left_left)) && (phase_right == phase && (!phase_invalide_right)))
                        {
                          assert(input_left_left != 0);
                          assert(phase_invalide_left || (input_left != 0.));
                          assert(phase_invalide_centre || (input_centre != 0.));
                          assert(input_right != 0);
                          flux_value = flux_4;
                        }
                      else
                        {
                          assert(phase_invalide_left || (input_left != 0.));
                          assert(phase_invalide_centre || (input_centre != 0.));
                          flux_value = flux_2;
                        }

                      if (n_centre < 0) // Si la cellule est pure
                        {
                          int index_ijk_per = 0;
                          while (index_ijk_per >= 0)
                            {
                              Int3 ijk = cut_cell_disc.ijk_per_of_index(i, j, k, index_ijk_per);
                              index_ijk_per = cut_cell_disc.next_index_ijk_per(i, j, k, index_ijk_per, 0, 1);

                              (*flux)(ijk[0],ijk[1],0) = flux_value;
                            }
                        }
                      else
                        {
                          diph_flux(n_centre) = flux_value;
                          if ((n_left < 0) && (phase == phase_left && (!phase_invalide_centre) && (!phase_invalide_left)))
                            {
                              int index_ijk_per = 0;
                              while (index_ijk_per >= 0)
                                {
                                  Int3 ijk = cut_cell_disc.ijk_per_of_index(i, j, k, index_ijk_per);
                                  index_ijk_per = cut_cell_disc.next_index_ijk_per(i, j, k, index_ijk_per, 0, 1);

                                  (*flux)(ijk[0],ijk[1],0) = flux_value;
                                }
                            }
                        }

                      // Si on est monophasique de la phase non-existante, le flux devrait etre nul
                      if (phase_invalide_centre || phase_invalide_left)
                        {
                          assert(flux_value == 0.);
                        }
                    }
                }
            }
        }
    }
}

#endif
