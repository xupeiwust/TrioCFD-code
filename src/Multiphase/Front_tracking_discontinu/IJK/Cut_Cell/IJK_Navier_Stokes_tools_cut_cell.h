/****************************************************************************
* Copyright (c) 2023, CEA
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

#ifndef IJK_Navier_Stokes_tools_cut_cell_included
#define IJK_Navier_Stokes_tools_cut_cell_included

#include <Cut_field.h>
#include <IJK_Interfaces.h>

void ijk_interpolate_cut_cell(bool next_time, int phase, const Cut_field_double& field, const DoubleTab& coordinates, ArrOfDouble& result);
void ijk_interpolate_cut_cell_skip_unknown_points(bool next_time, int phase, const Cut_field_double& field, const DoubleTab& coordinates, ArrOfDouble& result, const double value_for_bad_points);

void ijk_interpolate_cut_cell(bool next_time, int phase, const Cut_cell_FT_Disc& cut_cell_disc, const DoubleTab& coordinates, IntTabFT& signed_independent_index, DoubleTabFT& coefficient);
void ijk_interpolate_cut_cell_skip_unknown_points(bool next_time, int phase, const Cut_cell_FT_Disc& cut_cell_disc, const DoubleTab& coordinates, IntTabFT& signed_independent_index, DoubleTabFT& coefficient, const double value_for_bad_points);

double ijk_interpolate_cut_cell_using_interface(bool next_time, int phase, const IJK_Field_double& field_ft, const Cut_field_double& field, const ArrOfDouble& interfacial_temperature, const Vecteur3& coordinates, int tolerate_not_within_tetrahedron, int& status);

double ijk_interpolate_cut_cell_using_interface_skip_unknown_points(bool next_time, int phase, const IJK_Field_double& field_ft, const Cut_field_double& field, const ArrOfDouble& interfacial_temperature, const Vecteur3& coordinates, int tolerate_not_within_tetrahedron, const double value_for_bad_points, int& status);

void ijk_interpolate(bool next_time, int phase, const Cut_cell_FT_Disc& cut_cell_disc, const Vecteur3& coordinates, IntTabFT& signed_independent_index, DoubleTabFT& coefficient);
void ijk_interpolate_skip_unknown_points(bool next_time, int phase, const Cut_cell_FT_Disc& cut_cell_disc, const Vecteur3& coordinates, IntTabFT& signed_independent_index, DoubleTabFT& coefficient, const double value_for_bad_points);

void cut_cell_switch_field_time(Cut_field_double& v);

void euler_explicit_update_cut_cell_transport(double timestep, const Cut_field_double& dv, Cut_field_double& v);
void runge_kutta3_update_cut_cell_transport(const Cut_field_double& dv, Cut_field_double& F, Cut_field_double& v, const int step, double dt_tot, const Cut_field_int& cellule_rk_restreint);

void euler_explicit_update_cut_cell_notransport(double timestep, bool next_time, const Cut_field_double& dv, Cut_field_double& v);
void runge_kutta3_update_cut_cell_notransport(bool next_time, const Cut_field_double& dv, Cut_field_double& F, Cut_field_double& v, const int step, double dt_tot, const Cut_field_int& cellule_rk_restreint);

void runge_kutta3_update_surfacic_fluxes(Cut_field_double& dv, Cut_field_double& F, const int step, const int k_layer, const int dir, double dt_tot, const Cut_field_int& cellule_rk_restreint);

void add_flux_times_vol_over_dt_surface(double fractional_timestep, const Cut_field_vector3_double& cut_field_current_fluxes, Cut_field_vector3_double& cut_field_RK3_F_fluxes);
void set_rk_restreint(int rk_step, int rk_restriction_leniency, const Cut_cell_FT_Disc& cut_cell_disc, Cut_field_int& cellule_rk_restreint);

#endif /* IJK_Navier_Stokes_tools_cut_cell_included */
