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

#ifndef OpConvQuickScalarIJK_cut_cell_included
#define OpConvQuickScalarIJK_cut_cell_included

#include <Operateur_IJK_elem_conv_base.h>

class OpConvQuickIJKScalar_cut_cell_double : public Operateur_IJK_elem_conv_base_double
{
  Declare_instanciable_sans_constructeur(OpConvQuickIJKScalar_cut_cell_double);

public:
  OpConvQuickIJKScalar_cut_cell_double() : Operateur_IJK_elem_conv_base_double() { };

  void initialise_cut_cell(Cut_cell_conv_scheme cut_cell_conv_scheme,
                           const FixedVector<IJK_Field_vector3_double, 2>& temperature_face,
                           bool ignore_small_cells,
                           FixedVector<Cut_cell_double, 3>& cut_cell_flux,
                           IJK_Field_int& treatment_count,
                           int& new_treatment)
  {
    cut_cell_conv_scheme_ = cut_cell_conv_scheme;
    temperature_face_ = &temperature_face;
    ignore_small_cells_ = ignore_small_cells;
    cut_cell_flux_ = &cut_cell_flux;
    treatment_count_ = &treatment_count;
    new_treatment_ = &new_treatment;
  }

  void set_runge_kutta(int rk_step, double dt_tot, IJK_Field_vector3_double& current_fluxes, IJK_Field_vector3_double& RK3_F_fluxes, IJK_Field_int& cellule_rk_restreint_conv_main_l, IJK_Field_int& cellule_rk_restreint_conv_main_v)
  {
    if (rk_step == -1)
      {
        runge_kutta_flux_correction_ = false;
        rk_step_ = -1;
        dt_tot_ = 0;
        current_fluxes_ = nullptr;
        RK3_F_fluxes_ = nullptr;
        cellule_rk_restreint_conv_main_l_ = nullptr;
        cellule_rk_restreint_conv_main_v_ = nullptr;
      }
    else
      {
        runge_kutta_flux_correction_ = true;
        rk_step_ = rk_step;
        dt_tot_ = dt_tot;
        current_fluxes_ = &current_fluxes;
        RK3_F_fluxes_ = &RK3_F_fluxes;
        cellule_rk_restreint_conv_main_l_ = &cellule_rk_restreint_conv_main_l;
        cellule_rk_restreint_conv_main_v_ = &cellule_rk_restreint_conv_main_v;
      }
  };

  const Cut_cell_FT_Disc* get_cut_cell_disc()
  {
    assert(&(*cut_cell_flux_)[0].get_cut_cell_disc() == &(*cut_cell_flux_)[1].get_cut_cell_disc());
    assert(&(*cut_cell_flux_)[0].get_cut_cell_disc() == &(*cut_cell_flux_)[2].get_cut_cell_disc());
    return &(*cut_cell_flux_)[0].get_cut_cell_disc();
  }

  FixedVector<Cut_cell_double, 3>* get_cut_cell_flux()
  {
    return cut_cell_flux_;
  }

  void compute_cut_cell_divergence(const FixedVector<Cut_cell_double, 3>& cut_cell_flux,
                                   const IJK_Field_local_double& flux_x,
                                   const IJK_Field_local_double& flux_y,
                                   const IJK_Field_local_double& flux_zmin,
                                   const IJK_Field_local_double& flux_zmax,
                                   IJK_Field_double& resu, int k_layer, bool add);

  void Operator_IJK_div(const IJK_Field_local_double& flux_x, const IJK_Field_local_double& flux_y,
                        const IJK_Field_local_double& flux_zmin, const IJK_Field_local_double& flux_zmax,
                        IJK_Field_double& resu, int k_layer, bool add) override
  {
    Operateur_IJK_elem_base_double::Operator_IJK_div(flux_x, flux_y, flux_zmin, flux_zmax, resu, k_layer, add);

    FixedVector<Cut_cell_double, 3>& cut_cell_flux = *get_cut_cell_flux();
    compute_cut_cell_divergence(cut_cell_flux, flux_x, flux_y, flux_zmin, flux_zmax, resu, k_layer, add);
  }

protected:

  inline void compute_flux_x(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::X>(resu,k_layer);
  }
  inline void compute_flux_y(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::Y>(resu,k_layer);
  }
  inline void compute_flux_z(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::Z>(resu,k_layer);
  }

private:
  template <DIRECTION _DIR_>
  void compute_flux_(IJK_Field_local_double& resu, const int k_layer);

  template <DIRECTION _DIR_>
  double compute_flux_local_(int i, int j, int k);

  template <DIRECTION _DIR_>
  double compute_flux_local_(int k_layer, double delta_xyz, double surface, double velocity, double input_left_left, double input_left, double input_centre, double input_right);

  template <DIRECTION _DIR_>
  double compute_flux_local_(double surface, double velocity, double input);

  template <DIRECTION _DIR_>
  bool flux_determined_by_wall_(int k);

  template <DIRECTION _DIR_>
  Vecteur3 compute_curv_fram_local_(int k_layer, double input_left, double input_centre, double input_right);

  void correct_flux(IJK_Field_local_double *const flux, const int k_layer, const int dir) override;

  template <DIRECTION _DIR_>
  inline void correct_flux_(IJK_Field_local_double *const flux, const int k_layer);

  bool ignore_small_cells_;
  Cut_cell_conv_scheme cut_cell_conv_scheme_;

  const FixedVector<IJK_Field_vector3_double, 2> *temperature_face_;
  FixedVector<Cut_cell_double, 3> *cut_cell_flux_;
  IJK_Field_int *cellule_rk_restreint_conv_main_l_;
  IJK_Field_int *cellule_rk_restreint_conv_main_v_;

  IJK_Field_int *treatment_count_;
  int *new_treatment_;

};

#include <IJK_Field_vector.h>
#include <OpConvQuickIJKScalar_cut_cell.tpp>

#endif /* OpConvQuickScalarIJK_cut_cell_included */
