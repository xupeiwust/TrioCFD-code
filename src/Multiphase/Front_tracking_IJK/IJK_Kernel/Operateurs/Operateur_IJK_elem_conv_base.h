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

#ifndef Operateur_IJK_elem_conv_base_included
#define Operateur_IJK_elem_conv_base_included

#include <IJK_Splitting.h>
#include <Operateur_IJK_base.h>
#include <Cut_cell_FT_Disc.h>
#include <Cut_cell_convection_auxiliaire.h>

class Corrige_flux_FT;

class Operateur_IJK_elem_conv_base_double : public Operateur_IJK_elem_base_double
{
  Declare_base(Operateur_IJK_elem_conv_base_double);
public:
  void initialize(const IJK_Splitting& splitting) override;
  virtual void set_indicatrice(const IJK_Field_double& indicatrice) { indicatrice_= &indicatrice; };
  virtual void set_corrige_flux(Corrige_flux_FT& corrige_flux) { corrige_flux_ = &corrige_flux; };

  virtual void calculer(const IJK_Field_double& field,
                        const IJK_Field_double& vx,
                        const IJK_Field_double& vy,
                        const IJK_Field_double& vz,
                        IJK_Field_double& result);

  virtual void ajouter(const IJK_Field_double& field,
                       const IJK_Field_double& vx,
                       const IJK_Field_double& vy,
                       const IJK_Field_double& vz,
                       IJK_Field_double& result);

  virtual void calculer_cut_cell(bool ignore_small_cells,
                                 CUT_CELL_CONV_SCHEME cut_cell_conv_scheme,
                                 const Cut_field_scalar& field,
                                 const Cut_field_vector& v,
                                 const FixedVector<FixedVector<IJK_Field_double, 3>, 2>& temperature_face,
                                 Cut_cell_vector& cut_cell_flux,
                                 Cut_field_scalar& result);

  virtual void ajouter_cut_cell(bool ignore_small_cells,
                                CUT_CELL_CONV_SCHEME cut_cell_conv_scheme,
                                const Cut_field_scalar& field,
                                const Cut_field_vector& v,
                                const FixedVector<FixedVector<IJK_Field_double, 3>, 2>& temperature_face,
                                Cut_cell_vector& cut_cell_flux,
                                Cut_field_scalar& result);

protected:

  void compute_curv_fram(DIRECTION _DIR_, int k_layer);
  void shift_curv_fram(IJK_Field_local_double& tmp_curv_fram);
  inline const IJK_Field_local_double& get_input_velocity(DIRECTION _DIR_)
  {
    switch(_DIR_)
      {
      case DIRECTION::X:
        return *input_velocity_x_;
      case DIRECTION::Y:
        return *input_velocity_y_;
      case DIRECTION::Z:
        return *input_velocity_z_;
      default:
        Cerr << "Error in OpConvDiscIJKQuickScalar::get_input_velocity: wrong direction..." << finl;
        Process::exit();
      }
    //for compilation only...
    return *input_velocity_x_;
  }

  Operateur_IJK_data_channel channel_data_;

  // Pointers to input data (set by calculer, used by compute_flux_...)
  const IJK_Field_local_double *input_field_;

  bool *ignore_small_cells_;
  CUT_CELL_CONV_SCHEME *cut_cell_conv_scheme_;

  const Cut_field_scalar *input_cut_field_;
  const FixedVector<FixedVector<IJK_Field_double, 3>, 2> *temperature_face_;
  Cut_cell_vector *cut_cell_flux_;
  const Cut_field_vector *cut_field_velocity_;

  const IJK_Field_local_double *input_velocity_x_;
  const IJK_Field_local_double *input_velocity_y_;
  const IJK_Field_local_double *input_velocity_z_;
  bool perio_k_;

  // Temporary array to store curvature and fram coefficients
  // for the current computed flux.
  // layer k=0 and k=1 are used for "curv", k=2 and k=3 are used for "fram".
  // layer k=0 and k=2 store the previous values computed in direction "k" (which is used 2 times)
  IJK_Field_local_double tmp_curv_fram_;
  int stored_curv_fram_layer_z_; // which (local) layer is currently stored in layer 0 of the tmp array ?

  Corrige_flux_FT *corrige_flux_;
  const IJK_Field_local_double *indicatrice_;

  bool is_corrected_;
  bool is_grad_;

private:

  void compute_curv_fram_loop_(DIRECTION _DIR_, int iter, double factor12, double factor01, const ConstIJK_double_ptr& input_field, IJK_double_ptr& curv_values, IJK_double_ptr& fram_values );

};

inline Simd_double operator/(const Simd_double& x, const Simd_double& y)
{
  return SimdDivideMed(x, y);
}

#endif
