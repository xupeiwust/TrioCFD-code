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

#ifndef Operateur_IJK_elem_diff_base_H
#define Operateur_IJK_elem_diff_base_H

#include <Operateur_IJK_base.h>
#include <Cut_cell_FT_Disc.h>

enum class BOUNDARY_FLUX { NOT_DETERMINED_BY_BOUNDARY=0, KMIN_WALL=1, KMAX_WALL=2 };

class Operateur_IJK_elem_diff_base_double : public Operateur_IJK_elem_base_double
{
  Declare_base_sans_constructeur(Operateur_IJK_elem_diff_base_double);
public:
  Operateur_IJK_elem_diff_base_double();
  virtual void initialize(const IJK_Splitting& splitting) override;
  virtual void set_indicatrice(const IJK_Field_double& indicatrice) { indicatrice_= &indicatrice; };
  virtual void set_corrige_flux(OWN_PTR(Corrige_flux_FT_base)& corrige_flux) { corrige_flux_ = &corrige_flux; };
  virtual void calculer(const IJK_Field_double& field,
                        IJK_Field_double& result,
                        const IJK_Field_local_double& boundary_flux_kmin,
                        const IJK_Field_local_double& boundary_flux_kmax);
  virtual void ajouter(const IJK_Field_double& field,
                       IJK_Field_double& result,
                       const IJK_Field_local_double& boundary_flux_kmin,
                       const IJK_Field_local_double& boundary_flux_kmax);

  inline void set_uniform_lambda(const double& uniform_lambda) { uniform_lambda_ = &uniform_lambda; };
  inline void set_uniform_lambda_liquid(const double& uniform_lambda_liquid) { uniform_lambda_liquid_ = &uniform_lambda_liquid; };
  inline void set_uniform_lambda_vapour(const double& uniform_lambda_vapour) { uniform_lambda_vapour_ = &uniform_lambda_vapour; };

  inline void set_lambda(const IJK_Field_local_double& lambda) { lambda_ = &lambda; };

  inline void set_coeff_x_y_z(const IJK_Field_local_double& coeff_field_x,
                              const IJK_Field_local_double& coeff_field_y,
                              const IJK_Field_local_double& coeff_field_z)
  {
    coeff_field_x_ = &coeff_field_x;
    coeff_field_y_ = &coeff_field_y;
    coeff_field_z_ = &coeff_field_z;
  }

  inline double get_uniform_lambda()
  {
    if (uniform_lambda_ == nullptr)
      return 1.;
    else
      return *uniform_lambda_;
  }

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
  inline double compute_flux_local_x(int i, int j, int k) override
  {
    return compute_flux_local_<DIRECTION::X>(i,j,k);
  }
  inline double compute_flux_local_y(int i, int j, int k) override
  {
    return compute_flux_local_<DIRECTION::Y>(i,j,k);
  }
  inline double compute_flux_local_z(int i, int j, int k) override
  {
    return compute_flux_local_<DIRECTION::Z>(i,j,k);
  }

protected:
  template <DIRECTION _DIR_>
  inline void compute_flux_(IJK_Field_local_double& resu, const int k_layer);

  template <DIRECTION _DIR_>
  inline double compute_flux_local_(int i, int j, int k);

  template <DIRECTION _DIR_>
  inline double compute_flux_local_(double d0, double d1, double surface, double input_left, double input_centre, double lambda_left, double lambda_centre, double structural_model);

  template <DIRECTION _DIR_>
  inline BOUNDARY_FLUX flux_determined_by_boundary_condition_(int k);

  template <DIRECTION _DIR_>
  inline double compute_flux_local_boundary_condition_(BOUNDARY_FLUX type_boundary_flux, int i, int j);

  template <DIRECTION _DIR_>
  inline Vecteur3 compute_surface_d0_d1_(int k);

  const IJK_Field_local_double& get_model(DIRECTION _DIR_);

  Operateur_IJK_data_channel channel_data_;
  bool perio_k_;

  // Pointers to input data (set by calculer, used by compute_flux_...)
  const IJK_Field_double *input_field_;

  const IJK_Field_local_double *lambda_;
  const double *uniform_lambda_;
  const double *uniform_lambda_liquid_;
  const double *uniform_lambda_vapour_;

  const IJK_Field_local_double *coeff_field_x_;
  const IJK_Field_local_double *coeff_field_y_;
  const IJK_Field_local_double *coeff_field_z_;

  const IJK_Field_local_double *boundary_flux_kmin_;
  const IJK_Field_local_double *boundary_flux_kmax_;

  OWN_PTR(Corrige_flux_FT_base) *corrige_flux_;
  const IJK_Field_local_double *indicatrice_;

  bool is_anisotropic_;
  bool is_vectorial_;
  bool is_structural_;
  bool is_uniform_;
  bool is_corrected_;
  bool is_hess_;
  bool is_flux_;

};

class OpDiffUniformIJKScalar_double : public Operateur_IJK_elem_diff_base_double
{
  Declare_instanciable_sans_constructeur(OpDiffUniformIJKScalar_double);
public:
  OpDiffUniformIJKScalar_double() : Operateur_IJK_elem_diff_base_double() { is_uniform_ = true; }
};

class OpDiffUniformIJKScalarCorrection_double : public Operateur_IJK_elem_diff_base_double
{
  Declare_instanciable_sans_constructeur(OpDiffUniformIJKScalarCorrection_double);
public:
  OpDiffUniformIJKScalarCorrection_double() : Operateur_IJK_elem_diff_base_double() { is_uniform_ = true, is_corrected_ = true; }
private:
  void correct_flux(IJK_Field_local_double *const flux, const int k_layer, const int dir) override;
  void correct_flux_spherical(Simd_double& a, Simd_double& b, const int& i, const int& j, int k_layer, int dir) override;
};

class OpDiffIJKScalar_cut_cell_double : public Operateur_IJK_elem_diff_base_double
{
  Declare_instanciable_sans_constructeur(OpDiffIJKScalar_cut_cell_double);
public:
  OpDiffIJKScalar_cut_cell_double() : Operateur_IJK_elem_diff_base_double() {}
  void initialise_cut_cell(bool ignore_small_cells,
                           FixedVector<Cut_cell_double, 3>& cut_cell_flux,
                           IJK_Field_int& treatment_count,
                           int& new_treatment)
  {
    ignore_small_cells_ = ignore_small_cells;
    cut_cell_flux_ = &cut_cell_flux;
    treatment_count_ = &treatment_count;
    new_treatment_ = &new_treatment;
  }

  void set_runge_kutta(int rk_step, double dt_tot, IJK_Field_vector3_double& current_fluxes, IJK_Field_vector3_double& RK3_F_fluxes, IJK_Field_int& cellule_rk_restreint_diff_main_l, IJK_Field_int& cellule_rk_restreint_diff_main_v)
  {
    if (rk_step == -1)
      {
        runge_kutta_flux_correction_ = false;
        rk_step_ = -1;
        dt_tot_ = 0;
        current_fluxes_ = nullptr;
        RK3_F_fluxes_ = nullptr;
        cellule_rk_restreint_diff_main_l_ = nullptr;
        cellule_rk_restreint_diff_main_v_ = nullptr;
      }
    else
      {
        runge_kutta_flux_correction_ = true;
        rk_step_ = rk_step;
        dt_tot_ = dt_tot;
        current_fluxes_ = &current_fluxes;
        RK3_F_fluxes_ = &RK3_F_fluxes;
        cellule_rk_restreint_diff_main_l_ = &cellule_rk_restreint_diff_main_l;
        cellule_rk_restreint_diff_main_v_ = &cellule_rk_restreint_diff_main_v;
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

private:
  void correct_flux(IJK_Field_local_double *const flux, const int k_layer, const int dir) override;

  template <DIRECTION _DIR_>
  inline void correct_flux_(IJK_Field_local_double *const flux, const int k_layer);


  bool ignore_small_cells_;

  FixedVector<Cut_cell_double, 3> *cut_cell_flux_;
  IJK_Field_int *cellule_rk_restreint_diff_main_l_;
  IJK_Field_int *cellule_rk_restreint_diff_main_v_;

  IJK_Field_int *treatment_count_;
  int *new_treatment_;

};

class OpDiffIJKScalar_double : public Operateur_IJK_elem_diff_base_double
{
  Declare_instanciable_sans_constructeur(OpDiffIJKScalar_double);
public:
  OpDiffIJKScalar_double() : Operateur_IJK_elem_diff_base_double() {}
};

class OpDiffAnisotropicIJKScalar_double : public Operateur_IJK_elem_diff_base_double
{
  Declare_instanciable_sans_constructeur(OpDiffAnisotropicIJKScalar_double);
public:
  OpDiffAnisotropicIJKScalar_double() : Operateur_IJK_elem_diff_base_double() { is_anisotropic_ = true; }
};

class OpDiffVectorialIJKScalar_double : public Operateur_IJK_elem_diff_base_double
{
  Declare_instanciable_sans_constructeur(OpDiffVectorialIJKScalar_double);
public:
  OpDiffVectorialIJKScalar_double() : Operateur_IJK_elem_diff_base_double() { is_vectorial_ = true; }
};

class OpDiffVectorialAnisotropicIJKScalar_double : public Operateur_IJK_elem_diff_base_double
{
  Declare_instanciable_sans_constructeur(OpDiffVectorialAnisotropicIJKScalar_double);
public:
  OpDiffVectorialAnisotropicIJKScalar_double() : Operateur_IJK_elem_diff_base_double() { is_vectorial_ = true, is_anisotropic_ = true; }
};

class OpDiffStructuralOnlyIJKScalar_double : public Operateur_IJK_elem_diff_base_double
{
  Declare_instanciable_sans_constructeur(OpDiffStructuralOnlyIJKScalar_double);
public:
  OpDiffStructuralOnlyIJKScalar_double() : Operateur_IJK_elem_diff_base_double() { is_structural_ = true; }
};

#include <Operateur_IJK_elem_diff_base.tpp>

#endif
