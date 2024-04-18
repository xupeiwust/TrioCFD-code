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
  virtual void set_corrige_flux(Corrige_flux_FT& corrige_flux) { corrige_flux_ = &corrige_flux; };
  virtual void calculer(const IJK_Field_double& field,
                        IJK_Field_double& result,
                        const IJK_Field_local_double& boundary_flux_kmin,
                        const IJK_Field_local_double& boundary_flux_kmax);
  virtual void ajouter(const IJK_Field_double& field,
                       IJK_Field_double& result,
                       const IJK_Field_local_double& boundary_flux_kmin,
                       const IJK_Field_local_double& boundary_flux_kmax);
  virtual void calculer_cut_cell(bool ignore_small_cells,
                                 const Cut_field_scalar& field,
                                 Cut_cell_vector& cut_cell_flux,
                                 Cut_field_scalar& result,
                                 const IJK_Field_local_double& boundary_flux_kmin,
                                 const IJK_Field_local_double& boundary_flux_kmax);
  virtual void ajouter_cut_cell(bool ignore_small_cells,
                                const Cut_field_scalar& field,
                                Cut_cell_vector& cut_cell_flux,
                                Cut_field_scalar& result,
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
  DoubleTabFT_cut_cell* get_diph_flux(int phase) override
  {
    return (phase == 0) ? &cut_cell_flux_->diph_v_ : &cut_cell_flux_->diph_l_;
  }
  inline void compute_cut_cell_divergence(int phase, const DoubleTabFT_cut_cell& diph_flux,
                                          const IJK_Field_local_double& flux_x,
                                          const IJK_Field_local_double& flux_y,
                                          const IJK_Field_local_double& flux_zmin,
                                          const IJK_Field_local_double& flux_zmax,
                                          DoubleTabFT_cut_cell& diph_resu, int k_layer, bool add) override
  {
    const Cut_cell_FT_Disc& cut_cell_disc = cut_cell_flux_->get_cut_cell_disc();

    for (int index = cut_cell_disc.get_k_value_index(k_layer); index < cut_cell_disc.get_k_value_index(k_layer+1); index++)
      {
        int n = cut_cell_disc.get_n_from_k_index(index);
        Int3 ijk = cut_cell_disc.get_ijk(n);

        int i = ijk[0];
        int j = ijk[1];
        int k = ijk[2];

        if (!cut_cell_disc.within_ghost(i, j, k, 0, 0))
          continue;

        BOUNDARY_FLUX type_boundary_flux = flux_determined_by_boundary_condition_<DIRECTION::Z>(k);
        if (type_boundary_flux != BOUNDARY_FLUX::NOT_DETERMINED_BY_BOUNDARY)
          {
            Cerr << "Le cas d'une cellule diphasique avec flux condition aux limites n'est pas traite" << finl;
            Process::exit();
          }
        else
          {
            int n_ip1 = cut_cell_disc.get_n(i+1,j,k);
            int n_jp1 = cut_cell_disc.get_n(i,j+1,k);
            int n_kp1 = cut_cell_disc.get_n(i,j,k+1); // ???k

            double indicatrice_ip1 = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().I(i+1,j,k) : cut_cell_disc.get_interfaces().I(i+1,j,k);
            double indicatrice_jp1 = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().I(i,j+1,k) : cut_cell_disc.get_interfaces().I(i,j+1,k);
            double indicatrice_kp1 = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().I(i,j,k+1) : cut_cell_disc.get_interfaces().I(i,j,k+1);

            double fx_centre = diph_flux(n,0);
            double fy_centre = diph_flux(n,1);
            double fz_centre = diph_flux(n,2);

            double fx_right  = (n_ip1 < 0) ? indicatrice_ip1*flux_x(i+1,j,0)  : diph_flux(n_ip1,0);
            double fy_right  = (n_jp1 < 0) ? indicatrice_jp1*flux_y(i,j+1,0)  : diph_flux(n_jp1,1);
            double fz_right  = (n_kp1 < 0) ? indicatrice_kp1*flux_zmax(i,j,0) : diph_flux(n_kp1,2);

            double r = 0;
            r += fx_centre - fx_right;
            r += fy_centre - fy_right;
            r += fz_centre - fz_right;

            if(add)
              {
                r += diph_resu(n);
              }
            diph_resu(n) = r;
          }
      }
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

  template <DIRECTION _DIR_>
  inline void correct_flux_(IJK_Field_local_double *const flux, const int k_layer);

  const IJK_Field_local_double& get_model(DIRECTION _DIR_);

  Operateur_IJK_data_channel channel_data_;
  bool perio_k_;

  // Pointers to input data (set by calculer, used by compute_flux_...)
  const IJK_Field_local_double *input_field_;

  bool *ignore_small_cells_;

  const Cut_field_scalar *input_cut_field_;
  Cut_cell_vector *cut_cell_flux_;

  const IJK_Field_local_double *lambda_;
  const double *uniform_lambda_;
  const double *uniform_lambda_liquid_;
  const double *uniform_lambda_vapour_;

  const IJK_Field_local_double *coeff_field_x_;
  const IJK_Field_local_double *coeff_field_y_;
  const IJK_Field_local_double *coeff_field_z_;

  const IJK_Field_local_double *boundary_flux_kmin_;
  const IJK_Field_local_double *boundary_flux_kmax_;

  Corrige_flux_FT *corrige_flux_;
  const IJK_Field_local_double *indicatrice_;

  bool is_anisotropic_;
  bool is_vectorial_;
  bool is_structural_;
  bool is_uniform_;
  bool is_corrected_;
  bool is_hess_;

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
  void calculer(const IJK_Field_double& field,
                IJK_Field_double& result,
                const IJK_Field_local_double& boundary_flux_kmin,
                const IJK_Field_local_double& boundary_flux_kmax) override
  {
    Cerr << "The cut cell operators demand the use of calculer_cut_cell instead of calculer." << finl;
    Process::exit();
  }
  void ajouter(const IJK_Field_double& field,
               IJK_Field_double& result,
               const IJK_Field_local_double& boundary_flux_kmin,
               const IJK_Field_local_double& boundary_flux_kmax) override
  {
    Cerr << "The cut cell operators demand the use of ajouter_cut_cell instead of ajouter." << finl;
    Process::exit();
  }
private:
  void correct_flux(IJK_Field_local_double *const flux, const int k_layer, const int dir) override;
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
