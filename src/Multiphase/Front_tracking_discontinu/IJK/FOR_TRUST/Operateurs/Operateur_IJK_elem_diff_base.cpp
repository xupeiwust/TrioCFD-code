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

#include <Operateur_IJK_elem_diff_base.h>
#include <IJK_Navier_Stokes_tools_cut_cell.h>

Implemente_base_sans_constructeur(Operateur_IJK_elem_diff_base_double, "Operateur_IJK_elem_diff_base_double", Operateur_IJK_elem_base_double);

Operateur_IJK_elem_diff_base_double::Operateur_IJK_elem_diff_base_double()
{
  input_field_ = nullptr;
  uniform_lambda_ = nullptr;
  uniform_lambda_liquid_ = nullptr;
  uniform_lambda_vapour_ = nullptr;
  lambda_ = nullptr;

  boundary_flux_kmin_ = boundary_flux_kmax_ = nullptr;

  coeff_field_x_ = nullptr;
  coeff_field_y_ = nullptr;
  coeff_field_z_ = nullptr;

  is_anisotropic_ = false;
  is_vectorial_ = false;
  is_structural_ = false;
  is_uniform_ = false;
  is_corrected_ = false;

  perio_k_ = false;
  is_hess_ = false;
  is_flux_ = false;

  corrige_flux_ = nullptr;
  indicatrice_ = nullptr;
}

Sortie& Operateur_IJK_elem_diff_base_double::printOn(Sortie& os) const
{
  return os;
}

Entree& Operateur_IJK_elem_diff_base_double::readOn(Entree& is)
{
  return is;
}

const IJK_Field_local_double& Operateur_IJK_elem_diff_base_double::get_model(DIRECTION _DIR_)
{
  assert(is_vectorial_);
  switch(_DIR_)
    {
    case DIRECTION::X:
      return *coeff_field_x_;
    case DIRECTION::Y:
      return *coeff_field_y_;
    case DIRECTION::Z:
      return *coeff_field_z_;
    default:
      Cerr << "Error in Operateur_IJK_elem_diff_base_double::get_model: wrong direction..." << finl;
      Process::exit();
    }
  return *coeff_field_x_;
}

void Operateur_IJK_elem_diff_base_double::initialize(const Domaine_IJK& splitting)
{
  channel_data_.initialize(splitting);
  perio_k_= splitting.get_periodic_flag(DIRECTION_K);
}

void Operateur_IJK_elem_diff_base_double::calculer(const IJK_Field_double& field,
                                                   IJK_Field_double& result,
                                                   const IJK_Field_local_double& boundary_flux_kmin,
                                                   const IJK_Field_local_double& boundary_flux_kmax)
{
  // Cerr << "Uniform lambda: " << get_uniform_lambda() << finl;
  input_field_ = &field;
  boundary_flux_kmin_ = &boundary_flux_kmin;
  boundary_flux_kmax_ = &boundary_flux_kmax;
  compute_set(result);
  input_field_ = nullptr;
  lambda_ = nullptr; // TODO: Why reset to nullptr? we could attribute it once only at initialize and never change it later. What was the reason?
  coeff_field_x_ = nullptr;
  coeff_field_y_ = nullptr;
  coeff_field_z_ = nullptr;
  boundary_flux_kmin_ = boundary_flux_kmax_ = nullptr;
}

void Operateur_IJK_elem_diff_base_double::ajouter(const IJK_Field_double& field,
                                                  IJK_Field_double& result,
                                                  const IJK_Field_local_double& boundary_flux_kmin,
                                                  const IJK_Field_local_double& boundary_flux_kmax)
{
  input_field_ = &field;
  boundary_flux_kmin_ = &boundary_flux_kmin;
  boundary_flux_kmax_ = &boundary_flux_kmax;
  compute_add(result);
  input_field_ = nullptr;
  lambda_ = nullptr;
  coeff_field_x_ = nullptr;
  coeff_field_y_ = nullptr;
  coeff_field_z_ = nullptr;
  boundary_flux_kmin_ = boundary_flux_kmax_ = nullptr;
}

Implemente_instanciable_sans_constructeur(OpDiffUniformIJKScalar_double, "OpDiffUniformIJKScalar_double", Operateur_IJK_elem_diff_base_double);

Sortie& OpDiffUniformIJKScalar_double::printOn(Sortie& os) const
{
  //  Operateur_IJK_elem_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffUniformIJKScalar_double::readOn(Entree& is)
{
  //  Operateur_IJK_elem_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffUniformIJKScalarCorrection_double, "OpDiffUniformIJKScalarCorrection_double", Operateur_IJK_elem_diff_base_double);

Sortie& OpDiffUniformIJKScalarCorrection_double::printOn(Sortie& os) const
{
  //  Operateur_IJK_elem_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffUniformIJKScalarCorrection_double::readOn(Entree& is)
{
  //  Operateur_IJK_elem_diff_base_double::readOn(is);
  return is;
}

void OpDiffUniformIJKScalarCorrection_double::correct_flux(IJK_Field_local_double *const flux, const int k_layer, const int dir)
{
  corrige_flux_->valeur().corrige_flux_diff_faceIJ(flux, k_layer, dir);
}

void OpDiffUniformIJKScalarCorrection_double::correct_flux_spherical(Simd_double& a, Simd_double& b, const int& i, const int& j, int k_layer, int dir)
{
  corrige_flux_->valeur().correct_flux_spherical(a, b, i, j, k_layer, dir);
}

Implemente_instanciable_sans_constructeur(OpDiffIJKScalar_cut_cell_double, "OpDiffIJKScalar_cut_cell_double", Operateur_IJK_elem_diff_base_double);

Sortie& OpDiffIJKScalar_cut_cell_double::printOn(Sortie& os) const
{
  return os;
}

Entree& OpDiffIJKScalar_cut_cell_double::readOn(Entree& is)
{
  return is;
}

void OpDiffIJKScalar_cut_cell_double::correct_flux(IJK_Field_local_double *const flux, int k_layer, const int dir)
{
  if (dir == 0)
    {
      correct_flux_<DIRECTION::X>(flux, k_layer);
    }
  else if (dir == 1)
    {
      correct_flux_<DIRECTION::Y>(flux, k_layer);
    }
  else if (dir == 2)
    {
      correct_flux_<DIRECTION::Z>(flux, k_layer);
    }
  else
    {
      Cerr << "Unexpected value of dir in OpDiffIJKScalar_cut_cell_double::correct_flux" << finl;
      Process::exit();
    }

  // Fluxes are stored in the flux variable (IJK_Field_local_double) for pure cells
  // as well as within the cut_cell_flux_ variable for cut cells.
  // The two variables must agrees for the fluxes between pure and cut cells.
  assert((*cut_cell_flux_)[dir].verify_consistency_within_layer(dir, k_layer, *flux));



  if (runge_kutta_flux_correction_)
    {
      Cut_field_vector3_double& cut_field_current_fluxes = static_cast<Cut_field_vector3_double&>(*current_fluxes_);
      Cut_field_vector3_double& cut_field_RK3_F_fluxes = static_cast<Cut_field_vector3_double&>(*RK3_F_fluxes_);

      assert(&(*cut_cell_flux_)[0].get_cut_cell_disc() == &(*cut_cell_flux_)[1].get_cut_cell_disc());
      assert(&(*cut_cell_flux_)[0].get_cut_cell_disc() == &(*cut_cell_flux_)[2].get_cut_cell_disc());
      const Cut_cell_FT_Disc& cut_cell_disc = (*cut_cell_flux_)[0].get_cut_cell_disc();

      {
        int ni = (dir == 0) ? flux->ni() : flux->ni() - 1;
        int nj = (dir == 1) ? flux->nj() : flux->nj() - 1;
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                cut_field_current_fluxes[dir].pure_(i,j,k_layer) = (*flux)(i,j,0);
                int n = cut_cell_disc.get_n(i, j, k_layer);
                if (n >= 0)
                  {
                    cut_field_current_fluxes[dir].diph_l_(n) = (*cut_cell_flux_)[dir].diph_l_(n);
                    cut_field_current_fluxes[dir].diph_v_(n) = (*cut_cell_flux_)[dir].diph_v_(n);
                  }
              }
          }
      }

      {
        int ni = (dir == 0) ? flux->ni() : flux->ni() - 1;
        int nj = (dir == 1) ? flux->nj() : flux->nj() - 1;
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                int n = cut_cell_disc.get_n(i, j, k_layer);
                if (n >= 0)
                  {
                    const DoubleTabFT_cut_cell_vector3& indicatrice_surfacique = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face();
                    double indicatrice_surface = indicatrice_surfacique(n,dir);

                    cut_field_RK3_F_fluxes[dir].diph_l_(n) = cut_field_RK3_F_fluxes[dir].diph_l_(n)*indicatrice_surface;
                    cut_field_RK3_F_fluxes[dir].diph_v_(n) = cut_field_RK3_F_fluxes[dir].diph_v_(n)*(1 -indicatrice_surface);
                  }
              }
          }
      }

      runge_kutta3_update_surfacic_fluxes(cut_field_current_fluxes[dir], cut_field_RK3_F_fluxes[dir], rk_step_, k_layer, dir, dt_tot_, *cellule_rk_restreint_diff_main_);

      {
        int ni = (dir == 0) ? flux->ni() : flux->ni() - 1;
        int nj = (dir == 1) ? flux->nj() : flux->nj() - 1;
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                (*flux)(i,j,0) = cut_field_current_fluxes[dir].pure_(i,j,k_layer);
                int n = cut_cell_disc.get_n(i, j, k_layer);
                if (n >= 0)
                  {
                    (*cut_cell_flux_)[dir].diph_l_(n) = cut_field_current_fluxes[dir].diph_l_(n);
                    (*cut_cell_flux_)[dir].diph_v_(n) = cut_field_current_fluxes[dir].diph_v_(n);
                  }
              }
          }
      }

      {
        int ni = (dir == 0) ? flux->ni() : flux->ni() - 1;
        int nj = (dir == 1) ? flux->nj() : flux->nj() - 1;
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                int n = cut_cell_disc.get_n(i, j, k_layer);
                if (n >= 0)
                  {
                    const DoubleTabFT_cut_cell_vector3& indicatrice_surfacique = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face();
                    double indicatrice_surface = indicatrice_surfacique(n,dir);

                    cut_field_RK3_F_fluxes[dir].diph_l_(n) = (indicatrice_surface == 0.)       ? 0. : cut_field_RK3_F_fluxes[dir].diph_l_(n)/indicatrice_surface;
                    cut_field_RK3_F_fluxes[dir].diph_v_(n) = ((1 - indicatrice_surface) == 0.) ? 0. : cut_field_RK3_F_fluxes[dir].diph_v_(n)/(1 -indicatrice_surface);
                  }
              }
          }
      }
    }

  assert((*cut_cell_flux_)[dir].verify_consistency_within_layer(dir, k_layer, *flux));
}

void OpDiffIJKScalar_cut_cell_double::compute_cut_cell_divergence(const FixedVector<Cut_cell_double, 3>& cut_cell_flux,
                                                                  const IJK_Field_local_double& flux_x,
                                                                  const IJK_Field_local_double& flux_y,
                                                                  const IJK_Field_local_double& flux_zmin,
                                                                  const IJK_Field_local_double& flux_zmax,
                                                                  IJK_Field_double& resu, int k_layer, bool add)
{
  assert(&(cut_cell_flux[0].get_cut_cell_disc()) == &(cut_cell_flux[1].get_cut_cell_disc()));
  assert(&(cut_cell_flux[0].get_cut_cell_disc()) == &(cut_cell_flux[2].get_cut_cell_disc()));
  const Cut_cell_FT_Disc& cut_cell_disc = cut_cell_flux[0].get_cut_cell_disc();

  Cut_field_double& cut_field_resu = static_cast<Cut_field_double&>(resu);

  for (int index = cut_cell_disc.get_k_value_index(k_layer); index < cut_cell_disc.get_k_value_index(k_layer+1); index++)
    {
      int n = cut_cell_disc.get_n_from_k_index(index);
      Int3 ijk = cut_cell_disc.get_ijk(n);

      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      if (!cut_cell_disc.get_domaine().within_ghost(i, j, k, 0, 0))
        continue;

      BOUNDARY_FLUX type_boundary_flux = flux_determined_by_boundary_condition_<DIRECTION::Z>(k);
      if (type_boundary_flux != BOUNDARY_FLUX::NOT_DETERMINED_BY_BOUNDARY)
        {
          Cerr << "Le cas d'une cellule diphasique avec flux condition aux limites n'est pas traite" << finl;
          Process::exit();
        }
      else
        {
          for (int phase = 0; phase < 2; phase++)
            {
              const DoubleTabFT_cut_cell& diph_flux_x = (phase == 0) ? cut_cell_flux[0].diph_v_ : cut_cell_flux[0].diph_l_;
              const DoubleTabFT_cut_cell& diph_flux_y = (phase == 0) ? cut_cell_flux[1].diph_v_ : cut_cell_flux[1].diph_l_;
              const DoubleTabFT_cut_cell& diph_flux_z = (phase == 0) ? cut_cell_flux[2].diph_v_ : cut_cell_flux[2].diph_l_;
              DoubleTabFT_cut_cell& diph_resu = (phase == 0) ? cut_field_resu.diph_v_ : cut_field_resu.diph_l_;

              int n_ip1 = cut_cell_disc.get_n(i+1,j,k);
              int n_jp1 = cut_cell_disc.get_n(i,j+1,k);
              int n_kp1 = cut_cell_disc.get_n(i,j,k+1); // ???k

              double fx_centre = diph_flux_x(n);
              double fy_centre = diph_flux_y(n);
              double fz_centre = diph_flux_z(n);

              double fx_right  = (n_ip1 < 0) ? (cut_cell_disc.indic_pure(i+1,j,k) == phase)*flux_x(i+1,j,0)  : diph_flux_x(n_ip1);
              double fy_right  = (n_jp1 < 0) ? (cut_cell_disc.indic_pure(i,j+1,k) == phase)*flux_y(i,j+1,0)  : diph_flux_y(n_jp1);
              double fz_right  = (n_kp1 < 0) ? (cut_cell_disc.indic_pure(i,j,k+1) == phase)*flux_zmax(i,j,0) : diph_flux_z(n_kp1);

              double r = 0;
              r += fx_centre - fx_right;
              r += fy_centre - fy_right;
              r += fz_centre - fz_right;

              if (add)
                {
                  r += diph_resu(n);
                }
              diph_resu(n) = r;
            }
        }
    }
}

Implemente_instanciable_sans_constructeur(OpDiffIJKScalar_double, "OpDiffIJKScalar_double", Operateur_IJK_elem_diff_base_double);

Sortie& OpDiffIJKScalar_double::printOn(Sortie& os) const
{
  //  Operateur_IJK_elem_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffIJKScalar_double::readOn(Entree& is)
{
  //  Operateur_IJK_elem_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffAnisotropicIJKScalar_double, "OpDiffAnisotropicIJKScalar_double", Operateur_IJK_elem_diff_base_double);

Sortie& OpDiffAnisotropicIJKScalar_double::printOn(Sortie& os) const
{
  //  Operateur_IJK_elem_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffAnisotropicIJKScalar_double::readOn(Entree& is)
{
  //  Operateur_IJK_elem_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffVectorialIJKScalar_double, "OpDiffVectorialIJKScalar_double", Operateur_IJK_elem_diff_base_double);

Sortie& OpDiffVectorialIJKScalar_double::printOn(Sortie& os) const
{
  //  Operateur_IJK_elem_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffVectorialIJKScalar_double::readOn(Entree& is)
{
  //  Operateur_IJK_elem_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffVectorialAnisotropicIJKScalar_double, "OpDiffVectorialAnisotropicIJKScalar_double", Operateur_IJK_elem_diff_base_double);

Sortie& OpDiffVectorialAnisotropicIJKScalar_double::printOn(Sortie& os) const
{
  //  Operateur_IJK_elem_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffVectorialAnisotropicIJKScalar_double::readOn(Entree& is)
{
  //  Operateur_IJK_elem_diff_base_double::readOn(is);
  return is;
}

Implemente_instanciable_sans_constructeur(OpDiffStructuralOnlyIJKScalar_double, "OpDiffStructuralOnlyIJKScalar_double", Operateur_IJK_elem_diff_base_double);

Sortie& OpDiffStructuralOnlyIJKScalar_double::printOn(Sortie& os) const
{
  //  Operateur_IJK_elem_diff_base_double::printOn(os);
  return os;
}

Entree& OpDiffStructuralOnlyIJKScalar_double::readOn(Entree& is)
{
  //  Operateur_IJK_elem_diff_base_double::readOn(is);
  return is;
}
