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

#ifndef Operateur_IJK_elem_diff_base_TPP_included
#define Operateur_IJK_elem_diff_base_TPP_included

static void copy_boundary_condition(const IJK_Field_local_double& boundary_flux, IJK_Field_local_double& resu)
{
  assert(resu.ni() >= boundary_flux.ni());
  assert(resu.nj() >= boundary_flux.nj());
  // resu is the temporary array where all fluxes are stored before computing divergence,
  // they might have more place than ni and nj because there 1 more flux value dans velocity values
  // to compute divergence
  const int ni = boundary_flux.ni();
  const int nj = boundary_flux.nj();
  for (int j = 0; j < nj; j++)
    for (int i = 0; i < ni; i++)
      resu(i,j,0) = boundary_flux(i,j,0);
}

template <DIRECTION _DIR_>
void Operateur_IJK_elem_diff_base_double::compute_flux_(IJK_Field_local_double& resu, const int k_layer)
{
  const int nx = _DIR_ == DIRECTION::X ? input_field_->ni() + 1 : input_field_->ni() ;
  const int ny = _DIR_ == DIRECTION::Y ? input_field_->nj() + 1 : input_field_->nj();

  ConstIJK_double_ptr input_field(*input_field_, 0, 0, k_layer);
  Simd_double uniform_lambda(1.);
  Simd_double avg_lambda(1.);
  if (is_uniform_ and uniform_lambda_!=0)
    uniform_lambda = Simd_double(*uniform_lambda_);
  /*
   *  M.G: lambda point toward the input field just to initialise *structural_model without error
   *  May not work in further configurations (may be handled in IJK_Thermal classes) when operators
   *  are cast.
   */
  if (is_uniform_)
    lambda_=input_field_;


  /*
   * Gives lambda field as a dummy field (Avoid the creation of a IJK_Field_local_double
   * field in the current scope)
   */
  ConstIJK_double_ptr lambda(is_vectorial_? get_model(_DIR_) : *lambda_, 0, 0, k_layer);
  ConstIJK_double_ptr structural_model(is_structural_ ? get_model(_DIR_) : *lambda_, 0, 0, k_layer);

  IJK_double_ptr resu_ptr(resu, 0, 0, 0);

  if(_DIR_ == DIRECTION::Z)
    {
      // Are we on the wall ?
      const int global_k_layer = k_layer + channel_data_.offset_to_global_k_layer();
      // global index of the layer of flux of the wall
      //  (assume one walls at zmin and zmax)
      const int first_global_k_layer = 0; // index of k_layer when we are on the wall
      // Fluxes in direction k are on the faces, hence, index of last flux is equal to number of elements
      const int last_global_k_layer =  channel_data_.nb_elem_k_tot();

      if (!perio_k_)
        {
          if (global_k_layer == first_global_k_layer)
            {
              // We are on wall at kmin, copy boundary condition fluxes to "resu"
              if (boundary_flux_kmin_) // boundary condition is not zero flux
                copy_boundary_condition(*boundary_flux_kmin_, resu);
              else
                putzero(resu);
              return;
            }
          else if (global_k_layer == last_global_k_layer)
            {
              if (boundary_flux_kmax_) // boundary condition is not zero flux
                copy_boundary_condition(*boundary_flux_kmax_, resu);
              else
                putzero(resu);
              return;
            }
        }
    }

  double d0 = 0., d1 = 0.;
  double surface = 1.;
  switch(_DIR_)
    {
    case DIRECTION::X:
      if (!is_hess_)
        {
          d0 = channel_data_.get_delta_x() * 0.5;
          d1 = d0;
          surface = channel_data_.get_delta_y() * channel_data_.get_delta_z()[k_layer];
        }
      else
        surface = 1 / (channel_data_.get_delta_x() * channel_data_.get_delta_x());
      break;
    case DIRECTION::Y:
      if (!is_hess_)
        {
          d0 = channel_data_.get_delta_y() * 0.5;
          d1 = d0;
          surface = channel_data_.get_delta_x() * channel_data_.get_delta_z()[k_layer];
        }
      else
        surface = 1 / (channel_data_.get_delta_y() * channel_data_.get_delta_y());
      break;
    case DIRECTION::Z:
      if (!is_hess_)
        {
          d0 = channel_data_.get_delta_z()[k_layer-1] * 0.5;
          d1 = channel_data_.get_delta_z()[k_layer] * 0.5;
          surface = channel_data_.get_delta_x() * channel_data_.get_delta_y();
        }
      break;
    }

  Simd_double lambda_m1(uniform_lambda);
  Simd_double lambda_m2(uniform_lambda);
  const int imax = nx;
  const int jmax = ny;
  const int vsize = Simd_double::size();
  for (int j = 0; ; j++)
    {
      for (int i = 0; i < imax; i += vsize)
        {
          Simd_double flux = 0.;
          if (is_structural_)
            {
              Simd_double s_mo, s_mo_dummy;
              structural_model.get_left_center(_DIR_, i, s_mo_dummy, s_mo);
              flux = (-1.) * s_mo * surface;
            }
          else
            {
              Simd_double left_val, right_val;
              input_field.get_left_center(_DIR_, i, left_val, right_val);
              Simd_double zeroVec = 0.;
              Simd_double oneVec = 1.;
              Simd_double minVec = DMINFLOAT;
              Simd_double d = 1.;
              if (!is_hess_)
                {
                  // Fetch conductivity on neighbour cells:
                  if (!is_uniform_)
                    {
                      lambda.get_left_center(_DIR_, i, lambda_m1, lambda_m2);
                    }
                  // geometric avg: (d0+d1) / ( d0 / lambda_m1 + d1 / lambda_m2 ), optimized with only 1 division:
                  Simd_double dsabs = SimdSelect(zeroVec, d0 * lambda_m2 + d1 * lambda_m1, d0 * lambda_m2 + d1 * lambda_m1, (-1) * (d0 * lambda_m2 + d1 * lambda_m1));
                  Simd_double ds = SimdSelect(dsabs, minVec, oneVec, d0 * lambda_m2 + d1 * lambda_m1);
                  if(is_anisotropic_)
                    d = d0 + d1;
                  avg_lambda = SimdSelect(dsabs, minVec, zeroVec, SimdDivideMed(d * lambda_m1 * lambda_m2, ds));
                }
              // thermal flux is positive if going from left to right => -grad(T)
              flux = (left_val - right_val) * avg_lambda * surface;
            }
          resu_ptr.put_val(i, flux);
        }
      // do not execute end_iloop at last iteration (because of assert on valid j+1)
      if (j+1==jmax)
        break;
      input_field.next_j();
      if (!is_uniform_)
        { lambda.next_j(); }
      if(is_structural_)
        structural_model.next_j();
      resu_ptr.next_j();
    }
}

template <DIRECTION _DIR_>
double Operateur_IJK_elem_diff_base_double::compute_flux_local_(int i, int j, int k)
{
  const int dir_i = (_DIR_ == DIRECTION::X);
  const int dir_j = (_DIR_ == DIRECTION::Y);
  const int dir_k = (_DIR_ == DIRECTION::Z);

  IJK_Field_local_double input_field = *input_field_;
  /*
   *  M.G: lambda point toward the input field just to initialise *structural_model without error
   *  May not work in further configurations (may be handled in IJK_Thermal classes) when operators
   *  are cast.
   */
  if (is_uniform_)
    lambda_=input_field_;


  /*
   * Gives lambda field as a dummy field (Avoid the creation of a IJK_Field_local_double
   * field in the current scope)
   */
  IJK_Field_local_double lambda = is_vectorial_? get_model(_DIR_) : *lambda_;
  IJK_Field_local_double structural_model = is_structural_ ? get_model(_DIR_) : *lambda_;

  BOUNDARY_FLUX type_boundary_flux = flux_determined_by_boundary_condition_<_DIR_>(k);
  if (type_boundary_flux != BOUNDARY_FLUX::NOT_DETERMINED_BY_BOUNDARY)
    {
      double flux = compute_flux_local_boundary_condition_<_DIR_>(type_boundary_flux, i, j);
      return flux;
    }
  else
    {
      Vecteur3 surface_d0_d1 = compute_surface_d0_d1_<_DIR_>(k);
      double surface = surface_d0_d1[0];
      double d0 = surface_d0_d1[1];
      double d1 = surface_d0_d1[2];

      double input_left = input_field(i-dir_i,j-dir_j,k-dir_k);
      double input_centre = input_field(i,j,k);
      double lambda_left = lambda(i-dir_i,j-dir_j,k-dir_k);
      double lambda_centre = lambda(i,j,k);
      double struct_model = is_structural_ ? structural_model(i,j,k) : -1;
      double flux = Operateur_IJK_elem_diff_base_double::compute_flux_local_<_DIR_>(d0, d1, surface, input_left, input_centre, lambda_left, lambda_centre, struct_model);
      return flux;
    }
}

template <DIRECTION _DIR_>
double Operateur_IJK_elem_diff_base_double::compute_flux_local_(double d0, double d1, double surface, double input_left, double input_centre, double lambda_left, double lambda_centre, double structural_model)
{
  double uniform_lambda = 1.;
  double avg_lambda = 1.;
  if (is_uniform_ and uniform_lambda_!=0)
    uniform_lambda = *uniform_lambda_;

  double lambda_m1 = uniform_lambda;
  double lambda_m2 = uniform_lambda;
  double flux = 0.;
  if (is_structural_)
    {
      double s_mo = structural_model;
      flux = (-1.) * s_mo * surface;
    }
  else
    {
      double d = 1.;
      if (!is_hess_)
        {
          // Fetch conductivity on neighbour cells:
          if (!is_uniform_)
            {
              lambda_m1 = lambda_left;
              lambda_m2 = lambda_centre;
            }
          // geometric avg: (d0+d1) / ( d0 / lambda_m1 + d1 / lambda_m2 ), optimized with only 1 division:
          double dsabs = (0. < d0 * lambda_m2 + d1 * lambda_m1) ? d0 * lambda_m2 + d1 * lambda_m1 : (-1) * (d0 * lambda_m2 + d1 * lambda_m1);
          double ds = (dsabs < DMINFLOAT) ? 1. : d0 * lambda_m2 + d1 * lambda_m1;
          if(is_anisotropic_)
            d = d0 + d1;
          avg_lambda = (dsabs < DMINFLOAT) ? 0. : (d * lambda_m1 * lambda_m2)/ds;
        }
      // thermal flux is positive if going from left to right => -grad(T)
      flux = (input_left - input_centre) * avg_lambda * surface;
    }
  return flux;
}

template <DIRECTION _DIR_>
BOUNDARY_FLUX Operateur_IJK_elem_diff_base_double::flux_determined_by_boundary_condition_(int k)
{
  if(_DIR_ == DIRECTION::Z)
    {
      // Are we on the wall ?
      const int global_k_layer = k + channel_data_.offset_to_global_k_layer();
      // global index of the layer of flux of the wall
      //  (assume one walls at zmin and zmax)
      const int first_global_k_layer = 0; // index of k when we are on the wall
      // Fluxes in direction k are on the faces, hence, index of last flux is equal to number of elements
      const int last_global_k_layer =  channel_data_.nb_elem_k_tot();

      if (!perio_k_)
        {
          if (global_k_layer == first_global_k_layer)
            {
              return BOUNDARY_FLUX::KMIN_WALL;
            }
          else if (global_k_layer == last_global_k_layer)
            {
              return BOUNDARY_FLUX::KMAX_WALL;
            }
        }
    }
  return BOUNDARY_FLUX::NOT_DETERMINED_BY_BOUNDARY;
}

template <DIRECTION _DIR_>
double Operateur_IJK_elem_diff_base_double::compute_flux_local_boundary_condition_(BOUNDARY_FLUX type_boundary_flux, int i, int j)
{
  if (type_boundary_flux == BOUNDARY_FLUX::KMIN_WALL)
    {
      // We are on wall at kmin, copy boundary condition fluxes to "resu"
      if (boundary_flux_kmin_) // boundary condition is not zero flux
        return (*boundary_flux_kmin_)(i,j,0);
      else
        return 0.;
    }
  else if (type_boundary_flux == BOUNDARY_FLUX::KMAX_WALL)
    {
      if (boundary_flux_kmax_) // boundary condition is not zero flux
        return (*boundary_flux_kmax_)(i,j,0);
      else
        return 0.;
    }
  else
    {
      Cerr << "Unexpected situation in compute_flux_local_boundary_condition_" << finl;
      Process::exit();
      return -1;
    }
}

template <DIRECTION _DIR_>
Vecteur3 Operateur_IJK_elem_diff_base_double::compute_surface_d0_d1_(int k)
{
  double d0 = 0., d1 = 0.;
  double surface = 1.;
  switch(_DIR_)
    {
    case DIRECTION::X:
      if (!is_hess_)
        {
          d0 = channel_data_.get_delta_x() * 0.5;
          d1 = d0;
          surface = channel_data_.get_delta_y() * channel_data_.get_delta_z()[k];
        }
      else
        surface = 1 / (channel_data_.get_delta_x() * channel_data_.get_delta_x());
      break;
    case DIRECTION::Y:
      if (!is_hess_)
        {
          d0 = channel_data_.get_delta_y() * 0.5;
          d1 = d0;
          surface = channel_data_.get_delta_x() * channel_data_.get_delta_z()[k];
        }
      else
        surface = 1 / (channel_data_.get_delta_y() * channel_data_.get_delta_y());
      break;
    case DIRECTION::Z:
      if (!is_hess_)
        {
          d0 = channel_data_.get_delta_z()[k-1] * 0.5;
          d1 = channel_data_.get_delta_z()[k] * 0.5;
          surface = channel_data_.get_delta_x() * channel_data_.get_delta_y();
        }
      break;
    }
  return {surface, d0, d1};
}

template <DIRECTION _DIR_>
void Operateur_IJK_elem_diff_base_double::correct_flux_(IJK_Field_local_double *const flux, int k_layer)
{
  int dir = static_cast<int>(_DIR_);

  IJK_Field_local_double input_field = *input_field_;
  IJK_Field_local_double lambda = is_vectorial_? get_model(_DIR_) : *lambda_;
  IJK_Field_local_double structural_model = is_structural_ ? get_model(_DIR_) : *lambda_;
  Cut_cell_FT_Disc& cut_cell_disc = cut_cell_flux_->get_cut_cell_disc();

  IJK_Field_int& treatment_count = cut_cell_disc.get_treatment_count();
  int new_treatment = cut_cell_disc.new_treatment();

  int backward_receptive_stencil = 0;
  int forward_receptive_stencil = 1;
  assert(backward_receptive_stencil <= cut_cell_disc.get_ghost_size());
  assert(forward_receptive_stencil <= cut_cell_disc.get_ghost_size());

  if (_DIR_ == DIRECTION::Z)
    {
      if (cut_cell_disc.get_splitting().get_grid_geometry().get_periodic_flag(dir))
        {
          const int kmax = cut_cell_disc.get_interfaces().I().nk();
          int n_dir = cut_cell_disc.get_splitting().get_nb_elem_local(dir);
          int n_dir_tot = cut_cell_disc.get_splitting().get_grid_geometry().get_nb_elem_tot(dir);

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

              if (!cut_cell_disc.within_ghost_<_DIR_>(i, j, k, 0, 1))
                continue;

              if (treatment_count(i,j,k) == new_treatment)
                continue;

              {
                int index_ijk_per = 0;
                while (index_ijk_per >= 0)
                  {
                    Int3 ijk = cut_cell_disc.ijk_per_of_index(i, j, k, index_ijk_per);
                    index_ijk_per = cut_cell_disc.next_index_ijk_per(i, j, k, index_ijk_per, 0, 1);
                    assert((index_ijk_per - 1 >= 5 || index_ijk_per < 0) || (k_layer == ijk[2]));

                    treatment_count(ijk[0],ijk[1],ijk[2]) = new_treatment;
                  }
              }

              double indicatrice_centre = cut_cell_disc.get_interfaces().In(i,j,k);
              int n_centre = cut_cell_disc.get_n(i, j, k);

              BOUNDARY_FLUX type_boundary_flux = flux_determined_by_boundary_condition_<_DIR_>(k);
              if (type_boundary_flux != BOUNDARY_FLUX::NOT_DETERMINED_BY_BOUNDARY)
                {
                  // Le flux est dans ce cas determine par la condition aux limites
                  // Aucune correction n'est donc necessaire
                  assert(n_centre < 0); // Le cas d'une cellule diphasique avec flux condition aux limites n'est pas traite
                }
              else
                {
                  assert((n_centre >= 0) || ((indicatrice_centre == 0.) || (indicatrice_centre == 1.)));
                  int phase_min = (n_centre < 0) ? (int)indicatrice_centre : 0;
                  int phase_max = (n_centre < 0) ? (int)indicatrice_centre : 1;

                  for (int phase = phase_min ; phase <= phase_max ; phase++)
                    {
                      const DoubleTabFT_cut_cell& diph_input = (phase == 0) ? input_cut_field_->diph_v_ : input_cut_field_->diph_l_;
                      DoubleTabFT_cut_cell& diph_flux = (phase == 0) ? cut_cell_flux_->diph_v_ : cut_cell_flux_->diph_l_;

                      const int dir_i = (_DIR_ == DIRECTION::X);
                      const int dir_j = (_DIR_ == DIRECTION::Y);
                      const int dir_k = (_DIR_ == DIRECTION::Z);

                      int n_left = cut_cell_disc.get_n(i-dir_i,j-dir_j,k-dir_k);

                      double indicatrice_left = cut_cell_disc.get_interfaces().In(i-dir_i,j-dir_j,k-dir_k);

                      double bar_dir_left = cut_cell_disc.get_interfaces().get_barycentre_phase1_next()[dir](i-dir_i,j-dir_j,k-dir_k);
                      bar_dir_left = (phase == 0) ? IJK_Interfaces::opposing_barycentre(bar_dir_left, indicatrice_left) : bar_dir_left;
                      assert((n_left >= 0) || (bar_dir_left == .5));

                      double bar_dir_centre = cut_cell_disc.get_interfaces().get_barycentre_phase1_next()[dir](i,j,k);
                      bar_dir_centre = (phase == 0) ? IJK_Interfaces::opposing_barycentre(bar_dir_centre, indicatrice_centre) : bar_dir_centre;
                      assert((n_centre >= 0) || (bar_dir_centre == .5));

                      Vecteur3 surface_d0_d1 = compute_surface_d0_d1_<_DIR_>(k);

                      const DoubleTabFT_cut_cell_vector3& indicatrice_surfacique = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face();
                      double indicatrice_surface = (n_centre < 0) ? 1. : ((phase == 0) ? 1 - indicatrice_surfacique(n_centre,dir) : indicatrice_surfacique(n_centre,dir));
                      assert((n_centre >= 0) || (indicatrice_surface == 1.));

                      const double surface = surface_d0_d1[0] * indicatrice_surface;
                      double d0 = surface_d0_d1[1] * 2*(1 - bar_dir_left);
                      double d1 = surface_d0_d1[2] * 2*(bar_dir_centre);

                      int phase_left = (n_left < 0) ? (int)indicatrice_left : phase;
                      assert((phase_left == phase) || (indicatrice_surface == 0.));

                      double input_left = (n_left < 0) ? input_field(i-dir_i,j-dir_j,k-dir_k) : diph_input(n_left);
                      double input_centre = (n_centre < 0) ? input_field(i,j,k) : diph_input(n_centre);

                      double lambda_value = (phase == 0) ? *uniform_lambda_vapour_ : *uniform_lambda_liquid_;
                      assert((phase != phase_left) || (n_left >= 0) || (lambda_value == lambda(i-dir_i,j-dir_j,k-dir_k))); // La cellule est pure, lambda(i,j,k) doit donc etre celui de la phase
                      assert((phase != phase_left) || (n_centre >= 0) || (lambda_value == lambda(i,j,k))); // La cellule est pure, lambda(i,j,k) doit donc etre celui de la phase

                      double struct_model = is_structural_ ? structural_model(i,j,k) : -1;

                      double old_indicatrice_left = cut_cell_disc.get_interfaces().I(i-dir_i,j-dir_j,k-dir_k);
                      double old_indicatrice_centre = cut_cell_disc.get_interfaces().I(i,j,k);

                      double flux_value;
                      if ((cut_cell_disc.get_interfaces().devient_pure(old_indicatrice_centre, indicatrice_centre)) && (cut_cell_disc.get_interfaces().devient_pure(old_indicatrice_left, indicatrice_left)) && ((int)(1 - indicatrice_centre) == phase) && ((int)(1 - indicatrice_left) == phase))
                        {
                          flux_value = 0.;
                        }
                      else if ((cut_cell_disc.get_interfaces().devient_pure(old_indicatrice_centre, indicatrice_centre)) && ((int)(1 - indicatrice_centre) == phase))
                        {
                          flux_value = 0.;
                        }
                      else if ((cut_cell_disc.get_interfaces().devient_pure(old_indicatrice_left, indicatrice_left)) && ((int)(1 - indicatrice_left) == phase))
                        {
                          flux_value = 0.;
                        }
                      else if ((cut_cell_disc.get_interfaces().devient_diphasique(old_indicatrice_centre, indicatrice_centre)) && (cut_cell_disc.get_interfaces().devient_diphasique(old_indicatrice_left, indicatrice_left)) && ((int)(1 - old_indicatrice_centre) == phase) && ((int)(1 - old_indicatrice_left) == phase))
                        {
                          flux_value = 0.;
                        }
                      else if ((cut_cell_disc.get_interfaces().devient_diphasique(old_indicatrice_centre, indicatrice_centre)) && ((int)(1 - old_indicatrice_centre) == phase))
                        {
                          flux_value = 0.;
                        }
                      else if ((cut_cell_disc.get_interfaces().devient_diphasique(old_indicatrice_left, indicatrice_left)) && ((int)(1 - old_indicatrice_left) == phase))
                        {
                          flux_value = 0.;
                        }
                      else
                        {
                          flux_value = Operateur_IJK_elem_diff_base_double::compute_flux_local_<_DIR_>(d0, d1, surface, input_left, input_centre, lambda_value, lambda_value, struct_model);
                        }

                      if (n_centre < 0) // Si la cellule est pure
                        {
                          int index_ijk_per = 0;
                          while (index_ijk_per >= 0)
                            {
                              Int3 ijk = cut_cell_disc.ijk_per_of_index(i, j, k, index_ijk_per);
                              index_ijk_per = cut_cell_disc.next_index_ijk_per(i, j, k, index_ijk_per, 0, 1);
                              assert((index_ijk_per - 1 >= 5 || index_ijk_per < 0) || (k_layer == ijk[2]));

                              (*flux)(ijk[0],ijk[1],0) = flux_value;
                            }
                        }
                      else
                        {
                          diph_flux(n_centre, dir) = flux_value;
                          if ((n_left < 0) && (phase == phase_left))
                            {
                              int index_ijk_per = 0;
                              while (index_ijk_per >= 0)
                                {
                                  Int3 ijk = cut_cell_disc.ijk_per_of_index(i, j, k, index_ijk_per);
                                  index_ijk_per = cut_cell_disc.next_index_ijk_per(i, j, k, index_ijk_per, 0, 1);
                                  assert((index_ijk_per - 1 >= 5 || index_ijk_per < 0) || (k_layer == ijk[2]));

                                  (*flux)(ijk[0],ijk[1],0) = flux_value;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#endif

