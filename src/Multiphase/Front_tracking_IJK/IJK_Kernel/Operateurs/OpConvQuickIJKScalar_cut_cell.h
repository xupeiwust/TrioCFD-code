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

#ifndef OpConvQuickScalarIJK_included
#define OpConvQuickScalarIJK_included

#include <Operateur_IJK_elem_conv_base.h>

class OpConvQuickIJKScalar_cut_cell_double : public Operateur_IJK_elem_conv_base_double
{
  Declare_instanciable_sans_constructeur(OpConvQuickIJKScalar_cut_cell_double);

public:
  OpConvQuickIJKScalar_cut_cell_double() : Operateur_IJK_elem_conv_base_double() { };

  void calculer(const IJK_Field_double& field,
                const IJK_Field_double& vx,
                const IJK_Field_double& vy,
                const IJK_Field_double& vz,
                IJK_Field_double& result) override
  {
    Cerr << "The cut cell operators demand the use of calculer_cut_cell instead of calculer." << finl;
    Process::exit();
  }
  void ajouter(const IJK_Field_double& field,
               const IJK_Field_double& vx,
               const IJK_Field_double& vy,
               const IJK_Field_double& vz,
               IJK_Field_double& result) override
  {
    Cerr << "The cut cell operators demand the use of ajouter_cut_cell instead of ajouter." << finl;
    Process::exit();
  }

  const Cut_cell_FT_Disc* get_cut_cell_disc() override
  {
    return &cut_cell_flux_->get_cut_cell_disc();
  }
  DoubleTabFT_cut_cell* get_diph_flux(int phase) override
  {
    return (phase == 0) ? &cut_cell_flux_->diph_v_ : &cut_cell_flux_->diph_l_;
  }

  inline void compute_cut_cell_divergence(int phase, const DoubleTabFT_cut_cell& diph_flux,
                                          const DoubleTabFT_cut_cell* flux_interface_ptr,
                                          const IJK_Field_local_double& flux_x,
                                          const IJK_Field_local_double& flux_y,
                                          const IJK_Field_local_double& flux_zmin,
                                          const IJK_Field_local_double& flux_zmax,
                                          DoubleTabFT_cut_cell& diph_resu, int k_layer, bool add) override
  {
    const Cut_cell_FT_Disc& cut_cell_disc = *get_cut_cell_disc();

    for (int index = cut_cell_disc.get_k_value_index(k_layer); index < cut_cell_disc.get_k_value_index(k_layer+1); index++)
      {
        int n = cut_cell_disc.get_n_from_k_index(index);
        Int3 ijk = cut_cell_disc.get_ijk(n);

        int i = ijk[0];
        int j = ijk[1];
        int k = ijk[2];

        if (!cut_cell_disc.within_ghost(i, j, k, 0, 0))
          continue;

        if (flux_determined_by_wall_<DIRECTION::Z>(k))
          {
            Cerr << "Le cas d'une cellule diphasique avec flux de paroi n'est pas traite" << finl;
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

            if (flux_interface_ptr)
              {
                const DoubleTabFT_cut_cell& flux_interface = *flux_interface_ptr;
                r += (phase == 0) ? -flux_interface(n) : flux_interface(n);
              }

            if(add)
              {
                r += diph_resu(n);
              }
            diph_resu(n) = r;
          }
      }
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
};

#include <OpConvQuickIJKScalar_cut_cell.tpp>

#endif /* OpConvQuickScalarIJK_included */
