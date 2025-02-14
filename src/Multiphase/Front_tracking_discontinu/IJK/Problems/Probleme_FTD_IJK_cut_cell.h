/****************************************************************************
* Copyright (c) 2025, CEA
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

#ifndef Probleme_FTD_IJK_cut_cell_included
#define Probleme_FTD_IJK_cut_cell_included

#include <Probleme_FTD_IJK_base.h>
#include <Facettes_Interp_FT.h>
#include <Cut_cell_FT_Disc.h>

class Probleme_FTD_IJK_cut_cell : public Probleme_FTD_IJK_base
{
  Declare_instanciable(Probleme_FTD_IJK_cut_cell) ;
public:
  void set_param(Param& param) override;

  void deplacer_interfaces(const double timestep, const int rk_step, ArrOfDouble& var_volume_par_bulle, const int first_step_interface_smoothing) override;
  void deplacer_interfaces_rk3(const double timestep, const int rk_step, ArrOfDouble& var_volume_par_bulle) override;

  void update_indicator_field() override;
  void update_twice_indicator_field() override;

  Cut_cell_FT_Disc* get_cut_cell_disc() override { return &cut_cell_disc_; }

  const Cut_field_vector3_double& get_cut_field_velocity() const;

  /*! Getter des objets Facettes_Interp_FT, permettant d'acceder aux indices et coefficients des points d'interpolation a une certaine distance des facettes de l'interface
   */
  const Facettes_Interp_FT& get_cut_cell_facettes_interpolation() const
  {
    return cut_cell_facettes_interpolation_;
  }

  void cut_cell_perform_interpolation_facettes()
  {
    cut_cell_facettes_interpolation_.cut_cell_perform_interpolation_facettes_next(get_interface().next());
  }


protected:
  Cut_cell_FT_Disc cut_cell_disc_;
  DoubleTabFT_cut_cell_vector3 velocity_interface_;

  TYPE_SURFACE_EFFICACE_FACE type_surface_efficace_face_ = TYPE_SURFACE_EFFICACE_FACE::NON_INITIALISE;
  TYPE_SURFACE_EFFICACE_INTERFACE type_surface_efficace_interface_ = TYPE_SURFACE_EFFICACE_INTERFACE::NON_INITIALISE;

  double seuil_indicatrice_petite_fixe_ = -1;   // Valeur predefinie du seuil
  double seuil_indicatrice_petite_facsec_ = -1; // Valeur du seuil exprimee relativement au facsec

  // Stockage des indices et coefficients des points d'interpolation a une certaine distance des facettes de l'interface
  Facettes_Interp_FT cut_cell_facettes_interpolation_;

  void initialize() override;

};

#endif /* Probleme_FTD_IJK_cut_cell_included */
