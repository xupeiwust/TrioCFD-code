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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Cut_cell_diffusion_auxiliaire.h
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Cut_cell_diffusion_auxiliaire_included
#define Cut_cell_diffusion_auxiliaire_included

#include <FixedVector.h>
#include <IJK_Field.h>
#include <Param.h>
#include <IJK_Interfaces.h>
#include <Champ_diphasique.h>
#include <Cut_cell_correction_petites_cellules.h>
#include <Cut_cell_schema_auxiliaire.h>
#include <Maillage_FT_IJK.h>

class IJK_FT_cut_cell;

enum class METHODE_FLUX_INTERFACE : int
{
  NON_INITIALISE,  // Valeur invalide par defaut, pour forcer le choix
  INTERP_PURE,     // Methode d'interpolation n'utilisant pas les donnees cut-cell (cf. Aymeric)
  INTERP_CUT_CELL  // Methode d'interpolation utilisant les donnees cut-cell
};

struct Facettes_data
{
  DoubleTabFT centre; // Centre de la facette
  DoubleTabFT liqu_1; // A une certaine distance du cote liquide
  DoubleTabFT vap_1;  // A une certaine distance du cote vapeur
  DoubleTabFT liqu_2; // A une certaine distance (plus grande) du cote liquide
  DoubleTabFT vap_2;  // A une certaine distance (plus grande) du cote vapeur
};

class Cut_cell_diffusion_auxiliaire : public Cut_cell_schema_auxiliaire
{
  Declare_instanciable(Cut_cell_diffusion_auxiliaire);
public:
  FixedVector<IJK_Field_double, 2> flux_interface_ns_; // Note : contrairement au champs de IJK_Interfaces, on a toujours [0] = old et [1] = next pour le flux_interface
  FixedVector<IJK_Field_double, 2> flux_interface_ft_;
  DoubleTabFT_cut_cell_scalar flux_interface_efficace_;
  int deactivate_correction_petites_cellules_diffusion_;

public:
  void set_param(Param& param);

  double dying_cells_flux(int num_face, int phase, int n, const Cut_field_vector3_double& cut_field_total_velocity, const Cut_field_double& cut_field_temperature) override;
  double small_nascent_cells_flux(int num_face, int phase, int n, const Cut_field_vector3_double& cut_field_total_velocity, const Cut_field_double& cut_field_temperature) override;

  void calculer_flux_interface(bool next_time,
                               double lambda_liquid,
                               double lambda_vapour,
                               Facettes_data& coord_facettes,
                               Facettes_data& interfacial_temperature,
                               DoubleTabFT& interfacial_phin_ai,
                               const Cut_field_double& cut_field_temperature,
                               REF(IJK_FT_cut_cell)& ref_ijk_ft,
                               const IJK_Field_double& temperature_ns,
                               IJK_Field_double& temperature_ft);
  void calculer_flux_interface_next(double lambda_liquid,
                                    double lambda_vapour,
                                    Facettes_data& coord_facettes,
                                    Facettes_data& interfacial_temperature,
                                    DoubleTabFT& interfacial_phin_ai,
                                    const Cut_field_double& cut_field_temperature,
                                    REF(IJK_FT_cut_cell)& ref_ijk_ft,
                                    const IJK_Field_double& temperature_ns,
                                    IJK_Field_double& temperature_ft);
  void calculer_flux_interface_old(double lambda_liquid,
                                   double lambda_vapour,
                                   Facettes_data& coord_facettes,
                                   Facettes_data& interfacial_temperature,
                                   DoubleTabFT& interfacial_phin_ai,
                                   const Cut_field_double& cut_field_temperature,
                                   REF(IJK_FT_cut_cell)& ref_ijk_ft,
                                   const IJK_Field_double& temperature_ns,
                                   IJK_Field_double& temperature_ft);
  void calculer_flux_interface_efficace();
  void compute_interfacial_temperature2(bool next_time,
                                        double lambda_liquid,
                                        double lambda_vapour,
                                        const IJK_Field_double& temperature_ft,
                                        const IJK_Grid_Geometry& geom,
                                        const Maillage_FT_IJK& maillage,
                                        Facettes_data& coord_facettes,
                                        Facettes_data& interfacial_temperature,
                                        DoubleTabFT& flux_normal_interp);
  void compute_interfacial_temperature_cut_cell(bool next_time,
                                                double lambda_liquid,
                                                double lambda_vapour,
                                                const Cut_field_double& cut_field_temperature,
                                                const IJK_Grid_Geometry& geom,
                                                const Maillage_FT_IJK& maillage,
                                                Facettes_data& coord_facettes,
                                                Facettes_data& interfacial_temperature,
                                                DoubleTabFT& flux_normal_interp);

  void ajout_flux_interface_a_divergence_simple(Cut_field_double& cut_field_div_coeff_grad_T_volume);

  void calcul_temperature_flux_interface(const IJK_Field_double& temperature, const double ldal, const double ldav,
                                         const double dist, const DoubleTab& positions, const DoubleTab& normal_on_interf,
                                         DoubleTabFT& temperature_interp,
                                         DoubleTabFT& flux_normal_interp,
                                         DoubleTabFT& temp_liqu,
                                         DoubleTabFT& temp_vap,
                                         DoubleTab& coo_liqu,
                                         DoubleTab& coo_vap);
  void calcul_temperature_flux_interface_cut_cell(bool next_time,
                                                  const Cut_field_double& temperature, const double ldal, const double ldav,
                                                  const double dist, const DoubleTab& positions, const DoubleTab& normal_on_interf,
                                                  DoubleTabFT& temperature_interp,
                                                  DoubleTabFT& flux_normal_interp,
                                                  DoubleTabFT& temp_liqu,
                                                  DoubleTabFT& temp_vap,
                                                  DoubleTab& coo_liqu,
                                                  DoubleTab& coo_vap);
  void calcul_temperature_flux_interface_second_order(const IJK_Field_double& temperature, const double ldal, const double ldav,
                                                      const double dist_1, const double dist_2, const DoubleTab& positions, const DoubleTab& normal_on_interf,
                                                      DoubleTabFT& temperature_interp,
                                                      DoubleTabFT& flux_normal_interp,
                                                      DoubleTabFT& temp_liqu_1,
                                                      DoubleTabFT& temp_vap_1,
                                                      DoubleTabFT& temp_liqu_2,
                                                      DoubleTabFT& temp_vap_2,
                                                      DoubleTab& coo_liqu_1,
                                                      DoubleTab& coo_vap_1,
                                                      DoubleTab& coo_liqu_2,
                                                      DoubleTab& coo_vap_2);
  void calcul_temperature_flux_interface_cut_cell_second_order(bool next_time,
                                                               const Cut_field_double& temperature, const double ldal, const double ldav,
                                                               const double dist_1, const double dist_2, const DoubleTab& positions, const DoubleTab& normal_on_interf,
                                                               DoubleTabFT& temperature_interp,
                                                               DoubleTabFT& flux_normal_interp,
                                                               DoubleTabFT& temp_liqu_1,
                                                               DoubleTabFT& temp_vap_1,
                                                               DoubleTabFT& temp_liqu_2,
                                                               DoubleTabFT& temp_vap_2,
                                                               DoubleTab& coo_liqu_1,
                                                               DoubleTab& coo_vap_1,
                                                               DoubleTab& coo_liqu_2,
                                                               DoubleTab& coo_vap_2);

protected:
  int second_order_diffusion_interface_;
  double scaled_distance_flux_interface_; // Distance a l'interface de l'interpolation utilisee pour calculer le flux a l'interface
  double scaled_distance_second_point_flux_interface_; // Distance a l'interface de l'interpolation utilisee pour calculer le flux a l'interface
  METHODE_FLUX_INTERFACE methode_flux_interface_;
};

#endif /* Cut_cell_diffusion_auxiliaire_included */
