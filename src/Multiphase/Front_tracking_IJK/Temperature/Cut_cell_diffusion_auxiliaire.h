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
#include <IJK_Interfaces.h>
#include <Champ_diphasique.h>
#include <Cut_cell_correction_petites_cellules.h>
#include <Maillage_FT_IJK.h>

class IJK_FT_cut_cell;

enum class METHODE_FLUX_INTERFACE : int
{
  NON_INITIALISE,  // Valeur invalide par defaut, pour forcer le choix
  INTERP_PURE,     // Methode d'interpolation n'utilisant pas les donnees cut-cell (cf. Aymeric)
  INTERP_CUT_CELL, // Methode d'interpolation utilisant les donnees cut-cell
  LOCAL_CELLULE    // Methode utilisant les points de la cellule, en negligeant les variations tangentes
};

enum class ETALEMENT_DIFFUSION : int
{
  AUCUN_ETALEMENT,             // Choix normal, par defaut : pas d'etalement particulier
  FLUX_INTERFACE_UNIQUEMENT,   // Etalement de la contribution du flux a l'interface
  DIVERGENCE_UNIQUEMENT,       // Etalement de la divergence des flux diffusifs
  FLUX_INTERFACE_ET_DIVERGENCE // Etalement de la contribution du flux a l'interface et de la divergence  des flux diffusifs
};

class Cut_cell_diffusion_auxiliaire
{
public:
  IJK_Field_double flux_interface_ns_;
  IJK_Field_double flux_interface_ft_;
  DoubleTabFT_cut_cell_scalar flux_interface_;

public:
  void calculer_flux_interface(METHODE_FLUX_INTERFACE methode_flux_interface,
                               double scaled_distance,
                               double lambda_liquid,
                               double lambda_vapour,
                               ArrOfDouble& interfacial_temperature,
                               ArrOfDouble& interfacial_phin_ai,
                               Cut_field_scalar& cut_field_temperature,
                               REF(IJK_FT_cut_cell)& ref_ijk_ft,
                               const IJK_Field_double& temperature_ns,
                               IJK_Field_double& temperature_ft);
  void compute_interfacial_temperature2(double scaled_distance,
                                        double lambda_liquid,
                                        double lambda_vapour,
                                        const IJK_Field_double& temperature_ft,
                                        const IJK_Grid_Geometry& geom,
                                        const Maillage_FT_IJK& maillage,
                                        ArrOfDouble& interfacial_temperature,
                                        ArrOfDouble& flux_normal_interp);
  void compute_interfacial_temperature_cut_cell(double scaled_distance,
                                                double lambda_liquid,
                                                double lambda_vapour,
                                                Cut_field_scalar& cut_field_temperature,
                                                const IJK_Grid_Geometry& geom,
                                                const Maillage_FT_IJK& maillage,
                                                ArrOfDouble& interfacial_temperature,
                                                ArrOfDouble& flux_normal_interp);
  void compute_interfacial_temperature_local_normal(double lambda_liquid,
                                                    double lambda_vapour,
                                                    Cut_field_scalar& cut_field_temperature,
                                                    const Maillage_FT_IJK& maillage,
                                                    const IJK_Splitting& s,
                                                    const IJK_Interfaces& interfaces,
                                                    ArrOfDouble& interfacial_temperature,
                                                    ArrOfDouble& flux_normal_interp);

  void ajout_flux_interface_a_divergence_simple(Cut_field_scalar& cut_field_div_coeff_grad_T_volume);
  void ajout_flux_interface_a_divergence_etale(Cut_field_scalar& cut_field_div_coeff_grad_T_volume);
  void etalement_divergence_flux_diffusifs(Cut_field_scalar& cut_field_div_coeff_grad_T_volume, Cut_field_scalar& cut_field_div_coeff_grad_T_volume_temp);

  void add_diffusion_small_cells(CORRECTION_PETITES_CELLULES diffusion_petites_cellules, const Cut_field_scalar& cut_field_temperature_post_convection, Cut_field_scalar& cut_field_temperature);

  void calcul_temperature_flux_interface(const IJK_Field_double& temperature, const double ldal, const double ldav,
                                         const double dist, const DoubleTab& positions, const DoubleTab& normal_on_interf,
                                         ArrOfDouble& temperature_interp,
                                         ArrOfDouble& flux_normal_interp,
                                         ArrOfDouble& temp_liqu,
                                         ArrOfDouble& temp_vap,
                                         DoubleTab& coo_liqu,
                                         DoubleTab& coo_vap);
  void calcul_temperature_flux_interface_cut_cell(Cut_field_scalar& temperature, const double ldal, const double ldav,
                                                  const double dist, const DoubleTab& positions, const DoubleTab& normal_on_interf,
                                                  ArrOfDouble& temperature_interp,
                                                  ArrOfDouble& flux_normal_interp,
                                                  ArrOfDouble& temp_liqu,
                                                  ArrOfDouble& temp_vap,
                                                  DoubleTab& coo_liqu,
                                                  DoubleTab& coo_vap);
  void calcul_temperature_flux_interface_local_normal(Cut_field_scalar& temperature, const double ldal, const double ldav,
                                                      const Vecteur3& position_centre_cell, const Vecteur3& positions, const Vecteur3& normal_on_interf,
                                                      double& temperature_interp, double& flux_normal_interp);
};

#endif /* Cut_cell_diffusion_auxiliaire_included */
