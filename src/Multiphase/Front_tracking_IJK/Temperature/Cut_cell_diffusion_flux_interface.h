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

#ifndef Cut_cell_diffusion_flux_interface_included
#define Cut_cell_diffusion_flux_interface_included

#include <FixedVector.h>
#include <IJK_Field.h>
#include <Cut_field.h>
#include <Maillage_FT_IJK.h>
#include <Facettes_Interp_FT.h>

enum class METHODE_FLUX_INTERFACE : int
{
  NON_INITIALISE,               // Valeur invalide par defaut, pour forcer le choix
  INTERP_PURE,                  // Methode d'interpolation n'utilisant pas les donnees cut-cell (cf. Aymeric)
  INTERP_PURE_NO_JUMP,          //  * variante negligeant le saut du gradient a l'interface
  INTERP_CUT_CELL,              // Methode d'interpolation utilisant les donnees cut-cell
  INTERP_CUT_CELL_NO_JUMP       //  * variante negligeant le saut du gradient a l'interface
};

void calculer_flux_interface_next(METHODE_FLUX_INTERFACE methode_flux_interface,
                                  double lambda_liquid,
                                  double lambda_vapour,
                                  DoubleTabFT& interfacial_temperature,
                                  DoubleTabFT& interfacial_phin_ai,
                                  const Cut_field_double& cut_field_temperature,
                                  const Facettes_Interp_FT& cut_cell_facettes_interpolation,
                                  IJK_Field_double& flux_interface_ft);
void calculer_flux_interface_old(METHODE_FLUX_INTERFACE methode_flux_interface,
                                 double lambda_liquid,
                                 double lambda_vapour,
                                 DoubleTabFT& interfacial_temperature,
                                 DoubleTabFT& interfacial_phin_ai,
                                 const Cut_field_double& cut_field_temperature,
                                 const Facettes_Interp_FT& cut_cell_facettes_interpolation,
                                 IJK_Field_double& flux_interface_ft);

void calculer_flux_interface_efficace(const IJK_Field_double& flux_interface_ns_old, const IJK_Field_double& flux_interface_ns_next, DoubleTabFT_cut_cell_scalar& flux_interface_efficace);

void ajout_flux_interface_a_divergence(const DoubleTabFT_cut_cell_scalar& flux_interface_efficace, Cut_field_double& cut_field_d_temperature);

#endif /* Cut_cell_diffusion_flux_interface_included */
