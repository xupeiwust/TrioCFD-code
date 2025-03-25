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
#include <Cut_field.h>
#include <Cut_cell_correction_petites_cellules.h>
#include <Cut_cell_schema_auxiliaire.h>
#include <Maillage_FT_IJK.h>
#include <Facettes_Interp_FT.h>
#include <Cut_cell_diffusion_flux_interface.h>

class IJK_FT_cut_cell;

class Cut_cell_diffusion_auxiliaire : public Cut_cell_schema_auxiliaire
{
  Declare_instanciable(Cut_cell_diffusion_auxiliaire);

public:
  int deactivate_correction_petites_cellules_diffusion_;
  void associer(DoubleTabFT_cut_cell_scalar& flux_interface_efficace);

  DoubleTabFT_cut_cell_scalar *flux_interface_efficace_ptr_ = nullptr;

public:
  void set_param(Param& param);

protected:
  double dying_cells_flux(int num_face, int phase, int n, const Cut_field_vector3_double& cut_field_total_velocity, const Cut_field_double& cut_field) override;
  double small_nascent_cells_flux(int num_face, int phase, int n, const Cut_field_vector3_double& cut_field_total_velocity, const Cut_field_double& cut_field) override;

  const DoubleTabFT_cut_cell_scalar& select_flux_interface(int phase);
};

#endif /* Cut_cell_diffusion_auxiliaire_included */
