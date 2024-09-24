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
// File      : Cut_cell_convection_auxiliaire.h
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Cut_cell_convection_auxiliaire_included
#define Cut_cell_convection_auxiliaire_included

#include <IJK_Field_vector.h>
#include <IJK_Field.h>
#include <Param.h>
#include <Champ_diphasique.h>
#include <Cut_cell_correction_petites_cellules.h>
#include <Cut_cell_schema_auxiliaire.h>
#include <Maillage_FT_IJK.h>

class IJK_FT_cut_cell;

enum class CUT_CELL_SCHEMA_CONVECTION : int
{
  QUICK_OU_CENTRE2_STENCIL,                  // Utilise le schema quick si le stencil est disponible, le schema centre2 sinon
  QUICK_OU_CENTRE2_PERPENDICULAR_DISTANCE,   // Utilise le schema quick si le stencil est disponible et si la distance perpendiculaire n'est pas grande, le schema centre2 sinon
  QUICK_OU_LINEAIRE2_STENCIL,                // Utilise le schema quick si le stencil est disponible, le schema lineaire2 sinon
  QUICK_OU_LINEAIRE2_PERPENDICULAR_DISTANCE, // Utilise le schema quick si le stencil est disponible et si la distance perpendiculaire n'est pas grande, le schema lineaire2 sinon
  QUICK_OU_AMONT_STENCIL,                    // Utilise le schema quick si le stencil est disponible, le schema amont sinon
  QUICK_OU_AMONT_PERPENDICULAR_DISTANCE,     // Utilise le schema quick si le stencil est disponible et si la distance perpendiculaire n'est pas grande, le schema amont sinon
  CENTRE2,                                   // Utilise toujours le schema centre2
  LINEAIRE2,                                 // Utilise toujours le schema lineaire2
  AMONT                                      // Utilise toujours le schema amont
};

struct Cut_cell_conv_scheme
{
  CUT_CELL_SCHEMA_CONVECTION scheme;
};


class Cut_cell_convection_auxiliaire : public Cut_cell_schema_auxiliaire
{
  Declare_instanciable(Cut_cell_convection_auxiliaire);

public:
  void set_param(Param& param);

  double dying_cells_flux(int num_face, int phase, int n, const Cut_field_vector3_double& cut_field_total_velocity, const Cut_field_double& cut_field_temperature) override;
  double small_nascent_cells_flux(int num_face, int phase, int n, const Cut_field_vector3_double& cut_field_total_velocity, const Cut_field_double& cut_field_temperature) override;

protected:
};

#endif /* Cut_cell_convection_auxiliaire_included */
