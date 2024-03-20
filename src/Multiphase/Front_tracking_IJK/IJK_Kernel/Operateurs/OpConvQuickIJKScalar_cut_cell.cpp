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

#include <OpConvQuickIJKScalar_cut_cell.h>

Implemente_instanciable_sans_constructeur(OpConvQuickIJKScalar_cut_cell_double, "OpConvQuickIJKScalar_cut_cell_double", Operateur_IJK_elem_conv_base_double);

Sortie& OpConvQuickIJKScalar_cut_cell_double::printOn(Sortie& os) const
{
  return os;
}

Entree& OpConvQuickIJKScalar_cut_cell_double::readOn(Entree& is)
{
  return is;
}

void OpConvQuickIJKScalar_cut_cell_double::correct_flux(IJK_Field_local_double *const flux, const int k_layer, const int dir)
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
      Cerr << "Unexpected value of dir in OpConvQuickIJKScalar_cut_cell_double::correct_flux" << finl;
      Process::exit();
    }
}

