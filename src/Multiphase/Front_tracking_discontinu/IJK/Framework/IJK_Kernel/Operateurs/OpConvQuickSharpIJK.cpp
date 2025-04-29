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

#include <OpConvQuickSharpIJK.h>

Implemente_instanciable(OpConvQuickSharpIJK_double, "OpConvQuickSharpIJK_double", Operateur_IJK_faces_conv_base_double);

Sortie& OpConvQuickSharpIJK_double::printOn(Sortie& os) const
{
  return os;
}

Entree& OpConvQuickSharpIJK_double::readOn(Entree& is)
{
  return is;
}

void OpConvQuickSharpIJK_double::initialize(const IJK_Splitting& splitting)
{
  Operateur_IJK_faces_conv_base_double::initialize(splitting);
  delta_x_ = splitting.get_grid_geometry().get_constant_delta(DIRECTION_I);
  delta_y_ = splitting.get_grid_geometry().get_constant_delta(DIRECTION_J);
  delta_z_ = splitting.get_grid_geometry().get_constant_delta(DIRECTION_K);
}
