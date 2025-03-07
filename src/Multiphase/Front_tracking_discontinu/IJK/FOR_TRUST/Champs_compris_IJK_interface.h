/****************************************************************************
* Copyright (c) 2024, CEA
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

#ifndef Champs_compris_IJK_interface_included
#define Champs_compris_IJK_interface_included

#include <IJK_Field_forward.h>
#include <IJK_Field_vector.h>
#include <Champ_Generique_base.h>  // For Entity and Nature_du_champ

class Motcle;

/*! @brief Similar to Champs_compris_interface but for IJK scalar and vector fields
 */
class Champs_compris_IJK_interface
{
public :
  /** Name / Localisation (elem, face, node) / Nature (vector, scalar) / true=Lagrangian (interface), false=Eulerian
   */
  using FieldInfo_t = std::tuple<Motcle, Entity, Nature_du_champ, bool>;

  virtual inline ~Champs_compris_IJK_interface() {}

  virtual const IJK_Field_double& get_IJK_field(const Motcle& nom)=0;
  virtual const IJK_Field_vector3_double& get_IJK_field_vector(const Motcle& nom)=0;

  // Might be a repetition of what is in Champs_compris_interface, but this is a pure virtual, so OK
  virtual bool has_champ(const Motcle& nom) const=0;

  virtual bool has_champ_vectoriel(const Motcle& nom) const=0;

  // This is not possible in C++, but this method should be implemented every time too:
//  static void Fill_postprocessable_fields(std::vector<FieldInfo_t>& chps)=0;
};

#endif
