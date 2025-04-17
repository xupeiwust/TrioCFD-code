/****************************************************************************
* Copyright (c) 2023, CEA
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
/*! @brief class Op_Diff_K_Omega_VEF_Face Cette classe represente l'operateur de diffusion turbulente
 *
 *    La discretisation est VEF
 *   Les methodes pour l'implicite sont codees.
 *
 */
#ifndef Op_Diff_K_Omega_VEF_Face_included
#define Op_Diff_K_Omega_VEF_Face_included

#include <Op_Dift_VEF_Face_Gen.h>
#include <Op_Diff_K_Omega_VEF_base.h>


class Domaine_Cl_VEF;
class Champ_P1NC;
class Champ_Don_base;

//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Op_Diff_K_Omega_VEF_Face
//
//////////////////////////////////////////////////////////////////////////////

class Op_Diff_K_Omega_VEF_Face : public Op_Diff_K_Omega_VEF_base, public Op_Dift_VEF_Face_Gen<Op_Diff_K_Omega_VEF_Face>
{

  Declare_instanciable(Op_Diff_K_Omega_VEF_Face);

public:

  DoubleTab& ajouter(const DoubleTab&, DoubleTab&) const override; // pour l'explicite

  void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const override; // pour l'implicite

  void modifier_pour_Cl(Matrice_Morse& matrice, DoubleTab& secmem) const override;

};

#endif
