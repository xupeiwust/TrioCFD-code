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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Terme_Boussinesq_Sensibility_VEFPreP1B_Face.h
// Directory : $SENSITIVITY_ANALYSIS_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Terme_Boussinesq_Sensibility_VEFPreP1B_Face_included
#define Terme_Boussinesq_Sensibility_VEFPreP1B_Face_included

#include <Terme_Boussinesq_VEFPreP1B_Face.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Terme_Boussinesq_Sensibility_VEFPreP1B_Face
//
// <Description of class Terme_Boussinesq_Sensibility_VEFPreP1B_Face>
//
/////////////////////////////////////////////////////////////////////////////

class Terme_Boussinesq_Sensibility_VEFPreP1B_Face : public Terme_Boussinesq_VEFPreP1B_Face
{

  Declare_instanciable( Terme_Boussinesq_Sensibility_VEFPreP1B_Face ) ;

public :

  DoubleTab& ajouter(DoubleTab& ) const override ;

};

#endif /* Terme_Boussinesq_Sensibility_VEFPreP1B_Face_included */