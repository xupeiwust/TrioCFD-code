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
// File      : Extraire_surface_ALE.h
// Directory : $TRIOCFD_ROOT/src/Fluid_Structure_Interaction/ALE
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Extraire_surface_ALE_included
#define Extraire_surface_ALE_included

#include <Extraire_surface.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Extraire_surface_ALE
//
// <Description of class Extraire_surface_ALE>
//
/////////////////////////////////////////////////////////////////////////////

class Extraire_surface_ALE : public Extraire_surface
{

  Declare_instanciable( Extraire_surface_ALE ) ;

public :
  Entree& interpreter_(Entree&) override;
  void extraire_surface(Domaine& domaine_surfacique,const Domaine& domaine_volumique, const Nom& nom_domaine_surfacique, const Domaine_VF& domaine_vf, const Nom& expr_elements,const Nom& expr_faces, int avec_les_bords, const Noms& noms_des_bords) ;

protected :

};

#endif /* Extraire_surface_ALE_included */
