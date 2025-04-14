/****************************************************************************
* Copyright (c) 2025, CEA
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

#include <Param.h>
#include <Solid_Particle_base.h>

// XD Solid_Particle_base fluide_incompressible Solid_Particle_base -1 base particle type for collision model
Implemente_base_sans_constructeur(Solid_Particle_base,"Solid_Particle_base",Fluide_Incompressible);

Solid_Particle_base::Solid_Particle_base()
{}

Sortie& Solid_Particle_base::printOn(Sortie& os) const
{
  os << "{"                  << finl;
  os << "e_dry " << e_dry_   << finl;
  os << "}"                  << finl;
  return Fluide_Incompressible::printOn(os);
}

Entree& Solid_Particle_base::readOn(Entree& is)
{
  Fluide_Incompressible::readOn(is);
  return is;
}

void Solid_Particle_base::set_param(Param& param)
{
  Fluide_Incompressible::set_param(param);
  param.ajouter("e_dry", &e_dry_,Param::REQUIRED); // XD_ADD_P double dry coefficient
}
