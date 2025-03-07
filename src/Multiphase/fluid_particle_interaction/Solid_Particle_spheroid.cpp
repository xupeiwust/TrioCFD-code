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
*****************************************************************************/

#include <Solid_Particle_spheroid.h>
#include <Probleme_FT_Disc_gen.h>


Implemente_instanciable_sans_constructeur(Solid_Particle_spheroid,"Solid_Particle_spheroid",Solid_Particle_base);

Solid_Particle_spheroid::Solid_Particle_spheroid()
{
}

Entree& Solid_Particle_spheroid::readOn (Entree& is)
{
  Solid_Particle_base::readOn(is); // do to first
  set_equivalent_diameter(half_long_axis_spheroid_); // for impact activation distance and Verlet tables detection
  set_equivalent_radius(half_long_axis_spheroid_/2); // for impact activation distance and Verlet tables detection
  // set_diameter(); not coded yet
  // set_volume(); not coded yet
  // set_mass(); not coded yet
  return is;
}

Sortie& Solid_Particle_spheroid::printOn(Sortie& os) const
{
  Cerr << "Error::printOn is not implemented." << finl;
  Process::exit();
  return os;
}

void Solid_Particle_spheroid::set_param(Param& param)
{
  Solid_Particle_base::set_param(param);
  param.ajouter("half_small_axis", &half_small_axis_spheroid_, Param::REQUIRED);
  param.ajouter("half_long_axis", &half_long_axis_spheroid_, Param::REQUIRED);
}
