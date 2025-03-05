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

#include <Solid_Particle_sphere.h>
#include <Probleme_FT_Disc_gen.h>


Implemente_instanciable_sans_constructeur(Solid_Particle_sphere,"Solid_Particle_sphere",Solid_Particle_base);

Solid_Particle_sphere::Solid_Particle_sphere()
{
}

Entree& Solid_Particle_sphere::readOn (Entree& is)
{
  set_diameter(2*radius_);
  set_equivalent_diameter(diameter_);
  set_equivalent_radius(radius_);
  set_volume(4 * M_PI * pow(radius_, 3) / 3);
  Solid_Particle_base::readOn(is);
  return is;
}

Sortie& Solid_Particle_sphere::printOn(Sortie& os) const
{
  Cerr << "Error::printOn is not implemented." << finl;
  Process::exit();
  return os;
}

void Solid_Particle_sphere::set_param(Param& param)
{
  Solid_Particle_base::set_param(param);
  param.ajouter("radius", &radius_, Param::REQUIRED); // XD_ADD_P radius of a spherical particle in the case of a sphere flow;
}
