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

#include <Solid_Particle.h>
#include <Param.h>

Implemente_instanciable_sans_constructeur(Solid_Particle,"Solid_Particle",Fluide_Incompressible);


Solid_Particle::Solid_Particle()
{}

Sortie& Solid_Particle::printOn(Sortie& os) const
{
  os << "{"                  << finl;
  os << "e_dry " << e_dry_   << finl;
  os << "}"                  << finl;
  return Fluide_Incompressible::printOn(os);
}

Entree& Solid_Particle::readOn(Entree& is)
{
  Fluide_Incompressible::readOn(is);
  set_diameter_sphere(2*radius_sphere_);
  set_volume_sphere(4 * M_PI * pow(radius_sphere_, 3) / 3);
  const double solid_density = masse_volumique().valeurs()(0, 0);
  set_mass_sphere(volume_sphere_*solid_density);
  return is;
}

void Solid_Particle::set_param(Param& param)
{
  Fluide_Incompressible::set_param(param);
  param.ajouter("e_dry", &e_dry_,Param::REQUIRED); // XD_ADD_P dry coefficient;
  param.ajouter("radius_sphere", &radius_sphere_); // XD_ADD_P radius of a spherical particle in the case of a sphere flow;
  param.ajouter_non_std("spheroid", (this)); // XD_ADD_P definition of a revolution ellipsoid;
}

int Solid_Particle::lire_motcle_non_standard(const Motcle& word, Entree& is)
{
  if (word=="spheroid")
    {
      Motcles words;
      words.add("half_small_axis");
      words.add("half_long_axis");
      Motcle secondword;
      is >> secondword;
      Motcle openbrace ="{";
      Motcle closedbrace="}";
      if (secondword==openbrace)
        {
          is >> secondword;
          while (secondword != closedbrace)
            {
              int rang2 = words.search(secondword);
              switch(rang2)
                {
                case 0:
                  is >> half_small_axis_spheroid_;
                  break;
                case 1:
                  is >> half_long_axis_spheroid_;
                  break;
                default:
                  Cerr << "Solid_Particle::lire_motcle_non_standard\n"
                       << " options of spheroid are:\n"
                       << words;
                  exit();
                }
              is >> secondword;
            }
        }
      return 1;
    }
  else
    {
      Cerr << word << " is not a keyword understood by " << que_suis_je() << " in lire_motcle_non_standard"<< finl;
      exit();
    }
  return -1;
}
