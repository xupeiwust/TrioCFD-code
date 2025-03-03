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

#ifndef SOLID_PARTICLE_included
#define SOLID_PARTICLE_included

#include <Fluide_Incompressible.h>


/*! @brief : class Solid_Particle
 *
 *  Description: This class describes spherical solid particles.
 *  This class is used in the module fluid_particle_interaction.
 *  Solid particles must be modeled with a high viscosity ratio
 *  relative to the fluid.
 *  suffix sphere :for sphere flow, therefore all particles have
 *  the same physical properties.
 *  The radius and e_dry coeff are used for the computation of the
 *  solid-solid interactions. See Modele_Collision_FT
 *  Location: ${TrioCFD_ROOT}/src/Multiphase/fluid_particle_interaction
 */
class Solid_Particle_base : public Fluide_Incompressible
{
  Declare_base_sans_constructeur(Solid_Particle_base);

public:
  Solid_Particle_base();

  // override functions
  void set_param(Param& param) override;

  // Setters
  void set_volume(const double volume) { volume_=volume; }
  void set_mass(const double mass) { mass_=mass; }
  void set_equivalent_diameter(const double equivalent_diameter)
  { equivalent_diameter_=equivalent_diameter; }
  void set_equivalent_radius(const double equivalent_radius)
  { equivalent_radius_=equivalent_radius; }
  // getters
  const double& get_e_dry() const { return e_dry_; }
  const double& get_volume() const { return volume_; }
  const double& get_mass() const { return mass_; }
  const double& get_equivalent_radius() const { return equivalent_radius_; }
  const double& get_equivalent_diameter() const { return equivalent_diameter_; }


protected :
  // Variables declaration and initialization
  double e_dry_=0.; // dry restitution coefficient
  double volume_=0.;
  double mass_=0.;
  double equivalent_radius_=0;
  double equivalent_diameter_=0;
};

#endif
