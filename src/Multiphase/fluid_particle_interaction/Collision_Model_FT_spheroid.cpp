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

#include <Collision_Model_FT_spheroid.h>
#include <Probleme_FT_Disc_gen.h>


Implemente_instanciable_sans_constructeur(Collision_Model_FT_spheroid,"Collision_Model_FT_spheroid",Collision_Model_FT_base);

Collision_Model_FT_spheroid::Collision_Model_FT_spheroid()
{
}

Entree& Collision_Model_FT_spheroid::readOn (Entree& is)
{
  Collision_Model_FT_base::readOn(is);
  return is;
}

Sortie& Collision_Model_FT_spheroid::printOn(Sortie& os) const
{
  Cerr << "Error::printOn is not implemented." << finl;
  Process::exit();
  return os;
}


void Collision_Model_FT_spheroid::compute_lagrangian_contact_forces(const Fluide_Diphasique& two_phase_fluid,
                                                                    const DoubleTab& particles_position,
                                                                    const DoubleTab& particles_velocity,
                                                                    const double& deltat_simu)
{
  Process::exit("Collision_Model_FT_spheroid::compute_lagrangian_contact_forces not coded for "
                "spheroid particles");
}

void Collision_Model_FT_spheroid::discretize_contact_forces_eulerian_field(
  const DoubleTab& volumic_phase_indicator_function,
  const Domaine_VF& domain_vf,
  const IntTab& particles_eulerian_id_number,
  DoubleTab& contact_force_source_term)
{
  Process::exit("Collision_Model_FT_spheroid::discretize_contact_forces_eulerian_field "
                "not coded for spheroid particles");
}

