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

#ifndef Modele_Collision_FT_sphere_included
#define Modele_Collision_FT_sphere_included

#include <Collision_Model_FT_base.h>

/*! @brief : class Collision_Model_FT
 *
 *  Description: This class enables to compute solid-solid
 *  interactions for fpi module under the framework of
 *  soft-sphere collision model. Under this framework,
 *  multiple collisions can occurs at the same time (ie a
 *  particle can collide with 2 or more particles). The
 *  collision is spread out on multiple time steps. A slight
 *  overlap (less than the mesh grid size) occurs during the
 *  process.
 */

class Collision_Model_FT_sphere : public Collision_Model_FT_base
{
  Declare_instanciable_sans_constructeur(Collision_Model_FT_sphere);

public:

  Collision_Model_FT_sphere();
  void compute_lagrangian_contact_forces(const Fluide_Diphasique& two_phase_fluid,
                                         const DoubleTab& particles_position,
                                         const DoubleTab& particles_velocity,
                                         const double& deltat_simu) override;

  void discretize_contact_forces_eulerian_field(const DoubleTab& volumic_phase_indicator_function,
                                                const Domaine_VF& domain_vf,
                                                const IntTab& particles_eulerian_id_number,
                                                DoubleTab& contact_force_source_term) override;

};

#endif

