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

#include <Probleme_FTD_IJK.h>

Implemente_instanciable(Probleme_FTD_IJK, "Probleme_FTD_IJK", Probleme_FTD_IJK_base);

Sortie& Probleme_FTD_IJK::printOn(Sortie& os) const { return os; }
Entree& Probleme_FTD_IJK::readOn(Entree& is) { return Probleme_FTD_IJK_base::readOn(is); }

void Probleme_FTD_IJK::initialize()
{
  Cerr << "Probleme_FTD_IJK::initialize()" << finl;

  domaine_ijk_->get_local_mesh_delta(DIRECTION_K, 2 /* ghost cells */, delta_z_local_);

  thermal_probes_ghost_cells_ = 4;
  thermals_.compute_ghost_cell_numbers_for_subproblems(domaine_ijk_.valeur(), thermal_probes_ghost_cells_);
  thermal_probes_ghost_cells_ = thermals_.get_probes_ghost_cells(thermal_probes_ghost_cells_);

  if (IJK_Shear_Periodic_helpler::defilement_ == 1)
    allocate_velocity(velocity_, domaine_ijk_.valeur(), 2, boundary_conditions_.get_dU_perio(boundary_conditions_.get_resolution_u_prime_()));
  else
    allocate_velocity(velocity_, domaine_ijk_.valeur(), thermal_probes_ghost_cells_);

  Probleme_FTD_IJK_base::initialize();
}

// Pour creer une dilatation forcee (cas test champs FT)
void Probleme_FTD_IJK::create_forced_dilation()
{
  if (coeff_evol_volume_ != 0.)
    {
      vol_bulle_monodisperse_ = vol_bulle_monodisperse_ * (1. + coeff_evol_volume_ * schema_temps_ijk().get_timestep());
      const int nb_reelles = interfaces_.get_nb_bulles_reelles();
      for (int ib = 0; ib < nb_reelles; ib++)
        vol_bulles_[ib] = vol_bulle_monodisperse_;
    }
}
