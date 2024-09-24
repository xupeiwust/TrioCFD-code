/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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
// File      : IJK_FT_cut_cell.h
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_FT_cut_cell_included
#define IJK_FT_cut_cell_included

#include <IJK_FT_base.h>
#include <IJK_Field_vector.h>
#include <Cut_cell_FT_Disc.h>
#include <IJK_FT_Post.h>
#include <Champ_diphasique.h>


/*! @brief : class IJK_FT_cut_cell
 *
 *  <Description of class IJK_FT_cut_cell>
 *
 *
 *  La classe IJK_FT_cut_cell herite de la classe IJK_FT_base.
 *
 */
class IJK_FT_cut_cell : public IJK_FT_base
{
  friend class IJK_Thermique;
  friend class IJK_Thermique_cut_cell;
  friend class Statistiques_dns_ijk_FT;
  Declare_instanciable_sans_constructeur(IJK_FT_cut_cell) ;

public :
  IJK_FT_cut_cell();
  Entree& interpreter(Entree&) override;
  void run() override;
  void euler_time_step(ArrOfDouble& var_volume_par_bulle) override;
  void rk3_sub_step(const int rk_step, const double total_timestep, const double fractionnal_timestep, const double time) override;

  void deplacer_interfaces(const double timestep,
                           const int rk_step,
                           ArrOfDouble& var_volume_par_bulle,
                           const int first_step_interface_smoothing) override;
  void deplacer_interfaces_rk3(const double timestep, const int rk_step, ArrOfDouble& var_volume_par_bulle) override;

  void update_indicator_field() override;
  void update_twice_indicator_field() override;

  Cut_cell_FT_Disc* get_cut_cell_disc() override
  {
    return &cut_cell_disc_;
  }
  const Cut_field_vector3_double& get_cut_field_velocity() const
  {
    const Cut_field_vector3_double& cut_field_velocity = static_cast<const Cut_field_vector3_double&>(velocity_);
    return cut_field_velocity;
  }
  const Cut_field_vector3_double& get_cut_field_remeshing_velocity() const
  {
    const Cut_field_vector3_double& cut_field_remeshing_velocity = static_cast<const Cut_field_vector3_double&>(remeshing_velocity_);
    return cut_field_remeshing_velocity;
  }
  const Cut_field_vector3_double& get_cut_field_total_velocity() const
  {
    const Cut_field_vector3_double& cut_field_total_velocity = static_cast<const Cut_field_vector3_double&>(total_velocity_);
    return cut_field_total_velocity;
  }

protected :
  friend class IJK_FT_Post;
  Cut_cell_FT_Disc cut_cell_disc_;

  IJK_Field_vector3_double remeshing_velocity_;
  IJK_Field_vector3_double total_velocity_;

  DoubleTabFT_cut_cell_vector3 velocity_interface_;
};

#endif /* IJK_FT_cut_cell_included */
