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

#ifndef IJK_Composantes_Connex_included
#define IJK_Composantes_Connex_included

#include <Objet_U.h>
#include <IJK_Field.h>

#define DIRECTION_I 0
#define DIRECTION_J 1
#define DIRECTION_K 2
#define LIQUID_INDICATOR_TEST 1.-1.e-12
#define VAPOUR_INDICATOR_TEST 1.e-12
#define NEIGHBOURS_I {-1, 1, 0, 0, 0, 0}
#define NEIGHBOURS_J {0, 0, -1, 1, 0, 0}
#define NEIGHBOURS_K {0, 0, 0, 0, -1, 1}

class Probleme_FTD_IJK_base;
class IJK_Interfaces;
class IJK_Composantes_Connex : public Objet_U
{

  Declare_instanciable( IJK_Composantes_Connex ) ;

public :
  void initialize(IJK_Interfaces& interfaces,
                  const bool is_switch);
  void allocate_fields(const Domaine_IJK& splitting,
                       const int& compute_compo_fields);
  void associer(const Probleme_FTD_IJK_base& ijk_ft);
  void initialise_bubbles_params();
  void associate_rising_velocities_parameters(const Domaine_IJK& splitting,
                                              const int& compute_rising_velocities,
                                              const int& fill_rising_velocities,
                                              const int& use_bubbles_velocities_from_interface,
                                              const int& use_bubbles_velocities_from_barycentres);
  void compute_bounding_box_fill_compo_connex();
  void compute_compo_connex_from_interface();
  void compute_rising_velocities();

  const IJK_Field_double& get_eulerian_compo_connex_ft() const {  return eulerian_compo_connex_ft_;  }
  const IJK_Field_double& get_eulerian_compo_connex() const  {    return eulerian_compo_connex_ns_;   }
  const IJK_Field_double& get_eulerian_compo_connex_ghost_ft() const   {     return eulerian_compo_connex_ghost_ft_;   }
  const IJK_Field_double& get_eulerian_compo_connex_ghost() const    {     return eulerian_compo_connex_ghost_ns_;   }
  const IJK_Field_double& get_eulerian_compo_connex_from_interface_ft() const  {     return eulerian_compo_connex_from_interface_ft_;   }
  const IJK_Field_double& get_eulerian_compo_connex_from_interface_ghost_ft() const   {     return eulerian_compo_connex_from_interface_ghost_ft_;   }
  const IJK_Field_double& get_eulerian_compo_connex_from_interface_ns() const   {     return eulerian_compo_connex_from_interface_ns_;   }
  const IJK_Field_double& get_eulerian_compo_connex_from_interface_ghost_ns() const   {     return eulerian_compo_connex_from_interface_ghost_ns_;   }
  const IJK_Field_int& get_eulerian_compo_connex_int_from_interface_ns() const   {     return eulerian_compo_connex_from_interface_int_ns_;   }
  const IJK_Field_int& get_eulerian_compo_connex_int_from_interface_ghost_ns() const   {     return eulerian_compo_connex_from_interface_ghost_int_ns_;   }
  const DoubleTab& get_bounding_box() const   {     return bounding_box_;   }
  const DoubleTab& get_bubbles_barycentre() const   {     return bubbles_barycentre_;   }
  const ArrOfDouble& get_bubbles_volume() const   {     return bubbles_volume_;   }
  const IJK_Field_double& get_eulerian_rising_velocities() const   {     return eulerian_rising_velocities_;   }
  const ArrOfDouble& get_rising_velocities() const   {     return rising_velocities_;   }
  const DoubleTab& get_rising_vectors() const   {     return rising_vectors_;   }
  const Vecteur3& get_rising_velocity_overall() const   {     return rising_velocity_overall_;   }
  const Vecteur3& get_liquid_velocity() const   {     return liquid_velocity_;   }
  const DoubleTab& get_min_max_larger_box() const   {     return min_max_larger_box_;   }
  const int& get_compute_from_bounding_box() const   {     return compute_from_bounding_box_;   }
  const int& get_compute_compo_fields() const   {     return compute_compo_fields_;   }

protected :
  void fill_mixed_cell_compo();
  OBS_PTR(Probleme_FTD_IJK_base) ref_ijk_ft_;
  IJK_Interfaces * interfaces_ = nullptr;
  bool is_switch_=false;
  int compute_compo_fields_=0;
  int compute_from_bounding_box_=0;

  IJK_Field_double eulerian_compo_connex_ft_;
  IJK_Field_double eulerian_compo_connex_ns_;
  IJK_Field_double eulerian_compo_connex_ghost_ft_;
  IJK_Field_double eulerian_compo_connex_ghost_ns_;

  /*
   * TODO: write redistribute for IJK_Field_int
   */
  IJK_Field_double eulerian_compo_connex_from_interface_ft_;
  IJK_Field_double eulerian_compo_connex_from_interface_ns_;
  IJK_Field_double eulerian_compo_connex_from_interface_ghost_ft_;
  IJK_Field_double eulerian_compo_connex_from_interface_ghost_ns_;
  IJK_Field_int eulerian_compo_connex_from_interface_int_ns_;
  IJK_Field_int eulerian_compo_connex_from_interface_ghost_int_ns_;
  IJK_Field_int eulerian_compo_connex_valid_compo_field_;

  DoubleTab bounding_box_;
  DoubleTab bubbles_barycentre_;
  ArrOfDouble bubbles_volume_;
  DoubleTab min_max_larger_box_;

  IJK_Field_double eulerian_rising_velocities_;
  ArrOfDouble rising_velocities_;
  Vecteur3 rising_velocity_overall_ = {0., 0., 0.};
  DoubleTab rising_vectors_;
  Vecteur3 liquid_velocity_ = {0., 0., 0.};

  int compute_rising_velocities_ = 0;
  int fill_rising_velocities_ = 0;

  int use_bubbles_velocities_from_interface_ = 0;
  int use_bubbles_velocities_from_barycentres_ = 0;

  bool is_updated_ = false;
};

#endif /* IJK_Composantes_Connex_included */
