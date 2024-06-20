/****************************************************************************
* Copyright (c) 2024, CEA
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
// File      : IJK_One_Dimensional_Subproblems_Interfaces_Fields.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_One_Dimensional_Subproblems_Interfaces_Fields_included
#define IJK_One_Dimensional_Subproblems_Interfaces_Fields_included

#include <Objet_U.h>
#include <IJK_Lata_writer.h>
#include <Motcle.h>
#include <Intersections_Elem_Facettes_Data.h>

#define MIN_SURFACE 1.e-30
#define INVALID_TEST -1.e30
#define select(a,x,y,z) ((a==0)?(x):((a==1)?(y):(z)))

// #include <IJK_FT.h>
/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_One_Dimensional_Subproblems_Interfaces_Fields
//
// <Description of class IJK_One_Dimensional_Subproblems_Interfaces_Fields>
//
/////////////////////////////////////////////////////////////////////////////

class IJK_FT_base;
class IJK_One_Dimensional_Subproblems;
// class Intersections_Elem_Facettes;
// class IJK_Interfaces;

class IJK_One_Dimensional_Subproblems_Interfaces_Fields : public Objet_U
{

  Declare_instanciable( IJK_One_Dimensional_Subproblems_Interfaces_Fields ) ;

public :
  int initialise(const IJK_Splitting& splitting,
                 IJK_One_Dimensional_Subproblems& thermal_local_subproblems,
                 const int& debug);
  void associer(const IJK_FT_base& ijk_ft);
  void set_subproblems_interfaces_fields(const int& interface_field_type);
  void copy_previous_interface_state();
  void reset_flags();
  void posttraiter_tous_champs(Motcles& liste) const;
  int posttraiter_champs_instantanes(const Motcles& liste_post_instantanes,
                                     const char *lata_name,
                                     const int lata_step);
protected :
  REF(IJK_FT_base) ref_ijk_ft_;
  Motcles liste_post_instantanes_; // liste des champs instantanes a postraiter

  IJK_One_Dimensional_Subproblems * thermal_local_subproblems_ = nullptr;
  int interface_field_type_ = 0;
  int counter_call_ = 0;
  int debug_ = 0;

  DoubleTab ft_vertices_;
  ArrOfInt elem_crossed_;
  ArrOfDouble eulerian_values_;

  ArrOfDouble vertices_surface_;
  ArrOfDouble interfacial_heat_flux_;
  ArrOfDouble interfacial_heat_flux_sol_;
  ArrOfDouble velocity_magnitude_;
  ArrOfDouble temperature_;
  ArrOfDouble temperature_sol_;
  ArrOfDouble temperature_gradient_;
  ArrOfDouble temperature_gradient_sol_;

  IJK_Field_double tmp_field_val_;
  IJK_Field_double tmp_ft_field_val_;
  IJK_Field_double tmp_field_double_;
  IJK_Field_double tmp_ft_field_double_;
  IJK_Field_int tmp_ft_field_int_;

  /*
   * TODO : PB interface at time (n+1) only
   * Need a copy
   * There's a problem with Intersections_Elem_Facettes_Data which is a pointer
   * Need a usr-defined destructor for the current class...
   *
   * + WTF all the vertices are saved into the .lata file
   * Problem: there is no sharing between procs here
   */
  // Intersections_Elem_Facettes intersections_; // = maillage.intersections_elem_facettes();
  ArrOfDouble surface_facettes_; // = maillage.get_update_surface_facettes();
  IntTab facettes_; // = maillage.facettes();

  bool updated_interfacial_heat_flux_ = false;
  bool updated_interfacial_heat_flux_sol_ = false;
  bool updated_velocity_magnitude_ = false;
  bool updated_temperature_ = false;
  bool updated_temperature_sol_ = false;
  bool updated_temperature_gradient_ = false;
  bool updated_temperature_gradient_sol_ = false;

  bool retrieve_interfacial_surface_quantities(const Nom& surface_post_name);
  bool retrieve_interfacial_surface_quantity(ArrOfDouble& surface_quantity, bool& has_been_updated, const int& val_index);
  void averaged_by_vertex_surface(ArrOfDouble& surface_quantity);
  void get_surrounding_value(const int& i, const int& j, const int& k, double& eulerian_value);
};

#endif /* IJK_One_Dimensional_Subproblems_Interfaces_Fields_included */
