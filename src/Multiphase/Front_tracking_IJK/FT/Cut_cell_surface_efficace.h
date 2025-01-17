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
// File      : Cut_cell_surface_efficace.h
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Cut_cell_surface_efficace_included
#define Cut_cell_surface_efficace_included

#include <IJK_Field_vector.h>
#include <IJK_Field.h>
#include <Cut_field.h>

enum class TYPE_SURFACE_EFFICACE_FACE : int
{
  NON_INITIALISE, // Valeur non valide
  EXPLICITE,  // Calcul explicite de la surface efficace
  ALGEBRIQUE_SIMPLE,  // Calcul algrebrique simple de la surface efficace
  CONSERVATION_VOLUME_ITERATIF // Calcul de la surface efficace fonde sur la conservation du volume (methode iterative)
};

enum class TYPE_SURFACE_EFFICACE_INTERFACE : int
{
  NON_INITIALISE, // Valeur non valide
  EXPLICITE,  // Calcul explicite de la surface efficace
  ALGEBRIQUE_SIMPLE,  // Calcul algrebrique simple de la surface efficace
  CONSERVATION_VOLUME // Calcul de la surface efficace fonde sur la conservation du volume
};

class Cut_cell_surface_efficace
{
public:
  static void calcul_surface_interface_efficace_initiale(
    int methode_explicite,
    const IJK_Field_double& old_indicatrice_ns,
    const IJK_Field_double& next_indicatrice_ns,
    const IJK_Field_double& surface_interface_ns_old,
    const IJK_Field_double& surface_interface_ns_next,
    const IJK_Field_vector3_double& normal_of_interf_ns_old,
    const IJK_Field_vector3_double& normal_of_interf_ns_next,
    DoubleTabFT_cut_cell_vector3& normale_deplacement_interface,
    DoubleTabFT_cut_cell_scalar& surface_efficace_interface,
    DoubleTabFT_cut_cell_scalar& surface_efficace_interface_initial);

  static void calcul_vitesse_interface(
    const IJK_Field_vector3_double& velocity,
    const IJK_Field_double& old_indicatrice_ns,
    const IJK_Field_double& next_indicatrice_ns,
    const IJK_Field_vector3_double& barycentre_phase1_ns_old,
    const IJK_Field_vector3_double& barycentre_phase1_ns_next,
    DoubleTabFT_cut_cell_vector3& coord_deplacement_interface,
    DoubleTabFT_cut_cell_vector3& vitesse_deplacement_interface);

  static void calcul_surface_interface_efficace(
    double timestep,
    const IJK_Field_vector3_double& velocity,
    const IJK_Field_double& old_indicatrice_ns,
    const IJK_Field_double& next_indicatrice_ns,
    const DoubleTabFT_cut_cell_vector3& vitesse_deplacement_interface,
    const DoubleTabFT_cut_cell_vector3& normale_deplacement_interface,
    DoubleTabFT_cut_cell_scalar& surface_efficace_interface);

  // calcul_surface_face_efficace_initiale: version with DoubleTabFT_cut_cell output (use cut-cell structures)
  static void calcul_surface_face_efficace_initiale(
    int methode_explicite,
    const IJK_Field_vector3_double& old_indicatrice_surfacique_face_ns,
    const IJK_Field_vector3_double& next_indicatrice_surfacique_face_ns,
    DoubleTabFT_cut_cell_vector3& indicatrice_surfacique_efficace_face,
    DoubleTabFT_cut_cell_vector3& indicatrice_surfacique_efficace_face_initial);

  // calcul_surface_face_efficace_initiale: version with IJK_Field output (does not use cut-cell structures)
  static void calcul_surface_face_efficace_initiale(
    int methode_explicite,
    const IJK_Field_vector3_double& old_indicatrice_surfacique_face_ns,
    const IJK_Field_vector3_double& next_indicatrice_surfacique_face_ns,
    IJK_Field_vector3_double& indicatrice_surfacique_efficace_face,
    IJK_Field_vector3_double& indicatrice_surfacique_efficace_face_initial);

  static void calcul_surface_face_efficace_iteratif(
    int verbosite_surface_efficace_face,
    double timestep,
    const Cut_field_vector3_double& velocity,
    int& iteration_solver_surface_efficace_face,
    const IJK_Field_double& old_indicatrice_ns,
    const IJK_Field_double& next_indicatrice_ns,
    const IJK_Field_vector3_double& old_indicatrice_surfacique_face_ns,
    const IJK_Field_vector3_double& next_indicatrice_surfacique_face_ns,
    DoubleTabFT_cut_cell_vector3& indicatrice_surfacique_efficace_face,
    const DoubleTabFT_cut_cell_vector3& indicatrice_surfacique_efficace_face_initial,
    DoubleTabFT_cut_cell_vector6& indicatrice_surfacique_efficace_face_correction,
    DoubleTabFT_cut_cell_scalar& indicatrice_surfacique_efficace_face_absolute_error);


  static void imprimer_informations_surface_efficace_interface(
    int verbosite_surface_efficace_interface,
    double timestep,
    const IJK_Field_vector3_double& velocity,
    const IJK_Field_double& old_indicatrice_ns,
    const IJK_Field_double& next_indicatrice_ns,
    const DoubleTabFT_cut_cell_scalar& surface_efficace_interface,
    const DoubleTabFT_cut_cell_scalar& surface_efficace_interface_initial,
    const DoubleTabFT_cut_cell_vector3& normale_deplacement_interface,
    const DoubleTabFT_cut_cell_vector3& vitesse_deplacement_interface);

  static void imprimer_informations_surface_efficace_face(
    int verbosite_surface_efficace_face,
    int iteration_solver_surface_efficace_face,
    double timestep,
    const Cut_field_vector3_double& velocity,
    const IJK_Field_double& old_indicatrice_ns,
    const IJK_Field_double& next_indicatrice_ns,
    const DoubleTabFT_cut_cell_vector3& indicatrice_surfacique_efficace_face,
    const DoubleTabFT_cut_cell_vector3& indicatrice_surfacique_efficace_face_initial);

  static void calcul_vitesse_remaillage(double timestep,
                                        const IJK_Field_double& indicatrice_avant_remaillage,
                                        const IJK_Field_double& indicatrice_apres_remaillage,
                                        const IJK_Field_double& indicatrice_fin_pas_de_temps,
                                        DoubleTabFT_cut_cell_vector3& indicatrice_surfacique_efficace_remaillage_face,
                                        Cut_field_vector3_double& remeshing_velocity);

  static void calcul_delta_volume_theorique_bilan(int compo, const DoubleTab& bounding_box_bulles, double timestep,
                                                  const IJK_Field_double& indicatrice_avant_deformation,
                                                  const IJK_Field_double& indicatrice_apres_deformation,
                                                  const IJK_Field_vector3_double& indicatrice_surfacique_efficace_deformation_face,
                                                  const Cut_field_vector3_double& deformation_velocity,
                                                  IJK_Field_double& delta_volume_theorique_bilan);
};

#endif /* Cut_cell_surface_efficace_included */
