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

#ifndef IJK_subdivided_mesh_tools_included
#define IJK_subdivided_mesh_tools_included

#include <IJK_Field.h>
#include <IJK_Field_vector.h>
#include <Champ_diphasique.h>

// Ce ficher regroupe differentes fonctions utiles pour la creation d'un maillage
// deux fois plus fin que le maillage de la simulation dans chaque direction, et
// le regroupement de l'information sur ce maillage pour creer des champs sur le
// maillage de depart.
// Par exemple, cela permet si connait l'indicatrice sur un maillage deux fois plus fin,
// de deduire l'indicatrice sur le maillage de depart, sur le volume de controle des
// elements ou bien sur les volumes de controles decales aux centres des faces.

void build_subdivided_splitting(const IJK_Splitting& split_ext, IJK_Splitting& split_8x);

void extend_ft_to_8x(const IJK_Field_double& field_ft, IJK_Field_double& field_8x);
void extend_ft_to_8x(const IJK_Field_vector3_double& field_ft, IJK_Field_vector3_double& field_8x);

void moyenne_sous_elements_elem(const IJK_Field_double& field_8x, IJK_Field_double& field_ft);
void moyenne_indicatrice_sous_elements_elem(const IJK_Field_double& field_8x, IJK_Field_double& field_ft, double tolerance);
void moyenne_indicatrice_sous_elements_face(int dir, const IJK_Field_double& field_8x, IJK_Field_double& field_ft, double tolerance);
void somme_sous_elements_si_indicatrice_elem(const IJK_Field_double& field_8x, IJK_Field_double& field_ft, const IJK_Field_double& indicatrice_ft);
void somme_sous_elements_si_indicatrice_face(int dir, const IJK_Field_double& field_8x, IJK_Field_double& field_ft, const IJK_Field_double& indicatrice_ft);

void moyenne_indicatrice_surfacique_sous_faces_elem(int face_dir, const IJK_Field_double& field_8x, IJK_Field_double& field_ft, const IJK_Field_double& indicatrice_ft);
void moyenne_indicatrice_surfacique_sous_faces_face(int dir, int face_dir, const IJK_Field_double& field_8x, IJK_Field_double& field_ft, const IJK_Field_double& indicatrice_ft);

void moyenne_barycentre_phase1_sous_elements_elem(int bary_compo, const IJK_Field_double& indic_8x, const IJK_Field_double& field_8x, IJK_Field_double& field_ft, const IJK_Field_double& indicatrice_ft);
void moyenne_barycentre_phase1_sous_elements_face(int dir, int bary_compo, const IJK_Field_double& indic_8x, const IJK_Field_double& field_8x, IJK_Field_double& field_ft, const IJK_Field_double& indicatrice_ft);
void moyenne_barycentre_phase1_sous_faces_elem(int face_dir, int bary_compo, const IJK_Field_double& indic_8x, const IJK_Field_double& field_8x, IJK_Field_double& field_ft, const IJK_Field_double& indicatrice_ft);
void moyenne_barycentre_phase1_sous_faces_face(int dir, int face_dir, int bary_compo, const IJK_Field_double& indic_8x, const IJK_Field_double& field_8x, IJK_Field_double& field_ft, const IJK_Field_double& indicatrice_ft);

#endif /* IJK_subdivided_mesh_tools_included */
