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
// File      : Cut_cell_schema_auxiliaire.h
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Cut_cell_schema_auxiliaire_included
#define Cut_cell_schema_auxiliaire_included

#include <FixedVector.h>
#include <IJK_Field.h>
#include <Param.h>
#include <Cut_field.h>
#include <Cut_cell_correction_petites_cellules.h>
#include <Maillage_FT_IJK.h>
#include <Objet_U.h>

class Probleme_FTD_IJK_cut_cell;

enum class METHODE_TEMPERATURE_REMPLISSAGE : int
{
  NON_INITIALISE,     // Valeur invalide par defaut, pour forcer le choix
  COPIE_DIRECTE,      // Copie de la valeur actuelle de la temperature.
  PONDERATION_VOISIN, // Moyenne ponderee des voisins pour estimer la temperature
  PONDERATION_DIRECTIONNELLE_VOISIN, // Moyenne ponderee et directionnelle des voisins pour estimer la temperature
  SEMI_LAGRANGIEN,    // Approximation semi-lagrangienne du deplacement pour estimer la temperature
  SEMI_LAGRANGIEN_INTERPOLATE // Approximation semi-lagrangienne du deplacement pour estimer la temperature
};


class Cut_cell_schema_auxiliaire : public Objet_U
{
  Declare_base(Cut_cell_schema_auxiliaire);

public:
  void initialise(Cut_cell_FT_Disc& cut_cell_disc);
  void set_param(Param& param);

  virtual double dying_cells_flux(int num_face, int phase, int n, const Cut_field_vector3_double& cut_field_total_velocity, const Cut_field_double& cut_field_temperature) = 0;
  virtual double small_nascent_cells_flux(int num_face, int phase, int n, const Cut_field_vector3_double& cut_field_total_velocity, const Cut_field_double& cut_field_temperature) = 0;

  void compute_flux_dying_cells(const Cut_field_vector3_double& cut_field_total_velocity, const Cut_field_double& cut_field_temperature);

  void compute_flux_small_nascent_cells(const Cut_field_vector3_double& cut_field_total_velocity, Cut_field_double& cut_field_temperature);

  void calcule_valeur_remplissage(double timestep, double lambda_liquid, double lambda_vapour, const IJK_Field_double& flux_interface_ns, const ArrOfDouble& interfacial_temperature, const IJK_Field_double& temperature_ft, const Cut_field_vector3_double& cut_field_total_velocity, const Cut_field_double& cut_field_temperature);

  void add_dying_cells(const Cut_field_vector3_double& cut_field_total_velocity, Cut_field_double& cut_field_temperature, bool write_flux, Cut_field_vector3_double& cut_field_current_fluxes);
  void add_small_nascent_cells(const Cut_field_vector3_double& cut_field_total_velocity, Cut_field_double& cut_field_temperature, bool write_flux, Cut_field_vector3_double& cut_field_current_fluxes);

protected:
  void calcule_valeur_remplissage_copie_directe(const Cut_field_double& cut_field_temperature, DoubleTabFT_cut_cell_scalar& valeur_remplissage);
  void calcule_valeur_remplissage_ponderation_voisin(bool est_directionnel, const Cut_field_vector3_double& cut_field_total_velocity, const Cut_field_double& cut_field_temperature, DoubleTabFT_cut_cell_scalar& valeur_remplissage);
  void calcule_valeur_remplissage_semi_lagrangien(double timestep, double lambda_liquid, double lambda_vapour, const IJK_Field_double& flux_interface_ns, const Cut_field_double& cut_field_temperature, DoubleTabFT_cut_cell_scalar& valeur_remplissage);
  void calcule_valeur_remplissage_semi_lagrangien_interpolate(double timestep, const ArrOfDouble& interfacial_temperature, const IJK_Field_double& temperature_ft, const Cut_field_double& cut_field_temperature, DoubleTabFT_cut_cell_scalar& valeur_remplissage);

  CORRECTION_PETITES_CELLULES correction_petites_cellules_;

  DoubleTabFT_cut_cell_vector6 flux_naive_;
  DoubleTabFT_cut_cell_scalar temperature_remplissage_;
  METHODE_TEMPERATURE_REMPLISSAGE methode_valeur_remplissage_;

  bool no_static_update_; // Disable the correction if there is no velocity

  int tolerate_not_within_tetrahedron_; // Parameter to disable crashes if no suitable tetrahedron is found in the cut-cell interpolation (for debug purpose normally)
};

#endif /* Cut_cell_schema_auxiliaire_included */
