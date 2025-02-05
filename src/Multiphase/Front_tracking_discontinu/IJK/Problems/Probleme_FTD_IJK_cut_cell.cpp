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

#include <Probleme_FTD_IJK_cut_cell.h>
#include <Navier_Stokes_FTD_IJK.h>
#include <Cut_cell_tools.h>

Implemente_instanciable(Probleme_FTD_IJK_cut_cell, "Probleme_FTD_IJK_cut_cell", Probleme_FTD_IJK_base);

Sortie& Probleme_FTD_IJK_cut_cell::printOn(Sortie& os) const { return os; }

Entree& Probleme_FTD_IJK_cut_cell::readOn(Entree& is)
{
  Probleme_FTD_IJK_base::readOn(is);

  // Determination pour le seuil des petites cellules en cut-cell
  if ((seuil_indicatrice_petite_fixe_ == -1) && (seuil_indicatrice_petite_facsec_ == -1))
    seuil_indicatrice_petite_facsec_ = 0.125; // Default value = facsec/8. -- Note: The value facsec/20. = 0.01 would often work but not always be stable

  double seuil_indicatrice_petite;
  if (seuil_indicatrice_petite_facsec_ != -1)
    seuil_indicatrice_petite = schema_temps_ijk().get_timestep_facsec()*seuil_indicatrice_petite_facsec_;
  else
    seuil_indicatrice_petite = seuil_indicatrice_petite_fixe_;
  interfaces_.set_seuil_indicatrice_petite(seuil_indicatrice_petite);
  Cerr << "Le seuil pour l'indicatrice des petites cellules est : " << seuil_indicatrice_petite << finl;

  return is;
}

void Probleme_FTD_IJK_cut_cell::set_param(Param& param)
{
  Probleme_FTD_IJK_base::set_param(param);

  param.ajouter("seuil_indicatrice_petite_fixe", &seuil_indicatrice_petite_fixe_);
  param.ajouter("seuil_indicatrice_petite_facsec", &seuil_indicatrice_petite_facsec_);

  param.ajouter("type_surface_efficace_face", (int*)&type_surface_efficace_face_);
  param.dictionnaire("non_initialise",(int)TYPE_SURFACE_EFFICACE_FACE::NON_INITIALISE);
  param.dictionnaire("explicite",(int)TYPE_SURFACE_EFFICACE_FACE::EXPLICITE);
  param.dictionnaire("algebrique_simple",(int)TYPE_SURFACE_EFFICACE_FACE::ALGEBRIQUE_SIMPLE);
  param.dictionnaire("conservation_volume_iteratif", (int)TYPE_SURFACE_EFFICACE_FACE::CONSERVATION_VOLUME_ITERATIF);
  param.ajouter("type_surface_efficace_interface", (int*)&type_surface_efficace_interface_);
  param.dictionnaire("non_initialise",(int)TYPE_SURFACE_EFFICACE_INTERFACE::NON_INITIALISE);
  param.dictionnaire("explicite",(int)TYPE_SURFACE_EFFICACE_INTERFACE::EXPLICITE);
  param.dictionnaire("algebrique_simple",(int)TYPE_SURFACE_EFFICACE_INTERFACE::ALGEBRIQUE_SIMPLE);
  param.dictionnaire("conservation_volume", (int)TYPE_SURFACE_EFFICACE_INTERFACE::CONSERVATION_VOLUME);

  param.ajouter("facettes_interpolation", &cut_cell_facettes_interpolation_);
}

void Probleme_FTD_IJK_cut_cell::initialize()
{
  Cerr << "Probleme_FTD_IJK_cut_cell::initialize()" << finl;

  // Activation des champs cut-cell de post_ et interfaces_ (obligatoirement avant l'initialisation)
  cut_cell_disc_.initialise(interfaces_, domaine_ijk_.valeur(), Domaine_IJK::ELEM);
  post_.activate_cut_cell();
  interfaces_.activate_cut_cell();

  cut_cell_facettes_interpolation_.associer(interfaces_, cut_cell_disc_, domaine_ft_, interfaces_.maillage_ft_ijk(), interfaces_.old_maillage_ft_ijk());

  domaine_ijk_->get_local_mesh_delta(DIRECTION_K, 2 /* ghost cells */, delta_z_local_);

  thermal_probes_ghost_cells_ = 4;
  thermals_.compute_ghost_cell_numbers_for_subproblems(domaine_ijk_.valeur(), thermal_probes_ghost_cells_);
  thermal_probes_ghost_cells_ = thermals_.get_probes_ghost_cells(thermal_probes_ghost_cells_);

  // TODO : FIXME : faut boucler plus tard sur les equations IJK
  Navier_Stokes_FTD_IJK& eq_ns = ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur());
  const auto& bc = eq_ns.get_boundary_conditions();

  Cut_field_vector3_double& cut_field_velocity = static_cast<Cut_field_vector3_double&>(eq_ns.get_velocity());
  if (IJK_Shear_Periodic_helpler::defilement_ == 1)
    allocate_velocity_persistant(cut_cell_disc_, cut_field_velocity, domaine_ijk_.valeur(), 2, bc.get_dU_perio(bc.get_resolution_u_prime_()));
  else
    allocate_velocity_persistant(cut_cell_disc_, cut_field_velocity, domaine_ijk_.valeur(), thermal_probes_ghost_cells_);

  Probleme_FTD_IJK_base::initialize();
}

const Cut_field_vector3_double& Probleme_FTD_IJK_cut_cell::get_cut_field_velocity() const
{
  return static_cast<const Cut_field_vector3_double&>(ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur()).get_velocity());
}

void Probleme_FTD_IJK_cut_cell::update_indicator_field()
{
  // La suppression des cellules mortes est vraiment au tout dernier moment,
  // pour laisser la possibilite d'utiliser ces cellules lors des bilans
  // pour determiner l'indicatrice cible du remaillage.
  cut_cell_disc_.remove_dead_and_virtual_cells(interfaces_.In());

  Probleme_FTD_IJK_base::update_indicator_field();
}

void Probleme_FTD_IJK_cut_cell::update_twice_indicator_field()
{
  for(int i=0; i<2; i++)
    {
      update_indicator_field();
      update_old_intersections();
    }

  // Mise a jour des structures cut-cell
  cut_cell_disc_.update(interfaces_.I(), interfaces_.In());

  // Mise a jour des indices et coefficients des points d'interpolation a une certaine distance des facettes de l'interface
  cut_cell_perform_interpolation_facettes();

  // Calcul pour le temps old() egalement, de telle maniere a ce que les coefficients next() et old() sont initialises
  cut_cell_facettes_interpolation_.cut_cell_perform_interpolation_facettes_old(interfaces_.old());

  Cut_field_vector3_double& cut_field_velocity = static_cast<Cut_field_vector3_double&>(ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur()).get_velocity());
  cut_field_velocity[0].copie_pure_vers_diph_sans_interpolation();
  cut_field_velocity[1].copie_pure_vers_diph_sans_interpolation();
  cut_field_velocity[2].copie_pure_vers_diph_sans_interpolation();
}

void Probleme_FTD_IJK_cut_cell::deplacer_interfaces(const double timestep, const int rk_step,
                                                    ArrOfDouble& var_volume_par_bulle,
                                                    const int first_step_interface_smoothing)
{
  Cut_field_vector3_double& cut_field_velocity = static_cast<Cut_field_vector3_double&>(ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur()).get_velocity());

  thermals_.echange_diph_vers_pure_cellules_finalement_pures();
  thermals_.vide_phase_invalide_cellules_diphasiques();
  update_old_intersections(); // Pour conserver les donnees sur l'interface au temps t_{n} (en plus de t_{n+1})

  cut_field_velocity[0].copie_pure_vers_diph_sans_interpolation();
  cut_field_velocity[1].copie_pure_vers_diph_sans_interpolation();
  cut_field_velocity[2].copie_pure_vers_diph_sans_interpolation();

  Probleme_FTD_IJK_base::deplacer_interfaces(timestep, rk_step, var_volume_par_bulle, first_step_interface_smoothing);

  // Mise a jour des structures cut-cell
  cut_cell_disc_.update(interfaces_.I(), interfaces_.In());

  // Mise a jour des indices et coefficients des points d'interpolation a une certaine distance des facettes de l'interface
  cut_cell_perform_interpolation_facettes();

  thermals_.echange_pure_vers_diph_cellules_initialement_pures();
  cut_field_velocity[0].copie_pure_vers_diph_sans_interpolation();
  cut_field_velocity[1].copie_pure_vers_diph_sans_interpolation();
  cut_field_velocity[2].copie_pure_vers_diph_sans_interpolation();

  interfaces_.calcul_surface_efficace_face_initial(type_surface_efficace_face_);
  interfaces_.calcul_surface_efficace_interface_initial(type_surface_efficace_interface_);

  interfaces_.calcul_surface_efficace_face(type_surface_efficace_face_, schema_temps_ijk().get_timestep(), cut_field_velocity);
  interfaces_.calcul_surface_efficace_interface(type_surface_efficace_interface_, schema_temps_ijk().get_timestep(), cut_field_velocity);

  if (interfaces_.get_dt_impression_bilan_indicatrice() >= 0 && schema_temps_ijk().get_tstep() % interfaces_.get_dt_impression_bilan_indicatrice() == interfaces_.get_dt_impression_bilan_indicatrice() - 1)
    interfaces_.imprime_bilan_indicatrice();
}

void Probleme_FTD_IJK_cut_cell::deplacer_interfaces_rk3(const double timestep, const int rk_step,
                                                        ArrOfDouble& var_volume_par_bulle)
{
  Cut_field_vector3_double& cut_field_velocity = static_cast<Cut_field_vector3_double&>(ref_cast(Navier_Stokes_FTD_IJK, equations_.front().valeur()).get_velocity());

  thermals_.echange_diph_vers_pure_cellules_finalement_pures();
  thermals_.vide_phase_invalide_cellules_diphasiques();
  update_old_intersections(); // Pour conserver les donnees sur l'interface au temps t_{n} (en plus de t_{n+1})

  cut_field_velocity[0].copie_pure_vers_diph_sans_interpolation();
  cut_field_velocity[1].copie_pure_vers_diph_sans_interpolation();
  cut_field_velocity[2].copie_pure_vers_diph_sans_interpolation();

  Probleme_FTD_IJK_base::deplacer_interfaces_rk3(timestep, rk_step, var_volume_par_bulle);

  // Mise a jour des structures cut-cell
  cut_cell_disc_.update(interfaces_.I(), interfaces_.In());

  // Mise a jour des indices et coefficients des points d'interpolation a une certaine distance des facettes de l'interface
  cut_cell_perform_interpolation_facettes();

  thermals_.echange_pure_vers_diph_cellules_initialement_pures();
  cut_field_velocity[0].copie_pure_vers_diph_sans_interpolation();
  cut_field_velocity[1].copie_pure_vers_diph_sans_interpolation();
  cut_field_velocity[2].copie_pure_vers_diph_sans_interpolation();

  interfaces_.calcul_surface_efficace_face_initial(type_surface_efficace_face_);
  interfaces_.calcul_surface_efficace_interface_initial(type_surface_efficace_interface_);

  const double fractionnal_timestep = compute_fractionnal_timestep_rk3(timestep, ref_cast(Schema_RK3_IJK, schema_temps_ijk()).get_rk_step());

  interfaces_.calcul_surface_efficace_face(type_surface_efficace_face_, fractionnal_timestep, cut_field_velocity);
  interfaces_.calcul_surface_efficace_interface(type_surface_efficace_interface_, fractionnal_timestep, cut_field_velocity);

  if (interfaces_.get_dt_impression_bilan_indicatrice() >= 0 && schema_temps_ijk().get_tstep() % interfaces_.get_dt_impression_bilan_indicatrice() == interfaces_.get_dt_impression_bilan_indicatrice() - 1)
    interfaces_.imprime_bilan_indicatrice();
}

