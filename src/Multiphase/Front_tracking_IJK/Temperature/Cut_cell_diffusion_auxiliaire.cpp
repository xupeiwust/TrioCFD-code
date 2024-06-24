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
// File      : Cut_cell_diffusion_auxiliaire.cpp
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#include <Cut_cell_diffusion_auxiliaire.h>
#include <Cut_cell_FT_Disc.h>
#include <IJK_Thermal_base.h>
#include <IJK_FT_cut_cell.h>
#include <Process.h>
#include <stat_counters.h>
#include <IJK_Navier_Stokes_tools.h>
#include <IJK_Navier_Stokes_tools_cut_cell.h>
#include <Param.h>

Implemente_instanciable_sans_constructeur(Cut_cell_diffusion_auxiliaire, "Cut_cell_diffusion_auxiliaire", Cut_cell_schema_auxiliaire) ;

Cut_cell_diffusion_auxiliaire::Cut_cell_diffusion_auxiliaire()
{
  methode_temperature_remplissage_ = METHODE_TEMPERATURE_REMPLISSAGE::COPIE_DIRECTE;

  second_order_diffusion_interface_ = 0;
  methode_flux_interface_ = METHODE_FLUX_INTERFACE::INTERP_CUT_CELL;
  scaled_distance_flux_interface_ = 1.0;
  scaled_distance_second_point_flux_interface_ = 2.0;

  deactivate_correction_petites_cellules_diffusion_ = 0;
  correction_petites_cellules_ = CORRECTION_PETITES_CELLULES::DIRECTION_PRIVILEGIEE_AVEC_LIMITATION_2;

  no_static_update_ = false;
}

Sortie& Cut_cell_diffusion_auxiliaire::printOn(Sortie& os) const
{
  Objet_U::printOn(os);
  return os;
}

Entree& Cut_cell_diffusion_auxiliaire::readOn(Entree& is)
{
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades(is);
  return is;
}

void Cut_cell_diffusion_auxiliaire::set_param(Param& param)
{
  Cut_cell_schema_auxiliaire::set_param(param);

  param.ajouter_flag("second_order_diffusion_interface", &second_order_diffusion_interface_);
  param.ajouter("methode_flux_interface", (int*)&methode_flux_interface_);
  param.dictionnaire("non_initialise", (int)METHODE_FLUX_INTERFACE::NON_INITIALISE);
  param.dictionnaire("interp_pure", (int)METHODE_FLUX_INTERFACE::INTERP_PURE);
  param.dictionnaire("interp_cut_cell", (int)METHODE_FLUX_INTERFACE::INTERP_CUT_CELL);
  param.dictionnaire("local_cellule", (int)METHODE_FLUX_INTERFACE::LOCAL_CELLULE);
  param.ajouter("scaled_distance_flux_interface", &scaled_distance_flux_interface_);
  param.ajouter("scaled_distance_second_point_flux_interface", &scaled_distance_second_point_flux_interface_);

  param.ajouter_flag("deactivate_correction_petites_cellules_diffusion", &deactivate_correction_petites_cellules_diffusion_);
}

void Cut_cell_diffusion_auxiliaire::calculer_flux_interface(bool next_time, double lambda_liquid, double lambda_vapour, Facettes_data& coord_facettes, Facettes_data& interfacial_temperature, DoubleTabFT& interfacial_phin_ai, const Cut_field_double& cut_field_temperature, REF(IJK_FT_cut_cell)& ref_ijk_ft, const IJK_Field_double& temperature_ns, IJK_Field_double& temperature_ft)
{
  // Transfer the field to FT splitting, needed for compute_interfacial_temperature
  ref_ijk_ft->redistrib_to_ft_elem().redistribute(temperature_ns, temperature_ft);
  temperature_ft.echange_espace_virtuel(temperature_ft.ghost());

  if (methode_flux_interface_ == METHODE_FLUX_INTERFACE::INTERP_PURE)
    {
      compute_interfacial_temperature2(next_time, lambda_liquid, lambda_vapour, temperature_ft, ref_ijk_ft->get_geometry(), ref_ijk_ft->itfce().maillage_ft_ijk(), coord_facettes, interfacial_temperature, interfacial_phin_ai);
    }
  else if (methode_flux_interface_ == METHODE_FLUX_INTERFACE::INTERP_CUT_CELL)
    {
      compute_interfacial_temperature_cut_cell(next_time, lambda_liquid, lambda_vapour, cut_field_temperature, ref_ijk_ft->get_geometry(), ref_ijk_ft->itfce().maillage_ft_ijk(), coord_facettes, interfacial_temperature, interfacial_phin_ai);
    }
  else if (methode_flux_interface_ == METHODE_FLUX_INTERFACE::LOCAL_CELLULE)
    {
      assert(flux_interface_ft_[next_time].get_splitting().get_nb_elem_local(0) == flux_interface_ft_[next_time].ni());
      assert(flux_interface_ft_[next_time].get_splitting().get_nb_elem_local(1) == flux_interface_ft_[next_time].nj());
      assert(flux_interface_ft_[next_time].get_splitting().get_nb_elem_local(2) == flux_interface_ft_[next_time].nk());
      compute_interfacial_temperature_local_normal(next_time, lambda_liquid, lambda_vapour, cut_field_temperature, ref_ijk_ft->itfce().maillage_ft_ijk(), flux_interface_ft_[next_time].get_splitting(), ref_ijk_ft->itfce(), coord_facettes, interfacial_temperature, interfacial_phin_ai);
    }
  else
    {
      Cerr << "Methode non reconnue pour le calcul du flux a l'interface." << finl;
      Process::exit();
    }

  // Calcul des flux sur le maillage FT
  {
    const Cut_cell_FT_Disc& cut_cell_disc = flux_interface_efficace_.get_cut_cell_disc();
    const Maillage_FT_IJK& mesh = cut_cell_disc.get_interfaces().maillage_ft_ijk();
    const Intersections_Elem_Facettes& intersec = next_time ? mesh.intersections_elem_facettes() : mesh.intersections_elem_facettes_old();
    const IJK_Splitting& s = flux_interface_ft_[next_time].get_splitting();

    const int ni = flux_interface_ft_[next_time].ni();
    const int nj = flux_interface_ft_[next_time].nj();
    const int nk = flux_interface_ft_[next_time].nk();

    // Initialisation
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  flux_interface_ft_[next_time](i, j, k) = 0.;
                }
            }
        }
    }

    // Calcul pour les faces coupees par l'interface
    {
      const ArrOfInt& index_elem = intersec.index_elem();
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  if (next_time ? (!cut_cell_disc.get_interfaces().est_pure(cut_cell_disc.get_interfaces().In_ft(i,j,k))) : (!cut_cell_disc.get_interfaces().est_pure(cut_cell_disc.get_interfaces().I_ft(i,j,k))))
                    {
                      assert(mesh.ref_splitting().valeur() == s);
                      const int num_elem = s.convert_ijk_cell_to_packed(i,j,k);
                      int index = index_elem[num_elem];
                      double somme_contrib = 0.;
                      // Boucle sur les facettes qui traversent cet element
                      while (index >= 0)
                        {
                          const Intersections_Elem_Facettes_Data& data = intersec.data_intersection(index);
                          const int fa7 = data.numero_facette_;
                          somme_contrib += interfacial_phin_ai[fa7] * data.fraction_surface_intersection_;

                          index = data.index_facette_suivante_;
                        };

                      flux_interface_ft_[next_time](i,j,k) = somme_contrib;
                    }
                }
            }
        }
    }
  }

  flux_interface_ft_[next_time].echange_espace_virtuel(flux_interface_ft_[next_time].ghost());

  // Calcul du flux interface sur le domaine NS :
  ref_ijk_ft->redistrib_from_ft_elem().redistribute(flux_interface_ft_[next_time], flux_interface_ns_[next_time]);
  flux_interface_ns_[next_time].echange_espace_virtuel(flux_interface_ns_[next_time].ghost());
}

void Cut_cell_diffusion_auxiliaire::calculer_flux_interface_efficace()
{
  // Correction par la surface efficace
  const Cut_cell_FT_Disc& cut_cell_disc = flux_interface_efficace_.get_cut_cell_disc();
  const DoubleTabFT_cut_cell_scalar& surface_efficace_interface = cut_cell_disc.get_interfaces().get_surface_efficace_interface();
  const IJK_Field_double& old_surface_interface = cut_cell_disc.get_interfaces().get_surface_interface_old();
  const IJK_Field_double& next_surface_interface = cut_cell_disc.get_interfaces().get_surface_interface_next();
  for (int n = 0; n < cut_cell_disc.get_n_tot(); n++)
    {
      Int3 ijk = cut_cell_disc.get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      assert(flux_interface_ns_[0].ghost() == flux_interface_ns_[1].ghost());
      if (!cut_cell_disc.within_ghost(i, j, k, flux_interface_ns_[0].ghost(), flux_interface_ns_[0].ghost()))
        continue;

      double old_indicatrice = cut_cell_disc.get_interfaces().I(i,j,k);
      double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);

      // Note : contrairement au champs de IJK_Interfaces, on a toujours [0] = old et [1] = next pour le flux_interface
      if (IJK_Interfaces::devient_pure(old_indicatrice, next_indicatrice))
        {
          assert(old_surface_interface(i,j,k) != 0.);
          flux_interface_efficace_(n) = (flux_interface_ns_[0](i,j,k)/old_surface_interface(i,j,k)) * surface_efficace_interface(n);
        }
      else if (IJK_Interfaces::devient_diphasique(old_indicatrice, next_indicatrice))
        {
          assert(next_surface_interface(i,j,k) != 0.);
          flux_interface_efficace_(n) = (flux_interface_ns_[1](i,j,k)/next_surface_interface(i,j,k)) * surface_efficace_interface(n);
        }
      else
        {
          assert(old_surface_interface(i,j,k) != 0.);
          assert(next_surface_interface(i,j,k) != 0.);
          flux_interface_efficace_(n) = .5 * (flux_interface_ns_[0](i,j,k)/old_surface_interface(i,j,k) + flux_interface_ns_[1](i,j,k)/next_surface_interface(i,j,k)) * surface_efficace_interface(n);
        }
    }
}

void Cut_cell_diffusion_auxiliaire::calculer_flux_interface_next(double lambda_liquid, double lambda_vapour, Facettes_data& coord_facettes, Facettes_data& interfacial_temperature, DoubleTabFT& interfacial_phin_ai, const Cut_field_double& cut_field_temperature, REF(IJK_FT_cut_cell)& ref_ijk_ft, const IJK_Field_double& temperature_ns, IJK_Field_double& temperature_ft)
{
  calculer_flux_interface(true, lambda_liquid, lambda_vapour, coord_facettes, interfacial_temperature, interfacial_phin_ai, cut_field_temperature, ref_ijk_ft, temperature_ns, temperature_ft);
}

void Cut_cell_diffusion_auxiliaire::calculer_flux_interface_old(double lambda_liquid, double lambda_vapour, Facettes_data& coord_facettes, Facettes_data& interfacial_temperature, DoubleTabFT& interfacial_phin_ai, const Cut_field_double& cut_field_temperature, REF(IJK_FT_cut_cell)& ref_ijk_ft, const IJK_Field_double& temperature_ns, IJK_Field_double& temperature_ft)
{
  calculer_flux_interface(false, lambda_liquid, lambda_vapour, coord_facettes, interfacial_temperature, interfacial_phin_ai, cut_field_temperature, ref_ijk_ft, temperature_ns, temperature_ft);
}

// Copie de IJK_Thermal_base::compute_interfacial_temperature2
void Cut_cell_diffusion_auxiliaire::compute_interfacial_temperature2(bool next_time,
                                                                     double lambda_liquid,
                                                                     double lambda_vapour,
                                                                     const IJK_Field_double& temperature_ft,
                                                                     const IJK_Grid_Geometry& geom,
                                                                     const Maillage_FT_IJK& maillage,
                                                                     Facettes_data& coord_facettes,
                                                                     Facettes_data& interfacial_temperature,
                                                                     DoubleTabFT& flux_normal_interp)
{
  const double dist_1 = scaled_distance_flux_interface_/sqrt(3) * std::pow(std::pow(geom.get_constant_delta(0), 2.) +
                                                                           std::pow(geom.get_constant_delta(1), 2.) +
                                                                           std::pow(geom.get_constant_delta(2), 2.),
                                                                           0.5);
  const double dist_2 = scaled_distance_second_point_flux_interface_/sqrt(3) * std::pow(std::pow(geom.get_constant_delta(0), 2.) +
                                                                                        std::pow(geom.get_constant_delta(1), 2.) +
                                                                                        std::pow(geom.get_constant_delta(2), 2.),
                                                                                        0.5);
  const DoubleTab& normale_facettes = next_time ? maillage.get_update_normale_facettes() : maillage.get_normale_facettes_old();
  const ArrOfDouble& surface_facettes = next_time ? maillage.get_update_surface_facettes() : maillage.get_surface_facettes_old();
  const int nb_facettes = next_time ? maillage.nb_facettes() : maillage.nb_facettes_old();
  const IntTab& facettes = next_time ? maillage.facettes() : maillage.facettes_old();
  const DoubleTab& sommets = next_time ? maillage.sommets() : maillage.sommets_old();

  coord_facettes.centre.resize(nb_facettes, 3);
  coord_facettes.centre = 0.;

  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    for (int dir = 0; dir < 3; dir++)
      for (int som = 0; som < 3; som++)
        coord_facettes.centre(fa7, dir) += sommets(facettes(fa7, som), dir);

  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    for (int dir = 0; dir < 3; dir++)
      coord_facettes.centre(fa7, dir) /= 3.;

  if (second_order_diffusion_interface_)
    {
      calcul_temperature_flux_interface_second_order(temperature_ft,
                                                     lambda_liquid,
                                                     lambda_vapour,
                                                     dist_1,
                                                     dist_2,
                                                     coord_facettes.centre,
                                                     normale_facettes,
                                                     interfacial_temperature.centre,
                                                     flux_normal_interp,
                                                     interfacial_temperature.liqu_1,
                                                     interfacial_temperature.vap_1,
                                                     interfacial_temperature.liqu_2,
                                                     interfacial_temperature.vap_2,
                                                     coord_facettes.liqu_1,
                                                     coord_facettes.vap_1,
                                                     coord_facettes.liqu_2,
                                                     coord_facettes.vap_2);
    }
  else
    {
      calcul_temperature_flux_interface(temperature_ft,
                                        lambda_liquid,
                                        lambda_vapour,
                                        dist_1,
                                        coord_facettes.centre,
                                        normale_facettes,
                                        interfacial_temperature.centre,
                                        flux_normal_interp,
                                        interfacial_temperature.liqu_1,
                                        interfacial_temperature.vap_1,
                                        coord_facettes.liqu_1,
                                        coord_facettes.vap_1);
    }

  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    {
      flux_normal_interp(fa7) *= surface_facettes(fa7);
      interfacial_temperature.centre(fa7) *= surface_facettes(fa7);
      interfacial_temperature.liqu_1(fa7) *= surface_facettes(fa7);
      interfacial_temperature.vap_1(fa7) *= surface_facettes(fa7);
    }
}

void Cut_cell_diffusion_auxiliaire::compute_interfacial_temperature_cut_cell(bool next_time,
                                                                             double lambda_liquid,
                                                                             double lambda_vapour,
                                                                             const Cut_field_double& cut_field_temperature,
                                                                             const IJK_Grid_Geometry& geom,
                                                                             const Maillage_FT_IJK& maillage,
                                                                             Facettes_data& coord_facettes,
                                                                             Facettes_data& interfacial_temperature,
                                                                             DoubleTabFT& flux_normal_interp)
{
  const double dist_1 = scaled_distance_flux_interface_/sqrt(3) * std::pow(std::pow(geom.get_constant_delta(0), 2.) +
                                                                           std::pow(geom.get_constant_delta(1), 2.) +
                                                                           std::pow(geom.get_constant_delta(2), 2.),
                                                                           0.5);
  const double dist_2 = scaled_distance_second_point_flux_interface_/sqrt(3) * std::pow(std::pow(geom.get_constant_delta(0), 2.) +
                                                                                        std::pow(geom.get_constant_delta(1), 2.) +
                                                                                        std::pow(geom.get_constant_delta(2), 2.),
                                                                                        0.5);
  const DoubleTab& normale_facettes = next_time ? maillage.get_update_normale_facettes() : maillage.get_normale_facettes_old();
  const ArrOfDouble& surface_facettes = next_time ? maillage.get_update_surface_facettes() : maillage.get_surface_facettes_old();
  const int nb_facettes = next_time ? maillage.nb_facettes() : maillage.nb_facettes_old();
  const IntTab& facettes = next_time ? maillage.facettes() : maillage.facettes_old();
  const DoubleTab& sommets = next_time ? maillage.sommets() : maillage.sommets_old();

  coord_facettes.centre.resize(nb_facettes, 3);
  coord_facettes.centre = 0.;

  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    for (int dir = 0; dir < 3; dir++)
      for (int som = 0; som < 3; som++)
        coord_facettes.centre(fa7, dir) += sommets(facettes(fa7, som), dir);

  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    for (int dir = 0; dir < 3; dir++)
      coord_facettes.centre(fa7, dir) /= 3.;

  if (second_order_diffusion_interface_)
    {
      calcul_temperature_flux_interface_cut_cell_second_order(next_time,
                                                              cut_field_temperature,
                                                              lambda_liquid,
                                                              lambda_vapour,
                                                              dist_1,
                                                              dist_2,
                                                              coord_facettes.centre,
                                                              normale_facettes,
                                                              interfacial_temperature.centre,
                                                              flux_normal_interp,
                                                              interfacial_temperature.liqu_1,
                                                              interfacial_temperature.vap_1,
                                                              interfacial_temperature.liqu_2,
                                                              interfacial_temperature.vap_2,
                                                              coord_facettes.liqu_1,
                                                              coord_facettes.vap_1,
                                                              coord_facettes.liqu_2,
                                                              coord_facettes.vap_2);
    }
  else
    {
      calcul_temperature_flux_interface_cut_cell(next_time,
                                                 cut_field_temperature,
                                                 lambda_liquid,
                                                 lambda_vapour,
                                                 dist_1,
                                                 coord_facettes.centre,
                                                 normale_facettes,
                                                 interfacial_temperature.centre,
                                                 flux_normal_interp,
                                                 interfacial_temperature.liqu_1,
                                                 interfacial_temperature.vap_1,
                                                 coord_facettes.liqu_1,
                                                 coord_facettes.vap_1);
    }

  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    {
      flux_normal_interp(fa7) *= surface_facettes(fa7);
      interfacial_temperature.centre(fa7) *= surface_facettes(fa7);
      interfacial_temperature.liqu_1(fa7) *= surface_facettes(fa7);
      interfacial_temperature.vap_1(fa7) *= surface_facettes(fa7);
    }
}

void Cut_cell_diffusion_auxiliaire::compute_interfacial_temperature_local_normal(
  bool next_time,
  double lambda_liquid,
  double lambda_vapour,
  const Cut_field_double& cut_field_temperature,
  const Maillage_FT_IJK& maillage,
  const IJK_Splitting& s,
  const IJK_Interfaces& interfaces,
  Facettes_data& coord_facettes,
  Facettes_data& interfacial_temperature,
  DoubleTabFT& flux_normal_interp)
{
  const Intersections_Elem_Facettes& intersec = next_time ? maillage.intersections_elem_facettes() : maillage.intersections_elem_facettes_old();
  const DoubleTab& normale_facettes = next_time ? maillage.get_update_normale_facettes() : maillage.get_normale_facettes_old();
  const ArrOfDouble& surface_facettes = next_time ? maillage.get_update_surface_facettes() : maillage.get_surface_facettes_old();
  const int nb_facettes = next_time ? maillage.nb_facettes() : maillage.nb_facettes_old();
  const IntTab& facettes = next_time ? maillage.facettes() : maillage.facettes_old();
  const DoubleTab& sommets = next_time ? maillage.sommets() : maillage.sommets_old();

  coord_facettes.centre.resize(nb_facettes, 3);
  coord_facettes.centre = 0.;

  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    for (int dir = 0; dir < 3; dir++)
      for (int som = 0; som < 3; som++)
        coord_facettes.centre(fa7, dir) += sommets(facettes(fa7, som), dir);

  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    for (int dir = 0; dir < 3; dir++)
      coord_facettes.centre(fa7, dir) /= 3.;

  interfacial_temperature.centre.resize(nb_facettes);
  interfacial_temperature.centre = 0.;
  flux_normal_interp.resize(nb_facettes);
  flux_normal_interp = 0.;

  if ((lambda_liquid + lambda_vapour) == 0.)
    {
      Cerr << "Corrige flux used with no conductivity. Ti and Qi set to 0. " << finl;
      interfacial_temperature.centre = 0.;
      flux_normal_interp = 0.;
      return;
    }

  const int ni = s.get_nb_elem_local(0);
  const int nj = s.get_nb_elem_local(1);
  const int nk = s.get_nb_elem_local(2);
  const ArrOfInt& index_elem = intersec.index_elem();
  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              if (next_time ? (!interfaces.est_pure(interfaces.In_ft(i,j,k))) : (!interfaces.est_pure(interfaces.I_ft(i,j,k))))
                {
                  assert(maillage.ref_splitting().valeur() == s);
                  const int num_elem = s.convert_ijk_cell_to_packed(i,j,k);

                  double x_centre_cell = (i + s.get_offset_local(DIRECTION_I) + .5)*s.get_grid_geometry().get_constant_delta(DIRECTION_I) + s.get_grid_geometry().get_origin(DIRECTION_I);
                  double y_centre_cell = (j + s.get_offset_local(DIRECTION_J) + .5)*s.get_grid_geometry().get_constant_delta(DIRECTION_J) + s.get_grid_geometry().get_origin(DIRECTION_J);
                  double z_centre_cell = (k + s.get_offset_local(DIRECTION_K) + .5)*s.get_grid_geometry().get_constant_delta(DIRECTION_K) + s.get_grid_geometry().get_origin(DIRECTION_K);
                  Vecteur3 coord_centre_cell = {x_centre_cell, y_centre_cell, z_centre_cell};

                  int index = index_elem[num_elem];
                  // Boucle sur les facettes qui traversent cet element
                  while (index >= 0)
                    {
                      const Intersections_Elem_Facettes_Data& data = intersec.data_intersection(index);
                      const int fa7 = data.numero_facette_;

                      double Ti_intersection = 0.;
                      double flux_normal_intersection = 0.;
                      Vecteur3 coord = {coord_facettes.centre(fa7,0), coord_facettes.centre(fa7,1), coord_facettes.centre(fa7,2)};
                      Vecteur3 normale = {normale_facettes(fa7,0), normale_facettes(fa7,1), normale_facettes(fa7,2)};
                      calcul_temperature_flux_interface_local_normal(cut_field_temperature,
                                                                     lambda_liquid,
                                                                     lambda_vapour,
                                                                     coord_centre_cell,
                                                                     coord,
                                                                     normale,
                                                                     Ti_intersection,
                                                                     flux_normal_intersection);

                      interfacial_temperature.centre[fa7] += Ti_intersection*data.fraction_surface_intersection_;
                      flux_normal_interp[fa7] += flux_normal_intersection*data.fraction_surface_intersection_;

                      index = data.index_facette_suivante_;
                    };
                }
            }
        }
    }

  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    {
      flux_normal_interp(fa7) *= surface_facettes(fa7);
      interfacial_temperature.centre(fa7) *= surface_facettes(fa7);
    }
}

void Cut_cell_diffusion_auxiliaire::ajout_flux_interface_a_divergence_simple(Cut_field_double& cut_field_div_coeff_grad_T_volume)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_div_coeff_grad_T_volume.get_cut_cell_disc();
  for (int n = 0; n < cut_cell_disc.get_n_loc(); n++)
    {
      cut_field_div_coeff_grad_T_volume.diph_l_(n) -= flux_interface_efficace_(n);
      cut_field_div_coeff_grad_T_volume.diph_v_(n) += flux_interface_efficace_(n);
    }
}

void Cut_cell_diffusion_auxiliaire::ajout_flux_interface_a_divergence_etale(Cut_field_double& cut_field_div_coeff_grad_T_volume)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_div_coeff_grad_T_volume.get_cut_cell_disc();
  const IJK_Grid_Geometry& geom = cut_cell_disc.get_splitting().get_grid_geometry();

  for (int n = 0; n < cut_cell_disc.get_n_loc(); n++)
    {
      Int3 ijk = cut_cell_disc.get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      double x_interface = cut_cell_disc.get_interfaces().get_coord_deplacement_interface()(n,0);
      double y_interface = cut_cell_disc.get_interfaces().get_coord_deplacement_interface()(n,1);
      double z_interface = cut_cell_disc.get_interfaces().get_coord_deplacement_interface()(n,2);

      double i_interface = (x_interface - geom.get_origin(DIRECTION_I))/geom.get_constant_delta(DIRECTION_I);
      double j_interface = (y_interface - geom.get_origin(DIRECTION_J))/geom.get_constant_delta(DIRECTION_J);
      double k_interface = (z_interface - geom.get_origin(DIRECTION_K))/geom.get_constant_delta(DIRECTION_K);
      double xfact = i_interface - std::floor(i_interface);
      double yfact = j_interface - std::floor(j_interface);
      double zfact = k_interface - std::floor(k_interface);
      assert((int)std::floor(i_interface) == i);
      assert((int)std::floor(j_interface) == j);
      assert((int)std::floor(k_interface) == k);

      // Note: The neighbours also include the self
      const int number_of_neighbours = 27;
      Int3 neighbour_offset[number_of_neighbours] =
      {
        {-1,-1,-1}, {-1,-1, 0}, {-1,-1,+1},
        {-1, 0,-1}, {-1, 0, 0}, {-1, 0,+1},
        {-1,+1,-1}, {-1,+1, 0}, {-1,+1,+1},
        { 0,-1,-1}, { 0,-1, 0}, { 0,-1,+1},
        { 0, 0,-1}, { 0, 0, 0}, { 0, 0,+1},
        { 0,+1,-1}, { 0,+1, 0}, { 0,+1,+1},
        {+1,-1,-1}, {+1,-1, 0}, {+1,-1,+1},
        {+1, 0,-1}, {+1, 0, 0}, {+1, 0,+1},
        {+1,+1,-1}, {+1,+1, 0}, {+1,+1,+1}
      };

      double total_volume_l = 0.;
      double total_volume_v = 0.;
      double vol_l[number_of_neighbours];
      double vol_v[number_of_neighbours];
      for (int index_neighbour = 0; index_neighbour < number_of_neighbours; index_neighbour++)
        {
          int i_neighbour = i + neighbour_offset[index_neighbour][0];
          int j_neighbour = j + neighbour_offset[index_neighbour][1];
          int k_neighbour = k + neighbour_offset[index_neighbour][2];

          double indicatrice_neighbour = cut_cell_disc.get_interfaces().In(i_neighbour, j_neighbour, k_neighbour);

          double max_volume_contrib_x = (neighbour_offset[index_neighbour][0] == -1) ? (1-xfact) : ((neighbour_offset[index_neighbour][0] == 0) ? 1 : xfact);
          double max_volume_contrib_y = (neighbour_offset[index_neighbour][1] == -1) ? (1-yfact) : ((neighbour_offset[index_neighbour][1] == 0) ? 1 : yfact);
          double max_volume_contrib_z = (neighbour_offset[index_neighbour][2] == -1) ? (1-zfact) : ((neighbour_offset[index_neighbour][2] == 0) ? 1 : zfact);
          double max_volume = max_volume_contrib_x*max_volume_contrib_y*max_volume_contrib_z;

          vol_l[index_neighbour] = max_volume * indicatrice_neighbour;
          vol_v[index_neighbour] = max_volume * (1 - indicatrice_neighbour);
          total_volume_l += vol_l[index_neighbour];
          total_volume_v += vol_v[index_neighbour];
        }

      for (int index_neighbour = 0; index_neighbour < number_of_neighbours; index_neighbour++)
        {
          int i_neighbour = i + neighbour_offset[index_neighbour][0];
          int j_neighbour = j + neighbour_offset[index_neighbour][1];
          int k_neighbour = k + neighbour_offset[index_neighbour][2];

          int n_neighbour = cut_cell_disc.get_n(i_neighbour, j_neighbour, k_neighbour);

          if (n_neighbour >= 0)
            {
              cut_field_div_coeff_grad_T_volume.diph_l_(n) -= vol_l[index_neighbour]/total_volume_l*flux_interface_efficace_(n);
              cut_field_div_coeff_grad_T_volume.diph_v_(n) += vol_v[index_neighbour]/total_volume_v*flux_interface_efficace_(n);
            }
          else
            {
              double indicatrice_neighbour = cut_cell_disc.get_interfaces().In(i_neighbour, j_neighbour, k_neighbour);
              assert((indicatrice_neighbour == 0.) || (indicatrice_neighbour == 1.));

              int phase = (int)indicatrice_neighbour;
              if (phase == 0)
                {
                  cut_field_div_coeff_grad_T_volume.pure_(i_neighbour,j_neighbour,k_neighbour) += vol_v[index_neighbour]/total_volume_v*flux_interface_efficace_(n);
                }
              else
                {
                  cut_field_div_coeff_grad_T_volume.pure_(i_neighbour,j_neighbour,k_neighbour) -= vol_l[index_neighbour]/total_volume_l*flux_interface_efficace_(n);
                }
            }
        }
    }
}

void Cut_cell_diffusion_auxiliaire::etalement_divergence_flux_diffusifs(Cut_field_double& cut_field_div_coeff_grad_T_volume, Cut_field_double& cut_field_div_coeff_grad_T_volume_temp)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_div_coeff_grad_T_volume.get_cut_cell_disc();
  const IJK_Grid_Geometry& geom = cut_cell_disc.get_splitting().get_grid_geometry();

  cut_field_div_coeff_grad_T_volume_temp.set_to_uniform_value(0);

  for (int n = 0; n < cut_cell_disc.get_n_tot(); n++)
    {
      Int3 ijk = cut_cell_disc.get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      if (!cut_cell_disc.within_ghost(i, j, k, 1, 1))
        continue;

      // Valeur actuelle de la divergence, a etaler sur les cellules voisines
      double valeur_a_etale_l = cut_field_div_coeff_grad_T_volume.diph_l_(n);
      double valeur_a_etale_v = cut_field_div_coeff_grad_T_volume.diph_v_(n);

      // Reinitialisation de la divergence, pour remplir a nouveau plus tard
      cut_field_div_coeff_grad_T_volume.diph_l_(n) = 0.;
      cut_field_div_coeff_grad_T_volume.diph_v_(n) = 0.;

      double x_interface = cut_cell_disc.get_interfaces().get_coord_deplacement_interface()(n,0);
      double y_interface = cut_cell_disc.get_interfaces().get_coord_deplacement_interface()(n,1);
      double z_interface = cut_cell_disc.get_interfaces().get_coord_deplacement_interface()(n,2);

      double i_interface = (x_interface - geom.get_origin(DIRECTION_I))/geom.get_constant_delta(DIRECTION_I);
      double j_interface = (y_interface - geom.get_origin(DIRECTION_J))/geom.get_constant_delta(DIRECTION_J);
      double k_interface = (z_interface - geom.get_origin(DIRECTION_K))/geom.get_constant_delta(DIRECTION_K);
      double xfact = i_interface - std::floor(i_interface);
      double yfact = j_interface - std::floor(j_interface);
      double zfact = k_interface - std::floor(k_interface);
      assert((int)std::floor(i_interface) == i);
      assert((int)std::floor(j_interface) == j);
      assert((int)std::floor(k_interface) == k);

      // Note: The neighbours also include the self
      const int number_of_neighbours = 27;
      Int3 neighbour_offset[number_of_neighbours] =
      {
        {-1,-1,-1}, {-1,-1, 0}, {-1,-1,+1},
        {-1, 0,-1}, {-1, 0, 0}, {-1, 0,+1},
        {-1,+1,-1}, {-1,+1, 0}, {-1,+1,+1},
        { 0,-1,-1}, { 0,-1, 0}, { 0,-1,+1},
        { 0, 0,-1}, { 0, 0, 0}, { 0, 0,+1},
        { 0,+1,-1}, { 0,+1, 0}, { 0,+1,+1},
        {+1,-1,-1}, {+1,-1, 0}, {+1,-1,+1},
        {+1, 0,-1}, {+1, 0, 0}, {+1, 0,+1},
        {+1,+1,-1}, {+1,+1, 0}, {+1,+1,+1}
      };

      double total_volume_l = 0.;
      double total_volume_v = 0.;
      double vol_l[number_of_neighbours];
      double vol_v[number_of_neighbours];
      for (int index_neighbour = 0; index_neighbour < number_of_neighbours; index_neighbour++)
        {
          int i_neighbour = i + neighbour_offset[index_neighbour][0];
          int j_neighbour = j + neighbour_offset[index_neighbour][1];
          int k_neighbour = k + neighbour_offset[index_neighbour][2];

          double indicatrice_neighbour = cut_cell_disc.get_interfaces().In(i_neighbour, j_neighbour, k_neighbour);

          double max_volume_contrib_x = (neighbour_offset[index_neighbour][0] == -1) ? (1-xfact) : ((neighbour_offset[index_neighbour][0] == 0) ? 1 : xfact);
          double max_volume_contrib_y = (neighbour_offset[index_neighbour][1] == -1) ? (1-yfact) : ((neighbour_offset[index_neighbour][1] == 0) ? 1 : yfact);
          double max_volume_contrib_z = (neighbour_offset[index_neighbour][2] == -1) ? (1-zfact) : ((neighbour_offset[index_neighbour][2] == 0) ? 1 : zfact);
          double max_volume = max_volume_contrib_x*max_volume_contrib_y*max_volume_contrib_z;

          vol_l[index_neighbour] = max_volume * indicatrice_neighbour;
          vol_v[index_neighbour] = max_volume * (1 - indicatrice_neighbour);
          total_volume_l += vol_l[index_neighbour];
          total_volume_v += vol_v[index_neighbour];
        }

      for (int index_neighbour = 0; index_neighbour < number_of_neighbours; index_neighbour++)
        {
          int i_neighbour = i + neighbour_offset[index_neighbour][0];
          int j_neighbour = j + neighbour_offset[index_neighbour][1];
          int k_neighbour = k + neighbour_offset[index_neighbour][2];

          int n_neighbour = cut_cell_disc.get_n(i_neighbour, j_neighbour, k_neighbour);

          if (n_neighbour >= 0)
            {
              cut_field_div_coeff_grad_T_volume_temp.diph_l_(n) += vol_l[index_neighbour]/total_volume_l*valeur_a_etale_l;
              cut_field_div_coeff_grad_T_volume_temp.diph_v_(n) += vol_v[index_neighbour]/total_volume_v*valeur_a_etale_v;
            }
          else
            {
              double indicatrice_neighbour = cut_cell_disc.get_interfaces().In(i_neighbour, j_neighbour, k_neighbour);
              assert((indicatrice_neighbour == 0.) || (indicatrice_neighbour == 1.));

              int phase = (int)indicatrice_neighbour;
              if (phase == 0)
                {
                  cut_field_div_coeff_grad_T_volume_temp.pure_(i_neighbour,j_neighbour,k_neighbour) += vol_v[index_neighbour]/total_volume_v*valeur_a_etale_v;
                }
              else
                {
                  cut_field_div_coeff_grad_T_volume_temp.pure_(i_neighbour,j_neighbour,k_neighbour) += vol_l[index_neighbour]/total_volume_l*valeur_a_etale_l;
                }
            }
        }
    }

  cut_field_div_coeff_grad_T_volume.add_from(cut_field_div_coeff_grad_T_volume_temp);
}

double Cut_cell_diffusion_auxiliaire::dying_cells_flux(int num_face, int phase, int n, const FixedVector<Cut_field_double, 3>& cut_field_total_velocity, const Cut_field_double& cut_field_temperature)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();

  Int3 ijk = cut_cell_disc.get_ijk(n);
  int i = ijk[0];
  int j = ijk[1];
  int k = ijk[2];

  int dir = num_face%3;
  int decalage = num_face/3;
  int sign = decalage*2 -1;

  int di = decalage*(dir == 0);
  int dj = decalage*(dir == 1);
  int dk = decalage*(dir == 2);


  double normal_x = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,0);
  double normal_y = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,1);
  double normal_z = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,2);
  int sign_flux_interf = (flux_interface_efficace_(n) == 0.) ? 0 : (1 - 2*phase)*(2*(flux_interface_efficace_(n) > 0) - 1);

  double normal_to_face = sign*select(dir, normal_x, normal_y, normal_z);

  int n_face = cut_cell_disc.get_n_face(num_face, n, i, j, k);
  if (n_face >= 0)
    {
      double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir) : cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir);
      if (surface_efficace > 0)
        {
          return -sign*surface_efficace*normal_to_face*sign_flux_interf;
        }
      else
        {
          return 0.;
        }
    }
  else
    {
      double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().I(i+di,j+dj,k+dk) : cut_cell_disc.get_interfaces().I(i+di,j+dj,k+dk);
      assert((surface_efficace == 0) || (surface_efficace == 1));
      if (surface_efficace > 0)
        {
          return -sign*surface_efficace*normal_to_face*sign_flux_interf;
        }
      else
        {
          return 0.;
        }
    }
}

double Cut_cell_diffusion_auxiliaire::small_nascent_cells_flux(int num_face, int phase, int n, const FixedVector<Cut_field_double, 3>& cut_field_total_velocity, const Cut_field_double& cut_field_temperature)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();

  Int3 ijk = cut_cell_disc.get_ijk(n);
  int i = ijk[0];
  int j = ijk[1];
  int k = ijk[2];

  int dir = num_face%3;
  int decalage = num_face/3;
  int sign = decalage*2 -1;

  int di = decalage*(dir == 0);
  int dj = decalage*(dir == 1);
  int dk = decalage*(dir == 2);

  double normal_x = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,0);
  double normal_y = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,1);
  double normal_z = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,2);
  int sign_flux_interf = (flux_interface_efficace_(n) == 0.) ? 0 : (1 - 2*phase)*(2*(flux_interface_efficace_(n) > 0) - 1);

  double normal_to_face = sign*select(dir, normal_x, normal_y, normal_z);

  int n_face = cut_cell_disc.get_n_face(num_face, n, i, j, k);
  if (n_face >= 0)
    {
      double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir) : cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir);
      if (surface_efficace > 0)
        {
          return -sign*surface_efficace*normal_to_face*sign_flux_interf;
        }
      else
        {
          return 0.;
        }
    }
  else
    {
      double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().In(i+di,j+dj,k+dk) : cut_cell_disc.get_interfaces().In(i+di,j+dj,k+dk);
      assert((surface_efficace == 0) || (surface_efficace == 1));
      if (surface_efficace > 0)
        {
          return -sign*surface_efficace*normal_to_face*sign_flux_interf;
        }
      else
        {
          return 0.;
        }
    }
}

void Cut_cell_diffusion_auxiliaire::calcul_temperature_flux_interface(
  const IJK_Field_double& temperature, const double ldal, const double ldav,
  const double dist, const DoubleTab& positions, const DoubleTab& normal_on_interf,
  DoubleTabFT& temperature_interp,
  DoubleTabFT& flux_normal_interp,
  DoubleTabFT& temp_liqu,
  DoubleTabFT& temp_vap,
  DoubleTab& coo_liqu,
  DoubleTab& coo_vap)
{
  Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(
    positions, normal_on_interf, dist, coo_liqu);
  Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(
    positions, normal_on_interf, -dist, coo_vap);

  temp_liqu.resize(coo_liqu.dimension(0));
  ijk_interpolate_skip_unknown_points(temperature, coo_liqu, temp_liqu, 1.e31);
  temp_vap.resize(coo_vap.dimension(0));
  ijk_interpolate_skip_unknown_points(temperature, coo_vap, temp_vap, 1.e31);
  const int n_point_interp = positions.dimension(0);
  temperature_interp.resize(n_point_interp);
  flux_normal_interp.resize(n_point_interp);
  if ((ldal + ldav) == 0.)
    {
      Cerr << "Corrige flux used with no conductivity. Ti and Qi set to 0. " << finl;
      temperature_interp = 0.;
      flux_normal_interp = 0.;
      return;
    }
  for (int i_point_interp = 0; i_point_interp < n_point_interp; i_point_interp++)
    {
      const double Ti =
        (temp_liqu(i_point_interp) * ldal + temp_vap(i_point_interp) * ldav) / (ldal + ldav);

      temperature_interp(i_point_interp) = Ti;
      flux_normal_interp(i_point_interp) = ldav * (Ti - temp_vap(i_point_interp)) / dist;
    }
}

void Cut_cell_diffusion_auxiliaire::calcul_temperature_flux_interface_cut_cell(
  bool next_time, const Cut_field_double& cut_field_temperature, const double ldal, const double ldav,
  const double dist, const DoubleTab& positions, const DoubleTab& normal_on_interf,
  DoubleTabFT& temperature_interp,
  DoubleTabFT& flux_normal_interp,
  DoubleTabFT& temp_liqu,
  DoubleTabFT& temp_vap,
  DoubleTab& coo_liqu,
  DoubleTab& coo_vap)
{
  Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(
    positions, normal_on_interf, dist, coo_liqu);
  Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(
    positions, normal_on_interf, -dist, coo_vap);

  temp_liqu.resize(coo_liqu.dimension(0));
  ijk_interpolate_cut_cell_skip_unknown_points(next_time, 1, cut_field_temperature, coo_liqu, temp_liqu, 1.e31);
  temp_vap.resize(coo_vap.dimension(0));
  ijk_interpolate_cut_cell_skip_unknown_points(next_time, 0, cut_field_temperature, coo_vap, temp_vap, 1.e31);
  const int n_point_interp = positions.dimension(0);
  temperature_interp.resize(n_point_interp);
  flux_normal_interp.resize(n_point_interp);
  if ((ldal + ldav) == 0.)
    {
      Cerr << "Corrige flux used with no conductivity. Ti and Qi set to 0. " << finl;
      temperature_interp = 0.;
      flux_normal_interp = 0.;
      return;
    }
  for (int i_point_interp = 0; i_point_interp < n_point_interp; i_point_interp++)
    {
      const double Ti =
        (temp_liqu(i_point_interp) * ldal + temp_vap(i_point_interp) * ldav) / (ldal + ldav);

      temperature_interp(i_point_interp) = Ti;
      flux_normal_interp(i_point_interp) = ldav * (Ti - temp_vap(i_point_interp)) / dist;
    }
}

void Cut_cell_diffusion_auxiliaire::calcul_temperature_flux_interface_second_order(
  const IJK_Field_double& temperature, const double ldal, const double ldav,
  const double dist_1, const double dist_2, const DoubleTab& positions, const DoubleTab& normal_on_interf,
  DoubleTabFT& temperature_interp,
  DoubleTabFT& flux_normal_interp,
  DoubleTabFT& temp_liqu_1,
  DoubleTabFT& temp_vap_1,
  DoubleTabFT& temp_liqu_2,
  DoubleTabFT& temp_vap_2,
  DoubleTab& coo_liqu_1,
  DoubleTab& coo_vap_1,
  DoubleTab& coo_liqu_2,
  DoubleTab& coo_vap_2)
{
  Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(
    positions, normal_on_interf, dist_1, coo_liqu_1);
  Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(
    positions, normal_on_interf, -dist_1, coo_vap_1);
  Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(
    positions, normal_on_interf, dist_2, coo_liqu_2);
  Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(
    positions, normal_on_interf, -dist_2, coo_vap_2);

  temp_liqu_1.resize(coo_liqu_1.dimension(0));
  ijk_interpolate_skip_unknown_points(temperature, coo_liqu_1, temp_liqu_1, 1.e31);
  temp_vap_1.resize(coo_vap_1.dimension(0));
  ijk_interpolate_skip_unknown_points(temperature, coo_vap_1, temp_vap_1, 1.e31);
  temp_liqu_2.resize(coo_liqu_2.dimension(0));
  ijk_interpolate_skip_unknown_points(temperature, coo_liqu_2, temp_liqu_2, 1.e31);
  temp_vap_2.resize(coo_vap_2.dimension(0));
  ijk_interpolate_skip_unknown_points(temperature, coo_vap_2, temp_vap_2, 1.e31);
  const int n_point_interp = positions.dimension(0);
  temperature_interp.resize(n_point_interp);
  flux_normal_interp.resize(n_point_interp);
  if ((ldal + ldav) == 0.)
    {
      Cerr << "Corrige flux used with no conductivity. Ti and Qi set to 0. " << finl;
      temperature_interp = 0.;
      flux_normal_interp = 0.;
      return;
    }
  for (int i_point_interp = 0; i_point_interp < n_point_interp; i_point_interp++)
    {
      const double Ti = (ldal*temp_liqu_1(i_point_interp)/dist_1 + ldav*temp_vap_1(i_point_interp)/dist_1 + ldal*temp_liqu_2(i_point_interp)/dist_2 + ldav*temp_vap_2(i_point_interp)/dist_2 + ldal*(temp_liqu_1(i_point_interp) - temp_liqu_2(i_point_interp))/(dist_2 - dist_1) - ldav*(temp_vap_2(i_point_interp) - temp_vap_1(i_point_interp))/(dist_2 - dist_1))/(ldal/dist_1 + ldav/dist_1 + ldal/dist_2 + ldav/dist_2);

      temperature_interp(i_point_interp) = Ti;
      flux_normal_interp(i_point_interp) = - ldav * ((temp_vap_1(i_point_interp) - Ti)/dist_1 + (temp_vap_2(i_point_interp) - Ti)/dist_2 - (temp_vap_2(i_point_interp) - temp_vap_1(i_point_interp))/(dist_2 - dist_1));
    }
}

void Cut_cell_diffusion_auxiliaire::calcul_temperature_flux_interface_cut_cell_second_order(
  bool next_time, const Cut_field_double& cut_field_temperature, const double ldal, const double ldav,
  const double dist_1, const double dist_2, const DoubleTab& positions, const DoubleTab& normal_on_interf,
  DoubleTabFT& temperature_interp,
  DoubleTabFT& flux_normal_interp,
  DoubleTabFT& temp_liqu_1,
  DoubleTabFT& temp_vap_1,
  DoubleTabFT& temp_liqu_2,
  DoubleTabFT& temp_vap_2,
  DoubleTab& coo_liqu_1,
  DoubleTab& coo_vap_1,
  DoubleTab& coo_liqu_2,
  DoubleTab& coo_vap_2)
{
  Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(
    positions, normal_on_interf, dist_1, coo_liqu_1);
  Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(
    positions, normal_on_interf, -dist_1, coo_vap_1);
  Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(
    positions, normal_on_interf, dist_2, coo_liqu_2);
  Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(
    positions, normal_on_interf, -dist_2, coo_vap_2);

  temp_liqu_1.resize(coo_liqu_1.dimension(0));
  ijk_interpolate_cut_cell_skip_unknown_points(next_time, 1, cut_field_temperature, coo_liqu_1, temp_liqu_1, 1.e31);
  temp_vap_1.resize(coo_vap_1.dimension(0));
  ijk_interpolate_cut_cell_skip_unknown_points(next_time, 0, cut_field_temperature, coo_vap_1, temp_vap_1, 1.e31);
  temp_liqu_2.resize(coo_liqu_2.dimension(0));
  ijk_interpolate_cut_cell_skip_unknown_points(next_time, 1, cut_field_temperature, coo_liqu_2, temp_liqu_2, 1.e31);
  temp_vap_2.resize(coo_vap_2.dimension(0));
  ijk_interpolate_cut_cell_skip_unknown_points(next_time, 0, cut_field_temperature, coo_vap_2, temp_vap_2, 1.e31);
  const int n_point_interp = positions.dimension(0);
  temperature_interp.resize(n_point_interp);
  flux_normal_interp.resize(n_point_interp);
  if ((ldal + ldav) == 0.)
    {
      Cerr << "Corrige flux used with no conductivity. Ti and Qi set to 0. " << finl;
      temperature_interp = 0.;
      flux_normal_interp = 0.;
      return;
    }
  for (int i_point_interp = 0; i_point_interp < n_point_interp; i_point_interp++)
    {
      const double Ti = (ldal*temp_liqu_1(i_point_interp)/dist_1 + ldav*temp_vap_1(i_point_interp)/dist_1 + ldal*temp_liqu_2(i_point_interp)/dist_2 + ldav*temp_vap_2(i_point_interp)/dist_2 + ldal*(temp_liqu_1(i_point_interp) - temp_liqu_2(i_point_interp))/(dist_2 - dist_1) - ldav*(temp_vap_2(i_point_interp) - temp_vap_1(i_point_interp))/(dist_2 - dist_1))/(ldal/dist_1 + ldav/dist_1 + ldal/dist_2 + ldav/dist_2);

      temperature_interp(i_point_interp) = Ti;
      flux_normal_interp(i_point_interp) = - ldav * ((temp_vap_1(i_point_interp) - Ti)/dist_1 + (temp_vap_2(i_point_interp) - Ti)/dist_2 - (temp_vap_2(i_point_interp) - temp_vap_1(i_point_interp))/(dist_2 - dist_1));
    }
}

void Cut_cell_diffusion_auxiliaire::calcul_temperature_flux_interface_local_normal(
  const Cut_field_double& cut_field_temperature, const double ldal, const double ldav,
  const Vecteur3& position_centre_cell, const Vecteur3& positions, const Vecteur3& normal_on_interf,
  double& temperature_interp, double& flux_normal_interp)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();

  const IJK_Grid_Geometry& geom = cut_cell_disc.get_splitting().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);

  double normal_x = normal_on_interf[0];
  double normal_y = normal_on_interf[1];
  double normal_z = normal_on_interf[2];
  double norm_normal = sqrt(normal_x*normal_x + normal_y*normal_y + normal_z*normal_z);
  normal_x /= norm_normal;
  normal_y /= norm_normal;
  normal_z /= norm_normal;

  double x = positions[0];
  double y = positions[1];
  double z = positions[2];

  double x_centre_cell = position_centre_cell[0];
  double y_centre_cell = position_centre_cell[1];
  double z_centre_cell = position_centre_cell[2];

  int i = cut_cell_disc.get_i_selon_dir(0, x_centre_cell);
  int j = cut_cell_disc.get_i_selon_dir(1, y_centre_cell);
  int k = cut_cell_disc.get_i_selon_dir(2, z_centre_cell);

  int i_aperio = (int)std::floor((x_centre_cell - geom.get_origin(DIRECTION_I))/dx) - cut_cell_disc.get_splitting().get_offset_local(DIRECTION_I);
  int j_aperio = (int)std::floor((y_centre_cell - geom.get_origin(DIRECTION_J))/dy) - cut_cell_disc.get_splitting().get_offset_local(DIRECTION_J);
  int k_aperio = (int)std::floor((z_centre_cell - geom.get_origin(DIRECTION_K))/dz) - cut_cell_disc.get_splitting().get_offset_local(DIRECTION_K);
  if (cut_cell_disc.within_ghost(i_aperio, j_aperio, k_aperio, 0, 0))
    {
      assert(i == i_aperio);
      assert(j == j_aperio);
      assert(k == k_aperio);
    }

  int n = cut_cell_disc.get_n(i,j,k);
  assert(n >= 0);

  double temp_liqu = cut_field_temperature.diph_l_(n);
  double temp_vap = cut_field_temperature.diph_v_(n);

  double old_indicatrice = cut_cell_disc.get_interfaces().I(i,j,k);
  double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);

  double x_bord_cellule = geom.get_origin(DIRECTION_I) + dx*i_aperio;
  double y_bord_cellule = geom.get_origin(DIRECTION_J) + dx*j_aperio;
  double z_bord_cellule = geom.get_origin(DIRECTION_K) + dx*k_aperio;

  double liqu_x = x_bord_cellule + dx*cut_cell_disc.get_interfaces().get_barycentre(true, 0, 1, i, j, k, old_indicatrice, next_indicatrice);
  double liqu_y = y_bord_cellule + dy*cut_cell_disc.get_interfaces().get_barycentre(true, 1, 1, i, j, k, old_indicatrice, next_indicatrice);
  double liqu_z = z_bord_cellule + dz*cut_cell_disc.get_interfaces().get_barycentre(true, 2, 1, i, j, k, old_indicatrice, next_indicatrice);

  double vap_x = x_bord_cellule + dx*cut_cell_disc.get_interfaces().get_barycentre(true, 0, 0, i, j, k, old_indicatrice, next_indicatrice);
  double vap_y = y_bord_cellule + dy*cut_cell_disc.get_interfaces().get_barycentre(true, 1, 0, i, j, k, old_indicatrice, next_indicatrice);
  double vap_z = z_bord_cellule + dz*cut_cell_disc.get_interfaces().get_barycentre(true, 2, 0, i, j, k, old_indicatrice, next_indicatrice);

  double dist_liqu = std::abs(normal_x*(liqu_x - x) + normal_y*(liqu_y - y) + normal_z*(liqu_z - z));
  double dist_vap = std::abs(normal_x*(vap_x - x) + normal_y*(vap_y - y) + normal_z*(vap_z - z));

  const double Ti = (temp_liqu*ldal/dist_liqu + temp_vap*ldav/dist_vap) / (ldal/dist_liqu + ldav/dist_vap);


  temperature_interp = Ti;
  flux_normal_interp = ldav * (Ti - temp_vap) / dist_vap;
}

