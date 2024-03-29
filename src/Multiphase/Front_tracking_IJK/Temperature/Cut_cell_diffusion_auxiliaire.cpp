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
#include <Process.h>
#include <stat_counters.h>
#include <IJK_FT_cut_cell.h>
#include <IJK_Navier_Stokes_tools.h>
#include <IJK_Navier_Stokes_tools_cut_cell.h>

void Cut_cell_diffusion_auxiliaire::calculer_flux_interface(METHODE_FLUX_INTERFACE methode_flux_interface, double scaled_distance, double lambda_liquid, double lambda_vapour, Cut_field_scalar& cut_field_temperature, REF(IJK_FT_cut_cell)& ref_ijk_ft, const IJK_Field_double& temperature_ns, IJK_Field_double& temperature_ft)
{
  ArrOfDouble interfacial_temperature;
  ArrOfDouble interfacial_phin_ai;

  // Transfer the field to FT splitting, needed for compute_interfacial_temperature
  ref_ijk_ft->redistrib_to_ft_elem().redistribute(temperature_ns, temperature_ft);
  temperature_ft.echange_espace_virtuel(temperature_ft.ghost());

  if (methode_flux_interface == METHODE_FLUX_INTERFACE::INTERP_PURE)
    {
      compute_interfacial_temperature2(scaled_distance, lambda_liquid, lambda_vapour, temperature_ft, ref_ijk_ft->get_geometry(), ref_ijk_ft->itfce().maillage_ft_ijk(), interfacial_temperature, interfacial_phin_ai);
    }
  else if (methode_flux_interface == METHODE_FLUX_INTERFACE::INTERP_CUT_CELL)
    {
      compute_interfacial_temperature_cut_cell(scaled_distance, lambda_liquid, lambda_vapour, cut_field_temperature, ref_ijk_ft->get_geometry(), ref_ijk_ft->itfce().maillage_ft_ijk(), interfacial_temperature, interfacial_phin_ai);
    }
  else if (methode_flux_interface == METHODE_FLUX_INTERFACE::LOCAL_CELLULE)
    {
      assert(flux_interface_ft_.get_splitting().get_nb_elem_local(0) == flux_interface_ft_.ni());
      assert(flux_interface_ft_.get_splitting().get_nb_elem_local(1) == flux_interface_ft_.nj());
      assert(flux_interface_ft_.get_splitting().get_nb_elem_local(2) == flux_interface_ft_.nk());
      compute_interfacial_temperature_local_normal(lambda_liquid, lambda_vapour, cut_field_temperature, ref_ijk_ft->itfce().maillage_ft_ijk(), flux_interface_ft_.get_splitting(), ref_ijk_ft->itfce(), interfacial_temperature, interfacial_phin_ai);
    }
  else
    {
      Cerr << "Methode non reconnue pour le calcul du flux a l'interface." << finl;
      Process::exit();
    }

  // Calcul des flux sur le maillage FT
  {
    const Cut_cell_FT_Disc& cut_cell_disc = flux_interface_.get_cut_cell_disc();
    const Maillage_FT_IJK& mesh = cut_cell_disc.get_interfaces().maillage_ft_ijk();
    const Intersections_Elem_Facettes& intersec = mesh.intersections_elem_facettes();
    const IJK_Splitting& s = flux_interface_ft_.get_splitting();

    const int ni = flux_interface_ft_.ni();
    const int nj = flux_interface_ft_.nj();
    const int nk = flux_interface_ft_.nk();

    // Initialisation
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  flux_interface_ft_(i, j, k) = 0.;
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
                  if (!cut_cell_disc.get_interfaces().est_pure(cut_cell_disc.get_interfaces().I_ft(i,j,k)))
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

                      flux_interface_ft_(i,j,k) = somme_contrib;
                    }
                }
            }
        }
    }
  }

  flux_interface_ft_.echange_espace_virtuel(flux_interface_ft_.ghost());

  // Calcul du flux interface sur le domaine NS :
  ref_ijk_ft->redistrib_from_ft_elem().redistribute(flux_interface_ft_, flux_interface_ns_);
  flux_interface_ns_.echange_espace_virtuel(flux_interface_ns_.ghost());

  // Correction par la surface efficace
  const Cut_cell_FT_Disc& cut_cell_disc = flux_interface_.get_cut_cell_disc();
  const DoubleTabFT_cut_cell_scalar& surface_efficace_interface = cut_cell_disc.get_interfaces().get_surface_efficace_interface();
  const IJK_Field_double& surface_interface_next = cut_cell_disc.get_interfaces().get_surface_interface_next();
  for (int n = 0; n < cut_cell_disc.get_n_tot(); n++)
    {
      Int3 ijk = cut_cell_disc.get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      flux_interface_(n) = (surface_interface_next(i,j,k) == 0.) ? 0. : flux_interface_ns_(i,j,k) * surface_efficace_interface(n)/surface_interface_next(i,j,k);
    }
}

// Copie de IJK_Thermal_base::compute_interfacial_temperature2
void Cut_cell_diffusion_auxiliaire::compute_interfacial_temperature2(double scaled_distance,
                                                                     double lambda_liquid,
                                                                     double lambda_vapour,
                                                                     const IJK_Field_double& temperature_ft,
                                                                     const IJK_Grid_Geometry& geom,
                                                                     const Maillage_FT_IJK& maillage,
                                                                     ArrOfDouble& interfacial_temperature,
                                                                     ArrOfDouble& flux_normal_interp)
{
  const double dist = scaled_distance/sqrt(3) * std::pow(std::pow(geom.get_constant_delta(0), 2.) +
                                                         std::pow(geom.get_constant_delta(1), 2.) +
                                                         std::pow(geom.get_constant_delta(2), 2.),
                                                         0.5);
  const DoubleTab& normale_facettes = maillage.get_update_normale_facettes();
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  const int nb_facettes = maillage.nb_facettes();
  const IntTab& facettes = maillage.facettes();
  const DoubleTab& sommets = maillage.sommets();
  DoubleTab coord_facettes;
  coord_facettes.resize(nb_facettes, 3);
  coord_facettes = 0.;
  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    for (int dir = 0; dir < 3; dir++)
      for (int som = 0; som < 3; som++)
        coord_facettes(fa7, dir) += sommets(facettes(fa7, som), dir);

  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    for (int dir = 0; dir < 3; dir++)
      coord_facettes(fa7, dir) /= 3.;

  DoubleTab coo_liqu, coo_vap;
  ArrOfDouble temp_liqu, temp_vap;
  temp_liqu.set_smart_resize(1);
  temp_vap.set_smart_resize(1);
  coo_liqu.set_smart_resize(1);
  coo_vap.set_smart_resize(1);
  calcul_temperature_flux_interface(temperature_ft,
                                    lambda_liquid,
                                    lambda_vapour,
                                    dist,
                                    coord_facettes,
                                    normale_facettes,
                                    interfacial_temperature,
                                    flux_normal_interp,
                                    temp_liqu,
                                    temp_vap,
                                    coo_liqu,
                                    coo_vap);
  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    {
      flux_normal_interp(fa7) *= surface_facettes(fa7);
      interfacial_temperature(fa7) *= surface_facettes(fa7);
    }
}

void Cut_cell_diffusion_auxiliaire::compute_interfacial_temperature_cut_cell(double scaled_distance,
                                                                             double lambda_liquid,
                                                                             double lambda_vapour,
                                                                             Cut_field_scalar& cut_field_temperature,
                                                                             const IJK_Grid_Geometry& geom,
                                                                             const Maillage_FT_IJK& maillage,
                                                                             ArrOfDouble& interfacial_temperature,
                                                                             ArrOfDouble& flux_normal_interp)
{
  const double dist = scaled_distance/sqrt(3) * std::pow(std::pow(geom.get_constant_delta(0), 2.) +
                                                         std::pow(geom.get_constant_delta(1), 2.) +
                                                         std::pow(geom.get_constant_delta(2), 2.),
                                                         0.5);
  const DoubleTab& normale_facettes = maillage.get_update_normale_facettes();
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  const int nb_facettes = maillage.nb_facettes();
  const IntTab& facettes = maillage.facettes();
  const DoubleTab& sommets = maillage.sommets();
  DoubleTab coord_facettes;
  coord_facettes.resize(nb_facettes, 3);
  coord_facettes = 0.;
  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    for (int dir = 0; dir < 3; dir++)
      for (int som = 0; som < 3; som++)
        coord_facettes(fa7, dir) += sommets(facettes(fa7, som), dir);

  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    for (int dir = 0; dir < 3; dir++)
      coord_facettes(fa7, dir) /= 3.;

  DoubleTab coo_liqu, coo_vap;
  ArrOfDouble temp_liqu, temp_vap;
  temp_liqu.set_smart_resize(1);
  temp_vap.set_smart_resize(1);
  coo_liqu.set_smart_resize(1);
  coo_vap.set_smart_resize(1);
  calcul_temperature_flux_interface_cut_cell(cut_field_temperature,
                                             lambda_liquid,
                                             lambda_vapour,
                                             dist,
                                             coord_facettes,
                                             normale_facettes,
                                             interfacial_temperature,
                                             flux_normal_interp,
                                             temp_liqu,
                                             temp_vap,
                                             coo_liqu,
                                             coo_vap);
  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    {
      flux_normal_interp(fa7) *= surface_facettes(fa7);
      interfacial_temperature(fa7) *= surface_facettes(fa7);
    }
}

void Cut_cell_diffusion_auxiliaire::compute_interfacial_temperature_local_normal(
  double lambda_liquid,
  double lambda_vapour,
  Cut_field_scalar& cut_field_temperature,
  const Maillage_FT_IJK& maillage,
  const IJK_Splitting& s,
  const IJK_Interfaces& interfaces,
  ArrOfDouble& interfacial_temperature,
  ArrOfDouble& flux_normal_interp)
{
  const Intersections_Elem_Facettes& intersec = maillage.intersections_elem_facettes();
  const DoubleTab& normale_facettes = maillage.get_update_normale_facettes();
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  const int nb_facettes = maillage.nb_facettes();
  const IntTab& facettes = maillage.facettes();
  const DoubleTab& sommets = maillage.sommets();
  DoubleTab coord_facettes;
  coord_facettes.resize(nb_facettes, 3);

  coord_facettes = 0.;
  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    for (int dir = 0; dir < 3; dir++)
      for (int som = 0; som < 3; som++)
        coord_facettes(fa7, dir) += sommets(facettes(fa7, som), dir);

  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    for (int dir = 0; dir < 3; dir++)
      coord_facettes(fa7, dir) /= 3.;

  interfacial_temperature.resize_array(nb_facettes);
  interfacial_temperature = 0.;
  flux_normal_interp.resize_array(nb_facettes);
  flux_normal_interp = 0.;

  if ((lambda_liquid + lambda_vapour) == 0.)
    {
      Cerr << "Corrige flux used with no conductivity. Ti and Qi set to 0. " << finl;
      interfacial_temperature = 0.;
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
              if (!interfaces.est_pure(interfaces.I_ft(i,j,k)))
                {
                  assert(maillage.ref_splitting().valeur() == s);
                  const int num_elem = s.convert_ijk_cell_to_packed(i,j,k);

                  double x_centre_cell = (i + .5)*s.get_grid_geometry().get_constant_delta(DIRECTION_I) + s.get_grid_geometry().get_origin(DIRECTION_I);
                  double y_centre_cell = (j + .5)*s.get_grid_geometry().get_constant_delta(DIRECTION_J) + s.get_grid_geometry().get_origin(DIRECTION_J);
                  double z_centre_cell = (k + .5)*s.get_grid_geometry().get_constant_delta(DIRECTION_K) + s.get_grid_geometry().get_origin(DIRECTION_K);
                  Vecteur3 coord_centre_cell = {x_centre_cell, y_centre_cell, z_centre_cell};

                  int index = index_elem[num_elem];
                  // Boucle sur les facettes qui traversent cet element
                  while (index >= 0)
                    {
                      const Intersections_Elem_Facettes_Data& data = intersec.data_intersection(index);
                      const int fa7 = data.numero_facette_;

                      double Ti_intersection = 0.;
                      double flux_normal_intersection = 0.;
                      Vecteur3 coord = {coord_facettes(fa7,0), coord_facettes(fa7,1), coord_facettes(fa7,2)};
                      Vecteur3 normale = {normale_facettes(fa7,0), normale_facettes(fa7,1), normale_facettes(fa7,2)};
                      calcul_temperature_flux_interface_local_normal(cut_field_temperature,
                                                                     lambda_liquid,
                                                                     lambda_vapour,
                                                                     coord_centre_cell,
                                                                     coord,
                                                                     normale,
                                                                     Ti_intersection,
                                                                     flux_normal_intersection);

                      interfacial_temperature[fa7] += Ti_intersection*data.fraction_surface_intersection_;
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
      interfacial_temperature(fa7) *= surface_facettes(fa7);
    }
}

void Cut_cell_diffusion_auxiliaire::ajout_flux_interface_a_divergence_simple(Cut_field_scalar& cut_field_div_coeff_grad_T_volume)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_div_coeff_grad_T_volume.get_cut_cell_disc();
  for (int n = 0; n < cut_cell_disc.get_n_loc(); n++)
    {
      cut_field_div_coeff_grad_T_volume.diph_l_(n) -= flux_interface_(n);
      cut_field_div_coeff_grad_T_volume.diph_v_(n) += flux_interface_(n);
    }
}

void Cut_cell_diffusion_auxiliaire::ajout_flux_interface_a_divergence_etale(Cut_field_scalar& cut_field_div_coeff_grad_T_volume)
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
              cut_field_div_coeff_grad_T_volume.diph_l_(n) -= vol_l[index_neighbour]/total_volume_l*flux_interface_(n);
              cut_field_div_coeff_grad_T_volume.diph_v_(n) += vol_v[index_neighbour]/total_volume_v*flux_interface_(n);
            }
          else
            {
              double indicatrice_neighbour = cut_cell_disc.get_interfaces().In(i_neighbour, j_neighbour, k_neighbour);
              assert((indicatrice_neighbour == 0.) || (indicatrice_neighbour == 1.));

              int phase = (int)indicatrice_neighbour;
              if (phase == 0)
                {
                  cut_field_div_coeff_grad_T_volume.pure_(i_neighbour,j_neighbour,k_neighbour) += vol_v[index_neighbour]/total_volume_v*flux_interface_(n);
                }
              else
                {
                  cut_field_div_coeff_grad_T_volume.pure_(i_neighbour,j_neighbour,k_neighbour) -= vol_l[index_neighbour]/total_volume_l*flux_interface_(n);
                }
            }
        }
    }
}

void Cut_cell_diffusion_auxiliaire::etalement_divergence_flux_diffusifs(Cut_field_scalar& cut_field_div_coeff_grad_T_volume, Cut_field_scalar& cut_field_div_coeff_grad_T_volume_temp)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_div_coeff_grad_T_volume.get_cut_cell_disc();
  const IJK_Grid_Geometry& geom = cut_cell_disc.get_splitting().get_grid_geometry();

  cut_field_div_coeff_grad_T_volume_temp.pure_.data()=0;
  cut_field_div_coeff_grad_T_volume_temp.set_valeur_cellules_diphasiques(0);

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

void Cut_cell_diffusion_auxiliaire::calcul_temperature_flux_interface(
  const IJK_Field_double& temperature, const double ldal, const double ldav,
  const double dist, const DoubleTab& positions, const DoubleTab& normal_on_interf,
  ArrOfDouble& temperature_interp,
  ArrOfDouble& flux_normal_interp,
  ArrOfDouble& temp_liqu,
  ArrOfDouble& temp_vap,
  DoubleTab& coo_liqu,
  DoubleTab& coo_vap)
{
  Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(
    positions, normal_on_interf, dist, coo_liqu);
  Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(
    positions, normal_on_interf, -dist, coo_vap);

  temp_liqu.resize_array(coo_liqu.dimension(0));
  ijk_interpolate_skip_unknown_points(temperature, coo_liqu, temp_liqu, 1.e31);
  temp_vap.resize(coo_vap.dimension(0));
  ijk_interpolate_skip_unknown_points(temperature, coo_vap, temp_vap, 1.e31);
  const int n_point_interp = positions.dimension(0);
  temperature_interp.resize_array(n_point_interp);
  flux_normal_interp.resize_array(n_point_interp);
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

      if (temp_liqu(i_point_interp) > 1.e9)
        {
          Cerr << "Problem temperature interface" << finl;
          temp_liqu(i_point_interp) = 1.e9;
        }
      temperature_interp(i_point_interp) = Ti;
      flux_normal_interp(i_point_interp) = ldav * (Ti - temp_vap(i_point_interp)) / dist;
    }
}

void Cut_cell_diffusion_auxiliaire::calcul_temperature_flux_interface_cut_cell(
  Cut_field_scalar& cut_field_temperature, const double ldal, const double ldav,
  const double dist, const DoubleTab& positions, const DoubleTab& normal_on_interf,
  ArrOfDouble& temperature_interp,
  ArrOfDouble& flux_normal_interp,
  ArrOfDouble& temp_liqu,
  ArrOfDouble& temp_vap,
  DoubleTab& coo_liqu,
  DoubleTab& coo_vap)
{
  Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(
    positions, normal_on_interf, dist, coo_liqu);
  Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(
    positions, normal_on_interf, -dist, coo_vap);

  temp_liqu.resize_array(coo_liqu.dimension(0));
  ijk_interpolate_cut_cell(1, cut_field_temperature, coo_liqu, temp_liqu);
  temp_vap.resize(coo_vap.dimension(0));
  ijk_interpolate_cut_cell(0, cut_field_temperature, coo_vap, temp_vap);
  const int n_point_interp = positions.dimension(0);
  temperature_interp.resize_array(n_point_interp);
  flux_normal_interp.resize_array(n_point_interp);
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

      if (temp_liqu(i_point_interp) > 1.e9)
        {
          Cerr << "Problem cut_field_temperature interface" << finl;
          temp_liqu(i_point_interp) = 1.e9;
        }
      temperature_interp(i_point_interp) = Ti;
      flux_normal_interp(i_point_interp) = ldav * (Ti - temp_vap(i_point_interp)) / dist;
    }
}

void Cut_cell_diffusion_auxiliaire::calcul_temperature_flux_interface_local_normal(
  Cut_field_scalar& cut_field_temperature, const double ldal, const double ldav,
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

  int i = cut_cell_disc.get_i_selon_dir(0, x_centre_cell, 2, cut_cell_disc.get_splitting(), false);
  int j = cut_cell_disc.get_i_selon_dir(1, y_centre_cell, 2, cut_cell_disc.get_splitting(), false);
  int k = cut_cell_disc.get_i_selon_dir(2, z_centre_cell, 2, cut_cell_disc.get_splitting(), false);

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

  double liqu_x = x_bord_cellule + dx*cut_cell_disc.get_interfaces().get_barycentre_next(0, 1, i, j, k, old_indicatrice, next_indicatrice);
  double liqu_y = y_bord_cellule + dy*cut_cell_disc.get_interfaces().get_barycentre_next(1, 1, i, j, k, old_indicatrice, next_indicatrice);
  double liqu_z = z_bord_cellule + dz*cut_cell_disc.get_interfaces().get_barycentre_next(2, 1, i, j, k, old_indicatrice, next_indicatrice);

  double vap_x = x_bord_cellule + dx*cut_cell_disc.get_interfaces().get_barycentre_next(0, 0, i, j, k, old_indicatrice, next_indicatrice);
  double vap_y = y_bord_cellule + dy*cut_cell_disc.get_interfaces().get_barycentre_next(1, 0, i, j, k, old_indicatrice, next_indicatrice);
  double vap_z = z_bord_cellule + dz*cut_cell_disc.get_interfaces().get_barycentre_next(2, 0, i, j, k, old_indicatrice, next_indicatrice);

  double dist_liqu = std::abs(normal_x*(liqu_x - x) + normal_y*(liqu_y - y) + normal_z*(liqu_z - z));
  double dist_vap = std::abs(normal_x*(vap_x - x) + normal_y*(vap_y - y) + normal_z*(vap_z - z));

  const double Ti = (temp_liqu*ldal/dist_liqu + temp_vap*ldav/dist_vap) / (ldal/dist_liqu + ldav/dist_vap);

  if (temp_liqu > 1.e9)
    {
      Cerr << "Problem cut_field_temperature interface" << finl;
      temp_liqu = 1.e9;
    }

  temperature_interp = Ti;
  flux_normal_interp = ldav * (Ti - temp_vap) / dist_vap;
}

