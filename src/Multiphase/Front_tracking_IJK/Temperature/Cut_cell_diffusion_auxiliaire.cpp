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

void Cut_cell_diffusion_auxiliaire::ajout_flux_interface_a_divergence_simple(Cut_field_double& cut_field_div_coeff_grad_T_volume)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_div_coeff_grad_T_volume.get_cut_cell_disc();
  for (int n = 0; n < cut_cell_disc.get_n_loc(); n++)
    {
      cut_field_div_coeff_grad_T_volume.diph_l_(n) -= flux_interface_efficace_(n);
      cut_field_div_coeff_grad_T_volume.diph_v_(n) += flux_interface_efficace_(n);
    }
}

double Cut_cell_diffusion_auxiliaire::dying_cells_flux(int num_face, int phase, int n, const Cut_field_vector3_double& cut_field_total_velocity, const Cut_field_double& cut_field_temperature)
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

double Cut_cell_diffusion_auxiliaire::small_nascent_cells_flux(int num_face, int phase, int n, const Cut_field_vector3_double& cut_field_total_velocity, const Cut_field_double& cut_field_temperature)
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

