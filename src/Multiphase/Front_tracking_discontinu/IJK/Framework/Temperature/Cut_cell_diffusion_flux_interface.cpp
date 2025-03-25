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

#include <Cut_cell_diffusion_flux_interface.h>
#include <Cut_cell_FT_Disc.h>
#include <IJK_Thermal_base.h>
#include <Process.h>
#include <stat_counters.h>
#include <IJK_Navier_Stokes_tools.h>
#include <IJK_Navier_Stokes_tools_cut_cell.h>

void calculer_flux_interface_sur_facettes(METHODE_FLUX_INTERFACE methode_flux_interface,
                                          bool next_time,
                                          double coeff_liquid,
                                          double coeff_vapour,
                                          DoubleTabFT& interfacial_temperature,
                                          DoubleTabFT& interfacial_phin_ai,
                                          const Cut_field_double& cut_field_temperature,
                                          const Facettes_Interp_FT& cut_cell_facettes_interpolation);
void calculer_flux_interface_sur_maillage_ft(bool next_time,
                                             const DoubleTabFT& interfacial_phin_ai,
                                             const Cut_cell_FT_Disc& cut_cell_disc,
                                             IJK_Field_double& flux_interface_ft);

void ajout_flux_interface_a_divergence_simple(int phase, const DoubleTabFT_cut_cell_scalar& flux_interface_efficace, Cut_field_double& cut_field_d_field);

void calcul_temperature_flux_interface(bool next_time,
                                       bool cut_cell,
                                       bool no_jump,
                                       const Cut_field_double& cut_field,
                                       const Facettes_Interp_FT& cut_cell_facettes_interpolation,
                                       const double ldal,
                                       const double ldav,
                                       DoubleTabFT& temperature_interp,
                                       DoubleTabFT& flux_normal_interp);
void calcul_velocity_flux_interface(bool next_time,
                                    bool cut_cell,
                                    bool no_jump,
                                    bool no_transpose,
                                    const Cut_field_vector3_double& cut_field_velocity,
                                    const Facettes_Interp_FT& cut_cell_facettes_interpolation,
                                    const double ldal,
                                    const double ldav,
                                    FixedVector<DoubleTabFT, 3>& velocity_interp,
                                    FixedVector<FixedVector<DoubleTabFT, 3>, 2>& flux_normal_interp);

void calculer_flux_interface_implementation(bool next_time,
                                            METHODE_FLUX_INTERFACE methode_flux_interface,
                                            double lambda_liquid,
                                            double lambda_vapour,
                                            DoubleTabFT& interfacial_temperature,
                                            DoubleTabFT& interfacial_phin_ai,
                                            const Cut_field_double& cut_field_temperature,
                                            const Facettes_Interp_FT& cut_cell_facettes_interpolation,
                                            IJK_Field_double& flux_interface_ft);

void ajout_flux_interface_a_divergence_implementation(int phase, const DoubleTabFT_cut_cell_scalar& flux_interface_efficace, Cut_field_double& cut_field_d_field);



void calculer_flux_interface_sur_facettes(METHODE_FLUX_INTERFACE methode_flux_interface,
                                          bool next_time,
                                          double coeff_liquid,
                                          double coeff_vapour,
                                          DoubleTabFT& interfacial_temperature,
                                          DoubleTabFT& interfacial_phin_ai,
                                          const Cut_field_double& cut_field_temperature,
                                          const Facettes_Interp_FT& cut_cell_facettes_interpolation)
{
  if (methode_flux_interface == METHODE_FLUX_INTERFACE::INTERP_PURE || methode_flux_interface == METHODE_FLUX_INTERFACE::INTERP_CUT_CELL
      || methode_flux_interface == METHODE_FLUX_INTERFACE::INTERP_PURE_NO_JUMP || methode_flux_interface == METHODE_FLUX_INTERFACE::INTERP_CUT_CELL_NO_JUMP)
    {
      bool cut_cell = (methode_flux_interface == METHODE_FLUX_INTERFACE::INTERP_CUT_CELL || methode_flux_interface == METHODE_FLUX_INTERFACE::INTERP_CUT_CELL_NO_JUMP);
      bool no_jump = (methode_flux_interface == METHODE_FLUX_INTERFACE::INTERP_PURE_NO_JUMP || methode_flux_interface == METHODE_FLUX_INTERFACE::INTERP_CUT_CELL_NO_JUMP);

      if (!cut_cell)
        {
          if (next_time)
            {
              assert(cut_field_temperature.check_agreement_diph_pure_cellules_finalement_pures());
            }
          else
            {
              assert(cut_field_temperature.check_agreement_diph_pure_cellules_initialement_pures());
            }

          assert(cut_field_temperature.check_agreement_tableau_pure_cellules_diphasiques(next_time));
        }

      calcul_temperature_flux_interface(next_time,
                                        cut_cell,
                                        no_jump,
                                        cut_field_temperature,
                                        cut_cell_facettes_interpolation,
                                        coeff_liquid,
                                        coeff_vapour,
                                        interfacial_temperature,
                                        interfacial_phin_ai);
    }
  else
    {
      Cerr << "Methode non reconnue pour le calcul du flux a l'interface." << finl;
      Process::exit();
    }

  const Maillage_FT_IJK& maillage = next_time ? cut_cell_facettes_interpolation.maillage_ft_ijk() : cut_cell_facettes_interpolation.old_maillage_ft_ijk();
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  const int nb_facettes = maillage.nb_facettes();
  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    {
      interfacial_phin_ai(fa7) *= surface_facettes(fa7);
      interfacial_temperature(fa7) *= surface_facettes(fa7);
    }
}

void calculer_flux_interface_sur_maillage_ft(bool next_time,
                                             const DoubleTabFT& interfacial_phin_ai,
                                             const Cut_cell_FT_Disc& cut_cell_disc,
                                             IJK_Field_double& flux_interface_ft)
{
  // Calcul des flux sur le maillage FT
  {
    const Maillage_FT_IJK& mesh = next_time ? cut_cell_disc.get_interfaces().maillage_ft_ijk() : cut_cell_disc.get_interfaces().old_maillage_ft_ijk();
    const Intersections_Elem_Facettes& intersec = mesh.intersections_elem_facettes();

    const int ni = flux_interface_ft.ni();
    const int nj = flux_interface_ft.nj();
    const int nk = flux_interface_ft.nk();

    // Initialisation
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  flux_interface_ft(i, j, k) = 0.;
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
                  if (next_time ? (!cut_cell_disc.get_interfaces().est_pure(cut_cell_disc.get_interfaces().In(i,j,k))) : (!cut_cell_disc.get_interfaces().est_pure(cut_cell_disc.get_interfaces().I(i,j,k))))
                    {
                      double somme_contrib = 0.;
                      const int num_elem = mesh.ref_splitting()->convert_ijk_cell_to_packed(i,j,k);

                      int index = index_elem[num_elem];
                      // Boucle sur les facettes qui traversent cet element
                      while (index >= 0)
                        {
                          const Intersections_Elem_Facettes_Data& data = intersec.data_intersection(index);
                          const int fa7 = data.numero_facette_;
                          somme_contrib += interfacial_phin_ai[fa7] * data.fraction_surface_intersection_;

                          index = data.index_facette_suivante_;
                        };

                      flux_interface_ft(i,j,k) = somme_contrib;
                    }
                }
            }
        }
    }
  }

  flux_interface_ft.echange_espace_virtuel(flux_interface_ft.ghost());
}

void calculer_flux_interface_efficace(const IJK_Field_double& flux_interface_ns_old, const IJK_Field_double& flux_interface_ns_next, DoubleTabFT_cut_cell_scalar& flux_interface_efficace)
{
  // Correction par la surface efficace
  const Cut_cell_FT_Disc& cut_cell_disc = flux_interface_efficace.get_cut_cell_disc();
  const DoubleTabFT_cut_cell_scalar& surface_efficace_interface = cut_cell_disc.get_interfaces().get_surface_efficace_interface();
  const IJK_Field_double& old_surface_interface = cut_cell_disc.get_interfaces().get_surface_interface_old();
  const IJK_Field_double& next_surface_interface = cut_cell_disc.get_interfaces().get_surface_interface_next();
  for (int n = 0; n < cut_cell_disc.get_n_tot(); n++)
    {
      Int3 ijk = cut_cell_disc.get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      if (IJK_Interfaces::est_pure(.5*(cut_cell_disc.get_interfaces().I(i, j, k) + cut_cell_disc.get_interfaces().In(i, j, k))))
        {
          // Si la cellule est purement monophasique, il n'y a pas d'intersection avec l'interface
          flux_interface_efficace(n) = 0.;
        }
      else
        {
          assert(flux_interface_ns_old.ghost() == flux_interface_ns_next.ghost());
          if (!cut_cell_disc.get_splitting().within_ghost(i, j, k, flux_interface_ns_old.ghost(), flux_interface_ns_old.ghost()))
            continue;

          double old_indicatrice = cut_cell_disc.get_interfaces().I(i,j,k);
          double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);

          if (IJK_Interfaces::devient_pure(old_indicatrice, next_indicatrice))
            {
              assert(old_surface_interface(i,j,k) != 0.);
              flux_interface_efficace(n) = (flux_interface_ns_old(i,j,k)/old_surface_interface(i,j,k)) * surface_efficace_interface(n);
            }
          else if (IJK_Interfaces::devient_diphasique(old_indicatrice, next_indicatrice))
            {
              assert(next_surface_interface(i,j,k) != 0.);
              flux_interface_efficace(n) = (flux_interface_ns_next(i,j,k)/next_surface_interface(i,j,k)) * surface_efficace_interface(n);
            }
          else
            {
              assert(old_surface_interface(i,j,k) != 0.);
              assert(next_surface_interface(i,j,k) != 0.);
              flux_interface_efficace(n) = .5 * (flux_interface_ns_old(i,j,k)/old_surface_interface(i,j,k) + flux_interface_ns_next(i,j,k)/next_surface_interface(i,j,k)) * surface_efficace_interface(n);
            }
        }
    }
}

void calculer_flux_interface_implementation(bool next_time,
                                            METHODE_FLUX_INTERFACE methode_flux_interface,
                                            double lambda_liquid,
                                            double lambda_vapour,
                                            DoubleTabFT& interfacial_temperature,
                                            DoubleTabFT& interfacial_phin_ai,
                                            const Cut_field_double& cut_field_temperature,
                                            const Facettes_Interp_FT& cut_cell_facettes_interpolation,
                                            IJK_Field_double& flux_interface_ft)
{
  calculer_flux_interface_sur_facettes(methode_flux_interface, next_time, lambda_liquid, lambda_vapour, interfacial_temperature, interfacial_phin_ai, cut_field_temperature, cut_cell_facettes_interpolation);

  calculer_flux_interface_sur_maillage_ft(next_time, interfacial_phin_ai, cut_field_temperature.get_cut_cell_disc(), flux_interface_ft);
}

void calculer_flux_interface_next(METHODE_FLUX_INTERFACE methode_flux_interface,
                                  double lambda_liquid,
                                  double lambda_vapour,
                                  DoubleTabFT& interfacial_temperature,
                                  DoubleTabFT& interfacial_phin_ai,
                                  const Cut_field_double& cut_field_temperature,
                                  const Facettes_Interp_FT& cut_cell_facettes_interpolation,
                                  IJK_Field_double& flux_interface_ft)
{
  calculer_flux_interface_implementation(true /* next time */, methode_flux_interface, lambda_liquid, lambda_vapour, interfacial_temperature, interfacial_phin_ai, cut_field_temperature, cut_cell_facettes_interpolation, flux_interface_ft);
}

void calculer_flux_interface_old(METHODE_FLUX_INTERFACE methode_flux_interface,
                                 double lambda_liquid,
                                 double lambda_vapour,
                                 DoubleTabFT& interfacial_temperature,
                                 DoubleTabFT& interfacial_phin_ai,
                                 const Cut_field_double& cut_field_temperature,
                                 const Facettes_Interp_FT& cut_cell_facettes_interpolation,
                                 IJK_Field_double& flux_interface_ft)
{
  calculer_flux_interface_implementation(false /* old time */, methode_flux_interface, lambda_liquid, lambda_vapour, interfacial_temperature, interfacial_phin_ai, cut_field_temperature, cut_cell_facettes_interpolation, flux_interface_ft);
}

void ajout_flux_interface_a_divergence(const DoubleTabFT_cut_cell_scalar& flux_interface_efficace, Cut_field_double& cut_field_d_temperature)
{
  for (int phase = 0; phase < 2; phase++)
    {
      ajout_flux_interface_a_divergence_implementation(phase, flux_interface_efficace, cut_field_d_temperature);
    }
}

void ajout_flux_interface_a_divergence_implementation(int phase, const DoubleTabFT_cut_cell_scalar& flux_interface_efficace, Cut_field_double& cut_field_d_field)
{
  ajout_flux_interface_a_divergence_simple(phase, flux_interface_efficace, cut_field_d_field);
}

void ajout_flux_interface_a_divergence_simple(int phase, const DoubleTabFT_cut_cell_scalar& flux_interface_efficace, Cut_field_double& cut_field_d_field)
{
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_d_field.get_cut_cell_disc();
  DoubleTabFT_cut_cell& diph_d_field = (phase == 0) ? cut_field_d_field.diph_v_ : cut_field_d_field.diph_l_;

  for (int n = 0; n < cut_cell_disc.get_n_loc(); n++)
    {
      int sign = (phase == 0) ? +1 : -1;
      diph_d_field(n) += sign * flux_interface_efficace(n);
    }
}

void calcul_temperature_flux_interface(
  bool next_time,
  bool cut_cell,
  bool no_jump,
  const Cut_field_double& cut_field,
  const Facettes_Interp_FT& cut_cell_facettes_interpolation,
  const double ldal,
  const double ldav,
  DoubleTabFT& temperature_interp,
  DoubleTabFT& flux_normal_interp)
{
  double dist_1 = cut_cell_facettes_interpolation.get_distance_interpolation_1();
  double dist_2 = cut_cell_facettes_interpolation.get_distance_interpolation_2();

  const FixedVector<IntTabFT, 4>& interpolation_signed_independent_index = next_time ? cut_cell_facettes_interpolation.get_signed_independent_index_next() : cut_cell_facettes_interpolation.get_signed_independent_index_old();
  const FixedVector<DoubleTabFT, 4>& interpolation_coefficient = next_time ? cut_cell_facettes_interpolation.get_coefficient_next() : cut_cell_facettes_interpolation.get_coefficient_old();

  const Maillage_FT_IJK& maillage = next_time ? cut_cell_facettes_interpolation.maillage_ft_ijk() : cut_cell_facettes_interpolation.old_maillage_ft_ijk();
  const int nb_facettes = maillage.nb_facettes();
  const DoubleTab& normale_facettes = maillage.get_update_normale_facettes();

  int number_of_interpolation_points = cut_cell_facettes_interpolation.get_number_of_interpolation_points();

  temperature_interp.resize(nb_facettes);
  flux_normal_interp.resize(nb_facettes);
  if ((ldal + ldav) == 0.)
    {
      Cerr << "Corrige flux used with no conductivity. Ti and Qi set to 0. " << finl;
      temperature_interp = 0.;
      flux_normal_interp = 0.;
      return;
    }
  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    {
      double T_1v;
      double T_1l;
      double T_2v;
      double T_2l;

      if (!cut_cell)
        {
          Vecteur3 coords_fa7 = maillage.coords_fa7(fa7);
          Vecteur3 normal(normale_facettes, fa7);

          Vecteur3 coord_1v = Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(coords_fa7, normal, -dist_1);
          Vecteur3 coord_1l = Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(coords_fa7, normal, dist_1);

          T_1v = ijk_interpolate_skip_unknown_points(cut_field, coord_1v, 1.e31);
          T_1l = ijk_interpolate_skip_unknown_points(cut_field, coord_1l, 1.e31);

          if (number_of_interpolation_points > 1)
            {
              Vecteur3 coord_2v = Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(coords_fa7, normal, -dist_2);
              Vecteur3 coord_2l = Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(coords_fa7, normal, dist_2);

              T_2v = ijk_interpolate_skip_unknown_points(cut_field, coord_2v, 1.e31);
              T_2l = ijk_interpolate_skip_unknown_points(cut_field, coord_2l, 1.e31);
            }
          else
            {
              T_2v = 0.;
              T_2l = 0.;
            }
        }
      else
        {
          double field_0_liqu_1 = cut_field.from_signed_independent_index(interpolation_signed_independent_index[1](fa7, 0));
          double field_1_liqu_1 = cut_field.from_signed_independent_index(interpolation_signed_independent_index[1](fa7, 1));
          double field_2_liqu_1 = cut_field.from_signed_independent_index(interpolation_signed_independent_index[1](fa7, 2));
          double field_3_liqu_1 = cut_field.from_signed_independent_index(interpolation_signed_independent_index[1](fa7, 3));
          double field_0_vap_1  = cut_field.from_signed_independent_index(interpolation_signed_independent_index[0](fa7, 0));
          double field_1_vap_1  = cut_field.from_signed_independent_index(interpolation_signed_independent_index[0](fa7, 1));
          double field_2_vap_1  = cut_field.from_signed_independent_index(interpolation_signed_independent_index[0](fa7, 2));
          double field_3_vap_1  = cut_field.from_signed_independent_index(interpolation_signed_independent_index[0](fa7, 3));

          T_1l = field_0_liqu_1*interpolation_coefficient[1](fa7, 0) + field_1_liqu_1*interpolation_coefficient[1](fa7, 1) + field_2_liqu_1*interpolation_coefficient[1](fa7, 2) + field_3_liqu_1*interpolation_coefficient[1](fa7, 3);
          T_1v = field_0_vap_1*interpolation_coefficient[0](fa7, 0) + field_1_vap_1*interpolation_coefficient[0](fa7, 1) + field_2_vap_1*interpolation_coefficient[0](fa7, 2) + field_3_vap_1*interpolation_coefficient[0](fa7, 3);

          if (number_of_interpolation_points > 1)
            {
              double field_0_liqu_2 = cut_field.from_signed_independent_index(interpolation_signed_independent_index[3](fa7, 0));
              double field_1_liqu_2 = cut_field.from_signed_independent_index(interpolation_signed_independent_index[3](fa7, 1));
              double field_2_liqu_2 = cut_field.from_signed_independent_index(interpolation_signed_independent_index[3](fa7, 2));
              double field_3_liqu_2 = cut_field.from_signed_independent_index(interpolation_signed_independent_index[3](fa7, 3));
              double field_0_vap_2  = cut_field.from_signed_independent_index(interpolation_signed_independent_index[2](fa7, 0));
              double field_1_vap_2  = cut_field.from_signed_independent_index(interpolation_signed_independent_index[2](fa7, 1));
              double field_2_vap_2  = cut_field.from_signed_independent_index(interpolation_signed_independent_index[2](fa7, 2));
              double field_3_vap_2  = cut_field.from_signed_independent_index(interpolation_signed_independent_index[2](fa7, 3));

              T_2l = field_0_liqu_2*interpolation_coefficient[3](fa7, 0) + field_1_liqu_2*interpolation_coefficient[3](fa7, 1) + field_2_liqu_2*interpolation_coefficient[3](fa7, 2) + field_3_liqu_2*interpolation_coefficient[3](fa7, 3);
              T_2v = field_0_vap_2*interpolation_coefficient[2](fa7, 0) + field_1_vap_2*interpolation_coefficient[2](fa7, 1) + field_2_vap_2*interpolation_coefficient[2](fa7, 2) + field_3_vap_2*interpolation_coefficient[2](fa7, 3);
            }
          else
            {
              T_2v = 0.;
              T_2l = 0.;
            }
        }


      if (!no_jump)
        {
          if (number_of_interpolation_points == 1)
            {
              const double Ti = (T_1l * ldal + T_1v * ldav) / (ldal + ldav);

              temperature_interp(fa7) = Ti;
              flux_normal_interp(fa7) = ldav * (Ti - T_1v) / dist_1;
            }
          else if (number_of_interpolation_points == 2)
            {
              const double Ti = (ldal*T_1l/dist_1 + ldav*T_1v/dist_1 + ldal*T_2l/dist_2 + ldav*T_2v/dist_2 + ldal*(T_1l - T_2l)/(dist_2 - dist_1) - ldav*(T_2v - T_1v)/(dist_2 - dist_1))/(ldal/dist_1 + ldav/dist_1 + ldal/dist_2 + ldav/dist_2);

              temperature_interp(fa7) = Ti;
              flux_normal_interp(fa7) = - ldav * ((T_1v - Ti)/dist_1 + (T_2v - Ti)/dist_2 - (T_2v - T_1v)/(dist_2 - dist_1));
            }
          else
            {
              Cerr << "Nombre de points d'interpolation non reconnu." << finl;
              Process::exit();
            }
        }
      else
        {
          if (number_of_interpolation_points == 1)
            {
              const double Ti = .5*(T_1l + T_1v);

              temperature_interp(fa7) = Ti;
              flux_normal_interp(fa7) = 2./(1./ldav + 1./ldal) * (Ti - T_1v) / dist_1;
            }
          else
            {
              Cerr << "Nombre de points d'interpolation non reconnu pour le cas 'no jump'." << finl;
              Process::exit();
            }
        }
    }
}

