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
// File      : Cut_cell_convection_auxiliaire.cpp
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#include <Cut_cell_convection_auxiliaire.h>
#include <Cut_cell_FT_Disc.h>
#include <IJK_Thermal_base.h>
#include <Process.h>
#include <IJK_FT_cut_cell.h>
#include <stat_counters.h>
#include <IJK_Navier_Stokes_tools.h>
#include <IJK_Navier_Stokes_tools_cut_cell.h>
#include <Param.h>

Implemente_instanciable_sans_constructeur(Cut_cell_convection_auxiliaire, "Cut_cell_convection_auxiliaire", Cut_cell_schema_auxiliaire) ;

Cut_cell_convection_auxiliaire::Cut_cell_convection_auxiliaire()
{
  methode_temperature_remplissage_ = METHODE_TEMPERATURE_REMPLISSAGE::NON_INITIALISE;
  correction_petites_cellules_ = CORRECTION_PETITES_CELLULES::CORRECTION_SYMETRIQUE;

  no_static_update_ = true;
}

Sortie& Cut_cell_convection_auxiliaire::printOn(Sortie& os) const
{
  Objet_U::printOn(os);
  return os;
}

Entree& Cut_cell_convection_auxiliaire::readOn(Entree& is)
{
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades(is);
  return is;
}

void Cut_cell_convection_auxiliaire::set_param(Param& param)
{
  Cut_cell_schema_auxiliaire::set_param(param);
}

double Cut_cell_convection_auxiliaire::dying_cells_flux(int num_face, int phase, int n, const FixedVector<Cut_field_scalar, 3>& cut_field_total_velocity, const Cut_field_scalar& cut_field_temperature)
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
  int di_decale = sign*(dir == 0);
  int dj_decale = sign*(dir == 1);
  int dk_decale = sign*(dir == 2);

  double next_indicatrice_decale = cut_cell_disc.get_interfaces().In(i+di_decale,j+dj_decale,k+dk_decale);

  double temperature_centre = (phase == 0) ? cut_field_temperature.diph_v_(n) : cut_field_temperature.diph_l_(n);

  double temperature_decale = (phase == ((int)next_indicatrice_decale))*cut_field_temperature.pure_(i+di_decale,j+dj_decale,k+dk_decale);
  int n_decale = cut_cell_disc.get_n(i+di_decale, j+dj_decale, k+dk_decale);
  if (n_decale >= 0)
    {
      temperature_decale = (phase == 0) ? cut_field_temperature.diph_v_(n_decale) : cut_field_temperature.diph_l_(n_decale);
    }

  int n_face = cut_cell_disc.get_n_face(num_face, n, i, j, k);
  if (n_face >= 0)
    {
      double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir) : cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir);
      if (surface_efficace > 0)
        {
          double total_velocity = (phase == 0) ? cut_field_total_velocity[dir].diph_v_(n_face) : cut_field_total_velocity[dir].diph_l_(n_face);
          double temperature = (sign*total_velocity < 0) ? temperature_decale : temperature_centre;
          return -sign*surface_efficace*temperature*total_velocity;
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
          double total_velocity = cut_field_total_velocity[dir].pure_(i+di,j+dj,k+dk);
          double temperature = (sign*total_velocity < 0) ? temperature_decale : temperature_centre;
          return -sign*surface_efficace*temperature*total_velocity;
        }
      else
        {
          return 0.;
        }
    }
}

double Cut_cell_convection_auxiliaire::small_nascent_cells_flux(int num_face, int phase, int n, const FixedVector<Cut_field_scalar, 3>& cut_field_total_velocity, const Cut_field_scalar& cut_field_temperature)
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
  int di_decale = sign*(dir == 0);
  int dj_decale = sign*(dir == 1);
  int dk_decale = sign*(dir == 2);

  double next_indicatrice_decale = cut_cell_disc.get_interfaces().In(i+di_decale,j+dj_decale,k+dk_decale);

  double temperature_decale = (phase == ((int)next_indicatrice_decale))*cut_field_temperature.pure_(i+di_decale,j+dj_decale,k+dk_decale);
  int n_decale = cut_cell_disc.get_n(i+di_decale, j+dj_decale, k+dk_decale);
  if (n_decale >= 0)
    {
      temperature_decale = (phase == 0) ? cut_field_temperature.diph_v_(n_decale) : cut_field_temperature.diph_l_(n_decale);
    }

  int n_face = cut_cell_disc.get_n_face(num_face, n, i, j, k);
  if (n_face >= 0)
    {
      double surface_efficace = (phase == 0) ? 1 - cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir) : cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face()(n_face, dir);
      if (surface_efficace > 0)
        {
          double total_velocity = (phase == 0) ? cut_field_total_velocity[dir].diph_v_(n_face) : cut_field_total_velocity[dir].diph_l_(n_face);
          double temperature = temperature_decale;
          return -sign*surface_efficace*temperature*total_velocity;
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
          double total_velocity = cut_field_total_velocity[dir].pure_(i+di,j+dj,k+dk);
          double temperature = temperature_decale;
          return -sign*surface_efficace*temperature*total_velocity;
        }
      else
        {
          return 0.;
        }
    }
}


void Cut_cell_convection_auxiliaire::calcule_temperature_face_depuis_centre(double lambda_liquid, double lambda_vapour, const IJK_Field_double& flux_interface_ns, const Cut_field_scalar& cut_field_temperature, FixedVector<FixedVector<IJK_Field_double, 3>, 2>& temperature_face)
{
  const bool next_time = false;

  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();
  const IJK_Field_double& surface_interface = next_time ? cut_cell_disc.get_interfaces().get_surface_interface_next() : cut_cell_disc.get_interfaces().get_surface_interface_old();

  const FixedVector<IJK_Field_double, 3>& old_indicatrice_surfacique_face  = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_face_old();
  const FixedVector<IJK_Field_double, 3>& next_indicatrice_surfacique_face = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_face_next();

  double dx = cut_cell_disc.get_splitting().get_grid_geometry().get_constant_delta(0);
  double dy = cut_cell_disc.get_splitting().get_grid_geometry().get_constant_delta(1);
  double dz = cut_cell_disc.get_splitting().get_grid_geometry().get_constant_delta(2);

  temperature_face[0][0].data()=0;
  temperature_face[0][1].data()=0;
  temperature_face[0][2].data()=0;
  temperature_face[1][0].data()=0;
  temperature_face[1][1].data()=0;
  temperature_face[1][2].data()=0;

  for (int n = 0; n < cut_cell_disc.get_n_tot(); n++)
    {
      Int3 ijk = cut_cell_disc.get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      for (int dir = 0; dir < 3; dir++)
        {
          for (int phase = 0; phase < 2; phase++)
            {
              double old_indicatrice  = cut_cell_disc.get_interfaces().I(i,j,k);
              double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);

              double normal_x = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,0);
              double normal_y = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,1);
              double normal_z = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n,2);

              double centre_bary_x = dx*cut_cell_disc.get_interfaces().get_barycentre(next_time, 0, phase, i,j,k, old_indicatrice, next_indicatrice);
              double centre_bary_y = dy*cut_cell_disc.get_interfaces().get_barycentre(next_time, 1, phase, i,j,k, old_indicatrice, next_indicatrice);
              double centre_bary_z = dz*cut_cell_disc.get_interfaces().get_barycentre(next_time, 2, phase, i,j,k, old_indicatrice, next_indicatrice);

              double old_indicatrice_surfacique  = old_indicatrice_surfacique_face[dir](i,j,k);
              double next_indicatrice_surfacique = next_indicatrice_surfacique_face[dir](i,j,k);

              double face_bary_x = dx*cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 0, phase, i,j,k, old_indicatrice_surfacique, next_indicatrice_surfacique);
              double face_bary_y = dy*cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 1, phase, i,j,k, old_indicatrice_surfacique, next_indicatrice_surfacique);
              double face_bary_z = dz*cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 2, phase, i,j,k, old_indicatrice_surfacique, next_indicatrice_surfacique);

              double temperature_centre = (phase == 0) ? cut_field_temperature.diph_v_(n) : cut_field_temperature.diph_l_(n);

              assert((flux_interface_ns(i,j,k) == 0.) || (surface_interface(i,j,k) != 0.)); // La surface ne devrait pas etre nulle si il y a un flux
              double lambda = (phase == 0) ? lambda_vapour : lambda_liquid;
              double dTdn = (flux_interface_ns(i,j,k) == 0.) ? 0. : flux_interface_ns(i,j,k)/(lambda * surface_interface(i,j,k));
              assert((flux_interface_ns(i,j,k) == 0.) || (surface_interface(i,j,k) != 0));

              const DoubleTabFT_cut_cell_vector3& indicatrice_surfacique = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face();
              double indicatrice_surface = ((phase == 0) ? 1 - indicatrice_surfacique(n,dir) : indicatrice_surfacique(n,dir));

              temperature_face[phase][dir](i,j,k) = (indicatrice_surface == 0.) ? 0. : (temperature_centre + dTdn * ((face_bary_x - centre_bary_x)*normal_x + (face_bary_y - centre_bary_y)*normal_y + (face_bary_z - centre_bary_z)*normal_z));
            }
        }
    }
}

void Cut_cell_convection_auxiliaire::calcule_temperature_face_depuis_facette(double lambda_liquid, double lambda_vapour, const ArrOfDouble& interfacial_temperature, const ArrOfDouble& interfacial_phin_ai, const Cut_field_scalar& cut_field_temperature, REF(IJK_FT_cut_cell)& ref_ijk_ft, FixedVector<FixedVector<IJK_Field_double, 3>, 2>& temperature_face_ft, FixedVector<FixedVector<IJK_Field_double, 3>, 2>& temperature_face_ns)
{
  const bool next_time = false;

  // Calcul de T_face sur le maillage FT
  {
    const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();
    const Maillage_FT_IJK& mesh = cut_cell_disc.get_interfaces().maillage_ft_ijk();
    const IntTab& facettes = next_time ? mesh.facettes() : mesh.facettes_old();
    const DoubleTab& sommets = next_time ? mesh.sommets() : mesh.sommets_old();
    const DoubleTab& normale_facettes = next_time ? mesh.get_update_normale_facettes() : mesh.get_normale_facettes_old();
    const ArrOfDouble& surface_facettes = next_time ? mesh.get_update_surface_facettes() : mesh.get_surface_facettes_old();
    const Intersections_Elem_Facettes& intersec = next_time ? mesh.intersections_elem_facettes() : mesh.intersections_elem_facettes_old();
    const IJK_Splitting& s = temperature_face_ft.get_splitting();

    double dx = s.get_grid_geometry().get_constant_delta(DIRECTION_I);
    double dy = s.get_grid_geometry().get_constant_delta(DIRECTION_J);
    double dz = s.get_grid_geometry().get_constant_delta(DIRECTION_K);

    const IJK_Grid_Geometry& geom = s.get_grid_geometry();
    double origin_x = geom.get_origin(DIRECTION_I);
    double origin_y = geom.get_origin(DIRECTION_J);
    double origin_z = geom.get_origin(DIRECTION_K);

    int offset_x = s.get_offset_local(DIRECTION_I);
    int offset_y = s.get_offset_local(DIRECTION_J);
    int offset_z = s.get_offset_local(DIRECTION_K);

    const FixedVector<IJK_Field_double, 3>& old_indicatrice_surfacique_face  = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_face_old();
    const FixedVector<IJK_Field_double, 3>& next_indicatrice_surfacique_face = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_face_next();

    const int ni_ft = temperature_face_ft[0][0].ni();
    const int nj_ft = temperature_face_ft[0][0].nj();
    const int nk_ft = temperature_face_ft[0][0].nk();

    // Initialisation
    {
      for (int k_ft = 0; k_ft < nk_ft; k_ft++)
        {
          for (int j_ft = 0; j_ft < nj_ft; j_ft++)
            {
              for (int i_ft = 0; i_ft < ni_ft; i_ft++)
                {
                  for (int dir = 0; dir < 3; dir++)
                    {
                      for (int phase = 0; phase < 2; phase++)
                        {
                          temperature_face_ft[phase][dir](i_ft, j_ft, k_ft) = 0.;
                        }
                    }
                }
            }
        }
    }

    // Calcul pour les faces coupees par l'interface
    {
      const ArrOfInt& index_elem = intersec.index_elem();
      for (int k_ft = 0; k_ft < nk_ft; k_ft++)
        {
          for (int j_ft = 0; j_ft < nj_ft; j_ft++)
            {
              for (int i_ft = 0; i_ft < ni_ft; i_ft++)
                {
                  if (next_time ? (!cut_cell_disc.get_interfaces().est_pure(cut_cell_disc.get_interfaces().In_ft(i_ft,j_ft,k_ft))) : (!cut_cell_disc.get_interfaces().est_pure(cut_cell_disc.get_interfaces().I_ft(i_ft,j_ft,k_ft))))
                    {
                      double x_centre = (i_ft + offset_x + .5)*dx + origin_x;
                      double y_centre = (j_ft + offset_y + .5)*dy + origin_y;
                      double z_centre = (k_ft + offset_z + .5)*dz + origin_z;
                      int i_ns = cut_cell_disc.get_i_selon_dir(0, x_centre);
                      int j_ns = cut_cell_disc.get_i_selon_dir(1, y_centre);
                      int k_ns = cut_cell_disc.get_i_selon_dir(2, z_centre);

                      int n = cut_cell_disc.get_n(i_ns, j_ns, k_ns);

                      assert(mesh.ref_splitting().valeur() == s);
                      const int num_elem = s.convert_ijk_cell_to_packed(i_ft,j_ft,k_ft);
                      int index = index_elem[num_elem];

                      int facette_la_plus_proche[2][3] = {{-1, -1, -1}, {-1, -1, -1}};
                      double min_distance_facette[2][3] = {{DMAXFLOAT, DMAXFLOAT, DMAXFLOAT}, {DMAXFLOAT, DMAXFLOAT, DMAXFLOAT}};

                      // Boucle sur les facettes qui traversent cet element
                      while (index >= 0)
                        {
                          const Intersections_Elem_Facettes_Data& data = intersec.data_intersection(index);
                          const int fa7 = data.numero_facette_;

                          double coord_facettes[3] = {0};
                          for (int dir = 0; dir < 3; dir++)
                            for (int som = 0; som < 3; som++)
                              coord_facettes[dir] += sommets(facettes(fa7, som), dir);
                          for (int dir = 0; dir < 3; dir++)
                            coord_facettes[dir] /= 3.;

                          for (int dir = 0; dir < 3; dir++)
                            {
                              for (int phase = 0; phase < 2; phase++)
                                {
                                  double old_indicatrice_surfacique  = old_indicatrice_surfacique_face[dir](i_ns,j_ns,k_ns);
                                  double next_indicatrice_surfacique = next_indicatrice_surfacique_face[dir](i_ns,j_ns,k_ns);
                                  double face_bary_x = origin_x + dx*(i_ft + offset_x + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 0, phase, i_ns,j_ns,k_ns, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                  double face_bary_y = origin_y + dy*(j_ft + offset_y + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 1, phase, i_ns,j_ns,k_ns, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                  double face_bary_z = origin_z + dz*(k_ft + offset_z + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 2, phase, i_ns,j_ns,k_ns, old_indicatrice_surfacique, next_indicatrice_surfacique)));

                                  double distance_facette = sqrt((face_bary_x - coord_facettes[0])*(face_bary_x - coord_facettes[0]) + (face_bary_y - coord_facettes[1])*(face_bary_y - coord_facettes[1]) + (face_bary_z - coord_facettes[2])*(face_bary_z - coord_facettes[2]));
                                  // Ce code est un debut d'implementation. Le calcul de la distance a la facette n'est pas termine (distance au centre uniquement).
                                  Cerr << "distance_facette : calcul pas implementee." << finl;
                                  Process::exit();
                                  facette_la_plus_proche[phase][dir] = distance_facette < min_distance_facette[phase][dir] ? fa7 : facette_la_plus_proche[phase][dir];
                                  min_distance_facette[phase][dir] = distance_facette < min_distance_facette[phase][dir] ? distance_facette : min_distance_facette[phase][dir];
                                }
                            }

                          index = data.index_facette_suivante_;
                        };

                      for (int dir = 0; dir < 3; dir++)
                        {
                          for (int phase = 0; phase < 2; phase++)
                            {
                              double coord_facettes[3] = {0};
                              for (int dir2 = 0; dir2 < 3; dir2++)
                                for (int som = 0; som < 3; som++)
                                  coord_facettes[dir2] += sommets(facettes(facette_la_plus_proche[phase][dir2], som), dir2);
                              for (int dir2 = 0; dir2 < 3; dir2++)
                                coord_facettes[dir2] /= 3.;


                              double old_indicatrice_surfacique  = old_indicatrice_surfacique_face[dir](i_ns,j_ns,k_ns);
                              double next_indicatrice_surfacique = next_indicatrice_surfacique_face[dir](i_ns,j_ns,k_ns);
                              double face_bary_x = origin_x + dx*(i_ft + offset_x + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 0, phase, i_ns,j_ns,k_ns, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                              double face_bary_y = origin_y + dy*(j_ft + offset_y + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 1, phase, i_ns,j_ns,k_ns, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                              double face_bary_z = origin_z + dz*(k_ft + offset_z + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 2, phase, i_ns,j_ns,k_ns, old_indicatrice_surfacique, next_indicatrice_surfacique)));

                              const DoubleTabFT_cut_cell_vector3& indicatrice_surfacique = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face();
                              double indicatrice_surface = ((phase == 0) ? 1 - indicatrice_surfacique(n,dir) : indicatrice_surfacique(n,dir));

                              double lambda = (phase == 0) ? lambda_vapour : lambda_liquid;
                              double dTdn = interfacial_phin_ai(facette_la_plus_proche[phase][dir])/(lambda * surface_facettes(facette_la_plus_proche[phase][dir]));

                              double normal_x = normale_facettes(facette_la_plus_proche[phase][dir], 0);
                              double normal_y = normale_facettes(facette_la_plus_proche[phase][dir], 1);
                              double normal_z = normale_facettes(facette_la_plus_proche[phase][dir], 2);
                              double norm_normal = sqrt(normal_x*normal_x + normal_y*normal_y + normal_z*normal_z);
                              normal_x /= norm_normal;
                              normal_y /= norm_normal;
                              normal_z /= norm_normal;
                              temperature_face_ft[phase][dir](i_ft,j_ft,k_ft) = (indicatrice_surface == 0.) ? 0. : (interfacial_temperature(facette_la_plus_proche[phase][dir])/surface_facettes(facette_la_plus_proche[phase][dir]) + dTdn * ((face_bary_x - coord_facettes[0])*normal_x + (face_bary_y - coord_facettes[1])*normal_y + (face_bary_z - coord_facettes[2])*normal_z));
                            }
                        }
                    }
                }
            }
        }
    }
  }

  for (int dir = 0; dir < 3; dir++)
    {
      for (int phase = 0; phase < 2; phase++)
        {
          temperature_face_ft[phase][dir].echange_espace_virtuel(temperature_face_ft[phase][dir].ghost());

          // Calcul du flux interface sur le domaine NS :
          ref_ijk_ft->redistrib_from_ft_elem().redistribute(temperature_face_ft[phase][dir], temperature_face_ns[phase][dir]);
          temperature_face_ns[phase][dir].echange_espace_virtuel(temperature_face_ns[phase][dir].ghost());
        }
    }
}

void Cut_cell_convection_auxiliaire::calcule_temperature_face_depuis_facette_interpolate(CUT_CELL_CONV_FACE_INTERPOLATION face_interp, double timestep, const ArrOfDouble& interfacial_temperature, const IJK_Field_double& temperature_ft, const Cut_field_scalar& cut_field_temperature, REF(IJK_FT_cut_cell)& ref_ijk_ft, FixedVector<FixedVector<IJK_Field_double, 3>, 2>& temperature_face_ft, FixedVector<FixedVector<IJK_Field_double, 3>, 2>& temperature_face_ns, const FixedVector<Cut_field_scalar, 3>& cut_field_total_velocity)
{
  const bool next_time = false;

  int count_status_not_found = 0;
  int count_status_below_10 = 0;
  int count_status_below_100 = 0;
  int count_status_below_1000 = 0;
  int count_status_below_10000 = 0;
  int count_status_below_100000 = 0;
  int count_status_below_1000000 = 0;
  int count_status_above_1000000 = 0;

  // Calcul de T_face sur le maillage FT
  {
    const Cut_cell_FT_Disc& cut_cell_disc = cut_field_temperature.get_cut_cell_disc();
    const IJK_Splitting& s = temperature_face_ft.get_splitting();

    double dx = s.get_grid_geometry().get_constant_delta(DIRECTION_I);
    double dy = s.get_grid_geometry().get_constant_delta(DIRECTION_J);
    double dz = s.get_grid_geometry().get_constant_delta(DIRECTION_K);

    const IJK_Grid_Geometry& geom = s.get_grid_geometry();
    double origin_x = geom.get_origin(DIRECTION_I);
    double origin_y = geom.get_origin(DIRECTION_J);
    double origin_z = geom.get_origin(DIRECTION_K);

    int offset_x = s.get_offset_local(DIRECTION_I);
    int offset_y = s.get_offset_local(DIRECTION_J);
    int offset_z = s.get_offset_local(DIRECTION_K);

    const FixedVector<IJK_Field_double, 3>& old_indicatrice_surfacique_face  = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_face_old();
    const FixedVector<IJK_Field_double, 3>& next_indicatrice_surfacique_face = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_face_next();

    const int ni_ft = temperature_face_ft[0][0].ni();
    const int nj_ft = temperature_face_ft[0][0].nj();
    const int nk_ft = temperature_face_ft[0][0].nk();

    // Initialisation
    {
      for (int k_ft = 0; k_ft < nk_ft; k_ft++)
        {
          for (int j_ft = 0; j_ft < nj_ft; j_ft++)
            {
              for (int i_ft = 0; i_ft < ni_ft; i_ft++)
                {
                  for (int dir = 0; dir < 3; dir++)
                    {
                      for (int phase = 0; phase < 2; phase++)
                        {
                          temperature_face_ft[phase][dir](i_ft, j_ft, k_ft) = 0.;
                        }
                    }
                }
            }
        }
    }

    // Calcul pour les faces coupees par l'interface
    {
      for (int k_ft = 0; k_ft < nk_ft; k_ft++)
        {
          for (int j_ft = 0; j_ft < nj_ft; j_ft++)
            {
              for (int i_ft = 0; i_ft < ni_ft; i_ft++)
                {
                  if (next_time ? (!cut_cell_disc.get_interfaces().est_pure(cut_cell_disc.get_interfaces().In_ft(i_ft,j_ft,k_ft))) : (!cut_cell_disc.get_interfaces().est_pure(cut_cell_disc.get_interfaces().I_ft(i_ft,j_ft,k_ft))))
                    {
                      double x_centre = (i_ft + offset_x + .5)*dx + origin_x;
                      double y_centre = (j_ft + offset_y + .5)*dy + origin_y;
                      double z_centre = (k_ft + offset_z + .5)*dz + origin_z;
                      int i_ns = cut_cell_disc.get_i_selon_dir(0, x_centre);
                      int j_ns = cut_cell_disc.get_i_selon_dir(1, y_centre);
                      int k_ns = cut_cell_disc.get_i_selon_dir(2, z_centre);

                      if (!cut_cell_disc.within_ghost(i_ns, j_ns, k_ns, 1, 1))
                        continue;

                      int n_centre = cut_cell_disc.get_n(i_ns, j_ns, k_ns);

                      for (int num_face = 0; num_face < 6; num_face++)
                        {
                          int dir = num_face%3;
                          int decalage = num_face/3;
                          int sign = decalage*2 -1;

                          int di = decalage*(dir == 0);
                          int dj = decalage*(dir == 1);
                          int dk = decalage*(dir == 2);
                          int di_decale = sign*(dir == 0);
                          int dj_decale = sign*(dir == 1);
                          int dk_decale = sign*(dir == 2);

                          int n = cut_cell_disc.get_n(i_ns+di, j_ns+dj, k_ns+dk);
                          int n_decale = cut_cell_disc.get_n(i_ns+di_decale, j_ns+dj_decale, k_ns+dk_decale);

                          double normal_x = 0.;
                          double normal_y = 0.;
                          double normal_z = 0.;
                          if (n_centre >= 0 && n_decale >= 0)
                            {
                              double normal_x_centre = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_centre,0);
                              double normal_y_centre = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_centre,1);
                              double normal_z_centre = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_centre,2);
                              double normal_x_decale = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_decale,0);
                              double normal_y_decale = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_decale,1);
                              double normal_z_decale = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_decale,2);
                              normal_x = .5*(normal_x_centre + normal_x_decale);
                              normal_y = .5*(normal_y_centre + normal_y_decale);
                              normal_z = .5*(normal_z_centre + normal_z_decale);
                            }
                          else if (n_centre >= 0)
                            {
                              double normal_x_centre = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_centre,0);
                              double normal_y_centre = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_centre,1);
                              double normal_z_centre = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_centre,2);
                              normal_x = .5*(normal_x_centre);
                              normal_y = .5*(normal_y_centre);
                              normal_z = .5*(normal_z_centre);
                            }
                          else if (n_decale >= 0)
                            {
                              double normal_x_decale = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_decale,0);
                              double normal_y_decale = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_decale,1);
                              double normal_z_decale = cut_cell_disc.get_interfaces().get_normale_deplacement_interface()(n_decale,2);
                              normal_x = .5*(normal_x_decale);
                              normal_y = .5*(normal_y_decale);
                              normal_z = .5*(normal_z_decale);
                            }
                          else
                            {
                              Cerr << "Cut_cell_convection_auxiliaire::calcule_temperature_face_depuis_facette_interpolate: Les deux n sont pures ?" << finl;
                              Process::exit();
                            }
                          double norm_normal = sqrt(normal_x*normal_x + normal_y*normal_y + normal_z*normal_z);
                          normal_x /= norm_normal;
                          normal_y /= norm_normal;
                          normal_z /= norm_normal;


                          for (int phase = 0; phase < 2; phase++)
                            {
                              //int n = cut_cell_disc.get_n(i_ns+di, j_ns+dj, k_ns+dk);
                              //const DoubleTabFT_cut_cell_vector3& indicatrice_surfacique = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face();
                              //double indicatrice_surface = ((phase == 0) ? 1 - indicatrice_surfacique(n,dir) : indicatrice_surfacique(n,dir));
                              double old_indicatrice_surfacique  = old_indicatrice_surfacique_face[dir](i_ns+di,j_ns+dj,k_ns+dk);
                              double next_indicatrice_surfacique = next_indicatrice_surfacique_face[dir](i_ns+di,j_ns+dj,k_ns+dk);
                              const double indicatrice_surfacique = next_time ? next_indicatrice_surfacique : old_indicatrice_surfacique;
                              double indicatrice_surface = ((phase == 0) ? 1 - indicatrice_surfacique : indicatrice_surfacique);
                              if (indicatrice_surface == 0.)
                                {
                                  continue; // Pas utile de determiner une temperature
                                }
                              else
                                {
                                  double coordinates[3] = {};

                                  if (face_interp == CUT_CELL_CONV_FACE_INTERPOLATION::INTERPOLATE)
                                    {
                                      double face_bary_x = origin_x + dx*(i_ft + di + offset_x + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 0, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_y = origin_y + dy*(j_ft + dj + offset_y + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 1, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_z = origin_z + dz*(k_ft + dk + offset_z + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 2, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));

                                      coordinates[0] = face_bary_x;
                                      coordinates[1] = face_bary_y;
                                      coordinates[2] = face_bary_z;
                                    }
                                  else if (face_interp == CUT_CELL_CONV_FACE_INTERPOLATION::INTERPOLATEAMONT0)
                                    {
                                      double v_x = (n < 0) ? cut_field_total_velocity[0].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[0].diph_v_(n) : cut_field_total_velocity[0].diph_l_(n));
                                      double v_y = (n < 0) ? cut_field_total_velocity[1].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[1].diph_v_(n) : cut_field_total_velocity[1].diph_l_(n));
                                      double v_z = (n < 0) ? cut_field_total_velocity[2].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[2].diph_v_(n) : cut_field_total_velocity[2].diph_l_(n));
                                      double v_dir = select(dir, v_x, v_y, v_z);

                                      int i_ft_amont = i_ft + di - (dir == 0)*(v_dir>0);
                                      int j_ft_amont = j_ft + dj - (dir == 1)*(v_dir>0);
                                      int k_ft_amont = k_ft + dk - (dir == 2)*(v_dir>0);

                                      int i_ns_amont = i_ns + di - (dir == 0)*(v_dir>0);
                                      int j_ns_amont = j_ns + dj - (dir == 1)*(v_dir>0);
                                      int k_ns_amont = k_ns + dk - (dir == 2)*(v_dir>0);

                                      double old_indicatrice_amont  = cut_cell_disc.get_interfaces().I(i_ns_amont,j_ns_amont,k_ns_amont);
                                      double next_indicatrice_amont = cut_cell_disc.get_interfaces().In(i_ns_amont,j_ns_amont,k_ns_amont);

                                      double amont_bary_x = origin_x + dx*(i_ft_amont + offset_x + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 0, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_y = origin_y + dy*(j_ft_amont + offset_y + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 1, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_z = origin_z + dz*(k_ft_amont + offset_z + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 2, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));

                                      coordinates[0] = amont_bary_x;
                                      coordinates[1] = amont_bary_y;
                                      coordinates[2] = amont_bary_z;
                                    }
                                  else if (face_interp == CUT_CELL_CONV_FACE_INTERPOLATION::INTERPOLATEAMONT1)
                                    {
                                      double v_x = (n < 0) ? cut_field_total_velocity[0].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[0].diph_v_(n) : cut_field_total_velocity[0].diph_l_(n));
                                      double v_y = (n < 0) ? cut_field_total_velocity[1].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[1].diph_v_(n) : cut_field_total_velocity[1].diph_l_(n));
                                      double v_z = (n < 0) ? cut_field_total_velocity[2].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[2].diph_v_(n) : cut_field_total_velocity[2].diph_l_(n));
                                      double v_dir = select(dir, v_x, v_y, v_z);

                                      int i_ft_amont = i_ft + di - (dir == 0)*(v_dir>0);
                                      int j_ft_amont = j_ft + dj - (dir == 1)*(v_dir>0);
                                      int k_ft_amont = k_ft + dk - (dir == 2)*(v_dir>0);

                                      int i_ns_amont = i_ns + di - (dir == 0)*(v_dir>0);
                                      int j_ns_amont = j_ns + dj - (dir == 1)*(v_dir>0);
                                      int k_ns_amont = k_ns + dk - (dir == 2)*(v_dir>0);

                                      double old_indicatrice_amont  = cut_cell_disc.get_interfaces().I(i_ns_amont,j_ns_amont,k_ns_amont);
                                      double next_indicatrice_amont = cut_cell_disc.get_interfaces().In(i_ns_amont,j_ns_amont,k_ns_amont);

                                      double amont_bary_x = origin_x + dx*(i_ft_amont + offset_x + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 0, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_y = origin_y + dy*(j_ft_amont + offset_y + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 1, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_z = origin_z + dz*(k_ft_amont + offset_z + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 2, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));

                                      double face_bary_x = origin_x + dx*(i_ft + di + offset_x + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 0, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_y = origin_y + dy*(j_ft + dj + offset_y + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 1, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_z = origin_z + dz*(k_ft + dk + offset_z + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 2, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));

                                      coordinates[0] = (dir==0) ? amont_bary_x : face_bary_x;
                                      coordinates[1] = (dir==1) ? amont_bary_y : face_bary_y;
                                      coordinates[2] = (dir==2) ? amont_bary_z : face_bary_z;
                                    }
                                  else if (face_interp == CUT_CELL_CONV_FACE_INTERPOLATION::INTERPOLATEAMONT2)
                                    {
                                      double v_x = (n < 0) ? cut_field_total_velocity[0].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[0].diph_v_(n) : cut_field_total_velocity[0].diph_l_(n));
                                      double v_y = (n < 0) ? cut_field_total_velocity[1].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[1].diph_v_(n) : cut_field_total_velocity[1].diph_l_(n));
                                      double v_z = (n < 0) ? cut_field_total_velocity[2].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[2].diph_v_(n) : cut_field_total_velocity[2].diph_l_(n));

                                      double face_bary_x = origin_x + dx*(i_ft + di + offset_x + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 0, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_y = origin_y + dy*(j_ft + dj + offset_y + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 1, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_z = origin_z + dz*(k_ft + dk + offset_z + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 2, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));

                                      coordinates[0] = face_bary_x - .5*v_x*timestep;
                                      coordinates[1] = face_bary_y - .5*v_y*timestep;
                                      coordinates[2] = face_bary_z - .5*v_z*timestep;
                                    }
                                  else if (face_interp == CUT_CELL_CONV_FACE_INTERPOLATION::INTERPOLATENORMALE0)
                                    {
                                      double v_x = (n < 0) ? cut_field_total_velocity[0].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[0].diph_v_(n) : cut_field_total_velocity[0].diph_l_(n));
                                      double v_y = (n < 0) ? cut_field_total_velocity[1].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[1].diph_v_(n) : cut_field_total_velocity[1].diph_l_(n));
                                      double v_z = (n < 0) ? cut_field_total_velocity[2].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[2].diph_v_(n) : cut_field_total_velocity[2].diph_l_(n));
                                      double v_dir = select(dir, v_x, v_y, v_z);

                                      int i_ft_amont = i_ft + di - (dir == 0)*(v_dir>0);
                                      int j_ft_amont = j_ft + dj - (dir == 1)*(v_dir>0);
                                      int k_ft_amont = k_ft + dk - (dir == 2)*(v_dir>0);

                                      int i_ns_amont = i_ns + di - (dir == 0)*(v_dir>0);
                                      int j_ns_amont = j_ns + dj - (dir == 1)*(v_dir>0);
                                      int k_ns_amont = k_ns + dk - (dir == 2)*(v_dir>0);

                                      double old_indicatrice_amont  = cut_cell_disc.get_interfaces().I(i_ns_amont,j_ns_amont,k_ns_amont);
                                      double next_indicatrice_amont = cut_cell_disc.get_interfaces().In(i_ns_amont,j_ns_amont,k_ns_amont);

                                      double amont_bary_x = origin_x + dx*(i_ft_amont + offset_x + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 0, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_y = origin_y + dy*(j_ft_amont + offset_y + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 1, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_z = origin_z + dz*(k_ft_amont + offset_z + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 2, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));

                                      double face_bary_x = origin_x + dx*(i_ft + di + offset_x + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 0, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_y = origin_y + dy*(j_ft + dj + offset_y + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 1, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_z = origin_z + dz*(k_ft + dk + offset_z + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 2, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));

                                      double normal_dir = select(dir, normal_x, normal_y, normal_z);
                                      bool case_perpendicular = (std::abs(normal_dir) > 0.72);

                                      coordinates[0] = case_perpendicular ? amont_bary_x : face_bary_x;
                                      coordinates[1] = case_perpendicular ? amont_bary_y : face_bary_y;
                                      coordinates[2] = case_perpendicular ? amont_bary_z : face_bary_z;
                                    }
                                  else if (face_interp == CUT_CELL_CONV_FACE_INTERPOLATION::INTERPOLATENORMALE1)
                                    {
                                      double v_x = (n < 0) ? cut_field_total_velocity[0].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[0].diph_v_(n) : cut_field_total_velocity[0].diph_l_(n));
                                      double v_y = (n < 0) ? cut_field_total_velocity[1].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[1].diph_v_(n) : cut_field_total_velocity[1].diph_l_(n));
                                      double v_z = (n < 0) ? cut_field_total_velocity[2].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[2].diph_v_(n) : cut_field_total_velocity[2].diph_l_(n));
                                      double v_dir = select(dir, v_x, v_y, v_z);

                                      int i_ft_amont = i_ft + di - (dir == 0)*(v_dir>0);
                                      int j_ft_amont = j_ft + dj - (dir == 1)*(v_dir>0);
                                      int k_ft_amont = k_ft + dk - (dir == 2)*(v_dir>0);

                                      int i_ns_amont = i_ns + di - (dir == 0)*(v_dir>0);
                                      int j_ns_amont = j_ns + dj - (dir == 1)*(v_dir>0);
                                      int k_ns_amont = k_ns + dk - (dir == 2)*(v_dir>0);

                                      double old_indicatrice_amont  = cut_cell_disc.get_interfaces().I(i_ns_amont,j_ns_amont,k_ns_amont);
                                      double next_indicatrice_amont = cut_cell_disc.get_interfaces().In(i_ns_amont,j_ns_amont,k_ns_amont);

                                      double amont_bary_x = origin_x + dx*(i_ft_amont + offset_x + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 0, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_y = origin_y + dy*(j_ft_amont + offset_y + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 1, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_z = origin_z + dz*(k_ft_amont + offset_z + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 2, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));

                                      double face_bary_x = origin_x + dx*(i_ft + di + offset_x + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 0, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_y = origin_y + dy*(j_ft + dj + offset_y + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 1, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_z = origin_z + dz*(k_ft + dk + offset_z + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 2, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));

                                      double normal_dir = select(dir, normal_x, normal_y, normal_z);
                                      bool case_perpendicular = (std::abs(normal_dir) > 0.72);

                                      coordinates[0] = case_perpendicular ? face_bary_x : amont_bary_x;
                                      coordinates[1] = case_perpendicular ? face_bary_y : amont_bary_y;
                                      coordinates[2] = case_perpendicular ? face_bary_z : amont_bary_z;
                                    }
                                  else if (face_interp == CUT_CELL_CONV_FACE_INTERPOLATION::INTERPOLATEDISTANCE0)
                                    {
                                      double v_x = (n < 0) ? cut_field_total_velocity[0].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[0].diph_v_(n) : cut_field_total_velocity[0].diph_l_(n));
                                      double v_y = (n < 0) ? cut_field_total_velocity[1].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[1].diph_v_(n) : cut_field_total_velocity[1].diph_l_(n));
                                      double v_z = (n < 0) ? cut_field_total_velocity[2].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[2].diph_v_(n) : cut_field_total_velocity[2].diph_l_(n));
                                      double v_dir = select(dir, v_x, v_y, v_z);

                                      int i_ft_amont = i_ft + di - (dir == 0)*(v_dir>0);
                                      int j_ft_amont = j_ft + dj - (dir == 1)*(v_dir>0);
                                      int k_ft_amont = k_ft + dk - (dir == 2)*(v_dir>0);

                                      int i_ns_amont = i_ns + di - (dir == 0)*(v_dir>0);
                                      int j_ns_amont = j_ns + dj - (dir == 1)*(v_dir>0);
                                      int k_ns_amont = k_ns + dk - (dir == 2)*(v_dir>0);

                                      double old_indicatrice_amont  = cut_cell_disc.get_interfaces().I(i_ns_amont,j_ns_amont,k_ns_amont);
                                      double next_indicatrice_amont = cut_cell_disc.get_interfaces().In(i_ns_amont,j_ns_amont,k_ns_amont);

                                      double amont_bary_x = origin_x + dx*(i_ft_amont + offset_x + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 0, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_y = origin_y + dy*(j_ft_amont + offset_y + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 1, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_z = origin_z + dz*(k_ft_amont + offset_z + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 2, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));

                                      double face_bary_x = origin_x + dx*(i_ft + di + offset_x + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 0, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_y = origin_y + dy*(j_ft + dj + offset_y + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 1, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_z = origin_z + dz*(k_ft + dk + offset_z + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 2, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));

                                      double dist_x = (amont_bary_x - face_bary_x);
                                      double dist_y = (amont_bary_y - face_bary_y);
                                      double dist_z = (amont_bary_z - face_bary_z);
                                      double dist_face_amont = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
                                      bool close = (dist_face_amont < 0.1);

                                      coordinates[0] = close ? amont_bary_x : face_bary_x;
                                      coordinates[1] = close ? amont_bary_y : face_bary_y;
                                      coordinates[2] = close ? amont_bary_z : face_bary_z;
                                    }
                                  else if (face_interp == CUT_CELL_CONV_FACE_INTERPOLATION::INTERPOLATEDISTANCE1)
                                    {
                                      double v_x = (n < 0) ? cut_field_total_velocity[0].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[0].diph_v_(n) : cut_field_total_velocity[0].diph_l_(n));
                                      double v_y = (n < 0) ? cut_field_total_velocity[1].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[1].diph_v_(n) : cut_field_total_velocity[1].diph_l_(n));
                                      double v_z = (n < 0) ? cut_field_total_velocity[2].pure_(i_ns+di,j_ns+dj,k_ns+dk) : ((phase == 0) ? cut_field_total_velocity[2].diph_v_(n) : cut_field_total_velocity[2].diph_l_(n));
                                      double v_dir = select(dir, v_x, v_y, v_z);

                                      int i_ft_amont = i_ft + di - (dir == 0)*(v_dir>0);
                                      int j_ft_amont = j_ft + dj - (dir == 1)*(v_dir>0);
                                      int k_ft_amont = k_ft + dk - (dir == 2)*(v_dir>0);

                                      int i_ns_amont = i_ns + di - (dir == 0)*(v_dir>0);
                                      int j_ns_amont = j_ns + dj - (dir == 1)*(v_dir>0);
                                      int k_ns_amont = k_ns + dk - (dir == 2)*(v_dir>0);

                                      double old_indicatrice_amont  = cut_cell_disc.get_interfaces().I(i_ns_amont,j_ns_amont,k_ns_amont);
                                      double next_indicatrice_amont = cut_cell_disc.get_interfaces().In(i_ns_amont,j_ns_amont,k_ns_amont);

                                      double amont_bary_x = origin_x + dx*(i_ft_amont + offset_x + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 0, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_y = origin_y + dy*(j_ft_amont + offset_y + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 1, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));
                                      double amont_bary_z = origin_z + dz*(k_ft_amont + offset_z + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 2, phase, i_ns_amont,j_ns_amont,k_ns_amont, old_indicatrice_amont, next_indicatrice_amont)));

                                      double face_bary_x = origin_x + dx*(i_ft + di + offset_x + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 0, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_y = origin_y + dy*(j_ft + dj + offset_y + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 1, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));
                                      double face_bary_z = origin_z + dz*(k_ft + dk + offset_z + (cut_cell_disc.get_interfaces().get_barycentre_face(next_time, dir, 2, phase, i_ns+di,j_ns+dj,k_ns+dk, old_indicatrice_surfacique, next_indicatrice_surfacique)));

                                      double dist_x = (amont_bary_x - face_bary_x);
                                      double dist_y = (amont_bary_y - face_bary_y);
                                      double dist_z = (amont_bary_z - face_bary_z);
                                      double dist_face_amont = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
                                      bool close = (dist_face_amont < 0.1);

                                      coordinates[0] = close ? face_bary_x : amont_bary_x;
                                      coordinates[1] = close ? face_bary_y : amont_bary_y;
                                      coordinates[2] = close ? face_bary_z : amont_bary_z;
                                    }
                                  else
                                    {
                                      Cerr << "Cut_cell_convection_auxiliaire::calcule_temperature_face_depuis_facette_interpolate: face_interp inconnue." << finl;
                                      Process::exit();
                                    }

                                  int status = -2;
                                  const int tolerate_not_within_tetrahedron = std::max(2, tolerate_not_within_tetrahedron_);
                                  double temperature_interpolate = ijk_interpolate_cut_cell_using_interface(false, phase, temperature_ft, cut_field_temperature, interfacial_temperature, coordinates, tolerate_not_within_tetrahedron, status);

                                  if (status == -1)
                                    {
                                      temperature_interpolate = 0.;
                                    }

                                  // Note : Quand le statut est -1, la valeur de la temperature est souvent absolument pas pertinente
                                  temperature_face_ft[phase][dir](i_ft+di,j_ft+dj,k_ft+dk) = (status == -1 || indicatrice_surface == 0.) ? 0. : temperature_interpolate;

                                  assert(status > -2);
                                  if (status == -1)
                                    {
                                      count_status_not_found++;
                                    }
                                  else if (status < 10)
                                    {
                                      count_status_below_10++;
                                    }
                                  else if (status < 100)
                                    {
                                      count_status_below_100++;
                                    }
                                  else if (status < 1000)
                                    {
                                      count_status_below_1000++;
                                    }
                                  else if (status < 10000)
                                    {
                                      count_status_below_10000++;
                                    }
                                  else if (status < 100000)
                                    {
                                      count_status_below_100000++;
                                    }
                                  else if (status < 1000000)
                                    {
                                      count_status_below_1000000++;
                                    }
                                  else if (status >= 1000000)
                                    {
                                      count_status_above_1000000++;
                                    }
                                  else
                                    {
                                      Cerr << "Un tel statut n'est pas possible." << finl;
                                      Process::exit();
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
  }

  for (int dir = 0; dir < 3; dir++)
    {
      for (int phase = 0; phase < 2; phase++)
        {
          temperature_face_ft[phase][dir].echange_espace_virtuel(temperature_face_ft[phase][dir].ghost());

          // Calcul du flux interface sur le domaine NS :
          ref_ijk_ft->redistrib_from_ft_elem().redistribute(temperature_face_ft[phase][dir], temperature_face_ns[phase][dir]);
          temperature_face_ns[phase][dir].echange_espace_virtuel(temperature_face_ns[phase][dir].ghost());
        }
    }

  Cerr << "Bilan des statuts pour Cut_cell_convection_auxiliaire::calcule_temperature_face_depuis_facette_interpolate : " << finl;
  Cerr << "    -1      " << count_status_not_found  << finl;
  Cerr << "   <10      " << count_status_below_10   << finl;
  Cerr << "  <100      " << count_status_below_100  << finl;
  Cerr << " <1000      " << count_status_below_1000 << finl;
  Cerr << " <10000     " << count_status_below_10000 << finl;
  Cerr << " <100000    " << count_status_below_100000 << finl;
  Cerr << " <1000000   " << count_status_below_1000000 << finl;
  Cerr << ">=1000000   " << count_status_above_1000000 << finl;
}
