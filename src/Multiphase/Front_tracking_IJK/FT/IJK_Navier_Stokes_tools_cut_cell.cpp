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

#include <IJK_Navier_Stokes_tools_cut_cell.h>
#include <Champ_diphasique.h>
#include <Cut_cell_FT_Disc.h>

struct struct_int_double
{
  int index;
  double value;
};

int compare_value(const void *a, const void *b)
{
  struct_int_double *a1 = (struct_int_double *)a;
  struct_int_double *a2 = (struct_int_double *)b;
  if ((*a1).value < (*a2).value)
    return -1;
  else if ((*a1).value > (*a2).value)
    return 1;
  else
    return 0;
}

extern "C" {
  // LU decomoposition of a general matrix
  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

  // generate inverse of a matrix given its LU decomposition
  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}

void inverse(double* A, int N)
{
  int *IPIV = new int[N];
  int LWORK = N*N;
  double *WORK = new double[LWORK];
  int INFO;

  dgetrf_(&N,&N,A,&N,IPIV,&INFO);
  dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

  delete[] IPIV;
  delete[] WORK;
}

static double ijk_interpolate_cut_cell_for_given_index(int idx, int phase, Cut_field_scalar& field, const DoubleTab& coordinates, ArrOfDouble& result, int skip_unknown_points, double value_for_bad_points)
{
  const Cut_cell_FT_Disc& cut_cell_disc = field.get_cut_cell_disc();

  //const int ghost = field.pure_.ghost();
  const int ghost = cut_cell_disc.get_ghost_size();
  assert(field.pure_.ghost() >= ghost);

  const int ni = field.pure_.ni();
  const int nj = field.pure_.nj();
  const int nk = field.pure_.nk();

  const IJK_Splitting& splitting = field.pure_.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  //const IJK_Splitting::Localisation loc = field.pure_.get_localisation();
  // L'origine est sur un noeud. Donc que la premiere face en I est sur get_origin(DIRECTION_I)
  double origin_x = geom.get_origin(DIRECTION_I); //+ ((loc == IJK_Splitting::FACES_J || loc == IJK_Splitting::FACES_K || loc == IJK_Splitting::ELEM) ? (dx * 0.5) : 0.);
  double origin_y = geom.get_origin(DIRECTION_J); //+ ((loc == IJK_Splitting::FACES_K || loc == IJK_Splitting::FACES_I || loc == IJK_Splitting::ELEM) ? (dy * 0.5) : 0.);
  double origin_z = geom.get_origin(DIRECTION_K); //+ ((loc == IJK_Splitting::FACES_I || loc == IJK_Splitting::FACES_J || loc == IJK_Splitting::ELEM) ? (dz * 0.5) : 0.);

  const int offset_x = splitting.get_offset_local(DIRECTION_I);
  const int offset_y = splitting.get_offset_local(DIRECTION_J);
  const int offset_z = splitting.get_offset_local(DIRECTION_K);

  const double x = coordinates(idx, 0);
  const double y = coordinates(idx, 1);
  const double z = coordinates(idx, 2);

  const double x2 = (x - origin_x) / dx;
  const double y2 = (y - origin_y) / dy;
  const double z2 = (z - origin_z) / dz;

  // Coordonnes barycentriques du points dans la cellule :
  const double xfact = x2 - floor(x2);
  const double yfact = y2 - floor(y2);
  const double zfact = z2 - floor(z2);

  // On travaille sur le maillage NS, on va donc corrige les indices de la periodicite.
  // Note : on ne corrige que l'index et pas les coordonnees, car on n'utilise plus les coordonnees par la suite.
  const int index_i = cut_cell_disc.get_i_selon_dir(0, x, 2, splitting, false);
  const int index_j = cut_cell_disc.get_i_selon_dir(1, y, 2, splitting, false);
  const int index_k = cut_cell_disc.get_i_selon_dir(2, z, 2, splitting, false);

  // is point in the domain ? (ghost cells ok...)
  bool ok = (index_i >= -ghost && index_i < ni + ghost) && (index_j >= -ghost && index_j < nj + ghost) && (index_k >= -ghost && index_k < nk + ghost);
  bool close_to_edge = (index_i == -ghost || index_i == ni + ghost - 1) || (index_j == -ghost || index_j == nj + ghost - 1) || (index_k == -ghost || index_k == nk + ghost - 1);
  if (!ok)
    {
      if (skip_unknown_points)
        {
          return value_for_bad_points;
        }
      else
        {
          // Error!
          Cerr << "Error in ijk_interpolate_cut_cell_implementation: request cut-cell interpolation of point " << x << " " << y << " " << z << " which is outside of the domain on processor " << Process::me() << finl;
          Process::exit();
        }
    }

  const int max_number_of_candidates = 27;
  Int3 candidate_offset[max_number_of_candidates] =
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

  int number_of_candidates = 0;
  struct_int_double dist[max_number_of_candidates];
  for (int i = 0; i < max_number_of_candidates; i++)
    {
      dist[i].index = i;

      int i_candidate_aperio = index_i + candidate_offset[i][0];
      int j_candidate_aperio = index_j + candidate_offset[i][1];
      int k_candidate_aperio = index_k + candidate_offset[i][2];

      // Prise en compte de la periodicite
      double x_candidate_centred_aperio = (i_candidate_aperio + offset_x + .5)*dx + origin_x;
      double y_candidate_centred_aperio = (j_candidate_aperio + offset_y + .5)*dy + origin_y;
      double z_candidate_centred_aperio = (k_candidate_aperio + offset_z + .5)*dz + origin_z;
      int i_candidate = cut_cell_disc.get_i_selon_dir(0, x_candidate_centred_aperio, 2, splitting, false);
      int j_candidate = cut_cell_disc.get_i_selon_dir(1, y_candidate_centred_aperio, 2, splitting, false);
      int k_candidate = cut_cell_disc.get_i_selon_dir(2, z_candidate_centred_aperio, 2, splitting, false);
      assert((i_candidate_aperio == i_candidate) || (close_to_edge));
      assert((j_candidate_aperio == j_candidate) || (close_to_edge));
      assert((k_candidate_aperio == k_candidate) || (close_to_edge));

      double next_indicatrice = cut_cell_disc.get_interfaces().In(i_candidate, j_candidate, k_candidate);
      if ((phase == 0 && next_indicatrice == 1.) || (phase == 1 && next_indicatrice == 0.))
        {
          dist[i].value = DMAXFLOAT; // Point invalide
        }
      else
        {
          number_of_candidates += 1;
          double old_indicatrice = cut_cell_disc.get_interfaces().I(i_candidate, j_candidate, k_candidate);
          double candidate_x = (double)candidate_offset[i][0] + cut_cell_disc.get_interfaces().get_barycentre_next(0, phase, i_candidate, j_candidate, k_candidate, old_indicatrice, next_indicatrice);
          double candidate_y = (double)candidate_offset[i][1] + cut_cell_disc.get_interfaces().get_barycentre_next(1, phase, i_candidate, j_candidate, k_candidate, old_indicatrice, next_indicatrice);
          double candidate_z = (double)candidate_offset[i][2] + cut_cell_disc.get_interfaces().get_barycentre_next(2, phase, i_candidate, j_candidate, k_candidate, old_indicatrice, next_indicatrice);
          double decalage_x = candidate_x - xfact;
          double decalage_y = candidate_y - yfact;
          double decalage_z = candidate_z - zfact;
          dist[i].value = sqrt(decalage_x*decalage_x + decalage_y*decalage_y + decalage_z*decalage_z);
        }
    }
  assert(number_of_candidates <= max_number_of_candidates);
  qsort(dist, max_number_of_candidates, sizeof(struct_int_double), compare_value);

  // On boucle d'abord sur n_neighbours, le nombre de points que l'on s'autorise a chercher
  // dans la liste des voisins. Par exemple, n_neighbours=6 veut dire que l'on cherche a former
  // des tetraedres a partir des 6 premiers voisins.
  for (int n_neighbours = 4; n_neighbours < number_of_candidates; n_neighbours++)
    {
      int index_vertex0 = n_neighbours; // Le premier point est toujours fixe au dernier possible. Cela garanti que tous les tetraedres d'un nouveau n_neighbours sont nouveaux.
      int vertex0 = dist[index_vertex0].index;

      int i_0_aperio = index_i + candidate_offset[vertex0][0];
      int j_0_aperio = index_j + candidate_offset[vertex0][1];
      int k_0_aperio = index_k + candidate_offset[vertex0][2];

      // Prise en compte de la periodicite
      double x_0_centred_aperio = (i_0_aperio + offset_x + .5)*dx + origin_x;
      double y_0_centred_aperio = (j_0_aperio + offset_y + .5)*dy + origin_y;
      double z_0_centred_aperio = (k_0_aperio + offset_z + .5)*dz + origin_z;
      int i_0 = cut_cell_disc.get_i_selon_dir(0, x_0_centred_aperio, 2, splitting, false);
      int j_0 = cut_cell_disc.get_i_selon_dir(1, y_0_centred_aperio, 2, splitting, false);
      int k_0 = cut_cell_disc.get_i_selon_dir(2, z_0_centred_aperio, 2, splitting, false);
      assert((i_0_aperio == i_0) || (close_to_edge));
      assert((j_0_aperio == j_0) || (close_to_edge));
      assert((k_0_aperio == k_0) || (close_to_edge));

      double next_indicatrice_0 = cut_cell_disc.get_interfaces().In(i_0, j_0, k_0);
      double old_indicatrice_0 = cut_cell_disc.get_interfaces().I(i_0, j_0, k_0);
      double x_0 = (double)candidate_offset[vertex0][0] + cut_cell_disc.get_interfaces().get_barycentre_next(0, phase, i_0, j_0, k_0, old_indicatrice_0, next_indicatrice_0);
      double y_0 = (double)candidate_offset[vertex0][1] + cut_cell_disc.get_interfaces().get_barycentre_next(1, phase, i_0, j_0, k_0, old_indicatrice_0, next_indicatrice_0);
      double z_0 = (double)candidate_offset[vertex0][2] + cut_cell_disc.get_interfaces().get_barycentre_next(2, phase, i_0, j_0, k_0, old_indicatrice_0, next_indicatrice_0);
      //double dx_0 = 0.;
      //double dy_0 = 0.;
      //double dz_0 = 0.;

      for (int index_vertex1 = 0; index_vertex1 < index_vertex0-2; index_vertex1++)
        {
          int vertex1 = dist[index_vertex1].index;

          int i_1_aperio = index_i + candidate_offset[vertex1][0];
          int j_1_aperio = index_j + candidate_offset[vertex1][1];
          int k_1_aperio = index_k + candidate_offset[vertex1][2];

          // Prise en compte de la periodicite
          double x_1_centred_aperio = (i_1_aperio + offset_x + .5)*dx + origin_x;
          double y_1_centred_aperio = (j_1_aperio + offset_y + .5)*dy + origin_y;
          double z_1_centred_aperio = (k_1_aperio + offset_z + .5)*dz + origin_z;
          int i_1 = cut_cell_disc.get_i_selon_dir(0, x_1_centred_aperio, 2, splitting, false);
          int j_1 = cut_cell_disc.get_i_selon_dir(1, y_1_centred_aperio, 2, splitting, false);
          int k_1 = cut_cell_disc.get_i_selon_dir(2, z_1_centred_aperio, 2, splitting, false);
          assert((i_1_aperio == i_1) || (close_to_edge));
          assert((j_1_aperio == j_1) || (close_to_edge));
          assert((k_1_aperio == k_1) || (close_to_edge));

          double next_indicatrice_1 = cut_cell_disc.get_interfaces().In(i_1, j_1, k_1);
          double old_indicatrice_1 = cut_cell_disc.get_interfaces().I(i_1, j_1, k_1);
          double x_1 = (double)candidate_offset[vertex1][0] + cut_cell_disc.get_interfaces().get_barycentre_next(0, phase, i_1, j_1, k_1, old_indicatrice_1, next_indicatrice_1);
          double y_1 = (double)candidate_offset[vertex1][1] + cut_cell_disc.get_interfaces().get_barycentre_next(1, phase, i_1, j_1, k_1, old_indicatrice_1, next_indicatrice_1);
          double z_1 = (double)candidate_offset[vertex1][2] + cut_cell_disc.get_interfaces().get_barycentre_next(2, phase, i_1, j_1, k_1, old_indicatrice_1, next_indicatrice_1);
          double dx_1 = x_1 - x_0;
          double dy_1 = y_1 - y_0;
          double dz_1 = z_1 - z_0;

          for (int index_vertex2 = index_vertex1+1; index_vertex2 < index_vertex0-1; index_vertex2++)
            {
              int vertex2 = dist[index_vertex2].index;

              int i_2_aperio = index_i + candidate_offset[vertex2][0];
              int j_2_aperio = index_j + candidate_offset[vertex2][1];
              int k_2_aperio = index_k + candidate_offset[vertex2][2];

              // Prise en compte de la periodicite
              double x_2_centred_aperio = (i_2_aperio + offset_x + .5)*dx + origin_x;
              double y_2_centred_aperio = (j_2_aperio + offset_y + .5)*dy + origin_y;
              double z_2_centred_aperio = (k_2_aperio + offset_z + .5)*dz + origin_z;
              int i_2 = cut_cell_disc.get_i_selon_dir(0, x_2_centred_aperio, 2, splitting, false);
              int j_2 = cut_cell_disc.get_i_selon_dir(1, y_2_centred_aperio, 2, splitting, false);
              int k_2 = cut_cell_disc.get_i_selon_dir(2, z_2_centred_aperio, 2, splitting, false);
              assert((i_2_aperio == i_2) || (close_to_edge));
              assert((j_2_aperio == j_2) || (close_to_edge));
              assert((k_2_aperio == k_2) || (close_to_edge));

              double next_indicatrice_2 = cut_cell_disc.get_interfaces().In(i_2, j_2, k_2);
              double old_indicatrice_2 = cut_cell_disc.get_interfaces().I(i_2, j_2, k_2);
              double x_2 = (double)candidate_offset[vertex2][0] + cut_cell_disc.get_interfaces().get_barycentre_next(0, phase, i_2, j_2, k_2, old_indicatrice_2, next_indicatrice_2);
              double y_2 = (double)candidate_offset[vertex2][1] + cut_cell_disc.get_interfaces().get_barycentre_next(1, phase, i_2, j_2, k_2, old_indicatrice_2, next_indicatrice_2);
              double z_2 = (double)candidate_offset[vertex2][2] + cut_cell_disc.get_interfaces().get_barycentre_next(2, phase, i_2, j_2, k_2, old_indicatrice_2, next_indicatrice_2);

              double dx_2 = x_2 - x_0;
              double dy_2 = y_2 - y_0;
              double dz_2 = z_2 - z_0;

              for (int index_vertex3 = index_vertex2+1; index_vertex3 < index_vertex0; index_vertex3++)
                {
                  int vertex3 = dist[index_vertex3].index;

                  int i_3_aperio = index_i + candidate_offset[vertex3][0];
                  int j_3_aperio = index_j + candidate_offset[vertex3][1];
                  int k_3_aperio = index_k + candidate_offset[vertex3][2];

                  // Prise en compte de la periodicite
                  double x_3_centred_aperio = (i_3_aperio + offset_x + .5)*dx + origin_x;
                  double y_3_centred_aperio = (j_3_aperio + offset_y + .5)*dy + origin_y;
                  double z_3_centred_aperio = (k_3_aperio + offset_z + .5)*dz + origin_z;
                  int i_3 = cut_cell_disc.get_i_selon_dir(0, x_3_centred_aperio, 2, splitting, false);
                  int j_3 = cut_cell_disc.get_i_selon_dir(1, y_3_centred_aperio, 2, splitting, false);
                  int k_3 = cut_cell_disc.get_i_selon_dir(2, z_3_centred_aperio, 2, splitting, false);
                  assert((i_3_aperio == i_3) || (close_to_edge));
                  assert((j_3_aperio == j_3) || (close_to_edge));
                  assert((k_3_aperio == k_3) || (close_to_edge));

                  double next_indicatrice_3 = cut_cell_disc.get_interfaces().In(i_3, j_3, k_3);
                  double old_indicatrice_3 = cut_cell_disc.get_interfaces().I(i_3, j_3, k_3);
                  double x_3 = (double)candidate_offset[vertex3][0] + cut_cell_disc.get_interfaces().get_barycentre_next(0, phase, i_3, j_3, k_3, old_indicatrice_3, next_indicatrice_3);
                  double y_3 = (double)candidate_offset[vertex3][1] + cut_cell_disc.get_interfaces().get_barycentre_next(1, phase, i_3, j_3, k_3, old_indicatrice_3, next_indicatrice_3);
                  double z_3 = (double)candidate_offset[vertex3][2] + cut_cell_disc.get_interfaces().get_barycentre_next(2, phase, i_3, j_3, k_3, old_indicatrice_3, next_indicatrice_3);

                  double dx_3 = x_3 - x_0;
                  double dy_3 = y_3 - y_0;
                  double dz_3 = z_3 - z_0;

                  // En coordonnees barycentriques, puisque dX_0 = 0
                  //
                  // x_target_to_interpolate = dx_1 lambda_1 + dx_2 lamdba_2 + dx_3 lambda_3
                  // y_target_to_interpolate = dy_1 lambda_1 + dy_2 lamdba_2 + dy_3 lambda_3
                  // z_target_to_interpolate = dz_1 lambda_1 + dz_2 lamdba_2 + dz_3 lambda_3
                  // ---
                  // X_target_to_interpolate = Matrix * Lambda_vector
                  //
                  // Donc Lambda_vector = Matrix^-1 * X_target_to_interpolate
                  double Matrix[3*3] =
                  {
                    dx_1, dx_2, dx_3,
                    dy_1, dy_2, dy_3,
                    dz_1, dz_2, dz_3,
                  };

                  inverse(Matrix, 3);

                  double lambda_1 = Matrix[0] * (xfact - x_0) + Matrix[1] * (yfact - y_0) + Matrix[2] * (zfact - z_0);
                  double lambda_2 = Matrix[3] * (xfact - x_0) + Matrix[4] * (yfact - y_0) + Matrix[5] * (zfact - z_0);
                  double lambda_3 = Matrix[6] * (xfact - x_0) + Matrix[7] * (yfact - y_0) + Matrix[8] * (zfact - z_0);

                  double lambda_0 = 1 - lambda_1 - lambda_2 - lambda_3;

                  int lambda_0_within_bounds = ((lambda_0 >= 0) && (lambda_0 <= 1));
                  int lambda_1_within_bounds = ((lambda_1 >= 0) && (lambda_1 <= 1));
                  int lambda_2_within_bounds = ((lambda_2 >= 0) && (lambda_2 <= 1));
                  int lambda_3_within_bounds = ((lambda_3 >= 0) && (lambda_3 <= 1));

                  int point_within_tetrahedron = lambda_0_within_bounds && lambda_1_within_bounds && lambda_2_within_bounds && lambda_3_within_bounds;

                  if (point_within_tetrahedron)
                    {
                      double field_0;
                      double field_1;
                      double field_2;
                      double field_3;

                      int n_0 = cut_cell_disc.get_n(i_0, j_0, k_0);
                      if (n_0 >= 0)
                        {
                          field_0 = (phase == 0) ? field.diph_v_(n_0) : field.diph_l_(n_0);
                        }
                      else
                        {
                          assert(cut_cell_disc.get_interfaces().In(i_0,j_0,k_0) == (double)phase);
                          field_0 = field.pure_(i_0,j_0,k_0);
                        }

                      int n_1 = cut_cell_disc.get_n(i_1, j_1, k_1);
                      if (n_1 >= 0)
                        {
                          field_1 = (phase == 0) ? field.diph_v_(n_1) : field.diph_l_(n_1);
                        }
                      else
                        {
                          assert(cut_cell_disc.get_interfaces().In(i_1,j_1,k_1) == (double)phase);
                          field_1 = field.pure_(i_1,j_1,k_1);
                        }

                      int n_2 = cut_cell_disc.get_n(i_2, j_2, k_2);
                      if (n_2 >= 0)
                        {
                          field_2 = (phase == 0) ? field.diph_v_(n_2) : field.diph_l_(n_2);
                        }
                      else
                        {
                          assert(cut_cell_disc.get_interfaces().In(i_2,j_2,k_2) == (double)phase);
                          field_2 = field.pure_(i_2,j_2,k_2);
                        }

                      int n_3 = cut_cell_disc.get_n(i_3, j_3, k_3);
                      if (n_3 >= 0)
                        {
                          field_3 = (phase == 0) ? field.diph_v_(n_3) : field.diph_l_(n_3);
                        }
                      else
                        {
                          assert(cut_cell_disc.get_interfaces().In(i_3,j_3,k_3) == (double)phase);
                          field_3 = field.pure_(i_3,j_3,k_3);
                        }

                      double r = field_0*lambda_0 + field_1*lambda_1 + field_2*lambda_2 + field_3*lambda_3;

                      return r;
                    }
                }
            }
        }
    }

  // No suitable tetrahedron was found
  Cerr << "Value of close_to_edge: " << (int)close_to_edge << finl;
  Cerr << "Error in ijk_interpolate_cut_cell_for_given_index: no tetrahedron containing the point " << x << " " << y << " " << z << " on processor " << Process::me() << finl;
  Process::exit();
  return -1;
}

// Interpolate the "field" at the requested "coordinates" (array with 3 columns), and stores into "result"
static void ijk_interpolate_cut_cell_implementation(int phase, Cut_field_scalar& field, const DoubleTab& coordinates, ArrOfDouble& result, int skip_unknown_points, double value_for_bad_points)
{
  const int nb_coords = coordinates.dimension(0);
  result.resize_array(nb_coords);
  for (int idx = 0; idx < nb_coords; idx++)
    {
      double interpolated_value = ijk_interpolate_cut_cell_for_given_index(idx, phase, field, coordinates, result, skip_unknown_points, value_for_bad_points);
      result[idx] = interpolated_value;
    }
}
void ijk_interpolate_cut_cell_skip_unknown_points(int phase, Cut_field_scalar& field, const DoubleTab& coordinates, ArrOfDouble& result, const double value_for_bad_points)
{
  ijk_interpolate_cut_cell_implementation(phase, field, coordinates, result, 1 /* yes:skip unknown points */, value_for_bad_points);
}

void ijk_interpolate_cut_cell(int phase, Cut_field_scalar& field, const DoubleTab& coordinates, ArrOfDouble& result)
{
  ijk_interpolate_cut_cell_implementation(phase, field, coordinates, result, 0 /* skip unknown points=no */, 0.);
}
