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
#include <Cut_field.h>
#include <IJK_Interfaces.h>
#include <Cut_cell_FT_Disc.h>

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

struct IndexCoefficient
{
  int index;
  double coefficient;
};

struct Sommet
{
  int sommet;
  int fa7;
  double value;
  int count;
};

struct Candidate
{
  int index;
  double dist;
  Vecteur3 coord;
  int signed_independent_index;
};

struct Candidate_with_direct_value
{
  int index;
  double dist;
  Vecteur3 coord;
  double value;
};

static const int max_number_of_involved_sommet = 512; // Note: Pour ce maximum, les sommets sont comptes une fois pour chaque facette et pour chaque cellule contenant cette facette

static const int max_number_of_cell_candidates = 27;
static const Int3 candidate_offset[max_number_of_cell_candidates] =
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

static const int max_number_of_candidates = max_number_of_involved_sommet + max_number_of_cell_candidates;

template<typename T>
Vecteur3 compute_lambda(int index_vertex0, int index_vertex1, int index_vertex2, int index_vertex3, T candidates[max_number_of_candidates], double xfact, double yfact, double zfact)
{
  double x_0 = candidates[index_vertex0].coord[0];
  double y_0 = candidates[index_vertex0].coord[1];
  double z_0 = candidates[index_vertex0].coord[2];

  double x_1 = candidates[index_vertex1].coord[0];
  double y_1 = candidates[index_vertex1].coord[1];
  double z_1 = candidates[index_vertex1].coord[2];
  double dx_1 = x_1 - x_0;
  double dy_1 = y_1 - y_0;
  double dz_1 = z_1 - z_0;

  double x_2 = candidates[index_vertex2].coord[0];
  double y_2 = candidates[index_vertex2].coord[1];
  double z_2 = candidates[index_vertex2].coord[2];

  double dx_2 = x_2 - x_0;
  double dy_2 = y_2 - y_0;
  double dz_2 = z_2 - z_0;

  double x_3 = candidates[index_vertex3].coord[0];
  double y_3 = candidates[index_vertex3].coord[1];
  double z_3 = candidates[index_vertex3].coord[2];

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

  Vecteur3 lambda_vec = {lambda_1, lambda_2, lambda_3};
  return lambda_vec;
}

static void ijk_interpolate_cut_cell_for_given_index(bool next_time, int phase, const Cut_cell_FT_Disc& cut_cell_disc, const Vecteur3& coordinates, IndexCoefficient tetrahedra_index_coefficient[4], int tolerate_not_within_tetrahedron, int skip_unknown_points, double value_for_bad_points, int& status)
{
  const int ghost = cut_cell_disc.get_ghost_size();
  const int reduced_ghost = ghost - 1;

  const double x = coordinates[0];
  const double y = coordinates[1];
  const double z = coordinates[2];

  const IJK_Splitting& splitting = cut_cell_disc.get_splitting();
  const int ni = splitting.get_nb_items_local(IJK_Splitting::ELEM, 0);
  const int nj = splitting.get_nb_items_local(IJK_Splitting::ELEM, 1);
  const int nk = splitting.get_nb_items_local(IJK_Splitting::ELEM, 2);

  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  // L'origine est sur un noeud. Donc que la premiere face en I est sur get_origin(DIRECTION_I)
  double origin_x = geom.get_origin(DIRECTION_I);
  double origin_y = geom.get_origin(DIRECTION_J);
  double origin_z = geom.get_origin(DIRECTION_K);

  const double x2 = (x - origin_x) / dx;
  const double y2 = (y - origin_y) / dy;
  const double z2 = (z - origin_z) / dz;

  // Coordonnes barycentriques du points dans la cellule :
  const double xfact = x2 - floor(x2);
  const double yfact = y2 - floor(y2);
  const double zfact = z2 - floor(z2);

  // On travaille sur le maillage NS, on va donc corrige les indices de la periodicite.
  // Note : on ne corrige que l'index et pas les coordonnees, car on n'utilise plus les coordonnees par la suite.
  const int index_i = cut_cell_disc.get_splitting().get_i_along_dir_perio(DIRECTION_I, x, IJK_Splitting::ELEM);
  const int index_j = cut_cell_disc.get_splitting().get_i_along_dir_perio(DIRECTION_J, y, IJK_Splitting::ELEM);
  const int index_k = cut_cell_disc.get_splitting().get_i_along_dir_perio(DIRECTION_K, z, IJK_Splitting::ELEM);

  // is point in the domain ? (ghost cells ok...)
  bool ok = (index_i >= -reduced_ghost && index_i < ni + reduced_ghost) && (index_j >= -reduced_ghost && index_j < nj + reduced_ghost) && (index_k >= -reduced_ghost && index_k < nk + reduced_ghost);
  bool close_to_edge = (index_i == -reduced_ghost || index_i == ni + reduced_ghost - 1) || (index_j == -reduced_ghost || index_j == nj + reduced_ghost - 1) || (index_k == -reduced_ghost || index_k == nk + reduced_ghost - 1);
  if (!ok)
    {
      if (skip_unknown_points)
        {
          tetrahedra_index_coefficient[0] = {-1,value_for_bad_points};
          tetrahedra_index_coefficient[1] = {-1,value_for_bad_points};
          tetrahedra_index_coefficient[2] = {-1,value_for_bad_points};
          tetrahedra_index_coefficient[3] = {-1,value_for_bad_points};
          return;
        }
      else
        {
          // Error!
          Cerr << "Error in ijk_interpolate_cut_cell_implementation: request cut-cell interpolation of point " << x << " " << y << " " << z << " which is outside of the domain on processor " << Process::me() << finl;
          Process::exit();
        }
    }

  int number_of_candidates = 0;
  Candidate candidates[max_number_of_candidates];
  for (int i = 0; i < max_number_of_cell_candidates; i++)
    {
      int i_candidate_aperio = index_i + candidate_offset[i][0];
      int j_candidate_aperio = index_j + candidate_offset[i][1];
      int k_candidate_aperio = index_k + candidate_offset[i][2];

      // Prise en compte de la periodicite
      int i_candidate = cut_cell_disc.get_splitting().correct_perio_i_local(0, i_candidate_aperio);
      int j_candidate = cut_cell_disc.get_splitting().correct_perio_i_local(1, j_candidate_aperio);
      int k_candidate = cut_cell_disc.get_splitting().correct_perio_i_local(2, k_candidate_aperio);

      double old_indicatrice = cut_cell_disc.get_interfaces().I(i_candidate, j_candidate, k_candidate);
      double next_indicatrice = cut_cell_disc.get_interfaces().In(i_candidate, j_candidate, k_candidate);
      double indicatrice = next_time ? next_indicatrice : old_indicatrice;
      if ((phase == 0 && indicatrice == 1.) || (phase == 1 && indicatrice == 0.))
        {
          // Point invalide
        }
      else
        {
          double candidate_x = (double)candidate_offset[i][0] + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 0, phase, i_candidate, j_candidate, k_candidate));
          double candidate_y = (double)candidate_offset[i][1] + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 1, phase, i_candidate, j_candidate, k_candidate));
          double candidate_z = (double)candidate_offset[i][2] + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 2, phase, i_candidate, j_candidate, k_candidate));
          double decalage_x = candidate_x - xfact;
          double decalage_y = candidate_y - yfact;
          double decalage_z = candidate_z - zfact;
          candidates[number_of_candidates].index = number_of_candidates;
          candidates[number_of_candidates].dist = sqrt(decalage_x*decalage_x + decalage_y*decalage_y + decalage_z*decalage_z);
          candidates[number_of_candidates].coord[0] = candidate_x;
          candidates[number_of_candidates].coord[1] = candidate_y;
          candidates[number_of_candidates].coord[2] = candidate_z;

          int n_candidate = cut_cell_disc.get_n(i_candidate, j_candidate, k_candidate);
          if (n_candidate >= 0)
            {
              // Nothing to assert
            }
          else
            {
              assert((!next_time) || (cut_cell_disc.get_interfaces().In(i_candidate,j_candidate,k_candidate) == (double)phase));
              assert((next_time) || (cut_cell_disc.get_interfaces().I(i_candidate,j_candidate,k_candidate) == (double)phase));
            }
          candidates[number_of_candidates].signed_independent_index = cut_cell_disc.get_splitting().get_signed_independent_index(phase, i_candidate, j_candidate, k_candidate);
          number_of_candidates += 1;
        }

    }


  assert(number_of_candidates <= max_number_of_candidates);
  assert(number_of_candidates <= max_number_of_cell_candidates);
  std::sort(candidates, candidates + number_of_candidates, [&](const Candidate& c1, const Candidate& c2)
  {
    return (c1.dist < c2.dist);
  });

  // On boucle d'abord sur n_neighbours, le nombre de points que l'on s'autorise a chercher
  // dans la liste des voisins. Par exemple, n_neighbours=6 veut dire que l'on cherche a former
  // des tetraedres a partir des 6 premiers voisins.
  const bool limit_count = false;
  const int max_count = 100;
  int count = 0;
  int closest_tetrahedron[4] = {-1,-1,-1,-1};
  double closest_lambda_error = DMAXFLOAT;
  for (int n_neighbours = 4; (n_neighbours < number_of_candidates && ((!limit_count) || count < max_count)); n_neighbours++)
    {
      int index_vertex0 = n_neighbours; // Le premier point est toujours fixe au dernier possible. Cela garanti que tous les tetraedres d'un nouveau n_neighbours sont nouveaux.

      for (int index_vertex1 = 0; (index_vertex1 < index_vertex0-2 && ((!limit_count) || count < max_count)); index_vertex1++)
        {
          for (int index_vertex2 = index_vertex1+1; (index_vertex2 < index_vertex0-1 && ((!limit_count) || count < max_count)); index_vertex2++)
            {
              for (int index_vertex3 = index_vertex2+1; (index_vertex3 < index_vertex0 && ((!limit_count) || count < max_count)); index_vertex3++)
                {
                  Vecteur3 lambda_vec = compute_lambda(index_vertex0, index_vertex1, index_vertex2, index_vertex3, candidates, xfact, yfact, zfact);
                  double lambda_1 = lambda_vec[0];
                  double lambda_2 = lambda_vec[1];
                  double lambda_3 = lambda_vec[2];

                  double lambda_0 = 1 - lambda_1 - lambda_2 - lambda_3;

                  int lambda_0_within_bounds = ((lambda_0 >= 0) && (lambda_0 <= 1));
                  int lambda_1_within_bounds = ((lambda_1 >= 0) && (lambda_1 <= 1));
                  int lambda_2_within_bounds = ((lambda_2 >= 0) && (lambda_2 <= 1));
                  int lambda_3_within_bounds = ((lambda_3 >= 0) && (lambda_3 <= 1));

                  int point_within_tetrahedron = lambda_0_within_bounds && lambda_1_within_bounds && lambda_2_within_bounds && lambda_3_within_bounds;

                  count++;
                  if (point_within_tetrahedron)
                    {
                      int signed_independent_index_0 = candidates[index_vertex0].signed_independent_index;
                      int signed_independent_index_1 = candidates[index_vertex1].signed_independent_index;
                      int signed_independent_index_2 = candidates[index_vertex2].signed_independent_index;
                      int signed_independent_index_3 = candidates[index_vertex3].signed_independent_index;

                      tetrahedra_index_coefficient[0] = {signed_independent_index_0, lambda_0};
                      tetrahedra_index_coefficient[1] = {signed_independent_index_1, lambda_1};
                      tetrahedra_index_coefficient[2] = {signed_independent_index_2, lambda_2};
                      tetrahedra_index_coefficient[3] = {signed_independent_index_3, lambda_3};

                      status = count;
                      return;
                    }
                  else
                    {
                      double lambda_error_0 = std::max((lambda_0 < 0)*(-lambda_0), (lambda_0 > 1)*(lambda_0 - 1));
                      double lambda_error_1 = std::max((lambda_1 < 0)*(-lambda_1), (lambda_1 > 1)*(lambda_1 - 1));
                      double lambda_error_2 = std::max((lambda_2 < 0)*(-lambda_2), (lambda_2 > 1)*(lambda_2 - 1));
                      double lambda_error_3 = std::max((lambda_3 < 0)*(-lambda_3), (lambda_3 > 1)*(lambda_3 - 1));

                      double lambda_error = std::max(lambda_error_0, std::max(lambda_error_1, std::max(lambda_error_2, lambda_error_3)));
                      assert(lambda_error > 0);

                      if (closest_lambda_error > lambda_error)
                        {
                          closest_tetrahedron[0] = index_vertex0;
                          closest_tetrahedron[1] = index_vertex1;
                          closest_tetrahedron[2] = index_vertex2;
                          closest_tetrahedron[3] = index_vertex3;
                          closest_lambda_error = lambda_error;
                        }
                    }
                }
            }
        }
    }

  // Utilise le tetrahedre le plus proche selon les parametres de tolerance
  if ((tolerate_not_within_tetrahedron == 2) || (tolerate_not_within_tetrahedron == 1 && closest_lambda_error < 1e-9))
    {
      Vecteur3 lambda_vec = compute_lambda(closest_tetrahedron[0], closest_tetrahedron[1], closest_tetrahedron[2], closest_tetrahedron[3], candidates, xfact, yfact, zfact);
      double lambda_1 = lambda_vec[0];
      double lambda_2 = lambda_vec[1];
      double lambda_3 = lambda_vec[2];

      double lambda_0 = 1 - lambda_1 - lambda_2 - lambda_3;

      assert(!(((lambda_0 >= 0) && (lambda_0 <= 1)) && ((lambda_1 >= 0) && (lambda_1 <= 1)) && ((lambda_2 >= 0) && (lambda_2 <= 1)) && ((lambda_3 >= 0) && (lambda_3 <= 1))));

      int signed_independent_index_0 = candidates[closest_tetrahedron[0]].signed_independent_index;
      int signed_independent_index_1 = candidates[closest_tetrahedron[1]].signed_independent_index;
      int signed_independent_index_2 = candidates[closest_tetrahedron[2]].signed_independent_index;
      int signed_independent_index_3 = candidates[closest_tetrahedron[3]].signed_independent_index;

      tetrahedra_index_coefficient[0] = {signed_independent_index_0, lambda_0};
      tetrahedra_index_coefficient[1] = {signed_independent_index_1, lambda_1};
      tetrahedra_index_coefficient[2] = {signed_independent_index_2, lambda_2};
      tetrahedra_index_coefficient[3] = {signed_independent_index_3, lambda_3};

      status = -1;
      return;
    }
  else
    {
      Cerr << "Value of close_to_edge: " << (int)close_to_edge << finl;
      Cerr << "Error in ijk_interpolate_cut_cell_for_given_index: no tetrahedron containing the point " << x << " " << y << " " << z << " on processor " << Process::me() << finl;
      Cerr << "For information, closest_lambda_error=" << closest_lambda_error << finl;
      // Note: While maybe not most likely, I noticed this error could occur if
      //  * the number of cells per particle diameter is too small <=3
      //  * ijk_splitting_ft_extension is not large enough
      Process::exit();
      tetrahedra_index_coefficient[0] = {-1,-1.};
      tetrahedra_index_coefficient[1] = {-1,-1.};
      tetrahedra_index_coefficient[2] = {-1,-1.};
      tetrahedra_index_coefficient[3] = {-1,-1.};
      return;
    }
}

static double ijk_interpolate_cut_cell_using_interface_for_given_index(bool next_time, int phase, const IJK_Field_double field_ft, const Cut_field_double& field, const ArrOfDouble& interfacial_temperature, const Vecteur3& coordinates, int tolerate_not_within_tetrahedron, int skip_unknown_points, double value_for_bad_points, int& status)
{
  if (Process::me() > 0)
    {
      Cerr << "Error in ijk_interpolate_cut_cell_using_interface_for_given_index: Le calcul est parallele mais cette methode n'est pas correcte dans le cas parallele" << finl;
      Process::exit();
    }

  const Cut_cell_FT_Disc& cut_cell_disc = field.get_cut_cell_disc();

  const Maillage_FT_IJK& mesh = next_time ? cut_cell_disc.get_interfaces().maillage_ft_ijk() : cut_cell_disc.get_interfaces().old_maillage_ft_ijk();
  const IntTab& facettes = mesh.facettes();
  const DoubleTab& sommets = mesh.sommets();
  const ArrOfDouble& surface_facettes = mesh.get_update_surface_facettes();
  const Intersections_Elem_Facettes& intersec = mesh.intersections_elem_facettes();

  //const int ghost = field.ghost();
  const int ghost = cut_cell_disc.get_ghost_size();
  assert(field.ghost() >= ghost);

  const double x = coordinates[0];
  const double y = coordinates[1];
  const double z = coordinates[2];

  const int ni_ft = field_ft.ni();
  const int nj_ft = field_ft.nj();
  const int nk_ft = field_ft.nk();

  const IJK_Splitting& splitting_ft = field_ft.get_splitting();
  const IJK_Grid_Geometry& geom_ft = splitting_ft.get_grid_geometry();

  const double dx_ft = geom_ft.get_constant_delta(DIRECTION_I);
  const double dy_ft = geom_ft.get_constant_delta(DIRECTION_J);
  const double dz_ft = geom_ft.get_constant_delta(DIRECTION_K);
  double origin_x_ft = geom_ft.get_origin(DIRECTION_I);
  double origin_y_ft = geom_ft.get_origin(DIRECTION_J);
  double origin_z_ft = geom_ft.get_origin(DIRECTION_K);

  const double x2_ft = (x - origin_x_ft) / dx_ft;
  const double y2_ft = (y - origin_y_ft) / dy_ft;
  const double z2_ft = (z - origin_z_ft) / dz_ft;

  const int index_i_ft = (int)(std::floor(x2_ft)) - splitting_ft.get_offset_local(DIRECTION_I);
  const int index_j_ft = (int)(std::floor(y2_ft)) - splitting_ft.get_offset_local(DIRECTION_J);
  const int index_k_ft = (int)(std::floor(z2_ft)) - splitting_ft.get_offset_local(DIRECTION_K);

  const int ni = field.ni();
  const int nj = field.nj();
  const int nk = field.nk();

  const IJK_Splitting& splitting = field.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  //const IJK_Splitting::Localisation loc = field.pure_.get_localisation();
  // L'origine est sur un noeud. Donc que la premiere face en I est sur get_origin(DIRECTION_I)
  double origin_x = geom.get_origin(DIRECTION_I);
  double origin_y = geom.get_origin(DIRECTION_J);
  double origin_z = geom.get_origin(DIRECTION_K);

  const double x2 = (x - origin_x) / dx;
  const double y2 = (y - origin_y) / dy;
  const double z2 = (z - origin_z) / dz;

  // Coordonnes barycentriques du points dans la cellule :
  const double xfact = x2 - floor(x2);
  const double yfact = y2 - floor(y2);
  const double zfact = z2 - floor(z2);

  // On travaille sur le maillage NS, on va donc corrige les indices de la periodicite.
  // Note : on ne corrige que l'index et pas les coordonnees, car on n'utilise plus les coordonnees par la suite.
  const int index_i = cut_cell_disc.get_splitting().get_i_along_dir_perio(DIRECTION_I, x, IJK_Splitting::ELEM);
  const int index_j = cut_cell_disc.get_splitting().get_i_along_dir_perio(DIRECTION_J, y, IJK_Splitting::ELEM);
  const int index_k = cut_cell_disc.get_splitting().get_i_along_dir_perio(DIRECTION_K, z, IJK_Splitting::ELEM);

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

  Sommet involved_sommet[max_number_of_involved_sommet] = {};
  for (int i = 0; i < max_number_of_involved_sommet; i++)
    {
      involved_sommet[i].sommet = -1;
      involved_sommet[i].fa7 = -1;
      involved_sommet[i].value = 0.;
      involved_sommet[i].count = 0;
    }
  int number_of_involved_sommet = 0;

  int number_of_candidates = 0;
  Candidate_with_direct_value candidates[max_number_of_candidates];
  for (int i = 0; i < max_number_of_cell_candidates; i++)
    {
      int i_candidate_aperio = index_i + candidate_offset[i][0];
      int j_candidate_aperio = index_j + candidate_offset[i][1];
      int k_candidate_aperio = index_k + candidate_offset[i][2];

      // Prise en compte de la periodicite
      int i_candidate = cut_cell_disc.get_splitting().correct_perio_i_local(0, i_candidate_aperio);
      int j_candidate = cut_cell_disc.get_splitting().correct_perio_i_local(1, j_candidate_aperio);
      int k_candidate = cut_cell_disc.get_splitting().correct_perio_i_local(2, k_candidate_aperio);

      double old_indicatrice = cut_cell_disc.get_interfaces().I(i_candidate, j_candidate, k_candidate);
      double next_indicatrice = cut_cell_disc.get_interfaces().In(i_candidate, j_candidate, k_candidate);
      double indicatrice = next_time ? next_indicatrice : old_indicatrice;
      if ((phase == 0 && indicatrice == 1.) || (phase == 1 && indicatrice == 0.))
        {
          // Point invalide
        }
      else
        {
          double candidate_x = (double)candidate_offset[i][0] + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 0, phase, i_candidate, j_candidate, k_candidate));
          double candidate_y = (double)candidate_offset[i][1] + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 1, phase, i_candidate, j_candidate, k_candidate));
          double candidate_z = (double)candidate_offset[i][2] + (cut_cell_disc.get_interfaces().get_barycentre(next_time, 2, phase, i_candidate, j_candidate, k_candidate));
          double decalage_x = candidate_x - xfact;
          double decalage_y = candidate_y - yfact;
          double decalage_z = candidate_z - zfact;
          candidates[number_of_candidates].index = number_of_candidates;
          candidates[number_of_candidates].dist = sqrt(decalage_x*decalage_x + decalage_y*decalage_y + decalage_z*decalage_z);
          candidates[number_of_candidates].coord[0] = candidate_x;
          candidates[number_of_candidates].coord[1] = candidate_y;
          candidates[number_of_candidates].coord[2] = candidate_z;

          int n_candidate = cut_cell_disc.get_n(i_candidate, j_candidate, k_candidate);
          if (n_candidate >= 0)
            {
              candidates[number_of_candidates].value = (phase == 0) ? field.diph_v_(n_candidate) : field.diph_l_(n_candidate);
            }
          else
            {
              assert((!next_time) || (cut_cell_disc.get_interfaces().In(i_candidate,j_candidate,k_candidate) == (double)phase));
              assert((next_time) || (cut_cell_disc.get_interfaces().I(i_candidate,j_candidate,k_candidate) == (double)phase));
              candidates[number_of_candidates].value = field.pure_(i_candidate,j_candidate,k_candidate);
            }
          assert(candidates[number_of_candidates].value != 0); // Suggests a bug, but not necessarily implies so
          number_of_candidates += 1;
        }

      int i_candidate_ft = index_i_ft + candidate_offset[i][0];
      int j_candidate_ft = index_j_ft + candidate_offset[i][1];
      int k_candidate_ft = index_k_ft + candidate_offset[i][2];

      if (i_candidate_ft < 0 || j_candidate_ft < 0 || k_candidate_ft < 0 || i_candidate_ft >= ni_ft || j_candidate_ft >= nj_ft || k_candidate_ft >= nk_ft)
        {
          continue; // Ne fait rien, car ces conditions semblent impliquer que le point est inutile.
        }

      const ArrOfInt& index_elem = intersec.index_elem();
      assert(mesh.ref_splitting().valeur() == splitting_ft);
      const int num_elem = splitting_ft.convert_ijk_cell_to_packed(i_candidate_ft,j_candidate_ft,k_candidate_ft);
      int index = index_elem[num_elem];

      // Boucle sur les facettes qui traversent cet element
      while (index >= 0)
        {
          const Intersections_Elem_Facettes_Data& data = intersec.data_intersection(index);
          const int fa7 = data.numero_facette_;

          for (int som = 0; som < 3; som++)
            {
              // Note : Incomplet
              // On ne prend en compte que les sommets reels
              // car il semble difficile de prendre en compte les facettes que le PE ne connait pas
              // et je suppose que l'on connait toutes les facettes associees a un sommet reel.
              // Cela veut dire qu'il y aura une difference entre sequentiel et parallele.
              if (!mesh.sommet_virtuel(facettes(fa7, som)))
                {
                  assert(number_of_involved_sommet < max_number_of_involved_sommet);
                  involved_sommet[number_of_involved_sommet].sommet = facettes(fa7, som);
                  involved_sommet[number_of_involved_sommet].fa7 = fa7;
                  involved_sommet[number_of_involved_sommet].value = interfacial_temperature(fa7)/surface_facettes(fa7);
                  involved_sommet[number_of_involved_sommet].count = 1;
                  number_of_involved_sommet += 1;
                }
            }

          index = data.index_facette_suivante_;
        };
    }

  std::sort(involved_sommet, involved_sommet + number_of_involved_sommet, [&](const Sommet& s1, const Sommet& s2)
  {
    if(s1.sommet == s2.sommet)
      return (s1.fa7 > s2.fa7);
    else
      return (s1.sommet > s2.sommet);
  });

  int initial_number_of_involved_sommet = number_of_involved_sommet;

  int precedente_fa7 = -1;
  int precedent_sommet = -1;
  int premier_i_avec_ce_sommet = -1;
  for (int i = 0; i < initial_number_of_involved_sommet; i++)
    {
      int sommet = involved_sommet[i].sommet;
      int fa7 = involved_sommet[i].fa7;
      assert(sommet != -1);
      if (sommet == precedent_sommet)
        {
          assert(fa7 != -1);
          if (fa7 == precedente_fa7)
            {
              // This is a duplicate, thus the information is destroyed (below).
            }
          else
            {
              // This is the contribution of a different facet,
              // the information is moved to the first instance of the vertex, then destroyed (below).
              involved_sommet[premier_i_avec_ce_sommet].value += involved_sommet[i].value;
              involved_sommet[premier_i_avec_ce_sommet].count += 1;

              precedente_fa7 = fa7;
            }

          involved_sommet[i].sommet = -1;
          involved_sommet[i].fa7 = -1;
          involved_sommet[i].value = 0.;
          involved_sommet[i].count = 1;

          number_of_involved_sommet -= 1;

        }
      else
        {
          premier_i_avec_ce_sommet = i;
          precedent_sommet = sommet;
          precedente_fa7 = -1;
        }
    }

  for (int i = 0; i < initial_number_of_involved_sommet; i++)
    {
      involved_sommet[i].value /= involved_sommet[i].count;
      involved_sommet[i].count = 1;
    }

  std::sort(involved_sommet, involved_sommet + initial_number_of_involved_sommet, [&](const Sommet& s1, const Sommet& s2)
  {
    if(s1.sommet == s2.sommet)
      return (s1.fa7 > s2.fa7);
    else
      return (s1.sommet > s2.sommet);
  });


  for (int i = 0; i < number_of_involved_sommet; i++)
    {
      int sommet = involved_sommet[i].sommet;
      if (i >= 1)
        {
          assert(sommet != involved_sommet[i-1].sommet);
        }

      double decalage_x = (sommets(sommet, 0) - x)/dx;
      double decalage_y = (sommets(sommet, 1) - y)/dy;
      double decalage_z = (sommets(sommet, 2) - z)/dz;
      double candidate_x = decalage_x + xfact;
      double candidate_y = decalage_y + yfact;
      double candidate_z = decalage_z + zfact;

      candidates[number_of_candidates].dist = sqrt(decalage_x*decalage_x + decalage_y*decalage_y + decalage_z*decalage_z);
      candidates[number_of_candidates].coord[0] = candidate_x;
      candidates[number_of_candidates].coord[1] = candidate_y;
      candidates[number_of_candidates].coord[2] = candidate_z;
      candidates[number_of_candidates].value = involved_sommet[i].value;
      number_of_candidates += 1;
    }

  assert(number_of_candidates <= max_number_of_candidates);
  std::sort(candidates, candidates + number_of_candidates, [&](const Candidate_with_direct_value& c1, const Candidate_with_direct_value& c2)
  {
    return (c1.dist < c2.dist);
  });

  // On boucle d'abord sur n_neighbours, le nombre de points que l'on s'autorise a chercher
  // dans la liste des voisins. Par exemple, n_neighbours=6 veut dire que l'on cherche a former
  // des tetraedres a partir des 6 premiers voisins.
  const bool limit_count = false;
  const bool force_use_fist_neighbour = true; // Force l'utilisation du voisin le plus proche. Avantage : garanti de rester proche de ce point si le premier voisin est tres proche.
  const int max_count = 100;
  int count = 0;
  int closest_tetrahedron[4] = {-1,-1,-1,-1};
  double closest_lambda_error = DMAXFLOAT;
  for (int n_neighbours = 4; (n_neighbours < number_of_candidates && ((!limit_count) || count < max_count)); n_neighbours++)
    {
      int index_vertex0 = n_neighbours; // Le premier point est toujours fixe au dernier possible. Cela garanti que tous les tetraedres d'un nouveau n_neighbours sont nouveaux.

      for (int index_vertex1 = 0; (index_vertex1 < (force_use_fist_neighbour ? 1 : index_vertex0-2) && ((!limit_count) || count < max_count)); index_vertex1++)
        {
          for (int index_vertex2 = index_vertex1+1; (index_vertex2 < index_vertex0-1 && ((!limit_count) || count < max_count)); index_vertex2++)
            {
              for (int index_vertex3 = index_vertex2+1; (index_vertex3 < index_vertex0 && ((!limit_count) || count < max_count)); index_vertex3++)
                {
                  Vecteur3 lambda_vec = compute_lambda(index_vertex0, index_vertex1, index_vertex2, index_vertex3, candidates, xfact, yfact, zfact);
                  double lambda_1 = lambda_vec[0];
                  double lambda_2 = lambda_vec[1];
                  double lambda_3 = lambda_vec[2];

                  double lambda_0 = 1 - lambda_1 - lambda_2 - lambda_3;

                  int lambda_0_within_bounds = ((lambda_0 >= 0) && (lambda_0 <= 1));
                  int lambda_1_within_bounds = ((lambda_1 >= 0) && (lambda_1 <= 1));
                  int lambda_2_within_bounds = ((lambda_2 >= 0) && (lambda_2 <= 1));
                  int lambda_3_within_bounds = ((lambda_3 >= 0) && (lambda_3 <= 1));

                  int point_within_tetrahedron = lambda_0_within_bounds && lambda_1_within_bounds && lambda_2_within_bounds && lambda_3_within_bounds;

                  count++;
                  if (point_within_tetrahedron)
                    {
                      double field_0 = candidates[index_vertex0].value;
                      double field_1 = candidates[index_vertex1].value;
                      double field_2 = candidates[index_vertex2].value;
                      double field_3 = candidates[index_vertex3].value;
                      assert(field_0 != 0); // Suggests a bug, but not necessarily implies so
                      assert(field_1 != 0); //  .
                      assert(field_2 != 0); //  .
                      assert(field_3 != 0); //  .

                      double r = field_0*lambda_0 + field_1*lambda_1 + field_2*lambda_2 + field_3*lambda_3;

                      status = count;
                      return r;
                    }
                  else
                    {
                      double lambda_error_0 = std::max((lambda_0 < 0)*(-lambda_0), (lambda_0 > 1)*(lambda_0 - 1));
                      double lambda_error_1 = std::max((lambda_1 < 0)*(-lambda_1), (lambda_1 > 1)*(lambda_1 - 1));
                      double lambda_error_2 = std::max((lambda_2 < 0)*(-lambda_2), (lambda_2 > 1)*(lambda_2 - 1));
                      double lambda_error_3 = std::max((lambda_3 < 0)*(-lambda_3), (lambda_3 > 1)*(lambda_3 - 1));

                      double lambda_error = std::max(lambda_error_0, std::max(lambda_error_1, std::max(lambda_error_2, lambda_error_3)));
                      assert(lambda_error > 0);

                      if (closest_lambda_error > lambda_error)
                        {
                          closest_tetrahedron[0] = index_vertex0;
                          closest_tetrahedron[1] = index_vertex1;
                          closest_tetrahedron[2] = index_vertex2;
                          closest_tetrahedron[3] = index_vertex3;
                          closest_lambda_error = lambda_error;
                        }
                    }
                }
            }
        }
    }

  // Utilise le tetrahedre le plus proche selon les parametres de tolerance
  if ((tolerate_not_within_tetrahedron == 2) || (tolerate_not_within_tetrahedron == 1 && closest_lambda_error < 5e-2))
    {
      Vecteur3 lambda_vec = compute_lambda(closest_tetrahedron[0], closest_tetrahedron[1], closest_tetrahedron[2], closest_tetrahedron[3], candidates, xfact, yfact, zfact);
      double lambda_1 = lambda_vec[0];
      double lambda_2 = lambda_vec[1];
      double lambda_3 = lambda_vec[2];

      double lambda_0 = 1 - lambda_1 - lambda_2 - lambda_3;

      assert(!(((lambda_0 >= 0) && (lambda_0 <= 1)) && ((lambda_1 >= 0) && (lambda_1 <= 1)) && ((lambda_2 >= 0) && (lambda_2 <= 1)) && ((lambda_3 >= 0) && (lambda_3 <= 1))));

      double field_0 = candidates[closest_tetrahedron[0]].value;
      double field_1 = candidates[closest_tetrahedron[1]].value;
      double field_2 = candidates[closest_tetrahedron[2]].value;
      double field_3 = candidates[closest_tetrahedron[3]].value;

      double r = field_0*lambda_0 + field_1*lambda_1 + field_2*lambda_2 + field_3*lambda_3;

      status = -1;
      return r;
    }
  else
    {
      Cerr << "Value of close_to_edge: " << (int)close_to_edge << finl;
      Cerr << "Error in ijk_interpolate_cut_cell_using_interface_for_given_index: no tetrahedron containing the point " << x << " " << y << " " << z << " on processor " << Process::me() << finl;
      Cerr << "For information, closest_lambda_error=" << closest_lambda_error << finl;
      // Note: While maybe not most likely, I noticed this error could occur if
      //  * the number of cells per particle diameter is too small <=3
      //  * ijk_splitting_ft_extension is not large enough
      Process::exit();
      return -1;
    }
}

// Interpolate the "field" at the requested "coordinates" (array with 3 columns), and stores into "result"
static void ijk_interpolate_cut_cell_implementation(bool next_time, int phase, const Cut_field_double& field, const DoubleTab& coordinates, ArrOfDouble& result, int skip_unknown_points, double value_for_bad_points)
{
  const int nb_coords = coordinates.dimension(0);
  result.resize_array(nb_coords);
  for (int idx = 0; idx < nb_coords; idx++)
    {
      Vecteur3 coordinates_for_given_index(coordinates(idx,0), coordinates(idx,1), coordinates(idx,2));
      int tolerate_not_within_tetrahedron = 1;
      int status = -2;
      IndexCoefficient tetrahedra_index_coefficient[4];
      ijk_interpolate_cut_cell_for_given_index(next_time, phase, field.get_cut_cell_disc(), coordinates_for_given_index, tetrahedra_index_coefficient, tolerate_not_within_tetrahedron, skip_unknown_points, value_for_bad_points, status);

      double field_0 = field.from_signed_independent_index(tetrahedra_index_coefficient[0].index);
      double field_1 = field.from_signed_independent_index(tetrahedra_index_coefficient[1].index);
      double field_2 = field.from_signed_independent_index(tetrahedra_index_coefficient[2].index);
      double field_3 = field.from_signed_independent_index(tetrahedra_index_coefficient[3].index);

      double interpolated_value = field_0*tetrahedra_index_coefficient[0].coefficient + field_1*tetrahedra_index_coefficient[1].coefficient + field_2*tetrahedra_index_coefficient[2].coefficient + field_3*tetrahedra_index_coefficient[3].coefficient;
      result[idx] = interpolated_value;
    }
}

static void ijk_interpolate_cut_cell_implementation(bool next_time, int phase, const Cut_cell_FT_Disc& cut_cell_disc, const DoubleTab& coordinates, IntTabFT& signed_independent_index, DoubleTabFT& coefficient, int skip_unknown_points, double value_for_bad_points)
{
  const int nb_coords = coordinates.dimension(0);
  signed_independent_index.resize(nb_coords, 4);
  coefficient.resize(nb_coords, 4);
  for (int idx = 0; idx < nb_coords; idx++)
    {
      Vecteur3 coordinates_for_given_index(coordinates(idx,0), coordinates(idx,1), coordinates(idx,2));
      int tolerate_not_within_tetrahedron = 1;
      int status = -2;
      IndexCoefficient tetrahedra_index_coefficient[4];
      ijk_interpolate_cut_cell_for_given_index(next_time, phase, cut_cell_disc, coordinates_for_given_index, tetrahedra_index_coefficient, tolerate_not_within_tetrahedron, skip_unknown_points, value_for_bad_points, status);

      signed_independent_index(idx,0) = tetrahedra_index_coefficient[0].index;
      signed_independent_index(idx,1) = tetrahedra_index_coefficient[1].index;
      signed_independent_index(idx,2) = tetrahedra_index_coefficient[2].index;
      signed_independent_index(idx,3) = tetrahedra_index_coefficient[3].index;
      coefficient(idx,0) = tetrahedra_index_coefficient[0].coefficient;
      coefficient(idx,1) = tetrahedra_index_coefficient[1].coefficient;
      coefficient(idx,2) = tetrahedra_index_coefficient[2].coefficient;
      coefficient(idx,3) = tetrahedra_index_coefficient[3].coefficient;
    }
}

void ijk_interpolate_cut_cell_skip_unknown_points(bool next_time, int phase, const Cut_field_double& field, const DoubleTab& coordinates, ArrOfDouble& result, const double value_for_bad_points)
{
  ijk_interpolate_cut_cell_implementation(next_time, phase, field, coordinates, result, 1 /* yes:skip unknown points */, value_for_bad_points);
}

void ijk_interpolate_cut_cell(bool next_time, int phase, const Cut_field_double& field, const DoubleTab& coordinates, ArrOfDouble& result)
{
  ijk_interpolate_cut_cell_implementation(next_time, phase, field, coordinates, result, 0 /* skip unknown points=no */, 0.);
}

void ijk_interpolate_cut_cell_skip_unknown_points(bool next_time, int phase, const Cut_cell_FT_Disc& cut_cell_disc, const DoubleTab& coordinates, IntTabFT& signed_independent_index, DoubleTabFT& coefficient, const double value_for_bad_points)
{
  ijk_interpolate_cut_cell_implementation(next_time, phase, cut_cell_disc, coordinates, signed_independent_index, coefficient, 1 /* yes:skip unknown points */, value_for_bad_points);
}

void ijk_interpolate_cut_cell(bool next_time, int phase, const Cut_cell_FT_Disc& cut_cell_disc, const DoubleTab& coordinates, IntTabFT& signed_independent_index, DoubleTabFT& coefficient)
{
  ijk_interpolate_cut_cell_implementation(next_time, phase, cut_cell_disc, coordinates, signed_independent_index, coefficient, 0 /* skip unknown points=no */, 0.);
}

double ijk_interpolate_cut_cell_using_interface_skip_unknown_points(bool next_time, int phase, const IJK_Field_double field_ft, const Cut_field_double& field, const ArrOfDouble& interfacial_temperature, const Vecteur3& coordinates, int tolerate_not_within_tetrahedron, const double value_for_bad_points, int& status)
{
  double interpolated_value = ijk_interpolate_cut_cell_using_interface_for_given_index(next_time, phase, field_ft, field, interfacial_temperature, coordinates, tolerate_not_within_tetrahedron, 1, value_for_bad_points, status);
  return interpolated_value;
}

double ijk_interpolate_cut_cell_using_interface(bool next_time, int phase, const IJK_Field_double field_ft, const Cut_field_double& field, const ArrOfDouble& interfacial_temperature, const Vecteur3& coordinates, int tolerate_not_within_tetrahedron, int& status)
{
  double interpolated_value = ijk_interpolate_cut_cell_using_interface_for_given_index(next_time, phase, field_ft, field, interfacial_temperature, coordinates, tolerate_not_within_tetrahedron, 0, 0., status);
  return interpolated_value;
}

// Returns the indices and coefficients to interpolate the "field" at the requested coordinate
static void ijk_interpolate_one_value(bool next_time, const Cut_cell_FT_Disc& cut_cell_disc, const Vecteur3& coordinates, IntTabFT& signed_independent_index, DoubleTabFT& coefficient, int skip_unknown_points, double value_for_bad_points)
{
  const int ghost = cut_cell_disc.get_ghost_size();
  const IJK_Splitting& splitting = cut_cell_disc.get_splitting();
  const int ni = splitting.get_nb_items_local(IJK_Splitting::ELEM, 0);
  const int nj = splitting.get_nb_items_local(IJK_Splitting::ELEM, 1);
  const int nk = splitting.get_nb_items_local(IJK_Splitting::ELEM, 2);

  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  // L'origine est sur un noeud. Donc que la premiere face en I est sur get_origin(DIRECTION_I)
  double origin_x = geom.get_origin(DIRECTION_I);
  double origin_y = geom.get_origin(DIRECTION_J);
  double origin_z = geom.get_origin(DIRECTION_K);
  const double x = coordinates[0];
  const double y = coordinates[1];
  const double z = coordinates[2];
  const double x2 = (x - origin_x) / dx;
  const double y2 = (y - origin_y) / dy;
  const double z2 = (z - origin_z) / dz;
  const int index_i = (int) (floor(x2)) - splitting.get_offset_local(DIRECTION_I);
  const int index_j = (int) (floor(y2)) - splitting.get_offset_local(DIRECTION_J);
  const int index_k = (int) (floor(z2)) - splitting.get_offset_local(DIRECTION_K);
  // Coordonnes barycentriques du points dans la cellule :
  const double xfact = x2 - floor(x2);
  const double yfact = y2 - floor(y2);
  const double zfact = z2 - floor(z2);

  // is point in the domain ? (ghost cells ok...)
  bool ok = (index_i >= -ghost && index_i < ni + ghost - 1) && (index_j >= -ghost && index_j < nj + ghost - 1) && (index_k >= -ghost && index_k < nk + ghost - 1);
  if (!ok)
    {
      if (skip_unknown_points)
        {
        }
      else
        {
          // Error!
          Cerr << "Error in ijk_interpolate_implementation: request interpolation of point " << x << " " << y << " " << z << " which is outside of the domain on processor " << Process::me()
               << finl;
          Process::exit();
        }
    }

  double indic_0  = next_time ? (1 - cut_cell_disc.get_interfaces().In_ft(index_i,   index_j,   index_k))   : (1 - cut_cell_disc.get_interfaces().I_ft(index_i,   index_j,   index_k));
  double indic_1  = next_time ? (1 - cut_cell_disc.get_interfaces().In_ft(index_i+1, index_j,   index_k))   : (1 - cut_cell_disc.get_interfaces().I_ft(index_i+1, index_j,   index_k));
  double indic_2  = next_time ? (1 - cut_cell_disc.get_interfaces().In_ft(index_i,   index_j+1, index_k))   : (1 - cut_cell_disc.get_interfaces().I_ft(index_i,   index_j+1, index_k));
  double indic_3  = next_time ? (1 - cut_cell_disc.get_interfaces().In_ft(index_i+1, index_j+1, index_k))   : (1 - cut_cell_disc.get_interfaces().I_ft(index_i+1, index_j+1, index_k));
  double indic_4  = next_time ? (1 - cut_cell_disc.get_interfaces().In_ft(index_i,   index_j,   index_k+1)) : (1 - cut_cell_disc.get_interfaces().I_ft(index_i,   index_j,   index_k+1));
  double indic_5  = next_time ? (1 - cut_cell_disc.get_interfaces().In_ft(index_i+1, index_j,   index_k+1)) : (1 - cut_cell_disc.get_interfaces().I_ft(index_i+1, index_j,   index_k+1));
  double indic_6  = next_time ? (1 - cut_cell_disc.get_interfaces().In_ft(index_i,   index_j+1, index_k+1)) : (1 - cut_cell_disc.get_interfaces().I_ft(index_i,   index_j+1, index_k+1));
  double indic_7  = next_time ? (1 - cut_cell_disc.get_interfaces().In_ft(index_i+1, index_j+1, index_k+1)) : (1 - cut_cell_disc.get_interfaces().I_ft(index_i+1, index_j+1, index_k+1));
  double indic_8  = next_time ? cut_cell_disc.get_interfaces().In_ft(index_i,   index_j,   index_k)         : cut_cell_disc.get_interfaces().I_ft(index_i,   index_j,   index_k);
  double indic_9  = next_time ? cut_cell_disc.get_interfaces().In_ft(index_i+1, index_j,   index_k)         : cut_cell_disc.get_interfaces().I_ft(index_i+1, index_j,   index_k);
  double indic_10 = next_time ? cut_cell_disc.get_interfaces().In_ft(index_i,   index_j+1, index_k)         : cut_cell_disc.get_interfaces().I_ft(index_i,   index_j+1, index_k);
  double indic_11 = next_time ? cut_cell_disc.get_interfaces().In_ft(index_i+1, index_j+1, index_k)         : cut_cell_disc.get_interfaces().I_ft(index_i+1, index_j+1, index_k);
  double indic_12 = next_time ? cut_cell_disc.get_interfaces().In_ft(index_i,   index_j,   index_k+1)       : cut_cell_disc.get_interfaces().I_ft(index_i,   index_j,   index_k+1);
  double indic_13 = next_time ? cut_cell_disc.get_interfaces().In_ft(index_i+1, index_j,   index_k+1)       : cut_cell_disc.get_interfaces().I_ft(index_i+1, index_j,   index_k+1);
  double indic_14 = next_time ? cut_cell_disc.get_interfaces().In_ft(index_i,   index_j+1, index_k+1)       : cut_cell_disc.get_interfaces().I_ft(index_i,   index_j+1, index_k+1);
  double indic_15 = next_time ? cut_cell_disc.get_interfaces().In_ft(index_i+1, index_j+1, index_k+1)       : cut_cell_disc.get_interfaces().I_ft(index_i+1, index_j+1, index_k+1);

  int output_index = -1;

  if (indic_0 != 0)
    {
      output_index += 1;
      signed_independent_index(output_index)  = splitting.get_signed_independent_index(0, index_i,   index_j,   index_k);
      coefficient(output_index)               = indic_0  * (1.-xfact) * (1.-yfact) * (1.-zfact);
    }

  if (indic_1 != 0)
    {
      output_index += 1;
      signed_independent_index(output_index)  = splitting.get_signed_independent_index(0, index_i+1, index_j,   index_k);
      coefficient(output_index)               = indic_1  * (xfact)    * (1.-yfact) * (1.-zfact);
    }

  if (indic_2 != 0)
    {
      output_index += 1;
      signed_independent_index(output_index)  = splitting.get_signed_independent_index(0, index_i,   index_j+1, index_k);
      coefficient(output_index)               = indic_2  * (1.-xfact) * (yfact)    * (1.-zfact);
    }

  if (indic_3 != 0)
    {
      output_index += 1;
      signed_independent_index(output_index)  = splitting.get_signed_independent_index(0, index_i+1, index_j+1, index_k);
      coefficient(output_index)               = indic_3  * (xfact)    * (yfact)    * (1.-zfact);
    }

  if (indic_4 != 0)
    {
      output_index += 1;
      signed_independent_index(output_index)  = splitting.get_signed_independent_index(0, index_i,   index_j,   index_k+1);
      coefficient(output_index)               = indic_4  * (1.-xfact) * (1.-yfact) * (zfact);
    }

  if (indic_5 != 0)
    {
      output_index += 1;
      signed_independent_index(output_index)  = splitting.get_signed_independent_index(0, index_i+1, index_j,   index_k+1);
      coefficient(output_index)               = indic_5  * (xfact)    * (1.-yfact) * (zfact);
    }

  if (indic_6 != 0)
    {
      output_index += 1;
      signed_independent_index(output_index)  = splitting.get_signed_independent_index(0, index_i,   index_j+1, index_k+1);
      coefficient(output_index)               = indic_6  * (1.-xfact) * (yfact)    * (zfact);
    }

  if (indic_7 != 0)
    {
      output_index += 1;
      signed_independent_index(output_index)  = splitting.get_signed_independent_index(0, index_i+1, index_j+1, index_k+1);
      coefficient(output_index)               = indic_7  * (xfact)    * (yfact)    * (zfact);
    }

  if (indic_8 != 0)
    {
      output_index += 1;
      signed_independent_index(output_index)  = splitting.get_signed_independent_index(1, index_i,   index_j,   index_k);
      coefficient(output_index)               = indic_8  * (1.-xfact) * (1.-yfact) * (1.-zfact);
    }

  if (indic_9 != 0)
    {
      output_index += 1;
      signed_independent_index(output_index)  = splitting.get_signed_independent_index(1, index_i+1, index_j,   index_k);
      coefficient(output_index)               = indic_9  * (xfact)    * (1.-yfact) * (1.-zfact);
    }

  if (indic_10 != 0)
    {
      output_index += 1;
      signed_independent_index(output_index) = splitting.get_signed_independent_index(1, index_i,   index_j+1, index_k);
      coefficient(output_index)              = indic_10 * (1.-xfact) * (yfact)    * (1.-zfact);
    }

  if (indic_11 != 0)
    {
      output_index += 1;
      signed_independent_index(output_index) = splitting.get_signed_independent_index(1, index_i+1, index_j+1, index_k);
      coefficient(output_index)              = indic_11 * (xfact)    * (yfact)    * (1.-zfact);
    }

  if (indic_12 != 0)
    {
      output_index += 1;
      signed_independent_index(output_index) = splitting.get_signed_independent_index(1, index_i,   index_j,   index_k+1);
      coefficient(output_index)              = indic_12 * (1.-xfact) * (1.-yfact) * (zfact);
    }

  if (indic_13 != 0)
    {
      output_index += 1;
      signed_independent_index(output_index) = splitting.get_signed_independent_index(1, index_i+1, index_j,   index_k+1);
      coefficient(output_index)              = indic_13 * (xfact)    * (1.-yfact) * (zfact);
    }

  if (indic_14 != 0)
    {
      output_index += 1;
      signed_independent_index(output_index) = splitting.get_signed_independent_index(1, index_i,   index_j+1, index_k+1);
      coefficient(output_index)              = indic_14 * (1.-xfact) * (yfact)    * (zfact);
    }

  if (indic_15 != 0)
    {
      output_index += 1;
      signed_independent_index(output_index) = splitting.get_signed_independent_index(1, index_i+1, index_j+1, index_k+1);
      coefficient(output_index)              = indic_15 * (xfact)    * (yfact)    * (zfact);
    }

  signed_independent_index.resize(output_index+1);
  coefficient.resize(output_index+1);
}

void ijk_interpolate_skip_unknown_points(bool next_time, int phase, const Cut_cell_FT_Disc& cut_cell_disc, const Vecteur3& coordinates, IntTabFT& signed_independent_index, DoubleTabFT& coefficient, const double value_for_bad_points)
{
  signed_independent_index.resize(16);
  coefficient.resize(16);
  ijk_interpolate_one_value(next_time, cut_cell_disc, coordinates, signed_independent_index, coefficient, 1 /* yes:skip unknown points */, value_for_bad_points);
}

void ijk_interpolate(bool next_time, int phase, const Cut_cell_FT_Disc& cut_cell_disc, const Vecteur3& coordinates, IntTabFT& signed_independent_index, DoubleTabFT& coefficient)
{
  signed_independent_index.resize(16);
  coefficient.resize(16);
  ijk_interpolate_one_value(next_time, cut_cell_disc, coordinates, signed_independent_index, coefficient, 0 /* skip unknown points=no */, 0.);
}

void euler_explicit_update_cut_cell_notransport(double timestep, bool next_time, const Cut_field_double& dv, Cut_field_double& v)
{
  const Cut_cell_FT_Disc& cut_cell_disc = v.get_cut_cell_disc();
  const double delta_t = timestep;
  const int imax = v.ni();
  const int jmax = v.nj();
  const int kmax = v.nk();
  for (int k = 0; k < kmax; k++)
    {
      for (int j = 0; j < jmax; j++)
        {
          for (int i = 0; i < imax; i++)
            {
              double nonzero_indicatrice_l = next_time ? cut_cell_disc.get_interfaces().In_nonzero(1,i,j,k) : cut_cell_disc.get_interfaces().I_nonzero(1,i,j,k);
              double nonzero_indicatrice_v = next_time ? cut_cell_disc.get_interfaces().In_nonzero(0,i,j,k) : cut_cell_disc.get_interfaces().I_nonzero(0,i,j,k);

              int n = cut_cell_disc.get_n(i,j,k);
              if (n < 0)
                {
                  double x = dv.pure_(i,j,k);
                  double next_v_vol = v.pure_(i,j,k) + x * delta_t;
                  v.pure_(i,j,k) = next_v_vol;
                }
              else
                {
                  double x_l = dv.diph_l_(n);
                  double x_v = dv.diph_v_(n);

                  double next_v_vol_l = v.diph_l_(n)*nonzero_indicatrice_l + x_l * delta_t;
                  double next_v_vol_v = v.diph_v_(n)*nonzero_indicatrice_v + x_v * delta_t;
                  v.diph_l_(n) = next_v_vol_l/nonzero_indicatrice_l;
                  v.diph_v_(n) = next_v_vol_v/nonzero_indicatrice_v;
                }
            }
        }
    }
}

void runge_kutta3_update_cut_cell_notransport(bool next_time, const Cut_field_double& dv, Cut_field_double& F, Cut_field_double& v, const int step, double dt_tot, const Cut_field_int& cellule_rk_restreint)
{
  const double coeff_a[3] = { 0., -5. / 9., -153. / 128. };
  // Fk[0] = 1; Fk[i+1] = Fk[i] * a[i+1] + 1
  const double coeff_Fk[3] = { 1., 4. / 9., 15. / 32. };

  const double facteurF = coeff_a[step];
  const double intermediate_dt = compute_fractionnal_timestep_rk3(dt_tot, step);
  const double delta_t_divided_by_Fk = intermediate_dt / coeff_Fk[step];
  const int imax = v.ni();
  const int jmax = v.nj();
  const int kmax = v.nk();
  const Cut_cell_FT_Disc& cut_cell_disc = v.get_cut_cell_disc();
  switch(step)
    {
    case 0:
      // don't read initial value of F (no performance benefit because write to F causes the
      // processor to fetch the cache line, but we don't wand to use a potentially uninitialized value
      for (int k = 0; k < kmax; k++)
        {
          for (int j = 0; j < jmax; j++)
            {
              for (int i = 0; i < imax; i++)
                {
                  double nonzero_indicatrice_l = next_time ? cut_cell_disc.get_interfaces().In_nonzero(1,i,j,k) : cut_cell_disc.get_interfaces().I_nonzero(1,i,j,k);
                  double nonzero_indicatrice_v = next_time ? cut_cell_disc.get_interfaces().In_nonzero(0,i,j,k) : cut_cell_disc.get_interfaces().I_nonzero(0,i,j,k);

                  int n = cut_cell_disc.get_n(i,j,k);
                  if (n < 0)
                    {
                      double x = dv.pure_(i,j,k);
                      double next_v_vol = v.pure_(i,j,k) + x * delta_t_divided_by_Fk;
                      F.pure_(i,j,k) = x;
                      v.pure_(i,j,k) = next_v_vol;
                    }
                  else
                    {
                      double x_l = dv.diph_l_(n);
                      double x_v = dv.diph_v_(n);

                      double next_v_vol_l = v.diph_l_(n)*nonzero_indicatrice_l + x_l * delta_t_divided_by_Fk;
                      double next_v_vol_v = v.diph_v_(n)*nonzero_indicatrice_v + x_v * delta_t_divided_by_Fk;
                      F.diph_l_(n) = x_l;
                      F.diph_v_(n) = x_v;
                      v.diph_l_(n) = next_v_vol_l/nonzero_indicatrice_l;
                      v.diph_v_(n) = next_v_vol_v/nonzero_indicatrice_v;
                    }
                }
            }
        }
      break;
    case 1:
      // general case, read and write F
      for (int k = 0; k < kmax; k++)
        {
          for (int j = 0; j < jmax; j++)
            {
              for (int i = 0; i < imax; i++)
                {
                  double nonzero_indicatrice_l = next_time ? cut_cell_disc.get_interfaces().In_nonzero(1,i,j,k) : cut_cell_disc.get_interfaces().I_nonzero(1,i,j,k);
                  double nonzero_indicatrice_v = next_time ? cut_cell_disc.get_interfaces().In_nonzero(0,i,j,k) : cut_cell_disc.get_interfaces().I_nonzero(0,i,j,k);

                  int n = cut_cell_disc.get_n(i,j,k);
                  if (n < 0)
                    {
                      if (cellule_rk_restreint.pure_(i,j,k) == 0)
                        {
                          double x = F.pure_(i, j, k) * facteurF + dv.pure_(i,j,k);
                          double next_v_vol = v.pure_(i,j,k) + x * delta_t_divided_by_Fk;
                          F.pure_(i,j,k) = x;
                          v.pure_(i,j,k) = next_v_vol;
                        }
                      else
                        {
                          double x = dv.pure_(i,j,k);
                          double next_v_vol = v.pure_(i,j,k) + x * intermediate_dt;
                          F.pure_(i,j,k) = x;
                          v.pure_(i,j,k) = next_v_vol;
                        }
                    }
                  else
                    {
                      if (cellule_rk_restreint.diph_v_(n) == 0)
                        {
                          double x_v = F.diph_v_(n) * facteurF + dv.diph_v_(n);

                          double next_v_vol_v = v.diph_v_(n)*nonzero_indicatrice_v + x_v * delta_t_divided_by_Fk;
                          F.diph_v_(n) = x_v;
                          v.diph_v_(n) = next_v_vol_v/nonzero_indicatrice_v;
                        }
                      else
                        {
                          double x_v = dv.diph_v_(n);

                          double next_v_vol_v = v.diph_v_(n)*nonzero_indicatrice_v + x_v * intermediate_dt;
                          F.diph_v_(n) = x_v;
                          v.diph_v_(n) = next_v_vol_v/nonzero_indicatrice_v;
                        }

                      if (cellule_rk_restreint.diph_l_(n) == 0)
                        {
                          double x_l = F.diph_l_(n) * facteurF + dv.diph_l_(n);

                          double next_v_vol_l = v.diph_l_(n)*nonzero_indicatrice_l + x_l * delta_t_divided_by_Fk;
                          F.diph_l_(n) = x_l;
                          v.diph_l_(n) = next_v_vol_l/nonzero_indicatrice_l;
                        }
                      else
                        {
                          double x_l = dv.diph_l_(n);

                          double next_v_vol_l = v.diph_l_(n)*nonzero_indicatrice_l + x_l * intermediate_dt;
                          F.diph_l_(n) = x_l;
                          v.diph_l_(n) = next_v_vol_l/nonzero_indicatrice_l;
                        }
                    }
                }
            }
        }
      break;
    case 2:
      // do not write F
      for (int k = 0; k < kmax; k++)
        {
          for (int j = 0; j < jmax; j++)
            {
              for (int i = 0; i < imax; i++)
                {
                  double nonzero_indicatrice_l = next_time ? cut_cell_disc.get_interfaces().In_nonzero(1,i,j,k) : cut_cell_disc.get_interfaces().I_nonzero(1,i,j,k);
                  double nonzero_indicatrice_v = next_time ? cut_cell_disc.get_interfaces().In_nonzero(0,i,j,k) : cut_cell_disc.get_interfaces().I_nonzero(0,i,j,k);

                  int n = cut_cell_disc.get_n(i,j,k);
                  if (n < 0)
                    {
                      if (cellule_rk_restreint.pure_(i,j,k) == 0)
                        {
                          double x = F.pure_(i, j, k) * facteurF + dv.pure_(i,j,k);
                          double next_v_vol = v.pure_(i,j,k) + x * delta_t_divided_by_Fk;
                          v.pure_(i,j,k) = next_v_vol;
                        }
                      else
                        {
                          double x = dv.pure_(i,j,k);
                          double next_v_vol = v.pure_(i,j,k) + x * intermediate_dt;
                          v.pure_(i,j,k) = next_v_vol;
                        }
                    }
                  else
                    {
                      if (cellule_rk_restreint.diph_v_(n) == 0)
                        {
                          double x_v = F.diph_v_(n) * facteurF + dv.diph_v_(n);

                          double next_v_vol_v = v.diph_v_(n)*nonzero_indicatrice_v + x_v * delta_t_divided_by_Fk;
                          v.diph_v_(n) = next_v_vol_v/nonzero_indicatrice_v;
                        }
                      else
                        {
                          double x_v = dv.diph_v_(n);

                          double next_v_vol_v = v.diph_v_(n)*nonzero_indicatrice_v + x_v * intermediate_dt;
                          v.diph_v_(n) = next_v_vol_v/nonzero_indicatrice_v;
                        }

                      if (cellule_rk_restreint.diph_l_(n) == 0)
                        {
                          double x_l = F.diph_l_(n) * facteurF + dv.diph_l_(n);

                          double next_v_vol_l = v.diph_l_(n)*nonzero_indicatrice_l + x_l * delta_t_divided_by_Fk;
                          v.diph_l_(n) = next_v_vol_l/nonzero_indicatrice_l;
                        }
                      else
                        {
                          double x_l = dv.diph_l_(n);

                          double next_v_vol_l = v.diph_l_(n)*nonzero_indicatrice_l + x_l * intermediate_dt;
                          v.diph_l_(n) = next_v_vol_l/nonzero_indicatrice_l;
                        }
                    }
                }
            }
        }
      break;
    default:
      Cerr << "Error in runge_kutta_update: wrong step" << finl;
      Process::exit();
    };
}

void euler_explicit_update_cut_cell_transport(double timestep, const Cut_field_double& dv, Cut_field_double& v)
{
  const Cut_cell_FT_Disc& cut_cell_disc = v.get_cut_cell_disc();
  const double delta_t = timestep;
  const int imax = v.ni();
  const int jmax = v.nj();
  const int kmax = v.nk();
  for (int k = 0; k < kmax; k++)
    {
      for (int j = 0; j < jmax; j++)
        {
          for (int i = 0; i < imax; i++)
            {
              double old_nonzero_indicatrice_l  = cut_cell_disc.get_interfaces().I_nonzero(1,i,j,k);
              double next_nonzero_indicatrice_l = cut_cell_disc.get_interfaces().In_nonzero(1,i,j,k);
              double old_nonzero_indicatrice_v  = cut_cell_disc.get_interfaces().I_nonzero(0,i,j,k);
              double next_nonzero_indicatrice_v = cut_cell_disc.get_interfaces().In_nonzero(0,i,j,k);

              int n = cut_cell_disc.get_n(i,j,k);
              if (n < 0)
                {
                  double x = dv.pure_(i,j,k);
                  assert(old_nonzero_indicatrice_l == next_nonzero_indicatrice_l);
                  double next_v_vol = v.pure_(i,j,k) + x * delta_t;
                  v.pure_(i,j,k) = next_v_vol;
                }
              else
                {
                  double x_l = dv.diph_l_(n);
                  double x_v = dv.diph_v_(n);

                  double next_v_vol_l = v.diph_l_(n)*old_nonzero_indicatrice_l + x_l * delta_t;
                  double next_v_vol_v = v.diph_v_(n)*old_nonzero_indicatrice_v + x_v * delta_t;
                  v.diph_l_(n) = next_v_vol_l/next_nonzero_indicatrice_l;
                  v.diph_v_(n) = next_v_vol_v/next_nonzero_indicatrice_v;
                }
            }
        }
    }
}

void runge_kutta3_update_cut_cell_transport(const Cut_field_double& dv, Cut_field_double& F, Cut_field_double& v, const int step, double dt_tot, const Cut_field_int& cellule_rk_restreint)
{
  const double coeff_a[3] = { 0., -5. / 9., -153. / 128. };
  // Fk[0] = 1; Fk[i+1] = Fk[i] * a[i+1] + 1
  const double coeff_Fk[3] = { 1., 4. / 9., 15. / 32. };

  const double facteurF = coeff_a[step];
  const double intermediate_dt = compute_fractionnal_timestep_rk3(dt_tot, step);
  const double delta_t_divided_by_Fk = intermediate_dt / coeff_Fk[step];
  const int imax = v.ni();
  const int jmax = v.nj();
  const int kmax = v.nk();
  const Cut_cell_FT_Disc& cut_cell_disc = v.get_cut_cell_disc();
  switch(step)
    {
    case 0:
      // don't read initial value of F (no performance benefit because write to F causes the
      // processor to fetch the cache line, but we don't wand to use a potentially uninitialized value
      for (int k = 0; k < kmax; k++)
        {
          for (int j = 0; j < jmax; j++)
            {
              for (int i = 0; i < imax; i++)
                {
                  double old_nonzero_indicatrice_l  = cut_cell_disc.get_interfaces().I_nonzero(1,i,j,k);
                  double next_nonzero_indicatrice_l = cut_cell_disc.get_interfaces().In_nonzero(1,i,j,k);
                  double old_nonzero_indicatrice_v  = cut_cell_disc.get_interfaces().I_nonzero(0,i,j,k);
                  double next_nonzero_indicatrice_v = cut_cell_disc.get_interfaces().In_nonzero(0,i,j,k);

                  int n = cut_cell_disc.get_n(i,j,k);
                  if (n < 0)
                    {
                      double x = dv.pure_(i,j,k);
                      assert(old_nonzero_indicatrice_l == next_nonzero_indicatrice_l);
                      double next_v_vol = v.pure_(i,j,k) + x * delta_t_divided_by_Fk;
                      F.pure_(i,j,k) = x;
                      v.pure_(i,j,k) = next_v_vol;
                    }
                  else
                    {
                      double x_l = dv.diph_l_(n);
                      double x_v = dv.diph_v_(n);

                      double next_v_vol_l = v.diph_l_(n)*old_nonzero_indicatrice_l + x_l * delta_t_divided_by_Fk;
                      double next_v_vol_v = v.diph_v_(n)*old_nonzero_indicatrice_v + x_v * delta_t_divided_by_Fk;
                      F.diph_l_(n) = x_l;
                      F.diph_v_(n) = x_v;

                      v.diph_l_(n) = next_v_vol_l/next_nonzero_indicatrice_l;
                      v.diph_v_(n) = next_v_vol_v/next_nonzero_indicatrice_v;
                    }
                }
            }
        }
      break;
    case 1:
      // general case, read and write F
      for (int k = 0; k < kmax; k++)
        {
          for (int j = 0; j < jmax; j++)
            {
              for (int i = 0; i < imax; i++)
                {
                  double old_nonzero_indicatrice_l  = cut_cell_disc.get_interfaces().I_nonzero(1,i,j,k);
                  double next_nonzero_indicatrice_l = cut_cell_disc.get_interfaces().In_nonzero(1,i,j,k);
                  double old_nonzero_indicatrice_v  = cut_cell_disc.get_interfaces().I_nonzero(0,i,j,k);
                  double next_nonzero_indicatrice_v = cut_cell_disc.get_interfaces().In_nonzero(0,i,j,k);

                  int n = cut_cell_disc.get_n(i,j,k);
                  if (n < 0)
                    {
                      if (cellule_rk_restreint.pure_(i,j,k) == 0)
                        {
                          double x = F.pure_(i, j, k) * facteurF + dv.pure_(i,j,k);
                          assert(old_nonzero_indicatrice_l == next_nonzero_indicatrice_l);
                          double next_v_vol = v.pure_(i,j,k) + x * delta_t_divided_by_Fk;
                          F.pure_(i,j,k) = x;
                          v.pure_(i,j,k) = next_v_vol;
                        }
                      else
                        {
                          double x = dv.pure_(i,j,k);
                          assert(old_nonzero_indicatrice_l == next_nonzero_indicatrice_l);
                          double next_v_vol = v.pure_(i,j,k) + x * intermediate_dt;
                          F.pure_(i,j,k) = x;
                          v.pure_(i,j,k) = next_v_vol;
                        }
                    }
                  else
                    {
                      if (cellule_rk_restreint.diph_v_(n) == 0)
                        {
                          double x_v = F.diph_v_(n) * facteurF + dv.diph_v_(n);

                          double next_v_vol_v = v.diph_v_(n)*old_nonzero_indicatrice_v + x_v * delta_t_divided_by_Fk;
                          F.diph_v_(n) = x_v;

                          v.diph_v_(n) = next_v_vol_v/next_nonzero_indicatrice_v;
                        }
                      else
                        {
                          double x_v = dv.diph_v_(n);

                          double next_v_vol_v = v.diph_v_(n)*old_nonzero_indicatrice_v + x_v * intermediate_dt;
                          F.diph_v_(n) = x_v;

                          v.diph_v_(n) = next_v_vol_v/next_nonzero_indicatrice_v;
                        }

                      if (cellule_rk_restreint.diph_l_(n) == 0)
                        {
                          double x_l = F.diph_l_(n) * facteurF + dv.diph_l_(n);

                          double next_v_vol_l = v.diph_l_(n)*old_nonzero_indicatrice_l + x_l * delta_t_divided_by_Fk;
                          F.diph_l_(n) = x_l;

                          v.diph_l_(n) = next_v_vol_l/next_nonzero_indicatrice_l;
                        }
                      else
                        {
                          double x_l = dv.diph_l_(n);

                          double next_v_vol_l = v.diph_l_(n)*old_nonzero_indicatrice_l + x_l * intermediate_dt;
                          F.diph_l_(n) = x_l;

                          v.diph_l_(n) = next_v_vol_l/next_nonzero_indicatrice_l;
                        }
                    }
                }
            }
        }
      break;
    case 2:
      // do not write F
      for (int k = 0; k < kmax; k++)
        {
          for (int j = 0; j < jmax; j++)
            {
              for (int i = 0; i < imax; i++)
                {
                  double old_nonzero_indicatrice_l  = cut_cell_disc.get_interfaces().I_nonzero(1,i,j,k);
                  double next_nonzero_indicatrice_l = cut_cell_disc.get_interfaces().In_nonzero(1,i,j,k);
                  double old_nonzero_indicatrice_v  = cut_cell_disc.get_interfaces().I_nonzero(0,i,j,k);
                  double next_nonzero_indicatrice_v = cut_cell_disc.get_interfaces().In_nonzero(0,i,j,k);

                  int n = cut_cell_disc.get_n(i,j,k);
                  if (n < 0)
                    {
                      if (cellule_rk_restreint.pure_(i,j,k) == 0)
                        {
                          double x = F.pure_(i, j, k) * facteurF + dv.pure_(i,j,k);
                          assert(old_nonzero_indicatrice_l == next_nonzero_indicatrice_l);
                          double next_v_vol = v.pure_(i,j,k) + x * delta_t_divided_by_Fk;
                          v.pure_(i,j,k) = next_v_vol;
                        }
                      else
                        {
                          double x = dv.pure_(i,j,k);
                          assert(old_nonzero_indicatrice_l == next_nonzero_indicatrice_l);
                          double next_v_vol = v.pure_(i,j,k) + x * intermediate_dt;
                          v.pure_(i,j,k) = next_v_vol;
                        }
                    }
                  else
                    {
                      if (cellule_rk_restreint.diph_v_(n) == 0)
                        {
                          double x_v = F.diph_v_(n) * facteurF + dv.diph_v_(n);

                          double next_v_vol_v = v.diph_v_(n)*old_nonzero_indicatrice_v + x_v * delta_t_divided_by_Fk;

                          v.diph_v_(n) = next_v_vol_v/next_nonzero_indicatrice_v;
                        }
                      else
                        {
                          double x_v = dv.diph_v_(n);

                          double next_v_vol_v = v.diph_v_(n)*old_nonzero_indicatrice_v + x_v * intermediate_dt;

                          v.diph_v_(n) = next_v_vol_v/next_nonzero_indicatrice_v;
                        }

                      if (cellule_rk_restreint.diph_l_(n) == 0)
                        {
                          double x_l = F.diph_l_(n) * facteurF + dv.diph_l_(n);

                          double next_v_vol_l = v.diph_l_(n)*old_nonzero_indicatrice_l + x_l * delta_t_divided_by_Fk;

                          v.diph_l_(n) = next_v_vol_l/next_nonzero_indicatrice_l;
                        }
                      else
                        {
                          double x_l = dv.diph_l_(n);

                          double next_v_vol_l = v.diph_l_(n)*old_nonzero_indicatrice_l + x_l * intermediate_dt;

                          v.diph_l_(n) = next_v_vol_l/next_nonzero_indicatrice_l;
                        }
                    }
                }
            }
        }
      break;
    default:
      Cerr << "Error in runge_kutta_update: wrong step" << finl;
      Process::exit();
    };
}

void cut_cell_switch_field_time(Cut_field_double& v)
{
  const Cut_cell_FT_Disc& cut_cell_disc = v.get_cut_cell_disc();
  const int imax = v.ni();
  const int jmax = v.nj();
  const int kmax = v.nk();
  for (int k = 0; k < kmax; k++)
    {
      for (int j = 0; j < jmax; j++)
        {
          for (int i = 0; i < imax; i++)
            {
              double old_nonzero_indicatrice_l  = cut_cell_disc.get_interfaces().I_nonzero(1,i,j,k);
              double next_nonzero_indicatrice_l = cut_cell_disc.get_interfaces().In_nonzero(1,i,j,k);
              double old_nonzero_indicatrice_v  = cut_cell_disc.get_interfaces().I_nonzero(0,i,j,k);
              double next_nonzero_indicatrice_v = cut_cell_disc.get_interfaces().In_nonzero(0,i,j,k);

              int n = cut_cell_disc.get_n(i,j,k);
              if (n < 0)
                {
                  assert(old_nonzero_indicatrice_l == next_nonzero_indicatrice_l);
                  double next_v_vol = v.pure_(i,j,k);
                  v.pure_(i,j,k) = next_v_vol;
                }
              else
                {
                  double next_v_vol_l = v.diph_l_(n)*old_nonzero_indicatrice_l;
                  double next_v_vol_v = v.diph_v_(n)*old_nonzero_indicatrice_v;

                  v.diph_l_(n) = next_v_vol_l/next_nonzero_indicatrice_l;
                  v.diph_v_(n) = next_v_vol_v/next_nonzero_indicatrice_v;
                }
            }
        }
    }
}

void runge_kutta3_update_surfacic_fluxes(Cut_field_double& dv, Cut_field_double& F, const int step, const int k_layer, const int dir, double dt_tot, const Cut_field_int& cellule_rk_restreint)
{
  const double coeff_a[3] = { 0., -5. / 9., -153. / 128. };
  // Fk[0] = 1; Fk[i+1] = Fk[i] * a[i+1] + 1
  const double coeff_Fk[3] = { 1., 4. / 9., 15. / 32. };

  int di = (dir == 0) ? -1 : 0;
  int dj = (dir == 1) ? -1 : 0;
  int dk = (dir == 2) ? -1 : 0;

  const double facteurF = coeff_a[step];
  const double one_divided_by_Fk = 1. / coeff_Fk[step];
  const int imax = dv.ni();
  const int jmax = dv.nj();
  const int ghost = dv.ghost();
  const Cut_cell_FT_Disc& cut_cell_disc = dv.get_cut_cell_disc();
  switch(step)
    {
    case 0:
      // don't read initial value of F (no performance benefit because write to F causes the
      // processor to fetch the cache line, but we don't wand to use a potentially uninitialized value
      for (int j = -ghost; j < jmax+ghost; j++)
        {
          for (int i = -ghost; i < imax+ghost; i++)
            {
              int n = cut_cell_disc.get_n(i,j,k_layer);
              if (n < 0)
                {
                  double x = dv.pure_(i,j,k_layer);
                  dv.pure_(i,j,k_layer) = x * one_divided_by_Fk;
                  F.pure_(i,j,k_layer) = x;
                }
              else
                {
                  double x_l = dv.diph_l_(n);
                  double x_v = dv.diph_v_(n);

                  dv.diph_l_(n) = x_l * one_divided_by_Fk;
                  dv.diph_v_(n) = x_v * one_divided_by_Fk;
                  F.diph_l_(n) = x_l;
                  F.diph_v_(n) = x_v;

                  int n_decale = cut_cell_disc.get_n(i+di,j+dj,k_layer+dk);
                  if (n_decale < 0)
                    {
                      int phase = IJK_Interfaces::convert_indicatrice_to_phase(cut_cell_disc.indic_pure(i+di,j+dj,k_layer+dk));
                      dv.pure_(i,j,k_layer) = (phase == 0) ? dv.diph_v_(n) : dv.diph_l_(n);
                      F.pure_(i,j,k_layer) = (phase == 0) ? F.diph_v_(n) : F.diph_l_(n);
                    }
                }
            }
        }
      break;
    case 1:
      // general case, read and write F
      for (int j = -ghost; j < jmax+ghost; j++)
        {
          for (int i = -ghost; i < imax+ghost; i++)
            {
              int n = cut_cell_disc.get_n(i,j,k_layer);
              if (n < 0)
                {
                  int phase = IJK_Interfaces::convert_indicatrice_to_phase(cut_cell_disc.indic_pure(i,j,k_layer));
                  if ((cellule_rk_restreint.pure_(i,j,k_layer) == 0) && (cellule_rk_restreint.from_ijk_and_phase(i+di,j+dj,k_layer+dk,phase) == 0))
                    {
                      double x = F.pure_(i, j, k_layer) * facteurF + dv.pure_(i,j,k_layer);
                      dv.pure_(i,j,k_layer) = x * one_divided_by_Fk;
                      F.pure_(i,j,k_layer) = x;
                    }
                  else
                    {
                      double x = dv.pure_(i,j,k_layer);
                      dv.pure_(i,j,k_layer) = x;
                      F.pure_(i,j,k_layer) = x;
                    }
                }
              else
                {
                  if ((cellule_rk_restreint.diph_v_(n) == 0) && (cellule_rk_restreint.from_ijk_and_phase(i+di,j+dj,k_layer+dk,0) == 0))
                    {
                      double x_v = F.diph_v_(n) * facteurF + dv.diph_v_(n);

                      dv.diph_v_(n) = x_v * one_divided_by_Fk;
                      F.diph_v_(n) = x_v;
                    }
                  else
                    {
                      double x_v = dv.diph_v_(n);

                      dv.diph_v_(n) = x_v;
                      F.diph_v_(n) = x_v;
                    }

                  if ((cellule_rk_restreint.diph_l_(n) == 0) && (cellule_rk_restreint.from_ijk_and_phase(i+di,j+dj,k_layer+dk,1) == 0))
                    {
                      double x_l = F.diph_l_(n) * facteurF + dv.diph_l_(n);

                      dv.diph_l_(n) = x_l * one_divided_by_Fk;
                      F.diph_l_(n) = x_l;
                    }
                  else
                    {
                      double x_l = dv.diph_l_(n);

                      dv.diph_l_(n) = x_l;
                      F.diph_l_(n) = x_l;
                    }

                  int n_decale = cut_cell_disc.get_n(i+di,j+dj,k_layer+dk);
                  if (n_decale < 0)
                    {
                      int phase = IJK_Interfaces::convert_indicatrice_to_phase(cut_cell_disc.indic_pure(i+di,j+dj,k_layer+dk));
                      dv.pure_(i,j,k_layer) = (phase == 0) ? dv.diph_v_(n) : dv.diph_l_(n);
                      F.pure_(i,j,k_layer) = (phase == 0) ? F.diph_v_(n) : F.diph_l_(n);
                    }
                }
            }
        }
      break;
    case 2:
      // do not write F
      for (int j = -ghost; j < jmax+ghost; j++)
        {
          for (int i = -ghost; i < imax+ghost; i++)
            {
              int n = cut_cell_disc.get_n(i,j,k_layer);
              if (n < 0)
                {
                  int phase = IJK_Interfaces::convert_indicatrice_to_phase(cut_cell_disc.indic_pure(i,j,k_layer));
                  if ((cellule_rk_restreint.pure_(i,j,k_layer) == 0) && (cellule_rk_restreint.from_ijk_and_phase(i+di,j+dj,k_layer+dk,phase) == 0))
                    {
                      double x = F.pure_(i, j, k_layer) * facteurF + dv.pure_(i,j,k_layer);
                      dv.pure_(i,j,k_layer) = x * one_divided_by_Fk;
                    }
                  else
                    {
                      double x = dv.pure_(i,j,k_layer);
                      dv.pure_(i,j,k_layer) = x;
                    }
                }
              else
                {
                  if ((cellule_rk_restreint.diph_v_(n) == 0) && (cellule_rk_restreint.from_ijk_and_phase(i+di,j+dj,k_layer+dk,0) == 0))
                    {
                      double x_v = F.diph_v_(n) * facteurF + dv.diph_v_(n);

                      dv.diph_v_(n) = x_v * one_divided_by_Fk;
                    }
                  else
                    {
                      double x_v = dv.diph_v_(n);

                      dv.diph_v_(n) = x_v;
                    }

                  if ((cellule_rk_restreint.diph_l_(n) == 0) && (cellule_rk_restreint.from_ijk_and_phase(i+di,j+dj,k_layer+dk,1) == 0))
                    {
                      double x_l = F.diph_l_(n) * facteurF + dv.diph_l_(n);

                      dv.diph_l_(n) = x_l * one_divided_by_Fk;
                    }
                  else
                    {
                      double x_l = dv.diph_l_(n);

                      dv.diph_l_(n) = x_l;
                    }

                  int n_decale = cut_cell_disc.get_n(i+di,j+dj,k_layer+dk);
                  if (n_decale < 0)
                    {
                      int phase = IJK_Interfaces::convert_indicatrice_to_phase(cut_cell_disc.indic_pure(i+di,j+dj,k_layer+dk));
                      dv.pure_(i,j,k_layer) = (phase == 0) ? dv.diph_v_(n) : dv.diph_l_(n);
                      F.pure_(i,j,k_layer) = (phase == 0) ? F.diph_v_(n) : F.diph_l_(n);
                    }
                }
            }
        }
      break;
    default:
      Cerr << "Error in runge_kutta_update: wrong step" << finl;
      Process::exit();
    };
}

// Cette fonction sert ajouter la pente liee a la correction a la pente liee au schema principal
// Le facteur vol_over_dt_surface existe car le flux calcule dans les routines de correction n'a
// pas la bonne unite.
void add_flux_times_vol_over_dt_surface(double fractional_timestep, const Cut_field_vector3_double& cut_field_current_fluxes, Cut_field_vector3_double& cut_field_RK3_F_fluxes)
{
  const IJK_Grid_Geometry& geom = cut_field_current_fluxes[0].get_cut_cell_disc().get_splitting().get_grid_geometry();
  assert(cut_field_current_fluxes[0].get_cut_cell_disc().get_splitting().get_grid_geometry().is_uniform(0));
  assert(cut_field_current_fluxes[0].get_cut_cell_disc().get_splitting().get_grid_geometry().is_uniform(1));
  assert(cut_field_current_fluxes[0].get_cut_cell_disc().get_splitting().get_grid_geometry().is_uniform(2));
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  double vol = dx*dy*dz;

  for (int dir = 0; dir < 3; dir++)
    {
      const int ni = cut_field_current_fluxes[dir].ni();
      const int nj = cut_field_current_fluxes[dir].nj();
      const int nk = cut_field_current_fluxes[dir].nk();
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  double vol_over_dt = vol/fractional_timestep;

                  cut_field_RK3_F_fluxes[dir].pure_(i,j,k) += cut_field_current_fluxes[dir].pure_(i,j,k)*vol_over_dt;
                }
            }
        }

      const Cut_cell_FT_Disc& cut_cell_disc = cut_field_current_fluxes[dir].get_cut_cell_disc();
      for (int n = 0; n < cut_cell_disc.get_n_tot(); n++)
        {
          double vol_over_dt = vol/fractional_timestep;

          const DoubleTabFT_cut_cell_vector3& indicatrice_surfacique = cut_cell_disc.get_interfaces().get_indicatrice_surfacique_efficace_face();
          double indicatrice_surface = indicatrice_surfacique(n,dir);

          cut_field_RK3_F_fluxes[dir].diph_l_(n) += (indicatrice_surface == 0.)       ? 0. : cut_field_current_fluxes[dir].diph_l_(n)*vol_over_dt/indicatrice_surface;
          cut_field_RK3_F_fluxes[dir].diph_v_(n) += ((1 - indicatrice_surface) == 0.) ? 0. : cut_field_current_fluxes[dir].diph_v_(n)*vol_over_dt/(1 - indicatrice_surface);
        }
    }
}

void set_rk_restreint(int rk_step, int rk_restriction_leniency, const Cut_cell_FT_Disc& cut_cell_disc, Cut_field_int& cellule_rk_restreint)
{
  if (rk_step == 0)
    {
      cellule_rk_restreint.set_to_uniform_value(0.);
    }

  if (rk_restriction_leniency == 10)
    {
      cellule_rk_restreint.set_to_uniform_value(1.);
    }

  if (rk_restriction_leniency == 11)
    {
      for (int n = 0; n < cut_cell_disc.get_n_tot(); n++)
        {
          Int3 ijk = cut_cell_disc.get_ijk(n);
          int i = ijk[0];
          int j = ijk[1];
          int k = ijk[2];

          cellule_rk_restreint.pure_(i,j,k) = 1;
          cellule_rk_restreint.diph_v_(n) = 1;
          cellule_rk_restreint.diph_l_(n) = 1;
        }
    }

  // Boucle sur les cellules qui disparaissent lors de ce sous pas de temps
  {
    int statut_diphasique = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::MOURRANT);
    int index_min = cut_cell_disc.get_statut_diphasique_value_index(statut_diphasique);
    int index_max = cut_cell_disc.get_statut_diphasique_value_index(statut_diphasique+1);
    for (int index = index_min; index < index_max; index++)
      {
        int n = cut_cell_disc.get_n_from_statut_diphasique_index(index);

        Int3 ijk = cut_cell_disc.get_ijk(n);
        int i = ijk[0];
        int j = ijk[1];
        int k = ijk[2];

        if ((rk_restriction_leniency == 4) || (rk_restriction_leniency == 5))
          continue;

        if ((rk_restriction_leniency == 1) || (rk_restriction_leniency == 3))
          {
            double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);
            int phase = 1 - IJK_Interfaces::convert_indicatrice_to_phase(next_indicatrice); // phase de la cellule mourrante

            if (phase == 0)
              {
                cellule_rk_restreint.diph_v_(n) = 1;
              }
            else
              {
                cellule_rk_restreint.diph_l_(n) = 1;
              }
          }
        else if ((rk_restriction_leniency == 0) || (rk_restriction_leniency == 2))
          {
            cellule_rk_restreint.diph_v_(n) = 1;
            cellule_rk_restreint.diph_l_(n) = 1;
          }
      }
  }

  // Boucle sur les cellules qui apparaissent ou bien qui sont petites lors de ce sous pas de temps
  {
    int statut_diphasique_naissant = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::NAISSANT);
    int statut_diphasique_petit = static_cast<int>(Cut_cell_FT_Disc::STATUT_DIPHASIQUE::DESEQUILIBRE_FINAL);
    assert(statut_diphasique_petit == statut_diphasique_naissant + 1);
    int index_min = cut_cell_disc.get_statut_diphasique_value_index(statut_diphasique_naissant);
    int index_max = cut_cell_disc.get_statut_diphasique_value_index(statut_diphasique_petit+1);
    for (int index = index_min; index < index_max; index++)
      {
        int n = cut_cell_disc.get_n_from_statut_diphasique_index(index);

        Int3 ijk = cut_cell_disc.get_ijk(n);
        int i = ijk[0];
        int j = ijk[1];
        int k = ijk[2];

        if (rk_step == 0 && ((rk_restriction_leniency == 5) || (rk_restriction_leniency == 4) || (rk_restriction_leniency == 3) || (rk_restriction_leniency == 2)))
          continue;

        if ((rk_restriction_leniency == 1) || (rk_restriction_leniency == 3) || (rk_restriction_leniency == 5))
          {
            double old_indicatrice = cut_cell_disc.get_interfaces().I(i,j,k);
            double next_indicatrice = cut_cell_disc.get_interfaces().In(i,j,k);
            int est_naissant = cut_cell_disc.get_interfaces().est_pure(old_indicatrice);
            int phase = est_naissant ? 1 - IJK_Interfaces::convert_indicatrice_to_phase(old_indicatrice) : ((cut_cell_disc.get_interfaces().below_small_threshold(next_indicatrice)) ? 1 : 0); // phase de la cellule petite ou naissante

            if (phase == 0)
              {
                cellule_rk_restreint.diph_v_(n) = 1;
              }
            else
              {
                cellule_rk_restreint.diph_l_(n) = 1;
              }
          }
        else if ((rk_restriction_leniency == 0) || (rk_restriction_leniency == 2) || (rk_restriction_leniency == 4))
          {
            cellule_rk_restreint.diph_v_(n) = 1;
            cellule_rk_restreint.diph_l_(n) = 1;
          }
      }
  }
}

