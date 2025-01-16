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
/////////////////////////////////////////////////////////////////////////////
//
// File      : IJK_Ghost_Fluid_tools.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Ghost_Fluid_tools.h>
#include <IJK_Field_vector.h>
#include <Probleme_base.h>
#include <DebogIJK.h>
#include <stat_counters.h>
#include <Operateur_IJK_elem_diff.h>

static int decoder_numero_bulle(const int code)
{
  const int num_bulle = code >> 6;
  return num_bulle;
}

static void extrapolate_with_elem_faces_connectivity(const Domaine_VF& domaine_vf,
                                                     const IJK_Field_double& distance,
                                                     IJK_Field_double& field,
                                                     const int stencil_width)
{
  const double invalid_test = INVALID_TEST;
  const IntTab& elem_faces = domaine_vf.elem_faces();
  const IntTab& faces_elem = domaine_vf.face_voisins();
  const int nb_faces_elem = elem_faces.dimension(1);
  const int nb_elem = elem_faces.dimension(0);
  IJK_Field_double field_old;
  // Use to locate the initial non-zero values
  IJK_Field_double field_ini(field);
  const IJK_Splitting& splitting_distance = distance.get_splitting();
  /*
   * n_iterations = stencil_width is the minimum to get a propagation of information from the interface to the border
   * of the extrapolation. But doing more will lead to smoother values... And it probably costs close to nothing
   */
  const double n_iterations = 5 * stencil_width;
  for (int iteration = 0; iteration < n_iterations; iteration++)
    {
      // Necessary !!!
      // Copy the old field value as we do not want to use the current iteration values.
      field_old = field;
      // La valeur sur un element est la moyenne des valeurs sur les elements voisins
      for (int i_elem = 0; i_elem < nb_elem; i_elem++)
        {
          // Do not touch field in interfacial cells.
          // Iterate on other values.
          const Int3 num_elem_ijk = splitting_distance.convert_packed_to_ijk_cell(i_elem);
          const double d = distance(num_elem_ijk[DIRECTION_I],
                                    num_elem_ijk[DIRECTION_J],
                                    num_elem_ijk[DIRECTION_K]);
          //          const double interfacial_area_elem = interfacial_area(num_elem_ijk[DIRECTION_I],
          //                                                                num_elem_ijk[DIRECTION_J],
          //                                                                num_elem_ijk[DIRECTION_K]);
          // Need a value of distance but don't overwrite the first calculated value
          // if ((d > invalid_test) && (interfacial_area_elem < invalid_test))
          const double field_ini_val = field_ini(num_elem_ijk[DIRECTION_I],
                                                 num_elem_ijk[DIRECTION_J],
                                                 num_elem_ijk[DIRECTION_K]);
          if ((d > invalid_test) && (field_ini_val == 0))
            {
              double sum_field = 0.;
              double coeff = 0.;
              for (int i_face = 0; i_face < nb_faces_elem; i_face++)
                {
                  const int face = elem_faces(i_elem, i_face);
                  const int neighbour = faces_elem(face, 0) + faces_elem(face, 1) - i_elem;
                  if (neighbour >= 0)
                    {
                      // Not a boundary...
                      const Int3 num_elem_neighbour_ijk = splitting_distance.convert_packed_to_ijk_cell(neighbour);
                      double field_neighbour = field_old(num_elem_neighbour_ijk[DIRECTION_I],
                                                         num_elem_neighbour_ijk[DIRECTION_J],
                                                         num_elem_neighbour_ijk[DIRECTION_K]);
                      const double distance_neighbour = distance(num_elem_neighbour_ijk[DIRECTION_I],
                                                                 num_elem_neighbour_ijk[DIRECTION_J],
                                                                 num_elem_neighbour_ijk[DIRECTION_K]);
                      // Don't use zero values
                      // if ((distance_neighbour > invalid_test) && (field_neighbour != 0))
                      // Use zero_values grad_T_ decreasing with distance
                      if (distance_neighbour > invalid_test)
                        {
                          // Give more weight in the smoothing to values closer to the interface:
                          if (fabs(distance_neighbour) < INVALID_TEST)
                            {
                              Cerr << "Distance is very much at zero whereas interfacial_area is zero too... Pathological case to be looked into closely. " << finl;
                              Cerr << "Contact TRUST support." << finl;
                              Process::exit();
                            }
                          /*
                           * TODO: Check the difference between extrapoler_champ_elem
                           * and extrapolate from Convection_Diffusion_Temperature_FT_Disc.cpp
                           */
//                          const double inv_distance_squared = 1./ (distance_neighbour * distance_neighbour);
//                          sum_field += field_neighbour * inv_distance_squared;
//                          coeff += inv_distance_squared;
                          sum_field += field_neighbour;
                          coeff++;
                        }
                    }
                }
              if (coeff > 0.)
                {
                  field(num_elem_ijk[DIRECTION_I],
                        num_elem_ijk[DIRECTION_J],
                        num_elem_ijk[DIRECTION_K]) = sum_field / coeff;
                }
            }
        }
      field.echange_espace_virtuel(field.ghost());
    }
}

static void extrapolate_with_ijk_indices(const IJK_Field_double& distance,
                                         const IJK_Field_double& indicator,
                                         IJK_Field_double& field,
                                         const int& stencil_width,
                                         const int& recompute_field_ini,
                                         const int& zero_neighbour_value_mean,
                                         const int& vapour_mixed_only,
                                         const int& smooth_factor)
{
  static const Stat_Counter_Id stat_counter = statistiques().new_counter(3, "GFM - Extrapolate gfm temperature values");
  statistiques().begin_count(stat_counter);

  int neighbours_i[6] = NEIGHBOURS_I;
  int neighbours_j[6] = NEIGHBOURS_J;
  int neighbours_k[6] = NEIGHBOURS_K;
  const double invalid_test = INVALID_TEST;
  const double n_iterations = zero_neighbour_value_mean ? smooth_factor * stencil_width : stencil_width;
  const int ni = field.ni();
  const int nj = field.nj();
  const int nk = field.nk();
  IJK_Field_double field_old;
  IJK_Field_double field_ini = field;
  for (int iteration = 0; iteration < n_iterations; iteration++)
    {
      field_old = field;
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double indic = indicator(i,j,k);
              const int vapour_mixed = !vapour_mixed_only ? 1 : (indic < LIQUID_INDICATOR_TEST);
              if (vapour_mixed)
                {
                  const double d = distance(i,j,k);
                  // const double field_ini_val = field_ini(i,j,k);
                  const double field_ini_val = (!zero_neighbour_value_mean && !recompute_field_ini) ? field_old(i,j,k) : field_ini(i,j,k);
                  if ((d > invalid_test) && (field_ini_val == 0))
                    {
                      double sum_field = 0.;
                      double coeff = 0.;
                      for (int l=0; l<6; l++)
                        {
                          const int ii = neighbours_i[l];
                          const int jj = neighbours_j[l];
                          const int kk = neighbours_k[l];
                          const double distance_neighbour = distance(i+ii,j+jj,k+kk);
                          const double field_neighbour = field_old(i+ii,j+jj,k+kk);
                          const int neighbour_condition = zero_neighbour_value_mean ? 1 : (field_neighbour != 0);
                          // if (distance_neighbour > invalid_test)
                          // Don't use zero values
                          // if ((distance_neighbour > invalid_test) && field_neighbour != 0))
                          // Use zero_values grad_T_ decreasing with distance
                          if (distance_neighbour > invalid_test && neighbour_condition)
                            {
                              sum_field += field_neighbour;
                              coeff++;
                            }
                        }
                      if (coeff > 0.)
                        field(i,j,k) = sum_field / coeff;
                    }
                }
            }
      field.echange_espace_virtuel(field.ghost());
    }
  statistiques().end_count(stat_counter);
}

/*
 * TODO: Store the facets' barycentres
 */
void compute_eulerian_normal_distance_facet_barycentre_field(const IJK_Interfaces& interfaces, //  ref_problem_ft_disc,
                                                             IJK_Field_double& distance_field,
                                                             IJK_Field_vector3_double& normal_vect,
                                                             IJK_Field_vector3_double& facets_barycentre,
                                                             IJK_Field_vector3_double& tmp_old_vector_val,
                                                             IJK_Field_vector3_double& tmp_new_vector_val,
                                                             IJK_Field_double& tmp_old_val,
                                                             IJK_Field_double& tmp_new_val,
                                                             IJK_Field_int& tmp_interf_cells,
                                                             IJK_Field_int& tmp_propagated_cells,
                                                             FixedVector<ArrOfInt,3>& interf_cells_indices,
                                                             FixedVector<ArrOfInt,3>& gfm_first_cells_indices_,
                                                             FixedVector<ArrOfInt,3>& propagated_cells_indices,
                                                             const int& n_iter,
                                                             const int& avoid_gfm_parallel_calls)
{
  /*
   * Compute the normal distance to the interface
   */
  const bool use_ijk = true;
  static const Stat_Counter_Id stat_counter = statistiques().new_counter(2, "GFM - Compute Eulerian normal distance field");
  statistiques().begin_count(stat_counter);

  static const double invalid_distance_value = INVALID_TEST;
  const int dim = 3; // in IJK

  // Vertex coordinates of the eulerian domain
  const Domaine_dis_base& mon_dom_dis = interfaces.get_domaine_dis();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, mon_dom_dis);
  const DoubleTab& centre_element = domaine_vf.xp();

  // Approximation of the normal distance
  distance_field.data() = invalid_distance_value * 1.1;
  for (int dir = 0; dir < dim; dir++)
    normal_vect[dir].data() = 0.;

  const int nb_elem = mon_dom_dis.domaine().nb_elem();

  const Maillage_FT_IJK& maillage = interfaces.maillage_ft_ijk();

  // Same splitting for the normal vector field
  const IJK_Splitting& splitting_distance = distance_field.get_splitting();
  const IJK_Grid_Geometry& geom = splitting_distance.get_grid_geometry();

  for (int l=0; l<dim; l++)
    {
      interf_cells_indices[l].reset();
      gfm_first_cells_indices_[l].reset();
      propagated_cells_indices[l].reset();
    }
  tmp_interf_cells.data() = 0;
  tmp_propagated_cells.data() = 0;

  /*
   * M.G: Copy of B.M from Transport_Interfaces_FT_Disc::calculer_distance_interface
   * Distance calculation for the thickness 0 (vertices of the elements crossed by
   * the interface). For each element, we calculate the plane intersecting the
   * barycentre of the facets portions. The normal vector to the plane corresponds to
   * the average normal vector weighting by the surface portions. The distance
   * interface/element is the distance between this plane and the centre of the element.
   */
  {
    const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
    const ArrOfInt& index_elem = intersections.index_elem();
    const DoubleTab& normale_facettes = maillage.get_update_normale_facettes();
    const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
    const IntTab& facettes = maillage.facettes();
    const DoubleTab& sommets = maillage.sommets();
    // Loop on the elements

    for (int elem = 0; elem < nb_elem; elem++)
      {
        int index = index_elem[elem];
        // Moyenne ponderee des normales aux facettes qui traversent l'element
        double normale[3] = {0., 0., 0.};
        // Barycentre of the facets/element intersections
        double centre[3] = {0., 0., 0.};
        // Sum of the weights
        double surface_tot = 0.;
        // Loop on the facets which cross the element
        while (index >= 0)
          {
            const Intersections_Elem_Facettes_Data& data =
              intersections.data_intersection(index);

            const int num_facette = data.numero_facette_;
#ifdef AVEC_BUG_SURFACES
            const double surface = data.surface_intersection_;
#else
            const double surface = data.fraction_surface_intersection_ * surface_facettes[num_facette];
#endif
            surface_tot += surface;
            for (int i = 0; i < dim; i++)
              {
                normale[i] += surface * normale_facettes(num_facette, i);
                // Barycentre calculation of the facets/element intersections
                double g_i = 0.; // i-component of the barycentre coordinate
                for (int j = 0; j < dim; j++)
                  {
                    const int som   = facettes(num_facette, j);
                    const double coord = sommets(som, i);
                    const double coeff = data.barycentre_[j];
                    g_i += coord * coeff;
                  }
                centre[i] += surface * g_i;
              }
            index = data.index_facette_suivante_;
          }
        const Int3 num_elem_ijk = splitting_distance.convert_packed_to_ijk_cell(elem);
        if (surface_tot > 0.)
          {
            /*
             * The stored vector is not normed : norm ~ surface portion
             * centre = sum(centre[facet] * surface) / sum(surface)
             */
            const double inverse_surface_tot = 1. / surface_tot;
            double norme = 0.;
            int j;
            for (j = 0; j < dim; j++)
              {
                norme += normale[j] * normale[j];
                normal_vect[j](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) = normale[j];
                centre[j] *= inverse_surface_tot;
              }

            if (norme > 0)
              {
                double i_norme = 1./sqrt(norme);
                double distance = 0.;
                for (j = 0; j < dim; j++)
                  {
                    double n_j = normale[j] * i_norme; // normal vector normed
                    distance += (centre_element(elem, j) - centre[j]) * n_j;
                  }
                distance_field(num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) = distance;
              }

            if (avoid_gfm_parallel_calls)
              tmp_interf_cells(num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) = 1;

          }
        for (int j = 0; j < dim; j++)
          facets_barycentre[j](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) = centre[j];
      }
    distance_field.echange_espace_virtuel(distance_field.ghost());
    normal_vect.echange_espace_virtuel();
    tmp_interf_cells.echange_espace_virtuel(tmp_interf_cells.ghost());
    for (int dir = 0; dir < dim; dir++)
      DebogIJK::verifier("IJK_Ghost_Fluid_tools::compute_eulerian_normal_distance_field", normal_vect[dir]);
    DebogIJK::verifier("IJK_Ghost_Fluid_tools::compute_eulerian_normal_distance_field", distance_field);
  }

  if (avoid_gfm_parallel_calls)
    {
      const int ni = tmp_interf_cells.ni();
      const int nj = tmp_interf_cells.nj();
      const int nk = tmp_interf_cells.nk();
      const int nb_ghost = tmp_interf_cells.ghost();
      for (int k = - nb_ghost; k < nk + nb_ghost; k++)
        for (int j = - nb_ghost; j < nj + nb_ghost; j++)
          for (int i = - nb_ghost; i < ni + nb_ghost; i++)
            if (tmp_interf_cells(i,j,k))
              for (int l=0; l<dim; l++)
                {
                  const int index = select_dir(l, i, j, k);
                  interf_cells_indices[l].append_array(index);
                }
    }

  //  IJK_Field_vector3_double terme_src(normal_vect);
  //  IJK_Field_vector3_double tmp(normal_vect);
  IJK_Field_vector3_double& terme_src = tmp_old_vector_val;
  IJK_Field_vector3_double& tmp = tmp_new_vector_val;
  for (int l=0; l<dim; l++)
    {
      terme_src[l].data() = normal_vect[l].data();
      tmp[l].data() = normal_vect[l].data();
    }

  const IntTab& face_voisins = domaine_vf.face_voisins();
  const IntTab& elem_faces   = domaine_vf.elem_faces();
  const int nb_elem_voisins = elem_faces.line_size();

  // Normal vector calculation at the element location:
  /*
   * TODO: Check the how fast it is compared to using elem_faces matrix
   */
  tmp_propagated_cells.data() = tmp_interf_cells.data();
  propagated_cells_indices = interf_cells_indices;
  FixedVector<ArrOfInt,3> propagated_cells_indices_tmp;
  for (int l=0; l<dim; l++)
    propagated_cells_indices_tmp[l].resize_array(1);

  int iteration;
  const int n_iter_tmp = avoid_gfm_parallel_calls ? n_iter + 1: n_iter;

  if (use_ijk)
    {
      int neighbours_i[6] = NEIGHBOURS_I;
      int neighbours_j[6] = NEIGHBOURS_J;
      int neighbours_k[6] = NEIGHBOURS_K;
      const int ni = normal_vect[0].ni();
      const int nj = normal_vect[0].nj();
      const int nk = normal_vect[0].nk();
      const int nghost = normal_vect[0].ghost();
      int m;
      for (iteration = 0; iteration < n_iter_tmp; iteration++)
        {
          /*
           * Smoothing iterator, in theory:
           * normal = normal + (source_term - laplacian(normal)) * factor
           * converge towards laplacian(normal) = source_term
           * in practice:
           * normal = average(normal on neighbours) + source_term
           */

          const double un_sur_ncontrib = 1. / (1. + nb_elem_voisins);
          if (avoid_gfm_parallel_calls)
            {
              propagated_cells_indices_tmp = propagated_cells_indices;

              for (int ielem = 0; ielem < propagated_cells_indices_tmp[0].size_array(); ielem++)
                {
                  const int i = propagated_cells_indices_tmp[DIRECTION_I](ielem);
                  const int j = propagated_cells_indices_tmp[DIRECTION_J](ielem);
                  const int k = propagated_cells_indices_tmp[DIRECTION_K](ielem);

                  // Averaging the normal vector on the neighbours
                  double n[3] = {0., 0., 0.};
                  for (m = 0; m < dim; m++)
                    n[m] = normal_vect[m](i,j,k);
                  for (int l=0; l<6; l++)
                    {
                      const int ii = neighbours_i[l];
                      const int jj = neighbours_j[l];
                      const int kk = neighbours_k[l];

                      const int is_outside_proc = (i + ii < -nghost || j + jj < -nghost || k + kk < -nghost)
                                                  || (i + ii >= ni + nghost || j + jj >= nj + nghost || k + kk >= nk + nghost);
                      if (!is_outside_proc)
                        {
                          for (m = 0; m < dim; m++)
                            n[m] += normal_vect[m](i+ii,j+jj,k+kk);

                          if (!tmp_propagated_cells(i+ii,j+jj,k+kk))
                            {
                              const Int3 ijk_index(i+ii, j+jj, k+kk);
                              for (m=0; m<dim; m++)
                                propagated_cells_indices[m].append_array(ijk_index[m]);
                              tmp_propagated_cells(i+ii,j+jj,k+kk) = 1;
                              if (!iteration)
                                for (m=0; m<dim; m++)
                                  gfm_first_cells_indices_[m].append_array(ijk_index[m]);
                            }
                        }
                    }
                  for (m = 0; m < dim; m++)
                    tmp[m](i,j,k) = terme_src[m](i,j,k) + n[m] * un_sur_ncontrib;
                }
              for (int l=0; l<dim; l++)
                normal_vect[l].data() = tmp[l].data();
            }
          else
            {
              for (int k = 0; k < nk; k++)
                for (int j = 0; j < nj; j++)
                  for (int i = 0; i < ni; i++)
                    {
                      // Averaging the normal vector on the neighbours
                      double n[3] = {0., 0., 0.};
                      for (m = 0; m < dim; m++)
                        n[m] = normal_vect[m](i,j,k);
                      for (int l=0; l<6; l++)
                        {
                          const int ii = neighbours_i[l];
                          const int jj = neighbours_j[l];
                          const int kk = neighbours_k[l];
                          for (m = 0; m < dim; m++)
                            n[m] += normal_vect[m](i+ii,j+jj,k+kk);
                        }
                      for (m = 0; m < dim; m++)
                        tmp[m](i,j,k) = terme_src[m](i,j,k) + n[m] * un_sur_ncontrib;
                    }
              for (int l=0; l<dim; l++)
                normal_vect[l].data() = tmp[l].data();
              normal_vect.echange_espace_virtuel();
            }
        }
    }
  else
    {
      for (iteration = 0; iteration < n_iter; iteration++)
        {
          /*
           * Smoothing iterator, in theory:
           * normal = normal + (source_term - laplacian(normal)) * factor
           * converge towards laplacian(normal) = source_term
           * in practice:
           * normal = average(normal on neighbours) + source_term
           */
          const double un_sur_ncontrib = 1. / (1. + nb_elem_voisins);
          int elem, i, k;
          for (elem = 0; elem < nb_elem; elem++)
            {
              // Averaging the normal vector on the neighbours
              double n[3] = {0., 0., 0.};
              const Int3 num_elem_ijk = splitting_distance.convert_packed_to_ijk_cell(elem);
              for (i = 0; i < dim; i++)
                {
                  n[i] = normal_vect[i](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]);
                }
              for (k = 0; k < nb_elem_voisins; k++)
                {
                  // We look for the neighbour by face k
                  const int face = elem_faces(elem, k);
                  const int e_voisin = face_voisins(face, 0) + face_voisins(face, 1) - elem;
                  const Int3 num_elem_voisin_ijk = splitting_distance.convert_packed_to_ijk_cell(e_voisin);
                  if (e_voisin >= 0) // Not on a boundary
                    for (i = 0; i < dim; i++)
                      n[i] += normal_vect[i](num_elem_voisin_ijk[DIRECTION_I],num_elem_voisin_ijk[DIRECTION_J],num_elem_voisin_ijk[DIRECTION_K]);
                }
              for (i = 0; i < dim; i++)
                {
                  tmp[i](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) =
                    terme_src[i](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) + n[i] * un_sur_ncontrib;
                }
            }
          normal_vect = tmp;
          normal_vect.echange_espace_virtuel();
        }
    }
  /*
   * We normalise the normal vector and we create a list of elements
   * for whom the normal is known
   */
  ArrOfIntFT liste_elements;
  if (!use_ijk)
    {
      int elem;
      for (elem = 0; elem < nb_elem; elem++)
        {
          const Int3 num_elem_ijk = splitting_distance.convert_packed_to_ijk_cell(elem);
          double nx = normal_vect[0](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]);
          double ny = normal_vect[1](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]);
          double nz = normal_vect[2](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]);
          double norme2 = nx*nx + ny*ny + nz*nz;
          if (norme2 > 0.)
            {
              double i_norme = 1. / sqrt(norme2);
              normal_vect[0](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) = nx * i_norme;
              normal_vect[1](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) = ny * i_norme;
              normal_vect[2](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) = nz * i_norme;
              liste_elements.append_array(elem);
            }
        }
      normal_vect.echange_espace_virtuel(); // This swap is essential ?
    }
  else
    {
      if (!avoid_gfm_parallel_calls)
        normal_vect.echange_espace_virtuel(); // This swap is essential ?
    }
  // Distance calculation at the interface
  /*
   * TODO: Check the how fast it is compared to using elem_faces matrix
   */
  // IJK_Field_double terme_src_dist(distance_field);
  // IJK_Field_double tmp_dist(distance_field);
  IJK_Field_double& terme_src_dist = tmp_old_val;
  IJK_Field_double& tmp_dist = tmp_new_val;
  terme_src_dist.data() = distance_field.data();
  tmp_dist.data() = distance_field.data();

  tmp_propagated_cells.data() = tmp_interf_cells.data();
  propagated_cells_indices = gfm_first_cells_indices_;
  for (int l=0; l<dim; l++)
    propagated_cells_indices_tmp[l].reset();

  if (use_ijk)
    {
      int neighbours_i[6] = NEIGHBOURS_I;
      int neighbours_j[6] = NEIGHBOURS_J;
      int neighbours_k[6] = NEIGHBOURS_K;
      const int ni = normal_vect[0].ni();
      const int nj = normal_vect[0].nj();
      const int nk = normal_vect[0].nk();
      const int nghost = normal_vect[0].ghost();
      int m;
      for (iteration = 0; iteration < n_iter; iteration++)
        {
          if (avoid_gfm_parallel_calls)
            {
              propagated_cells_indices_tmp = propagated_cells_indices;
              for (int l=0; l<dim; l++)
                propagated_cells_indices[l].reset();
              for (int ielem = 0; ielem < propagated_cells_indices_tmp[0].size_array(); ielem++)
                {
                  const int i = propagated_cells_indices_tmp[DIRECTION_I](ielem);
                  const int j = propagated_cells_indices_tmp[DIRECTION_J](ielem);
                  const int k = propagated_cells_indices_tmp[DIRECTION_K](ielem);

                  // For all the element already crossed by the interface, the value is not computed again
                  if (terme_src_dist(i,j,k) > invalid_distance_value)
                    tmp_dist(i,j,k) = distance_field(i,j,k);
                  else
                    {
                      if (!tmp_propagated_cells(i,j,k))
                        tmp_propagated_cells(i,j,k) = 1;
                      // For the others, we compute a distance value per neighbour
                      double ncontrib = 0.;
                      double somme_distances = 0.;
                      for (int l = 0; l < 6; l++)
                        {
                          // Look for a neighbour
                          const int ii = neighbours_i[l];
                          const int jj = neighbours_j[l];
                          const int kk = neighbours_k[l];

                          const int is_outside_proc = (i + ii < -nghost || j + jj < -nghost || k + kk < -nghost)
                                                      || (i + ii >= ni + nghost || j + jj >= nj + nghost || k + kk >= nk + nghost);
                          if (!is_outside_proc)
                            {
                              const Int3 ijk_index(i+ii, j+jj, k+kk);
                              if (!tmp_propagated_cells(i+ii,j+jj,k+kk))
                                {
                                  for (m=0; m<dim; m++)
                                    propagated_cells_indices[m].append_array(ijk_index[m]);
                                  tmp_propagated_cells(i+ii,j+jj,k+kk) = 1;
                                }

                              const double distance_voisin = distance_field(i+ii, j+jj, k+kk);
                              if (distance_voisin > invalid_distance_value)
                                {
                                  // Average normal distance between an element and its neighbours
                                  double nx = normal_vect[0](i,j,k) + normal_vect[0](i+ii, j+jj, k+kk);
                                  double ny = normal_vect[1](i,j,k) + normal_vect[1](i+ii, j+jj, k+kk);
                                  double nz = normal_vect[2](i,j,k) + normal_vect[2](i+ii, j+jj, k+kk);
                                  double norm2 = nx*nx + ny*ny + nz*nz;
                                  if (norm2 > 0.)
                                    {
                                      double i_norm = 1./sqrt(norm2);
                                      nx *= i_norm;
                                      ny *= i_norm;
                                      nz *= i_norm;
                                    }
                                  // Element to neighbour vector calculation
                                  double dx = - geom.get_constant_delta(DIRECTION_I) * ii;
                                  double dy = - geom.get_constant_delta(DIRECTION_J) * jj;
                                  double dz = - geom.get_constant_delta(DIRECTION_K) * kk;
                                  double d = nx * dx + ny * dy + nz * dz + distance_voisin;
                                  somme_distances += d;
                                  ncontrib++;
                                }
                            }
                        }
                      // Averaging the distances obtained from neighbours
                      if (ncontrib > 0.)
                        {
                          double d = somme_distances / ncontrib;
                          tmp_dist(i,j,k) = d;
                        }
                    }
                }
              distance_field.data() = tmp_dist.data();
            }
          else
            {
              for (int k = 0; k < nk; k++)
                for (int j = 0; j < nj; j++)
                  for (int i = 0; i < ni; i++)
                    {
                      // For all the element already crossed by the interface, the value is not computed again
                      if (terme_src_dist(i,j,k) > invalid_distance_value)
                        tmp_dist(i,j,k) = distance_field(i,j,k);
                      else
                        {
                          // For the others, we compute a distance value per neighbour
                          double ncontrib = 0.;
                          double somme_distances = 0.;
                          for (int l = 0; l < 6; l++)
                            {
                              // Look for a neighbour
                              const int ii = neighbours_i[l];
                              const int jj = neighbours_j[l];
                              const int kk = neighbours_k[l];
                              const double distance_voisin = distance_field(i+ii, j+jj, k+kk);
                              if (distance_voisin > invalid_distance_value)
                                {
                                  // Average normal distance between an element and its neighbours
                                  double nx = normal_vect[0](i,j,k) + normal_vect[0](i+ii, j+jj, k+kk);
                                  double ny = normal_vect[1](i,j,k) + normal_vect[1](i+ii, j+jj, k+kk);
                                  double nz = normal_vect[2](i,j,k) + normal_vect[2](i+ii, j+jj, k+kk);
                                  double norm2 = nx*nx + ny*ny + nz*nz;
                                  if (norm2 > 0.)
                                    {
                                      double i_norm = 1./sqrt(norm2);
                                      nx *= i_norm;
                                      ny *= i_norm;
                                      nz *= i_norm;
                                    }
                                  // Element to neighbour vector calculation
                                  double dx = - geom.get_constant_delta(DIRECTION_I) * ii;
                                  double dy = - geom.get_constant_delta(DIRECTION_J) * jj;
                                  double dz = - geom.get_constant_delta(DIRECTION_K) * kk;
                                  double d = nx * dx + ny * dy + nz * dz + distance_voisin;
                                  somme_distances += d;
                                  ncontrib++;
                                }
                            }
                          // Averaging the distances obtained from neighbours
                          if (ncontrib > 0.)
                            {
                              double d = somme_distances / ncontrib;
                              tmp_dist(i,j,k) = d;
                            }
                        }
                    }
              distance_field.data() = tmp_dist.data();
              distance_field.echange_espace_virtuel(distance_field.ghost());
            }
        }
    }
  else
    {
      for (iteration = 0; iteration < n_iter; iteration++)
        {
          int i_elem, elem;
          const int liste_elem_size = liste_elements.size_array();
          for (i_elem = 0; i_elem < liste_elem_size; i_elem++)
            {
              elem = liste_elements[i_elem];
              const Int3 num_elem_ijk = splitting_distance.convert_packed_to_ijk_cell(elem);
              if (terme_src_dist(num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) > invalid_distance_value)
                {
                  // For all the element already crossed by the interface, the value is not computed again
                  tmp_dist(num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) = distance_field(num_elem_ijk[DIRECTION_I],
                                                                                                                             num_elem_ijk[DIRECTION_J],
                                                                                                                             num_elem_ijk[DIRECTION_K]);
                }
              else
                {
                  // For the others, we compute a distance value per neighbour
                  double ncontrib = 0.;
                  double somme_distances = 0.;
                  int k;
                  for (k = 0; k < nb_elem_voisins; k++)
                    {
                      // Look for a neighbour by the face k
                      const int face = elem_faces(elem, k);
                      const int e_voisin = face_voisins(face, 0) + face_voisins(face, 1) - elem;
                      const Int3 num_elem_voisin_ijk = splitting_distance.convert_packed_to_ijk_cell(e_voisin);
                      if (e_voisin >= 0) // Not on a boundary
                        {
                          const double distance_voisin = distance_field(num_elem_voisin_ijk[DIRECTION_I],
                                                                        num_elem_voisin_ijk[DIRECTION_J],
                                                                        num_elem_voisin_ijk[DIRECTION_K]);
                          if (distance_voisin > invalid_distance_value)
                            {
                              // Average normal distance between an element and its neighbours
                              double nx = normal_vect[0](num_elem_ijk[DIRECTION_I],
                                                         num_elem_ijk[DIRECTION_J],
                                                         num_elem_ijk[DIRECTION_K]) +
                                          normal_vect[0](num_elem_voisin_ijk[DIRECTION_I],
                                                         num_elem_voisin_ijk[DIRECTION_J],
                                                         num_elem_voisin_ijk[DIRECTION_K]);
                              double ny = normal_vect[1](num_elem_ijk[DIRECTION_I],
                                                         num_elem_ijk[DIRECTION_J],
                                                         num_elem_ijk[DIRECTION_K]) +
                                          normal_vect[1](num_elem_voisin_ijk[DIRECTION_I],
                                                         num_elem_voisin_ijk[DIRECTION_J],
                                                         num_elem_voisin_ijk[DIRECTION_K]);
                              double nz = normal_vect[2](num_elem_ijk[DIRECTION_I],
                                                         num_elem_ijk[DIRECTION_J],
                                                         num_elem_ijk[DIRECTION_K]) +
                                          normal_vect[2](num_elem_voisin_ijk[DIRECTION_I],
                                                         num_elem_voisin_ijk[DIRECTION_J],
                                                         num_elem_voisin_ijk[DIRECTION_K]);
                              double norm2 = nx*nx + ny*ny + nz*nz;
                              if (norm2 > 0.)
                                {
                                  double i_norm = 1./sqrt(norm2);
                                  nx *= i_norm;
                                  ny *= i_norm;
                                  nz *= i_norm;
                                }
                              // Element to neighbour vector calculation
                              double dx = centre_element(elem, 0) - centre_element(e_voisin, 0);
                              double dy = centre_element(elem, 1) - centre_element(e_voisin, 1);
                              double dz = centre_element(elem, 2) - centre_element(e_voisin, 2);
                              double d = nx * dx + ny * dy + nz * dz + distance_voisin;
                              somme_distances += d;
                              ncontrib++;
                            }
                        }
                    }
                  // Averaging the distances obtained from neighbours
                  if (ncontrib > 0.)
                    {
                      double d = somme_distances / ncontrib;
                      tmp_dist(num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) = d;
                    }
                }
            }
          distance_field = tmp_dist;
          distance_field.echange_espace_virtuel(distance_field.ghost());
        }
    }
  if (avoid_gfm_parallel_calls)
    {
      normal_vect.echange_espace_virtuel();
      distance_field.echange_espace_virtuel(distance_field.ghost());
    }

  statistiques().end_count(stat_counter);
}

void compute_eulerian_curvature_field_from_distance_field(const IJK_Field_double& distance,
                                                          IJK_Field_double& curvature,
                                                          const IJK_Field_local_double& boundary_flux_kmin,
                                                          const IJK_Field_local_double& boundary_flux_kmax)
{
  /*
   * Compute the divergence of the normal vector field or the laplacian of the eulerian distance field
   */
  static const Stat_Counter_Id stat_counter = statistiques().new_counter(2, "GFM - Compute Eulerian curvature dield from distance field");
  statistiques().begin_count(stat_counter);

  // Laplacian operator
  Operateur_IJK_elem_diff laplacian_distance;
  laplacian_distance.typer_diffusion_op("uniform");
  // Initialise with unit lambda
  laplacian_distance.initialize(distance.get_splitting());
  const double lambda = 1.;
  laplacian_distance->set_uniform_lambda(lambda);
  // Calculate Laplacian(dist)
  laplacian_distance->calculer(distance, curvature, boundary_flux_kmin, boundary_flux_kmax);
  const IJK_Grid_Geometry& geom = curvature.get_splitting().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  const double vol = dx * dy * dz;
  const int nx = curvature.ni();
  const int ny = curvature.nj();
  const int nz = curvature.nk();
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        curvature(i,j,k) /= vol;
  curvature.echange_espace_virtuel(curvature.ghost());

  statistiques().end_count(stat_counter);
}

void compute_eulerian_curvature_field_from_normal_vector_field(const IJK_Field_vector3_double& normal_vect,
                                                               IJK_Field_double& curvature)
{

}

void compute_eulerian_curvature_field_from_interface(const IJK_Field_vector3_double& normal_vect,
                                                     const IJK_Interfaces& interfaces,
                                                     IJK_Field_double& interfacial_area,
                                                     IJK_Field_double& curvature,
                                                     IJK_Field_double& tmp_old_val,
                                                     IJK_Field_double& tmp_new_val,
                                                     const int& n_iter,
                                                     const int igroup)
{
  /*
   * From IJK_Interfaces::calculer_normales_et_aires_interfaciales
   * Called in update_stat_ft IJK_FT through an instance of IJK_FT_Post !!!!!
   */
  const bool use_ijk = true;
  // Vertex coordinates of the eulerian domain
  interfacial_area.echange_espace_virtuel(interfacial_area.ghost());
  curvature.echange_espace_virtuel(curvature.ghost());

  static const Stat_Counter_Id stat_counter = statistiques().new_counter(2, "GFM - Compute Eulerian Curvature field from interface");
  statistiques().begin_count(stat_counter);

  static const double invalid_curvature_value = INVALID_TEST;

  interfacial_area.data() = invalid_curvature_value * 1.1;
  curvature.data() = invalid_curvature_value * 1.1;

  // Vertex coordinates of the eulerian domain
  const Domaine_dis_base& mon_dom_dis = interfaces.get_domaine_dis();
  const int nb_elem = mon_dom_dis.domaine().nb_elem();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, mon_dom_dis);
  const IJK_Splitting& splitting_curvature = normal_vect.get_splitting();

  const Maillage_FT_IJK& maillage = interfaces.maillage_ft_ijk();
  const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  const IntTab& facettes = maillage.facettes();
  const ArrOfDouble& courbure = maillage.get_update_courbure_sommets();

  const int n_fa7 = maillage.nb_facettes();
  // Calculate the curvature in the cells crossed by the interface
  const ArrOfInt& compo_facette = maillage.compo_connexe_facettes();
  const ArrOfInt& compo_to_group = interfaces.get_compo_to_group();
  for (int fa7 = 0; fa7 < n_fa7; fa7++)
    {
      int icompo = compo_facette[fa7];
      if (icompo<0)
        {
          // Portion d'interface ghost. On recherche le vrai numero
          icompo = decoder_numero_bulle(-icompo);
        }
      if ((compo_to_group[icompo] != igroup) && (igroup != -1))
        continue;
      const double sf = surface_facettes[fa7];
      int index = intersections.index_facette()[fa7];
      while (index >= 0)
        {
          const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
          const int num_elem = data.numero_element_;
          const Int3 ijk = splitting_curvature.convert_packed_to_ijk_cell(num_elem);
          const double surf = data.fraction_surface_intersection_ * sf;
          for (int isom = 0; isom < 3; isom++)
            {
              const int num_som = facettes(fa7, isom);
              const double kappa = courbure[num_som];
              // No volume consideration
              if (curvature(ijk[DIRECTION_I], ijk[DIRECTION_J], ijk[DIRECTION_K]) < invalid_curvature_value)
                curvature(ijk[DIRECTION_I], ijk[DIRECTION_J], ijk[DIRECTION_K]) = kappa * surf / 3.;
              else
                curvature(ijk[DIRECTION_I], ijk[DIRECTION_J], ijk[DIRECTION_K]) += kappa * surf / 3.;
            }
          if (interfacial_area(ijk[DIRECTION_I], ijk[DIRECTION_J], ijk[DIRECTION_K]) < invalid_curvature_value)
            interfacial_area(ijk[DIRECTION_I], ijk[DIRECTION_J], ijk[DIRECTION_K]) = surf;
          else
            interfacial_area(ijk[DIRECTION_I], ijk[DIRECTION_J], ijk[DIRECTION_K]) += surf;
          index = data.index_element_suivant_;
        }
    }
  interfacial_area.echange_espace_virtuel(interfacial_area.ghost());
  curvature.echange_espace_virtuel(curvature.ghost());

  {
    const int nx = curvature.ni();
    const int ny = curvature.nj();
    const int nz = curvature.nk();
    for (int k=0; k < nz ; k++)
      for (int j=0; j< ny; j++)
        for (int i=0; i < nx; i++)
          {
            const double kappa = curvature(i,j,k);
            const double ai = interfacial_area(i,j,k);
            // TODO: Why do I get a floating point exception ? Interface portion too small ?
            if ((kappa > invalid_curvature_value) && (ai > invalid_curvature_value))
              {
                if (fabs(ai) < DMINFLOAT)
                  {
                    Cerr << "Interfacial_area is very much at zero... Pathological case to be looked into closely. " << finl;
                    Cerr << "Curvature is set to invalid value to be overwritten by its neighbours" << finl;
                    // Be careful if the distance is not calculated well the spreading algorithm will not work
                    curvature(i,j,k) = invalid_curvature_value * 1.1;
                    // Process::exit();
                  }
                else
                  {
                    curvature(i,j,k) = kappa / ai;
                  }
              }
          }
  }
  interfacial_area.echange_espace_virtuel(interfacial_area.ghost());
  curvature.echange_espace_virtuel(curvature.ghost());

  // Get back the cells filled with non-zero normal vectors
  ArrOfIntFT liste_elements;
  if (!use_ijk)
    {
      for (int elem = 0; elem < nb_elem; elem++)
        {
          const Int3 num_elem_ijk = splitting_curvature.convert_packed_to_ijk_cell(elem);
          double nx = normal_vect[0](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]);
          double ny = normal_vect[1](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]);
          double nz = normal_vect[2](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]);
          double norme2 = nx*nx + ny*ny + nz*nz;
          if (norme2 > 0.)
            liste_elements.append_array(elem);
        }
    }

  // Curvature calculation at the interface
  // IJK_Field_double terme_src_curv(curvature);
  // IJK_Field_double tmp_curv(curvature);

  IJK_Field_double& terme_src_curv = tmp_old_val;
  IJK_Field_double& tmp_curv = tmp_new_val;
  terme_src_curv.data() = curvature.data();
  tmp_curv.data() = curvature.data();

  /*
   * TODO: Check the how fast it is compared to using elem_faces matrix
   */
  if (use_ijk)
    {
      int neighbours_i[6] = NEIGHBOURS_I;
      int neighbours_j[6] = NEIGHBOURS_J;
      int neighbours_k[6] = NEIGHBOURS_K;
      const int ni = normal_vect[0].ni();
      const int nj = normal_vect[0].nj();
      const int nk = normal_vect[0].nk();
      for (int iteration = 0; iteration < n_iter; iteration++)
        {
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              for (int i = 0; i < ni; i++)
                {
                  // For all the element already crossed by the interface, the value is not computed again
                  if (terme_src_curv(i,j,k) > invalid_curvature_value)
                    tmp_curv(i,j,k) = curvature(i,j,k);
                  else
                    {
                      // For the others, we compute a distance value per neighbour
                      double ncontrib = 0.;
                      double sum_kappa = 0.;
                      for (int l = 0; l < 6; l++)
                        {
                          // Look for a neighbour
                          const int ii = neighbours_i[l];
                          const int jj = neighbours_j[l];
                          const int kk = neighbours_k[l];
                          const double curvature_voisin = curvature(i+ii,j+jj,k+kk);
                          if (curvature_voisin > invalid_curvature_value)
                            {
                              // Average normal distance between an element and its neighbours
                              sum_kappa += curvature_voisin;
                              ncontrib++;
                            }
                        }
                      // Averaging the distances obtained from neighbours
                      if (ncontrib > 0.)
                        {
                          double kappa = sum_kappa / ncontrib;
                          tmp_curv(i,j,k) = kappa;
                        }
                    }
                }
          curvature.data() = tmp_curv.data();
          curvature.echange_espace_virtuel(curvature.ghost());
        }
    }
  else
    {
      const IntTab& face_voisins = domaine_vf.face_voisins();
      const IntTab& elem_faces   = domaine_vf.elem_faces();
      const int nb_elem_voisins = elem_faces.line_size();
      for (int iteration = 0; iteration < n_iter; iteration++)
        {
          int i_elem, elem;
          const int liste_elem_size = liste_elements.size_array();
          for (i_elem = 0; i_elem < liste_elem_size; i_elem++)
            {
              elem = liste_elements[i_elem];
              const Int3 num_elem_ijk = splitting_curvature.convert_packed_to_ijk_cell(elem);
              if (terme_src_curv(num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) > invalid_curvature_value)
                {
                  // For all the element already crossed by the interface, the value is not computed again
                  tmp_curv(num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) = curvature(num_elem_ijk[DIRECTION_I],
                                                                                                                        num_elem_ijk[DIRECTION_J],
                                                                                                                        num_elem_ijk[DIRECTION_K]);
                }
              else
                {
                  // For the others, we compute a distance value per neighbour
                  double ncontrib = 0.;
                  double sum_kappa = 0.;
                  int k;
                  for (k = 0; k < nb_elem_voisins; k++)
                    {
                      // Look for a neighbour by the face k
                      const int face = elem_faces(elem, k);
                      const int e_voisin = face_voisins(face, 0) + face_voisins(face, 1) - elem;
                      const Int3 num_elem_voisin_ijk = splitting_curvature.convert_packed_to_ijk_cell(e_voisin);
                      if (e_voisin >= 0) // Not on a boundary
                        {
                          const double curvature_voisin = curvature(num_elem_voisin_ijk[DIRECTION_I],
                                                                    num_elem_voisin_ijk[DIRECTION_J],
                                                                    num_elem_voisin_ijk[DIRECTION_K]);
                          if (curvature_voisin > invalid_curvature_value)
                            {
                              // Average normal distance between an element and its neighbours
                              sum_kappa += curvature_voisin;
                              ncontrib++;
                            }
                        }
                    }
                  // Averaging the distances obtained from neighbours
                  if (ncontrib > 0.)
                    {
                      double kappa = sum_kappa / ncontrib;
                      tmp_curv(num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) = kappa;
                    }
                }
            }
          curvature = tmp_curv;
          curvature.echange_espace_virtuel(curvature.ghost());
        }
    }
  statistiques().end_count(stat_counter);
}

void compute_eulerian_normal_temperature_gradient_interface(const IJK_Field_double& distance,
                                                            const IJK_Field_double& indicator,
                                                            const IJK_Field_double& interfacial_area,
                                                            const IJK_Field_double& curvature,
                                                            const	IJK_Field_double& temperature,
                                                            IJK_Field_double& grad_T_interface,
                                                            const int& spherical_approx,
                                                            const double& temperature_interf)
{
  /*
   * Compute the normal temperature gradient at the bubble interface
   * Write in the ijk manner !
   */
  static const Stat_Counter_Id stat_counter = statistiques().new_counter(3, "GFM - Compute Eulerian normal temperature gradient interface");
  statistiques().begin_count(stat_counter);

  int neighbours_i[6] = NEIGHBOURS_I;
  int neighbours_j[6] = NEIGHBOURS_J;
  int neighbours_k[6] = NEIGHBOURS_K;
  static const double invalid_value = INVALID_TEST;
  static const double liquid_indicator = LIQUID_INDICATOR_TEST;
  const int ni = temperature.ni();
  const int nj = temperature.nj();
  const int nk = temperature.nk();

  //  grad_T_interface.data() = 1.1 * invalid_value;
  grad_T_interface.data() = 0.;
  grad_T_interface.echange_espace_virtuel(grad_T_interface.ghost());

  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double ai = interfacial_area(i,j,k);
          if (ai > invalid_value)
            {
              for (int l=0; l < 6; l++)
                {
                  const int ii = neighbours_i[l];
                  const int jj = neighbours_j[l];
                  const int kk = neighbours_k[l];
                  const double d = distance(i+ii,j+jj,k+kk);
                  const double indic = indicator(i+ii,j+jj,k+kk);
                  // if ((indic > liquid_indicator) && (d > invalid_value) && grad_T_interface(i+ii,j+jj,k+kk) < invalid_value)
                  if ((indic > liquid_indicator) && (d > invalid_value) && grad_T_interface(i+ii,j+jj,k+kk) == 0)
                    {
                      const double temperature_liquid = temperature(i+ii,j+jj,k+kk);
                      const double second_order_gradient = (temperature_liquid - temperature_interf) / d;
                      const double kappa = curvature(i+ii,j+jj,k+kk);
                      double grad_T_modified = 0.;
                      // TODO: Check sign kappa
                      if (spherical_approx)
                        grad_T_modified = second_order_gradient * (1. - 0.5 * kappa * d);
                      else
                        {
                          const double kappa_non_zero = kappa + 1.e-16;
                          grad_T_modified = pow((d / 2 - 2 / kappa_non_zero),2) * second_order_gradient / (pow((0. - 2 / kappa_non_zero),2) + 1e-16);
                        }
                      grad_T_interface(i+ii,j+jj,k+kk) = grad_T_modified;
                    }
                }
            }
        }
  grad_T_interface.echange_espace_virtuel(grad_T_interface.ghost());
  statistiques().end_count(stat_counter);
  /*
   * Check if indicatrice of the neighbours is zero + interfacial_area to locate the mixed cells
   */
}

void propagate_eulerian_normal_temperature_gradient_interface(const IJK_Interfaces& interfaces,
                                                              const IJK_Field_double& distance,
                                                              IJK_Field_double& grad_T_interface,
                                                              const int& stencil_width,
                                                              const int& recompute_field_ini,
                                                              const int& zero_neighbour_value_mean,
                                                              const int& vapour_mixed_only,
                                                              const int& smooth_factor)
{
  /*
   * Propagate value of grad_T_int stored in pure liquid phase towards the vapour phase and mixed cells
   */
  const bool use_ijk = true;
  if (use_ijk)
    {
      // Using the ijk indices
      extrapolate_with_ijk_indices(distance,
                                   interfaces.I_ft(),
                                   grad_T_interface,
                                   stencil_width,
                                   recompute_field_ini,
                                   zero_neighbour_value_mean,
                                   vapour_mixed_only,
                                   smooth_factor);
    }
  else
    {
      // Using the elem_faces connectivity matrix of TrioCFD
      const Domaine_dis_base& mon_dom_dis = interfaces.get_domaine_dis();
      const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, mon_dom_dis);
      extrapolate_with_elem_faces_connectivity(domaine_vf, distance, grad_T_interface, stencil_width);
    }
}

void compute_eulerian_extended_temperature(const IJK_Field_double& indicator,
                                           const IJK_Field_double& distance,
                                           const IJK_Field_double& curvature,
                                           IJK_Field_double& grad_T_interface,
                                           IJK_Field_double& temperature,
                                           const int& spherical_approx,
                                           const double& temperature_interf)
{
  /*
   * Compute the extended temperature field using propagated values of the temperature gradient
   */
  static const Stat_Counter_Id stat_counter = statistiques().new_counter(3, "GFM - Compute Eulerian ghost fluid temperature extension");
  statistiques().begin_count(stat_counter);

  const double invalid_test = INVALID_TEST;
  const int ni = temperature.ni();
  const int nj = temperature.nj();
  const int nk = temperature.nk();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double indicator_vapour = fabs(1 - indicator(i,j,k));
          const double d = distance(i,j,k);
          const double grad_T = grad_T_interface(i,j,k);
          const double temperature_val = temperature(i,j,k);
          if ((d > invalid_test) && (indicator_vapour > VAPOUR_INDICATOR_TEST) && (grad_T != 0) && (temperature_val == 0))
            {
              const double kappa = curvature(i,j,k);
              double temperature_ghost = 0.;
              if (spherical_approx)
                temperature_ghost = temperature_interf + d * grad_T * (1. + 0.5 * kappa * d + kappa * kappa * d * d / 6.);
              else
                {
                  const double kappa_non_zero = kappa + 1.e-16;
                  temperature_ghost = temperature_interf + grad_T * (- 2 / kappa_non_zero) * (1 - (- 2 / kappa_non_zero) / ((d - (2 / kappa_non_zero)) + 1e-16));
                }
              temperature(i,j,k) = temperature_ghost;
            }
        }
  temperature.echange_espace_virtuel(temperature.ghost());

  statistiques().end_count(stat_counter);
}

void smooth_vector_field(IJK_Field_vector3_double& vector_field,
                         IJK_Field_vector3_double& vector_field_init,
                         const IJK_Field_vector3_double * eulerian_normal_vectors_ns_normed,
                         const IJK_Interfaces& interfaces,
                         const double (&direct_smoothing_factors) [7],
                         const double (&gaussian_smoothing_factors) [3][3][3],
                         const int& smooth_numbers,
                         const int& remove_normal_compo,
                         const int& direct_neighbours,
                         const int& use_field_init,
                         const int& use_unique_phase)
{
  for (int c=0; c<3; c++)
    {
      IJK_Field_double& field = vector_field[c];
      IJK_Field_double& field_init = vector_field_init[c];
      smooth_eulerian_field(field,
                            field_init,
                            c,
                            vector_field_init,
                            eulerian_normal_vectors_ns_normed,
                            interfaces,
                            direct_smoothing_factors,
                            gaussian_smoothing_factors,
                            smooth_numbers,
                            remove_normal_compo,
                            direct_neighbours,
                            use_field_init,
                            use_unique_phase);
    }
}

void smooth_eulerian_field(IJK_Field_double& field,
                           IJK_Field_double& field_init,
                           const int& dir,
                           IJK_Field_vector3_double& vector_field_init,
                           const IJK_Field_vector3_double * eulerian_normal_vectors_ns_normed,
                           const IJK_Interfaces& interfaces,
                           const double (&direct_smoothing_factors) [7],
                           const double (&gaussian_smoothing_factors) [3][3][3],
                           const int& smooth_numbers,
                           const int& remove_normal_compo,
                           const int& direct_neighbours,
                           const int& use_field_init,
                           const int& use_unique_phase)
{
  const IJK_Field_double& indicator = interfaces.I();
  const int smooth_numbers_end = (smooth_numbers < 1) ? 1: smooth_numbers;
  const int use_field_init_usr = use_field_init && (smooth_numbers_end == 1);
  IJK_Field_double field_copy;
  for (int m=0; m<smooth_numbers_end; m++)
    {
      const int ni = field.ni();
      const int nj = field.nj();
      const int nk = field.nk();
      IJK_Field_double& field_raw = use_field_init_usr ? field_init : field_copy;
      if (!use_field_init_usr)
        {
          if (m == 0)
            field_copy = field_init;
          else
            field_copy.data() = field.data();
          field_copy.echange_espace_virtuel(field_copy.ghost());
        }
      const int neighbours_i[6] = NEIGHBOURS_I;
      const int neighbours_j[6] = NEIGHBOURS_J;
      const int neighbours_k[6] = NEIGHBOURS_K;

      double sum_factors = 0;
      double sum_factors_phase = 0;
      double sum_direct_smoothing_factors = 0.;
      double sum_gaussian_smoothing_factors = 0.;

      for (int c=0; c<7; c++)
        sum_direct_smoothing_factors += direct_smoothing_factors[c];
      for (int c=0; c<3; c++)
        for (int l=0; l<3; l++)
          for (int n=0; n<3; n++)
            sum_gaussian_smoothing_factors += gaussian_smoothing_factors[n][l][c];

      sum_factors = (direct_neighbours) ? sum_direct_smoothing_factors : sum_gaussian_smoothing_factors;
      sum_factors_phase = sum_factors;
      Cerr << "Sum of smoothing factors: " << sum_factors << finl;
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              if (direct_neighbours)
                {
                  field(i,j,k) = direct_smoothing_factors[6] * field_raw(i,j,k);
                  for (int l=0; l<6; l++)
                    {
                      const int ii = neighbours_i[l];
                      const int jj = neighbours_j[l];
                      const int kk = neighbours_k[l];
                      if (!remove_normal_compo)
                        field(i,j,k) += direct_smoothing_factors[l] * field_raw(i+ii,j+jj,k+kk);
                      else
                        {
                          const double normal_vect = (*eulerian_normal_vectors_ns_normed)[dir](i,j,k);
                          const double normal_compo = (field_raw(i+ii,j+jj,k+kk) *
                                                       direct_smoothing_factors[l] *
                                                       normal_vect);
                          // const double normal_vect = (*eulerian_normal_vectors_ns_normed)[dir](i+ii,j+jj,k+kk);
                          // const double normal_compo = (field_raw(i+ii,j+jj,k+kk) *
                          //                              direct_smoothing_factors[l] *
                          //                              (*eulerian_normal_vectors_ns_normed)[dir](i+ii,j+jj,k+kk));
                          field(i,j,k) += (direct_smoothing_factors[l] * field_raw(i+ii,j+jj,k+kk) - normal_compo);
                        }
                    }
                }
              else
                {
                  sum_factors_phase = use_unique_phase ? 0 : sum_factors;
                  bool test_indicator;
                  field(i,j,k) = 0.;
                  // for (int c=0; c<3; c++)
                  for (int c=-1; c<=1; c++)
                    for (int l=-1; l<=1; l++)
                      for (int n=-1; n<=1; n++)
                        {
                          test_indicator=true;
                          // const int ii = select_dir(c, l, 0, 0);
                          // const int jj = select_dir(c, 0, l, 0);
                          // const int kk = select_dir(c, 0, 0, l);
                          const int ii = n;
                          const int jj = l;
                          const int kk = c;
                          const double indic = indicator(i+ii,j+jj,k+kk);
                          if (use_unique_phase && indic < VAPOUR_INDICATOR_TEST)
                            test_indicator = false;
                          if (test_indicator)
                            {
                              const double smoothing_factor_index = gaussian_smoothing_factors[n+1][l+1][c+1];
                              if (!remove_normal_compo)
                                field(i,j,k) += smoothing_factor_index * field_raw(i+ii,j+jj,k+kk);
                              else
                                {
                                  const double normal_vect = (*eulerian_normal_vectors_ns_normed)[dir](i,j,k);
//                                  const double normal_compo = (smoothing_factor_index *
//                                                               field_raw(i+ii,j+jj,k+kk) *
//                                                               normal_vect * normal_vect);
                                  double normal_compo_tot = 0;
                                  for (int cc=0; cc<3; cc++)
                                    {
                                      const double normal_vect_dir = (*eulerian_normal_vectors_ns_normed)[cc](i,j,k);
                                      normal_compo_tot += (smoothing_factor_index *
                                                           vector_field_init[cc](i+ii,j+jj,k+kk) *
                                                           normal_vect_dir * normal_vect);
                                    }
                                  // normal_compo_tot = normal_compo;
                                  // const double sign_projection = signbit(normal_vect) ? -1.: 1.;
                                  // const double normal_compo = (smoothing_factor_index *
                                  //                             field_raw(i+ii,j+jj,k+kk) *
                                  //                             (*eulerian_normal_vectors_ns_normed)[dir](i+ii,j+jj,k+kk));
                                  field(i,j,k) += (smoothing_factor_index * field_raw(i+ii,j+jj,k+kk) - normal_compo_tot);
                                  if (use_unique_phase)
                                    sum_factors_phase += smoothing_factor_index;
                                }
                            }
                        }
                }
              if (sum_factors_phase != 0)
                {
                  field(i,j,k) /= sum_factors_phase;

                  if (remove_normal_compo)
                    {
                      const double normal_vect = (*eulerian_normal_vectors_ns_normed)[dir](i,j,k);
                      double normal_compo_tot = 0;
                      for (int cc=0; cc<3; cc++)
                        {
                          const double normal_vect_dir = (*eulerian_normal_vectors_ns_normed)[cc](i,j,k);
                          normal_compo_tot += vector_field_init[cc](i,j,k) * normal_vect_dir;
                        }
                      // normal_compo_tot = field_raw(i,j,k) * normal_vect;
                      const double normal_compo = normal_compo_tot * normal_vect;
                      field(i,j,k) += normal_compo;
                    }
                }
            }
    }
}

void fill_tangential_gradient(const IJK_Field_vector3_double& vector_field,
                              const IJK_Field_vector3_double * eulerian_normal_vectors_ns_normed,
                              IJK_Field_vector3_double& tangential_vector_field)
{
  for (int c=0; c<3; c++)
    {
      IJK_Field_double& field = tangential_vector_field[c];
      fill_tangential_gradient_compo(vector_field,
                                     eulerian_normal_vectors_ns_normed,
                                     field, c);
    }
}

void fill_tangential_gradient_compo(const IJK_Field_vector3_double& vector_field,
                                    const IJK_Field_vector3_double * eulerian_normal_vectors_ns_normed,
                                    IJK_Field_double& tangential_field,
                                    const int& dir)
{
  const int ni = tangential_field.ni();
  const int nj = tangential_field.nj();
  const int nk = tangential_field.nk();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double normal_vect = (*eulerian_normal_vectors_ns_normed)[dir](i,j,k);
          double normal_compo_tot = 0;
          for (int cc=0; cc<3; cc++)
            {
              const double normal_vect_dir = (*eulerian_normal_vectors_ns_normed)[cc](i,j,k);
              normal_compo_tot += vector_field[cc](i,j,k) * normal_vect_dir * normal_vect;
            }
          tangential_field(i,j,k) = vector_field[dir](i,j,k) - normal_compo_tot;
        }
}
