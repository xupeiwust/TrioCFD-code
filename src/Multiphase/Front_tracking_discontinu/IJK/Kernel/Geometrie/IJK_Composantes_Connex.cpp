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
// File      : IJK_Composantes_Connex.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/FT
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Composantes_Connex.h>
#include <Probleme_FTD_IJK.h>
#include <IJK_Interfaces.h>
#include <IJK_Bubble_tools.h>

Implemente_instanciable( IJK_Composantes_Connex, "IJK_Composantes_Connex", Objet_U ) ;

static int decoder_numero_bulle(const int code)
{
  const int num_bulle = code >> 6;
  return num_bulle;
}

Sortie& IJK_Composantes_Connex::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  return os;
}

Entree& IJK_Composantes_Connex::readOn( Entree& is )
{
  Objet_U::readOn( is );
  return is;
}

void IJK_Composantes_Connex::initialize(IJK_Interfaces& interfaces,
                                        const bool is_switch)
{
  is_switch_ = is_switch;
  if (!is_switch)
    interfaces_ = &interfaces;
}

void IJK_Composantes_Connex::allocate_fields(const Domaine_IJK& splitting,
                                             const int& allocate_compo_fields)
{
  compute_compo_fields_ = allocate_compo_fields;
  if (!is_switch_ && allocate_compo_fields)
    {
      if (Process::nproc() == 1 && compute_from_bounding_box_)
        {
          eulerian_compo_connex_ft_.allocate(ref_ijk_ft_->get_domaine_ft(), Domaine_IJK::ELEM, 2);
          eulerian_compo_connex_ft_.data() = -1.;
          eulerian_compo_connex_ft_.echange_espace_virtuel(eulerian_compo_connex_ft_.ghost());

          eulerian_compo_connex_ghost_ft_.allocate(ref_ijk_ft_->get_domaine_ft(), Domaine_IJK::ELEM, 2);
          eulerian_compo_connex_ghost_ft_.data() = -1.;
          eulerian_compo_connex_ghost_ft_.echange_espace_virtuel(eulerian_compo_connex_ghost_ft_.ghost());

          eulerian_compo_connex_ns_.allocate(splitting, Domaine_IJK::ELEM, 0);
          eulerian_compo_connex_ns_.echange_espace_virtuel(eulerian_compo_connex_ns_.ghost());

          eulerian_compo_connex_ghost_ns_.allocate(splitting, Domaine_IJK::ELEM, 0);
          eulerian_compo_connex_ghost_ns_.echange_espace_virtuel(eulerian_compo_connex_ghost_ns_.ghost());
        }
      eulerian_compo_connex_from_interface_ft_.allocate(ref_ijk_ft_->get_domaine_ft(), Domaine_IJK::ELEM, 0);
      eulerian_compo_connex_from_interface_ft_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ft_.ghost());

      eulerian_compo_connex_from_interface_ns_.allocate(splitting, Domaine_IJK::ELEM, 0);
      eulerian_compo_connex_from_interface_ns_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ns_.ghost());

      eulerian_compo_connex_from_interface_ghost_ft_.allocate(ref_ijk_ft_->get_domaine_ft(), Domaine_IJK::ELEM, 0);
      eulerian_compo_connex_from_interface_ghost_ft_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ghost_ft_.ghost());

      eulerian_compo_connex_from_interface_ghost_ns_.allocate(splitting, Domaine_IJK::ELEM, 0);
      eulerian_compo_connex_from_interface_ghost_ns_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ghost_ns_.ghost());

      eulerian_compo_connex_from_interface_int_ns_.allocate(splitting, Domaine_IJK::ELEM, 1);
      eulerian_compo_connex_from_interface_int_ns_.data() = -1;
      eulerian_compo_connex_from_interface_int_ns_.echange_espace_virtuel(eulerian_compo_connex_from_interface_int_ns_.ghost());

      eulerian_compo_connex_from_interface_ghost_int_ns_.allocate(splitting, Domaine_IJK::ELEM, 1);
      eulerian_compo_connex_from_interface_ghost_int_ns_.data() = -1;
      eulerian_compo_connex_from_interface_ghost_int_ns_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ghost_int_ns_.ghost());

      eulerian_compo_connex_valid_compo_field_.allocate(splitting, Domaine_IJK::ELEM, 1);
      eulerian_compo_connex_valid_compo_field_.data() = 0;
      eulerian_compo_connex_valid_compo_field_.echange_espace_virtuel(eulerian_compo_connex_valid_compo_field_.ghost());
    }
}

void IJK_Composantes_Connex::associer(const Probleme_FTD_IJK_base& ijk_ft)
{
  ref_ijk_ft_ = ijk_ft;
}

void IJK_Composantes_Connex::initialise_bubbles_params()
{
  interfaces_->calculer_volume_bulles(bubbles_volume_, bubbles_barycentre_);
}

void IJK_Composantes_Connex::associate_rising_velocities_parameters(const Domaine_IJK& splitting,
                                                                    const int& compute_rising_velocities,
                                                                    const int& fill_rising_velocities,
                                                                    const int& use_bubbles_velocities_from_interface,
                                                                    const int& use_bubbles_velocities_from_barycentres)
{
  if (compute_compo_fields_)
    {
      compute_rising_velocities_ = compute_rising_velocities;
      fill_rising_velocities_ = fill_rising_velocities;
      use_bubbles_velocities_from_interface_ = use_bubbles_velocities_from_interface;
      use_bubbles_velocities_from_barycentres_ = use_bubbles_velocities_from_barycentres;
      if (use_bubbles_velocities_from_barycentres_)
        use_bubbles_velocities_from_interface_ = 0;
      if (fill_rising_velocities_)
        {
          eulerian_rising_velocities_.allocate(splitting, Domaine_IJK::ELEM, 0);
          eulerian_rising_velocities_.data() = 0;
          eulerian_rising_velocities_.echange_espace_virtuel(eulerian_rising_velocities_.ghost());
        }
    }
}

void IJK_Composantes_Connex::compute_bounding_box_fill_compo_connex()
{
  static Stat_Counter_Id cnt_compo_connex_bounding_box = statistiques().new_counter(2, "Compo Connex - Bounding Box");
  statistiques().begin_count(cnt_compo_connex_bounding_box);

  if (compute_compo_fields_)
    {
      if (Process::nproc() != 1 || !compute_from_bounding_box_)
        interfaces_->calculer_bounding_box_bulles(bounding_box_);
      else
        {
          compute_bounding_box_fill_compo(*interfaces_,
                                          bounding_box_,
                                          min_max_larger_box_,
                                          eulerian_compo_connex_ft_,
                                          eulerian_compo_connex_ghost_ft_,
                                          bubbles_barycentre_);
          eulerian_compo_connex_ft_.echange_espace_virtuel(eulerian_compo_connex_ft_.ghost());
          eulerian_compo_connex_ghost_ft_.echange_espace_virtuel(eulerian_compo_connex_ghost_ft_.ghost());

          eulerian_compo_connex_ns_.data() = -1;
          eulerian_compo_connex_ns_.echange_espace_virtuel(eulerian_compo_connex_ns_.ghost());
          ref_ijk_ft_->eq_ns().redistribute_from_splitting_ft_elem(eulerian_compo_connex_ft_, eulerian_compo_connex_ns_);
          eulerian_compo_connex_ns_.echange_espace_virtuel(eulerian_compo_connex_ns_.ghost());

          eulerian_compo_connex_ghost_ns_.data() = -1;
          eulerian_compo_connex_ghost_ns_.echange_espace_virtuel(eulerian_compo_connex_ghost_ns_.ghost());
          ref_ijk_ft_->eq_ns().redistribute_from_splitting_ft_elem(eulerian_compo_connex_ghost_ft_, eulerian_compo_connex_ghost_ns_);
          eulerian_compo_connex_ghost_ns_.echange_espace_virtuel(eulerian_compo_connex_ghost_ns_.ghost());
        }
    }

  statistiques().end_count(cnt_compo_connex_bounding_box);
}

void IJK_Composantes_Connex::compute_compo_connex_from_interface()
{
  static Stat_Counter_Id cnt_compo_connex_interface = statistiques().new_counter(2, "Compo Connex - From interface");
  statistiques().begin_count(cnt_compo_connex_interface);

  if (compute_compo_fields_)
    {
      // interfaces_->calculer_volume_bulles(bubbles_volume_, bubbles_barycentre_);
      interfaces_->compute_bubbles_volume_and_barycentres(bubbles_volume_, bubbles_barycentre_, 1);

      fill_mixed_cell_compo();

      const Domaine_IJK& splitting = eulerian_compo_connex_from_interface_int_ns_.get_domaine();
      int neighours_i[6] = NEIGHBOURS_I;
      int neighours_j[6] = NEIGHBOURS_J;
      int neighours_k[6] = NEIGHBOURS_K;

      const int nx = eulerian_compo_connex_from_interface_int_ns_.ni();
      const int ny = eulerian_compo_connex_from_interface_int_ns_.nj();
      const int nz = eulerian_compo_connex_from_interface_int_ns_.nk();
      ArrOfInt elems_valid;
      for (int k=0; k < nz ; k++)
        for (int j=0; j< ny; j++)
          for (int i=0; i < nx; i++)
            if(eulerian_compo_connex_valid_compo_field_(i,j,k))
              elems_valid.append_array(splitting.convert_ijk_cell_to_packed(i,j,k));


      int elems_valid_size = elems_valid.size_array();
      while (elems_valid_size > 0)
        {
          ArrOfInt elems_valid_copy = elems_valid;
          elems_valid.reset();
          for (int elem=0; elem<elems_valid_copy.size_array(); elem++)
            {
              const Int3 num_elem_ijk = splitting.convert_packed_to_ijk_cell(elems_valid_copy[elem]);
              const int i = num_elem_ijk[DIRECTION_I];
              const int j = num_elem_ijk[DIRECTION_J];
              const int k = num_elem_ijk[DIRECTION_K];
              const int num_compo = eulerian_compo_connex_from_interface_int_ns_(i,j,k);
              const int num_compo_ghost = eulerian_compo_connex_from_interface_ghost_int_ns_(i,j,k);
              for (int l = 0; l < 6; l++)
                {
                  const int ii = neighours_i[l];
                  const int jj = neighours_j[l];
                  const int kk = neighours_k[l];
                  if((i + ii < 0 || j + jj < 0 || k + kk < 0) || (i + ii >= nx || j + jj >= ny || k + kk >= nz))
                    break;
                  const int num = eulerian_compo_connex_from_interface_int_ns_(i + ii,j + jj,k + kk);
                  const double indic_neighbour =  interfaces_->In()(i + ii,j + jj,k + kk);
                  if (num == -1 && indic_neighbour < VAPOUR_INDICATOR_TEST)
                    {
                      const int num_elem = splitting.convert_ijk_cell_to_packed(i + ii,j + jj,k + kk);
                      elems_valid.append_array(num_elem);
                      eulerian_compo_connex_from_interface_int_ns_(i + ii,j + jj,k + kk) = num_compo;
                      eulerian_compo_connex_from_interface_ghost_int_ns_(i + ii,j + jj,k + kk) = num_compo_ghost;
                    }
                }
            }
          elems_valid_size = elems_valid.size_array();
        }
      eulerian_compo_connex_from_interface_int_ns_.echange_espace_virtuel(eulerian_compo_connex_from_interface_int_ns_.ghost());
      eulerian_compo_connex_from_interface_ghost_int_ns_.echange_espace_virtuel(eulerian_compo_connex_from_interface_int_ns_.ghost());
    }

  statistiques().end_count(cnt_compo_connex_interface);
}

void IJK_Composantes_Connex::fill_mixed_cell_compo()
{
  static Stat_Counter_Id cnt_fill_mixed_cell_compo = statistiques().new_counter(3, "Fill Compo connex");
  statistiques().begin_count(cnt_fill_mixed_cell_compo);

  const Domaine_dis_base& mon_dom_dis = interfaces_->get_domaine_dis();
  const int nb_elem = mon_dom_dis.domaine().nb_elem();
  const Maillage_FT_IJK& maillage = interfaces_->maillage_ft_ijk();

  // Same splitting for the normal vector field
  const Domaine_IJK& splitting_ft = ref_ijk_ft_->get_domaine_ft();

  const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
  const ArrOfInt& index_elem = intersections.index_elem();
  eulerian_compo_connex_from_interface_ft_.data() = -1;
  eulerian_compo_connex_from_interface_ft_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ft_.ghost());
  eulerian_compo_connex_from_interface_ghost_ft_.data() = -1;
  eulerian_compo_connex_from_interface_ghost_ft_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ghost_ft_.ghost());
  // Loop on the elements
  const ArrOfInt& compo_facettes = maillage.compo_connexe_facettes();
  ArrOfInt compo_per_cell;
  ArrOfInt compo_ghost_per_cell;
  ArrOfInt count_compo_per_cell;
  ArrOfInt count_compo_ghost_per_cell;
  int counter, counter_ghost;
  FixedVector<IJK_Field_double *,2> compo_connex_non_ghost_ghost;
  compo_connex_non_ghost_ghost[0] = &eulerian_compo_connex_from_interface_ft_;
  compo_connex_non_ghost_ghost[1] = &eulerian_compo_connex_from_interface_ghost_ft_;
  const int nbulles_reelles = interfaces_->get_nb_bulles_reelles();
  int l;
  for (int elem = 0; elem < nb_elem; elem++)
    {
      int index = index_elem[elem];
      compo_per_cell.reset();
      compo_ghost_per_cell.reset();
      count_compo_per_cell.reset();
      count_compo_ghost_per_cell.reset();
      // Loop on the facets which cross the element
      while (index >= 0)
        {
          const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
          const int num_facette = data.numero_facette_;
          // const int n = mesh.nb_facettes();
          const int compo_facet = compo_facettes[num_facette];
          const int compo_size_array = compo_per_cell.size_array();
          const int compo_ghost_size_array = compo_ghost_per_cell.size_array();
          int compo_real;
          int compo_ghost;
          if (compo_facet < 0)
            {
              compo_real = decoder_numero_bulle(-compo_facet);
              compo_ghost = compo_real + nbulles_reelles;
            }
          else
            {
              compo_real = compo_facet;
              compo_ghost = compo_facet;
            }
          counter = 0;
          counter_ghost = 0;
          for (l=0; l<compo_size_array; l++)
            {
              if (compo_real == compo_per_cell(l))
                {
                  count_compo_per_cell(l) += 1;
                  break;
                }
              counter++;
            }
          for (l=0; l<compo_ghost_size_array; l++)
            {
              if (compo_ghost == compo_ghost_per_cell(l))
                {
                  count_compo_ghost_per_cell(l) += 1;
                  break;
                }
              counter_ghost++;
            }
          if (counter == compo_size_array)
            {

              compo_per_cell.append_array(compo_real);
              count_compo_per_cell.append_array(1);
            }
          if (counter_ghost == compo_ghost_size_array)
            {

              compo_ghost_per_cell.append_array(compo_ghost);
              count_compo_ghost_per_cell.append_array(1);
            }
          index = data.index_facette_suivante_;
        }
      const int n = compo_per_cell.size_array();
      const int n_ghost = compo_ghost_per_cell.size_array();
      if (n > 0)
        {
          std::vector<int> indices(n);
          for (int j=0; j<n; j++)
            indices[j]=j;
          std::sort(indices.begin(), indices.end(), [&count_compo_per_cell](int i, int j) {return count_compo_per_cell[i] < count_compo_per_cell[j];});
          const int max_compo_per_cell = compo_per_cell(indices[n-1]);
          const Int3 num_elem_ijk = splitting_ft.convert_packed_to_ijk_cell(elem);
          eulerian_compo_connex_from_interface_ft_(num_elem_ijk[DIRECTION_I],num_elem_ijk[DIRECTION_J],num_elem_ijk[DIRECTION_K]) = (double) max_compo_per_cell;
        }
      if (n_ghost > 0)
        {
          std::vector<int> indices(n);
          for (int j=0; j<n; j++)
            indices[j]=j;
          std::sort(indices.begin(), indices.end(), [&count_compo_ghost_per_cell](int i, int j) {return count_compo_ghost_per_cell[i] < count_compo_ghost_per_cell[j];});
          const int max_compo_ghost_per_cell = compo_ghost_per_cell(indices[n-1]);
          const Int3 num_elem_ijk = splitting_ft.convert_packed_to_ijk_cell(elem);
          eulerian_compo_connex_from_interface_ghost_ft_(num_elem_ijk[DIRECTION_I],num_elem_ijk[DIRECTION_J],num_elem_ijk[DIRECTION_K]) = (double) max_compo_ghost_per_cell;
        }
    }
  eulerian_compo_connex_from_interface_ft_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ft_.ghost());
  eulerian_compo_connex_from_interface_ghost_ft_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ghost_ft_.ghost());
  eulerian_compo_connex_from_interface_ns_.data() = -1;
  eulerian_compo_connex_from_interface_ghost_ns_.data() = -1;
  ref_ijk_ft_->eq_ns().redistribute_from_splitting_ft_elem(eulerian_compo_connex_from_interface_ft_, eulerian_compo_connex_from_interface_ns_);
  ref_ijk_ft_->eq_ns().redistribute_from_splitting_ft_elem(eulerian_compo_connex_from_interface_ghost_ft_, eulerian_compo_connex_from_interface_ghost_ns_);
  eulerian_compo_connex_from_interface_ns_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ns_.ghost());
  eulerian_compo_connex_from_interface_ghost_ns_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ghost_ns_.ghost());
  const int nx = eulerian_compo_connex_from_interface_int_ns_.ni();
  const int ny = eulerian_compo_connex_from_interface_int_ns_.nj();
  const int nz = eulerian_compo_connex_from_interface_int_ns_.nk();
  eulerian_compo_connex_valid_compo_field_.data() = 0;
  eulerian_compo_connex_from_interface_int_ns_.data() = -1;
  eulerian_compo_connex_from_interface_ghost_int_ns_.data() = -1;
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        {
          eulerian_compo_connex_from_interface_int_ns_(i,j,k) = (int) eulerian_compo_connex_from_interface_ns_(i,j,k);
          eulerian_compo_connex_from_interface_ghost_int_ns_(i,j,k) = (int) eulerian_compo_connex_from_interface_ghost_ns_(i,j,k);
          if (eulerian_compo_connex_from_interface_int_ns_(i,j,k) >= 0)
            eulerian_compo_connex_valid_compo_field_(i,j,k) = 1;
        }
  eulerian_compo_connex_from_interface_int_ns_.echange_espace_virtuel(eulerian_compo_connex_from_interface_int_ns_.ghost());
  eulerian_compo_connex_from_interface_ghost_int_ns_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ghost_int_ns_.ghost());
  eulerian_compo_connex_valid_compo_field_.echange_espace_virtuel(eulerian_compo_connex_valid_compo_field_.ghost());

  statistiques().end_count(cnt_fill_mixed_cell_compo);
}

void IJK_Composantes_Connex::compute_rising_velocities()
{
  static Stat_Counter_Id cnt_compute_fill_rising_vel = statistiques().new_counter(2, "Compute and fill rising velocity");
  statistiques().begin_count(cnt_compute_fill_rising_vel);

  if (compute_rising_velocities_)
    {
      const DoubleTab& bubbles_velocities_from_interface = interfaces_->get_bubble_velocities_from_interface();
      const DoubleTab& bubbles_velocities_from_barycentres = interfaces_->get_bubble_velocities_from_barycentres();
      int nb_bubbles = ref_ijk_ft_->get_interface().get_nb_bulles_reelles();
      rising_velocities_ = ArrOfDouble(nb_bubbles);
      rising_vectors_ = DoubleTab(nb_bubbles, 3);

      int use_bubbles_velocities_from_interface_tmp = use_bubbles_velocities_from_interface_;
      if (bubbles_velocities_from_interface.size() == 0)
        use_bubbles_velocities_from_interface_tmp = 0;

      int use_bubbles_velocities_from_barycentres_tmp = use_bubbles_velocities_from_barycentres_;
      if (ref_ijk_ft_->schema_temps_ijk().get_tstep() == 0)
        use_bubbles_velocities_from_barycentres_tmp = 0;

      compute_rising_velocity(ref_ijk_ft_->eq_ns().get_velocity(),
                              ref_ijk_ft_->get_interface(),
                              eulerian_compo_connex_from_interface_int_ns_,
                              ref_ijk_ft_->milieu_ijk().get_direction_gravite(),
                              rising_velocities_,
                              rising_vectors_,
                              liquid_velocity_,
                              bubbles_velocities_from_interface,
                              use_bubbles_velocities_from_interface_tmp,
                              bubbles_velocities_from_barycentres,
                              use_bubbles_velocities_from_barycentres_tmp);

      compute_rising_velocity_overall(ref_ijk_ft_->get_interface(),
                                      rising_vectors_,
                                      rising_velocities_,
                                      bubbles_volume_,
                                      rising_velocity_overall_,
                                      bubbles_velocities_from_interface,
                                      use_bubbles_velocities_from_interface_tmp,
                                      bubbles_velocities_from_barycentres,
                                      use_bubbles_velocities_from_barycentres_tmp);
      if (fill_rising_velocities_)
        {
          eulerian_rising_velocities_.data() = 0.;
          eulerian_rising_velocities_.echange_espace_virtuel(eulerian_rising_velocities_.ghost());
          fill_rising_velocity_int(eulerian_compo_connex_from_interface_int_ns_, rising_velocities_, eulerian_rising_velocities_);
        }
    }
  else
    Cerr << "Don't compute the ghost temperature field" << finl;

  statistiques().end_count(cnt_compute_fill_rising_vel);
}
