/****************************************************************************
* Copyright (c) 2024, CEA
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
// File      : IJK_One_Dimensional_Subproblems_Interfaces_Fields.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_One_Dimensional_Subproblems_Interfaces_Fields.h>
#include <IJK_FT.h>
#include <IJK_One_Dimensional_Subproblems.h>

Implemente_instanciable_sans_constructeur( IJK_One_Dimensional_Subproblems_Interfaces_Fields, "IJK_One_Dimensional_Subproblems_Interfaces_Fields", Objet_U ) ;

IJK_One_Dimensional_Subproblems_Interfaces_Fields::IJK_One_Dimensional_Subproblems_Interfaces_Fields()
{

}

int IJK_One_Dimensional_Subproblems_Interfaces_Fields::initialise(const IJK_Splitting& splitting,
                                                                  IJK_One_Dimensional_Subproblems& thermal_local_subproblems,
                                                                  const int& debug)
{
  thermal_local_subproblems_ = &thermal_local_subproblems;
  debug_ = debug;

  int nalloc = 0;
  tmp_ft_field_val_.allocate(ref_ijk_ft_->get_splitting_ft(), IJK_Splitting::ELEM, 1);
  nalloc += 1;
  tmp_ft_field_val_.data() = 0.;
  tmp_ft_field_val_.echange_espace_virtuel(tmp_ft_field_val_.ghost());

  tmp_field_val_.allocate(splitting, IJK_Splitting::ELEM, 1);
  nalloc += 1;
  tmp_field_val_.data() = 0.;
  tmp_field_val_.echange_espace_virtuel(tmp_field_val_.ghost());

  return nalloc;
}

void IJK_One_Dimensional_Subproblems_Interfaces_Fields::set_subproblems_interfaces_fields(const int& interface_field_type)
{
  /*
   * interface_field_type_ = 0 - Do Nothing
   * interface_field_type_ = 1 - GFM
   * interface_field_type_ = 2 - Subres
   */
  interface_field_type_ = interface_field_type;
}

void IJK_One_Dimensional_Subproblems_Interfaces_Fields::copy_previous_interface_state()
{
  if (interface_field_type_)
    {
      /*
       * TODO : PB interface at time (n+1) only
       * Need a copy
       * There's a problem with Intersections_Elem_Facettes_Data which is a pointer
       * Need a usr-defined destructor for the current class...
       */
      // const IJK_Interfaces& interfaces = ref_ijk_ft_->itfce();
      // const Maillage_FT_IJK& maillage = interfaces.maillage_ft_ijk();
      // const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
      // const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
      // const IntTab& facettes = maillage.facettes();
      // intersections_ = intersections;
      // surface_facettes_ = surface_facettes;
      // facettes_ = facettes;
    }
}

void IJK_One_Dimensional_Subproblems_Interfaces_Fields::reset_flags()
{
  updated_interfacial_heat_flux_ = false;
  updated_interfacial_heat_flux_sol_ = false;
  updated_velocity_magnitude_ = false;
  updated_temperature_ = false;
  updated_temperature_sol_ = false;
  updated_temperature_gradient_ = false;
  updated_temperature_gradient_sol_ = false;
  counter_call_ = 0;
}

Sortie& IJK_One_Dimensional_Subproblems_Interfaces_Fields::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  return os;
}

Entree& IJK_One_Dimensional_Subproblems_Interfaces_Fields::readOn( Entree& is )
{
  Objet_U::readOn( is );
  return is;
}

void IJK_One_Dimensional_Subproblems_Interfaces_Fields::associer(const IJK_FT_double& ijk_ft)
{
  ref_ijk_ft_ = ijk_ft;
  liste_post_instantanes_ = ijk_ft.get_post().get_liste_post_instantanes();
}

void IJK_One_Dimensional_Subproblems_Interfaces_Fields::posttraiter_tous_champs(Motcles& liste) const
{
  liste.add("FT_INTERFACIAL_HEAT_FLUX");
  liste.add("FT_INTERFACIAL_VELOCITY_MAGNITUDE");
  liste.add("FT_INTERFACIAL_TEMPERATURE");
}


int IJK_One_Dimensional_Subproblems_Interfaces_Fields::posttraiter_champs_instantanes(const Motcles& liste_post_instantanes,
                                                                                      const char *lata_name,
                                                                                      const int lata_step)
{
  int n = 0;
  if (interface_field_type_)
    {
      const Nom interf_heat_flux("FT_INTERFACIAL_HEAT_FLUX_DENSITY");
      if (liste_post_instantanes.contient_(interf_heat_flux))
        {
          bool is_updated = retrieve_interfacial_surface_quantities(interf_heat_flux);
          if (is_updated)
            n++, dumplata_ft_field(lata_name, "INTERFACES", interf_heat_flux, "SOM",  interfacial_heat_flux_sol_, lata_step);
        }

      const Nom interf_heat_flux_interp("FT_INTERFACIAL_HEAT_FLUX_DENSITY_INTERP");
      if (liste_post_instantanes.contient_(interf_heat_flux_interp))
        {
          bool is_updated = retrieve_interfacial_surface_quantities(interf_heat_flux_interp);
          if (is_updated)
            n++, dumplata_ft_field(lata_name, "INTERFACES", interf_heat_flux_interp, "SOM",  interfacial_heat_flux_, lata_step);
        }

      const Nom interf_velocity_magnitude("FT_INTERFACIAL_VELOCITY_MAGNITUDE");
      if (liste_post_instantanes.contient_(interf_velocity_magnitude))
        {
          bool is_updated = retrieve_interfacial_surface_quantities(interf_velocity_magnitude);
          if (is_updated)
            n++, dumplata_ft_field(lata_name, "INTERFACES", interf_velocity_magnitude, "SOM",  velocity_magnitude_, lata_step);
        }

      const Nom interf_temperature("FT_INTERFACIAL_TEMPERATURE");
      if (liste_post_instantanes.contient_(interf_temperature))
        {
          bool is_updated = retrieve_interfacial_surface_quantities(interf_temperature);
          if (is_updated)
            n++, dumplata_ft_field(lata_name, "INTERFACES", interf_temperature, "SOM",  temperature_sol_, lata_step);
        }

      const Nom interf_temperature_gradient("FT_INTERFACIAL_TEMPERATURE_GRADIENT");
      if (liste_post_instantanes.contient_(interf_temperature_gradient))
        {
          bool is_updated = retrieve_interfacial_surface_quantities(interf_temperature_gradient);
          if (is_updated)
            n++, dumplata_ft_field(lata_name, "INTERFACES", interf_temperature_gradient, "SOM",  temperature_gradient_sol_, lata_step);
        }

      const Nom interf_temperature_interp("FT_INTERFACIAL_TEMPERATURE_INTERP");
      if (liste_post_instantanes.contient_(interf_temperature_interp))
        {
          bool is_updated = retrieve_interfacial_surface_quantities(interf_temperature_interp);
          if (is_updated)
            n++, dumplata_ft_field(lata_name, "INTERFACES", interf_temperature_interp, "SOM",  temperature_, lata_step);
        }

      const Nom interf_temperature_gradient_interp("FT_INTERFACIAL_TEMPERATURE_GRADIENT_INTERP");
      if (liste_post_instantanes.contient_(interf_temperature_gradient_interp))
        {
          bool is_updated = retrieve_interfacial_surface_quantities(interf_temperature_gradient_interp);
          if (is_updated)
            n++, dumplata_ft_field(lata_name, "INTERFACES", interf_temperature_gradient_interp, "SOM",  temperature_gradient_, lata_step);
        }

      reset_flags();
    }
  return n;
}

bool IJK_One_Dimensional_Subproblems_Interfaces_Fields::retrieve_interfacial_surface_quantities(const Nom& surface_post_name)
{
  if (surface_post_name == Nom("FT_INTERFACIAL_HEAT_FLUX_DENSITY"))
    return retrieve_interfacial_surface_quantity(interfacial_heat_flux_sol_, updated_interfacial_heat_flux_sol_, 0);
  else if (surface_post_name == Nom("FT_INTERFACIAL_VELOCITY_MAGNITUDE"))
    return retrieve_interfacial_surface_quantity(velocity_magnitude_, updated_velocity_magnitude_, 1);
  else if (surface_post_name == Nom("FT_INTERFACIAL_TEMPERATURE"))
    return retrieve_interfacial_surface_quantity(temperature_sol_, updated_temperature_sol_, 2);
  else if (surface_post_name == Nom("FT_INTERFACIAL_TEMPERATURE_GRADIENT"))
    return retrieve_interfacial_surface_quantity(temperature_gradient_sol_, updated_temperature_gradient_sol_, 3);
  else if (surface_post_name == Nom("FT_INTERFACIAL_TEMPERATURE_INTERP"))
    return retrieve_interfacial_surface_quantity(temperature_, updated_temperature_, 4);
  else if (surface_post_name == Nom("FT_INTERFACIAL_TEMPERATURE_GRADIENT_INTERP"))
    return retrieve_interfacial_surface_quantity(temperature_gradient_, updated_temperature_gradient_, 5);
  else if (surface_post_name == Nom("FT_INTERFACIAL_HEAT_FLUX_DENSITY_INTERP"))
    return retrieve_interfacial_surface_quantity(interfacial_heat_flux_, updated_interfacial_heat_flux_, 6);
  else
    return false;
}

bool IJK_One_Dimensional_Subproblems_Interfaces_Fields::retrieve_interfacial_surface_quantity(ArrOfDouble& surface_quantity, bool& has_been_updated, const int& val_index)
{
  int counter_val_elem = 0;
  int counter_val_interf = 0;
  static const double invalid_value = INVALID_TEST;
  if (!has_been_updated)
    {
      vertices_surface_.reset();
      eulerian_values_.reset();
      surface_quantity.reset();
      tmp_field_val_.data() = invalid_value * 1.1;
      tmp_ft_field_val_.data() = invalid_value * 1.1;
      if (!counter_call_)
        elem_crossed_.reset();
      const int dim = 3;
      const IJK_Interfaces& interfaces = ref_ijk_ft_->itfce();
      const IJK_Splitting& splitting = ref_ijk_ft_->get_splitting_ft();

      const Domaine_dis_base& mon_dom_dis = interfaces.get_domaine_dis().valeur();
      const Maillage_FT_IJK& maillage = interfaces.maillage_ft_ijk();
      const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
      const ArrOfInt& index_elem = intersections.index_elem();
      const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
      const ArrOfInt& sommet_elem = maillage.sommet_elem();
      const IntTab& facettes = maillage.facettes();
      const DoubleTab& sommets = maillage.sommets();

      /*
       * TODO : PB interface at time (n+1) only
       */
      // const Intersections_Elem_Facettes& intersections= intersections_;
      // const ArrOfInt& index_elem = intersections.index_elem();
      // const ArrOfDouble& surface_facettes = surface_facettes_;
      // const IntTab& facettes = facettes_;
      // const ArrOfInt& sommets_elems = maillage.sommet_elem();
      const int nb_vertices = sommets.dimension(0);
      const int nb_vertices_elem = sommet_elem.size_array();
      if (debug_)
        {
          Cerr << "Nb vertices: " << nb_vertices << finl;
          Cerr << "Nb vertices elem: " << nb_vertices_elem << finl;
        }
      vertices_surface_.resize(nb_vertices);
      surface_quantity.resize(nb_vertices);
      vertices_surface_ *= 0.;
      surface_quantity *= 0.;

      int index_i, index_j, index_k;
      const int nb_subproblems = thermal_local_subproblems_->get_effective_subproblems_counter();
      for (int m = 0; m < nb_subproblems; m++)
        {
          const double eulerian_value = thermal_local_subproblems_->get_thermal_subproblem_value_at_ijk_index(m, index_i, index_j, index_k, val_index);
          tmp_field_val_(index_i, index_j, index_k) = eulerian_value;
          eulerian_values_.append_array(eulerian_value);
          if (!counter_call_)
            {
              const int elem = splitting.convert_ijk_cell_to_packed(index_i, index_j, index_k);
              elem_crossed_.append_array(elem);
            }
        }

      ref_ijk_ft_->redistribute_to_splitting_ft_elem(tmp_field_val_, tmp_ft_field_val_);
      tmp_ft_field_val_.echange_espace_virtuel(tmp_ft_field_val_.ghost());

      const int nb_elem = mon_dom_dis.domaine().nb_elem();

      // const int nb_elem = (int) elem_crossed_.size_array();
      for (int elem = 0; elem < nb_elem; elem++)
        {
          // const int local_elem = elem_crossed_[elem];
          // int index = index_elem[local_elem];
          // const double eulerian_value = eulerian_values_[elem];
          int index = index_elem[elem];
          if (index >= 0)
            {
              double eulerian_value = 0.;
              const Int3 num_elem_ijk = splitting.convert_packed_to_ijk_cell(elem);
              const int& i = num_elem_ijk[DIRECTION_I];
              const int& j = num_elem_ijk[DIRECTION_J];
              const int& k = num_elem_ijk[DIRECTION_K];
              get_surrounding_value(i, j, k, eulerian_value);
              // Loop on the facets which cross the element
              while (index >= 0)
                {
                  const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
                  const int num_facette = data.numero_facette_;
                  const double surface = data.fraction_surface_intersection_ * surface_facettes[num_facette];
                  if (abs(surface) > MIN_SURFACE)
                    for (int l = 0; l < dim; l++)
                      {
                        const int som = facettes(num_facette, l);
                        vertices_surface_[som] += surface;
                        const double surface_value = eulerian_value;
                        surface_quantity[som] += (surface_value * surface);
                        counter_val_interf++;
                      }
                  else if (debug_)
                    Cerr << "Small surface intersection detected: " << surface << finl;
                  index = data.index_facette_suivante_;
                }
              counter_val_elem++;
            }
        }
      averaged_by_vertex_surface(surface_quantity);
      has_been_updated = true;
      counter_call_ += 1;
    }
  if (debug_)
    {
      Cerr << "Elems contributing to interfacial post-processing: " << counter_val_elem << finl;
      Cerr << "Facets contributing to interfacial post-processing: " << counter_val_interf << finl;
    }
  return true && has_been_updated;
}

void IJK_One_Dimensional_Subproblems_Interfaces_Fields::averaged_by_vertex_surface(ArrOfDouble& surface_quantity)
{
  const int size_array = (int) vertices_surface_.size_array();
  for (int m=0; m<size_array; m++)
    {
      const double surface_quantity_tmp = surface_quantity[m];
      const double surface_tot_som = vertices_surface_[m];
      if (abs(surface_tot_som) > MIN_SURFACE)
        surface_quantity[m] = surface_quantity_tmp / surface_tot_som;
      else if (debug_)
        {
          Cerr << "Small interface portion detected" << finl;
          Cerr << "Surface value of: " << surface_tot_som << finl;
        }
    }
}

void IJK_One_Dimensional_Subproblems_Interfaces_Fields::get_surrounding_value(const int& i, const int& j, const int& k, double& eulerian_value)
{
  static const double invalid_value = INVALID_TEST;
  eulerian_value = tmp_ft_field_val_(i,j,k);
  if (eulerian_value < invalid_value)
    {
      // double max_val = 0.;
      double avg_val = 0.;
      int nb_avg = 0;
      for (int c=0; c<3; c++)
        for (int l=-1; l<=1; l++)
          {
            const int ii = select(c, l, 0, 0);
            const int jj = select(c, 0, l, 0);
            const int kk = select(c, 0, 0, l);
            const double local_val = tmp_ft_field_val_(i+ii, j+jj, k+kk);
            if (local_val > INVALID_TEST)
              {
                // max_val = (abs(local_val) > abs(max_val) ? local_val : max_val);
                avg_val += local_val;
                nb_avg++;
              }
          }
      if (nb_avg)
        avg_val /= nb_avg;
      // tmp_ft_field_val_(i, j, k) = avg_val;
      eulerian_value = avg_val;
    }
}

