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
// File      : IJK_Ghost_Fluid_Fields.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Ghost_Fluid_Fields.h>
#include <Probleme_FTD_IJK_base.h>
#include <IJK_Navier_Stokes_tools.h>
#include <IJK_Ghost_Fluid_tools.h>
#include <IJK_Bubble_tools.h>


Implemente_instanciable_sans_constructeur( IJK_Ghost_Fluid_Fields, "IJK_Ghost_Fluid_Fields", Objet_U ) ;

IJK_Ghost_Fluid_Fields::IJK_Ghost_Fluid_Fields()
{

}

Sortie& IJK_Ghost_Fluid_Fields::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  return os;
}

Entree& IJK_Ghost_Fluid_Fields::readOn( Entree& is )
{
  Objet_U::readOn( is );
  return is;
}

void IJK_Ghost_Fluid_Fields::associer(const Probleme_FTD_IJK_base& ijk_ft)
{
  ref_ijk_ft_ = ijk_ft;
}

void IJK_Ghost_Fluid_Fields::initialize(const Domaine_IJK& splitting)
{
  if (compute_distance_)
    {
      /*
       * TODO: Move to IJK_Interfaces
       */
      // Laplacian(d) necessitates 2 ghost cells like temperature
      const int n_iter_base = use_n_iter_distance_ ? n_iter_distance_ : nb_cells_gfm_parallel_calls_;
      const int nb_ghost_parallel_calls = 2;
      const int dist_ghost = avoid_gfm_parallel_calls_ ? n_iter_base + nb_ghost_parallel_calls : 2;

      eulerian_distance_ft_.allocate(ref_ijk_ft_->get_domaine_ft(), Domaine_IJK::ELEM, dist_ghost);

      tmp_old_dist_val_ = eulerian_distance_ft_;
      tmp_new_dist_val_ = eulerian_distance_ft_;

      const int dist_tmp_ghost = avoid_gfm_parallel_calls_ ? n_iter_base + nb_ghost_parallel_calls : 0;
      tmp_interf_cells_.allocate(ref_ijk_ft_->get_domaine_ft(), Domaine_IJK::ELEM, dist_tmp_ghost);
      tmp_propagated_cells_.allocate(ref_ijk_ft_->get_domaine_ft(), Domaine_IJK::ELEM, dist_tmp_ghost);

      // grad(d) necessitates 1 ghost cell ?
      const int normal_ghost = avoid_gfm_parallel_calls_ ? n_iter_base + nb_ghost_parallel_calls : 1;
      allocate_cell_vector(eulerian_normal_vectors_ft_, ref_ijk_ft_->get_domaine_ft(), normal_ghost);

      tmp_old_vector_val_ = eulerian_normal_vectors_ft_;
      tmp_new_vector_val_ = eulerian_normal_vectors_ft_;
      // allocate_velocity(eulerian_normal_vectors_, ref_ijk_ft_->get_domaine_ft(), 1);

      allocate_cell_vector(eulerian_facets_barycentre_ft_, ref_ijk_ft_->get_domaine_ft(), 0);

      eulerian_distance_ft_.echange_espace_virtuel(eulerian_distance_ft_.ghost());
      eulerian_normal_vectors_ft_.echange_espace_virtuel();
      eulerian_facets_barycentre_ft_.echange_espace_virtuel();
      /*
       * TODO: This is already calculated in IJK_Interfaces
       * Keep it for now and clean later
       */
      eulerian_distance_ns_.allocate(splitting, Domaine_IJK::ELEM, 2);
      allocate_cell_vector(eulerian_normal_vectors_ns_, splitting, 1);
      allocate_cell_vector(eulerian_facets_barycentre_ns_, splitting, 0);
      eulerian_distance_ns_.echange_espace_virtuel(eulerian_distance_ns_.ghost());
      eulerian_normal_vectors_ns_.echange_espace_virtuel();
      eulerian_facets_barycentre_ns_.echange_espace_virtuel();
      allocate_cell_vector(eulerian_normal_vectors_ns_normed_, splitting, 1);
      eulerian_normal_vectors_ns_normed_.echange_espace_virtuel();

      if (avoid_gfm_parallel_calls_)
        {
          assert(eulerian_distance_ft_.ghost()==tmp_old_dist_val_.ghost());
          assert(eulerian_distance_ft_.ghost()==tmp_new_dist_val_.ghost());
          assert(eulerian_normal_vectors_ft_[0].ghost()==tmp_old_vector_val_[0].ghost());
          assert(eulerian_normal_vectors_ft_[0].ghost()==tmp_new_vector_val_[0].ghost());
        }
    }
  if (compute_curvature_)
    {
      // Laplacian(d) necessitates 0 ghost cells like div_lambda_grad_T
      // but if calculated using the neighbours maybe 1
      // const int curvature_ghost = avoid_gfm_parallel_calls_ ? n_iter_distance_ : 1;
      const int curvature_ghost = avoid_gfm_parallel_calls_ ? 1 : 1;
      eulerian_curvature_ft_.allocate(ref_ijk_ft_->get_domaine_ft(), Domaine_IJK::ELEM, curvature_ghost);

      tmp_old_curv_val_ = eulerian_curvature_ft_;
      tmp_new_curv_val_ = eulerian_curvature_ft_;

      eulerian_curvature_ft_.echange_espace_virtuel(eulerian_curvature_ft_.ghost());
      // Only calculated in the mixed cells ghost_cells = 0
      eulerian_interfacial_area_ft_.allocate(ref_ijk_ft_->get_domaine_ft(), Domaine_IJK::ELEM, 0);
      eulerian_interfacial_area_ft_.echange_espace_virtuel(eulerian_interfacial_area_ft_.ghost());
      /*
       * TODO: This is already calculated in IJK_Interfaces
       * Keep it for now and clean later
       */
      eulerian_curvature_ns_.allocate(splitting, Domaine_IJK::ELEM, 1);
      eulerian_curvature_ns_.echange_espace_virtuel(eulerian_curvature_ns_.ghost());
      eulerian_interfacial_area_ns_.allocate(splitting, Domaine_IJK::ELEM, 0);
      eulerian_interfacial_area_ns_.echange_espace_virtuel(eulerian_interfacial_area_ns_.ghost());
    }
}

void IJK_Ghost_Fluid_Fields::Fill_postprocessable_fields(std::vector<FieldInfo_t>& chps)
{
  std::vector<FieldInfo_t> c =
  {
    // Name     /     Localisation (elem, face, ...) /    Nature (scalare, vector)   /  Located on interface?
    { "TEMPERATURE", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "TEMPERATURE_ADIMENSIONNELLE_THETA", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "DIV_LAMBDA_GRAD_T_VOLUME", Entity::ELEMENT, Nature_du_champ::scalaire, false },

  };
  chps.insert(chps.end(), c.begin(), c.end());
}



void IJK_Ghost_Fluid_Fields::get_noms_champs_postraitables(Noms& noms,Option opt) const
{
  for (const auto& n : champs_compris_.liste_noms_compris())
    noms.add(n);
  for (const auto& n : champs_compris_.liste_noms_compris_vectoriel())
    noms.add(n);
}

bool IJK_Ghost_Fluid_Fields::has_champ(const Motcle& nom) const
{
  if (champs_compris_.liste_noms_compris().contient_(nom))
    return champs_compris_.has_champ(nom);
  else if (champs_compris_.liste_noms_compris_vectoriel().contient_(nom))
    return champs_compris_.has_champ_vectoriel(nom);
  else
    return false;
}

const IJK_Field_double& IJK_Ghost_Fluid_Fields::get_IJK_field(const Motcle& nom)
{
  if (champs_compris_.liste_noms_compris().contient_(nom))
    return champs_compris_.get_champ(nom);
  else
    {
      Cerr << "ERROR in IJK_Ghost_Fluid_Fields::get_IJK_field : " << finl;
      Cerr << "Requested field '" << nom << "' is not recognized by IJK_Ghost_Fluid_Fields::get_IJK_field()." << finl;
      throw;
    }
}

const IJK_Field_vector3_double& IJK_Ghost_Fluid_Fields::get_IJK_field_vector(const Motcle& nom)
{
  if (champs_compris_.liste_noms_compris_vectoriel().contient_(nom))
    return champs_compris_.get_champ_vectoriel(nom);
  else
    {
      Cerr << "ERROR in IJK_Ghost_Fluid_Fields::get_IJK_field_vector : " << finl;
      Cerr << "Requested field '" << nom << "' is not recognized by IJK_Ghost_Fluid_Fields::get_IJK_field_vector()." << finl;
      throw;
    }
}
static void enforce_zero_value_eulerian_field(IJK_Field_double& eulerian_field)
{
  const int nx = eulerian_field.ni();
  const int ny = eulerian_field.nj();
  const int nz = eulerian_field.nk();
  static const double invalid_distance_value = -1.e30;
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        if (eulerian_field(i,j,k) < invalid_distance_value)
          eulerian_field(i,j,k) = 0.;
}

static void enforce_max_value_eulerian_field(IJK_Field_double& eulerian_field)
{
  double eulerian_field_max = -1.e20;
  const int nx = eulerian_field.ni();
  const int ny = eulerian_field.nj();
  const int nz = eulerian_field.nk();
  static const double invalid_distance_value = -1.e30;
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        eulerian_field_max = std::max(eulerian_field_max, eulerian_field(i,j,k));
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        if (eulerian_field(i,j,k) < invalid_distance_value)
          eulerian_field(i,j,k) = eulerian_field_max;
}

//static void enforce_min_value_eulerian_field(IJK_Field_double& eulerian_field)
//{
//  double eulerian_field_min = 1.e20;
//  const int nx = eulerian_field.ni();
//  const int ny = eulerian_field.nj();
//  const int nz = eulerian_field.nk();
//  static const double invalid_distance_value = -1.e30;
//  for (int k=0; k < nz ; k++)
//    for (int j=0; j< ny; j++)
//      for (int i=0; i < nx; i++)
//        eulerian_field_min = std::min(eulerian_field_min, eulerian_field(i,j,k));
//  for (int k=0; k < nz ; k++)
//    for (int j=0; j< ny; j++)
//      for (int i=0; i < nx; i++)
//        if (eulerian_field(i,j,k) < invalid_distance_value)
//          eulerian_field(i,j,k) = eulerian_field_min;
//}

void IJK_Ghost_Fluid_Fields::compute_eulerian_distance()
{
  if (!has_computed_distance_)
    {
      if (compute_distance_)
        {
          Navier_Stokes_FTD_IJK& ns = ref_ijk_ft_->eq_ns();

          // TODO: Do we need to perform an echange_virtuel with interfaces ?
          compute_eulerian_normal_distance_facet_barycentre_field(ref_ijk_ft_->get_interface(),
                                                                  eulerian_distance_ft_,
                                                                  eulerian_normal_vectors_ft_,
                                                                  eulerian_facets_barycentre_ft_,
                                                                  tmp_old_vector_val_,
                                                                  tmp_new_vector_val_,
                                                                  tmp_old_dist_val_,
                                                                  tmp_new_dist_val_,
                                                                  tmp_interf_cells_,
                                                                  tmp_propagated_cells_,
                                                                  interf_cells_indices_,
                                                                  gfm_first_cells_indices_,
                                                                  propagated_cells_indices_,
                                                                  n_iter_distance_,
                                                                  avoid_gfm_parallel_calls_);
          eulerian_distance_ft_.echange_espace_virtuel(eulerian_distance_ft_.ghost());
          eulerian_distance_ns_.data() = 0.;
          eulerian_distance_ns_.echange_espace_virtuel(eulerian_distance_ns_.ghost());
          ns.redistribute_from_splitting_ft_elem(eulerian_distance_ft_, eulerian_distance_ns_);
          eulerian_distance_ns_.echange_espace_virtuel(eulerian_distance_ns_.ghost());
          for(int dir=0; dir<3; dir++)
            {
              ns.redistribute_from_splitting_ft_elem(eulerian_normal_vectors_ft_[dir], eulerian_normal_vectors_ns_[dir]);
              ns.redistribute_from_splitting_ft_elem(eulerian_facets_barycentre_ft_[dir], eulerian_facets_barycentre_ns_[dir]);
            }
          eulerian_normal_vectors_ns_normed_[0].data() = 0.;
          eulerian_normal_vectors_ns_normed_[1].data() = 0.;
          eulerian_normal_vectors_ns_normed_[2].data() = 0.;
          const int nx = eulerian_normal_vectors_ns_normed_[0].ni();
          const int ny = eulerian_normal_vectors_ns_normed_[0].nj();
          const int nz = eulerian_normal_vectors_ns_normed_[0].nk();
          for (int k=0; k < nz ; k++)
            for (int j=0; j< ny; j++)
              for (int i=0; i < nx; i++)
                {
                  double norm_x = eulerian_normal_vectors_ns_[0](i,j,k);
                  double norm_y = eulerian_normal_vectors_ns_[1](i,j,k);
                  double norm_z = eulerian_normal_vectors_ns_[2](i,j,k);
                  norm_x *= norm_x;
                  norm_y *= norm_y;
                  norm_z *= norm_z;
                  const double norm = norm_x + norm_y + norm_z;
                  if (norm > 0)
                    {
                      eulerian_normal_vectors_ns_normed_[0](i,j,k) = eulerian_normal_vectors_ns_[0](i,j,k) / sqrt(norm);
                      eulerian_normal_vectors_ns_normed_[1](i,j,k) = eulerian_normal_vectors_ns_[1](i,j,k) / sqrt(norm);
                      eulerian_normal_vectors_ns_normed_[2](i,j,k) = eulerian_normal_vectors_ns_[2](i,j,k) / sqrt(norm);
                    }
                }
          eulerian_normal_vectors_ns_normed_.echange_espace_virtuel();
          // has_computed_distance_ = true;
        }
      else
        Cerr << "Don't compute the eulerian distance field" << finl;
    }
}

void IJK_Ghost_Fluid_Fields::enforce_distance_curvature_values_for_post_processings()
{
  /*
   * For post-processings purposes
   */
  enforce_zero_value_eulerian_distance();
  enforce_max_value_eulerian_curvature();
  enforce_max_value_eulerian_field(eulerian_interfacial_area_ft_);
}

void IJK_Ghost_Fluid_Fields::enforce_zero_value_eulerian_distance()
{
  if (compute_distance_)
    {
      enforce_zero_value_eulerian_field(eulerian_distance_ft_);
      enforce_zero_value_eulerian_field(eulerian_distance_ns_);
    }
  else
    Cerr << "Eulerian distance has not been computed" << finl;
}

void IJK_Ghost_Fluid_Fields::compute_eulerian_curvature()
{
  if (compute_curvature_)
    {
      /*
       * Laplacian operator may not work properly with FT_field ?
       */
      eulerian_distance_ft_.echange_espace_virtuel(eulerian_distance_ft_.ghost());
      compute_eulerian_curvature_field_from_distance_field(eulerian_distance_ft_,
                                                           eulerian_curvature_ft_,
                                                           boundary_flux_kmin_,
                                                           boundary_flux_kmax_);
      eulerian_curvature_ft_.echange_espace_virtuel(eulerian_curvature_ft_.ghost());
      ref_ijk_ft_->eq_ns().redistribute_from_splitting_ft_elem(eulerian_curvature_ft_, eulerian_curvature_ns_);
    }
  else
    Cerr << "Don't compute the eulerian curvature field" << finl;
}

void IJK_Ghost_Fluid_Fields::compute_eulerian_curvature_from_interface()
{
  if (!has_computed_curvature_)
    {
      if (compute_curvature_)
        {
          eulerian_interfacial_area_ft_.echange_espace_virtuel(eulerian_interfacial_area_ft_.ghost());
          eulerian_normal_vectors_ft_.echange_espace_virtuel();
          int nb_groups = ref_ijk_ft_->get_interface().nb_groups();
          // Boucle debute a -1 pour faire l'indicatrice globale.
          // S'il n'y a pas de groupes de bulles (monophasique ou monodisperse), on passe exactement une fois dans la boucle
          if (nb_groups == 1)
            nb_groups = 0; // Quand il n'y a qu'un groupe, on ne posttraite pas les choses pour ce groupe unique puisque c'est identique au cas global
          for (int igroup = -1; igroup < nb_groups; igroup++)
            {
              compute_eulerian_curvature_field_from_interface(eulerian_normal_vectors_ft_,
                                                              ref_ijk_ft_->get_interface(),
                                                              eulerian_interfacial_area_ft_,
                                                              eulerian_curvature_ft_,
                                                              tmp_old_curv_val_,
                                                              tmp_new_curv_val_,
                                                              n_iter_distance_,
                                                              igroup);
            }
          eulerian_interfacial_area_ft_.echange_espace_virtuel(eulerian_interfacial_area_ft_.ghost());
          eulerian_curvature_ft_.echange_espace_virtuel(eulerian_curvature_ft_.ghost());
          ref_ijk_ft_->eq_ns().redistribute_from_splitting_ft_elem(eulerian_interfacial_area_ft_, eulerian_interfacial_area_ns_);
          ref_ijk_ft_->eq_ns().redistribute_from_splitting_ft_elem(eulerian_curvature_ft_, eulerian_curvature_ns_);
          // has_computed_curvature_ = true;
        }
      else
        Cerr << "Don't compute the eulerian curvature field" << finl;
    }
}

void IJK_Ghost_Fluid_Fields::enforce_zero_value_eulerian_curvature()
{
  if (compute_curvature_)
    enforce_zero_value_eulerian_field(eulerian_curvature_ft_);
  else
    Cerr << "Eulerian curvature has not been computed" << finl;
}

void IJK_Ghost_Fluid_Fields::enforce_max_value_eulerian_curvature()
{
  if (compute_curvature_)
    {
      enforce_max_value_eulerian_field(eulerian_curvature_ft_);
      enforce_max_value_eulerian_field(eulerian_curvature_ns_);
    }
  else
    Cerr << "Eulerian curvature has not been computed" << finl;
}
