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

#include <Cut_cell_FT_Disc.h>
#include <Maillage_FT_Disc.h>
#include <IJK_Field.h>
#include <Array_tools.h>
#include <IJK_Splitting.h>

Cut_cell_FT_Disc::Cut_cell_FT_Disc(IJK_Interfaces& interfaces, IJK_Splitting& splitting) :
  interfaces_(interfaces),
  splitting_(splitting),
  ghost_size_(4),
  processed_count_(0),
  n_loc_(0),
  n_tot_(0)
{
  linear_index_.associer_paresseux(*this);
  permutation_.associer_paresseux(*this);
  processed_.associer_paresseux(*this);

  coord_.associer_persistant(*this);

  index_sorted_by_k_.associer_paresseux(*this);
  k_value_index_.resize(0, 1);

  index_sorted_by_indicatrice_.associer_paresseux(*this);

  index_sorted_by_statut_diphasique_.associer_paresseux(*this);
  statut_diphasique_value_index_.resize(0, 1);
}

void Cut_cell_FT_Disc::initialise()
{
  /* Initialisation sans indicatrice (suppose aucune cellule diphasique) */

  // Initialisation des tableaux
  n_loc_ = 0;
  n_tot_ = 0;

  permutation_.resize(0, permutation_.dimension(1));

  processed_.resize(0, processed_.dimension(1));

  resize_data(0);

  // Initialisation du champ IJK_Field utilise pour le post-traitement
  write_buffer_.allocate(splitting_, IJK_Splitting::ELEM, 1);

  // Initialisation du champ IJK_Field des indices diphasiques
  indice_diphasique_.allocate(splitting_, IJK_Splitting::ELEM, ghost_size_);

  // Remplissage de certains champs a partir des valeurs sur la structure IJK_Field
  set_coord();

  // Communication des elements virtuels
  initialise_schema_comm();

  n_tot_ = initialise_communications();
  resize_data(n_tot_);

  echange_espace_virtuel();

  // Calcul de l'indice lineaire pour les elements virtuels
  linear_index_.resize(n_tot_, linear_index_.dimension(1));
  compute_virtual_linear_index();

  // Remplissage du tableau des indices diphasiques
  remplir_indice_diphasique();

  // Creation des tableaux des indices tries par k
  update_index_sorted_by_k();

  // Verifications
  assert(verifier_coherence_coord_linear_index());
  assert(verifier_taille_tableaux());
}


void Cut_cell_FT_Disc::initialise(const IJK_Field_double& old_indicatrice, const IJK_Field_double& next_indicatrice)
{
  // Initialisation des tableaux
  n_loc_ = initialise_linear_index(old_indicatrice, next_indicatrice);
  n_tot_ = n_loc_;

  permutation_.resize(n_loc_, permutation_.dimension(1));
  initialise_permutation();

  processed_.resize(n_loc_, processed_.dimension(1));
  initialise_processed();

  resize_data(n_loc_);

  // Initialisation du champ IJK_Field utilise pour le post-traitement
  write_buffer_.allocate(splitting_, IJK_Splitting::ELEM, 1);

  // Initialisation du champ IJK_Field des indices diphasiques
  indice_diphasique_.allocate(splitting_, IJK_Splitting::ELEM, ghost_size_);

  // Remplissage de certains champs a partir des valeurs sur la structure IJK_Field
  set_coord();

  // Communication des elements virtuels
  initialise_schema_comm();

  n_tot_ = initialise_communications();
  resize_data(n_tot_);

  echange_espace_virtuel();

  // Calcul de l'indice lineaire pour les elements virtuels
  linear_index_.resize(n_tot_, linear_index_.dimension(1));
  compute_virtual_linear_index();

  // Remplissage du tableau des indices diphasiques
  remplir_indice_diphasique();

  // Creation des tableaux des indices tries par k
  update_index_sorted_by_k();

  // Creation des tableaux des indices tries par l'indicatrice
  update_index_sorted_by_indicatrice(old_indicatrice, next_indicatrice);

  // Creation des tableaux des indices tries par le statut diphasique
  update_index_sorted_by_statut_diphasique(old_indicatrice, next_indicatrice);

  // Verifications
  assert(verifier_coherence_coord_linear_index());
  assert(verifier_taille_tableaux());
}

void Cut_cell_FT_Disc::update(const IJK_Field_double& old_indicatrice, const IJK_Field_double& next_indicatrice)
{
  // Increment du compteur des iterations
  processed_count_ += 1;

  // Supprime les elements virtuels
  linear_index_.resize(n_loc_, linear_index_.dimension(1));
  permutation_.resize(n_loc_, permutation_.dimension(1));
  processed_.resize(n_loc_, processed_.dimension(1));
  resize_data(n_loc_);
  n_tot_ = n_loc_;

  // Mise a jour des elements locaus
  n_loc_ = add_and_remove_local_elements(old_indicatrice, next_indicatrice);
  n_tot_ = n_loc_;

  // Remplissage de certains champs a partir des valeurs sur la structure IJK_Field
  set_coord();

  // Communication des elements virtuels
  desc_.reset();

  n_tot_ = initialise_communications();
  resize_data(n_tot_);

  echange_espace_virtuel();

  // Calcul de l'indice lineaire pour les elements virtuels
  linear_index_.resize(n_tot_, linear_index_.dimension(1));
  compute_virtual_linear_index();

  // Remplissage du tableau des indices diphasiques
  remplir_indice_diphasique();

  // Creation des tableaux des indices tries par k
  update_index_sorted_by_k();

  // Creation des tableaux des indices tries par l'indicatrice
  update_index_sorted_by_indicatrice(old_indicatrice, next_indicatrice);

  // Creation des tableaux des indices tries par le statut diphasique
  update_index_sorted_by_statut_diphasique(old_indicatrice, next_indicatrice);

  // Verifications
  assert(verifier_coherence_coord_linear_index());
  assert(verifier_taille_tableaux());

  Cerr << "(Cut_cell_FT_Disc::update) Iter: " << processed_count_-1 << " n_tot_=" << n_tot_  << " n_loc_=" << n_loc_ << finl;
}

// Initialise le tableau des indices a partir du champ de l'indicatrice
// et renvoie le nombre local d'elements diphasiques n_loc.
int Cut_cell_FT_Disc::initialise_linear_index(const IJK_Field_double& old_indicatrice, const IJK_Field_double& next_indicatrice)
{
  const int ni = splitting_.get_nb_elem_local(0);
  const int nj = splitting_.get_nb_elem_local(1);
  const int nk = splitting_.get_nb_elem_local(2);

  int n_loc = n_loc_;

  assert(n_loc == 0);
  assert(n_loc == n_tot_);

  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              if (!interfaces_.est_pure(0.5*(old_indicatrice(i,j,k) + next_indicatrice(i,j,k))))
                {
                  n_loc += 1;
                }
            }
        }
    }

  linear_index_.resize(n_loc, linear_index_.dimension(1));

  {
    int n = 0;
    for (int k = 0; k < nk; k++)
      {
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                if (!interfaces_.est_pure(0.5*(old_indicatrice(i,j,k) + next_indicatrice(i,j,k))))
                  {
                    linear_index_(n) = get_linear_index(i, j, k, ghost_size_, splitting_, true);

                    // Verification de la fonction function get_ijk_from_linear_index
                    assert(get_ijk_from_linear_index(linear_index_(n), ghost_size_, splitting_, true)[0] == i);
                    assert(get_ijk_from_linear_index(linear_index_(n), ghost_size_, splitting_, true)[1] == j);
                    assert(get_ijk_from_linear_index(linear_index_(n), ghost_size_, splitting_, true)[2] == k);

                    n += 1;
                  }
              }
          }
      }
    assert(n == n_loc);
  }

  return n_loc;
}

void Cut_cell_FT_Disc::initialise_permutation()
{
  for (int n = 0; n < n_loc_; n++)
    {
      permutation_(n) = n;
    }
}

void Cut_cell_FT_Disc::initialise_processed()
{
  processed_count_ = 1;
  for (int n = 0; n < n_loc_; n++)
    {
      processed_(n) = processed_count_;
    }
}

void Cut_cell_FT_Disc::set_coord()
{
  const int processed_sign = (processed_.dimension(0) == 0) ? 1 : (processed_(0) > 0 ? 1 : -1);
  assert(splitting_.get_grid_geometry().is_uniform(0));
  assert(splitting_.get_grid_geometry().is_uniform(1));
  assert(splitting_.get_grid_geometry().is_uniform(2));
  for (int n = 0; n < n_loc_; n++)
    {
      Int3 ijk = get_ijk_from_linear_index(linear_index_(n), ghost_size_, splitting_, true);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      // Pour les cellules anciennement diphasiques, il n'y a pas de sens
      // a fixer la valeur de coord_.
      if (processed_(n) * processed_sign >= processed_count_)
        {
          coord_(n,0) = (i + splitting_.get_offset_local(0) + .5)*splitting_.get_grid_geometry().get_constant_delta(0) + splitting_.get_grid_geometry().get_origin(0);
          coord_(n,1) = (j + splitting_.get_offset_local(1) + .5)*splitting_.get_grid_geometry().get_constant_delta(1) + splitting_.get_grid_geometry().get_origin(1);
          coord_(n,2) = (k + splitting_.get_offset_local(2) + .5)*splitting_.get_grid_geometry().get_constant_delta(2) + splitting_.get_grid_geometry().get_origin(2);

          // Verification de la fonction function get_ijk_from_linear_index
          assert(i == get_ijk_from_coord(coord_(n,0), coord_(n,1), coord_(n,2), ghost_size_, splitting_, true)[0]);
          assert(j == get_ijk_from_coord(coord_(n,0), coord_(n,1), coord_(n,2), ghost_size_, splitting_, true)[1]);
          assert(k == get_ijk_from_coord(coord_(n,0), coord_(n,1), coord_(n,2), ghost_size_, splitting_, true)[2]);
        }
    }
}

// Met a jour les tableaux des indices pour tenir compte des cellules
// nouvellement diphasiques et les cellules anciennement diphasiques a supprimer
// et renvoie le nombre local d'elements diphasiques n_loc.
int Cut_cell_FT_Disc::add_and_remove_local_elements(const IJK_Field_double& old_indicatrice, const IJK_Field_double& next_indicatrice)
{
  const int ni = splitting_.get_nb_elem_local(0);
  const int nj = splitting_.get_nb_elem_local(1);
  const int nk = splitting_.get_nb_elem_local(2);

  int n_loc = n_loc_;

  assert(n_loc == n_tot_);
  assert(verifier_si_ordonne(linear_index_, 0, n_loc));
  assert(verifier_coherence_coord_linear_index());

  // Calcule une borne superieure au nombre d'elements anciennement diphasiques a supprimer
  // Le nombre reeel pourrait etre plus faible si certains de ces elements redeviennent diphasiques.
  const int processed_sign = (processed_.dimension(0) == 0) ? 1 : (processed_(0) > 0 ? 1 : -1);

  // Mise a jour des tableaux linear_index_, permutation_ et processed_.
  // Procedure :
  // * Parcours du champ de l'indicatrice a la recherche de cellules diphasiques.
  //   Pour chaque cellule diphasique, verifier si elle deja dans la structure :
  //   les nouvelles cellules diphasiques sont ajoutes a la fin des tableaux,
  //   tout en notant l'endroit ou il aurait fallu les inserer dans permutation_.
  // * Permutation des tableaux pour retrouver l'ordre des indices croissants,
  //   en placant les anciennes cellules diphasiques a supprimer a la fin.
  //   Il suffit alors de diminuer la taille des tableaux pour les supprimer.

  int index_perm = 0;

  int n_initial = n_loc;
  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              if (!interfaces_.est_pure(0.5*(old_indicatrice(i,j,k) + next_indicatrice(i,j,k))))
                {
                  int linear_index = get_linear_index(i, j, k, ghost_size_, splitting_, true);

                  // Recherche de l'element dans la structure cut-cell
                  // La zone de recherche est restreinte car on sait que l'on va trouver les
                  // elements dans l'ordre du tableau
                  int m = indice_diphasique_(i,j,k);

                  if (m < 0)
                    {
                      // Nouvel element diphasique trouve
                      // Ajout du nouvel element a la fin de linear_index_ et processed_
                      linear_index_.resize(n_loc + 1, linear_index_.dimension(1));
                      permutation_.resize(n_loc + 1, permutation_.dimension(1));
                      processed_.resize(n_loc + 1, processed_.dimension(1));

                      linear_index_(n_loc) = linear_index;
                      processed_(n_loc) = processed_count_ * processed_sign;

                      // Ajout a permutation_ : nouvel element diphasique trouve
                      permutation_(index_perm) = n_loc;
                      index_perm += 1;
                      n_loc += 1;
                    }
                  else
                    {
                      // Ajout a permutation_ : element diphasique trouve
                      permutation_(index_perm) = m;
                      index_perm += 1;
                    }

                }
            }
        }
    }

  // Ajout a permutation_ : elements anciennement diphasiques a supprimer
  int n_deletion = 0;
  for (int n = 0; n < n_initial; n++)
    {
      Int3 ijk = get_ijk_from_coord(coord_(n,0), coord_(n,1), coord_(n,2), ghost_size_, splitting_, true);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      if (interfaces_.est_pure(0.5*(old_indicatrice(i,j,k) + next_indicatrice(i,j,k))))
        {
          permutation_(index_perm) = n;
          index_perm += 1;

          n_deletion += 1;
        }
    }

  assert(index_perm == n_loc);
  assert(verifier_valide_permutation(permutation_));

  // Application des differentes permutations
  apply_permutation(processed_, permutation_, processed_);
  apply_permutation(linear_index_, permutation_, processed_);

  assert(verifier_si_ordonne(linear_index_, 0, n_loc - n_deletion));

  assert(persistent_double_data_(0).valeur().dimension(0) <= n_loc);
  resize_data(n_loc);

  for (int i = 0; i < persistent_double_data_.size(); i++)
    {
      apply_permutation(persistent_double_data_(i).valeur(), permutation_, processed_);
    }

  for (int i = 0; i < persistent_int_data_.size(); i++)
    {
      apply_permutation(persistent_int_data_(i).valeur(), permutation_, processed_);
    }

  // Redimensionnement pour suppression des elements anciennement diphasiques

  resize_data(n_loc - n_deletion);

  linear_index_.resize(n_loc - n_deletion, linear_index_.dimension(1));
  permutation_.resize(n_loc - n_deletion, permutation_.dimension(1));
  processed_.resize(n_loc - n_deletion, permutation_.dimension(1));

  n_loc -= n_deletion;

  assert(verifier_si_ordonne(linear_index_, 0, n_loc));

  return n_loc;
}

void Cut_cell_FT_Disc::resize_data(int size)
{
  // Changement de taille des tableaux de double
  for (int i = 0; i < persistent_double_data_.size(); i++)
    {
      int n_initial = persistent_double_data_(i).valeur().dimension(0);

      persistent_double_data_(i).valeur().resize(size, persistent_double_data_(i).valeur().dimension(1));

      // Si la nouvelle taille est plus grande, initialise les valeurs
      for (int n = n_initial; n < size; n++)
        {
          for (int j = 0; j < persistent_double_data_(i).valeur().dimension(1); j++)
            {
              persistent_double_data_(i).valeur()(n,j) = 0;
            }
        }
    }

  for (int i = 0; i < transient_double_data_.size(); i++)
    {
      transient_double_data_(i).valeur().resize(size, transient_double_data_(i).valeur().dimension(1));

      // Initialise toutes les valeurs
      for (int n = 0; n < size; n++)
        {
          for (int j = 0; j < transient_double_data_(i).valeur().dimension(1); j++)
            {
              transient_double_data_(i).valeur()(n,j) = 0;
            }
        }
    }

  for (int i = 0; i < lazy_double_data_.size(); i++)
    {
      lazy_double_data_(i).valeur().resize(size, lazy_double_data_(i).valeur().dimension(1));
    }

  // Changement de taille des tableaux d'entier
  for (int i = 0; i < persistent_int_data_.size(); i++)
    {
      int n_initial = persistent_int_data_(i).valeur().dimension(0);

      persistent_int_data_(i).valeur().resize(size, persistent_int_data_(i).valeur().dimension(1));

      // Si la nouvelle taille est plus grande, initialise les valeurs
      for (int n = n_initial; n < size; n++)
        {
          for (int j = 0; j < persistent_int_data_(i).valeur().dimension(1); j++)
            {
              persistent_int_data_(i).valeur()(n,j) = 0;
            }
        }
    }

  for (int i = 0; i < transient_int_data_.size(); i++)
    {
      transient_int_data_(i).valeur().resize(size, transient_int_data_(i).valeur().dimension(1));

      // Initialise toutes les valeurs
      for (int n = 0; n < size; n++)
        {
          for (int j = 0; j < transient_int_data_(i).valeur().dimension(1); j++)
            {
              transient_int_data_(i).valeur()(n,j) = 0;
            }
        }
    }

  for (int i = 0; i < lazy_int_data_.size(); i++)
    {
      lazy_int_data_(i).valeur().resize(size, lazy_int_data_(i).valeur().dimension(1));
    }
}

void Cut_cell_FT_Disc::compute_virtual_linear_index()
{
  assert(linear_index_.dimension(0) == n_tot_);

  int linear_index;

  // Verification de l'indice lineaire pour les elements locaux
  for (int n = 0; n < n_loc_; n++)
    {
      Int3 ijk = get_ijk_from_coord(coord_(n,0), coord_(n,1), coord_(n,2), ghost_size_, splitting_, true);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      linear_index = get_linear_index(i, j, k, ghost_size_, splitting_, true);

      assert(linear_index_(n) == linear_index);
    }

  // Calcul de l'indice lineaire pour les elements virtuels
  for (int n = n_loc_; n < n_tot_; n++)
    {
      Int3 ijk = get_ijk_from_coord(coord_(n,0), coord_(n,1), coord_(n,2), ghost_size_, splitting_, false);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      linear_index = get_linear_index(i, j, k, ghost_size_, splitting_, false);

      linear_index_(n) = linear_index;
    }
}

void Cut_cell_FT_Disc::initialise_schema_comm()
{
  ArrOfIntFT pe_list;
  for (int next = 0; next < 2; next++)
    {
      for (int direction = 0; direction < 3; direction++)
        {
          for (int next_2 = 0; next_2 < 2; next_2++)
            {
              for (int direction_2 = direction + 1; direction_2 < 4; direction_2++) // Le case direction_2 = 3 correspond a aucun effet de next_2
                {
                  for (int next_3 = 0; next_3 < 2; next_3++)
                    {
                      for (int direction_3 = std::min(3, direction_2 + 1); direction_3 < 4; direction_3++) // Le case direction_3 = 3 correspond a aucun effet de next_3
                        {
                          if ((direction_2 == 3) && (next_2 == 1))
                            continue;
                          if ((direction_3 == 3) && (next_3 == 1))
                            continue;

                          int dest_pe_position_x = splitting_.get_local_slice_index(0) + (2*next-1)*(direction == 0) + (2*next_2-1)*(direction_2 == 0) + (2*next_3-1)*(direction_3 == 0);
                          int dest_pe_position_y = splitting_.get_local_slice_index(1) + (2*next-1)*(direction == 1) + (2*next_2-1)*(direction_2 == 1) + (2*next_3-1)*(direction_3 == 1);
                          int dest_pe_position_z = splitting_.get_local_slice_index(2) + (2*next-1)*(direction == 2) + (2*next_2-1)*(direction_2 == 2) + (2*next_3-1)*(direction_3 == 2);
                          int dest_pe = periodic_get_processor_by_ijk(dest_pe_position_x, dest_pe_position_y, dest_pe_position_z);
                          assert((direction_2 != 3 || direction_3 != 3) || (dest_pe == splitting_.get_neighbour_processor(next, direction)));
                          if ((dest_pe == splitting_.me()) || (dest_pe == -1))
                            continue;

                          pe_list.append_array(dest_pe);
                        }
                    }
                }
            }
        }
    }
  array_trier_retirer_doublons(pe_list);
  schema_comm_.set_send_recv_pe_list(pe_list, pe_list);
}

int Cut_cell_FT_Disc::initialise_communications()
{
  const int ni = splitting_.get_nb_elem_local(0);
  const int nj = splitting_.get_nb_elem_local(1);
  const int nk = splitting_.get_nb_elem_local(2);

  Descripteur_FT& espace_distant = desc_.espace_distant();

  schema_comm_.begin_comm();

  int n_loc = n_loc_;
  int n_tot = n_tot_;

  assert(n_tot == n_loc);

  // next indique la direction du processeur :
  //   0    voisin vers les petits indices
  //   1    voisin vers les grands indices
  for (int next = 0; next < 2; next++)
    {
      for (int direction = 0; direction < 3; direction++)
        {
          for (int next_2 = 0; next_2 < 2; next_2++)
            {
              for (int direction_2 = direction + 1; direction_2 < 4; direction_2++) // Le case direction_2 = 3 correspond a aucun effet de next_2
                {
                  for (int next_3 = 0; next_3 < 2; next_3++)
                    {
                      for (int direction_3 = std::min(3, direction_2 + 1); direction_3 < 4; direction_3++) // Le case direction_3 = 3 correspond a aucun effet de next_3
                        {
                          if ((direction_2 == 3) && (next_2 == 1))
                            continue;
                          if ((direction_3 == 3) && (next_3 == 1))
                            continue;

                          int dest_pe_position_x = splitting_.get_local_slice_index(0) + (2*next-1)*(direction == 0) + (2*next_2-1)*(direction_2 == 0) + (2*next_3-1)*(direction_3 == 0);
                          int dest_pe_position_y = splitting_.get_local_slice_index(1) + (2*next-1)*(direction == 1) + (2*next_2-1)*(direction_2 == 1) + (2*next_3-1)*(direction_3 == 1);
                          int dest_pe_position_z = splitting_.get_local_slice_index(2) + (2*next-1)*(direction == 2) + (2*next_2-1)*(direction_2 == 2) + (2*next_3-1)*(direction_3 == 2);
                          int dest_pe = periodic_get_processor_by_ijk(dest_pe_position_x, dest_pe_position_y, dest_pe_position_z);
                          assert((direction_2 != 3 || direction_3 != 3) || (dest_pe == splitting_.get_neighbour_processor(next, direction)));
                          if ((dest_pe == splitting_.me()) || (dest_pe == -1))
                            continue;

                          int ni_dir = select(direction, ni, nj, nk);
                          int ni_dir2 = (direction_2 == 3) ? 0 : select(direction_2, ni, nj, nk);
                          int ni_dir3 = (direction_3 == 3) ? 0 : select(direction_3, ni, nj, nk);

                          Sortie& send_buffer = schema_comm_.send_buffer(dest_pe);

                          for (int n = 0; n < n_loc; n++)
                            {
                              int i_selon_dir = get_i_selon_dir(direction, coord_(n,direction), ghost_size_, splitting_, true, false);
                              int i_selon_dir2 = (direction_2 == 3) ? 0 : get_i_selon_dir(direction_2, coord_(n,direction_2), ghost_size_, splitting_, true, false);
                              int i_selon_dir3 = (direction_3 == 3) ? 0 : get_i_selon_dir(direction_3, coord_(n,direction_3), ghost_size_, splitting_, true, false);
                              if ((!next) && (i_selon_dir >= ghost_size_))
                                continue;
                              if ((next) && (i_selon_dir < ni_dir - ghost_size_))
                                continue;
                              if ((direction_2 != 3) && (!next_2) && (i_selon_dir2 >= ghost_size_))
                                continue;
                              if ((direction_2 != 3) && (next_2) && (i_selon_dir2 < ni_dir2 - ghost_size_))
                                continue;
                              if ((direction_3 != 3) && (!next_3) && (i_selon_dir3 >= ghost_size_))
                                continue;
                              if ((direction_3 != 3) && (next_3) && (i_selon_dir3 < ni_dir3 - ghost_size_))
                                continue;

                              assert(! espace_distant.contient_element(dest_pe, n));

                              espace_distant.ajoute_element(dest_pe, n);
                              double x = coord_(n,0);
                              double y = coord_(n,1);
                              double z = coord_(n,2);
                              send_buffer << n << x << y << z;
                            }
                        }
                    }
                }
            }
        }
    }

  schema_comm_.echange_taille_et_messages();

  Descripteur_FT& espace_virtuel = desc_.espace_virtuel();

  const ArrOfInt& recv_pe_list = schema_comm_.get_recv_pe_list();
  for (int indice_pe = 0; indice_pe < recv_pe_list.size_array(); indice_pe++)
    {
      const int pe_source = recv_pe_list[indice_pe];
      Entree& recv_buffer = schema_comm_.recv_buffer(pe_source);

      do
        {
          int numero_sur_pe_source;
          double pos_x;
          double pos_y;
          double pos_z;
          recv_buffer >> numero_sur_pe_source >> pos_x >> pos_y >> pos_z;
          if (recv_buffer.eof())
            break;

          int nsom = coord_.dimension(0);
          n_tot += 1;

          coord_.append_line(pos_x, pos_y, pos_z);
          espace_virtuel.ajoute_element(pe_source, nsom);
        }
      while (1);
    }

  schema_comm_.end_comm();

  espace_distant.calcul_liste_pe_voisins();
  espace_virtuel.calcul_liste_pe_voisins();
  desc_.calcul_schema_comm(coord_.dimension(0));

  return n_tot;
}

void Cut_cell_FT_Disc::echange_espace_virtuel()
{
  // Echange des tableaux de double
  for (int i = 0; i < persistent_double_data_.size(); i++)
    {
      desc_.echange_espace_virtuel(persistent_double_data_(i).valeur());
    }

  // Echange des tableaux d'entier
  for (int i = 0; i < persistent_int_data_.size(); i++)
    {
      desc_.echange_espace_virtuel(persistent_int_data_(i).valeur());
    }
}

void Cut_cell_FT_Disc::echange_espace_virtuel(MD_Vector_tools::Operations_echange op)
{
  // Echange des tableaux de double
  for (int i = 0; i < persistent_double_data_.size(); i++)
    {
      desc_.echange_espace_virtuel(persistent_double_data_(i).valeur(), op);
    }

  // Echange des tableaux d'entier
  for (int i = 0; i < persistent_int_data_.size(); i++)
    {
      desc_.echange_espace_virtuel(persistent_int_data_(i).valeur(), op);
    }
}

void Cut_cell_FT_Disc::update_index_sorted_by_k()
{
  // Remplissage du tableau
  for (int n = 0; n < n_loc_; n++)
    {
      Int3 ijk = get_ijk_from_linear_index(linear_index_(n), ghost_size_, splitting_, true);
      int k = ijk[2];

      index_sorted_by_k_(n,0) = n;
      index_sorted_by_k_(n,1) = k;
    }

  for (int n = n_loc_; n < n_tot_; n++)
    {
      Int3 ijk = get_ijk_from_linear_index(linear_index_(n), ghost_size_, splitting_, false);
      int k = ijk[2];

      index_sorted_by_k_(n,0) = n;
      index_sorted_by_k_(n,1) = k;
    }

  // Trie selon la colonne des k
  index_sorted_by_k_.sort_tot(1);

  // Trouve les indices des differents k
  const int nk_tot = splitting_.get_nb_elem_local(2) + 2*ghost_size_;
  k_value_index_.resize(nk_tot+1, 1);

  int previous_k = -ghost_size_;
  k_value_index_(previous_k + ghost_size_) = 0;
  for (int index = 0; index < n_tot_; index++)
    {
      int k = index_sorted_by_k_(index,1);

      if (k > previous_k)
        {
          for (int k_intermediaire = previous_k+1; k_intermediaire <= k; k_intermediaire++)
            {
              k_value_index_(k_intermediaire + ghost_size_) = index;
            }
          previous_k = k;
        }
    }
  for (int k_intermediaire = previous_k+1; k_intermediaire <= splitting_.get_nb_elem_local(2) + ghost_size_; k_intermediaire++)
    {
      k_value_index_(k_intermediaire + ghost_size_) = n_tot_;
    }
}

void Cut_cell_FT_Disc::update_index_sorted_by_indicatrice(const IJK_Field_double& old_indicatrice, const IJK_Field_double& next_indicatrice)
{
  // Remplissage du tableau
  for (int n = 0; n < n_loc_; n++)
    {
      Int3 ijk = get_ijk_from_linear_index(linear_index_(n), ghost_size_, splitting_, true);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      union
      {
        double d;
        long long int i;
      } u;
      u.i = (long long int)n;

      index_sorted_by_indicatrice_(n,0) = u.d; // Storing the bits of the (int64) index of the cell into a double variable
      index_sorted_by_indicatrice_(n,1) = std::min(next_indicatrice(i,j,k), 1 - next_indicatrice(i,j,k));
      index_sorted_by_indicatrice_(n,2) = std::min(old_indicatrice(i,j,k), 1 - old_indicatrice(i,j,k));
    }

  for (int n = n_loc_; n < n_tot_; n++)
    {
      Int3 ijk = get_ijk_from_linear_index(linear_index_(n), ghost_size_, splitting_, false);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      union
      {
        double d;
        long long int i;
      } u;
      u.i = (long long int)n;

      index_sorted_by_indicatrice_(n,0) = u.d; // Storing the bits of the (int64) index of the cell into a double variable
      index_sorted_by_indicatrice_(n,1) = std::min(next_indicatrice(i,j,k), 1 - next_indicatrice(i,j,k));
      index_sorted_by_indicatrice_(n,2) = std::min(old_indicatrice(i,j,k), 1 - old_indicatrice(i,j,k));
    }

  // Trie selon la colonne des indicatrices
  index_sorted_by_indicatrice_.sort_tot(1,2);
}

void Cut_cell_FT_Disc::update_index_sorted_by_statut_diphasique(const IJK_Field_double& old_indicatrice, const IJK_Field_double& next_indicatrice)
{
  for (int index = 0; index < n_tot_; index++)
    {
      int n = get_n_from_indicatrice_index(index);

      Int3 ijk = get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      bool est_reguliere = interfaces_.est_reguliere(old_indicatrice(i,j,k), next_indicatrice(i,j,k));
      bool devient_pure = interfaces_.devient_pure(old_indicatrice(i,j,k), next_indicatrice(i,j,k));
      bool devient_diphasique = interfaces_.devient_diphasique(old_indicatrice(i,j,k), next_indicatrice(i,j,k));
      bool desequilibre_final = interfaces_.a_desequilibre_final(old_indicatrice(i,j,k), next_indicatrice(i,j,k));
      bool desequilibre_initial_uniquement = interfaces_.a_desequilibre_initial_uniquement(old_indicatrice(i,j,k), next_indicatrice(i,j,k));
      assert(est_reguliere + devient_pure + devient_diphasique + desequilibre_final + desequilibre_initial_uniquement == 1);

      index_sorted_by_statut_diphasique_(n,0) = n;
      index_sorted_by_statut_diphasique_(n,1) = est_reguliere*static_cast<int>(STATUT_DIPHASIQUE::REGULIER) + devient_pure*static_cast<int>(STATUT_DIPHASIQUE::MOURRANT) + devient_diphasique*static_cast<int>(STATUT_DIPHASIQUE::NAISSANT)
                                                + desequilibre_final*static_cast<int>(STATUT_DIPHASIQUE::DESEQUILIBRE_FINAL) + desequilibre_initial_uniquement*static_cast<int>(STATUT_DIPHASIQUE::DESEQUILIBRE_INITIAL);
      assert(index_sorted_by_statut_diphasique_(n,1) < STATUT_DIPHASIQUE::count);
      index_sorted_by_statut_diphasique_(n,2) = index;
    }

  // Trie selon la colonne des statut diphasiques
  // Le tri secondaire est selon l'indice du tableau trie selon l'indicatrice, c'est-a-dire
  // dans le sens de l'indicatrice (de la plus petite phase) croissante.
  index_sorted_by_statut_diphasique_.sort_tot(1,2);

  // Trouve les indices des differents statut dipashiques
  statut_diphasique_value_index_.resize(STATUT_DIPHASIQUE::count+1, 1);

  int previous_statut_diphasique = 0;
  statut_diphasique_value_index_(previous_statut_diphasique) = 0;
  for (int index = 0; index < n_tot_; index++)
    {
      int statut_diphasique = index_sorted_by_statut_diphasique_(index,1);

      if (statut_diphasique > previous_statut_diphasique)
        {
          for (int k_intermediaire = previous_statut_diphasique+1; k_intermediaire <= statut_diphasique; k_intermediaire++)
            {
              statut_diphasique_value_index_(k_intermediaire) = index;
            }
          previous_statut_diphasique = statut_diphasique;
        }
    }
  for (int k_intermediaire = previous_statut_diphasique+1; k_intermediaire <= static_cast<int>(STATUT_DIPHASIQUE::count); k_intermediaire++)
    {
      statut_diphasique_value_index_(k_intermediaire) = n_tot_;
    }
}

void Cut_cell_FT_Disc::remplir_indice_diphasique()
{
  const int ni = splitting_.get_nb_elem_local(0);
  const int nj = splitting_.get_nb_elem_local(1);
  const int nk = splitting_.get_nb_elem_local(2);

  for (int k = -ghost_size_; k < nk+ghost_size_; k++)
    {
      for (int j = -ghost_size_; j < nj+ghost_size_; j++)
        {
          for (int i = -ghost_size_; i < ni+ghost_size_; i++)
            {
              indice_diphasique_(i,j,k) = -1;
            }
        }
    }

  for (int n = 0; n < n_loc_; n++)
    {
      Int3 ijk_no_per = get_ijk_from_linear_index(linear_index_(n), ghost_size_, splitting_, true);
      int i = ijk_no_per[0];
      int j = ijk_no_per[1];
      int k = ijk_no_per[2];

      {
        int index_ijk_per = 0;
        while (index_ijk_per >= 0)
          {
            Int3 ijk = ijk_per_of_index(i, j, k, index_ijk_per);
            index_ijk_per = next_index_ijk_per(i, j, k, index_ijk_per, ghost_size_, ghost_size_);

            indice_diphasique_(ijk[0],ijk[1],ijk[2]) = n;
          }
      }
    }

  for (int n = n_loc_; n < n_tot_; n++)
    {
      Int3 ijk_no_per = get_ijk_from_linear_index(linear_index_(n), ghost_size_, splitting_, false);
      int i = ijk_no_per[0];
      int j = ijk_no_per[1];
      int k = ijk_no_per[2];

      {
        int index_ijk_per = 0;
        while (index_ijk_per >= 0)
          {
            Int3 ijk = ijk_per_of_index(i, j, k, index_ijk_per);
            index_ijk_per = next_index_ijk_per(i, j, k, index_ijk_per, ghost_size_, ghost_size_);

            indice_diphasique_(ijk[0],ijk[1],ijk[2]) = n;
          }
      }
    }
}

void Cut_cell_FT_Disc::remove_dead_and_virtual_cells(const IJK_Field_double& next_indicatrice)
{
  // Annule l'information de l'indice diphasique pour les cellules pures au pas suivant.
  // Motivation : Les cellules qui meurrent et reapparaissent au pas suivant doivent
  // etre supprimees puis reajoutees.

  linear_index_.resize(n_loc_, linear_index_.dimension(1));
  permutation_.resize(n_loc_, permutation_.dimension(1));
  processed_.resize(n_loc_, processed_.dimension(1));
  resize_data(n_loc_);
  n_tot_ = n_loc_;

  int index_perm = 0;
  int n_deletion = 0;

  for (int n = 0; n < n_loc_; n++)
    {
      Int3 ijk = get_ijk_from_linear_index(linear_index_(n), ghost_size_, splitting_, true);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      if ((next_indicatrice(i,j,k) == 0) || (next_indicatrice(i,j,k) == 1))
        {
          n_deletion += 1;
          permutation_(n_loc_ - n_deletion) = n;
        }
      else
        {
          permutation_(index_perm) = n;
          index_perm += 1;
        }
    }

  assert(index_perm + n_deletion  == n_loc_);
  assert(verifier_valide_permutation(permutation_));

  // Application des differentes permutations
  apply_permutation(processed_, permutation_, processed_);
  apply_permutation(linear_index_, permutation_, processed_);

  assert(verifier_si_ordonne(linear_index_, 0, n_loc_ - n_deletion));

  assert(persistent_double_data_(0).valeur().dimension(0) <= n_loc_);
  resize_data(n_loc_);

  for (int i = 0; i < persistent_double_data_.size(); i++)
    {
      apply_permutation(persistent_double_data_(i).valeur(), permutation_, processed_);
    }

  for (int i = 0; i < persistent_int_data_.size(); i++)
    {
      apply_permutation(persistent_int_data_(i).valeur(), permutation_, processed_);
    }

  // Redimensionnement pour suppression des elements anciennement diphasiques
  resize_data(n_loc_ - n_deletion);
  linear_index_.resize(n_loc_ - n_deletion, linear_index_.dimension(1));
  permutation_.resize(n_loc_ - n_deletion, permutation_.dimension(1));
  processed_.resize(n_loc_ - n_deletion, permutation_.dimension(1));

  n_loc_ -= n_deletion;

  assert(verifier_si_ordonne(linear_index_, 0, n_loc_));

  n_tot_ = n_loc_;

  remplir_indice_diphasique();
}

template<typename T>
void Cut_cell_FT_Disc::fill_buffer_with_variable(const TRUSTTabFT<T>& array, int component /* = 0 by default */) const
{
  assert(component < array.dimension(1));
  const int ni = splitting_.get_nb_elem_local(0);
  const int nj = splitting_.get_nb_elem_local(1);
  const int nk = splitting_.get_nb_elem_local(2);

  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              write_buffer_(i,j,k) = 0.;
            }
        }
    }

  for (int n = 0; n < n_loc_; n++)
    {
      Int3 ijk = get_ijk_from_linear_index(linear_index_(n), ghost_size_, splitting_, true);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      write_buffer_(i,j,k) = (double)(array(n,component));
    }
}

// Instantiation explicite pour IntTabFT et DoubleTabFT
template void Cut_cell_FT_Disc::fill_buffer_with_variable<int>(const TRUSTTabFT<int>&, int) const;
template void Cut_cell_FT_Disc::fill_buffer_with_variable<double>(const TRUSTTabFT<double>&, int) const;

template<typename T>
void Cut_cell_FT_Disc::fill_variable_with_buffer(TRUSTTabFT<T>& array, int component /* = 0 by default */) const
{
  assert(component < array.dimension(1));

  for (int n = 0; n < n_loc_; n++)
    {
      Int3 ijk = get_ijk_from_linear_index(linear_index_(n), ghost_size_, splitting_, true);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      array(n,component) = (T)write_buffer_(i,j,k);
    }
}

// Instantiation explicite pour IntTabFT et DoubleTabFT
template void Cut_cell_FT_Disc::fill_variable_with_buffer<int>(TRUSTTabFT<int>&, int) const;
template void Cut_cell_FT_Disc::fill_variable_with_buffer<double>(TRUSTTabFT<double>&, int) const;

bool Cut_cell_FT_Disc::verifier_coherence_coord_linear_index()
{
  // Verifie la coherence entre les tableaux coord_ et linear_index_
  for (int n = 0; n < n_loc_; n++)
    {
      Int3 ijk_1 = get_ijk_from_coord(coord_(n,0), coord_(n,1), coord_(n,2), ghost_size_, splitting_, true);
      int i_1 = ijk_1[0];
      int j_1 = ijk_1[1];
      int k_1 = ijk_1[2];
      Int3 ijk_2 = get_ijk_from_linear_index(linear_index_(n), ghost_size_, splitting_, true);
      int i_2 = ijk_2[0];
      int j_2 = ijk_2[1];
      int k_2 = ijk_2[2];

      if ((i_1 != i_2) || (j_1 != j_2) || (k_1 != k_2))
        return false;
    }

  for (int n = n_loc_; n < n_tot_; n++)
    {
      Int3 ijk_1 = get_ijk_from_coord(coord_(n,0), coord_(n,1), coord_(n,2), ghost_size_, splitting_, false);
      int i_1 = ijk_1[0];
      int j_1 = ijk_1[1];
      int k_1 = ijk_1[2];
      Int3 ijk_2 = get_ijk_from_linear_index(linear_index_(n), ghost_size_, splitting_, false);
      int i_2 = ijk_2[0];
      int j_2 = ijk_2[1];
      int k_2 = ijk_2[2];

      if ((i_1 != i_2) || (j_1 != j_2) || (k_1 != k_2))
        return false;
    }
  return true;
}

bool Cut_cell_FT_Disc::verifier_taille_tableaux()
{
  if (linear_index_.dimension(0) != n_tot_)
    return false;

  for (int i = 0; i < persistent_double_data_.size(); i++)
    {
      if (persistent_double_data_(i).valeur().dimension(0) != n_tot_)
        return false;
    }

  for (int i = 0; i < transient_double_data_.size(); i++)
    {
      if (transient_double_data_(i).valeur().dimension(0) != n_tot_)
        return false;
    }

  for (int i = 0; i < lazy_double_data_.size(); i++)
    {
      if (lazy_double_data_(i).valeur().dimension(0) != n_tot_)
        return false;
    }

  for (int i = 0; i < persistent_int_data_.size(); i++)
    {
      if (persistent_int_data_(i).valeur().dimension(0) != n_tot_)
        return false;
    }

  for (int i = 0; i < transient_int_data_.size(); i++)
    {
      if (transient_int_data_(i).valeur().dimension(0) != n_tot_)
        return false;
    }

  for (int i = 0; i < lazy_int_data_.size(); i++)
    {
      if (lazy_int_data_(i).valeur().dimension(0) != n_tot_)
        return false;
    }
  return true;
}

void Cut_cell_FT_Disc::imprime_elements_diphasiques()
{
  Cerr << "# Impression des elements diphasiques pour le PE#" << splitting_.me() << " [N IJK XYZ T]" << finl;
  for (int n = 0; n < n_loc_; n++)
    {
      Int3 ijk = get_ijk_from_linear_index(linear_index_(n), ghost_size_, splitting_, true);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      // symbole qualifie le type d'element a l'impression :
      // (vide) element local
      //   %    element virtuel
      Nom symbole = (n < n_loc_ ? "" : "%");

      Cerr << n << symbole << " " << i << " " << j << " " << k << " " << coord_(n,0) << " " << coord_(n,1) << " " << coord_(n,2) << finl;
    }
  for (int n = n_loc_; n < n_tot_; n++)
    {
      Int3 ijk = get_ijk_from_linear_index(linear_index_(n), ghost_size_, splitting_, false);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      // symbole qualifie le type d'element a l'impression :
      // (vide) element local
      //   %    element virtuel
      Nom symbole = (n < n_loc_ ? "" : "%");

      Cerr << n << symbole << " " << i << " " << j << " " << k << " " << coord_(n,0) << " " << coord_(n,1) << " " << coord_(n,2) << finl;
    }
  Cerr << "# Fin impression elements diphasiques" << finl;
}

void Cut_cell_FT_Disc::imprime_elements_distants()
{
  Cerr << "# Impression des elements distants [N PE IJK XYZ T]" << finl;
  const int ni = splitting_.get_nb_elem_local(0);
  const int nj = splitting_.get_nb_elem_local(1);
  const int nk = splitting_.get_nb_elem_local(2);

  // next indique la direction du processeur :
  //   0    voisin vers les petits indices
  //   1    voisin vers les grands indices
  for (int next = 0; next < 2; next++)
    {
      for (int direction = 0; direction < 3; direction++)
        {
          for (int next_2 = 0; next_2 < 2; next_2++)
            {
              for (int direction_2 = direction + 1; direction_2 < 4; direction_2++) // Le case direction_2 = 3 correspond a aucun effet de next_2
                {
                  for (int next_3 = 0; next_3 < 2; next_3++)
                    {
                      for (int direction_3 = std::min(3, direction_2 + 1); direction_3 < 4; direction_3++) // Le case direction_3 = 3 correspond a aucun effet de next_3
                        {
                          if ((direction_2 == 3) && (next_2 == 1))
                            continue;
                          if ((direction_3 == 3) && (next_3 == 1))
                            continue;

                          int dest_pe_position_x = splitting_.get_local_slice_index(0) + (2*next-1)*(direction == 0) + (2*next_2-1)*(direction_2 == 0) + (2*next_3-1)*(direction_3 == 0);
                          int dest_pe_position_y = splitting_.get_local_slice_index(1) + (2*next-1)*(direction == 1) + (2*next_2-1)*(direction_2 == 1) + (2*next_3-1)*(direction_3 == 1);
                          int dest_pe_position_z = splitting_.get_local_slice_index(2) + (2*next-1)*(direction == 2) + (2*next_2-1)*(direction_2 == 2) + (2*next_3-1)*(direction_3 == 2);
                          int dest_pe = periodic_get_processor_by_ijk(dest_pe_position_x, dest_pe_position_y, dest_pe_position_z);
                          assert((direction_2 != 3 || direction_3 != 3) || (dest_pe == splitting_.get_neighbour_processor(next, direction)));
                          if ((dest_pe == splitting_.me()) || (dest_pe == -1))
                            continue;

                          int ni_dir = select(direction, ni, nj, nk);
                          int ni_dir2 = (direction_2 == 3) ? 0 : select(direction_2, ni, nj, nk);
                          int ni_dir3 = (direction_3 == 3) ? 0 : select(direction_3, ni, nj, nk);

                          for (int n = 0; n < n_loc_; n++)
                            {
                              int i_selon_dir = get_i_selon_dir(direction, coord_(n,direction), ghost_size_, splitting_, true, false);
                              int i_selon_dir2 = (direction_2 == 3) ? 0 : get_i_selon_dir(direction_2, coord_(n,direction_2), ghost_size_, splitting_, true, false);
                              int i_selon_dir3 = (direction_3 == 3) ? 0 : get_i_selon_dir(direction_3, coord_(n,direction_3), ghost_size_, splitting_, true, false);
                              if ((!next) && (i_selon_dir >= ghost_size_))
                                continue;
                              if ((next) && (i_selon_dir < ni_dir - ghost_size_))
                                continue;
                              if ((direction_2 != 3) && (!next_2) && (i_selon_dir2 >= ghost_size_))
                                continue;
                              if ((direction_2 != 3) && (next_2) && (i_selon_dir2 < ni_dir2 - ghost_size_))
                                continue;
                              if ((direction_3 != 3) && (!next_3) && (i_selon_dir3 >= ghost_size_))
                                continue;
                              if ((direction_3 != 3) && (next_3) && (i_selon_dir3 < ni_dir3 - ghost_size_))
                                continue;

                              Int3 ijk = get_ijk_from_linear_index(linear_index_(n), ghost_size_, splitting_, true);
                              int i = ijk[0];
                              int j = ijk[1];
                              int k = ijk[2];

                              Cerr << n << " " << dest_pe << " " << i << " " << j << " " << k << " " << coord_(n,0) << " " << coord_(n,1) << " " << coord_(n,2) << finl;
                            }
                        }
                    }
                }
            }
        }
    }
  Cerr << "# Fin impression elements distants" << finl;
}

// Recherche d'une valeur dans un tableau non trie.
template<typename T>
int Cut_cell_FT_Disc::find_value_unsorted(T value, const TRUSTTabFT<T>& array, int imin, int imax)
{
  assert(array.dimension(1) == 1);
  assert(imax < array.dimension(0));

  for (int i = 0; i <= imax; i++)
    {
      if (array(i) == value)
        {
          return i;
        }
    }
  return -1;
}

// Instantiation explicite pour IntTabFT et DoubleTabFT
template int Cut_cell_FT_Disc::find_value_unsorted<int>(int, const TRUSTTabFT<int>&, int, int);
template int Cut_cell_FT_Disc::find_value_unsorted<double>(double, const TRUSTTabFT<double>&, int, int);

bool Cut_cell_FT_Disc::verifier_toujours_meme_signe_et_non_nul(const IntTabFT_cut_cell& array)
{
  assert(array.dimension(1) == 1);

  const int initial_sign = array(0) > 0 ? 1 : -1;
  for (int i = 0; i < array.dimension(0); i++)
    {
      const int sign = array(i) > 0 ? 1 : -1;
      if ((array(i) == 0) || (sign != initial_sign))
        return false;
    }
  return true;
}

bool Cut_cell_FT_Disc::verifier_tableau_jamais_nul(const DoubleTabFT_cut_cell& array)
{
  assert(array.dimension(1) == 3);

  for (int i = 0; i < array.dimension(0); i++)
    {
      if ((array(i,0) == 0) || (array(i,0) == 1) || (array(i,0) == 2))
        return false;
    }
  return true;
}

// Verifie que la permutation est valide, c'est-a-dire inclue tous les
// indices et une seule fois.
bool Cut_cell_FT_Disc::verifier_valide_permutation(const IntTabFT_cut_cell& array)
{
  assert(array.dimension(1) == 1);
  for (int i = 0; i < array.dimension(0); i++)
    {
      const int m = find_value_unsorted(i, array, 0, array.dimension(0) - 1);
      if (m < 0)
        return false;
    }
  return true;
}

// Verifie qu'il n'y a pas de doublons.
bool Cut_cell_FT_Disc::verifier_pas_de_doublons(const IntTabFT_cut_cell& array)
{
  assert(array.dimension(1) == 1);
  for (int i = 0; i < array.dimension(0); i++)
    {
      const int m = find_value_unsorted(array(i), array, 0, array.dimension(0) - 1);
      if (m != i)
        return false;
    }
  return true;
}

// Permutate le tableau array a partir des indices de permutation 'permutation'.
// Le tableau 'processed' permet de suivre les changements par son signe.
// La seule contrainte sur le tableau 'processed' est qu'il est toujours de meme signe et non nul.
// Le signe du tableau 'processed' est modifie par la fonction mais le tableau de permutation n'est pas modifie.
template<typename T>
void Cut_cell_FT_Disc::apply_permutation(TRUSTTabFT<T>& array, const IntTabFT_cut_cell& permutation, IntTabFT_cut_cell& processed)
{
  if (array.dimension(0) == 0)
    {
      assert(permutation.dimension(0) == 0);
      assert(processed.dimension(0) == 0);
    }
  else
    {
      assert(array.dimension(0) == permutation.dimension(0));
      assert(array.dimension(0) == processed.dimension(0));

      const int max_nb_columns = 6;
      assert(array.dimension(1) <= max_nb_columns);

      assert(verifier_toujours_meme_signe_et_non_nul(processed));
      int sign = processed(0) > 0 ? 1 : -1;

      for (int i = 0; i < array.dimension(0); i++)
        {
          if (sign*processed(i) < 0)
            {
              continue;
            }

          // L'element n'a pas deja ete traite, on stocke la valeur et on la remplace,
          // puis, on suit le tableau des permutations et on remplace chaque valeur
          // jusqu'a retrouver l'element initialement trouve.
          // Tous les elements traites sont marques dans 'processed', ce qui permet
          // de ne traiter par la suite que les elements d'un cycle different.

          T first_cycle_value[max_nb_columns];
          for (int j = 0; j < array.dimension(1); j++)
            {
              first_cycle_value[j] = array(i,j);
            }
          int current_index = i;
          int next_index = permutation[i];

          while (next_index != i)
            {
              for (int j = 0; j < array.dimension(1); j++)
                {
                  array(current_index,j) = array(next_index,j);
                }

              processed(current_index) = -processed(current_index);

              current_index = next_index;
              next_index = permutation[next_index];
            }

          for (int j = 0; j < array.dimension(1); j++)
            {
              array(current_index,j) = first_cycle_value[j];
            }
          processed(current_index) = -processed(current_index);
        }
    }
}

bool Cut_cell_FT_Disc::verifier_si_ordonne(const IntTabFT_cut_cell& array, int imin, int imax)
{
  int old_linear_index = -1;
  for (int n = imin; n < imax; n++)
    {
      int linear_index = array(n);
      if (old_linear_index >= linear_index)
        return false;
      old_linear_index = linear_index;
    }
  return true;
}

