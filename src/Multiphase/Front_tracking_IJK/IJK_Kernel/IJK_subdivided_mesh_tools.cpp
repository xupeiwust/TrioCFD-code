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

#include <IJK_subdivided_mesh_tools.h>

static void subdivide_array(const IJK_Grid_Geometry& geom1, const int direction, ArrOfDouble& delta, double& origin)
{
  delta = geom1.get_delta(direction);
  origin = geom1.get_origin(direction);

  // Modification du tableau delta
  // Le nombre de cellules est double et chaque cellule est deux foix plus petite
  const int n = delta.size_array(); // nombre de mailles initial
  delta.resize_array(2*n);
  int i;
  for (i = n - 1; i >= 0; i--)
    {
      double current_delta = delta[i]/2.;
      delta[2*i] = current_delta;
      delta[2*i + 1] = current_delta;
    }
}

// Division d'un splitting pour chaque maille dans chaque direction.
// split1 : Maillage sur le domaine etendu ou vivent les interfaces.
// split2 : Resultat subdivise utilise pour le parcours de l'interface.
void build_subdivided_splitting(const IJK_Splitting& split1, IJK_Splitting& split2)
{
  const IJK_Grid_Geometry& geom1 = split1.get_grid_geometry();

  double origin_x, origin_y, origin_z;
  ArrOfDouble dx, dy, dz;
  subdivide_array(geom1, DIRECTION_I, dx, origin_x);
  subdivide_array(geom1, DIRECTION_J, dy, origin_y);
  subdivide_array(geom1, DIRECTION_K, dz, origin_z);

  // Le domaine etendu n'est pas periodique: le champ n'est pas continu
  // entre les bords opposes du domaine etendu.
  IJK_Grid_Geometry geom2;
  Nom n(geom1.le_nom());
  geom2.nommer(n + "_X8");
  geom2.initialize_origin_deltas(origin_x, origin_y, origin_z, dx, dy, dz, geom1.get_periodic_flag(0), geom1.get_periodic_flag(1), geom1.get_periodic_flag(2));
  // Construction du decoupage parallele: on utilise les memes parametres
  // de decoupage que pour le maillage d'origine:
  split2.initialize(geom2, split1.get_nprocessor_per_direction(DIRECTION_I), split1.get_nprocessor_per_direction(DIRECTION_J), split1.get_nprocessor_per_direction(DIRECTION_K));
}

// Extension d'un champ ft vers 8x
void extend_ft_to_8x(const IJK_Field_double& field_ft, IJK_Field_double& field_8x)
{
  const int ni = field_ft.ni();
  const int nj = field_ft.nj();
  const int nk = field_ft.nk();
  const int ghost = field_ft.ghost();
  assert(field_8x.ghost() == 2*field_ft.ghost());
  for (int k = -ghost; k < nk+ghost; k++)
    {
      for (int j = -ghost; j < nj+ghost; j++)
        {
          for (int i = -ghost; i < ni+ghost; i++)
            {
              field_8x(2*i,2*j,2*k) = field_ft(i,j,k);
              field_8x(2*i+1,2*j,2*k) = field_ft(i,j,k);
              field_8x(2*i,2*j+1,2*k) = field_ft(i,j,k);
              field_8x(2*i+1,2*j+1,2*k) = field_ft(i,j,k);
              field_8x(2*i,2*j,2*k+1) = field_ft(i,j,k);
              field_8x(2*i+1,2*j,2*k+1) = field_ft(i,j,k);
              field_8x(2*i,2*j+1,2*k+1) = field_ft(i,j,k);
              field_8x(2*i+1,2*j+1,2*k+1) = field_ft(i,j,k);
            }
        }
    }
}

// Extension d'un champ ft vers 8x
void extend_ft_to_8x(const IJK_Field_vector3_double& field_ft, IJK_Field_vector3_double& field_8x)
{
  for (int dir=0; dir<3; dir++)
    {
      const int ni = field_ft[dir].ni();
      const int nj = field_ft[dir].nj();
      const int nk = field_ft[dir].nk();
      const int ghost = field_ft[dir].ghost();
      assert(field_8x[dir].ghost() == 2*field_ft[dir].ghost());
      for (int k = -ghost; k < nk+ghost; k++)
        {
          for (int j = -ghost; j < nj+ghost; j++)
            {
              for (int i = -ghost; i < ni+ghost; i++)
                {
                  field_8x[dir](2*i,2*j,2*k) = field_ft[dir](i,j,k);
                  field_8x[dir](2*i+1,2*j,2*k) = field_ft[dir](i,j,k);
                  field_8x[dir](2*i,2*j+1,2*k) = field_ft[dir](i,j,k);
                  field_8x[dir](2*i+1,2*j+1,2*k) = field_ft[dir](i,j,k);
                  field_8x[dir](2*i,2*j,2*k+1) = field_ft[dir](i,j,k);
                  field_8x[dir](2*i+1,2*j,2*k+1) = field_ft[dir](i,j,k);
                  field_8x[dir](2*i,2*j+1,2*k+1) = field_ft[dir](i,j,k);
                  field_8x[dir](2*i+1,2*j+1,2*k+1) = field_ft[dir](i,j,k);
                }
            }
        }
    }
}

// Combine les sous-elements deux-a-deux dans un element ft.
void moyenne_sous_elements_elem(const IJK_Field_double& field_8x, IJK_Field_double& field_ft)
{
  const int ni = field_ft.ni();
  const int nj = field_ft.nj();
  const int nk = field_ft.nk();
  const int ghost = field_ft.ghost();
  assert(field_8x.ghost() == 2*field_ft.ghost());
  for (int k = -ghost; k < nk+ghost; k++)
    {
      for (int j = -ghost; j < nj+ghost; j++)
        {
          for (int i = -ghost; i < ni+ghost; i++)
            {
              field_ft(i,j,k) = .125*(field_8x(2*i,2*j,2*k)       + field_8x(2*i+1,2*j,2*k)
                                      + field_8x(2*i,2*j+1,2*k)   + field_8x(2*i+1,2*j+1,2*k)
                                      + field_8x(2*i,2*j,2*k+1)   + field_8x(2*i+1,2*j,2*k+1)
                                      + field_8x(2*i,2*j+1,2*k+1) + field_8x(2*i+1,2*j+1,2*k+1));
            }
        }
    }
}

// Combine les sous-elements deux-a-deux avec une certaine tolerance.
// Utile pour l'indicatrice, car cela reproduit la tolerance du calcul de l'indicatrice sur le maillage ft.
void moyenne_indicatrice_sous_elements_elem(const IJK_Field_double& field_8x, IJK_Field_double& field_ft, double tolerance)
{
  const int ni = field_ft.ni();
  const int nj = field_ft.nj();
  const int nk = field_ft.nk();
  const int ghost = field_ft.ghost();
  assert(field_8x.ghost() == 2*field_ft.ghost());
  for (int k = -ghost; k < nk+ghost; k++)
    {
      for (int j = -ghost; j < nj+ghost; j++)
        {
          for (int i = -ghost; i < ni+ghost; i++)
            {
              double indic = .125*(field_8x(2*i,2*j,2*k)       + field_8x(2*i+1,2*j,2*k)
                                   + field_8x(2*i,2*j+1,2*k)   + field_8x(2*i+1,2*j+1,2*k)
                                   + field_8x(2*i,2*j,2*k+1)   + field_8x(2*i+1,2*j,2*k+1)
                                   + field_8x(2*i,2*j+1,2*k+1) + field_8x(2*i+1,2*j+1,2*k+1));
              if (indic<tolerance)
                {
                  field_ft(i,j,k) = 0.;
                }
              else if ((1.-indic)<tolerance)
                {
                  field_ft(i,j,k) = 1.;
                }
              else
                {
                  field_ft(i,j,k) = indic;
                }
            }
        }
    }
}

// Combine les sous-elements deux-a-deux avec une certaine tolerance, sur un volume de controle ft decale d'une demi-maille dans une direction (c'est-a-dire centre sur une face).
// Utile pour l'indicatrice, car cela reproduit la tolerance du calcul de l'indicatrice sur le maillage ft.
void moyenne_indicatrice_sous_elements_face(int dir, const IJK_Field_double& field_8x, IJK_Field_double& field_ft, double tolerance)
{
  const int ni = field_ft.ni();
  const int nj = field_ft.nj();
  const int nk = field_ft.nk();
  const int ghost = field_ft.ghost();
  assert(field_8x.ghost() == 2*field_ft.ghost());
  if (dir == 0)
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  double indic = .125*(field_8x(2*i,2*j,2*k)       + field_8x(2*i-1,2*j,2*k)
                                       + field_8x(2*i,2*j+1,2*k)   + field_8x(2*i-1,2*j+1,2*k)
                                       + field_8x(2*i,2*j,2*k+1)   + field_8x(2*i-1,2*j,2*k+1)
                                       + field_8x(2*i,2*j+1,2*k+1) + field_8x(2*i-1,2*j+1,2*k+1));
                  if (indic<tolerance)
                    {
                      field_ft(i,j,k) = 0.;
                    }
                  else if ((1.-indic)<tolerance)
                    {
                      field_ft(i,j,k) = 1.;
                    }
                  else
                    {
                      field_ft(i,j,k) = indic;
                    }
                }
            }
        }
    }
  else if (dir == 1)
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  double indic = .125*(field_8x(2*i,2*j,2*k)       + field_8x(2*i+1,2*j,2*k)
                                       + field_8x(2*i,2*j-1,2*k)   + field_8x(2*i+1,2*j-1,2*k)
                                       + field_8x(2*i,2*j,2*k+1)   + field_8x(2*i+1,2*j,2*k+1)
                                       + field_8x(2*i,2*j-1,2*k+1) + field_8x(2*i+1,2*j-1,2*k+1));
                  if (indic<tolerance)
                    {
                      field_ft(i,j,k) = 0.;
                    }
                  else if ((1.-indic)<tolerance)
                    {
                      field_ft(i,j,k) = 1.;
                    }
                  else
                    {
                      field_ft(i,j,k) = indic;
                    }
                }
            }
        }
    }
  else if (dir == 2)
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  double indic = .125*(field_8x(2*i,2*j,2*k)       + field_8x(2*i+1,2*j,2*k)
                                       + field_8x(2*i,2*j+1,2*k)   + field_8x(2*i+1,2*j+1,2*k)
                                       + field_8x(2*i,2*j,2*k-1)   + field_8x(2*i+1,2*j,2*k-1)
                                       + field_8x(2*i,2*j+1,2*k-1) + field_8x(2*i+1,2*j+1,2*k-1));
                  if (indic<tolerance)
                    {
                      field_ft(i,j,k) = 0.;
                    }
                  else if ((1.-indic)<tolerance)
                    {
                      field_ft(i,j,k) = 1.;
                    }
                  else
                    {
                      field_ft(i,j,k) = indic;
                    }
                }
            }
        }
    }
  else
    {
      Cerr << "moyenne_indicatrice_sous_elements_face: erroneous value of dir" << finl;
      Process::exit();
    }

  field_ft.echange_espace_virtuel(ghost);
}

// Combine les sous-elements deux-a-deux par une somme directe.
// Interdit une valeur non-nulle si le volume de controle ft est pure (selon l'indication de l'indicatrice).
// Utile pour la surface de l'interface.
void somme_sous_elements_si_indicatrice_elem(const IJK_Field_double& field_8x, IJK_Field_double& field_ft, const IJK_Field_double& indicatrice_ft)
{
  const int ni = field_ft.ni();
  const int nj = field_ft.nj();
  const int nk = field_ft.nk();
  const int ghost = field_ft.ghost();
  assert(field_8x.ghost() == 2*field_ft.ghost());
  for (int k = -ghost; k < nk+ghost; k++)
    {
      for (int j = -ghost; j < nj+ghost; j++)
        {
          for (int i = -ghost; i < ni+ghost; i++)
            {
              // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
              // qu'une cellule consideree pure soit traversee par l'interface.
              // Pour coherence, on considere les faces non coupees dans ce cas.
              if (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)))
                {
                  field_ft(i,j,k) =      (field_8x(2*i,2*j,2*k)       + field_8x(2*i+1,2*j,2*k)
                                          + field_8x(2*i,2*j+1,2*k)   + field_8x(2*i+1,2*j+1,2*k)
                                          + field_8x(2*i,2*j,2*k+1)   + field_8x(2*i+1,2*j,2*k+1)
                                          + field_8x(2*i,2*j+1,2*k+1) + field_8x(2*i+1,2*j+1,2*k+1));
                }
              else
                {
                  field_ft(i,j,k) = 0.;
                }
            }
        }
    }
}

// Combine les sous-elements deux-a-deux par une somme directe, sur un volume de controle ft decale d'une demi-maille dans une direction (c'est-a-dire centre sur une face).
// Interdit une valeur non-nulle si le volume de controle ft est pure (selon l'indication de l'indicatrice).
// Utile pour la surface de l'interface.
void somme_sous_elements_si_indicatrice_face(int dir, const IJK_Field_double& field_8x, IJK_Field_double& field_ft, const IJK_Field_double& indicatrice_ft)
{
  const int ni = field_ft.ni();
  const int nj = field_ft.nj();
  const int nk = field_ft.nk();
  const int ghost = field_ft.ghost();
  assert(field_8x.ghost() == 2*field_ft.ghost());
  if (dir == 0)
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                  // qu'une cellule consideree pure soit traversee par l'interface.
                  // Pour coherence, on considere les faces non coupees dans ce cas.
                  if (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)))
                    {
                      field_ft(i,j,k) =      (field_8x(2*i,2*j,2*k)       + field_8x(2*i-1,2*j,2*k)
                                              + field_8x(2*i,2*j+1,2*k)   + field_8x(2*i-1,2*j+1,2*k)
                                              + field_8x(2*i,2*j,2*k+1)   + field_8x(2*i-1,2*j,2*k+1)
                                              + field_8x(2*i,2*j+1,2*k+1) + field_8x(2*i-1,2*j+1,2*k+1));
                    }
                  else
                    {
                      field_ft(i,j,k) = 0.;
                    }
                }
            }
        }
    }
  else if (dir == 1)
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                  // qu'une cellule consideree pure soit traversee par l'interface.
                  // Pour coherence, on considere les faces non coupees dans ce cas.
                  if (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)))
                    {
                      field_ft(i,j,k) =      (field_8x(2*i,2*j,2*k)       + field_8x(2*i+1,2*j,2*k)
                                              + field_8x(2*i,2*j-1,2*k)   + field_8x(2*i+1,2*j-1,2*k)
                                              + field_8x(2*i,2*j,2*k+1)   + field_8x(2*i+1,2*j,2*k+1)
                                              + field_8x(2*i,2*j-1,2*k+1) + field_8x(2*i+1,2*j-1,2*k+1));
                    }
                  else
                    {
                      field_ft(i,j,k) = 0.;
                    }
                }
            }
        }
    }
  else if (dir == 2)
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                  // qu'une cellule consideree pure soit traversee par l'interface.
                  // Pour coherence, on considere les faces non coupees dans ce cas.
                  if (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)))
                    {
                      field_ft(i,j,k) =      (field_8x(2*i,2*j,2*k)       + field_8x(2*i+1,2*j,2*k)
                                              + field_8x(2*i,2*j+1,2*k)   + field_8x(2*i+1,2*j+1,2*k)
                                              + field_8x(2*i,2*j,2*k-1)   + field_8x(2*i+1,2*j,2*k-1)
                                              + field_8x(2*i,2*j+1,2*k-1) + field_8x(2*i+1,2*j+1,2*k-1));
                    }
                  else
                    {
                      field_ft(i,j,k) = 0.;
                    }
                }
            }
        }
    }
  else
    {
      Cerr << "somme_sous_elements_si_indicatrice_face: erroneous value of dir" << finl;
      Process::exit();
    }

  field_ft.echange_espace_virtuel(ghost);
}

// Combine les sous-faces deux-a-deux dans un element ft.
void moyenne_sous_faces(int face_dir, const IJK_Field_double& field_8x, IJK_Field_double& field_ft)
{
  const int ni = field_ft.ni();
  const int nj = field_ft.nj();
  const int nk = field_ft.nk();
  const int ghost = field_ft.ghost();
  assert(field_8x.ghost() == 2*field_ft.ghost());
  if (face_dir == 0)
    {
      for (int k = -ghost; k < nk+ghost; k++)
        {
          for (int j = -ghost; j < nj+ghost; j++)
            {
              for (int i = -ghost; i < ni+ghost; i++)
                {
                  field_ft(i,j,k) = .25*(field_8x(2*i,2*j,2*k) + field_8x(2*i,2*j+1,2*k)
                                         + field_8x(2*i,2*j,2*k+1) + field_8x(2*i,2*j+1,2*k+1));
                }
            }
        }
    }
  else if (face_dir == 1)
    {
      for (int k = -ghost; k < nk+ghost; k++)
        {
          for (int j = -ghost; j < nj+ghost; j++)
            {
              for (int i = -ghost; i < ni+ghost; i++)
                {
                  field_ft(i,j,k) = .25*(field_8x(2*i,2*j,2*k)       + field_8x(2*i+1,2*j,2*k)
                                         + field_8x(2*i,2*j,2*k+1)   + field_8x(2*i+1,2*j,2*k+1));
                }
            }
        }
    }
  else if (face_dir == 2)
    {
      for (int k = -ghost; k < nk+ghost; k++)
        {
          for (int j = -ghost; j < nj+ghost; j++)
            {
              for (int i = -ghost; i < ni+ghost; i++)
                {
                  field_ft(i,j,k) = .25*(field_8x(2*i,2*j,2*k)       + field_8x(2*i+1,2*j,2*k)
                                         + field_8x(2*i,2*j+1,2*k)   + field_8x(2*i+1,2*j+1,2*k));
                }
            }
        }
    }
  else
    {
      Cerr << "moyenne_sous_faces: erroneous value of face_dir" << finl;
      Process::exit();
    }
}

// Combine les sous-faces deux-a-deux dans un element ft.
// Utile pour l'indicatrice surfacique de la face : la fonction interdit une valeur non pure si l'un des deux volumes de controle contenant la face est pure.
void moyenne_indicatrice_surfacique_sous_faces_elem(int face_dir, const IJK_Field_double& field_8x, IJK_Field_double& field_ft, const IJK_Field_double& indicatrice_ft)
{
  const int ni = field_ft.ni();
  const int nj = field_ft.nj();
  const int nk = field_ft.nk();
  const int ghost = field_ft.ghost();
  assert(field_8x.ghost() == 2*field_ft.ghost());
  if (face_dir == 0)
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  double indicatrice_surfacique = .25*(field_8x(2*i,2*j,2*k) + field_8x(2*i,2*j+1,2*k)
                                                       + field_8x(2*i,2*j,2*k+1) + field_8x(2*i,2*j+1,2*k+1));
                  // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                  // qu'une cellule consideree pure soit traversee par l'interface.
                  // Pour coherence, on considere les faces non coupees dans ce cas.
                  //
                  // Une cellule adjacente pure suggere ou bien qu'il n'y a pas d'intersections
                  // ou bien qu'il y a des intersections mais que l'indicatrice resultante est
                  // trop faible et en dessous du seuil dans calculer_indicatrice.
                  // Pour coherence, on neglige la coupure de la surface egalement.
                  bool cellule_cote_positif_non_pure = (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)));
                  bool cellule_cote_negatif_non_pure = (!((indicatrice_ft(i-1, j, k) == 0.) || (indicatrice_ft(i-1, j, k) == 1.)));
                  if (cellule_cote_positif_non_pure && cellule_cote_negatif_non_pure)
                    {
                      field_ft(i,j,k) = indicatrice_surfacique;
                    }
                  else
                    {
                      field_ft(i,j,k) = (indicatrice_surfacique > .5) ? 1. : 0.;
                    }
                }
            }
        }
    }
  else if (face_dir == 1)
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  double indicatrice_surfacique = .25*(field_8x(2*i,2*j,2*k)       + field_8x(2*i+1,2*j,2*k)
                                                       + field_8x(2*i,2*j,2*k+1)   + field_8x(2*i+1,2*j,2*k+1));
                  // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                  // qu'une cellule consideree pure soit traversee par l'interface.
                  // Pour coherence, on considere les faces non coupees dans ce cas.
                  //
                  // Une cellule adjacente pure suggere ou bien qu'il n'y a pas d'intersections
                  // ou bien qu'il y a des intersections mais que l'indicatrice resultante est
                  // trop faible et en dessous du seuil dans calculer_indicatrice.
                  // Pour coherence, on neglige la coupure de la surface egalement.
                  bool cellule_cote_positif_non_pure = (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)));
                  bool cellule_cote_negatif_non_pure = (!((indicatrice_ft(i, j-1, k) == 0.) || (indicatrice_ft(i, j-1, k) == 1.)));
                  if (cellule_cote_positif_non_pure && cellule_cote_negatif_non_pure)
                    {
                      field_ft(i,j,k) = indicatrice_surfacique;
                    }
                  else
                    {
                      field_ft(i,j,k) = (indicatrice_surfacique > .5) ? 1. : 0.;
                    }
                }
            }
        }
    }
  else if (face_dir == 2)
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  double indicatrice_surfacique = .25*(field_8x(2*i,2*j,2*k)       + field_8x(2*i+1,2*j,2*k)
                                                       + field_8x(2*i,2*j+1,2*k)   + field_8x(2*i+1,2*j+1,2*k));
                  // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                  // qu'une cellule consideree pure soit traversee par l'interface.
                  // Pour coherence, on considere les faces non coupees dans ce cas.
                  //
                  // Une cellule adjacente pure suggere ou bien qu'il n'y a pas d'intersections
                  // ou bien qu'il y a des intersections mais que l'indicatrice resultante est
                  // trop faible et en dessous du seuil dans calculer_indicatrice.
                  // Pour coherence, on neglige la coupure de la surface egalement.
                  bool cellule_cote_positif_non_pure = (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)));
                  bool cellule_cote_negatif_non_pure = (!((indicatrice_ft(i, j, k-1) == 0.) || (indicatrice_ft(i, j, k-1) == 1.)));
                  if (cellule_cote_positif_non_pure && cellule_cote_negatif_non_pure)
                    {
                      field_ft(i,j,k) = indicatrice_surfacique;
                    }
                  else
                    {
                      field_ft(i,j,k) = (indicatrice_surfacique > .5) ? 1. : 0.;
                    }
                }
            }
        }
    }
  else
    {
      Cerr << "moyenne_indicatrice_surfacique_sous_faces_elem: erroneous value of face_dir" << finl;
      Process::exit();
    }

  field_ft.echange_espace_virtuel(ghost);
}

// Combine les sous-faces deux-a-deux dans un element ft, sur un volume de controle ft decale d'une demi-maille dans une direction (c'est-a-dire centre sur une face).
// Utile pour l'indicatrice surfacique de la face : la fonction interdit une valeur non pure si l'un des deux volumes de controle contenant la face est pure.
void moyenne_indicatrice_surfacique_sous_faces_face(int dir, int face_dir, const IJK_Field_double& field_8x, IJK_Field_double& field_ft, const IJK_Field_double& indicatrice_ft)
{
  const int ni = field_ft.ni();
  const int nj = field_ft.nj();
  const int nk = field_ft.nk();
  const int ghost = field_ft.ghost();
  assert(field_8x.ghost() == 2*field_ft.ghost());
  if (dir == 0)
    {
      if (face_dir == 0)
        {
          for (int k = 0; k < nk; k++)
            {
              for (int j = 0; j < nj; j++)
                {
                  for (int i = 0; i < ni; i++)
                    {
                      double indicatrice_surfacique = .25*(field_8x(2*i-1,2*j,2*k) + field_8x(2*i-1,2*j+1,2*k)
                                                           + field_8x(2*i-1,2*j,2*k+1) + field_8x(2*i-1,2*j+1,2*k+1));
                      // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                      // qu'une cellule consideree pure soit traversee par l'interface.
                      // Pour coherence, on considere les faces non coupees dans ce cas.
                      //
                      // Une cellule adjacente pure suggere ou bien qu'il n'y a pas d'intersections
                      // ou bien qu'il y a des intersections mais que l'indicatrice resultante est
                      // trop faible et en dessous du seuil dans calculer_indicatrice.
                      // Pour coherence, on neglige la coupure de la surface egalement.
                      bool cellule_cote_positif_non_pure = (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)));
                      bool cellule_cote_negatif_non_pure = (!((indicatrice_ft(i-1, j, k) == 0.) || (indicatrice_ft(i-1, j, k) == 1.)));
                      if (cellule_cote_positif_non_pure && cellule_cote_negatif_non_pure)
                        {
                          field_ft(i,j,k) = indicatrice_surfacique;
                        }
                      else
                        {
                          field_ft(i,j,k) = (indicatrice_surfacique > .5) ? 1. : 0.;
                        }
                    }
                }
            }
        }
      else if (face_dir == 1)
        {
          for (int k = 0; k < nk; k++)
            {
              for (int j = 0; j < nj; j++)
                {
                  for (int i = 0; i < ni; i++)
                    {
                      double indicatrice_surfacique = .25*(field_8x(2*i-1,2*j,2*k)       + field_8x(2*i,2*j,2*k)
                                                           + field_8x(2*i-1,2*j,2*k+1)   + field_8x(2*i,2*j,2*k+1));
                      // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                      // qu'une cellule consideree pure soit traversee par l'interface.
                      // Pour coherence, on considere les faces non coupees dans ce cas.
                      //
                      // Une cellule adjacente pure suggere ou bien qu'il n'y a pas d'intersections
                      // ou bien qu'il y a des intersections mais que l'indicatrice resultante est
                      // trop faible et en dessous du seuil dans calculer_indicatrice.
                      // Pour coherence, on neglige la coupure de la surface egalement.
                      bool cellule_cote_positif_non_pure = (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)));
                      bool cellule_cote_negatif_non_pure = (!((indicatrice_ft(i, j-1, k) == 0.) || (indicatrice_ft(i, j-1, k) == 1.)));
                      if (cellule_cote_positif_non_pure && cellule_cote_negatif_non_pure)
                        {
                          field_ft(i,j,k) = indicatrice_surfacique;
                        }
                      else
                        {
                          field_ft(i,j,k) = (indicatrice_surfacique > .5) ? 1. : 0.;
                        }
                    }
                }
            }
        }
      else if (face_dir == 2)
        {
          for (int k = 0; k < nk; k++)
            {
              for (int j = 0; j < nj; j++)
                {
                  for (int i = 0; i < ni; i++)
                    {
                      double indicatrice_surfacique = .25*(field_8x(2*i-1,2*j,2*k)       + field_8x(2*i,2*j,2*k)
                                                           + field_8x(2*i-1,2*j+1,2*k)   + field_8x(2*i,2*j+1,2*k));
                      // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                      // qu'une cellule consideree pure soit traversee par l'interface.
                      // Pour coherence, on considere les faces non coupees dans ce cas.
                      //
                      // Une cellule adjacente pure suggere ou bien qu'il n'y a pas d'intersections
                      // ou bien qu'il y a des intersections mais que l'indicatrice resultante est
                      // trop faible et en dessous du seuil dans calculer_indicatrice.
                      // Pour coherence, on neglige la coupure de la surface egalement.
                      bool cellule_cote_positif_non_pure = (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)));
                      bool cellule_cote_negatif_non_pure = (!((indicatrice_ft(i, j, k-1) == 0.) || (indicatrice_ft(i, j, k-1) == 1.)));
                      if (cellule_cote_positif_non_pure && cellule_cote_negatif_non_pure)
                        {
                          field_ft(i,j,k) = indicatrice_surfacique;
                        }
                      else
                        {
                          field_ft(i,j,k) = (indicatrice_surfacique > .5) ? 1. : 0.;
                        }
                    }
                }
            }
        }
      else
        {
          Cerr << "moyenne_indicatrice_surfacique_sous_faces_face: erroneous value of face_dir" << finl;
          Process::exit();
        }
    }
  else if (dir == 1)
    {
      if (face_dir == 0)
        {
          for (int k = 0; k < nk; k++)
            {
              for (int j = 0; j < nj; j++)
                {
                  for (int i = 0; i < ni; i++)
                    {
                      double indicatrice_surfacique = .25*(field_8x(2*i,2*j-1,2*k) + field_8x(2*i,2*j,2*k)
                                                           + field_8x(2*i,2*j-1,2*k+1) + field_8x(2*i,2*j,2*k+1));
                      // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                      // qu'une cellule consideree pure soit traversee par l'interface.
                      // Pour coherence, on considere les faces non coupees dans ce cas.
                      //
                      // Une cellule adjacente pure suggere ou bien qu'il n'y a pas d'intersections
                      // ou bien qu'il y a des intersections mais que l'indicatrice resultante est
                      // trop faible et en dessous du seuil dans calculer_indicatrice.
                      // Pour coherence, on neglige la coupure de la surface egalement.
                      bool cellule_cote_positif_non_pure = (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)));
                      bool cellule_cote_negatif_non_pure = (!((indicatrice_ft(i-1, j, k) == 0.) || (indicatrice_ft(i-1, j, k) == 1.)));
                      if (cellule_cote_positif_non_pure && cellule_cote_negatif_non_pure)
                        {
                          field_ft(i,j,k) = indicatrice_surfacique;
                        }
                      else
                        {
                          field_ft(i,j,k) = (indicatrice_surfacique > .5) ? 1. : 0.;
                        }
                    }
                }
            }
        }
      else if (face_dir == 1)
        {
          for (int k = 0; k < nk; k++)
            {
              for (int j = 0; j < nj; j++)
                {
                  for (int i = 0; i < ni; i++)
                    {
                      double indicatrice_surfacique = .25*(field_8x(2*i,2*j-1,2*k)       + field_8x(2*i+1,2*j-1,2*k)
                                                           + field_8x(2*i,2*j-1,2*k+1)   + field_8x(2*i+1,2*j-1,2*k+1));
                      // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                      // qu'une cellule consideree pure soit traversee par l'interface.
                      // Pour coherence, on considere les faces non coupees dans ce cas.
                      //
                      // Une cellule adjacente pure suggere ou bien qu'il n'y a pas d'intersections
                      // ou bien qu'il y a des intersections mais que l'indicatrice resultante est
                      // trop faible et en dessous du seuil dans calculer_indicatrice.
                      // Pour coherence, on neglige la coupure de la surface egalement.
                      bool cellule_cote_positif_non_pure = (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)));
                      bool cellule_cote_negatif_non_pure = (!((indicatrice_ft(i, j-1, k) == 0.) || (indicatrice_ft(i, j-1, k) == 1.)));
                      if (cellule_cote_positif_non_pure && cellule_cote_negatif_non_pure)
                        {
                          field_ft(i,j,k) = indicatrice_surfacique;
                        }
                      else
                        {
                          field_ft(i,j,k) = (indicatrice_surfacique > .5) ? 1. : 0.;
                        }
                    }
                }
            }
        }
      else if (face_dir == 2)
        {
          for (int k = 0; k < nk; k++)
            {
              for (int j = 0; j < nj; j++)
                {
                  for (int i = 0; i < ni; i++)
                    {
                      double indicatrice_surfacique = .25*(field_8x(2*i,2*j-1,2*k)       + field_8x(2*i+1,2*j-1,2*k)
                                                           + field_8x(2*i,2*j,2*k)   + field_8x(2*i+1,2*j,2*k));
                      // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                      // qu'une cellule consideree pure soit traversee par l'interface.
                      // Pour coherence, on considere les faces non coupees dans ce cas.
                      //
                      // Une cellule adjacente pure suggere ou bien qu'il n'y a pas d'intersections
                      // ou bien qu'il y a des intersections mais que l'indicatrice resultante est
                      // trop faible et en dessous du seuil dans calculer_indicatrice.
                      // Pour coherence, on neglige la coupure de la surface egalement.
                      bool cellule_cote_positif_non_pure = (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)));
                      bool cellule_cote_negatif_non_pure = (!((indicatrice_ft(i, j, k-1) == 0.) || (indicatrice_ft(i, j, k-1) == 1.)));
                      if (cellule_cote_positif_non_pure && cellule_cote_negatif_non_pure)
                        {
                          field_ft(i,j,k) = indicatrice_surfacique;
                        }
                      else
                        {
                          field_ft(i,j,k) = (indicatrice_surfacique > .5) ? 1. : 0.;
                        }
                    }
                }
            }
        }
      else
        {
          Cerr << "moyenne_indicatrice_surfacique_sous_faces_face: erroneous value of face_dir" << finl;
          Process::exit();
        }
    }
  else if (dir == 2)
    {
      if (face_dir == 0)
        {
          for (int k = 0; k < nk; k++)
            {
              for (int j = 0; j < nj; j++)
                {
                  for (int i = 0; i < ni; i++)
                    {
                      double indicatrice_surfacique = .25*(field_8x(2*i,2*j,2*k-1) + field_8x(2*i,2*j+1,2*k-1)
                                                           + field_8x(2*i,2*j,2*k) + field_8x(2*i,2*j+1,2*k));
                      // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                      // qu'une cellule consideree pure soit traversee par l'interface.
                      // Pour coherence, on considere les faces non coupees dans ce cas.
                      //
                      // Une cellule adjacente pure suggere ou bien qu'il n'y a pas d'intersections
                      // ou bien qu'il y a des intersections mais que l'indicatrice resultante est
                      // trop faible et en dessous du seuil dans calculer_indicatrice.
                      // Pour coherence, on neglige la coupure de la surface egalement.
                      bool cellule_cote_positif_non_pure = (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)));
                      bool cellule_cote_negatif_non_pure = (!((indicatrice_ft(i-1, j, k) == 0.) || (indicatrice_ft(i-1, j, k) == 1.)));
                      if (cellule_cote_positif_non_pure && cellule_cote_negatif_non_pure)
                        {
                          field_ft(i,j,k) = indicatrice_surfacique;
                        }
                      else
                        {
                          field_ft(i,j,k) = (indicatrice_surfacique > .5) ? 1. : 0.;
                        }
                    }
                }
            }
        }
      else if (face_dir == 1)
        {
          for (int k = 0; k < nk; k++)
            {
              for (int j = 0; j < nj; j++)
                {
                  for (int i = 0; i < ni; i++)
                    {
                      double indicatrice_surfacique = .25*(field_8x(2*i,2*j,2*k-1)       + field_8x(2*i+1,2*j,2*k-1)
                                                           + field_8x(2*i,2*j,2*k)   + field_8x(2*i+1,2*j,2*k));
                      // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                      // qu'une cellule consideree pure soit traversee par l'interface.
                      // Pour coherence, on considere les faces non coupees dans ce cas.
                      //
                      // Une cellule adjacente pure suggere ou bien qu'il n'y a pas d'intersections
                      // ou bien qu'il y a des intersections mais que l'indicatrice resultante est
                      // trop faible et en dessous du seuil dans calculer_indicatrice.
                      // Pour coherence, on neglige la coupure de la surface egalement.
                      bool cellule_cote_positif_non_pure = (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)));
                      bool cellule_cote_negatif_non_pure = (!((indicatrice_ft(i, j-1, k) == 0.) || (indicatrice_ft(i, j-1, k) == 1.)));
                      if (cellule_cote_positif_non_pure && cellule_cote_negatif_non_pure)
                        {
                          field_ft(i,j,k) = indicatrice_surfacique;
                        }
                      else
                        {
                          field_ft(i,j,k) = (indicatrice_surfacique > .5) ? 1. : 0.;
                        }
                    }
                }
            }
        }
      else if (face_dir == 2)
        {
          for (int k = 0; k < nk; k++)
            {
              for (int j = 0; j < nj; j++)
                {
                  for (int i = 0; i < ni; i++)
                    {
                      double indicatrice_surfacique = .25*(field_8x(2*i,2*j,2*k-1)       + field_8x(2*i+1,2*j,2*k-1)
                                                           + field_8x(2*i,2*j+1,2*k-1)   + field_8x(2*i+1,2*j+1,2*k-1));
                      // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                      // qu'une cellule consideree pure soit traversee par l'interface.
                      // Pour coherence, on considere les faces non coupees dans ce cas.
                      //
                      // Une cellule adjacente pure suggere ou bien qu'il n'y a pas d'intersections
                      // ou bien qu'il y a des intersections mais que l'indicatrice resultante est
                      // trop faible et en dessous du seuil dans calculer_indicatrice.
                      // Pour coherence, on neglige la coupure de la surface egalement.
                      bool cellule_cote_positif_non_pure = (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)));
                      bool cellule_cote_negatif_non_pure = (!((indicatrice_ft(i, j, k-1) == 0.) || (indicatrice_ft(i, j, k-1) == 1.)));
                      if (cellule_cote_positif_non_pure && cellule_cote_negatif_non_pure)
                        {
                          field_ft(i,j,k) = indicatrice_surfacique;
                        }
                      else
                        {
                          field_ft(i,j,k) = (indicatrice_surfacique > .5) ? 1. : 0.;
                        }
                    }
                }
            }
        }
      else
        {
          Cerr << "moyenne_indicatrice_surfacique_sous_faces_face: erroneous value of face_dir" << finl;
          Process::exit();
        }
    }
  else
    {
      Cerr << "moyenne_indicatrice_surfacique_sous_faces_face: erroneous value of dir" << finl;
      Process::exit();
    }

  field_ft.echange_espace_virtuel(ghost);
}

// Combine le barycentre des sous-elements dans un element ft.
void moyenne_barycentre_phase1_sous_elements_elem(int bary_compo, const IJK_Field_double& indic_8x, const IJK_Field_double& field_8x, IJK_Field_double& field_ft, const IJK_Field_double& indicatrice_ft)
{
  const int ni = field_ft.ni();
  const int nj = field_ft.nj();
  const int nk = field_ft.nk();
  const int ghost = field_ft.ghost();
  assert(field_8x.ghost() == 2*field_ft.ghost());
  for (int k = -ghost; k < nk+ghost; k++)
    {
      for (int j = -ghost; j < nj+ghost; j++)
        {
          for (int i = -ghost; i < ni+ghost; i++)
            {
              // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
              // qu'une cellule consideree pure soit traversee par l'interface.
              // Pour coherence, on considere les faces non coupees dans ce cas.
              if (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)))
                {
                  double indic = .125*(indic_8x(2*i,2*j,2*k)       + indic_8x(2*i+1,2*j,2*k)
                                       + indic_8x(2*i,2*j+1,2*k)   + indic_8x(2*i+1,2*j+1,2*k)
                                       + indic_8x(2*i,2*j,2*k+1)   + indic_8x(2*i+1,2*j,2*k+1)
                                       + indic_8x(2*i,2*j+1,2*k+1) + indic_8x(2*i+1,2*j+1,2*k+1));
                  field_ft(i,j,k) = (indic == 0 || indic == 1) ? .5 :
                                    (  indic_8x(2*i,2*j,2*k)       * (.5*field_8x(2*i,2*j,2*k)                                      ) + indic_8x(2*i+1,2*j,2*k)     * (.5*field_8x(2*i+1,2*j,2*k)     + (bary_compo==0                    )*.5)
                                       + indic_8x(2*i,2*j+1,2*k)   * (.5*field_8x(2*i,2*j+1,2*k)   + (          bary_compo==1          )*.5) + indic_8x(2*i+1,2*j+1,2*k)   * (.5*field_8x(2*i+1,2*j+1,2*k)   + (bary_compo==0 || bary_compo==1          )*.5)
                                       + indic_8x(2*i,2*j,2*k+1)   * (.5*field_8x(2*i,2*j,2*k+1)   + (                    bary_compo==2)*.5) + indic_8x(2*i+1,2*j,2*k+1)   * (.5*field_8x(2*i+1,2*j,2*k+1)   + (bary_compo==0 ||           bary_compo==2)*.5)
                                       + indic_8x(2*i,2*j+1,2*k+1) * (.5*field_8x(2*i,2*j+1,2*k+1) + (          bary_compo==1 || bary_compo==2)*.5) + indic_8x(2*i+1,2*j+1,2*k+1) * (.5*field_8x(2*i+1,2*j+1,2*k+1) + (bary_compo==0 || bary_compo==1 || bary_compo==2)*.5))/(8*indic);
                }
              else
                {
                  field_ft(i,j,k) = .5;
                }
            }
        }
    }
}

// Combine le barycentre des sous-elements dans un element ft, sur un volume de controle ft decale d'une demi-maille dans une direction (c'est-a-dire centre sur une face).
void moyenne_barycentre_phase1_sous_elements_face(int dir, int bary_compo, const IJK_Field_double& indic_8x, const IJK_Field_double& field_8x, IJK_Field_double& field_ft, const IJK_Field_double& indicatrice_ft)
{
  const int ni = field_ft.ni();
  const int nj = field_ft.nj();
  const int nk = field_ft.nk();
  const int ghost = field_ft.ghost();
  assert(field_8x.ghost() == 2*field_ft.ghost());
  if (dir == 0)
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                  // qu'une cellule consideree pure soit traversee par l'interface.
                  // Pour coherence, on considere les faces non coupees dans ce cas.
                  if (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)))
                    {
                      double indic = .125*(indic_8x(2*i,2*j,2*k)       + indic_8x(2*i-1,2*j,2*k)
                                           + indic_8x(2*i,2*j+1,2*k)   + indic_8x(2*i-1,2*j+1,2*k)
                                           + indic_8x(2*i,2*j,2*k+1)   + indic_8x(2*i-1,2*j,2*k+1)
                                           + indic_8x(2*i,2*j+1,2*k+1) + indic_8x(2*i-1,2*j+1,2*k+1));
                      field_ft(i,j,k) = (indic == 0 || indic == 1) ? .5 :
                                        (  indic_8x(2*i,2*j,2*k)       * (.5*field_8x(2*i,2*j,2*k)                                             ) + indic_8x(2*i-1,2*j,2*k)     * (.5*field_8x(2*i-1,2*j,2*k)     + (bary_compo==0                                  )*.5)
                                           + indic_8x(2*i,2*j+1,2*k)   * (.5*field_8x(2*i,2*j+1,2*k)   + (   bary_compo==1                 )*.5) + indic_8x(2*i-1,2*j+1,2*k)   * (.5*field_8x(2*i-1,2*j+1,2*k)   + (bary_compo==0 || bary_compo==1                 )*.5)
                                           + indic_8x(2*i,2*j,2*k+1)   * (.5*field_8x(2*i,2*j,2*k+1)   + (                    bary_compo==2)*.5) + indic_8x(2*i-1,2*j,2*k+1)   * (.5*field_8x(2*i-1,2*j,2*k+1)   + (bary_compo==0 ||                  bary_compo==2)*.5)
                                           + indic_8x(2*i,2*j+1,2*k+1) * (.5*field_8x(2*i,2*j+1,2*k+1) + (   bary_compo==1 || bary_compo==2)*.5) + indic_8x(2*i-1,2*j+1,2*k+1) * (.5*field_8x(2*i-1,2*j+1,2*k+1) + (bary_compo==0 || bary_compo==1 || bary_compo==2)*.5))/(8*indic);
                    }
                  else
                    {
                      field_ft(i,j,k) = .5;
                    }
                }
            }
        }
    }
  else if (dir == 1)
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                  // qu'une cellule consideree pure soit traversee par l'interface.
                  // Pour coherence, on considere les faces non coupees dans ce cas.
                  if (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)))
                    {
                      double indic = .125*(indic_8x(2*i,2*j,2*k)       + indic_8x(2*i+1,2*j,2*k)
                                           + indic_8x(2*i,2*j-1,2*k)   + indic_8x(2*i+1,2*j-1,2*k)
                                           + indic_8x(2*i,2*j,2*k+1)   + indic_8x(2*i+1,2*j,2*k+1)
                                           + indic_8x(2*i,2*j-1,2*k+1) + indic_8x(2*i+1,2*j-1,2*k+1));
                      field_ft(i,j,k) = (indic == 0 || indic == 1) ? .5 :
                                        (  indic_8x(2*i,2*j,2*k)       * (.5*field_8x(2*i,2*j,2*k)                                             ) + indic_8x(2*i+1,2*j,2*k)     * (.5*field_8x(2*i+1,2*j,2*k)     + (bary_compo==0                                  )*.5)
                                           + indic_8x(2*i,2*j-1,2*k)   * (.5*field_8x(2*i,2*j-1,2*k)   + (   bary_compo==1                 )*.5) + indic_8x(2*i+1,2*j-1,2*k)   * (.5*field_8x(2*i+1,2*j-1,2*k)   + (bary_compo==0 || bary_compo==1                 )*.5)
                                           + indic_8x(2*i,2*j,2*k+1)   * (.5*field_8x(2*i,2*j,2*k+1)   + (                    bary_compo==2)*.5) + indic_8x(2*i+1,2*j,2*k+1)   * (.5*field_8x(2*i+1,2*j,2*k+1)   + (bary_compo==0 ||                  bary_compo==2)*.5)
                                           + indic_8x(2*i,2*j-1,2*k+1) * (.5*field_8x(2*i,2*j-1,2*k+1) + (   bary_compo==1 || bary_compo==2)*.5) + indic_8x(2*i+1,2*j-1,2*k+1) * (.5*field_8x(2*i+1,2*j-1,2*k+1) + (bary_compo==0 || bary_compo==1 || bary_compo==2)*.5))/(8*indic);
                    }
                  else
                    {
                      field_ft(i,j,k) = .5;
                    }
                }
            }
        }
    }
  else if (dir == 2)
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                  // qu'une cellule consideree pure soit traversee par l'interface.
                  // Pour coherence, on considere les faces non coupees dans ce cas.
                  if (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)))
                    {
                      double indic = .125*(indic_8x(2*i,2*j,2*k)       + indic_8x(2*i+1,2*j,2*k)
                                           + indic_8x(2*i,2*j+1,2*k)   + indic_8x(2*i+1,2*j+1,2*k)
                                           + indic_8x(2*i,2*j,2*k-1)   + indic_8x(2*i+1,2*j,2*k-1)
                                           + indic_8x(2*i,2*j+1,2*k-1) + indic_8x(2*i+1,2*j+1,2*k-1));
                      field_ft(i,j,k) = (indic == 0 || indic == 1) ? .5 :
                                        (  indic_8x(2*i,2*j,2*k)       * (.5*field_8x(2*i,2*j,2*k)                                             ) + indic_8x(2*i+1,2*j,2*k)     * (.5*field_8x(2*i+1,2*j,2*k)     + (bary_compo==0                                  )*.5)
                                           + indic_8x(2*i,2*j+1,2*k)   * (.5*field_8x(2*i,2*j+1,2*k)   + (   bary_compo==1                 )*.5) + indic_8x(2*i+1,2*j+1,2*k)   * (.5*field_8x(2*i+1,2*j+1,2*k)   + (bary_compo==0 || bary_compo==1                 )*.5)
                                           + indic_8x(2*i,2*j,2*k-1)   * (.5*field_8x(2*i,2*j,2*k-1)   + (                    bary_compo==2)*.5) + indic_8x(2*i+1,2*j,2*k-1)   * (.5*field_8x(2*i+1,2*j,2*k-1)   + (bary_compo==0 ||                  bary_compo==2)*.5)
                                           + indic_8x(2*i,2*j+1,2*k-1) * (.5*field_8x(2*i,2*j+1,2*k-1) + (   bary_compo==1 || bary_compo==2)*.5) + indic_8x(2*i+1,2*j+1,2*k-1) * (.5*field_8x(2*i+1,2*j+1,2*k-1) + (bary_compo==0 || bary_compo==1 || bary_compo==2)*.5))/(8*indic);
                    }
                  else
                    {
                      field_ft(i,j,k) = .5;
                    }
                }
            }
        }
    }
  else
    {
      Cerr << "moyenne_barycentre_phase1_sous_elements_face: erroneous value of dir" << finl;
      Process::exit();
    }

  field_ft.echange_espace_virtuel(ghost);
}

// Combine le barycentre surfacique des sous-faces dans un element ft.
void moyenne_barycentre_phase1_sous_faces_elem(int face_dir, int bary_compo, const IJK_Field_double& indic_8x, const IJK_Field_double& field_8x, IJK_Field_double& field_ft, const IJK_Field_double& indicatrice_ft)
{
  // Pour la face x      dir0=z -> (2->0)   dir1=y -> (1->1)
  // Pour la face y      dir0=x -> (0->0)   dir1=z -> (2->1)
  // Pour la face z      dir0=y -> (1->0)   dir1=x -> (0->1)

  const int ni = field_ft.ni();
  const int nj = field_ft.nj();
  const int nk = field_ft.nk();
  const int ghost = field_ft.ghost();
  assert(field_8x.ghost() == 2*field_ft.ghost());
  if (face_dir == 0)
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                  // qu'une cellule consideree pure soit traversee par l'interface.
                  // Pour coherence, on considere les faces non coupees dans ce cas.
                  //
                  // Une cellule adjacente pure suggere ou bien qu'il n'y a pas d'intersections
                  // ou bien qu'il y a des intersections mais que l'indicatrice resultante est
                  // trop faible et en dessous du seuil dans calculer_indicatrice.
                  // Pour coherence, on neglige la coupure de la surface egalement.
                  bool cellule_cote_positif_non_pure = (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)));
                  bool cellule_cote_negatif_non_pure = (!((indicatrice_ft(i-1, j, k) == 0.) || (indicatrice_ft(i-1, j, k) == 1.)));
                  if (cellule_cote_positif_non_pure && cellule_cote_negatif_non_pure)
                    {
                      double indic = .25*(indic_8x(2*i,2*j,2*k) + indic_8x(2*i,2*j+1,2*k)
                                          + indic_8x(2*i,2*j,2*k+1) + indic_8x(2*i,2*j+1,2*k+1));
                      field_ft(i,j,k) = (indic == 0 || indic == 1) ? .5 :
                                        (  indic_8x(2*i,2*j,2*k)   * (.5*field_8x(2*i,2*j,2*k)                                             ) + indic_8x(2*i,2*j+1,2*k)   * (.5*field_8x(2*i,2*j+1,2*k)   + (          bary_compo==1                 )*.5)
                                           + indic_8x(2*i,2*j,2*k+1) * (.5*field_8x(2*i,2*j,2*k+1) + (                    bary_compo==0)*.5) + indic_8x(2*i,2*j+1,2*k+1) * (.5*field_8x(2*i,2*j+1,2*k+1) + (          bary_compo==1 || bary_compo==0)*.5))/(4*indic);
                    }
                  else
                    {
                      field_ft(i,j,k) = .5;
                    }
                }
            }
        }
    }
  else if (face_dir == 1)
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                  // qu'une cellule consideree pure soit traversee par l'interface.
                  // Pour coherence, on considere les faces non coupees dans ce cas.
                  //
                  // Une cellule adjacente pure suggere ou bien qu'il n'y a pas d'intersections
                  // ou bien qu'il y a des intersections mais que l'indicatrice resultante est
                  // trop faible et en dessous du seuil dans calculer_indicatrice.
                  // Pour coherence, on neglige la coupure de la surface egalement.
                  bool cellule_cote_positif_non_pure = (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)));
                  bool cellule_cote_negatif_non_pure = (!((indicatrice_ft(i, j-1, k) == 0.) || (indicatrice_ft(i, j-1, k) == 1.)));
                  if (cellule_cote_positif_non_pure && cellule_cote_negatif_non_pure)
                    {
                      double indic = .25*(indic_8x(2*i,2*j,2*k)       + indic_8x(2*i+1,2*j,2*k)
                                          + indic_8x(2*i,2*j,2*k+1)   + indic_8x(2*i+1,2*j,2*k+1));
                      field_ft(i,j,k) = (indic == 0 || indic == 1) ? .5 :
                                        (  indic_8x(2*i,2*j,2*k)   * (.5*field_8x(2*i,2*j,2*k)                                             )   + indic_8x(2*i+1,2*j,2*k)   * (.5*field_8x(2*i+1,2*j,2*k)   + (bary_compo==0                           )*.5)
                                           + indic_8x(2*i,2*j,2*k+1) * (.5*field_8x(2*i,2*j,2*k+1) + (                    bary_compo==1)*.5)   + indic_8x(2*i+1,2*j,2*k+1) * (.5*field_8x(2*i+1,2*j,2*k+1) + (bary_compo==0           || bary_compo==1)*.5))/(4*indic);
                    }
                  else
                    {
                      field_ft(i,j,k) = .5;
                    }
                }
            }
        }
    }
  else if (face_dir == 2)
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                  // qu'une cellule consideree pure soit traversee par l'interface.
                  // Pour coherence, on considere les faces non coupees dans ce cas.
                  //
                  // Une cellule adjacente pure suggere ou bien qu'il n'y a pas d'intersections
                  // ou bien qu'il y a des intersections mais que l'indicatrice resultante est
                  // trop faible et en dessous du seuil dans calculer_indicatrice.
                  // Pour coherence, on neglige la coupure de la surface egalement.
                  bool cellule_cote_positif_non_pure = (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)));
                  bool cellule_cote_negatif_non_pure = (!((indicatrice_ft(i, j, k-1) == 0.) || (indicatrice_ft(i, j, k-1) == 1.)));
                  if (cellule_cote_positif_non_pure && cellule_cote_negatif_non_pure)
                    {
                      double indic = .25*(indic_8x(2*i,2*j,2*k)       + indic_8x(2*i+1,2*j,2*k)
                                          + indic_8x(2*i,2*j+1,2*k)   + indic_8x(2*i+1,2*j+1,2*k));
                      field_ft(i,j,k) = (indic == 0 || indic == 1) ? .5 :
                                        (  indic_8x(2*i,2*j,2*k)   * (.5*field_8x(2*i,2*j,2*k)                                             )   + indic_8x(2*i+1,2*j,2*k)   * (.5*field_8x(2*i+1,2*j,2*k)   + (bary_compo==1                           )*.5)
                                           + indic_8x(2*i,2*j+1,2*k) * (.5*field_8x(2*i,2*j+1,2*k) + (          bary_compo==0          )*.5)   + indic_8x(2*i+1,2*j+1,2*k) * (.5*field_8x(2*i+1,2*j+1,2*k) + (bary_compo==1 || bary_compo==0          )*.5))/(4*indic);
                    }
                  else
                    {
                      field_ft(i,j,k) = .5;
                    }
                }
            }
        }
    }
  else
    {
      Cerr << "moyenne_barycentre_phase1_sous_faces_elem: erroneous value of face_dir" << finl;
      Process::exit();
    }

  field_ft.echange_espace_virtuel(ghost);
}

// Combine le barycentre surfacique des sous-faces dans un element ft, sur un volume de controle ft decale d'une demi-maille dans une direction (c'est-a-dire centre sur une face).
void moyenne_barycentre_phase1_sous_faces_face(int dir, int face_dir, int bary_compo, const IJK_Field_double& indic_8x, const IJK_Field_double& field_8x, IJK_Field_double& field_ft, const IJK_Field_double& indicatrice_ft)
{
  // Pour la face x      dir0=z -> (2->0)   dir1=y -> (1->1)
  // Pour la face y      dir0=x -> (0->0)   dir1=z -> (2->1)
  // Pour la face z      dir0=y -> (1->0)   dir1=x -> (0->1)

  const int ni = field_ft.ni();
  const int nj = field_ft.nj();
  const int nk = field_ft.nk();
  const int ghost = field_ft.ghost();
  assert(field_8x.ghost() == 2*field_ft.ghost());
  if (face_dir == 0)
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                  // qu'une cellule consideree pure soit traversee par l'interface.
                  // Pour coherence, on considere les faces non coupees dans ce cas.
                  //
                  // Une cellule adjacente pure suggere ou bien qu'il n'y a pas d'intersections
                  // ou bien qu'il y a des intersections mais que l'indicatrice resultante est
                  // trop faible et en dessous du seuil dans calculer_indicatrice.
                  // Pour coherence, on neglige la coupure de la surface egalement.
                  bool cellule_cote_positif_non_pure = (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)));
                  bool cellule_cote_negatif_non_pure = (!((indicatrice_ft(i-1, j, k) == 0.) || (indicatrice_ft(i-1, j, k) == 1.)));
                  if (cellule_cote_positif_non_pure && cellule_cote_negatif_non_pure)
                    {
                      double indic = .25*(indic_8x(2*i,2*j,2*k) + indic_8x(2*i,2*j+1,2*k)
                                          + indic_8x(2*i,2*j,2*k+1) + indic_8x(2*i,2*j+1,2*k+1));
                      field_ft(i,j,k) = (indic == 0 || indic == 1) ? .5 :
                                        (  indic_8x(2*i,2*j,2*k)   * (.5*field_8x(2*i,2*j,2*k)                                    ) + indic_8x(2*i,2*j+1,2*k)   * (.5*field_8x(2*i,2*j+1,2*k)    + (          bary_compo==1          )*.5)
                                           + indic_8x(2*i,2*j,2*k+1) * (.5*field_8x(2*i,2*j,2*k+1) + (                    bary_compo==0)*.5) + indic_8x(2*i,2*j+1,2*k+1) * (.5*field_8x(2*i,2*j+1,2*k+1) + (          bary_compo==1 || bary_compo==0)*.5))/(4*indic);
                    }
                  else
                    {
                      field_ft(i,j,k) = .5;
                    }
                }
            }
        }
    }
  else if (face_dir == 1)
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                  // qu'une cellule consideree pure soit traversee par l'interface.
                  // Pour coherence, on considere les faces non coupees dans ce cas.
                  //
                  // Une cellule adjacente pure suggere ou bien qu'il n'y a pas d'intersections
                  // ou bien qu'il y a des intersections mais que l'indicatrice resultante est
                  // trop faible et en dessous du seuil dans calculer_indicatrice.
                  // Pour coherence, on neglige la coupure de la surface egalement.
                  bool cellule_cote_positif_non_pure = (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)));
                  bool cellule_cote_negatif_non_pure = (!((indicatrice_ft(i, j-1, k) == 0.) || (indicatrice_ft(i, j-1, k) == 1.)));
                  if (cellule_cote_positif_non_pure && cellule_cote_negatif_non_pure)
                    {
                      double indic = .25*(indic_8x(2*i,2*j,2*k)       + indic_8x(2*i+1,2*j,2*k)
                                          + indic_8x(2*i,2*j,2*k+1)   + indic_8x(2*i+1,2*j,2*k+1));
                      field_ft(i,j,k) = (indic == 0 || indic == 1) ? .5 :
                                        (  indic_8x(2*i,2*j,2*k)   * (.5*field_8x(2*i,2*j,2*k)                                    )   + indic_8x(2*i+1,2*j,2*k)   * (.5*field_8x(2*i+1,2*j,2*k)    + (bary_compo==0                    )*.5)
                                           + indic_8x(2*i,2*j,2*k+1) * (.5*field_8x(2*i,2*j,2*k+1) + (                    bary_compo==1)*.5)   + indic_8x(2*i+1,2*j,2*k+1) * (.5*field_8x(2*i+1,2*j,2*k+1) + (bary_compo==0           || dir==1)*.5))/(4*indic);
                    }
                  else
                    {
                      field_ft(i,j,k) = .5;
                    }
                }
            }
        }
    }
  else if (face_dir == 2)
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                  // qu'une cellule consideree pure soit traversee par l'interface.
                  // Pour coherence, on considere les faces non coupees dans ce cas.
                  //
                  // Une cellule adjacente pure suggere ou bien qu'il n'y a pas d'intersections
                  // ou bien qu'il y a des intersections mais que l'indicatrice resultante est
                  // trop faible et en dessous du seuil dans calculer_indicatrice.
                  // Pour coherence, on neglige la coupure de la surface egalement.
                  bool cellule_cote_positif_non_pure = (!((indicatrice_ft(i, j, k) == 0.) || (indicatrice_ft(i, j, k) == 1.)));
                  bool cellule_cote_negatif_non_pure = (!((indicatrice_ft(i, j, k-1) == 0.) || (indicatrice_ft(i, j, k-1) == 1.)));
                  if (cellule_cote_positif_non_pure && cellule_cote_negatif_non_pure)
                    {
                      double indic = .25*(indic_8x(2*i,2*j,2*k)       + indic_8x(2*i+1,2*j,2*k)
                                          + indic_8x(2*i,2*j+1,2*k)   + indic_8x(2*i+1,2*j+1,2*k));
                      field_ft(i,j,k) = (indic == 0 || indic == 1) ? .5 :
                                        (  indic_8x(2*i,2*j,2*k)   * (.5*field_8x(2*i,2*j,2*k)                                    )   + indic_8x(2*i+1,2*j,2*k)   * (.5*field_8x(2*i+1,2*j,2*k)    + (dir==1                    )*.5)
                                           + indic_8x(2*i,2*j+1,2*k) * (.5*field_8x(2*i,2*j+1,2*k) + (          dir==0          )*.5)   + indic_8x(2*i+1,2*j+1,2*k) * (.5*field_8x(2*i+1,2*j+1,2*k) + (dir==1 || dir==0          )*.5))/(4*indic);
                    }
                  else
                    {
                      field_ft(i,j,k) = .5;
                    }
                }
            }
        }
    }
  else
    {
      Cerr << "moyenne_barycentre_phase1_sous_faces_face: erroneous value of face_dir" << finl;
      Process::exit();
    }

  field_ft.echange_espace_virtuel(ghost);
}

