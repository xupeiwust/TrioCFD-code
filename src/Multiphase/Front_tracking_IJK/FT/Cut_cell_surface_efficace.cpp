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
// File      : Cut_cell_surface_efficace.cpp
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#include <Cut_cell_surface_efficace.h>
#include <IJK_Field_vector.h>
#include <Cut_cell_FT_Disc.h>
#include <IJK_Thermal_base.h>
#include <Process.h>
#include <stat_counters.h>
#include <IJK_Navier_Stokes_tools.h>

#include <SolveurSys.h>
#include <Matrice_Morse.h>
#include <Matrice_Bloc.h>
#include <Matrice_Dense.h>
#include <EChaine.h>

extern "C" {
  void dgelsd_(int* M, int* N, int* NRHS, double* A, int* lda, double* B, int* ldb, double* S, double* RCOND, int* rank, double* WORK, int* lwork, int* iwork, int* INFO);
}

// Solve A.x = b using dgelsd_ from Lapack.
// A has dimension (M,N) but is in column-major order
// b has dimension (M)
// x has dimension (N)
// The result is overwritten within b.
void dgelsd_resoudre_systeme(double* A, double* B, int M, int N)
{
  int NRHS = 1;
  int INFO;
  int LDB = std::max(N,M);
  int rank;
  double RCOND = -1.;

  int NLVL = std::max(0, int(std::log2(std::min(N,M)/(2+1))) + 1);
  int LIWORK = std::max(1, 3*std::min(N,M)*NLVL + 11*std::min(N,M));
  int *IWORK = new int[LIWORK];

  int S_size = std::min(N,M);
  double *S = new double[S_size];

  int LWORK = -1;
  double LWORK_opt;

  // Query LWORK size
  dgelsd_(&M,&N,&NRHS,A,&M,B,&LDB,S,&RCOND,&rank,&LWORK_opt,&LWORK,IWORK,&INFO);

  LWORK = (int)LWORK_opt;
  double *WORK = new double[LWORK];

  // Solve
  dgelsd_(&M,&N,&NRHS,A,&M,B,&LDB,S,&RCOND,&rank,WORK,&LWORK,IWORK,&INFO);

  delete[] WORK;
}


void Cut_cell_surface_efficace::calcul_surface_interface_efficace_initiale(
  int methode_explicite,
  const IJK_Field_double& old_indicatrice_ns,
  const IJK_Field_double& next_indicatrice_ns,
  const IJK_Field_double& surface_interface_ns_old,
  const IJK_Field_double& surface_interface_ns_next,
  const IJK_Field_vector3_double& normal_of_interf_ns_old,
  const IJK_Field_vector3_double& normal_of_interf_ns_next,
  DoubleTabFT_cut_cell_vector3& normale_deplacement_interface,
  DoubleTabFT_cut_cell_scalar& surface_efficace_interface,
  DoubleTabFT_cut_cell_scalar& surface_efficace_interface_initial)
{
  const Cut_cell_FT_Disc& cut_cell_disc = surface_efficace_interface.get_cut_cell_disc();
  for (int n = 0; n < cut_cell_disc.get_n_tot(); n++)
    {
      Int3 ijk = cut_cell_disc.get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      if (IJK_Interfaces::est_pure(.5*(old_indicatrice_ns(i, j, k) + next_indicatrice_ns(i, j, k))))
        {
          // Si la cellule est purement monophasique, il n'y a pas d'intersection avec l'interface
          normale_deplacement_interface(n,0) = 0.;
          normale_deplacement_interface(n,1) = 0.;
          normale_deplacement_interface(n,2) = 0.;
          surface_efficace_interface(n) = 0.;
          surface_efficace_interface_initial(n) = 0.;
        }
      else
        {
          if (methode_explicite) // methode explicite : on se place a de S_n
            {
              if (IJK_Interfaces::devient_pure(old_indicatrice_ns(i, j, k), next_indicatrice_ns(i, j, k)))
                {
                  surface_efficace_interface(n) = surface_interface_ns_old(i, j, k);
                  for (int dir = 0 ; dir < 3 ; dir++)
                    {
                      normale_deplacement_interface(n,dir) = normal_of_interf_ns_old[dir](i, j, k);
                    }
                }
              else if (IJK_Interfaces::devient_diphasique(old_indicatrice_ns(i, j, k), next_indicatrice_ns(i, j, k)))
                {
                  surface_efficace_interface(n) = surface_interface_ns_next(i, j, k);
                  for (int dir = 0 ; dir < 3 ; dir++)
                    {
                      normale_deplacement_interface(n,dir) = normal_of_interf_ns_next[dir](i, j, k);
                    }
                }
              else
                {
                  surface_efficace_interface(n) = surface_interface_ns_old(i, j, k);
                  for (int dir = 0 ; dir < 3 ; dir++)
                    {
                      normale_deplacement_interface(n,dir) = normal_of_interf_ns_old[dir](i, j, k);
                    }
                }
            }
          else // methode algebrique : on ne place au milieu de S_n et S_n+1
            {
              if (IJK_Interfaces::devient_pure(old_indicatrice_ns(i, j, k), next_indicatrice_ns(i, j, k)))
                {
                  surface_efficace_interface(n) = (surface_interface_ns_old(i, j, k) + 3*surface_interface_ns_next(i, j, k)) * 0.25;
                  for (int dir = 0 ; dir < 3 ; dir++)
                    {
                      normale_deplacement_interface(n,dir) = normal_of_interf_ns_old[dir](i, j, k);
                    }
                }
              else if (IJK_Interfaces::devient_diphasique(old_indicatrice_ns(i, j, k), next_indicatrice_ns(i, j, k)))
                {
                  surface_efficace_interface(n) = (3*surface_interface_ns_old(i, j, k) + surface_interface_ns_next(i, j, k)) * 0.25;
                  for (int dir = 0 ; dir < 3 ; dir++)
                    {
                      normale_deplacement_interface(n,dir) = normal_of_interf_ns_next[dir](i, j, k);
                    }
                }
              else
                {
                  surface_efficace_interface(n) = (surface_interface_ns_old(i, j, k) + surface_interface_ns_next(i, j, k)) * 0.5;
                  for (int dir = 0 ; dir < 3 ; dir++)
                    {
                      normale_deplacement_interface(n,dir) = (normal_of_interf_ns_old[dir](i, j, k) + normal_of_interf_ns_next[dir](i, j, k)) * 0.5;
                    }
                }
            }

          double norm_normale = sqrt(normale_deplacement_interface(n,0)*normale_deplacement_interface(n,0) + normale_deplacement_interface(n,1)*normale_deplacement_interface(n,1) + normale_deplacement_interface(n,2)*normale_deplacement_interface(n,2));
          normale_deplacement_interface(n,0) /= norm_normale;
          normale_deplacement_interface(n,1) /= norm_normale;
          normale_deplacement_interface(n,2) /= norm_normale;

          surface_efficace_interface_initial(n) = surface_efficace_interface(n);
        }
    }
}

void Cut_cell_surface_efficace::calcul_vitesse_interface(
  const IJK_Field_vector3_double& velocity,
  const IJK_Field_double& old_indicatrice_ns,
  const IJK_Field_double& next_indicatrice_ns,
  const IJK_Field_vector3_double& barycentre_phase1_ns_old,
  const IJK_Field_vector3_double& barycentre_phase1_ns_next,
  DoubleTabFT_cut_cell_vector3& coord_deplacement_interface,
  DoubleTabFT_cut_cell_vector3& vitesse_deplacement_interface)
{
  const Cut_cell_FT_Disc& cut_cell_disc = vitesse_deplacement_interface.get_cut_cell_disc();
  const IJK_Grid_Geometry& geom = cut_cell_disc.get_splitting().get_grid_geometry();

  // Calculer la coordonnee correspondant au deplacement de l'interface :
  //    x_depl_interf = (I[new]*bary[new] - I[old]*bary[old])/(delta I)
  // Bien sur, celle-ci ne depend pas de la phase consideree.
  for (int n = 0; n < cut_cell_disc.get_n_loc(); n++)
    {
      Int3 ijk = cut_cell_disc.get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      for (int dir = 0; dir < 3; dir++)
        {
          double delta_I = next_indicatrice_ns(i,j,k) - old_indicatrice_ns(i,j,k);

          double bary_depl_interf;
          if (delta_I == 0.)
            {
              bary_depl_interf = barycentre_phase1_ns_next[dir](i,j,k);
            }
          else
            {
              bary_depl_interf = (next_indicatrice_ns(i,j,k)*barycentre_phase1_ns_next[dir](i,j,k) - old_indicatrice_ns(i,j,k)*barycentre_phase1_ns_old[dir](i,j,k))/delta_I;

              // Le calcul est pas correct (!), mais on evite de planter
              if (bary_depl_interf <= 0)
                {
                  bary_depl_interf = 1e-12;
                }
              else if (bary_depl_interf > 1)
                {
                  bary_depl_interf = 1 - 1e-12;
                }
            }
          assert((bary_depl_interf >= 0) && (bary_depl_interf <= 1));

          // Passage a un systeme de coordonnees dimensionel et absolue
          const int i_dir = select_dir(dir, i, j, k);
          const int offset_dir = cut_cell_disc.get_splitting().get_offset_local(dir);
          const double origin_dir = geom.get_origin(dir);
          const double delta_dir = geom.get_constant_delta(dir);
          const double coord_dir  = origin_dir + (i_dir + offset_dir + bary_depl_interf) * delta_dir;

          coord_deplacement_interface(n, dir) = coord_dir;
        }
    }

  // Interpolation de la vitesse sur la coordonnee du deplacement de l'interface
  ArrOfDouble vinterp_component(cut_cell_disc.get_n_loc());
  for (int dir = 0; dir < 3; dir++)
    {
      ijk_interpolate_skip_unknown_points(velocity[dir], coord_deplacement_interface, vinterp_component, 1.e5 /* value for unknown points */);
      for (int i = 0; i < cut_cell_disc.get_n_loc(); i++)
        {
          vitesse_deplacement_interface(i, dir) = vinterp_component[i];
        }
    }
  vitesse_deplacement_interface.echange_espace_virtuel();
}

void Cut_cell_surface_efficace::calcul_surface_interface_efficace(
  double timestep,
  const IJK_Field_vector3_double& velocity,
  const IJK_Field_double& old_indicatrice_ns,
  const IJK_Field_double& next_indicatrice_ns,
  const DoubleTabFT_cut_cell_vector3& vitesse_deplacement_interface,
  const DoubleTabFT_cut_cell_vector3& normale_deplacement_interface,
  DoubleTabFT_cut_cell_scalar& surface_efficace_interface)
{
  static Stat_Counter_Id calculer_surface_efficace_interface_counter_ =
    statistiques().new_counter(2, "cut_cell: calcul des surfaces efficaces interface");
  statistiques().begin_count(calculer_surface_efficace_interface_counter_);

  const Cut_cell_FT_Disc& cut_cell_disc = surface_efficace_interface.get_cut_cell_disc();
  const IJK_Grid_Geometry& geom = cut_cell_disc.get_splitting().get_grid_geometry();

  const double delta_x = geom.get_constant_delta(0);
  const double delta_y = geom.get_constant_delta(1);
  const int offset = cut_cell_disc.get_splitting().get_offset_local(DIRECTION_K);
  const ArrOfDouble& delta_z_all = geom.get_delta(DIRECTION_K);

  // Calcul de la correction de surface efficace de l'interface
  for (int n = 0; n < cut_cell_disc.get_n_loc(); n++)
    {
      Int3 ijk = cut_cell_disc.get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      const double delta_z = delta_z_all[k + offset];
      const double volume = delta_x * delta_y * delta_z;

      const double vitesse_interface = -vitesse_deplacement_interface(n,0)*normale_deplacement_interface(n,0) - vitesse_deplacement_interface(n,1)*normale_deplacement_interface(n,1) - vitesse_deplacement_interface(n,2)*normale_deplacement_interface(n,2);
      const double surface_interface = surface_efficace_interface(n);

      double delta_volume_cible = volume * (next_indicatrice_ns(i, j, k) - old_indicatrice_ns(i, j, k));
      double delta_volume_interface = vitesse_interface * surface_interface * timestep;

      double erreur_relative = delta_volume_cible == 0. ? delta_volume_interface/volume : std::fabs((delta_volume_interface - delta_volume_cible)/delta_volume_cible);

      if (erreur_relative == 0)
        {
        }
      else
        {
          surface_efficace_interface(n) = delta_volume_cible/(vitesse_interface * timestep);
        }
    }

  surface_efficace_interface.echange_espace_virtuel();

  statistiques().end_count(calculer_surface_efficace_interface_counter_);
}

// calcul_surface_face_efficace_initiale: version with DoubleTabFT_cut_cell output (use cut-cell structures)
void Cut_cell_surface_efficace::calcul_surface_face_efficace_initiale(
  int methode_explicite,
  const IJK_Field_vector3_double& old_indicatrice_surfacique_face_ns,
  const IJK_Field_vector3_double& next_indicatrice_surfacique_face_ns,
  DoubleTabFT_cut_cell_vector3& indicatrice_surfacique_efficace_face,
  DoubleTabFT_cut_cell_vector3& indicatrice_surfacique_efficace_face_initial)
{
  if (methode_explicite) // methode explicite : on se place a de S_n
    {
      const Cut_cell_FT_Disc& cut_cell_disc = indicatrice_surfacique_efficace_face.get_cut_cell_disc();
      for (int n = 0; n < cut_cell_disc.get_n_tot(); n++)
        {
          Int3 ijk = cut_cell_disc.get_ijk(n);
          int i = ijk[0];
          int j = ijk[1];
          int k = ijk[2];

          for (int dir = 0; dir < 3; dir++)
            {
              if ((old_indicatrice_surfacique_face_ns[dir](i, j, k) == 1) && (next_indicatrice_surfacique_face_ns[dir](i, j, k) == 0))
                {
                  indicatrice_surfacique_efficace_face(n, dir) = .5;
                }
              else if ((old_indicatrice_surfacique_face_ns[dir](i, j, k) == 0) && (next_indicatrice_surfacique_face_ns[dir](i, j, k) == 1))
                {
                  indicatrice_surfacique_efficace_face(n, dir) = .5;
                }
              else if (IJK_Interfaces::devient_pure(old_indicatrice_surfacique_face_ns[dir](i, j, k), next_indicatrice_surfacique_face_ns[dir](i, j, k)))
                {
                  indicatrice_surfacique_efficace_face(n, dir) = old_indicatrice_surfacique_face_ns[dir](i, j, k);
                }
              else if (IJK_Interfaces::devient_diphasique(old_indicatrice_surfacique_face_ns[dir](i, j, k), next_indicatrice_surfacique_face_ns[dir](i, j, k)))
                {
                  indicatrice_surfacique_efficace_face(n, dir) = next_indicatrice_surfacique_face_ns[dir](i, j, k);
                }
              else
                {
                  indicatrice_surfacique_efficace_face(n, dir) = old_indicatrice_surfacique_face_ns[dir](i, j, k);
                }
              indicatrice_surfacique_efficace_face_initial(n, dir) = indicatrice_surfacique_efficace_face(n, dir);
            }
        }
    }
  else // methode algebrique : on ne place au milieu de S_n et S_n+1
    {
      const Cut_cell_FT_Disc& cut_cell_disc = indicatrice_surfacique_efficace_face.get_cut_cell_disc();
      for (int n = 0; n < cut_cell_disc.get_n_tot(); n++)
        {
          Int3 ijk = cut_cell_disc.get_ijk(n);
          int i = ijk[0];
          int j = ijk[1];
          int k = ijk[2];

          for (int dir = 0; dir < 3; dir++)
            {
              if (IJK_Interfaces::devient_pure(old_indicatrice_surfacique_face_ns[dir](i, j, k), next_indicatrice_surfacique_face_ns[dir](i, j, k)))
                {
                  indicatrice_surfacique_efficace_face(n, dir) = (old_indicatrice_surfacique_face_ns[dir](i, j, k) + 3*next_indicatrice_surfacique_face_ns[dir](i, j, k)) * 0.25;
                }
              else if (IJK_Interfaces::devient_diphasique(old_indicatrice_surfacique_face_ns[dir](i, j, k), next_indicatrice_surfacique_face_ns[dir](i, j, k)))
                {
                  indicatrice_surfacique_efficace_face(n, dir) = (3*old_indicatrice_surfacique_face_ns[dir](i, j, k) + next_indicatrice_surfacique_face_ns[dir](i, j, k)) * 0.25;
                }
              else
                {
                  indicatrice_surfacique_efficace_face(n, dir) = (old_indicatrice_surfacique_face_ns[dir](i, j, k) + next_indicatrice_surfacique_face_ns[dir](i, j, k)) * 0.5;
                }
              indicatrice_surfacique_efficace_face_initial(n, dir) = indicatrice_surfacique_efficace_face(n, dir);
            }
        }
    }
}

// calcul_surface_face_efficace_initiale: version with IJK_Field output (does not use cut-cell structures)
void Cut_cell_surface_efficace::calcul_surface_face_efficace_initiale(
  int methode_explicite,
  const IJK_Field_vector3_double& old_indicatrice_surfacique_face_ns,
  const IJK_Field_vector3_double& next_indicatrice_surfacique_face_ns,
  IJK_Field_vector3_double& indicatrice_surfacique_efficace_face,
  IJK_Field_vector3_double& indicatrice_surfacique_efficace_face_initial)
{
  if (methode_explicite) // methode explicite : on se place a de S_n
    {
      for (int dir = 0; dir < 3; dir++)
        {
          const int ni = indicatrice_surfacique_efficace_face[dir].ni();
          const int nj = indicatrice_surfacique_efficace_face[dir].nj();
          const int nk = indicatrice_surfacique_efficace_face[dir].nk();
          const int ghost = indicatrice_surfacique_efficace_face[dir].ghost();

          for (int k = -ghost; k < nk+ghost; k++)
            {
              for (int j = -ghost; j < nj+ghost; j++)
                {
                  for (int i = -ghost; i < ni+ghost; i++)
                    {
                      if ((old_indicatrice_surfacique_face_ns[dir](i, j, k) == 1) && (next_indicatrice_surfacique_face_ns[dir](i, j, k) == 0))
                        {
                          indicatrice_surfacique_efficace_face[dir](i, j, k) = .5;
                        }
                      else if ((old_indicatrice_surfacique_face_ns[dir](i, j, k) == 0) && (next_indicatrice_surfacique_face_ns[dir](i, j, k) == 1))
                        {
                          indicatrice_surfacique_efficace_face[dir](i, j, k) = .5;
                        }
                      else if (IJK_Interfaces::devient_pure(old_indicatrice_surfacique_face_ns[dir](i, j, k), next_indicatrice_surfacique_face_ns[dir](i, j, k)))
                        {
                          indicatrice_surfacique_efficace_face[dir](i, j, k) = old_indicatrice_surfacique_face_ns[dir](i, j, k);
                        }
                      else if (IJK_Interfaces::devient_diphasique(old_indicatrice_surfacique_face_ns[dir](i, j, k), next_indicatrice_surfacique_face_ns[dir](i, j, k)))
                        {
                          indicatrice_surfacique_efficace_face[dir](i, j, k) = next_indicatrice_surfacique_face_ns[dir](i, j, k);
                        }
                      else
                        {
                          indicatrice_surfacique_efficace_face[dir](i, j, k) = old_indicatrice_surfacique_face_ns[dir](i, j, k);
                        }
                      indicatrice_surfacique_efficace_face_initial[dir](i, j, k) = indicatrice_surfacique_efficace_face[dir](i, j, k);
                    }
                }
            }
        }
    }
  else // methode algebrique : on ne place au milieu de S_n et S_n+1
    {
      for (int dir = 0; dir < 3; dir++)
        {
          const int ni = indicatrice_surfacique_efficace_face[dir].ni();
          const int nj = indicatrice_surfacique_efficace_face[dir].nj();
          const int nk = indicatrice_surfacique_efficace_face[dir].nk();
          const int ghost = indicatrice_surfacique_efficace_face[dir].ghost();

          for (int k = -ghost; k < nk+ghost; k++)
            {
              for (int j = -ghost; j < nj+ghost; j++)
                {
                  for (int i = -ghost; i < ni+ghost; i++)
                    {
                      if (IJK_Interfaces::devient_pure(old_indicatrice_surfacique_face_ns[dir](i, j, k), next_indicatrice_surfacique_face_ns[dir](i, j, k)))
                        {
                          indicatrice_surfacique_efficace_face[dir](i, j, k) = (old_indicatrice_surfacique_face_ns[dir](i, j, k) + 3*next_indicatrice_surfacique_face_ns[dir](i, j, k)) * 0.25;
                        }
                      else if (IJK_Interfaces::devient_diphasique(old_indicatrice_surfacique_face_ns[dir](i, j, k), next_indicatrice_surfacique_face_ns[dir](i, j, k)))
                        {
                          indicatrice_surfacique_efficace_face[dir](i, j, k) = (3*old_indicatrice_surfacique_face_ns[dir](i, j, k) + next_indicatrice_surfacique_face_ns[dir](i, j, k)) * 0.25;
                        }
                      else
                        {
                          indicatrice_surfacique_efficace_face[dir](i, j, k) = (old_indicatrice_surfacique_face_ns[dir](i, j, k) + next_indicatrice_surfacique_face_ns[dir](i, j, k)) * 0.5;
                        }
                      indicatrice_surfacique_efficace_face_initial[dir](i, j, k) = indicatrice_surfacique_efficace_face[dir](i, j, k);
                    }
                }
            }
        }
    }
}


void Cut_cell_surface_efficace::calcul_surface_face_efficace_iteratif(
  int verbosite_surface_efficace_face,
  double timestep,
  const Cut_field_vector3_double& velocity,
  int& iteration_solver_surface_efficace_face,
  const IJK_Field_double& old_indicatrice_ns,
  const IJK_Field_double& next_indicatrice_ns,
  const IJK_Field_vector3_double& old_indicatrice_surfacique_face_ns,
  const IJK_Field_vector3_double& next_indicatrice_surfacique_face_ns,
  DoubleTabFT_cut_cell_vector3& indicatrice_surfacique_efficace_face,
  const DoubleTabFT_cut_cell_vector3& indicatrice_surfacique_efficace_face_initial,
  DoubleTabFT_cut_cell_vector6& indicatrice_surfacique_efficace_face_correction,
  DoubleTabFT_cut_cell_scalar& indicatrice_surfacique_efficace_face_absolute_error)
{
  static Stat_Counter_Id calculer_surface_efficace_face_counter_ =
    statistiques().new_counter(2, "cut_cell: calcul des surfaces efficaces");
  statistiques().begin_count(calculer_surface_efficace_face_counter_);

  const Cut_cell_FT_Disc& cut_cell_disc = indicatrice_surfacique_efficace_face.get_cut_cell_disc();
  const IJK_Grid_Geometry& geom = cut_cell_disc.get_splitting().get_grid_geometry();
  const double delta_x = geom.get_constant_delta(0);
  const double delta_y = geom.get_constant_delta(1);
  const int offset = cut_cell_disc.get_splitting().get_offset_local(DIRECTION_K);
  const ArrOfDouble& delta_z_all = geom.get_delta(DIRECTION_K);

  for (int n = 0; n < cut_cell_disc.get_n_loc(); n++)
    {
      indicatrice_surfacique_efficace_face_absolute_error(n) = 0.;
    }

  // Iteration du solveur de la surface efficace
  const int maximum_iteration = 499;
  const int iteration_impression_intermediaire = 600;
  int solution_not_found = 0;
  do
    {
      iteration_solver_surface_efficace_face += 1;
      solution_not_found = 0;

      int solution_locally_not_found = 0;

      // Initialisation du tableau de correction de surface
      for (int n = 0; n < cut_cell_disc.get_n_loc(); n++)
        {
          Int3 ijk = cut_cell_disc.get_ijk(n);
          int i = ijk[0];
          int j = ijk[1];
          int k = ijk[2];

          for (int num_face = 0; num_face < 6; num_face++)
            {
              int dir = num_face%3;
              int decalage = num_face/3;

              int di = decalage*(dir == 0);
              int dj = decalage*(dir == 1);
              int dk = decalage*(dir == 2);

              int n_face = cut_cell_disc.get_n_face(num_face, n, i, j, k);
              if (n_face >= 0)
                {
                  indicatrice_surfacique_efficace_face_correction(n, num_face) = indicatrice_surfacique_efficace_face(n_face, dir);
                }
              else
                {
                  assert((next_indicatrice_ns(i+di,j+dj,k+dk) == 0) || (next_indicatrice_ns(i+di,j+dj,k+dk) == 1));
                  indicatrice_surfacique_efficace_face_correction(n, num_face) = next_indicatrice_ns(i+di,j+dj,k+dk);
                }
            }
        }

      // Calcul de la correction de surface
      // Note : Contrairement aux autres tableaux de surface, la correction est stockee sur
      // les six faces de chaque cellule diphasique. On a donc pour chaque face deux informations
      // de correction, qui correspondent a la correction calculee independamment pour les deux
      // cellules contenant la face. On calcule dans un deuxieme temps la surface efficace
      // a partir de ces deux corrections.
      for (int n = 0; n < cut_cell_disc.get_n_loc(); n++)
        {
          Int3 ijk = cut_cell_disc.get_ijk(n);
          int i = ijk[0];
          int j = ijk[1];
          int k = ijk[2];

          const double delta_z = delta_z_all[k + offset];
          const double volume = delta_x * delta_y * delta_z;

          const double fx = delta_y * delta_z * timestep;
          const double fy = delta_x * delta_z * timestep;
          const double fz = delta_x * delta_y * timestep;

          double surface_efficace[6] = {0};
          double surface_max[6] = {0};
          double surface_min[6] = {0};
          double flux_diphasique[6] = {0};
          double delta_surface_max[6] = {0};

          double delta_volume_total = 0;
          double delta_volume_diphasique = 0;
          for (int num_face = 0; num_face < 6; num_face++)
            {
              int dir = num_face%3;
              int decalage = num_face/3;
              int sign = decalage*2 -1;

              int di = decalage*(dir == 0);
              int dj = decalage*(dir == 1);
              int dk = decalage*(dir == 2);

              double f = select_dir(dir, fx, fy, fz);

              int n_face = cut_cell_disc.get_n_face(num_face, n, i, j, k);
              if (n_face >= 0)
                {
                  double surface_old = old_indicatrice_surfacique_face_ns[dir](i+di,j+dj,k+dk);
                  double surface_next = next_indicatrice_surfacique_face_ns[dir](i+di,j+dj,k+dk);
                  surface_max[num_face] = std::max(surface_old, surface_next);
                  surface_min[num_face] = std::min(surface_old, surface_next);
                  surface_efficace[num_face] = indicatrice_surfacique_efficace_face(n_face, dir);
                  flux_diphasique[num_face] = (indicatrice_surfacique_efficace_face(n_face, dir) == 1.) ? 0. : sign*f*indicatrice_surfacique_efficace_face(n_face, dir)*velocity[dir].diph_l_(n_face);
                  if (surface_efficace[num_face] == 0.)
                    {
                      assert(surface_max[num_face] == 0.);
                    }

                  delta_volume_total -= sign*f*surface_efficace[num_face]*velocity[dir].diph_l_(n_face);
                  delta_volume_diphasique -= flux_diphasique[num_face];

                }
              else
                {
                  delta_volume_total -= sign*f*next_indicatrice_ns(i+di,j+dj,k+dk)*velocity[dir].pure_(i+di,j+dj,k+dk);
                }

            }

          double delta_volume_cible = volume * (next_indicatrice_ns(i, j, k) - old_indicatrice_ns(i, j, k));
          double delta_volume_pure = delta_volume_total - delta_volume_diphasique;
          double delta_volume_diphasique_cible = delta_volume_cible - delta_volume_pure;

          double erreur_absolue = std::fabs(delta_volume_total - delta_volume_cible);
          double erreur_relative = delta_volume_cible == 0. ? delta_volume_total/volume : std::fabs((delta_volume_total - delta_volume_cible)/delta_volume_cible);
          double ecart_volume_diphasique = delta_volume_diphasique - delta_volume_diphasique_cible;

          indicatrice_surfacique_efficace_face_absolute_error(n) = erreur_absolue;

          if (ecart_volume_diphasique == 0)
            {
            }
          else
            {
              if (erreur_relative > 1e-5)
                {
                  solution_locally_not_found = 1;
                }

              double somme_ratio_max_sur_requis = 0;
              for (int num_face = 0; num_face < 6; num_face++)
                {
                  if (flux_diphasique[num_face] != 0)
                    {
                      bool flux_et_ecart_meme_signe = ((flux_diphasique[num_face] < 0) == (ecart_volume_diphasique < 0));

                      delta_surface_max[num_face] = flux_et_ecart_meme_signe ? surface_max[num_face] - surface_efficace[num_face] : surface_efficace[num_face] - surface_min[num_face];

                      // Note : Le calcul de delta_surface_requis n'est pas fait si la surface est tres faible
                      // En effet, cela suggere que le flux_diphasique[num_face] est tres faible egalement, ce qui peut provoquer une 'Arithmetic exception'.
                      // De plus, modifier davantage la surface efficace dans ce cas n'est pas utile.
                      double delta_surface_requis = (std::abs(surface_efficace[num_face]) < 1e-80) ? 0. : std::fabs(ecart_volume_diphasique/flux_diphasique[num_face])*surface_efficace[num_face];
                      somme_ratio_max_sur_requis += (delta_surface_requis == 0.) ? 0. : delta_surface_max[num_face]/delta_surface_requis;
                    }
                }

              for (int num_face = 0; num_face < 6; num_face++)
                {
                  if (flux_diphasique[num_face] != 0)
                    {
                      int dir = num_face%3;

                      int n_face = cut_cell_disc.get_n_face(num_face, n, i, j, k);
                      assert(n_face >= 0);

                      bool flux_et_ecart_meme_signe = ((flux_diphasique[num_face] < 0) == (ecart_volume_diphasique < 0));
                      int sign = 2*flux_et_ecart_meme_signe - 1;

                      // Note : La correction est choisi telle que chaque surface contribue la meme fraction de delta_S_max^i.
                      // Normalement, cette fraction est toujours inferieur a 100%.
                      // La correction est :
                      //     delta_S^i = delta_S_max^i / (somme delta_S_max/delta_S_requis),
                      // ce qui implique (delta_S^i/delta_S_max^i) est constant, et par ailleurs somme delta_S^i/delta_S_requis^i = 1
                      double correction = (somme_ratio_max_sur_requis == 0) ? 0. : sign*delta_surface_max[num_face]/somme_ratio_max_sur_requis;

                      indicatrice_surfacique_efficace_face_correction(n, num_face) = indicatrice_surfacique_efficace_face(n_face, dir) + correction;

                      {
                        double minimal_acceptable_surface = 1e-15;
                        if (indicatrice_surfacique_efficace_face_correction(n, num_face) < (surface_min[num_face] + minimal_acceptable_surface))
                          {
                            indicatrice_surfacique_efficace_face_correction(n, num_face) = std::min(surface_min[num_face] + minimal_acceptable_surface, .5*(surface_min[num_face] + surface_max[num_face]));
                          }
                        if (indicatrice_surfacique_efficace_face_correction(n, num_face) > (surface_max[num_face] - minimal_acceptable_surface))
                          {
                            indicatrice_surfacique_efficace_face_correction(n, num_face) = std::max(surface_max[num_face] - minimal_acceptable_surface, .5*(surface_min[num_face] + surface_max[num_face]));
                          }
                      }
                    }
                }

            }
        }

      indicatrice_surfacique_efficace_face_correction.echange_espace_virtuel();

      // Mise a jour de la surface efficace
      for (int n = 0; n < cut_cell_disc.get_n_loc(); n++)
        {
          Int3 ijk = cut_cell_disc.get_ijk(n);
          int i = ijk[0];
          int j = ijk[1];
          int k = ijk[2];

          for (int dir = 0; dir < 3; dir++)
            {
              int di = (-1)*(dir == 0);
              int dj = (-1)*(dir == 1);
              int dk = (-1)*(dir == 2);

              int n_face = cut_cell_disc.get_n(i+di, j+dj, k+dk);

              if (n_face >= 0)
                {
                  double new_surface = (indicatrice_surfacique_efficace_face_absolute_error(n_face) > indicatrice_surfacique_efficace_face_absolute_error(n)) ? indicatrice_surfacique_efficace_face_correction(n_face, dir+3) : indicatrice_surfacique_efficace_face_correction(n, dir);

                  // Si une des cellules disparait, on liu donne la priorite dans le calcul de la correction
                  bool centre_switch = (IJK_Interfaces::devient_pure(old_indicatrice_ns(i,j,k), next_indicatrice_ns(i,j,k)) || (IJK_Interfaces::devient_diphasique(old_indicatrice_ns(i,j,k), next_indicatrice_ns(i,j,k))));
                  bool decale_switch = (IJK_Interfaces::devient_pure(old_indicatrice_ns(i+di,j+dj,k+dk), next_indicatrice_ns(i+di,j+dj,k+dk)) || (IJK_Interfaces::devient_diphasique(old_indicatrice_ns(i+di,j+dj,k+dk), next_indicatrice_ns(i+di,j+dj,k+dk))));
                  if ((centre_switch) && (!decale_switch))
                    {
                      new_surface = indicatrice_surfacique_efficace_face_correction(n, dir);
                    }
                  else if ((!centre_switch) && (decale_switch))
                    {
                      new_surface = indicatrice_surfacique_efficace_face_correction(n_face, dir+3);
                    }

                  // Filtrage temporel de la correction
                  //indicatrice_surfacique_efficace_face(n, dir) = 0.25*(3*indicatrice_surfacique_efficace_face(n, dir) + new_surface);
                  indicatrice_surfacique_efficace_face(n, dir) = new_surface;
                }
            }
        }

      indicatrice_surfacique_efficace_face.echange_espace_virtuel();

      if (iteration_solver_surface_efficace_face%iteration_impression_intermediaire == 0)
        {
          imprimer_informations_surface_efficace_face(
            verbosite_surface_efficace_face,
            iteration_solver_surface_efficace_face,
            timestep,
            velocity,
            old_indicatrice_ns,
            next_indicatrice_ns,
            indicatrice_surfacique_efficace_face,
            indicatrice_surfacique_efficace_face_initial);
        }

      solution_not_found = Process::mp_max(solution_locally_not_found);
    }
  while (solution_not_found && iteration_solver_surface_efficace_face < maximum_iteration);

  statistiques().end_count(calculer_surface_efficace_face_counter_);
}

void Cut_cell_surface_efficace::imprimer_informations_surface_efficace_interface(
  int verbosite_surface_efficace_interface,
  double timestep,
  const IJK_Field_vector3_double& velocity,
  const IJK_Field_double& old_indicatrice_ns,
  const IJK_Field_double& next_indicatrice_ns,
  const DoubleTabFT_cut_cell_scalar& surface_efficace_interface,
  const DoubleTabFT_cut_cell_scalar& surface_efficace_interface_initial,
  const DoubleTabFT_cut_cell_vector3& normale_deplacement_interface,
  const DoubleTabFT_cut_cell_vector3& vitesse_deplacement_interface)
{
  const Cut_cell_FT_Disc& cut_cell_disc = surface_efficace_interface.get_cut_cell_disc();
  const IJK_Grid_Geometry& geom = cut_cell_disc.get_splitting().get_grid_geometry();
  const double delta_x = geom.get_constant_delta(0);
  const double delta_y = geom.get_constant_delta(1);
  const int offset = cut_cell_disc.get_splitting().get_offset_local(DIRECTION_K);
  const ArrOfDouble& delta_z_all = geom.get_delta(DIRECTION_K);

  // Impression des resultats
  int n_count = 0;
  int n_count_surface = 0;
  double moyenne_erreur_absolue = 0.;
  double moyenne_erreur_relative = 0.;
  double maximum_erreur_absolue = 0.;
  double maximum_erreur_relative = 0.;
  double moyenne_ecart_absolue_surface = 0.;
  double moyenne_ecart_relatif_surface = 0.;
  double maximum_ecart_absolue_surface = 0.;
  double maximum_ecart_relatif_surface = 0.;
  int switch_n_count = 0;
  int switch_n_count_surface = 0;
  double switch_moyenne_erreur_absolue = 0.;
  double switch_moyenne_erreur_relative = 0.;
  double switch_maximum_erreur_absolue = 0.;
  double switch_maximum_erreur_relative = 0.;
  double switch_moyenne_ecart_absolue_surface = 0.;
  double switch_moyenne_ecart_relatif_surface = 0.;
  double switch_maximum_ecart_absolue_surface = 0.;
  double switch_maximum_ecart_relatif_surface = 0.;
  for (int n = 0; n < cut_cell_disc.get_n_loc(); n++)
    {
      Int3 ijk = cut_cell_disc.get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      bool centre_switch = (IJK_Interfaces::devient_pure(old_indicatrice_ns(i,j,k), next_indicatrice_ns(i,j,k)) || (IJK_Interfaces::devient_diphasique(old_indicatrice_ns(i,j,k), next_indicatrice_ns(i,j,k))));

      const double delta_z = delta_z_all[k + offset];
      const double volume = delta_x * delta_y * delta_z;

      {
        double erreur_relative;
        if (surface_efficace_interface_initial(n) != 0.)
          {
            erreur_relative = std::fabs((surface_efficace_interface(n) - surface_efficace_interface_initial(n))/surface_efficace_interface_initial(n));
          }
        else
          {
            erreur_relative = 0.;
            assert(surface_efficace_interface(n) == 0.);
          }
        double erreur_absolue = std::fabs(surface_efficace_interface(n) - surface_efficace_interface_initial(n));

        moyenne_ecart_absolue_surface += erreur_absolue;
        moyenne_ecart_relatif_surface += erreur_relative;
        maximum_ecart_relatif_surface = erreur_relative > maximum_ecart_relatif_surface ? erreur_relative : maximum_ecart_relatif_surface;
        maximum_ecart_absolue_surface = erreur_absolue > maximum_ecart_absolue_surface ? erreur_absolue : maximum_ecart_absolue_surface;

        n_count_surface += 1;

        if (centre_switch)
          {
            switch_moyenne_ecart_absolue_surface += erreur_absolue;
            switch_moyenne_ecart_relatif_surface += erreur_relative;
            switch_maximum_ecart_relatif_surface = erreur_relative > switch_maximum_ecart_relatif_surface ? erreur_relative : switch_maximum_ecart_relatif_surface;
            switch_maximum_ecart_absolue_surface = erreur_absolue > switch_maximum_ecart_absolue_surface ? erreur_absolue : switch_maximum_ecart_absolue_surface;

            switch_n_count_surface += 1;
          }
      }

      const double vitesse_interface = -vitesse_deplacement_interface(n,0)*normale_deplacement_interface(n,0) - vitesse_deplacement_interface(n,1)*normale_deplacement_interface(n,1) - vitesse_deplacement_interface(n,2)*normale_deplacement_interface(n,2);
      const double surface_interface = surface_efficace_interface(n);

      double delta_volume_cible = volume * (next_indicatrice_ns(i, j, k) - old_indicatrice_ns(i, j, k));
      double delta_volume_interface = vitesse_interface * surface_interface * timestep;

      double erreur_relative = delta_volume_cible == 0. ? delta_volume_interface/volume : std::fabs((delta_volume_interface - delta_volume_cible)/delta_volume_cible);
      double erreur_absolue = std::fabs(delta_volume_interface - delta_volume_cible);

      if (IJK_Interfaces::est_pure(.5*(old_indicatrice_ns(i, j, k) + next_indicatrice_ns(i, j, k))))
        {
          // A failure of this assert suggests a non-zero velocity divergence
          //assert(std::abs(erreur_relative) < 1e-12);
        }
      else
        {
          n_count += 1;
          moyenne_erreur_relative += erreur_relative;
          moyenne_erreur_absolue += erreur_absolue;
          maximum_erreur_relative = erreur_relative > maximum_erreur_relative ? erreur_relative : maximum_erreur_relative;
          maximum_erreur_absolue = erreur_absolue > maximum_erreur_absolue ? erreur_absolue : maximum_erreur_absolue;

          if (centre_switch)
            {
              switch_n_count += 1;
              switch_moyenne_erreur_relative += erreur_relative;
              switch_moyenne_erreur_absolue += erreur_absolue;
              switch_maximum_erreur_relative = erreur_relative > switch_maximum_erreur_relative ? erreur_relative : switch_maximum_erreur_relative;
              switch_maximum_erreur_absolue = erreur_absolue > switch_maximum_erreur_absolue ? erreur_absolue : switch_maximum_erreur_absolue;
            }
        }
    }
  moyenne_erreur_absolue = Process::mp_sum(moyenne_erreur_absolue);
  moyenne_erreur_relative = Process::mp_sum(moyenne_erreur_relative);
  maximum_erreur_absolue = Process::mp_max(maximum_erreur_absolue);
  maximum_erreur_relative = Process::mp_max(maximum_erreur_relative);
  moyenne_ecart_absolue_surface = Process::mp_sum(moyenne_ecart_absolue_surface);
  moyenne_ecart_relatif_surface = Process::mp_sum(moyenne_ecart_relatif_surface);
  maximum_ecart_absolue_surface = Process::mp_max(maximum_ecart_absolue_surface);
  maximum_ecart_relatif_surface = Process::mp_max(maximum_ecart_relatif_surface);
  switch_moyenne_erreur_absolue = Process::mp_sum(switch_moyenne_erreur_absolue);
  switch_moyenne_erreur_relative = Process::mp_sum(switch_moyenne_erreur_relative);
  switch_maximum_erreur_absolue = Process::mp_max(switch_maximum_erreur_absolue);
  switch_maximum_erreur_relative = Process::mp_max(switch_maximum_erreur_relative);
  switch_moyenne_ecart_absolue_surface = Process::mp_sum(switch_moyenne_ecart_absolue_surface);
  switch_moyenne_ecart_relatif_surface = Process::mp_sum(switch_moyenne_ecart_relatif_surface);
  switch_maximum_ecart_absolue_surface = Process::mp_max(switch_maximum_ecart_absolue_surface);
  switch_maximum_ecart_relatif_surface = Process::mp_max(switch_maximum_ecart_relatif_surface);
  n_count = Process::check_int_overflow(Process::mp_sum(n_count));
  n_count_surface = Process::check_int_overflow(Process::mp_sum(n_count_surface));
  switch_n_count = Process::check_int_overflow(Process::mp_sum(switch_n_count));
  switch_n_count_surface = Process::check_int_overflow(Process::mp_sum(switch_n_count_surface));
  assert(n_count_surface == n_count);
  assert(switch_n_count_surface == switch_n_count);
  if (n_count >= 1)
    {
      moyenne_erreur_relative /= n_count;
      moyenne_erreur_absolue /= n_count;
      moyenne_ecart_absolue_surface /= (n_count_surface);
      moyenne_ecart_relatif_surface /= (n_count_surface);
    }
  if (switch_n_count >= 1)
    {
      switch_moyenne_erreur_relative /= switch_n_count;
      switch_moyenne_erreur_absolue /= switch_n_count;
      switch_moyenne_ecart_absolue_surface /= (switch_n_count_surface);
      switch_moyenne_ecart_relatif_surface /= (switch_n_count_surface);
    }

  if (verbosite_surface_efficace_interface == 2)
    {
      Cerr << "Calcul surface efficace interface:" << finl;
      Cerr << " moyenne_erreur: " << moyenne_erreur_absolue << " (relative: " << moyenne_erreur_relative << ") max_erreur: " << maximum_erreur_absolue << " (relative: " << maximum_erreur_relative << ")" << finl;
      Cerr << " moyenne_ecart_surface: " << moyenne_ecart_absolue_surface << " (relatif: " << moyenne_ecart_relatif_surface << ") max_ecart_surface: " << maximum_ecart_absolue_surface << " (relatif: " << maximum_ecart_relatif_surface << ")" << finl;
      Cerr << " switch moyenne_erreur: " << switch_moyenne_erreur_absolue << " (relative: " << switch_moyenne_erreur_relative << ") max_erreur: " << switch_maximum_erreur_absolue << " (relative: " << switch_maximum_erreur_relative << ")" << finl;
      Cerr << " switch moyenne_ecart_surface: " << switch_moyenne_ecart_absolue_surface << " (relatif: " << switch_moyenne_ecart_relatif_surface << ") max_ecart_surface: " << switch_maximum_ecart_absolue_surface << " (relatif: " << switch_maximum_ecart_relatif_surface << ")" << finl;
    }
  else if (verbosite_surface_efficace_interface == 1)
    {
      Cerr << "Calcul surface efficace interface:  mean_rel: " << moyenne_erreur_relative << " max_rel: " << maximum_erreur_relative << "   switch_mean_rel: " << switch_moyenne_erreur_relative << " switch_max_rel: " << switch_maximum_erreur_relative << finl;
    }
  else
    {
      Cerr << "Cut_cell_surface_efficace::imprimer_informations_surface_efficace_interface non reconnu." << finl;
      Process::exit();
    }
}

void Cut_cell_surface_efficace::imprimer_informations_surface_efficace_face(
  int verbosite_surface_efficace_face,
  int iteration_solver_surface_efficace_face,
  double timestep,
  const Cut_field_vector3_double& velocity,
  const IJK_Field_double& old_indicatrice_ns,
  const IJK_Field_double& next_indicatrice_ns,
  const DoubleTabFT_cut_cell_vector3& indicatrice_surfacique_efficace_face,
  const DoubleTabFT_cut_cell_vector3& indicatrice_surfacique_efficace_face_initial)
{
  const Cut_cell_FT_Disc& cut_cell_disc = indicatrice_surfacique_efficace_face.get_cut_cell_disc();
  const IJK_Grid_Geometry& geom = cut_cell_disc.get_splitting().get_grid_geometry();
  const double delta_x = geom.get_constant_delta(0);
  const double delta_y = geom.get_constant_delta(1);
  const int offset = cut_cell_disc.get_splitting().get_offset_local(DIRECTION_K);
  const ArrOfDouble& delta_z_all = geom.get_delta(DIRECTION_K);

  // Impression des resultats
  int n_count = 0;
  int n_count_surface = 0;
  double moyenne_erreur_absolue = 0.;
  double moyenne_erreur_relative = 0.;
  double maximum_erreur_absolue = 0.;
  double maximum_erreur_relative = 0.;
  double moyenne_ecart_absolue_surface = 0.;
  double moyenne_ecart_relatif_surface = 0.;
  double maximum_ecart_absolue_surface = 0.;
  double maximum_ecart_relatif_surface = 0.;
  int switch_n_count = 0;
  int switch_n_count_surface = 0;
  double switch_moyenne_erreur_absolue = 0.;
  double switch_moyenne_erreur_relative = 0.;
  double switch_maximum_erreur_absolue = 0.;
  double switch_maximum_erreur_relative = 0.;
  double switch_moyenne_ecart_absolue_surface = 0.;
  double switch_moyenne_ecart_relatif_surface = 0.;
  double switch_maximum_ecart_absolue_surface = 0.;
  double switch_maximum_ecart_relatif_surface = 0.;
  for (int n = 0; n < cut_cell_disc.get_n_loc(); n++)
    {
      Int3 ijk = cut_cell_disc.get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      bool centre_switch = (IJK_Interfaces::devient_pure(old_indicatrice_ns(i,j,k), next_indicatrice_ns(i,j,k)) || (IJK_Interfaces::devient_diphasique(old_indicatrice_ns(i,j,k), next_indicatrice_ns(i,j,k))));

      const double delta_z = delta_z_all[k + offset];
      const double volume = delta_x * delta_y * delta_z;

      const double fx = delta_y * delta_z * timestep;
      const double fy = delta_x * delta_z * timestep;
      const double fz = delta_x * delta_y * timestep;

      double delta_volume_total = 0;
      for (int num_face = 0; num_face < 6; num_face++)
        {
          int dir = num_face%3;
          int decalage = num_face/3;
          int sign = decalage*2 -1;

          int di = decalage*(dir == 0);
          int dj = decalage*(dir == 1);
          int dk = decalage*(dir == 2);

          double f = select_dir(dir, fx, fy, fz);

          int n_face = cut_cell_disc.get_n_face(num_face, n, i, j, k);
          if (n_face >= 0)
            {
              delta_volume_total -= sign*f*indicatrice_surfacique_efficace_face(n_face, dir)*velocity[dir].diph_l_(n_face);

              {
                double erreur_relative;
                if (indicatrice_surfacique_efficace_face_initial(n_face, dir) != 0.)
                  {
                    erreur_relative = std::fabs((indicatrice_surfacique_efficace_face(n_face, dir) - indicatrice_surfacique_efficace_face_initial(n_face, dir))/indicatrice_surfacique_efficace_face_initial(n_face, dir));
                  }
                else
                  {
                    erreur_relative = 0.;
                    assert(indicatrice_surfacique_efficace_face(n_face, dir) == 0.);
                  }
                double erreur_absolue = std::fabs(indicatrice_surfacique_efficace_face(n_face, dir) - indicatrice_surfacique_efficace_face_initial(n_face, dir));

                moyenne_ecart_absolue_surface += erreur_absolue;
                moyenne_ecart_relatif_surface += erreur_relative;
                maximum_ecart_relatif_surface = erreur_relative > maximum_ecart_relatif_surface ? erreur_relative : maximum_ecart_relatif_surface;
                maximum_ecart_absolue_surface = erreur_absolue > maximum_ecart_absolue_surface ? erreur_absolue : maximum_ecart_absolue_surface;

                //position_typique = (indicatrice_surfacique_efficace_face(n_face, dir) - surface_min)/(surface_max - surface_min);
                n_count_surface += 1;

                if (centre_switch)
                  {
                    switch_moyenne_ecart_absolue_surface += erreur_absolue;
                    switch_moyenne_ecart_relatif_surface += erreur_relative;
                    switch_maximum_ecart_relatif_surface = erreur_relative > switch_maximum_ecart_relatif_surface ? erreur_relative : switch_maximum_ecart_relatif_surface;
                    switch_maximum_ecart_absolue_surface = erreur_absolue > switch_maximum_ecart_absolue_surface ? erreur_absolue : switch_maximum_ecart_absolue_surface;

                    switch_n_count_surface += 1;
                  }
              }
            }
          else
            {
              assert((next_indicatrice_ns(i+di,j+dj,k+dk) == 0) || (next_indicatrice_ns(i+di,j+dj,k+dk) == 1));
              delta_volume_total -= sign*f*next_indicatrice_ns(i+di,j+dj,k+dk)*velocity[dir].pure_(i+di,j+dj,k+dk);
            }

        }

      double delta_volume_cible = volume * (next_indicatrice_ns(i, j, k) - old_indicatrice_ns(i, j, k));

      double erreur_relative = delta_volume_cible == 0. ? delta_volume_total/volume : std::fabs((delta_volume_total - delta_volume_cible)/delta_volume_cible);
      double erreur_absolue = std::fabs(delta_volume_total - delta_volume_cible);

      if (IJK_Interfaces::est_pure(.5*(old_indicatrice_ns(i, j, k) + next_indicatrice_ns(i, j, k))))
        {
          // A failure of this assert suggests a non-zero velocity divergence
          //assert(std::abs(erreur_relative) < 1e-12);
        }
      else
        {
          n_count += 1;
          moyenne_erreur_relative += erreur_relative;
          moyenne_erreur_absolue += erreur_absolue;
          maximum_erreur_relative = erreur_relative > maximum_erreur_relative ? erreur_relative : maximum_erreur_relative;
          maximum_erreur_absolue = erreur_absolue > maximum_erreur_absolue ? erreur_absolue : maximum_erreur_absolue;

          if (centre_switch)
            {
              switch_n_count += 1;
              switch_moyenne_erreur_relative += erreur_relative;
              switch_moyenne_erreur_absolue += erreur_absolue;
              switch_maximum_erreur_relative = erreur_relative > switch_maximum_erreur_relative ? erreur_relative : switch_maximum_erreur_relative;
              switch_maximum_erreur_absolue = erreur_absolue > switch_maximum_erreur_absolue ? erreur_absolue : switch_maximum_erreur_absolue;
            }
        }
    }
  moyenne_erreur_absolue = Process::mp_sum(moyenne_erreur_absolue);
  moyenne_erreur_relative = Process::mp_sum(moyenne_erreur_relative);
  maximum_erreur_absolue = Process::mp_max(maximum_erreur_absolue);
  maximum_erreur_relative = Process::mp_max(maximum_erreur_relative);
  moyenne_ecart_absolue_surface = Process::mp_sum(moyenne_ecart_absolue_surface);
  moyenne_ecart_relatif_surface = Process::mp_sum(moyenne_ecart_relatif_surface);
  maximum_ecart_absolue_surface = Process::mp_max(maximum_ecart_absolue_surface);
  maximum_ecart_relatif_surface = Process::mp_max(maximum_ecart_relatif_surface);
  switch_moyenne_erreur_absolue = Process::mp_sum(switch_moyenne_erreur_absolue);
  switch_moyenne_erreur_relative = Process::mp_sum(switch_moyenne_erreur_relative);
  switch_maximum_erreur_absolue = Process::mp_max(switch_maximum_erreur_absolue);
  switch_maximum_erreur_relative = Process::mp_max(switch_maximum_erreur_relative);
  switch_moyenne_ecart_absolue_surface = Process::mp_sum(switch_moyenne_ecart_absolue_surface);
  switch_moyenne_ecart_relatif_surface = Process::mp_sum(switch_moyenne_ecart_relatif_surface);
  switch_maximum_ecart_absolue_surface = Process::mp_max(switch_maximum_ecart_absolue_surface);
  switch_maximum_ecart_relatif_surface = Process::mp_max(switch_maximum_ecart_relatif_surface);
  n_count = Process::check_int_overflow(Process::mp_sum(n_count));
  n_count_surface = Process::check_int_overflow(Process::mp_sum(n_count_surface));
  switch_n_count = Process::check_int_overflow(Process::mp_sum(switch_n_count));
  switch_n_count_surface = Process::check_int_overflow(Process::mp_sum(switch_n_count_surface));
  if (n_count >= 1)
    {
      moyenne_erreur_relative /= n_count;
      moyenne_erreur_absolue /= n_count;
      moyenne_ecart_absolue_surface /= (n_count_surface);
      moyenne_ecart_relatif_surface /= (n_count_surface);
    }
  if (switch_n_count >= 1)
    {
      switch_moyenne_erreur_relative /= switch_n_count;
      switch_moyenne_erreur_absolue /= switch_n_count;
      switch_moyenne_ecart_absolue_surface /= (switch_n_count_surface);
      switch_moyenne_ecart_relatif_surface /= (switch_n_count_surface);
    }

  if (verbosite_surface_efficace_face == 2)
    {
      Cerr << "Calcul surface efficace face: #iteration: " << iteration_solver_surface_efficace_face << finl;
      Cerr << " moyenne_erreur: " << moyenne_erreur_absolue << " (relative: " << moyenne_erreur_relative << ") max_erreur: " << maximum_erreur_absolue << " (relative: " << maximum_erreur_relative << ")" << finl;
      Cerr << " moyenne_ecart_surface: " << moyenne_ecart_absolue_surface << " (relatif: " << moyenne_ecart_relatif_surface << ") max_ecart_surface: " << maximum_ecart_absolue_surface << " (relatif: " << maximum_ecart_relatif_surface << ")" << finl;
      Cerr << " switch moyenne_erreur: " << switch_moyenne_erreur_absolue << " (relative: " << switch_moyenne_erreur_relative << ") max_erreur: " << switch_maximum_erreur_absolue << " (relative: " << switch_maximum_erreur_relative << ")" << finl;
      Cerr << " switch moyenne_ecart_surface: " << switch_moyenne_ecart_absolue_surface << " (relatif: " << switch_moyenne_ecart_relatif_surface << ") max_ecart_surface: " << switch_maximum_ecart_absolue_surface << " (relatif: " << switch_maximum_ecart_relatif_surface << ")" << finl;
    }
  else if (verbosite_surface_efficace_face == 1)
    {
      Cerr << "Calcul surface efficace face: #" << iteration_solver_surface_efficace_face << " mean_rel: " << moyenne_erreur_relative << " max_rel: " << maximum_erreur_relative << "   switch_mean_rel: " << switch_moyenne_erreur_relative << " switch_max_rel: " << switch_maximum_erreur_relative << finl;
    }
  else
    {
      Cerr << "Cut_cell_surface_efficace::imprimer_informations_surface_efficace_face: verbosite_surface_efficace_face non reconnu." << finl;
      Process::exit();
    }
}

void Cut_cell_surface_efficace::calcul_vitesse_remaillage(double timestep,
                                                          const IJK_Field_double& indicatrice_avant_remaillage,
                                                          const IJK_Field_double& indicatrice_apres_remaillage,
                                                          const IJK_Field_double& indicatrice_fin_pas_de_temps,
                                                          DoubleTabFT_cut_cell_vector3& indicatrice_surfacique_efficace_remaillage_face,
                                                          Cut_field_vector3_double& remeshing_velocity)
{
  remeshing_velocity[0].set_to_uniform_value(0);
  remeshing_velocity[1].set_to_uniform_value(0);
  remeshing_velocity[2].set_to_uniform_value(0);

  assert(&remeshing_velocity[0].get_cut_cell_disc() == &remeshing_velocity[1].get_cut_cell_disc());
  assert(&remeshing_velocity[0].get_cut_cell_disc() == &remeshing_velocity[2].get_cut_cell_disc());
  const Cut_cell_FT_Disc& cut_cell_disc = remeshing_velocity[0].get_cut_cell_disc();
  const IJK_Grid_Geometry& geom = cut_cell_disc.get_splitting().get_grid_geometry();

  const double delta_x = geom.get_constant_delta(0);
  const double delta_y = geom.get_constant_delta(1);
  const double delta_z = geom.get_constant_delta(2);
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(0));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(1));
  assert(cut_cell_disc.get_splitting().get_grid_geometry().is_uniform(2));

  int max_number_of_pass = 10000;
  for (int pass = 0; pass < max_number_of_pass; pass++)
    {
      bool at_least_one_cell_has_changed = false;
      double max_erreur_absolue = 0.;

      // large_phase = 0 : we consider the smaller phase of the cell
      // large_phase = 1 : we consider the larger phase of the cell
      for (int large_phase = 0; large_phase <= 1; large_phase++)
        {
          for (int index = 0; index < cut_cell_disc.get_n_tot(); index++)
            {
              // Pour la grande phase, il est necessaire d'inverser le sens de parcours pour
              // forcer le sens croissant du volume de la phase.
              int index_reverted = (large_phase) ? cut_cell_disc.get_n_tot() - 1 - index : index;
              int n = cut_cell_disc.get_n_from_indicatrice_index(index_reverted);

              Int3 ijk = cut_cell_disc.get_ijk(n);
              int i = ijk[0];
              int j = ijk[1];
              int k = ijk[2];

              if (!cut_cell_disc.within_ghost(i, j, k, 1, 1))
                continue;


              int phase = large_phase ? (indicatrice_fin_pas_de_temps(i, j, k) >= .5) : (indicatrice_fin_pas_de_temps(i, j, k) < .5);

              const double volume = delta_x * delta_y * delta_z;

              const double fx = delta_y * delta_z * timestep;
              const double fy = delta_x * delta_z * timestep;
              const double fz = delta_x * delta_y * timestep;

              double indic_fin_pas_de_temps = (phase == 0) ? 1 - indicatrice_fin_pas_de_temps(i, j, k) : indicatrice_fin_pas_de_temps(i, j, k);

              double area_dt_free = 0;
              double area_dt_total = 0;
              double delta_volume_total = 0;
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

                  DoubleTabFT_cut_cell& remeshing_diph_velocity = (phase == 0) ? remeshing_velocity[dir].diph_v_ : remeshing_velocity[dir].diph_l_;

                  double indic_decale_fin_pas_de_temps = (phase == 0) ? 1 - indicatrice_fin_pas_de_temps(i+di_decale,j+dj_decale,k+dk_decale) : indicatrice_fin_pas_de_temps(i+di_decale,j+dj_decale,k+dk_decale);

                  double f = select_dir(dir, fx, fy, fz);

                  int n_face = cut_cell_disc.get_n_face(num_face, n, i, j, k);
                  int n_decale = cut_cell_disc.get_n(i+di_decale, j+dj_decale, k+dk_decale);
                  if (n_face >= 0)
                    {
                      double surface_efficace = (phase == 0) ? 1 - indicatrice_surfacique_efficace_remaillage_face(n_face, dir) : indicatrice_surfacique_efficace_remaillage_face(n_face, dir);

                      int decale_smaller = (n_decale >= 0) && (indic_decale_fin_pas_de_temps <= indic_fin_pas_de_temps);
                      area_dt_free += (decale_smaller) ? 0. : f*surface_efficace;
                      area_dt_total += f*surface_efficace;
                      delta_volume_total -= sign*f*surface_efficace*remeshing_diph_velocity(n_face);
                      if (!decale_smaller)
                        {
                          assert((pass != 0.) || (remeshing_diph_velocity(n_face) == 0.));
                        }
                    }
                  else
                    {
                      double surface_efficace = (phase == 0) ? 1 - indicatrice_fin_pas_de_temps(i+di,j+dj,k+dk) : indicatrice_fin_pas_de_temps(i+di,j+dj,k+dk);
                      assert((surface_efficace == 0) || (surface_efficace == 1));

                      area_dt_free += f*surface_efficace;
                      area_dt_total += f*surface_efficace;
                      delta_volume_total -= sign*f*surface_efficace*remeshing_velocity[dir].pure_(i+di,j+dj,k+dk);
                      if (surface_efficace > 0)
                        {
                          assert((pass != 0.) || (remeshing_velocity[dir].pure_(i+di,j+dj,k+dk) == 0.));
                        }
                    }
                }

              double indic_apres_remaillage = (phase == 0) ? 1 - indicatrice_apres_remaillage(i, j, k) : indicatrice_apres_remaillage(i, j, k);
              double indic_avant_remaillage = (phase == 0) ? 1 - indicatrice_avant_remaillage(i, j, k) : indicatrice_avant_remaillage(i, j, k);
              double delta_volume_cible = volume * (indic_apres_remaillage - indic_avant_remaillage);

              double erreur_absolue = std::abs((delta_volume_cible - delta_volume_total)/volume);
              max_erreur_absolue = std::max(max_erreur_absolue, erreur_absolue);

              if (erreur_absolue > 1e-12)
                {
                  at_least_one_cell_has_changed = true;

                  double vel;
                  if (area_dt_free == 0.)
                    {
                      Cerr << "Warning: in Cut_cell_surface_efficace::calcul_vitesse_remaillage, there is a cell without free surface to solve the problem. This is most probably due to the fact that all neighbouring cells are smaller. This proves the algorithm is not suitable. The error for this cell is " << erreur_absolue << finl;
                      vel = (delta_volume_cible - delta_volume_total)/area_dt_total; // Dans ce cas, tant pis on utilise toute la surface.
                    }
                  else
                    {
                      vel = (delta_volume_cible - delta_volume_total)/area_dt_free;
                    }

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

                      DoubleTabFT_cut_cell& remeshing_diph_velocity = (phase == 0) ? remeshing_velocity[dir].diph_v_ : remeshing_velocity[dir].diph_l_;

                      double indic_decale_fin_pas_de_temps = (phase == 0) ? 1 - indicatrice_fin_pas_de_temps(i+di_decale,j+dj_decale,k+dk_decale) : indicatrice_fin_pas_de_temps(i+di_decale,j+dj_decale,k+dk_decale);

                      int n_face = cut_cell_disc.get_n_face(num_face, n, i, j, k);
                      int n_decale = cut_cell_disc.get_n(i+di_decale, j+dj_decale, k+dk_decale);
                      if ((area_dt_free == 0.) || (n_decale < 0) || (indic_decale_fin_pas_de_temps > indic_fin_pas_de_temps))
                        {
                          if (n_face >= 0)
                            {
                              double surface_efficace = (phase == 0) ? 1 - indicatrice_surfacique_efficace_remaillage_face(n_face, dir) : indicatrice_surfacique_efficace_remaillage_face(n_face, dir);
                              if (surface_efficace > 0)
                                {
                                  assert((pass != 0.) || (remeshing_diph_velocity(n_face) == 0.));
                                  remeshing_diph_velocity(n_face) += -sign*vel;
                                }
                            }
                          else
                            {
                              double surface_efficace = (phase == 0) ? 1 - indicatrice_fin_pas_de_temps(i+di,j+dj,k+dk) : indicatrice_fin_pas_de_temps(i+di,j+dj,k+dk);
                              assert((surface_efficace == 0) || (surface_efficace == 1));
                              if (surface_efficace > 0)
                                {
                                  assert((pass != 0.) || (remeshing_velocity[dir].pure_(i+di,j+dj,k+dk) == 0.));
                                  remeshing_velocity[dir].pure_(i+di,j+dj,k+dk) += -sign*vel;
                                }
                            }
                        }
                    }
                }
            }
        }


      max_erreur_absolue = Process::mp_max(max_erreur_absolue);
      if (pass >= 5)
        {
          Cerr << "  (Debug) on pass " << pass << " the max absolute error is " << max_erreur_absolue << finl;
        }

      if (Process::mp_max(at_least_one_cell_has_changed))
        {
          remeshing_velocity[0].echange_espace_virtuel(remeshing_velocity[0].ghost());
          remeshing_velocity[1].echange_espace_virtuel(remeshing_velocity[1].ghost());
          remeshing_velocity[2].echange_espace_virtuel(remeshing_velocity[2].ghost());

          if (pass >= max_number_of_pass-1)
            {
              Cerr << "Cut_cell_surface_efficace::calcul_vitesse_remaillage: maximum number of pass reached, yet no convergence." << finl;
              Process::exit();
            }
        }
      else
        {
          Cerr << "Cut_cell_surface_efficace::calcul_vitesse_remaillage: the number of pass is " << pass << finl;
          break;
        }
    }
}

void Cut_cell_surface_efficace::calcul_delta_volume_theorique_bilan(int compo, const DoubleTab& bounding_box_bulles, double timestep,
                                                                    const IJK_Field_double& indicatrice_avant_deformation,
                                                                    const IJK_Field_double& indicatrice_apres_deformation,
                                                                    const IJK_Field_vector3_double& indicatrice_surfacique_efficace_deformation_face,
                                                                    const Cut_field_vector3_double& deformation_velocity,
                                                                    IJK_Field_double& delta_volume_theorique_bilan)
{
  assert(&deformation_velocity[0].get_cut_cell_disc() == &deformation_velocity[1].get_cut_cell_disc());
  assert(&deformation_velocity[0].get_cut_cell_disc() == &deformation_velocity[2].get_cut_cell_disc());
  const Cut_cell_FT_Disc& cut_cell_disc = deformation_velocity[0].get_cut_cell_disc();

  const int ni = delta_volume_theorique_bilan.ni();
  const int nj = delta_volume_theorique_bilan.nj();
  const int nk = delta_volume_theorique_bilan.nk();

  const IJK_Splitting& splitting = delta_volume_theorique_bilan.get_splitting();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  assert(geom.is_uniform(0));
  assert(geom.is_uniform(1));
  assert(geom.is_uniform(2));

  const double delta_x = geom.get_constant_delta(DIRECTION_I);
  const double delta_y = geom.get_constant_delta(DIRECTION_J);
  const double delta_z = geom.get_constant_delta(DIRECTION_K);

  double origin_x = geom.get_origin(DIRECTION_I);
  double origin_y = geom.get_origin(DIRECTION_J);
  double origin_z = geom.get_origin(DIRECTION_K);

  const int offset_x = splitting.get_offset_local(DIRECTION_I);
  const int offset_y = splitting.get_offset_local(DIRECTION_J);
  const int offset_z = splitting.get_offset_local(DIRECTION_K);

  int imin = std::max(0,  (int)((bounding_box_bulles(compo, 0, 0) - origin_x)/delta_x - offset_x - 2));
  int imax = std::min(ni, (int)((bounding_box_bulles(compo, 0, 1) - origin_x)/delta_x - offset_x + 2));
  int jmin = std::max(0,  (int)((bounding_box_bulles(compo, 1, 0) - origin_y)/delta_y - offset_y - 2));
  int jmax = std::min(nj, (int)((bounding_box_bulles(compo, 1, 1) - origin_y)/delta_y - offset_y + 2));
  int kmin = std::max(0,  (int)((bounding_box_bulles(compo, 2, 0) - origin_z)/delta_z - offset_z - 2));
  int kmax = std::min(nk, (int)((bounding_box_bulles(compo, 2, 1) - origin_z)/delta_z - offset_z + 2));

  for (int k = kmin; k < kmax; k++)
    {
      for (int j = jmin; j < jmax; j++)
        {
          for (int i = imin; i < imax; i++)
            {
              bool apres_deformation_pure = ((indicatrice_apres_deformation(i,j,k) == 0.) || (indicatrice_apres_deformation(i,j,k) == 1.));
              if (apres_deformation_pure)
                {
                  assert(delta_volume_theorique_bilan(i,j,k) == 0.);
                }
              else
                {
                  const double fx = delta_y * delta_z * timestep;
                  const double fy = delta_x * delta_z * timestep;
                  const double fz = delta_x * delta_y * timestep;

                  double delta_volume_total[2] = {0};

                  for (int phase = 0 ; phase < 2 ; phase++)
                    {
                      for (int num_face = 0; num_face < 6; num_face++)
                        {
                          int dir = num_face%3;
                          int decalage = num_face/3;
                          int sign = decalage*2 -1;

                          int di = decalage*(dir == 0);
                          int dj = decalage*(dir == 1);
                          int dk = decalage*(dir == 2);

                          double surface_efficace = (phase == 0) ? 1 - indicatrice_surfacique_efficace_deformation_face[dir](i+di,j+dj,k+dk) : indicatrice_surfacique_efficace_deformation_face[dir](i+di,j+dj,k+dk);

                          double f = select_dir(dir, fx, fy, fz);

                          int n = cut_cell_disc.get_n(i, j, k);
                          int n_face = cut_cell_disc.get_n_face(num_face, n, i, j, k);
                          if (n_face >= 0)
                            {
                              const DoubleTabFT_cut_cell& deformation_diph_velocity = (phase == 0) ? deformation_velocity[dir].diph_v_ : deformation_velocity[dir].diph_l_;

                              assert(deformation_diph_velocity(n_face) != 6.3e32);
                              delta_volume_total[phase] -= sign*f*surface_efficace*deformation_diph_velocity(n_face);
                            }
                          else
                            {
                              assert(deformation_velocity[dir].pure_(i+di,j+dj,k+dk) != 6.3e32);
                              delta_volume_total[phase] -= sign*f*surface_efficace*deformation_velocity[dir].pure_(i+di,j+dj,k+dk);
                            }
                        }
                    }


                  // moyenne entre les deux phases de la prediction sur la variation de volume
                  double delta_volume_cible = .5*(delta_volume_total[1] - delta_volume_total[0]); // la variation de la phase 0 est negative

                  assert(delta_volume_theorique_bilan(i,j,k) == 0.);
                  delta_volume_theorique_bilan(i,j,k) = -delta_volume_cible; // signe negatif car on souhaite la variation de la phase 0
                }
            }
        }
    }

  delta_volume_theorique_bilan.echange_espace_virtuel(delta_volume_theorique_bilan.ghost());
}
