/****************************************************************************
* Copyright (c) 2025, CEA
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

#include <Navier_Stokes_FTD_IJK_tools.h>
#include <IJK_Shear_Periodic_helpler.h>
#include <IJK_Navier_Stokes_tools.h>
#include <Boundary_Conditions.h>
#include <IJK_Interfaces.h>

void force_upstream_velocity(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz, double v_imposed, const IJK_Interfaces& interfaces,
                             double nb_diam, int upstream_dir, int gravity_dir, int upstream_stencil)
{
  int dir = 0;
  if (upstream_dir == -1)
    {
      dir = gravity_dir;
      if (dir == -1)
        dir = 0;
    }

  const Domaine_IJK& geom = vx.get_domaine();

  bool perio = geom.get_periodic_flag(dir);

  assert(interfaces.get_nb_bulles_reelles() == 1);
  DoubleTab bounding_box;
  Cerr << "Upstream Velocity - Compute Bounding box" << finl;
  interfaces.calculer_bounding_box_bulles(bounding_box);
  // Calcule la hauteur en x de la permiere bulle et la position de son cdg :
  const double Dbdir = bounding_box(0, dir, 1) - bounding_box(0, dir, 0);
  const double dirb = (bounding_box(0, dir, 1) + bounding_box(0, dir, 0)) / 2.;
  const double ldir = geom.get_domain_length(dir);
  if (nb_diam == 0.)
    nb_diam = (ldir / Dbdir) / 2;
  double dirobj = dirb + nb_diam * Dbdir;

  // L'origine est sur un noeud. Donc que la premiere face en I est sur get_origin(DIRECTION_I)
  const double ddir = geom.get_constant_delta(dir);
  const double origin_dir = geom.get_origin(dir);
  const int offset_dir = geom.get_offset_local(dir);

  // FIXME: If nb_diam is too large it will iterate a lot
  if (perio)
    {
      while (dirobj < origin_dir)
        dirobj += ldir;
      while (dirobj > origin_dir + ldir)
        dirobj -= ldir;
    }

  // On devrait avoir xobj dans le domaine, sinon, on a choisi nb_diam trop grand :
  assert(((dirobj >= origin_dir) && (dirobj <= origin_dir + ldir)));

  const double x2 = (dirobj - origin_dir) / ddir;
  int index_dir = (int) (floor(x2)) - offset_dir; // C'est l'index local, donc potentiellement negatif...
  int ndir;
  switch(dir)
    {
    case 0:
      ndir = vx.ni();
      break;
    case 1:
      ndir = vx.nj();
      break;
    case 2:
      ndir = vx.nk();
      break;
    default:
      ndir = vx.ni();
      break;
    }
  // Cerr << "index_dir " << index_dir << finl;
  if ((index_dir >= 0) && (index_dir < ndir))
    {
      // On est sur le bon proc...
      if (index_dir + upstream_stencil >= ndir)
        {
          // On ne veut pas s'embeter sur 2 procs...
          index_dir = ndir - upstream_stencil;
        }
    }
  else
    return;

  double imposed[3] = { 0., 0., 0. };
  imposed[dir] = v_imposed;
  for (int direction = 0; direction < 3; direction++)
    {
      IJK_Field_double& velocity = select_dir(direction, vx, vy, vz);
      int imin, jmin, kmin;
      int imax, jmax, kmax;
      switch(dir)
        {
        case 0:
          imin = index_dir;
          jmin = 0;
          kmin = 0;
          imax = imin + upstream_stencil;
          jmax = velocity.nj();
          kmax = velocity.nk();
          break;
        case 1:
          imin = 0;
          jmin = index_dir;
          kmin = 0;
          imax = velocity.ni();
          jmax = jmin + upstream_stencil;
          kmax = velocity.nk();
          break;
        case 2:
          imin = 0;
          jmin = 0;
          kmin = index_dir;
          imax = velocity.ni();
          jmax = velocity.nj();
          kmax = kmin + upstream_stencil;
          break;
        default:
          imin = index_dir;
          jmin = 0;
          kmin = 0;
          imax = imin + upstream_stencil;
          jmax = velocity.nj();
          kmax = velocity.nk();
          break;
        }
      for (int k = kmin; k < kmax; k++)
        for (int j = jmin; j < jmax; j++)
          for (int i = imin; i < imax; i++)
            velocity(i, j, k) = imposed[direction];
    }

  Cerr << "Upstream Velocity has been forced" << finl;
}

void force_upstream_velocity_shear_perio(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz, double v_imposed, const IJK_Interfaces& interfaces, double nb_diam, Boundary_Conditions& bc,
                                         double nb_diam_ortho_shear_perio, double Ux0, double Uy0, double Uz0, int epaisseur_maille)
{
  assert(interfaces.get_nb_bulles_reelles() == 1);
  DoubleTab bounding_box;
  Cerr << "Upstream Velocity - Compute Bounding box" << finl;
  interfaces.calculer_bounding_box_bulles(bounding_box);

  // Calcule la hauteur en x et z de la bulle et la position de son cdg :
  const double Dbx = bounding_box(0, 0, 1) - bounding_box(0, 0, 0);
  const double Dbz = bounding_box(0, 2, 1) - bounding_box(0, 2, 0);
  const double xb = (bounding_box(0, 0, 1) + bounding_box(0, 0, 0)) / 2.;
  const double zb = (bounding_box(0, 2, 1) + bounding_box(0, 2, 0)) / 2.;

  const Domaine_IJK& geom = vx.get_domaine();
  double origin_x = geom.get_origin(DIRECTION_I);
  double lx = geom.get_domain_length(DIRECTION_I);
  double origin_z = geom.get_origin(DIRECTION_K);
  double lz = geom.get_domain_length(DIRECTION_K);
  double z_min = zb - (Dbz * nb_diam_ortho_shear_perio) / 2.;
  double z_max = zb + (Dbz * nb_diam_ortho_shear_perio) / 2.;
  double z_min_modulo = z_min;
  double z_max_modulo = z_max;

  bool perio = geom.get_periodic_flag(DIRECTION_K);
  if (perio)
    {
      z_min_modulo = std::fmod(std::fmod(z_min - origin_z, lz) + lz, lz) + origin_z;
      z_max_modulo = std::fmod(std::fmod(z_max - origin_z, lz) + lz, lz) + origin_z;
    }
  // Calcule la hauteur du plan ou sera impose la vitesse par rapport a la bulle reelle
  double xobj = xb + nb_diam * Dbx;
  // Calcule la hauteur du plan ou sera impose la vitesse par rapport a la bulle shear periodique par le shear positif
  double xb_offsetp = xb + std::fmod(std::fmod(IJK_Shear_Periodic_helpler::shear_x_time_, lx) + lx, lx);
  double xobj_offsetp = xb_offsetp + nb_diam * Dbx;
  // Calcule la hauteur du plan ou sera impose la vitesse par rapport a la bulle shear periodique par le shear negatif
  double xb_offsetm = xb - std::fmod(std::fmod(IJK_Shear_Periodic_helpler::shear_x_time_, lx) + lx, lx);
  double xobj_offsetm = xb_offsetm + nb_diam * Dbx;

  perio = geom.get_periodic_flag(DIRECTION_I);
  // on s'assure que la position des plans en x reste dans le domaine physique
  if (perio)
    {
      xobj = std::fmod(std::fmod(xobj - origin_x, lx) + lx, lx) + origin_x;
      xobj_offsetp = std::fmod(std::fmod(xobj_offsetp - origin_x, lx) + lx, lx) + origin_x;
      xobj_offsetm = std::fmod(std::fmod(xobj_offsetm - origin_x, lx) + lx, lx) + origin_x;
    }

  double dx = geom.get_constant_delta(DIRECTION_I);
  double dz = geom.get_constant_delta(DIRECTION_K);
  int offset_i = geom.get_offset_local(DIRECTION_I);
  int offset_k = geom.get_offset_local(DIRECTION_K);

  // position des plans : conversion en indice du tableau NS
  // en shear perio, pas de decoupage sur x. les index_i sont compris entre 0 et ni_tot
  int ni = vy.ni();
  int index_i_offsetp = (int) (round((xobj_offsetp - origin_x) / dx)) - offset_i;
  int index_i_offsetm = (int) (round((xobj_offsetm - origin_x) / dx)) - offset_i;
  int index_i = (int) (round((xobj - origin_x) / dx)) - offset_i;
  // decoupage en z autorise. ATTENTION, indices locaux potentiellement negatifs
  int index_k_min = (int) (round((z_min_modulo - origin_z) / dz)) - offset_k;
  int index_k_max = (int) (((z_max_modulo - origin_z) / dz)) - offset_k;

  for (int direction = 0; direction < 3; direction++)
    {
      IJK_Field_double& velocity = select_dir(direction, vx, vy, vz);
      for (int k = 0; k < velocity.nk(); k++)
        {
          for (int j = 0; j < velocity.nj(); j++)
            {
              for (int i = 0; i < velocity.ni(); i++)
                {
                  bool go_i = false;
                  bool go_k = false;
                  int index_i_real = index_i;

                  // coord Z
                  if (z_min_modulo > z_max_modulo)
                    {
                      // la plan est a cheval sur la frontiere z
                      // Il y donc deux plans distincts ou imposer la vitesse
                      // un pour la bulle "de droite" (k=nk), un pour la bulle "de gauche" (k=0).
                      if (z_min < origin_z)
                        {
                          // la bulle reelle est en 0
                          // la bulle fantome est soumise a un shear positif
                          if (k > index_k_min and index_k_min < velocity.nk())
                            {
                              // indice i du plan pour la bulle fantome
                              index_i_real = index_i_offsetp;
                              go_k = true;

                            }
                          else if (k < index_k_max and index_k_max >= 0)
                            {
                              // indice i du plan pour la bulle reelle
                              index_i_real = index_i;
                              go_k = true;

                            }
                        }
                      else if (z_max > lz + origin_z)
                        {
                          // la bulle reelle est en nk
                          // la bulle fantome est soumise a un shear negatif
                          if (k > index_k_min and index_k_min < velocity.nk())
                            {
                              // indice i du plan pour la bulle reelle
                              index_i_real = index_i;
                              go_k = true;
                            }
                          else if (k < index_k_max and index_k_max >= 0)
                            {
                              // indice i du plan pour la bulle fantome
                              index_i_real = index_i_offsetm;
                              go_k = true;
                            }
                        }
                    }
                  else
                    {
                      // le plan ne traverse pas la frontiere z
                      // un seul plan defini pour la bulle entiere
                      if (z_max_modulo != z_max and z_min_modulo != z_min)
                        {
                          // le plan est entierement contenu dans le domaine etendu
                          // pas a cheval sur la frontiere dz
                          // lechange de compo n a pas encore ete fait
                          if (z_min < origin_z)
                            {
                              if ((index_k_min < velocity.nk() and k > index_k_min) and (k < index_k_max and index_k_max >= 0))
                                {
                                  index_i_real = index_i_offsetp;
                                  go_k = true;
                                }
                            }
                          else if (z_max > lz + origin_z)
                            {
                              if ((index_k_min < velocity.nk() and k > index_k_min) and (k < index_k_max and index_k_max >= 0))
                                {
                                  index_i_real = index_i_offsetm;
                                  go_k = true;
                                }
                            }

                        }
                      else
                        {
                          // la bulle reelle est entierement dans le domaine NS
                          if ((index_k_min < velocity.nk() and k > index_k_min) and (k < index_k_max and index_k_max >= 0))
                            {
                              index_i_real = index_i;
                              go_k = true;
                            }
                        }

                    }

                  // coord X
                  if (i == index_i_real % ni)
                    {
                      // On est sur la ligne du plan ou imposer la vitesse
                      go_i = true;
                    }

                  if (go_i && go_k)
                    {
                      // on impose la vitesse sur une couche de 3 mailles
                      // U(z) = DU_perio * z / Lz
                      // V = W = 0
                      // ATTENTION : doit etre utilisee avec corrections_qdm pour eviter toute deviation
                      // ATTENTION : s assurer que la condition de corrections_qdm est compatible avec ce profil de vitesse
                      if (direction == 0)
                        {
                          double z = dz / 2. + (k + offset_k) * dz;
                          for (int m = 0; m < epaisseur_maille; m++)
                            velocity((i + m) % ni, j, k) = Ux0 + bc.get_dU_perio(bc.get_resolution_u_prime_()) * z / lz;
                        }
                      else if (direction == 2)
                        {
                          for (int m = 0; m < epaisseur_maille; m++)
                            velocity((i + m) % ni, j, k) = Uz0;
                        }
                      else
                        {
                          for (int m = 0; m < epaisseur_maille; m++)
                            velocity((i + m) % ni, j, k) = Uy0;

                        }
                    }

                }
            }
        }
    }

  Cerr << "Upstream Velocity has been forced" << finl;
}
