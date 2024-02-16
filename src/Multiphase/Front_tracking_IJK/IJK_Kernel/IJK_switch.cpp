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

#include <IJK_switch.h>
#include <math.h>
#include <Probleme_base.h>
#include <SChaine.h>

Implemente_base(Switch_double, "Switch_double", Interprete);

Sortie& Switch_double::printOn(Sortie& s) const
{
  return s;
}
Entree& Switch_double::readOn(Entree& s)
{
  return s;
}

// Interpole source dans cible les champs sont aux faces.
void Switch_double::switch_vit(DoubleTab coeff_i, IntTab Indice_i,
                               DoubleTab coeff_j ,IntTab Indice_j,
                               DoubleTab coeff_k ,IntTab Indice_k,
                               const int dir)
{
  for (int k = 0; k < new_nk_; k++)
    for (int j = 0; j < new_nj_; j++)
      for (int i = 0; i < new_ni_; i++)
        {
          const int i2 = Indice_i[i];
          const int j2 = Indice_j[j];
          const int k2 = Indice_k[k];
          double x = 0.;
          for (int di = 0; di < 2; di++)
            for (int dj = 0; dj < 2; dj++)
              for (int dk = 0; dk < 2; dk++)
                x += coeff_i(i, di) * coeff_j(j, dj) * coeff_k(k, dk) * old_velocity_[dir](i2+di, j2+dj, k2+dk);
          new_velocity_[dir](i,j,k) = x;
        }
}

// Interpole source dans cible les champs sont aux faces.
void Switch_double::switch_vit_direct(SFichier& binary_file)
{
  DoubleTab coeff_i[3], coeff_j[3], coeff_k[3];
  IntTab Indice_i[3], Indice_j[3], Indice_k[3];
  for (int dir = 0; dir < 3; dir++)
    {
      calculer_coords_Vi(dir);
      calculer_coeff(coeff_i[dir],Indice_i[dir],
                     coeff_j[dir],Indice_j[dir],
                     coeff_k[dir],Indice_k[dir]);

    }
  for (int k = 0; k < new_nk_; k++)
    {
      ArrOfFloat tmp((new_ni_+1) * (new_nj_+1) * 3);
      tmp = 0.;
      for (int dir = 0; dir < 3; dir++)
        {
          for (int j = 0; j < new_nj_; j++)
            {
              for (int i = 0; i < new_ni_; i++)
                {
                  const int i2 = Indice_i[dir][i];
                  const int j2 = Indice_j[dir][j];
                  const int k2 = Indice_k[dir][k];
                  double x = 0.;
                  for (int di = 0; di < 2; di++)
                    for (int dj = 0; dj < 2; dj++)
                      for (int dk = 0; dk < 2; dk++)
                        x += coeff_i[dir](i, di) * coeff_j[dir](j, dj) * coeff_k[dir](k, dk) * old_velocity_[dir](i2+di, j2+dj, k2+dk);
                  tmp[(j * (new_ni_+1) + i) * 3 + dir] = (float)x;
                }
            }
        }
      Cerr << "Writing velocity, layer " << k << " / " << new_nk_ << endl;
      binary_file.put(tmp.addr(), tmp.size_array(), 1);
    }
  // Last layer, write zeros.
  // In z velocities are on the wall, x and y velocities are filling values
  // not used.
  {
    int k = new_nk_;
    ArrOfFloat tmp((new_ni_+1) * (new_nj_+1) * 3);
    tmp = 0.;
    Cerr << "Writing velocity, layer " << k << " / " << new_nk_ << endl;
    binary_file.put(tmp.addr(), tmp.size_array(), 1);
  }
}

void Switch_double::prepare_run()
{
  // Recuperation des donnees de maillage
  old_mesh_ = ref_cast(IJK_Splitting, Interprete_bloc::objet_global(old_mesh_name_));
  new_mesh_ = ref_cast(IJK_Splitting, Interprete_bloc::objet_global(new_mesh_name_));
  // On garde en memoire si on a un cas perio :
  perio_k_ = old_mesh_.get_grid_geometry().get_periodic_flag(DIRECTION_K);

  lire_fichier_reprise(nom_reprise_);
}

void Switch_double::set_param(Param& param)
{
  direct_write_=1;
  perio_k_=1;
  param.ajouter("old_ijk_splitting", &old_mesh_name_, Param::REQUIRED);
  param.ajouter("new_ijk_splitting", &new_mesh_name_, Param::REQUIRED);

  param.ajouter("nom_sauvegarde", &nom_sauvegarde_, Param::REQUIRED);
  param.ajouter("nom_reprise", &nom_reprise_, Param::REQUIRED);
  param.ajouter("direct_write", &direct_write_);
}

Entree& Switch_double::interpreter(Entree& is)
{
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades(is);

  prepare_run();

  run();
  return is;
}

void Switch_double::initialise()
{
  const Nom& oldgeomname = old_mesh_.get_grid_geometry().le_nom();

  Cout << "Lecture vitesse initiale dans fichier " << fichier_old_vitesse_ << " timestep= " << timestep_reprise_vitesse_ << finl;
  lire_dans_lata(fichier_old_vitesse_, timestep_reprise_vitesse_, oldgeomname, "VELOCITY",
                 old_velocity_[0], old_velocity_[1], old_velocity_[2]); // fonction qui lit un champ a partir d'un lata .

  old_ni_ = old_mesh_.get_nb_items_local(IJK_Splitting::ELEM, DIRECTION_I);
  old_nj_ = old_mesh_.get_nb_items_local(IJK_Splitting::ELEM, DIRECTION_J);
  old_nk_ = old_mesh_.get_nb_items_local(IJK_Splitting::ELEM, DIRECTION_K);

  new_ni_ = new_mesh_.get_nb_items_local(IJK_Splitting::ELEM, DIRECTION_I);
  new_nj_ = new_mesh_.get_nb_items_local(IJK_Splitting::ELEM, DIRECTION_J);
  new_nk_ = new_mesh_.get_nb_items_local(IJK_Splitting::ELEM, DIRECTION_K);
}

void Switch_double::set_param_reprise(Param& param)
{
  param.ajouter("tinit", &current_time_);
  param.ajouter("terme_acceleration_init", &terme_source_acceleration_);
  // GAB : gabriel.rmairez@cea.fr
  /* Voir reprendre probleme dans IJK_FT.cpp */
  /*
  param.ajouter("forcage", &forcage_);
  param.ajouter("reprise_qdm_source", &qdm_source_);
  param.ajouter("reprise_vap_velocity_tmoy", &vap_velocity_tmoy_);
  param.ajouter("reprise_liq_velocity_tmoy", &liq_velocity_tmoy_);
  param.ajouter("last_source_qdm_update_time", &last_source_qdm_update_time_);
  param.ajouter("offset_list_index_", &offset_list_index_);
  param.ajouter("reprise_size_listes_moyennes_glissantes", &size_listes_source_);
  liste_instants_.resize_array(size_listes_source_);
  liste_vap_dl_.resize_array(size_listes_source_);
  liste_liq_dl_.resize_array(size_listes_source_);
  param.ajouter("reprise_liste_instants", &liste_instants_);
  param.ajouter("reprise_liste_vap_dl", &liste_vap_dl_);
  param.ajouter("reprise_liste_liq_dl", &liste_liq_dl_);
  param.ajouter("reprise_v_target", &reprise_v_target_);
  */
  // fin GAB : gabriel.ramirez@cea.fr
  param.ajouter("fichier_reprise_vitesse", &fichier_old_vitesse_);
  param.ajouter("timestep_reprise_vitesse", &timestep_reprise_vitesse_);
}

void Switch_double::lire_fichier_reprise(const char *fichier_reprise)
{
  // Lecture par tous les processeurs, on retire les commentaires etc...
  LecFicDiffuse_JDD fichier(fichier_reprise);
  Param param(que_suis_je());
  set_param_reprise(param);
  param.lire_avec_accolades(fichier);
  //Cout << "hahahah" <<  forcage_.get_b_flt() << finl;
  //Cout << "hehehe" << forcage_.get_semi_gen() << finl;
  // Appeler ensuite initialize() pour lire les fichiers lata etc...
}


void Switch_double::calculer_coeff(DoubleTab& coeff_i, IntTab& Indice_i,
                                   DoubleTab& coeff_j ,IntTab& Indice_j,
                                   DoubleTab& coeff_k ,IntTab& Indice_k)
{
  coeff_i.resize(new_ni_,2);
  Indice_i.resize(new_ni_);
  for (int i = 0; i < new_ni_; i++)
    {
      // trouver le plus petit i2 tel que coord[i2] >= coord[i]
      int i2 = 0 ;
      for (i2 = -old_x_.ghost(); i2 < old_ni_ + old_x_.ghost() && old_x_[i2] < new_x_[i]; i2++)
        ;
      if (i2 > -old_x_.ghost())
        i2--;
      Indice_i[i] = i2;

      const double delta = old_x_[i2+1] - old_x_[i2];
      coeff_i(i,0) = (old_x_[i2+1] - new_x_[i]) / delta;
      coeff_i(i,1) = 1. - coeff_i(i,0);
    }

  coeff_j.resize(new_nj_,2);
  Indice_j.resize(new_nj_);
  for (int j = 0; j < new_nj_; j++)
    {
      // trouver le plus petit i2 tel que coord[i2] >= coord[i]
      int j2 = 0 ;
      for (j2 = -old_y_.ghost(); j2 < old_nj_ + old_y_.ghost() && old_y_[j2] < new_y_[j]; j2++)
        ;
      if (j2 > -old_y_.ghost())
        j2--;
      Indice_j[j] = j2;

      const double delta = old_y_[j2+1] - old_y_[j2];
      coeff_j(j,0) = (old_y_[j2+1] - new_y_[j]) / delta;
      coeff_j(j,1) = 1. - coeff_j(j,0);
    }
  coeff_k.resize(new_nk_,2);
  Indice_k.resize(new_nk_);
  for (int k = 0; k < new_nk_; k++)
    {
      // trouver le plus petit i2 tel que coord[i2] >= coord[i]
      int k2 = 0 ;
      for (k2 = -old_z_.ghost(); k2 < old_nk_ + old_z_.ghost() && old_z_[k2] < new_z_[k]; k2++)
        ;

      if (k2 > -old_z_.ghost())
        k2--;

      if (k2+1==old_nk_ + old_z_.ghost())
        k2--;
      double delta = old_z_[k2+1] - old_z_[k2];
      if (delta<DMINFLOAT)
        {
          // Les points (old_z_[k2+1] et old_z_[k2]) sont confondus (z_ghost et rempli avec la coord de la paroi...)
          // on se decalle d'un cran :
          k2++;
          delta = old_z_[k2+1] - old_z_[k2];
        }
      Indice_k[k] = k2;
      coeff_k(k,0) = (old_z_[k2+1] - new_z_[k]) / delta;
      coeff_k(k,1) = 1. - coeff_k(k,0);
    }

  // Quand ce n'est pas perio, il n'y a pas de cellules ghost pour le champ aux elem.
  // Donc impossible d'interroger ce champ en [-1] ou en [nk].
  // Dans l'absolu, il faudrait changer la technique d'interpolation dans ce cas, pour tenir compte de la CL.
  // Il est beaucoup plus simple de reduire la precision de la methode d'interpolation a l'ordre 0,
  // en mettant pour cela simplement tout le poids sur l'autre element.
  if (!perio_k_)
    {
      const int kmin = old_mesh_.get_offset_local(DIRECTION_K);
      const bool own_last = (kmin+old_nk_ == old_mesh_.get_grid_geometry().get_nb_elem_tot(DIRECTION_K));
      if (kmin == 0)
        {
          Indice_k[0] = 0;
          coeff_k(0,0) = 1. ; // tous le poids sur l'elem de bord.
          coeff_k(0,1) = 0.; // On va quand meme pas aller prendre le voisin dedans... Et pourquoi pas?
          //                    On pourrait construire une technique d'extrapolation basee sur l'evaluation du gradient entre cell[1] et cell[0].
          //                    Le coeff serait alors negatif.
        }
      if (own_last)
        {
          Indice_k[new_nk_-1] = old_nk_-2; // Je crains que si je ne mets pas n-2, pour le coeff[...,1] il utilisera la case +1 qui va sortir du tableau.
          //                    Cela ouvre aussi l'option a l'extrapolation.
          coeff_k(new_nk_-1,0) = 0;
          coeff_k(new_nk_-1,1) = 1.; // Du coup, c'est bien lui qui doit etre a 1.
        }
    }
}

void Switch_double::calculer_coords(const IJK_Splitting::Localisation loc)
{
  // ancien maillage
  const IJK_Grid_Geometry& old_geom = old_mesh_.get_grid_geometry();
  const double s_dx = old_geom.get_constant_delta(DIRECTION_I);
  const double s_dy = old_geom.get_constant_delta(DIRECTION_J);
  const ArrOfDouble& s_dz = old_geom.get_delta(DIRECTION_K);
  /*
  const int old_offset_i = old_geom.get_offset_local(DIRECTION_I);
  const int old_offset_j = old_geom.get_offset_local(DIRECTION_J);
  const int old_offset_k = old_geom.get_offset_local(DIRECTION_K);
  double s_origin_x = old_geom.get_origin(DIRECTION_I)
                    + ((loc==IJK_Splitting::FACES_J || loc==IJK_Splitting::FACES_K || loc==IJK_Splitting::ELEM) ? (s_dx * 0.5) : 0. ) ;
  double s_origin_y = old_geom.get_origin(DIRECTION_J)
                    + ((loc==IJK_Splitting::FACES_K || loc==IJK_Splitting::FACES_I || loc==IJK_Splitting::ELEM) ? (s_dy * 0.5) : 0. ) ;
  double s_origin_z = old_geom.get_origin(DIRECTION_K)
                    + ((loc==IJK_Splitting::FACES_I || loc==IJK_Splitting::FACES_J || loc==IJK_Splitting::ELEM) ? (s_dz[0] * 0.5) : 0. ) ;
  */
  // Coords of the local origin :
  Vecteur3 old_xyz0 = old_mesh_.get_coords_of_dof(0,0,0,loc);
  double s_origin_x = old_xyz0[DIRECTION_I] ;
  double s_origin_y = old_xyz0[DIRECTION_J] ;
  double s_origin_z = old_xyz0[DIRECTION_K] ;

  old_z_[0] = s_origin_z;
  for (int k=1; k < old_nk_ ; k++)
    old_z_[k] = old_z_[k-1] + 0.5*s_dz[k-1] + 0.5 * s_dz[k];

  // Il faut donner la position des ghosts :
  old_z_[-1]= old_z_[0]-s_dz[0];
  old_z_[old_nk_] = old_z_[old_nk_-1] +s_dz[old_nk_-1];
  if (!perio_k_)
    {
      // Il faut peut-etre corriger la position des ghosts :
      const int kmin = old_mesh_.get_offset_local(DIRECTION_K);
      const bool own_last = (kmin+old_nk_ == old_mesh_.get_grid_geometry().get_nb_elem_tot(DIRECTION_K));
      if (kmin == 0)
        old_z_[-1]= old_z_[0]-0.5*s_dz[0]; /* CL  K_min*/
      if (own_last)
        old_z_[old_nk_] = old_z_[old_nk_-1]+ 0.5*s_dz[old_nk_-1]; /* CL  K_max*/
    }

  /* DIR I */
  old_x_[-1]= s_origin_x-s_dx;
  for (int i =0 ; i <= old_ni_ ; i++)
    old_x_[i] = old_x_[i-1]+s_dx;

  /* DIR J */
  old_y_[-1]=s_origin_y-s_dy;
  for (int j =0 ; j <= old_nj_ ; j++)
    old_y_[j] = old_y_[j-1]+s_dy;

  // nouveau maillage
  const IJK_Grid_Geometry& new_geom = new_mesh_.get_grid_geometry();
  const double c_dx = new_geom.get_constant_delta(DIRECTION_I);
  const double c_dy = new_geom.get_constant_delta(DIRECTION_J);
  const ArrOfDouble& c_dz = new_geom.get_delta(DIRECTION_K);
  /*
  double c_origin_x = new_geom.get_origin(DIRECTION_I)
                    + ((loc==IJK_Splitting::FACES_J || loc==IJK_Splitting::FACES_K || loc==IJK_Splitting::ELEM) ? (c_dx * 0.5) : 0. ) ;
  double c_origin_y = new_geom.get_origin(DIRECTION_J)
                    + ((loc==IJK_Splitting::FACES_K || loc==IJK_Splitting::FACES_I || loc==IJK_Splitting::ELEM) ? (c_dy * 0.5) : 0. ) ;
  double c_origin_z = new_geom.get_origin(DIRECTION_K)
                    + ((loc==IJK_Splitting::FACES_I || loc==IJK_Splitting::FACES_J || loc==IJK_Splitting::ELEM) ? (c_dz[0] * 0.5) : 0. ) ;
  */
  // Coords of the local origin :
  Vecteur3 new_xyz0 = new_mesh_.get_coords_of_dof(0,0,0,loc);
  double c_origin_x = new_xyz0[DIRECTION_I] ;
  double c_origin_y = new_xyz0[DIRECTION_J] ;
  double c_origin_z = new_xyz0[DIRECTION_K] ;

  new_z_[0] = c_origin_z;
  for (int k=1; k < new_nk_ ; k++)
    new_z_[k] = new_z_[k-1] + 0.5*c_dz[k-1] + 0.5 * c_dz[k];

  if (perio_k_)
    {
      // new_z_[-1]= ...ori...-c_dz[0]; new_z_ n'a pas de ghost car pas necessaire !
      // new_z_[new_nk_] = new_z_[new_nk_-1] + c_dz[new_nk_-1];
    }
  else
    {
      // new_z_ n'a pas de ghost car pas necessaire !
      // new_z_[-1]= 0 ; /* CL  K_min*/
      // new_z_[new_nk_] = new_z_[new_nk_-1]+ 0.5*c_dz[new_nk_-1]; /* CL  K_max*/
    }

  /* DIR I */
  new_x_[0]=c_origin_x;
  for (int i =1 ; i < new_ni_ ; i++)
    new_x_[i] = new_x_[i-1]+c_dx;

  /* DIR J */
  new_y_[0]=c_origin_y;
  for (int j =1 ; j < new_nj_ ; j++)
    new_y_[j] = new_y_[j-1]+c_dy;
}

void Switch_double::calculer_coords_elem()
{
  //const IJK_Splitting::Localisation loc = field.get_localisation();
  const IJK_Splitting::Localisation loc=IJK_Splitting::ELEM;
  calculer_coords(loc);
  return;
}

void Switch_double::calculer_coords_Vi(const int dir)
{
  switch(dir)
    {
    case 0:
      calculer_coords(IJK_Splitting::FACES_I);
      break;
    case 1:
      calculer_coords(IJK_Splitting::FACES_J);
      break;
    case 2:
      calculer_coords(IJK_Splitting::FACES_K);
      break;
    default:
      Cerr << "Error in calculer_coords_Vi: wrong dir" << finl;
      Process::exit();
    };
  return;
}

void Switch_double::ecrire_header_lata(const Nom lata_name) // const
{
  Cout << "Dumping lata header and time into " << lata_name << finl;
  dumplata_header(lata_name, new_velocity_[0] /* on passe un champ pour ecrire la geometrie */);
  dumplata_newtime(lata_name, current_time_);
}

// Interpole source dans cible les champs sont aux elements
// Pas besoin de ghost dans new, donc pas de parcours des ghosts.
void Switch_double::switch_scalar_field(const IJK_Field_double& oldf, IJK_Field_double& newf,
                                        DoubleTab coeff_i, IntTab Indice_i,
                                        DoubleTab coeff_j ,IntTab Indice_j,
                                        DoubleTab coeff_k ,IntTab Indice_k) const
{
  const int ni = newf.ni();
  const int nj = newf.nj();
  const int nk = newf.nk();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const int i2 = Indice_i[i];
          const int j2 = Indice_j[j];
          const int k2 = Indice_k[k];
          double x = 0.;
          for (int di = 0; di < 2; di++)
            for (int dj = 0; dj < 2; dj++)
              for (int dk = 0; dk < 2; dk++)
                x += coeff_i(i, di) * coeff_j(j, dj) * coeff_k(k, dk) * oldf(i2+di, j2+dj, k2+dk);
          newf(i,j,k) = x;
        }
}


// Swap J et K et interpole, ecrit le resultat directement dans un fichier lata
void Switch_double::switch_scalar_field_direct(SFichier& binary_file,
                                               const IJK_Field_double& fld,
                                               DoubleTab coeff_i, IntTab Indice_i,
                                               DoubleTab coeff_j ,IntTab Indice_j,
                                               DoubleTab coeff_k ,IntTab Indice_k)
{
  for (int k = 0; k < new_nk_; k++)
    {
      ArrOfFloat tmp(new_ni_ * new_nj_);

      for (int j = 0; j < new_nj_; j++)
        for (int i = 0; i < new_ni_; i++)
          {
            const int i2 = Indice_i[i];
            const int j2 = Indice_j[j];
            const int k2 = Indice_k[k];
            double x = 0.;
            for (int di = 0; di < 2; di++)
              for (int dj = 0; dj < 2; dj++)
                for (int dk = 0; dk < 2; dk++)
                  x += coeff_i(i, di) * coeff_j(j, dj) * coeff_k(k, dk) * fld(i2+di, j2+dj, k2+dk);
            tmp[j * new_ni_ + i] = (float)x;
          }
      Cerr << "Writing field, layer " << k << " / " << new_nk_ << endl;
      binary_file.put(tmp.addr(), tmp.size_array(), 1);
    }
}

int Switch_double::allocate_fields(double& sz_arr)
{
  // Velocity
  if (!direct_write_)
    {
      allocate_velocity(new_velocity_, new_mesh_, 0);
      new_velocity_[0].data() = 0 ;
      new_velocity_[1].data() = 0 ;
      new_velocity_[2].data() = 0 ;
    }
  allocate_velocity(old_velocity_, old_mesh_, 2);

  // Fill with valid floating point data in walls and ghost cells:
  old_velocity_[0].data() = 0 ;
  old_velocity_[1].data() = 0 ;
  old_velocity_[2].data() = 0 ;
  int nb_allocated_arrays = 6; // it's a mix of old and new!

  nb_allocated_arrays += init_thermique();

  sz_arr = old_velocity_[0].data().size_array();
  return nb_allocated_arrays;
}

void Switch_double::write_velocity(const Nom lata_name) const
{
  Cout << "Adding velocities to " << lata_name << finl;
  if (!direct_write_)
    dumplata_vector(lata_name,"VELOCITY", new_velocity_[0], new_velocity_[1], new_velocity_[2], 0);
  else
    {
      // Ecrit a la main les lignes dans le fichier lata maitre:
      if (Process::je_suis_maitre())
        {
          SFichier f;
          const IJK_Grid_Geometry& geom = new_mesh_.get_grid_geometry();
          f.ouvrir(lata_name, ios::app);
          // Attention, peut ne pas tenir dans un int:
          long long n;
          char sz_string[100];
          n = ((long long) geom.get_nb_elem_tot(DIRECTION_I)+1)
              * ((long long) geom.get_nb_elem_tot(DIRECTION_J)+1)
              * ((long long) geom.get_nb_elem_tot(DIRECTION_K)+1);
          snprintf(sz_string, 100, "%lld", n); // Apparemment %lld est la bonne syntaxe pour les long long
          f << "Champ VELOCITY " << (lata_name + Nom(".VELOCITY.data")) <<  " geometrie=" << geom.le_nom() << " size=" << sz_string
            << " localisation=FACES composantes=3 nature=vector" << finl;

        }
    }
}


void Switch_double::run()
{
  Cerr << "IJK_problem_double::run()" << finl;

  // Field allocation:
  double sz_arr;

  int nb_allocated_arrays = allocate_fields(sz_arr);

  Cerr << " Allocating " << nb_allocated_arrays << " arrays, approx total size= "
       << sz_arr * sizeof(double) * nb_allocated_arrays * 9.537E-07 << " MB per core" << finl;


  // Intitialize other fields (velocity, interfaces, rho):
  initialise();

  // ghost cells not needed for new as we never do interpolation on it (or gradient or whatever...)
  old_mesh_.get_local_mesh_delta(DIRECTION_I, 1 /* ghost cells */, old_x_);
  new_mesh_.get_local_mesh_delta(DIRECTION_I, 0 /* ghost cells not needed for new */, new_x_);

  old_mesh_.get_local_mesh_delta(DIRECTION_J, 1 /* ghost cells */, old_y_);
  new_mesh_.get_local_mesh_delta(DIRECTION_J, 0 /* ghost cells */, new_y_);

  old_mesh_.get_local_mesh_delta(DIRECTION_K, 1 /* ghost cells */, old_z_);
  new_mesh_.get_local_mesh_delta(DIRECTION_K, 0 /* ghost cells */, new_z_);

  DoubleTab coeff_i(new_ni_,2);
  DoubleTab coeff_j(new_nj_,2);
  DoubleTab coeff_k(new_nk_,2);

  IntTab Indice_i(new_ni_);
  IntTab Indice_j(new_nj_);
  IntTab Indice_k(new_nk_);

  coeff_i =0;
  coeff_j =0;
  coeff_k =0;
  Indice_i=0;
  Indice_j=0;
  Indice_k=0;

  // We have to fill a single layer of ghost to get a better interpolation at procs boundaries (I'm not sure?)
  old_velocity_[0].echange_espace_virtuel(1 /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_I*/);
  old_velocity_[1].echange_espace_virtuel(1 /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_J*/);
  old_velocity_[2].echange_espace_virtuel(1/*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_K*/);
  // useless??
  //old_rho_.echange_espace_virtuel(1/*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_K*/);

  Nom lata_name = nom_sauvegarde_ + Nom(".lata");
  prepare_thermique(lata_name);

  // force les conditions de bords
  // remplir_gost(); -> I don't understand the use of this method
  if (!perio_k_)
    {
      Cout << "Forcing zero velocities at walls on the old velocity field" << finl;
      force_zero_on_walls(old_velocity_[DIRECTION_K]);
    }

  // Sortie des donnees (sauv + header du lata + interfaces dans le lata)
  ecrire_fichier_reprise(nom_sauvegarde_, false);
  ecrire_header_lata(lata_name);

  // Interpolate fields:
  if (!direct_write_)
    {
      Cout << "direct_write est nul ";
      // Compute and write velocity:
      for (int dir = 0 ; dir < 3 ; dir ++)
        {
          calculer_coords_Vi(dir);
          calculer_coeff(coeff_i,Indice_i,
                         coeff_j,Indice_j,
                         coeff_k,Indice_k);
          switch_vit(coeff_i,Indice_i,
                     coeff_j,Indice_j,
                     coeff_k,Indice_k,
                     dir);
          // Il faut mener les actions pour ajouter les modes ou que sais-ja a process_b et pour remplir correctement forcage
        }

      if (!perio_k_)
        {
          Cout << "Forcing zero velocities at walls on the new velocity field" << finl;
          force_zero_on_walls(new_velocity_[DIRECTION_K]);
        }
      write_velocity(lata_name);

      // Compute and write rho/thermic if needed:
      compute_and_write_extra_fields(lata_name, coeff_i,Indice_i,
                                     coeff_j,Indice_j,
                                     coeff_k,Indice_k);
    }
  else
    {
      SFichier file;
      file.set_bin(1);
      calculer_coords_elem();

      // Direct writing of velocity:
      file.ouvrir(nom_sauvegarde_ + Nom(".lata.VELOCITY.data"));
      switch_vit_direct(file);
      file.close();

      // Direct writing of extra fields:
      compute_and_write_extra_fields_direct(file, coeff_i,Indice_i,
                                            coeff_j,Indice_j,
                                            coeff_k,Indice_k);
    }
}


