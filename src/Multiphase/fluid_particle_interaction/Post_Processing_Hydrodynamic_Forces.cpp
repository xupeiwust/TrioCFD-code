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

#include <Post_Processing_Hydrodynamic_Forces.h>

#include <Convection_Diffusion_Temperature_FT_Disc.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Navier_Stokes_FT_Disc.h>
#include <Fluide_Incompressible.h>
#include <Connex_components_FT.h>
#include <Solid_Particle_base.h>
#include <Fluide_Diphasique.h>
#include <Maillage_FT_Disc.h>
#include <Matrice_Dense.h>
#include <Domaine_VDF.h>
#include <TRUSTTab.h>


Implemente_instanciable_sans_constructeur(Post_Processing_Hydrodynamic_Forces,
                                          "Post_Processing_Hydrodynamic_Forces",Objet_U);

Post_Processing_Hydrodynamic_Forces::Post_Processing_Hydrodynamic_Forces()
{
}

Entree& Post_Processing_Hydrodynamic_Forces::readOn(Entree& is)
{
  Param p(que_suis_je());
  set_param(p);
  p.lire_avec_accolades_depuis(is);
  return is;
}

Sortie& Post_Processing_Hydrodynamic_Forces::printOn(Sortie& os) const
{
  Cerr << "Error: Post_Processing_Hydrodynamic_Forces::printOn is not implemented." << finl;
  Process::exit();
  return os;
}

void Post_Processing_Hydrodynamic_Forces::set_param(Param& p)
{
  p.ajouter_flag("compute_hydrodynamic_forces", &is_compute_forces_);
  p.ajouter_flag("compute_heat_transfer", &is_compute_heat_transfer_);
  p.ajouter_flag("compute_stokes_theoretical_forces", &is_compute_stokes_theoretical_forces_);
  p.ajouter_flag("post_process_pressure_fa7", &is_post_process_pressure_fa7_);
  p.ajouter_flag("post_process_pressure_force_fa7", &is_post_process_pressure_force_fa7_);
  p.ajouter_flag("post_process_friction_fa7", &is_post_process_friction_force_fa7_);
  p.ajouter_flag("post_process_stress_tensor_fa7", &is_post_process_stress_tensor_fa7_);
  p.ajouter("interpolation_distance_pressure_P1", &interpolation_distance_pressure_P1_);
  p.ajouter("interpolation_distance_pressure_P2", &interpolation_distance_pressure_P2_);
  p.ajouter("interpolation_distance_temperature_P1", &interpolation_distance_temperature_P1_);
  p.ajouter("interpolation_distance_temperature_P2", &interpolation_distance_temperature_P2_);
  p.ajouter("interpolation_distance_gradU_P1", &interpolation_distance_gradU_P1_);
  p.ajouter("interpolation_distance_gradU_P2", &interpolation_distance_gradU_P2_);
  p.ajouter_non_std("method_pressure_force_computation", (this), Param::REQUIRED);
  p.ajouter_non_std("method_friction_force_computation", (this), Param::REQUIRED);
  p.ajouter_non_std("location_stress_tensor", (this));
}

int Post_Processing_Hydrodynamic_Forces::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{

  if (mot=="method_pressure_force_computation")
    {
      Motcles mots;
      mots.add("trilinear_linear");
      Motcle motbis;
      is >> motbis;
      Cerr << "Reading methode_calcul_force_pressure : " << motbis << finl;
      const int r = mots.search(motbis);
      switch(r)
        {
        case 0:
          method_pressure_force_computation_ = Method_pressure_force_computation::TRILINEAR_LINEAR;
          break;
        default:
          Cerr << "Error " << mots << "was expected whereas " << motbis <<" has been found."<< finl;
          Process::exit();
        }
      return 1;
    }
  else if (mot=="method_friction_force_computation")
    {
      Motcles mots;
      mots.add("trilinear_linear_complete_tensor");
      mots.add("trilinear_linear_projected_tensor");
      Motcle motbis;
      is >> motbis;
      Cerr << "Reading methode_calcul_force_frottements : " << motbis << finl;
      const int r = mots.search(motbis);
      switch(r)
        {
        case 0:
          method_friction_force_computation_ =
            Method_friction_force_computation::TRILINEAR_LINEAR_COMPLET_TENSOR;
          break;
        case 1:
          method_friction_force_computation_ =
            Method_friction_force_computation::TRILINEAR_LINEAR_PROJECTED_TENSOR;
          break;
        default:
          Cerr << "Error " << mots << "was expected whereas " << motbis
               << " has been found." << finl;
          Process::exit();
        }
      return 1;
    }
  else if (mot=="location_stress_tensor")
    {
      Motcles mots;
      mots.add("faces_normale_x");
      mots.add("elements");
      Motcle motbis;
      is >> motbis;
      Cerr << "Reading methode_interpolation : " << motbis << finl;
      const int r = mots.search(motbis);
      switch(r)
        {
        case 0:
          location_stress_tensor_ = Location_stress_tensor::FACES_NORMALE_X;
          break;
        case 1:
          location_stress_tensor_ = Location_stress_tensor::ELEMENTS;
          break;
        default:
          Cerr << "Error " << mots << "was expected whereas " << motbis <<" has been found."<< finl;
          Process::exit();
        }
      return 1;
    }
  else
    {
      Cerr << mot << " is not a keyword understood by " << que_suis_je() <<
           " in lire_motcle_non_standard"<< finl;
      Process::exit();
    }
  return -1;
}


void Post_Processing_Hydrodynamic_Forces::compute_hydrodynamic_forces()
{
  if (!flag_force_computation_)
    {
      Cerr << "Post_Processing_Hydrodynamic_Forces::compute_hydrodynamic_forces"  <<  finl;
      if (dimension!=3)
        Process::exit("Post_Processing_Hydrodynamic_Forces::compute_hydrodynamic_forces"
                      " only implemented for 3D configurations.");
      Transport_Interfaces_FT_Disc& eq_transport = ptr_eq_transport_.valeur();
      const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
      const Maillage_FT_Disc& mesh = eq_transport.maillage_interface_pour_post();
      const int nb_fa7 = mesh.nb_facettes();
      IntVect compo_connexes_fa7(nb_fa7); // Init a zero
      int n = search_connex_components_local_FT(mesh, compo_connexes_fa7);
      int nb_particles_tot=compute_global_connex_components_FT(mesh, compo_connexes_fa7, n);

      resize_and_init_tables(nb_particles_tot);

      if (nb_fa7>0)
        {
          resize_data_fa7(nb_fa7);
          resize_coord_neighbor_fluid_fa7(nb_fa7);

          const ArrOfDouble& fa7_surface = mesh.get_update_surface_facettes();
          const DoubleTab& tab_fa7_normal = mesh.get_update_normale_facettes();
          const DoubleTab& gravity_center_fa7=mesh.get_gravity_center_fa7();

          const OBS_PTR(Convection_Diffusion_Temperature_FT_Disc)& eq_thermal =
            eq_ns.variables_internes().ref_equation_mpoint_;
          bool is_discr_elem_diph=false;
          Convection_Diffusion_Temperature_FT_Disc::Thermal_correction_discretization_method
          thermal_correction_discretization_method=Convection_Diffusion_Temperature_FT_Disc::
                                                   Thermal_correction_discretization_method::P1;
          if (eq_thermal.non_nul())
            {
              thermal_correction_discretization_method=
                eq_thermal.valeur().get_thermal_correction_discretization_method();
              is_discr_elem_diph=thermal_correction_discretization_method==
                                 Convection_Diffusion_Temperature_FT_Disc::
                                 Thermal_correction_discretization_method::ELEM_DIPH;
            }

          if (is_discr_elem_diph)
            list_elem_diph_.resize(nb_fa7,2); // EB (j,0) : num_elem_diph, (j,1) : num compo associee

          const IntTab& particles_eulerian_id_number=eq_ns.get_particles_eulerian_id_number();
          compute_neighbors_coordinates_fluid_fa7(nb_fa7,
                                                  is_discr_elem_diph,
                                                  gravity_center_fa7,
                                                  mesh,
                                                  tab_fa7_normal,
                                                  particles_eulerian_id_number);

          // ----------------------------- Pressure force -----------------------------
          compute_pressure_force_trilinear_linear(nb_fa7,mesh,thermal_correction_discretization_method,
                                                  compo_connexes_fa7, fa7_surface, tab_fa7_normal);
          // ----------------------------- Friction force -----------------------------
          if (method_friction_force_computation_==Method_friction_force_computation::
              TRILINEAR_LINEAR_COMPLET_TENSOR)
            {
              compute_friction_force_complet_tensor(nb_fa7, mesh, compo_connexes_fa7,
                                                    fa7_surface, tab_fa7_normal);
            }
          else if (method_friction_force_computation_==Method_friction_force_computation::
                   TRILINEAR_LINEAR_PROJECTED_TENSOR) // see Butaye et al., Computers and Fluids, 2023.
            {
              compute_friction_force_projected_tensor(nb_fa7, mesh, compo_connexes_fa7,
                                                      fa7_surface, tab_fa7_normal);
            }
        }

      mp_sum_for_each_item(total_pressure_force_);
      mp_sum_for_each_item(total_friction_force_);
      mp_sum_for_each_item(proportion_fa7_ok_UP2_);
      mp_sum_for_each_item(total_surface_interf_);
      mp_sum_for_each_item(prop_P2_fluid_compo_);
      mp_sum_for_each_item(U_P2_moy_);
      mp_sum_for_each_item(Nb_fa7_tot_par_compo_);

      compute_U_P2_moy(nb_particles_tot);
      compute_proportion_fa7_ok_and_is_fluid_P2(nb_particles_tot);

      raise_the_flag();
    }
}

int Post_Processing_Hydrodynamic_Forces::elem_faces_for_interp(int num_elem, int i) const
{
  const Domaine_VF& domain_vf=ref_cast(Domaine_VF,ptr_eq_ns_.valeur().domaine_dis());
  if (num_elem<0 || num_elem>domain_vf.elem_faces().dimension_tot(0))
    {
      Journal() << "WARNING : IMPOSSIBLE to access the face  of the element " << num_elem <<
                " in the direction " << i << finl;
      return -1;
    }
  return domain_vf.elem_faces(num_elem, i);
}

int Post_Processing_Hydrodynamic_Forces::face_voisins_for_interp(int num_face,int i) const
{
  const Domaine_VF& domain_vf=ref_cast(Domaine_VF,ptr_eq_ns_.valeur().domaine_dis());

  if (num_face<0 || num_face > domain_vf.face_voisins().dimension_tot(0))
    {
      Journal() << "WARNING : IMPOSSIBLE to access the neighbor element of the face "
                << num_face << " in the direction " << i << finl;
      return -1;
    }
  return domain_vf.face_voisins(num_face,i);
}

int Post_Processing_Hydrodynamic_Forces::trilinear_interpolation_elem(
  const DoubleTab& valeurs_champ,
  DoubleTab& coord,
  DoubleTab& resu)
{
  return trilinear_interpolation_elem(valeurs_champ, coord, resu, 0,
                                      Convection_Diffusion_Temperature_FT_Disc::Thermal_correction_discretization_method::P1);
}

// EB
/*! @brief Interpolation trilineaire d'un champs "valeurs_champs" aux coordonnees coord a partir des 8 elements (4 en 2D) voisins les plus proches. Valeurs_chaamps contient
 * typiquement le champ de presssion ou le champ de temperature.
 * On utilise egalement cette fonction pour sauvegarder la liste des elments P1 ainsi que l'indiatrice des elements P2
 */
int Post_Processing_Hydrodynamic_Forces::trilinear_interpolation_elem(
  const DoubleTab& valeurs_champ,
  DoubleTab& coord,
  DoubleTab& resu,
  const int is_P2,
  const Convection_Diffusion_Temperature_FT_Disc::Thermal_correction_discretization_method
  thermal_correction_discr_method)
{

  const Domaine_VDF& domain_vdf = ref_cast(Domaine_VDF, ptr_eq_ns_.valeur().domaine_dis());
  int nb_voisins=8; // neighbouring elements
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  Transport_Interfaces_FT_Disc& eq_transport = ptr_eq_transport_.valeur();
  const Maillage_FT_Disc& mesh = eq_transport.maillage_interface();
  const DoubleTab& phase_indicator_function = eq_transport.get_update_indicatrice().valeurs();
  const int nb_fa7=coord.dimension(0);

  const IntTab& particles_eulerian_id_number=eq_ns.get_particles_eulerian_id_number();
  const int sauv_list_P1 = (thermal_correction_discr_method==
                            Convection_Diffusion_Temperature_FT_Disc::
                            Thermal_correction_discretization_method::P1) ? 1 : 0;
  if (sauv_list_P1)
    {
      list_elem_P1_.resize(nb_fa7,2); // In (j,0) we save element id, in (j,1) wa save the particle id
      list_elem_P1_=-1;
    }

  // if P1_ALL, we discretize over all elements used for the interpolation in P1
  const int sauv_list_P1_all = (thermal_correction_discr_method==
                                Convection_Diffusion_Temperature_FT_Disc::
                                Thermal_correction_discretization_method::P1_ALL) ? 1 : 0;
  int nb_elem_P1_all=0;

  if (sauv_list_P1_all)
    list_elem_P1_all_.resize(0,2); // In (j,0) we save element id, in (j,1) wa save the particle id

  const DoubleTab& gravity_center_fa7=mesh.get_gravity_center_fa7();

  for (int fa7=0; fa7<nb_fa7; fa7++)
    {
      if (!mesh.facette_virtuelle(fa7))
        {
          DoubleVect coord_elem_interp(dimension);
          for (int dim=0; dim<dimension; dim++)
            coord_elem_interp(dim)=coord(fa7,dim);
          IntVect neighboring_elements(nb_voisins);
          if (is_P2)
            {
              phase_indicator_function_P2_(fa7)=find_neighboring_elements(
                                                  coord_elem_interp,
                                                  neighboring_elements);
            }
          else
            {
              find_neighboring_elements(coord_elem_interp,neighboring_elements,
                                        sauv_list_P1,
                                        fa7);
              if (sauv_list_P1)
                {
                  int elem_diph=domain_vdf.domaine().chercher_elements(
                                  gravity_center_fa7(fa7,0),
                                  gravity_center_fa7(fa7,1),
                                  gravity_center_fa7(fa7,2));
                  int particle_eulerian_id_number = particles_eulerian_id_number(elem_diph);
                  list_elem_P1_(fa7,1)= particle_eulerian_id_number;
                }
            }
          int access_elem=1;
          for (int i=0; i<nb_voisins; i++)
            {
              if (neighboring_elements(i)<0)
                access_elem= 0;
              if (sauv_list_P1_all)
                {
                  int elem_voisin=neighboring_elements(i);
                  int elem_diph=domain_vdf.domaine().chercher_elements(
                                  gravity_center_fa7(fa7,0),
                                  gravity_center_fa7(fa7,1),
                                  gravity_center_fa7(fa7,2));
                  int particle_eulerian_id_number = particles_eulerian_id_number(elem_diph);
                  if (elem_voisin>=0 && phase_indicator_function(elem_voisin)==1)
                    {
                      list_elem_P1_all_.append_line(elem_voisin,particle_eulerian_id_number);
                      nb_elem_P1_all++;
                    }
                }
            }
          if (access_elem)
            {
              DoubleVect delta_i(dimension);
              delta_i(0) = fabs(domain_vdf.dist_elem(neighboring_elements(0),
                                                     neighboring_elements(1), 0));
              delta_i(1) = fabs(domain_vdf.dist_elem(neighboring_elements(0),
                                                     neighboring_elements(2), 1));
              delta_i(2) = fabs(domain_vdf.dist_elem(neighboring_elements(0),
                                                     neighboring_elements(4), 2));
              DoubleVect coord_elem_0(dimension);
              for (int dim=0; dim<dimension; dim++) coord_elem_0(dim)=domain_vdf.xp(
                                                                          neighboring_elements(0),dim);

              double xfact=fabs((coord_elem_interp(0)-coord_elem_0(0))/delta_i(0));
              double yfact=fabs((coord_elem_interp(1)-coord_elem_0(1))/delta_i(1));
              double zfact=fabs((coord_elem_interp(2)-coord_elem_0(2))/delta_i(2));

              resu(fa7)=(1-zfact)*(
                          (1-yfact)*(
                            (1-xfact)*(valeurs_champ(neighboring_elements(0))) +
                            xfact*(valeurs_champ(neighboring_elements(1)))
                          ) +
                          yfact*(
                            (1-xfact)*(valeurs_champ(neighboring_elements(2))) +
                            xfact*(valeurs_champ(neighboring_elements(3)))
                          )
                        ) +
                        zfact*(
                          (1-yfact)*(
                            (1-xfact)*(valeurs_champ(neighboring_elements(4))) +
                            xfact*(valeurs_champ(neighboring_elements(5)))
                          ) +
                          yfact*(
                            (1-xfact)*(valeurs_champ(neighboring_elements(6))) +
                            xfact*(valeurs_champ(neighboring_elements(7)))
                          )
                        );
            }
          else
            resu(fa7)=-1e15;
        }
    }
  if (sauv_list_P1_all)
    {
      list_elem_P1_all_.resize(nb_elem_P1_all,2);
    }
  return 1;
}


int Post_Processing_Hydrodynamic_Forces::trilinear_interpolation_face(
  const DoubleTab& valeurs_champ,
  DoubleTab& coord,
  DoubleTab& resu)
{
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  const Domaine_VDF& domain_vdf = ref_cast(Domaine_VDF, eq_ns.domaine_dis());
  const Transport_Interfaces_FT_Disc& eq_transport = ptr_eq_transport_.valeur();
  const Maillage_FT_Disc& mesh = eq_transport.maillage_interface();

  int nb_neighbors=8;
  for (int fa7=0; fa7<coord.dimension(0); fa7++)
    {
      if (!mesh.facette_virtuelle(fa7))
        {
          // On recupere les faces voisines : les faces euleriennes qui vont servir a l'interpolation

          DoubleVect coord_elem_interp(dimension); // coordinate of the fa7 gravity center
          for (int dim=0; dim<dimension; dim++)
            coord_elem_interp(dim)=coord(fa7,dim);
          IntTab faces_voisines(dimension,nb_neighbors); // 8 elements voisins au point de coordonnees coord. L'element elem est inclu dedans.
          find_neighboring_faces_xyz(coord_elem_interp,faces_voisines);
          int acces_faces=1;

          for (int dim=0; dim<dimension; dim++)
            {
              for (int i=0; i<nb_neighbors; i++)
                {
                  if (faces_voisines(dim,i)<0)
                    acces_faces= 0; // s'il y a eu un bug dans la recuperation des faces (pbm d'acces : zone de joint, bord), on renvoie 0 et la fa7 n'est pas prise en compte dans le calcul
                }
            }

          if (acces_faces)
            {
              DoubleVect xfact(dimension);
              DoubleVect yfact(dimension);
              DoubleVect zfact(dimension);
              DoubleTab Delta_(dimension);
              Delta_(0)=fabs(domain_vdf.dist_face(faces_voisines(0,0),faces_voisines(0,1),0));
              Delta_(1)=fabs(domain_vdf.dist_face(faces_voisines(0,0),faces_voisines(0,2),1));
              Delta_(2)=fabs(domain_vdf.dist_face(faces_voisines(0,0),faces_voisines(0,4),2));

              for (int dim=0; dim<dimension; dim++)
                {
                  xfact(dim)=fabs((coord(fa7,0)-domain_vdf.xv(faces_voisines(dim,0),0))/Delta_(0));
                  yfact(dim)=fabs((coord(fa7,1)-domain_vdf.xv(faces_voisines(dim,0),1))/Delta_(1));
                  zfact(dim)=fabs((coord(fa7,2)-domain_vdf.xv(faces_voisines(dim,0),2))/Delta_(2));
                }

              for (int dim=0; dim<dimension; dim++)
                {
                  if (xfact(dim)>1)
                    Cerr << "xfact > 1 pour le point " <<coord(fa7,0)<< " " << coord(fa7,1)<< " "
                         << coord(fa7,2) <<finl;
                  if (yfact(dim)>1)
                    Cerr << "yfact > 1 pour le point " <<coord(fa7,0)<< " " << coord(fa7,1)<< " "
                         << coord(fa7,2) <<finl;
                  if (zfact(dim)>1)
                    Cerr << "zfact > 1 pour le point " <<coord(fa7,0)<< " " << coord(fa7,1)<< " "
                         << coord(fa7,2) <<finl;
                }

              for (int dim=0; dim<dimension; dim++)
                {
                  resu(fa7,dim)=(1.-zfact(dim))*(
                                  (1.-yfact(dim))*(
                                    (1.-xfact(dim))*(valeurs_champ(faces_voisines(dim,0))) +
                                    xfact(dim)*(valeurs_champ(faces_voisines(dim,1)))
                                  ) +
                                  yfact(dim)*(
                                    (1.-xfact(dim))*(valeurs_champ(faces_voisines(dim,2))) +
                                    xfact(dim)*(valeurs_champ(faces_voisines(dim,3))))
                                ) +
                                zfact(dim)*(
                                  (1.-yfact(dim))*(
                                    (1.-xfact(dim))*(valeurs_champ(faces_voisines(dim,4))) +
                                    xfact(dim)*(valeurs_champ(faces_voisines(dim,5)))
                                  ) +
                                  yfact(dim)*(
                                    (1.-xfact(dim))*(valeurs_champ(faces_voisines(dim,6))) +
                                    xfact(dim)*(valeurs_champ(faces_voisines(dim,7))))
                                );
                }
            }
          else
            {
              for (int dim=0; dim<dimension; dim++)
                resu(fa7,dim)=-1e15; // de cette maniere, on ne calcule pas la force pour la fa7 pour laquelle on n'a pas acces a P2
            }
        }
    }
  return 1;
}

int Post_Processing_Hydrodynamic_Forces::trilinear_interpolation_gradU_face(
  const DoubleTab& valeurs_champ,
  DoubleTab& coord,
  DoubleTab& resu)
{
  const Transport_Interfaces_FT_Disc& eq_transport = ptr_eq_transport_.valeur();
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  const Maillage_FT_Disc& mesh = eq_transport.maillage_interface();
  const Domaine_VDF& domain_vdf = ref_cast(Domaine_VDF, eq_ns.domaine_dis());
  const Domaine& domain = domain_vdf.domaine();
  const DoubleTab& xv = domain_vdf.xv();
  IntVect faces_elem_interp(2*dimension);
  int nb_neighbors=8;

  for (int fa7=0; fa7<coord.dimension(0); fa7++)
    {
      if (!mesh.facette_virtuelle(fa7))
        {
          const int elem=domain.chercher_elements(coord(fa7,0), coord(fa7,1), coord(fa7,2));
          for (int dim=0; dim<dimension; dim++)
            {
              faces_elem_interp(dim)=elem_faces_for_interp(elem,dim);
              faces_elem_interp(dimension+dim)=elem_faces_for_interp(elem,dimension+dim);
              if (faces_elem_interp(dim)<0 || faces_elem_interp(dimension+dim)<0)
                return 0;
            }

          DoubleTab coord_face(2*dimension,dimension);
          for (int i=0; i<2*dimension; i++)
            {
              for (int dim=0; dim<dimension; dim++)
                {
                  coord_face(i,dim)=domain_vdf.xv(faces_elem_interp(i),dim);
                }
            }

          DoubleVect coord_elem_interp(dimension);
          for (int dim=0; dim<dimension; dim++) coord_elem_interp(dim)=coord(fa7,dim);
          IntVect faces_voisines(nb_neighbors); // 8 elements voisins au point de coordonnees coord. L'element elem est inclu dedans.
          find_neighboring_faces(coord_elem_interp,faces_voisines,0);

          for (int i=0; i<nb_neighbors; i++)
            {
              if (faces_voisines(i)<0)
                return 0;
            }

          DoubleTab gradUx(nb_neighbors, dimension, dimension); // le x signifie que l'on calcule le tenseur en une facette dont la composante de la vitesse est en x
          for (int face=0; face<8; face++) //on calcule le tenseur gradient de la vitesse pour chaque face voisine
            {
              // 														SCHEMA EN 2D
              //  										 ---- --- --- --- --- --- --- --- --- --- -
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //					 										 3
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //  										 ---- --- --- --- -7- -8- --- --- --- --- -
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //											  	  	  	 1 	 x 	 2
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //  									     ---- --- --- --- -5- -6- --- --- --- --- -	  y
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |	 |	  ^
              //					 										 4						  |
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |	 |	  |
              //  										 ---- --- --- --- --- --- --- --- --- --- -   o----->x
              //																					  z
              // Evaluation du gradient sur les faces de normale x (choix arbitraire mais peut tuer une certaine symetrie)
              // Les termes + designent les faces juxtaposees selon +z
              // Les termes - designent les faces juxtaposees selon -z
              // Les termes 57 et 86 pour w designent les faces de normale z secantes aux faces 5,7 et 8,6
              // 		   --					   													   	   										   													   --
              //		   |  u_2-u_1	 		   				  	  u_3-u_4				 						 		  u_x+ - u_x- 	 								   							|
              //		   | ----------								------------	        								---------------																|
              //		   | (x_2-x1)			 					(y_3-y_4)												  (z_x+ - z_x-)	  															|
              //		   |												    																														|
              //  grad U = | 1    v_7-v_8	 v_5 - v_6    			1	v_5-v_7	  v_6-v_8								  1	  v_7+ - v_7- 	     v_8+ - v_8-	   v_5+ - v_5-	      v_6+ - v_6-		|
              // 	  	   | - ( --------- + --------- ) 			- ( ------ + -------)								  - ( --------------- + --------------- + --------------- + --------------- )	|
              //		   | 2    x_7-x_8	 x_5 - x_6   			2	y_5-y_7	  y_6-y_8								  4	  (z_7+ - z_7-)	    (z_8+ - z_8-)	   (z_5+ - z_8-)	  (z_6+ - z_6-)	|
              //		   |																																											|
              // 	  	   | 1	 w57+ - w86+	 w57- - w86-		1	w27-w29	 	w28-w30	    w23-w25      w24-w26	  1   w_57+ - w_57-   w_86+ - w_86-												|
              //		   | - ( ----------- +	 ----------- )		- (	-------- + 	-------- + --------- +  --------- )   - ( ------------ + --------------)	    									|
              //  		   | 2	 x57+ - x86+	 x57- - x86-		4	y27-y29	 	y28-y30	    y23-y25		 y24-y26 	  2   z_57+ - z_57-   z_86+ - z_86-												|
              //		   --					   		   									              									   	   													   --
              //
              // Pour chaque face, il faut donc 30 faces voisines : 1-8, x+, x-, 5+, 5-, 6+, 6-, 7+, 7-, 8+, 8-, 57+, 57-, 86+, 86-
              // Pour plus de facilites, on numerote :  x+:9, x-:10, 5+:11, 5-:12, 6+:13, 6-:14, 7+:15, 7-:16, 8+:17, 8-:18, 57+:19, 57-:20, 86+:21, 86-:22
              // 23-30 : les faces de normales z autour de x. Numerotation de bas en haut, de devant a derriere, de gauche a droite

              int nb_voisinsx=30;
              IntVect voisinsx(nb_voisinsx);
              int num_facex=faces_voisines(face);
              int elem_gauche=face_voisins_for_interp(num_facex, 0);
              int elem_droite=face_voisins_for_interp(num_facex, 1);
              int elem_haut_gauche=face_voisins_for_interp(elem_faces_for_interp(elem_gauche,2+dimension), 1);
              int elem_haut_droite=face_voisins_for_interp(elem_faces_for_interp(elem_droite,2+dimension), 1);
              int elem_bas_gauche=face_voisins_for_interp(elem_faces_for_interp(elem_gauche,2), 0);
              int elem_bas_droite=face_voisins_for_interp(elem_faces_for_interp(elem_droite,2), 0);
              int elem_avant_gauche=face_voisins_for_interp(elem_faces_for_interp(elem_gauche,1+dimension),1);
              int elem_avant_droite=face_voisins_for_interp(elem_faces_for_interp(elem_droite,1+dimension),1);
              int elem_arriere_gauche=face_voisins_for_interp(elem_faces_for_interp(elem_gauche,1),0);
              int elem_arriere_droite=face_voisins_for_interp(elem_faces_for_interp(elem_droite,1),0);

              voisinsx(0)  = elem_faces_for_interp(elem_gauche,0);
              voisinsx(1)  = elem_faces_for_interp(elem_droite,0+dimension);
              voisinsx(2)  = elem_faces_for_interp(elem_avant_gauche,0+dimension);
              voisinsx(3)  = elem_faces_for_interp(elem_arriere_gauche,0+dimension);
              voisinsx(4)  = elem_faces_for_interp(elem_gauche,1);
              voisinsx(5)  = elem_faces_for_interp(elem_droite,1);
              voisinsx(6)  = elem_faces_for_interp(elem_gauche,1+dimension);
              voisinsx(7)  = elem_faces_for_interp(elem_droite,1+dimension);
              voisinsx(8)  = elem_faces_for_interp(elem_haut_gauche,0+dimension);
              voisinsx(9)  = elem_faces_for_interp(elem_bas_gauche,0+dimension);
              voisinsx(10) = elem_faces_for_interp(elem_haut_gauche,1);
              voisinsx(11) = elem_faces_for_interp(elem_bas_gauche,1);
              voisinsx(12) = elem_faces_for_interp(elem_haut_droite,1);
              voisinsx(13) = elem_faces_for_interp(elem_bas_droite,1);
              voisinsx(14) = elem_faces_for_interp(elem_haut_gauche,1+dimension);
              voisinsx(15) = elem_faces_for_interp(elem_bas_gauche,1+dimension);
              voisinsx(16) = elem_faces_for_interp(elem_haut_droite,1+dimension);
              voisinsx(17) = elem_faces_for_interp(elem_bas_droite,1+dimension);
              voisinsx(18) = elem_faces_for_interp(elem_gauche,2+dimension);
              voisinsx(19) = elem_faces_for_interp(elem_gauche,2);
              voisinsx(20) = elem_faces_for_interp(elem_droite,2+dimension);
              voisinsx(21) = elem_faces_for_interp(elem_droite,2);
              voisinsx(22) = elem_faces_for_interp(elem_arriere_gauche,2);
              voisinsx(23) = elem_faces_for_interp(elem_arriere_droite,2);
              voisinsx(24) = elem_faces_for_interp(elem_avant_gauche,2);
              voisinsx(25) = elem_faces_for_interp(elem_avant_droite,2);
              voisinsx(26) = elem_faces_for_interp(elem_arriere_gauche,2+dimension);
              voisinsx(27) = elem_faces_for_interp(elem_arriere_droite,2+dimension);
              voisinsx(28) = elem_faces_for_interp(elem_avant_gauche,2+dimension);
              voisinsx(29) = elem_faces_for_interp(elem_avant_droite,2+dimension);

              for (int i=0; i<nb_voisinsx; i++)
                {
                  if (voisinsx(i)<0)
                    return 0;
                }

              gradUx(face,0,0) = 	  (valeurs_champ(voisinsx(1))  - valeurs_champ(voisinsx(0)))  /
                                    (xv(voisinsx(1),0)  - xv(voisinsx(0),0));

              gradUx(face,0,1) =      (valeurs_champ(voisinsx(2))  - valeurs_champ(voisinsx(3)))  /
                                      (xv(voisinsx(2),1)  - xv(voisinsx(3),1));

              gradUx(face,0,2) =      (valeurs_champ(voisinsx(8))  - valeurs_champ(voisinsx(9)))  /
                                      (xv(voisinsx(8),2)  - xv(voisinsx(9),2));

              gradUx(face,1,0) = 1./2.*((valeurs_champ(voisinsx(6))  - valeurs_champ(voisinsx(7)))  /
                                        (xv(voisinsx(6),0)  - xv(voisinsx(7),0))  + (valeurs_champ(voisinsx(4))  -
                                                                                     valeurs_champ(voisinsx(5)))  / (xv(voisinsx(4),0)  - xv(voisinsx(5),0)));

              gradUx(face,1,1) = 1./2.*((valeurs_champ(voisinsx(4))  - valeurs_champ(voisinsx(6)))  /
                                        (xv(voisinsx(4),1)  - xv(voisinsx(6),1))  + (valeurs_champ(voisinsx(5))  -
                                                                                     valeurs_champ(voisinsx(7)))  / (xv(voisinsx(5),1)  - xv(voisinsx(7),1)));

              gradUx(face,1,2) = 1./4.*((valeurs_champ(voisinsx(14)) - valeurs_champ(voisinsx(15))) /
                                        (xv(voisinsx(14),2) - xv(voisinsx(15),2)) + (valeurs_champ(voisinsx(16)) -
                                                                                     valeurs_champ(voisinsx(17))) / (xv(voisinsx(16),2) - xv(voisinsx(17),2))  +
                                        (valeurs_champ(voisinsx(10)) - valeurs_champ(voisinsx(11))) / (xv(voisinsx(10),2)
                                                                                                       - xv(voisinsx(11),2)) + (valeurs_champ(voisinsx(12)) -
                                                                                                                                valeurs_champ(voisinsx(13))) / (xv(voisinsx(12),2) - xv(voisinsx(13),2)));

              gradUx(face,2,0) = 1./2.*((valeurs_champ(voisinsx(18)) - valeurs_champ(voisinsx(20))) /
                                        (xv(voisinsx(18),0) - xv(voisinsx(20),0)) + (valeurs_champ(voisinsx(19)) -
                                                                                     valeurs_champ(voisinsx(21))) / (xv(voisinsx(19),0) - xv(voisinsx(21),0)));

              gradUx(face,2,1) = 1./4.*((valeurs_champ(voisinsx(26)) - valeurs_champ(voisinsx(28))) /
                                        (xv(voisinsx(26),1) - xv(voisinsx(28),1)) + (valeurs_champ(voisinsx(27)) -
                                                                                     valeurs_champ(voisinsx(29))) / (xv(voisinsx(27),1) - xv(voisinsx(29),1))  +
                                        (valeurs_champ(voisinsx(22)) - valeurs_champ(voisinsx(24))) / (xv(voisinsx(22),1) -
                                                                                                       xv(voisinsx(24),1)) + (valeurs_champ(voisinsx(23)) - valeurs_champ(voisinsx(25))) /
                                        (xv(voisinsx(23),1) - xv(voisinsx(25),1)));

              gradUx(face,2,2) = 1./2.*((valeurs_champ(voisinsx(18)) - valeurs_champ(voisinsx(19))) /
                                        (xv(voisinsx(18),2) - xv(voisinsx(19),2)) + (valeurs_champ(voisinsx(20)) -
                                                                                     valeurs_champ(voisinsx(21))) / (xv(voisinsx(20),2) - xv(voisinsx(21),2)));
            }

          double xfact;
          double yfact;
          double zfact;

          DoubleTab Delta_x(dimension);

          Delta_x(0)=fabs(domain_vdf.dist_face(faces_voisines(0),faces_voisines(1),0));
          Delta_x(1)=fabs(domain_vdf.dist_face(faces_voisines(0),faces_voisines(2),1));
          Delta_x(2)=fabs(domain_vdf.dist_face(faces_voisines(0),faces_voisines(4),2));

          xfact=fabs((coord(fa7,0)-coord_face(0,0))/Delta_x(0));
          yfact=fabs((coord(fa7,1)-coord_face(0,1))/Delta_x(1));
          zfact=fabs((coord(fa7,2)-coord_face(0,2))/Delta_x(2));

          for (int i=0; i<dimension; i++)
            {
              for (int j=0; j<dimension; j++)
                {
                  resu(fa7,i,j)=(1-zfact)*(
                                  (1-yfact)*( (1-xfact)*(gradUx(0,i,j)) + xfact*(gradUx(1,i,j)) ) +
                                  yfact*( (1-xfact)*(gradUx(2,i,j)) + xfact*(gradUx(3,i,j)) )
                                ) +
                                zfact*(
                                  (1-yfact)*( (1-xfact)*(gradUx(4,i,j)) + xfact*(gradUx(5,i,j)) ) +
                                  yfact*( (1-xfact)*(gradUx(6,i,j)) + xfact*(gradUx(7,i,j))) );
                }
            }
        }
    }
  return 1;
}

int Post_Processing_Hydrodynamic_Forces::trilinear_interpolation_gradU_elem_P1(
  const DoubleTab& valeurs_champ,
  DoubleTab& coord,
  DoubleTab& resu)
{
  Transport_Interfaces_FT_Disc& eq_transport = ptr_eq_transport_.valeur();
  const DoubleTab& phase_indicator_function = eq_transport.get_update_indicatrice().valeurs();
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  const IntTab& particles_eulerian_id_number=eq_ns.get_particles_eulerian_id_number();
  const Fluide_Diphasique& two_phase_fluid = eq_ns.fluide_diphasique();
  const Maillage_FT_Disc& mesh = eq_transport.maillage_interface();
  const Domaine_VDF& domain_vdf = ref_cast(Domaine_VDF, eq_ns.domaine_dis());
  const DoubleTab& xv = domain_vdf.xv();
  int nb_voisins=8;
  const int id_fluid_phase=two_phase_fluid.get_id_fluid_phase();
  const double mu_f =two_phase_fluid.fluide_phase(id_fluid_phase).
                     viscosite_dynamique().valeurs()(0, 0);
  const double mu_p = two_phase_fluid.fluide_phase(1-id_fluid_phase).
                      viscosite_dynamique().valeurs()(0, 0);

  for (int fa7=0; fa7<coord.dimension(0); fa7++)
    {
      if (!mesh.facette_virtuelle(fa7))
        {
          // find the 8 neighboring elements where the stress tensor will be computed
          DoubleVect coord_elem_interp(dimension);
          for (int dim=0; dim<dimension; dim++) coord_elem_interp(dim)=coord(fa7,dim);
          IntVect elem_voisins(nb_voisins);
          find_neighboring_elements(coord_elem_interp,elem_voisins);

          for (int i=0; i<nb_voisins; i++)
            {
              if (elem_voisins(i)<0)
                return 0;
            }

          // compute distance between neighboring meshes in all direction to compute
          // interpolation coefficients
          DoubleVect coord_elem_0(dimension);
          DoubleVect delta_i(dimension);
          delta_i(0) = fabs(domain_vdf.dist_elem(elem_voisins(0), elem_voisins(1), 0));
          delta_i(1) = fabs(domain_vdf.dist_elem(elem_voisins(0), elem_voisins(3), 1));
          delta_i(2) = fabs(domain_vdf.dist_elem(elem_voisins(0), elem_voisins(5), 2));

          for (int dim=0; dim<dimension; dim++)
            coord_elem_0(dim)=domain_vdf.xp(elem_voisins(0),dim);

          double xfact=fabs((coord_elem_interp(0)-coord_elem_0(0))/delta_i(0));
          double yfact=fabs((coord_elem_interp(1)-coord_elem_0(1))/delta_i(1));
          double zfact=fabs((coord_elem_interp(2)-coord_elem_0(2))/delta_i(2));

          // compute the velocity gradient tensor for each neighboring element
          DoubleTab gradU(nb_voisins, dimension, dimension);
          for (int elem=0; elem<nb_voisins; elem++)
            {

              // 														SCHEMA EN 2D
              //  					FRONT				 ---- ---- --- --- --- --- --- --- ---     FRONT
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //														 9	 10
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //  										 ---- --- --- -5- -3- -6- --- --- ---
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //	  	  	    	  	LEFT  	  	  	  	  	  	     1 x 2              		RIGHT
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //  									     ---- --- --- -7- -4- -8- --- --- --- 	      			y
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |	 	  			^
              //												         11  12									|
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |	 	  			|
              //  					BACK				 ---- --- --- --- --- --- --- --- --- ---    BACK    o----->x
              //																					 			 z
              // Evaluation du gradient sur les faces de normale x (choix arbitraire mais peut tuer une certaine symetrie)
              // Les termes + designent les faces juxtaposees selon +z
              // Les termes - designent les faces juxtaposees selon -z
              // Les termes 57 et 86 pour w designent les faces de normale z secantes aux faces 5,7 et 8,6
              //
              //   ON CALCULE mu*gradU et pas gradU
              //   Pour plus de lisibilite, mu n'est pas reecrit pour chaque composante (mais sans ce mu_h local, la decomposition n'as plus d'interet)
              //   ON UTILISE LE MU HARMONIQUE LOCAL, IE : LE MU_HARMONIQUE AUX ARETES
              //
              // 		   --					   													   	   																															   								   --
              //		   |   u_2-u_1	 		   				                                   1  1  u_9-u_1	 u_1-u_11   1  u_10-u_2	   u_2-u_12	          								1     1   u_1+ - u_1    u_1 - u_1-	   1   u_2+ - u_2   u_2 - u_2-      |
              //		   | ----------							                                   -( - (--------- + --------)+ - (--------- + ---------) )	      								-	( - ( ----------- + ---------- ) + - ( ----------- + ---------- ) ) |
              //		   |   x_2-x1			 				                                   2  2  y_9-y_1	 y_1-y_11   2  y_10-y_2	   y_2-y_12			  								2     2   z_1+ - z_1	z_1 - z_1- 	   2   z_2+ - z_2   z_2 - z_2-      |
              //		   |												    																																														|
              //  grad U = | 1   1   v_5-v_3	v_3 - v_6    1  v_7-v_4	   v_4-v_8 	                       v_3-v_4								   				  								1     1   v_3+ - v_3    v_3 - v_3-	   1   v_4+ - v_4    v_4 - v_4-		|
              // 	  	   | - ( - ( -------- + ---------)+ - (-------- + --------) ) 			 	       ------- 						    					  								-	( - ( ----------- + ---------- ) + - ( ----------- + ---------- ) )	|
              //		   | 2   2   x_5-x_3    x_3 - x_6    2  x_7-x_4	   x_4-x_8	                       y_3-y_4						     					  								2     2   z_3+ - z_3	z_3 - z_3- 	   2   z_4+ - z_4    z_4 - z_4- 	|
              //		   |																																																											|
              // 	  	   | 1	 1  w57+ - w34+	 w34+ - w86+	 1 w57- - w34-	 w34- - w86-       1   1   w_1112+ - w_34+   w_34+ - w_910+	  1     w_1112- - w_34-   w_34- - w_910-	      									w_34+ - w_34-							|
              //		   | - ( - (---------- + ----------- ) + -(---------- + ------------) )    - ( - ( --------------- + -------------- ) + - ( --------------- + -------------- ) )					    				-------------	    					|
              //  		   | 2	 2  x57+ - x34+	 x34+ - x86+     2 x57- - x34-	 x34- - x86-       2   2   y_1112+ - y_34+	 y_34+ - y_910+   2     y_1112- - y_34-   y_34- - y_910- 	  		  	     					    z_34+ - z_34-							|
              //		   --					   		   									              									 																						   								   --
              //
              // Pour chaque face, il faut donc 30 faces voisines : 1-12, 1+, 1-, 2+, 2-, 3+, 3-, 4+, 4-, 34+, 34-, 57+, 57-, 86+, 86-, 1112+, 1112-, 910+, 910-
              // Pour plus de facilites, on numerote :  1+:13, 1-:14, 2+:15, 2-:16, 3+:17, 3-:18, 4+:19, 4-:20, 34+:21, 34-:22, 57+:23, 57-:24, 86+:25, 86-:26, 1112+:27, 1112-:28, 910+:29, 910-:30
              // On identifie les faces voisines pour le calcul des composantes du tenseur

              int nb_faces_voisines=30;
              IntVect les_faces_voisines(nb_faces_voisines);
              int elem_=elem_voisins(elem);
              int elem_gauche = face_voisins_for_interp(elem_faces_for_interp(elem_,0),0);
              int elem_droite = face_voisins_for_interp(elem_faces_for_interp(elem_,0+dimension),1);
              int elem_haut=face_voisins_for_interp(elem_faces_for_interp(elem_,2+dimension), 1);
              int elem_bas=face_voisins_for_interp(elem_faces_for_interp(elem_,2), 0);
              int elem_avant=face_voisins_for_interp(elem_faces_for_interp(elem_,1+dimension),1);
              int elem_arriere=face_voisins_for_interp(elem_faces_for_interp(elem_,1),0);

              les_faces_voisines(0)=elem_faces_for_interp(elem_,0);
              les_faces_voisines(1)=elem_faces_for_interp(elem_,0+dimension);
              les_faces_voisines(2)=elem_faces_for_interp(elem_,1);
              les_faces_voisines(3)=elem_faces_for_interp(elem_,1+dimension);
              les_faces_voisines(4)=elem_faces_for_interp(elem_gauche,1+dimension);
              les_faces_voisines(5)=elem_faces_for_interp(elem_droite,1+dimension);
              les_faces_voisines(6)=elem_faces_for_interp(elem_gauche,1);
              les_faces_voisines(7)=elem_faces_for_interp(elem_droite,1);
              les_faces_voisines(8)=elem_faces_for_interp(elem_avant,0);
              les_faces_voisines(9)=elem_faces_for_interp(elem_avant,0+dimension);
              les_faces_voisines(10)=elem_faces_for_interp(elem_arriere,0);
              les_faces_voisines(11)=elem_faces_for_interp(elem_arriere,0+dimension);
              les_faces_voisines(12)=elem_faces_for_interp(elem_haut,0);
              les_faces_voisines(13)=elem_faces_for_interp(elem_bas,0);
              les_faces_voisines(14)=elem_faces_for_interp(elem_haut,0+dimension);
              les_faces_voisines(15)=elem_faces_for_interp(elem_bas,0+dimension);
              les_faces_voisines(16)=elem_faces_for_interp(elem_haut,1+dimension);
              les_faces_voisines(17)=elem_faces_for_interp(elem_bas,1+dimension);
              les_faces_voisines(18)=elem_faces_for_interp(elem_haut,1);
              les_faces_voisines(19)=elem_faces_for_interp(elem_bas,1);
              les_faces_voisines(20)=elem_faces_for_interp(elem_,2);
              les_faces_voisines(21)=elem_faces_for_interp(elem_,2+dimension);
              les_faces_voisines(22)=elem_faces_for_interp(elem_gauche,2);
              les_faces_voisines(23)=elem_faces_for_interp(elem_gauche,2+dimension);
              les_faces_voisines(24)=elem_faces_for_interp(elem_droite,2);
              les_faces_voisines(25)=elem_faces_for_interp(elem_droite,2+dimension);
              les_faces_voisines(26)=elem_faces_for_interp(elem_arriere,2);
              les_faces_voisines(27)=elem_faces_for_interp(elem_arriere,2+dimension);
              les_faces_voisines(28)=elem_faces_for_interp(elem_avant,2);
              les_faces_voisines(29)=elem_faces_for_interp(elem_avant,2+dimension);

              for (int i=0; i<nb_faces_voisines; i++)
                {
                  if (les_faces_voisines(i)<0)
                    return 0;
                }

              int particle_id= particles_eulerian_id_number(elem);

              gradU(elem,0,0) = ( (mu_p*mu_f/(mu_f-phase_indicator_function(elem_voisins(elem))*
                                              (mu_f-mu_p)))*(valeurs_champ(les_faces_voisines(1))  -
                                                             valeurs_champ(les_faces_voisines(0)))  /
                                  (xv(les_faces_voisines(1),0)  - xv(les_faces_voisines(0),0)));

              gradU(elem,0,1) = 1./2.*  ( 1./2.*( compute_viscosity_edges_sphere(les_faces_voisines(8),
                                                                                 les_faces_voisines(0),particle_id)*(valeurs_champ(les_faces_voisines(8))  -
                                                                                                                     valeurs_champ(les_faces_voisines(0))) / (xv(les_faces_voisines(8),1)  -
                                                                                                                         xv(les_faces_voisines(0),1))   + compute_viscosity_edges_sphere(
                                                    les_faces_voisines(0),les_faces_voisines(10),particle_id)*(valeurs_champ(
                                                                                                                 les_faces_voisines(0))  - valeurs_champ(les_faces_voisines(10)))  /
                                                  (xv(les_faces_voisines(0),1) - xv(les_faces_voisines(10),1)))
                                          + 1./2.*( compute_viscosity_edges_sphere(les_faces_voisines(9),
                                                                                   les_faces_voisines(1),particle_id)*(valeurs_champ(les_faces_voisines(9))  -
                                                                                                                       valeurs_champ(les_faces_voisines(1))) / (xv(les_faces_voisines(9),1)  -
                                                                                                                           xv(les_faces_voisines(1),1))   + compute_viscosity_edges_sphere(
                                                      les_faces_voisines(1),les_faces_voisines(11),particle_id)*(valeurs_champ(
                                                                                                                   les_faces_voisines(1))  - valeurs_champ(les_faces_voisines(11)))  /
                                                    (xv(les_faces_voisines(1),1) - xv(les_faces_voisines(11),1))) );

              gradU(elem,0,2) = 1./2.*  ( 1./2.*( compute_viscosity_edges_sphere(les_faces_voisines(12),
                                                                                 les_faces_voisines(0),particle_id)*(valeurs_champ(les_faces_voisines(12)) -
                                                                                                                     valeurs_champ(les_faces_voisines(0))) / (xv(les_faces_voisines(12),2) -
                                                                                                                         xv(les_faces_voisines(0),2))   + compute_viscosity_edges_sphere(
                                                    les_faces_voisines(0),les_faces_voisines(13),particle_id)*(valeurs_champ(
                                                                                                                 les_faces_voisines(0)) - valeurs_champ(les_faces_voisines(13)))  /
                                                  (xv(les_faces_voisines(0),2)  - xv(les_faces_voisines(13),2)))
                                          + 1./2.*( compute_viscosity_edges_sphere(les_faces_voisines(14),
                                                                                   les_faces_voisines(1),particle_id)*(valeurs_champ(les_faces_voisines(14)) -
                                                                                                                       valeurs_champ(les_faces_voisines(1))) / (xv(les_faces_voisines(14),2) -
                                                                                                                           xv(les_faces_voisines(1),2))   + compute_viscosity_edges_sphere(
                                                      les_faces_voisines(1),les_faces_voisines(15),particle_id)*(valeurs_champ(
                                                                                                                   les_faces_voisines(1)) - valeurs_champ(les_faces_voisines(15)))  /
                                                    (xv(les_faces_voisines(1),2)  - xv(les_faces_voisines(15),2))) );

              gradU(elem,1,0) = 1./2.*  ( 1./2.*( compute_viscosity_edges_sphere(les_faces_voisines(4),
                                                                                 les_faces_voisines(2),particle_id)*(valeurs_champ(les_faces_voisines(4))  -
                                                                                                                     valeurs_champ(les_faces_voisines(2)))  / (xv(les_faces_voisines(4),0)  -
                                                                                                                         xv(les_faces_voisines(2),0))    + compute_viscosity_edges_sphere(
                                                    les_faces_voisines(3),les_faces_voisines(5),particle_id)*(valeurs_champ(
                                                                                                                les_faces_voisines(2))  - valeurs_champ(les_faces_voisines(5)))   /
                                                  (xv(les_faces_voisines(2),0)   - xv(les_faces_voisines(5),0)))
                                          + 1./2.*( compute_viscosity_edges_sphere(les_faces_voisines(6),
                                                                                   les_faces_voisines(3),particle_id)*(valeurs_champ(les_faces_voisines(6))  -
                                                                                                                       valeurs_champ(les_faces_voisines(3)))  / (xv(les_faces_voisines(6),0)  -
                                                                                                                           xv(les_faces_voisines(3),0))    + compute_viscosity_edges_sphere(
                                                      les_faces_voisines(3),les_faces_voisines(7),particle_id)*(valeurs_champ(
                                                                                                                  les_faces_voisines(3))  - valeurs_champ(les_faces_voisines(7)))   /
                                                    (xv(les_faces_voisines(3),0)   - xv(les_faces_voisines(7),0))) );

              gradU(elem,1,1) = ( (mu_p*mu_f/(mu_f-phase_indicator_function(elem_voisins(elem))*
                                              (mu_f-mu_p)))*(valeurs_champ(les_faces_voisines(2))  - valeurs_champ(
                                                               les_faces_voisines(3)))  / (xv(les_faces_voisines(2),1)  - xv(les_faces_voisines(3),1)));

              gradU(elem,1,2) = 1./2.*  ( 1./2.*( compute_viscosity_edges_sphere(
                                                    les_faces_voisines(16),les_faces_voisines(2),particle_id)*(valeurs_champ(
                                                                                                                 les_faces_voisines(16)) - valeurs_champ(les_faces_voisines(2))) / (
                                                    xv(les_faces_voisines(16),2) - xv(les_faces_voisines(2),2)) +
                                                  compute_viscosity_edges_sphere(les_faces_voisines(2),les_faces_voisines(17),
                                                                                 particle_id)*(valeurs_champ(les_faces_voisines(2)) - valeurs_champ(
                                                                                                 les_faces_voisines(17)))  / (xv(les_faces_voisines(2),2) -
                                                                                                                              xv(les_faces_voisines(17),2)))
                                          + 1./2.*( compute_viscosity_edges_sphere(les_faces_voisines(18),
                                                                                   les_faces_voisines(3),particle_id)*(valeurs_champ(les_faces_voisines(18))
                                                                                                                       - valeurs_champ(les_faces_voisines(3))) / (xv(les_faces_voisines(18),2) -
                                                                                                                           xv(les_faces_voisines(3),2))   + compute_viscosity_edges_sphere(
                                                      les_faces_voisines(3),les_faces_voisines(19),particle_id)*(valeurs_champ(
                                                                                                                   les_faces_voisines(3)) - valeurs_champ(les_faces_voisines(19)))  /
                                                    (xv(les_faces_voisines(3),2)  - xv(les_faces_voisines(19),2))));

              gradU(elem,2,0) = 1./2.* ( 1./2.*( compute_viscosity_edges_sphere(les_faces_voisines(22),
                                                                                les_faces_voisines(20),particle_id)*(valeurs_champ(les_faces_voisines(22)) -
                                                                                                                     valeurs_champ(les_faces_voisines(20))) / (xv(les_faces_voisines(22),0) -
                                                                                                                         xv(les_faces_voisines(20),0))   + compute_viscosity_edges_sphere(
                                                   les_faces_voisines(20),les_faces_voisines(24),particle_id)*(valeurs_champ(
                                                                                                                 les_faces_voisines(20)) - valeurs_champ(les_faces_voisines(24)))  /
                                                 (xv(les_faces_voisines(20),0)  - xv(les_faces_voisines(24),0)))
                                         + 1./2.*( compute_viscosity_edges_sphere(les_faces_voisines(23),
                                                                                  les_faces_voisines(21),particle_id)*(valeurs_champ(les_faces_voisines(23)) -
                                                                                                                       valeurs_champ(les_faces_voisines(21))) / (xv(les_faces_voisines(23),0) -
                                                                                                                           xv(les_faces_voisines(21),0)) + compute_viscosity_edges_sphere(
                                                     les_faces_voisines(21),les_faces_voisines(25),particle_id)*(valeurs_champ(
                                                                                                                   les_faces_voisines(21)) - valeurs_champ(les_faces_voisines(25)))  /
                                                   (xv(les_faces_voisines(21),0)  - xv(les_faces_voisines(25),0))) );

              gradU(elem,2,1) = 1./2.*  ( 1./2.*( compute_viscosity_edges_sphere(les_faces_voisines(26),
                                                                                 les_faces_voisines(20),particle_id)*(valeurs_champ(les_faces_voisines(26)) -
                                                                                                                      valeurs_champ(les_faces_voisines(20))) / (xv(les_faces_voisines(26),1) -
                                                                                                                          xv(les_faces_voisines(20),1))   + compute_viscosity_edges_sphere(
                                                    les_faces_voisines(20),les_faces_voisines(28),particle_id)*(valeurs_champ(
                                                                                                                  les_faces_voisines(20)) - valeurs_champ(les_faces_voisines(28)))  /
                                                  (xv(les_faces_voisines(20),1)  - xv(les_faces_voisines(28),1)))
                                          + 1./2.*( compute_viscosity_edges_sphere(les_faces_voisines(27),
                                                                                   les_faces_voisines(21),particle_id)*(valeurs_champ(les_faces_voisines(27)) -
                                                                                                                        valeurs_champ(les_faces_voisines(21))) / (xv(les_faces_voisines(27),1) -
                                                                                                                            xv(les_faces_voisines(21),1)) + compute_viscosity_edges_sphere(
                                                      les_faces_voisines(21),les_faces_voisines(29),particle_id)*(valeurs_champ(
                                                                                                                    les_faces_voisines(21)) - valeurs_champ(les_faces_voisines(29)))  /
                                                    (xv(les_faces_voisines(21),1)  - xv(les_faces_voisines(29),1))));

              gradU(elem,2,2) = ( (mu_p*mu_f/(mu_f-phase_indicator_function(elem_voisins(elem))*
                                              (mu_f-mu_p)))*(valeurs_champ(les_faces_voisines(20)) - valeurs_champ(
                                                               les_faces_voisines(21))) / (xv(les_faces_voisines(20),2) -
                                                                                           xv(les_faces_voisines(21),2)));
            }

          // trilinear interpolation of each component of the tensor
          for (int i=0; i<dimension; i++)
            {
              for (int j=0; j<dimension; j++)
                {
                  resu(fa7,i,j)=((1-zfact)*(
                                   (1-yfact)*(
                                     (1-xfact)*(gradU(0,i,j)) + xfact*(gradU(1,i,j))) +
                                   yfact*((1-xfact)*(gradU(2,i,j)) + xfact*(gradU(3,i,j)))
                                 ) +
                                 zfact*(
                                   (1-yfact)*((1-xfact)*(gradU(4,i,j)) + xfact*(gradU(5,i,j))) +
                                   yfact*((1-xfact)*(gradU(6,i,j)) + xfact*(gradU(7,i,j))))
                                );
                }
            }
        }
    }
  return 1;
}

int Post_Processing_Hydrodynamic_Forces::trilinear_interpolation_gradU_elem(
  const DoubleTab& valeurs_champ,
  DoubleTab& coord,
  DoubleTab& resu)
{
  Transport_Interfaces_FT_Disc& eq_transport = ptr_eq_transport_.valeur();
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  const Fluide_Diphasique& two_phase_fluid = eq_ns.fluide_diphasique();
  const Maillage_FT_Disc& mesh = eq_transport.maillage_interface();
  const Domaine_VDF& domain_vdf = ref_cast(Domaine_VDF, eq_ns.domaine_dis());
  const DoubleTab& xv = domain_vdf.xv();
  int nb_voisins=8;
  const int id_fluid_phase=two_phase_fluid.get_id_fluid_phase();
  const double mu_f =two_phase_fluid.fluide_phase(id_fluid_phase).
                     viscosite_dynamique().valeurs()(0, 0);

  for (int fa7=0; fa7<coord.dimension(0); fa7++)
    {
      if (!mesh.facette_virtuelle(fa7))
        {
          // find the 8 neighboring elements where the stress tensor will be computed
          DoubleVect coord_elem_interp(dimension);
          for (int dim=0; dim<dimension; dim++) coord_elem_interp(dim)=coord(fa7,dim);
          IntVect elem_voisins(nb_voisins);
          find_neighboring_elements(coord_elem_interp,elem_voisins);

          for (int i=0; i<nb_voisins; i++)
            {
              if (elem_voisins(i)<0)
                return 0;
            }

          // compute distance between neighboring meshes in all direction to compute
          // interpolation coefficients

          DoubleVect delta_i(dimension);
          delta_i(0) = fabs(domain_vdf.dist_elem(elem_voisins(0), elem_voisins(1), 0));
          delta_i(1) = fabs(domain_vdf.dist_elem(elem_voisins(0), elem_voisins(3), 1));
          delta_i(2) = fabs(domain_vdf.dist_elem(elem_voisins(0), elem_voisins(5), 2));
          DoubleVect coord_elem_0(dimension);
          for (int dim=0; dim<dimension; dim++) coord_elem_0(dim)=domain_vdf.xp(elem_voisins(0),dim);

          double xfact=fabs((coord_elem_interp(0)-coord_elem_0(0))/delta_i(0));
          double yfact=fabs((coord_elem_interp(1)-coord_elem_0(1))/delta_i(1));
          double zfact=fabs((coord_elem_interp(2)-coord_elem_0(2))/delta_i(2));

          // compute the velocity gradient tensor for each neighboring element
          DoubleTab gradU(nb_voisins, dimension, dimension);
          for (int elem=0; elem<nb_voisins; elem++)
            {

              // 														SCHEMA EN 2D
              //  					AVANT				 ---- ---- --- --- --- --- --- --- --- ---   AVANT
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //														 9	 10
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //  										 ---- --- --- -5- -3- -6- --- --- ---
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //	  	  	  	  	GAUCHE  	  	  	  	  	  	     1 x 2              		DROITE
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |
              //  									     ---- --- --- -7- -4- -8- --- --- --- 	      			y
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |	 	  			^
              //												         11  12									|
              // 										|	 |	 |	 |	 |	 |	 |	 |	 |	 |	 	  			|
              //  					ARRIERE					 ---- --- --- --- --- --- --- --- ---    ARRIERE    o----->x
              //																					 			 z
              // Evaluation du gradient sur les faces de normale x (choix arbitraire mais peut tuer une certaine symetrie)
              // Les termes + designent les faces juxtaposees selon +z
              // Les termes - designent les faces juxtaposees selon -z
              // Les termes 57 et 86 pour w designent les faces de normale z secantes aux faces 5,7 et 8,6
              // 		   --					   													   	   									--
              //		   |   u_2-u_1	 		   				  1	  u_9-u_11	   u_10-u_12		          1   u_1+ - u_1- 	 u_2+ - u_2-	|
              //		   | ----------							  -	(---------- + ----------- )	        	  -	( -----------  + ----------- )	|
              //		   |   x_2-x1			 				  2	  y_9-y_11	   y_10-y_12				  2   z_1+ - z_1-	 z_2+ - z_2- 	|
              //		   |												    																|
              //  grad U = | 1    v_5-v_6	 v_7 - v_8    			 	v_3-v_4								  1	  v_3+ -v_3-	v_4+ - v_4-     |
              // 	  	   | - ( --------- + --------- ) 			 	------- 						      - ( ----------- + ----------- )	|
              //		   | 2    x_5-x_6	 x_7 - x_8   			 	y_3-y_4						     	  2	  z_3+ - z_3-	z_4+ - z_4-		|
              //		   |																													|
              // 	  	   | 1	 w57+ - w86+	 w57- - w86-	  1  	w1112+ - w910+    w1112- -w910-  	      w_34+ - w_34-					|
              //		   | - ( ----------- +	 ----------- )	  - (	-------------- + -------------- )		  -------------	    			|
              //  		   | 2	 x57+ - x86+	 x57- - x86-  	  2	    y1112+ - y910+	  y1112--y910-	  		  z_34+ - z_34-					|
              //		   --					   		   									              									   --
              //
              // Pour chaque face, il faut donc 30 faces voisines : 1-12, 1+, 1-, 2+, 2-, 3+, 3-, 4+, 4-, 34+, 34-, 57+, 57-, 86+, 86-, 1112+, 1112-, 910+, 910-
              // Pour plus de facilites, on numerote :  1+:13, 1-:14, 2+:15, 2-:16, 3+:17, 3-:18, 4+:19, 4-:20, 34+:21, 34-:22, 57+:23, 57-:24, 86+:25, 86-:26, 1112+:27, 1112-:28, 910+:29, 910-:30

              // On identifie les faces voisines pour le calcul des composantes du tenseur
              int nb_faces_voisines=30;
              IntVect les_faces_voisines(nb_faces_voisines);
              int elem_=elem_voisins(elem);
              int elem_gauche = face_voisins_for_interp(elem_faces_for_interp(elem_,0),0);
              int elem_droite = face_voisins_for_interp(elem_faces_for_interp(elem_,0+dimension),1);
              int elem_haut=face_voisins_for_interp(elem_faces_for_interp(elem_,2+dimension), 1);
              int elem_bas=face_voisins_for_interp(elem_faces_for_interp(elem_,2), 0);
              int elem_avant=face_voisins_for_interp(elem_faces_for_interp(elem_,1+dimension),1);
              int elem_arriere=face_voisins_for_interp(elem_faces_for_interp(elem_,1),0);

              les_faces_voisines(0)=elem_faces_for_interp(elem_,0);
              les_faces_voisines(1)=elem_faces_for_interp(elem_,0+dimension);
              les_faces_voisines(2)=elem_faces_for_interp(elem_,1);
              les_faces_voisines(3)=elem_faces_for_interp(elem_,1+dimension);
              les_faces_voisines(4)=elem_faces_for_interp(elem_gauche,1+dimension);
              les_faces_voisines(5)=elem_faces_for_interp(elem_droite,1+dimension);
              les_faces_voisines(6)=elem_faces_for_interp(elem_gauche,1);
              les_faces_voisines(7)=elem_faces_for_interp(elem_droite,1);
              les_faces_voisines(8)=elem_faces_for_interp(elem_avant,0);
              les_faces_voisines(9)=elem_faces_for_interp(elem_avant,0+dimension);
              les_faces_voisines(10)=elem_faces_for_interp(elem_arriere,0);
              les_faces_voisines(11)=elem_faces_for_interp(elem_arriere,0+dimension);
              les_faces_voisines(12)=elem_faces_for_interp(elem_haut,0);
              les_faces_voisines(13)=elem_faces_for_interp(elem_bas,0);
              les_faces_voisines(14)=elem_faces_for_interp(elem_haut,0+dimension);
              les_faces_voisines(15)=elem_faces_for_interp(elem_bas,0+dimension);
              les_faces_voisines(16)=elem_faces_for_interp(elem_haut,1+dimension);
              les_faces_voisines(17)=elem_faces_for_interp(elem_bas,1+dimension);
              les_faces_voisines(18)=elem_faces_for_interp(elem_haut,1);
              les_faces_voisines(19)=elem_faces_for_interp(elem_bas,1);
              les_faces_voisines(20)=elem_faces_for_interp(elem_,2);
              les_faces_voisines(21)=elem_faces_for_interp(elem_,2+dimension);
              les_faces_voisines(22)=elem_faces_for_interp(elem_gauche,2);
              les_faces_voisines(23)=elem_faces_for_interp(elem_gauche,2+dimension);
              les_faces_voisines(24)=elem_faces_for_interp(elem_droite,2);
              les_faces_voisines(25)=elem_faces_for_interp(elem_droite,2+dimension);
              les_faces_voisines(26)=elem_faces_for_interp(elem_arriere,2);
              les_faces_voisines(27)=elem_faces_for_interp(elem_arriere,2+dimension);
              les_faces_voisines(28)=elem_faces_for_interp(elem_avant,2);
              les_faces_voisines(29)=elem_faces_for_interp(elem_avant,2+dimension);

              for (int i=0; i<nb_faces_voisines; i++)
                {
                  if (les_faces_voisines(i)<0)
                    return 0;
                }

              gradU(elem,0,0) = ( (valeurs_champ(les_faces_voisines(1))  - valeurs_champ(
                                     les_faces_voisines(0)))  / (xv(les_faces_voisines(1),0)  - xv(les_faces_voisines(0),0)));

              gradU(elem,0,1) = 1./2.*( (valeurs_champ(les_faces_voisines(8))  - valeurs_champ(
                                           les_faces_voisines(10))) / (xv(les_faces_voisines(8),1)  -
                                                                       xv(les_faces_voisines(10),1))   + (valeurs_champ(les_faces_voisines(9))  -
                                                                                                          valeurs_champ(les_faces_voisines(11)))  / (xv(les_faces_voisines(9),1)   -
                                                                                                              xv(les_faces_voisines(11),1)));

              gradU(elem,0,2) = 1./2.*( (valeurs_champ(les_faces_voisines(12)) - valeurs_champ(
                                           les_faces_voisines(13))) / (xv(les_faces_voisines(12),2) -
                                                                       xv(les_faces_voisines(13),2))   + (valeurs_champ(les_faces_voisines(14)) -
                                                                                                          valeurs_champ(les_faces_voisines(15)))  / (xv(les_faces_voisines(14),2)  -
                                                                                                              xv(les_faces_voisines(15),2)));

              gradU(elem,1,0) = 1./2.*( (valeurs_champ(les_faces_voisines(4)) - valeurs_champ(
                                           les_faces_voisines(5)))  / (xv(les_faces_voisines(4),0)  -
                                                                       xv(les_faces_voisines(5),0))    + (valeurs_champ(les_faces_voisines(6)) -
                                                                                                          valeurs_champ(les_faces_voisines(7)))   / (xv(les_faces_voisines(6),0) -
                                                                                                              xv(les_faces_voisines(7),0)));

              gradU(elem,1,1) = ( (valeurs_champ(les_faces_voisines(2))  - valeurs_champ(
                                     les_faces_voisines(3)))  / (xv(les_faces_voisines(2),1)  -
                                                                 xv(les_faces_voisines(3),1)));

              gradU(elem,1,2) = 1./2.*( (valeurs_champ(les_faces_voisines(16)) - valeurs_champ(
                                           les_faces_voisines(17))) / (xv(les_faces_voisines(16),2) -
                                                                       xv(les_faces_voisines(17),2))   + (valeurs_champ(les_faces_voisines(18)) -
                                                                                                          valeurs_champ(les_faces_voisines(19)))  / (xv(les_faces_voisines(18),2)  -
                                                                                                              xv(les_faces_voisines(19),2)));

              gradU(elem,2,0) = 1./2.*( (valeurs_champ(les_faces_voisines(22)) - valeurs_champ(
                                           les_faces_voisines(24))) / (xv(les_faces_voisines(22),0) -
                                                                       xv(les_faces_voisines(24),0))   + (valeurs_champ(les_faces_voisines(23)) -
                                                                                                          valeurs_champ(les_faces_voisines(25)))  / (xv(les_faces_voisines(23),0)  -
                                                                                                              xv(les_faces_voisines(25),0)));

              gradU(elem,2,1) = 1./2.*( (valeurs_champ(les_faces_voisines(26)) - valeurs_champ(
                                           les_faces_voisines(28))) / (xv(les_faces_voisines(26),1) -
                                                                       xv(les_faces_voisines(28),1))   + (valeurs_champ(les_faces_voisines(27)) -
                                                                                                          valeurs_champ(les_faces_voisines(29)))  / (xv(les_faces_voisines(27),1)  -
                                                                                                              xv(les_faces_voisines(29),1)));

              gradU(elem,2,2) = ( (valeurs_champ(les_faces_voisines(20)) - valeurs_champ(
                                     les_faces_voisines(21))) / (xv(les_faces_voisines(20),2) -
                                                                 xv(les_faces_voisines(21),2)));
            }

          // On fait une interpolation trilineaire de chaque composante du tenseur
          for (int i=0; i<dimension; i++)
            {
              for (int j=0; j<dimension; j++)
                {
                  resu(fa7,i,j)=mu_f*((1-zfact)*(
                                        (1-yfact)*(
                                          (1-xfact)*(gradU(0,i,j)) + xfact*(gradU(1,i,j))) +
                                        yfact*((1-xfact)*(gradU(2,i,j)) + xfact*(gradU(3,i,j)))
                                      ) +
                                      zfact*(
                                        (1-yfact)*((1-xfact)*(gradU(4,i,j)) + xfact*(gradU(5,i,j))) +
                                        yfact*((1-xfact)*(gradU(6,i,j)) + xfact*(gradU(7,i,j))))
                                     );
                }
            }
        }
    }
  return 1;
}

int Post_Processing_Hydrodynamic_Forces::trilinear_interpolation_face_sommets(
  const DoubleTab& valeurs_champ,
  DoubleTab& coord,
  DoubleTab& resu)
{
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  const Transport_Interfaces_FT_Disc& eq_transport = ptr_eq_transport_.valeur();
  const Domaine_VDF& domain_vdf = ref_cast(Domaine_VDF, eq_ns.domaine_dis());
  const Maillage_FT_Disc& mesh = eq_transport.maillage_interface();
  const ArrOfInt& elem = mesh.sommet_elem();
  const DoubleTab& pos = mesh.sommets();
  const int nb_pos_tot = pos.dimension(0);
  int nb_voisins=8;
  for (int som=0; som<nb_pos_tot; som++)
    {
      const int element = elem[som];
      if (element >= 0)   // real vertice ?
        {
          // find the neighboring faces that will be used for interpolation
          DoubleVect coord_elem_interp(dimension); // coordinate of the fa7 gravity center
          for (int dim=0; dim<dimension; dim++)
            coord_elem_interp(dim)=coord(som,dim);
          IntTab faces_voisines(dimension,nb_voisins);
          find_neighboring_faces_xyz(coord_elem_interp,faces_voisines);
          for (int dim=0; dim<dimension; dim++)
            {
              for (int i=0; i<nb_voisins; i++)
                {
                  if (faces_voisines(dim,i)<0)
                    return 0; // if one face could not be computed (access problem : joint zone, wall),
                  // we return 0 and the fa7 is not considered in the computation
                }
            }

          DoubleVect xfact(dimension);
          DoubleVect yfact(dimension);
          DoubleVect zfact(dimension);
          DoubleTab Delta_(dimension);
          Delta_(0)=fabs(domain_vdf.dist_face(faces_voisines(0,0),faces_voisines(0,1),0));
          Delta_(1)=fabs(domain_vdf.dist_face(faces_voisines(0,0),faces_voisines(0,2),1));
          Delta_(2)=fabs(domain_vdf.dist_face(faces_voisines(0,0),faces_voisines(0,4),2));

          for (int dim=0; dim<dimension; dim++)
            {
              xfact(dim)=fabs((coord(som,0)-domain_vdf.xv(faces_voisines(dim,0),0))/Delta_(0));
              yfact(dim)=fabs((coord(som,1)-domain_vdf.xv(faces_voisines(dim,0),1))/Delta_(1));
              zfact(dim)=fabs((coord(som,2)-domain_vdf.xv(faces_voisines(dim,0),2))/Delta_(2));
            }

          for (int dim=0; dim<dimension; dim++)
            {
              if (xfact(dim)>1) Cerr << "xfact > 1 pour le point " <<coord(som,0)<< " " << coord(som,1)<< " " << coord(som,2) <<finl;
              if (yfact(dim)>1) Cerr << "yfact > 1 pour le point " <<coord(som,0)<< " " << coord(som,1)<< " " << coord(som,2) <<finl;
              if (zfact(dim)>1) Cerr << "zfact > 1 pour le point " <<coord(som,0)<< " " << coord(som,1)<< " " << coord(som,2) <<finl;
            }
          for (int dim=0; dim<dimension; dim++)
            {
              resu(som,dim)=(1.-zfact(dim))*(
                              (1.-yfact(dim))*(
                                (1.-xfact(dim))*(valeurs_champ(faces_voisines(dim,0))) +
                                xfact(dim)*(valeurs_champ(faces_voisines(dim,1)))
                              ) +
                              yfact(dim)*(
                                (1.-xfact(dim))*(valeurs_champ(faces_voisines(dim,2))) +
                                xfact(dim)*(valeurs_champ(faces_voisines(dim,3))))
                            ) +
                            zfact(dim)*(
                              (1.-yfact(dim))*(
                                (1.-xfact(dim))*(valeurs_champ(faces_voisines(dim,4))) +
                                xfact(dim)*(valeurs_champ(faces_voisines(dim,5)))
                              ) +
                              yfact(dim)*(
                                (1.-xfact(dim))*(valeurs_champ(faces_voisines(dim,6))) +
                                xfact(dim)*(valeurs_champ(faces_voisines(dim,7))))
                            );
            }
        }
    }

  return 1;
}

double Post_Processing_Hydrodynamic_Forces::find_neighboring_elements(
  DoubleVect& coord_elem_interp,
  IntVect& neighboring_elements,
  const int sauv_list_P1,
  const int num_fa7)
{
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  Transport_Interfaces_FT_Disc& eq_transport = ptr_eq_transport_.valeur();
  const DoubleTab& phase_indicator_function = eq_transport.get_update_indicatrice().valeurs();
  const Domaine_VDF& domain_vdf = ref_cast(Domaine_VDF, eq_ns.domaine_dis());
  const Domaine& domain = domain_vdf.domaine();

  DoubleVect coord_elem_eulerian(dimension);
  int elem_eulerien=domain.chercher_elements(coord_elem_interp(0),
                                             coord_elem_interp(1),
                                             coord_elem_interp(2));
  if (sauv_list_P1)
    list_elem_P1_(num_fa7,0)=elem_eulerien;
  if (elem_eulerien<0)
    {
      neighboring_elements=-1;
      return -1;
    }

  for (int dim=0; dim<dimension; dim++)
    coord_elem_eulerian(dim)=domain_vdf.xp(elem_eulerien,dim);
  IntVect direction_interp(dimension); // Pour chaque direction, on regarde de quel cote de l'element le point se trouve. 0 : a gauche (pos_i_point<=pos_i_elem), 1 : a droite(pos_i_point>pos_i_elem)
  IntVect faces_elem_interp(2*dimension);
  for (int dim=0; dim<dimension; dim++)
    {
      faces_elem_interp(dim)=elem_faces_for_interp(elem_eulerien,dim);
      faces_elem_interp(dimension+dim)=elem_faces_for_interp(elem_eulerien,dimension+dim);

      if (coord_elem_interp(dim)<=coord_elem_eulerian(dim))
        direction_interp(dim)=0;
      else
        direction_interp(dim)=1;
    }

  // The 8 elements form a large cube with 2 cubes in each direction (1 cube = 1 elem).
  // Let's consider a cube with side 1. Each vertice represent a neighboring element.
  // We fill the neighboring list beginning with coordinates
  // z=0, y=0 then z=0, y=1, then z=1, y=0 then z=1, y=1
  // Thus, we have (0,0,0) ; (1,0,0) ; (0,1,0) ; (1,1,0) ; (0,0,+-1) ; (1,0,+-1) ; (0,1,+-1) ; (1,1,+-1)
  if (direction_interp(0)==1)
    {
      if (direction_interp(1)==1)
        {
          neighboring_elements(0)=elem_eulerien;
          neighboring_elements(1)=face_voisins_for_interp(faces_elem_interp(0+dimension),1);
          neighboring_elements(2)=face_voisins_for_interp(faces_elem_interp(1+dimension),1);
          neighboring_elements(3)=face_voisins_for_interp(elem_faces_for_interp(
                                                            neighboring_elements(1),1+dimension),1);
        }
      else
        {
          neighboring_elements(0)=face_voisins_for_interp(faces_elem_interp(1),0);
          neighboring_elements(1)=face_voisins_for_interp(elem_faces_for_interp(
                                                            neighboring_elements(0),0+dimension),1);
          neighboring_elements(2)=elem_eulerien;
          neighboring_elements(3)=face_voisins_for_interp(faces_elem_interp(0+dimension),1);
        }
    }
  else
    {
      if (direction_interp(1)==1)
        {
          neighboring_elements(0)=face_voisins_for_interp(faces_elem_interp(0),0);
          neighboring_elements(1)=elem_eulerien;
          neighboring_elements(2)=face_voisins_for_interp(elem_faces_for_interp(
                                                            neighboring_elements(0),1+dimension),1);
          neighboring_elements(3)=face_voisins_for_interp(faces_elem_interp(1+dimension),1);
        }
      else
        {
          neighboring_elements(1)=face_voisins_for_interp(faces_elem_interp(1),0);
          neighboring_elements(2)=face_voisins_for_interp(faces_elem_interp(0),0);
          neighboring_elements(3)=elem_eulerien;
          neighboring_elements(0)=face_voisins_for_interp(elem_faces_for_interp(
                                                            neighboring_elements(1),0),0);
        }
    }

  if (direction_interp(2)==1)
    {
      neighboring_elements(4)=face_voisins_for_interp(elem_faces_for_interp(
                                                        neighboring_elements(0),2+dimension),1);
      neighboring_elements(5)=face_voisins_for_interp(elem_faces_for_interp(
                                                        neighboring_elements(1),2+dimension),1);
      neighboring_elements(6)=face_voisins_for_interp(elem_faces_for_interp(
                                                        neighboring_elements(2),2+dimension),1);
      neighboring_elements(7)=face_voisins_for_interp(elem_faces_for_interp(
                                                        neighboring_elements(3),2+dimension),1);
    }
  else
    {
      neighboring_elements(4)=face_voisins_for_interp(elem_faces_for_interp(
                                                        neighboring_elements(0),2),0);
      neighboring_elements(5)=face_voisins_for_interp(elem_faces_for_interp(
                                                        neighboring_elements(1),2),0);
      neighboring_elements(6)=face_voisins_for_interp(elem_faces_for_interp(
                                                        neighboring_elements(2),2),0);
      neighboring_elements(7)=face_voisins_for_interp(elem_faces_for_interp(
                                                        neighboring_elements(3),2),0);

      for (int i=0; i<4; i++)
        {
          int tmp=neighboring_elements(i);
          neighboring_elements(i)=neighboring_elements(i+4);
          neighboring_elements(i+4)=tmp;
        }
    }
  return (phase_indicator_function(elem_eulerien));
}

void Post_Processing_Hydrodynamic_Forces::find_neighboring_faces(
  DoubleVect& coord_elem_interp,
  IntVect& neighboring_faces,
  int orientation) const
{
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  const Domaine_VDF& domain_vdf = ref_cast(Domaine_VDF, eq_ns.domaine_dis());

  DoubleVect coord_elem_eulerien(dimension);

  int elem_eulerien=domain_vdf.domaine().chercher_elements(coord_elem_interp(0),
                                                           coord_elem_interp(1),
                                                           coord_elem_interp(2));
  if (elem_eulerien<0)
    {
      neighboring_faces=-1;
      return;
    }
  for (int dim=0; dim<dimension; dim++)
    coord_elem_eulerien(dim)=domain_vdf.xp(elem_eulerien,dim);
  // In each direction, we look at which side of the element the point is.
  // 0 : left (pos_i_point<=pos_i_elem), 1 : right (pos_i_point>pos_i_elem)
  IntVect direction_interp(dimension);
  IntVect faces_elem_interp(2*dimension);

  for (int dim=0; dim<dimension; dim++)
    {
      faces_elem_interp(dim)=elem_faces_for_interp(elem_eulerien,dim);
      faces_elem_interp(dimension+dim)=elem_faces_for_interp(elem_eulerien,dimension+dim);

      if (coord_elem_interp(dim)<=coord_elem_eulerien(dim))
        direction_interp(dim)=0;
      else
        direction_interp(dim)=1;
    }

  if (orientation==0)
    {
      if (direction_interp(1)==1)
        {
          neighboring_faces(0)=faces_elem_interp(orientation);
          neighboring_faces(1)=faces_elem_interp(orientation+dimension);
          neighboring_faces(2)=elem_faces_for_interp(face_voisins_for_interp(
                                                       faces_elem_interp(1+dimension),1),orientation);
          neighboring_faces(3)=elem_faces_for_interp(face_voisins_for_interp(
                                                       faces_elem_interp(1+dimension),1),orientation+dimension);
        }
      else
        {
          neighboring_faces(0)=elem_faces_for_interp(face_voisins_for_interp(
                                                       faces_elem_interp(1),0),orientation);
          neighboring_faces(1)=elem_faces_for_interp(face_voisins_for_interp(
                                                       faces_elem_interp(1),0),orientation+dimension);
          neighboring_faces(2)=faces_elem_interp(orientation);
          neighboring_faces(3)=faces_elem_interp(orientation+dimension);
        }
    }
  if (orientation==1)
    {
      if (direction_interp(0)==1)
        {
          neighboring_faces(0)=faces_elem_interp(orientation);
          neighboring_faces(1)=elem_faces_for_interp(face_voisins_for_interp(
                                                       faces_elem_interp(0+dimension),1),orientation);
          neighboring_faces(2)=faces_elem_interp(orientation+dimension);
          neighboring_faces(3)=elem_faces_for_interp(face_voisins_for_interp(
                                                       faces_elem_interp(0+dimension),1),orientation+dimension);
        }
      else
        {
          neighboring_faces(0)=elem_faces_for_interp(face_voisins_for_interp(
                                                       faces_elem_interp(0),0),orientation);
          neighboring_faces(1)=faces_elem_interp(orientation);
          neighboring_faces(2)=elem_faces_for_interp(face_voisins_for_interp(
                                                       faces_elem_interp(0),0),orientation+dimension);
          neighboring_faces(3)=faces_elem_interp(orientation+dimension);
        }
    }
  if (direction_interp(2)==1 && orientation!=2)
    {
      neighboring_faces(4)=elem_faces_for_interp(face_voisins_for_interp(
                                                   elem_faces_for_interp(face_voisins_for_interp(neighboring_faces(0),1),
                                                                         2+dimension),1),orientation);
      neighboring_faces(5)=elem_faces_for_interp(face_voisins_for_interp(
                                                   elem_faces_for_interp(face_voisins_for_interp(neighboring_faces(1),1),
                                                                         2+dimension),1),orientation);
      neighboring_faces(6)=elem_faces_for_interp(face_voisins_for_interp(
                                                   elem_faces_for_interp(face_voisins_for_interp(neighboring_faces(2),1),
                                                                         2+dimension),1),orientation);
      neighboring_faces(7)=elem_faces_for_interp(face_voisins_for_interp(
                                                   elem_faces_for_interp(face_voisins_for_interp(neighboring_faces(3),1),
                                                                         2+dimension),1),orientation);
    }
  else if (direction_interp(2)==0 && orientation!=2)
    {
      neighboring_faces(4)=elem_faces_for_interp(face_voisins_for_interp(
                                                   elem_faces_for_interp(face_voisins_for_interp(neighboring_faces(0),1),
                                                                         2),0),orientation);
      neighboring_faces(5)=elem_faces_for_interp(face_voisins_for_interp(
                                                   elem_faces_for_interp(face_voisins_for_interp(neighboring_faces(1),1),
                                                                         2),0),orientation);
      neighboring_faces(6)=elem_faces_for_interp(face_voisins_for_interp(
                                                   elem_faces_for_interp(face_voisins_for_interp(neighboring_faces(2),1),
                                                                         2),0),orientation);
      neighboring_faces(7)=elem_faces_for_interp(face_voisins_for_interp(
                                                   elem_faces_for_interp(face_voisins_for_interp(neighboring_faces(3),1)
                                                                         ,2),0),orientation);

      for (int i=0; i<4; i++)
        {
          int tmp=neighboring_faces(i);
          neighboring_faces(i)=neighboring_faces(i+4);
          neighboring_faces(i+4)=tmp;
        }
    }

  if (orientation==2)
    {
      if (direction_interp(0)==1)
        {
          neighboring_faces(0)=faces_elem_interp(orientation);
          neighboring_faces(1)=elem_faces_for_interp(face_voisins_for_interp(
                                                       faces_elem_interp(0+dimension),1),orientation);
          neighboring_faces(4)=faces_elem_interp(orientation+dimension);
          neighboring_faces(5)=elem_faces_for_interp(face_voisins_for_interp(
                                                       faces_elem_interp(0+dimension),1),orientation+dimension);
        }
      else
        {
          neighboring_faces(0)=elem_faces_for_interp(face_voisins_for_interp(
                                                       faces_elem_interp(0),0),orientation);
          neighboring_faces(1)=faces_elem_interp(orientation);
          neighboring_faces(4)=elem_faces_for_interp(face_voisins_for_interp(
                                                       faces_elem_interp(0),0),orientation+dimension);
          neighboring_faces(5)=faces_elem_interp(orientation+dimension);

        }

      if (direction_interp(1)==1)
        {
          neighboring_faces(2)=elem_faces_for_interp(face_voisins_for_interp(
                                                       elem_faces_for_interp(face_voisins_for_interp(neighboring_faces(0),1),
                                                                             1+dimension),1),orientation);
          neighboring_faces(3)=elem_faces_for_interp(face_voisins_for_interp(
                                                       elem_faces_for_interp(face_voisins_for_interp(neighboring_faces(1),1),
                                                                             1+dimension),1),orientation);
          neighboring_faces(6)=elem_faces_for_interp(face_voisins_for_interp(
                                                       elem_faces_for_interp(face_voisins_for_interp(neighboring_faces(4),1),
                                                                             1+dimension),1),orientation);
          neighboring_faces(7)=elem_faces_for_interp(face_voisins_for_interp(
                                                       elem_faces_for_interp(face_voisins_for_interp(neighboring_faces(5),1),
                                                                             1+dimension),1),orientation);
        }
      else
        {
          neighboring_faces(2)=elem_faces_for_interp(face_voisins_for_interp(
                                                       elem_faces_for_interp(face_voisins_for_interp(neighboring_faces(0),1),
                                                                             1),0),orientation);
          neighboring_faces(3)=elem_faces_for_interp(face_voisins_for_interp(
                                                       elem_faces_for_interp(face_voisins_for_interp(neighboring_faces(1),1),
                                                                             1),0),orientation);
          neighboring_faces(6)=elem_faces_for_interp(face_voisins_for_interp(
                                                       elem_faces_for_interp(face_voisins_for_interp(neighboring_faces(4),1),
                                                                             1),0),orientation);
          neighboring_faces(7)=elem_faces_for_interp(face_voisins_for_interp(
                                                       elem_faces_for_interp(face_voisins_for_interp(neighboring_faces(5),1),
                                                                             1),0),orientation);

          int tmp0=neighboring_faces(0);
          int tmp1=neighboring_faces(1);
          int tmp4=neighboring_faces(4);
          int tmp5=neighboring_faces(5);

          neighboring_faces(0)=neighboring_faces(2);
          neighboring_faces(1)=neighboring_faces(3);
          neighboring_faces(4)=neighboring_faces(6);
          neighboring_faces(5)=neighboring_faces(7);

          neighboring_faces(2)=tmp0;
          neighboring_faces(3)=tmp1;
          neighboring_faces(6)=tmp4;
          neighboring_faces(7)=tmp5;
        }
    }
}

void Post_Processing_Hydrodynamic_Forces::find_neighboring_faces_xyz(
  DoubleVect& coord_elem_interp,
  IntTab& neighboring_faces) const
{
  int nb_neighboring_faces=8;
  IntVect neighboring_faces_x(nb_neighboring_faces);
  IntVect neighboring_faces_y(nb_neighboring_faces);
  IntVect neighboring_faces_z(nb_neighboring_faces);

  find_neighboring_faces(coord_elem_interp,neighboring_faces_x,0);
  find_neighboring_faces(coord_elem_interp,neighboring_faces_y,1);
  find_neighboring_faces(coord_elem_interp,neighboring_faces_z,2);

  for (int i=0; i<8; i++)
    {
      neighboring_faces(0,i)=neighboring_faces_x(i);
      neighboring_faces(1,i)=neighboring_faces_y(i);
      neighboring_faces(2,i)=neighboring_faces_z(i);
    }
}

// see Vincent et al., 2014
double Post_Processing_Hydrodynamic_Forces::compute_viscosity_edges_sphere(int face1, int face2, int particle_id)
{
  Transport_Interfaces_FT_Disc& eq_transport = ptr_eq_transport_.valeur();
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  const Fluide_Diphasique& two_phase_fluid = eq_ns.fluide_diphasique();

  const Domaine_VDF& domain_vdf = ref_cast(Domaine_VDF, eq_ns.domaine_dis());
  const DoubleTab& xp = domain_vdf.xp();
  const DoubleTab& particles_position=eq_transport.get_particles_position();

  const int id_fluid_phase=two_phase_fluid.get_id_fluid_phase();
  const Solid_Particle_base& solid_particle=ref_cast(Solid_Particle_base,two_phase_fluid.fluide_phase(id_fluid_phase));
  const double& effective_radius=solid_particle.get_equivalent_radius();
  const double mu_f =two_phase_fluid.fluide_phase(id_fluid_phase).
                     viscosite_dynamique().valeurs()(0, 0);
  const double mu_p = two_phase_fluid.fluide_phase(1-id_fluid_phase).
                      viscosite_dynamique().valeurs()(0, 0);

  int elem1 = face_voisins_for_interp(face1,0);
  int elem2 = face_voisins_for_interp(face1,1);
  int elem3 = face_voisins_for_interp(face2,0);
  int elem4 = face_voisins_for_interp(face2,1);

  double indicatrice_arete=0;
  int pvx=10;
  int pvy=10;
  int pvz=10;
  double x_min,y_min,z_min,x_max;
  double x_cg=particles_position(particle_id,0);
  double y_cg=particles_position(particle_id,1);
  double z_cg=particles_position(particle_id,2);

  x_min=std::min(xp(elem1,0),xp(elem2,0));
  x_min=std::min(x_min,xp(elem3,0));
  x_min=std::min(x_min,xp(elem4,0));
  x_max=std::max(xp(elem1,0),xp(elem2,0));
  x_max=std::max(x_min,xp(elem3,0));
  x_max=std::max(x_min,xp(elem4,0));
  y_min=std::min(xp(elem1,1),xp(elem2,1));
  y_min=std::min(y_min,xp(elem3,1));
  y_min=std::min(y_min,xp(elem4,1));
  z_min=std::min(xp(elem1,2),xp(elem2,2));
  z_min=std::min(z_min,xp(elem3,2));
  z_min=std::min(z_min,xp(elem4,2));

  double delta_x=x_max-x_min;

  for (int i=0; i<pvx; i++)
    {
      for (int j=0; j<pvy; j++)
        {
          for (int k=0; k<pvz; k++)
            {
              double x=x_min+i*delta_x/pvx;
              double y=y_min+i*delta_x/pvy;
              double z=z_min+i*delta_x/pvz;
              if (sqrt(pow(x-x_cg,2)+pow(y-y_cg,2)+pow(z-z_cg,2))>effective_radius)
                indicatrice_arete+=1;
            }
        }
    }
  indicatrice_arete/=(pvx*pvy*pvz);

  return(mu_p*mu_f/(mu_f-indicatrice_arete*(mu_f-mu_p)));
}

void Post_Processing_Hydrodynamic_Forces::resize_sigma(int nb_fa7)
{
  sigma_xx_fa7_.resize(nb_fa7);
  sigma_xy_fa7_.resize(nb_fa7);
  sigma_xz_fa7_.resize(nb_fa7);
  sigma_yx_fa7_.resize(nb_fa7);
  sigma_yy_fa7_.resize(nb_fa7);
  sigma_yz_fa7_.resize(nb_fa7);
  sigma_zx_fa7_.resize(nb_fa7);
  sigma_zy_fa7_.resize(nb_fa7);
  sigma_zz_fa7_.resize(nb_fa7);
}

void Post_Processing_Hydrodynamic_Forces::resize_gradU_P1(int nb_fa7)
{
  dUdx_P1_.resize(nb_fa7);
  dUdy_P1_.resize(nb_fa7);
  dUdz_P1_.resize(nb_fa7);
  dVdx_P1_.resize(nb_fa7);
  dVdy_P1_.resize(nb_fa7);
  dVdz_P1_.resize(nb_fa7);
  dWdx_P1_.resize(nb_fa7);
  dWdy_P1_.resize(nb_fa7);
  dWdz_P1_.resize(nb_fa7);
}

void Post_Processing_Hydrodynamic_Forces::resize_gradU_P2(int nb_fa7)
{
  dUdx_P2_.resize(nb_fa7);
  dUdy_P2_.resize(nb_fa7);
  dUdz_P2_.resize(nb_fa7);
  dVdx_P2_.resize(nb_fa7);
  dVdy_P2_.resize(nb_fa7);
  dVdz_P2_.resize(nb_fa7);
  dWdx_P2_.resize(nb_fa7);
  dWdy_P2_.resize(nb_fa7);
  dWdz_P2_.resize(nb_fa7);
}

void Post_Processing_Hydrodynamic_Forces::fill_gradU_P1(int fa7, DoubleTab gradU_P1)
{
  dUdx_P1_(fa7)=gradU_P1(fa7,0,0);
  dUdy_P1_(fa7)=gradU_P1(fa7,0,1);
  dUdz_P1_(fa7)=gradU_P1(fa7,0,2);
  dVdx_P1_(fa7)=gradU_P1(fa7,1,0);
  dVdy_P1_(fa7)=gradU_P1(fa7,1,1);
  dVdz_P1_(fa7)=gradU_P1(fa7,1,2);
  dWdx_P1_(fa7)=gradU_P1(fa7,2,0);
  dWdy_P1_(fa7)=gradU_P1(fa7,2,1);
  dWdz_P1_(fa7)=gradU_P1(fa7,2,2);
}

void Post_Processing_Hydrodynamic_Forces::fill_gradU_P2(int fa7, DoubleTab gradU_P2)
{
  dUdx_P2_(fa7)=gradU_P2(fa7,0,0);
  dUdy_P2_(fa7)=gradU_P2(fa7,0,1);
  dUdz_P2_(fa7)=gradU_P2(fa7,0,2);
  dVdx_P2_(fa7)=gradU_P2(fa7,1,0);
  dVdy_P2_(fa7)=gradU_P2(fa7,1,1);
  dVdz_P2_(fa7)=gradU_P2(fa7,1,2);
  dWdx_P2_(fa7)=gradU_P2(fa7,2,0);
  dWdy_P2_(fa7)=gradU_P2(fa7,2,1);
  dWdz_P2_(fa7)=gradU_P2(fa7,2,2);
}

void Post_Processing_Hydrodynamic_Forces::fill_sigma(int fa7, Matrice_Dense stress_tensor)
{
  sigma_xx_fa7_(fa7)=stress_tensor(0,0);
  sigma_xy_fa7_(fa7)=stress_tensor(0,1);
  sigma_xz_fa7_(fa7)=stress_tensor(0,2);
  sigma_yx_fa7_(fa7)=stress_tensor(1,0);
  sigma_yy_fa7_(fa7)=stress_tensor(1,1);
  sigma_yz_fa7_(fa7)=stress_tensor(1,2);
  sigma_zx_fa7_(fa7)=stress_tensor(2,0);
  sigma_zy_fa7_(fa7)=stress_tensor(2,1);
  sigma_zz_fa7_(fa7)=stress_tensor(2,2);
}

void Post_Processing_Hydrodynamic_Forces::resize_and_init_tables(int nb_particles_tot)
{
  total_pressure_force_.resize(nb_particles_tot,dimension);
  total_friction_force_.resize(nb_particles_tot,dimension);
  proportion_fa7_ok_UP2_.resize(nb_particles_tot,dimension);
  prop_P2_fluid_compo_.resize(nb_particles_tot);
  U_P2_moy_.resize(nb_particles_tot, dimension);
  total_surface_interf_.resize(nb_particles_tot);
  Nb_fa7_tot_par_compo_.resize(nb_particles_tot);
  total_heat_transfer_.resize(nb_particles_tot);

  total_pressure_force_=0;
  total_friction_force_=0;
  proportion_fa7_ok_UP2_=0;
  prop_P2_fluid_compo_=0;
  U_P2_moy_=0;
  total_surface_interf_=0;
  Nb_fa7_tot_par_compo_=0;
}

void Post_Processing_Hydrodynamic_Forces::resize_data_fa7(int nb_fa7)
{
  phase_indicator_function_P2_.resize(nb_fa7);

  if (is_post_process_pressure_fa7_)
    {
      pressure_fa7_.resize(nb_fa7);
      pressure_fa7_=3e15;
    }

  if (is_post_process_friction_force_fa7_)
    {
      friction_force_fa7_.resize(nb_fa7,dimension);
      friction_force_fa7_=3e15;
    }

  if (is_post_process_pressure_force_fa7_)
    {
      pressure_force_fa7_.resize(nb_fa7,dimension);
      pressure_force_fa7_=3e15;
    }

  if (is_post_process_stress_tensor_fa7_)
    {
      resize_sigma(nb_fa7);
      resize_gradU_P1(nb_fa7);
      resize_gradU_P2(nb_fa7);
    }
  U_P1_.resize(nb_fa7,dimension);
  U_P2_.resize(nb_fa7,dimension);
}

void Post_Processing_Hydrodynamic_Forces::resize_coord_neighbor_fluid_fa7(int nb_fa7)
{
  coord_neighbor_fluid_fa7_pressure_1_.resize(nb_fa7,dimension);
  coord_neighbor_fluid_fa7_pressure_2_.resize(nb_fa7,dimension);
  coord_neighbor_fluid_fa7_gradU_1_.resize(nb_fa7,dimension);
  coord_neighbor_fluid_fa7_gradU_2_.resize(nb_fa7,dimension);
}

void Post_Processing_Hydrodynamic_Forces::compute_friction_force_complet_tensor(int nb_fa7,
                                                                                const Maillage_FT_Disc& mesh,
                                                                                const IntVect& compo_connexes_fa7,
                                                                                const ArrOfDouble& fa7_surface,
                                                                                const DoubleTab& tab_fa7_normal)
{
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  DoubleTab grad_U_P1(nb_fa7, dimension, dimension);
  DoubleTab grad_U_P2(nb_fa7, dimension, dimension);
  grad_U_P1=-1e15;
  grad_U_P2=-1e30;
  int interp_gradU_P1_ok=0;
  int interp_gradU_P2_ok=0;
  if (location_stress_tensor_== Location_stress_tensor::FACES_NORMALE_X)
    {
      interp_gradU_P1_ok=trilinear_interpolation_gradU_face(eq_ns.la_vitesse->valeurs(),
                                                            coord_neighbor_fluid_fa7_gradU_1_, grad_U_P1);
      interp_gradU_P2_ok=trilinear_interpolation_gradU_face(eq_ns.la_vitesse->valeurs(),
                                                            coord_neighbor_fluid_fa7_gradU_2_, grad_U_P2);
    }
  else if (location_stress_tensor_== Location_stress_tensor::ELEMENTS)
    {
      interp_gradU_P1_ok=trilinear_interpolation_gradU_elem(eq_ns.la_vitesse->valeurs(),
                                                            coord_neighbor_fluid_fa7_gradU_1_,
                                                            grad_U_P1);
      interp_gradU_P2_ok=trilinear_interpolation_gradU_elem(eq_ns.la_vitesse->valeurs(),
                                                            coord_neighbor_fluid_fa7_gradU_2_,
                                                            grad_U_P2);
    }
  if ( interp_gradU_P1_ok==1 && interp_gradU_P2_ok==1 )
    {
      DoubleTab grad_U_extrapole(nb_fa7, dimension, dimension);
      grad_U_extrapole=1e20;

      for (int fa7=0; fa7<nb_fa7; fa7++)
        {
          fill_gradU_P1(fa7, grad_U_P1);
          fill_gradU_P2(fa7, grad_U_P2);

          for (int i=0; i<dimension; i++)
            {
              for (int j=0; j<dimension; j++)
                {
                  grad_U_extrapole(fa7,i,j) = grad_U_P2(fa7,i,j) -
                                              interpolation_distance_gradU_P2_*
                                              (grad_U_P2(fa7,i,j)-grad_U_P1(fa7,i,j))/
                                              (interpolation_distance_gradU_P2_-
                                               interpolation_distance_gradU_P1_);
                }
            }
        }
      for (int fa7=0; fa7<nb_fa7; fa7++)
        {
          int compo=compo_connexes_fa7(fa7);
          Matrice_Dense stress_tensor(dimension,dimension);
          for (int i=0; i<dimension; i++)
            {
              for (int j=0; j<dimension; j++)
                {
                  stress_tensor(i,j) = (grad_U_extrapole(fa7,i,j) +
                                        grad_U_extrapole(fa7,j,i));
                }
            }

          if (is_post_process_stress_tensor_fa7_)
            fill_sigma(fa7,stress_tensor);

          DoubleTab la_normale_fa7_x_surface(dimension);
          for (int dim=0; dim<dimension; dim++) la_normale_fa7_x_surface(dim) =
              fa7_surface(fa7)*tab_fa7_normal(fa7,dim);
          DoubleVect friction_force_fa7=stress_tensor*la_normale_fa7_x_surface;

          if (!mesh.facette_virtuelle(fa7))
            {
              if (is_post_process_friction_force_fa7_)
                {
                  for (int dim=0; dim<dimension; dim++)
                    friction_force_fa7_(fa7,dim)=friction_force_fa7(dim);
                }
              for (int dim=0; dim<dimension; dim++)
                total_friction_force_(compo,dim)+=friction_force_fa7(dim);
            }
        }
    }
}

void Post_Processing_Hydrodynamic_Forces::compute_friction_force_projected_tensor(int nb_fa7,
                                                                                  const Maillage_FT_Disc& mesh,
                                                                                  const IntVect& compo_connexes_fa7,
                                                                                  const ArrOfDouble& fa7_surface,
                                                                                  const DoubleTab& tab_fa7_normal)
{
  Transport_Interfaces_FT_Disc& eq_transport = ptr_eq_transport_.valeur();
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  const Domaine_VDF& domain_vdf = ref_cast(Domaine_VDF, eq_ns.domaine_dis());
  const Fluide_Diphasique& two_phase_fluid = eq_ns.fluide_diphasique();
  const Domaine& domain = domain_vdf.domaine();
  const int id_fluid_phase=two_phase_fluid.get_id_fluid_phase();
  const double mu_f =two_phase_fluid.fluide_phase(id_fluid_phase).
                     viscosite_dynamique().valeurs()(0, 0);
  const DoubleTab& gravity_center_fa7=mesh.get_gravity_center_fa7();

  int interp_U_P1_ok=0;
  int interp_U_P2_ok=0;
  U_P1_.resize(nb_fa7, dimension);
  U_P2_.resize(nb_fa7, dimension);

  DoubleTab U_P1_spherique(nb_fa7, dimension);
  DoubleTab U_P2_spherique(nb_fa7, dimension);
  DoubleTab U_cg_spherique(nb_fa7, dimension);
  DoubleTab Urr(nb_fa7);
  DoubleTab Uthetar(nb_fa7);
  DoubleTab Uphir(nb_fa7);

  U_P1_=-1e15;
  U_P2_=-1e30;
  U_P1_spherique=-1e15;
  U_P2_spherique=-1e30;
  U_cg_spherique=-1e20;
  Urr=1e8;
  Uthetar=1e12;
  Uphir=1e15;

  double theta=0;
  double phi=0;
  double distance_au_cg=0;
  const DoubleTab& positions_compo=eq_transport.get_particles_position();
  const DoubleTab& vitesses_compo = eq_transport.get_particles_velocity();

  if (location_stress_tensor_== Location_stress_tensor::ELEMENTS)
    {
      // 1. Trilinear interpolation in cartesian coordinates
      interp_U_P1_ok=trilinear_interpolation_face(eq_ns.la_vitesse->valeurs(),
                                                  coord_neighbor_fluid_fa7_gradU_1_,
                                                  U_P1_);
      interp_U_P2_ok=trilinear_interpolation_face(eq_ns.la_vitesse->valeurs(),
                                                  coord_neighbor_fluid_fa7_gradU_2_,
                                                  U_P2_);
    }
  if ( interp_U_P1_ok && interp_U_P2_ok )
    {
      // 2. Change of reference frame to spherical coordinates
      for (int fa7=0; fa7<nb_fa7; fa7++)
        {
          int compo=compo_connexes_fa7(fa7);
          if (!mesh.facette_virtuelle(fa7))
            {
              Nb_fa7_tot_par_compo_(compo)++;

              int cont=0; // if the velocity in P2 could not be computed
              // we do not compute the friction force for the fa7
              for (int dim=0; dim<dimension; dim++)
                {
                  if (U_P2_(fa7,dim)>-1e10)
                    {
                      U_P2_moy_(compo,dim)+= U_P2_(fa7,dim); // compute the mean velocity in P2
                      proportion_fa7_ok_UP2_(compo,dim)+=1;
                    }
                  else
                    cont=1;
                }
              if (cont)
                continue;

              DoubleVect distance_cg_vect(dimension);
              for (int i=0; i<dimension; i++) distance_cg_vect(i)=
                  coord_neighbor_fluid_fa7_gradU_1_(fa7,i)-positions_compo(compo,i);

              distance_au_cg=sqrt(local_carre_norme_vect(distance_cg_vect));

              if (fabs((coord_neighbor_fluid_fa7_gradU_1_(fa7,2)-positions_compo(compo,2))/
                       distance_au_cg)<=1)
                {
                  theta=acos((coord_neighbor_fluid_fa7_gradU_1_(fa7,2)-
                              positions_compo(compo,2))/distance_au_cg);
                }
              else if ((coord_neighbor_fluid_fa7_gradU_1_(fa7,2)-positions_compo(compo,2))/
                       distance_au_cg>1)
                {
                  theta=0;
                }
              else if ((coord_neighbor_fluid_fa7_gradU_1_(fa7,2)-positions_compo(compo,2))/
                       distance_au_cg<-1)
                {
                  theta=M_PI;
                }

              if ((coord_neighbor_fluid_fa7_gradU_1_(fa7,0)-positions_compo(compo,0))>0 &&
                  (coord_neighbor_fluid_fa7_gradU_1_(fa7,1)-positions_compo(compo,1))>=0)
                {
                  phi=atan((coord_neighbor_fluid_fa7_gradU_1_(fa7,1)-positions_compo(compo,1))/
                           (coord_neighbor_fluid_fa7_gradU_1_(fa7,0)-positions_compo(compo,0)));
                }
              else if ((coord_neighbor_fluid_fa7_gradU_1_(fa7,0)-positions_compo(compo,0))>0
                       && (coord_neighbor_fluid_fa7_gradU_1_(fa7,1)-positions_compo(compo,1))<0)
                {
                  phi=atan((coord_neighbor_fluid_fa7_gradU_1_(fa7,1)-positions_compo(compo,1))/
                           (coord_neighbor_fluid_fa7_gradU_1_(fa7,0)-positions_compo(compo,0)))
                      +2*M_PI;
                }
              else if ((coord_neighbor_fluid_fa7_gradU_1_(fa7,0)-positions_compo(compo,0))<0)
                {
                  phi=atan((coord_neighbor_fluid_fa7_gradU_1_(fa7,1)-positions_compo(compo,1))/
                           (coord_neighbor_fluid_fa7_gradU_1_(fa7,0)-positions_compo(compo,0)))
                      +M_PI;
                }
              else if ((coord_neighbor_fluid_fa7_gradU_1_(fa7,0)-positions_compo(compo,0))==0
                       && (coord_neighbor_fluid_fa7_gradU_1_(fa7,1)-positions_compo(compo,1))>0)
                {
                  phi=M_PI/2.;
                }
              else if ((coord_neighbor_fluid_fa7_gradU_1_(fa7,0)-positions_compo(compo,0))==0
                       && (coord_neighbor_fluid_fa7_gradU_1_(fa7,1)-positions_compo(compo,1))<0)
                {
                  phi=3.*M_PI/2.;
                }

              U_P1_spherique(fa7,0)=sin(theta)*cos(phi)*U_P1_(fa7,0)+sin(theta)*sin(phi)*
                                    U_P1_(fa7,1)+cos(theta)*U_P1_(fa7,2);
              U_P1_spherique(fa7,1)=cos(theta)*cos(phi)*U_P1_(fa7,0)+cos(theta)*sin(phi)*
                                    U_P1_(fa7,1)-sin(theta)*U_P1_(fa7,2);
              U_P1_spherique(fa7,2)=sin(phi)*U_P1_(fa7,0)+cos(phi)*U_P1_(fa7,1);

              U_P2_spherique(fa7,0)=sin(theta)*cos(phi)*U_P2_(fa7,0)+sin(theta)*sin(phi)*
                                    U_P2_(fa7,1)+cos(theta)*U_P2_(fa7,2);
              U_P2_spherique(fa7,1)=cos(theta)*cos(phi)*U_P2_(fa7,0)+cos(theta)*sin(phi)*
                                    U_P2_(fa7,1)-sin(theta)*U_P2_(fa7,2);
              U_P2_spherique(fa7,2)=sin(phi)*U_P2_(fa7,0)+cos(phi)*U_P2_(fa7,1);

              U_cg_spherique(fa7,0)=sin(theta)*cos(phi)*vitesses_compo(compo,0)+sin(theta)*
                                    sin(phi)*vitesses_compo(compo,1)+cos(theta)*vitesses_compo(compo,2);
              U_cg_spherique(fa7,1)=cos(theta)*cos(phi)*vitesses_compo(compo,0)+cos(theta)*
                                    sin(phi)*vitesses_compo(compo,1)-sin(theta)*vitesses_compo(compo,2);
              U_cg_spherique(fa7,2)=sin(phi)*vitesses_compo(compo,0)+cos(phi)*
                                    vitesses_compo(compo,1);

              // we compute delta again because we need epsilon
              int elem_diph=domain.chercher_elements(gravity_center_fa7(fa7,0),
                                                     gravity_center_fa7(fa7,1),
                                                     gravity_center_fa7(fa7,2));

              DoubleVect delta_i(dimension);
              int elem00=face_voisins_for_interp(elem_faces_for_interp(elem_diph, 0+dimension),1);
              int elem11=face_voisins_for_interp(elem_faces_for_interp(elem_diph, 1+dimension),1);
              int elem22=face_voisins_for_interp(elem_faces_for_interp(elem_diph, 2+dimension),1);
              int elem33=face_voisins_for_interp(elem_faces_for_interp(elem_diph, 2),0);

              delta_i(0) = elem00>=0 ? fabs(domain_vdf.dist_elem(elem_diph, elem00, 0)) : -1e15;
              delta_i(1) = elem11>=0 ? fabs(domain_vdf.dist_elem(elem_diph, elem11, 1)) : -1e15;

              if (tab_fa7_normal(fa7,2)>0)
                delta_i(2) = elem22>=0 ? fabs(domain_vdf.dist_elem(elem_diph, elem22, 2)) : -1e15;
              else
                delta_i(2) = elem33>=0 ? fabs(domain_vdf.dist_elem(elem_diph, elem33 , 2)) : -1e15;

              double epsilon=0;
              // Interpolation distance varies with mesh refinement
              for (int dim=0; dim<dimension; dim++)
                epsilon+= fabs(delta_i(dim)*fabs(tab_fa7_normal(fa7,dim)));

              // We compute the components of the friction force in spherical coordinates
              // afet simplification: ff=mu*(2*Urr, Uthetar, Uphir)
              Urr(fa7)=(-U_P2_spherique(fa7,0)+4.*U_P1_spherique(fa7,0)-3.*
                        U_cg_spherique(fa7,0))/(2.*epsilon);
              Uthetar(fa7)=(-U_P2_spherique(fa7,1)+4.*U_P1_spherique(fa7,1)-3.*
                            U_cg_spherique(fa7,1))/(2.*epsilon);
              Uphir(fa7)=(-U_P2_spherique(fa7,2)+4.*U_P1_spherique(fa7,2)-3.*
                          U_cg_spherique(fa7,2))/(2.*epsilon);

              DoubleVect ff(dimension);
              ff(0)=mu_f*fa7_surface(fa7)*(2.*sin(theta)*cos(phi)*Urr(fa7)+cos(theta)*
                                           cos(phi)*Uthetar(fa7)-sin(phi)*Uphir(fa7));
              ff(1)=mu_f*fa7_surface(fa7)*(2.*sin(theta)*sin(phi)*Urr(fa7)+cos(theta)*
                                           sin(phi)*Uthetar(fa7)+cos(phi)*Uphir(fa7));
              ff(2)=mu_f*fa7_surface(fa7)*(2.*cos(theta)*Urr(fa7)-sin(theta)*Uthetar(fa7));

              for (int dim=0; dim<dimension; dim++)
                total_friction_force_(compo,dim)+=ff(dim);

              if (is_post_process_friction_force_fa7_)
                {
                  for (int dim=0; dim<dimension; dim++)
                    friction_force_fa7_(fa7,dim)=ff(dim);
                }
            }
        }
    }
}

void Post_Processing_Hydrodynamic_Forces::fill_absurd_value_coord_neighbor_fluid_fa7(int fa7)
{
  for (int dim=0; dim<dimension; dim++)
    {
      coord_neighbor_fluid_fa7_pressure_1_(fa7,dim)=-1e15;
      coord_neighbor_fluid_fa7_pressure_2_(fa7,dim)=-1e15;
      coord_neighbor_fluid_fa7_gradU_1_(fa7,dim)=-1e15;
      coord_neighbor_fluid_fa7_gradU_2_(fa7,dim)=-1e15;
    }
}

void Post_Processing_Hydrodynamic_Forces::compute_pressure_force_trilinear_linear(int nb_fa7,
                                                                                  const Maillage_FT_Disc& mesh,
                                                                                  Convection_Diffusion_Temperature_FT_Disc::
                                                                                  Thermal_correction_discretization_method
                                                                                  thermal_correction_discretization_method,
                                                                                  const IntVect& compo_connexes_fa7,
                                                                                  const ArrOfDouble& fa7_surface,
                                                                                  const DoubleTab& tab_fa7_normal)
{
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  DoubleTab pressure_P1(nb_fa7);
  DoubleTab pressure_P2(nb_fa7);

  if (trilinear_interpolation_elem(eq_ns.la_pression->valeurs(),
                                   coord_neighbor_fluid_fa7_pressure_1_, pressure_P1, 0, thermal_correction_discretization_method) &&
      trilinear_interpolation_elem(eq_ns.la_pression->valeurs(),
                                   coord_neighbor_fluid_fa7_pressure_2_, pressure_P2,1,
                                   Convection_Diffusion_Temperature_FT_Disc::
                                   Thermal_correction_discretization_method::P1))
    {
      DoubleTab extrapolated_pressure(nb_fa7);
      prop_P2_fluid_compo_=0;
      for (int fa7=0; fa7<nb_fa7; fa7++)
        {
          if (!mesh.facette_virtuelle(fa7))
            {
              int compo = compo_connexes_fa7(fa7);
              if (phase_indicator_function_P2_(fa7)==1) prop_P2_fluid_compo_(compo)+=1;

              if (pressure_P2(fa7)<-1e10)
                {
                  if (is_post_process_pressure_fa7_) pressure_fa7_(fa7)=1e15;
                  extrapolated_pressure(fa7)=1e15;
                  continue;
                }

              // 3*pressure_P2(compo,fa7)-2*pressure_P1(compo,fa7); //
              // If one cannot interpolate in P1 and P2, then one do not compute the pressure force
              extrapolated_pressure(fa7) = pressure_P2(fa7)-interpolation_distance_pressure_P2_*
                                           (pressure_P2(fa7)-pressure_P1(fa7))/(interpolation_distance_pressure_P2_-
                                                                                interpolation_distance_pressure_P1_);
              if (is_post_process_pressure_fa7_)
                pressure_fa7_(fa7)=extrapolated_pressure(fa7);
            }
          else
            {
              if (is_post_process_pressure_fa7_) pressure_fa7_(fa7)=1e15;
              extrapolated_pressure(fa7)=-1e15;
            }
        }
      for (int fa7=0; fa7<nb_fa7; fa7++)
        {
          int compo = compo_connexes_fa7(fa7);
          if (extrapolated_pressure(fa7)>1e10)
            continue;
          double coeff=-extrapolated_pressure(fa7)*fa7_surface(fa7);
          DoubleVect pressure_force_fa7(dimension);
          for (int dim=0; dim<dimension; dim++)
            pressure_force_fa7(dim)=coeff*tab_fa7_normal(fa7,dim);
          if (!mesh.facette_virtuelle(fa7))
            {
              total_surface_interf_(compo)+=fa7_surface(fa7);
              if (is_post_process_pressure_force_fa7_)
                {
                  for (int dim=0; dim<dimension; dim++)
                    pressure_force_fa7_(fa7,dim)=pressure_force_fa7(dim);
                }
              for (int dim=0; dim<dimension; dim++)
                total_pressure_force_(compo,dim)+=pressure_force_fa7(dim);
            }
        }
    }
}

void Post_Processing_Hydrodynamic_Forces::compute_neighbors_coordinates_fluid_fa7(const int nb_fa7,
                                                                                  const int is_discr_elem_diph,
                                                                                  const DoubleTab& gravity_center_fa7,
                                                                                  const Maillage_FT_Disc& mesh,
                                                                                  const DoubleTab& tab_fa7_normal,
                                                                                  const IntTab& particles_eulerian_id_number)
{
  const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
  const Domaine_VDF& domain_vdf = ref_cast(Domaine_VDF, eq_ns.domaine_dis());
  const Domaine& domain = domain_vdf.domaine();

  for (int fa7 =0 ; fa7<nb_fa7 ; fa7++)
    {
      if (!mesh.facette_virtuelle(fa7))
        {
          DoubleVect fa7_normal(dimension);
          int elem_diph=domain.chercher_elements(gravity_center_fa7(fa7,0),
                                                 gravity_center_fa7(fa7,1),
                                                 gravity_center_fa7(fa7,2));
          if (is_discr_elem_diph)
            {
              list_elem_diph_(fa7,0) = elem_diph;
              list_elem_diph_(fa7,1) = particles_eulerian_id_number(elem_diph);
            }
          if (elem_diph>=0)
            {
              DoubleVect delta_i(dimension);
              // One compute eulerian mesh thicknesses in which are located the lagrangian facets
              // If we have access, one compute the thickness outside of the particle
              // Otherwise, one compute the thickness inside of the particle
              // It just comes to the selection of the mesh juxtaposed to the two-phase cell.
              int acces=1;
              for (int dim=0; dim<dimension; dim++)
                {
                  int elem_haut=face_voisins_for_interp(
                                  elem_faces_for_interp(elem_diph, dim+dimension),1);
                  int elem_bas=face_voisins_for_interp(
                                 elem_faces_for_interp(elem_diph, dim),0);
                  if (elem_bas<0 && elem_haut<0)
                    acces=0;
                  if (tab_fa7_normal(fa7,dim)>0)
                    {
                      delta_i(dim) =  (elem_haut>=0) ? fabs(domain_vdf.dist_elem(elem_diph,elem_haut, dim)) :
                                      fabs(domain_vdf.dist_elem(elem_diph,elem_bas, dim));
                    }
                  else
                    {
                      delta_i(dim) =  (elem_bas>=0) ? fabs(domain_vdf.dist_elem(elem_diph,elem_bas, dim)) :
                                      fabs(domain_vdf.dist_elem(elem_diph,elem_haut, dim)) ;
                    }
                }
              if (acces==0)
                fill_absurd_value_coord_neighbor_fluid_fa7(fa7);
              else
                {
                  double epsilon=0;
                  for (int dim=0; dim<dimension; dim++)
                    {
                      // Interpolation distance varies with mesh refinement
                      epsilon+= fabs(delta_i(dim)*fabs(tab_fa7_normal(fa7,dim)));
                    }
                  for (int dim=0; dim<dimension; dim++)
                    {
                      fa7_normal(dim)=tab_fa7_normal(fa7,dim);
                      coord_neighbor_fluid_fa7_pressure_1_(fa7,dim)=gravity_center_fa7(fa7,dim)+
                                                                    interpolation_distance_pressure_P1_*epsilon*fa7_normal(dim);
                      coord_neighbor_fluid_fa7_pressure_2_(fa7,dim)=gravity_center_fa7(fa7,dim)+
                                                                    interpolation_distance_pressure_P2_*epsilon*fa7_normal(dim);
                      coord_neighbor_fluid_fa7_gradU_1_(fa7,dim)=gravity_center_fa7(fa7,dim)+
                                                                 interpolation_distance_pressure_P1_*epsilon*fa7_normal(dim);
                      coord_neighbor_fluid_fa7_gradU_2_(fa7,dim)=gravity_center_fa7(fa7,dim)+
                                                                 interpolation_distance_pressure_P2_*epsilon*fa7_normal(dim);
                    }
                }
            }
          else
            fill_absurd_value_coord_neighbor_fluid_fa7(fa7);
        }
    }
}

void Post_Processing_Hydrodynamic_Forces::compute_U_P2_moy(const int nb_particles_tot)
{
  for (int compo=0; compo<nb_particles_tot; compo++)
    {
      for (int dim=0; dim<dimension; dim++)
        {
          if (proportion_fa7_ok_UP2_(compo,dim)>0)
            U_P2_moy_(compo,dim)/=proportion_fa7_ok_UP2_(compo,dim);
        }
    }
}

void Post_Processing_Hydrodynamic_Forces::compute_proportion_fa7_ok_and_is_fluid_P2
(const int nb_particles_tot)
{
  for (int compo=0; compo<nb_particles_tot; compo++)
    {
      for (int dim=0; dim<dimension; dim++)
        {
          if (Nb_fa7_tot_par_compo_(compo)>0)
            proportion_fa7_ok_UP2_(compo,dim)/=Nb_fa7_tot_par_compo_(compo);
        }
      if (Nb_fa7_tot_par_compo_(compo)>0)
        prop_P2_fluid_compo_(compo)/=Nb_fa7_tot_par_compo_(compo);

    }
}

void Post_Processing_Hydrodynamic_Forces::compute_heat_transfer()
{
  if(!flag_heat_transfer_computation_)
    {
      Cerr << "Post_Processing_Hydrodynamic_Forces::compute_heat_transfer"  <<  finl;
      const Navier_Stokes_FT_Disc& eq_ns = ptr_eq_ns_.valeur();
      Transport_Interfaces_FT_Disc& eq_transport = ptr_eq_transport_.valeur();
      Convection_Diffusion_Temperature_FT_Disc& eq_temp = ptr_eq_temp_.valeur();
      const DoubleTab& temperature = eq_temp.inconnue().valeurs();
      const Maillage_FT_Disc& mesh = eq_transport.maillage_interface();
      const Domaine_VDF& domain_vdf = ref_cast(Domaine_VDF, eq_ns.domaine_dis());
      const Domaine& domain = domain_vdf.domaine();
      const Fluide_Diphasique& mon_fluide = eq_ns.fluide_diphasique();
      double lambda_f=mon_fluide.fluide_phase(1).conductivite().valeurs()(0, 0);
      const int nb_fa7 = mesh.nb_facettes();

      IntVect compo_connexes_fa7(nb_fa7);
      int n = search_connex_components_local_FT(mesh, compo_connexes_fa7);
      int nb_compo_tot=compute_global_connex_components_FT(mesh, compo_connexes_fa7, n);

      if (total_heat_transfer_.size_array() != nb_compo_tot)
        {
          total_heat_transfer_.resize(nb_compo_tot);
          total_heat_transfer_=0;
        }

      if (nb_fa7>0)
        {
          const ArrOfDouble& les_surfaces_fa7 = mesh.get_update_surface_facettes();
          const DoubleTab& les_normales_fa7 = mesh.get_update_normale_facettes();
          heat_transfer_fa7_.resize(nb_fa7);
          heat_transfer_fa7_=1e15;

          const DoubleTab& les_cg_fa7=mesh.get_gravity_center_fa7();
          coord_neighbor_fluid_fa7_temp_1_.resize(nb_fa7,dimension);
          coord_neighbor_fluid_fa7_temp_2_.resize(nb_fa7,dimension);

          for (int fa7 =0 ; fa7<nb_fa7 ; fa7++)
            {
              if (!mesh.facette_virtuelle(fa7))
                {
                  DoubleVect normale_fa7(dimension);
                  int elem_diph=domain.chercher_elements(les_cg_fa7(fa7,0),
                                                         les_cg_fa7(fa7,1),les_cg_fa7(fa7,2));
                  DoubleVect delta_i(dimension);
                  for (int dim=0; dim<dimension; dim++)
                    {
                      int elem_haut=face_voisins_for_interp(elem_faces_for_interp(elem_diph,
                                                                                  dim+dimension),1);
                      int elem_bas=face_voisins_for_interp(elem_faces_for_interp(elem_diph, dim),0);
                      if (les_normales_fa7(fa7,dim)>0)
                        delta_i(dim) =  (elem_haut>=0) ? fabs(domain_vdf.dist_elem(elem_diph,
                                                                                   elem_haut, dim)) : fabs(domain_vdf.dist_elem(elem_diph,elem_bas, dim));
                      else delta_i(dim) =  (elem_bas>=0) ? fabs(domain_vdf.dist_elem(elem_diph,
                                                                                       elem_bas, dim)) : fabs(domain_vdf.dist_elem(elem_diph,elem_haut, dim));
                    }
                  double epsilon=0;
                  for (int dim=0; dim<dimension; dim++)
                    epsilon+= fabs(delta_i(dim)*fabs(les_normales_fa7(fa7,dim)));
                  for (int dim=0; dim<dimension; dim++)
                    {
                      normale_fa7(dim)=les_normales_fa7(fa7,dim);
                      coord_neighbor_fluid_fa7_temp_1_(fa7,dim)=les_cg_fa7(fa7,dim)+
                                                                interpolation_distance_temperature_P1_*epsilon*normale_fa7(dim);
                      coord_neighbor_fluid_fa7_temp_2_(fa7,dim)=les_cg_fa7(fa7,dim)+
                                                                interpolation_distance_temperature_P2_*epsilon*normale_fa7(dim);
                    }
                }
            }

          DoubleTab temp_P1(nb_fa7);
          DoubleTab temp_P2(nb_fa7);

          int interp_T_P1_ok=trilinear_interpolation_elem(temperature,
                                                          coord_neighbor_fluid_fa7_temp_1_, temp_P1);
          int interp_T_P2_ok=trilinear_interpolation_elem(temperature,
                                                          coord_neighbor_fluid_fa7_temp_2_, temp_P2);
          if (interp_T_P1_ok &&  interp_T_P2_ok)
            {
              for (int fa7=0; fa7<nb_fa7; fa7++)
                {
                  int compo=compo_connexes_fa7(fa7);
                  if (!mesh.facette_virtuelle(fa7))
                    {
                      int elem_diph=domain.chercher_elements(les_cg_fa7(fa7,0),
                                                             les_cg_fa7(fa7,1),les_cg_fa7(fa7,2));
                      DoubleVect delta_i(dimension);
                      int elem00=face_voisins_for_interp(elem_faces_for_interp(elem_diph, 0+dimension),1);
                      int elem11=face_voisins_for_interp(elem_faces_for_interp(elem_diph, 1+dimension),1);
                      int elem22=face_voisins_for_interp(elem_faces_for_interp(elem_diph, 2+dimension),1);
                      int elem33=face_voisins_for_interp(elem_faces_for_interp(elem_diph, 2),0);

                      delta_i(0) = elem00>=0 ? fabs(domain_vdf.dist_elem(elem_diph, elem00, 0)) : -1e15;
                      delta_i(1) = elem11>=0 ? fabs(domain_vdf.dist_elem(elem_diph, elem11, 1)) : -1e15;

                      if (les_normales_fa7(fa7,2)>0)
                        delta_i(2) = elem22>=0 ? fabs(domain_vdf.dist_elem(elem_diph, elem22, 2)) : -1e15;
                      else
                        delta_i(2) = elem33>=0 ? fabs(domain_vdf.dist_elem(elem_diph, elem33 , 2)) : -1e15;

                      double epsilon=0;
                      for (int dim=0; dim<dimension; dim++)
                        epsilon+= fabs(delta_i(dim)*fabs(les_normales_fa7(fa7,dim))); // la distance d'interpolation varie en fonction du raffinement du maillage
                      heat_transfer_fa7_(fa7)=lambda_f*(-temp_P2(fa7)+4.*temp_P1(fa7)-3.*eq_temp.get_tsat_constant())/(2.*epsilon)*les_surfaces_fa7(fa7); // schema decentre avant d'ordre 2
                      total_heat_transfer_(compo)+=heat_transfer_fa7_(fa7);
                    }
                }
            }
        }
      mp_sum_for_each_item(total_heat_transfer_);
      raise_the_flag_heat_transfer();
    }
}

