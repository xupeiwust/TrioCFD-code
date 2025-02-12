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
*****************************************************************************/

#include <Collision_Model_FT.h>
#include <EFichier.h>
#include <EcritureLectureSpecial.h>
#include <Probleme_FT_Disc_gen.h>

Implemente_instanciable_sans_constructeur(Collision_Model_FT,"Collision_Model_FT",Objet_U);

Collision_Model_FT::Collision_Model_FT()
{
  fictive_wall_coordinates_.resize(2*dimension);
}

Entree& Collision_Model_FT::readOn (Entree& is)
{
  Param p(que_suis_je());
  set_param(p);
  p.lire_avec_accolades_depuis(is);
  return is;
}

Sortie& Collision_Model_FT::printOn(Sortie& os) const
{
  Cerr << "Erreur : ::printOn n'est pas code." << finl;
  assert(0);
  return os;
}

void Collision_Model_FT::set_param(Param& p)
{
  p.ajouter_non_std("collision_model", (this),Param::REQUIRED);
  p.ajouter_non_std("collision_detection", (this),Param::REQUIRED);
  p.ajouter("collision_duration", &collision_duration_, Param::REQUIRED); // XD_ADD_P duration of the collision in seconds;
  p.ajouter("activate_collision_before_impact", &is_collision_activated_before_impact_, Param::REQUIRED);
  p.ajouter("activation_distance_percentage_diameter", &activation_distance_percentage_diameter_, Param::REQUIRED);
  p.ajouter_flag("force_on_two_phase_elem", &is_force_on_two_phase_elem_);
}

int Collision_Model_FT::lire_motcle_non_standard(const Motcle& word, Entree& is)
{
  if (word=="collision_model")
    {
      Motcles words;
      words.add("hybrid_esi");
      words.add("breugem");
      Motcle secondword;
      is >> secondword;
      Cerr << "Reading collision_model attributes: " << secondword << finl;
      const int r = words.search(secondword);
      switch(r)
        {
        case 0:
          collision_model_ = Collision_Model_FT::HYBRID_ESI;
          break;
        case 1:
          collision_model_ = Collision_Model_FT::BREUGEM;
          break;
        default:
          Cerr << "Error " << words << "was expected whereas " << secondword <<" has been found."<< finl;
          barrier();
          exit();
        }
      return 1;
    }
  else if (word=="collision_detection")
    {
      Motcles words;
      words.add("use_Verlet_tables");
      words.add("detection_thickness_Verlet");
      words.add("nb_pas_dt_max");
      words.add("activate_linked_cell");
      Motcle openbrace ="{";
      Motcle closedbrace="}";
      Motcle secondword;
      is >> secondword;
      if (secondword==openbrace)
        {
          is >> secondword;
          while (secondword != openbrace)
            {
              int rang2 = words.search(secondword);
              switch(rang2)
                {
                case 0:
                  is >> is_detection_Verlet_;
                  break;
                case 1:
                  is >> detection_thickness_Verlet_;
                  break;
                case 2:
                  is >> nb_pas_dt_max_Verlet_;
                  break;
                case 3:
                  is >> is_linked_cell_activated_;
                  break;
                default:
                  Cerr << "Collision_Model_FT::lire_motcle_non_standard\n"
                       << " options of collision_detection are:\n"
                       << words;
                  exit();
                }
              is >> secondword;
            }
        }
      return 1;
    }
  else
    {
      Cerr << word << " is not a keyword understood by " << que_suis_je() << " in lire_motcle_non_standard"<< finl;
      exit();
    }
  return -1;
}

void write_table_collision_model(Sortie& os, const DoubleTab& tab)
{
  const int dim0 = tab.dimension(0);
  if (Process::je_suis_maitre())
    os << dim0 << tspace << tab.dimension(1) << finl;
  os.put(tab.addr(), tab.size_array());
  os.syncfile();
}

void read_table_collision_model(Entree& is, DoubleTab& tab, Entree * fichier)
{
  int dim0;
  int dim1;
  (*fichier)  >> dim0  >> dim1;
  DoubleTab tmp;
  tmp.resize(dim0,dim1);
  fichier->get(tmp.addr(), tmp.size_array());
  tab=tmp;
}

void Collision_Model_FT::reset()
{
  int nb_compo_tot = nb_compo_tot_;
  int nb_boundaries=2*dimension;
  F_old_.resize(nb_compo_tot,nb_compo_tot+nb_boundaries);
  stiffness_.resize(nb_compo_tot,nb_compo_tot+nb_boundaries);
  e_eff_.resize(nb_compo_tot,nb_compo_tot+nb_boundaries);
}


void open_file_collision(SFichier& os, const int& flag, const Transport_Interfaces_FT_Disc& equation, Nom backup_file)
{
  // null flag, do not open the file
  if (flag==0)
    return ;

  const Schema_Temps_base& sch=equation.probleme().schema_temps();
  const int& precision=sch.precision_impr();

  os.ouvrir(backup_file,std::_S_out);
  os.precision(precision);
  os.setf(ios::scientific);
}

void open_file_collision(EFichier& os, const int& flag, const Transport_Interfaces_FT_Disc& equation, Nom restart_file)
{
  // null flag, do not open the file
  if (flag==0)
    return ;

  const Schema_Temps_base& sch=equation.probleme().schema_temps();
  const int& precision=sch.precision_impr();

  os.ouvrir(restart_file,ios::app);
  os.precision(precision);
  os.setf(ios::scientific);
}

int Collision_Model_FT::reprendre(Entree& is)
{
  Nom readword;
  const int format_xyz = EcritureLectureSpecial::is_lecture_special();
  reset();
  if (format_xyz)
    {
      if (Process::je_suis_maitre())
        {
          EFichier file_collision_model;
          open_file_collision(file_collision_model,1,refequation_transport_.valeur(),filename_data_fpi_);
          file_collision_model >> readword;
          if (readword != que_suis_je())
            {
              Cerr << "Error in Collision_Model_FT::reprendre\n";
              Cerr << "We was expecting " << que_suis_je();
              Cerr << "\n We found " << readword << finl;
              Process::exit();
            }
          file_collision_model >> readword;
          F_old_.lit(file_collision_model);
          file_collision_model >> readword;
          stiffness_.lit(file_collision_model);
          file_collision_model >> readword;
          e_eff_.lit(file_collision_model);
          file_collision_model.close();
        }
      envoyer_broadcast(F_old_,0);
      envoyer_broadcast(stiffness_,0);
      envoyer_broadcast(e_eff_,0);
      barrier();
      return 1;
    }
  else
    {
      is >> readword;
      if (readword != que_suis_je())
        {
          Cerr << "Error in Collision_Model_FT::reprendre\n";
          Cerr << "We was expecting " << que_suis_je();
          Cerr << "\n We found " << readword << finl;
          Process::exit();
        }
      is >> readword;
      F_old_.lit(is);
      is >> readword;
      stiffness_.lit(is);
      is >> readword;
      e_eff_.lit(is);
    }
  return 1;
}

int Collision_Model_FT::sauvegarder(Sortie& os) const
{
  int special, afaire;
  const int format_xyz = EcritureLectureSpecial::is_ecriture_special(special, afaire);

  if (format_xyz)
    {
      if (Process::je_suis_maitre())
        {
          SFichier file_collision_model;
          open_file_collision(file_collision_model,1,refequation_transport_.valeur(),filename_data_fpi_);
          file_collision_model << Nom(que_suis_je()) << finl;
          file_collision_model << "F_old" << finl;
          F_old_.ecrit(file_collision_model);
          file_collision_model << "stiffness" << finl;
          stiffness_.ecrit(file_collision_model);
          file_collision_model << "e_eff" << finl;
          e_eff_.ecrit(file_collision_model);
          file_collision_model.close();
        }
      return 0;
    }
  else
    {
      int bytes = 0;
      os << que_suis_je() << finl;
      os << "F_old" << finl;
      F_old_.ecrit(os);
      bytes += 8 * F_old_.size_array();
      os << "stiffness" << finl;
      stiffness_.ecrit(os);
      bytes += 8 * stiffness_.size_array();
      os << "e_eff" << finl;
      e_eff_.ecrit(os);
      bytes += 8 * e_eff_.size_array();
      return bytes;
    }
  return 0;
}


void  Collision_Model_FT::compute_contact_force(DoubleTab& force_contact, int& isFirstStepOfCollision, double& dist_int, double& next_dist_int, DoubleTab& norm, DoubleTab& dUn, double& masse_eff, int& compo, int& voisin, double& Stb, double& ed, double& vitesseRelNorm, double& dt, double& prod_scal)
{
  switch (collision_model_)
    {
    case Collision_Model_FT::HYBRID_ESI: // see Hamidi et al., IJMF, 2023.
      {
        DoubleTab& tab_stiffness=get_stiffness();
        DoubleTab& e_eff=get_e_eff();
        if (isFirstStepOfCollision) // we compute the stiffness based on the impact velocity
          {
            tab_stiffness(compo,voisin)= compute_stiffness_breugem(masse_eff,ed);
            e_eff(compo,voisin)=ed *compute_ewet_legendre(Stb);
          }
        int is_compression_step = prod_scal <= 0; // if true : the particle is approaching, if false, the particle is moving away
        double e_eff_part = is_compression_step ? 1 : e_eff(compo, voisin);
        double stiffness = tab_stiffness(compo, voisin);

        for (int d = 0; d < dimension; d++)
          force_contact(d)= -pow(e_eff_part,2) * stiffness * next_dist_int * norm(d);
      }
      break;
    case Collision_Model_FT::BREUGEM: // See. W-P. Breugem, 2010.
      {
        double stiffness = compute_stiffness_breugem(masse_eff,ed);
        double damper = compute_damper_breugem(masse_eff,ed);

        for (int d = 0; d < dimension; d++)
          force_contact(d)= -stiffness*next_dist_int*norm(d) -damper*dUn(d);
      }
      break;

    default:
      Cerr << "The method specified for modele_collision in not recognized. \n" << finl;
      Process::exit();
    }
}

void Collision_Model_FT::compute_fictive_wall_coordinates(const double& radius)
{
  DoubleVect offset_values(2*dimension);
  const double diameter=2*radius;

  switch(is_collision_activated_before_impact_)
    {
    case 0:
      offset_values=0;
      break;
    case 1:
      if (activation_distance_percentage_diameter_>0) offset_values=activation_distance_percentage_diameter_*diameter/100;
      break;
    default:
      Cerr << "Collision_Model_FT::compute_fictive_wall_coordinates error"  <<finl;
      Process::exit();
      break;
    }

  // an offset is computed for the walls to activate the collision process before the impact
  fictive_wall_coordinates_(0) = origin_(0) - radius + offset_values(0);
  fictive_wall_coordinates_(1) = origin_(1) - radius + offset_values(1);
  fictive_wall_coordinates_(2) = origin_(2) - radius + offset_values(2);
  fictive_wall_coordinates_(3) = origin_(0) + domain_dimensions_(0) + radius - offset_values(3);
  fictive_wall_coordinates_(4) = origin_(1) + domain_dimensions_(1) + radius - offset_values(4);
  fictive_wall_coordinates_(5) = origin_(2) + domain_dimensions_(2) + radius - offset_values(5);
}

void Collision_Model_FT::resize_geometric_parameters()
{
  domain_dimensions_.resize(dimension);
  nb_nodes_.resize(dimension);
  origin_.resize(dimension);
}

int Collision_Model_FT::check_for_duplicates(ArrOfInt& vector)
{
  int flag =0;
  ArrOfInt copy_vector(vector);
  const int size = copy_vector.size_array();
  copy_vector.ordonne_array();
  for (int i = 0; i < size-1; i++)
    {
      if (copy_vector(i)==copy_vector(i+1))
        {
          flag = 1;
          Cerr << copy_vector(i) << " is duplicate !!" << finl ;
        }
    }
  return flag;
}

/*! @brief Recover the geometric parameters of the domain:
 * number of nodes in each direction
 * origin of the domain
 * dimensions of the domain
 */
void Collision_Model_FT::set_geometric_parameters(Domaine_VDF& domaine_vdf)
{
  const Domaine& domain = domaine_vdf.domaine();
  DoubleTab BB=domain.getBoundingBox();
  const Bords& bords=domain.faces_bord();
  DoubleVect NiNj(dimension); // NiNj=(NyNz NxNz NxNy )
  resize_geometric_parameters();

  // 1. Number of nodes per direction
  NiNj=0;
  for (int i=0; i<bords.nb_bords(); i++)
    {
      int nb_boundary_faces = mp_sum(ref_cast(Frontiere,bords(i)).nb_faces());
      int nb_boundary_faces_local=ref_cast(Frontiere,bords(i)).nb_faces();
      if (nb_boundary_faces_local>0)
        {
          int face1=ref_cast(Frontiere,bords(i)).num_premiere_face();
          int orientation_face1=domaine_vdf.orientation(face1);
          NiNj(orientation_face1)=nb_boundary_faces;
        }
    }
  long long NxNy= static_cast<int>(mp_max(NiNj(2))) ;
  long long NxNz= static_cast<int>(mp_max(NiNj(1)));
  long long NyNz= static_cast<int>(mp_max(NiNj(0)));

  int Nx,Ny,Nz;

  Ny= NxNz>0 ? static_cast<int>(sqrt(NxNy*NyNz/NxNz)) : 0; // nb elem in the y-direction
  Nz= NxNy>0 ? static_cast<int>((NxNz*Ny/NxNy)) : 0; // nb elem in the z-direction, WARNING: operations order matter because Nz is an integer
  Nx= Ny>0 ? static_cast<int>(NxNy/Ny) : 0; // nb elem in the x-direction

  nb_nodes_(0)=Nx++; // nb nodes in the x-direction
  nb_nodes_(1)=Ny++; // nb nodes in the y-direction
  nb_nodes_(2)=Nz++; // nb nodes in the z-direction

  // 2. Origin and Domain_dimensions
  DoubleVect Origin(dimension);
  DoubleVect Domain_dimensions(dimension);
  Origin=0.;
  Domain_dimensions=0.;

  for (int j=0; j<dimension; j++)
    {
      double min_ = mp_min(BB(j,0));
      double max_ = mp_max(BB(j,1));

	  Origin(j)=min_;
	  Domain_dimensions(j)=max_-min_;
    }

  set_origin(Origin);
  set_domain_dimensions(Domain_dimensions);

  Cerr << "Origin " << Origin << finl;
  Cerr << "Domain length" << Domain_dimensions << finl;
  Cerr << "Nx Ny Nz " << nb_nodes_ << finl;
}

void Collision_Model_FT::associate_transport_equation(const Equation_base& equation)
{
  const Transport_Interfaces_FT_Disc& eq = ref_cast(Transport_Interfaces_FT_Disc,equation);
  refequation_transport_ = eq;
}


