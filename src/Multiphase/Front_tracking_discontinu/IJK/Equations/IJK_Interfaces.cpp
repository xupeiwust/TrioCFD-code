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

#include <Connex_components.h>
#include <IJK_Field_vector.h>
#include <Connex_components_FT.h>
#include <EcrFicPartageBin.h>
#include <Probleme_FTD_IJK_base.h>
#include <Probleme_FTD_IJK_cut_cell.h>
#include <Cut_cell_tools.h>
#include <IJK_Interfaces.h>
#include <IJK_Lata_writer.h>
#include <IJK_Navier_Stokes_tools.h>
#include <Cut_cell_tools.h>
#include <LataDB.h>
#include <LecFicDiffuse_JDD.h>
#include <Linear_algebra_tools.h>
#include <Octree_Double.h>
#include <Ouvrir_fichier.h>
#include <Param.h>
#include <SFichier.h>
#include <Domaine_VF.h>
#include <algorithm>
#include <communications.h>
#include <iostream>
#include <memory>
#include <stat_counters.h>
#include <Cut_field.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Navier_Stokes_FTD_IJK.h>
#include <vector>
#include <map>

// #define GB_VERBOSE

// permet de preciser la methode utilisee pour mettre a jour l'indicatrice
// 0: methode classique avec recherche des composantes connexes
// 1: methode optimisee
#define CLASSIC_METHOD 1
#define BORD -10000

Implemente_instanciable_sans_constructeur(IJK_Interfaces, "IJK_Interfaces", Equation_base);

IJK_Interfaces::IJK_Interfaces()
{
  positions_reference_.resize(0); // Par defaut, a dimensionner ensuite
  mean_force_.resize(0);          // Par defaut, a dimensionner ensuite
  compo_to_group_.resize(0); // Par defaut, a dimensionner ensuite par nb_bulles
}

static void FixedVector_to_DoubleTab(const FixedVector<ArrOfDouble,3> fixed_arr, DoubleTab& tab)
{
  int size_fixed_arr[3] = {0,0,0};
  for (int c=0; c<3; c++)
    size_fixed_arr[c] = fixed_arr[c].size_array();
  assert(size_fixed_arr[0] == size_fixed_arr[1]);
  assert(size_fixed_arr[0] == size_fixed_arr[2]);
  tab.reset();
  tab.resize(size_fixed_arr[0],3);
  for (int c=0; c<3; c++)
    for (int l=0; l<size_fixed_arr[0]; l++)
      tab(l,c) = fixed_arr[c](l);
}

void IJK_Interfaces::set_param_reprise_pb(Param& param)
{
  param.ajouter("interfaces", &(probleme_ijk().get_interface())); // on lit dans l'eq les params reprise
}

/*! Output the FT mesh information in the master LATA file.
 *
 * ASSUMPTION: for now the FT mesh always fits within 32b (in terms of nb of elems / vertices). This is checked.
 *
 * Ajoute ceci dans le fichier lata maitre:
 *  GEOM meshname type_elem=TRIANGLE_3D
 *  CHAMP SOMMETS  filename.step.meshname.SOMMETS geometry=meshname size=... composantes=3
 *  CHAMP ELEMENTS filename.step.meshname.ELEMENTS geometry=meshname size=... composantes=3 format=INT32|64
 */
void IJK_Interfaces::dumplata_ft_mesh(const char *filename, const char *meshname, int step) const
{
  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  Nom prefix = Nom(filename) + Nom(".") + Nom(step) + Nom(".") + Nom(meshname) + Nom(".");
  Nom fdsom = prefix + Nom("SOMMETS");
  Nom fdelem = prefix + Nom("ELEMENTS");
  Nom fdpelocal = prefix + Nom("FACETTE_PE_LOCAL");
  Nom fdpeowner = prefix + Nom("FACETTE_PE_OWNER");
  Nom fd_sommet_reel = prefix + Nom("INDEX_SOMMET_REEL");

  const int nb_som = mesh.nb_sommets();
  const int nb_elem = mesh.nb_facettes();
  const int nsomtot = Process::check_int_overflow(Process::mp_sum(nb_som));
  const int nelemtot = Process::check_int_overflow(Process::mp_sum(nb_elem));
  // valeur a ajouter a un indice local de sommet pour obtenir un indice global
  // de sommet
  const int offset_sommet = Process::check_int_overflow(Process::mppartial_sum(nb_som));
  FloatTab tmp(nb_som, 3);
  const DoubleTab& sommets = mesh.sommets();

  for (int i = 0; i < nb_som; i++)
    for (int j = 0; j < 3; j++)
      tmp(i, j) = (float)sommets(i, j);
  EcrFicPartageBin file;

  // See top doc above, we do everything in 32b
  file.set_64b(false);

  file.ouvrir(fdsom);
  file.put(tmp.addr(), tmp.size_array(), 3);
  file.syncfile();
  file.close();

  // Tableau des facettes locales ecrites dans le fichier: contient les indices
  // globaux des sommets
  //  (ajout de offset_sommet a l'indice du sommet local)
  IntTab tmp_facettes_renum(nb_elem, 3);
  const IntTab& facettes = mesh.facettes();
  for (int i = 0; i < nb_elem; i++)
    for (int j = 0; j < 3; j++)
      tmp_facettes_renum(i, j) = facettes(i, j) + offset_sommet;
  file.ouvrir(fdelem);
  file.put(tmp_facettes_renum.addr(), tmp_facettes_renum.size_array(), 3);
  file.syncfile();
  file.close();

  // Postraitement du numero du processeur local
  file.ouvrir(fdpelocal);
  ArrOfInt tmp2(nb_elem);
  tmp2 = Process::me();
  file.put(tmp2.addr(), nb_elem, 1);
  file.syncfile();
  file.close();

  // Postraitement du numero du processeur proprietaire
  file.ouvrir(fdpeowner);
  mesh.facette_PE_owner(tmp2);
  file.put(tmp2.addr(), nb_elem, 1);
  file.syncfile();
  file.close();

  // Postraitement de l'indice global du sommet reel associe a chaque sommet
  // postraite
  {
    ArrOfInt indice_global(nb_som);
    const int np = Process::nproc();
    // Pour chaque processeur, quel est l'offset des sommets de ce processeur:
    ArrOfInt offset_received(np);
    ArrOfInt offset_to_send(np);
    // On va faire un echange: chaque proc envoye a tous les processeurs sont
    // offset Je remplis le tableau avec ma valeur:
    offset_to_send = offset_sommet;
    // Je l'envoie a tout le monde:
    envoyer_all_to_all(offset_to_send, offset_received);
    // Maintenant offset_received[i] contient la valeur offset_sommet du
    // processeur i
    const ArrOfInt& sommet_num_owner = mesh.sommet_num_owner();
    const ArrOfInt& sommet_PE_owner = mesh.sommet_PE_owner();
    // Calcul de l'indice global du sommet reel correspondant a chaque sommet
    // local:
    for (int i = 0; i < nb_som; i++)
      {
        const int indice_local_sur_proc_proprietaire = sommet_num_owner[i];
        const int proc_proprietaire = sommet_PE_owner[i];
        const int offset_proc_proprietaire = offset_received[proc_proprietaire];
        indice_global[i] = indice_local_sur_proc_proprietaire + offset_proc_proprietaire;
      }
    file.ouvrir(fd_sommet_reel);
    file.put(indice_global.addr(), indice_global.size_array(), 1);
    file.syncfile();
    file.close();
  }

  const Nom format = "INT32";
  if (Process::je_suis_maitre())
    {
      SFichier master_file;
      master_file.ouvrir(filename, ios::app);
      master_file << "Geom " << meshname << " type_elem=TRIANGLE_3D" << finl;

      master_file << "Champ SOMMETS " << basename(fdsom) << " geometrie=" << meshname;
      master_file << " size=" << nsomtot << " composantes=3" << finl;

      master_file << "Champ ELEMENTS " << basename(fdelem) << " geometrie=" << meshname;
      master_file << " size=" << nelemtot << " composantes=3 format=" << format << finl;

      master_file << "Champ FACETTE_PE_LOCAL " << basename(fdpelocal) << " geometrie=" << meshname;
      master_file << " size=" << nelemtot << " composantes=1 format=" << format << " localisation=ELEM" << finl;

      master_file << "Champ FACETTE_PE_OWNER " << basename(fdpeowner) << " geometrie=" << meshname;
      master_file << " size=" << nelemtot << " composantes=1 format=" << format << " localisation=ELEM" << finl;

      master_file << "Champ INDEX_SOMMET_REEL " << basename(fd_sommet_reel) << " geometrie=" << meshname;
      master_file << " size=" << nsomtot << " composantes=1 format=" << format << " localisation=SOM" << finl;
    }
}

// Copie de la methode precedente mais pour les DoubleTab aux sommets :
void runge_kutta3_update(const DoubleTab& dvi, DoubleTab& G, DoubleTab& l, const int step, const double dt_tot,
                         const Maillage_FT_IJK& maillage)
{
  const double coeff_a[3] = {0., -5. / 9., -153. / 128.};
  // Gk[0] = 1; Gk[i+1] = Gk[i] * a[i+1] + 1
  const double coeff_Gk[3] = {1., 4. / 9., 15. / 32.};

  const double facteurG = coeff_a[step];
  const double intermediate_dt = compute_fractionnal_timestep_rk3(dt_tot, step);
  const double delta_t_divided_by_Gk = intermediate_dt / coeff_Gk[step];
  const int nbsom = maillage.nb_sommets();

  // Resize du tableau
  G.resize(nbsom, 3);

  switch (step)
    {
    case 0:
      // don't read initial value of G (no performance benefit because write to G
      // causes the processor to fetch the cache line, but we don't wand to use a
      // potentially uninitialized value
      for (int isom = 0; isom < nbsom; isom++)
        {
          for (int dir = 0; dir < 3; dir++)
            {
              if (maillage.sommet_virtuel(isom))
                {
                  G(isom, dir) = 111. / intermediate_dt;
                  l(isom, dir) = 111. / intermediate_dt;
                }
              double x = dvi(isom, dir);
              G(isom, dir) = x;
              l(isom, dir) += x * delta_t_divided_by_Gk;
            }
        }
      break;
    case 1:
      // general case, read and write G
      for (int isom = 0; isom < nbsom; isom++)
        {
          for (int dir = 0; dir < 3; dir++)
            {
              if (maillage.sommet_virtuel(isom))
                {
                  G(isom, dir) = 111. / intermediate_dt;
                  l(isom, dir) = 111. / intermediate_dt;
                }
              double x = G(isom, dir) * facteurG + dvi(isom, dir);
              G(isom, dir) = x;
              l(isom, dir) += x * delta_t_divided_by_Gk;
            }
        }
      break;
    case 2:
      // do not write G
      for (int isom = 0; isom < nbsom; isom++)
        {
          for (int dir = 0; dir < 3; dir++)
            {
              if (maillage.sommet_virtuel(isom))
                {
                  G(isom, dir) = 111. / intermediate_dt;
                  l(isom, dir) = 111. / intermediate_dt;
                }
              double x = G(isom, dir) * facteurG + dvi(isom, dir);
              l(isom, dir) += x * delta_t_divided_by_Gk;
            }
        }
      break;
    default:
      Cerr << "Error in runge_kutta_update: wrong step" << finl;
      Process::exit();
    };
}

Sortie& IJK_Interfaces::printOn(Sortie& os) const
{
  os << "{\n"
     << "   fichier_reprise_interface " << basename(fichier_sauvegarde_interface_) << "\n"
     << "   timestep_reprise_interface " << timestep_sauvegarde_interface_
     << "\n"
     //     << "   lata_meshname " << nom_par_defaut_interfaces << "\n"
     << "   remaillage_ft_ijk " << remaillage_ft_ijk_;
  if (terme_gravite_ == GRAVITE_RHO_G)
    os << "   terme_gravite rho_g \n";
  else
    os << "   terme_gravite grad_i \n";

  if (correction_gradient_potentiel_)
    os << "   correction_gradient_potentiel \n";
  if (delta_p_max_repulsion_ > 0.)
    {
      os << "   portee_force_repulsion " << portee_force_repulsion_ << "\n"
         << "   delta_p_max_repulsion " << delta_p_max_repulsion_ << "\n";
    }
  if (delta_p_wall_max_repulsion_ > 0.)
    {
      os << "   portee_wall_repulsion " << portee_wall_repulsion_ << "\n"
         << "   delta_p_wall_max_repulsion " << delta_p_wall_max_repulsion_ << "\n";
    }
  if (active_repulsion_paroi_)
    os << "   active_repulsion_paroi \n";
  if (frozen_)
    os << "   frozen \n" ;
  if (follow_colors_)
    {
      // Sauvegarde des couleurs des bulles (on imprime le tableau) :
      Cerr << "IJK_Interfaces::printOn -- couleurs des bulles stockees " << finl;
      os << "   follow_colors " <<  "\n"
         << "   reprise_colors " << through_yminus_  ;
    }
  if (compo_to_group_.size_array()>0)
    {
      Cerr << "IJK_Interfaces::printOn -- Groupes d'appartenance stockes " << finl;
      os << "  bubble_groups " << compo_to_group_;
    }
  if (flag_positions_reference_)
    os << "  positions_reference " << positions_reference_;

  if (parcours_.get_correction_parcours_thomas())
    os << "  parcours_interface { correction_parcours_thomas } " << "\n";
  if (parcours_.get_parcours_sans_tolerance())
    os << "  parcours_interface { parcours_sans_tolerance } " << "\n";

  if (use_barycentres_velocity_)
    os << "use_barycentres_velocity" << "\n";
  if (read_barycentres_velocity_)
    os << "use_barycentres_velocity" << "\n";

  double max_force_compo = 0.;
  if (mean_force_.size_array() > 0)
    max_force_compo = max_abs_array(mean_force_);
  if (max_force_compo > 1.e-16)
    {
      // On a quelque part une force
      Cerr << "IJK_Interfaces::printOn -- Storing mean into sauv for future "
           "restart. (max_force_compo= "
           << max_force_compo << " )." << finl;
      os << "  mean_force " << mean_force_;
    }
  os << " }\n";

  return os;
}
// XD interfaces interprete nul 1 not_set
Entree& IJK_Interfaces::readOn(Entree& is)
{
  Param param(que_suis_je());
  lata_interfaces_meshname_ = "INTERFACES"; // This line is necessary for reprendre_probleme in Probleme_FTD_IJK_base

  param.ajouter("bubble_groups", &compo_to_group_);
  param.ajouter("fichier_reprise_interface", &fichier_reprise_interface_, Param::REQUIRED); // XD_ADD_P  chaine not_set
  param.ajouter("timestep_reprise_interface", &timestep_reprise_interface_); // XD_ADD_P entier not_set
  param.ajouter("lata_meshname", &lata_interfaces_meshname_); // XD_ADD_P chaine not_set
  param.ajouter("remaillage_ft_ijk", &remaillage_ft_ijk_); // XD_ADD_P remaillage_ft_ijk not_set
  param.ajouter_flag("use_tryggvason_interfacial_source", &use_tryggvason_interfacial_source_); // XD_ADD_P remaillage_ft_ijk not_set
  param.ajouter("portee_force_repulsion", &portee_force_repulsion_);
  param.ajouter("delta_p_max_repulsion", &delta_p_max_repulsion_);
  param.ajouter("portee_wall_repulsion", &portee_wall_repulsion_);
  param.ajouter("delta_p_wall_max_repulsion", &delta_p_wall_max_repulsion_);
  param.ajouter_flag("active_repulsion_paroi", &active_repulsion_paroi_);
  param.ajouter("no_octree_method", &no_octree_method_);  // XD_ADD_P entier if the bubbles repel each other, what method should be used to compute relative velocities? Octree method by default, otherwise we used the IJK discretization
  param.ajouter_flag("follow_colors", &follow_colors_);
  param.ajouter_flag("compute_distance_autres_interfaces", &compute_distance_autres_interfaces_); // XD_ADD_P rien not_set
  param.ajouter("reprise_colors", &through_yminus_);
  param.ajouter_flag("correction_gradient_potentiel", &correction_gradient_potentiel_);
  param.ajouter_flag("avoid_duplicata", &avoid_duplicata_);
  param.ajouter("factor_length_duplicata", &factor_length_duplicata_);
  param.ajouter("ncells_forbidden", &ncells_forbidden_);
  param.ajouter("ncells_deleted", &ncells_deleted_);
  param.ajouter_flag("frozen", &frozen_);
  param.ajouter("positions_reference", &positions_reference_);
  param.ajouter("mean_force", &mean_force_);
  param.ajouter("parcours_interface",&parcours_);
  param.ajouter("dt_impression_bilan_indicatrice", &dt_impression_bilan_indicatrice_);
  param.ajouter("verbosite_surface_efficace_face", &verbosite_surface_efficace_face_);
  param.ajouter("verbosite_surface_efficace_interface", &verbosite_surface_efficace_interface_);
  param.ajouter("seuil_indicatrice_negligeable", &seuil_indicatrice_negligeable_);

  param.ajouter_flag("use_barycentres_velocity", &use_barycentres_velocity_);
  param.ajouter_flag("read_barycentres_velocity", &read_barycentres_velocity_);
  param.ajouter("maillage_ft_ijk", &maillage_ft_ijk_);


  // param.ajouter_non_std("terme_gravite",(this));
  param.ajouter("terme_gravite", &terme_gravite_); // XD_ADD_P chaine(into=["rho_g","grad_i"]) not_set
  param.dictionnaire("rho_g", GRAVITE_RHO_G);
  param.dictionnaire("grad_i", GRAVITE_GRAD_I);

  param.lire_avec_accolades(is);

  Cout << "IJK_Interfaces::readOn : Option gravite : " << terme_gravite_
       << " { " << (int)GRAVITE_RHO_G << " : GRAVITE_RHO_G, "
       << (int)GRAVITE_GRAD_I << " : GRAVITE_GRAD_I} "
       << finl;

  //    Cout << "IJK_Interfaces::readOn : Les options lues sont : " << finl;
  //    param.print(Cout);

  if (compo_to_group_.size_array() != 0)
    {
      nb_groups_ = max_array(compo_to_group_) + 1;
    }
  Cout << "IJK_Interfaces::readOn : il y a : " << nb_groups_
       << " classes/groupes de bulles suivis " << finl;

  if (positions_reference_.size_array() != 0)
    {
      flag_positions_reference_ = 1;
      Cout << "IJK_Interfaces::readOn : " << positions_reference_.dimension(0)
           << " bulles reelles reprises avec une position de reference. " << finl;
    }
  if (mean_force_.size_array() != 0)
    {
      Cout << "IJK_Interfaces::readOn : " << mean_force_.dimension(0)
           << " bulles reelles reprises avec une force_moyenne imposee. " << finl;
    }
  return is;
}

int IJK_Interfaces::lire_motcle_non_standard(const Motcle& un_mot, Entree& is)
{
  if (un_mot=="parcours_interface")
    {
      is >> parcours_;
      return 1;
    }
  else
    {
      Cerr << "Unknown Keyword " << un_mot << " in IJK_Interfaces::lire_motcle_non_standard" << finl;
      Process::exit();
    }
  return 1;
}

void IJK_Interfaces::activate_cut_cell()
{
  cut_cell_activated_ = 1;
  indicatrice_surfacique_efficace_face_.associer_ephemere(*ref_ijk_ft_->get_cut_cell_disc());
  indicatrice_surfacique_efficace_face_initial_.associer_ephemere(*ref_ijk_ft_->get_cut_cell_disc());
  indicatrice_surfacique_efficace_face_correction_.associer_ephemere(*ref_ijk_ft_->get_cut_cell_disc());
  indicatrice_surfacique_efficace_face_absolute_error_.associer_ephemere(*ref_ijk_ft_->get_cut_cell_disc());
  surface_efficace_interface_.associer_ephemere(*ref_ijk_ft_->get_cut_cell_disc());
  surface_efficace_interface_initial_.associer_ephemere(*ref_ijk_ft_->get_cut_cell_disc());
  coord_deplacement_interface_.associer_ephemere(*ref_ijk_ft_->get_cut_cell_disc());
  vitesse_deplacement_interface_.associer_ephemere(*ref_ijk_ft_->get_cut_cell_disc());
  normale_deplacement_interface_.associer_ephemere(*ref_ijk_ft_->get_cut_cell_disc());

  if (!parcours_.get_correction_parcours_thomas())
    {
      Cerr << "Warning: The cut-cell method overrides default or user parameters to force the use of correction_parcours_thomas in parcours_interface." << finl;
      parcours_.set_correction_parcours_thomas();
    }

  if (!parcours_.get_parcours_sans_tolerance())
    {
      Cerr << "Warning: The cut-cell method overrides default or user parameters to force the use of parcours_sans_tolerance in parcours_interface." << finl;
      parcours_.set_parcours_sans_tolerance();
    }
}

void IJK_Interfaces::imprime_bilan_indicatrice()
{
  int count_total = 0;
  int count_40pct = 0;
  int count_30pct = 0;
  int count_20pct = 0;
  int count_10pct = 0;
  int count_5pct = 0;
  int count_1pct = 0;
  int count_0pct = 0;
  int count_pure = 0;

  double min_indicatrice = 1.;

  const Cut_cell_FT_Disc& cut_cell_disc = surface_efficace_interface_.get_cut_cell_disc();
  for (int n = 0; n < cut_cell_disc.get_n_loc(); n++)
    {
      Int3 ijk = cut_cell_disc.get_ijk(n);
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      count_total += 1;

      double indicatrice = std::min(indicatrice_ns_[next()](i,j,k), 1 - indicatrice_ns_[next()](i,j,k));
      if (indicatrice > 0.)
        min_indicatrice = std::min(min_indicatrice, indicatrice);

      if (indicatrice <= .5 && indicatrice > .4)
        count_40pct += 1;
      else if (indicatrice <= .4 && indicatrice > .3)
        count_30pct += 1;
      else if (indicatrice <= .3 && indicatrice > .2)
        count_20pct += 1;
      else if (indicatrice <= .2 && indicatrice > .1)
        count_10pct += 1;
      else if (indicatrice <= .1 && indicatrice > .05)
        count_5pct += 1;
      else if (indicatrice <= .05 && indicatrice > .01)
        count_1pct += 1;
      else if (indicatrice <= .01 && indicatrice > .0)
        count_0pct += 1;
      else
        {
          assert(indicatrice == 0.);
          count_pure += 1;
        }
    }
  count_total = Process::check_int_overflow(Process::mp_sum(count_total));
  count_40pct = Process::check_int_overflow(mp_sum(count_40pct));
  count_30pct = Process::check_int_overflow(mp_sum(count_30pct));
  count_20pct = Process::check_int_overflow(mp_sum(count_20pct));
  count_10pct = Process::check_int_overflow(mp_sum(count_10pct));
  count_5pct  = Process::check_int_overflow(mp_sum(count_5pct));
  count_1pct  = Process::check_int_overflow(mp_sum(count_1pct));
  count_0pct  = Process::check_int_overflow(mp_sum(count_0pct));
  count_pure  = Process::check_int_overflow(mp_sum(count_pure));

  min_indicatrice = Process::mp_min(min_indicatrice);

  Cerr << "Bilan de l'indicatrice:" << finl;
  Cerr << "  Valeur minimale: " << min_indicatrice << finl;
  Cerr << "  Compte sur les cellules de la structure diphasique:" << finl;
  Cerr << "    40-50%    " << count_40pct << finl;
  Cerr << "    30-40%    " << count_30pct << finl;
  Cerr << "    20-30%    " << count_20pct << finl;
  Cerr << "    10-20%    " << count_10pct << finl;
  Cerr << "     5-10%    " << count_5pct  << finl;
  Cerr << "     1-5%     " << count_1pct  << finl;
  Cerr << "     0-1%     " << count_0pct  << finl;
  Cerr << "     pure     " << count_pure  << finl;
  Cerr << "    ------    " << finl;
  Cerr << "     total    " << count_total  << finl;
}

void IJK_Interfaces::calcul_surface_efficace_face(TYPE_SURFACE_EFFICACE_FACE type_surface_efficace_face, double timestep, const Cut_field_vector3_double& total_velocity)
{
  // Attention : suppose que indicatrice_surfacique_efficace_face_/initial_ est correctement initialisee.
  // (Appel prealable a calcul_surface_efficace_face_initial())

  int iteration_solver_surface_efficace_face = 0;

  if (type_surface_efficace_face == TYPE_SURFACE_EFFICACE_FACE::CONSERVATION_VOLUME_ITERATIF)
    {
      // Raffinement de l'estimation initiale la surface efficace des faces
      Cut_cell_surface_efficace::calcul_surface_face_efficace_iteratif(
        verbosite_surface_efficace_face_,
        timestep,
        total_velocity,
        iteration_solver_surface_efficace_face,
        indicatrice_ns_[old()],
        indicatrice_ns_[next()],
        indicatrice_surfacique_face_ns_[old()],
        indicatrice_surfacique_face_ns_[next()],
        indicatrice_surfacique_efficace_face_,
        indicatrice_surfacique_efficace_face_initial_,
        indicatrice_surfacique_efficace_face_correction_,
        indicatrice_surfacique_efficace_face_absolute_error_);
    }
  Cut_cell_surface_efficace::imprimer_informations_surface_efficace_face(
    verbosite_surface_efficace_face_,
    iteration_solver_surface_efficace_face,
    timestep,
    total_velocity,
    indicatrice_ns_[old()],
    indicatrice_ns_[next()],
    indicatrice_surfacique_efficace_face_,
    indicatrice_surfacique_efficace_face_initial_);
}

void IJK_Interfaces::calcul_surface_efficace_interface(TYPE_SURFACE_EFFICACE_INTERFACE type_surface_efficace_interface, double timestep, const Cut_field_vector3_double& velocity)
{
  // Attention : suppose que indicatrice_surfacique_efficace_interface_/initial_ est correctement initialisee.
  // (Appel prealable a calcul_surface_efficace_interface_initial())

  Cut_cell_surface_efficace::calcul_surface_interface_efficace_initiale(
    (type_surface_efficace_interface == TYPE_SURFACE_EFFICACE_INTERFACE::EXPLICITE),
    indicatrice_ns_[old()],
    indicatrice_ns_[next()],
    surface_interface_ns_[old()],
    surface_interface_ns_[next()],
    normal_of_interf_ns_[old()],
    normal_of_interf_ns_[next()],
    normale_deplacement_interface_,
    surface_efficace_interface_,
    surface_efficace_interface_initial_);

  Cut_cell_surface_efficace::calcul_vitesse_interface(
    velocity,
    indicatrice_ns_[old()],
    indicatrice_ns_[next()],
    barycentre_phase1_ns_[old()],
    barycentre_phase1_ns_[next()],
    coord_deplacement_interface_,
    vitesse_deplacement_interface_);

  if (type_surface_efficace_interface == TYPE_SURFACE_EFFICACE_INTERFACE::CONSERVATION_VOLUME)
    {
      // Raffinement de l'estimation initiale la surface efficace de l'interface
      Cut_cell_surface_efficace::calcul_surface_interface_efficace(
        timestep,
        velocity,
        indicatrice_ns_[old()],
        indicatrice_ns_[next()],
        vitesse_deplacement_interface_,
        normale_deplacement_interface_,
        surface_efficace_interface_);
    }
  Cut_cell_surface_efficace::imprimer_informations_surface_efficace_interface(
    verbosite_surface_efficace_interface_,
    timestep,
    velocity,
    indicatrice_ns_[old()],
    indicatrice_ns_[next()],
    surface_efficace_interface_,
    surface_efficace_interface_initial_,
    normale_deplacement_interface_,
    vitesse_deplacement_interface_);
}

void IJK_Interfaces::calcul_surface_efficace_face_initial(TYPE_SURFACE_EFFICACE_FACE type_surface_efficace_face)
{
  Cut_cell_surface_efficace::calcul_surface_face_efficace_initiale(
    (type_surface_efficace_face == TYPE_SURFACE_EFFICACE_FACE::EXPLICITE),
    indicatrice_surfacique_face_ns_[old()],
    indicatrice_surfacique_face_ns_[next()],
    indicatrice_surfacique_efficace_face_,
    indicatrice_surfacique_efficace_face_initial_);
}

void IJK_Interfaces::calcul_surface_efficace_interface_initial(TYPE_SURFACE_EFFICACE_INTERFACE type_surface_efficace_interface)
{
  Cut_cell_surface_efficace::calcul_surface_interface_efficace_initiale(
    (type_surface_efficace_interface == TYPE_SURFACE_EFFICACE_INTERFACE::EXPLICITE),
    indicatrice_ns_[old()],
    indicatrice_ns_[next()],
    surface_interface_ns_[old()],
    surface_interface_ns_[next()],
    normal_of_interf_ns_[old()],
    normal_of_interf_ns_[next()],
    normale_deplacement_interface_,
    surface_efficace_interface_,
    surface_efficace_interface_initial_);
}

void IJK_Interfaces::compute_vinterp()
{
  const IJK_Field_vector3_double& velocity_ft = ref_ijk_ft_->eq_ns().get_velocity_ft();
  Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const DoubleTab& sommets = mesh.sommets(); // Tableau des coordonnees des marqueurs.
  int nbsom = sommets.dimension(0);
  ArrOfDouble vinterp_component(nbsom);
  vinterp_.resize(nbsom, 3);
  for (int direction = 0; direction < 3; direction++)
    {
      // Interpolate the "field" at the requested "coordinates" (array with 3
      // columns), and stores into "result" void ijk_interpolate(const
      // IJK_Field_double & field, const DoubleTab &coordinates, ArrOfDouble &
      // result)
      ijk_interpolate_skip_unknown_points(velocity_ft[direction], sommets, vinterp_component,
                                          1.e5 /* value for unknown points */);
      for (int i = 0; i < nbsom; i++)
        vinterp_(i, direction) = vinterp_component[i];
    }
}

Probleme_FTD_IJK_base& IJK_Interfaces::probleme_ijk()
{
  return ref_cast(Probleme_FTD_IJK_base, mon_probleme.valeur());
}

const Probleme_FTD_IJK_base& IJK_Interfaces::probleme_ijk() const
{
  return ref_cast(Probleme_FTD_IJK_base, mon_probleme.valeur());
}

const IJK_Field_double& IJK_Interfaces::get_IJK_field(const Motcle& nom)
{
  if (nom=="COURBURE")
    {
      // TODO ABN may use ref_array to avoid copy ??
      IJK_Field_double& courb = scalar_post_fields_.at("COURBURE");
      courb.data() = maillage_ft_ijk_.get_update_courbure_sommets();
    }
  if (nom=="CONCENTRATION_INTERFACE")
    {
      IJK_Field_double& conc_interf = scalar_post_fields_.at("CONCENTRATION_INTERFACE");
      conc_interf.data() = maillage_ft_ijk_.Surfactant_facettes().get_FT_field_Array();
    }
  if (nom=="GRADX_CONCENTRATION_INTERFACE")
    {
      IJK_Field_double& dconc_interf_x = scalar_post_fields_.at("GRADX_CONCENTRATION_INTERFACE");
      dconc_interf_x.data() = maillage_ft_ijk_.Surfactant_facettes().get_Grad_FT_field_Array(0);
    }
  if (nom=="GRADY_CONCENTRATION_INTERFACE")
    {
      IJK_Field_double& dconc_interf_y = scalar_post_fields_.at("GRADY_CONCENTRATION_INTERFACE");
      dconc_interf_y.data() = maillage_ft_ijk_.Surfactant_facettes().get_Grad_FT_field_Array(1);
    }
  if (nom=="GRADZ_CONCENTRATION_INTERFACE")
    {
      IJK_Field_double& dconc_interf_z = scalar_post_fields_.at("GRADZ_CONCENTRATION_INTERFACE");
      dconc_interf_z.data() = maillage_ft_ijk_.Surfactant_facettes().get_Grad_FT_field_Array(2);
    }
  if (nom=="LAPLACIAN_CONCENTRATION_INTERFACE")
    {
      IJK_Field_double& lapla_interf = scalar_post_fields_.at("LAPLACIAN_CONCENTRATION_INTERFACE");
      lapla_interf.data() = maillage_ft_ijk_.Surfactant_facettes().get_Laplacian_FT_field_Array();
    }
  if (nom=="SIGMA")
    {
      IJK_Field_double& sigma = scalar_post_fields_.at("SIGMA");
      sigma.data() = maillage_ft_ijk_.Surfactant_facettes().get_sigma_sommets();
    }
  if (nom=="GRADX_SIGMA")
    {
      IJK_Field_double& dsigma_x = scalar_post_fields_.at("GRADX_SIGMA");
      dsigma_x.data() = maillage_ft_ijk_.Surfactant_facettes().get_interfacial_source_term_sommets(0);
    }
  if (nom=="GRADY_SIGMA")
    {
      IJK_Field_double& dsigma_y = scalar_post_fields_.at("GRADY_SIGMA");
      dsigma_y.data() = maillage_ft_ijk_.Surfactant_facettes().get_interfacial_source_term_sommets(1);
    }
  if (nom=="GRADZ_SIGMA")
    {
      IJK_Field_double& dsigma_z = scalar_post_fields_.at("GRADZ_SIGMA");
      dsigma_z.data() = maillage_ft_ijk_.Surfactant_facettes().get_interfacial_source_term_sommets(2);
    }
  if (nom=="DISTANCE_AUTRES_INTERFACES")
    {
      DoubleTab vr_to_closer; // The velocity of the closest neighbour
      IJK_Field_double& d = scalar_post_fields_.at("DISTANCE_AUTRES_INTERFACES");
      calculer_distance_autres_compo_connexe2(d.data(), vr_to_closer);
    }

  if(has_champ(nom))
    return champs_compris_.get_champ(nom);

  Cerr << "ERROR in IJK_Interfaces::get_IJK_field_vector : " << finl;
  Cerr << "Requested field '" << nom << "' is not recognized by IJK_Interfaces::get_IJK_field_vector()." << finl;
  throw;
}

const IJK_Field_vector3_double& IJK_Interfaces::get_IJK_field_vector(const Motcle& nom)
{
  if(has_champ_vectoriel(nom))
    return champs_compris_.get_champ_vectoriel(nom);

  Cerr << "ERROR in IJK_Interfaces::get_IJK_field_vector : " << finl;
  Cerr << "Requested field '" << nom << "' is not recognized by IJK_Interfaces::get_IJK_field_vector()." << finl;
  throw;
}

void IJK_Interfaces::Fill_postprocessable_fields(std::vector<FieldInfo_t>& chps)
{
  std::vector<FieldInfo_t> c =
  {
    // Name     /     Localisation (elem, face, ...) /    Nature (scalare, vector)   / Located on interface?
    { "INDICATRICE", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "INDICATRICE_FT", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "OLD_INDICATRICE", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "OLD_INDICATRICE_FT", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "REPULSION_FT", Entity::ELEMENT, Nature_du_champ::scalaire, false },
    { "CUT_FIELDS_BARY_L", Entity::ELEMENT, Nature_du_champ::vectoriel, false },

    { "COURBURE", Entity::NODE, Nature_du_champ::scalaire, true },
    { "CONCENTRATION_INTERFACE", Entity::ELEMENT, Nature_du_champ::scalaire, true },
    { "GRADX_CONCENTRATION_INTERFACE", Entity::NODE, Nature_du_champ::scalaire, true },
    { "GRADY_CONCENTRATION_INTERFACE", Entity::NODE, Nature_du_champ::scalaire, true },
    { "GRADZ_CONCENTRATION_INTERFACE", Entity::NODE, Nature_du_champ::scalaire, true },
    { "SIGMA", Entity::NODE, Nature_du_champ::scalaire, true },
    { "GRADX_SIGMA", Entity::NODE, Nature_du_champ::scalaire, true },
    { "GRADY_SIGMA", Entity::NODE, Nature_du_champ::scalaire, true },
    { "GRADZ_SIGMA", Entity::NODE, Nature_du_champ::scalaire, true },
    { "LAPLACIAN_CONCENTRATION_INTERFACE", Entity::ELEMENT, Nature_du_champ::scalaire, true },
    { "DISTANCE_AUTRES_INTERFACES", Entity::NODE, Nature_du_champ::scalaire, true }
  };

  chps.insert(chps.end(), c.begin(), c.end());
}

void IJK_Interfaces::get_noms_champs_postraitables(Noms& noms,Option opt) const
{
  for (const auto& n : champs_compris_.liste_noms_compris())
    noms.add(n);
  for (const auto& n : champs_compris_.liste_noms_compris_vectoriel())
    noms.add(n);
}

void IJK_Interfaces::update_indicatrice_variables_monofluides()
{
  Cerr << "Reset bubble rising velocity calculations" << finl;
  reset_flags_and_counters();
  Cerr << "Compute compo_connex from bounding box" << finl;
  compute_compo_connex_from_bounding_box();
  Cerr << "Compute compo_connex from interface compo in mixed cells" << finl;
  compute_compo_connex_from_interface();
  Cerr << "Compute rising velocity from compo connex (barycentre calc)" << finl;
  compute_rising_velocities_from_compo();
}

void IJK_Interfaces::initialize(const Domaine_IJK& domaine_FT,
                                const Domaine_IJK& domaine_NS,
                                const Domaine_dis_base& domaine_dis,
                                const int thermal_probes_ghost_cells,
                                const bool compute_vint,
                                const bool is_switch)
{
  Cerr << "Entree dans IJK_Interfaces::initialize" << finl;
  set_recompute_indicator(CLASSIC_METHOD);

  ref_domaine_ = domaine_FT;

  surface_vapeur_par_face_computation_.initialize(domaine_FT);
  val_par_compo_in_cell_computation_.initialize(domaine_FT, maillage_ft_ijk_);

  field_repulsion_.allocate(ref_domaine_, Domaine_IJK::ELEM, 0, "REPULSION_FT");
  champs_compris_.ajoute_champ(field_repulsion_);

  if ((not is_switch) || cut_cell_activated_)
    {
      const int nb_ghost_cells = std::max(thermal_probes_ghost_cells, (int) 4);


      // IMPORTANT
      // fields using the old/next syntax must be given a name each
      // starting with OLD_ for the one at old(). Same name after that
      // then these pointers are switched in champ_compris_ when old and next are swapped
      // see method Champs_compris_IJK::switch_ft_fields()

      indicatrice_ft_[old()].allocate(domaine_FT, Domaine_IJK::ELEM, 2, "OLD_INDICATRICE_FT");
      indicatrice_ft_[old()].data() = 1.;
      indicatrice_ft_[old()].echange_espace_virtuel(indicatrice_ft_[old()].ghost());
      indicatrice_ft_[next()].allocate(domaine_FT, Domaine_IJK::ELEM, 2, "INDICATRICE_FT");
      indicatrice_ft_[next()].data() = 1.;
      indicatrice_ft_[next()].echange_espace_virtuel(indicatrice_ft_[next()].ghost());
      // Register INDICATRICE_FT:
      champs_compris_.ajoute_champ(indicatrice_ft_[old()]);
      champs_compris_.ajoute_champ(indicatrice_ft_[next()]);

      indicatrice_ns_[old()].allocate(domaine_NS, Domaine_IJK::ELEM, nb_ghost_cells, "OLD_INDICATRICE");
      indicatrice_ns_[old()].data() = 1.;
      allocate_cell_vector(groups_indicatrice_ns_[old()], domaine_NS, 1);
      allocate_cell_vector(groups_indicatrice_ns_[next()], domaine_NS, 1);
      indicatrice_ns_[old()].echange_espace_virtuel(indicatrice_ns_[old()].ghost());
      indicatrice_ns_[next()].allocate(domaine_NS, Domaine_IJK::ELEM, nb_ghost_cells, "INDICATRICE");
      indicatrice_ns_[next()].data() = 1.;
      indicatrice_ns_[next()].echange_espace_virtuel(indicatrice_ns_[next()].ghost());
      // Register INDICATRICE:
      champs_compris_.ajoute_champ(indicatrice_ns_[old()]);
      champs_compris_.ajoute_champ(indicatrice_ns_[next()]);

      indicatrice_avant_remaillage_ft_.allocate(domaine_FT, Domaine_IJK::ELEM, 2);
      indicatrice_avant_remaillage_ft_.data() = 1.;
      indicatrice_avant_remaillage_ft_.echange_espace_virtuel(indicatrice_avant_remaillage_ft_.ghost());
      indicatrice_avant_remaillage_ns_.allocate(domaine_NS, Domaine_IJK::ELEM, nb_ghost_cells);
      indicatrice_avant_remaillage_ns_.data() = 1.;
      indicatrice_avant_remaillage_ns_.echange_espace_virtuel(indicatrice_avant_remaillage_ns_.ghost());
      indicatrice_apres_remaillage_ft_.allocate(domaine_FT, Domaine_IJK::ELEM, 2);
      indicatrice_apres_remaillage_ft_.data() = 1.;
      indicatrice_apres_remaillage_ft_.echange_espace_virtuel(indicatrice_apres_remaillage_ft_.ghost());
      indicatrice_apres_remaillage_ns_.allocate(domaine_NS, Domaine_IJK::ELEM, nb_ghost_cells);
      indicatrice_apres_remaillage_ns_.data() = 1.;
      indicatrice_apres_remaillage_ns_.echange_espace_virtuel(indicatrice_apres_remaillage_ns_.ghost());
      delta_volume_theorique_bilan_ns_.allocate(domaine_NS, Domaine_IJK::ELEM, nb_ghost_cells);
      delta_volume_theorique_bilan_ns_.data() = 0.;
      delta_volume_theorique_bilan_ns_.echange_espace_virtuel(delta_volume_theorique_bilan_ns_.ghost());
      allocate_cell_vector(groups_indicatrice_ft_[old()], domaine_FT, 1);
      allocate_cell_vector(groups_indicatrice_ft_[next()], domaine_FT, 1);
#if VERIF_INDIC
      indicatrice_ft_test_.allocate(domaine_FT, Domaine_IJK::ELEM, 1);
      allocate_cell_vector(groups_indicatrice_ft_test_, domaine_FT, 1);
      allocate_cell_vector(groups_indicatrice_ft_test_, domaine_FT, 1);
#endif
      nb_compo_traversante_[old()].allocate(domaine_FT, Domaine_IJK::ELEM, 0);
      nb_compo_traversante_[next()].allocate(domaine_FT, Domaine_IJK::ELEM, 0);
      for (int i = 0; i < max_authorized_nb_of_components_; i++)
        {
          compos_traversantes_[old()][i].allocate(domaine_FT, Domaine_IJK::ELEM, 1);
          surface_par_compo_[old()][i].allocate(domaine_FT, Domaine_IJK::ELEM, 1);
          indicatrice_par_compo_[old()][i].allocate(domaine_FT, Domaine_IJK::ELEM, 1);
          courbure_par_compo_[old()][i].allocate(domaine_FT, Domaine_IJK::ELEM, 1);
          phi_par_compo_[old()][i].allocate(domaine_FT, Domaine_IJK::ELEM, 1);
          surf_par_compo_[old()][i].allocate(domaine_FT, Domaine_IJK::ELEM, 1);
          for (int dir = 0; dir < 3; dir++)
            source_interf_par_compo_[old()][dir][i].allocate(domaine_FT, Domaine_IJK::ELEM, 1);
          repuls_par_compo_[old()][i].allocate(domaine_FT, Domaine_IJK::ELEM, 1);
          compos_traversantes_[next()][i].allocate(domaine_FT, Domaine_IJK::ELEM, 1);
          surface_par_compo_[next()][i].allocate(domaine_FT, Domaine_IJK::ELEM, 1);
          indicatrice_par_compo_[next()][i].allocate(domaine_FT, Domaine_IJK::ELEM, 1);
          phi_par_compo_[next()][i].allocate(domaine_FT, Domaine_IJK::ELEM, 1);
          surf_par_compo_[next()][i].allocate(domaine_FT, Domaine_IJK::ELEM, 1);
          for (int dir = 0; dir < 3; dir++)
            source_interf_par_compo_[next()][dir][i].allocate(domaine_FT, Domaine_IJK::ELEM, 1);
          repuls_par_compo_[next()][i].allocate(domaine_FT, Domaine_IJK::ELEM, 1);
          courbure_par_compo_[next()][i].allocate(domaine_FT, Domaine_IJK::ELEM, 1);
          // Et pour les vecteurs :
          for (int dir = 0; dir < 3; dir++)
            {
              const int idx = i * 3 + dir;
              normale_par_compo_[old()][idx].allocate(domaine_FT, Domaine_IJK::ELEM, 2);
              bary_par_compo_[old()][idx].allocate(domaine_FT, Domaine_IJK::ELEM, 1);
              normale_par_compo_[next()][idx].allocate(domaine_FT, Domaine_IJK::ELEM, 2);
              bary_par_compo_[next()][idx].allocate(domaine_FT, Domaine_IJK::ELEM, 1);
            }
        }

      surface_interface_ft_[old()].allocate(domaine_FT, Domaine_IJK::ELEM, 2);
      surface_interface_ft_[next()].allocate(domaine_FT, Domaine_IJK::ELEM, 2);
      surface_interface_ns_[old()].allocate(domaine_NS, Domaine_IJK::ELEM, nb_ghost_cells);
      surface_interface_ns_[next()].allocate(domaine_NS, Domaine_IJK::ELEM, nb_ghost_cells);

      surface_interface_ft_[old()].data() = 0.;
      surface_interface_ft_[next()].data() = 0.;
      surface_interface_ns_[old()].data() = 0.;
      surface_interface_ns_[next()].data() = 0.;

      allocate_cell_vector(barycentre_phase1_ft_[old()], domaine_FT, 2);
      allocate_cell_vector(barycentre_phase1_ft_[next()], domaine_FT, 2);
      allocate_cell_vector(barycentre_phase1_ns_[old()], domaine_NS, nb_ghost_cells, "OLD_CUT_FIELDS_BARY_L");
      allocate_cell_vector(barycentre_phase1_ns_[next()], domaine_NS, nb_ghost_cells, "CUT_FIELDS_BARY_L");
      champs_compris_.ajoute_champ_vectoriel(barycentre_phase1_ns_[old()]);
      champs_compris_.ajoute_champ_vectoriel(barycentre_phase1_ns_[next()]);

      for (int bary_compo = 0; bary_compo < 3; bary_compo++)
        {
          barycentre_phase1_ft_[old()][bary_compo].data() = 0.;
          barycentre_phase1_ft_[next()][bary_compo].data() = 0.;
          barycentre_phase1_ns_[old()][bary_compo].data() = 0.;
          barycentre_phase1_ns_[next()][bary_compo].data() = 0.;
        }

      allocate_velocity(indicatrice_surfacique_face_ft_[old()], domaine_FT, 2);
      allocate_velocity(indicatrice_surfacique_face_ft_[next()], domaine_FT, 2);
      allocate_velocity(indicatrice_surfacique_face_ns_[old()], domaine_NS, nb_ghost_cells);
      allocate_velocity(indicatrice_surfacique_face_ns_[next()], domaine_NS, nb_ghost_cells);

      for (int face_dir = 0; face_dir < 3; face_dir++)
        {
          indicatrice_surfacique_face_ft_[old()][face_dir].data() = 0.;
          indicatrice_surfacique_face_ft_[next()][face_dir].data() = 0.;
          indicatrice_surfacique_face_ns_[old()][face_dir].data() = 0.;
          indicatrice_surfacique_face_ns_[next()][face_dir].data() = 0.;
        }

      allocate_velocity(indicatrice_surfacique_avant_remaillage_face_ft_, domaine_FT, 2);
      allocate_velocity(indicatrice_surfacique_avant_remaillage_face_ns_, domaine_NS, nb_ghost_cells);
      allocate_velocity(indicatrice_surfacique_apres_remaillage_face_ft_, domaine_FT, 2);
      allocate_velocity(indicatrice_surfacique_apres_remaillage_face_ns_, domaine_NS, nb_ghost_cells);

      for (int face_dir = 0; face_dir < 3; face_dir++)
        {
          indicatrice_surfacique_avant_remaillage_face_ft_[face_dir].data() = 0.;
          indicatrice_surfacique_avant_remaillage_face_ns_[face_dir].data() = 0.;
          indicatrice_surfacique_apres_remaillage_face_ft_[face_dir].data() = 0.;
          indicatrice_surfacique_apres_remaillage_face_ns_[face_dir].data() = 0.;
        }

      for (int face_dir = 0; face_dir < 3; face_dir++)
        {
          for (int bary_compo = 0; bary_compo < 2; bary_compo++)
            {
              barycentre_phase1_face_ft_[old()][face_dir][bary_compo].allocate(domaine_FT, Domaine_IJK::ELEM, 2);
              barycentre_phase1_face_ft_[next()][face_dir][bary_compo].allocate(domaine_FT, Domaine_IJK::ELEM, 2);
              barycentre_phase1_face_ns_[old()][face_dir][bary_compo].allocate(domaine_NS, Domaine_IJK::ELEM, nb_ghost_cells);
              barycentre_phase1_face_ns_[next()][face_dir][bary_compo].allocate(domaine_NS, Domaine_IJK::ELEM, nb_ghost_cells);
            }
        }

      for (int face_dir = 0; face_dir < 3; face_dir++)
        {
          for (int bary_compo = 0; bary_compo < 2; bary_compo++)
            {
              barycentre_phase1_face_ft_[old()][face_dir][bary_compo].data() = 0.;
              barycentre_phase1_face_ft_[next()][face_dir][bary_compo].data() = 0.;
              barycentre_phase1_face_ns_[old()][face_dir][bary_compo].data() = 0.;
              barycentre_phase1_face_ns_[next()][face_dir][bary_compo].data() = 0.;
            }
        }

      allocate_velocity(normal_of_interf_[old()], domaine_FT, 2);
      allocate_velocity(normal_of_interf_[next()], domaine_FT, 2);
      allocate_velocity(normal_of_interf_ns_[old()], domaine_NS, nb_ghost_cells);
      allocate_velocity(normal_of_interf_ns_[next()], domaine_NS, nb_ghost_cells);

      allocate_velocity(bary_of_interf_[old()], domaine_FT, 1);
      allocate_velocity(bary_of_interf_[next()], domaine_FT, 1);
      allocate_velocity(bary_of_interf_ns_[old()], domaine_NS, 1);
      allocate_velocity(bary_of_interf_ns_[next()], domaine_NS, 1);

      allocate_velocity(surface_vapeur_par_face_[old()], domaine_FT, 1);
      allocate_velocity(surface_vapeur_par_face_[next()], domaine_FT, 1);
      allocate_velocity(surface_vapeur_par_face_ns_[old()], domaine_NS, 1);
      allocate_velocity(surface_vapeur_par_face_ns_[next()], domaine_NS, 1);

      for (int d = 0; d < 3; d++)
        {
          surface_vapeur_par_face_[old()][d].data() = 0.;
          surface_vapeur_par_face_[next()][d].data() = 0.;
          allocate_velocity(barycentre_vapeur_par_face_[old()][d], domaine_FT, 1);
          allocate_velocity(barycentre_vapeur_par_face_[next()][d], domaine_FT, 1);
          surface_vapeur_par_face_ns_[old()][d].data() = 0.;
          surface_vapeur_par_face_ns_[next()][d].data() = 0.;
          allocate_velocity(barycentre_vapeur_par_face_ns_[old()][d], domaine_NS, 1);
          allocate_velocity(barycentre_vapeur_par_face_ns_[next()][d], domaine_NS, 1);
        }

      if (cut_cell_activated_)
        {
          Cut_field_vector3_double& cut_field_deformation_velocity = static_cast<Cut_field_vector3_double&>(deformation_velocity_);
          allocate_velocity_ephemere(*ref_ijk_ft_->get_cut_cell_disc(), cut_field_deformation_velocity, domaine_NS, 2);
        }

      allocate_velocity(indicatrice_surfacique_efficace_deformation_face_, domaine_NS, nb_ghost_cells);
    }

  if (Option_IJK::DISABLE_DIPHASIQUE)
    return;

  refdomaine_dis_ = domaine_dis;

  // Calcul de la bounding box de Navier Stokes et stockage en memoire de la
  // perio.
  bounding_box_NS_domain_.resize(3, 2);
  // Calcul du domaine au-dela duquel une bulle doit etre repliquee.
  // Cette domaine est une peu plus petite que le domaine NS_domain, si la bulle
  // depasse de cette domaine elle est dupliquee.
  bounding_box_duplicate_criteria_.resize(3, 2);
  // Calcul du domaine au-dela duquel une bulle est detruite et remplacee
  // par son duplicata a l'autre extremite du domaine.
  // Cette domaine est un peu plus petite que le domaine etendu ou evolue le
  // maillage FT.
  bounding_box_forbidden_criteria_.resize(3, 2);
  for (int direction = 0; direction < 3; direction++)
    {
      const double ori = domaine_NS.get_origin(direction);
      const double len = domaine_NS.get_domain_length(direction);
      bounding_box_NS_domain_(direction, 0) = ori;
      bounding_box_NS_domain_(direction, 1) = ori + len;
      bool perio = domaine_NS.get_periodic_flag(direction);
      perio_NS_[direction] = perio;
      if (perio && !avoid_duplicata_)
        {
          // Les bulles qui entre dans les ncells_forbidden_ dernieres mailles
          // doivent etre deplacees.
          const double oriFT = domaine_FT.get_origin(direction);
          const double lenFT = domaine_FT.get_domain_length(direction);
          const double delta = domaine_FT.get_constant_delta(direction);
          // largeur du domaine, au bord du domaine navier stokes, au sein de
          // laquelle on declanche la duplication des bulles (si une bulle depasse a
          // l'exterieur du domaine, on duplique) Cette domaine tient compte du stencil
          // des forces de tension superficielle et de repulsion
          //duCluz  : 2 cest mieux pour les echanges espaces virtuels pour la condition de shear-periodicite
          double duplicate_stencil_width ;

          if(IJK_Shear_Periodic_helpler::defilement_ == 1)
            duplicate_stencil_width =
              std::max(2 * delta,
                       portee_force_repulsion_);
          else
            duplicate_stencil_width =
              std::max(delta,
                       portee_force_repulsion_);

          // GB2020.12.20 : avant c'etait 2. Est-ce
          // que la precaution etait necessaire? Elle
          // conduit a de plus gros cas tests comme
          // interfacial_temperature_and_flux
          bounding_box_duplicate_criteria_(direction, 0) = ori + duplicate_stencil_width;
          bounding_box_duplicate_criteria_(direction, 1) = ori + len - duplicate_stencil_width;
          bounding_box_forbidden_criteria_(direction, 0) = oriFT + ncells_forbidden_ * delta;
          bounding_box_forbidden_criteria_(direction, 1) = oriFT + lenFT - ncells_forbidden_ * delta;
        }
      else
        {
          if (avoid_duplicata_)
            {
              bounding_box_duplicate_criteria_(direction, 0) = ori - factor_length_duplicata_ * len;
              bounding_box_duplicate_criteria_(direction, 1) = ori + (factor_length_duplicata_ + 1.)*len;
              // Les bulles ne sont pas limitees au domaine NS :
              bounding_box_forbidden_criteria_(direction, 0) = bounding_box_duplicate_criteria_(direction, 0);
              bounding_box_forbidden_criteria_(direction, 1) = bounding_box_duplicate_criteria_(direction, 1);
            }
          else
            {
              bounding_box_duplicate_criteria_(direction, 0) = bounding_box_NS_domain_(direction, 0);
              bounding_box_duplicate_criteria_(direction, 1) = bounding_box_NS_domain_(direction, 1);
              // Les bulles sont limitees au domaine NS :
              bounding_box_forbidden_criteria_(direction, 0) = bounding_box_NS_domain_(direction, 0);
              bounding_box_forbidden_criteria_(direction, 1) = bounding_box_NS_domain_(direction, 1);
            }
        }
    }

  parcours_.associer_domaine_dis(domaine_dis);
  Domaine_VF& domaine_vf = ref_cast_non_const(Domaine_VF, domaine_dis);
  connectivite_frontieres_.associer_domaine_vf(domaine_vf);
  parcours_.associer_connectivite_frontieres(connectivite_frontieres_);

  maillage_ft_ijk_.initialize(domaine_FT, domaine_dis, parcours_, use_tryggvason_interfacial_source_);
  if (cut_cell_activated_)
    {
      old_maillage_ft_ijk_.initialize(domaine_FT, domaine_dis, parcours_);
    }

  // Lecture du maillage initial:
  maillage_ft_ijk_.lire_maillage_ft_dans_lata(fichier_reprise_interface_,
                                              timestep_reprise_interface_,
                                              lata_interfaces_meshname_);

  remaillage_ft_ijk_.associer_domaine(domaine_dis);

  // Si des bulles ghost sont lues, on les detruit
  // pour etre sur de les creer toutes.
  // Le compteur nb_bulles_ghost_ est remis a zero.
  supprimer_duplicata_bulles();

  // Initialisation du nombre de bulles :
  int loc_max;
  if (maillage_ft_ijk_.compo_connexe_facettes().size_array())
    loc_max = max_array(maillage_ft_ijk_.compo_connexe_facettes());
  else
    // on n'a pas une seule facette dans cette partie du decoupage :
    loc_max = -1; // Invalid value... Comme ca si on a -1 partout, on a 0 bulle.
  // ::mp_max force l'appel a la methode hors classe (version pour les ints)
  nb_bulles_reelles_ = Process::mp_max(loc_max) + 1; // car les composantes connexes commencent a 0.
  // On est en mesure de redimensionner le tableau (uniquement s'il n'est pas lu
  // en param:
  if (compo_to_group_.size_array() == 0)
    {
      compo_to_group_.resize(nb_bulles_reelles_);
      compo_to_group_ = 0;
    }
  assert(compo_to_group_.size_array() == nb_bulles_reelles_);

  // Si on commence de suivre les couleurs mais qu'on ne les reprends pas d'un
  // calcul precedent, il faut dimensionner le tableau et l'initaliser a 0 :
  if ((follow_colors_) && (through_yminus_.size_array() == 0))
    through_yminus_.resize(nb_bulles_reelles_);

  DoubleTab centre_gravite;
  ArrOfDouble volume_bulles;
  calculer_volume_bulles(volume_bulles, centre_gravite);
  maillage_ft_ijk_.Surfactant_facettes_non_const().initialize(maillage_ft_ijk_, centre_gravite);

  transferer_bulle_perio();
  // Si des interfaces sont lues en dehors du domaine de NS, il faut re-creer
  // leur ghost.
  creer_duplicata_bulles();

  // Initialisation du tableau de stockage des vitesses aux sommets en RK3 :
  RK3_G_store_vi_.resize(0, 3);
  distance_autres_interfaces_.resize_array(0);

// Maybe needed to post-pro initial condition :
  if (compute_vint)
    {
      if (compute_distance_autres_interfaces_ || (delta_p_max_repulsion_ > 0. && portee_force_repulsion_ > 0.))
        {
          DoubleTab vr_to_closer; // The velocity of the closest neighbour
          // Also calls to compute_vinterp.
          calculer_distance_autres_compo_connexe2(distance_autres_interfaces_, vr_to_closer);
        }
      else
        compute_vinterp();
    }
  // Mise en place des compteurs :

  force_time_n_.resize(nb_bulles_reelles_, 3);
  force_time_n_ = 0.;
  if (mean_force_.size_array() == 0)
    {
      mean_force_.resize(nb_bulles_reelles_, 3);
      mean_force_ = 0.;
    }

  nb_compo_in_num_compo_ = 0; // initially, waiting for the update when table is completed.
  if ((parser_ == 0) && (recompute_indicator_ == 0))
    {
      Cerr << "Error in option combination: invalid choice of optimized methods "
           << "for both the color_function in the forces computation "
           << "(parser_=0) and the indicator function calculation "
           << "(recompute_indicator_=0)"
           << finl;
      Cerr << "Maybe you are not using the forces computation?" << finl;
      Process::exit();
    }
  else if ((recompute_indicator_ == 1) || (parser_ == 0))
    {
      // Pour la methode historique (non-optim) de calcul de l'indicatrice ou pour
      // le calcul optimise de la force de rappel. Cree un tableau parallele
      // structure comme un tableau aux elements du maillage vdf, initialise a
      // zero.
      const Domaine& domaine = domaine_vf.domaine();
      domaine.creer_tableau_elements(num_compo_);
    }

  read_bubbles_barycentres_old_new(fichier_reprise_interface_);
  intersection_ijk_cell_.initialize(domaine_NS, *this);
  intersection_ijk_face_.initialize(domaine_NS, *this);
  ijk_compo_connex_.initialize(*this, is_switch);
}

void IJK_Interfaces::register_fields()
{
  // Register fields:
  scalar_post_fields_["COURBURE"] = IJK_Field_double();
  auto& courb = scalar_post_fields_.at("COURBURE");
  courb.nommer("COURBURE");
  champs_compris_.ajoute_champ(courb);

  scalar_post_fields_["CONCENTRATION_INTERFACE"] = IJK_Field_double();
  auto& conc_interf = scalar_post_fields_.at("CONCENTRATION_INTERFACE");
  conc_interf.nommer("CONCENTRATION_INTERFACE");
  champs_compris_.ajoute_champ(conc_interf);

  scalar_post_fields_["GRADX_CONCENTRATION_INTERFACE"] = IJK_Field_double();
  auto& dconc_interf_x = scalar_post_fields_.at("GRADX_CONCENTRATION_INTERFACE");
  dconc_interf_x.nommer("GRADX_CONCENTRATION_INTERFACE");
  champs_compris_.ajoute_champ(dconc_interf_x);

  scalar_post_fields_["GRADY_CONCENTRATION_INTERFACE"] = IJK_Field_double();
  auto& dconc_interf_y = scalar_post_fields_.at("GRADY_CONCENTRATION_INTERFACE");
  dconc_interf_y.nommer("GRADY_CONCENTRATION_INTERFACE");
  champs_compris_.ajoute_champ(dconc_interf_y);

  scalar_post_fields_["GRADZ_CONCENTRATION_INTERFACE"] = IJK_Field_double();
  auto& dconc_interf_z = scalar_post_fields_.at("GRADZ_CONCENTRATION_INTERFACE");
  dconc_interf_z.nommer("GRADZ_CONCENTRATION_INTERFACE");
  champs_compris_.ajoute_champ(dconc_interf_z);

  scalar_post_fields_["SIGMA"] = IJK_Field_double();
  auto& sigma = scalar_post_fields_.at("SIGMA");
  sigma.nommer("SIGMA");
  champs_compris_.ajoute_champ(sigma);

  scalar_post_fields_["GRADX_SIGMA"] = IJK_Field_double();
  auto& dsigma_x = scalar_post_fields_.at("GRADX_SIGMA");
  dsigma_x.nommer("GRADX_SIGMA");
  champs_compris_.ajoute_champ(dsigma_x);

  scalar_post_fields_["GRADY_SIGMA"] = IJK_Field_double();
  auto& dsigma_y = scalar_post_fields_.at("GRADY_SIGMA");
  dsigma_y.nommer("GRADY_SIGMA");
  champs_compris_.ajoute_champ(dsigma_y);

  scalar_post_fields_["GRADZ_SIGMA"] = IJK_Field_double();
  auto& dsigma_z = scalar_post_fields_.at("GRADZ_SIGMA");
  dsigma_z.nommer("GRADZ_SIGMA");
  champs_compris_.ajoute_champ(dsigma_z);

  scalar_post_fields_["LAPLACIAN_CONCENTRATION_INTERFACE"] = IJK_Field_double();
  auto& lapla_interf = scalar_post_fields_.at("LAPLACIAN_CONCENTRATION_INTERFACE");
  lapla_interf.nommer("LAPLACIAN_CONCENTRATION_INTERFACE");
  champs_compris_.ajoute_champ(lapla_interf);

  scalar_post_fields_["DISTANCE_AUTRES_INTERFACES"] = IJK_Field_double();
  auto& d = scalar_post_fields_.at("DISTANCE_AUTRES_INTERFACES");
  d.nommer("DISTANCE_AUTRES_INTERFACES");
  champs_compris_.ajoute_champ(d);
}

const Milieu_base& IJK_Interfaces::milieu() const
{
  Cerr << "IJK_Interfaces::milieu not coded ! access it from NS. " << finl;
  throw;
}

Milieu_base& IJK_Interfaces::milieu()
{
  Cerr << "IJK_Interfaces::milieu not coded ! access it from NS. " << finl;
  throw;
}

void IJK_Interfaces::associer_pb_base(const Probleme_base& pb)
{
  if (!sub_type(Probleme_FTD_IJK_base, pb))
    {
      Cerr << "Error for the method IJK_Interfaces::associer_pb_base\n";
      Cerr << " IJK_Interfaces equation must be associated to\n";
      Cerr << " a Probleme_FTD_IJK_base problem type\n";
      Process::exit();
    }
  mon_probleme = pb;
  if (nom_ == "??")
    {
      nom_ = pb.le_nom();
      nom_ += que_suis_je();
    }

  ref_ijk_ft_ = ref_cast(Probleme_FTD_IJK_base, pb);
  ijk_compo_connex_.associer(ref_ijk_ft_.valeur());
  // liste_post_instantanes_ = ijk_ft.post_.get_liste_post_instantanes();
}

void IJK_Interfaces::associer_switch(const Switch_FT_double& ijk_ft_switch)
{
  ref_ijk_ft_switch_ = ijk_ft_switch;
}

/*! Methode appelee lorsqu'on a mis "TOUS" dans la liste des champs a postraiter.
 * Elle ajoute a la liste tous les noms de champs postraitables par
 * IJK_Interfaces
 */
void IJK_Interfaces::posttraiter_tous_champs(Motcles& liste) const
{
  liste.add("INTERFACES");
  liste.add("COMPO_CONNEXE");
  liste.add("DISTANCE_AUTRES_INTERFACES");
  if (follow_colors_)
    liste.add("COLOR_Y");
  trustIdType size = RK3_G_store_vi_.size_array();
  size = Process::mp_sum(size);
  if (size)
    {
      liste.add("RK3_STORE_VI");
      liste.add("VI");
    }
  //  liste.add("NORMALE");
}

int IJK_Interfaces::posttraiter_champs_instantanes(const Motcles& liste_post_instantanes,
                                                   const char *lata_name,
                                                   const int lata_step) const
{
  throw; // THIS METHOD SHOULD NOT BE CALLED ANYMORE

  int n = 0; // nombre de champs postraites
  if (liste_post_instantanes.contient_("INTERFACES"))
    n++, dumplata_ft_mesh(lata_name, "INTERFACES", lata_step);
//  if (liste_post_instantanes.contient_("COMPO_CONNEXE"))
//    n++, dumplata_ft_field(lata_name, "INTERFACES", "COMPO_CONNEXE", "ELEM",  maillage_ft_ijk_.compo_connexe_facettes(), lata_step);
  if (liste_post_instantanes.contient_("COLOR_Y"))
    {
      if (!follow_colors_)
        {
          Cerr << "On demande le post-traitement du champ COLOR_Y alors que le flag follow_colors "
               << "n'a pas ete active... Ce champ est donc inconnu. Post-traitement refuse. " << finl;
          Process::exit();
        }
      ArrOfInt color;
      calculer_color(color);
      n++, dumplata_ft_field(lata_name, "INTERFACES", "COLOR_Y", "ELEM",  color, lata_step);
    }
  if (liste_post_instantanes.contient_("VINTERP"))
    {
      n++, dumplata_ft_field(lata_name, "INTERFACES", "VINTERP", "SOM",  RK3_G_store_vi_, lata_step);
    }
  if (liste_post_instantanes.contient_("VI"))
    {
      ArrOfDoubleFT vi;
      const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
      int nbsom = mesh.sommets().dimension(0);
      vi.resize_array(nbsom);
      // La methode est const! Un calcul de vinterp ici changerait l'objet interfaces_
      //compute_vinterp(); // Pour s'assurer que vinterp soit a jour avec le maillage (bon nb_som)
      //                      Ce n'est pas exactement celui utilise pour le transport puisqu'ici, la compo tangeante n'est pas retiree.
      for (int dir= 0; dir < 3; dir++)
        {
          if (vinterp_.dimension(0) == 0)
            {
              Cerr << "Posttraitement du champ VI vide : on stocke 0." << finl;
              vi = 0.;
            }
          else
            for (int i = 0; i < nbsom; i++)
              vi[i] = vinterp_(i,dir);

          if (dir==0)
            n++, dumplata_ft_field(lata_name, "INTERFACES", "VI_X", "SOM",  vi, lata_step);
          if (dir==1)
            n++, dumplata_ft_field(lata_name, "INTERFACES", "VI_Y", "SOM",  vi, lata_step);
          if (dir==2)
            n++, dumplata_ft_field(lata_name, "INTERFACES", "VI_Z", "SOM",  vi, lata_step);
        }

    }

  if (liste_post_instantanes.contient_("RK3_STORE_VI"))
    {
      ArrOfDoubleFT store_vi_compo;
      const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
      int nbsom = mesh.sommets().dimension(0);
      store_vi_compo.resize_array(nbsom);

      //    int size = RK3_G_store_vi_.size_array();
      //    size = Process::mp_sum(size);
      for (int dir= 0; dir < 3; dir++)
        {
          // Si le tableau n'est pas rempli, on retourne 0
          if (RK3_G_store_vi_.dimension(0) == 0)
            {
              Cerr << "Posttraitement du champ RK3_STORE_VI vide : on stocke 0." << finl;
              store_vi_compo = 0.;
            }
          else
            {
              for (int i = 0; i < nbsom; i++)
                {
                  const double val =  RK3_G_store_vi_(i,dir);
                  store_vi_compo[i] = val;
                }
            }
          if (dir==0)
            n++, dumplata_ft_field(lata_name, "INTERFACES", "RK3_STORE_VI_X", "SOM",  store_vi_compo, lata_step);
          if (dir==1)
            n++, dumplata_ft_field(lata_name, "INTERFACES", "RK3_STORE_VI_Y", "SOM",  store_vi_compo, lata_step);
          if (dir==2)
            n++, dumplata_ft_field(lata_name, "INTERFACES", "RK3_STORE_VI_Z", "SOM",  store_vi_compo, lata_step);
        }
    }
  return n;
}

// Suppression des bulles dans le domaine a eliminer pres des bords periodiques
// definie par ncells_deleted_ Dans ce cas, le maillage est modifie:
//  o Les ghosts sont supprimes de la sauvegarde.
//  o Les bulles qui sortent du domaine defini via ncells_deleted_ sont
//  supprimees aussi. o Le numero des autres bulles est re-attribue pour creer
//  une liste contigue (0, 1, 2, ...)
//    et supprimer les trous
void IJK_Interfaces::supprimer_certaines_bulles_reelles()
{
  Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  int flag = 0;
  if (Process::je_suis_maitre())
    {
      // Calcul du domaine dans lequelle une bulle est supprimee:
      bounding_box_delete_criteria_.resize(3, 2);
      const Domaine_IJK& geom_FT = ref_domaine_.valeur();
      for (int direction = 0; direction < 3; direction++)
        {
          if (perio_NS_[direction])
            {
              // Les bulles qui entre dans les ncells_forbidden_ dernieres mailles
              // doivent etre deplacees.
              const double oriFT = geom_FT.get_origin(direction);
              const double lenFT = geom_FT.get_domain_length(direction);
              const double delta = geom_FT.get_constant_delta(direction);
              bounding_box_delete_criteria_(direction, 0) = oriFT + ncells_deleted_ * delta;
              bounding_box_delete_criteria_(direction, 1) = oriFT + lenFT - ncells_deleted_ * delta;
            }
          else
            {
              // Les bulles sont limitees au domaine NS :
              bounding_box_delete_criteria_(direction, 0) = bounding_box_NS_domain_(direction, 0);
              bounding_box_delete_criteria_(direction, 1) = bounding_box_NS_domain_(direction, 1);
            }
        }
    }
  // broadcast to all mpi processes:
  envoyer_broadcast(bounding_box_delete_criteria_, 0);

  // Evaluation du cube contenant chaque bulle :
  DoubleTab bounding_box;
  calculer_bounding_box_bulles(bounding_box);

  // Pour chaque compo_connexe, remplir dans le tableau
  // masque_duplicata_pour_compo un encodage du deplacement maximal pour toutes
  // les bulles qui sortent de delete_criteria:
  ArrOfInt masque_delete_pour_compo;
  preparer_duplicata_bulles_masque_6bit(bounding_box, bounding_box_delete_criteria_, masque_delete_pour_compo);
  // Le masque reste a zero pour les bulles qui sont dans la
  // box_delete_criteria. Il faut donc supprimer les bulles en dehors, dont le
  // masque est different de 0.

  // Duplique et deplace les bulles de la liste :
  // dupliquer_bulle_perio(masque_delete_pour_compo);

  int nb_bulles_reelles_futur = nb_bulles_reelles_;
  int nb_bulles_reelles_deleted = 0;
  ArrOfInt old_to_new_compo;
  old_to_new_compo.resize_array(nb_bulles_reelles_);
  if (Process::je_suis_maitre())
    {
      int count = 0;
      for (int i = 0; i < masque_delete_pour_compo.size_array(); i++)
        {
          if (masque_delete_pour_compo[i] != 0)
            {
              flag = 1;
              nb_bulles_reelles_futur -= 1;
              nb_bulles_reelles_deleted += 1;
              old_to_new_compo[i] = -1;
            }
          else
            {
              old_to_new_compo[i] = count;
              count++;
            }
        }
      nb_bulles_reelles_futur = count;
    }
  // broadcast to all mpi processes:
  envoyer_broadcast(old_to_new_compo, 0);
  envoyer_broadcast(nb_bulles_reelles_futur, 0);
  envoyer_broadcast(flag, 0);
  envoyer_broadcast(nb_bulles_reelles_deleted, 0);

  if (flag)
    {
      Cerr << nb_bulles_reelles_deleted << " marked for deletion out of " << nb_bulles_reelles_
           << " so that only " << nb_bulles_reelles_futur << " will remain." << finl;
      const ArrOfInt& compo_connexe_facettes = mesh.compo_connexe_facettes();
      int icompo;
      // on remplace les compo a supprimer par -1
      for (int fa7 = 0; fa7 < mesh.nb_facettes(); fa7++)
        {
          icompo = compo_connexe_facettes[fa7];
          assert(icompo >= 0); // les duplicatas ne sont pas la.
          // imasque = masque_delete_pour_compo[icompo];
          // On lui donne son nouveau numero ou on met '-1' pour la marquer pour
          // suppression
          mesh.set_composante_connexe(fa7, old_to_new_compo[icompo]);
        }

      // On supprime les bulles que l'on vient de renumeroter -1:
      // (les duplicatas n'etaient pas presents).
      supprimer_duplicata_bulles();
      Cerr << "The number of bubbles has been reduced from " << nb_bulles_reelles_
           << " to " << nb_bulles_reelles_futur << finl;

      recompute_indicator_ = 1;

      if (Process::je_suis_maitre())
        {
          for (icompo = 0; icompo < nb_bulles_reelles_; icompo++)
            {
              int inew = old_to_new_compo[icompo];
              assert(inew <= icompo);
              // On ne fait quelque chose que pour les bulles conservees (pas marquees
              // -1)
              if (inew > -1)
                {
                  compo_to_group_[inew] = compo_to_group_[icompo];
                  if (through_yminus_.size_array())
                    through_yminus_[inew] = through_yminus_[icompo];
                  if (positions_reference_.size_array())
                    positions_reference_[inew] = positions_reference_[icompo];
                }
            }
          compo_to_group_.resize_array(nb_bulles_reelles_futur);
          if (through_yminus_.size_array())
            through_yminus_.resize_array(nb_bulles_reelles_futur);
          if (positions_reference_.size_array())
            positions_reference_.resize_array(nb_bulles_reelles_futur);
          Cerr << "The table of bubbles groups (compo_to_group_), "
               << " (as well as possibly through_yminus_ and positions_reference_ "
               << "if needed)"  << " have been updated in accordance." << finl;
        }
      envoyer_broadcast(compo_to_group_, 0);
      envoyer_broadcast(through_yminus_, 0);
      envoyer_broadcast(positions_reference_, 0);
      nb_bulles_reelles_ = nb_bulles_reelles_futur;
    }
}

// Attention a l'usage du mot cle 'ncells_deleted_' qui conduit a la suppression
// de bulles (reelles) et pas uniquement des ghosts.
void IJK_Interfaces::sauvegarder_interfaces(const char *lata_name, const Nom& interf_name) // const
{
  fichier_sauvegarde_interface_ = lata_name;
  timestep_sauvegarde_interface_ = 1;
  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;

  store_bubbles_barycentres(interf_name);

  // Suppression des bulles dans le domaine a eliminer pres des bords periodiques
  // lors de la sauvegarde.
  if (ncells_deleted_ > 0)
    supprimer_certaines_bulles_reelles();

  dumplata_ft_mesh(lata_name, "INTERFACES", 0);
  dumplata_ft_field(lata_name, "INTERFACES", "COMPO_CONNEXE", "ELEM", mesh.compo_connexe_facettes(), 0);
}

void IJK_Interfaces::postraiter_colors(Sortie& os, const double current_time) const
{
  if (!Process::je_suis_maitre())
    return;

  os << "# Impression des couleurs de bulles" << finl;
  os << "# colonne 1 : temps" << finl;
  os << "# colonne K : Couleur de la bulle K-2" << finl;
  const int n = through_yminus_.size_array();
  os << current_time << " ";
  for (int j = 0; j < n; j++)
    os << through_yminus_[j] << " ";
  os << finl;
}

void IJK_Interfaces::calculer_color(ArrOfInt& color) const
{
  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const int n = mesh.nb_facettes();
  color.resize_array(n);
  const ArrOfInt& compo_facettes = mesh.compo_connexe_facettes();
  for (int i = 0; i < n; i++)
    {
      //    if (mesh.facette_virtuelle(i)) {
      //      color[i] = -1000; // Valeur invalide pour les facettes virtuelles.
      //      continue;
      //    }
      const int compo = compo_facettes[i];
      // ignorer les bulles dupliquees
      if (compo < 0)
        {
          color[i] = -1000; // Valeur invalide pour les duplicatatas.
          continue;
        }
      const int value = through_yminus_[compo];
      color[i] = value;
    }
  //  mp_max_for_each_item(color);
}

// Sorte de dico renversant le tableau ghost_compo_converter_
// Cherche dans le tableau ghost_compo_converter_ la case qui a la composante
// fournie en argument.
// Retourne le numero de bulle ghost construit comme :
//      -1 - idx_case
// L'entier retourner prend donc une valeur entre -1 et -nb_bulles_ghost
// (inclus)
int IJK_Interfaces::get_ghost_number_from_compo(const int compo) const
{
  const int n = ghost_compo_converter_.size_array();
  for (int i = 0; i < n; i++)
    {
      if (ghost_compo_converter_[i] == compo)
        {
          return -1 - i;
        }
    }

  Cerr << "Exception dans IJK_Interfaces::get_ghost_number_from_compo(). "
       << " Compo demandee " << compo << " introuvable...";
  Process::exit();
  return -1000000;
}

void IJK_Interfaces::calculer_surface_bulles(ArrOfDouble& surfaces) const
{
  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const int n = mesh.nb_facettes();
  const int nbulles_reelles = get_nb_bulles_reelles();
  const int nbulles_ghost = get_nb_bulles_ghost();
  surfaces.resize_array(nbulles_reelles + nbulles_ghost, RESIZE_OPTIONS::NOCOPY_NOINIT);
  surfaces = 0.;
  const ArrOfDouble& surfaces_facettes = mesh.get_update_surface_facettes();
  const ArrOfInt& compo_facettes = mesh.compo_connexe_facettes();
  for (int i = 0; i < n; i++)
    {
      if (mesh.facette_virtuelle(i))
        continue;
      int compo = compo_facettes[i];
      // les bulles dupliquees a la fin :
      if (compo < 0)
        {
          // L'index de la bulle ghost est (entre -1 et -nbulles_ghost):
          const int idx_ghost = get_ghost_number_from_compo(compo);
          // On la place en fin de tableau :
          compo = nbulles_reelles - 1 - idx_ghost;
        }

      const double s = surfaces_facettes[i];
      surfaces[compo] += s;
    }
  mp_sum_for_each_item(surfaces);
}

// in : should be an array of field [Unit] for each facette, with value
// field*area [Unit*m2]
//      For instance, fields can be :
//                    interfacial_temperature(fa7) = Ti*surf;
//             interfacial_phin_ai[fa7] += phin*surf;
void IJK_Interfaces::compute_surface_average_per_bubble(const ArrOfDouble& surfaces, const ArrOfDouble& in,
                                                        ArrOfDouble& out) const
{
  const int nbulles = get_nb_bulles_reelles();
  out.resize_array(nbulles);
  out = 0.;
  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const int n = mesh.nb_facettes();
  const ArrOfInt& compo_facettes = mesh.compo_connexe_facettes();
  for (int i = 0; i < n; i++)
    {
      if (mesh.facette_virtuelle(i))
        continue;
      int compo = compo_facettes[i];
      // On passe les bulles ghosts...
      if (compo < 0)
        continue;

      const double s = in[i];
      out[compo] += s;
    }
  mp_sum_for_each_item(out);

  for (int c = 0; c < nbulles; c++)
    out[c] /= surfaces[c];
}
void IJK_Interfaces::update_surface_normale() const
{
  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const ArrOfDouble& surfaces_facettes = mesh.get_update_surface_facettes();
  const DoubleTab& normales_facettes = mesh.get_update_normale_facettes();
  if (surfaces_facettes.size_array() * normales_facettes.size_array() > 0)
    {
      // juste pour bien qu'elles soient utilisees
      Cerr << "interfacial surface and normale are sure to be up-to-date" << finl;
      return;
    }
  return;
}

void IJK_Interfaces::reset_flags_and_counters()
{
  has_computed_bubble_barycentres_ = false;
}

void IJK_Interfaces::read_bubbles_barycentres_old_new(const Nom& interf_name)
{
  if (reprise_ && use_barycentres_velocity_)
    read_barycentres_velocity_ = 1;
  if (read_barycentres_velocity_)
    {
      has_readen_barycentres_prev_ = true;
      bool tmp_readen_flag = false;
      Nom suffix_tmp = ".old";
      FixedVector<ArrOfDouble,3> bubbles_bary_old;
      tmp_readen_flag = read_bubbles_barycentres(interf_name, suffix_tmp, bubbles_bary_old);
      if (tmp_readen_flag)
        FixedVector_to_DoubleTab(bubbles_bary_old, bubbles_bary_old_);
      has_readen_barycentres_prev_ = has_readen_barycentres_prev_ && tmp_readen_flag;

      suffix_tmp = ".new";
      FixedVector<ArrOfDouble,3> bubbles_bary_new;
      tmp_readen_flag = read_bubbles_barycentres(interf_name, suffix_tmp, bubbles_bary_new);
      if (tmp_readen_flag)
        FixedVector_to_DoubleTab(bubbles_bary_new, bubbles_bary_new_);
      has_readen_barycentres_prev_ = has_readen_barycentres_prev_ && tmp_readen_flag;

      FixedVector<ArrOfDouble,3> bubbles_rising_dir;
      FixedVector<ArrOfDouble,3> bubbles_rising_vel;
      tmp_readen_flag = read_bubbles_barycentres_vel(interf_name, bubbles_rising_dir, bubbles_rising_vel, bubbles_velocities_bary_magnitude_);
      if (tmp_readen_flag)
        {
          FixedVector_to_DoubleTab(bubbles_rising_dir, bubbles_rising_vectors_bary_);
          FixedVector_to_DoubleTab(bubbles_rising_vel, bubbles_velocities_bary_);
        }
      has_readen_barycentres_prev_ = has_readen_barycentres_prev_ && tmp_readen_flag;
    }
}

bool IJK_Interfaces::read_bubbles_barycentres_vel(const Nom& interf_name,
                                                  FixedVector<ArrOfDouble,3>& bubbles_rising_dir,
                                                  FixedVector<ArrOfDouble,3>& bubbles_rising_vel,
                                                  ArrOfDouble& bubbles_rising_vel_mag)
{
  bool is_readen = false;
  int bubbles_bary_computed = ijk_compo_connex_.get_compute_compo_fields();
  bubbles_bary_computed = bubbles_bary_computed || use_barycentres_velocity_;
  bubbles_rising_vel_mag.reset();
  for (int c=0; c<3; c++)
    {
      bubbles_rising_dir[c].reset();
      bubbles_rising_vel[c].reset();
    }
  if (bubbles_bary_computed)
    {
      int line_counter = 0;
      int var_index = 0;
      int dir = 0;
      std::string line;
      const Nom interf_dir = dirname(interf_name);
      const Nom case_name = (Objet_U::nom_du_cas());
      const Nom suffix_bary = ".sauv.barycentres";
      const Nom file_folder = interf_dir + case_name + suffix_bary + ".vel";
      // creating an ifstream object named file
      ifstream read_tmp;
      read_tmp.open(file_folder);
      const bool read_file = read_tmp ? true : false;
      read_tmp.close();
      if (read_file)
        {
          EFichier fic_bary_vel(file_folder);
          ifstream& ifstream_bary_vel = fic_bary_vel.get_ifstream();
          const char delimiter = ' ';
          Cerr << "Read coordinates of bubbles barycentres from: " << file_folder << finl;
          while (std::getline(ifstream_bary_vel, line))
            {
              std::stringstream ssline(line);
              Cerr << "Line number: " << line_counter << finl;
              while (std::getline(ssline, line, delimiter))
                {
                  if(line_counter)
                    {
                      Cerr << "Param: " << line << finl;
                      switch(var_index)
                        {
                        case 0:
                          break;
                        case 1:
                        case 2:
                        case 3:
                          dir = var_index - 1;
                          bubbles_rising_vel[dir].append_array(std::stod(line));
                          break;
                        case 4:
                        case 5:
                        case 6:
                          dir = var_index - 4;
                          bubbles_rising_dir[dir].append_array(std::stod(line));
                          break;
                        case 7:
                          bubbles_rising_vel_mag.append_array(std::stod(line));
                          break;
                        default:
                          break;
                        }
                      var_index++;
                    }
                }
              var_index = 0;
              line_counter++;
            }
          fic_bary_vel.close();
          is_readen = true;
        }
    }
  return is_readen;
}

bool IJK_Interfaces::read_bubbles_barycentres(const Nom& interf_name, const Nom& suffix, FixedVector<ArrOfDouble,3>& bubbles_bary)
{
  bool is_readen = false;
  int bubbles_bary_computed = ijk_compo_connex_.get_compute_compo_fields();
  bubbles_bary_computed = bubbles_bary_computed || use_barycentres_velocity_;
  for (int c=0; c<3; c++)
    {
      bubbles_bary[c].reset();
    }
  if (bubbles_bary_computed)
    {
      int line_counter = 0;
      int var_index = 0;
      int dir = 0;
      std::string line;
      const Nom interf_dir = dirname(interf_name);
      const Nom case_name = (Objet_U::nom_du_cas());
      const Nom suffix_bary = ".sauv.barycentres";
      const Nom file_folder = interf_dir + case_name + suffix_bary + suffix;
      ifstream read_tmp;
      read_tmp.open(file_folder);
      const bool read_file = read_tmp ? true : false;
      read_tmp.close();
      if (read_file)
        {
          EFichier fic_bary(file_folder);
          ifstream& ifstream_bary_old = fic_bary.get_ifstream();
          const char delimiter = ' ';
          Cerr << "Read coordinates of bubbles barycentres from: " << file_folder << finl;
          while (std::getline(ifstream_bary_old, line))
            {
              std::stringstream ssline(line);
              Cerr << "Line number: " << line_counter << finl;
              while (std::getline(ssline, line, delimiter))
                {
                  if(line_counter)
                    {
                      Cerr << "Param: " << line << finl;
                      switch(var_index)
                        {
                        case 0:
                          break;
                        case 1:
                        case 2:
                        case 3:
                          dir = var_index - 1;
                          bubbles_bary[dir].append_array(std::stod(line));
                          break;
                        default:
                          break;
                        }
                      var_index++;
                    }
                }
              var_index = 0;
              line_counter++;
            }
          fic_bary.close();
          is_readen = true;
        }
    }
  return is_readen;
}

void IJK_Interfaces::store_bubbles_barycentres(const Nom& interf_name)
{
  const Nom end_space = " ";
  const Nom escape = "\n";
  int bubbles_bary_computed = ijk_compo_connex_.get_compute_compo_fields();
  bubbles_bary_computed = bubbles_bary_computed || use_barycentres_velocity_;
  if (Process::je_suis_maitre() && bubbles_bary_computed)
    {
      const Nom interf_dir = dirname(interf_name);
      const Nom suffix = ".sauv.barycentres";
      const int reset = 1;
      const Nom bary_header_old = "ibubble bary_old_x bary_old_y bary_old_z";
      const Nom bary_header_new = "ibubble bary_new_x bary_new_y bary_new_z";
      const Nom bary_vel_header = "ibubble bary_vel_x bary_vel_y bary_vel_z bary_vect_x bary_vect_y bary_vect_z bary_vel_val";
      const int nbulles_reelles = get_nb_bulles_reelles();
      SFichier fic_bary_old = Open_file_folder(interf_dir, suffix + ".old", bary_header_old, reset, 0);
      for (int ibubble=0; ibubble<nbulles_reelles; ibubble++)
        {
          fic_bary_old << ibubble << end_space;
          fic_bary_old << bubbles_bary_old_(ibubble, 0) << end_space;
          fic_bary_old << bubbles_bary_old_(ibubble, 1) << end_space;
          fic_bary_old << bubbles_bary_old_(ibubble, 2) << escape;
        }
      fic_bary_old.close();
      SFichier fic_bary_new = Open_file_folder(interf_dir, suffix + ".new", bary_header_new, reset, 0);
      for (int ibubble=0; ibubble<nbulles_reelles; ibubble++)
        {
          fic_bary_new << ibubble << end_space;
          fic_bary_new << bubbles_bary_new_(ibubble, 0) << end_space;
          fic_bary_new << bubbles_bary_new_(ibubble, 1) << end_space;
          fic_bary_new << bubbles_bary_new_(ibubble, 2) << escape;
        }
      fic_bary_new.close();
      SFichier fic_bary_vel = Open_file_folder(interf_dir, suffix + ".vel", bary_vel_header, reset, 0);
      for (int ibubble=0; ibubble<nbulles_reelles; ibubble++)
        {
          fic_bary_vel << ibubble << end_space;
          fic_bary_vel << bubbles_velocities_bary_(ibubble, 0) << end_space;
          fic_bary_vel << bubbles_velocities_bary_(ibubble, 1) << end_space;
          fic_bary_vel << bubbles_velocities_bary_(ibubble, 2) << end_space;
          fic_bary_vel << bubbles_rising_vectors_bary_(ibubble, 0) << end_space;
          fic_bary_vel << bubbles_rising_vectors_bary_(ibubble, 1) << end_space;
          fic_bary_vel << bubbles_rising_vectors_bary_(ibubble, 2) << end_space;
          fic_bary_vel << bubbles_velocities_bary_magnitude_(ibubble) << escape;
        }
      fic_bary_vel.close();
    }
}

void IJK_Interfaces::compute_bubbles_volume_and_barycentres(ArrOfDouble& volumes,
                                                            DoubleTab& barycentres,
                                                            const int& store_values)
{
  calculer_volume_bulles(volumes, barycentres);
  const Domaine_IJK& geom = I().get_domaine();
  if (store_values && !has_computed_bubble_barycentres_)
    {
      if (ref_ijk_ft_->schema_temps_ijk().get_tstep() == 0)
        {
          if (!has_readen_barycentres_prev_)
            bubbles_bary_old_ = barycentres;
          else
            bubbles_bary_old_ = bubbles_bary_new_;
        }
      else
        bubbles_bary_old_ = bubbles_bary_new_;
      bubbles_bary_new_ = barycentres;
      const int nbulles_reelles = get_nb_bulles_reelles();
      bubbles_velocities_bary_magnitude_.resize(nbulles_reelles);
      bubbles_velocities_bary_ = bubbles_bary_old_;
      for (int i = 0; i < nbulles_reelles; i++)
        {
          bubbles_velocities_bary_magnitude_(i) = 0.;
          for (int dir=0; dir<3; dir++)
            {
              const double ldir = geom.get_domain_length(dir);
              const double pos_old = bubbles_bary_old_(i, dir);
              const double pos_new = bubbles_bary_new_(i, dir);
              const int old_new_up_down = (pos_new - pos_old) < (-ldir / 2.);
              const int old_new_down_up = (pos_new - pos_old) > (ldir / 2.);
              if (old_new_up_down || old_new_down_up)
                {
                  if (old_new_up_down)
                    bubbles_bary_old_(i, dir) += (-ldir);
                  else
                    bubbles_bary_old_(i, dir) += (ldir);
                }
              bubbles_velocities_bary_(i, dir) = bubbles_bary_old_(i, dir);
              const double vel_old = bubbles_velocities_bary_(i, dir);
              bubbles_velocities_bary_(i, dir) = bubbles_bary_new_(i, dir) - vel_old;
              bubbles_velocities_bary_(i, dir) *= (1 / ref_ijk_ft_->schema_temps_ijk().get_timestep());
              bubbles_velocities_bary_magnitude_(i) += pow(bubbles_velocities_bary_(i, dir), 2);
            }
          bubbles_velocities_bary_magnitude_(i) = sqrt(bubbles_velocities_bary_magnitude_(i));
          bubbles_rising_vectors_bary_ = bubbles_velocities_bary_;
          for (int dir=0; dir<3; dir++)
            {
              if (abs(bubbles_velocities_bary_magnitude_(i)) > DMINFLOAT)
                bubbles_rising_vectors_bary_(i, dir) /= bubbles_velocities_bary_magnitude_(i);
              else
                bubbles_rising_vectors_bary_(i, dir) = 0.;
            }
        }
      has_computed_bubble_barycentres_ = true;
    }
}

void IJK_Interfaces::calculer_volume_bulles(ArrOfDouble& volumes,
                                            DoubleTab& centre_gravite) const
{
  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const int n = mesh.nb_facettes();
  const int nbulles_reelles = get_nb_bulles_reelles();
  const int nbulles_ghost = get_nb_bulles_ghost();
  const int nbulles_tot = nbulles_reelles + nbulles_ghost;
  volumes.resize_array(nbulles_tot, RESIZE_OPTIONS::NOCOPY_NOINIT);
  volumes = 0.;
  centre_gravite.resize(nbulles_tot, 3, RESIZE_OPTIONS::NOCOPY_NOINIT);
  centre_gravite = 0.;
  const ArrOfDouble& surfaces_facettes = mesh.get_update_surface_facettes();
  const DoubleTab& normales_facettes = mesh.get_update_normale_facettes();
  const IntTab& facettes = mesh.facettes();
  const DoubleTab& sommets = mesh.sommets();
  const ArrOfInt& compo_facettes = mesh.compo_connexe_facettes();
  for (int i = 0; i < n; i++)
    {
      if (mesh.facette_virtuelle(i))
        continue;
      int compo = compo_facettes[i];
      // les bulles dupliquees a la fin :
      if (compo < 0)
        {
          // L'index de la bulle ghost est (entre -1 et -nbulles_ghost):
          const int idx_ghost = get_ghost_number_from_compo(compo);
          // On la place en fin de tableau :
          compo = nbulles_reelles - 1 - idx_ghost;
        }
      const double s = surfaces_facettes[i];
      const double normale_scalaire_direction = normales_facettes(i, 0); // On projette sur x
      // Coordonnee du centre de gravite de la facette
      const int i0 = facettes(i, 0);
      const int i1 = facettes(i, 1);
      const int i2 = facettes(i, 2);
      const double coord_centre_gravite_i = (sommets(i0, 0) + sommets(i1, 0) + sommets(i2, 0)) / 3.;
      const double coord_centre_gravite_j = (sommets(i0, 1) + sommets(i1, 1) + sommets(i2, 1)) / 3.;
      const double coord_centre_gravite_k = (sommets(i0, 2) + sommets(i1, 2) + sommets(i2, 2)) / 3.;
      const double volume_prisme = coord_centre_gravite_i * s * normale_scalaire_direction;
      // centre de gravite du prisme pondere par son volume avec signe
      centre_gravite(compo, 0) += volume_prisme * (coord_centre_gravite_i * 0.5);
      centre_gravite(compo, 1) += volume_prisme * coord_centre_gravite_j;
      centre_gravite(compo, 2) += volume_prisme * coord_centre_gravite_k;
      volumes[compo] += volume_prisme;
    }
  mp_sum_for_each_item(volumes);
  mp_sum_for_each_item(centre_gravite);
  //Cerr << "volumes : " << volumes << finl;
  for (int i = 0; i < nbulles_tot; i++)
    {
      // const double x = 1./volumes[i];
      const double x = (volumes[i] == 0.) ? 0. : 1. / volumes[i];
      centre_gravite(i, 0) *= x;
      centre_gravite(i, 1) *= x;
      centre_gravite(i, 2) *= x;
    }
}


void IJK_Interfaces::calculer_aspect_ratio(ArrOfDouble& aspect_ratio) const
{
  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const int n = mesh.nb_facettes();
  const int nbulles_reelles = get_nb_bulles_reelles();
  const int nbulles_ghost = get_nb_bulles_ghost();
  const int nbulles_tot = nbulles_reelles + nbulles_ghost;
  const IntTab& facettes = mesh.facettes();
  const DoubleTab& sommets = mesh.sommets();
  const ArrOfInt& compo_facettes = mesh.compo_connexe_facettes();
  aspect_ratio.resize_array(nbulles_tot);

  ArrOfDouble volumes;
  DoubleTab centre_gravite;
  this->calculer_volume_bulles(volumes,centre_gravite);

  DoubleTab d_max(nbulles_tot);
  d_max = -1;
  DoubleTab d_min(nbulles_tot);
  d_min = 300;

  double d_imax;
  double d_imin;

  for (int i = 0; i < n; i++)
    {
      if (mesh.facette_virtuelle(i))
        continue;
      int compo = compo_facettes[i];
      // les bulles dupliquees a la fin :
      if (compo < 0)
        {
          // L'index de la bulle ghost est (entre -1 et -nbulles_ghost):
          const int idx_ghost = get_ghost_number_from_compo(compo);
          // On la place en fin de tableau :
          compo = nbulles_reelles - 1 - idx_ghost;
        }
      // Calcul des distances entre sommet_i et le centre de gravite
      const int i0 = facettes(i, 0);
      const int i1 = facettes(i, 1);
      const int i2 = facettes(i, 2);
      const double d_i_0 = sqrt( (sommets(i0, 0)-centre_gravite(compo,0))*(sommets(i0, 0)-centre_gravite(compo,0)) + (sommets(i0, 1)-centre_gravite(compo,1))*(sommets(i0, 1)-centre_gravite(compo,1)) + (sommets(i0, 2)-centre_gravite(compo,2))*(sommets(i0, 2)-centre_gravite(compo,2)) );
      const double d_i_1 = sqrt( (sommets(i1, 0)-centre_gravite(compo,0))*(sommets(i1, 0)-centre_gravite(compo,0)) + (sommets(i1, 1)-centre_gravite(compo,1))*(sommets(i1, 1)-centre_gravite(compo,1)) + (sommets(i1, 2)-centre_gravite(compo,2))*(sommets(i1, 2)-centre_gravite(compo,2)) );
      const double d_i_2 = sqrt( (sommets(i2, 0)-centre_gravite(compo,0))*(sommets(i2, 0)-centre_gravite(compo,0)) + (sommets(i2, 1)-centre_gravite(compo,1))*(sommets(i2, 1)-centre_gravite(compo,1)) + (sommets(i2, 2)-centre_gravite(compo,2))*(sommets(i2, 2)-centre_gravite(compo,2)) );

      // On recupere la plus grande distance et la plus petite distance parmi les 3 calculees
      d_imax = std::max(d_i_0,std::max(d_i_1,d_i_2));
      d_imin = std::min(d_i_0,std::min(d_i_1,d_i_2));

      // On met a jour le grand axe et le petit axe
      if (d_imax > d_max[compo])
        d_max[compo] = d_imax;

      if (d_imin < d_min[compo])
        d_min[compo] = d_imin;
    }

  mp_min_for_each_item(d_min);
  mp_max_for_each_item(d_max);
  for (int i = 0; i < nbulles_tot; i++)
    {
      aspect_ratio[i] = d_max[i]/d_min[i];
    }
}

void IJK_Interfaces::calculer_surfactant(ArrOfDouble& surfactant,ArrOfDouble& surfactant_min,ArrOfDouble& surfactant_max) const
{
  if (maillage_ft_ijk_.Surfactant_facettes().get_disable_surfactant())
    {
      return ;
    }
  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const int n = mesh.nb_facettes();
  const int nbulles_reelles = get_nb_bulles_reelles();
  const ArrOfInt& compo_facettes = mesh.compo_connexe_facettes();
  const ArrOfDouble& Surf = mesh.Surfactant_facettes().get_FT_field_Array();
  const ArrOfDouble& Sfa7 = mesh.get_update_surface_facettes();
  surfactant.resize_array(nbulles_reelles);
  surfactant_max.resize_array(nbulles_reelles);
  surfactant_min.resize_array(nbulles_reelles);
  for (int bulle = 0; bulle < nbulles_reelles; bulle++)
    {
      surfactant_min(bulle) = 1.e10;
      surfactant_max(bulle) = -1.e10;
    }
  for (int i = 0; i < n; i++)
    {
      if (mesh.facette_virtuelle(i) or compo_facettes(i)<0)
        continue;
      surfactant(compo_facettes(i))+=Surf(i)*Sfa7(i);
      if(Surf(i)>surfactant_max(compo_facettes(i)))
        {
          surfactant_max(compo_facettes(i)) = Surf(i);
        }
      if(Surf(i)<surfactant_min(compo_facettes(i)))
        {
          surfactant_min(compo_facettes(i)) = Surf(i);
        }
    }

  for (int bulle = 0; bulle < nbulles_reelles; bulle++)
    {
      surfactant(bulle)=Process::mp_sum(surfactant(bulle));
      surfactant_max(bulle)=Process::mp_max(surfactant_max(bulle));
      surfactant_min(bulle)=Process::mp_min(surfactant_min(bulle));
    }
}


void IJK_Interfaces::calculer_poussee_bulles(const DoubleTab& grav,
                                             DoubleTab& poussee) const
{
  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const int n = mesh.nb_facettes();
  const int nbulles_reelles = get_nb_bulles_reelles();
  const int nbulles_ghost = get_nb_bulles_ghost();
  const int nbulles_tot = nbulles_reelles + nbulles_ghost;
  poussee.resize(nbulles_tot, 3, RESIZE_OPTIONS::NOCOPY_NOINIT);
  poussee = 0.;
  const ArrOfDouble& surfaces_facettes = mesh.get_update_surface_facettes();
  const DoubleTab& normales_facettes = mesh.get_update_normale_facettes();
  const IntTab& facettes = mesh.facettes();
  const DoubleTab& sommets = mesh.sommets();
  const ArrOfInt& compo_facettes = mesh.compo_connexe_facettes();
  for (int i = 0; i < n; i++)
    {
      if (mesh.facette_virtuelle(i))
        continue;
      int compo = compo_facettes[i];
      // les bulles dupliquees a la fin :
      if (compo < 0)
        {
          // L'index de la bulle ghost est (entre -1 et -nbulles_ghost):
          const int idx_ghost = get_ghost_number_from_compo(compo);
          // On la place en fin de tableau :
          compo = nbulles_reelles - 1 - idx_ghost;
        }
      const double s = surfaces_facettes[i];
      // Coordonnee du centre de gravite de la facette
      const int i0 = facettes(i, 0);
      const int i1 = facettes(i, 1);
      const int i2 = facettes(i, 2);
      const double coord_centre_gravite_i = (sommets(i0, 0) + sommets(i1, 0) + sommets(i2, 0)) / 3.;
      const double coord_centre_gravite_j = (sommets(i0, 1) + sommets(i1, 1) + sommets(i2, 1)) / 3.;
      const double coord_centre_gravite_k = (sommets(i0, 2) + sommets(i1, 2) + sommets(i2, 2)) / 3.;
      const double grav_scalaire_position_fois_s =
        (grav(0,0) * coord_centre_gravite_i + grav(0,1) * coord_centre_gravite_j + grav(0,2) * coord_centre_gravite_k) * s;
      for (int dir = 0; dir < 3; dir++)
        poussee(compo, dir) += grav_scalaire_position_fois_s * normales_facettes(i, dir);
    }
  mp_sum_for_each_item(poussee);
}

void IJK_Interfaces::calculer_aire_interfaciale(IJK_Field_double& ai) const
{

  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const Intersections_Elem_Facettes& intersections = mesh.intersections_elem_facettes();
  const ArrOfDouble& surface_facettes = mesh.get_update_surface_facettes();

  const int n = mesh.nb_facettes();
  const Domaine_IJK& geom = ai.get_domaine();
  const double dxi = geom.get_constant_delta(DIRECTION_I);
  const double dxj = geom.get_constant_delta(DIRECTION_J);
  const double dxk = geom.get_constant_delta(DIRECTION_K);
  const double vol = dxi * dxj * dxk;

  // Remise a zero du champ.
  ai.data() = 0;

  for (int fa7 = 0; fa7 < n; fa7++)
    {

      // On compte aussi les facettes_virtuelles
      // if (maillage.facette_virtuelle(fa7))
      //  continue;

      const double sf = surface_facettes[fa7];
      int index = intersections.index_facette()[fa7];
      while (index >= 0)
        {
          const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
          const int num_elem = data.numero_element_;
          // Anciennement la methode etait portee par le mesh :
          //      const Int3 ijk = mesh.convert_packed_to_ijk_cell(num_elem);
          // A present, elle est dans le splitting :
          const Int3 ijk = geom.convert_packed_to_ijk_cell(num_elem);
          const double surf = data.fraction_surface_intersection_ * sf;
          ai(ijk[0], ijk[1], ijk[2]) += surf / vol;
          index = data.index_element_suivant_;
        }
    }
}

void IJK_Interfaces::calculer_aire_interfaciale_for_compo(IJK_Field_double& ai, const int compo) const
{

  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const Intersections_Elem_Facettes& intersections = mesh.intersections_elem_facettes();
  const ArrOfDouble& surface_facettes = mesh.get_update_surface_facettes();
  const ArrOfInt& compo_connex = mesh.compo_connexe_facettes();

  const int n = mesh.nb_facettes();
  const Domaine_IJK& geom = ai.get_domaine();
  const double dxi = geom.get_constant_delta(DIRECTION_I);
  const double dxj = geom.get_constant_delta(DIRECTION_J);
  const double dxk = geom.get_constant_delta(DIRECTION_K);
  const double vol = dxi * dxj * dxk;

  // Remise a zero du champ.
  ai.data() = 0;

  for (int fa7 = 0; fa7 < n; fa7++)
    {

      // On compte aussi les facettes_virtuelles
      // if (maillage.facette_virtuelle(fa7))
      //  continue;
      int compo_fa7 = compo_connex[fa7];

      // if compo reelle, on fait rien
      // si compo virtuelle, on stocke les ai des bulles virtuelles a la fin du tableau
      //if (compo_fa7 < 0)
      //  compo_fa7 = - compo_fa7 + nb_bulles_reelles_ - 1 ;

      if (compo_fa7 != compo)
        continue;


      const double sf = surface_facettes[fa7];
      int index = intersections.index_facette()[fa7];
      while (index >= 0)
        {
          const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
          const int num_elem = data.numero_element_;
          // Anciennement la methode etait portee par le mesh :
          //      const Int3 ijk = mesh.convert_packed_to_ijk_cell(num_elem);
          // A present, elle est dans le splitting :
          const Int3 ijk = geom.convert_packed_to_ijk_cell(num_elem);
          const double surf = data.fraction_surface_intersection_ * sf;
          ai(ijk[0], ijk[1], ijk[2]) += surf / vol;
          index = data.index_element_suivant_;
        }
    }
}

double IJK_Interfaces::calculer_aire_interfaciale_for_compo(const int compo, const int i_ref, const int j_ref, const int k_ref) const
{

  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const Intersections_Elem_Facettes& intersections = mesh.intersections_elem_facettes();
  const ArrOfDouble& surface_facettes = mesh.get_update_surface_facettes();
  const ArrOfInt& compo_connex = mesh.compo_connexe_facettes();

  const int n = mesh.nb_facettes();
  const Domaine_IJK& geom = mesh.get_domaine();
  const double dxi = geom.get_constant_delta(DIRECTION_I);
  const double dxj = geom.get_constant_delta(DIRECTION_J);
  const double dxk = geom.get_constant_delta(DIRECTION_K);
  const double vol = dxi * dxj * dxk;
  double ai = 0 ;
  for (int fa7 = 0; fa7 < n; fa7++)
    {

      // On compte aussi les facettes_virtuelles
      // if (maillage.facette_virtuelle(fa7))
      //  continue;

      const double sf = surface_facettes[fa7];
      const int compo_fa7 = compo_connex[fa7];

      if (compo_fa7 != compo)
        continue;

      int index = intersections.index_facette()[fa7];
      while (index >= 0)
        {
          const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
          const int num_elem = data.numero_element_;
          // Anciennement la methode etait portee par le mesh :
          //      const Int3 ijk = mesh.convert_packed_to_ijk_cell(num_elem);
          // A present, elle est dans le splitting :
          const Int3 ijk = geom.convert_packed_to_ijk_cell(num_elem);
          const double surf = data.fraction_surface_intersection_ * sf;
          if (ijk[0]==i_ref and ijk[1]==j_ref and ijk[2]==k_ref)
            {
              ai += surf / vol;
            }
          index = data.index_element_suivant_;
        }
    }
  return ai;
}

// Retourne le code contenant a la fois le numero de bulle et le code de
// deplacement.
static int encoder_compo(int num_bulle, int code_deplacement)
{
  // zzzzzzzzzzzzzzxxxyyy :
  // z = Numero de compo connexe.
  // x = bit de signe par direction (0: negatif ou 1:positif).
  // y = bit de mouvement par direction (0:pas mouvement ou 1: mouvement)
  return (num_bulle << 6) | (code_deplacement);
}

static int decoder_numero_bulle(const int code)
{
  const int num_bulle = code >> 6;
  return num_bulle;
}

// Calcule le champ de courbure moyenne sur le domaine etendu IJK_ft
// Utile pour gerer les conditions de shear-periodicite :: interpolation de la pression monofluide
void IJK_Interfaces::calculer_kappa_ft(IJK_Field_double& kappa_ft)
{
  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const Intersections_Elem_Facettes& intersections = mesh.intersections_elem_facettes();
  const ArrOfDouble& surface_facettes = mesh.get_update_surface_facettes();
  const IntTab& facettes = mesh.facettes();
  const ArrOfDouble& courbure = maillage_ft_ijk_.get_update_courbure_sommets();

  const int n = mesh.nb_facettes();
  const Domaine_IJK& s = kappa_ft.get_domaine();


  IJK_Field_double SI_ft;
  SI_ft.allocate(s, Domaine_IJK::ELEM, 0);
  SI_ft.data() = 0.;

  kappa_ft.data() = 0.;

  for (int fa7 = 0; fa7 < n; fa7++)
    {

      // On compte aussi les facettes_virtuelles
      //if (mesh.facette_virtuelle(fa7))
      //  continue;
      const double sf=surface_facettes[fa7];
      int index=intersections.index_facette()[fa7];
      while (index >= 0)
        {
          const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
          const int num_elem = data.numero_element_;
          const Int3 ijk = s.convert_packed_to_ijk_cell(num_elem);
          const double surf = data.fraction_surface_intersection_ * sf;
          for (int isom = 0; isom< 3; isom++)
            {
              const int num_som = facettes(fa7, isom);
              const double kappa = courbure[num_som];
              kappa_ft(ijk[0],ijk[1],ijk[2]) += kappa*surf/3.;
            }
          SI_ft(ijk[0],ijk[1],ijk[2]) += surf;
          index = data.index_element_suivant_;
        }
    }


  const int nx = kappa_ft.ni();
  const int ny = kappa_ft.nj();
  const int nz = kappa_ft.nk();

  for (int k = 0; k < nz; k++)
    for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++)
        {
          if (kappa_ft(i,j,k)!=0.)
            kappa_ft(i,j,k)/=SI_ft(i,j,k);
        }
}


/** Le champ de normale n'est pas sur une grille decallee.
 * Il doit etre a la meme localisation que "ai" : Domaine_IJK::ELEM
 * Le champ kappa_ai contient le produit de la courbure moyenne sur la cellule
 * eulerienne par l'aire interfaciale dans cette cellule,
 * divisee par le volume de la cellule.

 * Le calcul repose sur la conversion vdf -> ijk du numero : num_elem
 * Pour que cette conversion soit valide, il faut que le champ soit sur le
 * splitting_ft_ car le dom_vdf n'est pas construit pour le splitting ns. Par
 * definition, mettre igroup a -1 pour inclure toutes les bulles
 */
void IJK_Interfaces::calculer_normales_et_aires_interfaciales(IJK_Field_double& ai, IJK_Field_double& kappa_ai,
                                                              IJK_Field_vector3_double& normale_cell,
                                                              const int igroup) const
{
  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const Intersections_Elem_Facettes& intersections = mesh.intersections_elem_facettes();
  const ArrOfDouble& surface_facettes = mesh.get_update_surface_facettes();
  const DoubleTab& normale_facettes = mesh.get_update_normale_facettes();
  const IntTab& facettes = mesh.facettes();
  const ArrOfDouble& courbure = maillage_ft_ijk_.get_update_courbure_sommets();

  const int n = mesh.nb_facettes();
  const Domaine_IJK& geom = ai.get_domaine();
  const double dxi = geom.get_constant_delta(DIRECTION_I);
  const double dxj = geom.get_constant_delta(DIRECTION_J);
  const double dxk = geom.get_constant_delta(DIRECTION_K);
  const double vol = dxi * dxj * dxk;

  // Remise a zero du champ.
  ai.data() = 0.;
  kappa_ai.data() = 0.;
  for (int dir = 0; dir < 3; dir++)
    normale_cell[dir].data() = 0.;

  const ArrOfInt& compo_facette = mesh.compo_connexe_facettes();
  for (int fa7 = 0; fa7 < n; fa7++)
    {

      // On compte aussi les facettes_virtuelles
      //if (mesh.facette_virtuelle(fa7))
      //  continue;

      // On ne veut comptabiliser que les facettes des compo appartenant a igroup
      // (ou toutes les compo si igroup == -1)
      int icompo = compo_facette[fa7];
      if (icompo<0)
        {
          // Portion d'interface ghost. On recherche le vrai numero
          icompo = decoder_numero_bulle(-icompo);
        }
      if ((compo_to_group_[icompo] != igroup) && (igroup != -1))
        continue;

      const double sf=surface_facettes[fa7];
      //    double surface_tot = 0.;
      int index=intersections.index_facette()[fa7];
      while (index >= 0)
        {
          const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
          const int num_elem = data.numero_element_;
          const Int3 ijk = geom.convert_packed_to_ijk_cell(num_elem);
          const double surf = data.fraction_surface_intersection_ * sf;
          for (int dir = 0; dir< 3; dir++)
            {
              const double nx = normale_facettes(fa7,dir);
              normale_cell[dir](ijk[0],ijk[1],ijk[2]) += nx * surf;
            }
          const double fac = surf/(3.*vol);
          for (int isom = 0; isom< 3; isom++)
            {
              const int num_som = facettes(fa7, isom);
              const double kappa = courbure[num_som];
              kappa_ai(ijk[0],ijk[1],ijk[2]) += kappa*fac;
            }
          //      surface_tot +=surf;
          ai(ijk[0],ijk[1],ijk[2]) += surf/vol;
          index = data.index_element_suivant_;
        }
    }

  // Nombre de mailles du domaine NS :
  // Quelle que soit la compo, le champ normale_cell a le meme nombre de mailles
  // qu'ai puisqu'il est localise entierement a ::ELEM
  const int nx = ai.ni();
  const int ny = ai.nj();
  const int nz = ai.nk();

  for (int k = 0; k < nz; k++)
    for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++)
        {
          double norme_carre = 0.;
          for (int dir = 0; dir < 3; dir++)
            {
              const double dnx = normale_cell[dir](i, j, k);
              norme_carre += dnx * dnx;
            }

          if (norme_carre > 0.)
            for (int dir = 0; dir < 3; dir++)
              normale_cell[dir](i, j, k) *= 1. / sqrt(norme_carre);
        }
}

// Je ne peux plus conserver cette methode statique a cause de l'utilisation de
// get_ghost_number_from_compo
void IJK_Interfaces::calculer_vmoy_translation_composantes_connexes(const Maillage_FT_IJK& maillage,
                                                                    const ArrOfDouble& surface_facette,
                                                                    const ArrOfDouble& surface_par_bulle,
                                                                    const ArrOfInt& compo_connexes_facettes,
                                                                    const int nbulles_reelles,
                                                                    const int nbulles_ghost,
                                                                    const DoubleTab& vitesse_sommets,
                                                                    DoubleTab& vitesses_translation_bulles) const
{
  const int nbulles_tot = nbulles_reelles + nbulles_ghost;
  assert(vitesses_translation_bulles.dimension(0) == nbulles_tot);
  assert(surface_par_bulle.size_array() == nbulles_tot);
  assert(vitesses_translation_bulles.dimension(1) == 3);

  const IntTab& facettes = maillage.facettes();
  // Calcul de la vitesse de deplacement moyenne
  vitesses_translation_bulles = 0.;

  // Cette option supprime la partie translation du mouvement rigide
  if (disable_rigid_translation_)
    return;


  // calcul de la vitesse moyenne de deplacement de l'interface
  for (int fa7 = 0; fa7 < maillage.nb_facettes(); fa7++)
    {

      // Ne prendre que les facettes reelles (sinon on les compte plusieurs fois)
      if (maillage.facette_virtuelle(fa7))
        continue;

      int compo = compo_connexes_facettes[fa7];
      // On met les bulles dupliquees a la fin
      // La premiere (numero -1) est a la position nbulles_reelles :
      if (compo < 0)
        {
          // L'index de la bulle ghost est (entre -1 et -nbulles_ghost):
          const int idx_ghost = get_ghost_number_from_compo(compo);
          // On la place en fin de tableau :
          compo = nbulles_reelles - 1 - idx_ghost;
        }
      assert(compo >= 0);
      assert(compo < nbulles_tot);

      const double sf = surface_facette[fa7];

      // Boucle sur les directions :
      for (int dir = 0; dir < 3; dir++)
        {
          double v = 0.;
          // Boucle sur les sommets
          for (int j = 0; j < 3; j++)
            {
              const int isom = facettes(fa7, j);
              v += vitesse_sommets(isom, dir);
            }
          vitesses_translation_bulles(compo, dir) += v * sf / 3.; // il y a 3 sommets.
        }
    }

  mp_sum_for_each_item(vitesses_translation_bulles);

  for (int icompo = 0; icompo < nbulles_tot; icompo++)
    if (surface_par_bulle[icompo] > 0.)
      for (int dir = 0; dir < 3; dir++)
        vitesses_translation_bulles(icompo, dir) /= surface_par_bulle[icompo];
}

void IJK_Interfaces::calculer_vmoy_rotation_composantes_connexes(const Maillage_FT_IJK& maillage,
                                                                 const ArrOfDouble& surface_facette,
                                                                 const ArrOfDouble& surface_par_bulle,
                                                                 const ArrOfInt& compo_connexes_facettes,
                                                                 const int nbulles_reelles,
                                                                 const int nbulles_ghost,
                                                                 const DoubleTab& centre_gravite,
                                                                 const DoubleTab& vitesse_sommets,
                                                                 const DoubleTab& vitesse_translation_sommets,
                                                                 DoubleTab& mean_bubble_rotation_vector) const
{
  const int nbulles_tot = nbulles_reelles + nbulles_ghost;
  assert(mean_bubble_rotation_vector.dimension(0) == nbulles_tot);
  assert(surface_par_bulle.size_array() == nbulles_tot);
  assert(mean_bubble_rotation_vector.dimension(1) == 3);

  const IntTab& facettes = maillage.facettes();
  const DoubleTab& sommets = maillage.sommets();
  // Calcul de la vitesse de deplacement moyenne
  mean_bubble_rotation_vector = 0.;

  // Cette option supprime la partie rotation du mouvement rigide
  if (disable_rigid_rotation_)
    return;

  // calcul de la vitesse moyenne de deplacement de l'interface
  for (int fa7 = 0; fa7 < maillage.nb_facettes(); fa7++)
    {
      // Ne prendre que les facettes reelles (sinon on les compte plusieurs fois)
      if (maillage.facette_virtuelle(fa7))
        continue;

      int compo = compo_connexes_facettes[fa7];
      // On met les bulles dupliquees a la fin
      // La premiere (numero -1) est a la position nbulles_reelles :
      if (compo < 0)
        {
          // L'index de la bulle ghost est (entre -1 et -nbulles_ghost):
          const int idx_ghost = get_ghost_number_from_compo(compo);
          // On la place en fin de tableau :
          compo = nbulles_reelles - 1 - idx_ghost;
        }
      assert(compo >= 0);
      assert(compo < nbulles_tot);

      const double sf = surface_facette[fa7];

      // Boucle sur les sommets
      for (int j = 0; j < 3; j++)
        {
          const int isom = facettes(fa7, j);

          double vit_x = vitesse_sommets(isom, 0) - vitesse_translation_sommets(compo, 0);
          double vit_y = vitesse_sommets(isom, 1) - vitesse_translation_sommets(compo, 1);
          double vit_z = vitesse_sommets(isom, 2) - vitesse_translation_sommets(compo, 2);
          Vecteur3 vit = {vit_x, vit_y, vit_z};

          double x = sommets(isom, 0);
          double y = sommets(isom, 1);
          double z = sommets(isom, 2);
          Vecteur3 coord_sommet = {x, y, z};

          Vecteur3 coord_centre = {centre_gravite(compo,0), centre_gravite(compo,1), centre_gravite(compo,2)};

          Vecteur3 radial_position = coord_sommet - coord_centre;

          Vecteur3 rotation_vector;
          Vecteur3::produit_vectoriel(radial_position, vit, rotation_vector);

          // Boucle sur les directions :
          for (int dir = 0; dir < 3; dir++)
            {
              rotation_vector[dir] /= Vecteur3::produit_scalaire(radial_position, radial_position);

              mean_bubble_rotation_vector(compo, dir) += rotation_vector[dir] * sf / 3.; // il y a 3 sommets.
            }
        }
    }

  mp_sum_for_each_item(mean_bubble_rotation_vector);

  for (int icompo = 0; icompo < nbulles_tot; icompo++)
    if (surface_par_bulle[icompo] > 0.)
      for (int dir = 0; dir < 3; dir++)
        mean_bubble_rotation_vector(icompo, dir) /= surface_par_bulle[icompo];
}

/*! @brief calcul du vecteur normal a l'interface, aux sommets du maillage d'interface.
 *
 * Le tableau "normale" est efface et resize.
 *   La normale est la moyenne des normales des facettes voisines, ponderees par
 *   la surface de la facette.
 *   La norme du vecteur normal n'est pas unitaire !
 *   L'espace virtuel n'est pas a jour !
 *
 */
static void calculer_normale_sommets_interface(const Maillage_FT_IJK& maillage,
                                               DoubleTab& normale)
{
  const int nsom = maillage.nb_sommets();
  const int nfaces = maillage.nb_facettes();
  const int dim = maillage.sommets().dimension(1);

  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  const DoubleTab& normale_facettes = maillage.get_update_normale_facettes();
  const IntTab& facettes = maillage.facettes();
  const int nsommets_faces = facettes.dimension(1);

  normale.resize(nsom, dim);
  normale = 0.;
  double n[3] = {0., 0. , 0.};
  // La surface autour d'un sommet :
  ArrOfDouble surface_environnante(nsom);
  surface_environnante = 0.;

  for (int i = 0; i < nfaces; i++)
    {

      // Somme pour les faces reelles:
      if (maillage.facette_virtuelle(i))
        continue;

      const double surface = surface_facettes[i];

      for (int k = 0; k < dim; k++)
        n[k] = normale_facettes(i, k) * surface;

      for (int j = 0; j < nsommets_faces; j++)
        {
          const int sommet = facettes(i, j);
          surface_environnante[sommet] += surface;
          for (int k = 0; k < dim; k++)
            normale(sommet, k) += n[k];
        }
    }

  // Sommer les contributions pour les sommets partages par plusieurs
  // processeurs:
  maillage.desc_sommets().collecter_espace_virtuel(normale, MD_Vector_tools::EV_SOMME);
  maillage.desc_sommets().collecter_espace_virtuel(surface_environnante, MD_Vector_tools::EV_SOMME);

  double norme;
  int print = 0;
  int count[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  for (int isom = 0; isom < nsom; isom++)
    if (!maillage.sommet_virtuel(isom))
      {
        for (int dir = 0; dir < 3; dir++)
          {
            if (surface_environnante[isom])
              normale(isom, dir) /= surface_environnante[isom];
            else
              normale(isom, dir) = 0.;
          }

        norme = normale(isom, 0) * normale(isom, 0) + normale(isom, 1) * normale(isom, 1) +
                normale(isom, 2) * normale(isom, 2);
        if (norme < 0.9)
          {
            print = 1;
            if (norme < 0.8)
              {
                if (norme < 0.7)
                  {
                    if (norme < 0.6)
                      {
                        if (norme < 0.5)
                          {
                            if (norme < 0.4)
                              {
                                if (norme < 0.3)
                                  {
                                    if (norme < 0.2)
                                      {
                                        if (norme < 0.1)
                                          {
                                            count[0]++;
                                          }
                                        else
                                          {
                                            count[1]++;
                                          }
                                      }
                                    else
                                      {
                                        count[2]++;
                                      }
                                  }
                                else
                                  {
                                    count[3]++;
                                  }
                              }
                            else
                              {
                                count[4]++;
                              }
                          }
                        else
                          {
                            count[5]++;
                          }
                      }
                    else
                      {
                        count[6]++;
                      }
                  }
                else
                  {
                    count[7]++;
                  }
              }
            else
              {
                count[8]++;
              }
          }
        //      Cerr << "Som " << isom << " norme " << norme << finl;
      }
  if (print)
    {
      Cerr << "IJK_Interfaces.cpp:calculer_normale_sommets_interface : Calcul "
           "normale non-unitaire :  ";
      double sum = 0;
      for (int i = 0; i < 9; i++)
        {
          sum += count[i];
          Cerr << " " << count[i] / nsom;
        }
      Cerr << " " << 1 - sum / nsom << " (nsom = " << nsom << finl;
    }
}

void IJK_Interfaces::calculer_var_volume_remaillage(double timestep,
                                                    const DoubleTab& vitesses_translation_bulles,
                                                    const DoubleTab& mean_bubble_rotation_vector,
                                                    const DoubleTab& centre_gravite,
                                                    ArrOfDouble& var_volume)
{
  Cut_field_vector3_double& cut_field_deformation_velocity = static_cast<Cut_field_vector3_double&>(deformation_velocity_);

  Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const DoubleTab& sommets = mesh.sommets(); // Tableau des coordonnees des marqueurs.
  int nbsom = sommets.dimension(0);

  const int nbulles_reelles = get_nb_bulles_reelles();
  const int nbulles_ghost = get_nb_bulles_ghost();
  const int nbulles_tot = nbulles_reelles + nbulles_ghost;

  // Initialisations
  var_volume.resize(nbsom);
  var_volume = 0.;

  delta_volume_theorique_bilan_ns_.data() = 0.;

  Cut_cell_surface_efficace::calcul_surface_face_efficace_initiale(
    false,
    indicatrice_surfacique_face_ns_[next()], // Pour l'instant, next() est la fin du pas de temps precedent
    indicatrice_surfacique_avant_remaillage_face_ns_,
    indicatrice_surfacique_efficace_deformation_face_,
    indicatrice_surfacique_efficace_deformation_face_);

  for (int dir = 0 ; dir < 3 ; dir++)
    {
      cut_field_deformation_velocity[dir].set_to_uniform_value(6.3e32); // Valeur absurde pour assurer que la valeur par defaut n'est jamais utilisee
    }

  DoubleTab bounding_box;
  calculer_bounding_box_bulles(bounding_box);

  // Calcul de la variation de volume cible separement pour chaque bulle.
  // En effet, chaque bulle produit une vitesse de deformation differente sur le maillage eulerien.
  // Calculer un champ de vitesse de deformation unique devrait restreindre le rapprochement des bulles.
  for (int icompo = 0; icompo < nbulles_tot; icompo++)
    {
      const Probleme_FTD_IJK_cut_cell& ref_ijk_ft_cut_cell_ = ref_cast(Probleme_FTD_IJK_cut_cell, ref_ijk_ft_.valeur());
      calculer_vitesse_de_deformation(icompo, bounding_box, ref_ijk_ft_cut_cell_.get_cut_field_velocity(), vitesses_translation_bulles, mean_bubble_rotation_vector, centre_gravite);
      //cut_field_deformation_velocity.echange_espace_virtuel(ghost);

      Cut_cell_surface_efficace::calcul_delta_volume_theorique_bilan(icompo, bounding_box, timestep,
                                                                     indicatrice_ns_[next()],
                                                                     indicatrice_avant_remaillage_ns_,
                                                                     indicatrice_surfacique_efficace_deformation_face_,
                                                                     cut_field_deformation_velocity,
                                                                     delta_volume_theorique_bilan_ns_);
    }

  // Transfert des donnees sur la variation de volume cible, du maillage eulerien vers le maillage lagrangien.
  // On converti d'abord les donnees sur un DoubleVect du maillage VDF, pour pouvoir utiliser la fonction Transport_Interfaces_FT_Disc::transfert_conservatif_eulerien_vers_lagrangien_sommets.
  DoubleVect delta_volume_theorique_bilan_ns_vdf;

  {
    const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, refdomaine_dis_.valeur());
    const Domaine& domaine = domaine_vf.domaine();
    domaine.creer_tableau_elements(delta_volume_theorique_bilan_ns_vdf);

    const int ni = delta_volume_theorique_bilan_ns_.ni();
    const int nj = delta_volume_theorique_bilan_ns_.nj();
    const int nk = delta_volume_theorique_bilan_ns_.nk();
    assert(ni == I_ft().ni());
    assert(nj == I_ft().nj());
    assert(nk == I_ft().nk());

    for (int k = 0; k < nk; k++)
      {
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                const int num_elem = ref_domaine_->convert_ijk_cell_to_packed(i, j, k); // Note: ref_domaine_ is a ft_splitting
                delta_volume_theorique_bilan_ns_vdf[num_elem] = delta_volume_theorique_bilan_ns_(i,j,k);
              }
          }
      }
  }

  Transport_Interfaces_FT_Disc::transfert_conservatif_eulerien_vers_lagrangien_sommets(mesh, delta_volume_theorique_bilan_ns_vdf, var_volume);
}

void IJK_Interfaces::calculer_vecteurs_de_deplacement_rigide(DoubleTab& vitesses_translation_bulles,
                                                             DoubleTab& mean_bubble_rotation_vector,
                                                             DoubleTab& centre_gravite,
                                                             const int first_step_interface_smoothing)
{
  Maillage_FT_IJK& mesh = maillage_ft_ijk_;

  // compute_vinterp(); // to resize and fill vinterp_
  compute_vinterp();
  if (first_step_interface_smoothing)
    vinterp_ = 0.;

  // Les sommets virtuels sont peut-etre trop loin pour pouvoir interpoler leur
  // vitesse, il faut faire un echange espace virtuel pour avoir leur vitesse.
  mesh.desc_sommets().echange_espace_virtuel(vinterp_);

  // Calcul d'un deplacement preservant la distribution des noeuds sur les
  // bulles: Inspiree de :
  // Transport_Interfaces_FT_Disc::calculer_vitesse_repere_local

  const int nbulles_reelles = get_nb_bulles_reelles();
  const int nbulles_ghost = get_nb_bulles_ghost();
  const int nbulles_tot = nbulles_reelles + nbulles_ghost;
  //  dvol.resize_array(nbulles_tot);
  //  dvol = 0.;
  const ArrOfInt& compo_connex = mesh.compo_connexe_facettes();
  //  assert(compo_connex.size_array() == 0 || min_array(compo_connex) >=0); //
  //  Les duplicatas ne sont pas presents pendant le transport.
  // Nouveau depuis le 13/03/2014 : Les bulles ghost sont autorisees lors du
  // transport...
  ArrOfDouble surface_par_bulle;
  calculer_surface_bulles(surface_par_bulle);
  const ArrOfDouble& surface_facette = mesh.get_update_surface_facettes();
  //  const DoubleTab& normale_facettes = mesh.get_update_normale_facettes();
  ArrOfDouble volume_par_bulle(nbulles_tot);
  vitesses_translation_bulles.resize(nbulles_tot, 3);
  mean_bubble_rotation_vector.resize(nbulles_tot, 3);
  centre_gravite.resize(nbulles_tot, 3);
  if (use_barycentres_velocity_)
    compute_bubbles_volume_and_barycentres(volume_par_bulle, centre_gravite, 1);
  else
    calculer_volume_bulles(volume_par_bulle, centre_gravite);


  // Calcul de la vitesse moyenne de chaque composante connexe :
  calculer_vmoy_translation_composantes_connexes(mesh,
                                                 surface_facette,
                                                 surface_par_bulle,
                                                 compo_connex,
                                                 nbulles_reelles,
                                                 nbulles_ghost,
                                                 vinterp_,
                                                 vitesses_translation_bulles);
  bubbles_velocities_ = vitesses_translation_bulles;

  if (use_barycentres_velocity_ && ref_ijk_ft_->schema_temps_ijk().get_tstep())
    vitesses_translation_bulles = bubbles_velocities_bary_;

  // Calcul de la vitesse due a la rotation de chaque composante connexe :
  calculer_vmoy_rotation_composantes_connexes(mesh,
                                              surface_facette,
                                              surface_par_bulle,
                                              compo_connex,
                                              nbulles_reelles,
                                              nbulles_ghost,
                                              centre_gravite,
                                              vinterp_,
                                              vitesses_translation_bulles,
                                              mean_bubble_rotation_vector);

  if (first_step_interface_smoothing)
    {
      vitesses_translation_bulles = 0.;
      mean_bubble_rotation_vector = 0.;
    }

#ifdef GB_VERBOSE
  if (Process::je_suis_maitre())
    {
      // ofstream fout;
      // fout.open("composantes_connexes.txt",ios::app);
      // fout << "TEMPS: " << "xxxxtempsxxxxx" << endl;
      for (int bulle = 0; bulle < nbulles_tot; bulle++)
        {
          Cerr << "composante " << bulle << " vitesse_translation ";
          for (int i = 0; i < 3; i++)
            Cerr << " " << vitesses_translation_bulles(bulle, i);
          Cerr << " rotation_vector ";
          for (int i = 0; i < 3; i++)
            Cerr << " " << mean_bubble_rotation_vector(bulle, i);
          Cerr << endl;
        }
      //   fout.close();
    }
#endif
}

// Valeur par defaut : rk_step = -1 si schema temps different de rk3.
// dvol : la variation de volume par bulle au cours du pas de temps.
// Cette methode est a present capable de transporter aussi les bulles ghost.
// Pre-requis : il faut que le tableau dvol soit bien dimensionne a nbulles_tot
// La fonction transporter_maillage a ete separee en trois pour inserer un
// calcul de l'indicatrice entre les differentes etapes.
// Cette fonction realise la partie deformation uniquement, ou dans ce cadre
// la deformation est la partie residuelle du mouvement lorsque translation
// et rotation solide sont retranchees au mouvement total.
void IJK_Interfaces::transporter_maillage_deformation(const int correction_semi_locale_volume_bulle,
                                                      const DoubleTab& vitesses_translation_bulles,
                                                      const DoubleTab& mean_bubble_rotation_vector,
                                                      const DoubleTab& centre_gravite,
                                                      const double dt_tot,
                                                      ArrOfDouble& dvol,
                                                      const int rk_step = -1,
                                                      const int first_step_interface_smoothing)
{
  // nouvelle version:
  Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const DoubleTab& sommets = mesh.sommets(); // Tableau des coordonnees des marqueurs.
  int nbsom = sommets.dimension(0);
  DoubleTab deplacement(nbsom, 3);

  // On suppose que vinterp_ est deja rempli

  // Calcul d'un deplacement preservant la distribution des noeuds sur les
  // bulles: Inspiree de :
  // Transport_Interfaces_FT_Disc::calculer_vitesse_repere_local

  const int nbulles_reelles = get_nb_bulles_reelles();
  const int nbulles_ghost = get_nb_bulles_ghost();
  const int nbulles_tot = nbulles_reelles + nbulles_ghost;
  //  dvol.resize_array(nbulles_tot);
  //  dvol = 0.;
  //  assert(compo_connex.size_array() == 0 || min_array(compo_connex) >=0); //
  //  Les duplicatas ne sont pas presents pendant le transport.
  // Nouveau depuis le 13/03/2014 : Les bulles ghost sont autorisees lors du
  // transport...
  //  const DoubleTab& normale_facettes = mesh.get_update_normale_facettes();
  ArrOfIntFT compo_connex_som;
  mesh.calculer_compo_connexe_sommets(compo_connex_som);


  // Calcul des normales (non-unitaire, espace virtuel pas a jour) aux sommets :
  DoubleTabFT normale_sommet;
  calculer_normale_sommets_interface(mesh, normale_sommet);

  for (int som = 0; som < nbsom; som++)
    {
      if (mesh.sommet_virtuel(som))
        {
          // Valeur pipo pour dire qu'on n'initialise pas
          deplacement(som, 0) = 100.;
          deplacement(som, 1) = 100.;
          deplacement(som, 2) = 100.;
          vinterp_(som, 0) = 100. / dt_tot;
          vinterp_(som, 1) = 100. / dt_tot;
          vinterp_(som, 2) = 100. / dt_tot;
          continue;
        }

      int icompo = compo_connex_som[som];
      // Les bulles dupliquees sont a la fin
      // La premiere (numero -1) est a la position nbulles_reelles :
      if (icompo < 0)
        {
          // L'index de la bulle ghost est (entre -1 et -nbulles_ghost):
          const int idx_ghost = get_ghost_number_from_compo(icompo);
          // On la place en fin de tableau :
          icompo = nbulles_reelles - 1 - idx_ghost;
        }
      assert(icompo >= 0);
      assert(icompo < nbulles_tot);

      // (v-vmoy) doit etre normal a l'interface
      // Donc on fait v_corrige = v_initial -
      // composante_tangentielle_de(v_initial-vmoy) Demonstration que (v_corrige -
      // vmoy) est normal a l'interface : On note ct() =
      // composante_tangentielle_de() ct(v_corrige - vmoy)
      //  = ct(v_corrige) - ct(vmoy)    car ct est lineaire
      //  = ct(v_initial - ct(v_initial) + ct(vmoy)) - ct(vmoy)    car ct
      //  est lineaire = ct(v_initial) - ct(v_initial) + ct(vmoy) - ct(vmoy) cat
      //  ct(ct(x)) = ct(x) et linearite = 0

      // v_corrige = cn(v_initial - vmoy) + vmoy
      // v_corrige = cn(v_initial) - cn(vmoy) + vmoy
      // v_corrige = cn(v_initial) + ct(vmoy)
      // v_corrige = v_initial - ct(v_inital) + ct(vmoy)
      // v_corrige = v_initial - ct(v_initial - vmoy)

      double rot_x = mean_bubble_rotation_vector(icompo, 0);
      double rot_y = mean_bubble_rotation_vector(icompo, 1);
      double rot_z = mean_bubble_rotation_vector(icompo, 2);
      Vecteur3 rot = {rot_x, rot_y, rot_z};

      double x = sommets(som, 0);
      double y = sommets(som, 1);
      double z = sommets(som, 2);
      Vecteur3 coord_sommet = {x, y, z};

      Vecteur3 coord_centre = {centre_gravite(icompo,0), centre_gravite(icompo,1), centre_gravite(icompo,2)};

      Vecteur3 radial_position = coord_sommet - coord_centre;

      Vecteur3 tangential_velocity;
      Vecteur3::produit_vectoriel(rot, radial_position, tangential_velocity);

      double prodscal = 0.;
      double norme_carre = 0.;
      for (int direction = 0; direction < 3; direction++)
        {
          double vi = vinterp_(som, direction);
          double vcompo = vitesses_translation_bulles(icompo, direction) + tangential_velocity[direction];
          double n = normale_sommet(som, direction);
          prodscal += (vi - vcompo) * n;
          norme_carre += n * n;
        }
      if (norme_carre)
        prodscal /= norme_carre;
      else
        prodscal = 0.;

      for (int direction = 0; direction < 3; direction++)
        {
          double n = normale_sommet(som, direction);
          // On enregirste la vitesse avec laquelle on souhaite deplacer les
          // marqueurs :

          if (correction_semi_locale_volume_bulle)
            {
              // Correction semi-locale du volume : on effectue une deformation de la bulle uniquement (la partie rigide du deplacement est realisee apres le remaillage).
              vinterp_(som, direction) = n * prodscal;
            }
          else
            {
              // Comportement par defaut (sans la correction semi-locale du volume) : on effectue un deplacement complet des marqueurs lagrangiens.
              vinterp_(som, direction) = n * prodscal + vitesses_translation_bulles(icompo, direction);
            }
        }
    }

  // Transport conservant le volume des bulles:
  // On va transporter avec la vitesse interpolee,
  // calculer la variation de volume de chaque bulle engendree par le
  // deplacement, calculer une correction de volume a appliquer a chaque sommet
  // pour retrouver le volume initial appeler le barycentrage_lissage avec la
  // correction de volume a appliquer

  // On met zero dans le tableau RK3_G_store_vi_
  if (rk_step == 0)
    {
      RK3_G_store_vi_.resize(nbsom, 3);
      RK3_G_store_vi_ *= 0.;
      mesh.desc_sommets().echange_espace_virtuel(RK3_G_store_vi_);
    }

  // On remplit le tableau de deplacement
  if (rk_step >= 0)
    {
      // Attention vinterp_ n'est pas a jour sur les sommets virtuels car il y a
      // un continue dans la boucle precedente
      //   Ce n'est pas grave pour le tableau deplacement car il ne sert qu'aux
      //   sommets reels... Par contre, vinterp peut etre stockee dans
      //   RK3_G_store_vi_ qui ne sera donc pas a jour pour les sommets virtuels.
#ifdef GB_VERBOSE
      Journal() << "RKSTEP  " << rk_step << finl;
      for (int i = 0; i < nbsom; i++)
        {
          Journal() << " Som " << i << " virt  " << mesh.sommet_virtuel(i) << " depl ";
          for (int direction = 0; direction < 3; direction++)
            {
              Journal() << vinterp(i, direction) << " ";
            }
          Journal() << " RK3_G_store_vi_ ";
          for (int direction = 0; direction < 3; direction++)
            {
              Journal() << RK3_G_store_vi_(i, direction) << " ";
            }
          Journal() << finl;
        }
#endif
      // Pas necessaire de mettre a jour l'EV car le transport n'utilise pas les
      // sommets virt :
      // mesh.desc_sommets().echange_espace_virtuel(vinterp);
      runge_kutta3_update(vinterp_, RK3_G_store_vi_, deplacement, rk_step, dt_tot, mesh);

      // La mise a jour de l'espace virt me semble inutile car mesh.transporter ne
      // deplace que les sommets reels
      // mesh.desc_sommets().echange_espace_virtuel(RK3_G_store_vi_);

      // On vide l'espace virtuel avant le transport :
      mesh.preparer_tableau_avant_transport(RK3_G_store_vi_, mesh.desc_sommets());
    }
  else
    {
      // Schema Euler :
      for (int som = 0; som < nbsom; som++)
        for (int direction = 0; direction < 3; direction++)
          deplacement(som, direction) = vinterp_(som, direction) * dt_tot;
    }

  // MaJ de l'EV par precaution pour le post-pro (et peut-etre calcul distance?)
  mesh.desc_sommets().echange_espace_virtuel(vinterp_);

  // Copie des coordonnees des sommets avant deplacement
  DoubleTab coord_sommets_avant_deplacement(mesh.sommets());
  mesh.preparer_tableau_avant_transport(coord_sommets_avant_deplacement, mesh.desc_sommets());
  mesh.preparer_tableau_avant_transport(vinterp_, mesh.desc_sommets());

  if (!mesh.Surfactant_facettes().get_disable_surfactant())
    {
      //supprimer_duplicata_bulles();
      //maillage_ft_ijk_.parcourir_maillage();
      FT_Field& Surfactant = mesh.Surfactant_facettes_non_const();
      ArrOfDouble surfactant_avant_remaillage = Surfactant.check_conservation(mesh);
      Surfactant.passer_variable_intensive(mesh);

      mesh.transporter(deplacement);
      double dt = dt_tot;
      if (rk_step >= 0)
        {
          dt = compute_fractionnal_timestep_rk3(dt_tot, rk_step);
        }
      Surfactant.avancer_en_temps(mesh, dt);
      Surfactant.passer_variable_extensive(mesh);
      ArrOfDouble surfactant_apres_remaillage = Surfactant.check_conservation(mesh);
      Surfactant.correction_conservation_globale(mesh,  surfactant_avant_remaillage,  surfactant_apres_remaillage);
      //transferer_bulle_perio();
      //creer_duplicata_bulles();
      //maillage_ft_ijk_.parcourir_maillage();
    }
  else
    {
      mesh.transporter(deplacement);
    }


  nbsom = sommets.dimension(0); // Le deplacement a peut-etre cree des sommets...
  mesh.update_tableau_apres_transport(coord_sommets_avant_deplacement, nbsom, mesh.desc_sommets());
  mesh.update_tableau_apres_transport(vinterp_, nbsom, mesh.desc_sommets());

  if (rk_step >= 0)
    {
      // on fait du RK3 :
      mesh.update_tableau_apres_transport(RK3_G_store_vi_, nbsom, mesh.desc_sommets());
      // L'update fini par l'echange EV
      // mesh.desc_sommets().echange_espace_virtuel(RK3_G_store_vi_);
    }

  if (!mesh.Surfactant_facettes().get_disable_surfactant())
    {
      //mesh.nettoyer_maillage();
      //supprimer_duplicata_bulles();
      //mesh.parcourir_maillage();
      nbsom = sommets.dimension(0); // Le deplacement a peut-etre cree des sommets...
    }

  // Tableau aux sommets :
  var_volume_deformation_.resize(nbsom);
  var_volume_deformation_ = 0.;
  remaillage_ft_ijk_.calculer_variation_volume(mesh, coord_sommets_avant_deplacement, var_volume_deformation_);

  // Comportement par defaut (sans la correction semi-locale du volume) :
  // L'effet du transport sur le volume des bulles est ajoute a dvol, ce qui sera pris en compte lors de la correction globale du volume (au remaillage)
  // Avec la correction semi-locale du volume, seul la variation par sommet est utile (var_volume_deformation_), et ce bloc n'est donc pas execute.
  if (!correction_semi_locale_volume_bulle)
    {
      // Tableau par compo connexe (integrale sur chaque bulle de la var volume au
      // cours du sous pas de temps : (la variation totale integree au cours du pas
      // de temps est dans dvol)
      ArrOfDouble var_volume_par_bulle(nbulles_tot);

      // On recalcule certaines grandeurs apres le deplacement :
      mesh.calculer_compo_connexe_sommets(compo_connex_som);

      for (int isom = 0; isom < nbsom; isom++)
        {
          // Ne prendre que les sommets reels (sinon on les compte plusieurs fois)
          if (mesh.sommet_virtuel(isom))
            continue;

          int compo = compo_connex_som[isom];
          // Les bulles dupliquees sont a la fin
          // La premiere (numero -1) est a la position nbulles_reelles :
          if (compo < 0)
            {
              // L'index de la bulle ghost est (entre -1 et -nbulles_ghost):
              const int idx_ghost = get_ghost_number_from_compo(compo);
              // On la place en fin de tableau :
              compo = nbulles_reelles - 1 - idx_ghost;
            }
          assert(compo >= 0);
          assert(compo < nbulles_tot);

          const double v = var_volume_deformation_[isom];
          var_volume_par_bulle[compo] += v;
        }
      mp_sum_for_each_item(var_volume_par_bulle);
      // Mise a jour de la variation totale de volume au cours du pas de temps :
      // (necessaire meme en euler. Si euler, on y passe qu'une fois)
      //  if (rk_step >=0) {
      for (int icompo = 0; icompo < nbulles_tot; icompo++)
        {
          dvol[icompo] += var_volume_par_bulle[icompo];
        }
      //  }
    }
}

// La fonction transporter_maillage a ete separee en trois pour inserer un
// calcul de l'indicatrice entre les differentes etapes.
// Cette fonction realise le remaillage uniquement. Elle suppose que la
// fonction transporter_maillage_deformation a ete precedemment appelee pour
// realiser le deplacement et remplir dvol en entree.
void IJK_Interfaces::transporter_maillage_remaillage(int correction_semi_locale_volume_bulle,
                                                     const DoubleTab& vitesses_translation_bulles,
                                                     const DoubleTab& mean_bubble_rotation_vector,
                                                     const DoubleTab& centre_gravite,
                                                     double dt_tot,
                                                     ArrOfDouble& dvol,
                                                     const int rk_step = -1,
                                                     const double temps = -1 /*pas de remaillage*/)
{
  Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const DoubleTab& sommets = mesh.sommets(); // Tableau des coordonnees des marqueurs.
  int nbsom = sommets.dimension(0);

  if (correction_semi_locale_volume_bulle)
    {
      const double fractionnal_timestep = (rk_step == -1) ? dt_tot : compute_fractionnal_timestep_rk3(dt_tot, rk_step);
      calculer_var_volume_remaillage(fractionnal_timestep, vitesses_translation_bulles, mean_bubble_rotation_vector, centre_gravite, var_volume_remaillage_);

      // Peut-etre que des sommets virtuels ont ete ajoutees (par le parcours de l'interface), mais je pense que le nombre de sommets reels n'a pas change.
      // Mise-a-jour du tableau var_volume_deformation_ pour cette eventualite :
      var_volume_deformation_.resize(nbsom);
      maillage_ft_ijk_.desc_sommets().echange_espace_virtuel(var_volume_deformation_);

      var_volume_remaillage_ -= var_volume_deformation_;
    }
  else
    {
      var_volume_remaillage_.resize(nbsom);
      var_volume_remaillage_ = 0.;
    }

  var_volume_correction_globale_.resize(nbsom);
  var_volume_correction_globale_ = 0.;

  ArrOfDouble surface_par_bulle;

  const int nbulles_reelles = get_nb_bulles_reelles();

  const ArrOfInt& compo_connex_apres_transport = mesh.compo_connexe_facettes();
  calculer_surface_bulles(surface_par_bulle);
  const ArrOfDouble& surface_facette_apres_transport = mesh.get_update_surface_facettes();

  // Calculer une variation de volume a imposer sur chaque sommet pour conserver
  // le volume des bulles (necessaire si euler, cad rk_step = -1 ou si dernier
  // sous pas de rk3) On ne fait la correction qu'au dernier sous pas de temps
  // car on ne veut pas appliquer le lissage
  //   au cours des sous pas de temps intermediaires car cela detruirait
  //   RK3_G_store_vi_. En effet, le lissage barycentrage deplace les sommets,
  //   donc en parallele il y a des sommets qui changent de proc, donc creation
  //   destruction d'entrees dans les tableaux
  if ((rk_step == -1) || (rk_step == 2))
    {
      const int nb_facettes = mesh.nb_facettes();
      const IntTab& facettes = mesh.facettes();
      for (int i = 0; i < nb_facettes; i++)
        {
          // Ignorer les facettes virtuelles
          if (mesh.facette_virtuelle(i))
            continue;
          int compo = compo_connex_apres_transport[i];
          // Les bulles dupliquees sont a la fin
          if (compo < 0)
            {
              // L'index de la bulle ghost est (entre -1 et -nbulles_ghost):
              const int idx_ghost = get_ghost_number_from_compo(compo);
              // On la place en fin de tableau :
              compo = nbulles_reelles - 1 - idx_ghost;
            }
          const double var_volume_tot = dvol[compo];
          const double surface_bulle = surface_par_bulle[compo];
          const double s = surface_facette_apres_transport[i];
          const double dv = var_volume_tot * s / (surface_bulle * 3.);
          for (int j = 0; j < 3; j++)
            {
              const int isom = facettes(i, j);
              var_volume_correction_globale_[isom] -= dv; // Signe negatif car on veut annuler la variation de volume
              // engendree par le transport
            }
        }
      // Sommer les contributions mises sur les sommets virtuels:
      mesh.desc_sommets().collecter_espace_virtuel(var_volume_correction_globale_, MD_Vector_tools::EV_SOMME);
      mesh.desc_sommets().echange_espace_virtuel(var_volume_correction_globale_);
    }

  RK3_G_store_vi_.resize(mesh.nb_sommets(), 3);

  var_volume_remaillage_ += var_volume_correction_globale_;

  // Au dernier pas de temps rk3 ou au pas de temps euler : on lisse et on
  // nettoie :
  if ((rk_step == -1) || (rk_step == 2))
    {
      // ca passe meme en laissant les ghost_bulles lors du lissage...car le
      // Domaine_VF sous-jacent est sur DOM_EXT.
      //   Attention : s'il faut les supprimer, le tableau var_volume_remaillage_ n'est alors
      //   plus valide (nombre de sommet a change).
      remailler_interface(temps, mesh, var_volume_remaillage_, remaillage_ft_ijk_);
      // remaillage_ft_ijk_.barycentrer_lisser_systematique_ijk(mesh, var_volume_remaillage_);
      mesh.nettoyer_maillage();
      // A partir d'ici, le tableau n'est plus valide donc on le reduit a 0 :
      RK3_G_store_vi_.resize(0, 3);
    }

  if (!mesh.Surfactant_facettes().get_disable_surfactant())
    {
      //Apres remaillage, recalculer les gradient, laplacien FT
      //mesh.nettoyer_maillage();
      //creer_duplicata_bulles();
      //mesh.parcourir_maillage();
      mesh.update_gradient_laplacien_Surfactant();
      mesh.update_sigma_and_interfacial_source_term_sommet(ref_domaine_, false, use_tryggvason_interfacial_source_);

    }



  const int nbsom_end_transport = mesh.nb_sommets();
  const int size_store = RK3_G_store_vi_.dimension(0);
  if (!((nbsom_end_transport == size_store) || (0 == size_store)))
    {
      Cerr << "Une des tailles de tableau n'est pas bonne... "
           << " size_store = " << size_store
           << " nbsom_end_transport = " << nbsom_end_transport
           << finl;
      Process::exit();
    }

  assert((mesh.nb_sommets() == RK3_G_store_vi_.dimension(0)) || (0 == RK3_G_store_vi_.dimension(0)));
}

// La fonction transporter_maillage a ete separee en trois pour inserer un
// calcul de l'indicatrice entre les differentes etapes.
// Cette fonction realise le deplacement rigide (indeformable) uniquement.
// Il est normalement effectue en dernier, apres la deformation et le remaillage
// de la bulle.
void IJK_Interfaces::transporter_maillage_rigide(const double dt_tot,
                                                 const DoubleTab& vitesses_translation_bulles,
                                                 const DoubleTab& mean_bubble_rotation_vector,
                                                 const DoubleTab& centre_gravite,
                                                 const int rk_step,
                                                 const int first_step_interface_smoothing)
{
  // nouvelle version:
  Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const DoubleTab& sommets = mesh.sommets(); // Tableau des coordonnees des marqueurs.
  int nbsom = sommets.dimension(0);
  DoubleTab deplacement(nbsom, 3);

  const int nbulles_reelles = get_nb_bulles_reelles();

  ArrOfIntFT compo_connex_som;
  mesh.calculer_compo_connexe_sommets(compo_connex_som);

  for (int som = 0; som < nbsom; som++)
    {
      if (mesh.sommet_virtuel(som))
        {
          // Valeur pipo pour dire qu'on n'initialise pas
          deplacement(som, 0) = 100.;
          deplacement(som, 1) = 100.;
          deplacement(som, 2) = 100.;
          continue;
        }

      int icompo = compo_connex_som[som];
      // Les bulles dupliquees sont a la fin
      // La premiere (numero -1) est a la position nbulles_reelles :
      if (icompo < 0)
        {
          // L'index de la bulle ghost est (entre -1 et -nbulles_ghost):
          const int idx_ghost = get_ghost_number_from_compo(icompo);
          // On la place en fin de tableau :
          icompo = nbulles_reelles - 1 - idx_ghost;
        }
      assert(icompo >= 0);
      assert(icompo < nbulles_reelles + get_nb_bulles_ghost());

      double rot_x = mean_bubble_rotation_vector(icompo, 0);
      double rot_y = mean_bubble_rotation_vector(icompo, 1);
      double rot_z = mean_bubble_rotation_vector(icompo, 2);
      Vecteur3 rot = {rot_x, rot_y, rot_z};

      double x = sommets(som, 0);
      double y = sommets(som, 1);
      double z = sommets(som, 2);
      Vecteur3 coord_sommet = {x, y, z};

      Vecteur3 coord_centre = {centre_gravite(icompo,0), centre_gravite(icompo,1), centre_gravite(icompo,2)};

      Vecteur3 radial_position = coord_sommet - coord_centre;

      Vecteur3 tangential_velocity;
      Vecteur3::produit_vectoriel(rot, radial_position, tangential_velocity);

      for (int direction = 0; direction < 3; direction++)
        {
          const double fractionnal_timestep = (rk_step == -1) ? dt_tot : compute_fractionnal_timestep_rk3(dt_tot, rk_step);
          deplacement(som, direction) = (vitesses_translation_bulles(icompo, direction) + tangential_velocity[direction]) * fractionnal_timestep;
        }
    }

  mesh.transporter(deplacement);
}

// Voir Topologie_Maillage_FT::remailler_interface  pour l'original.
void IJK_Interfaces::remailler_interface(const double temps,
                                         Maillage_FT_IJK& maillage,
                                         ArrOfDouble& var_volume,
                                         Remaillage_FT_IJK& algo_remaillage_local)
{
  // Le remaillage qui suit le traitement des coalescences peut a nouveau
  //  produire un maillage impropre (avec des facettes qui se coupent).
  // Dans ce cas, on s'arrete apres le traitement des coalescences.

  // Note B.M. : le remaillag systematique peut retarder la necessite de
  // remailler localement (apparition de grandes ou petites aretes). Donc
  // je le fais d'abord, et ensuite je teste s'il faut faire un remaillage
  // local.
  // L'intervalle de temps entre deux lissages est-il ecoule ?
  //  if (algo_remaillage_local.a_lisser(temps))
  //    {
  // On a choisi de lisser systematique :
  if (!maillage_ft_ijk_.Surfactant_facettes().get_only_remaillage())
    algo_remaillage_local.barycentrer_lisser_systematique_ijk(maillage, var_volume);

  // L'intervalle de temps entre deux remaillages locaux est-il ecoule:
  if (algo_remaillage_local.a_remailler(temps, maillage))
    {
      // Declanchement d'un remaillage local.
      // Pour que ca marche bien, les parametres de barycentrage et lissage
      // "apres_remaillage" doivent etre suffisants pour ramener le maillage dans
      // un etat correct sans appliquer de lissage systematique apres
      algo_remaillage_local.remaillage_local_interface(temps, maillage);
    }
}

void IJK_Interfaces::calculer_vitesse_de_deformation(
  int compo,
  const DoubleTab& bounding_box_bulles,
  const Cut_field_vector3_double& cut_field_velocity,
  const DoubleTab& vitesses_translation_bulles,
  const DoubleTab& mean_bubble_rotation_vector,
  const DoubleTab& positions_bulles)
{
  Cut_field_vector3_double& cut_field_deformation_velocity = static_cast<Cut_field_vector3_double&>(deformation_velocity_);
  assert(&cut_field_deformation_velocity[0].get_cut_cell_disc() == &cut_field_deformation_velocity[1].get_cut_cell_disc());
  assert(&cut_field_deformation_velocity[0].get_cut_cell_disc() == &cut_field_deformation_velocity[2].get_cut_cell_disc());
  const Cut_cell_FT_Disc& cut_cell_disc = cut_field_deformation_velocity[0].get_cut_cell_disc();

  const int ni = cut_field_deformation_velocity[0].ni();
  const int nj = cut_field_deformation_velocity[0].nj();
  const int nk = cut_field_deformation_velocity[0].nk();
  const int ghost = cut_field_deformation_velocity[0].ghost();

  const Domaine_IJK& geom = cut_field_deformation_velocity[0].get_domaine();
  assert(geom.is_uniform(0));
  assert(geom.is_uniform(1));
  assert(geom.is_uniform(2));

  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);

  double origin_x = geom.get_origin(DIRECTION_I);
  double origin_y = geom.get_origin(DIRECTION_J);
  double origin_z = geom.get_origin(DIRECTION_K);

  const int offset_x = geom.get_offset_local(DIRECTION_I);
  const int offset_y = geom.get_offset_local(DIRECTION_J);
  const int offset_z = geom.get_offset_local(DIRECTION_K);

  int imin = std::max(-ghost,     (int)((bounding_box_bulles(compo, 0, 0) - origin_x)/dx - offset_x - 2));
  int imax = std::min(ni + ghost, (int)((bounding_box_bulles(compo, 0, 1) - origin_x)/dx - offset_x + 2));
  int jmin = std::max(-ghost,     (int)((bounding_box_bulles(compo, 1, 0) - origin_y)/dy - offset_y - 2));
  int jmax = std::min(nj + ghost, (int)((bounding_box_bulles(compo, 1, 1) - origin_y)/dy - offset_y + 2));
  int kmin = std::max(-ghost,     (int)((bounding_box_bulles(compo, 2, 0) - origin_z)/dz - offset_z - 2));
  int kmax = std::min(nk + ghost, (int)((bounding_box_bulles(compo, 2, 1) - origin_z)/dz - offset_z + 2));
  for (int k = kmin; k < kmax; k++)
    {
      for (int j = jmin; j < jmax; j++)
        {
          for (int i = imin; i < imax; i++)
            {
              const double x_centre_cell = (i + offset_x + .5)*dx + origin_x;
              const double y_centre_cell = (j + offset_y + .5)*dy + origin_y;
              const double z_centre_cell = (k + offset_z + .5)*dz + origin_z;

              double rot_x = mean_bubble_rotation_vector(compo, 0);
              double rot_y = mean_bubble_rotation_vector(compo, 1);
              double rot_z = mean_bubble_rotation_vector(compo, 2);
              Vecteur3 rot = {rot_x, rot_y, rot_z};

              int n = cut_cell_disc.get_n(i,j,k);
              if (n >= 0)
                {
                  for (int phase = 0 ; phase < 2 ; phase++)
                    {
                      double x_cut_cell = x_centre_cell - .5*dx + dx*cut_cell_disc.get_interfaces().get_barycentre(1, 0, phase, i,j,k);
                      double y_cut_cell = y_centre_cell - .5*dy + dy*cut_cell_disc.get_interfaces().get_barycentre(1, 1, phase, i,j,k);
                      double z_cut_cell = z_centre_cell - .5*dz + dz*cut_cell_disc.get_interfaces().get_barycentre(1, 2, phase, i,j,k);

                      Vecteur3 coord = {x_cut_cell, y_cut_cell, z_cut_cell};

                      Vecteur3 coord_centre = {positions_bulles(compo,0), positions_bulles(compo,1), positions_bulles(compo,2)};

                      Vecteur3 radial_position = coord - coord_centre;

                      Vecteur3 tangential_velocity;
                      Vecteur3::produit_vectoriel(rot, radial_position, tangential_velocity);

                      for (int direction = 0; direction < 3; direction++)
                        {
                          const DoubleTabFT_cut_cell& diph_velocity = (phase == 0) ? cut_field_velocity[direction].diph_v_ : cut_field_velocity[direction].diph_l_;
                          DoubleTabFT_cut_cell& deformation_diph_velocity = (phase == 0) ? cut_field_deformation_velocity[direction].diph_v_ : cut_field_deformation_velocity[direction].diph_l_;
                          deformation_diph_velocity(n) = diph_velocity(n) - vitesses_translation_bulles(compo, direction) - tangential_velocity[direction];
                        }
                    }
                }
              else
                {
                  Vecteur3 coord = {x_centre_cell, y_centre_cell, z_centre_cell};

                  Vecteur3 coord_centre = {positions_bulles(compo,0), positions_bulles(compo,1), positions_bulles(compo,2)};

                  Vecteur3 radial_position = coord - coord_centre;

                  Vecteur3 tangential_velocity;
                  Vecteur3::produit_vectoriel(rot, radial_position, tangential_velocity);

                  for (int direction = 0; direction < 3; direction++)
                    {
                      cut_field_deformation_velocity[direction].pure_(i,j,k) = cut_field_velocity[direction].pure_(i,j,k) - vitesses_translation_bulles(compo, direction) - tangential_velocity[direction];
                    }
                }
            }
        }
    }
}

// bounding_box(b,dir,m) :
//      b -> Numero de la composante connexe de la bulle.
//      dir -> Direction i,j ou k.
//      m   -> min (0) ou max (1)
void IJK_Interfaces::calculer_bounding_box_bulles(DoubleTab& bounding_box, int option_shear) const
{
  const int nbulles = get_nb_bulles_reelles();
  bounding_box.resize(nbulles, 3, 2);

  int direction;
  // Initialisation avec des valeurs inatteignables :
  const double inval = 1.e20;
  bounding_box = inval;

  // Puis on parcours le maillage :
  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  ArrOfIntFT compo_connex_som;
  mesh.calculer_compo_connexe_sommets(compo_connex_som);
  const DoubleTab& sommets = mesh.sommets(); // Tableau des coordonnees des marqueurs.
  const int nbsom = sommets.dimension(0);
  ArrOfDouble volume_reel;
  DoubleTab position;
  calculer_volume_bulles(volume_reel, position);

  ArrOfDouble position_xmax_compo;
  ArrOfDouble position_xmin_compo;
  if (option_shear != 0 && IJK_Shear_Periodic_helpler::defilement_ == 1)
    {
      position_xmax_compo.resize_array(nbulles, RESIZE_OPTIONS::NOCOPY_NOINIT);
      position_xmin_compo.resize_array(nbulles, RESIZE_OPTIONS::NOCOPY_NOINIT);
      position_xmax_compo = -10000.;
      position_xmin_compo = 10000.;

      for (int i_sommet2 = 0; i_sommet2 < nbsom; i_sommet2++)
        {
          int iconnex = compo_connex_som[i_sommet2];
          double coord = sommets(i_sommet2, 0);
          if (coord>position_xmax_compo[iconnex])
            position_xmax_compo[iconnex] = coord;
          if (coord<position_xmin_compo[iconnex])
            position_xmin_compo[iconnex] = coord;
        }
      mp_min_for_each_item(position_xmin_compo);
      mp_max_for_each_item(position_xmax_compo);
    }

  for (int i_sommet = 0; i_sommet < nbsom; i_sommet++)
    {
      for (direction=0; direction<3; direction++)
        {
          double coord = sommets(i_sommet, direction);
          int iconnex = compo_connex_som[i_sommet];

          if (direction==0 && option_shear != 0 && IJK_Shear_Periodic_helpler::defilement_ == 1)
            {
              // position du barycentre de la bulle de reference a laquelle appartient le sommet
              double pos_ref = position(iconnex,0);
              //const Domaine_IJK& split = ref_domaine_.valeur();
              double Lx =  IJK_Shear_Periodic_helpler::Lx_for_shear_perio;
              //double Lx =  split.get_domain_length(0) - (position_xmax_compo(iconnex)-pos_ref);
              double offset = option_shear * IJK_Shear_Periodic_helpler::shear_x_time_;
              // le barycentre de la bulle reelle (compo >0) est situe entre db et Lx + db (pas entre 0 et Lx)
              // vrai uniquement pour des bulles qui montent.
              // Ne fonctionnera pas pour des bulles descendantes...
              // Idem donc pour la bulle ghost
              double decallage_bulle_reel_ext_domaine_reel = 1.*(position_xmax_compo[iconnex]-position_xmin_compo[iconnex]); // verifier si cest ca la valeur
              // assure une position entre db et Lx + db --> vrai que pour une bulle qui monte....
              // Si la bulle descend, risque de ne pas fonctionner --> il faudrait assurer une position entre -db et Lx+db ?
              double pos = std::fmod(std::fmod(pos_ref + offset - decallage_bulle_reel_ext_domaine_reel, Lx) + Lx, Lx) + decallage_bulle_reel_ext_domaine_reel;
              // Tous les sommets d'une meme bulle deplaces de la meme maniere sur x
              coord += (pos - pos_ref);
            }



          if (iconnex >= 0)
            {
              double mmin = bounding_box(iconnex, direction, 0);
              double mmax = -bounding_box(iconnex, direction, 1);
              if (coord < mmin)
                bounding_box(iconnex, direction, 0) = coord; // Un nouveau min est trouve.
              if (coord > mmax)
                bounding_box(iconnex, direction, 1) = -coord; // Un nouveau max est trouve. On stock son oppose
            }
          else
            {
              // le sommet n'appartient a aucune facette sur ce processeur, il est pris en compte sur un autre proc,
              // on ne fait rien.
            }
        }
    }
  // Communication inter-proc pour trouver les vrais limites :
  mp_min_for_each_item(bounding_box);
  for (int ibulle = 0; ibulle < nbulles; ibulle++)
    for (direction = 0; direction < 3; direction++)
      bounding_box(ibulle, direction, 1) = -bounding_box(ibulle, direction, 1); // a present on a le max.

}

// Methode pour creer les duplication de bulles lorsquelles sorte du
// domaine NS et entrent dans le domaine geom_FT...
void IJK_Interfaces::creer_duplicata_bulles()
{
  // Evaluation du cube contenant chaque bulle :
  DoubleTab bounding_box;
  calculer_bounding_box_bulles(bounding_box);
  ArrOfInt masque_duplicata_pour_compo_reel;

  // Pour chaque compo_connexe, remplir dans le tableau
  // masque_duplicata_pour_compo_reel un encodage du deplacement maximal pour toutes
  // les bulles reelles qui sortent de NS: Le critere pour declancher la duplication des
  // bulles est bounding_box_duplicate_criteria_
  // Pour les conditions shear-periodic, besoin d'un autre
  // masque_duplicata_pour_compo_ghost, un encodage du deplacement maximal pour toutes
  // les bulles duplique/decalle par le cisaillement qui sortent de NS: Le critere pour declancher la duplication des
  // bulles est le meme que pour les bulles reelles.
  if (IJK_Shear_Periodic_helpler::defilement_ == 1)
    {
      // Evaluation du cube contenant chaque bulle offset par le shear positif
      DoubleTab bounding_box_offsetp;
      calculer_bounding_box_bulles(bounding_box_offsetp, 1);
      // Evaluation du cube contenant chaque bulle offset par le shear negatif
      DoubleTab bounding_box_offsetm;
      calculer_bounding_box_bulles(bounding_box_offsetm, -1);
      preparer_duplicata_bulles(bounding_box, bounding_box_offsetp, bounding_box_offsetm, bounding_box_duplicate_criteria_, masque_duplicata_pour_compo_reel);
    }
  else
    {
      preparer_duplicata_bulles_masque_6bit(bounding_box, bounding_box_duplicate_criteria_, masque_duplicata_pour_compo_reel);
    }

  // Duplique et deplace les bulles de la liste :
  dupliquer_bulle_perio(masque_duplicata_pour_compo_reel);
}

// Input :
//   - code : contient soit seulement l'encodage du deplacement, soit aussi
//   celui du numero de bulle.
//            Dans les 2 cas, on ne retiendra dans code_deplacement que
//            l'encodage pure du deplacement.
// Retourne 0, -1 ou +1 selon le deplacement necessaire dans la direction dir :
static int decoder_deplacement(const int code, const int dir, int& compo_bulle_reel)
{
  // decodage du deplacement :
  const int code_deplacement = code & 63;        // 63 = 111111b est le  masque ne conservant que les 6 derniers bits.
  int tmp = code_deplacement & (1 << (3 + dir)); // l'operateur & permet de ne reveler
  // que le bit en face de 3+dir
  //                 cad le code (0 ou 1) pour le signe (resp. negatif ou
  //                 positif) dans la direction dir.
  // Si tmp == 0, c'est que le code signe pour la direction 'dir' etait 0, donc
  // que le deplacement est negatif.
  int signe = (tmp == 0) ? -1 /*deplacement negatif*/ : 1 /*deplacement positif*/;
  tmp = code_deplacement & (1 << dir); // retourne le code deplacement pour la direction consideree
  // (sur le 'dir'ieme bit)
  int index = tmp >> dir;              // ramene ce bit en derniere position.
  //         index = 0 si pas de deplacement, 1 sinon.


  compo_bulle_reel = code >> 6;
  // Cerr << "Deplacer bulle " << num_bulle << " Direction: " << dir
  //     << " Code_deplacement : " << code_deplacement
  //     << " Move : " << signe*index << finl;
  return signe * index;
}

// Le code pour le deplacement est code dans la compo connexe. Il faut le
// recuperer et le decoder. Le maillage transmis doit avoir son tableau des
// composantes connexes a jour.
static void calculer_deplacement_from_code_compo_connexe(const Maillage_FT_IJK& m, const Domaine_IJK& split,
                                                         DoubleTab& deplacement,
                                                         DoubleTab& bounding_box_NS, DoubleTab position,
                                                         const int nbulles, const Maillage_FT_IJK& mesh)
{
  // Creation du tableau deplacement pour le maillage m :
  const int nbsom = m.nb_sommets();
  deplacement.resize(nbsom, 3);
  // sommets du maillage temporaire (avec les ghosts a deplacer)
  ArrOfIntFT compo_sommets;
  m.calculer_compo_connexe_sommets(compo_sommets);

  // sommets du maillage initial (avec juste les bulles reelles)
  ArrOfIntFT compo_sommets_init;
  mesh.calculer_compo_connexe_sommets(compo_sommets_init);
  const int nbsomreel = mesh.nb_sommets();

  ArrOfDouble position_xmax_compo;
  ArrOfDouble position_xmin_compo;
  if (IJK_Shear_Periodic_helpler::defilement_ == 1)
    {
      position_xmax_compo.resize_array(nbulles, RESIZE_OPTIONS::NOCOPY_NOINIT);
      position_xmin_compo.resize_array(nbulles, RESIZE_OPTIONS::NOCOPY_NOINIT);
      position_xmax_compo = -10000.;
      position_xmin_compo = 10000.;

      const DoubleTab& sommets = mesh.sommets(); // Tableau des coordonnees des marqueurs.

      for (int i_sommet = 0; i_sommet < nbsomreel; i_sommet++)
        {
          // decodage du deplacement :
          //const int code = compo_sommets[i_sommet];
          int iconnex = compo_sommets_init[i_sommet];
          //std::cout << " coord = " << sommets(i_sommet, 0) << "compo_sommets" << iconnex << std::endl;
          double coord = sommets(i_sommet, 0);
          if (coord>position_xmax_compo[iconnex])
            position_xmax_compo[iconnex] = coord;
          if (coord<position_xmin_compo[iconnex])
            position_xmin_compo[iconnex] = coord;
        }
      mp_min_for_each_item(position_xmin_compo);
      mp_max_for_each_item(position_xmax_compo);
    }

  for (int dir = 0; dir < 3; dir++)
    {
      for (int i_sommet = 0; i_sommet < nbsom; i_sommet++)
        {
          // On relit le code pour le deplacement :
          const int code = compo_sommets[i_sommet];
          // On remplit le deplacement des sommets de la facette:

          // decodage du deplacement
          // le vrai compo_bulle_reel est lu par decoder_deplacement
          int compo_bulle_reel = 0;
          //int iconnex = compo_connex_som[i_sommet];
          int decode = decoder_deplacement(code, dir, compo_bulle_reel);

          double depl = decode * (bounding_box_NS(dir, 1) - bounding_box_NS(dir, 0));
          deplacement(i_sommet, dir) = depl;
          double pos_ref = 0 ;
          double pos = 0;
          double decallage_bulle_reel_ext_domaine_reel = 0.;
          // si seulement on a traverser une frontiere shear periodique
          if (dir==2 && depl != 0. && IJK_Shear_Periodic_helpler::defilement_ == 1)
            {
              double Lx =  IJK_Shear_Periodic_helpler::Lx_for_shear_perio;
              double offset = decode * IJK_Shear_Periodic_helpler::shear_x_time_;
              // position du barycentre de la bulle de reference a laquelle appartient le sommet
              pos_ref = position(compo_bulle_reel,0);
              // on veut le barycentre de la bulle decallee dans le domaine reel
              decallage_bulle_reel_ext_domaine_reel = 1.*(position_xmax_compo[compo_bulle_reel]-position_xmin_compo[compo_bulle_reel]); // verifier si cest ca la valeur
              pos = std::fmod(std::fmod(deplacement(i_sommet, 0) + pos_ref + offset - decallage_bulle_reel_ext_domaine_reel, Lx) + Lx, Lx) + decallage_bulle_reel_ext_domaine_reel;
              // Tous les sommets d'une meme bulle deplaces de la meme maniere sur x
              deplacement(i_sommet, 0) += (pos - pos_ref);
            }


        }
    }
}

// Le code pour le deplacement est code dans la compo connexe. Il faut le
// recuperer et le decoder. Le maillage transmis doit avoir son tableau des
// composantes connexes a jour. cette methode ..._negatif est a utiliser dans
// l'algo lorsque le code des bulles ghost est negatif
static void calculer_deplacement_from_code_compo_connexe_negatif(const Maillage_FT_IJK& m,
                                                                 DoubleTab& deplacement,
                                                                 DoubleTab& bounding_box_NS)
{
  // Creation du tableau deplacement pour le maillage m :
  const int nbsom = m.nb_sommets();
  deplacement.resize(nbsom, 3);
  ArrOfIntFT compo_sommets;
  m.calculer_compo_connexe_sommets(compo_sommets);

  for (int i_sommet = 0; i_sommet < nbsom; i_sommet++)
    {
      // On relit le code pour le deplacement :
      int code = compo_sommets[i_sommet];
      if (code >= 0)
        {
          // La bulle est reelle. Code est son numero. Il n'y a pas de deplacement a
          // prevoir. On force la reinit du tableau :
          for (int dir = 0; dir < 3; dir++)
            deplacement(i_sommet, dir) = 0.;
        }
      else
        {
          // La bulle est ghost. Il faut l'oppose de code pour decoder son
          // deplacement.
          code *= -1;
          // On remplit le deplacement des sommets de la facette:
          for (int dir = 0; dir < 3; dir++)
            {
              // decodage du deplacement :
              int unused_variable = 0 ;
              int decode = decoder_deplacement(code, dir, unused_variable);
              double depl = decode * (bounding_box_NS(dir, 1) - bounding_box_NS(dir, 0));
              deplacement(i_sommet, dir) = depl;
            }
        }
    }
}

// Le code ci-dessous est bugge car il peut y avoir des sommets  reel
// n'appartenant a aucune facette reelle. Une bulle passe sur un processeur
// voisin, du processeur 1 au processeur 2. Le premier sommet qui passe le joint
// est forcement sans facette car les facettes sont affectees au processeur de
// rang le plus bas qui possede un des sommets.
#if 0
// Input :
//   - m : Le maillage_ft doit contenir dans la compo_connexe les numeros de bulles reeles.
//   - masque_array[icompo] : Contient le code_deplacement de la ieme bulle.
static void calculer_deplacement_from_masque_in_array(const Maillage_FT_IJK& m,
                                                      DoubleTab& deplacement,
                                                      const ArrOfInt& masque_array,
                                                      DoubleTab& bounding_box_NS)
{
  // Creation du tableau deplacement pour le maillage m :
  const int nf = m.nb_facettes();
  const int nbsom = m.nb_sommets();
  const IntTab& mes_facettes = m.facettes();
  deplacement.resize(nbsom,3);
  const ArrOfInt& compo_connexe_tempo = m.compo_connexe_facettes();
  for (int i_facette = 0; i_facette < nf; i_facette++)
    {
      const int icompo = compo_connexe_tempo[i_facette];
      // On relit le code pour le deplacement :
      const int code = masque_array[icompo];
      // On remplit le deplacement des sommets de la facette:
      for (int dir = 0; dir < 3; dir++)
        {
          // decodage du deplacement :
          int unused_variable = 0 ;
          int decode = decoder_deplacement(code, dir, unused_variable);
          double depl = decode * (bounding_box_NS(dir,1) - bounding_box_NS(dir,0));
          for (int isom = 0; isom < 3; isom++)
            {
              const int idx_global_som = mes_facettes(i_facette, isom);
              deplacement(idx_global_som, dir) = depl;
            }
        }
    }
}
#else
// Input :
//   - m : Le maillage_ft doit contenir dans la compo_connexe les numeros de
//   bulles reeles.
//   - masque_array[icompo] : Contient le code_deplacement de la ieme bulle.
static void calculer_deplacement_from_masque_in_array(const Maillage_FT_IJK& m,
                                                      DoubleTab& deplacement,
                                                      const ArrOfInt& masque_array,
                                                      DoubleTab& bounding_box_NS, DoubleTab position, const int nbulles)
{
  // Creation du tableau deplacement pour le maillage m :
  const int nbsom = m.nb_sommets();
  deplacement.resize(nbsom,3);
  ArrOfIntFT compo_connexe_som;
  m.calculer_compo_connexe_sommets(compo_connexe_som);

  ArrOfDouble position_xmax_compo;
  ArrOfDouble position_xmin_compo;
  if (IJK_Shear_Periodic_helpler::defilement_ == 1)
    {
      position_xmax_compo.resize_array(nbulles, RESIZE_OPTIONS::NOCOPY_NOINIT);
      position_xmin_compo.resize_array(nbulles, RESIZE_OPTIONS::NOCOPY_NOINIT);
      position_xmax_compo = -10000.;
      position_xmin_compo = 10000.;

      const DoubleTab& sommets = m.sommets(); // Tableau des coordonnees des marqueurs.

      for (int i_sommet = 0; i_sommet < nbsom; i_sommet++)
        {
          int iconnex = compo_connexe_som[i_sommet];
          double coord = sommets(i_sommet, 0);
          if (coord>position_xmax_compo[iconnex])
            position_xmax_compo[iconnex] = coord;
          if (coord<position_xmin_compo[iconnex])
            position_xmin_compo[iconnex] = coord;
        }
      mp_min_for_each_item(position_xmin_compo);
      mp_max_for_each_item(position_xmax_compo);
    }

  for (int i_sommet = 0; i_sommet < nbsom; i_sommet++)
    {
      // On relit le code pour le deplacement :
      const int icompo = compo_connexe_som[i_sommet];
      const int code = masque_array[icompo];
      // On remplit le deplacement des sommets de la facette:
      for (int dir = 0; dir < 3; dir++)
        {
          // decodage du deplacement :
          int compo_bulle_reel = 0;
          //int iconnex = compo_connex_som[i_sommet];
          int decode = decoder_deplacement(code, dir, compo_bulle_reel);

          double depl = decode * (bounding_box_NS(dir, 1) - bounding_box_NS(dir, 0));
          deplacement(i_sommet, dir) = depl;
          double pos_ref = 0 ;
          double pos = 0;
          double decallage_bulle_reel_ext_domaine_reel = 0.;
          // si seulement on a traverser une frontiere shear periodique
          if (dir==2 && depl != 0. && IJK_Shear_Periodic_helpler::defilement_ == 1)
            {
              double Lx =  IJK_Shear_Periodic_helpler::Lx_for_shear_perio;
              double offset = decode * IJK_Shear_Periodic_helpler::shear_x_time_;
              // position du barycentre de la bulle de reference a laquelle appartient le sommet
              pos_ref = position(compo_bulle_reel,0);
              // on veut le barycentre de la bulle decallee dans le domaine reel
              decallage_bulle_reel_ext_domaine_reel = 1.*(position_xmax_compo[compo_bulle_reel]-position_xmin_compo[compo_bulle_reel]); // verifier si cest ca la valeur
              pos = std::fmod(std::fmod(deplacement(i_sommet, 0) + pos_ref + offset - decallage_bulle_reel_ext_domaine_reel, Lx) + Lx, Lx) + decallage_bulle_reel_ext_domaine_reel;

              deplacement(i_sommet, 0) += (pos - pos_ref);
            }


        }
    }
}
#endif

// Duplique les bulles a partir du masque_duplicata_pour_compo.
// Ce masque contient pour chaque bulle reelle un encodage du deplacement
// maximal a partir duquel :
//   o  On calcul toutes les copies a creer (jusqu'a 7 par bulle).
//   o  on ajoute les duplicatas dans un nouveau maillage,
//   o  on les deplace
//   o  on met -(code) dans la composante connexe
//             (ou code contient l'info du numero de bulle et du deplacement
//             realise a ce moment)
//   o  on ajoute les duplicatas au maillage actuel.
void IJK_Interfaces::dupliquer_bulle_perio(ArrOfInt& masque_duplicata_pour_compo)
{
  // Algorithme alternatif: optimise pour eviter de copier le maillage tout int.
  Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  Maillage_FT_IJK maillage_temporaire; // Maillage ou on va cumuler les bulles
  // dupliquees deplacees.
  maillage_temporaire.initialize(ref_domaine_.valeur(), refdomaine_dis_.valeur(), parcours_);
  // dupliquer tout le maillage:
  Maillage_FT_IJK a_dupliquer;
  // On recopie tout le maillage dans a_dupliquer:
  a_dupliquer.recopie(mesh, Maillage_FT_Disc::MINIMAL); // copie l'etat MINIMAL du maillage et
  // initialise la compo_connexe a partir
  // du mesh fourni.

  // Transforme le numero de composante en encodage numero composante +
  // deplacement:
  const ArrOfInt& compo_connexe_facettes = a_dupliquer.compo_connexe_facettes();
  int icompo, index, signe;
  for (int i_facette = 0; i_facette < a_dupliquer.nb_facettes(); i_facette++)
    {
      icompo = compo_connexe_facettes[i_facette];
      a_dupliquer.set_composante_connexe(i_facette, encoder_compo(icompo /*numero bulle*/, 0 /* pas de deplacement */));
    }
  // Boucle sur les duplications:
  // Exemple :
  // Composante connexe 0 : masque sans les bit de signes = 001 : sortie sur z --> 1 ghost (pas vrai en shear-perio)
  // Composante connexe 1 : masque sans les bit de signes = 010 : sortie sur y --> 1 ghost (pas vrai en shear-perio)
  // Composante connexe 2 : masque sans les bit de signes = 011 : sortie sur z et z --> 3 ghost (pas vrai en shear-perio)
  // pour shear_periodic_conditions : si sortie en z --> nb de ghost depend de la position du point periodique en face
  // Premiere iteration, je vais dupliquer les compo
  //  0 -> direction 001
  //  1 -> direction 010
  //  2 -> direction 001
  // Deuxieme iteration:
  //  2 -> direction 010
  // Troisieme iteration:
  //  2 -> direction 011
  ArrOfInt liste_bulles_crees;
  const int nbulles = get_nb_bulles_reelles();
  int nbulles_crees = 0;
  for (int mon_numero_iteration = 1; ; mon_numero_iteration++)
    {
      // Determine les composantes connexes a copier et quelle copie est a faire
      // lors de cette iteration:
      ArrOfInt index_copie(nbulles);
      ArrOfInt index_signe(nbulles);
      // On determine aussi les composantes connexe dont on a plus besoin
      // et on les supprime...
      ArrOfInt liste_facettes_pour_suppression;
      int reste_a_faire = 0;


      //ducluzeau : determination de index_copie a modifier (boucle ci-dessous) pour le shear-perio
      for (icompo = 0; icompo < nbulles; icompo++)
        {
          // Si on a epuise toutes les copies a faire pour cette composante,
          // index_copie restera a -1:
          index_copie[icompo] = -1;
          if (IJK_Shear_Periodic_helpler::defilement_ != 1)
            {
              int compteur = 0;
              const int masque = masque_duplicata_pour_compo[icompo] & 7; // masque sans les bits de signe
              for (index = 1; index <= 7; index++)
                {
                  if ((masque & index) == index)
                    compteur++; // Cette copie est a creer

                  if (compteur == mon_numero_iteration)
                    {
                      // Ok, c'est cette copie qu'il faut faire lors de cette iteration
                      nbulles_crees++;
                      index_copie[icompo] = index;
                      // Journal() << "IJK_Interfaces::dupliquer_bulle_perio : Duplication de la composante " << icompo << finl;
                      // Journal() << "vers le domaine FT index=" << index << " a l'iteration " << mon_numero_iteration << finl;
                      // Il reste au moins une bulle a copier donc on continuera la boucle...
                      reste_a_faire = 1;
                      break;
                    }
                }
            }
          else
            {
              index_signe[icompo] = -1;
              ArrOfInt index_bulle(7);
              ArrOfInt signe_bulle(7);
              // perio_xxx dit si il y periodicite ou non
              const int perio_x_reel = masque_duplicata_pour_compo[icompo] & 1;
              const int perio_y_reel = masque_duplicata_pour_compo[icompo] & (1 << 1);
              const int perio_z_reel = masque_duplicata_pour_compo[icompo] & (1 << 2);
              const int perio_x_ghost = masque_duplicata_pour_compo[icompo] & (1 << 3);
              // signe_xx dit si cette periodicite est de droite ou de gauche
//          int signe_4bit  = masque_duplicata_pour_compo[icompo] & (15 << 4);
              int signe_3bit  = (masque_duplicata_pour_compo[icompo] & (7 << 4)) >> 1;
//          int signe_x = masque_duplicata_pour_compo[icompo] & 16;
//          int signe_y = masque_duplicata_pour_compo[icompo] & (16 << 1);
//          int signe_z = masque_duplicata_pour_compo[icompo] & (16 << 2);

              int bit_position_x_signe_6bit = 3;
              int bit_position_x_ghost_8bit = 7;
              int signe_x_ghost = (masque_duplicata_pour_compo[icompo] & (1 << bit_position_x_ghost_8bit)) >> bit_position_x_ghost_8bit;
              int mask = 1<<bit_position_x_signe_6bit;

              for (int nb_bulle_duplique = 0; nb_bulle_duplique < 7; nb_bulle_duplique++)
                {
                  index_bulle[nb_bulle_duplique] = -1;
                  signe_bulle[nb_bulle_duplique] = signe_3bit;
                }

              // condition perio classique
              //index_reel = 1 --> direction x : bit 0001 --> 1 ghost
              //index_reel = 2 --> direction y : bit 0010 --> 1 ghost
              //index_reel = 3 --> direction x & y : bit 0011 --> 3 ghost
              //index_reel = 4 --> direction z : bit 0100 --> 1 ghost
              //index_reel = 5 --> direction z & x : bit 0101 --> 2 ghost  (au lieu de 3 en perio classique)
              //index_reel = 6 --> direction z & y : bit 0110 --> 3 ghost
              //index_reel = 7 --> direction z & y & x : bit 0111 --> 5 ghost (au lieu de 7 en perio classique)

              // condition supplementaire pour perio shear
              //index_reel = 8 --> impossible (pas de ghost sans direction z) : bit 1000
              //index_reel = 9 --> impossible (pas de ghost sans direction z) : bit 1001
              //index_reel = 10 --> impossible (pas de ghost sans direction z): bit 1010
              //index_reel = 11 --> impossible (pas de ghost sans direction z): bit 1011
              //index_reel = 12 --> direction z & x ghost : bit 1100 --> 2 ghost (et pas 3)
              //index_reel = 13 --> direction z & x & x ghost : bit 1101 --> 3 ghost (et pas 7)
              //index_reel = 14 --> direction z & y & x ghost : bit 1110 --> 5 ghost (et pas 7)
              //index_reel = 15 --> direction z & y & x & x ghost : bit 1111 --> 7 ghost (et pas 14)

              if (perio_z_reel == 4) // possible shear-perio a gerer, pas le meme nombre de ghost a droite et a gauche
                {
                  if (perio_y_reel == 2 && perio_x_reel == 1 && perio_x_ghost == 8 )
                    {
                      // l index bulle na plus que 3 bit lui (pour les 3 dimensions despace)
                      index_bulle[0] = 1;
                      index_bulle[1] = 2;
                      index_bulle[2] = 3;
                      index_bulle[3] = 4;
                      index_bulle[4] = 5;
                      signe_bulle[4] &= (~mask);
                      signe_bulle[4] |= (signe_x_ghost << bit_position_x_signe_6bit); // change le signe pour les bulles ghost de bulle ghost
                      index_bulle[5] = 6;
                      index_bulle[6] = 7;
                      signe_bulle[6] &= (~mask);
                      signe_bulle[6] |= (signe_x_ghost << bit_position_x_signe_6bit);
                    }
                  else if (perio_y_reel == 2 && perio_x_reel == 1 && perio_x_ghost != 8 )
                    {
                      index_bulle[0] = 1;
                      index_bulle[1] = 2;
                      index_bulle[2] = 3;
                      index_bulle[3] = 4;
                      index_bulle[4] = 6;
                    }
                  else if (perio_y_reel == 2 && perio_x_reel != 1 && perio_x_ghost == 8 )
                    {
                      index_bulle[0] = 5;
                      signe_bulle[0] &= (~mask);
                      signe_bulle[0] |= (signe_x_ghost << bit_position_x_signe_6bit);
                      index_bulle[1] = 6;
                      index_bulle[2] = 7;
                      signe_bulle[2] &= (~mask);
                      signe_bulle[2] |= (signe_x_ghost << bit_position_x_signe_6bit);
                      index_bulle[3] = 4;
                      index_bulle[4] = 2;
                    }
                  else if (perio_y_reel != 2 && perio_x_reel == 1 && perio_x_ghost == 8 )
                    {
                      index_bulle[0] = 1;
                      index_bulle[1] = 4;
                      index_bulle[2] = 5;
                      signe_bulle[2] &= (~mask);
                      signe_bulle[2] |= (signe_x_ghost << bit_position_x_signe_6bit);

                    }
                  else if (perio_y_reel != 2 && perio_x_reel != 1 && perio_x_ghost != 8 )
                    {
                      index_bulle[0] = 4;
                    }
                  else if (perio_y_reel != 2 && perio_x_reel != 1 && perio_x_ghost == 8 )
                    {
                      index_bulle[0] = 5;
                      signe_bulle[0] &= (~mask);
                      signe_bulle[0] |= (signe_x_ghost << bit_position_x_signe_6bit);
                      index_bulle[1] = 4;
                    }
                  else if (perio_y_reel != 2 && perio_x_reel == 1 && perio_x_ghost != 8 )
                    {
                      index_bulle[0] = 4;
                      index_bulle[1] = 1;
                    }
                  else if (perio_y_reel == 2 && perio_x_reel != 1 && perio_x_ghost != 8 )
                    {
                      index_bulle[0] = 4;
                      index_bulle[1] = 2;
                      index_bulle[2] = 6;
                    }
                  else
                    {
                      Cerr << "Erreur dans la gestion des bulles ghost" << finl;
                      Process::exit();
                    }
                }
              else // pas de probleme a gerer avec le shear-perio
                {
                  if (perio_y_reel == 2 && perio_x_reel == 1)
                    {
                      index_bulle[0] = 1;
                      index_bulle[1] = 2;
                      index_bulle[2] = 3;
                    }
                  else if (perio_y_reel == 2 && perio_x_reel != 1)
                    {
                      index_bulle[0] = 2;
                    }
                  else if (perio_y_reel != 2 && perio_x_reel == 1)
                    {
                      index_bulle[0] = 1;
                    }
                  else if (perio_y_reel != 2 && perio_x_reel != 1)
                    {
                      index_bulle[0] = -1;
                    }
                  else
                    {
                      Cerr << "Erreur dans la gestion des bulles ghost" << finl;
                      Process::exit();
                    }
                }

              if (index_bulle[mon_numero_iteration-1] != -1 && mon_numero_iteration-1 <= 6)
                {
                  // alors il reste au moins une bulle a cree pour cette compo dont on stocke le deplacement dans index_copie
                  index_copie[icompo] = index_bulle[mon_numero_iteration-1];
                  index_signe[icompo] = signe_bulle[mon_numero_iteration-1];
                  nbulles_crees++;
                  reste_a_faire = 1;
                }
            }
        }
      reste_a_faire = Process::mp_max(reste_a_faire);
      if (reste_a_faire == 0)
        {
          // plus aucune bulle a dupliquer, on sort de la boucle
          // mon_numero_iteration
          break;
        }
      // Supprime les facettes qui ne sont plus a dupliquer
      const int nf = a_dupliquer.nb_facettes();
      for (int i_facette = 0; i_facette < nf; i_facette++)
        {
          icompo = decoder_numero_bulle(compo_connexe_facettes[i_facette]);
          // Pour cette composante, quelle est la prochaine copie a faire ?
          index = index_copie[icompo];
          if(IJK_Shear_Periodic_helpler::defilement_ == 1)
            signe = index_signe[icompo];

          if (index > 0)
            {
              // Cette facette est a copier. on encode le deplacement dans la
              // composante connexe: xxxyyy -> Code binaire. x = 0 pour un deplacement
              // negatif, 1 pour un deplacement positif.
              //                         y = 0 : pas de deplacement dans cette
              //                         direction, 1 : deplacement.
              // code_deplacement : de 1 a 63 : mais il n'y a que 26 possibilite.
              //                    Pour chaque y = 0, x n'a pas d'influence; Par
              //                    convention, on met x=0
              // Exemple dans le plan xy :
              // ________________
              // | 21 | 20 | 27 |
              // |___|_____|____|
              // |   |     |    |
              // | 1 | NS  | 9  |
              // |___|_____|____|
              // | 3 |  2  | 11 |
              // |___|_____|____|
              // Calcul du deplacement a faire,
              if(IJK_Shear_Periodic_helpler::defilement_ != 1)
                signe = masque_duplicata_pour_compo[icompo] & (7 << 3); // Recupere seulement le signe.
              const int code_deplacement = signe | index;                 // l'index donne les directions a deplacer lors de
              // cette iteration.
              //                                    Le signe donne le sens.
              // On change la valeur de la composante connexe associee a la facette
              // pour memoriser le deplacement qu'elle subira ensuite.
              a_dupliquer.set_composante_connexe(i_facette, encoder_compo(icompo, code_deplacement));
            }
          else
            {
              if (index == 0)
                {
                  Cerr << "Comment est-ce possible? " << finl;
                  Process::exit();
                }
              // index est reste a sa valeur initiale de -1. Cette facette est donc a
              // effacer:
              liste_facettes_pour_suppression.append_array(i_facette);
            }
        }
      // La suppression gere les compo_connex (surcharge dans
      // Maillage_FT_IJK::nettoyer_maillage)
      a_dupliquer.supprimer_facettes(liste_facettes_pour_suppression);

      // La methode Maillage_FT_IJK::ajouter_maillage gere aussi les
      // compo_connexe_facettes
      maillage_temporaire.Surfactant_facettes_non_const().set_disable_surfactant(mesh.Surfactant_facettes().get_disable_surfactant());
      maillage_temporaire.ajouter_maillage_IJK(a_dupliquer);

      // fin de l'iteration mon_numero_iteration.
    }

  // ducluzeau
  // a ce niveau, maillage_temporaire contient toutes les facettes duplique, mais pas transportee.
  // Le numero de composante connexe de chaque facette contient l information sur son deplacement (encode)

  // Si le maillage temporaire contient des facettes, on les deplace :
  int nbf = maillage_temporaire.nb_facettes();
  const int max_nbf = Process::mp_max(nbf);
  if (max_nbf)
    {
      DoubleTab deplacement;
      // Le maillage_temporaire transmis a son tableau des composantes connexes a
      // jour. le tableau contient l'encodage pour le deplacement que l'on va
      // decoder :
      const Domaine_IJK& split = ref_domaine_.valeur();
      ArrOfDouble volume_reel;
      DoubleTab position;
      calculer_volume_bulles(volume_reel, position);

      calculer_deplacement_from_code_compo_connexe(maillage_temporaire, split,
                                                   deplacement,
                                                   bounding_box_NS_domain_, position, nbulles, mesh);
      // La methode transporter gere les compo connexes.
      maillage_temporaire.transporter(deplacement);
      maillage_temporaire.nettoyer_maillage();

      // Lors du transport, des sommets peuvent changer de processeur.
      // Le nombre de facette peut varier (eg de 0 a qqch...) :
      nbf = maillage_temporaire.nb_facettes();
      // forcer les composantes connexes du maillage temporaire a -code
      const ArrOfInt& compo_connexe_facettes_tempo = maillage_temporaire.compo_connexe_facettes();
      for (int iface = 0; iface < nbf; iface++)
        {
          const int code_negatif = -compo_connexe_facettes_tempo[iface];
          maillage_temporaire.set_composante_connexe(iface, code_negatif);
          liste_bulles_crees.append_array(code_negatif);
        }
      // ajouter_maillage_IJK gere la compo_connexe correctement.
      // On pourrait aussi appeler Maillage_FT_IJK::ajouter_maillage.
      mesh.ajouter_maillage_IJK(maillage_temporaire);
    }

  // On allege chaque liste, mais il faut de toute facon le refaire a la fin :
  array_trier_retirer_doublons(liste_bulles_crees);
  assert(ghost_compo_converter_.size_array() == nb_bulles_ghost_);
  // Si moi ou un autre proc a cree des bulles, il faut mettre a jour
  // nb_bulles_ghost_ et ghost_compo_converter_ :
  if (Process::mp_sum(nbulles_crees))
    {
      if (Process::je_suis_maitre())
        {
          ArrOfInt prov;
          int nbproc = Process::nproc();
          int taille;
          for (int p = 1; p < nbproc; p++)
            {
              recevoir(prov, p, 0, 2001 + p);
              taille = liste_bulles_crees.size_array();
              // On injecte le tableau prov a partir de la case taille :
              liste_bulles_crees.resize_array(taille + prov.size_array());
              liste_bulles_crees.inject_array(prov, -1 /* tout*/,
                                              taille /* first_element_dest */,
                                              0 /* first_element_source */);
            }

          // On retire les doublons (duplicatas trouves sur plusieurs procs...);
          array_trier_retirer_doublons(liste_bulles_crees);

          // On a cree des bulles ghost :
          nbulles_crees = liste_bulles_crees.size_array();
          nb_bulles_ghost_ += nbulles_crees;
          if (nb_bulles_ghost_ != nb_bulles_ghost_before_)
            recompute_indicator_ = 1;

          if (ghost_compo_converter_.size_array())
            {
              Cerr << " On vient de creer des bulles ghost alors qu'il y en avait "
                   "deja car le "
                   << " tableau ghost_compo_converter_ n'est pas de taille nulle! " << finl;
              Process::exit();
              // On append le tableau liste_bulles_crees a la fin de
              // ghost_compo_converter_ :
              const int nb_ghost_before = ghost_compo_converter_.size_array();
              ghost_compo_converter_.resize_array(nb_ghost_before + nbulles_crees);
              ghost_compo_converter_.inject_array(liste_bulles_crees, -1 /* tout*/, nb_ghost_before /* first_element_dest */,
                                                  0 /* first_element_source */);
            }
          else
            {
              ghost_compo_converter_.copy_array(liste_bulles_crees); // L'oparateur "=" existe aussi. Je ne sais pas
              // ce qui est mieux.
            }

#ifdef GB_VERBOSE
          // Tests et commentaires :
          Cerr << "IJK_Interfaces.cpp::dupliquer_bulle_perio liste_bulles_crees : "
               << liste_bulles_crees << finl;
          const int old_size = ghost_compo_converter_.size_array();
          ghost_compo_converter_.array_trier_retirer_doublons();
          const int new_size = ghost_compo_converter_.size_array();
          if (new_size != old_size)
            {
              Cerr << "On dirait que des doublons etaient presents dans la liste ghost_compo_converter_ ... "
                   << "Ca signifie peut-etre que l'on vient de creer des bulles qui existaient deja..."
                   << "C'est etrange... on s'arrete... ";
              Process::exit();
            }
#endif

        }
      else
        {
          // Tous les proc envoient au maitre :
          envoyer(liste_bulles_crees, Process::me(), 0, 2001 + Process::me());
        }

      // Le proc maitre a fini de remplir le tableau. Il communique a tous les
      // procs :
      envoyer_broadcast(nb_bulles_ghost_, 0);
      envoyer_broadcast(ghost_compo_converter_, 0);
      envoyer_broadcast(recompute_indicator_, 0);
    }

  assert(ghost_compo_converter_.size_array() == nb_bulles_ghost_);
}

// Methode pour remplir le masque codant les duplications a creer
// a partir des bounding box de chaques bulles.
// La methode dupliquer_bulle_perio construira les maillages et les vecteurs deplacement
// a partir du masque.
// Parametres :
//   - bounding_box : Tableau des coordonnees limites des cubes contenant chaques bulles.
//   - authorized_bounding_box : Tableau des coordonnees limites a ne pas depasser
//                               (eg, le domaine NS ou le domaine_(FT - ncells_forbidden_)...
// Resultat :
//   - masque_duplicata_pour_compo : Encodage du deplacement a appliquer aux bulles sortant
//                                   du domaine autorise.
void IJK_Interfaces::preparer_duplicata_bulles(const DoubleTab& bounding_box,
                                               const DoubleTab& bounding_box_offsetp,
                                               const DoubleTab& bounding_box_offsetm,
                                               const DoubleTab& authorized_bounding_box,
                                               ArrOfInt& masque_duplicata_pour_compo)
{

  if (Process::je_suis_maitre())
    {
      const int nbulles = get_nb_bulles_reelles();
      masque_duplicata_pour_compo.resize_array(nbulles);
      for (int icompo = 0; icompo < nbulles; icompo++)
        {
          int masque_sortie_domaine_reel = 0;
          // Masque sortie est un chiffre binaire qui dit dans quelles directions
          // la bulle sort du domaine (et le domaine est periodique)

          // pour le shear_periodic,
          // on passe d'une variable 6 bit, a 8 bit, puisqu on a un degres de liberte supplementaire sur les ghost
          for (int direction = 0; direction < 3; direction++)
            {
              if (perio_NS_[direction])
                {
                  // Est-ce qu'on sort par la gauche ?
                  if (bounding_box(icompo, direction, 0) < authorized_bounding_box(direction, 0))
                    {
                      masque_sortie_domaine_reel |= (1 << direction); // met le bit "direction" a 1 dans le masque
                      masque_sortie_domaine_reel |= (16 << direction); // met le bit de signe a 1 dans le masque
                      // il faudra deplacer la copie vers la droite
                      if(direction==2 && IJK_Shear_Periodic_helpler::defilement_ == 1)
                        // on est sorti en z, est-ce que la bulle ghost depasse en x ? pour condition perio shear
                        // si sortie de la bulle en z, verifier la sortie en x de la bulle ghost
                        // ici, sortie par la gauche, donc shear positif dans bounding_box_offsetp
                        {
                          if (bounding_box_offsetp(icompo, 0, 0) < authorized_bounding_box(0, 0))
                            {
                              masque_sortie_domaine_reel |= (1 << 3); // met le bit "direction" a 1 dans le masque
                              masque_sortie_domaine_reel |= (16 << 3); // met le bit de signe a 1 dans le masque
                            }
                          if (bounding_box_offsetp(icompo, 0, 1) > authorized_bounding_box(0, 1))
                            {
                              masque_sortie_domaine_reel |= (1 << 3); // met le bit "direction" a 1 dans le masque
                              // le bit de signe reste a zero, qui signifie un deplacement vers
                              // les coord negatives.
                            }
                        }
                    }
                  if (bounding_box(icompo, direction, 1) > authorized_bounding_box(direction, 1))
                    {
                      // On sort par la droite, il faudra deplacer la copie vers la
                      // gauche.
                      masque_sortie_domaine_reel |= (1 << direction); // met le bit "direction" a 1 dans le masque
                      // le bit de signe reste a zero, qui signifie un deplacement vers
                      // les coord negatives.
                      if(direction==2 && IJK_Shear_Periodic_helpler::defilement_ == 1)
                        // on est sorti en z, est-ce que la bulle ghost depasse en x ? pour condition perio shear
                        // si sortie de la bulle en z, verifier la sortie en x de la bulle ghost
                        // ici, sortie par la droite, donc shear negatif dans bounding_box_offsetm
                        {
                          if (bounding_box_offsetm(icompo, 0, 0) < authorized_bounding_box(0, 0))
                            {
                              masque_sortie_domaine_reel |= (1 << 3); // met le bit "direction" a 1 dans le masque
                              masque_sortie_domaine_reel |= (16 << 3); // met le bit de signe a 1 dans le masque
                            }
                          if (bounding_box_offsetm(icompo, 0, 1) > authorized_bounding_box(0, 1))
                            {
                              masque_sortie_domaine_reel |= (1 << 3); // met le bit "direction" a 1 dans le masque
                              // le bit de signe reste a zero, qui signifie un deplacement vers
                              // les coord negatives.
                            }
                        }
                    }
                }
            }

          // Stocke le masque dans le tableau resultat :
          masque_duplicata_pour_compo[icompo] = masque_sortie_domaine_reel;
          // duCluzeau : pour shear-perio, il faut stocker un autre tableau avec le nombre de duplicata par bulles.
          // pour ca, il faut tester la position de la bulle ghost en z, et verifier si elle sort en x
        }
    }
  envoyer_broadcast(masque_duplicata_pour_compo, 0);

}

void IJK_Interfaces::preparer_duplicata_bulles_masque_6bit(const DoubleTab& bounding_box,
                                                           const DoubleTab& authorized_bounding_box,
                                                           ArrOfInt& masque_duplicata_pour_compo)
{

  if (Process::je_suis_maitre())
    {
      const int nbulles = get_nb_bulles_reelles();
      masque_duplicata_pour_compo.resize_array(nbulles);
      for (int icompo = 0; icompo < nbulles; icompo++)
        {
          int masque_sortie_domaine_reel = 0;
          // Masque sortie est un chiffre binaire qui dit dans quelles directions
          // la bulle sort du domaine (et le domaine est periodique)

          // pour le shear_periodic,
          // on passe d'une variable 6 bit, a 8 bit, puisqu on a un degres de liberte supplementaire sur les ghost
          for (int direction = 0; direction < 3; direction++)
            {
              if (perio_NS_[direction])
                {
                  // Est-ce qu'on sort par la gauche ?
                  if (bounding_box(icompo, direction, 0) < authorized_bounding_box(direction, 0))
                    {
                      masque_sortie_domaine_reel |= (1 << direction); // met le bit "direction" a 1 dans le masque
                      masque_sortie_domaine_reel |= (8 << direction); // met le bit de signe a 1 dans le masque
                      // il faudra deplacer la copie vers la droite
                    }
                  if (bounding_box(icompo, direction, 1) > authorized_bounding_box(direction, 1))
                    {
                      masque_sortie_domaine_reel |= (1 << direction); // met le bit "direction" a 1 dans le masque
                    }
                }
            }

          // Stocke le masque dans le tableau resultat :
          masque_duplicata_pour_compo[icompo] = masque_sortie_domaine_reel;
        }
    }
  // Le proc maitre a fini de remplir le tableau de masque.
  // Il les communiques a tous les procs.
  // broadcast to all mpi processes:
  envoyer_broadcast(masque_duplicata_pour_compo, 0);
}

// Methode supprimant toutes les bulles dupliquees (reconnues par compo_connex <
// 0);
void IJK_Interfaces::supprimer_duplicata_bulles()
{
  Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const ArrOfInt& compo_connex = mesh.compo_connexe_facettes();
  const int nbfacettes = mesh.nb_facettes();

  // On parcours toutes les facettes du maillage.
  // On marque celles dont la compo_connexe est negative:
  ArrOfInt liste_facettes_pour_suppression;
  for (int iface = 0; iface < nbfacettes; iface++)
    {
      const int icompo = compo_connex[iface];
      //  if ((icompo < 0) or (icompo>70)) //  Pour ne recuperer que les 70
      //  premieres bulles.
      if (icompo < 0)
        liste_facettes_pour_suppression.append_array(iface);
    }

  mesh.supprimer_facettes(liste_facettes_pour_suppression);
  // A partir d'ici, le tableau n'est plus valide donc on le reduit a 0 :
  RK3_G_store_vi_.resize(0, 3);

  // On remet le nombre de bulles ghost a 0 :
  nb_bulles_ghost_before_ = nb_bulles_ghost_;
  nb_bulles_ghost_ = 0;
  ghost_compo_converter_.resize_array(0);
}

void IJK_Interfaces::transferer_bulle_perio()
{
  // Evaluation du cube contenant chaque bulle :
  DoubleTab bounding_box;
  calculer_bounding_box_bulles(bounding_box);

  // Pour chaque compo_connexe, remplir dans le tableau
  // masque_duplicata_pour_compo un encodage du deplacement maximal :
  ArrOfInt masque;
  // Recherche des bulles qui sortent du domaine authorise et remplit
  // dans le tableau masque un encodage du deplacement :
  preparer_duplicata_bulles_masque_6bit(bounding_box, bounding_box_forbidden_criteria_, masque);

  // Deplace les bulles de la liste :
  deplacer_bulle_perio(masque);

  if (follow_colors_ && Process::je_suis_maitre())
    {
      const int nbulles = get_nb_bulles_reelles();
      assert(through_yminus_.size_array() == nbulles);
      assert(perio_NS_[DIRECTION_J]);
      for (int icompo = 0; icompo < nbulles; icompo++)
        {
          const int masque_sortie = masque[icompo];
          if (masque_sortie & 2)
            {
              // La bulle icompo contient 1 sur le bit de deplacement selon y (bit 2^1
              // = 2) Il faut donc changer la valeur stockee dans through_yminus :
              if (masque_sortie & 16)
                {
                  // La bulle icompo contient 1 sur le bit de sign du deplacement selon
                  // y (bit 2^(3+1)=16) On sort par la gauche et la bulle va etre
                  // deplacer vers la droite :
                  //   donc on traverse yminus et rentre par la droite.
                  through_yminus_[icompo] = -1;
                }
              else
                {
                  // On traverse yplus et rentre par la gauche.
                  through_yminus_[icompo] = 1;
                }
            }
        }
    }
  // Le maitre informe les autres des nouvelles valeurs :
  if (follow_colors_)
    envoyer_broadcast(through_yminus_, 0);
}

// Deplace les bulles a travers la frontiere perio a partir du masque.
// Ce masque contient pour chaque bulle reelle un encodage du deplacement
// maximal.
//   - On marque les facettes a deplacer.
//   - On les transporte.
// ducluzeau : deplace la bulle reelle quand elle est sortie de domain_ns
void IJK_Interfaces::deplacer_bulle_perio(const ArrOfInt& masque_deplacement_par_compo)
{
  // ducluzeau : masque_deplacement_par_compo doit etre en 6 bit, pas 8 !
  Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  DoubleTab deplacement;
  const int nbulles = get_nb_bulles_reelles();
  ArrOfDouble volume_reel;
  DoubleTab position;
  calculer_volume_bulles(volume_reel, position);
  calculer_deplacement_from_masque_in_array(mesh, deplacement, masque_deplacement_par_compo, bounding_box_NS_domain_, position, nbulles);

  // On transporte le maillage :
  mesh.transporter(deplacement);
  // Un petit message si on transporte :

  for (int ib = 0; ib < nbulles; ib++)
    {
      const int code_deplacement = masque_deplacement_par_compo[ib];
      if (code_deplacement != 0)
        {
          Journal() << "IJK_Interfaces::deplacer_bulle_perio : Un deplacement a eu "
                    << "lieu pour la composante " << ib << finl;
          Journal() << "Deplacement x y z : ";
          for (int dir = 0; dir < 3; dir++)
            {
              int unused_variable = 0 ;
              const int decode = decoder_deplacement(code_deplacement, dir, unused_variable);
              Journal() << decode << " ";
            }
          Journal() << finl;

          recompute_indicator_ = 1;
          envoyer_broadcast(recompute_indicator_, 0);
        }
    }
}

// retourne le nombre de compo_connexe presentes dans l'element.
// Pour cela, la methode utilise les Intersections_Elem_Facettes_Data
// Or, il est possible qu'on n'ait des facettes reconnues comme intersectant un
// element avec une surface vide (en fait, elles ne le coupent pas). C'est une
// astuce du parcours pour ne pas les parcourir plusieurs fois... Voir
// Parcours_interface::calcul_intersection_facelem_3D pour explication :
//       // Intersection de surface nulle. On l'enregistre pour ne pas traiter
//       // cette facette a nouveau lors du parcours.
//      maillage.intersections_elem_facettes_.ajoute_intersection(num_facette,
//          num_element,
//          0.,
//          0.,
//          0., 0., 0.);
//
//

// The method have to receive the extended field indic_ft because
// the splitting and the conversion "num_elem = s.convert_ijk_cell_to_packed(i,
// j, k);" are required for num_compo_ which is on domaineVDF which is on the
// extended domain
void IJK_Interfaces::calculer_indicatrice(IJK_Field_double& indic)
{
  static Stat_Counter_Id calculer_indicatrice_counter_ =
    statistiques().new_counter(2, "Calcul rho mu indicatrice: calcul de l'indicatrice");
  statistiques().begin_count(calculer_indicatrice_counter_);

  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, refdomaine_dis_.valeur());
  const IntTab& elem_faces = domaine_vf.elem_faces();
  const IntTab& faces_voisins = domaine_vf.face_voisins();
  const Intersections_Elem_Facettes& intersec = maillage_ft_ijk_.intersections_elem_facettes();
  const Domaine_IJK& s = indic.get_domaine();

  indic.data() = 1.; // Tout liquide pour commencer
  num_compo_ = 0;    // Re-initialize the table...
  const int ni = indic.ni();
  const int nj = indic.nj();
  const int nk = indic.nk();

  // La methode de maillage_ft_disc utilise l'equation de transport pour
  // recuperer la fonction distance et un bricolage pour trouver l'indicatrice
  // au voisinage de l'interface. On va faire autrement. Debut identique a la
  // methode VDF: calcul des fractions de volume de phase 1 dans les mailles
  // traversees par des interfaces:
  {
    const ArrOfInt& index_elem = intersec.index_elem();
    //    const int nb_elem = index_elem.size_array();
    // Boucle sur les elements euleriens
    for (int k = 0; k < nk; k++)
      {
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                // Anciennement la methode etait portee par le mesh :
                //    const int num_elem =
                // maillage_ft_ijk_.convert_ijk_cell_to_packed(i, j, k);
                // A present, elle est dans le splitting :
                const int num_elem = s.convert_ijk_cell_to_packed(i, j, k);
                int index = index_elem[num_elem];
                double somme_contrib = 0.;
                // Boucle sur les facettes qui traversent cet element
                while (index >= 0)
                  {
                    const Intersections_Elem_Facettes_Data& data = intersec.data_intersection(index);
                    somme_contrib += data.contrib_volume_phase1_;
                    index = data.index_facette_suivante_;
                  };

                int nb_increment_somme_contrib = 0;
                while (somme_contrib > 1.)
                  {
                    somme_contrib -= 1.;
                    nb_increment_somme_contrib++;
                  }
                while (somme_contrib < 0.)
                  {
                    somme_contrib += 1.;
                    nb_increment_somme_contrib++;
                  }

                if (nb_increment_somme_contrib > 1)
                  {
                    Cerr << "nb_increment_somme_contrib= " << nb_increment_somme_contrib << finl;
                    // Process::exit("Error in IJK_Interfaces::calculer_indicatrice !");
                  }
                // GB Fix 2022: tolerance play:
                // Si l'on est proche de 0 ou de 1, on ne sait pas vraiment si on a bien fait nos calculs
                // (les modulos et les sommes peuvent avoir conduit a une imprecision).
                // On fait le choix de le considerer comme non-traverser, et on laisse le soin au calcul
                // des compos connexes de determiner si c'est 0 ou 1 :
                if ((somme_contrib < seuil_indicatrice_negligeable_) || ((1.-somme_contrib) < seuil_indicatrice_negligeable_))
                  {
                    somme_contrib = 0.; // L'elem retourne dans la liste des elements "non traverses"
                  }
                if (somme_contrib > 0.)
                  {
                    // if(somme_contrib == 1.)
                    //  somme_contrib -= 1e-3;

                    indic(i, j, k) = somme_contrib;
                    num_compo_[num_elem] = -1; // marque l'element pour apres
                  }
              }
          }
      }
  }

  static Stat_Counter_Id search_connex_components_counter_ =
    statistiques().new_counter(2, "Calcul de l'indicatrice : recherche compo connexes");
  statistiques().begin_count(search_connex_components_counter_);
  num_compo_.echange_espace_virtuel();
  // Recherche des composantes connexes sur le maillage eulerien
  int nb_compo_locales = search_connex_components_local(elem_faces, faces_voisins, num_compo_);
  nb_compo_in_num_compo_ = compute_global_connex_components(num_compo_, nb_compo_locales);
  // num_compo_.echange_espace_virtuel(); // On n'aura pas besoin de son ev a
  // jour.

  // Il y a au moins une phase continue :
  assert(nb_compo_in_num_compo_ - (nb_bulles_reelles_ + nb_bulles_ghost_) > 0);
  statistiques().end_count(search_connex_components_counter_);

  // Calcul des drapeaux selon la methode la plus robuste :
  ArrOfInt drapeau(nb_compo_in_num_compo_);
  compute_drapeaux_vapeur_v4(num_compo_, drapeau);

#ifdef GB_VERBOSE
  Cerr << "GB-check-interf nb_compo_in_num_compo_ - (nb_reelles + nb_ghost)= " << nb_compo_in_num_compo_ << " - "
       << nb_bulles_reelles_ + nb_bulles_ghost_ << " = "
       << nb_compo_in_num_compo_ - (nb_bulles_reelles_ + nb_bulles_ghost_) << finl;
  Cerr << "nb_compo_in_num_compo_= " << nb_compo_in_num_compo_ << finl;
  Cerr << "IJK_Interfaces::calculer_indicatrice Evaluations drapeau vapeur " << drapeau << finl;
#endif

#define VERIF_DRAPEAU 0
#if VERIF_DRAPEAU
  // Calcul des drapeaux selon les autres methodes :
  {
    // Pour debug :
    //      maillage_ft_ijk_.nettoyer_maillage();
    //      maillage_ft_ijk_.parcourir_maillage();

    const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
    ArrOfInt drapeau1(nb_compo_in_num_compo_);
    ArrOfInt drapeau1b(nb_compo_in_num_compo_);
    ArrOfInt drapeau2(nb_compo_in_num_compo_);
    ArrOfInt drapeau3(nb_compo_in_num_compo_);

    // marquer les composantes interieures aux bulles:
    compute_drapeaux_vapeur_v1(mesh, elem_faces, faces_voisins, num_compo, true, drapeau1);
    compute_drapeaux_vapeur_v1(mesh, elem_faces, faces_voisins, num_compo, false, drapeau1b);
    compute_drapeaux_vapeur_v2(num_compo, drapeau2);

    const Domaine_IJK& split = ref_domaine_.valeur();
    compute_drapeaux_vapeur_v3(mesh, split, num_compo, drapeau3);

    // HACKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
    drapeau = drapeau3;

    Cerr << "IJK_Interfaces::calculer_indicatrice Evaluations drapeaux" << finl;
    Cerr << "nb_compo_in_num_compo_= " << nb_compo_in_num_compo_ << finl;
    Cerr << "drapeau v1 towards_liquid=true  : " << drapeau1 << finl;
    Cerr << "drapeau v1 towards_liquid=false : " << drapeau1b << finl;
    Cerr << "drapeau v2 :       " << drapeau2 << finl;
    Cerr << "drapeau v3 :       " << drapeau3 << finl;
    Cerr << "drapeau v4 :       " << drapeau << finl;

    for (int idx = 0; idx < nb_compo_in_num_compo_; idx++)
      {
        int sum = drapeau[idx] + drapeau1[idx] + drapeau1b[idx] + drapeau2[idx] + drapeau3[idx];
        if ((drapeau[idx] != drapeau1[idx]) || (drapeau[idx] != drapeau1b[idx]) || (drapeau[idx] != drapeau2[idx]) ||
            (drapeau[idx] != drapeau3[idx]))
          Cerr << "Les methodes different... "
               << "valeurs donnees par les methodes : " << drapeau[idx] << drapeau1[idx] << drapeau1b[idx]
               << drapeau2[idx] << drapeau3[idx] << " pour idx = " << idx << finl;

        if ((sum != 0) && (sum != 5))
          {
            Cerr << "Les methodes different... "
                 << "sum = " << sum << " pour idx = " << idx << finl;
            // Process::exit();
          }
      }
    assert(drapeau == drapeau1);
    assert(drapeau == drapeau1b);
    assert(drapeau == drapeau2);
    assert(drapeau == drapeau3);

    // Verifie qu'il n'y a qu'une phase liquide :
    check_somme_drapeau(drapeau);
    check_somme_drapeau(drapeau1);
    check_somme_drapeau(drapeau1b);
    check_somme_drapeau(drapeau2);
    check_somme_drapeau(drapeau3);
  }
#endif

#define DEBUG_INDIC 0
#if DEBUG_INDIC
  // Pour debugger: on met dans l'indicatrice le numero de la compo connexe
  // globale:
  {
    for (int k = 0; k < nk; k++)
      {
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                const int num_elem = s.convert_ijk_cell_to_packed(i, j, k);
                int compo = num_compo[num_elem];
                // Element dans une bulle:
                indic(i, j, k) = compo;
              }
          }
      }
  }
#else
  {
    // Mettre a zero l'indicatrice pour les elements qui sont dans des
    // composantes marquees:
    for (int k = 0; k < nk; k++)
      {
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                const int num_elem = s.convert_ijk_cell_to_packed(i, j, k);
                int compo = num_compo_[num_elem];
                // color(i,j,k) = compo;
                if (compo >= 0 && drapeau[compo] == 0)
                  {
                    // Element dans une bulle:
                    indic(i, j, k) = 0.;
                  }
              }
          }
      }
  }
#endif
  statistiques().end_count(calculer_indicatrice_counter_);
}

void IJK_Interfaces::calculer_indicatrice_optim(IJK_Field_double& indic)
{
  static Stat_Counter_Id calculer_indicatrice_counter_ =
    statistiques().new_counter(2, "Calcul rho mu indicatrice: calcul de l'indicatrice");
  statistiques().begin_count(calculer_indicatrice_counter_);

  const Intersections_Elem_Facettes& intersec = maillage_ft_ijk_.intersections_elem_facettes();
  const Domaine_IJK& s = indic.get_domaine();

  const int ni = indic.ni();
  const int nj = indic.nj();
  const int nk = indic.nk();

  // calcul des fractions de volume de phase 1 dans les mailles traversees par
  // des interfaces:
  {
    const ArrOfInt& index_elem = intersec.index_elem();
    // Boucle sur les elements euleriens
    for (int k = 0; k < nk; k++)
      {
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                const int num_elem = s.convert_ijk_cell_to_packed(i, j, k);
                int index = index_elem[num_elem];

                if (0. < indic(i, j, k) && indic(i, j, k) < 1.) // element traverse par l'interface
                  indic(i, j, k) = 2.;

                double somme_contrib = 0.;
                // Boucle sur les facettes qui traversent cet element
                while (index >= 0)
                  {
                    const Intersections_Elem_Facettes_Data& data = intersec.data_intersection(index);
                    somme_contrib += data.contrib_volume_phase1_;
                    index = data.index_facette_suivante_;
                  };

                int nb_increment_somme_contrib = 0;
                while (somme_contrib > 1.)
                  {
                    somme_contrib -= 1.;
                    nb_increment_somme_contrib++;
                  }
                while (somme_contrib < 0.)
                  {
                    somme_contrib += 1.;
                    nb_increment_somme_contrib++;
                  }

                if (nb_increment_somme_contrib > 1)
                  Process::exit("Error in IJK_Interfaces::calculer_indicatrice_optim !");

                if (somme_contrib > 0.)
                  {
                    // if(somme_contrib == 1.)
                    //   somme_contrib -= 1e-3;
                    indic(i, j, k) = somme_contrib;
                  }
              }
          }
      }
  }

  update_indicatrice(indic);
  statistiques().end_count(calculer_indicatrice_counter_);
}

void IJK_Interfaces::calculer_indicatrices(IJK_Field_vector3_double& indic)
{
  static Stat_Counter_Id calculer_indicatrice_counter_ =
    statistiques().new_counter(2, "Calcul rho mu indicatrice: calcul des indicatrices");
  statistiques().begin_count(calculer_indicatrice_counter_);

  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, refdomaine_dis_.valeur());
  const Domaine& domaine = domaine_vf.domaine();
  const IntTab& elem_faces = domaine_vf.elem_faces();
  const IntTab& faces_voisins = domaine_vf.face_voisins();
  const Intersections_Elem_Facettes& intersec = maillage_ft_ijk_.intersections_elem_facettes();
  const Domaine_IJK& s = indic.get_domaine();

  for (int igroup = 0; igroup < nb_groups_; igroup++)
    {
      indic[igroup].data() = 1.; // Tout liquide pour commencer
    }
  if (parser_ == 0)
    {
      Cerr << "This choice of options has to be validated. "
           << "It seems you are using a multi-groups of bubbles "
           << "and the forcing method with the color_function "
           << "which has not been tested yet!"
           << finl;
      Process::exit();
    }
  // ISO C++ forbids variable length array : IntVect num_compo[nb_groups_];
  IntVect num_compo[max_authorized_nb_of_groups_];

  // Cree un tableau parallele structure comme un tableau aux elements
  // du maillage vdf, initialise a zero.
  for (int i = 0; i < nb_groups_; i++)
    domaine.creer_tableau_elements(num_compo[i]);

  const int ni = indic[0].ni();
  const int nj = indic[0].nj();
  const int nk = indic[0].nk();

  assert(nb_groups_ <= max_authorized_nb_of_groups_);
  const ArrOfInt& compo_connexe = maillage_ft_ijk_.compo_connexe_facettes();

  // La methode de maillage_ft_disc utilise l'equation de transport pour
  // recuperer la fonction distance et un bricolage pour trouver l'indicatrice
  // au voisinage de l'interface. On va faire autrement. Debut identique a la
  // methode VDF: calcul des fractions de volume de phase 1 dans les mailles
  // traversees par des interfaces:
  {
    const ArrOfInt& index_elem = intersec.index_elem();
    // Boucle sur les elements euleriens

    for (int k = 0; k < nk; k++)
      {
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                // Anciennement la methode etait portee par le mesh :
                //    const int num_elem =
                // maillage_ft_ijk_.convert_ijk_cell_to_packed(i, j, k);
                // A present, elle est dans le splitting :
                const int num_elem = s.convert_ijk_cell_to_packed(i, j, k);
                int index = index_elem[num_elem];
                double somme_contrib[max_authorized_nb_of_groups_] = {0.};
                // Boucle sur les facettes qui traversent cet element
                while (index >= 0)
                  {
                    const Intersections_Elem_Facettes_Data& data = intersec.data_intersection(index);
                    const int facette = data.numero_facette_;
                    int icompo = compo_connexe[facette];
                    int m = icompo;
                    if (icompo < 0)
                      {
                        // Portion d'interface ghost.
                        // On recherche le vrai numero
                        m = decoder_numero_bulle(-icompo);
                      }
                    assert(0 <= m && m < nb_bulles_reelles_);
                    const int igroup = compo_to_group_[m];
                    somme_contrib[igroup] += data.contrib_volume_phase1_;

                    index = data.index_facette_suivante_;
                  };
                for (int igroup = 0; igroup < nb_groups_; igroup++)
                  {
                    while (somme_contrib[igroup] > 1.)
                      somme_contrib[igroup] -= 1.;
                    while (somme_contrib[igroup] < 0.)
                      somme_contrib[igroup] += 1.;
                    if (somme_contrib[igroup] > 0.)
                      {
                        // if(somme_contrib[igroup] == 1.)
                        // somme_contrib[igroup] -= 1e-3;
                        indic[igroup](i, j, k) = somme_contrib[igroup];
                        num_compo[igroup][num_elem] = -1; // marque l'element pour apres
                      }
                  }
              }
          }
      }
  }

  ArrOfInt drapeau[max_authorized_nb_of_groups_];
  static Stat_Counter_Id search_connex_components_counter_ =
    statistiques().new_counter(2, "Calcul de l'indicatrice : recherche compo connexes");
  for (int i = 0; i < nb_groups_; i++)
    {
      statistiques().begin_count(search_connex_components_counter_);
      num_compo[i].echange_espace_virtuel();
      // Recherche des composantes connexes sur le maillage eulerien
      int nb_compo_locales = search_connex_components_local(elem_faces, faces_voisins, num_compo[i]);
      nb_compo_in_num_compo_ = compute_global_connex_components(num_compo[i], nb_compo_locales);
      num_compo[i].echange_espace_virtuel();
      statistiques().end_count(search_connex_components_counter_);

      // Calcul des drapeaux selon la methode la plus robuste :
      drapeau[i].resize_array(nb_compo_in_num_compo_);
      compute_drapeaux_vapeur_v4(num_compo[i], drapeau[i]);
    }

  {
    // Mettre a zero l'indicatrice pour les elements qui sont dans des
    // composantes marquees:
    for (int k = 0; k < nk; k++)
      {
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                const int num_elem = s.convert_ijk_cell_to_packed(i, j, k);
                for (int igroup = 0; igroup < nb_groups_; igroup++)
                  {
                    int compo = num_compo[igroup][num_elem];
                    if (compo >= 0 && drapeau[igroup][compo] == 0)
                      {
                        // Element dans une bulle:
                        indic[igroup](i, j, k) = 0.;
                      }
                  }
              }
          }
      }
  }
  statistiques().end_count(calculer_indicatrice_counter_);
}

void IJK_Interfaces::calculer_indicatrices_optim(IJK_Field_vector3_double& indic)
{
  static Stat_Counter_Id calculer_indicatrice_counter_ =
    statistiques().new_counter(2, "Calcul rho mu indicatrice: calcul des indicatrices");
  statistiques().begin_count(calculer_indicatrice_counter_);

  const Intersections_Elem_Facettes& intersec = maillage_ft_ijk_.intersections_elem_facettes();
  const Domaine_IJK& s = indic.get_domaine();

  const int ni = indic[0].ni();
  const int nj = indic[0].nj();
  const int nk = indic[0].nk();

  assert(nb_groups_ <= max_authorized_nb_of_groups_);
  const ArrOfInt& compo_connexe = maillage_ft_ijk_.compo_connexe_facettes();

  // calcul des fractions de volume de phase 1 dans les mailles traversees par
  // des interfaces:
  {
    const ArrOfInt& index_elem = intersec.index_elem();
    // Boucle sur les elements euleriens

    for (int k = 0; k < nk; k++)
      {
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                const int num_elem = s.convert_ijk_cell_to_packed(i, j, k);
                int index = index_elem[num_elem];
                double somme_contrib[max_authorized_nb_of_groups_] = {0.};

                for (int igroup = 0; igroup < nb_groups_; igroup++)
                  {
                    if (0. < indic[igroup](i, j, k) && indic[igroup](i, j, k) < 1.) // element traverse par l'indicatrice
                      indic[igroup](i, j, k) = 2.;
                  }

                // Boucle sur les facettes qui traversent cet element
                while (index >= 0)
                  {
                    const Intersections_Elem_Facettes_Data& data = intersec.data_intersection(index);
                    const int facette = data.numero_facette_;
                    int icompo = compo_connexe[facette];
                    int m = icompo;
                    if (icompo < 0)
                      {
                        // Portion d'interface ghost.
                        // On recherche le vrai numero
                        m = decoder_numero_bulle(-icompo);
                      }
                    assert(0 <= m && m < nb_bulles_reelles_);
                    const int igroup = compo_to_group_[m];
                    somme_contrib[igroup] += data.contrib_volume_phase1_;
                    index = data.index_facette_suivante_;
                  };
                for (int igroup = 0; igroup < nb_groups_; igroup++)
                  {
                    while (somme_contrib[igroup] > 1.)
                      somme_contrib[igroup] -= 1.;
                    while (somme_contrib[igroup] < 0.)
                      somme_contrib[igroup] += 1.;
                    if (somme_contrib[igroup] > 0.)
                      {
                        indic[igroup](i, j, k) = somme_contrib[igroup];
                      }
                  }
              }
          }
      }
  }

  for (int igroup = 0; igroup < nb_groups_; igroup++)
    {
      update_indicatrice(indic[igroup]);
    }

  statistiques().end_count(calculer_indicatrice_counter_);
}

int IJK_Interfaces::update_indicatrice(IJK_Field_double& indic)
{
  const int ni = indic.ni();
  const int nj = indic.nj();
  const int nk = indic.nk();

  const Domaine_IJK& s = indic.get_domaine();
  const int imin = s.get_offset_local(DIRECTION_I);
  const int jmin = s.get_offset_local(DIRECTION_J);
  const int kmin = s.get_offset_local(DIRECTION_K);

  const Domaine_IJK::Localisation loc = indic.get_localisation();
  const int nitot = s.get_nb_items_global(loc, DIRECTION_I);
  const int njtot = s.get_nb_items_global(loc, DIRECTION_J);
  const int nktot = s.get_nb_items_global(loc, DIRECTION_K);

  indic.echange_espace_virtuel(indic.ghost());

  //(l'algo peut etre optimise si on stocke les elements a phase indeterminee
  // plutot que de parcourir tous les elements du processeur a chaque passe)
  int continuer;
  do
    {
      continuer = 0;
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  if (indic(i, j, k) == 2.) // element a phase indeterminee: on doit le mettre a jour
                    {
                      ArrOfInt voisins_num(6);
                      voisins_num[0] = i == 0 ? BORD : s.convert_ijk_cell_to_packed(i - 1, j, k);
                      voisins_num[1] = j == 0 ? BORD : s.convert_ijk_cell_to_packed(i, j - 1, k);
                      voisins_num[2] = k == 0 ? BORD : s.convert_ijk_cell_to_packed(i, j, k - 1);
                      voisins_num[3] = i == ni - 1 ? BORD : s.convert_ijk_cell_to_packed(i + 1, j, k);
                      voisins_num[4] = j == nj - 1 ? BORD : s.convert_ijk_cell_to_packed(i, j + 1, k);
                      voisins_num[5] = k == nk - 1 ? BORD : s.convert_ijk_cell_to_packed(i, j, k + 1);

                      ArrOfDouble voisins_valeur(6);
                      voisins_valeur[0] = imin + i == 0 ? BORD : indic(i - 1, j, k);
                      voisins_valeur[1] = jmin + j == 0 ? BORD : indic(i, j - 1, k);
                      voisins_valeur[2] = kmin + k == 0 ? BORD : indic(i, j, k - 1);
                      voisins_valeur[3] = imin + i == nitot ? BORD : indic(i + 1, j, k);
                      voisins_valeur[4] = jmin + j == njtot ? BORD : indic(i, j + 1, k);
                      voisins_valeur[5] = kmin + k == nktot ? BORD : indic(i, j, k + 1);

                      int found = 0; // indique si on a trouve la phase de l'element ijk

                      // on regarde si parmi nos voisins, il y en a un qui est liquide ou
                      // gazeux on regarde meme les voisins virtuels si oui, on adopte la
                      // meme phase que lui
                      for (int v = 0; v < 6; v++)
                        {
                          double valeur_indic_voisin = voisins_valeur[v];
                          if (valeur_indic_voisin == 0. || valeur_indic_voisin == 1.)
                            {
                              indic(i, j, k) = valeur_indic_voisin;
                              found = 1;
                              break;
                            }
                        }

                      // sinon, on regarde si parmi nos voisins, il y a des cellules
                      // coupees par l'interface on ne regarde que les voisins reels,
                      // parce qu'on aura ensuite besoin de regarder les facettes qui les
                      // traversent (pas possible avec les elements virtuels ? a voir...)
                      if (!found)
                        {
                          for (int v = 0; v < 6; v++)
                            {
                              double valeur_indic_voisin = voisins_valeur[v];
                              const int direction = v % 3;
                              const int face_plus = (v > 2) ? 1 : -1; // +1 si c'est notre voisin de droite

                              if (0. < valeur_indic_voisin && valeur_indic_voisin < 1.)
                                {
                                  int num_elem = voisins_num[v];
                                  if (num_elem != BORD)
                                    indic(i, j, k) =
                                      compute_cell_phase_with_interface_normal(num_elem,
                                                                               direction,
                                                                               -face_plus /* si le voisin courant est mon voisin de droite, alors notre face commune se trouve a sa gauche */);
                                  if (indic(i, j, k) != -1)
                                    {
                                      found = 1;
                                      break;
                                    }
                                }
                            }
                        }

                      if (!found)
                        {
                          // element sur le bord entoure uniquement de cellules a phase
                          // indeterminee ses voisins a phase indeterminee ont forcement des
                          // cellules voisines traversees par l'interface ils pourront donc
                          // trouver la valeur de leur indicatrice lors de la passe courante
                          // la phase de l'element actuel pourra donc aussi etre trouvee
                          // lors de la prochaine passe
                          continuer = 1;
                        }
                    }
                }
            }
        }
    }
  while (continuer);

  return 0;
}

// Calcul de l'indicatrice surfacique et du barycentre de la phase sur les faces euleriennes
void IJK_Interfaces::calculer_indicatrice_surfacique_barycentre_face(IJK_Field_vector3_double& indic_surfacique_face, FixedVector<FixedVector<IJK_Field_double, 2>, 3>& baric_face, IJK_Field_double& indic, IJK_Field_vector3_double& norme)
{
  static Stat_Counter_Id calculer_indicatrice_surfacique_face_counter_ =
    statistiques().new_counter(2, "Calcul rho mu indicatrice: calcul de l'indicatrice surface face");
  statistiques().begin_count(calculer_indicatrice_surfacique_face_counter_);

  const Intersections_Elem_Facettes& intersec = maillage_ft_ijk_.intersections_elem_facettes();
  const Domaine_IJK& s = indic_surfacique_face.get_domaine();

  const int ni = indic_surfacique_face[0].ni();
  const int nj = indic_surfacique_face[0].nj();
  const int nk = indic_surfacique_face[0].nk();

  // Initialisation
  {
    for (int k = 0; k < nk; k++)
      {
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                if (est_pure(indic(i, j, k)))
                  {
                    // Dans les cellules pures, on utilise l'indicatrice pour determiner la phase
                    indic_surfacique_face[0](i, j, k) = indic(i, j, k);
                    indic_surfacique_face[1](i, j, k) = indic(i, j, k);
                    indic_surfacique_face[2](i, j, k) = indic(i, j, k);
                  }
                else
                  {
                    // Dans les cellules diphasiques, on determine la phase a partir de la normale a l'interface
                    // Ce calcul est important si la maille est diphasique mais que l'interface ne coupe pas la face
                    indic_surfacique_face[0](i, j, k) = norme[0](i, j, k) == 0. ? 0. : (norme[0](i, j, k) > 0. ? 0. : 1.);
                    indic_surfacique_face[1](i, j, k) = norme[1](i, j, k) == 0. ? 0. : (norme[1](i, j, k) > 0. ? 0. : 1.);
                    indic_surfacique_face[2](i, j, k) = norme[2](i, j, k) == 0. ? 0. : (norme[2](i, j, k) > 0. ? 0. : 1.);

                    // Ce n'est pas toujours vrai, on corrige le cas ou une des cellules adjacentes est pure
                    // Note : Cette methode pourrait ne pas toujours fonctionner ?
                    if (est_pure(indic(i-1,j,k)))
                      {
                        indic_surfacique_face[0](i, j, k) = indic(i-1, j, k);
                      }
                    if (est_pure(indic(i,j-1,k)))
                      {
                        indic_surfacique_face[1](i, j, k) = indic(i, j-1, k);
                      }
                    if (est_pure(indic(i,j,k-1)))
                      {
                        indic_surfacique_face[2](i, j, k) = indic(i, j, k-1);
                      }
                  }

                // Le barycentre est par defaut au milieu de la maille
                baric_face[0][0](i, j, k) = 0.5;
                baric_face[0][1](i, j, k) = 0.5;
                baric_face[1][0](i, j, k) = 0.5;
                baric_face[1][1](i, j, k) = 0.5;
                baric_face[2][0](i, j, k) = 0.5;
                baric_face[2][1](i, j, k) = 0.5;
              }
          }
      }
  }

  // Correction pour les faces coupees par l'interface
  // Note : methode similaire a calculer_indicatrice
  {
    const ArrOfInt& index_elem = intersec.index_elem();
    //    const int nb_elem = index_elem.size_array();
    // Boucle sur les elements euleriens
    for (int k = 0; k < nk; k++)
      {
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                // qu'une cellule consideree pure soit traversee par l'interface.
                // Pour coherence, on considere les faces non coupees dans ce cas.
                if (!est_pure(indic(i, j, k)))
                  {
                    // Anciennement la methode etait portee par le mesh :
                    //    const int num_elem =
                    // maillage_ft_ijk_.convert_ijk_cell_to_packed(i, j, k);
                    // A present, elle est dans le splitting :
                    const int num_elem = s.convert_ijk_cell_to_packed(i, j, k);
                    int index = index_elem[num_elem];
                    double somme_contrib[3] = {0., 0., 0.};
                    double somme_contrib_baryc[3][2] = {{0., 0.}, {0., 0.}, {0., 0.}};
                    // Boucle sur les facettes qui traversent cet element
                    while (index >= 0)
                      {
                        const Intersections_Elem_Facettes_Data& data = intersec.data_intersection(index);
                        somme_contrib[0] += data.contrib_aire_faces_phase1_[0];
                        somme_contrib[1] += data.contrib_aire_faces_phase1_[1];
                        somme_contrib[2] += data.contrib_aire_faces_phase1_[2];

                        somme_contrib_baryc[0][0] += data.contrib_barycentre_faces_phase1_[0][0]*abs(data.contrib_aire_faces_phase1_[0]);
                        somme_contrib_baryc[0][1] += data.contrib_barycentre_faces_phase1_[0][1]*abs(data.contrib_aire_faces_phase1_[0]);
                        somme_contrib_baryc[1][0] += data.contrib_barycentre_faces_phase1_[1][0]*abs(data.contrib_aire_faces_phase1_[1]);
                        somme_contrib_baryc[1][1] += data.contrib_barycentre_faces_phase1_[1][1]*abs(data.contrib_aire_faces_phase1_[1]);
                        somme_contrib_baryc[2][0] += data.contrib_barycentre_faces_phase1_[2][0]*abs(data.contrib_aire_faces_phase1_[2]);
                        somme_contrib_baryc[2][1] += data.contrib_barycentre_faces_phase1_[2][1]*abs(data.contrib_aire_faces_phase1_[2]);

                        index = data.index_facette_suivante_;
                      };

                    for (int dir=0; dir<3; dir++)
                      {
                        const int dir_i = (dir == 0);
                        const int dir_j = (dir == 1);
                        const int dir_k = (dir == 2);

                        if (est_pure(indic(i-dir_i,j-dir_j,k-dir_k)))
                          {
                            // Une cellule adjacente pure suggere ou bien qu'il n'y a pas d'intersections
                            // ou bien qu'il y a des intersections mais que l'indicatrice resultante est
                            // trop faible et en dessous du seuil dans calculer_indicatrice.
                            // Pour coherence, on neglige la coupure de la surface egalement.
                          }
                        else
                          {
                            // Dans chaque direction, on ne touche qu'aux faces coupees
                            if (somme_contrib[dir]*somme_contrib[dir] > 0)
                              {
                                while (somme_contrib[dir] > 1.)
                                  {
                                    somme_contrib[dir] -= 1.;
                                    somme_contrib_baryc[dir][0] -= 1./2.;
                                    somme_contrib_baryc[dir][1] -= 1./2.;
                                  }
                                while (somme_contrib[dir] < 0.)
                                  {
                                    somme_contrib[dir] += 1.;
                                    somme_contrib_baryc[dir][0] += 1./2.;
                                    somme_contrib_baryc[dir][1] += 1./2.;
                                  }

                                indic_surfacique_face[dir](i, j, k) = somme_contrib[dir];

                                if (est_pure(somme_contrib[dir]))
                                  {
                                    baric_face[dir][0](i, j, k) = 1./2.;
                                    baric_face[dir][1](i, j, k) = 1./2.;
                                  }
                                else
                                  {
                                    baric_face[dir][0](i, j, k) = somme_contrib_baryc[dir][0]/abs(somme_contrib[dir]);
                                    baric_face[dir][1](i, j, k) = somme_contrib_baryc[dir][1]/abs(somme_contrib[dir]);

                                    // Un barycentre legerement negatif ou superieur a un est indicatif
                                    // d'une imprecision dans le calcul du barycentre.
                                    if ((baric_face[dir][0](i, j, k) <= 0) && (baric_face[dir][0](i, j, k) >= -1e-12 || somme_contrib[dir] < 1e-8))
                                      {
                                        baric_face[dir][0](i, j, k) = 1e-12;
                                      }
                                    else if ((baric_face[dir][0](i, j, k) >= 1) && (baric_face[dir][0](i, j, k) <= 1+1e-12 || somme_contrib[dir] < 1e-8))
                                      {
                                        baric_face[dir][0](i, j, k) = 1 - 1e-12;
                                      }
                                    if ((baric_face[dir][1](i, j, k) <= 0) && (baric_face[dir][1](i, j, k) >= -1e-12 || somme_contrib[dir] < 1e-8))
                                      {
                                        baric_face[dir][1](i, j, k) = 1e-12;
                                      }
                                    else if ((baric_face[dir][1](i, j, k) >= 1) && (baric_face[dir][1](i, j, k) <= 1+1e-12 || somme_contrib[dir] < 1e-8))
                                      {
                                        baric_face[dir][1](i, j, k) = 1 - 1e-12;
                                      }

                                    assert((baric_face[dir][0](i, j, k) >= 0) && (baric_face[dir][0](i, j, k) <= 1));
                                    assert((baric_face[dir][1](i, j, k) >= 0) && (baric_face[dir][1](i, j, k) <= 1));
                                    assert(((indic_surfacique_face[dir](i, j, k) < 0.) && (baric_face[dir][0](i, j, k) < 0.)) || ((indic_surfacique_face[dir](i, j, k) > 0.) && (baric_face[dir][0](i, j, k) > 0.)));
                                    assert(((indic_surfacique_face[dir](i, j, k) < 0.) && (baric_face[dir][1](i, j, k) < 0.)) || ((indic_surfacique_face[dir](i, j, k) > 0.) && (baric_face[dir][1](i, j, k) > 0.)));
                                  }
                              }
                          }
                      }
                  }
              }
          }
      }
  }
  statistiques().end_count(calculer_indicatrice_surfacique_face_counter_);
}

// Calcul de l'indicatrice surfacique sur les faces euleriennes (sans le barycentre)
void IJK_Interfaces::calculer_indicatrice_surfacique_face(IJK_Field_vector3_double& indic_surfacique_face, IJK_Field_double& indic, IJK_Field_vector3_double& norme)
{
  static Stat_Counter_Id calculer_indicatrice_surfacique_face_counter_ =
    statistiques().new_counter(2, "Calcul rho mu indicatrice: calcul de l'indicatrice surface face");
  statistiques().begin_count(calculer_indicatrice_surfacique_face_counter_);

  const Intersections_Elem_Facettes& intersec = maillage_ft_ijk_.intersections_elem_facettes();
  const Domaine_IJK& s = indic_surfacique_face.get_domaine();

  const int ni = indic_surfacique_face[0].ni();
  const int nj = indic_surfacique_face[0].nj();
  const int nk = indic_surfacique_face[0].nk();

  // Initialisation
  {
    for (int k = 0; k < nk; k++)
      {
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                if (est_pure(indic(i, j, k)))
                  {
                    // Dans les cellules pures, on utilise l'indicatrice pour determiner la phase
                    indic_surfacique_face[0](i, j, k) = indic(i, j, k);
                    indic_surfacique_face[1](i, j, k) = indic(i, j, k);
                    indic_surfacique_face[2](i, j, k) = indic(i, j, k);
                  }
                else
                  {
                    // Dans les cellules diphasiques, on determine la phase a partir de la normale a l'interface
                    // Ce calcul est important si la maille est diphasique mais que l'interface ne coupe pas la face
                    indic_surfacique_face[0](i, j, k) = norme[0](i, j, k) == 0. ? 0. : (norme[0](i, j, k) > 0. ? 0. : 1.);
                    indic_surfacique_face[1](i, j, k) = norme[1](i, j, k) == 0. ? 0. : (norme[1](i, j, k) > 0. ? 0. : 1.);
                    indic_surfacique_face[2](i, j, k) = norme[2](i, j, k) == 0. ? 0. : (norme[2](i, j, k) > 0. ? 0. : 1.);

                    // Ce n'est pas toujours vrai, on corrige le cas ou une des cellules adjacentes est pure
                    // Note : Cette methode pourrait ne pas toujours fonctionner ?
                    if (est_pure(indic(i-1,j,k)))
                      {
                        indic_surfacique_face[0](i, j, k) = indic(i-1, j, k);
                      }
                    if (est_pure(indic(i,j-1,k)))
                      {
                        indic_surfacique_face[1](i, j, k) = indic(i, j-1, k);
                      }
                    if (est_pure(indic(i,j,k-1)))
                      {
                        indic_surfacique_face[2](i, j, k) = indic(i, j, k-1);
                      }
                  }
              }
          }
      }
  }

  // Correction pour les faces coupees par l'interface
  // Note : methode similaire a calculer_indicatrice
  {
    const ArrOfInt& index_elem = intersec.index_elem();
    //    const int nb_elem = index_elem.size_array();
    // Boucle sur les elements euleriens
    for (int k = 0; k < nk; k++)
      {
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                // qu'une cellule consideree pure soit traversee par l'interface.
                // Pour coherence, on considere les faces non coupees dans ce cas.
                if (!est_pure(indic(i, j, k)))
                  {
                    // Anciennement la methode etait portee par le mesh :
                    //    const int num_elem =
                    // maillage_ft_ijk_.convert_ijk_cell_to_packed(i, j, k);
                    // A present, elle est dans le splitting :
                    const int num_elem = s.convert_ijk_cell_to_packed(i, j, k);
                    int index = index_elem[num_elem];
                    double somme_contrib[3] = {0., 0., 0.};
                    // Boucle sur les facettes qui traversent cet element
                    while (index >= 0)
                      {
                        const Intersections_Elem_Facettes_Data& data = intersec.data_intersection(index);
                        somme_contrib[0] += data.contrib_aire_faces_phase1_[0];
                        somme_contrib[1] += data.contrib_aire_faces_phase1_[1];
                        somme_contrib[2] += data.contrib_aire_faces_phase1_[2];

                        index = data.index_facette_suivante_;
                      };

                    for (int dir=0; dir<3; dir++)
                      {
                        const int dir_i = (dir == 0);
                        const int dir_j = (dir == 1);
                        const int dir_k = (dir == 2);

                        if (est_pure(indic(i-dir_i,j-dir_j,k-dir_k)))
                          {
                            // Une cellule adjacente pure suggere ou bien qu'il n'y a pas d'intersections
                            // ou bien qu'il y a des intersections mais que l'indicatrice resultante est
                            // trop faible et en dessous du seuil dans calculer_indicatrice.
                            // Pour coherence, on neglige la coupure de la surface egalement.
                          }
                        else
                          {
                            // Dans chaque direction, on ne touche qu'aux faces coupees
                            if (somme_contrib[dir]*somme_contrib[dir] > 0)
                              {
                                while (somme_contrib[dir] > 1.)
                                  {
                                    somme_contrib[dir] -= 1.;
                                  }
                                while (somme_contrib[dir] < 0.)
                                  {
                                    somme_contrib[dir] += 1.;
                                  }

                                indic_surfacique_face[dir](i, j, k) = somme_contrib[dir];
                              }
                          }
                      }
                  }
              }
          }
      }
  }
  statistiques().end_count(calculer_indicatrice_surfacique_face_counter_);
}


void IJK_Interfaces::calculer_surface_interface(IJK_Field_double& surf_interface, IJK_Field_double& indic)
{
  static Stat_Counter_Id calculer_surface_interface_counter_ =
    statistiques().new_counter(2, "Calcul rho mu indicatrice: calcul de la surface interfaciale");
  statistiques().begin_count(calculer_surface_interface_counter_);

  const Intersections_Elem_Facettes& intersec = maillage_ft_ijk_.intersections_elem_facettes();
  const Domaine_IJK& s = surf_interface.get_domaine();

  const ArrOfDouble& surface_facettes = maillage_ft_ijk_.get_update_surface_facettes();

  const int ni = surf_interface.ni();
  const int nj = surf_interface.nj();
  const int nk = surf_interface.nk();

  // Initialisation
  {
    for (int k = 0; k < nk; k++)
      {
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                surf_interface(i, j, k) = 0.;
              }
          }
      }
  }

  // Correction pour les faces coupees par l'interface
  // Note : methode similaire a calculer_indicatrice
  {
    const ArrOfInt& index_elem = intersec.index_elem();
    //    const int nb_elem = index_elem.size_array();
    // Boucle sur les elements euleriens
    for (int k = 0; k < nk; k++)
      {
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                // qu'une cellule consideree pure soit traversee par l'interface.
                if (!est_pure(indic(i, j, k)))
                  {
                    // Anciennement la methode etait portee par le mesh :
                    //    const int num_elem =
                    // maillage_ft_ijk_.convert_ijk_cell_to_packed(i, j, k);
                    // A present, elle est dans le splitting :
                    const int num_elem = s.convert_ijk_cell_to_packed(i, j, k);
                    int index = index_elem[num_elem];
                    double somme_contrib_surf = 0.;
                    // Boucle sur les facettes qui traversent cet element
                    while (index >= 0)
                      {
                        const Intersections_Elem_Facettes_Data& data = intersec.data_intersection(index);
                        const int fa7 = data.numero_facette_;
                        somme_contrib_surf += data.fraction_surface_intersection_ * surface_facettes[fa7];

                        index = data.index_facette_suivante_;
                      };

                    surf_interface(i, j, k) = somme_contrib_surf;
                    assert(somme_contrib_surf >= 0.);
                  }
              }
          }
      }
  }
  statistiques().end_count(calculer_surface_interface_counter_);
}

void IJK_Interfaces::calculer_barycentre(IJK_Field_vector3_double& baric, IJK_Field_double& indic)
{
  static Stat_Counter_Id calculer_barycentre_counter_ =
    statistiques().new_counter(2, "Calcul rho mu indicatrice: calcul du barycentre");
  statistiques().begin_count(calculer_barycentre_counter_);

  const Intersections_Elem_Facettes& intersec = maillage_ft_ijk_.intersections_elem_facettes();
  const Domaine_IJK& s = baric.get_domaine();

  const int ni = baric[0].ni();
  const int nj = baric[0].nj();
  const int nk = baric[0].nk();

  // Initialisation
  {
    for (int k = 0; k < nk; k++)
      {
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                // Le barycentre est par defaut au milieu de la maille
                baric[0](i, j, k) = 0.5;
                baric[1](i, j, k) = 0.5;
                baric[2](i, j, k) = 0.5;
              }
          }
      }
  }

  // Correction pour les faces coupees par l'interface
  // Note : methode similaire a calculer_indicatrice
  {
    const ArrOfInt& index_elem = intersec.index_elem();
    //    const int nb_elem = index_elem.size_array();
    // Boucle sur les elements euleriens
    for (int k = 0; k < nk; k++)
      {
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                // Puisqu'il y a une tolerance sur le calcul de l'indicatrice, il se peut
                // qu'une cellule consideree pure soit traversee par l'interface.
                if (!est_pure(indic(i, j, k)))
                  {
                    // Anciennement la methode etait portee par le mesh :
                    //    const int num_elem =
                    // maillage_ft_ijk_.convert_ijk_cell_to_packed(i, j, k);
                    // A present, elle est dans le splitting :
                    const int num_elem = s.convert_ijk_cell_to_packed(i, j, k);
                    int index = index_elem[num_elem];
                    double somme_contrib = 0.;
                    double somme_contrib_baryc[3] = {0., 0., 0.};
                    // Boucle sur les facettes qui traversent cet element
                    while (index >= 0)
                      {
                        const Intersections_Elem_Facettes_Data& data = intersec.data_intersection(index);
                        somme_contrib += data.contrib_volume_phase1_;

                        somme_contrib_baryc[0] += data.contrib_barycentre_phase1_[0]*abs(data.contrib_volume_phase1_);
                        somme_contrib_baryc[1] += data.contrib_barycentre_phase1_[1]*abs(data.contrib_volume_phase1_);
                        somme_contrib_baryc[2] += data.contrib_barycentre_phase1_[2]*abs(data.contrib_volume_phase1_);

                        index = data.index_facette_suivante_;
                      };

                    // int nb_increment_somme_contrib = 0;
                    while (somme_contrib > 1.)
                      {
                        somme_contrib -= 1.;
                        somme_contrib_baryc[0] -= 1./2.;
                        somme_contrib_baryc[1] -= 1./2.;
                        somme_contrib_baryc[2] -= 1./2.;
                        // nb_increment_somme_contrib++;
                      }
                    while (somme_contrib < 0.)
                      {
                        somme_contrib += 1.;
                        somme_contrib_baryc[0] += 1./2.;
                        somme_contrib_baryc[1] += 1./2.;
                        somme_contrib_baryc[2] += 1./2.;
                        // nb_increment_somme_contrib++;
                      }

                    // Note : On recalcule l'indicatrice, c'est un peu dommage
                    assert(indic(i, j, k) == somme_contrib);

                    baric[0](i, j, k) = somme_contrib_baryc[0]/abs(somme_contrib);
                    baric[1](i, j, k) = somme_contrib_baryc[1]/abs(somme_contrib);
                    baric[2](i, j, k) = somme_contrib_baryc[2]/abs(somme_contrib);

                    // Un barycentre legerement negatif ou superieur a un est indicatif
                    // d'une imprecision dans le calcul du barycentre.
                    if ((baric[0](i, j, k) <= 0) && (baric[0](i, j, k) >= -1e-12 || somme_contrib < 1e-8))
                      {
                        baric[0](i, j, k) = 1e-12;
                      }
                    else if ((baric[0](i, j, k) >= 1) && (baric[0](i, j, k) <= 1+1e-12 || somme_contrib < 1e-8))
                      {
                        baric[0](i, j, k) = 1 - 1e-12;
                      }
                    if ((baric[1](i, j, k) <= 0) && (baric[1](i, j, k) >= -1e-12 || somme_contrib < 1e-8))
                      {
                        baric[1](i, j, k) = 1e-12;
                      }
                    else if ((baric[1](i, j, k) >= 1) && (baric[1](i, j, k) <= 1+1e-12 || somme_contrib < 1e-8))
                      {
                        baric[1](i, j, k) = 1 - 1e-12;
                      }
                    if ((baric[2](i, j, k) <= 0) && (baric[2](i, j, k) >= -1e-12 || somme_contrib < 1e-8))
                      {
                        baric[2](i, j, k) = 1e-12;
                      }
                    else if ((baric[2](i, j, k) >= 1) && (baric[2](i, j, k) <= 1+1e-12 || somme_contrib < 1e-8))
                      {
                        baric[2](i, j, k) = 1 - 1e-12;
                      }

                    assert((baric[0](i, j, k) >= 0) && (baric[0](i, j, k) <= 1));
                    assert((baric[1](i, j, k) >= 0) && (baric[1](i, j, k) <= 1));
                    assert((baric[2](i, j, k) >= 0) && (baric[2](i, j, k) <= 1));
                    assert(((indic(i, j, k) < 0.) && (baric[0](i, j, k) < 0.)) || ((indic(i, j, k) > 0.) && (baric[0](i, j, k) > 0.)));
                    assert(((indic(i, j, k) < 0.) && (baric[1](i, j, k) < 0.)) || ((indic(i, j, k) > 0.) && (baric[1](i, j, k) > 0.)));
                    assert(((indic(i, j, k) < 0.) && (baric[2](i, j, k) < 0.)) || ((indic(i, j, k) > 0.) && (baric[2](i, j, k) > 0.)));
                  }
              }
          }
      }
  }
  statistiques().end_count(calculer_barycentre_counter_);
}

////////////////////////////////////////////////////////////////////////////////

void IJK_Interfaces::convert_to_IntVect(const ArrOfInt& in, IntVect& out) const
{
  // Cree un tableau parallele structure comme un tableau aux elements
  // du maillage vdf, initialise a zero.

  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, refdomaine_dis_.valeur());
  const Domaine& domaine = domaine_vf.domaine();
  domaine.creer_tableau_elements(out);

  const int nb_elem = domaine.nb_elem();
  const int nb_elem_tot = domaine.nb_elem_tot();
  Cerr << "nb_elem= " << nb_elem << " nb_elem_tot= " << nb_elem_tot;

  // Pour s'assurer que les tableaux ont la meme taille :
  assert(in.size_array() == out.size());

  // La copie veut un pointeur vers le md_vector non initialise.
  //  out.copy(in, /*RESIZE_OPTIONS opt*/ RESIZE_OPTIONS::COPY_INIT);

  // On utilise dont le inject_array
  // mais comment s'assurer qu'on a mis le tableau au bon endroit? et pas un peu
  // dans l'espace virt?
  out.inject_array(in);
  out.echange_espace_virtuel();
}

#define NUM_COMPO_INVALID (-2000000000)

// Ajout d'un terme en d(rho*v)/dt
// delta_rho = rho_liq - rho_vap
// la variable vpoint correspond a d(v)/dt, ou v correspond a nptquel valeur
// vrepul : OWN_PTR(Champ_base) etendu de potentiel de repulsion seul
// vabsrepul : OWN_PTR(Champ_base) etendu de la valeur absolue des repulsions.
void IJK_Interfaces::ajouter_terme_source_interfaces(
  IJK_Field_vector3_double& vpoint,
  IJK_Field_vector3_double& vrepul,
  IJK_Field_vector3_double& vabsrepul
) const
{
  statistiques().begin_count(source_counter_);

  const Domaine_IJK& geom = ref_domaine_.valeur();

  // deux options pour le calcul du terme source interfacial
  // 1 : use_tryggvason_interfacial_source  --> Voir Muradoglu et Tryggvason
  // 2 : !use_tryggvason_interfacial_source --> Discretisation originale de B.Mathieu pour les cas sigma = cte (p.92 et suivantes de sa thèse)
  // on peut aussi faire des mixte : B.Mathieu pour la partie en sigma + terme source de Marangoni avec la formulation de Tryggvason
  // Le terme source de Marangoni nest code que dans lapproche de tryggvason pour le moment

  if (use_tryggvason_interfacial_source_ or !maillage_ft_ijk_.Surfactant_facettes().get_disable_marangoni_source_term())
    {
      double interf1, interf2 ;
      double s1, s2 ;
      Int3 ijk_other_elem;
      for (int icol1 = 0; icol1 < max_authorized_nb_of_components_; icol1++)
        {
          for (int direction = 0; direction < 3; direction++)
            {
              for (int k = 0; k < vpoint[direction].nk(); k++)
                {
                  for (int j = 0; j < vpoint[direction].nj(); j++)
                    {
                      for (int i = 0; i < vpoint[direction].ni(); i++)
                        {
                          // le traitement du terme source de tension de surface doit etre entierement explicit
                          // on a donc besoin de tout stocker au pas de temps precedent dans next() avant le deplacement de linterface
                          // on utilise alors les champ old n (avant deplacement).
                          // les champs _par_compo sont calcules aux elements
                          // on interpole aux faces pour ajouter à vpoint
                          ijk_other_elem[0] = i;
                          ijk_other_elem[1] = j;
                          ijk_other_elem[2] = k;
                          ijk_other_elem[direction] -= 1;
                          s1 = surf_par_compo_[old()][icol1](i,j,k);
                          s2 = surf_par_compo_[old()][icol1](ijk_other_elem[0],ijk_other_elem[1],ijk_other_elem[2]);
                          interf1 = source_interf_par_compo_[old()][direction][icol1](i, j, k);
                          interf2 = source_interf_par_compo_[old()][direction][icol1](ijk_other_elem[0],ijk_other_elem[1],ijk_other_elem[2]);
                          vpoint[direction](i, j, k) += 1./2. * (s1 * interf1 + s2 * interf2);
                        }
                    }
                }
            }
        }
    }
  if (!use_tryggvason_interfacial_source_)
    {
      // calculer la courbure et le terme de gravite aux sommets du maillage
      // lagrangien On appelle ce terme "phi", potentiel aux sommets
      // ducluz : voir p.93 et suivantes dans MathieuPhD
      Int3 ijk_face;
      const int nkmax = std::max(vpoint[DIRECTION_I].nk(), std::max(vpoint[DIRECTION_J].nk(), vpoint[DIRECTION_K].nk()));
      for (int k = 0; k < nkmax; k++)
        {
          for (int direction = 0; direction < 3; direction++)
            {
              if (k >= vpoint[direction].nk())
                continue;

              const double delta_dir = geom.get_constant_delta(direction);
              const int offset = geom.get_offset_local(direction);
              const bool perio = geom.get_periodic_flag(direction);
              Domaine_IJK::Localisation loc = vpoint[direction].get_localisation();
              const int nb_items_tot = geom.get_nb_items_global(loc, direction);
              for (int j = 0; j < vpoint[direction].nj(); j++)
                {
                  for (int i = 0; i < vpoint[direction].ni(); i++)
                    {
                      ijk_face[0] = i;
                      ijk_face[1] = j;
                      ijk_face[2] = k;
                      const int global_face_position = ijk_face[direction] + offset;
                      if (!perio && (global_face_position == 0 || global_face_position == nb_items_tot))
                        continue; // on a wall...

                      Int3 ijk_droite = ijk_face; // l'element de droite a le meme num que la face
                      Int3 ijk_gauche = ijk_face; // l'element de gauche est a l'indice - 1 ...
                      ijk_gauche[direction]--;

                      // Boucle sur les elements a gauche et a droite de la face.
                      // Selon le tour, on appelle un des elem  : elem1
                      for (int gauche_droite = 0; gauche_droite <= 1; gauche_droite++)
                        {
                          // Boucle sur les colonnes de l'elem1 :
                          for (int icol1 = 0; icol1 < max_authorized_nb_of_components_; icol1++)
                            {
                              const Int3& elem1 = gauche_droite ? ijk_droite : ijk_gauche;
                              const Int3& elem2 = gauche_droite ? ijk_gauche : ijk_droite;
                              // Le signe pour le gradient :
                              int signe = gauche_droite ? -1 : 1;
                              int num_compo = compos_traversantes_[old()][icol1](elem1[0], elem1[1], elem1[2]);
                              if (num_compo == NUM_COMPO_INVALID)
                                break;
                              // cette composante est-elle presente sur elem2 et a quelle
                              // colonne ?
                              int icol2;
                              for (icol2 = 0; icol2 < max_authorized_nb_of_components_; icol2++)
                                if (compos_traversantes_[old()][icol2](elem2[0], elem2[1], elem2[2]) == num_compo)
                                  break;

                              if (icol2 < max_authorized_nb_of_components_ && gauche_droite == 1)
                                {
                                  // on a deja traite cette composante lorsqu'on a fait la boucle
                                  // pour gauche_droite = 0
                                  continue;
                                }

                              const double indic = indicatrice_par_compo_[old()][icol1](elem1[0], elem1[1], elem1[2]);


                              const double phi = phi_par_compo_[old()][icol1](elem1[0], elem1[1], elem1[2]);

                              const double surface = surface_par_compo_[old()][icol1](
                                                       elem1[0], elem1[1], elem1[2]
                                                     ); // la compo dans l'elem c'est icol1.

                              const double repul = repuls_par_compo_[old()][icol1](elem1[0], elem1[1], elem1[2]);
                              if (correction_gradient_potentiel_)
                                {
                                  double indic_voisin = 0., phi_voisin = 0.;
                                  double repul_voisin = 0.;
                                  if (icol2 == max_authorized_nb_of_components_)
                                    {
                                      // il n'y a aucune intersection (reelle) par num_compo dans
                                      // l'elem voisin (elem2).

                                      // Determine l'indicatrice dans l'element voisin lorsqu'il
                                      // n'est pas traverse. Pour cela, on se base sur une
                                      // comparaison de la position du centre de la face au plan
                                      // moyen defini par les facettes presentes dans l'element
                                      // elem.
                                      //   Point du plan :     bary_facettes_dans_elem
                                      //   Normale au plan :   normale  (non-unitaire)
                                      //
                                      // Equation du plan :
                                      //   F(X,Y,Z)=X*NX+Y*NY+Z*NZ+ CONSTANTE
                                      //   CONSTANTE = -(NX*PX+NY*PY+NZ*PZ) ou P est un point du
                                      //   plan,
                                      //                                    ici, c'est le cdg des
                                      //                                    facettes dans l'elem.
                                      //
                                      //   F(X,Y,Z) = (X-PX)*NX + (Y-PY)*NY + (Z-PZ)*NZ

                                      // F n'est pas une distance car normale non-unitaire.
                                      // Mais F est signee (positive si on va dans le sens de N, cad
                                      // vers le liquide)

                                      // Coordonnees centre face :
                                      Vecteur3 centre_face = geom.get_coords_of_dof(ijk_face[0], ijk_face[1], ijk_face[2], loc);

                                      Vecteur3 normale;
                                      Vecteur3 bary_facettes_dans_elem;
                                      for (int dir = 0; dir < 3; dir++)
                                        {
                                          const int idx = icol1 * 3 + dir;
                                          // TODO: aym attention, il a des chances que
                                          // normale_par_compo et bary_par_compo_[old()] ne soient pas
                                          // calcules au meme moment que indic_par_compo par exemple.
                                          // Ca risque de poser pb. Il faudra donc sirtir le calcul de
                                          // indic_par_compo de cette methode.
                                          normale[dir] = normale_par_compo_[old()][idx](elem1[0], elem1[1], elem1[2]);
                                          bary_facettes_dans_elem[dir] = bary_par_compo_[old()][idx](elem1[0], elem1[1], elem1[2]);
                                        }
                                      // Calcul du produit scalaire pour savoir de quel cote on est
                                      // :
                                      const double ps = Vecteur3::produit_scalaire(centre_face - bary_facettes_dans_elem, normale);

                                      // Si la fonction F est positive, le voisin est liquide :
                                      indic_voisin = (ps > 0 ? 1. : 0.);
                                      phi_voisin = phi; // On prend le phi de l'elem1 puisque le
                                      // voisin (elem2) n'est pas traverse.
                                      repul_voisin = repul;
                                    }
                                  else
                                    {
                                      // la composante est presente dans l'element voisin
                                      // Le voisin est aussi traverse par num_compo.
                                      phi_voisin = phi_par_compo_[old()][icol2](elem2[0], elem2[1],
                                                                                elem2[2]); // la compo dans le voisin c'est icol2
                                      repul_voisin = repuls_par_compo_[old()][icol2](elem2[0], elem2[1],
                                                                                     elem2[2]); // la compo dans le voisin c'est icol2
                                      indic_voisin = indicatrice_par_compo_[old()][icol2](elem2[0], elem2[1], elem2[2]);
                                    }

                                  // Calcul du gradient a la face :
                                  double gradient_phi_indic = (phi_voisin * indic_voisin - phi * indic) / delta_dir * signe;

                                  vpoint[direction](i, j, k) += gradient_phi_indic;
                                  double gradient_repul_indic = (repul_voisin * indic_voisin - repul * indic) / delta_dir * signe;
                                  vrepul[direction](i, j, k) += gradient_repul_indic;
                                  vabsrepul[direction](i, j, k) += fabs(gradient_repul_indic);
                                }
                              else
                                {
                                  double indic_voisin = 0., phi_face = 0.;
                                  double repul_face = 0.;
                                  if (icol2 == max_authorized_nb_of_components_)
                                    {
                                      // il n'y a aucune intersection (reelle) par num_compo dans
                                      // l'elem voisin (elem2).

                                      // Determine l'indicatrice dans l'element voisin lorsqu'il
                                      // n'est pas traverse. Pour cela, on se base sur une
                                      // comparaison de la position du centre de la face au plan
                                      // moyen defini par les facettes presentes dans l'element
                                      // elem.
                                      //   Point du plan :     bary_facettes_dans_elem
                                      //   Normale au plan :   normale  (non-unitaire)
                                      //
                                      // Equation du plan :
                                      //   F(X,Y,Z)=X*NX+Y*NY+Z*NZ+ CONSTANTE
                                      //   CONSTANTE = -(NX*PX+NY*PY+NZ*PZ) ou P est un point du
                                      //   plan,
                                      //                                    ici, c'est le cdg des
                                      //                                    facettes dans l'elem.
                                      //
                                      //   F(X,Y,Z) = (X-PX)*NX + (Y-PY)*NY + (Z-PZ)*NZ

                                      // F n'est pas une distance car normale non-unitaire.
                                      // Mais F est signee (positive si on va dans le sens de N, cad
                                      // vers le liquide)

                                      // Coordonnees centre face :
                                      Vecteur3 centre_face = geom.get_coords_of_dof(ijk_face[0], ijk_face[1], ijk_face[2], loc);

                                      Vecteur3 normale;
                                      Vecteur3 bary_facettes_dans_elem;
                                      for (int dir = 0; dir < 3; dir++)
                                        {
                                          const int idx = icol1 * 3 + dir;
                                          // TODO: AYM pareil attention a la synchro avec le maillage.
                                          normale[dir] = normale_par_compo_[old()][idx](elem1[0], elem1[1], elem1[2]);
                                          bary_facettes_dans_elem[dir] = bary_par_compo_[old()][idx](elem1[0], elem1[1], elem1[2]);
                                        }
                                      // Calcul du produit scalaire pour savoir de quel cote on est
                                      // :
                                      const double ps = Vecteur3::produit_scalaire(centre_face - bary_facettes_dans_elem, normale);

                                      // Si la fonction F est positive, le voisin est liquide :
                                      indic_voisin = (ps > 0 ? 1. : 0.);
                                      phi_face = phi; // On prend le phi de l'elem1 puisque le
                                      // voisin (elem2) n'est pas traverse.
                                      repul_face = repul;

                                    }
                                  else
                                    {
                                      // la composante est presente dans l'element voisin
                                      // Le voisin est aussi traverse par num_compo.
                                      const double phi_voisin = phi_par_compo_[old()][icol2](elem2[0], elem2[1],
                                                                                             elem2[2]); // la compo dans le voisin c'est icol2

                                      const double repul_voisin = repuls_par_compo_[old()][icol2](elem2[0], elem2[1],
                                                                                                  elem2[2]); // la compo dans le voisin c'est icol2
                                      indic_voisin = indicatrice_par_compo_[old()][icol2](elem2[0], elem2[1], elem2[2]);
#if 0
                                      // Seulement pour debug :
                                      field_indicatrice(elem2[0], elem2[1], elem2[2]) = indic_voisin;
#endif
                                      const double surface_voisin = surface_par_compo_[old()][icol2](
                                                                      elem2[0], elem2[1], elem2[2]
                                                                    ); // la compo dans le voisin c'est icol2
                                      // assert provisoir :

                                      assert(compos_traversantes_[old()][icol2](elem2[0], elem2[1], elem2[2]) == num_compo);
                                      assert(indic_voisin > -0.5);
                                      if (est_egal(surface+surface_voisin,0.))
                                        {
                                          phi_face = 0.;
                                          repul_face = 0.;
                                        }
                                      else
                                        {
                                          // Il faut calculer le phi moyen a la face :
                                          phi_face = (phi * surface + phi_voisin * surface_voisin) / (surface + surface_voisin);
                                          repul_face = (repul * surface + repul_voisin * surface_voisin) / (surface + surface_voisin);
                                        }
                                    }
                                  // Calcul du gradient a la face :
                                  double gradient_indic = (indic_voisin - indic) / delta_dir * signe;


                                  // terme de repulsion
                                  // parcourir les elements voisins jusqu'a une distance d (en
                                  // mailles) si un element voisin contient une composante connexe
                                  // differente, trouver la plus petite distance a cette
                                  // composante


                                  vpoint[direction](i, j, k) += phi_face * gradient_indic ;
                                  vrepul[direction](i, j, k) += repul_face * gradient_indic;
                                  vabsrepul[direction](i, j, k) += fabs(repul_face * gradient_indic);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
  statistiques().end_count(source_counter_);
}

static inline double determinant(const Vecteur3& v1, const Vecteur3& v2, const Vecteur3& v3)
{
  Vecteur3 tmp;
  Vecteur3::produit_vectoriel(v1, v2, tmp);
  return Vecteur3::produit_scalaire(tmp, v3);
}
// A et B sont les points du segment.
// C, D et E les sommets du triangle.
// Retourne false dans les cas pathologiques ou un determinant est nul.
static bool intersection_segment_triangle(const Vecteur3& A,
                                          const Vecteur3& B,
                                          const Vecteur3& C,
                                          const Vecteur3& D,
                                          const Vecteur3& E)
{
  // A et B sont-ils de part et d'autre du plan?
  // Le produit vectoriel det(CD, CE, CA) ne doit pas avoir le meme signe que
  // det(CD, CE, CB)
  Vecteur3 CD, CE, CA, CB;
  CD = D - C;
  CE = E - C;
  CA = A - C;
  CB = B - C;
  const double det1 = determinant(CD, CE, CA);
  const double det2 = determinant(CD, CE, CB);

  if (det1 * det2 >= 0)
    return false; // Le segment ne coupe pas le plan ABC.

  // L'intersection est-elle dans le triangle?
  // Tous les determinants pris dans le meme ordre doivent etre de meme signe :
  Vecteur3 AB, AC, AD, AE;
  AB = B - A;
  AC = C - A;
  AD = D - A;
  AE = E - A;
  const double det3 = determinant(AC, AD, AB);
  const double det4 = determinant(AD, AE, AB);
  const double det5 = determinant(AE, AC, AB);

  if ((det3 * det4 <= 0) || (det3 * det5 <= 0))
    return false; // Le segment coupe le plan ABC mais en dehors du triangle.

  return true;
}

int IJK_Interfaces::compute_cell_phase_with_interface_normal(int num_elem, int direction, int face_plus)
{
  const Domaine_IJK& split = ref_domaine_.valeur();
  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const Intersections_Elem_Facettes& intersections = mesh.intersections_elem_facettes();

  const IntTab& facettes = mesh.facettes();
  const DoubleTab& sommets = mesh.sommets();
  const Vecteur3 cell_size(split.get_constant_delta(DIRECTION_I),
                           split.get_constant_delta(DIRECTION_J),
                           split.get_constant_delta(DIRECTION_K));

  ArrOfInt facettes_traversantes;
  intersections.get_liste_facettes_traversantes(num_elem, facettes_traversantes);
  const int N = facettes_traversantes.size_array();

  if (N == 0)
    return -1;

  // on parcourt les facettes qui traversent l'element num_elem
  int index = 0;
  while (1)
    {
      int fa7 = facettes_traversantes[index];
      const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(fa7);

      // Comme dans compute_list_compo_connex_in_element, on ne parcourt pas les
      // facettes de fraction surface nulle. Note BM: je suis tombe sur un cas ou
      // l'interface a une fraction de 1e-9 dans la facette,
      //  ca produit une erreur ensuite a cause des arrondis, le resultat est
      //  faux, donc je mets
      // une limite a 1e3 pour ne pas utiliser les facettes qui ont une fraction
      // d'intersection trop faible
      if (data.fraction_surface_intersection_ < 1e-2)
        {
          index++;
          if (index == N) // on a parcouru toutes les facettes et aucune ne
            // convient: echec de la methode
            return -1;
          continue;
        }

      FixedVector<Vecteur3, 3> sommets_facette;
      for (int i = 0; i < 3; i++)
        sommets_facette[i] = Vecteur3(sommets, facettes(fa7, i));

      /* calcule du segment AB partant du centre de gravite
       * de l'intersection facette elem (A) vers la projection orthogonale de A
       * sur la face commune avec le voisin qu'on cherche a evaluer (B) */

      // Coordonnees de A (pas barycentriques) :
      Vecteur3 coordA(0., 0., 0.);
      for (int isom = 0; isom < 3; isom++)
        {
          const double bary_som = data.barycentre_[isom];
          for (int dir = 0; dir < 3; dir++)
            coordA[dir] += bary_som * sommets_facette[isom][dir];
        }
      // Coordonnees de B, au lieu de prendre la coordonnee du centre de la face,
      // on prend un point un peu plus loin dans l'element voisin, de toutes
      // facons l'element voisin n'est pas traverse par l'interface donc si on
      // trouve une intersection c'est forcement dans l'element courant.
      Vecteur3 coordB = coordA;
      coordB[direction] += cell_size[direction] * face_plus;

      // on parcourt toutes les facettes de l'element pour voir si l'une d'entre
      // elles intersecte le segment AB
      bool coupe = false;
      for (int index2 = intersections.index_elem()[num_elem]; index2 >= 0;
           index2 = intersections.data_intersection(index2).index_facette_suivante_)
        {
          const Intersections_Elem_Facettes_Data& data2 = intersections.data_intersection(index2);
          const int nvl_fa7 = data2.numero_facette_;
          if ((nvl_fa7 == fa7) || (data2.fraction_surface_intersection_ == 0))
            continue;
          // Les Coords des 3 sommets de la nouvelle facette :
          const Vecteur3 s0(sommets, facettes(nvl_fa7, 0));
          const Vecteur3 s1(sommets, facettes(nvl_fa7, 1));
          const Vecteur3 s2(sommets, facettes(nvl_fa7, 2));

          // Y-a-t-il interstion entre AB et cette facette?
          coupe = intersection_segment_triangle(coordA, coordB, s0, s1, s2);

          if (coupe)
            {
              // il y a une facette entre fa7 et la face de l'elem qui nous interesse
              // on saute toutes les fa7 jusqu'a nvl_fa7
              for (int j = 0; j < N; j++)
                {
                  if (facettes_traversantes[j] == nvl_fa7)
                    index = j;
                }
              break;
            }
        }
      // pas de coupe : on peut determiner la phase de l'element voisin
      // car il n'y a aucun obstacle entre la facete fa7 et l'element voisin
      // on evalue pour cela le vecteur normal a la facette
      if (!coupe)
        {
          Vecteur3 normale_facette;
          Vecteur3::produit_vectoriel(sommets_facette[1] - sommets_facette[0], sommets_facette[2] - sommets_facette[0],
                                      normale_facette);

          if (normale_facette[direction] * face_plus > 0.) // cellule liquide
            return 1;
          else // cellule gazeuse
            return 0;
        }
    }
}

// Ce define active une validation complete:
// pour chaque element pour lequel on peut determiner s'il est liquide ou
// vapeur, on incremente un compteur et on verifie a la fin que tous les
// compteurs sont coherents.

#define COMPUTE_DRAPEAUX_VALIDATION

void IJK_Interfaces::compute_drapeaux_vapeur_v4(const IntVect& vecteur_composantes,
                                                ArrOfInt& drapeau_vapeur) const
{
  static Stat_Counter_Id calculs_drapeaux_counter_ =
    statistiques().new_counter(2, "Calcul de l'indicatrice: calculs des drapeaux");
  statistiques().begin_count(calculs_drapeaux_counter_);

  const Domaine_IJK& split = ref_domaine_.valeur();
  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const Intersections_Elem_Facettes& intersections = mesh.intersections_elem_facettes();
  const ArrOfInt& index_elem = maillage_ft_ijk_.intersections_elem_facettes().index_elem();
  const ArrOfInt& index_facette = maillage_ft_ijk_.intersections_elem_facettes().index_facette();

  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, refdomaine_dis_.valeur());
  const IntTab& elem_faces = domaine_vf.elem_faces();
  const IntTab& faces_voisins = domaine_vf.face_voisins();

#ifdef COMPUTE_DRAPEAUX_VALIDATION
  ArrOfInt composantes_comptage_phase0(drapeau_vapeur.size_array());
  ArrOfInt composantes_comptage_phase1(drapeau_vapeur.size_array());
#endif

  drapeau_vapeur = 0; // Tout vapeur
  const Vecteur3 cell_size(split.get_constant_delta(DIRECTION_I),
                           split.get_constant_delta(DIRECTION_J),
                           split.get_constant_delta(DIRECTION_K));

  const int nb_facettes = mesh.nb_facettes();
  const IntTab& facettes = mesh.facettes();
  const DoubleTab& sommets = mesh.sommets();
  for (int fa7 = 0; fa7 < nb_facettes; fa7++)
    {
      // On parcourt aussi les facettes virtuelles
      // mais on ne parcourt pas les facettes de
      // data.fraction_surface_intersection_ == 0.
      FixedVector<Vecteur3, 3> sommets_facette;
      {
        for (int i = 0; i < 3; i++)
          sommets_facette[i] = Vecteur3(sommets, facettes(fa7, i));
      }
      int index = index_facette[fa7];
      while (index >= 0)
        {
          const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);

          // Comme dans compute_list_compo_connex_in_element, on ne parcours pas les
          // facettes de fraction surface nulle. Note BM: je suis tombe sur un cas
          // ou l'interface a une fraction de 1e-9 dans la facette,
          //  ca produit une erreur ensuite a cause des arrondis, le resultat est
          //  faux, donc je mets
          // une limite a 1e3 pour ne pas utiliser les facettes qui ont une fraction
          // d'intersection trop faible
          if (data.fraction_surface_intersection_ < 1e-2)
            {
              /*    Cerr << "compute_drapeaux_vapeur_v4 : elem= " <<
                 data.numero_element_
                   << " Intersection issue du parcours : numero_facette_= " <<
                 data.numero_facette_
                   << " numero_element_= "  << data.numero_element_
                   << " fraction= " << data.fraction_surface_intersection_  << finl;
               */
              index = data.index_element_suivant_;
              continue; // On zappe cette facette...
            }

          const int elem = data.numero_element_;
          /* boucle sur les faces de cet element */

          for (int iface = 0; iface < 6; iface++)
            {
              const int direction = iface % 3;
              const int face_plus = (iface > 2) ? 1 : -1; // +1 si on est sur la face de droite de l'element

              const int num_face = elem_faces(elem, iface);
              const int elem_voisin = faces_voisins(num_face, 0) + faces_voisins(num_face, 1) - elem;
              // Si on est au bord du domaine, ne pas faire la suite:
              if (elem_voisin < 0)
                {
                  continue;
                }
              // Si l'element voisin n'est pas traverse par une facette (sinon,
              // global_compo=-1) et que sa compo n'est pas deja marquee :
              const int global_compo_voisin = vecteur_composantes[elem_voisin];

              // Si l'element voisin est traverse par une interface, ne pas faire la
              // suite
              if (global_compo_voisin < 0)
                continue;

#ifndef COMPUTE_DRAPEAUX_VALIDATION
              // s'il est deja marque "1" ne pas faire la suite (sauf si validation de
              // l'algorihme)
              if (drapeau_vapeur[global_compo_voisin] == 1)
                continue;
#endif
              Vecteur3 normale_facette;
              Vecteur3::produit_vectoriel(sommets_facette[1] - sommets_facette[0],
                                          sommets_facette[2] - sommets_facette[0],
                                          normale_facette);

#ifndef COMPUTE_DRAPEAUX_VALIDATION
              // Si la normale a la facette est telle que l'element voisin doit etre
              // d'indicatrice 0, ne pas faire la suite (sauf si validation de
              // l'algorithme)
              if (normale_facette[direction] * face_plus <= 0.)
                continue;
#endif
              /* calcule du segment AB partant du centre de gravite
                 de l'intersection facette elem (A) vers la projection de A
                 sur la face orthogonalement a la face (B) */

              // Coordonnees de A (pas barycentriques) :
              Vecteur3 coordA(0., 0., 0.);
              for (int isom = 0; isom < 3; isom++)
                {
                  const double bary_som = data.barycentre_[isom];
                  for (int dir = 0; dir < 3; dir++)
                    coordA[dir] += bary_som * sommets_facette[isom][dir];
                }

              // Coordonnees de B, au lieu de prendre la coordonnee du centre de la
              // face, on prend un point un peu plus loin dans l'element voisin, de
              // toutes facons l'element voisin n'est pas traverse par l'interface
              // donc si on trouve une intersection c'est forcement dans l'element
              // courant.
              Vecteur3 coordB = coordA;
              coordB[direction] += cell_size[direction] * face_plus;

              /* si le segment AB coupe une autre facette dans l'element
                 ou s'il y a une ambiguite (un determinant nul dans le calcul)
                 on ne peut pas determiner le signe de l'indicatrice chez le voisin */
              // Boucle sur les facettes dans l'element :
              bool coupe = false;
              // Une facon differente d'ecrire la boucle while, permet de ne pas
              // dupliquer la ligne data2.index_facette_suivante_
              for (int index2 = index_elem[elem]; index2 >= 0;
                   index2 = intersections.data_intersection(index2).index_facette_suivante_)
                {
                  const Intersections_Elem_Facettes_Data& data2 = intersections.data_intersection(index2);
                  const int nvl_fa7 = data2.numero_facette_;
                  if ((nvl_fa7 == fa7) || (data2.fraction_surface_intersection_ == 0))
                    continue;
                  // Les Coords des 3 sommets de la nouvelle facette :
                  const Vecteur3 s0(sommets, facettes(nvl_fa7, 0));
                  const Vecteur3 s1(sommets, facettes(nvl_fa7, 1));
                  const Vecteur3 s2(sommets, facettes(nvl_fa7, 2));

                  // Y-a-t-il interstion entre AB et cette facette?
                  coupe = intersection_segment_triangle(coordA, coordB, s0, s1, s2);

                  if (coupe)
                    break; // Il y a une facette (nvl_fa7) entre fa7 et l'element
                  // voisin, on ne pourra pas determiner l'indicatrice
                }

              // Il n'y a pas d'autre facette entre "fa7" et l'element voisin, on peut
              // donc determiner le signe de l'indicatrice dans l'element voisin
              if (!coupe)
                {
                  if (normale_facette[direction] * face_plus > 0.)
                    drapeau_vapeur[global_compo_voisin] = 1;
#ifdef COMPUTE_DRAPEAUX_VALIDATION
                  // Incremente les compteurs pour la phase trouvee:
                  if (normale_facette[direction] * face_plus > 0.)
                    {
                      composantes_comptage_phase1[global_compo_voisin]++;
                      // Arrete tout de suite si c'est incoherent (on a trouve
                      // precedemment que ca doit etre phase 0)
                      // if (composantes_comptage_phase0[global_compo_voisin]>0) {
                      // Cerr << "Erreur ! " << finl;
                      //}
                      // assert(composantes_comptage_phase0[global_compo_voisin] == 0);
                    }
                  if (normale_facette[direction] * face_plus < 0.)
                    {
                      composantes_comptage_phase0[global_compo_voisin]++;
                      // Arrete tout de suite si c'est incoherent (on a trouve
                      // precedemment que ca doit etre phase 1)
                      // if (composantes_comptage_phase1[global_compo_voisin]>0) {
                      // Cerr << "Erreur ! " << finl;
                      //}
                      // assert(composantes_comptage_phase1[global_compo_voisin] == 0);
                    }
#endif
                }
            }
          // Boucle while externe.
          index = data.index_element_suivant_;
        }
    }

  // Synchroniser les marqueurs sur tous les processeurs
  mp_max_for_each_item(drapeau_vapeur);

  // Verifier que les comptages sont coherents entre les differents processeurs:
#ifdef COMPUTE_DRAPEAUX_VALIDATION
  {
    const int n = composantes_comptage_phase0.size_array();
    mp_sum_for_each_item(composantes_comptage_phase0);
    mp_sum_for_each_item(composantes_comptage_phase1);
    if (Process::je_suis_maitre())
      {
        int problem = 0;
        for (int i = 0; i < n; i++)
          {
            // si on a trouve a la fois phase 0 et phase 1 pour une composante:
            // probleme si on a trouve aucune phase pour une composante: probleme,
            // phase indeterminee
            if ((composantes_comptage_phase0[i] > 0 && composantes_comptage_phase1[i] > 0) ||
                (composantes_comptage_phase0[i] == 0 && composantes_comptage_phase1[i] == 0))
              {
                problem += 1;
                if (problem == 1)
                  {
                    Cerr << "[WARNING-Indic-Vote] compo/#phase0/#phase1 ";
                  }
                Cerr << i
                     << "/" << composantes_comptage_phase0[i]
                     << "/" << composantes_comptage_phase1[i]
                     << " ";

                // vote majoritaire:
                if (composantes_comptage_phase0[i] > composantes_comptage_phase1[i])
                  drapeau_vapeur[i] = 0;
                else
                  drapeau_vapeur[i] = 1;
              }
          }
        if (problem > 0)
          Cerr << "End." << finl;
      }
    envoyer_broadcast(drapeau_vapeur, 0);
  }
#endif
  statistiques().end_count(calculs_drapeaux_counter_);
}

static inline double norme_carre(const Vecteur3& x)
{
  return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
}

inline Vecteur3 operator*(double x, const Vecteur3& v)
{
  return Vecteur3(x * v[0], x * v[1], x * v[2]);
}
inline Vecteur3 operator+(const Vecteur3& v, const Vecteur3& w)
{
  return Vecteur3(v[0] + w[0], v[1] + w[1], v[2] + w[2]);
}

static inline double calculer_carre_distance_sommet_facette(const Vecteur3& coord,
                                                            const Maillage_FT_IJK& mesh,
                                                            int num_facette)
{
  const IntTab& facettes = mesh.facettes();
  const DoubleTab& sommets = mesh.sommets();
  Vecteur3 s0(sommets, facettes(num_facette, 0));
  Vecteur3 s1(sommets, facettes(num_facette, 1));
  Vecteur3 s2(sommets, facettes(num_facette, 2));

  // Il faut coder une algorithme d'Uzawa.
  // Met de cote, pour l'instant, je prends juste la distance min
  // entre coord et les trois sommets.

  double d = 1e10;
  for (double x = 0.; x < 1.01; x += 0.2)
    for (double y = 0.; y < 1. - x + 0.01; y += 0.2)
      {
        Vecteur3 s = s0 + x * (s1 - s0) + y * (s2 - s0);
        d = std::min(d, norme_carre(s - coord));
      }
  return d;
}

static void fill_relative_velocity(const DoubleTab& vinterp_tmp, const DoubleTab& vinterp, const IntTab& facettes, int id_facette, int som, DoubleTab& vr_to_other)
{
  const double un_tiers = 1. / 3.;
  if (id_facette == -1)
    {
      // Noone else is found in the given neighbourhood... default value for vr
      // is set to zero... why not?
      for (int idir = 0; idir < 3; idir++)
        vr_to_other(som, idir) = 0.;
    }
  else
    {
      // indexes of the 3 vertices of the facette:
      const int isom0 = facettes(id_facette, 0);
      const int isom1 = facettes(id_facette, 1);
      const int isom2 = facettes(id_facette, 2);
      for (int idir = 0; idir < 3; idir++)
        {
          // Carreful, one is an index in the list (i), whereas "isomN" are
          // indexes in the mesh!
          const double velocity_me = vinterp_tmp(som, idir);
          const double velocity_other =
            un_tiers  * (vinterp(isom0, idir) + vinterp(isom1, idir) + vinterp(isom2, idir));
          vr_to_other(som, idir) = velocity_me - velocity_other;
        }
    }
}

// Warning : sizes are for nb_sommets in the list sommets_a_tester, not for
// mesh!! vr_to_other : the relative velocity to the closest marker.
void IJK_Interfaces::calculer_distance_autres_compo_connexe_octree(const DoubleTab& sommets_a_tester,
                                                                   const ArrOfInt& compo_connexe_sommets,
                                                                   const DoubleTab& vinterp_tmp,
                                                                   const Maillage_FT_IJK& mesh,
                                                                   ArrOfDouble& distance,
                                                                   DoubleTab& vr_to_other,
                                                                   const double distmax)
{
  // Construction d'un octree avec les facettes du maillage
  Octree_Double octree;
  const IntTab& facettes = mesh.facettes();
  octree.build_elements(mesh.sommets(), facettes,
                        0. /* epsilon */, 0 /* utiliser dimension(0) pas dimension_tot(0),
                 ca n'a pas d'importance ici */);

  ArrOfInt liste_facettes;
  const ArrOfInt& compo_connexe_facettes = mesh.compo_connexe_facettes();

  const int nb_som = sommets_a_tester.dimension(0);
  distance.resize_array(nb_som, RESIZE_OPTIONS::NOCOPY_NOINIT);
  distance = distmax;
  vr_to_other.resize(nb_som, 3, RESIZE_OPTIONS::NOCOPY_NOINIT);
  vr_to_other = -1e5; // Invalid value
  assert(vinterp_tmp.dimension(0) == nb_som);
  for (int i = 0; i < nb_som; i++)
    {
      Vecteur3 coord(sommets_a_tester, i);
      const int compo_connexe_som = compo_connexe_sommets[i];

      double dmin = distmax * distmax;

      // Recherche dans l'octree des facettes proches de ce point:
      octree.search_elements_box(coord[0] - distmax, coord[1] - distmax, coord[2] - distmax, coord[0] + distmax,
                                 coord[1] + distmax, coord[2] + distmax, liste_facettes);

      // pour chaque facette, calcule la distance entre la facette et coord, puis
      // prend le min:
      const int nliste = liste_facettes.size_array();
      int idx_facette = -1;
      for (int j = 0; j < nliste; j++)
        {
          const int num_facette = liste_facettes[j];
          const int compo = compo_connexe_facettes[num_facette];
          if (compo == compo_connexe_som)
            continue; // cette facette appartient a la meme compo connexe que le
          // sommet,
          // ne m'interesse pas
          const double d = calculer_carre_distance_sommet_facette(coord, mesh, num_facette);
          if (d < dmin)
            {
              dmin = d;
              idx_facette = num_facette; // On retient le numero de la facette
            }
          // dmin = std::min(d, dmin);
        }
      distance[i] = std::min(distance[i], sqrt(dmin));
      fill_relative_velocity(vinterp_tmp,vinterp_, facettes, idx_facette, i, vr_to_other);

    }
}


// search for the center of mass of a neighbouring face among all cells at distance n from the current cell IJK (in a given direction)
static void check_neighbouring_layer_in_one_direction(int dir0, int dir1, int dir2, int n,
                                                      const Int3& nb_elem_loc, const Int3& ijk,
                                                      const std::map<std::array<int,3>, std::set<int>>& bary_ijk_loc, const DoubleTab& bary,
                                                      const ArrOfInt& compo_connexe_facettes, const int compo_connexe_som,
                                                      const double x, const double y, const double z,
                                                      double& distance, int& id_facette )
{

  for(int sens=0; sens<2; sens++)
    {
      int zero = 0;
      int a = sens == 0 ? std::max(ijk[dir0]-n,zero) : std::min(ijk[dir0]+n,nb_elem_loc[dir0]);
      for(int b=std::max(ijk[dir1]-n,zero); b<std::min(ijk[dir1]+n,nb_elem_loc[dir1]); b++)
        for(int c=std::max(ijk[dir2]-n,zero); c<std::min(ijk[dir2]+n,nb_elem_loc[dir2]); c++)
          {
            std::array<int,3> current_ijk;
            current_ijk[0] = a;
            current_ijk[1] = b;
            current_ijk[2] = c;
            if(bary_ijk_loc.find(current_ijk)!=bary_ijk_loc.end())
              {
                std::set<int> bary_in_current_ijk = bary_ijk_loc.at(current_ijk);
                for(const auto b_fa7: bary_in_current_ijk)
                  {
                    if(compo_connexe_facettes[b_fa7] != compo_connexe_som)
                      {
                        //computing distance between the center of mass of this face and the vertice I'm searching the closest neighbour of
                        double dist = (bary(b_fa7,0)-x)*(bary(b_fa7,0)-x) + (bary(b_fa7,1)-y)*(bary(b_fa7,1)-y) + (bary(b_fa7,2)-z)*(bary(b_fa7,2)-z);
                        distance = std::min(sqrt(dist),distance);
                        id_facette = b_fa7;
                      }
                  }
              }
          }
    }

}

/* @brief For each vertex of the front mesh, compute the distance and the relative velocity to the nearest face belonging to another bubble based on the IJK discretization :
 * starting from the IJK-cell containing my vertex, we progressively search the neighbouring layers until we find the center of mass of a face of another bubble
 * (or we reach the border of the domain)
 * Exemple :
 * checking layer 1:     checking layer 2:        checking layer 3 (and we stop here, as we found the barycenter of a neighbouring face)
 *                                                |x|x|x|x|x|x|x|
 *                        |x|x|x|x|x|             |x|_|_|_|_|_|x|
 *  |x|x|x|               |x|_|_|_|x|             |x|_|_|_|_|_|x|
 *  |x|o|x|               |x|_|o|_|x|             |x|_|_|o|_|_|x|
 *  |x|x|x|               |x|_|_|_|x|             |x|_|_|_|_|_|x|
 *                        |x|x|x|x|x|             |x|x|x|x|x|x|b|
 *
 * WARNING: here, we compute the distance between a vertice and a face as the distance between the vertice and the center of mass of the face.
 * In the method calculer_distance_autres_compo_connexe_octree, the distance taken into account is the shortest distance between the vertex
 * and a set of points on the surface of the face. So the results of the two algorithms may differ
 * WARNING: this algorithm is not validated and is slower than calculer_distance_autres_compo_connexe_octree. To be used only if you encounter a memory problem with the other method
 */
void IJK_Interfaces::calculer_distance_autres_compo_connexe_ijk(const DoubleTab& sommets_a_tester,
                                                                const ArrOfInt& compo_connexe_sommets,
                                                                const DoubleTab& vinterp_tmp,
                                                                const Maillage_FT_IJK& mesh,
                                                                ArrOfDouble& distance,
                                                                DoubleTab& vr_to_other,
                                                                const double distmax)
{
  const IntTab& facettes = mesh.facettes();
  const ArrOfInt& compo_connexe_facettes = mesh.compo_connexe_facettes();
  const int nb_som = sommets_a_tester.dimension(0);
  const int nb_fa7 = facettes.dimension(0);

  const Domaine_IJK& splitting = I_ft().get_domaine();
  Int3 ijk_glob, ijk_loc, useless;
  Int3 nb_elem_loc;
  for(int dir=0; dir<3; dir++)
    nb_elem_loc[dir] = splitting.get_nb_elem_local(dir);

  distance.resize_array(nb_som, RESIZE_OPTIONS::NOCOPY_NOINIT);
  distance = distmax;
  vr_to_other.resize(nb_som, 3, RESIZE_OPTIONS::NOCOPY_NOINIT);
  vr_to_other = -1e5; // Invalid value
  assert(vinterp_tmp.dimension(0) == nb_som);

  DoubleTab bary(nb_fa7,3);
  // for each i,j,k cell of my local domain, list of all the centers of mass that it contains
  std::map<std::array<int,3>, std::set<int>> bary_ijk_loc;
  for (int fa7=0; fa7<nb_fa7; fa7++)
    {
      // computing simple center of mass for each faces
      int s0 = facettes(fa7,0), s1 = facettes(fa7,1), s2 = facettes(fa7,2);
      for(int dir=0; dir<3; dir++)
        bary(fa7,dir) = (sommets_a_tester(s0,dir) + sommets_a_tester(s1,dir) + sommets_a_tester(s2,dir)) / 3.;

      splitting.search_elem(bary(fa7,0), bary(fa7,1), bary(fa7,2), ijk_glob, ijk_loc, useless);
      std::array<int,3> ijk;
      for(int dir=0; dir<3; dir++) ijk[dir] = ijk_loc[dir];
      bary_ijk_loc[ijk].insert(fa7);
    }

  for (int som=0; som<nb_som; som++)
    {
      double x=sommets_a_tester(som,0), y=sommets_a_tester(som,1), z=sommets_a_tester(som,2);
      splitting.search_elem(x, y, z, ijk_glob, ijk_loc, useless);
      const int compo_connexe_som = compo_connexe_sommets[som];

      bool local_domain_checked = false;
      int n_layer = 0;
      int id_facette = -1; //id of the closest face found
      double& dist = distance[som];

      // we stop searching for the closest neighbour if we found one, or if the entire local domain has been covered
      while(!local_domain_checked && id_facette==-1)
        {
          // checking left and right neighbourhood of my cell
          check_neighbouring_layer_in_one_direction(0 /*fixed dir*/, 1, 2, n_layer,
                                                    nb_elem_loc, ijk_loc,
                                                    bary_ijk_loc, bary,
                                                    compo_connexe_facettes, compo_connexe_som,
                                                    x, y, z,
                                                    dist, id_facette );

          // check for up and down neighbourhood
          check_neighbouring_layer_in_one_direction(1 /*fixed dir*/, 0, 2, n_layer,
                                                    nb_elem_loc, ijk_loc,
                                                    bary_ijk_loc, bary,
                                                    compo_connexe_facettes, compo_connexe_som,
                                                    x, y, z,
                                                    dist, id_facette );



          // check for front and back neighbourhood
          check_neighbouring_layer_in_one_direction(2 /*fixed dir*/, 0, 1, n_layer,
                                                    nb_elem_loc, ijk_loc,
                                                    bary_ijk_loc, bary,
                                                    compo_connexe_facettes, compo_connexe_som,
                                                    x, y, z,
                                                    dist, id_facette );


          // have we checked the whole local domain yet ?
          int i=ijk_loc[0], j=ijk_loc[1], k=ijk_loc[2];
          local_domain_checked = i-n_layer <= 0 && i+n_layer>= nb_elem_loc[0]
                                 && j-n_layer <= 0 && j+n_layer>= nb_elem_loc[1]
                                 && k-n_layer <= 0 && k+n_layer>= nb_elem_loc[2];

          // next layer
          n_layer++;
        }

      // Once the closest neighbour is found, fill velocity
      fill_relative_velocity(vinterp_tmp,vinterp_, facettes, id_facette, som, vr_to_other);

    }
}

// Cette fonction cherche dans les coordonnees les sommets qui sont au bord du
// domaine dans les directions
//    >= dir, les envoie aux voisins concernes et recupere la distance.
// Pour cela on
//  cherche les sommets a envoyer aux voisins, qu'on envoir et qu'on ajoute a la
//  liste coord_sommets appelle recursivement la fonction avec dir+1 recupere
//  les resultats des sommets envoyes aux voisins, et remet coord_sommets dans
//  son etat initial
// Si dir=3, on calcule la distance avec la liste de sommets coord_sommets, fin
// de la recursion Attention, recursion inhabituelle: la methode modifie le
// parametre passe par adresse et ajoute
//  des sommets a chaque appel recursif. Elle les retire a nouveau avant de
//  sortir.
// (on aurait pu passer un tableau "const" et le dupliquer pour avoir une copie
// de travail).
//
// Cette facon de faire permet de faire tous les echanges avec seulement 12
// messages echanges (6 a l'aller et 6 au retour). Si on veut envoyer
// directement tous les sommets a tous les processeurs, il y a 26 voisins
// (coins) donc 52 messages (aller et retour).
void IJK_Interfaces::recursive_calcul_distance_chez_voisin(DoubleTab& vinterp_tmp, int dir,
                                                           const Maillage_FT_IJK& mesh, DoubleTab& coord_sommets,
                                                           ArrOfInt& compo_sommet, ArrOfDouble& distance,
                                                           DoubleTab& vr_to_other, double distmax)
{
  const Domaine_IJK& splitting = ref_domaine_.valeur();
  if (dir == 3)
    {
      if(!no_octree_method_)
        calculer_distance_autres_compo_connexe_octree(coord_sommets, compo_sommet, vinterp_tmp, mesh, distance, vr_to_other, distmax);
      else
        calculer_distance_autres_compo_connexe_ijk(coord_sommets, compo_sommet, vinterp_tmp, mesh, distance, vr_to_other, distmax);
    }
  else
    {
      // schema pour envoyer recevoir des donnees de mes voisins gauche droite
      Schema_Comm schema;
      ArrOfInt pe_list;
      ArrOfInt flags(Process::nproc());
      flags = 0;
      double min_coord, max_coord;
      const int processor_at_left = splitting.get_neighbour_processor(0 /* previous */, dir);
      if (processor_at_left < 0)
        {
          min_coord = -1e10;
        }
      else
        {
          const ArrOfDouble& coord_nodes = splitting.get_node_coordinates(dir);
          const int offset = splitting.get_offset_local(dir);
          const double left_node = coord_nodes[offset];
          min_coord = left_node + distmax;
          if(!flags[processor_at_left])
            {
              flags[processor_at_left] = 1;
              pe_list.append_array(processor_at_left);
            }
        }
      const int processor_at_right = splitting.get_neighbour_processor(1 /* next */, dir);
      if (processor_at_right < 0)
        {
          max_coord = 1e10;
        }
      else
        {
          const ArrOfDouble& coord_nodes = splitting.get_node_coordinates(dir);
          const int offset = splitting.get_offset_local(dir);
          const int i_last_node = offset + splitting.get_nb_elem_local(dir);
          const double right_node = coord_nodes[i_last_node];
          max_coord = right_node - distmax;
          if(!flags[processor_at_right])
            {
              flags[processor_at_right] = 1;
              pe_list.append_array(processor_at_right);
            }
        }
      schema.set_send_recv_pe_list(pe_list, pe_list);

      schema.begin_comm();
      ArrOfInt index_sent_to_left;
      ArrOfInt index_sent_to_right;
      if (processor_at_left >= 0 || processor_at_right >= 0)
        {
          // Bypass in sequential or if no splitting in this direction
          const int nb_som = coord_sommets.dimension(0);
          for (int i_som = 0; i_som < nb_som; i_som++)
            {
              Vecteur3 coord(coord_sommets, i_som);
              if (coord[dir] < min_coord)
                {
                  // envoyer ce sommet au voisin de gauche
                  // je stocke en meme temps l'indice i_som de ce sommet dans le tableau d'origine
                  // (pour mettre a jour la distance avec ce que m'aura envoye le processeur voisin)
                  index_sent_to_left.append_array(i_som);
                  schema.send_buffer(processor_at_left) << compo_sommet[i_som]
                                                        << coord_sommets(i_som, 0)
                                                        << coord_sommets(i_som, 1)
                                                        << coord_sommets(i_som, 2)
                                                        << vinterp_tmp(i_som, 0)
                                                        << vinterp_tmp(i_som, 1)
                                                        << vinterp_tmp(i_som, 2);
                }
              if (coord[dir] > max_coord)
                {
                  // envoyer ce sommet au voisin de droite
                  index_sent_to_right.append_array(i_som);
                  schema.send_buffer(processor_at_right) << compo_sommet[i_som]
                                                         << coord_sommets(i_som, 0)
                                                         << coord_sommets(i_som, 1)
                                                         << coord_sommets(i_som, 2)
                                                         << vinterp_tmp(i_som, 0)
                                                         << vinterp_tmp(i_som, 1)
                                                         << vinterp_tmp(i_som, 2);
                }
            }
        }

      schema.echange_taille_et_messages();

      const int init_count = coord_sommets.dimension(0); // Nombre de sommets initialement recus par la fonction
      int count_left = 0; // nb sommets recus de gauche et de droite
      int count_right = 0;
      for (int left_right = 0; left_right < 2; left_right++)
        {
          const int pe_voisin = (left_right==0) ? processor_at_left : processor_at_right;
          int& count = (left_right==0) ? count_left : count_right;  // attention, c'est une reference, on va modifier la var.
          if (pe_voisin >= 0)
            {
              Entree& buf = schema.recv_buffer(pe_voisin);
              while(1)
                {
                  int num_compo;
                  buf >> num_compo;
                  if (buf.eof())
                    break;
                  double x, y, z, vx, vy, vz;
                  buf >> x >> y >> z >> vx >> vy >> vz;
                  coord_sommets.append_line(x, y, z);
                  vinterp_tmp.append_line(vx, vy, vz);
                  compo_sommet.append_array(num_compo);
                  count++;
                }
            }
        }
      // dans les tableaux on a maintenant init_count sommets venant de la fonction qui m'a appelee,
      // plus count_left sommets qui viennent de ma gauche, plus count_right sommets qui viennent de ma droite
      schema.end_comm();

      // appel recursif avec la direction suivante
      recursive_calcul_distance_chez_voisin(vinterp_tmp, dir + 1, mesh, coord_sommets, compo_sommet, distance, vr_to_other, distmax);

      // on recupere la distance pour tous les sommets du tableau,
      // il faut envoyer la distance aux voisins pour les sommets apres init_count:

      int index = init_count; // indice du premier sommet qui a ete recu d'un voisin
      schema.begin_comm();
      for (int left_right = 0; left_right < 2; left_right++)
        {
          const int pe_voisin = (left_right==0) ? processor_at_left : processor_at_right;
          const int count = (left_right==0) ? count_left : count_right;
          if (pe_voisin >= 0)
            {
              Sortie& buf = schema.send_buffer(pe_voisin);
              for (int i = 0; i < count; i++)
                buf << distance[index+i]
                    << vr_to_other(index+i,0)
                    << vr_to_other(index+i,1)
                    << vr_to_other(index+i,2);
              index += count;
            }
        }
      schema.echange_taille_et_messages();
      for (int left_right = 0; left_right < 2; left_right++)
        {
          int pe_voisin = (left_right==0) ? processor_at_left : processor_at_right;
          const ArrOfInt indices_sommets = (left_right==0) ? index_sent_to_left : index_sent_to_right;
          if (pe_voisin >= 0)
            {
              Entree& buf = schema.recv_buffer(pe_voisin);
              const int count = indices_sommets.size_array();
              for (int i = 0; i < count; i++)
                {
                  double d;
                  buf >> d;
                  double vx, vy, vz;
                  buf >> vx >> vy >> vz;
                  int idx = indices_sommets[i];
                  distance[idx] = std::min(distance[idx], d);
                  vr_to_other(idx,0) = vx;
                  vr_to_other(idx,1) = vy;
                  vr_to_other(idx,2) = vz;
                }
            }
        }
      schema.end_comm();

      // Remet les tableaux dans leur etat d'origine
      coord_sommets.resize(init_count, 3);
      vinterp_tmp.resize(init_count, 3); // GB. Heu.. je ne comprends pas l'utilite avec ma comprehension limitee de l'octree, mais bon.
      compo_sommet.resize_array(init_count);
      // On pourrait verifier avec un assert que coord_sommets est identique a la valeur qu'il avait a l'entree.
      distance.resize_array(init_count);
      vr_to_other.resize(init_count, 3);   // GB. Heu..
    }
}

// Remplit le tableau distance qui, pour chaque sommet du maillage mesh,
// donne la distance entre le point et la composante connexe differente la plus proche
// Si cette distance est superieure a distmax, on ne cherche pas, on prend distmax.

void IJK_Interfaces::calculer_distance_autres_compo_connexe2(ArrOfDouble& distance,
                                                             DoubleTab& v_closer)
{
  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const double distmax=	portee_force_repulsion_;
  //statistiques().begin_count(cnt_CalculerDistance);
  ArrOfIntFT compo_connexe_sommets;
  mesh.calculer_compo_connexe_sommets(compo_connexe_sommets);
  DoubleTab tmp_sommets = mesh.sommets();
  compute_vinterp(); // To be sure the velocity of markers is up-to-date

  // vinterp_tmp is a resizable list that will be as tmp_sommets, but containing velocities.
  DoubleTab vinterp_tmp(vinterp_);
  // Calcul de la distance avec les interfaces locales et autres processeurs
  recursive_calcul_distance_chez_voisin(vinterp_tmp,
                                        0, /* direction I, premiere etape de la recursion */
                                        mesh,
                                        tmp_sommets,
                                        compo_connexe_sommets,
                                        distance,
                                        v_closer,
                                        distmax);
  //statistiques().end_count(cnt_CalculerDistance);
}

// Methodes outils permettant depuis GDB d'ecrire des fichiers tracables dans gnuplot
// Le nom du fichier de sortie est mis en dur dans le code:
void dump_facette(const Maillage_FT_Disc& mesh, int fa7, int append = 0)
{
  SFichier f;
  if (append)
    {
      f.ouvrir("dump_facette.txt", ios::app);
      f << finl << finl;
    }
  else
    {
      f.ouvrir("dump_facette.txt");
    }
  for (int i = 0; i < 4; i++)
    {
      int som = mesh.facettes()(fa7, i%3);
      f << mesh.sommets()(som, 0) << " " << mesh.sommets()(som, 1) << " " << mesh.sommets()(som, 2) << finl;
    }
}

void dump_elem(const Domaine_IJK& split, int ii, int jj, int kk, int append = 0)
{
  SFichier f;
  if (append)
    {
      f.ouvrir("dump_elem.txt", ios::app);
      f << finl << finl;
    }
  else
    {
      f.ouvrir("dump_elem.txt");
    }
  Vecteur3 p[2];
  Int3 ijk;
  ijk[0] = ii;
  ijk[1] = jj;
  ijk[2] = kk;
  for (int i = 0; i < 3; i++)
    {
      int j = split.get_offset_local(i) + ijk[i];
      p[0][i] = split.get_node_coordinates(i)[j];
      p[1][i] = split.get_node_coordinates(i)[j+1];
    }
  for (int dir = 0; dir < 3; dir++)
    {
      for (int seg = 0; seg < 8; seg++)
        {
          int pstart = seg;
          int pend   = seg + (1 << dir);
          if (pend < 8)
            {
              f << p[pstart & 1][0] << " " << p[(pstart >> 1) & 1][1] << " " << p[(pstart >> 2) & 1][2] << finl;
              f << p[pend & 1][0]   << " " << p[(pend >> 1) & 1][1]   << " " << p[(pend >> 2) & 1][2] << finl;
              f << finl << finl;
            }
        }
    }
}

void dump_elem(const Domaine_VF& domaine, int elem, int append = 0)
{
  SFichier f;
  if (append)
    {
      f.ouvrir("dump_elem.txt", ios::app);
      f << finl << finl;
    }
  else
    {
      f.ouvrir("dump_elem.txt");
    }
  Vecteur3 p[2];
  for (int i = 0; i < 3; i++)
    {
      p[0][i] = domaine.xv(domaine.elem_faces(elem, i), i);
      p[1][i] = domaine.xv(domaine.elem_faces(elem, i+3), i);
    }
  for (int dir = 0; dir < 3; dir++)
    {
      for (int seg = 0; seg < 8; seg++)
        {
          int pstart = seg;
          int pend   = seg + (1 << dir);
          if (pend < 8)
            {
              f << p[pstart & 1][0] << " " << p[(pstart >> 1) & 1][1] << " " << p[(pstart >> 2) & 1][2] << finl;
              f << p[pend & 1][0]   << " " << p[(pend >> 1) & 1][1]   << " " << p[(pend >> 2) & 1][2] << finl;
              f << finl << finl;
            }
        }
    }
}


//Renvoie 1 si l option GRAVITE_RHO_G est activee 0 sinon
int IJK_Interfaces::is_terme_gravite_rhog() const
{
  if (terme_gravite_ == GRAVITE_RHO_G)
    return 1;
  else
    return 0;
}

void IJK_Interfaces::detecter_et_supprimer_rejeton(bool duplicatas_etaient_presents)
{
  // Selon ou cette methode est appelee, il faut qu'elle recree les duplicatas.

  if (duplicatas_etaient_presents)
    {
      supprimer_duplicata_bulles();
      maillage_ft_ijk_.parcourir_maillage();
    }
  Maillage_FT_IJK& maillage = maillage_ft_ijk_;
  const ArrOfInt& compo_std = maillage.compo_connexe_facettes();
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  const int nb_facettes=maillage.nb_facettes();
  //dumplata_newtime("DNS2.lata",0.);
  //dumplata_ft_mesh("DNS2.lata", "INTERFACES", maillage_ft_ijk_, 0);
  //dumplata_ft_field("DNS2.lata", "INTERFACES", "COMPO_CONNEXE_AVT", "ELEM",  maillage_ft_ijk_.compo_connexe_facettes(), 0);
  ArrOfInt compo_new(nb_facettes); // Init a zero
  int n = search_connex_components_local_FT(maillage, compo_new);
  //dumplata_ft_field("DNS2.lata", "INTERFACES", "COMPO_CONNEXE_MLX", "ELEM",  compo_new,  0);
  int nb_compo_tot=compute_global_connex_components_FT(maillage, compo_new, n);
  //dumplata_ft_field("DNS2.lata", "INTERFACES", "COMPO_CONNEXE_FIN", "ELEM",  compo_new, 0);
  if (nb_compo_tot>nb_bulles_reelles_)
    {
      int nb_rejeton = nb_compo_tot-nb_bulles_reelles_;
      ArrOfInt rejeton(nb_rejeton);
      DoubleTab rejetonArea(nb_bulles_reelles_, nb_compo_tot);
      // initialisation
      rejetonArea=0.0;
      rejeton=-1 ;

      // remplissage de la matrice
      for (int fa7 = 0; fa7 < nb_facettes; fa7++)
        rejetonArea(compo_std[fa7], compo_new[fa7])+=surface_facettes[fa7];

      // On somme les contributions de chaque processeur
      mp_sum_for_each_item(rejetonArea);

      // On parcoure la matrice pour connaitre les compo_new a supprimer
      for (int ibul = 0 ; ibul<nb_bulles_reelles_ ; ibul++)
        {
          // pour chaque compo on repere celle qui va rester
          int ibul_stay = 0 ;
          for (int ibul_new = 0 ; ibul_new<nb_compo_tot ; ibul_new++)
            if (rejetonArea(ibul, ibul_new)!=0.0 && rejetonArea(ibul, ibul_stay)<rejetonArea(ibul, ibul_new))
              {
                ibul_stay=ibul_new ;
                break;
              }

          // On envoie les autres dans rejetons a supprimer
          for (int ibul_new = 0 ; ibul_new<nb_compo_tot ; ibul_new++)
            if (rejetonArea(ibul, ibul_new)!=0.0 && ibul_new!=ibul_stay)
              for (int ii = 0; ii<nb_rejeton ; ii++)
                if (rejeton[ii]==-1)
                  {
                    rejeton[ii]=ibul_new;
                    break;
                  }
        }
      Cerr << "matrice de surface des compo" << finl;
      Cerr << rejetonArea << finl;
      Cerr << "Les rejetons a supprimer sont" << finl;
      Cerr <<  rejeton << finl;

      // on remplace les compo a supprimer par -1
      for (int fa7 = 0; fa7 < nb_facettes; fa7++)
        for (int ii = 0; ii<nb_rejeton ; ii++)
          if (compo_new[fa7] ==  rejeton[ii])
            maillage.set_composante_connexe(fa7, -1);

      // On supprime la bulle que l'on vient de renumeroter -1:
      supprimer_duplicata_bulles();
    }

  //Selon ou cette methode est appelee, il faut qu'elle recree les duplicatas.
  if (duplicatas_etaient_presents)
    {
      creer_duplicata_bulles();
      maillage_ft_ijk_.parcourir_maillage();
    }

//dumplata_newtime("DNS3.lata",0.);
//dumplata_ft_mesh("DNS3.lata", "INTERFACES", maillage_ft_ijk_, 0);
//dumplata_ft_field("DNS3.lata", "INTERFACES", "COMPO_CONNEXE_AVT", "ELEM",  maillage_ft_ijk_.compo_connexe_facettes(), 0);
}

// Rempli le champ de force de rappel pour les bulles fixes.
// coef_rayon_force_rappel : coef de taille du domaine de rappel const (attention aux superposition de bulles)
void IJK_Interfaces::compute_external_forces_(IJK_Field_vector3_double& rappel_ft,
                                              IJK_Field_vector3_double& rappel,
                                              const IJK_Field_vector3_double& vitesse,
                                              const IJK_Field_double& indic/*_ns*/,
                                              const IJK_Field_double& indic_ft,
                                              const double coef_immo,
                                              const int tstep,
                                              const double current_time,
                                              const double coef_ammortissement,
                                              const double coef_rayon_force_rappel,
                                              const double integration_time,
                                              const double coef_mean_force,
                                              const double coef_force_time_n)
{
  //////////////////////////////////// CALCUL DES VITESSE MOYENNE PAR BULLES ///////////////////////////////////////
  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  //const DoubleTab& sommets = mesh.sommets() ; // Tableau des coordonnees des marqueurs.
  //int nbsom = sommets.dimension(0);
  //DoubleTab deplacement(nbsom,3);

  compute_vinterp(); // to resize and fill vinterp_

  // Les sommets virtuels sont peut-etre trop loin pour pouvoir interpoler leur vitesse,
  // il faut faire un echange espace virtuel pour avoir leur vitesse.
  mesh.desc_sommets().echange_espace_virtuel(vinterp_);

  // Calcul d'un deplacement preservant la distribution des noeuds sur les bulles:
  // Inspiree de : Transport_Interfaces_FT_Disc::calculer_vitesse_repere_local

  const int nbulles_reelles = get_nb_bulles_reelles();
  const int nbulles_ghost = get_nb_bulles_ghost();
  const int nbulles_tot = nbulles_reelles + nbulles_ghost;
  const ArrOfInt& compo_connex = mesh.compo_connexe_facettes();
  //  assert(compo_connex.size_array() == 0 || min_array(compo_connex) >=0); // Les duplicatas ne sont pas presents pendant le transport.
  // Nouveau depuis le 13/03/2014 : Les bulles ghost sont autorisees lors du transport...
  ArrOfDouble surface_par_bulle;
  calculer_surface_bulles(surface_par_bulle);
  const ArrOfDouble& surface_facette = mesh.get_update_surface_facettes();
  ArrOfIntFT compo_connex_som;
  mesh.calculer_compo_connexe_sommets(compo_connex_som);

  DoubleTab vitesses_translation_bulles(nbulles_tot,3);
  calculer_vmoy_translation_composantes_connexes(mesh,
                                                 surface_facette,
                                                 surface_par_bulle,
                                                 compo_connex,
                                                 nbulles_reelles,
                                                 nbulles_ghost,
                                                 vinterp_,
                                                 vitesses_translation_bulles);

  //////////////////////////////////// FIN CALCUL DES VITESSE MOYENNE PAR BULLES ///////////////////////////////////////

  //////////////////////////////// CALCUL DES FORCES (algo inspired by Thomas & Bolotnov) ///////////////////////////////////
  ArrOfDouble volume_reel;
  DoubleTab position;
  calculer_volume_bulles(volume_reel, position);

  const double deltat = 1.; // Dirty code... but it's actually a time integral...
  DoubleTab individual_forces(nb_bulles_reelles_,3);
  for (int idir=0; idir < 3; idir++)
    {
      const double ldom = rappel.get_domaine().get_domain_length(idir);
      for (int ib=0; ib < nb_bulles_reelles_; ib++)
        {
          if (Process::je_suis_maitre())
            {
              double dx = positions_reference_(ib,idir)-position(ib,idir);
              // sort of modulo of domain length:
              while (dx> 0.8*ldom)
                dx -=ldom;

              while (dx< -0.8*ldom)
                dx +=ldom;

              individual_forces(ib,idir) = coef_mean_force*mean_force_(ib,idir);
              individual_forces(ib,idir) += (1.-coef_mean_force)*coef_immo*(dx);
              individual_forces(ib,idir) += (1.-coef_mean_force)*coef_ammortissement*vitesses_translation_bulles(ib,idir);
              individual_forces(ib,idir) += (1.-coef_mean_force)*coef_force_time_n*force_time_n_(ib,idir);
              /*
              mean_force_(ib,idir)*=integration_time;
              mean_force_(ib,idir)+=individual_forces(ib,idir)*deltat;
              mean_force_(ib,idir)/=(integration_time+deltat);
              force_time_n_(ib,idir)=individual_forces(ib,idir);
              */
            }
        }
    }
  if (Process::je_suis_maitre())
    {
      for (int idir=0; idir < 3; idir++)
        for (int ib = 0; ib < nb_bulles_reelles_; ib++)
          {
            mean_force_(ib,idir)*=integration_time;
            mean_force_(ib,idir)+=individual_forces(ib,idir)*deltat;
            mean_force_(ib,idir)/=(integration_time+deltat);
            force_time_n_(ib,idir)=individual_forces(ib,idir);
          }
    }
  envoyer_broadcast(individual_forces, 0);
  // envoyer_broadcast(mean_force_, 0); Normally, only proc 0 needs the mean_force_

  ////////////////////////////// FIN CALCUL DES FORCES (algo inspired by Thomas & Bolotnov) ////////////////////////////////

  // Choix de la methode implementee :
  const int parser = parser_;
  if (parser)
    {
      compute_external_forces_parser(rappel, indic, individual_forces, volume_reel, position, coef_rayon_force_rappel);
    }
  else
    {
      // The table "individual_forces" enters with the value for each force and is modified to contain
      // the integral of F over the volume where the bubble force is applied. And is finally divided by that volume.
      // Then, the value (that would be applied over the whole volume) is in the table
      compute_external_forces_color_function(rappel_ft,indic/*_ns*/, indic_ft,individual_forces, volume_reel, position);
      //Cerr << "t= "<< current_time << " Max-abs(individual_forces)= "<< max_abs_array(individual_forces) << finl;
    }

  // In the current state, every time-step is post-processed.
  if (Process::je_suis_maitre())
    {
      char s[1000];
      const char *nomcas = nom_du_cas();
      SFichier fic;
      int reset = (!reprise_) && (tstep==0);
      IOS_OPEN_MODE mode = (reset) ? ios::out : ios::app;
      for (True_int idir=0; idir < 3; idir++)
        {
          snprintf(s, 1000, "%s_bulles_external_force_every_%d.out", nomcas, idir);
          // Cerr << "Ecriture des donnees par bulles: fichier " << s << finl;
          fic.ouvrir(s, mode);
          if (reset)
            {
              // Header in the file:
              snprintf(s, 1000, "# Individual forces applied inside chiv[bubble_i]=1.\n# value=1./Vol_bubble \\int F_j(bubble_i) dv");
              fic << s;
              fic << finl;
            }
          snprintf(s, 1000, "%.16e ", current_time);
          fic << s;
          for (int ib = 0; ib < nb_bulles_reelles_; ib++)
            {
              snprintf(s, 1000, "%.16e ", individual_forces(ib,idir));
              fic << s;
            }
          fic << finl;
          fic.close();
        }
    }
}

// The method is based on an eulerian color function refering to each bubble
// BEWARE : individual_forces should contain the value of each force in inlet and
//          modify it to return the integrated value (homogeneous to "F*vol") at the end of the function
void IJK_Interfaces::compute_external_forces_color_function(IJK_Field_vector3_double& rappel_ft,
                                                            const IJK_Field_double& indic_ns,
                                                            const IJK_Field_double& indic_ft,
                                                            DoubleTab& individual_forces,
                                                            const ArrOfDouble& volume_reel,
                                                            const DoubleTab& position)
{
  assert(ghost_compo_converter_.size_array() == nb_bulles_ghost_);
  // Step 1. Conversion table.
  ArrOfInt decodeur_num_compo(nb_compo_in_num_compo_);
  decodeur_num_compo=NUM_COMPO_INVALID; // invalid value is large and negative.
  // The values in num_compo are not related to the compo numbers of the bubbles.
  // Some reordering is required.

  // Carefull : domaine_vdf is built on splitting_ft_, so num_compo_ refers to these numbers
  const Domaine_IJK& s_ft =indic_ft.get_domaine();
  Int3 ijk_global, ijk_local, ijk_processeur;
  for (int ib=0; ib < nb_bulles_reelles_+nb_bulles_ghost_; ib++)
    {
      //int num_proc = -1; // invalid initialization.
      double x = position(ib,0);
      double y = position(ib,1);
      double z = position(ib,2);
      ijk_processeur[0]=0;
      ijk_processeur[1]=0;
      ijk_processeur[2]=0;
      // Le centre des bulles est forcement dans le domaine ft (etendu) sur lequel est construit num_compo_
      // Donc les indices retournes doivent etre dans les bornes.
      s_ft.search_elem(x,y,z, ijk_global, ijk_local, ijk_processeur);

      if (ijk_processeur[0] * ijk_processeur[1] * ijk_processeur[2] == 1)
        {
          // num_elem_zvdf contient le numero de l'elem dans sa Domaine_VF sur le proc num_proc.
          // Normalement, c'est cette cellule qui est au centre de la bulle ib.
          const int num_elem_zvdf = s_ft.convert_ijk_cell_to_packed(ijk_local);
          const int icompo = num_compo_[num_elem_zvdf]; // la valeur dans le tableau eulerien (Domaine i,j,k ecrit en non structure VDF)
          assert(icompo>=0);
          decodeur_num_compo[icompo] = ib; // On remplit le decodeur avec le numero de la bulle.
          //Cerr << "Decodeur : ib=" << ib << " in icompo= " << icompo << " on proc: " << Process::me() << finl;
          //Journal() << "Decodeur : ib=" << ib << " in icompo= " << icompo << " on proc: " << Process::me() << finl;
        }
    }
  // To put every process with the correct information for de-cryption
  mp_max_for_each_item(decodeur_num_compo);

  // Si j'ai bien tout saisi, il en manque une qui n'est pas dans le decodeur, c'est le numero de la phase continue:
  // Recherche de l'indice invalid dans la liste :
  // decodeur_num_compo["num compo de la phase continue"]  vaut NUM_COMPO_INVALID et donc il suffirait d'une simple boucle pour le trouver
#if 0
  int phase_continue=-1000;
  {
    int i=0;
    for (i=0; i<nb_compo_in_num_compo_; i++)
      {
        if (decodeur_num_compo[i] == NUM_COMPO_INVALID)
          {
            phase_continue = i;

            break; // on a trouve le liquid.
          }
      }
    if (i==nb_compo_in_num_compo_)
      {
        Cerr << "Continuous phase has not been found in num_compo_ " <<finl;
        Process::exit();
      }
  }
#endif
  ArrOfInt list_continuous_phase;
  list_continuous_phase.resize_array(0);
  int phase_continue1=-1000;
  int phase_continue2=-1000;
  int phase_continue3=-1000;
  int phase_continue4=-1000;
  int phase_continue5=-1000;
  int phase_continue6=-1000;
  {
    for (int i=0; i<nb_compo_in_num_compo_; i++)
      {
        if (decodeur_num_compo[i] == NUM_COMPO_INVALID)
          {
            list_continuous_phase.append_array(i);
            if (phase_continue1<0)
              {
                phase_continue1 = i;
              }
            else if (phase_continue2<0)
              {
                phase_continue2 = i;
              }
            else if (phase_continue3<0)
              {
                phase_continue3 = i;
              }
            else if (phase_continue4<0)
              {
                phase_continue4 = i;
              }
            else if (phase_continue5<0)
              {
                phase_continue5 = i;
              }
            else if (phase_continue6<0)
              {
                phase_continue6 = i;
              }
            /*            else
                          {
                            Cerr << "Too many portions of continuous phase : " << list_continuous_phase.size_array() << ". Contact TRUST support." << finl;
                            Process::exit();
                          } */
          }
      }
  }
  const int nb_parts_continuous = list_continuous_phase.size_array();
  Cerr << "nb parts continuous phase : " << list_continuous_phase.size_array()
       << " continuous phases are : " << phase_continue1 << " " << phase_continue2 << " "
       << phase_continue3 << " " << phase_continue4 << " " << phase_continue5 << " " << phase_continue6 << finl;
  Cerr << "GB-check nb_compo - nb_continuous phase - nb_bulles_tot = "
       << nb_compo_in_num_compo_<< " - "  << nb_parts_continuous << " - " << nb_bulles_reelles_+nb_bulles_ghost_
       << " = " << nb_compo_in_num_compo_ - (nb_parts_continuous +nb_bulles_reelles_+nb_bulles_ghost_)<< finl;
  assert(nb_compo_in_num_compo_ - (nb_parts_continuous +nb_bulles_reelles_+nb_bulles_ghost_)==0);
  /*if (list_continuous_phase.size_array() >6)
  {
   Cerr << "Too many portions of continuous phase : " << list_continuous_phase.size_array() << ". Contact TRUST support." << finl;
   Process::exit();
  } else */if (nb_parts_continuous == 0)
    {
      Cerr << "No continuous phase found. " << finl;
      Process::exit();
    }
  // Step 2. On all processors.
  // Par pure commodite, on travaille sur les champs aux elements num_compo_ et indic pour decider si
  // si on applique la force de rappel sur la face (i,j,k) "de gauche". C'est donc legerement biaise.
  // Mais on n'a pas besoin d'une immense rigueur non-plus... donc on s'en accomode.
  //
  // Cette partie doit etre faite sur le domaine FT pour pouvoir interroger num_compo_[num_elem_zvdf]
  // qui doit etre construit sur le FT!

  DoubleTab integrated_forces(nb_bulles_reelles_/* For real bubbles only*/, 3/*dims*/);
  const double vol = s_ft.get_constant_delta(0)
                     * s_ft.get_constant_delta(1)
                     * s_ft.get_constant_delta(2);
  IntTab integration_cells_per_bubble(nb_bulles_reelles_/* For real bubbles only*/, 3/*dims because can be different due to MAC scheme */);
  integration_cells_per_bubble=0;
  integrated_forces =0.;
  for (int idir=0; idir < 3; idir++)
    {
      const int nx = rappel_ft[idir].ni();
      const int ny = rappel_ft[idir].nj();
      int nzdeb = 0;
      int nz = rappel_ft[idir].nk();
      if ((!s_ft.get_periodic_flag(DIRECTION_K)) &&
          (idir==DIRECTION_K))
        {
          force_zero_on_walls(rappel_ft[DIRECTION_K]);
          const int kmin = s_ft.get_offset_local(DIRECTION_K);
          const int nktot = s_ft.get_nb_items_global(Domaine_IJK::FACES_K, DIRECTION_K);
          if (kmin + nz == nktot)
            {
              // On the "last" proc (in dir==2), there is a last layer of faces in this direction. We cannot do convert_ijk_cell_to_packed
              // for them, so we should skip them (and set them to 0.
              nz-=1;
            }

          if (kmin == 0)
            {
              // On the "first" (in dir==2), we should skip the wall to the left...
              nzdeb = 1;
            }
        }
      for (int k=nzdeb; k < nz; k++)
        for (int j=0; j< ny; j++)
          for (int i=0; i < nx; i++)
            {
              // Probleme : s_ft est etendu, il faut fonc etre sur rappel_ft pour les indices (i,j,k).
              const int num_elem_zvdf = s_ft.convert_ijk_cell_to_packed(i,j,k);
              const int icompo = num_compo_[num_elem_zvdf]; // Ce ne sont pas les numeros qu'a le FT sur lui.
              bool continuous = false;
              for (int c=0; c<nb_parts_continuous; c++)
                {
                  if (icompo == list_continuous_phase[c])
                    {
                      continuous = true;
                      break;
                    }
                }
              if ((icompo==-1) || (continuous) )
                {
                  rappel_ft[idir](i,j,k) = 0.; // pas de force dans la phase continue, ni dans les mailles interfaciales.
                }
              else
                {
                  int id_bulle = decodeur_num_compo[icompo]; // On relit le numero de la bulle grace au decodeur.
                  // Si on obtient une bulle ghost, on cherche l'indice reel correspondant
                  // afin de trouver la valeur de la force que l'on doit y imposer
                  if (id_bulle>=nb_bulles_reelles_)
                    {
                      // Dealing with a ghost bubble:
                      const int ighost = ghost_compo_converter(id_bulle-nb_bulles_reelles_);
                      const int ibulle_reelle = decoder_numero_bulle(-ighost);
                      //Cerr <<"id_bulle vs ibulle_reelle : " <<  id_bulle << " " <<ibulle_reelle << finl;
                      //Cerr << "decodeur_num_compo : " << decodeur_num_compo << finl;
                      id_bulle = ibulle_reelle;
                    }
                  else
                    {
                      // There should not be negative indices here.
                      assert(id_bulle>=0);
                      // For real bubbles, we compute the integral over the extended domain:
                      const double f = individual_forces(id_bulle,idir);
                      integrated_forces(id_bulle,idir) += vol*f;
                      integration_cells_per_bubble(id_bulle,idir) +=1;
                    }
                  // We have taken necessary precautions for id_bulle to be the index of the real bubble:
                  assert((id_bulle<nb_bulles_reelles_) &&(id_bulle>=0));

                  // For all bubbles (real or ghost), we apply the force to their volume
                  //    (it will truely be effective in NS only, were velocity is really solved!)
                  const double f = individual_forces(id_bulle,idir);
                  rappel_ft[idir](i,j,k) = f;
                }
            }
    }
  mp_sum_for_each_item(integrated_forces);
  mp_sum_for_each_item(integration_cells_per_bubble);
  // Individual_forces has been used. It can be updated with the true value (weighted by the true volume where the force is applied).
// individual_forces=integrated_forces;
//  individual_forces.copy(integrated_forces, RESIZE_OPTIONS::COPY_INIT); // Return the integrated value in individual_forces.
//
  for (int idir=0; idir < 3; idir++)
    for (int ib = 0; ib < nb_bulles_reelles_; ib++)
      {
        // individual_forces(ib, idir) = integrated_forces(ib, idir) / volume_reel(ib); // WRONG : the force is not applied to the rigorous bubble volume..
        individual_forces(ib, idir) = integrated_forces(ib, idir) / (vol*integration_cells_per_bubble(ib,idir));
      }
}

// individual_forces : The instantaneous value of the force for each bubble.
// force_time_n_ : The same, but stored in the class for later use (only on master process).
// mean_force_   : Time average of the force for each bubble (only on master process).
void IJK_Interfaces::compute_external_forces_parser(IJK_Field_vector3_double& rappel,
                                                    const IJK_Field_double& indic, // ns
                                                    const DoubleTab& individual_forces,
                                                    const ArrOfDouble& volume_reel,
                                                    const DoubleTab& position,
                                                    const double coef_rayon_force_rappel)
{
  Noms noms_forces; // on attend trois expressions
  noms_forces.dimensionner_force(3);
  for (int idir=0; idir < 3; idir++)
    {
      noms_forces[idir] = "";
      for (int ib=0; ib < nb_bulles_reelles_; ib++)
        {
          double xir = position(ib,0);
          double yir = position(ib,1);
          double zir = position(ib,2);
          double r2 = coef_rayon_force_rappel*coef_rayon_force_rappel*pow((volume_reel[ib]*3)/(4.*M_PI), 2./3.);
          noms_forces[idir] += "+((";
          noms_forces[idir] += Nom("(X-(")+Nom(xir, "%g")+Nom("))*(X-(")+Nom(xir, "%g")+Nom("))");
          noms_forces[idir] += Nom("+(Y-(")+Nom(yir, "%g")+Nom("))*(Y-(")+Nom(yir, "%g")+Nom("))");
          noms_forces[idir] += Nom("+(Z-(")+Nom(zir, "%g")+Nom("))*(Z-(")+Nom(zir, "%g")+Nom("))");
          noms_forces[idir] += Nom("-")+Nom(r2, "%g");
          // 1-ff is the code name for chi_v in our case :
          // on veut > 0.9999 pour ne pas changer la gravite dans les mailles ou il y a une interface, ce qui rend la methode instable
          // d'apres Thomas & Bolotnov.
          noms_forces[idir] += ")_LT_0.)*(1.-ff)*((1.-ff)_GT_0.000001)*(" ;
          noms_forces[idir] += Nom(individual_forces(ib,idir), "%g");
          noms_forces[idir] += ")" ;
        }
    }

  for (int idir=0; idir < 3; idir++)
    {
      // Pour les bulles ghosts
      //Cerr << "Ghost_compo_converter : " << ghost_compo_converter_ << finl;
      assert(ghost_compo_converter_.size_array() == nb_bulles_ghost_);
      for (int ibg=0; ibg < nb_bulles_ghost_; ibg++)
        {
          const int ighost = ghost_compo_converter(ibg);
          const int ibulle_reelle = decoder_numero_bulle(-ighost);
          //Cerr << ighost << " " << ibulle_reelle << finl;

          double xir = position(nb_bulles_reelles_+ibg,0);
          double yir = position(nb_bulles_reelles_+ibg,1);
          double zir = position(nb_bulles_reelles_+ibg,2);
          double r2 = coef_rayon_force_rappel*coef_rayon_force_rappel*pow((volume_reel[ibg]*3)/(4.*M_PI), 2./3.);
          noms_forces[idir] += "+((";
          noms_forces[idir] += Nom("(X-(")+Nom(xir, "%g")+Nom("))*(X-(")+Nom(xir, "%g")+Nom("))");
          noms_forces[idir] += Nom("+(Y-(")+Nom(yir, "%g")+Nom("))*(Y-(")+Nom(yir, "%g")+Nom("))");
          noms_forces[idir] += Nom("+(Z-(")+Nom(zir, "%g")+Nom("))*(Z-(")+Nom(zir, "%g")+Nom("))");
          noms_forces[idir] += Nom("-")+Nom(r2, "%g");
          // 1-ff is the code name for chi_v in our case :
          // on veut > 0.9999 pour ne pas changer la gravite dans les mailles ou il y a une interface, ce qui rend la methode instable
          noms_forces[idir] += ")_LT_0.)*(1.-ff)*((1.-ff)_GT_0.000001)*(" ;
          noms_forces[idir] += Nom(individual_forces(ibulle_reelle,idir), "%g");
          noms_forces[idir] += ")" ;
        }
      // Cerr << "Setting force [dir=" << idir << "] : " << noms_forces[idir] << finl;
      // Cette methode parcours ni(), nj() et nk() et donc pas les ghost...
      set_field_data(rappel[idir], noms_forces[idir], indic, 0. /* could be time... */ );
    }
}

void IJK_Interfaces::compute_indicatrice_non_perturbe(IJK_Field_double& indic_np,
                                                      const IJK_Field_double& indic,
                                                      const ArrOfDouble& volume_reel,
                                                      const DoubleTab& position) const
{
  Nom nom_indicatrices_np;
  nom_indicatrices_np = "";
  for (int ib=0; ib < nb_bulles_reelles_; ib++)
    {
      double xir = positions_reference_(ib,0);
      double yir = positions_reference_(ib,1);
      double zir = positions_reference_(ib,2);
      // Une ellipse autour de la bulle :
      double r=pow((volume_reel[ib]*3)/(4.*M_PI), 1./3.);
      double a=5.0*r;
      double b=2.5*r;
      double c=b;
      double x0=xir-2*r;
      double y0=yir;
      double z0=zir;

      nom_indicatrices_np += "+((";
      nom_indicatrices_np += Nom("((X-(")+Nom(x0, "%g")+Nom("))*(X-(")+Nom(x0, "%g")+Nom("))/(")+Nom(a, "%g")+Nom("*")+Nom(a, "%g")+Nom("))");
      nom_indicatrices_np += Nom("+((Y-(")+Nom(y0, "%g")+Nom("))*(Y-(")+Nom(y0, "%g")+Nom("))/(")+Nom(b, "%g")+Nom("*")+Nom(b, "%g")+Nom("))");
      nom_indicatrices_np += Nom("+((Z-(")+Nom(z0, "%g")+Nom("))*(Z-(")+Nom(z0, "%g")+Nom("))/(")+Nom(c, "%g")+Nom("*")+Nom(c, "%g")+Nom("))");
      nom_indicatrices_np += Nom("-1.0");
      // 1-ff is the code name for chi_v in our case :
      nom_indicatrices_np += ")_LT_0.)*(1.0)" ;
    }
  // Pour les bulles ghosts
  //Cerr << "Ghost_compo_converter : " << ghost_compo_converter_ << finl;

  assert(ghost_compo_converter_.size_array() == nb_bulles_ghost_);
  for (int ibg=0; ibg < nb_bulles_ghost_; ibg++)
    {
      //const int ighost = ghost_compo_converter(ibg);
      //const int ibulle_reelle = decoder_numero_bulle(-ighost);
      //Cerr << ighost << " " << ibulle_reelle << finl;

      double xir = position(nb_bulles_reelles_+ibg,0);
      double yir = position(nb_bulles_reelles_+ibg,1);
      double zir = position(nb_bulles_reelles_+ibg,2);
      double r=pow((volume_reel[ibg]*3)/(4.*M_PI), 1./3.);
      double a=5*r;
      double b=2.5*r;
      double c=b;
      double x0=xir-2*r;
      double y0=yir;
      double z0=zir;

      nom_indicatrices_np += "+((";
      nom_indicatrices_np += Nom("((X-(")+Nom(x0, "%g")+Nom("))*(X-(")+Nom(x0, "%g")+Nom("))/(")+Nom(a, "%g")+Nom("*")+Nom(a, "%g")+Nom("))");
      nom_indicatrices_np += Nom("+((Y-(")+Nom(y0, "%g")+Nom("))*(Y-(")+Nom(y0, "%g")+Nom("))/(")+Nom(b, "%g")+Nom("*")+Nom(b, "%g")+Nom("))");
      nom_indicatrices_np += Nom("+((Z-(")+Nom(z0, "%g")+Nom("))*(Z-(")+Nom(z0, "%g")+Nom("))/(")+Nom(c, "%g")+Nom("*")+Nom(c, "%g")+Nom("))");
      nom_indicatrices_np += Nom("-1.0");
      // 1-ff is the code name for chi_v in our case :
      nom_indicatrices_np += ")_LT_0.)*(1.0)" ;

    }
  Cerr << "Setting indicatrice_non_perturbe :" << nom_indicatrices_np << finl;
  // Cette methode parcours ni(), nj() et nk() et donc pas les ghost...
  set_field_data(indic_np, nom_indicatrices_np, indic, 0. /* could be time... */ );
}

// Dans cette methode on calcule l'indicatrice_next_ en fonction de
// la variable interfaces_.
// /!\ Si l'interface n'a pas ete deplacee on recalcule la mm chose.
// A n'appeler qu'apres un deplacement d'interface donc ou pour
// initialisation.
void IJK_Interfaces::calculer_indicatrice_next(
  const DoubleTab& gravite,
  const double delta_rho,
  const double sigma,
  const double time,
  const int itstep,
  const bool parcourir
)
{
  // En monophasique, les champs sont a jours donc on zap :
  if (Option_IJK::DISABLE_DIPHASIQUE)
    return;

  Navier_Stokes_FTD_IJK& ns = ref_ijk_ft_->eq_ns();

  static Stat_Counter_Id calculer_indicatrice_next_counter_ = statistiques().new_counter(2, "Calcul Indicatrice Next");
  statistiques().begin_count(calculer_indicatrice_next_counter_);
  // En diphasique sans bulle (pour cas tests), on met tout a 1.
  if (get_nb_bulles_reelles() == 0)
    {
      indicatrice_ft_[next()].data() = 1.;
      indicatrice_ns_[next()].data() = 1.;
      indicatrice_ft_[next()].echange_espace_virtuel(indicatrice_ft_[next()].ghost());
      indicatrice_ns_[next()].echange_espace_virtuel(indicatrice_ns_[next()].ghost());

      if (parcourir)
        parcourir_maillage();
      return;
    }

  if (parcourir)
    parcourir_maillage();

  // Calcul de l'indicatrice sur le domaine etendu :
  // faut-il calculer les valeurs de l'indicatrice from scratch ? (avec methode
  // des composantes connexes)
  if (get_recompute_indicator())
    calculer_indicatrice(indicatrice_ft_[next()]);
  else
    calculer_indicatrice_optim(indicatrice_ft_[next()]);

  indicatrice_ft_[next()].echange_espace_virtuel(indicatrice_ft_[next()].ghost());

  // Calcul de l'indicatrice sur le domaine NS :
  ns.redistrib_from_ft_elem().redistribute(
    indicatrice_ft_[next()], indicatrice_ns_[next()]);
  indicatrice_ns_[next()].echange_espace_virtuel(indicatrice_ns_[next()].ghost());

  // Calcul des indicatrices s'il y a des groupes :
  const int nb_grps = IJK_Interfaces::nb_groups();
  if (nb_grps > 1)
    {
      if (get_recompute_indicator())
        calculer_indicatrices(groups_indicatrice_ft_[next()]);
      else
        calculer_indicatrices_optim(groups_indicatrice_ft_[next()]);

      groups_indicatrice_ft_[next()].echange_espace_virtuel();

      // Calcul de l'indicatrice sur le domaine NS :
      ns.redistrib_from_ft_elem().redistribute(groups_indicatrice_ft_[next()], groups_indicatrice_ns_[next()]);
      groups_indicatrice_ns_[next()].echange_espace_virtuel();
    }

  // Aux cellules diphasiques, calcule toutes les moyennes de l'interface
  // dans les cellules pour chaque compo. Le but est de le faire une fois
  // pour toute de maniere synchronisee (et pas au moment ou on calcule la
  // force par exemple).

  val_par_compo_in_cell_computation_.calculer_valeur_par_compo(
    time, itstep,
    nb_compo_traversante_[next()],
    compos_traversantes_[next()],
    normale_par_compo_[next()],
    bary_par_compo_[next()],
    indicatrice_par_compo_[next()],
    surface_par_compo_[next()],
    courbure_par_compo_[next()]);

  // calcul de la force de repulsion
  calculer_phi_repuls_par_compo(
    surf_par_compo_[next()],
    source_interf_par_compo_[next()],
    phi_par_compo_[next()],
    repuls_par_compo_[next()],
    gravite,
    delta_rho,
    sigma,
    time,
    itstep
  );

#if VERIF_INDIC
  verif_indic();
#endif

  // Calcul normale_of_interf_ bary_of_interf_ et passage a NS
  mean_over_compo(normale_par_compo_[next()], nb_compo_traversante_[next()], normal_of_interf_[next()]);
  mean_over_compo(bary_par_compo_[next()], nb_compo_traversante_[next()], bary_of_interf_[next()]);

  for (int c=0; c < 3; c++)
    {
      ns.redistrib_from_ft_elem().redistribute(
        normal_of_interf_[next()][c],
        normal_of_interf_ns_[next()][c]);
      ns.redistrib_from_ft_elem().redistribute(
        bary_of_interf_[next()][c],
        bary_of_interf_ns_[next()][c]);
    }
  normal_of_interf_ns_[next()].echange_espace_virtuel();
  bary_of_interf_ns_[next()].echange_espace_virtuel();

  // Aux faces mouillees
  n_faces_mouilles_[next()] = surface_vapeur_par_face_computation_.compute_surf_and_barys(
                                maillage_ft_ijk_,
                                indicatrice_ft_[next()],
                                normal_of_interf_[next()],
                                surface_vapeur_par_face_[next()],
                                barycentre_vapeur_par_face_[next()]);
  // Passage au domaine NS
  ns.get_redistribute_from_splitting_ft_faces(
    surface_vapeur_par_face_[next()],
    surface_vapeur_par_face_ns_[next()]);
  for (int c=0; c<3; c++)
    ns.get_redistribute_from_splitting_ft_faces(
      barycentre_vapeur_par_face_[next()][c],
      barycentre_vapeur_par_face_ns_[next()][c]);

  // Overwriting the MedCoupling computation of the face surfaces.
  //
  //{
  //  const Domaine_IJK& s = indicatrice_ns_[next()].get_domaine();
  //
  //  double dx = s.get_constant_delta(DIRECTION_I);
  //  double dy = s.get_constant_delta(DIRECTION_J);
  //  double dz = s.get_constant_delta(DIRECTION_K);
  //
  //  const Domaine_IJK& geom = s.get_grid_geometry();
  //  double origin_x = geom.get_origin(DIRECTION_I);
  //  double origin_y = geom.get_origin(DIRECTION_J);
  //  double origin_z = geom.get_origin(DIRECTION_K);
  //
  //  int offset_x = s.get_offset_local(DIRECTION_I);
  //  int offset_y = s.get_offset_local(DIRECTION_J);
  //  int offset_z = s.get_offset_local(DIRECTION_K);
  //
  //  const int ni = indicatrice_ns_[next()].ni();
  //  const int nj = indicatrice_ns_[next()].nj();
  //  const int nk = indicatrice_ns_[next()].nk();
  //  for (int k = 0; k < nk; k++)
  //    {
  //      for (int j = 0; j < nj; j++)
  //        {
  //          for (int i = 0; i < ni; i++)
  //            {
  //              int phase = 0;
  //
  //              for (int dir = 0; dir < 3; dir++)
  //                {
  //                  double face_bary_x = origin_x + dx*(i + offset_x + (get_barycentre_face(1, dir, 0, phase, i,j,k)));
  //                  double face_bary_y = origin_y + dy*(j + offset_y + (get_barycentre_face(1, dir, 1, phase, i,j,k)));
  //                  double face_bary_z = origin_z + dz*(k + offset_z + (get_barycentre_face(1, dir, 2, phase, i,j,k)));
  //
  //
  //                  double f_dir = select_dir(dir, dy*dz, dx*dz, dx*dy);
  //                  surface_vapeur_par_face_ns_[next()][dir](i,j,k) = (1 - next_indicatrice_surfacique)*f_dir;
  //
  //                  barycentre_vapeur_par_face_ns_[next()][dir][0](i,j,k) = face_bary_x;
  //                  barycentre_vapeur_par_face_ns_[next()][dir][1](i,j,k) = face_bary_y;
  //                  barycentre_vapeur_par_face_ns_[next()][dir][2](i,j,k) = face_bary_z;
  //                }
  //            }
  //        }
  //    }
  //}
  //
  //surface_vapeur_par_face_ns_[old()].echange_espace_virtuel();
  //for (int c = 0; c < 3; c++)
  //  {
  //    barycentre_vapeur_par_face_ns_[old()][c].echange_espace_virtuel();
  //  }
  //
  //// Passage au domaine FT
  //ref_ijk_ft_->redistrib_to_ft_elem().redistribute(surface_vapeur_par_face_ns_[next()], surface_vapeur_par_face_[next()]);
  //for (int c=0; c<3; c++)
  //  {
  //    ref_ijk_ft_->redistrib_to_ft_elem().redistribute(barycentre_vapeur_par_face_ns_[next()][c], barycentre_vapeur_par_face_[next()][c]);
  //  }
  //
  //surface_vapeur_par_face_[old()].echange_espace_virtuel();
  //for (int c = 0; c < 3; c++)
  //  {
  //    barycentre_vapeur_par_face_[old()][c].echange_espace_virtuel();
  //  }

  // Calcul de la surface interfaciale
  calculer_surface_interface(surface_interface_ft_[next()], indicatrice_ft_[next()]);
  surface_interface_ft_[next()].echange_espace_virtuel(surface_interface_ft_[next()].ghost());

  // Calcul du barycentre de la phase
  calculer_barycentre(barycentre_phase1_ft_[next()], indicatrice_ft_[next()]);
  for (int bary_compo = 0; bary_compo < 3; bary_compo++)
    {
      barycentre_phase1_ft_[next()][bary_compo].echange_espace_virtuel(barycentre_phase1_ft_[next()][bary_compo].ghost());
    }

  // Calcul de l'indicatrice surfacique et du barycentre correspondant
  calculer_indicatrice_surfacique_barycentre_face(indicatrice_surfacique_face_ft_[next()], barycentre_phase1_face_ft_[next()], indicatrice_ft_[next()], normal_of_interf_[next()]);
  indicatrice_surfacique_face_ft_[next()].echange_espace_virtuel();

  for (int face_dir = 0; face_dir < 3; face_dir++)
    {
      for (int bary_compo = 0; bary_compo < 2; bary_compo++)
        {
          barycentre_phase1_face_ft_[next()][face_dir][bary_compo].echange_espace_virtuel(barycentre_phase1_face_ft_[next()][face_dir][bary_compo].ghost());
        }
    }

  // Passage au domaine NS
  ns.redistrib_from_ft_elem().redistribute(surface_interface_ft_[next()], surface_interface_ns_[next()]);
  surface_interface_ns_[next()].echange_espace_virtuel(surface_interface_ns_[next()].ghost());

  for (int d = 0; d < 3; d++)
    {
      ns.redistrib_from_ft_elem().redistribute(barycentre_phase1_ft_[next()][d], barycentre_phase1_ns_[next()][d]);
      barycentre_phase1_ns_[next()][d].echange_espace_virtuel(barycentre_phase1_ns_[next()][d].ghost());
    }

  ns.get_redistribute_from_splitting_ft_faces(
    indicatrice_surfacique_face_ft_[next()],
    indicatrice_surfacique_face_ns_[next()]);
  indicatrice_surfacique_face_ns_[next()].echange_espace_virtuel();

  for (int face_dir = 0; face_dir < 3; face_dir++)
    {
      for (int bary_compo = 0; bary_compo < 2; bary_compo++)
        {
          ns.redistrib_from_ft_elem().redistribute(barycentre_phase1_face_ft_[next()][face_dir][bary_compo], barycentre_phase1_face_ns_[next()][face_dir][bary_compo]);
          barycentre_phase1_face_ns_[next()][face_dir][bary_compo].echange_espace_virtuel(barycentre_phase1_face_ns_[next()][face_dir][bary_compo].ghost());
        }
    }

  statistiques().end_count(calculer_indicatrice_next_counter_);
}

// Dans cette methode on calcule l'indicatrice_intermediaire_ en fonction de
// la variable interfaces_intermediaire_, qui doit correspondre a l'interface
// deplacee mais avant l'etape de remaillage.
void IJK_Interfaces::calculer_indicatrice_intermediaire(
  IJK_Field_double& indicatrice_intermediaire_ft,
  IJK_Field_double& indicatrice_intermediaire_ns,
  IJK_Field_vector3_double& indicatrice_surfacique_intermediaire_face_ft,
  IJK_Field_vector3_double& indicatrice_surfacique_intermediaire_face_ns,
  const bool parcourir
)
{
  // En monophasique, les champs sont a jours donc on zap :
  if (Option_IJK::DISABLE_DIPHASIQUE)
    return;

  static Stat_Counter_Id calculer_indicatrice_next_counter_ = statistiques().new_counter(2, "Calcul Indicatrice Next");
  statistiques().begin_count(calculer_indicatrice_next_counter_);
  // En diphasique sans bulle (pour cas tests), on met tout a 1.
  if (get_nb_bulles_reelles() == 0)
    {
      indicatrice_intermediaire_ft.data() = 1.;
      indicatrice_intermediaire_ns.data() = 1.;
      indicatrice_intermediaire_ft.echange_espace_virtuel(indicatrice_intermediaire_ft.ghost());
      indicatrice_intermediaire_ns.echange_espace_virtuel(indicatrice_intermediaire_ns.ghost());

      if (parcourir)
        parcourir_maillage();
      return;
    }

  if (parcourir)
    parcourir_maillage();

  // Calcul de l'indicatrice sur le domaine etendu :
  // faut-il calculer les valeurs de l'indicatrice from scratch ? (avec methode
  // des composantes connexes)
  if (get_recompute_indicator())
    calculer_indicatrice(indicatrice_intermediaire_ft);
  else
    calculer_indicatrice_optim(indicatrice_intermediaire_ft);

  indicatrice_intermediaire_ft.echange_espace_virtuel(indicatrice_intermediaire_ft.ghost());

  // Calcul de l'indicatrice sur le domaine NS :
  ref_ijk_ft_->eq_ns().redistrib_from_ft_elem().redistribute(
    indicatrice_intermediaire_ft, indicatrice_intermediaire_ns);
  indicatrice_intermediaire_ns.echange_espace_virtuel(indicatrice_intermediaire_ns.ghost());

  // Calcul de l'indicatrice surfacique correspondante
  calculer_indicatrice_surfacique_face(indicatrice_surfacique_intermediaire_face_ft, indicatrice_intermediaire_ft, normal_of_interf_[next()]);
  indicatrice_surfacique_intermediaire_face_ft.echange_espace_virtuel();

  ref_ijk_ft_->eq_ns().get_redistribute_from_splitting_ft_faces(
    indicatrice_surfacique_intermediaire_face_ft,
    indicatrice_surfacique_intermediaire_face_ns);
  indicatrice_surfacique_intermediaire_face_ns.echange_espace_virtuel();

  statistiques().end_count(calculer_indicatrice_next_counter_);
}

#if VERIF_INDIC
void IJK_Interfaces::verif_indic()
{
  calculer_indicatrice(indicatrice_ft_test_);
  indicatrice_ft_test_.echange_espace_virtuel(indicatrice_ft_test_.ghost());
  SChaine indic;
  if (nb_grps > 1)
    {
      calculer_indicatrices(groups_indicatrice_ft_test_);
      groups_indicatrice_ft_test_.echange_espace_virtuel();
    }
  SChaine group_indic;

  const int ni = indicatrice_ft_[next()].ni();
  const int nj = indicatrice_ft_[next()].nj();
  const int nk = indicatrice_ft_[next()].nk();

  int egalite = 1;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          if (indicatrice_ft_[next()](i, j, k) != indicatrice_ft_test_(i, j, k))
            {
              egalite = 0;
              indic << "(" << i << "," << j << "," << k << ") : " << indicatrice_ft_[next()](i, j, k) << " VS "
                    << indicatrice_ft_test_(i, j, k) << finl;
            }
          for (int igroup = 0; igroup < nb_grps; igroup++)
            {
              if (groups_indicatrice_ft_[next()][igroup](i, j, k) != groups_indicatrice_ft_test_[igroup](i, j, k))
                {
                  egalite = 0;
                  group_indic << "groupe = " << igroup << " (" << i << "," << j << "," << k
                              << ") : " << groups_indicatrice_ft_[next()][igroup](i, j, k) << " VS "
                              << groups_indicatrice_ft_next_test_[igroup](i, j, k) << finl;
                }
            }
        }

  Process::Journal() << indic.get_str() << group_indic.get_str() << finl;
  if (!egalite)
    {
      Cerr << "Probleme_FTD_IJK_base:: calcul de l'indicatrice faux ! (iteration " << tstep_ << " sur " << nb_timesteps_ << " )"
           << finl;
      Process::exit();
    }
}
#endif

double IJK_Interfaces::get_barycentre(bool next_time, int bary_compo, int phase, int i, int j, int k) const
{
  int current_time = next_time ? next() : old();
  int other_time = next_time ? old() : next();

  double current_indicatrice = next_time ? In(i,j,k) : I(i,j,k);
  double other_time_indicatrice = next_time ? I(i,j,k) : In(i,j,k);

  if ((I(i,j,k) == 0.) && (In(i,j,k) == 0.))
    {
      return .5;
    }
  else if ((I(i,j,k) == 1.) && (In(i,j,k) == 1.))
    {
      return .5;
    }
  else if (current_indicatrice == 0.)
    {
      assert(barycentre_phase1_ns_[current_time][bary_compo](i,j,k) == .5);
      if (phase == 0)
        {
          return .5;
        }
      else
        {
          double other_time_bary_compo = barycentre_phase1_ns_[other_time][bary_compo](i,j,k);
          if (other_time_bary_compo > .5)
            {
              return 1;
            }
          else
            {
              return 0;
            }
        }
    }
  else if (current_indicatrice == 1.)
    {
      assert(barycentre_phase1_ns_[current_time][bary_compo](i,j,k) == .5);
      if (phase == 1)
        {
          return .5;
        }
      else
        {
          double other_time_bary_compo = barycentre_phase1_ns_[other_time][bary_compo](i,j,k);
          other_time_bary_compo = opposing_barycentre(other_time_bary_compo, other_time_indicatrice);
          if (other_time_bary_compo > .5)
            {
              return 1;
            }
          else
            {
              return 0;
            }
        }
    }
  else
    {
      double bary = barycentre_phase1_ns_[current_time][bary_compo](i,j,k);

      if (phase == 0)
        {
          bary = opposing_barycentre(bary, current_indicatrice);
        }
      assert(bary >= 0 && bary <= 1.);
      return bary;
    }
}
double IJK_Interfaces::get_barycentre_face(bool next_time, int face_dir, int bary_compo, int phase, int i, int j, int k) const
{
  int current_time = next_time ? next() : old();
  int other_time = next_time ? old() : next();

  double old_indicatrice_surfacique  = get_indicatrice_surfacique_face_old()[face_dir](i,j,k);
  double next_indicatrice_surfacique = get_indicatrice_surfacique_face_next()[face_dir](i,j,k);
  double current_indicatrice_surfacique = next_time ? next_indicatrice_surfacique : old_indicatrice_surfacique;
  double other_time_indicatrice_surfacique = next_time ? old_indicatrice_surfacique : next_indicatrice_surfacique;

  if (face_dir == bary_compo)
    {
      return 0.;
    }
  else
    {
      // Pour la face x      dir1=z -> (2->0)   dir2=y -> (1->1)
      // Pour la face y      dir1=x -> (0->0)   dir2=z -> (2->1)
      // Pour la face z      dir1=y -> (1->0)   dir2=x -> (0->1)
      int compo2D = (face_dir == 0) ? ((bary_compo == 2) ? 0 : ((bary_compo == 1) ? 1 : -1)) :
                      ((face_dir == 1) ? ((bary_compo == 0) ? 0 : ((bary_compo == 2) ? 1 : -1)) :
                       ((face_dir == 2) ? ((bary_compo == 1) ? 0 : ((bary_compo == 0) ? 1 : -1)) :
                        -1));
      assert(compo2D >= 0);
      if ((old_indicatrice_surfacique == 0.) && (next_indicatrice_surfacique == 0.))
      {
          return .5;
        }
      else if ((old_indicatrice_surfacique == 1.) && (next_indicatrice_surfacique == 1.))
        {
          return .5;
        }
      else if (current_indicatrice_surfacique == 0.)
        {
          //assert(barycentre_phase1_face_ns_[current_time][face_dir][compo2D](i,j,k) == .5);
          if (phase == 0)
            {
              return .5;
            }
          else
            {
              double other_time_bary_compo = barycentre_phase1_face_ns_[other_time][face_dir][compo2D](i,j,k);
              if (other_time_bary_compo > .5)
                {
                  return 1;
                }
              else
                {
                  return 0;
                }
            }
        }
      else if (current_indicatrice_surfacique == 1.)
        {
          //assert(barycentre_phase1_face_ns_[current_time][face_dir][compo2D](i,j,k) == .5);
          if (phase == 1)
            {
              return .5;
            }
          else
            {
              double other_time_bary_compo = barycentre_phase1_face_ns_[other_time][face_dir][compo2D](i,j,k);
              other_time_bary_compo = opposing_barycentre(other_time_bary_compo, other_time_indicatrice_surfacique);
              if (other_time_bary_compo > .5)
                {
                  return 1;
                }
              else
                {
                  return 0;
                }
            }
        }
      else
        {
          double bary = barycentre_phase1_face_ns_[current_time][face_dir][compo2D](i,j,k);

          if (phase == 0)
            {
              bary = opposing_barycentre(bary, current_indicatrice_surfacique);
            }
          assert(bary >= 0 && bary <= 1.);
          return bary;
        }
    }
}

// Copie de l'interface et des intersections associees au pas de temps precedent
void IJK_Interfaces::update_old_intersections()
{
  old_maillage_ft_ijk_.recopie(maillage_ft_ijk_, Maillage_FT_Disc::COMPLET);
}

void IJK_Interfaces::switch_indicatrice_next_old()
{
  if (Option_IJK::DISABLE_DIPHASIQUE)
    return;

  old_en_premier_ = not old_en_premier_;

  champs_compris_.switch_ft_fields();

  // TODO: verifier la liste des echanges espace virtuels
  // TODO: il faut choisir, soit je les fait sur les next soit sur les old, mais
  // pas les deux.
  indicatrice_ft_[old()].echange_espace_virtuel(indicatrice_ft_[old()].ghost());
  groups_indicatrice_ft_[old()].echange_espace_virtuel();
  nb_compo_traversante_[old()].echange_espace_virtuel(nb_compo_traversante_[old()].ghost());
  normale_par_compo_[old()].echange_espace_virtuel();
  bary_par_compo_[old()].echange_espace_virtuel();
  surface_par_compo_[old()].echange_espace_virtuel();

  surface_interface_ft_[old()].echange_espace_virtuel(surface_interface_ft_[old()].ghost());
  surface_interface_ns_[old()].echange_espace_virtuel(surface_interface_ns_[old()].ghost());

  for (int bary_compo = 0; bary_compo < 3; bary_compo++)
    {
      barycentre_phase1_ft_[old()][bary_compo].echange_espace_virtuel(barycentre_phase1_ft_[old()][bary_compo].ghost());
      barycentre_phase1_ns_[old()][bary_compo].echange_espace_virtuel(barycentre_phase1_ns_[old()][bary_compo].ghost());
    }

  indicatrice_surfacique_face_ft_[old()].echange_espace_virtuel();
  indicatrice_surfacique_face_ns_[old()].echange_espace_virtuel();

  for (int face_dir = 0; face_dir < 3; face_dir++)
    {
      for (int bary_compo = 0; bary_compo < 2; bary_compo++)
        {
          barycentre_phase1_face_ft_[old()][face_dir][bary_compo].echange_espace_virtuel(barycentre_phase1_face_ft_[old()][face_dir][bary_compo].ghost());
          barycentre_phase1_face_ns_[old()][face_dir][bary_compo].echange_espace_virtuel(barycentre_phase1_face_ns_[old()][face_dir][bary_compo].ghost());
        }
    }

  normal_of_interf_[old()].echange_espace_virtuel();
  bary_of_interf_[old()].echange_espace_virtuel();
  normal_of_interf_ns_[old()].echange_espace_virtuel();
  bary_of_interf_ns_[old()].echange_espace_virtuel();

  surface_vapeur_par_face_[old()].echange_espace_virtuel();
  surface_vapeur_par_face_ns_[old()].echange_espace_virtuel();
  for (int c = 0; c < 3; c++)
    {
      barycentre_vapeur_par_face_[old()][c].echange_espace_virtuel();
      barycentre_vapeur_par_face_ns_[old()][c].echange_espace_virtuel();
    }

  indicatrice_ns_[old()].echange_espace_virtuel(indicatrice_ns_[old()].ghost());
  groups_indicatrice_ns_[old()].echange_espace_virtuel();
}


void IJK_Interfaces::calculer_phi_repuls_par_compo(
  FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& surf_par_compo,
  FixedVector<FixedVector<IJK_Field_double, max_authorized_nb_of_components_>, 3>& source_interf_par_compo,
  FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& phi_par_compo,
  FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& repuls_par_compo,
  const DoubleTab& gravite,
  const double delta_rho,
  const double sigma,
  const double time,
  const int itstep
)
{
  field_repulsion_.data() = -1.;

  ArrOfDouble potentiels_sommets;
  ArrOfDouble repulsions_sommets;
  calculer_phi_repuls_sommet(potentiels_sommets, repulsions_sommets, gravite, delta_rho, sigma, time, itstep);
  val_par_compo_in_cell_computation_.calculer_moy_field_sommet_par_compo(
    potentiels_sommets, phi_par_compo);


  if (!maillage_ft_ijk_.Surfactant_facettes().get_disable_surfactant() or (maillage_ft_ijk_.Surfactant_facettes().get_disable_surfactant() and use_tryggvason_interfacial_source_ ))
    {
      // on passe ici si sigma = variable
      // ou si sigma = cte, mais terme source interf = formulation Tryggvason
      DoubleTab interfacial_source_term_sommet = maillage_ft_ijk_.update_sigma_and_interfacial_source_term_sommet(ref_domaine_, true, use_tryggvason_interfacial_source_, sigma);
      ArrOfDouble interfacial_source_term_sommet_dir, unite;
      interfacial_source_term_sommet_dir.resize(maillage_ft_ijk_.nb_sommets());
      unite.resize(maillage_ft_ijk_.nb_sommets());
      for (int som = 0; som < maillage_ft_ijk_.nb_sommets(); som++)
        unite(som) = 1. ;
      val_par_compo_in_cell_computation_.calculer_somme_field_sommet_par_compo(unite, surf_par_compo);

      for (int dir = 0; dir < 3; dir++)
        {
          for (int som = 0; som < maillage_ft_ijk_.nb_sommets(); som++)
            interfacial_source_term_sommet_dir(som)=interfacial_source_term_sommet(som, dir);
          val_par_compo_in_cell_computation_.calculer_moy_field_sommet_par_compo(interfacial_source_term_sommet_dir, source_interf_par_compo[dir]);
        }
    }

  val_par_compo_in_cell_computation_.calculer_moy_field_sommet_par_compo(
    repulsions_sommets, repuls_par_compo);

  const Domaine_IJK& split = ref_domaine_;
  const int ni = split.get_nb_elem_local(DIRECTION_I);
  const int nj = split.get_nb_elem_local(DIRECTION_J);
  const int nk = split.get_nb_elem_local(DIRECTION_K);

  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          field_repulsion_(i, j, k) = repuls_par_compo[0](i, j, k);
        }
}


void IJK_Interfaces::calculer_phi_repuls_sommet(
  ArrOfDouble& potentiels_sommets,
  ArrOfDouble& repulsions_sommets,
  const DoubleTab& gravite,
  const double delta_rho,
  const double sigma,
  const double time,
  const int itstep
)
{

  // Initialisation forcee a -1 :

  const Domaine_IJK& geom = ref_domaine_.valeur();

  // calculer la courbure et le terme de gravite aux sommets du maillage
  // lagrangien On appelle ce terme "phi", potentiel aux sommets
  const Maillage_FT_IJK& mesh = maillage_ft_ijk_;
  const ArrOfDouble& courbure = mesh.get_update_courbure_sommets();
  potentiels_sommets = courbure;
  const double dxi = geom.get_constant_delta(DIRECTION_I);
  const double dxj = geom.get_constant_delta(DIRECTION_J);
  const double dxk = geom.get_constant_delta(DIRECTION_K);
  const double vol_cell = dxi * dxj * dxk;
  const DoubleTab& sommets = mesh.sommets();
  const int nb_som = potentiels_sommets.size_array();
  //ducluz : c'est ici qu'on multiplie le potentiel par sigma
  //Il faut quil passe en sigma variable ici dans le cas de surfactant
  if (!maillage_ft_ijk_.Surfactant_facettes().get_disable_surfactant())
    {
      //DoubleTab interfacial_source_term_sommet = maillage_ft_ijk_.update_sigma_and_interfacial_source_term_sommet(ref_domaine_, false, use_tryggvason_interfacial_source_);
      const ArrOfDouble& sigma_sommets = maillage_ft_ijk_.Surfactant_facettes().get_sigma_sommets();
      for (int i = 0; i < nb_som; i++)
        {
          potentiels_sommets[i] *= sigma_sommets[i];
        }
    }
  else
    {
      potentiels_sommets *= sigma;
    }


  // Terme source de gravite:
  //  on ajoute le potentiel de gravite = +/- delta_rho * (vecteur_g scalaire
  //  Ox), Ox est le vecteur qui va de l'origine au sommet de l'interface
  repulsions_sommets.resize_array(nb_som);

  // Ajout du terme de gravite:

  ArrOfIntFT compo_connexe_sommets;
  mesh.calculer_compo_connexe_sommets(compo_connexe_sommets);
  if (terme_gravite_ == GRAVITE_GRAD_I)
    {
      // Nouvelle version :
      DoubleTab deplacement;
      // le tableau contient l'encodage pour le deplacement que l'on va decoder :
      calculer_deplacement_from_code_compo_connexe_negatif(
        mesh, deplacement, bounding_box_NS_domain_
      );

      for (int i = 0; i < nb_som; i++)
        {
          double correction_potentiel_deplacement = 0.;
          if (compo_connexe_sommets[i] < 0)
            {
              // Bulle ghost : correction du deplacement:
              correction_potentiel_deplacement = -(
                                                   gravite(0,0) * deplacement(i, 0) +
                                                   gravite(0,1) * deplacement(i, 1) +
                                                   gravite(0,2) * deplacement(i, 2)
                                                 );
            }
          const Vecteur3 coord(sommets, i);
          double p = (
                       gravite(0,0) * coord[0] +
                       gravite(0,1) * coord[1] +
                       gravite(0,2) * coord[2] +
                       correction_potentiel_deplacement
                     );
          potentiels_sommets[i] -= p * delta_rho;
        }
    }

  // Terme source de repulsion des bulles:
  if (compute_distance_autres_interfaces_ || (delta_p_max_repulsion_ > 0. && portee_force_repulsion_ > 0.))
    {
      double vrx = DMAXFLOAT;
      double vry = DMAXFLOAT;
      double vrz = DMAXFLOAT;
      DoubleTab vr_to_closer; // The velocity of the closest neighbour
      calculer_distance_autres_compo_connexe2(distance_autres_interfaces_, vr_to_closer);
      if (delta_p_max_repulsion_ > 0. && portee_force_repulsion_ > 0.)
        {
          double dmin = portee_force_repulsion_;
          const ArrOfDouble& z_grid_nodes = geom.get_node_coordinates(DIRECTION_K);
          const int nznodes = z_grid_nodes.size_array();
          const double zmin = z_grid_nodes[0];
          const double zmax = z_grid_nodes[nznodes - 1];
          double zsmin = (zmax - zmin) / 2.;
          for (int i = 0; i < nb_som; i++)
            {
              double d = distance_autres_interfaces_[i];
              if (d < dmin)
                {
                  dmin = d;
                  vrx = vr_to_closer(i, 0);
                  vry = vr_to_closer(i, 1);
                  vrz = vr_to_closer(i, 2);
                }
              // dmin = std::min(dmin, d);
              double phi = 0.;
              if (active_repulsion_paroi_ && !geom.get_periodic_flag(DIRECTION_K))
                {
                  // Repulsion des parois haute et basse:
                  const double z = sommets(i, 2);
                  double dzs = std::min(z - zmin, zmax - z);
                  zsmin = std::min(zsmin, dzs);
                  if (delta_p_wall_max_repulsion_ > 0. && portee_wall_repulsion_ > 0.)
                    {
                      dzs = std::min(dzs, portee_wall_repulsion_);
                      // Ici, la repulsion paroi est donnee par une loi specifique :
                      phi = (
                              delta_p_wall_max_repulsion_ *
                              (portee_wall_repulsion_ - dzs) /
                              portee_wall_repulsion_
                            );
                    }
                  else
                    {
                      // Ici, la repulsion paroi est mise dans d, donc c'est la meme loi
                      // que pour l'inter-bulles :
                      d = std::min(d, dzs);
                    }
                }
              // Le "+=" sert dans le cas ou on a mis 2 lois differentes :
              phi += delta_p_max_repulsion_ * (portee_force_repulsion_ - d) / portee_force_repulsion_;
              potentiels_sommets[i] -= phi;
              repulsions_sommets[i] = -phi * vol_cell;
            }
          double dmin_local = dmin;
          dmin = Process::mp_min(dmin);
          zsmin = Process::mp_min(zsmin);

          int flag = 0;
          int iproc = Process::nproc();
          if (fabs(dmin_local - dmin) < 1e-12)
            {
              flag = 1;
              iproc = Process::me();
            }
          const int sumflag = Process::check_int_overflow(Process::mp_sum(flag));
          if (sumflag != 1)
            {
              Cerr << "Warning. There were equalities ("
                   << Process::mp_sum(flag)
                   << ") in dmin computation. "
                   << "2 equalities seems possible if it's between two markers "
                   "separated by a proc frontieer. "
                   << finl;
              // En cas d'egalite, il faut trancher pour que tout le monde soit
              // d'accord sur les envois/receptions :
              // on prend le petit proc:
              iproc = Process::mp_min(iproc);
            }
          else
            {
              // S'il n'y a pas d'egalite, il faut faire connaitre iproc a tous pour
              // l'envoi : C'est une sorte d'envoyer_broadcast(iproc, iproc), sauf
              // qu'evidemment, ca n'a pas de sens.
              iproc = Process::mp_min(iproc);
              // Voila. Au moins tout le monde connait le proc du min...
            }
          envoyer_broadcast(vrx, iproc);
          envoyer_broadcast(vry, iproc);
          envoyer_broadcast(vrz, iproc);

          // Impression dans le fichier _dmin.out :
          const double dx = geom.get_constant_delta(DIRECTION_I);
          const double dy = geom.get_constant_delta(DIRECTION_J);
          const double dz = geom.get_constant_delta(DIRECTION_K);
          if (Process::je_suis_maitre())
            {
              double norm_ur = std::sqrt(vrx * vrx + vry * vry + vrz * vrz);
              const double dt = ref_ijk_ft_->schema_temps_ijk().get_timestep();
              double cfl = dt * std::min(std::min(std::fabs(vrx) / dx, std::fabs(vry) / dy), std::fabs(vrz) / dz);
              int reset = (!reprise_) && (itstep == 0);
              SFichier fic = Ouvrir_fichier(
                               "_dmin.out", "tstep\ttime\t\tdmin\tzsmin\tur\tCFL_local", reset
                             );
              fic << itstep << " " << time << " " << dmin << " " << zsmin;
              fic << " " << norm_ur << " " << cfl << finl;
              fic.close();
            }
        }
    }
  potentiels_sommets *= dxi * dxj * dxk;

  // Au pire, on prend la plus grande largeur de maille eulerienne :
  double lg_euler = std::max(std::max(dxi, dxj), dxk);
  double lg_lagrange = mesh.minimum_longueur_arrete();
  const double sqrt_3 = 1.7320508075688772;
  bool ok = (lg_lagrange > lg_euler * sqrt_3);
  Cerr << "Test de la taille de maille eulerienne : lg_lagrange=" << lg_lagrange << " > (lg_euler=" << lg_euler
       << ")*1.7320... ? " << (ok ? "ok" : "ko") << " (ratio : " << lg_lagrange / (lg_euler * sqrt_3)
       << (ok ? " >" : " <") << "1 )" << finl;
}
