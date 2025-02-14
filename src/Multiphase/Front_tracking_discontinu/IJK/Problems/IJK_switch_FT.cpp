//TRUST_NO_INDENT
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

#include <IJK_switch_FT.h>
#include <init_forcage_THI.h>
#include <corrections_qdm.h>
#include <Force_sp.h>


Implemente_instanciable(Switch_FT_double, "Switch_FT_double", Switch_double);

Sortie & Switch_FT_double::printOn(Sortie&s) const
{
  return s;
}

Entree & Switch_FT_double::readOn(Entree&s)
{
  return s;
}

void Switch_FT_double::set_param(Param& param)
{
  Switch_double::set_param(param);

  // Parametres pour le FT:
  param.ajouter("interfaces", &interfaces_);
  param.ajouter("old_ijk_splitting_ft_extension", &old_ijk_splitting_ft_extension_);

  // Parametres pour la thermique:
  param.ajouter("thermals", &thermals_);

  /*
   * GAB : gabriel.ramirez@cea.fr
   * Parametres pour le forcage spectral
   * Voir reprendre probleme dans Probleme_FTD_IJK_base.cpp
   */
  param.ajouter("forcage", &old_forcage_);
  // Parametres pour la correction de qdm 
  param.ajouter("corrections_qdm", &old_qdm_corrections_);
  param.ajouter("reprise_vap_velocity_tmoy", &vap_velocity_tmoy_);
  param.ajouter("reprise_liq_velocity_tmoy", &liq_velocity_tmoy_);
  param.ajouter("vitesse_upstream_reprise", &vitesse_upstream_reprise_);
  param.ajouter("velocity_bubble_old", &velocity_bubble_old_);
}

void Switch_FT_double::prepare_run()
{
  Switch_double::prepare_run();
  Nom prefix = dirname(nom_reprise_);
  interfaces_.set_fichier_reprise(prefix + interfaces_.get_fichier_reprise());
  if (!thermals_.est_vide())
	  thermals_.set_fichier_reprise(prefix + thermals_.get_fichier_reprise());
  fichier_old_vitesse_ = prefix + fichier_old_vitesse_;
}


void Switch_FT_double::initialise()
{
  Cout << que_suis_je() <<"::initialise() Pb of type Probleme_FTD_IJK_base detected." << finl;

  // Probleme of type FT:
  // old_mesh_ and new_mesh_ are acutally splittings...
  const double dx = old_mesh_.get_constant_delta(DIRECTION_I);
  const double dy = old_mesh_.get_constant_delta(DIRECTION_J);

  const double new_dx = new_mesh_.get_constant_delta(DIRECTION_I);
  const double new_dy = new_mesh_.get_constant_delta(DIRECTION_J);
  const Nom& vdf_name = new_mesh_.le_nom()+Nom("_VDF");

  const double old_to_new_ratio = std::min(dx/new_dx, dy/new_dy);
  const int new_ijk_splitting_ft_extension = (int) std::lrint(std::ceil(old_ijk_splitting_ft_extension_* old_to_new_ratio));

  Cerr << "Extended splitting dimensions. old = " << old_ijk_splitting_ft_extension_
       << " new = " << new_ijk_splitting_ft_extension << finl;
  Cerr << "Construction du domaine VDF..." << finl;
  Domaine_IJK splitting_ft;
  build_extended_splitting(new_mesh_, splitting_ft, new_ijk_splitting_ft_extension);

  // Le probleme ft disc qui porte le maillage vdf pour les algorithmes front-tracking
  Probleme_base& refprobleme_ft_disc = creer_domaine_vdf(splitting_ft, vdf_name);
  const Domaine_dis_base& domaine_dis = refprobleme_ft_disc.domaine_dis();

  interfaces_.initialize(splitting_ft /* splitting_FT */,
                         new_mesh_ /* splitting_NS */,
                         domaine_dis,
                         0,
                         false,
                         true);
  //interfaces_.associer_switch(*this);
  interfaces_.set_reprise(1);
  interfaces_.lire_maillage_ft_dans_lata();
  
  Switch_double::initialise();

  // GAB
  const int nproc_tot = Process::nproc();
  old_forcage_.compute_initial_chouippe(nproc_tot,old_mesh_,old_ni_,old_nj_,old_nk_,splitting_ft,nom_sauvegarde_);
  new_forcage_.compute_initial_chouippe(nproc_tot,new_mesh_,new_ni_,new_nj_,new_nk_,splitting_ft,nom_sauvegarde_);
  Cout << "new_forcage_.get_semi_gen()" << new_forcage_.get_semi_gen() << finl;
  Cout << "old_forcage_.get_semi_gen()" << old_forcage_.get_semi_gen() << finl;
  Cout << "new_forcage_.get_b_flt()" << new_forcage_.get_b_flt() << finl;
}


int Switch_FT_double::init_thermique()
// The Thermique.initialize does both the allocate and the initialize:
{
  return 0;
}

int Switch_FT_double::init_thermals()
{
	// init_thermals() does both the allocate and the initialize()
	// thermals_.associer_switch(*this);
  return thermals_.init_switch_thermals(old_mesh_);
}

void Switch_FT_double::prepare_thermique(const Nom lata_name)
{
}

void Switch_FT_double::prepare_thermals(const Nom lata_name)
{
	thermals_.prepare_thermals(lata_name);
}

// flag and_lata to know if we also create the associated lata
void Switch_FT_double::ecrire_fichier_reprise(const char *fichier_sauvegarde, const bool and_lata)
{
  Nom lata_name(fichier_sauvegarde);
  lata_name += ".lata";

  if (Process::je_suis_maitre())
    {
      Cerr << "T= " << current_time_ << " Checkpointing dans le fichier " << fichier_sauvegarde << finl;
      SFichier fichier(fichier_sauvegarde);
      fichier.precision(17);
      fichier.setf(std::ios_base::scientific);
      fichier << "{\n"
          << " tinit " << current_time_ << "\n"
          << " terme_acceleration_init " << terme_source_acceleration_ << "\n"
          << " reprise_vap_velocity_tmoy " << vap_velocity_tmoy_ << "\n"
          << " reprise_liq_velocity_tmoy " << liq_velocity_tmoy_ << "\n"
		  << " velocity_bubble_old " << velocity_bubble_old_ << "\n"
		  << " vitesse_upstream_reprise " << vitesse_upstream_reprise_ << "\n"
          << " fichier_reprise_vitesse " << lata_name << "\n"
          << " timestep_reprise_vitesse 1\n";

      interfaces_.set_fichier_sauvegarde(lata_name);
      Cerr << "  saving interfaces... " << finl;
      fichier << " interfaces " << interfaces_  << "\n";
      fichier <<  " forcage " << new_forcage_;
      fichier <<  " corrections_qdm " << new_qdm_corrections_;

      Cerr << "  potentially saving temperature fields... " << finl;

      thermals_.ecrire_fichier_reprise(fichier, lata_name);

      fichier << "}\n";
    }

  if (and_lata)
    {
      ecrire_header_lata(lata_name);
      write_velocity(lata_name);
    }
}

void Switch_FT_double::ecrire_header_lata(const Nom lata_name) // const
{
  Switch_double::ecrire_header_lata(lata_name);

  Cout << "Adding interfaces into " << lata_name << finl;
  interfaces_.sauvegarder_interfaces(lata_name);
  // Le bloc suivant aurait pour consequence d'ecrire la temperature initiale (sur le maillage relu)
  // (car c'est pour lui qu'on a fait un allocate dans l'objet thermique pour le relire)
  // dans le fichier de sauvegarde des champs exrits pour la reprise. Ce n'est pas ce que l'on veut.
  // On avait pu se permettre cela pour l'interface car on ecrit directement celui qu'on a relu, meme
  // s'il n'a pas le bon support (le vieux geom).
}

void Switch_FT_double::set_param_reprise(Param& param)
{
  Switch_double::set_param_reprise(param);
  param.ajouter("interfaces", & interfaces_);
  param.ajouter("thermals", & thermals_);
  // GAB : gabriel.rmairez@cea.fr
  /* Voir reprendre probleme dans Probleme_FTD_IJK_base.cpp */
  param.ajouter("forcage", &new_forcage_);
  param.ajouter("corrections_qdm", &new_qdm_corrections_);
  param.ajouter("reprise_vap_velocity_tmoy", &vap_velocity_tmoy_);
  param.ajouter("reprise_liq_velocity_tmoy", &liq_velocity_tmoy_);
  // fin GAB : gabriel.ramirez@cea.fr
  param.ajouter("vitesse_upstream_reprise", &vitesse_upstream_reprise_);
  param.ajouter("velocity_bubble_old", &velocity_bubble_old_);
}

void Switch_FT_double::compute_and_write_extra_fields(const Nom& lata_name,
                                                      DoubleTab& coeff_i, IntTab Indice_i,
                                                      DoubleTab& coeff_j, IntTab Indice_j,
                                                      DoubleTab& coeff_k, IntTab Indice_k)
{
  thermals_.compute_new_thermal_field((*this),
  																		new_mesh_,
																			lata_name,
																			coeff_i,
																			Indice_i,
																			coeff_j,
																			Indice_j,
																			coeff_k,
																			Indice_k);
}

void Switch_FT_double::compute_and_write_extra_fields_direct(SFichier& file,
																														 DoubleTab& coeff_i, IntTab Indice_i,
																														 DoubleTab& coeff_j, IntTab Indice_j,
																														 DoubleTab& coeff_k, IntTab Indice_k)
{
	// Use writing_direct = 0 instead...
}



