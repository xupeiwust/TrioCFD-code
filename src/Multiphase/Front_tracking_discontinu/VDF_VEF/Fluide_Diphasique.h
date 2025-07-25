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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Fluide_Diphasique.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Fluide_Diphasique_included
#define Fluide_Diphasique_included

#include <Fluide_Incompressible.h>
#include <Milieu_base.h>
#include <Objet_U.h>
#include <Milieu_base.h>
#include <Solid_Particle_base.h>

class Fluide_Diphasique: public Milieu_base
{
  Declare_instanciable_sans_constructeur(Fluide_Diphasique);
public:

  Fluide_Diphasique()
  {
    indic_rayo_ = NONRAYO;
    formule_mu_ = "standard";
    is_solid_particle_ = false;
  }

  const Fluide_Incompressible& fluide_phase(int la_phase) const;
  double sigma() const;
  double chaleur_latente() const;
  int formule_mu() const;

  // Surcharge des methodes standard du milieu_base :
  void set_param(Param& param) override;
  void verifier_coherence_champs(int& err, Nom& message) override;
  int initialiser(const double temps) override;
  void mettre_a_jour(double temps) override;
  void discretiser(const Probleme_base& pb, const Discretisation_base& dis) override;

  // L'appel a ces methodes est invalide et genere une erreur
  const Champ_base& masse_volumique() const override { return invalid_<const Champ_base&>(__func__); }
  Champ_base& masse_volumique() override { return invalid_<Champ_base&>(__func__); }
  const Champ_Don_base& diffusivite() const override { return invalid_<const Champ_Don_base&>(__func__); }
  Champ_Don_base& diffusivite() override { return invalid_<Champ_Don_base&>(__func__); }
  const Champ_Don_base& conductivite() const override  { return invalid_<const Champ_Don_base&>(__func__); }
  Champ_Don_base& conductivite() override { return invalid_<Champ_Don_base&>(__func__); }
  const Champ_Don_base& capacite_calorifique() const override  { return invalid_<const Champ_Don_base&>(__func__); }
  Champ_Don_base& capacite_calorifique() override { return invalid_<Champ_Don_base&>(__func__); }
  const Champ_Don_base& beta_t() const override  { return invalid_<const Champ_Don_base&>(__func__); }
  Champ_Don_base& beta_t() override { return invalid_<Champ_Don_base&>(__func__); }

  const bool& get_is_solid_particle() const { return is_solid_particle_; }
  const int& get_id_fluid_phase() const { return id_fluid_phase_; }

protected:
  OWN_PTR(Milieu_base) phase0_, phase1_;
  OWN_PTR(Champ_Don_base) sigma_; // Tension de surface (J/m^2)
  OWN_PTR(Champ_Don_base) chaleur_latente_; // Enthalpie de changement de phase h(phase1_) - h(phase0_) (J/kg/K)
  Motcle formule_mu_; // Formule utilisee pour le calcul de la moyenne de mu
  bool is_solid_particle_; // True if phase0_ or phase1_ is a Solid_Particle
  int id_fluid_phase_; // number (0 or 1) of the Fluid_Incompressible phase
  template <typename RETURN_TYPE>
  RETURN_TYPE invalid_(const char * nom_funct) const
  {
    Cerr << "Invalid call to the method Fluide_Diphasique::" << nom_funct << " !!!" << finl;
    throw;
  }
};

#endif /* Fluide_Diphasique_included */
