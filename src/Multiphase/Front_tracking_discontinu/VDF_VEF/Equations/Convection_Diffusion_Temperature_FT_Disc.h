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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Convection_Diffusion_Temperature_FT_Disc.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/15
//
//////////////////////////////////////////////////////////////////////////////
#ifndef Convection_Diffusion_Temperature_FT_Disc_included
#define Convection_Diffusion_Temperature_FT_Disc_included

#include <Convection_Diffusion_Temperature.h>


#include <Assembleur_base.h>
#include <TRUST_Ref.h>

class Fluide_Diphasique;
class Navier_Stokes_std;
class Transport_Interfaces_FT_Disc;


class Convection_Diffusion_Temperature_FT_Disc: public Convection_Diffusion_Temperature
{
  Declare_instanciable_sans_constructeur(Convection_Diffusion_Temperature_FT_Disc);
public:
  friend class Post_Processing_Hydrodynamic_Forces;

  Convection_Diffusion_Temperature_FT_Disc();
  void set_param(Param& titi) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  virtual void preparer_pas_de_temps();
  virtual void corriger_pas_de_temps(double dt);
  void compute_divergence_free_velocity_extension();
  DoubleTab&   derivee_en_temps_inco(DoubleTab&) override;
  void         mettre_a_jour(double temps) override;
  int preparer_calcul() override;
  double get_flux_to_face(const int num_face) const;
  double get_Twall_at_face(const int num_face) const;
  double get_Twall_at_elem(const int elem) const;
  const double& get_tsat_constant() const;
  void get_flux_and_Twall(const int num_face,
                          double& flux, double& Twall) const;
  double get_Twall(const int num_face) const;
  virtual void suppression_interfaces(const IntVect& num_compo, const ArrOfInt& flags_compo_a_supprimer, int nouvelle_phase);
  void                associer_milieu_base(const Milieu_base& milieu) override;
  Milieu_base&        milieu() override;
  const Milieu_base& milieu() const override;
  const Champ_base& vitesse_pour_transport() const override;
  void    discretiser() override;

  int get_phase() const;
  void discretiser_assembleur_pression();
  void completer() override;


  const DoubleTab& get_mpoint() const
  {
    return mpoint_->valeurs();
  };

  const Transport_Interfaces_FT_Disc& get_eq_interface() const
  {
    return ref_eq_interface_.valeur();
  };
  Transport_Interfaces_FT_Disc& eq_interface()
  {
    return ref_eq_interface_.valeur();
  };

  ArrOfInt& mixed_elems()
  {
    return mixed_elems_;
  };
  ArrOfDouble& lost_fluxes()
  {
    return lost_fluxes_;
  }
  void calculer_grad_t();
  void calculer_mpoint(Champ_base& mpoint);
  void calculer_mpoint();

  enum class Thermal_correction_discretization_method { P1, ELEM_DIPH, P1_ALL };
  const Thermal_correction_discretization_method& get_thermal_correction_discretization_method()
  const  { return thermal_correction_discretization_method_; }


protected:
  void correct_mpoint();
  void set_is_solid_particle(const bool is_solid_particle) { is_solid_particle_=is_solid_particle; }

  // Quelle phase cette equation concerne-t-elle ? 0 ou 1
  int phase_;
  //GB : Ajout de variables :
  int stencil_width_;
  int correction_courbure_ordre_;
  int is_solid_particle_;
  Nom nom_sous_domaine_;
  double temp_moy_ini_;
  bool maintien_temperature_;
  bool is_prescribed_mpoint_;
  double prescribed_mpoint_;
  ArrOfInt correction_mpoint_diff_conv_energy_ ; // on attend trois flags 0 ou 1

  OBS_PTR(Fluide_Diphasique) fluide_dipha_;

  // Champs compris par le postraitement
  LIST(OBS_PTR(Champ_base)) liste_champs_compris_;

  // Reference a l'equation de transport de l'interface pour l'indicatrice de phase
  OBS_PTR(Transport_Interfaces_FT_Disc) ref_eq_interface_;
  // Reference a l'equation de navier_stokes pour le champ de vitesse (convection)
  OBS_PTR(Navier_Stokes_std) ref_eq_ns_;

  // Gradient normal de temperature a l'interface phase 0
  // (grad T scalaire normale a l'interface, normale dirigee
  //  vers la phase 1)
  OWN_PTR(Champ_Fonc_base)  grad_t_;
  OWN_PTR(Champ_Fonc_base)  mpoint_;
  OWN_PTR(Champ_Fonc_base)  mpoint_uncorrected_;
  OWN_PTR(Champ_Inc_base) vitesse_convection_;

  // To make a divergence-free velocity extension :
  int divergence_free_velocity_extension_;
  OWN_PTR(Assembleur_base) assembleur_pression_;
  OWN_PTR(Champ_Inc_base) la_pression; // Of course, it's a fake :D
  OWN_PTR(Champ_Inc_base) gradient_pression_;
  OWN_PTR(Champ_Inc_base) divergence_delta_U;
  SolveurSys solveur_pression_;
  Matrice matrice_pression_;
  OWN_PTR(Domaine_Cl_dis_base) zcl_fictitious_;
  Noms name_bc_opening_pressure_; // Liste de noms pour laisser sortir la source de div(delta u)

  ArrOfInt mixed_elems_;
  ArrOfDouble lost_fluxes_;
  ArrOfDouble derivee_energy_;
  ArrOfInt mixed_elems_diffu_;
  ArrOfDouble lost_fluxes_diffu_;
  ArrOfInt mixed_elems_conv_;
  ArrOfDouble lost_fluxes_conv_;

  Thermal_correction_discretization_method thermal_correction_discretization_method_=
    Thermal_correction_discretization_method::P1_ALL;

};
#endif
