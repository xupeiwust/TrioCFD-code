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
#ifndef Statistiques_dns_ijk_FT_H
#define Statistiques_dns_ijk_FT_H
#include <IJK_Field_vector.h>
#include <IJK_Field.h>
#include <Statistiques_dns_ijk.h>
#include <TRUSTArrays.h>
#include <TRUST_Ref.h>

class Probleme_FTD_IJK_base;
class Domaine_IJK;

class Statistiques_dns_ijk_FT : public Statistiques_dns_ijk
{
  Declare_instanciable_sans_constructeur(Statistiques_dns_ijk_FT);
public:
  Statistiques_dns_ijk_FT(); // Je ne sais pas compiler le Vect(Statistiques_dns_ijk_FT) sans lui...
  void associer_probleme(const Probleme_FTD_IJK_base&);
  using Statistiques_dns_ijk::initialize;
  void initialize(const Probleme_FTD_IJK_base& ijk_ft,const Domaine_IJK&);
  void initialize(const Probleme_FTD_IJK_base& ijk_ft,const Domaine_IJK& splitting,
                  const int check_stats);
  void alloc_fields();
  Sortie& completer_print(Sortie& os) const override;
  void completer_read(Param& param) override;
  int lire_motcle_non_standard(const Motcle& mot, Entree& is) override;

  void update_stat(Probleme_FTD_IJK_base& cas, const double dt);

  void postraiter(Sortie&, int flag_valeur_instantanee = 0) const;
  void postraiter_thermique(const double t) const;
  double compute_desequil_alpha(const Domaine_IJK& geom_NS,
                                const double portee_wall_repulsion) const;

  const IJK_Field_vector3_double& get_IJK_field_vector(const Motcle& nom) override;
  const IJK_Field_double& get_IJK_field(const Motcle& nom) override;

protected:
  int check_stats_;
  int nb_thermal_fields_; // Number of objects thermique_ in the list
  int nvalt_;             // Number of variables post-processed per field
  // Last instantaneous value of the space average (only on processor 0)
  DoubleTab moyenne_spatiale_instantanee_temperature_; // (i,j,k) ->  (nvalt_,nb_elem_k_tot,nb_thermal_fields_)
  // Temporal integral of statistics variables
  DoubleTab integrale_temporelle_temperature_; // (i,j,k) ->  (nvalt_,nb_elem_k_tot,nb_thermal_fields_)
  VECT(Nom) noms_moyennes_temperature_;
};
#endif
