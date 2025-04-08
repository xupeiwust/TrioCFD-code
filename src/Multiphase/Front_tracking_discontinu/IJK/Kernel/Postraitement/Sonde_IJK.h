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

#ifndef Sonde_IJK_included
#define Sonde_IJK_included

#include <Sonde.h>
#include <Postraitement.h>
#include <IJK_Field.h>
#include <TRUST_Ref.h>

class Probleme_FTD_IJK_base;

/*! @brief class Sonde_IJK
 */
class Sonde_IJK : public Sonde
{

  Declare_instanciable( Sonde_IJK ) ;

public :

  void completer() override;

protected :
  const Domaine& get_domaine_geom() const override;
  const Noms get_noms_champ() const override;
  int get_nb_compo_champ() const override;
  double get_temps_champ() const override;
  void validate_type(const Motcle& loc) const override;
  void validate_position() const override;
  void create_champ_generique(Entree& is, const Motcle& motlu) override { /* do nothing */ }
  void fix_probe_position() override;
  void fix_probe_position_grav() override;
  void fill_local_values() override;
  void update_source(double un_temps) override { /* do nothing */ }

  OBS_PTR(Probleme_FTD_IJK_base) ref_ijk_ft_;
  OBS_PTR(IJK_Field_double) ref_ijk_field_;
  IJK_Field_double tmp_storage_;

  int ncomp = 1;                           ///< Number of the component to handle
  int nbre_points_tot;

private:
  void fix_probe_position_generic(Domaine_IJK::Localisation loc);

  /** In IJK we need to store the initial name of the requested field because of the switch done on some fields
   *  (see Champs_compris_IJK::switch_ft_fields()). This switch might rename the field, but we want this to stay
   *  constant for a probe!
   */
  Nom field_name_;
};

#endif /* Sonde_IJK_included */
