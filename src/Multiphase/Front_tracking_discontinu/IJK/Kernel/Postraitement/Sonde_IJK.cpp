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

#include <Sonde_IJK.h>

#include <IJK_Field.h>
#include <Probleme_FT_Disc_gen.h>
#include <Domaine.h>
#include <Domaine_VF.h>
#include <Champ_Generique_Interpolation.h>
#include <communications.h>
#include <Probleme_FTD_IJK_base.h>
#include <IJK_tools.h>

Implemente_instanciable( Sonde_IJK, "Sonde_IJK", Sonde ) ;

Sortie& Sonde_IJK::printOn(Sortie& s ) const
{
  return s << que_suis_je();
}

void Sonde_IJK::completer()
{
  ref_ijk_ft_ = ref_cast(Probleme_FTD_IJK_base, mon_post->probleme());
  ref_ijk_field_ = ref_ijk_ft_->get_IJK_field(nom_champ_lu_);
  field_name_ = ref_ijk_field_->le_nom();

  // Make sure the cell-projected counterpart of the initial field is registered for postprocessing:
  if(grav && ref_ijk_field_->get_localisation() != Domaine_IJK::ELEM)
    ref_ijk_ft_->get_post().register_one_field(ref_ijk_field_->le_nom(), "ELEM");

  initialiser();
}

Entree& Sonde_IJK::readOn( Entree& is )
{
  return Sonde::readOn(is);
}

/** Override. Exit if invalid type for the probe.
 */
void Sonde_IJK::validate_type(const Motcle& loc) const
{
  // Valid IJK localisation, all lower case!
  static constexpr auto VALID = {"point", "points", "segment", "plan", "volume", "segmentxdx",
                                 "planxdxdy", "circle", "segmentpoints"
                                };

  std::string s = loc.getString(); // copy
  // from cppreference - convert to lowercase:
  std::transform(s.begin(), s.end(), s.begin(),
  [](unsigned char c) { return std::tolower(c); }
                );

  if (std::find(VALID.begin(), VALID.end(), s) == VALID.end())
    Process::exit("ERROR: Invalid localisation for IJK probe!!");
}

/** Override. Exit if invalid localisation for the probe.
 */
void Sonde_IJK::validate_position() const
{
  if(chsom || gravcl || nodes || som)
    Process::exit("ERROR: only 'grav' position modificator is supported for IJK probe!");
}

const Domaine& Sonde_IJK::get_domaine_geom() const
{
  const Domaine_IJK& geom_ft = ref_ijk_ft_->get_domaine_ft();
  return ref_ijk_ft_->probleme(geom_ft).domaine();
}

const Noms Sonde_IJK::get_noms_champ() const
{
  Noms ret;
  ret.add(ref_ijk_field_->le_nom());
  return ret;
}

int Sonde_IJK::get_nb_compo_champ() const
{
  return ref_ijk_field_->nb_comp();
}

// Here no difference between time of the field, and current time as given by the problem
double Sonde_IJK::get_temps_champ() const
{
  return ref_ijk_ft_->schema_temps().temps_courant();
}


/** In IJK, velocity components might be interpolated/moved to element center.
 * This is equivalent to having a ELEM field.
 */
void Sonde_IJK::fix_probe_position_grav()
{
  fix_probe_position_generic(Domaine_IJK::ELEM);
}

void Sonde_IJK::fix_probe_position_generic(Domaine_IJK::Localisation loc)
{
  int nbre_points = les_positions_sondes_.dimension(0);
  const Domaine_IJK& geom_ft = ref_ijk_ft_->get_domaine_ft();
  const Domaine_IJK& field_splitting = ref_ijk_field_->get_domaine();
  Vecteur3 origin(0., 0., 0.), delta(0., 0., 0.);
  delta[0] = geom_ft.get_constant_delta(DIRECTION_I);
  delta[1] = geom_ft.get_constant_delta(DIRECTION_J);
  delta[2] = geom_ft.get_constant_delta(DIRECTION_K);

  Int3 offset;
  offset[0] =  geom_ft.get_offset_local(DIRECTION_I);
  offset[1] =  geom_ft.get_offset_local(DIRECTION_J);
  offset[2] =  geom_ft.get_offset_local(DIRECTION_K);

  // L'origine est sur un noeud. Donc que la premiere face en I est sur get_origin(DIRECTION_I)
  origin[0] = geom_ft.get_origin(DIRECTION_I)
              + ((loc==Domaine_IJK::FACES_J || loc==Domaine_IJK::FACES_K || loc==Domaine_IJK::ELEM) ? (delta[DIRECTION_I] * 0.5) : 0. ) ;
  origin[1] = geom_ft.get_origin(DIRECTION_J)
              + ((loc==Domaine_IJK::FACES_K || loc==Domaine_IJK::FACES_I || loc==Domaine_IJK::ELEM) ? (delta[DIRECTION_J] * 0.5) : 0. ) ;
  origin[2] = geom_ft.get_origin(DIRECTION_K)
              + ((loc==Domaine_IJK::FACES_I || loc==Domaine_IJK::FACES_J || loc==Domaine_IJK::ELEM) ? (delta[DIRECTION_K] * 0.5) : 0. ) ;

  for (int idx =0; idx < nbre_points; idx++)
    {
      const int num_elem = elem_[idx];
      if (num_elem == -1) continue; // Handle only my processor

      // L'element appartient a ce proc. On peut donc trouver son ijk...
      // Les autres positions seront mises a jour par les autres procs...
      const Int3 ijk = geom_ft.convert_packed_to_ijk_cell(num_elem);
      Int3 new_ijk(ijk);

      // Corrige la position de la sonde :
      for (int i = 0; i < 3; i++)
        les_positions_sondes_(idx, i) = origin[i]+(ijk[i]+offset[i])*delta[i];
      Cerr << "Sonde " << nom_ << " Point " << idx
           << " x= " << les_positions_sondes_(idx, 0)
           << " y= " << les_positions_sondes_(idx, 1)
           << " z= " << les_positions_sondes_(idx, 2) << finl;

      //
      // TODO ABN to be removed once no more extended domain:
      // TODO: GB : So strange, geom_ft.ft_extension() = 0 whereas field_splitting.ft_extension() !=0
      if ((field_splitting != geom_ft) &&
          (field_splitting.ft_extension() != 0))
        //(geom_ft.ft_extension() != 0))
        {
          // Les 2 splittings ne sont pas identiques il faut changer l'elem :
          for (int i = 0; i < 3; i++)
            {
              //Cerr << " " << (field_geom.get_origin(i)-geom.get_origin(i))/delta[i] << finl;
              if (geom_ft.get_periodic_flag(i))
                {
                  //assert(int((field_geom.get_origin(i)-geom.get_origin(i))/delta[i]) == ref_ijk_ft_->get_splitting_extension());
                  new_ijk[i] -= field_splitting.ft_extension();
                }
            }
          const int new_num_elem =  field_splitting.convert_ijk_cell_to_packed(new_ijk[0], new_ijk[1], new_ijk[2]);
          elem_[idx] = new_num_elem;
        }
    }
}

/**  Fix probe position. In IJK, if no position modifier is specified, we systematically
 * move the probe to the center of the element for a ELEM field, and to the face corresponding to the field
 * component if this is a FACE field.
 */
void Sonde_IJK::fix_probe_position()
{
  const Domaine_IJK::Localisation loc = ref_ijk_field_->get_localisation();
  fix_probe_position_generic(loc);
}

void Sonde_IJK::fill_local_values()
{
  valeurs_locales.resize(elem_.size(),1);  // WARNING: must be dim 2 to be compatible with TRUST ...
  const int nb_pts = elem_.size();

  // Make sure to call get_IJK_field() to trigger the update of the field:
  const IJK_Field_double* ijk_field = &(ref_ijk_ft_->get_IJK_field(field_name_));
  const Domaine_IJK& dom = ijk_field->get_domaine();

  // If necessary, interpolate to cell center:
  if(grav && ref_ijk_field_->get_localisation() != Domaine_IJK::ELEM)
    {
      if( tmp_storage_.ni() == 0 )
        tmp_storage_ = *ijk_field; // deep copy
      interpolate_to_center_compo(tmp_storage_, *ijk_field);
      ijk_field = &tmp_storage_;
    }

  for (int idx =0; idx < nb_pts; idx++)
    {
      const int num_elem = elem_[idx];
      if (num_elem < 0)
        {
          Cerr << "ERROR - Sonde_IJK::fill_local_values() - the probe '" << nom_ << "' defines point #" << idx << " which falls outside the domain!!" << finl;
          Process::exit();
        }
      const Int3 ijk = dom.convert_packed_to_ijk_cell(num_elem);
      valeurs_locales[idx] = (*ijk_field)(ijk[0],ijk[1],ijk[2]);
    }
}

