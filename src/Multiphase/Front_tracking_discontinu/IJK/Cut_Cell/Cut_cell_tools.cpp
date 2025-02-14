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

#include <Cut_cell_tools.h>
#include <IJK_Field.h>
#include <Cut_field.h>
#include <Cut_cell_FT_Disc.h>
#include <IJK_Lata_writer.h>

void dumplata_scalar_cut_cell(int cut_cell_activated,
                              const char *filename, const char *fieldname,
                              const std::shared_ptr<IJK_Field_double>& f, int step)
{
  if (cut_cell_activated)
    {
      Cut_field_double& cut_field_f = static_cast<Cut_field_double&>(*f);
      cut_field_f.echange_diph_vers_pure_cellules_finalement_pures();
      cut_field_f.remplir_tableau_pure_cellules_diphasiques(true);
    }

  dumplata_scalar(filename, fieldname, *f, step);

  if (cut_cell_activated)
    {
      Cut_field_double& cut_field_f = static_cast<Cut_field_double&>(*f);

      cut_field_f.get_cut_cell_disc().fill_buffer_with_variable(cut_field_f.diph_l_);
      dumplata_scalar(filename, Nom(fieldname) + "_CUT_L", cut_field_f.get_cut_cell_disc().get_write_buffer(), step);

      cut_field_f.get_cut_cell_disc().fill_buffer_with_variable(cut_field_f.diph_v_);
      dumplata_scalar(filename, Nom(fieldname) + "_CUT_V", cut_field_f.get_cut_cell_disc().get_write_buffer(), step);
    }
}


void lire_dans_lata_cut_cell(int cut_cell_activated,
                             const char *filename_with_path, int tstep, const char *geometryname, const char *fieldname,
                             std::shared_ptr<IJK_Field_double>& f)
{
  lire_dans_lata(filename_with_path, tstep, geometryname, fieldname, *f);

  if (cut_cell_activated)
    {
      Cut_field_double& cut_field_f = static_cast<Cut_field_double&>(*f);

      lire_dans_lata(filename_with_path, tstep, geometryname, Nom(fieldname) + "_CUT_L", cut_field_f.get_cut_cell_disc().get_write_buffer());
      cut_field_f.get_cut_cell_disc().fill_variable_with_buffer(cut_field_f.diph_l_);

      lire_dans_lata(filename_with_path, tstep, geometryname, Nom(fieldname) + "_CUT_V", cut_field_f.get_cut_cell_disc().get_write_buffer());
      cut_field_f.get_cut_cell_disc().fill_variable_with_buffer(cut_field_f.diph_v_);

      cut_field_f.echange_pure_vers_diph_cellules_initialement_pures();
    }
}


void allocate_velocity(Cut_field_vector3_double& v, const Domaine_IJK& s, int ghost, double DU)
{
  v.get_ptr(0) = std::make_shared<Cut_field_double>();
  v.get_ptr(1) = std::make_shared<Cut_field_double>();
  v.get_ptr(2) = std::make_shared<Cut_field_double>();

  v[0].allocate(s, Domaine_IJK::FACES_I, ghost);
  v[1].allocate(s, Domaine_IJK::FACES_J, ghost);
  v[2].allocate(s, Domaine_IJK::FACES_K, ghost);
  v[0].get_shear_BC_helpler().set_dU_(DU);
  v[1].get_shear_BC_helpler().set_dU_(0.);
  v[2].get_shear_BC_helpler().set_dU_(0.);
}

void allocate_velocity_persistant(Cut_cell_FT_Disc& cut_cell_disc, Cut_field_vector3_double& v, const Domaine_IJK& s, int ghost, double DU)
{
  allocate_velocity(v, s, ghost, DU);

  v[0].associer_persistant(cut_cell_disc);
  v[1].associer_persistant(cut_cell_disc);
  v[2].associer_persistant(cut_cell_disc);
}

void allocate_velocity_ephemere(Cut_cell_FT_Disc& cut_cell_disc, Cut_field_vector3_double& v, const Domaine_IJK& s, int ghost, double DU)
{
  allocate_velocity(v, s, ghost, DU);

  v[0].associer_ephemere(cut_cell_disc);
  v[1].associer_ephemere(cut_cell_disc);
  v[2].associer_ephemere(cut_cell_disc);
}

void allocate_velocity_paresseux(Cut_cell_FT_Disc& cut_cell_disc, Cut_field_vector3_double& v, const Domaine_IJK& s, int ghost, double DU)
{
  allocate_velocity(v, s, ghost, DU);

  v[0].associer_paresseux(cut_cell_disc);
  v[1].associer_paresseux(cut_cell_disc);
  v[2].associer_paresseux(cut_cell_disc);
}

