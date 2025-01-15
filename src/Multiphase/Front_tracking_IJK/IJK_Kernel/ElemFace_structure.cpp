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
/////////////////////////////////////////////////////////////////////////////
//
// File      : ElemFace_structure.cpp
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Interfaces.h>
#include <Cut_cell_FT_Disc.h>
#include <ElemFace_structure.h>
#include <IJK_Interfaces.h>

void ElemFace_Cut_cell_FT_Disc::initialise(IJK_Interfaces& interfaces, IJK_Splitting& splitting)
{
  elem_data_.initialise(interfaces, splitting, IJK_Splitting::ELEM);
  face_data_[0].initialise(interfaces, splitting, IJK_Splitting::FACES_I);
  face_data_[1].initialise(interfaces, splitting, IJK_Splitting::FACES_J);
  face_data_[2].initialise(interfaces, splitting, IJK_Splitting::FACES_K);
}

void ElemFace_Cut_cell_FT_Disc::initialise(IJK_Interfaces& interfaces, IJK_Splitting& splitting, const ElemFace_IJK_Field_double& old_indicatrice, const ElemFace_IJK_Field_double& next_indicatrice)
{
  elem_data_.initialise(interfaces, splitting, IJK_Splitting::ELEM, old_indicatrice.elem_data(), next_indicatrice.elem_data());
  face_data_[0].initialise(interfaces, splitting, IJK_Splitting::FACES_I, old_indicatrice.face_data(0), next_indicatrice.face_data(0));
  face_data_[1].initialise(interfaces, splitting, IJK_Splitting::FACES_J, old_indicatrice.face_data(1), next_indicatrice.face_data(1));
  face_data_[2].initialise(interfaces, splitting, IJK_Splitting::FACES_K, old_indicatrice.face_data(2), next_indicatrice.face_data(2));
}

void ElemFace_Cut_cell_FT_Disc::update(const ElemFace_IJK_Field_double& old_indicatrice, const ElemFace_IJK_Field_double& next_indicatrice)
{
  elem_data_.update(old_indicatrice.elem_data(), next_indicatrice.elem_data());
  face_data_[0].update(old_indicatrice.face_data(0), next_indicatrice.face_data(0));
  face_data_[1].update(old_indicatrice.face_data(1), next_indicatrice.face_data(1));
  face_data_[2].update(old_indicatrice.face_data(2), next_indicatrice.face_data(2));
}

void ElemFace_Cut_cell_FT_Disc::remove_dead_and_virtual_cells(const ElemFace_IJK_Field_double& next_indicatrice)
{
  elem_data_.remove_dead_and_virtual_cells(next_indicatrice.elem_data());
  face_data_[0].remove_dead_and_virtual_cells(next_indicatrice.face_data(0));
  face_data_[1].remove_dead_and_virtual_cells(next_indicatrice.face_data(1));
  face_data_[2].remove_dead_and_virtual_cells(next_indicatrice.face_data(2));
}

void ElemFace_IJK_Field_double::set_to_uniform_value(double valeur)
{
  elem_data_->data() = valeur;
  face_data_[0]->data() = valeur;
  face_data_[1]->data() = valeur;
  face_data_[2]->data() = valeur;
}

void ElemFace_IJK_Field_vector3_double::set_to_uniform_value(double valeur)
{
  elem_data_[0].data() = valeur;
  elem_data_[1].data() = valeur;
  elem_data_[2].data() = valeur;
  face_data_[0][0].data() = valeur;
  face_data_[0][1].data() = valeur;
  face_data_[0][2].data() = valeur;
  face_data_[1][0].data() = valeur;
  face_data_[1][1].data() = valeur;
  face_data_[1][2].data() = valeur;
  face_data_[2][0].data() = valeur;
  face_data_[2][1].data() = valeur;
  face_data_[2][2].data() = valeur;
}

void ElemFace_IJK_Field_double::echange_espace_virtuel(int ghost)
{
  elem_data_->echange_espace_virtuel(ghost);
  face_data_[0]->echange_espace_virtuel(ghost);
  face_data_[1]->echange_espace_virtuel(ghost);
  face_data_[2]->echange_espace_virtuel(ghost);
}

void ElemFace_IJK_Field_vector3_double::echange_espace_virtuel(int ghost)
{
  elem_data_[0].echange_espace_virtuel(ghost);
  elem_data_[1].echange_espace_virtuel(ghost);
  elem_data_[2].echange_espace_virtuel(ghost);
  face_data_[0][0].echange_espace_virtuel(ghost);
  face_data_[0][1].echange_espace_virtuel(ghost);
  face_data_[0][2].echange_espace_virtuel(ghost);
  face_data_[1][0].echange_espace_virtuel(ghost);
  face_data_[1][1].echange_espace_virtuel(ghost);
  face_data_[1][2].echange_espace_virtuel(ghost);
  face_data_[2][0].echange_espace_virtuel(ghost);
  face_data_[2][1].echange_espace_virtuel(ghost);
  face_data_[2][2].echange_espace_virtuel(ghost);
}

void ElemFace_IJK_Field_vector3_double::echange_espace_virtuel()
{
  for (int i = 0; i < 3; i++)
    {
      elem_data_[i].echange_espace_virtuel(elem_data_[i].ghost());
    }

  for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
        {
          face_data_[i][j].echange_espace_virtuel(face_data_[i][j].ghost());
        }
    }
}

void ElemFace_DoubleTabFT_cut_cell_scalar::echange_espace_virtuel()
{
  elem_data_.echange_espace_virtuel();
  face_data_[0].echange_espace_virtuel();
  face_data_[1].echange_espace_virtuel();
  face_data_[2].echange_espace_virtuel();
}

void ElemFace_DoubleTabFT_cut_cell_vector3::echange_espace_virtuel()
{
  elem_data_.echange_espace_virtuel();
  face_data_[0].echange_espace_virtuel();
  face_data_[1].echange_espace_virtuel();
  face_data_[2].echange_espace_virtuel();
}

void ElemFace_DoubleTabFT_cut_cell_vector6::echange_espace_virtuel()
{
  elem_data_.echange_espace_virtuel();
  face_data_[0].echange_espace_virtuel();
  face_data_[1].echange_espace_virtuel();
  face_data_[2].echange_espace_virtuel();
}

void ElemFace_IntTabFT_cut_cell_scalar::echange_espace_virtuel()
{
  elem_data_.echange_espace_virtuel();
  face_data_[0].echange_espace_virtuel();
  face_data_[1].echange_espace_virtuel();
  face_data_[2].echange_espace_virtuel();
}

void ElemFace_IntTabFT_cut_cell_vector3::echange_espace_virtuel()
{
  elem_data_.echange_espace_virtuel();
  face_data_[0].echange_espace_virtuel();
  face_data_[1].echange_espace_virtuel();
  face_data_[2].echange_espace_virtuel();
}

void ElemFace_IntTabFT_cut_cell_vector6::echange_espace_virtuel()
{
  elem_data_.echange_espace_virtuel();
  face_data_[0].echange_espace_virtuel();
  face_data_[1].echange_espace_virtuel();
  face_data_[2].echange_espace_virtuel();
}

int ElemFace_IJK_Field_double::ghost()
{
  assert(face_data_[0]->ghost() == elem_data_->ghost());
  assert(face_data_[1]->ghost() == elem_data_->ghost());
  assert(face_data_[2]->ghost() == elem_data_->ghost());
  return elem_data_->ghost();
}

int ElemFace_IJK_Field_vector3_double::ghost()
{
  assert(elem_data_[1].ghost() == elem_data_[0].ghost());
  assert(elem_data_[2].ghost() == elem_data_[0].ghost());
  assert(face_data_[0][0].ghost() == elem_data_[0].ghost());
  assert(face_data_[0][1].ghost() == elem_data_[0].ghost());
  assert(face_data_[0][2].ghost() == elem_data_[0].ghost());
  assert(face_data_[1][0].ghost() == elem_data_[0].ghost());
  assert(face_data_[1][1].ghost() == elem_data_[0].ghost());
  assert(face_data_[1][2].ghost() == elem_data_[0].ghost());
  assert(face_data_[2][0].ghost() == elem_data_[0].ghost());
  assert(face_data_[2][1].ghost() == elem_data_[0].ghost());
  assert(face_data_[2][2].ghost() == elem_data_[0].ghost());
  return elem_data_[0].ghost();
}

void ElemFace_IJK_Field_double::allocate(const IJK_Splitting& splitting, int ghost_size)
{
  elem_data_ = std::make_shared<IJK_Field_double>();
  face_data_[0] = std::make_shared<IJK_Field_double>();
  face_data_[1] = std::make_shared<IJK_Field_double>();
  face_data_[2] = std::make_shared<IJK_Field_double>();

  elem_data_->allocate(splitting, IJK_Splitting::ELEM, ghost_size);
  face_data_[0]->allocate(splitting, IJK_Splitting::FACES_I, ghost_size);
  face_data_[1]->allocate(splitting, IJK_Splitting::FACES_J, ghost_size);
  face_data_[2]->allocate(splitting, IJK_Splitting::FACES_K, ghost_size);
}

void ElemFace_DoubleTabFT_cut_cell_scalar::associer_persistant(ElemFace_Cut_cell_FT_Disc& cut_cell_disc)
{
  elem_data_.associer_persistant(cut_cell_disc.elem_data());
  face_data_[0].associer_persistant(cut_cell_disc.face_data(0));
  face_data_[1].associer_persistant(cut_cell_disc.face_data(1));
  face_data_[2].associer_persistant(cut_cell_disc.face_data(2));
}

void ElemFace_DoubleTabFT_cut_cell_vector3::associer_persistant(ElemFace_Cut_cell_FT_Disc& cut_cell_disc)
{
  elem_data_.associer_persistant(cut_cell_disc.elem_data());
  face_data_[0].associer_persistant(cut_cell_disc.face_data(0));
  face_data_[1].associer_persistant(cut_cell_disc.face_data(1));
  face_data_[2].associer_persistant(cut_cell_disc.face_data(2));
}

void ElemFace_DoubleTabFT_cut_cell_vector6::associer_persistant(ElemFace_Cut_cell_FT_Disc& cut_cell_disc)
{
  elem_data_.associer_persistant(cut_cell_disc.elem_data());
  face_data_[0].associer_persistant(cut_cell_disc.face_data(0));
  face_data_[1].associer_persistant(cut_cell_disc.face_data(1));
  face_data_[2].associer_persistant(cut_cell_disc.face_data(2));
}

void ElemFace_DoubleTabFT_cut_cell_scalar::associer_ephemere(ElemFace_Cut_cell_FT_Disc& cut_cell_disc)
{
  elem_data_.associer_ephemere(cut_cell_disc.elem_data());
  face_data_[0].associer_ephemere(cut_cell_disc.face_data(0));
  face_data_[1].associer_ephemere(cut_cell_disc.face_data(1));
  face_data_[2].associer_ephemere(cut_cell_disc.face_data(2));
}

void ElemFace_DoubleTabFT_cut_cell_vector3::associer_ephemere(ElemFace_Cut_cell_FT_Disc& cut_cell_disc)
{
  elem_data_.associer_ephemere(cut_cell_disc.elem_data());
  face_data_[0].associer_ephemere(cut_cell_disc.face_data(0));
  face_data_[1].associer_ephemere(cut_cell_disc.face_data(1));
  face_data_[2].associer_ephemere(cut_cell_disc.face_data(2));
}

void ElemFace_DoubleTabFT_cut_cell_vector6::associer_ephemere(ElemFace_Cut_cell_FT_Disc& cut_cell_disc)
{
  elem_data_.associer_ephemere(cut_cell_disc.elem_data());
  face_data_[0].associer_ephemere(cut_cell_disc.face_data(0));
  face_data_[1].associer_ephemere(cut_cell_disc.face_data(1));
  face_data_[2].associer_ephemere(cut_cell_disc.face_data(2));
}

void ElemFace_DoubleTabFT_cut_cell_scalar::associer_paresseux(ElemFace_Cut_cell_FT_Disc& cut_cell_disc)
{
  elem_data_.associer_paresseux(cut_cell_disc.elem_data());
  face_data_[0].associer_paresseux(cut_cell_disc.face_data(0));
  face_data_[1].associer_paresseux(cut_cell_disc.face_data(1));
  face_data_[2].associer_paresseux(cut_cell_disc.face_data(2));
}

void ElemFace_DoubleTabFT_cut_cell_vector3::associer_paresseux(ElemFace_Cut_cell_FT_Disc& cut_cell_disc)
{
  elem_data_.associer_paresseux(cut_cell_disc.elem_data());
  face_data_[0].associer_paresseux(cut_cell_disc.face_data(0));
  face_data_[1].associer_paresseux(cut_cell_disc.face_data(1));
  face_data_[2].associer_paresseux(cut_cell_disc.face_data(2));
}

void ElemFace_DoubleTabFT_cut_cell_vector6::associer_paresseux(ElemFace_Cut_cell_FT_Disc& cut_cell_disc)
{
  elem_data_.associer_paresseux(cut_cell_disc.elem_data());
  face_data_[0].associer_paresseux(cut_cell_disc.face_data(0));
  face_data_[1].associer_paresseux(cut_cell_disc.face_data(1));
  face_data_[2].associer_paresseux(cut_cell_disc.face_data(2));
}

void ElemFace_IntTabFT_cut_cell_scalar::associer_persistant(ElemFace_Cut_cell_FT_Disc& cut_cell_disc)
{
  elem_data_.associer_persistant(cut_cell_disc.elem_data());
  face_data_[0].associer_persistant(cut_cell_disc.face_data(0));
  face_data_[1].associer_persistant(cut_cell_disc.face_data(1));
  face_data_[2].associer_persistant(cut_cell_disc.face_data(2));
}

void ElemFace_IntTabFT_cut_cell_vector3::associer_persistant(ElemFace_Cut_cell_FT_Disc& cut_cell_disc)
{
  elem_data_.associer_persistant(cut_cell_disc.elem_data());
  face_data_[0].associer_persistant(cut_cell_disc.face_data(0));
  face_data_[1].associer_persistant(cut_cell_disc.face_data(1));
  face_data_[2].associer_persistant(cut_cell_disc.face_data(2));
}

void ElemFace_IntTabFT_cut_cell_vector6::associer_persistant(ElemFace_Cut_cell_FT_Disc& cut_cell_disc)
{
  elem_data_.associer_persistant(cut_cell_disc.elem_data());
  face_data_[0].associer_persistant(cut_cell_disc.face_data(0));
  face_data_[1].associer_persistant(cut_cell_disc.face_data(1));
  face_data_[2].associer_persistant(cut_cell_disc.face_data(2));
}

void ElemFace_IntTabFT_cut_cell_scalar::associer_ephemere(ElemFace_Cut_cell_FT_Disc& cut_cell_disc)
{
  elem_data_.associer_ephemere(cut_cell_disc.elem_data());
  face_data_[0].associer_ephemere(cut_cell_disc.face_data(0));
  face_data_[1].associer_ephemere(cut_cell_disc.face_data(1));
  face_data_[2].associer_ephemere(cut_cell_disc.face_data(2));
}

void ElemFace_IntTabFT_cut_cell_vector3::associer_ephemere(ElemFace_Cut_cell_FT_Disc& cut_cell_disc)
{
  elem_data_.associer_ephemere(cut_cell_disc.elem_data());
  face_data_[0].associer_ephemere(cut_cell_disc.face_data(0));
  face_data_[1].associer_ephemere(cut_cell_disc.face_data(1));
  face_data_[2].associer_ephemere(cut_cell_disc.face_data(2));
}

void ElemFace_IntTabFT_cut_cell_vector6::associer_ephemere(ElemFace_Cut_cell_FT_Disc& cut_cell_disc)
{
  elem_data_.associer_ephemere(cut_cell_disc.elem_data());
  face_data_[0].associer_ephemere(cut_cell_disc.face_data(0));
  face_data_[1].associer_ephemere(cut_cell_disc.face_data(1));
  face_data_[2].associer_ephemere(cut_cell_disc.face_data(2));
}

void ElemFace_IntTabFT_cut_cell_scalar::associer_paresseux(ElemFace_Cut_cell_FT_Disc& cut_cell_disc)
{
  elem_data_.associer_paresseux(cut_cell_disc.elem_data());
  face_data_[0].associer_paresseux(cut_cell_disc.face_data(0));
  face_data_[1].associer_paresseux(cut_cell_disc.face_data(1));
  face_data_[2].associer_paresseux(cut_cell_disc.face_data(2));
}

void ElemFace_IntTabFT_cut_cell_vector3::associer_paresseux(ElemFace_Cut_cell_FT_Disc& cut_cell_disc)
{
  elem_data_.associer_paresseux(cut_cell_disc.elem_data());
  face_data_[0].associer_paresseux(cut_cell_disc.face_data(0));
  face_data_[1].associer_paresseux(cut_cell_disc.face_data(1));
  face_data_[2].associer_paresseux(cut_cell_disc.face_data(2));
}

void ElemFace_IntTabFT_cut_cell_vector6::associer_paresseux(ElemFace_Cut_cell_FT_Disc& cut_cell_disc)
{
  elem_data_.associer_paresseux(cut_cell_disc.elem_data());
  face_data_[0].associer_paresseux(cut_cell_disc.face_data(0));
  face_data_[1].associer_paresseux(cut_cell_disc.face_data(1));
  face_data_[2].associer_paresseux(cut_cell_disc.face_data(2));
}

void allocate_elem_face_velocity(ElemFace_IJK_Field_vector3_double& v, const IJK_Splitting& s, int ghost, double DU)
{
  v.elem_data().get_ptr(0) = std::make_shared<IJK_Field_double>();
  v.elem_data().get_ptr(1) = std::make_shared<IJK_Field_double>();
  v.elem_data().get_ptr(2) = std::make_shared<IJK_Field_double>();
  v.elem_data()[0].allocate(s, IJK_Splitting::FACES_I, ghost);
  v.elem_data()[1].allocate(s, IJK_Splitting::FACES_J, ghost);
  v.elem_data()[2].allocate(s, IJK_Splitting::FACES_K, ghost);
  v.elem_data()[0].get_shear_BC_helpler().set_dU_(DU);
  v.elem_data()[1].get_shear_BC_helpler().set_dU_(0.);
  v.elem_data()[2].get_shear_BC_helpler().set_dU_(0.);

  // Note : la vrai localisation est aux aretes, ce qui n'existe actuellement pas dans IJK_Splitting.
  v.face_data(0).get_ptr(0) = std::make_shared<IJK_Field_double>();
  v.face_data(0).get_ptr(1) = std::make_shared<IJK_Field_double>();
  v.face_data(0).get_ptr(2) = std::make_shared<IJK_Field_double>();
  v.face_data(1).get_ptr(0) = std::make_shared<IJK_Field_double>();
  v.face_data(1).get_ptr(1) = std::make_shared<IJK_Field_double>();
  v.face_data(1).get_ptr(2) = std::make_shared<IJK_Field_double>();
  v.face_data(2).get_ptr(0) = std::make_shared<IJK_Field_double>();
  v.face_data(2).get_ptr(1) = std::make_shared<IJK_Field_double>();
  v.face_data(2).get_ptr(2) = std::make_shared<IJK_Field_double>();
  v.face_data(0)[0].allocate(s, IJK_Splitting::FACES_I, ghost);
  v.face_data(0)[1].allocate(s, IJK_Splitting::FACES_I, ghost);
  v.face_data(0)[2].allocate(s, IJK_Splitting::FACES_I, ghost);
  v.face_data(1)[0].allocate(s, IJK_Splitting::FACES_J, ghost);
  v.face_data(1)[1].allocate(s, IJK_Splitting::FACES_J, ghost);
  v.face_data(1)[2].allocate(s, IJK_Splitting::FACES_J, ghost);
  v.face_data(2)[0].allocate(s, IJK_Splitting::FACES_K, ghost);
  v.face_data(2)[1].allocate(s, IJK_Splitting::FACES_K, ghost);
  v.face_data(2)[2].allocate(s, IJK_Splitting::FACES_K, ghost);
  v.face_data(0)[0].get_shear_BC_helpler().set_dU_(0.);
  v.face_data(0)[1].get_shear_BC_helpler().set_dU_(0.);
  v.face_data(0)[2].get_shear_BC_helpler().set_dU_(0.);
  v.face_data(1)[0].get_shear_BC_helpler().set_dU_(0.);
  v.face_data(1)[1].get_shear_BC_helpler().set_dU_(0.);
  v.face_data(1)[2].get_shear_BC_helpler().set_dU_(0.);
  v.face_data(2)[0].get_shear_BC_helpler().set_dU_(0.);
  v.face_data(2)[1].get_shear_BC_helpler().set_dU_(0.);
  v.face_data(2)[2].get_shear_BC_helpler().set_dU_(0.);
}

void allocate_elem_face_velocity(ElemFace_Cut_field_vector3_double& v, const IJK_Splitting& s, int ghost, double DU)
{
  v.elem_data().get_ptr(0) = std::make_shared<Cut_field_double>();
  v.elem_data().get_ptr(1) = std::make_shared<Cut_field_double>();
  v.elem_data().get_ptr(2) = std::make_shared<Cut_field_double>();
  v.elem_data()[0].allocate(s, IJK_Splitting::FACES_I, ghost);
  v.elem_data()[1].allocate(s, IJK_Splitting::FACES_J, ghost);
  v.elem_data()[2].allocate(s, IJK_Splitting::FACES_K, ghost);
  v.elem_data()[0].get_shear_BC_helpler().set_dU_(DU);
  v.elem_data()[1].get_shear_BC_helpler().set_dU_(0.);
  v.elem_data()[2].get_shear_BC_helpler().set_dU_(0.);

  // Note : la vrai localisation est aux aretes, ce qui n'existe actuellement pas dans IJK_Splitting.
  v.face_data(0).get_ptr(0) = std::make_shared<Cut_field_double>();
  v.face_data(0).get_ptr(1) = std::make_shared<Cut_field_double>();
  v.face_data(0).get_ptr(2) = std::make_shared<Cut_field_double>();
  v.face_data(1).get_ptr(0) = std::make_shared<Cut_field_double>();
  v.face_data(1).get_ptr(1) = std::make_shared<Cut_field_double>();
  v.face_data(1).get_ptr(2) = std::make_shared<Cut_field_double>();
  v.face_data(2).get_ptr(0) = std::make_shared<Cut_field_double>();
  v.face_data(2).get_ptr(1) = std::make_shared<Cut_field_double>();
  v.face_data(2).get_ptr(2) = std::make_shared<Cut_field_double>();
  v.face_data(0)[0].allocate(s, IJK_Splitting::FACES_I, ghost);
  v.face_data(0)[1].allocate(s, IJK_Splitting::FACES_I, ghost);
  v.face_data(0)[2].allocate(s, IJK_Splitting::FACES_I, ghost);
  v.face_data(1)[0].allocate(s, IJK_Splitting::FACES_J, ghost);
  v.face_data(1)[1].allocate(s, IJK_Splitting::FACES_J, ghost);
  v.face_data(1)[2].allocate(s, IJK_Splitting::FACES_J, ghost);
  v.face_data(2)[0].allocate(s, IJK_Splitting::FACES_K, ghost);
  v.face_data(2)[1].allocate(s, IJK_Splitting::FACES_K, ghost);
  v.face_data(2)[2].allocate(s, IJK_Splitting::FACES_K, ghost);
  v.face_data(0)[0].get_shear_BC_helpler().set_dU_(0.);
  v.face_data(0)[1].get_shear_BC_helpler().set_dU_(0.);
  v.face_data(0)[2].get_shear_BC_helpler().set_dU_(0.);
  v.face_data(1)[0].get_shear_BC_helpler().set_dU_(0.);
  v.face_data(1)[1].get_shear_BC_helpler().set_dU_(0.);
  v.face_data(1)[2].get_shear_BC_helpler().set_dU_(0.);
  v.face_data(2)[0].get_shear_BC_helpler().set_dU_(0.);
  v.face_data(2)[1].get_shear_BC_helpler().set_dU_(0.);
  v.face_data(2)[2].get_shear_BC_helpler().set_dU_(0.);
}

void allocate_elem_face_velocity_persistant(ElemFace_Cut_cell_FT_Disc& cut_cell_disc, ElemFace_Cut_field_vector3_double& v, const IJK_Splitting& s, int ghost, double DU)
{
  allocate_elem_face_velocity(v, s, ghost, DU);

  v.elem_data()[0].associer_persistant(cut_cell_disc.elem_data());
  v.elem_data()[1].associer_persistant(cut_cell_disc.elem_data());
  v.elem_data()[2].associer_persistant(cut_cell_disc.elem_data());
  v.face_data(0)[0].associer_persistant(cut_cell_disc.face_data(0));
  v.face_data(0)[1].associer_persistant(cut_cell_disc.face_data(0));
  v.face_data(0)[2].associer_persistant(cut_cell_disc.face_data(0));
  v.face_data(1)[0].associer_persistant(cut_cell_disc.face_data(1));
  v.face_data(1)[1].associer_persistant(cut_cell_disc.face_data(1));
  v.face_data(1)[2].associer_persistant(cut_cell_disc.face_data(1));
  v.face_data(2)[0].associer_persistant(cut_cell_disc.face_data(2));
  v.face_data(2)[1].associer_persistant(cut_cell_disc.face_data(2));
  v.face_data(2)[2].associer_persistant(cut_cell_disc.face_data(2));
}

void allocate_elem_face_velocity_ephemere(ElemFace_Cut_cell_FT_Disc& cut_cell_disc, ElemFace_Cut_field_vector3_double& v, const IJK_Splitting& s, int ghost, double DU)
{
  allocate_elem_face_velocity(v, s, ghost, DU);

  v.elem_data()[0].associer_ephemere(cut_cell_disc.elem_data());
  v.elem_data()[1].associer_ephemere(cut_cell_disc.elem_data());
  v.elem_data()[2].associer_ephemere(cut_cell_disc.elem_data());
  v.face_data(0)[0].associer_ephemere(cut_cell_disc.face_data(0));
  v.face_data(0)[1].associer_ephemere(cut_cell_disc.face_data(0));
  v.face_data(0)[2].associer_ephemere(cut_cell_disc.face_data(0));
  v.face_data(1)[0].associer_ephemere(cut_cell_disc.face_data(1));
  v.face_data(1)[1].associer_ephemere(cut_cell_disc.face_data(1));
  v.face_data(1)[2].associer_ephemere(cut_cell_disc.face_data(1));
  v.face_data(2)[0].associer_ephemere(cut_cell_disc.face_data(2));
  v.face_data(2)[1].associer_ephemere(cut_cell_disc.face_data(2));
  v.face_data(2)[2].associer_ephemere(cut_cell_disc.face_data(2));
}

void allocate_elem_face_velocity_paresseux(ElemFace_Cut_cell_FT_Disc& cut_cell_disc, ElemFace_Cut_field_vector3_double& v, const IJK_Splitting& s, int ghost, double DU)
{
  allocate_elem_face_velocity(v, s, ghost, DU);

  v.elem_data()[0].associer_paresseux(cut_cell_disc.elem_data());
  v.elem_data()[1].associer_paresseux(cut_cell_disc.elem_data());
  v.elem_data()[2].associer_paresseux(cut_cell_disc.elem_data());
  v.face_data(0)[0].associer_paresseux(cut_cell_disc.face_data(0));
  v.face_data(0)[1].associer_paresseux(cut_cell_disc.face_data(0));
  v.face_data(0)[2].associer_paresseux(cut_cell_disc.face_data(0));
  v.face_data(1)[0].associer_paresseux(cut_cell_disc.face_data(1));
  v.face_data(1)[1].associer_paresseux(cut_cell_disc.face_data(1));
  v.face_data(1)[2].associer_paresseux(cut_cell_disc.face_data(1));
  v.face_data(2)[0].associer_paresseux(cut_cell_disc.face_data(2));
  v.face_data(2)[1].associer_paresseux(cut_cell_disc.face_data(2));
  v.face_data(2)[2].associer_paresseux(cut_cell_disc.face_data(2));
}

void allocate_elem_face_cell_vector(ElemFace_IJK_Field_vector3_double& v, const IJK_Splitting& s, int ghost)
{
  v.elem_data().get_ptr(0) = std::make_shared<IJK_Field_double>();
  v.elem_data().get_ptr(1) = std::make_shared<IJK_Field_double>();
  v.elem_data().get_ptr(2) = std::make_shared<IJK_Field_double>();
  v.elem_data()[0].allocate(s, IJK_Splitting::ELEM, ghost);
  v.elem_data()[1].allocate(s, IJK_Splitting::ELEM, ghost);
  v.elem_data()[2].allocate(s, IJK_Splitting::ELEM, ghost);
  v.elem_data()[0].get_shear_BC_helpler().set_dU_(0.);
  v.elem_data()[1].get_shear_BC_helpler().set_dU_(0.);
  v.elem_data()[2].get_shear_BC_helpler().set_dU_(0.);

  for (int i=0; i<3 ; i++)
    {
      v.face_data(0).get_ptr(i) = std::make_shared<IJK_Field_double>();
      v.face_data(1).get_ptr(i) = std::make_shared<IJK_Field_double>();
      v.face_data(2).get_ptr(i) = std::make_shared<IJK_Field_double>();
      v.face_data(0)[i].allocate(s, IJK_Splitting::FACES_I, ghost);
      v.face_data(1)[i].allocate(s, IJK_Splitting::FACES_J, ghost);
      v.face_data(2)[i].allocate(s, IJK_Splitting::FACES_K, ghost);
      v.face_data(0)[i].get_shear_BC_helpler().set_dU_(0.);
      v.face_data(1)[i].get_shear_BC_helpler().set_dU_(0.);
      v.face_data(2)[i].get_shear_BC_helpler().set_dU_(0.);
    }
}

// Extension d'un champ elem vers les faces, imposant a la face la moyenne des valeurs non-nulles des elements contenant la face
void extend_elem_to_face(const IJK_Field_double& field_elem, IJK_Field_double& field_face_x, IJK_Field_double& field_face_y, IJK_Field_double& field_face_z)
{
  const int ni = field_elem.ni();
  const int nj = field_elem.nj();
  const int nk = field_elem.nk();
  const int ghost = field_elem.ghost();
  assert(field_elem.ghost() == field_face_x.ghost());
  assert(field_elem.ghost() == field_face_y.ghost());
  assert(field_elem.ghost() == field_face_z.ghost());
  assert(nk == field_face_x.nk());
  assert(nk == field_face_y.nk());
  assert(nk == field_face_z.nk());
  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              field_face_x(i,j,k) = (field_elem(i,j,k) == 0) ? field_elem(i-1,j,k) : ((field_elem(i-1,j,k) == 0) ? field_elem(i,j,k) : (.5*(field_elem(i,j,k) + field_elem(i-1,j,k))));
              field_face_y(i,j,k) = (field_elem(i,j,k) == 0) ? field_elem(i,j-1,k) : ((field_elem(i,j-1,k) == 0) ? field_elem(i,j,k) : (.5*(field_elem(i,j,k) + field_elem(i,j-1,k))));
              field_face_z(i,j,k) = (field_elem(i,j,k) == 0) ? field_elem(i,j,k-1) : ((field_elem(i,j,k-1) == 0) ? field_elem(i,j,k) : (.5*(field_elem(i,j,k) + field_elem(i,j,k-1))));
            }
        }
    }

  field_face_x.echange_espace_virtuel(ghost);
  field_face_y.echange_espace_virtuel(ghost);
  field_face_z.echange_espace_virtuel(ghost);
}

// Extension d'un champ elem vers les faces, imposant a la face la moyenne des valeurs non-nulles des elements contenant la face
void extend_elem_to_face(const IJK_Field_vector3_double& field_elem, IJK_Field_vector3_double& field_face_x, IJK_Field_vector3_double& field_face_y, IJK_Field_vector3_double& field_face_z)
{
  for (int dir = 0; dir < 3; dir++)
    {
      const int ni = field_elem[dir].ni();
      const int nj = field_elem[dir].nj();
      const int nk = field_elem[dir].nk();
      const int ghost = field_elem[dir].ghost();
      assert(field_elem[dir].ghost() == field_face_x[dir].ghost());
      assert(field_elem[dir].ghost() == field_face_y[dir].ghost());
      assert(field_elem[dir].ghost() == field_face_z[dir].ghost());
      assert(nk == field_face_x[dir].nk());
      assert(nk == field_face_y[dir].nk());
      assert(nk == field_face_z[dir].nk());
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  field_face_x[dir](i,j,k) = (field_elem[dir](i,j,k) == 0) ? field_elem[dir](i-1,j,k) : ((field_elem[dir](i-1,j,k) == 0) ? field_elem[dir](i,j,k) : (.5*(field_elem[dir](i,j,k) + field_elem[dir](i-1,j,k))));
                  field_face_y[dir](i,j,k) = (field_elem[dir](i,j,k) == 0) ? field_elem[dir](i,j-1,k) : ((field_elem[dir](i,j-1,k) == 0) ? field_elem[dir](i,j,k) : (.5*(field_elem[dir](i,j,k) + field_elem[dir](i,j-1,k))));
                  field_face_z[dir](i,j,k) = (field_elem[dir](i,j,k) == 0) ? field_elem[dir](i,j,k-1) : ((field_elem[dir](i,j,k-1) == 0) ? field_elem[dir](i,j,k) : (.5*(field_elem[dir](i,j,k) + field_elem[dir](i,j,k-1))));
                }
            }
        }

      field_face_x[dir].echange_espace_virtuel(ghost);
      field_face_y[dir].echange_espace_virtuel(ghost);
      field_face_z[dir].echange_espace_virtuel(ghost);
    }
}

