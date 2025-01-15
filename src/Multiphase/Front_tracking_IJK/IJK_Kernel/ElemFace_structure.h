/****************************************************************************
* Copyright (c) 2022, CEA
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
#ifndef ElemFace_structure_included
#define ElemFace_structure_included
#include <assert.h>
#include <FixedVector.h>
#include <IJK_Field_vector.h>
#include <IJK_Field.h>
#include <IJK_Navier_Stokes_tools.h>
#include <Champ_diphasique.h>
#include <Cut_cell_FT_Disc.h>

class ElemFace_IJK_Field_double;
class IJK_Interfaces;

/*! @brief : class ElemFace_structure<T>
 *
 *  Cette classe contient un champ sur le volume de controle des elements
 *  et trois champs sur les volumes de controle centres sur les face.
 *
 *  Pour l'utilisation sur des champs de donnees, des classes derivees sont
 *  disponibles, incluant des routines specifiques permettant d'operer
 *  certaines actions sur l'ensemble des champs.
 *
 *
 */
template<class T>
class ElemFace_structure
{
public:
  ElemFace_structure() { }

  T elem_data_;
  FixedVector<T, 3> face_data_;

  inline T& elem_data()
  {
    return this->elem_data_;
  }

  inline const T& elem_data() const
  {
    return this->elem_data_;
  }

  inline T& face_data(int dir)
  {
    return this->face_data_[dir];
  }

  inline const T& face_data(int dir) const
  {
    return this->face_data_[dir];
  }

  inline T& data_loc(IJK_Splitting::Localisation loc)
  {
    if (loc == IJK_Splitting::ELEM)
      {
        return this->elem_data_;
      }
    else if (loc == IJK_Splitting::FACES_I)
      {
        return this->face_data_[0];
      }
    else if (loc == IJK_Splitting::FACES_J)
      {
        return this->face_data_[1];
      }
    else if (loc == IJK_Splitting::FACES_K)
      {
        return this->face_data_[2];
      }
    else
      {
        Cerr << "Invalid localisation in ElemFace_structure<T>::data_loc." << finl;
        Process::exit();
        return this->elem_data_;
      }
  }

  inline const T& data_loc(IJK_Splitting::Localisation loc) const
  {
    if (loc == IJK_Splitting::ELEM)
      {
        return this->elem_data_;
      }
    else if (loc == IJK_Splitting::FACES_I)
      {
        return this->face_data_[0];
      }
    else if (loc == IJK_Splitting::FACES_J)
      {
        return this->face_data_[1];
      }
    else if (loc == IJK_Splitting::FACES_K)
      {
        return this->face_data_[2];
      }
    else
      {
        Cerr << "Invalid localisation in ElemFace_structure<T>::data_loc." << finl;
        Process::exit();
        return this->elem_data_;
      }
  }

protected:
};

/*! @brief : class ElemFace_ptr_structure<T>
 *
 *  En comparaison a ElemFace_structure<T>, cette classe ne contient pas
 *  les champs directement mais un std::shared_ptr<T>,
 *
 *
 */
template<class T>
class ElemFace_ptr_structure
{
public:
  ElemFace_ptr_structure() { }

  std::shared_ptr<T> elem_data_;
  FixedVector<std::shared_ptr<T>, 3> face_data_;

  inline T& elem_data()
  {
    return *this->elem_data_;
  }

  inline const T& elem_data() const
  {
    return *this->elem_data_;
  }

  inline T& face_data(int dir)
  {
    return *this->face_data_[dir];
  }

  inline const T& face_data(int dir) const
  {
    return *this->face_data_[dir];
  }

  inline T& data_loc(IJK_Splitting::Localisation loc)
  {
    if (loc == IJK_Splitting::ELEM)
      {
        return *this->elem_data_;
      }
    else if (loc == IJK_Splitting::FACES_I)
      {
        return *this->face_data_[0];
      }
    else if (loc == IJK_Splitting::FACES_J)
      {
        return *this->face_data_[1];
      }
    else if (loc == IJK_Splitting::FACES_K)
      {
        return *this->face_data_[2];
      }
    else
      {
        Cerr << "Invalid localisation in ElemFace_ptr_structure<T>::data_loc." << finl;
        Process::exit();
        return *this->elem_data_;
      }
  }

  inline const T& data_loc(IJK_Splitting::Localisation loc) const
  {
    if (loc == IJK_Splitting::ELEM)
      {
        return *this->elem_data_;
      }
    else if (loc == IJK_Splitting::FACES_I)
      {
        return *this->face_data_[0];
      }
    else if (loc == IJK_Splitting::FACES_J)
      {
        return *this->face_data_[1];
      }
    else if (loc == IJK_Splitting::FACES_K)
      {
        return *this->face_data_[2];
      }
    else
      {
        Cerr << "Invalid localisation in ElemFace_ptr_structure<T>::data_loc." << finl;
        Process::exit();
        return *this->elem_data_;
      }
  }

protected:
};

/*! @brief : class ElemFace_Cut_cell_FT_Disc
 *
 *  <Description of class ElemFace_Cut_cell_FT_Disc>
 *
 *
 */
class ElemFace_Cut_cell_FT_Disc : public ElemFace_structure<Cut_cell_FT_Disc>
{
public:
  ElemFace_Cut_cell_FT_Disc() { }

  void initialise(IJK_Interfaces& interfaces, IJK_Splitting& splitting);
  void initialise(IJK_Interfaces& interfaces, IJK_Splitting& splitting, const ElemFace_IJK_Field_double& old_indicatrice, const ElemFace_IJK_Field_double& next_indicatrice);
  void update(const ElemFace_IJK_Field_double& old_indicatrice, const ElemFace_IJK_Field_double& next_indicatrice);
  void remove_dead_and_virtual_cells(const ElemFace_IJK_Field_double& next_indicatrice);

protected:
};

/*! @brief : class ElemFace_IJK_Field_double
 *
 *  Note : on utilise ElemFace_ptr_structure pour permettre le polymorphisme
 *  entre IJK_Field_double et Cut_field_double.
 *
 */
class ElemFace_IJK_Field_double : public ElemFace_ptr_structure<IJK_Field_double>
{
public:
  ElemFace_IJK_Field_double() { }

  void set_to_uniform_value(double valeur);
  void echange_espace_virtuel(int ghost);
  int ghost();

  void allocate(const IJK_Splitting&, int ghost_size);

protected:
};

/*! @brief : class ElemFace_IJK_Field_vector3_double
 *
 *  <Description of class ElemFace_IJK_Field_vector3_double>
 *
 *
 */
class ElemFace_IJK_Field_vector3_double : public ElemFace_structure<IJK_Field_vector3_double>
{
public:
  ElemFace_IJK_Field_vector3_double() { }

  void set_to_uniform_value(double valeur);
  void echange_espace_virtuel(int ghost);
  void echange_espace_virtuel();
  int ghost();

protected:
};

/*! @brief : class ElemFace_DoubleTabFT_cut_cell_scalar
 *
 *  <Description of class ElemFace_DoubleTabFT_cut_cell_scalar>
 *
 *
 */
class ElemFace_DoubleTabFT_cut_cell_scalar : public ElemFace_structure<DoubleTabFT_cut_cell_scalar>
{
public:
  ElemFace_DoubleTabFT_cut_cell_scalar() { }

  void echange_espace_virtuel();

  void associer_persistant(ElemFace_Cut_cell_FT_Disc& cut_cell_disc);
  void associer_ephemere(ElemFace_Cut_cell_FT_Disc& cut_cell_disc);
  void associer_paresseux(ElemFace_Cut_cell_FT_Disc& cut_cell_disc);

protected:
};

/*! @brief : class ElemFace_DoubleTabFT_cut_cell_vector3
 *
 *  <Description of class ElemFace_DoubleTabFT_cut_cell_vector3>
 *
 *
 */
class ElemFace_DoubleTabFT_cut_cell_vector3 : public ElemFace_structure<DoubleTabFT_cut_cell_vector3>
{
public:
  ElemFace_DoubleTabFT_cut_cell_vector3() { }

  void echange_espace_virtuel();

  void associer_persistant(ElemFace_Cut_cell_FT_Disc& cut_cell_disc);
  void associer_ephemere(ElemFace_Cut_cell_FT_Disc& cut_cell_disc);
  void associer_paresseux(ElemFace_Cut_cell_FT_Disc& cut_cell_disc);

protected:
};

/*! @brief : class ElemFace_DoubleTabFT_cut_cell_vector6
 *
 *  <Description of class ElemFace_DoubleTabFT_cut_cell_vector6>
 *
 *
 */
class ElemFace_DoubleTabFT_cut_cell_vector6 : public ElemFace_structure<DoubleTabFT_cut_cell_vector6>
{
public:
  ElemFace_DoubleTabFT_cut_cell_vector6() { }

  void echange_espace_virtuel();

  void associer_persistant(ElemFace_Cut_cell_FT_Disc& cut_cell_disc);
  void associer_ephemere(ElemFace_Cut_cell_FT_Disc& cut_cell_disc);
  void associer_paresseux(ElemFace_Cut_cell_FT_Disc& cut_cell_disc);

protected:
};

/*! @brief : class ElemFace_IntTabFT_cut_cell_scalar
 *
 *  <Description of class ElemFace_IntTabFT_cut_cell_scalar>
 *
 *
 */
class ElemFace_IntTabFT_cut_cell_scalar : public ElemFace_structure<IntTabFT_cut_cell_scalar>
{
public:
  ElemFace_IntTabFT_cut_cell_scalar() { }

  void echange_espace_virtuel();

  void associer_persistant(ElemFace_Cut_cell_FT_Disc& cut_cell_disc);
  void associer_ephemere(ElemFace_Cut_cell_FT_Disc& cut_cell_disc);
  void associer_paresseux(ElemFace_Cut_cell_FT_Disc& cut_cell_disc);

protected:
};

/*! @brief : class ElemFace_IntTabFT_cut_cell_vector3
 *
 *  <Description of class ElemFace_IntTabFT_cut_cell_vector3>
 *
 *
 */
class ElemFace_IntTabFT_cut_cell_vector3 : public ElemFace_structure<IntTabFT_cut_cell_vector3>
{
public:
  ElemFace_IntTabFT_cut_cell_vector3() { }

  void echange_espace_virtuel();

  void associer_persistant(ElemFace_Cut_cell_FT_Disc& cut_cell_disc);
  void associer_ephemere(ElemFace_Cut_cell_FT_Disc& cut_cell_disc);
  void associer_paresseux(ElemFace_Cut_cell_FT_Disc& cut_cell_disc);

protected:
};

/*! @brief : class ElemFace_IntTabFT_cut_cell_vector6
 *
 *  <Description of class ElemFace_IntTabFT_cut_cell_vector6>
 *
 *
 */
class ElemFace_IntTabFT_cut_cell_vector6 : public ElemFace_structure<IntTabFT_cut_cell_vector6>
{
public:
  ElemFace_IntTabFT_cut_cell_vector6() { }

  void echange_espace_virtuel();

  void associer_persistant(ElemFace_Cut_cell_FT_Disc& cut_cell_disc);
  void associer_ephemere(ElemFace_Cut_cell_FT_Disc& cut_cell_disc);
  void associer_paresseux(ElemFace_Cut_cell_FT_Disc& cut_cell_disc);

protected:
};

/*! @brief : class ElemFace_Cut_field_double
 *
 *  <Description of class ElemFace_Cut_field_double>
 *
 *
 *
 */
class ElemFace_Cut_field_double : public ElemFace_IJK_Field_double
{
public :
  ElemFace_Cut_field_double() { }

  inline Cut_field_double& elem_data()
  {
    return static_cast<Cut_field_double&>(*this->elem_data_);
  }

  inline const Cut_field_double& elem_data() const
  {
    return static_cast<const Cut_field_double&>(*this->elem_data_);
  }

  inline Cut_field_double& face_data(int dir)
  {
    return static_cast<Cut_field_double&>(*this->face_data_[dir]);
  }

  inline const Cut_field_double& face_data(int dir) const
  {
    return static_cast<const Cut_field_double&>(*this->face_data_[dir]);
  }

  inline Cut_field_double& data_loc(IJK_Splitting::Localisation loc)
  {
    if (loc == IJK_Splitting::ELEM)
      {
        return static_cast<Cut_field_double&>(*this->elem_data_);
      }
    else if (loc == IJK_Splitting::FACES_I)
      {
        return static_cast<Cut_field_double&>(*this->face_data_[0]);
      }
    else if (loc == IJK_Splitting::FACES_J)
      {
        return static_cast<Cut_field_double&>(*this->face_data_[1]);
      }
    else if (loc == IJK_Splitting::FACES_K)
      {
        return static_cast<Cut_field_double&>(*this->face_data_[2]);
      }
    else
      {
        Cerr << "Invalid localisation in ElemFace_Cut_field_double<Cut_field_double>::data_loc." << finl;
        Process::exit();
        return static_cast<Cut_field_double&>(*this->elem_data_);
      }
  }

  inline const Cut_field_double& data_loc(IJK_Splitting::Localisation loc) const
  {
    if (loc == IJK_Splitting::ELEM)
      {
        return static_cast<const Cut_field_double&>(*this->elem_data_);
      }
    else if (loc == IJK_Splitting::FACES_I)
      {
        return static_cast<const Cut_field_double&>(*this->face_data_[0]);
      }
    else if (loc == IJK_Splitting::FACES_J)
      {
        return static_cast<const Cut_field_double&>(*this->face_data_[1]);
      }
    else if (loc == IJK_Splitting::FACES_K)
      {
        return static_cast<const Cut_field_double&>(*this->face_data_[2]);
      }
    else
      {
        Cerr << "Invalid localisation in ElemFace_Cut_field_double<Cut_field_double>::data_loc." << finl;
        Process::exit();
        return static_cast<const Cut_field_double&>(*this->elem_data_);
      }
  }

protected :
};

/*! @brief : class ElemFace_Cut_field_vector3_double
 *
 *  <Description of class ElemFace_Cut_field_vector3_double>
 *
 *
 *
 */
class ElemFace_Cut_field_vector3_double : public ElemFace_IJK_Field_vector3_double
{
public :
  ElemFace_Cut_field_vector3_double() { }

  inline Cut_field_vector3_double& elem_data()
  {
    return static_cast<Cut_field_vector3_double&>(this->elem_data_);
  }

  inline const Cut_field_vector3_double& elem_data() const
  {
    return static_cast<const Cut_field_vector3_double&>(this->elem_data_);
  }

  inline Cut_field_vector3_double& face_data(int dir)
  {
    return static_cast<Cut_field_vector3_double&>(this->face_data_[dir]);
  }

  inline const Cut_field_vector3_double& face_data(int dir) const
  {
    return static_cast<const Cut_field_vector3_double&>(this->face_data_[dir]);
  }

  inline Cut_field_vector3_double& data_loc(IJK_Splitting::Localisation loc)
  {
    if (loc == IJK_Splitting::ELEM)
      {
        return static_cast<Cut_field_vector3_double&>(this->elem_data_);
      }
    else if (loc == IJK_Splitting::FACES_I)
      {
        return static_cast<Cut_field_vector3_double&>(this->face_data_[0]);
      }
    else if (loc == IJK_Splitting::FACES_J)
      {
        return static_cast<Cut_field_vector3_double&>(this->face_data_[1]);
      }
    else if (loc == IJK_Splitting::FACES_K)
      {
        return static_cast<Cut_field_vector3_double&>(this->face_data_[2]);
      }
    else
      {
        Cerr << "Invalid localisation in ElemFace_Cut_field_vector3_double<Cut_field_vector3_double>::data_loc." << finl;
        Process::exit();
        return static_cast<Cut_field_vector3_double&>(this->elem_data_);
      }
  }

  inline const Cut_field_vector3_double& data_loc(IJK_Splitting::Localisation loc) const
  {
    if (loc == IJK_Splitting::ELEM)
      {
        return static_cast<const Cut_field_vector3_double&>(this->elem_data_);
      }
    else if (loc == IJK_Splitting::FACES_I)
      {
        return static_cast<const Cut_field_vector3_double&>(this->face_data_[0]);
      }
    else if (loc == IJK_Splitting::FACES_J)
      {
        return static_cast<const Cut_field_vector3_double&>(this->face_data_[1]);
      }
    else if (loc == IJK_Splitting::FACES_K)
      {
        return static_cast<const Cut_field_vector3_double&>(this->face_data_[2]);
      }
    else
      {
        Cerr << "Invalid localisation in ElemFace_Cut_field_vector3_double<Cut_field_vector3_double>::data_loc." << finl;
        Process::exit();
        return static_cast<const Cut_field_vector3_double&>(this->elem_data_);
      }
  }

protected :
};


//
// Routines utilitaires associees aux champs ElemFace :
//

void allocate_elem_face_velocity(ElemFace_IJK_Field_vector3_double& v, const IJK_Splitting& s, int ghost, double DU=0.);

void allocate_elem_face_velocity(ElemFace_Cut_field_vector3_double& v, const IJK_Splitting& s, int ghost, double DU=0.);
void allocate_elem_face_velocity_persistant(ElemFace_Cut_cell_FT_Disc& cut_cell_disc, ElemFace_Cut_field_vector3_double& v, const IJK_Splitting& s, int ghost, double DU=0.);
void allocate_elem_face_velocity_ephemere(ElemFace_Cut_cell_FT_Disc& cut_cell_disc,   ElemFace_Cut_field_vector3_double& v, const IJK_Splitting& s, int ghost, double DU=0.);
void allocate_elem_face_velocity_paresseux(ElemFace_Cut_cell_FT_Disc& cut_cell_disc,  ElemFace_Cut_field_vector3_double& v, const IJK_Splitting& s, int ghost, double DU=0.);

void allocate_elem_face_cell_vector(ElemFace_IJK_Field_vector3_double& v, const IJK_Splitting& s, int ghost);

void extend_elem_to_face(const IJK_Field_double& field_elem, IJK_Field_double& field_face_x, IJK_Field_double& field_face_y, IJK_Field_double& field_face_z);
void extend_elem_to_face(const IJK_Field_vector3_double& field_elem, IJK_Field_vector3_double& field_face_x, IJK_Field_vector3_double& field_face_y, IJK_Field_vector3_double& field_face_z);

#endif
