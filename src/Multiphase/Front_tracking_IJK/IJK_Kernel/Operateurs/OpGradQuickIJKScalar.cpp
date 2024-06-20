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
/////////////////////////////////////////////////////////////////////////////
//
// File      : OpGradQuickIJKScalar.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/IJK_Kernel/Operateurs
//
/////////////////////////////////////////////////////////////////////////////

#include <OpGradQuickIJKScalar.h>

Implemente_instanciable_sans_constructeur( OpGradQuickIJKScalar_double, "OpGradQuickIJKScalar_double", OpConvQuickIJKScalar_double ) ;

Sortie& OpGradQuickIJKScalar_double::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  return os;
}

Entree& OpGradQuickIJKScalar_double::readOn( Entree& is )
{
  Objet_U::readOn( is );
  return is;
}

void OpGradQuickIJKScalar_double::calculer_grad(const IJK_Field_double& field,
                                                FixedVector<IJK_Field_double, 3>& result)
{
  input_field_ = &field;
  const int ni = field.ni();
  const int nj = field.nj();
  tmp_curv_fram_.allocate(ni, nj, 4, 1);
  compute_grad(result);
  input_field_ = nullptr;
  input_velocity_x_ = nullptr;
  input_velocity_y_ = nullptr;
  input_velocity_z_ = nullptr;
}

void OpGradQuickIJKScalar_double::calculer_grad_x(const IJK_Field_double& field,
                                                  IJK_Field_double& result)
{
  input_field_ = &field;
  const int ni = field.ni();
  const int nj = field.nj();
  tmp_curv_fram_.allocate(ni, nj, 4, 1);
  compute_grad_x(result);
  input_field_ = nullptr;
  input_velocity_x_ = nullptr;
  input_velocity_y_ = nullptr;
  input_velocity_z_ = nullptr;
}

void OpGradQuickIJKScalar_double::calculer_grad_y(const IJK_Field_double& field,
                                                  IJK_Field_double& result)
{
  input_field_ = &field;
  const int ni = field.ni();
  const int nj = field.nj();
  tmp_curv_fram_.allocate(ni, nj, 4, 1);
  compute_grad_y(result);
  input_field_ = nullptr;
  input_velocity_x_ = nullptr;
  input_velocity_y_ = nullptr;
  input_velocity_z_ = nullptr;
}

void OpGradQuickIJKScalar_double::calculer_grad_z(const IJK_Field_double& field,
                                                  IJK_Field_double& result)
{
  input_field_ = &field;
  const int ni = field.ni();
  const int nj = field.nj();
  tmp_curv_fram_.allocate(ni, nj, 4, 1);
  compute_grad_z(result);
  input_field_ = nullptr;
  input_velocity_x_ = nullptr;
  input_velocity_y_ = nullptr;
  input_velocity_z_ = nullptr;
}

void OpGradQuickIJKScalar_double::fill_grad_field_x_y_(IJK_Field_local_double& flux, IJK_Field_double& resu, int k, int dir)
{
  int ni = resu.ni();
  int nj = resu.nj();
  switch(dir)
    {
    case 0:
      for (int i=0; i < ni; i++)
        for (int j=0; j < nj; j++)
          resu(i,j,k) = flux(i+1,j,0) - flux(i,j,0);
      break;
    case 1:
      for (int i=0; i < ni; i++)
        for (int j=0; j < nj; j++)
          resu(i,j,k) = flux(i,j+1,0) - flux(i,j,0);
      break;
    }
}

void OpGradQuickIJKScalar_double::fill_grad_field_z_(IJK_Field_local_double& flux_min, IJK_Field_local_double& flux_max, IJK_Field_double& resu, int k)
{
  int ni = resu.ni();
  int nj = resu.nj();
  for (int i=0; i < ni; i++)
    for (int j=0; j < nj; j++)
      {
        // TODO: What happen if dz is not constant ? The FD operator is no longer define
        double dz_inv = is_flux_ ? 1.: 1 / channel_data_.get_delta_z()[k];
        resu(i,j,k) = (flux_max(i,j,0) - flux_min(i,j,0)) * dz_inv;
      }
}

Implemente_instanciable_sans_constructeur( OpGradFluxQuickIJKScalar_double, "OpGradFluxQuickIJKScalar_double", OpGradQuickIJKScalar_double ) ;

Sortie& OpGradFluxQuickIJKScalar_double::printOn( Sortie& os ) const
{
  return os;
}

Entree& OpGradFluxQuickIJKScalar_double::readOn( Entree& is )
{
  return is;
}

void OpGradFluxQuickIJKScalar_double::calculer_grad_flux(const IJK_Field_double& field,
                                                         const IJK_Field_double& vx,
                                                         const IJK_Field_double& vy,
                                                         const IJK_Field_double& vz,
                                                         FixedVector<IJK_Field_double, 3>& result)
{
  input_velocity_x_ = &vx;
  input_velocity_y_ = &vy;
  input_velocity_z_ = &vz;
  calculer_grad(field, result);
  result.echange_espace_virtuel();
}

void OpGradFluxQuickIJKScalar_double::fill_grad_field_x_y_(IJK_Field_local_double& flux, IJK_Field_double& resu, int k, int dir)
{
  int ni = resu.ni();
  int nj = resu.nj();
  for (int i=0; i < ni; i++)
    for (int j=0; j < nj; j++)
      resu(i,j,k) = flux(i,j,0);
}

void OpGradFluxQuickIJKScalar_double::fill_grad_field_z_(IJK_Field_local_double& flux_min,
                                                         IJK_Field_local_double& flux_max,
                                                         IJK_Field_double& resu, int k)
{
  fill_grad_field_x_y_(flux_min, resu, k, 2);
}
