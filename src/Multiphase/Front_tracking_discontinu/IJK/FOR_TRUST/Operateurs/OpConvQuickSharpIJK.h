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

#ifndef OpConvQuickSharpIJK_included
#define OpConvQuickSharpIJK_included
#include <Operateur_IJK_faces_conv_base.h>


class OpConvQuickSharpIJK_double : public Operateur_IJK_faces_conv_base_double
{
  Declare_instanciable(OpConvQuickSharpIJK_double);
public:
  void initialize(const Domaine_IJK& splitting) override;

protected:
  inline void compute_flux_x_vx(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::X,DIRECTION::X>(resu,k_layer);
  }
  inline void compute_flux_x_vy(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::X,DIRECTION::Y>(resu,k_layer);
  }
  inline void compute_flux_x_vz(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::X,DIRECTION::Z>(resu,k_layer);
  }
  inline void compute_flux_y_vx(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::Y,DIRECTION::X>(resu,k_layer);
  }
  inline void compute_flux_y_vy(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::Y,DIRECTION::Y>(resu,k_layer);
  }
  inline void compute_flux_y_vz(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::Y,DIRECTION::Z>(resu,k_layer);
  }
  inline void compute_flux_z_vx(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::Z,DIRECTION::X>(resu,k_layer);
  }
  inline void compute_flux_z_vy(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::Z,DIRECTION::Y>(resu,k_layer);
  }
  inline void compute_flux_z_vz(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::Z,DIRECTION::Z>(resu,k_layer);
  }

  double delta_x_ = 0;
  double delta_y_ = 0;
  double delta_z_ = 0; // coded for uniform mesh

private:
  template <DIRECTION _DIR_, DIRECTION _VCOMPO_>
  void compute_flux_(IJK_Field_local_double& resu, const int k_layer);

  inline double get_delta(DIRECTION _DIR_) const
  {
    switch(_DIR_)
      {
      case DIRECTION::X:
        return delta_x_;
      case DIRECTION::Y:
        return delta_y_;
      case DIRECTION::Z:
        return delta_z_;
      default:
        Cerr << "Error in OpConvIJKQuickSsarp::get_delta: wrong direction..." << finl;
        Process::exit();
      }
    //for compilation only...
    return 0.0;
  }

};

#include <OpConvQuickSharpIJK.tpp>

#endif /* OpConvQuickSharpIJK_included */
