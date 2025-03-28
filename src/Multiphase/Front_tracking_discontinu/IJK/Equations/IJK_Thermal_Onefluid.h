/****************************************************************************
* Copyright (c) 2023, CEA
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
// File      : IJK_Thermal_Onefluid.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_Thermal_Onefluid_included
#define IJK_Thermal_Onefluid_included

#include <IJK_Thermal_base.h>
#include <IJK_Field_vector.h>
#include <IJK_Field.h>
#include <Boundary_Conditions_Thermique.h>
#include <Domaine_IJK.h>
#include <IJK_Field.h>
#include <Parser.h>
#include <IJK_Lata_writer.h>
#include <OpConvQuickIJKScalar.h>
#include <OpConvCentre2IJKScalar.h>
#include <Ouvrir_fichier.h>
#include <Corrige_flux_FT_base.h>
#include <TRUST_Ref.h>
#include <Operateur_IJK_elem_diff_base.h>
#include <OpConvAmontIJK.h>
#include <OpConvDiscQuickIJKScalar.h>
#include <OpConvCentre4IJK.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_Thermal_Onefluid
//
// <Description of class IJK_Thermal_Onefluid>
//
/////////////////////////////////////////////////////////////////////////////


class IJK_Thermal_Onefluid : public IJK_Thermal_base
{

  Declare_instanciable( IJK_Thermal_Onefluid ) ;

public :

  void initialize(const Domaine_IJK& splitting, const int idx) override;
  double compute_timestep(const double timestep,
                          const double rho_l, const double rho_v,
                          const double dxmin);
  void update_thermal_properties() override;
  void set_param( Param& param ) override;
  virtual void euler_rustine_step(const double timestep) override;
  virtual void rk3_rustine_sub_step(const int rk_step, const double total_timestep,
                                    const double fractionnal_timestep, const double time) override;
  void euler_rustine_step(const double timestep, const double dE);
  void rk3_rustine_sub_step(const int rk_step, const double total_timestep,
                            const double fractionnal_timestep, const double time, const double dE);
  void compute_temperature_convection_conservative(const IJK_Field_vector3_double& velocity)  override ;

protected :
  void compute_dT_rustine(const double dE);
  void compute_T_rust(const IJK_Field_vector3_double& velocity);

  void add_temperature_diffusion() override;
  void compute_diffusion_increment() override;
  /* correct_temperature_for_eulerian_fluxes() May be clearly overridden later */
  void correct_temperature_for_eulerian_fluxes() override { ; };
  double get_rho_cp_ijk(int i, int j, int k) const;
  double get_rho_cp_u_ijk(const IJK_Field_double& vx, int i, int j, int k) const override;
  double compute_rho_cp_u_mean(const IJK_Field_double& vx) override;
  double get_div_lambda_ijk(int i, int j, int k) const override;
  double compute_temperature_dimensionless_theta_mean(const IJK_Field_double& vx) override;

  //Rustine
  double E0_;//volumique

  int deprecated_rho_cp_;
  int rho_cp_moy_harmonic_;
  int lambda_moy_arith_;

  IJK_Field_double T_rust_;
  IJK_Field_double d_T_rustine_; // Temperature increment to conserve the energy.
  IJK_Field_double RK3_F_rustine_; // Temporary storage for substeps in the RK3 algorithm for the rustine calculation.
  IJK_Field_double div_rho_cp_T_;
  IJK_Field_double lambda_;
  IJK_Field_double cp_;
  IJK_Field_double rho_cp_inv_;
};

#endif /* IJK_Thermal_Onefluid_included */
