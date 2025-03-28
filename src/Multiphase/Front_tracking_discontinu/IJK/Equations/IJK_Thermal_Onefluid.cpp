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
// File      : IJK_Thermal_Onefluid.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Navier_Stokes_tools.h>
#include <IJK_Thermal_Onefluid.h>
#include <Probleme_FTD_IJK.h>
#include <Schema_RK3_IJK.h>
#include <DebogIJK.h>

Implemente_instanciable_sans_constructeur( IJK_Thermal_Onefluid, "IJK_Thermal_Onefluid", IJK_Thermal_base ) ;

IJK_Thermal_Onefluid::IJK_Thermal_Onefluid()
{
  needs_op_unform_ = 0;
  lambda_moy_arith_=0;
  type_temperature_convection_form_ = 0;  // Default value: 0 : non conservative
  conserv_energy_global_=0;
  rho_cp_moy_harmonic_=0;
  rho_cp_post_=0;
  E0_=0;
  deprecated_rho_cp_=0;
}

Sortie& IJK_Thermal_Onefluid::printOn( Sortie& os ) const
{
  IJK_Thermal_base::printOn( os );
  os<< "  {\n";

  os<< "    type_T_source " << type_T_source_ << "\n";

  if (rho_cp_moy_harmonic_)
    os<< "    rho_cp_moy_harmonic \n";
  if (lambda_variable_)
    os<< "    lambda_variable \n";
  if (lambda_moy_arith_)
    os<< "    lambda_moy_arith \n";
  if (deprecated_rho_cp_)
    os<< "    depracated_rho_cp \n";
  if (conserv_energy_global_)
    os<< "    conserv_energy_global \n";

  os<< "  \n}";
  return os;
}

Entree& IJK_Thermal_Onefluid::readOn( Entree& is )
{
  IJK_Thermal_base::readOn( is );
  return is;
}

void IJK_Thermal_Onefluid::set_param( Param& param )
{
  IJK_Thermal_base::set_param(param);
  param.ajouter_flag("conserv_energy_global", &conserv_energy_global_);
  param.ajouter_flag("lambda_moy_arith_", &lambda_moy_arith_);
  param.ajouter_flag("rho_cp_moy_harmonic", &rho_cp_moy_harmonic_);
  param.ajouter_flag("deprecated_rho_cp", &deprecated_rho_cp_);
}

void IJK_Thermal_Onefluid::initialize(const Domaine_IJK& splitting, const int idx)
{
  Cout << que_suis_je() << "::initialize()" << finl;
  IJK_Thermal_base::initialize(splitting, idx);
  temperature_diffusion_op_.set_conductivity_coefficient(uniform_lambda_, lambda_, *temperature_, *temperature_, *temperature_);
  lambda_.allocate(splitting, Domaine_IJK::ELEM, 1);

  if (rho_cp_moy_harmonic_)
    rho_cp_inv_.allocate(splitting, Domaine_IJK::ELEM, 2);
  else
    {
      if (deprecated_rho_cp_)
        cp_.allocate(splitting, Domaine_IJK::ELEM, 2);
      else
        {
          if (!rho_cp_post_ || (type_temperature_convection_form_==2))
            rho_cp_.allocate(splitting, Domaine_IJK::ELEM, 2);
        }
    }

  // Compute initial energy :
  if (conserv_energy_global_)
    {
      E0_ = compute_global_energy(*temperature_);
      d_T_rustine_.allocate(splitting, Domaine_IJK::ELEM, 1);
      if (ref_ijk_ft_.non_nul() && sub_type(Schema_RK3_IJK, ref_ijk_ft_->schema_temps_ijk()))
        RK3_F_rustine_.allocate(splitting, Domaine_IJK::ELEM, 0);
      Cout << "Initial energy at time t=" << ref_ijk_ft_->schema_temps_ijk().get_current_time() << " is " << E0_ << " [W.m-3]." << finl;
      Cerr << "Initial energy at time t=" << ref_ijk_ft_->schema_temps_ijk().get_current_time() << " is " << E0_ << " [W.m-3]." << finl;
    }
  if (type_temperature_convection_form_==1)
    {
      rho_cp_T_.allocate(splitting, Domaine_IJK::ELEM, 2);
      div_rho_cp_T_.allocate(splitting, Domaine_IJK::ELEM, 0);
    }

  const Nom nom_rust= get_field_name_with_rank("T_RUST");
  if (ref_ijk_ft_post_->is_post_required(nom_rust) || (type_temperature_convection_form_==2))
    {
      T_rust_.allocate(splitting, Domaine_IJK::ELEM, 0, nom_rust);
      champs_compris_.ajoute_champ(T_rust_);
      T_rust_.data() = 0.;
    }
  Cout << "End of " << que_suis_je() << "::initialize()" << finl;
}

void IJK_Thermal_Onefluid::update_thermal_properties()
{
  IJK_Thermal_base::update_thermal_properties();
  const IJK_Field_double& indic = ref_ijk_ft_->get_interface().I();
  // Nombre de mailles du domaine NS :
  const int nx = indic.ni();
  const int ny = indic.nj();
  const int nz = indic.nk();
  const double rho_l = ref_ijk_ft_->milieu_ijk().get_rho_liquid();
  const double rho_v = ref_ijk_ft_->milieu_ijk().get_rho_vapour();
  const bool harmonic_mean = ((!lambda_moy_arith_) and (lambda_liquid_ > DMINFLOAT) and (lambda_vapour_ >DMINFLOAT));
  const bool rho_cp_harmonic_mean = ((rho_cp_moy_harmonic_) and (rho_l*cp_liquid_ > DMINFLOAT) and (rho_v*cp_vapour_ >DMINFLOAT));
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        {
          double chi_l = indic(i,j,k);
          if (rho_cp_moy_harmonic_)
            if (rho_cp_harmonic_mean)
              rho_cp_inv_(i,j,k)= chi_l / (rho_l * cp_liquid_)  + (1 - chi_l) / (rho_v * cp_vapour_);
            else
              rho_cp_inv_(i,j,k) = rho_l * cp_liquid_ * chi_l + rho_v * cp_vapour_ * (1.- chi_l);
          else
            {
              if (deprecated_rho_cp_)
                cp_(i,j,k) = cp_liquid_ * chi_l + cp_vapour_ * (1.- chi_l);
              else if (!rho_cp_post_)
                rho_cp_(i,j,k) = rho_l * cp_liquid_ * chi_l + rho_v * cp_vapour_ * (1.- chi_l);
            }
          if (rho_cp_post_)
            rho_cp_(i,j,k) = rho_l * cp_liquid_ * chi_l + rho_v * cp_vapour_ * (1.- chi_l);
          if (harmonic_mean)
            lambda_(i,j,k) = lambda_liquid_ * lambda_vapour_ / ((1 - chi_l) * lambda_liquid_ + chi_l * lambda_vapour_);
          else
            lambda_(i,j,k) = lambda_liquid_  * chi_l + (1.- chi_l) * lambda_vapour_ ;
        }
  lambda_.echange_espace_virtuel(lambda_.ghost());
  if (rho_cp_post_)
    rho_cp_.echange_espace_virtuel(rho_cp_.ghost());
  if (rho_cp_moy_harmonic_)
    rho_cp_inv_.echange_espace_virtuel(rho_cp_inv_.ghost());
  else
    {
      if (deprecated_rho_cp_)
        cp_.echange_espace_virtuel(cp_.ghost());
      else
        {
          if (!rho_cp_post_)
            rho_cp_.echange_espace_virtuel(rho_cp_.ghost());
        }
    }
}

void IJK_Thermal_Onefluid::add_temperature_diffusion()
{
  lambda_.echange_espace_virtuel(lambda_.ghost());
  DebogIJK::verifier("lambda", lambda_);

  temperature_diffusion_op_->set_lambda(lambda_);
  IJK_Thermal_base::add_temperature_diffusion();
}

void IJK_Thermal_Onefluid::compute_diffusion_increment()
{
  const IJK_Field_double& div_coeff_grad_T_volume = *div_coeff_grad_T_volume_;
  IJK_Field_double& d_temperature                 = *d_temperature_;

  // Update d_temperature
  const int ni = d_temperature.ni();
  const int nj = d_temperature.nj();
  const int nk = d_temperature.nk();
  const double vol_inv = 1./vol_;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          if (rho_cp_moy_harmonic_)
            {
              const double rhocpV_inv = rho_cp_inv_(i,j,k) * vol_inv;
              const double ope = div_coeff_grad_T_volume(i,j,k);
              const double resu = ope*rhocpV_inv;
              d_temperature(i,j,k) +=resu ;
            }
          else
            {
              double rhocpV = 0;
              if (deprecated_rho_cp_)
                {
                  const double rho = ref_ijk_ft_->eq_ns().get_rho_field_ijk(i,j,k);
                  const double cp = cp_(i,j,k);
                  rhocpV = rho * cp * vol_;
                }
              else
                {
                  rhocpV = rho_cp_(i,j,k) * vol_;
                }
              const double ope = div_coeff_grad_T_volume(i,j,k);
              const double resu = ope/rhocpV;
              d_temperature(i,j,k) +=resu ;
            }
        }
}

double IJK_Thermal_Onefluid::compute_rho_cp_u_mean(const IJK_Field_double& vx)
{
  double rho_cp_u_mean = 0.;
  if (rho_cp_moy_harmonic_)
    {
      const double rho_l = ref_ijk_ft_->milieu_ijk().get_rho_liquid();
      const double rho_v = ref_ijk_ft_->milieu_ijk().get_rho_vapour();
      const bool rho_cp_harmonic_mean = ((rho_cp_moy_harmonic_) and (rho_l*cp_liquid_ > DMINFLOAT) and (rho_v*cp_vapour_ >DMINFLOAT));
      if (rho_cp_harmonic_mean)
        rho_cp_u_mean = calculer_rho_cp_u_moyen(vx, rho_cp_inv_, vx, 0., 3);
      else
        rho_cp_u_mean = calculer_rho_cp_u_moyen(vx, rho_cp_inv_, vx, 0., 1);
    }
  else
    {
      if (deprecated_rho_cp_)
        {
          rho_cp_u_mean = calculer_rho_cp_u_moyen(vx, cp_, ref_ijk_ft_->eq_ns().get_rho_field(), 0., 0);
        }
      else
        {
          rho_cp_u_mean = calculer_rho_cp_u_moyen(vx, rho_cp_, vx, 0., 1);
        }
    }
  return rho_cp_u_mean;
}

double IJK_Thermal_Onefluid::get_rho_cp_ijk(int i, int j, int k) const
{
  double rho_cp = 0.;
  if (rho_cp_moy_harmonic_)
    {
      const double rho_l = ref_ijk_ft_->milieu_ijk().get_rho_liquid();
      const double rho_v = ref_ijk_ft_->milieu_ijk().get_rho_vapour();
      const bool rho_cp_harmonic_mean = ((rho_cp_moy_harmonic_) and (rho_l*cp_liquid_ > DMINFLOAT) and (rho_v*cp_vapour_ >DMINFLOAT));
      if (rho_cp_harmonic_mean)
        rho_cp = 1 / rho_cp_inv_(i,j,k);
      else
        rho_cp = rho_cp_inv_(i,j,k);
    }
  else
    {
      if (deprecated_rho_cp_)
        {
          rho_cp = cp_(i,j,k) * (ref_ijk_ft_->eq_ns().get_rho_field_ijk(i,j,k));
        }
      else
        {
          rho_cp = rho_cp_(i,j,k);
        }
    }
  return rho_cp;
}

double IJK_Thermal_Onefluid::get_rho_cp_u_ijk(const IJK_Field_double& vx, int i, int j, int k) const
{
  return get_rho_cp_ijk(i,j,k) * vx(i,j,k);
}

double IJK_Thermal_Onefluid::get_div_lambda_ijk(int i, int j, int k) const
{
  return (lambda_(i+1,j,k)-lambda_(i-1,j,k));
}

double IJK_Thermal_Onefluid::compute_temperature_dimensionless_theta_mean(const IJK_Field_double& vx)
{
  double theta_adim_moy = 0.;
  if (rho_cp_moy_harmonic_)
    {
      const double rho_l = ref_ijk_ft_->milieu_ijk().get_rho_liquid();
      const double rho_v = ref_ijk_ft_->milieu_ijk().get_rho_vapour();
      const bool rho_cp_harmonic_mean = ((rho_cp_moy_harmonic_) and (rho_l*cp_liquid_ > DMINFLOAT) and (rho_v*cp_vapour_ >DMINFLOAT));
      if (rho_cp_harmonic_mean)
        theta_adim_moy = calculer_temperature_adimensionnelle_theta_moy(vx, temperature_adimensionnelle_theta_, rho_cp_inv_, vx, 0., 3);
      else
        theta_adim_moy = calculer_temperature_adimensionnelle_theta_moy(vx, temperature_adimensionnelle_theta_, rho_cp_inv_, vx, 0., 1);
    }
  else
    {
      if (deprecated_rho_cp_)
        {
          theta_adim_moy = calculer_temperature_adimensionnelle_theta_moy(vx, temperature_adimensionnelle_theta_, cp_, ref_ijk_ft_->eq_ns().get_rho_field(), 0., 0);
        }
      else
        {
          theta_adim_moy = calculer_temperature_adimensionnelle_theta_moy(vx, temperature_adimensionnelle_theta_, rho_cp_, vx, 0., 1);
        }
    }
  return theta_adim_moy;
}

void IJK_Thermal_Onefluid::compute_dT_rustine(const double dE)
{
  const int ni = T_rust_.ni();
  const int nj = T_rust_.nj();
  const int nk = T_rust_.nk();
  const double rho_l = ref_ijk_ft_->milieu_ijk().get_rho_liquid();
  const double rho_v = ref_ijk_ft_->milieu_ijk().get_rho_vapour();
  const IJK_Field_double& indic = ref_ijk_ft_->get_interface().I();
  double int_rhocpTrust = 0;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          int_rhocpTrust +=  (rho_l*cp_liquid_*indic(i,j,k) + rho_v*cp_vapour_*(1.-indic(i,j,k)))*T_rust_(i,j,k);
        }
  const int ntot = T_rust_.get_domaine().get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_I)
                   *T_rust_.get_domaine().get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_J)
                   *T_rust_.get_domaine().get_nb_items_global(Domaine_IJK::ELEM, DIRECTION_K);
  int_rhocpTrust = mp_sum(int_rhocpTrust)/(double)(ntot);
  if ((std::fabs(int_rhocpTrust) < DMINFLOAT) && (std::fabs(dE) > DMINFLOAT))
    {
      Cerr << "Trying to normalize the energy by an infinite coef : " << dE << " / " << int_rhocpTrust << finl;
      Process::exit();
    }
  else
    {
      Cerr << "Le coeff de manque d'energie dE/int_rhocpTrust vaut : " << dE/int_rhocpTrust << finl;
    }
  if (int_rhocpTrust)
    {
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              d_T_rustine_(i,j,k) = dE/ int_rhocpTrust * T_rust_(i,j,k);
            }
    }
}

void IJK_Thermal_Onefluid::rk3_rustine_sub_step(const int rk_step, const double total_timestep,
                                                const double fractionnal_timestep, const double time)
{
  update_thermal_properties();
  if (get_conserv_energy_global())
    {
      const double dE = get_E0() - compute_global_energy();
      rk3_rustine_sub_step(rk_step, total_timestep, fractionnal_timestep, time, dE);
    }
}

void IJK_Thermal_Onefluid::rk3_rustine_sub_step(const int rk_step, const double total_timestep,
                                                const double fractionnal_timestep, const double time, const double dE)
{
  compute_dT_rustine(dE);
  // Update the temperature :
  const int kmax = temperature_->nk();
  const double ene_ini = compute_global_energy();
  for (int k = 0; k < kmax; k++)
    {
      runge_kutta3_update(d_T_rustine_, RK3_F_rustine_, *temperature_, rk_step, k, total_timestep);
    }
  temperature_->echange_espace_virtuel(temperature_->ghost());
  const double ene_post = compute_global_energy();
  Cerr << "[Energy-Budget-T"<<rang_<<"RK3 rustine step "<<rk_step<<"] time t=" << ref_ijk_ft_->schema_temps_ijk().get_current_time()
       << " " << ene_ini
       << " " << ene_post << " [W.m-3]. [step"<< rk_step << "]" << finl;
  source_callback();
}

void IJK_Thermal_Onefluid::euler_rustine_step(const double timestep)
{
  update_thermal_properties();
  if (get_conserv_energy_global())
    {
      const double dE = get_E0() - compute_global_energy();
      euler_rustine_step(timestep, dE);
    }
}

void IJK_Thermal_Onefluid::euler_rustine_step(const double timestep, const double dE)
{
  compute_dT_rustine(dE);
  // Update the temperature :
  const int kmax = temperature_->nk();
  const double ene_ini = compute_global_energy();
  for (int k = 0; k < kmax; k++)
    ref_ijk_ft_->eq_ns().euler_explicit_update(d_T_rustine_, *temperature_, k);
  temperature_->echange_espace_virtuel(temperature_->ghost());
  const double ene_post = compute_global_energy();
  Cerr << "[Energy-Budget-T"<<rang_<<" euler rustine] time t=" << ref_ijk_ft_->schema_temps_ijk().get_current_time()
       << " " << ene_ini
       << " " << ene_post << " [W.m-3]."
       << " dE "<< dE
       << finl;
  source_callback();
}

void IJK_Thermal_Onefluid::compute_temperature_convection_conservative(const IJK_Field_vector3_double& velocity)
{
  static Stat_Counter_Id cnt_conv_temp = statistiques().new_counter(1, "FT convection rho");
  statistiques().begin_count(cnt_conv_temp);
  IJK_Field_double& d_temperature = *d_temperature_;
  if (conv_temperature_negligible_)
    {
      d_temperature.data()=0;
    }
  else
    {
      // TODO Mathis : Voir la nouvelle facon sans switch
      /*
        switch (type_temperature_convection_op_)
          {
          case 1:
            Cerr << "Operateur amont non implemente" << finl;
            break;
          case 2:
            Cerr << "Operateur centre non implemente" << finl;
            break;
          case 3:
            //          rho_cp_convection_op_quick_.calculer(rho_cp_T_, velocity[0], ref_ijk_ft_->itfce().I(), velocity[1], velocity[2], div_rho_cp_T_);
            rho_cp_convection_op_quick_.set_indicatrice(ref_ijk_ft_->itfce().I());
            rho_cp_convection_op_quick_.calculer(rho_cp_T_, velocity[0], velocity[1], velocity[2], div_rho_cp_T_);

            break;
          case 4:
            Cerr << "Operateur centre4 non implemente" << finl;
            break;

          default:
            Cerr << "Unknown convection operator for the temperature." << finl;
            Process::exit();
          } */
      const int ni = d_temperature.ni();
      const int nj = d_temperature.nj();
      const int nk = d_temperature.nk();
      const Domaine_IJK& geom = d_temperature.get_domaine();
      const double dx = geom.get_constant_delta(DIRECTION_I);
      const double dy = geom.get_constant_delta(DIRECTION_J);
      const double dz = geom.get_constant_delta(DIRECTION_K);
      const double vol = dx*dy*dz;
      const double rho_l = ref_ijk_ft_->milieu_ijk().get_rho_liquid();
      const double rho_v = ref_ijk_ft_->milieu_ijk().get_rho_vapour();
      const IJK_Field_double& indic = ref_ijk_ft_->get_interface().I();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              const double I = indic(i,j,k);
              // dT = 1/rho_cp_h * div( rho*cp*T*v )
              div_rho_cp_T_(i,j,k) /= vol;
              d_temperature(i,j,k) = 1/((I*rho_l*cp_liquid_) + (1.-I)*rho_v*cp_vapour_)*div_rho_cp_T_(i,j,k);
            }
      // on met a jour T_rust_
      compute_T_rust(velocity);
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            {
              // dT = 1/rho_cp_h * div( rho*cp*T*v ) - T/rho_cp_h *div( rho*cp*v)
              d_temperature(i,j,k) -= T_rust_(i,j,k) ;
            }
      d_temperature.echange_espace_virtuel(d_temperature.ghost());
    }
  statistiques().end_count(cnt_conv_temp);
  DebogIJK::verifier("op_conv(rho)", d_temperature);
  return;
}

void IJK_Thermal_Onefluid::compute_T_rust(const IJK_Field_vector3_double& velocity)
{
  static Stat_Counter_Id cnt_conv_temp = statistiques().new_counter(1, "FT convection rho cp");
  statistiques().begin_count(cnt_conv_temp);
  const Domaine_IJK& geom = T_rust_.get_domaine();
  // DONE: remplacer rho_cp par un champ rho_cp_ mis a jour dans update_thermal_properties. Necessaire pour que ca marche.
  //On calcule div(rho_cp*v) qu'on stocke dans T_rust
  // TODO Mathis check
  /* switch (type_temperature_convection_op_)
    {
    case 1:
      temperature_convection_op_amont_.calculer(rho_cp_, rho_cp_, rho_cp_, velocity[0], velocity[1], velocity[2], T_rust_, T_rust_, T_rust_);
      break;
    case 2:
      temperature_convection_op_centre2_.calculer(rho_cp_, velocity[0], velocity[1], velocity[2], T_rust_);
      break;
    case 3:
      temperature_convection_op_quick_.calculer(rho_cp_, velocity[0], velocity[1], velocity[2], T_rust_);
      break;
    case 4:
      temperature_convection_op_centre4_.calculer(rho_cp_,rho_cp_,rho_cp_, velocity[0], velocity[1], velocity[2], T_rust_, T_rust_, T_rust_);
      break;

    default:
      Cerr << "Unknown convection operator for the temperature." << finl;
      Process::exit();
    }
      */
  temperature_convection_op_->calculer(T_rust_, velocity[0], velocity[1], velocity[2], d_T_rustine_);

  // TODO : GB il manquait un truc comme ca? ou
  const int kmax = temperature_->nk();
  for (int k = 0; k < kmax; k++)
    ref_ijk_ft_->eq_ns().euler_explicit_update(d_T_rustine_, T_rust_, k);
  // ou bien :
//  const double timestep = ref_ijk_ft_->schema_temps_ijk().get_timestep();
//  euler_rustine_step(timestep, 0. /* dE */);
  //
  const int ni = T_rust_.ni();
  const int nj = T_rust_.nj();
  const int nk = T_rust_.nk();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  // To be sure we're on a regular mesh
  assert(dz >0);
  const double vol = dx*dy*dz;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          T_rust_(i,j,k) /=vol ;
        }
  //On met a jour T_rust_ on le multipliant par T/rho_cp_hamro
  const double rho_l = ref_ijk_ft_->milieu_ijk().get_rho_liquid();
  const double rho_v = ref_ijk_ft_->milieu_ijk().get_rho_vapour();
  const IJK_Field_double& indic = ref_ijk_ft_->get_interface().I();
  const IJK_Field_double& temperature = *temperature_;
  for (int k=0; k < nk ; k++)
    for (int j=0; j< nj; j++)
      for (int i=0; i < ni; i++)
        {
          const double chi_l = indic(i,j,k);
          T_rust_(i,j,k) *= temperature(i,j,k)/(chi_l*rho_l*cp_liquid_ + (1-chi_l)*rho_v*cp_vapour_);
        }
  statistiques().end_count(cnt_conv_temp);
  return;
}
