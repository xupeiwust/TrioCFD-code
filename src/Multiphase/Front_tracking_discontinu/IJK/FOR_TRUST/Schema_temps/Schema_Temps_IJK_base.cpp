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

#include <Schema_Temps_IJK_base.h>
#include <Probleme_FTD_IJK_base.h>
#include <EFichier.h>
#include <Param.h>

Implemente_base(Schema_Temps_IJK_base,"Schema_Temps_IJK_base",Schema_Temps_base);

Sortie& Schema_Temps_IJK_base::printOn(Sortie& os) const
{
  os << "dt " << dt_ << finl;
  os << "temps_courant " << temps_courant_ << finl ;
  os << "tinit " << tinit_ << finl;
  os << "timestep_facsec " << timestep_facsec_ << finl ;
  os << "nb_pas_dt_max " << nb_pas_dt_max_ << finl;
  os << "max_simu_time " << max_simu_time_ << finl ;
  os << "tstep_init " << tstep_init_ << finl;
  os << "use_tstep_init " << use_tstep_init_ << finl ;
  os << "dt_sauvegarde " << dt_sauvegarde_ << finl ;
  os << "cfl " << cfl_ << finl ;
  os << "fo " << fo_ << finl ;
  os << "oh " << oh_ << finl ;
  os << "fin " << finl;
  return os;
}

Entree& Schema_Temps_IJK_base::readOn(Entree& is)
{
  Cerr << "Reading of data for a " << que_suis_je() << " time scheme" << finl;
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);
  temps_courant_ = tinit_;
  lu_ = 1;
  stop_file_ = nom_du_cas() + Nom(".stop");
  return is;
}

void Schema_Temps_IJK_base::set_param(Param& param)
{
//  Schema_Temps_base::set_param(param);

  param.ajouter("tinit", &tinit_, Param::REQUIRED); // XD_ADD_P floattant initial time
  param.ajouter("timestep", &dt_, Param::REQUIRED); // XD_ADD_P floattant Upper limit of the timestep
  param.ajouter("timestep_facsec", &timestep_facsec_); // XD_ADD_P floattant Security factor on timestep
  param.ajouter("cfl", &cfl_); // XD_ADD_P floattant  To provide a value of the limiting CFL number used for setting the timestep
  param.ajouter("fo", &fo_); // XD_ADD_P floattant not_set
  param.ajouter("oh", &oh_); // XD_ADD_P floattant not_set
  param.ajouter("nb_pas_dt_max", &nb_pas_dt_max_, Param::REQUIRED); // XD_ADD_P entier maximum limit for the number of timesteps
  param.ajouter("max_simu_time", &max_simu_time_); // XD_ADD_P double maximum limit for the simulation time
  param.ajouter("tstep_init", &tstep_init_); // XD_ADD_P entier index first interation for recovery
  param.ajouter("use_tstep_init", &use_tstep_init_); // XD_ADD_P entier use tstep init for constant post-processing step
  param.ajouter("dt_sauvegarde", &dt_sauvegarde_); // XD_ADD_P entier saving frequency (writing files for computation restart)
  param.ajouter_non_std( "tcpumax",(this)); // XD_ADD_P double CPU time limit (must be specified in hours) for which the calculation is stopped (1e30s by default).

  param.ajouter_flag("enable_dt_oh_ideal_length_factor", &enable_dt_oh_ideal_length_factor_);
  param.ajouter_flag("first_step_interface_smoothing", &first_step_interface_smoothing_);
}

void Schema_Temps_IJK_base::set_param_reprise_pb(Param& param)
{
  param.ajouter("tinit", &temps_courant_);
  param.ajouter("tstep_init", &tstep_init_);
}

void Schema_Temps_IJK_base::completer()
{
  first_step_interface_smoothing_ = (first_step_interface_smoothing_ &&
                                     (!ref_cast(Probleme_FTD_IJK_base, mon_probleme.valeur()).get_reprise() && temps_courant_ == 0.));
  if (tstep_init_)
    use_tstep_init_ = 1;
}

// Methode de calcul du pas de temps max base sur CFL, Oh et Fo
// Pour les maillages uniformes uniquement.
// On choisi de mettre le facteur 0.5 dans dt_cfl et dt_fo
// pour que le calcul soit stable avec un CFL <=1.0 et Fo <= 1.0.
// Sinon, il faudrait recommander CFL <= 0.5 et Fo <=0.5 ce qui n'est pas la valeur par defaut...
double Schema_Temps_IJK_base::find_timestep(const double max_timestep, const double cfl, const double fo, const double oh)
{
  statistiques().begin_count(dt_counter_);

  // XXX pffff pas const car Thermals classe est tout xxxx ....
  Probleme_FTD_IJK_base& pb_ijk = ref_cast(Probleme_FTD_IJK_base, mon_probleme.valeur());

  dt_cfl_ = 1.e20;
  double dxmin = 1.e20;
  // On ne connait pas la longueur minimum du maillage lagrangien, mais on en prend une approximation :
  // lg \approx 1.7 delta
  double lg_cube = 1.;
  double lg_cube_raw = 1.;
  for (int dir = 0; dir < 3; dir++)
    {
      const IJK_Field_double& v = pb_ijk.eq_ns().get_velocity()[dir];
      double max_v = 1.e-20; // Pas zero pour la division a la fin...
      const int ni = v.ni();
      const int nj = v.nj();
      const int nk = v.nk();
      for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
          for (int i = 0; i < ni; i++)
            max_v = std::max(max_v, fabs(v(i, j, k)));
      max_v = Process::mp_max(max_v);
      const Domaine_IJK& geom = v.get_domaine();
#ifndef VARIABLE_DZ
      const double delta = geom.get_constant_delta(dir);
#else
      const ArrOfDouble& tab_dz=geom.get_delta(dir);
      const double delta = Process::mp_min(min_array(tab_dz));
#endif
      lg_cube *= 1.7 * delta;
      lg_cube_raw *= delta;
      dxmin = std::min(delta, dxmin);
      // QUESTION GAB : pourquoi on ne reprend pas dxmin ?
      if (max_v > 0)
        dt_cfl_ = std::min(dt_cfl_, delta / max_v * 0.5);
    }
  dt_cfl_ *= cfl;
  dt_cfl_liq_ = dt_cfl_;
  dt_cfl_vap_ = dt_cfl_;

  const double mu_l = pb_ijk.milieu_ijk().get_mu_liquid(), rho_l = pb_ijk.milieu_ijk().get_rho_liquid(), mu_v = pb_ijk.milieu_ijk().get_mu_vapour(), rho_v = pb_ijk.milieu_ijk().get_rho_vapour();

  dt_fo_liq_ = dxmin * dxmin / ((mu_l / rho_l) + 1.e-20) * fo * 0.125;
  dt_fo_vap_ = dxmin * dxmin / ((mu_v / rho_v) + 1.e-20) * fo * 0.125;
  dt_fo_ = std::min(dt_fo_liq_, dt_fo_vap_);
  if (pb_ijk.eq_ns().get_disable_diffusion_qdm())
    dt_fo_ = 1.e20;

  /*
   * Popinet et.al 2018 (review surface tension calculation)
   * Au cas ou sigma = 0, on utilise (sigma + 1e-20)
   */
  double ideal_length_factor = 1.7;

  if (pb_ijk.has_interface())
    ideal_length_factor = pb_ijk.get_remaillage_ft_ijk().get_facteur_longueur_ideale();

  double max_sigma = -1e10;
  /*if (!interfaces_.maillage_ft_ijk().Surfactant_facettes().get_disable_surfactant())
   {
   // alors on a une tension de surface variable
   // Il faut prendre le max de la tension de surface dans le critere de stabilite ?
   const int nb_som=interfaces_.maillage_ft_ijk().nb_sommets();
   const ArrOfDouble& sigma_sommets = interfaces_.maillage_ft_ijk().Surfactant_facettes().get_sigma_sommets();
   for (int i = 0; i < nb_som; i++)
   {
   if (max_sigma<sigma_sommets[i])
   max_sigma= sigma_sommets[i];
   }
   max_sigma=Process::mp_max(max_sigma);
   }
   else*/
  max_sigma = pb_ijk.milieu_ijk().sigma();

  if (enable_dt_oh_ideal_length_factor_)
    dt_oh_ = sqrt((rho_l + rho_v) / (2 * M_PI) * lg_cube_raw * pow(ideal_length_factor, 3) / (max_sigma + 1e-20)) * oh;
  else
    dt_oh_ = sqrt((rho_l + rho_v) / 2. * lg_cube / (max_sigma + 1e-20)) * oh;
  if (pb_ijk.eq_ns().get_disable_convection_qdm())
    dt_oh_ = 1.e20;

  // Seems underestimated !
  //  const double dt_oh  = sqrt((rho_liquide_+rho_vapeur_)/2. * lg_cube/(sigma_+1e-20) ) * oh;

  const double dt_eq_velocity = 1. / (1. / dt_cfl_ + 1. / dt_fo_ + 1. / dt_oh_);


  double dt_thermals = 1.e20;

  if (pb_ijk.has_thermals())
    pb_ijk.get_ijk_thermals().compute_timestep(dt_thermals, dxmin);

  const double dt = std::min(max_timestep, timestep_facsec_ * std::min(dt_eq_velocity, dt_thermals));

  if (Process::je_suis_maitre())
    {
      int reset = (!pb_ijk.get_reprise()) && (nb_pas_dt_ == 0);
      SFichier fic = Ouvrir_fichier(".dt_ev", "tstep\ttime\ttimestep\tdt_cfl\tdt_fo\tdt_oh\tdt_diff_th", reset);
      fic << nb_pas_dt_ << " " << temps_courant_ << " " << dt;
      fic << " " << dt_cfl_ << " " << dt_fo_ << " " << dt_oh_;
      fic << " " << dt_thermals; // If no thermal equation, value will be large.
      fic << finl;
      fic.close();
    }
  statistiques().end_count(dt_counter_);

  /* a bouger un jour dans mettre_a_jour ... existe pas pour le moment */
  if (!ind_temps_cpu_max_atteint)
    ind_temps_cpu_max_atteint = (temps_cpu_ecoule_ >= tcpumax_);

  if (je_suis_maitre())
    temps_cpu_ecoule_ = statistiques().last_time(temps_total_execution_counter_);

  envoyer_broadcast(temps_cpu_ecoule_,0);

  assert(dt > 0);

  return dt;
}

void Schema_Temps_IJK_base::check_stop_criteria(bool& stop) const
{
  stop = false;
  int stop_i = 0;

  // 1. verification du fichier stop
  if (!get_disable_stop())
    {
      if (je_suis_maitre())
        {
          EFichier f;
          stop_i = f.ouvrir(stop_file_);
          if (stop_i) // file exists, check if it contains 1
            f >> stop_i;
        }
      envoyer_broadcast(stop_i, 0);

      if (stop_i)
        {
          Cerr << "---------------------------------------------------------" << finl
               << "The problem " << pb_base().le_nom() << " wants to stop : stop file detected" << finl << finl;
        }
    }

  // 2. verification nb_pas_dt_max
  if (get_tstep() == get_nb_timesteps() - 1)
    {
      stop_i = 1;
      Cerr << "---------------------------------------------------------" << finl
           << "The problem " << pb_base().le_nom() << " wants to stop : the maximum number of time steps reached" << finl << finl;
    }

  // 3. verification tmax
  if (get_current_time() >= get_max_simu_time())
    {
      stop_i = 1;
      Cerr << "---------------------------------------------------------" << finl
           << "The problem " << pb_base().le_nom() << " wants to stop : final time reached" << finl << finl;
    }

  // 4. verification tcpu_max
  if (ind_temps_cpu_max_atteint)
    {
      stop_i = 1;
      Cerr << "---------------------------------------------------------" << finl
           << "The problem " << pb_base().le_nom() << " wants to stop : max cpu time reached" << finl << finl;
    }

  stop = stop_i;
}

double Schema_Temps_IJK_base::computeTimeStep(bool& stop) const
{
  const Probleme_FTD_IJK_base& pb_ijk = ref_cast(Probleme_FTD_IJK_base, mon_probleme.valeur());

  if (timestep_facsec_ > 0. && !stop)
    {
      Schema_Temps_IJK_base& sh_non_cst = const_cast<Schema_Temps_IJK_base&>(*this); // FIXME
      double max_post_simu_timestep = pb_ijk.get_post().get_timestep_simu_post(get_current_time(), max_simu_time_);
      // WTF - find_timestep non const ?!!
      sh_non_cst.dt_ = sh_non_cst.find_timestep(std::min(max_timestep_, max_post_simu_timestep), cfl_, fo_, oh_);
    }

  return dt_;
}
