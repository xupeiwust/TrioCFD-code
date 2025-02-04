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

#ifndef Schema_Temps_IJK_base_included
#define Schema_Temps_IJK_base_included

#include <Schema_Temps_base.h>

class Schema_Temps_IJK_base: public Schema_Temps_base
{
  Declare_base(Schema_Temps_IJK_base);
public :
  // Renvoie le nombre de valeurs temporelles a conserver. Ici : n et n+1, donc 2.
  int nb_valeurs_temporelles() const override { return 2 ; }

  // Renvoie le nombre de valeurs temporelles futures. Ici : n+1, donc 1.
  int nb_valeurs_futures() const override { return 1 ; }

  // Renvoie le le temps a la i-eme valeur future. Ici : t(n+1)
  double temps_futur(int i) const override
  {
    assert(i == 1);
    return temps_courant() + pas_de_temps();
  }

  // Renvoie le le temps le temps que doivent rendre les champs a l'appel de valeurs(). Ici : t(n+1)
  double temps_defaut() const override { return temps_courant() + pas_de_temps(); }

  // a surcharger si utile
  void completer() override;
  double computeTimeStep(bool& stop) const override;

  void set_param(Param& ) override;
  void set_param_reprise_pb(Param& );

  double get_dt_cfl() const { return dt_cfl_; }
  double get_dt_fo() const { return dt_fo_; }
  double get_dt_oh() const { return dt_oh_; }
  double get_dt_cfl_liq() const { return dt_cfl_liq_; }
  double get_dt_cfl_vap_() const { return dt_cfl_vap_; }
  double get_dt_fo_liq() const { return dt_fo_liq_; }
  double get_dt_fo_vap_() const { return dt_fo_vap_; }
  double get_timestep_facsec() const { return timestep_facsec_; }
  double get_modified_time_ini() const { return modified_time_ini_; }
  double get_max_simu_time() const { return max_simu_time_; }
  double get_current_time() const { return temps_courant(); }
  double get_timestep() const { return pas_de_temps(); } // en double
  double& set_timestep() { return set_dt(); } // en double

  int get_nb_timesteps() const { return nb_pas_dt_max() ; }
  int get_tstep() const { return nb_pas_dt(); }
  int& get_tstep() { return nb_pas_dt_; }
  int get_tstep_init() const { return tstep_init_; }
  int get_tstep_sauv() const { return tstep_sauv_; }
  int get_use_tstep_init() const { return use_tstep_init_; }
  int& get_first_step_interface_smoothing() { return first_step_interface_smoothing_; }
  int get_first_step_interface_smoothing() const { return first_step_interface_smoothing_; }
  int get_enable_dt_oh_ideal_length_factor() const { return enable_dt_oh_ideal_length_factor_; }

  double find_timestep(const double max_timestep, const double cfl, const double fo, const double oh);

  void set_modified_time_ini(const double t) { modified_time_ini_ = t; }
  void set_max_timestep(const double t) { max_timestep_ = t; }
  void set_current_time(const double t) { changer_temps_courant(t); }
  void set_tstep_sauv(const int ts) { tstep_sauv_ = ts; }

  void check_stop_criteria(bool& stop) const;

protected:
  double dt_cfl_ = 1.e20, dt_fo_ = 1.e20, dt_oh_ = 1.e20, dt_fo_liq_ = 1.e20;
  double dt_fo_vap_ = 1.e20, dt_cfl_liq_ = 1.e20, dt_cfl_vap_ = 1.e20;
  double timestep_facsec_ = 1., max_simu_time_ = 1e6;
  double modified_time_ini_ = 0., max_timestep_ = -123.;
  double cfl_ = 1., fo_ = 1., oh_ = 1.;

  int enable_dt_oh_ideal_length_factor_ = 0, first_step_interface_smoothing_ = 0;
  int tstep_sauv_ = 0, tstep_init_ = 0, use_tstep_init_ = 0;
};

#endif /* Schema_Temps_IJK_base_included */
