/****************************************************************************
* Copyright (c) 2019, CEA
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
// File      : Force_ph.h
// Directory : $IJK_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef force_ph_included
#define force_ph_included


#include <iostream>
#include <string>
#include <vector>
#include <Force_sp.h>
#include <IJK_Field_vector.h>
// #include <fftw3.h>
#include <IJK_Splitting.h>
#include <communications.h>
#include <Objet_U.h>

#include <IJK_Field.h>

class Force_ph : public Objet_U
{

  Declare_instanciable( Force_ph ) ;

public:

  // void initialise(int nproc_tot, int ni, int nj, int nk, int nl,int nm,int nn,
  // 		double Lx, double Ly, double Lz, double Ox,double Oy,double Oz, int momin, int momax, double kmin, double kmax,
  // 	  std::string nom_fichier, const IJK_Splitting& splitting
  // 	);
  void initialise(int nproc_tot, int ni, int nj, int nk, int nl,int nm,int nn,
                  double Lx, double Ly, double Lz, double Ox,double Oy,double Oz, int momin, int momax, double kmin, double kmax,
                  std::string nom_fichier, const IJK_Splitting& splitting,
                  int a_i_offset, int a_j_offset, int a_k_offset
                 );
  void from_spect_to_phys2(const std:: vector <double >& coeff_force);
  void from_spect_to_phys_opti2(ArrOfDouble& coeff_force);
  void from_spect_to_phys_opti2_advection(ArrOfDouble& coeff_force, const ArrOfDouble& advection_length);
  void from_spect_to_phys_bis(const std::vector< std::vector< std:: vector <double >>>& coeff_force);
  void set_zero();
  void cheat_function();
  // void from_spect_to_phys_fftw(fftw_complex* in);
  void write(std::string nom_fichier_sortie, double t);
  void write_separate(std::string nom_fichier_sortie, double t);
  void write_offset_index_position( const IJK_Splitting& my_splitting);
  void compute_energie();
  double get_energie();
  IJK_Field_vector3_double get_force_attribute();
  IJK_Field_vector3_double& get_force_attribute2();

  void gbz_gather(IJK_Field_vector3_double force_);

private:
  int nproc_tot = 0;

  int ni = 0;
  int nj = 0;
  int nk = 0;
  int n_ijk = 0;
  int nl = 0;
  int nm = 0;
  int nn = 0;
  int n_lmn = 0;
  double Lx = 0.;
  double Ly = 0.;
  double Lz = 0.;
  double Ox = 0.;
  double Oy = 0.;
  double Oz = 0.;
  double kmin = 0.;
  double kmax = 0.;
  int momin = 0;
  int momax = 0;
  IJK_Field_vector3_double force_;
  std::vector<std::vector< std:: vector < double > > > force;
  double energie = 0.;

  int nproc_i = 0;
  int nproc_j = 0;
  int nproc_k = 0;
  int i_offset = 0;
  int j_offset = 0;
  int k_offset = 0;


};

std::vector< std::vector< std:: vector <double >>> set_dimensions(std::vector< std::vector< std:: vector <double >>> the_vector, int dim_one, int dim_two, int dim_three);

// IJK_Field_vector3_double set_to_zero(IJK_Field_vector3_double vector)
// {
// 		for (int dir=0; dir<3; dir++)
// 			for (int i=0; i<vector[0].ni(); i++)
// 				for (int j=0; j<vector[1].nj(); j++)
// 					for (int k=0; k<vector[2].nk(); k++)
// 						vector[dir](i,j,k) = 0.;
// 		return vector;
// }

#endif
