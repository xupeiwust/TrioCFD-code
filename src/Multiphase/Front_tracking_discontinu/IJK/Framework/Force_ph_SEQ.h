// #ifndef DEF_FORCE_PH
// #define DEF_FORCE_PH
//
// #include <iostream>
// #include <string>
// #include <vector>
// #include <Force_sp.h>
// #include <IJK_Field_vector.h>
// #include <IJK_Field.h>
// #include <fftw3.h>
// #include <communications.h>
//
// class Force_ph
// {
// 	public:
//
// 	Force_ph();
// 	~Force_ph();
// 	void initialise(int nproc_tot, int ni, int nj, int nk, int nl,int nm,int nn,
// 			double Lx, double Ly, double Lz, double Ox,double Oy,double Oz, int momin, int momax, double kmin, double kmax,
// 		  std::string nom_fichier, const Domaine_IJK& splitting);
// 	// void initialise(int ni, int nj, int nk, int nl,int nm,int nn,
// 	// 		double Lx, double Ly, double Lz, double kmin, double kmax);
// 	void from_spect_to_phys(const std::vector< std::vector< std:: vector <double >>>& coeff_force);
// 	void from_spect_to_phys2(const std:: vector <double >& coeff_force);
// 	void from_spect_to_phys_opti(const std::vector< std::vector< std:: vector <double >>>& coeff_force);
// 	void from_spect_to_phys_opti2(const std:: vector <double >& coeff_force);
// 	void from_spect_to_phys_bis(const std::vector< std::vector< std:: vector <double >>>& coeff_force);
// 	void cheat_function();
// 	// void from_spect_to_phys_fftw(fftw_complex* in);
// 	void write(std::string nom_fichier_sortie, double t);
// 	void write_separate(std::string nom_fichier_sortie, double t);
// 	void write_offset_index_position( const Domaine_IJK& my_splitting);
// 	void compute_energie();
// 	double get_energie();
// 	IJK_Field_vector3_double get_force_attribute();
//
// 	void gbz_gather(IJK_Field_vector3_double force_);
//
// 	private:
// 		int nproc_tot;
//
// 	int ni,nj,nk,n_ijk;
// 	int nl,nm,nn,n_lmn;
// 	double Lx, Ly, Lz;
// 	double Ox, Oy, Oz;
// 	double kmin,kmax;
// 	int momin,momax;
// 	IJK_Field_vector3_double force_;
// 	std::vector<std::vector< std:: vector < double > > > force;
// 	double energie;
//
// 	int nproc_i,nproc_j,nproc_k;
//
//
// };
//
// std::vector< std::vector< std:: vector <double >>> set_dimensions(std::vector< std::vector< std:: vector <double >>> the_vector, int dim_one, int dim_two, int dim_three);
//
// // IJK_Field_vector3_double set_to_zero(IJK_Field_vector3_double vector)
// // {
// // 		for (int dir=0; dir<3; dir++)
// // 			for (int i=0; i<vector[0].ni(); i++)
// // 				for (int j=0; j<vector[1].nj(); j++)
// // 					for (int k=0; k<vector[2].nk(); k++)
// // 						vector[dir](i,j,k) = 0.;
// // 		return vector;
// // }
//
// #endif
