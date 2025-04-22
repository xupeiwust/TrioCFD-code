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

#ifndef FT_Field_included
#define FT_Field_included


#include <TRUST_Ref.h>
#include <TRUSTVect.h>
#include <Domaine_IJK.h>
//#include <Maillage_FT_IJK.h>
#include <Operator_FT_Disc.h>
#include <Descripteur_FT.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <variant>
#include <Process.h>
using namespace std;

/*! @brief : class FT_Field
 */
class Maillage_FT_IJK;

class Point2D
{
public:
  double x, y;
  static double tol;
  Point2D() : x(0), y(0) {}
  Point2D(double xVal, double yVal) : x(xVal), y(yVal) {}
  Point2D(const Point2D& pt): x(pt.x) , y(pt.y) {}

  void operator = (const Point2D& pt)
  {
    x = pt.x ;
    y = pt.y ;
  }
  double get_x()const
  {
    return x;
  }
  double get_y()const
  {
    return y;
  }
  /* void set(double xVal, double yVal)
   {
     x = xVal;
     y = yVal;
   }*/
  void print() const
  {
    std::cout << "Point2D(" << x << ", " << y << ")" << std::endl;
  }
  bool operator==(const Point2D& rhs)
  {
    return abs(x-rhs.x)<tol && abs(y-rhs.y)<tol;
  }
};


class Point3D
{
public:
  static double tol;
  double x, y, z;
  Point3D(const Point3D& pt): x(pt.x) , y(pt.y), z(pt.z) {}
  //Point3D(const Point3D pt): x(pt.x) , y(pt.y), z(pt.z) {}
  Point3D() : x(0), y(0), z(0) {}
  Point3D(double xVal, double yVal, double zVal) : x(xVal), y(yVal), z(zVal) {}

  double norme()
  {
    return std::sqrt(x*x + y*y + z*z) ;
  }
  void operator = (const Point3D& pt)
  {
    x = pt.x ;
    y = pt.y ;
    z = pt.z ;
  }

  Point3D operator-(const Point3D& other) const
  {
    return Point3D(x - other.x, y - other.y, z - other.z);
  }
  Point3D operator+(const Point3D& other) const
  {
    return Point3D(x + other.x, y + other.y, z + other.z);
  }
  Point3D operator/(const double& other) const
  {
    return Point3D(x/other, y/other, z/other);
  }
  double get_x() const
  {
    return x;
  }
  double get_y() const
  {
    return y;
  }
  double get_z() const
  {
    return z;
  }
  void print() const
  {
    std::cout << "Point3D(" << x << ", " << y <<  ", " << z << ")" << std::endl;
  }
  // Comparison operator for sorting Point3D
  bool operator==(const Point3D& rhs) const
  {
    return abs(x-rhs.x)<tol && abs(y-rhs.y)<tol && abs(z-rhs.z)<tol;
  }
  bool operator<(const Point3D& rhs) const
  {
    if (x != rhs.x) return x < rhs.x;
    if (y != rhs.y) return y < rhs.y;
    return z < rhs.z;
  }
};



class triangle
{
public:
  Point3D p1, p2, p3;
  //triangle() : p1, p2, p3 {}
  triangle(Point3D xVal, Point3D yVal, Point3D zVal) : p1(xVal), p2(yVal), p3(zVal) {}
  /*void set(Point3D xVal, Point3D yVal, Point3D zVal)
  {
    p1 = xVal;
    p2 = yVal;
    p3 = zVal;
  }*/
  bool operator==(const triangle& rhs)
  {
    return p1 == rhs.p1 && p2 == rhs.p2 && p3 == rhs.p3;
  }
};

class FT_Field : public Objet_U
{
  Declare_instanciable_sans_constructeur( FT_Field ) ;
private:
  int nb_sommet_change_ = 0;
  DoubleTab sommet_bouge_ ;

  //IntTab index_fa7_structure_ ;
  //IntTab facette_par_structure_;
  IntTab indice_global_to_local_final_;
  // sommet_change contient les indices des sommets principaux
  // et le nombre de fa7 initiales impliquees sur ce sommet
  // et le nombre de fa7 finales impliquees sur ce sommet
  DoubleTab triangle_initiaux_ ;
  DoubleTab normale_facette_initiale_;
  DoubleTab normale_facette_finale_;
  DoubleTab triangle_finaux_ ;
  DoubleTab Surfactant_facette_initiale_ ;
  IntTab indice_facette_finaux_ ;
  bool variable_intensive_=false;
  bool disable_surfactant_=true;
  int Taylor_test_ = 0;
  int disable_marangoni_source_term_ = 0;
  int print_debug_surfactant_ = 0;
  int only_remaillage_ = 0;
  int patch_conservation_surfactant_locale_ = 0;
  int patch_conservation_surfactant_globale_ = 0;
  int check_triangle_duplicata_ = 0;
  double Diff_coeff_surfactant_ = 0. ;
  double Concentration_surfactant_init_ = 0. ;
  double Surfactant_theoric_case_=0. ;


  double mean_surfactant_ = 0.;
  double mean_surface_ = 0.;
  DoubleTab facettes_sommets_full_compo_;
  DoubleTab liste_sommets_et_deplacements_;
  IntTab compo_transmises_a_envoyer_;
  ArrOfInt compo_connexe_sommets_;
  ArrOfInt proc_numero_;
  ArrOfInt proc_deja_echange_;
  ArrOfInt sorted_index_;
  int nb_compo_a_envoyer_max_ = 0;
  int nb_proc_echange_=0;
  int index_local_Ft_field_ = 0;

public:
  FT_Field();
  void initialize(const Maillage_FT_IJK& mesh, const DoubleTab& centre_mass);
  void update_gradient_laplacien_FT(const Maillage_FT_IJK& mesh);
  void update_sigma_grad_sigma(const Maillage_FT_IJK& mesh, const Domaine_IJK& splitting);
  void dimensionner_remaillage_FT_Field(Maillage_FT_IJK& mesh, const ArrOfIntFT& table_old_new);
  void sauvegarder_triangle(const Maillage_FT_IJK& mesh, const int i, const int avant_apres_remaillage);
  void echanger_triangles(Maillage_FT_IJK& mesh);
  void remailler_FT_Field(Maillage_FT_IJK& mesh);
  ArrOfDouble check_conservation(const Maillage_FT_IJK& mesh);
  void correction_conservation_globale(const Maillage_FT_IJK& mesh, const ArrOfDouble& surfactant_avant_remaillage, const ArrOfDouble& surfactant_apres_remaillage);
  ArrOfDouble FT_field_Array_;
  ArrOfDouble FT_field_Array_sommets_;
  DoubleTab  Grad_FT_field_Array_;
  ArrOfDouble Laplacian_FT_field_Array_;
  Operator_FT_Disc OpFTDisc_;
  double sigma0_;
  double R_;
  double T_;
  double Gamma_inf_;
  ArrOfDouble sigma_sommets_;
  ArrOfDouble sigma_facettes_;
  DoubleTab Grad_sigma_sommets_;
  void avancer_en_temps(const Maillage_FT_IJK& mesh, const double time_step);
  void passer_variable_extensive(const Maillage_FT_IJK& mesh);
  void passer_variable_intensive(const Maillage_FT_IJK& mesh);
  void preparer_tableau_avant_transport();
  void update_tableau_apres_transport();
  void nettoyer_espace_virtuel_facette(const Maillage_FT_IJK& mesh);
  void set_field_facettes(ArrOfDouble field);
  void set_field_sommets(ArrOfDouble field);
  void update_Field_sommets(const Maillage_FT_IJK& FTmesh, const ArrOfDouble& Field_facettes, ArrOfDouble& field_sommet);

  void exchange_data(int pe_send_, /* processor to send data */
                     int pe_recv_ /* processor to received data */,
                     const Maillage_FT_IJK& mesh, const DoubleTab& liste_sommets_avant_deplacement, const DoubleTab& liste_sommets_apres_deplacement, const ArrOfInt& compo_connexe_sommets_deplace);

  void champ_sommet_from_facettes(const ArrOfInt& compo_connexe_facettes, const Maillage_FT_IJK& mesh);

  void exchange_compo_connexe(int pe_send_,int pe_recv_, const Maillage_FT_IJK& mesh);
  void update_FT_Field_local_from_full_compo(const Maillage_FT_IJK& mesh);
  void completer_compo_connexe_partielle(const Maillage_FT_IJK& mesh, const Domaine_IJK& splitting, const DoubleTab& liste_sommets_apres_deplacement, const DoubleTab& liste_sommets_avant_deplacement, const ArrOfInt& compo_connexe_sommets_deplace);

  bool sauv_num_pe_echange(int pe);
  bool is_compo_in_proc(const int compo_connexe, const int pe_send);
  void calculer_volume_bulles(ArrOfDouble& volumes, DoubleTab& centre_gravite, const Maillage_FT_IJK& mesh) const;
  // Equality operator for Point3D
  /*bool operator==(const Point3D& lhs, const Point3D& rhs)
  {
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
  }

  // Comparison operator for sorting Point3D
  bool operator<(const Point3D& lhs, const Point3D& rhs)
  {
    if (lhs.x != rhs.x) return lhs.x < rhs.x;
    if (lhs.y != rhs.y) return lhs.y < rhs.y;
    return lhs.z < rhs.z;
  }*/

  void sortAndTrackIndices(const std::vector<double>& arr, std::vector<size_t>& indices);
  double det(const Point2D& a, const Point2D& b, const Point2D& c);
  bool isPointInTriangle(const Point2D& pt, const Point2D& v1, const Point2D& v2, const Point2D& v3);
  void Calculate_Facette_Intersection_Area(DoubleTab& Surface_fa7init, DoubleTab& Surface_fa7fin, DoubleTab& Surface_intersection, vector<Point3D> points_fa7_originale, vector<Point3D> points_fa7_finale, IntTab points_triangle_originaux, IntTab points_triangle_finaux, IntTab& normale_triangle_originaux, IntTab& normale_triangle_finaux);
  bool lineIntersection(const Point2D& a, const Point2D& b, const Point2D& c, const Point2D& d, Point2D& intersection);
  double polygonArea(const std::vector<Point2D>& vertices);
  double intersectionArea(Point2D t1[3], Point2D t2[3]);
  double norme(const Point3D& pt);
  std::vector<Point3D> removeDuplicates(std::vector<Point3D>& points);
  Point3D computeCentroid(const vector<Point3D>& points);
  void computeCovarianceMatrix(const vector<Point3D>& points, const Point3D& centroid, double cov[3][3]);
  void powerIteration(const double cov[3][3], double eigenVector[3], double& eigenValue);
  Point2D projectPointToPlane(const Point3D& point, const Point3D& centroid, const array<double, 3>& eigenVector1, const array<double, 3>& eigenVector2);
  int orientation_triangle(const Point3D& normale, const array<double, 3>& eigenVector1, const array<double, 3>& eigenVector2);
  vector<pair<double, array<double, 3>>> Main_2D_plane_eigenvectors(vector<Point3D> points);
  double triangleArea(const Point2D& p1, const Point2D& p2, const Point2D& p3);
  Point3D crossProduct(const Point3D& u, const Point3D& v);
  double scalarProduct(const Point3D& u, const Point3D& v);
  double magnitude(const Point3D& v);
  double triangleArea3D(const Point3D& A, const Point3D& B, const Point3D& C);
  Point3D calculer_normale_apres_deplacement(const int fa, const int somfa7, const Vecteur3 pos_apres_dep);



  void copy_FT_Field(FT_Field copy)
  {
    FT_field_Array_.copy_array(copy.FT_field_Array_);
    FT_field_Array_sommets_.copy_array(copy.FT_field_Array_sommets_);
    Grad_FT_field_Array_.copy_array(copy.Grad_FT_field_Array_);
    Laplacian_FT_field_Array_.copy_array(copy.Laplacian_FT_field_Array_);
    sigma_facettes_.copy_array(copy.sigma_facettes_);
    Grad_sigma_sommets_.copy_array(copy.Grad_sigma_sommets_);
    copy.OpFTDisc_=OpFTDisc_;
    copy.disable_surfactant_=disable_surfactant_;
  };
  void inject_array(const FT_Field& source, int nb_elements, int first_element_dest, int first_element_source)
  {
    FT_field_Array_.inject_array(source.FT_field_Array_, nb_elements, /* tous les elements de la source */
                                 first_element_dest, /* indice destination */
                                 first_element_source); /* indice source */
  };
  void resize_array(int index)
  {
    FT_field_Array_.resize_array(index);
  }
  void resize(int index)
  {
    FT_field_Array_.resize(index);
  }
  double& operator()(int index)
  {
    return FT_field_Array_(index);
  }
  double& operator[](int index)
  {
    return FT_field_Array_[index];
  }

  void echange_espace_virtuel(const Maillage_FT_Disc& mesh);

  int size() const
  {
    return FT_field_Array_.size_array();
  };
  int size_array() const
  {
    return FT_field_Array_.size_array();
  };

  int size_sommets() const
  {
    return FT_field_Array_sommets_.size_array();
  };

  int get_only_remaillage() const
  {
    return only_remaillage_;
  }

  Operator_FT_Disc get_OpFTDisc() const
  {
    return OpFTDisc_;
  };
  ArrOfDouble get_FT_field_Array() const
  {
    return FT_field_Array_;
  };
  bool get_disable_surfactant() const
  {
    return disable_surfactant_;
  };
  int get_disable_marangoni_source_term() const
  {
    return disable_marangoni_source_term_;
  };

  void set_disable_surfactant(bool disable_surfactant)
  {
    disable_surfactant_ = disable_surfactant;
  };

  ArrOfInt get_compo_connexe_sommets() const
  {
    return compo_connexe_sommets_;
  };
  ArrOfDouble get_FT_field_Array_non_const()
  {
    return FT_field_Array_;
  };
  ArrOfDouble get_FT_field_Array_sommets() const
  {
    return FT_field_Array_sommets_;
  };
  ArrOfDouble get_Grad_FT_field_Array(int dir) const
  {
    int nbsom = Grad_FT_field_Array_.dimension(0);
    ArrOfDouble  Grad_dir_FT_field_Array;
    Grad_dir_FT_field_Array.resize(nbsom);
    for (int som=0 ; som<nbsom ; som++)
      {
        Grad_dir_FT_field_Array(som)=Grad_FT_field_Array_(som, dir);
      }
    return Grad_dir_FT_field_Array;
  };
  ArrOfDouble get_sigma_facettes() const
  {
    return sigma_facettes_;
  };
  ArrOfDouble get_sigma_sommets() const
  {
    return sigma_sommets_;
  };
  DoubleTab get_grad_sigma_sommets() const
  {
    return Grad_sigma_sommets_;
  };
  ArrOfDouble get_grad_sigma_sommets(int dir) const
  {
    int nbsom = Grad_sigma_sommets_.dimension(0);
    ArrOfDouble  Grad_dir;
    Grad_dir.resize(nbsom);
    for (int som=0 ; som<nbsom ; som++)
      {
        Grad_dir(som)=Grad_sigma_sommets_(som, dir);
      }
    return Grad_dir;
  };

  ArrOfDouble get_grad_sigma_sommets_non_const(int dir)
  {
    int nbsom = Grad_sigma_sommets_.dimension(0);
    ArrOfDouble  Grad_dir;
    Grad_dir.resize(nbsom);
    for (int som=0 ; som<nbsom ; som++)
      {
        Grad_dir(som)=Grad_sigma_sommets_(som, dir);
      }
    return Grad_dir;
  };
  ArrOfDouble get_Laplacian_FT_field_Array() const
  {
    return Laplacian_FT_field_Array_;
  };
  DoubleTab get_sommet_bouge() const
  {
    return sommet_bouge_;
  };
  DoubleTab& get_facettes_sommets_full_compo_non_const()
  {
    return facettes_sommets_full_compo_;
  };
  DoubleTab& get_liste_sommets_et_deplacements_non_const()
  {
    return liste_sommets_et_deplacements_;
  };
  const ArrOfInt get_sorted_index()
  {
    return sorted_index_;
  };

  int get_nb_triangle_finaux()
  {
    return triangle_finaux_.dimension(0);
  }
  int get_nb_triangle_initiaux()
  {
    return triangle_initiaux_.dimension(0);
  }
  void reinit_remeshing_table()
  {
    normale_facette_initiale_.resize(0, 3, RESIZE_OPTIONS::NOCOPY_NOINIT);
    normale_facette_finale_.resize(0, 3, RESIZE_OPTIONS::NOCOPY_NOINIT);
    triangle_initiaux_.resize(0, 3, 3, RESIZE_OPTIONS::NOCOPY_NOINIT);
    triangle_finaux_.resize(0, 3, 3, RESIZE_OPTIONS::NOCOPY_NOINIT);
    Surfactant_facette_initiale_.resize(0, RESIZE_OPTIONS::NOCOPY_NOINIT);
    indice_facette_finaux_.resize(0, RESIZE_OPTIONS::NOCOPY_NOINIT);
  }

  void set_tolerance_point_identique(double longueur_cara_fa7)
  {
    Point3D::tol = 1.e-3*longueur_cara_fa7;
    Point2D::tol = 1.e-3*longueur_cara_fa7;
  }
  double get_tolerance_point_identique()
  {
    return Point3D::tol;
  }

};

#endif /* FT_Field_included */
