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

//#include <Maillage_FT_IJK.h>

#include <FT_Field.h>
#include <Maillage_FT_IJK.h>
#include <Param.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <variant>
#include <IJK_communications.h>
#include <Process.h>
using namespace std;

Implemente_instanciable_sans_constructeur( FT_Field, "FT_Field", Objet_U ) ;

double Point3D::tol = 1.e-10;
double Point2D::tol = 1.e-10;

void FT_Field::initialize(const Maillage_FT_IJK& mesh, const DoubleTab& centre_mass)
{
  if (Surfactant_theoric_case_==0. and Concentration_surfactant_init_==0.)
    {
      disable_surfactant_ = true;
      return;
    }
  else
    {
      disable_surfactant_ = false;
    }
  int nbfa7=mesh.nb_facettes();
  FT_field_Array_.resize(nbfa7);
  mean_surfactant_= 0.;
  mean_surface_ = 0.;

  const DoubleTab& sommets=mesh.sommets();
  const ArrOfDouble& Sfa7 = mesh.get_update_surface_facettes();
  const IntTab& facettes=mesh.facettes();
  const ArrOfInt& compo_connex = mesh.compo_connexe_facettes();


  if (Surfactant_theoric_case_==1.)
    {
      // on initialise un champ pour solution analytique si option activee
      for (int fa7=0 ; fa7<nbfa7 ; fa7++)
        {
          int compo = compo_connex(fa7);
          Point3D xg_bulle = {centre_mass(compo, 0), centre_mass(compo, 1), centre_mass(compo, 2)};
          int indice_sommet1 = facettes(fa7,0);
          int indice_sommet2 = facettes(fa7,1);
          int indice_sommet3 = facettes(fa7,2);
          Point3D x0 = {sommets(indice_sommet1,0),sommets(indice_sommet1,1),sommets(indice_sommet1,2)};
          Point3D x1 = {sommets(indice_sommet2,0),sommets(indice_sommet2,1),sommets(indice_sommet2,2)};
          Point3D x2 = {sommets(indice_sommet3,0),sommets(indice_sommet3,1),sommets(indice_sommet3,2)};
          Point3D xg_fa7 = (x0 + x1 + x2)/3.;
          Point3D AB = xg_fa7-xg_bulle;

          double costheta = AB.x/norme(AB);
          FT_field_Array_[fa7]= 0.5 * (1.-costheta) ;
          mean_surfactant_+=FT_field_Array_(fa7)*Sfa7(fa7);
          mean_surface_+=Sfa7(fa7);
        }

    }
  else
    {
      // sinon on initialise avec la concentration initiale renseignee dans le jdd
      for (int fa7=0 ; fa7<nbfa7 ; fa7++)
        {
          FT_field_Array_[fa7]= Concentration_surfactant_init_ ;
          mean_surfactant_+=FT_field_Array_(fa7)*Sfa7(fa7);
          mean_surface_+=Sfa7(fa7);
        }
    }
  mean_surfactant_=Process::mp_sum(mean_surfactant_);
  mean_surface_=Process::mp_sum(mean_surface_);
  mean_surfactant_/=mean_surface_;
  Concentration_surfactant_init_=mean_surfactant_;
  mesh.desc_facettes().echange_espace_virtuel(FT_field_Array_);
  update_Field_sommets(mesh, FT_field_Array_, FT_field_Array_sommets_);
}

FT_Field::FT_Field()
{
}


Entree& FT_Field::readOn( Entree& is )
{
  // Objet_U::readOn( is );
  print_debug_surfactant_=0;
  only_remaillage_=0;
  patch_conservation_surfactant_locale_=0;
  patch_conservation_surfactant_globale_=0;
  check_triangle_duplicata_=0;
  Diff_coeff_surfactant_=0.;
  Taylor_test_=0 ;
  disable_marangoni_source_term_ = 0;
  sigma0_ = 0.;
  R_= 8.314; // constante des gaz parfaits
  T_ = 290 ; // temperature absolue
  Gamma_inf_ = 0. ;
  Param param(que_suis_je());
  param.ajouter("print_debug_surfactant", &print_debug_surfactant_);
  param.ajouter("only_remaillage", &only_remaillage_);
  param.ajouter("patch_conservation_surfactant_locale", &patch_conservation_surfactant_locale_);
  param.ajouter("patch_conservation_surfactant_globale", &patch_conservation_surfactant_globale_);
  param.ajouter("check_triangle_duplicata", &check_triangle_duplicata_);
  param.ajouter("Diff_coeff_surfactant", &Diff_coeff_surfactant_);
  param.ajouter("Surfactant_theoric_case", &Surfactant_theoric_case_);
  param.ajouter("Concentration_surfactant_init", &Concentration_surfactant_init_);
  param.ajouter("sigma0", &sigma0_);
  param.ajouter("Taylor_test", &Taylor_test_);
  param.ajouter("disable_marangoni_source_term", &disable_marangoni_source_term_);
  param.ajouter("R", &R_);
  param.ajouter("T", &T_);
  param.ajouter("Gamma_inf", &Gamma_inf_);
  param.lire_avec_accolades(is);

  if (Surfactant_theoric_case_==0. and Concentration_surfactant_init_==0.)
    {
      disable_surfactant_ = true;
    }
  else
    {
      disable_surfactant_ = false;
    }
  return is;
}

Sortie& FT_Field::printOn( Sortie& os ) const
{
  //Objet_U::printOn(os);
  os  << "{\n"
      << "   print_debug_surfactant " << print_debug_surfactant_ << "\n"
      << "   only_remaillage " << only_remaillage_ << "\n"
      << "   patch_conservation_surfactant_locale " << patch_conservation_surfactant_locale_ << "\n"
      << "   patch_conservation_surfactant_globale " << patch_conservation_surfactant_globale_ << "\n"
      << "   Diff_coeff_surfactant " << Diff_coeff_surfactant_ << "\n"
      << "   Surfactant_theoric_case " << Surfactant_theoric_case_ << "\n"
      << "   Concentration_surfactant_init " << Concentration_surfactant_init_ << "\n";
  return os;
}


void FT_Field::avancer_en_temps(const Maillage_FT_IJK& mesh, const double time_step)
{
// au pas de temps N+1, je vais calculer S_{n+1}gamma_{n+1}=S_{n}gamma_{n}+DT*Dgamma_{n+1}+DT*Source
// S_{n}gamma_{n} est la version intensive de la variable de concentration de Surfactant


  if (!variable_intensive_)
    {
      std::cout << "la variable doit etre intensive a l'appel de cette fonction" << std::endl;
      Process::exit();
    }
  // on met a jour le Laplacien interfacial apres transport
  passer_variable_extensive(mesh);
  update_gradient_laplacien_FT(mesh);
  passer_variable_intensive(mesh);

  // ajout du Laplacien
  for (int fa7=0 ; fa7<FT_field_Array_.size_array() ; fa7++)
    {
      if(! mesh.facette_virtuelle(fa7))
        {
          FT_field_Array_[fa7]+=time_step*Diff_coeff_surfactant_*Laplacian_FT_field_Array_(fa7);
        }
    }
  mesh.desc_facettes().echange_espace_virtuel(FT_field_Array_);
  // TODO :: ajout du terme source pour les cas surfactant soluble


}

void FT_Field::update_gradient_laplacien_FT(const Maillage_FT_IJK& mesh)
{
  OpFTDisc_.Operator_Laplacian_FT_element(FT_field_Array_,mesh, Laplacian_FT_field_Array_, Grad_FT_field_Array_);
}

void FT_Field::update_sigma_grad_sigma(const Maillage_FT_IJK& mesh, const Domaine_IJK& splitting)
{
  // on calcule les variations de sigma associees aux variations de FT_field_Array_
  const int nbsom=mesh.nb_sommets();
  const int nbfa7=mesh.nb_facettes();
  sigma_sommets_.resize(nbsom);
  sigma_facettes_.resize(nbfa7);

  // on exprime la tension de surface aux facettes, comme pour FT_field_Array_
  if (!Taylor_test_)
    {
      for (int fa=0 ; fa<nbfa7 ; fa++)
        {
          if(! mesh.facette_virtuelle(fa))
            {
              sigma_facettes_[fa]= max(1.e-16,sigma0_ - R_ * T_ * FT_field_Array_[fa]);
              // max(1.e-5,sigma0_+R_*T_*Gamma_inf_*log(max(1.-FT_field_Array_[fa]/Gamma_inf_, 1.e-5)));
            }
        }
    }
  else
    {
      // on teste le cas analytique de Young (voir Muradoglu & Trygvason 2008)
      const DoubleTab& sommets=mesh.sommets();
      const IntTab& facettes=mesh.facettes();
      const ArrOfInt& compo_connex = mesh.compo_connexe_facettes();
      ArrOfDouble volume_reel;
      DoubleTab position;
      calculer_volume_bulles(volume_reel, position, mesh);
      double Lx =  splitting.get_domain_length(0);
      for (int fa=0 ; fa<nbfa7 ; fa++)
        {
          if(! mesh.facette_virtuelle(fa))
            {
              int compo = compo_connex(fa);
              Point3D xg_bulle = {position(compo, 0), position(compo, 1), position(compo, 2)};
              int indice_sommet1 = facettes(fa,0);
              int indice_sommet2 = facettes(fa,1);
              int indice_sommet3 = facettes(fa,2);
              Point3D x0 = {sommets(indice_sommet1,0),sommets(indice_sommet1,1),sommets(indice_sommet1,2)};
              Point3D x1 = {sommets(indice_sommet2,0),sommets(indice_sommet2,1),sommets(indice_sommet2,2)};
              Point3D x2 = {sommets(indice_sommet3,0),sommets(indice_sommet3,1),sommets(indice_sommet3,2)};
              Point3D xg_fa7 = (x0 + x1 + x2)/3.;
              // il faut faire comme si la bulle etait fixe au centre du domaine (voir these kalyani)
              double x_centre = Lx/2. + (xg_fa7.x-xg_bulle.x);
              // on prend beta = 1.
              //double pos = std::fmod(std::fmod(pos_ref + offset - decallage_bulle_reel_ext_domaine_reel, Lx) + Lx, Lx) + decallage_bulle_reel_ext_domaine_reel;
              sigma_facettes_[fa]= sigma0_ *(1.- 1. * x_centre/Lx);
            }
        }
    }
  mesh.desc_facettes().echange_espace_virtuel(sigma_facettes_);
  // calcule du gradient de sigma
  OpFTDisc_.Operator_Gradient_FT_sommets(sigma_facettes_, mesh, Grad_sigma_sommets_);//, true);

  // on extrapole la valeur de la tension de surface aux sommets pour l'algo de IJK_Interfaces_
  update_Field_sommets(mesh, sigma_facettes_, sigma_sommets_);
}



void FT_Field::passer_variable_intensive(const Maillage_FT_IJK& mesh)
{
  // Le maillage ne doit pas avoir de doublon de facettes a l'appel de cette fonction
  // Utilisation de nettoyer maillage en amont recommandé
  if (!variable_intensive_)
    {
      const ArrOfDouble& Sfa7 = mesh.get_update_surface_facettes();
      for (int fa7=0 ; fa7<FT_field_Array_.size_array() ; fa7++)
        {
          if(! mesh.facette_virtuelle(fa7))
            {
              FT_field_Array_[fa7]*=Sfa7(fa7);
            }
        }
      mesh.desc_facettes().echange_espace_virtuel(FT_field_Array_);
      variable_intensive_=true;
    }
}
void FT_Field::passer_variable_extensive(const Maillage_FT_IJK& mesh)
{
  // Le maillage ne doit pas avoir de doublon de facettes a l'appel de cette fonction
  // Utilisation de nettoyer maillage en amont recommandé
  if (variable_intensive_)
    {
      const ArrOfDouble& Sfa7 = mesh.get_update_surface_facettes();
      for (int fa7=0 ; fa7<FT_field_Array_.size_array() ; fa7++)
        {
          if(! mesh.facette_virtuelle(fa7))
            {
              if (Sfa7(fa7)!=0.)
                {
                  FT_field_Array_[fa7]/=Sfa7(fa7);
                }
              else
                {
                  FT_field_Array_[fa7] = 0.;
                }
            }
        }
      mesh.desc_facettes().echange_espace_virtuel(FT_field_Array_);
      variable_intensive_=false;
    }
}


void FT_Field::nettoyer_espace_virtuel_facette(const Maillage_FT_IJK& mesh)
{
  // On decale toutes les facettes reelles pour remplir lespace libere par les facettes virtuelles
  {
    const int nbfacettes = mesh.facettes().dimension(0);
    int n = 0;
    int i;
    for (i = 0; i < nbfacettes; i++)
      {
        const int virtuelle = mesh.facette_virtuelle(i);
        if (!virtuelle)
          {
            FT_field_Array_(n) = FT_field_Array_(i);
            n++;
          }
      }
    FT_field_Array_.resize(n);
  }
}


void FT_Field::update_Field_sommets(const Maillage_FT_IJK& FTmesh, const ArrOfDouble& Field_facettes, ArrOfDouble& field_sommet)
{
  const int nbfa7=FTmesh.nb_facettes();
  const int nbsom=FTmesh.nb_sommets();
  const IntTab& facettes=FTmesh.facettes();
  const ArrOfDouble& Sfa7 = FTmesh.get_surface_facettes();

  ArrOfDouble Surface_sommet;
  field_sommet.resize(nbsom);
  Surface_sommet.resize(nbsom);

  for (int som=0 ; som<nbsom ; som++)
    {
      field_sommet[som]=0.;
      Surface_sommet[som]=0.;
    }

  // interpolation des valeur de Phi et de la normale aux sommets en moyennant les contributions adjacentes

  for (int fa7=0 ; fa7<nbfa7 ; fa7++)
    {
      if(! FTmesh.facette_virtuelle(fa7))
        {
          for (int sommet_fa7=0 ; sommet_fa7<3 ; sommet_fa7++)
            {
              int indice_sommet = facettes(fa7,sommet_fa7);
              field_sommet[indice_sommet]+=Field_facettes[fa7]*Sfa7[fa7];
              Surface_sommet[indice_sommet]+=Sfa7[fa7];
            }
        }
    }
  // On a calcule la contribution de chaque facette reelle aux differents sommets.
  // Certaines contributions ont ete ajoutees a des sommets virtuels, il
  // faut recuperer ces contributions sur le sommet reel.
  const Desc_Structure_FT& desc_sommets = FTmesh.desc_sommets();
  desc_sommets.collecter_espace_virtuel(field_sommet, MD_Vector_tools::EV_SOMME);
  desc_sommets.collecter_espace_virtuel(Surface_sommet, MD_Vector_tools::EV_SOMME);

  for (int som=0 ; som<nbsom ; som++)
    {
      if(! FTmesh.sommet_virtuel(som))
        {
          field_sommet[som]/= Surface_sommet[som];
        }
      Surface_sommet[som] = 0.;
    }
  // les sommets reels sont mis a jour
  // il reste a mettre a jour les sommets virtuels
  desc_sommets.echange_espace_virtuel(field_sommet);
}

void FT_Field::echange_espace_virtuel(const Maillage_FT_Disc& mesh)
{
  mesh.desc_facettes().echange_espace_virtuel(FT_field_Array_);
}




// Toute les fonction qui suivent servent a calculer les intersections entre deux triangles dans un espace 2D.
// Cela sert a redistribuer de maniere conservative et en diffusant le moins possible la quantite de surfactant
// lorsque lon supprime les petites cellules.

// Comparator function to sort indices based on the values in the array

bool compare(const std::pair<size_t, double>& a, const std::pair<size_t, double>& b)
{
  return a.second < b.second;
}

void FT_Field::sortAndTrackIndices(const std::vector<double>& arr, std::vector<size_t>& indices)
{
  size_t n = arr.size();

  // Create a vector of pairs where the first element is the index and the second is the value
  std::vector<std::pair<size_t, double>> indexedArray(n);
  for (size_t i = 0; i < n; ++i)
    {
      indexedArray[i] = std::make_pair(i, arr[i]);
    }

  // Sort the indexedArray using the custom comparator
  std::sort(indexedArray.begin(), indexedArray.end(), compare);

  // Extract the sorted indices
  for (size_t i = 0; i < n; ++i)
    {
      indices[i] = indexedArray[i].first;
    }
}


// Function to calculate the determinant
double FT_Field::det(const Point2D& a, const Point2D& b, const Point2D& c)
{
  return a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y);
}

// Function to check if a point is inside a triangle
bool FT_Field::isPointInTriangle(const Point2D& pt, const Point2D& v1, const Point2D& v2, const Point2D& v3)
{
  double d1 = det(pt, v1, v2);
  double d2 = det(pt, v2, v3);
  double d3 = det(pt, v3, v1);
  bool has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
  bool has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);
  return !(has_neg && has_pos);
}

// Function to find the intersection of two lines
bool FT_Field::lineIntersection(const Point2D& a, const Point2D& b, const Point2D& c, const Point2D& d, Point2D& intersection)
{
  double a1 = b.y - a.y;
  double b1 = a.x - b.x;
  double c1 = a1 * a.x + b1 * a.y;

  double a2 = d.y - c.y;
  double b2 = c.x - d.x;
  double c2 = a2 * c.x + b2 * c.y;

  double determinant = a1 * b2 - a2 * b1;

  if (std::abs(determinant) < 1e-10)
    {
      return false; // The lines are parallel
    }

  intersection.x = (b2 * c1 - b1 * c2) / determinant;
  intersection.y = (a1 * c2 - a2 * c1) / determinant;

  if (std::min(a.x, b.x) <= intersection.x && intersection.x <= std::max(a.x, b.x) &&
      std::min(a.y, b.y) <= intersection.y && intersection.y <= std::max(a.y, b.y) &&
      std::min(c.x, d.x) <= intersection.x && intersection.x <= std::max(c.x, d.x) &&
      std::min(c.y, d.y) <= intersection.y && intersection.y <= std::max(c.y, d.y))
    {
      return true;
    }

  return false;
}

// Function to calculate the area using the Shoelace formula
double FT_Field::polygonArea(const std::vector<Point2D>& vertices)
{
  double area = 0;
  int n = static_cast<int>(vertices.size());

  for (int i = 0; i < n; ++i)
    {
      area += (vertices[i].x * vertices[(i + 1) % n].y) - (vertices[(i + 1) % n].x * vertices[i].y);
    }
  return std::abs(area) / 2.0;
}

double FT_Field::intersectionArea(Point2D t1[3], Point2D t2[3])
{
  vector<Point2D> intersectionPoints;

  // Find intersection points of edges
  for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
        {
          Point2D intersect;
          if (lineIntersection(t1[i], t1[(i + 1) % 3], t2[j], t2[(j + 1) % 3], intersect))
            {
              intersectionPoints.push_back(intersect);
            }
        }
    }

  // Add vertices of t1 that are inside t2
  for (int i = 0; i < 3; ++i)
    {
      if (isPointInTriangle(t1[i], t2[0], t2[1], t2[2]))
        {
          intersectionPoints.push_back(t1[i]);
        }
    }

  // Add vertices of t2 that are inside t1
  for (int i = 0; i < 3; ++i)
    {
      if (isPointInTriangle(t2[i], t1[0], t1[1], t1[2]))
        {
          intersectionPoints.push_back(t2[i]);
        }
    }

  // Remove duplicate points
  sort(intersectionPoints.begin(), intersectionPoints.end(), [](const Point2D &a, const Point2D &b)
  {
    return a.x == b.x ? a.y < b.y : a.x < b.x;
  });
  intersectionPoints.erase(unique(intersectionPoints.begin(), intersectionPoints.end(), [](const Point2D &a, const Point2D &b)
  {
    return std::abs(a.x - b.x) < 1e-10 && std::abs(a.y - b.y) < 1e-10;
  }), intersectionPoints.end());

  // Sort points to form a closed polygon
  if (!intersectionPoints.empty())
    {
      Point2D centroid = {0, 0};
      for (const auto &pt : intersectionPoints)
        {
          centroid.x += pt.x;
          centroid.y += pt.y;
        }
      centroid.x /= static_cast<double>(intersectionPoints.size());
      centroid.y /= static_cast<double>(intersectionPoints.size());

      sort(intersectionPoints.begin(), intersectionPoints.end(), [&centroid](const Point2D &a, const Point2D &b)
      {
        double angleA = atan2(a.y - centroid.y, a.x - centroid.x);
        double angleB = atan2(b.y - centroid.y, b.x - centroid.x);
        return angleA < angleB;
      });
    }

  // Calculate the area of the intersection polygon
  double intersectionArea = polygonArea(intersectionPoints);

  return intersectionArea;
}


double FT_Field::norme(const Point3D& pt)
{
  return std::sqrt(pt.x*pt.x + pt.y*pt.y + pt.z*pt.z) ;
}
// Function to remove duplicates from a vector of Point3D
std::vector<Point3D> FT_Field::removeDuplicates(std::vector<Point3D>& points)
{
  // Sort the vector of points
  std::sort(points.begin(), points.end());
  // Erase duplicates using unique algorithm
  points.erase(std::unique(points.begin(), points.end()), points.end());
  return points;
}



// Function to compute the centroid of the points
Point3D FT_Field::computeCentroid(const vector<Point3D>& points)
{
  Point3D centroid = {0.0, 0.0, 0.0};

  for (const auto &p : points)
    {
      centroid.x += p.x;
      centroid.y += p.y;
      centroid.z += p.z;
    }

  centroid.x /= static_cast<double>(points.size());
  centroid.y /= static_cast<double>(points.size());
  centroid.z /= static_cast<double>(points.size());

  return centroid;
}

// Function to compute the covariance matrix
void FT_Field::computeCovarianceMatrix(const vector<Point3D>& points, const Point3D& centroid, double cov[3][3])
{
  // Initialize covariance matrix to zero
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      cov[i][j] = 0.0;

  for (const auto &p : points)
    {
      double x = p.x - centroid.x;
      double y = p.y - centroid.y;
      double z = p.z - centroid.z;
      cov[0][0] += x * x;
      cov[0][1] += x * y;
      cov[0][2] += x * z;
      cov[1][0] += y * x;
      cov[1][1] += y * y;
      cov[1][2] += y * z;
      cov[2][0] += z * x;
      cov[2][1] += z * y;
      cov[2][2] += z * z;
    }

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      cov[i][j] /=  static_cast<double>(points.size());
}

// Function to compute eigenvalues and eigenvectors of a 3x3 matrix using the power iteration method
void FT_Field::powerIteration(const double cov[3][3], double eigenVector[3], double& eigenValue)
{
  double b_k[3] = {1.0, 1.0, 1.0};
  double b_k1[3];

  const int maxIterations = 100;
  const double tolerance = 1e-10;

  for (int iter = 0; iter < maxIterations; ++iter)
    {
      // Multiply cov by b_k: b_k1 = cov * b_k
      for (int i = 0; i < 3; ++i)
        {
          b_k1[i] = 0.0;
          for (int j = 0; j < 3; ++j)
            {
              b_k1[i] += cov[i][j] * b_k[j];
            }
        }

      // Normalize b_k1
      double norm = sqrt(b_k1[0] * b_k1[0] + b_k1[1] * b_k1[1] + b_k1[2] * b_k1[2]);
      for (int i = 0; i < 3; ++i)
        {
          b_k1[i] /= norm;
        }

      // Check for convergence
      double diff = sqrt((b_k1[0] - b_k[0]) * (b_k1[0] - b_k[0]) +
                         (b_k1[1] - b_k[1]) * (b_k1[1] - b_k[1]) +
                         (b_k1[2] - b_k[2]) * (b_k1[2] - b_k[2]));
      if (diff < tolerance)
        {
          break;
        }

      // Update b_k
      for (int i = 0; i < 3; ++i)
        {
          b_k[i] = b_k1[i];
        }
    }

  // Compute the eigenvalue
  eigenValue = 0.0;
  for (int i = 0; i < 3; ++i)
    {
      double temp = 0.0;
      for (int j = 0; j < 3; ++j)
        {
          temp += cov[i][j] * b_k1[j];
        }
      eigenValue += temp * b_k1[i];
    }

  // Copy the eigenvector
  for (int i = 0; i < 3; ++i)
    {
      eigenVector[i] = b_k1[i];
    }
}

// Function to project a 3D point onto the 2D plane defined by two eigenvectors
Point2D FT_Field::projectPointToPlane(const Point3D& point, const Point3D& centroid, const array<double, 3>& eigenVector1, const array<double, 3>& eigenVector2)
{
  double x = point.x - centroid.x;
  double y = point.y - centroid.y;
  double z = point.z - centroid.z;

  Point2D projectedPoint;
  projectedPoint.x = x * eigenVector1[0] + y * eigenVector1[1] + z * eigenVector1[2];
  projectedPoint.y = x * eigenVector2[0] + y * eigenVector2[1] + z * eigenVector2[2];
  return projectedPoint;
}


// Function to project a 3D point onto the 2D plane defined by two eigenvectors
int FT_Field::orientation_triangle(const Point3D& normale, const array<double, 3>& eigenVector1, const array<double, 3>& eigenVector2)
{
  Point3D eigenVector1_pt {eigenVector1[0],eigenVector1[1],eigenVector1[2]};
  Point3D eigenVector2_pt {eigenVector2[0],eigenVector2[1],eigenVector2[2]};
  double norm = scalarProduct(crossProduct(normale, eigenVector1_pt), eigenVector2_pt);
  if (norm >= 0.)
    return 1 ;
  else
    return -1 ;
}


vector<pair<double, array<double, 3>>> FT_Field::Main_2D_plane_eigenvectors(vector<Point3D> points)
{
  // Step 1: Compute the centroid
  Point3D centroid(computeCentroid(points));

  // Step 2: Compute the covariance matrix
  double cov[3][3];
  computeCovarianceMatrix(points, centroid, cov);

  // Step 3: Compute the eigenvalues and eigenvectors using power iteration
  double eigenVector1[3], eigenVector2[3], eigenVector3[3];
  double eigenValue1, eigenValue2, eigenValue3;

  powerIteration(cov, eigenVector1, eigenValue1);

  // Deflate the covariance matrix by removing the component of the first eigenvector
  for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
        {
          cov[i][j] -= eigenValue1 * eigenVector1[i] * eigenVector1[j];
        }
    }

  powerIteration(cov, eigenVector2, eigenValue2);

  // Deflate the covariance matrix by removing the component of the second eigenvector
  for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
        {
          cov[i][j] -= eigenValue2 * eigenVector2[i] * eigenVector2[j];
        }
    }

  powerIteration(cov, eigenVector3, eigenValue3);

  // Sort eigenvalues and corresponding eigenvectors in descending order
  vector<pair<double, array<double, 3>>> eigenPairs =
  {
    {eigenValue1, {eigenVector1[0], eigenVector1[1], eigenVector1[2]}},
    {eigenValue2, {eigenVector2[0], eigenVector2[1], eigenVector2[2]}},
    {eigenValue3, {eigenVector3[0], eigenVector3[1], eigenVector3[2]}}
  };
  sort(eigenPairs.rbegin(), eigenPairs.rend());

  return eigenPairs;
}

double FT_Field::triangleArea(const Point2D& p1, const Point2D& p2, const Point2D& p3)
{
  return 0.5 * std::abs((p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y)));
}

// Function to calculate the cross product of two vectors
Point3D FT_Field::crossProduct(const Point3D& u, const Point3D& v)
{
  Point3D cross;
  cross.x = u.y * v.z - u.z * v.y;
  cross.y = u.z * v.x - u.x * v.z;
  cross.z = u.x * v.y - u.y * v.x;
  return cross;
}
// Function to calculate the cross product of two vectors
double FT_Field::scalarProduct(const Point3D& u, const Point3D& v)
{
  return u.x*v.x + u.y*v.y + u.z*v.z;
}

// Function to calculate the magnitude of a vector
double FT_Field::magnitude(const Point3D& v)
{
  return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

// Function to calculate the area of a triangle given its vertices
double FT_Field::triangleArea3D(const Point3D& A, const Point3D& B, const Point3D& C)
{
  // Create vectors AB and AC
  Point3D AB = {B.x - A.x, B.y - A.y, B.z - A.z};
  Point3D AC = {C.x - A.x, C.y - A.y, C.z - A.z};

  // Calculate the cross product of AB and AC
  Point3D cross = crossProduct(AB, AC);

  // Calculate the magnitude of the cross product
  double area = magnitude(cross) / 2.0;

  return area;
}

void FT_Field::Calculate_Facette_Intersection_Area(DoubleTab& Surface_fa7init, DoubleTab& Surface_fa7fin, DoubleTab& Surface_intersection, vector<Point3D> points_fa7_originale, vector<Point3D> points_fa7_finale, IntTab points_triangle_originaux, IntTab points_triangle_finaux, IntTab& normale_triangle_originaux, IntTab& normale_triangle_finaux)
{
  // Les grandes lignes de l'algo :
  //  1 : on calcule les 2 vecteurs propres qui definissent le plan moyen de la structure formee par les triangles initiaux
  //  2 : on gere le changement eventuel d'orientation des facettes dans certains cas particuliers de repliement
  //  3 : on projette les points 3D (sommets des facettes) sur le plan 2D defini par les vecteurs propres
  //  4 : on calcule les intersections de tous les triangles initiaux avec tous les triangles finaux dans ce plan 2D

  if(points_fa7_originale.size()==0 and points_fa7_finale.size()==0)
    {
      return;
    }
  // Output the initial points
  int nbtriangle_originaux = points_triangle_originaux.dimension(0);
  int nbtriangle_finaux = points_triangle_finaux.dimension(0);
  // Step 1: Compute the centroid and The first two eigenvectors of the mean plane described by points_fa7_originale

  Point3D centroid = computeCentroid(points_fa7_originale);
  vector<pair<double, array<double, 3>>> eigenPairs = Main_2D_plane_eigenvectors(points_fa7_originale);

  array<double, 3> planeVector1 = eigenPairs[0].second;
  array<double, 3> planeVector2 = eigenPairs[1].second;

  // Project each 3D point onto the 2D plane
  vector<Point2D> projectedPoints_fa7_originale;
  for (const auto &point : points_fa7_originale)
    {
      projectedPoints_fa7_originale.push_back(projectPointToPlane(point, centroid, planeVector1, planeVector2));
    }
  vector<Point2D> projectedPoints_fa7_finale;
  for (const auto &point : points_fa7_finale)
    {
      projectedPoints_fa7_finale.push_back(projectPointToPlane(point, centroid, planeVector1, planeVector2));
    }

  // teste des mailles repliee sur elle-même
  // on calcul l'orientation + ou - de la normale du triangle

  for (int to = 0; to < nbtriangle_originaux; to++)
    {
      const Point3D normale = {normale_facette_initiale_(to, 0),normale_facette_initiale_(to, 1),normale_facette_initiale_(to, 2)};
      normale_triangle_originaux(to) = orientation_triangle(normale, planeVector1, planeVector2);
    }

  for (int to = 0; to < nbtriangle_finaux; to++)
    {
      // Pour les fa7 finale, il faut checker si la normale ne change pas d'orientation pendant le deplacement.
      // Pour ca, on calcule les deux normale a la fa7 avant et apres
      // si ça change, on change le signe de la normale de reference
      Point3D AB = {triangle_initiaux_(to, 1, 0)-triangle_initiaux_(to, 0, 0),
                    triangle_initiaux_(to, 1, 1)-triangle_initiaux_(to, 0, 1),
                    triangle_initiaux_(to, 1, 2)-triangle_initiaux_(to, 0, 2)
                   };
      Point3D AC = {triangle_initiaux_(to, 2, 0)-triangle_initiaux_(to, 0, 0),
                    triangle_initiaux_(to, 2, 1)-triangle_initiaux_(to, 0, 1),
                    triangle_initiaux_(to, 2, 2)-triangle_initiaux_(to, 0, 2)
                   };
      Point3D normale = crossProduct(AB, AC);
      int orientation_normale_avant_deplacement = orientation_triangle(normale, planeVector1, planeVector2);
      AB = {triangle_finaux_(to, 1, 0)-triangle_finaux_(to, 0, 0),
            triangle_finaux_(to, 1, 1)-triangle_finaux_(to, 0, 1),
            triangle_finaux_(to, 1, 2)-triangle_finaux_(to, 0, 2)
           };
      AC = {triangle_finaux_(to, 2, 0)-triangle_finaux_(to, 0, 0),
            triangle_finaux_(to, 2, 1)-triangle_finaux_(to, 0, 1),
            triangle_finaux_(to, 2, 2)-triangle_finaux_(to, 0, 2)
           };
      normale = crossProduct(AB, AC);
      int orientation_normale_apres_deplacement = orientation_triangle(normale, planeVector1, planeVector2);

      const Point3D normale_ref = {normale_facette_initiale_(to, 0),normale_facette_initiale_(to, 1),normale_facette_initiale_(to, 2)};
      normale_triangle_finaux(to) = orientation_triangle(normale_ref, planeVector1, planeVector2)*orientation_normale_apres_deplacement*orientation_normale_avant_deplacement;
    }

  // on projette les points 3D sur le plan 2D defini par les vecteurs propres
  // Puis on calcule les intersections de tous les triangles initiaux avec tous les triangles finaux
  for (int tf = 0; tf < nbtriangle_finaux; tf++)
    {
      Point2D p1f, p2f, p3f;
      p1f = projectedPoints_fa7_finale[points_triangle_finaux(tf,0)];
      p2f = projectedPoints_fa7_finale[points_triangle_finaux(tf,1)];
      p3f = projectedPoints_fa7_finale[points_triangle_finaux(tf,2)];

      Point2D t1[3] = {{p1f.x, p1f.y},{p2f.x, p2f.y},{p3f.x, p3f.y}};
      Surface_fa7fin[tf] = triangleArea(t1[0],t1[1],t1[2]);

      for (int to = 0; to < nbtriangle_originaux; to++)
        {
          Point2D p1o, p2o, p3o;
          p1o = projectedPoints_fa7_originale[points_triangle_originaux(to,0)];
          p2o = projectedPoints_fa7_originale[points_triangle_originaux(to,1)];
          p3o = projectedPoints_fa7_originale[points_triangle_originaux(to,2)];
          Point2D t2[3] = {{p1o.x, p1o.y},{p2o.x, p2o.y},{p3o.x, p3o.y}};
          Surface_fa7init[to] = triangleArea(t2[0],t2[1],t2[2]);
          Surface_intersection(tf,to) = intersectionArea(t1, t2);
        }
    }
  return;
}

void FT_Field::sauvegarder_triangle(const Maillage_FT_IJK& mesh, const int i, const int avant_apres_remaillage)
{
  bool triangle_already_sauv = false;
  if (avant_apres_remaillage==0)
    {
      /* Normalement pas necessaire de supprimer les duplicata avec la maniere dont est geré le parallelisme
       * On nest pas sense en avoir a ce stade
       * Il y en a malgre tout parfois dans les premiere iteration ou le remaillage est violent
       * Raison non identifiee --> paliatif
      */
      if (check_triangle_duplicata_)
        {
          for (int triangle = 0; triangle < triangle_initiaux_.dimension(0) ; triangle++)
            {
              Point3D p1new = {facettes_sommets_full_compo_(i, 0),facettes_sommets_full_compo_(i, 1),facettes_sommets_full_compo_(i, 2)};
              Point3D p2new = {facettes_sommets_full_compo_(i, 3),facettes_sommets_full_compo_(i, 4),facettes_sommets_full_compo_(i, 5)};
              Point3D p3new = {facettes_sommets_full_compo_(i, 6),facettes_sommets_full_compo_(i, 7),facettes_sommets_full_compo_(i, 8)};
              Point3D p1ref = {triangle_initiaux_(triangle, 0, 0),triangle_initiaux_(triangle, 0, 1),triangle_initiaux_(triangle, 0, 2)};
              Point3D p2ref = {triangle_initiaux_(triangle, 1, 0),triangle_initiaux_(triangle, 1, 1),triangle_initiaux_(triangle, 1, 2)};
              Point3D p3ref = {triangle_initiaux_(triangle, 2, 0),triangle_initiaux_(triangle, 2, 1),triangle_initiaux_(triangle, 2, 2)};
              if (p1new == p1ref and p2new == p2ref and p3new == p3ref)
                triangle_already_sauv = true;
            }
        }
      if(!triangle_already_sauv)
        {
          int index_triangle = triangle_initiaux_.dimension(0);
          triangle_initiaux_.resize(index_triangle+1, 3, 3);
          normale_facette_initiale_.resize(index_triangle+1, 3);
          Surfactant_facette_initiale_.resize(index_triangle+1);
          Surfactant_facette_initiale_(index_triangle) = facettes_sommets_full_compo_(i, 9);
          for (int i_som = 0; i_som < 3; i_som++)
            for (int dir = 0; dir < 3; dir++)
              {
                triangle_initiaux_(index_triangle, i_som, dir) = facettes_sommets_full_compo_(i, 3*i_som+dir);
              }
          for (int dir = 0; dir < 3; dir++)
            {
              normale_facette_initiale_(index_triangle, dir) = facettes_sommets_full_compo_(i, 12+dir);
            }
        }
    }
  if (avant_apres_remaillage==1)
    {
      /* Normalement pas necessaire de supprimer les duplicata avec la maniere dont est geré le parallelisme
       * On nest pas sense en avoir a ce stade
       * Il y en a malgre tout parfois dans les premiere iteration ou le remaillage est violent
       * Raison non identifiee --> paliatif
      */
      if (check_triangle_duplicata_)
        {
          for (int triangle = 0; triangle < triangle_finaux_.dimension(0) ; triangle++)
            {
              Point3D p1new = {facettes_sommets_full_compo_(i, 0),facettes_sommets_full_compo_(i, 1),facettes_sommets_full_compo_(i, 2)};
              Point3D p2new = {facettes_sommets_full_compo_(i, 3),facettes_sommets_full_compo_(i, 4),facettes_sommets_full_compo_(i, 5)};
              Point3D p3new = {facettes_sommets_full_compo_(i, 6),facettes_sommets_full_compo_(i, 7),facettes_sommets_full_compo_(i, 8)};
              Point3D p1ref = {triangle_finaux_(triangle, 0, 0),triangle_finaux_(triangle, 0, 1),triangle_finaux_(triangle, 0, 2)};
              Point3D p2ref = {triangle_finaux_(triangle, 1, 0),triangle_finaux_(triangle, 1, 1),triangle_finaux_(triangle, 1, 2)};
              Point3D p3ref = {triangle_finaux_(triangle, 2, 0),triangle_finaux_(triangle, 2, 1),triangle_finaux_(triangle, 2, 2)};
              if (p1new == p1ref and p2new == p2ref and p3new == p3ref)
                triangle_already_sauv = true;
            }
        }
      if(!triangle_already_sauv)
        {
          int index_triangle = triangle_finaux_.dimension(0);
          triangle_finaux_.resize(index_triangle+1, 3, 3);
          normale_facette_finale_.resize(index_triangle+1, 3);
          indice_facette_finaux_.resize(index_triangle+1);
          indice_facette_finaux_(index_triangle)=i;
          for (int i_som = 0; i_som < 3; i_som++)
            for (int dir = 0; dir < 3; dir++)
              triangle_finaux_(index_triangle, i_som, dir) = facettes_sommets_full_compo_(i, 3*i_som+dir);
          for (int dir = 0; dir < 3; dir++)
            {
              normale_facette_finale_(index_triangle, dir) = facettes_sommets_full_compo_(i, 12+dir);
            }
        }
    }
  return;
}


void FT_Field::remailler_FT_Field(Maillage_FT_IJK& mesh)
{

  int nb_triangle_fin = triangle_finaux_.dimension(0);
  int nb_triangle_init = triangle_initiaux_.dimension(0);

  if (nb_triangle_fin == 0 or nb_triangle_init==0)
    {
      return;
    }
  if ((nb_triangle_fin == 0 or nb_triangle_init==0) and nb_triangle_fin!=nb_triangle_init)
    {
      std::cout << "Erreur : nb_triangle_fin = " << nb_triangle_fin << " et triangle_initiaux_= " << nb_triangle_init << std::endl;
      Process::exit();
    }
  if (!variable_intensive_)
    {
      std::cout << "Erreur : la variable doit etre intensive ici" ;
      Process::exit();
    }

  ArrOfDouble surfactant_copy ;
  vector<Point3D> points_fa7_originale((nb_triangle_init)*3);
  vector<Point3D> points_fa7_finale((nb_triangle_fin)*3);

  //on remplit tous les points, il y a les doublons car sommets communs a plusieurs fa7
  for (int ifaforsom = 0; ifaforsom < nb_triangle_init; ifaforsom++)
    {
      int i_sommet ;
      for (i_sommet = 0; i_sommet < 3; i_sommet++)
        points_fa7_originale[3*ifaforsom+i_sommet]= {triangle_initiaux_(ifaforsom, i_sommet, 0),triangle_initiaux_(ifaforsom, i_sommet, 1),triangle_initiaux_(ifaforsom, i_sommet, 2)};
    }
  for (int ifaforsom= 0; ifaforsom < nb_triangle_fin; ifaforsom++)
    {
      int i_sommet ;
      for (i_sommet = 0; i_sommet < 3; i_sommet++)
        points_fa7_finale[3*ifaforsom+i_sommet]= {triangle_finaux_(ifaforsom, i_sommet, 0),triangle_finaux_(ifaforsom, i_sommet, 1),triangle_finaux_(ifaforsom, i_sommet, 2)};
    }

  // Supprimer les points dupliques
  vector<Point3D> points_fa7_finale_unique = removeDuplicates(points_fa7_finale);
  vector<Point3D> points_fa7_originale_unique = removeDuplicates(points_fa7_originale);

  // On remplit ensuite les triangles
  IntTab points_triangle_originaux(nb_triangle_init,3);
  IntTab points_triangle_finaux(nb_triangle_fin,3);

  for (int ifaforsom = 0; ifaforsom < nb_triangle_fin; ifaforsom++)
    {
      int i_sommet ;
      int index_point ;
      for (i_sommet = 0; i_sommet < 3; i_sommet++)
        {
          int n = static_cast<int>(points_fa7_finale_unique.size());
          for (index_point = 0; index_point < n; index_point++)
            {
              Point3D point_liste = points_fa7_finale_unique[index_point];
              Point3D point_triangle = {triangle_finaux_(ifaforsom, i_sommet, 0),triangle_finaux_(ifaforsom, i_sommet, 1),triangle_finaux_(ifaforsom, i_sommet, 2)};
              if(point_liste==point_triangle)
                {
                  points_triangle_finaux(ifaforsom, i_sommet)=index_point;
                }
            }
        }
    }

  for (int ifaforsom = 0; ifaforsom < nb_triangle_init; ifaforsom++)
    {
      int i_sommet ;
      int index_point ;

      for (i_sommet = 0; i_sommet < 3; i_sommet++)
        {
          int n = static_cast<int>(points_fa7_originale_unique.size());
          for (index_point = 0; index_point < n; index_point++)
            {
              Point3D point_liste = points_fa7_originale_unique[index_point];
              Point3D point_triangle = {triangle_initiaux_(ifaforsom, i_sommet, 0),triangle_initiaux_(ifaforsom, i_sommet, 1),triangle_initiaux_(ifaforsom, i_sommet, 2)};
              if(point_liste==point_triangle)
                {
                  points_triangle_originaux(ifaforsom, i_sommet)=index_point;
                }
            }
        }
    }

  // Check si la structure initiale de triangle est complete
  // Pour ca on verifie que chaque points est partage par au moins 2 triangles
  // Sinon cela veut dire quil y a des trous (cela peut arriver en // rarement). La source de lerreur est encore a chercher
  // Dans le cas dune structure incomplete, on interpole rien car ca peut donner nimporte quoi.
  // On impose la concentration moyenne de la structure partout.
  bool structure_complete = true;
  int n = static_cast<int>(points_fa7_originale_unique.size());
  for (int index_point = 0; index_point < n; index_point++)
    {
      int nb_triangle = 0;
      for (int ifaforsom = 0; ifaforsom < nb_triangle_init; ifaforsom++)
        {
          for (int i_sommet = 0; i_sommet < 3; i_sommet++)
            {
              if (points_triangle_originaux(ifaforsom, i_sommet)==index_point)
                {
                  nb_triangle+=1;
                }
            }
        }
      if (nb_triangle<2)
        {
          structure_complete=false;
        }
    }
  if(!structure_complete)
    {
      std::cout << "WARNING : structure incomplete de fa7 pour le remaillage surfactant" << std::endl;
      std::cout << "la concentration initiale moyenne est appliquee a toute la structure en question" << std::endl;
      std::cout << "Cela peut provoquer des erreurs de conservation" << std::endl;
    }

  // on calcule ensuite les intersections des triangles
  DoubleTab Surface_intersection(nb_triangle_fin, nb_triangle_init);
  DoubleTab Surface_fa7init(nb_triangle_init);
  DoubleTab Surface_fa7fin(nb_triangle_fin);
  IntTab normale_triangle_originaux(nb_triangle_init);
  IntTab normale_triangle_finaux(nb_triangle_fin);
  Calculate_Facette_Intersection_Area(Surface_fa7init, Surface_fa7fin, Surface_intersection, points_fa7_originale_unique, points_fa7_finale_unique, points_triangle_originaux, points_triangle_finaux,normale_triangle_originaux,normale_triangle_finaux);


  // on calcule les surfactants a ajouter au fa7 finale de maniere conservative

  surfactant_copy.resize(facettes_sommets_full_compo_.dimension(0));
  for (int fa = 0; fa < facettes_sommets_full_compo_.dimension(0); fa++)
    surfactant_copy(fa)=0.;


  // on calcule l'orientation principale (normale exterieur)

  int orientation_principale = 0;
  for (int tf = 0; tf < nb_triangle_fin; tf++)
    {
      orientation_principale+=normale_triangle_finaux(tf);
    }
  if (orientation_principale>=0)
    orientation_principale = 1;
  else
    orientation_principale = -1;


  // on calcule la quantite de surfactant sur la structure finale a partir des intersection calculee
  ArrOfDouble surfactant_final(nb_triangle_fin);


  if (structure_complete)
    {
      for (int tf = 0; tf < nb_triangle_fin; tf++)
        {
          int fa7 = indice_facette_finaux_(tf);

          double surfactant_facette_finale=0.;
          for (int to = 0; to < nb_triangle_init; to++)
            {
              if(normale_triangle_finaux(tf) == orientation_principale)
                {
                  surfactant_facette_finale+=Surfactant_facette_initiale_(to)*Surface_intersection(tf, to)/Surface_fa7init[to];
                }
            }
          surfactant_copy(fa7) = surfactant_facette_finale ;
          surfactant_final(tf) = surfactant_facette_finale ;
        }

      double reste_surfactant = 0.;
      double tot_surfactant = 0.;
      for (int to = 0; to < nb_triangle_init; to++)
        {
          tot_surfactant+=Surfactant_facette_initiale_(to);
        }
      reste_surfactant = tot_surfactant;
      for (int tf = 0; tf < nb_triangle_fin; tf++)
        {
          reste_surfactant-=surfactant_final(tf);
        }

      // on calcule les surfaces pour lesquelles ont a pas trouve dintersection (et qui explique le bilan non conservatif de surfactant)
      double surface_finale_manquante = 0. ;
      double surface_structure_finale = 0.;
      for (int tf = 0; tf < Surface_fa7fin.size_array(); tf++)
        {
          double conservation = 0. ;
          if(normale_triangle_finaux(tf) == orientation_principale)
            {
              surface_structure_finale+=Surface_fa7fin[tf];
            }
          for (int to = 0; to < Surface_fa7init.size_array(); to++)
            {
              if(normale_triangle_finaux(tf) == orientation_principale)
                conservation += Surface_intersection(tf,to);
            }
          conservation/=Surface_fa7fin[tf];
          if (abs(conservation-1.)>1.e-8)
            {
              if(normale_triangle_finaux(tf) == orientation_principale)
                surface_finale_manquante += (1.-conservation)*Surface_fa7fin[tf];
            }
        }


      for (int tf = 0; tf < Surface_fa7fin.size_array(); tf++)
        {
          double conservation = 0. ;
          for (int to = 0; to < Surface_fa7init.size_array(); to++)
            {
              if(normale_triangle_finaux(tf) == orientation_principale)
                conservation += Surface_intersection(tf,to);
            }
          conservation/=Surface_fa7fin[tf];

          // si reste_surfactant!= 0. et surface_finale_manquante != 0.
          // alors on est dans le cas d'une structure particuliere non parfaitement convexe
          // Les intersections peuvent donc ne pas etre complete
          // Normalement, la quantite de surfactant qui doit combler les intersections manquantes
          // vient de structures voisines (tres complique a mettre en place pour des evenements qui ne se produisent quasiment jamais)
          // Pour simplifier, on redistribue ce qui manque sur ces surfaces manquantes
          // Cela permet dassurer la conservation
          if (surface_finale_manquante!=0.)
            {
              if(normale_triangle_finaux(tf) == orientation_principale)
                {
                  surfactant_final(tf)+=reste_surfactant * ((1.-conservation) * Surface_fa7fin[tf]/surface_finale_manquante) ;
                  int fa7 = indice_facette_finaux_(tf);
                  surfactant_copy(fa7) = surfactant_final(tf);
                }
            }
          else
            {
              // Il peut encore rester des surfactants non distribues meme si conservation == 1
              // Dans le cas d'une maille retournee initiale rebarycentree
              // La structure initiale peut alors avoir une surface projetee plus grande que la structure finale (inverse du cas précedent de la structure convexe pour lequel Sfin > Sinit)
              // La totalite de la surface des mailles finale trouvent une intersection
              // Mais des bouts de surface initiale peuvent disparaitre
              // On redistribue ces bouts-là sur les mailles de la structure finale
              if(normale_triangle_finaux(tf) == orientation_principale)
                {
                  surfactant_final(tf)+=reste_surfactant * Surface_fa7fin[tf] / surface_structure_finale ;
                  int fa7 = indice_facette_finaux_(tf);
                  surfactant_copy(fa7) = surfactant_final(tf);
                }
            }
        }
    }
  else
    {
      // Paliatif pour structure incomplete (quil faudrait debug)
      // on impose la concentration moyenne de la structure initiale partout
      // Cela peut generer des erreurs de conservation
      double quantite_surf = 0.;
      double surf_tot = 0.;
      for (int to = 0; to < nb_triangle_init; to++)
        {
          if (nb_triangle_init>2)
            {
              if(normale_triangle_finaux(to) == orientation_principale)
                {
                  quantite_surf += Surfactant_facette_initiale_(to);
                  surf_tot += Surface_fa7init[to];
                }
            }
          else
            {
              quantite_surf += Surfactant_facette_initiale_(to);
              surf_tot += Surface_fa7init[to];
            }
        }

      for (int tf = 0; tf < nb_triangle_fin; tf++)
        {
          int fa7 = indice_facette_finaux_(tf);
          surfactant_copy(fa7) = (quantite_surf/surf_tot)*Surface_fa7fin[tf] ;
          surfactant_final(tf) = (quantite_surf/surf_tot)*Surface_fa7fin[tf] ;
        }
    }

  for (int tf = 0; tf < nb_triangle_fin; tf++)
    {
      int fa7 = indice_facette_finaux_(tf);
      facettes_sommets_full_compo_(fa7,9)=surfactant_copy(fa7);
    }

  // option d ecriture dans le cout pour debugage si necessaire
  if (print_debug_surfactant_)
    {
      double min_surfactant_init = 1.e+10;
      double max_surfactant_init = -1.e+10;
      for (int ifaforsom = 0; ifaforsom < nb_triangle_init; ifaforsom++)
        {
          if (normale_triangle_originaux(ifaforsom)==orientation_principale)
            {
              Point3D p1 = {triangle_initiaux_(ifaforsom, 0, 0), triangle_initiaux_(ifaforsom, 0, 1), triangle_initiaux_(ifaforsom, 0, 2)};
              Point3D p2 = {triangle_initiaux_(ifaforsom, 1, 0), triangle_initiaux_(ifaforsom, 1, 1), triangle_initiaux_(ifaforsom, 1, 2)};
              Point3D p3 = {triangle_initiaux_(ifaforsom, 2, 0), triangle_initiaux_(ifaforsom, 2, 1), triangle_initiaux_(ifaforsom, 2, 2)};
              double S = triangleArea3D( p1, p2, p3);
              if (Surfactant_facette_initiale_(ifaforsom)/S > max_surfactant_init)
                max_surfactant_init = Surfactant_facette_initiale_(ifaforsom)/S;
              if (Surfactant_facette_initiale_(ifaforsom)/S < min_surfactant_init)
                min_surfactant_init = Surfactant_facette_initiale_(ifaforsom)/S;
            }
        }

      bool bizarre = false;
      for (int ifaforsom = 0; ifaforsom < nb_triangle_fin; ifaforsom++)
        {
          if (normale_triangle_finaux(ifaforsom)==orientation_principale)
            {
              Point3D p1 = {triangle_finaux_(ifaforsom, 0, 0), triangle_finaux_(ifaforsom, 0, 1), triangle_finaux_(ifaforsom, 0, 2)};
              Point3D p2 = {triangle_finaux_(ifaforsom, 1, 0), triangle_finaux_(ifaforsom, 1, 1), triangle_finaux_(ifaforsom, 1, 2)};
              Point3D p3 = {triangle_finaux_(ifaforsom, 2, 0), triangle_finaux_(ifaforsom, 2, 1), triangle_finaux_(ifaforsom, 2, 2)};
              double S = triangleArea3D( p1, p2, p3);
              if(surfactant_final(ifaforsom)/S > max_surfactant_init*1.5)
                bizarre = true;
            }
        }

      if (bizarre)
        {
          std::cout << " valeur anormale de surfactant " <<std::endl;

          std::cout << "normale initiaux : " <<std::endl;
          for (int ifaforsom = 0; ifaforsom < nb_triangle_init; ifaforsom++)
            {
              std::cout << normale_triangle_originaux(ifaforsom) << std::endl;
            }
          std::cout <<std::endl;
          std::cout << "normale finale : " <<std::endl;
          for (int ifaforsom = 0; ifaforsom < nb_triangle_fin; ifaforsom++)
            {
              std::cout << normale_triangle_finaux(ifaforsom) << std::endl;
            }
          std::cout <<std::endl;

          std::cout << "triangle initiaux : "<<std::endl;
          for (int ifaforsom = 0; ifaforsom < nb_triangle_init; ifaforsom++)
            {
              //indice_facette_initiaux(i, ifaforsom);
              int i_sommet ;
              std::cout << "[" ;
              for (i_sommet = 0; i_sommet < 3; i_sommet++)
                std::cout << "[ " << triangle_initiaux_(ifaforsom, i_sommet, 0) << "," << triangle_initiaux_(ifaforsom, i_sommet, 1) << "," << triangle_initiaux_(ifaforsom, i_sommet, 2) << "],";
              std::cout << "],"<<std::endl;
            }
          std::cout << "triangle finaux : " <<std::endl;
          for (int ifaforsom = 0; ifaforsom < nb_triangle_fin; ifaforsom++)
            {
              //indice_facette_initiaux(i, ifaforsom);
              int i_sommet ;
              std::cout << "[" ;
              for (i_sommet = 0; i_sommet < 3; i_sommet++)
                std::cout << "[ " << triangle_finaux_(ifaforsom, i_sommet, 0) << "," << triangle_finaux_(ifaforsom, i_sommet, 1) << "," << triangle_finaux_(ifaforsom, i_sommet, 2)<< "],";
              std::cout << "],"<<std::endl;
            }
          std::cout << "Surfactant initiaux : " <<std::endl;
          for (int ifaforsom = 0; ifaforsom < nb_triangle_init; ifaforsom++)
            {
              std::cout << Surfactant_facette_initiale_(ifaforsom) << std::endl;
            }
          std::cout <<std::endl;

          std::cout << "Surfactant initiaux / Sfa7: " <<std::endl;
          for (int ifaforsom = 0; ifaforsom < nb_triangle_init; ifaforsom++)
            {
              Point3D p1 = {triangle_initiaux_(ifaforsom, 0, 0), triangle_initiaux_(ifaforsom, 0, 1), triangle_initiaux_(ifaforsom, 0, 2)};
              Point3D p2 = {triangle_initiaux_(ifaforsom, 1, 0), triangle_initiaux_(ifaforsom, 1, 1), triangle_initiaux_(ifaforsom, 1, 2)};
              Point3D p3 = {triangle_initiaux_(ifaforsom, 2, 0), triangle_initiaux_(ifaforsom, 2, 1), triangle_initiaux_(ifaforsom, 2, 2)};
              double S = triangleArea3D( p1, p2, p3);
              std::cout << Surfactant_facette_initiale_(ifaforsom)/S << std::endl;
            }
          std::cout <<std::endl;

          std::cout << "Surfactant finaux / Sfa7: " <<std::endl;
          for (int ifaforsom = 0; ifaforsom < nb_triangle_fin; ifaforsom++)
            {
              Point3D p1 = {triangle_finaux_(ifaforsom, 0, 0), triangle_finaux_(ifaforsom, 0, 1), triangle_finaux_(ifaforsom, 0, 2)};
              Point3D p2 = {triangle_finaux_(ifaforsom, 1, 0), triangle_finaux_(ifaforsom, 1, 1), triangle_finaux_(ifaforsom, 1, 2)};
              Point3D p3 = {triangle_finaux_(ifaforsom, 2, 0), triangle_finaux_(ifaforsom, 2, 1), triangle_finaux_(ifaforsom, 2, 2)};
              double S = triangleArea3D( p1, p2, p3);
              std::cout << surfactant_final(ifaforsom)/S << std::endl;
            }
          std::cout <<std::endl;

          std::cout << "indice facette finaux : " <<std::endl;
          for (int ifaforsom = 0; ifaforsom < nb_triangle_fin; ifaforsom++)
            {
              std::cout << indice_facette_finaux_(ifaforsom) << std::endl;
            }
          std::cout <<std::endl;

          for (const auto &point : Surface_fa7fin)
            {
              cout << "Surface fa7 fin = " << point << endl;
            }
          cout << endl;
          for (const auto &point : Surface_fa7init)
            {
              cout << "Surface fa7 init = " << point << endl;
            }
          cout << endl;

          for (int tf = 0; tf < Surface_fa7fin.size_array(); tf++)
            {
              for (int to = 0; to < Surface_fa7init.size_array(); to++)
                {
                  cout << "Intersection normalisee " << tf << " / " << to << " = " <<  Surface_intersection(tf,to)/Surface_fa7fin[tf] << std::endl;
                }
            }
          cout << endl;

          for (int tf = 0; tf < Surface_fa7fin.size_array(); tf++)
            {
              double conservation = 0. ;
              for (int to = 0; to < Surface_fa7init.size_array(); to++)
                {
                  if(normale_triangle_finaux(tf) == orientation_principale)
                    conservation += Surface_intersection(tf,to);
                }
              conservation/=Surface_fa7fin[tf];
              Process::Journal() << "conservation de la surface apres calcul des intersections remaillage pour le triangle final " << tf << " vaut " << conservation << std::endl;
              std::cout << "conservation de la surface apres calcul des intersections remaillage pour le triangle final " << tf << " vaut " << conservation << std::endl;


            }
        }
    }
}

ArrOfDouble FT_Field::check_conservation(const Maillage_FT_IJK& mesh)
{
  // calcule la quantite totale de surfactant pour chaque compo_connexe
  const int nbfa7=mesh.nb_facettes();
  const ArrOfInt& compo_connex = mesh.compo_connexe_facettes();
  const ArrOfDouble& Sfa7 = mesh.get_update_surface_facettes();

  int nb_bulles_reelles = 0 ;
  if (compo_connex.size_array()!=0)
    nb_bulles_reelles = max_array(compo_connex);
  nb_bulles_reelles = Process::mp_max(nb_bulles_reelles) + 1;
  ArrOfDouble Surfactant_par_compo(nb_bulles_reelles);
  for (int fa7=0 ; fa7<nbfa7 ; fa7++)
    {
      if(! mesh.facette_virtuelle(fa7) and compo_connex(fa7)>=0)
        {
          Surfactant_par_compo(compo_connex(fa7))+=FT_field_Array_(fa7)*Sfa7(fa7);
        }
    }
  for (int bulle=0 ; bulle<nb_bulles_reelles ; bulle++)
    {
      Surfactant_par_compo(bulle) = Process::mp_sum(Surfactant_par_compo(bulle));
    }
  return Surfactant_par_compo;
}

void FT_Field::correction_conservation_globale(const Maillage_FT_IJK& mesh, const ArrOfDouble& surfactant_avant_remaillage, const ArrOfDouble& surfactant_apres_remaillage)
{
  // redistribue l'erreur de conservation sur l'ensemble des facettes (par compo)
  // A eviter dutiliser. Cela signifie qu'une erreur de conservation locale a ete faite
  if(patch_conservation_surfactant_globale_)
    {
      const int nbfa7=mesh.nb_facettes();
      const ArrOfInt& compo_connex = mesh.compo_connexe_facettes();
      const ArrOfDouble& Sfa7 = mesh.get_surface_facettes();

      int nb_bulles_reelles = 0 ;
      if (compo_connex.size_array()!=0)
        nb_bulles_reelles = max_array(compo_connex);
      nb_bulles_reelles = Process::mp_max(nb_bulles_reelles) + 1;
      ArrOfDouble Surface_par_compo(nb_bulles_reelles);

      for (int fa7=0 ; fa7<nbfa7 ; fa7++)
        {
          if(! mesh.facette_virtuelle(fa7) and compo_connex(fa7)>=0)
            {
              Surface_par_compo(compo_connex(fa7))+=Sfa7(fa7);
            }
        }
      for (int bulle=0 ; bulle<nb_bulles_reelles ; bulle++)
        {
          Surface_par_compo(bulle) = Process::mp_sum(Surface_par_compo(bulle));
        }

      for (int fa7=0 ; fa7<nbfa7 ; fa7++)
        {
          if(! mesh.facette_virtuelle(fa7) and compo_connex(fa7)>=0)
            {
              FT_field_Array_(fa7) -= (surfactant_apres_remaillage(compo_connex(fa7))-surfactant_avant_remaillage(compo_connex(fa7)))/max(1.e-8,Surface_par_compo(compo_connex(fa7)));
            }
        }

      ArrOfDouble apres_correction = check_conservation(mesh);
      std::cout << "Apres correction" << std::endl;
      for (int bulle=0 ; bulle<nb_bulles_reelles ; bulle++)
        {
          double percent_error = 100.*(apres_correction(bulle)-surfactant_avant_remaillage(bulle))/max(1.e-8,apres_correction(bulle)) ;
          std::cout << " compo connex " << bulle << " = " <<  percent_error <<" %" << std::endl;
        }
    }
}

bool FT_Field::is_compo_in_proc(const int compo_connexe, const int pe_send)
{
// renvoie true sur la compo est dans le proc pe_send, false sinon
  int index_proc = -1 ;
  for (int proc = 0; proc < proc_numero_.size_array(); proc++)
    {
      if (proc_numero_(proc)==pe_send)
        {
          index_proc = proc ;
        }
    }
  if (index_proc==-1)
    {
      std::cout << "ERREUR : ne trouve pas avec quel proc echanger des compo_connexe" << std::endl;
      Process::exit();
    }

  for (int compo = 0; compo < compo_transmises_a_envoyer_.dimension(1); compo++)
    {
      if(compo_transmises_a_envoyer_(index_proc, compo)-1 == compo_connexe)
        return true;
    }
  return false;
}


void FT_Field::champ_sommet_from_facettes(const ArrOfInt& compo_connexe_facettes, const Maillage_FT_IJK& mesh)
{
  compo_connexe_sommets_.resize(mesh.nb_sommets());
  const int nbfa = compo_connexe_facettes.size_array();
  const IntTab& facettes=mesh.facettes();
  for (int fa7=0 ; fa7<nbfa ; fa7++)
    {
      if(! mesh.facette_virtuelle(fa7))
        {
          for (int sommet_fa7=0 ; sommet_fa7<3 ; sommet_fa7++)
            {
              int indice_sommet = facettes(fa7,sommet_fa7);
              compo_connexe_sommets_[indice_sommet] = compo_connexe_facettes(fa7);
            }
        }
    }
  const Desc_Structure_FT& desc_sommets = mesh.desc_sommets();
  desc_sommets.collecter_espace_virtuel(compo_connexe_sommets_, MD_Vector_tools::EV_MAX);
}


bool FT_Field::sauv_num_pe_echange(int pe)
{
  for (int proc_deja_sauv = 0; proc_deja_sauv < proc_deja_echange_.size_array(); proc_deja_sauv++)
    {
      if (proc_deja_echange_(proc_deja_sauv)==pe)
        {
          return true;
        }
    }
  nb_proc_echange_++;
  proc_deja_echange_.resize(nb_proc_echange_);
  proc_deja_echange_(nb_proc_echange_-1) = pe;
  return false;
}

void FT_Field::update_FT_Field_local_from_full_compo(const Maillage_FT_IJK& mesh)
{
  for (int fa = 0; fa < mesh.nb_facettes(); fa++)
    FT_field_Array_(fa) = -123.;

  int index_init = index_local_Ft_field_;
  for (int fa = 0; fa < mesh.nb_facettes(); fa++)
    {
      index_init++;
      const int indice_fa_locale = int(facettes_sommets_full_compo_(index_init, 10));
      FT_field_Array_(indice_fa_locale) = facettes_sommets_full_compo_(index_init, 9);
    }
}


void FT_Field::completer_compo_connexe_partielle(const Maillage_FT_IJK& mesh, const Domaine_IJK& splitting, const DoubleTab& liste_sommets_apres_deplacement, const DoubleTab& liste_sommets_avant_deplacement, const ArrOfInt& compo_connexe_sommets_deplace)
{
// Version de la methode qui prend en argument :
// la liste des coordonnees de sommets avant deplacement
// la liste des coordonnees de sommets apres deplacement
// la liste contenant la compo_connexe associee
// permet de mutualiser la methode pour le remaillage (supprimer petites arretes) et pour le barycentrage.
// Cette methode complete les tableaux  facettes_sommets_full_compo_ et liste_sommets_et_deplacements_.
// Les tableaux contiennent alors la totalite des informations de toutes les composantes connexes qui traversent le proc (position des sommets, deplacement, facette).
// Mais il n'y a plus d'information sur l'indexation locale

  if (!((liste_sommets_apres_deplacement.dimension(0) == liste_sommets_avant_deplacement.dimension(0)) and (liste_sommets_avant_deplacement.dimension(0) == compo_connexe_sommets_deplace.size_array())))
    {
      std::cout << "Erreur dans completer_compo_connexe_partielle : les tableaux en arguments ne font pas tous la meme taille" << std::endl;
      Process::exit();
    }
  const ArrOfInt& compo_connex = mesh.compo_connexe_facettes();
  champ_sommet_from_facettes(compo_connex, mesh);

  IntVect slice_3D(3) ;
  IntVect N_3D(3) ;
  for (int dir = 0; dir < 3; dir++)
    {
      slice_3D(dir)= splitting.get_local_slice_index(dir);
      N_3D(dir)=splitting.get_nprocessor_per_direction(dir);
    }

  nb_proc_echange_=0;
  nb_compo_a_envoyer_max_=0;
  proc_deja_echange_.resize(0);
  proc_numero_.resize(0);
  compo_transmises_a_envoyer_.resize(0, 0);
  // on boucle sur l'ensemble des proc_voisins dans chaque direction (voisin face/arrete/coin)
  // compatible tri-periodique
  // on echange les compo_connexe concernees
  for (int iproc = -1; iproc < 2; iproc++)
    {
      for (int jproc = -1; jproc < 2; jproc++)
        {
          for (int kproc = -1; kproc < 2; kproc++)
            {
              int pe_send = splitting.get_processor_by_ijk(((slice_3D(0)-iproc)%N_3D(0)+N_3D(0))%N_3D(0), ((slice_3D(1)-jproc)%N_3D(1)+N_3D(1))%N_3D(1), ((slice_3D(2)-kproc)%N_3D(2)+N_3D(2))%N_3D(2));
              int pe_rcdv = splitting.get_processor_by_ijk(((slice_3D(0)+iproc)%N_3D(0)+N_3D(0))%N_3D(0), ((slice_3D(1)+jproc)%N_3D(1)+N_3D(1))%N_3D(1), ((slice_3D(2)+kproc)%N_3D(2)+N_3D(2))%N_3D(2));
              exchange_compo_connexe(pe_send, pe_rcdv, mesh);
            }
        }
    }

  // On connait les compo connex a echanger
  // On peut maintenant echanger les donnees
  // on remet le tableau des echange realises a 0
  nb_proc_echange_=0;
  proc_deja_echange_.resize(0);
  facettes_sommets_full_compo_.resize(0, 0);
  liste_sommets_et_deplacements_.resize(0,0);
  int moi = splitting.get_processor_by_ijk(slice_3D(0), slice_3D(1), slice_3D(2));
  exchange_data(moi, moi, mesh, liste_sommets_avant_deplacement, liste_sommets_apres_deplacement, compo_connexe_sommets_deplace);

  for (int iproc = -1; iproc < 2; iproc++)
    {
      for (int jproc = -1; jproc < 2; jproc++)
        {
          for (int kproc = -1; kproc < 2; kproc++)
            {
              int pe_send = splitting.get_processor_by_ijk(((slice_3D(0)-iproc)%N_3D(0)+N_3D(0))%N_3D(0), ((slice_3D(1)-jproc)%N_3D(1)+N_3D(1))%N_3D(1), ((slice_3D(2)-kproc)%N_3D(2)+N_3D(2))%N_3D(2));
              int pe_rcdv = splitting.get_processor_by_ijk(((slice_3D(0)+iproc)%N_3D(0)+N_3D(0))%N_3D(0), ((slice_3D(1)+jproc)%N_3D(1)+N_3D(1))%N_3D(1), ((slice_3D(2)+kproc)%N_3D(2)+N_3D(2))%N_3D(2));
              if (!(pe_send==moi and pe_rcdv==moi))
                exchange_data(pe_send, pe_rcdv, mesh, liste_sommets_avant_deplacement, liste_sommets_apres_deplacement, compo_connexe_sommets_deplace);
            }
        }
    }

  // on supprime les eventuels doublons pour que chaque proc ait bien la meme chose
  // necessaire car on a echange certaines valeurs virtuelles pour etre sur de perdre personne
  // on demarre a partir de nb_facettes() pour etre sûr que les donnes du proc courant ne seront pas supprimees
  /* double tol_sauv = get_tolerance_point_identique();
  set_tolerance_point_identique(tol_sauv/100.);

  std::cout << "avant suppression doublons = " << facettes_sommets_full_compo_.dimension(0) << std::endl;
  int size_init = facettes_sommets_full_compo_.dimension(0);
  {
    for (int index = mesh.nb_facettes() ; index < size_init; index++)
      {
        for (int index2 = index + 1 ; index2 < size_init; )
          {
            Point3D p1ref= {facettes_sommets_full_compo_(index, 0),facettes_sommets_full_compo_(index, 1),facettes_sommets_full_compo_(index, 2)};
            Point3D p2ref= {facettes_sommets_full_compo_(index, 3),facettes_sommets_full_compo_(index, 4),facettes_sommets_full_compo_(index, 5)};
            Point3D p3ref= {facettes_sommets_full_compo_(index, 6),facettes_sommets_full_compo_(index, 7),facettes_sommets_full_compo_(index, 8)};

            Point3D p1ref2= {facettes_sommets_full_compo_(index2, 0),facettes_sommets_full_compo_(index2, 1),facettes_sommets_full_compo_(index2, 2)};
            Point3D p2ref2= {facettes_sommets_full_compo_(index2, 3),facettes_sommets_full_compo_(index2, 4),facettes_sommets_full_compo_(index2, 5)};
            Point3D p3ref2= {facettes_sommets_full_compo_(index2, 6),facettes_sommets_full_compo_(index2, 7),facettes_sommets_full_compo_(index2, 8)};
            double faref = facettes_sommets_full_compo_(index, 9);
            double faref2 = facettes_sommets_full_compo_(index2, 9);
            if(p1ref==p1ref2 and p2ref==p2ref2 and p3ref==p3ref2 and abs(faref-faref2)<1.e-10)
              {
                std::cout << "triangle sup = " << " " << facettes_sommets_full_compo_(index, 0) << " " << facettes_sommets_full_compo_(index, 1) << " " << facettes_sommets_full_compo_(index, 2) << " " << facettes_sommets_full_compo_(index, 3) << " " << facettes_sommets_full_compo_(index, 4) << facettes_sommets_full_compo_(index, 5) << " " << facettes_sommets_full_compo_(index, 6)<< " " <<facettes_sommets_full_compo_(index, 7)<< " " <<facettes_sommets_full_compo_(index, 8) << " " << facettes_sommets_full_compo_(index, 9) << std::endl ;
                for (int k = index2; k < size_init - 1; ++k)
                  for (int item = 0 ; item < facettes_sommets_full_compo_.dimension(1); item++)
                    facettes_sommets_full_compo_(k, item)=facettes_sommets_full_compo_(k+1, item);

                --size_init;
              }
            else
              {
                ++index2;
              }
          }
      }
  }
  facettes_sommets_full_compo_.resize(size_init, facettes_sommets_full_compo_.dimension(1));
  std::cout << "apres suppression doublons = " << facettes_sommets_full_compo_.dimension(0) << std::endl;

  size_init = liste_sommets_et_deplacements_.dimension(0);
  for (int index = 0 ; index < size_init; index++)
    {
      for (int index2 = index + 1 ; index2 < size_init; )
        {
          Point3D p1ref= {liste_sommets_et_deplacements_(index, 0),liste_sommets_et_deplacements_(index, 1),liste_sommets_et_deplacements_(index, 2)};
          Point3D p1ref2= {liste_sommets_et_deplacements_(index2, 0),liste_sommets_et_deplacements_(index2, 1),liste_sommets_et_deplacements_(index2, 2)};
          if(p1ref==p1ref2)
            {
              for (int k = index2; k < size_init - 1; ++k)
                for (int item = 0 ; item < liste_sommets_et_deplacements_.dimension(1); item++)
                  liste_sommets_et_deplacements_(k, item)=liste_sommets_et_deplacements_(k+1, item);

              --size_init;
            }
          else
            {
              ++index2;
            }
        }
    }
  liste_sommets_et_deplacements_.resize(size_init, liste_sommets_et_deplacements_.dimension(1));
  set_tolerance_point_identique(tol_sauv);

  */
  // si on veut etre sur que ca marche bien en //, il faut que chaque proc fasse les operations dans le meme ordre !
  // on peut trier liste_sommets_et_deplacements_full_compo selon par exemple, les normes des points deplaces

  int n = liste_sommets_et_deplacements_.dimension(0);
  std::vector<double> arr(n);
  sorted_index_.resize(n);
  for (int item = 0; item < n; item++)
    {
      arr[item]=norme({liste_sommets_et_deplacements_(item,0),liste_sommets_et_deplacements_(item,1),liste_sommets_et_deplacements_(item,2)});
    }
  // Vector to store the sorted indices
  std::vector<long unsigned int> indices(n);
  // Sort the array and track indices
  sortAndTrackIndices(arr, indices);
  for (int item = 0; item < n; item++)
    {
      sorted_index_(item)=int(indices[item]);
    }

  if (print_debug_surfactant_)
    {
      Process::barrier();
      int nb_proc = splitting.get_nprocessor_per_direction(0)*splitting.get_nprocessor_per_direction(1)*splitting.get_nprocessor_per_direction(2);
      int ny = splitting.get_nprocessor_per_direction(1);
      int nz = splitting.get_nprocessor_per_direction(2);
      int x = splitting.get_local_slice_index(0);
      int y = splitting.get_local_slice_index(1);
      int z = splitting.get_local_slice_index(2);
      int proc_index = z + y * nz + x * ny * nz ;
      int index_global = 0;
      for (int proc = 0; proc < nb_proc; proc++)
        {
          index_global = Process::mp_max(index_global);
          Process::barrier();
          if (proc == proc_index)
            {
              std::cout << "proc = " << proc_index << " total facette connues = " << facettes_sommets_full_compo_.dimension(0) << std::endl;
              std::cout << "proc = " << proc_index << " nb facette_local = " << mesh.nb_facettes() <<  std::endl;
              std::cout << "proc = " << proc_index << " total sommets connus = " << liste_sommets_et_deplacements_.dimension(0) <<  std::endl;
              std::cout << "proc = " << proc_index << " nb liste_sommets local = " << liste_sommets_apres_deplacement.dimension(0) <<  std::endl;
              std::cout << "proc = " << proc_index << " nb sommets local = " << mesh.nb_sommets() <<  std::endl;
            }
        }
      Process::barrier();
    }
}

Point3D FT_Field::calculer_normale_apres_deplacement(const int fa, const int somfa7, const Vecteur3 pos_apres_dep)
{
// si n = AB x AC
// alors AC' = n x AB doit etre dans la meme direction que AC original, soit AC'.AC > 0
// Soit AB le vecteur qui joint les deux premiers sommets
  Point3D AB = {facettes_sommets_full_compo_(fa, 3*1+0)-facettes_sommets_full_compo_(fa, 0),
                facettes_sommets_full_compo_(fa, 3*1+1)-facettes_sommets_full_compo_(fa, 1),
                facettes_sommets_full_compo_(fa, 3*1+2)-facettes_sommets_full_compo_(fa, 2)
               };

  Point3D AC = {facettes_sommets_full_compo_(fa, 3*2+0)-facettes_sommets_full_compo_(fa, 0),
                facettes_sommets_full_compo_(fa, 3*2+1)-facettes_sommets_full_compo_(fa, 1),
                facettes_sommets_full_compo_(fa, 3*2+2)-facettes_sommets_full_compo_(fa, 2)
               };
  Point3D nrecon = crossProduct(AB, AC);
  Point3D nreel = {facettes_sommets_full_compo_(fa, 12),facettes_sommets_full_compo_(fa, 13),facettes_sommets_full_compo_(fa, 14)};

  int A,B,C ;
  if (scalarProduct(nrecon,nreel)>=0.)
    {
      // alors on a la bonne definition du produit vectoriel qui donne n_reel
      A = 0;
      B = 1;
      C = 2;
    }
  else
    {
      // alors on a la definition inverse du produit vectoriel qui donne n_reel
      // il suffit d'inverse un vecteur pour obtenir la bonne
      A = 0;
      B = 2;
      C = 1;
    }

// une fois obtenu, on peut reconstruire la nouvelle normale apres deplacement

  if (somfa7==A)// sommet A qui bouge
    {
      AB = {facettes_sommets_full_compo_(fa, 3*B+0)-pos_apres_dep[0],
            facettes_sommets_full_compo_(fa, 3*B+1)-pos_apres_dep[1],
            facettes_sommets_full_compo_(fa, 3*B+2)-pos_apres_dep[2]
           };
      AC = {facettes_sommets_full_compo_(fa, 3*C+0)-pos_apres_dep[0],
            facettes_sommets_full_compo_(fa, 3*C+1)-pos_apres_dep[1],
            facettes_sommets_full_compo_(fa, 3*C+2)-pos_apres_dep[2]
           };
    }
  else if (somfa7==B)
    {
      AB = {pos_apres_dep[0]-facettes_sommets_full_compo_(fa, 3*A+0),
            pos_apres_dep[1]-facettes_sommets_full_compo_(fa, 3*A+1),
            pos_apres_dep[2]-facettes_sommets_full_compo_(fa, 3*A+2)
           };
      AC = {facettes_sommets_full_compo_(fa, 3*C+0)-facettes_sommets_full_compo_(fa, 3*A+0),
            facettes_sommets_full_compo_(fa, 3*C+1)-facettes_sommets_full_compo_(fa, 3*A+1),
            facettes_sommets_full_compo_(fa, 3*C+2)-facettes_sommets_full_compo_(fa, 3*A+2)
           };
    }
  else if (somfa7==C)
    {
      AB = {facettes_sommets_full_compo_(fa, 3*B+0)-facettes_sommets_full_compo_(fa, 3*A+0),
            facettes_sommets_full_compo_(fa, 3*B+1)-facettes_sommets_full_compo_(fa, 3*A+1),
            facettes_sommets_full_compo_(fa, 3*B+2)-facettes_sommets_full_compo_(fa, 3*A+2)
           };
      AC = {pos_apres_dep[0]-facettes_sommets_full_compo_(fa, 3*A+0),
            pos_apres_dep[1]-facettes_sommets_full_compo_(fa, 3*A+1),
            pos_apres_dep[2]-facettes_sommets_full_compo_(fa, 3*A+2)
           };
    }
  return crossProduct(AB, AC);
}

void FT_Field::exchange_compo_connexe(int pe_send_, /* processor to send data */
                                      int pe_recv_ /* processor to received data */,
                                      const Maillage_FT_IJK& mesh)
{
  const ArrOfInt& compo_connexe_facettes=mesh.compo_connexe_facettes();
  if (pe_send_ == Process::me() && pe_recv_ == Process::me())
    {
      // Self (periodicity on same processor)
      return;
    }

  bool echange_already_made = sauv_num_pe_echange(pe_send_);
  if (echange_already_made)
    {
      // ATTENTION : Dans le cas de calculs a faible nb de proc, les memes echanges risquent de se realiser plusieurs fois par periodicite des procs
      // Il faut s'assurer qu'on ne le fait qu'une seule fois
      return;
    }

  // Le volume de donnees echangees ne sont pas connues a l'avance
  // On est oblige de faire l'echange en plusieurs etapes
  // 1 - Le proc partage le nb de compo_connexe qui le traverse
  // 2 - Le proc partage les valeurs des compo_connexe qui le traverse

  const int int_size = sizeof(int);
  int *send_buffer_nb_compo_connexe = 0;
  int *recv_buffer_nb_compo_connexe = 0;
  int *send_buffer_compo_connexe = 0;
  int *recv_buffer_compo_connexe = 0;
  IntTab compo_to_send(1);
  compo_to_send(0)=-1000;
  int nb_compo_a_envoyer = 0;
  int nb_compo_recu = 0;

  ////////////////////////////////// Stockage des compo connexe a envoyer
  if (pe_send_ >= 0)
    {
      for (int fa = 0; fa < mesh.nb_facettes(); fa++)
        {
          bool new_compo = true ;
          if (!mesh.facette_virtuelle(fa))
            {
              for (int compo_deja_sauv = 0; compo_deja_sauv < compo_to_send.size_array(); compo_deja_sauv++)
                {
                  if(compo_connexe_facettes[fa]==compo_to_send(compo_deja_sauv))
                    new_compo = false;
                  break;
                }
              if(new_compo)
                {
                  nb_compo_a_envoyer++;
                  compo_to_send.resize(nb_compo_a_envoyer);
                  compo_to_send(nb_compo_a_envoyer-1) = compo_connexe_facettes[fa];
                }
            }
        }
    }

  if (pe_send_ >= 0)
    {
      ////////////////////// envoi du nb de compo connexe /////////////////////
      send_buffer_nb_compo_connexe = new int[1];
      int *buf_nb_compo_connexe = send_buffer_nb_compo_connexe;
      *buf_nb_compo_connexe = nb_compo_a_envoyer ;
      ////////////////////// envoi des compo connexe elle-meme  /////////////////////
      send_buffer_compo_connexe = new int[nb_compo_a_envoyer];
      int *buf_compo_connexe = send_buffer_compo_connexe;
      for (int compo = 0; compo < nb_compo_a_envoyer; compo++, buf_compo_connexe++)
        *buf_compo_connexe = compo_to_send(compo);

    }

  ////////////////////// reception du nb de compo connexe /////////////////////
  if (pe_recv_ >= 0)
    {
      recv_buffer_nb_compo_connexe = new int[1];
    }
  ::envoyer_recevoir(send_buffer_nb_compo_connexe, int_size, pe_send_, recv_buffer_nb_compo_connexe, int_size, pe_recv_);

  if (pe_recv_ >= 0)
    {
      int *buf_nb_compo_connexe = recv_buffer_nb_compo_connexe;
      nb_compo_recu = *buf_nb_compo_connexe;
      if (nb_compo_recu>nb_compo_a_envoyer_max_)
        nb_compo_a_envoyer_max_ = nb_compo_recu;
    }

  ////////////////////// reception des compo connexe elle-memes /////////////////////
  if (pe_recv_ >= 0)
    {
      recv_buffer_compo_connexe = new int[nb_compo_recu];
    }
  ::envoyer_recevoir(send_buffer_compo_connexe, nb_compo_a_envoyer*int_size, pe_send_, recv_buffer_compo_connexe, nb_compo_recu*int_size, pe_recv_);

  if (pe_recv_ >= 0)
    {
      // on stocke les compo+1 pour que les compo commence a 1 et pas 0
      // Les valeurs initalisees non changees reste donc a 0
      proc_numero_.resize(nb_proc_echange_);
      proc_numero_(nb_proc_echange_-1)=pe_recv_;
      compo_transmises_a_envoyer_.resize(nb_proc_echange_, nb_compo_a_envoyer_max_);
      int *buf_compo_connexe = recv_buffer_compo_connexe;
      for (int compo = 0; compo < nb_compo_recu; compo++, buf_compo_connexe++)
        compo_transmises_a_envoyer_(nb_proc_echange_-1, compo) = *buf_compo_connexe+1;
    }

  delete[] send_buffer_nb_compo_connexe;
  delete[] recv_buffer_nb_compo_connexe;
  delete[] send_buffer_compo_connexe;
  delete[] recv_buffer_compo_connexe;
}

void FT_Field::exchange_data(int pe_send_, /* processor to send data */
                             int pe_recv_ /* processor to received data */,
                             const Maillage_FT_IJK& mesh, const DoubleTab& liste_sommets_avant_deplacement, const DoubleTab& liste_sommets_apres_deplacement, const ArrOfInt& compo_connexe_sommets_deplace)
{

  const DoubleTab& sommets=mesh.sommets();
  const IntTab& facettes=mesh.facettes();
  const ArrOfInt& compo_connexe_facettes=mesh.compo_connexe_facettes();
  const DoubleTab& normale_facette=mesh.get_update_normale_facettes();
  bool echange_already_made = sauv_num_pe_echange(pe_send_);
  if (echange_already_made)
    {
      // ATTENTION : Dans le cas de calculs a faible nb de proc, les memes echanges risquent de se realiser plusieurs fois par periodicite des procs
      // Il faut s'assurer qu'on ne le fait qu'une seule fois
      return;
    }

  int dimensions_tab = 15;
  if (pe_send_ == Process::me() && pe_recv_ == Process::me())
    {
      // echange avec soi-meme, je remplis le tableau avec les valeurs locales deja connues
      int size_avant_echange = facettes_sommets_full_compo_.dimension(0);
      int new_size = size_avant_echange + mesh.nb_facettes() ;
      facettes_sommets_full_compo_.resize(new_size, dimensions_tab);
      int index = size_avant_echange-1;
      index_local_Ft_field_=index;

      for (int fa = 0; fa < mesh.nb_facettes(); fa++)
        {
          //if (!mesh.facette_virtuelle(fa))
          {
            index++;
            facettes_sommets_full_compo_(index, 0) = sommets(facettes(fa, 0), 0) ;
            facettes_sommets_full_compo_(index, 1) = sommets(facettes(fa, 0), 1) ;
            facettes_sommets_full_compo_(index, 2) = sommets(facettes(fa, 0), 2) ;
            facettes_sommets_full_compo_(index, 3) = sommets(facettes(fa, 1), 0) ;
            facettes_sommets_full_compo_(index, 4) = sommets(facettes(fa, 1), 1) ;
            facettes_sommets_full_compo_(index, 5) = sommets(facettes(fa, 1), 2) ;
            facettes_sommets_full_compo_(index, 6) = sommets(facettes(fa, 2), 0) ;
            facettes_sommets_full_compo_(index, 7) = sommets(facettes(fa, 2), 1) ;
            facettes_sommets_full_compo_(index, 8) = sommets(facettes(fa, 2), 2) ;
            facettes_sommets_full_compo_(index, 9) = FT_field_Array_(fa);
            facettes_sommets_full_compo_(index, 10) = double(fa);
            if (!mesh.facette_virtuelle(fa))
              {
                facettes_sommets_full_compo_(index, 11) = 0.;
              }
            else
              {
                facettes_sommets_full_compo_(index, 11) = 1.;
              }
            facettes_sommets_full_compo_(index, 12) = normale_facette(fa, 0);
            facettes_sommets_full_compo_(index, 13) = normale_facette(fa, 1);
            facettes_sommets_full_compo_(index, 14) = normale_facette(fa, 2);
          }
        }


      size_avant_echange = liste_sommets_et_deplacements_.dimension(0);
      new_size = size_avant_echange + liste_sommets_avant_deplacement.dimension(0) ;
      liste_sommets_et_deplacements_.resize(new_size, 6);
      index = size_avant_echange-1;
      for (int som = 0; som < liste_sommets_avant_deplacement.dimension(0); som++)
        {
          index++;
          liste_sommets_et_deplacements_(index, 0) = liste_sommets_avant_deplacement(som, 0);
          liste_sommets_et_deplacements_(index, 1) = liste_sommets_avant_deplacement(som, 1);
          liste_sommets_et_deplacements_(index, 2) = liste_sommets_avant_deplacement(som, 2);
          liste_sommets_et_deplacements_(index, 3) = liste_sommets_apres_deplacement(som, 0);
          liste_sommets_et_deplacements_(index, 4) = liste_sommets_apres_deplacement(som, 1);
          liste_sommets_et_deplacements_(index, 5) = liste_sommets_apres_deplacement(som, 2);
        }

      return;
    }

  // Le volume de donnees echangees ne sont pas connues a l'avance
  // On est oblige de faire l'echange en plusieurs etapes
  // 1 - Une fois ces compo_connexe_connu, on repere les donnes a echanger
  // 2 - on transmet au proc le volume de donnees qu'il va recevoir
  // 3 - on transmet au proc les donnees elles-même

  const int double_size = sizeof(double);
  const int int_size = sizeof(int);
  int *send_buffer_volume_facette = 0;
  int *recv_buffer_volume_facette = 0;
  int *send_buffer_volume_sommet = 0;
  int *recv_buffer_volume_sommet = 0;
  double *send_buffer_data_facette = 0;
  double *recv_buffer_data_facette = 0;
  double *send_buffer_data_sommet = 0;
  double *recv_buffer_data_sommet = 0;
  int nb_fa7_a_envoyer = 0;
  int nb_som_a_envoyer = 0;


  DoubleTab item_to_send_facette;
  DoubleTab item_to_send_sommet;
  int volume_donnees_facette_envoyees = 0;
  int volume_donnees_sommet_envoyees = 0;
  int volume_donnees_facette_recue = 0;
  int volume_donnees_sommet_recue = 0;

  if (pe_send_ >= 0)
    {
      // A partir des compo_connexe recue, on cherche maintenant les donnes a envoyer
      // Cest a dire toute les donnees liees aux facettes dont la compo_connexe est partagee avec le proc a qui on partage
      // position des sommets qui compose la fa7 et quantite de surfactant associe a la fa7
      for (int fa = 0; fa < mesh.nb_facettes(); fa++)
        {
          if (!mesh.facette_virtuelle(fa) and is_compo_in_proc(compo_connexe_facettes[fa], pe_send_))
            {
              // alors il faut envoyer l'info de cette fa7
              nb_fa7_a_envoyer++;
              item_to_send_facette.resize(nb_fa7_a_envoyer, dimensions_tab);
              item_to_send_facette(nb_fa7_a_envoyer-1, 0) = sommets(facettes(fa, 0), 0) ;
              item_to_send_facette(nb_fa7_a_envoyer-1, 1) = sommets(facettes(fa, 0), 1) ;
              item_to_send_facette(nb_fa7_a_envoyer-1, 2) = sommets(facettes(fa, 0), 2) ;
              item_to_send_facette(nb_fa7_a_envoyer-1, 3) = sommets(facettes(fa, 1), 0) ;
              item_to_send_facette(nb_fa7_a_envoyer-1, 4) = sommets(facettes(fa, 1), 1) ;
              item_to_send_facette(nb_fa7_a_envoyer-1, 5) = sommets(facettes(fa, 1), 2) ;
              item_to_send_facette(nb_fa7_a_envoyer-1, 6) = sommets(facettes(fa, 2), 0) ;
              item_to_send_facette(nb_fa7_a_envoyer-1, 7) = sommets(facettes(fa, 2), 1) ;
              item_to_send_facette(nb_fa7_a_envoyer-1, 8) = sommets(facettes(fa, 2), 2) ;
              item_to_send_facette(nb_fa7_a_envoyer-1, 9) = FT_field_Array_(fa);
              item_to_send_facette(nb_fa7_a_envoyer-1, 10) = -1.;
              item_to_send_facette(nb_fa7_a_envoyer-1, 11) = 0.;
              item_to_send_facette(nb_fa7_a_envoyer-1, 12) = normale_facette(fa, 0);
              item_to_send_facette(nb_fa7_a_envoyer-1, 13) = normale_facette(fa, 1);
              item_to_send_facette(nb_fa7_a_envoyer-1, 14) = normale_facette(fa, 2);

            }
        }
      volume_donnees_facette_envoyees = nb_fa7_a_envoyer;

      for (int som = 0; som < liste_sommets_avant_deplacement.dimension(0); som++)
        {
          if (is_compo_in_proc(compo_connexe_sommets_deplace[som], pe_send_))
            {
              nb_som_a_envoyer++;
              item_to_send_sommet.resize(nb_som_a_envoyer, 6);
              item_to_send_sommet(nb_som_a_envoyer-1, 0) = liste_sommets_avant_deplacement(som, 0);
              item_to_send_sommet(nb_som_a_envoyer-1, 1) = liste_sommets_avant_deplacement(som, 1);
              item_to_send_sommet(nb_som_a_envoyer-1, 2) = liste_sommets_avant_deplacement(som, 2);
              item_to_send_sommet(nb_som_a_envoyer-1, 3) = liste_sommets_apres_deplacement(som, 0);
              item_to_send_sommet(nb_som_a_envoyer-1, 4) = liste_sommets_apres_deplacement(som, 1);
              item_to_send_sommet(nb_som_a_envoyer-1, 5) = liste_sommets_apres_deplacement(som, 2);
            }
        }
      volume_donnees_sommet_envoyees = nb_som_a_envoyer;
    }


  if (pe_send_ >= 0)
    {
      ////////////////////// envoi du volume de donnnes /////////////////////
      send_buffer_volume_facette = new int[1];
      int *buf_volume_facette = send_buffer_volume_facette;
      *buf_volume_facette = volume_donnees_facette_envoyees ;

      ////////////////////// envoi des donnnes elles-memes /////////////////////
      send_buffer_data_facette = new double[volume_donnees_facette_envoyees*dimensions_tab];
      double *buf_data_facette = send_buffer_data_facette;

      for (int i = 0; i < nb_fa7_a_envoyer; i++)
        for (int item = 0; item < dimensions_tab; item++, buf_data_facette++)
          *buf_data_facette = item_to_send_facette(i, item);


      send_buffer_volume_sommet = new int[1];
      int *buf_volume_sommet = send_buffer_volume_sommet;
      *buf_volume_sommet = volume_donnees_sommet_envoyees ;

      ////////////////////// envoi des donnnes elles-memes /////////////////////
      send_buffer_data_sommet = new double[volume_donnees_sommet_envoyees*6];
      double *buf_data_sommet = send_buffer_data_sommet;

      for (int i = 0; i < nb_som_a_envoyer; i++)
        for (int item = 0; item < 6; item++, buf_data_sommet++)
          *buf_data_sommet = item_to_send_sommet(i, item);
    }

  ////////////////////// reception du volume de donnees/////////////////////
  if (pe_recv_ >= 0)
    {
      recv_buffer_volume_facette = new int[1];
      recv_buffer_volume_sommet = new int[1];
    }
  ::envoyer_recevoir(send_buffer_volume_facette, int_size, pe_send_, recv_buffer_volume_facette, int_size, pe_recv_);
  ::envoyer_recevoir(send_buffer_volume_sommet, int_size, pe_send_, recv_buffer_volume_sommet, int_size, pe_recv_);
  if (pe_recv_ >= 0)
    {
      int *buf_volume_facette = recv_buffer_volume_facette;
      volume_donnees_facette_recue = *buf_volume_facette;

      int *buf_volume_sommet = recv_buffer_volume_sommet;
      volume_donnees_sommet_recue = *buf_volume_sommet;
    }

  ////////////////////// reception des donnees elles-memes/////////////////////

  if (pe_recv_ >= 0)
    {
      recv_buffer_data_facette = new double[volume_donnees_facette_recue*dimensions_tab];
      recv_buffer_data_sommet = new double[volume_donnees_sommet_recue*6];
    }
  ::envoyer_recevoir(send_buffer_data_facette, volume_donnees_facette_envoyees*dimensions_tab * double_size, pe_send_, recv_buffer_data_facette, volume_donnees_facette_recue * dimensions_tab * double_size, pe_recv_);
  ::envoyer_recevoir(send_buffer_data_sommet, volume_donnees_sommet_envoyees*6 * double_size, pe_send_, recv_buffer_data_sommet, volume_donnees_sommet_recue * 6 * double_size, pe_recv_);

  if (pe_recv_ >= 0)
    {
      int size_avant_echange = facettes_sommets_full_compo_.dimension(0);
      int new_size = size_avant_echange + volume_donnees_facette_recue ;
      facettes_sommets_full_compo_.resize(new_size, dimensions_tab);
      double *buf_data_facette = recv_buffer_data_facette;
      for (int i = size_avant_echange; i < new_size; i++)
        for (int item = 0; item < dimensions_tab; item++, buf_data_facette++)
          facettes_sommets_full_compo_(i, item) = *buf_data_facette;

      size_avant_echange = liste_sommets_et_deplacements_.dimension(0);
      new_size = size_avant_echange + volume_donnees_sommet_recue ;
      liste_sommets_et_deplacements_.resize(new_size, 6);
      double *buf_data_sommet = recv_buffer_data_sommet;
      for (int i = size_avant_echange; i < new_size; i++)
        for (int item = 0; item < 6; item++, buf_data_sommet++)
          liste_sommets_et_deplacements_(i, item) = *buf_data_sommet;
    }

  delete[] send_buffer_volume_facette;
  delete[] recv_buffer_volume_facette;
  delete[] send_buffer_data_facette;
  delete[] recv_buffer_data_facette;
  delete[] send_buffer_volume_sommet;
  delete[] recv_buffer_volume_sommet;
  delete[] send_buffer_data_sommet;
  delete[] recv_buffer_data_sommet;
}


void FT_Field::calculer_volume_bulles(ArrOfDouble& volumes, DoubleTab& centre_gravite, const Maillage_FT_IJK& mesh) const
{
  const int n = mesh.nb_facettes();
  const ArrOfInt& compo_connex = mesh.compo_connexe_facettes();
  int nb_bulles_reelles = 0 ;
  if (compo_connex.size_array()!=0)
    nb_bulles_reelles = max_array(compo_connex);
  nb_bulles_reelles = Process::mp_max(nb_bulles_reelles) + 1;

  const int nbulles_tot = nb_bulles_reelles ;
  volumes.resize_array(nbulles_tot, RESIZE_OPTIONS::NOCOPY_NOINIT);
  volumes = 0.;
  centre_gravite.resize(nbulles_tot, 3, RESIZE_OPTIONS::NOCOPY_NOINIT);
  centre_gravite = 0.;
  const ArrOfDouble& surfaces_facettes = mesh.get_update_surface_facettes();
  const DoubleTab& normales_facettes = mesh.get_update_normale_facettes();
  const IntTab& facettes = mesh.facettes();
  const DoubleTab& sommets = mesh.sommets();
  const ArrOfInt& compo_facettes = mesh.compo_connexe_facettes();
  for (int i = 0; i < n; i++)
    {
      if (mesh.facette_virtuelle(i))
        continue;
      int compo = compo_facettes[i];
      // les bulles dupliquees a la fin :
      const double s = surfaces_facettes[i];
      const double normale_scalaire_direction = normales_facettes(i, 0); // On projette sur x
      // Coordonnee du centre de gravite de la facette
      const int i0 = facettes(i, 0);
      const int i1 = facettes(i, 1);
      const int i2 = facettes(i, 2);
      const double coord_centre_gravite_i = (sommets(i0, 0) + sommets(i1, 0) + sommets(i2, 0)) / 3.;
      const double coord_centre_gravite_j = (sommets(i0, 1) + sommets(i1, 1) + sommets(i2, 1)) / 3.;
      const double coord_centre_gravite_k = (sommets(i0, 2) + sommets(i1, 2) + sommets(i2, 2)) / 3.;
      const double volume_prisme = coord_centre_gravite_i * s * normale_scalaire_direction;
      // centre de gravite du prisme pondere par son volume avec signe
      centre_gravite(compo, 0) += volume_prisme * (coord_centre_gravite_i * 0.5);
      centre_gravite(compo, 1) += volume_prisme * coord_centre_gravite_j;
      centre_gravite(compo, 2) += volume_prisme * coord_centre_gravite_k;
      volumes[compo] += volume_prisme;
    }
  mp_sum_for_each_item(volumes);
  mp_sum_for_each_item(centre_gravite);
  //Cerr << "volumes : " << volumes << finl;
  for (int i = 0; i < nbulles_tot; i++)
    {
      // const double x = 1./volumes[i];
      const double x = (volumes[i] == 0.) ? 0. : 1. / volumes[i];
      centre_gravite(i, 0) *= x;
      centre_gravite(i, 1) *= x;
      centre_gravite(i, 2) *= x;
    }
}
