/****************************************************************************
* Copyright (c) 2025, CEA
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

#include <Probleme_FTD_IJK_tools.h>
#include <Probleme_FTD_IJK_base.h>

void copy_field_values(IJK_Field_double& field, const IJK_Field_double& field_to_copy)
{
  const int ni = field.ni(), nj = field.nj(), nk = field.nk();

  const int ni_to_copy = field_to_copy.ni(), nj_to_copy = field_to_copy.nj(), nk_to_copy = field_to_copy.nk();

  const bool bool_dim = (ni == ni_to_copy && nj == nj_to_copy && nk == nk_to_copy);

  if (!bool_dim)
    Process::exit();

  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        field(i,j,k) = field_to_copy(i,j,k);

  field.echange_espace_virtuel(field.ghost());
}

//  PRODUITS DE CHAMPS
IJK_Field_double scalar_product(const Probleme_FTD_IJK_base& pb, const IJK_Field_vector3_double& V1, const IJK_Field_vector3_double& V2)
{
  /*
   * * ATTENTION : valide pour un maillage cartesien, de maille cubiques uniquement !
   */
  IJK_Field_double resu;
  resu.allocate(pb.domaine_ijk(), Domaine_IJK::ELEM, 3);

  const int nk = V1[0].nk();
  if (nk != V2[0].nk())
    Cerr << "scalar product of fields with different dimensions (nk)" << finl;

  const int nj = V1[0].nj();
  if (nj != V2[0].nj())
    Cerr << "scalar product of fields with different dimensions (nj)" << finl;

  const int ni = V1[0].ni();
  if (ni != V2[0].ni())
    Cerr << "scalar product of fields with different dimensions (ni)" << finl;

  for (int k = 0; k < nk; ++k)
    for (int j = 0; j < nj; ++j)
      for (int i = 0; i < ni; ++i)
        {
          resu(i, j, k) = 0.25 * ((V1[0](i, j, k) + V1[0](i + 1, j, k)) * (V2[0](i, j, k) + V2[0](i + 1, j, k)) +
                                  (V1[1](i, j, k) + V1[1](i, j + 1, k)) * (V2[1](i, j, k) + V2[1](i, j + 1, k))
                                  + (V1[2](i, j, k) + V1[2](i, j, k + 1)) * (V2[2](i, j, k) + V2[2](i, j, k + 1)));
        }
  // Communication avec tous les process ?
  return resu;
}

IJK_Field_vector3_double scalar_times_vector(const Probleme_FTD_IJK_base& pb, const IJK_Field_double& Sca, const IJK_Field_vector3_double& Vec)
{
  /*
   * Produit d'un champ scalaire (Sca) par un champ de vecteur (Vec).
   * Le champ scalaire est aux centre des elements, le champ de vecteur est aux faces
   * Le resultat reste localise au meme endroit que le champ de vecteur passe en entree.
   * ATTENTION : valide pour un maillage cartesien, de maille cubiques uniquement !
   */

  IJK_Field_vector3_double resu;
  allocate_velocity(resu, pb.domaine_ijk(), 3); // j'ai besoin de mettre des cellules ghost ? non, je ne pense pas

  const int nk = Vec[0].nk();
  if (nk != Sca.nk())
    Cerr << "scalar fields has different dimension from vector field  (nk)" << finl;

  const int nj = Vec[0].nj();
  if (nj != Sca.nj())
    Cerr << "scalar fields has different dimension from vector field  (nj)" << finl;

  const int ni = Vec[0].ni();
  if (ni != Sca.ni())
    Cerr << "scalar fields has different dimension from vector field  (ni)" << finl;

  for (int k = 0; k < nk; ++k)
    for (int j = 0; j < nj; ++j)
      for (int i = 0; i < ni; ++i)
        {
          resu[0](i, j, k) = 0.5 * (Sca(i - 1, j, k) + Sca(i, j, k)) * Vec[0](i, j, k);
          resu[1](i, j, k) = 0.5 * (Sca(i, j - 1, k) + Sca(i, j, k)) * Vec[1](i, j, k);
          resu[2](i, j, k) = 0.5 * (Sca(i, j, k - 1) + Sca(i, j, k)) * Vec[2](i, j, k);
        }
  // Communication avec tous les process ?
  return resu;
}

IJK_Field_double scalar_fields_product(const Probleme_FTD_IJK_base& pb, const IJK_Field_double& S1, const IJK_Field_double& S2, int dir)
{
  /*
   * Produit d'un champ scalaire aux centres (S1) par une des composantes d'un champ de vecteur (S2).
   * Le resultat est localise au meme endroit que le champ de vecteur dont est issu S2.
   * ATTENTION : valide pour un maillage cartesien, de maille cubiques uniquement !
   */
  IJK_Field_double resu;
  resu.allocate(pb.domaine_ijk(), Domaine_IJK::ELEM, 3);

  const int nk = S1.nk();
  if (nk != S2.nk())
    Cerr << "scalar fields have different dimensions for the product (nk)" << finl;

  const int nj = S1.nj();
  if (nj != S2.nj())
    Cerr << "scalar fields have different dimensions for the product (nj)" << finl;

  const int ni = S1.ni();
  if (ni != S2.ni())
    Cerr << "scalar fields have different dimensions for the product (ni)" << finl;

  for (int k = 0; k < nk; ++k)
    for (int j = 0; j < nj; ++j)
      for (int i = 0; i < ni; ++i)
        {
          if (dir == 0)
            resu(i, j, k) = 0.5 * (S1(i - 1, j, k) + S1(i, j, k)) * S2(i, j, k);

          if (dir == 1)
            resu(i, j, k) = 0.5 * (S1(i, j - 1, k) + S1(i, j, k)) * S2(i, j, k);

          if (dir == 2)
            resu(i, j, k) = 0.5 * (S1(i, j, k - 1) + S1(i, j, k)) * S2(i, j, k);
        }
  // Communication avec tous les process ?
  return resu;
}

void ecrire_donnees(const Probleme_FTD_IJK_base& pb, const IJK_Field_vector3_double& f3compo, SFichier& le_fichier, const int compo, bool binary)
{
  const IJK_Field_double& f =  f3compo[compo];

  ArrOfDouble coord_i, coord_j, coord_k;
  build_local_coords(f, coord_i, coord_j, coord_k);

  const int ni = f.ni();
  const int nj = f.nj();
  const int nk = f.nk();

  int cnt = 0;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          // le_fichier << coord_i[i] << Separateur::SPACE << coord_j[j] << Separateur::SPACE << coord_k[k] << f(i,j,k) << Separateur::SPACE;
          le_fichier << coord_i[i] << coord_j[j] << coord_k[k] << f(i,j,k);
          cnt++;
        }

  const Domaine_IJK& geom = pb.get_domaine();
  const int idx_min = f.get_domaine().get_offset_local(compo);
  if ((idx_min == 0) &&  (geom.get_periodic_flag(compo)))
    {
      double l = geom.get_domain_length(compo) + geom.get_origin(compo);
      if (compo == 0)
        {
          Cerr << "compo 0: " << l << " " << coord_i[0] << finl;
          for (int k = 0; k < nk; k++)
            for (int j = 0; j < nj; j++)
              {
                le_fichier << l << coord_j[j] << coord_k[k] << f(0,j,k);
                cnt++;
              }
        }
      if (compo == 1)
        {
          Cerr << "compo 1: " << l << " " << coord_j[0] << finl;
          for (int k = 0; k < nk; k++)
            for (int i = 0; i < ni; i++)
              {
                le_fichier << coord_i[i] << l << coord_k[k] << f(i,0,k);
                cnt++;
              }
        }
      if (compo == 2)
        {
          for (int j = 0; j < nj; j++)
            for (int i = 0; i < ni; i++)
              {
                le_fichier << coord_i[i] << coord_j[j] << l << f(i,j,0);
                cnt++;
              }
        }
    }
  le_fichier.flush();
  Cerr << "Written " << cnt << " items with (ni, nj, nk)." << ni << " " << nj << " " << nk << finl;
}

// Initialize field with specified string expression (must be understood by Parser class)
void dumpxyz_vector(const Probleme_FTD_IJK_base& pb, const IJK_Field_vector3_double& f3compo, const char * filename, bool binary)
{
  int np = Process::nproc();
  int rank = Process::me();

  Process::barrier();
  int token = 1;
  if (Process::je_suis_maitre())
    {
      // Write and send token to rank+1
      SFichier le_fichier;
      le_fichier.set_bin(1);
      le_fichier.ouvrir(filename);
      for (unsigned compo=0; compo<3; compo++)
        ecrire_donnees(pb, f3compo, le_fichier, compo, binary);
      if (np > 1)
        envoyer(token, rank, rank+1, 2345);
    }
  else
    {
      int rcv;
      recevoir(rcv, rank-1, rank, 2345);

      SFichier le_fichier;
      le_fichier.set_bin(1);
      le_fichier.ouvrir(filename, ios::app);  // in append mode!
      for (unsigned compo=0; compo<3; compo++)
        ecrire_donnees(pb, f3compo, le_fichier, compo, binary);

      if (rank != np-1)
        envoyer(token, rank, rank+1, 2345);
    }

  Process::barrier();
  Cerr << "Fin de l ecriture dans le fichier XYZ: " << filename << finl;
}
