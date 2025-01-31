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
  if (ni != Sca.nk())
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
