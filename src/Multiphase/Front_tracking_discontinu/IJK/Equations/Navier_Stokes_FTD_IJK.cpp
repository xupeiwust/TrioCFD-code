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

#include <Navier_Stokes_FTD_IJK.h>
#include <Probleme_FTD_IJK_base.h>
#include <Fluide_Diphasique_IJK.h>

Implemente_instanciable(Navier_Stokes_FTD_IJK, "Navier_Stokes_FTD_IJK", Equation_base);


Sortie& Navier_Stokes_FTD_IJK::printOn(Sortie& os) const
{
  return os<< que_suis_je() << " " << le_nom();
}

Entree& Navier_Stokes_FTD_IJK::readOn(Entree& is)
{
  return is;
}

void Navier_Stokes_FTD_IJK::associer_pb_base(const Probleme_base& pb)
{
  if (!sub_type(Probleme_FTD_IJK_base, pb))
    {
      Cerr << "Error for the method Navier_Stokes_FTD_IJK::associer_pb_base\n";
      Cerr << " Navier_Stokes_FTD_IJK equation must be associated to\n";
      Cerr << " a Navier_Stokes_FTD_IJK problem type\n";
      Process::exit();
    }
  mon_probleme = pb;
  if (nom_ == "??")
    {
      nom_ = pb.le_nom();
      nom_ += que_suis_je();
    }
}

const Milieu_base& Navier_Stokes_FTD_IJK::milieu() const
{
  if (le_fluide.est_nul())
    {
      Cerr << "You forgot to associate a fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return le_fluide.valeur();
}

Milieu_base& Navier_Stokes_FTD_IJK::milieu()
{
  if (le_fluide.est_nul())
    {
      Cerr << "You forgot to associate a fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return le_fluide.valeur();
}

void Navier_Stokes_FTD_IJK::associer_milieu_base(const Milieu_base& un_milieu)
{
  if (sub_type(Fluide_Diphasique_IJK, un_milieu))
    {
      const Fluide_base& un_fluide = ref_cast(Fluide_base, un_milieu);
      le_fluide = un_fluide;
    }
  else
    {
      Cerr << "Error of fluid type for the method Navier_Stokes_FTD_IJK::associer_milieu_base" << finl;
      Process::exit();
    }
}
