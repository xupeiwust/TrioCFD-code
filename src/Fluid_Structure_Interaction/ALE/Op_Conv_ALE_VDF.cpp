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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Op_Conv_ALE_VDF.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/ALE/src
// Version:     /main/7
//
//////////////////////////////////////////////////////////////////////////////

#include <Op_Conv_ALE_VDF.h>

Implemente_instanciable(Op_Conv_ALE_VDF,"Op_Conv_ALE_VDF",Op_Conv_ALE);

Sortie& Op_Conv_ALE_VDF::printOn(Sortie& os) const
{
  return os;
}


Entree& Op_Conv_ALE_VDF::readOn(Entree& is)
{

  Cerr<<" The VDF discrtisation is not working with the ALE module"<<finl;
  Cerr<<" Please contact the TrioCFD team"<<finl;
  Process::exit();
  return is>>op_conv;
}
DoubleTab& Op_Conv_ALE_VDF::ajouterALE(const DoubleTab& inco, DoubleTab& resu) const
{
  Cerr<<" The VDF discrtisation is not working with the ALE module"<<finl;
  Cerr<<" Please contact the TrioCFD team"<<finl;
  Process::exit();
  return resu;
}
