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

#ifndef POST_PROCESSING_HYDRODYNAMIC_FORCES_H_
#define POST_PROCESSING_HYDRODYNAMIC_FORCES_H_


class Post_Processing_Hydrodynamic_Forces: public Objet_U
{
  Declare_instanciable_sans_constructeur(Post_Processing_Hydrodynamic_Forces);

public:
  friend class Navier_Stokes_FT_Disc;
  friend class Transport_Interfaces_FT_Disc;

  Post_Processing_Hydrodynamic_Forces();

  void associate_transport_equation(OBS_PTR(Transport_Interfaces_FT_Disc)& ptr_eq_transport)
  { ptr_eq_transport_=ptr_eq_transport; }
  void associate_ns_equation(OBS_PTR(Navier_Stokes_FT_Disc)& ptr_eq_ns)
  { ptr_eq_ns_=ptr_eq_ns; }

  void compute_hydrodynamic_forces();

private:
  OBS_PTR(Transport_Interfaces_FT_Disc) ptr_eq_transport_;
  OBS_PTR(Navier_Stokes_FT_Disc) ptr_eq_ns_;
};



#endif
