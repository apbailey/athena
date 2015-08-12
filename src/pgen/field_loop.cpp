//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================

// Primary header
#include "../mesh.hpp" 

// C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena headers
#include "../athena.hpp"           // enums, Real
#include "../athena_arrays.hpp"    // AthenaArray
#include "../parameter_input.hpp"  // ParameterInput
#include "../fluid/eos/eos.hpp"    // eos
#include "../fluid/fluid.hpp"      // Fluid
#include "../field/field.hpp"      // Field
#include "../coordinates/coordinates.hpp" // Coordinates

//======================================================================================
//! \file field_loop.c
//  \brief Problem generator for advection of a field loop test. 
//
// Can only be run in 2D or 3D.  Input parameters are:
//   -  problem/rad   = radius of field loop
//   -  problem/amp   = amplitude of vector potential (and therefore B)
//   -  problem/vflow = flow velocity
//   -  problem/drat  = density ratio in loop, to test density advection and conduction
// The flow is automatically set to run along the diagonal. 
//
// Various test cases are possible:
//   - (iprob=1): field loop in x1-x2 plane (cylinder in 3D)
//   - (iprob=2): field loop in x2-x3 plane (cylinder in 3D)
//   - (iprob=3): field loop in x3-x1 plane (cylinder in 3D) 
//   - (iprob=4): rotated cylindrical field loop in 3D.
//   - (iprob=5): spherical field loop in rotated plane
//
// REFERENCE: T. Gardiner & J.M. Stone, "An unsplit Godunov method for ideal MHD via
// constrined transport", JCP, 205, 509 (2005)
//======================================================================================

//======================================================================================
//! \fn ProblemGenerator
//  \brief field loop advection problem generator for 2D/3D problems.
//======================================================================================

void Mesh::ProblemGenerator(Fluid *pfl, Field *pfd, ParameterInput *pin)
{
  MeshBlock *pmb = pfl->pmy_block;
  Coordinates *pco = pmb->pcoord;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  Real gm1 = (pfl->pf_eos->GetGamma() - 1.0);

  AthenaArray<Real> ax,ay,az;
  int nx1 = (ie-is)+1 + 2*(NGHOST);
  int nx2 = (je-js)+1 + 2*(NGHOST);
  int nx3 = (ke-ks)+1 + 2*(NGHOST);
  ax.NewAthenaArray(nx3,nx2,nx1);
  ay.NewAthenaArray(nx3,nx2,nx1);
  az.NewAthenaArray(nx3,nx2,nx1);

// Read initial conditions, diffusion coefficients (if needed)

  Real rad = pin->GetReal("problem","rad");
  Real amp = pin->GetReal("problem","amp");
  Real vflow = pin->GetReal("problem","vflow");
  Real drat = pin->GetOrAddReal("problem","drat",1.0);
  int iprob = pin->GetInteger("problem","iprob");
  Real ang_2,cos_a2,sin_a2,lambda;

// For (iprob=4) -- rotated cylinder in 3D -- set up rotation angle and wavelength

  if(iprob == 4){
    Real x1size = pmb->pmy_mesh->mesh_size.x1max - pmb->pmy_mesh->mesh_size.x1min;
    Real x3size = pmb->pmy_mesh->mesh_size.x3max - pmb->pmy_mesh->mesh_size.x3min;

// We put 1 wavelength in each direction.  Hence the wavelength
//     lambda = x1size*cos_a;
//     AND   lambda = x3size*sin_a;  are both satisfied.

    if(x1size == x3size){
      ang_2 = PI/4.0;
      cos_a2 = sin_a2 = sqrt(0.5);
    } else{
      ang_2 = atan(x1size/x3size);
      sin_a2 = sin(ang_2);
      cos_a2 = cos(ang_2);
    }
    // Use the larger angle to determine the wavelength
    if (cos_a2 >= sin_a2) {
      lambda = x1size*cos_a2;
    } else {
      lambda = x3size*sin_a2;
    }
  }

// Use vector potential to initialize field loop

  for (int k=ks; k<=ke+1; k++) {
  for (int j=js; j<=je+1; j++) {
  for (int i=is; i<=ie+1; i++) {
     
    // (iprob=1): field loop in x1-x2 plane (cylinder in 3D) */
    if(iprob==1) {  
      ax(k,j,i) = 0.0;
      ay(k,j,i) = 0.0;
      if ((SQR(pco->x1f(i)) + SQR(pco->x2f(j))) < rad*rad) {
        az(k,j,i) = amp*(rad - sqrt(SQR(pco->x1f(i)) + SQR(pco->x2f(j))));
      } else {
        az(k,j,i) = 0.0;
      }
    }

    // (iprob=2): field loop in x2-x3 plane (cylinder in 3D)
    if(iprob==2) {  
      if ((SQR(pco->x2f(j)) + SQR(pco->x3f(k))) < rad*rad) {
        ax(k,j,i) = amp*(rad - sqrt(SQR(pco->x2f(j)) + SQR(pco->x3f(k))));
      } else {
        ax(k,j,i) = 0.0;
      }
      ay(k,j,i) = 0.0;
      az(k,j,i) = 0.0;
    }

    // (iprob=3): field loop in x3-x1 plane (cylinder in 3D)
    if(iprob==3) {  
      if ((SQR(pco->x1f(i)) + SQR(pco->x3f(k))) < rad*rad) {
        ay(k,j,i) = amp*(rad - sqrt(SQR(pco->x1f(i)) + SQR(pco->x3f(k))));
      } else {
        ay(k,j,i) = 0.0;
      }
      ax(k,j,i) = 0.0;
      az(k,j,i) = 0.0;
    }

    // (iprob=4): rotated cylindrical field loop in 3D.  Similar to iprob=1 with a 
    // rotation about the x2-axis.  Define coordinate systems (x1,x2,x3) and (x,y,z)
    // with the following transformation rules:
    //    x =  x1*cos(ang_2) + x3*sin(ang_2)
    //    y =  x2
    //    z = -x1*sin(ang_2) + x3*cos(ang_2)
    // This inverts to:
    //    x1  = x*cos(ang_2) - z*sin(ang_2)
    //    x2  = y
    //    x3  = x*sin(ang_2) + z*cos(ang_2)
    if(iprob==4) {
      Real x = pco->x1v(i)*cos_a2 + pco->x3f(k)*sin_a2;
      Real y = pco->x2f(j);
      // shift x back to the domain -0.5*lambda <= x <= 0.5*lambda
      while(x >  0.5*lambda) x -= lambda;
      while(x < -0.5*lambda) x += lambda;
      if ((x*x + y*y) < rad*rad) {
        ax(k,j,i) = amp*(rad - sqrt(x*x + y*y))*(-sin_a2);
      } else {
        ax(k,j,i) = 0.0;
      }
      ay(k,j,i) = 0.0;

      x = pco->x1f(i)*cos_a2 + pco->x3v(k)*sin_a2;
      y = pco->x2f(j);
      // shift x back to the domain -0.5*lambda <= x <= 0.5*lambda
      while(x >  0.5*lambda) x -= lambda;
      while(x < -0.5*lambda) x += lambda;
      if ((x*x + y*y) < rad*rad) {
        az(k,j,i) = amp*(rad - sqrt(x*x + y*y))*(cos_a2);
      } else {
        az(k,j,i) = 0.0;
      }
    }

    // (iprob=5): spherical field loop in rotated plane
    if(iprob==5) { 
      ax(k,j,i) = 0.0;
      if ((SQR(pco->x1f(i)) + SQR(pco->x2v(j)) + SQR(pco->x3f(k))) < rad*rad) {
        ay(k,j,i) = amp*(rad-sqrt(SQR(pco->x1f(i))+SQR(pco->x2v(j))+SQR(pco->x3f(k))));
      } else {
        ay(k,j,i) = 0.0;
      }
      if ((SQR(pco->x1f(i)) + SQR(pco->x2f(j)) + SQR(pco->x3v(k))) < rad*rad) {
        az(k,j,i) = amp*(rad-sqrt(SQR(pco->x1f(i))+SQR(pco->x2f(j))+SQR(pco->x3v(k))));
      } else {
        az(k,j,i) = 0.0;
      }
    }

  }}}

// Initialize density and momenta.  If drat != 1, then density and temperature will be
// different inside loop than background values

  Real x1size = pmb->pmy_mesh->mesh_size.x1max - pmb->pmy_mesh->mesh_size.x1min;
  Real x2size = pmb->pmy_mesh->mesh_size.x2max - pmb->pmy_mesh->mesh_size.x2min;
  Real x3size = pmb->pmy_mesh->mesh_size.x3max - pmb->pmy_mesh->mesh_size.x3min;
  Real diag = sqrt(x1size*x1size + x2size*x2size + x3size*x3size);
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
     pfl->u(IDN,k,j,i) = 1.0;
     pfl->u(IM1,k,j,i) = pfl->u(IDN,k,j,i)*vflow*x1size/diag;
     pfl->u(IM2,k,j,i) = pfl->u(IDN,k,j,i)*vflow*x2size/diag;
     pfl->u(IM3,k,j,i) = pfl->u(IDN,k,j,i)*vflow*x3size/diag;

     if ((SQR(pco->x1v(i)) + SQR(pco->x2v(j)) + SQR(pco->x3v(k))) < rad*rad) {
       pfl->u(IDN,k,j,i) = drat;
       pfl->u(IM1,k,j,i) = pfl->u(IDN,k,j,i)*vflow*x1size/diag;
       pfl->u(IM2,k,j,i) = pfl->u(IDN,k,j,i)*vflow*x2size/diag;
       pfl->u(IM3,k,j,i) = pfl->u(IDN,k,j,i)*vflow*x3size/diag;
     }
  }}}

// initialize interface B

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie+1; i++) {
    pfd->b.x1f(k,j,i) = (az(k,j+1,i) - az(k,j,i))/pco->dx2f(j) -
                        (ay(k+1,j,i) - ay(k,j,i))/pco->dx3f(k);
  }}}
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je+1; j++) {
  for (int i=is; i<=ie; i++) {
    pfd->b.x2f(k,j,i) = (ax(k+1,j,i) - ax(k,j,i))/pco->dx3f(k) -
                        (az(k,j,i+1) - az(k,j,i))/pco->dx1f(i);
  }}}
  for (int k=ks; k<=ke+1; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    pfd->b.x3f(k,j,i) = (ay(k,j,i+1) - ay(k,j,i))/pco->dx1f(i) -
                        (ax(k,j+1,i) - ax(k,j,i))/pco->dx2f(j);
  }}}

// initialize total energy

  if (NON_BAROTROPIC_EOS) {
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        pfl->u(IEN,k,j,i) = 1.0/gm1 +
          0.5*(SQR(0.5*(pfd->b.x1f(k,j,i) + pfd->b.x1f(k,j,i+1))) +
               SQR(0.5*(pfd->b.x2f(k,j,i) + pfd->b.x2f(k,j+1,i))) +
               SQR(0.5*(pfd->b.x3f(k,j,i) + pfd->b.x3f(k+1,j,i)))) + (0.5)*
          (SQR(pfl->u(IM1,k,j,i)) + SQR(pfl->u(IM2,k,j,i)) + SQR(pfl->u(IM3,k,j,i)))
          /pfl->u(IDN,k,j,i);
      }
    }}
  }

  ax.DeleteAthenaArray();
  ay.DeleteAthenaArray();
  az.DeleteAthenaArray();

  return;
}