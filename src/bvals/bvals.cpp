
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
#include "bvals.hpp"

// C++ headers
#include <iostream>   // endl
#include <iomanip>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cstring>    // memcpy

// Athena headers
#include "../athena.hpp"          // Real
#include "../athena_arrays.hpp"   // AthenaArray
#include "../mesh.hpp"            // MeshBlock
#include "../fluid/fluid.hpp"     // Fluid
#include "../parameter_input.hpp" // ParameterInput

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


// arrays of start and end points, created in InitBoundaryBuffer
int fluid_send_se_[6][6];
int fluid_recv_se_[6][6];
int field_send_se_[6][3][6];
int field_recv_se_[6][3][6];
int fluid_bufsize_[6];
int field_bufsize_[6];
int eflux_bufsize_[6];

//======================================================================================
//! \file bvals.cpp
//  \brief implements functions that initialize/apply BCs on each dir
//======================================================================================

// BoundaryValues constructor - sets functions for the appropriate
// boundary conditions at each of the 6 dirs of a MeshBlock

BoundaryValues::BoundaryValues(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_mblock_ = pmb;
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

// Set BC functions for each of the 6 boundaries in turn -------------------------------
// Inner x1
  switch(pmb->block_bcs[inner_x1]){
    case 1:
      FluidBoundary_[inner_x1] = ReflectInnerX1;
      FieldBoundary_[inner_x1] = ReflectInnerX1;
      EFluxBoundary_[inner_x1] = DefaultEFluxInnerX1;
      break;
    case 2:
      FluidBoundary_[inner_x1] = OutflowInnerX1;
      FieldBoundary_[inner_x1] = OutflowInnerX1;
      EFluxBoundary_[inner_x1] = DefaultEFluxInnerX1;
      break;
    case -1: // block boundary
    case 3: // do nothing, useful for user-enrolled BCs
    case 4: // periodic boundary
      FluidBoundary_[inner_x1] = NULL;
      FieldBoundary_[inner_x1] = NULL;
      EFluxBoundary_[inner_x1] = NULL;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "Flag ix1_bc=" << pmb->block_bcs[inner_x1] << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
   }

// Outer x1
  switch(pmb->block_bcs[outer_x1]){
    case 1:
      FluidBoundary_[outer_x1] = ReflectOuterX1;
      FieldBoundary_[outer_x1] = ReflectOuterX1;
      EFluxBoundary_[outer_x1] = DefaultEFluxOuterX1;
      break;
    case 2:
      FluidBoundary_[outer_x1] = OutflowOuterX1;
      FieldBoundary_[outer_x1] = OutflowOuterX1;
      EFluxBoundary_[outer_x1] = DefaultEFluxOuterX1;
      break;
    case -1: // block boundary
    case 3: // do nothing, useful for user-enrolled BCs
    case 4: // periodic boundary
      FluidBoundary_[outer_x1] = NULL;
      FieldBoundary_[outer_x1] = NULL;
      EFluxBoundary_[outer_x1] = NULL;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
          << "Flag ox1_bc=" << pmb->block_bcs[outer_x1] << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }

// Inner x2
  if (pmb->block_size.nx2 > 1) {
    switch(pmb->block_bcs[inner_x2]){
      case 1:
        FluidBoundary_[inner_x2] = ReflectInnerX2;
        FieldBoundary_[inner_x2] = ReflectInnerX2;
        EFluxBoundary_[inner_x2] = DefaultEFluxInnerX2;
        break;
      case 2:
        FluidBoundary_[inner_x2] = OutflowInnerX2;
        FieldBoundary_[inner_x2] = OutflowInnerX2;
        EFluxBoundary_[inner_x2] = DefaultEFluxInnerX2;
        break;
      case -1: // block boundary
      case 3: // do nothing, useful for user-enrolled BCs
      case 4: // periodic boundary
        FluidBoundary_[inner_x2] = NULL;
        FieldBoundary_[inner_x2] = NULL;
        EFluxBoundary_[inner_x2] = NULL;
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ix2_bc=" << pmb->block_bcs[inner_x2] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
     }

// Outer x2
    switch(pmb->block_bcs[outer_x2]){
      case 1:
        FluidBoundary_[outer_x2] = ReflectOuterX2;
        FieldBoundary_[outer_x2] = ReflectOuterX2;
        EFluxBoundary_[outer_x2] = DefaultEFluxOuterX2;
        break;
      case 2:
        FluidBoundary_[outer_x2] = OutflowOuterX2;
        FieldBoundary_[outer_x2] = OutflowOuterX2;
        EFluxBoundary_[outer_x2] = DefaultEFluxOuterX2;
        break;
      case -1: // block boundary
      case 3: // do nothing, useful for user-enrolled BCs
      case 4: // periodic boundary
        FluidBoundary_[outer_x2] = NULL;
        FieldBoundary_[outer_x2] = NULL;
        EFluxBoundary_[outer_x2] = NULL;
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ox2_bc=" << pmb->block_bcs[outer_x2] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
    }
  }

// Inner x3
  if (pmb->block_size.nx3 > 1) {
    switch(pmb->block_bcs[inner_x3]){
      case 1:
        FluidBoundary_[inner_x3] = ReflectInnerX3;
        FieldBoundary_[inner_x3] = ReflectInnerX3;
        EFluxBoundary_[inner_x3] = DefaultEFluxInnerX3;
        break;
      case 2:
        FluidBoundary_[inner_x3] = OutflowInnerX3;
        FieldBoundary_[inner_x3] = OutflowInnerX3;
        EFluxBoundary_[inner_x3] = DefaultEFluxInnerX3;
        break;
      case -1: // block boundary
      case 3: // do nothing, useful for user-enrolled BCs
      case 4: // periodic boundary
        FluidBoundary_[inner_x3] = NULL;
        FieldBoundary_[inner_x3] = NULL;
        EFluxBoundary_[inner_x3] = NULL;
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ix3_bc=" << pmb->block_bcs[inner_x3] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
     }

// Outer x3
    switch(pmb->block_bcs[outer_x3]){
      case 1:
        FluidBoundary_[outer_x3] = ReflectOuterX3;
        FieldBoundary_[outer_x3] = ReflectOuterX3;
        EFluxBoundary_[outer_x3] = DefaultEFluxOuterX3;
        break;
      case 2:
        FluidBoundary_[outer_x3] = OutflowOuterX3;
        FieldBoundary_[outer_x3] = OutflowOuterX3;
        EFluxBoundary_[outer_x3] = DefaultEFluxOuterX3;
        break;
      case -1: // block boundary
      case 3: // do nothing, useful for user-enrolled BCs
      case 4: // periodic boundary
        FluidBoundary_[outer_x3] = NULL;
        FieldBoundary_[outer_x3] = NULL;
        EFluxBoundary_[outer_x3] = NULL;
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in BoundaryValues constructor" << std::endl
            << "Flag ox3_bc=" << pmb->block_bcs[outer_x3] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
    }
  }

  // Allocate Buffers
  int r=2;
  if(pmb->block_size.nx2 > 1) r=4;
  if(pmb->block_size.nx3 > 1) r=6;
  for(int i=0;i<r;i++) {
    fluid_send_[i]=new Real[fluid_bufsize_[i]];
    fluid_recv_[i]=new Real[fluid_bufsize_[i]];
  }

  if (MAGNETIC_FIELDS_ENABLED) {
    for(int i=0;i<r;i++) {
      field_send_[i]=new Real[field_bufsize_[i]];
      field_recv_[i]=new Real[field_bufsize_[i]];
    }
    if(pmb->block_size.nx2>1) { // 2D or 3D
      for(int i=0;i<r;i++) {
        eflux_send_[i]=new Real[eflux_bufsize_[i]];
        eflux_recv_[i]=new Real[eflux_bufsize_[i]];
      }
    }
  }

  // initialize flags
  for(int i=0;i<6;i++) {
    fluid_flag_[i][0][0]=false;
    fluid_flag_[i][0][1]=false;
    fluid_flag_[i][1][0]=false;
    fluid_flag_[i][1][1]=false;
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int i=0;i<6;i++) {
      field_flag_[i][0][0]=false;
      field_flag_[i][0][1]=false;
      field_flag_[i][1][0]=false;
      field_flag_[i][1][1]=false;
      eflux_flag_[i][0][0]=false;
      eflux_flag_[i][0][1]=false;
      eflux_flag_[i][1][0]=false;
      eflux_flag_[i][1][1]=false;
    }
  }
}

// destructor

BoundaryValues::~BoundaryValues()
{
  int r=2;
  if(pmy_mblock_->block_size.nx2 > 1) r=4;
  if(pmy_mblock_->block_size.nx3 > 1) r=6;
  for(int i=0;i<r;i++) {
    delete [] fluid_send_[i];
    delete [] fluid_recv_[i];
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int i=0;i<r;i++) { 
      delete [] field_send_[i];
      delete [] field_recv_[i];
    }
    if(pmy_mblock_->block_size.nx2 > 1) {
      for(int i=0;i<r;i++) { 
        delete [] eflux_send_[i];
        delete [] eflux_recv_[i];
      }
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::EnrollFluidBoundaryFunction(enum direction dir,
//                                                       BValFluid_t my_bc)
//  \brief Enroll a user-defined boundary function for fluid

void BoundaryValues::EnrollFluidBoundaryFunction(enum direction dir, BValFluid_t my_bc)
{
  std::stringstream msg;
  if(dir<0 || dir>5)
  {
    msg << "### FATAL ERROR in EnrollFluidBoundaryCondition function" << std::endl
        << "dirName = " << dir << " not valid" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(pmy_mblock_->pmy_mesh->mesh_bcs[dir]!=3) {
    msg << "### FATAL ERROR in EnrollFluidBoundaryCondition function" << std::endl
        << "A user-defined boundary condition flag (3) must be specified "
        << "in the input file to use a user-defined boundary function." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(pmy_mblock_->neighbor[dir][0][0].gid==-1)
    FluidBoundary_[dir]=my_bc;
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::EnrollFieldBoundaryFunction(enum direction dir,
//                                                       BValField_t my_bc)
//  \brief Enroll a user-defined boundary function for magnetic fields

void BoundaryValues::EnrollFieldBoundaryFunction(enum direction dir,BValField_t my_bc)
{
  std::stringstream msg;
  if(dir<0 || dir>5)
  {
    msg << "### FATAL ERROR in EnrollFieldBoundaryCondition function" << std::endl
        << "dirName = " << dir << " is not valid" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(pmy_mblock_->pmy_mesh->mesh_bcs[dir]!=3) {
    msg << "### FATAL ERROR in EnrollFieldBoundaryCondition function" << std::endl
        << "A user-defined boundary condition flag (3) must be specified "
        << "in the input file to use a user-defined boundary function." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(pmy_mblock_->neighbor[dir][0][0].gid==-1)
    FieldBoundary_[dir]=my_bc;
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::EnrollEFluxBoundaryFunction(enum direction dir,
//                                                       BValEFlux_t my_bc)
//  \brief Enroll a user-defined boundary function for electric fields

void BoundaryValues::EnrollEFluxBoundaryFunction(enum direction dir,BValEFlux_t my_bc)
{
  std::stringstream msg;
  if(dir<0 || dir>5)
  {
    msg << "### FATAL ERROR in EnrollEFluxBoundaryCondition function" << std::endl
        << "dirName = " << dir << " is not valid" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(pmy_mblock_->pmy_mesh->mesh_bcs[dir]!=3) {
    msg << "### FATAL ERROR in EnrollEFluxBoundaryCondition function" << std::endl
        << "A user-defined boundary condition flag (3) must be specified "
        << "in the input file to use a user-defined boundary function." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if(pmy_mblock_->neighbor[dir][0][0].gid==-1)
    EFluxBoundary_[dir]=my_bc;
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::StartReceivingFluid(int flag)
//  \brief initiate MPI_Irecv for fluid
void BoundaryValues::StartReceivingFluid(int flag)
{
#ifdef MPI_PARALLEL
  MeshBlock *pmb=pmy_mblock_;
  int tag;
  for(int i=0;i<6;i++) {
    if(pmb->neighbor[i][0][0].gid!=-1 && pmb->neighbor[i][0][0].rank!=myrank) { 
      tag=CreateMPITag(pmb->lid, flag, i, tag_fluid, 0, 0);
      MPI_Irecv(fluid_recv_[i],fluid_bufsize_[i],MPI_ATHENA_REAL,
        pmb->neighbor[i][0][0].rank,tag,MPI_COMM_WORLD,&req_fluid_recv_[i][0][0]);
    }
  }
#endif
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::StartReceivingField(int flag)
//  \brief initiate MPI_Irecv for field
void BoundaryValues::StartReceivingField(int flag)
{
#ifdef MPI_PARALLEL
  MeshBlock *pmb=pmy_mblock_;
  int tag;
  for(int i=0;i<6;i++) {
    if(pmb->neighbor[i][0][0].gid!=-1 && pmb->neighbor[i][0][0].rank!=myrank) {
      tag=CreateMPITag(pmb->lid, flag, i, tag_field, 0, 0);
      MPI_Irecv(field_recv_[i],field_bufsize_[i],MPI_ATHENA_REAL,
          pmb->neighbor[i][0][0].rank,tag,MPI_COMM_WORLD,&req_field_recv_[i][0][0]);
    }
  }
#endif
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::StartReceivingEFlux(int flag)
//  \brief initiate MPI_Irecv for field
void BoundaryValues::StartReceivingEFlux(int flag)
{
  MeshBlock *pmb=pmy_mblock_;
  int tag, ndir;
  if(pmb->block_size.nx2==1)
    return; // 1D
  ndir=4; // 2D
  if(pmb->block_size.nx3>1)
    ndir=6; // 3D
  for(int i=0;i<ndir;i++) {
    eflux_flag_[i][0][0]=0;
#ifdef MPI_PARALLEL
    if(pmb->neighbor[i][0][0].gid!=-1 && pmb->neighbor[i][0][0].rank!=myrank) {
      tag=CreateMPITag(pmb->lid, flag, i, tag_eflux, 0, 0);
      MPI_Irecv(eflux_recv_[i],eflux_bufsize_[i],MPI_ATHENA_REAL,
          pmb->neighbor[i][0][0].rank,tag,MPI_COMM_WORLD,&req_eflux_recv_[i][0][0]);
    }
#endif
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::LoadAndSendFluidBoundaryBuffer
//                          (enum direction dir, AthenaArray<Real> &src, int flag)
//  \brief Set boundary buffer for x1 direction using boundary functions
//  note: some geometric boundaries (e.g. origin and pole) are not implemented yet
void BoundaryValues::LoadAndSendFluidBoundaryBuffer
                     (enum direction dir, AthenaArray<Real> &src, int flag)
{
  MeshBlock *pmb=pmy_mblock_;
  MeshBlock *pbl=pmb->pmy_mesh->pblock;
  int oside;
  Real *sendbuf=fluid_send_[dir];
  int si, sj, sk, ei, ej, ek, mylevel;
#ifdef MPI_PARALLEL
  int tag;
#endif

  if(pmb->neighbor[dir][0][0].gid==-1)
    return; // do nothing for physical boundary

  si=fluid_send_se_[dir][0];
  ei=fluid_send_se_[dir][1];
  sj=fluid_send_se_[dir][2];
  ej=fluid_send_se_[dir][3];
  sk=fluid_send_se_[dir][4];
  ek=fluid_send_se_[dir][5];

  if(dir%2==0)
    oside=dir+1;
  else
    oside=dir-1;

  // Set buffers
  int p=0;
  for (int n=0; n<(NFLUID); ++n) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma simd
        for (int i=si; i<=ei; ++i) {
          // buffer is always fully packed
          sendbuf[p++]=src(n,k,j,i);
        }
      }
    }
  }

  // Send the buffer; modify this for MPI and AMR
  if(pmb->neighbor[dir][0][0].rank == myrank) // myrank
  {
    while(pbl!=NULL)
    {
      if(pbl->gid==pmb->neighbor[dir][0][0].gid)
        break;
      pbl=pbl->next;
    }
    std::memcpy(pbl->pbval->fluid_recv_[oside], sendbuf,
                fluid_bufsize_[dir]*sizeof(Real));
    pbl->pbval->fluid_flag_[oside][0][0]=true; // the other side
  }
  else // MPI
  {
#ifdef MPI_PARALLEL
    // on the same level
    tag=CreateMPITag(pmb->neighbor[dir][0][0].lid, flag, oside, tag_fluid, 0, 0);
    MPI_Isend(sendbuf,fluid_bufsize_[dir],MPI_ATHENA_REAL,
      pmb->neighbor[dir][0][0].rank,tag,MPI_COMM_WORLD,&req_fluid_send_[dir][0][0]);
#endif
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveAndSetFluidBoundary(enum direction dir,
//                                                      AthenaArray<Real> &dst)
//  \brief load boundary buffer for x1 direction into the array
bool BoundaryValues::ReceiveAndSetFluidBoundary(enum direction dir,
                                                AthenaArray<Real> &dst)
{
  MeshBlock *pmb=pmy_mblock_;
  Real *recvbuf=fluid_recv_[dir];
  int si, sj, sk, ei, ej, ek;

  if(pmb->neighbor[dir][0][0].gid==-1) // physical boundary
    FluidBoundary_[dir](pmb,dst);
  else // block boundary
  {
#ifdef MPI_PARALLEL
    if(fluid_flag_[dir][0][0] == false)
    {
      if(pmb->neighbor[dir][0][0].rank!=myrank) // MPI boundary
        MPI_Wait(&req_fluid_recv_[dir][0][0],MPI_STATUS_IGNORE);
    }
#endif

    si=fluid_recv_se_[dir][0];
    ei=fluid_recv_se_[dir][1];
    sj=fluid_recv_se_[dir][2];
    ej=fluid_recv_se_[dir][3];
    sk=fluid_recv_se_[dir][4];
    ek=fluid_recv_se_[dir][5];

    int p=0;
    for (int n=0; n<(NFLUID); ++n) {
      for (int k=sk; k<=ek; ++k) {
        for (int j=sj; j<=ej; ++j) {
#pragma simd
          for (int i=si; i<=ei; ++i) {
            // buffer is always fully packed
            dst(n,k,j,i) = recvbuf[p++];
          }
        }
      }
    }
    fluid_flag_[dir][0][0] = false; // clear the flag
  }
  return true;
}

//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::WaitSendFluid(enum direction dir)
//  \brief wait until MPI_Isend for fluid completes
void BoundaryValues::WaitSendFluid(enum direction dir)
{
#ifdef MPI_PARALLEL
  MeshBlock *pmb=pmy_mblock_;
  if(pmb->neighbor[dir][0][0].rank!=-1 && pmb->neighbor[dir][0][0].rank!=myrank)
    MPI_Wait(&req_fluid_send_[dir][0][0],MPI_STATUS_IGNORE);
#endif
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::LoadAndSendFieldBoundaryBuffer
//                           (enum direction dir, InterfaceField &src, int flag)
//  \brief Set boundary buffer for x1 direction using boundary functions
//  note: some geometric boundaries (e.g. origin and pole) are not implemented yet
void BoundaryValues::LoadAndSendFieldBoundaryBuffer(enum direction dir,
                                                    InterfaceField &src, int flag)
{
  MeshBlock *pmb=pmy_mblock_;
  MeshBlock *pbl=pmb->pmy_mesh->pblock;
  int oside;
  Real *sendbuf=field_send_[dir];
  AthenaArray<Real>& x1src=src.x1f;
  AthenaArray<Real>& x2src=src.x2f;
  AthenaArray<Real>& x3src=src.x3f;
  int si, sj, sk, ei, ej, ek;
#ifdef MPI_PARALLEL
  int tag;
#endif

  if(pmb->neighbor[dir][0][0].gid==-1)
    return; // do nothing for physical boundary

  if(dir%2==0)
    oside=dir+1;
  else
    oside=dir-1;

  // Set buffers; x1f
  int p=0;
  si=field_send_se_[dir][x1face][0];
  ei=field_send_se_[dir][x1face][1];
  sj=field_send_se_[dir][x1face][2];
  ej=field_send_se_[dir][x1face][3];
  sk=field_send_se_[dir][x1face][4];
  ek=field_send_se_[dir][x1face][5];
  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i) {
        // buffer is always fully packed
        sendbuf[p++]=x1src(k,j,i);
      }
    }
  }
  // Set buffers; x2f
  si=field_send_se_[dir][x2face][0];
  ei=field_send_se_[dir][x2face][1];
  sj=field_send_se_[dir][x2face][2];
  ej=field_send_se_[dir][x2face][3];
  sk=field_send_se_[dir][x2face][4];
  ek=field_send_se_[dir][x2face][5];
  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i) {
        // buffer is always fully packed
        sendbuf[p++]=x2src(k,j,i);
      }
    }
  }
  // Set buffers; x3f
  si=field_send_se_[dir][x3face][0];
  ei=field_send_se_[dir][x3face][1];
  sj=field_send_se_[dir][x3face][2];
  ej=field_send_se_[dir][x3face][3];
  sk=field_send_se_[dir][x3face][4];
  ek=field_send_se_[dir][x3face][5];
  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma simd
      for (int i=si; i<=ei; ++i) {
        // buffer is always fully packed
        sendbuf[p++]=x3src(k,j,i);
      }
    }
  }

  // Send the buffer; modify this for MPI and AMR
  if(pmb->neighbor[dir][0][0].rank == myrank) // myrank
  {
    while(pbl!=NULL)
    {
      if(pbl->gid==pmb->neighbor[dir][0][0].gid)
        break;
      pbl=pbl->next;
    }
    std::memcpy(pbl->pbval->field_recv_[oside], field_send_[dir],
                field_bufsize_[dir]*sizeof(Real));
    pbl->pbval->field_flag_[oside][0][0]=true; // the other side
  }
  else // MPI
  {
#ifdef MPI_PARALLEL
    // on the same level
    tag=CreateMPITag(pmb->neighbor[dir][0][0].lid, flag, oside, tag_field, 0, 0);
    MPI_Isend(sendbuf,field_bufsize_[dir],MPI_ATHENA_REAL,
        pmb->neighbor[dir][0][0].rank,tag,MPI_COMM_WORLD,&req_field_send_[dir][0][0]);
#endif
  }
  return;
}

//--------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveAndSetFieldBoundary(enum direction dir,
//                                                      InterfaceField &dst)
//  \brief load boundary buffer for x1 direction into the array
bool BoundaryValues::ReceiveAndSetFieldBoundary(enum direction dir, InterfaceField &dst)
{
  MeshBlock *pmb=pmy_mblock_;
  Real *recvbuf=field_recv_[dir];
  AthenaArray<Real>& x1dst=dst.x1f;
  AthenaArray<Real>& x2dst=dst.x2f;
  AthenaArray<Real>& x3dst=dst.x3f;
  int si, sj, sk, ei, ej, ek;

  if(pmb->neighbor[dir][0][0].gid==-1) // physical boundary
    FieldBoundary_[dir](pmb,dst);
  else // block boundary
  {
#ifdef MPI_PARALLEL
    if(field_flag_[dir][0][0] == false)
    {
      if(pmb->neighbor[dir][0][0].rank!=myrank) { // MPI boundary
        // temporary: wait the communication
        MPI_Wait(&req_field_recv_[dir][0][0],MPI_STATUS_IGNORE);

        // for the future interleaving: check if it is ready
        //   return false; return if it is not ready yet ; for task implementation
      }
    }
#endif

    // Load buffers; x1f
    int p=0;
    si=field_recv_se_[dir][x1face][0];
    ei=field_recv_se_[dir][x1face][1];
    sj=field_recv_se_[dir][x1face][2];
    ej=field_recv_se_[dir][x1face][3];
    sk=field_recv_se_[dir][x1face][4];
    ek=field_recv_se_[dir][x1face][5];
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma simd
        for (int i=si; i<=ei; ++i) {
          // buffer is always fully packed
          x1dst(k,j,i)=recvbuf[p++];
        }
      }
    }
    // Load buffers; x2f
    si=field_recv_se_[dir][x2face][0];
    ei=field_recv_se_[dir][x2face][1];
    sj=field_recv_se_[dir][x2face][2];
    ej=field_recv_se_[dir][x2face][3];
    sk=field_recv_se_[dir][x2face][4];
    ek=field_recv_se_[dir][x2face][5];
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma simd
        for (int i=si; i<=ei; ++i) {
          // buffer is always fully packed
          x2dst(k,j,i)=recvbuf[p++];
        }
      }
    }
    // Load buffers; x3f
    si=field_recv_se_[dir][x3face][0];
    ei=field_recv_se_[dir][x3face][1];
    sj=field_recv_se_[dir][x3face][2];
    ej=field_recv_se_[dir][x3face][3];
    sk=field_recv_se_[dir][x3face][4];
    ek=field_recv_se_[dir][x3face][5];
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma simd
        for (int i=si; i<=ei; ++i) {
          // buffer is always fully packed
          x3dst(k,j,i)=recvbuf[p++];
        }
      }
    }
    field_flag_[dir][0][0] = false; // clear the flag
  }

  return true;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::WaitSendField(enum direction dir)
//  \brief wait until MPI_Isend for magnetic fields completes
void BoundaryValues::WaitSendField(enum direction dir)
{
#ifdef MPI_PARALLEL
  MeshBlock *pmb=pmy_mblock_;
  if(pmb->neighbor[dir][0][0].rank!=-1 && pmb->neighbor[dir][0][0].rank!=myrank)
    MPI_Wait(&req_field_send_[dir][0][0],MPI_STATUS_IGNORE);
#endif
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::LoadAndSendEFluxBoundaryBuffer(InterfaceField &fsrc,
//                                                      InterfaceField &wsrc, int flag)
//  \brief Set boundary buffer for x1 direction using boundary functions
//  note: some geometric boundaries (e.g. origin and pole) are not implemented yet
void BoundaryValues::LoadAndSendEFluxBoundaryBuffer(InterfaceField &fsrc,
                                                    InterfaceField &wsrc, int flag)
{
  MeshBlock *pmb=pmy_mblock_;
  int oside, p, dir, ndir;
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;
#ifdef MPI_PARALLEL
  int tag;
#endif

  if(pmb->block_size.nx2==1)
    return; // 1D
  if(pmb->block_size.nx3==1) { // 2D
    if(pmb->neighbor[inner_x1][0][0].gid!=-1) {
      p=0;
      for(int j=js; j<=je+1; j++)
        eflux_send_[inner_x1][p++]=fsrc.x2f(X2E3,ks,j,is);
      for(int j=js; j<=je+1; j++)
        eflux_send_[inner_x1][p++]=wsrc.x2f(ks,j,is);
    }
    if(pmb->neighbor[outer_x1][0][0].gid!=-1) {
      p=0;
      for(int j=js; j<=je+1; j++)
        eflux_send_[outer_x1][p++]=fsrc.x2f(X2E3,ks,j,ie);
      for(int j=js; j<=je+1; j++)
        eflux_send_[outer_x1][p++]=wsrc.x2f(ks,j,ie);
    }
    if(pmb->neighbor[inner_x2][0][0].gid!=-1) {
      p=0;
      for(int i=is; i<=ie+1; i++)
        eflux_send_[inner_x2][p++]=fsrc.x1f(X1E3,ks,js,i);
      for(int i=is; i<=ie+1; i++)
        eflux_send_[inner_x2][p++]=wsrc.x1f(ks,js,i);
    }
    if(pmb->neighbor[outer_x2][0][0].gid!=-1) {
      p=0;
      for(int i=is; i<=ie+1; i++)
        eflux_send_[outer_x2][p++]=fsrc.x1f(X1E3,ks,je,i);
      for(int i=is; i<=ie+1; i++)
        eflux_send_[outer_x2][p++]=wsrc.x1f(ks,je,i);
    }
  }
  else {  // 3D
    if(pmb->neighbor[inner_x1][0][0].gid!=-1) {
      p=0;
      for(int k=ks; k<=ke; k++) {
        for(int j=js; j<=je+1; j++)
          eflux_send_[inner_x1][p++]=fsrc.x2f(X2E3,k,j,is);
      }
      for(int k=ks; k<=ke; k++) {
        for(int j=js; j<=je+1; j++)
          eflux_send_[inner_x1][p++]=wsrc.x2f(k,j,is);
      }
      for(int k=ks; k<=ke+1; k++) {
        for(int j=js; j<=je; j++)
          eflux_send_[inner_x1][p++]=fsrc.x3f(X3E2,k,j,is);
      }
      for(int k=ks; k<=ke+1; k++) {
        for(int j=js; j<=je; j++)
          eflux_send_[inner_x1][p++]=wsrc.x3f(k,j,is);
      }
    }
    if(pmb->neighbor[outer_x1][0][0].gid!=-1) {
      p=0;
      for(int k=ks; k<=ke; k++) {
        for(int j=js; j<=je+1; j++)
          eflux_send_[outer_x1][p++]=fsrc.x2f(X2E3,k,j,ie);
      }
      for(int k=ks; k<=ke; k++) {
        for(int j=js; j<=je+1; j++)
          eflux_send_[outer_x1][p++]=wsrc.x2f(k,j,ie);
      }
      for(int k=ks; k<=ke+1; k++) {
        for(int j=js; j<=je; j++)
          eflux_send_[outer_x1][p++]=fsrc.x3f(X3E2,k,j,ie);
      }
      for(int k=ks; k<=ke+1; k++) {
        for(int j=js; j<=je; j++)
          eflux_send_[outer_x1][p++]=wsrc.x3f(k,j,ie);
      }
    }
    if(pmb->neighbor[inner_x2][0][0].gid!=-1) {
      p=0;
      for(int k=ks; k<=ke; k++) {
        for(int i=is; i<=ie+1; i++)
          eflux_send_[inner_x2][p++]=fsrc.x1f(X1E3,k,js,i);
      }
      for(int k=ks; k<=ke; k++) {
        for(int i=is; i<=ie+1; i++)
          eflux_send_[inner_x2][p++]=wsrc.x1f(k,js,i);
      }
      for(int k=ks; k<=ke+1; k++) {
        for(int i=is; i<=ie; i++)
          eflux_send_[inner_x2][p++]=fsrc.x3f(X3E1,k,js,i);
      }
      for(int k=ks; k<=ke+1; k++) {
        for(int i=is; i<=ie; i++)
          eflux_send_[inner_x2][p++]=wsrc.x3f(k,js,i);
      }
    }
    if(pmb->neighbor[outer_x2][0][0].gid!=-1) {
      p=0;
      for(int k=ks; k<=ke; k++) {
        for(int i=is; i<=ie+1; i++)
          eflux_send_[outer_x2][p++]=fsrc.x1f(X1E3,k,je,i);
      }
      for(int k=ks; k<=ke; k++) {
        for(int i=is; i<=ie+1; i++)
          eflux_send_[outer_x2][p++]=wsrc.x1f(k,je,i);
      }
      for(int k=ks; k<=ke+1; k++) {
        for(int i=is; i<=ie; i++)
          eflux_send_[outer_x2][p++]=fsrc.x3f(X3E1,k,je,i);
      }
      for(int k=ks; k<=ke+1; k++) {
        for(int i=is; i<=ie; i++)
          eflux_send_[outer_x2][p++]=wsrc.x3f(k,je,i);
      }
    }
    if(pmb->neighbor[inner_x3][0][0].gid!=-1) {
      p=0;
      for(int j=js; j<=je; j++) {
        for(int i=is; i<=ie+1; i++)
          eflux_send_[inner_x3][p++]=fsrc.x1f(X1E2,ks,j,i);
      }
      for(int j=js; j<=je; j++) {
        for(int i=is; i<=ie+1; i++)
          eflux_send_[inner_x3][p++]=wsrc.x1f(ks,j,i);
      }
      for(int j=js; j<=je+1; j++) {
        for(int i=is; i<=ie; i++)
          eflux_send_[inner_x3][p++]=fsrc.x2f(X2E1,ks,j,i);
      }
      for(int j=js; j<=je+1; j++) {
        for(int i=is; i<=ie; i++)
          eflux_send_[inner_x3][p++]=wsrc.x2f(ks,j,i);
      }
    }
    if(pmb->neighbor[outer_x3][0][0].gid!=-1) {
      p=0;
      for(int j=js; j<=je; j++) {
        for(int i=is; i<=ie+1; i++)
          eflux_send_[outer_x3][p++]=fsrc.x1f(X1E2,ke,j,i);
      }
      for(int j=js; j<=je; j++) {
        for(int i=is; i<=ie+1; i++)
          eflux_send_[outer_x3][p++]=wsrc.x1f(ke,j,i);
      }
      for(int j=js; j<=je+1; j++) {
        for(int i=is; i<=ie; i++)
          eflux_send_[outer_x3][p++]=fsrc.x2f(X2E1,ke,j,i);
      }
      for(int j=js; j<=je+1; j++) {
        for(int i=is; i<=ie; i++)
          eflux_send_[outer_x3][p++]=wsrc.x2f(ke,j,i);
      }
    }
  }

  ndir=4;
  if(pmb->block_size.nx3>1) // 3D
    ndir=6;
  for(dir=0;dir<ndir;dir++)
  {
    if(pmb->neighbor[dir][0][0].gid==-1)
      continue; // do nothing for physical boundary
    if(dir%2==0)
      oside=dir+1;
    else
      oside=dir-1;
    // Send the buffer; modify this for MPI and AMR
    if(pmb->neighbor[dir][0][0].rank == myrank) // myrank
    {
      MeshBlock *pbl=pmb->pmy_mesh->pblock;
      while(pbl!=NULL)
      {
        if(pbl->gid==pmb->neighbor[dir][0][0].gid)
          break;
        pbl=pbl->next;
      }
      std::memcpy(pbl->pbval->eflux_recv_[oside], eflux_send_[dir],
                  eflux_bufsize_[dir]*sizeof(Real));
      pbl->pbval->eflux_flag_[oside][0][0]=true; // the other side
    }
    else // MPI
    {
#ifdef MPI_PARALLEL
      // on the same level
      tag=CreateMPITag(pmb->neighbor[dir][0][0].lid, flag, oside, tag_eflux, 0, 0);
      MPI_Isend(eflux_send_[dir],eflux_bufsize_[dir],MPI_ATHENA_REAL,
          pmb->neighbor[dir][0][0].rank,tag,MPI_COMM_WORLD,&req_eflux_send_[dir][0][0]);
#endif
    }
  }
  return;
}

//--------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveAndSetEFluxBoundary(InterfaceField &fdst,
//                                                      InterfaceField &wdst)
//  \brief load boundary buffer for x1 direction into the array
bool BoundaryValues::ReceiveAndSetEFluxBoundary(InterfaceField &fdst,
                                                InterfaceField &wdst)
{
  MeshBlock *pmb=pmy_mblock_;
  int dir, ndir, p;
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  if(pmb->block_size.nx2==1) // 1D
    return true;
#ifdef MPI_PARALLEL
  ndir=4;
  if(pmb->block_size.nx3>1) // 3D
    ndir=6;
  for(dir=0;dir<ndir;dir++) {
    if(pmb->neighbor[dir][0][0].gid!=-1 && pmb->neighbor[dir][0][0].rank!=myrank) {
      MPI_Wait(&req_eflux_recv_[dir][0][0],MPI_STATUS_IGNORE);
      eflux_flag_[dir][0][0]=true;
    }
  }
#endif
  // at this point, all the buffers should be completed

  if(pmb->block_size.nx3==1) { // 2D
    if(pmb->neighbor[inner_x1][0][0].gid==-1) // physical boundary
      EFluxBoundary_[inner_x1](pmb,fdst,wdst);
    else {
      p=0;
      for(int j=js; j<=je+1; j++)
        fdst.x2f(X2E3,ks,j,is-1)=eflux_recv_[inner_x1][p++];
      for(int j=js; j<=je+1; j++)
        wdst.x2f(ks,j,is-1)=eflux_recv_[inner_x1][p++];
    }
    eflux_flag_[inner_x1][0][0] = false; // clear the flag
    if(pmb->neighbor[outer_x1][0][0].gid==-1) // physical boundary
      EFluxBoundary_[outer_x1](pmb,fdst,wdst);
    else {
      p=0;
      for(int j=js; j<=je+1; j++)
        fdst.x2f(X2E3,ks,j,ie+1)=eflux_recv_[outer_x1][p++];
      for(int j=js; j<=je+1; j++)
        wdst.x2f(ks,j,ie+1)=eflux_recv_[outer_x1][p++];
    }
    eflux_flag_[outer_x1][0][0] = false; // clear the flag
    if(pmb->neighbor[inner_x2][0][0].gid==-1) // physical boundary
      EFluxBoundary_[inner_x2](pmb,fdst,wdst);
    else {
      p=0;
      for(int i=is; i<=ie+1; i++)
        fdst.x1f(X1E3,ks,js-1,i)=eflux_recv_[inner_x2][p++];
      for(int i=is; i<=ie+1; i++)
        wdst.x1f(ks,js-1,i)=eflux_recv_[inner_x2][p++];
    }
    eflux_flag_[inner_x2][0][0] = false; // clear the flag
    if(pmb->neighbor[outer_x2][0][0].gid==-1) // physical boundary
      EFluxBoundary_[outer_x2](pmb,fdst,wdst);
    else {
      p=0;
      for(int i=is; i<=ie+1; i++)
        fdst.x1f(X1E3,ks,je+1,i)=eflux_recv_[outer_x2][p++];
      for(int i=is; i<=ie+1; i++)
        wdst.x1f(ks,je+1,i)=eflux_recv_[outer_x2][p++];
    }
    eflux_flag_[outer_x2][0][0] = false; // clear the flag
  }
  else { // 3D
    if(pmb->neighbor[inner_x1][0][0].gid==-1) // physical boundary
      EFluxBoundary_[inner_x1](pmb,fdst,wdst);
    else {
      p=0;
      for(int k=ks; k<=ke; k++) {
        for(int j=js; j<=je+1; j++)
          fdst.x2f(X2E3,k,j,is-1)=eflux_recv_[inner_x1][p++];
      }
      for(int k=ks; k<=ke; k++) {
        for(int j=js; j<=je+1; j++)
          wdst.x2f(k,j,is-1)=eflux_recv_[inner_x1][p++];
      }
      for(int k=ks; k<=ke+1; k++) {
        for(int j=js; j<=je; j++)
          fdst.x3f(X3E2,k,j,is-1)=eflux_recv_[inner_x1][p++];
      }
      for(int k=ks; k<=ke+1; k++) {
        for(int j=js; j<=je; j++)
          wdst.x3f(k,j,is-1)=eflux_recv_[inner_x1][p++];
      }
    }
    eflux_flag_[inner_x1][0][0] = false; // clear the flag
    if(pmb->neighbor[outer_x1][0][0].gid==-1) // physical boundary
      EFluxBoundary_[outer_x1](pmb,fdst,wdst);
    else {
      p=0;
      for(int k=ks; k<=ke; k++) {
        for(int j=js; j<=je+1; j++)
          fdst.x2f(X2E3,k,j,ie+1)=eflux_recv_[outer_x1][p++];
      }
      for(int k=ks; k<=ke; k++) {
        for(int j=js; j<=je+1; j++)
          wdst.x2f(k,j,ie+1)=eflux_recv_[outer_x1][p++];
      }
      for(int k=ks; k<=ke+1; k++) {
        for(int j=js; j<=je; j++)
          fdst.x3f(X3E2,k,j,ie+1)=eflux_recv_[outer_x1][p++];
      }
      for(int k=ks; k<=ke+1; k++) {
        for(int j=js; j<=je; j++)
          wdst.x3f(k,j,ie+1)=eflux_recv_[outer_x1][p++];
      }
    }
    eflux_flag_[outer_x1][0][0] = false; // clear the flag
    if(pmb->neighbor[inner_x2][0][0].gid==-1) // physical boundary
      EFluxBoundary_[inner_x2](pmb,fdst,wdst);
    else {
      p=0;
      for(int k=ks; k<=ke; k++) {
        for(int i=is; i<=ie+1; i++)
          fdst.x1f(X1E3,k,js-1,i)=eflux_recv_[inner_x2][p++];
      }
      for(int k=ks; k<=ke; k++) {
        for(int i=is; i<=ie+1; i++)
          wdst.x1f(k,js-1,i)=eflux_recv_[inner_x2][p++];
      }
      for(int k=ks; k<=ke+1; k++) {
        for(int i=is; i<=ie; i++)
          fdst.x3f(X3E1,k,js-1,i)=eflux_recv_[inner_x2][p++];
      }
      for(int k=ks; k<=ke+1; k++) {
        for(int i=is; i<=ie; i++)
          wdst.x3f(k,js-1,i)=eflux_recv_[inner_x2][p++];
      }
    }
    eflux_flag_[inner_x2][0][0] = false; // clear the flag
    if(pmb->neighbor[outer_x2][0][0].gid==-1) // physical boundary
      EFluxBoundary_[outer_x2](pmb,fdst,wdst);
    else {
      p=0;
      for(int k=ks; k<=ke; k++) {
        for(int i=is; i<=ie+1; i++)
          fdst.x1f(X1E3,k,je+1,i)=eflux_recv_[outer_x2][p++];
      }
      for(int k=ks; k<=ke; k++) {
        for(int i=is; i<=ie+1; i++)
          wdst.x1f(k,je+1,i)=eflux_recv_[outer_x2][p++];
      }
      for(int k=ks; k<=ke+1; k++) {
        for(int i=is; i<=ie; i++)
          fdst.x3f(X3E1,k,je+1,i)=eflux_recv_[outer_x2][p++];
      }
      for(int k=ks; k<=ke+1; k++) {
        for(int i=is; i<=ie; i++)
          wdst.x3f(k,je+1,i)=eflux_recv_[outer_x2][p++];
      }
    }
    eflux_flag_[outer_x2][0][0] = false; // clear the flag
    if(pmb->neighbor[inner_x3][0][0].gid==-1) // physical boundary
      EFluxBoundary_[inner_x3](pmb,fdst,wdst);
    else {
      p=0;
      for(int j=js; j<=je; j++) {
        for(int i=is; i<=ie+1; i++)
          fdst.x1f(X1E2,ks-1,j,i)=eflux_recv_[inner_x3][p++];
      }
      for(int j=js; j<=je; j++) {
        for(int i=is; i<=ie+1; i++)
          wdst.x1f(ks-1,j,i)=eflux_recv_[inner_x3][p++];
      }
      for(int j=js; j<=je+1; j++) {
        for(int i=is; i<=ie; i++)
          fdst.x2f(X2E1,ks-1,j,i)=eflux_recv_[inner_x3][p++];
      }
      for(int j=js; j<=je+1; j++) {
        for(int i=is; i<=ie; i++)
          wdst.x2f(ks-1,j,i)=eflux_recv_[inner_x3][p++];
      }
    }
    eflux_flag_[inner_x3][0][0] = false; // clear the flag
    if(pmb->neighbor[outer_x3][0][0].gid==-1) // physical boundary
      EFluxBoundary_[outer_x3](pmb,fdst,wdst);
    else {
      p=0;
      for(int j=js; j<=je; j++) {
        for(int i=is; i<=ie+1; i++)
          fdst.x1f(X1E2,ke+1,j,i)=eflux_recv_[outer_x3][p++];
      }
      for(int j=js; j<=je; j++) {
        for(int i=is; i<=ie+1; i++)
          wdst.x1f(ke+1,j,i)=eflux_recv_[outer_x3][p++];
      }
      for(int j=js; j<=je+1; j++) {
        for(int i=is; i<=ie; i++)
          fdst.x2f(X2E1,ke+1,j,i)=eflux_recv_[outer_x3][p++];
      }
      for(int j=js; j<=je+1; j++) {
        for(int i=is; i<=ie; i++)
          wdst.x2f(ke+1,j,i)=eflux_recv_[outer_x3][p++];
      }
    }
    eflux_flag_[outer_x3][0][0] = false; // clear the flag
  }

  return true;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::WaitSendEFlux(void)
//  \brief wait until MPI_Isend completes for EFlux
void BoundaryValues::WaitSendEFlux(void) /***!!!something is wrong here!!!***/
{
#ifdef MPI_PARALLEL
  MeshBlock *pmb=pmy_mblock_;
  int ndir=4, dir;
  if(pmb->block_size.nx2==1) return; // 1D
  if(pmb->block_size.nx3>1) // 3D
    ndir=6;
  for(dir=0;dir<ndir;dir++) {
    if(pmb->neighbor[dir][0][0].rank!=-1 && pmb->neighbor[dir][0][0].rank!=myrank)
      MPI_Wait(&req_eflux_send_[dir][0][0],MPI_STATUS_IGNORE);
  }
#endif
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void InitBoundaryBuffer(int nx1, int nx2, int nx3)
//  \brief creates a list of the sizes and offsets of boundary buffers
void InitBoundaryBuffer(int nx1, int nx2, int nx3)
{
  int is, ie, js, je, ks, ke;

  is = NGHOST;
  ie = is + nx1 - 1;

  if (nx2 > 1) {
    js = NGHOST;
    je = js + nx2 - 1;
  } else {
    js = je = 0;
  }

  if (nx3 > 1) {
    ks = NGHOST;
    ke = ks + nx3 - 1;
  } else {
    ks = ke = 0;
  }
  fluid_send_se_[inner_x1][0]=is;
  fluid_send_se_[inner_x1][1]=is+NGHOST-1;
  fluid_send_se_[inner_x1][2]=js;
  fluid_send_se_[inner_x1][3]=je;
  fluid_send_se_[inner_x1][4]=ks;
  fluid_send_se_[inner_x1][5]=ke;

  fluid_send_se_[outer_x1][0]=ie-NGHOST+1;
  fluid_send_se_[outer_x1][1]=ie;
  fluid_send_se_[outer_x1][2]=js;
  fluid_send_se_[outer_x1][3]=je;
  fluid_send_se_[outer_x1][4]=ks;
  fluid_send_se_[outer_x1][5]=ke;

  fluid_send_se_[inner_x2][0]=0;
  fluid_send_se_[inner_x2][1]=ie+NGHOST;
  fluid_send_se_[inner_x2][2]=js;
  fluid_send_se_[inner_x2][3]=js+NGHOST-1;
  fluid_send_se_[inner_x2][4]=ks;
  fluid_send_se_[inner_x2][5]=ke;

  fluid_send_se_[outer_x2][0]=0;
  fluid_send_se_[outer_x2][1]=ie+NGHOST;
  fluid_send_se_[outer_x2][2]=je-NGHOST+1;
  fluid_send_se_[outer_x2][3]=je;
  fluid_send_se_[outer_x2][4]=ks;
  fluid_send_se_[outer_x2][5]=ke;

  fluid_send_se_[inner_x3][0]=0;
  fluid_send_se_[inner_x3][1]=ie+NGHOST;
  fluid_send_se_[inner_x3][2]=0;
  fluid_send_se_[inner_x3][3]=je+NGHOST;
  fluid_send_se_[inner_x3][4]=ks;
  fluid_send_se_[inner_x3][5]=ks+NGHOST-1;

  fluid_send_se_[outer_x3][0]=0;
  fluid_send_se_[outer_x3][1]=ie+NGHOST;
  fluid_send_se_[outer_x3][2]=0;
  fluid_send_se_[outer_x3][3]=je+NGHOST;
  fluid_send_se_[outer_x3][4]=ke-NGHOST+1;
  fluid_send_se_[outer_x3][5]=ke;


  fluid_recv_se_[inner_x1][0]=is-NGHOST;
  fluid_recv_se_[inner_x1][1]=is-1;
  fluid_recv_se_[inner_x1][2]=js;
  fluid_recv_se_[inner_x1][3]=je;
  fluid_recv_se_[inner_x1][4]=ks;
  fluid_recv_se_[inner_x1][5]=ke;

  fluid_recv_se_[outer_x1][0]=ie+1;
  fluid_recv_se_[outer_x1][1]=ie+NGHOST;
  fluid_recv_se_[outer_x1][2]=js;
  fluid_recv_se_[outer_x1][3]=je;
  fluid_recv_se_[outer_x1][4]=ks;
  fluid_recv_se_[outer_x1][5]=ke;

  fluid_recv_se_[inner_x2][0]=0;
  fluid_recv_se_[inner_x2][1]=ie+NGHOST;
  fluid_recv_se_[inner_x2][2]=js-NGHOST;
  fluid_recv_se_[inner_x2][3]=js-1;
  fluid_recv_se_[inner_x2][4]=ks;
  fluid_recv_se_[inner_x2][5]=ke;

  fluid_recv_se_[outer_x2][0]=0;
  fluid_recv_se_[outer_x2][1]=ie+NGHOST;
  fluid_recv_se_[outer_x2][2]=je+1;
  fluid_recv_se_[outer_x2][3]=je+NGHOST;
  fluid_recv_se_[outer_x2][4]=ks;
  fluid_recv_se_[outer_x2][5]=ke;

  fluid_recv_se_[inner_x3][0]=0;
  fluid_recv_se_[inner_x3][1]=ie+NGHOST;
  fluid_recv_se_[inner_x3][2]=0;
  fluid_recv_se_[inner_x3][3]=je+NGHOST;
  fluid_recv_se_[inner_x3][4]=ks-NGHOST;
  fluid_recv_se_[inner_x3][5]=ks-1;

  fluid_recv_se_[outer_x3][0]=0;
  fluid_recv_se_[outer_x3][1]=ie+NGHOST;
  fluid_recv_se_[outer_x3][2]=0;
  fluid_recv_se_[outer_x3][3]=je+NGHOST;
  fluid_recv_se_[outer_x3][4]=ke+1;
  fluid_recv_se_[outer_x3][5]=ke+NGHOST;

  fluid_bufsize_[inner_x1]=NGHOST*nx2*nx3*NFLUID;
  fluid_bufsize_[outer_x1]=NGHOST*nx2*nx3*NFLUID;
  fluid_bufsize_[inner_x2]=(nx1+2*NGHOST)*NGHOST*nx3*NFLUID;
  fluid_bufsize_[outer_x2]=(nx1+2*NGHOST)*NGHOST*nx3*NFLUID;
  fluid_bufsize_[inner_x3]=(nx1+2*NGHOST)*(nx2+2*NGHOST)*NGHOST*NFLUID;
  fluid_bufsize_[outer_x3]=(nx1+2*NGHOST)*(nx2+2*NGHOST)*NGHOST*NFLUID;

  if (MAGNETIC_FIELDS_ENABLED) {
    field_send_se_[inner_x1][x1face][0]=is+1;
    field_send_se_[inner_x1][x1face][1]=is+NGHOST;
    field_send_se_[inner_x1][x1face][2]=js;
    field_send_se_[inner_x1][x1face][3]=je;
    field_send_se_[inner_x1][x1face][4]=ks;
    field_send_se_[inner_x1][x1face][5]=ke;

    field_send_se_[inner_x1][x2face][0]=is;
    field_send_se_[inner_x1][x2face][1]=is+NGHOST-1;
    field_send_se_[inner_x1][x2face][2]=js;
    field_send_se_[inner_x1][x2face][3]=je+1;
    field_send_se_[inner_x1][x2face][4]=ks;
    field_send_se_[inner_x1][x2face][5]=ke;

    field_send_se_[inner_x1][x3face][0]=is;
    field_send_se_[inner_x1][x3face][1]=is+NGHOST-1;
    field_send_se_[inner_x1][x3face][2]=js;
    field_send_se_[inner_x1][x3face][3]=je;
    field_send_se_[inner_x1][x3face][4]=ks;
    field_send_se_[inner_x1][x3face][5]=ke+1;

    field_send_se_[outer_x1][x1face][0]=ie-NGHOST+1;
    field_send_se_[outer_x1][x1face][1]=ie;
    field_send_se_[outer_x1][x1face][2]=js;
    field_send_se_[outer_x1][x1face][3]=je;
    field_send_se_[outer_x1][x1face][4]=ks;
    field_send_se_[outer_x1][x1face][5]=ke;

    field_send_se_[outer_x1][x2face][0]=ie-NGHOST+1;
    field_send_se_[outer_x1][x2face][1]=ie;
    field_send_se_[outer_x1][x2face][2]=js;
    field_send_se_[outer_x1][x2face][3]=je+1;
    field_send_se_[outer_x1][x2face][4]=ks;
    field_send_se_[outer_x1][x2face][5]=ke;

    field_send_se_[outer_x1][x3face][0]=ie-NGHOST+1;
    field_send_se_[outer_x1][x3face][1]=ie;
    field_send_se_[outer_x1][x3face][2]=js;
    field_send_se_[outer_x1][x3face][3]=je;
    field_send_se_[outer_x1][x3face][4]=ks;
    field_send_se_[outer_x1][x3face][5]=ke+1;

    field_send_se_[inner_x2][x1face][0]=0;
    field_send_se_[inner_x2][x1face][1]=ie+NGHOST+1;
    field_send_se_[inner_x2][x1face][2]=js;
    field_send_se_[inner_x2][x1face][3]=js+NGHOST-1;
    field_send_se_[inner_x2][x1face][4]=ks;
    field_send_se_[inner_x2][x1face][5]=ke;

    field_send_se_[inner_x2][x2face][0]=0;
    field_send_se_[inner_x2][x2face][1]=ie+NGHOST;
    field_send_se_[inner_x2][x2face][2]=js+1;
    field_send_se_[inner_x2][x2face][3]=js+NGHOST;
    field_send_se_[inner_x2][x2face][4]=ks;
    field_send_se_[inner_x2][x2face][5]=ke;

    field_send_se_[inner_x2][x3face][0]=0;
    field_send_se_[inner_x2][x3face][1]=ie+NGHOST;
    field_send_se_[inner_x2][x3face][2]=js;
    field_send_se_[inner_x2][x3face][3]=js+NGHOST-1;
    field_send_se_[inner_x2][x3face][4]=ks;
    field_send_se_[inner_x2][x3face][5]=ke+1;

    field_send_se_[outer_x2][x1face][0]=0;
    field_send_se_[outer_x2][x1face][1]=ie+NGHOST+1;
    field_send_se_[outer_x2][x1face][2]=je-NGHOST+1;
    field_send_se_[outer_x2][x1face][3]=je;
    field_send_se_[outer_x2][x1face][4]=ks;
    field_send_se_[outer_x2][x1face][5]=ke;

    field_send_se_[outer_x2][x2face][0]=0;
    field_send_se_[outer_x2][x2face][1]=ie+NGHOST;
    field_send_se_[outer_x2][x2face][2]=je-NGHOST+1;
    field_send_se_[outer_x2][x2face][3]=je;
    field_send_se_[outer_x2][x2face][4]=ks;
    field_send_se_[outer_x2][x2face][5]=ke;

    field_send_se_[outer_x2][x3face][0]=0;
    field_send_se_[outer_x2][x3face][1]=ie+NGHOST;
    field_send_se_[outer_x2][x3face][2]=je-NGHOST+1;
    field_send_se_[outer_x2][x3face][3]=je;
    field_send_se_[outer_x2][x3face][4]=ks;
    field_send_se_[outer_x2][x3face][5]=ke+1;

    field_send_se_[inner_x3][x1face][0]=0;
    field_send_se_[inner_x3][x1face][1]=ie+NGHOST+1;
    field_send_se_[inner_x3][x1face][2]=0;
    field_send_se_[inner_x3][x1face][3]=je+NGHOST;
    field_send_se_[inner_x3][x1face][4]=ks;
    field_send_se_[inner_x3][x1face][5]=ks+NGHOST-1;

    field_send_se_[inner_x3][x2face][0]=0;
    field_send_se_[inner_x3][x2face][1]=ie+NGHOST;
    field_send_se_[inner_x3][x2face][2]=0;
    field_send_se_[inner_x3][x2face][3]=je+NGHOST+1;
    field_send_se_[inner_x3][x2face][4]=ks;
    field_send_se_[inner_x3][x2face][5]=ks+NGHOST-1;

    field_send_se_[inner_x3][x3face][0]=0;
    field_send_se_[inner_x3][x3face][1]=ie+NGHOST;
    field_send_se_[inner_x3][x3face][2]=0;
    field_send_se_[inner_x3][x3face][3]=je+NGHOST;
    field_send_se_[inner_x3][x3face][4]=ks+1;
    field_send_se_[inner_x3][x3face][5]=ks+NGHOST;

    field_send_se_[outer_x3][x1face][0]=0;
    field_send_se_[outer_x3][x1face][1]=ie+NGHOST+1;
    field_send_se_[outer_x3][x1face][2]=0;
    field_send_se_[outer_x3][x1face][3]=je+NGHOST;
    field_send_se_[outer_x3][x1face][4]=ke-NGHOST+1;
    field_send_se_[outer_x3][x1face][5]=ke;

    field_send_se_[outer_x3][x2face][0]=0;
    field_send_se_[outer_x3][x2face][1]=ie+NGHOST;
    field_send_se_[outer_x3][x2face][2]=0;
    field_send_se_[outer_x3][x2face][3]=je+NGHOST+1;
    field_send_se_[outer_x3][x2face][4]=ke-NGHOST+1;
    field_send_se_[outer_x3][x2face][5]=ke;

    field_send_se_[outer_x3][x3face][0]=0;
    field_send_se_[outer_x3][x3face][1]=ie+NGHOST;
    field_send_se_[outer_x3][x3face][2]=0;
    field_send_se_[outer_x3][x3face][3]=je+NGHOST;
    field_send_se_[outer_x3][x3face][4]=ke-NGHOST+1;
    field_send_se_[outer_x3][x3face][5]=ke;

    field_recv_se_[inner_x1][x1face][0]=is-NGHOST;
    field_recv_se_[inner_x1][x1face][1]=is-1;
    field_recv_se_[inner_x1][x1face][2]=js;
    field_recv_se_[inner_x1][x1face][3]=je;
    field_recv_se_[inner_x1][x1face][4]=ks;
    field_recv_se_[inner_x1][x1face][5]=ke;

    field_recv_se_[inner_x1][x2face][0]=is-NGHOST;
    field_recv_se_[inner_x1][x2face][1]=is-1;
    field_recv_se_[inner_x1][x2face][2]=js;
    field_recv_se_[inner_x1][x2face][3]=je+1;
    field_recv_se_[inner_x1][x2face][4]=ks;
    field_recv_se_[inner_x1][x2face][5]=ke;

    field_recv_se_[inner_x1][x3face][0]=is-NGHOST;
    field_recv_se_[inner_x1][x3face][1]=is-1;
    field_recv_se_[inner_x1][x3face][2]=js;
    field_recv_se_[inner_x1][x3face][3]=je;
    field_recv_se_[inner_x1][x3face][4]=ks;
    field_recv_se_[inner_x1][x3face][5]=ke+1;

    field_recv_se_[outer_x1][x1face][0]=ie+2;
    field_recv_se_[outer_x1][x1face][1]=ie+NGHOST+1;
    field_recv_se_[outer_x1][x1face][2]=js;
    field_recv_se_[outer_x1][x1face][3]=je;
    field_recv_se_[outer_x1][x1face][4]=ks;
    field_recv_se_[outer_x1][x1face][5]=ke;

    field_recv_se_[outer_x1][x2face][0]=ie+1;
    field_recv_se_[outer_x1][x2face][1]=ie+NGHOST;
    field_recv_se_[outer_x1][x2face][2]=js;
    field_recv_se_[outer_x1][x2face][3]=je+1;
    field_recv_se_[outer_x1][x2face][4]=ks;
    field_recv_se_[outer_x1][x2face][5]=ke;

    field_recv_se_[outer_x1][x3face][0]=ie+1;
    field_recv_se_[outer_x1][x3face][1]=ie+NGHOST;
    field_recv_se_[outer_x1][x3face][2]=js;
    field_recv_se_[outer_x1][x3face][3]=je;
    field_recv_se_[outer_x1][x3face][4]=ks;
    field_recv_se_[outer_x1][x3face][5]=ke+1;

    field_recv_se_[inner_x2][x1face][0]=0;
    field_recv_se_[inner_x2][x1face][1]=ie+NGHOST+1;
    field_recv_se_[inner_x2][x1face][2]=js-NGHOST;
    field_recv_se_[inner_x2][x1face][3]=js-1;
    field_recv_se_[inner_x2][x1face][4]=ks;
    field_recv_se_[inner_x2][x1face][5]=ke;

    field_recv_se_[inner_x2][x2face][0]=0;
    field_recv_se_[inner_x2][x2face][1]=ie+NGHOST;
    field_recv_se_[inner_x2][x2face][2]=js-NGHOST;
    field_recv_se_[inner_x2][x2face][3]=js-1;
    field_recv_se_[inner_x2][x2face][4]=ks;
    field_recv_se_[inner_x2][x2face][5]=ke;

    field_recv_se_[inner_x2][x3face][0]=0;
    field_recv_se_[inner_x2][x3face][1]=ie+NGHOST;
    field_recv_se_[inner_x2][x3face][2]=js-NGHOST;
    field_recv_se_[inner_x2][x3face][3]=js-1;
    field_recv_se_[inner_x2][x3face][4]=ks;
    field_recv_se_[inner_x2][x3face][5]=ke+1;

    field_recv_se_[outer_x2][x1face][0]=0;
    field_recv_se_[outer_x2][x1face][1]=ie+NGHOST+1;
    field_recv_se_[outer_x2][x1face][2]=je+1;
    field_recv_se_[outer_x2][x1face][3]=je+NGHOST;
    field_recv_se_[outer_x2][x1face][4]=ks;
    field_recv_se_[outer_x2][x1face][5]=ke;

    field_recv_se_[outer_x2][x2face][0]=0;
    field_recv_se_[outer_x2][x2face][1]=ie+NGHOST;
    field_recv_se_[outer_x2][x2face][2]=je+2;
    field_recv_se_[outer_x2][x2face][3]=je+NGHOST+1;
    field_recv_se_[outer_x2][x2face][4]=ks;
    field_recv_se_[outer_x2][x2face][5]=ke;

    field_recv_se_[outer_x2][x3face][0]=0;
    field_recv_se_[outer_x2][x3face][1]=ie+NGHOST;
    field_recv_se_[outer_x2][x3face][2]=je+1;
    field_recv_se_[outer_x2][x3face][3]=je+NGHOST;
    field_recv_se_[outer_x2][x3face][4]=ks;
    field_recv_se_[outer_x2][x3face][5]=ke+1;

    field_recv_se_[inner_x3][x1face][0]=0;
    field_recv_se_[inner_x3][x1face][1]=ie+NGHOST+1;
    field_recv_se_[inner_x3][x1face][2]=0;
    field_recv_se_[inner_x3][x1face][3]=je+NGHOST;
    field_recv_se_[inner_x3][x1face][4]=ks-NGHOST;
    field_recv_se_[inner_x3][x1face][5]=ks-1;

    field_recv_se_[inner_x3][x2face][0]=0;
    field_recv_se_[inner_x3][x2face][1]=ie+NGHOST;
    field_recv_se_[inner_x3][x2face][2]=0;
    field_recv_se_[inner_x3][x2face][3]=je+NGHOST+1;
    field_recv_se_[inner_x3][x2face][4]=ks-NGHOST;
    field_recv_se_[inner_x3][x2face][5]=ks-1;

    field_recv_se_[inner_x3][x3face][0]=0;
    field_recv_se_[inner_x3][x3face][1]=ie+NGHOST;
    field_recv_se_[inner_x3][x3face][2]=0;
    field_recv_se_[inner_x3][x3face][3]=je+NGHOST;
    field_recv_se_[inner_x3][x3face][4]=ks-NGHOST;
    field_recv_se_[inner_x3][x3face][5]=ks-1;

    field_recv_se_[outer_x3][x1face][0]=0;
    field_recv_se_[outer_x3][x1face][1]=ie+NGHOST+1;
    field_recv_se_[outer_x3][x1face][2]=0;
    field_recv_se_[outer_x3][x1face][3]=je+NGHOST;
    field_recv_se_[outer_x3][x1face][4]=ke+1;
    field_recv_se_[outer_x3][x1face][5]=ke+NGHOST;

    field_recv_se_[outer_x3][x2face][0]=0;
    field_recv_se_[outer_x3][x2face][1]=ie+NGHOST;
    field_recv_se_[outer_x3][x2face][2]=0;
    field_recv_se_[outer_x3][x2face][3]=je+NGHOST+1;
    field_recv_se_[outer_x3][x2face][4]=ke+1;
    field_recv_se_[outer_x3][x2face][5]=ke+NGHOST;

    field_recv_se_[outer_x3][x3face][0]=0;
    field_recv_se_[outer_x3][x3face][1]=ie+NGHOST;
    field_recv_se_[outer_x3][x3face][2]=0;
    field_recv_se_[outer_x3][x3face][3]=je+NGHOST;
    field_recv_se_[outer_x3][x3face][4]=ke+2;
    field_recv_se_[outer_x3][x3face][5]=ke+NGHOST+1;

    field_bufsize_[inner_x1]=field_bufsize_[outer_x1]
                            =NGHOST*(nx2*nx3+(nx2+1)*nx3+nx2*(nx3+1));
    field_bufsize_[inner_x2]=field_bufsize_[outer_x2]=NGHOST*((nx1+2*NGHOST)*nx3
                            +(nx1+2*NGHOST+1)*nx3+(nx1+2*NGHOST)*(nx3+1));
    field_bufsize_[inner_x3]=field_bufsize_[outer_x3]
                  =NGHOST*((nx1+2*NGHOST+1)*(nx2+2*NGHOST)
                  +(nx1+2*NGHOST)*(nx2+2*NGHOST+1)+(nx1+2*NGHOST)*(nx2+2*NGHOST));
    if(nx2==1) return; // 1D
    else {
      if(nx3==1) { // 2D
        eflux_bufsize_[inner_x1]=eflux_bufsize_[outer_x1]=(nx2+1)*2;
        eflux_bufsize_[inner_x2]=eflux_bufsize_[outer_x2]=(nx1+1)*2;
      }
      else { // 3D
        eflux_bufsize_[inner_x1]=eflux_bufsize_[outer_x1]=(nx2+1)*nx3*2+nx2*(nx3+1)*2;
        eflux_bufsize_[inner_x2]=eflux_bufsize_[outer_x2]=(nx1+1)*nx3*2+nx1*(nx3+1)*2;
        eflux_bufsize_[inner_x3]=eflux_bufsize_[outer_x3]=(nx1+1)*nx2*2+nx1*(nx2+1)*2;
      }
    }
  }
  return;
}

