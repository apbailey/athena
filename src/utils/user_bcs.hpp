#ifndef UTILS_USER_BCS_HPP_
#define UTILS_USER_BCS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file user_bcs.hpp
//! \brief Runtime-selectable user boundary conditions:
//!   - "diode" : outflow + clamp v_perp to 0 on attempted inflow
//!   - "vacuum": small (rho, P) ghosts with v = 0

#include "../athena.hpp"
#include "../bvals/bvals_interfaces.hpp"

class Mesh;
class MeshBlock;
class Coordinates;
class ParameterInput;

// --- Diode BCs (12 face directions split across the 6 names below) ---
void DiodeInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int is, int ie, int js, int je, int ks, int ke, int ngh);
void DiodeOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int is, int ie, int js, int je, int ks, int ke, int ngh);
void DiodeInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int is, int ie, int js, int je, int ks, int ke, int ngh);
void DiodeOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int is, int ie, int js, int je, int ks, int ke, int ngh);
void DiodeInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int is, int ie, int js, int je, int ks, int ke, int ngh);
void DiodeOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int is, int ie, int js, int je, int ks, int ke, int ngh);

// --- Vacuum BCs ---
void VacuumInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt,
                   int is, int ie, int js, int je, int ks, int ke, int ngh);
void VacuumOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt,
                   int is, int ie, int js, int je, int ks, int ke, int ngh);
void VacuumInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt,
                   int is, int ie, int js, int je, int ks, int ke, int ngh);
void VacuumOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt,
                   int is, int ie, int js, int je, int ks, int ke, int ngh);
void VacuumInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt,
                   int is, int ie, int js, int je, int ks, int ke, int ngh);
void VacuumOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt,
                   int is, int ie, int js, int je, int ks, int ke, int ngh);

// Override vacuum floors.  Defaults match the EOS floors
// (<hydro>/dfloor, <hydro>/pfloor) when ResolveUserBC is invoked.
void SetVacuumFloors(Real rho_floor, Real press_floor);

// Read <problem>/<param> string and return the corresponding BC function
// pointer for `face`.  Names:
//   "diode"  -> Diode*X*
//   "vacuum" -> Vacuum*X*
//   ""       -> nullptr (no user BC requested; caller should skip enrolling)
//   else     -> ATHENA_ERROR
//
// Also reads <problem>/vacuum_rho_floor and <problem>/vacuum_press_floor
// (defaults 1e-6) so any later vacuum call uses the configured values.
//
// Mesh::EnrollUserBoundaryFunction is private, so the caller (which lives
// in Mesh::InitUserMeshData and thus has access) must do the enrollment:
//
//   BValFunc fn = ResolveUserBC(pin, "ix3_user_bc", BoundaryFace::inner_x3);
//   if (fn) EnrollUserBoundaryFunction(BoundaryFace::inner_x3, fn);
BValFunc ResolveUserBC(ParameterInput *pin, const char *param,
                       BoundaryFace face);

#endif  // UTILS_USER_BCS_HPP_
