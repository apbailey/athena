//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file disk_planet.cpp
//! \brief Global star-centered disk (spherical_polar, rotating frame) with an embedded
//! planet — the global equivalent of the local cartesian Hill setup
//! (bondi_accretion_cartesian.cpp), built on disk.cpp.
//!
//! Frame/geometry: the star sits at the origin of a spherical_polar grid; the whole
//! domain corotates with the planet at Omega0 = sqrt((GM+GM_p)/r0^3) about the z-axis.
//! The planet is FIXED at (r, theta, phi) = (r0, pi/2, 0). Locally, refined cells near
//! the planet (dr ~ r0 dtheta ~ r0 dphi) form a quasi-cartesian patch.
//!
//! Source terms:
//!   built-in  : star point-mass gravity  (<problem> GM, pointmass.cpp)
//!               Coriolis + centrifugal   (<orbital_advection> Omega0, set here;
//!                                         rotating_system_srcterms.cpp)
//!   this pgen : planet gravity  a = -GM_p d/(|d|^2+eps^2)^(3/2), d = x - x_p,
//!               with the same sin^2 ramp over tinj as the local pgen;
//!               indirect term  a = -GM_p/r0^2 xhat  (acceleration of the star-centered,
//!               non-inertial origin; with it the expansion about the planet reproduces
//!               the Hill terms 3 Omega^2 x_loc, -Omega^2 z_loc exactly);
//!               vacuum sink    every stage, cells with |x-x_p| < r_sink are reset to
//!               (sink_density, zero momentum in the rotating frame); for adiabatic the
//!               energy is set to sink_density*cs0^2/(gamma-1). Identical prescription
//!               to bondi_accretion_cartesian.cpp.
//!
//! Units and the map to the local Hill runs: with GM = r0 = 1 the planet's orbit has
//! Omega_K = 1 and the local unit length H0 = cs/Omega_K = iso_sound_speed (= the disk
//! aspect ratio h0). The dimensionless inputs carry over from the local runs unchanged:
//!   qthermal              GM_p = qthermal * cs^3/Omega_K   (same number as local)
//!   r_sink_h, epsilon_h   sink radius / softening in units of H0 (local code units)
//! so e.g. the local q20 run (qthermal=19.58, r_sink=0.1) maps to qthermal=19.58,
//! r_sink_h=0.1 here. R_H = r0 (GM_p/3GM)^(1/3) = (qthermal/3)^(1/3) h0 r0, matching the
//! local R_H in H0 units.
//!
//! Boundary conditions: disk.cpp-style user BCs pin all ghost zones to the rotating-frame
//! equilibrium background. theta faces may instead use reflecting (e.g. ox2 reflecting at
//! the midplane to evolve only z>0, as in the local runs).
//!
//! Optional AMR: with <mesh> refinement=adaptive, <problem> amr_r_hill (in Hill radii)
//! refines blocks within that distance of the planet up to <mesh> numlevel.

// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>      // sqrt, sin, cos, pow
#include <cstring>    // strcmp()
#include <iostream>   // endl, cout
#include <limits>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../parameter_input.hpp"
#include "../utils/user_bcs.hpp"

namespace {
void GetCylCoord(Coordinates *pco, Real &rad, Real &phi, Real &z, int i, int j, int k);
Real DenProfileCyl(const Real rad, const Real phi, const Real z);
Real PoverR(const Real rad, const Real phi, const Real z);
Real VelProfileCyl(const Real rad, const Real phi, const Real z);
Real PlanetPotential(Real x, Real y, Real z, Real time);
Real LocalIsoCs(Real x1, Real x2, Real x3);
// disk background parameters (as in disk.cpp)
Real gm0, r0, rho0, dslope, p0_over_r0, pslope, gamma_gas;
Real dfloor;
Real Omega0;
// planet parameters
Real gmp;            // GM_planet = qthermal * cs0^3 / Omega_K
Real qthermal;       // dimensionless thermal mass (input)
Real eps_soft;       // gravitational softening (code units; input in H0 units)
Real tinj;           // planet-gravity ramp time (code units)
Real r_sink;         // sink radius (code units; input in H0 units)
Real sink_density;   // density inside the sink
Real cs0;            // sound speed at r0
Real h0;             // H0 = cs0/Omega_K (local unit length in code units)
Real r_hill;         // Hill radius (code units)
bool indirect_term;  // include -GM_p/r0^2 xhat
Real amr_dist;       // refine blocks within this distance of the planet (code units)
} // namespace

void PlanetSource(MeshBlock *pmb, const Real time, const Real dt,
                  const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                  AthenaArray<Real> &cons_scalar);
int RefinementCondition(MeshBlock *pmb);

// User-defined boundary conditions (background pin, from disk.cpp)
void DiskInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
// Azimuthal background pin, for runs on a PATCH of the disk rather than the
// full 2pi annulus.  With x3 spanning the whole circle these are unused and
// <mesh>/ix3_bc,ox3_bc stay 'periodic'; set both to 'user' to pin instead.
void DiskInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief read parameters, derive planet/frame quantities, enroll everything
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in disk_planet::InitUserMeshData" << std::endl
        << "Coordinate system must be 'spherical_polar'; got '"
        << COORDINATE_SYSTEM << "'." << std::endl;
    ATHENA_ERROR(msg);
  }
  if (mesh_size.nx2 == 1 || mesh_size.nx3 == 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in disk_planet::InitUserMeshData" << std::endl
        << "This pgen is 3D-only; got nx2=" << mesh_size.nx2
        << " nx3=" << mesh_size.nx3 << "." << std::endl;
    ATHENA_ERROR(msg);
  }

  // ---- star + disk background (disk.cpp parameters) ----
  gm0 = pin->GetOrAddReal("problem", "GM", 1.0);
  r0 = pin->GetOrAddReal("problem", "r0", 1.0);
  rho0 = pin->GetReal("problem", "rho0");
  dslope = pin->GetOrAddReal("problem", "dslope", 0.0);
  if (NON_BAROTROPIC_EOS) {
    p0_over_r0 = pin->GetOrAddReal("problem", "p0_over_r0", 0.0025);
    pslope = pin->GetOrAddReal("problem", "pslope", 0.0);
    gamma_gas = pin->GetReal("hydro", "gamma");
  } else {
    p0_over_r0 = SQR(pin->GetReal("hydro", "iso_sound_speed"));
    // locally isothermal: P = rho c_s^2(R) with c_s^2 = p0_over_r0 (R/r0)^pslope
    // (R = cylindrical radius; vertically isothermal). pslope = 0 recovers the
    // globally isothermal disk. Requires the locally isothermal EOS extension.
    pslope = pin->GetOrAddReal("problem", "pslope", 0.0);
  }
  Real float_min = std::numeric_limits<float>::min();
  dfloor = pin->GetOrAddReal("hydro", "dfloor", (1024*(float_min)));

  // ---- planet: dimensionless inputs, local-run conventions ----
  // Guard against the local-pgen parameter names being used with global units.
  if (pin->DoesParameterExist("problem", "r_sink") ||
      pin->DoesParameterExist("problem", "epsilon")) {
    std::stringstream msg;
    msg << "### FATAL ERROR in disk_planet::InitUserMeshData" << std::endl
        << "Use r_sink_h / epsilon_h (in units of H0 = cs/Omega_K, i.e. the local"
        << std::endl
        << "Hill-run code units), not r_sink / epsilon." << std::endl;
    ATHENA_ERROR(msg);
  }
  qthermal = pin->GetOrAddReal("problem", "qthermal", 0.0);
  Real epsilon_h = pin->GetOrAddReal("problem", "epsilon_h", 0.0);
  Real r_sink_h = pin->GetOrAddReal("problem", "r_sink_h", 0.0);
  sink_density = pin->GetOrAddReal("problem", "sink_density", 1.0e-10);
  tinj = pin->GetOrAddReal("problem", "tinj", 0.0);
  indirect_term = pin->GetOrAddBoolean("problem", "indirect_term", true);

  cs0 = std::sqrt(p0_over_r0);                       // sound speed at r0
  Real omega_k = std::sqrt(gm0/(r0*r0*r0));          // star-only Keplerian
  h0 = cs0/omega_k;
  gmp = qthermal*cs0*cs0*cs0/omega_k;
  eps_soft = epsilon_h*h0;
  r_sink = r_sink_h*h0;
  r_hill = r0*std::cbrt(gmp/(3.0*gm0));

  // Frame rotation: planet on a circular orbit of the two-body problem. Set Omega0 in
  // the input so the built-in rotating-frame source terms and the BCs use the same
  // value; warn if the input supplied something else.
  Real omega_p = std::sqrt((gm0 + gmp)/(r0*r0*r0));
  Real omega_in = pin->GetOrAddReal("orbital_advection", "Omega0", 0.0);
  if (omega_in != 0.0 && std::abs(omega_in - omega_p) > 1.0e-12*omega_p
      && Globals::my_rank == 0) {
    std::cout << "### WARNING in disk_planet: <orbital_advection>/Omega0 = " << omega_in
              << " overridden by the corotation value " << omega_p << std::endl;
  }
  pin->SetReal("orbital_advection", "Omega0", omega_p);
  Omega0 = omega_p;

  // Warn if EOS floors would override the vacuum-sink clamp (as in the local pgen).
  if (r_sink > 0.0 && Globals::my_rank == 0) {
    Real eos_dfloor = pin->GetOrAddReal(
        "hydro", "dfloor", std::sqrt(1024.0*std::numeric_limits<Real>::min()));
    if (eos_dfloor > sink_density) {
      std::cout << "### WARNING in disk_planet: <hydro>/dfloor = " << eos_dfloor
                << " > problem/sink_density = " << sink_density << std::endl
                << "  the vacuum sink clamp is not honored;"
                << " lower <hydro>/dfloor below sink_density." << std::endl;
    }
    if (NON_BAROTROPIC_EOS) {
      Real eos_pfloor = pin->GetOrAddReal(
          "hydro", "pfloor", std::sqrt(1024.0*std::numeric_limits<Real>::min()));
      Real target_press = sink_density*p0_over_r0;
      if (eos_pfloor > target_press) {
        std::cout << "### WARNING in disk_planet: <hydro>/pfloor = " << eos_pfloor
                  << " > sink pressure target = " << target_press << std::endl
                  << "  the vacuum sink pressure is not honored;"
                  << " lower <hydro>/pfloor." << std::endl;
      }
    }
  }

  if (Globals::my_rank == 0) {
    std::cout << "disk_planet: qthermal=" << qthermal << " GM_p=" << gmp
              << " (GM=" << gm0 << ", r0=" << r0 << ")" << std::endl
              << "  h0=H0/r0=" << h0/r0 << "  R_H=" << r_hill << " (=" << r_hill/h0
              << " H0)  Omega0=" << Omega0 << std::endl
              << "  r_sink=" << r_sink << " (=" << r_sink_h << " H0, "
              << (r_hill > 0.0 ? r_sink/r_hill : 0.0) << " R_H)  eps=" << eps_soft
              << "  tinj=" << tinj << "  sink_density=" << sink_density << std::endl
              << "  indirect_term=" << (indirect_term ? "true" : "false")
              << "  EOS=" << (NON_BAROTROPIC_EOS ? "adiabatic" : "isothermal")
              << std::endl;
  }

  if (!NON_BAROTROPIC_EOS && pslope != 0.0) {
    EquationOfState::EnrollIsoSoundSpeed(LocalIsoCs);
    if (Globals::my_rank == 0) {
      std::cout << "disk_planet: locally isothermal, c_s = " << cs0
                << " (R/r0)^" << 0.5*pslope << std::endl;
    }
  }

  EnrollUserExplicitSourceFunction(PlanetSource);

  // Optional distance-based AMR around the planet.
  Real amr_r_hill = pin->GetOrAddReal("problem", "amr_r_hill", 0.0);
  amr_dist = amr_r_hill*r_hill;
  if (adaptive) {
    EnrollUserRefinementCondition(RefinementCondition);
  }

  // Background-pin boundary conditions on faces flagged 'user'.
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user"))
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, DiskInnerX1);
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user"))
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, DiskOuterX1);
  if (mesh_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("user"))
    EnrollUserBoundaryFunction(BoundaryFace::inner_x2, DiskInnerX2);
  // outer_x2 is the MIDPLANE (x2max = pi/2 in these half-domain runs).  The
  // physical condition there is 'reflecting': for a symmetric disk the
  // midplane is a true symmetry plane.  Setting <problem>/ox2_user_bc=diode
  // instead lets gas leave through the midplane but never re-enter, which is
  // deliberately UNPHYSICAL -- it removes the material the mirror hemisphere
  // would have supplied.  It exists only as a boundary-sensitivity test: it
  // asks how much of the polar funnel survives without the midplane bounce.
  // The local Hill-box runs have the same switch (ix3_user_bc=diode there).
  // Unset, this falls back to the background pin exactly as before.
  if (mesh_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("user")) {
    BValFunc fn = ResolveUserBC(pin, "ox2_user_bc", BoundaryFace::outer_x2);
    EnrollUserBoundaryFunction(BoundaryFace::outer_x2, fn ? fn : DiskOuterX2);
  }
  if (mesh_bcs[BoundaryFace::inner_x3] == GetBoundaryFlag("user"))
    EnrollUserBoundaryFunction(BoundaryFace::inner_x3, DiskInnerX3);
  if (mesh_bcs[BoundaryFace::outer_x3] == GetBoundaryFlag("user"))
    EnrollUserBoundaryFunction(BoundaryFace::outer_x3, DiskOuterX3);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief disk.cpp initial condition (rotating-frame velocities via Omega0)
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real den, vel;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        GetCylCoord(pcoord, rad, phi, z, i, j, k);
        den = DenProfileCyl(rad, phi, z);
        vel = VelProfileCyl(rad, phi, z);
        phydro->u(IDN,k,j,i) = den;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = den*vel;
        if (NON_BAROTROPIC_EOS) {
          Real p_over_r = PoverR(rad, phi, z);
          phydro->u(IEN,k,j,i) = p_over_r*phydro->u(IDN,k,j,i)/(gamma_gas - 1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }
      }
    }
  }
  return;
}

namespace {
//----------------------------------------------------------------------------------------
//! transform to cylindrical coordinates (spherical_polar host grid)

void GetCylCoord(Coordinates *pco, Real &rad, Real &phi, Real &z, int i, int j, int k) {
  rad = std::abs(pco->x1v(i)*std::sin(pco->x2v(j)));
  phi = pco->x3v(k);
  z = pco->x1v(i)*std::cos(pco->x2v(j));
  return;
}

//----------------------------------------------------------------------------------------
//! background density (vertical hydrostatic equilibrium about the star)

Real DenProfileCyl(const Real rad, const Real phi, const Real z) {
  Real den;
  Real p_over_r = PoverR(rad, phi, z);
  Real denmid = rho0*std::pow(rad/r0, dslope);
  Real dentem = denmid*std::exp(gm0/p_over_r*(1./std::sqrt(SQR(rad)+SQR(z))-1./rad));
  den = dentem;
  return std::max(den, dfloor);
}

//----------------------------------------------------------------------------------------
//! P/rho profile

Real PoverR(const Real rad, const Real phi, const Real z) {
  return p0_over_r0*std::pow(rad/r0, pslope);
}

//----------------------------------------------------------------------------------------
//! background rotational velocity, rotating frame (subtracts rad*Omega0)

Real VelProfileCyl(const Real rad, const Real phi, const Real z) {
  Real p_over_r = PoverR(rad, phi, z);
  Real vel = (dslope+pslope)*p_over_r/(gm0/rad) + (1.0+pslope)
             - pslope*rad/std::sqrt(rad*rad+z*z);
  vel = std::sqrt(gm0/rad)*std::sqrt(vel) - rad*Omega0;
  return vel;
}

//----------------------------------------------------------------------------------------
//! locally isothermal sound speed c_s(R_cyl) on the spherical_polar grid
//! (enrolled with the EOS when pslope != 0; vertically isothermal)

Real LocalIsoCs(Real x1, Real x2, Real x3) {
  Real rad = std::abs(x1*std::sin(x2));
  return cs0*std::pow(rad/r0, 0.5*pslope);
}

//----------------------------------------------------------------------------------------
//! planet + indirect potential at cartesian (x, y, z), with the tinj ramp.
//! Used for the conservative flux-form energy source (adiabatic only).

Real PlanetPotential(Real x, Real y, Real z, Real time) {
  Real dx = x - r0;
  Real d = std::sqrt(dx*dx + y*y + z*z + eps_soft*eps_soft);
  Real pot = -gmp/d;
  if (indirect_term) pot += gmp*x/(r0*r0);
  if (time < tinj)
    pot *= SQR(std::sin(PI*time/tinj/2.0));
  return pot;
}
} // namespace

//----------------------------------------------------------------------------------------
//! PlanetSource: planet gravity + indirect term + vacuum sink.
//! Star gravity and Coriolis/centrifugal are handled by the built-in source terms.

void PlanetSource(MeshBlock *pmb, const Real time, const Real dt,
                  const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                  AthenaArray<Real> &cons_scalar) {
  const Real r_sink2 = r_sink*r_sink;
  Real ramp = 1.0;
  if (time < tinj) ramp = SQR(std::sin(PI*time/tinj/2.0));
  const Real a_ind = indirect_term ? -ramp*gmp/(r0*r0) : 0.0;  // along -xhat

  AthenaArray<Real> x1area, x2aream, x2areap, x3aream, x3areap, vol;
  if (NON_BAROTROPIC_EOS) {
    x1area.NewAthenaArray(pmb->ncells1 + 1);
    vol.NewAthenaArray(pmb->ncells1);
    x2aream.NewAthenaArray(pmb->ncells1);
    x2areap.NewAthenaArray(pmb->ncells1);
    x3aream.NewAthenaArray(pmb->ncells1);
    x3areap.NewAthenaArray(pmb->ncells1);
  }

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    Real phv = pmb->pcoord->x3v(k);
    Real sp = std::sin(phv), cp = std::cos(phv);
    for (int j=pmb->js; j<=pmb->je; ++j) {
      Real thv = pmb->pcoord->x2v(j);
      Real st = std::sin(thv), ct = std::cos(thv);
      if (NON_BAROTROPIC_EOS) {
        pmb->pcoord->Face1Area(k, j, pmb->is, pmb->ie + 1, x1area);
        pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
        pmb->pcoord->Face2Area(k, j,     pmb->is, pmb->ie, x2aream);
        pmb->pcoord->Face2Area(k, j + 1, pmb->is, pmb->ie, x2areap);
        pmb->pcoord->Face3Area(k,     j, pmb->is, pmb->ie, x3aream);
        pmb->pcoord->Face3Area(k + 1, j, pmb->is, pmb->ie, x3areap);
      }
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real rv = pmb->pcoord->x1v(i);
        // cell center in star-centered cartesian; planet at (r0, 0, 0)
        Real xc = rv*st*cp, yc = rv*st*sp, zc = rv*ct;
        Real dx = xc - r0, dy = yc, dz = zc;
        Real d2 = dx*dx + dy*dy + dz*dz;

        // ---------- planet gravity + indirect term (cartesian, then project) ----------
        Real g_scale = -ramp*gmp/std::pow(d2 + eps_soft*eps_soft, 1.5);
        Real ax = g_scale*dx + a_ind;
        Real ay = g_scale*dy;
        Real az = g_scale*dz;
        // spherical components: rhat, thhat, phhat at (theta, phi)
        Real a_r  = ax*st*cp + ay*st*sp + az*ct;
        Real a_th = ax*ct*cp + ay*ct*sp - az*st;
        Real a_ph = -ax*sp + ay*cp;
        Real rho_loc = prim(IDN,k,j,i);
        cons(IM1,k,j,i) += rho_loc*a_r*dt;
        cons(IM2,k,j,i) += rho_loc*a_th*dt;
        cons(IM3,k,j,i) += rho_loc*a_ph*dt;

        // ---------- conservative energy source via flux-form grad(Potential) ----------
        if (NON_BAROTROPIC_EOS) {
          Real rm = pmb->pcoord->x1f(i), rp = pmb->pcoord->x1f(i+1);
          Real thm = pmb->pcoord->x2f(j), thp = pmb->pcoord->x2f(j+1);
          Real phm = pmb->pcoord->x3f(k), php = pmb->pcoord->x3f(k+1);
          Real phic = PlanetPotential(xc, yc, zc, time);
          Real phil, phir;

          phil = PlanetPotential(rm*st*cp, rm*st*sp, rm*ct, time);
          phir = PlanetPotential(rp*st*cp, rp*st*sp, rp*ct, time);
          cons(IEN,k,j,i) -=
              dt*(pmb->phydro->flux[X1DIR](IDN,k,j,i+1)*x1area(i+1)*(phir - phic)
                  + pmb->phydro->flux[X1DIR](IDN,k,j,i)*x1area(i)*(phic - phil))/vol(i);

          Real stm = std::sin(thm), ctm = std::cos(thm);
          Real stp = std::sin(thp), ctp = std::cos(thp);
          phil = PlanetPotential(rv*stm*cp, rv*stm*sp, rv*ctm, time);
          phir = PlanetPotential(rv*stp*cp, rv*stp*sp, rv*ctp, time);
          cons(IEN,k,j,i) -=
              dt*(pmb->phydro->flux[X2DIR](IDN,k,j+1,i)*x2areap(i)*(phir - phic)
                  + pmb->phydro->flux[X2DIR](IDN,k,j,i)*x2aream(i)*(phic - phil))/vol(i);

          Real spm = std::sin(phm), cpm = std::cos(phm);
          Real spp = std::sin(php), cpp_ = std::cos(php);
          phil = PlanetPotential(rv*st*cpm, rv*st*spm, rv*ct, time);
          phir = PlanetPotential(rv*st*cpp_, rv*st*spp, rv*ct, time);
          cons(IEN,k,j,i) -=
              dt*(pmb->phydro->flux[X3DIR](IDN,k+1,j,i)*x3areap(i)*(phir - phic)
                  + pmb->phydro->flux[X3DIR](IDN,k,j,i)*x3aream(i)*(phic - phil))/vol(i);
        }

        // ---------- vacuum sink: clamp cells inside r_sink of the planet ----------
        if (r_sink > 0.0 && d2 < r_sink2) {
          cons(IDN,k,j,i) = sink_density;
          cons(IM1,k,j,i) = 0.0;
          cons(IM2,k,j,i) = 0.0;
          cons(IM3,k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            cons(IEN,k,j,i) = sink_density*p0_over_r0/(gamma_gas - 1.0);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! Distance-based AMR: refine blocks within amr_dist of the planet, derefine far away.

int RefinementCondition(MeshBlock *pmb) {
  if (amr_dist <= 0.0) return 0;
  Real d2min = std::numeric_limits<Real>::max();
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    Real phv = pmb->pcoord->x3v(k);
    Real sp = std::sin(phv), cp = std::cos(phv);
    for (int j=pmb->js; j<=pmb->je; ++j) {
      Real thv = pmb->pcoord->x2v(j);
      Real st = std::sin(thv), ct = std::cos(thv);
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real rv = pmb->pcoord->x1v(i);
        Real dx = rv*st*cp - r0, dy = rv*st*sp, dz = rv*ct;
        Real d2 = dx*dx + dy*dy + dz*dz;
        if (d2 < d2min) d2min = d2;
      }
    }
  }
  if (d2min < SQR(amr_dist)) return 1;
  if (d2min > SQR(2.5*amr_dist)) return -1;
  return 0;
}

//----------------------------------------------------------------------------------------
//! User BCs: pin ghost zones to the rotating-frame background (disk.cpp, spherical only)

namespace {
inline void FillDiskGhost(Coordinates *pco, AthenaArray<Real> &prim,
                          int k, int j, int i) {
  Real rad(0.0), phi(0.0), z(0.0);
  GetCylCoord(pco, rad, phi, z, i, j, k);
  prim(IDN,k,j,i) = DenProfileCyl(rad, phi, z);
  prim(IM1,k,j,i) = 0.0;
  prim(IM2,k,j,i) = 0.0;
  prim(IM3,k,j,i) = VelProfileCyl(rad, phi, z);
  if (NON_BAROTROPIC_EOS)
    prim(IEN,k,j,i) = PoverR(rad, phi, z)*prim(IDN,k,j,i);
}
} // namespace

void DiskInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=kl; k<=ku; ++k)
    for (int j=jl; j<=ju; ++j)
      for (int i=1; i<=ngh; ++i)
        FillDiskGhost(pco, prim, k, j, il-i);
}

void DiskOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=kl; k<=ku; ++k)
    for (int j=jl; j<=ju; ++j)
      for (int i=1; i<=ngh; ++i)
        FillDiskGhost(pco, prim, k, j, iu+i);
}

void DiskInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=kl; k<=ku; ++k)
    for (int j=1; j<=ngh; ++j)
      for (int i=il; i<=iu; ++i)
        FillDiskGhost(pco, prim, k, jl-j, i);
}

void DiskOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=kl; k<=ku; ++k)
    for (int j=1; j<=ngh; ++j)
      for (int i=il; i<=iu; ++i)
        FillDiskGhost(pco, prim, k, ju+j, i);
}

// Azimuthal faces.  Same background pin as x1/x2: the ghost zones hold the
// unperturbed disk, with VelProfileCyl already in the rotating frame (it
// subtracts rad*Omega0), so the inflowing edge carries the correct Keplerian
// shear past the planet.
//
// NOTE both faces are pinned, upstream and downstream alike.  That is what
// the local Hill-box runs do on their own side boundaries, and it is the
// point of these runs -- a pinned patch cannot deepen a gap.  The cost is
// that the planetary wake is truncated where it meets the downstream face,
// so the azimuthal extent must be wide enough that the wake has left the
// region of interest before it gets there.
void DiskInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=1; k<=ngh; ++k)
    for (int j=jl; j<=ju; ++j)
      for (int i=il; i<=iu; ++i)
        FillDiskGhost(pco, prim, kl-k, j, i);
}

void DiskOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=1; k<=ngh; ++k)
    for (int j=jl; j<=ju; ++j)
      for (int i=il; i<=iu; ++i)
        FillDiskGhost(pco, prim, ku+k, j, i);
}
