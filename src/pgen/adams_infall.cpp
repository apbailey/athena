//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file adams_infall.cpp
//! \brief 2D problem generator for the Adams & Batygin (2022) infall model.
//!
//! Ballistic, pressure-free collapse of gas from a circumstellar disk onto a young
//! planet inside its Hill sphere, set up in spherical_polar coordinates (r, theta).
//! The analytic infall solution (eqs. 5, 9-12 of the paper) is used as the initial
//! condition everywhere and as a Dirichlet inflow boundary at x1-outer. Other faces:
//!   x1-inner: outflow (planet accretes)
//!   x2-inner: polar_wedge (pole)
//!   x2-outer: outflow (material leaves the polar wedge toward the disk midplane)
//!
//! Code units: H = c_s = rho_0 = 1, Omega_planet = 1. The user supplies the
//! dimensionless thermal mass qthermal = GM_p directly in <problem>.

// C headers

// C++ headers
#include <algorithm>  // max
#include <cmath>      // cos, sin, sqrt, cbrt, pow, atan
#include <iostream>
#include <sstream>
#include <string>

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
#include "../parameter_input.hpp"

namespace {
Real gm_p;       // GM_p in code units (== qthermal)
Real r_hill;     // R_H
Real r_cent;     // R_C = lambda^2 * R_H / 3  (lambda=1 -> Adams 2022)
Real rho_norm;   // density prefactor: gives rho=1 at (r_outer, pole) for ISOTROPIC
Real gamma_idx;  // adiabatic gamma (set only if NON_BAROTROPIC_EOS)
Real lambda_J;   // Adams 2025 angular-momentum bias: J = lambda * Omega * R_H^2
bool hill_source = false;        // Coriolis + tidal (-1.5 R^2 <cos^2 phi>) in the rotating frame
bool vertical_gravity = false;   // tidal vertical compression (+0.5 z^2)

// Taylor & Adams (2024) Icarus 415, 116044 eq. (1): five inflow geometries
// parameterized by mu_0 = cos(theta_0). All five are normalized so that the
// total mass inflow integrated over the full 4 pi solid angle equals Mdot,
// which is what rho_norm absorbs. So changing flux_kind redistributes ρ
// across mu_0 but keeps the total injected Mdot the same.
enum FluxKind { FLUX_POLAR, FLUX_QUASIPOLAR, FLUX_ISOTROPIC,
                FLUX_QUASIEQUATORIAL, FLUX_EQUATORIAL, FLUX_GAUSSIAN,
                FLUX_GAUSSIAN_NORMALIZED };
FluxKind flux_kind = FLUX_ISOTROPIC;
Real r_out_boundary;     // = x1max; set in InitUserMeshData for Gaussian profile
Real gaussian_norm = 1.0;// integral of raw f_eff_gaussian over mu_0 in [0,1]

constexpr Real kTiny = 1.0e-14;

// f_i(mu0) -- the inflow asymmetry function from T&A 2024 eq. (1).
// Isotropic returns 1 (reproduces Adams & Batygin 2022).
static Real FluxFactor(Real mu0) {
  Real m2 = mu0 * mu0;
  switch (flux_kind) {
    case FLUX_POLAR:
      return 3.0 * m2;
    case FLUX_QUASIPOLAR:
      return 2.0 * mu0;                       // mu0 >= 0 in upper hemisphere
    case FLUX_ISOTROPIC:
      return 1.0;
    case FLUX_QUASIEQUATORIAL:
      return (4.0 / PI) * std::sqrt(std::max(1.0 - m2, 0.0));
    case FLUX_EQUATORIAL:
      return 1.5 * (1.0 - m2);
    case FLUX_GAUSSIAN:
    case FLUX_GAUSSIAN_NORMALIZED: {
      // Effective f for an isothermal-disk vertical Gaussian rho = exp(-z^2 / (2 H^2))
      // imposed at the outer boundary r = r_out_boundary. With c_s = H = 1 in code units.
      // _NORMALIZED variant divides by integral over mu0 in [0,1] so total Mdot
      // matches the T&A geometries.
      Real zeta_out = r_cent / r_out_boundary;
      Real omu0sq = 1.0 - m2;
      Real mu_at_RH = mu0 * (1.0 - zeta_out * omu0sq);
      Real z_at_RH = r_out_boundary * mu_at_RH;
      Real vr_RH = std::sqrt(gm_p / r_out_boundary * (2.0 - zeta_out * omu0sq));
      Real P2_mu0 = 0.5 * (3.0 * m2 - 1.0);
      Real factor = r_out_boundary * r_out_boundary * vr_RH
                    * (1.0 + 2.0 * zeta_out * P2_mu0);
      Real f_raw = std::exp(-0.5 * z_at_RH * z_at_RH) * factor / rho_norm;
      if (flux_kind == FLUX_GAUSSIAN_NORMALIZED) f_raw /= gaussian_norm;
      return f_raw;
    }
  }
  return 1.0;  // unreachable
}

// Closed-form Cardano solution to the orbit-equation cubic
//   zeta * mu0^3 + (1 - zeta) * mu0 - mu = 0
// Ported from analytic.py (functions mu_crit, Q, K, mu0).
//   - zeta < 1, or zeta >= 1 with mu^2 >= mu_crit^2: single-real-root radical branch.
//   - zeta >= 1 and mu^2 < mu_crit^2: trigonometric branch (three real roots).
static Real SolveMu0(Real mu, Real zeta) {
  // zeta -> 0 limit (e.g., lambda = 0): cubic reduces to mu0 = mu (no deflection).
  if (zeta < kTiny) return mu;
  Real mu_crit_sq = (zeta >= 1.0)
      ? 4.0 * std::pow(zeta - 1.0, 3) / (27.0 * zeta)
      : -1.0;
  if (zeta < 1.0 || (mu*mu) >= mu_crit_sq) {
    Real inside = mu*mu - 4.0 * std::pow(zeta - 1.0, 3) / (27.0 * zeta);
    Real sqrt_inside = std::sqrt(std::max(inside, 0.0));
    Real cube_arg = mu + sqrt_inside;
    Real Q = 3.0 * std::pow(zeta, 2.0/3.0) * std::cbrt(cube_arg);
    if (std::abs(Q) < kTiny) return 0.0;
    return Q / (3.0 * std::cbrt(2.0) * zeta)
           + std::cbrt(2.0) * (zeta - 1.0) / Q;
  } else {
    if (std::abs(mu) < kTiny) {
      return std::sqrt((zeta - 1.0) / zeta);
    }
    Real K = std::sqrt(108.0 * std::pow(zeta - 1.0, 3) / zeta - 729.0 * mu * mu)
             / (27.0 * mu);
    return 2.0 * std::sqrt((zeta - 1.0) / (3.0 * zeta))
           * std::cos(std::atan(K) / 3.0);
  }
}

// Evaluate the Adams-Batygin infall state (rho, v_r, v_theta, v_phi, press)
// at a point (r, theta). Pressure is set to rho (sound speed = 1 in code units)
// when NON_BAROTROPIC_EOS; otherwise the press output is unused.
static void InfallState(Real r, Real theta,
                        Real &rho, Real &vr, Real &vtheta, Real &vphi,
                        Real &press) {
  Real mu = std::cos(theta);
  Real zeta = r_cent / r;
  Real mu0 = SolveMu0(mu, zeta);
  Real v0 = std::sqrt(gm_p / r);
  Real one_minus_mu0sq = 1.0 - mu0 * mu0;
  Real one_minus_musq  = 1.0 - mu * mu;

  // eq. 9 (radial)
  Real radicand = 2.0 - zeta * one_minus_mu0sq;
  vr = -v0 * std::sqrt(std::max(radicand, 0.0));

  // eqs. 10, 11 (polar, azimuthal). NOTE: eq. 10 as printed in Adams & Batygin
  // (2022) is MISSING a factor of zeta inside the square root; with the
  // correction included here the three velocity components satisfy
  // v_r^2 + v_theta^2 + v_phi^2 = 2 v_0^2 (parabolic E=0). Verified by
  // re-deriving v_theta from the Cassen-Moosman formula via the orbit eq.
  if (one_minus_musq < kTiny) {
    vtheta = 0.0;
    vphi   = 0.0;
  } else {
    Real diff = mu0 * mu0 - mu * mu;       // > 0 for mu0 > mu (the physical case)
    vtheta = v0 * std::sqrt(std::max(zeta * one_minus_mu0sq * diff / one_minus_musq, 0.0));
    // v_phi sign follows sign(lambda_J): retrograde (lambda<0) gives v_phi<0.
    // sqrt(zeta) absorbs |lambda|, the sign is restored here.
    Real sign_lam = (lambda_J > 0.0) - (lambda_J < 0.0);
    vphi   = sign_lam * v0 * one_minus_mu0sq / std::sqrt(one_minus_musq) * std::sqrt(zeta);
  }

  // density from T&A 2024 eq. (7): rho = (Mdot/(4 pi)) f_i(mu0) / (r^2 |v_r|)
  //                                       / [1 + 2 zeta P_2(mu0)]
  // The Mdot/(4 pi) factor is absorbed into rho_norm (calibrated for isotropic
  // f=1; the same value of rho_norm gives the SAME total Mdot for every flux
  // function, since each f_i integrates to 1 over mu0 in [0,1]).
  Real P2 = 0.5 * (3.0 * mu0 * mu0 - 1.0);
  Real denom = r * r * std::abs(vr) * (1.0 + 2.0 * zeta * P2);
  rho = (denom > 0.0) ? rho_norm * FluxFactor(mu0) / denom : 0.0;

  if (NON_BAROTROPIC_EOS) {
    press = rho;  // c_s = 1, isothermal-equivalent pressure
  } else {
    press = 0.0;
  }
}
}  // namespace

// Forward declarations (file-scope linkage).
void OuterX1Infall(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt,
                   int is, int ie, int js, int je, int ks, int ke, int ngh);
void InnerX1Diode(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int is, int ie, int js, int je, int ks, int ke, int ngh);
void OuterX2Diode(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int is, int ie, int js, int je, int ks, int ke, int ngh);

void PointMassGravity(MeshBlock *pmb, const Real time, const Real dt,
                      const AthenaArray<Real> &prim,
                      const AthenaArray<Real> &prim_scalar,
                      const AthenaArray<Real> &bcc,
                      AthenaArray<Real> &cons,
                      AthenaArray<Real> &cons_scalar);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  gm_p = pin->GetReal("problem", "qthermal");
  r_hill = std::cbrt(gm_p / 3.0);
  // Adams 2025 angular momentum bias factor: J = lambda * Omega * R_H^2.
  // R_C = J^2 / GM_p = lambda^2 * R_H / 3. lambda=1 reproduces Adams 2022.
  lambda_J = pin->GetOrAddReal("problem", "lambda", 1.0);
  r_cent = lambda_J * lambda_J * r_hill / 3.0;
  hill_source      = pin->GetOrAddBoolean("problem", "hill_source", false);
  vertical_gravity = pin->GetOrAddBoolean("problem", "vertical_gravity", false);

  // Taylor & Adams 2024 inflow geometry. Default = isotropic (reproduces
  // Adams & Batygin 2022).
  std::string fkind = pin->GetOrAddString("problem", "flux_function", "isotropic");
  if      (fkind == "polar")            flux_kind = FLUX_POLAR;
  else if (fkind == "quasipolar")       flux_kind = FLUX_QUASIPOLAR;
  else if (fkind == "isotropic")        flux_kind = FLUX_ISOTROPIC;
  else if (fkind == "quasiequatorial")  flux_kind = FLUX_QUASIEQUATORIAL;
  else if (fkind == "equatorial")       flux_kind = FLUX_EQUATORIAL;
  else if (fkind == "gaussian")             flux_kind = FLUX_GAUSSIAN;
  else if (fkind == "gaussian_normalized")  flux_kind = FLUX_GAUSSIAN_NORMALIZED;
  else {
    std::stringstream msg;
    msg << "### FATAL ERROR in adams_infall::InitUserMeshData" << std::endl
        << "Unknown problem/flux_function = \"" << fkind << "\"" << std::endl
        << "Valid: polar, quasipolar, isotropic, quasiequatorial, equatorial,"
        << " gaussian, gaussian_normalized"
        << std::endl;
    ATHENA_ERROR(msg);
  }

  if (NON_BAROTROPIC_EOS) {
    gamma_idx = pin->GetReal("hydro", "gamma");
  } else {
    gamma_idx = 0.0;
  }

  // Normalization: at (r_outer, theta = 0 -> mu = 1), the cubic has the exact
  // root mu0 = 1, giving v_r = -v0 sqrt(2) and density factor 1/(1 + 2 zeta_out).
  // Set rho_norm so that rho(r_outer, pole) = 1.
  Real r_out = mesh_size.x1max;
  r_out_boundary = r_out;     // needed by FluxFactor for the gaussian case
  Real v0_out = std::sqrt(gm_p / r_out);
  Real zeta_out = r_cent / r_out;
  rho_norm = r_out * r_out * v0_out * std::sqrt(2.0) * (1.0 + 2.0 * zeta_out);

  // Numerically integrate the raw f_eff_gaussian over mu0 in [0,1] so the
  // gaussian_normalized variant matches T&A normalization (integral = 1).
  // Use composite Simpson's rule over a temporary set to FLUX_GAUSSIAN.
  {
    FluxKind saved = flux_kind;
    flux_kind = FLUX_GAUSSIAN;
    const int N = 4001;
    Real h = 1.0 / (N - 1);
    Real integral = 0.0;
    for (int k = 0; k < N; ++k) {
      Real mu0 = k * h;
      Real w = (k == 0 || k == N - 1) ? 1.0 : ((k % 2 == 1) ? 4.0 : 2.0);
      integral += w * FluxFactor(mu0);
    }
    gaussian_norm = integral * h / 3.0;
    flux_kind = saved;
  }

  if (Globals::my_rank == 0) {
    std::cout << "adams_infall: qthermal = " << gm_p
              << ", R_H = " << r_hill
              << ", R_C = " << r_cent
              << ", zeta_outer = " << zeta_out
              << ", rho_norm = " << rho_norm
              << ", flux_function = " << fkind
              << ", gaussian_norm = " << gaussian_norm
              << ", lambda = " << lambda_J
              << ", hill_source = " << (hill_source ? "true" : "false")
              << ", vertical_gravity = " << (vertical_gravity ? "true" : "false")
              << std::endl;
  }

  EnrollUserExplicitSourceFunction(PointMassGravity);

  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, OuterX1Infall);
  }
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, InnerX1Diode);
  }
  if (mesh_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x2, OuterX2Diode);
  }
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initialize each cell with the analytic Adams-Batygin infall state.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      Real theta = pcoord->x2v(j);
      for (int i = is; i <= ie; ++i) {
        Real r = pcoord->x1v(i);
        Real rho, vr, vth, vph, press;
        InfallState(r, theta, rho, vr, vth, vph, press);
        phydro->u(IDN, k, j, i) = rho;
        phydro->u(IM1, k, j, i) = rho * vr;
        phydro->u(IM2, k, j, i) = rho * vth;
        phydro->u(IM3, k, j, i) = rho * vph;
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN, k, j, i) = press / (gamma_idx - 1.0)
              + 0.5 * rho * (vr*vr + vth*vth + vph*vph);
        }
      }
    }
  }
  return;
}

//========================================================================================
//! \fn void OuterX1Infall(...)
//! \brief Dirichlet inflow at the outer radial boundary, set from the analytic
//!        infall solution evaluated at each ghost cell.
//========================================================================================

void OuterX1Infall(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt,
                   int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      Real theta = pco->x2v(j);
      for (int i = 1; i <= ngh; ++i) {
        Real r = pco->x1v(ie + i);
        Real rho, vr, vth, vph, press;
        InfallState(r, theta, rho, vr, vth, vph, press);
        prim(IDN, k, j, ie + i) = rho;
        prim(IVX, k, j, ie + i) = vr;
        prim(IVY, k, j, ie + i) = vth;
        prim(IVZ, k, j, ie + i) = vph;
        if (NON_BAROTROPIC_EOS) {
          prim(IPR, k, j, ie + i) = press;
        }
      }
    }
  }
  return;
}

//========================================================================================
//! \fn Real Potential(Real r, Real theta)
//! \brief Total gravitational potential per unit mass used by PointMassGravity.
//!   Includes the planet point-mass and (if hill_source) the axisymmetric Hill
//!   tidal and vertical terms. Omega = 1 in code units. sq = <cos^2 phi> = 1/2.
//========================================================================================

namespace {
inline Real Potential(Real r, Real theta) {
  Real pot = -gm_p / r;
  if (hill_source) {
    const Real sq = 0.5;
    Real R = r * std::sin(theta);
    pot += -1.5 * R * R * sq;
  }
  if (vertical_gravity) {
    Real z = r * std::cos(theta);
    pot += 0.5 * z * z;
  }
  return pot;
}
}  // namespace

//========================================================================================
//! \fn void PointMassGravity(...)
//! \brief Momentum sources from -rho grad(Potential) plus (if hill_source)
//!        Coriolis. Energy sources are computed CONSERVATIVELY via the
//!        flux-divergence form -div(rho v Phi), matching bondi_accretion_unified.
//========================================================================================

void PointMassGravity(MeshBlock *pmb, const Real time, const Real dt,
                      const AthenaArray<Real> &prim,
                      const AthenaArray<Real> &prim_scalar,
                      const AthenaArray<Real> &bcc,
                      AthenaArray<Real> &cons,
                      AthenaArray<Real> &cons_scalar) {
  AthenaArray<Real> x1area, vol, x2aream, x2areap;
  x1area.NewAthenaArray(pmb->ncells1 + 1);
  vol.NewAthenaArray(pmb->ncells1);
  x2aream.NewAthenaArray(pmb->ncells1);
  x2areap.NewAthenaArray(pmb->ncells1);

  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      Real x2  = pmb->pcoord->x2v(j);
      Real x2m = pmb->pcoord->x2f(j);
      Real x2p = pmb->pcoord->x2f(j+1);
      Real sth = std::sin(x2), cth = std::cos(x2);
      pmb->pcoord->Face1Area(k, j, pmb->is, pmb->ie + 1, x1area);
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Face2Area(k, j,     pmb->is, pmb->ie, x2aream);
        pmb->pcoord->Face2Area(k, j + 1, pmb->is, pmb->ie, x2areap);
      }
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real x1  = pmb->pcoord->x1v(i);
        Real x1m = pmb->pcoord->x1f(i);
        Real x1p = pmb->pcoord->x1f(i+1);
        Real rho = prim(IDN, k, j, i);

        // Conservative energy update from grad(Potential), per direction.
        if (NON_BAROTROPIC_EOS) {
          Real phil, phic, phir;
          phil = Potential(x1m, x2);
          phic = Potential(x1,  x2);
          phir = Potential(x1p, x2);
          cons(IEN, k, j, i) -= dt * (
              pmb->phydro->flux[X1DIR](IDN, k, j, i+1) * x1area(i+1) * (phir - phic)
            + pmb->phydro->flux[X1DIR](IDN, k, j, i)   * x1area(i)   * (phic - phil)
          ) / vol(i);
          if (pmb->block_size.nx2 > 1) {
            phil = Potential(x1, x2m);
            phic = Potential(x1, x2);
            phir = Potential(x1, x2p);
            cons(IEN, k, j, i) -= dt * (
                pmb->phydro->flux[X2DIR](IDN, k, j+1, i) * x2areap(i) * (phir - phic)
              + pmb->phydro->flux[X2DIR](IDN, k, j, i)   * x2aream(i) * (phic - phil)
            ) / vol(i);
          }
        }

        // Momentum sources -rho grad(Potential).
        // Point-mass radial:
        Real g_r = -gm_p / (x1 * x1);
        cons(IM1, k, j, i) += dt * rho * g_r;

        // Coriolis + tidal-radial momentum sources (axisymmetric average).
        // Coriolis is velocity-dependent and does no work.
        if (hill_source) {
          const Real sq = 0.5;
          Real vr  = prim(IVX, k, j, i);
          Real vth = prim(IVY, k, j, i);
          Real vph = prim(IVZ, k, j, i);
          // Coriolis (Omega = z-hat):
          cons(IM1, k, j, i) += dt * 2.0 * sth * rho * vph;
          cons(IM2, k, j, i) += dt * 2.0 * cth * rho * vph;
          cons(IM3, k, j, i) -= dt * (2.0 * sth * rho * vr + 2.0 * cth * rho * vth);
          // Tidal radial component from -1.5 R^2 <cos^2 phi>:
          cons(IM1, k, j, i) += dt * 3.0 * rho * x1 * sth * sth * sq;
          cons(IM2, k, j, i) += dt * 3.0 * rho * x1 * sth * cth * sq;
        }
        // Vertical-compressional tidal gravity: g_z = -z, potential +0.5 z^2.
        if (vertical_gravity) {
          Real z = x1 * cth;
          cons(IM1, k, j, i) -= dt * rho * z * cth;
          cons(IM2, k, j, i) += dt * rho * z * sth;
        }
      }
    }
  }
  vol.DeleteAthenaArray();
  x1area.DeleteAthenaArray();
  x2aream.DeleteAthenaArray();
  x2areap.DeleteAthenaArray();
  return;
}

//========================================================================================
//! \fn void InnerX1Diode(...)
//! \brief Outflow-only BC at inner-x1 (planet boundary): copy interior, clamp
//!        radial velocity to v_r <= 0 so no gas can flow back outward.
//========================================================================================

void InnerX1Diode(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int n = 0; n < NHYDRO; ++n) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = 1; i <= ngh; ++i) {
          prim(n, k, j, is - i) = prim(n, k, j, is);
        }
      }
    }
  }
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = 1; i <= ngh; ++i) {
        if (prim(IVX, k, j, is - i) > 0.0) prim(IVX, k, j, is - i) = 0.0;
      }
    }
  }
  return;
}

//========================================================================================
//! \fn void OuterX2Diode(...)
//! \brief Outflow-only BC at outer-x2 (equatorial side of wedge): copy interior,
//!        clamp polar velocity to v_theta >= 0 so no gas can flow back toward pole.
//========================================================================================

void OuterX2Diode(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int n = 0; n < NHYDRO; ++n) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = 1; j <= ngh; ++j) {
        for (int i = is; i <= ie; ++i) {
          prim(n, k, je + j, i) = prim(n, k, je, i);
        }
      }
    }
  }
  for (int k = ks; k <= ke; ++k) {
    for (int j = 1; j <= ngh; ++j) {
      for (int i = is; i <= ie; ++i) {
        if (prim(IVY, k, je + j, i) < 0.0) prim(IVY, k, je + j, i) = 0.0;
      }
    }
  }
  return;
}
