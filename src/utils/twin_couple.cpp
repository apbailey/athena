//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file twin_couple.cpp
//! \brief One-way coupling of two Mesh objects over a planet-centred spherical shell.

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "twin_couple.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//----------------------------------------------------------------------------------------
//! \brief Read the shell geometry, build the cell list and verify its thickness.

TwinCoupler::TwinCoupler(Mesh *psrc, Mesh *pdst, ParameterInput *pin)
    : psrc_(psrc), pdst_(pdst) {
  rmin_ = pin->GetReal("twin", "couple_rmin");
  rmax_ = pin->GetReal("twin", "couple_rmax");
  // The shell is centred on the planet, which the global disk pgens place at
  // (r0, 0, 0) in star-centred cartesian coordinates.
  Real r0 = pin->GetOrAddReal("problem", "r0", 1.0);
  cx_ = pin->GetOrAddReal("twin", "couple_x0", r0);
  cy_ = pin->GetOrAddReal("twin", "couple_y0", 0.0);
  cz_ = pin->GetOrAddReal("twin", "couple_z0", 0.0);

  if (rmax_ <= rmin_ || rmin_ < 0.0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in TwinCoupler" << std::endl
        << "need 0 <= <twin>/couple_rmin < couple_rmax, got " << rmin_
        << " and " << rmax_ << std::endl;
    ATHENA_ERROR(msg);
  }
  BuildIndex();
  CheckStencil();
}

//----------------------------------------------------------------------------------------
//! \brief Distance of a cell centre from the shell centre.  Ghost indices are valid:
//! their coordinates coincide with the neighbouring block's active cells.

Real TwinCoupler::PlanetRadius(MeshBlock *pmb, int k, int j, int i) const {
  Real r = pmb->pcoord->x1v(i);
  Real th = pmb->pcoord->x2v(j);
  Real ph = pmb->pcoord->x3v(k);
  Real st = std::sin(th), ct = std::cos(th);
  Real dx = r*st*std::cos(ph) - cx_;
  Real dy = r*st*std::sin(ph) - cy_;
  Real dz = r*ct - cz_;
  return std::sqrt(dx*dx + dy*dy + dz*dz);
}

//----------------------------------------------------------------------------------------
//! \brief Flag the active cells of every local block that fall inside the shell.

void TwinCoupler::BuildIndex() {
  idx_.resize(psrc_->nblocal);
  std::int64_t nforced = 0, ninner = 0;
  for (int b = 0; b < psrc_->nblocal; ++b) {
    MeshBlock *pmb = psrc_->my_blocks(b);
    idx_[b].clear();
    for (int k = pmb->ks; k <= pmb->ke; ++k) {
      for (int j = pmb->js; j <= pmb->je; ++j) {
        for (int i = pmb->is; i <= pmb->ie; ++i) {
          Real d = PlanetRadius(pmb, k, j, i);
          if (d < rmin_) {
            ninner++;
          } else if (d < rmax_) {
            idx_[b].push_back((k*pmb->ncells2 + j)*pmb->ncells1 + i);
            nforced++;
          }
        }
      }
    }
  }
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, &nforced, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &ninner, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (Globals::my_rank == 0) {
    std::cout << "TwinCoupler: shell " << rmin_ << " <= r_p < " << rmax_
              << " centred on (" << cx_ << ", " << cy_ << ", " << cz_ << ")"
              << std::endl
              << "  forced cells = " << nforced
              << ",  free cells inside r_min = " << ninner << std::endl;
  }
  if (nforced == 0 || ninner == 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in TwinCoupler" << std::endl
        << "the shell encloses " << ninner << " free cells and forces " << nforced
        << "; both must be nonzero." << std::endl;
    ATHENA_ERROR(msg);
  }
}

//----------------------------------------------------------------------------------------
//! \brief Verify the shell is thick enough to act as ghost zones for the free interior.
//!
//! Reconstruction reaches NGHOST cells along each coordinate direction, so no free cell
//! (r_p < rmin_) may have a stencil neighbour at r_p >= rmax_: such a neighbour evolves
//! freely in the destination Mesh and would leak the destination's own solution inward.

void TwinCoupler::CheckStencil() const {
  Real worst = 0.0;          // largest r_p reached by any free cell's stencil
  std::int64_t nbad = 0;
  for (int b = 0; b < psrc_->nblocal; ++b) {
    MeshBlock *pmb = psrc_->my_blocks(b);
    for (int k = pmb->ks; k <= pmb->ke; ++k) {
      for (int j = pmb->js; j <= pmb->je; ++j) {
        for (int i = pmb->is; i <= pmb->ie; ++i) {
          if (PlanetRadius(pmb, k, j, i) >= rmin_) continue;   // not a free cell
          Real reach = 0.0;
          for (int dk = -NGHOST; dk <= NGHOST; ++dk) {
            int kk = std::min(std::max(k + dk, 0), pmb->ncells3 - 1);
            for (int dj = -NGHOST; dj <= NGHOST; ++dj) {
              int jj = std::min(std::max(j + dj, 0), pmb->ncells2 - 1);
              for (int di = -NGHOST; di <= NGHOST; ++di) {
                int ii = std::min(std::max(i + di, 0), pmb->ncells1 - 1);
                reach = std::max(reach, PlanetRadius(pmb, kk, jj, ii));
              }
            }
          }
          worst = std::max(worst, reach);
          if (reach >= rmax_) nbad++;
        }
      }
    }
  }
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, &worst, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &nbad, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (Globals::my_rank == 0) {
    std::cout << "  NGHOST=" << NGHOST << " stencil of the free interior reaches r_p = "
              << worst << " (shell outer edge " << rmax_ << ")" << std::endl;
  }
  if (nbad > 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in TwinCoupler" << std::endl
        << nbad << " free cells have a reconstruction stencil reaching r_p = " << worst
        << " >= couple_rmax = " << rmax_ << "." << std::endl
        << "The forced shell is thinner than NGHOST=" << NGHOST
        << " cells somewhere; widen <twin>/couple_rmax." << std::endl;
    ATHENA_ERROR(msg);
  }
}

//----------------------------------------------------------------------------------------
//! \brief Overwrite the destination's shell cells with the source's current state.
//!
//! Both conserved and primitive variables are written: the next cycle reconstructs from
//! w and updates u, so leaving either stale would make the shell inconsistent.

void TwinCoupler::Apply() {
  for (int b = 0; b < psrc_->nblocal; ++b) {
    MeshBlock *ps = psrc_->my_blocks(b);
    MeshBlock *pd = pdst_->my_blocks(b);
    AthenaArray<Real> &us = ps->phydro->u, &ud = pd->phydro->u;
    AthenaArray<Real> &ws = ps->phydro->w, &wd = pd->phydro->w;
    const int nc1 = ps->ncells1, nc2 = ps->ncells2;
    for (std::size_t m = 0; m < idx_[b].size(); ++m) {
      int f = idx_[b][m];
      int i = f % nc1;
      int j = (f / nc1) % nc2;
      int k = f / (nc1*nc2);
      for (int n = 0; n < NHYDRO; ++n) {
        ud(n, k, j, i) = us(n, k, j, i);
        wd(n, k, j, i) = ws(n, k, j, i);
      }
    }
  }
}
