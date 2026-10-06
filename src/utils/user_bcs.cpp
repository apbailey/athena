//========================================================================================
// Athena++ astrophysical MHD code
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file user_bcs.cpp
//! \brief Runtime-selectable user boundary conditions: diode and vacuum.
//!
//! - "diode": outflow + clamp v_perp in the ghost zone to 0 on attempted inflow.
//!   Soft block; produces no spurious wave reflection.
//! - "vacuum": ghost cells set to (rho_floor, P_floor, v=0).  Magnetic field
//!   (if any) uses pure outflow so that div B = 0 is preserved at the face.
//!
//! Selection happens in the pgen via EnrollUserBCByName, which reads the
//! per-face BC name from <problem>.

// C++ headers
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>

// Athena++ headers
#include "user_bcs.hpp"

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals_interfaces.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

// File-scope vacuum-floor values (set by SetVacuumFloors).
namespace {
Real vacuum_rho_floor_   = 1.0e-6;
Real vacuum_press_floor_ = 1.0e-6;
constexpr Real kZero = 0.0;
}  // namespace

void SetVacuumFloors(Real rho_floor, Real press_floor) {
  vacuum_rho_floor_   = rho_floor;
  vacuum_press_floor_ = press_floor;
}

//========================================================================================
// Diode boundary conditions
//========================================================================================
// Pattern (per face):
//   1. Outflow copy of all NHYDRO primitive components into the ghost zone.
//   2. Clamp the perpendicular velocity in the ghost zone so that the boundary
//      cannot induce inflow.
//   3. If MHD is on, pure outflow on the three face-centered B components
//      (mirrors bvals/fc/outflow_fc.cpp exactly).
//========================================================================================

void DiodeInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int n=0; n<NHYDRO; ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          prim(n,k,j,is-i) = prim(n,k,j,is);
        }
      }
    }
  }
  // Inflow at inner_x1 means v1 > 0; clamp ghost v1 to <= 0.
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        prim(IVX,k,j,is-i) = std::min(prim(IVX,k,j,is), kZero);
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x1f(k,j,(is-i)) = b.x1f(k,j,is);
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x2f(k,j,(is-i)) = b.x2f(k,j,is);
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x3f(k,j,(is-i)) = b.x3f(k,j,is);
        }
      }
    }
  }
}

void DiodeOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int n=0; n<NHYDRO; ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          prim(n,k,j,ie+i) = prim(n,k,j,ie);
        }
      }
    }
  }
  // Inflow at outer_x1 means v1 < 0; clamp ghost v1 to >= 0.
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        prim(IVX,k,j,ie+i) = std::max(prim(IVX,k,j,ie), kZero);
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x1f(k,j,(ie+i+1)) = b.x1f(k,j,(ie+1));
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x2f(k,j,(ie+i)) = b.x2f(k,j,ie);
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x3f(k,j,(ie+i)) = b.x3f(k,j,ie);
        }
      }
    }
  }
}

void DiodeInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int n=0; n<NHYDRO; ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(n,k,js-j,i) = prim(n,k,js,i);
        }
      }
    }
  }
  // Inflow at inner_x2 means v2 > 0; clamp ghost v2 to <= 0.
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        prim(IVY,k,js-j,i) = std::min(prim(IVY,k,js,i), kZero);
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          b.x1f(k,(js-j),i) = b.x1f(k,js,i);
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x2f(k,(js-j),i) = b.x2f(k,js,i);
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x3f(k,(js-j),i) = b.x3f(k,js,i);
        }
      }
    }
  }
}

void DiodeOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int n=0; n<NHYDRO; ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(n,k,je+j,i) = prim(n,k,je,i);
        }
      }
    }
  }
  // Inflow at outer_x2 means v2 < 0; clamp ghost v2 to >= 0.
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        prim(IVY,k,je+j,i) = std::max(prim(IVY,k,je,i), kZero);
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          b.x1f(k,(je+j),i) = b.x1f(k,je,i);
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x2f(k,(je+j+1),i) = b.x2f(k,(je+1),i);
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x3f(k,(je+j),i) = b.x3f(k,je,i);
        }
      }
    }
  }
}

void DiodeInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int n=0; n<NHYDRO; ++n) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(n,ks-k,j,i) = prim(n,ks,j,i);
        }
      }
    }
  }
  // Inflow at inner_x3 means v3 > 0; clamp ghost v3 to <= 0.
  for (int k=1; k<=ngh; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        prim(IVZ,ks-k,j,i) = std::min(prim(IVZ,ks,j,i), kZero);
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          b.x1f((ks-k),j,i) = b.x1f(ks,j,i);
        }
      }
    }
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x2f((ks-k),j,i) = b.x2f(ks,j,i);
        }
      }
    }
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x3f((ks-k),j,i) = b.x3f(ks,j,i);
        }
      }
    }
  }
}

void DiodeOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                  FaceField &b, Real time, Real dt,
                  int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int n=0; n<NHYDRO; ++n) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(n,ke+k,j,i) = prim(n,ke,j,i);
        }
      }
    }
  }
  // Inflow at outer_x3 means v3 < 0; clamp ghost v3 to >= 0.
  for (int k=1; k<=ngh; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        prim(IVZ,ke+k,j,i) = std::max(prim(IVZ,ke,j,i), kZero);
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          b.x1f((ke+k),j,i) = b.x1f(ke,j,i);
        }
      }
    }
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x2f((ke+k),j,i) = b.x2f(ke,j,i);
        }
      }
    }
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x3f((ke+k+1),j,i) = b.x3f((ke+1),j,i);
        }
      }
    }
  }
}

//========================================================================================
// Vacuum boundary conditions
//========================================================================================
// Ghost cells get (rho_floor, 0, 0, 0, P_floor); face-centered B uses pure
// outflow (preserves div B = 0 at the face).
//========================================================================================

namespace {

void FillVacuumPrim(AthenaArray<Real> &prim,
                    int k, int j, int i) {
  prim(IDN,k,j,i) = vacuum_rho_floor_;
  prim(IVX,k,j,i) = 0.0;
  prim(IVY,k,j,i) = 0.0;
  prim(IVZ,k,j,i) = 0.0;
  if (NON_BAROTROPIC_EOS) {
    prim(IPR,k,j,i) = vacuum_press_floor_;
  }
}

}  // namespace

void VacuumInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt,
                   int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        FillVacuumPrim(prim, k, j, is-i);
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x1f(k,j,(is-i)) = b.x1f(k,j,is);
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x2f(k,j,(is-i)) = b.x2f(k,j,is);
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x3f(k,j,(is-i)) = b.x3f(k,j,is);
        }
      }
    }
  }
}

void VacuumOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt,
                   int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        FillVacuumPrim(prim, k, j, ie+i);
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x1f(k,j,(ie+i+1)) = b.x1f(k,j,(ie+1));
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x2f(k,j,(ie+i)) = b.x2f(k,j,ie);
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x3f(k,j,(ie+i)) = b.x3f(k,j,ie);
        }
      }
    }
  }
}

void VacuumInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt,
                   int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=is; i<=ie; ++i) {
        FillVacuumPrim(prim, k, js-j, i);
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          b.x1f(k,(js-j),i) = b.x1f(k,js,i);
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x2f(k,(js-j),i) = b.x2f(k,js,i);
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x3f(k,(js-j),i) = b.x3f(k,js,i);
        }
      }
    }
  }
}

void VacuumOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt,
                   int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=is; i<=ie; ++i) {
        FillVacuumPrim(prim, k, je+j, i);
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          b.x1f(k,(je+j),i) = b.x1f(k,je,i);
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x2f(k,(je+j+1),i) = b.x2f(k,(je+1),i);
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x3f(k,(je+j),i) = b.x3f(k,je,i);
        }
      }
    }
  }
}

void VacuumInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt,
                   int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int k=1; k<=ngh; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        FillVacuumPrim(prim, ks-k, j, i);
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          b.x1f((ks-k),j,i) = b.x1f(ks,j,i);
        }
      }
    }
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x2f((ks-k),j,i) = b.x2f(ks,j,i);
        }
      }
    }
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x3f((ks-k),j,i) = b.x3f(ks,j,i);
        }
      }
    }
  }
}

void VacuumOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt,
                   int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int k=1; k<=ngh; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        FillVacuumPrim(prim, ke+k, j, i);
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          b.x1f((ke+k),j,i) = b.x1f(ke,j,i);
        }
      }
    }
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x2f((ke+k),j,i) = b.x2f(ke,j,i);
        }
      }
    }
    for (int k=1; k<=ngh; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b.x3f((ke+k+1),j,i) = b.x3f((ke+1),j,i);
        }
      }
    }
  }
}

//========================================================================================
// Dispatcher: read <problem>/<param> and enroll the matching BC.
//========================================================================================

namespace {

BValFunc DiodeForFace(BoundaryFace face) {
  switch (face) {
    case BoundaryFace::inner_x1: return DiodeInnerX1;
    case BoundaryFace::outer_x1: return DiodeOuterX1;
    case BoundaryFace::inner_x2: return DiodeInnerX2;
    case BoundaryFace::outer_x2: return DiodeOuterX2;
    case BoundaryFace::inner_x3: return DiodeInnerX3;
    case BoundaryFace::outer_x3: return DiodeOuterX3;
    default:                     return nullptr;
  }
}

BValFunc VacuumForFace(BoundaryFace face) {
  switch (face) {
    case BoundaryFace::inner_x1: return VacuumInnerX1;
    case BoundaryFace::outer_x1: return VacuumOuterX1;
    case BoundaryFace::inner_x2: return VacuumInnerX2;
    case BoundaryFace::outer_x2: return VacuumOuterX2;
    case BoundaryFace::inner_x3: return VacuumInnerX3;
    case BoundaryFace::outer_x3: return VacuumOuterX3;
    default:                     return nullptr;
  }
}

}  // namespace

BValFunc ResolveUserBC(ParameterInput *pin, const char *param,
                       BoundaryFace face) {
  std::string name = pin->GetOrAddString("problem", param, "");
  if (name.empty()) return nullptr;

  // Default vacuum floors track the EOS hydro floors (<hydro>/dfloor,
  // <hydro>/pfloor).  Match the expression EOS uses (src/eos/adiabatic_hydro.cpp,
  // src/eos/isothermal_hydro.cpp) so the parameter is added with an identical
  // default if not present, and EOS subsequently reads back the same value.
  const Real float_min  = std::numeric_limits<float>::min();
  const Real hydro_default_floor = std::sqrt(1024.0 * float_min);
  Real rho_floor   = pin->GetOrAddReal("hydro", "dfloor", hydro_default_floor);
  Real press_floor = pin->GetOrAddReal("hydro", "pfloor", hydro_default_floor);
  // <problem>/vacuum_*_floor overrides if present.
  if (pin->DoesParameterExist("problem", "vacuum_rho_floor"))
    rho_floor   = pin->GetReal("problem", "vacuum_rho_floor");
  if (pin->DoesParameterExist("problem", "vacuum_press_floor"))
    press_floor = pin->GetReal("problem", "vacuum_press_floor");
  SetVacuumFloors(rho_floor, press_floor);

  if (name == "diode")  return DiodeForFace(face);
  if (name == "vacuum") return VacuumForFace(face);

  std::stringstream msg;
  msg << "[ResolveUserBC] unknown user BC name '" << name
      << "' for parameter '" << param << "'." << std::endl;
  ATHENA_ERROR(msg);
  return nullptr;
}
