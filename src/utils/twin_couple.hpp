#ifndef UTILS_TWIN_COUPLE_HPP_
#define UTILS_TWIN_COUPLE_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file twin_couple.hpp
//! \brief One-way coupling of two Mesh objects over a planet-centred spherical shell.
//!
//! Used by the twin-mesh experiment: a destination Mesh (e.g. one carrying a diode
//! midplane) has its cells in \f$r_{\rm min}\le r_p<r_{\rm max}\f$ overwritten every
//! cycle with the corresponding cells of a source Mesh (the unmodified run), so the
//! interior \f$r_p<r_{\rm min}\f$ evolves under exactly the source's boundary state.
//!
//! The shell is a Dirichlet layer, so it plays the role of ghost zones for the free
//! interior and must be at least NGHOST cells thick along every coordinate direction.
//! That is checked at construction, not assumed.

#include <vector>

#include "../athena.hpp"

class Mesh;
class ParameterInput;

//! \class TwinCoupler
//! \brief Copies hydro state from one Mesh to another over a fixed set of cells.

class TwinCoupler {
 public:
  TwinCoupler(Mesh *psrc, Mesh *pdst, ParameterInput *pin);
  //! Overwrite the destination's shell cells with the source's current values.
  void Apply();

 private:
  Mesh *psrc_, *pdst_;
  Real rmin_, rmax_;          // planet-centred radii bounding the forced shell
  Real cx_, cy_, cz_;         // shell centre in star-centred cartesian
  std::vector<std::vector<int>> idx_;   // per local block: flattened cell indices
  Real PlanetRadius(MeshBlock *pmb, int k, int j, int i) const;
  void BuildIndex();
  void CheckStencil() const;
};

#endif // UTILS_TWIN_COUPLE_HPP_
