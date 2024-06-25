//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file exact.cpp
//! \brief Exact solver for the Riemann problem
//!
//!  Computes 1D fluxes by solving the exact problem iteratively. Currently only for
//!  isothermal setup.
//!
//! REFERENCES:
//! - R.J. LeVeque, "Numerical Methods for Conservation Laws", 2nd ed.,
//!   Birkhauser Verlag, Basel, (1992).
//! - E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics", 2nd ed.,
//!   Springer-Verlag, Berlin, (1999) chpt. 10.

// C headers

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../eos/eos.hpp"
#include "../../hydro.hpp"

static void srder(double cs, double dm, double vl, double vr, double dmin, double dmax, double &y, double &dydx);
static double rtsafe(double cs, double x1, double x2, double xacc, double vl, double vr, double dmin, double dmax);

//----------------------------------------------------------------------------------------
//! \fn void Hydro::RiemannSolver
//! \brief The exact Riemann solver for hydrodynamics (isothermal)

void Hydro::RiemannSolver(const int k, const int j, const int il, const int iu,
  const int ivx, AthenaArray<Real> &wl,
  AthenaArray<Real> &wr, AthenaArray<Real> &flx,
  const AthenaArray<Real> &dxw) {
    int ivy = IVX + ((ivx-IVX)+1)%3;
    int ivz = IVX + ((ivx-IVX)+2)%3;
    Real wli[(NHYDRO)],wri[(NHYDRO)];
    Real flxi[(NHYDRO)];
    Real iso_cs = pmy_block->peos->GetIsoSoundSpeed();
    Real tmp,zm,dm,zl,zr,dmin,dmax,vxm,mxm,sl,sr,hdl,tll,hdr,tlr;
    int soln;
    bool fassign; // True = flux has been assigned, False flux has not been assigned

    for (int i=il; i<=iu; ++i) {
      fassign = false;
      //--- Step 1.  Load L/R states into local variables
      wli[IDN]=wl(IDN,i);
      wli[IVX]=wl(ivx,i);
      wli[IVY]=wl(ivy,i);
      wli[IVZ]=wl(ivz,i);

      wri[IDN]=wr(IDN,i);
      wri[IVX]=wr(ivx,i);
      wri[IVY]=wr(ivy,i);
      wri[IVZ]=wr(ivz,i);

      //--- Step 2. Compute the density and momentum of the intermediate state
      zl = std::sqrt(wli[IDN]);
      zr = std::sqrt(wri[IDN]);

      // 1-shock and 2-shock
      soln = 0;

      // start by finding density if shocks on both left and right.
      // This will only be the case if dm > wli[IDN] and dm > wri[IDN]
      tmp = zl*zr*(wli[IVX] - wri[IVX])/(2.0*iso_cs*(zl + zr));
      zm = tmp + std::sqrt(tmp*tmp + zl*zr);
      dm = zm*zm;

      // Get velocity from 1-shock formula
      vxm = wli[IVX] - iso_cs*(dm - wli[IDN])/(zm*zl);

      // If left or right density is greater than intermediate density,
      // then at least one side has rarefaction instead of shock
      dmin = std::min(wli[IDN], wri[IDN]);
      dmax = std::max(wli[IDN], wri[IDN]);
      if (dm < dmax){
        // 1-rarefaction and 2-rarefaction
        soln = 3;

        // Try rarefactions on both left and right, since it's a quicker
        // calculation than 1-shock+2-rarefaction or 1-raref+2-shock
        dm = zl*zr*std::exp((wli[IVX]-wri[IVX])/(2.0*iso_cs));

        // Get velocity from 1-rarefaction formula
        vxm = wli[IVX] - iso_cs*std::log(dm/wli[IDN]);

        // If left or right density is small than intermediate density,
        // then we must instead have a combination of shock and rarefaction
        if (dm > dmin){
          // EITHER 1-rarefaction and 2-shock
          // OR     1-shock and 2-rarefaction

          // Solve iteratively equation for shock and rarefaction
          // If wli[IDN] > wri[IDN] ==> 1-rarefaction and 2-shock
          // If wri[IDN] > wli[IDN] ==> 1-shock and 2-rarefaction
          if (wli[IDN] > wri[IDN])
          soln = 2;
          else
          soln = 1;

          dm = rtsafe(iso_cs, dmin, dmax, 2.0*std::numeric_limits<double>::epsilon(), wli[IVX], wri[IVX], dmin, dmax);

          // Don't be foolish enough to take ln of zero
          if ((dm > dmin) && (dm <= dmax)){
            if (wli[IDN] > wri[IDN]){
              // Get velocity from 1-rarefaction formula
              vxm = wli[IVX] - iso_cs*std::log(dm/wli[IDN]);
            }
            else{
              // Get velocity from 2-rarefaction formula
              vxm = wri[IVX] + iso_cs*std::log(dm/wri[IDN]);
            }
          }
          else{
            // DEFAULT 1-rarefaction and 2-rarefaction
            soln = 3;

            // In the event that the intermediate density fails to fall between
            // the left and right densities (should only happen when left and
            // right densities differ only slightly and intermediate density
            // calculated in any step has significant truncation and/or roundoff
            // errors), default to rarefactions on both left and right
            dm = zl*zr*std::exp((wli[IVX]-wri[IVX])/(2.0*iso_cs));

            // Get velocity from 1-rarefaction formula
            vxm = wli[IVX] - iso_cs*std::log(dm/wli[IDN]);
          }
        }
      }

      //--- Step 2. Calculate the Interface Flux if the wave speeds are such
      //            that we aren't actually in the intermediate state

      if (soln & 2){ // left rarefaction, cases 2 & 3
        // The L-going rarefaction head/tail velocity
        hdl = wli[IVX] - iso_cs;
        tll = vxm - iso_cs;

        if (hdl >= 0.0){
          // To left of rarefaction
          flxi[IDN] = wli[IDN]*wli[IVX];
          flxi[IVX] = wli[IDN]*wli[IVX]*wli[IVX] + wli[IDN]*iso_cs*iso_cs;
          flxi[IVY] = wli[IDN]*wli[IVY]*wli[IVX];
          flxi[IVZ] = wli[IDN]*wli[IVZ]*wli[IVX];
          fassign = true;
        }
        else if (tll >= 0.0){
          // Inside rarefaction fan
          dm = wli[IDN]*std::exp(hdl/iso_cs);
          mxm = wli[IDN]*iso_cs*std::exp(hdl/iso_cs);
          vxm = (dm == 0.0 ? 0.0 : mxm/dm);

          flxi[IDN] = mxm;
          flxi[IVX] = mxm*vxm + dm*iso_cs*iso_cs;
          flxi[IVY] = mxm*wli[IVY];
          flxi[IVZ] = mxm*wli[IVZ];
          fassign = true;
        }
      }
      else{ // left shock
        // The L-going shock velocity
        sl = wli[IVX] - iso_cs*std::sqrt(dm)/zl;

        if (sl >= 0.0){
          flxi[IDN] = wli[IDN]*wli[IVX];
          flxi[IVX] = wli[IDN]*wli[IVX]*wli[IVX] + wli[IDN]*iso_cs*iso_cs;
          flxi[IVY] = wli[IDN]*wli[IVY]*wli[IVX];
          flxi[IVZ] = wli[IDN]*wli[IVZ]*wli[IVX];
          fassign = true;
        }
      }
      if (soln & 1){ // right rarefaction, cases 1 & 3
        hdr = wri[IVX] + iso_cs;
        tlr = vxm + iso_cs;

        if (hdr <= 0.0){
          flxi[IDN] = wri[IDN]*wri[IVX];
          flxi[IVX] = wri[IDN]*wri[IVX]*wri[IVX] + wri[IDN]*iso_cs*iso_cs;
          flxi[IVY] = wri[IDN]*wri[IVY]*wri[IVX];
          flxi[IVZ] = wri[IDN]*wri[IVZ]*wri[IVX];
          fassign = true;
        }
        else if (tlr <= 0.0){
          // Inside rarefaction fan
          tmp = dm;
          dm = tmp*std::exp(-tlr/iso_cs);
          mxm = -tmp*iso_cs*std::exp(-tlr/iso_cs);
          vxm = (dm == 0.0 ? 0.0 : mxm/dm);

          flxi[IDN] = mxm;
          flxi[IVX] = mxm*vxm + dm*iso_cs*iso_cs;
          flxi[IVY] = mxm*wri[IVY];
          flxi[IVZ] = mxm*wri[IVZ];
          fassign = true;
        }
      }
      else{ // right shock
        // The R-going shock velocity
        sr = wri[IVX] + iso_cs*std::sqrt(dm)/zr;

        if (sr <= 0.0){
          // To right of shock
          flxi[IDN] = wri[IDN]*wri[IVX];
          flxi[IVX] = wri[IDN]*wri[IVX]*wri[IVX] + wri[IDN]*iso_cs*iso_cs;
          flxi[IVY] = wri[IDN]*wri[IVX]*wri[IVY];
          flxi[IVZ] = wri[IDN]*wri[IVX]*wri[IVZ];
          fassign = true;
        }
      }

      //--- Step 3. If the flux has not been assigned at this point,
      //            we're in the intermediate state and now assign the flux
      if (!fassign){
        if (vxm >= 0.0){
          flxi[IDN] = dm*vxm;
          flxi[IVX] = dm*vxm*vxm + dm*iso_cs*iso_cs;
          flxi[IVY] = dm*vxm*wli[IVY];
          flxi[IVZ] = dm*vxm*wli[IVZ];
        }
        else{
          flxi[IDN] = dm*vxm;
          flxi[IVX] = dm*vxm*vxm + dm*iso_cs*iso_cs;
          flxi[IVY] = dm*vxm*wri[IVY];
          flxi[IVZ] = dm*vxm*wri[IVZ];
        }
      }
      flx(IDN,k,j,i) = flxi[IDN];
      flx(ivx,k,j,i) = flxi[IVX];
      flx(ivy,k,j,i) = flxi[IVY];
      flx(ivz,k,j,i) = flxi[IVZ];
    }
    return;
  }

  static void srder(double cs, double dm, double vl, double vr, double dmin, double dmax, double &y, double &dydx){
    y = (vr - vl) + cs*(std::log(dm/dmax) + (dm-dmin)/std::sqrt(dm*dmin));
    dydx = cs/dm*(1.0 + 0.5*(dm+dmin)/std::sqrt(dm*dmin));
    return;
  }

  static double rtsafe(double cs, double x1, double x2, double xacc, double vl, double vr, double dmin, double dmax){
    int j;
    double df,dx,dxold,f,fh,fl;
    double temp,xh,xl,rts;
    int maxit = 100;

    srder(cs,x1,vl,vr,dmin,dmax,fl,df);
    srder(cs,x2,vl,vr,dmin,dmax,fh,df);
    if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)){
      return 0.0;
    }
    if (fl == 0.0) return x1;
    if (fh == 0.0) return x2;
    if (fl < 0.0){
      xl = x1;
      xh = x2;
    }
    else{
      xh = x1;
      xl = x2;
    }
    rts = 0.5*(x1+x2);
    dxold = std::fabs(x2-x1);
    dx = dxold;
    srder(cs,rts,vl,vr,dmin,dmax,f,df);
    for (j=1;j<=maxit;j++){
      if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) ||
      (std::fabs(2.0*f) > std::fabs(dxold*df))){
        dxold = dx;
        dx = 0.5*(xh-xl);
        rts = xl+dx;
        if (xl == rts)
        return rts;
      }
      else{
        dxold = dx;
        dx = f/df;
        temp=rts;
        rts -= dx;
        if (temp == rts)
        return rts;
      }
      if (std::fabs(dx) < xacc)
      return rts;
      srder(cs,rts,vl,vr,dmin,dmax,f,df);
      if (f < 0.0)
      xl = rts;
      else
      xh = rts;
    }
    return rts;
  }
