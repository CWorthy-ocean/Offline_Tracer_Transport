
// The concentration of a tracer \psi within a grid cell is represented as a smoothly-varying
// polynomial that can be expressed as (Bott, 1989; Easter, 1993)
//
// \psi_{j,l}(x) = \sum_{k=0}^l a_{j,k}*(x-x_j)/\Delta_x.  (1)
//
// We know the average value of this polynomial within a cell j (i.e. the integral of (1) divided by
// the cell width) is equal to \psi_j; therefore, we can solve for the coefficients by solving a matrix
// equation (the matrix is n x n, where n=l+1). (White and Adcroft, 2008)
//
// We then evaluate the polynomial at the location of the cell edge to get the interpolated value of \psi.
// Bott (1989) gives coefficients of the approximating polynomial for uniform grids, and White and
// Adcroft give expressions for the actual edge values for uniform grids (note that the expressions in these
// papers are not equal!)
//
// On the C-grid the flux is evaluated by multiplying the interpolated tracer concentration by the volume flux.

// The "constant grid flux form" (CGF) method uses slope ratios, defined as
//
// r_{j+1/2} = (psi_j - \psi_{j-1}) / (\psi_{j+1} - \psi_j)
//
// The flux limiter \phi works to blend the contributions of a low- and high-order flux, so that
//
// F_{j+1/2} = F_L_{j+1/2} + \phi_{j+1/2}(r) * (F_H_{j+1/2} - F^L_{j+1/2})
//
// We will use third-order estimates here for the high-order flux (Table 1 of Easter, 1993). For evenly spaced grids,
// the coefficients of the second-order approximating polynomial (l=2, third-order accurate) in \psi are
// a0 = (1/24) * (-\psi_{j+1} + 26*\psi_j - \psi_{j-1})
// a1 = 0.5 * (\psi_{j+1} - \psi_{j-1})
// a2 = 0.5 * (\psi_{j+1} - 2*\psi_j + \psi_{j-1}).
//
// For the low-order flux we will use a first-order approximation, which assumes a constant value in the cell.
//
// We will use the Superbee limiter, defined so that
//
// \phi(r) = max[0, min(1,2r), min(r,2)].
//

use AllLocalesBarriers;

use domains;
use dynamics;
use INPUTS;
use tracers;

proc calc_horizontal_fluxes_U(ref U, ref tmp_U, ref arr, const t) {

  /////////////////////////////////////////
  //              U-fluxes               //
  /////////////////////////////////////////

  forall (i,j,k) in D3_int_u.localSubdomain() {

    // Slope ratio
      var r = (arr.localAccess[i,j,t,k] - arr.localAccess[i-1,j,t,k]) / (arr.localAccess[i+1,j,t,k] - arr.localAccess[i,j,t,k] + eps);

    // Superbee limiter
      var sb = max(0, min(1, 2*r), min(r,2));

    // 3rd order, upstream-biased parabolic interpolation: SM05, after 4.13
      var tmp_h = 0.5*(arr.localAccess[i,j,t,k] + arr.localAccess[i+1,j,t,k]) * U.localAccess[i,j,k] - mask_rho.localAccess[i+2,j]*
                           (one_sixth * max(U.localAccess[i,j,k], 0.0) * (arr.localAccess[i+1,j,t,k] - arr.localAccess[i,j,t,k])
                           +  one_sixth * min(U.localAccess[i,j,k], 0.0) * (arr.localAccess[i+2,j,t,k] - 2*arr.localAccess[i+1,j,t,k] + arr.localAccess[i,j,t,k]));
    // 1st-order interpolation (constant in the cell)
      var tmp_l = max(U.localAccess[i,j,k], 0.0) * arr.localAccess[i,j,t,k] + min(U.localAccess[i,j,k], 0.0) * arr.localAccess[i+1,j,t,k];

    // Limited flux
      tmp_U.localAccess[i,j,k] = tmp_l + sb * (tmp_h - tmp_l);

  }

  forall (i,j,k) in D3_edge_u_w.localSubdomain() {

    // 1st-order interpolation (constant in the cell)
      var tmp_l = max(U.localAccess[i,j,k], 0.0) * arr.localAccess[t,i,j,k] + min(U.localAccess[i,j,k], 0.0) * arr.localAccess[t,i+1,j,k];

    // Limited flux
      tmp_U.localAccess[i,j,k] = tmp_l;

  }

  forall (i,j,k) in D3_edge_u_e.localSubdomain() {

    // 1st-order interpolation (constant in the cell)
      var tmp_l = max(U.localAccess[i,j,k], 0.0) * arr.localAccess[t,i,j,k] + min(U.localAccess[i,j,k], 0.0) * arr.localAccess[t,i+1,j,k];

    // Limited flux
      tmp_U.localAccess[i,j,k] = tmp_l;

  }

  update_halos(tmp_U);

}

proc calc_horizontal_fluxes_V(ref V, ref tmp_V, ref arr, const t) {

  /////////////////////////////////////////
  //              V-fluxes               //
  /////////////////////////////////////////

  forall (i,j,k) in D3_int_v.localSubdomain() {

    // Slope ratio
      var r = (arr.localAccess[i,j,t,k] - arr.localAccess[i,j-1,t,k]) / (arr.localAccess[i,j+1,t,k] - arr.localAccess[i,j,t,k] + eps);

    // Superbee limiter
      var sb = max(0, min(1, 2*r), min(r,2));

    // 3rd order, upstream-biased parabolic interpolation: SM05, after 4.13
      var tmp_h = 0.5*(arr.localAccess[i,j,t,k] + arr.localAccess[i,j+1,t,k]) * V.localAccess[i,j,k] - mask_rho.localAccess[i,j+2]*
                           (one_sixth * max(V.localAccess[i,j,k], 0.0) * (arr.localAccess[i,j+1,t,k] - arr.localAccess[i,j,t,k])
                           +  one_sixth * min(V.localAccess[i,j,k], 0.0) * (arr.localAccess[i,j+2,t,k] - 2*arr.localAccess[i,j+1,t,k] + arr.localAccess[i,j,t,k]));
    // 1st-order interpolation (constant in the cell)
      var tmp_l = max(V.localAccess[i,j,k], 0.0) * arr.localAccess[i,j,t,k] + min(V.localAccess[i,j,k], 0.0) * arr.localAccess[i,j+1,t,k];

    // Limited flux
      tmp_V.localAccess[i,j,k] = tmp_l + sb * (tmp_h - tmp_l);

  }

  forall (i,j,k) in D3_edge_v_s.localSubdomain() {

    // 1st-order interpolation (constant in the cell)
      var tmp_l = max(V.localAccess[i,j,k], 0.0) * arr.localAccess[t,i,j,k] + min(V.localAccess[i,j,k], 0.0) * arr.localAccess[t,i,j+1,k];

    // Limited flux
      tmp_V.localAccess[i,j,k] = tmp_l;

  }

  forall (i,j,k) in D3_edge_v_n.localSubdomain() {

    // 1st-order interpolation (constant in the cell)
      var tmp_l = max(V.localAccess[i,j,k], 0.0) * arr.localAccess[t,i,j,k] + min(V.localAccess[i,j,k], 0.0) * arr.localAccess[t,i,j+1,k];

    // Limited flux
      tmp_V.localAccess[i,j,k] = tmp_l;

  }

  update_halos(tmp_V);

}

