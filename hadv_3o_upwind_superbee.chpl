
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

proc calc_horizontal_fluxes_U(ref U, ref tmp_U, ref arr) {

  /////////////////////////////////////////
  //              U-fluxes               //
  /////////////////////////////////////////

  forall (k,j,i) in D3_u.localSubdomain() {

    // Slope ratio
      var r = (arr[k,j,i] - arr[k,j,i-1]) / (arr[k,j,i+1] - arr[k,j,i] + eps);

    // Superbee limiter
      var sb = max(0, min(1, 2*r), min(r,2));

    // 3rd order, upstream-biased parabolic interpolation: SM05, after 4.13
      var tmp_h = 0.5*(arr[k,j,i] + arr[k,j,i+1]) * U[k,j,i] - mask_rho[j,i+2]*
                           (one_sixth * max(U[k,j,i], 0.0) * (arr[k,j,i+1] - arr[k,j,i])
                           +  one_sixth * min(U[k,j,i], 0.0) * (arr[k,j,i+2] - 2*arr[k,j,i+1] + arr[k,j,i]));
    // 1st-order interpolation (constant in the cell)
      var tmp_l = max(U[k,j,i], 0.0) * arr[k,j,i] + min(U[k,j,i], 0.0) * arr[k,j,i+1];

    // Limited flux
    tmp_U[k,j,i] = tmp_l + sb * (tmp_h - tmp_l);
  }

  update_halos(tmp_U);

}

proc calc_horizontal_fluxes_V(ref V, ref tmp_V, ref arr) {

  /////////////////////////////////////////
  //              V-fluxes               //
  /////////////////////////////////////////

  forall (k,j,i) in D3_v.localSubdomain() {

    // Slope ratio
      var r = (arr[k,j,i] - arr[k,j-1,i]) / (arr[k,j+1,i] - arr[k,j,i] + eps);

    // Superbee limiter
      var sb = max(0, min(1, 2*r), min(r,2));

    // 3rd order, upstream-biased parabolic interpolation: SM05, after 4.13
      var tmp_h = 0.5*(arr[k,j,i] + arr[k,j+1,i]) * V[k,j,i] - mask_rho[j+2,i]*
                           (one_sixth * max(V[k,j,i], 0.0) * (arr[k,j+1,i] - arr[k,j,i])
                           +  one_sixth * min(V[k,j,i], 0.0) * (arr[k,j+2,i] - 2*arr[k,j+1,i] + arr[k,j,i]));
    // 1st-order interpolation (constant in the cell)
      var tmp_l = max(V[k,j,i], 0.0) * arr[k,j,i] + min(V[k,j,i], 0.0) * arr[k,j+1,i];

    // Limited flux
    tmp_V[k,j,i] = tmp_l + sb * (tmp_h - tmp_l);
  }

  update_halos(tmp_V);

}

