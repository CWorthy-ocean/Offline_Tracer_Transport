use AllLocalesBarriers;

use domains;
use dynamics;
use INPUTS;
use NetCDF_IO;
use tracers;


// This function will apply the piecewise parabolic method (PPM) as described in
// White and Adcroft (2008).

proc Polyfit(ref tracer_n, ref tracer_dagger, num_tracers) {

  if (num_tracers > 0) {

  var D3_loc = D3_int_3D.localSubdomain();
  forall (i,j) in {D3_loc.dim[0], D3_loc.dim[1]} {

  if (mask_rho.localAccess[i,j] == 1) {

    var interface_values : [1..num_tracers,0..Nz] real;

    interface_values = calc_interface_values(i, j, tracer_dagger, H_dagger, num_tracers);

    // Create the coefficients for the polynomial in each cell
    // Going to treat "left" as equivalent to "bottom", and "right" as equivalent to "top"

    var a0 : [1..num_tracers,0..<Nz] real;
    var a1 : [1..num_tracers,0..<Nz] real;
    var a2 : [1..num_tracers,0..<Nz] real;

    for (t,k) in {1..num_tracers,0..<Nz} {
      a0[t,k] = interface_values[t,k];
      a1[t,k] = 6*tracer_dagger.localAccess[i,j,t,k] - 4*interface_values[t,k] - 2*interface_values[t,k+1];
      a2[t,k] = 3*(interface_values[t,k] + interface_values[t,k+1] - 2*tracer_dagger.localAccess[i,j,t,k]);
    }

    // These are vector copies of each column.
      var H_orig : [0..<Nz] real;
      var H_new  : [0..<Nz] real;
      for kk in 0..<Nz {
        H_orig[kk] = H_dagger.localAccess[i,j,kk];
        H_new[kk]  = H_tmp2.localAccess[i,j,kk];
      }

      var total_tracer_orig : [1..num_tracers] real = 0;
      for (t,k) in {1..num_tracers,0..<Nz} {
        total_tracer_orig[t] += H_orig[k] * tracer_dagger.localAccess[i,j,t,k];
      }

    // Going to normalize the thicknesses to make the forthcoming loop logic and
    // tracer conservation a bit easier. This is just for bookkeeping here and will not
    // affect the actual thicknesses.
      var H_new_tot = + reduce H_new;
      var H_orig_tot = + reduce H_orig;

      H_orig = H_orig * H_new_tot / H_orig_tot;

    // The division and multiplication in the above step may not make the total column thicknesses
    // match exactly, so one more really small correction will be added to the top layer.
      H_orig_tot = + reduce H_orig;
      H_orig[Nz-1] += H_new_tot - H_orig_tot;

    // Also keep track of the original interface depths.
      var z_orig_interfaces : [0..Nz] real;
      for k in 0..<Nz {
        z_orig_interfaces[k+1] = z_orig_interfaces[k] + H_orig[k];
      }

    // For each interface in H_new, figure out in which cell of H_orig it lies
    var curr_src_index = 0;
    var tgt_index_start : [0..<Nz] int = 0;
    var tgt_index_end : [0..<Nz] int = 0;

    // Keep track of how far up we are in the water column
    var curr_tgt_height = H_new[0];
    var curr_src_height = H_orig[0];

    // Keep track of where in each H_orig cell the interfaces of H_new are
    var tgt_frac_start : [0..<Nz] real = 0;
    var tgt_frac_end : [0..<Nz] real = 0;

    for k_new in 0..<(Nz-1) {

      while (curr_tgt_height > curr_src_height) {
        curr_src_index += 1;
        curr_src_height += H_orig[curr_src_index];
      }

      tgt_index_end[k_new] = curr_src_index;
      tgt_index_start[k_new+1] = curr_src_index;

      tgt_frac_end[k_new] = (curr_tgt_height - z_orig_interfaces[curr_src_index]) / H_orig[curr_src_index];
      tgt_frac_start[k_new+1] = tgt_frac_end[k_new];

      curr_tgt_height += H_new[k_new+1];

    }
    tgt_index_end[Nz-1] = Nz-1;
    tgt_frac_end[Nz-1] = 1;

    // This holds the integrated_values
    var definite_integral : [1..num_tracers, 0..<Nz] real = 0;

    for (t,k) in {1..num_tracers,0..<(Nz-1)} {

      var tot_dist : real = 0;
      for idx in tgt_index_start[k]..tgt_index_end[k] {

        if (tgt_index_start[k] == tgt_index_end[k]) {
          definite_integral[t,k] = integrate(a0[t,idx], a1[t,idx], a2[t,idx], tgt_frac_start[k], tgt_frac_end[k]);
          tot_dist = tgt_frac_end[k] - tgt_frac_start[k];
        }
        else if ((idx == tgt_index_start[k]) && (idx < tgt_index_end[k])) {
          definite_integral[t,k] = integrate(a0[t,idx], a1[t,idx], a2[t,idx], tgt_frac_start[k], 1);
          tot_dist = 1 - tgt_frac_start[k];
        }
        else if ((idx < tgt_index_end[k]) && (idx > tgt_index_start[k])) {
          definite_integral[t,k] = definite_integral[t,k] + integrate(a0[t,idx], a1[t,idx], a2[t,idx], 0, 1);
          tot_dist = tot_dist + 1;
        }
        else {
          definite_integral[t,k] = definite_integral[t,k] + integrate(a0[t,idx], a1[t,idx], a2[t,idx], 0, tgt_frac_end[k]);
          tot_dist = tot_dist + tgt_frac_end[k];
        }

      }
      tracer_n.localAccess[i,j,t,k] = definite_integral[t,k] / tot_dist;
    }

    var total_tracer_new : [1..num_tracers] real = 0;
    for (t,k) in {1..num_tracers,0..<(Nz-1)} {
      total_tracer_new[t] += H_new[k] * tracer_n.localAccess[i,j,t,k];
    }
    // Can I just take a guess at this, instead of doing this whole calculation?
    for t in 1..num_tracers {
      tracer_n.localAccess[i,j,t,Nz-1] = (total_tracer_orig[t] - total_tracer_new[t]) / H_new[Nz-1];
    }

    total_tracer_new = 0;
    for (t,k) in {1..num_tracers,0..<Nz} {
      total_tracer_new[t] += H_new[k] * tracer_n.localAccess[i,j,t,k];
    }
    tracer_n.localAccess[i,j,1,Nz-1] += total_tracer_orig[1] - total_tracer_new[1];

    } // mask_rho

  } // end of forall

} // if num_tracers > 0

} // end of subroutine


proc integrate(a0, a1, a2, z0, z1) {

  var definite_integral = a0 * (z1 - z0) + 0.5*a1*(z1**2 - z0**2) + one_third*a2*(z1**3 - z0**3);

  return definite_integral;
}

proc calc_interface_values(i, j, ref arr, ref H, num_tracers) {

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                     Get tracer values at layer interfaces                               //
  //  Calculated with Implicit Fourth-order scheme using Thomas algorithm:  White and Adcroft, 2008, Eq. 46  //
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // First get the tracer values at the top and bottom boundaries using conservative
  // piecewise polynomial reconstruction: White and Adcroft, 2008, Eq. 7

    var Dp : domain(1) = {0..Nz};
    var DpDp : domain(2) = {0..Nz, 0..Nz};

  // Bottom boundary extrapolation
  // Calculated using a polynomial of order "ord" over the bottom (ord+1) cells

    var B : [1..num_tracers,0..ord] real;
    var M : [0..ord,0..ord] real;

    var h_b : real = 0;
    var h_t : real = H.localAccess[i,j,0];
    var iH = 1.0 / (h_t - h_b);

    var Ts_bot : [1..num_tracers] real;
    var Ts_top : [1..num_tracers] real;

    for k in 0..ord {
      for kk in 0..ord {
        var kkp1 = kk+1;
        M[k,kk] = (1.0 / kkp1)*iH*(h_t**(kkp1) - h_b**(kkp1));
      }

      h_b = h_b + H.localAccess[i,j,k];
      h_t = h_t + H.localAccess[i,j,k+1];
      iH = 1.0/(h_t - h_b);

      for t in 1..num_tracers {
        B[t,k] = arr.localAccess[i,j,t,k];
      }
    }

    for t in 1..num_tracers {
      Ts_bot[t] = gauss(M, B[t,..]);
    }

    // Top boundary extrapolation
    h_b = 0;
    h_t = H.localAccess[i,j,Nz-1];
    iH = 1.0/(h_t - h_b);
    for k in 0..ord {
      for kk in 0..ord {
        var kkp1 = kk+1;
        M[k,kk] = (1.0 / kkp1)*iH*(h_t**kkp1 - h_b**kkp1);
      }

      h_b = h_b + H.localAccess[i,j,Nz-1-k];
      h_t = h_t + H.localAccess[i,j,Nz-2-k];
      iH = 1.0/(h_t - h_b);

      for t in 1..num_tracers {
        B[t,k] = arr.localAccess[i,j,t,Nz-1-k];
      }
    }

    for t in 1..num_tracers {
      Ts_top[t] = gauss(M, B[t,..]);
    }

    /////////////////////////////

    var interface_values : [1..num_tracers,0..Nz] real;
    for t in 1..num_tracers {
      interface_values[t,..] = thomas_PPM(t,i,j,Ts_bot, Ts_top, arr, H);
    }

  return interface_values;

}

proc gauss(ref M, ref b) {

    // This routine uses Gaussian Elimination to reduce the matrix M to row-echelon form (lower triangular).
    // It only returns the first element of the solution vector.

    // Starting from the last row, we're going to work upwards to
    // make M a lower triangular matrix
    for i in 1..ord by -1 {

      // Iterate over all rows above i
      for j in 0..(i-1) {
        var ratio = M[j,i] / M[i,i];

        // Iterate over all columns up to i.
        // This operation basically is multiplying the ith row by "ratio" and adding it to the jth row
        for k in 0..i {
          M[j,k] -= ratio*M[i,k];
        }
        b[j] -= ratio*b[i];
      }
    }

    return (b[0] / M[0,0]);

}

proc thomas_PPM(t, i, j, Ts_bot, Ts_top, ref arr, ref H) {

  //////////////////////////////////////////////////////////////////////////////////////
  //                   Get tracer values at layer interfaces                          //
  //  Calculated with Implicit Fourth-order scheme:  White and Adcroft, 2008, Eq. 46  //
  //////////////////////////////////////////////////////////////////////////////////////

    var n = Nz+1;
    var Dp : domain(1) = {1..n};

    var a : [Dp] real;
    var b : [Dp] real;
    var c : [Dp] real;
    var d : [Dp] real;

    var cp : [Dp] real;
    var dp : [Dp] real;
    var x  : [Dp] real;

    b[1] = 1.0;
    b[n] = 1.0;
    d[1] = Ts_bot[t];
    d[n] = Ts_top[t];

    for k in 1..(n-2) {
      var h0 = H.localAccess[i,j,k-1];
      var h1 = H.localAccess[i,j,k];

      var alpha = (h1**2) / ((h0 + h1)**2);
      var beta = (h0**2) / ((h0 + h1)**2);
      var d1 = 2*(h1**2)*(h1**2 + 2*h0**2 + 3*h0*h1) / ((h0+h1)**4);
      var d2 = 2*(h0**2)*(h0**2 + 2*h1**2 + 3*h0*h1) / ((h0+h1)**4);

      a[k+1] = alpha;
      b[k+1] = 1.0;
      c[k+1] = beta;

      d[k+1] = d1*arr.localAccess[i,j,t,k-1] + d2*arr.localAccess[i,j,t,k];
    }

    cp[1] = c[1] / b[1];
    dp[1] = d[1] / b[1];

    for k in 2..(n-1) {
      cp[k] = c[k] / (b[k] - a[k]*cp[k-1]);
      dp[k] = (d[k] - a[k]*dp[k-1]) / (b[k] - a[k]*cp[k-1]);
    }

    dp[n] = (d[n] - a[n]*dp[n-1]) / (b[n] - a[n]*cp[n-1]);

    x[n] = dp[n];
    for k in 1..(n-1) by -1 {
      x[k] = dp[k] - cp[k]*x[k+1];
    }

    return x;
}

