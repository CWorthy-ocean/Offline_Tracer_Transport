use AllLocalesBarriers;

use domains;
use INPUTS;
use tracers;


proc calc_vertical_diffusion(ref arr, ref H) {

  // This will update tracer_dagger with an implicit timestep for the vertical diffusion

  var D3_loc = D3_tr.localSubdomain();
  forall (t,i,j) in {D3_loc.dim[0], D3_loc.dim[1], D3_loc.dim[2]} {

    var tmmp = thomas_diff(t,i,j, arr, H);

    for kk in 0..<Nz {
      arr[t,i,j,kk] = tmmp[kk+1];
    }
  }

  allLocalesBarrier.barrier();

}

proc thomas_diff(t, i, j, ref arr, ref H) {

  var n = Nz;
  var Dp : domain(1) = {1..n};

  var a : [Dp] real;
  var b : [Dp] real;
  var c : [Dp] real;
  var d : [Dp] real;

  var cp : [Dp] real;
  var dp : [Dp] real;
  var x  : [Dp] real;

  a[1] = 0;
  c[1] = -(2*dt*kappa_v[i,j,1]) / (H[i,j,0] * (H[i,j,1]+H[i,j,0]));
  b[1] = 1 - a[1] - c[1];
  d[1] = arr[t,i,j,0];

  a[n] = -(2*dt*kappa_v[i,j,n-1]) / (H[i,j,n-1] * (H[i,j,n-1]+H[i,j,n-2]));
  c[n] = 0;
  b[n] = 1 - a[n] - c[n];
  d[n] = arr[t,i,j,n-1];

  for k in 2..<n {
    a[k] = -(2*dt*kappa_v[i,j,k-1]) / (H[i,j,k-1] * (H[i,j,k-1]+H[i,j,k-2]));
    c[k] = -(2*dt*kappa_v[i,j,k]) / (H[i,j,k-1] * (H[i,j,k]+H[i,j,k-1]));
    b[k] = 1 - a[k] - c[k];

    d[k] = arr[t,i,j,k-1];
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
