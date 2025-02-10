use Math;

use INPUTS;

proc Cs(ref sigma: [?D] real) {

  var C = (1 - cosh(theta_s * sigma)) / (cosh(theta_s) - 1);
  var C2 = (exp(theta_b * C) - 1) / (1 - exp(-theta_b));

  return C2;
}

proc get_H0(ref h: [?D] real) {

  var k_w : [0..Nz] int;

  for idx in k_w.domain {
    k_w[idx] = idx;
  }

  var sigma_w = (k_w - Nz) / (1.0*Nz);

  var Cs_w = Cs(sigma_w);

  var z_w : [D.dim[0], D.dim[1], 0..Nz] real;
  var H0 : [D.dim[0], D.dim[1], 0..<Nz] real;

  forall (i,j,k) in z_w.domain {
    z_w[i,j,k] = h[i,j] * (hc*sigma_w[k] + h[i,j]*Cs_w[k]) / (hc + h[i,j]);
  }

  forall (i,j,k) in H0.domain {
    H0[i,j,k] = z_w[i,j,k+1] - z_w[i,j,k];
  }

  return H0;
}

