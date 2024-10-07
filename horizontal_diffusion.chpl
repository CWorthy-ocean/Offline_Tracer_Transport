use AllLocalesBarriers;

use INPUTS;
use domains;
//use sigma_coordinate;
use tracers;

/*
record Diffusion {

  var u0, u1, u2 = 1..0;
  var v0, v1, v2 = 1..0;

  var U_n : [u0, u1, u2] real;
  var U_np1h : [u0, u1, u2] real;
  var U_np1 : [u0, u1, u2] real;
  var tmp_U : [u0, u1, u2] real;

  var V_n : [v0, v1, v2] real;
  var V_np1h : [v0, v1, v2] real;
  var V_np1 : [v0, v1, v2] real;
  var tmp_V : [v0, v1, v2] real;

  proc init(arg: Domains) {
    this.u0 = arg.u_3D.dim[0];
    this.u1 = arg.u_3D.dim[1];
    this.u2 = arg.u_3D.dim[2];

    this.v0 = arg.v_3D.dim[0];
    this.v1 = arg.v_3D.dim[1];
    this.v2 = arg.v_3D.dim[2];

  }

}
*/

// Apply Laplacian diffusion in the sponge layer near the domain boundary

proc calc_diffusive_fluxes(ref U, ref V, ref arr, ref H) {

  /////////////////////////////////////////
  //             Zonal fluxes            //
  /////////////////////////////////////////

  forall (k,j,i) in D3_u.localSubdomain() {

    U[k,j,i] = 0.5*(sponge[k,j,i] + sponge[k,j,i+1]) * (arr[k,j,i+1] - arr[k,j,i])
                         * dy * 0.5 * (H[k,j,i] + H[k,j,i+1]) / dx;
  }

  /////////////////////////////////////////
  //         Meridional fluxes           //
  /////////////////////////////////////////

  forall (k,j,i) in D3_v.localSubdomain() {

    V[k,j,i] = 0.5*(sponge[k,j,i] + sponge[k,j+1,i]) * (arr[k,j+1,i] - arr[k,j,i])
                         * dx * 0.5 * (H[k,j,i] + H[k,j+1,i]) / dy;
  }

  allLocalesBarrier.barrier();

}

// Initialize the boundary sponge

proc initialize_sponge() {

  coforall loc in Locales with (ref H_n) do on loc {
    forall (k,j,i) in D3.localSubdomain() {
    //forall (k,j,i) in D.rho_3D {
      var dist_from_x = min(i, Nx - 1 - i);
      var dist_from_y = min(j, Ny - 1 - j);
      var min_dist = max(0.0, min(dist_from_x, dist_from_y) - 1);
      var amp = max(0.0, (sponge_width - min_dist) / sponge_width);
      sponge[k,j,i] = v_sponge * amp;
    }
  }
}


