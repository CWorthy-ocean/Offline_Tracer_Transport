use AllLocalesBarriers;

use domains;
use dynamics;
use INPUTS;
use tracers;


// Apply Laplacian diffusion in the sponge layer near the domain boundary

proc calc_diffusive_fluxes_U(ref U, ref arr, ref H) {

  /////////////////////////////////////////
  //             Zonal fluxes            //
  /////////////////////////////////////////

  forall (k,j,i) in D3_u.localSubdomain() {

    U[k,j,i] = 0.5*(visc[k,j,i] + visc[k,j,i+1]) * (arr[k,j,i+1] - arr[k,j,i])
                         * dy * 0.5 * (H[k,j,i] + H[k,j,i+1]) / dx;
  }

  update_halos(U);

}

proc calc_diffusive_fluxes_V(ref V, ref arr, ref H) {

  /////////////////////////////////////////
  //         Meridional fluxes           //
  /////////////////////////////////////////

  forall (k,j,i) in D3_v.localSubdomain() {

    V[k,j,i] = 0.5*(visc[k,j,i] + visc[k,j+1,i]) * (arr[k,j+1,i] - arr[k,j,i])
                         * dx * 0.5 * (H[k,j,i] + H[k,j+1,i]) / dy;
  }

  update_halos(V);

}

// Initialize the boundary sponge

proc initialize_sponge() {

  coforall loc in Locales with (ref H_n) do on loc {
    forall (k,j,i) in D3.localSubdomain() {
      var dist_from_x = min(i, Nx - 1 - i);
      var dist_from_y = min(j, Ny - 1 - j);
      var min_dist = max(0.0, min(dist_from_x, dist_from_y) - 1);
      var amp = max(0.0, (sponge_width - min_dist) / sponge_width);
      bry_sponge[k,j,i] = v_sponge * amp;
    }
  }
}


