use AllLocalesBarriers;

use domains;
use dynamics;
use INPUTS;
use tracers;
use updates;

// Apply Laplacian diffusion in the sponge layer near the domain boundary

proc calc_diffusive_fluxes_U(ref U, ref arr, ref H, const t) {

  /////////////////////////////////////////
  //             Zonal fluxes            //
  /////////////////////////////////////////

  forall (i,j,k) in D3_u.localSubdomain() {

    U[i,j,k] = 0.5*(visc[i,j,k] + visc[i+1,j,k]) * (arr[i+1,j,t,k] - arr[i,j,t,k])
                         * dy * 0.5 * (H[i,j,k] + H[i+1,j,k]) / dx;
  }

  update_halos(U);

}

proc calc_diffusive_fluxes_V(ref V, ref arr, ref H, const t) {

  /////////////////////////////////////////
  //         Meridional fluxes           //
  /////////////////////////////////////////

  forall (i,j,k) in D3_v.localSubdomain() {

    V[i,j,k] = 0.5*(visc[i,j,k] + visc[i,j+1,k]) * (arr[i,j+1,t,k] - arr[i,j,t,k])
                         * dx * 0.5 * (H[i,j,k] + H[i,j+1,k]) / dy;
  }

  update_halos(V);

}

// Initialize the boundary sponge

proc initialize_sponge() {

  forall (i,j,k) in D3.localSubdomain() {
    var dist_from_x = min(Nx - 1 - i,i);
    var dist_from_y = min(Ny - 1 - j,j);
    var min_dist = max(0.0, min(dist_from_x, dist_from_y) - 1);
    var amp = max(0.0, (sponge_width - min_dist) / sponge_width);
    visc[i,j,k] = v_sponge * amp;
  }

  update_halos(visc);

}


