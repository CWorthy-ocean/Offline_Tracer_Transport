use AllLocalesBarriers;
use Zarr;

use domains;
use dynamics;
use INPUTS;
use NetCDF_IO;
use tracers;


proc prepare_to_timestep(step : int) {

  // Load velocity fields for this timestep
    get_var(velfiles[step], 'u', u_n, D3_u);
    get_var(velfiles[step], 'v', v_n, D3_v);

  // Load viscosity for this timestep
    get_var(velfiles[step], 'Aks', kappa_v, D3);
    update_halos(kappa_v);

  calc_volumetric_fluxes(u_n, v_n, U_n, V_n, H_n);

  // Load thickness for future timestep (used for vertical remapping)
    update_thickness(zeta_np1, H_np1, H0, h, step+1);

}

proc prepare_next_timestep(step : int) {

  // Copy next thickness to current thickness
    forall (i,j,k) in D3.localSubdomain() {
      H_n[i,j,k] = H_np1[i,j,k];
    }
    update_halos(H_n);

  // Load the boundary data for the upcoming timestep
    var bryloc = bryfiles[step+1];
    for t in 1..num_tracers {
      set_bry(bryloc, "dye", tracer_n, D3.localSubdomain(), t);
    }

  allLocalesBarrier.barrier();
  update_halos(tracer_n);

}

proc update_thickness(ref zeta, ref H, ref H0, ref h, step : int) {

  // Read in SSH, update thicknesses for (n+1)
    get_var(velfiles[Nt_start], "zeta", zeta, D2);

  // From SM09, Eq. 2.13
    forall (i,j,k) in D3.localSubdomain() {
      H[i,j,k] = H0[i,j,k] * (1 + zeta[i,j] / h[i,j]);
    }

  allLocalesBarrier.barrier();
  update_halos(H);

}

proc update_halos(ref arr) {

  // Update the halos. Only need to issue the updateFluff
  // method from one Locale.

  allLocalesBarrier.barrier();
  if (here.id == 0) {
    arr.updateFluff();
  }
  allLocalesBarrier.barrier();

}
