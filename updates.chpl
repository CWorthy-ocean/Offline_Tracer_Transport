use AllLocalesBarriers;
//use Zarr;

use domains;
use dynamics;
use INPUTS;
use NetCDF_IO;
use tracers;


proc prepare_to_timestep(step : int) {

  // Load velocity fields for this timestep
    get_var(velfiles[step], 'u', u_n, D3_u);
    get_var(velfiles[step], 'v', v_n, D3_v);

  calc_volumetric_fluxes(u_n, v_n, U_n, V_n, H_n);

  // Load thickness for future timestep (used for vertical remapping)
    update_thickness(zeta_np1, H_np1, H0, h, step+1);

}

proc prepare_next_timestep(step : int) {

  // Copy next thickness to current thickness
    forall (i,j,k) in D3.localSubdomain() {
      H_n.localAccess[i,j,k] = H_np1.localAccess[i,j,k];
    }
    update_halos(H_n);

  // Load the boundary data for the upcoming timestep
    var bryloc = bryfiles[step+1];
    for t in 1..num_ts_tracers {
      set_bry(bryloc, ts_namelist[t-1], tracers_ts_n, D3.localSubdomain(), t);
    }
    for t in 1..num_marbl_tracers {
      set_bry(bryloc, marbl_namelist[t-1], tracers_marbl_n, D3.localSubdomain(), t);
    }
    for t in 1..num_other_tracers {
      set_bry(bryloc, other_namelist[t-1], tracers_other_n, D3.localSubdomain(), t);
    }

  allLocalesBarrier.barrier();
  update_halos(tracers_ts_n);
  update_halos(tracers_marbl_n);
  update_halos(tracers_other_n);

}

proc update_thickness(ref zeta, ref H, ref H0, ref h, step : int) {

  // Read in SSH, update thicknesses for (n+1)
    get_var(velfiles[Nt_start], "zeta", zeta, D2);

  // From SM09, Eq. 2.13
    forall (i,j,k) in D3.localSubdomain() {
      H.localAccess[i,j,k] = H0.localAccess[i,j,k] * (1 + zeta.localAccess[i,j] / h.localAccess[i,j]);
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
