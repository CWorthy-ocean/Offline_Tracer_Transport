use AllLocalesBarriers;
use Zarr;

use domains;
use dynamics;
use INPUTS;
use NetCDF_IO;
use tracers;


proc prepare_to_timestep(step : int) {

  // Load velocity fields for this timestep
    u_n = readZarrArray(velfiles[step] + '/u/', real(32), 3, targetLocales=myTargetLocales3D);
    v_n = readZarrArray(velfiles[step] + '/v/', real(32), 3, targetLocales=myTargetLocales3D);

    coforall loc in Locales do on loc {
      calc_volumetric_fluxes(u_n, v_n, U_n, V_n, H_n);
    }

  // Load thickness for future timestep (used for vertical remapping)
    update_thickness(zeta_np1, H_np1, H0, h, step+1);

}

proc prepare_next_timestep(step : int) {

  // Copy next thickness to current thickness
    forall (k,j,i) in D3.localSubdomain() {
      H_n[k,j,i] = H_np1[k,j,i];
    }
    update_halos(H_n);

  // Load the boundary data for the upcoming timestep
    var bryloc = bryfiles[step+1];
    set_bry(bryloc, "temp", tracer_n, D3.localSubdomain());

    allLocalesBarrier.barrier();
    update_halos(tracer_n);

}

proc update_thickness(ref zeta, ref H, ref H0, ref h, step : int) {

  // Read in SSH, update thicknesses for (n+1)
    zeta = readZarrArray(velfiles[step] + '/zeta/', real(32), 2, targetLocales=myTargetLocales2D);

  // From SM09, Eq. 2.13
  coforall loc in Locales do on loc {
    forall (k,j,i) in D3.localSubdomain() {
      H[k,j,i] = H0[k,j,i] * (1 + zeta[j,i] / h[j,i]);
    }
    allLocalesBarrier.barrier();
    update_halos(H);
  }

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
