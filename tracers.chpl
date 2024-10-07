use INPUTS;
use domains;
//use params;
//use dynamics;
use sigma_coordinate;
use updates;
use NetCDF_IO;

use Zarr;
use StencilDist;
use AllLocalesBarriers;
use FileSystem;

// Notation here will be consistent with the description in
// https://adcroft.github.io/assets/pdf/ALE_workshop_NCWCP_2016.pdf

// For RK3
  var k1 : [D3] real;
  var k2 : [D3] real;
  var k3 : [D3] real;
  var ktmp : [D3] real;

  var tracer_n : [D3] real;
//  var tracer_np1h : [D3] real;
//  var tracer_np1 : [D3] real;
  var tracer_dagger : [D3] real;

  var mask_rho : [D2] real;
  var h : [D2] real;
  var H0 : [D3] real;
  var H_n : [D3] real;
  var H_np1h : [D3] real;
  var H_np1 : [D3] real;
  var H_dagger : [D3] real;

  var zeta_n : [D2] real;
  var zeta_np1h : [D2] real;
  var zeta_np1 : [D2] real;

  var sponge : [D3] real;

  var kappa_v : [D3] real;

  var velfiles = glob(velocity_files);
  var bryfiles = glob(boundary_files);

proc initialize_tr() {

    mask_rho = readZarrArray(maskfile, real(64), 2, targetLocales=myTargetLocales2D);
    h = readZarrArray(hfile, real(64), 2, targetLocales=myTargetLocales2D);

    coforall loc in Locales do on loc {
      var D2_loc = D2.localSubdomain();
      var D3_loc = D3.localSubdomain();
      H0[D3_loc] = get_H0(h[D2_loc]);
    }

    kappa_v = readZarrArray(velfiles[Nt_start] + '/Akt/', real(32), 3, targetLocales=myTargetLocales3D);

    // Initialize tracer fields
    tracer_n = readZarrArray(velfiles[Nt_start] + '/temp/', real(32), 3, targetLocales=myTargetLocales3D);

    // Initialize zeta and thicknesses
      update_thickness(zeta_n, H_n, H0, h, Nt_start);

    // Update the halos
      update_halos(mask_rho);
      update_halos(h);
      update_halos(H0);
      update_halos(tracer_n);
//        mask_rho.updateFluff();
//        h.updateFluff();
//        H0.updateFluff();
//        tracer_n.updateFluff();

//    // Initialize thickness for next timestep
//      update_thickness(zeta_np1, H_np1, H0, h, Nt_start+1);
}

/*
proc update_halos(ref arr) {

    // Update the halos. Only need to issue the updateFluff
    // method from one Locale.

    allLocalesBarrier.barrier();
    if (here.id == 0) {
      arr.updateFluff();
    }
    allLocalesBarrier.barrier();

}

proc update_thickness(ref zeta, ref H, ref H0, ref h, step : int) {

  // Read in SSH, update thicknesses for (n+1)
//    zeta[D.rho_2D] = get_var(P.velfiles[step], "zeta", D.rho_2D);
//    zeta[D.rho_2D] = readZarrArrayLocal(P.velfiles[step] + (here.id : string) + '/zeta/', real(32), 2);
    zeta = readZarrArray(velfiles[step] + '/zeta/', real(32), 2, targetLocales=myTargetLocales2D);

  // From SM09, Eq. 2.13
  coforall loc in Locales do on loc {
    forall (k,j,i) in D3.localSubdomain() {
      H[k,j,i] = H0[k,j,i] * (1 + zeta[j,i] / h[j,i]);
    }
  }

//  allLocalesBarrier.barrier();
//  if (here.id == 0) {
    H.updateFluff();
//  }
//  allLocalesBarrier.barrier();

}

*/

proc calc_half_step_tr() {

    forall (j,i) in D2.localSubdomain() {
      zeta_np1h[j,i] = 0.5 * (zeta_n[j,i] + zeta_np1[j,i]);
    }

  // From SM09, Eq. 2.13
    forall (k,j,i) in D3.localSubdomain() {
      H_np1h[k,j,i] = H0[k,j,i] * (1 + zeta_np1h[j,i] / h[j,i]);
    }

    update_halos(H_np1h);

/*
    allLocalesBarrier.barrier();
    if (here.id == 0) {
      H_np1h.updateFluff();
    }
    allLocalesBarrier.barrier();
*/
}

