use INPUTS;
use domains;
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
  var ktmp : [D3] real;

  var tracer_n : [D3] real;
  var tracer_tilde : [D3] real;
  var tracer_dagger : [D3] real;

  var mask_rho : [D2] real;
  var h : [D2] real;
  var H0 : [D3] real;
  var H_n : [D3] real;
  var H_np1 : [D3] real;
  var H_tilde : [D3] real;
  var H_dagger : [D3] real;

  var zeta_n : [D2] real;
  var zeta_np1 : [D2] real;

  var bry_sponge : [D3] real;
  var visc : [D3] real;

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
    coforall loc in Locales do on loc {
      update_halos(mask_rho);
      update_halos(h);
      update_halos(H0);
      update_halos(tracer_n);
    }
}
