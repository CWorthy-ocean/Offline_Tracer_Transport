use INPUTS;
use domains;
use sigma_coordinate;
use updates;
use NetCDF_IO;
use Marbl;

use StencilDist;
use AllLocalesBarriers;
use FileSystem;

// Notation here will be consistent with the description in
// https://adcroft.github.io/assets/pdf/ALE_workshop_NCWCP_2016.pdf

// For RK3
  var ktmp : [D3] real;
  var ktmp_tr : [D3] real;

  var tracers_ts_n : [D3_ts] real;
  var tracers_ts_tilde : [D3_ts] real;
  var tracers_ts_dagger : [D3_ts] real;

  var tracers_marbl_n : [D3_marbl] real;
  var tracers_marbl_tilde : [D3_marbl] real;
  var tracers_marbl_dagger : [D3_marbl] real;

  var tracers_other_n : [D3_other] real;
  var tracers_other_tilde : [D3_other] real;
  var tracers_other_dagger : [D3_other] real;

  var mask_rho : [D2] real;
  var h : [D2] real;
  var H0 : [D3] real;
  var H_n : [D3] real;
  var H_np1 : [D3] real;
  var H_tilde : [D3] real;
  var H_dagger : [D3] real;

  var zeta_n : [D2] real;
  var zeta_np1 : [D2] real;

  var visc : [D3] real;

  var kappa_v : [D3] real;

  var velfiles = glob(velocity_files);
  var bryfiles = glob(boundary_files);

proc initialize_tr() {

  var D3_loc = D3.localSubdomain();
  var D2_loc = D2.localSubdomain();

  get_var(maskfile, "mask_rho", mask_rho, D2);
  get_var(hfile, "h", h, D2);

  H0[D3_loc] = get_H0(h[D2_loc]);

//  get_var(velfiles[Nt_start], 'Akt', kappa_v, D3);

  if (restart == 1) {
//    for t in 1..num_ts_tracers {
//      get_var(restart_file, ts_namelist[t], ktmp, D3);
//      tracers_ts_n[D3_loc.dim[0], D3_loc.dim[1], t, D3_loc.dim[2]] = ktmp[D3_loc];
//    }
//    for t in 1..num_marbl_tracers {
//      get_var(restart_file, marbl_namelist[t], ktmp, D3);
//      tracers_marbl_n[D3_loc.dim[0], D3_loc.dim[1], t, D3_loc.dim[2]] = ktmp[D3_loc];
//    }
    for t in 1..num_other_tracers {
      get_var(restart_file, other_namelist[t], ktmp, D3);
      tracers_other_n[D3_loc.dim[0], D3_loc.dim[1], t, D3_loc.dim[2]] = ktmp[D3_loc];
    }
  }
  else {
//    for t in 1..num_ts_tracers {
//      get_var(velfiles[Nt_start], ts_namelist[t], ktmp, D3);
//      tracers_ts_n[D3_loc.dim[0], D3_loc.dim[1], t, D3_loc.dim[2]] = ktmp[D3_loc];
//    }
//    for t in 1..num_marbl_tracers {
//      get_var(velfiles[Nt_start], marbl_namelist[t], ktmp, D3);
//      tracers_marbl_n[D3_loc.dim[0], D3_loc.dim[1], t, D3_loc.dim[2]] = ktmp[D3_loc];
//    }
    for t in 1..num_other_tracers {
      get_var(velfiles[Nt_start], other_namelist[t], ktmp, D3);
      tracers_other_n[D3_loc.dim[0], D3_loc.dim[1], t, D3_loc.dim[2]] = ktmp[D3_loc];
    }
  }

  // Initialize zeta and thicknesses
    update_thickness(zeta_n, H_n, H0, h, Nt_start);

  // Update the halos
    update_halos(mask_rho);
    update_halos(h);
    update_halos(H0);
    update_halos(kappa_v);
    update_halos(tracers_ts_n);
    update_halos(tracers_marbl_n);
    update_halos(tracers_other_n);
    
}
