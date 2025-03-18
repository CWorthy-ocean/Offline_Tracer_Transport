use AllLocalesBarriers;
//use Zarr;

use domains;
use dynamics;
use INPUTS;
use NetCDF_IO;
use tracers;
use forcings;
use IO;

proc prepare_to_timestep(step : int) {

  var curr_idx = (step / read_freq) : int;

  if (((step % read_freq) == 0) || restart) {
    restart = false;

  // Load velocity fields for current time
    get_var(velfiles[curr_idx], 'u', u_n, D3_u);
    get_var(velfiles[curr_idx], 'v', v_n, D3_v);
  // Load velocity fields for future time
    get_var(velfiles[curr_idx+1], 'u', u_np1, D3_u);
    get_var(velfiles[curr_idx+1], 'v', v_np1, D3_v);

  // Load SSH field for current time
    get_var(velfiles[curr_idx], "zeta", zeta_n, D2);
  // Load SSH field for future time
    get_var(velfiles[curr_idx+1], "zeta", zeta_np1, D2);

  // Load vertical diffusivity for current time
    get_var(velfiles[curr_idx], 'Aks', kappa_v_n, D3);
  // Load vertical diffusivity for future time
    get_var(velfiles[curr_idx+1], 'Aks', kappa_v_np1, D3);
  }


  // Interpolate to obtain velocity, SSH, and diffusivity fields for this timestep
  var frac = (step % read_freq) / read_freq;
  var frac2 = frac + (1/read_freq);

  forall (i,j,k) in D3_u.localSubdomain() {
    u_tmp.localAccess[i,j,k] = (1-frac) * u_n.localAccess[i,j,k] + frac * u_np1.localAccess[i,j,k];
  }
  forall (i,j,k) in D3_v.localSubdomain() {
    v_tmp.localAccess[i,j,k] = (1-frac) * v_n.localAccess[i,j,k] + frac * v_np1.localAccess[i,j,k];
  }
  forall (i,j) in D2.localSubdomain() {
    zeta_tmp1.localAccess[i,j] = (1-frac) * zeta_n.localAccess[i,j] + frac * zeta_np1.localAccess[i,j];
    zeta_tmp2.localAccess[i,j] = (1-frac2) * zeta_n.localAccess[i,j] + frac2 * zeta_np1.localAccess[i,j];
  }
  forall (i,j,k) in D3.localSubdomain() {
    kappa_v_tmp.localAccess[i,j,k] = (1-frac) * kappa_v_n.localAccess[i,j,k] + frac * kappa_v_np1.localAccess[i,j,k];
  }

  // From SM09, Eq. 2.13
    forall (i,j,k) in D3.localSubdomain() {
      H_tmp1.localAccess[i,j,k] = H0.localAccess[i,j,k] * (1 + zeta_tmp1.localAccess[i,j] / h.localAccess[i,j]);
      H_tmp2.localAccess[i,j,k] = H0.localAccess[i,j,k] * (1 + zeta_tmp2.localAccess[i,j] / h.localAccess[i,j]);
    }

/*
  writeln("Initializing marbl wrapper cell geometry");
  stdout.flush();
  var littleDomain = {0..10, 0..10};
  for (i,j) in littleDomain {
    ref localH = H_tmp1.localSlice(H_tmp1.domain.localSubdomain());
    var thicknesses_reversed: [0..<Nz] real;
    for k in 0..<Nz do thicknesses_reversed[k] = localH[i,j,Nz-1-k];

    var depths_reversed: [0..<Nz] real = + scan thicknesses_reversed[..];
    var midpoints_reversed = depths_reversed - 0.5 * thicknesses_reversed[..];

    marbl_wrappers[i,j].initMarblInstance(Nz, numParSubcols,
      numElementsSurfaceFlux, thicknesses_reversed, depths_reversed, midpoints_reversed, Nz);
  }
  writeln("Done initializing");
*/

  allLocalesBarrier.barrier();
  update_halos(H_tmp1);
  update_halos(H_tmp2);



  calc_volumetric_fluxes(u_tmp, v_tmp, U_n, V_n, H_tmp1);

  // Load MARBL forcings
    get_var(surf_forcing_files[0], 'sss', sss, D2);
    get_var(surf_forcing_files[0], 'sst', sst, D2);
    get_var(surf_forcing_files[0], 'u10_sqr', u10_sqr, D2);
    get_var(surf_forcing_files[0], 'atm_pressure', atm_pressure, D2);
    get_var(surf_forcing_files[0], 'xco2', xco2, D2);
    get_var(surf_forcing_files[0], 'xco2_alt_co2', xco2_alt_co2, D2);
    get_var(surf_forcing_files[0], 'dust_flux_surf', dust_flux_surf, D2);
    get_var(surf_forcing_files[0], 'iron_flux', iron_flux, D2);
    get_var(surf_forcing_files[0], 'nox_flux', nox_flux, D2);
    get_var(surf_forcing_files[0], 'nhy_flux', nhy_flux, D2);
    get_var(surf_forcing_files[0], 'ice_fraction', ice_fraction, D2);

    get_var(int_forcing_files[0], 'dust_flux_int', dust_flux_int, D3);
    get_var(int_forcing_files[0], 'pot_temp', pot_temp, D3);
    get_var(int_forcing_files[0], 'salinity', salinity, D3);
    get_var(int_forcing_files[0], 'surface_sw', surface_sw, D3);
    get_var(int_forcing_files[0], 'par_col_frac', par_col_frac, D3);
    get_var(int_forcing_files[0], 'pressure', pressure, D3);
    get_var(int_forcing_files[0], 'o2_scale_fac', o2_scale_fac, D3);
    get_var(int_forcing_files[0], 'iron_sedflux', iron_sedflux, D3);

}

proc prepare_next_timestep(step : int) {

  // Load the boundary data for the upcoming timestep
    var bryloc = bryfiles[0];
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

proc update_halos(ref arr) {

  // Update the halos. Only need to issue the updateFluff
  // method from one Locale.

  allLocalesBarrier.barrier();
  if (here.id == 0) {
    arr.updateFluff();
  }
  allLocalesBarrier.barrier();

}
