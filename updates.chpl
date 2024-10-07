use INPUTS;
use domains;
use tracers;
use dynamics;

proc prepare_to_timestep(step : int) {

  // Load the boundary data for the current timestep
  coforall loc in Locales do on loc {
    set_bry(bryfiles[step], "temp", tracer_n, D3.localSubdomain());
  }
    update_halos(tracer_n);

  // Load velocity fields for the next timestep
    update_thickness(zeta_np1, H_np1, H0, h, step+1);
    update_dynamics(u_np1, v_np1, U_np1, V_np1, H_np1, step+1);

}

proc prepare_next_timestep(step : int) {

  // Copy next thickness to current thickness
    forall (k,j,i) in D3.localSubdomain() {
      H_n[k,j,i] = H_np1[k,j,i];
    }
//    set_bry(P, P.bryfiles[step+1], "temp", tracer_n, D.rho_3D);

    update_halos(H_n);
//    update_halos(tracer_n);

  // Copy next velocities to current velocities
    forall (k,j,i) in D3_u.localSubdomain() {
      U_n[k,j,i] = U_np1[k,j,i];
    }
    update_halos(U_n);

    forall (k,j,i) in D3_v.localSubdomain() {
      V_n[k,j,i] = V_np1[k,j,i];
    }
    update_halos(V_n);

//    allLocalesBarrier.barrier();

//  // Load velocity fields for the next timestep
//    update_thickness(zeta_np1, H_np1, H0, h, step+2);
//    update_dynamics(u_np1, v_np1, U_np1, V_np1, H_np1, step+2);

}

proc update_dynamics(ref u, ref v, ref U, ref V, ref H, step : int) {

  // Read in u and v

//    u = get_var(P.velfiles[step], "u", D.u_3D);
//    v = get_var(P.velfiles[step], "v", D.v_3D);
//    u = readZarrArrayLocal(P.velfiles[step] + (here.id : string) + '/u/', real(32), 3);
//    v = readZarrArrayLocal(P.velfiles[step] + (here.id : string) + '/v/', real(32), 3);
    u = readZarrArray(velfiles[step] + '/u/', real(32), 3, targetLocales=myTargetLocales3D);
    v = readZarrArray(velfiles[step] + '/v/', real(32), 3, targetLocales=myTargetLocales3D);

  coforall loc in Locales do on loc {
    calc_volumetric_fluxes(u, v, U, V, H);
  }
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

  update_halos(H);

//  allLocalesBarrier.barrier();
//  if (here.id == 0) {
//    H.updateFluff();
//  }
//  allLocalesBarrier.barrier();

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
