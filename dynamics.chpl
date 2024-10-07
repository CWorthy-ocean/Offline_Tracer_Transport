//use params;
use INPUTS;
use domains;
use tracers;
//use sigma_coordinate;
use NetCDF_IO;

use Zarr;
use AllLocalesBarriers;

  var U_n : [D3_u] real;
  var U_np1h : [D3_u] real;
  var U_np1 : [D3_u] real;
  var tmp_U : [D3_u] real;

  var V_n : [D3_v] real;
  var V_np1h : [D3_v] real;
  var V_np1 : [D3_v] real;
  var tmp_V : [D3_v] real;

  var u_n : [D3_u] real;
  var u_np1h : [D3_u] real;
  var u_np1 : [D3_u] real;

  var v_n : [D3_v] real;
  var v_np1h : [D3_v] real;
  var v_np1 : [D3_v] real;

/*
record Dynamics {

  var r0, r1, r2 = 1..0;
  var u0, u1, u2 = 1..0;
  var v0, v1, v2 = 1..0;
  var w0, w1, w2 = 1..0;

  var U_n : [u0, u1, u2] real;
  var U_np1h : [u0, u1, u2] real;
  var U_np1 : [u0, u1, u2] real;
  var tmp_U : [u0, u1, u2] real;

  var V_n : [v0, v1, v2] real;
  var V_np1h : [v0, v1, v2] real;
  var V_np1 : [v0, v1, v2] real;
  var tmp_V : [v0, v1, v2] real;

  var u_n : [u0, u1, u2] real;
  var u_np1h : [u0, u1, u2] real;
  var u_np1 : [u0, u1, u2] real;

  var v_n : [v0, v1, v2] real;
  var v_np1h : [v0, v1, v2] real;
  var v_np1 : [v0, v1, v2] real;


  proc init(arg: Domains) {

    this.r0 = arg.rho_3D.dim[0];
    this.r1 = arg.rho_3D.dim[1];
    this.r2 = arg.rho_3D.dim[2];

    this.u0 = arg.u_3D.dim[0];
    this.u1 = arg.u_3D.dim[1];
    this.u2 = arg.u_3D.dim[2];

    this.v0 = arg.v_3D.dim[0];
    this.v1 = arg.v_3D.dim[1];
    this.v2 = arg.v_3D.dim[2];

  }

}
*/

proc initialize_dynamics() {

  // Read in u and v

//    u = get_var(P.velfiles[step], "u", D.u_3D);
//    v = get_var(P.velfiles[step], "v", D.v_3D);
    u_n = readZarrArray(velfiles[Nt_start] + '/u/', real(32), 3, targetLocales=myTargetLocales3D);
    v_n = readZarrArray(velfiles[Nt_start] + '/v/', real(32), 3, targetLocales=myTargetLocales3D);
//    u_np1 = readZarrArray(velfiles[Nt_start+1] + '/u/', real(32), 3, targetLocales=myTargetLocales3D);
//    v_np1 = readZarrArray(velfiles[Nt_start+1] + '/v/', real(32), 3, targetLocales=myTargetLocales3D);

//    u_n.updateFluff();
//    v_n.updateFluff();
//    u_np1.updateFluff();
//    v_np1.updateFluff();

    coforall loc in Locales do on loc {
      calc_volumetric_fluxes(u_n, v_n, U_n, V_n, H_n);
      calc_volumetric_fluxes(u_np1, v_np1, U_np1, V_np1, H_np1);
    }

/*
  coforall loc in Locales do on loc {
    WriteOutput(u_n, "after", "stuff", 1);
  }
*/

}

/*
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
*/

proc calc_volumetric_fluxes(ref u, ref v, ref U, ref V, ref H) {

//    var fluffDom_u = D3_u.localSubdomain().expand(fluff3D);
//    forall (k,j,i) in fluffDom_u {
    forall (k,j,i) in D3_u.localSubdomain() {
      U[k,j,i] = u[k,j,i] * 0.5 * (H[k,j,i] + H[k,j,i+1]) * dy;
    }

//    var fluffDom_v = D3_v.localSubdomain().expand(fluff3D);
//    forall (k,j,i) in fluffDom_v {
    forall (k,j,i) in D3_v.localSubdomain() {
      V[k,j,i] = v[k,j,i] * 0.5 * (H[k,j,i] + H[k,j+1,i]) * dx;
    }

/*
    allLocalesBarrier.barrier();
    if (here.id == 0) {
      U.updateFluff();
      V.updateFluff();
    }
    allLocalesBarrier.barrier();
*/
    update_halos(U);
    update_halos(V);

}

proc calc_half_step_dyn() {

//  var fluffDom_u = D3_u.localSubdomain().expand(fluff3D);
//  forall (k,j,i) in fluffDom_u {
  forall (k,j,i) in D3_u.localSubdomain() {
    u_np1h[k,j,i] = 0.5 * (u_n[k,j,i] + u_np1[k,j,i]);
  }

//  var fluffDom_v = D3_v.localSubdomain().expand(fluff3D);
//  forall (k,j,i) in fluffDom_v {
  forall (k,j,i) in D3_v.localSubdomain() {
    v_np1h[k,j,i] = 0.5 * (v_n[k,j,i] + v_np1[k,j,i]);
  }

//  allLocalesBarrier.barrier();
//  if (here.id == 0) {
//    u_np1h.updateFluff();
//    v_np1h.updateFluff();
//  }
//  allLocalesBarrier.barrier();

}
