use AllLocalesBarriers;
//use Zarr;

use domains;
use INPUTS;
use NetCDF_IO;
use tracers;
use updates;

  var U_n : [D3_u] real;
  var tmp_U_adv : [D3_u] real;
  var tmp_U_diff : [D3_u] real;

  var V_n : [D3_v] real;
  var tmp_V_adv : [D3_v] real;
  var tmp_V_diff : [D3_v] real;

  var u_n : [D3_u] real;
  var v_n : [D3_v] real;

proc calc_volumetric_fluxes(ref u, ref v, ref U, ref V, ref H) {

  forall (i,j,k) in D3_u.localSubdomain() {
    U.localAccess[i,j,k] = u.localAccess[i,j,k] * 0.5 * (H.localAccess[i,j,k] + H.localAccess[i+1,j,k]) * dy;
  }

  forall (i,j,k) in D3_v.localSubdomain() {
    V.localAccess[i,j,k] = v.localAccess[i,j,k] * 0.5 * (H.localAccess[i,j,k] + H.localAccess[i,j+1,k]) * dx;
  }

  allLocalesBarrier.barrier();

  update_halos(U);
  update_halos(V);

}

