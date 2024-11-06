use AllLocalesBarriers;
use Zarr;

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

    forall (k,j,i) in D3_u.localSubdomain() {
      U[k,j,i] = u[k,j,i] * 0.5 * (H[k,j,i] + H[k,j,i+1]) * dy;
    }

    forall (k,j,i) in D3_v.localSubdomain() {
      V[k,j,i] = v[k,j,i] * 0.5 * (H[k,j,i] + H[k,j+1,i]) * dx;
    }

    allLocalesBarrier.barrier();

    update_halos(U);
    update_halos(V);

}

