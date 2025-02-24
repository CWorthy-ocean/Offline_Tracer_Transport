use AllLocalesBarriers;
use Time;

use domains;
use dynamics;
use hadv_3o_upwind_superbee;
use horizontal_diffusion;
use INPUTS;
use NetCDF_IO;
use tracers;
use updates;
use vertical_diffusion;

proc Explicit_TimeStep_H() {

  // Update in x-direction

    RHS_H_U(ktmp, U_n);
    forall (i,j,k) in D3_int_3D.localSubdomain() {
      H_tilde[i,j,k] = H_n[i,j,k] + dt*ktmp[i,j,k];
    }

  // Update in y-direction

    RHS_H_V(ktmp, V_n);
    forall (i,j,k) in D3_int_3D.localSubdomain() {
      H_dagger[i,j,k] = H_tilde[i,j,k] + dt*ktmp[i,j,k];
    }
}

proc Explicit_TimeStep_tr(ref tracer_n, ref tracer_tilde, ref tracer_dagger, num_tracers, step : int) {

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //     Step forward the tracer using Easter's Psuedo-compressibility and Strang splitting.          //
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  if (num_tracers > 0) {

  // Update in x-direction

    for t in 1..num_tracers {
      calc_horizontal_fluxes_U(U_n, tmp_U_adv, tracer_n, t);
      calc_diffusive_fluxes_U(tmp_U_diff, tracer_n, H_n, t);

      RHS_tr_U(ktmp_tr, tmp_U_adv, tmp_U_diff);

      forall (i,j,k) in D3_int_3D.localSubdomain() {
        tracer_tilde[i,j,t,k] =  (tracer_n[i,j,t,k]*H_n[i,j,k] + dt*ktmp_tr[i,j,k]) / H_tilde[i,j,k];
      }
    }

    update_halos(tracer_tilde);

  // Update in y-direction

    for t in 1..num_tracers {
      calc_horizontal_fluxes_V(V_n, tmp_V_adv, tracer_tilde, t);
      calc_diffusive_fluxes_V(tmp_V_diff, tracer_tilde, H_tilde, t);

      RHS_tr_V(ktmp_tr, tmp_V_adv, tmp_V_diff);

      forall (i,j,k) in D3_int_3D.localSubdomain() {
        tracer_dagger[i,j,t,k] =  (tracer_tilde[i,j,t,k]*H_tilde[i,j,k] + dt*ktmp_tr[i,j,k]) / H_dagger[i,j,k];
      }
    }

    allLocalesBarrier.barrier();

  } // num_tracers > 0

}

proc Implicit_TimeStep(ref tracer_dagger, num_tracers) {

  if (num_tracers > 0) {

    calc_vertical_diffusion(tracer_dagger, H_dagger);

  } // num_tracers > 0

}

proc RHS_H_U(ref tmp, ref U) {

  forall (i,j,k) in D3_int_3D.localSubdomain() {
    tmp[i,j,k] = -iarea * (U[i,j,k] - U[i-1,j,k]);
  }

}

proc RHS_H_V(ref tmp, ref V) {

  forall (i,j,k) in D3_int_3D.localSubdomain() {
    tmp[i,j,k] = -iarea * (V[i,j,k] - V[i,j-1,k]);
  }

}

proc RHS_tr_U(ref tmp, ref adv_U, ref diff_U) {

  forall (i,j,k) in D3_int_3D.localSubdomain() {
    tmp[i,j,k] = - iarea * (  (adv_U[i,j,k] - adv_U[i-1,j,k])
                            - (diff_U[i,j,k] - diff_U[i-1,j,k]) );
  }

}

proc RHS_tr_V(ref tmp, ref adv_V, ref diff_V) {

  forall (i,j,k) in D3_int_3D.localSubdomain() {
      tmp[i,j,k] = - iarea * (  (adv_V[i,j,k] - adv_V[i,j-1,k])
                              - (diff_V[i,j,k] - diff_V[i,j-1,k]) );
  }

}
