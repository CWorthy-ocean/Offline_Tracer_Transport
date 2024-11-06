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


proc Explicit_TimeStep(step : int) {

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //     Step forward the tracer using Easter's Psuedo-compressibility and Strang splitting.          //
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  // Update in x-direction

    RHS_H_U(ktmp, U_n);
    forall (k,j,i) in D3.localSubdomain() {
      H_tilde[k,j,i] = H_n[k,j,i] + dt*ktmp[k,j,i];
    }

    calc_horizontal_fluxes_U(U_n, tmp_U_adv, tracer_n);
    calc_diffusive_fluxes_U(tmp_U_diff, tracer_n, H_n);

    RHS_tr_U(ktmp, tmp_U_adv, tmp_U_diff);

    forall (k,j,i) in D3.localSubdomain() {
      tracer_tilde[k,j,i] =  (tracer_n[k,j,i]*H_n[k,j,i] + dt*ktmp[k,j,i]) / H_tilde[k,j,i];
    }
    update_halos(tracer_tilde);

  // Update in y-direction

    RHS_H_V(ktmp, V_n);
    forall (k,j,i) in D3.localSubdomain() {
      H_dagger[k,j,i] = H_tilde[k,j,i] + dt*ktmp[k,j,i];
    }

    calc_horizontal_fluxes_V(V_n, tmp_V_adv, tracer_tilde);
    calc_diffusive_fluxes_V(tmp_V_diff, tracer_tilde, H_tilde);

    RHS_tr_V(ktmp, tmp_V_adv, tmp_V_diff);

    forall (k,j,i) in D3.localSubdomain() {
      tracer_dagger[k,j,i] =  (tracer_tilde[k,j,i]*H_tilde[k,j,i] + dt*ktmp[k,j,i]) / H_dagger[k,j,i];
    }
    allLocalesBarrier.barrier();

}

proc Implicit_TimeStep(step : int) {

    calc_vertical_diffusion(tracer_dagger, H_dagger);

}

proc RHS_H_U(ref tmp, ref U) {

    forall (k,j,i) in D3.localSubdomain() {
      tmp[k,j,i] = -iarea * (U[k,j,i] - U[k,j,i-1]);
    }

}

proc RHS_H_V(ref tmp, ref V) {

    forall (k,j,i) in D3.localSubdomain() {
      tmp[k,j,i] = -iarea * (V[k,j,i] - V[k,j-1,i]);
    }

}

proc RHS_tr_U(ref tmp, ref U, ref diff_U) {

  forall (k,j,i) in D3.localSubdomain() {
      tmp[k,j,i] = - iarea * (  (U[k,j,i] - U[k,j,i-1])
                                  - (diff_U[k,j,i] - diff_U[k,j,i-1]) );
  }

}

proc RHS_tr_V(ref tmp, ref V, ref diff_V) {

  forall (k,j,i) in D3.localSubdomain() {
      tmp[k,j,i] = - iarea * (  (V[k,j,i] - V[k,j-1,i])
                                  - (diff_V[k,j,i] - diff_V[k,j-1,i]) );
  }

}
