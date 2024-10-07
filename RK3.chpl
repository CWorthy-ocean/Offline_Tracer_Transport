use INPUTS;
use domains;
use tracers;
use dynamics;
use horizontal_diffusion;
use vertical_diffusion;
//use params;
use hadv_3o_upwind;
use updates;
use NetCDF_IO;
//use utils;

use Math;
use AllLocalesBarriers;
use Time;


////////////////////////////////////////////////////////
//                                                    //
//               Runge-Kutta 3th order                //
//                                                    //
//         y_(n+1)=y_n+(1/6)*h*(k_1+4*k_2+k_3)        //
//                      where                         //
//                  k_1=f(x_n,y_n)                    //
//          k2=f(x_n+(1/2)*h,y_n+(1/2)*h*k_1)         //
//            k3=f(x_n+h,y_n+2*h*k_2−h*k_1)           //
//                                                    //
////////////////////////////////////////////////////////


proc Explicit_TimeStep(step : int) {

  // Calculate horizontal velocities at the (n+1/2) time step.
    calc_half_step_dyn();

  // Calculate thickness at the (n+1/2) time step. This is only used for the diffusion module.
    calc_half_step_tr();

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                       Update H using RK3.                                        //
  //             Here we are following "ALE algorithm, flavor 2", which is on slide 4 of              //
  //                 https://adcroft.github.io/assets/pdf/ALE_workshop_NCWCP_2016.pdf                 //
  //          We are pretending that v_dagger is obtained by reading in the velocity fields.          //
  //   We then need to calculate h_dagger, which is then used to calculate h_dagger * theta_dagger    //
  //                                     in a consistent manner.                                      //
  //////////////////////////////////////////////////////////////////////////////////////////////////////

    RHS_H(k1, U_n, V_n);

    forall (k,j,i) in D3.localSubdomain() {
      ktmp[k,j,i] =  H_n[k,j,i] + 0.5*dt*k1[k,j,i];
    }
    update_halos(ktmp);

    calc_volumetric_fluxes(u_np1h, v_np1h, U_np1h, V_np1h, ktmp);
    RHS_H(k2, U_np1h, V_np1h);

    forall (k,j,i) in D3.localSubdomain() {
      ktmp[k,j,i] =  H_n[k,j,i] - dt*k1[k,j,i] + 2*dt*k2[k,j,i];
    }
    update_halos(ktmp);

    calc_volumetric_fluxes(u_np1, v_np1, U_np1, V_np1, ktmp);

    RHS_H(k3, U_np1, V_np1);

    forall (k,j,i) in D3.localSubdomain() {
      H_dagger[k,j,i] = H_n[k,j,i] + one_sixth*dt*
                           (k1[k,j,i] + 4*k2[k,j,i] + k3[k,j,i]);
    }

  //////////////////////////////////////////
  //       Update tracer using RK3        //
  //////////////////////////////////////////

    calc_horizontal_fluxes(U_n, V_n, tmp_U, tmp_V, tracer_n);
    calc_diffusive_fluxes(tmp_U, tmp_V, tracer_n, H_n);

    RHS_tr(k1, tmp_U, tmp_V, tmp_U, tmp_V);

    forall (k,j,i) in D3.localSubdomain() {
      ktmp[k,j,i] =  (tracer_n[k,j,i]*H_n[k,j,i] + 0.5*dt*k1[k,j,i]) / H_np1h[k,j,i];
    }
    update_halos(ktmp);

    calc_horizontal_fluxes(U_np1h, V_np1h, tmp_U, tmp_V, ktmp);
    calc_diffusive_fluxes(tmp_U, tmp_V, ktmp, H_np1h);

    RHS_tr(k2, tmp_U, tmp_V, tmp_U, tmp_V);

    forall (k,j,i) in D3.localSubdomain() {
      ktmp[k,j,i] =  (tracer_n[k,j,i]*H_np1h[k,j,i] - dt*k1[k,j,i] + 2*dt*k2[k,j,i]) / H_np1[k,j,i];
    }
    update_halos(ktmp);

    calc_horizontal_fluxes(U_np1, V_np1, tmp_U, tmp_V, ktmp);
    calc_diffusive_fluxes(tmp_U, tmp_V, ktmp, H_np1);

    RHS_tr(k3, tmp_U, tmp_V, tmp_U, tmp_V);

    forall (k,j,i) in D3.localSubdomain() {
      tracer_dagger[k,j,i] = (tracer_n[k,j,i]*H_n[k,j,i] + one_sixth*dt*
                           (k1[k,j,i] + 4*k2[k,j,i] + k3[k,j,i])) / H_dagger[k,j,i];
    }

    allLocalesBarrier.barrier();

//WriteOutput(H_n, "H_n", "stuff", step);
//WriteOutput(H_np1h, "H_np1h", "stuff", step);
//WriteOutput(H_np1, "H_np1", "stuff", step);
//WriteOutput(H_dagger, "H_dagger", "stuff", step);
//WriteOutput(tracer_n, "tracer_n", "stuff", step);
//WriteOutput(tracer_dagger, "tracer_dagger", "stuff", step);
//WriteOutput(u_n, "u_n", "stuff", step);
//WriteOutput(u_np1h, "u_np1h", "stuff", step);
//WriteOutput(u_np1, "u_np1", "stuff", step);
//WriteOutput(v_n, "v_n", "stuff", step);
//WriteOutput(v_np1h, "v_np1h", "stuff", step);
//WriteOutput(v_np1, "v_np1", "stuff", step);

}

proc Implicit_TimeStep(step : int) {

    calc_vertical_diffusion(tracer_dagger, H_dagger);

}

/*
proc update_fields(step : int) {

    forall (k,j,i) in D3.localSubdomain() {
      H_n[k,j,i] = H_np1[k,j,i];
    }
    set_bry(P, P.bryfiles[step+1], "temp", tracer_n, D.rho_3D);

    update_halos(H_n);
    update_halos(tracer_n);

    U_n = U_np1;
    V_n = V_np1;
    allLocalesBarrier.barrier();

  // Load velocity fields for the next timestep
    update_thickness(zeta_np1, H_np1, H0, h, step+2);
    update_dynamics(u_np1, v_np1, U_np1, V_np1, H_np1, step+2);

}
*/

proc RHS_H(ref tmp, ref U, ref V) {

  /////////////////////////////////////////
  //  Calculate tracer field at (n+1) timestep  //
  /////////////////////////////////////////

/*
  if (here.id == 0) {
    forall (k,j,i) in {D.rho_3D.dim[0], (D.rho_3D.first[1]+1)..(D.rho_3D.last[1]-1), (D.rho_3D.first[2]+1)..D.rho_3D.last[2]}  {
      tmp[k,j,i] = - P.iarea * (  (U[k,j,i] - U[k,j,i-1])
                        + (V[k,j,i] - V[k,j-1,i])  );
    }
  }
  else if (here.id == (Locales.size-1)) {
    forall (k,j,i) in {D.rho_3D.dim[0], (D.rho_3D.first[1]+1)..(D.rho_3D.last[1]-1), D.rho_3D.first[2]..(D.rho_3D.last[2]-1)}  {
      tmp[k,j,i] = - P.iarea * (  (U[k,j,i] - U[k,j,i-1])
                        + (V[k,j,i] - V[k,j-1,i])  );
    }
  }
  else {
    forall (k,j,i) in {D.rho_3D.dim[0], (D.rho_3D.first[1]+1)..(D.rho_3D.last[1]-1), D.rho_3D.first[2]..D.rho_3D.last[2]}  {
      tmp[k,j,i] = - P.iarea *(  (U[k,j,i] - U[k,j,i-1])
                        + (V[k,j,i] - V[k,j-1,i])  );
    }
  }

  allLocalesBarrier.barrier();
*/

    forall (k,j,i) in D3.localSubdomain() {
      tmp[k,j,i] = -iarea * (  (U[k,j,i] - U[k,j,i-1])
                        + (V[k,j,i] - V[k,j-1,i])  );
    }

}




proc RHS_tr(ref tmp, ref U, ref V, ref diff_U, ref diff_V) {

  /////////////////////////////////////////
  //  Calculate tracer field at (n+1) timestep  //
  /////////////////////////////////////////
/*
  if (here.id == 0) {
    forall (k,j,i) in {D.rho_3D.dim[0], (D.rho_3D.first[1]+1)..(D.rho_3D.last[1]-1), (D.rho_3D.first[2]+1)..D.rho_3D.last[2]}  {
      tmp[k,j,i] = - P.iarea * (  (U[k,j,i] - U[k,j,i-1])
                                  + (V[k,j,i] - V[k,j-1,i])
                                  - (diff_U[k,j,i] - diff_U[k,j,i-1])
                                  - (diff_V[k,j,i] - diff_V[k,j-1,i]) );
    }
  }
  else if (here.id == (Locales.size-1)) {
    forall (k,j,i) in {D.rho_3D.dim[0], (D.rho_3D.first[1]+1)..(D.rho_3D.last[1]-1), D.rho_3D.first[2]..(D.rho_3D.last[2]-1)}  {
      tmp[k,j,i] = - P.iarea * (  (U[k,j,i] - U[k,j,i-1])
                                  + (V[k,j,i] - V[k,j-1,i])
                                  - (diff_U[k,j,i] - diff_U[k,j,i-1])
                                  - (diff_V[k,j,i] - diff_V[k,j-1,i]) );
    }
  }
  else {
    forall (k,j,i) in {D.rho_3D.dim[0], (D.rho_3D.first[1]+1)..(D.rho_3D.last[1]-1), D.rho_3D.first[2]..D.rho_3D.last[2]}  {
      tmp[k,j,i] = - P.iarea * (  (U[k,j,i] - U[k,j,i-1])
                                  + (V[k,j,i] - V[k,j-1,i])
                                  - (diff_U[k,j,i] - diff_U[k,j,i-1])
                                  - (diff_V[k,j,i] - diff_V[k,j-1,i]) );
    }
  }

*/

  forall (k,j,i) in D3.localSubdomain() {
      tmp[k,j,i] = - iarea * (  (U[k,j,i] - U[k,j,i-1])
                                  + (V[k,j,i] - V[k,j-1,i])
                                  - (diff_U[k,j,i] - diff_U[k,j,i-1])
                                  - (diff_V[k,j,i] - diff_V[k,j-1,i]) );
  }

  allLocalesBarrier.barrier();

}
