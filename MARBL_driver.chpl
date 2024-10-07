/*
-----------------------------------------------------------------------
     MODULE: marbl_driver

     DESCRIPTION:
     This module provides an interface between MARBL and OTT.
     It configures an array of MARBL "instances", each of which has a particular number of tracers
     and diagnostics, and each of which is assigned to a certain subdomain on each node (Locale).
     At each time step and gridpoint the instance receives tracer/forcing
     information from OTT and sends that information to MARBL (Fortran) via a shared object, which
     is used to calculate interior tendencies of tracers and surface fluxes. OTT then updates its tracer and diagnostics
     arrays based on these calculations.

     NOTES:
     Configured to run with MARBL 0.45 on the "development" branch:
     https://github.com/marbl-ecosys/MARBL/releases/tag/marbl0.45.0

     AUTHOR: SB
     UPDATED: 2024-09-03
-----------------------------------------------------------------------

     KEY ELEMENTS:
     Chapel access to key phases of MARBL:
       init: pass in settings to configure MARBL
       surface_flux_compute: compute the tracer surface fluxes
       interior_tendency_compute: compute the interior source/sink terms
       shutdown: close the simulation
     Additional requirements:
       Parse interior_tendency_diags and surface_flux_diags to write to output file. (+ metadata)
       Handle "saved state", e.g., surface_flux_saved_state to enable exact restart. (+ metadata)
       Support MARBL timers and logging
       (low priority) support running means/global average
*/

use IO;

use INPUTS;
use domains;
use tracers;


var u10_sqr_ind, sss_ind, sst_ind, ifrac_ind, dust_dep_ind : int;
var fe_dep_ind, nox_flux_ind, nhy_flux_ind, atmpress_ind, xco2_ind, xco2_al_ind : int;

var dustflux_ind, PAR_col_frac_ind, surf_shortwave_ind : int;
var potemp_ind, salinity_ind, pressure_ind : int;
var fesedflux_ind, o2_scalef_ind, remin_scalef_ind : int;

/*------------------------------------------------*/
/*------------------------------------------------*/


class marbl_shared_object {


 proc init() { }
  // Calls init (line 196, marbl_interface.F90)
  // !!! Something goes here
 }

/*------------------------------------------------*/

  proc put_settings() {
    // Reads in the parameter settings from marbl_in,
    // then calls put_settings_file_line (line 562, marbl_interface.F90)

    // Read marbl_in into buffer
      var infile = open("marbl_in", ioMode.r);
      var reader = infile.reader();

    // Line-by-line read
      for line in reader.lines() {
        // Call marbl_instance%put_setting(namelist_line)
        // !!! Something goes here
      }

      writeln("Successfully read marbl_in to end of file.");
   }

/*------------------------------------------------*/

  proc surface_flux_compute() {
  // Calls surface_flux_compute (line 999, marbl_interface.F90)
  // !!! Something goes here
  }

/*------------------------------------------------*/

  proc interior_tendency_compute() {
  // Calls interior_tendency_compute (line 940, marbl_interface.F90)
  // !!! Something goes here
  }

/*------------------------------------------------*/
/*------------------------------------------------*/

  proc do_column_physics(istr : int,iend : int,jstr : int,jend : int, tracer_array : ?) {
    // Calculate surface fluxes, interior tracer tendencies,
    // saved state variables, and diagnostics using MARBL,
    // then update values in OTT.

    // Populate MARBL instance surface flux forcing values:
      for j in jstr..jend {
        for i in istr..iend {

          if (mask_rho[i,j] == 1) {
            // sea surface salinity
            if (sss_ind > 0) {
// !!!        marbl_instance%surface_flux_forcings(sss_ind)%field_0d(1)
                   = tracer_array(i,j,nz,nnew,isalt) ! (psu)
            }
            // sea surface temperature
            if (sst_ind > 0) {
// !!!        marbl_instance%surface_flux_forcings(sst_ind)%field_0d(1)
                   = tracer_array(i,j,nz,nnew,itemp) ! (degC)
            }
            // ice fraction (unused)
            if (ifrac_ind > 0) {
// !!!        marbl_instance%surface_flux_forcings(ifrac_ind)%field_0d(1)=0
            }
            // squared 10m windspeed
            if (u10_sqr_ind > 0) {
// !!!        marbl_instance%surface_flux_forcings(u10_sqr_ind)%field_0d(1)=
                    (uwnd(i,j)**2) + (vwnd(i,j)**2) ! (m2/s2)
            }
            // atmospheric pressure (constant)
            if (atmpress_ind > 0) {
              marbl_instance%surface_flux_forcings(atmpress_ind)%field_0d(1)=1.
            }
            // pco2
            if (xco2_ind > 0) {
               marbl_instance%surface_flux_forcings(xco2_ind)%field_0d(1)=
                   pco2air(i,j) !(ppm)
            }
            // pco2 (alternative co2 forcing)
            if (xco2_alt_ind > 0) {
               marbl_instance%surface_flux_forcings(xco2_alt_ind)%field_0d(1)=
                   pco2air_alt(i,j) !(ppm)
            }
            // dust
            if (dust_dep_ind > 0) {
               marbl_instance%surface_flux_forcings(dust_dep_ind)%field_0d(1)
                   = dust(i,j) !(kg/m2/s)
            }
            // iron deposition
            if (fe_dep_ind > 0) {
               marbl_instance%surface_flux_forcings(fe_dep_ind)%field_0d(1)=
                   iron(i,j)*0.01 ! nmol/cm2/s (ROMS) -> mmol/m2/s (MARBL MKS)
            }
            // nox
            if (nox_flux_ind > 0) {
               marbl_instance%surface_flux_forcings(nox_flux_ind)%field_0d(1)=
                   nox(i,j)*71394.200220751 ! kg(N)/m2/s (ROMS) -> mmol/m2/s (MARBL MKS)
            }
            // nhy
            if (nhy_flux_ind > 0) {
               marbl_instance%surface_flux_forcings(nhy_flux_ind)%field_0d(1)=
                   nhy(i,j)*71394.200220751 ! kg(N)/m2/s (ROMS) -> mmol/m2/s (MARBL MKS)
            }

            // Surface tracer values
              for m in 1..nt_marbl {
                marbl_instance%tracers_at_surface(1,m)=
                   tracer_array(i,j,nz,nnew, m+(NT-nt_marbl))
              }

            //// Compute surface fluxes using MARBL
              // !!! call marbl_instance%surface_flux_compute()

            // -------------------------------------------------------

            // Populate MARBL instance interior tendency forcing values:
              // First, update domain in MARBL instance to local geometry:

// !!!                MARBL_instance%domain%zw(:)      = -z_w(i,j,nz-1:0:-1) ! bottom interface depth
// !!!                MARBL_instance%domain%zt(:)      = -z_r(i,j,nz  :1:-1) ! centre depth
// !!!                MARBL_instance%domain%delta_z(:) = Hz(  i,j,nz  :1:-1) ! thickness
// !!!                MARBL_instance%domain%kmt        = nz                  ! number of active levels

              // dust flux
              if (dustflux_ind > 0) {
                MARBL_instance%interior_tendency_forcings(dustflux_ind)%field_0d(1)=
                  dust(i,j)
              }
              // PAR subcolumns (unused as the dimension of this field is the number of ice categories and we have none)
              if (PAR_col_frac_ind > 0) {
                continue
              }
              // surface shortwave
              if (surf_shortwave_ind > 0) {
                MARBL_instance%interior_tendency_forcings(surf_shortwave_ind)%field_1d(1,1)=
                   srflx(i,j) * Cp * rho0
              }
              // potential temperature
              if (potemp_ind > 0) {
                MARBL_instance%interior_tendency_forcings(potemp_ind)%field_1d(1,:)
                   = tracer_array(i,j,nz:1:-1,nnew,itemp)
              }
              // salinity
              if (salinity_ind > 0) {
                MARBL_instance%interior_tendency_forcings(salinity_ind)%field_1d(1,:)
                   =tracer_array(i,j,nz:1:-1,nnew,isalt)
              }
              // pressure (currently calculating from depth)
              if (pressure_ind > 0) {
                MARBL_instance%interior_tendency_forcings(pressure_ind)%field_1d(1,:)=
                   -z_r(i,j,nz:1:-1)*0.1
              }
              // iron flux
              if (fesedflux_ind > 0) {
                MARBL_instance%interior_tendency_forcings(fesedflux_ind)%field_1d(1,:)=0.
              }
              // column tracers
              for m in 1..nt_marbl {
                marbl_instance%tracers(m,:)=
                   tracer_array(i,j,nz:1:-1,nnew, m+(NT-nt_marbl) )
              }

              // Populate bottom fluxes (0 for now)
//                marbl_instance%bot_flux_to_tend(:)=0.
//                marbl_instance%bot_flux_to_tend(Nz)=
//                1./MARBL_instance%domain%delta_z(Nz)


              //// Compute interior tendencies using MARBL:
                // !!! call marbl_instance%interior_tendency_compute()


              //// Apply calculated increments to ROMS tracer array:
              for m in 1..nt_marbl {
                tracer_array(i,j,:,nnew, m+(NT-nt_marbl) ) =
                   tracer_array(i,j,:,nnew, m+(NT-nt_marbl) ) +
                   marbl_instance%interior_tendencies(m,nz:1:-1)*dt
              }
        }
      }

}

/*------------------------------------------------*/
/*------------------------------------------------*/

proc init_marbl () {

  coforall loc in Locales do on loc {

    // !!! Pretty certain that marbl_so will not persist in memory after this coforall ends.  How to fix?
    var marbl_so : [1..here.maxTaskPar] owned marbl_shared_object?;
    forall i in 1..here.maxTaskPar do
      marbl_so[i] = create_marbl_instance();

      // Apply the settings from marbl_in to the MARBL instance
      marbl_so[i].put_settings();

      // Initialize the MARBL instance
      marbl_so[i].init();

      if (here.id == 0) {
        // Confirm how many tracers MARBL is using
          // !!! Can we access/call the MARBL instance directly?
          // nt_marbl=size(marbl_instance%tracer_metadata)

        // Print MARBL internal log
          // !!! call print_marbl_log(marbl_instance%StatusLog)
          // !!! call marbl_instance%StatusLog%erase()

        // Write info about surface forcings
          writeln("Here are the SF forcings MARBL requested:");
          // !!! Can we call the marbl instance to get the surface_flux_forcing objects?
          // for idx in 1..size(marbl_instance%surface_flux_forcings) {
          //  print *, 'var: ', (trim(MARBL_instance%surface_flux_forcings(idx)%metadata%varname))
          //  print *, 'units:',(trim(MARBL_instance%surface_flux_forcings(idx)%metadata%field_units))

        // Write info about interior tendency forcings
          writeln("Here are the int. tend. forcings MARBL requested:");
          // !!! Can we call the marbl instance to get the interior_tendency_forcings objects?
          // for idx in 1..size(marbl_instance%interior_tendency_forcings) {
          //  print *, 'var: ', (trim(MARBL_instance%interior_tendency_forcings(idx)%metadata%varname))
          //  print *, 'units:',(trim(MARBL_instance%interior_tendency_forcings(idx)%metadata%field_units))

      }

      // Determine indices of MARBL surface forcings
      for idx in 1..size(marbl_instance%surface_flux_forcings) {

        select (trim(MARBL_instance%surface_flux_forcings(idx)%metadata%varname)) {
          when 'u10_sqr'              {u10_sqr_ind = idx;}
          when 'sss'                  {sss_ind = idx;}
          when 'sst'                  {sst_ind = idx;}
          when 'Ice Fraction'         {ifrac_ind = idx;}
          when 'Dust Flux'            {dust_dep_ind = idx;}
          when 'Iron Flux'            {fe_dep_ind = idx;}
          when 'NOx Flux'             {nox_flux_ind = idx;}
          when 'NHy Flux'             {nhy_flux_ind = idx;}
          when 'Atmospheric Pressure' {atmpress_ind = idx;}
          when 'xco2'                 {xco2_ind = idx;}
          when 'xco2_alt_co2'         {xco2_alt_ind = idx;}
          otherwise {
            writeln("Additional forcing requested but not indexed: ",
              trim(MARBL_instance%surface_flux_forcings(idx)
              %metadata%varname));
          }
        }
      }

      // Determine indices of MARBL interior tendency forcings
      for idx in 1..size(marbl_instance%surface_flux_forcings) {

        select (trim(MARBL_instance%interior_tendency_forcings(idx)%metadata%varname)) {
          when 'Dust Flux'                        {dustflux_ind = idx;}
          when 'PAR Column Fraction'              {PAR_col_frac_ind = idx;}
          when 'Surface Shortwave'                {surf_shortwave_ind = idx;}
          when 'Potential Temperature'            {potemp_ind = idx;}
          when 'Salinity'                         {salinity_ind = idx;}
          when 'Pressure'                         {pressure_ind = idx;}
          when 'Iron Sediment Flux'               {fesedflux_ind = idx;}
          when 'O2 Consumption Scale Factor'      {o2_scalef_ind = idx;}
          when 'Particulate Remin Scale Factor'   {remin_scalef_ind = idx;}
          otherwise {
            writeln("Additional forcing requested but not indexed: ",
              trim(MARBL_instance%interior_tendency_forcings(idx)
              %metadata%varname));
          }
        }
      }
  } // coforall
}   // init_marbl

/*----------------------------------------------------------------------- */

proc create_marbl_instance() {

  // Allocate memory for a MARBL instance
    !!! Something goes here
}


