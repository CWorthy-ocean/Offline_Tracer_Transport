use IO;
use BlockDist;
use StencilDist;
use Time;
use AutoMath;
use LinearAlgebra;
use IO.FormattedIO;
use Math;
use AllLocalesBarriers;
use Zarr;

use INPUTS;
use domains;

use dynamics;
use horizontal_diffusion;
use tracers;
use updates;

use NetCDF_IO;
use RK3;
use PPM;
//use params;


proc main() {

  var t : stopwatch;
  t.start();

  initialize_tr();
  initialize_sponge();
  initialize_dynamics();

  // timestepping loop
    for step in (Nt_start)..(Nt_start+Nt) {

      prepare_to_timestep(step);

      coforall loc in Locales do on loc {
          var t1 : stopwatch;
          var t2 : stopwatch;
          var t3 : stopwatch;
          var t4 : stopwatch;

          // Step forward

          t1.start();
          Explicit_TimeStep(step);
          t1.stop();

          t2.start();
          Implicit_TimeStep(step);
          t2.stop();

          // Create polynomial fit to current grid

          t3.start();
          Polyfit();
          t3.stop();

          // Update fields to prepare for next time step

          t4.start();
          prepare_next_timestep(step);
          t4.stop();

// NEED TO WRITE ZARR OUTSIDE OF COFORALL LOOP
          WriteOutput(tracer_n, "after", "stuff", step);
          allLocalesBarrier.barrier();

          writeln("Locale ", here.id, " time for explicit step: ", t1.elapsed());
          writeln("Locale ", here.id, " time for implicit step: ", t2.elapsed());
          writeln("Locale ", here.id, " time for polyfit: ", t3.elapsed());
//          writeln("Locale ", here.id, " time for reading: ", t4.elapsed());

      } // coforall loop

//  var (maxVal, maxLoc) = maxloc reduce zip(tracer_dagger, tracer_dagger.domain);
//  var (minVal, minLoc) = minloc reduce zip(tracer_dagger, tracer_dagger.domain);

//  writeln("Max of v is ", maxVal, ' at ', maxLoc);
//  writeln("Min of v is ", minVal, ' at ', minLoc);

    } // timestepping loop

  t.stop();
  writeln("Program finished in ", t.elapsed(), " seconds.");


} // end program
