use AllLocalesBarriers;
use Time;
use Zarr;

use domains;
use dynamics;
use horizontal_diffusion;
use INPUTS;
use NetCDF_IO;
use PPM;
use forward_step;
use tracers;
use updates;

proc main() {

  initialize_tr();
  initialize_sponge();

  // timestepping loop
    for step in (Nt_start)..(Nt_start+Nt) {

      var t0 : stopwatch;

      prepare_to_timestep(step);

      coforall loc in Locales do on loc {
          var t1 : stopwatch;
          var t2 : stopwatch;
          var t3 : stopwatch;
          var t4 : stopwatch;
          var t5 : stopwatch;

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

          allLocalesBarrier.barrier();

          t5.start();
          WriteOutput(tracer_n, "tracer", "stuff", step);
          t5.stop();

          writeln("Locale ", here.id, " time for explicit step: ", t1.elapsed());
          writeln("Locale ", here.id, " time for implicit step: ", t2.elapsed());
          writeln("Locale ", here.id, " time for polyfit: ", t3.elapsed());
          writeln("Locale ", here.id, " time for reading: ", t4.elapsed());
          writeln("Locale ", here.id, " time for writing NetCDF: ", t5.elapsed());

      } // coforall loop

    } // timestepping loop

} // end program
