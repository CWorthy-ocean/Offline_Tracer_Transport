use AllLocalesBarriers;
use Time;
//use Zarr;

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

  coforall loc in Locales do on loc {

    initialize_tr();
    initialize_sponge();

    // timestepping loop
      for step in (Nt_start)..(Nt_start+Nt) {

        var t0  : stopwatch;
        var t1a : stopwatch;
        var t1b : stopwatch;
        var t1c : stopwatch;
        var t2a  : stopwatch;
        var t2b  : stopwatch;
        var t2c  : stopwatch;
        var t3a  : stopwatch;
        var t3b  : stopwatch;
        var t3c  : stopwatch;
        var t4  : stopwatch;
        var t5  : stopwatch;

        prepare_to_timestep(step);


        // Step forward thickness
        t0.start();
        Explicit_TimeStep_H();
        t0.stop();

        // Step forward tracers
        t1a.start();
        Explicit_TimeStep_tr(tracers_ts_n, tracers_ts_tilde, tracers_ts_dagger, num_ts_tracers, step);
        t1a.stop();

        t1b.start();
        Explicit_TimeStep_tr(tracers_marbl_n, tracers_marbl_tilde, tracers_marbl_dagger, num_marbl_tracers, step);
        t1b.stop();

        t1c.start();
        Explicit_TimeStep_tr(tracers_other_n, tracers_other_tilde, tracers_other_dagger, num_other_tracers, step);
        t1c.stop();


        t2a.start();
        Implicit_TimeStep(tracers_ts_dagger, num_ts_tracers);
        t2a.stop();

        t2b.start();
        Implicit_TimeStep(tracers_marbl_dagger, num_marbl_tracers);
        t2b.stop();

        t2c.start();
        Implicit_TimeStep(tracers_other_dagger, num_other_tracers);
        t2c.stop();


        // Create polynomial fit to current grid
        t3a.start();
        Polyfit(tracers_ts_n, tracers_ts_dagger, num_ts_tracers);
        t3a.stop();

        t3b.start();
        Polyfit(tracers_marbl_n, tracers_marbl_dagger, num_marbl_tracers);
        t3b.stop();

        t3c.start();
        Polyfit(tracers_other_n, tracers_other_dagger, num_other_tracers);
        t3c.stop();


        // Update fields to prepare for next time step
        t4.start();
        prepare_next_timestep(step);
        t4.stop();

        allLocalesBarrier.barrier();

        t5.start();
        if ((step % output_freq) == 0) {
          WriteOutput(tracers_other_n, D3, "tracer", "stuff", step, 1);
	}
        t5.stop();

        writeln("Locale ", here.id, " time for explicit step (ts): ", t1a.elapsed());
        writeln("Locale ", here.id, " time for explicit step (marbl): ", t1b.elapsed());
        writeln("Locale ", here.id, " time for explicit step (other): ", t1c.elapsed());
        writeln("Locale ", here.id, " time for implicit step (ts): ", t2a.elapsed());
        writeln("Locale ", here.id, " time for implicit step (marbl): ", t2b.elapsed());
        writeln("Locale ", here.id, " time for implicit step (other): ", t2c.elapsed());
        writeln("Locale ", here.id, " time for polyfit (ts): ", t3a.elapsed());
        writeln("Locale ", here.id, " time for polyfit (marbl): ", t3b.elapsed());
        writeln("Locale ", here.id, " time for polyfit (other): ", t3c.elapsed());
        writeln("Locale ", here.id, " time for reading: ", t4.elapsed());
        writeln("Locale ", here.id, " time for writing NetCDF: ", t5.elapsed());
        writeln();


    } // timestepping loop

  } // coforall loop

} // end program
