use tracers;
use INPUTS;
use domains;
use forcings;
use Marbl;


proc step_marbl_wrappers(step) {
  writeln("stepping marbl wrappers for step ", step);
  for (i,j) in marbl_wrappers.domain.localSubdomain() {

    ref marbl_wrapper = marbl_wrappers[i,j];
    // Populate surface flux forcing values
    marbl_wrapper.setSurfaceFluxForcingValue("sss", sss[i,j]);

    // Populate surface tracers

    // Call surface flux compute


    // Populate interior tendency forcing values

    // Populate interior tracers

    // Call interior tendency compute
  }
}
