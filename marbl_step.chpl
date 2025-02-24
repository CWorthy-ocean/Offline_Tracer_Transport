use tracers;
use INPUTS;
use domains;

use Marbl;


proc step_marbl_wrappers(step) {
  writeln("stepping marbl wrappers for step ", step);
  for (i,j) in marbl_wrappers.domain.localSubdomain() {
    // Populate surface flux forcing values

    // Populate surface tracers

    // Call surface flux compute

    
    // Populate interior tendency forcing values

    // Populate interior tracers

    // Call interior tendency compute
  }
}