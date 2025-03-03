use INPUTS;
use domains;

use StencilDist;
use FileSystem;



var bryfiles = glob(boundary_files);
var surf_forcing_files = glob(surface_forcing_files);
var int_forcing_files = glob(interior_forcing_files);

// Surface forcings
var sss            : [D2] real;
var sst            : [D2] real;
var u10_sqr        : [D2] real;
var atm_pressure   : [D2] real;
var xco2           : [D2] real;
var xco2_alt_co2   : [D2] real;
var dust_flux_surf : [D2] real;
var iron_flux      : [D2] real;
var nox_flux       : [D2] real;
var nhy_flux       : [D2] real;
var ice_fraction   : [D2] real;

// Interior forcings
var dust_flux_int  : [D3] real;
var pot_temp       : [D3] real;
var salinity       : [D3] real;
var surface_sw     : [D3] real;
var par_col_frac   : [D3] real;
var pressure       : [D3] real;
var o2_scale_fac   : [D3] real;
var iron_sedflux   : [D3] real;
