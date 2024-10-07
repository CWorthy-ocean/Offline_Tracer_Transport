use FileSystem;

/* INPUT PARAMETERS */

//config const gridfile       = '/glade/derecho/scratch/bachman/UCLA-ROMS/Work/Iceland1/INPUT/Iceland1_grd.nc';
//config const velocity_files = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/Iceland1/AVG/Iceland1_avg.??????????????.nc';
//config const boundary_files = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/Iceland1/AVG/Iceland1_bry.??????????????.nc';
//config const forcing_files  = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/Iceland1/AVG/Iceland1_bry.??????????????.nc';

config const maskfile       = '/glade/derecho/scratch/bachman/chapel_experiments/offline_BGC/remove_time/INPUTS_FULL/Iceland1_grd.zarr/mask/';
config const hfile          = '/glade/derecho/scratch/bachman/chapel_experiments/offline_BGC/remove_time/INPUTS_FULL/Iceland1_grd.zarr/h/';
config const velocity_files = '/glade/derecho/scratch/bachman/chapel_experiments/offline_BGC/remove_time/INPUTS_FULL/Iceland1_avg.??????????????.zarr';
//config const velocity_files = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/Iceland1/AVG/Iceland1_avg.??????????????.nc';
//config const boundary_files = '/glade/derecho/scratch/bachman/chapel_experiments/offline_BGC/remove_time/INPUTS_FULL/Iceland1_avg.??????????????_bry.zarr';
config const boundary_files = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/Iceland1/AVG/Iceland1_bry.??????????????.nc';
config const forcing_files  = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/Iceland1/AVG/Iceland1_bry.??????????????.nc';

config const Nx = 66;
config const Ny = 34;
config const Nz = 100;

/* Sigma coordinate parameters */
config const theta_s : real = 5.0;
config const theta_b : real = 2.0;
config const hc      : real = 300.0;

config const dx : real = 4000;
config const dy : real = 4000;
const area = dx * dy;
const iarea = 1.0 / area;

config const dt : real = 60.0;

/* Timestepping */  // do these need to be config vars?
config const Nt_start : int = 0;
config const Nt : int = 100;

// For LF-AM3 scheme
config const gamma = 0.0833333333333;
config const us = 0.16666666666666;

// For AB3 scheme
config const beta = 5.0/12.0;

// For RK4 scheme
config const one_sixth = 1.0 / 6.0;

// For PPM scheme
config const one_third = 1.0 / 3.0;

// For sponge
config const v_sponge : real = 300;
config const sponge_width : real = 15;

// Order of polynomial for boundary value extrapolation
config const ord : int = 3;
