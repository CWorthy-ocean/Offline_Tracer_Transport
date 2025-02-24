use FileSystem;

/* INPUT PARAMETERS */

config const maskfile = '/pscratch/sd/s/sbachman/DYE/dye_grd.nc_perm';
config const hfile = '/pscratch/sd/s/sbachman/DYE/dye_grd.nc_perm';
config const velocity_files = '/pscratch/sd/s/sbachman/DYE/dye_rnd.??????????????.nc??_perm';
config const boundary_files = '/pscratch/sd/s/sbachman/DYE/dye_bry.??????????????.nc??_perm';
config const forcing_files  = '/glade/derecho/scratch/bachman/UCLA-ROMS/run/Iceland1/AVG/Iceland1_bry.??????????????.nc';

const num_ts_tracers = 0;
const num_marbl_tracers = 0;
const num_other_tracers = 1;

config const Nx = 66;
config const Ny = 34;
config const Nz = 100;

/* Restart? */
config const restart = 0;
config const restart_file = '/glade/derecho/scratch/bachman/chapel_experiments/offline_BGC/remove_time/tracer.0000003510.nc';
config const Nt_start : int = 0;
config const Nt : int = 100;

/* Sigma coordinate parameters */
config const theta_s : real = 5.0;
config const theta_b : real = 2.0;
config const hc      : real = 300.0;

config const dx : real = 4000;
config const dy : real = 4000;
const area = dx * dy;
const iarea = 1.0 / area;

config const dt : real = 60.0;

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

// Horizontal viscosity
config const grid_Pe : real = 10;

// Order of polynomial for boundary value extrapolation
config const ord : int = 3;

// A really small number;
config const eps : real = 1e-16;

// Output frequency (in timesteps);
config const output_freq : int = 6;

// Tracers namelists
var marbl_namelist = ["PO4", "NO3", "SiO3", "NH4", "Fe", "Lig", "O2", "DIC", "DIC_ALT_CO2",
                        "ALK", "ALK_ALT_CO2", "DOC", "DON", "DOP", "DOPr", "DONr", "DOCr",
                        "zooC", "spChl", "spC", "spP", "spFe", "spCaCO3", "diatChl", "diatC",
                        "diatP", "diatFe", "diatSi", "diazChl", "diazC", "diazP", "diazFe"];

var ts_namelist = ["temp", "salt"];

var other_namelist = ["dye"];
