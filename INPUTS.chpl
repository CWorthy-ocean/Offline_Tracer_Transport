use FileSystem;

/* INPUT PARAMETERS */

config const maskfile = '/glade/derecho/scratch/bachman/roms_marbl_setup_assistant/cases/dye_Atlantic/INPUTS/dye_grd.nc_perm';
config const hfile = '/glade/derecho/scratch/bachman/roms_marbl_setup_assistant/cases/dye_Atlantic/INPUTS/dye_grd.nc_perm';
config const initial_file = '/glade/derecho/scratch/bachman/roms_marbl_setup_assistant/cases/dye_Atlantic32/INPUTS/spinup_rst32.20120301000000.nc00_perm';
config const velocity_files = '/glade/derecho/scratch/bachman/roms_marbl_setup_assistant/cases/dye_Atlantic32/run_OTT/OUTPUT/HIS/dye_his.??????????????.nc??_perm';
config const boundary_files = '/glade/derecho/scratch/bachman/roms_marbl_setup_assistant/cases/dye_Atlantic32/run_OTT/OUTPUT/HIS/dye_bry.??????????????.nc??_perm';
config const surface_forcing_files  = '/glade/derecho/scratch/bachman/chapel_experiments/offline_BGC/Offline_Tracer_Transport/BGC_inputs/surface_forcing.nc_perm';
config const interior_forcing_files  = '/glade/derecho/scratch/bachman/chapel_experiments/offline_BGC/Offline_Tracer_Transport/BGC_inputs/interior_forcing.nc_perm';

const num_ts_tracers = 0;
const num_marbl_tracers = 0;
const num_other_tracers = 1;

config const Nx = 1082;
config const Ny = 1082;
config const Nz = 100;

/* Restart? */
config var restart = false;
config const restart_file = '/glade/derecho/scratch/bachman/chapel_experiments/offline_BGC/remove_time/tracer.0000003510.nc';
config const Nt_start : int = 0;
config const Nt : int = 4;

/* Sigma coordinate parameters */
config const theta_s : real = 5.0;
config const theta_b : real = 2.0;
config const hc      : real = 300.0;

config const dx : real = 1000;
config const dy : real = 1000;
const area = dx * dy;
const iarea = 1.0 / area;

/* Timesteps for OTT model and inputs */
config const dt : real = 300.0;
config const input_dt : real = 3600;
const read_freq = (input_dt / dt) : int;

// For 3rd-order upstream advection scheme
config const one_sixth = 1.0 / 6.0;

// For PPM scheme
config const one_third = 1.0 / 3.0;

// For sponge
config const v_sponge : real = 100;
config const sponge_width : real = 15;

// Order of polynomial for boundary value extrapolation
config const ord : int = 3;

// A really small number;
config const eps : real = 1e-16;

// I/O frequency (in timesteps);
config const write_freq : int = 0;
config const restart_freq : int = 120;
config const report_freq : int = 100;

// Tracers namelists
var marbl_namelist = ["PO4", "NO3", "SiO3", "NH4", "Fe", "Lig", "O2", "DIC", "DIC_ALT_CO2",
                        "ALK", "ALK_ALT_CO2", "DOC", "DON", "DOP", "DOPr", "DONr", "DOCr",
                        "zooC", "spChl", "spC", "spP", "spFe", "spCaCO3", "diatChl", "diatC",
                        "diatP", "diatFe", "diatSi", "diazChl", "diazC", "diazP", "diazFe"];

var ts_namelist = ["temp", "salt"];

var other_namelist = ["dye1", "dye2", "dye3", "dye4", "dye5", "dye6", "dye7", "dye8",
                      "dye9", "dye10", "dye11", "dye12", "dye13", "dye14", "dye15", "dye16",
                      "dye17", "dye18", "dye19", "dye20", "dye21", "dye22", "dye23", "dye24",
                      "dye25", "dye26", "dye27", "dye28", "dye29", "dye30", "dye31", "dye32"];
var numParSubcols = 1;
var numElementsSurfaceFlux = 5;
