use StencilDist;

use INPUTS;

// Creating a singleton first dimension to store the time
  const FullDomain_2D = {0..<Nx, 0..<Ny};
  const FullDomain_3D = {0..<Nx, 0..<Ny, 0..<Nz};
  const FullDomain_ts = {0..<Nx, 0..<Ny, 1..num_ts_tracers, 0..<Nz};
  const FullDomain_marbl = {0..<Nx, 0..<Ny, 1..num_marbl_tracers, 0..<Nz};
  const FullDomain_other = {0..<Nx, 0..<Ny, 1..num_other_tracers, 0..<Nz};
  const FullDomain_u = {0..<(Nx-1), 0..<Ny, 0..<Nz};
  const FullDomain_v = {0..<Nx, 0..<(Ny-1), 0..<Nz};

  const FullDomain_n = {0..<Nx, (Ny-1)..(Ny-1), 0..<Nz};
  const FullDomain_s = {0..<Nx, 0..0, 0..<Nz};
  const FullDomain_w = {0..0, 0..<Ny, 0..<Nz};
  const FullDomain_e = {(Nx-1)..(Nx-1), 0..<Ny, 0..<Nz};

  var myTargetLocales_2D = reshape(Locales, {1..Locales.size, 1..1});
  var myTargetLocales_3D = reshape(Locales, {1..Locales.size, 1..1, 1..1});
  var myTargetLocales_tr = reshape(Locales, {1..Locales.size, 1..1, 1..1, 1..1});
  var myTargetLocales_u = reshape(Locales, {1..Locales.size, 1..1, 1..1});
  var myTargetLocales_v = reshape(Locales, {1..Locales.size, 1..1, 1..1});

  const fluff_2D = (3,0);
  const fluff_3D = (3,0,0);
  const fluff_tr = (3,0,0,0);
  const fluff_u = (3,0,0);
  const fluff_v = (3,0,0);

// Create stencilDist arrays for the tracers. This is necessary to do the halo updates efficiently.
  const stencil_2D = new stencilDist(boundingBox=FullDomain_2D, targetLocales=myTargetLocales_2D, fluff=fluff_2D);
  const stencil_3D = new stencilDist(boundingBox=FullDomain_3D, targetLocales=myTargetLocales_3D, fluff=fluff_3D);
  const stencil_ts = new stencilDist(boundingBox=FullDomain_ts, targetLocales=myTargetLocales_tr, fluff=fluff_tr);
  const stencil_marbl = new stencilDist(boundingBox=FullDomain_marbl, targetLocales=myTargetLocales_tr, fluff=fluff_tr);
  const stencil_other = new stencilDist(boundingBox=FullDomain_other, targetLocales=myTargetLocales_tr, fluff=fluff_tr);
  const stencil_u = new stencilDist(boundingBox=FullDomain_u, targetLocales=myTargetLocales_u, fluff=fluff_u);
  const stencil_v = new stencilDist(boundingBox=FullDomain_v, targetLocales=myTargetLocales_v, fluff=fluff_v);

  const D2 = stencil_2D.createDomain(FullDomain_2D);
  const D3 = stencil_3D.createDomain(FullDomain_3D);
  const D3_ts = stencil_ts.createDomain(FullDomain_ts);
  const D3_marbl = stencil_marbl.createDomain(FullDomain_marbl);
  const D3_other = stencil_other.createDomain(FullDomain_other);
  const D3_u = stencil_u.createDomain(FullDomain_u);
  const D3_v = stencil_v.createDomain(FullDomain_v);

