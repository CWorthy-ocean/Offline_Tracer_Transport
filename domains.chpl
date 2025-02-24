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

  const IntDomain_3D = {1..<(Nx-1), 1..<(Ny-1), 0..<Nz};
  const IntDomain_ts = {1..<(Nx-1), 1..<(Ny-1), 1..num_ts_tracers, 0..<Nz};
  const IntDomain_marbl = {1..<(Nx-1), 1..<(Ny-1), 1..num_marbl_tracers, 0..<Nz};
  const IntDomain_other = {1..<(Nx-1), 1..<(Ny-1), 1..num_other_tracers, 0..<Nz};
  const IntDomain_u = {1..<(Nx-2), 0..<Ny, 0..<Nz};
  const IntDomain_v = {0..<Nx, 1..<(Ny-2), 0..<Nz};

  const EdgeDomain_v_n = {0..<Nx,(Ny-2)..(Ny-2), 0..<Nz};
  const EdgeDomain_v_s = {0..<Nx, 0..0, 0..<Nz};
  const EdgeDomain_u_e = {(Nx-2)..(Nx-2), 0..<Ny, 0..<Nz};
  const EdgeDomain_u_w = {0..0, 0..<Ny, 0..<Nz};

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

  const stencil_int_3D = new stencilDist(boundingBox=IntDomain_3D, targetLocales=myTargetLocales_3D, fluff=fluff_3D);
  const stencil_int_ts = new stencilDist(boundingBox=IntDomain_ts, targetLocales=myTargetLocales_tr, fluff=fluff_tr);
  const stencil_int_marbl = new stencilDist(boundingBox=IntDomain_marbl, targetLocales=myTargetLocales_tr, fluff=fluff_tr);
  const stencil_int_other = new stencilDist(boundingBox=IntDomain_other, targetLocales=myTargetLocales_tr, fluff=fluff_tr);
  const stencil_int_u = new stencilDist(boundingBox=IntDomain_u, targetLocales=myTargetLocales_u, fluff=fluff_u);
  const stencil_int_v = new stencilDist(boundingBox=IntDomain_v, targetLocales=myTargetLocales_v, fluff=fluff_v);

  const stencil_edge_v_n = new stencilDist(boundingBox=EdgeDomain_v_n, targetLocales=myTargetLocales_v, fluff=fluff_v);
  const stencil_edge_v_s = new stencilDist(boundingBox=EdgeDomain_v_s, targetLocales=myTargetLocales_v, fluff=fluff_v);
  const stencil_edge_u_e = new stencilDist(boundingBox=EdgeDomain_u_e, targetLocales=myTargetLocales_u, fluff=fluff_u);
  const stencil_edge_u_w = new stencilDist(boundingBox=EdgeDomain_u_w, targetLocales=myTargetLocales_u, fluff=fluff_u);



  const D2 = stencil_2D.createDomain(FullDomain_2D);
  const D3 = stencil_3D.createDomain(FullDomain_3D);
  const D3_ts = stencil_ts.createDomain(FullDomain_ts);
  const D3_marbl = stencil_marbl.createDomain(FullDomain_marbl);
  const D3_other = stencil_other.createDomain(FullDomain_other);
  const D3_u = stencil_u.createDomain(FullDomain_u);
  const D3_v = stencil_v.createDomain(FullDomain_v);

  const D3_int_3D = stencil_int_3D.createDomain(IntDomain_3D);
  const D3_int_ts = stencil_int_ts.createDomain(IntDomain_ts);
  const D3_int_marbl = stencil_int_marbl.createDomain(IntDomain_marbl);
  const D3_int_other = stencil_int_other.createDomain(IntDomain_other);
  const D3_int_u = stencil_int_u.createDomain(IntDomain_u);
  const D3_int_v = stencil_int_v.createDomain(IntDomain_v);

  const D3_edge_v_n = stencil_edge_v_n.createDomain(EdgeDomain_v_n);
  const D3_edge_v_s = stencil_edge_v_s.createDomain(EdgeDomain_v_s);
  const D3_edge_u_e = stencil_edge_u_e.createDomain(EdgeDomain_u_e);
  const D3_edge_u_w = stencil_edge_u_w.createDomain(EdgeDomain_u_w);
