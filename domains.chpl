use StencilDist;

use INPUTS;

// Creating a singleton first dimension to store the time
  const FullDomain2D = {0..<Ny, 0..<Nx};
  const FullDomain3D = {0..<Nz, 0..<Ny, 0..<Nx};
  const FullDomain2D_h = {0..<Ny, 0..<Nx};
  const FullDomain3D_h = {0..<Nz, 0..<Ny, 0..<Nx};
  const FullDomain2D_u = {0..<Ny, 0..<(Nx-1)};
  const FullDomain3D_u = {0..<Nz, 0..<Ny, 0..<(Nx-1)};
  const FullDomain2D_v = {0..<(Ny-1), 0..<Nx};
  const FullDomain3D_v = {0..<Nz, 0..<(Ny-1), 0..<Nx};
  const FullDomain_n = {0..<Nz, (Ny-1)..(Ny-1), 0..<Nx};
  const FullDomain_s = {0..<Nz, 0..0, 0..<Nx};
  const FullDomain_w = {0..<Nz, 0..<Ny, 0..0};
  const FullDomain_e = {0..<Nz, 0..<Ny, (Nx-1)..(Nx-1)};

  var myTargetLocales2D = reshape(Locales, {1..1, 1..Locales.size});
  var myTargetLocales3D = reshape(Locales, {1..1, 1..1, 1..Locales.size});

  const fluff2D = (0,3);
  const fluff3D = (0,0,3);

// Create stencilDist arrays for the tracers. This is necessary to do the halo updates efficiently.
  const stencil2D = new stencilDist(boundingBox=FullDomain2D, targetLocales=myTargetLocales2D, fluff=fluff2D);
  const stencil3D = new stencilDist(boundingBox=FullDomain3D, targetLocales=myTargetLocales3D, fluff=fluff3D);
  const stencil2D_h = new stencilDist(boundingBox=FullDomain2D_h, targetLocales=myTargetLocales2D, fluff=fluff2D);
  const stencil3D_h = new stencilDist(boundingBox=FullDomain3D_h, targetLocales=myTargetLocales3D, fluff=fluff3D);
  const stencil2D_u = new stencilDist(boundingBox=FullDomain2D_u, targetLocales=myTargetLocales2D, fluff=fluff2D);
  const stencil3D_u = new stencilDist(boundingBox=FullDomain3D_u, targetLocales=myTargetLocales3D, fluff=fluff3D);
  const stencil2D_v = new stencilDist(boundingBox=FullDomain2D_v, targetLocales=myTargetLocales2D, fluff=fluff2D);
  const stencil3D_v = new stencilDist(boundingBox=FullDomain3D_v, targetLocales=myTargetLocales3D, fluff=fluff3D);

  const D2 = stencil2D.createDomain(FullDomain2D);
  const D3 = stencil3D.createDomain(FullDomain3D);
  const D2_h = stencil2D_h.createDomain(FullDomain2D_h);
  const D3_h = stencil3D_h.createDomain(FullDomain3D_h);
  const D2_u = stencil2D_u.createDomain(FullDomain2D_u);
  const D3_u = stencil3D_u.createDomain(FullDomain3D_u);
  const D2_v = stencil2D_v.createDomain(FullDomain2D_v);
  const D3_v = stencil3D_v.createDomain(FullDomain3D_v);

