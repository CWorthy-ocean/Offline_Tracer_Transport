use INPUTS;
use StencilDist;

// Creating a singleton first dimension to store the time
  const FullDomain2D = {0..<Ny, 0..<Nx};
  const FullDomain3D = {0..<Nz, 0..<Ny, 0..<Nx};
  const FullDomain2D_h = {0..<Ny, 0..<Nx};
  const FullDomain3D_h = {0..<Nz, 0..<Ny, 0..<Nx};
  const FullDomain2D_u = {0..<Ny, 0..<(Nx-1)};
  const FullDomain3D_u = {0..<Nz, 0..<Ny, 0..<(Nx-1)};
  const FullDomain2D_v = {0..<(Ny-1), 0..<Nx};
  const FullDomain3D_v = {0..<Nz, 0..<(Ny-1), 0..<Nx};

  var myTargetLocales2D = reshape(Locales, {1..1, 1..Locales.size});
  var myTargetLocales3D = reshape(Locales, {1..1, 1..1, 1..Locales.size});

  const fluff2D = (0,2);
  const fluff3D = (0,0,2);

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
  const D2_u = stencil2D_h.createDomain(FullDomain2D_u);
  const D3_u = stencil3D_h.createDomain(FullDomain3D_u);
  const D2_v = stencil2D_h.createDomain(FullDomain2D_v);
  const D3_v = stencil3D_h.createDomain(FullDomain3D_v);

/*
class Domains {
  var grid : domain(2);
  var rho_2D : domain(2);
  var rho_3D : domain(3);
  var u_3D : domain(3);
  var v_3D : domain(3);
  var w_3D : domain(3);
}

proc set_domains(arg: Domains, const h_domain: ?, const u_domain: ?, const v_domain: ?, const w_domain: ?, const grid_domain: ?) {

  arg.rho_3D = h_domain.localSubdomain();

  arg.grid = {arg.rho_3D.dim[1], arg.rho_3D.dim[2]};

  arg.rho_2D = {arg.rho_3D.dim[1], arg.rho_3D.dim[2]};

//  arg.u_3D = u_domain.localSubdomain();
  arg.u_3D = get_u(arg.rho_3D);

//  arg.v_3D = v_domain.localSubdomain();
  arg.v_3D = get_v(arg.rho_3D);

//  arg.w_3D = w_domain.localSubdomain();
  arg.w_3D = get_w(arg.rho_3D);

}

proc get_u(ref D: domain(3)) {

  var D_u : domain(3);

  if (here.id == 0) {
    D_u = {D.dim[0], D.dim[1], 0..D.last[2]};
  }
  else if (here.id == (Locales.size-1)) {
    D_u = {D.dim[0], D.dim[1], (D.first[2]-1)..(D.last[2]-1)};
  }
  else {
    D_u = {D.dim[0], D.dim[1], (D.first[2]-1)..(D.last[2])};
  }

  return D_u;
}


proc get_v(ref D: domain(3)) {

  var D_v = {D.dim[0], 0..(D.last[1]-1), D.dim[2]};

  return D_v;
}

proc get_w(ref D: domain(3)) {

  var D_w = {-1..D.last[0], D.dim[1], D.dim[2]};

  return D_w;
}

*/
