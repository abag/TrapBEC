load('./data/dims.log');
global Psi nproc nmeshx nmeshy nmeshz Nx Ny Nz xx yy zz Lx Ly Lz vel qpress vortex_pos
Nx=dims(1);
Ny=dims(2);
Nz=dims(3);
nprocx=dims(4);
nprocy=dims(5);
nmeshx=dims(6);
nmeshy=dims(7);
nmeshz=dims(8);
nproc=nprocx*nprocy;
Psi(Nz,Ny,Nx)=0.;
vel(Nz,Ny,Nx,3)=0.;
qpress(Nz,Ny,Nx,3)=0.;
xx(Nx)=0;
yy(Ny)=0;
zz(Nz)=0;
Lx=[]; Ly=[]; Lz=[];
clear dims ;
vortex_pos=[];
