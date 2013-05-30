function var_plot
global Psi nproc nmeshx nmeshy nmeshz Nx Ny Nz xx yy zz Lx Ly Lz t
[x,y,z]=meshgrid(1:Nx,1:Ny,1:Nz);
x=x*Lx/Nx-Lx/2.;
y=y*Lx/Nx-Lx/2.;
z=z*Lx/Nx-Lx/2.;
xi=x*2.;
yi=y*2.;
zi=z*2.;
Psi2=interp3(x,y,z,Psi,xi,yi,zi);
Psi2(isnan(Psi2))=0.;
A1=real(Psi2);
A2=imag(Psi2);
B(1,:,:,:)=A1;
B(2,:,:,:)=A2;
size(B)
%slice(xi,yi,zi,abs(Psi2),0,0,0)
%shading interp
