function var_1d
global Psi nproc nmeshx nmeshy nmeshz Nx Ny Nz xx yy zz Lx Ly Lz t
for i=Nx/2-15:Nx/2+15 
plot(xx,abs(squeeze(Psi(:,Ny/2,i))),'-o')
hold on
plot(xx,abs(squeeze(Psi(i,Ny/2,:))),'-ro')
plot(xx,abs(squeeze(Psi(Nz/2,:,i))),'-go')
plot(xx,abs(squeeze(Psi(:,i,Nx/2))),'-mo')
dum_x=[0 0]
dum_y=[0 1]
plot(dum_x,dum_y)
axis([-Lx/2 Lx/2 0 inf])
hold off
pause
end