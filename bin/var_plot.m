function var_plot
global Psi nproc nmeshx nmeshy nmeshz Nx Ny Nz xx yy zz Lx Ly Lz t
figure('Name',sprintf('slice in xy plane at z=%f, t=%f',Lz/2,t))
 pcolor(xx,yy,squeeze(abs(Psi(Nz/2,:,:))))
 shading interp
 xlabel('x','FontSize',16)
 ylabel('y','FontSize',16)
c=colorbar;
 set(gca,'FontSize',16)
ylabel(c,'|\Psi|','FontSize',16,'rot',0)
daspect([1 1 1])
figure('Name',sprintf('slice in yz plane at x=%f, t=%f',Lx/2,t))
 pcolor(yy,zz,squeeze(abs(Psi(:,:,Nx/2))))
 shading interp
 xlabel('y','FontSize',16)
 ylabel('z','FontSize',16)
c=colorbar;
 set(gca,'FontSize',16)
ylabel(c,'|\Psi|','FontSize',16,'rot',0)
daspect([1 1 1])
figure('Name',sprintf('slice in xz plane at y=%f, t=%f',Ly/2,t))
 pcolor(xx,zz,squeeze(abs(Psi(:,Ny/2,:))))
 shading interp
 xlabel('x','FontSize',16)
 ylabel('z','FontSize',16)
 c=colorbar;
 set(gca,'FontSize',16)
ylabel(c,'|\Psi|','FontSize',16,'rot',0)
daspect([1 1 1])

figure('Name','slices')
xslice = [0.0]; yslice = 0.; zslice = [0.,0.];
temp_Psi=permute(Psi,[3 2 1]);
slice(xx,yy,zz,abs(temp_Psi),xslice,yslice,zslice)
%colormap hsv
shading interp
 xlabel('x','FontSize',16)
 ylabel('y','FontSize',16)
 zlabel('z','FontSize',16)
 c=colorbar;
 set(gca,'FontSize',16)
 daspect([1 1 1])

