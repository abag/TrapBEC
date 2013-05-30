function qpress_plot
global Psi qpress Nx Ny Nz xx yy zz Lx Ly Lz t
qpress2=sqrt(qpress(:,:,:,1).^2+qpress(:,:,:,2).^2+qpress(:,:,:,3).^2);
figure('Name',sprintf('slice in xy plane at z=%f, t=%f',Lz/2,t))
 pcolor(xx,yy,squeeze(qpress2(Nz/2,:,:)))
 shading interp
 xlabel('x','FontSize',16)
 ylabel('y','FontSize',16)
c=colorbar;
 set(gca,'FontSize',16)
ylabel(c,'|\Psi|','FontSize',16,'rot',0)
figure('Name',sprintf('slice in yz plane at x=%f, t=%f',Lx/2,t))
 pcolor(yy,zz,squeeze(qpress2(:,:,Nx/2)))
 shading interp
 xlabel('y','FontSize',16)
 ylabel('z','FontSize',16)
c=colorbar;
 set(gca,'FontSize',16)
ylabel(c,'|\Psi|','FontSize',16,'rot',0)
figure('Name',sprintf('slice in xz plane at y=%f, t=%f',Ly/2,t))
 pcolor(xx,yy,squeeze(qpress2(:,Ny/2,:)))
 shading interp
 xlabel('x','FontSize',16)
 ylabel('z','FontSize',16)
 c=colorbar;
 set(gca,'FontSize',16)
ylabel(c,'|\Psi|','FontSize',16,'rot',0)

