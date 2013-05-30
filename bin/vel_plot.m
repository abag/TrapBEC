global Psi nproc nmeshx nmeshy nmeshz Nx Ny Nz xx yy zz Lx Ly Lz t vel
vel2=sqrt(vel(:,:,:,1).^2+vel(:,:,:,2).^2+vel(:,:,:,3).^2);
figure
xslice = [0.05*Lx,0.9*Lx]; yslice = [0.4*Ly,0.9*Ly]; zslice = [0.1*Lz,0.75*Lz];
slice(xx,yy,zz,vel2(:,:,:),xslice,yslice,zslice)
colormap hot
shading interp
colorbar
set(gca,'FontSize',16)
return
figure
temp_Psi=smooth3(Psi);
temp_Psi=permute(temp_Psi,[3 2 1]);
temp_vel2=smooth3(vel2);
temp_vel2=permute(temp_vel2,[3 2 1]);
func=abs(temp_Psi).*temp_vel2;
p = patch(isosurface(func,1.));
isonormals(func,p)
set(p,'FaceColor','m','EdgeColor','none');
camlight 
lighting gouraud
daspect([1,1,1])
view(3); 
alpha(0.3)
axis([ 1 Nx 1 Ny 1 Nz])
axis off
set(gca,'FontSize',16)
zoom(1.5)
