function BECexpand(filenumbers)
global Psi nproc nmeshx nmeshy nmeshz Nx Ny Nz xx yy zz Lx Ly Lz t
iso=0;
cmap=jet(filenumbers(end));
cmap2=jet(filenumbers(end));
for i=filenumbers
gather(i)
figure(1)
plot(xx,abs(squeeze(Psi(Nz/2,Ny/2,:))),'color',cmap(i,:))
hold on
figure(2)
plot(yy,abs(squeeze(Psi(Nz/2,:,Nx/2))),'color',cmap2(i,:))
hold on
if iso==1
figure(3)
temp_Psi=smooth3(Psi);
temp_Psi=permute(temp_Psi,[3 2 1]);
p = patch(isosurface(1-abs(temp_Psi.^2),.8));
isonormals(1-abs(Psi),p)
set(p,'FaceColor',cmap(i,:),'EdgeColor','none');
camlight 
lighting gouraud
daspect([1,1,1])
view(3); 
alpha(0.3)
end
end
if iso==1
figure(3)
axis([ 1 Nx 1 Ny 1 Nz])
axis off
set(gca,'FontSize',16)
zoom(1.5)
end