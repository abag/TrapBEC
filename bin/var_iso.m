function var_iso
global Nx Ny Nz Psi Lx Ly Lz 
whitebg('w')
temp_Psi=smooth3(Psi);
temp_Psi=permute(temp_Psi,[3 2 1]);
p = patch(isosurface(1-abs(temp_Psi.^2),.7));
isonormals(1-abs(Psi),p)
set(p,'FaceColor','g','EdgeColor','none');
camlight 
lighting gouraud
daspect([1,1,1])
view(3); 
alpha(0.3)
axis([ 1 Nx 1 Ny 1 Nz])
axis off
set(gca,'FontSize',16)
zoom(1.5)