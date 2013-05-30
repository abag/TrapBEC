%function var_vor_iso
global Nx Ny Nz Psi Lx Ly Lz vortex_pos
whitebg('w')
temp_Psi=smooth3(Psi);
temp_Psi=permute(temp_Psi,[3 2 1]);
p = patch(isosurface(1-abs(temp_Psi.^2),.8));
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
hold on
size(vortex_pos)
plot3(vortex_pos(:,2),vortex_pos(:,1),vortex_pos(:,3),'o','MarkerFaceColor','b', 'MarkerSize',20)
%zoom(1.5)