function var_iso
global Nx Ny Nz Psi Lx Ly Lz
whitebg('w')
test=angle(Psi);
[velx vely velz]=gradient(test);
mom=(Psi.*(conj(Psi))); %.*sqrt(velx.^2+vely.^2+velz.^2);%%
temp_Psi=smooth3(Psi);
temp_Psi=permute(temp_Psi,[3 2 1]);
temp_Psi=abs(temp_Psi.^2);
p = patch(isosurface(1-temp_Psi,0.7));
isonormals(1-temp_Psi,p)
set(p,'FaceColor','b','EdgeColor','none');
camlight 
lighting gouraud
daspect([1,1,1])
view(3); 
alpha(0.8)
axis([ 1 Ny 1 Nx 1 Nz])
axis on
box on
set(gca,'FontSize',16)
camproj('perspective')
xlabel('$y$','FontSize',18,'Interpreter','latex')
ylabel('$x$','FontSize',18,'Interpreter','latex')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
set(gca,'zticklabel',[])