global Nx Ny Nz Psi
for k=1:Nz
    for j=1:Ny
        for i=1:Nx
            %Z(k,j,i)=angle(Psi(k,j,i));
            Z(k,j,i)=atan2(imag(Psi(k,j,i)),real(Psi(k,j,i)));
        end
    end
end
temp_Psi=(Psi);
Z=permute(Z,[3 2 1]);
temp_Psi=permute(temp_Psi,[3 2 1]);
p = patch(isosurface(1-abs(temp_Psi.^2),.72));
isonormals(1-abs(temp_Psi.^2),p)
isocolors(Z,p);
set(p,'FaceColor','interp','EdgeColor','none')
colormap('jet')
camlight 
lighting gouraud
daspect([1,1,1])
view(3); 
alpha(1)
axis([ 1 Ny 1 Nx 1 Nz])
axis on
box on
set(gca,'FontSize',16)
camproj('perspective')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
set(gca,'zticklabel',[])
view([-90 90])