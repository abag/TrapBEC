global Nx Ny Nz Psi
for k=1:Nz
    for j=1:Ny
        for i=1:Nx
            %Z(k,j,i)=angle(Psi(k,j,i));
            Z(k,j,i)=atan2(imag(Psi(k,j,i)),real(Psi(k,j,i)));
        end
    end
end
temp_Psi=smooth3(Psi);
idx=find((1-abs(temp_Psi.^2))<0.6);
[rows cols pags] = ind2sub(size(temp_Psi),idx);
for i=1:size(rows,1)
   Z(rows(i),cols(i),pags(i))=0.;
end
whitebg('w')
Z=permute(Z,[3 2 1]);
p = patch(isosurface(Z,0.7*pi));
isonormals(Z,p)
set(p,'FaceColor','k','EdgeColor','none');
camlight 
lighting gouraud
daspect([1,1,1])
view(3); 
axis([ 1 Nx 1 Ny 1 Nz])
axis on
box on
set(gca,'FontSize',16)
camproj('perspective')
hold on
whitebg('w')
temp_Psi=permute(temp_Psi,[3 2 1]);
p = patch(isosurface(1-abs(temp_Psi.^2),.7));
isonormals(1-abs(Psi),p)
set(p,'FaceColor','g','EdgeColor','none');
%isocolors(Z,p);
%set(p,'FaceColor','interp','EdgeColor','none')
colormap('hsv')
camlight 
lighting gouraud
daspect([1,1,1])
view(3); 
alpha(.7)
axis([ 1 Nx 1 Ny 1 Nz])
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