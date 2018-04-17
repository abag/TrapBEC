function var_iso
global Nx Ny Nz Psi Lx Ly Lz xx
harm_wx=7.5 ;  harm_wy=7.5 ;  harm_wz=7.5 ; 
for ii=1:Nx
    for jj=1:Ny
        for kk=1:Nz
            TF(ii,jj,kk)=1-harm_wx*((xx(ii))/Lx)^2-harm_wy*((xx(jj))/Ly)^2-harm_wz*((xx(kk))/Lz)^2;
            if TF(ii,jj,kk)<0
                TF(ii,jj,kk)=0;
            else
                TF(ii,jj,kk)=1.;
            end
        end
    end
end
whitebg('w')
test=angle(Psi);
test2=gradient(test);
size(test2)
temp_Psi=smooth3(Psi);
temp_Psi=permute(temp_Psi,[3 2 1]);
reldense=TF-Psi.*(conj(Psi));
reldense2=reldense(38:90,38:90,38:90);
p = patch(isosurface(reldense2,0.85));
isonormals(reldense2,p)
set(p,'FaceColor','g','EdgeColor','none');
camlight 
lighting gouraud
daspect([1,1,1])
view(3); 
alpha(0.4)
p.FaceColor = 'red';
p.EdgeColor = 'none';
%axis([ 1 Nx 1 Ny 1 Nz])
axis on
box on
set(gca,'FontSize',16)
camproj('perspective')
xlabel('$x$','FontSize',18,'Interpreter','latex')
ylabel('$y$','FontSize',18,'Interpreter','latex')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
set(gca,'zticklabel',[])