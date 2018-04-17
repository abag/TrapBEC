global Nx Ny Nz Psi
for k=1:Nz
    for j=1:Ny
        for i=1:Nx
            %Z(k,j,i)=angle(Psi(k,j,i));
            Z(k,j,i)=atan2(imag(Psi(k,j,i)),real(Psi(k,j,i)));
        end
    end
end
if Nz==1
    imagesc(squeeze((Z(1,:,:))))
else
    imagesc(squeeze((Z(Nz/2,:,:))))
end
shading interp
colorbar
colormap(fireprint)
set(gca,'FontSize',16)
daspect([1 1 1])
figure 
xslice = [0.0]; yslice = 0.; zslice = [0.,0.];
temp_Z=permute(Z,[3 2 1]);
slice(xx,yy,zz,abs(temp_Z),xslice,yslice,zslice)
shading interp


