global Nx Ny Nz Psi
for k=1:Nz
    for j=1:Ny
        for i=1:Nx
            %Z(k,j,i)=angle(Psi(k,j,i));
            Z(k,j,i)=atan2(imag(Psi(k,j,i)),real(Psi(k,j,i)));
        end
    end
end
imagesc(squeeze((Z(Nz/2,:,:))))
shading interp
colorbar
colormap(fireprint)
set(gca,'FontSize',16)


