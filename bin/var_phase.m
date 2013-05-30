global Nx Ny Nz Psi
for k=1:Nz
    for j=1:Ny
        for i=1:Nx
            Z(k,j,i)=phase(Psi(k,j,i));
        end
    end
end
pcolor(squeeze(abs(Z(Nz/2,:,:))))
shading interp
colorbar
set(gca,'FontSize',16)
xlabel('x','FontSize',16)
ylabel('y','FontSize',16)