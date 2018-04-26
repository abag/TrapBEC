function var2para_TF(filenumbers)
global Nx Ny Nz Psi Lx Ly Lz xx
harm_wx=5 ;  harm_wy=5;  harm_wz=5 ; 
for ii=1:Nx
    for jj=1:Ny
        for kk=1:Nz
            TF(ii,jj,kk)=1-harm_wx*((xx(ii))/Lx)^2-harm_wy*((xx(jj))/Ly)^2-harm_wz*((xx(kk))/Lz)^2;
            if TF(ii,jj,kk)<0
                TF(ii,jj,kk)=0;
            end
        end
    end
end
for i=filenumbers
    gather(i)
    rho=(1-TF)+abs(Psi.^2);
    rho=permute(rho,[3 2 1]);
    fOUT=sprintf('./data/var_TF%04d.vtk',i);
    writevtk(rho,fOUT)
end