function GeorgeFind(filenumber)
global xx yy Psi Nx Ny
%for i=filenumber
    filenumber
    gather(filenumber)
    for j=1:Ny
        for k=1:Nx
            phase(j,k)=angle(Psi(1,j,k));
            potential(j,k)=0.;
        end
    end
    [xlocs,ylocs,pol] = gpeget2dvort_homg(squeeze(abs(Psi)).^2,phase,xx,yy,potential)
    

    