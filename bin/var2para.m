function var2para(filenumbers)
global Psi
for i=filenumbers
    gather(i)
    rho=abs(Psi.^2);
    rho=permute(rho,[3 2 1]);
    fOUT=sprintf('./data/var%04d.vtk',i);
    writevtk(rho,fOUT)
end