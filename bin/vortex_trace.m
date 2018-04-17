global Nx Ny Nz Psi

[field1x field1y field1z]=gradient(real(Psi),0.3265);
[field2x field2y field2z]=gradient(imag(Psi),0.3265);
field3(Nz,Ny,Nx,1:3)=0.;
for k=1:Nz
    k
    for j=1:Ny
        for i=1:Nx
            vect1=[field1x(k,j,i) field1y(k,j,i) field1z(k,j,i)];
            vect2=[field2x(k,j,i) field2y(k,j,i) field2z(k,j,i)];
            field3(k,j,i,:)=cross(vect1,vect2);
        end
    end
end