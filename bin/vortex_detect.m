function vortex_detect
global Psi nproc nmeshx nmeshy nmeshz Nx Ny Nz xx yy zz Lx Ly Lz t vortex_pos
%fix z, create slices in xy-plane
%loop over z
vect1(1:8)=0.; 
absPsi=abs(Psi);
full_counter=0;
modPsitol=0.3;
phasetol=0.7;
for k=20:Nz-20
   Pslice=squeeze(Psi(k,:,:));
   absPslice=squeeze(absPsi(k,:,:));
   for i=1:Nx ; for j=1:Ny
     Pphase(j,i)=phase(Pslice(j,i));
   end ; end
   %loop over and find points which are the local minimum
   for i=20:Nx-20
       for j=20:Ny-20
           if absPslice(j,i)<absPslice(j+1,i) && absPslice(j,i)<absPslice(j-1,i)...
               && absPslice(j,i)<absPslice(j,i+1) && absPslice(j,i)<absPslice(j,i-1)...
               && absPslice(j,i)<absPslice(j+1,i+1) && absPslice(j,i)<absPslice(j+1,i-1)...
               && absPslice(j,i)<absPslice(j-1,i+1) && absPslice(j,i)<absPslice(j-1,i-1)...
               && absPslice(j,i)<modPsitol
               %put local points into a vector
               count=0;
               for l1=-1:1 ; for l2=-1:1 
                 if l1==0 && l2==0 ; continue ; end
                 count=count+1;
                 vect1(count)=Pphase(j+l1,i+l2);
               end ; end
               if max(vect1)-min(vect1)>phasetol*2*pi
                   pos(j,i)=1;
                   full_counter=full_counter+1;
                   %vortex_pos(full_counter,1)=xx(i);
                   %vortex_pos(full_counter,2)=yy(j);
                   %vortex_pos(full_counter,3)=zz(k);
                   vortex_pos(full_counter,1)=i;
                   vortex_pos(full_counter,2)=j;
                   vortex_pos(full_counter,3)=k;
               end
           else
               pos(j,i)=0;
           end
       end
   end
   %figure
   %pcolor(pos)
   %colormap('jet')
   %figure
   %pcolor(Pphase)
   %hold on
   %contour(Pphase)
   %return
end
for j=20:Ny-20
   Pslice=squeeze(Psi(:,j,:));
   absPslice=squeeze(absPsi(:,j,:));
   for i=1:Nx ; for k=1:Nz
     Pphase(k,i)=phase(Pslice(k,i));
   end ; end
   %loop over and find points which are the local minimum
   for i=20:Nx-20
       for k=20:Nz-20
           if absPslice(k,i)<absPslice(k+1,i) && absPslice(k,i)<absPslice(k-1,i)...
               && absPslice(k,i)<absPslice(k,i+1) && absPslice(k,i)<absPslice(k,i-1)...
               && absPslice(k,i)<absPslice(k+1,i+1) && absPslice(k,i)<absPslice(k+1,i-1)...
               && absPslice(k,i)<absPslice(k-1,i+1) && absPslice(k,i)<absPslice(k-1,i-1)...
               && absPslice(k,i)<modPsitol
               %put local points into a vector
               count=0;
               for l1=-1:1 ; for l2=-1:1 
                 if l1==0 && l2==0 ; continue ; end
                 count=count+1;
                 vect1(count)=Pphase(k+l1,i+l2);
               end ; end
               if max(vect1)-min(vect1)>phasetol*2*pi
                   full_counter=full_counter+1;
                   %vortex_pos(full_counter,1)=xx(i);
                   %vortex_pos(full_counter,2)=yy(j);
                   %vortex_pos(full_counter,3)=zz(k);
                   vortex_pos(full_counter,1)=i;
                   vortex_pos(full_counter,2)=j;
                   vortex_pos(full_counter,3)=k;
               end
           end
       end
   end
end
for i=20:Nx-20
   Pslice=squeeze(Psi(:,:,i));
   absPslice=squeeze(absPsi(:,:,i));
   for j=1:Ny ; for k=1:Nz
     Pphase(k,j)=phase(Pslice(k,j));
   end ; end
   %loop over and find points which are the local minimum
   for j=20:Ny-20
       for k=20:Nz-20
           if absPslice(k,j)<absPslice(k+1,j) && absPslice(k,j)<absPslice(k-1,j)...
               && absPslice(k,j)<absPslice(k,j+1) && absPslice(k,j)<absPslice(k,j-1)...
               && absPslice(k,j)<absPslice(k+1,j+1) && absPslice(k,j)<absPslice(k+1,j-1)...
               && absPslice(k,j)<absPslice(k-1,j+1) && absPslice(k,j)<absPslice(k-1,j-1)...
               && absPslice(k,j)<modPsitol
               %put local points into a vector
               count=0;
               for l1=-1:1 ; for l2=-1:1 
                 if l1==0 && l2==0 ; continue ; end
                 count=count+1;
                 vect1(count)=Pphase(k+l1,j+l2);
               end ; end
               if max(vect1)-min(vect1)>phasetol*2*pi
                   full_counter=full_counter+1;
                   %vortex_pos(full_counter,1)=xx(i);
                   %vortex_pos(full_counter,2)=yy(j);
                   %vortex_pos(full_counter,3)=zz(k);
                   vortex_pos(full_counter,1)=i;
                   vortex_pos(full_counter,2)=j;
                   vortex_pos(full_counter,3)=k;
               end
           end
       end
   end
end