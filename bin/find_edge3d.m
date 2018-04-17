function [selecteddensity] = find_edge3d(perc,psi)
% input perc and psi
% output selecteddensity, a matrix which is zero outside the condensate
% edge and 1 inside the condensate edge. 
% note, we could also output sphere edge

%clear all

dims = size(psi);
%realphi_1 = real(psi);
%imagphi_1 = imag(psi);
denstot = psi.*conj(psi);

% set up the parameters

mm = dims(1);

lx = dims(1);
ly = dims(2);
lz = dims(3);

% 3D density slice
dens = denstot(:,:,:);

comp = lx*ly*lz; 

% find the maximum density

sum1 = max(dens);
sum2 = max(sum1);
densmax = max(sum2);

% define a cutoff density

denscut = perc*densmax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We want to determine the edge of a condensate with vorticies, as defined
% by a certain density value.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

boundaryleft=zeros(1,mm/2);
boundaryright=zeros(1,mm/2);

Arrayleft=zeros(mm,mm,1);
Arrayright=zeros(mm,mm,1);

selecteddensity=zeros(mm,mm,mm);

sphereedge=zeros(mm,mm,mm);

%left half [i,j,k] corresponds to [z,y,x]

for kk = 1:1:lz
for jj = 1:1:ly
for ii = 1:1:lx/2
    if dens(ii,jj,kk) <= denscut
   boundaryleft(ii)=comp; % just needs to be a number greater than the value of lx/2 here
    else
   boundaryleft(ii)=ii;
    end        
end
Arrayleft(jj,kk)=min(boundaryleft);
boundaryleft=zeros(1,mm/2);
for ii = 1:1:lx/2
    if (Arrayleft(jj,kk)~=(comp) && ii>=Arrayleft(jj,kk))
        selecteddensity(ii,jj,kk)=1;
    end
    if (Arrayleft(jj,kk)~=(comp) && ii==Arrayleft(jj,kk))
     sphereedge(ii,jj,kk)=2.0;
    end
end
end
end


% right half

for kk = 1:1:lz
for jj = 1:1:ly
    for ii = (lx/2+1):1:lx
        if dens(ii,jj,kk) <= denscut
            boundaryright(ii-lx/2)=0;
        else
            boundaryright(ii-lx/2)=ii;
        end
    end
Arrayright(jj,kk)= max(boundaryright);
boundaryright=zeros(1,mm/2);   
for ii = (lx/2+1):1:lx
    if (Arrayright(jj,kk)~=(0) && ii<=Arrayright(jj,kk))
        selecteddensity(ii,jj,kk)=1;
    end
    if (Arrayright(jj,kk)~=(0) && ii==Arrayright(jj,kk))
     sphereedge(ii,jj,kk)=2.0;
    end
end
end
end
%end % function

%figure(1); pcolor(selecteddensity(:,:,floor(dims(3)/2))); shading interp
%figure(2); pcolor(sphereedge(:,:,floor(dims(3)/2))); shading interp

%figure(3); pcolor(squeeze(selecteddensity(:,floor(dims(3)/2),:))); shading interp
%figure(4); pcolor(squeeze(sphereedge(:,floor(dims(3)/2),:))); shading interp

%figure(5); pcolor(squeeze(selecteddensity(floor(dims(3)/2),:,:))); shading interp
%figure(6); pcolor(squeeze(sphereedge(floor(dims(3)/2),:,:))); shading interp
end
