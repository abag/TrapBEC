% Outputs compressible and incompressible kinetic energy density
% Now assumes we are working in the case of box only
% AND 3d.
function [varargout] = kineticenergydecomposition3d(psik,dx)
% Input arguments: 
% psik - wavefunction in k-space. 
% dx - spacing of x, y & z 
% assumes dx = dy = dz and same gridsize in x, y and z

gridDim = size(psik);
% get the gridsize
mm = gridDim(1);


%--------------------------------------------------------------------------
% A: Calculating Velocities
%--------------------------------------------------------------------------

% spacing in kspace
dk = 2*pi/(dx*mm);

% setting up kx ky and kz
ki = [linspace(0,(mm/2-1)*dk,mm/2) linspace(-mm/2*dk,-dk,mm/2)];
kj = [linspace(0,(mm/2-1)*dk,mm/2) linspace(-mm/2*dk,-dk,mm/2)];
kz = [linspace(0,(mm/2-1)*dk,mm/2) linspace(-mm/2*dk,-dk,mm/2)];

% Get the wavefunction and its complex conjugate in real-space
psi2 = ifftn(psik);
conjpsi2 = conj(psi2);

% construct density 
dens2 = psi2.*conjpsi2;

% complex conjugate wfn in k-space
psiconjk = (fftn(conjpsi2));

% initialize
kxpsi= zeros(mm,mm,mm);
kxpsiconj = zeros(mm,mm,mm);
kypsi = zeros(mm,mm,mm);
kypsiconj = zeros(mm,mm,mm);
kzpsi = zeros(mm,mm,mm);
kzpsiconj = zeros(mm,mm,mm);

% calculate derivatives
for ii = 1:1:length(ki)
  for jj = 1:1:length(kj)
      for kk = 1:1:length(kz)
          kxpsi(ii,jj,kk)=(1i*ki(ii))*psik(ii,jj,kk);
          kxpsiconj(ii,jj,kk)=(1i*ki(ii))*psiconjk(ii,jj,kk);
          kypsi(ii,jj,kk)=(1i*kj(jj))*psik(ii,jj,kk);
          kypsiconj(ii,jj,kk)=(1i*kj(jj))*psiconjk(ii,jj,kk);
          kzpsi(ii,jj,kk)=(1i*kz(kk))*psik(ii,jj,kk);
          kzpsiconj(ii,jj,kk)=(1i*kz(kk))*psiconjk(ii,jj,kk);
      end
  end
end

dxpsi=(ifftn(kxpsi));
dypsi=(ifftn(kypsi));
dzpsi=(ifftn(kzpsi));
dxpsiconj=(ifftn(kxpsiconj));
dypsiconj=(ifftn(kypsiconj));
dzpsiconj=(ifftn(kzpsiconj));

clear kxpsi; clear kypsi; clear kzpsi;
clear kxpsiconj; clear kypsiconj; clear kzpsiconj;
clear psik; clear conjpsik; 

% calculate velocity components
velx = real(-0.5.*1i.*(conjpsi2(:,:,:).*dxpsi(:,:,:)-psi2(:,:,:).*dxpsiconj(:,:,:))./dens2(:,:,:));
vely = real(-0.5.*1i.*(conjpsi2(:,:,:).*dypsi(:,:,:)-psi2(:,:,:).*dypsiconj(:,:,:))./dens2(:,:,:));
velz = real(-0.5.*1i.*(conjpsi2(:,:,:).*dzpsi(:,:,:)-psi2(:,:,:).*dzpsiconj(:,:,:))./dens2(:,:,:));

clear conjpsi2; clear psi2; 

%--------------------------------------------------------------------------
% B: Calculating the incompressible and compressible components
%--------------------------------------------------------------------------

% construct density weighted velocity field
omegax = sqrt(dens2(:,:,:)).*(velx(:,:,:));
omegay = sqrt(dens2(:,:,:)).*(vely(:,:,:));
omegaz = sqrt(dens2(:,:,:)).*(velz(:,:,:));

clear dens2; 

% density weighted velocity field in k-space
omegax_kx = (fftn(omegax));
omegay_ky = (fftn(omegay));
omegaz_kz = (fftn(omegaz));

% initialize
komegac_kx = zeros(mm,mm,mm);
komegac_ky = zeros(mm,mm,mm);
komegac_kz = zeros(mm,mm,mm);

komegai_kx = zeros(mm,mm,mm);
komegai_ky = zeros(mm,mm,mm);
komegai_kz = zeros(mm,mm,mm);

absk = zeros(mm,mm);

% construct compressible and incompressible components
for ii = 1:1:length(ki)
  for jj = 1:1:length(kj)
    for kk = 1:1:length(kz)
        absk(ii,jj,kk) = ki(ii)*ki(ii)+kj(jj)*kj(jj)+kz(kk)*kz(kk);
        komegac_ky(ii,jj,kk) = (kj(jj)*(ki(ii)*omegax_kx(ii,jj,kk)+kz(kk)*omegaz_kz(ii,jj,kk))+kj(jj)*kj(jj)*omegay_ky(ii,jj,kk))/(absk(ii,jj,kk));
        komegac_kx(ii,jj,kk) = (ki(ii)*ki(ii)*omegax_kx(ii,jj,kk)+ki(ii)*(kj(jj)*omegay_ky(ii,jj,kk)+kz(kk)*omegaz_kz(ii,jj,kk)))/(absk(ii,jj,kk));
        komegac_kz(ii,jj,kk) = (kz(kk)*(ki(ii)*omegax_kx(ii,jj,kk)+kj(jj)*omegay_ky(ii,jj,kk))+kz(kk)*kz(kk)*omegaz_kz(ii,jj,kk))/(absk(ii,jj,kk));
        komegai_ky(ii,jj,kk) = omegay_ky(ii,jj,kk) - komegac_ky(ii,jj,kk);
        komegai_kx(ii,jj,kk) = omegax_kx(ii,jj,kk) - komegac_kx(ii,jj,kk);
        komegai_kz(ii,jj,kk) = omegaz_kz(ii,jj,kk) - komegac_kz(ii,jj,kk);
    end
  end
end

% deal with nan's... 
komegac_ky(find(isnan(komegac_ky))) = 0;
komegac_kx(find(isnan(komegac_kx))) = 0;
komegac_kz(find(isnan(komegac_kz))) = 0;
komegai_ky(find(isnan(komegai_ky))) = 0;
komegai_kx(find(isnan(komegai_kx))) = 0;
komegai_kz(find(isnan(komegai_kz))) = 0;

omegac_x = real(ifftn(komegac_kx));
omegac_y = real(ifftn(komegac_ky));
omegac_z = real(ifftn(komegac_kz));

omegai_x = real(ifftn(komegai_kx));
omegai_y = real(ifftn(komegai_ky));
omegai_z = real(ifftn(komegai_kz));

% compressible and incompressible energy components
Ekinsq_c = 0.5.*((omegac_x.^2+omegac_y.^2+omegac_z.^2));
Ekinsq_i = 0.5.*((omegai_x.^2+omegai_y.^2+omegai_z.^2));

if nargout == 2
 varargout = {Ekinsq_c, Ekinsq_i};
else
 varargout = {Ekinsq_i};
end
end