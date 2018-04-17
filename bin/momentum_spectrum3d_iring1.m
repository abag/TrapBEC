dens = Psi.*conj(Psi);
dens=Ekinsq_i;
dims = size(Psi);
mm = dims(1); 
L = -2*min(xx);
dx = xx(2)-xx(1);
dk = 2*pi/(dx*mm);

% plot condensate density;
%cd /data/acw/ringfigures 
% figure, close
% %figure(1); 
% %set(gcf,'Visible', 'off');
%figure(1); 
%isosurface(dens(:,:,:),0.0003); alpha('0.2'); axis equal; 
%disp('here1')
% %saveas(gcf, sprintf('densiso_iring1g8000fine_t%d.jpg',tt),'jpg'); 
% print(gcf,'-depsc','densiso_iring1.jpg')
% close(gcf);


% densx(:,:) =  dens(100,:,:);
% densy(:,:) =  dens(:,100,:);
% densz(:,:) =  dens(:,:,100);
% 
% figure(2); pcolor(densx); shading interp; axis equal; 
% figure(3); pcolor(densy); shading interp; axis equal; 
% figure(4); pcolor(densz); shading interp; axis equal;

% setting up grids for fourier space
ki = [linspace(0,(mm/2-1)*dk,mm/2) linspace(-mm/2*dk,-dk,mm/2)];
kj = [linspace(0,(mm/2-1)*dk,mm/2) linspace(-mm/2*dk,-dk,mm/2)];
kz = [linspace(0,(mm/2-1)*dk,mm/2) linspace(-mm/2*dk,-dk,mm/2)];
[kxgrid,kygrid,kzgrid]=ndgrid(ki,kj,kz);

% |k|
magk = sqrt(kxgrid.^2+kygrid.^2+kzgrid.^2);

% angle averaging k  
maxk = max(max(max(ki))); 

% \rho(kx,ky,kz)
densk = (fftn(dens));

mmkbin = floor(maxk/dk);  
kbin = 0:dk:(mmkbin+1)*dk;  
kind = length(kbin); 
 
%nkselect = zeros(1,kind); 
nkselectangleav = zeros(1,kind); 
nkselectcount = zeros(1,kind); 
 disp('here2')
for ii = 1:1:length(ki) 
    for jj = 1:1:length(kj) 
        for kk = 1:1:length(kz) 
            for ks = 1:1:kind-1 
                if (magk(ii,jj,kk)<kbin(ks+1) && magk(ii,jj,kk)>=kbin(ks)) 
                    nkselectangleav(ks) = nkselectangleav(ks) + abs(densk(ii,kk,jj)).*(magk(ii,jj,kk).*magk(ii,jj,kk)); 
                    nkselectcount(ks) = nkselectcount(ks)+1;  
                end % if 
            end 
        end 
    end 
end 
disp('here3')
nkangleav = nkselectangleav*4*pi; 
nkangleavcount = dx.*dx.*dx.*((nkangleav)./nkselectcount);

 
%figure, close 
%set(gcf,'Visible','off');
figure(2);
hold 
plot((log10(kbin)),log10(nkangleavcount),'-','LineWidth',1) 
hold 

%xlabel('Log_{10}(|k|)'); 
%ylabel('Log_{10}(n(k))');

%cd /data/acw/ringfigures
%print(gcf,'-depsc','momdist_iring1.jpg');

%saveas(gcf, sprintf('momdist_iring1_Radius1p8_t%d.fig',tt),'fig'); 
%close(gcf);

%save('momdist_iring15_Radius1p6.mat','nkangleavcount','kbin');%,'dens');

%exit; 
