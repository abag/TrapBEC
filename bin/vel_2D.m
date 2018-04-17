Z=squeeze(angle(Psi));
[gradpsix, gradpsiy]=gradient(squeeze(Psi));
[gradhatpsix, gradhatpsiy]=gradient(conj(squeeze(Psi)));
velx=conj(squeeze(Psi)).*(gradpsix)-(squeeze(Psi)).*(gradhatpsix);
vely=conj(squeeze(Psi)).*(gradpsiy)-(squeeze(Psi)).*(gradhatpsiy);
% for i=1:Nx
%     for j=1:Ny
%         if abs(Psi(1,j,i))>0.01
%           velx(j,i)=velx(j,i)/abs((Psi(1,j,i)))^2;
%           vely(j,i)=vely(j,i)/abs((Psi(1,j,i)))^2;
%         else
%           velx(j,i)=0.;
%           vely(j,i)=0.; 
%         end
%     end
% end
velx=-squeeze(velx)*1i;
vely=-squeeze(vely)*1i;
vel2=sqrt(velx(:,:).^2+vely(:,:).^2);
imagesc(xx,yy,(vel2))
set(gca,'YDir','normal');
shading interp
colorbar
hold on
skip=10
quiver(xx(1:skip:Nx),yy(1:skip:Ny),velx(1:skip:end,1:skip:end),vely(1:skip:end,1:skip:end),'w','LineWidth',2)
daspect([1 1 1])
set(gca,'FontSize',16)
ylabel('$y$','Interpreter','LaTex','FontSize',20','rot',0);
    xlabel('$x$','Interpreter','LaTex','FontSize',20);
    c=colorbar;
    ylabel(c,'$mv$','Interpreter','LaTex','FontSize',20','rot',0)

