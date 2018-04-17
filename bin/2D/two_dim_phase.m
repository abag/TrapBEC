
function two_dim_phase(filenumber)
global xx yy Nx Ny Nz Psi
%for i=filenumber
    gather(filenumber)
    %colormap(hot)
    for j=1:Ny
        for k=1:Nx
            phase(j,k)=angle(Psi(ceil(Nz/2),j,k));
        end
    end
    colormap(fireprint)
    imagesc(xx,yy,phase) 
    set(gca,'YDir','normal');
    axis square
    colorbar
    set(gca,'FontSize',16)
    ylabel('$y$','Interpreter','LaTex','FontSize',20','rot',0);
    xlabel('$x$','Interpreter','LaTex','FontSize',20);
    c=colorbar;
    ylabel(c,'$\theta$','Interpreter','LaTex','FontSize',20','rot',0)
    daspect([1 1 1])
%end