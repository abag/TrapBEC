
function two_dim_phase(filenumber)
global xx yy Nx Ny Nz Psi
    gather(filenumber)
    %colormap(hot)
    for j=1:Ny
        for k=1:Nx
            if abs(Psi(1,j,k))>0.5
              phase(j,k)=angle(Psi(ceil(Nz/2),j,k));
            else
              phase(j,k)=NaN;
            end
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