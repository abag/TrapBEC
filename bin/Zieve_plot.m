function Zieve_plot(filenumber)
global xx yy Psi Nx Ny Nz
for i=filenumber
    filenumber
    gather(i)
    %colormap(hot)
    colormap(gray)
    imagesc(xx,yy,interp2(squeeze(abs(Psi(ceil(Nz/2),:,:))).^2,1))
    set(gca,'YDir','normal');
    daspect([1 1 1])
    colorbar
    set(gca,'FontSize',16)
    ylabel('$y$','Interpreter','LaTex','FontSize',20','rot',0);
    xlabel('$x$','Interpreter','LaTex','FontSize',20);
    c=colorbar;
    ylabel(c,'$|\Psi|^2$','Interpreter','LaTex','FontSize',20','rot',0)
    %caxis([0. 1.2])
    for j=1:Ny
        for k=1:Nx
            phase(j,k)=angle(Psi(1,j,k));
            rr=sqrt(xx(k).^2+yy(j).^2);
            if rr>50
                potential(j,k)=1.1;
            else
                potential(j,k)=0.;
            end
        end
    end
    [xlocs,ylocs,pol] = gpeget2dvort_homg(squeeze(abs(Psi)).^2,phase,xx,yy,potential);
    hold on
    g = gscatter(xlocs,ylocs,pol,['b','r'],['^','o'],10,'off');
    if(length(g)==1 && pol(1)==1)
        set(g(1), 'MarkerFaceColor', 'r')
        set(g(1),'Marker','o');
        set(g(1),'MarkerEdgeColor','none');
    end
    if(length(g)==1 && pol(1)==-1)
        set(g(1), 'MarkerFaceColor', 'b')
        set(g(1),'Marker','^');
        set(g(1),'MarkerEdgeColor','none');
    end
    if(length(g)>1)
        set(g(1),'MarkerEdgeColor','none');
        set(g(1), 'MarkerFaceColor', 'b')
        set(g(2),'MarkerEdgeColor','none');
        set(g(2), 'MarkerFaceColor', 'r')
    end
    hold off
    pause(0.1)
    fOUT=sprintf('./data/out2d%02d.jpeg',i)
    print('-djpeg',fOUT)
end