function two_dim_plot(filenumber)
global xx yy Psi Nz
for i=filenumber
    filenumber
    gather(i)
    %colormap(hot)
    colormap(fireprint)
    imagesc(xx,yy,interp2(squeeze(abs(Psi(ceil(Nz/2),:,:))),1))
    set(gca,'YDir','normal');
    daspect([1 1 1])
    colorbar
    set(gca,'FontSize',16)
    ylabel('$y$','Interpreter','LaTex','FontSize',20','rot',0);
    xlabel('$x$','Interpreter','LaTex','FontSize',20);
    c=colorbar;
    ylabel(c,'$|\Psi|$','Interpreter','LaTex','FontSize',20','rot',0)
    %caxis([0. 1.2])
    pause(0.1)
    fOUT=sprintf('./data/out2d%02d.jpeg',i)
    print('-djpeg',fOUT)
end