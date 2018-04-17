function two_dim_plot(filenumber)
global xx yy zz Psi Nz
for i=filenumber
    gather(i)
    colormap(fireprint)
    xslice = [0.0]; yslice = 0.; zslice = [0.,0.];
    temp_Psi=permute(Psi,[3 2 1]);
    slice(xx,yy,zz,abs(temp_Psi).^2,xslice,yslice,zslice)
    daspect([1 1 1])
    colorbar
    set(gca,'FontSize',16)
    ylabel('$y$','Interpreter','LaTex','FontSize',20','rot',0);
    xlabel('$x$','Interpreter','LaTex','FontSize',20);
    c=colorbar;
    ylabel(c,'$|\Psi|$','Interpreter','LaTex','FontSize',20','rot',0)
    caxis([0. 1.2])
    shading interp
    axis tight
    pause(0.1)
    fOUT=sprintf('./data/out2d%02d.jpeg',i)
    print('-djpeg',fOUT)
end