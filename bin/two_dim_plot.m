function two_dim_plot(filenumber)
global xx yy Psi
for i=filenumber
    gather(i)
    %colormap(hot)
    colormap(fireprint)
    imagesc(xx,yy,squeeze(abs(Psi(1,:,:)))) 
    axis square
    fOUT=sprintf('./data/out2d%02d.jpeg',i)
    print('-djpeg',fOUT)
end