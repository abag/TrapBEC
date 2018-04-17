function two_dim_identify(filenumber)
global xx yy Nx Ny Psi
TFR=100.;
for i=filenumber
    gather(i)
    vcount=0;
    for j=1:Nx
        for k=1:Ny
            if abs(Psi(1,j,k))^2<0.2            
                if sqrt(xx(j)^2+yy(k)^2)<TFR
                  vcount=vcount+1;
                  vortex(vcount,1)=xx(k);
                  vortex(vcount,2)=yy(j);
                end
            end
        end
    end
    %colormap(hot)
    colormap(fireprint)
    imagesc(xx,yy,squeeze(abs(Psi(1,:,:)))) 
    axis square
    hold on
    plot(vortex(:,1),vortex(:,2),'o')
    fOUT=sprintf('./data/out2d%02d.jpeg',i)
    print('-djpeg',fOUT)
end