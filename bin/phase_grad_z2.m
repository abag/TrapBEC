function phase_grad_z2(filenumber)
for fi=filenumber
    gather(fi)
    global Nx Ny Nz Psi xx yy zz
    Z=angle(Psi);
    for i=1:Ny
        for j=1:Nx
            vel(:,i,j)=diff(squeeze(Z(:,i,j)));
            for k=1:Nz-1
                if abs(vel(k,i,j))>4
                    if k==1 || k==(Nz-1)
                        vel(k,i,j)=NaN;
                    else
                        vel(k,i,j)=0.5*(vel(k-1,i,j)+vel(k+1,i,j));
                    end
                end
            end
        end
    end
    A=nanmean(vel(:,:,[1:Nx/2-4,(Nx/2+4):Nx]),3);
    B=nanmean(A(:,[1:Ny/2-4,(Ny/2+4):Ny]),2);
    plot(zz(1:end-1),B,'LineWidth',2)
    xlabel('$z$','FontSize',20,'Interpreter','LaTex')
    ylabel('$<u_z>$','FontSize',20,'Interpreter','LaTex')
    set(gca,'FontSize',16)
    axis([-70 70 -0.05 0.15])
    fOUT=sprintf('./data/uz%03d.jpeg',fi)
    print('-djpeg',fOUT)
end