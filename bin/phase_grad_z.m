function phase_grad_z(filenumber)
gather(filenumber)
global Nx Ny Nz Psi xx yy zz
Z=angle(Psi);
grad_phase=diff(squeeze(Z(:,40,64)));
for i=1:length(grad_phase)
    if abs(grad_phase(i))>2
        if i==1 || i==(Nz-1)
            grad_phase(i)=NaN;
        else
            grad_phase(i)=0.5*(grad_phase(i-1)+grad_phase(i+1));
        end
    end
end
subplot(2,1,1)
plot(zz(1:end),squeeze(Z(:,40,64)),'LineWidth',2)
hold on
subplot(2,1,2)
plot(zz(1:end-1),grad_phase,'LineWidth',2)
hold on
