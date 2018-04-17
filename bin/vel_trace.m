global Nx Ny Nz Psi Z
[velx vely velz]=gradient(Z,0.3265);
vel2=sqrt(velx.^2+vely.^2+velz.^2);
vel2=vel2.*(abs(Psi).^2);
velx=velx.*(abs(Psi).^2);
vely=vely.*(abs(Psi).^2);
velz=velz.*(abs(Psi).^2);
[sx sy sz] = meshgrid(40:20:160, 40:20:160, 40:20:200)
streamline(velx,vely,velz,sx,sy,sz)