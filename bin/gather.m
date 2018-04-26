function gather(var_number)
global Psi nproc nmeshx nmeshy nmeshz Nx Ny Nz xx yy zz Lx Ly Lz t
cluster=0 % if set to 1 means double precision endians
t(1:nproc)=0.;
for i=0:(nproc-1)
  filename=sprintf('./data/proc%03d/var%04d.dat',i,var_number);
  fid=fopen(filename,'r');
  if (fid<0)
      disp('file does not exist-exiting script')
      return
  end
  frewind(fid);
  if (cluster==1)
    endian1=fread(fid,1,'float64');
  else
    endian1=fread(fid,1,'float32');
  end
  coordx=fread(fid,1,'int');
  coordy=fread(fid,1,'int');
  Lx=fread(fid,1,'float32');
  Ly=fread(fid,1,'float32');
  Lz=fread(fid,1,'float32');
  dumt=fread(fid,1,'float32');
  t(i+1)=dumt;
  A=fread(fid,[2*nmeshz*nmeshy*nmeshx],'float32');
  B=reshape(A,2,nmeshz,nmeshy,nmeshx);
  Psi(1:nmeshz,coordy*nmeshy+1:(coordy+1)*nmeshy,coordx*nmeshx+1:(coordx+1)*nmeshx)=B(1,:,:,:)+sqrt(-1.)*B(2,:,:,:);
  endian2=fread(fid,1,'float32');
  fclose(fid);
end
fclose('all');
%sanity check - is t consistent?
if std(t)>0. ; disp('t not consistent across procs') ; end
t=[] ; t=dumt ; %if std=0 OK to take one value
for i=1:Nx
  xx(i)=Lx*((2*i-1)/(2*Nx))-Lx/2.;
end
for i=1:Ny
  yy(i)=Ly*((2*i-1)/(2*Ny))-Ly/2.;
end
for i=1:Nz
  zz(i)=Lz*((2*i-1)/(2*Nz))-Lz/2.;
end
