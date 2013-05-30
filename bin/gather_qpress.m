function gather_qpress(var_number)
global qpress nproc nmeshx nmeshy nmeshz 
cluster=0 % if set to 1 means double precision endians
for i=0:(nproc-1)
  filename=sprintf('./data/proc%02d/qpress%02d.dat',i,var_number);
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
  C=fread(fid,[3*nmeshz*nmeshy*nmeshx],'float32');
  qpress(1:nmeshz,coordy*nmeshy+1:(coordy+1)*nmeshy,coordx*nmeshx+1:(coordx+1)*nmeshx,:)=reshape(C,nmeshz,nmeshy,nmeshx,3);
  endian2=fread(fid,1,'float32');
  fclose(fid);
end
fclose('all');

                   
