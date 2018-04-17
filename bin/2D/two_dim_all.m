function two_dim_all(filenumbers)
global xx yy Nx Ny Nz Psi
  for iloop=filenumbers;
      
      close all
  subplot(2,2,1)
    two_dim_plot(iloop)
  subplot(2,2,2)
    two_dim_phase(iloop)
  subplot(2,2,3)
    two_dim_phase2(iloop)
  if Nz==1
  subplot(2,2,4)
    vel_2D
  end
  fout=sprintf('./plot%03d.jpeg',iloop)
  print('-djpeg',fout)
  end