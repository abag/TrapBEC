module output
  use cdata
  use mpi_var
  implicit none
  character (len=40) :: print_file !for printing to file
  character (len=40) :: print_save !for printing to file
  contains
  !---------------------------------------------
  subroutine initial_print
    if (rank==0) then
      write(*,*) 'TrapBEC code initialised'
      write(*,'(a,i4.4,a,i4.4,a,i4.4,a)') ' Running a ',Nx,'x',Ny,'x',Nz,' GPE simulation'
      write(*,'(a,f6.2,a,f6.2,a,f6.2,a,f6.2)') ' Dimensions: Lx=', Lx, ' Ly=',Ly,' Lz=',Lz
      write(*,'(a, f7.3)') ' Timestep dt:', dt
      write(*,*) 'Boundaries: ', boundaries
      write(*,*) 'Initial condition: ', init_cond
      write(*,*) 'Potential: ', potential
      select case(potential)
        case('zero')
        case('harmonic')
          write(*,'(a,f5.2,a,f5.2,a,f5.2)') ' trapping freqs. wx: ', harm_wx, ' wy: ', harm_wy, ' wz ', harm_wz
        case default
          call emergency_stop('error incorrect potential set')
      end select
      write(*,'(a,i6.6,a,i4.4)') ' Running for ', nsteps, ' timesteps, file print every ', shots
      write(*,'(a,i5.5,a)') ' Propagating in imaginary time for maximum of ', imag_nsteps, ' timesteps'
      if (potential_off<1E10) then
        write(*,*) 'potential switched off at t= ', potential_off
      end if
    end if
  end subroutine
  !---------------------------------------------
  subroutine dims_print
   implicit none
   open(unit=87,file='./data/dims.log')
     write(87,*) Nx
     write(87,*) Ny
     write(87,*) Nz
     write(87,*) nprocx
     write(87,*) nprocy
     write(87,*) nmeshx
     write(87,*) nmeshy
     write(87,*) nmeshz
   close(87)
  end subroutine
  !---------------------------------------------
  subroutine var_print
    implicit none
    write(unit=print_file,fmt="(a,i2.2,a,i3.3,a)")"./data/proc",rank,"/var",itime/shots, ".dat"
    open(unit=98,file=print_file,status='replace',form='unformatted')
      write(98) coords, Lx, Ly, Lz, t, Psi(1:nmeshz,1:nmeshy,1:nmeshx), vel
    close(98)
    write(unit=print_file,fmt="(a,i2.2,a,i2.2,a)")"./data/proc",rank,"/qpress",itime/shots, ".dat"
    open(unit=98,file=print_file,status='replace',form='unformatted')
      write(98) coords, qpressure
    close(98)
    write(unit=print_save,fmt="(a,i2.2,a)")"./data/proc",rank,"/var.dat"
    open(unit=98,file=print_save,status='replace',form='unformatted')
      write(98) itime+1
      write(98) t
      write(98) Psi(1:nmeshz,1:nmeshy,1:nmeshx)
    close(98)
  end subroutine
  !---------------------------------------------
  subroutine print_timing
   implicit none
   open(unit=87,file='./data/timing.log')
     write(87,*) Nx*Ny*Nz, nprocx*nprocy, time_per_timestep
   close(87)
  end subroutine
end module
