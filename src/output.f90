module output
  use cdata
  use mpi_var
  implicit none
  character (len=40) :: print_file !for printing to file
  character (len=40) :: print_save !for printing to file
  contains
  !---------------------------------------------
  subroutine initial_print
    implicit none
    if (rank==0) then
      write(*,*) 'TrapBEC code initialised'
      if (two_dim) write(*,*) ''//achar(27)//'[31mWARNING: Running in 2D mode '//achar(27)//'[0m.'
      write(*,'(a,i4.4,a,i4.4,a,i4.4,a)') ' Running a ',Nx,'x',Ny,'x',Nz,' GPE simulation'
      write(*,'(a,f6.2,a,f6.2,a,f6.2,a,f6.2)') ' Dimensions: Lx=', Lx, ' Ly=',Ly,' Lz=',Lz
      write(*,'(a, f7.3)') ' Timestep dt:', dt
      write(*,*) 'Boundaries: ', boundaries
      write(*,*) 'Initial condition: ', init_cond
      write(*,*) 'Potential: ', potential
      select case(potential)
        case('zero','sphere','barrier','channel','mobius','box_trap')
        case('harmonic')
          write(*,'(a,f5.2,a,f5.2,a,f5.2)') ' trapping freqs. wx: ', harm_wx, ' wy: ', harm_wy, ' wz ', harm_wz
        case('harmonic_ring','harmonic_ring2')
          write(*,'(a,f5.2)') ' trapping freq: ', harm_wx
        case('quartic')
              write(*,'(a,f5.2,a,f5.2,a,f5.2)') ' trapping freqs. wx: ', harm_wx, ' wy: ', harm_wy, ' wz ', harm_wz
        case default
          call emergency_stop('error incorrect potential set')
      end select
      write(*,'(a,i6.6,a,i4.4)') ' Running for ', nsteps, ' timesteps, file print every ', shots
      write(*,'(a,i5.5,a)') ' Propagating in imaginary time for maximum of ', imag_nsteps, ' timesteps'
      if (potential_off<1E10) then
        write(*,*) 'potential switched off at t= ', potential_off
      end if
      write(*,*) 'Dissipation, gamma= ', gamma
      if (abs(advective_vx)>epsilon(0.)) then
         write(*,*) 'Moving frame, vx=', advective_vx
      else if (abs(omega_z)>epsilon(0.)) then
         write(*,*) 'Rotating frame, Omega=', omega_z
         if (t_spindown<1E10) then
           write(*,*) 'rotation switched off at t= ', t_spindown
         end if
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
  subroutine imaginary_print
    implicit none
    integer :: i, j
      write(unit=print_file,fmt="(a,i2.2,a,i4.4,a)")"./data/proc",rank,"/imag_var",itime/shots, ".dat"
      open(unit=98,file=print_file,status='replace',form='unformatted')
        write(98) coords, Lx, Ly, Lz, t, Psi(1:nmeshz,1:nmeshy,1:nmeshx)
      close(98)
   end subroutine
  !---------------------------------------------
  subroutine var_print
    implicit none
    integer :: i, j
    if (full_binary_output) then
      write(unit=print_file,fmt="(a,i2.2,a,i4.4,a)")"./data/proc",rank,"/var",itime/shots, ".dat"
      open(unit=98,file=print_file,status='replace',form='unformatted')
        write(98) coords, Lx, Ly, Lz, t, Psi(1:nmeshz,1:nmeshy,1:nmeshx)
      close(98)
    end if
    if (vel_print) then
      write(unit=print_file,fmt="(a,i2.2,a,i2.2,a)")"./data/proc",rank,"/vel",itime/shots, ".dat"
      open(unit=98,file=print_file,status='replace',form='unformatted')
        write(98) coords, vel
      close(98)
    end if
    if (qpress_print) then
      write(unit=print_file,fmt="(a,i2.2,a,i2.2,a)")"./data/proc",rank,"/qpress",itime/shots, ".dat"
      open(unit=98,file=print_file,status='replace',form='unformatted')
        write(98) coords, qpressure
      close(98)
    end if
    if (gnuplot_2D_print) then
      write(unit=print_file,fmt="(a,i2.2,a,i2.2,a)")"./data/proc",rank,"/2D_slice",itime/shots, ".log"
      open(unit=98,file=print_file,status='replace')
      do i=1, nmeshx
        do j=1, nmeshy
          write(98,*) xx(i), yy(j), abs(Psi(1,j,i))**2
        end do
        write(98,*)
      end do
      close(98)
    end if
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
     write(87,*) '---------communication timing-----------'
     write(87,*) tot_comm_time
   close(87)
  end subroutine
end module
