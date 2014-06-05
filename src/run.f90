program run
   use cdata
   use mpi_var
   use output
   use initial
   use ghost
   use deriv
   use diagnostic
   use finite_diff
   use statistics
   use force
   implicit none
   call mpi_setup !mpi_var.f90
   call dims_print !output.f90
   call initial_print !output.f90
   call mesh_init !init.f90
   call ghost_init !ghost.f90
   if (forcing) call initialise_tempory_force_array !force.f90
   if ((can_restart.eqv..false.).and.(imaginary_time.eqv..true.)) then
     do itime=1, imag_nsteps
        call get_rhs
        Psi(1:nmeshz,1:nmeshy,1:nmeshx)=Psi(1:nmeshz,1:nmeshy,1:nmeshx)+dt*0.5*rhs
        call Psi_norm ; imag_Psiratio=intPsi2/intPsi2_old ; intPsi2_old=intPsi2
        !Psi=Psi*(intPsi2_old/intPsi2)
        !could put something in here to leave if ratio is close to unity
        if (mod(itime,shots)==0.and.rank==0) then
           print*, itime/shots, imag_Psiratio
        end if
     end do
     if (rank==0) write(*,*) 'finished imaginary time propagation'
   end if
   call MPI_BARRIER(comm2D,ierror)
   do itime=nstart,nsteps
     !timestep using RK3 - low storage
     CALL CPU_TIME(timing_t1)
     call get_rhs !finite_diff.f90
     w=eye*0.5*dt*rhs
     Psi(1:nmeshz,1:nmeshy,1:nmeshx)=Psi(1:nmeshz,1:nmeshy,1:nmeshx)+(w/3.)
     call get_rhs
     w=(-2./3.)*w+eye*0.5*dt*rhs
     Psi(1:nmeshz,1:nmeshy,1:nmeshx)=Psi(1:nmeshz,1:nmeshy,1:nmeshx)+w
     call get_rhs
     w=-w+eye*0.5*dt*rhs
     Psi(1:nmeshz,1:nmeshy,1:nmeshx)=Psi(1:nmeshz,1:nmeshy,1:nmeshx)+0.5*w
     t=t+dt
     CALL CPU_TIME(timing_t2)
     time_per_timestep=time_per_timestep+timing_t2-timing_t1
     !-------------output---------
     if (mod(itime,shots)==0) then
       call Psi_info ; call Psi_norm
       if (rank==0) then
         print*, itime/shots, t, maxPsi, minPsi, meanPsi, intPsi2
       end if
       call ghost_comm !ghost.f90
       !call get_phase !diagnostic.f90
       !call get_qpressure !diagnostic.f90
       call var_print !output.f90
     end if
     !-------------forcing---------
     if (forcing) then
       if (mod(itime,force_freq)==0) then
         call force_random_points
       end if
     end if
     !-----------------------------


   end do !end itime loop
   time_per_timestep=time_per_timestep/nsteps
   call print_timing !output.f90
   call MPI_FINALIZE(ierror)
end program

