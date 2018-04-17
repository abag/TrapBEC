module cdata
  implicit none
  include 'cparam.local'
  include 'mpif.h'
  real :: hx, hy, hz !resolution in x, y, z
  real, allocatable :: xx(:),yy(:),zz(:) !local physical coords
  complex, allocatable :: Psi(:,:,:) !local wavefunction
  real, allocatable :: Trap(:,:,:) !trapping potential 
  complex, allocatable :: rhs(:,:,:), w(:,:,:) !local helpers
  real, allocatable :: phase(:,:,:), vel(:,:,:,:)
  real, allocatable :: density(:,:,:), laplace_rho(:,:,:), qpressure(:,:,:,:)
  integer :: nmeshz, nmeshy, nmeshx !local mesh sizes
  real, parameter :: pi=3.14159265359
  complex, parameter :: eye=(0.,1.)
  !--------------------------------------
  integer :: itime, nstart=1
  real :: t=0.
  logical :: can_restart=.false.
  real :: timing_t1, timing_t2, time_per_timestep=0.
  real :: comm_timing_t1, comm_timing_t2, tot_comm_time=0.
  !--------------------------------------
  real :: imag_Psiratio
end module
