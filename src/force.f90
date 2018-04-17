module force
  use cdata
  use statistics
  use mpi_var
  complex, allocatable, private :: ftmp(:,:,:) !helper array
  contains
    subroutine initialise_tempory_force_array
      implicit none
      allocate(ftmp(1:nmeshz,1:nmeshy,1:nmeshx))
    end subroutine
    !----------------------------------------------------------------
    subroutine force_random_points
    implicit none
    integer :: i, j, k
    integer :: v1
    real :: vort_pos(2) !position of vortex in xy-plane
    real :: r, theta, func !to insert a vortex
    real :: dir
        if (rank==0) then
          vort_pos(1)=runif(-Lx/4.,Lx/4.)
          vort_pos(2)=runif(-Ly/4.,Ly/4.)
          dir=runif(0.,1.)
          if (dir<0.5) then
            dir=-1.
          else
            dir=1.
          end if
        end if
        !now broadcast these quantities across procs
        call MPI_BCAST(vort_pos, 2, MPI_REAL, 0,comm2d,ierror)
        call MPI_BCAST(dir, 1, MPI_REAL, 0,comm2d,ierror)
        do i=1, nmeshx
          do j=1, nmeshy
            do k=1, nmeshz
                r=sqrt((xx(i)-vort_pos(1))**2+(yy(j)-vort_pos(2))**2)
                theta=atan2(yy(j)-vort_pos(2),xx(i)-vort_pos(1))
                !func=1.-exp(-0.7*r**1.15)
                func=pade_density(r)
                ftmp(k,j,i)=func*exp(dir*eye*theta)
                Psi(k,j,i)=Psi(k,j,i)*ftmp(k,j,i)
            end do
          end do
        end do
  end subroutine
  real function pade_density(r)
    implicit none
    real, intent(IN) :: r
    real :: c1, c2, c3
    c1=11./32.
    c3=(5.-32.*c1)/(48.-192*c1)
    c2=c1*(c3-0.25)
    pade_density=sqrt((r**2)*(c1+c2*r**2)/(1+c3*r**2+c2*r**4))
  end function
end module
