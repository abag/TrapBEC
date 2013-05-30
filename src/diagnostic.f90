module diagnostic
  use cdata
  use mpi_var
  implicit none
  real :: maxPsi, minPsi, meanPsi, intPsi2, intPsi2_old=1.
  contains
  subroutine Psi_info
    implicit none
    real :: maxPsi_loc, minPsi_loc, meanPsi_loc
    maxPsi_loc=maxval(abs(Psi(1:nmeshz,1:nmeshy,1:nmeshx)))
    call MPI_REDUCE(maxPsi_loc,maxPsi,1,MPI_REAL,MPI_MAX,0,comm2d,ierror)
    minPsi_loc=minval(abs(Psi(1:nmeshz,1:nmeshy,1:nmeshx)))
    call MPI_REDUCE(minPsi_loc,minPsi,1,MPI_REAL,MPI_MIN,0,comm2d,ierror)
    meanPsi_loc=sum(abs(Psi(1:nmeshz,1:nmeshy,1:nmeshx)))
    call MPI_REDUCE(meanPsi_loc,meanPsi,1,MPI_REAL,MPI_SUM,0,comm2d,ierror)
    meanPsi=meanPsi/(Nx*Ny*Nz)
  end subroutine
  !------------------------------------------------------------
  subroutine Psi_norm
    implicit none
    real :: intPsi2_loc
    intPsi2_loc=sum(abs(Psi(1:nmeshz,1:nmeshy,1:nmeshx))**2)
    call MPI_REDUCE(intPsi2_loc,intPsi2,1,MPI_REAL,MPI_SUM,0,comm2d,ierror)
    call MPI_BCAST(intPsi2, 1, MPI_REAL, 0,comm2d,ierror)
  end subroutine
  !------------------------------------------------------------
  subroutine get_phase
    use deriv
    implicit none
    integer :: i, j, k
    do i=-2, nmeshx+3
      do j=-2, nmeshy+3
        do k=-2, nmeshz+3
          phase(k,j,i)=atan2(aimag(Psi(k,j,i)),real(Psi(k,j,i)))
        end do
      end do
    end do
    !now get the velocity field
    do i=1, nmeshx
      do j=1, nmeshy
        do k=1, nmeshz
          vel(k,j,i,1)=first_deriv(phase(k,j,(i-3):(i+3)),hx)
          vel(k,j,i,2)=first_deriv(phase(k,(j-3):(j+3),i),hy)
          vel(k,j,i,3)=first_deriv(phase((k-3):(k+3),j,i),hz)
        end do
      end do
    end do
  end subroutine
  !------------------------------------------------------------
  subroutine get_qpressure
    use deriv
    implicit none
    real :: rho
    integer :: i, j, k
    !get the  density field
    do i=-2, nmeshx+3
      do j=-2, nmeshy+3
        do k=-2, nmeshz+3
          density(k,j,i)=abs(Psi(k,j,i))
        end do
      end do
    end do
    !take the laplacian
    do i=0, nmeshx+1
      do j=0, nmeshy+1
        do k=0, nmeshz+1
          laplace_rho(k,j,i)=real_second_deriv(density(k,j,(i-2):(i+2)),hx)+&
                             real_second_deriv(density(k,(j-2):(j+2),i),hy)+&
                             real_second_deriv(density((k-2):(k+2),j,i),hz)
        end do
      end do
    end do
    do i=1, nmeshx
      do j=1, nmeshy
        do k=1, nmeshz
          qpressure(k,j,i,1)=low_first_deriv(laplace_rho(k,j,(i-1):(i+1)),hx)
          qpressure(k,j,i,2)=low_first_deriv(laplace_rho(k,(j-1):(j+1),i),hy)
          qpressure(k,j,i,3)=low_first_deriv(laplace_rho((k-1):(k+1),j,i),hz)
          !qpressure(k,j,i,1)=laplace_rho(k,j,i)
          !qpressure(k,j,i,2:3)=0.
        end do
      end do
    end do
  end subroutine
end module
