module finite_diff
  use cdata
  use ghost
  use deriv
  use mpi_var
  contains
  !-----------------------------------
  subroutine get_rhs()
  implicit none
  integer :: i , j , k
  real :: pot
  call ghost_comm
  do i=1, nmeshx
    do j=1, nmeshy
      do k=1, nmeshz
      !Laplacian
      rhs(k,j,i)=second_deriv(Psi(k,j,(i-3):(i+3)),hx)
      rhs(k,j,i)=rhs(k,j,i)+second_deriv(Psi(k,(j-3):(j+3),i),hy)
      if (two_dim.eqv..false.) then
        rhs(k,j,i)=rhs(k,j,i)+second_deriv(Psi((k-3):(k+3),j,i),hz)
      end if
      !evaluate potential
      pot=0.
      select case(potential)
        case('zero')
          !do nothing but allow the possibility
        case('harmonic')
          pot=harm_wx*((xx(i))/Lx)**2+harm_wy*((yy(j))/Ly)**2+harm_wz*((zz(k))/Lz)**2
      end select
      if (t>potential_off) pot=0.
      rhs(k,j,i)=rhs(k,j,i)+(1.-abs(Psi(k,j,i)**2)-pot)*Psi(k,j,i)
      end do
    end do
  end do
  end subroutine
end module
