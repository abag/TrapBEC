  module finite_diff
  use cdata
  use ghost
  use deriv
  use mpi_var
  contains
  !-----------------------------------
  subroutine get_rhs()
  implicit none
  integer :: i , j , k, l1
  real :: pot, r, r1,r2, ring_loc1, ring_loc2, pow
  real :: theta, phi
  real :: las_x, las_y
  real :: las_freq
  CALL CPU_TIME(comm_timing_t1)
  call ghost_comm
  CALL CPU_TIME(comm_timing_t2)
  tot_comm_time=tot_comm_time+comm_timing_t2-comm_timing_t1
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
        case('harmonic_ring')
          r=sqrt((xx(i)/Lx)**2+(yy(j)/Lx)**2)
          pot=harm_wx*(r-0.3)**2
          if ((abs(r-0.3)<0.01)) then
             pot=max(0.,5.-0.1*t)
           end if
        case('harmonic_ring2')
          !ring_loc1=0.5
          !r1=(abs(xx(i)/(1.3*Lx))**(10)+abs(yy(j)/(.5*Ly)-ring_loc1)**(10))**(1./10)
          !r2=(abs(xx(i)/(1.3*Lx))**(10)+abs(yy(j)/(.5*Ly)+ring_loc1)**(10))**(1./10)
          !pot=5-(5*exp(-4000*((r1-0.25)/10)**2)+5*exp(-4000*((r2-0.25)/10)**2));
          pow=6
          r1=abs(0.25**pow-0.15*(xx(i)/Lx)**pow-4*(yy(j)/Ly-0.2)**pow)**(1./pow)
          r2=abs(0.25**pow-0.15*(xx(i)/Lx)**pow-4*(yy(j)/Ly+0.2)**pow)**(1./pow)
          if (r1<0.2495) then
            r1=0.
          else
            r1=5.
          end if
          if (r2<0.2495) then
            r2=0.
          else
            r2=5.
          end if
        pot=r1*r2
          if (pot<0.) pot=0.
          if ((abs(yy(j)/Ly)<0.01).and.(pot<0.01)) then
            pot=max(0.,5.-0.1*t)
          end if
        case('barrier')
          pot=0
          if (abs(yy(j)/Ly)>0.3) pot=5
          if (abs(yy(j)/Ly)<0.01) then
            pot=max(0.,5.-0.1*t)
          end if
        case('box_trap')
          pot=10
          r=sqrt(xx(i)**2+yy(j)**2)
          if ((r<Lx/3).and.(abs(zz(k))<Lz/2.5)) then
            pot=0
          end if
          pot=pot+5*sin(0.1*t)*zz(k)/Lz
        case('channel')
          pot=0
          if ((abs(yy(j))>0.75*Ly/2).or.(abs(zz(k))>0.75*Lz/2)) then
            pot=5.
          end if
          do l1=-1,1
          if (abs(xx(i))<Lx/Nx) then
            if (abs(yy(j)-l1*Ly/5)<Ly/Ny) then
              pot=7.5
            end if
          end if
          end do
          do l1=-1,1
          if (abs(xx(i))<Lx/Nx) then
          if (abs(zz(k)-l1*Lz/5)<Ly/Ny) then
           pot=7.5
          end if
          end if
          end do
        case('quartic')
          r=sqrt(((xx(i)/Lx)**2)+(yy(j)/Ly)**2)
          pot=harm_wx*r**8
        case('sphere')
          theta=atan2(yy(j),xx(i))
          r1=sqrt((xx(i)/Lx)**2+(yy(j)/Ly)**2+(zz(k)/Lz)**2)
          phi=acos((zz(k)/Lz)/r1)
          if (r1<0.2+0.01*sin(20*phi)+0.01*cos(30*theta+0.03)) then
            pot=10
          else
            pot=0
          end if
          !r1=sqrt((xx(i)/Lx)**2+(yy(j)/Ly)**2)
          !if (r1>0.45) then
          !  pot=5
          !end if
        case('mobius')
          pot=Trap(k,j,i)
          if ((abs(yy(j))<1.).and.(abs(xx(i)-25.)<1.)) then
            pot=pot+10
          end if
      end select
      if (lbeam) then
         las_freq=0.0
         las_x=70*cos(las_freq*t)*(1.-sin(las_freq*t))
         las_y=30*sin(las_freq*t)*cos(las_freq*t)
         pot=pot+20*exp(-((xx(i)-las_x)/wbeam)**2-((yy(j)-las_y)/wbeam)**2)
      end if
      if (t>potential_off) pot=0.
      !advective potential
      if (abs(advective_vx)>epsilon(0.)) then
        rhs(k,j,i)=rhs(k,j,i)+2*eye*advective_vx*cmplx_first_deriv(Psi(k,j,(i-3):(i+3)),hx)
      else if ((abs(omega_z)>epsilon(0.)).and.(t<t_spindown)) then
        rhs(k,j,i)=rhs(k,j,i)+omega_z*eye*(xx(i)*cmplx_first_deriv(Psi(k,(j-3):(j+3),i),hy) &
        -yy(j)*cmplx_first_deriv(Psi(k,j,(i-3):(i+3)),hx))
      end if
      rhs(k,j,i)=rhs(k,j,i)+(1.-abs(Psi(k,j,i)**2)-pot)*Psi(k,j,i)
      end do
    end do
  end do
  end subroutine
end module
