module initial
  use cdata
  use statistics
  use mpi_var
  complex, allocatable :: tmp(:,:,:) !helper array 
  contains
  subroutine mesh_init
    implicit none
    integer :: i, j, k
    real :: vort_pos(2), vort_pos2(2) !position of vortex in xy-plane
    real :: r, theta, func !to insert a vortex
    character(len=40) :: restart_file
    logical :: can_restart_loc=.false.
    complex :: ring_p1, ring_p2

    allocate(Psi(-2:nmeshz+3,-2:nmeshy+3,-2:nmeshx+3))
    allocate(phase(-2:nmeshz+3,-2:nmeshy+3,-2:nmeshx+3))
    allocate(density(-2:nmeshz+3,-2:nmeshy+3,-2:nmeshx+3))
    allocate(laplace_rho(0:nmeshz+1,0:nmeshy+1,0:nmeshx+1))
    allocate(qpressure(1:nmeshz,1:nmeshy,1:nmeshx,3))
    allocate(vel(1:nmeshz,1:nmeshy,1:nmeshx,3))
    allocate(xx(1:nmeshx),yy(1:nmeshy),zz(1:nmeshz))
    allocate(rhs(1:nmeshz,1:nmeshy,1:nmeshx),w(1:nmeshz,1:nmeshy,1:nmeshx))
    allocate(tmp(1:nmeshz,1:nmeshy,1:nmeshx))
    hx=Lx/Nx ; hy=Lx/Ny ; hz=Lz/Nz
    Psi=(0.,0.) !zero for now
    !how about the timestep 
    if (dt>(1./8.)*(min(hx,hy,hz)**2)) then
      if (rank==0) print*, 'WARNING TIMESTEP IS PROBABLY TOO LARGE!'
    else
      if (rank==0) print*, 'dt is below upper-limit of', (1./8.)*(min(hx,hy,hz)**2)
    end if
    do k=1,nmeshz
      zz(k)=real((2*k-1)/(2.*nmeshz))*Lz-Lz/2.
    end do
    do j=1,nmeshy
      yy(j)=real((2*j-1)/(2.*nmeshy))*Ly/dims(2)+coords(2)*Ly/dims(2)-Ly/2.
    end do
    do i=1,nmeshx
      xx(i)=real((2*i-1)/(2.*nmeshx))*Lx/dims(1)+coords(1)*Lx/dims(1)-Lx/2.
    end do
    !initialise the wavefunction
    write(unit=restart_file,fmt="(a,i2.2,a)")"./data/proc",rank,"/var.dat"
    inquire(file=restart_file, exist=can_restart_loc)
    call MPI_REDUCE(can_restart_loc,can_restart,1,MPI_LOGICAL,MPI_LAND,0,comm2d,ierror)
    if (rank==0) then
      if (can_restart) then
        write(*,'(a,i3.3,a)'), ' attempting to restart from ', nprocs, ' processes, this may fail!'
      end if
    end if
    !let the other processes know if we can restart
    call MPI_BCAST(can_restart, 1, MPI_LOGICAL, 0,comm2d,ierror)
    if (can_restart) then
      call reload_data
      return ! go no further
    end if
    select case(init_cond)
    case('constant')
      do i=1, nmeshx
        do j=1, nmeshy
          do k=1, nmeshz
            Psi(k,j,i)=cmplx(1./sqrt(2.),1./sqrt(2.))
          end do
        end do
      end do
    case('random_phase')
      do i=1, nmeshx
        do j=1, nmeshy
          do k=1, nmeshz
            Psi(k,j,i)=exp(eye*runif(0.,2.*pi))
          end do
        end do
      end do
    case('crow')
      vort_pos(2)=-Ly/4. ; vort_pos2(2)=-Ly/4.
      do i=1, nmeshx
        do j=1, nmeshy
          do k=1, nmeshz
              vort_pos(1)=-4+cos(pi*zz(k)/Lz)
              r=sqrt((xx(i)-vort_pos(1))**2+(yy(j)-vort_pos(2))**2)
              theta=atan2(yy(j)-vort_pos(2),xx(i)-vort_pos(1))
              func=1.-exp(-0.7*r**1.15)
              Psi(k,j,i)=func*exp(eye*theta)
              vort_pos2(1)=4-cos(pi*zz(k)/Lz)
              r=sqrt((xx(i)-vort_pos2(1))**2+(yy(j)-vort_pos2(2))**2)
              theta=atan2(yy(j)-vort_pos2(2),xx(i)-vort_pos2(1))
              func=1.-exp(-0.7*r**1.15)
              tmp(k,j,i)=func*exp(-eye*theta)
              Psi(k,j,i)=Psi(k,j,i)*tmp(k,j,i)
          end do
        end do
      end do
    case('single_vortex')
      vort_pos(1)=0.
      vort_pos(2)=0.
      do i=1, nmeshx
        do j=1, nmeshy
          do k=1, nmeshz
              r=sqrt((xx(i)-vort_pos(1))**2+(yy(j)-vort_pos(2))**2)
              theta=atan2(yy(j)-vort_pos(2),xx(i)-vort_pos(1))
              func=1.-exp(-0.7*r**1.15)
              Psi(k,j,i)=func*exp(eye*theta)
          end do
        end do
      end do
    case('sausage')
      vort_pos(1)=0.
      vort_pos(2)=0.
      do i=1, nmeshx
        do j=1, nmeshy
          do k=1, nmeshz
              r=sqrt((xx(i)-vort_pos(1))**2+(yy(j)-vort_pos(2))**2)
              theta=atan2(yy(j)-vort_pos(2),xx(i)-vort_pos(1))
              if (abs(zz(k))<4.) then
                func=1.-exp(-0.1*r**1.15)
              else
                func=1.-exp(-0.7*r**1.15)
              end if
              Psi(k,j,i)=func*exp(eye*theta)
          end do
        end do
      end do
    case('random_points')
      call setup_random_points
    case('vortex_rings')
      call setup_vortex_rings
    case('line_of_lines')
      call setup_line_of_lines
    case('linked_rings')
      call setup_linked_rings
    case('lattice')
      call setup_lattice
    case('lattice2')
      call setup_lattice2
    case default
      call emergency_stop('wrong initial conditions set')
    end select
    deallocate(tmp)
  end subroutine
  !----------------------------------------------------------
  subroutine setup_vortex_rings
    implicit none
    real :: r, theta
    real :: TF_dense
    real :: anglex, angley, anglez
    real :: temp_x, temp_y, temp_z
    real :: translatex,translatey,translatez
    complex :: ring_p1, ring_p2
    integer :: i, j, k, v
      do i=1, nmeshx
        do j=1, nmeshy
          do k=1, nmeshz
            TF_dense=1.-(0.5*harm_wx*1.5/Lx**2)*(xx(i))**2 &
                         -(0.5*harm_wy*1.5/Ly**2)*(yy(j))**2 &
                         -(0.5*harm_wz*1.5/Lz**2)*(zz(k))**2
            TF_dense=TF_dense/sqrt(2.)
            if (TF_dense<0.) TF_dense=0.
            Psi(k,j,i)=cmplx(TF_dense,TF_dense)
            Psi(k,j,i)=(1.,1.)
          end do
        end do
      end do
      do v=1,vor_ring_count
      if (rank==0) then
        anglex=0.
        angley=runif(0.,2.*pi)
        anglez=runif(0.,2.*pi)
        translatex=runif(-Lx/4,Lx/4)
        translatey=runif(-Ly/4,Ly/4)
        translatez=runif(-Lz/4,Lz/4)
      end if
      !now broadcast these angles across procs
      call MPI_BCAST(anglex, 1, MPI_REAL, 0,comm2d,ierror)
      call MPI_BCAST(angley, 1, MPI_REAL, 0,comm2d,ierror)
      call MPI_BCAST(anglez, 1, MPI_REAL, 0,comm2d,ierror)
      call MPI_BCAST(translatex, 1, MPI_REAL, 0,comm2d,ierror)
      call MPI_BCAST(translatey, 1, MPI_REAL, 0,comm2d,ierror)
      call MPI_BCAST(translatez, 1, MPI_REAL, 0,comm2d,ierror)
      do i=1, nmeshx
        do j=1, nmeshy
          do k=1, nmeshz
              temp_x=cos(angley)*cos(anglez)*(xx(i))+&
                    (cos(anglex)*sin(anglez)+sin(anglex)*sin(angley)*cos(anglez))*(yy(j))+&
                    (sin(anglex)*sin(anglez)-cos(anglex)*sin(angley)*cos(anglez))*(zz(k))
              temp_y=-cos(angley)*sin(anglez)*(xx(i))+&
                     (cos(anglex)*cos(anglez)-sin(anglex)*sin(angley)*sin(anglez))*(yy(j))+&
                     (sin(anglex)*cos(angley)+cos(anglex)*sin(angley)*sin(anglez))*(zz(k))
              temp_z=sin(angley)*(xx(i))-sin(anglex)*cos(angley)*(yy(j))+cos(anglex)*cos(angley)*(zz(k))
              r=sqrt((temp_x-translatex)**2+(sqrt((temp_y-translatey)**2+(temp_z-translatez)**2)+ring_rad)**2)
              theta=atan2((temp_z-translatez),(temp_y-translatey))
              ring_p1=sqrt((r**2)*(0.347+0.0286*r**2)/(1+0.3333*r**2+0.0286*r**4))
              ring_p1=ring_p1*exp(eye*theta)
              r=sqrt((temp_x-translatex)**2+(sqrt((temp_y-translatey)**2+(temp_z-translatez)**2)-ring_rad)**2)
              ring_p2=sqrt((r**2)*(0.347+0.0286*r**2)/(1+0.3333*r**2+0.0286*r**4))
              ring_p2=ring_p2*exp(eye*theta)
              tmp(k,j,i)=ring_p1*conjg(ring_p2)
              Psi(k,j,i)=Psi(k,j,i)*tmp(k,j,i)
              if (isnan(real(Psi(k,j,i)))) then
                print*, i,j,k
              end if
              if (isnan(imag(Psi(k,j,i)))) then
                print*, i,j,k
              end if
          end do
        end do
      end do
      call MPI_BARRIER(comm2D,ierror)
      end do
  end subroutine
  !----------------------------------------------------------
  subroutine setup_random_points
    integer :: i, j, k
    integer :: v1
    real :: vort_pos(2) !position of vortex in xy-plane
    real :: r, theta, func !to insert a vortex
    real :: dir
      Psi(:,:,:)=(1.,1.)
      do v1=1,100
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
                func=1.-exp(-0.7*r**1.15)
                tmp(k,j,i)=func*exp(dir*eye*theta)
                Psi(k,j,i)=Psi(k,j,i)*tmp(k,j,i)
            end do
          end do
        end do
      end do
  end subroutine
  !----------------------------------------------------------
  subroutine setup_lattice
    integer :: i, j, k
    integer :: v1, v2
    real :: vort_pos(2) !position of vortex in xy-plane
    real :: r, theta, func !to insert a vortex
      Psi(:,:,:)=(1.,1.)
      do v1=1,3
      do v2=1,3
      vort_pos(1)=v1*(Lx/4.)
      vort_pos(2)=v2*(Ly/4.)
      do i=1, nmeshx
        do j=1, nmeshy
          do k=1, nmeshz
              r=sqrt((xx(i)-vort_pos(1))**2+(yy(j)-vort_pos(2))**2)
              theta=atan2(yy(j)-vort_pos(2),xx(i)-vort_pos(1))
              func=1.-exp(-0.7*r**1.15)
              tmp(k,j,i)=func*exp(eye*theta)
              Psi(k,j,i)=Psi(k,j,i)*tmp(k,j,i)
          end do
        end do
      end do
      end do ; end do
  end subroutine
  !----------------------------------------------------------
  subroutine setup_lattice2
    implicit none
    integer :: i, j, k
    integer :: v1, v2
    real :: vort_pos(2) !position of vortex in xy-plane
    real :: r, theta, func !to insert a vortex
      Psi(:,:,:)=(1.,1.)
      do v1=1,3
      do v2=1,3
      vort_pos(1)=v1*(Lx/4.)-Lx/2.
      vort_pos(2)=v2*(Ly/4.)-Ly/2.
      do i=1, nmeshx
        do j=1, nmeshy
          do k=1, nmeshz
              r=sqrt((xx(i)-vort_pos(1))**2+(yy(j)-vort_pos(2))**2)
              theta=atan2(yy(j)-vort_pos(2),xx(i)-vort_pos(1))
              func=1.-exp(-0.7*r**1.15)
              tmp(k,j,i)=func*exp(eye*theta)
              Psi(k,j,i)=Psi(k,j,i)*tmp(k,j,i)
          end do
        end do
      end do
      end do ; end do
      do v1=1,3
      do v2=1,3
      vort_pos(1)=v1*(Lx/4.)+10-Lx/2.
      vort_pos(2)=v2*(Lz/4.)+10-Lz/2.
      do i=1, nmeshx
        do j=1, nmeshy
          do k=1, nmeshz
              r=sqrt((xx(i)-vort_pos(1))**2+(zz(k)-vort_pos(2))**2)
              theta=atan2(zz(k)-vort_pos(2),xx(i)-vort_pos(1))
              func=1.-exp(-0.7*r**1.15)
              tmp(k,j,i)=func*exp(eye*theta)
              Psi(k,j,i)=Psi(k,j,i)*tmp(k,j,i)
          end do
        end do
      end do
      end do ; end do
  end subroutine
  !----------------------------------------------------------
  subroutine setup_line_of_lines
    implicit none
    integer :: i, j, k
    integer :: v1
    real :: translatex,translatey,translatez
    complex :: ring_p1, ring_p2
    real :: vort_pos(2) !position of vortex in xy-plane
    real :: r, theta, func !to insert a vortex
    Psi=cmplx(1/sqrt(2.),1/sqrt(2.))
      do v1=1,vor_line_count
      vort_pos(1)=-Lx/2.+(real(v1)/(vor_line_count+1))*Lx
      vort_pos(2)=0.
      do i=1, nmeshx
        do j=1, nmeshy
          do k=1, nmeshz
              r=sqrt((xx(i)-vort_pos(1))**2+(yy(j)-vort_pos(2))**2)
              theta=atan2(yy(j)-vort_pos(2),xx(i)-vort_pos(1))
              func=1.-exp(-0.7*r**1.15)
              if (mod(v1,2)==0) then
               tmp(k,j,i)=func*exp(eye*theta)
              else
               tmp(k,j,i)=func*exp(-eye*theta)
              end if
              Psi(k,j,i)=Psi(k,j,i)*tmp(k,j,i)
          end do
        end do
      end do
      end do
      translatex=-Lx/2.+5
      translatey=0.
      translatez=0.
      do i=1, nmeshx
        do j=1, nmeshy
          do k=1, nmeshz
              theta=atan2(xx(i)-translatex,sqrt(yy(j)**2+zz(k)**2)+ring_rad)
              r=sqrt((xx(i)-translatex)**2+(sqrt(yy(j)**2+zz(k)**2)+ring_rad)**2)
              ring_p1=ring_density(r)*exp(eye*theta)
              theta=atan2(xx(i)-translatex,sqrt(yy(j)**2+zz(k)**2)-ring_rad)
              r=sqrt((xx(i)-translatex)**2+(sqrt(yy(j)**2+zz(k)**2)-ring_rad)**2)
              ring_p2=ring_density(r)*exp(eye*theta)
              Psi(k,j,i)=Psi(k,j,i)*ring_p1*conjg(ring_p2)
          end do
        end do
      end do
  end subroutine
  !----------------------------------------------------------
  subroutine setup_linked_rings
    implicit none
    integer :: i, j, k
    integer :: v1
    real :: translatex,translatey,translatez
    complex :: ring_p1, ring_p2
    real :: vort_pos(2) !position of vortex in xy-plane
    real :: r, theta, func !to insert a vortex
      translatex=0.
      translatey=0.
      translatez=ring_rad*1.5
      do i=1, nmeshx
        do j=1, nmeshy
          do k=1, nmeshz
              theta=atan2(xx(i)-translatex,sqrt(yy(j)**2+zz(k)**2)+ring_rad)
              r=sqrt((xx(i)-translatex)**2+(sqrt(yy(j)**2+zz(k)**2)+ring_rad)**2)
              ring_p1=ring_density(r)*exp(eye*theta)
              theta=atan2(xx(i)-translatex,sqrt(yy(j)**2+zz(k)**2)-ring_rad)
              r=sqrt((xx(i)-translatex)**2+(sqrt(yy(j)**2+zz(k)**2)-ring_rad)**2)
              ring_p2=ring_density(r)*exp(eye*theta)
              Psi(k,j,i)=ring_p1*conjg(ring_p2)
          end do
        end do
      end do
      do i=1, nmeshx
        do j=1, nmeshy
          do k=1, nmeshz
              theta=atan2(yy(j)-translatey,sqrt(xx(i)**2+(zz(k)-translatez)**2)+ring_rad)
              r=sqrt((yy(j)-translatey)**2+(sqrt(xx(i)**2+(zz(k)-translatez)**2)+ring_rad)**2)
              ring_p1=ring_density(r)*exp(eye*theta)
              theta=atan2(yy(j)-translatey,sqrt(xx(i)**2+(zz(k)-translatez)**2)-ring_rad)
              r=sqrt((yy(j)-translatey)**2+(sqrt(xx(i)**2+(zz(k)-translatez)**2)-ring_rad)**2)
              ring_p2=ring_density(r)*exp(eye*theta)
              Psi(k,j,i)=Psi(k,j,i)*ring_p1*conjg(ring_p2)
          end do
        end do
      end do
  end subroutine
!----------------------------------------------------------
  real function ring_density(r)
    real, intent(IN) :: r
    real :: c1, c2, c3
    c1=11./32.
    c3=(5.-32.*c1)/(48.-192*c1)
    c2=c1*(c3-0.25)
    ring_density=sqrt((r**2)*(c1+c2*r**2)/(1+c3*r**2+c2*r**4))
  end function
  !----------------------------------------------------------
  subroutine reload_data
    implicit none
    character(len=40) :: restart_file
    integer :: itime_min, itime_max
    write(unit=restart_file,fmt="(a,i2.2,a)")"./data/proc",rank,"/var.dat"
    open(unit=98,file=restart_file,form='unformatted')
    read(98) itime
    read(98) t
    read(98) Psi(1:nmeshz,1:nmeshy,1:nmeshx)
    close(98)
    call MPI_REDUCE(itime,itime_min,1,MPI_INTEGER,MPI_MIN,0,comm2d,ierror)
    call MPI_REDUCE(itime,itime_max,1,MPI_INTEGER,MPI_MAX,0,comm2d,ierror)
    if (rank==0) then
      !check itime is consistent, on root proc
      if ((itime_min/=itime).or.(itime_max/=itime)) then
        !houston we have...
        call emergency_stop('itime not consistent in reloaded data')
      end if
    end if
    nstart=itime
    call MPI_BARRIER(comm2D,ierror)
    if (rank==0) write(*,*) 'succsesful restart, t=',t
  end subroutine
end module
