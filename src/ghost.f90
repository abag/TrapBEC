module ghost
  use cdata
  use mpi_var
  implicit none
  complex, allocatable :: ghost_l1(:,:), ghost_l2(:,:), ghost_l3(:,:) 
  complex, allocatable :: ghost_l1_rec(:,:), ghost_l2_rec(:,:), ghost_l3_rec(:,:)
  complex, allocatable :: ghost_r1(:,:), ghost_r2(:,:), ghost_r3(:,:)
  complex, allocatable :: ghost_r1_rec(:,:), ghost_r2_rec(:,:), ghost_r3_rec(:,:)
  complex, allocatable :: ghost_u1(:,:), ghost_u2(:,:), ghost_u3(:,:)
  complex, allocatable :: ghost_u1_rec(:,:), ghost_u2_rec(:,:), ghost_u3_rec(:,:)
  complex, allocatable :: ghost_d1(:,:), ghost_d2(:,:), ghost_d3(:,:)
  complex, allocatable :: ghost_d1_rec(:,:), ghost_d2_rec(:,:), ghost_d3_rec(:,:)
  integer :: lsend(6), rsend(6), lrec(6), rrec(6)
  integer :: usend(6), dsend(6), urec(6), drec(6)
  contains
  subroutine ghost_init
    implicit none
     allocate(ghost_l1(nmeshz,nmeshy),ghost_l2(nmeshz,nmeshy),ghost_l3(nmeshz,nmeshy))
     allocate(ghost_l1_rec(nmeshz,nmeshy),ghost_l2_rec(nmeshz,nmeshy),ghost_l3_rec(nmeshz,nmeshy))
     allocate(ghost_r1(nmeshz,nmeshy),ghost_r2(nmeshz,nmeshy),ghost_r3(nmeshz,nmeshy))
     allocate(ghost_r1_rec(nmeshz,nmeshy),ghost_r2_rec(nmeshz,nmeshy),ghost_r3_rec(nmeshz,nmeshy))
     allocate(ghost_u1(nmeshz,nmeshx),ghost_u2(nmeshz,nmeshx),ghost_u3(nmeshz,nmeshx))
     allocate(ghost_u1_rec(nmeshz,nmeshx),ghost_u2_rec(nmeshz,nmeshx),ghost_u3_rec(nmeshz,nmeshx))
     allocate(ghost_d1(nmeshz,nmeshx),ghost_d2(nmeshz,nmeshx),ghost_d3(nmeshz,nmeshx))
     allocate(ghost_d1_rec(nmeshz,nmeshx),ghost_d2_rec(nmeshz,nmeshx),ghost_d3_rec(nmeshz,nmeshx))
  end subroutine
  subroutine ghost_comm
    implicit none
    if (left/=MPI_PROC_NULL) then
      !----------------send to the left-----------------
      ghost_l1(1:nmeshz,1:nmeshy)=Psi(1:nmeshz,1:nmeshy,1)
      ghost_l2(1:nmeshz,1:nmeshy)=Psi(1:nmeshz,1:nmeshy,2)
      ghost_l3(1:nmeshz,1:nmeshy)=Psi(1:nmeshz,1:nmeshy,3)
      call MPI_ISEND(ghost_l1,nmeshz*nmeshy,MPI_COMPLEX,left,tag(1),comm2D,lsend(1),ierror)
      call MPI_ISEND(ghost_l2,nmeshz*nmeshy,MPI_COMPLEX,left,tag(2),comm2D,lsend(2),ierror)
      call MPI_ISEND(ghost_l3,nmeshz*nmeshy,MPI_COMPLEX,left,tag(3),comm2D,lsend(3),ierror)
    end if
    if (right/=MPI_PROC_NULL) then
      call MPI_IRECV(ghost_l1_rec,nmeshz*nmeshy,MPI_COMPLEX,right,tag(1),comm2D,lrec(1),ierror)
      call MPI_IRECV(ghost_l2_rec,nmeshz*nmeshy,MPI_COMPLEX,right,tag(2),comm2D,lrec(2),ierror)
      call MPI_IRECV(ghost_l3_rec,nmeshz*nmeshy,MPI_COMPLEX,right,tag(3),comm2D,lrec(3),ierror)
    end if
    if (left/=MPI_PROC_NULL) then
      call MPI_WAIT(lsend(1),stat,ierror) ; call MPI_WAIT(lsend(2),stat,ierror) ; call MPI_WAIT(lsend(3),stat,ierror)
    end if
    if (right/=MPI_PROC_NULL) then
      call MPI_WAIT(lrec(1),stat,ierror) ; call MPI_WAIT(lrec(2),stat,ierror) ;  call MPI_WAIT(lrec(3),stat,ierror)
      !now implement these ghostzones
      Psi(1:nmeshz,1:nmeshy,nmeshx+1)=ghost_l1_rec(:,:)
      Psi(1:nmeshz,1:nmeshy,nmeshx+2)=ghost_l2_rec(:,:)
      Psi(1:nmeshz,1:nmeshy,nmeshx+3)=ghost_l3_rec(:,:)
    else
      select case(boundaries)
      case('zero')
      Psi(1:nmeshz,1:nmeshy,nmeshx+1:nmeshx+3)=(0.,0.) !zero boundaries
      case('solid')
      Psi(1:nmeshz,1:nmeshy,nmeshx+1)=Psi(1:nmeshz,1:nmeshy,nmeshx-1) !solid
      Psi(1:nmeshz,1:nmeshy,nmeshx+2)=Psi(1:nmeshz,1:nmeshy,nmeshx-2)
      Psi(1:nmeshz,1:nmeshy,nmeshx+3)=Psi(1:nmeshz,1:nmeshy,nmeshx-3)
      end select
    end if
    !----------------send to the right-----------------
    if (right/=MPI_PROC_NULL) then
      ghost_r1(1:nmeshz,1:nmeshy)=Psi(1:nmeshz,1:nmeshy,nmeshx)
      ghost_r2(1:nmeshz,1:nmeshy)=Psi(1:nmeshz,1:nmeshy,nmeshx-1)
      ghost_r3(1:nmeshz,1:nmeshy)=Psi(1:nmeshz,1:nmeshy,nmeshx-2)
      call MPI_ISEND(ghost_r1,nmeshz*nmeshy,MPI_COMPLEX,right,tag(1),comm2D,rsend(1),ierror)
      call MPI_ISEND(ghost_r2,nmeshz*nmeshy,MPI_COMPLEX,right,tag(2),comm2D,rsend(2),ierror)
      call MPI_ISEND(ghost_r3,nmeshz*nmeshy,MPI_COMPLEX,right,tag(3),comm2D,rsend(3),ierror)
    end if
    if (left/=MPI_PROC_NULL) then
      call MPI_IRECV(ghost_r1_rec,nmeshz*nmeshy,MPI_COMPLEX,left,tag(1),comm2D,rrec(1),ierror)
      call MPI_IRECV(ghost_r2_rec,nmeshz*nmeshy,MPI_COMPLEX,left,tag(2),comm2D,rrec(2),ierror)
      call MPI_IRECV(ghost_r3_rec,nmeshz*nmeshy,MPI_COMPLEX,left,tag(3),comm2D,rrec(3),ierror)
    end if
    if (right/=MPI_PROC_NULL) then
      call MPI_WAIT(rsend(1),stat,ierror) ; call MPI_WAIT(rsend(2),stat,ierror) ; call MPI_WAIT(rsend(3),stat,ierror)
    end if
    if (left/=MPI_PROC_NULL) then
      call MPI_WAIT(rrec(1),stat,ierror) ; call MPI_WAIT(rrec(2),stat,ierror) ;  call MPI_WAIT(rrec(3),stat,ierror)
      !now implement these ghostzones
      Psi(1:nmeshz,1:nmeshy,0)=ghost_r1_rec(:,:)
      Psi(1:nmeshz,1:nmeshy,-1)=ghost_r2_rec(:,:)
      Psi(1:nmeshz,1:nmeshy,-2)=ghost_r3_rec(:,:)
    else
      select case(boundaries)
      case('zero')
      Psi(1:nmeshz,1:nmeshy,-2:0)=(0.,0.) !zero boundaries
      case('solid')
      Psi(1:nmeshz,1:nmeshy,0)=Psi(1:nmeshz,1:nmeshy,2) !solid
      Psi(1:nmeshz,1:nmeshy,-1)=Psi(1:nmeshz,1:nmeshy,3)
      Psi(1:nmeshz,1:nmeshy,-2)=Psi(1:nmeshz,1:nmeshy,4)
      end select
    end if
    !----------------send upwards---------------------
    if (top/=MPI_PROC_NULL) then
      ghost_u1(1:nmeshz,1:nmeshx)=Psi(1:nmeshz,1,1:nmeshx)
      ghost_u2(1:nmeshz,1:nmeshx)=Psi(1:nmeshz,2,1:nmeshx)
      ghost_u3(1:nmeshz,1:nmeshx)=Psi(1:nmeshz,3,1:nmeshx)
      call MPI_ISEND(ghost_u1,nmeshz*nmeshx,MPI_COMPLEX,top,tag(1),comm2D,usend(1),ierror)
      call MPI_ISEND(ghost_u2,nmeshz*nmeshx,MPI_COMPLEX,top,tag(2),comm2D,usend(2),ierror)
      call MPI_ISEND(ghost_u3,nmeshz*nmeshx,MPI_COMPLEX,top,tag(3),comm2D,usend(3),ierror)
    end if
    if (bot/=MPI_PROC_NULL) then
      call MPI_IRECV(ghost_u1_rec,nmeshz*nmeshx,MPI_COMPLEX,bot,tag(1),comm2D,urec(1),ierror)
      call MPI_IRECV(ghost_u2_rec,nmeshz*nmeshx,MPI_COMPLEX,bot,tag(2),comm2D,urec(2),ierror)
      call MPI_IRECV(ghost_u3_rec,nmeshz*nmeshx,MPI_COMPLEX,bot,tag(3),comm2D,urec(3),ierror)
    end if
    if (top/=MPI_PROC_NULL) then
      call MPI_WAIT(usend(1),stat,ierror) ; call MPI_WAIT(usend(2),stat,ierror) ; call MPI_WAIT(usend(3),stat,ierror)
    end if
    if (bot/=MPI_PROC_NULL) then
      call MPI_WAIT(urec(1),stat,ierror) ; call MPI_WAIT(urec(2),stat,ierror) ;  call MPI_WAIT(urec(3),stat,ierror)
      !now implement these ghostzones
      Psi(1:nmeshz,nmeshy+1,1:nmeshx)=ghost_u1_rec(:,:)
      Psi(1:nmeshz,nmeshy+2,1:nmeshx)=ghost_u2_rec(:,:)
      Psi(1:nmeshz,nmeshy+3,1:nmeshx)=ghost_u3_rec(:,:)
    else
      select case(boundaries)
      case('zero')
      Psi(1:nmeshz,nmeshy+1:nmeshy+3,1:nmeshx)=(0.,0.) !zero boundaries
      case('solid')
      Psi(1:nmeshz,nmeshy+1,1:nmeshx)=Psi(1:nmeshz,nmeshy-1,1:nmeshx)
      Psi(1:nmeshz,nmeshy+2,1:nmeshx)=Psi(1:nmeshz,nmeshy-2,1:nmeshx)
      Psi(1:nmeshz,nmeshy+3,1:nmeshx)=Psi(1:nmeshz,nmeshy-3,1:nmeshx)
      end select
    end if
    !----------------send down---------------------
    if (bot/=MPI_PROC_NULL) then
      ghost_d1(1:nmeshz,1:nmeshx)=Psi(1:nmeshz,nmeshy,1:nmeshx)
      ghost_d2(1:nmeshz,1:nmeshx)=Psi(1:nmeshz,nmeshy-1,1:nmeshx)
      ghost_d3(1:nmeshz,1:nmeshx)=Psi(1:nmeshz,nmeshy-2,1:nmeshx)
      call MPI_ISEND(ghost_d1,nmeshz*nmeshx,MPI_COMPLEX,bot,tag(1),comm2D,dsend(1),ierror)
      call MPI_ISEND(ghost_d2,nmeshz*nmeshx,MPI_COMPLEX,bot,tag(2),comm2D,dsend(3),ierror)
      call MPI_ISEND(ghost_d3,nmeshz*nmeshx,MPI_COMPLEX,bot,tag(3),comm2D,dsend(2),ierror)
    end if
    if (top/=MPI_PROC_NULL) then
      call MPI_IRECV(ghost_d1_rec,nmeshz*nmeshx,MPI_COMPLEX,top,tag(1),comm2D,drec(1),ierror)
      call MPI_IRECV(ghost_d2_rec,nmeshz*nmeshx,MPI_COMPLEX,top,tag(2),comm2D,drec(2),ierror)
      call MPI_IRECV(ghost_d3_rec,nmeshz*nmeshx,MPI_COMPLEX,top,tag(3),comm2D,drec(3),ierror)
    end if
    if (bot/=MPI_PROC_NULL) then
      call MPI_WAIT(dsend(1),stat,ierror) ; call MPI_WAIT(dsend(2),stat,ierror) ; call MPI_WAIT(dsend(3),stat,ierror)
    end if
    if (top/=MPI_PROC_NULL) then
      call MPI_WAIT(drec(1),stat,ierror) ; call MPI_WAIT(drec(2),stat,ierror) ;  call MPI_WAIT(drec(3),stat,ierror)
      !now implement these ghostzones
      Psi(1:nmeshz,0,1:nmeshx)=ghost_d1_rec(:,:)
      Psi(1:nmeshz,-1,1:nmeshx)=ghost_d2_rec(:,:)
      Psi(1:nmeshz,-2,1:nmeshx)=ghost_d3_rec(:,:)
    else
      select case(boundaries)
      case('zero')
      Psi(1:nmeshz,-2:0,1:nmeshx)=(0.,0.) !zero boundaries
      case('solid')
      Psi(1:nmeshz,0,1:nmeshx)=Psi(1:nmeshz,2,1:nmeshx)
      Psi(1:nmeshz,-1,1:nmeshx)=Psi(1:nmeshz,3,1:nmeshx)
      Psi(1:nmeshz,-2,1:nmeshx)=Psi(1:nmeshz,4,1:nmeshx)
      end select
    end if
    call MPI_BARRIER(comm2D,ierror)
    !finally must account for boundaries in z
    if (periodicz) then
      Psi(nmeshz+1,1:nmeshy,1:nmeshx)=Psi(1,1:nmeshy,1:nmeshx)
      Psi(nmeshz+2,1:nmeshy,1:nmeshx)=Psi(2,1:nmeshy,1:nmeshx)
      Psi(nmeshz+3,1:nmeshy,1:nmeshx)=Psi(3,1:nmeshy,1:nmeshx)
      Psi(0,1:nmeshy,1:nmeshx)=Psi(nmeshz,1:nmeshy,1:nmeshx)
      Psi(-1,1:nmeshy,1:nmeshx)=Psi(nmeshz-1,1:nmeshy,1:nmeshx)
      Psi(-2,1:nmeshy,1:nmeshx)=Psi(nmeshz-2,1:nmeshy,1:nmeshx)
    else
    select case(boundaries)
      case('periodic')
        Psi(nmeshz+1,1:nmeshy,1:nmeshx)=Psi(1,1:nmeshy,1:nmeshx)
        Psi(nmeshz+2,1:nmeshy,1:nmeshx)=Psi(2,1:nmeshy,1:nmeshx)
        Psi(nmeshz+3,1:nmeshy,1:nmeshx)=Psi(3,1:nmeshy,1:nmeshx)
        Psi(0,1:nmeshy,1:nmeshx)=Psi(nmeshz,1:nmeshy,1:nmeshx)
        Psi(-1,1:nmeshy,1:nmeshx)=Psi(nmeshz-1,1:nmeshy,1:nmeshx)
        Psi(-2,1:nmeshy,1:nmeshx)=Psi(nmeshz-2,1:nmeshy,1:nmeshx)
       case('zero')
        Psi(nmeshz+1:nmeshz+3,1:nmeshy,1:nmeshx)=0.
        Psi(-2:0,1:nmeshy,1:nmeshx)=0.
       case('solid')
        !call emergency_stop('still to do!')
        Psi(0,1:nmeshy,1:nmeshx)=Psi(2,1:nmeshy,1:nmeshx)
        Psi(-1,1:nmeshy,1:nmeshx)=Psi(3,1:nmeshy,1:nmeshx)
        Psi(-2,1:nmeshy,1:nmeshx)=Psi(4,1:nmeshy,1:nmeshx)
        Psi(nmeshz+1,1:nmeshy,1:nmeshx)=Psi(nmeshz-1,1:nmeshy,1:nmeshx)
        Psi(nmeshz+2,1:nmeshy,1:nmeshx)=Psi(nmeshz-2,1:nmeshy,1:nmeshx)
        Psi(nmeshz+3,1:nmeshy,1:nmeshx)=Psi(nmeshz-3,1:nmeshy,1:nmeshx)
     end select
     end if
  end subroutine
end module
