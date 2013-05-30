module remesh
  use mpi_var
  contains
  subroutine remesh_Psi
    complex, allocatable :: Psi_rem(:,:,:), Psi_tmp(:,:,:)
    integer :: i, j, k
    integer :: coords_in(2)
    integer :: coords_send, coords_recv
    integer :: chunk_send, chunk_recv
    allocate(Psi_rem(1:nmeshz/2,1:nmeshy/2,1:nmeshx/2))
    do i=1, nmeshx/2
      do j=1, nmeshy/2
        do k=1, nmeshz/2
          !downsample taking the mean of 8 points
          !Psi_rem(k,j,i)=(Psi(2*k,2*j,2*i-1)+Psi(2*k,2*j,2*i)+&
          !               Psi(2*k,2*j-1,2*i-1)+Psi(2*k,2*j-1,2*i)+&
          !               Psi(2*k-1,2*j,2*i-1)+Psi(2*k-1,2*j,2*i)+&
          !               Psi(2*k-1,2*j-1,2*i-1)+Psi(2*k-1,2*j-1,2*i))/8.
          Psi_rem(k,j,i)=Psi(2*k-1,2*j-1,2*i-1)
        end do
      end do
    end do
    if (rank==0) then
      !allocate main array on root node
      allocate(Psi_tmp(Nz,Ny,Nx))
      Psi_tmp=(0.,0.);
      !put the root chunk into the array
      Psi_tmp(Nz/4+1:(Nz/4+Nz/2),(Ny/4+coords(2)*nmeshy/2+1):(Ny/4+(coords(2)+1)*nmeshy/2), &
            (Nx/4+coords(1)*nmeshx/2+1):(Nx/4+(coords(1)+1)*nmeshx/2))=Psi_rem(:,:,:)
    end if
    do i=1, nprocs-1
      if (rank==i) then
        !now get the coords of the other processes
        call MPI_ISEND(coords,2,MPI_INTEGER,0,tag(1),comm2D,coords_send,ierror)
        call MPI_WAIT(coords_send,stat,ierror)
        !send the Psi_rem chunk
        call MPI_ISEND(Psi_rem,nmeshz*nmeshy*nmeshx/8,MPI_COMPLEX,0,tag(2),comm2D,chunk_send,ierror)
        call MPI_WAIT(chunk_send,stat,ierror)
      end if
      if (rank==0) then
        call MPI_IRECV(coords_in,2,MPI_INTEGER,i,tag(1),comm2D,coords_recv,ierror)
        call MPI_WAIT(coords_recv,stat,ierror)
        call MPI_IRECV(Psi_rem,nmeshz*nmeshy*nmeshx/8,MPI_COMPLEX,i,tag(2),comm2D,chunk_recv,ierror)
        call MPI_WAIT(chunk_recv,stat,ierror)
        Psi_tmp(Nz/4+1:(Nz/4+Nz/2),(Ny/4+coords_in(2)*nmeshy/2+1):(Ny/4+(coords_in(2)+1)*nmeshy/2), &
               (Nx/4+coords_in(1)*nmeshx/2+1):(Nx/4+(coords_in(1)+1)*nmeshx/2))=Psi_rem(:,:,:)
      end if
      call MPI_BARRIER(comm2D,ierror)
    end do
    if (rank==0) then
      open(unit=98,file='./data/remesh.dat',status='replace',form='unformatted')
        write(98) Psi_tmp
      close(98)
    end if
    !now distribute back out
    do i=1, nprocs-1
      if (rank==i) then
        !now get the coords of the other processes
        call MPI_ISEND(coords,2,MPI_INTEGER,0,tag(1),comm2D,coords_send,ierror)
        call MPI_WAIT(coords_send,stat,ierror)
      end if
      if (rank==0) then
        call MPI_IRECV(coords_in,2,MPI_INTEGER,i,tag(1),comm2D,coords_recv,ierror)
        call MPI_WAIT(coords_recv,stat,ierror)
        rhs(1:nmeshz,1:nmeshy,1:nmeshx)=Psi_tmp(1:nmeshz,(coords_in(2)*nmeshy+1):((coords_in(2)+1)*nmeshy), &
                                               (coords_in(1)*nmeshx+1):((coords_in(1)+1)*nmeshx))
        call MPI_ISEND(rhs,nmeshz*nmeshy*nmeshx,MPI_COMPLEX,i,tag(3),comm2D,chunk_send,ierror)
        call MPI_WAIT(chunk_send,stat,ierror)
      end if
      if (rank==i) then
        call MPI_IRECV(rhs,nmeshz*nmeshy*nmeshx,MPI_COMPLEX,0,tag(3),comm2D,chunk_recv,ierror)
        call MPI_WAIT(chunk_recv,stat,ierror)
        Psi(1:nmeshz,1:nmeshy,1:nmeshx)=rhs(:,:,:)
      end if
      call MPI_BARRIER(comm2D,ierror)
    end do
    !and finally sort out roots section
    if (rank==0) then
      Psi(1:nmeshz,1:nmeshy,1:nmeshx)=Psi_tmp(1:nmeshz,(coords(2)*nmeshy+1):((coords(2)+1)*nmeshy), &
                                             (coords(1)*nmeshx+1):((coords(1)+1)*nmeshx))
      deallocate(Psi_tmp)
    end if
    deallocate(Psi_rem)
    hx=hx*2 ; hy=hy*2 ; hz=hz*2
    Lx=Lx*2 ; Ly=Ly*2 ; Lz=Lz*2
    xx=xx*2 ; yy=yy*2 ; zz=zz*2
    if (rank==0) then
      write(*,*) 'successfully remeshed code, box size now: ', Lx, 'x', Ly, 'x', Lz
    end if
  end subroutine
end module
