module mpi_var
  use cdata
  implicit none
  logical :: periodic(2)=.true. !are the boundaries periodic or not?
  integer :: rank, nprocs, ierror
  integer :: nprocx, nprocy !number of processes in x/y dir
  integer, dimension(MPI_STATUS_SIZE) :: stat
  integer :: comm2D !new MPI comm
  integer :: tag(3)=(/100,101,102/) !for message sending
  integer :: left, right, top, bot !cartesian array of processes
  integer :: dims(2)=0, coords(2) ! cartesian array of processes
  contains
  !----------------------------------------------------------------
  subroutine mpi_setup
    implicit none
      select case(boundaries)
        case('periodic')
          periodic=.true.
        case('zero','solid')
          periodic=.false.
        case default
          call emergency_stop('wrong boundaries set')
      end select
      call MPI_INIT(ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
      !divide up the processors to x and y co-ords
      call MPI_DIMS_CREATE(nprocs,2,dims,ierror)
      !YOU NEED TO CHECK OUT COMM2D VS MPI_COMM_WORLD
      call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periodic,.true.,comm2D,ierror)
      call MPI_CART_COORDS(comm2D, rank, 2, coords, ierror)
      call MPI_CART_SHIFT(comm2D,0,1,left,right,ierror)
      call MPI_CART_SHIFT(comm2D,1,1,top,bot,ierror)
      !now get the number of processes in the x and y direction
      nprocx=dims(1) ; nprocy=dims(2) !note it is an xy-plane
      !create local meshes
      nmeshx=Nx/nprocx ; nmeshy=Ny/nprocy ; nmeshz=Nz
      if (nmeshx*nprocx/=Nx.or.nmeshy*nprocy/=Ny) then
         call emergency_stop("unable to evenly subdivide grid check meshsize and number of processors")
         stop
      end if
      call init_random_seed
  end subroutine
  !----------------------------------------------------------------
  subroutine emergency_stop(err_string)
   implicit none
   character(*), intent(in) :: err_string

   if (rank == 0) then
     ! Print error to screen, write error to ERROR.
     print*, '*****FATAL ERROR ABORTING CODE*****'
     print*, err_string
     open (96, file = 'ERROR')
       write (96, *) err_string
     close (96)
   end if
   ! Stop the MPI process grid.
   call MPI_BARRIER(MPI_COMM_WORLD, ierror)
   call MPI_FINALIZE(ierror)
   ! Emergency stop.
   stop
  end subroutine emergency_stop
  !----------------------------------------------------------------
  SUBROUTINE init_random_seed()
    real :: test
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
  
    clock=seedgen(rank)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)
    DEALLOCATE(seed)
 !   call random_number(test) ; print*, rank, test
  END SUBROUTINE
  !----------------------------------------------------------------
  function seedgen(pid)
    !use iso_fortran_env
    implicit none

    integer :: seedgen
    integer, intent(IN) :: pid
    integer :: s

    call system_clock(s)
    seedgen = abs( mod((s*181)*((pid-83)*359), 104729) )
  end function seedgen
end module
