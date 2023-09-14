! 
! module created to handle hdf5 file i/o
! it can create, append and read hdf5 files needed for 2DMC
! inspiration for h5md fortran package for molecular data.
!
!
module hdf5output
  use types,       only: swarm, list_of_swarms
  use parameters,  only: con2, matdens, vfrag, nbins, smallr, path
  use constants
  use discstruct,  only: densg, sigmag, cs, omegaK, Pg, alpha
  use grid,        only: g
  use initproblem, only: m0, mswarm, nord, nmonom0
  use hdf5
  implicit none

  private
  public      :: hdf5_file_write, hdf5_file_t, hdf5_file_close, hdf5_file_read
   
  type hdf5_file_t
    integer(HID_T) :: id
    integer(HID_T) :: swarms
    integer(HID_T) :: times
    integer(HID_T) :: memtype
    integer(HID_T) :: dset_id
    integer(HID_T) :: dspace_id
    integer(HID_T) :: prop_list
    integer(HID_T) :: t_dset_id
    integer(HID_T) :: t_dspace_id
    integer(HID_T) :: t_prop_list
    integer :: error

  end type hdf5_file_t
  contains
  subroutine hdf5_file_write (this, swrm, time, mode, outnr)
    type(hdf5_file_t), intent (out) :: this
    character(len=15) :: filename
    character(len=200) :: filepath
    !character(len=41) ::path
    type (swarm), dimension(:), intent(in) :: swrm
    character(len=*), intent (in) :: mode
    integer, intent (in) :: outnr
    real, intent (in) :: time
    real :: timeout
    integer :: error
    !filename = 'swarms.h5'
    write(filename,'(a7,i5.5,a3)') 'swarms-',outnr,'.h5'
    !path = '/scratch/vaikundaraman/2DMC/hdf5_v2/data/'
    filepath = trim(path)//"/"//trim(filename)
    write(*,*) filepath
    !call h5open_f(error)
    if (mode .eq. 'create') then
      call h5open_f(error)
      call hdf5_file_create(this, filepath)                   ! creates file
      call hdf5_create_dataspace(this, swrm)                  ! creates dataspace and adds swarm  data
      timeout = time/year
      call hdf5_time_create(this, timeout)                    ! creates dataspace for time outputs
      call hdf5_file_close(this)
      call h5close_f(error)    
    else
      write(*,*) 'ERROR : only create mode works for this version'
      stop
    endif
    !call hdf5_file_close(this)
    !call h5close_f(error)
    return
  end subroutine hdf5_file_write

  subroutine swarm_unit_conversion(swrm, tempswarm)
    type(swarm), dimension(:), intent(in)  :: swrm
    type(swarm), dimension(:), intent(out) :: tempswarm
    integer :: i
    do i=1,size(swrm)
      tempswarm(i)%idnr = swrm(i)%idnr
      tempswarm(i)%mass = swrm(i)%mass
      tempswarm(i)%rdis = swrm(i)%rdis/AU
      tempswarm(i)%zdis = swrm(i)%zdis/AU
      tempswarm(i)%stnr = swrm(i)%stnr
    end do

    return
  end subroutine swarm_unit_conversion

  subroutine hdf5_create_dataspace (this, swrm)
    type(hdf5_file_t), intent (inout) :: this
    type (swarm), dimension(:), intent(in) :: swrm
    type (swarm), dimension(size(swrm)), TARGET :: tempswrm
    type(C_PTR) :: f_ptr
    integer, parameter :: rank = 2
     
    integer(HSIZE_T) :: dims(rank), maxdims(rank), chunkdims(rank)
    dims(1:rank-1) = shape(swrm)
    dims(rank) = 1
    chunkdims(1:rank-1) = shape(swrm)
    chunkdims(rank) = 1
    call swarm_unit_conversion(swrm, tempswrm)                   ! convert data to the units we want
    maxdims(1:rank-1) = shape(swrm)
    maxdims(rank) = H5S_UNLIMITED_F
    call create_swarm_type(this, swrm) ! create compund datatype of swarm
    call h5screate_simple_f(rank, dims, this%dspace_id, this% error, maxdims) ! create chunked dataspace
    call h5pcreate_f(H5P_DATASET_CREATE_F, this%prop_list, this%error) ! create chunkable dataset
    call h5pset_chunk_f(this%prop_list, rank, chunkdims, this%error) ! chunk the dataset
    
   !write data to the dataset
    call h5dcreate_f(this%swarms, 'swarmsout', this%memtype, this%dspace_id, this%dset_id, this% error, this%prop_list)
    f_ptr = C_LOC(tempswrm(1))
    call h5dwrite_f(this%dset_id, this%memtype, f_ptr, this% error)
    call h5pclose_f(this%prop_list, this%error)
    call h5dclose_f(this%dset_id, this% error)
    call h5sclose_f(this%dspace_id, this% error)
    call h5gclose_f(this%swarms, this%error)
    !write(*,*) this%dset_id, this%dspace_id

    return
  end subroutine hdf5_create_dataspace

! subroutine to create dataset for writing output time
  subroutine hdf5_time_create(this, times)
    real, intent(in), target :: times
    type(hdf5_file_t), intent(inout) :: this
    integer, parameter :: rank_t = 0
    type(C_PTR) :: f_ptr 
   !integer(HSIZE_T) :: time_dims(rank_t), time_maxdims(rank_t), time_chunkdims(rank_t)
   integer(HSIZE_T), dimension(1) :: time_dims
   time_dims(1) = 0

    f_ptr = C_LOC(times)
    call h5screate_simple_f(rank_t, time_dims, this%t_dspace_id, this%error)
    !call h5screate_simple_f(rank_t, this%t_dspace_id, this%error)
    ! write data to the dataset
    call h5dcreate_f(this%times, 'timesout', H5T_NATIVE_DOUBLE, this%t_dspace_id, this%t_dset_id, this% error, this%t_prop_list)
    call h5dwrite_f(this%t_dset_id, H5T_NATIVE_DOUBLE, f_ptr, this% error)

    ! close dependencies
    call h5pclose_f(this%t_prop_list, this%error)
    call h5dclose_f(this%t_dset_id, this% error)
    call h5sclose_f(this%t_dspace_id, this% error)
    call h5gclose_f(this%times, this%error)
    return
  end subroutine hdf5_time_create

! subroutine to add the swarm object as a compound datatype to the HDF5 file
  subroutine create_swarm_type(this,tempswrm)
    
    type(hdf5_file_t), intent(inout)  :: this
    type(swarm), dimension(:), intent(in)        :: tempswrm
    type(swarm), dimension(size(tempswrm)), TARGET        :: swrm
    !integer(HID_T) :: memtype
    swrm = tempswrm
    call h5tcreate_f(H5T_COMPOUND_F, H5OFFSETOF(C_LOC(swrm(1)),C_LOC(swrm(2))), & 
                                                                  this%memtype, this%error)
    call h5tinsert_f(this%memtype, "id number", H5OFFSETOF(C_LOC(swrm(1)),C_LOC(swrm(1)%idnr)), &
                                                                                  H5T_NATIVE_INTEGER, this%error)
    call h5tinsert_f(this%memtype, "mass of a particle [g]", H5OFFSETOF(C_LOC(swrm(1)),C_LOC(swrm(1)%mass)), &
                                                                                  H5T_NATIVE_DOUBLE, this%error)
    call h5tinsert_f(this%memtype, "cylindrical radius [AU]", H5OFFSETOF(C_LOC(swrm(1)),C_LOC(swrm(1)%rdis)), &
                                                                                  H5T_NATIVE_DOUBLE, this%error)
    call h5tinsert_f(this%memtype, "height above midplane [AU]", H5OFFSETOF(C_LOC(swrm(1)),C_LOC(swrm(1)%zdis)), &
                                                                                  H5T_NATIVE_DOUBLE, this%error)
    call h5tinsert_f(this%memtype, "Stokes number", H5OFFSETOF(C_LOC(swrm(1)),C_LOC(swrm(1)%stnr)), &
                                                                                  H5T_NATIVE_DOUBLE, this%error)
    return
  end subroutine create_swarm_type
  
  ! subroutine to create a hdf5 file 
  subroutine hdf5_file_create(this, filename)
    type(hdf5_file_t), intent(out) :: this
    character(len=*), intent(in) :: filename
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, this% id, this% error) ! creates file
    call hdf5_write_attribute(this%id, '2DMC', 'swarm properties at each output time')
    call hdf5_write_attribute(this%id, 'version', '0.1')
    call hdf5_write_attribute(this%id, 'author', 'Dr Joanna Drazkowska') 
    
    call h5gcreate_f(this% id, 'swarms', this% swarms, this% error)    ! creates group for swarm data to be written into
    call h5gcreate_f(this% id, 'times', this% times, this% error)      ! creates group for timesout
  
    return
  end subroutine hdf5_file_create

  !subroutine to write attributes
  subroutine hdf5_write_attribute(id, name, value)
    integer(HID_T), intent(inout) :: id
    character(len=*), intent(in)  :: name
    character(len=*), intent(in)  :: value

    integer(HID_T) :: a, s, t
    integer :: error
    integer(HSIZE_T) :: dims(1)
    integer(SIZE_T) :: length

    call h5screate_f(H5S_SCALAR_F, s, error)
    dims(1) = len(value)
    length = len(value)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, t, error)
    call h5tset_size_f(t, length, error)
    call h5tset_strpad_f(t, H5T_STR_NULLTERM_F, error)
    call h5acreate_f(id, name, t, s, a, error)
    call h5awrite_f(a, t, value, dims, error)
    call h5aclose_f(a, error)
    call h5tclose_f(t, error)
    call h5sclose_f(s, error)
    
    return
  end subroutine hdf5_write_attribute

  ! subroutine to read the last row 
  subroutine hdf5_file_read(Ntot, swrm)
    type(hdf5_file_t) :: this
    type(swarm), dimension(:), allocatable, TARGET ::swrm
    type(C_PTR) :: f_ptr
    integer, intent(in) :: Ntot
    integer, parameter :: rank = 1
    integer(HSIZE_T) :: dims(rank), maxdims(rank), offset(rank), counts(rank)
    integer(HID_T) :: mem_space
    integer :: error, i
    if (.not.allocated(swrm)) allocate( swrm(Ntot) )

    ! open file and group 
    call h5open_f(error)
    call h5fopen_f('restart.h5', H5F_ACC_RDONLY_F, this% id, this% error)
    call check_error(this% error, 'error opening the file restart.h5')
    call h5gopen_f(this% id, 'swarms', this%swarms, this% error)
    call check_error(this% error, 'group opening error')

    ! open dataset, dataspace and get data dimensions
    call h5dopen_f (this%swarms, 'swarmsout', this%dset_id, this%error)
    CALL h5dget_space_f(this%dset_id, this%dspace_id, error)
    call h5sget_simple_extent_dims_f(this%dspace_id, dims, maxdims, this%error)
    
    ! create swarm datatype
    call create_swarm_type(this,swrm)
    ! read 
    f_ptr = C_LOC(swrm(1))
    call h5dread_f(this%dset_id, this%memtype, f_ptr, this%error)
    ! close dependencies
    call h5dclose_f(this%dset_id, this%error)
    call h5sclose_f(this%dspace_id, this%error)
    !call h5fclose_f(this%id, this%error)
    call hdf5_file_close(this)
    call h5close_f(error)

    do i = 1, Ntot
         swrm(i)%rdis = swrm(i)%rdis * AU
         swrm(i)%zdis = swrm(i)%zdis * AU
    enddo

    return
  end subroutine hdf5_file_read
  
  ! subroutine to close hdf5 file
  subroutine hdf5_file_close(this)
    type(hdf5_file_t), intent(inout) :: this

    logical :: valid

    call h5iis_valid_f(this% swarms, valid, this% error)
    if (valid) then
        call h5gclose_f(this% swarms, this% error)
    end if
    call h5fclose_f(this% id, this% error)
    call check_error(this% error, 'file close error')

    return
  end subroutine hdf5_file_close


  ! routine to check error
  subroutine check_error(e, msg)
    integer, intent(in) :: e
    character(len=*), intent(in) :: msg

    if (e /= 0) then
       write(*,*) 'error code ', e
       write(*,*) msg
       stop
    end if

    return
  end subroutine check_error

end module hdf5output

   
