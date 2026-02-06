! --------------------------------------------------------------------------------
! Module created to handle hdf5 file i/o
! Create, append and read hdf5 files needed for mcdust
! --------------------------------------------------------------------------------
! File structure
! --------------------------------------------------------------------------------
! Attributes
! code - mcdust v1.0
! author - Dr Joanna Drazkowska, Vignesh Vaikundaraman, Nerea Gurrutxaga
! mass_of_swarm - mswarm - mass of a swarm
! output_number - nout - output number of file
! number_of_particles_per_cell - ncell 
! number_or_radial_zones - nr
! number_of_vertical_zones - nz
! steps_between_outputs - fout
! time_between_outputs - dtime
! maximum_time_of_simulation - tend
! minimum_radius_[AU] - minrad
! maximum_radius_[AU] - maxrad
! monomer_radius_[cm] - a0
! material_density_[g/cm3] - matdens
! minimum_mass_difference - dmmax                             
! evaporation_radius_[AU] - smallr 
! dust_to_gas_ratio  - dtg
! restart  - restart
! fragmentation_velocity_[cm/s] - vfrag
! alpha   - alpha_t                       
! sigma_gas_[g/cm2] - sigmag0
! temperature_[K]  - t0
! erosion_mass_ratio - erosion_mass_ratio
! data_directory - datadir
! --------------------------------------------------------------------------------
! Data
! /swarms/swarmsout - data of the swarms (id, r, z, St, m, velr, velz)
! /times/timesout - snapshot time of the data
! /times/timestep - the latest advection timestep before output
! --------------------------------------------------------------------------------
module hdf5output
  use types,       only: swarm, list_of_swarms
  use constants
  use hdf5
  use parameters
  implicit none

  private
  public      :: hdf5_file_write, hdf5_file_t, hdf5_file_close, hdf5_file_read, hdf5_write_att

  interface hdf5_write_att
    module procedure hdf5_write_att_str
    module procedure hdf5_write_att_int
    module procedure hdf5_write_att_real
  end interface hdf5_write_att
   
  type hdf5_file_t
    integer(HID_T) :: id
    integer(HID_T) :: swarms
    integer(HID_T) :: times
    integer(HID_T) :: params
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
  subroutine hdf5_file_write (this, swrm, time, outnr,mswarm, resdt)

    implicit none
    type(hdf5_file_t), intent (out) :: this
    character(len=15) :: filename
    character(len=200) :: filepath
    type (swarm), dimension(:), intent(in) :: swrm
    integer, intent (in) :: outnr
    real, intent(in) :: mswarm
    real, intent (in) :: time, resdt
    real :: timeout, resdtout
    integer :: error
    write(filename,'(a7,i5.5,a3)') 'swarms-',outnr,'.h5'
    filepath = trim(path)//"/"//trim(filename)
    call h5open_f(error)
    resdtout = resdt/year
    call hdf5_file_create(this, filepath, outnr , mswarm, resdtout)     ! creates file
    call hdf5_create_dataspace(this, swrm)                  ! creates dataspace and adds swarm  data
    timeout = time/year
    call hdf5_time_create(this, timeout)                    ! creates dataspace for time outputs
    call hdf5_file_close(this,this%swarms)
    call h5close_f(error)    
    return
  end subroutine hdf5_file_write

  subroutine swarm_unit_conversion(swrm, tempswarm)
    implicit none
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
    implicit none
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
    write(*,*) 'create and write successfull'
    call h5pclose_f(this%prop_list, this%error)
    call h5dclose_f(this%dset_id, this% error)
    call h5sclose_f(this%dspace_id, this% error)
    call h5gclose_f(this%swarms, this%error)

    return
  end subroutine hdf5_create_dataspace

! subroutine to create dataset for writing output time
  subroutine hdf5_time_create(this, times)
    implicit none
    real, intent(in), target :: times
    type(hdf5_file_t), intent(inout) :: this
    integer, parameter :: rank_t = 1
    type(C_PTR) :: f_ptr
    integer(HSIZE_T) :: time_dims(rank_t), time_maxdims(rank_t), time_chunkdims(rank_t)

    f_ptr = C_LOC(times)
    time_dims(1) = 1
    time_chunkdims(1) = 10
    time_maxdims(1) = H5S_UNLIMITED_F
    ! write time of the snapshot
    call h5screate_simple_f(rank_t, time_dims, this%t_dspace_id, this%error, time_maxdims) ! create chunked dataspace
    call h5pcreate_f(H5P_DATASET_CREATE_F, this%t_prop_list, this%error) ! create chunkable dataset
    call h5pset_chunk_f(this%t_prop_list, rank_t, time_chunkdims, this%error) ! chunk the dataset

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
    implicit none
    
    type(hdf5_file_t), intent(inout)  :: this
    type(swarm), dimension(:), intent(in)        :: tempswrm
    type(swarm), dimension(size(tempswrm)), TARGET        :: swrm
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
    call h5tinsert_f(this%memtype, "Radial velocity v_r [cm/s]", H5OFFSETOF(C_LOC(swrm(1)),C_LOC(swrm(1)%velr)), &
                                                                                  H5T_NATIVE_DOUBLE, this%error)                                                                              
    call h5tinsert_f(this%memtype, "Vertical velocity v_z [cm/s]", H5OFFSETOF(C_LOC(swrm(1)),C_LOC(swrm(1)%velz)), &
                                                                                  H5T_NATIVE_DOUBLE, this%error)
    return
  end subroutine create_swarm_type
  
  ! subroutine to create a hdf5 file 
  subroutine hdf5_file_create(this, filename, nout, mswarm, resdtout)
    implicit none
    type(hdf5_file_t), intent(out) :: this
    character(len=*), intent(in) :: filename
    real :: mswarm, resdtout
    integer :: nout
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, this% id, this% error) ! creates file
    call hdf5_write_att(this%id, 'code', 'mcdust   v1.0')
    call hdf5_write_att(this%id, 'authors', 'Dr Joanna Drazkowska, Vignesh Vaikundaraman, Nerea Gurrutuxaga')
    call hdf5_write_att(this%id, 'output_number', nout)
    call hdf5_write_att(this%id, 'mass_of_swarm[g]', mswarm)
    call hdf5_write_att(this%id, 'number_of_particles', Ntot)
    call hdf5_write_att(this%id, 'number_of_particles_per_cell', ncell)
    call hdf5_write_att(this%id, 'number_of_radial_zones', nr)
    call hdf5_write_att(this%id, 'number_of_vertical_zones', nz)
    call hdf5_write_att(this%id, 'steps_between_outputs', fout)
    call hdf5_write_att(this%id, 'time_between_outputs[yr]', dtime)
    call hdf5_write_att(this%id, 'maximum_time_of_simulation[yr]', tend)
    call hdf5_write_att(this%id, 'minimum_radius_[AU]', minrad0)
    call hdf5_write_att(this%id, 'maximum_radius_[AU]', maxrad0)
    call hdf5_write_att(this%id, 'monomer_radius_[cm]', r0)
    call hdf5_write_att(this%id, 'material_density_[g/cm3]', matdens)
    call hdf5_write_att(this%id, 'dmmax', dmmax)
    call hdf5_write_att(this%id, 'evaporation_radius_[AU]', smallr)
    call hdf5_write_att(this%id, 'dust_to_gas_ratio', dtg)  
    call hdf5_write_att(this%id, 'fragmentation_velocity_[cm/s]', vfrag)
    call hdf5_write_att(this%id, 'alpha_t', alpha_t)
    call hdf5_write_att(this%id, 'sigma_gas_[g/cm2]', sigmag0)
    call hdf5_write_att(this%id, 'temperature_[K]', temperature)
    call hdf5_write_att(this%id, 'pressure_gradient_eta', eta)
    call hdf5_write_att(this%id, 'data_directory', datadir) 
    call hdf5_write_att(this%id, 'final_timestep[yr]', resdtout) 
#ifdef EROSION
    call hdf5_write_att(this%id, 'erosion_mass_ratio', erosion_mass_ratio)
#endif
    call h5gcreate_f(this% id, 'swarms', this% swarms, this% error)    ! creates group for swarm data to be written into
    call h5gcreate_f(this% id, 'times', this% times, this% error)      ! creates group for timesout
    return
  end subroutine hdf5_file_create

  !subroutine to write attributes
  subroutine hdf5_write_att_str(id, name, value)
    implicit none
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
  end subroutine hdf5_write_att_str

  subroutine hdf5_write_att_real(id, name, value)
    implicit none
    integer(hid_t) :: id
    integer(hid_t) :: attr_id
    integer(hid_t) :: aspace_id
    integer(size_t), dimension(1) :: adims = (/0/)
    character(len=*), intent(in)  :: name
    real :: value
    integer     :: ierr

    call h5screate_f(H5S_SCALAR_F, aspace_id, ierr)
    call h5acreate_f(id, name, H5T_NATIVE_REAL, &
                     aspace_id, attr_id, ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, value, adims, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5sclose_f(aspace_id, ierr)
  end subroutine hdf5_write_att_real

  subroutine hdf5_write_att_int(id, name, value)
    implicit none
    integer(hid_t) :: id
    integer(hid_t) :: attr_id
    integer(hid_t) :: aspace_id
    integer(size_t), dimension(1) :: adims = (/0/)
    character(len=*), intent(in)  :: name
    integer :: value
    integer :: ierr

    call h5screate_f(H5S_SCALAR_F, aspace_id, ierr)
    call h5acreate_f(id, name, H5T_NATIVE_REAL, &
                     aspace_id, attr_id, ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, value, adims, ierr)
    call h5aclose_f(attr_id, ierr)
    call h5sclose_f(aspace_id, ierr)
  end subroutine hdf5_write_att_int

  ! subroutine to read the last row 
  subroutine hdf5_file_read(Ntot, swrm, nout, mswarm, time, resdt)
    implicit none
    type(hdf5_file_t) :: this
    type(swarm), dimension(:), allocatable, TARGET ::swrm
    type(C_PTR) :: f_ptr
    integer, intent(in) :: Ntot
    integer, parameter :: rank = 2
    integer(HSIZE_T) :: dims(rank), maxdims(rank), offset(rank), counts(rank)
    integer(HID_T) :: mem_space
    integer :: error, i, nout
    real :: mswarm, time, resdt
    if (.not.allocated(swrm)) allocate( swrm(Ntot) )

    ! open file and group 
    call h5open_f(error)
    call h5fopen_f('restart.h5', H5F_ACC_RDONLY_F, this% id, this% error)
    call check_error(this% error, 'error opening the file restart.h5')
    call h5gopen_f(this% id, 'swarms', this%swarms, this% error)
    call check_error(this% error, 'swarms group opening error')

    ! open dataset, dataspace and get data dimensions
    call h5dopen_f (this%swarms, 'swarmsout', this%dset_id, this%error)
    CALL h5dget_space_f(this%dset_id, this%dspace_id, error)
    call h5sget_simple_extent_dims_f(this%dspace_id, dims, maxdims, this%error)
    
    ! create swarm datatype
    call create_swarm_type(this,swrm)
    ! initiate offset and counts to read the last row of the file
    offset(1:rank-1) = 0
    offset(rank) = dims(rank) - 1
    counts(1:rank-1) = shape(swrm)
    counts(rank) = 1
    
    ! create dataspace for the one row
    call h5screate_simple_f(rank-1, counts(1), mem_space, this% error)
    ! define hyperslab to be read in the file
    call h5sselect_hyperslab_f(this%dspace_id, H5S_SELECT_SET_F, offset, counts, this%error)
    ! read using defined hyperslab
    f_ptr = C_LOC(swrm(1))
    call h5dread_f(this%dset_id, this%memtype, f_ptr, this%error, mem_space, this%dspace_id)
    ! close dependencies
    call h5dclose_f(this%dset_id, this%error)
    call h5sclose_f(this%dspace_id, this%error)
    call hdf5_file_close(this, this%swarms)
    
    do i = 1, Ntot
         swrm(i)%rdis = swrm(i)%rdis * AU
         swrm(i)%zdis = swrm(i)%zdis * AU
    enddo
    resdt = resdt*year
    ! now reading snapshotime
    call hdf5_read_time(time)
    time = time * year
    call hdf5_read_att(nout,mswarm,resdt)
    call h5close_f(error)

    return
  end subroutine hdf5_file_read
  ! subroutine to read snapshottime
  subroutine hdf5_read_time(time)
    implicit none
    type(hdf5_file_t) :: this 
    type(C_PTR) :: f_ptr
    real, TARGET :: time
    integer(hsize_t), dimension(1) :: dims, maxdims

    ! reading in snapshot time
    call h5fopen_f('restart.h5', H5F_ACC_RDONLY_F, this% id, this% error)
    call check_error(this% error, 'error opening the file restart.h5')
    call h5gopen_f(this% id, 'times', this%times, this% error)
    call check_error(this% error, 'times group opening error')

    ! open dataset, dataspace and get data dimensions
    call h5dopen_f (this%times, 'timesout', this%dset_id, this%error)
    CALL h5dget_space_f(this%dset_id, this%dspace_id, this%error)
    call h5sget_simple_extent_dims_f(this%dspace_id, dims, maxdims, this%error)

    f_ptr = C_LOC(time)
    call h5dread_f(this%dset_id, H5T_NATIVE_DOUBLE, f_ptr, this%error, this%dspace_id)
    ! close dependencies
    call h5dclose_f(this%dset_id, this%error)
    call h5sclose_f(this%dspace_id, this%error)
    call hdf5_file_close(this, this%times)

  end subroutine hdf5_read_time

  subroutine hdf5_read_att(nout, mswarm, resdt)
    implicit none
    type(hdf5_file_t) :: this
    type(C_PTR) :: f_ptr
    real, target :: outnr, mswarm, resdt
    integer(HSIZE_T), dimension(1) :: dims
    integer(HID_T) :: att1, att2, att3
    integer :: nout

    dims=1

    call h5fopen_f('restart.h5', H5F_ACC_RDONLY_F, this%id, this%error)
    call check_error(this% error, 'error opening the file restart.h5')

    CALL h5aopen_name_f(this%id, 'mass_of_swarm[g]', att1, this%error)
    call check_error(this% error, 'error opening the attribute m_swarm')

    f_ptr = C_LOC(mswarm)
    CALL h5aread_f(att1, H5T_NATIVE_DOUBLE, mswarm, dims, this%error)
    CALL h5aclose_f(att1, this%error)                  

    CALL h5aopen_name_f(this%id, 'output_number', att2, this%error)
    call check_error(this% error, 'error opening the attribute output_number')

    f_ptr = C_LOC(outnr)
    CALL h5aread_f(att2, H5T_NATIVE_INTEGER, nout, dims, this%error)
    CALL h5aclose_f(att2, this%error)

    CALL h5aopen_name_f(this%id, 'final_timestep[yr]', att3, this%error)
    call check_error(this% error, 'error opening the attribute final_timestep[yr]')

    f_ptr = C_LOC(resdt)
    CALL h5aread_f(att3, H5T_NATIVE_INTEGER, resdt, dims, this%error)
    CALL h5aclose_f(att3, this%error)

    
    call h5fclose_f(this%id, this%error)
    call check_error(this% error, 'error closing the file restart.h5')
  end subroutine hdf5_read_att

  ! subroutine to close hdf5 file
  subroutine hdf5_file_close(this,group)
    implicit none
    type(hdf5_file_t), intent(inout) :: this
    integer(HID_T) :: group
    logical :: valid

    call h5iis_valid_f(group, valid, this% error)
    if (valid) then
        call h5gclose_f(group, this% error)
    end if
    call h5fclose_f(this% id, this% error)
    call check_error(this% error, 'file close error')

    return
  end subroutine hdf5_file_close


  ! routine to check error
  subroutine check_error(e, msg)
    implicit none
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

   
