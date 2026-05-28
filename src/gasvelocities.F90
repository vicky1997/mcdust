!subroutines to read velocity data
module gasvelocities
    use constants
    use parameters
    use utils
    implicit none
    private
#ifdef READGASVELOCITES_R
    public :: read_data_r, r_data_vgr, zh_data_vgr, vgr_data 
#endif
#ifdef READGASVELOCITES_Z
    public :: read_data_z, r_data_vgz, zh_data_vgz, vgz_data
#endif

    real, dimension(:), allocatable :: r_data, zh_data
    real, dimension(:,:),allocatable :: vgr_data, vgz_data
    
    contains
#ifdef READGASVELOCITES_R
! to read radial gas velocity data
   subroutine read_data_r()
      implicit none
      real, dimension(:), allocatable :: v_cs_data, rr_data !values from the data file which is scale free model
      integer :: i, nr !loop variable
      real :: nrr
      open(1,file=velfile_r,action='read')
      read(1,*) nrr
      nr = int(nrr)
      allocate(r_data_vgr(nr),v_cs_data(nr*nr),zh_data_vgr(nr),rr_data(nr*nr))
      read(1,*)
      do i=1,nr
         read(1,*) r_data_vgr(i), zh_data_vgr(i)
      enddo
      read(1,*)
      do i=1,nr*nr
         read(1,*) rr_data(i), v_cs_data(i)
      enddo
      close(1)
      
      r_data_vgr(:) = r_data_vgr(:)*AU 
   
      allocate(vgr_data(nr,nr))
      vgr_data = reshape(v_cs_data,(/nr,nr/))
      deallocate(vr_cs_data)
      return
   end subroutine read_data_r
#endif
   
#ifdef READGASVELOCITES_Z
! to read vertical gas velocity data
subroutine read_data_z()
   implicit none
   real, dimension(:), allocatable :: v_cs_data, rr_data !values from the data file which is scale free model
   integer :: i, nr !loop variable
   real :: nrr
   open(1,file=velfile_r,action='read')
   read(1,*) nrr
   nr = int(nrr)
   allocate(r_data_vgz(nr),v_cs_data(nr*nr),zh_data_vgz(nr),rr_data(nr*nr))
   read(1,*)
   do i=1,nr
      read(1,*) r_data_vgz(i), zh_data_vgz(i)
   enddo
   read(1,*)
   do i=1,nr*nr
      read(1,*) rr_data(i), v_cs_data(i)
   enddo
   close(1)
   
   r_data_vgz(:) = r_data_vgz(:)*AU 

   allocate(vgz_data(nr,nr))
   vgz_data = reshape(v_cs_data,(/nr,nr/))
   deallocate(vr_cs_data)
   return
end subroutine read_data_z
#endif
   
end module gasvelocities
