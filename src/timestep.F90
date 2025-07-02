! this module estimates maximum time step of the code taking into account Courant-like condition for advection
! if you implement other processes requiring time step limitation, put corresponding subroutine here
! and take the result into account in the time_step routine
module timestep

   use grid,       only: g
   use types,      only: swarm, list_of_swarms
   use constants,  only: AU
   private
   public      :: time_step

   contains

   ! gives back timestep as minimal timestep of all processes
   subroutine time_step(bin, ncolls, dtout, final_dtime)
      implicit none
      type(list_of_swarms), dimension(:,:), allocatable, target :: bin
      integer, dimension(:,:), allocatable           :: ncolls
      real, intent(in)                               :: dtout
      real, intent(inout)                            :: final_dtime   ! final timestep
      real                                           :: old_dt
      integer                                        :: i, j, c1=10, c2=1 !c1,c2 - timestep enlargement factors

      old_dt = final_dtime
      final_dtime = max(1.,dtout)

      do i = 1, size(ncolls(:,:),dim=1)
        do j = 1, size(ncolls(:,:),dim=2)
          final_dtime = min(final_dtime, c1*size(bin(i,j)%p) * old_dt / real(max(1,ncolls(i,j)))) !we don't want more collisions than ~ nr of particles in a zone per timestep
#ifdef VERSETTLING
          final_dtime = min(final_dtime, min(AU,c2*g%dz(i,j)/maxval(abs(bin(i,j)%p(:)%velz)))) ! don't let particles jump over cells in z
#endif
#ifdef RADRIFT
          final_dtime = min(final_dtime, c2*g%dr(i)/maxval(abs(bin(i,j)%p(:)%velr))) ! don't let particles jump over cells in r
#endif
         enddo
      enddo

      
      return
   end subroutine time_step

end
