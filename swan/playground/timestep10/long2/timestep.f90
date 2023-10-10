! this module estimates maximum time step of the code taking into account Courant-like condition for advection
! if you implement other processes requiring time step limitation, put corresponding subroutine here
! and take the result into account in the time_step routine
module timestep

   use constants,  only: pi, third, mH2, AH2, year, Ggrav, Msun
   use discstruct, only: omegaK
   use grid,       only: g
   use parameters, only: smallr, vertsett, rdrift
   use types,      only: swarm, list_of_swarms
   use advection,  only: stokesnr, vel_rd, vel_vn

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
      integer                                        :: i, j, c=10

      old_dt = final_dtime
      final_dtime = dtout
      
      
      !final_dtime = min(final_dtime, size(bin(1,1)%p) * old_dt / real(maxval(ncolls(:,:)))) !we don't want more collisions than ~ nr of particles in a zone per timestep
      !if (vertsett) final_dtime = min(final_dtime, g%dz(1,1)/maxval(abs(bin(:,:)%p(:)%velz))) ! don't let particles jump over cells in z
      !if (rdrift) final_dtime = min(final_dtime, g%dr(1)/maxval(abs(bin(:,:)%p(:)%velr))) ! don't let particles jump over cells in r
      

      do i = 1, size(ncolls(:,:),dim=1)
        do j = 1, size(ncolls(:,:),dim=2)
          final_dtime = min(final_dtime, c*size(bin(i,j)%p) * old_dt / real(max(1,ncolls(i,j)))) !we don't want more collisions than ~ nr of particles in a zone per timestep
          if (vertsett) final_dtime = min(final_dtime, c*g%dz(i,j)/maxval(abs(bin(i,j)%p(:)%velz))) ! don't let particles jump over cells in z
          if (rdrift) final_dtime = min(final_dtime, c*g%dr(i)/maxval(abs(bin(i,j)%p(:)%velr))) ! don't let particles jump over cells in r
        enddo
      enddo

      
      return
   end subroutine time_step

end
