! this module estimates maximum time step of the code taking into account Courant-like condition for advection
! if you implement other processes requiring time step limitation, put corresponding subroutine here
! and take the result into account in the time_step routine
module timestep

   use constants,  only: pi, third, mH2, AH2, year
   use discstruct, only: omegaK
   use grid,       only: g
   use parameters, only: smallr, vertsett, eta
   use types,      only: swarm
   use advection,  only: stokesnr, vel_rd, vel_vn

   private
   public      :: time_step

   contains

   ! gives back timestep as minimal timestep of all processes
   subroutine time_step(swrm, final_dtime)
      implicit none
      type(swarm), dimension(:), allocatable, target :: swrm
      real, intent(inout)                            :: final_dtime   ! final timestep
      real                                           :: ver_dtime     ! maximum timestep for vertical settling
      real, parameter                                :: increase = 1.5

      call vertical(swrm, ver_dtime)
      final_dtime = max(final_dtime * increase,  1. / (eta * omegaK(swrm(1)%rdis))) ! estimate of the collisional timescale
      final_dtime = min(final_dtime, ver_dtime)

      return
   end subroutine time_step

   ! estimation of dt because of vertical settling (it's always faster than radial drift)
   subroutine vertical(swrm, dtime) ! TODO
      implicit none
      type(swarm), dimension(:), allocatable, target :: swrm
      real, dimension(:), allocatable                :: stk            ! Stokes number
      real, intent(out)                              :: dtime          ! resulting time step
      real, parameter                                :: dtparam = 0.4  ! arbitrary
      integer                                        :: i

      allocate(stk(size(swrm)))

      dtime = 1. / (eta * omegaK(swrm(1)%rdis))

      do i = 1, size(swrm)
         if (swrm(i)%rdis > smallr) dtime = min(dtime, dtparam / (omegaK(swrm(i)%rdis) * swrm(i)%stnr / (1. + swrm(i)%stnr**2.)) )
      enddo

      deallocate(stk)

      return
   end subroutine vertical

end
