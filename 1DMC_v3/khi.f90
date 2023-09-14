! this module calculates the effect of the shear instability that is then added in advection
module khi
   use constants,   only: pi, mH2, AH2, third
   use grid,        only: g
   use types,       only: swarm, list_of_swarms
   use initproblem, only: mswarm
   use discstruct,  only: sigmag, omegaK, densg, pg
   
   implicit none
   
   private
   public      :: alpha_KH
   
   contains
   
   subroutine alpha_KH(nr, rbin, alphaKH, realtime)
      integer, intent(in)                                     :: nr         ! which radial zone are we in?
      type(list_of_swarms), dimension(:), allocatable         :: rbin       ! all the particles we have in this zone are rbin(nr)
      real, dimension(:), allocatable                         :: alphaKH    ! alpha coming from the sher insta for this zone 
      real                                                    :: Z, eta, meanSt
      real, intent(in)                                        :: realtime
      integer                                                 :: i
      integer                                                 :: nparts
      real, parameter                                         :: C1 = 1.0
      real, parameter                                         :: C2 = 1.6
      real, parameter                                         :: Ceff = 0.19
      
      nparts = size(rbin(nr)%p(:))
      meanSt = 0.
      
      do i = 1, nparts
         meanSt = meanSt + rbin(nr)%p(i)%stnr
      enddo
      meanSt = meanSt / real(nparts)
      meanSt = min(meanSt,0.5)
      
      eta = 0.25 * (Pg(g%rce(nr)+1., 0.0, realtime) - Pg(g%rce(nr)-1.,0.0, realtime)) / &
           densg(g%rce(nr),0.0,realtime) / omegaK(g%rce(nr))**2. / g%rce(nr)
      Z = real(nparts)*mswarm / (pi*g%rup(nr)**2-pi*g%rlo(nr)**2) / sigmag(g%rce(nr),realtime)
      
      alphaKH(nr) = ((C1 * Ceff * abs(eta) * Z)**(-2.*third) + (C2 * Ceff * abs(eta) / Z)**(-2.))**(-1.) * meanSt
   
      return   
   end subroutine alpha_KH

end
