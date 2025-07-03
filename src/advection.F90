! this module performes advection of every representative particle due to radial drift, vertical settling and turbulent
! diffusion
! in this module all the quantities should be calculated locally (at the particle location) therefore we don't need the
! grid to perform advection
module advection

   use constants,  only: pi, AU, third, AH2, mH2, year
   use discstruct, only: Pg, densg, omegaK, cs, diffcoefgas, vgas, ddensgdr, ddensgdz, alpha
   use types,      only: swarm
   use parameters, only: matdens, smallr, con2

   implicit none

   private
   public      :: mc_advection, stokesnr, update_St
#ifdef RADRIFT
   public      :: vel_rd, vel_vn
#endif

#ifdef VERSETTLING
   public      :: vel_vs
#endif
   contains

   ! routine performes radial and vertical advection due to radial drift, vertical settling and turbulent diffusion
   ! for every particle in the simulation
   subroutine mc_advection(swrm, dtime, realtime)
      implicit none
      type(swarm), dimension(:), allocatable, target :: swrm
      type(swarm), pointer                           :: particle
      real, intent(in)                               :: dtime
      real, intent(in)                               :: realtime
      real                                           :: vs, vr, velr, velv, vn
      integer                                        :: i

      ! loop over all the particles: calculating advection velocities
      
      do i = 1, size(swrm)
         if (swrm(i)%rdis < smallr) cycle ! if particle is inside photoevaporation radius, no advection

         particle => swrm(i)

         particle%stnr = stokesnr(particle,realtime)
#ifdef RADRIFT
         call vel_vn(particle, vn, realtime)
#endif
         
#ifdef VERSETTLING
         call vel_vs(particle, vs)
         call vel_ver(particle, vs, velv, dtime, realtime)
         particle%zdis = particle%zdis + velv * dtime
#else
         call vertical_redistr(particle, realtime)
#endif
#ifdef RADRIFT
         call vel_rd(particle, vr, vn, realtime)
         call vel_rad(particle, velr, vr, dtime, realtime)
         particle%rdis = particle%rdis + velr * dtime
         if (particle%rdis > smallr) then
            particle%stnr = stokesnr(particle,realtime)
         else
            particle%rdis = 0.99 * smallr
            particle%velr = 1.e-20
            particle%velz = 1.e-20
         endif
#endif
      enddo
      

   end subroutine mc_advection

   ! the maximum radial drift speed due to the gas pressure gradient
   subroutine vel_vn(particle, vn, realtime)
      implicit none
      type(swarm)                                    :: particle
      real, intent(in)                               :: realtime
      real                                           :: vn

      vn = 0.25 * (Pg(particle%rdis+1., particle%zdis, realtime) - Pg(particle%rdis-1.,particle%zdis, realtime)) / &
              densg(particle%rdis,particle%zdis,realtime) / omegaK(particle%rdis)

      return
   end subroutine vel_vn

#ifdef VERSETTLING
   ! the complete vertical velocity
   subroutine vel_ver(particle, vs, velv, dtime, realtime)
      implicit none
      type(swarm)                                     :: particle
      real, intent(in)                                :: dtime, realtime
      real, intent(in)                                :: vs
      real                                            :: velv
      real                                            :: dz
      real                                            :: Ldiff
      real, dimension(2)                              :: rand
      real                                            :: vD1, vD2


      Ldiff = sqrt(2. * dtime * diffcoefgas(particle%rdis)  / (1. + particle%stnr**2))
      call random_number(rand)
      dz = Ldiff/(sqrt(2.*log(2.))) * sqrt(-2.*log(rand(1))) * cos(2.*pi*rand(2))
      vD1 = dz / dtime
      vD2 = (diffcoefgas(particle%rdis) / (1. + particle%stnr**2) / densg(particle%rdis,particle%zdis, realtime)) * &
         ddensgdz(particle%rdis,particle%zdis, realtime)
      velv = - vs + vD1 + vD2
   end subroutine vel_ver
#endif

#ifdef RADRIFT
   ! the complete radial velocity
   subroutine vel_rad(particle, velr, vr, dtime, realtime)
      implicit none
      type(swarm)                                     :: particle
      real, intent(in)                                :: dtime, realtime
      real, intent(in)                                :: vr
      real                                            :: velr
      real                                            :: Ldiff
      real, dimension(2)                              :: rand
      real                                            :: dr
      real                                            :: vD1, vD2


      Ldiff = sqrt(2. * dtime * diffcoefgas(particle%rdis) / (1. + particle%stnr**2))
      call random_number(rand)
      dr = Ldiff/(sqrt(2.*log(2.))) * sqrt(-2.*log(rand(1))) * cos(2.*pi*rand(2))
      vD1 = dr / dtime
      vD2 = (diffcoefgas(particle%rdis) / (1. + particle%stnr**2) / densg(particle%rdis,particle%zdis, realtime)) * &
         (ddensgdr(particle%rdis,particle%zdis, realtime))
      velr = vr + vD1 + vD2


      return
   end subroutine vel_rad
    ! vertical redistribution of particles to a theoretical Gaussian profile
   subroutine vertical_redistr(particle, realtime)
      implicit none
      type(swarm)                                     :: particle
      real                                            :: Hd
      real                                            :: rand1, rand2
      real                                            :: Hg
      real, intent(in)                                :: realtime


      Hg = cs(particle%rdis) / omegaK(particle%rdis)
      Hd = Hg * (alpha(particle%rdis) / (alpha(particle%rdis)+particle%stnr))**(0.5)

      call random_number(rand1)
      call random_number(rand2)
      particle%zdis = Hd * sqrt(-2.*log(rand1)) * cos(2.*pi*rand2)

      return
   end subroutine vertical_redistr
#endif
   ! calculates the Stokes numbers of particles locally
   real function stokesnr(particle, realtime)
      implicit none
      type(swarm)                                     :: particle
      real                                            :: lmfp
      real                                            :: rad
      real                                            :: gasdens, css
      real, intent(in)                                :: realtime

      gasdens = densg(particle%rdis,particle%zdis,realtime)
      css = cs(particle%rdis)
      lmfp = mH2 / ( gasdens * AH2 )
      rad = con2 * particle%mass**third

      if (rad > 2.25 * lmfp) then ! Stokes regime
         stokesnr = sqrt(2.*pi) * rad**2. * matdens * omegaK(particle%rdis) * AH2 / (9. * css * mH2)
      else                 ! Epstein regime
         stokesnr = rad * matdens / (sqrt(8./pi) * css * gasdens) * omegaK(particle%rdis)
      endif
      if(isnan(stokesnr)) then
         write(*,*)'in St w id',particle%idnr,'r',particle%rdis/AU,'z', particle%zdis/AU,'mass', particle%mass
         write(*,*)'dens', gasdens, 'omegaK', omegaK(particle%rdis),'size', rad,  'cs', css, 'con2', con2
         write(*,*)'npar', particle%npar
         stop
      endif 
      return
   end function stokesnr

#ifdef VERSETTLING
   ! vertical settling velocity
   subroutine vel_vs(particle, vs)
      implicit none
      type(swarm)                                     :: particle
      real                                            :: vs

      vs =particle%zdis * omegaK(particle%rdis) * particle%stnr / (1. + particle%stnr**2.)
      particle%velz = vs
      return
   end subroutine vel_vs
#endif

#ifdef RADRIFT
   ! velocity of the radial drift
   subroutine vel_rd(particle, vr, vn, realtime)
      implicit none
      type(swarm)                                     :: particle
      real, intent(in)                                :: realtime
      real                                            :: vr
      real, intent(in)                                :: vn

      vr = 2. * vn / (particle%stnr + 1./particle%stnr) + vgas(particle%rdis, realtime) / (1. + particle%stnr* particle%stnr)
      particle%velr = vr

      return
   end subroutine vel_rd
#endif

subroutine update_St(swrm, realtime)
   implicit none
   type(swarm), dimension(:), allocatable, target :: swrm
   type(swarm), pointer                           :: particle
   real, intent(in)                               :: realtime
   integer                                        :: i


   do i = 1, size(swrm)
      if (swrm(i)%rdis < smallr) cycle
      particle => swrm(i)
      particle%stnr = stokesnr(particle,realtime)
   enddo

   return
end subroutine

end
