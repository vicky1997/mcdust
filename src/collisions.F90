! this module performes collisions between representative bodies
! in this module all the quantities should be computed at the center of current cell
! in the matrices: first index -> representative particle, 2nd -> physical
module collisions

   use constants,  only: pi, third, mH2, AH2
   use discstruct, only: cs, omegaK, densg, Pg, vgas, alpha
   use grid,       only: g
   use initproblem,only: m0, mswarm
   use parameters, only: matdens, dmmax, con1, con2, vfrag 
#ifdef EROSION
   use parameters, only: erosion_mass_ratio
#endif
   use types,      only: swarm, list_of_swarms

   implicit none

   private
   public :: mc_collisions

#ifdef TESTCOLLISIONS
   public :: col_rates, rel_vels, col_rates_r, vel_vs_centr, vel_rd_centr, col_accel, stokes_nr_centr
#endif

   contains

   ! the routine performes collisional evolution on swarms located in the cell nr,ni with the aid of the MC algorithm
   subroutine mc_collisions(nr, ni, bin, swrm, dtime, realtime, kk, first_idx)
      implicit none
      type(swarm), dimension(:), allocatable, target            :: swrm
      type(list_of_swarms), dimension(:,:), allocatable, target :: bin
      type(swarm), dimension(:), pointer              :: swarms      ! local rps array
      integer, intent(in)                             :: nr, ni      ! indices of the cell
      real, intent(in)                                :: dtime       ! time step
      real, intent(in)                                :: realtime    ! physical time
      integer, intent(in)                             :: first_idx   ! first index of swrm in rbin
      integer                                         :: nsws        ! number of representative particles in given cell
      integer                                         :: nri, nrk    ! indices of physical and representative particles choosen to the next collision
      real, dimension(:,:), allocatable               :: colrates    ! collision rates matrix
      real, dimension(:,:), allocatable               :: relvels     ! relative velocities matrix
      real, dimension(:,:), allocatable               :: accelncol   ! coagulation acceleration matrix (in the case of high mass ratio,
                                                                     ! instead of performing every collision separately, we group the collisions)
      real, dimension(:), allocatable                 :: stokesnr    ! stokes numbers of particles
      real, dimension(:), allocatable                 :: vs, vr      ! vertical settilng and radial drift velocities
      real, dimension(:), allocatable                 :: colri       ! collision rates for rps
      real                                            :: local_time, dt ! physical time spent in the cell, time step between two collisions
      real                                            :: rand        ! random number
      real                                            :: totr        ! total collision rate
      real                                            :: vn          ! maximum radial velocity from the pressure gradient
      real                                            :: Reynolds, v0, Vg2, veta, tL, teta ! relative velocities stuff
      real                                            :: lmfp, gasdens ! mean free path, gas density
      integer                                         :: i, k, l
      integer, intent(out)                            :: kk ! collisions counter
      real                                            :: dustdens
      real                                            :: Kepler_freq, cs_speed
      real:: deltar, deltaz   ! radial and vertical size of grid (for adaptative dmmax)
      
      deltar = g%rup(nr) - g%rlo(nr)
      deltaz = g%zup(nr,ni) - g%zlo(nr,ni)


      ! calculation of some values needed for calculations
      gasdens = densg(g%rce(nr),g%zce(nr,ni),realtime)         ! gas density in the center of cell
      lmfp = mH2 / ( gasdens * AH2 )                           ! mean free path in gas in the center of cell
      Kepler_freq = omegaK(g%rce(nr))                          ! keplerian frequency at the radial centre of cell
      cs_speed = cs(g%rce(nr))
      Reynolds = sqrt(0.5 * pi) * alpha(g%rce(nr)) * cs_speed / (Kepler_freq * lmfp) ! Reynolds number
      v0 = sqrt(alpha(g%rce(nr))) * cs_speed                                 ! velocity of the largest eddy
      veta = v0 * Reynolds**(-0.25)                                               ! velocity of the smallest eddy
      tL = 1. / Kepler_freq
      teta = Reynolds**(-0.5) * tL                                                ! overturn time of the smallest eddy
      vn = 0.25 * (Pg(g%rce(nr)+1., g%zce(nr,ni),realtime) - Pg(g%rce(nr)-1., g%zce(nr,ni), realtime)) / &  ! maximum radial velocity from the pressure gradient
              gasdens / Kepler_freq
      Vg2 = 1.5 *v0**2.0                                                          ! turbulent velocity of gas squared

      !-------------------------------------------------------------------------------------

      swarms => bin(nr,ni)%p       ! points to the swarms array that contains only swarms located in the current cell nr,ni
      swarms(:)%coll_f = 0         ! setting collision flag to zero, will be updated to 1 if collisions happen
      nsws = size(swarms)          ! number of swarms in the current cell

      dustdens =  real(nsws) * mswarm / g%vol(nr,ni)

      allocate (colrates(nsws,nsws), accelncol(nsws,nsws), relvels(nsws,nsws), stokesnr(nsws), vs(nsws), vr(nsws), colri(nsws))

      ! calculates initial Stokes number for particles and their settling and radial drift velocities
      do i = 1, nsws
         call stokes_nr_centr(nr, i, swarms, stokesnr, lmfp, gasdens, Kepler_freq, cs_speed)
         call vel_vs_centr(nr, ni, i, stokesnr, Kepler_freq, vs)
#ifdef RADRIFT
         call vel_rd_centr(nr, i, stokesnr, vr, vn, realtime)
#endif
      enddo

      ! calculates relative velocities and collision rates matrix
      do i = 1, nsws
         call rel_vels(nr, i, swarms, stokesnr, vr, vs, relvels, vn, Reynolds, veta, Vg2, tL, teta, Kepler_freq, cs_speed)
         call col_rates(i, swarms, relvels, colrates, accelncol, g%vol(nr,ni), deltar, deltaz)
      enddo

      ! collision rate for representative particles
      colri(:) = sum(colrates(:,:),dim=2)

      !------------ MAIN COLLISIONS LOOP ----------------------------------------------------
      local_time = 0.0
      kk = 0
      do while (local_time < dtime)
         totr = sum(colri)                 ! total collision rate
         call random_number(rand)
         dt = -1. * log(rand) / totr       ! time step between this and the next collision
         if (dt > dtime) then ! 0 or 1 collisions, decided by a random number
           call random_number(rand)
           if (rand> dtime/dt) then
             local_time = dtime
             cycle
           endif
         endif
         local_time = local_time + dt      ! update of the local time
         call choose_swarms(nri,nrk,colrates,colri,totr)
         call collision(nri,nrk,swarms,relvels, accelncol)
         call stokes_nr_centr(nr, nri, swarms, stokesnr, lmfp, gasdens, Kepler_freq, cs_speed)
         call vel_vs_centr(nr, ni, nri, stokesnr, Kepler_freq, vs)
#ifdef RADRIFT
         call vel_rd_centr(nr, nri, stokesnr, vr, vn, realtime)
#endif
         call rel_vels(nr, nri, swarms, stokesnr, vr, vs, relvels, vn, Reynolds, veta, Vg2, tL, teta, Kepler_freq,cs_speed)
         relvels(:,nri) = relvels(nri,:)
         colri(:) = colri(:) - colrates(:,nri)
         call col_rates(nri, swarms, relvels, colrates, accelncol, g%vol(nr,ni), deltar, deltaz)
         call col_rates_r(nri, swarms, relvels, colrates, accelncol, g%vol(nr,ni), deltar, deltaz)
         colri(:) = colri(:) + colrates(:,nri)
         colri(nri) = sum(colrates(nri,:))
         kk = kk + 1
      enddo
      !-------------------------------------------------------------------------------------
      write(*,*) '       collisions in zone',nr,ni,'done: total',kk,'collisions'

      deallocate (colrates, accelncol, relvels, stokesnr, vs, vr, colri)

      write(*,*) '      Updating swrm:'            ! TODO: update only the modified swarms
      do k = 1, nsws
         if(swarms(k)%coll_f /= 0) then
            l = first_idx
            l = FINDLOC(swrm(:)%idnr,swarms(k)%idnr,dim=1)
            swrm(l) = swarms(k)
         else
            cycle
         endif   
      enddo
      write(*,*) '       swrm updated!'
      nullify(swarms)

      return
   end subroutine mc_collisions

   ! calculating Stokes numbers of particle "i" IN THE CENTER OF CELL nr,ni
   subroutine stokes_nr_centr(nr, i, swarms, stokesnr, lmfp, gasdens, Kepler_freq, cs_speed)
      implicit none
      type(swarm), dimension(:)                       :: swarms
      integer, intent(in)                             :: nr
      integer, intent(in)                             :: i
      real, dimension(:), allocatable                 :: stokesnr
      real, intent(in)                                :: lmfp, gasdens, Kepler_freq, cs_speed
      real                                            :: rad

      rad = con2 * swarms(i)%mass**third                       ! particle radius

      if (rad > 2.25 * lmfp) then ! Stokes regime
         stokesnr(i) = sqrt(2.*pi) * matdens * AH2 * rad**2. * Kepler_freq / (9. * mH2 * cs_speed )
      else                        ! Epstein regime
         stokesnr(i) = rad * matdens / (sqrt(8./pi) * cs_speed  * gasdens) * Kepler_freq
      endif

      return
   end subroutine stokes_nr_centr

   ! calculation of relative velocities of bodies:
   ! we take 5 sources:
   ! Brownian motion vB
   ! turbulence vT (Ormel & Cuzzi 2007), implementation stolen from Til Birnstiel
   ! radial drift vr
   ! vertical settling vs
   ! azimuthal drift vtan
   subroutine rel_vels(nr, nl, swarms, stokesnr, vr, vs, relvels, vn, Reynolds, veta, Vg2, tL, teta, Kepler_freq, cs_speed)
      implicit none
      type(swarm), dimension(:), pointer              :: swarms
      integer, intent(in)                             :: nr
      integer, intent(in)                             :: nl
      real, dimension(:,:), allocatable               :: relvels
      real, dimension(:), allocatable                 :: stokesnr, vr, vs
      real, intent(in)                                :: vn
      real, dimension(size(swarms))                   :: vB2, vT2, vtan
      real                                            :: gts, lts
      real, intent(in)                                :: tL, teta
      real, intent(in)                                :: Reynolds
      real, intent(in)                                :: veta, Vg2
      integer                                         :: i
      real                                            :: St1, St2
      real                                            :: y_star, c1, c2, c3, c0, eps, hulp1, hulp2
      real, parameter                                 :: ya = 1.6
      real, intent(in)                                :: Kepler_freq, cs_speed

      ! Brownian motions
      vB2(:) = 8.* mH2 * cs_speed**2 * (swarms(nl)%mass+swarms(:)%mass) / (pi*swarms(nl)%mass*swarms(:)%mass)

      ! turbulence
      do i = 1, size(swarms)
         if (stokesnr(i) > stokesnr(nl)) then
            gts = stokesnr(i) / Kepler_freq
            lts = stokesnr(nl) / Kepler_freq
            St1 = stokesnr(i)
            St2 = stokesnr(nl)
         else
            gts = stokesnr(nl) / Kepler_freq
            lts = stokesnr(i) / Kepler_freq
            St1 = stokesnr(nl)
            St2 = stokesnr(i)
         endif

         if (gts < 0.2*teta) then
             vT2(i) = 1.5 *(veta/teta *(gts - lts))**2.0

         elseif (gts < teta/ya) then
             vT2(i) = Vg2 *(St1-St2)/(St1+St2)*(St1**2.0/(St1+Reynolds**(-0.5)) - St2**2.0/(St2+Reynolds**(-0.5)))
         elseif (gts < 5.0*teta) then
            !Eq. 17 of OC07. The second term with St_i**2.0 is negligible (assuming !Re>>1)
            !hulp1 = Eq. 17; hulp2 = Eq. 18

            hulp1 = ( (St1-St2)/(St1+St2) * (St1**2.0/(St1+ya*St1) - St2**2.0/(St2+ya*St1)) )!note the -sign
            hulp2 = 2.0*(ya*St1-Reynolds**(-0.5)) + St1**2.0/(ya*St1+St1) - St1**2.0/(St1+Reynolds**(-0.5)) +&
                                           St2**2.0/(ya*St1+St2) - St2**2.0/(St2+Reynolds**(-0.5))
            vT2(i) = Vg2 *(hulp1 + hulp2)

         elseif (gts < tL*0.2)  then
            eps=St2/St1!stopping time ratio
            vT2(i) = Vg2 *( St1*(2.0*ya - (1.0+eps) + 2.0/(1.0+eps) *(1.0/(1.0+ya) + eps**3.0/(ya+eps) )) )

         elseif (gts < tL) then
           !now y* lies between 1.6 (St1 << 1) and 1.0 (St1>=1). The fit below fits ystar to less than 1%
           c3 =-0.29847604
           c2 = 0.32938936
           c1 =-0.63119577
           c0 = 1.6015125
           y_star = c0 + c1*St1 + c2*St1**2.0 + c3*St1**3.0
           !we can then employ the same formula as before
           eps=St2/St1
           vT2(i) = Vg2 *( St1*(2.0*y_star - (1.0+eps) + 2.0/(1.0+eps) *(1.0/(1.0+y_star) + eps**3.0/(y_star+eps) )) )

         else
            vT2(i) = Vg2 *( 1.0/(1.0+St1) + 1.0/(1.0+St2) )
         endif

      enddo

      ! tangential
      vtan(:) = vn * ( 1. /(1.+stokesnr(:)**2.) - 1. / (1.+stokesnr(nl)**2.) )

      ! total
#ifdef RADRIFT
         relvels(nl,:) = sqrt(vB2(:)  + vT2(:) + (vs(:) - vs(nl))**2 + (vr(:) - vr(nl))**2 + vtan(:)**2)
#else
         relvels(nl,:) = sqrt(vB2(:) + (vs(:) - vs(nl))**2 + (vr(:) - vr(nl))**2 + vtan(:)**2)
#endif

      return
   end subroutine rel_vels

   ! calculation of the collision rates between representative particle nl and all physical particles
   subroutine col_rates(nl, swarms, relvels, colrates, accelncol, vol, deltar, deltaz)
      implicit none
      integer, intent(in)                             :: nl       ! number of line of colrates to update
      type(swarm), dimension(:), pointer              :: swarms
      real, dimension(:,:), allocatable               :: colrates
      real, dimension(:,:), allocatable               :: accelncol
      real, dimension(:,:), allocatable               :: relvels
      real, intent(in)                                :: vol      ! volume of the cell
      real, intent(in) :: deltar, deltaz ! size of radial and vertical grid

      ! basic eq. (see e.g. Drazkowska et al 2013, Eqs 19-20)
      colrates(nl,:) = swarms(:)%npar * relvels(nl,:) * con1 * (swarms(nl)%mass**third + swarms(:)%mass**third)**2./vol
      ! adaptative dmmax
      where(swarms(:)%mass/swarms(nl)%mass < dmmax)
         accelncol(nl,:) = max(colrates(nl,:)*min(deltar/abs(swarms(nl)%velr),deltaz/abs(swarms(nl)%velz)), 1.)
      elsewhere
         accelncol(nl,:) = 1. 
      endwhere
      colrates(nl,:) = colrates(nl,:) / accelncol(nl,:)
      ! if the representative particle represents less than 1 particle, the method is not valid anymore,
      ! so the collision rate is supressed
      where(swarms(:)%npar <= 1.0) colrates(nl,:) = 0.0

      return
   end subroutine col_rates

   ! for updating the collision rates after collision: calculating the collsion rates between physical particle nl
   ! and all the representative particles
   subroutine col_rates_r(nl, swarms, relvels, colrates, accelncol, vol, deltar, deltaz)
      implicit none
      integer, intent(in)                             :: nl       ! number of column of colrates to update
      type(swarm), dimension(:), pointer              :: swarms
      real, dimension(:,:), allocatable               :: colrates
      real, dimension(:,:), allocatable               :: accelncol
      real, dimension(:,:), allocatable               :: relvels
      real, intent(in)                                :: vol      ! volume of the cell
      real, intent(in) :: deltar, deltaz ! size of radial and vertical grid

      colrates(:,nl) = swarms(nl)%npar * relvels(nl,:) * con1 * (swarms(nl)%mass**third + swarms(:)%mass**third)**2./vol
      
      ! adaptative dmmax
      where(swarms(nl)%mass/swarms(:)%mass < dmmax)
         accelncol(:,nl) = max(colrates(:,nl) * min(deltar/abs(swarms(:)%velr),deltaz/abs(swarms(:)%velz)), 1.)
      elsewhere
         accelncol(:,nl) = 1. 
      endwhere
      colrates(:,nl) = colrates(:,nl) / accelncol(:,nl)
      where(swarms(:)%npar <= 1.0) colrates(:,nl) = 0.0

      return
   end subroutine col_rates_r

   ! calculating the number of collisions in the case of collision acceleration
   ! only in the case when the representative particle is much larger than the physical particle
   subroutine col_accel(nri, swarms, accelncol)
      implicit none
      integer, intent(in)                             :: nri
      real, dimension(:), pointer                     :: accelncol_c
      type(swarm), dimension(:), pointer              :: swarms
      real, allocatable, dimension(:,:), target       :: accelncol

      accelncol_c => accelncol(nri,:)

      accelncol_c(:) = 1.0
      where ((swarms(:)%mass / swarms(nri)%mass) < dmmax)
         accelncol_c(:) = swarms(nri)%mass * dmmax / swarms(:)%mass
      endwhere

      return
   end subroutine col_accel

   subroutine col_accel_r(nri, swarms, accelncol)
      implicit none
      integer, intent(in)                             :: nri
      real, dimension(:), pointer                     :: accelncol_r
      type(swarm), dimension(:), pointer              :: swarms
      real, allocatable, dimension(:,:), target       :: accelncol

      accelncol_r => accelncol(:,nri)

      accelncol_r(:) = 1.0
      where ((swarms(nri)%mass / swarms(:)%mass) < dmmax)
         accelncol_r(:) = swarms(:)%mass * dmmax / swarms(nri)%mass
      endwhere

      return
   end subroutine col_accel_r

   ! choosing particles to the next collision
   ! nri -> representative
   ! nrk -> physical
   subroutine choose_swarms(nri, nrk, colrates, ri, totrate)
      implicit none
      integer, intent(out)                            :: nri, nrk
      real, dimension(:,:), allocatable, intent(in)   :: colrates
      integer                                         :: i
      real, dimension(:), allocatable, intent(in)     :: ri
      real, intent(in)                                :: totrate
      real, dimension(2)                              :: rand
      real                                            :: fin

      call random_number(rand) ! drawing 2 random numbers

      ! choosing the representative particle: the higher ri, the higher chance that it's particle "i"
      rand(1) = rand(1) * totrate
      i = 1
      fin = ri(1)
      do while(rand(1) > fin)
         fin = fin + ri(i+1)
         i = i + 1
      enddo
      nri = i

      ! choosing the physical particle
      i = 1
      rand(2) = rand(2) * ri(nri)
      fin = colrates(nri,1)
      do while (rand(2) > fin)
         fin = fin + colrates(nri,i+1)
         i = i + 1
      enddo
      nrk = i

      return
   end subroutine choose_swarms

   ! performing the collision: deciding the collision outcome - put your collision model here
   ! only the representative particle is updated
   subroutine collision(nri,nrk,swarms,relvels, accelncol)
      implicit none
      integer, intent(in)                             :: nri, nrk
      type(swarm), dimension(:), pointer              :: swarms
      real, dimension(:,:), allocatable               :: accelncol
      real, dimension(:,:), allocatable               :: relvels
      real                                            :: rvel

      rvel = relvels(nri,nrk)

      if (rvel < vfrag) then
         call hit_and_stick(nri,nrk,swarms, accelncol)
#ifdef EROSION
      else if (swarms(nri)%mass/swarms(nrk)%mass .ge. erosion_mass_ratio) then
         call erosion(nri,nrk,swarms,accelncol)
#endif
      else
         call fragmentation(nri, swarms)
      endif
      swarms(nri)%npar = mswarm / swarms(nri)%mass
      swarms(nri)%coll_f =  1
      return
   end subroutine collision

   ! sticking collision
   subroutine hit_and_stick(nri,nrk,swarms, accelncol)
      implicit none
      integer, intent(in)                             :: nri, nrk
      type(swarm), dimension(:), pointer              :: swarms
      real, dimension(:,:), allocatable               :: accelncol
      
      !if (accelncol(nri, nrk) > 0.0) then
      !   swarms(nri)%mass = swarms(nri)%mass + accelncol(nri, nrk) * swarms(nrk)%mass
      !else
      !   swarms(nri)%mass = swarms(nri)%mass + swarms(nrk)%mass
      !endif
      swarms(nri)%mass = swarms(nri)%mass + accelncol(nri, nrk) * swarms(nrk)%mass
      return
   end subroutine hit_and_stick

   ! fragmentation collision
   ! put your fragment size distribution here:
   ! n(m) ~ m^(kappa - 2)
   subroutine fragmentation(nri,swarms)
      implicit none
      integer, intent(in)                             :: nri
      real                                            :: ran
      type(swarm), dimension(:), pointer              :: swarms
      real, parameter                                 :: kappa = 1./6.  ! n(m) ~ m^(kappa - 2)

      call random_number(ran)
      swarms(nri)%mass = (ran * (swarms(nri)%mass**kappa - m0**kappa) +  m0**kappa )**(1./kappa)
      swarms(nri)%mass = max(swarms(nri)%mass, m0)

      return
   end subroutine fragmentation

#ifdef EROSION
   !erosion collision as in dustpy
   !when small aggregate hits large one, it fragments
   !and big aggregate loses a chunk mass same as target particle

   subroutine erosion(nri, nrk, swarms, accelncol)
      implicit none
      integer, intent(in)                             :: nri, nrk
      type(swarm), dimension(:), pointer              :: swarms
      real, dimension(:,:), allocatable               :: accelncol
      real                                            :: ran, p
      

      p = swarms(nrk)%mass/swarms(nri)%mass
      call random_number(ran)
      if (ran .ge. p) then
            !if(accelncol(nri,nrk) > 0.0) then
            !   swarms(nri)%mass = swarms(nri)%mass - accelncol(nri, nrk) * swarms(nrk)%mass
            !else
            !   swarms(nri)%mass = swarms(nri)%mass - swarms(nrk)%mass
            !endif
            swarms(nri)%mass = swarms(nri)%mass - accelncol(nri, nrk) * swarms(nrk)%mass
      else
         swarms(nri)%mass = swarms(nrk)%mass
      endif
      return
   end subroutine erosion
#endif

   ! vertical settling velocity
   subroutine vel_vs_centr(nr, ni, i, stokesnr, Kepler_freq, vs)
      implicit none
      integer, intent(in)                             :: nr, ni
      real, dimension(:), allocatable, intent(in)     :: stokesnr
      real, dimension(:), allocatable                 :: vs
      integer, intent(in)                             :: i       ! index of particle
      real, intent (in)                               :: Kepler_freq
      vs(i) = g%zce(nr,ni) * Kepler_freq * min(stokesnr(i), 0.5)

      return
   end subroutine vel_vs_centr

#ifdef RADRIFT
   ! velocity of radial drift
   subroutine vel_rd_centr(nr, i, stokesnr, vr, vn, realtime)
      implicit none
      integer, intent(in)                             :: nr
      real, dimension(:), allocatable                 :: stokesnr
      real, dimension(:), allocatable                 :: vr
      real, intent(in)                                :: vn
      integer, intent(in)                             :: i      ! index of particle
      real, intent(in)                                :: realtime

      vr(i) = 2. * vn * stokesnr(i) / (stokesnr(i)**2. + 1.) + &
              vgas(g%rce(nr), realtime) / (1. + stokesnr(i) * stokesnr(i))

      return
   end subroutine vel_rd_centr
#endif

end
