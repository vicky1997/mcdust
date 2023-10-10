! this module should contain all the initial conditions for dust
module initproblem

   use constants,  only: pi, third, AU
   use discstruct, only: gasmass, cs, omegaK
   use parameters, only: dtg, minrad0, maxrad0, r0, matdens
   use types,      only: swarm
   use advection,  only: stokesnr
   
   implicit none
   
   private
   public   :: init_swarms, init_random_seed, mswarm, m0, nord, nmonom0
   
   real            :: mswarm, m0, nord, nmonom0
   
   contains
   
   ! initializing the swrm array and some variables
   subroutine init_swarms(Ntot, swrm)
      implicit none
      type(swarm), dimension(:), allocatable, target   :: swrm        ! list of all swarms in the simulation
      integer, intent(in)                              :: Ntot        ! number of all swarms
      real                                             :: mdust       ! mass of the dust in the simulatated domain
      real, dimension(2)                               :: rand        ! random numbers
      real                                             :: Hg          ! pressure height scale of the gas
      integer                                          :: i
      real, parameter                                  :: s = 0.0   ! to initialize the particles density r slope (for MMSN s = -0.25)

      ! total mass of dust =  dust to gas ratio x mass of the gas               
      mdust = dtg * gasmass(minrad0*AU,maxrad0*AU,0.0)
      ! mass of one swarm
      mswarm = mdust / real(Ntot)
      ! monomer mass
      m0 = 4. * third * pi * r0**3 * matdens
      ! orders of magnitude in mass possible to cover in the simulation
      nord = (log10(mswarm/m0))
      ! initial number of monomers
      nmonom0 = mdust/m0
      
      if (.not.allocated(swrm)) allocate( swrm(Ntot) )

      ! initializing the particles
      do i = 1, Ntot
         swrm(i)%idnr = i
         swrm(i)%mass = m0 
         swrm(i)%npar = mswarm / swrm(i)%mass
         call random_number(rand)
        ! swrm(i)%rdis = minrad0 * AU + (maxrad0 - minrad0) * AU * (real(i)-0.5) / real(Ntot)   ! BE CAREFUL WITH 1D vertical column tests!!!!!!; you don't want the s then
        ! swrm(i)%rdis = minrad0 * AU + (maxrad0 - minrad0) * AU * ((s+1.)*rand(1))**(1./(s+1.)) 
         swrm(i)%rdis = (((maxrad0*AU)**(s+1.) - (minrad0*AU)**(s+1.))*rand(1) + (minrad0*AU)**(s+1.))**(1./(s+1.))
         Hg = cs(swrm(i)%rdis) / omegaK(swrm(i)%rdis)
         call random_number(rand)
         swrm(i)%zdis = Hg * sqrt(-2.*log(rand(1))) * cos(2.*pi*rand(2))
         swrm(i)%stnr = stokesnr(swrm(i), 0.0)
         swrm(i)%velr = 0.0
         swrm(i)%velz = 0.0
         swrm(i)%coll_f = 0
      enddo
          
      return
   end subroutine init_swarms
   
   ! initialize the random number generator
   subroutine init_random_seed
      implicit none
      integer                            :: i, n, clock
      integer, dimension(:), allocatable :: seed

      call random_seed(size = n)
      allocate(seed(n))

      call system_clock(count=clock)

      seed = 37 * [(i - 1, i = 1, n)]
      seed = seed + clock

      call random_seed(put = seed)
      deallocate(seed)

   end subroutine init_random_seed

end
