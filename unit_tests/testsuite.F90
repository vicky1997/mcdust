program test_advection
    use advection
    use constants
    use parameters
    use types
    use collisions
    use discstruct
    use grid, only: make_grid, g

    implicit none


#ifdef TESTCOLLISIONS
    type(swarm), dimension(:), allocatable, target    :: swrm        ! list of all swarms in the simulation
    type(swarm), dimension(:), pointer                :: swarms        ! list of all swarms in the simulation
    type(list_of_swarms), dimension(:,:), allocatable, target       :: bin         ! swarms binned into cells
    type(list_of_swarms), dimension(:), allocatable                 :: rbin        ! swarms binned into radial zones
    integer                                          :: Ntest = 256  ! number of all swarms
    real, dimension(:,:), allocatable               :: colrates    ! collision rates matrix
    real, dimension(:,:), allocatable               :: relvels     ! relative velocities matrix
    real, dimension(:,:), allocatable               :: accelncol   ! coagulation acceleration matrix (in the case of high mass ratio,
                                                                    ! instead of performing every collision separately, we group the collisions)
    real, dimension(:), allocatable                 :: stokesno    ! stokes numbers of particles
    real, dimension(:), allocatable                 :: vs, vr      ! vertical settilng and radial drift velocities
    real, dimension(:), allocatable                 :: colri       ! collision rates for rps
    real                                            :: vn          ! maximum radial velocity from the pressure gradient
    real                                            :: Reynolds, v0, Vg2, veta, tL, teta ! relative velocities stuff
    real                                            :: lmfp, gasdens ! mean free path, gas density
    integer                                         :: i
    real                                            :: Kepler_freq, cs_speed
    integer                                         :: ni=1, nrtest=1
    real                                            :: deltar, deltaz
    real                                            :: totmass
    integer, dimension(:,:), allocatable            :: ncolls
    real, dimension(5)                              :: zarray = (/-0.02, -0.01, 0., 0.01, 0.02/)
    real, dimension(5)                              :: rarray = (/9.98, 9.99, 10.00, 10.01, 10.02/)
#endif
 
#ifdef TESTADVECTION
    type(swarm)                                      :: particle, particle2
    real                                             :: vs1, vn1, vr1, realtime, vs2, vn2, vr2
    real                                             :: St1, St2, Stratio
#endif


    real                                             :: mdust       ! mass of the dust in the simulatated domain
    real                                             :: mswarm
    real                                             :: m0
    character(len=100)  :: ctrl_file           ! parameter file
    
    call get_command_argument(1, ctrl_file)
    write(*,*) 'mcdust v1.0 Testing module'
    write(*,*) '------------------------------------------------------------------'
    write(*,*) 'Reading parameters...'
    call read_parameters(ctrl_file)
    write(*,*) '------------------------------------------------------------------'
    write(*,*) 'Initializing representative bodies...'
    
    mdust = dtg * gasmass(minrad0*AU,maxrad0*AU,0.0)

    m0 = 4. * third * pi * r0**3 * matdens


#ifdef TESTCOLLISIONS

    if (.not.allocated(swrm)) allocate( swrm(Ntest) )
    mswarm = mdust / real(Ntot)
      ! initializing the particles
    do i = 1, Ntest
         swrm(i)%idnr = i
         swrm(i)%mass = m0 
         swrm(i)%npar = mswarm / swrm(i)%mass
         swrm(i)%rdis = rarray(mod(i,5)+1)*AU      
         swrm(i)%zdis= zarray(mod(i,5)+1)*AU
         swrm(i)%stnr = stokesnr(swrm(i), 0.0)
         swrm(i)%velr = 1.0
         swrm(i)%velz = 1.0
         swrm(i)%coll_f = 0
    enddo
    swarms => swrm
#endif
#ifdef TESTADVECTION
    particle%idnr = 1
    particle%mass = m0
    particle%rdis = 10.*AU
    particle%zdis = 0.1*AU
    particle%stnr = stokesnr(particle, 0.0)
    particle%velr = 0.0
    particle%velz = 0.0
    particle%coll_f = 0

    write(*,*) '------------------------'
    write(*,*) 'Testing advection module'
    write(*,*) '------------------------'
    write(*,*) 'Calculating velocities for particle with the following properties to compare and check'
    write(*,*) 'Particle 1'
    write(*,*) ' | mass (g): ', particle%mass, ' | Stokes Number (St):', particle%stnr, ' |'
    write(*,*) ' | Distance from star (AU):', particle%rdis/AU, ' | vertical height (AU)', particle%zdis/AU,' |'

   

    call vel_vs(particle, vs1)
    call vel_vn(particle, vn1, realtime)
    call vel_rd(particle, vr1, vn1, realtime)
    
    particle2%idnr = 2
    particle2%mass = 10*m0
    particle2%rdis = 10.*AU
    particle2%zdis = 0.1*AU
    particle2%stnr = stokesnr(particle2, 0.0)
    particle2%velr = 0.0
    particle2%velz = 0.0
    particle2%coll_f = 0
    St1 = particle%stnr
    St2 = particle2%stnr
    Stratio = St2*(1+St1**2)/(St1*(1+St2**2))

    write(*,*) 'Particle 2'
    write(*,*) ' | mass (g): ', particle2%mass, ' | Stokes Number (St):', particle2%stnr, ' |'
    write(*,*) ' | Distance from star (AU):', particle2%rdis/AU, ' | vertical height (AU)', particle2%zdis/AU,' |'

    call vel_vs(particle2, vs2)    
    call vel_vn(particle2, vn2, realtime)
    call vel_rd(particle2, vr2, vn2, realtime)
    
    if((abs(vs2/vs1 - Stratio)<smallv))  then
        write(*,*)'vertical setting velocity works'
    else
        write(*,*)'error in vertical settling velocity routine vel_vs'
    endif

    if((abs(vn1-vn2)<smallv)) then
        write(*,*)'maximal drift velocity vel_vn works'
    else
        write(*,*)'error in maximal drift velocity vel_vn'
    endif

    if((abs(vr2/vr1 - Stratio)<smallv)) then
        write(*,*)'total radial velocity works'
    else
        write(*,*)'error in total radial velocity routine vel_rd'
    endif

    write(*,*)'---------------------------'
    write(*,*)'Testing advection completed'
    write(*,*)'---------------------------'

#endif

#ifdef TESTCOLLISIONS
    write(*,*)'-------------------------'
    write(*,*)'Testing collisions module'
    write(*,*)'-------------------------'


    ! making a single grid
    write(*,*) 'making grid...'
    call make_grid(swrm, bin, rbin, nrtest, 1, smallr,totmass, ncolls)
    write(*,*)' Grid done! with no. of particles in the cell =', size(swrm)

    ! calculation of some values needed for calculations
    gasdens = densg(g%rce(nrtest),g%zce(nrtest,ni),realtime)         ! gas density in the center of cell
    lmfp = mH2 / ( gasdens * AH2 )                           ! mean free path in gas in the center of cell
    Kepler_freq = omegaK(g%rce(nrtest))                          ! keplerian frequency at the radial centre of cell
    cs_speed = cs(g%rce(nrtest))
    Reynolds = sqrt(0.5 * pi) * alpha(g%rce(nrtest)) * cs_speed / (Kepler_freq * lmfp) ! Reynolds number
    v0 = sqrt(alpha(g%rce(nrtest))) * cs_speed                                 ! velocity of the largest eddy
    veta = v0 * Reynolds**(-0.25)                                               ! velocity of the smallest eddy
    tL = 1. / Kepler_freq
    teta = Reynolds**(-0.5) * tL                                                ! overturn time of the smallest eddy
    vn = 0.25 * (Pg(g%rce(nrtest)+1., g%zce(nrtest,ni),realtime) - Pg(g%rce(nrtest)-1., g%zce(nrtest,ni), realtime)) / &  ! maximum radial velocity from the pressure gradient
            gasdens / Kepler_freq
    Vg2 = 1.5 *v0**2.0
    deltar = g%rup(nrtest) - g%rlo(nrtest)
    deltaz = g%zup(nrtest,ni) - g%zlo(nrtest,ni)        
    allocate (colrates(Ntest,Ntest), accelncol(Ntest,Ntest), relvels(Ntest,Ntest), stokesno(Ntest), vs(Ntest),&
                     vr(Ntest), colri(Ntest))
    do i = 1, Ntest
        call stokes_nr_centr(i, swarms, stokesno, lmfp, gasdens, Kepler_freq, cs_speed)
        call vel_vs_centr(nrtest, ni, i, stokesno, Kepler_freq, vs)
#ifdef RADRIFT
        call vel_rd_centr(nrtest, i, stokesno, vr, vn, realtime)
#endif
    enddo

    ! calculates relative velocities and collision rates matrix
    do i = 1, ntest
        call rel_vels(i, swarms, stokesno, vr, vs, relvels, vn, Reynolds, veta, Vg2, tL, teta, Kepler_freq, cs_speed)

        call col_rates(i, swarms, relvels, colrates, accelncol, g%vol(nrtest,ni), deltar, deltaz)
    enddo

    colri(:) = sum(colrates(:,:),dim=2)
     
    if(abs(sum(relvels)/size(relvels) - relvels(1,1)) < smallv) then
        write(*,*) 'Relative velocities tested successfully'
    else
        write(*,*) 'Error in relative velocities subroutine rel_vels()'
    endif

    if(abs(sum(colri)/size(colri) - colri(1)) < smallv) then
        write(*,*) 'Collision rates tested successfully'
    else
        write(*,*) 'Error in collision rates subroutine col_rates()'
    endif
    
    write(*,*) '-------------------------'
    write(*,*) 'Collision test completed'
    write(*,*) '-------------------------'

    deallocate(colrates, accelncol, relvels, stokesno, vs, vr, colri)

#endif

    write (*,*) ' End of test suite'

end