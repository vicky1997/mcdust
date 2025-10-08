! this module reads parameters from params.par file
module parameters
    
    use constants,                only: year, AU, pi, third

    implicit none

    private
    public :: read_parameters, Ntot, ncell, nr, nz, fout, dtime, tend, smallr, restart, restime, minrad0, maxrad0, &
              r0, matdens, dmmax, eta, dtg, vfrag, con1, con2, db_data, &
              datadir, path, sigmag0, temperature, alpha_t, nbins
#ifdef EROSION
    public :: erosion_mass_ratio
#endif
#ifdef LOGTIME
    public :: tstart, ntimeout
#endif
    integer                             :: Ntot     ! total number of representative particles in the simulation
    integer                             :: ncell    ! number of particles per cell
    integer                             :: nr       ! nr of radial zones
    integer                             :: nz       ! nr of zones in z
    real                                :: dtime    ! do output at least every dtime time
    real                                :: tend     ! time to finish the simulation
    real                                :: smallr   ! smallest distance (evaporation radius)
    real                                :: dtg      ! dust to gas ratio
    logical                             :: restart  ! is the simulation run from restart?
    real                                :: restime  ! time of restart
    real                                :: minrad0  ! minimum distance from the star of initial rim
    real                                :: maxrad0  ! maximum distance from the star of initial rim
    real                                :: r0       ! monomer radius
    real, protected                     :: matdens  ! material density
    real                                :: dmmax    ! MC-acceleration parameter
    real                                :: vfrag    ! fragmentation treshold velocity
    integer                             :: nbins    ! number of bins for mass histograms
    integer                             :: fout     ! steps between outputs
    real                                :: alpha_t, sigmag0, temperature, eta ! gas disk properties
#ifdef EROSION
    integer                             :: erosion_mass_ratio ! mass ratio between bodies to trigger erosion
#endif
    logical                             :: db_data  ! write out files for sanity checks?
    character(len=100)                  :: datadir  ! data directory to write the data into.
    character(len=100)                  :: path     ! absolute path of data directory
    real, protected  :: con1  ! optimization
    real, protected  :: con2  ! optimization
#ifdef LOGTIME
    real                                :: tstart
    integer                             :: ntimeout
#endif
    contains

    subroutine read_parameters(ctrl_file)
    
        implicit none
        character(len=100), intent(in)   :: ctrl_file
        character(len=100)               :: buffer, label
        integer                          :: pos
        integer                          :: ios = 0
        integer                          :: line = 0
        integer, parameter               :: fh = 1
        character(len=100)               :: command

        ! setting default values
        minrad0 = 3.0
        maxrad0 = 5.0
        r0 = 1.e-2
        ncell = 200
        nr   = 10
        nz = 1
        fout = 5
        dtime = 100. * year
        tend = 1000. * year
        matdens = 1.
        dmmax = 0.001
        smallr = 3.0
        dtg = 0.01
        vfrag = 100.
        alpha_t = 1.e-3
        sigmag0 = 800.
        temperature = 280.
        eta = 0.05
        restart = .false.
        restime = 0.0
        nbins = 200
#ifdef EROSION
        erosion_mass_ratio = 10
#endif
        db_data = .true.
        datadir = "data"

        ! reading the parameter file
      open(fh,file=ctrl_file,action='read')

      ! ios is negative if an end of record condition is encountered or if
      ! an endfile condition was detected.  It is positive if an error was
      ! detected.  ios is zero otherwise.

      do while (ios == 0)
            read(fh, '(A)', iostat=ios) buffer
            if (ios == 0) then
                line = line + 1

                pos = scan(buffer, ' 	')
                label = buffer(1:pos)
                buffer = buffer(pos+1:)

                select case (label)
                case ('minimum_radius_[AU]')
                    read(buffer, *, iostat=ios) minrad0
                    print *, 'Read minimum radius: ', minrad0
                case ('maximum_radius_[AU]')
                    read(buffer, *, iostat=ios) maxrad0
                    print *, 'Read maximum radius: ', maxrad0
                case ('monomer_radius_[cm]')
                    read(buffer, *, iostat=ios) r0
                    print *, 'Read monomer radius: ', r0
                case('number_of_particles_per_cell')
                    read(buffer, *, iostat=ios) ncell
                    print *, 'Read number of particles per cell: ', ncell
                case('number_of_radial_zones')
                    read(buffer, *, iostat=ios) nr
                    print *, 'Read number of radial zones: ', nr
                case('number_of_vertical_zones')
                    read(buffer, *, iostat=ios) nz
                    print *, 'Read number of vertical zones: ', nz
                case('time_between_outputs_[yrs]')
                    read(buffer, *, iostat=ios) dtime
                    print *, 'Read time between outputs: ', dtime
                    dtime = dtime * year
                case('steps_between_outputs')
                    read(buffer, *, iostat=ios) fout
                    print *, 'Read number of iterations between outputs: ', fout
#ifdef LOGTIME
                case('no_of_outputs')
                    read(buffer, *, iostat=ios) ntimeout
                    print *, 'Number of outputs: ', ntimeout
                
                case('output_start_time_[yr]')
                    read(buffer, *, iostat=ios) tstart
                    print *, 'Output start time: ', tstart
                    tstart = tstart*year
#endif
                case('maximum_time_of_simulation_[yrs]')
                    read(buffer, *, iostat=ios) tend
                    print *, 'Read maximum time of simulation: ', tend
                    tend = tend * year
                case('material_density_[g/cm3]')
                    read(buffer, *, iostat=ios) matdens
                    print *, 'Read material density: ', matdens
                case('dmmax')
                    read(buffer, *, iostat=ios) dmmax
                    print *, 'Read dmmax: ', dmmax
                case('evaporation_radius_[AU]')
                    read(buffer, *, iostat=ios) smallr
                    print *, 'Read evaporation radius: ', smallr
                    smallr = smallr * AU
                case('dust_to_gas_ratio')
                    read(buffer, *, iostat=ios) dtg
                    print *, 'Read dust to gas ratio: ', dtg
                case('fragmentation_velocity_[cm/s]')
                    read(buffer, *, iostat=ios) vfrag
                    print *, 'Read fragmentation velocity: ', vfrag
                case('restart')
                    read(buffer, *, iostat=ios) restart
                    print *, 'Restart? ',restart
                case('Restart_time_[yrs]')
                    read(buffer, *, iostat=ios) restime
                    print *, 'Read restart time: ', restime
                case('alpha')
                    read(buffer, *, iostat=ios) alpha_t
                    print *, 'alpha ', alpha_t
                case('sigma_gas_[g/cm2]')
                    read(buffer, *, iostat=ios) sigmag0
                    print *, 'gas surface density ', sigmag0
                case('temperature_[K]')
                    read(buffer, *, iostat=ios) temperature
                    print *, 'gas temperature ', temperature
                case('eta')
                    read(buffer, *, iostat=ios) eta
                    print *, 'pressure gradient parameter eta ', eta
#ifdef EROSION
                case('erosion_mass_ratio')
                    read(buffer, *,iostat=ios) erosion_mass_ratio
                    print *, 'erosion?', erosion_mass_ratio
#endif
                case('generate_debug_files?')
                    read(buffer, *,iostat=ios) db_data
                    print *, 'generate files for sanity check?', db_data
                case('data_directory')
                    read(buffer, *,iostat=ios) datadir
                    print *, 'data output directory', datadir

                ! if you want other parameters add another cases here and in the setup file

                case default
                    print *, 'Skipping invalid label at line', line
                end select
            end if
        end do
        Ntot = ncell * nr * nz
        con1 = pi**third * (0.75 / matdens)**(2. * third)
        con2 = (0.75 / pi / matdens)**third
        
        close(fh)

        write(command,*) './directory.sh '//trim(datadir)
        CALL SYSTEM(command)
        open(unit=2,file='outputs/path.txt',action='read') !read path from file
        read(2,'(A)') path
        close(2)
        CALL SYSTEM('rm -rf outputs/path.txt')


        return
    end subroutine read_parameters

end module parameters
