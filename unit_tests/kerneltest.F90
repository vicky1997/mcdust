program testkernel

    use types
    use hdf5output
    use hdf5
    use parameters, only: path, datadir, read_parameters

    implicit none
    !  ------- set parameters here --------------------------------------------------------------------------------------
   integer, parameter             :: ntot = 1000      ! total number of representative particles (swarms)
   real, parameter                :: dens = 1.       ! dust volume mass density
   real, parameter                :: mswrm = 10.e20  ! mass of the swarm
#ifdef FRAGTEST
   real, parameter                :: m0 = 1e-6
   real, parameter                :: gamma = 1e4
#else
   real, parameter                :: m0 = 1.         ! mass of the monomer (initial mass of the particles)
#endif
   real, parameter                :: tend = 20.      ! maximum time of the simulation
   real, parameter                :: dtout = 4.      ! time step for the output
   integer, parameter             :: nbins = 250     ! number of bins for output histograms
   integer, parameter             :: repeat = 10     ! number of times you want to repeat a simulation
   real, parameter                :: dmmax = 0.001
   ! -------------------------------------------------------------------------------------------------------------------

   type(swarm), dimension(ntot)      :: swarms             ! list of the swarms
   integer                        :: nri               ! index of representative partcile to undergo the next collision
   integer                        :: nrk               ! index of physical partcile to undergo the next collision
   real                           :: vol               ! volume of the cell
   real                           :: time = 0.0        ! time
   real                           :: dt             ! time step between the subsequent collisions
   real                           :: tout = dtout      ! time of the next output
   real, dimension(ntot,ntot)     :: colrates          ! matrix of the collision rates between the particles
   real, dimension(ntot)          :: colrates_rp       ! probability of choosing given representative particles
   real                           :: totrate           ! total collision rate
   real                           :: rand              ! random number
   real                           :: fin               ! temporary value
   character(len=11)              :: fname             ! output file name
   real, dimension(ntot)          :: colrates_old_col  ! for optimization reasons
   real, dimension(ntot)          :: accelncol ! collision optimisation
   integer                        :: i, j, cont_time, loops, kk
   character(len=100)               :: loops_str
   character(len=100)               :: command
   character(len=100)               :: ctrl_file
   character(len=200)               :: filepath
   type(hdf5_file_t)                :: file
   real, dimension(:), allocatable :: t_arr
   real, dimension(nbins)         :: m2fm              ! mass distribution function
   real, dimension(nbins+1)       :: mgrid             ! mass grid
   integer                        :: k = 1
   real                           :: ll, lll, nord

  datadir = 'kernels'
  call get_command_argument(1, ctrl_file)
  write(*,*) 'Reading parameters...'
  call read_parameters(ctrl_file)

#ifdef LINEARKERNEL
   !dtime = 20.
   allocate(t_arr(5))
   t_arr = (/4., 8., 12., 16., 20./)
#elif PRODUCTKERNEL
   !dtime = 0.9
   allocate(t_arr(3))
   t_arr = (/0.4, 0.7, 0.9/)
#else
   !dtime = 100000.
   allocate(t_arr(6))
   t_arr = (/1., 10., 100., 1000., 10000., 100000./)
#endif

accelncol(:) = 1.
! calculate the cell volume
vol = ntot * mswrm / dens

! initialize mass bins for histograms
mgrid(1) = m0
nord = log10(mswrm * ntot / m0)  ! how many orders of magnitude in mass should the histogram go through?
ll = real(nord) / real(nbins)
lll = ll
do i = 2, nbins+1
   mgrid(i) = m0 * 10.0 ** lll
   lll = ll * i
enddo

    do loops=1, repeat
        time = 0.
        cont_time = 1
        kk = 0
        swarms(:)%npar = mswrm / m0
        swarms(:)%mass = m0
        swarms(:)%rdis = 1.
        swarms(:)%zdis = 1.
        swarms(:)%velr = 0.
        swarms(:)%velz = 0.
        swarms(:)%coll_f = 0.
        swarms(:)%stnr = 0.

#ifdef LINEARKERNEL
        ! calculate the collision rates: this is the linear kernel
        do i = 1, ntot
            where((swarms(:)%mass / swarms(i)%mass) < dmmax)
               accelncol(:) = swarms(i)%mass * dmmax / swarms(:)%mass
            endwhere
            colrates(i,:) = swarms(:)%npar * 0.5 * (swarms(i)%mass + swarms(:)%mass) / vol / accelncol(:)
        enddo
#elif PRODUCTKERNEL
        do i = 1, ntot
            where((swarms(:)%mass / swarms(i)%mass) < dmmax)
               accelncol(:) = swarms(i)%mass * dmmax / swarms(:)%mass
            endwhere
            colrates(i,:) = swarms(:)%npar * (swarms(i)%mass * swarms(:)%mass) / vol / accelncol(:)
        enddo
#else        
        do i = 1, ntot
            where((swarms(:)%mass / swarms(i)%mass) < dmmax)
               accelncol(:) = swarms(i)%mass * dmmax / swarms(:)%mass
            endwhere
            colrates(i,:) = swarms(:)%npar * 1. / vol
        enddo
#endif
        colrates_rp = sum(colrates, dim=2)
        write(*,*) 'going inside main loop'
 !--------------- MAIN LOOP ------------------------------------------------------------------------------------------
        do
            


            ! calculate a total rate for collisions
            totrate = sum(colrates_rp)
    
            ! determine the time step for the next collision
            call random_number(rand)
            dt = - 1. / totrate * log(rand)  ! based on the Poisson distribution
    
            ! update time
            time = time + dt
    
            ! select representative particle nri to undergo the collision
            call random_number(rand)
            rand = rand * totrate
            j = 1
    
            fin = colrates_rp(1)
            do while (rand > fin)
            fin = fin + colrates_rp(j+1)
            j = j + 1
            enddo
            nri = j
    
            ! select physical particle nrk to undergo the collision
            call random_number(rand)
            rand = rand * colrates_rp(nri)
            j = 1
    
            fin = colrates(nri,1)
            do while (rand > fin)
            fin = fin + colrates(nri,j+1)
            j = j + 1
            enddo
            nrk = j
    
            ! perform the collision and update the representative particle nri (here: only sticking)
#if !defined(FRAGTEST)
            swarms(nri)%mass = swarms(nri)%mass + swarms(nrk)%mass
            swarms(nri)%npar = mswrm / swarms(nri)%mass

#else
            call random_number(rand)
            w = rand*(exp(-1*gamma*(swarms(nri)%mass + swarms(nrk)%mass))- exp(-1*gamma*m0)) &
                                                               + exp(-1*gamma*m0)
            swarns(nri)%mass = -1*log(w)/gamma
#endif
    
            ! remember the old column (optimization)
            colrates_old_col(:) = colrates(:,nri)
    
            ! update the collision rates matrix
            accelncol(:) = 1.
#ifdef LINEARKERNEL
            where((swarms(:)%mass / swarms(nri)%mass) < dmmax)
                  accelncol(:) = swarms(nri)%mass * dmmax / swarms(:)%mass
            endwhere
            colrates(nri,:) = swarms(:)%npar * 0.5 * (swarms(nri)%mass + swarms(:)%mass) / vol /  accelncol(:)
            where((swarms(nri)%mass / swarms(:)%mass) < dmmax)
                  accelncol(:) = swarms(:)%mass * dmmax / swarms(nri)%mass
            endwhere
            colrates(:,nri) = swarms(nri)%npar * 0.5 * (swarms(nri)%mass + swarms(:)%mass) / vol
            
#elif PRODUCTKERNEL
            where((swarms(:)%mass / swarms(nri)%mass) < dmmax)
               accelncol(:) = swarms(nri)%mass * dmmax / swarms(:)%mass
            endwhere
            colrates(nri,:) = swarms(:)%npar * (swarms(nri)%mass * swarms(:)%mass) / vol
            where((swarms(nri)%mass / swarms(:)%mass) < dmmax)
                  accelncol(:) = swarms(:)%mass * dmmax / swarms(nri)%mass
            endwhere
            colrates(:,nri) = swarms(nri)%npar * (swarms(nri)%mass * swarms(:)%mass) / vol

#else
            where((swarms(:)%mass / swarms(nri)%mass) < dmmax)
               accelncol(:) = swarms(nri)%mass * dmmax / swarms(:)%mass
            endwhere
            colrates(nri,:) = swarms(:)%npar * 1. / vol
            where((swarms(nri)%mass / swarms(:)%mass) < dmmax)
                  accelncol(:) = swarms(:)%mass * dmmax / swarms(nri)%mass
            endwhere
            colrates(:,nri) = swarms(nri)%npar * 1. / vol
#endif

            colrates_rp(:) = colrates_rp(:) + (colrates(:,nri) - colrates_old_col(:))           ! one could just write colrates_rp = sum(colrates, dim=2) here
            colrates_rp(nri) = sum(colrates(nri,:))                                             ! but this way is much faster
    
            if (t_arr(cont_time)<time) then
               write(loops_str, '(I0)') loops
               write(command, '(A,A,A)') './directory.sh ', trim(datadir), trim(loops_str)
               CALL SYSTEM(command)
               open(unit=2, file='outputs/path.txt', action='read')
               read(2,'(A)') path
               close(2)
               CALL SYSTEM('rm -rf outputs/path.txt')
               !call hdf5_file_write(file, swarms, time, cont_time-1, mswrm, time)
               !write(*,*) swarms(1)%npar
               !write(*,*) 'time: ',time,'producing output ',cont_time

               ! make histogram
               m2fm(:) = 0.0
               do i = 1, ntot
                  j = 1
                  do while (swarms(i)%mass >= mgrid(j))
                     j = j + 1
                  enddo
                  m2fm(j) = m2fm(j) + ((swarms(i)%npar) * swarms(i)%mass**2) / ((mgrid(j+1) - mgrid(j)) * mswrm * ntot)
               enddo
      
               write(fname,'(a4,i3.3,a4)') 'out-',cont_time,'.dat'
               filepath = trim(path)//"/"//trim(fname)
               open(cont_time,file=filepath)
                  do j = 1, nbins
                     write(cont_time, *) sqrt(mgrid(j+1)*mgrid(j)), m2fm(j)
                  enddo
               close(cont_time)
               write(fname,'(a4,i3.3,a4)') 'mgridout-',cont_time,'.dat'
               filepath = trim(path)//"/"//trim(fname)
               open(cont_time,file=filepath)
                  do j = 1,nbins+1
                     write(cont_time, *) mgrid(j)
                  enddo
               open(cont_time,file=trim(path)//trim('/timesout.dat'),status='unknown',position='append')
                   write(cont_time, *) time
               close(cont_time)

               cont_time = cont_time + 1
               if (cont_time>size(t_arr)) exit
            endif
    
        enddo
        !------------------- END OF THE MAIN LOOP ---------------------------------------------------------------------------
    
        write(*,*) 'tend exceeded, finishing simulation no.',  loops

    enddo
    deallocate(t_arr)

end program