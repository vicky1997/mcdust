! This code performes 2D simulation of dust evolution in protoplanetary disk.
! The gas disk is treated 100% analytically (discstruct module)
! The dust is treated as representative particles (RPs) undergoing advection (advection module)
! as well as collisions performed with Monte Carlo algorithm (collisions module)
! To perform collisions the RPs are binned using an adaptive grid (grid module)
!
! citation: Drążkowska, Windmark & Dullemond (2013) A&A 556, A37
!
! Author: Joanna Drążkowska, Heidelberg University, ZAH, ITA
! Albert-Ueberle-Str. 2, 69120 Heidelberg, Germany
! mailto: asiadrazkowska@gmail.com
!
program main
   use constants
   use advection,    only: mc_advection
   use collisions,   only: mc_collisions
   use grid,         only: g, make_grid, deallocate_grid
   use initproblem,  only: init_swarms, init_random_seed, m0, mswarm, nord
   use discstruct,   only: alpha, cs, omegaK, sigmag, densg, Pg, vgas, gasmass, dlogPg
   use parameters,   only: read_parameters, Ntot, nz, nr, dtime, fout, tend, smallr, restart, restime, minrad0, &
                           maxrad0, matdens, r0, dtg, db_data
   use output,       only: write_output, read_restart
   use timestep,     only: time_step
   use types
   use hdf5
   use hdf5output, only: hdf5_file_write, hdf5_file_t

   implicit none

   ! the array of the representative particles (swarms) is declared here:
   type(swarm), dimension(:), allocatable, target                  :: swrm        ! list of all swarms in the simulation
   type(list_of_swarms), dimension(:,:), allocatable, target       :: bin         ! swarms binned into cells
   type(list_of_swarms), dimension(:), allocatable                 :: rbin        ! swarms binned into radial zones
   type(hdf5_file_t)                                               :: file
   real                       :: total
   real(kind=4), dimension(2) :: elapsed

   integer             :: i, j, iter
   real                :: time = 0.0          ! physical time
   real                :: timeofnextout = 0.0 ! time of next output
   real                :: resdt = 0.1*year               ! resulting physical time step
   integer             :: nout = 0            ! number of the next output
   real                :: totmass             ! total mass of dust beyond evaporation line
   character(len=100)  :: ctrl_file           ! parameter file
   real                :: mdust               ! mass of dust ! TODO: can be merged with totmass?
   integer, dimension(:,:), allocatable :: ncolls
   !real, dimension(:), allocatable  :: mgrid  ! mass grid

   ! random number generator initialization
   call init_random_seed

   ! reading parameters and initializing the simulation
   call get_command_argument(1, ctrl_file)
   write(*,*) 'Reading parameters...'
   call read_parameters(ctrl_file)
   write(*,*) '------------------------------------------------------------------'
   write(*,*) 'Initializing representative bodies...'
   if (restart) then
      write(*,*) ' Reading restart...'
      call read_restart(Ntot, swrm)
      !call hdf5_file_read(Ntot, swrm)
      write(*,*) '  restart read!'
      time = restime * year
      nout = nint(time/dtime)
      timeofnextout = nout * dtime
      mdust = dtg * gasmass(minrad0*AU,maxrad0*AU,0.0)  ! TODO
      mswarm = mdust / real(Ntot)
      m0 = 4. * third * pi * r0**3 * matdens
      nord = (log10(mswarm/m0))
   else
      call init_swarms(Ntot,swrm)
   endif
   write(*,*) 'succeed'

   write(*,*) 'Initial disk mass: ', gasmass(0.1*AU,maxrad0*AU,0.0)/Msun

   write(*,*) ' Making grid for the first time...'
   call make_grid(swrm, bin, rbin, nr, nz, smallr,totmass, ncolls)
   write(*,*) '  grid done'
   
   !allocate(ncolls(nr,nz))
   ncolls(:,:) = 1
   
   write(*,*) 'going into the main loop...'

   iter = 0
   ! ------------------- MAIN LOOP -------------------------------------------------------------------------------------
   do while (time < tend)

      ! producing output
      if (modulo(iter,fout) == 0 .or. time>=timeofnextout) then
         !call write_output(swrm, nout)
         call hdf5_file_write(file, swrm, time, 'create', nout)
         write(*,*) 'Time: ', time/year, 'produced output: ',nout
         open(23,file='timesout.dat',status='unknown',position='append')
         write(23,*) 'time: ', time/year, 'produced output: ',nout
         close(23)
         timeofnextout = time+dtime
         nout = nout + 1
      endif

      iter = iter + 1

      ! writing max mass value for each timestep for bug fixes
      if(db_data) then
         open(123,file='mmax.dat',position='append')
         write(123,*) time/year, maxval(swrm(:)%mass)
         close(123)
      endif


      ! determining the time step
      if (iter == 1) then 
         resdt = 1./omegaK(minval(swrm(:)%rdis))
      else
         call time_step(bin, ncolls, timeofnextout-time, resdt)
      
      end if
      write(*,*) ' Performing advection: timestep',resdt/year,'yrs'
      if (db_data) then
         open(23,file='timestep.dat',status='unknown',position='append')
         write(23,*) time/year,resdt/year
         close(23)
      endif

      ! performing advection
      call mc_advection(swrm, resdt, time)
      write(*,*) '  advection done'


      ! removing old grid and building new one
      call deallocate_grid
      if (allocated(bin))  deallocate(bin)
      if (allocated(rbin)) deallocate(rbin)
      if (allocated(ncolls)) deallocate(ncolls)
      write(*,*) '    Making grid...'
      call make_grid(swrm, bin, rbin, nr, nz, smallr,totmass, ncolls)
      write(*,*) '     grid done'

      !deallocate(ncolls)
      deallocate(swrm)
      !allocate(ncolls(nr,nz))

      ! performing collisions
      write(*,*) '   Performing collisions...'

      !$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(DYNAMIC)
      do i = 1, size(g%rce)
         do j = 1, size(g%zce,dim=2)
            if (.not.allocated(bin(i, j)%p)) cycle
            write(*,*) '    entering zone',i,j,'including ',size(bin(i, j)%p),' rbs','.....'
            call mc_collisions(i, j, bin, resdt, time, ncolls(i,j))
         enddo
      enddo
      !$OMP END PARALLEL DO

      write(*,*) '    collisions done!'

      ! build the new swrm list after collisions
      allocate(swrm(0))
      do i = 1, size(g%rce)
        do j = 1, size(g%zce,dim=2)
          swrm = [swrm, bin(i, j)%p(:)]
        enddo
      enddo

      time = time + resdt

   enddo
   ! ---- END OF THE MAIN LOOP -----------------------------------------------------------------------------------------
   !nout = nout+1
   write(*,*) 'time: ', time/year, 'produced output: ',nout
   open(23,file='timesout.dat',status='unknown',position='append')
   !call write_output(swrm, nout)
   write(23,*) 'time: ', time/year, 'produced output: ',nout
   close(23)

   call hdf5_file_write(file, swrm, time, 'create', nout)
   
   
   deallocate(bin)
   deallocate(rbin)

   write(*,*) '------------------------------------------------------------------'
   write(*,*) 'tend exceeded, finishing simulation...'

   ! this causes problems when used with intel compilers, so just remove it in case you want to use one
   total = etime(elapsed)
   write(*,*) 'Elapsed time [s]: ', total, ' user:', elapsed(1), ' system:', elapsed(2)

end
