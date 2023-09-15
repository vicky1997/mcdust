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
   use parameters,   only: read_parameters, Ntot, nz, nr, fout, tend, smallr, restart, restime, minrad0, &
                           maxrad0, matdens, r0, eta, db_data
   use output,       only: write_output, read_restart
   use timestep,     only: time_step
   use types
   use hdf5
   use hdf5output!, only: hdf5_file_write, hdf5_file_t, hdf5_file_close

   implicit none

   ! the array of the representative particles (swarms) is declared here:
   type(swarm), dimension(:), allocatable, target                  :: swrm        ! list of all swarms in the simulation
   type(list_of_swarms), dimension(:,:), allocatable, target       :: bin         ! swarms binned into cells
   type(list_of_swarms), dimension(:), allocatable                 :: rbin        ! swarms binned into radial zones
   type(hdf5_file_t)                                               :: file
   real                       :: total
   real(kind=4), dimension(2) :: elapsed

   integer             :: i, j, iter, error
   real                :: time = 0.0          ! physical time
   real                :: timeofnextout = 0.0 ! time of next output
   real                :: resdt               ! resulting physical time step
   integer             :: nout = 0            ! number of the nect output
   real                :: totmass             ! total mass of dust beyond evaporation line
   character(len=100)  :: ctrl_file           ! parameter file
   real                :: mdust               ! mass of dust ! TODO: can be merged with totmass?
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
      !call read_restart(Ntot, swrm)
      call hdf5_file_read(Ntot, swrm)
      write(*,*) '  restart read!'
      time = restime * year
      !nout = 49 ! TODO
      timeofnextout = restime * year
      mdust = eta * gasmass(minrad0*AU,maxrad0*AU,0.0)  ! TODO
      mswarm = mdust / real(Ntot)
      m0 = 4. * third * pi * r0**3 * matdens
      nord = (log10(mswarm/m0))
   else
      call init_swarms(Ntot,swrm)
   endif
   write(*,*) 'succeed'

   write(*,*) 'Initial disk mass: ', gasmass(0.1*AU,maxrad0*AU,0.0)/Msun

   write(*,*) ' Making grid for the first time...'
   call make_grid(swrm, bin, rbin, nr, nz, smallr,totmass)
   write(*,*) '  grid done'

   write(*,*) 'going into the main loop...'

   iter = 0
   !call h5open_f(error)
   call hdf5_file_write (file, swrm, 0., 'create',nout)
   !call hdf5_file_close(file)
   !call h5close_f(error)
   !stop
   ! ------------------- MAIN LOOP -------------------------------------------------------------------------------------
   do while (time < tend)

      ! determining the time step
      call time_step(swrm, resdt)
      !write(*,*) ' Performing advection: timestep',resdt/year,'yrs'
      if (db_data) then
         open(23,file='timestep.dat',status='unknown',position='append')
         write(23,*) time/year,resdt/year
         close(23)
      endif

      ! performing advection
      call mc_advection(swrm, resdt, time)
      write(*,*) '  advection done'

      ! writing max mass value for each timestep for bug fixes
      if(db_data) then
         open(123,file='mmax.dat',position='append')
         write(123,*) time/year, maxval(swrm(:)%mass)
         close(123)
      endif

      ! removing old grid and building new one
      call deallocate_grid
      if (allocated(bin))  deallocate(bin)
      if (allocated(rbin)) deallocate(rbin)
      write(*,*) '    Making grid...'
      call make_grid(swrm, bin, rbin, nr, nz, smallr,totmass)
      write(*,*) '     grid done'

      
      ! performing collisions
      write(*,*) '   Performing collisions...'

      !$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(DYNAMIC)
      do i = 1, size(g%rce)
         do j = 1, size(g%zce,dim=2)
            if (.not.allocated(bin(i, j)%p)) cycle
            !write(*,*) '    entering zone',i,j,'including ',size(bin(i, j)%p),' rbs','.....'
            call mc_collisions(i, j, bin, swrm, resdt, time)
         enddo
      enddo
      !$OMP END PARALLEL DO

      write(*,*) '    collisions done!'

      ! producing output
      if (modulo(iter,fout) == 0) then
         !call write_output(swrm, nout)
         nout = nout + 1
         !call h5open_f(error)
         call hdf5_file_write(file, swrm, time, 'create', nout)
         !call hdf5_file_close(file)
         !call h5close_f(error)
         write(*,*) 'Time: ', time/year, 'produced output: ',nout
         !open(23,file='timesout.dat',status='unknown',position='append')
         !write(23,*) 'time: ', time/year, 'produced output: ',nout
         !close(23)
         !nout = nout + 1
      endif

      iter = iter + 1


      time = time + resdt

   enddo
   ! ---- END OF THE MAIN LOOP -----------------------------------------------------------------------------------------
   nout = nout+1
   write(*,*) 'time: ', time/year, 'produced output: ',nout
   !open(23,file='timesout.dat',status='unknown',position='append')
   !call write_output(swrm, nout)
   !write(23,*) 'time: ', time/year, 'produced output: ',nout
   !close(23)
   
   !call h5open_f(error)
   call hdf5_file_write(file, swrm, time, 'create', nout)
   !call hdf5_file_close(file)
   !call h5close_f(error)

   
   deallocate(bin)
   deallocate(rbin)

   write(*,*) '------------------------------------------------------------------'
   write(*,*) 'tend exceeded, finishing simulation...'

   ! this causes problems when used with intel compilers, so just remove it in case you want to use one
   total = etime(elapsed)
   write(*,*) 'Elapsed time [s]: ', total, ' user:', elapsed(1), ' system:', elapsed(2)

end
