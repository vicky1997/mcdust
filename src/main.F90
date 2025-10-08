! Authors: Joanna Drążkowska, Vignesh Vaikundaraman, Nerea Gurrutxaga
! Max Planck Institute for Solar System Research, Göttingen, Germany
!
! This code performs a 2D simulation of dust evolution in a protoplanetary disk.
! The gas disk is an input and is not evolved by mcdust (discstruct module)
! The dust is treated as representative particles (RPs) undergoing advection (advection module)
! as well as collisions performed with Monte Carlo algorithm (collisions module)
! To perform collisions the RPs are binned using an adaptive grid (grid module)
!
! citations: Drążkowska, Windmark & Dullemond (2013) A&A 556, A37, Vaikundaraman et al. (2025, submitted to JOSS)
!
! This code is licensed under the GNU General Public License v3
!
program main
   use constants
   use advection,    only: mc_advection, update_St
   use collisions,   only: mc_collisions
   use grid,         only: g, make_grid, deallocate_grid
   use initproblem,  only: init_swarms, init_random_seed, m0, mswarm, nord
   use discstruct,   only: cs, omegaK, gasmass
   use parameters,   only: read_parameters, Ntot, nz, nr, dtime, fout, tend, smallr, restart, &
                           maxrad0, matdens, r0, db_data, path
#ifdef LOGTIME
   use parameters,   only: read_parameters, Ntot, nz, nr, dtime, fout, tend, smallr, restart, &
                           maxrad0, matdens, r0, db_data, path, tstart, ntimeout
#endif
   use timestep,     only: time_step
   use types
   use hdf5
   use hdf5output, only: hdf5_file_write, hdf5_file_t, hdf5_file_read

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
#ifdef LOGTIME
   real, dimension(:),allocatable  :: timearray
   real                       :: logdt
#endif
   real                :: timeofnextout = 0.0 ! time of next output
   real                :: resdt = 0.1*year               ! resulting physical time step
   integer             :: nout = 0            ! number of the next output
   real                :: totmass             ! total mass of dust beyond evaporation line
   character(len=100)  :: ctrl_file           ! parameter file
   real                :: mdust               ! mass of dust ! TODO: can be merged with totmass?
   integer, dimension(:,:), allocatable :: ncolls

   ! random number generator initialization
   call init_random_seed

   ! reading parameters and initializing the simulation
   call get_command_argument(1, ctrl_file)
   write(*,*) 'mcdust v1.0'
   write(*,*) '------------------------------------------------------------------'
   write(*,*) 'Reading parameters...'
   call read_parameters(ctrl_file)
   write(*,*) '------------------------------------------------------------------'
   write(*,*) 'Initializing representative bodies...'
#ifdef RESTART
      write(*,*) ' Reading restart...'
      call hdf5_file_read(Ntot, swrm, nout, mswarm, time, resdt)
      write(*,*) time, nout, mswarm
      write(*,*) '  restart read!'
      timeofnextout = time + dtime
      mdust = mswarm * real(Ntot)
      m0 = 4. * third * pi * r0**3 * matdens
      nord = (log10(mswarm/m0))
#else
      call init_swarms(Ntot,swrm)
#endif
#ifdef LOGTIME
      allocate(timearray(ntimeout+1))
      logdt = log(tend/tstart)/ntimeout
      do i = 1,ntimeout+1
         timearray(i) = tstart*exp((i-1)*logdt)
      enddo
      open(33,file=trim(path)//trim('/timescheck.dat'),status='unknown',position='append')
      do i = 1,size(timearray)
         write(33,*) timearray(i)/year
      enddo
      close(33)
      timeofnextout = timearray(nout+1)
      stop
#endif
   write(*,*) 'succeed'

   write(*,*) 'Initial disk mass: ', gasmass(0.1*AU,maxrad0*AU,0.0)/Msun

   write(*,*) ' Making grid for the first time...'
   call make_grid(swrm, bin, rbin, nr, nz, smallr,totmass, ncolls)
   write(*,*) '  grid done'
   
   ncolls(:,:) = 1
   
   write(*,*) 'going into the main loop...'

   iter = 0
   ! ------------------- MAIN LOOP -------------------------------------------------------------------------------------
   do while (time < tend)

      ! determining the time step
      if (iter == 0) then 
         resdt = 1./omegaK(minval(swrm(:)%rdis))
      else
         call time_step(bin, ncolls, timeofnextout-time, resdt)
      end if

#ifdef LOGTIME
      if(time>=timeofnextout) then
         call update_St(swrm, time) ! update Stokes number; it will later be update during advection too
         call hdf5_file_write(file, swrm, time, nout, mswarm,resdt)
         write(*,*) 'Time: ', time/year, 'produced output: ',nout
         open(23,file=trim(path)//trim('/timesout.dat'),status='unknown',position='append')
         write(23,*) 'time: ', time/year, 'produced output: ',nout
         close(23)
         nout = nout + 1
         timeofnextout = timearray(nout+1)
      endif
#else

      ! producing output
      if (modulo(iter,fout) == 0 .or. time>=timeofnextout) then
         call update_St(swrm, time) ! update Stokes number; it will later be update during advection too
         call hdf5_file_write(file, swrm, time, nout, mswarm,resdt)
         write(*,*) 'Time: ', time/year, 'produced output: ',nout
         open(23,file=trim(path)//trim('/timesout.dat'),status='unknown',position='append')
         write(23,*) 'time: ', time/year, 'produced output: ',nout
         close(23)
         timeofnextout = time+dtime
         nout = nout + 1
      endif
#endif
      iter = iter + 1

      ! writing max mass value for each timestep for bug fixes
#ifdef AUXDATA
         open(123,file=trim(path)//trim('/mmax.dat'),position='append')
         write(123,*) time/year, maxval(swrm(:)%mass)
         close(123)
#endif

#ifdef TRANSPORT
      write(*,*) ' Performing advection: timestep',resdt/year,'yrs'
#ifdef AUXDATA
         open(23,file=trim(path)//trim('/timestep.dat'),status='unknown',position='append')
         write(23,*) time/year,resdt/year
         close(23)
#endif

      ! performing advection
      call mc_advection(swrm, resdt, time)
      write(*,*) '  advection done'
#endif

#ifdef COLLISIONS
      ! removing old grid and building new one
      call deallocate_grid
      if (allocated(bin))  deallocate(bin)
      if (allocated(rbin)) deallocate(rbin)
      if (allocated(ncolls)) deallocate(ncolls)
      write(*,*) '    Making grid...'
      call make_grid(swrm, bin, rbin, nr, nz, smallr,totmass, ncolls)
      write(*,*) '     grid done'

      ! performing collisions
      write(*,*) '   Performing collisions...'

      !$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(DYNAMIC)
      do i = 1, size(g%rce)
         do j = 1, size(g%zce,dim=2)
            if (.not.allocated(bin(i, j)%p)) cycle
            write(*,*) '    entering zone',i,j,'including ',size(bin(i, j)%p),' rbs','.....'
            call mc_collisions(i, j, bin, swrm, resdt, time, ncolls(i,j), rbin(i)%first_idx)
         enddo
      enddo
      !$OMP END PARALLEL DO

      write(*,*) '    collisions done!'
#endif
      time = time + resdt

   enddo
   ! ---- END OF THE MAIN LOOP -----------------------------------------------------------------------------------------

   write(*,*) 'time: ', time/year, 'produced output: ',nout
   open(23,file=trim(path)//trim('/timesout.dat'),status='unknown',position='append')
   write(23,*) 'time: ', time/year, 'produced output: ',nout
   close(23)
   call update_St(swrm, time)
   call hdf5_file_write(file, swrm, time, nout, mswarm, resdt)
   
   
   deallocate(bin)
   deallocate(rbin)
#ifdef LOGTIME
   deallocate(timearray)
#endif
   write(*,*) '------------------------------------------------------------------'
   write(*,*) 'tend exceeded, finishing simulation...'

   ! this causes problems when used with intel compilers, so just remove it in case you want to use one
   total = etime(elapsed)
   write(*,*) 'Elapsed time [s]: ', total, ' user:', elapsed(1), ' system:', elapsed(2)

end
