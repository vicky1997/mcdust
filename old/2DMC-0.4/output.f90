! module gathers output producing routines
module output

   use types,       only: swarm, list_of_swarms
   use parameters,  only: con2, matdens, vfrag, nbins, smallr
   use constants,   only: AU, pi, third, mH2, AH2
   use discstruct,  only: densg, sigmag, cs, omegaK, Pg, alpha
   use grid,        only: g
   use initproblem, only: m0, mswarm, nord, nmonom0

   implicit none

   private
   public      :: make_histogram, write_output, read_restart, make_mass_grid

   contains

   ! make the mass grid for histrograms
   subroutine make_mass_grid(mgrid)
      implicit none
      real, dimension(nbins+1)               :: mgrid
      real                                   :: mmin
      real                                   :: ll, lll
      integer                                :: i

      mmin = m0
      mgrid(1) = mmin
      ll = real(nord) / real(nbins)
      lll = ll
      do i = 2, nbins+1
         mgrid(i) = mmin * 10.0 ** lll
         lll = ll * i
      enddo

      return
   end subroutine make_mass_grid

   ! this routine makes f(m) for given list of swarms
   subroutine make_histogram(list,nout,mgrid)
      implicit none
      type(swarm), dimension(:), allocatable :: list
      integer, intent(in)                    :: nout
      character(len=16)                      :: filename
      integer                                :: i, j
      integer                                :: ncur
      real, dimension(nbins+1)               :: mgrid
      real, dimension(nbins)                 :: mf2m

      ncur = size(list)

      mf2m = 0.0

      do i = 1, nbins-1
         do j = 1, ncur
            if ((list(j)%mass.ge.mgrid(i)).and.(list(j)%mass.lt.(mgrid(i+1))).and.(list(j)%rdis.ge.smallr)) &
            mf2m(i)=mf2m(i) + ((list(j)%npar) * list(j)%mass**2) / ((mgrid(i+1) - mgrid(i)) * m0)
        enddo
      enddo

      i=nbins
      do j = 1, ncur
         if (list(j)%mass.ge.(mgrid(i)).and.(list(j)%rdis.ge.smallr)) mf2m(nbins) = mf2m(nbins) + &
         ((list(j)%npar) * list(j)%mass**2) / ((mgrid(i) - mgrid(i-1)) * m0)
      enddo

      mf2m(:) = mf2m(:) / nmonom0

      write(filename,'(a7,i5.5,a4)') 'm2fbin-',nout,'.dat'
      open(88,file=filename,position='append')
      do i = 1, nbins
         write(88,*) 0.5 * (mgrid(i+1)+mgrid(i))/m0, con2*(0.5*(mgrid(i+1)+mgrid(i)))**third, mf2m(i)
      enddo
      close(88)

      return
   end subroutine make_histogram

   ! this routine writes down all the swarms' properties
   subroutine write_output(list,outnr)
      implicit none
      type(swarm), dimension(:), allocatable :: list
      integer, intent(in)                    :: outnr
      integer                                :: i
      character(len=16)                      :: filename

      write(filename,'(a7,i5.5,a4)') 'swarms-',outnr,'.out'
      open(11,file=filename,status='unknown')

      do i = 1, size(list)
         write(11,*) list(i)%idnr, list(i)%npar, list(i)%mass, list(i)%rdis/AU, list(i)%zdis/AU, &
                    list(i)%stnr, list(i)%velr, list(i)%velz
      enddo

      close(11)

      return
   end subroutine write_output

   ! reads restart
   subroutine read_restart(Ntot, swrm)
      use constants, only: third, pi, AU
      implicit none
      type(swarm), dimension(:), allocatable :: swrm
      integer, intent(in)                    :: Ntot
      integer                                :: i

      if (.not.allocated(swrm)) allocate( swrm(Ntot) )

      open(1,file='restart',action='read')
      do i = 1, Ntot
         read(1,*) swrm(i)%idnr, swrm(i)%npar, swrm(i)%mass, swrm(i)%rdis, swrm(i)%zdis, swrm(i)%stnr, &
         swrm(i)%velr, swrm(i)%velz
         swrm(i)%rdis = swrm(i)%rdis * AU
         swrm(i)%zdis = swrm(i)%zdis * AU
      enddo
      close(1)

      return
   end subroutine read_restart

end
