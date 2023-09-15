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
   public      :: write_output, read_restart

   contains

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
         write(11,*) list(i)%idnr, list(i)%npar, list(i)%mass, list(i)%rdis/AU, list(i)%zdis/AU, list(i)%stnr
      enddo

      close(11)

      return
   end subroutine write_output

   ! reads restart
   subroutine read_restart(Ntot, swrm)
      use constants, only: AU
      implicit none
      type(swarm), dimension(:), allocatable :: swrm
      integer, intent(in)                    :: Ntot
      integer                                :: i

      if (.not.allocated(swrm)) allocate( swrm(Ntot) )

      open(1,file='restart',action='read')
      do i = 1, Ntot
         read(1,*) swrm(i)%idnr, swrm(i)%npar, swrm(i)%mass, swrm(i)%rdis, swrm(i)%zdis, swrm(i)%stnr
         swrm(i)%rdis = swrm(i)%rdis * AU
         swrm(i)%zdis = swrm(i)%zdis * AU
      enddo
      close(1)

      return
   end subroutine read_restart

end
