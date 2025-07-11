! the module sets the grid to bin particles for performing collisions
! the number of particles per cell is at least equal to the one set in params.par
! if you want to refer to any value related to grid, you have to use this module
! see Drazkowska et al 2013, Sect 2.4, Fig 1
module grid
   use types
   use parallel_sort
   use constants ! all?
   implicit none

   private
   public           :: g, make_grid, deallocate_grid

   type(grid_container)          :: g

   contains

   ! find location of the walls in radial dimension
   subroutine make_grid_r(swrm, rbin, nrad, smallr,totmass)
      use constants,  only: smallv
      use initproblem, only: mswarm
      implicit none
      type(swarm), dimension(:), allocatable :: swrm
      type(list_of_swarms), dimension(:), allocatable :: rbin
      integer                             :: nrad    ! nr of zones in r
      integer, dimension(:), allocatable  :: nzone_arr ! nr of representative particles per r zone
      real                                :: minr, maxr ! innermost and outermost wall
      integer                             :: nrad0   ! initial number of particles per radial zone
      real, intent(in)                    :: smallr  ! "evaporation radius"
      real, intent(out)                   :: totmass ! total dust mass ouside of the "evaporation radius"
      integer                             :: i, k, nin
      integer, dimension(size(swrm))      :: r_order ! rank for sorting wrt radius
      real                                :: x       ! random variable
      integer                             :: xi      ! random index
      integer                             :: i_last 

      

      nrad0 = nrad

      ! making grid in radial direction: adjusting to the radial distances of particles
      minr = minval(swrm(:)%rdis)-smallv
      maxr = maxval(swrm(:)%rdis)+smallv
      if (minr < smallr)  minr = smallr - smallv
      if (maxr < smallr)  then
         write(*,*) 'END OF SIMULATION: all the dust is inside the evaporation radius...'
         stop
      endif

      ! counting particles outside smallr:
      k = count(swrm(:)%rdis > smallr)
      totmass = k * mswarm  ! mass of dust outside smallr

      nin = size(swrm) - k  ! nr of particles inside the small r

      ! how many zones in r?
      do while ( (real(k)/real(nrad)) < (real(size(swrm))/real(nrad0)) )
         nrad = nrad - 1
      enddo
      if (nrad < 1) nrad = 1
      write(*,*) '    -- Number of radial zones: ', nrad


      ! how many particles per zone?

      allocate(nzone_arr(nrad))
      nzone_arr(:) = floor(real(k)/real(nrad))

      ! redistribute the rest of particles (k-nzone*nrad) randomly
      do i = 1, k-nzone_arr(1)*nrad
         call random_number(x)
         xi = 1 + floor(x*nrad) ! find random index in (1, nrand) range
         nzone_arr(xi) = nzone_arr(xi) + 1
      enddo

      allocate ( g%rlo(nrad), g%rup(nrad), g%rce(nrad), g%dr(nrad) )
      ! sorting all the swarms by r distance
      if (nrad > 1) then
         if(size(swrm) .ge. 2097152 ) then
         call rad_parallel_sort(swrm, r_order)
         swrm = swrm(r_order)
         else 
           call shell_sort_r(swrm)
         endif
      endif

      allocate(rbin(nrad))
      
      ! save first index in swrm for every radial grid: useful in collision subroutine
      rbin(1)%first_idx = 1 + nin
      do i = 2, nrad
         rbin(i)%first_idx = rbin(i-1)%first_idx + nzone_arr(i-1)
      enddo

      ! establishing walls of r zones & putting particles there
      g%rlo(1) = minr
      if (nrad > 1) then
         g%rup(1) = 0.5*(swrm(nin + nzone_arr(1))%rdis + swrm(nin + nzone_arr(1) + 1)%rdis)
         g%rce(1) = sqrt(g%rlo(1) * g%rup(1))
         allocate(rbin(1)%p(nzone_arr(1)))
         rbin(1)%p = swrm(nin+1:nin + nzone_arr(1))
      endif
      do i = 2, nrad - 1
         g%rlo(i) = g%rup(i - 1)
         g%rup(i) = sqrt(swrm(nin + sum(nzone_arr(:i)))%rdis * swrm(nin + sum(nzone_arr(:i)) + 1)%rdis)
         g%rce(i) = sqrt(g%rlo(i) * g%rup(i))
         allocate(rbin(i)%p(nzone_arr(i)))
         rbin(i)%p = swrm(nin+sum(nzone_arr(:(i-1)))+1:nin + sum(nzone_arr(:i)))
      enddo
      if (nrad > 1) then
         g%rlo(nrad) = g%rup(nrad - 1)
         i_last = sum(nzone_arr(:nrad-1))
      else
         i_last = 0
      endif
      g%rup(nrad) = maxr
      g%rce(nrad) = sqrt(g%rlo(nrad) * g%rup(nrad))
      allocate(rbin(nrad)%p(nzone_arr(nrad)))
      rbin(nrad)%p = swrm(nin+i_last+1:)

      g%dr(:) = g%rup(:) - g%rlo(:)
      deallocate(nzone_arr)
      return
   end subroutine make_grid_r

   ! find location of the walls in vertical dimension for given radial zone
   subroutine make_grid_z(rbin, bin, k, nz)
      use constants,  only: smallv
      implicit none
      type(list_of_swarms), dimension(:),   allocatable               :: rbin
      type(list_of_swarms), dimension(:,:), allocatable, target       :: bin
      integer, intent(in)     :: k      ! nr of current r zone
      integer                 :: nz     ! nr of zones in z
      integer, dimension(:), allocatable  :: nzzone_arr ! nr of representative particles per z zone
      real                    :: minz   ! min z of swarms
      real                    :: maxz   ! max z of swarms
      integer                 :: i
      real                    :: x      ! random variable
      integer                 :: xi     ! random index
      integer                 :: i_last ! position of particle right before last grid nz

      ! adjusting the bottom and the top wall
      minz = minval(rbin(k)%p(:)%zdis)-smallv
      maxz = maxval(rbin(k)%p(:)%zdis)+smallv

      ! sorting the particles by z
      if (nz > 1) call sort_swrm_z(rbin(k)%p)

      ! how many particles per zone?
      allocate(nzzone_arr(nz))
      nzzone_arr(:) = floor(real(size(rbin(k)%p)) / real(nz))

      do i =1, k-nzzone_arr(1)*nz
         call random_number(x)
         xi = 1 + floor(x*nz)
         nzzone_arr(xi) = nzzone_arr(xi) + 1
      enddo

      ! setting the walls and putting particles into bins
      g%zlo(k,1) = minz
      if (nz > 1) then
         g%zup(k,1) = 0.5 * (rbin(k)%p(nzzone_arr(1))%zdis + rbin(k)%p(nzzone_arr(1) + 1)%zdis)
         g%zce(k,1) = 0.5 * (g%zlo(k,1) + g%zup(k,1))
         allocate(bin(k,1)%p(nzzone_arr(1)))
         bin(k,1)%p = rbin(k)%p(:nzzone_arr(1))
      endif
      do i = 2, nz - 1
         g%zlo(k,i) = g%zup(k,i - 1)
         g%zup(k,i) = 0.5 * (rbin(k)%p(sum(nzzone_arr(:i)))%zdis + rbin(k)%p(sum(nzzone_arr(:i)) + 1)%zdis)
         g%zce(k,i) = 0.5 * (g%zlo(k,i) + g%zup(k,i))
         allocate(bin(k,i)%p(nzzone_arr(i)))
         bin(k,i)%p = rbin(k)%p(sum(nzzone_arr(:(i-1)))+1:sum(nzzone_arr(:i)))
      enddo
      if (nz > 1) then
         g%zlo(k,nz) = g%zup(k,nz - 1)
         i_last = sum(nzzone_arr(:nz-1))
      else
         i_last = 0
      endif
      g%zup(k,nz) = maxz
      g%zce(k,nz) = 0.5 * (g%zlo(k,nz) + g%zup(k,nz))
      allocate(bin(k,nz)%p(nzzone_arr(nz)))
      bin(k,nz)%p = rbin(k)%p(i_last+1:size(rbin(k)%p))

      g%dz(k,:) = g%zup(k,:) - g%zlo(k,:)

      deallocate(nzzone_arr)
      return
   end subroutine make_grid_z

   ! makes 2D grid using cylindrical coordinates: r and z
   subroutine make_grid(swrm, bin, rbin, nr, nz, smallr,totmass,ncolls)
      use constants,  only: pi
      implicit none
      type(swarm), dimension(:), allocatable :: swrm
      type(list_of_swarms), dimension(:,:), allocatable, target       :: bin
      type(list_of_swarms), dimension(:),   allocatable               :: rbin
      integer, dimension(:,:), allocatable                            :: ncolls
      integer                 :: nz      ! nr of zones in z
      integer                 :: nr      ! nominal nr of radial zones
      real, intent(in)        :: smallr  ! "evaporation radius" - inner edge of the simulation
      real, intent(out)       :: totmass ! total dust mass outside of evaporation radius
      integer                 :: nrad    ! present nr of radial zones
      integer                 :: i, k

      ! we start with the number of radial zones read from params.par
      nrad = nr
      ! making the grid in r
      call make_grid_r(swrm, rbin, nrad, smallr,totmass)
      ! now we know how many r zones do we actually have, so we can allocate arrays
      allocate( bin(nrad,nz) )
      allocate( ncolls(nrad,nz) )
      allocate( g%zlo(nrad,nz), g%zup(nrad,nz), g%zce(nrad,nz), g%dz(nrad,nz) )
      ! making the grid in z: separately for every radial zone
      !$OMP PARALLEL DO SCHEDULE(DYNAMIC)
      do k = 1, nrad
         call make_grid_z(rbin, bin, k, nz)
      enddo
      !$OMP END PARALLEL DO
      ! calculating the volumes of cells and total volume of simulation
      allocate(g%vol(nrad,nz))
      do i = 1, nrad
         g%vol(i,:) = pi * (g%rup(i)**2 - g%rlo(i)**2) * g%dz(i,:)
      enddo
      g%totvol = sum(g%vol(:,:))
      return
   end subroutine make_grid

   subroutine deallocate_grid
      implicit none

      if (allocated(g%rlo))      deallocate(g%rlo)
      if (allocated(g%rce))      deallocate(g%rce)
      if (allocated(g%rup))      deallocate(g%rup)
      if (allocated(g%dr))       deallocate(g%dr)
      if (allocated(g%zlo))      deallocate(g%zlo)
      if (allocated(g%zup))      deallocate(g%zup)
      if (allocated(g%zce))      deallocate(g%zce)
      if (allocated(g%dz))       deallocate(g%dz)
      if (allocated(g%vol))      deallocate(g%vol)

      return
   end subroutine deallocate_grid

   ! sorting particles by radial distance using shell sort algorithm
   subroutine shell_sort_r(swrm)
      implicit none
      type(swarm), dimension(:), allocatable :: swrm
      integer                                :: i, j, increment
      type(swarm)                            :: tempswarm

      increment = size(swrm) / 2
      do while (increment > 0)
         do i = increment+1, size(swrm)
            j = i
            tempswarm = swrm(j)
            do 
               if (j < increment+1) exit
               if (swrm(j-increment)%rdis <= tempswarm%rdis) exit 
               swrm(j) = swrm(j-increment)
               j = j - increment
            enddo
            swrm(j) = tempswarm
         enddo
         if (increment == 2) then
            increment = 1
         else
            increment = increment * 5 / 11
         endif
      enddo

      return
   end subroutine shell_sort_r

   ! sorting swrm table by z using insertion sort algorithm
   subroutine sort_swrm_z(lista)
      implicit none
      type(swarm), dimension(:), allocatable :: lista
      real                                   :: temp
      type(swarm)                            :: tempswarm
      integer                                :: i, j

      ! in this sort method the first element has to be already the smallest
      i = minloc(lista(:)%zdis, dim=1)
      tempswarm = lista(i)
      lista(i) = lista(1)
      lista(1) = tempswarm

      do i = 2, size(lista)
         j = i - 1
         temp = lista(i)%zdis
         tempswarm = lista(i)
         do while (j>=1 .and. lista(j)%zdis>temp)
            lista(j+1) = lista(j)
            j = j - 1
         end do
         lista(j+1) = tempswarm
      enddo

      return
   end subroutine sort_swrm_z

end
