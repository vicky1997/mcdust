! this module gathers types declarations that can be used by other modules
module types
   implicit none
   public

   type :: swarm
      integer        :: idnr     ! id of the swarm
      real           :: npar     ! number of bodies
      real           :: mass     ! mass of one body
      real           :: rdis     ! cylindrical radius
      real           :: zdis     ! distance from the midplane
      real           :: velr     ! radial velocity (advection)
      real           :: velz     ! vertical velocity (advection)
      real           :: stnr     ! Stokes number of the particle
   end type swarm

   type :: list_of_swarms
      type(swarm), dimension(:), allocatable :: p
   end type list_of_swarms

   type :: grid_container
      real, dimension(:,:), allocatable :: zlo       ! lower boundary of bin in z
      real, dimension(:,:), allocatable :: zup       ! upper boundary of bin in z
      real, dimension(:,:), allocatable :: zce       ! central point of bin in z
      real, dimension(:), allocatable :: rlo         ! lower boundary of bin in r
      real, dimension(:), allocatable :: rup         ! upper boundary of bin in r
      real, dimension(:), allocatable :: rce         ! central point of bin in r
      real, dimension(:), allocatable :: dr          ! radial span of the bin
      real, dimension(:,:), allocatable :: dz        ! vertical span of the bin
      real, dimension(:,:), allocatable   :: vol     ! volume of the zone
      real                            :: totvol      ! total volume
   end type grid_container

end
