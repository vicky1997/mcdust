!this module performs parallel sorting of a given array
! an array of swarms here specifically.
! routines taken from https://github.com/cphyc/Fortran-parallel-sort
module parallel_sort

    use omp_lib
    use mrgrnk
    use types

    implicit none

    public:: z_parallel_sort, rad_parallel_sort

    contains

 subroutine Z_parallel_sort (A, order)
    type(swarm), intent(in), dimension(:) :: A
    integer, intent(out), dimension(:) :: order 

    integer :: ilen, from, middle, ito, nthreads, thread, chunk, chunk2, i, iremainder, extraThread

    ilen      = size(A)
    nthreads = omp_get_max_threads()
    chunk    = ilen / nthreads

    iremainder = mod(ilen, nthreads)
    if (iremainder /= 0) then
        extraThread = 1
    else
        extraThread = 0
    endif
    !----------------------------------------
    ! Initialize order
    !----------------------------------------
    !$OMP parallel do shared(order)
    do i = 1, ilen
       order(i) = i
    end do
    !$OMP end parallel do

    !----------------------------------------
    ! Sort each chunk
    !----------------------------------------
    !$OMP parallel do default(shared) private(thread, from, ito) schedule(guided)
    do thread = 0, nthreads-1 + extraThread
       from = thread*chunk + 1
       ito  = min((thread + 1)*chunk, ilen)

       call z_mrgrnk(A(from:ito), order(from:ito))
       order(from:ito) = order(from:ito) + from - 1
    end do
    !$OMP end parallel do

    !----------------------------------------
    ! Merge pieces together
    !----------------------------------------
    i = 1
    chunk2 = chunk
    do while (chunk2 < size(A))

       !$OMP parallel do default(shared) private(thread, from, middle, ito)
       do thread = 0, ceiling(.5 * size(A) / chunk2)
          from   = thread*2     * chunk2 + 1
          middle = (thread*2 + 1) * chunk2
          ito     = (thread*2 + 2) * chunk2

          middle = min(middle, size(A))
          ito     = min(ito, size(A))
          if (from < ito) then
             call z_merge(A, order, from, middle, ito)
          end if
       end do
       !$OMP end parallel do

       chunk2 = chunk2 * 2
       i = i + 1
    end do
  end subroutine Z_parallel_sort

  !> Merge two parts of A, ordered by order from left to right
  !! around middle.
  subroutine Z_merge (A, order, left, middle, right)
    type(swarm), intent(in), dimension(:) :: A
    integer, intent(out), dimension(:) :: order    

    integer, intent(in) :: left, middle, right

    integer :: leftA, rightA, leftB, rightB
    integer :: iA, iB, i

    integer, dimension(left    :middle) :: orderA
    integer, dimension(middle+1:right ) :: orderB

    ! copy order
    orderA = order(left    :middle)
    orderB = order(middle+1:right)

    ! more explicit variables
    leftA  = left
    rightA = middle
    leftB  = middle+1
    rightB = right

    ! initialize iA, iB to their leftmost position
    iA = leftA
    iB = leftB

    i = leftA

    do while ((iA <= rightA) .and. (iB <= rightB))
       if (A(orderA(iA))%zdis <= A(orderB(iB))%zdis) then
          order(i) = orderA(iA)
          iA = iA + 1
       else
          order(i) = orderB(iB)
          iB = iB + 1
       end if

       i = i + 1
    end do

    ! either A or B still have elements, append them to the new order
    do while (iA <= rightA)
       order(i) = orderA(iA)
       iA = iA + 1

       i  = i + 1

    end do
    do while (iB <= rightB)
       order(i) = orderB(iB)
       iB = iB + 1

       i  = i + 1
    end do

  end subroutine z_merge
  subroutine rad_parallel_sort (A, order)
    !real(8), intent(in),  dimension(:) :: A
    !integer, intent(out), dimension(size(A)) :: order
    type(swarm), intent(in), dimension(:) :: A
    integer, intent(out), dimension(size(A)) :: order 

    integer :: ilen, from, middle, ito, nthreads, thread, chunk, chunk2, i, iremainder, extraThread
    
    !allocate(order(size(A)))
    ilen      = size(A)
    nthreads = omp_get_max_threads()
    chunk    = ilen / nthreads

    iremainder = mod(ilen, nthreads)
    if (iremainder /= 0) then
        extraThread = 1
    else
        extraThread = 0
    endif
    !----------------------------------------
    ! Initialize order
    !----------------------------------------
    !$OMP parallel do shared(order)
    do i = 1, ilen
       order(i) = i
    end do
    !$OMP end parallel do
    !----------------------------------------
    ! Sort each chunk
    !----------------------------------------
    !$OMP parallel do default(shared) private(thread, from, ito) schedule(guided)
    do thread = 0, nthreads-1 + extraThread
       from = thread*chunk + 1
       ito  = min((thread + 1)*chunk, ilen)

       call rad_mrgrnk(A(from:ito), order(from:ito))
       order(from:ito) = order(from:ito) + from - 1
    end do
    !$OMP end parallel do
    !----------------------------------------
    ! Merge pieces together
    !----------------------------------------
    i = 1
    chunk2 = chunk
    do while (chunk2 < size(A))

       !$OMP parallel do default(shared) private(thread, from, middle, ito)
       do thread = 0, ceiling(.5 * size(A) / chunk2)
          from   = thread*2     * chunk2 + 1
          middle = (thread*2 + 1) * chunk2
          ito     = (thread*2 + 2) * chunk2

          middle = min(middle, size(A))
          ito     = min(ito, size(A))
          if (from < ito) then
             call rad_merge(A, order, from, middle, ito)
          end if
       end do
       !$OMP end parallel do

       chunk2 = chunk2 * 2
       i = i + 1
    end do
  end subroutine rad_parallel_sort

  !> Merge two parts of A, ordered by order from left to right
  !! around middle.
  subroutine rad_merge (A, order, left, middle, right)
    type(swarm), intent(in), dimension(:) :: A
    integer, intent(out), dimension(:) :: order    

    integer, intent(in) :: left, middle, right

    integer :: leftA, rightA, leftB, rightB
    integer :: iA, iB, i

    integer, dimension(left    :middle) :: orderA
    integer, dimension(middle+1:right ) :: orderB
    ! copy order
    orderA = order(left    :middle)
    orderB = order(middle+1:right)

    ! more explicit variables
    leftA  = left
    rightA = middle
    leftB  = middle+1
    rightB = right

    ! initialize iA, iB to their leftmost position
    iA = leftA
    iB = leftB

    i = leftA

    do while ((iA <= rightA) .and. (iB <= rightB))
       if (A(orderA(iA))%rdis <= A(orderB(iB))%rdis) then
          order(i) = orderA(iA)
          iA = iA + 1
       else
          order(i) = orderB(iB)
          iB = iB + 1
       end if

       i = i + 1
    end do

    ! either A or B still have elements, append them to the new order
    do while (iA <= rightA)
       order(i) = orderA(iA)
       iA = iA + 1

       i  = i + 1

    end do
    do while (iB <= rightB)
       order(i) = orderB(iB)
       iB = iB + 1

       i  = i + 1
    end do

  end subroutine rad_merge
 
end module parallel_sort   