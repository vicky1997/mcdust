module utils
! code from astrofrog/fortranlib
    use constants
    implicit none

    public :: interp2d, interp1d

    !interface interp2d
    !    module procedure interp2d
    !end interface interp2d

    contains

    function interp2d(x,y,array,x0,y0,bounds_error,fill_value) result(value)
        ! Bilinar interpolation of array = f(x,y) at (x0,y0)

        implicit none
        real,intent(in) :: x(:),y(:),array(:,:),x0,y0

        logical,intent(in),optional :: bounds_error
        ! whether to raise an out of bounds error

        real,intent(in),optional :: fill_value
        ! value for out of bounds if bounds_error is .false.

        real :: value,norm
        integer :: i1,i2,j1,j2

        logical :: bounds_error_tmp
        real :: fill_value_tmp

        if(present(bounds_error)) then
        bounds_error_tmp = bounds_error
        else
        bounds_error_tmp = .true.
        end if

        if(.not.bounds_error_tmp) then
        if(present(fill_value)) then
            fill_value_tmp = fill_value
        else
            fill_value_tmp = 0.
        end if
        end if

        !write(*,*) x0/AU
        !write(*,*) y0/AU

        if(size(x).ne.size(array,1)) stop "x does not match array"
        if(size(y).ne.size(array,2)) stop "y does not match array"

        i1 = locate(x,x0) ; i2 = i1 + 1
        j1 = locate(y,y0) ; j2 = j1 + 1

        if(i1==-1) then
        if(bounds_error_tmp) then
            write(0,'("ERROR: Interpolation out of bounds : ",ES11.4," in [",ES11.4,":",ES11.4,"]")') x0,x(1),x(size(x))
        else
            value = fill_value_tmp
            return
        end if
        end if

        if(j1==-1) then
        if(bounds_error_tmp) then
            write(0,'("ERROR: Interpolation out of bounds : ",ES11.4," in [",ES11.4,":",ES11.4,"]")') y0,y(1),y(size(y))
            !stop
        else
            value = fill_value_tmp
            return
        end if
        end if

        norm = 1. / (x(i2) - x(i1)) / (y(j2)-y(j1))

        value =   array(i1,j1) * (x(i2)-x0) * (y(j2)-y0) * norm &
            &       + array(i2,j1) * (x0-x(i1)) * (y(j2)-y0) * norm &
            &       + array(i1,j2) * (x(i2)-x0) * (y0-y(j1)) * norm &
            &       + array(i2,j2) * (x0-x(i1)) * (y0-y(j1)) * norm

    end function interp2d
    
    
    
	real function interp1d(x,y,xval,bounds_error,fill_value) result(yval)
		implicit none
		real,dimension(:),intent(in) :: x,y
		real,intent(in) :: xval
		logical,intent(in),optional :: bounds_error
		real,intent(in),optional :: fill_value
		yval = interp1d_general(x,y,xval,interp1d_single,bounds_error,fill_value)
	end function interp1d
    
    
    
	real function interp1d_single(x1,y1,x2,y2,xval) result(yval)
		real,intent(in) :: x1,y1,x2,y2,xval
		real :: frac
		frac = ( xval - x1 ) / ( x2 - x1 )
		yval = y1 + frac * ( y2 - y1 )
	end function interp1d_single
    
    


	real function interp1d_general(x,y,xval,f_single,bounds_error,fill_value) result(yval)

		implicit none

		integer :: n
		! the size of the array

		real,dimension(:),intent(in) :: x,y
		! the x and y arrays

		real,intent(in) :: xval
		! the value at which to interpolate y

		interface
		   real function f_single(x1,y1,x2,y2,xval)
		     implicit none
		     real,intent(in) :: x1,y1,x2,y2,xval
		   end function f_single
		end interface

		logical,intent(in),optional :: bounds_error
		! whether to raise an out of bounds error

		real,intent(in),optional :: fill_value
		! value for out of bounds if bounds_error is .false.

		integer :: ipos
		! position of x value in x array

		logical :: bounds_error_tmp
		real :: fill_value_tmp

		if(present(bounds_error)) then
		   bounds_error_tmp = bounds_error
		else
		   bounds_error_tmp = .true.
		end if

		if(.not.bounds_error_tmp) then
		   if(present(fill_value)) then
		      fill_value_tmp = fill_value
		   else
		      fill_value_tmp = 0.
		   end if
		end if

		n = size(x)

		ipos = locate(x,xval)

		! --- First some error checking --- !

		if(ipos == -1) then
		   if(bounds_error_tmp) then
		      write(0,'("ERROR: Interpolation out of bounds : ",ES11.4," in [",ES11.4,":",ES11.4,"]")') xval,x(1),x(n)
		      stop
		   else
		      yval = fill_value_tmp
		      return
		   end if
		end if

		if( ipos < n .and. ipos > 0) then
		   yval = f_single(x(ipos), y(ipos), x(ipos+1), y(ipos+1), xval)
		else if(ipos == n) then
		   yval = y(n)
		else if(ipos == 0) then
		   yval = y(1)
		else
		   write(0,'("ERROR: Unexpected value of ipos : ",I0)') ipos
		   stop
		end if

	end function interp1d_general
  
  
  
    
    integer function locate(xx,x)
    ! Locate a value in a sorted array

    implicit none
    real, dimension(:), intent(in) :: xx
    real, intent(in) :: x
    integer :: n,jl,jm,ju
    logical :: ascnd
    n=size(xx)
    ascnd = (xx(n) >= xx(1))
    jl=0
    ju=n+1
    do
       if (ju-jl <= 1) exit
       jm=(ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl=jm
       else
          ju=jm
       end if
    end do

    if (x == xx(1)) then
       locate = 1
    else if (x == xx(n)) then
       locate = n-1
    else if(ascnd.and. (x > xx(n) .or. x < xx(1))) then
       locate = -1
    else if(.not.ascnd.and. (x < xx(n) .or. x > xx(1))) then
       locate = -1
    else
       locate = jl
    end if

  end function locate
  

end module utils
