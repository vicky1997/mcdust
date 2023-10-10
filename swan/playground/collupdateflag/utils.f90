module utils
    implicit none

    public :: interp2d

    interface interp2d
        module procedure interp2d_dp
    end interface interp2d

    contains

    function interp2d_dp(x,y,array,x0,y0,bounds_error,fill_value) result(value)
        ! Bilinar interpolation of array = f(x,y) at (x0,y0)

        implicit none

        real(dp),intent(in) :: x(:),y(:),array(:,:),x0,y0

        logical,intent(in),optional :: bounds_error
        ! whether to raise an out of bounds error

        real(dp),intent(in),optional :: fill_value
        ! value for out of bounds if bounds_error is .false.

        real(dp) :: value,norm
        integer :: i1,i2,j1,j2

        logical :: bounds_error_tmp
        real(dp) :: fill_value_tmp

        if(present(bounds_error)) then
        bounds_error_tmp = bounds_error
        else
        bounds_error_tmp = .true.
        end if

        if(.not.bounds_error_tmp) then
        if(present(fill_value)) then
            fill_value_tmp = fill_value
        else
            fill_value_tmp = 0._dp
        end if
        end if

        if(size(x).ne.size(array,1)) stop "x does not match array"
        if(size(y).ne.size(array,2)) stop "y does not match array"

        i1 = locate(x,x0) ; i2 = i1 + 1
        j1 = locate(y,y0) ; j2 = j1 + 1

        if(i1==-1) then
        if(bounds_error_tmp) then
            write(0,'("ERROR: Interpolation out of bounds : ",ES11.4," in [",ES11.4,":",ES11.4,"]")') x0,x(1),x(size(x))
            stop
        else
            value = fill_value_tmp
            return
        end if
        end if

        if(j1==-1) then
        if(bounds_error_tmp) then
            write(0,'("ERROR: Interpolation out of bounds : ",ES11.4," in [",ES11.4,":",ES11.4,"]")') y0,y(1),y(size(y))
            stop
        else
            value = fill_value_tmp
            return
        end if
        end if

        norm = 1._dp / (x(i2) - x(i1)) / (y(j2)-y(j1))

        value =   array(i1,j1) * (x(i2)-x0) * (y(j2)-y0) * norm &
            &       + array(i2,j1) * (x0-x(i1)) * (y(j2)-y0) * norm &
            &       + array(i1,j2) * (x(i2)-x0) * (y0-y(j1)) * norm &
            &       + array(i2,j2) * (x0-x(i1)) * (y0-y(j1)) * norm

    end function interp2d_dp
end module utils