! put your gas disk model here
! x is distance from the star (in cm)
module discstruct
    use constants, only: AU, kB, mH2, Msun, year
    use parameters, only: alphaMRI, sigmag0, temperature, eta
    implicit none

    private
    public   :: alpha, cs, omegaK, sigmag, densg, Pg, vgas, ddensgdz, ddensgdr, gasmass, dlogPg

    interface alpha
        module procedure alpha_ds
        module procedure alpha_d1
    end interface alpha

    interface cs
        module procedure cs_ds
        module procedure cs_d1
    end interface cs

    interface omegaK
        module procedure omegaK_ds
        module procedure omegaK_d1
    end interface omegaK

    interface sigmag
        module procedure sigmag_ds
        module procedure sigmag_d1
    end interface sigmag

    interface densg
        module procedure densg_ds
        module procedure densg_d1
    end interface densg

    interface Pg
        module procedure Pg_ds
        module procedure Pg_d1
    end interface Pg

    interface vgas
        module procedure vgas_ds
        module procedure vgas_d1
    end interface vgas

    interface ddensgdz
        module procedure ddensgdz_ds
        module procedure ddensgdz_d1
    end interface ddensgdz

    interface ddensgdr
        module procedure ddensgdr_ds
        module procedure ddensgdr_d1
    end interface ddensgdr

    interface gasmass
        module procedure gasmass_ds
        module procedure gasmass_d1
    end interface gasmass

    interface dlogPg
        module procedure dlogPg_ds
        module procedure dlogPg_d1
    end interface dlogPg

    !real, parameter             :: alphaMRI = 1.e-3                     ! Shakura-Sunyaev's turbulence parameter; MRI turbulence ONLY!
    real, parameter             :: q       = 0.5                        ! T propto r^-q
    real, parameter             :: p       = 1.                         ! gas surface density propto r^-p
    !real, parameter             :: T0      = 280.0                      ! temperature at 1 AU [K]
    real, parameter             :: cs0     = sqrt(kB * temperature / mH2)        ! sound speed at 1 AU
    real, parameter             :: sigmag1AU = 800.                     ! gas surface density at 1AU [g cm^-2]
    real, parameter             :: dM      = 1.e-9 * Msun / year        ! gas accretion rate
    contains

    ! Shakura-Sunyaev's turbulence parameter
    real function alpha_ds(x)
        implicit none
        real, intent(in)  :: x

        alpha_ds = alphaMRI ! this just means constant alpha

        return
    end function

    function alpha_d1(x)
        implicit none
        real, dimension(:), intent(in)  :: x
        real, dimension(size(x)) :: alpha_d1

        alpha_d1 = alphaMRI ! this just means constant alpha

        return
    end function

    ! speed of the sound in gas
    real function cs(x)
        implicit none
        real, intent(in)  :: x

        cs = cs0 * (x/(1.*AU))**(-0.5 * q)

        return
    end function

    ! Keplerian frequency
    real function omegaK(x)
        use constants, only: Ggrav, Msun
        implicit none
        real, intent(in)  :: x

        omegaK = sqrt(Ggrav * Msun / x**3)

        return
    end function

    ! gas surface density
    real function sigmag(x, time)
        use constants, only: smallv
        implicit none
        real, intent(in) :: x
        real, intent(in) :: time ! in seconds

        sigmag = sigmag1AU * (x/AU)**(-1.*p)
        sigmag = max(sigmag, smallv)
        
        return
    end function

    
    ! gas volume density
    real function densg(x,z, time)
        use constants, only: pi
        implicit none
        real, intent(in)  :: x, z
        real  :: Hg
        real, intent(in)  :: time

        Hg = cs(x) / omegaK(x) ! gas disk scaleheight
        densg = (sigmag(x, time) / (sqrt(2.*pi) * Hg)) * exp(-0.5 * (z / Hg)**2)

        return
    end function

    ! gas pressure
    real function Pg(x, z,time)
        implicit none
        real, intent(in)  :: x, z, time

        Pg = densg(x,z, time) * cs(x)**2.

        return
    end function

    ! radial gas velocity
    real function vgas(x, time)
        use constants, only: pi
        implicit none
        real, intent(in)  :: x, time

        !vgas = -1. * dM / (2. * pi * sigmag(x, time) * x)
        vgas = 0.0

        return
    end function

    ! turbulent diffusion coefficient
    real function diffcoefgas(x)
        implicit none
        real, intent(in)  :: x

        diffcoefgas = alpha(x) * cs(x)**2 / omegaK(x)

        return
    end function

    ! partial derivative of gas density with respect to vertical height
    real function ddensgdz(x,z, time)
        use constants, only: pi
        implicit none
        real, intent(in)  :: x, z, time
        real              :: Hg

        Hg = cs(x) / omegaK(x)
        ddensgdz = -z * sigmag(x, time) * exp(-0.5 * (z / Hg)**2) / (sqrt(2.*pi) * Hg**3)

        return
    end function

    ! partial derivative of gas density with respect to r
    real function ddensgdr(x,z, time)
        implicit none
        real, intent(in)  :: x, z, time
        real              :: dr

        dr = 0.01 * x
        ddensgdr = ( densg(x+dr,z, time) - densg(x-dr,z, time) ) / (2. * dr)

        return
    end function

    real function dlogPg(x, z, time)
        implicit none
        real, intent(in)  :: x, time,z
        real :: Deltar = 1.e4

        dlogPg = 0.5 * x * (Pg(x+Deltar, z, time) - Pg(x-Deltar,z,time)) / Pg(x,z,time) / Deltar

        return
    end function

    ! mass of the gas included between xmin and xmax
    real function gasmass(xmin,xmax,time)
        use constants, only: pi
        implicit none
        real, intent(in)  :: xmin,xmax, time
        real              :: x
        integer           :: i
        integer, parameter:: N = 10000
        real              :: dx

        dx = (xmax - xmin) / real(N)

        gasmass = 0.0
        x = xmin + 0.5 * dx
        do i = 1, N-1
            gasmass = gasmass + sigmag(x,time) * 2.*pi*x*dx
            x = x + dx
        enddo

        return
    end function
end module discstruct
