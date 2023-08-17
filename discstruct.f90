! put your gas disk model here
! x is distance from the star (in cm)
module discstruct
    use constants, only: AU, kB, mH2, Msun, year
    use parameters, only: alphaMRI, sigmag0, temperature, eta
    implicit none

    private
    public   :: alpha, cs, omegaK, sigmag, densg, Pg, vgas, ddensgdz, ddensgdr, gasmass, dlogPg

    !real, parameter             :: alphaMRI = 1.e-3                     ! Shakura-Sunyaev's turbulence parameter; MRI turbulence ONLY!
    real, parameter             :: q       = 0.5                        ! T propto r^-q
    real, parameter             :: p       = 1.                         ! gas surface density propto r^-p
    !real, parameter             :: T0      = 280.0                      ! temperature at 1 AU [K]
    real, parameter             :: cs0     = sqrt(kB * temperature / mH2)        ! sound speed at 1 AU
    real, parameter             :: sigmag1AU = 800.                     ! gas surface density at 1AU [g cm^-2]
    real, parameter             :: dM      = 1.e-9 * Msun / year        ! gas accretion rate
    contains

    ! Shakura-Sunyaev's turbulence parameter
    real function alpha(x)
        implicit none
        real, intent(in)  :: x

        alpha = alphaMRI ! this just means constant alpha

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

    
