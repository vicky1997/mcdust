! put your gas disk model here
! x is distance from the star (in cm)
module discstruct
   use constants, only: AU, kB, mH2, Msun, year
   use parameters, only: alphaMRI, sigmag0, temperature, eta
   implicit none

   private
   public   :: alpha, cs, omegaK, sigmag, densg, Pg, vgas, ddensgdz, ddensgdr, gasmass, dlogPg

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

      cs = sqrt(kB * temperature / mH2) !* (x/(1.*AU))**(-0.5 * q)

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
      use constants, only: pi, smallv
      implicit none
      real, intent(in) :: x
      real, intent(in) :: time ! in seconds

      sigmag = sigmag0

      return
   end function

   ! gas disk scaleheight
   real function Hg(x)
      implicit none
      real, intent(in)  :: x

      Hg = cs(x) / omegaK(x)

      return
   end function

   ! gas volume density
   real function densg(x,z, time)
      use constants, only: pi
      implicit none
      real, intent(in)  :: x, z
      real, intent(in)  :: time

      densg = (sigmag(x, time) / (sqrt(2.*pi) * Hg(x))) * exp(-0.5 * (z / Hg(x))**2)

      return
   end function

   ! gas pressure
   real function Pg(x, z,time)
      use constants, only: pi
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

      vgas = 0.0

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
end
