import numpy as np
import os as os
from astropy import units as u
from astropy import constants as c
from tqdm import tqdm

au = u.au.to(u.cm)
k_b = c.k_B.cgs.value
m_p = c.m_p.cgs.value
Grav = c.G.cgs.value
Msun = c.M_sun.cgs.value

class Disk_func:

    def omega_K (self, rdis) :
        """Function to calculate the orbital frequency at a given radial locatin

        Args:
            rdis (real): radial location

        Returns:
            omega_K(real): orbital frequency 
        """    
        keplerfreq = np.sqrt(Grav*Msun/(rdis)**3)
        return keplerfreq
    def temperature (self,rdis,q,t0) :
        """Function to calculate the temperature at a given radial location


        Args:
            rdis (real): radial location
            q (real): temperature power law exponent
            t0 (real): temperature at 1 AU in Kelvin.

        Returns:
            temp(real): temperature in Kelvin
        """    
        temp = t0 * (rdis/au)**(-1*q)
        return temp
    def cs_speed (self, rdis,q,t0) :
        """Function to calculate the sound speed at a given radial location

        Args:
            rdis (real): radial location
            q (real): temperature power law exponent
            t0 (real): temperature at 1 AU in Kelvin.

        Returns:
            cs(real): sound speed
        """    
        cs0 = np.sqrt(k_b*t0/(2.3*m_p))
        cs = cs0 *(rdis/au)**(-0.5*q)
        return cs
    def H_g (self, rdis,q,t0) :
        """Function to calculate the gas scale height at a given radial location

        Args:
            rdis (real): radial location
            q (real): temperature power law exponent
            t0 (real): temperature at 1 AU in Kelvin.

        Returns:
            Hg(real): scale height
        """    
        omegaK = self.omega_K(rdis)
        cs = self.cs_speed(rdis,q,t0)
        Hg = cs/omegaK
        return Hg
    def sigma_g (self,rdis,p,sigmag0,r0=au) :
        """Function to compute gas surface density at a given
        radial location.

        Args:
            rdis (real): radial location
            p (real): surface density power law exponent
            sigmag0 (real): gas surface density at radius r0
            r0 (real, optional): radial location where sigma=sigmag0. Defaults to au.

        Returns:
            sigmagas(real): gas surface density
        """    
        sigmagas  = sigmag0*(rdis/r0)**(-1*p)
        return sigmagas
    def densg (self,rdis,zdis,p,q,sigmag0,t0) :
        """Function to compute the gas volume density at a given
        radial location and vertical height.

        Args:
            rdis (real): radial location
            zdis (_type_): _description_
            p (real): surface density power law exponent
            q (real): temperature power law exponent
            sigmag0 (real): gas surface density at radius r0
            t0 (real): temperature at 1 AU in Kelvin.

        Returns:
            rho_g(real): gas volume density
        """    
        Hgas = self.H_g(rdis,q,t0)
        sigmagas = self.sigma_g(rdis,p,sigmag0)
        exponent = (zdis)**2/(2 * Hgas**2)
        rho_g = sigmagas / np.sqrt(2 * np.pi) / Hgas * np.exp(-1 * exponent) #gas density
        return rho_g
    def pressure(self,rdis,zdis,p,q,sigmag0,t0) :
        """Function to compute the gas pressure at a given
        radial location and vertical height.


        Args:
            rdis (real): radial location
            zdis (real): vertical height
            p (real): surface density power law exponent
            q (real): temperature power law exponent
            sigmag0 (real): gas surface density at radius r0
            t0 (real): temperature at 1 AU in Kelvin.
            
        Returns:
            pg(real): gas pressure
        """    
        rho_g = self.densg(rdis,zdis,p,q,sigmag0,t0)
        cs = self.cs_speed(rdis,q,t0)
        pg = rho_g*(cs**2)
        return pg

    def dlogP(self,rdis, zdis, p, q, sigmag0, t0):
        """Function to calculate the pressure gratient at a
        given radial location and vertical height

        Args:
            rdis (real): radial location
            zdis (real): vertical height
            p (real): surface density power law exponent
            q (real): temperature power law exponent
            sigmag0 (real): gas surface density at radius r0
            t0 (real): temperature at 1 AU in Kelvin.
            
        Returns:
            dlogPg(real): gas pressure gradient
        """    
        x = rdis
        z = zdis
        Deltar = 1.e4
        dlogPg = 0.5 * x * (self.pressure(x+Deltar, z,p,q,sigmag0,t0) - self.pressure(x-Deltar,z,p,q,sigmag0,t0)) / pressure(x,z,p,q,sigmag0,t0) / Deltar
        return dlogPg
 

class Diskbuild():
    """Class to build a basic static analytical disk

    """    
    def __init__(self,rbins=100,zbins=100,rmin=1,rmax=100) :
        self.r = np.linspace(rmin,rmax,rbins+1)*au
        self.rcents = 0.5*(self.r[1:]+self.r[:-1])
        self.drdis = self.r[1:]-self.r[:-1]
        self.zh = np.linspace(0,5,zbins)
        self.zwalls = np.zeros(zbins)
        self.rlog = np.logspace(np.log10(rmin),np.log10(rmax),rbins+1)*au
        self.agrid = np.logspace(-4,4,100)
        self.sigmag = np.zeros(rbins)
        self.sigmad = np.zeros(rbins)
        self.cs = np.zeros(rbins)
        self.T = np.zeros(rbins)
        self.Hg = np.zeros(rbins)
        self.omegaK = np.zeros(rbins)
        self.rhog = np.zeros((zbins,rbins))
        self.Pgas = np.zeros((zbins,rbins))
        self.alpha = None
        self.vfrag = None
        self.Stfrag = np.zeros(rbins)
        self.afrag = np.zeros(rbins)
        self.Stdrift = np.zeros(rbins)
        self.Stdf = np.zeros(rbins)
        self.adrift = np.zeros(rbins)
        self.adf = np.zeros(rbins)
        self.taugrowth = np.zeros(rbins)
        self.tgrowth = np.zeros(rbins)
        self.dtg = 0.01
        self.matdens = None
    def build(self,sigmag0=800.,p=1.,q=0.5,t0=280.,alpha=1e-3, vfrag=100.,matdens=1.):
        """Function that builds a static disk with the given parameters

        Args:
            sigmag0 (real, optional): gas surface density at 1 AU. Defaults to 800..
            p (real, optional): gas surface density power law exponent. Defaults to 1.
            q (real, optional): temperature power law exponent. Defaults to 0.5.
            t0 (real, optional): temperature at 1 AU in Kelvin. Defaults to 280..
            alpha (real, optional): alpha viscosity/turbulence strength. Defaults to 1e-3.
            vfrag (real, optional): fragmentation velocity [cm/s]. Defaults to 100..
            matdens (real, optional): material density of the particles. Defaults to 1..
        """        
        diskfunc = Disk_func()
        self.alpha = alpha
        self.vfrag = vfrag
        self.matdens = matdens
        for i in range(np.size(self.rcents)) :
            rdis = self.rcents[i]
            self.omegaK[i] = diskfunc.omega_K(rdis)
            self.sigmag[i] = diskfunc.sigma_g(rdis,p,sigmag0,au)
            self.sigmad[i] = self.dtg * self.sigmag[i]
            self.cs[i]     = diskfunc.cs_speed(rdis,q,t0)
            self.Hg[i]     = diskfunc.H_g(rdis,q,t0)
            self.T[i]      = diskfunc.temperature(rdis,q,t0)
            
            for j in range(np.size(self.zh)) :
                zdis = self.zh[j]*self.Hg[i]
                self.rhog[j,i] = diskfunc.densg(rdis,zdis,p,q,sigmag0,t0)
                self.Pgas[j,i] = diskfunc.pressure(rdis,zdis,p,q,sigmag0,t0)
            self.taugrowth = 1/self.dtg/self.omegaK[:]
            self.tmix = 1./(1e-3 * self.omegaK[:])
            self.rhog_mid  = self.sigmag/np.sqrt(2*np.pi)/self.Hg
    def vel_vn(self,p,q,sigmag0,t0):
        """_summary_

        Args:
            p (_type_): _description_
            q (_type_): _description_
            sigmag0 (_type_): _description_
            t0 (_type_): _description_

        Returns:
            _type_: _description_
        """        
        vn = np.zeros_like(self.Pgas)
        for i in range(np.size(self.rcents)):
            rdis = self.rcents[i]
            for j in range(np.size(self.z)):
                zdis = self.zh[j]*Hg[i]
                
                Pg2 = diskfunc.pressure(rdis+1,zdis,p,q,sigmag0,t0)
                Pg1 = diskfunc.pressure(rdis-1,zdis,p,q,sigmag0,t0)
                vn[j,i] = 0.25 * (Pg2-Pg1) / dens[j,i] / omega[i]
        return vn
    def calc_growthbarriers(self):
        """Function to calculate growth barriers
        """        
        tgrowth = 1/self.dtg/self.omegaK[:]
        self.Stfrag[:]  = self.vfrag**2/self.alpha/self.cs[:]**2
        gamma = 11/4
        self.Stdf[:] = self.vfrag*self.rcents[:]*self.omegaK[:]/self.cs[:]**2/0.5/gamma
        self.Stdrift[:] = self.dtg*(self.rcents[:]*self.omegaK[:])**2/self.cs[:]**2/gamma
        self.afrag[:] = self.Stfrag[:] * self.rhog_mid[:] * self.cs[:] * np.sqrt(8. / np.pi) / (self.matdens * self.omegaK[:])
        self.adf[:]   = self.Stdf[:] * self.rhog_mid[:] * self.cs[:] * np.sqrt(8. / np.pi) / (self.matdens * self.omegaK[:])
        self.adrift[:] = self.Stdrift[:] * self.rhog_mid[:] * self.cs[:] * np.sqrt(8./ np.pi) / (self.matdens * self.omegaK[:])
        a1 = np.minimum(self.afrag,self.adrift)
        self.tgrowth = self.taugrowth[:]*np.log(a1[:]/1.e-4)