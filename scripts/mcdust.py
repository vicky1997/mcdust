import numpy as np
import os as os
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.constants import golden
from astropy import units as u
from astropy import constants as c
import glob
from tqdm import tqdm
import pandas as pd
from mcdisk import Diskbuild
import warnings

au = u.au.to(u.cm)
k_b = c.k_B.cgs.value
m_p = c.m_p.cgs.value
Grav = c.G.cgs.value
Msun = c.M_sun.cgs.value
aH2 = 2e-15
plt.rcParams.update({'font.size': 12})

class Swarm:
    """Class for reading/handling data produced by mcdust
    """    
    def __init__(self, par):
        """Function to initialize variables of the class

        Args:
            par (Params): object containing the parameters of the simulation
        """        
        self.idnr = np.zeros((par.ntime,par.nswarms))
        self.mass = np.zeros_like(self.idnr)
        self.rdis = np.zeros_like(self.idnr)
        self.zdis = np.zeros_like(self.idnr)
        self.St = np.zeros_like(self.idnr)
        self.grain_size = np.zeros_like(self.idnr)
        self.indens = np.ones_like(self.idnr)
        self.velr = np.zeros_like(self.idnr)
        self.velz = np.zeros_like(self.idnr)
        self.mswarm = None
        self.snapt = np.zeros(par.ntime)
        self.sigmad = None
        self.rwalls = None
        self.rcents = None
    def read_hdf5(self, par) :
        """Function 

        Args:
            par (Params): object containing the parameters of the simulation
        """        
        snapshottime = np.zeros(par.ntime)
        f1 = os.path.join(par.datadir+'swarms-00000.h5')
        with h5py.File(f1,"r") as f:
            a1 = f.attrs['authors'].decode('UTF-8')
            a2 = f.attrs['code'].decode('UTF-8')
            a3 = f.attrs['mass_of_swarm[g]']
            a4 = f.attrs['output_number']
            print(str(a2))
            print("Authors : " + str(a1))
            self.mswarm = a3
            f.close()
        print("Please cite Drążkowska, Windmark & Dullemond (2013) A&A 556, A37")
        print("Reading data ...")
        for iout in tqdm(range(par.ntime)):
            string="%05d"%(iout)
            fname = os.path.join(par.datadir+'swarms-'+string+'.h5')
            with h5py.File(fname,"r") as f:
                dset = f['swarms/swarmsout']
                dset1 = f['times/timesout']
                swarmlist = dset[...]
                timesout = dset1[...]
                f.close()
            self.idnr[iout,:] = np.array(int(swarmlist['id number']))
            self.mass[iout,:] = np.array(swarmlist['mass of a particle [g]'])
            self.rdis[iout,:] = np.array(swarmlist['cylindrical radius [AU]'])
            self.zdis[iout,:] = np.array(swarmlist['height above midplane [AU]'])
            self.St[iout,:]   = np.array(swarmlist['Stokes number'])
            self.velr[iout,:] = np.array(swarmlist['Radial velocity v_r [cm/s]'])
            self.velz[iout,:] = np.array(swarmlist['Vertical velocity v_z [cm/s]'])
            snapshottime[iout]= timesout
        self.indens = self.indens*par.rho_s
        con = np.zeros_like(self.St)
        con = (0.75/np.pi/self.indens)**(1./3.)
        mass_to_one_third = np.zeros_like(con)
        mass_to_one_third = self.mass[:,:]**(1./3.)
        self.grain_size = np.multiply(mass_to_one_third,con)
        self.snapt = snapshottime
        print("Done!")
    def calculate_properties(self,pars,time=-1):
        """Function to calculate dust surface density

        Args:
            pars (Params): object containing the parameters of the simulation
            time (int, optional): time of the simumation. Defaults to -1.
        """        
        rbins = 100
        rwalls = np.linspace(0.99*pars.minr,pars.maxr+0.1, rbins)*au
        rcents = 0.5*(rwalls[1:]+rwalls[:-1])
        drdis = rwalls[1:]-rwalls[:-1]
        if(time==-1):
            time = self.snapt[-1]
        it = self.snapt.searchsorted(time)
        rdis = self.rdis[it,:]
        rcounts,binsr = np.histogram(rdis,bins=rwalls)
        self.sigmad = rcounts[:-1] *self.mswarm / (2 * np.pi * rcents*drdis)
        self.rwalls = rwalls
        self.rcents = rcents
    def describe(self,time=-1):
        """Function that prints the mean, median, max and min values for the
           properties of the swarms.

        Args:
            time (int, optional): Time of the simulation. Defaults to -1.
        """        
        if(time==-1):
            it = -1
        else :
            it = self.snapt.searchsorted(time)
        frame = {"mass [g]":self.mass[it,:],"Stokes Nr":self.St[it,:], "grain_size [cm]":self.grain_size[it,:]}
        df = pd.DataFrame(frame)
        print(df.describe()) 
    

class Params :
    """Class to handle the parameters for a simulation in mcdust
    """    
    def __init__(self):
        """Function to initialize the parameters to be handled.
        """        
        self.alpha_t = None
        self.tgas = None
        self.sigmagas = None
        self.minr = None
        self.maxr = None
        self.a0 = None
        self.vfrag = None  
        self.rho_s = None
        self.dtg = None
        self.erosion_m_ratio = None
        self.t_end = None
        self.nr = None
        self.nz = None
        self.n_cell = None
        self.nswarms = None
        self.ntime = None
        self.eta = None
    def read_params(self, path, parfile='setup.par'):
        """Function to read parameters from the specified path and
           and parameter file"

        Args:
            path (str): path to the directory of the par file
            parfile (str, optional): Name of the par file. Defaults to 'setup.par'.
        """        
        parfile = os.path.join(path+parfile)
        p = open(parfile, 'r')
        for line in p:
            s = line.split()
            label = s[0]
            if label == 'alpha':
                self.alpha = float(s[1])
            elif label == 'temperature_[K]':
                self.tgas = float(s[1])
            elif label == 'sigma_gas_[g/cm2]':
                self.sigmagas = float(s[1])
            elif label == 'minimum_radius_[AU]':
                self.minr = float(s[1])
            elif label == 'maximum_radius_[AU]':
                self.maxr = float(s[1])
            elif label == 'monomer_radius_[cm]':
                self.a0 = float(s[1])
            elif label == 'fragmentation_velocity_[cm/s]':
                self.vfrag = float(s[1])
            elif label == 'material_density_[g/cm3]':
                self.rho_s = float(s[1])
            elif label == 'dust_to_gas_ratio':
                self.dtg = float(s[1])
            elif label == 'number_of_particles_per_cell':
                self.n_cell = int(s[1])
            elif label == 'number_of_vertical_zones':
                self.nz = int(s[1])
            elif label == 'number_of_radial_zones':
                self.nr = int(s[1])
            elif label == 'erosion_mass_ratio':
                self.erosion_m_ratio = float(s[1])
            elif label == 'data_directory':
                self.datadir = s[1]
                self.datadir = self.datadir[1:-1]
                
        self.nswarms = self.nr*self.nz*self.n_cell
        self.datadir = path + self.datadir + '/'
        nroutputs = len(glob.glob1(self.datadir+'/',"*.h5"))   
        self.ntime = nroutputs
    def read_params_hdf5(self,path):
        """Read parameters directly from the hdf5 file

        Args:
            path (str): path to the hdf5 data files.
        """    
        self.datadir = os.path.join(path + '/')
        path_file0 = os.path.join(path + 'swarms-00000.h5')
        print(' Reading paramters file from data files in '+ path )
        with h5py.File(path_file0,"r") as f:
            self.nswarms = int(f.attrs['number_of_particles'])
            self.n_cell = int(f.attrs['number_of_particles_per_cell'])
            self.nr = int(f.attrs['number_of_radial_zones'])
            self.nz = int(f.attrs['number_of_vertical_zones'])
            self.t_end = f.attrs['maximum_time_of_simulation[yr]']
            self.minr = f.attrs['minimum_radius_[AU]']
            self.maxr = f.attrs['maximum_radius_[AU]']
            self.a0 = f.attrs['monomer_radius_[cm]']
            self.rho_s = f.attrs['material_density_[g/cm3]']
            self.dtg = f.attrs['dust_to_gas_ratio']
            self.vfrag = f.attrs['fragmentation_velocity_[cm/s]']
            self.alpha_t = f.attrs['alpha_t']
            self.sigmagas = f.attrs['sigma_gas_[g/cm2]']
            self.tgas = f.attrs['temperature_[K]']
            self.eta = f.attrs['pressure_gradient_eta']
            self.erosion_m_ratio = f.attrs['erosion_mass_ratio']
            self.alpha = f.attrs['alpha_t']
            f.close()
        self.ntime = len(glob.glob1(self.datadir+'/',"*.h5"))    
#TODO: function to write parameter files needed to start the code


class Simulation:
    """Class simulation to handle the data and the parameters of the simulation along with a simple
       disk model built with the parameters of the simulation.
    """
    def __init__(self):
        self.swarms = None
        self.pars = None
        self.disk = None
    def read(self,path):
        """Function to read the data files, parameter files and create a simple disk model.
           The parameters are read from the hdf5 files. If you want it to be read from the setup.par
           file in the setups directory, it has be to modified.

        Args:
            path (str): path to data directory
        """        
        #self.pars.read_params(path)
        #self.swarms.read_hdf5(self.pars)
        self.pars = Params()
        self.pars.read_params_hdf5(path)
        self.swarms = Swarm(self.pars)
        self.swarms.read_hdf5(self.pars)
        self.disk = Diskbuild(rmin=self.pars.minr,rmax=self.pars.maxr)
        self.disk.build(sigmag0=self.pars.sigmagas,t0=self.pars.tgas,alpha=self.pars.alpha, matdens=self.pars.rho_s)


#Plotting routines

class Plots:
    def rz_scatter(swarms, time=-1, scatter=np.zeros(100), step=1, cmap = 'viridis', cbarlabel='') :
        """Function to produce a scatter plot optionally with a quantity to be shown
           with color.

        Args:
            swarms (Swarms): Object to store swarm data
            time (int, optional): Time of the simulation. Defaults to -1.
            scatter (real, optional): Quantity to be denoted with the color. Defaults to np.zeros(100).
            step (int, optional): Plot every nth element. Defaults to 1.
            cmap (str, optional): Colomarp. Defaults to 'viridis'.
            cbarlabel (str, optional): Colorbar label. Defaults to ''.
        """        
        if(time==-1):
            time = swarms.snapt[-1]
        it = swarms.snapt.searchsorted(time)

        f,ax = plt.subplots()
        if scatter.all() ==0 :
            ax.scatter(swarms.rdis[it,0::step], swarms.zdis[it,0::step])

        else:
            sc = ax.scatter(swarms.rdis[it,0::step], swarms.zdis[it,0::step],c=scatter[it,0::step],norm=LogNorm())
            cbar = plt.colorbar(sc, cmap=cmap)   
            cbar.ax.set_ylabel(cbarlabel)
            ax.set_xscale('log')
        ax.set_xlabel('distance from central star [au]')
        ax.set_ylabel('height above midplane [au]')
        ax.set_title('Time = %d yrs'%(swarms.snapt[it])) 
    def scatter_all(swarms, time=-1, step=1) :
        """Function to produce scatter plots for Stokes Number, mass,
           grain size, material density.

        Args:
            swarms (Swarms): Object to store swarm data
            time (int, optional): Time of the simulation. Defaults to -1.
            step (int, optional): Plot every nth element. Defaults to 1.
        """        
        if(time==-1):
            time = swarms.snapt[-1]
        it = swarms.snapt.searchsorted(time)
        cmaps = ['magma_r','afmhot_r','viridis','viridis']
        cbarlabel = ['Stokes Number','Mass [g]', 'grain size [cm]', 'material density [g/cm3]']
        scatter = [swarms.St, swarms.mass, swarms.grain_size, swarms.indens]
        layout = """AB;CD"""
        f,axd = plt.subplot_mosaic(layout,figsize=((10,8)))
        ax = [axd["A"],  axd["B"], axd["C"], axd["D"]]
        for axs,cm,cb,scat in zip(ax,cmaps,cbarlabel,scatter):
            sc = axs.scatter(swarms.rdis[it,0::step], swarms.zdis[it,0::step],c=scat[it,0::step], norm=LogNorm())
            cbar = plt.colorbar(sc, cmap=cm,ax=axs)
            cbar.ax.set_ylabel(cb)
            axs.set_xscale('log')
        axd["D"].set_xlabel('distance from star [au]')
        axd["C"].set_xlabel('distance from star [au]')
        axd["A"].set_ylabel('z [au]')
        axd["C"].set_ylabel('z [au]')
        plt.suptitle('Time = %d yrs'%(swarms.snapt[it]))
        plt.tight_layout()
    def radial_scatter_all(swarms, time=-1, step=1) :
        """Function to produce radial scatter plots for Stokes Number
           Mass and grain size.

        Args:
            swarms (Swarms): Object to store swarm data
            time (int, optional): Time of the simulation. Defaults to -1.
            step (int, optional): Plot every nth element. Defaults to 1.
        """        
        if(time==-1):
            time = swarms.snapt[-1]
        it = swarms.snapt.searchsorted(time)
        ylabel = ['Stokes Number','Mass [g]', 'grain size [cm]']
        scatter = [swarms.St, swarms.mass, swarms.grain_size]
        layout = """ABC"""
        f,axd = plt.subplot_mosaic(layout,figsize=((9,3)))
        ax = [axd["A"],  axd["B"], axd["C"],]
        for axs,yl,scat in zip(ax,ylabel,scatter):
            axs.scatter(swarms.rdis[it,0::step], scat[it,0::step])
            axs.set_ylabel(yl)
            axs.set_xlabel('distance from star [au]')
            axs.set_xscale('log')
            axs.set_yscale('log')
        plt.suptitle('Time = %d yrs'%(swarms.snapt[it]))
        plt.tight_layout()
    def sigma_d(swarms,pars,disk,time=-1):
        """Function to produce plot of dust surface density vs radius

        Args:
            swarms (Swarms): Object to store swarm data
            pars (Params): object containing the parameters of the simulation
            disk (Diskbuild): object containing the static disk model
            time (int, optional): time of the simulation. Defaults to -1.
        """        
        #rbins = 100
        #rwalls = np.linspace(0.99*pars.minr,pars.maxr+0.1, rbins)*au
        #rcents = 0.5*(rwalls[1:]+rwalls[:-1])
        #drdis = rwalls[1:]-rwalls[:-1]
        if(time==-1):
            time = swarms.snapt[-1]
        it = swarms.snapt.searchsorted(time)
        width = 8.
        dpi=150.
        f,ax = plt.subplots(figsize=(width/2., width/2),dpi=dpi)
        ax.loglog(swarms.rcents/au,swarms.sigmad,label='2DMC')
        ax.loglog(disk.rcents/au,0.01*disk.sigmag,label='analytical')
        ax.legend()
        ax.set_xlabel('r [au]',fontsize='x-large')
        ax.set_ylabel('$\Sigma_d$',fontsize='x-large')
        plt.show()
    def mass_dens_radial(swarms,pars,time=-1,nmbins=32,rbins=32,show_hist=False):
        """Function to produce density plots of mass vs radius

        Args:
            swarms (Swarms): Object to store swarm data
            pars (Params): object containing the parameters of the simulation
            time (int, optional): time of the simulation. Defaults to -1.
            nmbins (int, optional): number of mass bins. Defaults to 100.
            rbins (int, optional): number of radial bins. Defaults to 100.
            show_hist (bool, optional): _description_. Defaults to True.
        """        
        if(time==-1):
            time = swarms.snapt[-1]
        it = swarms.snapt.searchsorted(time)
        mwalls = np.logspace(np.log10(np.min(swarms.mass)),np.log10(np.max(swarms.mass)),nmbins+1)
        mcents = np.sqrt(mwalls[1:]*mwalls[:-1])
        rlogwalls = np.logspace(np.log10(pars.minr),np.log10(pars.maxr),rbins+1)*au
        rlogcents = np.sqrt(rlogwalls[1:]*rlogwalls[:-1])
        drlog = rlogwalls[1:]-rlogwalls[:-1]
        width=7.
        it=-1
        mdens = np.zeros((pars.ntime,rbins,nmbins))
        mdens1 = np.zeros((rbins,nmbins))
        f,ax = plt.subplots(figsize=(width/2., width/2), dpi=300)
        ax.set_xlabel('distance from star [au]')
        ax.set_ylabel('particle mass [g]')
        ax.set_xscale('log')
        ax.set_yscale('log')
        counts,xedges,yedges,im=ax.hist2d(swarms.rdis[it,:],swarms.mass[it,:],bins=(rlogwalls/au,mwalls),norm=LogNorm(),cmap='magma_r')
        for i in range(np.size(rlogwalls)-1):
            mdens1[i,:] = counts[i,:]* swarms.mswarm /(2.*np.pi*rlogcents[i]*drlog[i])
        f.colorbar(im)
        time = int(swarms.snapt[it])
        ax.set_title('t = {:.1f} yr'.format(swarms.snapt[it]))
        if(show_hist==True):
            plt.show()
        f,ax=plt.subplots(figsize=(width/1.5, width/2.), dpi=200)
        ax.set_ylim(np.min(swarms.mass),np.max(swarms.mass)+5)
        ax.set_xlim(0.99,101.)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylabel("particle mass [g]")
        ax.set_xlabel("distance from star [au]")

        sigMax = 3.
        sigMin = -7.
        levels = np.arange(sigMin, sigMax, 1)
        pcf = ax.contourf(rlogcents/au, mcents, np.log10(mdens1.T), levels=levels, extend="both", cmap="magma_r")
        cbar = plt.colorbar(pcf, ax=ax)
        cbar.ax.set_ylabel(r"$\log\ \Sigma_\mathrm{d}$ [g/cm²]")
        ax.set_title('t = {:.1f} yr'.format(swarms.snapt[it]))
        ax.tick_params(width=1.,which='both')
        plt.show()
        
        

#todo : write routine to make video