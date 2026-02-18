#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 11:46:14 2026

@author: gurrutxaga
"""

# script to reproduce Figure 4 in the paper of global disk

import numpy as np
from Plotting import init_plot
import matplotlib.pyplot as plt
import glob
import os as os
import h5py
from astropy import units as u



time = 10e4 # years
nmbins = 100
rbins = 100

au = u.au.to(u.cm)
width=30.

sigMax = 2.
sigMin = -5.
levels = np.arange(sigMin, sigMax, 1)

size_walls = np.logspace(np.log10(10**(-4)),
                     np.log10(20), nmbins+1)
size_cents = np.sqrt(size_walls[1:]*size_walls[:-1])


# Radial bins
rlogwalls = np.logspace(np.log10(1.), 
                        np.log10(100.), rbins+1)*au
rlogcents = np.sqrt(rlogwalls[1:]*rlogwalls[:-1])
drlog = rlogwalls[1:]-rlogwalls[:-1]

path_list = ['../outputs/global_disk_mcdust/',
             '../outputs/global_disk_X1/',
             '../outputs/global_disk_X2/',
             '../outputs/global_disk_X3/']



fig,axx = plt.subplots(nrows=2, ncols=2, figsize=(width, 0.8*width), dpi=300,
                     sharex=True, sharey=True)


for i_plot, (path, ax) in enumerate(zip(path_list,axx.flatten())):
    init_plot(ax)

    # read last snapshot
    all_files = sorted(glob.glob(os.path.join(path, "*.h5")))


    with h5py.File(all_files[-1], 'r') as f:
        dset = f['swarms/swarmsout']
        swarmlist = dset[...]  # Load the whole array
        dset1 = f['times/timesout']
        rdis = np.array(swarmlist['cylindrical radius [AU]'])[0,:]
        mass = np.array(swarmlist['mass of a particle [g]'])[0,:]

        # internal density of dust is 1g/cm**3
        grain_size = (0.75/np.pi*mass)**(1./3.)
        # Read compound dataset from swarms/swarmsout
        if i_plot == 0:
            mswarm = f.attrs['mass_of_swarm[g]']
        else:
            mswarm = np.array(swarmlist['mass of swarm [g]'])[0,:]
            


    mdens = np.zeros((rbins,nmbins))
    mdens1 = np.zeros((rbins,nmbins))
    

    if i_plot!=0:
        counts,xedges,yedges=np.histogram2d(rdis, grain_size,
                                            bins=(rlogwalls/au,size_walls),
                                            weights=mswarm)
    else:
        counts,xedges,yedges=np.histogram2d(rdis, grain_size,
                                            bins=(rlogwalls/au,size_walls))
        counts *= mswarm
   # Surface density
    mdens1 = counts / (2. * np.pi * rlogcents[:, None] * drlog[:, None])
    mdens1[mdens1 == 0] = 10**(-14)


    # Contour plot for second subplot

    pcf = ax.contourf(rlogcents / au, size_cents, np.log10(mdens1.T),
                      levels=levels, extend="both", cmap="magma_r")

    if i_plot==0:
        # Axes settings
        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.set_ylabel("dust size [cm]")
        ax.set_xlim(0.95, 100.01)
        ax.set_ylim(1e-4, 20)

        ax.tick_params(width=1., which='both')
        ax.set_title('mcdust')
        
    elif i_plot==1:
        ax.set_title('X=0.1')
    elif i_plot==2:
        ax.set_ylabel("dust size [cm]")
        ax.set_title('X=0.01')
        ax.set_xlabel("distance from star [au]")
    else:
        ax.set_title('X=0.001')
        ax.set_xlabel("distance from star [au]")

cbar = fig.colorbar(pcf, ax=axx)
cbar.set_label(r'$\log \left(\frac{\sigma_\mathrm{d}(a)}{g/cm^{2}}\right)$')#/ [g/cm²]

fig.subplots_adjust(wspace=0.1, hspace=0.15, top=0.85, right=0.77)
fig.savefig('mcdust.pdf', dpi=400, bbox_inches='tight', format='pdf')
