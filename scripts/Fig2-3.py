#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 16 20:44:21 2026

@author: gurrutxaga
"""


# script to reproduce Figures 2 and 3 in the paper
# fig 3, high_res=True
# fig 4, high_res=False

import numpy as np
from scipy.special import gammaln
from Plotting import init_plot
import matplotlib.pyplot as plt
import glob
import os as os
import h5py



high_res = True # if False, low resolution plot




class KernelType:
    
    def __init__(self, function, res=True, mswrm=10.e20, m0=1.):
        self.function = function
        
        self.mswrm = mswrm
        self.m0 = m0
        self.m2fm = None
        self.res = res

        add_ = '_hr' if self.res else ''
        
        # read number of particles and how many times we repeated test
        parfile = '../setups/' + 'kernel_'+function+ add_ +'/setup.par'
        p = open(parfile, 'r')
        for line in p:
            s = line.split()
            label = s[0]
            if label == 'number_of_particles_per_cell':
                self.ntot = float(s[1])
            elif label == 'repeat_test':
                self.sim_num = float(s[1])
                self.sim_num = int(self.sim_num)
        
        self.file_names = []
        
        [self.file_names.append('kernel_'+function+ add_ + f'{i}') 
                 for i in range(1, self.sim_num+1)]


        if function == 'constant':
            self.analytical_kernel = self._analytical_constant_kernel
            self.t_arr = np.array([1, 10, 100, 1_000, 10_000, 100_000])
            self.xlim = (1, 10**6)
            self.nbins = 250 if self.res else 120
            self.nord = np.log10(self.mswrm * self.ntot / self.m0)
            self.mgrid = np.empty((self.nbins+1,))
            self.mgrid[0] = self.m0
            self.nord = np.log10(self.mswrm * self.ntot / self.m0)  # how many orders of magnitude in mass should the histogram go through?
            ll = self.nord / self.nbins
            lll = ll
            for i in range(1, self.nbins+1):
               self.mgrid[i] = self.m0 * 10.0 ** lll
               lll = ll * i
               
        elif function == 'linear':
            self.analytical_kernel = self._analytical_linear_kernel
            self.t_arr = np.array([4, 8, 12, 16, 20])
            self.xlim = (1, 10**10)
            self.nbins = 150 if self.res else 50
            self.nord = np.log10(self.mswrm * self.ntot / self.m0)
            self.mgrid = np.empty((self.nbins+1,))
            self.mgrid[0] = self.m0

            self.nord = np.log10(self.mswrm * self.ntot / self.m0)  # how many orders of magnitude in mass should the histogram go through?

            ll = self.nord / self.nbins
            lll = ll
            for i in range(1, self.nbins+1):
               self.mgrid[i] = self.m0 * 10.0 ** lll
               lll = ll * i
        elif function == 'product':
            self.analytical_kernel = self._analytical_product_kernel
            self.t_arr = np.array([0.4, 0.7, 0.9])
            self.xlim = (1, 10**3)
            self.nbins = 30  if self.res else 10
            self.nord = 3
            self.mgrid = np.logspace(np.log10(self.m0), self.nord,
                                     self.nbins+1).astype(int)
        else:
            raise ValueError(f'Unknown kernel type: {function}')
        
        self.mgrid_mid = np.sqrt(self.mgrid[:-1] * self.mgrid[1:] )

    # from Tanaka & Nakazawa 1994, their eq. 2.1
    def _analytical_constant_kernel(self, t, m, a=1.):
        m0 = m[0]
        N0 = 1./m0
        N = N0/m0*4./(a*N0*t)**2 *np.exp((1.-m/m0)*2/(a*N0*t))
        return N*m**2

    # from Tanaka & Nakazawa 1994, their eq. 2.2
    def _analytical_linear_kernel(self, t, m, a=0.5):
        m0 = m[0]
        N0 = 1./m0**2
        g = np.exp(-a*t)
        k = (m / m0)
        
        # compute it with logarithm because large numbers give numerical problems
        logN = np.log(N0*g) - k*(1.-g) + (k-1.)*np.log(k*(1.-g)) - gammaln(k+1.)
        N = np.exp(logN) 
        return N*m**2

    
    # from Tanaka & Nakazawa 1994, their eq. 2.3
    def _analytical_product_kernel(self, t, m):
        # compute it with logarithm because large numbers give numerical problems
        logN = (m-1)*np.log(m*t) - m*t - gammaln(m+1.) - np.log(m)
        N = np.exp(logN) 
        return N*m**2
    
    def compute_m2fm(self):
        if len(self.file_names) == 0:
            raise ValueError('Initiate file_names')
        self.m2fm = np.zeros((self.nbins, len(self.t_arr),
                         len(self.file_names)))
        for sim, filename in enumerate(self.file_names):
            # read files
            output_dir = '../outputs/' + filename
            all_files = sorted(glob.glob(os.path.join(output_dir, "*.h5")))
            for file_id, t in enumerate(all_files):
                with h5py.File(all_files[file_id], 'r') as f:
                    dset = f['swarms/swarmsout']
                    swarmlist = dset[...]  # Load the whole array
                    
                    # Read compound dataset from swarms/swarmsout
                    npar_arr = swarmlist['npar'][...][0]
                    mass_arr = swarmlist['mass of a particle [g]'][0]
    
                    print(file_id)
                    for i in range(0, mass_arr.shape[0]):
                    
                        m = mass_arr[i]
                        w = npar_arr[i]
                    
                        # find bin: mgrid[j] <= m < mgrid[j+1]
                        j = np.searchsorted(self.mgrid, m, side="right") - 1
                        # skip particles outside the grid
                        if j>self.mgrid.shape[0]-2: continue

                        dm = self.mgrid[j+1] - self.mgrid[j]
                        self.m2fm[j, file_id, sim] += w * m**2 / dm  
        self.m2fm /= (self.mswrm*self.ntot)
    


colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

sim_list = ['constant', 'linear', 'product']

kernel_list = []
[kernel_list.append(KernelType(sim, res=high_res)) for i, sim in enumerate(sim_list)]
    

fig, ax = plt.subplots(1, 3, figsize=(50,15), sharey=True)
lw = 3

for i, axx in enumerate(ax): 
    ylabel = r'm$^{2}$f(m)' if i==0 else None
    kernel = kernel_list[i]
    init_plot(axx, xlabel='m', ylabel=ylabel, title=kernel.function.capitalize())
    if high_res:
        axx.set_ylim(10**-5, 1)
    else:
        axx.set_ylim(10**-3, 1)
        
    axx.set_xlim(kernel.xlim)
    
    axx.set_yscale('log')
    axx.set_xscale('log')


    kernel.compute_m2fm()
    start = 1
    for it, time in enumerate(kernel.t_arr):
        c = colors[it % len(colors)]
        axx.errorbar(kernel.mgrid_mid[start:], 
                     np.mean(kernel.m2fm[start:, it,:], axis=1),
                     yerr=np.std(kernel.m2fm[start:, it, :], axis=1),
                     fmt='none', marker='o', c=c, capsize=3, ls=None)
        axx.scatter(kernel.mgrid_mid[start:], 
                     np.mean(kernel.m2fm[start:, it,:], axis=1),
                     marker='o', c=c)
        axx.plot(kernel.mgrid, kernel.analytical_kernel(time, kernel.mgrid),
                 c=c, ls='--')



fig.subplots_adjust(wspace=0.15, hspace=0.1, top=0.85, right=0.77)

if high_res:
    fig.savefig('kernel_highres.pdf', dpi=400, bbox_inches='tight', format='pdf')
else:
    fig.savefig('kernel_lowres.pdf', dpi=400, bbox_inches='tight', format='pdf')

