#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 11:16:13 2026

@author: gurrutxaga
"""
# script to reproduce figure 6 in the paper
import numpy as np
from Plotting import init_plot
import matplotlib.pyplot as plt
import glob
import os as os
import h5py


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 15))
init_plot(ax, None, ylabel='mean ice fraction', xlabel='time [yr]')


path_list = ['../outputs/water/',
             '../outputs/water_mcdust/']
ntot = 1000
for i, (path, label) in enumerate(zip(path_list, ['This work', 'Z&D08'])):
    all_files = sorted(glob.glob(os.path.join(path, "*.h5")))
    with h5py.File(all_files[0], 'r') as f:
        ntot = int(f.attrs['number_of_particles_per_cell'])
    fw_arr = np.empty((len(all_files),ntot))
    snapt = np.empty((len(all_files),))

    if i == 0:
        mswarm_arr = np.empty((len(all_files),ntot))
    for file_id, t in enumerate(all_files):
        with h5py.File(all_files[file_id], 'r') as f:
            dset = f['swarms/swarmsout']
            swarmlist = dset[...]  # Load the whole array
            dset1 = f['times/timesout']
            snapt[file_id] = dset1[0]
            fw_arr[file_id,:] = np.array(swarmlist['fw'])[0,:]
            # Read compound dataset from swarms/swarmsout
            if i == 0:
                mswarm_arr[file_id,:] = np.array(swarmlist['mass of swarm [g]'])[0,:]

    if i == 0:
        M0 = np.mean(mswarm_arr[0,:])*mswarm_arr.shape[1]
        ax.plot(snapt, 
                 np.sum(fw_arr*mswarm_arr, axis=1)/M0,
                 label=label, lw=5)

    if i == 1:
        ax.plot(snapt, np.mean(fw_arr, axis=1), label=label, lw=5, ls='--')

ax.legend(loc='lower right', fontsize=40,  framealpha=1)


ax.set_ylim(0, 1)

fig.savefig('water_conservation.pdf', dpi=400, bbox_inches='tight', format='pdf')
