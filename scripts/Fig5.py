#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 11:16:13 2026

@author: gurrutxaga
"""
# script to reproduce Figure 5 in the paper
import numpy as np
from Plotting import init_plot
import matplotlib.pyplot as plt
import glob
import os as os
import h5py


"""
l--------- l ------ l ------ l
l          l        l        l
l          l   ax1  l   ax2  l
l          l -------l--------l
l          l        l        l
l    ax0   l   ax3  l   ax4  l
l          l -------l--------l
l          l        l        l
l          l   ax5  l   ax6  l
l--------- l ------ l ------ l
"""
from matplotlib.gridspec import GridSpec
fig = plt.figure(figsize=(30, 15))

gs = GridSpec(3, 4, figure=fig, width_ratios=[1, 1, 1, 1])  


ax0 = fig.add_subplot(gs[:,:2])
ax1 = fig.add_subplot(gs[0,2])
ax2 = fig.add_subplot(gs[0,3])
ax3 = fig.add_subplot(gs[1,2])
ax4 = fig.add_subplot(gs[1,3])
ax5 = fig.add_subplot(gs[2,2])
ax6 = fig.add_subplot(gs[2,3])
plt.subplots_adjust(hspace=0.2)

#fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15, 15))
init_plot(ax0, None, ylabel='mean ice fraction',
          xlabel='time [yr]')


rhow, rhos = 1., 3. # densities of water and silicates
a_th = 0.1 # size threshold in cm
lw = 10 # line width
marker_size = 100 
file_id_list = [900, 1010, 1500] # which snapshots to plot
path_list = ['../outputs/water2/', '../outputs/water_mcdust/']
# get list of default colors in matplotlib
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

for i, (path, label) in enumerate(zip(path_list, ['This work', 'Z&D08'])):
    all_files = sorted(glob.glob(os.path.join(path, "*.h5")))
    with h5py.File(all_files[0], 'r') as f:
        ntot = int(f.attrs['number_of_particles_per_cell'])
    fw_arr = np.empty((len(all_files),ntot))
    m_arr = np.empty((len(all_files),ntot))
    grain_arr = np.empty((len(all_files),ntot))
    fw_mean = np.empty((len(all_files),2))
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
            w = fw_arr[file_id,:]/(1-fw_arr[file_id,:])
            
            # Read compound dataset from swarms/swarmsout
            if i == 0:
                mswarm_arr[file_id,:] = np.array(swarmlist['mass of swarm [g]'])[0,:]
                m_arr[file_id,:] = np.array(swarmlist['mass of a particle [g]'])[0,:]
                rhoi = 1/(fw_arr[file_id,:]/rhow + (1-fw_arr[file_id,:])/rhos)
                grain_arr[file_id,:] = (m_arr[file_id,:]/(4*np.pi/3*rhoi))**(1/3)
                bool_a = grain_arr[file_id,:]>=a_th
                fw_mean[file_id,0] = np.sum(w[bool_a]*mswarm_arr[file_id,bool_a])/np.sum(mswarm_arr[file_id,bool_a]*(1.+w[bool_a]))
                fw_mean[file_id,1] = np.sum(w[~bool_a]*mswarm_arr[file_id,~bool_a])/np.sum((1+w[~bool_a])*mswarm_arr[file_id,~bool_a])
            else:
                m_arr[file_id,:] = np.array(swarmlist['mass of a particle [g]'])[0,:]
                rhoi = 1/(fw_arr[file_id,:]/rhow + (1-fw_arr[file_id,:])/rhos)
                grain_arr[file_id,:] = (m_arr[file_id,:]/(4*np.pi/3*rhoi))**(1/3)
                bool_a = grain_arr[file_id,:]>=a_th
                fw_mean[file_id,0] = np.sum(w[bool_a])/np.sum(1.+w[bool_a])
                fw_mean[file_id,1] = np.sum(w[~bool_a])/np.sum(1+w[~bool_a])


    if i == 0:
        M0 = np.mean(mswarm_arr[:,:]/(1-fw_arr[:,:]), axis=1)
        ax0.plot(snapt, 
                 np.mean(fw_arr/(1-fw_arr)*mswarm_arr, axis=1)/M0,
                 label=label, lw=lw, c=colors[i])
        label = [None, None]
        ax0.plot(snapt, fw_mean[:,0], label=label[0], lw=lw/2, c=colors[i],
                ls='--', zorder=-1, alpha=0.5)

        ax0.plot(snapt, fw_mean[:,1], label=label[1], lw=lw/2, ls=':',
                c=colors[i], zorder=-1, alpha=0.5)

        # mark times
        for ww in range(0, len(file_id_list)):
            x_timestep = snapt[file_id_list[ww]]
            ax0.plot(x_timestep, 1., 'o',
                        markersize=20,
                        #markerfacecolor='none',   # hollow circle
                        c='blueviolet',
                        #markeredgewidth=lw,
                        transform=ax0.get_xaxis_transform(),
                        clip_on=False)

    if i == 1:
        ax0.plot(snapt, 
                 np.mean(fw_arr/(1-fw_arr), axis=1)/np.mean(1/(1-fw_arr[:,:]), axis=1),
                 label=label, lw=lw, c=colors[i])
        ax0.plot(snapt, fw_mean[:,0], lw=lw/2, c=colors[i],
                ls='--', zorder=-1, alpha=0.5)

        ax0.plot(snapt, fw_mean[:,1], lw=lw/2, ls=':',
                c=colors[i], zorder=-1, alpha=0.5)

ax0.plot([],[], c='k', ls='-', label='all sizes', lw=lw)
ax0.plot([],[], c='k', ls=':', label=f'a$\leq${a_th}$\,$cm', lw=lw/2)
ax0.plot([],[], c='k', ls='--', label=f'a>{a_th}$\,$cm', lw=lw/2)
#ax.hlines(y=0.5, xmin=0, xmax=3000)
ax0.legend(loc='upper left', fontsize=40,  framealpha=1)

ax0.set_xlim(0, 1950)
ax0.set_ylim(-0.1, 1.1)






    #color = ['#609D9B','orange', '#543B87']
for i, (path, label, ax_list) in enumerate(zip(path_list, ['This work', 'Z&D08'], [[ax1, ax3, ax5], [ax2, ax4, ax6]])):
    all_files = sorted(glob.glob(os.path.join(path, "*.h5")))
    with h5py.File(all_files[0], 'r') as f:
        ntot = int(f.attrs['number_of_particles_per_cell'])
    snapt = np.empty((len(all_files),))


    for j, (file_id, ax) in enumerate(zip(file_id_list, ax_list)):
        t = float(file_id)
        init_plot(ax)
        with h5py.File(all_files[file_id], 'r') as f:
            dset = f['swarms/swarmsout']
            swarmlist = dset[...]  # Load the whole array
            dset1 = f['times/timesout']
            snapt[file_id] = dset1[0]
            fw_arr = np.array(swarmlist['fw'])[0]
            mass_arr = swarmlist['mass of a particle [g]'][0]
            # Read compound dataset from swarms/swarmsout
    
            rhoi = 1/(fw_arr/rhow + (1-fw_arr)/rhos)
            grain_size = (mass_arr/(4*np.pi/3*rhoi))**(1/3)  

            ax.scatter(grain_size, fw_arr, label=f'{snapt[file_id]-1000:.0f} yr',
                       s=marker_size, alpha=0.5, c=colors[i])
 
        # horizontal line as a reference
        ax.axhline(y=0.5, xmin=0.5e-4, xmax=10, zorder=-1, color='lightgray',
                   lw=lw/2, ls='--')
        ax.set_xscale('log')
        ax.set_xlim(0.5e-4, 10.)
        ax.set_ylim(-0.1, 1.1)
        
        if i==1:
            round_val = round(snapt[file_id_list[j]] / 10) * 10
            ax.text(0.1, 0.8, s=f'{round_val:.0f}$\,$yr', color='blueviolet',
                    fontsize=33)


# adjusment of the plots 
ax0.tick_params(top=True, labeltop=True, bottom=True, labelbottom=True)


ax1.tick_params(left=True, labelleft=True, right=True, labelright=False,
                top=True, labeltop=True, bottom=True, labelbottom=False)
ax3.tick_params(left=True, labelleft=True, right=True, labelright=False,
                top=True, labeltop=False, bottom=True, labelbottom=False)
ax5.tick_params(left=True, labelleft=True, right=True, labelright=False,
                top=True, labeltop=False, bottom=True, labelbottom=True)


ax2.tick_params(left=True, labelleft=False, right=True, labelright=True,
                top=True, labeltop=True, bottom=True, labelbottom=False)
ax4.tick_params(left=True, labelleft=False, right=True, labelright=True,
                top=True, labeltop=False, bottom=True, labelbottom=False)
ax6.tick_params(left=True, labelleft=False, right=True, labelright=True,
                top=True, labeltop=False, bottom=True, labelbottom=True)




# label at the bottom of subplots
ax5.set_xlabel('dust size [cm]')

# 
ax5.xaxis.set_label_coords(1.1, -0.25)

# label on the left size of subplots

pos2 = ax2.get_position()
pos4 = ax4.get_position()
pos6 = ax6.get_position()


y_center = (pos2.y1 + pos6.y0) / 2

x_rigth = min(pos2.x1, pos4.x1, pos6.x1)

# Add shared y-label
fig.text(x_rigth + 0.02, y_center, 'ice fraction',
         va='center', ha='left', rotation=90)

fig.savefig('water_conservation.pdf', dpi=400, bbox_inches='tight', format='pdf')
