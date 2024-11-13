#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: trcn27
"""

import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import matplotlib
from scipy.io import netcdf_file
from fltrace import trace_fieldlines



class Grid():
    def __init__(self):
        self.x0 = xs[0]; self.x1 = xs[-1]
        self.y0 = ys[0]; self.y1 = ys[-1]
        self.z0 = zs[0]; self.z1 = zs[-1]
        self.nx = nx ; self.ny = ny; self.nz = nz

i = 0

fig, axs = plt.subplots(2,2, figsize = (10,7))

nruns = 6
nsnaps = 150
do_nulls = True

angles = [0.0,0.01,0.05,0.1,0.5,1.0,5.0,10.0]
cs = ['red','blue','green','orange','brown', 'grey','purple','black']

for pressure in [0,1]:
    for run in [0,1,2,3,4,5,6,7]:
        run = pressure*10 + run
        backfield_angle = angles[run%10] #Angle of background field in degrees

        fname = './diagnostics/run%02d.nc' % (run)

        try:
            data = netcdf_file(fname, 'r', mmap=False)
            paras = np.loadtxt('parameters/variables%03d.txt' % run)
            print('Diagnostics from run', run)
        except:
            continue

        nx = int(paras[1])
        ny = int(paras[2])
        nz = int(paras[3])

        nsnaps = int(paras[5])

        remote_flag = int(paras[22])

        xs = np.linspace(-12,12, nx+1)
        ys = np.linspace(-12,12, ny+1)
        zs = np.linspace(0,24, nz+1)

        xc = 0.5*(xs[1:] + xs[:-1])
        yc = 0.5*(ys[1:] + ys[:-1])
        zc = 0.5*(zs[1:] + zs[:-1])

        dx = xs[1] - xs[0]
        dy = ys[1] - ys[0]
        dz = zs[1] - zs[0]

        xc = 0.5*(xs[1:] + xs[:-1])
        yc = 0.5*(ys[1:] + ys[:-1])
        zc = 0.5*(zs[1:] + zs[:-1])

        dx = xs[1] - xs[0]
        dy = ys[1] - ys[0]
        dz = zs[1] - zs[0]

        if do_nulls:
            nulls = []
            nullts = []
            def null_height(bz):
                #Determine(roughly) the height of the null point, up to the eruption
                zslice = bz[slice_index,slice_index,:]
                flip_index = np.where(np.sign(zslice[1:])*np.sign(zslice[:-1]) < 0)
                return 0.5*(zs[flip_index[-1]] + zs[flip_index[-1]+1])

            #Calculate null height to have a look at eruption. Can't do this in diagnostics because dim
            if remote_flag:
                data_directory = '/nobackup/trcn27/mf3d0/%03d/' % run
            else:
                data_directory = '/extra/tmp/trcn27/mf3d/%03d/' % run

            slice_index = ny//2

            for snap in range(nsnaps):
                wait = 0
                fname = '%s%04d.nc' % (data_directory, snap)
                try:
                    data = netcdf_file(fname, 'r', mmap=False)
                except:
                    break
                bz = np.zeros((nx+2,ny+2,nz+1))
                bz[1:-1,1:-1,:] = np.swapaxes(data.variables['bz'][:],0,2)
                null = null_height(bz)
                if len(null) > 0:
                    nulls.append(null[-1])
                else:
                    nulls.append(np.nan)
                nullts.append(750*snap/150)

        fname = './diagnostics/run%02d.nc' % (run)

        data = netcdf_file(fname, 'r', mmap=False)

        time = data.variables['time'][:]
        openflux = data.variables['openflux'][:]
        avgcurrent = data.variables['avgcurrent'][:]
        energy = data.variables['energy'][:]

        end = np.sum(energy > 0.)

        if pressure == 0:
            im = axs[0,0].plot(time[:end], openflux[:end], label = 'Background tilt = %.2f' % (backfield_angle), color = cs[run%10])
            axs[0,0].set_title('Open Flux')
            im = axs[0,1].plot(time[:end], avgcurrent[:end], color = cs[run%10])
            axs[0,1].set_title('Avg. Current')

            im = axs[1,0].plot(time[:end], energy[:end], color = cs[run%10])
            axs[1,0].set_title('Magnetic Energy')

            if do_nulls:
                im = axs[1,1].plot(nullts, nulls, color = cs[run%10])
                axs[1,1].set_title('Null Height')
                plt.ylim(ymin = 0.0)

        else:
            im = axs[0,0].plot(time[:end], openflux[:end], color = cs[run%10], linestyle = 'dashed')
            axs[0,0].set_title('Open Flux')
            im = axs[0,1].plot(time[:end], avgcurrent[:end], color = cs[run%10], linestyle = 'dashed')
            axs[0,1].set_title('Avg. Current')

            im = axs[1,0].plot(time[:end], energy[:end], color = cs[run%10], linestyle = 'dashed')
            axs[1,0].set_title('Magnetic Energy')

            if do_nulls:
                im = axs[1,1].plot(nullts, nulls, color = cs[run%10], linestyle = 'dashed')
                axs[1,1].set_title('Null Height')
                plt.ylim(ymin = 0.0)


axs[0,0].legend()
plt.savefig('backtilt_nopressure.png')
plt.show()














