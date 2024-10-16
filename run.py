#Shell script to run the 2D magnetofrictional code. Initially NOT using MPI as it's hard and should be quick enough without. Just a nice test of my abilities.

import os
import shutil
import numpy as np
import sys
from numpy import random
import time
import matplotlib.pyplot as plt

from init import compute_initial_condition
from fltrace import trace_fieldlines

if len(sys.argv) > 1:
    run = int(sys.argv[1])
else:
    run = 0

if len(sys.argv) > 2:
    nprocs = int(sys.argv[2])
else:
    nprocs = 1

hflag = 0
#DYNAMIC SYSTEM PARAMETERS
#-------------------------------------
voutfact = 0.0
shearfact = 5e-4
eta0 = 0.0

tmax = 600.0

nx = 64
ny = 64
nz = 64

nplots = 600
ndiags = 10

nu0 = 1.0
eta = 1e-4

x0 = -12.0; x1 = 12.0
y0 = -12.0; y1 = 12.0
z0 = -0.0; z1 = 24.0

#Variables for the pressure term
a = 0.0; b = 1.0; deltaz = 0.05; zstar = 0.0

#INITIAL LOWER BOUNDARY CONDITION
#-------------------------------------

def lbound_fn(x,y):
    #Outputs lower boundary radial magnetic field as a function of position x
    lbound_fn = np.zeros((len(x),len(y)))
    #lbound_fn[:,:] = -np.sin(3*x[:,np.newaxis])*np.sin(y[np.newaxis,:]) + 0.01*random.rand(len(x),len(y))
    d = 0.25; z0 = 0.5
    for i, xi in enumerate(x):
        for j,yj in enumerate(y):
            lbound_fn[i,j] = z0/(((xi-d)**2 + yj**2 + z0**2)**1.5) - z0/(((xi+d)**2 + yj**2 + z0**2)**1.5)

    return lbound_fn
    #return random.rand(len(x))

def lbound_pariat(x,y):
    #Outputs lower boundary radial magnetic field as a function of position x
    lbound_fn = np.zeros((len(x),len(y)))
    sf = x1/12.0

    dipole_mag = 25.0; zstar = z1*1.5/24.0
    for i, xi in enumerate(x):
        for j,yj in enumerate(y):
            lbound_fn[i,j] = sf**3*dipole_mag*(2*(zstar)**2 - ((xi)**2 + (yj)**2))/(((xi)**2 + (yj)**2 + (zstar)**2)**2.5)

    print('Max lbound', np.max(lbound_fn))
    print('Lbound flux', np.sum(np.abs(lbound_fn)))
    #plt.pcolormesh(lbound_fn)
    #plt.show()
    return lbound_fn

def surface_flow(grid, bz_lbound):
    br = 13; bl = 0.1; kb = 15.0
    #Finds the velocity field for the Pariat jet simulations, using the contours of bz.
    bzdy = (bz_lbound[:,1:] - bz_lbound[:,:-1])/grid.dy
    bzdx = (bz_lbound[1:,:] - bz_lbound[:-1,:])/grid.dx

    bzdy = 0.5*(bzdy[1:,:] + bzdy[:-1,:])
    bzdx = 0.5*(bzdx[:,1:] + bzdx[:,:-1])

    fact = kb*(br-bl)/bz_lbound*np.tanh(kb*(bz_lbound - bl*np.ones((nx+2,ny+2)))/(br-bl))

    fact = 0.25*(fact[:-1,:-1] + fact[:-1,1:] + fact[1:,:-1] + fact[:-1,:-1])

    return -fact*bzdy,fact*bzdx

#SAVE VARIABLES TO BE READ BY FORTRAN
#-------------------------------------

variables = np.zeros((30))
#Basics
variables[0] = run
variables[1] = nx
variables[2] = ny
variables[3] = nz
variables[4] = tmax
#Outputs
variables[5] = nplots
variables[6] = ndiags
#Parameters
variables[7] = voutfact
variables[8] = shearfact
variables[9] = eta
variables[10] = nu0
variables[11] = eta0
#Grid things
variables[12] = x0
variables[13] = x1
variables[14] = y0
variables[15] = y1
variables[16] = z0
variables[17] = z1

#Pressure things
variables[18] = a
variables[19] = b
variables[20] = deltaz
variables[21] = zstar

variables[22] = hflag

np.savetxt('parameters/variables%03d.txt' % run, variables)   #variables numbered based on run number (up to 1000)

#SOME FOLDER ADMIN
#-------------------------------------


data_directory = '/extra/tmp/trcn27/mf3d/%03d/' % run

if os.path.isdir(data_directory):
    for i in range(1000):
        if os.path.isfile('%s%04d.nc' % (data_directory, i)):
            os.remove('%s%04d.nc' % (data_directory, i))
else:
    os.mkdir(data_directory)


#Create initial condition using new init.py (potential field with arbitrary lower boundary and domain dimensions)
#-------------------------------------

class Grid():
    def __init__(self):
        self.x0 = x0; self.x1 = x1
        self.y0 = y0; self.y1 = y1
        self.z0 = z0; self.z1 = z1
        self.nx = nx ; self.ny = ny; self.nz = nz
        self.dx = (self.x1 - self.x0)/nx
        self.dy = (self.y1 - self.y0)/ny
        self.dz = (self.z1 - self.z0)/nz
        self.xs = np.linspace(self.x0,self.x1,self.nx+1)
        self.ys = np.linspace(self.y0,self.y1,self.ny+1)
        self.zs = np.linspace(self.z0,self.z1,self.nz+1)

if True:
    grid= Grid()
    init = compute_initial_condition(grid, lbound_pariat, run, background_strength = 1.0, boundary_error_limit = 1e-6, init_filename = './inits/init%03d.nc' % run)

    bx = init.bx; by = init.by; bz = init.bz

    vx_surf, vy_surf = surface_flow(grid,bz[:,:,0])

    if False:
        plt.pcolormesh(grid.xs,grid.ys,bz[1:-1,1:-1,0].T,cmap='plasma')
        plt.colorbar()

        plt.quiver(grid.xs,grid.ys,vx_surf.T,vy_surf.T)
        plt.show()
    #trace_fieldlines(Grid(),bx,by,bz)

#RUN CODE
#-------------------------------------

if True:
    os.system('make')
    print('Using output directory "%s"' % (data_directory))
    if nprocs <= 4:
        os.system('/usr/lib64/openmpi/bin/mpiexec -np %d ./bin/mf3d %d' % (nprocs, run))
    else:
        os.system('/usr/lib64/openmpi/bin/mpiexec -np %d --oversubscribe ./bin/mf3d %d' % (nprocs, run))
#os.system('python diagnostics.py %d' % run)

#os.system('rm parameters/variables%03d.txt' % run)





























