#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 16:49:08 2024

@author: trcn27

Coding a simple fast fourier transform so I understand how it works
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft2
import random

class Grid():
    """In the interest of doing it properly, put grid parameters in here"""
    def __init__(self):
        self.x0 = -1.0; self.x1 = 1.0
        self.y0 = -1.0; self.y1 = 1.0
        self.z0 = -0.0; self.z1 = 2.0
        
        self.nx = 192; self.ny = 192; self.nz = 128
        
        self.xs = np.linspace(self.x0, self.x1, self.nx+1)
        self.ys = np.linspace(self.y0, self.y1, self.ny+1)
        self.zs = np.linspace(self.z0, self.z1, self.nz+1)
        
        self.dx = self.xs[1] - self.xs[0]
        self.dy = self.ys[1] - self.ys[0]
        self.dz = self.zs[1] - self.zs[0]
        
        self.xc = np.linspace(self.x0 - self.dx/2, self.x1 + self.dx/2, self.nx+2)
        self.yc = np.linspace(self.y0 - self.dy/2, self.y1 + self.dy/2, self.ny+2)
        self.zc = np.linspace(self.z0 - self.dz/2, self.z1 + self.dz/2, self.nz+2)

        self.dx = self.xs[1] - self.xs[0]
        self.dy = self.ys[1] - self.ys[0]
        self.dz = self.zs[1] - self.zs[0]
        self.dzdx = self.dz/self.dx; self.dzdy = self.dz/self.dy

    def laplacian(self, phi):
        """Calculates the Laplacian of the array phi (on grid centres)"""
        lap = (phi[2:,1:-1,1:-1] + phi[:-2,1:-1,1:-1] - 2*phi[1:-1,1:-1,1:-1])/self.dx**2
        lap = lap + (phi[1:-1,2:,1:-1,] + phi[1:-1,:-2,1:-1,] - 2*phi[1:-1,1:-1,1:-1,])/self.dy**2
        lap = lap + (phi[1:-1,1:-1,2:] + phi[1:-1,1:-1,:-2] - 2*phi[1:-1,1:-1,1:-1])/self.dz**2
        return lap
    
class FT():
    """Creates various quantities related to the Fourier transform (mutliplication matrix etc.)"""
    def __init__(self,grid,lbound):
        l = np.arange(grid.nx+2)
        m = np.arange(grid.ny+2)
        j = l.reshape((grid.nx+2, 1))
        k = m.reshape((grid.ny+2, 1))
        
        self.xbasis = np.exp(-2j * np.pi * l * j / (grid.nx+2))
        self.ybasis = np.exp(-2j * np.pi * m * k / (grid.ny+2))

        #Each basis function can then be obtained with self.M[nb,:] (or the other way around, they're completely equivalent)
        self.aNMs = ifft2(lbound)
        self.d2Xdx2 = 2*np.cos(2*np.pi*l/(grid.nx+2)) - 2
        self.d2Ydy2 = 2*np.cos(2*np.pi*m/(grid.ny+2)) - 2

        #Check everything adds up with normalisation etc.
        lbound_test = np.zeros((grid.nx+2,grid.ny+2))
        for lb in range(grid.nx+2):
            for mb in range(grid.ny+2):
                lbound_test += (self.aNMs[lb,mb]*self.xbasis[lb][:,np.newaxis]*self.ybasis[mb][np.newaxis,:]).real      
        print('Fourier transform sucessful?', np.allclose(lbound_test.real, lbound))
       
    def find_zfns(self, grid, l, m):
        """Finds the basis functions in the z direction. Have to use direct integration I think, with no snazziness"""
        zfn = np.zeros((grid.nz+2))
        zfn[-1] = 1e-6; zfn[-2] = -1e-6
        for k in range(grid.nz-1,-1,-1):
            zfn[k] = 2*zfn[k+1] - zfn[k+2] - zfn[k+1]*(grid.dzdx**2)*self.d2Xdx2[l] - zfn[k+1]*(grid.dzdy**2)*self.d2Ydy2[m]
        zfn[:] = zfn[:]*grid.dz/(zfn[1] - zfn[0])
        return zfn
    
    def find_phi(self,grid,mode_cutoff=1e-6):
        phi = np.zeros((grid.nx+2,grid.ny+2,grid.nz+2))
        for l in range(grid.nx+2):
            for m in range(grid.ny+2):
                if np.abs(self.aNMs[l,m]) > mode_cutoff:
                    print(l, m, np.abs(self.aNMs[l,m]))
                    phi = phi + (self.aNMs[l,m]*self.xbasis[l,:][:,np.newaxis,np.newaxis]*self.ybasis[m,:][np.newaxis,:,np.newaxis]).real*self.find_zfns(grid,l,m)[:][np.newaxis,np.newaxis,:]
        lbound_test = (phi[:,:,1]- phi[:,:,0])/grid.dz
        plt.pcolormesh(lbound_test.T)
        plt.show()
        print('Lbound match?', np.allclose(lbound,lbound_test))
        print('Max. Lbound error =', np.max(np.abs(lbound-lbound_test)))
        print('Max Laplacian =', np.max(np.abs(grid.laplacian(phi))))
        return phi
    
    
def generate_lbound(grid):
    """Generates the lower boundary function. To be normalised such that this is just the X basis"""
    sf = grid.x1/12.0

    dipole_mag = 25.0; zstar = grid.z1*1.5/24.0

    lbound = 0.0*grid.xc[:,np.newaxis]*grid.yc[np.newaxis,:]
    for ic, ix in enumerate(grid.xc[:]):
        for jc, jy in enumerate(grid.yc[:]):
            lbound[ic,jc] = sf**3*dipole_mag*(2*(zstar)**2 - ((ix)**2 + (jy)**2))/(((ix)**2 + (jy)**2 + (zstar)**2)**2.5)
        
    #Enforce boundary conditions
    lbound[-1,:] = lbound[-2,:]
    lbound[0,:]  = lbound[1,:]
    lbound[:,-1] = lbound[:,-2]
    lbound[:,0]  = lbound[:,1]

    
    return lbound
    
grid = Grid()
lbound = generate_lbound(grid) 
ft = FT(grid, lbound)
phi = ft.find_phi(grid,mode_cutoff=1e-10)
