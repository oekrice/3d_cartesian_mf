#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 16:49:08 2024

@author: trcn27

Coding a simple fast fourier transform so I understand how it works
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft2, fft2, ifft
import random

class Grid():
    """In the interest of doing it properly, put grid parameters in here"""
    def __init__(self):
        self.x0 = -np.pi; self.x1 = np.pi
        self.y0 = -1.0; self.y1 = 1.0
        
        self.nx = 128; self.ny = 128
        
        self.xs = np.linspace(self.x0, self.x1, self.nx+1)
        self.ys = np.linspace(self.y0, self.y1, self.ny+1)
        
        self.dx = self.xs[1] - self.xs[0]
        self.dy = self.ys[1] - self.ys[0]
        
        self.xc = np.linspace(self.x0 - self.dx/2, self.x1 + self.dx/2, self.nx+2)
        self.yc = np.linspace(self.y0 - self.dy/2, self.y1 + self.dy/2, self.ny+2)

        self.dx = self.xs[1] - self.xs[0]
        self.dy = self.ys[1] - self.ys[0]
        
    def lap_xribs(self, xribs):
        """Calculates laplacian of a quantity saved on the x ribs, using staggered as appropriate"""
        lap_x = np.zeros((self.nx+2, self.ny+1))
        #x direction
        lap_x[1:-1,1:-1] += (xribs[2:,1:-1] -  2*xribs[1:-1,1:-1] + xribs[:-2,1:-1])/self.dx**2
        lap_x[1:-1,1:-1] += (xribs[1:-1,2:] -  2*xribs[1:-1,1:-1] + xribs[1:-1,:-2])/self.dy**2
        return lap_x

    def lap_yribs(self, yribs):
        """Calculates laplacian of a quantity saved on the y ribs, using staggered as appropriate"""
        lap_y = np.zeros((self.nx+1, self.ny+2))
        #x direction
        lap_y[1:-1,1:-1] += (yribs[2:,1:-1] -  2*yribs[1:-1,1:-1] + yribs[:-2,1:-1])/self.dx**2
        lap_y[1:-1,1:-1] += (yribs[1:-1,2:] -  2*yribs[1:-1,1:-1] + yribs[1:-1,:-2])/self.dy**2
        return lap_y


    def curl_inplane(self, C):
        """Calculates the in-plane curl of the quantity C(x,y)"""
        """Outputs on the respective ribs, hopefully"""
        
        curl_x = np.zeros((self.nx+2, self.ny+1))
        curl_y = np.zeros((self.nx+1, self.ny+2))

        curl_x = (C[:,1:] - C[:,:-1])/self.dy
        curl_y = -(C[1:,:] - C[:-1,:])/self.dx
        
        return curl_x, curl_y

    def curl_E(self, E_x, E_y):
        """Returns the in-plane curl of the vectors E"""
        curl = np.zeros((grid.nx+2, grid.ny+2))
        curl[1:-1,1:-1] += (E_x[1:-1,1:] - E_x[1:-1,:-1])/self.dy
        curl[1:-1,1:-1] -= (E_y[1:,1:-1] - E_y[:-1,1:-1])/self.dx

        return curl
        

    def lap_points(self, points):
        """Calculates laplacian of a quantity saved on the y ribs, using staggered as appropriate"""
        lap_p = np.zeros((self.nx+1, self.ny+1))
        #x direction
        lap_p[1:-1,1:-1] += (points[2:,1:-1] -  2*points[1:-1,1:-1] + points[:-2,1:-1])/self.dx**2
        lap_p[1:-1,1:-1] += (points[1:-1,2:] -  2*points[1:-1,1:-1] + points[1:-1,:-2])/self.dy**2
        return lap_p

    def lap_centres(self, points):
        """Calculates laplacian of a quantity saved on the y ribs, using staggered as appropriate"""
        lap_p = np.zeros((self.nx+2, self.ny+2))
        #x direction
        lap_p[1:-1,1:-1] += (points[2:,1:-1] -  2*points[1:-1,1:-1] + points[:-2,1:-1])/self.dx**2
        lap_p[1:-1,1:-1] += (points[1:-1,2:] -  2*points[1:-1,1:-1] + points[1:-1,:-2])/self.dy**2
        return lap_p

    def lap_points_eh(self, points):
        """Calculates laplacian of a quantity saved on the y ribs, using staggered as appropriate"""
        
        grad_x = np.zeros((self.nx+2, self.ny+1))
        grad_y = np.zeros((self.nx+1, self.ny+2))
        grad_x[1:-1,:] = (points[1:,:] - points[:-1,:])/self.dx
        grad_y[:,1:-1] = (points[:,1:] - points[:,:-1])/self.dy

        div = np.zeros((grid.nx+1, grid.ny+1))
        div[:,:] += (grad_x[1:,:] - grad_x[:-1,:])/self.dx
        div[:,:] += (grad_y[:,1:] - grad_y[:,:-1])/self.dy

        lap_p = np.zeros((self.nx+1, self.ny+1))
        #x direction
        lap_p[1:-1,1:-1] += (points[2:,1:-1] -  2*points[1:-1,1:-1] + points[:-2,1:-1])/self.dx**2
        lap_p[1:-1,1:-1] += (points[1:-1,2:] -  2*points[1:-1,1:-1] + points[1:-1,:-2])/self.dy**2
        return div


    def div_E(self, E_x, E_y):
        """Returns the in-plane curl of the vectors E"""
        div = np.zeros((grid.nx+1, grid.ny+1))
        div[:,:] += (E_x[1:,:] - E_x[:-1,:])/self.dx
        div[:,:] += (E_y[:,1:] - E_y[:,:-1])/self.dy

        return div
    
    def grad(self, phi):
        """Returns gradients of phi (at grid points)"""
        grad_x = np.zeros((self.nx+2, self.ny+1))
        grad_y = np.zeros((self.nx+1, self.ny+2))
        grad_x[1:-1,:] = (phi[1:,:] - phi[:-1,:])/self.dx
        grad_y[:,1:-1] = (phi[:,1:] - phi[:,:-1])/self.dy
        
        return grad_x, grad_y
        
        
class FT():
    
    """Creates various quantities related to the Fourier transform (mutliplication matrix etc.)"""
    def __init__(self,grid):
        pass
        
    def test_1d(self, grid):
        """1D test to find periodic solutions with given RHS double derivative. Using x direction"""
        rhs = np.sin(grid.xc*6)         #Double x derivative
        rhs_transform = ifft(rhs)  #k twiddle
        j = np.arange(grid.nx+2)
        xbasis = np.exp(-2j * np.pi * j / (grid.nx+2))

        d2Xdx2 = (2*np.cos(2*np.pi*j/(grid.nx+2)) - 2)/grid.dx**2
        a_twiddle = xbasis*0.
        a_twiddle[1:] = rhs_transform[1:]/d2Xdx2[1:]
        
        E = fft(a_twiddle)
            
        d2EdE2 = (E[2:] + E[:-2] - 2*E[1:-1])/grid.dx**2
        plt.plot(rhs[2:-2])
        #plt.plot(E[2:-2])
        plt.plot(d2EdE2[1:-1])
        plt.show()
        
    def x_transform(self, rhs_x):
        """2D test in the x rib direction"""
        rhs_transform = ifft2(rhs_x)  
        
        j = np.arange(grid.nx+2)
        k = np.arange(grid.ny+1)

        d2Xdx2 = (2*np.cos(2*np.pi*j/(grid.nx+2)) - 2)/grid.dx**2
        d2Ydy2 = (2*np.cos(2*np.pi*k/(grid.ny+1)) - 2)/grid.dy**2
        
        d2 = d2Xdx2.reshape((grid.nx+2,1)) + d2Ydy2
        
        d2[0,0] = 1.0   #Avoid problems with the zero value. All this does is add a constant anyway so set to zero in the transform
        
        a_twiddle = rhs_transform/d2
        
        a_twiddle[0,0] = 0.0
        
        E_x = (fft2(a_twiddle)).real
        
        return E_x
           
    def y_transform(self, rhs_y):
        """2D test in the y rib direction"""
        rhs_transform = ifft2(rhs_y)  
        
        j = np.arange(grid.nx+1)
        k = np.arange(grid.ny+2)

        d2Xdx2 = (2*np.cos(2*np.pi*j/(grid.nx+1)) - 2)/grid.dx**2
        d2Ydy2 = (2*np.cos(2*np.pi*k/(grid.ny+2)) - 2)/grid.dy**2
        
        d2 = d2Xdx2.reshape((grid.nx+1,1)) + d2Ydy2
        
        d2[0,0] = 1.0   #Avoid problems with the zero value. All this does is add a constant anyway so set to zero in the transform
        
        a_twiddle = rhs_transform/d2
        
        a_twiddle[0,0] = 0.0
        
        E_y = (fft2(a_twiddle)).real
        
        return E_y

    def point_transform(self, rhs_points):
        """2D transform in at the grid points"""
        rhs_transform = ifft2(rhs_points)  
        
        j = np.arange(grid.nx+1)
        k = np.arange(grid.ny+1)

        d2Xdx2 = (2*np.cos(2*np.pi*j/(grid.nx+1)) - 2)/grid.dx**2
        d2Ydy2 = (2*np.cos(2*np.pi*k/(grid.ny+1)) - 2)/grid.dy**2
        
        d2 = d2Xdx2.reshape((grid.nx+1,1)) + d2Ydy2
        
        d2[0,0] = 1.0   #Avoid problems with the zero value. All this does is add a constant anyway so set to zero in the transform
        
        a_twiddle = rhs_transform/d2
        
        a_twiddle[0,0] = 0.0
        
        phi = (fft2(a_twiddle)).real
        
        return phi

    def centre_transform(self, rhs_centres):
        """2D transform in at the grid centres"""
        rhs_transform = ifft2(rhs_centres)  
        
        j = np.arange(grid.nx+2)
        k = np.arange(grid.ny+2)

        d2Xdx2 = (2*np.cos(2*np.pi*j/(grid.nx+2)) - 2)/grid.dx**2
        d2Ydy2 = (2*np.cos(2*np.pi*k/(grid.ny+2)) - 2)/grid.dy**2
        
        d2 = d2Xdx2.reshape((grid.nx+2,1)) + d2Ydy2
        
        d2[0,0] = 1.0   #Avoid problems with the zero value. All this does is add a constant anyway so set to zero in the transform
        
        a_twiddle = rhs_transform/d2
        
        a_twiddle[0,0] = 0.0
        
        G = (fft2(a_twiddle)).real
        
        return G

    def test_2d(self, grid, e_x):
        """2D test in the x rib direction"""
        rhs = e_x#np.sin(grid.xc[:,np.newaxis]**2*grid.ys[np.newaxis,:] + 2.0)         #Double x derivative
        rhs_transform = ifft2(rhs)  #k twiddle
        
        j = np.arange(grid.nx+2)
        k = np.arange(grid.ny+1)

        xbasis = np.exp(-2j * np.pi * j / (grid.nx+2))
        ybasis = np.exp(-2j * np.pi * k / (grid.ny+1))

        d2Xdx2 = (2*np.cos(2*np.pi*j/(grid.nx+2)) - 2)/grid.dx**2
        d2Ydy2 = (2*np.cos(2*np.pi*k/(grid.ny+1)) - 2)/grid.dy**2
        
        d2 = d2Xdx2.reshape((grid.nx+2,1)) + d2Ydy2
        
        d2[0,0] = 1.0   #Avoid problems with the zero value. All this does is add a constant anyway so set to zero in the transform
        
        a_twiddle = rhs_transform*0.
        a_twiddle = rhs_transform[:]/d2
        
        a_twiddle[0,0] = rhs_transform[0,0]
        
        E = (fft2(a_twiddle)).real
        
        Ediff = grid.lap_xribs(E)
    
        
        print('ediffsum', np.sum(Ediff[1:-1,1:-1]), np.sum(rhs[1:-1,1:-1]))
        #plt.plot(E[2:-2])
        fig, ax = plt.subplots(3,1, figsize = (10,7))
        ax[0].pcolormesh(Ediff[1:-1,1:-1])
        ax[1].pcolormesh(rhs[1:-1,1:-1])
        im = ax[2].pcolormesh(Ediff[1:-1,1:-1] - rhs[1:-1,1:-1])
        
        fig.colorbar(im)
        
        plt.tight_layout()
        plt.show()
        
        diff = np.max(np.abs(Ediff[1:-1,1:-1] - rhs[1:-1,1:-1]))
        print('Zero test', diff, np.sum(Ediff[1:-1,1:-1]), np.sum(rhs[1:-1,1:-1]))

    
def generate_lbound(grid, offset = 0):
    """Generates the lower boundary function. To be normalised such that this is just the X basis"""
    sf = grid.x1/12.0

    dipole_mag = 25.0; zstar = 2.*1.5/24.0

    lbound = 0.0*grid.xc[:,np.newaxis]*grid.yc[np.newaxis,:]
    for ic, ix in enumerate(grid.xc[:]):
        for jc, jy in enumerate(grid.yc[:] - offset):
            lbound[ic,jc] = sf**3*dipole_mag*(2*(zstar)**2 - ((ix)**2 + (jy)**2))/(((ix)**2 + (jy)**2 + (zstar)**2)**2.5)

    #Enforce boundary conditions
    lbound[-1,:] = lbound[-2,:]
    lbound[0,:]  = lbound[1,:]
    lbound[:,-1] = lbound[:,-2]
    lbound[:,0]  = lbound[:,1]

    #Remove net flux

    lbound -= np.sum(lbound)/np.size(lbound)

    plt.pcolormesh(lbound.T)
    plt.show()

    return lbound
    
grid = Grid()
lbound1 = generate_lbound(grid, offset = 0.0) 
lbound2 = generate_lbound(grid, offset = 0.2) 

C = lbound2 - lbound1    #Difference in the boundaries -- this is the important RHS bit
#curl_x, curl_y = grid.curl_inplane(C)
#Assume zero in-plane divergence for now, for E? Can add stuff easily later I hope. 

ft = FT(grid)
G = ft.centre_transform(-C)
e_x, e_y = grid.curl_inplane(G)
curl_test = grid.curl_E(e_x, e_y)

#e_x = ft.x_transform(-curl_x)
#e_y = ft.y_transform(-curl_y)

#lap_x = grid.lap_xribs(e_x)
#lap_y = grid.lap_yribs(e_y)


print('Curl test', np.max(np.abs(curl_test[1:-1,1:-1] + C[1:-1,1:-1])))
#Determine the curl-free part of the electric field and remove it (for now -- can specify it later instead)
div_test = grid.div_E(e_x, e_y)
phi = ft.point_transform(-div_test)
correct_x, correct_y = grid.grad(phi)
e_x += correct_x
e_y += correct_y
div_test = grid.div_E(e_x, e_y)
print('Div Test', np.max(np.abs(div_test[1:-1,1:-1])))

curl_test = grid.curl_E(e_x, e_y)

print('Overall', np.max(np.abs(curl_test[1:-1,1:-1] + C[1:-1,1:-1])))
#print('xs', np.max(np.abs(lap_x[1:-1,1:-1]+curl_x[1:-1,1:-1])))
#print('ys', np.max(np.abs(lap_y[1:-1,1:-1]+curl_y[1:-1,1:-1])))
#ft.test_1d(grid)
#ft.test_2d(grid, curl_x)
#phi = ft.find_phi(grid,mode_cutoff=1e-10)





















