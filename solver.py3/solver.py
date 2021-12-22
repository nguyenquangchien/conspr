#!/usr/bin/env python
""" Module solver.py
Class Flow and Domain
Flow solver for a two-dimensional domain
on a uniform rectangular grid
using finite volume method
Written by Nguyen Quang Chien, 2009

Revised 2021

Numerical method based on XBeach, http://xbeach.org/
License GNU GPL, http://www.gnu.org/copyleft/gpl.html """

import numpy as np
from copy import deepcopy

class Param(object):
    """ Parameters related to the flow """
    def __init__(self):
        # hydraulic parameters
        self.Ch = 65        # Chezy coefficient, m^0.5 s^-1 
        self.rho = 1000     # water density, kg m^-3
        self.eps = 0.1      # diffusivity, m^2 s-1
        self.g = 9.81       # acceleration of gravity, m^2 s^-1
        self.hmin = 0.001   # minimum water depth, m
        self.umin = 0.1     # minimum flow velocity, m s^-1
        self.CFL = 0.4      # Courant-Fredricht-Lewy number
        self.zs0 = 0        # initial water level, m
        self.v0 = 1.0       # initial velocity, m s^-1
        self.i = 1.E-4      # longitudinal slope
        
        # timing parameters
        self.tstart = 0     # start time of simulation, s
        self.t = self.tstart
        self.tint = 1       # time interval, s
        self.tnext = 0 + self.tint
        self.tstop = 100    # end time of simulation, s
        
        # contaminant parameters
        self.k = 0          # decay rate, s^-1
        self.Elong = 0      # dispersion rate in longitudinal direction, m^2 s^-1
        self.Etran = 0      # dispersion rate in transversal direction, m^2 s^-1


class Domain(object):
    def __init__(self, flowpar, type='river', ext=None):
        self.makeBathy(flowpar, type, 4, 0.01, 0.5)
        
        # Forcing etc.
        domainDim = np.shape(self.zb)
        self.Fx = np.zeros(domainDim)
        self.Fy = np.zeros(domainDim)
        self.W = np.zeros(domainDim)
        
        # Init flow condition
        zs0 = flowpar.zs0
        zs = np.zeros(domainDim) - flowpar.i * self.y
        self.zs = np.maximum(zs, self.zb)
        self.h = self.zs - self.zb
        self.dzsdt = np.zeros(domainDim)
        self.dzbdt = np.zeros(domainDim)
        self.uu = np.zeros(domainDim)
        self.vv = np.zeros(domainDim)
        self.vv = flowpar.Ch * np.sqrt(self.h * flowpar.i)
        self.vv = self.vv * (self.h > 0)
        self.qx = np.zeros(domainDim)
        self.qy = np.zeros(domainDim)
        self.c = np.zeros(domainDim)


    def makeBathy(self, flowpar, kind, *arg):
        """ Creates the bathymetry for the model.

            Parameters
            ----------
            `flowpar`: instance of `Param`
                parameters of the flow

            `kind`: `str`
                the kind of bathymetry
                valid values include `'river'` and
                `'San Francisco bay'`
            
            `*arg`: `tuple`
                providing additional info in case `kind == 'river'`
                * `h0` the maximum river depth
                * `slope_t` the transversal slope
                * `bank_height` the river bank height
        """
        if kind.lower() == 'river':
            if not 'ext' in locals():
                ext = (40, 40, 10, 40)
            (self.nx, self.ny, self.dx, self.dy) = ext
            domainDim = (self.nx + 1, self.ny + 1)
            self.x = np.zeros(domainDim)
            self.y = np.zeros(domainDim)
            self.zb = np.zeros(domainDim)
            h0 = arg[0]         # max. river depth
            slope_t = arg[1]    # transversal slope
            bank_height = arg[2]     # river bank
            for j in range(self.ny + 1):
                for i in range(self.nx + 1):
                    self.x[i,j] = i * self.dx
                    self.y[i,j] = j * self.dy
                    self.zb[i,j] = -h0 + self.x[i,j] * slope_t
            self.zb[30:, :] = bank_height 
            self.zb -= flowpar.i * self.y     # Longitudinal slope i
        elif kind.lower() == 'san francisco bay':
            self.nx = self.ny = 199
            self.dx = 100
            self.dy = 200
            
            self.x, self.y = np.meshgrid(
                    np.arange(0, 20000, self.dx), np.arange(0, 40000, self.dy))
            self.zb = np.fromfile('newgrid.txt', sep=' ').reshape(self.x.shape)
            print(self.x.shape, self.y.shape, self.zb.shape)


    def solve(self, flowpar):
        x = deepcopy(self.x)
        y = deepcopy(self.y)
        dx = deepcopy(self.dx)
        dy = deepcopy(self.dy)
        uu  = deepcopy(self.uu)
        vv = deepcopy(self.vv)
        qx = deepcopy(self.qx)
        qy = deepcopy(self.qy)
        zb = deepcopy(self.zb)
        zs = deepcopy(self.zs)
        Fx = deepcopy(self.Fx)
        Fy = deepcopy(self.Fy)
        c = deepcopy(self.c)
        W = deepcopy(self.W)

        g = flowpar.g
        Ch = flowpar.Ch
        rho = flowpar.rho
        umin = flowpar.umin
        hmin = flowpar.hmin
        eps = flowpar.eps
        t = flowpar.t
        tint = flowpar.tint
        tnext = flowpar.tnext
        CFL = flowpar.CFL
        Elong = flowpar.Elong
        Etran = flowpar.Etran

        domainDim = np.shape(x)
        vu      = np.zeros(domainDim)
        vsu     = np.zeros(domainDim)
        usu     = np.zeros(domainDim)
        vsv     = np.zeros(domainDim)
        usv     = np.zeros(domainDim)
        uv      = np.zeros(domainDim)
        vmagu   = np.zeros(domainDim)
        vmague  = np.zeros(domainDim)
        vmagv   = np.zeros(domainDim)
        vmagve  = np.zeros(domainDim)
        veu     = np.zeros(domainDim)
        ueu     = np.zeros(domainDim)
        vev     = np.zeros(domainDim)
        uev     = np.zeros(domainDim)

        u       = np.zeros(domainDim)
        v       = np.zeros(domainDim)
        dzsdx   = np.zeros(domainDim)
        dzsdy   = np.zeros(domainDim)
        dzsdt   = np.zeros(domainDim)
        dcdt    = np.zeros(domainDim)
        mx      = np.zeros(domainDim)
        mx1     = np.zeros(domainDim)
        mx2     = np.zeros(domainDim)
        my      = np.zeros(domainDim)
        my1     = np.zeros(domainDim)
        my2     = np.zeros(domainDim)
        f1      = np.zeros(domainDim)
        f2      = np.zeros(domainDim)
        ududx   = np.zeros(domainDim)
        vdudy   = np.zeros(domainDim)
        udvdx   = np.zeros(domainDim)
        vdvdy   = np.zeros(domainDim)
        qxm     = np.zeros(domainDim)
        qym     = np.zeros(domainDim)
        us      = np.zeros(domainDim)
        vs      = np.zeros(domainDim)
        ue      = np.zeros(domainDim)
        ve      = np.zeros(domainDim)
        Ecx     = np.zeros(domainDim)
        Ecy     = np.zeros(domainDim)

        # V-velocities @ u-points
        vu[:-1, 1:-1] = 0.25 * (vv[:-1, :-2] + vv[:-1, 1:-1] +
								vv[1:, :-2] + vv[1:, 1:-1])
        # at boundaries
        vu[:,0] = vu[:,1]
        vu[:, -1] = vu[:, -2]
        
        # V-stokes velocities at U point
        vsu[:-2, 1:-2] = 0
        vsu[:, 0] = vsu[:, 1]
        vsu[:, -1] = vsu[:, -2]
        # U-stokes velocities at U point
        usu[:-2, 1:-2] = 0
        usu[:, 0] = usu[:, 1]
        usu[:, -1] = usu[:, -2]
        
        veu = vu - vsu                  # V-euler velocities @ u-point
        ueu = uu - usu                  # U-euler velocties @ u-point
        vmagu = np.sqrt(uu**2 + vu**2)  # velocity magnitude
        vmageu = np.sqrt(ueu**2 + veu**2)    # Eulerian velocity magnitude

        # U-velocities at v-points
        uv[1:-1, :-1] = 0.25 * (uu[:-2, :-1] + uu[1:-1, :-1] +
                				uu[:-2, 1:] + uu[1:-1, 1:])

        uv[:, -1] = uv[:, -2]   # at the boundary

        # V-stokes velocities at V point
        vsv[1:-1, :-1] = 0
        vsv[:, 0] = vsv[:, 1]
        vsv[:, -1] = vsv[:, -2]
        # U-stokes velocities at V point
        usv[1:-1, :-1] = 0
        usv[:, 0] = usv[:, 1]
        usv[:, -1] = usv[:, -2]

        vev = vv - vsv  # V-euler velocities at V-point
        uev = uv - usv  # U-euler velocties at V-point
        vmagv = np.sqrt(uv**2 + vv**2)  # Velocity magnitude
        vmagev = np.sqrt(uev**2 + vev**2)   # Eulerian velocity magnitude

        # Water level slopes
        dzsdx[:-1, :] = (zs[1:, :] - zs[:-1, :]) / dx
        dzsdy[:, :-1] = (zs[:, 1:] - zs[:, :-1]) / dy

        # Upwind method implemented through weight factors mx1 and mx2
        
        h = np.maximum(zs - zb, hmin)    # water depth
        weth = h > hmin
        # X-direction
        mx = np.sign(uu)
        mx[np.abs(uu) < umin] = 0 
        mx1 = (mx + 1) / 2
        mx2 = -(mx - 1) / 2
        f1[:-1, :] = h[:-1, :]
        f2[:-1, :] = h[1:, :]
        
        hu = mx1 * f1 + mx2 * f2    # Water depth @ u-points for cont. eq.: upwind
        hum = np.max(0.5 * (f1 + f2), hmin) # Water depth @ u-points for mom. eq.: mean
        # Advection terms (momentum conserving method)
        f1.fill(0)
        f2.fill(0)
        f1[1:-1, :] = (qx[1:-1, :] + qx[:-2,:]) * (uu[1:-1, :] - uu[:-2,:]) / (2 * dx)
        f2[:-2, :] = f1[1:-1,:]
        ududx = 1 / hum * (mx1 * f1 + mx2 * f2)
        vdudy[:, 1:-1] = vu[:, 1:-1] * (uu[:, 2:] - uu[:, :-2]) / (2 * dy)
        
        wetu = hu > eps     # wetting & drying crit. (only for mom. balance)
        ur = uu[1,:]        # store velocity @ seaward boundary (given by BC)
        # Compute automatic timestep
        dt = CFL * np.min(dx / (np.sqrt(g*h) 
                            + np.maximum(vmagu,vmageu)))
        dt = np.maximum(dt, 0.05)
        t += dt
        if t >= tnext:
            dt -= t - tnext
            t = tnext
            tnext += tint
        
        # Explicit Euler step momentum u-direction
        uu[wetu] = uu[wetu] - dt * ( ududx[wetu] + vdudy[wetu]
                                     + g * dzsdx[wetu]
                                     + g / (Ch**2) / hu[wetu] * vmageu[wetu] * ueu[wetu]  # GLM approach
                                     - Fx[wetu] / rho / hu[wetu] )

        uu[not wetu] = 0

        uu[-1, :] = 0
        uu[-2, :] = 0   # reflection at the last grid line (usually dry)
        uu[1,:] = ur    # restore the seaward boundary condition
        qx = uu * hu    # flux @ u-point
        
        # Y-direction
        my = np.sign(vv)
        my[np.abs(vv) < umin] = 0
        my1 = (my + 1) / 2
        my2 = -(my - 1) / 2
        f1.fill(0)
        f2.fill(0)
        f1[:, :-1] = h[:, :-1]
        f2[:, :-1] = h[:, 1:]
        
        hv = my1 * f1 + my2 * f2    # h @ v-points for cont. eq.: upwind
        hvm = 0.5 * (f1 + f2)       # h @ v-points for mom. eq.: mean
        
        # Advection terms (momentum conserving method)
        udvdx[1:-1, :] = uv[1:-1,:] * (vv[2:, :] - vv[:-2, :]) / (2 * dx)
        f1.fill(0)
        f2.fill(0)
        f1[:, 1:-1] = 0.5 * (qy[:, 1:-1] + qy[:, :-2]) * (vv[:, 1:-1] - vv[:, :-2]) / dy
        f2[:, :-2] = f1[:, 1:-1]
        vdvdy = my1 * f1 + my2 * f2
        
        wetv = hv > eps  # wetting & drying crit. (only for mom. balance)
        
        # Explicit Euler step momentum v-direction
        vv[wetv] = vv[wetv] - dt * ( udvdx[wetv] + vdvdy[wetv] 
                                     + g * dzsdy[wetv] 
                                     + g / Ch**2 / hv[wetv] * vmagev[wetv] * vev[wetv] 
                                     - Fy[wetv] / rho / hv[wetv] )
        vv[not wetv] = 0
        
        qy = vv * hv    # flux @ v-points
        
        #  U and V @ cell centre
        u[1:-1, :] = 0.5 * (uu[:-2, :] + uu[1:-1, :])
        u[0, :] = uu[0, :]
        v[:, 1:-1] = 0.5 * (vv[:, :-2] + vv[:, 1:-1])
        v[:, 0] = vv[:, 0]
        
        # Ue and Ve @ cell centre
        ue[1:-1, :] = 0.5 * (ueu[:-2, :] + ueu[1:-1,:])
        ue[0, :] = ueu[0, :]
        ve[:, 1:-1] = 0.5 * (vev[:, :-2] + vev[:, 1:-1])
        ve[:, 0] = vev[:, 0]
        
        # Update water level using continuity eq.
        dzsdt[1:-1, 1:-1] = -(qx[1:-1, 1:-1] - qx[:-2, 1:-1]) / dx \
                - (qy[1:-1, 1:-1] - qy[1:-1, :-2]) / dy
        
        zs[1:-1, 1:-1] += dzsdt[1:-1, 1:-1] * dt
        
        # Transport of contaminants
        k = flowpar.k
        # @ boundaries
        c[0, :] = c[1, :]
        c[-1, :] = c[-2, :]
        c[:, 0] = c[:, 1]
        c[:, -1] = c[:, -2]
        
        # X-direction
        f1.fill(0)
        f2.fill(0)
        f1[:-1, :] = c[:-1, :]
        f2[:-1, :] = c[1:, :]   # water depth @ u-points
        cu = mx1 * f1 + mx2 * f2
        Su = cu * hu * uu
        Ecx[1:-1, :] = Etran * ( 
                0.5 * (h[:-2, :] + h[1:-1, :]) * (c[:-2, :] - c[1:-1, :]) * wetu[:-2, :] * wetu[1:-1, :]
                + 0.5 * (h[1:-1, :] + h[2:, :])  * (c[2:, :] - c[1:-1, :]) * wetu[2:, :] * wetu[1:-1, :]
                ) / (dx*dx) / h[1:-1, :]

        # Y-direction
        f1.fill(0)
        f2.fill(0)
        f1[:, :-1] = c[:, :-1]
        f2[:, :-1] = c[:, 1:]
        cv = my * f1 + my2 * f2
        Sv = cv * hv * vv
        Ecy[:, 1:-1] = Elong * (  0.5  * wetv[:, :-2] * wetv[:, 1:-1] * 
									(h[:, :-2] + h[:, 1:-1]) *	(c[:, :-2] - c[:, 1:-1])
								+ 0.5  * wetv[:, 2:] * wetv[:, 1:-1] *
									(h[:, 1:-1] + h[:, 2:]) * (c[:, 2:] - c[:, 1:-1])
								) / (dy*dy) / h[:, 1:-1]
        
        # Update using mass balance
        dcdt[1:-1, 1:-1] =  \
                - ((Su[1:-1, 1:-1] - Su[:-2, 1:-1]) / dx + \
                   (Sv[1:-1, 1:-1] - Sv[1:-1, :-2]) / dy ) / h[1:-1, 1:-1] \
                + Ecx[1:-1, 1:-1] \
                + Ecy[1:-1, 1:-1] \
                + W[1:-1, 1:-1] / (h[1:-1, 1:-1] * dx * dy) \
                - k * c[1:-1, 1:-1]

        c[1:-1, 1:-1] += dcdt[1:-1, 1:-1] * dt
        c[not wetv] = 0
        c[not wetu] = 0
        
        # Output
        self.uu = uu
        self.vv  = vv 
        self.ueu = ueu
        self.vev = vev 
        self.qx = qx
        self.qy = qy 
        self.vmagu = vmagu
        self.vmagv = vmagv
        self.u = u
        self.v = v 
        self.ue = ue
        self.ve = ve 
        self.zs = zs
        self.hold = h
        self.h = np.maximum(zs - zb, hmin)
        self.wetu = wetu
        self.wetv = wetv
        self.hu = hu
        self.hv = hv
        self.Su = Su
        self.Sv = Sv
        self.c = c

        flowpar.dt = dt
        flowpar.t = t
        flowpar.tnext = tnext
