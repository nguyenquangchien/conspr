#!/usr/bin/env python
# Module solver.py

# Class Flow and Domain
# Flow solver for 2DH domain on a uniform rectangular grid
# using finite volume method

# Written by Nguyen Quang Chien, 2009
# based on XBeach, http://xbeach.org/

# License GNU GPL, http://www.gnu.org/copyleft/gpl.html
 
import numpy
import copy

class Param(object):
    """ Parameters related to flow """
    def __init__(self):
        # hydraulic parameters
        self.Ch = 65
        self.rho = 1000
        self.eps = 0.1
        self.g = 9.81
        self.hmin = 0.001
        self.umin = 0.1
        self.CFL = 0.4
        self.zs0 = 0
        self.v0 = 1.0
        self.i = 1.E-4      # longitudinal slope
        
        # timing parameters
        self.tstart =0
        self.t = self.tstart
        self.tint = 1
        self.tnext = 0 + self.tint
        self.tstop = 100
        
        # contaminant parameters
        self.k = 0      # decay rate
        self.Elong = 0      # dispersion rate
        self.Etran = 0      # dispersion rate


class Domain(object):
    def __init__(self, flowpar, type='river', ext=None):
        self.makeBathy(flowpar, type, 4, 0.01, 0.5)
        
        # Forcing etc.
        domainDim = numpy.shape(self.zb)
        self.Fx = numpy.zeros(domainDim)  # need specified
        self.Fy = numpy.zeros(domainDim)
        self.W = numpy.zeros(domainDim)
        
        # Init flow condition
        zs0 = flowpar.zs0
        zs = numpy.zeros(domainDim) - flowpar.i * self.y
        self.zs = numpy.maximum(zs, self.zb)
        self.h = self.zs - self.zb
        self.dzsdt=numpy.zeros(domainDim)
        self.dzbdt=numpy.zeros(domainDim)
        self.uu=numpy.zeros(domainDim)
        self.vv=numpy.zeros(domainDim)
        self.vv = flowpar.Ch * numpy.sqrt(self.h * flowpar.i)
        self.vv = self.vv * (self.h > 0)
        self.qx = numpy.zeros(domainDim)
        self.qy = numpy.zeros(domainDim)
        self.c = numpy.zeros(domainDim)
    
    def makeBathy(self, flowpar, type, *arg):
        if type == 'river':
            if not 'ext' in locals():
                ext = (40, 40, 10, 40)
            (self.nx, self.ny, self.dx, self.dy) = ext
            domainDim = (self.nx+1, self.ny+1)
            self.x = numpy.zeros(domainDim)
            self.y = numpy.zeros(domainDim)
            self.zb = numpy.zeros(domainDim)
            h0 = arg[0]         # max. river depth
            slope_t = arg[1]    # transversal slope
            bank_height = arg[2]     # river bank
            for j in range(self.ny+1):      # 0 .. ny
                for i in range(self.nx+1):
                    self.x[i,j] = i * self.dx
                    self.y[i,j] = j * self.dy
                    self.zb[i,j] = -h0 + self.x[i,j] * slope_t
            self.zb[30:, :] = bank_height 
            self.zb -= flowpar.i * self.y     # Longitudinal slope i
        elif type == 'San Francisco bay':
            self.nx = self.ny = 199
            self.dx = 100
            self.dy = 200
            
            self.x, self.y = numpy.meshgrid(
                    numpy.arange(0,20000,self.dx), numpy.arange(0,40000,self.dy))
            self.zb = numpy.fromfile('newgrid.txt', sep=' ').reshape(self.x.shape)
            print(self.x.shape, self.y.shape, self.zb.shape)
            
    def solve(self, flowpar):
        x = copy.deepcopy(self.x)
        y = copy.deepcopy(self.y)
        dx = copy.deepcopy(self.dx)
        dy = copy.deepcopy(self.dy)
        uu  = copy.deepcopy(self.uu)
        vv = copy.deepcopy(self.vv)
        qx = copy.deepcopy(self.qx)
        qy = copy.deepcopy(self.qy)
        zb = copy.deepcopy(self.zb)
        zs = copy.deepcopy(self.zs)
        Fx = copy.deepcopy(self.Fx)
        Fy = copy.deepcopy(self.Fy)
        c = copy.deepcopy(self.c)
        W = copy.deepcopy(self.W)
        
        g  = flowpar.g
        Ch  = flowpar.Ch
        rho  = flowpar.rho
        umin = flowpar.umin
        hmin = flowpar.hmin
        eps  = flowpar.eps
        t     = flowpar.t
        tint = flowpar.tint
        tnext   =flowpar.tnext
        CFL     =flowpar.CFL
        Elong = flowpar.Elong
        Etran = flowpar.Etran
        
        domainDim = numpy.shape(x)
        vu      = numpy.zeros(domainDim)
        vsu     = numpy.zeros(domainDim)
        usu     = numpy.zeros(domainDim)
        vsv     = numpy.zeros(domainDim)
        usv     = numpy.zeros(domainDim)
        uv      = numpy.zeros(domainDim)
        vmagu   = numpy.zeros(domainDim)
        vmague  = numpy.zeros(domainDim)
        vmagv   = numpy.zeros(domainDim)
        vmagve  = numpy.zeros(domainDim)
        veu     = numpy.zeros(domainDim)
        ueu     = numpy.zeros(domainDim)
        vev     = numpy.zeros(domainDim)
        uev     = numpy.zeros(domainDim)

        u       = numpy.zeros(domainDim)
        v       = numpy.zeros(domainDim)
        dzsdx   = numpy.zeros(domainDim)
        dzsdy   = numpy.zeros(domainDim)
        dzsdt   = numpy.zeros(domainDim)
        dcdt   = numpy.zeros(domainDim)
        mx      = numpy.zeros(domainDim)
        mx1     = numpy.zeros(domainDim)
        mx2     = numpy.zeros(domainDim)
        my      = numpy.zeros(domainDim)
        my1     = numpy.zeros(domainDim)
        my2     = numpy.zeros(domainDim)
        f1      = numpy.zeros(domainDim)
        f2      = numpy.zeros(domainDim)
        ududx   = numpy.zeros(domainDim)
        vdudy   = numpy.zeros(domainDim)
        udvdx   = numpy.zeros(domainDim)
        vdvdy   = numpy.zeros(domainDim)
        qxm     = numpy.zeros(domainDim)
        qym     = numpy.zeros(domainDim)
        us      = numpy.zeros(domainDim)
        vs      = numpy.zeros(domainDim)
        ue      = numpy.zeros(domainDim)
        ve      = numpy.zeros(domainDim)
        Ecx      = numpy.zeros(domainDim)
        Ecy      = numpy.zeros(domainDim)
        
        # V-velocities at u-points
        vu[:-1, 1:-1] = 0.25 * (vv[:-1, :-2] 
                + vv[:-1, 1:-1] + vv[1:, :-2] + vv[1:, 1:-1])
        
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
        # V-euler velocities at u-point
        veu = vu - vsu
        # U-euler velocties at u-point
        ueu = uu - usu
        # Velocity magnitude at u-points
        vmagu = numpy.sqrt(uu**2 + vu**2)
        # Eulerian velocity magnitude at u-points
        vmageu = numpy.sqrt(ueu**2 + veu**2)

        # U-velocities at v-points
        uv[1:-1, :-1] = 0.25 * (uu[:-2, :-1] + uu[1:-1, :-1] 
                + uu[:-2, 1:] + uu[1:-1, 1:])

        # boundaries?
        uv[:, -1] = uv[:, -2]
        # V-stokes velocities at V point
        vsv[1:-1, :-1] = 0
        vsv[:, 0] = vsv[:, 1]
        vsv[:, -1] = vsv[:, -2]
        # U-stokes velocities at V point
        usv[1:-1, :-1] = 0
        usv[:, 0] = usv[:, 1]
        usv[:, -1] = usv[:, -2]

        # V-euler velocities at V-point
        vev = vv - vsv
        # U-euler velocties at V-point
        uev = uv - usv
        # Velocity magnitude at v-points
        vmagv = numpy.sqrt(uv**2 + vv**2)
        # Eulerian velocity magnitude at v-points
        vmagev = numpy.sqrt(uev**2 + vev**2)
        # Water level slopes
        dzsdx[:-1, :] = (zs[1:, :] - zs[:-1, :]) / dx
        dzsdy[:, :-1] = (zs[:, 1:] - zs[:, :-1]) / dy
        
        # Upwind method implemented through weight factors mx1 and mx2
        
        # Water depth
        h = zs - zb
        h = numpy.maximum(h, hmin)
        weth = h > hmin
        # X-direction
        mx = numpy.sign(uu)
        mx[numpy.abs(uu) < umin] = 0 
        mx1 = (mx + 1) / 2
        mx2 = -(mx - 1) / 2
        f1[:-1, :] = h[:-1, :]
        f2[:-1, :] = h[1:, :]
        # Water depth in u-points for continuity equation: upwind
        hu = mx1 * f1 + mx2 * f2
        # Water depth in u-points for momentum equation: mean
        hum = numpy.max(0.5 * (f1 + f2), hmin)
        # Advection terms (momentum conserving method)
        f1.fill(0)
        f2.fill(0)
        f1[1:-1, :] = (qx[1:-1, :] + qx[:-2,:]) * (uu[1:-1, :] - uu[:-2,:]) / (2 * dx)
        f2[:-2, :] = f1[1:-1,:]
        ududx = 1 / hum * (mx1 * f1 + mx2 * f2)
        vdudy[:, 1:-1] = vu[:, 1:-1] * (uu[:, 2:] - uu[:, :-2]) / (2 * dy)
        # Wetting and drying criterion (only for momentum balance)
        wetu = hu > eps
        # Store velocity at seaward boundary (given by boundary condition)
        ur = uu[1,:]
        # Compute automatic timestep
        dt = CFL * numpy.min(dx / (numpy.sqrt(g*h) 
                + numpy.maximum(vmagu,vmageu)))
        dt = numpy.maximum(dt, 0.05)
        t += dt
        if t >= tnext:
            dt -= t - tnext
            t = tnext
            tnext += tint
        
        # Explicit Euler step momentum u-direction
        uu[wetu] = uu[wetu] - dt * ( ududx[wetu] + vdudy[wetu]
                             + g * dzsdx[wetu]
                             + g / (Ch**2) / hu[wetu] * vmageu[wetu] * ueu[wetu]  # GLM approach + g/C^2./hu(wetu).*vmagu(wetu).*uu(wetu) ... 
                             - Fx[wetu] / rho / hu[wetu] )
        
        uu[wetu==False] = 0
        
        uu[-1, :] = 0
        uu[-2, :] = 0   # reflection at the last grid line (usually dry)
        # Restore seaward boundary condition
        uu[1,:] = ur
        # Flux in u-point
        qx = uu * hu
        
        # Y-direction
        my = numpy.sign(vv)
        my[numpy.abs(vv) < umin] = 0
        my1 = (my + 1) / 2
        my2 = -(my - 1) / 2
        f1.fill(0)
        f2.fill(0)
        f1[:, :-1] = h[:, :-1]
        f2[:, :-1] = h[:, 1:]
        # Water depth in v-points for continuity equation: upwind
        hv = my1 * f1 + my2 * f2
        # Water depth in v-points for momentum equation: mean
        hvm = 0.5 * (f1 + f2)
        # Advection terms (momentum conserving method)
        udvdx[1:-1, :] = uv[1:-1,:] * (vv[2:, :] - vv[:-2, :]) / (2 * dx)
        
        f1.fill(0)
        f2.fill(0)
        f1[:, 1:-1] = 0.5 * (qy[:, 1:-1] + qy[:, :-2]) * (vv[:, 1:-1] - vv[:, :-2]) / dy
        f2[:, :-2] = f1[:, 1:-1]
        vdvdy = my1 * f1 + my2 * f2
        # Wetting and drying criterion (only for momentum balance)
        wetv = hv > eps
        
        # Explicit Euler step momentum v-direction
        vv[wetv] = vv[wetv] - dt * ( udvdx[wetv] + vdvdy[wetv] 
                              + g * dzsdy[wetv] 
                              + g / Ch**2 / hv[wetv] * vmagev[wetv] * vev[wetv] 
                              - Fy[wetv] / rho / hv[wetv] )
        vv[wetv==False] = 0
        
        # Flux in v-points
        qy = vv * hv
        
        #  U and V in cell centre
        u[1:-1, :] = 0.5 * (uu[:-2, :] + uu[1:-1, :])
        u[0, :] = uu[0, :]
        v[:, 1:-1] = 0.5 * (vv[:, :-2] + vv[:, 1:-1])
        v[:, 0] = vv[:, 0]
        
        # Ue and Ve in cell centre
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
        # at boundaries
        c[0, :] = c[1, :]
        c[-1, :] = c[-2, :]
        c[:, 0] = c[:, 1]
        c[:, -1] = c[:, -2]
        
        # X-direction
        f1.fill(0)
        f2.fill(0)
        f1[:-1, :] = c[:-1, :]
        f2[:-1, :] = c[1:, :]
        # water depth in u-points
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
        Ecy[:, 1:-1] = Elong * ( 
                0.5 * (h[:, :-2] + h[:, 1:-1]) * (c[:, :-2] - c[:, 1:-1]) * wetv[:, :-2] * wetv[:, 1:-1]
                + 0.5 * (h[:, 1:-1] + h[:, 2:])  * (c[:, 2:] - c[:, 1:-1]) * wetv[:, 2:] * wetv[:, 1:-1]
                ) / (dy*dy) / h[:, 1:-1]
        
        # Update using mass balance
        dcdt[1:-1, 1:-1] =  \
                - ((Su[1:-1, 1:-1] - Su[:-2, 1:-1]) / dx + (Sv[1:-1, 1:-1] - Sv[1:-1, :-2]) / dy ) / h[1:-1, 1:-1] \
                + Ecx[1:-1, 1:-1] \
                + Ecy[1:-1, 1:-1] \
                + W[1:-1, 1:-1] / (h[1:-1, 1:-1] * dx * dy) \
                - k * c[1:-1, 1:-1]
        
        c[1:-1, 1:-1] += dcdt[1:-1, 1:-1] * dt
        c[wetv==False] = 0
        c[wetu==False] = 0
        
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
        self.h = numpy.maximum(zs - zb, hmin)
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
    