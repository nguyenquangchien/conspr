import numpy

K1 = 0.51444       # knot to m/s


def to_sec(hourmin):
    return (hourmin // 100) * 3600 + (hourmin % 100) * 60


# ATTN: water level in Feet above MLLW
ts1 = [(-51, 7.3), (622, -1.8), (1340, 5.0), (1757, 2.8), (2400, 7.4)]
ts1 = [(to_sec(datapt[0]), datapt[1] * 0.3048) for datapt in ts1]
ts2 = [(-32, 9.4), (700, -1.9), (1359, 6.4), (1835, 2.9), (2419, 9.5)]  # San Mateo
ts2 = [(to_sec(datapt[0]), datapt[1] * 0.3048) for datapt in ts2]
ts3 = [(-26, 10.4), (704, -1.9), (1405, 7.1), (1839, 2.9), (2425, 10.6)]  # Dumbarton Bridge
ts3 = [(to_sec(datapt[0]), datapt[1] * 0.3048) for datapt in ts3]


def init(par, dom):
    par.Elong = 20
    par.Etran = 20
    par.hmin = 0.005
    # set up a linear water slope
    zsSea = interp(ts1, par.t)
    zsRiv = interp(ts3, par.t) 
    slope = (zsSea - zsRiv) / float(dom.y[-1,0] - dom.y[0,0])
    print(slope)
    dom.zs = zsRiv + dom.y * slope
    dom.zs = numpy.maximum(dom.zs, dom.zb + par.hmin)
    dom.uu.fill(0)  # still water
    dom.vv.fill(0)


def interp(dataXY, xp):
    xList, yList = list(zip(*dataXY)[0]), list(zip(*dataXY)[1])
    if xp < xList[0]:
        return yList[0]
    elif xp > xList[-1]:
        return yList[-1]
    else:
        for i in range(len(xList) - 1):
            if xList[i] <= xp <= xList[i+1]:
                return yList[i] + ( (yList[i+1] - yList[i]) 
                        * (xp - xList[i]) / (xList[i+1] - xList[i]) )


def applyBC(par, dom):
    t = par.t
    
    # North boundary - Yerba Buena Island (29d)
    vts1 = [ (-30, 0), (213, -1.1), (543, 0), (1051, 1.3),
            (1353, 0), (1545, -0.5), (1758, 0), (2137, 1.0),
            (2423, 0), ]
    vts1 = [(to_sec(datapt[0]), datapt[1] * K1 / 1.1) for datapt in vts1]
    vel1 = interp(vts1, t)
    dom.vv[-1] = vel1 * (dom.h[0] > par.hmin)
    dom.uu[-1,:] = 0
    
    wl1 = interp(ts1, t)
    dom.zs[-1] = numpy.maximum(wl1, dom.zb[-1] + par.hmin)
    
    # South boundary - Dumbarton bridge
    vts2 = [ (-123, 1.6), (21, 0), (343, -1.9), (805, 0),
            (1052, 1.8), (1429, 0), (1712, -0.7), (1911, 0),
            (2227, 1.6), (2508, 0),]
    vts2 = [(to_sec(datapt[0]), datapt[1] * K1 / 1.1) for datapt in vts2]
    vel2 = interp(vts2, t)
    dom.vv[0] = vel2 * (dom.h[0] > par.hmin)
    dom.uu[0,:] = 0

    wl3 = interp(ts3, t)
    dom.zs[0] = numpy.maximum(wl3, dom.zb[0] + par.hmin)

    # Block lateral boundaries
    dom.uu[:,0] = 0
    dom.uu[:,-1] = 0
    dom.vv[:,0] = 0
    dom.vv[:,-1] = 0
