from __future__ import division, print_function, absolute_import

import sys
import os
import subprocess
from subprocess import DEVNULL
import struct

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Script for pretex to build a mesh for a given radius gap ratio
PRETEXSCRIPT = """\
~/Projects/Nek5000/bin/pretex << EOF
newmesh
   1 READ PREVIOUS PARAMETERS
original
   1 BUILD FROM FILE
original
   3 GLOBAL REFINE
  11 STRETCH
   7 Stretch Outside Circle
.5  Input protected radius:
 -{}  {} Input new xminxmax (0,0 to scale): old -x; old +x
 -{}  {} Input new yminymax (0,0 to scale): old -y; old +y
   1 UP MENU
   1 END GLOBAL REFINE
   4 CURVE SIDES
   8 Convert Midside to Circle
500  Input maxium radius:
   1 BUILD MENU
   1 END    ELEMENTS
   1 ACCEPT MATL,QVOL
   1 ACCEPT B.C.s
   1 ACCEPT B.C.s
   1 EXIT
EOF
"""

# Script for n2to3 to extrude the 2d mesh. The script answers the
# following questions from n2to3.
#
# Input old (source) file name:
# Input new (output) file name:
# For ASCII output only: 0; For .rea/.re2: 1
# input number of levels: (1, 2, 3,...; < 0 for circular sweep.):
# input z min:
# input z max:
# input gain (0=custom,1=uniform,other=geometric spacing):
# input file containing ASCENDING z-values:
# This is for CEM: yes or no:
# Enter Z (5) boundary condition (P,PEC,PMC,PML):
# Enter Z (6) boundary condition (PEC,PMC,PML):
#
N2TO3SCRIPT = """\
~/Projects/Nek5000/bin/n2to3 << EOF
newmesh
{}
0
{}
0
{}
0
layers
yes
PML
PML
EOF
"""

# Script to create a map file. The script answers the following
# questions from genmap.
#
# Input .rea / .re2 name:
# Input mesh tolerance
GENMAPSCRIPT = """\
~/Projects/Nek5000/bin/genmap << EOF
{}
0.2
EOF
"""


def config_layers(rad, gap, ll, lb):
    # number of pml layers
    pml = 2
    # number of total layers in mesh
    num_lev = 2*pml + ll + 2*lb
    # wavelength of incoming wave
    lmbda = 633.0
    # height of the lens
    h = 155.0
    # z-coordinate of the first layer
    lay = [0]
    for i in range(1, num_lev):
        if i < pml+1 or i > num_lev-pml: # PML layers
            lay.append(lay[i-1]+(((lmbda/2)/pml)/rad)*0.5)
        elif ((i > pml and i < pml + lb + 1) or
              (i > num_lev - pml - lb and i < num_lev - pml + 1)): # Buffers
            lay.append(lay[i-1]+((lmbda/lb)/rad)*0.5)
        else: # lens
            lay.append(lay[i-1]+((h/ll)/rad)*0.5)
    zmax = (((3*lmbda)+h)/rad)*0.5 # total height of mesh
    lay.append(zmax)
    with open('layers', 'w') as f:
        for i in range(num_lev+1):
            f.write('{}\n'.format(lay[i]))
    return num_lev, zmax


def create_mesh(rad, gap, ll, lb):
    reastem = '{}_{}'.format(rad, gap)
    # Stretch the mesh in `original.rea'
    xy = 0.5*(rad + 0.5*gap)/rad
    script = PRETEXSCRIPT.format(xy, xy, xy, xy)
    subprocess.check_call(['bash', '-c', script], stdout=DEVNULL,
                          stderr=DEVNULL)
    # Extrude the mesh
    num_lev, zmax = config_layers(rad, gap, ll, lb)
    script = N2TO3SCRIPT.format(reastem, num_lev, zmax)
    subprocess.check_call(['bash', '-c', script], stdout=DEVNULL,
                           stderr=DEVNULL)
    # Generate the map file
    script = GENMAPSCRIPT.format(reastem)
    subprocess.check_call(['bash', '-c', script], stdout=DEVNULL,
                          stderr=DEVNULL)
    # Turn the Nek5000 read file into a NekCEM one
    reaname = '{}.rea'.format(reastem)
    with open('rea_start.txt', 'r') as f, open(reaname, 'r') as g:
        start = f.readlines()
        start[73] = '  {}\n'.format(rad*1e-9)
        end = g.readlines()[141:]
        readfile = start + end
    with open(reaname, 'w') as g:
        g.write(''.join(readfile))
    return reastem


def load_buf(filename):
    with open(filename, mode='rb') as f:
        content = f.read()
        data = [x for x in struct.iter_unpack('d', content)]
    return np.array(data)


def get_timeseries(index, firststep, laststep, step):
    series = []
    for i in range(firststep, laststep, step):
        filename = 'vtk/hx{:07d}.dat'.format(i)
        hx = load_buf(filename)
        series.append(hx[index])
    return np.squeeze(np.asarray(series))


class FittedSine():
    def __init__(self, tdata, hdata, kz, z, omega):
        self.tdata, self.hdata = tdata, hdata
        self.kz, self.z, self.omega = kz, z, omega
        self.fit()

    def f(self, t, A, phi):
        return A*np.sin(-self.kz*self.z - self.omega*t - phi)

    def jac(self, t, A, phi):
        out = np.empty((t.size, 2))
        arg = -self.kz*self.z - self.omega*t - phi
        out[:,0] = np.sin(arg)
        out[:,1] = -A*np.cos(arg)
        return out

    def fit(self):
        bounds = ([0, -np.inf], [np.inf, np.inf])
        popt, _ = curve_fit(self.f, self.tdata, self.hdata,
                            bounds=bounds)
        self.A, self.phi = popt

    def __call__(self, t):
        return self.f(t, self.A, self.phi)

    def __repr__(self):
        info = []
        info.append('A = {}'.format(self.A))
        info.append('phi = {}'.format(self.phi))
        return '\n'.join(info)


def get_amplitude_and_phi(index, firststep, laststep, step):
    dt = 8.49e-3
    z = -7.104
    omega = 0.993
    eps = 1.0
    mu = 1.0
    kz = np.sqrt(eps*mu)*omega

    hx = get_timeseries(index, firststep, laststep, step)
    times = dt*np.arange(firststep, laststep, step)
    hxinc = np.sin(-kz*z - omega*times)
    fit = FittedSine(times, hx, kz, z, omega)

    # plt.plot(times, hxinc, 'b', label='incident')
    # plt.plot(times, hx, 'r', label='scattered')
    # plt.plot(times, fit(times), 'g', label='fit')
    # plt.legend()
    # plt.show()

    return fit.A, fit.phi


def main():
    ll = 3
    lb = 3
    gaps = [46]
    radii = range(130, 181)

    firststep = 9000
    laststep = 10575
    step = 25

    x = load_buf('vtk/xcoordinates.dat')
    y = load_buf('vtk/ycoordinates.dat')
    index = x.size//2

    subprocess.call(['make', '-j4'], stdout=DEVNULL, stderr=DEVNULL)
    phis = []
    for gap in gaps:
        for rad in radii:
            reastem = create_mesh(rad, gap, ll, lb)
            subprocess.check_call(['nek', reastem, '4'])
            A, phi = get_amplitude_and_phi(index, firststep,
                                           laststep, step)
            os.remove('{}.rea'.format(reastem))
            os.remove('{}.map'.format(reastem))
            os.remove('{}.np=4.output'.format(reastem))
            phis.append(str(phi))
    with open('phi.txt', 'w') as f:
        f.write('\n'.join(phis))


if __name__ == '__main__':
    main()
