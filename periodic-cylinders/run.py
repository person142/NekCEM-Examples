# python makemesh.py rad gap lev_lens lev_buff
# This script uses original.rea. If you want to use a different original mesh,
# you can change original to your .rea name in lines 17 and 19
# The output mesh file is rad_gap.rea.
from __future__ import division, print_function, absolute_import

import sys
import subprocess
import analysis

# Script for pretex to build a mesh for a given radius gap ratio
PRETEXSCRIPT = """\
pretex << EOF
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

# Script for n2to3 to extrude the 2d mesh
# THE FOLLOWING ARE THE QUESTIONS THAT N2TO3SCRIPT GIVES ANSWERS TO	 
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
N2TO3SCRIPT = """\
n2to3 << EOF
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

# Script to create .map file
# Input .rea / .re2 name:
# Input mesh tolerance
GENMAPSCRIPT = """\
genmap << EOF
{}
0.2
EOF
"""

MAKENEKSCRIPT = "./makenek cylinder"

NEKSCRIPT = "./nek {} 4"


def main():
    run_make()
    ll = 3
    lb = 3
    gaps = [46]
    radii = [180]
    for gap in gaps:
        for rad in radii:
            reaname = create_mesh(rad, gap, ll, lb)
            run_nek(reaname)
            phi = analysis.main()
            with open('PHI', 'a+') as f:
                f.write(str(rad) + '_' + str(gap)+' : ' + str(phi)+"\n") 


def run_make():
    script = MAKENEKSCRIPT
    subprocess.call(['bash', '-c', script])


def run_nek(reaname):
    script = NEKSCRIPT.format(reaname.rstrip(".rea"))
    subprocess.call(['bash', '-c', script])


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
        if i < pml+1 or i > num_lev-pml: #PML layers
            lay.append(lay[i-1]+(((lmbda/2)/pml)/rad)*0.5)
        elif ((i > pml and i < pml + lb + 1) or
              (i > num_lev - pml - lb and i < num_lev - pml + 1)): #Buffers
            lay.append(lay[i-1]+((lmbda/lb)/rad)*0.5)
        else: #lens
            lay.append(lay[i-1]+((h/ll)/rad)*0.5)
    zmax = (((3*lmbda)+h)/rad)*0.5 #total height of mesh
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
    subprocess.call(['bash', '-c', script])
    # Extrude the mesh
    num_lev, zmax = config_layers(rad, gap, ll, lb)
    script = N2TO3SCRIPT.format(reastem, num_lev, zmax)
    subprocess.call(['bash', '-c', script])
    # Generate the map file
    script = GENMAPSCRIPT.format(reastem)
    subprocess.call(['bash', '-c', script])
    # Turn the Nek5000 read file into a NekCEM one
    reaname = '{}.rea'.format(reastem)
    with open('rea_start.txt', 'r') as f, open(reaname, 'r') as g:
        start = f.readlines()
        start[73] = '  {}\n'.format(rad*1e-9)
        end = g.readlines()[141:]
        readfile = start + end
    with open(reaname, 'w') as g:
        g.write(''.join(readfile))
    return reaname


if __name__ == '__main__':
    main()
