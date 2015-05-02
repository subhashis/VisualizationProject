from sdfpy import load_sdf
from thingking import loadtxt
import numpy as np
import matplotlib.pylab as pl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as colors
import matplotlib.cm as cmx


import vtk
from vtk import *

def isInSphere(cx,cy,cz,radius,px,py,pz):
    distance = (cx - px)**2 + (cy-py)**2 + (cz-pz)**2
    if distance <= radius**2:
        return True
    else:
        return False

particles = load_sdf("/home/subhashis/VisData/contestData/2015/ds14_scivis_0128_e4_dt04_1.0000")

# Now we want to convert the proper kpc of the particle position to comoving
# Mpc/h, a common unit used in computational cosmology in general, but
# specifically is used as the output unit in the merger tree halo list loaded
# in above. First we get the Hubble parameter, here stored as 'h_100' in the
# SDF parameters. Then we load the simulation width, L0, which is also in
# proper kpc. Finally we load the scale factor, a, which for this particular
# snapshot is equal to 1 since we are loading the final snapshot from the
# simulation. 
h_100 = particles.parameters['h_100']
width = particles.parameters['L0']
cosmo_a = particles.parameters['a']


kpc_to_Mpc = 1./1000
sl = slice(0,None)

# Define a simple function to convert proper to comoving Mpc/h.
convert_to_cMpc = lambda proper: (proper + width/2.) * h_100 * kpc_to_Mpc / cosmo_a


#pl.scatter(convert_to_cMpc(particles['x'][sl]),convert_to_cMpc(particles['y'][sl]), color='b', s=1.0, alpha=0.05)

print convert_to_cMpc(particles['x'][2097151])
nop = len(particles['x'])
print nop
print max(convert_to_cMpc(particles['x']))
print min(convert_to_cMpc(particles['x']))
print max(convert_to_cMpc(particles['y']))
print min(convert_to_cMpc(particles['y']))
print max(convert_to_cMpc(particles['z']))
print min(convert_to_cMpc(particles['z']))


#print isInSphere(4,33,1,2,5.785,33.0123,1)

grid = np.zeros((4,4,4))

#for z in range(0,64):
#    for y in range(0,64):
#        for x in range(0,64):
#            count = 0
#            for p in range(0,nop):
#                px = convert_to_cMpc(particles['x'][p])
#                py = convert_to_cMpc(particles['y'][p])
#                pz = convert_to_cMpc(particles['z'][p])
#                if(isInSphere(x,y,z,1.0,px,py,pz)):
#                    count += 1;
#            grid[z][y][x] = count;


#for z in range(0,2):
#    for y in range(0,2):
#        for z in range(0,2):
#            l = x - 1.5 #left
#            r = x + 1.5 #right
#            u = y - 1.5 #up
#            d = y + 1.5 #down
#            f = z - 1.5 #front
#            b = z + 1.5 #back
#            count = 0
#            for p in range(0,nop):
#                px = convert_to_cMpc(particles['x'][p])
#                py = convert_to_cMpc(particles['y'][p])
#                pz = convert_to_cMpc(particles['z'][p])
#                if px >= l and px <= r and py >= u and py <= d and pz <= f and pz >=b:
#                    count += 1
#            print "completed : " + str(count)

particlePos = np.zeros((nop*3))
j =0;
for i in range(0,nop):
    px = convert_to_cMpc(particles['x'][i])
    py = convert_to_cMpc(particles['y'][i])
    pz = convert_to_cMpc(particles['z'][i])
    particlePos[j] = px
    particlePos[j +1] = py
    particlePos[j +2] = pz
    j += 3;
       
            

#print grid
#grid.tofile("particlegrid.raw")
print particlePos[0]
print particlePos[1]
print particlePos[2]
#particlePos.tofile("particlePostion.raw")     


