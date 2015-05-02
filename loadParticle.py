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
convert_to_cMpc = lambda proper: (proper + width/2.) * h_100 * kpc_to_Mpc


pl.scatter(convert_to_cMpc(particles['x'][sl]),convert_to_cMpc(particles['y'][sl]), color='b', s=1.0, alpha=0.05)

print convert_to_cMpc(particles['x'][2097151])
nop = len(particles['x'])

Points = vtk.vtkPoints()
phi_array = vtk.vtkDoubleArray()
#phi_array.SetName("phi")

for i in range(0,nop):
    Points.InsertNextPoint(convert_to_cMpc(particles['x'][i]), convert_to_cMpc(particles['y'][i]), convert_to_cMpc(particles['z'][i]))
    #phi_array.InsertNextTuple1(particles['phi'][i])

polydata = vtk.vtkPolyData()
polydata.SetPoints(Points)
#polydata.GetPointData().AddArray(phi_array)

if vtk.VTK_MAJOR_VERSION <= 5:
    polydata.Update()
    
#outputFile = "/home/subhashis/HaloTS_" + str(timeslice) + ".vtp"
outputFile = "/home/subhashis/VisData/cosmicParticle.vtp"
    
writer = vtk.vtkXMLPolyDataWriter();
writer.SetFileName(outputFile);
if vtk.VTK_MAJOR_VERSION <= 5:
    writer.SetInput(polydata)
else:
    writer.SetInputData(polydata)
writer.Write()