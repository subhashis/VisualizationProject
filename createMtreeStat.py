# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import vtk
from vtk import *
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerLine2D
import matplotlib.cm as cm

file_name = "/home/subhashis/VisData/merger_trees/firsttreewithUmass.vtp"
 
# Read the source file.
reader = vtk.vtkXMLPolyDataReader()
reader.SetFileName(file_name)
reader.Update() # Needed because of GetScalarRange
#output = reader.GetOutput()
#scalar_range = output.GetScalarRange()
#pt_array = vtk_to_numpy(reader.GetOutput.GetPoint())
haloid = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray('haloid'))
mvir = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray('mvir'))
particleCount = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray('particleCount'))
snap_num = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray('Snap_num'))
totalParMass = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray('totalParMass'))
uMass = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray('uMass'))


noh = len(haloid)
xaxis = range(0,noh)
#area = [3.14159, 12.56636, 28.27431, 50.26544, 78.53975, 113.09724]
#plt.plot(xaxis, np.sort(mvir))
#plt.plot(xaxis, np.sort(uMass))
#plt.show()

virial_mass = np.zeros(89)

for i in range(0,89):
    count = 0
    for j in range(0,noh):
        if snap_num[j] == i:
            count += 1
            virial_mass[i] += mvir[j]
    #virial_mass[i] /= count
            
    #print count

#print virial_mass

totalParticle_mass = np.zeros(89)

for i in range(0,89):
    count = 0
    for j in range(0,noh):
        if snap_num[j] == i:
            count += 1
            totalParticle_mass[i] += totalParMass[j]
    #totalParticle_mass[i] /= count
            
u_mass = np.zeros(89)

for i in range(0,89):
    count = 0
    for j in range(0,noh):
        if snap_num[j] == i:
            count += 1
            u_mass[i] += uMass[j]
    #u_mass[i] /= count
            
    

plt.plot(range(0,89),virial_mass,label='Virial Mass')
plt.plot(range(0,89),totalParticle_mass,label='DM Particle Mass')
#plt.plot(range(0,89),u_mass)
plt.xlabel('Timestep')
plt.ylabel('Mass 1e14 Msun/h')

plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)
plt.show()

