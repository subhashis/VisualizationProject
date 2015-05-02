# -*- coding: utf-8 -*-
from sdfpy import load_sdf
from thingking import loadtxt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import time, sys
import csv


import vtk
from vtk import *

from vtk.util.numpy_support import vtk_to_numpy

def getParticleMass(ts):
    if ts == 88:
        suffix = "1.0000"
    else:
        i = ts + 12
        suffix = "0." + str(i) + "00"
    
    path = "/home/subhashis/VisData/contestData/2015/ds14_scivis_0128_e4_dt04_" + suffix
    
    particles = load_sdf(path)
    
    h_100 = particles.parameters['h_100']
    width = particles.parameters['L0']
    cosmo_a = particles.parameters['a']
    pm = particles.parameters['particle_mass']
    pm = pm / cosmo_a
    return pm

def getParticleCount(hid,ts):
    f_path = "/home/subhashis/VisData/merger_trees/particleList/time" + str(ts) + ".vtp"
    # Read the source file.
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(f_path)
    reader.Update()
    hid_array = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray('haloid'))
    totalNum = len(hid_array)
    count = 0
    for i in range(0,totalNum):
        if hid == hid_array[i]:
            count += 1
    return count
    
    

prefix = "/home/subhashis/VisData/merger_trees/test.dat"
scale, id, desc_scale, desc_id, num_prog, pid, upid, desc_pid, phantom, \
    sam_mvir, mvir, rvir, rs, vrms, mmp, scale_of_last_MM, vmax, x, y, z, \
    vx, vy, vz, Jx, Jy, Jz, Spin, Breadth_first_ID, Depth_first_ID, \
    Tree_root_ID, Orig_halo_ID, Snap_num, Next_coprogenitor_depthfirst_ID, \
    Last_progenitor_depthfirst_ID, Rs_Klypin, M_all, M200b, M200c, M500c, \
    M2500c, Xoff, Voff, Spin_Bullock, b_to_a, c_to_a, A_x, A_y, A_z, \
    b_to_a_500c, c_to_a_500c, A_x_500c, A_y_500c, A_z_500c, T_over_U, \
    M_pe_Behroozi, M_pe_Diemer = \
    loadtxt(prefix, skiprows=0, unpack=True)

int_snap_num = Snap_num.astype(int)
noh = len(id)
#haloID = id[0]
#timeSnap = int_snap_num[0]
#particleMass = getParticleMass(timeSnap)
#particleCount = getParticleCount(haloID,timeSnap)
#print particleMass
#print particleCount

realcount = 0
for i in range(0,noh):
    if int_snap_num[i] >= 23:
        realcount += 1
print realcount

particleMass = np.zeros(noh)
particleCount = np.zeros(noh)
totalParMass = np.zeros(noh)
uncertainMass = np.zeros(noh)

for i in range(0,noh):
    haloID = id[i]
    timeSnap = int_snap_num[i]
    particleMass[i] = getParticleMass(timeSnap) * 10**10
    particleCount[i] = getParticleCount(haloID,timeSnap)
    totalParMass[i] = particleMass[i]*particleCount[i]
    uncertainMass[i] = mvir[i] - totalParMass[i]


#plotting graphs
#xaxis = range(0,noh)
#plt.plot(xaxis, mvir)
#plt.show()

Points = vtk.vtkPoints()
id_array = vtk.vtkIntArray()
#id_array.SetNumberofComponents(1)
id_array.SetName("haloid")

mvir_array = vtk.vtkDoubleArray()
mvir_array.SetName("mvir")

rvir_array = vtk.vtkDoubleArray()
rvir_array.SetName("rvir")

pid_array = vtk.vtkIntArray()
pid_array.SetName("pid")

velocity_array = vtk.vtkDoubleArray()
velocity_array.SetNumberOfComponents(3)
velocity_array.SetName("v")

snapID_array = vtk.vtkIntArray()
snapID_array.SetName("Snap_num")

ba_array = vtk.vtkDoubleArray()
ba_array.SetName("ba")

ca_array = vtk.vtkDoubleArray()
ca_array.SetName("ca")

axis_array = vtk.vtkDoubleArray()
axis_array.SetNumberOfComponents(3)
axis_array.SetName("Axis")

particleCount_array = vtk.vtkDoubleArray()
particleCount_array.SetName("particleCount")
totalParMass_array = vtk.vtkDoubleArray()
totalParMass_array.SetName("totalParMass")
uMass_array = vtk.vtkDoubleArray()
uMass_array.SetName("uMass")

for i in range(0,noh):
    Points.InsertNextPoint(x[i],y[i],z[i])
    id_array.InsertNextTuple1(id[i])
    pid_array.InsertNextTuple1(pid[i])
    mvir_array.InsertNextTuple1(mvir[i])
    rvir_array.InsertNextTuple1(rvir[i]/1000)
    snapID_array.InsertNextTuple1(Snap_num[i])
    ba_array.InsertNextTuple1(b_to_a[i])
    ca_array.InsertNextTuple1(c_to_a[i])
    #velo = [vx[i],vy[i],vz[i]]
    velocity_array.InsertNextTuple3(vx[i],vy[i],vz[i])
    axis_array.InsertNextTuple3(A_x[i],A_y[i],A_z[i])
    particleCount_array.InsertNextTuple1(particleCount[i])
    totalParMass_array.InsertNextTuple1(totalParMass[i])
    uMass_array.InsertNextTuple1(uncertainMass[i])
    
    

  
polydata = vtk.vtkPolyData()
polydata.SetPoints(Points)
polydata.GetPointData().AddArray(id_array)
polydata.GetPointData().AddArray(pid_array)
polydata.GetPointData().AddArray(mvir_array)
polydata.GetPointData().AddArray(rvir_array)
polydata.GetPointData().AddArray(snapID_array)
polydata.GetPointData().AddArray(ba_array)
polydata.GetPointData().AddArray(ca_array)
polydata.GetPointData().AddArray(particleCount_array)
polydata.GetPointData().AddArray(totalParMass_array)
polydata.GetPointData().AddArray(uMass_array)
#polydata.GetPointData().SetScalars(hostHaloId)
polydata.GetPointData().AddArray(velocity_array)
polydata.GetPointData().SetVectors(velocity_array)
polydata.GetPointData().AddArray(axis_array)

if vtk.VTK_MAJOR_VERSION <= 5:
    polydata.Update()
    
#outputFile = "/home/subhashis/HaloTS_" + str(timeslice) + ".vtp"
outputFile = "/home/subhashis/VisData/merger_trees/firsttreewithUmass.vtp"
    
writer = vtk.vtkXMLPolyDataWriter();
writer.SetFileName(outputFile);
if vtk.VTK_MAJOR_VERSION <= 5:
    writer.SetInput(polydata)
else:
    writer.SetInputData(polydata)
writer.Write()
    
