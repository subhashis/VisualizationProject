from sdfpy import load_sdf
from thingking import loadtxt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import time, sys


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
convert_to_cMpc = lambda proper: (proper + (width*cosmo_a)/2.) * h_100 * kpc_to_Mpc / cosmo_a

nop = len(particles['x'])
print nop


#prefix = "/home/subhashis/VisData/contestData/2015/rockstar/trees/tree_0_0_0.dat"
#prefix = "/home/subhashis/VisData/merger_trees/tree_0_0_0.dat"
#loaded just one tree in test.dat
prefix = "/home/subhashis/VisData/merger_trees/test.dat"
# Load the a=1 Rockstar hlist file. The header of the file lists the useful
# units/information.
scale, id, desc_scale, desc_id, num_prog, pid, upid, desc_pid, phantom, \
    sam_mvir, mvir, rvir, rs, vrms, mmp, scale_of_last_MM, vmax, x, y, z, \
    vx, vy, vz, Jx, Jy, Jz, Spin, Breadth_first_ID, Depth_first_ID, \
    Tree_root_ID, Orig_halo_ID, Snap_num, Next_coprogenitor_depthfirst_ID, \
    Last_progenitor_depthfirst_ID, Rs_Klypin, M_all, M200b, M200c, M500c, \
    M2500c, Xoff, Voff, Spin_Bullock, b_to_a, c_to_a, A_x, A_y, A_z, \
    b_to_a_500c, c_to_a_500c, A_x_500c, A_y_500c, A_z_500c, T_over_U, \
    M_pe_Behroozi, M_pe_Diemer = \
    loadtxt(prefix, skiprows=0, unpack=True)
    
print len(id)   
print max(pid)
print min(pid)
print len(Tree_root_ID)
print max(Tree_root_ID)
print min(Tree_root_ID)
maxMass = max(mvir)

for i in range(0,len(x)):
    if maxMass == mvir[i]:
        print id[i]


#fig = plt.figure()
#ax = fig.gca(projection='3d')
#numPoints = len(x)
#for i in range(0,numPoints-1):
#    ax.plot([x[i],x[i+1]],[y[i],y[i+1]],[z[i],z[i+1]],'r-',lw=1.5)
#ax.set_xlim3d(0,70)
#ax.set_ylim3d(0,70)
#ax.set_zlim3d(0,70)
#plt.show()


    
Points = vtk.vtkPoints()
id_array = vtk.vtkIntArray()
#id_array.SetNumberofComponents(1)
id_array.SetName("haloid")

phi_array = vtk.vtkDoubleArray()
phi_array.SetName("phi")
#
#rvir_array = vtk.vtkDoubleArray()
#rvir_array.SetName("rvir")
#
#pid_array = vtk.vtkIntArray()
#pid_array.SetName("pid")
#
#velocity_array = vtk.vtkDoubleArray()
#velocity_array.SetNumberOfComponents(3)
#velocity_array.SetName("v")

halo_cx = x[0]
halo_cy = y[0]
halo_cz = z[0]
halo_radius = rvir[i]/1000
hid = id[0]
count = 0
for i in range(0,nop):
    px = convert_to_cMpc(particles['x'][i])
    py = convert_to_cMpc(particles['y'][i])
    pz = convert_to_cMpc(particles['z'][i])
    if isInSphere(halo_cx,halo_cy,halo_cz,halo_radius,px,py,pz):
        Points.InsertNextPoint(px,py,pz)
        id_array.InsertNextTuple1(hid)
        phi_array.InsertNextTuple1(particles['phi'][i])
        count += 1
print "count=" + str(count)

 
#for i in range(0,len(x)):
#    Points.InsertNextPoint(x[i],y[i],z[i])
#    id_array.InsertNextTuple1(id[i])
#    pid_array.InsertNextTuple1(pid[i])
#    mvir_array.InsertNextTuple1(mvir[i])
#    rvir_array.InsertNextTuple1(rvir[i])
#    #velo = [vx[i],vy[i],vz[i]]
#    velocity_array.InsertNextTuple3(vx[i],vy[i],vz[i])
    
#draw_3d_lines(x,y,z)

#for i in range(0,len(x)):
#    if id[i] == 679582:
#        Points.InsertNextPoint(x[i],y[i],z[i])
#    elif pid[i] == 679582:
#        Points.InsertNextPoint(x[i],y[i],z[i])
#    else:
#        ran=2
  
polydata = vtk.vtkPolyData()
polydata.SetPoints(Points)
polydata.GetPointData().AddArray(id_array)
#polydata.GetPointData().AddArray(pid_array)
polydata.GetPointData().AddArray(phi_array)
#polydata.GetPointData().AddArray(rvir_array)
##polydata.GetPointData().SetScalars(hostHaloId)
#polydata.GetPointData().AddArray(velocity_array)
#polydata.GetPointData().SetVectors(velocity_array)

if vtk.VTK_MAJOR_VERSION <= 5:
    polydata.Update()
    
#outputFile = "/home/subhashis/HaloTS_" + str(timeslice) + ".vtp"
outputFile = "/home/subhashis/VisData/merger_trees/haloParticle.vtp"
    
writer = vtk.vtkXMLPolyDataWriter();
writer.SetFileName(outputFile);
if vtk.VTK_MAJOR_VERSION <= 5:
    writer.SetInput(polydata)
else:
    writer.SetInputData(polydata)
writer.Write()

 



