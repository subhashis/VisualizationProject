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


import vtk
from vtk import *

def isInEllipsoid(cx,cy,cz,ax,ay,az,ba,ca,px,py,pz):
    CT = np.matrix([cx,cy,cz])
    XT = np.matrix([px,py,pz])
    C = CT.getT()
    X = XT.getT()
    XT_CT = XT - CT
    X_C = X - C
    a = np.sqrt([ax**2 + ay**2 + az**2])
    uax = ax/a
    uay = ay/a
    uaz = az/a
    b = ba*a
    c = ca*a
    A = np.matrix([[1.0/a**2,0.0,0.0],[0.0,1.0/b**2,0.0],[0.0,0.0,1.0/c**2]])
    #calculate rotation matrix
    unitX = np.matrix([1.0,0.0,0.0])
    ua = np.matrix([uax[0],uay[0],uaz[0]])
    #print "-------"
    #print unitX.shape
    #print ua.shape
    v = np.cross(unitX,ua)
    sine = np.linalg.norm(v)
    #cosine = np.dot(unitX,ua)
    cosine = uax[0]
    v1 = v[0,0]
    v2 = v[0,1]
    v3 = v[0,2]
    vx = np.matrix([[0.0,-v3,v2],[v3,0.0,-v1],[-v2,v1,0.0]])
    #print vx
    vx2 = vx*vx
    #print vx2
    I = np.matrix([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
    R = I + vx + vx2*(1-cosine)/sine**2
    #print R
    #R = np.matrix(np.array(R))
    RT = R.getT()
    #ellipsoid equation
    eq = XT_CT*RT*A*R*X_C
    
    e1=eq[0,0][0,0]
    #print e1
    if e1 <= 1.0:
        return True
    else:
        return False
        

def extractParticleinTime(timeStep,noh):
    
    if timeStep == 88:
        suffix = "1.0000"
    else:
        i = timeStep + 12
        suffix = "0." + str(i) + "00"
    
    path = "/home/subhashis/VisData/contestData/2015/ds14_scivis_0128_e4_dt04_" + suffix
    
    particles = load_sdf(path)
    
    h_100 = particles.parameters['h_100']
    width = particles.parameters['L0']
    cosmo_a = particles.parameters['a']
    
    kpc_to_Mpc = 1./1000
    
    # Define a simple function to convert proper to comoving Mpc/h.
    convert_to_cMpc = lambda proper: (proper + (width*cosmo_a)/2.) * h_100 * kpc_to_Mpc / cosmo_a
    
    nop = len(particles['x'])
    
    
        
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
    print "For %d" %timeStep
    for i in range(0,noh):
        if int_snap_num[i] == timeStep:
            halo_cx = x[i]
            halo_cy = y[i]
            halo_cz = z[i]
            #calculate the bounding box for faster computation
            xmax = halo_cx + 3.0
            xmin = halo_cx - 3.0
            ymax = halo_cy + 3.0
            ymin = halo_cy - 3.0
            zmax = halo_cz + 3.0
            zmin = halo_cz - 3.0
            
            #halo_radius = rvir[i]/1000
            hid = id[i]
            halo_ax = A_x[i]/1000
            halo_ay = A_y[i]/1000
            halo_az = A_z[i]/1000
            ba = b_to_a[i]
            ca = c_to_a[i]
            count = 0
            qcount = 0
            
            #print isInEllipsoid(halo_cx,halo_cy,halo_cz,halo_ax,halo_ay,halo_az,ba,ca,px,py,pz)
            for j in range(0,nop):
                px = convert_to_cMpc(particles['x'][j])
                py = convert_to_cMpc(particles['y'][j])
                pz = convert_to_cMpc(particles['z'][j])
                if px <= xmax and px >= xmin and py <= ymax and py >= ymin and pz <= zmax and pz >= zmin:
                    qcount += 1
                    if isInEllipsoid(halo_cx,halo_cy,halo_cz,halo_ax,halo_ay,halo_az,ba,ca,px,py,pz):
                        Points.InsertNextPoint(px,py,pz)
                        id_array.InsertNextTuple1(hid)
                        phi_array.InsertNextTuple1(particles['phi'][j])
                        count += 1
                #print j
            print "count=" + str(count) 
            print "qcount=" + str(qcount)
     
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
    #outputFile = "/home/subhashis/VisData/merger_trees/haloParticle.vtp"
    #outputFile = "/home/subhashis/VisData/merger_trees/haloParticleEllipsoid1.vtp"
    outputFile = "/home/subhashis/VisData/merger_trees/particleList/time" + str(timeStep) + ".vtp" 
        
    writer = vtk.vtkXMLPolyDataWriter();
    writer.SetFileName(outputFile);
    if vtk.VTK_MAJOR_VERSION <= 5:
        writer.SetInput(polydata)
    else:
        writer.SetInputData(polydata)
    writer.Write()
    print "Done generating output for time %d" %timeStep

def justLikeThat():
    print "jlt"
    
justLikeThat()


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
for k in range(0,87):
    extractParticleinTime(86-k,noh)

#extractParticleinTime(87,noh)

