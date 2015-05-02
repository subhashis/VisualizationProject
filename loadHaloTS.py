# -*- coding: utf-8 -*-

from sdfpy import load_sdf
from thingking import loadtxt
import numpy as np


import vtk
from vtk import *

def createVTKfile(timeslice):
    
    fileNo = "/home/subhashis/VisData/hlists/" + str(timeslice) + ".list"
    
    scale, id, desc_scale, desc_id, num_prog, pid, upid, desc_pid, phantom, \
    sam_mvir, mvir, rvir, rs, vrms, mmp, scale_of_last_MM, vmax, x, y, z, \
    vx, vy, vz, Jx, Jy, Jz, Spin, Breadth_first_ID, Depth_first_ID, \
    Tree_root_ID, Orig_halo_ID, Snap_num, Next_coprogenitor_depthfirst_ID, \
    Last_progenitor_depthfirst_ID, Rs_Klypin, M_all, M200b, M200c, M500c, \
    M2500c, Xoff, Voff, Spin_Bullock, b_to_a, c_to_a, A_x, A_y, A_z, \
    b_to_a_500c, c_to_a_500c, A_x_500c, A_y_500c, A_z_500c, T_over_U, \
    M_pe_Behroozi, M_pe_Diemer, Macc, Mpeak, Vacc, Vpeak, Halfmass_Scale, \
    Acc_Rate_Inst, Acc_Rate_100Myr, Acc_Rate_Tdyn = \
    loadtxt(fileNo, unpack=True)
    
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
    
    
    for i in range(0,len(x)):
        Points.InsertNextPoint(x[i],y[i],z[i])
        id_array.InsertNextTuple1(id[i])
        pid_array.InsertNextTuple1(pid[i])
        mvir_array.InsertNextTuple1(mvir[i])
        rvir_array.InsertNextTuple1(rvir[i])
        #velo = [vx[i],vy[i],vz[i]]
        velocity_array.InsertNextTuple3(vx[i],vy[i],vz[i])
    
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
    polydata.GetPointData().AddArray(pid_array)
    polydata.GetPointData().AddArray(mvir_array)
    polydata.GetPointData().AddArray(rvir_array)
    #polydata.GetPointData().SetScalars(hostHaloId)
    polydata.GetPointData().AddArray(velocity_array)
    polydata.GetPointData().SetVectors(velocity_array)
    
    if vtk.VTK_MAJOR_VERSION <= 5:
        polydata.Update()
        
    #outputFile = "/home/subhashis/HaloTS_" + str(timeslice) + ".vtp"
    outputFile = "/home/subhashis/VisData/HaloTS/HaloTS_" + str(timeslice) + ".vtp"
        
    writer = vtk.vtkXMLPolyDataWriter();
    writer.SetFileName(outputFile);
    if vtk.VTK_MAJOR_VERSION <= 5:
        writer.SetInput(polydata)
    else:
        writer.SetInputData(polydata)
    writer.Write()



for i in range(0,88):
    createVTKfile(i+1)
#createVTKfile(89)
#readData(88)