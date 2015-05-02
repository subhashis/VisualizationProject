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

def readParticle(ts):
    f_path = "/home/subhashis/VisData/merger_trees/particleList/time" + str(ts) + ".vtp"
    # Read the source file.
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(f_path)
    reader.Update()
    hid_arr = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray('haloid'))
    phi_arr = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray('phi'))
    nop = len(hid_arr)
    x = np.zeros(nop)
    y = np.zeros(nop)
    z = np.zeros(nop)
    for index in range(nop):
        pt = [0,0,0]
        reader.GetOutput().GetPoint(index, point)
        x[index] = pt[0]
        y[index] = pt[1]
        z[index] = pt[2]
    
    return x,y,z,hid_arr,phi_arr

def createTrail(ts):
    Points = vtk.vtkPoints()
    id_array = vtk.vtkIntArray()
    #id_array.SetNumberofComponents(1)
    id_array.SetName("haloid")
    
    phi_array = vtk.vtkDoubleArray()
    phi_array.SetName("phi")
    for i in range(0,ts+1):
        px,py,pz,phid,pphi = readParticle(i)
        for j in range(0,len(px)):
            Points.InsertNextPoint(px[j],py[j],pz[j])
            id_array.InsertNextTuple1(phid[j])
            phi_array.InsertNextTuple1(pphi[j])
    
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(Points)
    polydata.GetPointData().AddArray(id_array)
    polydata.GetPointData().AddArray(phi_array)
    
    if vtk.VTK_MAJOR_VERSION <= 5:
        polydata.Update()
        
    
    outputFile = "/home/subhashis/VisData/merger_trees/particleTrail/time" + str(ts) + ".vtp" 
        
    writer = vtk.vtkXMLPolyDataWriter();
    writer.SetFileName(outputFile);
    if vtk.VTK_MAJOR_VERSION <= 5:
        writer.SetInput(polydata)
    else:
        writer.SetInputData(polydata)
    writer.Write()
    print "Done generating output for time %d" %ts
       


#px,py,pz,phid,pphi = readParticle(4)

#print px.shape
for k in range(0,89):
    createTrail(k)

    
    


