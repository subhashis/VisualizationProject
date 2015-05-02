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

noh = len(id)    
print len(id)   
#print max(pid)
#print min(pid)
#print len(Tree_root_ID)
#print max(Tree_root_ID)
#print min(Tree_root_ID)
#maxMass = max(mvir)
#
#for i in range(0,len(x)):
#    if maxMass == mvir[i]:
#        print id[i]
num = Snap_num.astype(int)
print max(num)
print len(num)
binTemp = range(0,90)


hist, bin_edges = np.histogram(num, bins=binTemp)

print hist
print bin_edges

print hist.shape
print bin_edges.shape

count = 0;

for i in range(0,noh):
    if num[i] == 87:
        count += 1
print "count = %d" %count
print "total = %d" %np.sum(hist)
print hist[num[2]]