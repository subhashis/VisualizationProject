# -*- coding: utf-8 -*-
import numpy as np

timeStep = 0

if timeStep == 88:
    suffix = "1.0000"
else:
    i = timeStep + 12
    suffix = "0." + str(i) + "00"

path = "/home/subhashis/VisData/contestData/2015/ds14_scivis_0128_e4_dt04_" + suffix

print path


