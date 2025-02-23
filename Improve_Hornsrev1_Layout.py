# Install PyWake if needed
import py_wake

import numpy as np
import matplotlib.pyplot as plt

#importing the properties of Hornsrev1, which are already stored in PyWake
from py_wake.examples.data.hornsrev1 import V80
from py_wake.examples.data.hornsrev1 import Hornsrev1Site
from py_wake.examples.data.hornsrev1 import wt_x, wt_y

# BastankhahGaussian combines the engineering wind farm model, `PropagateDownwind` with
# the `BastankhahGaussianDeficit` wake deficit model and the `SquaredSum` super position model
from py_wake.literature.gaussian_models import Bastankhah_PorteAgel_2014

# After we import the objects we instatiate them:
site = Hornsrev1Site()
wt = V80()
windFarmModel = Bastankhah_PorteAgel_2014(site, wt, k=0.0324555)

site.plot_wd_distribution(n_wd=12)

# Original layout
plt.figure()
wt.plot(wt_x, wt_y)
plt.xlabel('x [m]')
plt.ylabel('y [m]')

# Original AEP
aep_ref = windFarmModel(wt_x,wt_y).aep().sum()
print ('Original AEP: %f GWh'%aep_ref)

# Here we define a function to print and plot your new layout. No need to change anything here
def add_offset_plot_and_print(row_offset_x, row_offset_y, col_offset_x, col_offset_y):
    x,y = wt_x, wt_y
    y = np.reshape(y,(10,8)).astype(float)
    x = np.reshape(x,(10,8)).astype(float)

    x+= np.array(row_offset_x)
    y+= np.array(row_offset_y)
    x+= np.array(col_offset_x)[:,np.newaxis]
    y+= np.array(col_offset_y)[:,np.newaxis]
    y = np.maximum(min(wt_y), np.minimum(max(wt_y), y.flatten()))
    x = np.maximum(min(wt_x), np.minimum(max(wt_x), x.flatten()))

    plt.plot()
    plt.plot(wt_x, wt_y,'b.')
    wt.plot(x, y)
    aep = windFarmModel(x,y).aep().sum()
    print ("AEP ref", aep_ref.values)
    print ("AEP", aep.values)
    print ("Increase: %f %%"%((aep-aep_ref)/aep_ref*100))
    
# =======================================
# Specify offsets
# =======================================

row_offset_x = np.linspace(0,1,8)* -500
row_offset_y = np.linspace(0,1,8) * 0

col_offset_x = np.linspace(0,1,10) * 0
col_offset_y = np.linspace(0,1,10) * 0

add_offset_plot_and_print(row_offset_x, row_offset_y, col_offset_x, col_offset_y)