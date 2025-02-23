# Install PyWake if needed
import py_wake

import numpy as np
import matplotlib.pyplot as plt
import random as rnd

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
# Specify random offsets
# =======================================

# random - número entre 0 e 1
numero_aleatorio_rox_1 = (rnd.random()*2)-1
numero_aleatorio_rox_2 = (rnd.random()*2)-1
numero_aleatorio_rox_3 = (rnd.random()*2)-1
numero_aleatorio_rox_3 = (rnd.random()*2)-1
numero_aleatorio_rox_4 = (rnd.random()*2)-1
numero_aleatorio_rox_5 = (rnd.random()*2)-1
numero_aleatorio_rox_6 = (rnd.random()*2)-1
numero_aleatorio_rox_7 = (rnd.random()*2)-1
numero_aleatorio_rox_8 = (rnd.random()*2)-1
numero_aleatorio_roy_1 = (rnd.random()*2)-1
numero_aleatorio_roy_2 = (rnd.random()*2)-1
numero_aleatorio_roy_3 = (rnd.random()*2)-1
numero_aleatorio_roy_3 = (rnd.random()*2)-1
numero_aleatorio_roy_4 = (rnd.random()*2)-1
numero_aleatorio_roy_5 = (rnd.random()*2)-1
numero_aleatorio_roy_6 = (rnd.random()*2)-1
numero_aleatorio_roy_7 = (rnd.random()*2)-1
numero_aleatorio_roy_8 = (rnd.random()*2)-1
numero_aleatorio_cox_1 = (rnd.random()*2)-1
numero_aleatorio_cox_2 = (rnd.random()*2)-1
numero_aleatorio_cox_3 = (rnd.random()*2)-1
numero_aleatorio_cox_4 = (rnd.random()*2)-1
numero_aleatorio_cox_5 = (rnd.random()*2)-1
numero_aleatorio_cox_6 = (rnd.random()*2)-1
numero_aleatorio_cox_7 = (rnd.random()*2)-1
numero_aleatorio_cox_8 = (rnd.random()*2)-1
numero_aleatorio_cox_9 = (rnd.random()*2)-1
numero_aleatorio_cox_10 = (rnd.random()*2)-1
numero_aleatorio_coy_1 = (rnd.random()*2)-1
numero_aleatorio_coy_2 = (rnd.random()*2)-1
numero_aleatorio_coy_3 = (rnd.random()*2)-1
numero_aleatorio_coy_4 = (rnd.random()*2)-1
numero_aleatorio_coy_5 = (rnd.random()*2)-1
numero_aleatorio_coy_6 = (rnd.random()*2)-1
numero_aleatorio_coy_7 = (rnd.random()*2)-1
numero_aleatorio_coy_8 = (rnd.random()*2)-1
numero_aleatorio_coy_9 = (rnd.random()*2)-1
numero_aleatorio_coy_10 = (rnd.random()*2)-1
scale = 250
row_offset_x = (numero_aleatorio_rox_1*scale,numero_aleatorio_rox_2*scale,numero_aleatorio_rox_3*scale,numero_aleatorio_rox_4*scale,numero_aleatorio_rox_5*scale,numero_aleatorio_rox_6*scale,numero_aleatorio_rox_7*scale,numero_aleatorio_rox_8*scale)
row_offset_y = (numero_aleatorio_roy_1*scale,numero_aleatorio_roy_2*scale,numero_aleatorio_roy_3*scale,numero_aleatorio_roy_4*scale,numero_aleatorio_roy_5*scale,numero_aleatorio_roy_6*scale,numero_aleatorio_roy_7*scale,numero_aleatorio_roy_8*scale)
col_offset_x = (numero_aleatorio_cox_1*scale,numero_aleatorio_cox_2*scale,numero_aleatorio_cox_3*scale,numero_aleatorio_cox_4*scale,numero_aleatorio_cox_5*scale,numero_aleatorio_cox_6*scale,numero_aleatorio_cox_7*scale,numero_aleatorio_cox_8*scale,numero_aleatorio_cox_9*scale, numero_aleatorio_cox_10*scale)
col_offset_y = (numero_aleatorio_coy_1*scale,numero_aleatorio_coy_2*scale,numero_aleatorio_coy_3*scale,numero_aleatorio_coy_4*scale,numero_aleatorio_coy_5*scale,numero_aleatorio_coy_6*scale,numero_aleatorio_coy_7*scale,numero_aleatorio_coy_8*scale,numero_aleatorio_coy_9*scale, numero_aleatorio_coy_10*scale)

add_offset_plot_and_print(row_offset_x, row_offset_y, col_offset_x, col_offset_y)