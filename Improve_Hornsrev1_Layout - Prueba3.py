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
# plt.figure()
# wt.plot(wt_x, wt_y)
# plt.xlabel('x [m]')
# plt.ylabel('y [m]')

x_max,y_max = wt_x, wt_y

# Original AEP
aep_ref = windFarmModel(wt_x,wt_y).aep().sum()
aep_max = 0
print ('Original AEP: %f GWh'%aep_ref)

# Here we define a function to print and plot your new layout. No need to change anything here
def add_offset_plot_and_print(row_offset_x_1, row_offset_y_1, col_offset_x_1, col_offset_y_1,row_offset_x_2, row_offset_y_2, col_offset_x_2, col_offset_y_2):
    global aep_rnd
    global x_rnd,y_rnd
    x,y = wt_x, wt_y
    y = np.reshape(y,(10,8)).astype(float)
    x = np.reshape(x,(10,8)).astype(float)
  
    print (y)    
    print (x)

    x+= np.array(row_offset_x_1)
    y+= np.array(row_offset_y_1)
    x+= np.array(col_offset_x_1)[:,np.newaxis]
    y+= np.array(col_offset_y_1)[:,np.newaxis]

    x+= np.array(row_offset_x_2)
    y+= np.array(row_offset_y_2)
    x+= np.array(col_offset_x_2)[:,np.newaxis]
    y+= np.array(col_offset_y_2)[:,np.newaxis] 

    Delta_x = [500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500]
    Delta_y = [500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500]
    Delta_x_neg = [-500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500]
    Delta_y_neg = [-500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500, -500]
   
    # Delta_y = np.reshape(Delta_y,(10,8)).astype(float)
    # Delta_x = np.reshape(Delta_x,(10,8)).astype(float)

    wt_y_wide = wt_y + Delta_y
    wt_x_wide = wt_x + Delta_x
    wt_y_narrow = wt_y + Delta_y_neg
    wt_x_narrow = wt_x + Delta_x_neg

    y = np.maximum(min(wt_y_narrow), np.minimum(max(wt_y_wide), y.flatten()))
    x = np.maximum(min(wt_x_narrow), np.minimum(max(wt_x_wide), x.flatten()))

    x_rnd,y_rnd = x,y
    aep_rnd = windFarmModel(x,y).aep().sum()
    
# =======================================
# Specify random offsets
# =======================================
scale = 100

# Offset#1:
row_offset_x_1 = np.linspace(0,1.3,8)* (+1.2*scale)
row_offset_y_1 = np.linspace(0,1.3,8) * (-1.2*scale)
col_offset_x_1 = np.linspace(0,1.3,10) * (+1.2*scale)
col_offset_y_1 = np.linspace(0,1.3,10) * (-1.2*scale)

# Offset#2: Random in points located after offset#1 
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

row_offset_x_2 = (numero_aleatorio_rox_1*scale,numero_aleatorio_rox_2*scale,numero_aleatorio_rox_3*scale,numero_aleatorio_rox_4*scale,numero_aleatorio_rox_5*scale,numero_aleatorio_rox_6*scale,numero_aleatorio_rox_7*scale,numero_aleatorio_rox_8*scale)
row_offset_y_2 = (numero_aleatorio_roy_1*scale,numero_aleatorio_roy_2*scale,numero_aleatorio_roy_3*scale,numero_aleatorio_roy_4*scale,numero_aleatorio_roy_5*scale,numero_aleatorio_roy_6*scale,numero_aleatorio_roy_7*scale,numero_aleatorio_roy_8*scale)
col_offset_x_2 = (numero_aleatorio_cox_1*scale,numero_aleatorio_cox_2*scale,numero_aleatorio_cox_3*scale,numero_aleatorio_cox_4*scale,numero_aleatorio_cox_5*scale,numero_aleatorio_cox_6*scale,numero_aleatorio_cox_7*scale,numero_aleatorio_cox_8*scale,numero_aleatorio_cox_9*scale, numero_aleatorio_cox_10*scale)
col_offset_y_2 = (numero_aleatorio_coy_1*scale,numero_aleatorio_coy_2*scale,numero_aleatorio_coy_3*scale,numero_aleatorio_coy_4*scale,numero_aleatorio_coy_5*scale,numero_aleatorio_coy_6*scale,numero_aleatorio_coy_7*scale,numero_aleatorio_coy_8*scale,numero_aleatorio_coy_9*scale, numero_aleatorio_coy_10*scale)

for n in range (1,200):
    add_offset_plot_and_print(row_offset_x_1, row_offset_y_1, col_offset_x_1, col_offset_y_1,row_offset_x_2, row_offset_y_2, col_offset_x_2, col_offset_y_2)
    if aep_rnd > aep_max:
        aep_max = aep_rnd
        x_max,y_max = x_rnd,y_rnd
    n += 1

plt.figure()
plt.plot(wt_x, wt_y,'b.')
wt.plot(x_max, y_max)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.show()

print ("AEP ref", aep_ref.values)
print ("AEP_Max", aep_max.values)
print ("Increase: %f %%"%((aep_max-aep_ref)/aep_ref*100))
    