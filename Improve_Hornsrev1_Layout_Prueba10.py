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

posit_x_max = 429492
posit_y_max = 6151447
posit_x_min = 423974
posit_y_min = 6147556

# Original AEP
aep_ref = windFarmModel(wt_x,wt_y).aep().sum()
aep_max = 0
print ('Original AEP: %f GWh'%aep_ref)

# Define the problem constants
WT_Num = 80
farm_widht = posit_x_max - posit_x_min
Farm_height = posit_y_max - posit_y_min
WT_Rad = 40
Gen_Num = 5

def fitness(s):
    penalty = 0
    # Calculate the total energy output and penalize overlapping turbines
    for i in range(0,(WT_Num-1)):
        x1,y1 = s[0][i], s[1][i]
        for j in range(i+1, (WT_Num-1)):
            x2,y2 = s[0][j], s[1][j]
            distance = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
            if distance > 5 * WT_Rad:
                penalty += 0
            else:
                penalty += -(distance - 5 * WT_Rad)
    optim_func = float(windFarmModel(s[0],s[1]).aep().sum()) - penalty # Simplified energy calculation
    return optim_func

# Generate Solutions
solutions = []
for s in range(100):
    x_sol = []
    y_sol = []
    x_rnd = []
    y_rnd = []
    for i in range(0,(WT_Num-1)):
        x_rnd.append(rnd.uniform(-25,+25))
        y_rnd.append(rnd.uniform(-25,+25))
    x_sol = wt_x + x_rnd
    y_sol = wt_y + y_rnd
    solutions.append( (x_sol,y_sol) )
for i in range(Gen_Num):
    rankedsolutions = []
    for s in solutions:
        rankedsolutions.append( (fitness(s),s) )
    rankedsolutions.sort(key=lambda a: a[0])   
    bestsolution = rankedsolutions[:1]
    print("Best solution:")
    print(bestsolution)
    
    NewGen = []
    for _ in range(100):
        NewGen_x = []
        NewGen_y = []
        for i in range(0,(WT_Num-1)):
            x_elements = []
            y_elements = [] 
            for s in bestsolution:
                new_gen_x = 0
                new_gen_y = 0
                good_new_gen_x = False
                good_new_gen_y = False
                while good_new_gen_x == False:
                    new_gen_x = (s[1][0][i]*rnd.uniform(0.99,1.01))
                    if ((new_gen_x <= posit_x_max) and (new_gen_x >= posit_x_min)):
                        good_new_gen_x = True
                while good_new_gen_y == False:
                    new_gen_y = (s[1][1][i]*rnd.uniform(0.99,1.01))
                    if ((new_gen_y <= posit_y_max) and (new_gen_y >= posit_y_min)):
                        good_new_gen_y = True
                x_elements.append(new_gen_x)
                y_elements.append(new_gen_y)
            NewGen_x.append((sum(x_elements)/len(x_elements)))
            NewGen_y.append((sum(y_elements)/len(y_elements)))
        NewGen.append( (NewGen_x,NewGen_y) )
    solutions = NewGen

finalsolution = []
for s in bestsolution:
    finalsolution.append( (s[1][0],s[1][1]) )
finalsolution = finalsolution[:1]
print("finalsolution")
print(finalsolution)
x_max,y_max = (finalsolution[0][0],finalsolution[0][1])
aep_ref = float(windFarmModel(wt_x,wt_y).aep().sum())
print("aep_ref:",aep_ref)
aep_max = float(windFarmModel(x_max,y_max).aep().sum())
print("aep_max:",aep_max)

plt.figure()
plt.plot(wt_x, wt_y,'b.')
wt.plot(x_max, y_max)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.show()

         
    