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

posit_x_max = max(wt_x)
posit_y_max = max(wt_y)
posit_x_min = min(wt_x)
posit_y_min = min(wt_y)

# Original AEP
aep_ref = windFarmModel(wt_x,wt_y).aep().sum()
aep_max = 0
print ('Original AEP: %f GWh'%aep_ref)

# Define the problem constants
WT_Num = 10
farm_widht = posit_x_max - posit_x_min
Farm_height = posit_y_max - posit_y_min
WT_Rad = 50

def fitness(s):
    penalty = 0
    # Calculate the total energy output and penalize overlapping turbines
    for i in range(0,9):
        x1,y1 = s[0][i], s[1][i]
        for j in range(i+1, 9):
            x2, y2 = s[0][j], s[1][j]
            distance = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
            if distance > 5 * WT_Rad:
                penalty += 0
            else:
                penalty += -(distance - 5 * WT_Rad)
    optim_func = float(windFarmModel(s[0],s[1]).aep().sum()) - penalty # Simplified energy calculation
    return optim_func

# Generate Solutions
solutions = []
for s in range(1000):
    solutions.append( ((rnd.uniform(posit_x_max,posit_x_min),rnd.uniform(posit_x_max,posit_x_min),rnd.uniform(posit_x_max,posit_x_min),rnd.uniform(posit_x_max,posit_x_min),rnd.uniform(posit_x_max,posit_x_min),rnd.uniform(posit_x_max,posit_x_min),rnd.uniform(posit_x_max,posit_x_min),rnd.uniform(posit_x_max,posit_x_min),rnd.uniform(posit_x_max,posit_x_min),rnd.uniform(posit_x_max,posit_x_min)),(rnd.uniform(posit_y_max,posit_y_min),rnd.uniform(posit_y_max,posit_y_min),rnd.uniform(posit_y_max,posit_y_min),rnd.uniform(posit_y_max,posit_y_min),rnd.uniform(posit_y_max,posit_y_min),rnd.uniform(posit_y_max,posit_y_min),rnd.uniform(posit_y_max,posit_y_min),rnd.uniform(posit_y_max,posit_y_min),rnd.uniform(posit_y_max,posit_y_min),rnd.uniform(posit_y_max,posit_y_min))) )

global bestsoluions
bestsolutions =  []
for i in range(10):
    rankedsolutions = []
    for s in solutions:
        rankedsolutions.append( (fitness(s),s) )
        rankedsolutions.sort()
        rankedsolutions.reverse()    
    bestsolutions = rankedsolutions[:15]
    print("Best solution:")
    print(bestsolutions[:5])
    
    elements = []
    e10 = []
    e11 = []
    e12 = []
    e13 = []
    e14 = []
    e15 = []
    e16 = []
    e17 = []
    e18 = []
    e19 = []
    e20 = []
    e21 = []
    e22 = []
    e23 = []
    e24 = []
    e25 = [] 
    e26 = []
    e27 = []
    e28 = []
    e29 = []
    for s in bestsolutions:
        e10.append(s[1][0][0]*rnd.uniform(0.99,1.01))
        e11.append(s[1][0][1]*rnd.uniform(0.99,1.01))
        e12.append(s[1][0][2]*rnd.uniform(0.99,1.01))
        e13.append(s[1][0][3]*rnd.uniform(0.99,1.01))
        e14.append(s[1][0][4]*rnd.uniform(0.99,1.01))
        e15.append(s[1][0][5]*rnd.uniform(0.99,1.01))
        e16.append(s[1][0][6]*rnd.uniform(0.99,1.01))
        e17.append(s[1][0][7]*rnd.uniform(0.99,1.01))
        e18.append(s[1][0][8]*rnd.uniform(0.99,1.01))
        e19.append(s[1][0][9]*rnd.uniform(0.99,1.01))
        e20.append(s[1][1][0]*rnd.uniform(0.99,1.01))
        e21.append(s[1][1][1]*rnd.uniform(0.99,1.01))
        e22.append(s[1][1][2]*rnd.uniform(0.99,1.01))
        e23.append(s[1][1][3]*rnd.uniform(0.99,1.01))
        e24.append(s[1][1][4]*rnd.uniform(0.99,1.01))
        e25.append(s[1][1][5]*rnd.uniform(0.99,1.01))
        e26.append(s[1][1][6]*rnd.uniform(0.99,1.01))
        e27.append(s[1][1][7]*rnd.uniform(0.99,1.01))
        e28.append(s[1][1][8]*rnd.uniform(0.99,1.01))
        e29.append(s[1][1][9]*rnd.uniform(0.99,1.01))
    NewGen = []
    for _ in range(1000):
        NewGen.append( ((rnd.choice(e10),rnd.choice(e11),rnd.choice(e12),rnd.choice(e13),rnd.choice(e14),rnd.choice(e15),rnd.choice(e16),rnd.choice(e17),rnd.choice(e18),rnd.choice(e19)),(rnd.choice(e20),rnd.choice(e21),rnd.choice(e22),rnd.choice(e23),rnd.choice(e24),rnd.choice(e25),rnd.choice(e26),rnd.choice(e27),rnd.choice(e28),rnd.choice(e29))) )
    solutions = NewGen

finalsolution = []
for s in bestsolutions:
    finalsolution.append( (s[1][0],s[1][1]) )
print("finalsolution")
print(finalsolution)
finalsolution = finalsolution[:1]
print("finalsolution")
print(finalsolution)
x_max,y_max = (finalsolution[0][0],finalsolution[0][1])
aep_ref = windFarmModel(wt_x,wt_y).aep().sum()
print("aep_ref:",aep_ref)
aep_max = float(windFarmModel(x_max,y_max).aep().sum())
print("aep_max:",aep_max)

plt.figure()
# plt.plot(wt_x, wt_y,'b.')
wt.plot(x_max, y_max)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.show()

         
    