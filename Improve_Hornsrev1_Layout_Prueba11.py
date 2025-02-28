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

posit_x_max = 429492 + 300
posit_y_max = 6151447 + 300
posit_x_min = 423974 - 300
posit_y_min = 6147556 - 300

# Original AEP
aep_ref = windFarmModel(wt_x,wt_y).aep().sum()
aep_max = 0
print ('Original AEP: %f GWh'%aep_ref)

# Define the problem constants
WT_Num = 80
farm_widht = posit_x_max - posit_x_min
Farm_height = posit_y_max - posit_y_min
WT_Rad = 40
Gen_Num = 10

def fitness(s):
    penalty = 0
    # Calculate the total energy output and penalize overlapping turbines
    for i in range(0,(WT_Num-1)):
        x1,y1 = s[0][i], s[1][i]
        for j in range(i+1,(WT_Num-1)):
            x2,y2 = s[0][j], s[1][j]
            distance = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
            if distance > 5 * WT_Rad:
                penalty += 0
            else:
                penalty += -(distance - 5 * WT_Rad)
    optim_func = float(windFarmModel(s[0],s[1]).aep().sum()) - penalty # Simplified energy calculation
    return optim_func

# Generate Solutions
def all_elements_greater_than(lst, value):
    return all(x > value for x in lst)
def all_elements_smaller_than(lst, value):
    return all(x < value for x in lst)

def list_randon_values(lst):
    lst_rnd = []
    for x in lst:
        lst_rnd.append(x * rnd.uniform(0.99999,1.00001))
    return lst_rnd
    
solutions = []
for s in range(200):
    x_rnd = []
    y_rnd = []
    good_x_sol = False
    while good_x_sol == False:
        x_sol = wt_x
        x_rnd = []
        for i in range(0,(WT_Num)):
            x_rnd.append(rnd.uniform(-10,+10))
        x_sol += np.array(x_rnd)
        if (all_elements_smaller_than(x_sol, posit_x_max)) and (all_elements_greater_than(x_sol, posit_x_min)):
            good_x_sol = True
    good_y_sol = False
    while good_y_sol == False:
        y_sol = wt_y
        y_rnd = []
        for i in range(0,(WT_Num)):
            y_rnd.append(rnd.uniform(-25,+25))
        y_sol += np.array(y_rnd)   
        if (all_elements_smaller_than(y_sol, posit_y_max)) and (all_elements_greater_than(y_sol, posit_y_min)):
            good_y_sol = True
    solutions.append( (x_sol,y_sol) )

print("First solutions:")
print(solutions[:3])

for i in range(Gen_Num):
    rankedsolutions = []
    for s in solutions:
        rankedsolutions.append( (fitness(s),s) )
    rankedsolutions.sort(key=lambda a: a[0])   
    bestsolution = rankedsolutions[199]
    print("Best solution:")
    print(bestsolution)
    bestsolution = ( (bestsolution[1][0],bestsolution[1][1]) )
    print("Best solution:")
    print(bestsolution)
    x_best_sol,y_best_sol = (bestsolution[0],bestsolution[1])
    aep_best_sol = float(windFarmModel(x_best_sol,y_best_sol).aep().sum())
    print(f"aep_best_sol {i}:",aep_best_sol)

    NewGen = []
    for _ in range(200):
        new_gen_x = []
        new_gen_y = []
        good_new_gen_x = False
        good_new_gen_y = False
        while good_new_gen_x == False:
            for s in bestsolution:
                new_gen_x = bestsolution[0]
                new_gen_x_rnd = []
                new_gen_x_rnd = list_randon_values(new_gen_x)
                if (all_elements_smaller_than(new_gen_x_rnd, posit_x_max)) and (all_elements_greater_than(new_gen_x_rnd, posit_x_min)):
                    good_new_gen_x = True
        while good_new_gen_y == False:
            for s in bestsolution:
                new_gen_y = bestsolution[1]
                new_gen_y_rnd = []
                new_gen_y_rnd = list_randon_values(new_gen_y)
                if (all_elements_smaller_than(new_gen_y_rnd, posit_y_max)) and (all_elements_greater_than(new_gen_y_rnd, posit_y_min)):
                    good_new_gen_y = True
        x_sol = new_gen_x_rnd
        y_sol = new_gen_y_rnd
        NewGen.append( (x_sol,y_sol) )
    solutions = NewGen
    print(f"New Gen {i}")
    print(NewGen[:4])

finalsolution = []
finalsolution = bestsolution
print("finalsolution")
print(finalsolution)
x_max,y_max = (finalsolution[0],finalsolution[1])
aep_ref = float(windFarmModel(wt_x,wt_y).aep().sum())
print("aep_ref:",aep_ref)
aep_max = float(windFarmModel(x_max,y_max).aep().sum())
print("aep_max:",aep_max)

plt.figure()
plt.plot(wt_x, wt_y,'b.')
plt.plot(x_sol, y_sol,'r.')
wt.plot(x_max, y_max)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.show()

         
    