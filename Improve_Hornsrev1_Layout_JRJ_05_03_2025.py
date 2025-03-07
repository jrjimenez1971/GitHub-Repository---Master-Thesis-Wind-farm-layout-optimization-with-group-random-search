# Juan R Jimenez DTU Wind Energy Master
# Master Thesis Wind farm layout optimization with group random search

# Install PyWake if needed
import py_wake
import numpy as np
import matplotlib.pyplot as plt
import random as rnd

#importing the properties of Hornsrev1, which are already stored in PyWake
from py_wake.examples.data.hornsrev1 import V80
from py_wake.examples.data.hornsrev1 import Hornsrev1Site
from py_wake.examples.data.iea37 import IEA37Site
from py_wake.site import UniformWeibullSite

# Alternative sites that can be selected by the user, besides Hornsrev1Site
sites = {"IEA37": IEA37Site(n_wt=16),
         "Hornsrev1": Hornsrev1Site(),
         "UniformSite": UniformWeibullSite(p_wd = [.20,.25,.35,.25], a = [9.176929,  9.782334,  9.531809,  9.909545], k = [2.392578, 2.447266, 2.412109, 2.591797], ti = 0.1)}

from py_wake.examples.data.hornsrev1 import wt_x, wt_y

# BastankhahGaussian combines the engineering wind farm model, `PropagateDownwind` with
# the `BastankhahGaussianDeficit` wake deficit model and the `SquaredSum` super position model
from py_wake.literature.gaussian_models import Bastankhah_PorteAgel_2014

def Data_Imputs ():
    global site
    global Site_Elect
    global Iter_Lenght
    global Iter_Num
    global Gen_Num
    Site_Elec = int(input ("Enter Site to be used 1-IEA37, 2-Hornsrev1, 3-UniforSite: "))
    Iter_Lenght = int(input ("Enter Iteration Lenght in meters: "))
    Iter_Num = int(input ("Enter Number of Iterations on each Generation: ")) 
    Gen_Num = int(input ("Enter Number of Genetic Generations: "))
    if Site_Elec == 1:
        site = sites["IEA37"]
    elif Site_Elec == 2:
        site = sites["Hornsrev1"]              
    else:
        site = sites["UniformSite"]        
    return site
    return Iter_Lenght
    return Iter_Num
    return Gen_Num

Data_Imputs ()

# After we import the objects we instatiate them:
# site = Hornsrev1Site()
wt = V80()
windFarmModel = Bastankhah_PorteAgel_2014(site, wt, k=0.0324555)

site.plot_wd_distribution(n_wd=12)

posit_x_max = 429492 + 25
posit_y_max = 6151447 + 25
posit_x_min = 423974 - 25
posit_y_min = 6147556 - 25

# Original AEP
aep_ref = windFarmModel(wt_x,wt_y).aep().sum()
aep_max = 0
print ('Original AEP: %f GWh'%aep_ref)

# Define the problem constants
WT_Num = 80
farm_widht = posit_x_max - posit_x_min
Farm_height = posit_y_max - posit_y_min
WT_Rad = 40
# Gen_Num = 10

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

def list_randon_values(lst, val_max, val_min):
    lst_rnd = []
    for x in lst:
        lst_rnd_val = 0
        good_lst_rnd_val = False 
        while good_lst_rnd_val == False:
            lst_rnd_val = (x + (Iter_Lenght * rnd.uniform(-1,1)))
            if (lst_rnd_val < val_max) and (lst_rnd_val > val_min):
                good_lst_rnd_val = True
        lst_rnd.append(lst_rnd_val)
    return lst_rnd

solutions = []
for s in range(Iter_Num):
    x_rnd = []
    y_rnd = []
    x_sol = wt_x
    x_rnd = list_randon_values(x_sol, posit_x_max, posit_x_min)
    x_sol = np.array(x_rnd)
    y_sol = wt_y
    y_rnd = list_randon_values(y_sol, posit_y_max, posit_y_min)
    y_sol = np.array(y_rnd)   
    solutions.append( (x_sol,y_sol) )

rankedsolutions = []
for s in solutions:
    rankedsolutions.append( (fitness(s),s) )
rankedsolutions.sort(key=lambda a: a[0])   
bestsolution = rankedsolutions[(Iter_Num-1)]
bestsolution = ( (bestsolution[1][0],bestsolution[1][1]) )
print("Best solution on reference layout randomized Generation 0:")
print(bestsolution)
x_best_sol,y_best_sol = (bestsolution[0],bestsolution[1])
aep_best_sol = float(windFarmModel(x_best_sol,y_best_sol).aep().sum())
print("aep_best_sol on reference layout randomized Generation 0:",aep_best_sol)

global Gen
Gen = []
NewGen_i_best = []
aep_best_sol_i = []
Delta_aep_best_sol_i = []
for i in range(1,Gen_Num+1):
    Gen.append(i)
    NewGen = []
    for _ in range(Iter_Num):
        new_gen_x = []
        new_gen_y = []
        for s in bestsolution:
            new_gen_x = bestsolution[0]
            new_gen_x_rnd = []
            new_gen_x_rnd = list_randon_values(new_gen_x, posit_x_max , posit_x_min)
        for s in bestsolution:
            new_gen_y = bestsolution[1]
            new_gen_y_rnd = []
            new_gen_y_rnd = list_randon_values(new_gen_y, posit_y_max, posit_y_min)
        x_sol = new_gen_x_rnd
        y_sol = new_gen_y_rnd
        NewGen.append( (x_sol,y_sol) )

    solutions = NewGen
    rankedsolutions = []
    for s in solutions:
        rankedsolutions.append( (fitness(s),s) )
    rankedsolutions.sort(key=lambda a: a[0])   
    bestsolution = rankedsolutions[(Iter_Num-1)]
    bestsolution = ( (bestsolution[1][0],bestsolution[1][1]) )
    print(f"Best solution on Genetic Generation {i}:")
    print(bestsolution)
    x_best_sol,y_best_sol = (bestsolution[0],bestsolution[1])
    NewGen_i_best.append( (x_best_sol,y_best_sol) )
    aep_best_sol = float(windFarmModel(x_best_sol,y_best_sol).aep().sum())
    aep_best_sol_i.append(aep_best_sol)
    print(f"aep_best_sol on Genetic Generation {i}:",aep_best_sol)

aep_best_sol_prev = aep_ref
for s in aep_best_sol_i:
    Delta_aep_best_sol_i_elem = s - aep_best_sol_prev
    aep_best_sol_prev = s
    Delta_aep_best_sol_i.append(Delta_aep_best_sol_i_elem)
    
Delta_aep_best_sol_i
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
plt.title('Original blue and evolutions of Generation Layouts yellow')
plt.plot(wt_x, wt_y,'b.')
for s in NewGen_i_best:
    plt.plot(s[0], s[1],'y.')     
plt.plot(x_max, y_max, 'r.')
plt.xlabel('Evolution of positions x [m]')
plt.ylabel('Evolution of positions y [m]')
plt.show()

plt.figure()
plt.title('Original blue and Final black Layouts')
plt.plot(wt_x, wt_y,'b.')
wt.plot(x_max, y_max)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.show()

plt.figure()
plt.title('Evolution of improvement throught genetic evolutions')
plt.plot(Gen, aep_best_sol_i, 'r.')
plt.xlabel('Genetic Generation i')
plt.ylabel('Best aep Generation i GWh')
plt.show()

plt.figure()
plt.title('Evolution of diferencial improvement throught genetic evolutions')
plt.plot(Gen, Delta_aep_best_sol_i, 'r.')
plt.xlabel('Genetic Generation i')
plt.ylabel('Best Delta aep Generation i GWh')
plt.show()
    