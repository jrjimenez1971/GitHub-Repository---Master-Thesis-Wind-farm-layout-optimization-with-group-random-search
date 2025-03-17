# Juan R Jimenez DTU Wind Energy Master
# Master Thesis Wind farm layout optimization with group random search

# Glabal variables definition and initiation
global bestsolution
global solutions_dic_var_iter_length_x
global bestfitness_dic_var_iter_length_x
global aep_best_sol_iter_length_x
global best_length
    
global genbestsolution
global Gen
global NewGen_i_best
global aep_best_gensol_i
global Delta_aep_best_sol_i
global Delta_aep_best_sol_i_elem
global gensolutions_dic_var_iter_length_x
global genbestfitness_dic_var_iter_length_x
global aep_best_gensol_iter_length_x
global genbest_length_i
global Gen_Iter_Lenght_Pack_x
global Gen_Generations_Pack_i
global aep_best_gensol_calc_i
    
Gen = []
NewGen_i_best = []
aep_best_sol_i = []
Delta_aep_best_sol_i = []
Delta_aep_best_sol_i_elem = []
genbest_length_i = []
Gen_Iter_Lenght_Pack_x = []
Gen_Generations_Pack_i = []
aep_best_gensol_calc_i = []
solutions_dic_var_iter_length_x = []

# Install PyWake if needed
import pandas as pd
import py_wake
import numpy as np
import matplotlib.pyplot as plt
import random as rnd
from matplotlib.path import Path

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
    global Iter_Length
    global Iter_Num
    global Gen_Num
    global Border_Margin
    Site_Elec = int(input ("Enter Site to be used 1-IEA37, 2-Hornsrev1, 3-UniforSite: "))
    Iter_Length = int(input ("Enter Iteration Lenght in meters: "))
    Iter_Num = int(input ("Enter Number of Iterations on each Generation: ")) 
    Gen_Num = int(input ("Enter Number of Genetic Generations: "))
    Border_Margin = int(input ("Enter border aditional distance from actual laylout in m: "))
    if Site_Elec == 1:
        site = sites["IEA37"]
    elif Site_Elec == 2:
        site = sites["Hornsrev1"]              
    else:
        site = sites["UniformSite"] 
    return site
    return Iter_Length
    return Iter_Num
    return Gen_Num
    return Border_Margin

Data_Imputs ()

iter_lengths = {"0":Iter_Length,"1":Iter_Length*0.7, "2":Iter_Length*1.3}

# After we import the objects we instatiate them:
wt = V80()
windFarmModel = Bastankhah_PorteAgel_2014(site, wt, k=0.0324555)

site.plot_wd_distribution(n_wd=12)

border_points = [((wt_x[0]-Border_Margin), (wt_y[0]+Border_Margin)), ((wt_x[7]-Border_Margin), (wt_y[7]-Border_Margin)), ((wt_x[79]+Border_Margin), (wt_y[79]-Border_Margin)), ((wt_x[72]+Border_Margin, wt_y[72]+Border_Margin)), ((wt_x[0]-Border_Margin, wt_y[0]+Border_Margin))]

# Creation of Mesh Border to Check all New Positions are whithin 
def create_mesh_border(points):
    codes = [Path.MOVETO] + [Path.LINETO] * (len(points) - 1)
    return Path(points, codes)

mesh_border = create_mesh_border(border_points)

# Original AEP
aep_ref = round(float(windFarmModel(wt_x,wt_y).aep().sum()),4)
aep_max = 0
print ("Original AEP: ",aep_ref,"GWh")

# Define the problem constants of Hornsrev1 windfarm
WT_Num = 80
WT_Rad = 40

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
                penalty += -10*(distance - 5 * WT_Rad)
    optim_func = float(windFarmModel(s[0],s[1]).aep().sum()) - penalty # Simplified energy calculation
    return optim_func

def filter_points_inside_border(border_path,point):
    if border_path.contains_point((point[0],point[1])):
        return True
    else:
        return False

def convertir_listas_a_lista_de_listas(lista_x, lista_y):
    if len(lista_x) != len(lista_y):
        raise ValueError("Las listas deben tener la misma longitud")
    tupple_de_listas = []
    for i in range(len(lista_x)):
        tupple_de_listas.append( (lista_x[i], lista_y[i]) )
    return tupple_de_listas

def list_randon_values(lst1,lst2,Length,border):
    lst1_rnd = []
    lst1_rnd_val = 0
    lst2_rnd = []
    lst2_rnd_val = 0
    lst = convertir_listas_a_lista_de_listas(lst1, lst2)
    for s in lst:
        good_lst_rnd_val = False
        sx, sy = [s[0],s[1]]
        while good_lst_rnd_val == False:
            lst1_rnd_val = sx + (Length * rnd.uniform(-1,1))
            lst2_rnd_val = sy + (Length * rnd.uniform(-1,1)) 
            new_point = (lst1_rnd_val,lst2_rnd_val)
            if filter_points_inside_border(border, new_point) == True:
                good_lst_rnd_val = True              
        lst1_rnd.append(lst1_rnd_val)
        lst2_rnd.append(lst2_rnd_val)      
    return lst1_rnd,lst2_rnd

def getSolutions():
    solutions_dic_var_iter_length_x = {"0":[],"1":[], "2":[]}
    bestfitness_dic_var_iter_length_x = {"0":[],"1":[], "2":[]}
    aep_best_sol_iter_length_x = {"0":[],"1":[], "2":[]}
    rankedsolutions_x = {"0":[],"1":[], "2":[]}
    for x in iter_lengths.keys():
        Iter_Length = iter_lengths[x]
        solutions = []
        for s in range(Iter_Num):
            x_sol,y_sol = list_randon_values(wt_x,wt_y,Iter_Length,mesh_border)
            solutions.append( (x_sol,y_sol) )    
        bestsolution = ranksolutions(solutions)
        rnk = ranksolutions(solutions)     
        bestsolution = rnk[0]
        rankedsoltions = rnk[1]
        bestfitness_dic_var_iter_length_x[x] = bestsolution[0]
        bestsolution = ( (bestsolution[1][0],bestsolution[1][1]) )
        solutions_dic_var_iter_length_x[x] = bestsolution  
        aep_best_sol_iter_length_x[x] = float(windFarmModel(bestsolution[0],bestsolution[1]).aep().sum())
        rankedsolutions_x[x] = rankedsoltions
    best_length = max( aep_best_sol_iter_length_x, key= aep_best_sol_iter_length_x.get)
    bestsolution = solutions_dic_var_iter_length_x[best_length] 
    x_best_sol,y_best_sol = (bestsolution[0],bestsolution[1])
    aep_best_sol_calc = float(windFarmModel(x_best_sol,y_best_sol).aep().sum())
    
    if aep_best_sol_calc < aep_ref:
        bestsolution = ( (wt_x,wt_y) )
        aep_best_sol = round(aep_ref,4)
    else:
        aep_best_sol = round(aep_best_sol_calc,4)
        bestsolution =  x_best_sol,y_best_sol 
    
    print("Best solution on reference layout randomized Generation 0:")
    print(np.round(bestsolution,4))
    print("aep_best_sol on reference layout randomized Generation 0:",aep_best_sol,"GWh")
    return bestsolution, aep_best_sol,solutions_dic_var_iter_length_x,aep_best_sol_iter_length_x,best_length,rankedsolutions_x
    
def getGenSolutions(imputsolution,aep_best_sol_prev):
    NewGen = imputsolution
    NewGen_i_best = []
    NewGen_i_best_prev = []
    aep_best_gensol_i = []
    genrankedsolutions_i = []
    aep_best_gensol_prev = 0
    for i in range(1,Gen_Num+1):
        Gen.append(i)
        genbest_length = []
        Gen_Iter_Lenght_Pack_x = []
        gensolutions_dic_var_iter_length_x = {"0":[],"1":[], "2":[]}
        genbestfitness_dic_var_iter_length_x = {"0":[],"1":[], "2":[]}
        aep_best_gensol_iter_length_x = {"0":[],"1":[], "2":[]}
        genrankedsolutions_x = {"0":[],"1":[], "2":[]}
        for x in iter_lengths.keys():
            Iter_Length = iter_lengths[x]
            gensolutions = []
            for _ in range(Iter_Num):
                new_gen_x = []
                new_gen_y = []
                for s in NewGen:
                    new_gen_x = NewGen[0]
                    new_gen_y = NewGen[1]
                    new_gen_x_rnd = []
                    new_gen_y_rnd = []
                    new_gen_x_rnd,new_gen_y_rnd = list_randon_values(new_gen_x,new_gen_y,Iter_Length,mesh_border)
                x_sol = new_gen_x_rnd
                y_sol = new_gen_y_rnd                
                gensolutions.append( (x_sol,y_sol) )
            rnk = ranksolutions(gensolutions)     
            genbestsolution = rnk[0]
            genrankedsoltions = rnk[1]
            genbestsolution = ( (genbestsolution[1][0],genbestsolution[1][1]) )
            gensolutions_dic_var_iter_length_x[x] = genbestsolution
            aep_best_gensol_iter_length_x[x] = float(windFarmModel(genbestsolution[0],genbestsolution[1]).aep().sum())
            genbestfitness_dic_var_iter_length_x[x] = genbestsolution
            genrankedsolutions_x[x] = genrankedsoltions
        Gen_Iter_Lenght_Pack_x.append((genbestfitness_dic_var_iter_length_x, aep_best_gensol_iter_length_x))
        genbest_length = max(aep_best_gensol_iter_length_x, key=aep_best_gensol_iter_length_x.get)
        genbestsolution = gensolutions_dic_var_iter_length_x[genbest_length]          
        aep_best_gensol_calc = float(windFarmModel(genbestsolution[0],genbestsolution[1]).aep().sum())
        Gen_Generations_Pack_i.append(Gen_Iter_Lenght_Pack_x)
        aep_best_gensol_calc_i.append(aep_best_gensol_calc)
        genbest_length_i.append(genbest_length)
        genrankedsolutions_i.append(genrankedsolutions_x)
        
        if i == 1:
            if aep_best_gensol_calc <= aep_best_sol_prev:
                genbestsolution = bestsolution_zero
                aep_best_gensol = round(aep_best_sol_prev,4)
            else:
                aep_best_gensol = round(aep_best_gensol_calc,4)
        else:
            if aep_best_gensol_calc <= aep_best_sol_prev:
                genbestsolution = NewGen_i_best_prev
                aep_best_gensol = round(aep_best_gensol_prev,4)
            else:
                aep_best_gensol = round(aep_best_gensol_calc,4)
        
        NewGen = genbestsolution
        NewGen_i_best.append(genbestsolution)
        aep_best_gensol_i.append(aep_best_gensol)
        aep_best_gensol_prev = aep_best_gensol
        NewGen_i_best_prev = NewGen
        print(f"Best solution on Genetic Generation {i}:")
        print(np.round(genbestsolution,4))
        print(f"aep_best_sol on Genetic Generation {i}:",aep_best_gensol,"GWh")
    return genbestsolution, Gen, NewGen_i_best, aep_best_gensol_i,aep_best_gensol_calc_i,genbestfitness_dic_var_iter_length_x,Gen_Generations_Pack_i,genbest_length,genrankedsolutions_i

def ranksolutions (solutions):
    rankedsolutions = []
    for s in solutions:
        rankedsolutions.append( (fitness(s),s) )
    rankedsolutions.sort(key=lambda a: a[0])   
    bestsolution = rankedsolutions[(Iter_Num-1)]
    return bestsolution,rankedsolutions

def crear_dataframe_con_arrays(arrays, nombres_columnas):
    # Verifica que el número de arrays coincida con el número de nombres de columnas
    if len(arrays) != len(nombres_columnas):
        raise ValueError("El número de arrays debe coincidir con el número de nombres de columnas.")
    # Crea un diccionario para almacenar los datos de cada columna
    datos = {}
    # Asigna los arrays a las columnas correspondientes en el diccionario de datos
    for i, nombre_columna in enumerate(nombres_columnas):
        datos[nombre_columna] = arrays[i]
    # Crea el DataFrame a partir del diccionario de datos
    df = pd.DataFrame(datos)
    return df    

# Getting the fist Solutions set
solutions = []
bestsolution = []
sol = getSolutions()
bestsolution_zero = sol[0]
aep_best_sol_zero = sol[1]
solutions_dic_var_iter_lenght_x = sol[2]
aep_best_sol_iter_length_x =sol[3]
best_length = sol[4]
rankedsolutions_x = sol[5]

# Solutions arrays
array_N = []
array_iter_length = []
n = 0
for i in range (3):
    for j in range (Iter_Num):
        n += 1
        array_N.append(n)

for i in range (3):
    for j in range (Iter_Num):
        array_iter_length.append(i)            

array_x_y = []
array_Fitness = []

for x in rankedsolutions_x.keys():
    for a in rankedsolutions_x[x]:
        array_x_y.append(a[1])
        array_Fitness.append(a[0])     

# Solutions Dataframe creation
df_Solutions = pd.DataFrame()
arrays_sol = [array_iter_length,array_N, array_x_y, array_Fitness]
nombres_columnas_sol = ['Iteration Length','Iterarion','xy', 'Fitness Function_xy_sol']
# Crea el DataFrame
df_Solutions = crear_dataframe_con_arrays(arrays_sol, nombres_columnas_sol)

print(np.round(df_Solutions),4)

# Getting the succesive Genetic Generations
SolGen = getGenSolutions(bestsolution_zero,aep_best_sol_zero)
genbestsolution = SolGen[0]
Gen = SolGen[1]
NewGen_i_best = SolGen[2] 
aep_best_gensol_i = SolGen[3]
aep_best_gensol_calc_i = SolGen[4]
genbestfitness_dic_var_iter_lenght_i = SolGen[5]
Gen_Generations_Pack_i = SolGen[6]
genbest_length = SolGen[7]
genrankedsolutions_i = SolGen[8]

# Succesive Genetic Generation solutions improvement comparison  
Delta_aep_best_gensol_i = []
aep_best_gensol_prev = aep_ref
for s in aep_best_gensol_i:
    Delta_aep_best_gensol_i_elem = s - aep_best_gensol_prev
    aep_best_gensol_prev = s
    Delta_aep_best_gensol_i.append(Delta_aep_best_gensol_i_elem)

# GenSolutions arrays
array_GenN = []
array_gen_iter_length = []
array_i = []
n = 0
for i in range (Gen_Num):
    for j in range (3):
        for k in range (Iter_Num):
            n += 1
            array_GenN.append(n)
for i in range (Gen_Num):
    for j in range (3):
        for k in range (Iter_Num):
            array_gen_iter_length.append(j)            
for i in range (Gen_Num):
    for j in range (3):
        for k in range (Iter_Num):
            array_i.append(i)
            
array_xy_i_x = []
array_Fitness_xy_i_x = []
for i in genrankedsolutions_i:
    for x in i.keys():
        for a in i[x]:
            array_xy_i_x.append(a[1])
            array_Fitness_xy_i_x.append(a[0])
   
# GenSolutions Dataframe creation
df_GenSolutions_i_x = pd.DataFrame()
arrays_gensol = [array_GenN,array_gen_iter_length,array_i, array_xy_i_x, array_Fitness_xy_i_x]
nombres_columnas_gensol = ['Generation','Iteration Length','Iterarion', 'xy', 'Fitness Function_xy_gensol']
# Crea el DataFrame
df_GenSolutions_i_x = crear_dataframe_con_arrays(arrays_gensol, nombres_columnas_gensol)
    
print(np.round(df_GenSolutions_i_x),4)

finalsolution = []
finalsolution = genbestsolution
print("finalsolution")
print(np.round(finalsolution),4)
x_max,y_max = (finalsolution[0],finalsolution[1])
aep_ref = round(float(windFarmModel(wt_x,wt_y).aep().sum()),4)
print("aep_ref:",aep_ref,"GWh")
aep_max = round(float(windFarmModel(x_max,y_max).aep().sum()),4)
print("aep_max:",aep_max,"GWh")

border_x, border_y = zip(*border_points)
plt.figure()
plt.title('Original blue and evolutions of Generation Layouts yellow')
plt.plot(wt_x, wt_y,'b.')
for s in NewGen_i_best:
    plt.plot(s[0], s[1],'y.')     
plt.plot(x_max, y_max, 'g.')
plt.plot(border_x, border_y, 'r-', label='Border')
plt.xlabel('Evolution of positions x [m]')
plt.ylabel('Evolution of positions y [m]')
plt.show()

plt.figure()
plt.title('Original blue and Final black Layouts')
plt.plot(wt_x, wt_y,'b.')
plt.plot(x_max, y_max,'g.')
# wt.plot(x_max, y_max)
plt.plot(border_x, border_y, 'r-', label='Border')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.show()

plt.figure()
plt.title('Evolution of improvement throught genetic evolutions')
plt.plot(Gen, aep_best_gensol_i, 'r.')
plt.xlabel('Genetic Generation i')
plt.ylabel('Best aep Generation i GWh')
plt.show()

plt.figure()
plt.title('Evolution of diferencial improvement throught genetic evolutions')
plt.plot(Gen, Delta_aep_best_gensol_i, 'r.')
plt.xlabel('Genetic Generation i')
plt.ylabel('Best Delta aep Generation i GWh')
plt.show()
    