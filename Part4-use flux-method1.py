
# coding: utf-8

# In[31]:


import cobra
from math import *
import scipy.io
import numpy
import zlib
import pickle
import pandas as pd
import pytfa
from pytfa.io import import_matlab_model, load_thermoDB


def fva(model, soln):
    res = {}
    bmax = soln.objective_value
    oldobj = model.objective
    olddir = model.objective_direction
    maxcons = model.problem.Constraint(oldobj.expression,lb=bmax,ub=bmax)
    model.add_cons_vars(maxcons)
    progre = 0
    for met in model.metabolites:

        if (progre % 10 == 0):
            print(progre)
        progre += 1
        
        up_b = tfa_model.log_concentration.get_by_id(met._id).variable.ub=1
        lo_b = tfa_model.log_concentration.get_by_id(met._id).variable.lb=-1
        
        res[met._id] = [0,0]
        s_objective = oldobj.expression
        for i, direc in enumerate(['min', 'max']):
            for rxn in model.reactions:
                if met in rxn.metabolites:
                    s= rxn.get_coefficient(met._id)
                else:
                    s=0             
                s_objective = s_objective + s * rxn.flux_expression
      
            model.objective = model.problem.Objective(s_objective - oldobj.expression, direction=direc)

            optsoln = model.optimize()
            
            res[met._id][i] = optsoln.objective_value
        tfa_model.log_concentration.get_by_id(met._id).variable.ub = up_b;
        tfa_model.log_concentration.get_by_id(met._id).variable.lb = lo_b;
            
    model.objective = oldobj
    model.objective_direction = olddir
    model.remove_cons_vars(maxcons)
    return pd.DataFrame.from_dict(data=res, orient='index', columns=['minimum','maximum'])





cobra_model = import_matlab_model('C:/users/vicci/Desktop/vcpre/small_ecoli.mat')
fba_solution = cobra_model.optimize()



thermo_data = load_thermoDB('C:/users/vicci/Desktop/vcpre/thermo_data.thermodb')
tfa_model = pytfa.ThermoModel(thermo_data, cobra_model)
tfa_model.solver = 'optlang-glpk'



tfa_model.prepare()
tfa_model.convert()


## Info on the model
tfa_model.print_info()

## Optimality
tfa_solution = tfa_model.optimize()



# integrating metabolomics data (comment out this part if not integrating)


df = pd.read_csv('metabolomics_data.csv', sep=';')
df['model ID'] = df['model ID'].astype(str)
display(df.head())


df_1 = pd.DataFrame(tfa_model.metabolites)
df_1 = df_1.reset_index()
df_1['ID'] = df_1.iloc[:,1]

list_i = []
for df_i in df_1['ID']:
    #print(str(df_i))
    list_i.append(str(df_i))
print(list_i[11])


df_final = pd.DataFrame(index=index, columns=df.columns.values.tolist())
df = df.fillna(0)
for i, ID in enumerate(list_i):
    #print(ID)
    df_i = df[df['model ID']==ID]
    if (not df_i.empty):
        df_final = pd.concat([df_final, df_i], axis=0)
        element = df_i.iloc[:,0].tolist()[0]
        up_b = float(df_i.iloc[:,2].tolist()[0]) + float(abs(df_i.iloc[:,3].tolist()[0]))
        
        lo_b = float(df_i.iloc[:,2].tolist()[0]) - float(abs(df_i.iloc[:,3].tolist()[0]))
        if lo_b<=0:
            lo_b = 0.0000001
        tfa_model.log_concentration.get_by_id(element).variable.ub = log(up_b)
        tfa_model.log_concentration.get_by_id(element).variable.lb = log(lo_b)

# This is the end of the metabolics data integration





# !!!!!!!!!!Please make sure that the new generated csv file won't cover the one that you already generated

#fva_solns = []


#Aerobic_ TFA
tfa_model.reactions.DM_o2_e.bounds=(-20,0)

tfa_model.reactions.DM_glc_e.bounds=(0,10)
tfa_model.reactions.DM_etoh_e.bounds=(0,10)
tfa_model.reactions.DM_ac_e.bounds=(0,10)
tfa_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)

#
tfa_model.reactions.DM_glc_e.bounds=(-10,10)
tfa_solution1 = tfa_model.optimize()

fva_solns1=(fva(tfa_model, tfa_solution1))
fva_solns1.to_csv('fva_solns1.csv',encoding='utf-8')

tfa_model.reactions.DM_glc_e.bounds=(0,10)
print(tfa_solution1)


#
tfa_model.reactions.DM_etoh_e.bounds=(-10,10)
tfa_solution2 = tfa_model.optimize()

fva_solns2=(fva(tfa_model, tfa_solution2))
fva_solns2.to_csv('fva_solns2.csv',encoding='utf-8')

tfa_model.reactions.DM_etoh_e.bounds=(0,10)
print(tfa_solution2)

#
tfa_model.reactions.DM_ac_e.bounds=(-10,10)
tfa_solution3 = tfa_model.optimize()

fva_solns3=(fva(tfa_model, tfa_solution3))
fva_solns3.to_csv('fva_solns3.csv',encoding='utf-8')

tfa_model.reactions.DM_ac_e.bounds=(0,10)
print(tfa_solution3)

#
tfa_model.reactions.get_by_id('DM_lac-D_e').bounds=(-10,10)
tfa_solution4 = tfa_model.optimize()

fva_solns4=(fva(tfa_model, tfa_solution4))
fva_solns4.to_csv('fva_solns4.csv',encoding='utf-8')

tfa_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)
print(tfa_solution4)





#Anaerobic_ TFA
tfa_model.reactions.DM_o2_e.bounds=(0,0)

tfa_model.reactions.DM_glc_e.bounds=(0,10)
tfa_model.reactions.DM_etoh_e.bounds=(0,10)
tfa_model.reactions.DM_ac_e.bounds=(0,10)
tfa_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)

#
tfa_model.reactions.DM_glc_e.bounds=(-10,10)
tfa_solution9 = tfa_model.optimize()
fva_solns9=(fva(tfa_model, tfa_solution9))
fva_solns9.to_csv('fva_solns9.csv',encoding='utf-8')
tfa_model.reactions.DM_glc_e.bounds=(0,10)
print(tfa_solution9)


#
tfa_model.reactions.DM_etoh_e.bounds=(-10,10)
tfa_solution10 = tfa_model.optimize()
fva_solns10=(fva(tfa_model, tfa_solution10))
fva_solns10.to_csv('fva_solns10.csv',encoding='utf-8')
tfa_model.reactions.DM_etoh_e.bounds=(0,10)
print(tfa_solution10)

#
tfa_model.reactions.DM_ac_e.bounds=(-10,10)
tfa_solution11 = tfa_model.optimize()
fva_solns11=(fva(tfa_model, tfa_solution11))
fva_solns11.to_csv('fva_solns11.csv',encoding='utf-8')
tfa_model.reactions.DM_ac_e.bounds=(0,10)
print(tfa_solution11)

#
tfa_model.reactions.get_by_id('DM_lac-D_e').bounds=(-10,10)
tfa_solution12 = tfa_model.optimize()
fva_solns12=(fva(tfa_model, tfa_solution12))
fva_solns12.to_csv('fva_solns12.csv',encoding='utf-8')
tfa_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)
print(tfa_solution12)


