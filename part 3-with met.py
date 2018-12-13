
# coding: utf-8

# In[ ]:



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
    maxcons = model.problem.Constraint(oldobj.expression,lb=bmax)
    #maxcons = model.problem.Constraint(oldobj.expression,lb=bmax,ub=bmax)
    model.add_cons_vars(maxcons)
    progre = 0
    for rxn in model.reactions:
        if (progre % 10 == 0):
            print(progre)
        progre += 1
        res[rxn._id] = [0,0]
        for i, direction in enumerate(['min', 'max']):
            model.objective = rxn
            model.objective_direction = direction
            optsoln = model.optimize()
            res[rxn._id][i] = optsoln.objective_value
    model.objective = oldobj
    model.objective_direction = olddir
    model.remove_cons_vars(maxcons)
    
    
    #sol = pd.DataFrame.from_dict(data=res, orient='index', columns=['minimum','maximum'])
    #sol.to_csv('sol.csv',encoding='utf-8')
    
    
    return pd.DataFrame.from_dict(data=res, orient='index', columns=['minimum','maximum'])

# In[32]:



cobra_model = import_matlab_model('C:/Users/farza/Documents/Master_3/Systemes biology/small_ecoli.mat')
fba_solution = cobra_model.optimize()



thermo_data = load_thermoDB('C:/Users/farza/Documents/Master_3/Systemes biology/thermo_data.thermodb')
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


#df_final = pd.DataFrame(index=index, columns=df.columns.values.tolist())
#df_final = pd.DataFrame(index = df.index,columns=df.columns.values.tolist())
df_final = pd.DataFrame(columns=df.columns.values.tolist())
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








# Please make sure that the new generated csv file won't cover the one that you already generated

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




#Aerobic_ FBA
cobra_model.reactions.DM_o2_e.bounds=(-20,0)

cobra_model.reactions.DM_glc_e.bounds=(0,10)
cobra_model.reactions.DM_etoh_e.bounds=(0,10)
cobra_model.reactions.DM_ac_e.bounds=(0,10)
cobra_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)

#
cobra_model.reactions.DM_glc_e.bounds=(-10,10)
fba_solution5 = cobra_model.optimize()

fva_solns5=(fva(cobra_model, fba_solution5))
fva_solns5.to_csv('fva_solns5.csv',encoding='utf-8')

cobra_model.reactions.DM_glc_e.bounds=(0,10)
print(fba_solution5)


#
cobra_model.reactions.DM_etoh_e.bounds=(-10,10)
fba_solution6 = cobra_model.optimize()

fva_solns6=(fva(cobra_model, fba_solution6))
fva_solns6.to_csv('fva_solns6.csv',encoding='utf-8')

cobra_model.reactions.DM_etoh_e.bounds=(0,10)
print(fba_solution6)

#
cobra_model.reactions.DM_ac_e.bounds=(-10,10)
fba_solution7 = cobra_model.optimize()

fva_solns7=(fva(cobra_model, fba_solution7))
fva_solns7.to_csv('fva_solns7.csv',encoding='utf-8')

cobra_model.reactions.DM_ac_e.bounds=(0,10)
print(fba_solution7)

#
cobra_model.reactions.get_by_id('DM_lac-D_e').bounds=(-10,10)
fba_solution8 = cobra_model.optimize()

fva_solns8=(fva(cobra_model, fba_solution8))
fva_solns8.to_csv('fva_solns8.csv',encoding='utf-8')

cobra_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)
print(fba_solution8)


# In[43]:


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


# In[278]:


#Anaerobic_ FBA
cobra_model.reactions.DM_o2_e.bounds=(0,0)

cobra_model.reactions.DM_glc_e.bounds=(0,10)
cobra_model.reactions.DM_etoh_e.bounds=(0,10)
cobra_model.reactions.DM_ac_e.bounds=(0,10)
cobra_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)

#
cobra_model.reactions.DM_glc_e.bounds=(-10,10)
fba_solution13 = cobra_model.optimize()

fva_solns13=(fva(cobra_model, fba_solution13))
fva_solns13.to_csv('fva_solns13.csv',encoding='utf-8')

cobra_model.reactions.DM_glc_e.bounds=(0,10)
print(fba_solution13)


#
cobra_model.reactions.DM_etoh_e.bounds=(-10,10)
fba_solution14 = cobra_model.optimize()
fva_solns14=(fva(cobra_model, fba_solution14))
fva_solns14.to_csv('fva_solns14.csv',encoding='utf-8')
cobra_model.reactions.DM_etoh_e.bounds=(0,10)
print(fba_solution14)

#
cobra_model.reactions.DM_ac_e.bounds=(-10,10)
fba_solution15 = cobra_model.optimize()
fva_solns15=(fva(cobra_model, fba_solution15))
fva_solns15.to_csv('fva_solns15.csv',encoding='utf-8')
cobra_model.reactions.DM_ac_e.bounds=(0,10)
print(fba_solution15)

#
cobra_model.reactions.get_by_id('DM_lac-D_e').bounds=(-10,10)
fba_solution16 = cobra_model.optimize()
fva_solns16=(fva(cobra_model, fba_solution16))
fva_solns16.to_csv('fva_solns16.csv',encoding='utf-8')
cobra_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)
print(fba_solution16)

