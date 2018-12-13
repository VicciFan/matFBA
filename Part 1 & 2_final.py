
# coding: utf-8

# In[1]:


import cobra
from math import *
import scipy.io
import numpy
import zlib
import pickle
import pytfa
from pytfa.io import import_matlab_model, load_thermoDB


# In[2]:



cobra_model = import_matlab_model('C:/Users/farza/Documents/Master_3/Systemes biology/small_ecoli.mat')


# In[3]:



# Run FBA on the model 
fba_solution = cobra_model.optimize()


# In[4]:


thermo_data = load_thermoDB('C:/Users/farza/Documents/Master_3/Systemes biology/thermo_data.thermodb')
tfa_model = pytfa.ThermoModel(thermo_data, cobra_model)
tfa_model.solver = 'optlang-glpk'


# In[5]:


## TFA conversion
tfa_model.prepare()
tfa_model.convert()

## Info on the model
tfa_model.print_info()

## Optimality
tfa_solution = tfa_model.optimize()


# In[6]:


#Aerobic_ TFA
tfa_model.reactions.DM_o2_e.bounds=(-20,0)

tfa_model.reactions.DM_glc_e.bounds=(0,10)
tfa_model.reactions.DM_etoh_e.bounds=(0,10)
tfa_model.reactions.DM_ac_e.bounds=(0,10)
tfa_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)

#
tfa_model.reactions.DM_glc_e.bounds=(-10,10)
tfa_solution1 = tfa_model.optimize()
tfa_model.reactions.DM_glc_e.bounds=(0,10)
print(tfa_solution1)


#
tfa_model.reactions.DM_etoh_e.bounds=(-10,10)
tfa_solution2 = tfa_model.optimize()
tfa_model.reactions.DM_etoh_e.bounds=(0,10)
print(tfa_solution2)

#
tfa_model.reactions.DM_ac_e.bounds=(-10,10)
tfa_solution3 = tfa_model.optimize()
tfa_model.reactions.DM_ac_e.bounds=(0,10)
print(tfa_solution3)

#
tfa_model.reactions.get_by_id('DM_lac-D_e').bounds=(-10,10)
tfa_solution4 = tfa_model.optimize()
tfa_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)
print(tfa_solution4)


# In[7]:


#Aerobic_ FBA
cobra_model.reactions.DM_o2_e.bounds=(-20,0)

cobra_model.reactions.DM_glc_e.bounds=(0,10)
cobra_model.reactions.DM_etoh_e.bounds=(0,10)
cobra_model.reactions.DM_ac_e.bounds=(0,10)
cobra_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)

#
cobra_model.reactions.DM_glc_e.bounds=(-10,10)
fba_solution = cobra_model.optimize()
cobra_model.reactions.DM_glc_e.bounds=(0,10)
print(fba_solution)


#
cobra_model.reactions.DM_etoh_e.bounds=(-10,10)
fba_solution = cobra_model.optimize()
cobra_model.reactions.DM_etoh_e.bounds=(0,10)
print(fba_solution)

#
cobra_model.reactions.DM_ac_e.bounds=(-10,10)
fba_solution = cobra_model.optimize()
cobra_model.reactions.DM_ac_e.bounds=(0,10)
print(fba_solution)

#
cobra_model.reactions.get_by_id('DM_lac-D_e').bounds=(-10,10)
fba_solution = cobra_model.optimize()
cobra_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)
print(fba_solution)


# In[8]:


#Anaerobic_ TFA
tfa_model.reactions.DM_o2_e.bounds=(0,0)

tfa_model.reactions.DM_glc_e.bounds=(0,10)
tfa_model.reactions.DM_etoh_e.bounds=(0,10)
tfa_model.reactions.DM_ac_e.bounds=(0,10)
tfa_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)

#
tfa_model.reactions.DM_glc_e.bounds=(-10,10)
tfa_solution1 = tfa_model.optimize()
tfa_model.reactions.DM_glc_e.bounds=(0,10)
print(tfa_solution1)


#
tfa_model.reactions.DM_etoh_e.bounds=(-10,10)
tfa_solution2 = tfa_model.optimize()
tfa_model.reactions.DM_etoh_e.bounds=(0,10)
print(tfa_solution2)

#
tfa_model.reactions.DM_ac_e.bounds=(-10,10)
tfa_solution3 = tfa_model.optimize()
tfa_model.reactions.DM_ac_e.bounds=(0,10)
print(tfa_solution3)

#
tfa_model.reactions.get_by_id('DM_lac-D_e').bounds=(-10,10)
tfa_solution4 = tfa_model.optimize()
tfa_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)
print(tfa_solution4)


# In[9]:


#Anaerobic_ FBA
cobra_model.reactions.DM_o2_e.bounds=(0,0)

cobra_model.reactions.DM_glc_e.bounds=(0,10)
cobra_model.reactions.DM_etoh_e.bounds=(0,10)
cobra_model.reactions.DM_ac_e.bounds=(0,10)
cobra_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)

#
cobra_model.reactions.DM_glc_e.bounds=(-10,10)
fba_solution = cobra_model.optimize()
cobra_model.reactions.DM_glc_e.bounds=(0,10)
print(fba_solution)


#
cobra_model.reactions.DM_etoh_e.bounds=(-10,10)
fba_solution = cobra_model.optimize()
cobra_model.reactions.DM_etoh_e.bounds=(0,10)
print(fba_solution)

#
cobra_model.reactions.DM_ac_e.bounds=(-10,10)
fba_solution = cobra_model.optimize()
cobra_model.reactions.DM_ac_e.bounds=(0,10)
print(fba_solution)

#
cobra_model.reactions.get_by_id('DM_lac-D_e').bounds=(-10,10)
fba_solution = cobra_model.optimize()
cobra_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)
print(fba_solution)


# In[6]:


print(type(tfa_model.metabolites))


# In[7]:


import pandas as pd
df = pd.read_csv('metabolomics_data.csv', sep=';')
df['model ID'] = df['model ID'].astype(str)
display(df.head())


# In[8]:


df_1 = pd.DataFrame(tfa_model.metabolites)
df_1 = df_1.reset_index()
df_1['ID'] = df_1.iloc[:,1]
#display(df_1)


# In[9]:


list_i = []
for df_i in df_1['ID']:
    #print(str(df_i))
    list_i.append(str(df_i))
print(list_i[11])


# In[10]:


df.loc[df['model ID'].isin(list_i)]


# In[11]:


# Integrate some conncentration data
df_final = pd.DataFrame(columns=df.columns.values.tolist())
df = df.fillna(0)


# In[17]:


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
        
display(df_final.reset_index(drop=True))

tfa_solution_w_concentration = tfa_model.optimize()


# In[110]:


#Aerobic_ TFA
tfa_model.reactions.DM_o2_e.bounds=(-20,0)

tfa_model.reactions.DM_glc_e.bounds=(0,10)
tfa_model.reactions.DM_etoh_e.bounds=(0,10)
tfa_model.reactions.DM_ac_e.bounds=(0,10)
tfa_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)

#
tfa_model.reactions.DM_glc_e.bounds=(-10,10)
tfa_solution1 = tfa_model.optimize()
tfa_model.reactions.DM_glc_e.bounds=(0,10)
print(tfa_solution1)


#
tfa_model.reactions.DM_etoh_e.bounds=(-10,10)
tfa_solution2 = tfa_model.optimize()
tfa_model.reactions.DM_etoh_e.bounds=(0,10)
print(tfa_solution2)

#
tfa_model.reactions.DM_ac_e.bounds=(-10,10)
tfa_solution3 = tfa_model.optimize()
tfa_model.reactions.DM_ac_e.bounds=(0,10)
print(tfa_solution3)

#
tfa_model.reactions.get_by_id('DM_lac-D_e').bounds=(-10,10)
tfa_solution4 = tfa_model.optimize()
tfa_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)
print(tfa_solution4)


# In[111]:


#Aerobic_ FBA
cobra_model.reactions.DM_o2_e.bounds=(-20,0)

cobra_model.reactions.DM_glc_e.bounds=(0,10)
cobra_model.reactions.DM_etoh_e.bounds=(0,10)
cobra_model.reactions.DM_ac_e.bounds=(0,10)
cobra_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)

#
cobra_model.reactions.DM_glc_e.bounds=(-10,10)
fba_solution = cobra_model.optimize()
cobra_model.reactions.DM_glc_e.bounds=(0,10)
print(fba_solution)


#
cobra_model.reactions.DM_etoh_e.bounds=(-10,10)
fba_solution = cobra_model.optimize()
cobra_model.reactions.DM_etoh_e.bounds=(0,10)
print(fba_solution)

#
cobra_model.reactions.DM_ac_e.bounds=(-10,10)
fba_solution = cobra_model.optimize()
cobra_model.reactions.DM_ac_e.bounds=(0,10)
print(fba_solution)

#
cobra_model.reactions.get_by_id('DM_lac-D_e').bounds=(-10,10)
fba_solution = cobra_model.optimize()
cobra_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)
print(fba_solution)


# In[112]:


#Anaerobic_ TFA
tfa_model.reactions.DM_o2_e.bounds=(0,0)

tfa_model.reactions.DM_glc_e.bounds=(0,10)
tfa_model.reactions.DM_etoh_e.bounds=(0,10)
tfa_model.reactions.DM_ac_e.bounds=(0,10)
tfa_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)

#
tfa_model.reactions.DM_glc_e.bounds=(-10,10)
tfa_solution1 = tfa_model.optimize()
tfa_model.reactions.DM_glc_e.bounds=(0,10)
print(tfa_solution1)


#
tfa_model.reactions.DM_etoh_e.bounds=(-10,10)
tfa_solution2 = tfa_model.optimize()
tfa_model.reactions.DM_etoh_e.bounds=(0,10)
print(tfa_solution2)

#
tfa_model.reactions.DM_ac_e.bounds=(-10,10)
tfa_solution3 = tfa_model.optimize()
tfa_model.reactions.DM_ac_e.bounds=(0,10)
print(tfa_solution3)

#
tfa_model.reactions.get_by_id('DM_lac-D_e').bounds=(-10,10)
tfa_solution4 = tfa_model.optimize()
tfa_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)
print(tfa_solution4)


# In[33]:


#Anaerobic_ FBA
cobra_model.reactions.DM_o2_e.bounds=(0,0)

cobra_model.reactions.DM_glc_e.bounds=(0,10)
cobra_model.reactions.DM_etoh_e.bounds=(0,10)
cobra_model.reactions.DM_ac_e.bounds=(0,10)
cobra_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)

#
cobra_model.reactions.DM_glc_e.bounds=(-10,10)
fba_solution = cobra_model.optimize()
cobra_model.reactions.DM_glc_e.bounds=(0,10)
print(fba_solution)


#
cobra_model.reactions.DM_etoh_e.bounds=(-10,10)
fba_solution = cobra_model.optimize()
cobra_model.reactions.DM_etoh_e.bounds=(0,10)
print(fba_solution)

#
cobra_model.reactions.DM_ac_e.bounds=(-10,10)
fba_solution = cobra_model.optimize()
cobra_model.reactions.DM_ac_e.bounds=(0,10)
print(fba_solution)

#
cobra_model.reactions.get_by_id('DM_lac-D_e').bounds=(-10,10)
fba_solution = cobra_model.optimize()
cobra_model.reactions.get_by_id('DM_lac-D_e').bounds=(0,10)
print(fba_solution)

