
# coding: utf-8

# In[1]:

from __future__ import division # ensures no rounding errors from division involving integers
import pandas as pd # gives us the dataframe concept
pd.options.display.max_rows = 9


# In[2]:

aod = pd.read_table('psu_20070731.aod', skiprows=[0,1,2,3,4], sep=" ", skipinitialspace=True)


# In[3]:

aod


# In[4]:

def flag_check(good, aod):
    return aod if good==0 else None


# In[5]:

good_aod501_w_none = [flag_check(good, aod501) for good, aod501 in zip(aod['0=good'],aod['AOD501'])]
good_aod501 = [x for x in good_aod501_w_none if x is not None]
avg_aod501 = sum(good_aod501)/len(good_aod501)
print "Average AOD501 is: " + str(avg_aod501)


# In[6]:

good_aod416_w_none = [flag_check(good, aod416) for good, aod416 in zip(aod['0=good'],aod['AOD416'])]
good_aod416 = [x for x in good_aod416_w_none if x is not None]
avg_aod416 = sum(good_aod416)/len(good_aod416)
print "Average AOD416 is: " + str(avg_aod416)


# In[7]:

good_pressure_w_none = [flag_check(good, press) for good, press in zip(aod['0=good'],aod['p_mb'])]
good_pressure = [x for x in good_pressure_w_none if x is not None]
avg_pressure = sum(good_pressure)/len(good_pressure)
print "Average Pressure (mb) is: " + str(avg_pressure)


# In[ ]:



