# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 15:52:15 2021

@author: elian
"""

import numpy as np
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt

'''
Functions
'''
#Convert time string to seconds 
#Time should have notation hh:mm:ss ,decimals in seconds allowed
def t_sec(time):
    h = int(time[:2])
    m = int(time[3:5])
    s = float(time[6:])
    
    t_s = 3600*h + 60*m + s
    return t_s

#Find index of data corresponding to t in array named times [seconds]
def i_time(t,times):
    times = [int(time) for time in times]
    print(times)
    try: 
        i = times.index(int(t))
    except:
        # Try for a second before, might have missed a second in measurements.
        try: 
            i = times.index(int(t)-1)
        except:   
            sys.exit("Time indexing failed")
            
    return i

#Using end time to retrieve the corresponding measurements in the manholes
#subsequently take average of CH4:CO2 ratio
def data_retr(df,t_end):
    #df needs to have column names TIMES with the measurement times, notation hh:mm:ss.ss
    times = [t_sec(time) for time in df["TIME"]]  #In seconds from start of day 
    
    #Find index corresponding to t_end
    t = t_sec(t_end)
    i_start = i_time((t_sec(t_end)-60),times)
    i_end = i_time(t_sec(t_end),times)
    
    print(i_end)
    
    #Compute ratio as timeseries for manhole measurement
    ratio = df.CH4[i_start:i_end+1]/df.CO2[i_start:i_end+1]
    ratio_dry = df.CH4_dry[i_start:i_end+1]/df.CO2_dry[i_start:i_end+1]
    
    print(ratio)
    


#Script runs if directories with data are places in same folder as this python file
path = os.path.dirname(__file__)

#Open relevant G2301 data as dataframe
# df1 = pd.read_table(os.path.join(path, "G2301\\15\\CFADS2394-20211015-130036Z-DataLog_User.dat"), 
#                     sep="\s+", 
#                     usecols=['DATE', 'TIME','FRAC_DAYS_SINCE_JAN1','FRAC_HRS_SINCE_JAN1','JULIAN_DAYS','EPOCH_TIME','ALARM_STATUS','INST_STATUS','CavityPressure','CavityTemp','DasTemp','EtalonTemp','WarmBoxTemp','species','MPVPosition','OutletValve','solenoid_valves','CO2','CO2_dry','CH4','CH4_dry','H2O'])
# df2 = pd.read_table(os.path.join(path,"G2301\\15\\CFADS2394-20211015-140042Z-DataLog_User.dat"), sep="\s+", usecols=['DATE', 'TIME','FRAC_DAYS_SINCE_JAN1','FRAC_HRS_SINCE_JAN1','JULIAN_DAYS','EPOCH_TIME','ALARM_STATUS','INST_STATUS','CavityPressure','CavityTemp','DasTemp','EtalonTemp','WarmBoxTemp','species','MPVPosition','OutletValve','solenoid_valves','CO2','CO2_dry','CH4','CH4_dry','H2O'])
# df=np.concatenate((df1,df2))

df1 = pd.read_csv(os.path.join(path, "G2301\\15\\CFADS2394-20211015-130036Z-DataLog_User.dat"),
                  sep="\s+")
df2 = pd.read_csv(os.path.join(path, "G2301\\15\\CFADS2394-20211015-140042Z-DataLog_User.dat"),
                  sep="\s+")
df3 = pd.read_csv(os.path.join(path, "G2301\\15\\CFADS2394-20211015-150048Z-DataLog_User.dat"),
                  sep="\s+")
df4 = pd.read_csv(os.path.join(path, "G2301\\15\\CFADS2394-20211015-155955Z-DataLog_User.dat"),
                  sep="\s+")

df = pd.concat([df1, df2, df3, df4], ignore_index=True)
print(df["TIME"])

# The manholes numbering/naming
manholes = np.arange(1,21)

# End times of each measurement to retrieve the data coresponding to the manholes
# Use notation 13:10:05 for 5 seconds past 10 minutes past 13h. 
t_end = '15:30:00'

plt.plot(df["TIME"])

#data_retr(df,t_end)



#data = pd.DataFrame(data={'manholes':manholes,'end times':t_end},index=manholes)



