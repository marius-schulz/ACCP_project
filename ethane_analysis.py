# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 15:04:09 2021

@author: vanlo
"""

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
def data_retr(df,t_end,):
    #df needs to have column names TIMES with the measurement times, notation hh:mm:ss.ss
    times = [t_sec(time)+7225 for time in df["TIME"]]  #In seconds from start of day in local dutch time
    
    #Find index corresponding to t_end
    i_start = i_time((t_sec(t_end)-60),times)
    i_end = i_time(t_sec(t_end),times)
    
    #Compute ratio as timeseries for manhole measurement
    ratio = df.CH4[i_start:i_end]/df.CO2[i_start:i_end]
    ratio_dry = df.CH4_dry[i_start:i_end]/df.CO2_dry[i_start:i_end]
    
    ratio_avg = np.average(ratio.to_numpy())
    ratio_dry_avg = np.average(ratio_dry.to_numpy())

    return ratio_avg, ratio_dry_avg    

    

#%%Import data

#Script runs if directories with data are places in same folder as this python file
path = os.path.dirname(__file__)

df1 = pd.read_csv(os.path.join(path, "G4302\\15\\NOMAD-20211015-125852Z-DataLog_User_Minimal.dat"),
                  sep="\s+")
df2 = pd.read_csv(os.path.join(path, "G4302\\15\\NOMAD-20211015-135856Z-DataLog_User_Minimal.dat"),
                  sep="\s+")
df3 = pd.read_csv(os.path.join(path, "G4302\\15\\NOMAD-20211015-145900Z-DataLog_User_Minimal.dat"),
                  sep="\s+")


df = pd.concat([df1,df2,df3], ignore_index=True)

# The manholes numbering/naming
manholes = np.arange(1,20)

#Import time stamps
timestamps = pd.read_csv(os.path.join(path,"measurement_times.csv"))

# End times of each measurement to retrieve the data coresponding to the manholes
# Use notation 13:10:05 for 5 seconds past 10 minutes past 13h. 
t_end = timestamps["END"]  

#Import coordinates
coor = pd.read_csv(os.path.join(path,"coordinates.csv"),delimiter=',')

#Background levels CH4 and CO2 [ppm] for the dry measurements
CH4_bg = 1.99
CO2_bg = 413

u_CH4_bg = 0.002
u_CO2_bg = 8



#%%Make plots of all times of interest (ToI), times at which the data is selected
import matplotlib
matplotlib.use('Agg')
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)

overshoot = 20 #Amount of datapoints overshoot before and after time of interest

#Times in seconds, local time, corresponding to data
times = [t_sec(time)+7200 for time in df["TIME"]] 

#Make plot of analysed data for every ToI
for i in range(len(t_end)):
    #Select ToI
        #In terms of index in all data
    i_start = i_time(t_sec(t_end[i])-50,times)
    i_end = i_time(t_sec(t_end[i]),times)
        #In terms of time
    t_min = times[i_time(t_sec(t_end[i])-50,times)]
    t_max = times[(i_time(t_sec(t_end[i]),times))]
    
    #Plotting aesthetics
    fig,axs = plt.subplots(2,1,sharex=True)
    
    axs[0].xaxis.set_major_locator(MultipleLocator(10))
    axs[1].xaxis.set_major_locator(MultipleLocator(10))
    # axs[0].xaxis.set_minor_locator(MultipleLocator(2))
    # axs[1].xaxis.set_minor_locator(MultipleLocator(2))
    axs[0].set_xticklabels([])
    axs[1].set_xticklabels([])
    
    axs[0].tick_params(labelrotation=45)
    axs[1].tick_params(labelrotation=45)
    
    axs[0].set_ylabel('CH4 [ppm]')
    axs[1].set_ylabel('C2H6 [ppm]')
    axs[1].set_xlabel('Time [10 s]')
    fig.suptitle('Manhole %i'%(i+1))
    
    #Plotting shaded area of ToI
    area_CH4 = axs[0].axvspan(t_min,t_max,color='darksalmon',ec='red')
    area_CO2 = axs[1].axvspan(t_min,t_max,color='darksalmon',ec='red')

    #Plot data little before and after ToI
    
    times_plot = times[int(i_start-overshoot):int(i_end+overshoot)] #Timespan that is plotted
    CH4 = df['CH4_dry'][int(i_start-overshoot):int(i_end+overshoot)]-CH4_bg #Select data
    CO2 = df['C2H6_dry'][int(i_start-overshoot):int(i_end+overshoot)] #Select data
    
    CH4_plot = axs[0].plot(times_plot,CH4)
    C02_plot = axs[1].plot(times_plot,CO2)
    
    fig.savefig(os.path.join(path,'figures_ToI_ethane',"manhole_%i.jpg"%(i+1)))
    
    
#%%Make plot of ToI of specific manhole for further examination
import matplotlib
matplotlib.use('Agg')
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)


#Select manhole:
manhole = 1
i = manhole-1

#Select ToI
    #In terms of index in all data
i_start = i_time(t_sec(t_end[i])-50,times)
i_end = i_time(t_sec(t_end[i]),times)
    #In terms of time
t_min = times[i_time(t_sec(t_end[i])-50,times)]
t_max = times[(i_time(t_sec(t_end[i]),times))]

#Plotting aesthetics
fig,axs = plt.subplots(2,1)

axs[0].xaxis.set_major_locator(MultipleLocator(60))
axs[1].xaxis.set_major_locator(MultipleLocator(60))
axs[0].xaxis.set_minor_locator(MultipleLocator(10))
axs[1].xaxis.set_minor_locator(MultipleLocator(10))

axs[0].tick_params(labelrotation=45)
axs[1].tick_params(labelrotation=45)
axs[0].set_xticklabels([])
axs[1].set_xticklabels([])
    
axs[0].set_ylabel('CH4 [ppm]')
axs[1].set_ylabel('CO2 [ppm]')
axs[1].set_xlabel('Time [10 s]')
fig.suptitle('Manhole %i'%(i+1))

#Plotting shaded area of ToI
area_CH4 = axs[0].axvspan(t_min,t_max,color='darksalmon',ec='red')
area_CO2 = axs[1].axvspan(t_min,t_max,color='darksalmon',ec='red')

#Plot data little before and after ToI
undershoot = 20 #Amount of datapoints overshoot before and after time of interest
overshoot = 20 #Amount of datapoints overshoot before and after time of interest

times_plot = times[int(i_start-undershoot):int(i_end+overshoot)] #Timespan that is plotted
CH4 = df['CH4'][int(i_start-undershoot):int(i_end+overshoot)]-CH4_bg #Select data
CO2 = df['CO2'][int(i_start-undershoot):int(i_end+overshoot)]-CO2_bg #Select data

CH4_plot = axs[0].plot(times_plot,CH4)
C02_plot = axs[1].plot(times_plot,CO2)

fig.savefig(os.path.join(path,"figures_individual_measurements","manhole_%i.jpg"%(i+1)))
print(MultipleLocator(10))



