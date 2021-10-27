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


df = pd.concat([df1, df2], ignore_index=True)

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



#%%Make Regression Fit

from scipy.optimize import curve_fit
def func(x, a,b): #Define function for Plotting
    return a*x+b
def lin_reg(i,plot): #Defining function that gives fit parameters for manhole i
    times = [t_sec(time)+7225 for time in df["TIME"]] #Get time corresponding to manhole measurement
    i_start = i_time((t_sec(t_end[i])-50),times)
    i_end = i_time(t_sec(t_end[i]),times)
    xdata = df.CO2_dry[i_start:i_end]-CO2_bg #Calculate measured CO2 excess concentrations
    ydata = df.CH4_dry[i_start:i_end]-CH4_bg #Calculate measured CH4 excess concentrations
    popt, pcov = curve_fit(func, xdata, ydata) #make fit
    if plot: #optional: Scatter plot concentrations and fitted function
        plt.scatter(xdata,ydata)
        plt.plot(xdata,func(xdata,popt[0],popt[1]))
        plt.xlabel('CO2 [ppm]')
        plt.ylabel('CH4 [ppm]')
    return popt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
fig, axs = plt.subplots(4,5,figsize=(18,14))
for i in range(len(t_end)):
    times = [t_sec(time)+7225 for time in df["TIME"]] #Get time corresponding to manhole measurement
    i_start = i_time((t_sec(t_end[i])-50),times)
    i_end = i_time(t_sec(t_end[i]),times)
    xdata = df.CO2_dry[i_start:i_end]-403 #Calculate measured CO2 excess concentrations
    ydata = df.CH4_dry[i_start:i_end]-1.97 #Calculate measured CH4 excess concentrations
    popt, pcov = curve_fit(func, xdata, ydata) #make fit
    k=int(i/5)
    l=i%5
    popt, pcov = curve_fit(func, xdata, ydata)
    axs[k, l].scatter(xdata,ydata,s=2)
    axs[k, l].plot(xdata,func(xdata,popt[0],popt[1]))
    axs[k, l].set_title('Manhole {0}'.format(i+1))
    axs[k,l].tick_params(axis='both', which='major', labelsize=9)
    axs[k,l].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    axs[k,0].set_ylabel('c(CH4) [ppm]')
    axs[3,l].set_xlabel('c(CO2) [ppm]')
plt.savefig('fit_subplots.pdf')
    
#print(lin_reg(5,plot=True))

params=[[],[]]
for i in range(len(t_end)):
    popt=lin_reg(i, plot=False)
    params[0].append(popt[0])
    params[1].append(popt[1])
print(params[0])#print regression slope values for all manholes

#plt.scatter(range(len(params[0])),params[0])
#plt.yscale('log')
#%%Combined Slope plot
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
data = pd.read_csv(os.path.join(path,"data.csv"))
print(data.keys())
xmax=1000 #ppm
x=np.linspace(0,xmax,100)
colorcode={'gas':'green','sewage':'red','rain':'blue','undefined':'black'}
for i in range(len(data['manholes'])):
    plt.plot(x,data['ratio'][i]*x+data['intercept'][i],linewidth=1,color=colorcode[data['pipe type'][i]])
l1,l2,l3,l4,p1 = Line2D([0], [0], label='rain', color='blue'),Line2D([0], [0], label='sewage', color='red'),Line2D([0], [0], label='gas', color='green'),Line2D([0], [0], label='undefined', color='black'),mpatches.Patch(color='lightsteelblue', label='combustion regime') 
plt.ylim([-1000,8000])
plt.fill_between(x,20*x,color='lightsteelblue')
plt.ylabel('CH4 [ppb]')
plt.xlabel('CO2 [ppm]')
plt.legend(handles=[l1,l2,l3,l4,p1])
plt.savefig('combined_slopes.pdf')
#%%Plotting data

from mpl_toolkits.axes_grid1 import make_axes_locatable

#Load map and setting coordinates
lat_max =52.0881
lat_min =52.0847   

lon_max = 5.1775
lon_min = 5.1632

BBox = ((lon_min,lon_max,lat_min,lat_max))
im = plt.imread(os.path.join(path,"map_52.0881_5.1775_52.0847_5.1632.png"))

#Dictionary for coding different pipe types on map
pipetype = {'sewege':'v',
            'rain':'o',
            'gas':'s'}

fig,ax = plt.subplots(figsize=(12,8))

vmin = np.min(params[0])
vmax = np.max(params[0])
for i in range(len(params[0])):
    
    plot = ax.scatter(coor['lon'][i],coor['lat'][i],
                      c=params[0][i],
                      zorder=1,s=40, 
                      cmap='plasma',vmin=vmin,vmax=vmax,alpha=0.8,
                      marker=pipetype[coor["pipetype"][i]])
    
    
ax.set_xlim(lon_min,lon_max)
ax.set_ylim(lat_min,lat_max)

ax.imshow(im,zorder=0,extent=BBox,aspect='equal',alpha=0.5)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.05)
cbar = plt.colorbar(plot, cax=cax)
cbar.set_label("ppm/ppm")


#%%Save data to csv file

#Build dataframe with data
data = pd.DataFrame(data={'manholes':manholes,
                          'end times':t_end,
                          'ratio':params[0],
                          'intercept':params[1],
                          'lat':coor["lat"],
                          'lon':coor["lon"],
                          'pipe type':coor["pipetype"]})

print(data)

data.to_csv(os.path.join(path,'data.csv'),index=False)



#%%Make plots of all times of interest (ToI), times at which the data is selected
import matplotlib
matplotlib.use('Agg')
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)

#Times in seconds, local time, corresponding to data
times = [t_sec(time)+7225 for time in df["TIME"]] 

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
    axs[0].xaxis.set_minor_locator(MultipleLocator(2))
    axs[1].xaxis.set_minor_locator(MultipleLocator(2))
    
    axs[0].tick_params(labelrotation=45)
    axs[1].tick_params(labelrotation=45)
    
    axs[0].set_ylabel('CH4 [ppm]')
    axs[1].set_ylabel('CO2 [ppm]')
    axs[1].set_xlabel('Time (s)')
    fig.suptitle('Manhole %i'%(i+1))
    
    
    #Plotting shaded area of ToI
    area_CH4 = axs[0].axvspan(t_min,t_max,color='darksalmon',ec='red')
    area_CO2 = axs[1].axvspan(t_min,t_max,color='darksalmon',ec='red')

    #Plot data little before and after ToI
    overshoot = 10 #Amount of datapoints overshoot before and after time of interest
    times_plot = times[int(i_start-overshoot):int(i_end+overshoot)] #Timespan that is plotted
    CH4 = df['CH4_dry'][int(i_start-overshoot):int(i_end+overshoot)]-CH4_bg #Select data
    CO2 = df['CO2_dry'][int(i_start-overshoot):int(i_end+overshoot)]-CO2_bg #Select data
    
    CH4_plot = axs[0].plot(times_plot,CH4)
    C02_plot = axs[1].plot(times_plot,CO2)
    
    fig.savefig(os.path.join(path,'figures_ToI',"manhole_%i.jpg"%(i+1)))
    
    
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

axs[0].set_ylabel('CH4 [ppm]')
axs[1].set_ylabel('CO2 [ppm]')
axs[1].set_xlabel('Time (s)')
fig.suptitle('Manhole %i'%(i+1))

#Plotting shaded area of ToI
area_CH4 = axs[0].axvspan(t_min,t_max,color='darksalmon',ec='red')
area_CO2 = axs[1].axvspan(t_min,t_max,color='darksalmon',ec='red')

#Plot data little before and after ToI
undershoot = 500 #Amount of datapoints overshoot before and after time of interest
overshoot = 10 #Amount of datapoints overshoot before and after time of interest

times_plot = times[int(i_start-undershoot):int(i_end+overshoot)] #Timespan that is plotted
CH4 = df['CH4'][int(i_start-undershoot):int(i_end+overshoot)]-CH4_bg #Select data
CO2 = df['CO2'][int(i_start-undershoot):int(i_end+overshoot)]-CO2_bg #Select data

CH4_plot = axs[0].plot(times_plot,CH4)
C02_plot = axs[1].plot(times_plot,CO2)

fig.savefig(os.path.join(path,"figures_individual_measurements","manhole_%i.jpg"%(i+1)))
print(MultipleLocator(10))

#%%


