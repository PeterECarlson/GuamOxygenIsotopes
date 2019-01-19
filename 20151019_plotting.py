# -*- coding: utf-8 -*-
"""
Created on Monday, October 19, 2015 

@author: anoronha
"""
import os
import datetime
from pylab import *
import pandas as pd
import numpy as np
from data_processing import plot_all

path = os.getcwd()
today = datetime.datetime.now().strftime("%Y%m%d")

FieldCO2Continous = pd.read_pickle('FieldCO2Continous.pickle')
FieldCO2Continous = FieldCO2Continous.sort_values('Datetime')
FieldCO2Continous = FieldCO2Continous.set_index(pd.DatetimeIndex(FieldCO2Continous['Datetime']))
FieldCO2_Daily = FieldCO2Continous.resample('D').mean()
FieldCO2_Daily = FieldCO2_Daily.reset_index(drop = False)
FieldCO2_Daily.rename(columns={'index': 'Date'}, inplace=True)
FieldCO2_Daily = FieldCO2_Daily[np.isfinite(FieldCO2_Daily['CO2'])]

fig_name = 'tmp'

savefig('%s/%s_%s.png' %(path, today, fig_name), format='png', dpi=1000)

Fig1 = plt.figure(num=1, figsize = (11,8.5), 
                  facecolor='w', edgecolor='k')
ax1 = Fig1.add_subplot(111)
ax1.plot(FieldCO2_Daily.Datetime,FieldCO2_Daily.CO2,'o', 
     markersize = 6, markeredgewidth=0.5, markeredgecolor='k',
     markerfacecolor = 'w' ,linestyle = '-', color = 'k')
ax1.ylabel(r'pCO2')
ax1.xlabel(r'Date')
ax1.title('Stumpy')

#JinapsanTrips = FieldTrip[FieldTrip.Location == 'Jinapsan Cave']
#for i in range(0, len(JinapsanTrips)):
#    ax1.axvline(JinapsanTrips.BeginFieldTrip.iloc[i], color = 'r')


plot_all(JinapsanCaveResults, 'Stumpy','Midpoint','1000lnalpha','w',1, 0.3)


d18OWaterResults = pd.read_csv('reprocessed_d18OWaterResults.csv')
