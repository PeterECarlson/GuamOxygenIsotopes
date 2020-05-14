# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 13:10:44 2016

@author: anoronha
"""

import pandas as pd
from pandas import DataFrame
import numpy as np
import os
from pylab import *
import datetime
import time
from scipy import stats, signal

molar_mass_Ca = 40.08
molar_mass_calcite = 100.09
sec_per_min = 60.0
min_per_hour = 60.0
hour_per_day = 24.0
seconds_per_day = sec_per_min * min_per_hour * hour_per_day
start_plot = datetime.datetime(2008, 6, 1)
end_plot = datetime.datetime(2016, 8, 1)

DripwaterData = pd.read_pickle('DripwaterData.pickle')
RainwaterData = pd.read_pickle('RainwaterData.pickle')
DripRate = pd.read_pickle('DripRate.pickle')
CationResults = pd.read_pickle('CationResults.pickle')
FieldCO2Continous = pd.read_pickle('FieldCO2Continous.pickle')
CalciteData = pd.read_pickle('CalciteData.pickle')
DripChem = pd.read_pickle('DripChem.pickle')
Climate = pd.read_pickle('Climate.pickle')
GroundwaterSampleDetails = pd.read_pickle('GroundwaterSampleDetails.pickle')
#d18OWaterResults_unfiltered = pd.read_pickle('d18OWaterResults_unfiltered.pickle')
FieldTrip = pd.read_pickle('FieldTrip.pickle')
PlateMass = pd.read_pickle('PlateMass_all.pickle')
Alkalinity = pd.read_pickle('Alkalinity.pickle')
JinapsanCaveResults = pd.read_pickle('JinapsanCaveResults.pickle')
dripwater_measured = pd.read_pickle('dripwater_measured.pickle')
ModeledvsMeasGrowth = pd.read_pickle('ModeledvsMeasGrowth.pickle')



DripNames = ['Flatman','Station 2', 'Station 1','Stumpy']

def rgb_convert(palette):
    for i in range(len(palette)):
        r, g, b = palette[i]
        palette[i] = (r / 255., g / 255., b / 255.)
    return palette

plot_table = DataFrame(columns = ['DripName','Color','Marker'])
plot_table.DripName = DripNames
plot_table.Color = rgb_convert([(27,158,119), (217,95,2), (117,112,179), (231,41,138)])
plot_table.Marker = ['o','D','v','^']
plot_table = plot_table.sort_values('DripName')


guam_data = pd.concat([JinapsanCaveResults.SampleName,
                       JinapsanCaveResults.SiteName,
                       JinapsanCaveResults.DripInterval,
                       JinapsanCaveResults.Temp,
                       JinapsanCaveResults.Ca,
                       JinapsanCaveResults.CO2,
                       JinapsanCaveResults.d18OWater], axis=1)
guam_data = guam_data[guam_data.SiteName.isin(DripNames)]
guam_data['DripInterval'] = [int(x) for x in guam_data.DripInterval]
guam_data = guam_data[np.isnan(guam_data.Ca) == False]
guam_data = guam_data[np.isnan(guam_data.d18OWater) == False]
guam_data = guam_data.reset_index(drop=True)
isolution_key = guam_data.SampleName
isolution_key  = isolution_key.reset_index(drop=False)
guam_data = guam_data.drop(['SampleName','SiteName'],1)
guam_data.to_csv('guam_data.csv',header=False)

isolution_output_affek = pd.read_csv('isolution_output_affek.csv',names=['index','d18OCalcite_Modeled_Affek'])
isolution_output_kim = pd.read_csv('isolution_output_kim.csv',names=['index','d18OCalcite_Modeled_Kim'])
isolution_output_coplen = pd.read_csv('isolution_output_coplen.csv',names=['index','d18OCalcite_Modeled_Coplen'])
isolution_output_affek = pd.merge(isolution_output_affek, isolution_key, how='left',on='index')
isolution_output_kim = pd.merge(isolution_output_kim, isolution_key, how='left',on='index')
isolution_output_coplen = pd.merge(isolution_output_coplen, isolution_key, how='left',on='index')
JinapsanCaveResults = pd.merge(JinapsanCaveResults, isolution_output_kim, how='left',on='SampleName')
JinapsanCaveResults = pd.merge(JinapsanCaveResults, isolution_output_affek, how='left',on='SampleName')
JinapsanCaveResults = pd.merge(JinapsanCaveResults, isolution_output_coplen, how='left',on='SampleName')
holder = []
for i in range(0,len(JinapsanCaveResults)):
#    equilib = (exp(17.4/(tmp.Temp.iloc[i]+273.15)-0.0286))*(1000+tmp.d18OWater.iloc[i])-1000
    equilib_affek = (exp(15.63/(JinapsanCaveResults.Temp.iloc[i]+273.15)-0.02329))*(1000+JinapsanCaveResults.d18OWater.iloc[i])-1000
    equilib_affek = (equilib_affek - 30.91)/1.03091
    holder.append(equilib_affek)
JinapsanCaveResults['equilib_affek'] = holder
holder = []
for i in range(0,len(JinapsanCaveResults)):
#    equilib = (exp(17.4/(tmp.Temp.iloc[i]+273.15)-0.0286))*(1000+tmp.d18OWater.iloc[i])-1000
    equilib_kim = (exp(18.03/(JinapsanCaveResults.Temp.iloc[i]+273.15)-0.03242))*(1000+JinapsanCaveResults.d18OWater.iloc[i])-1000
    equilib_kim = (equilib_kim - 30.91)/1.03091
    holder.append(equilib_kim)
JinapsanCaveResults['equilib_kim'] = holder
holder = []
for i in range(0,len(JinapsanCaveResults)):
#    equilib = (exp(17.4/(tmp.Temp.iloc[i]+273.15)-0.0286))*(1000+tmp.d18OWater.iloc[i])-1000
    equilib_coplen = (exp(17.4/(JinapsanCaveResults.Temp.iloc[i]+273.15)-0.0286))*(1000+JinapsanCaveResults.d18OWater.iloc[i])-1000
    equilib_coplen = (equilib_coplen - 30.91)/1.03091
    holder.append(equilib_coplen)
JinapsanCaveResults['equilib_coplen'] = holder
#JinapsanCaveResults['mismatch_equilib'] = tmp.equilib - tmp.d18O_VPDB
JinapsanCaveResults['mismatch'] = JinapsanCaveResults.d18OCalcite_Modeled_Affek - JinapsanCaveResults.d18O_VPDB

#%%Fig 1
plot1_min = min(JinapsanCaveResults['1000lnalpha'])
plot1_max = max(JinapsanCaveResults['1000lnalpha'])

print(min(JinapsanCaveResults['1000lnalpha']))
print(max(JinapsanCaveResults['1000lnalpha']))

fig, axes = plt.subplots(nrows=2,ncols=2,sharex=True, figsize=(12, 6))
for i, ax in enumerate(axes.flatten()):
#        for j in range(0,len(grid_lines)):
#            ax.axvline(grid_lines[j], color='0.75', linestyle='--', lw=1)
    plotname = plot_table.DripName.iloc[i]
    plotmarker = plot_table.Marker.iloc[i]
    plotcolor = plot_table.Color.iloc[i]
    set_to_plot = JinapsanCaveResults[JinapsanCaveResults.SiteName == plotname]
    set_to_plot = set_to_plot.sort_values('Midpoint')
    ax.plot(set_to_plot['Midpoint'], set_to_plot['1000lnalpha'],
                 linestyle = '-',
                 color = plotcolor,
                 marker = plotmarker,
                 markeredgewidth = 1,
                 linewidth = 1,
                 markeredgecolor = plotcolor,
                 label = plotname,
                 markersize = 5)
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    for label in ax.yaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    ax.set_title(plotname)
    plot_ylim = ax.get_ylim()
    if plot_ylim[0] >=0:
        plot_ylim = (plot_ylim[0] - .1*plot_ylim[1], plot_ylim[1])
    ax.set_ylim([plot1_min, plot1_max])
    ax.set_xlim([start_plot, end_plot])
#ax0.set_xlim([start_plot, end_plot])
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    if i == 1 or i == 3:
        for label in ax.yaxis.get_ticklabels():
            label.set_visible(False)
    else:
        ax.set_ylabel(r'$1000\ln\alpha$', fontsize=15)          
plt.savefig('1000lnalpha_agu.png', format='png', dpi=600,bbox_inches='tight')

#%%Fig 2 

DripwaterData = DripwaterData[DripwaterData.SiteName.isin(DripNames)]
DripwaterData = DripwaterData[DripwaterData.d18OWater < -4]
DripwaterData = pd.merge(DripwaterData, FieldTrip, how='left',on='idFieldTrip')

plot2_min = -10#min(DripwaterData['d18OWater'])
plot2_max = -4#max(DripwaterData['d18OWater'])

plot3_min = -10#min(CalciteData['d18O_VPDB'])
plot3_max = -4#max(CalciteData['d18O_VPDB'])


fig, axes = plt.subplots(nrows=2,ncols=4,sharex=True, figsize=(10, 5.63))
for i, ax in enumerate(axes.flatten()):
    if i < 4:
        plotname = plot_table.DripName.iloc[i]
        plotmarker = plot_table.Marker.iloc[i]
        plotcolor = plot_table.Color.iloc[i]
        set_to_plot = DripwaterData[DripwaterData.SiteName == plotname]
        set_to_plot = set_to_plot.sort_values('EndFieldTrip')
        ax.plot(set_to_plot['EndFieldTrip'], set_to_plot['d18OWater'],
                     linestyle = '-',
                     color = plotcolor,
                     marker = plotmarker,
                     markeredgewidth = 1,
                     linewidth = 1,
                     markeredgecolor = plotcolor,
                     label = plotname,
                     markersize = 5)
        for label in ax.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
    #    for label in ax.yaxis.get_ticklabels()[::2]:
    #        label.set_visible(False)
        ax.set_title(plotname)
        plot_ylim = ax.get_ylim()
    #    if plot_ylim[0] >=0:
    #        plot_ylim = (plot_ylim[0] - .1*plot_ylim[1], plot_ylim[1])
        ax.set_ylim([plot2_min, plot2_max])
        ax.set_xlim([start_plot, end_plot])
    #ax0.set_xlim([start_plot, end_plot])
        for label in ax.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        if i > 0 and i < 4:
            for label in ax.yaxis.get_ticklabels():
                label.set_visible(False)   
        if i == 0:
            ax.set_ylabel(u'$\delta^{18}$O Water \n ($\u2030$ VSMOW)', fontsize=15)
    else:
        i = i - 4
        plotname = plot_table.DripName.iloc[i]
        plotmarker = plot_table.Marker.iloc[i]
        plotcolor = plot_table.Color.iloc[i]
        set_to_plot = CalciteData[CalciteData.SiteName == plotname]
        set_to_plot = set_to_plot.sort_values('Midpoint')
        ax.plot(set_to_plot['Midpoint'], set_to_plot['d18O_VPDB'],
                     linestyle = '-',
                     color = plotcolor,
                     marker = plotmarker,
                     markeredgewidth = 1,
                     linewidth = 1,
                     markeredgecolor = plotcolor,
                     label = plotname,
                     markersize = 5)
        for label in ax.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
    #    for label in ax.yaxis.get_ticklabels()[::2]:
    #        label.set_visible(False)
        plot_ylim = ax.get_ylim()
    #    if plot_ylim[0] >=0:
    #        plot_ylim = (plot_ylim[0] - .1*plot_ylim[1], plot_ylim[1])
        ax.set_ylim([plot3_min, plot3_max])
        ax.set_xlim([start_plot, end_plot])
    #ax0.set_xlim([start_plot, end_plot])
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)
        for label in ax.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        if i > 0 and i < 4:
            for label in ax.yaxis.get_ticklabels():
                label.set_visible(False)    
        if i == 0:
            ax.set_ylabel(u'$\delta^{18}O$ Calcite$ \n ($\u2030$ VPDB)', fontsize=15)
            plt.subplots_adjust(wspace=0.1, hspace=0.1,bottom=0.1, top=0.85)
plt.savefig('measured_d18o_agu.png', format='png', dpi=600,bbox_inches='tight')

#%%Fig 3 1000lnalpha vs drip rate and deposition rate.
fig4 = plt.figure(figsize=(10,5.6))
ax2 = fig4.add_subplot(121)      
to_regress = JinapsanCaveResults[(JinapsanCaveResults['1000lnalpha'].isnull()==False) & (JinapsanCaveResults.DripsPerMin.isnull() == False)]
text_y = max(to_regress['1000lnalpha'])
text_x = max(to_regress['DripsPerMin'])*0.5
slope, intercept, r_value, p_value, std_err = stats.linregress(to_regress.DripsPerMin, to_regress['1000lnalpha'])
yp = polyval([slope,intercept], to_regress.DripsPerMin)
for i in range(0,len(plot_table)):
    plotname = plot_table.DripName.iloc[i]
    plotmarker = plot_table.Marker.iloc[i]
    plotcolor = plot_table.Color.iloc[i]
    set_to_plot = to_regress[to_regress.SiteName == plotname]
    ax2.plot(set_to_plot['DripsPerMin'], set_to_plot['1000lnalpha'],
                 linestyle = 'none',
                 color = plotcolor,
                 marker = plotmarker,
                 markeredgewidth = 1,
                 linewidth = 1,
                 markeredgecolor = plotcolor,
                 label = plotname,
                 markersize = 5)
#for label in ax2.xaxis.get_ticklabels()[::2]:
#    label.set_visible(False)
ax2.tick_params(axis='both', which='major', labelsize=15)
plot(to_regress.DripsPerMin, yp, 'k-')
xlabel(u'Drip rate  \n(drip/min)', fontsize=15, labelpad=15)
ylabel(r'$1000\ln\alpha$', fontsize=15)
text(text_x,text_y-0.05,'$r^{2} =$ %5.2f' %(r_value**2), fontsize=15, ha='center')
text(text_x,text_y-0.2,'$p =$ %.1E' %(p_value), fontsize=15, ha='center')
ax2.set_xlim([0,7])
tmp_ticks_labels = [10, 20, 30, 60, 300]
tmp_ticks_pos = [60.0/x for x in tmp_ticks_labels]
ax3 = ax2.twiny()
ax3.set_xlim([0,7])
ax3.set_xticks(tmp_ticks_pos)
ax3.set_xticklabels(tmp_ticks_labels)
ax3.tick_params(axis='x', labelsize=12)
ax3.set_xlabel(u'Drip interval\n (seconds)', fontsize=15)
#plt.savefig('dripate_alpha_agu.png', format='png', dpi=600, bbox_inches='tight')
ax4 = fig4.add_subplot(122)      
to_regress = JinapsanCaveResults[(JinapsanCaveResults['1000lnalpha'].isnull()==False) & (JinapsanCaveResults.GrowthRate.isnull() == False)]
to_regress = to_regress[to_regress.GrowthRate > 0]
to_regress['GrowthRate'] = to_regress.GrowthRate * 1000
text_y = max(to_regress['1000lnalpha'])
text_x = max(to_regress['GrowthRate'])*0.5
slope, intercept, r_value, p_value, std_err = stats.linregress(to_regress.GrowthRate, to_regress['1000lnalpha'])
yp = polyval([slope,intercept], to_regress.GrowthRate)
for i in range(0,len(plot_table)):
    plotname = plot_table.DripName.iloc[i]
    plotmarker = plot_table.Marker.iloc[i]
    plotcolor = plot_table.Color.iloc[i]
    set_to_plot = to_regress[to_regress.SiteName == plotname]
    ax4.plot(set_to_plot['GrowthRate'], set_to_plot['1000lnalpha'],
                 linestyle = 'none',
                 color = plotcolor,
                 marker = plotmarker,
                 markeredgewidth = 1,
                 linewidth = 1,
                 markeredgecolor = plotcolor,
                 label = plotname,
                 markersize = 5)
for label in ax4.yaxis.get_ticklabels():
    label.set_visible(False)
ax4.tick_params(axis='both', which='major', labelsize=15)
plot(to_regress.GrowthRate, yp, 'k-')
xlabel(u'Calcite Deposition Rate  \n(mg/day)', fontsize=15, labelpad=15)
#ylabel(r'$1000\ln\alpha$', fontsize=15)
text(text_x,text_y-0.05,'$r^{2} =$ %5.2f' %(r_value**2), fontsize=15, ha='center')
text(text_x,text_y-0.2,'$p =$ %.1E' %(p_value), fontsize=15, ha='center')
plt.subplots_adjust(wspace=0.1, hspace=0.1,bottom=0.1, top=0.85)
plt.savefig('dripgrowth_alpha_agu.png', format='png', dpi=600, bbox_inches='tight')

#%% Fig 4 drip d18O vs Drip rate
fig4 = plt.figure(figsize=(7,5.6))
ax2 = fig4.add_subplot(111)      
to_regress = JinapsanCaveResults[(JinapsanCaveResults['d18OWater'].isnull()==False) & (JinapsanCaveResults.DripsPerMin.isnull() == False)]
to_regress = to_regress[to_regress.DripsPerMin > 0]
text_y = max(to_regress['d18OWater'])
text_x = max(to_regress['DripsPerMin'])*0.5
slope, intercept, r_value, p_value, std_err = stats.linregress(to_regress.DripsPerMin, to_regress['d18OWater'])
yp = polyval([slope,intercept], to_regress.DripsPerMin)
for i in range(0,len(plot_table)):
    plotname = plot_table.DripName.iloc[i]
    plotmarker = plot_table.Marker.iloc[i]
    plotcolor = plot_table.Color.iloc[i]
    set_to_plot = to_regress[to_regress.SiteName == plotname]
    ax2.plot(set_to_plot['DripsPerMin'], set_to_plot['d18OWater'],
                 linestyle = 'none',
                 color = plotcolor,
                 marker = plotmarker,
                 markeredgewidth = 1,
                 linewidth = 1,
                 markeredgecolor = plotcolor,
                 label = plotname,
                 markersize = 5)
#for label in ax2.xaxis.get_ticklabels()[::2]:
#    label.set_visible(False)
ax3.tick_params(axis='both', which='major', labelsize=15)
for label in ax3.yaxis.get_ticklabels():
    label.set_visible(False)
plot(to_regress.DripsPerMin, yp, 'k-')
xlabel(u'Drip Rate  \n(mg/day)', fontsize=15, labelpad=15)
ylabel(r'$Drip Water \delta^{18}O (VSMOW)$', fontsize=15)
text(text_x,text_y-0.05,'$r^{2} =$ %5.2f' %(r_value**2), fontsize=15, ha='center')
text(text_x,text_y-0.2,'$p =$ %.1E' %(p_value), fontsize=15, ha='center')
plt.savefig('Driprated18OCorrelation.png', format='png', dpi=600, bbox_inches='tight')

#%%Fig 5 Drip rate timeseries
fig5 = plt.figure(figsize=(10,5.63))
ax = fig5.add_subplot(111)      
for i in range(0,len(plot_table)):
    plotname = plot_table.DripName.iloc[i]
    plotmarker = plot_table.Marker.iloc[i]
    plotcolor = plot_table.Color.iloc[i]
    set_to_plot = DripRate[DripRate.SiteName == plotname]
    set_to_plot = set_to_plot.sort_values('EndFieldTrip')
    ax.plot(set_to_plot['EndFieldTrip'], set_to_plot['DripsPerMin'],
                 linestyle = '-',
                 color = plotcolor,
                 marker = plotmarker,
                 markeredgewidth = 1,
                 linewidth = 1,
                 markeredgecolor = plotcolor,
                 label = plotname,
                 markersize = 5)
    ax.set_xlim([start_plot, end_plot])
    ax.tick_params(axis='both', which='major', labelsize=15)
for label in ax.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)
ax.set_ylim([0,7])
legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), frameon=False, fancybox=False, shadow=False, ncol=4,numpoints = 1,prop={'size':15})
tmp_ticks_labels = [10, 20, 30, 60, 120, 300]
tmp_ticks_pos = [60.0/x for x in tmp_ticks_labels]
ax2 = ax.twinx()
ax2.set_ylim([0,7])
ax2.set_yticks(tmp_ticks_pos)
ax2.set_yticklabels(tmp_ticks_labels)
ax2.tick_params(axis='y', labelsize=12)
ax2.set_ylabel(u'Drip interval\n (seconds)', fontsize=15)
ax.set_ylabel(u'Drip rate \n(drips/min)', fontsize=15)
#text(-7,-5,'$r^{2} = %5.2f$' %(r_value), fontsize=15)
#plt.subplots_adjust(wspace=0, hspace=0.1,bottom=0.1, top=0.85)
plt.savefig('driprate_agu.png', format='png', dpi=600)

#%%Fig 6
fig6_tick = [-10,-9,-8,-7, -6, -5]
plot_lim = [-10,-5]
fig6 = plt.figure(figsize=(10,5.5))
ax = fig6.add_subplot(231)      
for i in range(0,len(plot_table)):
    plotname = plot_table.DripName.iloc[i]
    plotmarker = plot_table.Marker.iloc[i]
    plotcolor = plot_table.Color.iloc[i]
    set_to_plot = JinapsanCaveResults[JinapsanCaveResults.SiteName == plotname]
#    set_to_plot = set_to_plot.sort_values('EndFieldTrip')
    ax.plot(set_to_plot['d18O_VPDB'], set_to_plot['equilib_kim'],
                 linestyle = 'none',
                 color = plotcolor,
                 marker = plotmarker,
                 markeredgewidth = 1,
                 linewidth = 1,
                 markeredgecolor = plotcolor,
                 label = plotname,
                 markersize = 5)
    ax.tick_params(axis='both', which='major', labelsize=15)
ax.plot([-10,-5],[-10,-5],'k-')
ax.set_xlim(plot_lim)
ax.set_ylim(plot_lim)
ax.set_xticks(fig6_tick)
ax.set_yticks(fig6_tick)
ax.set_yticklabels(fig6_tick)
#for label in ax.yaxis.get_ticklabels()[::2]:
#    label.set_visible(False)
for label in ax.xaxis.get_ticklabels():
    label.set_visible(False)
ax.set_title('Kim and O\'Neil (1997)', fontsize=15)
#legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), frameon=False, fancybox=False, shadow=False, ncol=4,numpoints = 1,prop={'size':15})
ax.set_ylabel(u'$\delta^{18}$O Modeled \nEqulibirum', fontsize=15)
ax1 = fig6.add_subplot(232)      
for i in range(0,len(plot_table)):
    plotname = plot_table.DripName.iloc[i]
    plotmarker = plot_table.Marker.iloc[i]
    plotcolor = plot_table.Color.iloc[i]
    set_to_plot = JinapsanCaveResults[JinapsanCaveResults.SiteName == plotname]
    ax1.plot(set_to_plot['d18O_VPDB'], set_to_plot['equilib_affek'],
                 linestyle = 'none',
                 color = plotcolor,
                 marker = plotmarker,
                 markeredgewidth = 1,
                 linewidth = 1,
                 markeredgecolor = plotcolor,
                 label = plotname,
                 markersize = 5)
    ax1.tick_params(axis='both', which='major', labelsize=15)
ax1.set_title('Affek and Zaarur (2014)', fontsize=15)
ax1.plot([-10,-5],[-10,-5],'k-')
ax1.set_xticks(fig6_tick)
ax1.set_yticks(fig6_tick)
ax1.set_xlim(plot_lim)
ax1.set_ylim(plot_lim)
for label in ax1.xaxis.get_ticklabels():
    label.set_visible(False)
for label in ax1.yaxis.get_ticklabels():
    label.set_visible(False)
ax2 = fig6.add_subplot(233)      
for i in range(0,len(plot_table)):
    plotname = plot_table.DripName.iloc[i]
    plotmarker = plot_table.Marker.iloc[i]
    plotcolor = plot_table.Color.iloc[i]
    set_to_plot = JinapsanCaveResults[JinapsanCaveResults.SiteName == plotname]
    ax2.plot(set_to_plot['d18O_VPDB'], set_to_plot['equilib_coplen'],
                 linestyle = 'none',
                 color = plotcolor,
                 marker = plotmarker,
                 markeredgewidth = 1,
                 linewidth = 1,
                 markeredgecolor = plotcolor,
                 label = plotname,
                 markersize = 5)
    ax2.tick_params(axis='both', which='major', labelsize=15)
ax2.set_title('Coplen (2002)', fontsize=15)
ax2.plot([-10,-5],[-10,-5],'k-')
ax2.set_xticks(fig6_tick)
ax2.set_yticks(fig6_tick)
ax2.set_xlim(plot_lim)
ax2.set_ylim(plot_lim)
for label in ax2.xaxis.get_ticklabels():
    label.set_visible(False)
for label in ax2.yaxis.get_ticklabels():
    label.set_visible(False)    

ax3 = fig6.add_subplot(234)      
for i in range(0,len(plot_table)):
    plotname = plot_table.DripName.iloc[i]
    plotmarker = plot_table.Marker.iloc[i]
    plotcolor = plot_table.Color.iloc[i]
    set_to_plot = JinapsanCaveResults[JinapsanCaveResults.SiteName == plotname]
    ax3.plot(set_to_plot['d18O_VPDB'], set_to_plot['d18OCalcite_Modeled_Kim'],
                 linestyle = 'none',
                 color = plotcolor,
                 marker = plotmarker,
                 markeredgewidth = 1,
                 linewidth = 1,
                 markeredgecolor = plotcolor,
                 label = plotname,
                 markersize = 5)
    ax3.tick_params(axis='both', which='major', labelsize=15)
ax3.plot([-10,-5],[-10,-5],'k-')
ax3.set_xlim(plot_lim)
ax3.set_ylim(plot_lim)
ax3.set_xticks(fig6_tick)
ax3.set_xticklabels(fig6_tick)
ax3.set_yticks(fig6_tick)
ax3.set_yticklabels(fig6_tick)
#for label in ax1.yaxis.get_ticklabels()[::2]:
#    label.set_visible(False)
#for label in ax1.xaxis.get_ticklabels()[::2]:
#    label.set_visible(False)
ax3.set_ylabel(u'$\delta^{18}$O Modeled \nRayleigh', fontsize=15)

ax4 = fig6.add_subplot(235)      
for i in range(0,len(plot_table)):
    plotname = plot_table.DripName.iloc[i]
    plotmarker = plot_table.Marker.iloc[i]
    plotcolor = plot_table.Color.iloc[i]
    set_to_plot = JinapsanCaveResults[JinapsanCaveResults.SiteName == plotname]
    ax4.plot(set_to_plot['d18O_VPDB'], set_to_plot['d18OCalcite_Modeled_Affek'],
                 linestyle = 'none',
                 color = plotcolor,
                 marker = plotmarker,
                 markeredgewidth = 1,
                 linewidth = 1,
                 markeredgecolor = plotcolor,
                 \
                 label = plotname,
                 markersize = 5)
    ax4.tick_params(axis='both', which='major', labelsize=15)
ax4.plot([-10,-5],[-10,-5],'k-')
ax4.set_xticks(fig6_tick)
ax4.set_yticks(fig6_tick)
ax4.set_xlim(plot_lim)
ax4.set_ylim(plot_lim)
#ax3.xticks(fig6_tick, fig6_tick)
ax4.set_xticks(fig6_tick)
ax4.set_xticklabels(fig6_tick)
for label in ax4.yaxis.get_ticklabels():
    label.set_visible(False)

ax4 = fig6.add_subplot(236)      
for i in range(0,len(plot_table)):
    plotname = plot_table.DripName.iloc[i]
    plotmarker = plot_table.Marker.iloc[i]
    plotcolor = plot_table.Color.iloc[i]
    set_to_plot = JinapsanCaveResults[JinapsanCaveResults.SiteName == plotname]
    ax4.plot(set_to_plot['d18O_VPDB'], set_to_plot['d18OCalcite_Modeled_Coplen'],
                 linestyle = 'none',
                 color = plotcolor,
                 marker = plotmarker,
                 markeredgewidth = 1,
                 linewidth = 1,
                 markeredgecolor = plotcolor,
                 label = plotname,
                 markersize = 5)
    ax4.tick_params(axis='both', which='major', labelsize=15)
ax4.plot([-10,-5],[-10,-5],'k-')
ax4.set_xlim(plot_lim)
ax4.set_ylim(plot_lim)
#ax3.xticks(fig6_tick, fig6_tick)
ax4.set_xticks(fig6_tick)
ax4.set_yticks(fig6_tick)
ax4.set_xticklabels(fig6_tick)
for label in ax4.yaxis.get_ticklabels():
    label.set_visible(False)    
fig6.text(0.5, -0.02, u'$\delta^{18}$O Calcite Measured', fontsize=15, ha='center')
plt.subplots_adjust(wspace=0.2, hspace=0.2,bottom=0.1, top=0.85)
plt.savefig('model_results_agu.png', format='png', dpi=600, bbox_inches='tight')

#%%Fig 7

#fig7 = plt.figure(figsize = (10,5.63), facecolor='w', edgecolor='k')
#ax1 = fig7.add_subplot(111)
#ax1.set_ylim([25,38])
#ax1.set_xlim([3.32,3.35])
#for i in range(0,len(plot_table)):
#    plotname = plot_table.DripName.iloc[i]
#    plotmarker = plot_table.Marker.iloc[i]
#    plotcolor = plot_table.Color.iloc[i]
#    set_to_plot = JinapsanCaveResults[JinapsanCaveResults.SiteName == plotname]
#    ax1.plot(1000/(set_to_plot['Temp'] + 273.15), set_to_plot['1000lnalpha'],
#             linestyle = 'none',
#             color = plotcolor,
#             marker = plotmarker,
#             markeredgewidth = 1,
#             linewidth = 1,
#             markeredgecolor = plotcolor,
#             label = plotname,
#             markersize = 5)
#Temp_range = [-5,50]
#Temp_range = [1000/(i+273.15) for i in Temp_range]
#Coplen_alpha = [17.4*(i)-28.6 for i in Temp_range]
#Affek_alpha = [15.63*(i)-23.29 for i in Temp_range]
#KimONeil_alpha = [18.03*(i)-32.42 for i in Temp_range]
#ax1.plot(Temp_range, 
#         KimONeil_alpha, 
#         linestyle=':', 
#         linewidth=3,
#         color = '#e41a1c',
#         label = 'Kim and O\'Neil (1997)')
#ax1.plot(Temp_range, 
#         Coplen_alpha, 
#         linestyle='-', 
#         linewidth=3,
#         color='#377eb8',
#         label = 'Coplen (2002)')
#ax1.plot(Temp_range, 
#         Affek_alpha, 
#         linestyle='--',
#         linewidth=3,
#         color='#4daf4a',
#         label = 'Affek and Zaarur (2014)')
#ax1.set_xlim([3.32,3.35])
#ax1.set_ylim([27,36])
#ax1.tick_params(axis='x', labelsize=20)
#ax1.tick_params(axis='y', labelsize=20)
#for label in ax1.yaxis.get_ticklabels()[::2]:
#    label.set_visible(False)
#xlabel(r'1000/T (K)', fontsize=20)
##    box = ax1.get_position()
##    ax1.set_position([box.x0, box.y0, box.width * 0.5, box.height])
#ax1.legend(loc='upper left',numpoints = 1, frameon=False)
#ax2 = ax1.twiny()
##    ax2.set_position([box.x0, box.y0, box.width * 0.5, box.height])
#tmp_ticks_labels = [27.0, 26.5, 26.0]
#tmp_ticks_pos = [1000.0/(x+273.15) for x in tmp_ticks_labels]
#ax2.set_xlim([3.32,3.35])
##    tmp_pos = np.array([3.2, 3.4, 3.5, 3.6])
#ax2.set_xticks(tmp_ticks_pos)
#ax2.set_xticklabels(tmp_ticks_labels)
#ax2.tick_params(axis='x', labelsize=20)
#ax2.set_xlabel(u'T (\u00b0C)', fontsize=20)
#ax1.set_ylabel(r'$1000\ln\alpha$', fontsize=20)
#savefig('AGU2016_Fig2.png', format='png', dpi=1000, bbox_inches='tight')


#fig6 = plt.figure(figsize=(10,5.63))
#ax = fig6.add_subplot(111)      
#for i in range(0,len(plot_table)):
#    plotname = plot_table.DripName.iloc[i]
#    plotmarker = plot_table.Marker.iloc[i]
#    plotcolor = plot_table.Color.iloc[i]
#    set_to_plot = dripwater_measured[dripwater_measured.SiteName == plotname]
#    set_to_plot = set_to_plot.sort_values('EndFieldTrip')
#    ax.plot(set_to_plot['EndFieldTrip'], set_to_plot['Ca']/(1e3 * molar_mass_Ca),
#                 linestyle = '-',
#                 color = plotcolor,
#                 marker = plotmarker,
#                 markeredgewidth = 1,
#                 linewidth = 1,
#                 markeredgecolor = plotcolor,
#                 label = plotname,
#                 markersize = 5)
#    ax.set_xlim([start_plot, end_plot])
#    ax.tick_params(axis='both', which='major', labelsize=15)
#for label in ax.xaxis.get_ticklabels()[::2]:
#    label.set_visible(False)
#legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), frameon=False, fancybox=False, shadow=False, ncol=4,numpoints = 1,prop={'size':15})
#ylabel(u'[Ca$^{2+}$] (mmol/L)', fontsize=15)
#plt.savefig('calcium_agu.png', format='png', dpi=600)


#to_regress = to_regress[to_regress.GrowthRate > 0]
#to_regress['SA'] = [(x*0.27)/2.0 for x in to_regress.Growth]
#to_regress['Time_sec'] = to_regress.DaysDeploy * seconds_per_day
#to_regress['Mols'] = [(x/molar_mass_calcite) for x in to_regress.Growth]
#to_regress['calc_R'] = (to_regress.Mols/to_regress.SA)/to_regress.Time_sec
#to_regress['log_R'] = [log10(x) for x in to_regress.calc_R]
#text_x = min(to_regress['log_R'])
#to_regress['log_R'] = [log10(x) for x in to_regress.]

JinapsanCaveResults['d18OWaterPDB'] = (JinapsanCaveResults.d18OWater - 30.91)/1.03091 
JinapsanCaveResults['Modeled_Appar_Frac'] = 1000*(log(((JinapsanCaveResults.d18OCalcite_Modeled_Affek/1000)+1)/((JinapsanCaveResults.d18OWaterPDB/1000)+1)))
JinapsanCaveResults['Expect_Frac'] = 1000*(log(((JinapsanCaveResults.equilib_affek/1000)+1)/((JinapsanCaveResults.d18OWaterPDB/1000)+1)))

plot7_min = 28#min(pd.concat(JinapsanCaveResults['Modeled_Appar_Frac'], JinapsanCaveResults['Expect_Frac'], JinapsanCaveResults['1000lnalpha']))
plot7_max = 31#min(pd.concat(JinapsanCaveResults['Modeled_Appar_Frac'], JinapsanCaveResults['Expect_Frac'], JinapsanCaveResults['1000lnalpha']))

fig, axes = plt.subplots(nrows=2,ncols=2,sharex=True, figsize=(12, 6))
for i, ax in enumerate(axes.flatten()):
#        for j in range(0,len(grid_lines)):
#            ax.axvline(grid_lines[j], color='0.75', linestyle='--', lw=1)
    plotname = plot_table.DripName.iloc[i]
    plotmarker = plot_table.Marker.iloc[i]
    plotcolor = plot_table.Color.iloc[i]
    set_to_plot = JinapsanCaveResults[JinapsanCaveResults.SiteName == plotname]
    set_to_plot = set_to_plot.sort_values('Midpoint')

    ax.plot(set_to_plot['Midpoint'], set_to_plot['1000lnalpha'],
                 linestyle = '-',
                 color = 'tab:gray',
                 marker = plotmarker,
                 markeredgewidth = 1,
                 linewidth = 1,
                 markeredgecolor = 'tab:gray',
                 markerfacecolor='None',
                 label = plotname,
                 markersize = 5)
    ax.plot(set_to_plot['Midpoint'], set_to_plot['Modeled_Appar_Frac'],
                 linestyle = '-',
                 color = plotcolor,
                 marker = plotmarker,
                 markeredgewidth = 1,
                 linewidth = 1,
                 markeredgecolor = plotcolor,
                 label = plotname,
                 markersize = 5)
    ax.plot(set_to_plot['Midpoint'], set_to_plot['Expect_Frac'],
                 linestyle = '-',
                 color = 'k',
                 marker = '',
                 markeredgewidth = 1,
                 linewidth = 1,
                 markeredgecolor = plotcolor,
                 label = plotname,
                 markersize = 5)
    print (min(set_to_plot['Expect_Frac']))
    print (max(set_to_plot['Expect_Frac']))
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    for label in ax.yaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    ax.set_title(plotname)
    plot_ylim = ax.get_ylim()
    if plot_ylim[0] >=0:
        plot_ylim = (plot_ylim[0] - .1*plot_ylim[1], plot_ylim[1])
    ax.set_ylim([plot7_min, plot7_max])
    ax.set_xlim([start_plot, end_plot])
#ax0.set_xlim([start_plot, end_plot])
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    if i == 1 or i ==3:
        for label in ax.yaxis.get_ticklabels():
            label.set_visible(False)
    else:
        ax.set_ylabel(r'$1000\ln\alpha$', fontsize=15)          
plt.savefig('apparent_modeled_1000lnalpha.png', format='png', dpi=600,bbox_inches='tight')

DripwaterData = DripwaterData[DripwaterData.SiteName.isin(DripNames)]
DripwaterData = DripwaterData[DripwaterData.d18OWater < -4]
DripwaterData = pd.merge(DripwaterData, FieldTrip, how='left',on='idFieldTrip')

#%%Fig 8

fig8 = plt.figure(figsize=(7,5.6))
ax = fig8.add_subplot(111)      
to_regress = JinapsanCaveResults[(JinapsanCaveResults['mismatch'].isnull()==False) & (JinapsanCaveResults.DripsPerMin.isnull() == False)]
to_regress = to_regress[to_regress.DripsPerMin > 0]

text_y = max(to_regress['mismatch'])
text_x = max(to_regress['DripsPerMin'])*0.5
slope, intercept, r_value, p_value, std_err = stats.linregress(to_regress.DripsPerMin, to_regress['mismatch'])
yp = polyval([slope,intercept], to_regress.DripsPerMin)
for i in range(0,len(plot_table)):
    plotname = plot_table.DripName.iloc[i]
    plotmarker = plot_table.Marker.iloc[i]
    plotcolor = plot_table.Color.iloc[i]
    set_to_plot = to_regress[to_regress.SiteName == plotname]
    ax.plot(set_to_plot['DripsPerMin'], set_to_plot['mismatch'],
                 linestyle = 'none',
                 color = plotcolor,
                 marker = plotmarker,
                 markeredgewidth = 1,
                 linewidth = 1,
                 markeredgecolor = plotcolor,
                 label = plotname,
                 markersize = 5)
#for label in ax2.xaxis.get_ticklabels()[::2]:
#    label.set_visible(False)

ax.tick_params(axis='both', which='major', labelsize=15)
for label in ax3.yaxis.get_ticklabels():
    label.set_visible(False)
plot(to_regress.DripsPerMin, yp, 'k-')
tmp_ticks_labels = [10, 20, 30, 60, 120, 300]
tmp_ticks_pos = [60.0/y for y in tmp_ticks_labels]
ax2 = ax.twiny()
ax2.set_xlim([0,7])
ax2.set_xticks(tmp_ticks_pos)
ax2.set_xticklabels(tmp_ticks_labels)
ax2.tick_params(axis='x', labelsize=11)
ax2.set_xlabel(u'Drip interval\n (seconds)', fontsize=15)
ax.set_xlabel(u'Drip Rate  \n(drip/min)', fontsize=15, labelpad=15)
ax.set_ylabel(r'$1000\ln\alpha$ error', fontsize=15)
text(text_x,text_y-0.05,'$r^{2} =$ %5.2f' %(r_value**2), fontsize=15, ha='center')
text(text_x,text_y-0.2,'$p =$ %.1E' %(p_value), fontsize=15, ha='center')
plt.savefig('mismatchDriprate.png', format='png', dpi=600, bbox_inches='tight')

#%%Fig 9 1000lnalpha error

plot9_min = -2#min(pd.concat(JinapsanCaveResults['Modeled_Appar_Frac'], JinapsanCaveResults['Expect_Frac'], JinapsanCaveResults['1000lnalpha']))
plot9_max = 1#min(pd.concat(JinapsanCaveResults['Modeled_Appar_Frac'], JinapsanCaveResults['Expect_Frac'], JinapsanCaveResults['1000lnalpha']))

fig, axes = plt.subplots(nrows=2,ncols=2,sharex=True, figsize=(12, 6)) #figsize=(12, 2.7)
for i, ax in enumerate(axes.flatten()):
#        for j in range(0,len(grid_lines)):
#            ax.axvline(grid_lines[j], color='0.75', linestyle='--', lw=1)
    plotname = plot_table.DripName.iloc[i]
    plotmarker = plot_table.Marker.iloc[i]
    plotcolor = plot_table.Color.iloc[i]
    set_to_plot = JinapsanCaveResults[JinapsanCaveResults.SiteName == plotname]
    set_to_plot = set_to_plot.sort_values('Midpoint')


    ax.plot(set_to_plot['Midpoint'], set_to_plot['mismatch'],
                 linestyle = '-',
                 color = plotcolor,
                 marker = plotmarker,
                 markeredgewidth = 1,
                 linewidth = 1,
                 markeredgecolor = plotcolor,
                 label = plotname,
                 markersize = 5)

    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    for label in ax.yaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    ax.set_title(plotname)
    plot_ylim = ax.get_ylim()
    if plot_ylim[0] >=0:
        plot_ylim = (plot_ylim[0] - .1*plot_ylim[1], plot_ylim[1])
    ax.set_ylim([plot9_min, plot9_max])
    ax.set_xlim([start_plot, end_plot])
#ax0.set_xlim([start_plot, end_plot])
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    if i == 1 or i == 3:
        for label in ax.yaxis.get_ticklabels():
            label.set_visible(False)
    else:
        ax.set_ylabel(r'$1000\ln\alpha$ error', fontsize=15)          
plt.savefig('mismatchTime.png', format='png', dpi=600,bbox_inches='tight')

#%% Fig 10

DripwaterData = DripwaterData[DripwaterData.SiteName.isin(DripNames)]
DripwaterData = DripwaterData[DripwaterData.d18OWater < -4]
DripwaterData = pd.merge(DripwaterData, FieldTrip, how='left',on='idFieldTrip')


fig10 = plt.figure(figsize=(12,5.6))
ax1 = fig10.add_subplot(121)      
ax3 = fig10.add_subplot(122)  

for i in range(0,len(plot_table)):
    plotname = plot_table.DripName.iloc[i]
    plotmarker = plot_table.Marker.iloc[i]
    plotcolor = plot_table.Color.iloc[i]
    set_to_plot = to_regress[to_regress.SiteName == plotname]
    AveDripRate = mean(set_to_plot['DripsPerMin'])
    DripNormVar = cov(set_to_plot['DripsPerMin']/AveDripRate)
    ax1.plot(AveDripRate, DripNormVar,
                 linestyle = 'none',
                 color = plotcolor,
                 marker = plotmarker,
                 markeredgewidth = 1,
                 linewidth = 1,
                 markeredgecolor = plotcolor,
                 label = plotname,
                 markersize = 5)



    AveDripd18O = mean(set_to_plot['d18OWater'])
    NormVar = cov(set_to_plot['d18OWater']/AveDripd18O)
    slope, intercept, r_value, p_value, std_err = stats.linregress(set_to_plot['Midpoint'].map(datetime.datetime.toordinal), set_to_plot.DripsPerMin)
    yp = polyval([slope,intercept], set_to_plot['Midpoint'].map(datetime.datetime.toordinal))
    DetrendNormVar = cov([set_to_plot.DripsPerMin-yp])/AveDripRate    
    print(plotname, DripNormVar)
    print(plotname, p_value)
    print(plotname, DetrendNormVar)
    print(plotname, min(set_to_plot['DripsPerMin']))
    print(plotname, AveDripRate)
    print(plotname, max(set_to_plot['DripsPerMin']))
    ax3.plot(AveDripRate, NormVar,
                 linestyle = 'none',
                 color = plotcolor,
                 marker = plotmarker,
                 markeredgewidth = 1,
                 linewidth = 1,
                 markeredgecolor = plotcolor,
                 label = plotname,
                 markersize = 5)
#for label in ax2.xaxis.get_ticklabels()[::2]:
#    label.set_visible(False)
#for label in ax2.xaxis.get_ticklabels()[::2]:
#    label.set_visible(False)
ax1.tick_params(axis='both', which='major', labelsize=15)
tmp_ticks_labels = [10, 20, 30, 60, 120, 300]
tmp_ticks_pos = [60.0/y for y in tmp_ticks_labels]
ax2 = ax1.twiny()
ax2.set_xlim([0,7])
ax2.set_xticks(tmp_ticks_pos)
ax2.set_xticklabels(tmp_ticks_labels)
ax2.tick_params(axis='x', labelsize=11)
ax2.set_xlabel(u'Drip interval\n (seconds)', fontsize=15)
ax1.set_xlabel(u'Mean Drip Rate  \n(drips/min)', fontsize=15, labelpad=15)
ax1.set_ylabel(u'Normalized Drip Rate Variance', fontsize=15)
ax3.tick_params(axis='both', which='major', labelsize=15)
tmp_ticks_labels = [10, 20, 30, 60, 120, 300]
tmp_ticks_pos = [60.0/y for y in tmp_ticks_labels]
ax4 = ax3.twiny()
ax4.set_xlim([0,7])
ax4.set_xticks(tmp_ticks_pos)
ax4.set_xticklabels(tmp_ticks_labels)
ax4.tick_params(axis='x', labelsize=11)
ax4.set_xlabel(u'Drip interval\n (seconds)', fontsize=15)
ax3.set_xlabel(u'Mean Drip Rate  \n(drips/min)', fontsize=15, labelpad=15)
ax3.set_ylabel(r'$Normalized$ Drip$ \delta^{18}O$ Variance', fontsize=15)
plt.savefig('DripVar.png', format='png', dpi=600, bbox_inches='tight')