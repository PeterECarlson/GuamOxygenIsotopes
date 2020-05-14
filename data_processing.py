# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 13:28 2015

@author: anoronha
"""

import sqlite3
import datetime
from pylab import *
import pandas as pd
import numpy as np
from math import *
from pandas import DataFrame
import matplotlib.pyplot
from matplotlib import cm
from scipy import stats, signal
from growth_rate_model import calc_co2, generate_input_files

molar_mass_Ca = 40.08
molar_mass_Ca = 40.08

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

#raingauge = pd.read_csv('raingauge.csv')
#raingauge['Time'] = pd.to_datetime(raingauge['Time'])
#plot(raingauge.Time, raingauge.ALN,'g-')
#plot(raingauge.Timee, raingauge.BFH,'k-')

affek = pd.read_csv('affek2014.csv')

DripNames = ['Flatman','Station 2', 'Station 1','Stumpy','Trinity','Amidala','Stumpys Brother']

def rgb_convert(palette):
    for i in range(len(palette)):
        r, g, b = palette[i]
        palette[i] = (r / 255., g / 255., b / 255.)
    return palette

plot_table = DataFrame(columns = ['DripName','Color','Marker'])
plot_table.DripName = DripNames
plot_table.Color = rgb_convert([(27,158,119), (217,95,2), (117,112,179), (231,41,138), (102,166,30), (230,171,2), (166,118,29)])
plot_table.Marker = ['o','D','v','^','s','<','>']
plot_table = plot_table.sort_values('DripName')

def plot_all(series, to_plot, x_val, y_val, fill,line, axis_num, y_error_bar,y_direct):
    fig1 = plt.figure(num=1, figsize = (11,8.5), facecolor='w', edgecolor='k')
    ax1 = fig1.add_subplot(111)

    if axis_num == 2:
        ax1 = ax1.twinx()

    if y_direct == 'reverse':
        ax1.invert_yaxis()

    plots = plot_table[plot_table.DripName == to_plot].reset_index(drop=True)
    if to_plot == 'all':
        plots = plot_table

    for i in range(0,len(plots)):
        plotname = plots.DripName[i]
        plotmarker = plots.Marker[i]
        plotcolor = plots.Color[i]
        set_to_plot = series[series.SiteName == plotname]
        set_to_plot = set_to_plot.sort_values(x_val)
        ax1.plot(set_to_plot[x_val].values, set_to_plot[y_val].values,
                 linestyle = line,
                 color = plotcolor,
                 marker = plotmarker,
                 markeredgewidth = 2,
                 markeredgecolor = plotcolor,
                 label = plotname,
                 markerfacecolor = fill,
                 markersize = 7)
        if type(y_error_bar) == float:
            ax1.errorbar(set_to_plot[x_val], set_to_plot[y_val],
                         yerr = y_error_bar,
                         ecolor = plotcolor,
                         linestyle = 'none')
    xlabel(x_val)
    ylabel(y_val)
    legend(loc = 'best', numpoints = 1)



def match_calcite_water(calcite, water, param):
    corr_trip_rec = []
    water_rec = []
    param_mean_rec = []
    for i in range(0, len(calcite)):
        test_site = calcite.SiteName.iloc[i]
        corr_trip = list(range(int(calcite.DeployFieldTrip.iloc[i]), int(calcite.CollectFieldTrip.iloc[i])))
        corr_trip.append(int(calcite.CollectFieldTrip.iloc[i]))
        corr_trip_rec.append(list(corr_trip))
        water_values = water[water['idFieldTrip'].isin(corr_trip) & ((water.SiteName == test_site) == True)]
        water_rec.append(list(water_values.idFieldTrip))
        corr_water_param = float(np.mean(water_values[param]))
        param_mean_rec.append(corr_water_param)
    calcite[param] = param_mean_rec
    return calcite


def lnalpha(series):
    results = []
    for i in range(0,len(series)):
        if (isnan(series.d18OCalcite.iloc[i]) == False) & (isnan(series.d18OWater.iloc[i]) == False):
            result = 1000*(log(((series.d18OCalcite.iloc[i]/1000)+1)/((series.d18OWater.iloc[i]/1000)+1)))
            results.append(result)
        else:
            results.append(float('Nan'))
    series['1000lnalpha'] = results
    return series

def fill_temp(series,DripChem):
    results = []
    DripChem_mean = []
    trips_list = list(set(DripChem.idFieldTrip))
    for i in range(0,len(trips_list)):
        corr_trip_set = DripChem[DripChem.idFieldTrip == trips_list[i]]
        DripChem_mean.append(np.mean(corr_trip_set.Temp))
    trips_list = pd.Series(trips_list, name='idFieldTrip')
    DripChem_mean = pd.Series(DripChem_mean, name='AverageWaterTemp')
    DripChem_mean = pd.concat([trips_list, DripChem_mean], axis = 1)
    for i in range(0,len(series)):
        if (isnan(series.Temp.iloc[i]) == True):
            corr_trip = list(range(int(series.DeployFieldTrip.iloc[i]), int(series.CollectFieldTrip.iloc[i])))
            corr_trip.append(int(series.CollectFieldTrip.iloc[i]))
            corr_mean_temp = DripChem_mean[DripChem_mean['idFieldTrip'].isin(corr_trip)]
            results.append(list(corr_mean_temp.idFieldTrip))
            series.Temp.iloc[i] = float(np.mean(corr_mean_temp.Temp))
        else:
            results.append(float('Nan'))
    series['WaterTempFlag'] = results
    return series

CalciteData = match_calcite_water(CalciteData, DripwaterData, 'd18OWater')
CalciteData = match_calcite_water(CalciteData, DripwaterData, 'Temp')
CalciteData = match_calcite_water(CalciteData, DripwaterData, 'pH')
CalciteData = match_calcite_water(CalciteData, DripRate, 'DripInterval')
CalciteData = match_calcite_water(CalciteData, DripRate, 'DripsPerMin')
CalciteData = match_calcite_water(CalciteData, DripRate, 'MlPerMin')
CalciteData = match_calcite_water(CalciteData, CationResults, 'Ca')

holder = []
for i in range(0,len(CalciteData)):
    if np.isnan(CalciteData.Temp.iloc[i]):
        holder.append(26.5)
    else:
        holder.append(CalciteData.Temp.iloc[i])
CalciteData['Temp'] = holder
#CalciteData = avg_co2(CalciteData, FieldCO2Continous)
JinapsanCaveResults = lnalpha(CalciteData)

JinapsanCaveResults = JinapsanCaveResults.reset_index(drop = False)
input_files, index_list = generate_input_files(JinapsanCaveResults)
results_table = calc_co2(input_files, index_list, JinapsanCaveResults)
results_table['index']=results_table['index'].apply(int)
JinapsanCaveResults = pd.merge(JinapsanCaveResults, results_table, how='left', on='index')

JinapsanCaveResults['Ca'] = [x/(1e6 * molar_mass_Ca) for x in JinapsanCaveResults.Ca]
JinapsanCaveResults['Delta_Ca'] = JinapsanCaveResults.Ca - JinapsanCaveResults.Ca_eq
#JinapsanCaveResults['ln_interval'] = [log(x) for x in JinapsanCaveResults.DripInterval]
JinapsanCaveResults.to_pickle('JinapsanCaveResults.pickle')

def match_modeled_meas(meas, modeled, param):
    corr_trip_rec = []
    modeled_rec = []
    param_mean_rec = []
    for i in range(0, len(meas)):
        test_site = meas.SiteName.iloc[i]
        corr_trip = list(range(int(meas.DeployFieldTrip.iloc[i]),
                               int(meas.CollectFieldTrip.iloc[i])))
        corr_trip.append(int(meas.CollectFieldTrip.iloc[i]))
        corr_trip_rec.append(list(corr_trip))
        modeled_values = modeled[modeled['idFieldTrip'].isin(corr_trip) &\
                               ((modeled.SiteName == test_site) == True)]
        modeled_rec.append(list(modeled_values.idFieldTrip))
        corr_mod_param = float(np.mean(modeled_values[param]))
        param_mean_rec.append(corr_mod_param)
    meas[param] = param_mean_rec
    return meas


FieldCO2Continous = FieldCO2Continous.sort_values('Datetime')
FieldCO2Continous = FieldCO2Continous.set_index(pd.DatetimeIndex(FieldCO2Continous['Datetime']))
FieldCO2_Daily = FieldCO2Continous.resample('D').mean()
FieldCO2_Daily = FieldCO2_Daily.reset_index(drop = False)
FieldCO2_Daily.rename(columns={'index': 'Date'}, inplace=True)
FieldCO2_Daily = FieldCO2_Daily[np.isfinite(FieldCO2_Daily['CO2'])]

Climate = Climate.sort_values('DATE')
Climate = Climate.set_index(pd.DatetimeIndex(Climate['DATE']))
Climate_Monthly = Climate.resample('M').mean()
Climate_Monthly = Climate_Monthly.reset_index(drop = False)
Climate_Monthly.rename(columns={'index': 'DATE'}, inplace=True)


FieldCO2_Monthly = FieldCO2Continous.set_index(pd.DatetimeIndex(FieldCO2Continous['Datetime']))
FieldCO2_Monthly = FieldCO2_Monthly.resample('M').mean()
FieldCO2_Monthly = FieldCO2_Monthly.reset_index(drop = False)
FieldCO2_Monthly.rename(columns={'index': 'Date'}, inplace=True)

ExternalMSResults = pd.read_csv('externalms.csv')
ExternalMSResults.rename(columns={'1000lnalpha': '1000lnalpha_pub'}, inplace=True)
ExternalMSResults = lnalpha(ExternalMSResults).sort_values('Cave')
ExternalMSResults = ExternalMSResults[ExternalMSResults.Category == 'Cave Calcite']
list_ExternalCaves = list(set(ExternalMSResults.Cave))
list_ExternalCaves.sort()


fig1 = plt.figure(figsize = (10,5.63), facecolor='w', edgecolor='k')
ax1 = fig1.add_subplot(111)
ax1.set_ylim([25,38])
ax1.set_xlim([3.2,3.65])
ax1.plot((1000/((JinapsanCaveResults.Temp)+273.15)), JinapsanCaveResults['1000lnalpha'],
         linestyle='None', 
         markersize = 5, 
         color = '0.5', 
         marker='o',
         label='This Study')
ax1.plot((1000/((ExternalMSResults.WaterTemp)+273.15)),ExternalMSResults['1000lnalpha'],
     marker = 'o',
     markersize = 7,
     markeredgecolor = 'k',
     markerfacecolor = 'none',
     markeredgewidth = 1,
     linestyle = 'none',
     label='Published Measurements')
#    ax1.plot((1000/((affek.Temp)+273.0)), affek['1000lnalpha'],'ko')
Temp_range = [-5,50]
Temp_range = [1000/(i+273.15) for i in Temp_range]
Coplen_alpha = [17.4*(i)-28.6 for i in Temp_range]
Affek_alpha = [15.63*(i)-23.29 for i in Temp_range]
KimONeil_alpha = [18.03*(i)-32.42 for i in Temp_range]
ax1.plot(Temp_range, 
         KimONeil_alpha, 
         linestyle=':', 
         linewidth=3,
         color = '#e41a1c',
         label = 'Kim and O\'Neil (1997)')
ax1.plot(Temp_range, 
         Coplen_alpha, 
         linestyle='-', 
         linewidth=3,
         color='#377eb8',
         label = 'Coplen (2002)')
ax1.plot(Temp_range, 
         Affek_alpha, 
         linestyle='--',
         linewidth=3,
         color='#4daf4a',
         label = 'Affek and Zaarur (2014)')
ax1.set_xlim([3.3,3.65])
ax1.set_ylim([27,36])
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=20)
for label in ax1.yaxis.get_ticklabels()[::2]:
    label.set_visible(False)
xlabel(r'1000/T (K)', fontsize=20)
#    box = ax1.get_position()
#    ax1.set_position([box.x0, box.y0, box.width * 0.5, box.height])
ax1.legend(loc='upper left',numpoints = 1, frameon=False)
ax2 = ax1.twiny()
#    ax2.set_position([box.x0, box.y0, box.width * 0.5, box.height])
tmp_ticks_labels = [25, 20, 15, 10, 5, 0]
tmp_ticks_pos = [1000.0/(x+273.0) for x in tmp_ticks_labels]
ax2.set_xlim([3.3,3.65])
#    tmp_pos = np.array([3.2, 3.4, 3.5, 3.6])
ax2.set_xticks(tmp_ticks_pos)
ax2.set_xticklabels(tmp_ticks_labels)
ax2.tick_params(axis='x', labelsize=20)
ax2.set_xlabel(u'T (\u00b0C)', fontsize=20)
ax1.set_ylabel(r'$1000\ln\alpha$', fontsize=20)
savefig('AGU2016_Fig1.png', format='png', dpi=1000, bbox_inches='tight')


def plot_cross(series, site, symbol, param1, param2):
    if (site == 'all') == False:
        series = series[series.SiteName == site]
    to_regress = series[(series[param1].isnull()==False) & (series[param2].isnull() == False)]
    slope, intercept, r_value, p_value, std_err = stats.linregress(to_regress[param1], to_regress[param2])
    yp = polyval([slope,intercept], to_regress[param1])
    plot(series[param1], series[param2], symbol)
    plot(to_regress[param1], yp, 'r-')
    title(site)
    xlabel(param1)
    ylabel(param2)
    x_cord = min(to_regress[param1])
    y_cord = max(to_regress[param2])
    text(x_cord,y_cord,'y = %5.2fx+%5.2f, r^2 = %5.2f, pval = %5.5f' %(slope, intercept, r_value, p_value))


#plot_query = input("Plot the feng growth rate stuff? : ")
#if plot_query == 'y':
#    feng = ExternalMSResults[(ExternalMSResults.PublicationKey == 'Feng2014') | (ExternalMSResults.PublicationKey == 'Feng2012')]
#
#    holder = []
#    for i in range(0,len(feng)):
#        if feng.GrowthRate.iloc[i]:
#            holder.append(0.000001)
#        else:
#            holder.append(feng.GrowthRate.iloc[i])
#    feng['Growth_tmp'] = holder
#    feng['SA'] = ((feng.GrowthRate*0.27)/2.0)
#    #feng['Rate_tmp'] = feng.Growth_tmp/86400.0
#    feng['mols_calcite'] = feng.GrowthRate*(1.0/100.09)
#    feng['R'] = (feng['mols_calcite']/feng['SA'])/86400.0
#
#    def take_log(series, param):
#        holder = []
#        for i in range(0,len(series)):
#            if (isnan(series[param].iloc[i])) == False and (series[param].iloc[i] > 0):
#                holder.append(log10(series[param].iloc[i]))
#            else:
#                holder.append(float('nan'))
#        col_label = 'log_%s' %(param)
#        series[col_label] = holder
#        return series
#
#
#    from decimal import *
#    dummy_guam = pd.concat([JinapsanCaveResults.Temp, JinapsanCaveResults['1000lnalpha'], JinapsanCaveResults.GrowthRate], axis=1)
#    #dummy_guam = take_log(dummy_guam, 'GrowthRate')
#    dummy_guam['source'] = 1
#    dummy_tx = pd.concat([feng.WaterTemp, feng['1000lnalpha'], feng.GrowthRate], axis=1)
#    #dummy_tx = take_log(dummy_tx, 'GrowthRate')
#    dummy_tx['source'] = 0
#    dummy_tx.rename(columns={'WaterTemp': 'Temp'}, inplace=True)
#    dummy = dummy_guam.append(dummy_tx)
#    dummy = dummy[np.isfinite(dummy['GrowthRate'])]
#    dummy = dummy[dummy.GrowthRate > 0]
#    dummy['Growthtmp'] = (dummy.GrowthRate*(1.0/100.09))/86400.0
#    dummy = dummy[np.isfinite(dummy['Temp'])]
#    dummy['SA'] = (dummy.Growthtmp*0.27)/2.0
#    dummy['mols_calcite'] = dummy.Growthtmp*(1.0/100.09)
#    holder = []
#    for i in range(0,len(dummy)):
#        tmp = (Decimal(dummy.mols_calcite.iloc[i])/(Decimal(dummy.SA.iloc[i])))/Decimal(8640000.0)
#        hodler = holder.append(tmp)
#    dummy['R'] = holder
#    dummy = take_log(dummy, 'R')
#    dummy = take_log(dummy, 'Growthtmp')
#
#    tmp = ((1000/(dummy.Temp+273.0))*17.4)-28.6
#    dummy['alpha_detrended'] = dummy['1000lnalpha']-tmp
#
#    #dummy = dummy[dummy.source == 1]
#
#    fig2 = plt.figure(num=1, figsize = (11,8.5), facecolor='w', edgecolor='k')
#    ax1 = fig2.add_subplot(111)
#    ax = fig2.gca()
#    ax.scatter(dummy.log_Growthtmp,dummy['alpha_detrended'], c = dummy['Temp'], s=100)
#    #Temp_range = [-5,50]
#    #Temp_range = [1000/(i+273.0) for i in Temp_range]
#    #Coplen_alpha = [17.4*(i)-28.6 for i in Temp_range]
#    #Tremaine_alpha = [16.1*(i)-24.6 for i in Temp_range]
#    #KimONeil_alpha = [18.03*(i)-32.17 for i in Temp_range]
#    #ax.plot(Temp_range, Coplen_alpha, 'k-', label = 'Coplen Line')
#    ax.set_xlim([3.2,3.6])
#    ##ax.plot(Temp_range, Tremaine_alpha, 'k--', label = 'Tremaine Line')
#    #ax.plot(Temp_range, KimONeil_alpha, 'k:', label = 'Kim and ONeil Line')
#    xlabel(r'lnGrowthRate')
#    ax.set_ylabel(r'$1000\ln\alpha$')
#    box = ax.get_position()
#    ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
#    ax2 = ax.twiny()
#    ax2.set_xlim([3.2,3.6])
#    ax2.set_position([box.x0, box.y0, box.width * 0.75, box.height])
#    tmp_ticks_labels = [30, 20, 10, 0]
#    tmp_ticks_pos = [1000.0/(x+273.0) for x in tmp_ticks_labels]
#    tmp_pos = np.array([3.2, 3.4, 3.5, 3.6])
#    ax2.set_xticks(tmp_ticks_pos)
#    ax2.set_xticklabels(tmp_ticks_labels)
#    ax2.set_xlabel(u'T (\u00b0C)')
#
#    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=False, shadow=False, ncol=2,numpoints = 1,prop={'size':10})
#
#    ax.scatter(1000/(feng['WaterTemp']+273.0),feng['1000lnalpha'], c = feng.GrowthRate, s=100)
#    #s.set_clim([cmin,cmax])
##    cb = fig.colorbar(s)
##    cb.set_label('log(Growth Rate (g/day))')
#    del plot_query
#else:
#    del plot_query





