# -*- coding: utf-8 -*-
"""
Created on Thu May 12 17:19:24 2016

@author: anoronha

modified from Matlab function written by mdeininger
"""

import pandas as pd
from pandas import DataFrame
import numpy as np
import os
from pylab import *
import datetime
import phreeqpy.iphreeqc.phreeqc_dll as phreeqc_mod
from my_functions import plot_all


datpath = r"\Program Files\USGS\IPhreeqcCOM 3.4.0-12927\database"

DreybrodtModel = pd.read_pickle('DreybrodtModel.pickle')
#DreybrodtModel['CO2'] = 500
PlateMass = pd.read_pickle('PlateMass_all.pickle')

def phreeqc_equil_co2_input_files(drip_data):
    input_files = []
    index_list = []
    selected_output = 'SELECTED_OUTPUT\n-molalities Ca+2 HCO3-\n-activities H+ Ca+2 CO2 HCO3- CO3-2'
    drip_data = drip_data[(np.isnan(drip_data.Temp) == False) & \
                          (np.isnan(drip_data.CO2) == False)]
    for i in range(0,len(drip_data)):
        index_rec = drip_data['index'].iloc[i]
        solution = 'SOLUTION %s\ntemp %5.2f\nunits ppm\n     8.0\nCa %5.2f\nNa %5.2f\nMg %5.2f' \
                    %(i, drip_data.Temp.iloc[i], drip_data.Na.iloc[i]/1000, drip_data.Ca.iloc[i]/1000, drip_data.Mg.iloc[i]/1000)
#        solution = 'SOLUTION %s\ntemp %5.2f\npH 8.0' %(i, drip_data.Temp.iloc[i])
        eql_phase = 'EQUILIBRIUM_PHASES\nCO2(g) %5.2f\nCalcite 0.0' \
                    %(log10(drip_data.CO2.iloc[i]/(10**6)))
        ind_file = '%s\n%s\n%s\nEND' %(solution, eql_phase, selected_output)
        input_files.append(ind_file)
        index_list.append(index_rec)
    return input_files, index_list

def phreeqc_calc_co2(input_files, index_list):
    def selected_array(db_path, input_string):
        dbase = phreeqc_mod.IPhreeqc(r"\Users\peter\Documents\PYTHON\vs2018_64_dll\Release\IPhreeqc.dll")
        dbase.load_database(db_path)
        dbase.run_string(input_string)
        return dbase.get_selected_output_array()
    results_table = pd.DataFrame(columns = ['index', 'Ca_eq'])
    for i in range(0,len(input_files)):
        result = selected_array(os.path.join(datpath, 'phreeqc.dat'), input_files[i])
        dicts = {'index': index_list[i], 'Ca_eq': result[2][8]}
        results_table = results_table.append([dicts])
    return results_table


def growth_rate(drip_interval, drip_vol, temp, drip_ca, ca_eq, ca_ex_flag):
    #g/L/(g/mol) = mol/L
    drip_ca = (drip_ca/40.08)
    #(mol/L - mol/L) * 1 L/0.001 m^3 = mol/m^3
    ca_ex = (drip_ca - ca_eq) * 1e3
    if ca_ex < 0:
        ca_ex = 0
    if ca_ex_flag == 'Y':
        ca_ex = 0.2
    #m, film thickness
    delta = 8e-4
    #ml * 1 m^3/1e6 mL
    drip_vol = drip_vol * 1e-6
    #m/s, Equation from Hansen2013 pg 244
    alpha = (0.52 + 0.04 * temp + 0.004 * temp**2) * 1e-3
    #m/(m/s) = s
    Z = delta / alpha
    #g/mol * m/s * (1 - exp(-s/s)) * mol/m^3 = g/(m^2*s), Growth rate
    W_0 = 100.09 * delta / drip_interval * (1 - exp(-drip_interval / Z)) * ca_ex
    #h/day * min/day * s/day = s
    seconds_per_day = 24 * 60 * 60
    #m^3/m = m^2
    area = drip_vol/delta
    #g/(m^2*s) * m^2 * s/day = g/day
    growth_rate = W_0 * area * seconds_per_day
    return growth_rate
 
DreybrodtModel = DreybrodtModel.reset_index(drop = False)
#tmp = DreybrodtModel.copy()
input_files, index_list = phreeqc_equil_co2_input_files(DreybrodtModel)
results_table = phreeqc_calc_co2(input_files, index_list)
results_table['index']=results_table['index'].apply(int)
DreybrodtModel = pd.merge(DreybrodtModel, results_table, how='left',on='index')
holder = []
for i in range(0,len(DreybrodtModel)):
    output = growth_rate(DreybrodtModel.DripInterval.iloc[i],
                              DreybrodtModel.DripVol.iloc[i],
                              DreybrodtModel.Temp.iloc[i],
                              DreybrodtModel.Ca.iloc[i]*1e-6,
                              DreybrodtModel.Ca_eq.iloc[i],'N')
    holder.append(output)
DreybrodtModel['GrowthRate_Modeled'] = holder

holder = []
for i in range(0,len(DreybrodtModel)):
    output = growth_rate(DreybrodtModel.DripInterval.iloc[i],
                              DreybrodtModel.DripVol.iloc[i],
                              DreybrodtModel.Temp.iloc[i],
                              DreybrodtModel.Ca.iloc[i]*1e-6,
                              DreybrodtModel.Ca_eq.iloc[i],'Y')
    holder.append(output)
DreybrodtModel['GrowthRate_Modeled_Cons'] = holder
#tmp['CO2'] = 500
#input_files, index_list = phreeqc_equil_co2_input_files(tmp)
#results_table = phreeqc_calc_co2(input_files, index_list)
#tmp = pd.merge(tmp, results_table, how='left',on='index')
#holder = []
#for i in range(0,len(tmp)):
#    output = growth_rate(tmp.DripInterval.iloc[i],
#                              tmp.DripVol.iloc[i],
#                              tmp.Temp.iloc[i],
#                              tmp.Ca.iloc[i]*1e-6,
#                              tmp.Ca_eq.iloc[i])
#    holder.append(output)
#DreybrodtModel['GrowthRate_Const_CO2'] = holder

#DreybrodtModel = DreybrodtModel[DreybrodtModel.GrowthRate_Modeled > 0]
PlateMass = PlateMass[PlateMass.GrowthRate > 0]


def match_modeled_meas(meas, modeled, param):
    param_mean_rec = []
    for i in range(0, len(meas)):
        corr_trip = []
        modeled_values =[]
        test_site = meas.SiteName.iloc[i]
        corr_trip = list(range(int(meas.DeployFieldTrip.iloc[i]),
                               int(meas.CollectFieldTrip.iloc[i])))
        corr_trip.append(int(meas.CollectFieldTrip.iloc[i]))
        modeled_values = modeled[modeled['idFieldTrip'].isin(corr_trip) & (( modeled.SiteName == test_site) == True)]
        corr_mod_param = float(np.mean(modeled_values[param]))
        param_mean_rec.append(corr_mod_param)
    meas[param] = param_mean_rec
    return meas
    
def standardize(series, param):
    sites = list(set(series.SiteName))
    holder = DataFrame(columns = ['SampleName', 'Standardized'])
    for i in range(0,len(sites)):
        site_set = series[series.SiteName == sites[i]]
        site_avg = np.mean(site_set[param])
        dev_holder = []
        for j in range(0,len(site_set)):
            dev_holder.append(site_avg - site_set[param].iloc[j])
        tmp = []
        tmp = DataFrame(columns = ['SampleName', 'Standardized'])
        tmp['SampleName'] = site_set.SampleName
        tmp['Standardized'] = dev_holder
        holder = pd.concat([holder, tmp])
    return holder
        
    

Growth = match_modeled_meas(PlateMass, DreybrodtModel, 'GrowthRate_Modeled')
Growth = match_modeled_meas(Growth, DreybrodtModel, 'Ca')
Growth = match_modeled_meas(Growth, DreybrodtModel, 'DripInterval')


#model_growth_standardized = standardize(DreybrodtModel, 'GrowthRate_Modeled')
#DreybrodtModel = pd.merge(DreybrodtModel, model_growth_standardized, how='left',on='SampleName')
#meas_growth_standardized = standardize(PlateMass, 'GrowthRate')
#PlateMass = pd.merge(PlateMass, meas_growth_standardized, how='left',on='SampleName')
#Growth = Growth[(Growth.Site == 'Trinity') == False]
#DreybrodtModel = DreybrodtModel[(DreybrodtModel.Site == 'Trinity') == False]
#PlateMass = PlateMass[(PlateMass.Site == 'Trinity') == False]
#Growth['Offset'] = Growth.GrowthRate_Modeled - Growth.GrowthRate


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
    
#plot_all(Growth, 'all', 'GrowthRate', 'GrowthRate_Modeled', 'w','none', 1, 'none','norm')
#plot([0,0.05], [0,0.05], 'k-')

Growth.to_pickle('ModeledvsMeasGrowth.pickle')
DreybrodtModel.to_pickle('ModeledGrowth.pickle')
Fig1 = plt.figure(num=1, figsize = (11,8.5), facecolor='w', edgecolor='k')
plot_all(DreybrodtModel, to_plot = 'Station 2', x_val = 'EndFieldTrip', y_val = 'GrowthRate_Modeled', fill ='k', line = '-',axis_num = 1, y_error_bar = 'none', y_direct ='norm')
Fig2 = plt.figure(num=1, figsize = (11,8.5), facecolor='w', edgecolor='k')
plot_all(PlateMass, 'Station 2', 'Midpoint', 'GrowthRate', 'w','-', 1, '-','norm')
