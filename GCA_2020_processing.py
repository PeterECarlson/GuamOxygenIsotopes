# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 20:05:12 2020

@author: peter
"""

import pandas as pd
import numpy as np
import datetime

import oxygen_isotope_stats_functions as oxystats

calcite_data = oxystats.load_cave_data()
drip_data = oxystats.load_drip_data()
external_data = oxystats.load_external_data()

cap18O = '${\Delta}^{18}O$'
d18O = '${\delta}^{18}O$'
d13C = '${\delta}^{13}C$'
permil = '$\u2030$'
# %% Initialize
    

FieldCO2Continous = pd.read_pickle('FieldCO2Continous.pickle')
FieldCO2Spot = pd.read_pickle('FieldCO2Spot.pickle')
DripNames = ['Flatman','Station 2', 'Station 1','Stumpy']

def t_to_k(t, t_unit='K'):
    assert(t_unit.lower() in ('k', 'c'))
    if t_unit.lower() == 'k':
        pass
    elif t_unit.lower() == 'c':
        t = t + 273.15
    return t

def affek_1000lna(t, t_unit='K'):
    t = t_to_k(t, t_unit)
    frac = 15630 / t - 23.29
    return frac

def coplen_1000lna(t, t_unit='K'):
    t = t_to_k(t, t_unit)
    frac = 17400 / t - 28.6
    return frac

def kim_1000lna(t, t_unit='K'):
    t = t_to_k(t, t_unit) 
    frac = 18030 / t - 32.42
    return frac

def tremaine_1000lna(t, t_unit='K'):
    t = t_to_k(t, t_unit)  
    frac = 16100 / t - 24.6
    return frac
      
def rgb_convert(palette):
    for i in range(len(palette)):
        r, g, b = palette[i]
        palette[i] = (r / 255., g / 255., b / 255.)
    return palette


dripintervals = []
dripinterval_errs = []
dripinterval_ns = []
drip_d18O_agg_mean = []
drip_d18O_agg_err = []
drip_d18O_err = []
for i, row in calcite_data.iterrows():
    site = row['SiteName']
    start = row['DateDeploy'] - pd.to_timedelta(2, unit='d')
    end = row['DateCollect'] + pd.to_timedelta(2, unit='d')
    drip = drip_data[((drip_data['SiteName'] == site)
                      & (drip_data['Date'] >= start)
                      & (drip_data['Date'] <= end))]
    est_s_per_drip = 60 / (drip['BottleFillRate'] / 0.07)
    drip['DripInterval'] = drip['DripInterval'].mask(
        drip['DripInterval'].isnull(),
        est_s_per_drip)    
    dripintervals.append(drip['DripInterval'].mean())
    dripinterval_errs.append(drip['DripInterval'].std())
    dripinterval_ns.append(len(drip['DripInterval']))
    drip_d18O_agg_mean.append(drip['δ18OWaterVSMOW'].mean())
    drip_d18O_agg_err.append(drip['δ18OWaterVSMOW'].std())
    drip_d18O_err.append(drip['δ18OWaterErr'].mean())
    
calcite_data['DripInterval'] = dripintervals
calcite_data['DripInterval_err'] = dripinterval_errs
calcite_data['DripInterval_n'] = dripinterval_ns
calcite_data['DripsPerMin'] = 60 / calcite_data['DripInterval']
calcite_data['DripRate_err'] = (60*calcite_data['DripInterval_err']
                                / calcite_data['DripInterval'].pow(2))
calcite_data['δ18OWaterVSMOW'] = drip_d18O_agg_mean
calcite_data['δ18OWater_Err_VSMOW'] = drip_d18O_agg_err
calcite_data['δ18OWater_Meas_Err_VSMOW'] = drip_d18O_err
calcite_data['δ18OWater_Err_VSMOW'] = calcite_data['δ18OWater_Err_VSMOW'].mask(
    calcite_data['δ18OWater_Err_VSMOW'].isnull(),
    calcite_data['δ18OWater_Meas_Err_VSMOW'])
print(calcite_data['δ18OWater_Err_VSMOW'])
calcite_data['δ18OWaterPDB'] = (calcite_data['δ18OWaterVSMOW'] - 30.91) / 1.03091
calcite_data['δ18OWater_Err_PDB'] = calcite_data['δ18OWater_Err_VSMOW'] / 1.03091
calcite_data['Δ18OPDB'] = (calcite_data['Calcite_δ18O_VPDB']
                           - calcite_data['δ18OWaterPDB'])
calcite_data['Δ18OPDB_Err'] = (calcite_data['Calcite_δ18O_Err_VPDB'].pow(2)
                               + calcite_data['δ18OWater_Err_PDB'].pow(2)
                               ).pow(0.5)
calcite_data['frac_kim'] = kim_1000lna(calcite_data['WaterTemp'], 'C')
calcite_data['frac_affek'] = affek_1000lna(calcite_data['WaterTemp'], 'C')
calcite_data['frac_coplen'] = coplen_1000lna(calcite_data['WaterTemp'], 'C')
calcite_data['alpha_kim'] = np.exp(calcite_data['frac_kim'] / 1000)
calcite_data['alpha_affek'] = np.exp(calcite_data['frac_affek'] / 1000)
calcite_data['alpha_coplen'] = np.exp(calcite_data['frac_coplen'] / 1000)
calcite_data['equilib_calcite_δ18O_kim'] = ((calcite_data['δ18OWaterPDB'] + 1000)
                                            * calcite_data['alpha_kim'] - 1000)
calcite_data['equilib_calcite_δ18O_affek'] = ((calcite_data['δ18OWaterPDB'] + 1000)
                                            * calcite_data['alpha_affek'] - 1000)
calcite_data['equilib_calcite_δ18O_coplen'] = ((calcite_data['δ18OWaterPDB'] + 1000)
                                            * calcite_data['alpha_coplen'] - 1000)


calcite_data.to_csv('calcite_reprocessed.csv')
FieldCO2Spot.to_csv('field_co2_spot.csv')
FieldCO2Continous.to_csv('field_co2_continous.csv')
