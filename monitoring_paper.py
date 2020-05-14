# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 11:19:39 2015

@author: alexandranoronha
"""

import datetime
import pylab 
import pandas as pd
import numpy as np
import numbers as nb
from math import *
from matplotlib.dates import date2num
from pandas import DataFrame
import collections
import matplotlib.pyplot as plt
import random

#load data as csv files, outputs from queries on database
samplename = pd.read_csv('samplename.csv')
dripwater = pd.read_csv('dripdata.csv')
dripwater = dripwater.sort_values(['Site','Trip']).reset_index(drop = True)
dripinterval = pd.read_csv('dripinterval.csv')
dripinterval['Datetime'] = pd.to_datetime(dripinterval['Datetime'])
bottledata = pd.read_csv('bottledata.csv')
bottledata['DeployTime'] = pd.to_datetime(bottledata.DeployTime)
bottledata['CollectTime'] = pd.to_datetime(bottledata.CollectTime)
bottledata = pd.merge(bottledata, samplename, how = 'left', on = 'SampleName')
climate = pd.read_csv('GUM_Daily_NOAACDO.csv')
climatedate = [datetime.datetime.strptime(str(i), '%Y%m%d') for i in climate['DATE']]
climate['DATE'] = pd.to_datetime(climatedate)
climateprecip = [i if i >= 0 else float('nan') for i in climate.PRCP ]
climate['PRCP'] = climateprecip

cdioxide = pd.read_csv('caveco2test.csv', names = ['idFieldCO2','Site','FieldInstrumentName','LogNumber','Datetime','CO2'])
cdioxide['Datetime'] = pd.to_datetime(cdioxide['Datetime'])
daily_cdioxide = cdioxide.set_index(pd.DatetimeIndex(cdioxide['Datetime']))
daily_cdioxide = daily_cdioxide.resample('D').mean()

plate_mass = pd.read_csv('platemass.csv')
plate_name = pd.read_csv('platedata.csv')
plate_isotope = pd.read_csv('plateisotope.csv')
trip = pd.read_csv('trip.csv')
ultrameter = pd.read_csv('ultrameter.csv').sort_values(['Site','Trip'])
ultrameter.rename(columns={'Temp': 'WaterTemp'}, inplace=True)
external = pd.read_csv('externalms.csv')
external.rename(columns={'1000lnalpha': '1000lnalpha_pub'}, inplace=True)


bottledata['MinutesDeployed'] = (bottledata.CollectTime-bottledata.DeployTime)
bottledata['MinutesDeployed'] = [(i.days*86400 + i.seconds)/60.0 for i in bottledata.MinutesDeployed]
bottledata['FlowRate'] = (bottledata.FinalMass - bottledata.InitialMass)/bottledata.MinutesDeployed

#group by PlateName to identify and average duplicate isotope measurements
plate_isotope = plate_isotope.groupby(['PlateName'])
singles = DataFrame()
platename = []
average = []
for name,group in plate_isotope:
    if group['PlateName'].count() > 1:
        platename.append(name)
        average.append(np.mean(group['Result']))
    else:
        singles = singles.append(group)
processed_dupes = DataFrame({'PlateName': platename, 'd18O_VPDB': average})
processed_singles = DataFrame({'PlateName': singles.PlateName, 'd18O_VPDB': singles.Result})
frames = [processed_dupes, processed_singles]
result = pd.concat(frames)
plate_isotope = result.reset_index(drop=True)
plate_isotope = pd.merge(plate_isotope, plate_name, how = 'left', on = 'PlateName')
plate_isotope['d18O_Calcite'] = (plate_isotope.d18O_VPDB*1.03091)+30.91 

#format dates in the plate table correctly
plate_mass['DateDeploy'] = pd.to_datetime(plate_mass['DateDeploy'])
plate_mass['DateRemove'] = pd.to_datetime(plate_mass['DateRemove'])
plate_mass['DateReplace'] = pd.to_datetime(plate_mass['DateReplace'])
plate_mass['DateCollect'] = pd.to_datetime(plate_mass['DateCollect'])
trip['BeginTrip'] = pd.to_datetime(trip['BeginTrip'])
trip['EndTrip'] = pd.to_datetime(trip['EndTrip'])

#remove plates that were either lost or are still in the field
plate_mass = plate_mass[plate_mass.DateDeploy.isnull()==False]
plate_mass = plate_mass[plate_mass.DateCollect.isnull()==False]
plate_mass = plate_mass.reset_index(drop=True)

#calculate the days the plate was in the field, convert to a decimal day
DaysDeploy = []
DaysDeployDecimal = []
remove_null = plate_mass.DateRemove.isnull()
for i in range(0, len(plate_mass)):
    if remove_null[i] == True:
        DaysDeploy.append(plate_mass.DateCollect[i]-plate_mass.DateDeploy[i])
    else:
        DaysDeploy.append(plate_mass.DateCollect[i]-plate_mass.DateDeploy[i] - (plate_mass.DateReplace[i]-plate_mass.DateRemove[i]))
    DaysDeployDecimal.append(DaysDeploy[i].days+DaysDeploy[i].seconds/86400.0)
plate_mass['DaysDeploy'] = DaysDeployDecimal

#calcluate the growth and growth rate in mg/day
plate_mass['Growth'] = (plate_mass.FinalWeight-plate_mass.StandardWeight_Final)-(plate_mass.InitialWeight-plate_mass.StandardWeight_Initial)
plate_mass['GrowthRate'] = (plate_mass.Growth/plate_mass.DaysDeploy)*1000
plate_mass['Midpoint'] = ((plate_mass.DateCollect-plate_mass.DateDeploy)/2)+plate_mass.DateDeploy
plate_mass = plate_mass[((plate_mass.Site == 'Flatman') == True) | ((plate_mass.Site == 'Station 1') == True) | ((plate_mass.Site == 'Station 2') == True) | ((plate_mass.Site == 'Stumpy') == True)]

#assemble table with mass and isotope data
calcite_data = pd.merge(plate_mass, plate_isotope, how = 'outer', left_on = 'SampleName', right_on = 'PlateName')
calcite_data = calcite_data.drop('Site_y', 1)
calcite_data.rename(columns={'Site_x': 'Site'}, inplace=True)

#format dripwater input tables
dripwater['BeginTrip'] = pd.to_datetime(dripwater['BeginTrip'])
dripwater['EndTrip'] = pd.to_datetime(dripwater['EndTrip'])
dripwater['DateAnalyzed'] = pd.to_datetime(dripwater['DateAnalyzed'])
dripwater = dripwater.reset_index(drop = False).sort_values(['Trip'])
dripwater.rename(columns={'index': 'index_record'}, inplace=True)
dripwater_full = dripwater.copy()

#find and remove all samples run in Breecker's lab that are replicates of samples run in Quinn's lab
last_quinn_index = len(dripwater)-1-int(list(dripwater['Lab'][::-1]).index('Quinn'))
last_quinn_trip = dripwater.Trip.iloc[last_quinn_index]
breecker_keep = dripwater[((dripwater['Trip'] > last_quinn_trip) & (dripwater['Lab'] == 'Breecker'))]
quinn_all = dripwater[dripwater['Lab'] == 'Quinn']
dripwater = breecker_keep.append(quinn_all)
dripwater = dripwater.sort_values(['Site', 'Trip'])

#find duplicates and multiples in the dripwater table and split them up by type
all_grouped = dripwater.groupby(['Site','Trip'])
processed = DataFrame()
singlets = DataFrame()
average_duplicates = DataFrame()
screen_duplicates = DataFrame()
screen_multiples = DataFrame()
left_overs = DataFrame()

tol = 0.2
    
for name, group in all_grouped:
    if group['d18O_Water'].count() == 1:
        singlets = singlets.append(group)
    elif group['d18O_Water'].count() > 1 and group['d18O_Water'].std() < tol:
        average_duplicates = average_duplicates.append(group)
    elif group['d18O_Water'].count() == 2 and group['d18O_Water'].std() > tol:
        screen_duplicates = screen_duplicates.append(group)
    elif group['d18O_Water'].count() > 2 and group['d18O_Water'].std() > tol:
        screen_multiples = screen_multiples.append(group)
    else:
        left_overs = left_overs.append(group)

#processed is a DataFrame that holds all isotope measurements that have been dealt with and how they were handled
processed = pd.concat([singlets.Site, singlets.Trip, singlets.d18O_Water,  singlets.d18O_Water_Error, singlets.index_record], axis = 1)

#makes a tuple from a DataFrame grouped by Site and Trip, this is probably totally unncessary
def make_tuple(c):
    try:
       dupes = [c.iloc[0]]
       curr_trip = dupes[0].Trip
       curr_site = dupes[0].Site
    
       ret = []
       for index, row in c[1:].iterrows():
           if row.Trip == curr_trip and row.Site == curr_site:
               dupes.append(row)
           else:
               ret.append(dupes)
               dupes = [row]
               curr_trip = row.Trip
               curr_site = row.Site
               ret.append(dupes)
    except IndexError:
        ret = []
    finally:
        return ret

#averages duplicate measurements if they are within tolerance, takes a tuple made by the make_tuple function
def average_duplicates_function(series):
    results = DataFrame(columns = ['Site','Trip','d18O_Water','d18O_Water_Error','index_record'])
    for i in range(0,len(series)):
        unpacked_series = pd.DataFrame(list(series[i]))
        dicts = {'Site' : unpacked_series.Site.iloc[0], 'Trip': unpacked_series.Trip.iloc[0], 'd18O_Water': np.mean(unpacked_series.d18O_Water), 'd18O_Water_Error':  np.std(unpacked_series.d18O_Water), 'index_record' : 'average of %s' %(list(unpacked_series.index_record))}         
        results = results.append([dicts])
    results = results.reset_index(drop = True)
    return results
    
tuple_average_duplicates = make_tuple(average_duplicates)
average_duplicates_processed = average_duplicates_function(tuple_average_duplicates)
processed = processed.append(average_duplicates_processed).sort_values(['Site','Trip']).reset_index(drop = True)

#in the case of duplicates that are not within tolerance, compares samples with nearest neighbors
def compare_neighbors(series, dripwater):
    site = str(series.Site.iloc[0])
    trip = int(series.Trip.iloc[0])
    index_site_first = list(dripwater['Site']).index(site)
    index_site_last_dummy = int(list(dripwater['Site'][::-1]).index(site))
    index_site_last = len(dripwater)-1-index_site_last_dummy
    index_sitetrip_first = list(dripwater['Trip']).index(trip, index_site_first)
    index_sitetrip_last_dummy = int(list(dripwater['Trip'][::-1]).index(trip, index_site_last_dummy))
    index_sitetrip_last = len(dripwater)-1-index_sitetrip_last_dummy
    if dripwater.Site.iloc[index_sitetrip_first-1] == site:
        prev_trip = dripwater.Trip.iloc[index_sitetrip_first-1]
    else:
        prev_trip = -1
    if dripwater.Site.iloc[index_sitetrip_last+1] == site:
        next_trip = dripwater.Trip.iloc[index_sitetrip_last+1]
    else:
        next_trip = -1
    return prev_trip, next_trip, site

tuple_screen_duplicates = make_tuple(screen_duplicates)

#processing duplicates that are not within tolerance
screen_duplicates_results = DataFrame(columns = ['Site','Trip','d18O_Water','d18O_Water_Error','index_record'])
for i in range(0,len(tuple_screen_duplicates)):
    unpacked_series = pd.DataFrame(list(tuple_screen_duplicates[i]))
    dicts = {'Site' : unpacked_series.Site.iloc[0], 'Trip':unpacked_series.Trip.iloc[0]}
    prev_trip, next_trip, site = compare_neighbors(unpacked_series, dripwater)
    if ((processed['Site'] == site) & (processed['Trip'] == prev_trip)).any() and ((processed['Site'] == site) & (processed['Trip'] == next_trip)).any() == True:
        prev_trip_value = processed[(processed['Site'] == site) & (processed['Trip'] == prev_trip)]
        next_trip_value = processed[(processed['Site'] == site) & (processed['Trip'] == next_trip)]
        compare_to = np.mean([prev_trip_value.d18O_Water.iloc[0], next_trip_value.d18O_Water.iloc[0]])
        if np.abs(unpacked_series.d18O_Water.iloc[0] - compare_to) < np.abs(unpacked_series.d18O_Water.iloc[1] - compare_to):
           dicts['d18O_Water'] = float(unpacked_series.d18O_Water.iloc[0])
           dicts['d18O_Water_Error'] = float(unpacked_series.d18O_Water_Error.iloc[0])
           dicts['index_record'] = 'selected %s over %s by comparison with nearest neighbors' %(unpacked_series.index_record.iloc[0], unpacked_series.index_record.iloc[1])
        else: 
           dicts['d18O_Water'] = float(unpacked_series.d18O_Water.iloc[1])
           dicts['d18O_Water_Error'] = float(unpacked_series.d18O_Water_Error.iloc[1])
           dicts['index_record'] = 'selected %s over %s by comparison with nearest neighbors' %(unpacked_series.index_record.iloc[1], unpacked_series.index_record.iloc[0])
    else:
        #need to revisit this shit eventually, but right now this case doesn't occur
        print('doot mr. skeletal, this is the end')
        dicts['d18O_Water'] = 0.0
        dicts['d18O_Water_Error'] = 0.0
        dicts['index_record'] = 'these werent processed'
    screen_duplicates_results = screen_duplicates_results.append([dicts])
screen_duplicates_processed = screen_duplicates_results.reset_index(drop = True)
processed = processed.append(screen_duplicates_processed).sort_values(['Site','Trip']).reset_index(drop = True)

tuple_screen_multiples = make_tuple(screen_multiples)

#this isn't totally working, the popping will take two out.  Should probably try popping both and then comarping to see which one does better.  
screen_multiples_results = DataFrame(columns = ['Site','Trip','d18O_Water','d18O_Water_Error','index_record'])
for i in range(0,len(tuple_screen_multiples)):
    unpacked_series = pd.DataFrame(list(tuple_screen_multiples[i])).sort_values(['d18O_Water'])
    dicts = {'Site' : unpacked_series.Site.iloc[0], 'Trip':unpacked_series.Trip.iloc[0]}
    to_exclude = []
    to_keep = []
    if len(unpacked_series) < 4:
        test = [0,-1]
        for j in range(0,len(test)):
            d18O_Water_list = list(unpacked_series['d18O_Water'])
            popped = d18O_Water_list.pop(test[j])
            if  np.abs(popped - np.mean(d18O_Water_list)) > tol and np.std(d18O_Water_list) < tol:
                to_exclude.append(unpacked_series['index_record'].iloc[test[j]])
    exclude_records = unpacked_series.loc[to_exclude]
    keep_records = unpacked_series.loc[list(set(unpacked_series.index_record) - set(to_exclude))]
    dicts['d18O_Water'] = np.mean(keep_records.d18O_Water)
    dicts['d18O_Water_Error'] = np.std(keep_records.d18O_Water)
    dicts['index_record'] = 'excluded %s from set of %s' %(list(exclude_records['index_record']), list(unpacked_series['index_record']))
    screen_multiples_results = screen_multiples_results.append([dicts])


screen_multiples_results = screen_multiples_results.reset_index(drop = True)
processed = processed.append(screen_multiples_results).sort_values(['Site','Trip']).reset_index(drop = True)

trip.rename(columns={'idTrip': 'Trip'}, inplace=True)

processed['Trip'] = processed.Trip.astype('int64')


processed = pd.merge(processed, trip, how = 'left', on = 'Trip')

water_data = pd.merge(processed, ultrameter, how = 'left',on = ['Site','Trip'])

pool_data = ultrameter[ultrameter.Site == 'Pool']
ultrameter_drips = ultrameter.loc[list(set(ultrameter.idUltrameterData) - set(pool_data.idUltrameterData))]
ultrameter_drips = pd.merge(ultrameter_drips, trip, how = 'left', on = 'Trip')

dripwater_temp_mean = np.mean(ultrameter_drips.WaterTemp)
dripwater_temp_max = np.max(ultrameter_drips.WaterTemp)
dripwater_temp_min = np.min(ultrameter_drips.WaterTemp)
dripwater_temp_mean_error = np.std(ultrameter_drips.WaterTemp)


dripaverage = []
for i in range(0,len(dripinterval)):
    dripaverage.append(np.mean([dripinterval.Count1[i], dripinterval.Count2[i], dripinterval.Count3[i]]))

dripinterval['DripIntervalAverage'] = dripaverage
    
dripinterval = pd.concat([dripinterval.Site, dripinterval.FieldTrip, dripinterval.DripIntervalAverage], axis = 1)

dripinterval = dripinterval.groupby(['Site', 'FieldTrip']).mean().reset_index()

dripinterval.rename(columns={'FieldTrip': 'idFieldTrip'}, inplace=True)
water_data.rename(columns={'Trip': 'idFieldTrip'}, inplace=True)

water_data = pd.merge(water_data,dripinterval, how = 'left', on = ('Site', 'idFieldTrip'))
water_data.rename(columns={'idFieldTrip': 'Trip'}, inplace=True)


#need to modify this funciton so that it's using the water temp, 
#this means I also need to find a way to deal with water temps 
#since they're at a different interval than the plates and they're not continous

def average_corr_water(calcite,water):
    corr_trip_rec = []
    water_rec = []
    d18O_water_mean_rec = []
    temp_water_mean_rec = []
    pH_water_mean_rec = []
    dripinterval_mean_rec = []
    for i in range(0, len(calcite)):
        test_site = calcite.Site[i]
        corr_trip = list(range(int(calcite.DeployFieldTrip.iloc[i]), int(calcite.CollectFieldTrip.iloc[i])))
        corr_trip.append(int(calcite.CollectFieldTrip.iloc[i]))
        corr_trip_rec.append(list(corr_trip))
        water_values = water[water['Trip'].isin(corr_trip) & ((water.Site == test_site) == True)]
        water_rec.append(list(water_values.Trip))
        corr_water_d18O = float(np.mean(water_values.d18O_Water))
        d18O_water_mean_rec.append(corr_water_d18O)
        corr_water_temp = float(np.mean(water_values.WaterTemp))
        temp_water_mean_rec.append(corr_water_temp)
        corr_water_pH = float(np.mean(water_values.pH))
        pH_water_mean_rec.append(corr_water_pH)
        corr_water_interval = float(np.mean(water_values.DripIntervalAverage))
        dripinterval_mean_rec.append(corr_water_interval)   
    calcite['d18O_Water'] = d18O_water_mean_rec
    calcite['Trips Averaged'] = water_rec
    calcite['WaterTemp'] = temp_water_mean_rec
    calcite['pH'] = pH_water_mean_rec
    calcite['DripInterval'] = dripinterval_mean_rec
    return calcite
    
def lnalpha(series):
    results = []
    for i in range(0,len(series)):
        if (isnan(series.d18OCalcite.iloc[i]) == False) & (isnan(series.d18OWater.iloc[i]) == False):
            result = 1000*(log((series.d18OCalcite.iloc[i] + 1000)/(series.d18OWater.iloc[i] + 1000)))
            results.append(result)
        else:
            results.append(float('Nan'))
    series['1000lnalpha'] = results
    return series
    
def fill_temp(series,ultrameter):
    results = []
    corr_trip_rec = []
    trinity_temp = ultrameter[ultrameter.Site == 'Trinity']
    for i in range(0,len(series)):
        if (isnan(series.WaterTemp.iloc[i]) == True):
            corr_trip = list(range(int(series.DeployFieldTrip.iloc[i]), int(series.CollectFieldTrip.iloc[i])))
            corr_trip.append(int(series.CollectFieldTrip.iloc[i]))
            corr_trip_rec.append('Trinity')
            corr_trip_rec.append(list(corr_trip))
            corr_trinity_temp = trinity_temp[trinity_temp['Trip'].isin(corr_trip)]
            results.append(list(corr_trinity_temp.Trip))
            series.WaterTemp.iloc[i] = float(np.mean(corr_trinity_temp.WaterTemp))
        else:
            results.append(float('Nan'))
    series['WaterTempFlag'] = results
    return series
    
fill_water_temp = average_corr_water(calcite_data, water_data)
fill_water_temp = fill_temp(fill_water_temp,ultrameter)
external_alpha = lnalpha(external)
fill_water_temp.rename(columns={'d18O_Water': 'd18OWater', 'd18O_Calcite':'d18OCalcite'}, inplace=True)
internal_alpha = lnalpha(fill_water_temp)

cutoff = datetime.datetime(2011,3,1)
filtered_internal = internal_alpha[internal_alpha['Midpoint'] > cutoff].copy()
plt.plot(1000/(internal_alpha.WaterTemp+273),internal_alpha['1000lnalpha'], 'wo', label = 'This publication')

external_records = list(set(external_alpha.PublicationKey))

marker_list = ['o','v','s','D']

color_palette = [(240,163,255),(0,117,220),(153,63,0),(76,0,92),(25,25,25),(0,92,49),(43,206,72),(255,204,153),(128,128,128),(148,255,181),(143,124,0),(157,204,0),(194,0,136),(0,51,128),(255,164,5),(255,168,187),(66,102,0),(255,0,16),(94,241,242),(0,153,143),(224,255,102),(116,10,255),(153,0,0),(255,255,128),(255,80,5)]    

for i in range(len(color_palette)):    
    r, g, b = color_palette[i]    
    color_palette[i] = (r / 255., g / 255., b / 255.)    

for i in range(0,len(external_records)):
    indv_record = external_alpha[external_alpha.PublicationKey == external_records[i]]
    plt.plot(1000/(indv_record.WaterTemp+273),indv_record['1000lnalpha'], linestyle='None', color = color_palette[random.randint(0,len(color_palette)-1)], marker = marker_list[random.randint(0,3)], label = external_records[i])

T = [1000.0/(273-10), 1000.0/(273+75)]
Coplen_alpha = [17.4*(i)-28.6 for i in T]
#Tremaine_alpha = [16.1*(i)-24.6 for i in T]
Affek_alpha = [15.63*(i)-23.29 for i in T]
KimONeil_alpha = [18.04*(i)-32.42 for i in T]
plt.plot(T, Coplen_alpha, 'k--', label = 'Coplen Line')
plt.plot(T, Affek_alpha, 'k-',label = 'Affek and Zaarur Line')
plt.plot(T, KimONeil_alpha, 'k:', label = 'Kim and ONeil Line')


plt.ylabel('1000ln(alpha)')
plt.xlabel('1000/T (K)')
plt.legend(loc ='lower right', numpoints = 1, ncol =4)



#internal_alpha['Month'] = [i.month for i in internal_alpha.Midpoint]
#cutoff = datetime.date(2011,3,1)
#filtered_internal = internal_alpha[internal_alpha.Midpoint > cutoff]
#stumpy_results = filtered_internal[filtered_internal.Site == 'Stumpy']
#station1_results = filtered_internal[filtered_internal.Site == 'Station 1']
#station2_results = filtered_internal[filtered_internal.Site == 'Station 2']
#flatman_results = filtered_internal[filtered_internal.Site == 'Flatman']
#plot(stumpy_results.Month,stumpy_results['1000lnalpha'], 'ro', label = 'Stumpy')
#plot(flatman_results.Month,flatman_results['1000lnalpha'], 'go', label = 'Flatman')
#plot(station1_results.Month ,station1_results['1000lnalpha'], 'mo', label = 'Station 1')
#plot(station2_results.Month ,station2_results['1000lnalpha'], 'yo', label = 'Station 2')
#plot(1000/(true_temp.WaterTemp+273), true_temp['1000lnalpha'], 'go')



def plot_cross(series, site, symbol, param1, param2):
    if (site == 'all') == False:
        series = series[series.Site == site]
    to_regress = series[(series[param1].isnull()==False) & (series[param2].isnull() == False)]
    (m,b) = np.polyfit(to_regress[param1], to_regress[param2], 1)
    yp = polyval([m,b], to_regress[param1])
    coeff = np.corrcoef(to_regress[param1], to_regress[param2])[0, 1]
    plot(series[param1], series[param2], symbol)
    plot(to_regress[param1], yp, 'r-')
    title(site)
    xlabel(param1)
    ylabel(param2)
    #x_cord = ((abs(max(to_regress[param1]))-abs(min(to_regress[param1])))/2)+min(to_regress[param1])
    x_cord = min(to_regress[param1])
    y_cord = max(to_regress[param2])
    #y_cord = ((abs(max(to_regress[param2]))-abs(min(to_regress[param2])))/2)+min(to_regress[param2])
    text(x_cord,y_cord,'y = %5.2fx+%5.2f, r^2 = %5.2f' %(m, b, coeff))

both_iso = pd.read_csv('plateisotope2.csv')
both_iso = pd.merge(both_iso,plate_name, how = 'left', left_on = 'SampleName', right_on = 'PlateName')
both_iso = pd.merge(both_iso,plate_mass, how = 'left', left_on = 'SampleName', right_on = 'SampleName')
both_iso.rename(columns={'Site_x': 'Site'}, inplace=True)


cutoff = datetime.datetime(2012,1,1)

both_iso = both_iso[both_iso.Midpoint > cutoff]

stumpy_internal = internal_alpha[internal_alpha.Site == 'Stumpy']
#stumpy_internal['TremaineEquilb'] = [exp(16.1/(i+273.0)-0.0246)) for i in stumpy_internal.WaterTemp]
stumpy_internal['TremaineEquilb']  = (e**(16.1/(stumpy_internal.WaterTemp+273.15)-0.0246))*(1000+stumpy_internal.d18OWater)-1000
stumpy_internal['TremaineDeparture'] = stumpy_internal.d18OCalcite-stumpy_internal.TremaineEquilb
stumpy_internal = stumpy_internal.sort_values('Midpoint').copy()

cdioxide = pd.read_csv('caveco2test.csv', names = ['idFieldCO2','Site','FieldInstrumentName','LogNumber','Datetime','CO2'])
cdioxide['Datetime'] = pd.to_datetime(cdioxide['Datetime'])
daily_cdioxide = cdioxide.set_index(pd.DatetimeIndex(cdioxide['Datetime']))
daily_cdioxide = daily_cdioxide.resample('D').mean()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(daily_cdioxide.index.to_pydatetime(), daily_cdioxide.CO2,'ko-',markersize = 4)
ax.set_ylim(0,5000)
ax.set_ylabel('pCO2')
ax2 = ax.twinx()
ax2.plot(stumpy_internal.Midpoint, stumpy_internal.TremaineDeparture, 'ro-')
ax2.set_ylabel ('Departure from Equilibrium', color ='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')
plt.xlabel('Date')
plt.title('Jinapsan Cave Daily Average CO2 and Stumpy Departure from Equilibrium')



#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.plot(station2.Midpoint, station2.d18O_Calcite, 'ko-', label ='d18O Calcite')
#ax.set_ylabel ('d18O Calcite',color ='k')
#ax2 = ax.twinx()
#ax2.plot(station2.Midpoint, station2.d18O_Water, 'ro-',label ='d18O Water')
#ax2.set_ylabel ('d18O Water', color ='r')
#for tl in ax2.get_yticklabels():
#    tl.set_color('r')


#def plot_cross(series, symbol):
#    series = series[series.GrowthRate > 0]
#    series = series[np.abs(series.GrowthRate) < np.mean(series.GrowthRate)+(np.std(series.GrowthRate))*3]
#    to_regress = series[(series.GrowthRate.isnull()==False) & (series['Coplen departure'].isnull() == False)]
#    (m,b) = np.polyfit(to_regress.GrowthRate, to_regress['Coplen departure'], 1)
#    yp = polyval([m,b], to_regress.GrowthRate)
#    plot(series.GrowthRate, series['Coplen departure'], symbol)
#    plot(to_regress.GrowthRate, yp, 'r-')
#    title(series.Site[0])
#    xlabel('Growth Rate (mg/day)')
#    ylabel('Coplen Departure')
#    x_cord = (max(to_regress.GrowthRate)-min(to_regress.GrowthRate))/2
#    y_cord = (max(to_regress['Coplen departure'])-min(to_regress['Coplen departure']))/2
#    text(x_cord,y_cord,'y = %5.2fx+%5.2f' %(m, b))
##
#def plot_diff(series, symbol):
#    plot(series.Midpoint, series['Coplen departure'], symbol)
#    title(series.Site[0])
#    xlabel('Date')
#    ylabel('Coplen Departure')        

#f, axarr = plt.subplots(2,2)
#axarr[0, 0].plot(station1_calcite_isotope['Midpoint'],station1_calcite_isotope['Tremaine departure'],'ko-')
#axarr[0, 0].set_title('Station 1')
#axarr[0, 0].plot(station1_water_isotope['EndTrip'],station1_water_isotope['Temp'],'ro')
#
#ax1 = plt.subplot(412)
#plt.title('Station 1')
#ax1.plot(station1_calcite_isotope['Midpoint'],station1_calcite_isotope['Tremaine departure'],'ko-')
##ax1.set_xlabel('Date')
#ax1.set_ylabel('d18Occ-e', color='k')
#for tl in ax1.get_yticklabels():
#    tl.set_color('k')
#
#ax2 = ax1.twinx()
#ax2.plot(station1_water_isotope['EndTrip'],station1_water_isotope['Temp'],'ro')
#ax2.set_ylabel('Water Temp', color='r')
#for tl in ax2.get_yticklabels():
#    tl.set_color('r')
#plt.show()
#
#ax3 = plt.subplot(422)
#plt.title('Flatman')
#ax3.plot(flatman_calcite_isotope['Midpoint'],flatman_calcite_isotope['Tremaine departure'],'ko-')
#ax3.set_xlabel('Date')
## Make the y-axis label and tick labels match the line color.
#ax3.set_ylabel('d18Occ-e', color='k')
#for tl in ax3.get_yticklabels():
#    tl.set_color('k')
#
#ax4 = ax3.twinx()
#ax4.plot(flatman_water_isotope['EndTrip'],flatman_water_isotope['Temp'],'ro')
#ax4.set_ylabel('Water Temp', color='r')
#for tl in ax4.get_yticklabels():
#    tl.set_color('r')
#plt.show()
#
#
#ax5 = plt.subplots(421)
#plt.title('Station 2')
#ax5.plot(station2_calcite_isotope['Midpoint'],station2_calcite_isotope['Tremaine departure'],'ko-')
#ax5.set_xlabel('Date')
## Make the y-axis label and tick labels match the line color.
#ax5.set_ylabel('d18Occ-e', color='k')
#for tl in ax5.get_yticklabels():
#    tl.set_color('k')
#
#ax6 = ax5.twinx()
#ax6.plot(station2_water_isotope['EndTrip'],station2_water_isotope['Temp'],'ro')
#ax6.set_ylabel('Water Temp', color='r')
#for tl in ax6.get_yticklabels():
#    tl.set_color('r')
#plt.show()
#
#
#ax7 = plt.subplots(422)
#plt.title('Stumpy')
#ax7.plot(stumpy_calcite_isotope['Midpoint'],stumpy_calcite_isotope['Tremaine departure'],'ko-')
#ax7.set_xlabel('Date')
## Make the y-axis label and tick labels match the line color.
#ax7.set_ylabel('d18Occ-e', color='k')
#for tl in ax7.get_yticklabels():
#    tl.set_color('k')
#
#ax8 = ax7.twinx()
#ax8.plot(stumpy_water_isotope['EndTrip'],stumpy_water_isotope['Temp'],'ro-')
#ax8.set_ylabel('Water Temp', color='r')
#for tl in ax8.get_yticklabels():
#    tl.set_color('r')
#plt.show()
#
