# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 13:10:44 2016

@author: Peter Carlson
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime
import matplotlib.dates as mdates

import oxygen_isotope_stats_functions as oxystats

calcite_data = oxystats.load_cave_data()
drip_data = oxystats.load_drip_data()
external_data = oxystats.load_external_data()
FieldCO2Continous = oxystats.load_continuous_co2()
FieldCO2Spot = oxystats.load_spot_co2()


cap18O = '${\Delta}^{18}O$'
d18O = '${\delta}^{18}O$'
d13C = '${\delta}^{13}C$'
permil = '$\u2030$'
# %% Initialize

start_plot = datetime.datetime(2008, 6, 1)
co2_start_plot = datetime.datetime(2012, 1, 1)
end_plot = datetime.datetime(2016, 8, 1)
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

plot_table = pd.DataFrame(columns = ['DripName','Color','Marker'])
plot_table.DripName = DripNames
plot_table.Color = rgb_convert([
    (117,112,179),  #orange
    (27,158,119), #green?
    (217,95,2), #pink
    (231,41,138), #purple
    ])
plot_table.Marker = ['o','D','v','^']
plot_table = plot_table.sort_values('DripName').reset_index(drop=True, 
                                                            inplace=False)


# %% Figure 1

## Drip Rate timeseries
print(drip_data.columns)
fig = plt.figure(figsize=(10,5.63))
ax = fig.add_subplot(111)      
for i, site in plot_table.iterrows():
    plotname = site['DripName']
    plotmarker = site['Marker']
    plotcolor = site['Color']
    set_to_plot = drip_data[drip_data['SiteName'] == plotname]
    set_to_plot = set_to_plot.sort_values('Date')
    set_to_plot['mdate'] =mdates.date2num(set_to_plot['Date'])
    est_drips_per_min = (set_to_plot['BottleFillRate']
                         / set_to_plot['DripVolume'].mean())
    set_to_plot['DripsPerMin'] = set_to_plot['DripsPerMin'].mask(
        set_to_plot['DripsPerMin'].isnull(),
        est_drips_per_min)
    ax.plot(set_to_plot['mdate'], 
            set_to_plot['DripsPerMin'],
            linestyle='-',
            color=plotcolor,
            marker=plotmarker,
            markeredgewidth=1,
            linewidth=1,
            markeredgecolor=plotcolor,
            label=plotname,
            markersize=5)
    
    ax.set_xlim([start_plot, end_plot])
    ax.tick_params(axis='both', which='major', labelsize=15)
for label in ax.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)
ax.set_ylim([0,7])
ax.legend(loc='upper center', 
          bbox_to_anchor=(0.5, 1.1), 
          frameon=False, 
          fancybox=False, 
          shadow=False, 
          ncol=4,
          numpoints=1,
          prop={'size':15})
tmp_ticks_labels = [10, 20, 30, 60, 120, 300]
tmp_ticks_pos = [60.0/x for x in tmp_ticks_labels]
ax2 = ax.twinx()
ax2.set_ylim([0,7])
ax2.set_yticks(tmp_ticks_pos)
ax2.set_yticklabels(tmp_ticks_labels)
ax2.tick_params(axis='y', labelsize=12)
ax2.set_ylabel(u'Drip interval\n (seconds)', fontsize=15)
ax.set_ylabel(u'Drip rate \n(drips/min)', fontsize=15)
plt.tight_layout()
plt.show()
plt.savefig('Carlson2020resubmit_fig01.png', format='png', dpi=600)

# %% Figure 2
# Drip and calcite d18O
plot2_min = -10
plot2_max = -3
fig, axes = plt.subplots(nrows=2,ncols=4,sharex=True, figsize=(10, 5.63))
for i, ax in enumerate(axes.flatten()):
    if i < 4:
        plotname = plot_table.DripName.iloc[i]
        plotmarker = plot_table.Marker.iloc[i]
        plotcolor = plot_table.Color.iloc[i]
        set_to_plot = drip_data[drip_data.SiteName == plotname]
        set_to_plot = set_to_plot.sort_values('Date')
        ax.errorbar(x=set_to_plot['Date'], 
                    y=set_to_plot['δ18OWaterVSMOW'],
                    yerr=set_to_plot['δ18OWaterErr'],
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
        ax.set_title(plotname)
        plot_ylim = ax.get_ylim()
        ax.set_ylim([plot2_min, plot2_max])
        ax.set_xlim([start_plot, end_plot])
        for label in ax.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        if i > 0 and i < 4:
            for label in ax.yaxis.get_ticklabels():
                label.set_visible(False)   
        if i == 0:
            ax.set_ylabel(f'{d18O} Water \n ({permil} VSMOW)', fontsize=15)
    else:
        i = i - 4
        plotname = plot_table.DripName.iloc[i]
        plotmarker = plot_table.Marker.iloc[i]
        plotcolor = plot_table.Color.iloc[i]
        set_to_plot = calcite_data[calcite_data.SiteName == plotname]
        set_to_plot = set_to_plot.sort_values('MidPoint')
        ax.errorbar(x=set_to_plot['MidPoint'], 
                    y=set_to_plot['Calcite_δ18O_VPDB'],
                    yerr=set_to_plot['Calcite_δ18O_Err_VPDB'],
                    linestyle = '-',
                    color = plotcolor,
                    marker = plotmarker,
                    markeredgewidth = 1,
                    linewidth = 1,
                    markeredgecolor = plotcolor,
                    label = plotname,
                    markersize = 5)
        plot_ylim = ax.get_ylim()
    #    if plot_ylim[0] >=0:
    #        plot_ylim = (plot_ylim[0] - .1*plot_ylim[1], plot_ylim[1])
        ax.set_ylim([plot2_min, plot2_max])
        ax.set_xlim([start_plot, end_plot])
    #ax0.set_xlim([start_plot, end_plot])
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
        for label in ax.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        if i > 0 and i < 4:
            for label in ax.yaxis.get_ticklabels():
                label.set_visible(False)    
        if i == 0:
            ax.set_ylabel(f'{d18O} Calcite \n ({permil} VPDB)', fontsize=15)
            plt.subplots_adjust(wspace=0.1, hspace=0.1, bottom=0.1, top=0.85)
plt.tight_layout()
plt.show()
plt.savefig('Carlson2020resubmit_fig02.png', format='png', dpi=600)
# %% Figure 3           
# Cap Delta vs. Time

fig8, axes = plt.subplots(2,2, figsize=(10,5.5))
for i, site in plot_table.iterrows():
    plotname = site['DripName']
    plotmarker = site['Marker']
    plotcolor = site['Color']
    ax = axes[i//2,i%2]
    set_to_plot = calcite_data[calcite_data.SiteName == plotname]
    set_to_plot = set_to_plot.sort_values('MidPoint')
    set_to_plot['mdate'] =mdates.date2num(set_to_plot['MidPoint'])
    ax.set_title(f'{plotname}')
    ax.errorbar(x=set_to_plot['mdate'],
                y=set_to_plot['Δ18OPDB'],
                yerr=set_to_plot['Δ18OPDB_Err'],
                linestyle='-',
                color=plotcolor,
                marker=plotmarker,
                markersize=4,
                linewidth=1,
                label='Measured')
    equil_names = ['Coplen', 'Affek', 'Kim']
    
    ax.legend(fontsize=8)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y')) 
    ax.set_xlim([start_plot, end_plot])
    ax.set_ylim([26, 31])
    ax.tick_params(axis='both', which='major', labelsize=10)
    
    plt.setp(ax.get_xticklabels(), rotation=0, ha='center')
    plt.savefig('PCO2.png', format='png', dpi=600)
    if i%2 != 0:
        for label in ax.yaxis.get_ticklabels():
            label.set_visible(False)
    else:
        ax.set_ylabel(f'{cap18O} (PDB)')
    if i//2 == 0:
        for label in ax.xaxis.get_ticklabels():
            label.set_visible(False)
plt.tight_layout()
plt.show()
plt.savefig('Carlson2020resubmit_fig03.png', format='png', dpi=600)
# %% Figure 4
# Published 1000lnalpha values
fig, ax = plt.subplots(1, 1, figsize=(5.6,5.6))
external_data['inv_TK'] = 1000 / t_to_k(external_data['WaterTemp'].values, 'C')
calcite_data['inv_TK'] = 1000 / t_to_k(calcite_data['WaterTemp'].values, 'C')

ax.plot(external_data['inv_TK'],
        external_data['1000lnα'],
        linestyle='none',
        color='xkcd:gray',
        marker='o',
        markersize=3,
        label='Other Publications',
        )
ax.plot(calcite_data['inv_TK'],
        calcite_data['meas_1000lnα'],
        linestyle='none',
        color='xkcd:dark gray',
        marker='o',
        markersize=4,
        label='This Study',
        )
xlim = ax.get_xlim()
ax.plot(xlim,
        affek_1000lna(1000/np.array(xlim)),
        linestyle='-',
        color='k',
        marker=None,
        label='Affek and Zaarur Line',
        )
ax.plot(xlim,
        coplen_1000lna(1000/np.array(xlim)),
        linestyle='--',
        color='k',
        marker=None,
        label='Coplen Line',
        )
ax.plot(xlim,
        kim_1000lna(1000/np.array(xlim)),
        linestyle=':',
        color='k',
        marker=None,
        label="Kim and O'Neil Line",
        )
ax.plot(xlim,
        tremaine_1000lna(1000/np.array(xlim)),
        linestyle='-.',
        color='k',
        marker=None,
        label="Tremaine Line",
        )
ax.set_xlim(xlim)
ax.set_ylabel('$1000ln{\alpha}$')
ax.set_xlabel('$1000/T [K]$')
ax.legend()
ax2 = ax.twiny()
ax2.set_xlabel('$T [{\degree}C]$')
ax2.set_xlim(xlim)
ax2ticks = [1000/t_to_k(t, 'C')
            for t 
            in range(0,40,5)]

ax2.set_xticks(ax2ticks)
ax2.set_xticklabels([str(t) for t in range(0,40,5)])
plt.tight_layout()
plt.show()
plt.savefig('Carlson2020resubmit_fig04.png', format='png', dpi=600)
# %% Figure 5
# Cap18O vs drip rate and growth rate
fig, axes = plt.subplots(1, 3, figsize=(10,5.6))

to_regress = calcite_data[(calcite_data['Δ18OPDB'].notnull()
                          & calcite_data['DripsPerMin'].notnull())]
text_y = max(to_regress['Δ18OPDB'])
text_x = max(to_regress['DripsPerMin'])*0.5
ax1 = axes[0]
ax2 = axes[1]
ax3 = axes[2]

fig_5_ylim = [27,31]

for i, site in plot_table.iterrows():
    plotname = site['DripName']
    plotmarker = site['Marker']
    plotcolor = site['Color']
    set_to_plot = to_regress[to_regress['SiteName'] == plotname]
    ax1.errorbar(set_to_plot['DripsPerMin'],
             set_to_plot['Δ18OPDB'],
             xerr=set_to_plot['DripRate_Err'],
             yerr=set_to_plot['Δ18OPDB_Err'],
             linestyle = 'none',
             color = plotcolor,
             marker = plotmarker,
             markeredgewidth = 1,
             linewidth = 1,
             markeredgecolor = plotcolor,
             label = plotname,
             markersize = 5)
    ax2.errorbar(set_to_plot['DripsPerMin'],
             set_to_plot['Δ18OPDB'],
             xerr=set_to_plot['DripRate_Err'],
             yerr=set_to_plot['Δ18OPDB_Err'],
             linestyle = 'none',
             color = plotcolor,
             marker = plotmarker,
             markeredgewidth = 1,
             linewidth = 1,
             markeredgecolor = plotcolor,
             label = plotname,
             markersize = 5)  
    ax3.errorbar(set_to_plot['Log_Rc'], 
             set_to_plot['Δ18OPDB'],
             xerr=set_to_plot['Log_Rc_Err'],
             yerr=set_to_plot['Δ18OPDB_Err'],
             linestyle = 'none',
             color = plotcolor,               
             marker = plotmarker,
             markeredgewidth = 1,             
             linewidth = 1,
             markeredgecolor = plotcolor,             
             label = plotname,
             markersize = 5)    
x = to_regress['DripsPerMin'].values
y = to_regress['Δ18OPDB'].values
func = lambda x, a, b: a * x + b
band_corr_dict = oxystats.correlate_with_uncertainty_band(x,
                                                          y,
                                                          func,
                                                          alpha=0.1,
                                                          band_type='prediction')

m = band_corr_dict['popt'][0]
b = band_corr_dict['popt'][1]
m_s = band_corr_dict['sigmas'][0]
b_s = band_corr_dict['sigmas'][1]
f_1 = lambda x: m * x + b
regr_dict = oxystats.correlate_with_ar_correction(f_1(x), y)
print(regr_dict)
ax1.fill_between(band_corr_dict['x_hat'],
                 y1=band_corr_dict['uncertainty_band'][0],
                 y2=band_corr_dict['uncertainty_band'][1],
                 color='k',
                 alpha=0.2)
ax1.plot(band_corr_dict['x_hat'], 
         band_corr_dict['x_hat'] * m + b,
         color='k',
         linestyle='-')
ax2.fill_between(band_corr_dict['x_hat'],
                 y1=band_corr_dict['uncertainty_band'][0],
                 y2=band_corr_dict['uncertainty_band'][1],
                 color='k',
                 alpha=0.2)
ax2.plot(band_corr_dict['x_hat'], 
         band_corr_dict['x_hat'] * m + b,
         color='k',
         linestyle='-')
ax2.text(text_x+1,
         text_y+0.15,
         f'$r^{{{2}}} =$ {regr_dict["r_squared"]:5.2f}',
         fontsize=12,
         ha='center')
ax2.text(text_x+1,
         text_y-0.15,
         f'$p =$ {regr_dict["pval"]:5.1e}',
         fontsize=12,
         ha='center')
ax2.vlines(x=1.5,
           ymin=fig_5_ylim[0],
           ymax=fig_5_ylim[1],
           linestyle='--',
           color='k',
           )
for label in ax2.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)
ax1.tick_params(axis='both', which='major', labelsize=15)
ax1.set_xlabel(u'Drip rate  \n[drip/min]', fontsize=15, labelpad=15)
ax1.set_ylabel(f'{cap18O} (PDB)', fontsize=15)
ax1.set_ylim(fig_5_ylim)
ax1.set_xlim([0, 1.5])
tmp_ticks_labels = [60, 120, 300]
tmp_ticks_pos = [60.0/x for x in tmp_ticks_labels]

ax1a = ax1.twiny()
ax1a.set_xlim([0,1.5])
ax1a.set_xlabel(u'Drip interval\n [seconds[)]', fontsize=15)
ax1a.set_xticks(tmp_ticks_pos)
ax1a.set_xticklabels(tmp_ticks_labels)
ax1a.tick_params(axis='x', labelsize=12)
ax1a.set_xlabel(u'Drip interval\n [seconds]', fontsize=15)

tmp_ticks_labels = [10, 20, 30, 60, 300]
tmp_ticks_pos = [60.0/x for x in tmp_ticks_labels]
ax2.set_xlim([0,7])
ax2.set_ylim(fig_5_ylim)
ax2.tick_params(axis='both', which='major', labelsize=15)
ax2.set_xlabel(u'Drip rate  \n[drip/min]', fontsize=15, labelpad=15)

for label in ax2.yaxis.get_ticklabels():
    label.set_visible(False)
ax2.set_ylim(fig_5_ylim)
ax2a = ax2.twiny()
ax2a.set_xlim([0,7])
ax2a.set_xticks(tmp_ticks_pos)
ax2a.set_xticklabels(tmp_ticks_labels)
ax2a.tick_params(axis='x', labelsize=12)
ax2a.set_xlabel(u'Drip interval\n [seconds]', fontsize=15)

# Growth Rate Plot
to_regress = calcite_data[(calcite_data['Δ18OPDB'].notnull() 
                          & calcite_data['Log_Rc'].notnull())]
x = to_regress['Log_Rc'].values
y = to_regress['Δ18OPDB'].values
func = lambda x, a, b: a * x + b
band_corr_dict = oxystats.correlate_with_uncertainty_band(x,
                                                          y,
                                                          func,
                                                          alpha=0.1,
                                                          band_type='prediction')

m = band_corr_dict['popt'][0]
b = band_corr_dict['popt'][1]
m_s = band_corr_dict['sigmas'][0]
b_s = band_corr_dict['sigmas'][1]
f_1 = lambda x: m * x + b
regr_dict = oxystats.correlate_with_ar_correction(f_1(x), y)
print(regr_dict)
ax3.fill_between(band_corr_dict['x_hat'],
                 y1=band_corr_dict['uncertainty_band'][0],
                 y2=band_corr_dict['uncertainty_band'][1],
                 color='k',
                 alpha=0.2)
ax3.plot(band_corr_dict['x_hat'], 
         band_corr_dict['x_hat'] * m + b,
         color='k',
         linestyle='-')

text_y = max(to_regress['Δ18OPDB'])
text_x = max(to_regress['Log_Rc'])*0.5
ax3.text(text_x+1,
         text_y+0.15,
         f'$r^{{{2}}} =$ {regr_dict["r_squared"]:5.2f}',
         fontsize=12,
         ha='center')
ax3.text(text_x+1,
         text_y-0.15,
         f'$p =$ {regr_dict["pval"]:5.1e}',
         fontsize=12,
         ha='center')

for label in ax3.yaxis.get_ticklabels():
    label.set_visible(False)
ax3.set_ylim(fig_5_ylim)
ax3.tick_params(axis='both', which='major', labelsize=15)
ax3.set_xlabel(u'Calcite Deposition Rate\n$[log({\mu}mol/m^2/hr)]$', 
               fontsize=15, 
               labelpad=15)
ax3.legend(loc='lower center', bbox_to_anchor=(0.15,1,0.7,0.2), ncol=2)
plt.tight_layout()
plt.show()
plt.savefig('Carlson2020resubmit_fig05.png', format='png', dpi=600)

# %% Figure 5 Alt.
# Cap18O vs drip rate and growth rate
fig, axes = plt.subplots(1, 3, figsize=(10,5.6))

to_regress = calcite_data[(calcite_data['Δ18OPDB'].notnull()
                          & calcite_data['DripsPerMin'].notnull())]
text_y = max(to_regress['Δ18OPDB'])
text_x = max(to_regress['DripsPerMin'])*0.5
ax1 = axes[0]
ax2 = axes[1]
ax3 = axes[2]

fig_5_ylim = [27,31]

for i, site in plot_table.iterrows():
    plotname = site['DripName']
    plotmarker = site['Marker']
    plotcolor = site['Color']
    set_to_plot = to_regress[to_regress['SiteName'] == plotname]
    ax1.errorbar(set_to_plot['DripsPerMin'],
             set_to_plot['Δ18OPDB'],
             xerr=set_to_plot['DripRate_Err'],
             yerr=set_to_plot['Δ18OPDB_Err'],
             linestyle = 'none',
             color = plotcolor,
             marker = plotmarker,
             markeredgewidth = 1,
             linewidth = 1,
             markeredgecolor = plotcolor,
             label = plotname,
             markersize = 5)
    ax2.errorbar(set_to_plot['DripsPerMin'],
             set_to_plot['Δ18OPDB'],
             xerr=set_to_plot['DripRate_Err'],
             yerr=set_to_plot['Δ18OPDB_Err'],
             linestyle = 'none',
             color = plotcolor,
             marker = plotmarker,
             markeredgewidth = 1,
             linewidth = 1,
             markeredgecolor = plotcolor,
             label = plotname,
             markersize = 5)  
    ax3.errorbar(set_to_plot['Log_Rc'], 
             set_to_plot['Δ18OPDB'],
             xerr=set_to_plot['Log_Rc_Err'],
             yerr=set_to_plot['Δ18OPDB_Err'],
             linestyle = 'none',
             color = plotcolor,               
             marker = plotmarker,
             markeredgewidth = 1,             
             linewidth = 1,
             markeredgecolor = plotcolor,             
             label = plotname,
             markersize = 5)    
x = to_regress['DripsPerMin'].values
y = to_regress['Δ18OPDB'].values
func = lambda x, a, b: a * x + b
band_corr_dict = oxystats.correlate_vs_driprate(calcite_data)
yp = oxystats.piecewise_linear_continuous(x=band_corr_dict['x_hat'],
                                 x_b=band_corr_dict['x_b'],
                                 m_1=band_corr_dict['m_1'],
                                 b_1=band_corr_dict['b_1'],
                                 m_2=band_corr_dict['m_2'])
ax1.fill_between(band_corr_dict['x_hat'],
                 y1=band_corr_dict['uncertainty_band'][0],
                 y2=band_corr_dict['uncertainty_band'][1],
                 color='k',
                 alpha=0.2)
ax1.plot(band_corr_dict['x_hat'], 
         yp,
         color='k',
         linestyle='-')
ax2.fill_between(band_corr_dict['x_hat'],
                 y1=band_corr_dict['uncertainty_band'][0],
                 y2=band_corr_dict['uncertainty_band'][1],
                 color='k',
                 alpha=0.2)

print(band_corr_dict)
ax2.plot(band_corr_dict['x_hat'], 
         yp,
         color='k',
         linestyle='-')
ax2.text(text_x+1,
         text_y+0.15,
         f'$r^{{{2}}} =$ {regr_dict["r_squared"]:5.2f}',
         fontsize=12,
         ha='center')
ax2.text(text_x+1,
         text_y-0.15,
         f'$p =$ {regr_dict["pval"]:5.1e}',
         fontsize=12,
         ha='center')
ax2.vlines(x=1.5,
           ymin=fig_5_ylim[0],
           ymax=fig_5_ylim[1],
           linestyle='--',
           color='k',
           )
for label in ax2.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)
ax1.tick_params(axis='both', which='major', labelsize=15)
ax1.set_xlabel(u'Drip rate  \n[drip/min]', fontsize=15, labelpad=15)
ax1.set_ylabel(f'{cap18O} (PDB)', fontsize=15)
ax1.set_ylim(fig_5_ylim)
ax1.set_xlim([0, 1.5])
tmp_ticks_labels = [60, 120, 300]
tmp_ticks_pos = [60.0/x for x in tmp_ticks_labels]

ax1a = ax1.twiny()
ax1a.set_xlim([0,1.5])
ax1a.set_xlabel(u'Drip interval\n [seconds[)]', fontsize=15)
ax1a.set_xticks(tmp_ticks_pos)
ax1a.set_xticklabels(tmp_ticks_labels)
ax1a.tick_params(axis='x', labelsize=12)
ax1a.set_xlabel(u'Drip interval\n [seconds]', fontsize=15)

tmp_ticks_labels = [10, 20, 30, 60, 300]
tmp_ticks_pos = [60.0/x for x in tmp_ticks_labels]
ax2.set_xlim([0,7])
ax2.set_ylim(fig_5_ylim)
ax2.tick_params(axis='both', which='major', labelsize=15)
ax2.set_xlabel(u'Drip rate  \n[drip/min]', fontsize=15, labelpad=15)

for label in ax2.yaxis.get_ticklabels():
    label.set_visible(False)
ax2.set_ylim(fig_5_ylim)
ax2a = ax2.twiny()
ax2a.set_xlim([0,7])
ax2a.set_xticks(tmp_ticks_pos)
ax2a.set_xticklabels(tmp_ticks_labels)
ax2a.tick_params(axis='x', labelsize=12)
ax2a.set_xlabel(u'Drip interval\n [seconds]', fontsize=15)

# Growth Rate Plot
to_regress = calcite_data[(calcite_data['Δ18OPDB'].notnull() 
                          & calcite_data['Log_Rc'].notnull())]
x = to_regress['Log_Rc'].values
y = to_regress['Δ18OPDB'].values
func = lambda x, a, b: a * x + b
band_corr_dict = oxystats.correlate_with_uncertainty_band(x,
                                                          y,
                                                          func,
                                                          alpha=0.1,
                                                          band_type='confidence')

m = band_corr_dict['popt'][0]
b = band_corr_dict['popt'][1]
m_s = band_corr_dict['sigmas'][0]
b_s = band_corr_dict['sigmas'][1]
f_1 = lambda x: m * x + b
regr_dict = oxystats.correlate_with_ar_correction(f_1(x), y)
print(regr_dict)
ax3.fill_between(band_corr_dict['x_hat'],
                 y1=band_corr_dict['uncertainty_band'][0],
                 y2=band_corr_dict['uncertainty_band'][1],
                 color='k',
                 alpha=0.2)
ax3.plot(band_corr_dict['x_hat'], 
         band_corr_dict['x_hat'] * m + b,
         color='k',
         linestyle='-')

text_y = max(to_regress['Δ18OPDB'])
text_x = max(to_regress['Log_Rc'])*0.5
ax3.text(text_x+1,
         text_y+0.15,
         f'$r^{{{2}}} =$ {regr_dict["r_squared"]:5.2f}',
         fontsize=12,
         ha='center')
ax3.text(text_x+1,
         text_y-0.15,
         f'$p =$ {regr_dict["pval"]:5.1e}',
         fontsize=12,
         ha='center')

for label in ax3.yaxis.get_ticklabels():
    label.set_visible(False)
ax3.set_ylim(fig_5_ylim)
ax3.tick_params(axis='both', which='major', labelsize=15)
ax3.set_xlabel(u'Calcite Deposition Rate\n$[log({\mu}mol/m^2/hr)]$', 
               fontsize=15, 
               labelpad=15)
ax3.legend(loc='lower center', bbox_to_anchor=(0.15,1,0.7,0.2), ncol=2)
plt.tight_layout()
plt.show()
plt.savefig('Carlson2020resubmit_fig05a.png', format='png', dpi=600)

# %% Figure 6
# Model Comparison
fig6_tick = [-10,-9,-8,-7, -6, -5]
plot_lim = [-10,-5]
fig6, axes = plt.subplots(2,3, figsize=(10,5.5))

models = ['equilib_calcite_δ18O_kim', 
          'equilib_calcite_δ18O_affek',
          'equilib_calcite_δ18O_coplen',
          'd18OCalcite_Modeled_Kim',
          'd18OCalcite_Modeled_Affek',
          'd18OCalcite_Modeled_Coplen'
          ]   

for j, model in enumerate(models):  
    to_regress = calcite_data[(calcite_data['Calcite_δ18O_VPDB'].notnull()
                               & calcite_data[model].notnull())].copy()
    to_regress = to_regress[to_regress[model].notnull()]  
    x = to_regress['Calcite_δ18O_VPDB']
    y = to_regress[model]
    func = lambda x, a, b: a * x + b
    band_corr_dict = oxystats.correlate_with_uncertainty_band(x,
                                                              y,
                                                              func,
                                                              alpha=0.1,
                                                              band_type='confidence')
    m = band_corr_dict['popt'][0]
    b = band_corr_dict['popt'][1]
    m_s = band_corr_dict['sigmas'][0]
    b_s = band_corr_dict['sigmas'][1]
    f_1 = lambda x: m * x + b
    regr_dict = oxystats.correlate_with_ar_correction(f_1(x), y)
    print(model)
    print(band_corr_dict['popt'])
    print(band_corr_dict['sigmas'])
    print(regr_dict)
    ax = axes[j//3,j%3]
    ax.plot([-10,-5],[-10,-5],'k-')
    for i in range(0,len(plot_table)):
        plotname = plot_table.DripName.iloc[i]
        plotmarker = plot_table.Marker.iloc[i]
        plotcolor = plot_table.Color.iloc[i]   
        set_to_plot = to_regress[to_regress.SiteName == plotname]
        ax.plot(set_to_plot['Calcite_δ18O_VPDB'], 
                set_to_plot[model],
                linestyle='none',
                color=plotcolor,
                marker=plotmarker,
                markeredgewidth=1,
                linewidth=1,
                markeredgecolor=plotcolor,
                label=plotname,
                markersize=5)
    ax.plot([x.min(), x.max()],
            [np.min(f_1(x)), np.max(f_1(x))],
            color='xkcd:dark gray',
            linewidth=2,
            linestyle='--')
    txt = f'''
    $p$={regr_dict['pval']:.2}
    $r^2$={regr_dict['r_squared']**2:.2}
    '''
    ax.text(-10, -6.5, txt)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_xlim(plot_lim)
    ax.set_ylim(plot_lim)
    ax.set_yticks(fig6_tick)
    ax.set_xticks(fig6_tick)
    ax.set_yticklabels(fig6_tick)
    if j<3:
        for label in ax.xaxis.get_ticklabels():
                label.set_visible(False)
    if j%3 != 0:
        for label in ax.yaxis.get_ticklabels():
            label.set_visible(False)

axes[0,0].set_title('Kim and O\'Neil (1997)', fontsize=15)
axes[0,0].set_ylabel(u'$\delta^{18}$O Modeled \nEqulibirum', fontsize=15)

axes[0,1].set_title('Affek and Zaarur (2014)', fontsize=15)
axes[0,2].set_title('Coplen (2002)', fontsize=15)
axes[1,0].set_ylabel(f'{d18O} Modeled \nRayleigh', fontsize=15)

fig6.text(0.5, -0.02, f'{d18O} Calcite Measured', fontsize=15, ha='center')
plt.tight_layout()
plt.show()
plt.savefig('Carlson2020resubmit_fig06.png', format='png', dpi=600)
# %% Figure 7
# Model Comparison timeseries
mod_affek_D18OPDB = (calcite_data['d18OCalcite_Modeled_Affek']
                     - calcite_data['δ18OWaterPDB'])
calcite_data['Mod_Affek_Δ18OPDB'] = mod_affek_D18OPDB

fig, axes = plt.subplots(2,2, figsize=(10,5.5))
for i in range(0,len(plot_table)):
    ax = axes[i//2,i%2]
    plotname = plot_table.DripName.iloc[i]
    plotmarker = plot_table.Marker.iloc[i]
    plotcolor = plot_table.Color.iloc[i]
    set_to_plot = calcite_data[calcite_data.SiteName == plotname]
    set_to_plot = set_to_plot.sort_values('MidPoint')
    set_to_plot['mdate'] =mdates.date2num(set_to_plot['MidPoint'])
    ax.set_title(f'{plotname}')
    ax.errorbar(x=set_to_plot['mdate'],
                y=set_to_plot['Δ18OPDB'],
                yerr=set_to_plot['Δ18OPDB_Err'],
                linestyle='-',
                color='xkcd:cement',
                marker=plotmarker,
                markersize=4,
                linewidth=1,
                label='Measured')
    equil_names = ['Coplen', 'Affek', 'Kim']
    
    for eq_name in equil_names:
        eq_D18OPDB = (set_to_plot[f'equilib_calcite_δ18O_{eq_name.lower()}']
                     - set_to_plot['δ18OWaterPDB'])      
        eq_col = f'Equilib_{eq_name}_Δ18OPDB'
        set_to_plot[eq_col] = eq_D18OPDB

        ax.plot(set_to_plot['mdate'],
                set_to_plot[eq_col],
                linestyle='-',
                color='xkcd:light gray',
                marker=plotmarker,\
                    
                markersize=2,
                linewidth=1,
                label=None)
        ax.text(set_to_plot['mdate'].iloc[-1] + 60,
                set_to_plot[eq_col].iloc[-1],
                eq_name,
                color='xkcd:cement',
                ha='left',
                va='center')
    ax.plot(set_to_plot['mdate'],
            set_to_plot['Mod_Affek_Δ18OPDB'],
            linestyle='-',
            color=plotcolor,
            marker=plotmarker,
            markersize=4,
            linewidth=1,
            label='ISOLUTION')
    rmse = (set_to_plot['Mod_Affek_Δ18OPDB']
            - set_to_plot['Δ18OPDB']).pow(2).mean()**0.5
    ax.text(np.mean(mdates.date2num(start_plot)) + 60,
            30.7,
            f'RMSE = {rmse:00.02f}',
            ha='left',
            va='center')
    ax.legend(fontsize=8)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y')) 
    ax.set_xlim([start_plot, end_plot])
    ax.set_ylim([26, 31])
    ax.tick_params(axis='both', which='major', labelsize=10)
    
    plt.setp(ax.get_xticklabels(), rotation=0, ha='center')
    if i%2 != 0:
        for label in ax.yaxis.get_ticklabels():
            label.set_visible(False)
    else:
        ax.set_ylabel(f'{cap18O} (PDB)')
    if i//2 == 0:
        for label in ax.xaxis.get_ticklabels():
            label.set_visible(False)
plt.tight_layout()
plt.show()
plt.savefig('Carlson2020resubmit_fig07.png', format='png', dpi=600)
# %% Figure 8
# d13C vs d18O


fig, ax = plt.subplots(1,1, figsize=(5.5,5.5))
for i in range(0,len(plot_table)):
    plotname = plot_table.DripName.iloc[i]
    plotmarker = plot_table.Marker.iloc[i]
    plotcolor = plot_table.Color.iloc[i]
    set_to_plot = calcite_data[calcite_data.SiteName == plotname]
    
    x = set_to_plot['Calcite_δ18O_VPDB']
    y = set_to_plot['Calcite_δ13C_VPDB']
    x = x[x.notnull() & y.notnull()]
    y = y[x.notnull() & y.notnull()]
    func = lambda x, a, b: a * x + b
    band_corr_dict = oxystats.correlate_with_uncertainty_band(x,
                                                              y,
                                                              func,
                                                              alpha=0.1,
                                                              band_type='confidence')
    m = band_corr_dict['popt'][0]
    b = band_corr_dict['popt'][1]
    m_s = band_corr_dict['sigmas'][0]
    b_s = band_corr_dict['sigmas'][1]
    f_1 = lambda x: m * x + b
    regr_dict = oxystats.correlate_with_ar_correction(f_1(x), y)
    print(model)
    print(band_corr_dict['popt'])
    print(band_corr_dict['sigmas'])
    print(regr_dict)
    ftext = f'{d13C}$={m:0.2f}*${d18O}${b:+0.2f}$'
    rtext = f'$r^2={regr_dict["r_squared"]:0.2f}$'
    ptext = f'$p={regr_dict["pval"]:0.2e}$'
    label_text = f'{plotname}:    {ftext}    {rtext}    {ptext}'
    ax.plot(x,
            f_1(x),
            linestyle='-',
            color=plotcolor,
            marker='None',
            markersize=4,
            linewidth=1,
            label=label_text)    
    ax.errorbar(x=set_to_plot['Calcite_δ18O_VPDB'],
                y=set_to_plot['Calcite_δ13C_VPDB'],
                xerr=set_to_plot['Calcite_δ18O_Err_VPDB'],
                yerr=set_to_plot['Calcite_δ13C_Err_VPDB'],
                linestyle='None',
                color=plotcolor,
                marker=plotmarker,
                markersize=4,
                linewidth=1,
                label=None)
    ax.legend(fontsize=8, loc='upper center', bbox_to_anchor = (0.5,-0.15))
    #ax.set_xlim([-13, -5.5])
    #ax.set_ylim([-13, -5.5])
    ax.tick_params(axis='both', which='major', labelsize=10)

    plt.setp(ax.get_xticklabels(), rotation=0, ha='center')
    plt.savefig('PCO2.png', format='png', dpi=600)
    ax.set_xlabel(f'Calcite {d18O} (PDB)')
    ax.set_ylabel(f'Calcite {d13C} (PDB)')
plt.tight_layout()
plt.show()
plt.savefig('Carlson2020resubmit_fig08.png', format='png', dpi=600)
# %% Figure 9
# d18O Departure versus CO2
fig7, axes = plt.subplots(2,2, figsize=(10,5.5))
co2_to_plot = FieldCO2Continous['2012-01-01':]
co2_to_plot['mdate'] = mdates.date2num(co2_to_plot['Datetime'])

for i in range(0,len(plot_table)):
    ax = axes[i//2,i%2]
    plotname = plot_table.DripName.iloc[i]
    plotmarker = plot_table.Marker.iloc[i]
    plotcolor = plot_table.Color.iloc[i]
    set_to_plot = calcite_data[calcite_data.SiteName == plotname]
    set_to_plot.set_index('MidPoint', inplace=True)
    set_to_plot = set_to_plot['2012-01-01':]
    set_to_plot.reset_index(inplace=True, drop=False)
    set_to_plot = set_to_plot.sort_values('MidPoint')
    set_to_plot['mdate'] =mdates.date2num(set_to_plot['MidPoint'])
    
    set_to_plot['departure'] = (set_to_plot['Calcite_δ18O_VPDB'] 
                                - set_to_plot['equilib_calcite_δ18O_affek'])
    ax.set_title(f'{plotname}')
    ax.plot(co2_to_plot['mdate'],
            co2_to_plot['CO2'],
            linestyle = '-',
            color = 'k',
            marker = None,
            linewidth = 1,
            label = None)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%Y')) 
    ax.set_xlim([co2_start_plot, end_plot])
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax2 = ax.twinx()
    ax2.plot(set_to_plot['mdate'], 
            set_to_plot['departure'],
            linestyle='-',
            color=plotcolor,
            marker=plotmarker,
            markeredgewidth=1,
            linewidth=1,
            markeredgecolor=plotcolor,
            label=plotname,
            markersize=5)    
    ax2.set_ylim([-1, 2.5])
    ax2.invert_yaxis()
    plt.setp(ax2.get_yticklabels(),color=plotcolor)

    plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
    plt.savefig('PCO2.png', format='png', dpi=600)
    y2str =  f'''
    Departure from 
    Affek and Zaarur 
    Equilibrium ({permil} PDB)
    '''
    if i%2 != 0:
        ax2.set_ylabel(y2str)
    else:
        ax.set_ylabel('$pCO_2$ (ppmv)')
    if i//2 == 0:
        for label in ax.xaxis.get_ticklabels():
            label.set_visible(False)
plt.tight_layout()
plt.show()
plt.savefig('Carlson2020resubmit_fig09.png', format='png', dpi=600)
# %% Figure 10
# d13C vs CO2
fig7, axes = plt.subplots(2,2, figsize=(10,5.5))
co2_to_plot = FieldCO2Continous['2012-01-01':]
co2_to_plot['mdate'] = mdates.date2num(co2_to_plot['Datetime'])

for i in range(0,len(plot_table)):
    ax = axes[i//2,i%2]
    plotname = plot_table.DripName.iloc[i]
    plotmarker = plot_table.Marker.iloc[i]
    plotcolor = plot_table.Color.iloc[i]
    set_to_plot = calcite_data[calcite_data.SiteName == plotname]
    set_to_plot.set_index('MidPoint', inplace=True)
    set_to_plot = set_to_plot['2012-01-01':]
    set_to_plot.reset_index(inplace=True, drop=False)
    set_to_plot = set_to_plot.sort_values('MidPoint')
    set_to_plot['mdate'] =mdates.date2num(set_to_plot['MidPoint'])

    ax.set_title(f'{plotname}')
    ax.plot(co2_to_plot['mdate'],
            co2_to_plot['CO2'],
            linestyle = '-',
            color = 'k',
            marker = None,
            linewidth = 1,
            label = None)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%Y')) 
    ax.set_xlim([co2_start_plot, end_plot])
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax2 = ax.twinx()
    ax2.plot(set_to_plot['mdate'], 
            set_to_plot['Calcite_δ13C_VPDB'],
            linestyle='-',
            color=plotcolor,
            marker=plotmarker,
            markeredgewidth=1,
            linewidth=1,
            markeredgecolor=plotcolor,
            label=plotname,
            markersize=5)    
    ax2.set_ylim([-15.5, -6.5])
    ax2.invert_yaxis()
    plt.setp(ax2.get_yticklabels(),color=plotcolor)

    plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
    plt.savefig('PCO2.png', format='png', dpi=600)
    y2str =  f'''
    Calcite{d13C}
    ({permil} PDB)
    '''
    if i%2 != 0:
        ax2.set_ylabel(y2str)
    else:
        ax.set_ylabel('$pCO_2$ (ppmv)')
    if i//2 == 0:
        for label in ax.xaxis.get_ticklabels():
            label.set_visible(False)
plt.tight_layout()
plt.show()
plt.savefig('Carlson2020resubmit_fig10.png', format='png', dpi=600)
# %% Figure A.1
# Published 1000lnalpha values
fig, ax = plt.subplots(1, 1, figsize=(5.6,5.6))
external_data['inv_TK'] = 1000 / t_to_k(external_data['WaterTemp'].values, 'C')
calcite_data['inv_TK'] = 1000 / t_to_k(calcite_data['WaterTemp'].values, 'C')
markers = ['o','D','v','^', '>','<', '*', '+', 'x']
studies = external_data.groupby('PublicationKey')
for i, (study, study_data) in enumerate(studies):
    ax.plot(study_data['inv_TK'],
            study_data['1000lnα'],
            linestyle='none',
#            color='xkcd:gray',
            marker=markers[i//len(markers)],
            markersize=3,
            label=study,
            )
ax.plot(calcite_data['inv_TK'],
        calcite_data['meas_1000lnα'],
        linestyle='none',
        color='xkcd:dark gray',
        marker='o',
        markersize=4,
        label='This Study',
        )
xlim = ax.get_xlim()
ax.plot(xlim,
        affek_1000lna(1000/np.array(xlim)),
        linestyle='-',
        color='k',
        marker=None,
        label='Affek and Zaarur Line',
        )
ax.plot(xlim,
        coplen_1000lna(1000/np.array(xlim)),
        linestyle='--',
        color='k',
        marker=None,
        label='Coplen Line',
        )
ax.plot(xlim,
        kim_1000lna(1000/np.array(xlim)),
        linestyle=':',
        color='k',
        marker=None,
        label="Kim and O'Neil Line",
        )
ax.plot(xlim,
        tremaine_1000lna(1000/np.array(xlim)),
        linestyle='-.',
        color='k',
        marker=None,
        label="Tremaine Line",
        )
ax.set_xlim(xlim)
ax.set_ylabel('$1000ln{\alpha}$')
ax.set_xlabel('$1000/T [K]$')
ax.legend(fontsize=8,
          loc='center left', 
          bbox_to_anchor=(1.1,0.5),
          ncol=2)
ax2 = ax.twiny()
ax2.set_xlabel('$T [{\degree}C]$')
ax2.set_xlim(xlim)
ax2ticks = [1000/t_to_k(t, 'C')
            for t 
            in range(0,40,5)]

ax2.set_xticks(ax2ticks)
ax2.set_xticklabels([str(t) for t in range(0,40,5)])         
plt.tight_layout()
plt.show()
plt.savefig('Carlson2020resubmit_figA1.png', format='png', dpi=600) 
# %% Figure A.2
rooms = ['Midslide', 'Shakey Room', 'Shakey Room', 'Stumpy Room']
fig7, axes = plt.subplots(2,2, figsize=(10,5.5))
co2_to_plot = FieldCO2Continous['2012-01-01':]
co2_to_plot['mdate'] = mdates.date2num(co2_to_plot['Datetime'])
for i in range(0,len(plot_table)):
    ax = axes[i//2,i%2]
    plotname = plot_table.DripName.iloc[i]
    plotmarker = plot_table.Marker.iloc[i]
    plotcolor = plot_table.Color.iloc[i]
    room = rooms[i]
    set_to_plot = FieldCO2Spot[FieldCO2Spot.SiteName == room]
    set_to_plot.set_index('Datetime', inplace=True)
    set_to_plot = set_to_plot['2012-01-01':]
    set_to_plot.reset_index(inplace=True, drop=False)
    set_to_plot = set_to_plot.sort_values('Datetime')
    set_to_plot['mdate'] =mdates.date2num(set_to_plot['Datetime'])
    ax.set_title(f'{plotname}: {room}')
    ax.plot(co2_to_plot['mdate'],
            co2_to_plot['CO2'],
            linestyle = '-',
            color = 'k',
            marker = None,
            linewidth = 1,
            label = None)
    
    if i != 3:
        ax.plot(set_to_plot['mdate'], 
                set_to_plot['CO2'],
                linestyle = '-',
                color = plotcolor,
                marker = plotmarker,
                markeredgewidth = 1,
                linewidth = 1,
                markeredgecolor = plotcolor,
                label = plotname,
                markersize = 5)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%Y')) 
    ax.set_xlim([co2_start_plot, end_plot])
    ax.tick_params(axis='both', which='major', labelsize=10)
    
    plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
    plt.savefig('PCO2.png', format='png', dpi=600)
    if i%2 != 0:
        for label in ax.yaxis.get_ticklabels():
            label.set_visible(False)
    else:
        ax.set_ylabel('$pCO_2$ (ppmv)')
    if i//2 == 0:
        for label in ax.xaxis.get_ticklabels():
            label.set_visible(False)
#plt.tight_layout()
plt.show()
plt.savefig('Carlson2020resubmit_figA2.png', format='png', dpi=600)