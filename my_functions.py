# -*- coding: utf-8 -*-
"""
Created on Thu May 12 18:13:04 2016

@author: anoronha
"""

import pandas as pd
from pandas import DataFrame
import numpy as np
import os
from pylab import *
import datetime


def plot_all(series, to_plot, x_val, y_val, fill,line, axis_num, y_error_bar,y_direct):

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