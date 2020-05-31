# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 20:06:38 2020
Statistical and Data Analysis Functions for Guam 2020 Oxygen isotopes paper.
Run in Python 3
@author: peter carlson
"""
import pandas as pd
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from statsmodels.tsa.arima_model import ARMA


# Set Constants
MOLAR_MASS_CA_G_MOL = 40.08
MOLAR_MASS_CALCITE_G_MOL = 100.09 
SECONDS_PER_DAY = 60 * 60 * 24
PLATE_AREA_CM2 = 100

#Set Uncertainties (1 sigma)
PLATE_AREA_ERR_CM2 = 6
PLATE_WEIGHT_ERR_MG = 0.9
STD_PLATE_G = 55.38745
SAMPLE_TIME_ERR_MIN = 5


#Constant Conversions
PLATE_AREA_M2 = PLATE_AREA_CM2 / 1E4
MOLAR_MASS_CALCITE_G_UMOL = MOLAR_MASS_CALCITE_G_MOL / 1E6

# Scale Factor
#g calcite per day per plate => umol calcite per m2 per hour
RC_SCALE = 1/(24. * PLATE_AREA_M2 * MOLAR_MASS_CALCITE_G_UMOL)

# Error Conversions
PLATE_AREA_ERR_M2 = PLATE_AREA_ERR_CM2 / 1E4
PLATE_WEIGHT_ERR_G = PLATE_WEIGHT_ERR_MG / 1E3
SAMPLE_TIME_ERR_DAYS = SAMPLE_TIME_ERR_MIN * 60 / SECONDS_PER_DAY

# Error Propagation
CALCITE_WEIGHT_ERR_G = (2*(PLATE_WEIGHT_ERR_G**2))**(0.5)
DEPLOY_TIME_ERR_DAYS = (2*(SAMPLE_TIME_ERR_DAYS**2))**(0.5)

def load_cave_data():
    """
    Reads and processes cave data from Carlson et al 2020 Supplemental Data.
    Shows error propagation for derived properties
    
    Returns
    -------
    cave_results_df : pandas dataframe
    """
    # Read Data 
    cave_results_df = pd.read_excel(
            './Carlson2020_Supplement.xlsx',
            sheet_name='TableA2_Calcite_&_Model_Results',
            skiprows=[0,2]
    )
    print(cave_results_df.columns)
    cave_results_df = add_drip_columns(cave_results_df)
    cave_results_df = add_calcite_growth_columns(cave_results_df)
    
    
    return cave_results_df

def add_calcite_growth_columns(cave_results_df):
    """
    Calculates Calcite Growth, Rc, and log Rc
    Calculates Growth Rate Error, Rc Error, and log(Rc) Error
    Adds all as columns ad returns input df.

    Parameters
    ----------
    cave_results_df : pandas dataframe
    
    Returns
    -------
    cave_results_df : pandas dataframe
    """
    # Data Conversions
    cave_results_df['Calcite_Growth'] = (cave_results_df['GrowthRate'] *
                   cave_results_df['DaysDeploy'])
    cave_results_df['Rc'] = RC_SCALE*cave_results_df['GrowthRate'] # Eq. 11
    cave_results_df['logRc'] = np.log10(cave_results_df['Rc'])
    
    # Error Propagation
    cave_results_df['GrowthRate_Err'] = cave_results_df['GrowthRate'] * (
            (CALCITE_WEIGHT_ERR_G / cave_results_df['Calcite_Growth']).pow(2)
            + (SAMPLE_TIME_ERR_DAYS / cave_results_df['DaysDeploy']).pow(2)
            ).pow(0.5)
    cave_results_df['Rc_err'] = cave_results_df['Rc']*(
            (cave_results_df['GrowthRate_Err'] / cave_results_df['GrowthRate']).pow(2)
            + (PLATE_AREA_ERR_CM2 / PLATE_AREA_CM2)**2
            ).pow(0.5)
    cave_results_df['logRc_err'] = cave_results_df['Rc_err'] / cave_results_df['Rc']
    
    return cave_results_df

def add_drip_columns(cave_results_df):
    """
    Calculates Drips per min
    Calculates Drip Rate Error
    Adds all as columns ad returns input df.

    Parameters
    ----------
    cave_results_df : pandas dataframe

    Returns
    -------
    cave_results_df : pandas dataframe
    """
    # Data Conversions
    cave_results_df['DripsPerMin'] = 60.0/cave_results_df['DripInterval']

    # Error Propagation
    cave_results_df['DripRate_Err'] = (60*cave_results_df['DripInterval_err'] / 
                   cave_results_df['DripInterval'].pow(2))
    
    return cave_results_df

def correlate_vs_driprate(cave_results_df):
    """
    Calculates Drips per min
    Calculates Drip Rate Error
    Adds all as columns ad returns input df.

    Parameters
    ----------
    cave_results_df : pandas dataframe

    """
    # Filter data for valid drip rate and fractionation values
    to_regress = cave_results_df[
            (cave_results_df['Δ18OPDB'].notnull()) 
            & (cave_results_df['DripsPerMin'].notnull())
            & (np.isfinite(cave_results_df['DripRate_Err']))
            ].copy()
    
    # Set correlation function
    correlation_function = lambda x, x_b, m_1, b_1, m_2: (
            piecewise_linear_continuous(x, x_b, m_1, b_1, m_2))
    x = to_regress['DripsPerMin'].values
    y = to_regress['Δ18OPDB'].values
    band_corr_dict = correlate_with_uncertainty_band(
            x,
            y,
            correlation_function,
            alpha=0.1,
            band_type = 'prediction'
            )
    
    x_b = band_corr_dict['popt'][0]
    m_1 = band_corr_dict['popt'][1]
    b_1 = band_corr_dict['popt'][2]
    m_2 = band_corr_dict['popt'][3]
    dm = (m_1 - m_2)
    b_2 = dm*x_b + b_1
    
    x_b_s = band_corr_dict['sigmas'][0]
    m_1_s = band_corr_dict['sigmas'][1]
    b_1_s = band_corr_dict['sigmas'][2]
    m_2_s = band_corr_dict['sigmas'][3]
    
    dm_s = (m_1_s**2 - m_2_s**2)**(0.5)
    
    b_2_s = (((dm*x_b)**2)*(
                    (dm_s/dm)**2 + (x_b_s/x_b)**2)
             + b_1_s**2)**(0.5)    
    print(f'''
    ###################################################################
    # CALCITE-WATER OXYGEN ISOTOPE FRACTIONATION FACTORS v. DRIP RATE #
    ###################################################################
        
        Piecewise linear correlation:
        Break Point = {x_b: 0.2f} +/- {x_b_s: 0.2f} drips per minute
    
        f(x) = {{
                ({m_1:0.2f} +/- {m_1_s:0.2f})*x + ({b_1:0.2f} +/- {b_1_s:0.2f}),\
    x <= {x_b:0.2f}
                ({m_2:0.2f} +/- {m_2_s:0.2f})*x + ({b_2:0.2f} +/- {b_2_s:0.2f}),\
    x >  {x_b:0.2f}
               }}''')
    
    x_1 = to_regress[to_regress['DripsPerMin'] <= x_b]['DripsPerMin'].values
    y_1 = to_regress[to_regress['DripsPerMin'] <= x_b]['Δ18OPDB'].values
    f_1 = lambda x: m_1*x + b_1

    ar_corr_dict= correlate_with_ar_correction(f_1(x_1),y_1)

    print(f'''
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~          
    Below Break Point:
        
        f_low(x) = ({m_1:0.2f} +/- {m_1_s:0.2f})*x + ({b_1:0.2f} +/- {b_1_s:0.2f})
    
        r_squared = {ar_corr_dict['r_squared']:0.2f}
        p = {ar_corr_dict['pval']:0.2e}''')
    
    
    x_2 = to_regress[to_regress['DripsPerMin'] > x_b]['DripsPerMin'].values
    y_2 = to_regress[to_regress['DripsPerMin'] > x_b]['Δ18OPDB'].values
    f_2 = lambda x: m_2*x + b_2
    
    ar_corr_dict = correlate_with_ar_correction(f_2(x_2),y_2)
    
    print(f'''
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~          
    Above Break Point:
        
        f_high(x) = ({m_2:0.2f} +/- {m_2_s:0.2f})*x + ({b_2:0.2f} +/- {b_2_s:0.2f})
    
        r_squared = {ar_corr_dict['r_squared']:0.2f}
        p = {ar_corr_dict['pval']:0.2e}''')
        
    return None
        
def correlate_vs_rc(cave_results_df):
    """
    Calculates Drips per min
    Calculates Drip Rate Error
    Adds all as columns ad returns input df.

    Parameters
    ----------
    cave_results_df : pandas dataframe
    """    
    to_regress = cave_results_df[
            (pd.notnull(cave_results_df['Δ18OPDB']))
            & (pd.notnull(cave_results_df.GrowthRate))
            & (cave_results_df.GrowthRate > 0)
            ].copy()
    x = to_regress['logRc'].values
    y = to_regress['Δ18OPDB'].values
    func = lambda x, a, b: a*x + b
    band_corr_dict = correlate_with_uncertainty_band(
            x,
            y,
            func,
            alpha=0.1,
            band_type = 'confidence')
    
    m = band_corr_dict['popt'][0]
    b = band_corr_dict['popt'][1]
    m_s = band_corr_dict['sigmas'][0]
    b_s = band_corr_dict['sigmas'][1]
    
    f_1 = lambda x: m*x + b
    ar_corr_dict= correlate_with_ar_correction(f_1(x),y)

    print(f'''
    ###################################################################
    #  CALCITE-WATER OXYGEN ISOTOPE FRACTIONATION FACTORS v. log(Rc)  #
    ###################################################################
        
        Linear correlation:
            
        f(x) = ({m:0.2f} +/- {m_s:0.2f})*x + ({b:0.2f} +/- {b_s:0.2f})
        
        r_squared = {ar_corr_dict['r_squared']:0.2f}
        p = {ar_corr_dict['pval']:0.2e}''')
        
    return None

    
def piecewise_linear_continuous(x, x_b, m_1, b_1, m_2):
    """
    Eq. 10. 2-Part Piecewise Linear Continuous Function
    Used to fit D18O_cc-w fractionation versus drip rate. 
    Describes two lines that intersect at x_b.
    
    Parameters
    ----------
    x : numpy array
        independent variable
    x_b : float
        breakpoint for the piecewise linear function.
    m_1 : float
        slope of the line left of x_b
    b_1 : float
        intercept of the line left of x_b
    m_2 : float
        slope of the line right of x_b
    
    Returns
    -------
    y : numpy array
        dependent variable
    """
    
    y=x*0 # initialize y
    
    # Function Left of x_b
    y[x<=x_b] = m_1*x[x<=x_b] + b_1
    
    # Function right of x_b
    b_2 = (m_1 - m_2)*(x_b)
    y[x>x_b] = m_2*x[x>x_b] + b_1 + b_2
    
    return y

def correlate_with_uncertainty_band(x, y, func, alpha=0.05, x_hat_res=200,
                                    band_type='confidence'):
    """
    Creates best-fit parameters and upper and lower confidence bands
    for a given set of x and y data and a user-supplied function.
    
    Parameters
    ----------
    x : numpy array
        independent variable data or model
    y : numpy array
        dependent variable data
    func : function
        anonymous function to use to correlate data.
    alpha : float, optional
        statistical significance threshold value. Symmetrical (0.05=0.95)
        defaul: 0.05
    x_hat_res : int, optional
        number of subdivisions in x for output band data
        default: 200
    band_type : string, optional
        type of uncertainty band. Must be 'confidence' or 'prediction'
        default: 'confidence'
    
    Returns
    -------
    dict 
    keys:
        popt : numpy array
            array of optimal parameter values for the given function
        sigmas : numpy array
            1-sigma uncertainty associated with each parameter in popt
        p_upper-popt : 
        x_hat : numpy array
            
        uncertainty_band : numpy array. 
            a (x_hat_res,2) array. 
            Uncertainty_band[:,0] is lower uncertainty.
            Uncertainty_band[:,1] is upper uncertainty.
        r_squared : float
            Pearson's r^2. Note that this is only useful if func is linear.
    """
    alpha = min([alpha, 1-alpha])   
    
    popt, covar = curve_fit(func, x, y)
    
    sigmas = np.sqrt(np.diagonal(covar))
    n = len(x)
    xbar = np.mean(x)
    ybar = np.mean(y)
    
    sum_of_squares = lambda x : np.dot(x,x)
    SSxx = sum_of_squares(x - xbar) # sum of squares in x
    SSyy = sum_of_squares(y - ybar) # sum of squares in y
    SSe = sum_of_squares(y - func(x, *popt)) # sum of squares of the error
    r_squared = 1 - SSe/SSyy
    
    
    t = stats.t.ppf(1-alpha/2,n-len(popt), )
    
    x_hat = np.linspace(np.min(x), np.max(x), x_hat_res)
    y_hat = func(x_hat, *popt)
    
    
    if band_type.lower() == 'prediction':
        sigma_xy = np.sqrt(SSe/(n-len(popt)))
        error_func = lambda x: t * sigma_xy * np.sqrt(
                1 + 1/n + (np.square(x - xbar))/SSxx)        
        p_lower = np.empty_like(popt)
        p_upper = np.empty_like(popt)
        bound_lower = y_hat - error_func(x_hat)
        bound_upper = y_hat + error_func(x_hat)
        
    elif band_type.lower() == 'confidence':
        p_lower = popt - sigmas*t
        p_upper = popt + sigmas*t
        bound_lower = func(x_hat, *p_lower)
        bound_upper = func(x_hat, *p_upper)
        
    else:
        raise ValueError(f'Unknown uncertainty band type: \"{band_type}\"')
    
    uncertainty_band = (bound_lower,bound_upper)
    parameter_bounds = (p_lower,p_upper)
    return_dict =  {
            'popt': popt,
            'sigmas': sigmas,
            'parameter_bounds': parameter_bounds,
            'x_hat': x_hat,
            'uncertainty_band': uncertainty_band,
            'r_squared': r_squared
            }
    return return_dict

def correlate_with_ar_correction(x,y):
    """
    linear correlation, with a correction for autoregressive parameter phi.
    For use with timeseries data.
    
    After Hu et al., 2017:
    Hu, J., Emile-Geay, J., Partin, J. 2017. Correlation-based interpretations 
    of paleoclimate data—where statistics meet past climates. 
    Earth and Plan. Sci. Let. 459:362-371.
    
    Parameters
    ----------
    x : numpy array
        independent variable data or model
    y : numpy array
        dependent variable data
     
    Returns
    -------
    dict
    keys:
        r_squared : float
            Pearson's r^2
        pval : float
            signficance, corrected for autoregressive characteristics of x and y
        phi_vals : tuple
            phi_vals[0] is the first-order AR coefficient for X.
            phi_vals[1] is the same for y
    """

    n = len(x)
    
    sum_of_squares = lambda x : np.dot(x,x)
    
    xbar = np.mean(x)
    ybar = np.mean(y)
    
    SSxx = sum_of_squares(x - xbar) # sum of squares in x
    SSyy = sum_of_squares(y - ybar) # sum of squares in y
    
    r =  np.dot(x - xbar,y - ybar)/np.sqrt(SSxx*SSyy)

    # Fit AR(1) model to y and to x
    xmod = ARMA(x, order=(1,0,0))
    x_coef = xmod.fit().arparams
    ymod = ARMA(y, order=(1,0,0))
    y_coef = ymod.fit().arparams
    
    # For effective degrees of freedom
    neff = n*(1-x_coef*y_coef)/(1+x_coef*y_coef)
    if neff <3:
        neff = 3   
    phi_vals = (x_coef[0], y_coef[0])
    tval = r/np.sqrt(1-r**2)*np.sqrt(neff-2)   
    pval = stats.t.sf(abs(tval),neff-2)*2
    
    return {'r_squared': r**2, 'pval': pval[0], 'phi_vals': phi_vals}

if __name__ == '__main__':
    cave_results_df = load_cave_data()
    correlate_vs_driprate(cave_results_df)
    correlate_vs_rc(cave_results_df)