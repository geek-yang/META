# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Statistical Operator for Climate Data
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.07.26
Last Update     : 2018.08.14
Contributor     :
Description     : This module provides several methods to perform statistical
                  analysis on MET and all kinds of fields.
Return Values   : pngs
Caveat!         :
"""

import numpy as np
import scipy
from scipy import stats
import cartopy
import os
import iris
#import matplotlib
import seaborn as sns
#import matplotlib.pyplot as plt

class operator:
    def __init__(self, var):
        """
        Statistical operations on climate data.
        param var: imput time series
        param outpath: the path for the output files
        """
        self.var = var

    def anomaly(self, Dim_month=True):
        """
        Remove seasonal cycling for monthly data.
        param Dim_month: there are two modes for removing the seasonal cycling
        -True (default) input time series have month dimension [year,month,...]
        -False input time series have only 1 dimension for time
        param white_var: time series without seasonal cycling
        """
        # white refers to the time series without seasonal cycling
        white_var = np.zeros(self.var.shape, dtype=float)
        #switch mode
        if Dim_month == True:
            print ('The input data has the dimension of month.')
            # check the dimension of input
            seansonal_cycle_var = np.mean(self.var, axis=0)
            if self.var.ndim == 2:
                t, m = white_var.shape
                for i in np.arange(t):
                        white_var[i,:] = self.var[i,:] - seansonal_cycle_var[:]
                # re-arrange into single time series - without month dimension
                white_var = white_var.reshape(t*m)
            elif self.var.ndim == 3:
                t, m, y = white_var.shape
                for i in np.arange(t):
                        white_var[i,:,:] = self.var[i,:,:] - seansonal_cycle_var[:]
                # re-arrange into single time series - without month dimension
                white_var = white_var.reshape(t*m,y)
            elif self.var.ndim == 4:
                t, m, y, x = white_var.shape
                for i in np.arange(t):
                        white_var[i,:,:,:] = self.var[i,:,:,:] - seansonal_cycle_var[:]
                # re-arrange into single time series - without month dimension
                white_var = white_var.reshape(t*m,y,x)
            else:
                raise IOError("This module can not work with any array with a \
                              dimension higher than 4!")
        else:
            print ('The input data does not have the dimension of month.')
            if self.var.ndim == 1:
                for i in np.arange(12):
                    seansonal_cycle_var = np.mean(self.var[i::12],axis=0)
                    white_var[i::12] = self.var[i::12] - seansonal_cycle_var
            elif self.var.ndim == 2:
                for i in np.arange(12):
                    seansonal_cycle_var = np.mean(self.var[i::12,:],axis=0)
                    white_var[i::12,:] = self.var[i::12,:] - seansonal_cycle_var
            elif self.var.ndim == 3:
                for i in np.arange(12):
                    seansonal_cycle_var = np.mean(self.var[i::12,:,:],axis=0)
                    white_var[i::12,:,:] = self.var[i::12,:,:] - seansonal_cycle_var
            else:
                raise IOError("This module can not work with any array with a \
                              dimension higher than 3!")
        self._anomaly = white_var
        
        print ("The output anomaly time series only contains one dimension for time!")
        
        return self._anomaly

    def detrend(self, order=2, obj='anomaly'):
        """
        Detrend time series through polynomial fit.
        param series: input time series, either 1D or 2/3D
        param order: order of polynomial for fitting
        param obj: objects for detrending, two options available
        -'anomaly' (default) the time series of anomaly will be detrended
        -'original' the original input time series will be detrended
        """
        if obj == 'anomaly':
            series = self._anomaly
        elif obj == 'original':
            print ("Make sure that the input time series has only 1 dimension for time!")
            series = self.var
        else:
            raise IOError("Please choose the right input mode for detrending!")
        # check the dimension of input
        if series.ndim == 1:
            polynomial = np.polyfit(np.arange(len(series)), series, order)
            poly_fit = np.poly1d(polynomial)
            poly_fit_var = poly_fit(np.arange(len(series)))
        elif series.ndim == 2:
            poly_fit_var = np.zeros(series.shape, dtype=float)
            t, y = poly_fit_var.shape
            for i in np.arange(y):
                polynomial = np.polyfit(np.arange(t), series[:,i], order)
                poly_fit = np.poly1d(polynomial)
                poly_fit_var[:,i] = poly_fit(np.arange(t))
        elif series.ndim == 3:
            poly_fit_var = np.zeros(series.shape, dtype=float)
            t, y, x = poly_fit_var.shape
            for i in np.arange(y):
                for j in np.arange(x):
                    polynomial = np.polyfit(np.arange(t), series[:,i,j], order)
                    poly_fit = np.poly1d(polynomial)
                    poly_fit_var[:,i,j] = poly_fit(np.arange(t))
        else:
            raise IOError("This module can not work with any array with a \
                            dimension higher than 3!")
        self._polyfit = poly_fit_var
        self._detrend = series - self._polyfit

        return self._detrend
    
    def lowpass(self, window=60, obj='anomaly'):
        """
        Apply low pass filter to the time series. The function gives running mean
        for the point AT The End Of The Window!!
        param series: input time series, either 1D or 2/3D
        param window: time span for the running mean
        param obj: object for detrending, two options available
        -'anomaly' (default) apply low pass filter to the time series of anomaly
        -'original' apply lowpass filter to the original input time series
        -'detrend' apply lowpass filter to the detrended time series
        """
        if obj == 'anomaly':
            series = self._anomaly
        elif obj == 'original':
            series = self.var
        elif obj == 'detrend':
            series = self._detrend        
        # check the dimension of input
        if series.ndim == 1:
            t = len(series)
            running_mean = np.zeros(t-window+1, dtype=float)
            for i in np.arange(t-window+1):
                running_mean[i] = np.mean(series[i:i+window])
        elif series.ndim == 2:
            t, y  = series.shape
            running_mean = np.zeros((t-window+1, y), dtype=float)
            for i in np.arange(t-window+1):
                running_mean[i,:] = np.mean(series[i:i+window,:],0)
        elif series.ndim == 3:
            t, y, x = series.shape
            running_mean = np.zeros((t-window+1, y, x), dtype=float)
            for i in np.arange(t-window+1):
                running_mean[i,:,:] = np.mean(series[i:i+window,:,:],0)
        else:
            raise IOError("This module can not work with any array with a \
                            dimension higher than 3!")
        self._lowpass = running_mean
        
        return self._lowpass
    
    @staticmethod
    def linearRegress(var_x, var_y, lag=0):
        """
        Linear regression of input time series. Lead/lag regression can also be performed.
        param var_x: input time series, either 1D or 2D
        param var_y: input time series as the regression target, either 1D or 3D
        param lag: time unit for lead / lag regression, lag must be an integer
        """
        # check the dimensions of input time series
        if var_x.shape == var_y.shape:
            print("One time series is regressed on another.")
            if var_y.ndim == 2:
                if lag == 0:
                    t, y  = var_y.shape
                    slope = np.zeros(y, dtype=float)
                    r_value = np.zeros(y, dtype=float)
                    p_value = np.zeros(y, dtype=float)
                    for i in np.arange(y):
                        slope[i], _, r_value[i], p_value[i], _ = stats.linregress(var_x[:,i], var_y[:,i])
                elif lag > 0:
                    print ("This a regression with lead/lag analysis.")
                    t, y  = var_y.shape
                    lag_index = np.arange(-lag,lag+1,1)
                    slope = np.zeros((len(lag_index),y), dtype=float)
                    r_value = np.zeros((len(lag_index),y), dtype=float)
                    p_value = np.zeros((len(lag_index),y), dtype=float)
                    # regress
                    for i in np.arange(len(lag_index)):
                        for j in np.arange(y):
                            if lag_index[i]<0: # var_x lead var_y
                                slope[i,j], _, r_value[i,j], p_value[i,j], _ = stats.linregress(
                                               var_x[:lag_index[i],j], var_y[-lag_index[i]:,j])
                            elif lag_index[i]>0: # var_y lead var_x
                                slope[i,j], _, r_value[i,j], p_value[i,j], _ = stats.linregress(
                                               var_x[lag_index[i]:,j], var_y[:-lag_index[i],j])
                            else:
                                slope[i,j], _, r_value[i,j], p_value[i,j], _ = stats.linregress(
                                               var_x[:,j], var_y[:,j])
                else:
                    IOError("The lead / lag coefficient should be positive integers.")
            else:
                raise IOError("The input time series must have dimensions @time,latitude@!")
        elif var_y.ndim == 2 and var_x.ndim == 1:
            print("A time series is regressed on a 1D field.")
            if lag == 0:
                t, y  = var_y.shape
                slope = np.zeros(y, dtype=float)
                r_value = np.zeros(y, dtype=float)
                p_value = np.zeros(y, dtype=float)
                for i in np.arange(y):
                    slope[i], _, r_value[i], p_value[i], _ = stats.linregress(var_x[:], var_y[:,i])
            elif lag > 0:
                print ("This a regression with lead/lag analysis.")
                t, y  = var_y.shape
                lag_index = np.arange(-lag,lag+1,1)
                slope = np.zeros((len(lag_index),y), dtype=float)
                r_value = np.zeros((len(lag_index),y), dtype=float)
                p_value = np.zeros((len(lag_index),y), dtype=float)
                # regress
                for i in np.arange(len(lag_index)):
                    for j in np.arange(y):
                        if lag_index[i]<0: # var_x lead var_y
                            slope[i,j], _, r_value[i,j], p_value[i,j], _ = stats.linregress(
                                            var_x[:lag_index[i]], var_y[-lag_index[i]:,j])
                        elif lag_index[i]>0: # var_y lead var_x
                            slope[i,j], _, r_value[i,j], p_value[i,j], _ = stats.linregress(
                                            var_x[lag_index[i]:], var_y[:-lag_index[i],j])
                        else:
                            slope[i,j], _, r_value[i,j], p_value[i,j], _ = stats.linregress(
                                            var_x[:], var_y[:,j])                
        elif var_y.ndim == 3 and var_x.ndim == 1:
            print("A time series is regressed on a 2D field.")
            if lag == 0:
                t, y, x  = var_y.shape
                slope = np.zeros((y, x), dtype=float)
                r_value = np.zeros((y, x), dtype=float)
                p_value = np.zeros((y, x), dtype=float)
                for i in np.arange(y):
                    for j in np.arange(x):
                        slope[i,j], _, r_value[i,j], p_value[i,j], _ = stats.linregress(var_x, var_y[:,i,j])
            elif type(lag) == int:
                print ("This a regression with lead/lag analysis.")
                print ("Positive lag means 2nd input leads 1st, vice versa.")
                t, y, x  = var_y.shape
                slope = np.zeros((y, x), dtype=float)
                r_value = np.zeros((y, x), dtype=float)
                p_value = np.zeros((y, x), dtype=float)
                for i in np.arange(y):
                    for j in np.arange(x):
                        if lag > 0:
                            slope[i,j], _, r_value[i,j], p_value[i,j], _ = stats.linregress(var_x[lag:],
                                                                                            var_y[:-lag,i,j])
                        elif lag < 0:
                            slope[i,j], _, r_value[i,j], p_value[i,j], _ = stats.linregress(var_x[:lag],
                                                                                            var_y[-lag:,i,j])                           
            else:
                IOError("The lead / lag coefficient should be integers.")
        else:
            IOError("The dimensons of input time series are not supported.")
        
        return slope, r_value, p_value

    @staticmethod
    def seasons(series, Dim_month=False):
        """
        Extract time series for summer / winter from given series.
        The given time series should include the time series of all seasons.
        Here summer is June, July and August, and winter is December, January and
        February.
        param series: input time series containing the data for all seasons.
        param Dim_month: A check whether the time series include the dimension of month.
        """
        if Dim_month == True:
            if series.ndim == 2:
                t, m = series.shape
                series = series.reshape(t*m)
            elif series.ndim == 3:
                t, m, y = series.shape
                series = series.reshape(t*m, y)
            elif series.ndim == 4:
                t, m, y, x = series.shape
                series = series.reshape(t*m, y, x)
            else:
                raise IOError("This module can not work with any array with a \
                              dimension higher than 4!")
        else:
            pass
        # seperate summer and winter from the rest of the months
        if series.ndim == 1:
            t = len(series)
            series_summer = np.zeros(t/4,dtype=float)
            series_winter = np.zeros(t/4,dtype=float)
            # obtain summer time series
            series_summer[0::3] = series[5::12] #June
            series_summer[1::3] = series[6::12] #July
            series_summer[2::3] = series[7::12] #August
            # obtain winter time series
            series_winter[2::3] = series[11::12]#December
            series_winter[0::3] = series[0::12] #January
            series_winter[1::3] = series[1::12] #February
        elif series.ndim == 2:
            t, y = series.shape
            series_summer = np.zeros((t/4,y),dtype=float)
            series_winter = np.zeros((t/4,y),dtype=float)
              # obtain summer time series
            series_summer[0::3,:] = series[5::12,:] #June
            series_summer[1::3,:] = series[6::12,:] #July
            series_summer[2::3,:] = series[7::12,:] #August
            # obtain winter time series
            series_winter[2::3,:] = series[11::12,:]#December
            series_winter[0::3,:] = series[0::12,:] #January
            series_winter[1::3,:] = series[1::12,:] #February              
        elif series.ndim == 3:
            t, y, x = series.shape
            series_summer = np.zeros((t/4,y,x),dtype=float)
            series_winter = np.zeros((t/4,y,x),dtype=float)
              # obtain summer time series
            series_summer[0::3,:,:] = series[5::12,:,:] #June
            series_summer[1::3,:,:] = series[6::12,:,:] #July
            series_summer[2::3,:,:] = series[7::12,:,:] #August
            # obtain winter time series
            series_winter[2::3,:,:] = series[11::12,:,:]#December
            series_winter[0::3,:,:] = series[0::12,:,:] #January
            series_winter[1::3,:,:] = series[1::12,:,:] #February
        else:
            raise IOError("This module can not work with any array with a \
                           dimension higher than 3!")            
        
        return series_summer, series_winter
        
    @staticmethod
    def interpolation(series, lat_nav, lat_tar, interp_kind='slinear',Dim_month=True):
        """
        Interpolate a time series onto certain latitudes for the coupling
        and comparison of atmosphere and ocean fields.
        It is recommended to interpolate time series of oceanic fields on
        the latitude of the certain atmospheric fields to avoid the data
        out of range issues.
        param series: input time series
        param lat_nav: original latitude for the input data
        param lat_tar: target latitude for interpolation
        param interp_kind: the methods for interpolation, it includse
        -linear
        -nearest
        -(spline)slinear (default) / quadratic / cubic
        param Dim_month: there are two modes for removing the seasonal cycling
        -True (default) input time series have month dimension [year,month,...]
        -False input time series have only 1 dimension for time
        """
        if series.ndim > 3:
            raise IOError("This module can not work with any array with a \
                           dimension higher than 3!")
        else:
            if Dim_month == True:
                t, m, y = series.shape
                interp_series = np.zeros((t, m, len(lat_tar)), dtype=float)
                for i in np.arange(t):
                    for j in np.arange(m):
                        # call the data attribute in case it has mask
                        ius = scipy.interpolate.interp1d(lat_nav.data, series[i,j,:], kind=interp_kind,
                                                         bounds_error=False, fill_value=0.0)
                        interp_series[i,j,:] = ius(lat_tar.data)
            else:
                t, y = series.shape
                interp_series = np.zeros((t, len(lat_tar)), dtype=float)
                for i in np.arange(t):
                    ius = scipy.interpolate.interp1d(lat_nav.data, series[i,:], kind=interp_kind,
                                                     bounds_error=False, fill_value=0.0)
                    interp_series[i,:] = ius(lat_tar.data)
            
        return interp_series
