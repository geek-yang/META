# -*- coding: utf-8 -*-
"""
Copyright Netherlands eScience Center
Function        : Plots generator for visualization
Author          : Yang Liu (y.liu@esciencecenter.nl)
First Built     : 2018.08.13
Last Update     : 2018.08.13
Contributor     :
Description     : This module provides several methods to perform statistical
                  analysis on MET and all kinds of fields.
Return Values   : pngs
Caveat!         :
"""

import numpy as np
import scipy
#from scipy import stats
import cartopy
import os
import iris
import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt

class plots:
    @staticmethod
    def linearRegress(xaxis, corr, figname='./LinearRegression.png'):
        """
        This module will make a x-y plot to display the correlation coefficient
        got from the linear regression.
        
        param 
        return: Figures
        rtype: png
        """
        print ("Create x-y plot of correlation coefficient.")
        fig = plt.figure()
        plt.plot(xaxis, corr)
        plt.xlabel("Latitude")
        #plt.xticks(np.linspace(20, 90, 11))
        plt.ylabel("Correlation Coefficient")
        plt.show()
        fig.savefig(figname,dpi=400)
        
    @staticmethod    
    def leadlagRegress(yaxis, corr, lag, figname='./LeadLagRegression.png'):
        """
        This module will make a contour plot to display the correlation coefficient
        got from the lead/lag regression.
        
        param 
        return: Figures
        rtype: png
        """
        print ("Create contour plot of correlation coefficient.")
        # calculate the lead/lag index as x axis
        lag_index = np.arange(-lag,lag+1,1)
        xaxis = lag_index / 12
        # make plots
        fig = plt.figure()
        contour_level = np.array([-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6])
        #contour_level = np.array([-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8])
        cs = plt.contour(xaxis, yaxis, corr.transpose(),
                         contour_level, colors='k')
        plt.clabel(cs, inline=1, fontsize=10)
        plt.xlabel("Time Lag (year)")
        #lead_year = ['-15','-12','-9','-6','-3','0','3','6','9','12','15']
        plt.ylabel("Latitude")
        plt.show()
        fig.savefig(figname,dpi=400)