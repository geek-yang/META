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
import os
import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.ticker as mticker
import iris
import iris.plot as iplt
import cartopy
import cartopy.crs as ccrs

class plots:
    @staticmethod
    def linearRegress(xaxis, corr, figname='./LinearRegression.png'):
        """
        This module will make a x-y plot to display the correlation coefficient
        got from the linear regression.
        
        param xaxis: latitude for the plot as x axis
        param corr: the correlation coefficient
        param figname: name and output path of figure
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
        plt.close(fig)
        
    @staticmethod    
    def leadlagRegress(yaxis, corr, lag, figname='./LeadLagRegression.png'):
        """
        This module will make a contour plot to display the correlation coefficient
        got from the lead/lag regression.
        
        param yaxis: latitude for the plot as y axis
        param corr: the correlation coefficient
        param lag: the maximum lag time for plot as x axis
        param figname: name and output path of figure
        return: Figures
        rtype: png
        """
        print ("Create contour plot of correlation coefficient.")
        # calculate the lead/lag index as x axis
        lag_index = np.arange(-lag,lag+1,1)
        xaxis = lag_index / 12
        # make plots
        fig = plt.figure()
        #contour_level = np.array([-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6])
        contour_level = np.array([-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8])
        cs = plt.contour(xaxis, yaxis, corr.transpose(),
                         contour_level, colors='k')
        plt.clabel(cs, inline=1, fontsize=10)
        plt.xlabel("Time Lag (year)")
        #lead_year = ['-15','-12','-9','-6','-3','0','3','6','9','12','15']
        plt.ylabel("Latitude")
        plt.show()
        fig.savefig(figname,dpi=400)
        plt.close(fig)
    
    @staticmethod
    def geograph(latitude, longitude, field, p_value, label, ticks,
                 figname='./NorthPolar.png', gridtype='geographical'):
        """
        This module will make a geographical plot to give a spatial view of fields.
        This module is built on iris and cartopy for the visualization of fields on
        both geographical and curvilinear grid.
        
        param lat: latitude coordinate for plot
        param lon: longitude coordinate for plot
        param field: input field for visualization
        param p_value: the significance level from t test. The significance level is set to be 99.5%.
        param gridtype: type of input spatial fields, it has two options
        - geographical (default) the coordinate is geographical, normally applied to atmosphere reanalysis
        - curvilinear the coordinate is curvilinear, normally applied to ocean reanalysis
        param figname: name and output path of figure
        return: figures
        rtype: png
        """
        print ("Create a NorthPolarStereo view of input fields.")
        if gridtype == 'geographical':
            print ("The input fields are originally on geographical grid")
            # first construct iris coordinate
            lat_iris = iris.coords.DimCoord(latitude, standard_name='latitude', long_name='latitude',
                                            var_name='lat', units='degrees')
            lon_iris = iris.coords.DimCoord(longitude, standard_name='longitude', long_name='longitude',
                                            var_name='lon', units='degrees')
            # assembly the cube
            cube_iris = iris.cube.Cube(field, long_name='geographical field', var_name='field', 
                                       units='1', dim_coords_and_dims=[(lat_iris, 0), (lon_iris, 1)])
            fig = plt.figure()
            ax = plt.axes(projection=ccrs.NorthPolarStereo())
            ax.set_extent([-180,180,20,90],ccrs.PlateCarree())
            ax.set_aspect('1')
            ax.coastlines()
            gl = ax.gridlines(linewidth=1, color='gray', alpha=0.5, linestyle='--')
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)
            ax.set_boundary(circle, transform=ax.transAxes)
            cs = iplt.pcolormesh(cube_iris, cmap='coolwarm', vmin=ticks[0], vmax=ticks[-1])
            cbar = fig.colorbar(cs,extend='both', orientation='horizontal',
                                 shrink =0.8, pad=0.05, format="%.1f")
            cbar.set_label(label,size = 8)
            cbar.set_ticks(ticks)
            cbar.ax.tick_params(labelsize = 6)
            ii, jj = np.where(p_value<=0.005)
            #ax.scatter(longitude[jj], latitude[ii], transform=ccrs.Geodetic(),
            #           s=0.1, c='g',alpha=0.3)
            iplt.show()
            fig.savefig(figname, dpi=300)
            plt.close(fig)
        elif gridtype == 'curvilinear':
            print ("The input fields are originally on curvilinear grid")
            # first construct iris coordinate
            lat_iris = iris.coords.AuxCoord(latitude, standard_name='latitude', long_name='latitude',
                                            var_name='lat', units='degrees')
            lon_iris = iris.coords.AuxCoord(longitude, standard_name='longitude', long_name='longitude',
                                            var_name='lon', units='degrees')
            # assembly the cube
            cube_iris = iris.cube.Cube(field, long_name='curvilinear field', var_name='field', 
                                       units='1', aux_coords_and_dims=[(lat_iris, (0,1)), (lon_iris, (0,1))])
            coord_sys = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
            cube_iris.coord('latitude').coord_system = coord_sys
            cube_iris.coord('longitude').coord_system = coord_sys
            # determine nx and ny for interpolation
            jj, ii = latitude.shape
            if ii > 1000:
                nx = 1440
                ny = 350
            else:
                nx = 720
                ny = 140
            cube_regrid, extent = iris.analysis.cartography.project(cube_iris, ccrs.PlateCarree(), nx, ny)
            # make plots
            fig = plt.figure()
            ax = plt.axes(projection=ccrs.NorthPolarStereo())
            ax.set_extent([-180,180,20,90],ccrs.PlateCarree())
            ax.set_aspect('1')
            ax.coastlines()
            gl = ax.gridlines(linewidth=1, color='gray', alpha=0.5, linestyle='--')
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)
            ax.set_boundary(circle, transform=ax.transAxes)
            cs = iplt.pcolormesh(cube_regrid, cmap='coolwarm', vmin=ticks[0], vmax=ticks[-1])
            cbar = fig.colorbar(cs,extend='both', orientation='horizontal',
                                shrink =0.8, pad=0.05, format="%.1f")
            cbar.set_label(label,size = 8)
            cbar.set_ticks(ticks)
            cbar.ax.tick_params(labelsize = 6)
            ii, jj = np.where(p_value<=0.005)
            #ax.scatter(longitude[jj], latitude[ii], transform=ccrs.Geodetic(),
            #           s=0.1, c='g',alpha=0.3)
            iplt.show()
            fig.savefig(figname, dpi=300)
            plt.close(fig)            
        else:
            raise IOError("This module only support fields on geographical or curvilinear grid!")