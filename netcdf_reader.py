# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 22:37:52 2021

@author: starg
"""


import netCDF4
import numpy as np
import glob
import re
import os
Mmol = 28.97 #[g/mol] src: https://www.cs.mcgill.ca/~rwest/wikispeedia/wpcd/wp/e/Earth%2527s_atmosphere.htm

C_d=1 # really doubtful!!!
R = 8.31446261815324 #[J⋅K^−1⋅mol^−1]


aij_list=glob.glob(r"EarthRotationRate/AIJ/*.nc")
oij_list=glob.glob(r"EarthRotationRate/OIJ/*.nc")
oijl_list=glob.glob(r"EarthRotationRate/OIJL/*.nc")

#periods=list(filter(reg.search,aij_list)) 
periods= [1,4,8,16,32,64,128,243,256]
for p in periods:
    reg=re.compile('_'+str(p)).search
    print(p,list(filter(reg,aij_list))[0],list(filter(reg,oij_list))[0],list(filter(reg,oijl_list))[0])
    atmo = netCDF4.Dataset(list(filter(reg,aij_list))[0])
    ocean = netCDF4.Dataset(list(filter(reg,oij_list))[0])
    oijl = netCDF4.Dataset(list(filter(reg,oijl_list))[0])
    atmo_vars=atmo.variables
    oc_vars=ocean.variables
    oijl_vars = oijl.variables
    lon=atmo_vars['lon']
    lat=atmo_vars['lat']
    wind_u=atmo_vars['usurf']
    wind_v=atmo_vars['vsurf']
    wind_vel=np.sqrt(np.power(wind_u[:],2)+np.power(wind_v[:],2))

    swnow_cov_per=atmo_vars['snowicefr']
    temp_surf=atmo_vars['tsurf']
    p_surf=atmo_vars['prsurf']

    

    oc_mixl=oc_vars['oij_mld']
    rho_surf = p_surf[:]/(temp_surf[:]*R/Mmol)
    tau=C_d*rho_surf*np.power(wind_vel[:],2)
'''
atmo = netCDF4.Dataset('EarthRotationRate/AIJ/ANN0940-0949.aijP211eoDOFP3Od_4daysg.nc')
ocean=netCDF4.Dataset('EarthRotationRate/OIJ/ANN0940-0949.oijP211eoDOFP3Od_4daysg.nc')
oijl=netCDF4.Dataset('EarthRotationRate/OIJL/ANN0940-0949.oijlP211eoDOFP3Od_4daysg.nc')


atmo256 = netCDF4.Dataset('EarthRotationRate/AIJ/ANN0940-0949.aijP211eoDOFP3Od_4daysg.nc')
ocean256=netCDF4.Dataset('EarthRotationRate/OIJ/ANN0940-0949.oijP211eoDOFP3Od_4daysg.nc')
oijl_256=netCDF4.Dataset('EarthRotationRate/OIJL/ANN0750-0799.oijlP211eoDOFP3Od_256daysg.nc')
atmo4_vars=atmo.variables
keys=atmo4_vars.keys()
#print(keys)
'''
'''

Hi Mike! 

That sounds like a fun project! Happy to help. 

The mixed layer depth (oij_mld) 
was unmodified and should be straightforward. 
Sea ice was only slightly modified. 
I just calculated sea ice cover with snowicefr (aij) 
after applying a land mask. 
I suppose oicefr (also aij) would be easier if you don’t mind some lakes. 

The others panels in Figure 7 required some simple calculations that use: 
    surface windspeed (wsurf, aij), potential density (pot_dens, oijl), 
    and downward vertical velocity (w, oijl). 
The equations are provided in Section 2.4 and for the calculations are fairly straightforward, 
the possible exception being upwelling. 

For the upwelling calculation, I wrote a silly script that looped through the i,j grid and found the ocean level that contained the base of the mixed layer using values from oij_mld. Then the script checked the downward vertical velocity (w in oijl). If w(i,j,l) was positive (downward) the script moved on; if w(i,j,l) was negative (upward), then w(i,j,l) was multiplied by the area of the cell and added to a running sum. This sum was ultimately normalized to sum of baseline Earth scenario and expressed as a positive number. You could normalize to one of your own simulations, doesn’t really matter. 

An easier way to do this would be just to do an area-weighted sum of the upward velocities at a fixed depth (e.g., 200 m). We got a qualitatively similar result doing it that way. (Obviously must still filter only upward velocities, otherwise the area-weighted sum should be zero...). 


'''
'''

lon=atmo4_vars['lon']
lat=atmo4_vars['lat']
wind_u=atmo4_vars['usurf']
wind_v=atmo4_vars['vsurf']
wind_vel=np.sqrt(np.power(wind_u[:],2)+np.power(wind_v[:],2))

swnow_cov_per=atmo4_vars['snowicefr']
temp_surf=atmo4_vars['tsurf']
p_surf=atmo4_vars['prsurf']

oc4_vars=ocean.variables

oc_mixl=oc4_vars['oij_mld']
rho_surf = p_surf[:]/(temp_surf[:]*R/Mmol)
tau=C_d*rho_surf*np.power(wind_vel[:],2)
'''
