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
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

Mmol = 28.97 #[g/mol] src: https://www.cs.mcgill.ca/~rwest/wikispeedia/wpcd/wp/e/Earth%2527s_atmosphere.htm

C_d=1 # really doubtful!!!
R = 8.31446261815324 #[J⋅K^−1⋅mol^−1]



aij_list=glob.glob(r"EarthRotationRate/AIJ/*.nc")
oij_list=glob.glob(r"EarthRotationRate/OIJ/*.nc")
oijl_list=glob.glob(r"EarthRotationRate/OIJL/*.nc")

atmo_1 = netCDF4.Dataset(list(filter(re.compile('_'+str(1)).search,aij_list))[0])
ocean_1 = netCDF4.Dataset(list(filter(re.compile('_'+str(1)).search,oij_list))[0])
oijl_1 = netCDF4.Dataset(list(filter(re.compile('_'+str(1)).search,oijl_list))[0])

atmo_vars_1=atmo_1.variables
oc_vars_1=ocean_1.variables
oijl_vars_1 = oijl_1.variables

lons=atmo_vars_1['lon'][:]
lats=atmo_vars_1['lat'][:]



downward_v_1=oijl_vars_1['w'][:]
oceanfr_1=atmo_vars_1['ocnfr'] #  ocean fraction
ocfrac_land_1=oc_vars_1['oxyp'] 
cellarea_1=ocfrac_land_1[:]/oceanfr_1[:] # grid cell area
downward_v_1[downward_v_1>0]=0
upwelling_1=cellarea_1*downward_v_1[0]
upwell_sum_1=np.sum(upwelling_1)

tau_1=oc_vars_1['oij_ustar']
oc_mixl_1=np.mean(oc_vars_1['oij_mld'][:])

land_1=atmo_vars_1['frac_land']
snow_cov_per_1=atmo_vars_1['snowicefr']
ice_1=snow_cov_per_1[:]-land_1[:]
ice_1[ice_1<0]=0
icecov_fr_1=np.mean(ice_1)

fig = plt.figure(600)

ax = fig.add_subplot(941)

ax.set_title("Ice coverage")
m = Basemap(projection='moll',lat_0=0,lon_0=0)

lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)
m.drawcoastlines()
cs = m.pcolor(xi,yi,ice_1[:],cmap='jet')
cbar = m.colorbar(cs, location='bottom', pad="10%")


ax = fig.add_subplot(942)
ax.set_title("Upwelling")
m = Basemap(projection='moll',lat_0=0,lon_0=0)
m.drawcoastlines()
lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.pcolor(xi,yi,np.squeeze(upwelling_1[:]),cmap='jet')
cbar = m.colorbar(cs, location='bottom', pad="10%")

ax = fig.add_subplot(943)
ax.set_title("Wind stress")
m = Basemap(projection='moll',lat_0=0,lon_0=0)
m.drawcoastlines()
lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.pcolor(xi,yi,np.squeeze(tau_1[:]),cmap='jet')
cbar = m.colorbar(cs, location='bottom', pad="10%")

ax = fig.add_subplot(944)
ax.set_title("Mixing layer depth")
m = Basemap(projection='moll',lat_0=0,lon_0=0)
m.drawcoastlines()
lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.pcolor(xi,yi,oc_vars_1['oij_mld'][:],cmap='jet')
cbar = m.colorbar(cs, location='bottom', pad="10%")


plt.show()


upwell_arr=[1]
tau_arr=[1]
mld=[oc_mixl_1]
icecov=[icecov_fr_1]

#periods=list(filter(reg.search,aij_list)) 
periods= [4,8,16,32,64,128,243,256]
i=4
for p in periods:
    reg=re.compile('_'+str(p)).search #Regular Expression looking for string matching period duration in .nc filenames
    atmo = netCDF4.Dataset(list(filter(reg,aij_list))[0]) #read AIJ .nc file with given period
    ocean = netCDF4.Dataset(list(filter(reg,oij_list))[0]) #read OIJ .nc file with given period
    oijl = netCDF4.Dataset(list(filter(reg,oijl_list))[0]) #read OIJL .nc file with given period
    atmo_vars=atmo.variables # dict of AIJ table variables
    oc_vars=ocean.variables # dict of OIJ table variables
    oijl_vars = oijl.variables# dict of OIJL table variables
    oclayers=np.linspace(0,5000,13) # ocean layers - obsolete thus far
    
   
    wind_u=atmo_vars['usurf'] #u component of surface wind velocity
    wind_v=atmo_vars['vsurf'] #v component of surface wind velocity
    wind_vel=np.sqrt(np.power(wind_u[:],2)+np.power(wind_v[:],2)) # wind velocity as sqrt(u^2 + v^2)

    snow_cov_per=atmo_vars['snowicefr'] #snow ince fraction
    land=atmo_vars['frac_land'] #land fraction
    
    temp_surf=atmo_vars['tsurf'] #surface temperature - obsolete
    p_surf=atmo_vars['prsurf'] #surface pressure - obsolete
    ice=snow_cov_per[:]-land[:] # masking out land from snow coverage map
    ice[ice<0]=0
    icecov_fr=np.mean(ice) # global average ice coverage fraction
    icecov.append(icecov_fr)

    oc_mixl=oc_vars['oij_mld'][:] # ocean mixing layer depth [m]
    # rho_surf = p_surf[:]/(temp_surf[:]*R/Mmol)
    # tau=C_d*rho_surf*np.power(wind_vel[:],2) -- obsolete substitution of ideal gas law into wind stress formula
    tau=oc_vars['oij_ustar'] # windstress
    
    pot_dens=oijl_vars['pot_dens']
    downward_v=oijl_vars['w'][:]
    oceanfr=atmo_vars['ocnfr'] #  ocean fraction
    ocfrac_land=oc_vars['oxyp'] 
    cellarea=ocfrac_land[:]/oceanfr[:] # grid cell area
    downward_v[downward_v>0]=0
    upwelling=cellarea*downward_v[0]
    upwell_sum=np.sum(upwelling)/upwell_sum_1
    upwell_arr.append(upwell_sum)
    tau_arr.append(np.mean(tau)/np.mean(tau_1))
    oc_mixl=np.mean(oc_vars['oij_mld'][:])
    mld.append(oc_mixl)
    i=i+1
    ax = fig.add_subplot(9,4,4+i)
    

