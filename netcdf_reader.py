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

fig, axs = plt.subplots(11,4, figsize=(16, 18))
#fig.subplots_adjust(hspace = .25, wspace=.1)
#ax = fig.add_subplot(941)
axs = axs.ravel()


pad = 5 # in points


axs[0].annotate("P = {} d".format(str(1)),xy=(0, 0.5), xytext=(-axs[0].yaxis.labelpad-pad,0), xycoords=axs[0].yaxis.label, textcoords='offset points',size='large', ha='right', va='center')
axs[0].set_title("Ice coverage")
m = Basemap(projection='moll',lat_0=0,lon_0=0,ax=axs[0])

lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)
m.drawcoastlines()
m.fillcontinents(color='white')
cs = m.pcolor(xi,yi,ice_1[:],cmap='jet')
cbar = m.colorbar(cs, location='bottom',size="5%", pad='2%')


#axs[2] = fig.add_subplot(942)
axs[1].set_title("Upwelling")
m = Basemap(projection='moll',lat_0=0,lon_0=0,ax=axs[1])
m.drawcoastlines()
m.fillcontinents(color='white')
lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.pcolor(xi,yi,np.squeeze(upwelling_1[:]),cmap='jet')
cbar = m.colorbar(cs, location='bottom',size="5%", pad='2%')

#ax = fig.add_subplot(943)
axs[2].set_title("Wind stress")
m = Basemap(projection='moll',lat_0=0,lon_0=0,ax=axs[2])
m.drawcoastlines()
m.fillcontinents(color='white')
lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.pcolor(xi,yi,np.squeeze(tau_1[:]),cmap='jet')
cbar = m.colorbar(cs, location='bottom',size="5%", pad='2%')

#ax = fig.add_subplot(944)
axs[3].set_title("Mixing layer depth")
m = Basemap(projection='moll',lat_0=0,lon_0=0,ax=axs[3])
m.drawcoastlines()
m.fillcontinents(color='white')
lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.pcolor(xi,yi,oc_vars_1['oij_mld'][:],cmap='jet')
cbar = m.colorbar(cs, location='bottom',size="5%", pad='2%')


upwell_arr=[1]
tau_arr=[1]
mld=[oc_mixl_1]
icecov=[icecov_fr_1]

#periods=list(filter(reg.search,aij_list)) 
periods= [4,8,16,32,64,128,243,256]
col=3

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
    
    axs[col+1].annotate("P = {} d".format(str(p)),xy=(0, 0.5), xytext=(-axs[col].yaxis.labelpad-pad,0), xycoords=axs[col+1].yaxis.label, textcoords='offset points',size='large', ha='right', va='center')
    m = Basemap(projection='moll',lat_0=0,lon_0=0,ax=axs[col+1])
    lon, lat = np.meshgrid(lons, lats)
    xi, yi = m(lon, lat)
    m.drawcoastlines()
    m.fillcontinents(color='white')
    cs = m.pcolor(xi,yi,ice[:],cmap='jet')
    cbar = m.colorbar(cs, location='bottom',size="5%", pad='2%')
    

    m = Basemap(projection='moll',lat_0=0,lon_0=0,ax=axs[col+2])
    m.drawcoastlines()
    m.fillcontinents(color='white')
    lon, lat = np.meshgrid(lons, lats)
    xi, yi = m(lon, lat)

    cs = m.pcolor(xi,yi,np.squeeze(upwelling[:]),cmap='jet')
    cbar = m.colorbar(cs, location='bottom',size="5%", pad='2%')


    m = Basemap(projection='moll',lat_0=0,lon_0=0,ax=axs[col+3])
    m.drawcoastlines()
    m.fillcontinents(color='white')
    lon, lat = np.meshgrid(lons, lats)
    xi, yi = m(lon, lat)

    cs = m.pcolor(xi,yi,np.squeeze(tau[:]),cmap='jet')
    cbar = m.colorbar(cs, location='bottom',size="5%", pad='2%')

    
    m = Basemap(projection='moll',lat_0=0,lon_0=0,ax=axs[col+4])
    m.drawcoastlines()
    m.fillcontinents(color='white')
    lon, lat = np.meshgrid(lons, lats)
    xi, yi = m(lon, lat)

    cs = m.pcolor(xi,yi,oc_vars['oij_mld'][:],cmap='jet')
    cbar = m.colorbar(cs, location='bottom',size="5%", pad='2%')
    
    col=col+4
'''
tidal locking - P = 365 d
'''    
atmo_365 = netCDF4.Dataset(list(filter(re.compile('_'+str(365)).search,aij_list))[0])
ocean_365 = netCDF4.Dataset(list(filter(re.compile('_'+str(365)).search,oij_list))[0])
oijl_365 = netCDF4.Dataset(list(filter(re.compile('_'+str(365)).search,oijl_list))[0])

atmo_vars_365=atmo_365.variables
oc_vars_365=ocean_365.variables
oijl_vars_365 = oijl_365.variables

lons=atmo_vars_365['lon'][:]
lats=atmo_vars_365['lat'][:]



downward_v_365=oijl_vars_365['w'][:]
oceanfr_365=atmo_vars_365['ocnfr'] #  ocean fraction
ocfrac_land_365=oc_vars_365['oxyp'] 
cellarea_365=ocfrac_land_365[:]/oceanfr_365[:] # grid cell area
downward_v_365[downward_v_365>0]=0
upwelling_365=cellarea_365*downward_v_365[0]
upwell_sum_365=np.sum(upwelling_365)

tau_365=oc_vars_365['oij_ustar']
oc_mixl_365=np.mean(oc_vars_365['oij_mld'][:])

land_365=atmo_vars_365['frac_land']
snow_cov_per_365=atmo_vars_365['snowicefr']
ice_365=snow_cov_per_365[:]-land_365[:]
ice_365[ice_365<0]=0
icecov_fr_365=np.mean(ice_365)
icecov.append(icecov_fr_365)

upwell_sum=np.sum(upwelling_365)/upwell_sum_1
upwell_arr.append(upwell_sum)
tau_arr.append(np.mean(tau_365)/np.mean(tau_1))
oc_mixl=np.mean(oc_vars_365['oij_mld'][:])
mld.append(oc_mixl)



axs[col+1].annotate("P = {} d".format(str(365)),xy=(0, 0.5), xytext=(-axs[col+1].yaxis.labelpad-pad,0), xycoords=axs[col+1].yaxis.label, textcoords='offset points',size='large', ha='right', va='center')

m = Basemap(projection='moll',lat_0=0,lon_0=0,ax=axs[col+1])

lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)
m.drawcoastlines()
m.fillcontinents(color='white')
cs = m.pcolor(xi,yi,ice_365[:],cmap='jet')
cbar = m.colorbar(cs, location='bottom',size="5%", pad='2%')


m = Basemap(projection='moll',lat_0=0,lon_0=0,ax=axs[col+2])
m.drawcoastlines()
m.fillcontinents(color='white')
lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.pcolor(xi,yi,np.squeeze(upwelling_365[:]),cmap='jet')
cbar = m.colorbar(cs, location='bottom',size="5%", pad='2%')


m = Basemap(projection='moll',lat_0=0,lon_0=0,ax=axs[col+3])
m.drawcoastlines()
m.fillcontinents(color='white')
lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.pcolor(xi,yi,np.squeeze(tau_365[:]),cmap='jet')
cbar = m.colorbar(cs, location='bottom',size="5%", pad='2%')


m = Basemap(projection='moll',lat_0=0,lon_0=0,ax=axs[col+4])
m.drawcoastlines()
m.fillcontinents(color='white')
lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.pcolor(xi,yi,oc_vars_365['oij_mld'][:],cmap='jet')
cbar = m.colorbar(cs, location='bottom',size="5%", pad='2%')

'''
tidal locking - P = 365 d, 0 obliquity
'''    
atmo_365_0 = netCDF4.Dataset(list(filter(re.compile('_'+str(365)+'daysg0o').search,aij_list))[0])
ocean_365_0 = netCDF4.Dataset(list(filter(re.compile('_'+str(365)+'daysg0o').search,oij_list))[0])
oijl_365_0 = netCDF4.Dataset(list(filter(re.compile('_'+str(365)+'daysg0o').search,oijl_list))[0])

atmo_vars_365_0=atmo_365_0.variables
oc_vars_365_0=ocean_365_0.variables
oijl_vars_365_0 = oijl_365_0.variables

lons=atmo_vars_365_0['lon'][:]
lats=atmo_vars_365_0['lat'][:]



downward_v_365_0=oijl_vars_365_0['w'][:]
oceanfr_365_0=atmo_vars_365_0['ocnfr'] #  ocean fraction
ocfrac_land_365_0=oc_vars_365_0['oxyp'] 
cellarea_365_0=ocfrac_land_365_0[:]/oceanfr_365_0[:] # grid cell area
downward_v_365_0[downward_v_365_0>0]=0
upwelling_365_0=cellarea_365_0*downward_v_365_0[0]
upwell_sum_365_0=np.sum(upwelling_365_0)

tau_365_0=oc_vars_365_0['oij_ustar']
oc_mixl_365_0=np.mean(oc_vars_365_0['oij_mld'][:])

land_365_0=atmo_vars_365_0['frac_land']
snow_cov_per_365_0=atmo_vars_365_0['snowicefr']
ice_365_0=snow_cov_per_365_0[:]-land_365_0[:]
ice_365_0[ice_365_0<0]=0
icecov_fr_365_0=np.mean(ice_365_0)


icecov.append(icecov_fr_365_0)

upwell_sum=np.sum(upwelling_365_0)/upwell_sum_1
upwell_arr.append(upwell_sum)
tau_arr.append(np.mean(tau_365_0)/np.mean(tau_1))
oc_mixl=np.mean(oc_vars_365_0['oij_mld'][:])
mld.append(oc_mixl)



axs[col+5].annotate("P = {} d, $\epsilon=0$".format(str(365)),xy=(0, 0.5), xytext=(-axs[col+5].yaxis.labelpad-pad,0), xycoords=axs[col+5].yaxis.label, textcoords='offset points',size='large', ha='right', va='center')

m = Basemap(projection='moll',lat_0=0,lon_0=0,ax=axs[col+5])

lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)
m.drawcoastlines()
m.fillcontinents(color='white')
cs = m.pcolor(xi,yi,ice_365_0[:],cmap='jet')
cbar = m.colorbar(cs, location='bottom',size="5%", pad='2%')


m = Basemap(projection='moll',lat_0=0,lon_0=0,ax=axs[col+6])
m.drawcoastlines()
m.fillcontinents(color='white')
lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.pcolor(xi,yi,np.squeeze(upwelling_365_0[:]),cmap='jet')
cbar = m.colorbar(cs, location='bottom',size="5%", pad='2%')


m = Basemap(projection='moll',lat_0=0,lon_0=0,ax=axs[col+7])
m.drawcoastlines()
m.fillcontinents(color='white')
lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.pcolor(xi,yi,np.squeeze(tau_365_0[:]),cmap='jet')
cbar = m.colorbar(cs, location='bottom',size="5%", pad='2%')

m = Basemap(projection='moll',lat_0=0,lon_0=0,ax=axs[col+8])
m.drawcoastlines()
m.fillcontinents(color='white')
lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.pcolor(xi,yi,oc_vars_365_0['oij_mld'][:],cmap='jet')
cbar = m.colorbar(cs, location='bottom',size="5%", pad='2%')
#plt.savefig("ParameterGlobalMaps.png",dpi=1200)    
plt.show()

'''
mean values' scatterplots
'''
periods1=[1]+periods+[365]
fig, axs = plt.subplots(2,2, figsize=(10,8))
fig.subplots_adjust(hspace = .5, wspace=.25)
axs[0,0].plot(periods1,upwell_arr[:-1],c='mediumseagreen',marker='o',linestyle='dashed',linewidth=1, markersize=8)
axs[0,0].set_title("Mean upwelling")
axs[0,0].set_xlabel("Rotation Period [d]")



axs[0,1].plot(periods1,icecov[:-1],c='grey',marker='o',linestyle='dashed',linewidth=1, markersize=8)
axs[0,1].set_title("Mean ice coverage")
axs[0,1].set_xlabel("Rotation Period [d]")
axs[0,1].set_ylabel("Ice coverage $[\%]$")



axs[1,0].plot(periods1,tau_arr[:-1],c='aqua',marker='o',linestyle='dashed',linewidth=1, markersize=8)
axs[1,0].set_title("Mean wind stress")
axs[1,0].set_xlabel("Rotation Period [d]")
axs[1,0].set_ylabel(r"$\tau [g m^{-1} s^{-2}]$")


axs[1,1].plot(periods1,mld[:-1],c='navy',marker='o',linestyle='dashed',linewidth=1, markersize=8)
axs[1,1].set_title("Mean mixing layer depth")
axs[1,1].set_xlabel("Rotation Period [d]")
axs[1,1].set_ylabel("Mixing layer depth [m]")
plt.show()