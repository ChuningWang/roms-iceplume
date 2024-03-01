#!/usr/bin/env python
# coding: utf-8
import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from netCDF4 import Dataset
from matplotlib.font_manager import FontProperties
import sys

#myfile='benguela-008x008_052'
#myncfile='croco_grd.nc'

if len(sys.argv) != 3:
	print ("Usage : mpp_plot.py GRID_FILE COVDTA_FILE")
	sys.exit (1)

myncfile=sys.argv[1]
myfile=sys.argv[2]

ncid   = Dataset(myncfile)
lat2d  = ncid.variables['lat_rho'    ][  :,:].squeeze()
lon2d  = ncid.variables['lon_rho'    ][  :,:].squeeze()
msk    = ncid.variables['mask_rho'][:,:].squeeze()
bathy  = ncid.variables['h'][:,:].squeeze()
ncid.close()

xmax=np.size(lat2d,1)-1
ymax=np.size(lat2d,0)-1

bathy=np.ma.masked_equal(bathy*msk,0. )

fig, ax = plt.subplots()

plt.contour(lon2d,lat2d,bathy,colors='k',linestyles='dashed',linewidths=0.5)
vlevel=np.arange(0,6000,100)

plt.contourf(lon2d,lat2d,bathy,levels=vlevel)

plt.title('Bathymetry (m)')
plt.ylabel('Latitude',fontsize=14)
plt.xlabel('Longitude',fontsize=14)
plt.colorbar()

NPTS=1
font0 = FontProperties()
alignment = {'horizontalalignment': 'center', 'verticalalignment': 'baseline','color': 'r'}
# Show family options
families = ['serif', 'sans-serif', 'cursive', 'fantasy', 'monospace']
font1 = font0.copy()
font1.set_weight('bold')


Path = mpath.Path
file = open(myfile,"r")
L0=file.readline()
L0=int(L0.replace('#',''))
for i in range(L0):
	L00=file.readline().split()

	L1=file.readline().split()
	LON1=lon2d[max(int(L1[1])-NPTS,0),max(int(L1[0])-NPTS,0)]
	LAT1=lat2d[max(int(L1[1])-NPTS,0),max(int(L1[0])-NPTS,0)]    

	L2=file.readline().split()
	LON2=lon2d[min(int(L2[1])+NPTS,ymax),min(int(L2[0])+NPTS,xmax)]
	LAT2=lat2d[max(int(L2[1])-NPTS,0),max(int(L2[0])-NPTS,0)]    

	L3=file.readline().split()
	LON3=lon2d[min(int(L3[1])+NPTS,ymax),min(int(L3[0])+NPTS,xmax)]
	LAT3=lat2d[min(int(L3[1])+NPTS,ymax),min(int(L3[0])+NPTS,xmax)]    

	L4=file.readline().split()
	LON4=lon2d[max(int(L4[1])-NPTS,0),max(int(L4[0])-NPTS,0)]
	LAT4=lat2d[min(int(L4[1])+NPTS,ymax),min(int(L4[0])+NPTS,xmax)]    

	L5=file.readline().split()
	LON5=lon2d[max(int(L5[1])-NPTS,0),max(int(L5[0])-NPTS,0)]
	LAT5=lat2d[max(int(L5[1])-NPTS,0),max(int(L5[0])-NPTS,0)]    

	if i==0:
		toto=(LAT3-LAT2)/(lat2d[-1,0]-lat2d[0,0])
		font1.set_size(70*toto)

	L000=file.readline()

	path_data = [
    (Path.MOVETO, (LON1,LAT1)),
    (Path.LINETO, (LON2,LAT2)),
    (Path.LINETO, (LON3,LAT3)),
    (Path.LINETO, (LON4,LAT4)),
    (Path.CLOSEPOLY, (LON5,LAT5))]

	codes, verts = zip(*path_data)
	path = mpath.Path(verts, codes)
	patch = mpatches.PathPatch(path, facecolor='b', alpha=0.2)
	ax.add_patch(patch)		
	ax.add_patch(patch)

	mynum=int(L00[1]) -1
	t = plt.text((LON1+LON2)/2., (LAT1+LAT3)/2.,str(mynum), fontproperties=font1,
             **alignment)	

iend=i

L0=file.readline()
L0=int(L0.replace('# vides:',''))
for i in range(L0):
	L00=file.readline()
	L1=file.readline().split()
	LON1=lon2d[max(int(L1[1])-NPTS,0),max(int(L1[0])-NPTS,0)]
	LAT1=lat2d[max(int(L1[1])-NPTS,0),max(int(L1[0])-NPTS,0)]    

	L2=file.readline().split()
	LON2=lon2d[min(int(L2[1])+NPTS,ymax),min(int(L2[0])+NPTS,xmax)]
	LAT2=lat2d[max(int(L2[1])-NPTS,0),max(int(L2[0])-NPTS,0)]    

	L3=file.readline().split()
	LON3=lon2d[min(int(L3[1])+NPTS,ymax),min(int(L3[0])+NPTS,xmax)]
	LAT3=lat2d[min(int(L3[1])+NPTS,ymax),min(int(L3[0])+NPTS,xmax)]    

	L4=file.readline().split()
	LON4=lon2d[max(int(L4[1])-NPTS,0),max(int(L4[0])-NPTS,0)]
	LAT4=lat2d[min(int(L4[1])+NPTS,ymax),min(int(L4[0])+NPTS,xmax)]    

	L5=file.readline().split()
	LON5=lon2d[max(int(L5[1])-NPTS,0),max(int(L5[0])-NPTS,0)]
	LAT5=lat2d[max(int(L5[1])-NPTS,0),max(int(L5[0])-NPTS,0)]    


	L000=file.readline()
	L000=file.readline()
	L000=file.readline()
	L000=file.readline()

	path_data = [
    (Path.MOVETO, (LON1,LAT1)),
    (Path.LINETO, (LON2,LAT2)),
    (Path.LINETO, (LON3,LAT3)),
    (Path.LINETO, (LON4,LAT4)),
    (Path.CLOSEPOLY, (LON5,LAT5))]

	codes, verts = zip(*path_data)
	path = mpath.Path(verts, codes)
	patch = mpatches.PathPatch(path, facecolor='r', alpha=0.2)
	ax.add_patch(patch)	
#	t = plt.text((LON1+LON2)/2., (LAT1+LAT3)/2., str(i+iend), fontproperties=font1,
#             **alignment)	

file.close()

plt.show()
