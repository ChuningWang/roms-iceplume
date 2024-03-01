%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Add the tides to the ROMS forcing file for tidal forcing
%  by the boundary conditions.
%
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2003-2006 by Patrick Marchesiello and Meinte Blass
%
%  Updated   1-Sep-2006 by Pierrick Penven (generalisation of romstools_param.m)
%  Updated   3-Oct-2006 by Pierrick Penven (cleaning + phase lag for Yorig time)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
% Common parameters
%
romstools_param
%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%
%
% Get start time of simulation in fractional mjd for nodal correction
%
date_mjd=0.; %mjd(Ymin,Mmin,Dmin);
deg=180.0/pi;
%
% Add a phase correction to be consistent with the 'Yorig' time
%
t0=0.; %t0-24*mjd(Yorig,1,1);
%
%  Read in ROMS grid.
%
disp('Reading ROMS grid parameters ...');
nc=netcdf(grdname,'r');
lonr=nc{'lon_rho'}(:);
latr=nc{'lat_rho'}(:);
lonu=nc{'lon_u'}(:);
latu=nc{'lat_u'}(:);
lonv=nc{'lon_v'}(:);
latv=nc{'lat_v'}(:);
rangle=nc{'angle'}(:); % devrait etre utilise....
h=nc{'h'}(:);
rmask=nc{'mask_rho'}(:);
close(nc)
%
% Create the forcing file and add tides
%
disp(' ')
disp(' Create the forcing file...')
create_forcing(frcname,grdname,ROMS_title,...
               coads_time,coads_time,coads_time,...
               coads_time,coads_time,coads_time,...
               coads_cycle,coads_cycle,coads_cycle,...
               coads_cycle,coads_cycle,coads_cycle)
components='S2';
disp(['Tidal components : ',components])
nc_add_tides(frcname,Ntides,date_mjd,components)
ncfrc=netcdf(frcname,'write');
%
% Loop on periods
%
for itide=1:Ntides;
  ncfrc{'tide_period'}(itide)=12.;
%
% Process surface elevation
%
  disp('  ssh...')
  fname='IGW_FILES/amp_ssh_S2.cdf';
  nc=netcdf(fname);
  ssh_amp=squeeze(nc{'amp_ssh_S2'}(1,3:7,:));
  close(nc)
  fname='IGW_FILES/pha_ssh_S2.cdf';
  nc=netcdf(fname);
  ssh_pha=squeeze(nc{'pha_ssh_S2'}(1,3:7,:));
  close(nc)
  ncfrc{'tide_Ephase'}(itide,:,:)=ssh_pha;     
  ncfrc{'tide_Eamp'}(itide,:,:)=ssh_amp;
%
% Process U
%
  disp('  u...')
  fname='IGW_FILES/amp_u_S2.cdf';
  nc=netcdf(fname);
  uamp=squeeze(nc{'amp_u_S2'}(1,3:7,2:end));
  close(nc)
  fname='IGW_FILES/pha_u_S2.cdf';
  nc=netcdf(fname);
  upha=squeeze(nc{'pha_u_S2'}(1,3:7,2:end));;
  close(nc)
  uamp=u2rho_2d(uamp);
  upha=u2rho_2d(upha);
%
% Process V
%
  disp('  v...')
  fname='IGW_FILES/amp_v_S2.cdf';
  nc=netcdf(fname);
  vamp=squeeze(nc{'amp_v_S2'}(1,4:7,:));
  close(nc)
  fname='IGW_FILES/pha_v_S2.cdf';
  nc=netcdf(fname);
  vpha=squeeze(nc{'pha_v_S2'}(1,4:7,:));;
  close(nc)
  vamp=v2rho_2d(vamp);
  vpha=v2rho_2d(vpha);
%
% Convert to tidal ellipses
%
  disp('  Convert to tidal ellipse parameters...')
  [major,eccentricity,inclination,phase]=ap2ep(uamp,upha,vamp,vpha);
  ncfrc{'tide_Cmin'}(itide,:,:)=major.*eccentricity;
  ncfrc{'tide_Cmax'}(itide,:,:)=major;;
  ncfrc{'tide_Cangle'}(itide,:,:)=inclination;
  ncfrc{'tide_Cphase'}(itide,:,:)=phase;
%
end
%
% Close the file
%
close(ncfrc)
%
% Plot
%
if makeplot==1
  disp('  Make a few plots')
  figure('position',[300 300 700 300])
  plot_tide(grdname,frcname,1,0.1,100,coastfileplot)
end
