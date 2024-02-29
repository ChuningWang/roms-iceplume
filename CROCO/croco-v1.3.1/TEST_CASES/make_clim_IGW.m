%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Build a ROMS climatology file
%
%  Extrapole and interpole temperature and salinity from a
%  Climatology to get boundary and initial conditions for
%  ROMS (climatology and initial netcdf files) .
%  Get the velocities and sea surface elevation via a 
%  geostrophic computation.
%
%  Data input format (netcdf):
%     temperature(T, Z, Y, X)
%     T : time [Months]
%     Z : Depth [m]
%     Y : Latitude [degree north]
%     X : Longitude [degree east]
%
%  Data source : IRI/LDEO Climate Data Library (World Ocean Atlas 1998)
%    http://ingrid.ldgo.columbia.edu/
%    http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NODC/.WOA98/
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
%  Copyright (c) 2002-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%  Contributions of P. Marchesiello (IRD)
%
%  Updated    1-Sep-2006 by Pierrick Penven
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
% Common parameters
%
romstools_param
%
%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%

%
% Title
%
disp(' ')
disp([' Making the clim: ',clmname])
disp(' ')
disp([' Title: ',ROMS_title])

%
% Read in the grid
%
disp(' ')
disp(' Read in the grid...')
nc=netcdf(grdname,'r');
Lp=length(nc('xi_rho'));
Mp=length(nc('eta_rho'));
hmax=max(max(nc{'h'}(:)));
close(nc);

%----------------------------------------------------------------------------
% Create the climatology file
%----------------------------------------------------------------------------

if (makeclim)
  disp(' ')
  disp(' Create the climatology file...')
  if  ~exist('vtransform')
      vtransform=1; %Old Vtransform
      disp([' NO VTRANSFORM parameter found'])
      disp([' USE VTRANSFORM default value vtransform = 1'])
  end
  create_climfile(clmname,grdname,ROMS_title,...
                  theta_s,theta_b,hc,N,...
                  woa_time,woa_cycle,'clobber',vtransform);
end

if (makeclim)
  nc=netcdf(clmname,'write');
  nc{'temp'}(:)=zeros(N,Mp,Lp);
  nc{'salt'}(:)=zeros(N,Mp,Lp);
  nc{'u'}(:)=zeros(N,Mp,Lp-1);
  nc{'v'}(:)=zeros(N,Mp-1,Lp);
  nc{'ubar'}(:)=zeros(Mp,Lp-1);
  nc{'vbar'}(:)=zeros(Mp-1,Lp);
  nc{'zeta'}(:)=zeros(Mp,Lp);
  nc{'SSH'}(:)=zeros(Mp,Lp);
  close(nc)
end

%----------------------------------------------------------------------------
% Initial file
%----------------------------------------------------------------------------
if (makeini)
  create_inifile(ininame,grdname,ROMS_title,...
                 theta_s,theta_b,hc,N,...
                 tini,'clobber',vtransform);
  nc=netcdf(ininame,'write');
  nc{'temp'}(:)=zeros(N,Mp,Lp);
  nc{'salt'}(:)=zeros(N,Mp,Lp);
  nc{'u'}(:)=zeros(N,Mp,Lp-1);
  nc{'v'}(:)=zeros(N,Mp-1,Lp);
  nc{'ubar'}(:)=zeros(Mp,Lp-1);
  nc{'vbar'}(:)=zeros(Mp-1,Lp);
  nc{'zeta'}(:)=zeros(Mp,Lp);
  close(nc)
end


%
% End
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
