%======================================================================
%
%     ---               TFLAT2DV Test Case                ---     
%
%
%  This file is part of CROCOTOOLS
%
%  CROCOTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  CROCOTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
% 
%  Patrick Marchesiello - 2012
%======================================================================
clear all
close all
%================== User defined parameters ===========================

fname='tidal_flat_his.nc';
makepdf   = 0;                     % make pdf file

%======================================================================
%
% Read data
%

idy=1;
idz=1;

nc=netcdf(fname);

h=squeeze(nc{'h'}(idy,:));
Dcrit=squeeze(nc{'Dcrit'}(idy,:));
X=squeeze(nc{'x_rho'}(idy,:))/1000;
T=squeeze(nc{'scrum_time'}(:,:))/86400;
sand=squeeze(nc{'SAND'}(:,idz,idy,:));
zeta=squeeze(nc{'zeta'}(:,idy,:));
ubar=u2rho_2d(squeeze(nc{'ubar'}(:,idy,:)));

close(nc);


T1=T(1);
T=T-T1;
D=zeta+repmat(h,size(zeta,1),1);
sand(D<Dcrit+0.01)=NaN;
zeta(D<Dcrit+0.01)=NaN;
ubar(D<Dcrit+0.01)=NaN;

%----------------------------------------------------------
%  Plot Hovmoller
%----------------------------------------------------------
figure('Position',[1 1 800 800])

subplot(2,1,1)
pcolor(X,T,ubar), shading flat
caxis([-3 3]);
%cbr=colorbar('fontsize',13);
c=colorbar();
set(c, 'YLim', [-3 3]);
newmap = jet(64);
ncol = size(newmap,1);          
zpos = 1 + floor(ncol);    
newmap(zpos,:) = [1 1 1];       
colormap(newmap);
xlabel('X (km)','fontsize',12);
ylabel('Time (days)','fontsize',12);
title(['Tidal Flat / Ubar (m/s) hovmoller'],'fontsize',14)
set(gca,'fontsize',15);

subplot(2,1,2)
pcolor(X,T,sand), shading flat
caxis([0 1]);
%cbr=colorbar('fontsize',13);
c=colorbar();
set(c, 'YLim', [0 1]);
newmap = jet(64);
ncol = size(newmap,1);          
zpos = 1 + floor(ncol);    
newmap(zpos,:) = [1 1 1];       
colormap(newmap);
xlabel('X (km)','fontsize',12);
ylabel('Time (days)','fontsize',12);
title(['Tidal Flat / Sand concentration (Kg/m3) hovmoller'],'fontsize',14)
set(gca,'fontsize',15);

set(gcf,'PaperPositionMode','auto');
export_fig -transparent tidal_flat.pdf

return

