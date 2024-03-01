%======================================================================
%
%     ---               Dune3d Test Case                ---     
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

% Choose usgs (1) or mustang (0) model output file
usgs=1

fname='dune3d_his.nc';
tdays=2; 

%======================================================================
%
% Read data
%
tndx=tdays+1; 

nc=netcdf(fname);

h=squeeze(nc{'h'}(:,:));
X=squeeze(nc{'x_rho'}(:,:));
Y=squeeze(nc{'y_rho'}(:,:));
%

if usgs ~= 1
    model='MUSTANG';
    var_hmorph='Hm';
    hmorph0  =squeeze(nc{var_hmorph}(1,:,:));
    hmorph   =squeeze(nc{var_hmorph}(tndx,:,:));
else
    model='USGS';
    var_hmorph='hmorph';
    hmorph0  =squeeze(nc{var_hmorph}(1,:,:));
    hmorph  =squeeze(nc{var_hmorph}(tndx,:,:));
end

close(nc);

%----------------------------------------------------------
%  Plot 2D bed evolution 
%----------------------------------------------------------
figure('Position',[1 1 1000 500])


% plot at t=0  
ax1=subplot(1,2,1);

pcolor(X,Y,-squeeze(hmorph0(:,:))), shading flat
caxis([-4 -2]);
cbr=colorbar('fontsize',13);
%cbr.Label.String = 'Hm';
xlabel('X (m)','fontsize',12);
ylabel('Y (m)','fontsize',12);

hold on
[C,h]=contour(X,Y,-squeeze(hmorph0(:,:)),[-3.8:0.4:-2],'k');
v = [-3.8:0.4:-2];
clabel(C,h,v,'FontSize',7,'Color','red');

%grid on
title(['Hm (t0=0)'],'fontsize',14)
set(gca,'fontsize',15);

% plot at t
ax2=subplot(1,2,2);

pcolor(X,Y,-squeeze(hmorph(:,:))), shading flat
caxis([-4 -2]);
cbr=colorbar('fontsize',13);
cbr.Label.String = 'Hm';
xlabel('X (m)','fontsize',12);

hold on
[C,h]=contour(X,Y,-squeeze(hmorph(:,:)),[-3.8:0.4:-2],'k');
v = [-3.8:0.4:-2];
clabel(C,h,v,'FontSize',7,'Color','red');
%grid on

title(['Hm (t= t0 + ',num2str(tdays),' days) / ', model],'fontsize',14)

set(gca,'fontsize',15);
set(gcf,'PaperPositionMode','auto');
hold off

if usgs ~= 1
 export_fig -transparent dune3d_mustang.pdf
else
 export_fig -transparent dune3d_usgs.pdf
end


return

