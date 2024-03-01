%======================================================================
%
%     ---               Internal Test Case                ---     
%     ---  Internal Gravity Wave solution over a ridge    ---
%
% Reference:
% ----------
% Di Lorenzo, E, W.R. Young and S.L. Smith, 2006: Numerical and 
% anlytical estimates of M2 tidal conversion at steep oceanic ridges, 
% J. Phys. Oceanogr., 36, 1072-1084.
%
%  Further Information:  
%  http://www.croco-ocean.org
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

fname='internal_his.nc';

makepdf=0;
%======================================================================
%
% Read data
%
j=3;
nc=netcdf(fname);
time=nc{'scrum_time'}(:)./86400;
tndx=length(time);
disp([ 'tndx = ',num2str(tndx), ...
       ' - Time = ',num2str(time(end)*24/12.4),' M2 cycles' ])
h=squeeze(nc{'h'}(j,:));
x=squeeze(nc{'x_rho'}(j,:));
zeta=squeeze(nc{'zeta'}(tndx,j,:));
L=length(zeta);
%
t=squeeze(nc{'rho'}(tndx,:,j,:));
t0=squeeze(nc{'rho'}(1,:,j,:));
[N,L]=size(t);
u=squeeze(nc{'u'}(tndx,:,j,:));
v=squeeze(nc{'v'}(tndx,:,j,:));
theta_s=nc.theta_s(:);
theta_b=nc.theta_b(:);
hc=nc.hc(:);
close(nc);
%
zr = zlevs(h,zeta,theta_s,theta_b,hc,N,'r',1);
zr=squeeze(zr);
xr=reshape(x,1,L);
xr=repmat(xr,[N 1])/1000;
%
%  Plot rho anomaly
%
figure('position',[500 500 800 400])
contourf(xr,zr,t-t0,40,'linestyle','none')
hold on
line(xr,-h,'color','k','Linewidth',2)
map=colormap(jet(40));
map(20:21,:)=[1 1 1; 1 1 1];
colormap(map)
caxis([-0.01 0.01])
colorbar
title(['Internal case - rho anomaly at 12 M2 cycles'])
set(gca,'fontsize',15);
set(gcf,'PaperPositionMode','auto');
hold off

if makepdf,
 export_fig -transparent internal_rho.pdf
end

return

