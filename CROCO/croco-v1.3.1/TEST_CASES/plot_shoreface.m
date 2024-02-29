%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make plot from the results of the SHOREFACE test case
% 
%  Further Information:  
%  http://www.crocoagrif.org/
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
%  Ref: Penven, P., L. Debreu, P. Marchesiello and J.C. McWilliams,
%       Application of the ROMS embedding procedure for the Central 
%      California Upwelling System,  Ocean Modelling, 2006.
%
%  Patrick Marchesiello, IRD 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
fname     = 'shoreface_his.nc';    % croco file name
yindex    = 1;                     % y index
makepdf   = 0;                     % make pdf file 
%
%======================================================================

if exist('subplot_tight')==2,
 mysubplot = str2func('subplot_tight');
else
 mysubplot = str2func('subplot');
end

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

nc=netcdf(fname,'r');
tindex=length(nc{'scrum_time'}(:)); % reads last record

%
% horizontal grid
hr=squeeze(nc{'h'}(yindex,:));
xindex=1;
hr=hr(xindex:end);
L=length(hr);
xr=squeeze(nc{'x_rho'}(yindex,xindex:end));
yr=squeeze(nc{'y_rho'}(yindex,xindex:end));
dx=xr(2)-xr(1);
%
% vertical grid
N=length(nc('s_rho'));
theta_s=nc.theta_s(:); 
theta_b=nc.theta_b(:); 
hc=nc.hc(:); 
zeta=squeeze(nc{'zeta'}(tindex,yindex,xindex:end));
zr=squeeze(zlevs(hr,zeta,theta_s,theta_b,hc,N,'r',2));
dzr=zr(2:end,:)-zr(1:end-1,:);               % ---> zw(2:N,:)
zru=0.5*(zr(:,1:end-1)+zr(:,2:end));
dzru=zru(2:end,:)-zru(1:end-1,:);            % ---> zwu(2:N,:)
zw=squeeze(zlevs(hr,zeta,theta_s,theta_b,hc,N,'w',2));
dzw=zw(2:end,:)-zw(1:end-1,:);               % ---> zr
zwu=0.5*(zw(:,1:end-1)+zw(:,2:end));
dzwu=zwu(2:end,:)-zwu(1:end-1,:);            % ---> zru
%
xr2d=repmat(xr,[N 1]);
xw2d=repmat(xr,[N+1 1]);
D=hr+zeta;
D2d=repmat(D,[N 1]);

% ---------------------------------------------------------------------
% --- read/compute numerical model fields (index 1) ---
% --------------------------------------------------------------------
time=nc{'scrum_time'}(tindex)/86400;

zeta1=zeta;

% ... zonal velocity ...                         ---> xu,zru
u1=squeeze(nc{'u'}(tindex,:,yindex,xindex:end));

% ... meridional velocity ...                    ---> xr,zr
v1=squeeze(nc{'v'}(tindex,:,yindex,xindex:end));

% ... vertical velocity ...                      ---> xr,zw
w1=squeeze(nc{'w'}(tindex,:,yindex,xindex:end));

% ... temperature ...                            ---> xr,zr
t1=squeeze(nc{'temp'}(tindex,:,yindex,xindex:end));

% ... vertical viscosity/diffusivity
Akv=squeeze(nc{'AKv'}(tindex,:,yindex,xindex:end));
Akt=squeeze(nc{'AKt'}(tindex,:,yindex,xindex:end));


% ... wave setup ...  
sup=squeeze(nc{'sup'}(tindex,yindex,xindex:end));

% ... zonal Stokes dritf                         ---> xu,zru
ust=squeeze(nc{'ust'}(tindex,:,yindex,xindex:end));

% ... meridional Stokes drift ...                ---> xr,zr
vst=squeeze(nc{'vst'}(tindex,:,yindex,xindex:end));

% ... vertical Stokes drift ...                  ---> xr,zw
wst=squeeze(nc{'wst'}(tindex,:,yindex,xindex:end));

% eddy viscosity due to depth-induced wave breaking
Akb=squeeze(nc{'Akb'}(tindex,:,yindex,xindex:end));

% eddy diffusivity due to primary waves
Akw=squeeze(nc{'Akw'}(tindex,:,yindex,xindex:end));

close(nc)

%============================================================
% --- plot ---
%=============================================================
Dcrit=0.2;

xr=1.e-3*xr-1;
xr2d=1.e-3*xr2d-1;
xw2d=1.e-3*xw2d-1;
u1(:,L)=u1(:,L-1);
zeta1(D<Dcrit)=NaN;
u1(D2d<Dcrit)=NaN;
v1(D2d<Dcrit)=NaN;
ust(:,L)=ust(:,L-1);
Akv(:,1)=NaN;
Akt(:,1)=NaN;
Akb(:,1)=NaN;
Akw(:,1)=NaN;

%-----------------------------------
% Eulerian velocities u,v,w
%
figure('Position',[100 150 800 600])

mysubplot(3,3,1)
cmin=-0.5; cmax=-cmin; nbcol=20;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
colormap(map);
contourf(xr2d,zr,u1,[cmin:cint:cmax],'linestyle','none'); 
hold on
colorbar;
plot(xr,-hr,'color','k','LineWidth',3);
plot(xr,zeta1,'color','g','LineWidth',3);
axis([-1 0 -12 1])
caxis([cmin cmax])
ylabel('Z [m]')
thour=floor(time*24);
text(-0.05,-11,'U [m/s]','HorizontalAlignment','right')
hold off

mysubplot(3,3,4)
cmin=-1.2; cmax=-cmin; nbcol=20;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
colormap(map);
contourf(xr2d,zr,v1,[cmin:cint:cmax],'linestyle','none'); 
hold on
colorbar;
plot(xr,-hr,'color','k','LineWidth',3);
plot(xr,zeta1,'color','g','LineWidth',3);
axis([-1 0 -12 1])
caxis([cmin cmax])
ylabel('Z [m]')
text(-0.05,-11,'V [m/s]','HorizontalAlignment','right')
hold off

mysubplot(3,3,7)
cmin=-7.e-3; cmax=-cmin; nbcol=20;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
colormap(map);
contourf(xr2d,zr,w1,[cmin:cint:cmax],'linestyle','none'); 
hold on
colorbar;
plot(xr,-hr,'color','k','LineWidth',3);
plot(xr,zeta1,'color','g','LineWidth',3);
axis([-1 0 -12 1])
xlabel('Distance to shore [km]')
caxis([cmin cmax])
ylabel('Z [m]')
text(-0.05,-11,'W [m/s]','HorizontalAlignment','right')
hold off

%-----------------------------------
% Stokes drift ust,vst,wst
%
mysubplot(3,3,2)
cmin=-0.2; cmax=0.2; nbcol=20;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
colormap(map);
contourf(xr2d,zr,ust,[cmin:cint:cmax],'linestyle','none'); 
hold on
colorbar;
plot(xr,-hr,'color','k','LineWidth',3);
plot(xr,zeta1,'color','g','LineWidth',3);
axis([-1 0 -12 1])
caxis([cmin cmax])
text(-0.05,-11,'U Stokes [m/s]','HorizontalAlignment','right')
title(['SHOREFACE - Time = ',num2str(thour),' h'],'fontsize',15)
hold off

mysubplot(3,3,5)
cmin=-2.e-2; cmax=2.e-2; nbcol=20;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
colormap(map);
contourf(xr2d,zr,vst,[cmin:cint:cmax],'linestyle','none'); 
hold on
colorbar;
plot(xr,-hr,'color','k','LineWidth',3);
plot(xr,zeta1,'color','g','LineWidth',3);
axis([-1 0 -12 1])
caxis([cmin cmax])
text(-0.05,-11,'V Stokes [m/s]','HorizontalAlignment','right')
hold off

mysubplot(3,3,8)
cmin=-2.e-3; cmax=2.e-3; nbcol=20;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
colormap(map);
contourf(xr2d,zr,wst,[cmin:cint:cmax],'linestyle','none'); 
hold on
colorbar;
plot(xr,-hr,'color','k','LineWidth',3);
plot(xr,zeta1,'color','g','LineWidth',3);
axis([-1 0 -12 1])
xlabel('Distance to shore [km]')
caxis([cmin cmax])
text(-0.05,-11,'W Stokes [m/s]','HorizontalAlignment','right')
hold off

%-------------------------------------------------
% Turbulent and wave-induced Viscosity/diffusivity
%
mysubplot(3,3,3)
cmin=-0.05; cmax=-cmin; nbcol=20;
cint=(cmax-cmin)/nbcol;
contourf(xw2d,zw,Akv,[cmin:cint:cmax],'linestyle','none'); 
hold on
colorbar;
plot(xr,-hr,'color','k','LineWidth',3);
axis([-1 0 -12 1])
ylabel('Z [m]')
caxis([cmin cmax])
text(-0.05,-11,'Akv (total visc) [m^2/s]','HorizontalAlignment','right')
hold off

mysubplot(3,3,6)
cmin=-0.05; cmax=-cmin; nbcol=20;
cint=(cmax-cmin)/nbcol;
contourf(xw2d,zw,Akb,[cmin:cint:cmax],'linestyle','none'); 
hold on
colorbar;
plot(xr,-hr,'color','k','LineWidth',3);
axis([-1 0 -12 1])
caxis([cmin cmax])
ylabel('Z [m]')
text(-0.05,-11,'Akb (break. wave visc.)','HorizontalAlignment','right')
hold off

mysubplot(3,3,9)
cmin=-3.e-3; cmax=-cmin; nbcol=20;
cint=(cmax-cmin)/nbcol;
contourf(xw2d,zw,Akw,[cmin:cint:cmax],'linestyle','none'); 
hold on
colorbar;
plot(xr,-hr,'color','k','LineWidth',3);
axis([-1 0 -12 1])
xlabel('Distance to shore [km]')
caxis([cmin cmax])
ylabel('Z [m]')
text(-0.05,-11,'Akw (non-break. wave diff.)','HorizontalAlignment','right')
hold off

if makepdf
 export_fig -transparent shoreface.pdf
end

return

%--------------------------------------------------------------------
%
% Surface elevation + wave setup
% 
figure('Position',[100 150 500 750])
xmin=-1; xmax=-0; zmin=-0.1; zmax=0.25;
plot(xr,zeta1+sup,'color','b','LineWidth',3);
grid on
axis([xmin xmax zmin zmax])
thour=floor(time*24);
title(['SHOREFACE: Wave Setup at ',num2str(thour),' h'])
hold off
if makepdf
 export_fig -transparent shoreface_setup.pdf
end


