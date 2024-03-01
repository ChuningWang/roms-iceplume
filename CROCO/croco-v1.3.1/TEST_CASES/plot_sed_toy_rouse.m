%======================================================================
%
%     ---               Sed_toy Rouse Test Case                ---     
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

fname='sed_toy_rouse_his.nc';
% Choose usgs (1) or mustang (0) sediment model output file
usgs=1


%======================================================================
%
% Read data
%

% model parameters
dt=1  ; nz=100 ; it=1;

% filter
nt1=12 ; idy=4 ; idx=4;

nc=netcdf(fname);

h=squeeze(nc{'h'}(:,:));
x=squeeze(nc{'x_rho'}(:,:));
y=squeeze(nc{'y_rho'}(:,:));
zeta=squeeze(nc{'zeta'}(nt1,idy,idx));
L=length(zeta);
if usgs == 1
  model='Usgs';
  tauskin=squeeze(nc{'bostr'}(nt1,idy,idx));
else
  model='Mustang';
  tauskin=squeeze(nc{'TAUSKIN'}(nt1,idy,idx));
end
u=squeeze(nc{'u'}(nt1,:,idy,idx));
N=size(u,1);
%
theta_s=nc.theta_s(:);
theta_b=nc.theta_b(:);
hc=nc.hc(:);
Vtrans=nc{'Vtransform'}(:);
%
depth=h(idy,idx);
zr=zlevs(depth,zeta,theta_s,theta_b,hc,N,'r',Vtrans);
zr=squeeze(zr);
zu=0.5*(zr(1:end-1)+zr(2:end));
%
h=depth+zeta;
z=-zr;

% Constant values 
ws=[0.001, 0.01, 0.02, 0.04, 0.08, 0.1]; % setling velocity (m/s)
if usgs == 1
  sandstr=["mud_01", "mud_02", "mud_03", "mud_04", "mud_05", "mud_06"];
else
  sandstr=["SED1", "SED2", "SED3", "SED4", "SED5", "SED6"];
end
rho=1030.; % 1025 ds croco.in
vk=0.41;
vk_hydro=0.41;
z0_val=0.0001;
nsand=6;
k_aref=1;
a=z(1);

%----------------------------------------------------------
% Compute rouse value
%----------------------------------------------------------
uet=sqrt(tauskin/rho);
rouse=ws/(vk*uet);

%----------------------------------------------------------
% Concentration tracer (Model/Rouse)
%----------------------------------------------------------
c_ana=zeros(nsand,nz);
c_sand=zeros(nsand,nz);

for isand=1:nsand
 sd=sandstr(isand);
 sandn=squeeze(nc{sd{1}}(nt1,:,idy,idx));
 c_sand(isand,:)=sandn(:,it);
 for k=1:nz
  c_ana(isand,k)=c_sand(isand,k_aref)*( ((h-z(k))/z(k))*(a./(h-a)) ).^rouse(isand);
 end
end

close(nc);

%----------------------------------------------------------
% Normalize Concentration (Model/Rouse)
%----------------------------------------------------------
rapC_ana=zeros(nsand,nz);
rapC_sand=zeros(nsand,nz);

for isand=1:nsand
 rapC_sand(isand,:)=c_sand(isand,:)/c_sand(isand,1);
 rapC_ana(isand,:)=fliplr(c_ana(isand,:)/c_ana(isand,100));
end

%----------------------------------------------------------
%  Plot Vertical profile (model:dt=1s / nz=100)
%----------------------------------------------------------

figure('Position',[1 1 1500 500])

title0=['Vertical profile of sand concentrations at equilibrium / Rouse profile (model:dt=',num2str(dt),'s,nz=',num2str(nz),')'];

for isand=1:nsand
	ax=subplot(1,nsand,isand);
    plot(fliplr(rapC_sand(isand,:)),z, 'g')
	hold on
    plot(fliplr(rapC_ana(isand,:)),z, '--')
	grid on
    xlabel('C/C[k=1]','fontsize',12)
	if isand == 1
	   ylabel('Water level [m]','fontsize',12)
	end
    title(['WS=',num2str(ws(isand)),' m/s'])
    xlim([0 1])
	hold off
end

legend(model,'Rouse')
set(gcf,'PaperPositionMode','auto');
set(gcf,'Name',title0)

export_fig -transparent sed_toy.pdf


return

