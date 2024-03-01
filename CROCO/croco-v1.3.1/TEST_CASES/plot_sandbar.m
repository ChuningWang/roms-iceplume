%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make plot from the results of the SANDBAR test case
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
%  Patrick Marchesiello, IRD 2017,2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- LIP experiment --- 
%
%prompt = 'Erosion (1B) or Accretion (1C) LIP experiment? [1B]: ';
%mycase = input(prompt,'s');
%if isempty(mycase)
    mycase = '1B';
%end
%
% --- model params ---
%
fname  = 'sandbar_his.nc'; % croco file name

if mycase == '1B',
 morph_fac = 18;      % morphological factor (from sediment.in)
else
 morph_fac = 13; 
end

morph_cpl  = 1;       % feedback to currents
makepdf    = 0;       % make pdf file
%
%======================================================================
%
% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------
%
yindex = 2; % Mm=1 with NS no-slip conditions

nc=netcdf(fname,'r');
tindex  =length(nc{'scrum_time'}(:)); % reads last record

time=morph_fac*nc{'scrum_time'}(:)/3600; % time in hours
if mycase == '1B',
 [d,tindex0]=min(abs(time-4));     %  0-8h (1B)
else
 [d,tindex0]=min(abs(time-3));     %  0-7h (1C)
end

% horizontal grid
hr=squeeze(nc{'h'}(yindex,:));
xr=squeeze(nc{'x_rho'}(yindex,:));
hu=0.5*(hr(1:end-1)+hr(2:end));
xu=0.5*(xr(1:end-1)+xr(2:end));
L=length(hr);

% new bathy from last record
if morph_cpl,
 hnew=squeeze(nc{'hmorph'}(tindex,yindex,:));
 h=hnew;
 h0=squeeze(nc{'hmorph'}(tindex0,yindex,:));
 if isempty(hnew),
  h=hr;
  hnew=hr;
  h0=hr;
 end
else
 h=hr;
 hnew=hr;
 h0=hr;
end

% vertical grid
N=length(nc('s_rho'));
theta_s=nc.theta_s(:); 
theta_b=nc.theta_b(:); 
hc=nc.hc(:); 
zeta=squeeze(nc{'zeta'}(tindex,yindex,:));
Dcrit=1.1*nc{'Dcrit'}(:);
zeta(h<Dcrit)=zeta(h<Dcrit)-h(h<Dcrit); % add land topo
zr=squeeze(zlevs(h,zeta,theta_s,theta_b,hc,N,'r',2));
zw=squeeze(zlevs(h,zeta,theta_s,theta_b,hc,N,'w',2));
zru=0.5*(zr(:,1:end-1)+zr(:,2:end));
zwu=0.5*(zw(:,1:end-1)+zw(:,2:end));
dz1=zr(1,:)-zw(1,:);
dzu1=zru(1,:)-zwu(1,:);
dzu3=zru(3,:)-zwu(1,:);
%
xr2d=repmat(xr,[N 1]);
xu2d=repmat(xu,[N 1]);
xw2d=repmat(xr,[N+1 1]);
D   =zw(N+1,:)-zw(1,:);
D2d =repmat(D,[N 1]);
Du  =zwu(N+1,:)-zwu(1,:);
Du2d=repmat(Du,[N 1]);

% ---------------------------------------------------------------------
% --- read/compute model fields (tindex) ---
% --------------------------------------------------------------------
time=morph_fac/86400*(nc{'scrum_time'}(tindex)- ...
                      nc{'scrum_time'}(1));

% ... zonal velocity ...                         ---> xu,zu
u=squeeze(nc{'u'}(tindex,:,yindex,:));

% ... vertical velocity ...                      ---> xr,zw
w=squeeze(nc{'w'}(tindex,:,yindex,:));

% ... sediment concentration ...                 ---> xr,zr
C=2*squeeze(nc{'sand_01'}(tindex,:,yindex,:));

% ... total viscosity
Akv=squeeze(nc{'AKv'}(tindex,:,yindex,:));

% ... viscosity due to wave breaking ...
Akb=squeeze(nc{'Akb'}(tindex,:,yindex,:));

% ... wave setup ...  
sup=squeeze(nc{'zeta'}(tindex0,yindex,:)); % mid exp time
sup(hr<0)=sup(hr<0)+hr(hr<0)-Dcrit;

% ---------------------------------------------------------------------
% --- read/compute bottom model fields (tindx0) ---
% --------------------------------------------------------------------

% ... u undertow ...
zref=0.1; 
u0=squeeze(nc{'u'}(tindex0,:,yindex,:));
hu=0.5*(h0(1:end-1)+h0(2:end));
L=length(hu);
hu2d=repmat(hu,[N 1]);
z0=squeeze(nc{'zeta'}(tindex0,yindex,:));
z0(h0<Dcrit)=z0(h0<Dcrit)-h0(h0<Dcrit); % add land topo
zr0=squeeze(zlevs(h0,z0,theta_s,theta_b,hc,N,'r',2));
zru0=0.5*(zr0(:,1:end-1)+zr0(:,2:end));
zzu=hu2d+zru0;
for ix=1:L
 nn=min(find(zzu(:,ix)>zref));
 if isempty(nn), nn=N; end
 ubot(ix)=squeeze(u0(nn,ix));
end

% ... sediment concentration ...                 ---> xr,zr
C0=2*squeeze(nc{'sand_01'}(tindex0,:,yindex,:));
zref=0.05;
L=length(h0);
h2d=repmat(h0,[N 1]);
zz=h2d+zr0;
for ix=1:L
 nn=min(find(zz(:,ix)>zref));
 if isempty(nn), nn=N; end
 Cbot(ix)=squeeze(C0(nn,ix));
end

% ... hrms ...  
hrms =squeeze(nc{'hrm'}(tindex0,yindex,:)); % init time
close(nc)

zeta(D<Dcrit)=NaN;
u(Du2d<Dcrit)=NaN;
sup(D<=max(0.1,Dcrit))=NaN;
sup(hr<0)=NaN;
%
%======================================================================
%     FLUME DATA (LIP: ROELVINK & RENIER 1995)
%======================================================================
%
xd=linspace(0,200,64);

% --- LIP-1B bathy at 18h on xd grid ---
hd1B = -[ ...
-4.0935  -4.0931  -4.0926  -4.0921  -4.0917  -4.0912  -4.0217  -3.9788 ...
-3.8282  -3.6638  -3.5145  -3.3757  -3.2330  -3.0776  -2.9134  -2.7556 ...
-2.6207  -2.4878  -2.3855  -2.3182  -2.2756  -2.2480  -2.2272  -2.2064 ...
-2.1806  -2.1467  -2.1033  -2.0508  -1.9910  -1.9269  -1.8621  -1.8006 ...
-1.7456  -1.6995  -1.6625  -1.6319  -1.6011  -1.5583  -1.4985  -1.4276 ...
-1.3233  -1.0958  -0.8628  -0.8772  -0.9906  -1.0213  -1.0502  -1.0420 ...
-0.9964  -0.8379  -0.5619  -0.4863  -0.4529  -0.4259  -0.3967  -0.3527 ...
-0.2812  -0.1732  -0.0251   0.1585   0.3641   0.5696   0.7438   0.8490];

x_hrms_d1B=[20  65  100 115 130 138 145 152 160 170];      % --- Hs 8h
hrms_d1B  =[.85 .80 .72 .65 .58 .52 .39 .36 .36 .25];

x_ubot_d1B= [65   102  115  130  138  145  152  160  170]; % --- Ub 8h
  ubot_d1B=-[12.3 13.0 12.9 17.2 30.3 30.4 17.4 15.2 13.6];

x_Cbot_d1B= [65  102 115 130 138  145  152  170];          % --- Cb 8h
  Cbot_d1B= [.40 .23 .32 .77 3.08 1.73 .64  0.87];

% --- LIP-1C bathy at 13h on xd grid ---
hd1C = -[ ...
-4.0935  -4.0931  -4.0926  -4.0921  -4.0917  -4.0912  -4.0733  -3.9439 ...
-3.8105  -3.6717  -3.5289  -3.3835  -3.2361  -3.0864  -2.9344  -2.7828 ...
-2.6401  -2.4939  -2.3813  -2.3146  -2.2765  -2.2534  -2.2351  -2.2143 ...
-2.1864  -2.1494  -2.1028  -2.0481  -1.9878  -1.9250  -1.8632  -1.8055 ...
-1.7545  -1.7112  -1.6749  -1.6427  -1.6085  -1.5627  -1.4999  -1.4240 ...
-1.3120  -1.1622  -0.9755  -0.7739  -0.7177  -1.0002  -1.0464  -1.0356 ...
-0.9975  -0.8333  -0.5603  -0.4902  -0.4558  -0.4270  -0.4016  -0.3659 ...
-0.3035  -0.2013  -0.0536   0.1341   0.3454   0.5532   0.7225   0.8159];

x_hrms_d1C=[20  40  65 100 115 130 132 138 145 152 160 170]; % --- Hs 7h
hrms_d1C  =[.4 .41 .43 .44 .43 .43 .46 .43 .35 .33 .32 .22];

x_ubot_d1C= [65 102 115 125 130 134 152 160];                % --- Ub 7h
  ubot_d1C=-[ 1  1   1   2   2   3  13  11];

x_Cbot_d1C= [65  102 115 125 130  134 152 160];              % --- Cb 7h
  Cbot_d1C= [0.1 0.1 0.3 0.2 0.35 0.5 0.3 0.8];

if mycase=='1B',
 hd=hd1B;
 x_hrms_d=x_hrms_d1B;
 hrms_d=hrms_d1B;
 x_ubot_d=x_ubot_d1B;
 ubot_d=ubot_d1B;
 x_Cbot_d=x_Cbot_d1B;
 Cbot_d=Cbot_d1B;
else
 hm=interp1(xr,hr,xd);
 hd=hd1C;
 x_hrms_d=x_hrms_d1C;
 hrms_d=hrms_d1C;
 x_ubot_d=x_ubot_d1C;
 ubot_d=ubot_d1C;
 x_Cbot_d=x_Cbot_d1C;
 Cbot_d=Cbot_d1C;
end
%============================================================
% --- plot ---
%=============================================================
%
figure('Position',[50 50 600 800])
xmin=60; xmax=190;
%
% section u ...
%
h1=subplot(4,1,1);
h1.Position = h1.Position + [0 -0.12 0 0.15];
zmin=-2.5; zmax=0.5;
cmin=-0.5; cmax=0.5; nbcol=20;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
colormap(map);
contourf(xu2d,zru,u,[cmin:cint:cmax],'linestyle','-'); 
hold on
colorbar('h');
leg(1)=plot(xr,-hr,  'k:','LineWidth',3);
leg(2)=plot(xr,-hnew,'k', 'LineWidth',3);
leg(3)=plot(xd,-hd  ,'r', 'LineWidth',3);
plot(xr,zeta,'color','g', 'LineWidth',3);
legend(leg(1:3),'Initial','Final Model','Final data', ...
                'location','southeast');
ylabel('Depth [m]','Fontsize',15)
grid on
axis([xmin xmax zmin zmax])
caxis([cmin cmax])
thour=floor(time*24);
if mycase=='1B',
 title(['SANDBAR EROSION   LIP-1B - U at Time ',num2str(thour),' hour'])
else
 title(['SANDBAR ACCRETION LIP-1C - U at Time ',num2str(thour),' hour'])
end
set(gca,'Fontsize',15)
hold off
%
% hrms ...
%
h2=subplot(4,1,2);
h2.Position = h2.Position + [0 -0.05 0 -0.05];
zmin=0; zmax=1;
plot(xr,hrms,'b',x_hrms_d,hrms_d,'b*','LineWidth',2); hold on
legend('Model','Flume','Location','SouthWest')
grid on
axis([xmin xmax zmin zmax])
thour=floor(time*24);
ylabel('Hrms [m]','Fontsize',15)
set(gca,'Fontsize',15)
hold off
%
% Undertow ...
%
h3=subplot(4,1,3);
h3.Position = h3.Position + [0 -0.025 0 -0.05];
zmin=-40; zmax=0;
plot(xu,100*ubot,'b',x_ubot_d,ubot_d,'b*','LineWidth',2); hold on;
legend('Model','Flume','Location','SouthWest')
grid on
axis([xmin xmax zmin zmax])
thour=floor(time*24);
ylabel('U [cm/s]','Fontsize',15)
set(gca,'Fontsize',15)
hold off
%
% Sand Concentration ...
%
h4=subplot(4,1,4);
h4.Position = h4.Position + [0 0.0 0 -0.05];
zmin=0; zmax=4;
plot(xr,Cbot,'b',x_Cbot_d,Cbot_d,'b*','LineWidth',2); hold on;
legend('Model','Flume','Location','SouthWest')
grid on
axis([xmin xmax zmin zmax])
thour=floor(time*24);
xlabel('X [m]','Fontsize',15)
ylabel('C [g/l]','Fontsize',15)
set(gca,'Fontsize',15)
hold off

if makepdf
 if mycase=='1B',
  export_fig -transparent sandbar_LIP_1B.pdf
 else 
  export_fig -transparent sandbar_LIP_1C.pdf
 end
end


