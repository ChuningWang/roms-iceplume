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
fname     = 'thacker_his.nc';      % croco file name
x0        = 101;                   % x and y origins
y0        = 2;                     %
makemovie = 0;                     % make movie using QTWriter
makepdf   = 0;                     % make pdf file
%
%======================================================================

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

nc=netcdf(fname);
tindex=length(nc{'scrum_time'}(:)); % reads last record

if makemovie,
 movObj = QTWriter('thacker.mov');
 tstr=1;
 tend=tindex;
else,
 tstr=tindex;
 tend=tstr;
end

hf = figure
axis tight; set(hf,'DoubleBuffer','on');
set(gca,'nextplot','replacechildren');

for tindex=tstr:tend % ---------------------------------------------

%
% horizontal grid
 hr=squeeze(nc{'h'}(y0,:));
 xindex=1;
 hr=hr(xindex:end);
 L=length(hr);
 xr=squeeze(nc{'x_rho'}(y0,xindex:end));
 yr=squeeze(nc{'y_rho'}(y0,xindex:end));
 dx=xr(2)-xr(1);
 Dcrit=nc{'Dcrit'}(:);
%
% vertical grid
 N=length(nc('s_rho'));
 theta_s=nc.theta_s(:); 
 theta_b=nc.theta_b(:); 
 hc=nc.hc(:); 
 zeta=squeeze(nc{'zeta'}(tindex,y0,xindex:end));
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
%
 if zeta(1)<Dcrit+0.1, % check if zeta needs redef on dry land
  zeta(hr<Dcrit)=zeta(hr<Dcrit)-hr(hr<Dcrit);
 end
 D=hr+zeta;
 D2d=repmat(D,[N 1]);

 % ---------------------------------------------------------------------
 % --- read/compute numerical model fields (index 1) ---
 % ---------------------------------------------------------------------
 time=nc{'scrum_time'}(tindex);

 zeta1=zeta;

 % ... num zonal velocity ...                         ---> xu,zru
 u1=squeeze(nc{'u'}(tindex,:,y0,xindex:end));

 % ... num meridional velocity ...                    ---> xr,zr
 v1=squeeze(nc{'v'}(tindex,:,y0,xindex:end));

 % ... num vertical velocity ...                      ---> xr,zw
 w1=squeeze(nc{'w'}(tindex,:,y0,xindex:end));

 % ... num temperature ...                            ---> xr,zr
 t1=squeeze(nc{'temp'}(tindex,:,y0,xindex:end));

 % ---------------------------------------------------------------------
 % --- compute analytical solutions (index 2) ---
 % ---------------------------------------------------------------------

 eta = 0.1;                % --> nondimensional periodic amplitude
 D0  = 10;                 % --> max depth at rest
 Lt  = 80.e3;              % --> distance at psi points
 f   = nc{'f'}(y0,1);      % --> 0 (2D case)
 g   = 9.81;               % gravity acceleration (m^2/s)

 omega=sqrt(f^2 +2*g*D0/Lt^2);

 u2    = -eta*omega*Lt*sin(omega*time)*ones(size(u1));
 v2    = -eta*omega*Lt*cos(omega*time)*ones(size(v1));
 zeta2 =  2*eta*D0/Lt*(xr.*cos(omega*time)-0.5*eta) ...
                                    .*ones(size(zeta1));

 %============================================================
 % --- plot ---
 %=============================================================

 xr=1.e-3*xr;
 xr2d=1.e-3*xr2d;
 u1(:,L)=u1(:,L-1);
 u2(:,L)=u2(:,L-1);
 zeta1(D<Dcrit+0.01)=NaN;
 u1(D2d<Dcrit)=NaN;
 zeta2(zeta2<-hr)=NaN;

 cmin=-100; cmax=-cmin; nbcol=20;
 cint=(cmax-cmin)/nbcol;
 map=colormap(jet(nbcol));
 map(nbcol/2  ,:)=[1 1 1];
 map(nbcol/2+1,:)=[1 1 1];
 colormap(map);

 uerr=100*(u1-u2)./u2;
 contourf(xr2d,zr,uerr,[cmin:cint:cmax],'linestyle','none');  
 hold on
 colorbar;
 plot(xr,-hr,'color','k','LineWidth',3);
 h=plot(xr,zeta2,'g',xr,zeta1,'r','LineWidth',3);
 legend(h,'Analytical','Numerical')
 hold off
 axis([-100 100 -10 5])
 xlabel('X [km]')
 ylabel('Z [m]')
 caxis([cmin cmax])
 grid on
 thour=floor(time/3600);
 set(gca,'fontsize',15)
 title(['THACKER: \eta [m] and U Error [%] at ',num2str(thour),' hour'])

 if makemovie,  
  % Write each frame to the file
  movObj.FrameRate =  5;
  writeMovie(movObj,getframe(hf));
  clf('reset')
 end
%----------------------------------

end % time loop

%if makepdf
% export_fig -transparent thacker_72h.pdf
%end

if makemovie,  
 movObj.Loop = 'loop'; % Set looping flag
 close(movObj);        % Finish writing movie and close file
end

%============================================================
% --- plot sea level extremes at 6, 9, 12 h ---
%============================================================

if tindex<64*2; 
 close(nc);
 return
end

xr=xr*1.e3;
tstr=1;

tindex=tstr+60*2; 
time=nc{'scrum_time'}(tindex);
zm=squeeze(nc{'zeta'}(tindex,y0,xindex:end));
za=2*eta*D0/Lt*(xr.*cos(omega*time)-0.5*eta) ...
                                 .*ones(size(zeta1));
if zm(1)<Dcrit+0.1, % check if zeta needs redef on dry land
 zm(hr<Dcrit)=zm(hr<Dcrit)-hr(hr<Dcrit);
end
D=hr+zm;
zm(D<Dcrit+0.01)=NaN;
za(za<-hr)=NaN;
zm_60h=zm; za_60h=za;

tindex=tstr+62*2;
time=nc{'scrum_time'}(tindex);
zm=squeeze(nc{'zeta'}(tindex,y0,xindex:end));
za=2*eta*D0/Lt*(xr.*cos(omega*time)-0.5*eta) ...
                                 .*ones(size(zeta1));
if zm(1)<Dcrit+0.1, % check if zeta needs redef on dry land
 zm(hr<Dcrit)=zm(hr<Dcrit)-hr(hr<Dcrit);
end
D=hr+zm;
zm(D<Dcrit+0.01)=NaN;
za(za<-hr)=NaN;
zm_62h=zm; za_62h=za;

tindex=tstr+64*2;
time=nc{'scrum_time'}(tindex);
zm=squeeze(nc{'zeta'}(tindex,y0,xindex:end));
za=2*eta*D0/Lt*(xr.*cos(omega*time)-0.5*eta) ...
                                 .*ones(size(zeta1));
if zm(1)<Dcrit+0.1, % check if zeta needs redef on dry land
 zm(hr<Dcrit)=zm(hr<Dcrit)-hr(hr<Dcrit);
end
D=hr+zm;
zm(D<Dcrit+0.01)=NaN;
za(za<-hr)=NaN;
zm_64h=zm; za_64h=za;

close(nc);

xr=xr*1.e-3;
figure
plot(xr,-hr,'color','k','LineWidth',4); hold on;
h1=plot(xr,za_60h,'g-',xr,zm_60h,'r-','LineWidth',3);
text(60,zm_60h(160)+0.3,'60 h','fontsize',15);
h2=plot(xr,za_62h,'g-',xr,zm_62h,'r-','LineWidth',3);
text(60,zm_62h(160)+0.2,'62 h','fontsize',15);
h3=plot(xr,za_64h,'g-',xr,zm_64h,'r-','LineWidth',3);
text(60,zm_64h(160)+0.1,'64 h','fontsize',15);
legend(h3,'Analytical','Numerical','location','Northwest')
axis([20 100 -3 3])
grid on
hold off
xlabel('X [km]')
ylabel('Z [m]')
set(gca,'fontsize',15)
title(['THACKER: sea level at various times'])

if makepdf
 export_fig -transparent thacker_zcomp.pdf
end

return

%============================================================
% --- plot u time series at center point ---
%=============================================================
nc=netcdf(fname);
t0=nc{'scrum_time'}(1:tindex);
u10=squeeze(nc{'u'}(1:tindex,2,y0,x0));
u20=-eta*omega*Lt*sin(omega*t0);

figure
t0=t0/86400;
h=plot(t0,u20,'g',t0,u10,'r');
set(h,'linewidth',2)
legend('Analytical','Numerical')
title('U timeseries at center point')
grid on
if makepdf
 export_fig -transparent thacker_Useries.pdf
end
close(nc);




