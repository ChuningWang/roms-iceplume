%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make plot from the results of the SHOREFACE test case
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
%  Ref: Penven, P., L. Debreu, P. Marchesiello and J.C. McWilliams,
%       Application of the CROCO embedding procedure for the Central 
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
fname     = 'swash_his.nc';   % croco history file name
varname   = 'u';              % var name [ 'u' 'w' ]

makemovie = 0;                % make movie using QTWriter
makepdf   = 0;                % make pdf file
%
%======================================================================

yindex = 2;
g = 9.81;

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

nc=netcdf(fname);
tindex=length(nc{'scrum_time'}(:)); % reads last record

if makemovie,
 movObj = QTWriter('swash.mov');
 tstr=1;
 tend=tindex;
else,
 tstr=tindex;
 tend=tstr;
end

hf = figure('position',[1000 500 800 400]);
axis tight; set(hf,'DoubleBuffer','on');
set(gca,'nextplot','replacechildren');

for tindex=tstr:tend % ---------------------------------------------
%
% horizontal grid
 hr=squeeze(nc{'h'}(yindex,:));
 L=length(hr);
 xr=squeeze(nc{'x_rho'}(yindex,:));
 yr=squeeze(nc{'y_rho'}(yindex,:));
 dx=xr(2)-xr(1);
 xmin=min(xr);
 xmax=max(xr);
 Dcrit=nc{'Dcrit'}(:);
%
% vertical grid
 N=length(nc('s_rho'));
 theta_s=nc.theta_s(:); 
 theta_b=nc.theta_b(:); 
 hc=nc.hc(:); 
 zeta=squeeze(nc{'zeta'}(tindex,yindex,:));
 zr=squeeze(zlevs(hr,zeta,theta_s,theta_b,hc,N,'r',2));
 zw=squeeze(zlevs(hr,zeta,theta_s,theta_b,hc,N,'w',2));
 dzr=zr(2:end,:)-zr(1:end-1,:);         % ---> zw(2:N,:)
 zru=0.5*(zr(:,1:end-1)+zr(:,2:end));
 dzru=zru(2:end,:)-zru(1:end-1,:);      % ---> zwu(2:N,:)
 dzw=zw(2:end,:)-zw(1:end-1,:);         % ---> zr
 zwu=0.5*(zw(:,1:end-1)+zw(:,2:end));
 dzwu=zwu(2:end,:)-zwu(1:end-1,:);      % ---> zru
%
 xr2d=repmat(xr,[N 1]);
 D=hr+zeta;
 D2d=repmat(D,[N 1]);
%
 xmin= min(xr); xmax=90;
 zmin=-max(hr); zmax=0.4*max(hr);

 % ---------------------------------------------------------------------
 % --- read/compute numerical model fields (index 1) ---
 % ---------------------------------------------------------------------
 time=nc{'scrum_time'}(tindex);

 % ... zonal velocity ...                         ---> xu,zru
 u=squeeze(nc{'u'}(tindex,:,yindex,:));
 u(:,2:L-1)=0.5*(u(:,1:L-2)+u(:,2:L-1));
 u(:,1)=u(:,2);u(:,L)=u(:,L-1);

 % ... vertical velocity ...                      ---> xr,zw
 w=squeeze(nc{'w'}(tindex,2:end,yindex,:));

 %=============================================================
 % --- plot ---
 %=============================================================
 
 if varname(:,1)=='u', var=u;   end
 if varname(:,1)=='w', var=w;   end

 if varname(:,1)=='u', cmin=-1.50; cmax=-cmin; nbcol=20; end % u
 if varname(:,1)=='w', cmin=-0.70; cmax=-cmin; nbcol=20; end % w
 cint=(cmax-cmin)/nbcol;

 Dcrit=Dcrit+1.e-6;
 var(D2d<Dcrit)=NaN;

 ztop=squeeze(zw(end,:,:));
 ztop(D<Dcrit+0.01)=NaN;

 map=colormap(jet(nbcol));
 map(nbcol/2  ,:)=[1 1 1];
 map(nbcol/2+1,:)=[1 1 1];
 colormap(map);
 
 contourf(xr2d,zr,var,[cmin:cint:cmax],'LineStyle','none'); hold on
 colorbar;
 plot(xr,-hr,'color','k','LineWidth',3);
 hn=plot(xr,ztop,'color','r','LineWidth',2);
 grid on
 axis([xmin xmax zmin zmax])
 caxis([cmin cmax])
 tmin=floor(time/60);
 title(['SWASH: ',varname,' at ',num2str(time),' sec'])
 hold off
 set(gca,'fontsize',15);
 set(gcf,'PaperPositionMode','auto');

 if makemovie,  
  % Write each frame to the file
  movObj.FrameRate =  5;
  writeMovie(movObj,getframe(hf));
  clf('reset')
 end
%----------------------------------

end % time loop

if makemovie,  
    movObj.Loop = 'loop'; % Set looping flag
    close(movObj);        % Finish writing movie and close file
end

if makepdf
 export_fig -transparent swash.pdf
end

return




