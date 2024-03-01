%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make plot from the results of the RIP test case
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
%  Patrick Marchesiello, IRD 2013,2017
%   uses QTwriter.m for creating .mov video files if desired
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
fname     = 'rip_his.nc';  % croco file name
makemovie = 0;             % make movie using QTWriter (for tidal case)
makepdf   = 0;             % make pdf file
pltvort   = 0;             % plot vort (else plot u)
%
%======================================================================

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

nc=netcdf(fname);
tindex=length(nc{'scrum_time'}(:)); % read last record
tstr=2;
tend=tindex;

if makemovie,
 movObj = QTWriter('rip.mov');
else
 tstr=tindex;
end

h=nc{'h'}(:);
x=nc{'x_rho'}(:);
y=nc{'y_rho'}(:);
pm=nc{'pm'}(:);
pn=nc{'pn'}(:);
N=length(nc('s_rho'));

% min depth in wetting-drying scheme
Dcrit=nc{'Dcrit'}(:);

% initialize plot
hf = figure;
axis tight; set(hf,'DoubleBuffer','on');
set(gca,'nextplot','replacechildren');

% ---------------------------------------------------------------------
% --- get variables and plot in time loop ---
% ---------------------------------------------------------------------

for tindex=tstr:tend  % ---

 disp(['Tindex = ',num2str(tindex)])

 time=nc{'scrum_time'}(tindex)/86400;
 zeta=nc{'zeta'}(tindex,:,:);
 u=nc{'u'}(tindex,N,:,:);
 v=nc{'v'}(tindex,N,:,:);
 ub=nc{'ubar'}(tindex,:,:)+nc{'ust2d'}(tindex,:,:);
 vb=nc{'vbar'}(tindex,:,:)+nc{'vst2d'}(tindex,:,:);

 mask=ones(size(zeta));
 mask((h+zeta)<=Dcrit)=NaN;
 ur=u2rho_2d(ub);
 vr=v2rho_2d(vb);
 if pltvort,
  vort=mask.*psi2rho(vorticity(u,v,pm,pn));
 else
  speed=mask.*ur;  
 end

 %======================================
 % --- plot ---
 %======================================

 if pltvort,
  cmin=-0.01; cmax=0.01; nbcol=10;
 else
  cmin=-0.25; cmax=0.25; nbcol=10;
 end
 cint=(cmax-cmin)/nbcol;
 map=colormap(cool(nbcol));
 colormap(map);
 set(gcf,'color','w');
 hold on

 if pltvort,
  [C,hh]=contourf(x,y,vort,[cmin:cint:cmax],'linestyle','none');
 else
  [C,hh]=contourf(x,y,speed,[cmin:cint:cmax],'linestyle','none');
 end
 colorbar;

 I = ~(sqrt(ur.^2+vr.^2)<0.05);
 quiver(x(I),y(I),ur(I),vr(I)); hold off
 axis([100 760 10 760])
 caxis([cmin cmax])
 thour=floor(time*24);
 title(['Rip Test Case - time = ',num2str(thour),' h'],'fontsize',14)
 xlabel('Cross-shore distance [m]','fontsize',14)
 ylabel('Along-shore distance [m]','fontsize',14)
 set(gca,'fontsize',14)
 hold off

 if makemovie,  
  % Write each frame to the file
  movObj.FrameRate =  5;
  writeMovie(movObj,getframe(hf));
  clf('reset')
 end

end % time loop -------------------------------------------------------

close(nc);

if makemovie,  
    movObj.Loop = 'loop'; % Set looping flag
    close(movObj);        % Finish writing movie and close file
end

if makepdf
 export_fig -transparent rip.pdf
end

return





