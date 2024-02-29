%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make plot for Kelvin-Helmholtz Instability test case (KH_INST)
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
warning off
%================== User defined parameters ===========================
%
% --- model params ---
%
fname     = 'khinst_his.nc';  % croco file name

makemovie = 0;             % make movie using QTWriter
makepdf   = 0;             % make pdf file

%======================================================================

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

nc=netcdf(fname,'r');
tindex=length(nc{'scrum_time'}(:)); % reads last record

if makemovie,
 tstr=1;
 tend=tindex;
 movObj = QTWriter('khinst.mov');
else
 tstr=tindex;
 tend=tstr;
end

hf = figure('position',[500 500 700 500]);
set(gca,'FontSize',15)
axis tight; set(hf,'DoubleBuffer','on');
%set(gca,'nextplot','replacechildren');

cmin=15.5; cmax=18.5; nbcol=100;

for tindex=tstr:tend % ---------------------------------------------

 if tindex==tstr
  h2d=nc{'h'}(:,:);
  h=squeeze(nc{'h'}(2,:));
  xl=nc{'xl'}(:);
  x=squeeze(nc{'x_rho'}(2,:));
  pm=squeeze(nc{'pm'}(2,:));
  N=length(nc('s_rho'));
  zeta2d=squeeze(nc{'zeta'}(tindex,:,:));
  rho=squeeze(nc{'rho'}(tindex,:,2,:));

  [N,M]=size(rho);
  theta_s=nc.theta_s(:);
  theta_b=nc.theta_b(:);
  hc=nc.hc(:);
  z=zlevs(h2d,zeta2d,theta_s,theta_b,hc,N,'r',2);
  z=squeeze(z(:,2,:));
  x=reshape(x,1,M);
  x=repmat(x,[N 1]);
 end

 time=nc{'scrum_time'}(tindex)/60;
 rho=squeeze(nc{'rho'}(tindex,:,2,:));

 %============================================================
 % --- plot ---
 %=============================================================
 cint=(cmax-cmin)/nbcol;
 colormap(jet)
 set(gcf,'color','w');
 set(gca,'fontsize',15);

 [C,hh]=contourf(x,z,rho,[cmin:cint:cmax],'LineStyle','none'); 
 %pcolor(x,z,rho); shading flat
 shading flat; colorbar;
 axis([0 256 -256 0])
 caxis([cmin cmax])
 thour=floor(time/60);
 tmin=floor(time)-60*thour;
 clock=[num2str(thour),' h ',num2str(tmin),' min '];
 title(['Time: ',clock],'fontsize',15)
 xlabel('X [m]')
 ylabel('Z [m]')
 disp(['Plot for time: ',clock,' (Record: ',num2str(tindex),')'])
 set(gcf,'PaperPositionMode','auto');

 if makemovie,  
  % Write each frame to the file
  movObj.FrameRate =  5;
  writeMovie(movObj,getframe(hf));
  clf('reset')
 end

%----------------------------------

end % time loop

close(nc);

if makemovie,  
    movObj.Loop = 'loop'; % Set looping flag
    close(movObj);        % Finish writing movie and close file
end

if makepdf, 
    export_fig -transparent khinst.pdf
end

return





