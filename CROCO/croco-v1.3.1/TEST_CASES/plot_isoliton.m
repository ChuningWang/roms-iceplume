%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make 1 plot from the results of the GRAV_ADJ test case
% 
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
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
%  Copyright (c) 2002-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
fname     = 'isoliton_his.nc';   % croco output file
tindex    = [30;50;70];          % record indices to plot
makepdf   = 0;                   % make pdf file
%
%======================================================================
nplot=length(tindex);
hFig = figure;
colormap('jet')
set(hFig, 'Position', [200 500 700 150*nplot])

for it=1:length(tindex); % --------------- time loop

tndx=tindex(it);

% -------------------------------------
%  --- Read model data ---
% -------------------------------------
nc=netcdf(fname,'r');
tndx=min(tndx,length(nc{'scrum_time'}(:)));
disp(['tndx = ',num2str(tndx)'']);
h=nc{'h'}(:);
x=squeeze(nc{'x_rho'}(2,:));
zeta=squeeze(nc{'zeta'}(tndx,:,:));
t=squeeze(nc{'rho'}(tndx,:,2,:));
w=1000*squeeze(nc{'w'}(tndx,:,2,:));
[N,M]=size(t);
theta_s=nc.theta_s(:);
theta_b=nc.theta_b(:);
hc=nc.hc(:);
close(nc);

zr = zlevs(h,zeta,theta_s,theta_b,hc,N,'r',2);
zr=squeeze(zr(:,2,:));
xr=reshape(x,1,M);
xr=repmat(xr,[N 1]);
t(t==0)=NaN;

% -------------------------------------
%  --- Plot ---
% -------------------------------------

subplot(length(tindex),1,it)
[C,h] = contourf(xr,zr,t,[-40:1:40]); hold on
set(h,'LineColor','none')
xlabel('X [m]')
ylabel('Depth [m]')
axis([-Inf Inf -0.25 -0.15])
caxis([20 30])
colorbar
if it==1,
 title('Internal Soliton')
end

end % ----------------------- time loop

if makepdf
 export_fig -transparent -pdf isoliton.pdf
end




