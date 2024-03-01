%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make 1 plot from the results of the SOLITON test case
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
%  Copyright (c) 2002-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%====================================================================
% This test problem considers the propagation of a Rossby soliton on 
% an equatorial beta-plane, for which an asymptotic solution exists to
% the inviscid, nonlinear shallow water equations. In principle, the 
% soliton should propagate westwards at fixed phase speed, without 
% change of shape. Since the uniform propagation and shape preservation 
% of the soliton are achieved through a delicate balance between linear 
% wave dynamics and nonlinearity, this is a good context in which to 
% look for erroneous wave dispersion and/or numerical damping.
%
% Problem is nondimensionalized with:
% H = 40 cm, L=295 km, T = 1.71 days and U=L/T=1.981 m/s.
%
% Theorical propagation speed is 0.4 (0.395 actually) so that at t=120, 
% the soliton should be back to its initial position after crossing 
% the periodic canal of length 48.
%
% Reference: 
% Boyd J.P., 1980: Equatorial solitary waves. Part I: Rossby solitons.
% JPO, 10, 1699-1717
%====================================================================

tndx=16;
nc=netcdf('soliton_his.nc');
makepdf=0;

time=(nc{'scrum_time'}(tndx));
x=nc{'x_rho'}(:);
y=nc{'y_rho'}(:);
z1=squeeze(nc{'zeta'}(tndx,:,:));
z0=squeeze(nc{'zeta'}(1,:,:));
close(nc);

H = 1; % 40;    % cm 
L = 1; % 295;   % km
T = 1; % 1.71;  % days
U = 1; % 1.981; % m/s.
z0=z0*H;
z1=z1*H;
time=time*T;
x=x*L;
y=y*L;

figure
cnt=[-0.05:0.02:0.2]*H;
contourf(x,y,z1,cnt,'linestyle','none')
axis image
colorbar('h')
hold on
contour(x,y,z0,cnt,'k')
hold off
title(['SOLITON - zeta at t=',num2str(time),' [non dim]'])
set(gca,'fontsize',15)

if makepdf,
 export_fig -transparent soliton.pdf
end

maxzi=max(max(z0));
maxzf=max(max(z1));
Xi=mean(x(z0==maxzi));
Xf=mean(x(z1==maxzf));
X0=48*L;
X_croco=X0-(Xf-Xi);

disp(['Final apmplitude (vs. theory) : ', num2str(maxzf), ...
                                     ' (',num2str(maxzi),')'])
disp(['Final position   (vs. theory) : ', num2str(X_croco), ...
                                     ' (',num2str(X0),')'])


