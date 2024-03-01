%======================================================================
%
%     ---            Analytical Dune Test Case                ---     
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

% Choose usgs (1) or mustang (0) model output file
usgs=1
    
fname='ana_dune_his.nc';

j=2;
tndx=2; % half-hourly outputs

%======================================================================

%----------------------------------------------------------
% Read Model data
%----------------------------------------------------------
%
nc=netcdf(fname);
t0=nc{'scrum_time'}(:);

nt=size(t0);
t=zeros(nt);
for i=2:nt(1)
  t(i)=t0(i)-t0(1);
end
tndx=min(tndx,length(t));
disp([ 'tndx = ',num2str(tndx), ...
    ' - Time = ',num2str(t(tndx)/3600),' hr' ])
x=squeeze(nc{'x_rho'}(j,:));

if usgs ~= 1
    var_hmorph='Hm';
else
    var_hmorph='hmorph';
end
hm=squeeze(nc{var_hmorph}(:,j,:));
close(nc);

%----------------------------------------------------------
% Analytical solution
%----------------------------------------------------------
% From Marrieu Phd thesis 2007 and 
% Long, Wen, James T. Kirby, and Zhiyu Shao, A numerical scheme for
% morphological bed level calculations.
% Coastal Engineering, 2007, doi :10.1016/j.coastaleng.2007.09.009.

% Model parameters
alpha=0.001; gamma=0.01; beta=3;

h0=6; xc=150; Q=10; poros=0.4;

% Initial bathymetry
h= h0 - 2*exp(-gamma*(x-xc).^2);

% Initial phase speed of bedform 
%   dq/dh*1/(1-poros), with beload q=alpha*u^beta
C=alpha*beta*Q^beta./(-h).^(beta+1)./(1-poros);

% Set array of final bathy
ha=NaN.*ones(length(t),length(x));
ha(1,:)=h;

% Mass budget under bedform
sed_dune=sum(h0-h); 

for k=2:length(t)

 % Find location of bedform: hai et hbi
 % -----------------------------------
 % hai and hbi are the 2 possible solutions of the equation
 %
 % At any time, a location can originate from a position
 % located both upstream and downstream of the bedform; the peak 
 % travels faster and overtakes the downstream section of the dune;
 % we use mass conservation to choose the best of 2 solutions

 hai=NaN.*ones(length(x));
 hbi=NaN.*ones(length(x));

 % Position of first half of bedform at t
 %   d : position of given point on initial bedform
 %   hai : new bathy for the first half of bedform
 for i=1:ceil(length(x)/2)
  d(i)=x(i)+C(i)*t(k); 
  idx=find(x>d(i),1,'first');
  hai(find(x>d(i),1,'first'))=h(i); 
 end
 % Position of second half of bedform at t
 %   hbi : new bathy for the second half of bedform
 for i=ceil(length(x)/2)+1:length(x)
  d(i)=x(i)+C(i)*t(k);
  idx=find(x>d(i),1,'first');
  hbi(find(x>d(i),1,'first'))=h(i); 
 end

% Interpolate bathy
% -----------------
 idnan=isnan(hai);
 hai=interp1(x(~idnan),hai(~idnan),x);
 idnan=isnan(hbi);
 hbi=interp1(x(~idnan),hbi(~idnan),x);

 % Locate shock : where hai and hbi solutions should merge
 % -------------------------------------------------------
 % hai and hbi are the solutions before and after the shock

 % initialize with the first point of the second half of bedform
 id_shock=find(~isnan(hbi),1,'first')-1; 

 % look for the solution (hai or hbi) that maintain mass budget
 mass_sed=0;
 mass_sedp1=0;
 while (mass_sed<=sed_dune & mass_sed<=mass_sedp1)
  id_shock=id_shock+1;
  mass_sed  =nansum(h0-hai(1:id_shock-1))+ ...
             nansum(h0-hbi(id_shock  :end)); 
  mass_sedp1=nansum(h0-hai(1:id_shock  ))+ ...
             nansum(h0-hbi(id_shock+1:end));
 end

 ha(k,:)=[hai(1:id_shock-1) hbi(id_shock:end)]; % final solution 

end

%----------------------------------------------------------
%  Plot analytical and numerical bed evolution 
%----------------------------------------------------------
figure('Position',[1 1 920 550])
hold on
for n=1:3:length(t)/2
 line(x,-ha(n,:),'color','r','Linewidth',2)
 line(x,-hm(n,:),'color','k','Linewidth',2)
end
text(120,-3.6,'Numerical', ...
              'color','k','fontsize',20)
text(120,-3.8,'Analytical', ...
	      'color','r','fontsize',20)
hold off
grid on
axis([120 200 -7 -3])
if usgs ~= 1
   title(['ANA-DUNE Test Case  (MUSTANG) - Bed evolution' ])
else
   title(['ANA-DUNE Test Case  (USGS) - Bed evolution'])
end


set(gca,'fontsize',15);
set(gcf,'PaperPositionMode','auto');

if usgs ~= 1
 export_fig -transparent ana_dune_mustang.pdf
else
 export_fig -transparent ana_dune_usgs.pdf
end

return

