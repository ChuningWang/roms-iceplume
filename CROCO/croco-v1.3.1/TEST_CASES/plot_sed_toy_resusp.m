%======================================================================
%
%     ---               Sed toy resuspension Test Case                ---     
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

fname='sed_toy_resusp_his.nc';

%======================================================================
%
% Read data
%

idy=1;
idx=1;
depth=20;       % 20m water column depth
beddepth=0.041; %  4cm sediment depth
it=61;          % 5days

nc=netcdf(fname);

h=squeeze(nc{'h'}(idy,idx));
T=squeeze(nc{'scrum_time'}(:,:))/86400;
[nt,t0]=size(T);
bostr=squeeze(nc{'bostr'}(:,idy,idx));
ALT=squeeze(nc{'act_thick'}(:,idy,idx));

s1=squeeze(nc{'sand_01'}(:,idy,idx));
s2=squeeze(nc{'sand_02'}(:,idy,idx));

m1=squeeze(nc{'mud_01'}(:,idy,idx));
m2=squeeze(nc{'mud_02'}(:,idy,idx));

sand1=s1*depth; % Kg/m2
sand2=s2*depth;
mud1=m1*depth;
mud2=m2*depth;


% for stratigraphy
fracs1=squeeze(nc{'bed_frac_sand_01'}(it,:,idy,idx));
fracs2=squeeze(nc{'bed_frac_sand_02'}(it,:,idy,idx));

fracm1=squeeze(nc{'bed_frac_mud_01'}(it,:,idy,idx));
fracm2=squeeze(nc{'bed_frac_mud_02'}(it,:,idy,idx));

bedthick=squeeze(nc{'bed_thick'}(it,:,idy,idx));
NL=size(bedthick,1);

zbed=[-beddepth;-beddepth+cumsum(flipud(bedthick(2:end,:)),1)];
zbed=flipud(zbed)*100; %cm

close(nc);

Yconc=zeros(4,nt); % 4 : 2 sands + 2 muds
for i=1:nt
  conc=[sand1(i) sand2(i) mud1(i) mud2(i)];
  Yconc(:,i)=conc;
end

nl=size(fracs1,1);
Yfrac=zeros(4,nl); % 4 : 2 sands + 2 muds
for k=1:nl
  frac=[fracs1(k) fracs2(k) fracm1(k) fracm2(k)];
  Yfrac(:,k)=frac;
end

%----------------------------------------------------------
%  Plot
%----------------------------------------------------------
figure('Position',[1 1 1050 800])

title('Two successive erosion–deposition events lasting 5 days')

subplot(3,2,1)
% Mass of sediment in suspension
area(T,transpose(Yconc))
ylabel('Mass of sediment (Kg/m2)','fontsize',12);
legend({'Sand1','Sand2 ', 'Mud1', 'Mud2'})
hold off
grid on

subplot(3,2,3)
% bostr
plot(T,bostr,'linewidth',1)
ylabel('Bottom stress (Pa)','fontsize',12);
hold off
grid on

subplot(3,2,5)
% Active Layer Thickness (cm)
plot(T,ALT*100,'linewidth',1)
ylabel('Active Layer thickness (cm)','fontsize',12);
hold off
grid on

subplot(3,2,[2,4,6])
% stratigraphy
width=1;
%barh(zbed(:,1),Yfrac,width,'stacked')
area(zbed(:,1),transpose(Yfrac))
camroll(-270)
set(gca,'YDir','reverse');
xline(zbed(:,1),'--')
xlim([-4.1 0])
ylim([0 1])

xlabel('Z (cm)','fontsize',12);
ylabel([ 'fraction of sediment ' , num2str(round(T(it))), ' Days'] ,'fontsize',12);
hold off
grid on

set(gcf,'PaperPositionMode','auto');


export_fig -transparent sed_toy_resusp.pdf


return

