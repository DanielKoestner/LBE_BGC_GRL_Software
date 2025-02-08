% LBE analysis, must extract only dates when we have floats.


close all
clear all
curdir=cd;

curdir=cd;

addpath(genpath([curdir '/Scripts']))
addpath(genpath([curdir '/Data']))
addpath(genpath([curdir '/OneArgo']))
addpath(genpath([curdir '/aux']))

print_flag=1;
save_flag=0;
fs=14;
lw=1.5;
set(0, 'DefaultAxesFontName', 'Times');
set(0, 'DefaultTextFontName', 'Times');


%% Get data
S=ncread('META3.2_DT_allsat_Anticyclonic_long_LBE.nc','speed_average');
R=ncread('META3.2_DT_allsat_Anticyclonic_long_LBE.nc','effective_radius');
R=R/1000;
A=ncread('META3.2_DT_allsat_Anticyclonic_long_LBE.nc','effective_area');
r=sqrt(A/pi)/1000; %effective area converted to equivalent circular radius in km

temp=ncread('META3.2_DT_allsat_Anticyclonic_long_LBE.nc','time');
%convert to DD-MM-YYYY
LBE_dates=datetime(temp*24*60*60,'ConvertFrom','epochtime','Epoch','1950-01-01')';
LBE_dnum=datenum(LBE_dates);


% get float dates

load('LBE_Average_Eddy.mat','dnums');

for i = 1:length(dnums)
% find closest eddy location date
[~,ind]=min(abs(dnums(i)-LBE_dnum));

r_lb(i)=R(ind)*1.5;
s_lb(i)=S(ind);
clear ind
end


%% plot

hf=figure();
set(hf,'Units','inches','Position', [5 5 6.5 4], 'PaperPosition', [0 0 6.5 4], 'PaperSize', [6.5 4]);
ha1=iSubplot(1,2, 'Gap', [0.02 0.04], 'Min', [0.07 0.105], 'Max', [0.98 0.96], 'XTickL', 'All', 'YTickL', 'All');

axes(ha1(1))
hold on
grid on
xlabel('Effective Radius [km]')
ylabel('counts')
set(gca,'yminortick','on')
ylim([0 150]);
xlim([20 160]);
set(gca,'fontsize',fs-2)
set(gca,'ticklength',[0.025 0.025])
set(gca,'xtick',[0:25:150])
psr=histfit(r_lb);
psr(1).FaceColor=[0.7 0.7 0.9];
psr(2).Color=[0.3 0.3 0.3];
text(0.05,0.95,'(a)','units','normalized')
text(0.76,0.95,sprintf('\\mu = %1.0f km',mean(r_lb)),'units','normalized')
text(0.76,0.9,sprintf('\\sigma = %1.0f km',std(r_lb)),'units','normalized')

axes(ha1(2))
hold on
grid on
ylim([0 150]);
xlim([0 0.6]);
set(gca,'yticklabel',{''})
set(gca,'yminortick','on')
xlabel('Average Rotational Speed [m/s]')
set(gca,'fontsize',fs-2)
set(gca,'ticklength',[0.025 0.025])
set(gca,'xtick',[0:.15:0.6])
pss=histfit(s_lb);
pss(1).FaceColor=[0.7 0.7 0.9];
pss(2).Color=[0.3 0.3 0.3];
text(0.05,0.95,'(b)','units','normalized')
text(0.71,0.95,sprintf('\\mu = %0.2f m/s',mean(s_lb)),'units','normalized')
text(0.71,0.9,sprintf('\\sigma = %0.2f m/s',std(s_lb)),'units','normalized')


%% print

if print_flag==1
    print('Figures/V8/SI/LB22_EddyStats','-dpdf','-r800')
end
