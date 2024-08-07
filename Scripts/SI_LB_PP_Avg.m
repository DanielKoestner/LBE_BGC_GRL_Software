%% Use satellite PP estimates to compare inside vs. outside
% 1. determine average PP data for each month 2010–2022 within bounds
% 2. derive spatial weighted mean and standard deviation for in and out
% eddy using eddy heatmap to derive in and out eddy PP per month
% 3. two panel plot, errobar of PP per month and second panel for E

close all
clear all
curdir=cd;
addpath(genpath([curdir '/Scripts']))
addpath(genpath([curdir '/Data']))

print_flag=1;
save_flag=0;
fs=14;
lw=1.5;
set(0, 'DefaultAxesFontName', 'Times');
set(0, 'DefaultTextFontName', 'Times');
alp=0.65;
mksz=7;
%% load data
load('LBE_Average_Eddy.mat')
load('LB_GlobColour_PP_2010_2024_15-Apr-2024.mat');

%% Average monthly satellite data for region
for i = 1:12
    inds=find(Mn==i & Yr<2023);
    PP_mn{i}=zeros(480,1080);
    for ii= 1:length(inds)
        tmp=PP{inds(ii)}';
        PP_mn{i}=PP_mn{i}+tmp;
        clear tmp
    end
    PP_mn{i}=PP_mn{i}/ii;
end


% interpolate to LB bounds used for floats
for i = 1:12
PP_mn_LB{i}=interp2(Plg,Plt,PP_mn{i},X_lb,Y_lb);
end

%% eddy weighting

EDDY=EDDY-10;
% wt_LBE=EDDY/max(max(EDDY));
% wt_LB=1-wt_LBE;

eddy=interp2(X_eddy,Y_eddy,EDDY,X_lb,Y_lb);
wt_LBE=eddy/max(max(eddy));
wt_LBE(isnan(wt_LBE))=0;
wt_LB=1-wt_LBE;

for i = 1:12
    
    tmp=PP_mn_LB{i};
    tmp(isnan(tmp))=0;
    pp_lbe(i)=sum(sum(tmp.*wt_LBE))./sum(sum(wt_LBE));
    pp_lbe_sd(i)=sqrt(sum(sum(wt_LBE.*((tmp-pp_lbe(i))).^2))/(sum(sum(wt_LBE))));
    clear tmp

    tmp=PP_mn_LB{i};
    tmp(isnan(tmp))=0;
    pp_lb(i)=sum(sum(tmp.*wt_LB))./sum(sum(wt_LB));
    pp_lb_sd(i)=sqrt(sum(sum(wt_LB.*((tmp-pp_lb(i))).^2))/(sum(sum(wt_LB))));
    clear tmp
end


% Flux (E) in difference between months (units are already average Flux?)
P_in=diff(pp_lbe); P_in(12)=0;
P_out=diff(pp_lb); P_out(12)=0;

%% Plot
colors2=crameri('batlow',8);
colors(1,:)=colors2(7,:);
colors(2,:)=colors2(3,:);
clear colors2

hf=figure();
set(hf,'Units','inches','Position', [5 5 8 2], 'PaperPosition', [0 0 8 2], 'PaperSize', [8 2]);
ha1=iSubplot(1,1, 'Gap', [0 0.05], 'Min', [0.055 0.1], 'Max', [0.99 0.96], 'XTickL', 'All', 'YTickL', 'All');

axes(ha1(1))

hold on

p1=errorbar((1:12)-.1,pp_lbe,pp_lbe_sd,'color',colors(1,:),'linestyle','none','marker','o','MarkerFaceColor',colors(1,:),'MarkerSize',5);
p2=errorbar((1:12)+.1,pp_lb,pp_lb_sd,'color',colors(2,:),'linestyle','none','marker','s','MarkerFaceColor',colors(2,:),'MarkerSize',5);

set(gca,'ticklength',[0.01 0.005])
box on
grid on
ha1(1).YMinorTick='on';
% set(gca,'ytick',[0:0.5:2])
set(gca,'fontsize',fs-2)
xlim([0 13])
ylim([0 700])
ylabel('PP [mg m^{-2} d^{-1}]')
set(gca,'xtick',[1:12],'xticklabel',{''});


[hl,lines]=legendflex([p1 p2],{'inside','outside'}, 'ref', gca, ...
'anchor', {'nw', 'nw'}, 'buffer', [5 -2], 'xscale', 0.5, 'ncol',1,'nrow',3, ...
'box', 'off', 'FontSize', fs-3);



% axes(ha1(2))
% 
% hold on
% 
% p1=plot((1:12),P_in,'color',colors(1,:),'linestyle','--','marker','o','MarkerFaceColor',colors(1,:),'MarkerSize',5);
% p2=plot((1:12),P_out,'color',colors(2,:),'linestyle','--','marker','s','MarkerFaceColor',colors(2,:),'MarkerSize',5);
% 
% set(gca,'ticklength',[0.01 0.005])
% box on
% grid on
% ha1(2).YMinorTick='on';
% % set(gca,'ytick',[0:0.5:2])
% set(gca,'fontsize',fs-2)
% xlim([0 13])
% ylim([-300 300])
% set(gca,'ytick',(-300:150:300))
% ylabel('\itE\rm [mg m^{-2} d^{-1}]')
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});

%% print

if print_flag==1
    print('Figures/V3/SI/LB22_SatPP','-dpdf','-r800')
end
