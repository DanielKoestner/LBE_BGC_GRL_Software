%% Use satellite Chla estimates to compare inside vs. outside
% 1. determine average CHL data for each month 2010–2022 within bounds
% 2. derive spatial weighted mean and standard deviation for in and out
% eddy using eddy heatmap to derive in and out eddy CHL per month
% 3. two panel plot, errobar of CHL per month and second panel for E

close all
clear all
curdir=cd;
addpath(genpath([curdir '/Scripts']))
addpath(genpath([curdir '/Data']))

print_flag=1;
save_flag=0;
scale_flag='A'; % A for 1.5r and B for 2r scaling of LBEZ

fs=14;
lw=1.5;
set(0, 'DefaultAxesFontName', 'Times');
set(0, 'DefaultTextFontName', 'Times');
alp=0.65;
mksz=7;
%% load data
load('LBE_Average_Eddy.mat')
load('LB_CCI_CHL_2010_2022_27-Dec-2024.mat');

% load('LBE_BGC_POC_2010_2022_MonthlyProfiles_30-Dec-2024','CHL*','Z')
load(['LBE_BGC_POC_2010_2022_MonthlyProfiles_' scale_flag '_06-Feb-2025'],'CHL*','Z')

%% Average monthly satellite data for region

% need to fix to create an appropriate counter for each pixel, i.e., how
% many actual "januarys" are available for each pixel 
for i = 1:12
    inds=find(Mn==i & Yr<2023);
    CHL_mn{i}=zeros(480,1080);
    cnt{i}=zeros(480,1080);
    for ii= 1:length(inds)
        tmp=CHL{inds(ii)}';
        cnt{i}=cnt{i}+~isnan(tmp);
        tmp(isnan(tmp))=0;
        CHL_mn{i}=CHL_mn{i}+tmp;
        clear tmp
    end
    CHL_mn{i}=CHL_mn{i}./cnt{i};
end


% interpolate to LB bounds used for floats
for i = 1:12
CHL_mn_LB{i}=interp2(Plg,Plt,CHL_mn{i},X_lb,Y_lb);
end

% indices 1–3 is average of upper 15 m
for i = 1:12
    chl_lbe(i)=mean(CHLA_in(1:3,i));
    chl_lbe_sd(i)=mean(CHLA_in_sd(1:3,i));
    chl_lbe_ci(i)=(mean(CHLA_in_hi(1:3,i)) - mean(CHLA_in_low(1:3,i)))/2;

    chl_lb(i)=mean(CHLA_out(1:3,i));
    chl_lb_sd(i)=mean(CHLA_out_sd(1:3,i));
    chl_lb_ci(i)=(mean(CHLA_out_hi(1:3,i)) - mean(CHLA_out_low(1:3,i)))/2;
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
    
    tmp=CHL_mn_LB{i};
    tmp(isnan(tmp))=0;
    CHL_lbe(i)=sum(sum(tmp.*wt_LBE))./sum(sum(wt_LBE));
    CHL_lbe_sd(i)=sqrt(sum(sum(wt_LBE.*((tmp-CHL_lbe(i))).^2))/(sum(sum(wt_LBE))));
    clear tmp

    tmp=CHL_mn_LB{i};
    tmp(isnan(tmp))=0;
    CHL_lb(i)=sum(sum(tmp.*wt_LB))./sum(sum(wt_LB));
    CHL_lb_sd(i)=sqrt(sum(sum(wt_LB.*((tmp-CHL_lb(i))).^2))/(sum(sum(wt_LB))));
    clear tmp
end



%% Plot
colors2=crameri('batlow',8);
colors(1,:)=colors2(7,:);
colors(2,:)=colors2(3,:);
clear colors2

hf=figure();
set(hf,'Units','inches','Position', [5 5 8 3], 'PaperPosition', [0 0 8 3], 'PaperSize', [8 3]);
ha1=iSubplot(1,1, 'Gap', [0 0.05], 'Min', [0.055 0.1], 'Max', [0.99 0.96], 'XTickL', 'All', 'YTickL', 'All');

axes(ha1(1))

hold on

p3=errorbar((1:12)-.2,CHL_lbe,CHL_lbe_sd,'color',colors(1,:),'linestyle','none','marker','*','MarkerFaceColor',colors(1,:),'MarkerSize',5,'CapSize',0);
p4=errorbar((1:12)+.2,CHL_lb,CHL_lb_sd,'color',colors(2,:),'linestyle','none','marker','*','MarkerFaceColor',colors(2,:),'MarkerSize',5,'CapSize',0);

p1=errorbar((1:12)-.1,chl_lbe,chl_lbe_ci,'color',colors(1,:),'linestyle','none','marker','d','MarkerFaceColor',colors(1,:),'MarkerSize',5);
p2=errorbar((1:12)+.1,chl_lb,chl_lb_ci,'color',colors(2,:),'linestyle','none','marker','d','MarkerFaceColor',colors(2,:),'MarkerSize',5);


set(gca,'ticklength',[0.01 0.005])
box on
grid on
ha1(1).YMinorTick='on';
% set(gca,'ytick',[0:0.5:2])
set(gca,'fontsize',fs-2)
xlim([0.5 12.5])
ylim([0 3])
set(gca,'yscale','lin')
ylabel('Chla [mg m^{-3}]')
set(gca,'xtick',[1:12],'xticklabel',{''});


[hl,lines]=legendflex([p1 p3 p2 p4],{'inside BGC-Argo','inside OC CCI','outside BGC-Argo','outside OC CCI'}, 'ref', gca, ...
'anchor', {'nw', 'nw'}, 'buffer', [5 -2], 'xscale', 0.5, 'ncol',2,'nrow',2, ...
'box', 'off', 'FontSize', fs-3);



set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});

%% print

if print_flag==1
    print(['Figures/V8/SI/LB22_SatCHL_' scale_flag],'-dpdf','-r800')
end
