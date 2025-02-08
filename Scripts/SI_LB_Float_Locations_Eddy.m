%% set up

% draw a polar plot with relative float locations to eddy center
% UPDATE ONLY FOR QC data...
% look into histogram or other statistical look at dist


close all
clear all

curdir=cd;

addpath(genpath([curdir '/Scripts']))
addpath(genpath([curdir '/Data']))

addpath(genpath([curdir '/aux']))

print_flag=1;


scale_flag='A'; % A for 1.5r and B for 2r scaling of LBEZ


fs=17;
lw=1.5;
set(0, 'DefaultAxesFontName', 'Times');
set(0, 'DefaultTextFontName', 'Times');
alp=0.65;
mksz=7;

% month colors
c=crameri('roma',12);
colorsm=[c(9:12,:); c(1:8,:)];

%% data

% additional qc data, what is used to compute iPOC data/fluxes/TE
load('LBE_BGC_POC_2010_2022_MonthlyParameters_A_06-Feb-2025_200.mat','R','TH','SP','days'); % July 30 used 200 m, and trap integration
mns=month(days);



for i = 1:length(R)
    polarplot(TH(i),log10(R(i)),'o','MarkerFaceColor',colorsm(mns(i),:),'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',SP(i)*18 + 2);
    hold on
end

thetaticks(0:45:360);
thetaticklabels({'E','', 'N','', 'W','', 'S',''})
ax = gca;
rlim([0 log10(200)])
rticks([log10(10) log10(50) log10(100) log10(200)])
rticklabels({'10','50','100','200'})

hf=gcf;

set(hf,'Units','inches','Position', [5 5 5 5], 'PaperPosition', [0 0 5 5], 'PaperSize', [5 5]);

%% legend

for i = 1:12
    s(i)=polarplot(0,300,'o','MarkerFaceColor',colorsm(i,:),'MarkerEdgeColor',[0.5 0.5 0.5]);
end

lgd_text={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
[hl,lines]=legendflex([s],lgd_text, 'ref', gca, ...
'anchor', {'s', 's'}, 'buffer', [0 -40], 'xscale', 0.5, 'ncol',6,'nrow',2, ...
'box', 'off', 'FontSize', fs-6);


s2(1)=polarplot(0,300,'o','MarkerFaceColor','none','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',(0.11)*18 + 2);
s2(2)=polarplot(0,300,'o','MarkerFaceColor','none','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',(0.46)*18 + 2);

[hl,lines]=legendflex([s2],{'0.11 m/s','0.46 m/s'}, 'ref', gca, ...
'anchor', {'ne', 'ne'}, 'buffer', [20 15], 'xscale', 0.5, 'ncol',1,'nrow',2, ...
'box', 'off', 'FontSize', fs-6);

text(0.01,1.03,'(a)','units','normalized')
%% violins

for i = 1:12
    inds=ismember(mns,i);
    r{i}=R(inds);
    sp{i}=SP(inds);
    
end

hf2=figure();
set(hf2,'Units','inches','Position', [5 5 5 5], 'PaperPosition', [0 0 5 5], 'PaperSize', [5 5]);
ha1=iSubplot(2,1, 'Gap', [0 0.03], 'Min', [0.08 0.05], 'Max', [.98 .96], 'XTickL', 'All', 'YTickL', 'All');

axes(ha1(1));
[h,L,MX,MED,bw]=violin(r,'x',[1:12],'facecolor',[0.993090729481532	0.692220292760577	0.667205149785594],'edgecolor',[0.6 0.6 0.6],'mc','k','medc','k-.');
ylim([0 200])
yticks([0:25:200])
set(gca,'xtick',[1:12],'xticklabel',{''});
ylabel('distance to eddy center [km]')
set(gca,'XTickLabelRotation',45)
grid on
text(0.01,1.04,'(b)','units','normalized')

axes(ha1(2));
[h,L,MX,MED,bw]=violin(sp,'x',[1:12],'facecolor',[0.993090729481532	0.692220292760577	0.667205149785594],'edgecolor',[0.6 0.6 0.6],'mc','k','medc','k-.');
ylim([0.1 0.5])
% yticks([0:25:200])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
ylabel('average rotational speed [m/s]')
set(gca,'XTickLabelRotation',45)
grid on
text(0.01,1.04,'(c)','units','normalized')
L.Visible='off';
%% print

if print_flag==1
    figure(hf)
    print(['Figures/V8/SI/LB22_Float_Locations_Eddy_' scale_flag],'-dpdf','-r800')
    figure(hf2)
    print(['Figures/V8/SI/LB22_Float_Locations_Eddy_Violin_' scale_flag],'-dpdf','-r800')
end