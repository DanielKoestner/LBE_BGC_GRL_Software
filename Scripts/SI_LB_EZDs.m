%% Comparison of EZDs and MLDs

% violin plots of three EZD metrics and MLD (only for when PAR exists)
% or violin plots of all metrics (regardless of PAR existence) (by month?)

%scatter plot of EZD(par) and EZD_chla w/ slope and 1:1 line, cor



%% set up

close all
clear all

curdir=cd;

addpath(genpath([curdir '/Scripts']))
addpath(genpath([curdir '/Data']))

addpath(genpath([curdir '/aux']))

print_flag=1;

scale_flag='A'; % A for 1.5r and B for 2r scaling of LBEZ


fs=16;
lw=1.5;
set(0, 'DefaultAxesFontName', 'Times');
set(0, 'DefaultTextFontName', 'Times');
alp=0.65;
mksz=7;





mk{1}='s';
mk{2}='^';
mk{3}='d';
mk{4}='v';
mk{5}='p';
mk{6}='h';
mk{7}='>';
mk{8}='<';

% month colors
c=crameri('roma',12);
colorsm=[c(9:12,:); c(1:8,:)];

warning off

thrs=0.702; %exclude profiles with less than 60 5-m bin resolved measures

%% data
load(['LBE_BGC_POC_2010_2022_' scale_flag '_06-Feb-2025.mat'],'in','out','lb'); 

colors2=crameri('batlow',8);
colsin=colors2(7,:);
colsout=colors2(3,:);
sat=0.4;

% in
ezd=[];ezd_c=[];ezd_iso=[]; ezd_p=[];
mld=[];
zprod=[];
mns=[]; %month for color?
for i = 1:length(in)
    inds=sum(isnan(in{i}.zbin.poc_s))<(thrs)*length(in{i}.zbin.z);
    ezd=[ezd in{i}.ezd(inds)];
    ezd_p=[ezd_p in{i}.ezd_p(inds)];
    ezd_c=[ezd_c in{i}.ezd_c(inds)];
    ezd_iso=[ezd_iso in{i}.ezd_iso(inds)];
    mld=[mld in{i}.mld(inds)];
    zprod=[zprod in{i}.zprod(inds)];
    mns=[mns month(datetime(in{i}.dnum(inds),'ConvertFrom','datenum'))];
end


for i = 1:12
    inds=ismember(mns,i) & zprod<900;
    EZD{i}=ezd(inds);
    EZD_c{i}=ezd_c(inds);
    EZD_iso{i}=ezd_iso(inds);
    EZD_p{i}=ezd_p(inds);
    MLD{i}=mld(inds);
end



%out

% in
ezd2=[];ezd2_c=[];ezd2_iso=[]; ezd2_p=[];
mld2=[];
zprod2=[];
mns2=[]; %month for color?
for i = 1:length(out)
    inds=sum(isnan(out{i}.zbin.poc_s))<(thrs)*length(out{i}.zbin.z);
    ezd2=[ezd2 out{i}.ezd(inds)];
    ezd2_p=[ezd2_p out{i}.ezd_p(inds)];
    ezd2_c=[ezd2_c out{i}.ezd_c(inds)];
    ezd2_iso=[ezd2_iso out{i}.ezd_iso(inds)];
    mld2=[mld2 out{i}.mld(inds)];
    zprod2=[zprod2 out{i}.zprod(inds)];
    mns2=[mns2 month(datetime(out{i}.dnum(inds),'ConvertFrom','datenum'))];
end


for i = 1:12
    inds=ismember(mns2,i) & zprod2<900;
    EZD2{i}=ezd2(inds);
    EZD2_c{i}=ezd2_c(inds);
    EZD2_iso{i}=ezd2_iso(inds);
    EZD2_p{i}=ezd2_p(inds);
    MLD2{i}=mld2(inds);

    [~,pezd(i),~,~]=ttest2(EZD2{i},EZD{i});
    [~,pmld(i),~,~]=ttest2(MLD2{i},MLD{i});
end



%% Start with monthly violins for (max) EZD and MLD

cols=crameri('imola',5);


%
hf1=figure();
set(hf1,'Units','inches','Position', [5 5 9 6], 'PaperPosition', [0 0 9 6], 'PaperSize', [9 6]);
ha=iSubplot(2,1, 'Gap', [0 0.015], 'Min', [0.07 0.04], 'Max', [.98 .97], 'XTickL', 'All', 'YTickL', 'All');


axes(ha(1));
text(-0.1,1,'(a)','units','normalized','fontsize',fs-3)
grid on
plot([3.5 3.5],[0 250],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([6.5 6.5],[0 250],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([9.5 9.5],[0 250],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
% plot([0.5 12.5],[0.5 0.5],'k:','linewidth',1.25)
hold on

% text(1,1.04,'EZD','units','normalized','fontsize',fs-4,'fontweight','bold','HorizontalAlignment','right')

[h2,L2,MX2,MED2,bw2]=violin(EZD2,'x',[1.2:1:12.2],'facecolor',colsout,'edgecolor',[0.8 0.8 0.8],'mc','k','medc','k-.');
[h,L,MX,MED,bw]=violin(EZD,'x',[0.8:1:11.8],'facecolor',colsin,'edgecolor',[0.8 0.8 0.8],'mc','k','medc','k-.');

L2.Visible='off';
L.Visible='off';
ylim([0 250])

text(0.01,0.94,'inside LBEZ','units','normalized','fontsize',fs-4,'FontWeight','bold','Color',colsin);
text(0.01,0.87,'outside LBEZ','units','normalized','fontsize',fs-4,'FontWeight','bold','Color',colsout);
% text(0.01,0.75,'Lower Twilight Zone','units','normalized','fontsize',fs-4,'FontWeight','bold','Color',coltz2);
%
w=1-pezd;

%alpha weighting by significance of difference in mean
for i = 1:12
    if pezd(i)<0.05
        h(i).FaceAlpha=0.8;
        h2(i).FaceAlpha=0.8;
    else
        h(i).FaceAlpha=0.4;
        h2(i).FaceAlpha=0.4;
    end
end



set(gca,'fontsize',fs)
set(gca,'ticklength',[0.01 0.01])
set(gca,'xtick',1:12)
set(gca,'xticklabel',{''})
set(gca,'xlim',[0.5 12.5])
set(gca,'yminortick','on')
% ylabel('iPOC^{sub}\it_s\rm or iPOC^{mes}\it_s\rm / iPOC\it_s\rm');
ylabel('EZD [m]');

% AX(1).YColor=colors2(2,:);
grid on

axes(ha(2));
text(-0.1,1,'(b)','units','normalized','fontsize',fs-3)

grid on
plot([3.5 3.5],[0 1000],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([6.5 6.5],[0 1000],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([9.5 9.5],[0 1000],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
% plot([0.5 12.5],[0.5 0.5],'k:','linewidth',1.25)
hold on

% text(1,1.04,'MLD','units','normalized','fontsize',fs-4,'fontweight','bold','HorizontalAlignment','right')

[h2,L2,MX2,MED2,bw2]=violin(MLD2,'x',[1.2:1:12.2],'facecolor',colsout,'edgecolor',[0.8 0.8 0.8],'mc','k','medc','k-.');
[h,L,MX,MED,bw]=violin(MLD,'x',[0.8:1:11.8],'facecolor',colsin,'edgecolor',[0.8 0.8 0.8],'mc','k','medc','k-.');
% L2.Visible='off';

ylim([0 1000])
% L2.Visible='off';
% text(0.01,0.95,'Productive Zone','units','normalized','fontsize',fs-4,'FontWeight','bold','Color',colpz);
% text(0.01,0.85,'Twilight Zone','units','normalized','fontsize',fs-4,'FontWeight','bold','Color',coltz);
% text(0.01,0.75,'Lower Twilight Zone','units','normalized','fontsize',fs-4,'FontWeight','bold','Color',coltz2);
%
w=1-pezd;

%alpha weighting by significance of difference in mean
for i = 1:12
    if pmld(i)<0.05
        h(i).FaceAlpha=0.8;
        h2(i).FaceAlpha=0.8;
    else
        h(i).FaceAlpha=0.4;
        h2(i).FaceAlpha=0.4;
    end
end



set(gca,'fontsize',fs)
set(gca,'ticklength',[0.01 0.01])
set(gca,'xtick',1:12)
set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'xlim',[0.5 12.5])
set(gca,'yminortick','on')
% ylabel('iPOC^{sub}\it_s\rm or iPOC^{mes}\it_s\rm / iPOC\it_s\rm');
ylabel('MLD [m]');

% AX(1).YColor=colors2(2,:);
grid on







%% intercompare different ezds
clearvars -except lb thrs fs hf1 print_flag scale_flag
ezd_p=[];ezd_c=[];ezd_iso=[];
mld=[];
mns=[]; %month for color?

for i = 1:length(lb)
     inds=sum(isnan(lb{i}.zbin.poc_s))<(thrs)*length(lb{i}.zbin.z) & lb{i}.zprod<900;
   
    ezd_p=[ezd_p lb{i}.ezd_p(inds)];
    ezd_c=[ezd_c lb{i}.ezd_c(inds)];
    ezd_iso=[ezd_iso lb{i}.ezd_iso(inds)];
    mld=[mld lb{i}.mld(inds)];
    mns=[mns month(lb{i}.dnum(inds))];
end


inds1=~isnan(ezd_c)&~isnan(ezd_p);
inds2=~isnan(ezd_c)&~isnan(ezd_iso);
inds3=~isnan(ezd_p)&~isnan(ezd_iso);
%% plot
% month colors
c=crameri('roma',12);
colorsm=[c(9:12,:); c(1:8,:)];

hf2=figure();
set(hf2,'Units','inches','Position', [5 5 5 5], 'PaperPosition', [0 0 5 5], 'PaperSize', [5 5]);
ha1=iSubplot(1,1, 'Gap', [0 0], 'Min', [0.13 0.13], 'Max', [0.96  0.96], 'XTickL', 'All', 'YTickL', 'All');

axes(ha1(1))
box on
grid on
axis square
xlim([0 200])
ylim([0 200])
plot(0:200,0:200,'k:')
s1=scatter(ezd_c,ezd_p,40,[0.6 0.3 0.1],'filled','s','markeredgecolor',[0.33 0.33 0.33],'linewidth',0.25,'MarkerFaceAlpha',0.85);
s2=scatter(ezd_c,ezd_iso,20,[0.6 0.7 0.7],'filled','^','markeredgecolor',[0.65 0.65 0.65],'linewidth',0.25,'MarkerFaceAlpha',0.6);

p1=polyfit(ezd_c(inds1),ezd_p(inds1),1);
r1=corr(ezd_c(inds1)',ezd_p(inds1)');
mr1=median(ezd_p(inds1)./ezd_c(inds1));
plot(0:200,p1(1)*(0:200)+p1(2),'color',[0.6 0.3 0.1])

p2=polyfit(ezd_c(inds2),ezd_iso(inds2),1);
r2=corr(ezd_c(inds2)',ezd_iso(inds2)');
mr2=median(ezd_iso(inds2)./ezd_c(inds2));
plot(0:200,p2(1)*(0:200)+p2(2),'color',[0.55 0.66 0.66])

mr3=median(ezd_p(inds3)./ezd_iso(inds3));

xlabel('EZD_{10%-Chla_{max}} [m]')
ylabel('EZD_{1%-PAR_0} or EZD_{iso} [m]')


lgd_text={'EZD_{1%-PAR_0}','EZD_{iso}'};
[hl,lines]=legendflex([s1 s2],lgd_text, 'ref', gca, ...
'anchor', {'nw', 'nw'}, 'buffer', [10 -10], 'xscale', 0.8, 'ncol',6,'nrow',2, ...
'box', 'off', 'FontSize', fs-3);

text(0.4,0.93,sprintf('y = %1.2fx + %1.1f m; MdR = %1.2f',p1(1),p1(2),mr1),'units','normalized','fontsize',fs-4,'color',[0.6 0.3 0.1],'fontweight','bold')
text(0.4,0.84,sprintf('y = %1.2fx + %1.1f m; MdR = %1.2f',p2(1),p2(2),mr2),'units','normalized','fontsize',fs-4,'color',[0.55 0.66 0.66],'fontweight','bold')
text(0.08,0.75,sprintf('MdR = %1.2f',mr3),'units','normalized','FontSize',fs-5,'FontWeight','bold')

set(gca,'fontsize',fs-3)

%% print
if print_flag==1

    figure(hf1);
    print(['Figures/V8/SI/LB22_EZD_MLD_Months_' scale_flag],'-dpdf','-r800')

    figure(hf2);
    print(['Figures/V8/SI/LB22_EZD_Compare_' scale_flag],'-dpdf','-r800')
end