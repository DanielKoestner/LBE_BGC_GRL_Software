%% iPOC violin plots to see spread for each month.
% work on total iPOC with top panel for inside eddy and bottom with outside
% use 3 sets of violin plots next to each other, PZ, TZ, and TZ+200
%% set up


close all
clear all

curdir=cd;

addpath(genpath([curdir '/Scripts']))
addpath(genpath([curdir '/Data']))

addpath(genpath([curdir '/aux']))

print_flag=0;

% size flag, 0 for total, 1 for small, 2 for large
size_flag=0;
scale_flag='A'; % A for 1.5r and B for 2r scaling of LBEZ


fs=17;
lw=1.5;
set(0, 'DefaultAxesFontName', 'Times');
set(0, 'DefaultTextFontName', 'Times');
alp=0.65;
mksz=7;

f=12; % LBE
colors=crameri('acton',f+1);
% colors=colors(2:end-2,:);
gcol=[0.2 0.45 0.85];

f=23; %LB
colors2=crameri('bamako',f+1);
colors2=colors2(2:end,:);

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
%% load data

% load('LBE_BGC_POC_2010_2022_MonthlyParameters_04-Aug-2024_200.mat')
% load('LBE_BGC_POC_2010_2022_MonthlyParameters_B_31-Dec-2024_200.mat'); % July 30 used 200 m, and trap integration
load(['LBE_BGC_POC_2010_2022_MonthlyParameters_' scale_flag '_06-Feb-2025_200.mat']); % July 30 used 200 m, and trap integration


colors=crameri('davos',7);

coltz=colors(4,:);
coltz2=colors(2,:);
colpz=[0.4 0.75 0.4];


if size_flag==0
    %set up total ipoc data
    for i = 1:12
        ipoc_pz_in{i}=(ipoc_eups_all{i}+ipoc_eups_l_all{i})./1000; % in g/m2
        ipoc_tz_in{i}=(ipoc_subs_all{i}+ipoc_subs_l_all{i})./1000;
        ipoc_tz2_in{i}=(ipoc_subbs_all{i}+ipoc_subbs_l_all{i})./1000;

        ipoc_pz_out{i}=(ipoc_eup2s_all{i}+ipoc_eup2s_l_all{i})./1000;
        ipoc_tz_out{i}=(ipoc_sub2s_all{i}+ipoc_sub2s_l_all{i})./1000;
        ipoc_tz2_out{i}=(ipoc_subb2s_all{i}+ipoc_subb2s_l_all{i})./1000;
    end
    lims=[0 20];

elseif size_flag==1

    %set up small ipoc data
    for i = 1:12
        ipoc_pz_in{i}=(ipoc_eups_all{i})./1000; % in g/m2
        ipoc_tz_in{i}=(ipoc_subs_all{i})./1000;
        ipoc_tz2_in{i}=(ipoc_subbs_all{i})./1000;

        ipoc_pz_out{i}=(ipoc_eup2s_all{i})./1000;
        ipoc_tz_out{i}=(ipoc_sub2s_all{i})./1000;
        ipoc_tz2_out{i}=(ipoc_subb2s_all{i})./1000;
    end
    lims=[0 15];

elseif size_flag==2
    %set up large ipoc data
    for i = 1:12
        ipoc_pz_in{i}=(ipoc_eups_l_all{i})./1000; % in g/m2
        ipoc_tz_in{i}=(ipoc_subs_l_all{i})./1000;
        ipoc_tz2_in{i}=(ipoc_subbs_l_all{i})./1000;

        ipoc_pz_out{i}=(ipoc_eup2s_l_all{i})./1000;
        ipoc_tz_out{i}=(ipoc_sub2s_l_all{i})./1000;
        ipoc_tz2_out{i}=(ipoc_subb2s_l_all{i})./1000;

    end
    lims=[0 8];
end


% remove nans and test significance of differences

for i = 1:12
ipoc_pz_in{i}=ipoc_pz_in{i}(~isnan(ipoc_pz_in{i}));
ipoc_tz_in{i}=ipoc_tz_in{i}(~isnan(ipoc_tz_in{i}));
ipoc_tz2_in{i}=ipoc_tz2_in{i}(~isnan(ipoc_tz2_in{i}));

ipoc_pz_out{i}=ipoc_pz_out{i}(~isnan(ipoc_pz_out{i}));
ipoc_tz_out{i}=ipoc_tz_out{i}(~isnan(ipoc_tz_out{i}));
ipoc_tz2_out{i}=ipoc_tz2_out{i}(~isnan(ipoc_tz2_out{i}));

[h_pz(i),p_pz(i),ci_pz(i,:),stats_pz{i}] = ttest2(ipoc_pz_in{i},ipoc_pz_out{i},'Vartype','unequal','alpha',0.05);
[h_tz(i),p_tz(i),ci_tz(i,:),stats_tz{i}] = ttest2(ipoc_tz_in{i},ipoc_tz_out{i},'Vartype','unequal','alpha',0.05);
[h_tz2(i),p_tz2(i),ci_tz2(i,:),stats_tz2{i}] = ttest2(ipoc_tz2_in{i},ipoc_tz2_out{i},'Vartype','unequal','alpha',0.05);

p2_pz(i) = ranksum(ipoc_pz_in{i},ipoc_pz_out{i});
p2_tz(i)= ranksum(ipoc_tz_in{i},ipoc_tz_out{i});
p2_tz2(i)= ranksum(ipoc_tz2_in{i},ipoc_tz2_out{i});


end

%%
hf1=figure();
set(hf1,'Units','inches','Position', [5 5 14 9], 'PaperPosition', [0 0 14 9], 'PaperSize', [14 9]);
ha=iSubplot(3,1, 'Gap', [0 0.015], 'Min', [0.03 0.04], 'Max', [.98 .97], 'XTickL', 'All', 'YTickL', 'All');


axes(ha(1));
text(-0.05,1,'(a)','units','normalized','fontsize',fs-3)
grid on
plot([3.5 3.5],[0 20],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([6.5 6.5],[0 20],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([9.5 9.5],[0 20],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
% plot([0.5 12.5],[0.5 0.5],'k:','linewidth',1.25)
hold on

text(1,1.04,'inside LBEZ','units','normalized','fontsize',fs-4,'fontweight','bold','HorizontalAlignment','right')

[h3,L3,MX3,MED3,bw3]=violin(ipoc_tz2_in,'x',[1.25:1:12.25],'facecolor',coltz2,'edgecolor','none','mc','k','medc','k-.');
[h,L,MX,MED,bw]=violin(ipoc_pz_in,'x',[0.75:1:11.75],'facecolor',colpz,'edgecolor','none','mc','k','medc','k-.');
[h2,L2,MX2,MED2,bw2]=violin(ipoc_tz_in,'x',[1:1:12],'facecolor',coltz,'edgecolor','none','mc','k','medc','k-.');

L2.Visible='off';
text(0.01,0.95,'Productive Zone','units','normalized','fontsize',fs-4,'FontWeight','bold','Color',colpz);
text(0.01,0.85,'Twilight Zone','units','normalized','fontsize',fs-4,'FontWeight','bold','Color',coltz);
text(0.01,0.75,'Lower Twilight Zone','units','normalized','fontsize',fs-4,'FontWeight','bold','Color',coltz2);

w=1-p_pz;
w2=1-p_tz;
w3=1-p_tz2;

%alpha weighting by significance of difference in mean
for i = 1:12
    
    h(i).FaceAlpha=w(i)*0.6+0.2;
    h2(i).FaceAlpha=w2(i)*0.6+0.2;
    h3(i).FaceAlpha=w3(i)*0.6+0.2;
end

ylim(lims)

set(gca,'fontsize',fs)
set(gca,'ticklength',[0.01 0.01])
set(gca,'xtick',1:12)
set(gca,'xticklabel',{''})
set(gca,'xlim',[0.5 12.5])
set(gca,'yminortick','on')
% ylabel('iPOC^{sub}\it_s\rm or iPOC^{mes}\it_s\rm / iPOC\it_s\rm');
ylabel('iPOC [g m^{-3}]');

% AX(1).YColor=colors2(2,:);
grid on

axes(ha(2));
text(-0.05,1,'(b)','units','normalized','fontsize',fs-3)
grid on
plot([3.5 3.5],[0 20],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([6.5 6.5],[0 20],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([9.5 9.5],[0 20],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
hold on


text(1,1.04,'outside LBEZ','units','normalized','fontsize',fs-4,'fontweight','bold','HorizontalAlignment','right')

[h3,L3,MX3,MED3,bw3]=violin(ipoc_tz2_out,'x',[1.25:1:12.25],'facecolor',coltz2,'edgecolor','none','mc','k','medc','k-.');
[h,L,MX,MED,bw]=violin(ipoc_pz_out,'x',[0.75:1:11.75],'facecolor',colpz,'edgecolor','none','mc','k','medc','k-.');
[h2,L2,MX2,MED2,bw2]=violin(ipoc_tz_out,'x',[1:1:12],'facecolor',coltz,'edgecolor','none','mc','k','medc','k-.');

L2.Location='northeast';
% text(0.01,0.25,'inside eddy','units','normalized','fontsize',fs-4,'FontWeight','bold','Color',colors(10,:));
% text(0.01,0.15,'outside eddy','units','normalized','fontsize',fs-4,'FontWeight','bold','Color',colors2(3,:)*1.1);

%alpha weighting by significance of difference in mean
for i = 1:12
    
    h(i).FaceAlpha=w(i)*0.6+0.2;
    h2(i).FaceAlpha=w2(i)*0.6+0.2;
    h3(i).FaceAlpha=w3(i)*0.6+0.2;
end

ylim(lims)

set(gca,'fontsize',fs)
set(gca,'ticklength',[0.01 0.01])
set(gca,'xtick',1:12)
set(gca,'yminortick','on')
set(gca,'xlim',[0.5 12.5])

% ylabel('iPOC^{sub}\it_s\rm or iPOC^{mes}\it_s\rm / iPOC\it_s\rm');
ylabel('iPOC [g m^{-3}]');

% AX(1).YColor=colors2(2,:);
grid on



set(gca,'xticklabel',{''})



% Plot p-values
axes(ha(3))
text(-0.05,1,'(c)','units','normalized','fontsize',fs-3)
b=bar([1-p_pz;1-p_tz;1-p_tz2]');
b(1).FaceColor=colpz;
b(2).FaceColor=coltz;
b(3).FaceColor=coltz2;
ylim([0 1])
ylabel('1 – \itp\rm-value')
set(gca,'ytick',[0:0.2:1])
plot([3.5 3.5],[0 20],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([6.5 6.5],[0 20],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([9.5 9.5],[0 20],'--','color',[0.4 0.4 0.4],'linewidth',0.9)


pgon=polyshape([0 0 12.5 12.5],[0.95 0 0 0.95]);
pg=plot(pgon);
pg.FaceColor=[0.6 0.6 0.6];
pg.EdgeColor='none';

% pgon=polyshape([0 0 12.5 12.5],[0.5 0 0 0.5]);
% pg=plot(pgon);
% pg.FaceColor=[0.4 0.4 0.4];
% pg.EdgeColor='none';

set(gca,'YMinorTick','on')
grid on
set(gca,'fontsize',fs)
set(gca,'ticklength',[0.01 0.01])
set(gca,'xtick',1:12)

set(gca,'xlim',[0.5 12.5])
set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
text(1,1.04,'two sample \itt\rm\bf-test ','units','normalized','fontsize',fs-4,'fontweight','bold','HorizontalAlignment','right')


%% Small to large POC ratios for each month
fs=13;

for i = 1:12
    in_t{i}=ipoc_tots_all{i}./ipoc_tots_l_all{i};
    out_t{i}=ipoc_tot2s_all{i}./ipoc_tot2s_l_all{i};

    in_pz{i}=ipoc_eups_all{i}./ipoc_eups_l_all{i};
    out_pz{i}=ipoc_eup2s_all{i}./ipoc_eup2s_l_all{i};

    in_tz{i}=ipoc_subs_all{i}./ipoc_subs_l_all{i};
    out_tz{i}=ipoc_sub2s_all{i}./ipoc_sub2s_l_all{i};
end

% tmp=out_pz{1};
% tmp(tmp==Inf)=NaN;
% out_pz{1}=tmp;

for i = 1:12
    tmp=in_pz{i};
    tmp(tmp>100)=NaN;
    in_pz{i}=tmp;

    tmp=out_pz{i};
    tmp(tmp>100)=NaN;
    out_pz{i}=tmp;

end

hf2=figure();
set(hf2,'Units','inches','Position', [5 5 8 3.5], 'PaperPosition', [0 0 8 3.5], 'PaperSize', [8 3.5]);
ha2=iSubplot(1,3, 'Gap', [0.01 0], 'Min', [0.035 0], 'Max', [0.98 1], 'XTickL', 'All', 'YTickL', 'All');

axes(ha2(1))
set(gca,'ticklength',[0.02 0.02])
box on
grid on
axis square
% set(gca,'yminortick','on')
set(gca,'fontsize',fs)
title('Productive Zone')
xlim([0 20])
ylim([0 20])
set(gca,'ytick',0:4:20)
set(gca,'xtick',0:4:20)
lims=xlim;
text(0.04,0.96,'(a)','units','normalized','fontsize',fs-3)

plot(lims(1):lims(end),lims(1):lims(end),'k:')
xlabel('\it s \rm: \it l \rm outside LBEZ')
ylabel('\it s \rm: \it l \rm inside LBEZ')

for i = 1:12
    [~,~,in_ci,~]=ttest(in_pz{i});
    [~,~,out_ci,~]=ttest(out_pz{i});

    errorbar(nanmean(out_pz{i}),nanmean(in_pz{i}),diff(in_ci)/2,diff(in_ci)/2,diff(out_ci)/2,diff(out_ci)/2,'color',colorsm(i,:),'CapSize',0)
    clearvars in_ci out_ci
end

for i = 1:12
    s(i)=scatter(nanmean(out_pz{i}),nanmean(in_pz{i}),((ipoc_eups(i,1)+ipoc_eups_l(i,1))./10000)*150,colorsm(i,:),'filled','markeredgecolor',[0.5 0.5 0.5],'markerfacealpha',0.8); 
end



lgd_text={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
[hl,lines]=legendflex([s],lgd_text, 'ref', gca, ...
'anchor', {'nw', 'nw'}, 'buffer', [0 4], 'xscale', 0.5, 'ncol',12,'nrow',1, ...
'box', 'off', 'FontSize', fs);

axes(ha2(2))
text(0.04,0.96,'(b)','units','normalized','fontsize',fs-3)

set(gca,'ticklength',[0.02 0.02])
box on
grid on
axis square
% set(gca,'yminortick','on')
set(gca,'fontsize',fs)
title('Twilight Zone')
xlim([0 5])
ylim([0 5])
set(gca,'ytick',0:1:5)
set(gca,'xtick',0:1:5)
lims=xlim;

plot(lims(1):lims(end),lims(1):lims(end),'k:')
xlabel('\it s \rm: \it l \rm outside LBEZ')

for i = 1:12
    [~,~,in_ci,~]=ttest(in_tz{i});
    [~,~,out_ci,~]=ttest(out_tz{i});

    errorbar(nanmean(out_tz{i}),nanmean(in_tz{i}),diff(in_ci)/2,diff(in_ci)/2,diff(out_ci)/2,diff(out_ci)/2,'color',colorsm(i,:),'CapSize',0)
    clearvars in_ci out_ci
end

for i = 1:12
    s(i)=scatter(nanmean(out_tz{i}),nanmean(in_tz{i}),((ipoc_subs(i,1)+ipoc_subs_l(i,1))./10000)*150,colorsm(i,:),'filled','markeredgecolor',[0.5 0.5 0.5],'markerfacealpha',0.8); 
end


s1(1)=scatter(10,10,(1000/10000)*150,'o','markeredgecolor',[0.5 0.5 0.5]);
s1(2)=scatter(100,100,(10000/10000)*150,'o','markeredgecolor',[0.5 0.5 0.5]);

[hl,lines]=legendflex([s1(1) s1(2)],[{'1 g m^{-2}','10 g m^{-2}'}], 'ref', gca, ...
'anchor', {'se', 'se'}, 'buffer', [0 50], 'xscale', 0.5, 'ncol',1,'nrow',2, ...
'box', 'off', 'FontSize', fs-1);

tmp=findobj(hl, 'type', 'patch');
tmp(2).MarkerSize=sqrt(15/pi)*2;
tmp(1).MarkerSize=sqrt(150/pi)*2;

text(0.745,0.3,'iPOC','units','normalized','fontsize',fs-1)
axes(ha2(3))
text(0.04,0.96,'(c)','units','normalized','fontsize',fs-3)
set(gca,'ticklength',[0.02 0.02])
box on
grid on
axis square
% set(gca,'yminortick','on')
set(gca,'fontsize',fs)
title('Total Water Column')
xlim([0 5])
ylim([0 5])
set(gca,'ytick',0:1:5)
set(gca,'xtick',0:1:5)
lims=xlim;

plot(lims(1):lims(end),lims(1):lims(end),'k:')
% ylabel('iPOC_s / iPOC_l outside eddy')
xlabel('\it s \rm: \it l \rm outside LBEZ')

for i = 1:12
    [~,~,in_ci,~]=ttest(in_t{i});
    [~,~,out_ci,~]=ttest(out_t{i});

    errorbar(nanmean(out_t{i}),nanmean(in_t{i}),diff(in_ci)/2,diff(in_ci)/2,diff(out_ci)/2,diff(out_ci)/2,'color',colorsm(i,:),'CapSize',0)
    clearvars in_ci out_ci
end

for i = 1:12
    s(i)=scatter(nanmean(out_t{i}),nanmean(in_t{i}),((ipoc_tots(i,1)+ipoc_tots_l(i,1))./10000)*150,colorsm(i,:),'filled','markeredgecolor',[0.5 0.5 0.5],'markerfacealpha',0.8); 
end



%% print

if print_flag==1

    figure(hf1)

    if size_flag==0
        print(['Figures/V8/SI/LB22_iPOC_Violins_Total_' scale_flag],'-dpdf','-r800')
    elseif size_flag==1
        print(['Figures/V8/SI/LB22_iPOC_Violins_Small_' scale_flag],'-dpdf','-r800')
    elseif size_flag==2
        print(['Figures/V8/SI/LB22_iPOC_Violins_Large_' scale_flag],'-dpdf','-r800')
    end

    figure(hf2)
    print(['Figures/V8/SI/LB22_iPOC_SizeRatios_INvOUT_' scale_flag],'-dpdf','-r800')
end

