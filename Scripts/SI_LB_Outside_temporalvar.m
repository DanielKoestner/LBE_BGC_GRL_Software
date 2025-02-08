% Look into temporal variability in iPOC

% Make a total and TZ version


close all
clear all

curdir=cd;

addpath(genpath([curdir '/Scripts']))
addpath(genpath([curdir '/Data']))


print_flag=1;
save_flag=0;
scale_flag='A'; % A for 1.5r and B for 2r scaling of LBEZ

fs=15;
lw=1.5;
set(0, 'DefaultAxesFontName', 'Times');
set(0, 'DefaultTextFontName', 'Times');

%% Get data

% load('LBE_BGC_POC_2010_2022_MonthlyParameters_04-Aug-2024_200.mat');
load(['LBE_BGC_POC_2010_2022_MonthlyParameters_' scale_flag '_06-Feb-2025_200.mat']); % July 30 used 200 m, and trap integration

% need to get years, already have months
for i = 1:length(days)
    years(i)=year(days(i));
end
for i = 1:length(days2)
    years2(i)=year(days2(i));
end

txt={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
%% outside
hf=figure();
set(hf,'Units','inches','Position', [5 5 8 9], 'PaperPosition', [0 0 8 9], 'PaperSize', [8 9]);
ha1=iSubplot(3,4, 'Gap', [0 0.08], 'Min', [0.045 0.04], 'Max', [0.99 0.94], 'XTickL', 'All', 'YTickL', 'All');



colors=crameri('lapaz',13);
nyrs=2010:2022;
for ii = 1:12
    inds=ismember(months2,ii);
    indsIN=ismember(months,ii); %for this month, index for inside lbez
    uni_years_in=unique(years(indsIN)); %unique years inside lbez for this month

    tmpyrs2=years2(inds);
    tmppoc2=(ipoc_sub2(inds)+ipoc_sub2_l(inds))/1000;

    r=corr(tmpyrs2(~isnan(tmppoc2))',tmppoc2(~isnan(tmppoc2)));

    axes(ha1(ii));
    title(txt(ii))
    text(0.7,1.05,sprintf('\\itr\\rm = %1.2f',r),'units','normalized')
    if ttest2(tmppoc2(ismember(tmpyrs2,uni_years_in)),tmppoc2)==1
        text(0.6,1.1,'*','units','normalized')
    end

    xlim([2009 2023])
    ylim([0 20])
    set(gca,'ticklength',[0.02 0.02])
    set(gca,'xminortick','off')
    set(gca,'yticklabel',{''})
    box on
    grid on
    set(gca,'fontsize',fs-2)
    set(gca,'xtick',2010:2:2022)

    hold on

    xx=[2009 2023 2023 2009];
    
%     % % % % USING "spread" in data
%     yy=[ipoc_sub2s(ii,4)+ipoc_sub2s_l(ii,4) ipoc_sub2s(ii,4)+ipoc_sub2s_l(ii,4) ipoc_sub2s(ii,5)+ipoc_sub2s_l(ii,5) ipoc_sub2s(ii,5)+ipoc_sub2s_l(ii,5) ];
%     yy2=[(ipoc_sub2s(ii,1)-ipoc_sub2s(ii,2))+(ipoc_sub2s_l(ii,1)-ipoc_sub2s_l(ii,2)) (ipoc_sub2s(ii,1)-ipoc_sub2s(ii,2))+(ipoc_sub2s_l(ii,1)-ipoc_sub2s_l(ii,2)) (ipoc_sub2s(ii,1)+ipoc_sub2s(ii,2))+(ipoc_sub2s_l(ii,1)+ipoc_sub2s_l(ii,2)) (ipoc_sub2s(ii,1)+ipoc_sub2s(ii,2))+(ipoc_sub2s_l(ii,1)+ipoc_sub2s_l(ii,2))];
% 
%     s2=fill(xx,yy2/1000,[0.3 0.3 0.3]);
%     s2.FaceAlpha=0.25;
%     s2.EdgeColor='none';
% 
%     s1=fill(xx,yy/1000,[0.6 0.6 0.6]);
%     s1.FaceAlpha=0.25;
%     s1.EdgeColor=[0.6 0.6 0.6];
%     s1.LineStyle=':';
%     s1.LineWidth=1;


    % % % using confidence intervals for both populations
    [~,~,ci,~]=ttest(tmppoc2(ismember(tmpyrs2,uni_years_in)));
    yy=[ci(1) ci(1) ci(2) ci(2)];
    yy2=[ipoc_sub2s(ii,6)+ipoc_sub2s_l(ii,6) ipoc_sub2s(ii,6)+ipoc_sub2s_l(ii,6) ipoc_sub2s(ii,7)+ipoc_sub2s_l(ii,7) ipoc_sub2s(ii,7)+ipoc_sub2s_l(ii,7)];
     s2=fill(xx,yy2/1000,[0.3 0.3 0.3]);
    s2.FaceAlpha=0.25;
    s2.EdgeColor='none';

    s1=fill(xx,yy,[0.6 0.6 0.6]);
    s1.FaceAlpha=0.25;
    s1.EdgeColor=[0.6 0.6 0.6];
    s1.LineStyle=':';
    s1.LineWidth=1;


    plot([2009 2023],[ipoc_sub2s(ii,1)+ipoc_sub2s_l(ii,1) ipoc_sub2s(ii,1)+ipoc_sub2s_l(ii,1)]/1000,'-','color',[0.1 0.1 0.1],'linewidth',1.5)
    plot([2009 2023],[nanmean(tmppoc2(ismember(tmpyrs2,uni_years_in))) nanmean(tmppoc2(ismember(tmpyrs2,uni_years_in)))],'--','color',[0.3 0.3 0.3],'linewidth',1.3)
    for i = 1:length(nyrs)
        inds=ismember(tmpyrs2,nyrs(i));
        if sum(ismember(uni_years_in,nyrs(i)))==0
            errorbar(nyrs(i),nanmean(tmppoc2(inds)),nanstd(tmppoc2(inds)),'marker','s','MarkerFaceColor',colors(i,:),'color',colors(i,:),'MarkerEdgeColor','none','linewidth',1)
        elseif sum(ismember(uni_years_in,nyrs(i)))==1
            errorbar(nyrs(i),nanmean(tmppoc2(inds)),nanstd(tmppoc2(inds)),'marker','s','MarkerFaceColor',colors(i,:),'color',colors(i,:),'MarkerEdgeColor',[0.3 0.3 0.3],'linewidth',1.3)
        end
    end

    clearvars r tmp*
end

for i = 1:4:12
axes(ha1(i))
ylabel('iPOC^{TZ} [g m^{-2}]')
set(gca,'yticklabel',num2cell(0:5:25))
end

axes(ha1(1))
text(0,1.2,'outside LBEZ','fontsize',fs-3,'units','normalized','fontweight','bold')

%% inside

hf2=figure();
set(hf2,'Units','inches','Position', [5 5 8 9], 'PaperPosition', [0 0 8 9], 'PaperSize', [8 9]);
ha1=iSubplot(3,4, 'Gap', [0 0.08], 'Min', [0.045 0.04], 'Max', [0.99 0.94], 'XTickL', 'All', 'YTickL', 'All');


for ii = 1:12
    inds=ismember(months,ii);

    tmpyrs=years(inds);
    tmppoc=(ipoc_sub(inds)+ipoc_sub_l(inds))/1000;

    r=corr(tmpyrs(~isnan(tmppoc))',tmppoc(~isnan(tmppoc)));

    axes(ha1(ii));
    title(txt(ii))
    text(0.7,1.05,sprintf('\\itr\\rm = %1.2f',r),'units','normalized')
    xlim([2009 2023])
    ylim([0 20])
    set(gca,'ticklength',[0.02 0.02])
    set(gca,'xminortick','off')
    set(gca,'yticklabel',{''})
    box on
    grid on
    set(gca,'fontsize',fs-2)
    set(gca,'xtick',2010:2:2022)

    hold on
    
    xx=[2009 2023 2023 2009];

    % % % % USING "spread" in data
%     yy=[ipoc_subs(ii,4)+ipoc_subs_l(ii,4) ipoc_subs(ii,4)+ipoc_subs_l(ii,4) ipoc_subs(ii,5)+ipoc_subs_l(ii,5) ipoc_subs(ii,5)+ipoc_subs_l(ii,5) ];
%     yy2=[(ipoc_subs(ii,1)-ipoc_subs(ii,2))+(ipoc_subs_l(ii,1)-ipoc_subs_l(ii,2)) (ipoc_subs(ii,1)-ipoc_subs(ii,2))+(ipoc_subs_l(ii,1)-ipoc_subs_l(ii,2)) (ipoc_subs(ii,1)+ipoc_subs(ii,2))+(ipoc_subs_l(ii,1)+ipoc_subs_l(ii,2)) (ipoc_subs(ii,1)+ipoc_subs(ii,2))+(ipoc_subs_l(ii,1)+ipoc_subs_l(ii,2))];
% 
%     s2=fill(xx,yy2/1000,[0.3 0.3 0.3]);
%     s2.FaceAlpha=0.25;
%     s2.EdgeColor='none';
% 
%     s1=fill(xx,yy/1000,[0.6 0.6 0.6]);
%     s1.FaceAlpha=0.25;
%     s1.EdgeColor=[0.6 0.6 0.6];
%     s1.LineStyle=':';
%     s1.LineWidth=1;

        % % % using confidence intervals for both populations
    yy2=[ipoc_subs(ii,6)+ipoc_subs_l(ii,6) ipoc_subs(ii,6)+ipoc_subs_l(ii,6) ipoc_subs(ii,7)+ipoc_subs_l(ii,7) ipoc_subs(ii,7)+ipoc_subs_l(ii,7)];
     s2=fill(xx,yy2/1000,[0.3 0.3 0.3]);
    s2.FaceAlpha=0.25;
    s2.EdgeColor='none';


    plot([2009 2023],[ipoc_subs(ii,1)+ipoc_subs_l(ii,1) ipoc_subs(ii,1)+ipoc_subs_l(ii,1)]/1000,'--','color',[0.3 0.3 0.3],'linewidth',1.3)
   
    for i = 1:length(nyrs)
        inds=ismember(tmpyrs,nyrs(i));
        errorbar(nyrs(i),nanmean(tmppoc(inds)),nanstd(tmppoc(inds)),'marker','o','MarkerFaceColor',colors(i,:),'color',colors(i,:),'MarkerEdgeColor',[0.3 0.3 0.3],'linewidth',1.3)
    end

    clearvars r tmp*
end

for i = 1:4:12
axes(ha1(i))
ylabel('iPOC^{TZ} [g m^{-2}]')
set(gca,'yticklabel',num2cell(0:5:25))
end
axes(ha1(1))
text(0,1.2,'inside LBEZ','fontsize',fs-3,'units','normalized','fontweight','bold')


%% print

if print_flag==1
    figure(hf)
    print(['Figures/V8/SI/LB22_iPOCtz_years_outside_' scale_flag],'-dpdf','-r800')

    figure(hf2)
    print(['Figures/V8/SI/LB22_iPOCtz_years_inside_' scale_flag],'-dpdf','-r800')
end