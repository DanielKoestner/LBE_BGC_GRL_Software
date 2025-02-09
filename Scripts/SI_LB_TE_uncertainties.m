% monthly average TEs w/ error bars
% TE +/- 1 sd results using
% surf plot feature. Do it as 2x2, inside left, outside right, 200 m top,
% 500 m bottom.

% repeat with std and quantiles (15 and 85)
%% set up
close all
clear all

curdir=cd;

addpath(genpath([curdir '/Scripts']))
addpath(genpath([curdir '/Data']))

addpath(genpath([curdir '/aux']))

print_flag=0;
save_flag=0;
err_flag=3; %1 for std, 2 for quantiles (15 and 85 percentile), 3 for 95% confidence interval
scale_flag='A'; % A for 1.5r and B for 2r scaling of LBEZ

fs=14;
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
%% Try with profile by profile calculated iPOCsub

load(['LBE_BGC_POC_2010_2022_MonthlyParameters_' scale_flag '_06-Feb-2025_200.mat']); % July 30 used 200 m, and trap integration

%1 for mean, 2 for std, 3 for median, 4 for iqr
a=1;
b=2;

% % % % % IN
iPOC_ez_in=ipoc_eups(:,a)+ipoc_eups_l(:,a);
iPOC_tz_in=ipoc_subs(:,a)+ipoc_subs_l(:,a);
iPOC_tz2_in=ipoc_subbs(:,a)+ipoc_subbs_l(:,a);

iPOCs_ez_in=ipoc_eups(:,a);
iPOCs_tz_in=ipoc_subs(:,a);
iPOCs_tz2_in=ipoc_subbs(:,a);

iPOCl_ez_in=ipoc_eups_l(:,a);
iPOCl_tz_in=ipoc_subs_l(:,a);
iPOCl_tz2_in=ipoc_subbs_l(:,a);

iPOC_ez_in_md=ipoc_eups(:,3)+ipoc_eups_l(:,3);
iPOC_tz_in_md=ipoc_subs(:,3)+ipoc_subs_l(:,3);
iPOC_tz2_in_md=ipoc_subbs(:,3)+ipoc_subbs_l(:,3);

if err_flag==1 % use std deviation
    iPOC_ez_in_err(1,:)=(ipoc_eups(:,a)-ipoc_eups(:,b))+(ipoc_eups_l(:,a)-ipoc_eups_l(:,b));
    iPOC_ez_in_err(2,:)=(ipoc_eups(:,a)+ipoc_eups(:,b))+(ipoc_eups_l(:,a)+ipoc_eups_l(:,b));
    iPOC_tz_in_err(1,:)=(ipoc_subs(:,a)-ipoc_subs(:,b))+(ipoc_subs_l(:,a)-ipoc_subs_l(:,b));
    iPOC_tz_in_err(2,:)=(ipoc_subs(:,a)+ipoc_subs(:,b))+(ipoc_subs_l(:,a)+ipoc_subs_l(:,b));
    iPOC_tz2_in_err(1,:)=(ipoc_subbs(:,a)-ipoc_subbs(:,b))+(ipoc_subbs_l(:,a)-ipoc_subbs_l(:,b));
    iPOC_tz2_in_err(2,:)=(ipoc_subbs(:,a)+ipoc_subbs(:,b))+(ipoc_subbs_l(:,a)+ipoc_subbs_l(:,b));

elseif err_flag==2 % use 15th and 85th percentiles
    iPOC_tz_in_err(1,:)=ipoc_subs(:,4)+ipoc_subs_l(:,4);
    iPOC_tz_in_err(2,:)=ipoc_subs(:,5)+ipoc_subs_l(:,5);
    iPOC_ez_in_err(1,:)=ipoc_eups(:,4)+ipoc_eups_l(:,4); 
    iPOC_ez_in_err(2,:)=ipoc_eups(:,5)+ipoc_eups_l(:,5); 
    iPOC_tz2_in_err(1,:)=ipoc_subbs(:,4)+ipoc_subbs_l(:,4);
    iPOC_tz2_in_err(2,:)=ipoc_subbs(:,5)+ipoc_subbs_l(:,5);

elseif err_flag==3 % use CI
    iPOC_tz_in_err(1,:)=ipoc_subs(:,6)+ipoc_subs_l(:,6);
    iPOC_tz_in_err(2,:)=ipoc_subs(:,7)+ipoc_subs_l(:,7);
    iPOC_ez_in_err(1,:)=ipoc_eups(:,6)+ipoc_eups_l(:,6);
    iPOC_ez_in_err(2,:)=ipoc_eups(:,7)+ipoc_eups_l(:,7);
    iPOC_tz2_in_err(1,:)=ipoc_subbs(:,6)+ipoc_subbs_l(:,6);
    iPOC_tz2_in_err(2,:)=ipoc_subbs(:,7)+ipoc_subbs_l(:,7);

end


% % % % % OUT
iPOC_ez_out=ipoc_eup2s(:,a)+ipoc_eup2s_l(:,a);
iPOC_tz_out=ipoc_sub2s(:,a)+ipoc_sub2s_l(:,a);
iPOC_tz2_out=ipoc_subb2s(:,a)+ipoc_subb2s_l(:,a);

iPOCs_ez_out=ipoc_eup2s(:,a);
iPOCs_tz_out=ipoc_sub2s(:,a);
iPOCs_tz2_out=ipoc_subb2s(:,a);

iPOCl_ez_out=ipoc_eup2s_l(:,a);
iPOCl_tz_out=ipoc_sub2s_l(:,a);
iPOCl_tz2_out=ipoc_subb2s_l(:,a);

iPOC_ez_out_md=ipoc_eup2s(:,3)+ipoc_eup2s_l(:,3);
iPOC_tz_out_md=ipoc_sub2s(:,3)+ipoc_sub2s_l(:,3);
iPOC_tz2_out_md=ipoc_subb2s(:,3)+ipoc_subb2s_l(:,3);

if err_flag==1 % use std deviation
    iPOC_ez_out_err(1,:)=(ipoc_eup2s(:,a)-ipoc_eup2s(:,b))+(ipoc_eup2s_l(:,a)-ipoc_eup2s_l(:,b));
    iPOC_ez_out_err(2,:)=(ipoc_eup2s(:,a)+ipoc_eup2s(:,b))+(ipoc_eup2s_l(:,a)+ipoc_eup2s_l(:,b));
    iPOC_tz_out_err(1,:)=(ipoc_sub2s(:,a)-ipoc_sub2s(:,b))+(ipoc_sub2s_l(:,a)-ipoc_sub2s_l(:,b));
    iPOC_tz_out_err(2,:)=(ipoc_sub2s(:,a)+ipoc_sub2s(:,b))+(ipoc_sub2s_l(:,a)+ipoc_sub2s_l(:,b));
    iPOC_tz2_out_err(1,:)=(ipoc_subb2s(:,a)-ipoc_subb2s(:,b))+(ipoc_subb2s_l(:,a)-ipoc_subb2s_l(:,b));
    iPOC_tz2_out_err(2,:)=(ipoc_subb2s(:,a)+ipoc_subb2s(:,b))+(ipoc_subb2s_l(:,a)+ipoc_subb2s_l(:,b));

elseif err_flag==2 % use 15th and 85th percentiles
    iPOC_tz2_out_err(1,:)=ipoc_subb2s(:,4)+ipoc_subb2s_l(:,4);
    iPOC_tz2_out_err(2,:)=ipoc_subb2s(:,5)+ipoc_subb2s_l(:,5);
    iPOC_tz_out_err(1,:)=ipoc_sub2s(:,4)+ipoc_sub2s_l(:,4);
    iPOC_tz_out_err(2,:)=ipoc_sub2s(:,5)+ipoc_sub2s_l(:,5);
    iPOC_ez_out_err(1,:)=ipoc_eup2s(:,4)+ipoc_eup2s_l(:,4);
    iPOC_ez_out_err(2,:)=ipoc_eup2s(:,5)+ipoc_eup2s_l(:,5);

elseif err_flag==3 % use CI
    iPOC_tz2_out_err(1,:)=ipoc_subb2s(:,6)+ipoc_subb2s_l(:,6);
    iPOC_tz2_out_err(2,:)=ipoc_subb2s(:,7)+ipoc_subb2s_l(:,7);
    iPOC_tz_out_err(1,:)=ipoc_sub2s(:,6)+ipoc_sub2s_l(:,6);
    iPOC_tz_out_err(2,:)=ipoc_sub2s(:,7)+ipoc_sub2s_l(:,7);
    iPOC_ez_out_err(1,:)=ipoc_eup2s(:,6)+ipoc_eup2s_l(:,6);
    iPOC_ez_out_err(2,:)=ipoc_eup2s(:,7)+ipoc_eup2s_l(:,7);

end

%% Rates
% remove any negative iPOC error
iPOC_ez_in_err(iPOC_ez_in_err<0)=0;
iPOC_tz_in_err(iPOC_tz_in_err<0)=0;
iPOC_tz2_in_err(iPOC_tz2_in_err<0)=0;

iPOC_ez_out_err(iPOC_ez_out_err<0)=0;
iPOC_tz_out_err(iPOC_tz_out_err<0)=0;
iPOC_tz2_out_err(iPOC_tz2_out_err<0)=0;

% % % %  rates, P and E
% dt=[31,28,31,30,31,30,31,31,30,31,30,31]; %days per month
dt=ones(1,12)*(356/12);
% dt=[29.5,29.5,30.5,30.5,30.5,30.5,31,30.5,30.5,30.5,30.5,31]; %avg days for each month (e.g., jan=(jan+feb)/2)
for i = 1:12
    if i <12
        % in
        P_in(i)= (iPOC_ez_in(i+1)-iPOC_ez_in(i))/dt(i);
        Ez_in(i)= (iPOC_tz_in(i+1)-iPOC_tz_in(i))/dt(i);
        Ez2_in(i)= (iPOC_tz2_in(i+1)-iPOC_tz2_in(i))/dt(i);

        P_s_in(i)= (iPOCs_ez_in(i+1)-iPOCs_ez_in(i))/dt(i);
        Ez_s_in(i)= (iPOCs_tz_in(i+1)-iPOCs_tz_in(i))/dt(i);
        Ez2_s_in(i)= (iPOCs_tz2_in(i+1)-iPOCs_tz2_in(i))/dt(i);

        P_l_in(i)= (iPOCl_ez_in(i+1)-iPOCl_ez_in(i))/dt(i);
        Ez_l_in(i)= (iPOCl_tz_in(i+1)-iPOCl_tz_in(i))/dt(i);
        Ez2_l_in(i)= (iPOCl_tz2_in(i+1)-iPOCl_tz2_in(i))/dt(i);


        Ez_in_md(i)= (iPOC_tz_in_md(i+1)-iPOC_tz_in_md(i))/dt(i);
        Ez2_in_md(i)= (iPOC_tz2_in_md(i+1)-iPOC_tz2_in_md(i))/dt(i);

        P_in_err(:,i)= (iPOC_ez_in_err(:,i+1)-iPOC_ez_in_err(:,i))./dt(i);
        Ez_in_err(:,i)= (iPOC_tz_in_err(:,i+1)-iPOC_tz_in_err(:,i))./dt(i);
        Ez2_in_err(:,i)= (iPOC_tz2_in_err(:,i+1)-iPOC_tz2_in_err(:,i))./dt(i);

        % out
        P_out(i)= (iPOC_ez_out(i+1)-iPOC_ez_out(i))/dt(i);
        Ez_out(i)= (iPOC_tz_out(i+1)-iPOC_tz_out(i))/dt(i);
        Ez2_out(i)= (iPOC_tz2_out(i+1)-iPOC_tz2_out(i))/dt(i);

        P_s_out(i)= (iPOCs_ez_out(i+1)-iPOCs_ez_out(i))/dt(i);
        Ez_s_out(i)= (iPOCs_tz_out(i+1)-iPOCs_tz_out(i))/dt(i);
        Ez2_s_out(i)= (iPOCs_tz2_out(i+1)-iPOCs_tz2_out(i))/dt(i);

        P_l_out(i)= (iPOCl_ez_out(i+1)-iPOCl_ez_out(i))/dt(i);
        Ez_l_out(i)= (iPOCl_tz_out(i+1)-iPOCl_tz_out(i))/dt(i);
        Ez2_l_out(i)= (iPOCl_tz2_out(i+1)-iPOCl_tz2_out(i))/dt(i);

        Ez_out_md(i)= (iPOC_tz_out_md(i+1)-iPOC_tz_out_md(i))/dt(i);
        Ez2_out_md(i)= (iPOC_tz2_out_md(i+1)-iPOC_tz2_out_md(i))/dt(i);

        P_out_err(:,i)= (iPOC_ez_out_err(:,i+1)-iPOC_ez_out_err(:,i))./dt(i);
        Ez_out_err(:,i)= (iPOC_tz_out_err(:,i+1)-iPOC_tz_out_err(:,i))./dt(i);
        Ez2_out_err(:,i)= (iPOC_tz2_out_err(:,i+1)-iPOC_tz2_out_err(:,i))./dt(i);

    else
% in
        P_in(i)= (iPOC_ez_in(1)-iPOC_ez_in(i))/dt(i);
        Ez_in(i)= (iPOC_tz_in(1)-iPOC_tz_in(i))/dt(i);
        Ez2_in(i)= (iPOC_tz2_in(1)-iPOC_tz2_in(i))/dt(i);

        P_s_in(i)= (iPOCs_ez_in(1)-iPOCs_ez_in(i))/dt(i);
        Ez_s_in(i)= (iPOCs_tz_in(1)-iPOCs_tz_in(i))/dt(i);
        Ez2_s_in(i)= (iPOCs_tz2_in(1)-iPOCs_tz2_in(i))/dt(i);

        P_l_in(i)= (iPOCl_ez_in(1)-iPOCl_ez_in(i))/dt(i);
        Ez_l_in(i)= (iPOCl_tz_in(1)-iPOCl_tz_in(i))/dt(i);
        Ez2_l_in(i)= (iPOCl_tz2_in(1)-iPOCl_tz2_in(i))/dt(i);

        Ez_in_md(i)= (iPOC_tz_in_md(1)-iPOC_tz_in_md(i))/dt(i);
        Ez2_in_md(i)= (iPOC_tz2_in_md(1)-iPOC_tz2_in_md(i))/dt(i);

        P_in_err(:,i)= (iPOC_ez_in_err(:,1)-iPOC_ez_in_err(:,i))./dt(i);
        Ez_in_err(:,i)= (iPOC_tz_in_err(:,1)-iPOC_tz_in_err(:,i))./dt(i);
        Ez2_in_err(:,i)= (iPOC_tz2_in_err(:,1)-iPOC_tz2_in_err(:,i))./dt(i);

% out
        P_out(i)= (iPOC_ez_out(1)-iPOC_ez_out(i))/dt(i);
        Ez_out(i)= (iPOC_tz_out(1)-iPOC_tz_out(i))/dt(i);
        Ez2_out(i)= (iPOC_tz2_out(1)-iPOC_tz2_out(i))/dt(i);

        P_s_out(i)= (iPOCs_ez_out(1)-iPOCs_ez_out(i))/dt(i);
        Ez_s_out(i)= (iPOCs_tz_out(1)-iPOCs_tz_out(i))/dt(i);
        Ez2_s_out(i)= (iPOCs_tz2_out(1)-iPOCs_tz2_out(i))/dt(i);

        P_l_out(i)= (iPOCl_ez_out(1)-iPOCl_ez_out(i))/dt(i);
        Ez_l_out(i)= (iPOCl_tz_out(1)-iPOCl_tz_out(i))/dt(i);
        Ez2_l_out(i)= (iPOCl_tz2_out(1)-iPOCl_tz2_out(i))/dt(i);

        Ez_out_md(i)= (iPOC_tz_out_md(1)-iPOC_tz_out_md(i))/dt(i);
        Ez2_out_md(i)= (iPOC_tz2_out_md(1)-iPOC_tz2_out_md(i))/dt(i);

        P_out_err(:,i)= (iPOC_ez_out_err(:,1)-iPOC_ez_out_err(:,i))./dt(i);
        Ez_out_err(:,i)= (iPOC_tz_out_err(:,1)-iPOC_tz_out_err(:,i))./dt(i);
        Ez2_out_err(:,i)= (iPOC_tz2_out_err(:,1)-iPOC_tz2_out_err(:,i))./dt(i);

    end
end


TE_in=Ez2_in./Ez_in;
TE_s_in=Ez2_s_in./Ez_s_in;
TE_l_in=Ez2_l_in./Ez_l_in;
TE_in_md=Ez2_in_md./Ez_in_md;
TE_in_err=Ez2_in_err./Ez_in_err;

TE_out=Ez2_out./Ez_out;
TE_s_out=Ez2_s_out./Ez_s_out;
TE_l_out=Ez2_l_out./Ez_l_out;
TE_out_md=Ez2_out_md./Ez_out_md;
TE_out_err=Ez2_out_err./Ez_out_err;




%% Plot

hf=figure();
set(hf,'Units','inches','Position', [5 5 8 6], 'PaperPosition', [0 0 8 6], 'PaperSize', [8 6]);
ha1=iSubplot(2,2, 'Gap', [0 0.05], 'Min', [0.065 0.05], 'Max', [0.98 0.96], 'XTickL', 'All', 'YTickL', 'All');

axes(ha1(1))
colors2=crameri('batlow',8);
text(0.015,0.92,'(a)','units','normalized','fontsize',fs-4)
text(0.75,1.04,'inside LBEZ','units','normalized','fontsize',fs-3,'fontweight','bold')
hold on
% plot([0 13],[0 0],'k-')
x=1:12;
y=TE_in; % set to shallowest depth horizon (100 or 200)
xx = [x;x];
yy = [min(TE_in_err);max(TE_in_err)];
alps=[ones(1,12)*0.4;ones(1,12)*0.4];
s1=surf(xx,yy,xx*0,...
    'alphadata',alps,...
    'facealpha','interp',...
    'edgecolor','none','FaceColor',colors2(7,:));


% p1=plot(1:12,TE_in,'color',colors2(7,:),'linewidth',0.75,'marker','o','MarkerFaceColor',colors2(7,:),'MarkerSize',1,'MarkerEdgeColor','none');
ps=plot(1:12,TE_s_in,'--','color',colors2(7,:),'linewidth',1.5);
pl=plot(1:12,TE_l_in,':','color',colors2(7,:),'linewidth',1.5);

% scatter(1:12,TE_in,iPOC_tz_in/115,'o','filled','MarkerFaceColor',colors2(7,:),'MarkerFaceAlpha',0.9);
s=scatter(1:12,TE_in,abs(Ez_in),'o','filled','MarkerFaceColor',colors2(7,:),'MarkerFaceAlpha',0.9);
% scatter(1:12,TE_in_md,50,'*','filled','MarkerEdgeColor',colors2(7,:),'MarkerFaceAlpha',0.9);

set(gca,'ticklength',[0.01 0.005])
box on
% grid on
ha1(1).YGrid='on';
ha1(1).YMinorTick='on';
set(gca,'ytick',[0:0.25:2])
set(gca,'fontsize',fs-2)
xlim([0.4 12.6])
ylim([0 1.5])
ylabel('TE^{200}')
set(gca,'xtick',[0.5:11.5],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});


[hl,lines]=legendflex([s(1) ps(1) pl(1)],{'total','small','large'}, 'ref', gca, ...
'anchor', {'ne', 'ne'}, 'buffer', [-5 -2], 'xscale', 0.5, 'ncol',1,'nrow',3, ...
'box', 'off', 'FontSize', fs-3);
% lines(5).Color='k';
% lines(6).Color='k';
% lines(7).Color='k';


%outside

axes(ha1(2))
colors2=crameri('batlow',8);
text(0.015,0.92,'(b)','units','normalized','fontsize',fs-4)
text(0.75,1.04,'outside LBEZ','units','normalized','fontsize',fs-3,'fontweight','bold')
hold on
% plot([0 13],[0 0],'k-')
x=1:12;
y=TE_out; % set to shallowest depth horizon (100 or 200)
xx = [x;x];
yy = [min(TE_out_err);max(TE_out_err)];
alps=[ones(1,12)*0.4;ones(1,12)*0.4];
s1=surf(xx,yy,xx*0,...
    'alphadata',alps,...
    'facealpha','interp',...
    'edgecolor','none','FaceColor',colors2(3,:));

% p1=plot(1:12,TE_out,'color',colors2(3,:),'linewidth',0.75,'marker','o','MarkerFaceColor',colors2(3,:),'MarkerSize',1,'MarkerEdgeColor','none');
ps=plot(1:12,TE_s_out,'--','color',colors2(3,:),'linewidth',1.5);
pl=plot(1:12,TE_l_out,':','color',colors2(3,:),'linewidth',1.5);

% scatter(1:12,TE_out,iPOC_tz_out/115,'s','filled','MarkerFaceColor',colors2(3,:),'MarkerFaceAlpha',0.9);
s=scatter(1:12,TE_out,abs(Ez_out),'s','filled','MarkerFaceColor',colors2(3,:),'MarkerFaceAlpha',0.9);
% scatter(1:12,TE_out_md,50,'*','filled','MarkerEdgeColor',colors2(3,:),'MarkerFaceAlpha',0.9);

set(gca,'ticklength',[0.01 0.005])
box on
% grid on
ha1(2).YGrid='on';
ha1(2).YMinorTick='on';
set(gca,'ytick',[0:0.25:2])
set(gca,'yticklabel',{''})
set(gca,'fontsize',fs-2)
xlim([0.4 12.6])
ylim([0 1.5])
% ylabel('TE')
set(gca,'xtick',[0.5:11.5],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});

[hl,lines]=legendflex([ s(1) ps(1) pl(1)],{'total','small','large'}, 'ref', gca, ...
'anchor', {'nw', 'nw'}, 'buffer', [20 -2], 'xscale', 0.5, 'ncol',1,'nrow',3, ...
'box', 'off', 'FontSize', fs-3);
% lines(5).Color='k';
% lines(6).Color='k';
% lines(7).Color='k';
%% deep

load(['LBE_BGC_POC_2010_2022_MonthlyParameters_' scale_flag '_06-Feb-2025_500.mat']); % July 30 used 200 m, and trap integration

% % % % % IN
iPOC_ez_in=ipoc_eups(:,a)+ipoc_eups_l(:,a);
iPOC_tz_in=ipoc_subs(:,a)+ipoc_subs_l(:,a);
iPOC_tz2_in=ipoc_subbs(:,a)+ipoc_subbs_l(:,a);

iPOCs_ez_in=ipoc_eups(:,a);
iPOCs_tz_in=ipoc_subs(:,a);
iPOCs_tz2_in=ipoc_subbs(:,a);

iPOCl_ez_in=ipoc_eups_l(:,a);
iPOCl_tz_in=ipoc_subs_l(:,a);
iPOCl_tz2_in=ipoc_subbs_l(:,a);

iPOC_ez_in_md=ipoc_eups(:,3)+ipoc_eups_l(:,3);
iPOC_tz_in_md=ipoc_subs(:,3)+ipoc_subs_l(:,3);
iPOC_tz2_in_md=ipoc_subbs(:,3)+ipoc_subbs_l(:,3);

if err_flag==1 % use std deviation
    iPOC_ez_in_err(1,:)=(ipoc_eups(:,a)-ipoc_eups(:,b))+(ipoc_eups_l(:,a)-ipoc_eups_l(:,b));
    iPOC_ez_in_err(2,:)=(ipoc_eups(:,a)+ipoc_eups(:,b))+(ipoc_eups_l(:,a)+ipoc_eups_l(:,b));
    iPOC_tz_in_err(1,:)=(ipoc_subs(:,a)-ipoc_subs(:,b))+(ipoc_subs_l(:,a)-ipoc_subs_l(:,b));
    iPOC_tz_in_err(2,:)=(ipoc_subs(:,a)+ipoc_subs(:,b))+(ipoc_subs_l(:,a)+ipoc_subs_l(:,b));
    iPOC_tz2_in_err(1,:)=(ipoc_subbs(:,a)-ipoc_subbs(:,b))+(ipoc_subbs_l(:,a)-ipoc_subbs_l(:,b));
    iPOC_tz2_in_err(2,:)=(ipoc_subbs(:,a)+ipoc_subbs(:,b))+(ipoc_subbs_l(:,a)+ipoc_subbs_l(:,b));

elseif err_flag==2 % use 15th and 85th percentiles
    iPOC_tz_in_err(1,:)=ipoc_subs(:,4)+ipoc_subs_l(:,4);
    iPOC_tz_in_err(2,:)=ipoc_subs(:,5)+ipoc_subs_l(:,5);
    iPOC_ez_in_err(1,:)=ipoc_eups(:,4)+ipoc_eups_l(:,4); 
    iPOC_ez_in_err(2,:)=ipoc_eups(:,5)+ipoc_eups_l(:,5); 
    iPOC_tz2_in_err(1,:)=ipoc_subbs(:,4)+ipoc_subbs_l(:,4);
    iPOC_tz2_in_err(2,:)=ipoc_subbs(:,5)+ipoc_subbs_l(:,5);

elseif err_flag==3 % use CI
    iPOC_tz_in_err(1,:)=ipoc_subs(:,6)+ipoc_subs_l(:,6);
    iPOC_tz_in_err(2,:)=ipoc_subs(:,7)+ipoc_subs_l(:,7);
    iPOC_ez_in_err(1,:)=ipoc_eups(:,6)+ipoc_eups_l(:,6);
    iPOC_ez_in_err(2,:)=ipoc_eups(:,7)+ipoc_eups_l(:,7);
    iPOC_tz2_in_err(1,:)=ipoc_subbs(:,6)+ipoc_subbs_l(:,6);
    iPOC_tz2_in_err(2,:)=ipoc_subbs(:,7)+ipoc_subbs_l(:,7);

end


% % % % % OUT
iPOC_ez_out=ipoc_eup2s(:,a)+ipoc_eup2s_l(:,a);
iPOC_tz_out=ipoc_sub2s(:,a)+ipoc_sub2s_l(:,a);
iPOC_tz2_out=ipoc_subb2s(:,a)+ipoc_subb2s_l(:,a);

iPOCs_ez_out=ipoc_eup2s(:,a);
iPOCs_tz_out=ipoc_sub2s(:,a);
iPOCs_tz2_out=ipoc_subb2s(:,a);

iPOCl_ez_out=ipoc_eup2s_l(:,a);
iPOCl_tz_out=ipoc_sub2s_l(:,a);
iPOCl_tz2_out=ipoc_subb2s_l(:,a);

iPOC_ez_out_md=ipoc_eup2s(:,3)+ipoc_eup2s_l(:,3);
iPOC_tz_out_md=ipoc_sub2s(:,3)+ipoc_sub2s_l(:,3);
iPOC_tz2_out_md=ipoc_subb2s(:,3)+ipoc_subb2s_l(:,3);

if err_flag==1 % use std deviation
    iPOC_ez_out_err(1,:)=(ipoc_eup2s(:,a)-ipoc_eup2s(:,b))+(ipoc_eup2s_l(:,a)-ipoc_eup2s_l(:,b));
    iPOC_ez_out_err(2,:)=(ipoc_eup2s(:,a)+ipoc_eup2s(:,b))+(ipoc_eup2s_l(:,a)+ipoc_eup2s_l(:,b));
    iPOC_tz_out_err(1,:)=(ipoc_sub2s(:,a)-ipoc_sub2s(:,b))+(ipoc_sub2s_l(:,a)-ipoc_sub2s_l(:,b));
    iPOC_tz_out_err(2,:)=(ipoc_sub2s(:,a)+ipoc_sub2s(:,b))+(ipoc_sub2s_l(:,a)+ipoc_sub2s_l(:,b));
    iPOC_tz2_out_err(1,:)=(ipoc_subb2s(:,a)-ipoc_subb2s(:,b))+(ipoc_subb2s_l(:,a)-ipoc_subb2s_l(:,b));
    iPOC_tz2_out_err(2,:)=(ipoc_subb2s(:,a)+ipoc_subb2s(:,b))+(ipoc_subb2s_l(:,a)+ipoc_subb2s_l(:,b));

elseif err_flag==2 % use 15th and 85th percentiles
    iPOC_tz2_out_err(1,:)=ipoc_subb2s(:,4)+ipoc_subb2s_l(:,4);
    iPOC_tz2_out_err(2,:)=ipoc_subb2s(:,5)+ipoc_subb2s_l(:,5);
    iPOC_tz_out_err(1,:)=ipoc_sub2s(:,4)+ipoc_sub2s_l(:,4);
    iPOC_tz_out_err(2,:)=ipoc_sub2s(:,5)+ipoc_sub2s_l(:,5);
    iPOC_ez_out_err(1,:)=ipoc_eup2s(:,4)+ipoc_eup2s_l(:,4);
    iPOC_ez_out_err(2,:)=ipoc_eup2s(:,5)+ipoc_eup2s_l(:,5);

elseif err_flag==3 % use CI
    iPOC_tz2_out_err(1,:)=ipoc_subb2s(:,6)+ipoc_subb2s_l(:,6);
    iPOC_tz2_out_err(2,:)=ipoc_subb2s(:,7)+ipoc_subb2s_l(:,7);
    iPOC_tz_out_err(1,:)=ipoc_sub2s(:,6)+ipoc_sub2s_l(:,6);
    iPOC_tz_out_err(2,:)=ipoc_sub2s(:,7)+ipoc_sub2s_l(:,7);
    iPOC_ez_out_err(1,:)=ipoc_eup2s(:,6)+ipoc_eup2s_l(:,6);
    iPOC_ez_out_err(2,:)=ipoc_eup2s(:,7)+ipoc_eup2s_l(:,7);

end

%% Rates
% remove any negative iPOC error
iPOC_ez_in_err(iPOC_ez_in_err<0)=0;
iPOC_tz_in_err(iPOC_tz_in_err<0)=0;
iPOC_tz2_in_err(iPOC_tz2_in_err<0)=0;

iPOC_ez_out_err(iPOC_ez_out_err<0)=0;
iPOC_tz_out_err(iPOC_tz_out_err<0)=0;
iPOC_tz2_out_err(iPOC_tz2_out_err<0)=0;

% % % %  rates, P and E
% dt=[31,28,31,30,31,30,31,31,30,31,30,31]; %days per month
dt=ones(1,12)*(356/12);
% dt=[29.5,29.5,30.5,30.5,30.5,30.5,31,30.5,30.5,30.5,30.5,31]; %avg days for each month (e.g., jan=(jan+feb)/2)
for i = 1:12
    if i <12
        % in
        P_in(i)= (iPOC_ez_in(i+1)-iPOC_ez_in(i))/dt(i);
        Ez_in(i)= (iPOC_tz_in(i+1)-iPOC_tz_in(i))/dt(i);
        Ez2_in(i)= (iPOC_tz2_in(i+1)-iPOC_tz2_in(i))/dt(i);

        P_s_in(i)= (iPOCs_ez_in(i+1)-iPOCs_ez_in(i))/dt(i);
        Ez_s_in(i)= (iPOCs_tz_in(i+1)-iPOCs_tz_in(i))/dt(i);
        Ez2_s_in(i)= (iPOCs_tz2_in(i+1)-iPOCs_tz2_in(i))/dt(i);

        P_l_in(i)= (iPOCl_ez_in(i+1)-iPOCl_ez_in(i))/dt(i);
        Ez_l_in(i)= (iPOCl_tz_in(i+1)-iPOCl_tz_in(i))/dt(i);
        Ez2_l_in(i)= (iPOCl_tz2_in(i+1)-iPOCl_tz2_in(i))/dt(i);


        Ez_in_md(i)= (iPOC_tz_in_md(i+1)-iPOC_tz_in_md(i))/dt(i);
        Ez2_in_md(i)= (iPOC_tz2_in_md(i+1)-iPOC_tz2_in_md(i))/dt(i);

        P_in_err(:,i)= (iPOC_ez_in_err(:,i+1)-iPOC_ez_in_err(:,i))./dt(i);
        Ez_in_err(:,i)= (iPOC_tz_in_err(:,i+1)-iPOC_tz_in_err(:,i))./dt(i);
        Ez2_in_err(:,i)= (iPOC_tz2_in_err(:,i+1)-iPOC_tz2_in_err(:,i))./dt(i);

        % out
        P_out(i)= (iPOC_ez_out(i+1)-iPOC_ez_out(i))/dt(i);
        Ez_out(i)= (iPOC_tz_out(i+1)-iPOC_tz_out(i))/dt(i);
        Ez2_out(i)= (iPOC_tz2_out(i+1)-iPOC_tz2_out(i))/dt(i);

        P_s_out(i)= (iPOCs_ez_out(i+1)-iPOCs_ez_out(i))/dt(i);
        Ez_s_out(i)= (iPOCs_tz_out(i+1)-iPOCs_tz_out(i))/dt(i);
        Ez2_s_out(i)= (iPOCs_tz2_out(i+1)-iPOCs_tz2_out(i))/dt(i);

        P_l_out(i)= (iPOCl_ez_out(i+1)-iPOCl_ez_out(i))/dt(i);
        Ez_l_out(i)= (iPOCl_tz_out(i+1)-iPOCl_tz_out(i))/dt(i);
        Ez2_l_out(i)= (iPOCl_tz2_out(i+1)-iPOCl_tz2_out(i))/dt(i);

        Ez_out_md(i)= (iPOC_tz_out_md(i+1)-iPOC_tz_out_md(i))/dt(i);
        Ez2_out_md(i)= (iPOC_tz2_out_md(i+1)-iPOC_tz2_out_md(i))/dt(i);

        P_out_err(:,i)= (iPOC_ez_out_err(:,i+1)-iPOC_ez_out_err(:,i))./dt(i);
        Ez_out_err(:,i)= (iPOC_tz_out_err(:,i+1)-iPOC_tz_out_err(:,i))./dt(i);
        Ez2_out_err(:,i)= (iPOC_tz2_out_err(:,i+1)-iPOC_tz2_out_err(:,i))./dt(i);

    else
% in
        P_in(i)= (iPOC_ez_in(1)-iPOC_ez_in(i))/dt(i);
        Ez_in(i)= (iPOC_tz_in(1)-iPOC_tz_in(i))/dt(i);
        Ez2_in(i)= (iPOC_tz2_in(1)-iPOC_tz2_in(i))/dt(i);

        P_s_in(i)= (iPOCs_ez_in(1)-iPOCs_ez_in(i))/dt(i);
        Ez_s_in(i)= (iPOCs_tz_in(1)-iPOCs_tz_in(i))/dt(i);
        Ez2_s_in(i)= (iPOCs_tz2_in(1)-iPOCs_tz2_in(i))/dt(i);

        P_l_in(i)= (iPOCl_ez_in(1)-iPOCl_ez_in(i))/dt(i);
        Ez_l_in(i)= (iPOCl_tz_in(1)-iPOCl_tz_in(i))/dt(i);
        Ez2_l_in(i)= (iPOCl_tz2_in(1)-iPOCl_tz2_in(i))/dt(i);

        Ez_in_md(i)= (iPOC_tz_in_md(1)-iPOC_tz_in_md(i))/dt(i);
        Ez2_in_md(i)= (iPOC_tz2_in_md(1)-iPOC_tz2_in_md(i))/dt(i);

        P_in_err(:,i)= (iPOC_ez_in_err(:,1)-iPOC_ez_in_err(:,i))./dt(i);
        Ez_in_err(:,i)= (iPOC_tz_in_err(:,1)-iPOC_tz_in_err(:,i))./dt(i);
        Ez2_in_err(:,i)= (iPOC_tz2_in_err(:,1)-iPOC_tz2_in_err(:,i))./dt(i);

% out
        P_out(i)= (iPOC_ez_out(1)-iPOC_ez_out(i))/dt(i);
        Ez_out(i)= (iPOC_tz_out(1)-iPOC_tz_out(i))/dt(i);
        Ez2_out(i)= (iPOC_tz2_out(1)-iPOC_tz2_out(i))/dt(i);

        P_s_out(i)= (iPOCs_ez_out(1)-iPOCs_ez_out(i))/dt(i);
        Ez_s_out(i)= (iPOCs_tz_out(1)-iPOCs_tz_out(i))/dt(i);
        Ez2_s_out(i)= (iPOCs_tz2_out(1)-iPOCs_tz2_out(i))/dt(i);

        P_l_out(i)= (iPOCl_ez_out(1)-iPOCl_ez_out(i))/dt(i);
        Ez_l_out(i)= (iPOCl_tz_out(1)-iPOCl_tz_out(i))/dt(i);
        Ez2_l_out(i)= (iPOCl_tz2_out(1)-iPOCl_tz2_out(i))/dt(i);

        Ez_out_md(i)= (iPOC_tz_out_md(1)-iPOC_tz_out_md(i))/dt(i);
        Ez2_out_md(i)= (iPOC_tz2_out_md(1)-iPOC_tz2_out_md(i))/dt(i);

        P_out_err(:,i)= (iPOC_ez_out_err(:,1)-iPOC_ez_out_err(:,i))./dt(i);
        Ez_out_err(:,i)= (iPOC_tz_out_err(:,1)-iPOC_tz_out_err(:,i))./dt(i);
        Ez2_out_err(:,i)= (iPOC_tz2_out_err(:,1)-iPOC_tz2_out_err(:,i))./dt(i);

    end
end


TE_in=Ez2_in./Ez_in;
TE_s_in=Ez2_s_in./Ez_s_in;
TE_l_in=Ez2_l_in./Ez_l_in;
TE_in_md=Ez2_in_md./Ez_in_md;
TE_in_err=Ez2_in_err./Ez_in_err;

TE_out=Ez2_out./Ez_out;
TE_s_out=Ez2_s_out./Ez_s_out;
TE_l_out=Ez2_l_out./Ez_l_out;
TE_out_md=Ez2_out_md./Ez_out_md;
TE_out_err=Ez2_out_err./Ez_out_err;






%% Plot

axes(ha1(3))
text(0.015,0.92,'(c)','units','normalized','fontsize',fs-4)
hold on
% plot([0 13],[0 0],'k-')
x=1:12;
y=TE_in; % set to shallowest depth horizon (100 or 200)
xx = [x;x];
yy = [min(TE_in_err);max(TE_in_err)];
alps=[ones(1,12)*0.4;ones(1,12)*0.4];
s1=surf(xx,yy,xx*0,...
    'alphadata',alps,...
    'facealpha','interp',...
    'edgecolor','none','FaceColor',colors2(7,:));

% p1=plot(1:12,TE_in,'color',colors2(7,:),'linewidth',0.75,'marker','o','MarkerFaceColor',colors2(7,:),'MarkerSize',1,'MarkerEdgeColor','none');
ps=plot(1:12,TE_s_in,'--','color',colors2(7,:),'linewidth',1.5);
pl=plot(1:12,TE_l_in,':','color',colors2(7,:),'linewidth',1.5);

% scatter(1:12,TE_in,iPOC_tz_in/115,'o','filled','MarkerFaceColor',colors2(7,:),'MarkerFaceAlpha',0.9);
scatter(1:12,TE_in,abs(Ez_in),'o','filled','MarkerFaceColor',colors2(7,:),'MarkerFaceAlpha',0.9);
% scatter(1:12,TE_in_md,50,'*','filled','MarkerEdgeColor',colors2(7,:),'MarkerFaceAlpha',0.9);

set(gca,'ticklength',[0.01 0.005])
box on
% grid on
ha1(3).YGrid='on';
ha1(3).YMinorTick='on';
set(gca,'ytick',[0:0.25:2])
set(gca,'fontsize',fs-2)
xlim([0.4 12.6])
ylim([0 1.5])
ylabel('TE^{500}')
set(gca,'xtick',[0.5:11.5],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});


%outside

axes(ha1(4))
colors2=crameri('batlow',8);
text(0.015,0.92,'(d)','units','normalized','fontsize',fs-4)

hold on
% plot([0 13],[0 0],'k-')
x=1:12;
y=TE_out; % set to shallowest depth horizon (100 or 200)
xx = [x;x];
yy = [min(TE_out_err);max(TE_out_err)];
alps=[ones(1,12)*0.4;ones(1,12)*0.4];
s1=surf(xx,yy,xx*0,...
    'alphadata',alps,...
    'facealpha','interp',...
    'edgecolor','none','FaceColor',colors2(3,:));

% p1=plot(1:12,TE_out,'color',colors2(3,:),'linewidth',0.75,'marker','o','MarkerFaceColor',colors2(3,:),'MarkerSize',1,'MarkerEdgeColor','none');
ps=plot(1:12,TE_s_out,'--','color',colors2(3,:),'linewidth',1.5);
pl=plot(1:12,TE_l_out,':','color',colors2(3,:),'linewidth',1.5);

% scatter(1:12,TE_out,iPOC_tz_out/115,'o','filled','MarkerFaceColor',colors2(3,:),'MarkerFaceAlpha',0.9);
scatter(1:12,TE_out,abs(Ez_out),'s','filled','MarkerFaceColor',colors2(3,:),'MarkerFaceAlpha',0.9);
% scatter(1:12,TE_out_md,50,'*','filled','MarkerEdgeColor',colors2(3,:),'MarkerFaceAlpha',0.9);

set(gca,'ticklength',[0.01 0.005])
box on
% grid on
ha1(4).YGrid='on';
ha1(4).YMinorTick='on';
set(gca,'ytick',[0:0.25:2])
set(gca,'fontsize',fs-2)
set(gca,'yticklabel',{''})
xlim([0.4 12.6])
ylim([0 1.5])
% ylabel('TE')
set(gca,'xtick',[0.5:11.5],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});




%% Accumulation in euphotic zone


colp=[80/255 175/255 200/255];
coltz2=[18/255 22/255 20/255];
coltz=[50/255 75/255 70/255]*0.8;

hf2=figure();
set(hf2,'Units','inches','Position', [5 5 8 4.25], 'PaperPosition', [0 0 8 4.25], 'PaperSize', [8 4.25]);
ha1=iSubplot(2,1, 'Gap', [0 0.045], 'Min', [0.065 0.03], 'Max', [0.90 0.95], 'XTickL', 'All', 'YTickL', 'All');

axes(ha1(1))
hold on
text(0.96,0.92,'(a)','units','normalized','fontsize',fs-4)

b1=bar(Ez_in,'FaceColor',coltz,'FaceAlpha',0.15,'edgecolor',[0.5 0.5 0.7],'linewidth',1.5);
errorbar(1:12,Ez_in,Ez_in-Ez_in_err(1,:),Ez_in_err(2,:)-Ez_in,'Marker','none','Color',coltz,'LineStyle','none')

x=1:12;
y=TE_out; % set to shallowest depth horizon (100 or 200)
xx = [x;x];
yy = [min(P_in_err);max(P_in_err)];
alps=[ones(1,12)*0.2;ones(1,12)*0.2];
s1=surf(xx,yy,xx*0,...
    'alphadata',alps,...
    'facealpha','interp',...
    'edgecolor','none','FaceColor',colp);


p1=plot(1:12,P_s_in,'--','Color',colp,'linewidth',2);
p2=plot(1:12,P_l_in,':','Color',colp,'linewidth',2);
p3=plot(1:12,P_in,'-','Color',colp,'linewidth',1.5);


ha1(1).YGrid='on';
set(gca,'fontsize',fs-2)
% xlim([0 13])
xlim([0.4 12.6])
ylim([-100 220])
set(gca,'ytick',[-80:40:240])
ylabel('\itE \rm[mg m^{-2} d^{-1}]')
set(gca,'xtick',[0.5:11.5],'xticklabel',{''});
ha1(1).YMinorTick='on';
title('inside LBEZ','fontsize',fs-3)


[hl,lines]=legendflex([p3 p1 p2 b1 ],{'PZ total','PZ small','PZ large','TZ total'}, 'ref', gca, ...
'anchor', {'nw', 'nw'}, 'buffer', [5 -4], 'xscale', 0.5, 'ncol',2,'nrow',3, ...
'box', 'off', 'FontSize', fs-3);
PatchInLegend = findobj(hl, 'type', 'patch');
PatchInLegend(1).FaceAlpha=0.15;

axes(ha1(2))
hold on
text(0.96,0.92,'(b)','units','normalized','fontsize',fs-4)

b1=bar(Ez_out,'FaceColor',coltz,'FaceAlpha',0.15,'edgecolor',[0.5 0.5 0.7],'linewidth',1.5);
errorbar(1:12,Ez_out,Ez_out-Ez_out_err(1,:),Ez_out_err(2,:)-Ez_out,'Marker','none','Color',coltz,'LineStyle','none')


x=1:12;
y=TE_out; % set to shallowest depth horizon (100 or 200)
xx = [x;x];
yy = [min(P_out_err);max(P_out_err)];
alps=[ones(1,12)*0.2;ones(1,12)*0.2];
s1=surf(xx,yy,xx*0,...
    'alphadata',alps,...
    'facealpha','interp',...
    'edgecolor','none','FaceColor',colp);


p1=plot(1:12,P_s_out,'--','Color',colp,'linewidth',2);
p2=plot(1:12,P_l_out,':','Color',colp,'linewidth',2);
p3=plot(1:12,P_out,'-','Color',colp,'linewidth',1.5);


ha1(2).YGrid='on';
set(gca,'fontsize',fs-2)
% xlim([0 13])
xlim([0.4 12.6])
ylim([-100 220])
set(gca,'ytick',[-80:40:240])
ylabel('\itE \rm[mg m^{-2} d^{-1}]')
set(gca,'xtick',[0.5:11.5],'xticklabel',{''});
ha1(2).YMinorTick='on';
title('outside LBEZ','fontsize',fs-3)

set(gca,'xtick',[0.5:11.5],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});

%% print

if print_flag==1
    if err_flag==1
        figure(hf)
        print(['Figures/V8/SI/LB22_TEerror_200_500_std_' scale_flag],'-dpdf','-r800')
        figure(hf2)
        print(['Figures/V8/SI/LB22_Pz_Ez_error_std_' scale_flag],'-dpdf','-r800')
    elseif err_flag==2
        figure(hf)
        print(['Figures/V8/SI/LB22_TEerror_200_500_q15_q85_' scale_flag],'-dpdf','-r800')
        figure(hf2)
        print(['Figures/V8/SI/LB22_Pz_Ez_error_q15_q85_' scale_flag],'-dpdf','-r800')
    elseif err_flag==3
        figure(hf)
        print(['Figures/V8/SI/LB22_TEerror_200_500_ci_' scale_flag],'-dpdf','-r800')
        figure(hf2)
        print(['Figures/V8/SI/LB22_Pz_Ez_error_ci_' scale_flag],'-dpdf','-r800')

    end
end
