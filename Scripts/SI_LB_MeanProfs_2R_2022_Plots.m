% Mean profiles
% Seasonal Average profiles in/out (jan-mar, april-june, july-sep,oct-dec)
% MAYBE WE CAN JUST DO A SINGLE LINE FOR EACH MONTH, using colors2, maybe
% try shaded error bar ersion
% POC, T or density, oxygen

% FLOATS ONLY
% code only calculates and saves for plotting later

% Add dotted line for EZD
%% set up


close all
clear all

curdir=cd;

addpath(genpath([curdir '/Scripts']))
addpath(genpath([curdir '/Data']))

addpath(genpath([curdir '/aux']))

print_flag=0;
save_flag=0;
scale_flag='A'; % A for 1.5r and B for 2r scaling of LBEZ
err_flag=2; %1 for std, 2 for CI

fs=13;
lw=2.3;
set(0, 'DefaultAxesFontName', 'Times');
set(0, 'DefaultTextFontName', 'Times');
alp=0.65;
mksz=7;


% month colors
c=crameri('roma',12);
colorsm=[c(9:12,:); c(1:8,:)];

tletter=['abcdefgh'];
pletter=['ijklmnop'];
%% load data
% load('LBE_BGC_POC_2010_2022_MonthlyProfiles_07-Jun-2024.mat')
% load('LBE_BGC_POC_2010_2022_MonthlyProfiles_30-Dec-2024.mat')
% load(['LBE_BGC_POC_2010_2022_MonthlyProfiles_' scale_flag '_09-Jan-2025'])
load(['LBE_BGC_POC_2010_2022_MonthlyProfiles_A_28-Jan-2025'])

rho_in=gsw_rho_t_exact(S_in,T_in,Z');
rho_out=gsw_rho_t_exact(S_out,T_out,Z');

tlims=[-1 12];
olims=[260 320];
aolims=[-20 40];
nlims=[2 20];
zlims=[-1000 0];

logpoc=1;

if logpoc==1
    poclims=[1 200];
    pocllims=[1 50];
else
    poclims=[0 150];
    pocllims=[0 45];
end

%% Plot!, 2 x 4 (each panel for season!)


sat=0.4;

%% T and POC
hf1=figure();
set(hf1,'Units','inches','Position', [5 5 8 10], 'PaperPosition', [0 0 8 10], 'PaperSize', [8 10]);
ha1=iSubplot(4,4, 'Gap', [0 0.015], 'Min', [0.045 0.04], 'Max', [0.98 0.96], 'XTickL', 'All', 'YTickL', 'All');

% T first
for i = 1:12

    %in first
    if i < 4
        axes(ha1(1))
        xlim(tlims)
        set(gca,'ytick',[-1000:200:0],'yticklabel',{'1000','800','600','400','200','0'})
        set(gca,'xticklabel',{''});
        ylabel('\itz \rm[m]')
        title('Winter')
        text(0.05,0.95,'(a)','units','normalized','fontsize',fs-3);
    elseif i < 7
        axes(ha1(2))
        xlim(tlims)
        set(gca,'xticklabel',{''});
        set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
        title('Spring')
        text(0.05,0.95,'(b)','units','normalized','fontsize',fs-3);
    elseif i < 10
        axes(ha1(3))
        xlim(tlims)
        set(gca,'xticklabel',{''});
        set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
        title('Summer')
        text(0.05,0.95,'(c)','units','normalized','fontsize',fs-3);
    else
        axes(ha1(4))
        xlim(tlims)
        set(gca,'xticklabel',{''});
        set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
        title('Autumn')
        text(0.05,0.95,'(d)','units','normalized','fontsize',fs-3);
    end

    set(gca,'ticklength',[0.02 0.02])
    grid on
    set(gca,'yminortick','on')
    set(gca,'xminortick','on')
    set(gca,'fontsize',fs)
    if err_flag==1
        h=shadedErrorBar_x(T_in(:,i),-Z,T_in_sd(:,i),'patchSaturation',sat,'lineProps',{'linewidth',lw,'color',colorsm(i,:)});
    else
        h=shadedErrorBar_x(T_in(:,i),-Z,[T_in_hi(:,i)-T_in(:,i) T_in(:,i)-T_in_low(:,i)]','patchSaturation',sat,'lineProps',{'linewidth',lw,'color',colorsm(i,:)});
    end

    h.edge(1).LineStyle='none';
    h.edge(2).LineStyle='none';
    xs=xlim;
    if err_flag==1
        h2=shadedErrorBar([10 xs(2)],-[MLD_in(:,i); MLD_in(:,i)],-[MLD_in_sd(:,i); MLD_in_sd(:,i)],'patchSaturation',sat-.2,'lineProps',{'linestyle','--','linewidth',lw-.2,'color',colorsm(i,:)});
    else
        h2=shadedErrorBar([10 xs(2)],-[MLD_in(:,i); MLD_in(:,i)],-[MLD_in_hi(:,i)-MLD_in(:,i); MLD_in(:,i)-MLD_in_low(:,i)],'patchSaturation',sat-.2,'lineProps',{'linestyle','--','linewidth',lw-.2,'color',colorsm(i,:)});
    end

    h2.edge(1).LineStyle='none';
    h2.edge(2).LineStyle='none';
    uistack(h.mainLine,'top')
    ylim(zlims);

    %out
    if i < 4
        axes(ha1(5))
        xlim(tlims)
        % set(gca,'xtick',[0:10:30],'xticklabel',{'0','10','20','30'});
        set(gca,'ytick',[-1000:200:0],'yticklabel',{'1000','800','600','400','200','0'})
        ylabel('\itz \rm[m]')
        xlabel('T [\circC]')
        text(0.05,0.95,'(e)','units','normalized','fontsize',fs-3);
    elseif i < 7
        axes(ha1(6))
        xlim(tlims)
        % set(gca,'xtick',[0:50:150],'xticklabel',{'0','50','100','150'});
        set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
        xlabel('T [\circC]')
        text(0.05,0.95,'(f)','units','normalized','fontsize',fs-3);
    elseif i < 10
        axes(ha1(7))
        xlim(tlims)
        % set(gca,'xtick',[0:50:150],'xticklabel',{'0','50','100','150'});
        set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
        xlabel('T [\circC]')
        text(0.05,0.95,'(g)','units','normalized','fontsize',fs-3);
    else
        axes(ha1(8))
        xlim(tlims)
        % set(gca,'xtick',[0:50:150],'xticklabel',{'0','50','100','150'});
        set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
        xlabel('T [\circC]')
        text(0.05,0.95,'(h)','units','normalized','fontsize',fs-3);
    end

    set(gca,'ticklength',[0.02 0.02])
    grid on
    set(gca,'yminortick','on')
    set(gca,'xminortick','on')
    set(gca,'fontsize',fs)
    if err_flag==1
        h=shadedErrorBar_x(T_out(:,i),-Z,T_out_sd(:,i),'patchSaturation',sat,'lineProps',{'linewidth',lw,'color',colorsm(i,:)});
    else
        h=shadedErrorBar_x(T_out(:,i),-Z,[T_out_hi(:,i)-T_out(:,i) T_out(:,i)-T_out_low(:,i)]','patchSaturation',sat,'lineProps',{'linewidth',lw,'color',colorsm(i,:)});
    end

    h.edge(1).LineStyle='none';
    h.edge(2).LineStyle='none';
    xs=xlim;
    if err_flag==1
        h2=shadedErrorBar([10 xs(2)],-[MLD_out(:,i); MLD_out(:,i)],-[MLD_out_sd(:,i); MLD_out_sd(:,i)],'patchSaturation',sat-.2,'lineProps',{'linestyle','--','linewidth',lw-.2,'color',colorsm(i,:)});
    else
        h2=shadedErrorBar([10 xs(2)],-[MLD_out(:,i); MLD_out(:,i)],-[MLD_out_hi(:,i)-MLD_out(:,i); MLD_out(:,i)-MLD_out_low(:,i)],'patchSaturation',sat-.2,'lineProps',{'linestyle','--','linewidth',lw-.2,'color',colorsm(i,:)});
    end

    h2.edge(1).LineStyle='none';
    h2.edge(2).LineStyle='none';
    uistack(h.mainLine,'top')
    ylim(zlims);
    clear xs
end

axes(ha1(1));
text(-0.25,1.1,'inside LBEZ','units','normalized','fontsize',fs-3,'fontweight','bold')
hl=plot(-1,-1,'k--','linewidth',lw);
[hl,lines]=legendflex(hl,{'Mean MLD'}, 'ref', gca, ...
    'anchor', {'se', 'se'}, 'buffer', [0 -25], 'xscale', 0.5, 'ncol',1,'nrow',1, ...
    'box', 'off', 'FontSize', fs-3);

axes(ha1(5));
text(-0.25,1.1,'outside LBEZ','units','normalized','fontsize',fs-3,'fontweight','bold')


%% Repeat for 2R

load(['LBE_BGC_POC_2010_2022_MonthlyProfiles_B_29-Jan-2025'])
% % % % % % % T
for i = 1:12

    %in first
    if i < 4
        axes(ha1(9))
        xlim(tlims)
        set(gca,'ytick',[-1000:200:0],'yticklabel',{'1000','800','600','400','200','0'})
        set(gca,'xticklabel',{''});
        ylabel('\itz \rm[m]')
        text(0.88,0.07,'(i)','units','normalized','fontsize',fs-3);
        % title('Winter')
    elseif i < 7
        axes(ha1(10))
        xlim(tlims)
        set(gca,'xticklabel',{''});
        set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
        % title('Spring')
        text(0.88,0.07,'(j)','units','normalized','fontsize',fs-3);
    elseif i < 10
        axes(ha1(11))
        xlim(tlims)
        set(gca,'xticklabel',{''});
        set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
        % title('Summer')
        text(0.88,0.07,'(k)','units','normalized','fontsize',fs-3);
    else
        axes(ha1(12))
        xlim(tlims)

        set(gca,'xticklabel',{''});
        set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
        % title('Autumn')
        text(0.88,0.07,'(l)','units','normalized','fontsize',fs-3);
    end

    if err_flag==1
        h=shadedErrorBar_x(T_in(:,i),-Z,T_in_sd(:,i),'patchSaturation',sat,'lineProps',{'linewidth',lw,'color',colorsm(i,:)});
    else
        h=shadedErrorBar_x(T_in(:,i),-Z,[T_in_hi(:,i)-T_in(:,i) T_in(:,i)-T_in_low(:,i)]','patchSaturation',sat,'lineProps',{'linewidth',lw,'color',colorsm(i,:)});
    end

    h.edge(1).LineStyle='none';
    h.edge(2).LineStyle='none';
    xs=xlim;
    if err_flag==1
        h2=shadedErrorBar([10 xs(2)],-[MLD_in(:,i); MLD_in(:,i)],-[MLD_in_sd(:,i); MLD_in_sd(:,i)],'patchSaturation',sat-.2,'lineProps',{'linestyle','--','linewidth',lw-.2,'color',colorsm(i,:)});
    else
        h2=shadedErrorBar([10 xs(2)],-[MLD_in(:,i); MLD_in(:,i)],-[MLD_in_hi(:,i)-MLD_in(:,i); MLD_in(:,i)-MLD_in_low(:,i)],'patchSaturation',sat-.2,'lineProps',{'linestyle','--','linewidth',lw-.2,'color',colorsm(i,:)});
    end

    h2.edge(1).LineStyle='none';
    h2.edge(2).LineStyle='none';
    uistack(h.mainLine,'top')
    ylim(zlims);

    %out
    if i < 4
        axes(ha1(13))
        xlim(tlims)

        set(gca,'ytick',[-1000:200:0],'yticklabel',{'1000','800','600','400','200','0'})
        ylabel('\itz \rm[m]')
        xlabel('T [\circC]')
        text(0.88,0.07,'(m)','units','normalized','fontsize',fs-3);
    elseif i < 7
        axes(ha1(14))
        xlim(tlims)


        set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
        xlabel('T [\circC]')
        text(0.88,0.07,'(n)','units','normalized','fontsize',fs-3);
    elseif i < 10
        axes(ha1(15))
        xlim(tlims)


        set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
        xlabel('T [\circC]')
        text(0.88,0.07,'(o)','units','normalized','fontsize',fs-3);
    else
        axes(ha1(16))
        xlim(tlims)


        set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
        xlabel('T [\circC]')
        text(0.88,0.07,'(p)','units','normalized','fontsize',fs-3);
    end

    set(gca,'ticklength',[0.02 0.02])
    grid on
    set(gca,'yminortick','on')
    set(gca,'xminortick','on')
    set(gca,'fontsize',fs)
    if err_flag==1
        h=shadedErrorBar_x(T_out(:,i),-Z,T_out_sd(:,i),'patchSaturation',sat,'lineProps',{'linewidth',lw,'color',colorsm(i,:)});
    else
        h=shadedErrorBar_x(T_out(:,i),-Z,[T_out_hi(:,i)-T_out(:,i) T_out(:,i)-T_out_low(:,i)]','patchSaturation',sat,'lineProps',{'linewidth',lw,'color',colorsm(i,:)});
    end

    h.edge(1).LineStyle='none';
    h.edge(2).LineStyle='none';
    xs=xlim;
    if err_flag==1
        h2=shadedErrorBar([10 xs(2)],-[MLD_out(:,i); MLD_out(:,i)],-[MLD_out_sd(:,i); MLD_out_sd(:,i)],'patchSaturation',sat-.2,'lineProps',{'linestyle','--','linewidth',lw-.2,'color',colorsm(i,:)});
    else
        h2=shadedErrorBar([10 xs(2)],-[MLD_out(:,i); MLD_out(:,i)],-[MLD_out_hi(:,i)-MLD_out(:,i); MLD_out(:,i)-MLD_out_low(:,i)],'patchSaturation',sat-.2,'lineProps',{'linestyle','--','linewidth',lw-.2,'color',colorsm(i,:)});
    end

    h2.edge(1).LineStyle='none';
    h2.edge(2).LineStyle='none';
    uistack(h.mainLine,'top')
    ylim(zlims);
    clear xs
end

axes(ha1(9));
text(-0.25,1.1,'inside LBEZ','units','normalized','fontsize',fs-3,'fontweight','bold')

hl=plot(-1,-1,'k:','linewidth',lw);
[hl,lines]=legendflex(hl,{'Mean EZD'}, 'ref', gca, ...
    'anchor', {'se', 'se'}, 'buffer', [0 -25], 'xscale', 0.5, 'ncol',1,'nrow',1, ...
    'box', 'off', 'FontSize', fs-3);

axes(ha1(13));
text(-0.25,1.1,'outside LBEZ','units','normalized','fontsize',fs-3,'fontweight','bold')

for i = [1,5]
yl=get(ha1(i),'ylabel');
yl.Position(1)=yl.Position(1)-yl.Position(1)/5;
clear yl
end

for i = [9,13]
yl=get(ha1(i),'ylabel');
yl.Position(1)=yl.Position(1)+yl.Position(1)/3;
clear yl
end

for i = 5:8
xl=get(ha1(i),'xlabel');
xl.Position(2)=xl.Position(2)-xl.Position(2)/25;
clear xl
end

for i = 13:16
xl=get(ha1(i),'xlabel');
xl.Position(2)=xl.Position(2)-xl.Position(2)/25;
clear xl
end

for i = 1:8
    ha1(i).Position(2)=ha1(i).Position(2)+0.015;
end
for i = 9:16
    ha1(i).Position(2)=ha1(i).Position(2)-0.015;
end


axes(ha1(7));
for i = 1:12
        s(i)=scatter(i,i,50,colorsm(i,:),'filled');
end


lgd_text={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
[hl,lines]=legendflex([s],lgd_text, 'ref', gca, ...
'anchor', {'nw', 'nw'}, 'buffer', [20 30], 'xscale', 0.5, 'ncol',6,'nrow',2, ...
'box', 'off', 'FontSize', fs-3);


%% print
if print_flag==1
    %     exportgraphics(hf1,['Figures/OSM24_out_MeanProfs_POC.pdf'],'ContentType','vector');
    %     exportgraphics(hf2,['Figures/OSM24_out_MeanProfs_Temp.pdf'],'ContentType','vector');
    %     exportgraphics(hf3,['Figures/OSM24_out_MeanProfs_DOXY.pdf'],'ContentType','vector');
    %     exportgraphics(hf4,['Figures/OSM24_out_MeanProfs_Nitrate.pdf'],'ContentType','vector');
    if err_flag==1
        figure(hf1);
        print(['Figures/V7/2_LB22_MonthlyProfiles_T_POCs_sd_' scale_flag],'-dpdf','-r800')
    elseif err_flag==2
        figure(hf1);
        print(['Figures/V7/2_LB22_MonthlyProfiles_T_POCs_ci_' scale_flag],'-dpdf','-r800')
    end

    figure(hf2);
    print(['Figures/V7/3_LB22_TimeSections_' scale_flag],'-dpdf','-r800')

    figure(hf4);
    print(['Figures/V7/SI/LB22_TimeSections_PDif_' scale_flag],'-dpdf','-r800')

    % figure(hf3);
    % print('Figures/LB22_TimeSections_ttest','-dpdf','-r800')
end
