% Mean profiles
% Seasonal Average profiles in/out (jan-mar, april-june, july-sep,oct-dec)
% MAYBE WE CAN JUST DO A SINGLE LINE FOR EACH MONTH, using colors2, maybe
% try shaded error bar ersion
% POC, T or density, oxygen

% FLOATS ONLY
% code only calculates and saves for plotting later

%% set up


close all
clear all

curdir=cd;

addpath(genpath([curdir '/Scripts']))
addpath(genpath([curdir '/Data']))

addpath(genpath([curdir '/aux']))

print_flag=0;
save_flag=0;
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
load('LBE_BGC_POC_2010_2022_MonthlyProfiles_07-Jun-2024.mat')


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
    h=shadedErrorBar_x(T_in(:,i),-Z,T_in_sd(:,i),'patchSaturation',sat,'lineProps',{'linewidth',lw,'color',colorsm(i,:)});
    h.edge(1).LineStyle='none';
    h.edge(2).LineStyle='none';
    xs=xlim;
    h2=shadedErrorBar([10 xs(2)],-[MLD_in(:,i); MLD_in(:,i)],-[MLD_in_sd(:,i); MLD_in_sd(:,i)],'patchSaturation',sat-.2,'lineProps',{'linestyle','--','linewidth',lw-.2,'color',colorsm(i,:)});
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
    h=shadedErrorBar_x(T_out(:,i),-Z,T_out_sd(:,i),'patchSaturation',sat,'lineProps',{'linewidth',lw,'color',colorsm(i,:)});
    h.edge(1).LineStyle='none';
    h.edge(2).LineStyle='none';
    xs=xlim;
    h2=shadedErrorBar([10 xs(2)],-[MLD_out(:,i); MLD_out(:,i)],-[MLD_out_sd(:,i); MLD_out_sd(:,i)],'patchSaturation',sat-.2,'lineProps',{'linestyle','--','linewidth',lw-.2,'color',colorsm(i,:)});
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


% % % % % % % POC
for i = 1:12

    %in first
    if i < 4
        axes(ha1(9))
        xlim(poclims)
        set(gca,'ytick',[-1000:200:0],'yticklabel',{'1000','800','600','400','200','0'})
        if logpoc==1
            set(gca,'xtick',[1 10 100],'xticklabel',{''},'xscale','log','XMinorGrid','off');
        else
            set(gca,'xtick',[0:50:150],'xticklabel',{''});
        end
        ylabel('\itz \rm[m]')
        text(0.88,0.07,'(i)','units','normalized','fontsize',fs-3);
        % title('Winter')
    elseif i < 7
        axes(ha1(10))
        xlim(poclims)
        if logpoc==1
            set(gca,'xtick',[1 10 100],'xticklabel',{''},'xscale','log','XMinorGrid','off');
        else
            set(gca,'xtick',[0:50:150],'xticklabel',{''});
        end
        set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
        % title('Spring')
        text(0.88,0.07,'(j)','units','normalized','fontsize',fs-3);
    elseif i < 10
        axes(ha1(11))
        xlim(poclims)
        if logpoc==1
            set(gca,'xtick',[1 10 100],'xticklabel',{''},'xscale','log','XMinorGrid','off');
        else
            set(gca,'xtick',[0:50:150],'xticklabel',{''});
        end
        set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
        % title('Summer')
        text(0.88,0.07,'(k)','units','normalized','fontsize',fs-3);
    else
        axes(ha1(12))
        xlim(poclims)
        if logpoc==1
            set(gca,'xtick',[1 10 100],'xticklabel',{''},'xscale','log','XMinorGrid','off');
        else
            set(gca,'xtick',[0:50:150],'xticklabel',{''});
        end
        set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
        % title('Autumn')
        text(0.88,0.07,'(l)','units','normalized','fontsize',fs-3);
    end

    set(gca,'ticklength',[0.02 0.02])
    grid on
    set(gca,'yminortick','on')
    set(gca,'fontsize',fs)
    h=shadedErrorBar_x_log(POC_in(:,i),-Z,POC_in_sd(:,i),'patchSaturation',sat,'lineProps',{'linewidth',lw,'color',colorsm(i,:)});
    h.edge(1).LineStyle='none';
    h.edge(2).LineStyle='none';
    xs=xlim;
    % h2=shadedErrorBar([xs(2)*2/3 xs(2)],-[EZD_in(:,i); EZD_in(:,i)],-[EZD_in_sd(:,i); EZD_in_sd(:,i)],'patchSaturation',sat-.2,'lineProps',{'linestyle',':','linewidth',lw-.2,'color',colorsm(i,:)});
    h2=shadedErrorBar([100 xs(2)],-[EZD_in(:,i); EZD_in(:,i)],-[EZD_in_sd(:,i); EZD_in_sd(:,i)],'patchSaturation',sat-.2,'lineProps',{'linestyle',':','linewidth',lw-.2,'color',colorsm(i,:)});
    h2.edge(1).LineStyle='none';
    h2.edge(2).LineStyle='none';
    uistack(h.mainLine,'top')
    ylim(zlims);

    %out
    if i < 4
        axes(ha1(13))
        xlim(poclims)
        if logpoc==1
            set(gca,'xtick',[1 10 100],'xticklabel',{'1','10','100'},'xscale','log','XMinorGrid','off');
        else
            set(gca,'xtick',[0:50:150],'xticklabel',{'0','50','100','150'});
        end

        set(gca,'ytick',[-1000:200:0],'yticklabel',{'1000','800','600','400','200','0'})
        ylabel('\itz \rm[m]')
        xlabel('POC\it_s\rm [mg m^{-3}]')
text(0.88,0.07,'(m)','units','normalized','fontsize',fs-3);
    elseif i < 7
        axes(ha1(14))
        xlim(poclims)
         if logpoc==1
            set(gca,'xtick',[1 10 100],'xticklabel',{'1','10','100'},'xscale','log','XMinorGrid','off');
        else
            set(gca,'xtick',[0:50:150],'xticklabel',{'0','50','100','150'});
        end

        set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
        xlabel('POC\it_s\rm [mg m^{-3}]')
text(0.88,0.07,'(n)','units','normalized','fontsize',fs-3);
    elseif i < 10
        axes(ha1(15))
        xlim(poclims)
        if logpoc==1
            set(gca,'xtick',[1 10 100],'xticklabel',{'1','10','100'},'xscale','log','XMinorGrid','off');
        else
            set(gca,'xtick',[0:50:150],'xticklabel',{'0','50','100','150'});
        end

        set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
        xlabel('POC\it_s\rm [mg m^{-3}]')
text(0.88,0.07,'(o)','units','normalized','fontsize',fs-3);
    else
        axes(ha1(16))
        xlim(poclims)
        if logpoc==1
            set(gca,'xtick',[1 10 100],'xticklabel',{'1','10','100'},'xscale','log','XMinorGrid','off');
        else
            set(gca,'xtick',[0:50:150],'xticklabel',{'0','50','100','150'});
        end

        set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
        xlabel('POC\it_s\rm [mg m^{-3}]')
text(0.88,0.07,'(p)','units','normalized','fontsize',fs-3);
    end

    set(gca,'ticklength',[0.02 0.02])
    grid on
    set(gca,'yminortick','on')
    set(gca,'fontsize',fs)
    h=shadedErrorBar_x_log(POC_out(:,i),-Z,POC_out_sd(:,i),'patchSaturation',sat,'lineProps',{'linewidth',lw,'color',colorsm(i,:)});
    h.edge(1).LineStyle='none';
    h.edge(2).LineStyle='none';
    xs=xlim;
    % h2=shadedErrorBar([xs(2)*2/3 xs(2)],-[EZD_out(:,i); EZD_out(:,i)],-[EZD_out_sd(:,i); EZD_out_sd(:,i)],'patchSaturation',sat-.2,'lineProps',{'linestyle',':','linewidth',lw-.2,'color',colorsm(i,:)});
    h2=shadedErrorBar([100 xs(2)],-[EZD_out(:,i); EZD_out(:,i)],-[EZD_out_sd(:,i); EZD_out_sd(:,i)],'patchSaturation',sat-.2,'lineProps',{'linestyle',':','linewidth',lw-.2,'color',colorsm(i,:)});
    h2.edge(1).LineStyle='none';
    h2.edge(2).LineStyle='none';
    uistack(h.mainLine,'top')
    ylim(zlims);
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
%% Section plots, 3 x 1 plots for in, out, and difference

% aoucols=crameri('bam',240);
% aoucols=(aoucols(61:end,:));
% aoucols=crameri('bam',28);
% aoucols=(aoucols(8:end,:));

aoucols=crameri('bam',20);
aoucols=(aoucols(6:end,:));

deltacols=crameri('vik',16);
cols=crameri('batlow',14);
colst=crameri('lapaz',13);

sat=-0.2; %brightness scale
aoucols=brighten(aoucols,sat);
deltacols=brighten(deltacols,sat);
cols=brighten(cols,0.1);
fs=12;


% POC
hf2=figure();
set(hf2,'Units','inches','Position', [5 5 8 11], 'PaperPosition', [0 0 8 11], 'PaperSize', [8 11]);
ha1=iSubplot(4,3, 'Gap', [0 0.01], 'Min', [0.04 0.01], 'Max', [0.98 1], 'XTickL', 'All', 'YTickL', 'All');

axes(ha1(4))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{'1000','800','600','400','200','0'})
ylabel('\itz \rm[m]')
set(gca,'XTickLabelRotation',45)

s=surface(1:12,-Z,log10(POC_in));
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,log10(poclims));
colormap(ha1(4),cols)
c=colorbar;
c.Location='northoutside';
text(0.55,1.4,'POC_{\its\rm} inside [mg m^{-3}]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(d)','units','normalized','fontsize',fs-2)

c.Ticks=log10([1 10 100]);
c.TickLabels={'1','10','100'};
plot3(1:12,-zPROD_in,max(log10(POC_in)),'-.','linewidth',1,'color',[0.7 0.7 0.7])
c.FontSize=fs-2;
% text(1,-zPROD_in(1),'PZD','units','normalized','fontsize',fs-2,'color',[0.5 0.5 0.5])

axes(ha1(5))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{''})

s=surface(1:12,-Z,log10(POC_out));
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,log10(poclims));
colormap(ha1(5),cols)
c=colorbar;
c.Location='northoutside';
c.Ticks=log10([1 10 100]);
c.TickLabels={'1','10','100'};
text(0.55,1.4,'POC_{\its\rm} outside [mg m^{-3}]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(e)','units','normalized','fontsize',fs-2)

set(gca,'XTickLabelRotation',45)
c.FontSize=fs-2;

plot3(1:12,-zPROD_out,max(log10(POC_out)),'-.','linewidth',1,'color',[0.7 0.7 0.7])
axes(ha1(6))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{''})

s=surface(1:12,-Z,POC_in-POC_out);
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
% caxis(gca,[-25 25]);
caxis(gca,[-20 20]);
colormap(ha1(6),deltacols)
c=colorbar;
c.Ticks=-15:5:15;
c.Location='northoutside';
text(0.55,1.4,'\Delta POC_{\its\rm} [mg m^{-3}]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(f)','units','normalized','fontsize',fs-2)

set(gca,'XTickLabelRotation',45)
c.FontSize=fs-2;

plot3(1:12,-zPROD_in,max(POC_in-POC_out),'-.','linewidth',1,'color',[0.4 0.4 0.4])

% % % Temp
% hf7=figure();
% set(hf7,'Units','inches','Position', [5 5 12 6], 'PaperPosition', [0 0 12 6], 'PaperSize', [12 6]);
% ha1=iSubplot(1,3, 'Gap', [0 0.02], 'Min', [0.04 0.05], 'Max', [0.99 0.96], 'XTickL', 'All', 'YTickL', 'All');

axes(ha1(1))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{'1000','800','600','400','200','0'})
ylabel('\itz \rm[m]')

s=surface(1:12,-Z,T_in);
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,tlims);
colormap(ha1(1),colst)
c=colorbar;
c.Location='northoutside';
text(0.55,1.4,'T inside [\circC]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(a)','units','normalized','fontsize',fs-2)

plot3(1:12,-MLD_in,max(T_in),'--','linewidth',1,'color',[0.7 0.7 0.7])
set(gca,'XTickLabelRotation',45)
c.FontSize=fs-2;
text(1,-MLD_in(1),'MLD','units','normalized','fontsize',fs-2,'color',[0.4 0.4 0.4])

axes(ha1(2))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{''})

s=surface(1:12,-Z,T_out);
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,tlims);
colormap(ha1(2),colst)
c=colorbar;
c.Location='northoutside';
text(0.55,1.4,'T outside [\circC]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(b)','units','normalized','fontsize',fs-2)

plot3(1:12,-MLD_out,max(T_out),'--','linewidth',1,'color',[0.7 0.7 0.7])
set(gca,'XTickLabelRotation',45)
c.FontSize=fs-2;

axes(ha1(3))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{''})

s=surface(1:12,-Z,T_in-T_out);
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,[-4 4]);
colormap(ha1(3),deltacols)
c=colorbar;
c.Ticks=-3:1:3;
c.Location='northoutside';
text(0.55,1.4,'\Delta T [\circC]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(c)','units','normalized','fontsize',fs-2)

plot3(1:12,-MLD_in,max(T_in-T_out),'--','linewidth',1,'color',[0.4 0.4 0.4])
c.FontSize=fs-2;

set(gca,'XTickLabelRotation',45)


% % % % AOU
% hf8=figure();
% set(hf8,'Units','inches','Position', [5 5 12 6], 'PaperPosition', [0 0 12 6], 'PaperSize', [12 6]);
% ha1=iSubplot(1,3, 'Gap', [0 0.02], 'Min', [0.04 0.05], 'Max', [0.99 0.96], 'XTickL', 'All', 'YTickL', 'All');
% 


axes(ha1(10))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{'1000','800','600','400','200','0'})
ylabel('\itz \rm[m]')
set(gca,'XTickLabelRotation',45)

s=surface(1:12,-Z,AOU_in);
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,aolims);
colormap(ha1(10),aoucols)
c=colorbar;
c.Location='northoutside';
text(0.55,1.4,'AOU inside [\muM]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(j)','units','normalized','fontsize',fs-2)

plot3(1:12,-zPROD_in,max(AOU_in),'-.','linewidth',1,'color',[0.4 0.4 0.4])
c.FontSize=fs-2;
text(1,-zPROD_in(1),'PZD','units','normalized','fontsize',fs-2,'color',[0.5 0.5 0.5])

axes(ha1(11))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
set(gca,'XTickLabelRotation',45)

s=surface(1:12,-Z,AOU_out);
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,aolims);
colormap(ha1(11),aoucols)
c=colorbar;
c.Location='northoutside';
text(0.55,1.4,'AOU outside [\muM]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(k)','units','normalized','fontsize',fs-2)

plot3(1:12,-zPROD_out,max(AOU_out),'-.','linewidth',1,'color',[0.4 0.4 0.4])
c.FontSize=fs-2;

axes(ha1(12))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
set(gca,'XTickLabelRotation',45)

s=surface(1:12,-Z,AOU_in-AOU_out);
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,[-20 20]);
colormap(ha1(12),deltacols)
c=colorbar;
c.Ticks=-15:5:15;
c.Location='northoutside';
text(0.55,1.4,'\Delta AOU [\muM]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(l)','units','normalized','fontsize',fs-2)

plot3(1:12,-zPROD_in,max(AOU_in-AOU_out),'-.','linewidth',1,'color',[0.4 0.4 0.4])
% plot3(1:12,-zPROD_in-200,max(AOU_in-AOU_out),'-','linewidth',0.2,'color',[0.4 0.4 0.4])

c.FontSize=fs-2;


% hf9=figure();
% set(hf9,'Units','inches','Position', [5 5 12 6], 'PaperPosition', [0 0 12 6], 'PaperSize', [12 6]);
% ha1=iSubplot(1,3, 'Gap', [0 0.02], 'Min', [0.04 0.05], 'Max', [0.99 0.96], 'XTickL', 'All', 'YTickL', 'All');

axes(ha1(7))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{'1000','800','600','400','200','0'})
ylabel('\itz \rm[m]')
set(gca,'XTickLabelRotation',45)

s=surface(1:12,-Z,log10(POCl_in));
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,log10(pocllims));
colormap(ha1(7),cols)
c=colorbar;
c.Location='northoutside';
c.Ticks=log10([1 5 50]);
c.TickLabels={'1','5','50'};
text(0.55,1.4,'POC_{\itl\rm} inside [mg m^{-3}]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(g)','units','normalized','fontsize',fs-2)

c.FontSize=fs-2;
text(1,-zPROD_in(1),'PZD','units','normalized','fontsize',fs-2,'color',[0.5 0.5 0.5])

plot3(1:12,-zPROD_in,max(log10(POCl_in)),'-.','linewidth',1,'color',[0.7 0.7 0.7])


axes(ha1(8))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
set(gca,'XTickLabelRotation',45)

s=surface(1:12,-Z,log10(POCl_out));
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,log10(pocllims));
% colormap(cols)
colormap(ha1(8),cols)
c=colorbar;
c.Location='northoutside';
c.Ticks=log10([1 5 50]);
c.TickLabels={'1','5','50'};
text(0.55,1.4,'POC_{\itl\rm} outside [mg m^{-3}]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(h)','units','normalized','fontsize',fs-2)

plot3(1:12,-zPROD_out,max(log10(POCl_out)),'-.','linewidth',1,'color',[0.7 0.7 0.7])
c.FontSize=fs-2;


axes(ha1(9))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
set(gca,'XTickLabelRotation',45)

s=surface(1:12,-Z,POCl_in-POCl_out);
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
% caxis(gca,[-20 20]);
caxis(gca,[-12 12]);
colormap(ha1(9),deltacols)
c=colorbar;
% c.Ticks=-15:7.5:15;
c.Ticks=-9:3:9;
c.Location='northoutside';
text(0.55,1.4,'\Delta POC_{\itl\rm} [mg m^{-3}]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(i)','units','normalized','fontsize',fs-2)

c.FontSize=fs-2;

plot3(1:12,-zPROD_in,max(POCl_in-POCl_out),'-.','linewidth',1,'color',[0.4 0.4 0.4])
% plot3(1:12,-zPROD_in-200,max(POCl_in-POCl_out),'-','linewidth',0.2,'color',[0.4 0.4 0.4])


for i = [1,4,7,10]
yl=get(ha1(i),'ylabel');
yl.Position(1)=yl.Position(1)-yl.Position(1)/5;
clear yl
end

for i = 1:12
    ha1(i).Position(4)=0.135;
end

for i = 1:3
    ha1(i).Position(2)=ha1(i).Position(2)+0.01;
end

for i = 4:6
    ha1(i).Position(2)=ha1(i).Position(2)+0.005;
end

for i = [1 4 7 10]
    c1=get(ha1(i),'colorbar');
    c2=get(ha1(i+1),'colorbar');
    c1.Position(2:3)=c2.Position(2:3);
end
%% Month legend plot

hf11=figure();
set(hf11,'Units','inches','Position', [5 5 4 4], 'PaperPosition', [0 0 4 4], 'PaperSize', [4 4]);
ha11=iSubplot(1,1, 'Gap', [0.02 0.00], 'Min', [0.07 0.02], 'Max', [0.98 0.98], 'XTickL', 'All', 'YTickL', 'All');

hold on
for i = 1:12
        s(i)=scatter(i,i,50,colorsm(i,:),'filled');
end
xlim([i i+1])
ylim([i i+1])


lgd_text={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
[hl,lines]=legendflex([s],lgd_text, 'ref', gca, ...
'anchor', {'w', 'w'}, 'buffer', [30 0], 'xscale', 0.5, 'ncol',4,'nrow',3, ...
'box', 'off', 'FontSize', fs);

%% statistical signifigance
for i = 1:12
for ii = 1:80
[~,p_aou(ii,i),~,~]=ttest2(AOU_in_all{i}(ii,:),AOU_out_all{i}(ii,:),'Vartype','unequal');
[~,p_t(ii,i),~,~]=ttest2(T_in_all{i}(ii,:),T_out_all{i}(ii,:),'Vartype','unequal');
[~,p_pocs(ii,i),~,~]=ttest2(POC_in_all{i}(ii,:),POC_out_all{i}(ii,:),'Vartype','unequal');
[~,p_pocl(ii,i),~,~]=ttest2(POCl_in_all{i}(ii,:),POCl_out_all{i}(ii,:),'Vartype','unequal');

end
end


cols=flip(crameri('lajolla'));
hf3=figure();
set(hf3,'Units','inches','Position', [5 5 6 7], 'PaperPosition', [0 0 6 7], 'PaperSize', [6 7]);
ha1=iSubplot(2,2, 'Gap', [0 0.03], 'Min', [0.07 0.03], 'Max', [0.98 1], 'XTickL', 'All', 'YTickL', 'All');

axes(ha1(1))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{'1000','800','600','400','200','0'})
ylabel('\itz \rm[m]')
set(gca,'XTickLabelRotation',45)

s=surface(1:12,-Z,1-p_t);
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,[0 1]);
colormap(ha1(1),cols)
c=colorbar;
c.Location='northoutside';
title('T')

axes(ha1(2))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
% ylabel('\itz \rm[m]')
set(gca,'XTickLabelRotation',45)

s=surface(1:12,-Z,1-p_aou);
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,[0 1]);
colormap(ha1(2),cols)
c=colorbar;
c.Location='northoutside';
title('AOU')

axes(ha1(3))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{'1000','800','600','400','200','0'})
ylabel('\itz \rm[m]')
set(gca,'XTickLabelRotation',45)

s=surface(1:12,-Z,1-p_pocs);
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,[0 1]);
colormap(ha1(3),cols)
c=colorbar;
c.Location='northoutside';
title('POCs')

axes(ha1(4))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{''})
% ylabel('\itz \rm[m]')
set(gca,'XTickLabelRotation',45)

s=surface(1:12,-Z,1-p_pocl);
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,[0 1]);
colormap(ha1(4),cols)
c=colorbar;
c.Location='northoutside';
title('POCl')

axes(ha1(1))
text(1.4,1.1,'two sample \itt\rm\bf-test, 1–\itp ','units','normalized','fontsize',fs-2,'fontweight','bold','HorizontalAlignment','right')


%% Section plots just percent difference

% aoucols=crameri('bam',240);
% aoucols=(aoucols(61:end,:));
aoucols=crameri('bam',28);
aoucols=(aoucols(8:end,:));

deltacols=crameri('vik',16);
cols=crameri('batlow',16);
colst=crameri('lapaz',16);

% sat=-0.2; %brightness scale
% aoucols=brighten(aoucols,sat);
% deltacols=brighten(deltacols,sat);
% cols=brighten(cols,-0.1);
fs=12;


hf4=figure();
set(hf4,'Units','inches','Position', [5 5 4 9.5], 'PaperPosition', [0 0 4 9.5], 'PaperSize', [4 9.5]);
ha1=iSubplot(4,1, 'Gap', [0 0.05], 'Min', [0.1 0.01], 'Max', [0.98 0.98], 'XTickL', 'All', 'YTickL', 'All');

axes(ha1(2))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{'1000','800','600','400','200','0'})
ylabel('\itz \rm[m]')

s=surface(1:12,-Z,100*(POC_in-POC_out)./POC_out);
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,[-100 100]);
colormap(ha1(2),deltacols)
c=colorbar;
c.Location='eastoutside';
title('\Delta POC_{\its\rm} [%]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
% text(0,1.4,'(f)','units','normalized','fontsize',fs-2)

set(gca,'XTickLabelRotation',45)
c.FontSize=fs-2;

plot3(1:12,-zPROD_in,max(100*(POC_in-POC_out)./POC_out),'-.','linewidth',1,'color',[0.4 0.4 0.4])

% % % Temp

axes(ha1(1))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{'1000','800','600','400','200','0'})
ylabel('\itz \rm[m]')

s=surface(1:12,-Z,100*(T_in-T_out)./T_out);
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,[-100 100]);
colormap(ha1(1),deltacols)
c=colorbar;
% c.Ticks=-3:1.5:3;
c.Location='eastoutside';
title('\Delta T [%]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
% text(0,1.4,'(c)','units','normalized','fontsize',fs-2)

plot3(1:12,-MLD_in,max(100*(T_in-T_out)./T_out),'--','linewidth',1,'color',[0.4 0.4 0.4])
c.FontSize=fs-2;

set(gca,'XTickLabelRotation',45)


% % % % AOU

axes(ha1(4))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{'1000','800','600','400','200','0'})
ylabel('\itz \rm[m]')
set(gca,'XTickLabelRotation',45)

s=surface(1:12,-Z,100*(AOU_in-AOU_out)./AOU_out);
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,[-100 100]);
colormap(ha1(4),deltacols)
c=colorbar;
% c.Ticks=-20:10:20;
c.Location='eastoutside';
title('\Delta AOU [%]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
% text(0,1.4,'(l)','units','normalized','fontsize',fs-2)

plot3(1:12,-zPROD_in,max(100*(AOU_in-AOU_out)./AOU_out),'-.','linewidth',1,'color',[0.4 0.4 0.4])
% plot3(1:12,-zPROD_in-200,max(AOU_in-AOU_out),'-','linewidth',0.2,'color',[0.4 0.4 0.4])

c.FontSize=fs-2;


% % % % %POC l
axes(ha1(3))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{'1000','800','600','400','200','0'})
ylabel('\itz \rm[m]')
set(gca,'XTickLabelRotation',45)

s=surface(1:12,-Z,100*(POCl_in-POCl_out)./POCl_out);
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,[-100 100]);
colormap(ha1(3),deltacols)
c=colorbar;
% c.Ticks=-15:7.5:15;
c.Location='eastoutside';
title('\Delta POC_{\itl\rm} [%]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
% text(0,1.4,'(i)','units','normalized','fontsize',fs-2)

c.FontSize=fs-2;

plot3(1:12,-zPROD_in,max(100*(POCl_in-POCl_out)./POCl_out),'-.','linewidth',1,'color',[0.4 0.4 0.4])
% plot3(1:12,-zPROD_in-200,max(POCl_in-POCl_out),'-','linewidth',0.2,'color',[0.4 0.4 0.4])



%% print
if print_flag==1
    %     exportgraphics(hf1,['Figures/OSM24_out_MeanProfs_POC.pdf'],'ContentType','vector');
    %     exportgraphics(hf2,['Figures/OSM24_out_MeanProfs_Temp.pdf'],'ContentType','vector');
    %     exportgraphics(hf3,['Figures/OSM24_out_MeanProfs_DOXY.pdf'],'ContentType','vector');
    %     exportgraphics(hf4,['Figures/OSM24_out_MeanProfs_Nitrate.pdf'],'ContentType','vector');

    figure(hf1);
    print('Figures/V4/2_LB22_MonthlyProfiles_T_POCs','-dpdf','-r800')
    figure(hf2);
    print('Figures/V4/3_LB22_TimeSections','-dpdf','-r800')

    figure(hf4);
    print('Figures/LB22_TimeSections_PDif','-dpdf','-r800')

    % figure(hf3);
    % print('Figures/LB22_TimeSections_ttest','-dpdf','-r800')
end
