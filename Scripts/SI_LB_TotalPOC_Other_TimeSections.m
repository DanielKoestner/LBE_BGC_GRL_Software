% SI figure
% maybe density time section (in/out/dif)
% total POC time section (in/out/dif)

% Salinity and total POC
% Add another plot for Nitrate and composition, in surface and during
% productive months only (Apr–Sep)
%% set up


close all
clear all

curdir=cd;

addpath(genpath([curdir '/Scripts']))
addpath(genpath([curdir '/Data']))

addpath(genpath([curdir '/aux']))

print_flag=1;
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
load(['LBE_BGC_POC_2010_2022_MonthlyProfiles_' scale_flag '_06-Feb-2025'])

rho_in=gsw_rho_t_exact(S_in,T_in,Z');
rho_out=gsw_rho_t_exact(S_out,T_out,Z');

rholims=[1027 1033];
slims=[34.9 35.2];
olims=[260 320];
aolims=[-20 40];
nlims=[0 12];
zlims=[-1000 0];
complims=[20 2000];
logpoc=1;

if logpoc==1
    poclims=[1 250];
    pocllims=[1 50];
else
    poclims=[0 150];
    pocllims=[0 45];
end

COMP_in=CHLA_in./BBP_in;
COMP_out=CHLA_out./BBP_out;

% fix comp stuff
% for i=1:12
%     tmp=COMP_in(:,i)==inf;
%     tmp2=COMP_out(:,i)==inf;
%     
%     COMP_in(tmp,i)=min(COMP_in(Z<200,i));
%     COMP_out(tmp2,i)=min(COMP_out(Z<200,i));
% end
%% Section plots, 3 x 1 plots for in, out, and difference

% aoucols=crameri('bam',240);
% aoucols=(aoucols(61:end,:));
% aoucols=crameri('bam',28);
% aoucols=(aoucols(8:end,:));

aoucols=crameri('bam',20);
aoucols=(aoucols(6:end,:));

deltacols=crameri('vik',12);
cols=crameri('batlow',12);
colst=crameri('lapaz',12);

sat=-0.1; %brightness scale
aoucols=brighten(aoucols,sat);
deltacols=brighten(deltacols,sat);
cols=brighten(cols,0.08);
fs=12;

hf2=figure();
set(hf2,'Units','inches','Position', [5 5 8 11], 'PaperPosition', [0 0 8 11], 'PaperSize', [8 11]);
ha1=iSubplot(4,3, 'Gap', [0 0.01], 'Min', [0.04 0.01], 'Max', [0.98 1], 'XTickL', 'All', 'YTickL', 'All');

%% POC

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

s=surface(1:12,-Z,log10(POC_in+POCl_in));
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,log10(poclims));
colormap(ha1(1),cols)
c=colorbar;
c.Location='northoutside';
text(0.55,1.4,'POC inside [mg m^{-3}]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(a)','units','normalized','fontsize',fs-2)

c.Ticks=log10([1 10 100]);
c.TickLabels={'1','10','100'};
plot3(1:12,-EZD_in,max(log10(POC_in+POCl_in)),':','linewidth',1,'color',[0.7 0.7 0.7])
plot3(1:12,-MLD_in,max(log10(POC_in+POCl_in)),'--','linewidth',1,'color',[0.7 0.7 0.7])



c.FontSize=fs-2;
% text(1,-zPROD_in(1),'PZD','units','normalized','fontsize',fs-2,'color',[0.5 0.5 0.5])

axes(ha1(2))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{''})

s=surface(1:12,-Z,log10(POC_out+POCl_out));
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,log10(poclims));
colormap(ha1(2),cols)
c=colorbar;
c.Location='northoutside';
c.Ticks=log10([1 10 100]);
c.TickLabels={'1','10','100'};
text(0.55,1.4,'POC outside [mg m^{-3}]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(b)','units','normalized','fontsize',fs-2)

set(gca,'XTickLabelRotation',45)
c.FontSize=fs-2;

plot3(1:12,-EZD_out,max(log10(POC_out+POCl_out)),':','linewidth',1,'color',[0.7 0.7 0.7])
plot3(1:12,-MLD_out,max(log10(POC_out+POCl_out)),'--','linewidth',1,'color',[0.7 0.7 0.7])

axes(ha1(3))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-1000 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-1000:200:0],'yticklabel',{''})

s=surface(1:12,-Z,(POC_in+POCl_in)-(POC_out+POCl_out));
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
% caxis(gca,[-25 25]);
caxis(gca,[-24 24]);
colormap(gca,deltacols)
c=colorbar;
c.Ticks=-16:8:16;
c.Location='northoutside';
text(0.55,1.4,'\Delta POC [mg m^{-3}]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(c)','units','normalized','fontsize',fs-2)

set(gca,'XTickLabelRotation',45)
c.FontSize=fs-2;

plot3(1:12,-zPROD_in,max((POC_in+POCl_in)-(POC_out+POCl_out)),'-.','linewidth',1,'color',[0.4 0.4 0.4])


%% POC zoom

axes(ha1(4))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-200 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-200:40:0],'yticklabel',{'200','160','120','80','40','0'})
ylabel('\itz \rm[m]')
set(gca,'XTickLabelRotation',45)

s=surface(1:12,-Z,log10(POC_in+POCl_in));
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,log10(poclims));
colormap(gca,cols)
c=colorbar;
c.Location='northoutside';
text(0.55,1.4,'POC inside [mg m^{-3}]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(d)','units','normalized','fontsize',fs-2)

c.Ticks=log10([1 10 100]);
c.TickLabels={'1','10','100'};
plot3(1:12,-EZD_in,max(log10(POC_in+POCl_in)),':','linewidth',1,'color',[0.7 0.7 0.7])
plot3(1:12,-MLD_in,max(log10(POC_in+POCl_in)),'--','linewidth',1,'color',[0.7 0.7 0.7])



c.FontSize=fs-2;
% text(1,-zPROD_in(1),'PZD','units','normalized','fontsize',fs-2,'color',[0.5 0.5 0.5])

axes(ha1(5))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-200 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-200:40:0],'yticklabel',{''})

s=surface(1:12,-Z,log10(POC_out+POCl_out));
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,log10(poclims));
colormap(gca,cols)
c=colorbar;
c.Location='northoutside';
c.Ticks=log10([1 10 100]);
c.TickLabels={'1','10','100'};
text(0.55,1.4,'POC outside [mg m^{-3}]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(e)','units','normalized','fontsize',fs-2)

set(gca,'XTickLabelRotation',45)
c.FontSize=fs-2;

plot3(1:12,-EZD_out,max(log10(POC_out+POCl_out)),':','linewidth',1,'color',[0.7 0.7 0.7])
plot3(1:12,-MLD_out,max(log10(POC_out+POCl_out)),'--','linewidth',1,'color',[0.7 0.7 0.7])

axes(ha1(6))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-200 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-200:40:0],'yticklabel',{''})

s=surface(1:12,-Z,(POC_in+POCl_in)-(POC_out+POCl_out));
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
% caxis(gca,[-25 25]);
caxis(gca,[-24 24]);
colormap(gca,deltacols)
c=colorbar;
c.Ticks=-16:8:16;
c.Location='northoutside';
text(0.55,1.4,'\Delta POC [mg m^{-3}]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(f)','units','normalized','fontsize',fs-2)

set(gca,'XTickLabelRotation',45)
c.FontSize=fs-2;

plot3(1:12,-zPROD_in,max((POC_in+POCl_in)-(POC_out+POCl_out)),'-.','linewidth',1,'color',[0.4 0.4 0.4])



%% COMP
axes(ha1(7))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-200 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-200:40:0],'yticklabel',{'200','160','120','80','40','0'})
ylabel('\itz \rm[m]')

s=surface(1:12,-Z,log10(COMP_in));
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,log10(complims));
colormap(gca,colst)
c=colorbar;
c.Ticks=log10([20 100 1000]);
c.TickLabels={'20' '100' '1000'};
c.Location='northoutside';
text(0.55,1.4,'\varsigma inside [mg m^{-2}]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(g)','units','normalized','fontsize',fs-2)

plot3(1:12,-MLD_in,ones(size(zPROD_in))*4,'--','linewidth',1,'color',[0.4 0.4 0.4])
plot3(1:12,-EZD_in,ones(size(zPROD_in))*4,':','linewidth',1,'color',[0.4 0.4 0.4])

set(gca,'XTickLabelRotation',45)
c.FontSize=fs-2;
text(1,-MLD_in(1),'MLD','units','normalized','fontsize',fs-2,'color',[0.4 0.4 0.4])

axes(ha1(8))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-200 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-200:40:0],'yticklabel',{''})


s=surface(1:12,-Z,log10(COMP_out));
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,log10(complims));
colormap(gca,colst)
c=colorbar;
c.Ticks=log10([20 100 1000]);
c.TickLabels={'20' '100' '1000'};
c.Location='northoutside';
text(0.55,1.4,'\varsigma outside [mg m^{-2}]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(h)','units','normalized','fontsize',fs-2)

plot3(1:12,-MLD_out,ones(size(zPROD_in))*4,'--','linewidth',1,'color',[0.4 0.4 0.4])
plot3(1:12,-EZD_out,ones(size(zPROD_in))*4,':','linewidth',1,'color',[0.4 0.4 0.4])

set(gca,'XTickLabelRotation',45)
c.FontSize=fs-2;

axes(ha1(9))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-200 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-200:40:0],'yticklabel',{''})

s=surface(1:12,-Z,COMP_in-COMP_out);
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,[-420 420]);
colormap(gca,deltacols)
c=colorbar;
c.Ticks=-280:140:280;
c.Location='northoutside';
text(0.55,1.4,'\Delta \varsigma [mg m^{-2}]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(i)','units','normalized','fontsize',fs-2)

plot3(1:12,-zPROD_in,ones(size(zPROD_in))*500,'-.','linewidth',1,'color',[0.4 0.4 0.4])
c.FontSize=fs-2;

set(gca,'XTickLabelRotation',45)


%% Nitrate

axes(ha1(10))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-200 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-200:40:0],'yticklabel',{'200','160','120','80','40','0'})
ylabel('\itz \rm[m]')

s=surface(1:12,-Z,N_in);
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,nlims);
colormap(gca,colst)
c=colorbar;
c.Ticks=[0:2:12];
c.Location='northoutside';
text(0.55,1.4,'NO_3^- [\mumol kg^{-1}]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(j)','units','normalized','fontsize',fs-2)

plot3(1:12,-MLD_in,ones(size(zPROD_in))*20,'--','linewidth',1,'color',[0.7 0.7 0.7])
plot3(1:12,-EZD_in,ones(size(zPROD_in))*20,':','linewidth',1,'color',[0.7 0.7 0.7])

set(gca,'XTickLabelRotation',45)
c.FontSize=fs-2;
text(1,-MLD_in(1),'MLD','units','normalized','fontsize',fs-2,'color',[0.4 0.4 0.4])

axes(ha1(11))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-200 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-200:40:0],'yticklabel',{''})


s=surface(1:12,-Z,N_out);
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,nlims);
colormap(gca,colst)
c=colorbar;
c.Ticks=[0:2:12];
c.Location='northoutside';
text(0.55,1.4,'NO_3^- [\mumol kg^{-1}]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(k)','units','normalized','fontsize',fs-2)

plot3(1:12,-MLD_out,ones(size(zPROD_in))*20,'--','linewidth',1,'color',[0.7 0.7 0.7])
plot3(1:12,-EZD_out,ones(size(zPROD_in))*20,':','linewidth',1,'color',[0.7 0.7 0.7])

set(gca,'XTickLabelRotation',45)
c.FontSize=fs-2;

axes(ha1(12))
set(gca,'ticklength',[0.02 0.02])
box on
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.6 12])
ylim([-200 0])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'ytick',[-200:40:0],'yticklabel',{''})

s=surface(1:12,-Z,N_in-N_out);
s.MeshStyle='column';
s.EdgeAlpha=0.15;
s.EdgeColor=[0.4 0.4 0.4];
s.FaceColor='interp';
caxis(gca,[-6 6]);
colormap(gca,deltacols)
c=colorbar;
c.Ticks=-4:2:4;
c.Location='northoutside';
text(0.55,1.4,'\Delta NO_3^- [\mumol kg^{-1}]','units','normalized','fontsize',fs-0.5,'fontweight','bold','horizontalalignment','center')
text(0,1.4,'(l)','units','normalized','fontsize',fs-2)

plot3(1:12,-zPROD_in,ones(size(zPROD_in))*20,'-.','linewidth',1,'color',[0.4 0.4 0.4])
c.FontSize=fs-2;

set(gca,'XTickLabelRotation',45)



%%

axes(ha1(9));
p1=plot(20:21,20:21,'--','linewidth',1,'color',[0.4 0.4 0.4]);
p2=plot(20:21,20:21,':','linewidth',1,'color',[0.4 0.4 0.4]);
p3=plot(20:21,20:21,'-.','linewidth',1,'color',[0.4 0.4 0.4]);

[hl,lines]=legendflex([p1 p2 p3],{'MLD','EZD','PZD'}, 'ref', gca, ...
'anchor', {'ne', 'ne'}, 'buffer', [5 68], 'xscale', 0.5, 'ncol',3,'nrow',1, ...
'box', 'off', 'FontSize', fs-2);

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


%% print
if print_flag==1

    figure(hf2);
    print(['Figures/V8/SI/LB22_TimeSections_Extra_' scale_flag],'-dpdf','-r800')


end
