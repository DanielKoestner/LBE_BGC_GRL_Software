%% Get floats, plot locations, calculate poc

close all
clear all
curdir=cd;

curdir=cd;

addpath(genpath([curdir '/Scripts']))
addpath(genpath([curdir '/Data']))
addpath(genpath([curdir '/OneArgo']))
addpath(genpath([curdir '/aux']))

print_flag=1;
save_flag=0;
scale_flag='A'; % A for 1.5r and B for 2r scaling of LBEZ


fs=16;
lw=1.5;
set(0, 'DefaultAxesFontName', 'Times');
set(0, 'DefaultTextFontName', 'Times');


alp=0.65;
alin = 0.9;
warning off
%%
hf1=figure();
set(hf1,'Units','inches','Position', [5 5 7.5 9], 'PaperPosition', [0 0 7.5 9], 'PaperSize', [7.5 9]);
ha1=iSubplot(1,1, 'Gap', [0 0], 'Min', [0.04 0.2], 'Max', [1 1], 'XTickL', 'All', 'YTickL', 'All');


    % LATLIMS=[60 80];
    LATLIMS=[63 76];
    LONLIMS=[-20 25];
m_proj('lambert','long',LONLIMS,'lat',LATLIMS);
% m_coast('patch',[1 .85 .7]);
%     m_gshhs_l('patch',[.5 .6 .5]);

axes(ha1(1)); hold on

% [CS,CH]=m_elev('contourf',[-4000:100:0],'edgecolor','none');
[CS,CH]=m_elev('contourf',[-4000 -2000],'edgecolor','k','ShowText','on');
% m_elev('contour',[-4000 -2000],'edgecolor','k','ShowText','on');
if print_flag==1
    m_gshhs_h('patch',[.85 0.95 .85]);
else
    m_gshhs_l('patch',[.85 0.95 .85]);
end
m_grid('box','fancy','fontsize',fs,'linestyle','none');
cols=m_colmap('blue',400);
cols=cols(200:end,:);
cols2=crameri('lajolla',100);
cols2(1,:)=cols(1,:);
cols2(2,:)=cols(1,:);
cols2(3,:)=cols(1,:);
cols2(4,:)=cols(1,:);
cols2(5,:)=cols(1,:);
COL=[cols(120,:);cols(150:170,:)];
colormap([cols;cols2])
caxis([-4000 2000]);
% colormap([cols(110:180,:)])

% title(['Lofoten Basin, 2010 – 2022'])
text(0.35,0.99,'\bfLofoten Basin, 2010 – 2022','units','normalized','fontsize',fs-2)

% [ax,h]=m_contfbar([.27 .72],0.9,CS,CH,'endpiece','no','axfrac',.02,'edgecolor','none');
ax=colorbar('north');

% title(ax,'\it f');
text(0.89,0.815,'\itn','units','normalized','fontsize',fs-2)
ax.FontSize=fs-5;
ax.Position=[0.85 0.8 0.12 0.015];
ax.Limits=[206 2000];
ax.Ticks=[210 1010 2000];
ax.TickLabels={'1','800','1800'};
% ax.XTickLabel={'4000','2000','0'};

% cb=colorbar;
% % cb.Ticks=[0:.2:1];
% % cb.TickLabels={'','30','60','90','120','150'};
% cb.Location='southoutside';
% cb.FontSize=fs-2;
% cb.Position=[0.33 0.13 0.35 0.025];

% text(0.7,-0.10,'\itz\rm [m]','units','normalized','fontsize',fs-1);

% bndry_lon=[-7 -7 15 15 -7];
% bndry_lat=[67.5 72.5 72.5 67.5 67.5];
% m_line(bndry_lon,bndry_lat,'linewi',1.5,'color','k');     % Area outline ...
% 

m_line(-7:15,67.5*ones(length(-7:15),1)','linewi',1,'color','k');     % Area outline ...
m_line(-7:15,72.5*ones(length(-7:15),1)','linewi',1,'color','k');     % Area outline ...
m_line(-7*ones(length(67.5:0.1:72.5),1)',67.5:0.1:72.5,'linewi',1,'color','k');     % Area outline ...
m_line(15*ones(length(67.5:0.1:72.5),1)',67.5:0.1:72.5,'linewi',1,'color','k');     % Area outline ...
 
%% add floats!
% valid as of 22 August 2023, used OneArgo_Lofoten code first


% load(['LBE_BGC_POC_2010_2022_07-Jun-2024.mat']);
load(['LBE_BGC_POC_2010_2022_' scale_flag '_06-Feb-2025'])


%% markers


mk{1}='s';
mk{2}='^';
mk{3}='d';
mk{4}='v';
mk{5}='p';
mk{6}='h';
mk{7}='>';
mk{8}='<';
mk=[mk mk mk mk];

dnums=[];
cnt=1;
while cnt<=length(lb)
    for i = 1:length(lb{cnt}.dnum)
        dnums=[dnums; lb{cnt}.dnum(i)];
    end
    cnt=cnt+1;
end
 %% plot eddy
EDDY=[];
  %setup a grid for query in/out
  lat_box=[68:0.05:71.5];
  lon_box=[-1:0.05:9];
  [X,Y]=meshgrid(lon_box,lat_box);

  load('LBE_locations.mat');

  for i = 1:length(dnums)
      % find closest eddy location date
      [~,ind]=min(abs(dnums(i)-LBE_dnum));
      ind=ind;
      eddy_lats=LBE_lats(:,ind);

      % if no eddy data, use next date
      if eddy_lats(1)==0
          ind=ind+1;
      end

      eddy_lats=LBE_lats(:,ind);
      eddy_lons=LBE_lons(:,ind);


      % polygon approach
      pgon{i}=polyshape(eddy_lons,eddy_lats);
      pgon2{i}=scale(pgon{i},1.5,LBE_center(:,ind)');
      center(i,:)=LBE_center(:,ind)';

      clear ind

      %check box indices which are inside eddy
      for ii = 1:length(lat_box)
          e(ii,:)=isinterior(pgon2{i},X(ii,:),Y(ii,:));

      end
      if i == 1
      EDDY=e;
      else
          EDDY=EDDY+e;
      end
      clear e
  end



  % 
  % ax2=gca;
  % ax2.Visible='off';
  % axes(ax2);
  EDDY=EDDY+10;
  m_image(lon_box,lat_box,EDDY)

%% Plot floats
f=length(lb);
colors=crameri('acton',f+4);
colors=colors(2:end-3,:);
% colors=flip(colors);
% colors=gray(f+2);


mns_in=[];
mns_out=[];
yrs_in=[];
yrs_out=[];
cnt=1;
thrs=0.702;
dnums2=[];
while cnt<=length(lb)

    floatcol=colors(cnt,:);

    for i = 1:length(lb{cnt}.dnum)
        if            sum(isnan(lb{cnt}.zbin.poc_s(:,i)))<(thrs)*length(lb{cnt}.zbin.z) & lb{cnt}.zprod(i)<900
            dnums2=[dnums2; lb{cnt}.dnum(i)];
            [eddy_id]=in_eddy(lb{cnt}.dnum(i),lb{cnt}.lat(i),lb{cnt}.lon(i),'LBE_locations.mat',sf);
            if eddy_id==1 %& month(datetime(lb{cnt}.dnum(i),'ConvertFrom','datenum'))==8

                m_plot(lb{cnt}.lon(i),lb{cnt}.lat(i),'lines','none','marker','o','markerfacecolor',floatcol,'markeredgecolor','none','markersize',2,'linewidth',0.25);
                mns_in=[mns_in; month(datetime(lb{cnt}.dnum(i),'ConvertFrom','datenum'))];
                yrs_in=[yrs_in; year(datetime(lb{cnt}.dnum(i),'ConvertFrom','datenum'))];

            else
                m_plot(lb{cnt}.lon(i),lb{cnt}.lat(i),'lines','none','marker','o','markerfacecolor',floatcol,'markeredgecolor','none','markersize',2,'linewidth',0.25);
                mns_out=[mns_out; month(datetime(lb{cnt}.dnum(i),'ConvertFrom','datenum'))];
                yrs_out=[yrs_out; year(datetime(lb{cnt}.dnum(i),'ConvertFrom','datenum'))];
            end
        else
        end
    end



    cnt=cnt+1;
end


cnt=1;
while cnt<=length(lb)

    floatcol=colors(cnt,:);
    m_plot(lb{cnt}.lon(end),lb{cnt}.lat(end),'marker','o','markerfacecolor',floatcol,'markersize',6,'markeredgecolor',[0 0 0]);

    cnt=cnt+1;
end

  % m_grid('linewi',2,'tickdir','out');
  
%% Legend
  for i = 1:f
    floatcol=colors(i,:);

    lgd_text{i}=sprintf('%1.0f',lb{i}.wmo);

    m(i)=scatter(2000,2000,50,'o','markerfacecolor',floatcol,'markeredgecolor',floatcol,'markeredgecolor',[0 0 0]);
end


[hl,lines]=legendflex([m],lgd_text, 'ref', gca, ...
    'anchor', {'n', 'n'}, 'buffer', [35 -10], 'xscale', .5, 'ncol',6,'nrow',4, ...
    'box', 'off', 'FontSize', fs-5);


%% Monthly stats

for i = 1:12
inds=ismember(mns_in,i);
yrs_in_cnts(i,:)=histcounts(yrs_in(inds),2009.5:2022.5);
clear inds

inds=ismember(mns_out,i);
yrs_out_cnts(i,:)=histcounts(yrs_out(inds),2009.5:2022.5);
clear inds
end

hf2=figure();
set(hf2,'Units','inches','Position', [5 5 7.5 2.75], 'PaperPosition', [0 0 7.5 2.75], 'PaperSize', [7.5 2.75]);
ha2=iSubplot(2,1, 'Gap', [0 0.09], 'Min', [0.05 0.06], 'Max', [0.98 0.9], 'XTickL', 'All', 'YTickL', 'All');
fs=13;

cols=crameri('lapaz',13);
% [map,num] = tab20(13);
% cols=map;

% cols=gray(13);
% inds=randperm(13);
inds=[7	4	9	13	11	8	6	12	1	5	2	10	3];
% cols=cols(inds,:);

axes(ha2(1))
set(gca,'ticklength',[0.01 0.01])
box on
ha2(1).YGrid='on';
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.5 12.5])
ylim([0 30])
% set(gca,'ytick',[0:10:40])
set(gca,'xtick',[1:12],'xticklabel',{''});
ylabel('\itN')
title('inside LBEZ')
% h=histogram(mns_in,0:13,'facecolor',[0.7 0.7 0.7]);
plot([3.5 3.5],[0 1000],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([6.5 6.5],[0 1000],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([9.5 9.5],[0 1000],'--','color',[0.4 0.4 0.4],'linewidth',0.9)

b1=bar(yrs_in_cnts,'stacked');


axes(ha2(2))
ylabel('\itN')
title('outside LBEZ')
set(gca,'ticklength',[0.01 0.01])
box on
ha2(2).YGrid='on';
set(gca,'yminortick','on')
set(gca,'fontsize',fs)
xlim([0.5 12.5])
ylim([0 150])
% set(gca,'ytick',[0:60:240])
set(gca,'xtick',[1:12],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'XTickLabelRotation',0)
% h=histogram(mns_out,0:13,'facecolor',[0.3 0.3 0.3]);
plot([3.5 3.5],[0 1000],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([6.5 6.5],[0 1000],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([9.5 9.5],[0 1000],'--','color',[0.4 0.4 0.4],'linewidth',0.9)

b2=bar(yrs_out_cnts,'stacked');
colormap(ha2(2),cols);
c2=colorbar;
c2.Position(4)=0.8;
for i = 1:13
b1(i).FaceColor=cols(i,:);
b2(i).FaceColor=cols(i,:);
end

ha2(1).Position(3)=0.82;
ha2(2).Position(3)=0.82;
c2.Ticks=[1/26:1/13:1];
c2.TickLabels=string(2010:2022);
c2.TickLength=0;

%%

hf3=figure();
set(hf3,'Units','inches','Position', [5 5 3 3], 'PaperPosition', [0 0 3 3], 'PaperSize', [3 3]);
ha3=iSubplot(1,1, 'Gap', [0 0], 'Min', [0.04 0.2], 'Max', [1 1], 'XTickL', 'All', 'YTickL', 'All');

cols=m_colmap('blue',400);
cols=cols(200:end,:);

m_proj('orthographic','lat',45,'long',0);
% m_proj('Mollweide','lat',[-90 90],'lon',[-90 90],'rot',45);

    %date of observation
ice_date='20230310';

%Near-Real-Time DMSP SSMIS Daily Polar Gridded Sea Ice Concentrations,
%Version 2 - Microwave observations 25 km

ice_lat=ncread('NSIDC0771_LatLon_PS_N25km_v1.0.nc','latitude');
ice_lon=ncread('NSIDC0771_LatLon_PS_N25km_v1.0.nc','longitude');

ice=ncread('NSIDC0081_SEAICE_PS_N25km_20230101_v2.0.nc','F18_ICECON');
ice=ice*100;
% ice(ice==255)=105;   % Put ocean at top of indices
ice(ice>100)=100;

m_pcolor(ice_lon,ice_lat,ice);
cols3=(crameri('oslo',1000));
cols3=cols3(600:end,:);
cols3(1,:)=cols(50,:);
colormap(cols3);

m_grid('backcolor',cols(50,:),'xtick',[-180:45:180],'xticklabels',[''],'ytick',[-90:30:90],'yticklabels',['']);

% m_gshhs('lc','patch',[.85 0.95 .85],'color','k')
m_coast('patch',[.85 0.95 .85]);

bndry_lon=[-7 -7:15 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0 -1 -2 -3 -4 -5 -6 -7];
bndry_lat=[67.5 72.5*ones(1,length(-7:15)) 67.5*ones(1,length(-7:15))];


% m_line(bndry_lon,bndry_lat,'linewi',1.5,'color','k');     % Area outline ...
m_line(-7:15,67.5*ones(length(-7:15),1)','linewi',1.5,'color','k');     % Area outline ...
m_line(-7:15,72.5*ones(length(-7:15),1)','linewi',1.5,'color','k');     % Area outline ...
m_line(-7*ones(length(67.5:0.1:72.5),1)',67.5:0.1:72.5,'linewi',1.5,'color','k');     % Area outline ...
m_line(15*ones(length(67.5:0.1:72.5),1)',67.5:0.1:72.5,'linewi',1.5,'color','k');     % Area outline ...
 
m_hatch(bndry_lon,bndry_lat,'cross','color',[0.5 0.5 0.5]);


%% Print
if print_flag==1
    figure(hf1);
    print(['Figures/V8/LB_FloatMap_2022_' scale_flag],'-dpdf','-r800')
    figure(hf2);
    print(['Figures/V8/LB_MonthlyYrsStats_2022_' scale_flag],'-dpdf','-r800')
%     exportgraphics(hf1,['Figures/OSM24_LB_Map.pdf'],'ContentType','vector');

    figure(hf3);
    print(['Figures/V8/LB_GlobalMap_2022_' scale_flag],'-dpdf','-r800')
end
