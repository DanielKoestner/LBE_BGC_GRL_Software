%% Get floats, make graphic

close all
clear all
curdir=cd;

curdir=cd;

addpath(genpath([curdir '/Scripts']))
addpath(genpath([curdir '/Data']))
addpath(genpath([curdir '/OneArgo']))
addpath(genpath([curdir '/aux']))

print_flag=0;
save_flag=0;
fs=16;
lw=1.5;
set(0, 'DefaultAxesFontName', 'Times');
set(0, 'DefaultTextFontName', 'Times');

alp=0.65;
alin = 0.9;

weights=[ones(1,5) 2*ones(1,5)]; %weights for linear fit
warning('off')
tic
%% Load data
load('LBE_BGC_POC_2010_2022_08-May-2024.mat');
load('LB_GlobColour_PP_2010_2024_15-Apr-2024.mat');
load('LBE_locations.mat');

%%
dates=734289:5:738561; % 31 May 2010 to Feb 2022
matrix=zeros(length(dates),length(lb));

% find indices of each float which correspond to date window, store in
% matrix with rows for date and columns for float
for i = 1:length(dates)-1
    tmpdate=[dates(i) dates(i+1)];
    for ii=1:length(lb)
        tmpfld=lb{ii}.dnum;
        ind=find(tmpfld>tmpdate(1) & tmpfld<tmpdate(2),1,'first');
        if ~isempty(ind)
        matrix(i,ii)=(ind);
        end
    end
end

% for i = 1:length(dates)
%     flts=find(matrix(i,:)>0);
%     for ii = 1:length(flts)
%         lats{i,flts(ii)}=lb{flts(ii)}.lat(1:matrix(i,flts(ii)));
%         lons{i,flts(ii)}=lb{flts(ii)}.lon(1:matrix(i,flts(ii)));
%         TEs(i,flts(ii))=lb{flts(ii)}.TE(matrix(i,flts(ii)));
%     end
%     clear flts
% end

%try only storing last location
for i = 1:length(dates)
    flts=find(matrix(i,:)>0);
    for ii = 1:length(flts)
        lats(i,flts(ii))=lb{flts(ii)}.lat(matrix(i,flts(ii)));
        lons(i,flts(ii))=lb{flts(ii)}.lon(matrix(i,flts(ii)));
        MES(i,flts(ii))=lb{flts(ii)}.ipoc_mes(matrix(i,flts(ii)));
        MES2(i,flts(ii))=lb{flts(ii)}.ipoc_mes2(matrix(i,flts(ii)));
        
    end
    clear flts
end

lats(lats==0)=NaN;
lons(lons==0)=NaN;

MES(MES==Inf)=NaN;
MES(MES==0)=NaN;
MES2(MES2==Inf)=NaN;
MES2(MES2==0)=NaN;

TEb=MES2./MES;

% TEs(TEs>1)=1; % assume higher than 1 is due to mixed layer and/or horizontal advection

colors=crameri('lapaz',330);
colors=colors(11:310,:);
colors=flip(colors);


for i = 3:length(dates)-4 %628 starts 2019
    time=datetime(dates(i),'convertfrom','datenum');
    hf1=figure();
    set(hf1,'Units','inches','Position', [5 5 7 7.5], 'PaperPosition', [0 0 7 7.5], 'PaperSize', [7 7.5]);
    ha1=iSubplot(1,1, 'Gap', [0 0], 'Min', [0.06 0], 'Max', [1 1], 'XTickL', 'All', 'YTickL', 'All');

    LATLIMS=[60 80];
    LONLIMS=[-20 25];

    m_proj('lambert','long',LONLIMS,'lat',LATLIMS);
    %     [CS,CH]=m_elev('contourf',[-4000:50:0],'edgecolor','none');


    % % % % Primary Production
    ind_pp=find(Mn==month(time)&Yr==year(time));
    m_pcolor(Plg,Plt,log10(PP{ind_pp})'); shading flat
    caxis([log10(100) log10(5000)]);

    m_gshhs_l('patch',[.9 0.95 .9]);
    m_grid('box','fancy','fontsize',fs,'linestyle','none');
    %     colormap(m_colmap('blue',256));
    %     colormap(gray(256)); %b/w bathy
    %     colormap(crameri('lapaz',256))

    % % % % Eddy
    % find closest eddy location date
    [~,inde]=min(abs(dates(i)-LBE_dnum));
    eddy_lats=LBE_lats(:,inde);

    % if no eddy data, use next date
    if eddy_lats(1)==0
        inde=inde+1;
    end


    pgon=polyshape(LBE_lons(:,inde),LBE_lats(:,inde));
    pgon2=scale(pgon,1.50,LBE_center(:,inde)');
    eddy_coords=pgon.Vertices;
    eddy_coords(end+1,:)=eddy_coords(1,:);
    eddy2_coords=pgon2.Vertices;
    eddy2_coords(end+1,:)=eddy2_coords(1,:);
    m_plot(LBE_center(1,inde),LBE_center(2,inde),'lines','none','color',[0.2 0.2 0.2],'marker','.'); 
    m_plot(eddy_coords(:,1),eddy_coords(:,2),'lines','-','color',[0 0 0],'marker','none','linewidth',2); 
    m_plot(eddy2_coords(:,1),eddy2_coords(:,2),'lines','--','color',[0.2 0.2 0.2],'marker','none','linewidth',1);
    clear inde eddy_coords eddy2_coords

    colormap(flip(m_colmap('green',256))); % chla

    ax2=axes;
    ax2.Visible='off';
    colormap(ax2,colors);
    c=colorbar(ax2,'southoutside');
    c.Ticks=[ 0.2 0.8 ];
    c.TickLabels={'Low ~TE','High ~TE'};
    c.FontSize=12;
    c.TickLength=0;
c.Position(2)=0.05;
    axes(ha1(1)); hold on

    title(['Nordic Seas, ' string(time)],'FontSize',fs+2)

    for fn=1:length(lb)
        if i<=10
            lat=lats(1:i,fn);
            lon=lons(1:i,fn);
            mes=MES(1:i,fn);
            mes2=MES2(1:i,fn);
            t=dates(1:i);
            te=TEb(1:i,fn);

            p=glmfit(t,mes,'normal','weights',weights(1:i));
            p2=glmfit(t,mes2,'normal','weights',weights(1:i));

            % p=polyfit(t,mes,1);
            % p2=polyfit(t,mes2,1);
           
            if sum(~isnan(mes))>3
                TE(i,fn)=p2(2)./p(2);
                Ez(i,fn)=p(2);
                Ez500(i,fn)=p2(2);
            else
                TE(i,fn)=NaN;
                Ez(i,fn)=NaN;
                Ez500(i,fn)=NaN;
            end


% for "real" TE
            % if ~isnan(TE(i,fn)) 
            %     tcol=ceil(TE(i,fn)*300);
            %     tcol(tcol>300)=300;
            %     tcol(tcol<1)=1;
            % 
            %     floatcol=colors(tcol,:);
            %     clear tcol
            % else
            %     floatcol='none';
            % end
% For quasi TE
            if ~isnan(te(end))
                tcol=ceil(te(end)*300);
                tcol(tcol>300)=300;
                tcol(tcol<1)=1;
                floatcol=colors(tcol,:);
                clear tcol
            else
                floatcol='none';
            end

          

            tmk=(mes/1000)/2;
            tmk=tmk+5;
            tmk(tmk>11)=11;
            if isnan(tmk(end))
                mk=nanmean(tmk(end-2:end));
            else
                mk=tmk(end);
            end
            if isnan(mk)
                mk=8;
            end


            % plot tail
            lat=lat(~isnan(lat));
            lon=lon(~isnan(lon));
            if ~isempty(lat)
                m_plot(lon(1:end),lat(1:end),'lines','-','color',[0.5 0.5 0.6],'marker','none');
                m_plot(lon(end),lat(end),'lines','-','marker','o','markerfacecolor',floatcol,'markeredgecolor',[0.5 0.5 0.6],'markersize',10);
            end
        elseif i>10
            lat=lats(i-9:i,fn);
            lon=lons(i-9:i,fn);
            mes=MES(i-9:i,fn);
            mes2=MES2(i-9:i,fn);
            t=dates(i-9:i);
            p=glmfit(t,mes,'normal','weights',weights);
            p2=glmfit(t,mes2,'normal','weights',weights);

            %             p=polyfit(t(~isnan(mes)),mes(~isnan(mes)),1);
            %             p2=polyfit(t(~isnan(mes)),mes2~(isnan(mes)),1);
            if sum(~isnan(mes))>4
                TE(i,fn)=p2(2)./p(2);
                Ez(i,fn)=p(2);
                Ez500(i,fn)=p2(2);
            else
                TE(i,fn)=NaN;
                Ez(i,fn)=NaN;
                Ez500(i,fn)=NaN;
            end

            % % for ~Real TE
            % if ~isnan(TE(i,fn))
            %     tcol=ceil(TE(i,fn)*300);
            %     tcol(tcol>300)=300;
            %     tcol(tcol<1)=1;
            % 
            %     floatcol=colors(tcol,:);
            %     clear tcol
            % else
            %     floatcol='none';
            % end

% For quasi TE
            if ~isnan(TEb(i,fn))
                tcol=ceil(TEb(i,fn)*300);
            else
                tcol=ceil(nanmean(TEb(i-3:i+3,fn))*300);  
            end

            if ~isnan(tcol)
                tcol(tcol>300)=300;
                tcol(tcol<1)=1;
                floatcol=colors(tcol,:);
                
            else
                floatcol='none';
            end
            clear tcol

                tmk=(mes/1000)/2;
                tmk=tmk+5;
                tmk(tmk>11)=11;
                if isnan(tmk(end))
                    mk=nanmean(tmk(end-5:end));
                else
                    mk=tmk(end);
                end
                if isnan(mk)
                    mk=8;
                end
               
            % plot tail
            lat=lat(~isnan(lon));
            lon=lon(~isnan(lon));
            if ~isempty(lat)
                m_plot(lon(1:end),lat(1:end),'lines','-','color',[0.5 0.5 0.6],'marker','none');
                m_plot(lon(end),lat(end),'lines','-','marker','o','markerfacecolor',floatcol,'markeredgecolor',[0.5 0.5 0.6],'markersize',mk);
            end
        clear tmk mk
        
        end
    end
print(['Figures/animation/NorBGCArgoAnimation_' num2str(i,'%04.f') '.jpeg'],'-djpeg','-r200');
close all
end



   toc
    %     print('Figures/OSM24_LB_Map.jpg','-jpg','-r300')
    %     exportgraphics(hf1,['Figures/OSM24_LB_Map.pdf'],'ContentType','vector');
%% make video
addpath(genpath('Figures'))
v = VideoWriter('NorBGCArgoAnimation_PP_2022','MPEG-4');
v.Quality=100;
v.FrameRate=16;

cd([curdir '/Figures/animation'])
d=dir('*.jpeg');
cd(curdir);

open(v)
for i = 1:length(d)
A=imread(d(i).name);
writeVideo(v,A);
end

close(v)



%% 
hf1=figure();
set(hf1,'Units','inches','Position', [5 5 7 7.5], 'PaperPosition', [0 0 7 7.5], 'PaperSize', [7 7.5]);
ha1=iSubplot(1,1, 'Gap', [0 0], 'Min', [0.06 0], 'Max', [0.8 1], 'XTickL', 'All', 'YTickL', 'All');


% % % % DEPTH
% LATLIMS=[60 80];
% LONLIMS=[-20 25];
% 
% m_proj('lambert','long',LONLIMS,'lat',LATLIMS);
% [CS,CH]=m_elev('contourf',[-4000:50:0],'edgecolor','none');
% m_gshhs_l('patch',[.9 0.95 .9]);
% m_grid('box','fancy','fontsize',fs,'linestyle','none');
% %     colormap(m_colmap('blue',256));
% colormap(gray(256));
% %     colormap(crameri('lapaz',256))
% 
%   h=colorbar;
%   set(get(h,'ylabel'),'String','seafloor depth [m]');
%   h.Location='eastoutside';
%     set(h,'tickdir','out');
% h.Ticks=[-3000 -2000 -1000 0];
% h.TickLabels={'3000','2000','1000','0'};
% 
% print(['Figures/NorBGCArgoAnimation_DepthLegend.jpeg'],'-djpeg','-r600')
% 


% % % PP
LATLIMS=[60 80];
LONLIMS=[-20 25];

m_proj('lambert','long',LONLIMS,'lat',LATLIMS);

m_pcolor(Plg,Plt,log10(PP{ind_pp})'); shading flat
caxis([log10(100) log10(5000)]);

m_gshhs_l('patch',[.9 0.95 .9]);
m_grid('box','fancy','fontsize',fs,'linestyle','none');
colormap(flip(m_colmap('green',256))); % chla


h=colorbar;
set(get(h,'ylabel'),'String','PP [mg C m^{-2} d^{-1}]');
h.Location='eastoutside';
set(h,'tickdir','out');
% h.Ticks=[log10(100) log10(1000) log10(5000)];
% h.TickLabels={'100','1000','5000'};

h.Ticks=[log10(100) log10(5000)];
h.TickLabels={'100','5000'};

ms=m_plot(0,70,'lines','none','marker','o','markerfacecolor','none','markeredgecolor',[0.5 0.5 0.6],'markersize',6);
ml=m_plot(10,70,'lines','none','marker','o','markerfacecolor','none','markeredgecolor',[0.5 0.5 0.6],'markersize',11);

[hl,lines]=legendflex([ms,ml],{'Low iPOC','High iPOC'}, 'ref', gca, ...
    'anchor', {'n', 'n'}, 'buffer', [0 -30], 'xscale', 1, 'ncol',1,'nrow',2, ...
    'box', 'off', 'FontSize', fs-1);
print(['Figures/NorBGCArgoAnimation_PPLegend.jpeg'],'-djpeg','-r600')