%% BGC_Argo_POC_Parameters2
% derive parameters for POC depth profiles from BGC Argo
% inputs are depth, poc, min and max depths, and a euphotic zone depth

% note that depth data will be negative but use positive values for min and
% max depth. Min depth is only used for sigmoid fitting.

% updated in V2 includes integration below euphotic zone to 200 m and
% redefine mesopelagic zone as 200 to 1000

% V3: DK 20230930
% Fixed depth poc, (epi <2000 and meso >200)

% update w/ comments and fewer variables
% also store best fit slope as one structure/vector
% add optional plot feature, requires comp values and some ID, can be cycle
% number, date, or just running number is ok

%input ID can be either a datenumber (Matlab standard from 0 jan 0000) or a
%cycle number in first row and WMO ID in second row

function [euphotic_poc,subeuphotic_poc,epipelagic_poc,mesopelagic_poc,total_poc,zp,slope,stat] = POC_zParameters(depth,poc,euphoticzone,TEzs,plot_flag,id,comp)

if nargin==3
    plot_flag=0;
    id=NaN(2,1); %not needed if not plotting
    comp=NaN(size(poc)); %not needed if not plotting
    TEzs=[];
elseif nargin==6
    comp=ones(size(poc)); %filler
end

% force depth to be positive values
depth=abs(depth);

% can be updated as needed, but min depth only potentialy impacts sigmoid
% fitting and max depth should usually be 1000
mindepth=10;
maxdepth=1000;


% Remove any NaN data for calculations and plotting
depth=depth(~isnan(poc));
comp=comp(~isnan(poc));
poc=poc(~isnan(poc));

% use a deepest depth
comp=comp(depth<=maxdepth);
poc=poc(depth<=maxdepth);
depth=depth(depth<=maxdepth);

% added DK 20240129, if few comp values over 0 exist, fix minimum value to
% be 10. This could be updated in future, but 14 seems to be a more common
% minimum. Not sure it makes a huge difference.
if sum(comp>0)>10
    comp(comp==0)=min(comp(comp>0)); % fix all 0 chla values to minimum comp value
else
    comp(comp==0)=10;
end


% determine cummulative integral of POC up to max depth
for ii = 2:sum(depth>-maxdepth)
    tpoc=poc(1:ii);
    z=depth(1:ii);
    if sum(~isnan(tpoc))>1
        int(ii)=trapz(z,tpoc);
    else
        int(ii)=0;
    end
end

intpoc=trapz(depth(depth<=maxdepth),poc(depth<=maxdepth));
int=int./intpoc; % as fraction of total poc integral
total_poc=intpoc;

% find indices nearest percentiles of interest
[~,ind10]=min(abs(int-0.1));
[~,ind25]=min(abs(int-0.25));
[~,ind50]=min(abs(int-0.5));
[~,ind75]=min(abs(int-0.75));
[~,ind90]=min(abs(int-0.9));


% find max and min POC  (min is lowest value over 0)
[maxp ind1]=max(poc);
[minp ind2]=min(poc(poc>0));


% find sigmoid function fit 
ind=find(depth>=mindepth&depth<maxdepth);
[param,stat]=sigm_fit(depth(ind),poc(ind),[],[minp,maxp,depth(ind50),0.05],0);

tnum=sum((poc(ind)-stat.ypred).^2);
tdem=sum((poc(ind)-nanmean(poc(ind))).^2);

R2=1-((length(stat.ypred)-1)/(length(stat.ypred)-length(param)))*(tnum/tdem);
clear tnum tdem


% poc integrated within various zones
indsub=depth>=euphoticzone;
if ~isnan(euphoticzone) & sum(indsub)>5 & sum(indsub==0)>1
    euphotic_poc=trapz(depth(depth<euphoticzone),poc(depth<euphoticzone));
    subeuphotic_poc=trapz(depth(indsub),poc(indsub));
    for i = 1:length(TEzs)
        indsub2=depth>(euphoticzone+TEzs(i));
        if sum(indsub2)>5
            subeuphotic_poc=[subeuphotic_poc trapz(depth(indsub2),poc(indsub2))];
        else
            subeuphotic_poc=[subeuphotic_poc NaN];
        end
        clear indsub2
    end
else
    euphotic_poc=NaN;
    subeuphotic_poc=NaN(1,length(TEzs)+1);
end



indepi=depth<200;
epipelagic_poc=trapz(depth(indepi),poc(indepi));

indme=depth>=200;
mesopelagic_poc=trapz(depth(indme),poc(indme));
    
% store variables
zp=[depth(ind10) depth(ind25) depth(ind50) depth(ind75) depth(ind90)];

ypred=stat.ypred;
slope=param(4);

pocmax=[maxp depth(ind1)];
pocmin=[minp depth(ind2)];

if plot_flag==1
    hf=figure();
    set(hf,'Units','inches','Position', [5 5 4 6], 'PaperPosition', [0 0 4 6], 'PaperSize', [4 6]);
    ha=iSubplot(1,1, 'Gap', [0 0.01], 'Min', [0.14 0.06], 'Max', [0.96 0.92], 'XTickL', 'All', 'YTickL', 'All');

    ylabel('\itz \rm[m]')
    set(gca,'ytick',[-1000:200:0],'yticklabel',{'1000','800','600','400','200','0'})
    xlabel('POC [mg m^{-3}]')
    grid on
    set(gca,'TickLength',[0.025 0.025])
    set(gca,'fontsize',12)
    set(0, 'DefaultAxesFontName', 'Times');
    set(0, 'DefaultTextFontName', 'Times');
    set(gca,'YMinorTick','on')
    set(gca,'XMinorTick','on')
    colormap(crameri('batlow',256));
    ylim([-1000 0]);
%     xlim([0 200])

    c=colorbar;
    c.Title.String='\varsigma';

    scatter(poc,-z,20,comp,'filled') %plot data
    xl=xlim;
    plot(ypred,-z(ind),'-','color',[0.2 0.2 0.2],'linewidth',1.25) %plot fit
    plot(xl,[-zp(1) -zp(1)],'--','color',[0.8 0.8 0.8]) % plot percentiles
    plot(xl,[-zp(2) -zp(2)],'--','color',[0.7 0.7 0.7]) % plot percentiles
    plot(xl,[-zp(3) -zp(3)],'--','color',[0.6 0.6 0.6]) % plot percentiles
    plot(xl,[-zp(4) -zp(4)],'--','color',[0.5 0.5 0.5]) % plot percentiles
    plot(xl,[-zp(5) -zp(5)],'--','color',[0.4 0.4 0.4]) % plot percentiles

    if id(1) > 1000
        title(string(datetime(id(1),'ConvertFrom','datenum')))
    else
        title(['Cycle ' num2str(id)])
    end

    text(-0.07,1.07,num2str(id(2)),'units','normalized','fontweight','bold')
end
end