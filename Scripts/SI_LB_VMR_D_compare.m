% LBE analysis, must extract only dates when we have floats.


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
fs=16;
lw=1.75;
set(0, 'DefaultAxesFontName', 'Times');
set(0, 'DefaultTextFontName', 'Times');


%% Get data
load('LBE_BGC_POC_2010_2022_22-Jul-2024.mat','in','out')

% note that vmr_100 is median diameter (µm) in upper 100 m for each
% profile, under some assumptions of Qbb and volume.
vmr_in=[];
for i = 1:length(in)
    vmr_in=[vmr_in in{i}.vmr_100];
end

vmr_out=[];
for i = 1:length(out)
    vmr_out=[vmr_out out{i}.vmr_100];
end


[f_in,x_in,flo_in,fup_in]=ecdf(vmr_in,'bounds','on');
[f_out,x_out,flo_out,fup_out]=ecdf(vmr_out,'bounds','on');
%% plot
colors2=crameri('batlow',8);
colsin=colors2(7,:);
colsout=colors2(3,:);
sat=0.4;

hf=figure();
set(hf,'Units','inches','Position', [5 5 5 5], 'PaperPosition', [0 0 5 5], 'PaperSize', [5 5]);
ha1=iSubplot(1,1, 'Gap', [0.02 0.04], 'Min', [0.09 0.09], 'Max', [0.95 0.96], 'XTickL', 'All', 'YTickL', 'All');

axes(ha1(1))
hold on
grid on
xlabel('Median diameter in upper 100 m [µm]')
ylabel('CDF')
set(gca,'yminortick','on')
ylim([0 1]);
xlim([10 1000]);
set(gca,'ytick',0:0.1:1)
set(gca,'xtick',[10 100 1000],'xticklabel',{'10','100','1000'})
set(gca,'fontsize',fs-2)
set(gca,'ticklength',[0.025 0.025])
% set(gca,'xtick',[0:25:100])
set(gca,'xscale','log')

h=shadedErrorBar(x_in,f_in,[f_in-flo_in fup_in-f_in],'patchSaturation',sat,'lineProps',{'linewidth',lw,'color',colsin});
h.edge(1).LineStyle='none';
h.edge(2).LineStyle='none';

h2=shadedErrorBar(x_out,f_out,[f_out-flo_out fup_out-f_out],'patchSaturation',sat,'lineProps',{'linewidth',lw,'color',colsout});
h2.edge(1).LineStyle='none';
h2.edge(2).LineStyle='none';

% text(0.76,0.95,sprintf('\\mu = %1.0f km',mean(r_lb)),'units','normalized')
% text(0.76,0.9,sprintf('\\sigma = %1.0f km',std(r_lb)),'units','normalized')
lgd_txt{1}=sprintf('inside LBEZ, \\itN \\rm= %1.0f',length(f_in));
lgd_txt{2}=sprintf('outside LBEZ, \\itN \\rm= %1.0f',length(f_out));


[hl,lines]=legendflex([h.mainLine h2.mainLine],lgd_txt, 'ref', gca, ...
    'anchor', {'nw', 'nw'}, 'buffer', [5 -5], 'xscale', .5, 'ncol',6,'nrow',4, ...
    'box', 'off', 'FontSize', fs-5);

%% print

if print_flag==1
    print(['Figures/V4/SI/LB22_VMR_Compare_' date],'-dpdf','-r800')
end
