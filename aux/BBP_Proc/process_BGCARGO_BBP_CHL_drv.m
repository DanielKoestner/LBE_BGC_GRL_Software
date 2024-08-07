%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A driver script for processing BBP and CHL small
% and large particle data. Includes steps to apply
% Real-time quality control to BBP, and then
% process small and large BBP and CHL for data
% wehre QC flags are < 3 or 8.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc


do_plot = 0;
prof=108;
do_noisy_esd = 1;
do_median_out = 0;
do_noisy_bin = 0
flag_samplesize = 0;
alpha_sf = 0.01;
alpha_subsf = 0.05;
alpha_meso = 0.15;
%load BGC-ARGo dataset

path_to_data = '/discover/nobackup/plerner/GO-BGC-ARGO/BGC_ARGO_toolbox/Project_Lofoten/matfiles/';
filename = 'BBP_CHL_BGC-ARGO.mat'
%fileout_ls = 'BBP_CHL_BGC-ARGO_processed_float_DallOlmonoise.mat';
%fileout_RTQC = 'BBP_CHL_BGC-ARGO_RTQC_DallOlmonoise.mat';
fileout_ls = 'BBP_CHL_BGC-ARGO_processed_float_DallOlmonoise_ESD_ameso015.mat';
fileout_RTQC = 'BBP_CHL_BGC-ARGO_RTQC_DallOlmonoise_ESD_ameso015.mat';
load([path_to_data,filename]);
%fieldselect = {'F6903590'};
fieldselect = {'F6900799','F6902547','F6903553','F6903554','F6903568','F6903577','F6903578','F6903590'}
fnames = fieldnames(Data_BBP);
for i = 1:length(fnames);
  if any(strcmp(fnames{i},fieldselect)) == 1
    mask(i) = 1;
  else
    mask(i) = 0;
  end
end

keep = logical(mask);
D_bbp = struct2cell(Data_BBP);
D_chl = struct2cell(Data_qc_Chl);
Data_BBP =  cell2struct(D_bbp(keep),fnames(keep));
Data_qc_Chl =   cell2struct(D_chl(keep),fnames(keep));
% path to OneArgo package
path_oneargo = '/discover/nobackup/plerner/GO-BGC-ARGO/BGC_ARGO_toolbox/Project_Lofoten/OneArgo/';

% real time quality control data
[Data_BBP_RTQC] = RTQC_BGCARGO_BBP(Data_BBP,path_oneargo,'B_RES_THRESHOLD',0.00025,'B_FRAC_OUTLIER',0.05,'do_noisy_bin',do_noisy_bin,'do_noisy_esd',do_noisy_esd,'do_median_out',do_median_out,'flag_samplesize',flag_samplesize,...
'alpha_sf',alpha_sf,'alpha_subsf',alpha_subsf,'alpha_meso',alpha_meso);
save([path_to_data,fileout_RTQC],'Data_BBP_RTQC');

% Apply qc flags: qc for 3,4, and 9 are Nanned.

fnames = fieldnames(Data_BBP_RTQC);

for i =1:length(fnames);%24
  Float = Data_BBP_RTQC.(fnames{i});
  PRES = Float.PRES;
  BBP700_tmp = Float.BBP700;
  BBP700_flag_tmp = Float.BBP700_QC;
  for iflag = [0,3,4,5,9];
    BBP700_tmp(find(BBP700_flag_tmp==iflag)) = NaN;
  end
  BBP700_RTQC = BBP700_tmp;
  Data_BBP_RTQC.(fnames{i}).BBP700_RTQC = BBP700_RTQC; % naming BBP700 Nanned where qc flags are 3,4,9 "BBP700_RTQC". 
  if do_plot == 1;
    BBP700_tmp = Float.BBP700;
    hf5 = figure();
    ylim([-2000 -5]);

    ylabel('\itp \rm[dbar]')
    set(gca,'ytick',[-2000 -1000 -100 -5],'yticklabel',{'','1000','100','5'})
    set(gca,'ytick',[-2000 -1000 -100 -5])
    set(gca,'yticklabel',{''})

    xlabel('bbp(700) [m^{-1}]')

    s=plot(BBP700_tmp(:,prof),PRES(:,prof).*-1,'ok','markerfacecolor','k');hold on;
    s=plot(BBP700_RTQC(:,prof),PRES(:,prof).*-1,'or','markerfacecolor','r');
    title([str2num(fnames{i}(2:8)) '; ' datestr(Float.TIME(1,prof),'yyyymmdd')])
    set(gca,'yscale','linear')
    set(gca,'TickLength',[0.04 0.04])
        
    bbplims=([0.00001 0.03]);
    set(gca,'xscale','log')
    set(gca,'xtick',[0.0001 0.001 0.01])
  end   
end



% processing of chl and BBP into large and small particles

[Data_bbp_proc,Data_chl_proc] = process_BGCARGO_BBP_CHL(Data_qc_Chl,Data_BBP_RTQC,path_oneargo,'var_bbp','BBP700_RTQC');

save([path_to_data,fileout_ls],'Data_bbp_proc','Data_chl_proc');
