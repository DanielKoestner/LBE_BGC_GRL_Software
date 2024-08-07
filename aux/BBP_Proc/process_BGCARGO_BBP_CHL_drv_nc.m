%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A driver script for processing BBP and CHL small
% and large particle data. Includes steps to apply
% Real-time quality control to BBP, and then
% process small and large BBP and CHL for data
% where QC flags are < 3 or 8.

% The difference between this script and 
% process_BGCARGO_BBP_CHL_drv.m" is that it uses as
% input the .nc files in the "Profiles" directory
% that are extracted using the OneArgo utilities,
% (e.g., "load_float_data"), and outputs only
% a revised "BBP700_QC" array with revised QC
% flags based on the Dal'Olmo processing routines.
% whereas the other takes as input data structures
% that the user must extact from the netcdf file in 
% "Profiles" (e.g., using the qc_filter function
% from OneArgo), and output the revised QC flags 
% in the same matlab data structure.

% As an added backup, the function "RTQC_BGCARGO_BBP_nc"
% still returns a data structure (here called "Data_BBP_RTQC")
%, which contains information on which QC tests failed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc


do_plot = 0;
prof=108;
do_noisy_esd = 0;
do_median_out = 0;
do_noisy_bin = 0
flag_samplesize = 0;
alpha_sf = 0.01;
alpha_subsf = 0.05;
alpha_meso = 0.10;
%load BGC-ARGo dataset

path_to_data = '/discover/nobackup/plerner/GO-BGC-ARGO/BGC_ARGO_toolbox/Project_Lofoten/matfiles/';
path_to_ncfiles = '/discover/nobackup/plerner/GO-BGC-ARGO/BGC_ARGO_toolbox/Project_Lofoten/Scripts/Profiles';
path_meta = '/discover/nobackup/plerner/GO-BGC-ARGO/BGC_ARGO_toolbox/Project_Lofoten/Scripts/Meta';
fileout_RTQC = 'BBP_CHL_BGC-ARGO_RTQC_DallOlmonoise.mat';
%fileout_RTQC = 'BBP_CHL_BGC-ARGO_RTQC_DallOlmonoise_ESD.mat';
%fieldselect = {'F6903590'};
floats = {'6900799','6902547','6903553','6903554','6903577','6903578','6903590','6903591'}
%floats = {'6902547'}
path_oneargo = '/discover/nobackup/plerner/GO-BGC-ARGO/BGC_ARGO_toolbox/Project_Lofoten/OneArgo/';

% real time quality control data
%[Data_BBP_RTQC] = RTQC_BGCARGO_BBP_nc(floats,path_to_ncfiles,path_meta,'B_RES_THRESHOLD',0.00025,'B_FRAC_OUTLIER',0.05,'do_noisy_bin',do_noisy_bin,'do_noisy_esd',do_noisy_esd,'do_median_out',do_median_out,'flag_samplesize',...
%flag_samplesize,'alpha_sf',alpha_sf,'alpha_subsf',alpha_subsf,'alpha_meso',alpha_meso);
[Data_BBP_RTQC] = RTQC_BGCARGO_BBP_nc(floats,path_to_ncfiles,path_meta);
save([path_to_data,fileout_RTQC],'Data_BBP_RTQC');

% Apply qc flags: qc for 3,4, and 9 are Nanned.

fnames = fieldnames(Data_BBP_RTQC);

for i =1:length(fnames);%24
  Float = Data_BBP_RTQC.(fnames{i});
  PRES = Float.PRES;
  BBP700_tmp = Float.BBP700;
  BBP700_flag_tmp = Float.BBP700_QC;
  for iflag = [0,3,4,9];
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



