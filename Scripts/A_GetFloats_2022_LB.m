% Get floats for Lofoten basin Project, 2022–2023 
% updated 2024.03.28
% Do not seperate floats here, keep floats which have been inside the
% domain of interest at some point (approx eddy region)

% clean up workspace
clear; close all; clc; 
p=genpath('aux');
addpath(genpath('Scripts'));
addpath(p); clear p
curdir=cd;


%% ============== User Inputs ================================

% lat and lon bounds 
% lat_LB = [68,  71.5];
% lon_LB = [-1, 10];

lat_LB = [67.5,  72.5];
lon_LB = [-7, 15];

% lat_LB = [62,  75];
% lon_LB = [-10, 15];

% date range 
t1 = [2010, 1, 1]; 
t2 = [2022, 2, 9];



initialize_tf = false ; % boolean to decide whether to initialzie or not. If it has been initialized recently, then make this false. 

save_flag=0;
proc_flag=0;
%% ============================================================
% main starts here

% initialize OneArgo toolbox. this may take a few minutes depending on internet speed

initialize_argo(); 

global Sprof Float Settings

%% get float ID and profiles that have BBP in LB general region
fprintf('\nSelecting BBP floats and profiles...\n\n');
% [LB_BBP_flts, LB_BBP_profs] = select_profiles(lon_LB, lat_LB, t1, t2, ...
%     'Sensor', 'BBP700','outside','space'); % keep profiles outside spatial domain

[LB_BBP_flts, LB_BBP_profs] = select_profiles(lon_LB, lat_LB, t1, t2, ...
    'Sensor', 'BBP700','outside','none','direction','a'); %switch to outside none to only keep profiles inside larger box
fprintf('\nFinished selecting LB BBP floats and profiles...\n\n');

fprintf('\nDownloading metadata and running BBP RTQC...\n\n');
download_meta_files(LB_BBP_flts);
[Data_BBP_RTQC] = RTQC_BGCARGO_BBP_nc(LB_BBP_flts,[curdir '/Profiles'],[curdir '/Meta']);
fprintf('\nFinished BBP RTQC...\n\n');

%not enough bbp data for meaningful conclusion, determined manually
bad_ids=[1902602; 1902603; 2903794; 6900876; 6902546; 6903567; 7901028];
inds=~ismember(LB_BBP_flts,bad_ids);
LB_BBP_profs=LB_BBP_profs(inds);
LB_BBP_flts=LB_BBP_flts(inds);


%% make plot of location
% fprintf('\nPlotting all BBP float locations...\n');
show_trajectories(LB_BBP_flts,...
   'color','multiple','float_profs',LB_BBP_profs, ... % this plots different floats in different colors
   'title', 'LBE BBP floats');

%% VARS TO GET
var2get={'TEMP', 'PSAL',...
    'BBP700','CHLA','DOXY','PRES','NITRATE'};
% var2get={'ALL'};

% choose which QFs you want to extract
QF2use = [1 2 5 8]; % good and probably good data. Need to add 8 because for some Provor floats, salinity is interpolated onto DOXY so it gets a 8 flag. 

%% NOT LBE
% NEED "RAW" BBP, PAR, AND CDOM BUT "ADJUSTED" CHLA, DOXY, AND NITRATE, SPLIT INTO TWO FILES
[data_LBbbp, meta_LBbbp] = load_float_data(LB_BBP_flts,{'PRES','BBP700','DOWNWELLING_PAR','CDOM','PSAL','TEMP'},LB_BBP_profs);
data_LBbbpqc = qc_filter(data_LBbbp,{'PRES','BBP700','DOWNWELLING_PAR','CDOM','PSAL','TEMP'}, QF2use, 'raw', 'y');

[data_LB, meta_LB] = load_float_data(LB_BBP_flts,'ALL',LB_BBP_profs);
data_LBqc = qc_filter(data_LB,var2get , QF2use);

data_LBbbpqc = calc_mld(data_LBbbpqc, 'calc_mld_temp', 1, 'calc_mld_dens',1);
%% add step to remove floats with all (or primarily) NaN BBP data


%% save
if save_flag==1
    save(['Data/LBE_BGCDATA_2010_2022_' date '.mat'],'data*','lat_LB','lon_LB','t1','t2')
end

%% Process/partition bbp and chla into small/large
if proc_flag==1
        % partition
    [LB_bbp,LB_chl] = process_BGCARGO_BBP_CHL(data_LBqc,data_LBbbpqc,[curdir '/aux/OneArgo']);
    save(['Data/LBE_BGCDATA_2010_2022_proc_' date '.mat'],'LB*','lon_LB','lat_LB')
end
