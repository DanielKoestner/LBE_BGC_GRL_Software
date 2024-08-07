%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [Data_bbp_proc,Data_chl_proc] = process_BGCARGO_BBP_CHL(Data_CHL,Data_BBP,path_oneargo,varargin);

% Script to Process BGC-ARGO BBP700 and chla from a .mat file, converting
% raw BBP700 and adjusted Chla to BBP700 and Chla associated with small/large particles

%% Output

% Data_bbp_proc: a structure containing as field float ids. Each float id is itself a substurcture
% that contains 

% i. BBP700_l: BBP700 associated with large particles (m^-1)
% ii BBP700_s: BBP700 associated with small particles (m^-1)
% iii BBP700_despike: despiked BBP700 that is the raw BBP700 < the threshold (thr; see code below for thr value) (m^-1)
% iv BBP700_noise: instrument noise compomenent, also considered "blank" of the spike signal (see Briggs et al., 2020) (m^-1).
% v. BBP700_r: Refractory pool of BBP700, constant for each float (m^-1)

% Data_chl_proc: same as Data_bbp_proc, but for Chla (mg/m3):

% i. chl_l: chla associated with large particles (mg/m3)
% ii chl_s: chla associated with small particles (mg/m3)
% iii CHLA_ADJUSTED_despike: chla at same pressure levels as despiked BBP700 (mg/m3)
% iv chl_noise: instrument noise compomenent, also considered "blank" of the spike signal (see Briggs et al., 2020) (mg/m3).
% v. chl_r: Refractory pool of chla, constant for each float (mg/m3)


%% Input


% i. Data: A matfile that contains two data structures, one for chla and one for BBP700. The data structures should contain the
%    following fields:

%    Data_bbp, Data_chl contain:
%    WMOID1 [1x1 structure]
%    WMOID2 [1x1 structure]
%    ...
%    WMOIDN [1x1 structure]

%    where WMOIDs are the float ids for the floats obtained from the BGC-ARGO repository
%    Each WMOID is itself a structure, and must contain the following field for Data_bbp:

%    BBP700 [NxM] (m^-1)
%    CHLA_ADJUSTED [NxM] (mg/m^3)
%    PRES_ADJSUTED [NxM] (dbar)
%    Time (serial date number, in days since Jan 0, 0000)

%    where N is the number of pressure levels, and M is the number of profiles

% ii. Data_chl: A string that is the name of the Chla data structure in "Data"

% iii. Data_bbp: A string that is the name of the BBP data structure in "Data" 

%% Optional Arguments

% i. path_OneArgo: path to OneArgo package. Default assumes that OneArgo is a directory in the level up the "scripts" directory

% ii. path_matfiles: path to Data matfile. Default assumes that matfiles is a directory in the level up the "scripts" directory

% iii. bckgnoise: how the background noise is processed. 'float' computes a single value for the entire float's lifetime, 'annual' computes anually varying background noise, 'monthy' computes monthly varying background noise,
%, and 'ann_winter/spring/summer/fall' computes anually varying background noise, but computed only from data for the winter, spring, summer, or fall of that year. 'ann_winter/spring/summer/fall' follows the protocol
% of Biff et al., 2023, Journal of Marine Systems


% The methodology to process chl/BBP data follows Briggs et al., 2020, Science 
% Paul Lerner, Sep 24, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% background noise: default 'float' is to compute one value for the each float. Can also compute noise monthly or yearly (yearly based on wintertime data)

bckgnoise = 'float';

% default names of chl, bbp700, and pres arrays

var_bbp = 'BBP700';
var_chl = 'CHLA_ADJUSTED';
var_pres = 'PRES';
win = 11;
% Parse optional arguments, 20230902 DK Returned error for length(varargin)
% (comment out for now)

for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'bckg_noise')
        bckgnoise = varargin{i+1}
    elseif strcmpi(varargin{i}, 'var_bbp')
        var_bbp = varargin{i+1};
    elseif strcmpi(varargin{i}, 'var_chl')
        var_chl = varargin{i+1};
    elseif strcmpi(varargin{i}, 'var_pres')
        var_pres = varargin{i+1};
    elseif strcmpi(varargin{i}, 'window')
        win = varargin{i+1};
    end
end


addpath(path_oneargo);

% load dataset: discontinued for now, but might add back in as opt-in at a later point

%load([path_matfiles,'/',Data,'.mat']);

% have global setting variables to easily share among functions

global proc_settings

% percentile of bbp_sr (small + "refractory" particles -> Briggs et al., 2020) between 850-900 m below which bbp is considerered to all be "refractory"
% For Briggs et al, this was at 25%
pctile_sr = 25; % percentage of all "small+refractory" bbp 
pctile_chl = 25;  % percentage of all "small+refractory" chl 

% threshold above which to remove bbp data
% For Brigs et al., the threshold is 0.008 

thr = 0.01 ; % in m^-1

thr_z = 25; % in dbar, this is the pressure range around which all data is removed if BBP700 greater than "thr".


% Pressure level below which BBP_resid is bin-averaged, concatenated (between profiles), and used to find instrument noise
% Briggs et al., 2020 use 300 m

noise_z = 300; 


% size of bins used for bin-averaging (Briggs et al., 2020 use 50 m).

noise_dz = 50; 

% upper and lower bounds wihtin which refractory compoment of small+refractory particle-associated bbp700 is computed (dbar)
% Brigss et al., 2020 use 850 and 900 m. In the presentation we used all depths below 850 (using 3000 dbar here since floats never go
% that deep);

refr_z1 = 850;
refr_z2  = 900;

% months for computing noise using only seasonal data
if strcmp(bckgnoise,'ann_winter') == 1;
  months = [12 1 2]; % winter months DEC, JAN, FEB
  seas_name = 'winter'
elseif strcmp(bckgnoise,'ann_spring') == 1;
  months = [3 4 5]; % spring months MAR, APR, MAY
  seas_name = 'spring'
elseif strcmp(bckgnoise,'ann_summer') == 1;
  months = [6 7 8]; % summer months JUN, JUL, AUG
  seas_name = 'summer'
elseif strcmp(bckgnoise,'ann_fall') == 1;
  months = [9 10 11]; % fall months SEP, OCT, NOV
  seas_name = 'fall'
else
  months = ['NaN'];
end     



% set global settings for fixed values

proc_settings.pctile_sr = pctile_sr; % should be input
proc_settings.pctile_chl = pctile_chl; % should be input
proc_settings.thr = thr; % should be input
proc_settings.thr_z = thr_z; % should be input
proc_settings.noise_z = noise_z;
proc_settings.noise_dz = noise_dz;
proc_settings.refer_z1 = refr_z1;
proc_settings.refer_z2 = refr_z2;
proc_settings.months = months;

% format for conversion from serial to datetime

formatOut_ymd = 'yyyy-MM-dd'; % for annual/ann_season noise computation
formatOut_ym = 'yyyyMM'; % for monthly noise computation


if contains(bckgnoise,'ann_');
  bckgnoise = 'ann_season';
end

% loop over floats

fieldnms = fieldnames(Data_BBP); % these are the float ids

for i=1:length(fieldnms);
  fieldnm = fieldnms(i); 
  fprintf(['\n Starting Float',char(fieldnm),'...\n']);

% Extract the float structure  for float i from the datasets
 
  Float_Struc_BBP = eval(['Data_BBP.',char(fieldnm)]);
  Float_Struc_Chl = eval(['Data_CHL.',char(fieldnm)]);

% assign variable "PRES" to chl data structure
  Float_Struc_Chl.PRES = Float_Struc_BBP.PRES;
% "despike" refers to BB700 < thr 
  Float_Struc_BBP.BBP700_despike = NaN.*ones(size(Float_Struc_BBP.(var_bbp)));
  Float_Struc_Chl.chl_despike = NaN.*ones(size(Float_Struc_Chl.(var_chl)));

% arrays containing background noise for variables

  Float_Struc_BBP.BBP700_noise = NaN.*ones(size(Float_Struc_BBP.(var_bbp)));
  Float_Struc_Chl.chl_noise = NaN.*ones(size(Float_Struc_Chl.(var_chl)));

% check if BBP and chl arrays are populated: if not, through warning

  BBP_check = Float_Struc_BBP.(var_bbp);
  chl_check = Float_Struc_Chl.(var_chl);
  
  warn_BBP = 0;
  warn_chl = 0;
  if sum(~isnan(BBP_check(:))) == 0; 
     warning('The float contains no BBP700 data. Please double-check float arrays')
     warn_BBP = 1;
  end

  if sum(~isnan(chl_check(:))) == 0;                                 
     warning('The float contains no chl-a data. Please double-check float arrays')
     warn_chl = 1;
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% First obtain the background noise  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% extract some arrays common to all noise processsing methods
  timearr = datetime(Float_Struc_BBP.TIME,'ConvertFrom','datenum','Format',formatOut_ymd); % Time in datetime (year, month, day) format
  timevec = timearr(1,:); % time along pressure dimension is uniform, only need time per cycle

  year_un = unique(year(timevec)); % find all unique years


%% if using annual noise based on the winter of that year

  switch bckgnoise
    case 'ann_season';
      fprintf(['processing noise anually using ',seas_name,' data\n'])
      [BBP700_noise_asea,chl_noise_asea] = extract_aseanoise(timevec,year_un,Float_Struc_BBP.BBP700,Float_Struc_Chl.(var_chl),Float_Struc_Chl.(var_pres));
      Float_Struc_BBP.BBP700_noise = BBP700_noise_asea;
      Float_Struc_Chl.chl_noise = chl_noise_asea;


%% if using anually varying background noise, considers dataset over entire year
    case 'annual'

      fprintf('processing noise anually using dataset over entire year\n')
      [BBP700_noise_ann,chl_noise_ann] = extract_annoise(timevec,year_un,Float_Struc_BBP.(var_bbp),Float_Struc_Chl.(var_chl),Float_Struc_Chl.(var_pres));
      Float_Struc_BBP.BBP700_noise = BBP700_noise_ann;
      Float_Struc_Chl.chl_noise = chl_noise_ann;

%% if using a monthly varying background noise, considers dataset for each month

    case 'monthly'
      fprintf('processing noise monthly \n')
      timearr_mon = datetime(Float_Struc_BBP.TIME,'ConvertFrom','datenum','Format',formatOut_ym); % Time in datetime (yearmonth) format
      timevec_mon = timearr_mon(1,:); % time along pressure dimension is uniform, only need time per cycle
      timevec_mon = year(timevec_mon).*100 + month(timevec_mon); % convert time-date to double
      yearmonth_un = unique(timevec_mon); % extract unique yearmonths
      [BBP700_noise_mon,chl_noise_mon] = extract_monnoise(timevec_mon,yearmonth_un,Float_Struc_BBP.(var_bbp),Float_Struc_Chl.(var_chl),Float_Struc_Chl.(var_pres));
      Float_Struc_BBP.BBP700_noise = BBP700_noise_mon;
      Float_Struc_Chl.chl_noise = chl_noise_mon;
%% if using a constant backgroiund noise for the entire float lifetime

    otherwise
      fprintf('processing noise for entire float \n')

% take variables for entire float and pass them to function to compute background noise

      BBP700 = Float_Struc_BBP.(var_bbp);
      chl = Float_Struc_Chl.(var_chl);
      Pres = Float_Struc_Chl.(var_pres);

      [BBP700_noise,chl_noise] = calc_bckg_noise(Pres,BBP700,chl,win);

      Float_Struc_BBP.BBP700_noise(:,:) = BBP700_noise;
      Float_Struc_Chl.chl_noise(:,:) = chl_noise;

    end % end over which type of backgroudn noise to compute   
 
% now we filter profiles again and add the instrument noisei

  BBP_sr_tot = 0;  % an array of BBP for each float, used to find the 5th percentile of BBP_sr between refr_z2 dbar;
  chl_sr_tot = 0;  % an array of chl for each float, used to find the 5th percentile of chl_sr between refr_z1 and dbar;

% preallocated arrays in structures

  Float_Struc_BBP.BBP700_l = NaN.*ones(size(Float_Struc_BBP.(var_bbp)));
  Float_Struc_BBP.BBP700_s = NaN.*ones(size(Float_Struc_BBP.(var_bbp)));
  Float_Struc_Chl.chl_l = NaN.*ones(size(Float_Struc_Chl.(var_chl)));
  Float_Struc_Chl.chl_s = NaN.*ones(size(Float_Struc_Chl.(var_chl)));


%% Second for loop is used here after the background noise for the float is computed, to be applied to the filtered BBP_tot profiles %%%%%%%%%%
  for j=1:length(Float_Struc_Chl.(var_pres)(1,:));
    BBP_tot_tmp = Float_Struc_BBP.(var_bbp)(:,j);
    Pres_tmp = Float_Struc_Chl.(var_pres)(:,j);
    chl_tot_tmp = Float_Struc_Chl.(var_chl)(:,j);

% determine window size for moving minimum and maximum filters


% zero out negative chl values

    chl_tot_tmp(find(chl_tot_tmp<0)) = 0;

% remove within 25 m of BBP_tot > thr 
    
    ind_hi = find(BBP_tot_tmp>thr);
    for k = 1:length(BBP_tot_tmp);
      if any(Pres_tmp(k) >Pres_tmp(ind_hi) -thr_z) && any(Pres_tmp(k) <Pres_tmp(ind_hi) +thr_z);
         BBP_tot_tmp(k) = NaN;
         Pres_tmp(k) = NaN;
         chl_tot_tmp(k) = NaN;
      end
    end

% find indices where both BBP700 and chl are real

    indnonan = find(~isnan(BBP_tot_tmp) & ~isnan(chl_tot_tmp));

% remove NaN values
    Pres_tmp = Pres_tmp(indnonan);
    BBP_tot_tmp = BBP_tot_tmp(indnonan);
    chl_tot_tmp = chl_tot_tmp(indnonan);


% extract noise from arrays, and remove values overlapping with NaNs above

    BBP700_noise_tmp = Float_Struc_BBP.BBP700_noise(:,j);
    chl_noise_tmp = Float_Struc_Chl.chl_noise(:,j);
    BBP700_noise_tmp = BBP700_noise_tmp(indnonan);
    chl_noise_tmp = chl_noise_tmp(indnonan);


% create arrays of BBP/chl after the initial despiking. These are stricly to compare with profiles after processing, to ensure only a small loss of data
% Note that despiking sets BBP of values in vcinity of thr to NaN, AS WELL AS Chl and Pressure

    Float_Struc_BBP.BBP700_despike(indnonan,j) = BBP_tot_tmp;
    Float_Struc_Chl.chl_despike(indnonan,j) = chl_tot_tmp;

% break out of loop if all NaNs
    if sum(~isnan(Pres_tmp)) == 0;
       continue
    end

% same min/max filtering as in prevous loop

    BBP_fil_tmp = movmin(BBP_tot_tmp,win); % 11 pt moving minimum filter
    BBP_fil_tmp = movmax(BBP_fil_tmp,win); % 11 pt moving minimum filter

    chl_fil_tmp = movmin(chl_tot_tmp,win); % 11 pt moving minimum filter
    chl_fil_tmp = movmax(chl_fil_tmp,win); % 11 pt moving minimum filter

    BBP_sr = BBP_fil_tmp+BBP700_noise_tmp; % refractory + small particles
    BBP_l =  BBP_tot_tmp - BBP_sr; % large particels
    BBP_l(find(BBP_l<10^-5)) = 0; % zero out negative BBP_l values. This may be possible in the case that the BBP700 signal is smaller than the instrument nosie of the float

    chl_sr = chl_fil_tmp+chl_noise_tmp; % refractory + small particles
    chl_l =  chl_tot_tmp - chl_sr; % large particles
    chl_l(find(chl_l<10^-5)) = 0; % as with BBP_l values, NaN out negative chl_l values

% concatenate chl_rs and BBP_sr values between 850-900m: to be used to compute refractory copmonent
    chl_sr_850to900 = chl_sr(find(Pres_tmp>refr_z1 & Pres_tmp<refr_z2));
    chl_sr_tot = vertcat(chl_sr_tot,chl_sr_850to900);

    BBP_sr_850to900 = BBP_sr(find(Pres_tmp>refr_z1 & Pres_tmp<refr_z2));
    BBP_sr_tot = vertcat(BBP_sr_tot,BBP_sr_850to900);

% make sure to insert processed data a pressure levels corresponding with real chla/bbp values
    Float_Struc_BBP.BBP700_l(indnonan,j) = BBP_l;
    Float_Struc_BBP.BBP700_s(indnonan,j) = BBP_sr; 
    Float_Struc_Chl.chl_l(indnonan,j) = chl_l;
    Float_Struc_Chl.chl_s(indnonan,j) = chl_sr;
  end

  BBP_sr_tot = BBP_sr_tot(2:end);
  chl_sr_tot = chl_sr_tot(2:end);


% calcualte pctile_sr of BBP_sr between 850-900 m. Use last float value if no deep data exist
  if sum(~isnan(BBP_sr_tot)) ~= 0;
    P_bbp = prctile(BBP_sr_tot,pctile_sr); 
    P_chl = prctile(chl_sr_tot,pctile_chl);
  else % use minimum values if not enought deep data 
    P_bbp = min(Float_Struc_BBP.BBP700_s(find(Float_Struc_BBP.BBP700_s>0)));
    P_chl = min(Float_Struc_Chl.chl_s(find(Float_Struc_Chl.chl_s>0)));
  end

% if there is no float data for chl or BBP700 array, set refractory value to NaN
  if warn_BBP == 1 || warn_chl == 1;
    P_bbp = NaN
    P_chl = NaN
  end
% apply correction to obtain BBP_s;
  Float_Struc_BBP.BBP700_s = Float_Struc_BBP.BBP700_s - P_bbp;  
  Float_Struc_BBP.BBP700_s(find(Float_Struc_BBP.BBP700_s<10^-5)) = 0; 
  Float_Struc_BBP.BBP700_r = P_bbp; % refractory material is a single value for the entire float: like instrument noise
  Float_Struc_Chl.chl_s = Float_Struc_Chl.chl_s - P_chl;
  Float_Struc_Chl.chl_s(find(Float_Struc_Chl.chl_s<10^-5)) = 0;
  Float_Struc_Chl.chl_r = P_chl; % refractory material is a single value for the entire float: like instrument noise
  Float_Struc_BBP.(var_pres) = Float_Struc_Chl.(var_pres);

  eval(['Data_bbp_proc.',char(fieldnm),'= Float_Struc_BBP;']);    
  eval(['Data_chl_proc.',char(fieldnm),'= Float_Struc_Chl;']);

end

%% apply QC flags to BBP in small and large particles (for BBP only)
%% To be consistent with the rest of the processing, chla values are also
%% set to NaN when the qc tests for BBP fail
%fprintf('QCing large and small bbp700 data\n');
%[Data_BBP_RTQC_s] = RTQC_BGCARGO_BBP(Data_bbp_proc,path_oneargo,'var_bbp','BBP700_s','var_pres',var_pres,'do_noisy_bin',1);
%[Data_BBP_RTQC_l] = RTQC_BGCARGO_BBP(Data_bbp_proc,path_oneargo,'var_bbp','BBP700_l','var_pres',var_pres,'B_RES_THRESHOLD',0.0002,'B_FRAC_OUTLIER',0.05,'do_noisy_bin',1);

%for i=1:length(fieldnms);
%  Float_BBP = Data_BBP_RTQC_s.(fieldnms{i});
%  Float_CHL = Data_chl_proc.(fieldnms{i});
%  BBP700_tmp = Float_BBP.BBP700_s;
%  CHL_tmp = Float_CHL.chl_s;
%  BBP700_flag_tmp = Float_BBP.BBP700_s_RTQC_Flags;
%  for iflag = [0,3,4,9];
%    BBP700_tmp(find(BBP700_flag_tmp==iflag)) = NaN;
%    CHL_tmp(find(BBP700_flag_tmp==iflag)) = NaN;
%  end

%  BBP700_s = BBP700_tmp;
%  CHL_s = CHL_tmp;
%  Data_bbp_proc.(fieldnms{i}).BBP700_s_RTQC = BBP700_s; % replace BBP700_s with BBP700_s nann'd where iflag = 3,4, or 9
%  Data_chl_proc.(fieldnms{i}).chl_s_RTQC = CHL_s;
%  Data_bbp_proc.(fieldnms{i}).BBP700_s_RTQC_Flags = BBP700_flag_tmp;
%  Data_bbp_proc.(fieldnms{i}).BBP700_s_RTQC_failed_tests = Data_BBP_RTQC_s.(fieldnms{i}).BBP700_s_RTQC_failed_tests;

%  Float_BBP = Data_BBP_RTQC_l.(fieldnms{i});
%  Float_CHL = Data_chl_proc.(fieldnms{i});
%  BBP700_tmp = Float_BBP.BBP700_l;
%  Float_BBP2 = Data_bbp_proc.(fieldnms{i});
%  BBP700_tmp2 = Float_BBP2.BBP700_l;
%  CHL_tmp = Float_CHL.chl_l;
%  BBP700_flag_tmp = Float_BBP.BBP700_l_RTQC_Flags;
%  for iflag = [0,3,4,9];
%    BBP700_tmp(find(BBP700_flag_tmp==iflag)) = NaN;
%    CHL_tmp(find(BBP700_flag_tmp==iflag)) = NaN;
%  end
      
%  BBP700_l = BBP700_tmp;
%  CHL_l = CHL_tmp;
%  Data_bbp_proc.(fieldnms{i}).BBP700_l_RTQC = BBP700_l; % replace BBP700_l with BBP700_l nann'd where iflag = 3,4, or 9
%  Data_bbp_proc.(fieldnms{i}).BBP700_l_RTQC_failed_tests = Data_BBP_RTQC_l.(fieldnms{i}).BBP700_l_RTQC_failed_tests;
%  Data_chl_proc.(fieldnms{i}).chl_l_RTQC = CHL_l;
%  Data_bbp_proc.(fieldnms{i}).BBP700_l_RTQC_Flags = BBP700_flag_tmp;

%end
end
