%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [Data_bbp_proc,Data_chl_proc] = process_BGCARGO_BBP_CHL(Data,Data_chl,Data_bbp,varargin);

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

%    BBP700 [NxM]
%    CHLA_ADJUSTED [NxM]
%    PRES_ADJSUTED [NxM]

%    where N is the number of pressure levels, and M is the number of profiles

% ii. Data_chl: A string that is the name of the Chla data structure in "Data"

% iii. Data_bbp: A string that is the name of the BBP data structure in "Data" 

%% Optional Arguments

% i. path_OneArgo: path to OneArgo package. Default assumes that OneArgo is a directory in the level up the "scripts" directory

% ii. path_matfiles: path to Data matfile. Default assumes that matfiles is a directory in the level up the "scripts" directory

% The methodology to process chl/BBP data follows Briggs et al., 2020, Science 
% Paul Lerner, Aug 22, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% defaults of optional arguments: assume OneArgo and matfiles are directories in the level up the "scripts" directory

path_OneArgo = pwd;
path_OneArgo = erase(path_OneArgo,'Scripts');
path_OneArgo = [path_OneArgo,'OneArgo'];

path_matfiles = pwd;
path_matfiles = erase(path_matfiles,'Scripts');
path_matfiles = [path_matfiles,'/Data'];



% Parse optional arguments, 20230902 DK Returned error for length(varargin)
% (comment out for now)

% for i = 1:2:length(varargin)-1
%     if strcmpi(varargin{i}, 'path_OneArgo')
%         path_OneArgo = varargin{i+1};
%     elseif strcmpi(varargin{i}, 'path_matfiles')
%         path_matfiles = varargin{i+1};
%     end
% end

% load dataset

addpath(path_OneArgo); 

load([path_matfiles,'/',Data,'.mat']);


% Assign Data BBP and CHL data structures from function arguments
Data_BBP = eval(Data_bbp); 
Data_CHL = eval(Data_chl);

% percentile of bbp_sr (small + "refractory" particles -> Briggs et al., 2020) between 850-900 m below which bbp is considerered to all be "refractory"
% For Briggs et al, this was at 25%
pctile_sr = 25; % percentage of all "small+refractory" bbp 
pctile_chl = 25;  % percentage of all "small+refractory" chl 

% threshold above which to remove bbp data
% For Brigs et al., the threshold is 0.008!

thr = 0.008 ; % in m^-1 % was 0.02

thr_z = 25; % in dbar, this is the pressure range around which all data is removed if BBP700 greater than "thr".

% window for moving min/max

win = 11;

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

% loop over floats

fieldnms = fieldnames(Data_BBP); % these are the float ids

for i=1:length(fieldnames(Data_BBP)); %updated 20230902 DK, removed "_qc" from name
  fieldnm = fieldnms(i); 
  fprintf(['\n Starting Float',char(fieldnm),'...\n']);

% Extract the float structure  for float i from the datasets
 
  Float_Struc_BBP = eval(['Data_BBP.',char(fieldnm)]);
  Float_Struc_Chl = eval(['Data_CHL.',char(fieldnm)]);

% loop over profiles, first to get instrument noise for each float. Arrays are initialized as "0" here so there is something to concatenate

  BBP_resid_bin_tot = 0; % total vector of all binned BBP_resid profiles, needed to compute instrucment noise
  Pres_BBP_tmp_bin_tot = 0;

  chl_resid_bin_tot = 0; % total vector of all binned chl profiles, needed to compute instrucment noise
  Pres_chl_tmp_bin_tot = 0;

% "despike" refers to BB700 < thr 
  Float_Struc_BBP.BBP700_despike = NaN.*ones(size(Float_Struc_BBP.BBP700));
  Float_Struc_Chl.CHLA_ADJUSTED_despike = NaN.*ones(size(Float_Struc_Chl.CHLA_ADJUSTED));

%%%%%% First for loop over is used to obtained background noise for the entire float (same noise applied to every prifle in the second loop)  %%%%%%%%%%

  for j=1:length(Float_Struc_Chl.PRES_ADJUSTED(1,:));


    
% filter BBP770

    BBP_tot_tmp = Float_Struc_BBP.BBP700(:,j);
    Pres_tmp = Float_Struc_Chl.PRES_ADJUSTED(:,j);
    chl_tot_tmp = Float_Struc_Chl.CHLA_ADJUSTED(:,j);
  
% zero out negative chl values

    chl_tot_tmp(find(chl_tot_tmp<0)) = 0;

% remove within 25 m of BBP_tot > thr 
    ind_hi = find(BBP_tot_tmp>thr);
    for k = 1:length(BBP_tot_tmp);
      if any(Pres_tmp(k) >Pres_tmp(ind_hi) - thr_z) && any(Pres_tmp(k) <Pres_tmp(ind_hi) + thr_z);
         BBP_tot_tmp(k) = NaN;
         chl_tot_tmp(k) = NaN;
         Pres_tmp(k) = NaN;
      end
    end

     
% remove NaN values
    indnonan = find(~isnan(BBP_tot_tmp) & ~isnan(chl_tot_tmp));
    
    Pres_tmp = Pres_tmp(indnonan);
    BBP_tot_tmp = BBP_tot_tmp(indnonan);
    chl_tot_tmp = chl_tot_tmp(indnonan);
   
 
% break out of loop if all NaNs
    if sum(~isnan(Pres_tmp)) == 0;
       continue
    end


% create arrays of BBP/chl after the initial despiking. These are stricly to compare with profiles after processing, to ensure only a small loss of data
% Note that despiking sets BBP of values in vcinity of thr to NaN, AS WELL AS Chl and Pressure

    Float_Struc_BBP.BBP700_despike(indnonan,j) = BBP_tot_tmp;
    Float_Struc_Chl.CHLA_ADJUSTED_despike(indnonan,j) = chl_tot_tmp;

% filter BBP.chl profiles. BBP_fil_tmp should be refractory + small particles
    BBP_fil_tmp = movmin(BBP_tot_tmp,win); % 11 pt moving minimum filter
    BBP_fil_tmp = movmax(BBP_fil_tmp,win); % 11 pt moving minimum filter

    chl_fil_tmp = movmin(chl_tot_tmp,win); % 11 pt moving minimum filter
    chl_fil_tmp = movmax(chl_fil_tmp,win); % 11 pt moving minimum filter
% isolate BBP700 residual

    BBP_resid = BBP_tot_tmp - BBP_fil_tmp; % this is large particels + background of instrument noise, or residual spikes in Briggs et al., 2020
    chl_resid = chl_tot_tmp - chl_fil_tmp; % this is large particels + background of instrument noise, or residual spikes in Briggs et al., 2020

% bin average BBP_resid level noise_z in order to find instrument noise
    if max(Pres_tmp) > noise_z;
      edges = (noise_z:noise_dz:max(Pres_tmp));
      [~,~,loc]=histcounts(Pres_tmp,edges); % number of depths per bin
      Pres_b300 = Pres_tmp;
      BBP_resid_b300 = BBP_resid;
      BBP_resid_b300(find(loc==0)) = []; % remove any bins where there are no pressure values
      Pres_b300(find(loc==0)) = []; % remove any bins where there are no pressure values
      loc(find(loc==0)) = [];% remove any bins where there are no pressure values
   
      Pres_BBP_tmp_bin = accumarray(loc(:),Pres_b300(:),[],@mean);
      BBP_resid_bin = accumarray(loc(:),BBP_resid_b300(:),[],@mean);% depth-binned mean
      BBP_resid_bin_tot = vertcat(BBP_resid_bin_tot,BBP_resid_bin(:));
      Pres_BBP_tmp_bin_tot = vertcat(Pres_BBP_tmp_bin_tot,Pres_BBP_tmp_bin(:));

      [~,~,loc]=histcounts(Pres_tmp,edges); % number of depths per bin
      Pres_b300 = Pres_tmp;
      chl_resid_b300 = chl_resid;
      chl_resid_b300(find(loc==0)) = []; % remove any bins where there are no pressure values
      Pres_b300(find(loc==0)) = []; % remove any bins where there are no pressure values
      loc(find(loc==0)) = [];% remove any bins where there are no pressure values

      Pres_chl_tmp_bin = accumarray(loc(:),Pres_b300(:),[],@mean);
      chl_resid_bin = accumarray(loc(:),chl_resid_b300(:),[],@mean);% depth-binned mean
      chl_resid_bin_tot = vertcat(chl_resid_bin_tot,chl_resid_bin(:));
      Pres_chl_tmp_bin_tot = vertcat(Pres_chl_tmp_bin_tot,Pres_chl_tmp_bin(:));
    end

  end

% calcualte median of each composite vertical profile

  Pres_BBP_tmp_bin_tot = Pres_BBP_tmp_bin_tot(2:end);
  BBP_resid_bin_tot = BBP_resid_bin_tot(2:end);
  chl_resid_bin_tot = chl_resid_bin_tot(2:end);

  median_resid_bbp = median(BBP_resid_bin_tot);
  median_resid_chl = median(chl_resid_bin_tot);

% remove all bins 2x> than the median

  Pres_BBP_tmp_bin_tot = Pres_BBP_tmp_bin_tot(find(BBP_resid_bin_tot<= 2*median_resid_bbp));
  BBP_resid_bin_tot = BBP_resid_bin_tot(find(BBP_resid_bin_tot<= 2*median_resid_bbp));
  Pres_chl_tmp_bin_tot = Pres_chl_tmp_bin_tot(find(chl_resid_bin_tot<= 2*median_resid_chl));
  chl_resid_bin_tot = chl_resid_bin_tot(find(chl_resid_bin_tot<= 2*median_resid_chl));

  Float_Struc_BBP.BBP700_noise = median(BBP_resid_bin_tot); % Instrument noise for this float
  Float_Struc_Chl.chl_noise = median(chl_resid_bin_tot); % Instrument noise for this float

% now we filter profiles again and add the instrument noisei

  BBP_sr_tot = 0;  % an array of BBP for each float, used to find the 5th percentile of BBP_sr between refr_z2 dbar;
  chl_sr_tot = 0;  % an array of chl for each float, used to find the 5th percentile of chl_sr between refr_z1 and dbar;

% preallocated arrays in structures

  Float_Struc_BBP.BBP700_l = NaN.*ones(size(Float_Struc_BBP.BBP700));
  Float_Struc_BBP.BBP700_s = NaN.*ones(size(Float_Struc_BBP.BBP700));
  Float_Struc_Chl.chl_l = NaN.*ones(size(Float_Struc_Chl.CHLA_ADJUSTED));
  Float_Struc_Chl.chl_s = NaN.*ones(size(Float_Struc_Chl.CHLA_ADJUSTED));


%% Second for loop is used here after the background noise for the float is computed, to be applied to the filtered BBP_tot profiles %%%%%%%%%%
  for j=1:length(Float_Struc_Chl.PRES_ADJUSTED(1,:));
    BBP_tot_tmp = Float_Struc_BBP.BBP700(:,j);
    Pres_tmp = Float_Struc_Chl.PRES_ADJUSTED(:,j);
    chl_tot_tmp = Float_Struc_Chl.CHLA_ADJUSTED(:,j);
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

% break out of loop if all NaNs
    if sum(~isnan(Pres_tmp)) == 0;
       continue
    end

% same min/max filtering as in prevous loop

    BBP_fil_tmp = movmin(BBP_tot_tmp,win); % 11 pt moving minimum filter
    BBP_fil_tmp = movmax(BBP_fil_tmp,win); % 11 pt moving minimum filter

    chl_fil_tmp = movmin(chl_tot_tmp,win); % 11 pt moving minimum filter
    chl_fil_tmp = movmax(chl_fil_tmp,win); % 11 pt moving minimum filter

    BBP_sr = BBP_fil_tmp+Float_Struc_BBP.BBP700_noise; % refractory + small particles
    BBP_l =  BBP_tot_tmp - BBP_sr; % large particels
    BBP_l(find(BBP_l<0)) = 0; % zero out negative BBP_l values. This may be possible in the case that the BBP700 signal is smaller than the instrument nosie of the float

    chl_sr = chl_fil_tmp+Float_Struc_Chl.chl_noise; % refractory + small particles
    chl_l =  chl_tot_tmp - chl_sr; % large particles
    chl_l(find(chl_l<0)) = 0; % as with BBP_l values, zero out negative chl_l values

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
  end
% apply correction to obtain BBP_s;
  Float_Struc_BBP.BBP700_s = Float_Struc_BBP.BBP700_s - P_bbp;  
  Float_Struc_BBP.BBP700_s(find(Float_Struc_BBP.BBP700_s<0)) = 0; 
  Float_Struc_BBP.BBP700_r = P_bbp; % refractory material is a single value for the entire float: like instrument noise
  Float_Struc_Chl.chl_s = Float_Struc_Chl.chl_s - P_chl;
  Float_Struc_Chl.chl_s(find(Float_Struc_Chl.chl_s<0)) = 0; 
  Float_Struc_Chl.chl_r = P_chl; % refractory material is a single value for the entire float: like instrument noise
  Float_Struc_BBP.PRES_ADJUSTED = Float_Struc_Chl.PRES_ADJUSTED;

  eval(['Data_bbp_proc.',char(fieldnm),'= Float_Struc_BBP;']);    
  eval(['Data_chl_proc.',char(fieldnm),'= Float_Struc_Chl;']);

end
%save([oneargo,'/matfiles/',filename,'_processed.mat'],[Data_bbp,'_proc'],[Data_chl,'_proc']); 

end

  
