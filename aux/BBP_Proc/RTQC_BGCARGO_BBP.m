%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [Data_BBP_RTQC] = RTQC_BGCARGO_BBP(Data_BBP,path_oneargo,varargin);

% Conversion from Dall'Olmo's script (https://github.com/euroargodev/BBP_RTQC) to assign QC flags to BBP700 from BGC-ARGO floats 

%% Output

% Data_BBP_RTQC: a structure containing as field float ids. Each float id is itself a substurcture
% that contains 

% i. BBP700: raw BBP700 data from float (m^-1)
% ii. PRES: pressure data from float (dbar)
% iii. BBP700_QC: QC flags both originated from float data and derived from the RTQC processing. Flags include
%      1: good data (all QC tests passed)
%      2: probably good data (should be used with caution)
%      3: probably bad data (do not use until and expert has checked these data)
%      4: bad data (do not use these data)
%      5: value changed
%      8: interpolated data
%      9: missing data (Dal'Ollmo assigns this only when the entire profile has no data. When profiles have only a fraction of missing data points, or missing pressure levels, they retain QC=0 as provided)
%      Dal'Ollmo test assigns QC 2,3,4, or 9 if tests failed. The remaining flags are from the float data, or have QC 0 assigned based of if pressure is available
% iv. PARK_PRESS: parking pressure used for parking hook test
% v.  BBP700_RTQC_failed_tests: An pxmxn dimensional array, where p is the number of tests, m is the number of pressure levels, and n is the number of profiles. An element is assigned a value equal to the index of the
%     test code if a test fails, otherwise it is assigned 0. The sequence of test codes is described in the portion of the code above the variable "QC_test_codes". order mattesr, so e.g., in the mxn array for p=1, each element
%     is assigned a 1 if the negative value test fails, or 0 if it passes. Or for the p = 8, each element in the mxn will be assigned an 8 if the test fails, or a 0 if the test passes.


%% Input


% i. Data_BBP: A data structures containing the following fields:
%    following fields:

%    WMOID1 [1x1 structure]
%    WMOID2 [1x1 structure]
%    ...
%    WMOIDN [1x1 structure]

%    where WMOIDs are the float ids for the floats obtained from the BGC-ARGO repository
%    Each WMOID is itself a structure, and must contain the following field for Data_bbp:

%    BBP700 [NxM] (m^-1): bbp700 can be named anything, but the default names is BBP700
%    PRES [NxM] (dbar): pressure can also be named anything, but default name is PRES

%    where N is the number of pressure levels, and M is the number of profiles

% ii. path_OneArgo: path to OneArgo package.

%% Optional Arguments

% i.  A_MIN_BBP700: minimum acceptable BBP700 value 
% ii. A_MAX_FRAC_BAD: fraction of BBP700 values that need to be negative before the profile is flagged with QC=4
% iii.B_RES_THRESHOLD: threshold for residuals relative to median fitlered value before a datapoint is flagged as noisy
% iv. B_FRAC_OUTLIER: fraction of data that are considered noisy above which profile is flagged with QC=3
% v.  B_PRES_THRESH: threshold depth below which noisy profile test is applied
% vi. C_DEPTH_THRESH: depth threshold below which deep value test is applied
% vii.C_DEEP_BBP700_THRESH: threshod BBP700 above which a datapoint is flagged in deep value test
% viii.C_N_DEEP_POINTS: number of datapoints needed to apply deep value test
% viv.G_DELTAPRES1: difference in depth from parking pressure over which parking hook test is applied
% x.G_DELTAPRES2: difference in depth from parking pressure which is used to compute baseline BBP700 for parking hook test. Baseline computed from median between max_pres + G_DELTAPRES1 and max_pres + G_DELAPRES2
% xi. G_DEV: deviation in BBP700 from median of baseline in parking test above which a datapoint is flagged with QC=4
% xii. G_DELTAPRES0: deviation between maximum and parking pressure above which parking hook test will not be applied
% xiii. E_MIN_N_PERBIN: minimum number of data points per bin below which missing data flags bin with QC=3 or 4
% xiv. E_MAXPRES: maximum pressure which if not exceeded by a profile, it is considered a shallow profile and given QC=3
% xv. var_bbp: variable name of BBP700 contained in input data structure
% xvi. var_pres: variable name of pressure contained in input data structure
% xvii. do_noisy_bin: set to 1 if noisy profile test is applied in 100-m bins, and to 0 if not. Default is 0
% xvii. do_noisy_esd: set to 1 if Generalized Extreme Studentized Deviate test is applied in 100-m bins, and to 0 if not. Default is 0
% xviii. alpha_sf: alpha for for ESD test in upper 100dbar. Default is 0.01
% xviv. alpha_subsf: alpha for for ESD test in 100-300dbar interval. Default is 0.05
% xx. alpha_subsf: alpha for for ESD test below 300 dbar Default is 0.10


 
% The methodology to process BBP data follows Dal'Ollmo et al., 2022, Open Research Europe 
% Paul Lerner, Dec 7, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% defaults of optional arguments: assume OneArgo and matfiles are directories in the level up the "scripts" directory

% path to OneArgo package
addpath(path_oneargo);



% Parse optional arguments for global parameters

global QC_settings

% Negative BBP
A_MIN_BBP700 = 0;    % [1/m]
%A_MAX_BBP700 = 0.01 % [1/m] REVISED VALUE (very conservative estimate based on histograms in fig 2 of Bisson et al., 2019, 10.1364/OE.27.030191)
%A_MAX_BBP700 = 0.03 #%[1/m] REVISED VALUE (very conservative estimate based on histograms in fig 2 of Bisson et al., 2019, 10.1364/OE.27.030191)
A_MAX_FRACTION_OF_BAD_POINTS = 0.1;


% Noisy profile
B_RES_THRESHOLD = 0.0005;                    % [1/m] threshold for relative residuals
B_FRACTION_OF_PROFILE_THAT_IS_OUTLIER = .1; % fraction of profile with relative residuals above RES_THRESHOLD: set default to 5%, while Dal'Olmo 2022 is 10%
B_PRES_THRESH = 100;                         % [dbar] this is to avoid flagging profiles with spikes in surface data (likely good data)

% High-Deep Value
C_DEPTH_THRESH = 700;          % [dbar] pressure threshold below which the test acts
C_DEEP_BBP700_THRESH = 0.0005; % [1/m] threshold for bbp at depth
C_N_DEEP_POINTS = 5;        % number of points deeper than C_DEPTH_THRESH required for the test to proceed

% Parking hook
G_DELTAPRES1 = 50;  % [dbar] difference in PRES from parking pressure over which the test is implemented
G_DELTAPRES2 = 20;  % [dbar] difference in PRES from parking pressure use to compute test baseline
G_DEV = 0.0002;     % [1/m] deviation from baseline that identifies anomalous data points
G_DELTAPRES0 = 100; % [dbar] define how close PARK_PRES has to be to max(PRES)
                   %        for the test to proceed

% Missing Data
E_MIN_N_PERBIN = 1; % [-] minimum number of data points per bin
E_MAXPRES = 200;   % [dbar] pressure below which the profile is considered shallow

% variable to go through QC procedure. For BGC-ARGO arrays, this will be the BBP700 field of the data structure 
% in "Data_bbp". However, the user may want to specify a different field name: for example, if the data has been
% processed into small and large particles, then they may wish to QC the separet BBP700-particle classes.

var_bbp = 'BBP700';
var_pres = 'PRES';
do_noisy_bin = 0;
do_noisy_esd = 0;
do_med_out = 0;
flag_samplesize = 0;
alpha_sf = 0.01;
alpha_subsf = 0.05;
alpha_meso = 0.10;
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'A_MIN_BBP700')
        A_MIN_BBP700 = A_MIN_BBP700{i+1};
    elseif strcmpi(varargin{i}, 'A_MAX_FRAC_BAD');
        A_MAX_FRACTION_OF_BAD_POINTS = varargin{i+1};
    elseif strcmpi(varargin{i}, 'B_RES_THRESHOLD')
        B_RES_THRESHOLD = varargin{i+1}
    elseif strcmpi(varargin{i}, 'B_FRAC_OUTLIER')
        B_FRACTION_OF_PROFILE_THAT_IS_OUTLIER = varargin{i+1}
    elseif strcmpi(varargin{i}, 'B_PRES_THRESH')
        B_PRES_THRESH = varargin{i+1}
    elseif strcmpi(varargin{i}, 'C_DEPTH_THRESH')
        C_DEPTH_THRESH= varargin{i+1}
    elseif strcmpi(varargin{i}, 'C_DEEP_BBP700_THRESH')
        C_DEEP_BBP700_THRESH= varargin{i+1}
    elseif strcmpi(varargin{i}, 'C_N_DEEP_POINTS')
        C_N_DEEP_POINTS= varargin{i+1}
    elseif strcmpi(varargin{i}, 'G_DELTAPRES1')
        G_DELTAPRES1= varargin{i+1}
    elseif strcmpi(varargin{i}, 'G_DELTAPRES2')
        G_DELTAPRES2= varargin{i+1}
    elseif strcmpi(varargin{i}, 'G_DEV')
        G_DEV = varargin{i+1}
    elseif strcmpi(varargin{i}, 'G_DELTAPRES0')
        G_DELTAPRES0 = varargin{i+1}
    elseif strcmpi(varargin{i}, 'E_MIN_N_PERBIN')
        E_MIN_N_PERBIN = varargin{i+1}
    elseif strcmpi(varargin{i}, 'E_MAXPRES')
        E_MAXPRES = varargin{i+1}
    elseif strcmpi(varargin{i}, 'var_bbp')
        var_bbp = varargin{i+1};
    elseif strcmpi(varargin{i}, 'var_pres')
        var_pres = varargin{i+1};
    elseif strcmpi(varargin{i}, 'do_noisy_bin')
        do_noisy_bin = varargin{i+1};
    elseif strcmpi(varargin{i}, 'do_noisy_esd')
        do_noisy_esd = varargin{i+1};
    elseif strcmpi(varargin{i}, 'do_median_out')
        do_med_out = varargin{i+1};
    elseif strcmpi(varargin{i}, 'flag_samplesize')
        flag_samplesize = varargin{i+1}
    elseif strcmpi(varargin{i}, 'alpha_sf')
        alpha_sf = varargin{i+1}
    elseif strcmpi(varargin{i}, 'alpha_subsf')
        alpha_subsf = varargin{i+1}
    elseif strcmpi(varargin{i}, 'alpha_meso')
        alpha_meso = varargin{i+1}
    end
end
% Test codes

% tests are:
% A: negative values above 5 dbar
% A2: negative values below 5 dbar
% B: noisy profile test
% C: High deep value test
% E: missing data test
% G: parking hook test
% H: noisy profile test for individual bins
% I: outlier test using generalized extreme studentized deviate test 
QC_test_codes = {'A','A2','B','C','E','G','H','I'}; % indices for test are 1:8, where 1 is for negative values test above 5 dbar, and 8 is for outlier detection using generalized ESD test

% load dataset: discontinued for now, but maybe should be added back in as optional argument

%load([path_matfiles,'/',Data,'.mat']);

% assign medfilt class object

medfilt_obj = medfilt;

% have global setting variables to easily share among functions

% set global settings for fixed values
QC_settings.A_MIN_BBP700 = A_MIN_BBP700;
QC_settings.A_MAX_FRACTION_OF_BAD_POINTS = A_MAX_FRACTION_OF_BAD_POINTS;
QC_settings.B_RES_THRESHOLD = B_RES_THRESHOLD;
QC_settings.B_FRACTION_OF_PROFILE_THAT_IS_OUTLIER = B_FRACTION_OF_PROFILE_THAT_IS_OUTLIER;
QC_settings.B_PRES_THRESH = B_PRES_THRESH;
QC_settings.C_DEPTH_THRESH = C_DEPTH_THRESH;
QC_settings.C_DEEP_BBP700_THRESH = C_DEEP_BBP700_THRESH;
QC_settings.C_N_DEEP_POINTS = C_N_DEEP_POINTS;
QC_settings.G_DELTAPRES1 = G_DELTAPRES1;
QC_settings.G_DELTAPRES2 = G_DELTAPRES2;
QC_settings.G_DEV = G_DEV;
QC_settings.G_DELTAPRES0 = G_DELTAPRES0;
QC_settings.E_MIN_N_PERBIN = E_MIN_N_PERBIN;
QC_settings.E_MAXPRES = E_MAXPRES;
QC_settings.QC_TEST_CODES = QC_test_codes;
QC_settings.do_med_out = do_med_out;
QC_settings.do_noisy_esd = do_noisy_esd;
QC_settings.flag_samplesize = flag_samplesize; % must be <10
QC_settings.alpha_sf = alpha_sf;
QC_settings.alpha_subsf = alpha_subsf;
QC_settings.alpha_meso = alpha_meso;

%% loop over floats

fieldnms = fieldnames(Data_BBP); % these are the float ids


% find integer for missing data test, which is used to flag all NaN values

N_QC_MISSING = 0; % this is to flag data that come from GDAC as NaN 
for i=1:length(fieldnms);%24
  fieldnm = fieldnms(i);
  fieldnm 
  fprintf(['\n Starting Float',char(fieldnm),'...\n']);



% Extract the float structure  for float i from the datasets
 
  Float_Struc_BBP = eval(['Data_BBP.',char(fieldnm)]);


% Initialize arrays of QC flag values

%% loop over profiles

  for j=1:length(Float_Struc_BBP.(var_pres)(1,:));%108
    BBP_tot_tmp = eval(['Float_Struc_BBP.',var_bbp,'(:,j);']);
    Pres_tmp = Float_Struc_BBP.(var_pres)(:,j);
    eval(['QC_Flags = Float_Struc_BBP.',var_bbp,'_QC(:,j);']);
% as a first step, NaN out BBP where Pres < 0. Set QC=0, since this is QC for values from float that are NaN
    Pres_tmp(find(Pres_tmp<0)) = NaN;
    BBP_tot_tmp(find(Pres_tmp<0)) = NaN;
    QC_Flags(find(Pres_tmp<0)) = 0;
    indnonan = find(~isnan(BBP_tot_tmp) & ~isnan(Pres_tmp));
    indnan = find(isnan(BBP_tot_tmp) | isnan(Pres_tmp));
    QC_Flags(indnan) = 0; % set QC flag of BBP where pressure is not defined to zero
    Pres_tmp = Pres_tmp(indnonan);
    BBP_tot_tmp = BBP_tot_tmp(indnonan);
    BBP_QC_failed_test = zeros(length(QC_test_codes),length(Pres_tmp));
  
    QC_Flags2 = QC_Flags; % this array is used to populate "BBP700_QC" field at the end, because now this has been modified based on criteria above 
    QC_Flags = QC_Flags(indnonan); % for QC_Flags passed to tests, only take values where BBP700 and PRES are not NaN
    QC_Flags3 = QC_Flags; % the initial QC_Flag array must be set to 1 for tests to work. This array will then replace indices where all tests passed, but QC was originally 8 (but are now temporarily 1)
    QC_Flags(find(QC_Flags>0)) = 1; % Aside from QC 8, we all QC flags are to be reasssessed here: so we start by "assuming" all data are "good", i.e. QC=1

% this array stores elements that indicate whether a test was passed or failed. For a given member of dimension p, whose length equals the number of tests, we assign 0 if a test is passed and the value X if the test failed, where X is
% the index for that test
    eval(['Float_Struc_BBP.',var_bbp,'_RTQC_failed_tests(:,:,j) = zeros(length(QC_test_codes),length(Float_Struc_BBP.(var_pres)(:,j)),1);']);

% compute median filter BBP
    BBPmf1 = medfilt_obj.adaptive_medfilt1(Pres_tmp, BBP_tot_tmp);
% find parking pressure for parking hook test
    park_press = Float_Struc_BBP.PARK_PRESS(j);


% first do ESD test, if applicable
    if do_noisy_esd == 1
      [QC_Flags, BBP_QC_failed_test] = BBP_Noisy_profile_test_esd(BBP_tot_tmp, BBPmf1, Pres_tmp, QC_Flags, BBP_QC_failed_test);
    end

% next to noisy profile test
% this test is applied to the entire profile: if >10% of the profile is "noisy", it is assumed the sensor for the cast is compromised
    [QC_Flags, BBP_QC_failed_test,tmp] = BBP_Noisy_profile_test(BBP_tot_tmp, BBPmf1, Pres_tmp, QC_Flags, BBP_QC_failed_test);
% this test is applied to 100m bins. In any inidivual bin, if 10% of the data is noisy, the data in that bin is assumed questionable
    if do_noisy_bin == 1 
      [QC_Flags, BBP_QC_failed_test,tmp] = BBP_Noisy_profile_bin_test(BBP_tot_tmp, BBPmf1, Pres_tmp, QC_Flags, BBP_QC_failed_test);
    end

% first do negative value test

    [QC_Flags, BBP_QC_failed_test] = BBP_Negative_BBP_test(BBP_tot_tmp, Pres_tmp, QC_Flags, BBP_QC_failed_test);
    
% next to parking hook test
    [QC_Flags, BBP_QC_failed_test] = BBP_Parking_hook_test(BBP_tot_tmp, BBPmf1, Pres_tmp, max(Pres_tmp), park_press,QC_Flags,BBP_QC_failed_test);


% next do high deep value test
    [QC_Flags, BBP_QC_failed_test] = BBP_High_Deep_Values_test(BBPmf1, Pres_tmp, QC_Flags, BBP_QC_failed_test);
   
%  do missing data test
    [QC_Flags,BBP_QC_failed_test] = BBP_Missing_Data_test(BBP_tot_tmp, Pres_tmp, max(Pres_tmp), QC_Flags, BBP_QC_failed_test);
% now we replace indices where QC_Flag was originally 8 and where all tests are passed
    %eval(['Float_Struc_BBP.',var_bbp,'_RTQC_Flags(indnonan,j) = QC_Flags;']);
    eval(['Float_Struc_BBP.',var_bbp,'_QC(:,j) = QC_Flags2;']); % replace QC flags with those for no pressure or negative pressure layres   
    eval(['Float_Struc_BBP.',var_bbp,'_QC(indnonan,j) = QC_Flags;']);
    eval(['Float_Struc_BBP.',var_bbp,'_RTQC_failed_tests(:,indnonan,j) = BBP_QC_failed_test;']);
  end 
  eval(['Data_BBP_RTQC.',char(fieldnm),'= Float_Struc_BBP;']);  
 
end
%save([oneargo,'/matfiles/',filename,'_processed.mat'],[Data_bbp,'_proc'],[Data_chl,'_proc']); 
end

  
