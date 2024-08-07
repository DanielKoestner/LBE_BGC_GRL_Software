function [QC_Flags,QC_1st_failed_test] = BBP_Noisy_profile_test(BBP, BBPmf1, PRES, QC_Flags, QC_1st_failed_test);
% BBP: nparray with all BBP data
% BBPmf1: smooth BBP array (medfilt1(BBP700, 31)
% QC_Flags: array with QC flags
% QC_flag_1st_failed_test: array with info on which test failed QC_TEST_CODE
                                                                         
% Objective: To flag profiles that are affected by noisy data. This noise could 
% indicate sensor malfunctioning, some animal spikes, or other anomalous conditions.

% Implementation: The absolute residuals between the median filtered BBP and the 
% raw BBP values are computed below a pressure threshold B_PRES_THRESH = 100 dbar 
% (this is to avoid surface data, where spikes are more common and generate false 
% positives). The test fails if residuals with values above B_RES_THRESHOLD = 0.0005 m-1 
% occur in at least B_FRACTION_OF_PROFILE_THAT_IS_OUTLIER = 10% of the profile. 
% These threshold values were selected after visual inspection of flagged profiles.
% Flagging: If the test fails, a QC flag of 3 is assigned to the entire profile.                                                                   
% __________________________________________________________________________________________

    global QC_settings

    B_PRES_THRESH = QC_settings.B_PRES_THRESH;
    B_RES_THRESHOLD = QC_settings.B_RES_THRESHOLD;
    B_FRACTION_OF_PROFILE_THAT_IS_OUTLIER = QC_settings.B_FRACTION_OF_PROFILE_THAT_IS_OUTLIER;
    flag_samplesize = QC_settings.flag_samplesize;
    alpha_sf = QC_settings.alpha_sf;
    alpha_subsf = QC_settings.alpha_subsf;
    alpha_meso = QC_settings.alpha_meso;

                                                                     
    FAILED = false(1);
    
    QC = 3; % flag to apply if the result of the test is true
    QC_TEST_CODE = 'I';
    ISBAD = [];  % index of where flags should be applied in the profile
    ibad = [];
    bins = linspace(0, 1000, 11); % create 10 bins between B_PRES_THRESH and 1000 dbars
    bin_counts = zeros(length(bins)-1,1); % initialise array with number of counts in each bin
    for i = 2:length(bins);
      bin_counts(i-1) = length(find((PRES >= bins(i-1)) & (PRES < bins(i))));
      iPRES_bin = find((PRES >= bins(i-1)) & (PRES < bins(i))); % this is for flagging all datapoints in bin
% set alpha depending on depth interval considered. Higher alpha at deeper depths because we expect a smaller BBP spread at depth, rejecting values the fall outside a smaller confidence interval as "noisy"
      if bins(i) <= 100;
        alpha = alpha_sf;
      elseif bins(i) <= 300;
        alpha = alpha_subsf;
      else
        alpha = alpha_meso;
      end 
      if length(iPRES_bin) > 10;
% Rosner, B. (1983). "Percentage points for a generalized ESD many-outlier
%    procedure". Technometrics 25.2, pp. 165-172.
% since values >0, a lognormal distribtuion is assmued
          BBP_esd = log(BBP(iPRES_bin));
          N = length(BBP_esd);
          k = 0.1*N; % maximum number of outliers is 10% of data
          X = 0*(1:k);
          G = 0*(1:k);
          R = 0*(1:k);
          for ii = 1:k
            BBPm = nanmean(BBP_esd);
            s = nanstd(BBP_esd);
            r = abs(BBP_esd-BBPm)/s;
            [rmax,index] = nanmax(r);
            R(ii) = rmax;
            X(ii) = BBP_esd(index);
            BBP_esd(index) = [];
            L = ii-1;
            p = 1-alpha/(2*(N-L));
            df = N-L-2;
            t = tinv(p,df);
            gamma = ((N-L-1)*t)/sqrt((df+t^2)*(N-L));
            G(ii) = gamma;
          end
          II = 1:k;
          P = max(II(R > G));
          if P == 0
            I = [];
          else
            [~,I,~] = intersect(log(BBP(iPRES_bin)),X(1:P));
          end
          ioutliers = I; 
          ibad = [ibad;iPRES_bin(ioutliers)]; % flag only those values that are identified as outliers
      elseif flag_samplesize ~= 0;
        if length(iPRES_bin) <= flag_samplesize;
          ibad = [ibad;iPRES_bin]; % if there is not enough data, flag entire bin as 'outlier'
        end
      end  %length(iPRES_bin)>10
    end
    ISBAD = ibad; 
     
    % update QC_Flags to 3 when bad profiles are found
    if size(ISBAD) ~= 0;
        FAILED = true(1);
        % apply flag
          [QC_Flags, QC_1st_failed_test] = apply_qc(QC_Flags, ISBAD, QC, QC_1st_failed_test, QC_TEST_CODE);
    end



end % function

