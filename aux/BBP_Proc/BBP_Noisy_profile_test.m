function [QC_Flags,QC_1st_failed_test,res] = BBP_Noisy_profile_test(BBP, BBPmf1, PRES, QC_Flags, QC_1st_failed_test);
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
    do_med_out = QC_settings.do_med_out;
    FAILED = false(1);
    
    QC = 3; % flag to apply if the result of the test is true
    QC_TEST_CODE = 'B';
    ISBAD = [];  % index of where flags should be applied in the profile
 
    res = NaN.*ones(size(BBPmf1));
     
    innan = find(~isnan(BBPmf1) & ~isnan(BBP)); % due to adapative_medfilt1 function, it is possible that BBPmf1 contains NaNs in indices where BBP does not
    if length(innan)>10; % if we have at least 10 points in the profile
        
        %new constraint: noise should be below 100 dbars
        iPRES = find(PRES(innan) > B_PRES_THRESH);
        PRES_i = PRES(iPRES);
        if any(iPRES); % if there are enough points, then continue wiht the test
            res(innan) = abs(BBP(innan)-BBPmf1(innan));
            res_innan = res(innan);
            res_innan_iPRES = res_innan(iPRES);
            ioutliers = find(abs(res_innan(iPRES)) > B_RES_THRESHOLD); % index of where the rel res are greater than the threshold
            ibad = innan(iPRES(ioutliers));
            if length(ioutliers)/length(iPRES) >= B_FRACTION_OF_PROFILE_THAT_IS_OUTLIER; % this is the actual test: is there more than a certain fraction of points that are noisy?
            %if length(ioutliers) >= 1; % remove all noisy points
              ISBAD = ioutliers; 
            end
            if do_med_out == 1; % apply QC flag to all outliers, regardless what percentatge of data
              if length(ioutliers) >= 1;
                ISBAD = innan(iPRES(ioutliers));
              end
            end
        end% any(iPres)
    end %length(innan)>10;
    % update QC_Flags to 3 when bad profiles are found
    if size(ISBAD) ~= 0;
        FAILED = true(1);
        % apply flag
        %if length(ioutliers)/length(iPRES) >= B_FRACTION_OF_PROFILE_THAT_IS_OUTLIER; 
        if length(ioutliers)/length(iPRES) >= B_FRACTION_OF_PROFILE_THAT_IS_OUTLIER;
          [QC_Flags, QC_1st_failed_test] =  apply_qc(QC_Flags, find(~isnan(BBP)), QC, QC_1st_failed_test, QC_TEST_CODE);% ~isnan(BBP) is used to flag the entire profile
        end
        if do_med_out == 1;
          [QC_Flags, QC_1st_failed_test] = apply_qc(QC_Flags, ISBAD, QC, QC_1st_failed_test, QC_TEST_CODE);
        end
         
    end



end % function

