function [QC_Flags, QC_1st_failed_test] = BBP_High_Deep_Values_test(BBPmf1, PRES, QC_Flags, QC_1st_failed_test);
% BBP: nparray with all BBP data
% BBPmf1: smooth BBP array (medfilt1(BBP700, 31)
% QC_Flags: array with QC flags
% QC_flag_1st_failed_test: array with info on which test failed QC_TEST_CODE
% fn: name of corresponding B-file
% PLOT: flag to plot results
% SAVEPLOT: flag to save plot
% VERBOSE: flag to display verbose output

% Objective: To flag profiles with anomalously high BBP values at depth. 
% These high values at depth could indicate a variety of problems, including 
% biofouling, incorrect calibration coefficients, sensor malfunctioning, etc. 
% A threshold value of 0.0005 m-1 was selected that is half of the value typical
% for surface BBP in the oligotrophic ocean: median-filtered BBP data at depth
% are expected to be considerably lower than this threshold value.

% Implementation: This tests fails if the BBP profile has at 
% least a certain number (C_N_DEEP_POINTS = 5) of points 
% below a threshold depth (C_DEPTH_THRESH = 700 dbar) and if the median 
% of these deep median-filtered BBP values is above C_DEEP_BBP700_THRESH. 
% Note that this test can only be implemented 
% if the profile reaches a maximum pressure greater than 700 dbar. 
    
% Flagging: If the test fails, a QC flag of 3 is applied to the entire profile. 
% High deep BBP values can result from a variety of reasons, including natural causes, 
% in which case data might be set to good quality during DMQC. Therefore, we decided
% to use QC=3 and to revise these profiles during DMQC.
%
% __________________________________________________________________________________________

    global QC_settings

    C_DEPTH_THRESH = QC_settings.C_DEPTH_THRESH;
    C_N_DEEP_POINTS = QC_settings.C_N_DEEP_POINTS;
    C_DEEP_BBP700_THRESH = QC_settings.C_DEEP_BBP700_THRESH;
    FAILED = false(1);

    QC = 3; % flag to apply if the result of the test is true
    QC_TEST_CODE = 'C';
    ISBAD = []; % flag for noisy profile
    % this is the test
    iDEEP = find(PRES > C_DEPTH_THRESH); % find deep part of the profile
    nPointsBelow700dbar = length(iDEEP); % number of points below C_DEPTH_THRESH 
    if nPointsBelow700dbar >= C_N_DEEP_POINTS; % check if we have enough points at depth  
        if nanmedian(BBPmf1(iDEEP)) > C_DEEP_BBP700_THRESH; % check if median value at depth is greater than threshold
            ISBAD = find(~isnan(PRES)); % flag entire profile
        end
    end
    if size(ISBAD) ~= 0; % if ISBAD, then apply QC_flag=3
        FAILED = true(1);
        % apply flag
        [QC_Flags, QC_1st_failed_test] = apply_qc(QC_Flags, ISBAD, QC, QC_1st_failed_test, QC_TEST_CODE);
    end
end % function

