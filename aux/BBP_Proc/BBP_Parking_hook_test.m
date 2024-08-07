function [QC_Flags,QC_1st_failed_test] = BBP_Parking_hook_test(BBP, BBPmf1, PRES, maxPRES, parkPRES, QC_Flags, QC_1st_failed_test);
% BBP: nparray with all BBP data
% BBPmf1: nparray with medfilt BBP data
% maxPRES: maximum pressure recorded in this profile
% parkPRES: programmed parking pressure for this profile%
% QC_Flags: array with QC flags
% QC_flag_1st_failed_test: array with info on which test failed QC_TEST_CODE
   
% Objective: To flag data points near the parking pressure with 
% anomalously high values, when the parking pressure is close to the 
% maximum pressure of the profile. This could indicate that particles 
% have accumulated on the sensor or the float and that are released when the float starts ascending.

% Implementation: First the parking pressure (parkPRES) is extracted 
% from the metadata file. Then, we verify that the vertical resolution of
% the data near parkPRES is greater than G_DELTAPRES2 = 20 dbar: if it is not, 
% the test cannot be applied to this profile. If the vertical resolution is sufficient, 
% we verify that the maximum pressure of the profile is less than G_DELTAPRES0 (100 dbar) different 
% from parkPRES (i.e., that parkPRES ~= max(PRES)), i.e., that the profile starts 
% from the parking pressure. If it does, a pressure range iPRESmed 
% (max(PRES) - G_DELTAPRES2 > PRES >= max(PRES) - G_DELTAPRES1, with G_DELTAPRES1 = 50 dbar) 
% is defined over which the baseline for the test will be calculated. This baseline 
% is computed as the median + G_DEV (with G_DEV =  0.0002 m-1). The test is implemented 
% in the pressure range iPREStest (where PRES>= maxPRES - G_DELTAPRES1). The test fails 
% if BBP within iPREStest is greater than the baseline.
 
% Flagging: A QC flag of 4 is applied to the points that fail the test.
 
% __________________________________________________________________________________________
    global QC_settings

    G_DELTAPRES1 = QC_settings.G_DELTAPRES1;  % [dbar] difference in PRES from parking pressure over which the test is implemented
    G_DELTAPRES2 = QC_settings.G_DELTAPRES2;  % [dbar] difference in PRES from parking pressure use to compute test baseline
    G_DEV = QC_settings.G_DEV;     % [1/m] deviation from baseline that identifies anomalous data points
    G_DELTAPRES0 = QC_settings.G_DELTAPRES0; % [dbar] define how close PARK_PRES has to be to max(PRES)

    FAILED = false(1);

    QC = 4;
    QC_TEST_CODE = 'G';
    ISBAD = []; % flag for noisy profile

    if isnan(maxPRES) | isnan(parkPRES);
      fprintf('no parking pressure data; parking hook test will not be applied\n');
      return
    end

    if length(PRES) < 2
      fprintf('not enough data for parking hook test\n');
      return
    end

    % check that there are enough data to run the test
    imaxPRES = find(PRES == maxPRES);
    deltaPRES = PRES(imaxPRES) - PRES(imaxPRES-1);
    if deltaPRES > G_DELTAPRES2;
        fprintf('vertical resolution is too low to check for Parking Hook; test will not be applied\n')
        return
    end

    % check if max PRES is 'close' (i.e., within 100 dbars) to parkPRES
    if abs(maxPRES - parkPRES) >= G_DELTAPRES0;
      return
    end

    % define PRES range over which to compute the baseline for the test
    iPRESmed = find((PRES >= maxPRES - G_DELTAPRES1 ) & (PRES < maxPRES - G_DELTAPRES2) );
    % define PRES range over which to apply the test
    iPREStest = find((PRES >= maxPRES - G_DELTAPRES1 ));

    % compute parameters to define baseline above which test fails
    medBBP = nanmedian(BBP(iPRESmed));
    baseline = medBBP + G_DEV;

    % this is the test
    ibad = find(BBP(iPREStest) > baseline);
    ISBAD = iPREStest(ibad);

    if size(ISBAD) ~= 0; % If ISBAD is not empty
        FAILED = true(1);
        % apply flag
        [QC_Flags, QC_1st_failed_test] = apply_qc(QC_Flags, ISBAD, QC, QC_1st_failed_test, QC_TEST_CODE);
    end
end %function
