function [QC_Flags,QC_1st_failed_test] = BBP_Missing_Data_test(BBP, PRES, maxPRES, QC_Flags, QC_1st_failed_test)
% BBP: nparray with all BBP data
% QC_Flags: array with QC flags
% QC_flag_1st_failed_test: array with info on which test failed QC_TEST_CODE
% fn: name of corresponding B-file
    
% Objective: To detect and flag profiles that have a large fraction of 
% missing data. Missing data could indicate shallow or incomplete profiles.

% Implementation: The upper 1000 dbar of the profile are divided into 10 
% pressure bins with the following lower boundaries (all in dbar): 
% 50, 156, 261, 367, 472, 578,  683, 789, 894, 1000. For example, 
% the first bin covers the pressure range [0, 50), the second [51, 156), 
% etc. The test fails if any of the bins contains fewer data points than MIN_N_PERBIN = 1.
 
% Flagging: Different flags are assigned depending on how many bins are empty.
% If only one bin contains data or the profile has no data at all, a QC flag of 4 
% is applied to the entire profile. This condition may indicate a malfunctioning 
% sensor or a profile that is so shallow that it is too difficult to quality control in real time.
% If there are bins with missing data, but the number of bins with data is 
% greater than one, then a QC flag of 3 is assigned to the entire profile.
% __________________________________________________________________________________________

    global QC_settings

    E_MIN_N_PERBIN = QC_settings.E_MIN_N_PERBIN;
    E_MAXPRES = QC_settings.E_MAXPRES;
    FAILED = false(1);

    QC_all = NaN.*ones(4,1);
    QC_all(1) = 3; % flag to apply if shallow profile
    QC_all(2) = 4; % flag to apply if the result of the test is true only in one bin
    QC_all(3) = 3; % flag to apply if the result of the test is true elsewhere
    QC_all(4) = 9; % flag to apply if there are no data at all

    QC_TEST_CODE = 'E';
    ISBAD = [];  % index of where flags should be applied in the profile


    % bin the profile into 100-dbars bins
    bins = linspace(0, 1000, 11); % create 10 bins between 0 and 1000 dbars
    bin_counts = zeros(length(bins)-1,1); % initialise array with number of counts in each bin
    for i = 2:length(bins);
      bin_counts(i-1) = length(find((PRES >= bins(i-1)) & (PRES < bins(i))));
    end
    % check if there are bins with missing data
    if any(bin_counts < E_MIN_N_PERBIN);
        isbad4plot = find(bin_counts < E_MIN_N_PERBIN);
        ISBAD = find(~isnan(BBP)); % flag the entire profile

        % find which bins contain data
        nonempty = find(bin_counts > 0); % index of bins that contain data points

        % select which flag to use
        if size(nonempty) ~= 0;

            % if shallow profile
            if (maxPRES < E_MAXPRES) & (length(find(bin_counts > E_MIN_N_PERBIN)) > 1);
                QC = QC_all(1);
            % if there is only one bin with data then
            elseif length(find(bin_counts > E_MIN_N_PERBIN)) == 1;  % with test
                QC = QC_all(2);

            % if missing data profile
            else;
                QC = QC_all(3);
                
            end
        else % this is for when we have no data at all, then
            QC = QC_all(4);

        end
    end
    if size(ISBAD) ~= 0; % if ISBAD, then apply QC_flag
        FAILED = true(1);
        % apply flag
        [QC_Flags, QC_1st_failed_test] = apply_qc(QC_Flags, ISBAD, QC, QC_1st_failed_test, QC_TEST_CODE);% np.where(BBP)[0] is used to flag the entire profile

    end

end % function

