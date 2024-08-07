function [QC_Flags, QC_1st_failed_test] = BBP_Negative_BBP_test(BBP, PRES, QC_Flags, QC_1st_failed_test);
% BBP: nparray with all BBP data
% QC_Flags: array with QC flags
% QC_flag_1st_failed_test: array with info on which test failed QC_TEST_CODE

% Objective: To flag data points or profiles outside an expected range of BBP values.
% The expected range is defined by two extrema: A_MIN_BBP700 = 0 m-1 and A_MAX_BBP700 = 0.01 m-1.
% A_MIN_BBP700 is defined to flag negative values, while A_MAX_BBP700 is a conservative
% estimate of the maximum BBP to be expected in the open ocean, based on statistics of
% satellite and BGC-Argo data (Bisson et al., 2019).
% Implementation: The test is implemented on data that have been median filtered (to remove spikes).
% Flagging: A QC flag of 3 is assigned to data points that fall above A_MAX_BBP700, while the
% entire profile is flagged with QC = 3 if any data point falls below A_MIN_BBP700
% (this is to reflect the more serious condition of having negative median filtered data in a profile).
% __________________________________________________________________________________________
%

    global QC_settings
    A_MIN_BBP700 = QC_settings.A_MIN_BBP700;
    A_MAX_FRACTION_OF_BAD_POINTS = QC_settings.A_MAX_FRACTION_OF_BAD_POINTS;

    FAILED = false(1);

    QC_TEST_CODE = 'A'; % or 'A2' if negative medfilt1 value is found

    % this is the test
    iLT5dbar = find( PRES < 5 ); % index for data shallower than 5 dbar
    i_ge5dbar = find( PRES >= 5 );% index for data deeper than or at 5 dbar
    ISBAD = find( BBP < A_MIN_BBP700 ); % first fill in all ISBAD indices where BBP < threshold
    ISBAD_ge5dbar = ISBAD(~ismember(ISBAD,iLT5dbar)); %  select only ISBAD indices deeper or equal than 5 dbar
    ISBAD_lt5dbar = ISBAD(ismember(ISBAD,iLT5dbar)); % select only ISBAD indices shallower than 5 dbar

    if length(ISBAD_ge5dbar) ~= 0;% If ISBAD_gt5dbar is not empty
        FAILED = true(1);
        QC_TEST_CODE = 'A2';
        % flag based on fraction of bad points
        fraction_of_bad_points = length(ISBAD_ge5dbar)/length(i_ge5dbar);
        if fraction_of_bad_points > A_MAX_FRACTION_OF_BAD_POINTS;
            QC = 4;
        else;
            QC = 3;
        end
        ISBAD = find(~isnan(BBP));  % flag entire profile

        % apply flag
        [QC_Flags, QC_1st_failed_test] = apply_qc(QC_Flags, ISBAD, QC, QC_1st_failed_test, QC_TEST_CODE);
    end



    %if (len(ISBAD_lt5dbar) > 0) & (len(ISBAD_ge5dbar)==0):  # if there are bad points only at PRES <5 dbar
    if length(ISBAD_lt5dbar) > 0;  % if there are bad points at PRES <5 dbar
        FAILED = true(1);
        QC = 4;
        QC_TEST_CODE = 'A';
        ISBAD = ISBAD_lt5dbar; % flag only negative values shallower than 5 dbar
        % apply flag
        [QC_Flags, QC_1st_failed_test] = apply_qc(QC_Flags, ISBAD, QC, QC_1st_failed_test, QC_TEST_CODE);
    end


end % function
