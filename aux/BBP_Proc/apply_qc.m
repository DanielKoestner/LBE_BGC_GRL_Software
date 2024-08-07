% apply QC flag to array with all flags
function [QC_Flags, QC_1st_failed_test] = apply_qc(QC_Flags, ISBAD, QC, QC_1st_failed_test, QC_TEST_CODE);

global QC_settings
    % find which part of the QC_Flag[ISBAD] array needs to be updated with the new flag
    i2flag = find( QC_Flags(ISBAD) < QC );  % find where the existing flag is lower than the new flag (cannot lower existing flags)
    % apply flag
    QC_Flags(ISBAD(i2flag)) = QC;
   %  record which test changed the flag
    N_QC_TEST_CODE = find(strcmp(QC_TEST_CODE,QC_settings.QC_TEST_CODES));
    QC_1st_failed_test(N_QC_TEST_CODE,ISBAD) = N_QC_TEST_CODE;
    if strcmp(QC_TEST_CODE,'B') == 1 && QC_settings.do_noisy_esd == 1; % if using noisy profiles test with esd test
      N_QC_TEST_CODE_ESD = find(strcmp('I',QC_settings.QC_TEST_CODES));% index for ESD test
      ind_passed_ESD =  find(QC_1st_failed_test(N_QC_TEST_CODE_ESD,ISBAD) == 0);
      QC_Flags(ind_passed_ESD) = 2; % set QC flag for measurements that pass ESD test but fail noisy profile to 2
    end 
end
