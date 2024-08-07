function [BBP700_noise_mon,chl_noise_mon] = extract_monnoise(time,yearmonths,bbp700,chl,pres,win);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to extract monthly varying background noise for BGC-ARGO floats

%% Output

% iii. BBP700_noise_,p:  NxM array of BBP700 background noise estimates, where N is the number of pressure levels and M is the number of profiles in bbp700 (m^-1)
% iv.  chl_noise_on: NxM array of chl-a background noise estimates, where N is the number of pressure levels and M is the number of profiles in chl (mg/m^3)


%% Input

% i.    time: a 1xM datetime array (YYYY-MM-DD), which contains the dates during which profiles are sampled
% ii.   yearmonths: 1xMr array, where Mr is the number of non-repeating months of the year (YYYYMM) of the float's lifetime
% iii.  pres: NxM array of pressure levels, where is the number of pressure levels (same for each profile, but pressure for each profile can vary from j=1,2..M), and M is the number of profiles for pres (dbar)
% iv.   bbp700:  NxM array of pressure levels, where is the number of pressure levels, and M is the number of profiles for BBP700 (m^-1)
% v.    chl:  NxM array of pressure levels, where is the number of pressure levels, and M is the number of profiles for chl-a (mg/m^3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global proc_settings

BBP700_noise_arr = NaN.*ones(size(bbp700));
chl_noise_arr = NaN.*ones(size(chl));

yearmon_flag = NaN.*ones(length(yearmonths),1); % flag for if data exists or not
% if there is no data, remove output
if sum(~isnan(bbp700(:))) == 0 || sum(~isnan(chl(:))) == 0;
  fprintf('Not enough chla or BBP700 values found to compute noise for this float. Noise Values set to NaN \n')
else
% first find months of the year with enough data to obtain noise estimate 
  for k =1:length(yearmonths);
    indmon = (find(time == yearmonths(k))); %  % extract all occurances current year
    bbp700_mon = bbp700(:,indmon); % subsample bbp700 at current year
    chl_mon = chl(:,indmon); % subsample chl-a at current month of the year
    pres_mon = pres(:,indmon); % subsample pressure at current month of the year

% compute background noise for the current month of the current year

    [BBP700_noise,chl_noise] = calc_bckg_noise(pres_mon,bbp700_mon,chl_mon,win);

% fill in noise arrays if noise was computed    
    if ~isnan(BBP700_noise(:))>0 && ~isnan(chl_noise(:))>0
      BBP700_noise_arr(:,indmon) = BBP700_noise;
      chl_noise_arr(:,indmon) = chl_noise;
      BBP700_noise_all(k) = BBP700_noise;
      chl_noise_all(k) = chl_noise;
      yearmon_flag(k) = 1;  % set flag for current month of the year to true (1) if noise estimates are available
    else

% otherwise set noise for current year to NaN
      BBP700_noise_all(k) = NaN;
      chl_noise_all(k) = NaN;
    end
  end
 
   
  if all(isnan(BBP700_noise_all));
    fprintf('Not enough chla or BBP700 values found to compute noise for this float. Noise values set to NaN. \n') 
  else

    yearmon_data = yearmonths(find(~isnan(yearmon_flag))); % find years where noise estimates exist

% fill out remaining months  with noise estimates of closest year
    for k =1:length(yearmonths);
      if isnan(yearmon_flag(k));
        ind_data = find(abs(yearmonths(k) - yearmon_data) == min(abs(yearmonths(k) - yearmon_data)));% find index of closest month with available noise estimate
        idx=randperm(length(ind_data),1);
        ind_data = ind_data(idx); % if two such months are available, pick one at random
        indmon_noise = (find(yearmonths == yearmon_data(ind_data))); % find the index for when the full yearmonth array equals the array containing month of real noise estiamtes
        indmon = (find(time == yearmonths(k))); % find all instances of the datetime array corresponding to the current month
        BBP700_noise_arr(:,indmon) = BBP700_noise_all(indmon_noise);% set bbp700 noise arr for profiles at the current month equal to the noise estimate available from the closets month
        chl_noise_arr(:,indmon) = chl_noise_all(indmon_noise); % set chl-a noise arr for profiles at the current month equal to the noise estimate available from the closets month
      end
    end

  end % if all BBP700_noise values are NaN
end %if all data are NaN


BBP700_noise_mon = BBP700_noise_arr;
chl_noise_mon = chl_noise_arr;
end
