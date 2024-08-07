function [BBP700_noise_ann,chl_noise_ann] = extract_annoise(time,years,bbp700,chl,pres,win);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to extract annually varying background noise for BGC-ARGO floats

%% Output

% iii. BBP700_noise_ann:  NxM array of BBP700 background noise estimates, where N is the number of pressure levels and M is the number of profiles in bbp700 (m^-1)
% iv.  chl_noise_ann: NxM array of chl-a background noise estimates, where N is the number of pressure levels and M is the number of profiles in chl (mg/m^3)


%% Input

% i.    time: a 1xM datetime array (YYYY-MM-DD), which contains the dates during which profiles are sampled
% ii.   years: 1xMr array, where Mr is the number of non-repeating years of the float's lifetime
% iii.  pres: NxM array of pressure levels, where is the number of pressure levels (same for each profile, but pressure for each profile can vary from j=1,2..M), and M is the number of profiles for pres (dbar)
% iv.   bbp700:  NxM array of pressure levels, where is the number of pressure levels, and M is the number of profiles for BBP700 (m^-1)
% v.    chl:  NxM array of pressure levels, where is the number of pressure levels, and M is the number of profiles for chl-a (mg/m^3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global proc_settings

BBP700_noise_arr = NaN.*ones(size(bbp700));
chl_noise_arr = NaN.*ones(size(chl));

year_flag = NaN.*ones(length(years),1); % flag for if data exists or not
% if there is no data, remove output
if sum(~isnan(bbp700(:))) == 0 || sum(~isnan(chl(:))) == 0;
  fprintf('Not enough chla or BBP700 values found to compute noise for this float. Noise Values set to NaN \n')
else
% first find years with enough data to obtain noise estimate 
  for k =1:length(years);
    indyr = (find(year(time) == years(k))); % extract all occurances current year
    bbp700_yr = bbp700(:,indyr); % subsample bbp700 at current year
    chl_yr = chl(:,indyr); % subsample chl-a at current year
    pres_yr = pres(:,indyr); % subsample pressure at current year

% compute background noise for the current year
    [BBP700_noise,chl_noise] = calc_bckg_noise(pres_yr,bbp700_yr,chl_yr,win);

% fill in noise arrays if noise was computed    
    if ~isnan(BBP700_noise(:))>0 && ~isnan(chl_noise(:))>0
      BBP700_noise_arr(:,indyr) = BBP700_noise;
      chl_noise_arr(:,indyr) = chl_noise;
      BBP700_noise_all(k) = BBP700_noise;
      chl_noise_all(k) = chl_noise;
      year_flag(k) = 1; % set flag for current year to true (1) if noise estimates are available
    else

% otherwise set noise for current year to NaN
      BBP700_noise_all(k) = NaN;
      chl_noise_all(k) = NaN;
    end
  end

  if all(isnan(BBP700_noise_all));
    fprintf('Not enough chla or BBP700 values found to compute noise for this float. Noise values set to NaN. \n') 
  else
    year_data = years(find(~isnan(year_flag))); % find years where noise estimates exist

% fill out remaining years with noise estimates of closest year
    for k =1:length(years);
      if isnan(year_flag(k));
        ind_data = find(abs(years(k) - year_data) == min(abs(years(k) - year_data))); % find index of closest year with available noise estimate
        idx=randperm(length(ind_data),1);
        ind_data = ind_data(idx); % if two such years are available, pick one at random
        indyr_noise = (find(years == year_data(ind_data))); % find the index for when the full year array equals the array containing years of real noise estiamtes
        indyr = (find(year(time) == years(k))); % find all instances of the datetime array corresponding to the current year
        BBP700_noise_arr(:,indyr) = BBP700_noise_all(indyr_noise); % set bbp700 noise arr for profiles at the current year equal to the noise estimate available from the closets year
        chl_noise_arr(:,indyr) = chl_noise_all(indyr_noise); % set chl noise arr for profiles at the current year equal to the noise estimate available from the closets year
      end
    end

  end % if all BBP700_noise values are NaN
end %if all data are NaN


BBP700_noise_ann = BBP700_noise_arr;
chl_noise_ann = chl_noise_arr;
end
