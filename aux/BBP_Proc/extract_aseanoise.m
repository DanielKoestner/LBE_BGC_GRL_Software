function [BBP700_noise_asea,chl_noise_asea] = extract_aseanoise(time,years,bbp700,chl,pres,win);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function return BBP700 and chl annually varying background noise, computed using only data for a specific season (Dec-Feb,Mar-May,June-Aug, or Sep-Nov).

%% Output

% iii. BBP700_noise_asea:  NxM array of BBP700 background noise estimates, where N is the number of pressure levels and M is the number of profiles in bbp700 (m^-1)
% iv.  chl_noise_asea: NxM array of chl-a background noise estimates, where N is the number of pressure levels and M is the number of profiles in chl (mg/m^3)


%% Input

% i.    time: a 1xM datetime array (YYYY-MM-DD), which contains the dates during which profiles are sampled
% ii.   years: 1xMr array, where Mr is the number of non-repeating years of the float's lifetime
% iii.  pres: NxM array of pressure levels, where is the number of pressure levels (same for each profile, but pressure for each profile can vary from j=1,2..M), and M is the number of profiles for pres (dbar)
% iv.   bbp700:  NxM array of pressure levels, where is the number of pressure levels, and M is the number of profiles for BBP700 (m^-1)
% v.    chl:  NxM array of pressure levels, where is the number of pressure levels, and M is the number of profiles for chl-a (mg/m^3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global proc_settings

months = proc_settings.months; % read in months from global array

% preallocate bbp700 and chl noise arrays

BBP700_noise_asea = NaN.*ones(size(bbp700));
chl_noise_asea = NaN.*ones(size(chl));

year_flag = NaN.*ones(length(years),1); % flag for if data exists or not
indwin = arrayfun(@(x) find(month(time)==x), months,'un',0); % find matching months in time
indwin = [indwin{:}]; % unpack indices

timevecwin = time(indwin); % find time corresponding to months in season
yearvecwin = year(timevecwin); % find years corresponding to months in season
if any(months == 12); % treat December of year Y-1 as part of winter of year Y. Only applicable if wintertime data are used for noise computation
  yearvecwin(month(timevecwin) == 12) = yearvecwin(month(timevecwin) == 12) + 1; 
end
yearvecwin_un = unique(yearvecwin); % find unique years with seasonal data to loop over
yearwin_flag = NaN.*ones(length(yearvecwin_un),1); % flag for determinig if seasonal data exists

% If the float does not contain data for the chosen season, compute the annually varying noise instead. This treats the noise computation as though the 'annual' option is used
if isempty(indwin);
  fprintf('no wintertime data found. Processing noise anually using all data for each year\n')
  [BBP700_noise_ann,chl_noise_ann] = extract_annoise(time,years,bbp700,chl,pres);
  BBP700_noise_asea = BBP700_noise_ann;
  chl_noise_asea = chl_noise_ann;
else
  yrs_nowin = 0; % for concatenating yrs without seasonal data
  inds_retr = 0; % for concatenating indices of yrs without seasonal data
  BBP700_noise_all = [];
  chl_noise_all = [];
  for j = 1:length(yearvecwin_un);
    indwintmp = indwin(find(yearvecwin == yearvecwin_un(j))); % extract only seasons of current year
% extract cycle of variables corresponding to season of current year
    BBP700_win = bbp700(:,indwintmp); 
    chl_win = chl(:,indwintmp);
    Pres_win = pres(:,indwintmp);
     
% compute noise for extracted arrays
    if sum(~isnan(BBP700_win(:)))>0 && sum(~isnan(chl_win(:)))>0; % only compute noise if there is data for this season of this year
      [BBP700_noise,chl_noise] = calc_bckg_noise(Pres_win,BBP700_win,chl_win,win);
      BBP700_noise_all(j) = BBP700_noise; % array to keep track of the noise value computed for season of year (j)
      chl_noise_all(j) = chl_noise; % array to keep track of the noise value computed for season of year (j)
      yearwin_flag(j) = 1; % if there is data, set flag for that season of that year to true (1)
      if any(years == yearvecwin_un(j));
        year_flag(find(years == yearvecwin_un(j))) = 1; % Also if the float was active during the year within which the season is contained, set the flag for that year to true (1)
      end
 
    else
% set to NaN if the data to compute the noise for season of year (j) is not available
      BBP700_noise_all(j) = NaN;
      chl_noise_all(j) = NaN; 
     end
  end % end over j=yearvecwin_un

% if we cannot compute the noise for any season of any year, then use the 'annual' method
  if all(isnan(BBP700_noise_all));
    fprintf('no wintertime data found. Processing noise anually using all data for each year\n')
    [BBP700_noise_ann,chl_noise_ann] = extract_annoise(time,years,bbp700,chl,pres);
    BBP700_noise_asea = BBP700_noise_ann;
    chl_noise_asea = chl_noise_ann;
  else
    for k =1:length(years);
      if ~isnan(year_flag(k)); 
        indyr = find(year(time) == years(k)); % find all indices corresponding to the current year
        indyr_noise = find(years(k) == yearvecwin_un); % find all indices of unique seasons corresponding to the year (k)
        BBP700_noise_asea(:,indyr) = BBP700_noise_all(indyr_noise);
        chl_noise_asea(:,indyr) = chl_noise_all(indyr_noise);

      else

% keep track of years without any seasonal float data: these will take the seasonal noise computed from the closest available year. If two equidistant years are found, the subsequent year is chosen
        yearvecwin_data = yearvecwin_un(~isnan(yearwin_flag)); % only take years with seasonal data
        ind_retr = find(abs(years(k) - yearvecwin_data) == min(abs(years(k) - yearvecwin_data))); % find the seasonal data-containing year that is closest to year with no seasonal data
        ind_retr = max(ind_retr); % if two years are the same distance from year with no seasonal data, take the subsequent year 

% These arrays of for looping over yrs without data 
        inds_retr = vertcat(inds_retr,ind_retr); 
        yrs_nowin = vertcat(yrs_nowin,years(k)); 
      end
    end % end over k=years
    yrs_nowin = yrs_nowin(2:end);
    inds_retr = inds_retr(2:end);

% Now loop over yrs without seasonal data and take noise from indices corresponding to closest available seasonal data
    for k = 1:length(yrs_nowin);
      indyr = find(year(time) == yrs_nowin(k)); % find all instances of the datetime array corresponding to the current year
      ind_retr = inds_retr(k); % find corresponding index of the closets season with an avaialbe noise estiamte
      indyr_noise = (find(yearvecwin_un == yearvecwin_data(ind_retr))); % find the index for when the full seasonal array (yrs where seasonal data is available) equals the array containing seasons of real noise estiamtes
      BBP700_noise_asea(:,indyr)  = BBP700_noise_all(indyr_noise); % set bbp700 noise arr for profiles at the current year equal to the closest available seasonal noise estimate
      chl_noise_asea(:,indyr) = chl_noise_all(indyr_noise);  % set chl-a noise arr for profiles at the current year equal to the closest available seasonal noise estimate
    end 
  end % end over if all BBP700 noise is NaN
end  % end over there is no seaosnal data

end % end function
