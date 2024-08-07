function [BBP700_noise,chl_noise] = calc_bckg_noise(pres,bbp700,chl,win);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to compute background noise for BGC-ARGO floats

%% Output

% iii. BBP700_noise:  NxM array of BBP700 background noise estimates, where N is the number of pressure levels and M is the number of profiles in bbp700 (m^&-1)
% iv. chl_noise: NxM array of chl-a background noise estimates, where N is the number of pressure levels and M is the number of profiles in chl (mg/m^3)


%% Input

% i.  pres: NxM array of pressure levels, where is the number of pressure levels (same for each profile, but pressure for each profile can vary from j=1,2..M), and M is the number of profiles for pres (dbar)
% ii. bbp700:  NxM array of pressure levels, where is the number of pressure levels, and M is the number of profiles for BBP700 (m^-1)
% iii.chl:  NxM array of pressure levels, where is the number of pressure levels, and M is the number of profiles for chl-a (mg/m^3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global proc_settings

thr = proc_settings.thr;
thr_z = proc_settings.thr_z;
noise_z = proc_settings.noise_z;
noise_dz = proc_settings.noise_dz;

% Arrays are initialized as "0" here so there is something to concatenate

  BBP_resid_bin_tot = 0; % total vector of all binned BBP_resid profiles, needed to compute instrument noise
  Pres_BBP_tmp_bin_tot = 0;

  chl_resid_bin_tot = 0; % total vector of all binned chl profiles, needed to compute instrument noise
  Pres_chl_tmp_bin_tot = 0;


for j=1:length(pres(1,:));

    BBP_tot_tmp = bbp700(:,j);
    Pres_tmp = pres(:,j);
    chl_tot_tmp = chl(:,j);
 
 
% zero out negative chl values. 

    chl_tot_tmp(find(chl_tot_tmp<0)) = 0;

% remove within thr_z m of BBP_tot > thr 
    ind_hi = find(BBP_tot_tmp>thr);
    
    for k = 1:length(BBP_tot_tmp);
      if any(Pres_tmp(k) >Pres_tmp(ind_hi) - thr_z) && any(Pres_tmp(k) <Pres_tmp(ind_hi) + thr_z);
         BBP_tot_tmp(k) = NaN;
         chl_tot_tmp(k) = NaN;
         Pres_tmp(k) = NaN;
      end
    end
% remove NaN values
    indnonan = find(~isnan(BBP_tot_tmp) & ~isnan(chl_tot_tmp));

    Pres_tmp = Pres_tmp(indnonan);
    BBP_tot_tmp = BBP_tot_tmp(indnonan);
    chl_tot_tmp = chl_tot_tmp(indnonan);

% break out of loop if all NaNs
    if sum(~isnan(Pres_tmp)) == 0;
       continue
    end

   
% determine window size for moving minimum and maximum filters

% filter BBP.chl profiles. BBP_fil_tmp should be refractory + small particles
    BBP_fil_tmp = movmin(BBP_tot_tmp,win); % 11 pt moving minimum filter
    BBP_fil_tmp = movmax(BBP_fil_tmp,win); % 11 pt moving maximum filter
    chl_fil_tmp = movmin(chl_tot_tmp,win); % 11 pt moving minimum filter
    chl_fil_tmp = movmax(chl_fil_tmp,win); % 11 pt moving maximum filter
% isolate BBP700 residual

    BBP_resid = BBP_tot_tmp - BBP_fil_tmp; % this is large particels + background of instrument noise, or residual spikes in Briggs et al., 2020
    chl_resid = chl_tot_tmp - chl_fil_tmp; % this is large particels + background of instrument noise, or residual spikes in Briggs et al., 2020
% bin average BBP_resid level noise_z in order to find instrument noise
    if max(Pres_tmp) > noise_z;
      edges = (noise_z:noise_dz:max(Pres_tmp));
      [~,~,loc]=histcounts(Pres_tmp,edges); % number of depths per bin
      Pres_b300 = Pres_tmp;
      BBP_resid_b300 = BBP_resid;
      BBP_resid_b300(find(loc==0)) = []; % remove any bins where there are no pressure values
      Pres_b300(find(loc==0)) = []; % remove any bins where there are no pressure values
      loc(find(loc==0)) = [];% remove any bins where there are no pressure values
   
      Pres_BBP_tmp_bin = accumarray(loc(:),Pres_b300(:),[],@mean);
      BBP_resid_bin = accumarray(loc(:),BBP_resid_b300(:),[],@mean);% depth-binned mean
      BBP_resid_bin_tot = vertcat(BBP_resid_bin_tot,BBP_resid_bin(:));
      Pres_BBP_tmp_bin_tot = vertcat(Pres_BBP_tmp_bin_tot,Pres_BBP_tmp_bin(:));

      [~,~,loc]=histcounts(Pres_tmp,edges); % number of depths per bin
      Pres_b300 = Pres_tmp;
      chl_resid_b300 = chl_resid;
      chl_resid_b300(find(loc==0)) = []; % remove any bins where there are no pressure values
      Pres_b300(find(loc==0)) = []; % remove any bins where there are no pressure values
      loc(find(loc==0)) = [];% remove any bins where there are no pressure values


      Pres_chl_tmp_bin = accumarray(loc(:),Pres_b300(:),[],@mean);
      chl_resid_bin = accumarray(loc(:),chl_resid_b300(:),[],@mean);% depth-binned mean
      chl_resid_bin_tot = vertcat(chl_resid_bin_tot,chl_resid_bin(:));
      Pres_chl_tmp_bin_tot = vertcat(Pres_chl_tmp_bin_tot,Pres_chl_tmp_bin(:));
    end

  end

% calcualte median of each composite vertical profile

  Pres_BBP_tmp_bin_tot = Pres_BBP_tmp_bin_tot(2:end);
  BBP_resid_bin_tot = BBP_resid_bin_tot(2:end);
  chl_resid_bin_tot = chl_resid_bin_tot(2:end);
  median_resid_bbp = median(BBP_resid_bin_tot);
  median_resid_chl = median(chl_resid_bin_tot);
% remove all bins 2x> than the median

  Pres_BBP_tmp_bin_tot = Pres_BBP_tmp_bin_tot(find(BBP_resid_bin_tot<= 2*median_resid_bbp));
  BBP_resid_bin_tot = BBP_resid_bin_tot(find(BBP_resid_bin_tot<= 2*median_resid_bbp));
  Pres_chl_tmp_bin_tot = Pres_chl_tmp_bin_tot(find(chl_resid_bin_tot<= 2*median_resid_chl));
  chl_resid_bin_tot = chl_resid_bin_tot(find(chl_resid_bin_tot<= 2*median_resid_chl));

  BBP700_noise = median(BBP_resid_bin_tot); % Instrument noise for this float
  chl_noise = median(chl_resid_bin_tot); % Instrument noise for this float

end % end function calc_bckg_noise
