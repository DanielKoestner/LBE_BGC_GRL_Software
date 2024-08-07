% Mean profiles
% Seasonal Average profiles in/out (jan-mar, april-june, july-sep,oct-dec)
% MAYBE WE CAN JUST DO A SINGLE LINE FOR EACH MONTH, using colors2, maybe
% try shaded error bar ersion
% POC, T or density, oxygen

% FLOATS ONLY
% code only calculates and saves for plotting later

%% set up


close all
clear all

curdir=cd;

addpath(genpath([curdir '/Scripts']))
addpath(genpath([curdir '/Data']))

addpath(genpath([curdir '/aux']))

save_flag=1;


%% load data
load('LBE_BGC_POC_2010_2022_07-Jun-2024')

%% try rebinning as defined earlier
% Rebin
bins=[0:5:195 200:20:980];
z=out{1}.zbin.z;
mids=[2.5:5:200 210:20:990];

y=discretize(z,bins);

% out
for i = 1:length(out)
    mn_out{i}=month(datetime(out{i}.dnum,'ConvertFrom','datenum'));
    for ii = 1:length(bins)
        inds=ismember(y,ii);
        [~, sz]=size(out{i}.zbin.poc_s);
        if sum(inds)>1
            poc_out{i}(ii,:)=nanmean(out{i}.zbin.poc_s(inds,:));
            pocl_out{i}(ii,:)=nanmean(out{i}.zbin.poc_l(inds,:));
            comp_out{i}(ii,:)=nanmean(out{i}.zbin.chla_s(inds,:)./out{i}.zbin.bbp700_s(inds,:));
            t_out{i}(ii,:)=nanmean(out{i}.zbin.t(inds,:));
            s_out{i}(ii,:)=nanmean(out{i}.zbin.s(inds,:));
            oxy_out{i}(ii,:)=nanmean(out{i}.zbin.oxy(inds,:));
            aou_out{i}(ii,:)=nanmean(out{i}.zbin.aou(inds,:));
            nitrate_out{i}(ii,:)=nanmean(out{i}.zbin.nitrate(inds,:));
        elseif sum(inds)==1
            poc_out{i}(ii,:)=out{i}.zbin.poc_s(inds,:);
            pocl_out{i}(ii,:)=out{i}.zbin.poc_l(inds,:);
            comp_out{i}(ii,:)=out{i}.zbin.chla_s(inds,:)./out{i}.zbin.bbp700_s(inds,:);
            t_out{i}(ii,:)=out{i}.zbin.t(inds,:);
            s_out{i}(ii,:)=out{i}.zbin.s(inds,:);
            oxy_out{i}(ii,:)=out{i}.zbin.oxy(inds,:);
            aou_out{i}(ii,:)=out{i}.zbin.aou(inds,:);
            nitrate_out{i}(ii,:)=out{i}.zbin.nitrate(inds,:);
        else
            poc_out{i}(ii,:)=NaN(sz,1);
            pocl_out{i}(ii,:)=NaN(sz,1);
            comp_out{i}(ii,:)=NaN(sz,1);
            t_out{i}(ii,:)=NaN(sz,1);
            s_out{i}(ii,:)=NaN(sz,1);
            oxy_out{i}(ii,:)=NaN(sz,1);
            aou_out{i}(ii,:)=NaN(sz,1);
            nitrate_out{i}(ii,:)=NaN(sz,1);
        end
    end
end

% in
for i = 1:length(in)
    mn_in{i}=month(datetime(in{i}.dnum,'ConvertFrom','datenum'));

    for ii = 1:length(bins)
        inds=ismember(y,ii);
        [~, sz]=size(in{i}.zbin.poc_s);
        if sum(inds)>1
            poc_in{i}(ii,:)=nanmean(in{i}.zbin.poc_s(inds,:));
            pocl_in{i}(ii,:)=nanmean(in{i}.zbin.poc_l(inds,:));
            comp_in{i}(ii,:)=nanmean(in{i}.zbin.chla_s(inds,:)./in{i}.zbin.bbp700_s(inds,:));
            t_in{i}(ii,:)=nanmean(in{i}.zbin.t(inds,:));
            s_in{i}(ii,:)=nanmean(in{i}.zbin.s(inds,:));
            oxy_in{i}(ii,:)=nanmean(in{i}.zbin.oxy(inds,:));
            aou_in{i}(ii,:)=nanmean(in{i}.zbin.aou(inds,:));
            nitrate_in{i}(ii,:)=nanmean(in{i}.zbin.nitrate(inds,:));
        elseif sum(inds)==1
            poc_in{i}(ii,:)=in{i}.zbin.poc_s(inds,:);
            pocl_in{i}(ii,:)=in{i}.zbin.poc_l(inds,:);
            comp_in{i}(ii,:)=in{i}.zbin.chla_s(inds,:)./in{i}.zbin.bbp700_s(inds,:);
            t_in{i}(ii,:)=in{i}.zbin.t(inds,:);
            s_in{i}(ii,:)=in{i}.zbin.s(inds,:);
            oxy_in{i}(ii,:)=in{i}.zbin.oxy(inds,:);
            aou_in{i}(ii,:)=in{i}.zbin.aou(inds,:);
            nitrate_in{i}(ii,:)=in{i}.zbin.nitrate(inds,:);
        else
            poc_in{i}(ii,:)=NaN(sz,1);
            pocl_in{i}(ii,:)=NaN(sz,1);
            comp_in{i}(ii,:)=NaN(sz,1);
            t_in{i}(ii,:)=NaN(sz,1);
            s_in{i}(ii,:)=NaN(sz,1);
            oxy_in{i}(ii,:)=NaN(sz,1);
            aou_in{i}(ii,:)=NaN(sz,1);
            nitrate_in{i}(ii,:)=NaN(sz,1);
        end
    end
end

%% Combine all data into one large matrix with columns for profiles and a known month
% out
mns_out=[];
pocs_out=[];
pocsl_out=[];
comps_out=[];
ss_out=[];
ts_out=[];
oxys_out=[];
aous_out=[];
nitrates_out=[];
mlds_out=[];
ezds_out=[];
zprods_out=[];
for i = 1:length(mn_out)
    mns_out=[mns_out mn_out{i}];
    mlds_out=[mlds_out out{i}.mld];
    ezds_out=[ezds_out out{i}.ezd];
    zprods_out=[zprods_out out{i}.zprod];

    pocs_out=[pocs_out poc_out{i}];
    pocsl_out=[pocsl_out pocl_out{i}];
    comps_out=[comps_out comp_out{i}];
    ts_out=[ts_out t_out{i}];
    ss_out=[ss_out s_out{i}];
    oxys_out=[oxys_out oxy_out{i}];
    aous_out=[aous_out aou_out{i}];
    nitrates_out=[nitrates_out nitrate_out{i}];
end

% zprod_out=max([mlds_out;ezds_out]);

%months
for i = 1:12
    inds=ismember(mns_out,i);
    POC_out(:,i)=nanmean(pocs_out(:,inds),2);
    POC_out_all{i}=pocs_out(:,inds);
    POC_out_sd(:,i)=nanstd(pocs_out(:,inds)');
    POC_out_q10(:,i)=quantile(pocs_out(:,inds)',0.10);
    POC_out_q90(:,i)=quantile(pocs_out(:,inds)',0.90);

    POCl_out(:,i)=nanmean(pocsl_out(:,inds),2);
    POCl_out_sd(:,i)=nanstd(pocsl_out(:,inds)');
    POCl_out_all{i}=pocsl_out(:,inds);

    COMP_out(:,i)=nanmean(comps_out(:,inds),2);
    COMP_out_sd(:,i)=nanstd(comps_out(:,inds)');
    COMP_out_all{i}=comps_out(:,inds);

    S_out(:,i)=nanmean(ss_out(:,inds),2);
    S_out_sd(:,i)=nanstd(ss_out(:,inds)');
    S_out_all{i}=ss_out(:,inds);

    T_out(:,i)=nanmean(ts_out(:,inds),2);
    T_out_sd(:,i)=nanstd(ts_out(:,inds)');
    T_out_all{i}=ts_out(:,inds);

    OXY_out(:,i)=nanmean(oxys_out(:,inds),2);
    OXY_out_sd(:,i)=nanstd(oxys_out(:,inds)');

    AOU_out(:,i)=nanmean(aous_out(:,inds),2);
    AOU_out_sd(:,i)=nanstd(aous_out(:,inds)');
    AOU_out_all{i}=aous_out(:,inds);

    N_out(:,i)=nanmean(nitrates_out(:,inds),2);
    N_out_sd(:,i)=nanstd(nitrates_out(:,inds)');

    MLD_out(:,i)=nanmean(mlds_out(inds));
    MLD_out_sd(:,i)=nanstd(mlds_out(inds));

    EZD_out(:,i)=nanmean(ezds_out(inds));
    EZD_out_sd(:,i)=nanstd(ezds_out(inds));
    
    zPROD_out(:,i)=nanmean(zprods_out(inds));
    zPROD_out_sd(:,i)=nanstd(zprods_out(inds));
    
end

% in
mns_in=[];
pocs_in=[];
pocsl_in=[];
comps_in=[];
ss_in=[];
ts_in=[];
oxys_in=[];
aous_in=[];
nitrates_in=[];
mlds_in=[];
ezds_in=[];
zprods_in=[];
for i = 1:length(mn_in)
    mns_in=[mns_in mn_in{i}];
    mlds_in=[mlds_in in{i}.mld];
    ezds_in=[ezds_in in{i}.ezd];
    zprods_in=[zprods_in in{i}.zprod];

    pocs_in=[pocs_in poc_in{i}];
    pocsl_in=[pocsl_in pocl_in{i}];
    comps_in=[comps_in comp_in{i}];
    ss_in=[ss_in s_in{i}];
    ts_in=[ts_in t_in{i}];
    oxys_in=[oxys_in oxy_in{i}];
    aous_in=[aous_in aou_in{i}];
    nitrates_in=[nitrates_in nitrate_in{i}];
end

% zprod_in=max([mlds_in;ezds_in]);
%months
for i = 1:12
    inds=ismember(mns_in,i);
    POC_in(:,i)=nanmean(pocs_in(:,inds),2);
    POC_in_all{i}=pocs_in(:,inds);
    POC_in_sd(:,i)=nanstd(pocs_in(:,inds)');
    POC_in_q10(:,i)=quantile(pocs_in(:,inds)',0.10);
    POC_in_q90(:,i)=quantile(pocs_in(:,inds)',0.90);

    POCl_in(:,i)=nanmean(pocsl_in(:,inds),2);
    POCl_in_sd(:,i)=nanstd(pocsl_in(:,inds)');
    POCl_in_all{i}=pocsl_in(:,inds);

    COMP_in(:,i)=nanmean(comps_in(:,inds),2);
    COMP_in_sd(:,i)=nanstd(comps_in(:,inds)');
    COMP_in_all{i}=comps_in(:,inds);

    S_in(:,i)=nanmean(ss_in(:,inds),2);
    S_in_sd(:,i)=nanstd(ss_in(:,inds)');
    S_in_all{i}=ss_in(:,inds);

    T_in(:,i)=nanmean(ts_in(:,inds),2);
    T_in_sd(:,i)=nanstd(ts_in(:,inds)');
    T_in_all{i}=ts_in(:,inds);

    OXY_in(:,i)=nanmean(oxys_in(:,inds),2);
    OXY_in_sd(:,i)=nanstd(oxys_in(:,inds)');

    AOU_in(:,i)=nanmean(aous_in(:,inds),2);
    AOU_in_sd(:,i)=nanstd(aous_in(:,inds)');
    AOU_in_all{i}=aous_in(:,inds);

    N_in(:,i)=nanmean(nitrates_in(:,inds),2);
    N_in_sd(:,i)=nanstd(nitrates_in(:,inds)');

    MLD_in(:,i)=nanmean(mlds_in(inds));
    MLD_in_sd(:,i)=nanstd(mlds_in(inds));

    EZD_in(:,i)=nanmean(ezds_in(inds));
    EZD_in_sd(:,i)=nanstd(ezds_in(inds));

    zPROD_in(:,i)=nanmean(zprods_in(inds));
    zPROD_in_sd(:,i)=nanstd(zprods_in(inds));
end



%% save
if save_flag==1
    Z=mids;
    save(['Data/LBE_BGC_POC_2010_2022_MonthlyProfiles_' date],'POC_in','POC_in_sd','POCl_in','POCl_in_sd','AOU_in','AOU_in_sd','zPROD_in','zPROD_in_sd','POC_out','POC_out_sd','POCl_out','POCl_out_sd','AOU_out','AOU_out_sd','zPROD_out','zPROD_out_sd','Z','AOU_in_all','AOU_out_all','T_in','T_out','T_in_sd','T_out_sd','T_in_all','T_out_all','MLD_in','MLD_in_sd','MLD_out','MLD_out_sd','EZD_in','EZD_in_sd','EZD_out','EZD_out_sd','POC_in_all','POC_out_all','POCl_in_all','POCl_out_all','S_in','S_out','S_in_sd','S_out_sd','S_in_all','S_out_all','COMP_in','COMP_out','COMP_in_sd','COMP_out_sd','COMP_in_all','COMP_out_all')
end