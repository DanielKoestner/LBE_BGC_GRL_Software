% Mean profiles
% Seasonal Average profiles in/out (jan-mar, april-june, july-sep,oct-dec)
% MAYBE WE CAN JUST DO A SINGLE LINE FOR EACH MONTH, using colors2, maybe
% try shaded error bar ersion
% POC, T or density, oxygen

% FLOATS ONLY
% code only calculates and saves for plotting later


% add chla

%fix storing of EZD as (maybe) max of PAR and Chla metric

% 22 jan 2025, add removal of profiles with threshold measurements and
% mld>900
%% set up


close all
clear all

curdir=cd;

addpath(genpath([curdir '/Scripts']))
addpath(genpath([curdir '/Data']))

addpath(genpath([curdir '/aux']))

save_flag=1;
scale_flag='A'; % A for 1.5r and B for 2r scaling of LBEZ

thrs=0.702; %exclude profiles with less than 60 5-m bin resolved measures
%% load data
load(['LBE_BGC_POC_2010_2022_' scale_flag '_06-Feb-2025'])

%% try rebinning as defined earlier
% Rebin
bins=[0:5:195 200:20:1000];
z=out{1}.zbin.z;
mids=[2.5:5:200 210:20:990];

y=discretize(z,bins);

% out
for i = 1:length(out)
    inds2=sum(isnan(out{i}.zbin.poc_s))<(thrs)*length(z);
    mn_out{i}=month(datetime(out{i}.dnum(inds2),'ConvertFrom','datenum'));
    mld_out{i}=out{i}.mld(inds2);
    ezd_out{i}=out{i}.ezd(inds2);
    zprod_out{i}=out{i}.zprod(inds2);

    for ii = 1:length(mids)
        inds=ismember(y,ii);
        [~, sz]=size(out{i}.zbin.poc_s(:,inds2));
        if sum(inds)>1
            poc_out{i}(ii,:)=nanmean(out{i}.zbin.poc_s(inds,inds2));
            pocl_out{i}(ii,:)=nanmean(out{i}.zbin.poc_l(inds,inds2));
            bbp_out{i}(ii,:)=nanmean(out{i}.zbin.bbp700_s(inds,inds2));
            chla_out{i}(ii,:)=nanmean(out{i}.zbin.chla_s(inds,inds2)+out{i}.zbin.chla_l(inds,inds2));
            t_out{i}(ii,:)=nanmean(out{i}.zbin.t(inds,inds2));
            s_out{i}(ii,:)=nanmean(out{i}.zbin.s(inds,inds2));
            oxy_out{i}(ii,:)=nanmean(out{i}.zbin.oxy(inds,inds2));
            aou_out{i}(ii,:)=nanmean(out{i}.zbin.aou(inds,inds2));
            nitrate_out{i}(ii,:)=nanmean(out{i}.zbin.nitrate(inds,inds2));
        elseif sum(inds)==1
            poc_out{i}(ii,:)=out{i}.zbin.poc_s(inds,inds2);
            pocl_out{i}(ii,:)=out{i}.zbin.poc_l(inds,inds2);
            bbp_out{i}(ii,:)=out{i}.zbin.bbp700_s(inds,inds2);
            chla_out{i}(ii,:)=out{i}.zbin.chla_s(inds,inds2)+out{i}.zbin.chla_l(inds,inds2);
            t_out{i}(ii,:)=out{i}.zbin.t(inds,inds2);
            s_out{i}(ii,:)=out{i}.zbin.s(inds,inds2);
            oxy_out{i}(ii,:)=out{i}.zbin.oxy(inds,inds2);
            aou_out{i}(ii,:)=out{i}.zbin.aou(inds,inds2);
            nitrate_out{i}(ii,:)=out{i}.zbin.nitrate(inds,inds2);
        else
            poc_out{i}(ii,:)=NaN(sz,1);
            pocl_out{i}(ii,:)=NaN(sz,1);
            bbp_out{i}(ii,:)=NaN(sz,1);
            chla_out{i}(ii,:)=NaN(sz,1);
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
    inds2=sum(isnan(in{i}.zbin.poc_s))<(thrs)*length(z);
    mn_in{i}=month(datetime(in{i}.dnum(inds2),'ConvertFrom','datenum'));
    mld_in{i}=in{i}.mld(inds2);
    ezd_in{i}=in{i}.ezd(inds2);
    zprod_in{i}=in{i}.zprod(inds2);

    for ii = 1:length(mids)
        inds=ismember(y,ii);
        [~, sz]=size(in{i}.zbin.poc_s(:,inds2));

        if sum(inds)>1
            poc_in{i}(ii,:)=nanmean(in{i}.zbin.poc_s(inds,inds2));
            pocl_in{i}(ii,:)=nanmean(in{i}.zbin.poc_l(inds,inds2));
            bbp_in{i}(ii,:)=nanmean(in{i}.zbin.bbp700_s(inds,inds2));
            chla_in{i}(ii,:)=nanmean(in{i}.zbin.chla_s(inds,inds2)+in{i}.zbin.chla_l(inds,inds2));
            t_in{i}(ii,:)=nanmean(in{i}.zbin.t(inds,inds2));
            s_in{i}(ii,:)=nanmean(in{i}.zbin.s(inds,inds2));
            oxy_in{i}(ii,:)=nanmean(in{i}.zbin.oxy(inds,inds2));
            aou_in{i}(ii,:)=nanmean(in{i}.zbin.aou(inds,inds2));
            nitrate_in{i}(ii,:)=nanmean(in{i}.zbin.nitrate(inds,inds2));
        elseif sum(inds)==1
            poc_in{i}(ii,:)=in{i}.zbin.poc_s(inds,inds2);
            pocl_in{i}(ii,:)=in{i}.zbin.poc_l(inds,inds2);
            bbp_in{i}(ii,:)=in{i}.zbin.bbp700_s(inds,inds2);
            chla_in{i}(ii,:)=in{i}.zbin.chla_s(inds,inds2)+in{i}.zbin.chla_l(inds,inds2);
            t_in{i}(ii,:)=in{i}.zbin.t(inds,inds2);
            s_in{i}(ii,:)=in{i}.zbin.s(inds,inds2);
            oxy_in{i}(ii,:)=in{i}.zbin.oxy(inds,inds2);
            aou_in{i}(ii,:)=in{i}.zbin.aou(inds,inds2);
            nitrate_in{i}(ii,:)=in{i}.zbin.nitrate(inds,inds2);
        else
            poc_in{i}(ii,:)=NaN(sz,1);
            pocl_in{i}(ii,:)=NaN(sz,1);
            bbp_in{i}(ii,:)=NaN(sz,1);
            chla_in{i}(ii,:)=NaN(sz,1);
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
bbps_out=[];
chlas_out=[];
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
    mlds_out=[mlds_out mld_out{i}];
    ezds_out=[ezds_out ezd_out{i}];
    zprods_out=[zprods_out zprod_out{i}];

    pocs_out=[pocs_out poc_out{i}];
    pocsl_out=[pocsl_out pocl_out{i}];
    bbps_out=[bbps_out bbp_out{i}];
    chlas_out=[chlas_out chla_out{i}];
    ts_out=[ts_out t_out{i}];
    ss_out=[ss_out s_out{i}];
    oxys_out=[oxys_out oxy_out{i}];
    aous_out=[aous_out aou_out{i}];
    nitrates_out=[nitrates_out nitrate_out{i}];
end

% zprod_out=max([mlds_out;ezds_out]);
% mlds_out(mlds_out>900)=NaN;
%months
for i = 1:12
%     inds=ismember(mns_out,i);
    inds=ismember(mns_out,i)&zprods_out<900; %exclude profiles with MLD>900 m or no MLD
    POC_out(:,i)=nanmean(pocs_out(:,inds),2);
    POC_out_all{i}=pocs_out(:,inds);
    POC_out_sd(:,i)=nanstd(pocs_out(:,inds)');
    POC_out_q10(:,i)=quantile(pocs_out(:,inds)',0.10);
    POC_out_q90(:,i)=quantile(pocs_out(:,inds)',0.90);
    [~,~,ci,~]=ttest(pocs_out(:,inds)');
    POC_out_low(:,i)=ci(1,:);
    POC_out_hi(:,i)=ci(2,:);
    clear ci

    POCl_out(:,i)=nanmean(pocsl_out(:,inds),2);
    POCl_out_sd(:,i)=nanstd(pocsl_out(:,inds)');
    POCl_out_all{i}=pocsl_out(:,inds);
    [~,~,ci,~]=ttest(pocsl_out(:,inds)');
    POCl_out_low(:,i)=ci(1,:);
    POCl_out_hi(:,i)=ci(2,:);
    clear ci

    BBP_out(:,i)=nanmean(bbps_out(:,inds),2);
    BBP_out_sd(:,i)=nanstd(bbps_out(:,inds)');
    BBP_out_all{i}=bbps_out(:,inds);
    [~,~,ci,~]=ttest(bbps_out(:,inds)');
    BBP_out_low(:,i)=ci(1,:);
    BBP_out_hi(:,i)=ci(2,:);
    clear ci

    CHLA_out(:,i)=nanmean(chlas_out(:,inds),2);
    CHLA_out_sd(:,i)=nanstd(chlas_out(:,inds)');
    CHLA_out_all{i}=chlas_out(:,inds);
    [~,~,ci,~]=ttest(chlas_out(:,inds)');
    CHLA_out_low(:,i)=ci(1,:);
    CHLA_out_hi(:,i)=ci(2,:);
    clear ci

    S_out(:,i)=nanmean(ss_out(:,inds),2);
    S_out_sd(:,i)=nanstd(ss_out(:,inds)');
    S_out_all{i}=ss_out(:,inds);
    [~,~,ci,~]=ttest(ss_out(:,inds)');
    S_out_low(:,i)=ci(1,:);
    S_out_hi(:,i)=ci(2,:);
    clear ci

    T_out(:,i)=nanmean(ts_out(:,inds),2);
    T_out_sd(:,i)=nanstd(ts_out(:,inds)');
    T_out_all{i}=ts_out(:,inds);
    [~,~,ci,~]=ttest(ts_out(:,inds)');
    T_out_low(:,i)=ci(1,:);
    T_out_hi(:,i)=ci(2,:);
    clear ci

    OXY_out(:,i)=nanmean(oxys_out(:,inds),2);
    OXY_out_sd(:,i)=nanstd(oxys_out(:,inds)');
    [~,~,ci,~]=ttest(oxys_out(:,inds)');
    OXY_out_low(:,i)=ci(1,:);
    OXY_out_hi(:,i)=ci(2,:);
    clear ci

    AOU_out(:,i)=nanmean(aous_out(:,inds),2);
    AOU_out_sd(:,i)=nanstd(aous_out(:,inds)');
    AOU_out_all{i}=aous_out(:,inds);
    [~,~,ci,~]=ttest(aous_out(:,inds)');
    AOU_out_low(:,i)=ci(1,:);
    AOU_out_hi(:,i)=ci(2,:);
    clear ci

    N_out(:,i)=nanmean(nitrates_out(:,inds),2);
    N_out_sd(:,i)=nanstd(nitrates_out(:,inds)');
    [~,~,ci,~]=ttest(nitrates_out(:,inds)');
    N_out_low(:,i)=ci(1,:);
    N_out_hi(:,i)=ci(2,:);
    clear ci

    MLD_out(:,i)=nanmean(mlds_out(inds));
    MLD_out_sd(:,i)=nanstd(mlds_out(inds));
    [~,~,ci,~]=ttest(mlds_out(:,inds)');
    MLD_out_low(:,i)=ci(1,:);
    MLD_out_hi(:,i)=ci(2,:);
    clear ci

    EZD_out(:,i)=nanmean(ezds_out(inds));
    EZD_out_sd(:,i)=nanstd(ezds_out(inds));
    [~,~,ci,~]=ttest(ezds_out(:,inds)');
    EZD_out_low(:,i)=ci(1,:);
    EZD_out_hi(:,i)=ci(2,:);
    clear ci

    zPROD_out(:,i)=nanmean(zprods_out(inds));
    zPROD_out_sd(:,i)=nanstd(zprods_out(inds));
    [~,~,ci,~]=ttest(zprods_out(:,inds)');
    zPROD_out_low(:,i)=ci(1,:);
    zPROD_out_hi(:,i)=ci(2,:);
    clear ci
end

% in
mns_in=[];
pocs_in=[];
pocsl_in=[];
bbps_in=[];
chlas_in=[];
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
    mlds_in=[mlds_in mld_in{i}];
    ezds_in=[ezds_in ezd_in{i}];
    zprods_in=[zprods_in zprod_in{i}];

    pocs_in=[pocs_in poc_in{i}];
    pocsl_in=[pocsl_in pocl_in{i}];
    bbps_in=[bbps_in bbp_in{i}];
    chlas_in=[chlas_in chla_in{i}];
    ss_in=[ss_in s_in{i}];
    ts_in=[ts_in t_in{i}];
    oxys_in=[oxys_in oxy_in{i}];
    aous_in=[aous_in aou_in{i}];
    nitrates_in=[nitrates_in nitrate_in{i}];
end

% zprod_in=max([mlds_in;ezds_in]);
%months
% mlds_in(mlds_in>900)=NaN;
for i = 1:12
%     inds=ismember(mns_in,i);
    inds=ismember(mns_in,i)&zprods_in<900; %exclude profiles with MLD>900 m or no MLD
    POC_in(:,i)=nanmean(pocs_in(:,inds),2);
    POC_in_all{i}=pocs_in(:,inds);
    POC_in_sd(:,i)=nanstd(pocs_in(:,inds)');
    POC_in_q10(:,i)=quantile(pocs_in(:,inds)',0.10);
    POC_in_q90(:,i)=quantile(pocs_in(:,inds)',0.90);
    [~,~,ci,~]=ttest(pocs_in(:,inds)');
    POC_in_low(:,i)=ci(1,:);
    POC_in_hi(:,i)=ci(2,:);
    clear ci

    POCl_in(:,i)=nanmean(pocsl_in(:,inds),2);
    POCl_in_sd(:,i)=nanstd(pocsl_in(:,inds)');
    POCl_in_all{i}=pocsl_in(:,inds);
    [~,~,ci,~]=ttest(pocsl_in(:,inds)');
    POCl_in_low(:,i)=ci(1,:);
    POCl_in_hi(:,i)=ci(2,:);
    clear ci

    BBP_in(:,i)=nanmean(bbps_in(:,inds),2);
    BBP_in_sd(:,i)=nanstd(bbps_in(:,inds)');
    BBP_in_all{i}=bbps_in(:,inds);
    [~,~,ci,~]=ttest(bbps_in(:,inds)');
    BBP_in_low(:,i)=ci(1,:);
    BBP_in_hi(:,i)=ci(2,:);
    clear ci

    CHLA_in(:,i)=nanmean(chlas_in(:,inds),2);
    CHLA_in_sd(:,i)=nanstd(chlas_in(:,inds)');
    CHLA_in_all{i}=chlas_in(:,inds);
    [~,~,ci,~]=ttest(chlas_in(:,inds)');
    CHLA_in_low(:,i)=ci(1,:);
    CHLA_in_hi(:,i)=ci(2,:);
    clear ci

    S_in(:,i)=nanmean(ss_in(:,inds),2);
    S_in_sd(:,i)=nanstd(ss_in(:,inds)');
    S_in_all{i}=ss_in(:,inds);
    [~,~,ci,~]=ttest(ss_in(:,inds)');
    S_in_low(:,i)=ci(1,:);
    S_in_hi(:,i)=ci(2,:);
    clear ci

    T_in(:,i)=nanmean(ts_in(:,inds),2);
    T_in_sd(:,i)=nanstd(ts_in(:,inds)');
    T_in_all{i}=ts_in(:,inds);
    [~,~,ci,~]=ttest(ts_in(:,inds)');
    T_in_low(:,i)=ci(1,:);
    T_in_hi(:,i)=ci(2,:);
    clear ci

    OXY_in(:,i)=nanmean(oxys_in(:,inds),2);
    OXY_in_sd(:,i)=nanstd(oxys_in(:,inds)');
    [~,~,ci,~]=ttest(oxys_in(:,inds)');
    OXY_in_low(:,i)=ci(1,:);
    OXY_in_hi(:,i)=ci(2,:);
    clear ci

    AOU_in(:,i)=nanmean(aous_in(:,inds),2);
    AOU_in_sd(:,i)=nanstd(aous_in(:,inds)');
    AOU_in_all{i}=aous_in(:,inds);
    [~,~,ci,~]=ttest(aous_in(:,inds)');
    AOU_in_low(:,i)=ci(1,:);
    AOU_in_hi(:,i)=ci(2,:);
    clear ci

    N_in(:,i)=nanmean(nitrates_in(:,inds),2);
    N_in_sd(:,i)=nanstd(nitrates_in(:,inds)');
    [~,~,ci,~]=ttest(nitrates_in(:,inds)');
    N_in_low(:,i)=ci(1,:);
    N_in_hi(:,i)=ci(2,:);
    clear ci

    MLD_in(:,i)=nanmean(mlds_in(inds));
    MLD_in_sd(:,i)=nanstd(mlds_in(inds));
    [~,~,ci,~]=ttest(mlds_in(:,inds)');
    MLD_in_low(:,i)=ci(1,:);
    MLD_in_hi(:,i)=ci(2,:);
    clear ci

    EZD_in(:,i)=nanmean(ezds_in(inds));
    EZD_in_sd(:,i)=nanstd(ezds_in(inds));
    [~,~,ci,~]=ttest(ezds_in(:,inds)');
    EZD_in_low(:,i)=ci(1,:);
    EZD_in_hi(:,i)=ci(2,:);
    clear ci

    zPROD_in(:,i)=nanmean(zprods_in(inds));
    zPROD_in_sd(:,i)=nanstd(zprods_in(inds));
    [~,~,ci,~]=ttest(zprods_in(:,inds)');
    zPROD_in_low(:,i)=ci(1,:);
    zPROD_in_hi(:,i)=ci(2,:);
    clear ci
end



%% save
if save_flag==1
    Z=mids;
    if sf==1.5
        save(['Data/LBE_BGC_POC_2010_2022_MonthlyProfiles_A_' date],'sf','POC*','CHLA*','MLD*','N*','OXY*','AOU*','zPROD*','EZD*','T*','S*','BBP*','Z')
    elseif sf==2
         save(['Data/LBE_BGC_POC_2010_2022_MonthlyProfiles_B_' date],'sf','POC*','CHLA*','MLD*','N*','OXY*','AOU*','zPROD*','EZD*','T*','S*','BBP*','Z')
       
%         save(['Data/LBE_BGC_POC_2010_2022_MonthlyProfiles_B_' date],'sf','POC_in','POC_in_sd','POCl_in','POCl_in_sd','AOU_in','AOU_in_sd','zPROD_in','zPROD_in_sd','POC_out','POC_out_sd','POCl_out','POCl_out_sd','AOU_out','AOU_out_sd','zPROD_out','zPROD_out_sd','Z','AOU_in_all','AOU_out_all','T_in','T_out','T_in_sd','T_out_sd','T_in_all','T_out_all','MLD_in','MLD_in_sd','MLD_out','MLD_out_sd','EZD_in','EZD_in_sd','EZD_out','EZD_out_sd','POC_in_all','POC_out_all','POCl_in_all','POCl_out_all','S_in','S_out','S_in_sd','S_out_sd','S_in_all','S_out_all','COMP_in','COMP_out','COMP_in_sd','COMP_out_sd','COMP_in_all','COMP_out_all','CHLA_in','CHLA_out','CHLA_in_sd','CHLA_out_sd','CHLA_in_all','CHLA_out_all')
    end
end