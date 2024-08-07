% script to sort data and remake an "lb" and "lbe" clear all

% copy lb and lbe, float index by float index, then reset data for each
% eddy_id

clear all
close all

save_flag=1;

% load('LBE_BGC_POC_2010_2022_07-Jun-2024.mat')
load('LBE_BGC_POC_2010_2022_22-Jul-2024.mat') %w/ vmr

out=lb;
in=lb;
inid=[];

for flt=1:length(lb)

    for i = 1:length(lb{flt}.dnum)

        argo_lat=lb{flt}.lat(i);
        argo_lon=lb{flt}.lon(i);
        argo_date=lb{flt}.dnum(i);

        [eddy_id(i)]=in_eddy(argo_date,argo_lat,argo_lon,'LBE_locations.mat');
        clear argo*
    end


    fields=fieldnames(lb{flt});

    % write only outside eddy data
    out{flt}.bbp700=lb{flt}.bbp700(:,eddy_id==0);
    out{flt}.bbp700_s=lb{flt}.bbp700_s(:,eddy_id==0);
    out{flt}.bbp700_l=lb{flt}.bbp700_l(:,eddy_id==0);
    if sum(ismember('par',fields))>0
        out{flt}.par=lb{flt}.par(:,eddy_id==0);
        out{flt}.ezd=lb{flt}.ezd(:,eddy_id==0);
        out{flt}.par0=lb{flt}.par0(:,eddy_id==0);
    end
    out{flt}.t=lb{flt}.t(:,eddy_id==0);
    out{flt}.s=lb{flt}.s(:,eddy_id==0);
    out{flt}.mld=lb{flt}.mld(:,eddy_id==0);
    out{flt}.ezd_c=lb{flt}.ezd_c(:,eddy_id==0);
    out{flt}.zprod=lb{flt}.zprod(:,eddy_id==0);
    out{flt}.vmr_100=lb{flt}.vmr_100(:,eddy_id==0);
    out{flt}.chla=lb{flt}.chla(:,eddy_id==0);
    out{flt}.chla_s=lb{flt}.chla_s(:,eddy_id==0);
    out{flt}.chla_l=lb{flt}.chla_l(:,eddy_id==0);
    out{flt}.pres=lb{flt}.pres(:,eddy_id==0);
    out{flt}.lat=lb{flt}.lat(:,eddy_id==0);
    out{flt}.lon=lb{flt}.lon(:,eddy_id==0);
    out{flt}.dnum=lb{flt}.dnum(:,eddy_id==0);
    if sum(ismember('oxy',fields))>0
        out{flt}.oxy=lb{flt}.oxy(:,eddy_id==0);
    end
    if sum(ismember('nitrate',fields))>0
        out{flt}.nitrate=lb{flt}.nitrate(:,eddy_id==0);
    end
    out{flt}.z=lb{flt}.z(:,eddy_id==0);
    out{flt}.ipoc_epi=lb{flt}.ipoc_epi(:,eddy_id==0);
    out{flt}.ipoc_mes=lb{flt}.ipoc_mes(:,eddy_id==0);
    out{flt}.ipoc_mes2=lb{flt}.ipoc_mes2(:,eddy_id==0);
    out{flt}.zbin.bbp700_s=lb{flt}.zbin.bbp700_s(:,eddy_id==0);
    out{flt}.zbin.bbp700_l=lb{flt}.zbin.bbp700_l(:,eddy_id==0);
    out{flt}.zbin.chla_s=lb{flt}.zbin.chla_s(:,eddy_id==0);
    out{flt}.zbin.chla_l=lb{flt}.zbin.chla_l(:,eddy_id==0);
    out{flt}.zbin.t=lb{flt}.zbin.t(:,eddy_id==0);
    out{flt}.zbin.s=lb{flt}.zbin.s(:,eddy_id==0);
    out{flt}.zbin.par=lb{flt}.zbin.par(:,eddy_id==0);
    out{flt}.zbin.oxy=lb{flt}.zbin.oxy(:,eddy_id==0);
    out{flt}.zbin.aou=lb{flt}.zbin.aou(:,eddy_id==0);
    out{flt}.zbin.nitrate=lb{flt}.zbin.nitrate(:,eddy_id==0);
    %     out{flt}.zbin.z=lb{flt}.zbin.z(:,eddy_id==0);
    out{flt}.zbin.poc_s=lb{flt}.zbin.poc_s(:,eddy_id==0);
    out{flt}.zbin.poc_l=lb{flt}.zbin.poc_l(:,eddy_id==0);
    out{flt}.zbin.pi_s=lb{flt}.zbin.pi_s{:,eddy_id==0};
    out{flt}.zbin.pi_l=lb{flt}.zbin.pi_l{:,eddy_id==0};



    % write only inside eddy data
    if sum(eddy_id)>1

        in{flt}.bbp700=lb{flt}.bbp700(:,eddy_id==1);
        in{flt}.bbp700_s=lb{flt}.bbp700_s(:,eddy_id==1);
        in{flt}.bbp700_l=lb{flt}.bbp700_l(:,eddy_id==1);
        if sum(ismember('par',fields))>0
            in{flt}.par=lb{flt}.par(:,eddy_id==1);
            in{flt}.ezd=lb{flt}.ezd(:,eddy_id==1);
            in{flt}.par0=lb{flt}.par0(:,eddy_id==1);
        end
        in{flt}.t=lb{flt}.t(:,eddy_id==1);
        in{flt}.s=lb{flt}.s(:,eddy_id==1);
        in{flt}.mld=lb{flt}.mld(:,eddy_id==1);
        in{flt}.ezd_c=lb{flt}.ezd_c(:,eddy_id==1);
        in{flt}.zprod=lb{flt}.zprod(:,eddy_id==1);
        in{flt}.vmr_100=lb{flt}.vmr_100(:,eddy_id==1);
        in{flt}.chla=lb{flt}.chla(:,eddy_id==1);
        in{flt}.chla_s=lb{flt}.chla_s(:,eddy_id==1);
        in{flt}.chla_l=lb{flt}.chla_l(:,eddy_id==1);
        in{flt}.pres=lb{flt}.pres(:,eddy_id==1);
        in{flt}.lat=lb{flt}.lat(:,eddy_id==1);
        in{flt}.lon=lb{flt}.lon(:,eddy_id==1);
        in{flt}.dnum=lb{flt}.dnum(:,eddy_id==1);
        if sum(ismember('oxy',fields))>0
            in{flt}.oxy=lb{flt}.oxy(:,eddy_id==1);
        end
        if sum(ismember('nitrate',fields))>0
            in{flt}.nitrate=lb{flt}.nitrate(:,eddy_id==1);
        end
        in{flt}.z=lb{flt}.z(:,eddy_id==1);
        in{flt}.ipoc_epi=lb{flt}.ipoc_epi(:,eddy_id==1);
        in{flt}.ipoc_mes=lb{flt}.ipoc_mes(:,eddy_id==1);
        in{flt}.ipoc_mes2=lb{flt}.ipoc_mes2(:,eddy_id==1);
        in{flt}.zbin.bbp700_s=lb{flt}.zbin.bbp700_s(:,eddy_id==1);
        in{flt}.zbin.bbp700_l=lb{flt}.zbin.bbp700_l(:,eddy_id==1);
        in{flt}.zbin.chla_s=lb{flt}.zbin.chla_s(:,eddy_id==1);
        in{flt}.zbin.chla_l=lb{flt}.zbin.chla_l(:,eddy_id==1);
        in{flt}.zbin.t=lb{flt}.zbin.t(:,eddy_id==1);
        in{flt}.zbin.s=lb{flt}.zbin.s(:,eddy_id==1);
        in{flt}.zbin.par=lb{flt}.zbin.par(:,eddy_id==1);
        in{flt}.zbin.oxy=lb{flt}.zbin.oxy(:,eddy_id==1);
        in{flt}.zbin.aou=lb{flt}.zbin.aou(:,eddy_id==1);
        in{flt}.zbin.nitrate=lb{flt}.zbin.nitrate(:,eddy_id==1);
        %         in{flt}.zbin.z=lb{flt}.zbin.z(:,eddy_id==1);
        in{flt}.zbin.poc_s=lb{flt}.zbin.poc_s(:,eddy_id==1);
        in{flt}.zbin.poc_l=lb{flt}.zbin.poc_l(:,eddy_id==1);
        in{flt}.zbin.pi_s=lb{flt}.zbin.pi_s{:,eddy_id==1};
        in{flt}.zbin.pi_l=lb{flt}.zbin.pi_l{:,eddy_id==1};

        inid=[inid; flt];
    else
        in{flt}=NaN;

    end

    clear eddy_id
end


% only store inside eddy floats with data
for i = 1:length(inid)
    in2{i}=in{inid(i)};
end
clear in
in=in2;

%%
if save_flag==1
    save(['Data/LBE_BGC_POC_2010_2022_' date '.mat'],'lb','out','in')
end