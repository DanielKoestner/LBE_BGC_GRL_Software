%% Get floats, calculate poc

close all
clear all
curdir=cd;

addpath(genpath([curdir '/Scripts']))
addpath(genpath([curdir '/Data']))
addpath(genpath([curdir '/OneArgo']))
addpath(genpath([curdir '/aux']))

print_flag=0;
save_flag=1;
fs=16;
lw=1.5;
set(0, 'DefaultAxesFontName', 'Times');
set(0, 'DefaultTextFontName', 'Times');

alp=0.65;
alin = 0.9;



%% add floats!


load(['Data/LBE_BGCDATA_2010_2022_proc_02-May-2024.mat']);

f=length(fieldnames(LB_bbp));
names=fieldnames(LB_bbp);

for i = 1:f
    temp=eval(['LB_bbp.',char(names{i})]);
    fields=fieldnames(temp);

    %only use profiles with sampling >4 days
    tdiff=diff(temp.TIME(1,:));
    inds=find(tdiff>4);

    % inds=~isnan(temp.BBP700_ADJUSTED) | ~isnan(temp.CHLA_ADJUSTED);
    floats{i}.bbp700=temp.BBP700_despike(:,inds);
    floats{i}.bbp700_s = temp.BBP700_s(:,inds);
    floats{i}.bbp700_l = temp.BBP700_l(:,inds);
    floats{i}.bbp700_r = temp.BBP700_r;



    if sum(ismember(fields,'DOWNWELLING_PAR'))>0
        floats{i}.par=temp.DOWNWELLING_PAR(:,inds);
    else
        floats{i}.par=NaN(size(temp.BBP700));
    end

    if sum(ismember(fields,'TEMP_ADJUSTED'))>0
        floats{i}.t=temp.TEMP_ADJUSTED(:,inds);
    else
        floats{i}.t=temp.TEMP(:,inds);
    end
   
    if sum(ismember(fields,'PSAL_ADJUSTED'))>0
        floats{i}.s=temp.PSAL_ADJUSTED(:,inds);
    else
        floats{i}.s=temp.PSAL(:,inds);
    end

    floats{i}.mld=temp.MLD_DENS(:,inds);
%     floats{i}.mld(floats{i}.mld>900)=NaN; % add jan 9 2025
    clear temp

    temp=eval(['LB_chl.',char(names{i})]);
    fields=fieldnames(temp);
    floats{i}.chla=temp.chl_despike(:,inds);
    floats{i}.chla_s = temp.chl_s(:,inds);
    floats{i}.chla_l = temp.chl_l(:,inds);
    floats{i}.chla_r = temp.chl_r;
    if sum(ismember(fields,'PRES_ADJUSTED'))>0
        floats{i}.pres=temp.PRES_ADJUSTED(:,inds);
    else
        floats{i}.pres=temp.PRES(:,inds);
    end

    floats{i}.lat=temp.LATITUDE(1,inds);
    floats{i}.lon=temp.LONGITUDE(1,inds);
    floats{i}.dnum=temp.TIME(1,inds);
    floats{i}.wmo=str2num(names{i}(2:8));

    if sum(ismember(fields,'TEMP_ADJUSTED'))>0
        floats{i}.t=temp.TEMP_ADJUSTED(:,inds);
    end
    if sum(ismember(fields,'PSAL_ADJUSTED'))>0
        floats{i}.s=temp.PSAL_ADJUSTED(:,inds);
    end



    if sum(ismember(fields,'DOXY_ADJUSTED'))>0 
        floats{i}.oxy=temp.DOXY_ADJUSTED(:,inds);
    end

    if sum(ismember(fields,'NITRATE_ADJUSTED'))>0
        floats{i}.nitrate=temp.NITRATE_ADJUSTED(:,inds);
    end


    floats{i}.z=gsw_z_from_p(floats{i}.pres,floats{i}.lat);




    clear temp inds
end


% % % % % % ADD VMR CALCULATION % % % % %
for i = 1:length(floats)
    tbbp=floats{i}.bbp700;
    tz_in=floats{i}.z;

    n=length(tz_in(1,:)); % number of profiles/columns
    for ii = 1:n
        [tvmr,tz_out]=vmr_bbp(tbbp(:,ii),tz_in(:,ii));
        vmr_100(1,ii)=nanmedian(tvmr(tz_out>-100));
        clear tvmr tz_out
    end
    floats{i}.vmr_100=vmr_100*1e6; % correct
clear vmr_100;
end

%% Calculate poc + aou for floats

lb=floats;

peu=0.01; %fraction for EZD calculation, 0.5 or 1 %, decide

zbins=[0:5:1000];


% LB
for fn = 1:length(lb)
    % make single matrix with binned depths for Chla, bbp, oxy, t, and s
    
    fields=fieldnames(lb{fn});
    for i = 1:length(lb{fn}.lat)

        bs=lb{fn}.bbp700_s(:,i);
        bl=lb{fn}.bbp700_l(:,i);
        cs=lb{fn}.chla_s(:,i);
        cl=lb{fn}.chla_l(:,i);
       
        te=lb{fn}.t(:,i);
        sa=lb{fn}.s(:,i);
        z=lb{fn}.z(:,i)*-1;

        if sum(ismember(fields,'oxy'))>0
            o=lb{fn}.oxy(:,i);
            taou=f_aou(sa,te,o);
        else
            o=NaN(size(bs));
            taou=NaN(size(bs));
        end
        
        if sum(ismember(fields,'nitrate'))>0
        n=lb{fn}.nitrate(:,i);
        else
            n=NaN(size(bs));
        end

        if sum(ismember(fields,'par'))>0
        p=lb{fn}.par(:,i);
        else
            p=NaN(size(bs));
        end


        y=discretize(z,zbins);
        for ii = 1:length(zbins)
            inds=ismember(y,ii);
            if sum(inds)>1
                bbps(ii,i)=nanmean(bs(inds));
                bbpl(ii,i)=nanmean(bl(inds));
                chlas(ii,i)=nanmean(cs(inds));
                chlal(ii,i)=nanmean(cl(inds));
                par(ii,i)=nanmean(p(inds));
                oxy(ii,i)=nanmean(o(inds));
                aou(ii,i)=nanmean(taou(inds));
                nitrate(ii,i)=nanmean(n(inds));
                t(ii,i)=nanmean(te(inds));
                s(ii,i)=nanmean(sa(inds));
            elseif sum(inds)==1
                bbps(ii,i)=(bs(inds));
                bbpl(ii,i)=(bl(inds));
                chlas(ii,i)=(cs(inds));
                chlal(ii,i)=(cl(inds));
                par(ii,i)=(p(inds));
                oxy(ii,i)=(o(inds));
                aou(ii,i)=(taou(inds));
                nitrate(ii,i)=(n(inds));
                t(ii,i)=(te(inds));
                s(ii,i)=(sa(inds));
            else
                bbps(ii,i)=NaN;
                bbpl(ii,i)=NaN;
                chlas(ii,i)=NaN;
                chlal(ii,i)=NaN;
                par(ii,i)=NaN;
                oxy(ii,i)=NaN;
                aou(ii,i)=NaN;
                nitrate(ii,i)=NaN;
                t(ii,i)=NaN;
                s(ii,i)=NaN;
            end


        end

        % EZD
        mxp=max(p);
        ind0=find(p>mxp*peu,1,'last');
        par0(1,i)=mxp;
        if ~isempty(ind0)
           ezd_p(1,i)=z(ind0); % first pass at ezd;
        else
            ezd_p(1,i)=NaN;
        end

        if ezd_p(1,i)>200
            ezd_p(1,i)=NaN;
        end
        clear ind0 mxp

        % ISOLUME version, 0.415 mol quanta / m2 / day, find daylight
        dl=day_length(lb{fn}.lat(i),lb{fn}.dnum(i));
        mxp=(0.415*1e6)/(dl*60*60); % determine equivalent µmol quanta / m2 / s, isolume instantenous PAR threshold
        ind0=find(p>mxp,1,'last');
        if ~isempty(ind0)
           ezd_iso(1,i)=z(ind0); % first pass at ezd;
        else
            ezd_iso(1,i)=NaN;
        end

        if ezd_iso(1,i)>200
            ezd_iso(1,i)=NaN;
        end

        clear ind0 mxp dl

        % Chla EZD
        c=cs+cl;
        thrs=0.1*max(c(z<100)); %ensure max c value is for upper water column
        if ~isempty(thrs)
        indz=find(c>thrs,1,'last');
        else
            indz=find(c>nan);
        end

        if ~isempty(indz)
            ezd_c(1,i)=z(indz);
        else
            ezd_c(1,i)=NaN;
        end
        if ezd_c(1,i)>200
            ezd_c(1,i)=NaN;
        end

        ezd(1,i)=max([ezd_c(1,i),ezd_p(1,i),ezd_iso(1,i)]);

        [POCs(:,i),PIs{i}]=Koest23_modelB_700p(bbps(:,i),chlas(:,i),0.9);
        [POCl(:,i),PIl{i}]=Koest23_modelB_700p(bbpl(:,i),chlal(:,i),0.9);
        clear bl bs cl cs c o z p thrs taou n te sa indz 
    end
    lb{fn}.zbin.bbp700_s=bbps;
    lb{fn}.zbin.bbp700_l=bbpl;
    lb{fn}.zbin.chla_s=chlas;
    lb{fn}.zbin.chla_l=chlal;
    lb{fn}.zbin.par=par;
    lb{fn}.zbin.t=t;
    lb{fn}.zbin.s=s;
    lb{fn}.zbin.oxy=oxy;
    lb{fn}.zbin.aou=aou;
    lb{fn}.zbin.nitrate=nitrate;
    lb{fn}.zbin.z=zbins;
    lb{fn}.zbin.poc_s=POCs;
    lb{fn}.zbin.poc_l=POCl;
    lb{fn}.zbin.pi_s=PIs;
    lb{fn}.zbin.pi_l=PIl;

    lb{fn}.ezd=ezd;
    lb{fn}.par0=par0;
    lb{fn}.ezd_p=ezd_p;
    lb{fn}.ezd_c=ezd_c;
    lb{fn}.ezd_iso=ezd_iso;
    
clear ezd par0 POCs POCl PIs PIl bbps bbpl chlas chlal par t s oxy nitrate aou ezd_c ezd_iso

end


%% iPOC fractions

poclims=[1 200];

thrs=0.702; %if over 80% values in profile are nan, exclude

% Calculate iPOC
% fixed depth epi (<200 m) vs meso (>200)
% < MLD/EZD (whichever is greater) vs. > MLD/EZD


z=lb{1}.zbin.z';
for i = 1:length(lb)
    for ii = 1:length(lb{i}.dnum)

        pocs=lb{i}.zbin.poc_s(:,ii);
        pocl=lb{i}.zbin.poc_l(:,ii);
        pocl(isnan(pocl))=0;
        poc=pocs+pocl;

        ezd=lb{i}.ezd(:,ii); % max of all ezd methods
        mld=lb{i}.mld(:,ii); % mld
        %         mld(mld>900)=NaN;

        %         hzn(ii)=max([ezd_mn,mld]);
        hzn(ii)=max([ezd,mld]);
        tmphzn=hzn(ii);
        tmphzn(tmphzn>900)=NaN;

        depth=z;
        %         pocl=lb{i}.zbin.poc_l(:,ii);


        if sum(isnan(poc))<(thrs)*length(depth) %if over % samples are nan


            % Remove any NaN data for calculations and plotting
            depth=depth(~isnan(poc));
            poc=poc(~isnan(poc));

            indme2=depth>=tmphzn+200;
            if sum(indme2)>1
            mesopelagic200_poc(ii)=trapz(depth(indme2),poc(indme2));
            else
                mesopelagic200_poc(ii)=NaN;
            end

            indme=depth>=tmphzn;
            if sum(indme)>1
                mesopelagic_poc(ii)=trapz(depth(indme),poc(indme));
            else
                mesopelagic_poc(ii)=NaN;
            end
            indep=depth<tmphzn;
            if sum(indep)>1
                epipelagic_poc(ii)=trapz(depth(indep),poc(indep));
            else
                epipelagic_poc(ii)=NaN;
            end
        else
            mesopelagic200_poc(ii)=NaN;
            mesopelagic_poc(ii)=NaN;
            epipelagic_poc(ii)=NaN;

        end
        clear poc tmphzn
    end

    lb{i}.ipoc_epi=epipelagic_poc;
    lb{i}.ipoc_mes=mesopelagic_poc;
    lb{i}.ipoc_mes2=mesopelagic200_poc;
    lb{i}.zprod=hzn;
    clear  *_poc *_poc_l *_poc_s hzn
end



%% Print


if save_flag==1
    save(['Data/LBE_BGC_POC_2010_2022_' date '.mat'],'lb')
end
