

% iPOC VS TIME (EPI/MESO and EZDorMLD/SUB) inside and outside
% Glider "leaving" in May, 
%% set up


close all
clear all

curdir=cd;

addpath(genpath([curdir '/Scripts']))
addpath(genpath([curdir '/Data']))

addpath(genpath([curdir '/aux']))

print_flag=0;
save_flag=1;
scale_flag='A'; % A for 1.5r and B for 2r scaling of LBEZ

fs=15;
lw=1.5;
set(0, 'DefaultAxesFontName', 'Times');
set(0, 'DefaultTextFontName', 'Times');
alp=0.65;
mksz=7;

TEzs=[100 500];
%quantile values for upper lower limits, 
q1=0.15;
q3=0.85;

f=12; % LBE
colors=crameri('acton',f+1);
% colors=colors(2:end-2,:);
gcol=[0.2 0.45 0.85];

f=23; %LB
colors2=crameri('bamako',f+1);
colors2=colors2(2:end,:);

mk{1}='s';
mk{2}='^';
mk{3}='d';
mk{4}='v';
mk{5}='p';
mk{6}='h';
mk{7}='>';
mk{8}='<';

% month colors
c=crameri('roma',12);
colorsm=[c(9:12,:); c(1:8,:)];

warning off
%% load data

% load('LBE_BGC_POC_2010_2022_07-Jun-2024')
load(['LBE_BGC_POC_2010_2022_' scale_flag '_06-Feb-2025'])
load(['LBE_BGC_POC_2010_2022_MonthlyProfiles_' scale_flag '_06-Feb-2025'],'zPROD_in','zPROD_out')

% 
poclims=[1 200];

thrs=0.702; %if over 70% values in profile are nan, exclude. This means at least 60 5-m binned measurements are required

%% Calculate iPOC
% fixed depth epi (<200 m) vs meso (>200)
% < MLD/EZD (whichever is greater) vs. > MLD/EZD

%% INSIDE
days=[];
Z=[];
ipoc_eup=[];
ipoc_sub=[];
ipoc_subb=[];
ipoc_epi=[];
ipoc_mes=[];
ipoc_tot=[];
Z50=[];
R=[];
TH=[];
SP=[];

ipoc_eup_l=[];
ipoc_sub_l=[];
ipoc_subb_l=[];
ipoc_epi_l=[];
ipoc_mes_l=[];
ipoc_tot_l=[];
Z50_l=[];

z=-in{1}.zbin.z';
for i = 1:length(in)
    for ii = 1:length(in{i}.dnum)

        poc=in{i}.zbin.poc_s(:,ii);
        pocl=in{i}.zbin.poc_l(:,ii);
        comp=in{i}.zbin.chla_s(:,ii)./in{i}.zbin.bbp700_s(:,ii);
        compl=in{i}.zbin.chla_l(:,ii)./in{i}.zbin.bbp700_l(:,ii);


        if sum(isnan(poc))<(thrs)*length(z) %if over % samples are nan and zero
            %         if sum(isnan(poc))<(thrs)*length(z) %if over % samples are nan and zero

            days=[days; in{i}.dnum(ii)];

            zt=in{i}.zprod(ii);
            if zt>900 %| isnan(zt)
                %                 tmp=month(in{i}.dnum(ii));
                %                 zt=zPROD_in(tmp);
                zt=NaN;
            end


            [euphotic_poc,subeuphotic_poc,epipelagic_poc,mesopelagic_poc,total_poc,zp,slope,~] = POC_zParameters(z,poc,zt,TEzs,0,[in{i}.dnum(ii); in{i}.wmo],comp);
            [euphotic_poc_l,subeuphotic_poc_l,epipelagic_poc_l,mesopelagic_poc_l,total_poc_l,zp_l,~,~] = POC_zParameters(z,pocl,zt,TEzs,0,[in{i}.dnum(ii); in{i}.wmo],compl);
            %store
            Z=[Z; zt];
            Z50=[Z50; zp(3)];
            ipoc_eup=[ipoc_eup; euphotic_poc];
            ipoc_sub=[ipoc_sub; subeuphotic_poc(1)];
            ipoc_subb=[ipoc_subb; subeuphotic_poc(3)];
            ipoc_epi=[ipoc_epi; epipelagic_poc];
            ipoc_mes=[ipoc_mes; mesopelagic_poc];
            ipoc_tot=[ipoc_tot; total_poc];

            ipoc_eup_l=[ipoc_eup_l; euphotic_poc_l];
            ipoc_sub_l=[ipoc_sub_l; subeuphotic_poc_l(1)];     
            ipoc_subb_l=[ipoc_subb_l; subeuphotic_poc_l(3)];     
            ipoc_epi_l=[ipoc_epi_l; epipelagic_poc_l];
            ipoc_mes_l=[ipoc_mes_l; mesopelagic_poc_l];
            ipoc_tot_l=[ipoc_tot_l; total_poc_l];
            Z50_l=[Z50_l; zp_l(3)];

            R=[R; in{i}.relative_location(2,ii)];
            TH=[TH; in{i}.relative_location(1,ii)];
            SP=[SP; in{i}.lbe_speed(ii)];
            

        else
        end
        clear poc comp zp slope *_poc *_poc_l
    end
end



%% LB
days2=[];
Z2=[];
ipoc_eup2=[];
ipoc_sub2=[];
ipoc_subb2=[];
ipoc_epi2=[];
ipoc_mes2=[];
ipoc_tot2=[];
Z502=[];

ipoc_eup2_l=[];
ipoc_sub2_l=[];
ipoc_subb2_l=[];
ipoc_epi2_l=[];
ipoc_mes2_l=[];
ipoc_tot2_l=[];
Z502_l=[];

z=-out{1}.zbin.z';
for i = 1:length(out)
    for ii = 1:length(out{i}.dnum);

        poc=out{i}.zbin.poc_s(:,ii);
        pocl=out{i}.zbin.poc_l(:,ii);
        comp=out{i}.zbin.chla_s(:,ii)./out{i}.zbin.bbp700_s(:,ii);

        %         if sum(isnan(poc))<(thrs)*length(z) & sum(poc>0)>10 %if over % samples are nan or only 10 POC
        if sum(isnan(poc))<(thrs)*length(z) %if over % samples are nan and zero
            days2=[days2; out{i}.dnum(ii)];


            zt=out{i}.zprod(ii);
            if zt>900 % | isnan(zt)
                %                 tmp=month(out{i}.dnum(ii));
                %                 zt=zPROD_out(tmp);
                zt=NaN;
            end


           [euphotic_poc,subeuphotic_poc,epipelagic_poc,mesopelagic_poc,total_poc,zp,slope,~] = POC_zParameters(z,poc,zt,TEzs,0,[out{i}.dnum(ii); out{i}.wmo],comp);
           [euphotic_poc_l,subeuphotic_poc_l,epipelagic_poc_l,mesopelagic_poc_l,total_poc_l,zp_l,~,~] = POC_zParameters(z,pocl,zt,TEzs,0,[out{i}.dnum(ii); out{i}.wmo],comp);
           %store
            Z2=[Z2; zt];
            Z502=[Z502; zp(3)];
            ipoc_eup2=[ipoc_eup2; euphotic_poc];
            ipoc_sub2=[ipoc_sub2; subeuphotic_poc(1)];
            ipoc_subb2=[ipoc_subb2; subeuphotic_poc(3)];
            ipoc_epi2=[ipoc_epi2; epipelagic_poc];
            ipoc_mes2=[ipoc_mes2; mesopelagic_poc];
            ipoc_tot2=[ipoc_tot2; total_poc];

            Z502_l=[Z502_l; zp_l(3)];
            ipoc_eup2_l=[ipoc_eup2_l; euphotic_poc_l];
            ipoc_sub2_l=[ipoc_sub2_l; subeuphotic_poc_l(1)];
            ipoc_subb2_l=[ipoc_subb2_l; subeuphotic_poc_l(3)];
            ipoc_epi2_l=[ipoc_epi2_l; epipelagic_poc_l];
            ipoc_mes2_l=[ipoc_mes2_l; mesopelagic_poc_l];
            ipoc_tot2_l=[ipoc_tot2_l; total_poc_l];

        else
        end
        clear poc comp zp slope *_poc *_poc_l
    end
end




%% average data into each month

months=month(datetime(days,'convertfrom','datenum'));
months2=month(datetime(days2,'convertfrom','datenum'));

mns=unique(months);
for i = 1:length(mns)
    inds=ismember(months,mns(i));
    Zs(i,1)=nanmean(Z(inds));
    Z50s(i,1)=nanmean(Z50(inds));
    ipoc_tots(i,1)=nanmean(ipoc_tot(inds));
    ipoc_eups(i,1)=nanmean(ipoc_eup(inds));
    ipoc_subs(i,1)=nanmean(ipoc_sub(inds));
    ipoc_subbs(i,1)=nanmean(ipoc_subb(inds));
    ipoc_epis(i,1)=nanmean(ipoc_epi(inds));
    ipoc_mess(i,1)=nanmean(ipoc_mes(inds));

    Zs_all{i}=Z(inds);
    Z50s_all{i}=Z50(inds);
    ipoc_tots_all{i}=ipoc_tot(inds);
    ipoc_eups_all{i}=ipoc_eup(inds);
    ipoc_subs_all{i}=ipoc_sub(inds);
    ipoc_subbs_all{i}=ipoc_subb(inds);
    ipoc_epis_all{i}=ipoc_epi(inds);
    ipoc_mess_all{i}=ipoc_mes(inds);
    ratsa{i}=ipoc_eup(inds)./ipoc_tot(inds); %ratio of euphotic zone poc to total
    ratsb{i}=ipoc_epi(inds)./ipoc_tot(inds); %ratio of epipelagic poc to total
    ratsc{i}=(ipoc_epi(inds)./(ipoc_epi(inds)+ipoc_epi_l(inds))); % fraction of small to small+large in epi
    ratsd{i}=(ipoc_mes(inds)./(ipoc_mes(inds)+ipoc_mes_l(inds))); % fraction of small to small+large in meso
    
    

    Zs(i,2)=nanstd(Z(inds));
    Z50s(i,2)=nanstd(Z50(inds));
    ipoc_tots(i,2)=nanstd(ipoc_tot(inds));
    ipoc_eups(i,2)=nanstd(ipoc_eup(inds));
    ipoc_subs(i,2)=nanstd(ipoc_sub(inds));
    ipoc_subbs(i,2)=nanstd(ipoc_subb(inds));
    ipoc_epis(i,2)=nanstd(ipoc_epi(inds));
    ipoc_mess(i,2)=nanstd(ipoc_mes(inds));

    Zs(i,3)=nanmedian(Z(inds));
    Z50s(i,3)=nanmedian(Z50(inds));
    ipoc_tots(i,3)=nanmedian(ipoc_tot(inds));
    ipoc_eups(i,3)=nanmedian(ipoc_eup(inds));
    ipoc_subs(i,3)=nanmedian(ipoc_sub(inds));
    ipoc_subbs(i,3)=nanmedian(ipoc_subb(inds));
    ipoc_epis(i,3)=nanmedian(ipoc_epi(inds));
    ipoc_mess(i,3)=nanmedian(ipoc_mes(inds));


    Zs(i,4)=quantile(Z(inds),q1);
    Z50s(i,4)=quantile(Z50(inds),q1);
    ipoc_tots(i,4)=quantile(ipoc_tot(inds),q1);
    ipoc_eups(i,4)=quantile(ipoc_eup(inds),q1);
    ipoc_subs(i,4)=quantile(ipoc_sub(inds),q1);
    ipoc_subbs(i,4)=quantile(ipoc_subb(inds),q1);
    ipoc_epis(i,4)=quantile(ipoc_epi(inds),q1);
    ipoc_mess(i,4)=quantile(ipoc_mes(inds),q1);

    Zs(i,5)=quantile(Z(inds),q3);
    Z50s(i,5)=quantile(Z50(inds),q3);
    ipoc_tots(i,5)=quantile(ipoc_tot(inds),q3);
    ipoc_eups(i,5)=quantile(ipoc_eup(inds),q3);
    ipoc_subs(i,5)=quantile(ipoc_sub(inds),q3);
    ipoc_subbs(i,5)=quantile(ipoc_subb(inds),q3);
    ipoc_epis(i,5)=quantile(ipoc_epi(inds),q3);
    ipoc_mess(i,5)=quantile(ipoc_mes(inds),q3);


    [~,~,Zs(i,6:7),~]=ttest(Z(inds));
    [~,~,Z50s(i,6:7),~]=ttest(Z50(inds));
    [~,~,ipoc_tots(i,6:7),~]=ttest(ipoc_tot(inds));
    [~,~,ipoc_eups(i,6:7),~]=ttest(ipoc_eup(inds));
    [~,~,ipoc_subs(i,6:7),~]=ttest(ipoc_sub(inds));
    [~,~,ipoc_subbs(i,6:7),~]=ttest(ipoc_subb(inds));


    % start large
    Z50s_l(i,1)=nanmean(Z50_l(inds));
    ipoc_tots_l(i,1)=nanmean(ipoc_tot_l(inds));
    ipoc_eups_l(i,1)=nanmean(ipoc_eup_l(inds));
    ipoc_subs_l(i,1)=nanmean(ipoc_sub_l(inds));
    ipoc_subbs_l(i,1)=nanmean(ipoc_subb_l(inds));
    ipoc_epis_l(i,1)=nanmean(ipoc_epi_l(inds));
    ipoc_mess_l(i,1)=nanmean(ipoc_mes_l(inds));

    Z50s_l_all{i}=Z50_l(inds);
    ipoc_tots_l_all{i}=ipoc_tot_l(inds);
    ipoc_eups_l_all{i}=ipoc_eup_l(inds);
    ipoc_subs_l_all{i}=ipoc_sub_l(inds);
    ipoc_subbs_l_all{i}=ipoc_subb_l(inds);
    ipoc_epis_l_all{i}=ipoc_epi_l(inds);
    ipoc_mess_l_all{i}=ipoc_mes_l(inds);
    ratla{i}=ipoc_eup_l(inds)./ipoc_tot_l(inds);
    ratlb{i}=ipoc_epi_l(inds)./ipoc_tot_l(inds);

    Z50s_l(i,2)=nanstd(Z50_l(inds));
    ipoc_tots_l(i,2)=nanstd(ipoc_tot_l(inds));
    ipoc_eups_l(i,2)=nanstd(ipoc_eup_l(inds));
    ipoc_subs_l(i,2)=nanstd(ipoc_sub_l(inds));
    ipoc_subbs_l(i,2)=nanstd(ipoc_subb_l(inds));
    ipoc_epis_l(i,2)=nanstd(ipoc_epi_l(inds));
    ipoc_mess_l(i,2)=nanstd(ipoc_mes_l(inds));

    Z50s_l(i,3)=nanmedian(Z50_l(inds));
    ipoc_tots_l(i,3)=nanmedian(ipoc_tot_l(inds));
    ipoc_eups_l(i,3)=nanmedian(ipoc_eup_l(inds));
    ipoc_subs_l(i,3)=nanmedian(ipoc_sub_l(inds));
    ipoc_subbs_l(i,3)=nanmedian(ipoc_subb_l(inds));
    ipoc_epis_l(i,3)=nanmedian(ipoc_epi_l(inds));
    ipoc_mess_l(i,3)=nanmedian(ipoc_mes_l(inds));

    Z50s_l(i,4)=quantile(Z50_l(inds),q1);
    ipoc_tots_l(i,4)=quantile(ipoc_tot_l(inds),q1);
    ipoc_eups_l(i,4)=quantile(ipoc_eup_l(inds),q1);
    ipoc_subs_l(i,4)=quantile(ipoc_sub_l(inds),q1);
    ipoc_subbs_l(i,4)=quantile(ipoc_subb_l(inds),q1);
    ipoc_epis_l(i,4)=quantile(ipoc_epi_l(inds),q1);
    ipoc_mess_l(i,4)=quantile(ipoc_mes_l(inds),q1);

    Z50s_l(i,5)=quantile(Z50_l(inds),q3);
    ipoc_tots_l(i,5)=quantile(ipoc_tot_l(inds),q3);
    ipoc_eups_l(i,5)=quantile(ipoc_eup_l(inds),q3);
    ipoc_subs_l(i,5)=quantile(ipoc_sub_l(inds),q3);
    ipoc_subbs_l(i,5)=quantile(ipoc_subb_l(inds),q3);
    ipoc_epis_l(i,5)=quantile(ipoc_epi_l(inds),q3);
    ipoc_mess_l(i,5)=quantile(ipoc_mes_l(inds),q3);

    [~,~,Z50s_l(i,6:7),~]=ttest(Z50_l(inds));
    [~,~,ipoc_tots_l(i,6:7),~]=ttest(ipoc_tot_l(inds));
    [~,~,ipoc_eups_l(i,6:7),~]=ttest(ipoc_eup_l(inds));
    [~,~,ipoc_subs_l(i,6:7),~]=ttest(ipoc_sub_l(inds));
    [~,~,ipoc_subbs_l(i,6:7),~]=ttest(ipoc_subb_l(inds));

%     % ADD SMALL + LARGE BEFORE COMPUTING STATS
%     ipoc_eupsT(i,1)=nanmean(ipoc_eup(inds)+ipoc_eup_l(inds));
%     ipoc_eupsT(i,2)=nanstd(ipoc_eup(inds)+ipoc_eup_l(inds));
%     ipoc_eupsT(i,3)=nanmedian(ipoc_eup(inds)+ipoc_eup_l(inds));
%     ipoc_eupsT(i,4)=quantile(ipoc_eup(inds)+ipoc_eup_l(inds),q1);
%     ipoc_eupsT(i,5)=quantile(ipoc_eup(inds)+ipoc_eup_l(inds),q3);
% 
%     ipoc_subsT(i,1)=nanmean(ipoc_sub(inds)+ipoc_sub_l(inds));
%     ipoc_subsT(i,2)=nanstd(ipoc_sub(inds)+ipoc_sub_l(inds));
%     ipoc_subsT(i,3)=nanmedian(ipoc_sub(inds)+ipoc_sub_l(inds));
%     ipoc_subsT(i,4)=quantile(ipoc_sub(inds)+ipoc_sub_l(inds),q1);
%     ipoc_subsT(i,5)=quantile(ipoc_sub(inds)+ipoc_sub_l(inds),q3);
% 
%     ipoc_subbsT(i,1)=nanmean(ipoc_subb(inds)+ipoc_subb_l(inds));
%     ipoc_subbsT(i,2)=nanstd(ipoc_subb(inds)+ipoc_subb_l(inds));
%     ipoc_subbsT(i,3)=nanmedian(ipoc_subb(inds)+ipoc_subb_l(inds));
%     ipoc_subbsT(i,4)=quantile(ipoc_subb(inds)+ipoc_subb_l(inds),q1);
%     ipoc_subbsT(i,5)=quantile(ipoc_subb(inds)+ipoc_subb_l(inds),q3);

    
    
end


for i = 1:length(mns)
    inds=ismember(months2,mns(i));
    Z2s(i,1)=nanmean(Z2(inds));
    Z502s(i,1)=nanmean(Z502(inds));
    ipoc_tot2s(i,1)=nanmean(ipoc_tot2(inds));
    ipoc_eup2s(i,1)=nanmean(ipoc_eup2(inds));
    ipoc_sub2s(i,1)=nanmean(ipoc_sub2(inds));
    ipoc_subb2s(i,1)=nanmean(ipoc_subb2(inds));
    ipoc_epi2s(i,1)=nanmean(ipoc_epi2(inds));
    ipoc_mes2s(i,1)=nanmean(ipoc_mes2(inds));

    Z2s_all{i}=Z2(inds);
    Z502s_all{i}=Z502(inds);
    ipoc_tot2s_all{i}=ipoc_tot2(inds);
    ipoc_eup2s_all{i}=ipoc_eup2(inds);
    ipoc_sub2s_all{i}=ipoc_sub2(inds);
    ipoc_subb2s_all{i}=ipoc_subb2(inds);
    ipoc_epi2s_all{i}=ipoc_epi2(inds);
    ipoc_mes2s_all{i}=ipoc_mes2(inds);
    rat2sa{i}=ipoc_eup2(inds)./ipoc_tot2(inds);
    rat2sb{i}=ipoc_epi2(inds)./ipoc_tot2(inds);

    rat2sc{i}=(ipoc_epi2(inds)./(ipoc_epi2(inds)+ipoc_epi2_l(inds))); % fraction of small to small+large in epi
    rat2sd{i}=(ipoc_mes2(inds)./(ipoc_mes2(inds)+ipoc_mes2_l(inds))); % fraction of small to small+large in meso

    Z2s(i,2)=nanstd(Z2(inds));
    Z502s(i,2)=nanstd(Z502(inds));
    ipoc_tot2s(i,2)=nanstd(ipoc_tot2(inds));
    ipoc_eup2s(i,2)=nanstd(ipoc_eup2(inds));
    ipoc_sub2s(i,2)=nanstd(ipoc_sub2(inds));
    ipoc_subb2s(i,2)=nanstd(ipoc_subb2(inds));
    ipoc_epi2s(i,2)=nanstd(ipoc_epi2(inds));
    ipoc_mes2s(i,2)=nanstd(ipoc_mes2(inds));

    Z2s(i,3)=nanmedian(Z2(inds));
    Z502s(i,3)=nanmedian(Z502(inds));
    ipoc_tot2s(i,3)=nanmedian(ipoc_tot2(inds));
    ipoc_eup2s(i,3)=nanmedian(ipoc_eup2(inds));
    ipoc_sub2s(i,3)=nanmedian(ipoc_sub2(inds));
    ipoc_subb2s(i,3)=nanmedian(ipoc_subb2(inds));
    ipoc_epi2s(i,3)=nanmedian(ipoc_epi2(inds));
    ipoc_mes2s(i,3)=nanmedian(ipoc_mes2(inds));

    Z2s(i,4)=quantile(Z2(inds),q1);
    Z502s(i,4)=quantile(Z502(inds),q1);
    ipoc_tot2s(i,4)=quantile(ipoc_tot2(inds),q1);
    ipoc_eup2s(i,4)=quantile(ipoc_eup2(inds),q1);
    ipoc_sub2s(i,4)=quantile(ipoc_sub2(inds),q1);
    ipoc_subb2s(i,4)=quantile(ipoc_subb2(inds),q1);
    ipoc_epi2s(i,4)=quantile(ipoc_epi2(inds),q1);
    ipoc_mes2s(i,4)=quantile(ipoc_mes2(inds),q1);

    Z2s(i,5)=quantile(Z2(inds),q3);
    Z502s(i,5)=quantile(Z502(inds),q3);
    ipoc_tot2s(i,5)=quantile(ipoc_tot2(inds),q3);
    ipoc_eup2s(i,5)=quantile(ipoc_eup2(inds),q3);
    ipoc_sub2s(i,5)=quantile(ipoc_sub2(inds),q3);
    ipoc_subb2s(i,5)=quantile(ipoc_subb2(inds),q3);
    ipoc_epi2s(i,5)=quantile(ipoc_epi2(inds),q3);
    ipoc_mes2s(i,5)=quantile(ipoc_mes2(inds),q3);

    [~,~,Z2s(i,6:7),~]=ttest(Z2(inds));
    [~,~,Z502s(i,6:7),~]=ttest(Z502(inds));
    [~,~,ipoc_tot2s(i,6:7),~]=ttest(ipoc_tot2(inds));
    [~,~,ipoc_eup2s(i,6:7),~]=ttest(ipoc_eup2(inds));
    [~,~,ipoc_sub2s(i,6:7),~]=ttest(ipoc_sub2(inds));
    [~,~,ipoc_subb2s(i,6:7),~]=ttest(ipoc_subb2(inds));

    
    % large 
    Z502s_l(i,1)=nanmean(Z502_l(inds));
    ipoc_tot2s_l(i,1)=nanmean(ipoc_tot2_l(inds));
    ipoc_eup2s_l(i,1)=nanmean(ipoc_eup2_l(inds));
    ipoc_sub2s_l(i,1)=nanmean(ipoc_sub2_l(inds));
    ipoc_subb2s_l(i,1)=nanmean(ipoc_subb2_l(inds));
    ipoc_epi2s_l(i,1)=nanmean(ipoc_epi2_l(inds));
    ipoc_mes2s_l(i,1)=nanmean(ipoc_mes2_l(inds));

    Z502s_l_all{i}=Z502_l(inds);
    ipoc_tot2s_l_all{i}=ipoc_tot2_l(inds);
    ipoc_eup2s_l_all{i}=ipoc_eup2_l(inds);
    ipoc_sub2s_l_all{i}=ipoc_sub2_l(inds);
    ipoc_subb2s_l_all{i}=ipoc_subb2_l(inds);
    ipoc_epi2s_l_all{i}=ipoc_epi2_l(inds);
    ipoc_mes2s_l_all{i}=ipoc_mes2_l(inds);
    rat2la{i}=ipoc_eup2_l(inds)./ipoc_tot2_l(inds);
    rat2lb{i}=ipoc_epi2_l(inds)./ipoc_tot2_l(inds);

    Z502s_l(i,2)=nanstd(Z502_l(inds));
    ipoc_tot2s_l(i,2)=nanstd(ipoc_tot2_l(inds));
    ipoc_eup2s_l(i,2)=nanstd(ipoc_eup2_l(inds));
    ipoc_sub2s_l(i,2)=nanstd(ipoc_sub2_l(inds));
    ipoc_subb2s_l(i,2)=nanstd(ipoc_subb2_l(inds));
    ipoc_epi2s_l(i,2)=nanstd(ipoc_epi2_l(inds));
    ipoc_mes2s_l(i,2)=nanstd(ipoc_mes2_l(inds));

    Z502s_l(i,3)=nanmedian(Z502_l(inds));
    ipoc_tot2s_l(i,3)=nanmedian(ipoc_tot2_l(inds));
    ipoc_eup2s_l(i,3)=nanmedian(ipoc_eup2_l(inds));
    ipoc_sub2s_l(i,3)=nanmedian(ipoc_sub2_l(inds));
    ipoc_subb2s_l(i,3)=nanmedian(ipoc_subb2_l(inds));
    ipoc_epi2s_l(i,3)=nanmedian(ipoc_epi2_l(inds));
    ipoc_mes2s_l(i,3)=nanmedian(ipoc_mes2_l(inds));

    Z502s_l(i,4)=quantile(Z502_l(inds),q1);
    ipoc_tot2s_l(i,4)=quantile(ipoc_tot2_l(inds),q1);
    ipoc_eup2s_l(i,4)=quantile(ipoc_eup2_l(inds),q1);
    ipoc_sub2s_l(i,4)=quantile(ipoc_sub2_l(inds),q1);
    ipoc_subb2s_l(i,4)=quantile(ipoc_subb2_l(inds),q1);
    ipoc_epi2s_l(i,4)=quantile(ipoc_epi2_l(inds),q1);
    ipoc_mes2s_l(i,4)=quantile(ipoc_mes2_l(inds),q1);

    Z502s_l(i,5)=quantile(Z502_l(inds),q3);
    ipoc_tot2s_l(i,5)=quantile(ipoc_tot2_l(inds),q3);
    ipoc_eup2s_l(i,5)=quantile(ipoc_eup2_l(inds),q3);
    ipoc_sub2s_l(i,5)=quantile(ipoc_sub2_l(inds),q3);
    ipoc_subb2s_l(i,5)=quantile(ipoc_subb2_l(inds),q3);
    ipoc_epi2s_l(i,5)=quantile(ipoc_epi2_l(inds),q3);
    ipoc_mes2s_l(i,5)=quantile(ipoc_mes2_l(inds),q3);

    [~,~,Z502s_l(i,6:7),~]=ttest(Z502_l(inds));
    [~,~,ipoc_tot2s_l(i,6:7),~]=ttest(ipoc_tot2_l(inds));
    [~,~,ipoc_eup2s_l(i,6:7),~]=ttest(ipoc_eup2_l(inds));
    [~,~,ipoc_sub2s_l(i,6:7),~]=ttest(ipoc_sub2_l(inds));
    [~,~,ipoc_subb2s_l(i,6:7),~]=ttest(ipoc_subb2_l(inds));

%     % ADD SMALL + LARGE BEFORE COMPUTING STATS
%     ipoc_eup2sT(i,1)=nanmean(ipoc_eup2(inds)+ipoc_eup2_l(inds));
%     ipoc_eup2sT(i,2)=nanstd(ipoc_eup2(inds)+ipoc_eup2_l(inds));
%     ipoc_eup2sT(i,3)=nanmedian(ipoc_eup2(inds)+ipoc_eup2_l(inds));
%     ipoc_eup2sT(i,4)=quantile(ipoc_eup2(inds)+ipoc_eup2_l(inds),q1);
%     ipoc_eup2sT(i,5)=quantile(ipoc_eup2(inds)+ipoc_eup2_l(inds),q3);
% 
%     ipoc_sub2sT(i,1)=nanmean(ipoc_sub2(inds)+ipoc_sub2_l(inds));
%     ipoc_sub2sT(i,2)=nanstd(ipoc_sub2(inds)+ipoc_sub2_l(inds));
%     ipoc_sub2sT(i,3)=nanmedian(ipoc_sub2(inds)+ipoc_sub2_l(inds));
%     ipoc_sub2sT(i,4)=quantile(ipoc_sub2(inds)+ipoc_sub2_l(inds),q1);
%     ipoc_sub2sT(i,5)=quantile(ipoc_sub2(inds)+ipoc_sub2_l(inds),q3);
% 
%     ipoc_subb2sT(i,1)=nanmean(ipoc_subb2(inds)+ipoc_subb2_l(inds));
%     ipoc_subb2sT(i,2)=nanstd(ipoc_subb2(inds)+ipoc_subb2_l(inds));
%     ipoc_subb2sT(i,3)=nanmedian(ipoc_subb2(inds)+ipoc_subb2_l(inds));
%     ipoc_subb2sT(i,4)=quantile(ipoc_subb2(inds)+ipoc_subb2_l(inds),q1);
%     ipoc_subb2sT(i,5)=quantile(ipoc_subb2(inds)+ipoc_subb2_l(inds),q3);

    

end


%%
alpha=0.6; % face alpha for bars
bw=0.9; %bar width
fs2=11;

hf1=figure();
set(hf1,'Units','inches','Position', [5 5 11.5 7], 'PaperPosition', [0 0 11.5 8], 'PaperSize', [11.5 8]);
ha1=iSubplot(2,1, 'Gap', [0 0.02], 'Min', [0.065 0.04], 'Max', [.9 .98], 'XTickL', 'All', 'YTickL', 'All');

colors2=crameri('oslo',3);
colors2(1,:)=[0.0037 0.0051 0.0025];
colors2(3,:)=[0.3250 0.4652 0.5210];
colors2(2,:)=[0.7646 0.7629 0.5230];

colors4=crameri('nuuk',2);
% colors2=copper(3);
% dayrange=days(1):30:days(end);


axes(ha1(1));
grid on
[AX,H1,H2]=plotyy(mns',[ipoc_tots(:,1)'/1000; ipoc_tot2s(:,1)'/1000],mns',[-Z50s(:,1) -Z502s(:,1)],@bar,@bar);
axes(AX(1));hold on
box off
H3=bar(mns,[ipoc_mess(:,1)'/1000;ipoc_mes2s(:,1)'/1000]);

axes(AX(2)); hold on
box off
plot([3.5 3.5],[-1000 0],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([6.5 6.5],[-1000 0],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([9.5 9.5],[-1000 0],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
hold off

AX(2).YLabel.String='\itz^{50}_s\rm [m]';
AX(2).YLim=[-1000 0];
AX(2).YTick=[-1000:200:0];
AX(2).YTickLabel={'1000','800','600','400','200','0'};
AX(2).YMinorTick='on';
AX(1).YLabel.String='iPOC\it_s\rm [g m^{-2}]';
AX(1).YLim=[0 15];
AX(1).YTick=[0:3:15];
AX(1).YMinorTick='on';


H2(1).FaceColor=[0 0 0];
H2(1).BarWidth=bw;
H2(1).FaceAlpha=0.8;
H2(2).FaceColor=[0.2 0.2 0.2];
H2(2).BarWidth=bw;
H2(2).FaceAlpha=alpha;

AX(2).YColor=[0 0 0];

H1(1).FaceColor=colors2(2,:)*1.2;
H1(1).BarWidth=bw;
H1(1).FaceAlpha=alpha;
H1(2).FaceColor=colors2(2,:)/1.2;
H1(2).BarWidth=bw;
H1(2).FaceAlpha=alpha;
AX(1).YColor=colors2(2,:);

H3(1).FaceColor=colors2(3,:)*1.2;
H3(1).BarWidth=bw;
H3(1).FaceAlpha=alpha;
H3(2).FaceColor=colors2(3,:)/1.2;
H3(2).BarWidth=bw;
H3(2).FaceAlpha=alpha;

set(AX,'fontsize',fs)
set(AX,'ticklength',[0.015 0.015])
set(AX,'xtick',1:12)
set(AX,'xticklabel',{''})
set(AX,'xlim',[0.5 12.5])


axes(AX(1));
% in vs out
xtips1 = H1(1).XEndPoints(1);
ytips1 = H1(1).YEndPoints(1);
labels1 = 'In';
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','fontsize',fs2)

xtips2 = H1(2).XEndPoints(1);
ytips2 = H1(2).YEndPoints(1);
labels2 = 'Out';
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','fontsize',fs2,'Color',[0.2 0.2 0.2])

% epi vs meso
xtips3 = H3(1).XEndPoints(1);
ytips3 = H3(1).YEndPoints(1);
labels3 = 'epi';
text(xtips3*.97,ytips3+((ytips1-ytips3)/2.5),labels3,'HorizontalAlignment','center',...
    'VerticalAlignment','middle','Rotation',90,'fontsize',fs2)

labels4 = 'meso';
text(xtips3*.97,ytips3/2.5,labels4,'HorizontalAlignment','center',...
    'VerticalAlignment','middle','Rotation',90,'fontsize',fs2)

axes(AX(2))


axes(ha1(2));
grid on
[AX,H1,H2]=plotyy(mns',[ipoc_tots_l(:,1)'/1000; ipoc_tot2s_l(:,1)'/1000],mns',[-Z50s_l(:,1) -Z502s_l(:,1)],@bar,@bar);
axes(AX(1));hold on
H3=bar(mns,[ipoc_mess_l(:,1)'/1000;ipoc_mes2s_l(:,1)'/1000]);

axes(AX(2)); hold on
plot([3.5 3.5],[-1000 0],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([6.5 6.5],[-1000 0],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([9.5 9.5],[-1000 0],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
hold off

AX(2).YLabel.String='\itz^{50}_l\rm [m]';
AX(2).YLim=[-1000 0];
AX(2).YTick=[-1000:200:0];
AX(2).YTickLabel={'1000','800','600','400','200','0'};
AX(2).YMinorTick='on';
AX(1).YLabel.String='iPOC\it_l\rm [g m^{-2}]';
AX(1).YLim=[0 10];
AX(1).YTick=[0:2:10];
AX(1).YMinorTick='on';


H2(1).FaceColor=[0 0 0];
H2(1).BarWidth=bw;
H2(1).FaceAlpha=0.8;
H2(2).FaceColor=[0.2 0.2 0.2];
H2(2).BarWidth=bw;
H2(2).FaceAlpha=alpha;

AX(2).YColor=[0 0 0];

H1(1).FaceColor=colors2(2,:)*1.2;
H1(1).BarWidth=bw;
H1(1).FaceAlpha=alpha;
H1(2).FaceColor=colors2(2,:)/1.2;
H1(2).BarWidth=bw;
H1(2).FaceAlpha=alpha;
AX(1).YColor=colors2(2,:);

H3(1).FaceColor=colors2(3,:)*1.2;
H3(1).BarWidth=bw;
H3(1).FaceAlpha=alpha;
H3(2).FaceColor=colors2(3,:)/1.2;
H3(2).BarWidth=bw;
H3(2).FaceAlpha=alpha;

set(AX,'fontsize',fs)
set(AX,'ticklength',[0.015 0.015])
set(AX,'xtick',1:12)
set(AX,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(AX,'xlim',[0.5 12.5])


%%
alpha=0.6; % face alpha for bars
bw=0.9; %bar width
fs2=9;

hf5=figure();
set(hf5,'Units','inches','Position', [5 5 8 5], 'PaperPosition', [0 0 8 5], 'PaperSize', [8 5]);
ha5=iSubplot(2,1, 'Gap', [0 0.03], 'Min', [0.065 0.04], 'Max', [.9 .98], 'XTickL', 'All', 'YTickL', 'All');

colors2=crameri('oslo',3);
colors2(1,:)=[0.0037 0.0051 0.0025];
colors2(3,:)=[0.3250 0.4652 0.5210];
colors2(2,:)=[0.7646 0.7629 0.5230];

colors4=crameri('nuuk',2);
% colors2=copper(3);
% dayrange=days(1):30:days(end);


axes(ha5(1));
grid on
[AX,H1,H2]=plotyy(mns',[ipoc_tots(:,1)'/1000; ipoc_tot2s(:,1)'/1000],mns',[-Z50s(:,1) -Z502s(:,1)],@bar,@bar);
axes(AX(1));hold on
box off
H3=bar(mns,[ipoc_subs(:,1)'/1000;ipoc_sub2s(:,1)'/1000]);

axes(AX(2)); hold on
box off
% plot([3.5 3.5],[-1000 0],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
% plot([6.5 6.5],[-1000 0],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
% plot([9.5 9.5],[-1000 0],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
hold off

text(-0.075,1,'(a)','units','normalized','fontsize',fs-4)

AX(2).YLabel.String='\itz^{50}_s\rm [m]';
AX(2).YLim=[-1000 0];
AX(2).YTick=[-1000:200:0];
AX(2).YTickLabel={'1000','800','600','400','200','0'};
AX(2).YMinorTick='on';
AX(1).YLabel.String='iPOC\it_s\rm [g m^{-2}]';
AX(1).YLim=[0 15];
AX(1).YTick=[0:3:15];
AX(1).YMinorTick='on';


H2(1).FaceColor=[0 0 0];
H2(1).BarWidth=bw;
H2(1).FaceAlpha=0.8;
H2(2).FaceColor=[0.2 0.2 0.2];
H2(2).BarWidth=bw;
H2(2).FaceAlpha=alpha;

AX(2).YColor=[0 0 0];

H1(1).FaceColor=colors2(2,:)*1.2;
H1(1).BarWidth=bw;
H1(1).FaceAlpha=alpha;
H1(2).FaceColor=colors2(2,:)/1.2;
H1(2).BarWidth=bw;
H1(2).FaceAlpha=alpha;
AX(1).YColor=colors2(3,:);

H3(1).FaceColor=colors2(3,:)*1.2;
H3(1).BarWidth=bw;
H3(1).FaceAlpha=alpha;
H3(2).FaceColor=colors2(3,:)/1.2;
H3(2).BarWidth=bw;
H3(2).FaceAlpha=alpha;

set(AX,'fontsize',fs-2)
set(AX,'ticklength',[0.015 0.015])
set(AX,'xtick',1:12)
set(AX,'xticklabel',{''})
set(AX,'xlim',[0.4 12.6])


axes(AX(1));
% in vs out
xtips1 = H1(1).XEndPoints(1);
ytips1 = H1(1).YEndPoints(1);
labels1 = 'in';
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','fontsize',fs2)

xtips2 = H1(2).XEndPoints(1);
ytips2 = H1(2).YEndPoints(1);
labels2 = 'out';
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','fontsize',fs2,'Color',[0.2 0.2 0.2])

% epi vs meso
xtips3 = H3(1).XEndPoints(1);
ytips3 = H3(1).YEndPoints(1);
labels3 = 'PZ';
text(xtips3,ytips3+((ytips1-ytips3)/2),labels3,'HorizontalAlignment','center',...
    'VerticalAlignment','middle','Rotation',0,'fontsize',fs2)

labels4 = 'TZ';
text(xtips3,ytips3/2,labels4,'HorizontalAlignment','center',...
    'VerticalAlignment','middle','Rotation',0,'fontsize',fs2)

axes(AX(2))
 AX(1).XAxis.TickLength=[0 0];
 AX(2).XAxis.TickLength=[0 0];
tx=AX(2).YLabel.Position(1);
AX(2).YLabel.Position(1)=tx-(tx/75);

axes(ha5(2));
grid on
[AX,H1,H2]=plotyy(mns',[ipoc_tots_l(:,1)'/1000; ipoc_tot2s_l(:,1)'/1000],mns',[-Z50s_l(:,1) -Z502s_l(:,1)],@bar,@bar);
axes(AX(1));hold on; box off
H3=bar(mns,[ipoc_subs_l(:,1)'/1000;ipoc_sub2s_l(:,1)'/1000]);

text(-0.075,1,'(b)','units','normalized','fontsize',fs-4)

axes(AX(2)); hold on
% plot([3.5 3.5],[-1000 0],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
% plot([6.5 6.5],[-1000 0],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
% plot([9.5 9.5],[-1000 0],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
hold off

AX(2).YLabel.String='\itz^{50}_l\rm [m]';
AX(2).YLim=[-1000 0];
AX(2).YTick=[-1000:200:0];
AX(2).YTickLabel={'1000','800','600','400','200','0'};
AX(2).YMinorTick='on';
AX(1).YLabel.String='iPOC\it_l\rm [g m^{-2}]';
AX(1).YLim=[0 15];
AX(1).YTick=[0:3:15];
AX(1).YMinorTick='on';


H2(1).FaceColor=[0 0 0];
H2(1).BarWidth=bw;
H2(1).FaceAlpha=0.8;
H2(2).FaceColor=[0.2 0.2 0.2];
H2(2).BarWidth=bw;
H2(2).FaceAlpha=alpha;

AX(2).YColor=[0 0 0];

H1(1).FaceColor=colors2(2,:)*1.2;
H1(1).BarWidth=bw;
H1(1).FaceAlpha=alpha;
H1(2).FaceColor=colors2(2,:)/1.2;
H1(2).BarWidth=bw;
H1(2).FaceAlpha=alpha;
AX(1).YColor=colors2(3,:);

H3(1).FaceColor=colors2(3,:)*1.2;
H3(1).BarWidth=bw;
H3(1).FaceAlpha=alpha;
H3(2).FaceColor=colors2(3,:)/1.2;
H3(2).BarWidth=bw;
H3(2).FaceAlpha=alpha;


set(AX,'fontsize',fs-2)
set(AX,'ticklength',[0.015 0.015])
set(AX,'xtick',1:12)
set(AX,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(AX,'xlim',[0.4 12.6])
 AX(1).XAxis.TickLength=[0 0];
 AX(2).XAxis.TickLength=[0 0];

tx=AX(2).YLabel.Position(1);
AX(2).YLabel.Position(1)=tx-(tx/75);
%% Ratio plots
% in/out eddy small/large only meso fraction

% redefine colors



hf2=figure();
set(hf2,'Units','inches','Position', [5 5 12 6], 'PaperPosition', [0 0 12 6], 'PaperSize', [12 6]);
ha2=iSubplot(2,1, 'Gap', [0 0.015], 'Min', [0.05 0.04], 'Max', [.98 .97], 'XTickL', 'All', 'YTickL', 'All');


axes(ha2(1));
grid on
plot([3.5 3.5],[0 1],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([6.5 6.5],[0 1],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([9.5 9.5],[0 1],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([0.5 12.5],[0.5 0.5],'k:','linewidth',1.25)
hold on

text(1,1.04,'epipelagic (\itz\rm\bf < 200 m) fraction of small iPOC','units','normalized','fontsize',fs-4,'fontweight','bold','HorizontalAlignment','right')

[h2,L2,MX2,MED2,bw2]=violin(rat2sb,'x',[1.18:1:12.18],'facecolor',colors2(3,:)*1.1,'edgecolor','none','mc','k','medc','k-.');
[h,L,MX,MED,bw]=violin(ratsb,'x',[0.82:1:11.82],'facecolor',colors(10,:),'edgecolor','none','mc','k','medc','k-.');
L2.Visible='off';
text(0.01,0.25,'inside eddy','units','normalized','fontsize',fs-4,'FontWeight','bold','Color',colors(10,:));
text(0.01,0.15,'outside eddy','units','normalized','fontsize',fs-4,'FontWeight','bold','Color',colors2(3,:)*1.1);

%alpha weighting
for i = 1:12
    w=ipoc_tots(i,1)/max(ipoc_tots(:,1));
    h(i).FaceAlpha=w*0.6+0.2;
    h2(i).FaceAlpha=w*0.6+0.2;
end

ylim([0 1]);
set(gca,'fontsize',fs)
set(gca,'ticklength',[0.01 0.01])
set(gca,'xtick',1:12)
set(gca,'xticklabel',{''})
set(gca,'xlim',[0.5 12.5])

% ylabel('iPOC^{sub}\it_s\rm or iPOC^{mes}\it_s\rm / iPOC\it_s\rm');
ylabel('epi fraction of  iPOC\it_s\rm');

% AX(1).YColor=colors2(2,:);
grid on

axes(ha2(2));
grid on
plot([3.5 3.5],[0 1],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([6.5 6.5],[0 1],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([9.5 9.5],[0 1],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([0.5 12.5],[0.5 0.5],'k:','linewidth',1.25)
hold on


[h4,L4,MX4,MED4,bw4]=violin(rat2lb,'x',[1.18:1:12.18],'facecolor',colors2(3,:)*1.1,'edgecolor','none','mc','k','medc','k-.');
[h3,L3,MX3,MED3,bw3]=violin(ratlb,'x',[0.82:1:11.82],'facecolor',colors(10,:),'edgecolor','none','mc','k','medc','k-.');
L4.Location='northeast';
text(1,1.04,'epipelagic (\itz\rm\bf < 200 m) fraction of large iPOC','units','normalized','fontsize',fs-4,'fontweight','bold','HorizontalAlignment','right')


%alpha weighting
for i = 1:12
    w=ipoc_tots_l(i,1)/max(ipoc_tots_l(:,1));
    h3(i).FaceAlpha=w*0.6+0.2;
    h4(i).FaceAlpha=w*0.6+0.2;
end

ylim([0 1]);
set(gca,'fontsize',fs)
set(gca,'ticklength',[0.01 0.01])
set(gca,'xtick',1:12)
set(gca,'xticklabel',{''})
set(gca,'xlim',[0.5 12.5])
% ylabel('iPOC^{mes}\it_s\rm / iPOC\it_s\rm');
ylabel('epi fraction of  iPOC\it_l\rm');

% AX(1).YColor=colors2(2,:);
grid on

set(gca,'fontsize',fs)
set(gca,'ticklength',[0.01 0.01])
set(gca,'xtick',1:12)
set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'xlim',[0.5 12.5])



%% Ratio plots
% outisde eddy

hf3=figure();
set(hf3,'Units','inches','Position', [5 5 12 6], 'PaperPosition', [0 0 12 6], 'PaperSize', [12 6]);
ha3=iSubplot(2,1, 'Gap', [0 0.015], 'Min', [0.05 0.04], 'Max', [.98 .97], 'XTickL', 'All', 'YTickL', 'All');


axes(ha3(1));
grid on
plot([3.5 3.5],[0 1],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([6.5 6.5],[0 1],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([9.5 9.5],[0 1],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([0.5 12.5],[0.5 0.5],'k:','linewidth',1.25)
hold on

text(1,1.04,'small particle fraction of total iPOC in epipelagic (\itz\rm\bf < 200 m)','units','normalized','fontsize',fs-4,'fontweight','bold','HorizontalAlignment','right')

[h2,L2,MX2,MED2,bw2]=violin(rat2sc,'x',[1.18:1:12.18],'facecolor',colors2(3,:)*1.1,'edgecolor','none','mc','k','medc','k-.');
[h,L,MX,MED,bw]=violin(ratsc,'x',[0.82:1:11.82],'facecolor',colors(10,:),'edgecolor','none','mc','k','medc','k-.');
L2.Visible='off';
text(0.26,0.25,'inside eddy','units','normalized','fontsize',fs-4,'FontWeight','bold','Color',colors(10,:));
text(0.26,0.15,'outside eddy','units','normalized','fontsize',fs-4,'FontWeight','bold','Color',colors2(3,:)*1.1);

%alpha weighting
for i = 1:12
    w=(ipoc_epis(i,1)+ipoc_epis_l(i,1))/max(ipoc_epis(:,1)+ipoc_epis_l(:,1));
    h(i).FaceAlpha=w*0.6+0.2;
    h2(i).FaceAlpha=w*0.6+0.2;
end

ylim([0 1]);
set(gca,'fontsize',fs)
set(gca,'ticklength',[0.01 0.01])
set(gca,'xtick',1:12)
set(gca,'xticklabel',{''})
set(gca,'xlim',[0.5 12.5])

ylabel('iPOC\it_s\rm / iPOC');
% ylabel('meso fraction of  iPOC\it_s\rm');

% AX(1).YColor=colors2(2,:);
grid on

axes(ha3(2));
grid on
plot([3.5 3.5],[0 1],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([6.5 6.5],[0 1],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([9.5 9.5],[0 1],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([0.5 12.5],[0.5 0.5],'k:','linewidth',1.25)
hold on


[h4,L4,MX4,MED4,bw4]=violin(rat2sd,'x',[1.18:1:12.18],'facecolor',colors2(3,:)*1.1,'edgecolor','none','mc','k','medc','k-.');
[h3,L3,MX3,MED3,bw3]=violin(ratsd,'x',[0.82:1:11.82],'facecolor',colors(10,:),'edgecolor','none','mc','k','medc','k-.');
L4.Location='southeast';

text(1,1.04,'small particle fraction of total iPOC in mesopelagic (\itz\rm\bf > 200 m)','units','normalized','fontsize',fs-4,'fontweight','bold','HorizontalAlignment','right')

%alpha weighting
for i = 1:12
    w=(ipoc_mess(i,1)+ipoc_mess_l(i,1))/max(ipoc_mess(:,1)+ipoc_mess_l(:,1));
    h3(i).FaceAlpha=w*0.6+0.2;
    h4(i).FaceAlpha=w*0.6+0.2;
end

ylim([0 1]);
set(gca,'fontsize',fs)
set(gca,'ticklength',[0.01 0.01])
set(gca,'xtick',1:12)
set(gca,'xticklabel',{''})
set(gca,'xlim',[0.5 12.5])
% ylabel('iPOC^{mes}\it_s\rm / iPOC\it_s\rm');
ylabel('iPOC\it_s\rm / iPOC');

% AX(1).YColor=colors2(2,:);
grid on

set(gca,'fontsize',fs)
set(gca,'ticklength',[0.01 0.01])
set(gca,'xtick',1:12)
set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'xlim',[0.5 12.5])


%% Ratio plots
% in/out eddy small/large only meso fraction

% redefine colors



hf4=figure();
set(hf4,'Units','inches','Position', [5 5 12 6], 'PaperPosition', [0 0 12 6], 'PaperSize', [12 6]);
ha4=iSubplot(2,1, 'Gap', [0 0.015], 'Min', [0.05 0.04], 'Max', [.98 .97], 'XTickL', 'All', 'YTickL', 'All');


axes(ha4(1));
grid on
plot([3.5 3.5],[0 1],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([6.5 6.5],[0 1],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([9.5 9.5],[0 1],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([0.5 12.5],[0.5 0.5],'k:','linewidth',1.25)
hold on

text(1,1.04,'productive (\itz\rm\bf < MLD/EZD) fraction of small iPOC','units','normalized','fontsize',fs-4,'fontweight','bold','HorizontalAlignment','right')

[h2,L2,MX2,MED2,bw2]=violin(rat2sa,'x',[1.18:1:12.18],'facecolor',colors2(3,:)*1.1,'edgecolor','none','mc','k','medc','k-.');
[h,L,MX,MED,bw]=violin(ratsa,'x',[0.82:1:11.82],'facecolor',colors(10,:),'edgecolor','none','mc','k','medc','k-.');
L2.Visible='off';
text(0.01,0.25,'inside eddy','units','normalized','fontsize',fs-4,'FontWeight','bold','Color',colors(10,:));
text(0.01,0.15,'outside eddy','units','normalized','fontsize',fs-4,'FontWeight','bold','Color',colors2(3,:)*1.1);

%alpha weighting
for i = 1:12
    w=ipoc_tots(i,1)/max(ipoc_tots(:,1));
    h(i).FaceAlpha=w*0.6+0.2;
    h2(i).FaceAlpha=w*0.6+0.2;
end

ylim([0 1]);
set(gca,'fontsize',fs)
set(gca,'ticklength',[0.01 0.01])
set(gca,'xtick',1:12)
set(gca,'xticklabel',{''})
set(gca,'xlim',[0.5 12.5])

% ylabel('iPOC^{sub}\it_s\rm or iPOC^{mes}\it_s\rm / iPOC\it_s\rm');
ylabel('fraction of  iPOC\it_s\rm');

% AX(1).YColor=colors2(2,:);
grid on

axes(ha4(2));
grid on
plot([3.5 3.5],[0 1],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([6.5 6.5],[0 1],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([9.5 9.5],[0 1],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
plot([0.5 12.5],[0.5 0.5],'k:','linewidth',1.25)
hold on


[h4,L4,MX4,MED4,bw4]=violin(rat2la,'x',[1.18:1:12.18],'facecolor',colors2(3,:)*1.1,'edgecolor','none','mc','k','medc','k-.');
[h3,L3,MX3,MED3,bw3]=violin(ratla,'x',[0.82:1:11.82],'facecolor',colors(10,:),'edgecolor','none','mc','k','medc','k-.');
L4.Location='northeast';
text(1,1.04,'productive (\itz\rm\bf < MLD/EZD) fraction of large iPOC','units','normalized','fontsize',fs-4,'fontweight','bold','HorizontalAlignment','right')


%alpha weighting
for i = 1:12
    w=ipoc_tots_l(i,1)/max(ipoc_tots_l(:,1));
    h3(i).FaceAlpha=w*0.6+0.2;
    h4(i).FaceAlpha=w*0.6+0.2;
end

ylim([0 1]);
set(gca,'fontsize',fs)
set(gca,'ticklength',[0.01 0.01])
set(gca,'xtick',1:12)
set(gca,'xticklabel',{''})
set(gca,'xlim',[0.5 12.5])
% ylabel('iPOC^{mes}\it_s\rm / iPOC\it_s\rm');
ylabel(' fraction of  iPOC\it_l\rm');

% AX(1).YColor=colors2(2,:);
grid on

set(gca,'fontsize',fs)
set(gca,'ticklength',[0.01 0.01])
set(gca,'xtick',1:12)
set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'xlim',[0.5 12.5])




%% print

if print_flag==1
%     exportgraphics(hf1,['Figures/OSM24_LB_iPOC_TimeSeries.pdf'],'ContentType','vector');
%     exportgraphics(hf2,['Figures/OSM24_IN_iPOC_Ratios.pdf'],'ContentType','vector');
%     exportgraphics(hf3,['Figures/OSM24_OUT_iPOC_Ratios.pdf'],'ContentType','vector');

% figure(hf1)
% print('Figures/LB22_iPOC_Monthly_Summary_Epi','-dpdf','-r800')
% 
% figure(hf2)
% print('Figures/LB22_iPOC_Vertical_Epi','-dpdf','-r800')
% 
% figure(hf3)
% print('Figures/LB22_iPOC_Size','-dpdf','-r800')
% 
% figure(hf4)
% print('Figures/LB22_iPOC_Vertical_Prod','-dpdf','-r800')

figure(hf5)
print(['Figures/V8/4_LB22_iPOC_Monthly_Summary_Prod_' scale_flag],'-dpdf','-r800')
end

%% save
if save_flag==1
    if sf==1.5
        save(['Data/LBE_BGC_POC_2010_2022_MonthlyParameters_A_' date '_' num2str(TEzs(2))],'sf','Z*','ipoc*','days*','months*','rat*','R','TH','SP')
    elseif sf==2
        save(['Data/LBE_BGC_POC_2010_2022_MonthlyParameters_B_' date '_' num2str(TEzs(2))],'sf','Z*','ipoc*','days*','months*','rat*','R','TH','SP')
    end
end