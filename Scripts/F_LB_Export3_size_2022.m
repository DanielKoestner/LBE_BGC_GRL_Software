% calculate monthly average production, export, and transport efficiency.

% From Dall'olmo + Buessler papers

% 1. production rate is change in POC integrated EZ to surface with time (here
% will be monthly)

% 2. export rate is change in POC integrated EZ to ~900 m with time (TZ)

% 3. TE: Export rate for increasing depths below EZ is then claculated for
% Transport efficiency. (standard is 100 m and 500 m, start with only 200)

% We will do this for small particles only, large particles only, and total
% (small + large).



% Here we can define EZ as average max of (EZD or MLD) for floats with PAR.
% We can also look into 10% max Chla.

% What do we do about changing depth of productive layer? For example,
% iPOC(EZ) is calculated for two different depths so it will likely be
% larger for a larger EZD. The same goes for mesopelagic. Perhaps this is
% not an issue.


% command to add to end of D_LB_MeanProfs
% Z=mids;
% save('Data/LBE_BGC_POC_Monthly','POC_in','POC_in_sd','POCl_in','POCl_in_sd','AOU_in','AOU_in_sd','zPROD_in','zPROD_in_sd','POC_out','POC_out_sd','POCl_out','POCl_out_sd','AOU_out','AOU_out_sd','zPROD_out','zPROD_out_sd','Z')

% for revisions, explore a version of this figure which showcases all,
% small, large, and total information.
%% set up
close all
clear all

curdir=cd;

addpath(genpath([curdir '/Scripts']))
addpath(genpath([curdir '/Data']))

addpath(genpath([curdir '/aux']))

print_flag=1;
save_flag=0;
size_flag=0; %0 for total, 1 for small particles
scale_flag='A'; % A for 1.5r and B for 2r scaling of LBEZ

fs=15;
lw=1.5;
set(0, 'DefaultAxesFontName', 'Times');
set(0, 'DefaultTextFontName', 'Times');
alp=0.65;
mksz=7;

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
%% Try with profile by profile calculated iPOCsub

% load('LBE_BGC_POC_2010_2022_MonthlyParameters_B_31-Dec-2024_200.mat'); % July 30 used 200 m, and trap integration
% load('LBE_BGC_POC_2010_2022_MonthlyParameters_04-Aug-2024_200.mat')

load(['LBE_BGC_POC_2010_2022_MonthlyParameters_' scale_flag '_06-Feb-2025_200.mat']); 


%1 for mean, 2 for std, 3 for median, 4 for iqr
a=1;
b=2;
% % % % % IN
iPOCs_ez_in=ipoc_eups(:,a); % 1 for mean, 2 for std
iPOCl_ez_in=ipoc_eups_l(:,a); % 1 for mean, 2 for std
iPOC_ez_in=iPOCs_ez_in+iPOCl_ez_in;
iPOC_ez_in_err(1,:)=(ipoc_eups(:,a)-ipoc_eups(:,b))+(ipoc_eups_l(:,a)-ipoc_eups_l(:,b));
iPOC_ez_in_err(2,:)=(ipoc_eups(:,a)+ipoc_eups(:,b))+(ipoc_eups_l(:,a)+ipoc_eups_l(:,b));

iPOCs_tz_in=ipoc_subs(:,a); % 1 for mean
iPOCl_tz_in=ipoc_subs_l(:,a); % 1 for mean
iPOC_tz_in=iPOCs_tz_in+iPOCl_tz_in;
iPOC_tz_in_err(1,:)=(ipoc_subs(:,a)-ipoc_subs(:,b))+(ipoc_subs_l(:,a)-ipoc_subs_l(:,b));
iPOC_tz_in_err(2,:)=(ipoc_subs(:,a)+ipoc_subs(:,b))+(ipoc_subs_l(:,a)+ipoc_subs_l(:,b));

iPOCs_tz2_in=ipoc_subbs(:,a); % 1 for mean
iPOCl_tz2_in=ipoc_subbs_l(:,a); % 1 for mean
iPOC_tz2_in=iPOCs_tz2_in+iPOCl_tz2_in;
iPOC_tz2_in_err(1,:)=(ipoc_subbs(:,a)-ipoc_subbs(:,b))+(ipoc_subbs_l(:,a)-ipoc_subbs_l(:,b));
iPOC_tz2_in_err(2,:)=(ipoc_subbs(:,a)+ipoc_subbs(:,b))+(ipoc_subbs_l(:,a)+ipoc_subbs_l(:,b));

%a=3;
% % % % % OUT
iPOCs_ez_out=ipoc_eup2s(:,a); % 1 for mean
iPOCl_ez_out=ipoc_eup2s_l(:,a); % 1 for mean
iPOC_ez_out=iPOCs_ez_out+iPOCl_ez_out;
iPOC_ez_out_err(1,:)=(ipoc_eup2s(:,a)-ipoc_eup2s(:,b))+(ipoc_eup2s_l(:,a)-ipoc_eup2s_l(:,b));
iPOC_ez_out_err(2,:)=(ipoc_eup2s(:,a)+ipoc_eup2s(:,b))+(ipoc_eup2s_l(:,a)+ipoc_eup2s_l(:,b));

iPOCs_tz_out=ipoc_sub2s(:,a); % 1 for mean
iPOCl_tz_out=ipoc_sub2s_l(:,a); % 1 for mean
iPOC_tz_out=iPOCs_tz_out+iPOCl_tz_out;
iPOC_tz_out_err(1,:)=(ipoc_sub2s(:,a)-ipoc_sub2s(:,b))+(ipoc_sub2s_l(:,a)-ipoc_sub2s_l(:,b));
iPOC_tz_out_err(2,:)=(ipoc_sub2s(:,a)+ipoc_sub2s(:,b))+(ipoc_sub2s_l(:,a)+ipoc_sub2s_l(:,b));

iPOCs_tz2_out=ipoc_subb2s(:,a); % 1 for mean
iPOCl_tz2_out=ipoc_subb2s_l(:,a); % 1 for mean
iPOC_tz2_out=iPOCs_tz2_out+iPOCl_tz2_out;
iPOC_tz2_out_err(1,:)=(ipoc_subb2s(:,a)-ipoc_subb2s(:,b))+(ipoc_subb2s_l(:,a)-ipoc_subb2s_l(:,b));
iPOC_tz2_out_err(2,:)=(ipoc_subb2s(:,a)+ipoc_subb2s(:,b))+(ipoc_subb2s_l(:,a)+ipoc_subb2s_l(:,b));

%% With Monthly profiles and average EZDs
% from monthly profiles and average EZDs, doesn't work well I think because
% of EZDs
% load('LBE_BGC_POC_2010_2022_MonthlyProfiles_08-May-2024.mat')
% % % % % IN Calculate iPOCs for each month
% for i = 1:12
%     inds=Z<=zPROD_in(i);
%     iPOCs_ez_in(i)=trapz(Z(inds),POC_in(inds,i));
%     iPOCl_ez_in(i)=trapz(Z(inds),POCl_in(inds,i));
%     iPOC_ez_in(i)=trapz(Z(inds),POC_in(inds,i)+POCl_in(inds,i));
%     iPOC_ez_in_err(1,i)=trapz(Z(inds),(POC_in(inds,i)-POC_in_sd(inds,i))+(POCl_in(inds,i)-nanmedian(POCl_in_sd(inds,i))));
%     iPOC_ez_in_err(2,i)=trapz(Z(inds),(POC_in(inds,i)+POC_in_sd(inds,i))+(POCl_in(inds,i)+nanmedian(POCl_in_sd(inds,i))));
% 
%     iAOU_ez_in(i)=trapz(Z(inds),AOU_in(inds,i));
%     clear inds
% 
%     inds=Z>zPROD_in(i) & Z<900;
%     iPOCs_tz_in(i)=trapz(Z(inds),POC_in(inds,i));
%     iPOCl_tz_in(i)=trapz(Z(inds),POCl_in(inds,i));
%     iPOC_tz_in(i)=trapz(Z(inds),POC_in(inds,i)+POCl_in(inds,i));
%     iPOC_tz_in_err(1,i)=trapz(Z(inds),(POC_in(inds,i)-POC_in_sd(inds,i))+(POCl_in(inds,i)-nanmedian(POCl_in_sd(inds,i))));
%     iPOC_tz_in_err(2,i)=trapz(Z(inds),(POC_in(inds,i)+POC_in_sd(inds,i))+(POCl_in(inds,i)+nanmedian(POCl_in_sd(inds,i))));
% 
%     iAOU_tz_in(i)=trapz(Z(inds),AOU_in(inds,i));
%     clear inds
% 
%     inds=Z>(zPROD_in(i)+200) & Z<900;
%     iPOCs_tz2_in(i)=trapz(Z(inds),POC_in(inds,i));
%     iPOCl_tz2_in(i)=trapz(Z(inds),POCl_in(inds,i));
%     iPOC_tz2_in(i)=trapz(Z(inds),POC_in(inds,i)+POCl_in(inds,i));
%     iPOC_tz2_in_err(1,i)=trapz(Z(inds),(POC_in(inds,i)-POC_in_sd(inds,i))+(POCl_in(inds,i)-nanmedian(POCl_in_sd(inds,i))));
%     iPOC_tz2_in_err(2,i)=trapz(Z(inds),(POC_in(inds,i)+POC_in_sd(inds,i))+(POCl_in(inds,i)+nanmedian(POCl_in_sd(inds,i))));
%    
%     iAOU_tz2_in(i)=trapz(Z(inds),AOU_in(inds,i));
%     clear inds
% end
%     
% % % % % OUT Calculate iPOCs for each month 
% for i = 1:12
%     inds=Z<=zPROD_out(i);
%     iPOCs_ez_out(i)=trapz(Z(inds),POC_out(inds,i));
%     iPOCl_ez_out(i)=trapz(Z(inds),POCl_out(inds,i));
%     iPOC_ez_out(i)=trapz(Z(inds),POC_out(inds,i)+POCl_out(inds,i));
%     iPOC_ez_out_err(1,i)=trapz(Z(inds),(POC_out(inds,i)-POC_out_sd(inds,i))+(POCl_out(inds,i)-nanmedian(POCl_out_sd(inds,i))));
%     iPOC_ez_out_err(2,i)=trapz(Z(inds),(POC_out(inds,i)+POC_out_sd(inds,i))+(POCl_out(inds,i)+nanmedian(POCl_out_sd(inds,i))));
% 
%     iAOU_ez_out(i)=trapz(Z(inds),AOU_out(inds,i));
%     clear inds
% 
%     inds=Z>zPROD_out(i) & Z<900;
%     iPOCs_tz_out(i)=trapz(Z(inds),POC_out(inds,i));
%     iPOCl_tz_out(i)=trapz(Z(inds),POCl_out(inds,i));
%     iPOC_tz_out(i)=trapz(Z(inds),POC_out(inds,i)+POCl_out(inds,i));
%     iPOC_tz_out_err(1,i)=trapz(Z(inds),(POC_out(inds,i)-POC_out_sd(inds,i))+(POCl_out(inds,i)-nanmedian(POCl_out_sd(inds,i))));
%     iPOC_tz_out_err(2,i)=trapz(Z(inds),(POC_out(inds,i)+POC_out_sd(inds,i))+(POCl_out(inds,i)+nanmedian(POCl_out_sd(inds,i))));
% 
%     iAOU_tz_out(i)=trapz(Z(inds),AOU_out(inds,i));
%     clear inds
% 
%     inds=Z>(zPROD_out(i)+200) & Z<900;
%     iPOCs_tz2_out(i)=trapz(Z(inds),POC_out(inds,i));
%     iPOCl_tz2_out(i)=trapz(Z(inds),POCl_out(inds,i));
%     iPOC_tz2_out(i)=trapz(Z(inds),POC_out(inds,i)+POCl_out(inds,i));
%     iPOC_tz2_out_err(1,i)=trapz(Z(inds),(POC_out(inds,i)-POC_out_sd(inds,i))+(POCl_out(inds,i)-nanmedian(POCl_out_sd(inds,i))));
%     iPOC_tz2_out_err(2,i)=trapz(Z(inds),(POC_out(inds,i)+POC_out_sd(inds,i))+(POCl_out(inds,i)+nanmedian(POCl_out_sd(inds,i))));
% 
%     iAOU_tz2_out(i)=trapz(Z(inds),AOU_out(inds,i));
%     clear inds
% end

%% Rates
% remove any negative iPOC error
iPOC_ez_in_err(iPOC_ez_in_err<0)=0;
iPOC_tz_in_err(iPOC_tz_in_err<0)=0;
iPOC_tz2_in_err(iPOC_tz2_in_err<0)=0;

iPOC_ez_out_err(iPOC_ez_out_err<0)=0;
iPOC_tz_out_err(iPOC_tz_out_err<0)=0;
iPOC_tz2_out_err(iPOC_tz2_out_err<0)=0;

% % % %  rates, P and E
% dt=[31,28,31,30,31,30,31,31,30,31,30,31]; %days per month
dt=ones(1,12)*(356/12);
% dt=[29.5,29.5,30.5,30.5,30.5,30.5,31,30.5,30.5,30.5,30.5,31]; %avg days for each month (e.g., jan=(jan+feb)/2)
for i = 1:12
    if i <12
        % in
        P_in(i)= (iPOC_ez_in(i+1)-iPOC_ez_in(i))/dt(i);
        Ez_in(i)= (iPOC_tz_in(i+1)-iPOC_tz_in(i))/dt(i);
        Ez2_in(i)= (iPOC_tz2_in(i+1)-iPOC_tz2_in(i))/dt(i);

        P_in_err(:,i)= (iPOC_ez_in_err(:,i+1)-iPOC_ez_in_err(:,i))./dt(i);
        Ez_in_err(:,i)= (iPOC_tz_in_err(:,i+1)-iPOC_tz_in_err(:,i))./dt(i);
        Ez2_in_err(:,i)= (iPOC_tz2_in_err(:,i+1)-iPOC_tz2_in_err(:,i))./dt(i);

        P_s_in(i)= (iPOCs_ez_in(i+1)-iPOCs_ez_in(i))/dt(i);
        Ez_s_in(i)= (iPOCs_tz_in(i+1)-iPOCs_tz_in(i))/dt(i);
        Ez2_s_in(i)= (iPOCs_tz2_in(i+1)-iPOCs_tz2_in(i))/dt(i);
        P_l_in(i)= (iPOCl_ez_in(i+1)-iPOCl_ez_in(i))/dt(i);
        Ez_l_in(i)= (iPOCl_tz_in(i+1)-iPOCl_tz_in(i))/dt(i);
        Ez2_l_in(i)= (iPOCl_tz2_in(i+1)-iPOCl_tz2_in(i))/dt(i);
        % out
        P_out(i)= (iPOC_ez_out(i+1)-iPOC_ez_out(i))/dt(i);
        Ez_out(i)= (iPOC_tz_out(i+1)-iPOC_tz_out(i))/dt(i);
        Ez2_out(i)= (iPOC_tz2_out(i+1)-iPOC_tz2_out(i))/dt(i);

        P_out_err(:,i)= (iPOC_ez_out_err(:,i+1)-iPOC_ez_out_err(:,i))./dt(i);
        Ez_out_err(:,i)= (iPOC_tz_out_err(:,i+1)-iPOC_tz_out_err(:,i))./dt(i);
        Ez2_out_err(:,i)= (iPOC_tz2_out_err(:,i+1)-iPOC_tz2_out_err(:,i))./dt(i);

        P_s_out(i)= (iPOCs_ez_out(i+1)-iPOCs_ez_out(i))/dt(i);
        Ez_s_out(i)= (iPOCs_tz_out(i+1)-iPOCs_tz_out(i))/dt(i);
        Ez2_s_out(i)= (iPOCs_tz2_out(i+1)-iPOCs_tz2_out(i))/dt(i);
        P_l_out(i)= (iPOCl_ez_out(i+1)-iPOCl_ez_out(i))/dt(i);
        Ez_l_out(i)= (iPOCl_tz_out(i+1)-iPOCl_tz_out(i))/dt(i);
        Ez2_l_out(i)= (iPOCl_tz2_out(i+1)-iPOCl_tz2_out(i))/dt(i);
    else
% in
        P_in(i)= (iPOC_ez_in(1)-iPOC_ez_in(i))/dt(i);
        Ez_in(i)= (iPOC_tz_in(1)-iPOC_tz_in(i))/dt(i);
        Ez2_in(i)= (iPOC_tz2_in(1)-iPOC_tz2_in(i))/dt(i);

        P_in_err(:,i)= (iPOC_ez_in_err(:,1)-iPOC_ez_in_err(:,i))./dt(i);
        Ez_in_err(:,i)= (iPOC_tz_in_err(:,1)-iPOC_tz_in_err(:,i))./dt(i);
        Ez2_in_err(:,i)= (iPOC_tz2_in_err(:,1)-iPOC_tz2_in_err(:,i))./dt(i);

        P_s_in(i)= (iPOCs_ez_in(1)-iPOCs_ez_in(i))/dt(i);
        Ez_s_in(i)= (iPOCs_tz_in(1)-iPOCs_tz_in(i))/dt(i);
        Ez2_s_in(i)= (iPOCs_tz2_in(1)-iPOCs_tz2_in(i))/dt(i);
        P_l_in(i)= (iPOCl_ez_in(1)-iPOCl_ez_in(i))/dt(i);
        Ez_l_in(i)= (iPOCl_tz_in(1)-iPOCl_tz_in(i))/dt(i);
        Ez2_l_in(i)= (iPOCl_tz2_in(1)-iPOCl_tz2_in(i))/dt(i);
% out
        P_out(i)= (iPOC_ez_out(1)-iPOC_ez_out(i))/dt(i);
        Ez_out(i)= (iPOC_tz_out(1)-iPOC_tz_out(i))/dt(i);
        Ez2_out(i)= (iPOC_tz2_out(1)-iPOC_tz2_out(i))/dt(i);

        P_out_err(:,i)= (iPOC_ez_out_err(:,1)-iPOC_ez_out_err(:,i))./dt(i);
        Ez_out_err(:,i)= (iPOC_tz_out_err(:,1)-iPOC_tz_out_err(:,i))./dt(i);
        Ez2_out_err(:,i)= (iPOC_tz2_out_err(:,1)-iPOC_tz2_out_err(:,i))./dt(i);

        P_s_out(i)= (iPOCs_ez_out(1)-iPOCs_ez_out(i))/dt(i);
        Ez_s_out(i)= (iPOCs_tz_out(1)-iPOCs_tz_out(i))/dt(i);
        Ez2_s_out(i)= (iPOCs_tz2_out(1)-iPOCs_tz2_out(i))/dt(i);
        P_l_out(i)= (iPOCl_ez_out(1)-iPOCl_ez_out(i))/dt(i);
        Ez_l_out(i)= (iPOCl_tz_out(1)-iPOCl_tz_out(i))/dt(i);
        Ez2_l_out(i)= (iPOCl_tz2_out(1)-iPOCl_tz2_out(i))/dt(i);
    end
end

% Ez2_out(7)=0;
% Ez2_out(12)=0;

TE_in=Ez2_in./Ez_in;
TE_in_err=Ez2_in_err./Ez_in_err;
TE_s_in=Ez2_s_in./Ez_s_in;
TE_l_in=Ez2_l_in./Ez_l_in;

TE_out=Ez2_out./Ez_out;
TE_out_err=Ez2_out_err./Ez_out_err;
TE_s_out=Ez2_s_out./Ez_s_out;
TE_l_out=Ez2_l_out./Ez_l_out;


%% Deepest depth horizon
% load('LBE_BGC_POC_2010_2022_MonthlyParameters_04-Aug-2024_500.mat')
load(['LBE_BGC_POC_2010_2022_MonthlyParameters_' scale_flag '_06-Feb-2025_500.mat']); 

%1 for mean, 2 for std, 3 for median, 4 for iqr
% a=1;
% b=2;
% % % % % IN
iPOCs_tz_in_deep=ipoc_subs(:,a);
iPOCs_tz2_in_deep=ipoc_subbs(:,a);
iPOCl_tz_in_deep=ipoc_subs_l(:,a);
iPOCl_tz2_in_deep=ipoc_subbs_l(:,a);
iPOC_tz_in_deep=ipoc_subs(:,a)+ipoc_subs_l(:,a);
iPOC_tz2_in_deep=ipoc_subbs(:,a)+ipoc_subbs_l(:,a);

%a=3;
% % % % % OUT
iPOCs_tz_out_deep=ipoc_sub2s(:,a);
iPOCs_tz2_out_deep=ipoc_subb2s(:,a);
iPOCl_tz_out_deep=ipoc_sub2s_l(:,a);
iPOCl_tz2_out_deep=ipoc_subb2s_l(:,a);
iPOC_tz_out_deep=ipoc_sub2s(:,a)+ipoc_sub2s_l(:,a);
iPOC_tz2_out_deep=ipoc_subb2s(:,a)+ipoc_subb2s_l(:,a);

for i = 1:12
    if i <12
        % in
        Ez_s_in_deep(i)= (iPOCs_tz_in_deep(i+1)-iPOCs_tz_in_deep(i))/dt(i);
        Ez2_s_in_deep(i)= (iPOCs_tz2_in_deep(i+1)-iPOCs_tz2_in_deep(i))/dt(i);

        Ez_l_in_deep(i)= (iPOCl_tz_in_deep(i+1)-iPOCl_tz_in_deep(i))/dt(i);
        Ez2_l_in_deep(i)= (iPOCl_tz2_in_deep(i+1)-iPOCl_tz2_in_deep(i))/dt(i);

        Ez_in_deep(i)= (iPOC_tz_in_deep(i+1)-iPOC_tz_in_deep(i))/dt(i);
        Ez2_in_deep(i)= (iPOC_tz2_in_deep(i+1)-iPOC_tz2_in_deep(i))/dt(i);

        % out
        Ez_s_out_deep(i)= (iPOCs_tz_out_deep(i+1)-iPOCs_tz_out_deep(i))/dt(i);
        Ez2_s_out_deep(i)= (iPOCs_tz2_out_deep(i+1)-iPOCs_tz2_out_deep(i))/dt(i);
        
        Ez_l_out_deep(i)= (iPOCl_tz_out_deep(i+1)-iPOCl_tz_out_deep(i))/dt(i);
        Ez2_l_out_deep(i)= (iPOCl_tz2_out_deep(i+1)-iPOCl_tz2_out_deep(i))/dt(i);

        Ez_out_deep(i)= (iPOC_tz_out_deep(i+1)-iPOC_tz_out_deep(i))/dt(i);
        Ez2_out_deep(i)= (iPOC_tz2_out_deep(i+1)-iPOC_tz2_out_deep(i))/dt(i);

    else
% in
        Ez_s_in_deep(i)= (iPOCs_tz_in_deep(1)-iPOCs_tz_in_deep(i))/dt(i);
        Ez2_s_in_deep(i)= (iPOCs_tz2_in_deep(1)-iPOCs_tz2_in_deep(i))/dt(i);
        
        Ez_l_in_deep(i)= (iPOCl_tz_in_deep(1)-iPOCl_tz_in_deep(i))/dt(i);
        Ez2_l_in_deep(i)= (iPOCl_tz2_in_deep(1)-iPOCl_tz2_in_deep(i))/dt(i);

        Ez_in_deep(i)= (iPOC_tz_in_deep(1)-iPOC_tz_in_deep(i))/dt(i);
        Ez2_in_deep(i)= (iPOC_tz2_in_deep(1)-iPOC_tz2_in_deep(i))/dt(i);
% out
        Ez_s_out_deep(i)= (iPOCs_tz_out_deep(1)-iPOCs_tz_out_deep(i))/dt(i);
        Ez2_s_out_deep(i)= (iPOCs_tz2_out_deep(1)-iPOCs_tz2_out_deep(i))/dt(i);

        Ez_l_out_deep(i)= (iPOCl_tz_out_deep(1)-iPOCl_tz_out_deep(i))/dt(i);
        Ez2_l_out_deep(i)= (iPOCl_tz2_out_deep(1)-iPOCl_tz2_out_deep(i))/dt(i);

        Ez_out_deep(i)= (iPOC_tz_out_deep(1)-iPOC_tz_out_deep(i))/dt(i);
        Ez2_out_deep(i)= (iPOC_tz2_out_deep(1)-iPOC_tz2_out_deep(i))/dt(i);

    end
end

TE_s_in_deep=Ez2_s_in_deep./Ez_s_in_deep;
TE_s_out_deep=Ez2_s_out_deep./Ez_s_out_deep;

TE_l_in_deep=Ez2_l_in_deep./Ez_l_in_deep;
TE_l_out_deep=Ez2_l_out_deep./Ez_l_out_deep;

TE_in_deep=Ez2_in_deep./Ez_in_deep;
TE_out_deep=Ez2_out_deep./Ez_out_deep;



%% Plot

% reset variables for plotting for smal particles only if desired.
if size_flag==1
    Ez2_in=Ez2_s_in;
    Ez_in=Ez_s_in;
    P_in=P_s_in;
    TE_in=TE_s_in;
    TE_in_deep=TE_s_in_deep;
    iPOC_tz_in=iPOCs_tz_in;


    Ez2_out=Ez2_s_out;
    Ez_out=Ez_s_out;
    P_out=P_s_out;
    TE_out=TE_s_out;
    TE_out_deep=TE_s_out_deep;
    iPOC_tz_out=iPOCs_tz_out;

end

colors=crameri('davos',7);
% colp=[0.4 0.75 0.4];
% coltz2=colors(2,:);
coltz=colors(3,:);

colp=[80/255 175/255 200/255];
coltz2=[18/255 22/255 20/255];
coltz=[50/255 75/255 70/255]*0.8;

col2tz2=[52/255   57/255    55/255];
col2tz=[124/255   139/255    136/255];

hf=figure();
set(hf,'Units','inches','Position', [5 5 8 5.25], 'PaperPosition', [0 0 8 5.25], 'PaperSize', [8 5.25]);
ha1=iSubplot(3,1, 'Gap', [0 0.04], 'Min', [0.065 0.03], 'Max', [0.90 0.96], 'XTickL', 'All', 'YTickL', 'All');

axes(ha1(1))
hold on
% plot([3.5 3.5],[-1000 1000],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
% plot([6.5 6.5],[-1000 1000],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
% plot([9.5 9.5],[-1000 1000],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
text(0.96,0.92,'(c)','units','normalized','fontsize',fs-4)
% plot(1:12,P_in_err,'-','Color',colp,'linewidth',0.5)
% plot(1:12,P_l_in,'-','Color',colp,'linewidth',0.5)

% p2=plot(1:12,P_s_in,'-','Color',colp,'linewidth',2);
% p3=plot(1:12,P_l_in,':','Color',colp,'linewidth',2);

b2=bar(Ez2_in,'FaceColor',coltz2,'FaceAlpha',0.8,'linestyle','-','edgecolor',[0.7 0.6 0.3],'linewidth',1.25);
b1=bar(Ez_in,'FaceColor',coltz,'FaceAlpha',0.55,'edgecolor',[0.5 0.5 0.7],'linewidth',1.5);

bar(1:12,[Ez_s_in;Ez_l_in],0.6,'FaceColor',col2tz,'FaceAlpha',1,'edgecolor',[0.7579    0.7579    0.8670],'linewidth',1);
bar(1:12,[Ez2_s_in;Ez2_l_in],0.6,'FaceColor',col2tz2,'FaceAlpha',1,'edgecolor',[0 0 0],'linestyle','-','edgecolor',[0.7 0.6 0.3],'linewidth',1);



p1=plot(1:12,P_in,'--','Color',colp,'linewidth',2.5);

set(gca,'ticklength',[0.01 0.005])
box on

ha1(1).YGrid='on';
set(gca,'fontsize',fs-2)
% xlim([0 13])
xlim([0.4 12.6])
ylim([-80 160])
set(gca,'ytick',[-80:40:160])
ylabel('\itE \rm[mg m^{-2} d^{-1}]')
set(gca,'xtick',[0.5:11.5],'xticklabel',{''});
ha1(1).YMinorTick='on';
title('inside LBEZ','fontsize',fs-3)

[hl,lines]=legendflex([p1 b1 b2],{'PZ','TZ','lower-TZ'}, 'ref', gca, ...
'anchor', {'nw', 'nw'}, 'buffer', [5 -2], 'xscale', 0.5, 'ncol',1,'nrow',3, ...
'box', 'off', 'FontSize', fs-3);
PatchInLegend = findobj(hl, 'type', 'patch');
PatchInLegend(2).FaceAlpha=0.55;
PatchInLegend(1).FaceAlpha=0.8;

% [hl,lines]=legendflex([p1 p2 p3],{'PZ','PZ small','PZ large'}, 'ref', gca, ...
% 'anchor', {'nw', 'nw'}, 'buffer', [5 -2], 'xscale', 0.5, 'ncol',1,'nrow',3, ...
% 'box', 'off', 'FontSize', fs-3);




axes(ha1(2))
% plot([3.5 3.5],[-1000 1000],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
% plot([6.5 6.5],[-1000 1000],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
% plot([9.5 9.5],[-1000 1000],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
text(0.96,0.92,'(d)','units','normalized','fontsize',fs-4)

hold on
ha1(2).YGrid='on';

% plot(1:12,P_out_err,'-','Color',colp,'linewidth',0.5)
% plot(1:12,P_l_out,'-','Color',colp,'linewidth',0.5)
% p2=plot(1:12,P_s_out,'-','Color',colp,'linewidth',2);
% p3=plot(1:12,P_l_out,':','Color',colp,'linewidth',2);

b2=bar(Ez2_out,'FaceColor',coltz2,'FaceAlpha',0.8,'linestyle','-','edgecolor',[0.7 0.6 0.3],'linewidth',1.25);
b1=bar(Ez_out,'FaceColor',coltz,'FaceAlpha',0.55,'edgecolor',[0.5 0.5 0.7],'linewidth',1.5);

bar(1:12,[Ez_s_out;Ez_l_out],0.6,'FaceColor',col2tz,'FaceAlpha',1,'edgecolor',[0.7579    0.7579    0.8670],'linewidth',1);
bar(1:12,[Ez2_s_out;Ez2_l_out],0.6,'FaceColor',col2tz2,'FaceAlpha',1,'edgecolor',[0 0 0],'linestyle','-','edgecolor',[0.7 0.6 0.3],'linewidth',1);

p1=plot(1:12,P_out,'--','Color',colp,'linewidth',2.5);

xlim([0.4 12.6])
ylim([-80 160])
set(gca,'ytick',[-80:40:160])
box on
ha1(2).YMinorTick='on';
set(gca,'fontsize',fs-2)
ylabel('\itE\rm [mg m^{-2} d^{-1}]')
set(gca,'xtick',[0.5:11.5],'xticklabel',{''});

title('outside LBEZ','fontsize',fs-3)

% text(0.7,0.87,sprintf('\\Sigma \\itE^{PZ}\\rm = %1.0f (%1.0f – %1.0f) mg m^{-2} y^{-1}',sum(P_out)*365,sum(P_out_err(1,:))*365,sum(P_out_err(2,:))*365),'units','normalized','color',[0.4 0.75 0.4],'fontsize',fs-4)
% text(0.7,0.72,sprintf('\\Sigma \\itE^{TZ}\\rm = %1.0f (%1.0f – %1.0f) mg m^{-2} y^{-1}',sum(Ez_out)*365,sum(Ez_out_err(1,:))*365,sum(Ez_out_err(2,:))*365),'units','normalized','color',colors(4,:),'fontsize',fs-4)
% text(0.76,0.695,sprintf('%1.2f (%1.1f – %1.1f) g m^{-2} yr^{-1}',sum(Ez2_out)*365/1000,sum(Ez2_out_err(1,:))*365/1000,sum(Ez2_out_err(2,:))*365/1000),'units','normalized','color',colors(2,:),'fontsize',fs-4)


axes(ha1(3))
colors2=crameri('batlow',8);
% plot([3.5 3.5],[-1000 1000],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
% plot([6.5 6.5],[-1000 1000],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
% plot([9.5 9.5],[-1000 1000],'--','color',[0.4 0.4 0.4],'linewidth',0.9)
text(0.96,0.92,'(e)','units','normalized','fontsize',fs-4)

hold on
% plot([0 13],[0 0],'k-')

% % % % % % % % SHADING - NEED TO SAVE AND LOAD ALL DEPTH HORIZONS

% out first, then do in
x=1:11;
y=TE_out(1:11); % set to shallowest depth horizon (100 or 200)
ydeep=TE_out_deep(1:11); %set to deepest depth horizon (400 or 500)
% ydeep(12)=y(12);
xx = [x;x];
yy = [y;ydeep];
% alps=yy*0.5; alps(alps<0)=0; alps(alps>0.4)=0.4;
% alps(1,12)=0.1;
alps=[ones(1,11)*0.4;ones(1,11)*0.4];
s2=surf(xx,yy,xx*0,...
    'alphadata',alps,...
    'facealpha','interp',...
    'edgecolor','none','FaceColor',colors2(3,:));

% x=8:11;
% y=TE_out(8:11); % set to shallowest depth horizon (100 or 200)
% ydeep=TE_out_deep(8:11); %set to deepest depth horizon (400 or 500)
% % ydeep(12)=y(12);
% xx = [x;x];
% yy = [y;ydeep];
% alps=yy*0.5; alps(alps<0)=0; alps(alps>0.4)=0.4;
% % alps(1,12)=0.1;
% % alps=[ones(1,12)*0.4;ones(1,12)*0.3];
% s2=surf(xx,yy,xx*0,...
%     'alphadata',alps,...
%     'facealpha','interp',...
%     'edgecolor','none','FaceColor',colors2(3,:));

% % line([0 x x(end) 0],[0 y 0 0],'linew',2)
p2=plot(1:11,TE_out(1:11),'color',colors2(3,:),'linewidth',0.75,'marker','s','MarkerFaceColor',colors2(3,:),'MarkerSize',1,'MarkerEdgeColor','none');
% plot([6 8],[(TE_out(6)+TE_out_deep(6))/2 (TE_out(8)+TE_out_deep(8))/2],'color',colors2(3,:),'linewidth',0.75,'linestyle','--')
% p2=plot(8:11,TE_out(8:11),'color',colors2(3,:),'linewidth',0.75,'marker','s','MarkerFaceColor',colors2(3,:),'MarkerSize',1,'MarkerEdgeColor','none');

% % % % % % 

% p2=plot(1:12,TE_out,'color',colors2(3,:),'linewidth',0.75,'marker','s','MarkerFaceColor',colors2(3,:),'MarkerSize',1,'MarkerEdgeColor','none');
% plot(1:12,TE_s_out,':','color',colors2(3,:),'linewidth',0.75,'marker','^','MarkerFaceColor',colors2(3,:),'MarkerSize',1,'MarkerEdgeColor','none');
% plot(1:12,TE_l_out,'--','color',colors2(3,:),'linewidth',0.75,'marker','v','MarkerFaceColor',colors2(3,:),'MarkerSize',1,'MarkerEdgeColor','none');


% for i = 1:12
%     errorbar(i,TE_out(i),TE_out_err(1,i)-TE_out(i),TE_out(i)-TE_out_err(2,i),'color',colors2(3,:),'CapSize',1,'markersize',0.1,'linewidth',0.25)
% end

coldeep=brighten(colors2(3,:),0.65);
plot(1:11,TE_out_deep(1:11),'color',coldeep,'linewidth',0.1);
% plot(8:11,TE_out_deep(8:11),'color',coldeep,'linewidth',0.1);
scatter(1:12,TE_out(1:12),iPOC_tz_out(1:12)/115,'s','filled','MarkerFaceColor',colors2(3,:),'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceAlpha',0.9);
% scatter(8:11,TE_out(8:11),iPOC_tz_out(8:11)/115,'s','filled','MarkerFaceColor',colors2(3,:),'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceAlpha',0.9);

x=0:12;
y=[TE_in(12) TE_in(1:12)]; % set to shallowest depth horizon (100 or 200)
ydeep=[TE_in_deep(12) TE_in_deep(1:12)]; %set to deepest depth horizon (400 or 500)
xx = [x;x];
yy = [y;ydeep];
% alps=yy*0.5; alps(alps<0)=0; alps(alps>0.4)=0.4;
alps=[ones(1,13)*0.4;ones(1,13)*0.4];
s1=surf(xx,yy,xx*0,...
    'alphadata',alps,...
    'facealpha','interp',...
    'edgecolor','none','FaceColor',colors2(7,:));

x=12:13;
y=[TE_in(12) TE_in(1)]; % set to shallowest depth horizon (100 or 200)
ydeep=[TE_in_deep(12) TE_in_deep(1)]; %set to deepest depth horizon (400 or 500)
xx = [x;x];
yy = [y;ydeep];
% alps=yy*0.5; alps(alps<0)=0; alps(alps>0.4)=0.4;
alps=[ones(1,2)*0.4;ones(1,2)*0.4];
s2=surf(xx,yy,xx*0,...
    'alphadata',alps,...
    'facealpha','interp',...
    'edgecolor','none','FaceColor',colors2(7,:));

coldeep=brighten(colors2(7,:),0.65);
p1=plot(1:12,TE_in(1:12),'color',colors2(7,:),'linewidth',0.75,'marker','o','MarkerFaceColor',colors2(7,:),'MarkerSize',1,'MarkerEdgeColor','none');
plot(1:12,TE_in_deep(1:12),'color',coldeep,'linewidth',0.1);

% p1=plot(12,TE_in(12),'color',colors2(7,:),'linewidth',0.75,'marker','o','MarkerFaceColor',colors2(7,:),'MarkerSize',1,'MarkerEdgeColor','none');
% plot(12,TE_in_deep(12),'color',coldeep,'linewidth',0.1);

% plot(1:12,TE_s_in,':','color',colors2(7,:),'linewidth',0.75,'marker','^','MarkerFaceColor',colors2(3,:),'MarkerSize',1,'MarkerEdgeColor','none');
% plot(1:12,TE_l_in,'--','color',colors2(7,:),'linewidth',0.75,'marker','v','MarkerFaceColor',colors2(3,:),'MarkerSize',1,'MarkerEdgeColor','none');


% for i = 1:12
%     errorbar(i,TE_in(i),TE_in_err(1,i)-TE_in(i),TE_in(i)-TE_in_err(2,i),'color',colors2(7,:),'CapSize',0,'markersize',0.1,'linewidth',0.25)
% end

scatter(1:12,TE_in,iPOC_tz_in/115,'o','filled','MarkerFaceColor',colors2(7,:),'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceAlpha',0.9);

set(gca,'ticklength',[0.01 0.005])
box on
% grid on
ha1(3).YGrid='on';
ha1(3).YMinorTick='on';
set(gca,'ytick',[0:0.25:2])
set(gca,'fontsize',fs-2)
xlim([0.4 12.6])
ylim([0 1.25])
ylabel('TE')
set(gca,'xtick',[0.5:11.5],'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});



title('Transfer Efficiency','fontsize',fs-3)

[hl,lines]=legendflex([p1 p2],{'inside','outside'}, 'ref', gca, ...
'anchor', {'nw', 'nw'}, 'buffer', [5 -4], 'xscale', 0.5, 'ncol',1,'nrow',3, ...
'box', 'off', 'FontSize', fs-3);
lines(4).MarkerSize=4;
lines(6).MarkerSize=4;
%% print

if print_flag==1
    if size_flag==0
    print(['Figures/V8/5_LB22_MonthlyExport_200_500_' scale_flag],'-dpdf','-r800')
%     elseif size_flag==1
%     print(['Figures/V5/5_LB22_MonthlyExport_200_500_' scale_flag '_small'],'-dpdf','-r800')
    end    
end
