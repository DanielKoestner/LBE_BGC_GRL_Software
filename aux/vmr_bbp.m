% function for variance-to-mean ratio calculation, following Briggs et al
% 2013 AO

% input: one profile or time-series of bbp (despiked with some threshold
% ~0.01), correspoding depth (or time) bins, and a window (default is 11)

% other inputs include backscattering efficiency factor (Qbb), volume of
% integration in ml (V_ml), and volume spreading factor (alp)

% output: diameter (µm) and corresponding depth/time steps

% Version: 20240723 – DK
function [D,steps_out]=vmr_bbp(bbp,steps_in,window,Qbb,V_ml,alp)

% % from Rembauville et al. 2017
% tres=0.1; %seconds
% tsample=1; %seconds
% tau=tres/tsample;
% if tau>=1
%     alp=1-(1/(3*tau));
% elseif tau<1
%     alp=tau-((tau^2)/3);
% end


if nargin==2
    window=11;
    Qbb=0.024;
    V_ml=10; %standard for bgc-argo?
    alp=1; %no correction for volume spreading
elseif nargin==3
    Qbb=0.024;
    V_ml=10; %standard for bgc-argo?
    alp=1; %no correction for volume spreading
end


V_m3=V_ml*1e-6;

% remove NaN data
steps_out=steps_in(~isnan(bbp));
bbp=bbp(~isnan(bbp));

% calculate vmr
bbp_mean=movmean(bbp,window); % moving mean for denominator 
detrend=movmean(movmedian(bbp,window),window); %detrending filter
bbp_detrend=bbp-detrend; % bbp, detrended
bbp_var=movvar(bbp_detrend,window); % moving variance of detrended bbp

vmr=bbp_var./bbp_mean;

A=vmr*(V_m3/Qbb)*(1/alp); % area calculation
D_m=2*sqrt(A/pi); % equivalent spherical diameter in m
D=D_m*1e6; % in µm
