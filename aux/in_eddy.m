%% function to find in/out eddy

% argo date in matlab datenum format (days since 0-0-0000)
% eddy_id = 1 if inside eddy, 0 if outside
% outputs also relative_location in polar coordinates, theta (radians) and
% r (km), and rotational speed in m/s
function [eddy_id,relative_location,speed]=in_eddy(argo_date,argo_lat,argo_lon,LBE_locations,factor)

load(LBE_locations);

% find closest eddy location date
[~,ind]=min(abs(argo_date-LBE_dnum));
ind=ind;
eddy_lats=LBE_lats(:,ind);

% if no eddy data, use next date
if eddy_lats(1)==0
    ind=ind+1;
end

eddy_lats=LBE_lats(:,ind);
eddy_lons=LBE_lons(:,ind);

speed=LBE_speed(ind);

% polygon approach
pgon=polyshape(eddy_lons,eddy_lats);
pgon2=scale(pgon,factor,LBE_center(:,ind)');
eddy_id=isinterior(pgon2,argo_lon,argo_lat);

elat=0.25; %assume 50% error in lat height (~1º, *0.5/2
elon=0.75; %assume 50% error in lon width (~3º *0.5/2)


% % % % % %  relative location in km

% if eddy_id==1

% Haversine formula to calculate distance
R = 6371; % Earth's radius in kilometers
dlat = deg2rad(argo_lat - LBE_center(2,ind));
dlon = deg2rad(argo_lon - LBE_center(1,ind));
a = sin(dlat/2).^2 + cos(deg2rad(LBE_center(2,ind))) .* cos(deg2rad(argo_lat)) .* sin(dlon/2).^2;
c = asin(sqrt(a));
rho = 2*R * c; % Distance in kilometers

% Convert lat/lon to polar coordinates relative to the fixed center
theta=atan2(dlat,dlon); % in radians

relative_location(1,1)=theta;
relative_location(2,1)=rho;

% else
%     relative_location=NaN(2,1);
% end

% elat=0; 
% elon=0;

% % simple rectangle approach?
% % is float within bounds?
% if argo_lat > min(eddy_lats)-elat && argo_lat < max(eddy_lats)+elat
%     if argo_lon > min(eddy_lons)-elon && argo_lon < max(eddy_lons)+elon
%         eddy_id=1; % INSIDE eddy
%     else
%         eddy_id=0; % OUTSIDE eddy
%     end
% else
%     eddy_id=0; % OUTSIDE eddy
% end

eddy=pgon;
eddy2=pgon2;

% % % plot test
% %
% figure(); box on
% title(string(LBE_dates(ind)));
% hold on
% ylabel('lat')
% xlabel('lon')
% plot(pgon)
% plot(pgon2)
% scatter(LBE_center(1,ind),LBE_center(2,ind),'kx')
% 
% [x,y]=meshgrid([min(eddy_lons) max(eddy_lons)],[min(eddy_lats) max(eddy_lats)]);
% plot(x,y,'b--');
% plot(x',y','b--');
% 
% [x,y]=meshgrid([min(eddy_lons)-elon max(eddy_lons)+elon],[min(eddy_lats)-elat max(eddy_lats)+elat]);
% plot(x,y,'r--');
% plot(x',y','r--');


end