%% function to find in/out eddy

% argo date in matlab datenum format (days since 0-0-0000)
% eddy_id = 1 if inside eddy, 0 if outside
function [eddy_id]=in_eddy(argo_date,argo_lat,argo_lon,LBE_locations)

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


% polygon approach
pgon=polyshape(eddy_lons,eddy_lats);
pgon2=scale(pgon,1.5,LBE_center(:,ind)');
eddy_id=isinterior(pgon2,argo_lon,argo_lat);

elat=0.25; %assume 50% error in lat height (~1º, *0.5/2
elon=0.75; %assume 50% error in lon width (~3º *0.5/2)

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