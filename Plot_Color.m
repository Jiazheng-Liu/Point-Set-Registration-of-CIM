
function Plot_Color(NewObserveData0,Latitude,Longitude,Inflected_latitudes1,Inflected_longitudes1,TruthLatitude,TruthLongitude,CentLat, CentLon)


R=2;
fanwei.Lat_min=max(floor(CentLat)-R,-89.9999);
fanwei.Lat_max=min(ceil(CentLat)+R,89.9999);
fanwei.Lon_min=max(floor(CentLon)-R,-179.9999);
fanwei.Lon_max=min(ceil(CentLon)+R,179.9999);
Scale=[-49, -46, -69, -65 ]

fanwei.Lat_min=Scale(1);
fanwei.Lat_max=Scale(2);
fanwei.Lon_min=Scale(3);
fanwei.Lon_max=Scale(4);

[X,Y] = meshgrid(1:size(NewObserveData0,2),1:size(NewObserveData0,
1));
%Adjust the size of the interpolation
Step=1/(2*4);
%Step=1/1;
[X1,Y1] = meshgrid(1:Step:size(NewObserveData0,2),1:Step:size(NewObserveData0,1));

NewObserveData0([find(Latitude<fanwei.Lat_min-2); find(Latitude>fanwei.Lat_max+2);...
find(Longitude<fanwei.Lon_min-2); find(Longitude>fanwei.Lon_max+2)]) = NaN;
Latitude([find(Latitude<fanwei.Lat_min-2); find(Latitude>fanwei.Lat_max+2)]) = NaN;%NaN
Longitude([find(Longitude<fanwei.Lon_min-2); find(Longitude>fanwei.Lon_max+2)]) = NaN;

NewLatitude = interp2(X, Y, Latitude, X1, Y1,'linear');
NewLongitude = interp2(X, Y, Longitude, X1, Y1,'linear');
NewObserveData = interp2(X, Y, NewObserveData0, X1, Y1,'cubic');
NewLatitude(:,end-(1/Step):end)=[];
NewLongitude(:,end-(1/Step):end)=[];
NewObserveData(:,end-(1/Step):end)=[];

%Start drawing

figure('color','white');
latlim = [fanwei.Lat_min fanwei.Lat_max];
lonlim = [fanwei.Lon_min fanwei.Lon_max];

axesm('mercator','MapLatLimit',latlim,'MapLonLimit',lonlim, ...
      'Frame','on', 'MeridianLabel','on','ParallelLabel','on'); 
axis off;
setm(gca,'MLabelLocation',1);%
setm(gca,'PLabelLocation',1);%

Threshold=5;% 10 previous
dlevels = floor(min(min(NewObserveData))/Threshold)*Threshold : Threshold: ceil(max(max(NewObserveData))/Threshold)*Threshold;
for k = 1:length(dlevels) - 1
   NewObserveData(find(NewObserveData>dlevels(k) & NewObserveData<=dlevels(k+1))) = k;
end
NewObserveData(find(NewObserveData==dlevels(1))) = 1;
cmap = colormap(jet(length(dlevels) - 1));
colormap(cmap);
caxis([0 length(dlevels)-1]);
cbar = colorbar;
set(cbar,'Ticks',0:1:length(dlevels)-1,'TickLabels',dlevels);

pcolorm(NewLatitude, NewLongitude, NewObserveData);
geoshow(TruthLatitude+0.05,TruthLongitude-0.10,'Color','black');%%%


hold on;

setm(gca, 'fontsize', 16);
set(gcf, 'outerposition', get(0, 'screensize'));
hold on;


end