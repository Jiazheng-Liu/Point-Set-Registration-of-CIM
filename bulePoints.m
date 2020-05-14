%Extract the coordinate information of the GSHHS corresponding to the detected coastline from the  .hdf into .mat 
%%
clc;clear;close all;warning off;
data=hdf5info('FY3C_MWRIA_GBAL_L1_20180904_1818_010KM_MS.hdf');
EARTH_OBSERVE_BT_10_to_89GHz=hdf5read(data.GroupHierarchy.Groups(1).Datasets(2));%data observation
Daycnt=hdf5read(data.GroupHierarchy.Groups(1).Datasets(5));%Day count
Mscnt=hdf5read(data.GroupHierarchy.Groups(1).Datasets(6));%Seconds count
%%
Output.Daycnt=Daycnt;
Output.Mscnt=Mscnt;
%%
LandSeaMask=hdf5read(data.GroupHierarchy.Groups(1).Datasets(4));%Surface type
LandCover=hdf5read(data.GroupHierarchy.Groups(1 ).Datasets(3));
%%
ObserveData=EARTH_OBSERVE_BT_10_to_89GHz(:,:,10);%Extract 89GHz data


ObserveData1=EARTH_OBSERVE_BT_10_to_89GHz(:,:,9);%Extract 89GHz data
ObserveData=double(ObserveData)*0.01+327.68;%Change to positive value
%%
Latitude=double(hdf5read(data.GroupHierarchy.Groups(2).Datasets(1)));%latitude
Longitude=double(hdf5read(data.GroupHierarchy.Groups(2).Datasets(2)));%longitude

LatCenter=(Latitude(127,:)+Latitude(128,:))/2;
LonCenter=(Longitude(127,:)+Longitude(128,:))/2;
TrackDir=[diff(LatCenter'),diff(LonCenter')];
HeadingAngle=acos(abs(TrackDir(:,1))./sqrt(TrackDir(:,1).^2+TrackDir(:,2).^2));
HeadingAngle=[HeadingAngle;HeadingAngle(end)];
HeadingAngle1=HeadingAngle*180/pi;

%%Calculate data between latitude -75 to +75 degrees
ZeroPosition=find(abs(Latitude)>125);
Latitude(ZeroPosition)=NaN;
Longitude(ZeroPosition)=NaN;

%%
x=double(reshape(Longitude,1,[]));
y=double(reshape(Latitude,1,[]));

x11=double(reshape(Longitude,1,[]));
y11=double(reshape(Latitude,1,[]));

z=double(reshape(ObserveData,1,[]));
l=double(reshape(LandSeaMask,1,[]));
figure;worldmap world;geoshow(y11,x11,'DisplayType','point');
xInt=round(x);
yInt=round(y);
ObserveXY=[yInt' xInt'];

%%
load gshhs_land_f.mat  %GSHHS coastline truth value


TruthLatitudeInt=round(TruthLatitude);
TruthLongitudeInt=round(TruthLongitude);
%figure;worldmap world;
hold on;geoshow(TruthLatitude,TruthLongitude,'Color','blue');

Truth=[TruthLatitudeInt' TruthLongitudeInt'];
[C,~,~]=intersect(ObserveXY,Truth,'rows');
tf=ismember(Truth,C,'rows');
subTF = 1:length(TruthLatitudeInt);
subTF = subTF(tf);
TruthLon=TruthLongitude(subTF);
TruthLat=TruthLatitude(subTF);
figure;scatter(TruthLon,TruthLat,'.b');
tf1=ismember(ObserveXY,C,'rows');
subTF1 = 1:length(x);
subTF1 = subTF1(~tf1);
x1=x;y1=y;z1=z;l1=l;
x1(subTF1)=NaN;%%Set the disjoint places to NaN to reduce the calculation time
y1(subTF1)=NaN;
z1(subTF1)=NaN;
l1(subTF1)=NaN;
hold on;plot(x1,y1,'r+');
figure;worldmap world;geoshow(y1,x1,'DisplayType','point');
hold on;geoshow(TruthLatitude,TruthLongitude,'Color','blue');
x2=reshape(x1,[size(Longitude,1) size(Longitude,2)]);%%Do not count disjoint places
y2=reshape(y1,[size(Latitude,1) size(Latitude,2)]);
z2=reshape(z1,[size(ObserveData,1) size(ObserveData,2)]);
l2=reshape(l1,[size(LandSeaMask,1) size(LandSeaMask,2)]);
count=0;
Information={};
%The truth map is regarded as continuous, and it is not necessary to find all the coastlines. 
%It only needs to fix the latitude (it can be considered that the latitude has not changed in a small range), and interpolate the longitude.
AllCount=0;
InCount=0;

name=['bulePoint_FY3C_',data.Filename(24:32),'.mat'];
%save('bulePoint_FY3C_0405_0249','TruthLon','TruthLat');
save(name,'TruthLon','TruthLat');