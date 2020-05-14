%Extract the coordinate information of the detected coastline from the .hdf into .mat
%%
clc;clear;close all;warning off;
data=hdf5info('FY3C_MWRIA_GBAL_L1_20180904_1818_010KM_MS1.hdf');
% a=double(hdf5read(data.GroupHierarchy.Groups(2).Datasets(1)));
% b=double(hdf5read(data.GroupHierarchy.Groups(4).Datasets(5)));

%全球地图
Inflected_latitudes=hdf5read(data.GroupHierarchy.Groups.Datasets(20));%12  6
Inflected_longitudes=hdf5read(data.GroupHierarchy.Groups.Datasets(21));%13  7
% figure;worldmap world;geoshow(Inflected_latitudes,Inflected_longitudes,'DisplayType','point');
% load GSHHS.mat;
% hold on;geoshow(TruthLatitude,TruthLongitude,'Color','blue');
% 
% %%
Cross_Tracks=hdf5read(data.GroupHierarchy.Groups.Datasets(8));
Along_Tracks=hdf5read(data.GroupHierarchy.Groups.Datasets(1));

% Threshold=Cross_Tracks.^2+Along_Tracks.^2>0.01;
% Cross_Tracks(Threshold)=[];Along_Tracks(Threshold)=[];

% Error_latitudes=hdf5read(data.GroupHierarchy.Groups.Datasets(9));
% Error_longitudes=hdf5read(data.GroupHierarchy.Groups.Datasets(10));
% %%查看均值
mean_CrossTracks=mean(Cross_Tracks);
mean_AlongTracks=mean(Along_Tracks);
% 
figure;scatter(Cross_Tracks,Along_Tracks,'b.');hold on;
scatter(mean(Cross_Tracks),mean(Along_Tracks),'r*','LineWidth',1.5);
xlabel('Cross-Track error (degree)');ylabel('In-Track error (degree)');
%title('20140520_0142');
%axis([-0.085 0.085 -0.085 0.085]);
axis([-0.16 0.16 -0.16 0.16]);
grid on;

% load 610150.mat;
% findLat=Y1(:,1);
% findLon=Y1(:,2);
% figure;worldmap world;geoshow(TruthLatitude,TruthLongitude,'Color','blue');
% geoshow(findLat,findLon,'DisplayType','point');

name=['redPoint_',data.Filename(24:32),'.mat'];
%save('redPoint_0405_0249.mat','Inflected_longitudes','Inflected_latitudes');
save(name,'Inflected_longitudes','Inflected_latitudes');