%The main function of drawing the local brightness temperature map
%%
clc;clear;close all;warning off;
currentFolder = pwd;
addpath(genpath(currentFolder));

datachar='FY3C_MWRIA_GBAL_L1_20180411_0259_010KM_MS.HDF';

dirName = datachar(20:27);
cdDir = [currentFolder, '\', dirName];

data=hdf5info(datachar);
Earth_Obs_BT_10_To_89=hdf5read(data.GroupHierarchy.Groups(1).Datasets(2));
ObserveData=double(Earth_Obs_BT_10_To_89(:,:,10))*0.01+327.68;
Latitude=double(hdf5read(data.GroupHierarchy.Groups(2).Datasets(1)));
Longitude=double(hdf5read(data.GroupHierarchy.Groups(2).Datasets(2)));


Latitude2=double(hdf5read(data.GroupHierarchy.Groups(4).Datasets(5)));
Longitude2=double(hdf5read(data.GroupHierarchy.Groups(4).Datasets(6)));

name1=[datachar(20:32) 'old'];
name2=[datachar(20:32) 'new'];

Scale=[58, 66, -86, -72]

LatScale=1;%30
LonScale=1;%90

load gshhs_land_f;
Inflected_latitudes1=TruthLatitude;
Inflected_longitudes1=TruthLongitude;
load GSHHS;
Plot_Color(ObserveData,Latitude,Longitude,Latitude2,Longitude2,TruthLatitude,TruthLongitude,-48, -68);
Plot_Color(ObserveData,Latitude2,Longitude2,Inflected_latitudes1,Inflected_longitudes1,TruthLatitude,TruthLongitude,-48, -68);