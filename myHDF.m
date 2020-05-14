%% ����ΪHDF�ļ�
% HDF.name=dataName;              %�ļ���
% HDF.Point_row_indices=int32(index1(index3));          %�к�
% HDF.Choosen_columns=index2(index3)+INDEX5(index3);    %�к�
% HDF.Closest_latitudes=Closest_point.latitudes;        %�����γ��
% HDF.Closest_longitudes=Closest_point.longitudes;      %����㾭��
% HDF.Inflected_latitudes=Solving_point.Latitude;       %��ֵ��γ��
% HDF.Inflected_longitudes=Solving_point.Longitude;     %��ֵ�㾭��
% hdf_create(HDF)
clc;clear;close all;warning off;
data=hdf5info('FY3C_MWRIA_GBAL_L1_20180725_0137_010KM_MS1.hdf');
Inflected_latitudes=hdf5read(data.GroupHierarchy.Groups.Datasets(20));%20
Inflected_longitudes=hdf5read(data.GroupHierarchy.Groups.Datasets(21));%21

%%
load pipei_0725_0137.mat;

row=hdf5read(data.GroupHierarchy.Groups.Datasets(15));%20
column=hdf5read(data.GroupHierarchy.Groups.Datasets(2));%21

%HDF.name='FY3C_MWRIA_GBAL_L1_20180405_0249_010KM_MS_MS1';              %�ļ���
HDF.name=[data.Filename(1:end-5),'_MS1']              %�ļ���
HDF.Point_row_indices=row;          %�к�
HDF.Choosen_columns=column;    %�к�
if (HDF.name(10)=='A')
    HDF.Closest_latitudes=Y1(:,1)-0.08 %-0.8 %-0.08    ;%+0.08;        %�����γ��
    HDF.Closest_longitudes=Y1(:,2)%+2.2 %+0.2   ;%-0.08;      %����㾭��
else
    HDF.Closest_latitudes=Y1(:,1)+0.08 %+0.8 %+0.08    ;%-0.08;        %�����γ��
    HDF.Closest_longitudes=Y1(:,2)%-2.2 %-0.2    ; %+0.08;      %����㾭��
end
HDF.Inflected_latitudes=Inflected_latitudes;       %��ֵ��γ��
HDF.Inflected_longitudes=Inflected_longitudes;     %��ֵ�㾭��
hdf_create(HDF)
% 
%data2=hdf5info('FY3C_MWRIA_GBAL_L1_20180403_1330_010KM_MS2.hdf');

