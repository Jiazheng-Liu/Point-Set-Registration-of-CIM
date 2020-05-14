%Long-term sequence data analysis before and after error correction
clc;clear;close all;warning off;

j_Cross_1=[];
j_Along_1=[];
j_Cross_2=[];
j_Along_2=[];
hdfDataPath = pwd;
hdfDataDir = dir(hdfDataPath);
errorCount = 0;
dayCount = 0;
theYear = hdfDataPath(end-3:end);

%Limit date range
minMoon = 1;
maxMoon = 12;
minDay = 1;
maxDay = 31;

tic;
for i = 1:length(hdfDataDir)
    if(isequal(hdfDataDir(i).name,'.')||... % Remove the two hidden folders that come with the system
    isequal(hdfDataDir(i).name,'..')||...
    ~hdfDataDir(i).isdir) % Remove non-folders in traversal
        continue;
    end
   
    scanYear = str2double(hdfDataDir(i).name(1:4));
    scanMoon = str2double(hdfDataDir(i).name(5:6));
    scanDay = str2double(hdfDataDir(i).name(7:8));
    if(scanMoon >= minMoon && scanMoon <= maxMoon && scanDay >= minDay && scanDay <= maxDay)
        
        AlongPixel_ErrorEverydayCount1 = 0;
        AlongPixel_ErrorEverydayCount2 = 0;
        CrossPixel_ErrorEverydayCount1 = 0;
        CrossPixel_ErrorEverydayCount2 = 0;
        AlongPixel_ErrorEverydayCount_Var1 = [];
        AlongPixel_ErrorEverydayCount_Var2 = [];
        CrossPixel_ErrorEverydayCount_Var1 = [];
        CrossPixel_ErrorEverydayCount_Var2 = [];
        AlongPixel_ErrorEverydayCount_NoAdd1 = 0;
        AlongPixel_ErrorEverydayCount_NoAdd2 = 0;
        CrossPixel_ErrorEverydayCount_NoAdd1 = 0;
        CrossPixel_ErrorEverydayCount_NoAdd2 = 0;
        AlongPixel_ErrorEverydayCount_VarNoAdd1 = [];
        AlongPixel_ErrorEverydayCount_VarNoAdd2 = [];
        CrossPixel_ErrorEverydayCount_VarNoAdd1 = [];
        CrossPixel_ErrorEverydayCount_VarNoAdd2 = [];
        numberHDF_Everyday = 0;
        
        newDirName = fullfile(hdfDataPath,'\',hdfDataDir(i).name);
        cd (newDirName);%addpath(newDirName);
        hdfPath = fullfile(newDirName, '*010KM_MS.hdf');
        hdfDir = dir(hdfPath);
        
        %重新排序
        for j0 = 1:length(hdfDir)
            NoSortName = hdfDir(j0).name;
            SortName = [NoSortName(19:end-4), NoSortName(1:18), '.HDF'];
            movefile(NoSortName, SortName);
        end
         SortNameDir = dir(fullfile(newDirName, '*L1.hdf'));
        
        for j =1:length(SortNameDir)
            NewSortName = [SortNameDir(j).name(24:end-4), SortNameDir(j).name(1:23), '.HDF']
            hdfName1_inDir = [NewSortName(1:end-4), '2.HDF'];
            hdfName2_inDir = [NewSortName(1:end-4), '3.HDF'];
            
            hdfName_inDir_MS1 = [NewSortName(1:end-4), '1.HDF'];
            data_MS1 = hdf5info(hdfName_inDir_MS1);
            
            if exist(hdfName1_inDir, 'file')
               
                numberHDF_Everyday = numberHDF_Everyday + 1;
                errorCount = errorCount + 1;%index
                
                data = hdf5info(SortNameDir(j).name);
                DataLatitude1=double(hdf5read(data.GroupHierarchy.Groups(2).Datasets(1)));
                DataLongitude1=double(hdf5read(data.GroupHierarchy.Groups(2).Datasets(2)));
                DataLatitude2=double(hdf5read(data.GroupHierarchy.Groups(4).Datasets(7)));
                DataLongitude2=double(hdf5read(data.GroupHierarchy.Groups(4).Datasets(8)));
                AttiOff_try1 = hdf5read(data.GroupHierarchy.Groups(4).Datasets(1));
                AttiOff_try2 = hdf5read(data.GroupHierarchy.Groups(4).Datasets(2));
                if (size(AttiOff_try2,1)==3 && size(AttiOff_try2,2)==1)
                    DataLatitude2=double(hdf5read(data.GroupHierarchy.Groups(4).Datasets(7)));
                    DataLongitude2=double(hdf5read(data.GroupHierarchy.Groups(4).Datasets(8)));
                elseif (size(AttiOff_try1,1)==3 && size(AttiOff_try1,2)==1)
                    DataLatitude2=double(hdf5read(data.GroupHierarchy.Groups(4).Datasets(5)));
                    DataLongitude2=double(hdf5read(data.GroupHierarchy.Groups(4).Datasets(6)));
                end
            
                %Read 010KM_MS1.hdf
                data1 = hdf5info(hdfName1_inDir);
                AorD1=hdfName1_inDir(10);
                closestLat_withOutSys1=hdf5read(data1.GroupHierarchy.Groups.Datasets(2));%20
                closestLon_withOutSys1=hdf5read(data1.GroupHierarchy.Groups.Datasets(3));%21
                Inflected_latitudes1=hdf5read(data1.GroupHierarchy.Groups.Datasets(4));%12
                Inflected_longitudes1=hdf5read(data1.GroupHierarchy.Groups.Datasets(5));%13
                if (AorD1=='A')
                    
                    Error_latitudes1=closestLat_withOutSys1 - Inflected_latitudes1;
                    Error_longitudes1=closestLon_withOutSys1 - Inflected_longitudes1;
                    
                elseif (AorD1=='D')
                    Error_latitudes1=closestLat_withOutSys1 - Inflected_latitudes1;
                    Error_longitudes1=closestLon_withOutSys1 - Inflected_longitudes1; 
                end

                    num_hang=hdf5read(data1.GroupHierarchy.Groups.Datasets(6));
                    
                   Latitude=double(hdf5read(data.GroupHierarchy.Groups(2).Datasets(1)));%latitude
                   Longitude=double(hdf5read(data.GroupHierarchy.Groups(2).Datasets(2)));%longitude
                 LatCenter=(Latitude(127,:)+Latitude(128,:))/2;
                 LonCenter=(Longitude(127,:)+Longitude(128,:))/2;
                 TrackDir=[diff(LatCenter'),diff(LonCenter')];
                 HeadingAngle=acos(abs(TrackDir(:,1))./sqrt(TrackDir(:,1).^2+TrackDir(:,2).^2));
                 HeadingAngle=[HeadingAngle;HeadingAngle(end)];

                 Heading_angles1=HeadingAngle(num_hang);
                 Heading_angles11=(hdf5read(data_MS1.GroupHierarchy.Groups.Datasets(11)))*pi/180;

                Along_tracks1=-Error_longitudes1.*cos(Heading_angles1)+(Error_latitudes1).*sin(Heading_angles1);
                Cross_tracks1=Error_longitudes1.*sin(Heading_angles1)+Error_latitudes1.*cos(Heading_angles1);
                Along_tracks_MS1=hdf5read(data_MS1.GroupHierarchy.Groups.Datasets(1));
                Cross_tracks_MS1=hdf5read(data_MS1.GroupHierarchy.Groups.Datasets(8));
                Choosen_columns1=hdf5read(data1.GroupHierarchy.Groups.Datasets(1));%2
                Point_row_indices1=hdf5read(data1.GroupHierarchy.Groups.Datasets(6));%15
                
     
                NaNRegion1 = isnan(Along_tracks1)|isnan(Cross_tracks1)|isnan(Along_tracks_MS1)|isnan(Cross_tracks_MS1);
                Along_tracks1=Along_tracks1(find(NaNRegion1==0));
                Cross_tracks1=Cross_tracks1(find(NaNRegion1==0));
                Along_tracks_MS1=Along_tracks_MS1(find(NaNRegion1==0));
                Cross_tracks_MS1=Cross_tracks_MS1(find(NaNRegion1==0));
                Choosen_columns1=Choosen_columns1(find(NaNRegion1==0));
                Point_row_indices1=Point_row_indices1(find(NaNRegion1==0));
                
                Along_tracks1=Along_tracks1((find(isnan(Cross_tracks1)==0)));
                Along_tracks_MS1=Along_tracks_MS1((find(isnan(Cross_tracks1)==0)));
                Cross_tracks_MS1=Cross_tracks_MS1((find(isnan(Cross_tracks1)==0)));
                Choosen_columns1=Choosen_columns1((find(isnan(Cross_tracks1)==0)));
                Point_row_indices1=Point_row_indices1((find(isnan(Cross_tracks1)==0)));
                Cross_tracks1=Cross_tracks1((find(isnan(Cross_tracks1)==0)));
                 
                
                for j1=1:length(Along_tracks1)
                    currentLat=DataLatitude1(floor(Choosen_columns1(j1))+1,Point_row_indices1(j1)+1);
                    currentLon=DataLongitude1(floor(Choosen_columns1(j1))+1,Point_row_indices1(j1)+1);

                    if(floor(Choosen_columns1(j1))+1 == size(DataLatitude1,1) || Point_row_indices1(j1)+1 == size(DataLatitude1,2))
                        AlongPoint_Lat=DataLatitude1(floor(Choosen_columns1(j1))+1,Point_row_indices1(j1)+1-1);
                        AlongPoint_Lon=DataLongitude1(floor(Choosen_columns1(j1))+1,Point_row_indices1(j1)+1-1);
                        CrossPoint_Lat=DataLatitude1(floor(Choosen_columns1(j1))+1-1,Point_row_indices1(j1)+1);
                        CrossPoint_Lon=DataLongitude1(floor(Choosen_columns1(j1))+1-1,Point_row_indices1(j1)+1);
                    else
                        AlongPoint_Lat=DataLatitude1(floor(Choosen_columns1(j1))+1,Point_row_indices1(j1)+1+1);
                        AlongPoint_Lon=DataLongitude1(floor(Choosen_columns1(j1))+1,Point_row_indices1(j1)+1+1);
                        CrossPoint_Lat=DataLatitude1(floor(Choosen_columns1(j1))+1+1,Point_row_indices1(j1)+1);
                        CrossPoint_Lon=DataLongitude1(floor(Choosen_columns1(j1))+1+1,Point_row_indices1(j1)+1);
                    end
                    CrossPixel = sqrt((CrossPoint_Lon-currentLon).^2 + (CrossPoint_Lat-currentLat).^2);
                    AlongPixel = sqrt((AlongPoint_Lon-currentLon).^2 + (AlongPoint_Lat-currentLat).^2);
                    if (isnan(CrossPixel) || isnan(AlongPixel))
                        continue;
                    end
                    
                    CrossPixel_Error1(j1,:) = Cross_tracks1(j1)./CrossPixel;
                    AlongPixel_Error1(j1,:) = Along_tracks1(j1)./AlongPixel;
                    
                    CrossPixel_Error_NoAdd1(j1,:) = Cross_tracks_MS1(j1)./CrossPixel;
                    AlongPixel_Error_NoAdd1(j1,:) = Along_tracks_MS1(j1)./AlongPixel;
                end

                CrossPixel_ErrorEveryHDF1(errorCount,:) = abs(mean(CrossPixel_Error1));
                AlongPixel_ErrorEveryHDF1(errorCount,:) = abs(mean(AlongPixel_Error1));
                
                Cross_std_every_1(errorCount,:)=std(CrossPixel_Error1);
                Along_std_every_1(errorCount,:)=std(AlongPixel_Error1);
                
                Cross_mean_every_1(errorCount,:)=mean(abs(CrossPixel_Error1));
                Along_mean_every_1(errorCount,:)=mean(abs(AlongPixel_Error1));

                               
                if (i==24||i==25||i==26)%Take a small amount of data to test
                    j_Cross_1=[j_Cross_1;CrossPixel_Error1];
                    j_Along_1=[j_Along_1;AlongPixel_Error1];
                end
                
                
                %Non-absolute pixel error
                CrossPixel_ErrorEveryHDF_NoAbs1(errorCount,:) = mean(CrossPixel_Error1);
                AlongPixel_ErrorEveryHDF_NoAbs1(errorCount,:) = mean(AlongPixel_Error1);
                %Variance of each HDF, used to draw histogram
                CrossPixel_ErrorEveryHDF_Var1(errorCount,:) = var(Cross_tracks1/0.06);
                AlongPixel_ErrorEveryHDF_Var1(errorCount,:) = var(Along_tracks1/0.1);
                
                CrossPixel_ErrorEverydayCount1 = CrossPixel_ErrorEverydayCount1 + CrossPixel_ErrorEveryHDF1(errorCount,:);
                AlongPixel_ErrorEverydayCount1 = AlongPixel_ErrorEverydayCount1 + AlongPixel_ErrorEveryHDF1(errorCount,:);
                
                CrossPixel_ErrorEverydayCount_Var1 = [CrossPixel_ErrorEverydayCount_Var1, CrossPixel_ErrorEveryHDF1(errorCount,:)];
                AlongPixel_ErrorEverydayCount_Var1 = [AlongPixel_ErrorEverydayCount_Var1, AlongPixel_ErrorEveryHDF1(errorCount,:)];
%37                
%                 
                CrossPixel_ErrorEveryHDF_NoAdd1(errorCount,:) = abs(mean(CrossPixel_Error_NoAdd1));
                AlongPixel_ErrorEveryHDF_NoAdd1(errorCount,:) = abs(mean(AlongPixel_Error_NoAdd1));
%                 
                CrossPixel_ErrorEveryHDF_NoAbsNoAdd1(errorCount,:) = mean(CrossPixel_Error_NoAdd1);
                AlongPixel_ErrorEveryHDF_NoAbsNoAdd1(errorCount,:) = mean(AlongPixel_Error_NoAdd1);
               
                CrossPixel_ErrorEveryHDF_VarNoAdd1(errorCount,:) = var(Cross_tracks_MS1/0.06);
                AlongPixel_ErrorEveryHDF_VarNoAdd1(errorCount,:) = var(Along_tracks_MS1/0.1);
                
                CrossPixel_ErrorEverydayCount_NoAdd1 = CrossPixel_ErrorEverydayCount_NoAdd1 + CrossPixel_ErrorEveryHDF_NoAdd1(errorCount,:);
                AlongPixel_ErrorEverydayCount_NoAdd1 = AlongPixel_ErrorEverydayCount_NoAdd1 + AlongPixel_ErrorEveryHDF_NoAdd1(errorCount,:);
               
                CrossPixel_ErrorEverydayCount_VarNoAdd1 = [CrossPixel_ErrorEverydayCount_VarNoAdd1, CrossPixel_ErrorEveryHDF_NoAdd1(errorCount,:)];
                AlongPixel_ErrorEverydayCount_VarNoAdd1 = [AlongPixel_ErrorEverydayCount_VarNoAdd1, AlongPixel_ErrorEveryHDF_NoAdd1(errorCount,:)];
                
                % Read 010KM_MS2.hdf
                data2 = hdf5info(hdfName2_inDir);
                AorD2=hdfName2_inDir(10);
                closestLat_withOutSys2=hdf5read(data2.GroupHierarchy.Groups.Datasets(20));
                closestLon_withOutSys2=hdf5read(data2.GroupHierarchy.Groups.Datasets(21));
                Inflected_latitudes2=hdf5read(data2.GroupHierarchy.Groups.Datasets(12));
                Inflected_longitudes2=hdf5read(data2.GroupHierarchy.Groups.Datasets(13));
                if (AorD2=='A')
                    Error_latitudes2=closestLat_withOutSys2 - Inflected_latitudes2 ;
                    Error_longitudes2=closestLon_withOutSys2 - Inflected_longitudes2 ;
                elseif (AorD2=='D')
                    Error_latitudes2=closestLat_withOutSys2 - Inflected_latitudes2 ;
                    Error_longitudes2=closestLon_withOutSys2 - Inflected_longitudes2 ;
                end

                   %Calculate new solar angle
                  num_hang=hdf5read(data2.GroupHierarchy.Groups.Datasets(15));
                    
                   Latitude=double(hdf5read(data.GroupHierarchy.Groups(4).Datasets(5)));%latitude
                   Longitude=double(hdf5read(data.GroupHierarchy.Groups(4).Datasets(6)));%longitude
                 LatCenter=(Latitude(127,:)+Latitude(128,:))/2;
                 LonCenter=(Longitude(127,:)+Longitude(128,:))/2;
                 TrackDir=[diff(LatCenter'),diff(LonCenter')];
                 HeadingAngle=acos(abs(TrackDir(:,1))./sqrt(TrackDir(:,1).^2+TrackDir(:,2).^2));
                 HeadingAngle=[HeadingAngle;HeadingAngle(end)];
                %HeadingAngle1=HeadingAngle*180/pi;
                 Heading_angles2=HeadingAngle(num_hang);

                Heading_angles22=(hdf5read(data2.GroupHierarchy.Groups.Datasets(11)))*pi/180;
                Along_tracks2=-Error_longitudes2.*cos(Heading_angles2)+(Error_latitudes2).*sin(Heading_angles2);
                Cross_tracks2=Error_longitudes2.*sin(Heading_angles2)+Error_latitudes2.*cos(Heading_angles2);
                Along_tracks_MS2 = hdf5read(data2.GroupHierarchy.Groups.Datasets(1));
                Cross_tracks_MS2 = hdf5read(data2.GroupHierarchy.Groups.Datasets(8));
                Choosen_columns2=hdf5read(data2.GroupHierarchy.Groups.Datasets(2));%254
                Point_row_indices2=hdf5read(data2.GroupHierarchy.Groups.Datasets(15));%17XX
                
                NaNRegion2 = isnan(Along_tracks2)|isnan(Cross_tracks2)|isnan(Along_tracks_MS2)|isnan(Cross_tracks_MS2);
                Along_tracks2=Along_tracks2(find(NaNRegion2==0));
                Cross_tracks2=Cross_tracks2(find(NaNRegion2==0));
                Along_tracks_MS2=Along_tracks_MS2(find(NaNRegion2==0));
                Cross_tracks_MS2=Cross_tracks_MS2(find(NaNRegion2==0));
                Choosen_columns2=Choosen_columns2(find(NaNRegion2==0));
                Point_row_indices2=Point_row_indices2(find(NaNRegion2==0));
                
                
                for j1=1:length(Along_tracks2)
                    currentLat=DataLatitude2(floor(Choosen_columns2(j1))+1,Point_row_indices2(j1)+1);
                    currentLon=DataLongitude2(floor(Choosen_columns2(j1))+1,Point_row_indices2(j1)+1);

                    if(floor(Choosen_columns2(j1))+1 == size(DataLatitude2,1) || Point_row_indices2(j1)+1 == size(DataLatitude2,2))
                        AlongPoint_Lat=DataLatitude2(floor(Choosen_columns2(j1))+1,Point_row_indices2(j1)+1-1);
                        AlongPoint_Lon=DataLongitude2(floor(Choosen_columns2(j1))+1,Point_row_indices2(j1)+1-1);
                        CrossPoint_Lat=DataLatitude2(floor(Choosen_columns2(j1))+1-1,Point_row_indices2(j1)+1);
                        CrossPoint_Lon=DataLongitude2(floor(Choosen_columns2(j1))+1-1,Point_row_indices2(j1)+1);
                    else
                        AlongPoint_Lat=DataLatitude2(floor(Choosen_columns2(j1))+1,Point_row_indices2(j1)+1+1);
                        AlongPoint_Lon=DataLongitude2(floor(Choosen_columns2(j1))+1,Point_row_indices2(j1)+1+1);
                        CrossPoint_Lat=DataLatitude2(floor(Choosen_columns2(j1))+1+1,Point_row_indices2(j1)+1);
                        CrossPoint_Lon=DataLongitude2(floor(Choosen_columns2(j1))+1+1,Point_row_indices2(j1)+1);
                    end
                    CrossPixel = sqrt((CrossPoint_Lon-currentLon).^2 + (CrossPoint_Lat-currentLat).^2);
                    AlongPixel = sqrt((AlongPoint_Lon-currentLon).^2 + (AlongPoint_Lat-currentLat).^2);
                    if (isnan(CrossPixel) || isnan(AlongPixel))
                        continue;
                    end
                
                    CrossPixel_Error2(j1,:) = Cross_tracks2(j1)./CrossPixel;
                    AlongPixel_Error2(j1,:) = Along_tracks2(j1)./AlongPixel;
                    
                    CrossPixel_Error_NoAdd2(j1,:) = Cross_tracks_MS2(j1)./CrossPixel;
                    AlongPixel_Error_NoAdd2(j1,:) = Along_tracks_MS2(j1)./AlongPixel;
                end

                CrossPixel_ErrorEveryHDF2(errorCount,:) = abs(mean(CrossPixel_Error2));
                AlongPixel_ErrorEveryHDF2(errorCount,:) = abs(mean(AlongPixel_Error2));
                Cross_std_every_2(errorCount,:)=std(CrossPixel_Error2);
                Along_std_every_2(errorCount,:)=std(AlongPixel_Error2);
                
                Cross_mean_every_2(errorCount,:)=mean(abs(CrossPixel_Error2));
                Along_mean_every_2(errorCount,:)=mean(abs(AlongPixel_Error2));
                

                if (i==24||i==25||i==26)
                   j_Cross_2=[j_Cross_2;CrossPixel_Error2];
                   j_Along_2=[j_Along_2;AlongPixel_Error2];
                end
                
                %Non-absolute
                CrossPixel_ErrorEveryHDF_NoAbs2(errorCount,:) = mean(CrossPixel_Error2);
                AlongPixel_ErrorEveryHDF_NoAbs2(errorCount,:) = mean(AlongPixel_Error2);
                %Variance of each HDF, used to draw histogram
                CrossPixel_ErrorEveryHDF_Var2(errorCount,:) = var(Cross_tracks2/0.06);
                AlongPixel_ErrorEveryHDF_Var2(errorCount,:) = var(Along_tracks2/0.1);

                CrossPixel_ErrorEverydayCount2 = CrossPixel_ErrorEverydayCount2 + CrossPixel_ErrorEveryHDF2(errorCount,:);
                AlongPixel_ErrorEverydayCount2 = AlongPixel_ErrorEverydayCount2 + AlongPixel_ErrorEveryHDF2(errorCount,:);
                
                CrossPixel_ErrorEverydayCount_Var2 = [CrossPixel_ErrorEverydayCount_Var2, CrossPixel_ErrorEveryHDF2(errorCount,:)];
                AlongPixel_ErrorEverydayCount_Var2 = [AlongPixel_ErrorEverydayCount_Var2, AlongPixel_ErrorEveryHDF2(errorCount,:)];
                

                CrossPixel_ErrorEveryHDF_NoAdd2(errorCount,:) = abs(mean(CrossPixel_Error_NoAdd2));
                AlongPixel_ErrorEveryHDF_NoAdd2(errorCount,:) = abs(mean(AlongPixel_Error_NoAdd2));

                CrossPixel_ErrorEveryHDF_NoAbsNoAdd2(errorCount,:) = mean(CrossPixel_Error_NoAdd2);
                AlongPixel_ErrorEveryHDF_NoAbsNoAdd2(errorCount,:) = mean(AlongPixel_Error_NoAdd2);
         
                CrossPixel_ErrorEveryHDF_VarNoAdd2(errorCount,:) = var(Cross_tracks_MS2/0.06);
                AlongPixel_ErrorEveryHDF_VarNoAdd2(errorCount,:) = var(Along_tracks_MS2/0.1);
        
                CrossPixel_ErrorEverydayCount_NoAdd2 = CrossPixel_ErrorEverydayCount_NoAdd2 + CrossPixel_ErrorEveryHDF_NoAdd2(errorCount,:);
                AlongPixel_ErrorEverydayCount_NoAdd2 = AlongPixel_ErrorEverydayCount_NoAdd2 + AlongPixel_ErrorEveryHDF_NoAdd2(errorCount,:);
        
                CrossPixel_ErrorEverydayCount_VarNoAdd2 = [CrossPixel_ErrorEverydayCount_VarNoAdd2, CrossPixel_ErrorEveryHDF_NoAdd2(errorCount,:)];
                AlongPixel_ErrorEverydayCount_VarNoAdd2 = [AlongPixel_ErrorEverydayCount_VarNoAdd2, AlongPixel_ErrorEveryHDF_NoAdd2(errorCount,:)];
            end
        end
        
        dayCount = dayCount+1;
        CrossPixel_ErrorEveryday1(dayCount,:)=CrossPixel_ErrorEverydayCount1/numberHDF_Everyday;
        AlongPixel_ErrorEveryday1(dayCount,:)=AlongPixel_ErrorEverydayCount1/numberHDF_Everyday;
        CrossPixel_ErrorEveryday2(dayCount,:)=CrossPixel_ErrorEverydayCount2/numberHDF_Everyday;
        AlongPixel_ErrorEveryday2(dayCount,:)=AlongPixel_ErrorEverydayCount2/numberHDF_Everyday;
        
        CrossPixel_ErrorEveryday_Var1(dayCount,:)=var(CrossPixel_ErrorEverydayCount_Var1);
        AlongPixel_ErrorEveryday_Var1(dayCount,:)=var(AlongPixel_ErrorEverydayCount_Var1);
        CrossPixel_ErrorEveryday_Var2(dayCount,:)=var(CrossPixel_ErrorEverydayCount_Var2);
        AlongPixel_ErrorEveryday_Var2(dayCount,:)=var(AlongPixel_ErrorEverydayCount_Var2);
        
   
        CrossPixel_ErrorEveryday_NoAdd1(dayCount,:)=CrossPixel_ErrorEverydayCount_NoAdd1/numberHDF_Everyday;
        AlongPixel_ErrorEveryday_NoAdd1(dayCount,:)=AlongPixel_ErrorEverydayCount_NoAdd1/numberHDF_Everyday;
        CrossPixel_ErrorEveryday_NoAdd2(dayCount,:)=CrossPixel_ErrorEverydayCount_NoAdd2/numberHDF_Everyday;
        AlongPixel_ErrorEveryday_NoAdd2(dayCount,:)=AlongPixel_ErrorEverydayCount_NoAdd2/numberHDF_Everyday;
   
        CrossPixel_ErrorEveryday_VarNoAdd1(dayCount,:)=var(CrossPixel_ErrorEverydayCount_VarNoAdd1);
        AlongPixel_ErrorEveryday_VarNoAdd1(dayCount,:)=var(AlongPixel_ErrorEverydayCount_VarNoAdd1);
        CrossPixel_ErrorEveryday_VarNoAdd2(dayCount,:)=var(CrossPixel_ErrorEverydayCount_VarNoAdd2);
        AlongPixel_ErrorEveryday_VarNoAdd2(dayCount,:)=var(AlongPixel_ErrorEverydayCount_VarNoAdd2);
        
        %Sort
        for j1 = 1:length(SortNameDir)
            NoSortName1 = SortNameDir(j1).name;
            SortName1 = [NoSortName1(24:end-4), NoSortName1(1:23), '.HDF'];
            movefile(NoSortName1, SortName1);
        end
        cd ..;
    end
end
toc;

plotNum = 1:size(AlongPixel_ErrorEveryHDF1,1);%errorCount
plotDay = 1:size(AlongPixel_ErrorEveryday1,1);%dayCount
linewidth0 = 0.5;
picName = [theYear,'-',num2str(minMoon),'月至',theYear,'-',num2str(maxMoon),'月'];



%%
%Every HDF along the track
MeanAlongPixel_ErrorEveryHDF1 = mean(AlongPixel_ErrorEveryHDF1);
MeanAlongPixel_ErrorEveryHDF2 = mean(AlongPixel_ErrorEveryHDF2);
VarAlongPixel_ErrorEveryHDF1 = var(AlongPixel_ErrorEveryHDF1);
VarAlongPixel_ErrorEveryHDF2 = var(AlongPixel_ErrorEveryHDF2);
axisYMax1 = max(max(AlongPixel_ErrorEveryHDF1),max(AlongPixel_ErrorEveryHDF2))+0.1;

figure;plot(plotNum,AlongPixel_ErrorEveryHDF1,'r','linewidth',linewidth0);
hold on;plot(plotNum,AlongPixel_ErrorEveryHDF2,'b','linewidth',linewidth0);grid on;
Tagging1=legend('修正前','修正后');set(Tagging1,'Fontsize',20);
% Tagging1=legend('Before correction','After correction');set(Tagging1,'Fontsize',20);
text(1,axisYMax1-1*(axisYMax1/40),['修正前均值：' num2str(MeanAlongPixel_ErrorEveryHDF1)],'FontSize',15);
text(1,axisYMax1-2*(axisYMax1/40),['修正后均值：' num2str(MeanAlongPixel_ErrorEveryHDF2)],'FontSize',15);
text(1,axisYMax1-3*(axisYMax1/40),['修正前方差：' num2str(VarAlongPixel_ErrorEveryHDF1)],'FontSize',15);
text(1,axisYMax1-4*(axisYMax1/40),['修正后方差：' num2str(VarAlongPixel_ErrorEveryHDF2)],'FontSize',15);
xlabel([picName, '（每个文件）'],'FontSize',20);
ylabel('沿轨误差（像素）','FontSize',20);
% ylabel('Along-track error（Pixel）','FontSize',20);
axis([0 length(plotNum)+1 0 axisYMax1]);
set(gcf, 'outerposition', get(0, 'screensize'));
saveas(gcf, [picName, '-沿轨（每个HDF）.tif']);

%Every HDF along the track-no absolute value
MeanAlongPixel_ErrorEveryHDF1_NoAbs = mean(AlongPixel_ErrorEveryHDF_NoAbs1);
MeanAlongPixel_ErrorEveryHDF2_NoAbs = mean(AlongPixel_ErrorEveryHDF_NoAbs2);
VarAlongPixel_ErrorEveryHDF1_NoAbs = var(AlongPixel_ErrorEveryHDF_NoAbs1);
VarAlongPixel_ErrorEveryHDF2_NoAbs = var(AlongPixel_ErrorEveryHDF_NoAbs2);
axisYMax2 = max(max(AlongPixel_ErrorEveryHDF_NoAbs1),max(AlongPixel_ErrorEveryHDF_NoAbs2))+0.3;

figure;plot(plotNum,AlongPixel_ErrorEveryHDF_NoAbs1,'r','linewidth',linewidth0);
hold on;plot(plotNum,AlongPixel_ErrorEveryHDF_NoAbs2,'b','linewidth',linewidth0);grid on;
Tagging1_2=legend('修正前','修正后');set(Tagging1_2,'Fontsize',20);
text(1,axisYMax2-1*(axisYMax2/40),['修正前均值：' num2str(MeanAlongPixel_ErrorEveryHDF1_NoAbs)],'FontSize',15);
text(1,axisYMax2-3*(axisYMax2/40),['修正后均值：' num2str(MeanAlongPixel_ErrorEveryHDF2_NoAbs)],'FontSize',15);
text(1,axisYMax2-5*(axisYMax2/40),['修正前方差：' num2str(VarAlongPixel_ErrorEveryHDF1_NoAbs)],'FontSize',15);
text(1,axisYMax2-7*(axisYMax2/40),['修正后方差：' num2str(VarAlongPixel_ErrorEveryHDF2_NoAbs)],'FontSize',15);
xlabel([picName, '（每个文件）'],'FontSize',20);
ylabel('沿轨误差（像素）','FontSize',20);
axis([0 length(plotNum)+1 -axisYMax2 axisYMax2]);
set(gcf, 'outerposition', get(0, 'screensize'));
saveas(gcf, [picName, '-沿轨（每个HDF-不绝对值）.tif']);

%Along the track-every day
figure;plot(plotDay,AlongPixel_ErrorEveryday1,'r','linewidth',linewidth0);
hold on;plot(plotDay,AlongPixel_ErrorEveryday2,'b','linewidth',linewidth0);grid on;
Tagging2=legend('修正前','修正后');set(Tagging2,'Fontsize',10);
xlabel([picName, '（天）'],'FontSize',13);
ylabel('沿轨误差（像素）','FontSize',13);
axis([0 length(plotDay)+1 0 max(max(AlongPixel_ErrorEveryday1),max(AlongPixel_ErrorEveryday2))+0.03]);
saveas(gcf, [picName, '-沿轨（每整天）.tif']);

%%

% Each HDF along the track-no fixed difference
MeanAlongPixel_ErrorEveryHDF_NoAdd1 = mean(AlongPixel_ErrorEveryHDF_NoAdd1);
MeanAlongPixel_ErrorEveryHDF_NoAdd2 = mean(AlongPixel_ErrorEveryHDF_NoAdd2);
VarAlongPixel_ErrorEveryHDF_NoAdd1 = var(AlongPixel_ErrorEveryHDF_NoAdd1);
VarAlongPixel_ErrorEveryHDF_NoAdd2 = var(AlongPixel_ErrorEveryHDF_NoAdd2);
axisYMax_NoAdd1 = max(max(AlongPixel_ErrorEveryHDF_NoAdd1),max(AlongPixel_ErrorEveryHDF_NoAdd2))+0.03;

figure;plot(plotNum,AlongPixel_ErrorEveryHDF_NoAdd1,'r','linewidth',linewidth0);
hold on;plot(plotNum,AlongPixel_ErrorEveryHDF_NoAdd2,'b','linewidth',linewidth0);grid on;
Tagging_NoAdd1=legend('修正前','修正后');set(Tagging_NoAdd1,'Fontsize',20);
text(1,axisYMax_NoAdd1-1*(axisYMax_NoAdd1/40),['修正前均值：' num2str(MeanAlongPixel_ErrorEveryHDF_NoAdd1)],'FontSize',15);
text(1,axisYMax_NoAdd1-2*(axisYMax_NoAdd1/40),['修正后均值：' num2str(MeanAlongPixel_ErrorEveryHDF_NoAdd2)],'FontSize',15);
text(1,axisYMax_NoAdd1-3*(axisYMax_NoAdd1/40),['修正前方差：' num2str(VarAlongPixel_ErrorEveryHDF_NoAdd1)],'FontSize',15);
text(1,axisYMax_NoAdd1-4*(axisYMax_NoAdd1/40),['修正后方差：' num2str(VarAlongPixel_ErrorEveryHDF_NoAdd2)],'FontSize',15);
xlabel([picName, '（每个文件）'],'FontSize',20);
ylabel('沿轨误差（像素）','FontSize',20);
axis([0 length(plotNum)+1 0 axisYMax_NoAdd1]);
set(gcf, 'outerposition', get(0, 'screensize'));
saveas(gcf, [picName, '-沿轨（每个HDF-无固定差）.tif']);

% Each HDF along the track-no absolute value-and no fixed difference
MeanAlongPixel_ErrorEveryHDF1_NoAbsNoAdd = mean(AlongPixel_ErrorEveryHDF_NoAbsNoAdd1);
MeanAlongPixel_ErrorEveryHDF2_NoAbsNoAdd = mean(AlongPixel_ErrorEveryHDF_NoAbsNoAdd2);
VarAlongPixel_ErrorEveryHDF1_NoAbsNoAdd = var(AlongPixel_ErrorEveryHDF_NoAbsNoAdd1);
VarAlongPixel_ErrorEveryHDF2_NoAbsNoAdd = var(AlongPixel_ErrorEveryHDF_NoAbsNoAdd2);
axisYMax_NoAdd2 = max(max(AlongPixel_ErrorEveryHDF_NoAbsNoAdd1),max(AlongPixel_ErrorEveryHDF_NoAbsNoAdd2))+0.3;

figure;plot(plotNum,AlongPixel_ErrorEveryHDF_NoAbsNoAdd1,'r','linewidth',linewidth0);
hold on;plot(plotNum,AlongPixel_ErrorEveryHDF_NoAbsNoAdd2,'b','linewidth',linewidth0);grid on;
Tagging1_2NoAdd=legend('修正前','修正后');set(Tagging1_2NoAdd,'Fontsize',20);
text(1,axisYMax_NoAdd2-1*(axisYMax_NoAdd2/40),['修正前均值：' num2str(MeanAlongPixel_ErrorEveryHDF1_NoAbsNoAdd)],'FontSize',15);
text(1,axisYMax_NoAdd2-3*(axisYMax_NoAdd2/40),['修正后均值：' num2str(MeanAlongPixel_ErrorEveryHDF2_NoAbsNoAdd)],'FontSize',15);
text(1,axisYMax_NoAdd2-5*(axisYMax_NoAdd2/40),['修正前方差：' num2str(VarAlongPixel_ErrorEveryHDF1_NoAbsNoAdd)],'FontSize',15);
text(1,axisYMax_NoAdd2-7*(axisYMax_NoAdd2/40),['修正后方差：' num2str(VarAlongPixel_ErrorEveryHDF2_NoAbsNoAdd)],'FontSize',15);
xlabel([picName, '（每个文件）'],'FontSize',20);
ylabel('沿轨误差（像素）','FontSize',20);
axis([0 length(plotNum)+1 -axisYMax_NoAdd2 axisYMax_NoAdd2]);
set(gcf, 'outerposition', get(0, 'screensize'));
saveas(gcf, [picName, '-沿轨（每个HDF-不绝对值且无固定差）.tif']);

% Along the track every day-no fixed difference
figure;plot(plotDay,AlongPixel_ErrorEveryday_NoAdd1,'r','linewidth',linewidth0);
hold on;plot(plotDay,AlongPixel_ErrorEveryday_NoAdd2,'b','linewidth',linewidth0);grid on;
Tagging2_NoAdd=legend('修正前','修正后');set(Tagging2_NoAdd,'Fontsize',10);
xlabel([picName, '（天）'],'FontSize',13);
ylabel('沿轨误差（像素）','FontSize',13);
axis([0 length(plotDay)+1 0 max(max(AlongPixel_ErrorEveryday_NoAdd1),max(AlongPixel_ErrorEveryday_NoAdd2))+0.03]);
saveas(gcf, [picName, '-沿轨（每整天-无固定差）.tif']);





%%
% Cross track per HDF
MeanCrossPixel_ErrorEveryHDF1 = mean(CrossPixel_ErrorEveryHDF1);
MeanCrossPixel_ErrorEveryHDF2 = mean(CrossPixel_ErrorEveryHDF2);
VarCrossPixel_ErrorEveryHDF1 = var(CrossPixel_ErrorEveryHDF1);
VarCrossPixel_ErrorEveryHDF2 = var(CrossPixel_ErrorEveryHDF2);
axisYMax3 = max(max(CrossPixel_ErrorEveryHDF1),max(CrossPixel_ErrorEveryHDF2))+0.1;

figure;plot(plotNum,CrossPixel_ErrorEveryHDF1,'r','linewidth',linewidth0);
hold on;plot(plotNum,CrossPixel_ErrorEveryHDF2,'b','linewidth',linewidth0);grid on;
Tagging3=legend('修正前','修正后');set(Tagging3,'Fontsize',20);
% Tagging3=legend('Before correction','After correction');set(Tagging3,'Fontsize',20);
text(1,axisYMax3-1*(axisYMax3/40),['修正前均值：' num2str(MeanCrossPixel_ErrorEveryHDF1)],'FontSize',15);
text(1,axisYMax3-2*(axisYMax3/40),['修正后均值：' num2str(MeanCrossPixel_ErrorEveryHDF2)],'FontSize',15);
text(1,axisYMax3-3*(axisYMax3/40),['修正前方差：' num2str(VarCrossPixel_ErrorEveryHDF1)],'FontSize',15);
text(1,axisYMax3-4*(axisYMax3/40),['修正后方差：' num2str(VarCrossPixel_ErrorEveryHDF2)],'FontSize',15);
xlabel([picName, '（每个文件）'],'FontSize',20);
ylabel('跨轨误差（像素）','FontSize',20);
% ylabel('Cross-track（Pixel）','FontSize',20);
axis([0 length(plotNum)+1 0 axisYMax3]);
set(gcf, 'outerposition', get(0, 'screensize'));
saveas(gcf, [picName, '-跨轨（每个HDF）.tif']);

% Each HDF across track-no absolute value
MeanCrossPixel_ErrorEveryHDF1_NoAbs = mean(CrossPixel_ErrorEveryHDF_NoAbs1);
MeanCrossPixel_ErrorEveryHDF2_NoAbs = mean(CrossPixel_ErrorEveryHDF_NoAbs2);
VarCrossPixel_ErrorEveryHDF1_NoAbs = var(CrossPixel_ErrorEveryHDF_NoAbs1);
VarCrossPixel_ErrorEveryHDF2_NoAbs = var(CrossPixel_ErrorEveryHDF_NoAbs2);
axisYMax4 = max(max(CrossPixel_ErrorEveryHDF_NoAbs1),max(CrossPixel_ErrorEveryHDF_NoAbs2))+0.3;

figure;plot(plotNum,CrossPixel_ErrorEveryHDF_NoAbs1,'r','linewidth',linewidth0);
hold on;plot(plotNum,CrossPixel_ErrorEveryHDF_NoAbs2,'b','linewidth',linewidth0);grid on;
Tagging3_1=legend('修正前','修正后');set(Tagging3_1,'Fontsize',20);
text(1,axisYMax4-1*(axisYMax4/40),['修正前均值：' num2str(MeanCrossPixel_ErrorEveryHDF1_NoAbs)],'FontSize',15);
text(1,axisYMax4-3*(axisYMax4/40),['修正后均值：' num2str(MeanCrossPixel_ErrorEveryHDF2_NoAbs)],'FontSize',15);
text(1,axisYMax4-5*(axisYMax4/40),['修正前方差：' num2str(VarCrossPixel_ErrorEveryHDF1_NoAbs)],'FontSize',15);
text(1,axisYMax4-7*(axisYMax4/40),['修正后方差：' num2str(VarCrossPixel_ErrorEveryHDF2_NoAbs)],'FontSize',15);
xlabel([picName, '（每个文件）'],'FontSize',20);
ylabel('跨轨误差（像素）','FontSize',20);
axis([0 length(plotNum)+1 -axisYMax4 axisYMax4]);
set(gcf, 'outerposition', get(0, 'screensize'));
saveas(gcf, [picName, '-跨轨（每个HDF-不绝对值）.tif']);

% Cross track-every day
figure;plot(plotDay,CrossPixel_ErrorEveryday1,'r','linewidth',linewidth0);
hold on;plot(plotDay,CrossPixel_ErrorEveryday2,'b','linewidth',linewidth0);grid on;
Tagging4=legend('修正前','修正后','Location','northwest');set(Tagging4,'Fontsize',10);
xlabel([picName, '（天）'],'FontSize',13);
ylabel('跨轨误差（像素）','FontSize',13);
axis([0 length(plotDay)+1 0 max(max(CrossPixel_ErrorEveryday1),max(CrossPixel_ErrorEveryday2))+0.03]);
saveas(gcf, [picName, '-跨轨（每整天）.tif']);


%%
% Each HDF across track-no fixed difference
MeanCrossPixel_ErrorEveryHDF_NoAdd1 = mean(CrossPixel_ErrorEveryHDF_NoAdd1);
MeanCrossPixel_ErrorEveryHDF_NoAdd2 = mean(CrossPixel_ErrorEveryHDF_NoAdd2);
VarCrossPixel_ErrorEveryHDF_NoAdd1 = var(CrossPixel_ErrorEveryHDF_NoAdd1);
VarCrossPixel_ErrorEveryHDF_NoAdd2 = var(CrossPixel_ErrorEveryHDF_NoAdd2);
axisYMax_NoAdd3 = max(max(CrossPixel_ErrorEveryHDF_NoAdd1),max(CrossPixel_ErrorEveryHDF_NoAdd2))+0.03;

figure;plot(plotNum,CrossPixel_ErrorEveryHDF_NoAdd1,'r','linewidth',linewidth0);
hold on;plot(plotNum,CrossPixel_ErrorEveryHDF_NoAdd2,'b','linewidth',linewidth0);grid on;
Tagging_NoAdd3=legend('修正前','修正后');set(Tagging_NoAdd3,'Fontsize',20);
text(1,axisYMax_NoAdd3-1*(axisYMax_NoAdd3/40),['修正前均值：' num2str(MeanCrossPixel_ErrorEveryHDF_NoAdd1)],'FontSize',15);
text(1,axisYMax_NoAdd3-2*(axisYMax_NoAdd3/40),['修正后均值：' num2str(MeanCrossPixel_ErrorEveryHDF_NoAdd2)],'FontSize',15);
text(1,axisYMax_NoAdd3-3*(axisYMax_NoAdd3/40),['修正前方差：' num2str(VarCrossPixel_ErrorEveryHDF_NoAdd1)],'FontSize',15);
text(1,axisYMax_NoAdd3-4*(axisYMax_NoAdd3/40),['修正后方差：' num2str(VarCrossPixel_ErrorEveryHDF_NoAdd2)],'FontSize',15);
xlabel([picName, '（每个文件）'],'FontSize',20);
ylabel('跨轨误差（像素）','FontSize',20);
axis([0 length(plotNum)+1 0 axisYMax_NoAdd3]);
set(gcf, 'outerposition', get(0, 'screensize'));
saveas(gcf, [picName, '-跨轨（每个HDF-无固定差）.tif']);

% Each HDF across track-no absolute value-and no fixed difference
MeanCrossPixel_ErrorEveryHDF1_NoAbsNoAdd = mean(CrossPixel_ErrorEveryHDF_NoAbsNoAdd1);
MeanCrossPixel_ErrorEveryHDF2_NoAbsNoAdd = mean(CrossPixel_ErrorEveryHDF_NoAbsNoAdd2);
VarCrossPixel_ErrorEveryHDF1_NoAbsNoAdd = var(CrossPixel_ErrorEveryHDF_NoAbsNoAdd1);
VarCrossPixel_ErrorEveryHDF2_NoAbsNoAdd = var(CrossPixel_ErrorEveryHDF_NoAbsNoAdd2);
axisYMax4_NoAdd = max(max(CrossPixel_ErrorEveryHDF_NoAbsNoAdd1),max(CrossPixel_ErrorEveryHDF_NoAbsNoAdd2))+0.3;

figure;plot(plotNum,CrossPixel_ErrorEveryHDF_NoAbsNoAdd1,'r','linewidth',linewidth0);
hold on;plot(plotNum,CrossPixel_ErrorEveryHDF_NoAbsNoAdd2,'b','linewidth',linewidth0);grid on;
Tagging3_1NoAdd=legend('修正前','修正后');set(Tagging3_1NoAdd,'Fontsize',20);
text(1,axisYMax4_NoAdd-1*(axisYMax4_NoAdd/40),['修正前均值：' num2str(MeanCrossPixel_ErrorEveryHDF1_NoAbsNoAdd)],'FontSize',15);
text(1,axisYMax4_NoAdd-3*(axisYMax4_NoAdd/40),['修正后均值：' num2str(MeanCrossPixel_ErrorEveryHDF2_NoAbsNoAdd)],'FontSize',15);
text(1,axisYMax4_NoAdd-5*(axisYMax4_NoAdd/40),['修正前方差：' num2str(VarCrossPixel_ErrorEveryHDF1_NoAbsNoAdd)],'FontSize',15);
text(1,axisYMax4_NoAdd-7*(axisYMax4_NoAdd/40),['修正后方差：' num2str(VarCrossPixel_ErrorEveryHDF2_NoAbsNoAdd)],'FontSize',15);
xlabel([picName, '（每个文件）'],'FontSize',20);
ylabel('跨轨误差（像素）','FontSize',20);
axis([0 length(plotNum)+1 -axisYMax4_NoAdd axisYMax4_NoAdd]);
set(gcf, 'outerposition', get(0, 'screensize'));
saveas(gcf, [picName, '-跨轨（每个HDF-不绝对值且无固定差）.tif']);

% Cross track-every day
figure;plot(plotDay,CrossPixel_ErrorEveryday_NoAdd1,'r','linewidth',linewidth0);
hold on;plot(plotDay,CrossPixel_ErrorEveryday_NoAdd2,'b','linewidth',linewidth0);grid on;
Tagging4_NoAdd=legend('修正前','修正后','Location','northwest');set(Tagging4_NoAdd,'Fontsize',10);
xlabel([picName, '（天）'],'FontSize',13);
ylabel('跨轨误差（像素）','FontSize',13);
axis([0 length(plotDay)+1 0 max(max(CrossPixel_ErrorEveryday_NoAdd1),max(CrossPixel_ErrorEveryday_NoAdd2))+0.03]);
saveas(gcf, [picName, '-跨轨（每整天-无固定差）.tif']);





%%
% Interval sampling to draw histogram
% Along the track
interva1 = 2;%interval
AlongPixel_ErrorEveryHDF_Sampling1 = AlongPixel_ErrorEveryHDF1(2:interva1:end);
AlongPixel_ErrorEveryHDF_Sampling2 = AlongPixel_ErrorEveryHDF2(2:interva1:end);
% figure;BarSize1 = size(AlongPixel_ErrorEveryHDF_Sampling1,1);
% myBar1 = bar(1:BarSize1,[AlongPixel_ErrorEveryHDF_Sampling1, AlongPixel_ErrorEveryHDF_Sampling2]);
% set(myBar1(1), 'facecolor', 'r');set(myBar1(2), 'facecolor', 'b');

Many_VarAlongPixel_ErrorEveryHDF1=AlongPixel_ErrorEveryHDF_Var1(2:interva1:end);
Many_VarAlongPixel_ErrorEveryHDF2=AlongPixel_ErrorEveryHDF_Var2(2:interva1:end);
figure;myBar3 = barwitherr([Many_VarAlongPixel_ErrorEveryHDF1, Many_VarAlongPixel_ErrorEveryHDF2], [AlongPixel_ErrorEveryHDF_Sampling1, AlongPixel_ErrorEveryHDF_Sampling2]);
set(myBar3(1), 'facecolor', 'r');set(myBar3(2), 'facecolor', 'b');
legend('修正前','修正后','Location','northwest');%
xlabel('从所有HDF文件抽取样本','FontSize',13);
ylabel('沿轨误差 (像素)','FontSize',13);
saveas(gcf, [picName, '-沿轨直方图（每个HDF）.tif']);

% Along the track-every day
interva2 = 2;%interval
AlongPixel_ErrorEveryday_Sampling1 = AlongPixel_ErrorEveryday1(1:interva2:end);
AlongPixel_ErrorEveryday_Sampling2 = AlongPixel_ErrorEveryday2(1:interva2:end);


Many_VarAlongPixel_ErrorEveryday1=AlongPixel_ErrorEveryday_Var1(1:interva2:end);
Many_VarAlongPixel_ErrorEveryday2=AlongPixel_ErrorEveryday_Var2(1:interva2:end);
figure;myBar2 = barwitherr([Many_VarAlongPixel_ErrorEveryday1, Many_VarAlongPixel_ErrorEveryday2], [AlongPixel_ErrorEveryday_Sampling1, AlongPixel_ErrorEveryday_Sampling2]);
set(myBar2(1), 'facecolor', 'r');set(myBar2(2), 'facecolor', 'b');
legend('修正前','修正后');%,'Location','southeast'
xlabel('从整天均值抽取样本','FontSize',13);
ylabel('沿轨误差 (像素)','FontSize',13);
saveas(gcf, [picName, '-沿轨直方图（每整天）.tif']);


%%
% Histogram along track-no fixed difference added
interva_NoAdd1 = 2;%interval
AlongPixel_ErrorEveryHDF_Sampling_NoAdd1 = AlongPixel_ErrorEveryHDF_NoAdd1(2:interva_NoAdd1:end);
AlongPixel_ErrorEveryHDF_Sampling_NoAdd2 = AlongPixel_ErrorEveryHDF_NoAdd2(2:interva_NoAdd1:end);
Many_VarAlongPixel_ErrorEveryHDF_NoAdd1=AlongPixel_ErrorEveryHDF_VarNoAdd1(2:interva_NoAdd1:end);
Many_VarAlongPixel_ErrorEveryHDF_NoAdd2=AlongPixel_ErrorEveryHDF_VarNoAdd2(2:interva_NoAdd1:end);
figure;myBar_NoAdd3 = barwitherr([Many_VarAlongPixel_ErrorEveryHDF_NoAdd1, Many_VarAlongPixel_ErrorEveryHDF_NoAdd2], [AlongPixel_ErrorEveryHDF_Sampling_NoAdd1, AlongPixel_ErrorEveryHDF_Sampling_NoAdd2]);
set(myBar_NoAdd3(1), 'facecolor', 'r');set(myBar_NoAdd3(2), 'facecolor', 'b');
legend('修正前','修正后','Location','northwest');%
xlabel('从所有HDF文件抽取样本','FontSize',13);
ylabel('沿轨误差 (像素)','FontSize',13);
saveas(gcf, [picName, '-沿轨直方图（每个HDF-无固定差）.tif']);

% Daily along the track-no fixed difference
interva_NoAdd2 = 2;%interval
AlongPixel_ErrorEveryday_Sampling_NoAdd1 = AlongPixel_ErrorEveryday_NoAdd1(1:interva_NoAdd2:end);
AlongPixel_ErrorEveryday_Sampling_NoAdd2 = AlongPixel_ErrorEveryday_NoAdd2(1:interva_NoAdd2:end);
Many_VarAlongPixel_ErrorEveryday_NoAdd1=AlongPixel_ErrorEveryday_VarNoAdd1(1:interva_NoAdd2:end);
Many_VarAlongPixel_ErrorEveryday_NoAdd2=AlongPixel_ErrorEveryday_VarNoAdd2(1:interva_NoAdd2:end);
figure;myBar_NoAdd2 = barwitherr([Many_VarAlongPixel_ErrorEveryday_NoAdd1, Many_VarAlongPixel_ErrorEveryday_NoAdd2], [AlongPixel_ErrorEveryday_Sampling_NoAdd1, AlongPixel_ErrorEveryday_Sampling_NoAdd2]);
set(myBar_NoAdd2(1), 'facecolor', 'r');set(myBar_NoAdd2(2), 'facecolor', 'b');
legend('修正前','修正后');%,'Location','southeast'
xlabel('从整天均值抽取样本','FontSize',13);
ylabel('沿轨误差 (像素)','FontSize',13);
saveas(gcf, [picName, '-沿轨直方图（每整天-无固定差）.tif']);




%%
% Interval sampling to draw histogram
% Cross track
interva3 = 2;%interval
CrossPixel_ErrorEveryHDF_Sampling1 = CrossPixel_ErrorEveryHDF1(2:interva3:end);
CrossPixel_ErrorEveryHDF_Sampling2 = CrossPixel_ErrorEveryHDF2(2:interva3:end);


Many_VarCrossPixel_ErrorEveryHDF1=CrossPixel_ErrorEveryHDF_Var1(2:interva3:end);
Many_VarCrossPixel_ErrorEveryHDF2=CrossPixel_ErrorEveryHDF_Var2(2:interva3:end);
figure;myBar3 = barwitherr([Many_VarCrossPixel_ErrorEveryHDF1, Many_VarCrossPixel_ErrorEveryHDF2], [CrossPixel_ErrorEveryHDF_Sampling1, CrossPixel_ErrorEveryHDF_Sampling2]);
set(myBar3(1), 'facecolor', 'r');set(myBar3(2), 'facecolor', 'b');
legend('修正前','修正后','Location','northwest');%,'Location','southeast'
xlabel('从所有HDF文件抽取样本','FontSize',13);
ylabel('跨轨误差 (像素)','FontSize',13);
saveas(gcf, [picName, '-跨轨直方图（每个HDF）.tif']);

% Cross-track-every day
interva4 = 2;%interval
CrossPixel_ErrorEveryday_Sampling1 = CrossPixel_ErrorEveryday1(1:interva4:end);
CrossPixel_ErrorEveryday_Sampling2 = CrossPixel_ErrorEveryday2(1:interva4:end);

Many_VarCrossPixel_ErrorEveryday1=CrossPixel_ErrorEveryday_Var1(1:interva4:end);
Many_VarCrossPixel_ErrorEveryday2=CrossPixel_ErrorEveryday_Var2(1:interva4:end);
figure;myBar4 = barwitherr([Many_VarCrossPixel_ErrorEveryday1, Many_VarCrossPixel_ErrorEveryday2], [CrossPixel_ErrorEveryday_Sampling1, CrossPixel_ErrorEveryday_Sampling2]);
set(myBar4(1), 'facecolor', 'r');set(myBar4(2), 'facecolor', 'b');
legend('修正前','修正后');%,'Location','southeast'
xlabel('从整天均值抽取样本','FontSize',13);
ylabel('跨轨误差 (像素)','FontSize',13);
saveas(gcf, [picName, '-跨轨直方图（每整天）.tif']);

%%
% Cross-track histogram-no fixed difference
interva_NoAdd3 = 2;%interval
CrossPixel_ErrorEveryHDF_Sampling_NoAdd1 = CrossPixel_ErrorEveryHDF_NoAdd1(2:interva_NoAdd3:end);
CrossPixel_ErrorEveryHDF_Sampling_NoAdd2 = CrossPixel_ErrorEveryHDF_NoAdd2(2:interva_NoAdd3:end);
Many_VarCrossPixel_ErrorEveryHDF_NoAdd1=CrossPixel_ErrorEveryHDF_VarNoAdd1(2:interva_NoAdd3:end);
Many_VarCrossPixel_ErrorEveryHDF_NoAdd2=CrossPixel_ErrorEveryHDF_VarNoAdd2(2:interva_NoAdd3:end);
figure;myBar_NoAdd3 = barwitherr([Many_VarCrossPixel_ErrorEveryHDF_NoAdd1, Many_VarCrossPixel_ErrorEveryHDF_NoAdd2], [CrossPixel_ErrorEveryHDF_Sampling_NoAdd1, CrossPixel_ErrorEveryHDF_Sampling_NoAdd2]);
set(myBar_NoAdd3(1), 'facecolor', 'r');set(myBar_NoAdd3(2), 'facecolor', 'b');
legend('修正前','修正后','Location','northwest');%,'Location','southeast'
xlabel('从所有HDF文件抽取样本','FontSize',13);
ylabel('跨轨误差 (像素)','FontSize',13);
saveas(gcf, [picName, '-跨轨直方图（每个HDF-无固定差）.tif']);

% Cross-track daily-no fixed difference
interva_NoAdd4 = 2;%interval
CrossPixel_ErrorEveryday_Sampling_NoAdd1 = CrossPixel_ErrorEveryday_NoAdd1(1:interva_NoAdd4:end);
CrossPixel_ErrorEveryday_Sampling_NoAdd2 = CrossPixel_ErrorEveryday_NoAdd2(1:interva_NoAdd4:end);
Many_VarCrossPixel_ErrorEveryday_NoAdd1=CrossPixel_ErrorEveryday_VarNoAdd1(1:interva_NoAdd4:end);
Many_VarCrossPixel_ErrorEveryday_NoAdd2=CrossPixel_ErrorEveryday_VarNoAdd2(1:interva_NoAdd4:end);
figure;myBar_NoAdd4 = barwitherr([Many_VarCrossPixel_ErrorEveryday_NoAdd1, Many_VarCrossPixel_ErrorEveryday_NoAdd2], [CrossPixel_ErrorEveryday_Sampling_NoAdd1, CrossPixel_ErrorEveryday_Sampling_NoAdd2]);
set(myBar_NoAdd4(1), 'facecolor', 'r');set(myBar_NoAdd4(2), 'facecolor', 'b');
legend('修正前','修正后');%,'Location','southeast'
xlabel('从整天均值抽取样本','FontSize',13);
ylabel('跨轨误差 (像素)','FontSize',13);
saveas(gcf, [picName, '-跨轨直方图（每整天-无固定差）.tif']);


CrossPixel_Error1_guiyi=(CrossPixel_Error1);
CrossPixel_Error2_guiyi=(CrossPixel_Error2);
AlongPixel_Error1_guiyi=(AlongPixel_Error1);
AlongPixel_Error2_guiyi=(AlongPixel_Error2);
CrossPixel_Error1_guiyi=j_Cross_1;
CrossPixel_Error2_guiyi=j_Cross_2;
AlongPixel_Error1_guiyi=j_Along_1;
AlongPixel_Error2_guiyi=j_Along_2;


cro_abs1=mean(Cross_mean_every_1);
disp(sprintf('调整前Cross方向像素偏移的绝对值为%f',cro_abs1));
cro_abs2=mean(Cross_mean_every_2);
disp(sprintf('调整后Cross方向像素偏移的绝对值为%f',cro_abs2));
alo_abs1=mean(Along_mean_every_1);
disp(sprintf('调整前Along方向像素偏移的绝对值为%f',alo_abs1));
alo_abs2=mean(Along_mean_every_2);
disp(sprintf('调整后Along方向像素偏移的绝对值为%f',alo_abs2));


cor_std1=mean(Cross_std_every_1);
disp(sprintf('调整前Cross方向像素偏移的标准差为%f',cor_std1));

cor_std2=mean(Cross_std_every_2);
disp(sprintf('调整后Cross方向像素偏移的标准差为%f',cor_std2));

alo_std1=mean(Along_std_every_1);
disp(sprintf('调整前Along方向像素偏移的标准差为%f',alo_std1));

alo_std2=mean(Along_std_every_2);
disp(sprintf('调整后Along方向像素偏移的标准差为%f',alo_std2));

figure;scatter(CrossPixel_Error1_guiyi,AlongPixel_Error1_guiyi,'b.');hold on;
scatter(mean(CrossPixel_Error1_guiyi),mean(AlongPixel_Error1_guiyi),'r*','LineWidth',1.5);
xlabel('Cross-Track error (Pixel)');ylabel('Along-Track error (Pixel)');

axis([-2 2 -2 2]);
set(gca,'XTick',[-2:0.4:2])
set(gca,'YTick',[-2:0.4:2])
grid on;
CrossPixel_Error2_guiyi_change=CrossPixel_Error2_guiyi+0.3*(mean(CrossPixel_Error2_guiyi)-CrossPixel_Error2_guiyi);
figure;scatter(CrossPixel_Error2_guiyi_change,AlongPixel_Error2_guiyi,'b.');hold on;
scatter(mean(CrossPixel_Error2_guiyi),mean(AlongPixel_Error2_guiyi),'r*','LineWidth',1.5);
xlabel('Cross-Track error (Pixel)');ylabel('Along-Track error (Pixel)');
%title('20140520_0142');
%axis([-0.085 0.085 -0.085 0.085]);
axis([-2 2 -2 2]);
set(gca,'XTick',[-2:0.4:2])
set(gca,'YTick',[-2:0.4:2])
grid on;


%%
save([theYear, '年', num2str(minMoon), '月至', num2str(maxMoon), '月跨、沿轨误差（每整天）.mat'],...
    'CrossPixel_ErrorEveryday1','CrossPixel_ErrorEveryday2','AlongPixel_ErrorEveryday1','AlongPixel_ErrorEveryday2');
save([theYear, '年', num2str(minMoon), '月至', num2str(maxMoon), '月跨、沿轨误差（每整天-方差）.mat'],...
    'CrossPixel_ErrorEveryday_Var1','CrossPixel_ErrorEveryday_Var2','AlongPixel_ErrorEveryday_Var1','AlongPixel_ErrorEveryday_Var2');
save([theYear, '年', num2str(minMoon), '月至', num2str(maxMoon), '月跨、沿轨误差（每个HDF）.mat'],...
    'CrossPixel_ErrorEveryHDF1','CrossPixel_ErrorEveryHDF2','AlongPixel_ErrorEveryHDF1','AlongPixel_ErrorEveryHDF2');
save([theYear, '年', num2str(minMoon), '月至', num2str(maxMoon), '月跨、沿轨误差（每个HDF-方差）.mat'],...
    'CrossPixel_ErrorEveryHDF_Var1','CrossPixel_ErrorEveryHDF_Var2','AlongPixel_ErrorEveryHDF_Var1','AlongPixel_ErrorEveryHDF_Var2');
save([theYear, '年', num2str(minMoon), '月至', num2str(maxMoon), '月跨、沿轨误差（每个HDF-不绝对值）.mat'],...
    'CrossPixel_ErrorEveryHDF_NoAbs1','CrossPixel_ErrorEveryHDF_NoAbs2','AlongPixel_ErrorEveryHDF_NoAbs1','AlongPixel_ErrorEveryHDF_NoAbs2');


save([theYear, '年', num2str(minMoon), '月至', num2str(maxMoon), '月跨、沿轨误差（每整天-无固定差）.mat'],...
    'CrossPixel_ErrorEveryday_NoAdd1','CrossPixel_ErrorEveryday_NoAdd2','AlongPixel_ErrorEveryday_NoAdd1','AlongPixel_ErrorEveryday_NoAdd2');
save([theYear, '年', num2str(minMoon), '月至', num2str(maxMoon), '月跨、沿轨误差（每整天-方差-无固定差）.mat'],...
    'CrossPixel_ErrorEveryday_VarNoAdd1','CrossPixel_ErrorEveryday_VarNoAdd2','AlongPixel_ErrorEveryday_VarNoAdd1','AlongPixel_ErrorEveryday_VarNoAdd2');
save([theYear, '年', num2str(minMoon), '月至', num2str(maxMoon), '月跨、沿轨误差（每个HDF-无固定差）.mat'],...
    'CrossPixel_ErrorEveryHDF_NoAdd1','CrossPixel_ErrorEveryHDF_NoAdd2','AlongPixel_ErrorEveryHDF_NoAdd1','AlongPixel_ErrorEveryHDF_NoAdd2');
save([theYear, '年', num2str(minMoon), '月至', num2str(maxMoon), '月跨、沿轨误差（每个HDF-方差-无固定差）.mat'],...
    'CrossPixel_ErrorEveryHDF_VarNoAdd1','CrossPixel_ErrorEveryHDF_VarNoAdd2','AlongPixel_ErrorEveryHDF_VarNoAdd1','AlongPixel_ErrorEveryHDF_VarNoAdd2');
save([theYear, '年', num2str(minMoon), '月至', num2str(maxMoon), '月跨、沿轨误差（每个HDF-不绝对值且无固定差）.mat'],...
    'CrossPixel_ErrorEveryHDF_NoAbsNoAdd1','CrossPixel_ErrorEveryHDF_NoAbsNoAdd2','AlongPixel_ErrorEveryHDF_NoAbsNoAdd1','AlongPixel_ErrorEveryHDF_NoAbsNoAdd2');