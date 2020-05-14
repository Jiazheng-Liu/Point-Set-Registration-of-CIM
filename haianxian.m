%Match the red point set (detected coastline) and blue point set (Coastline truth point) obtained previously
%save the matching data in the target domain
%%
clc;clear;
warning off all;

load redPoint_0403_1836.mat;
load bulePoint_FY3C_0403_1836.mat;

p=[Inflected_latitudes,Inflected_longitudes];
q=[TruthLat',TruthLon'];

normalize = 1;
conf.lambda1 = 3;
conf.lambda2 = 0.05;
conf.lambda3 = 0.05;
conf.r = 0.08;
conf.a = 10;
conf.M = floor(size(Inflected_latitudes,1)/2);
%conf.beta = 0.1

show_result = 1;

X = p;    %Raw data
Y = q;
tic;

conf = LLT_init(conf);

[idx, dist] = knnsearch(Y,X,'dist','euclidean','k',1);
Y1 = Y(idx,:);  %Initialize Y1

ind = 1:size(p);

normal.xm=0; normal.ym=0;
normal.xscale=1; normal.yscale=1;

if normalize
    [nX, nY, normal]=norm_ind(X,Y1,ind);
end

if ~exist('conf'), conf = []; end

VecFld=MR(nX, nY, conf,ind);

for iii = 1:30
    if normalize,VecFld.TX=(VecFld.TX)*normal.yscale+repmat(normal.ym,size(VecFld.TX,1),1);end
    Xk = VecFld.TX;
    
   [idx, dist] = knnsearch(Y,Xk,'dist','euclidean','k',1); 
   Y1 = Y(idx,:); 

    if normalize
        [nX, nY, normal]=norm_ind(Xk,Y1,ind);
    end
    VecFld=MR(nX, nY, conf,ind);
end

if normalize,VecFld.TX=(VecFld.TX)*normal.yscale+repmat(normal.ym,size(VecFld.TX,1),1);end
Xk = VecFld.TX;
toc

  [idx, dist] = knnsearch(Y,Xk,'dist','euclidean','k',1); 
   Y1 = Y(idx,:); 
plot(X(:,1),X(:,2),'r.');  hold on; plot(Y1(:,1),Y1(:,2),'b.'); 

data=hdf5info('FY3C_MWRIA_GBAL_L1_20180403_1836_010KM_MS1.hdf');
Inflected_latitudes=hdf5read(data.GroupHierarchy.Groups.Datasets(20));%20 12
Inflected_longitudes=hdf5read(data.GroupHierarchy.Groups.Datasets(21));%21 13
figure;worldmap world;geoshow(Inflected_latitudes,Inflected_longitudes,'DisplayType','point');
load gshhs_land_f.mat;
hold on;geoshow(TruthLatitude,TruthLongitude,'Color','blue');
 
figure; 
[idx2, dist] = knnsearch(Y,X,'dist','euclidean','k',1); 

save('pipei_0403_1836','Y1');



