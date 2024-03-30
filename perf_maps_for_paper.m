%%
clear all;clc;close all;


files = dir('**/*_v2.mat');
testSig = zeros(50,1);
%%
segs = dir('*.nii.gz');
segNames = extractfield(segs,'name');

outDir = ['kinetic_movies_nifti_exports' ,date] ;
mkdir(outDir)

%% 

e17 = files(7);
e14 = files(4) ;


%% doing e`4 first

load(fullfile(e14.folder,e14.name));

segDat = niftiread(fullfile(segs(4).folder,segs(4).name));

perfMap=(steepest_slope/AIFmax*100);%, [2 1 3]);

segDat = segDat >0;

regions_applied = perfMap(:,:,13) .* segDat(:,:,13);

unique(segDat)

 [splitmap,combmap,~,~] = combineRGB_nD(steepest_slope(:,:,13),regions_applied,'test','test','test');

 
 figure;
 subplot(1,2,1)
 imshow(splitmap);
 subplot(1,2,2);
 imshow(combmap);
 
 
 %% now doing e17
 
 load(fullfile(e17.folder,e17.name));

segDat = niftiread(fullfile(segs(7).folder,segs(7).name));

perfMap=(steepest_slope/AIFmax*100);%, [2 1 3]);

segDat = segDat >0; 
regions_applied = perfMap(:,:,14) .* segDat(:,:,14);

unique(segDat)

 [splitmap,combmap,~,~] = combineRGB_nD(steepest_slope(:,:,13),regions_applied,'test','test','test');

 figure;
 subplot(1,2,1)
 imshow(splitmap);
 subplot(1,2,2);
 imshow(combmap);
 