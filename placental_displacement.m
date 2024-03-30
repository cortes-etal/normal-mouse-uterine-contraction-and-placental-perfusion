

% purpose: create kinetic movies for each animal
%load each scan as a separate image, create movies across the different
%slices to create kinetic movie maps

clear all;close all;clc; % start pretty


%%

files = dir('opticflow_movies17-May-2022\**\*.mat');
testSig = zeros(50,1);
%%
segs = dir('*.nii.gz');
segNames = extractfield(segs,'name');

outDir = ['opticflow_movies_and_displacement' ,date] ;
mkdir(outDir)


%%

for ii = 1:numel(files)
    fname_p = fullfile(files(ii).folder,files(ii).name);
    load(fname_p);
    [X,Y] = meshgrid(0:1:255, 0:1:159);
    
    for jj = 1:size(Vxs,3)
        
        disp_field(:,:,jj) = divergence(X,Y,Vxs(:,:,jj), Vys(:,:,jj));
        
        split_fname = split(files(ii).name,'.');
        vid_out = [split_fname{1}, '_displacement_.mp4'];
    end
    
    imageToMP4(disp_field,fullfile(files(ii).folder,vid_out),jet(256), [0 prctile(disp_field(:), 95)],1);
    
end

%%

%end
