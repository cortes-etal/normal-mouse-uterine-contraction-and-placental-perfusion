

% purpose: create kinetic movies for each animal
%load each scan as a separate image, create movies across the different
%slices to create kinetic movie maps

%%
numAnimals = 6;
clear all;
close all;
%for control = 1:numAnimals
folders = uigetfile_n_dir;

mainpath='C:\Users\ddblabwu\Documents\MATLAB\MATLAB\code\code\Anthony\Anthony'; %PUT PATH TO THIS FOLDER HERE
addpath(fullfile(mainpath,'bruker_io'));

%     for ii = 1:numel(folders)
%         temp = circshift(folders, [0 -ii]);
%         orderedDirs(ii).cell=temp;
%     end
%
k=1;

for j = 1:numel(folders)
    j
    % k=1;
    r3=@(x)reshape(x,1,1,[]);
    foldername = folders{j};
    % for j=1:numel(orderedFolders)
    % foldername=orderedFolders{j};
    acqp{j}=read_bruker_acqp(fullfile(foldername,'acqp'));
    try
        reco=read_bruker_acqp(fullfile(foldername,'pdata/1/reco'));
    catch %try again in case Box just isn't responsive the first time
        reco=read_bruker_acqp(fullfile(foldername,'pdata/1/reco'));
    end
    im=load_bruker_2dseq(reco.RECO_size(1),fullfile(foldername,'pdata/1/2dseq'),'int32');
    im=reshape(im,reco.RECO_size(2),reco.RECO_size(1),[]);
    ims(:,:,1:size(im,3),k)=bsxfun(@rdivide,im,r3(reco.RECO_map_slope)); %store the image, correcting for the storage scaling factor
    alphas(1,k)=acqp{j}.ACQ_flip_angle; %also store the flip angle
    try
        TRs(1,k)=acqp{j}.ACQ_repetition_time*1e-3; %also store the repetition time
    catch
        acqR = acqp{j}.ACQ_repetition_time(1)*1e-3;
        TRs(1,k) = acqR;
    end
    k=k+1;
    % end
    
    
   % Nims = numel(allIms);
    
    prts = regexp(folders{1},'\','split');
    maindir= prts{5};
    
    
    
    % for ii =1:numel(allIms)
%     ii
%     I = allIms(ii).im;
  %  I = ims;
   % enhancement=(im - mean(im(:)))/(std(im(:)));
   % enhancement =(im-max(im(:)))/(max(im(:)) - min(im(:)));
        % enhancement = max(I,[],4)-mean(I(:,:,:,1:4),4);
    enhancement=permute(im, [2 1 3]);
    
    %     figure;
    %    imshow3D(enhancement)
    slicesOne(:,:,j) = enhancement(:,:,1);
    slicesTwo(:,:,j) = enhancement(:,:,2);
    slicesThree(:,:,j) = enhancement(:,:,3);
    slicesFour(:,:,j) = enhancement(:,:,4);
    slicesFive(:,:,j) = enhancement(:,:,5);
    
    
end

%%
mkdir(fullfile('kineticMovies_16Aug',maindir));

imageToMP4(slicesOne,fullfile('kineticMovies_16Aug',maindir,'sliceOne.mp4'));
imageToMP4(slicesTwo,fullfile('kineticMovies_16Aug',maindir,'slicesTwo.mp4'));
imageToMP4(slicesThree,fullfile('kineticMovies_16Aug',maindir,'slicesThree.mp4'));
imageToMP4(slicesFour, fullfile('kineticMovies_16Aug',maindir,'slicesFour.mp4'));
imageToMP4(slicesFive,fullfile('kineticMovies_16Aug',maindir, 'slicesFive.mp4'));
%end
