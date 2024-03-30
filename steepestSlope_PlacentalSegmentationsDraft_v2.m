clear all;
close all hidden;

mainpath='C:\Users\ddblabwu\Documents\MATLAB\MATLAB\code\code\Anthony\Anthony'; %PUT PATH TO THIS FOLDER HERE
addpath(fullfile(mainpath,'bruker_io'));

%% first loading R1 post maps and doing placental segmentations
% Devin Cortes
% Dr Yijen Wu
% purpose of this code is to use FA 75 for manual segmentation of placentas
% and further DCE analysis.

clear all;clc;close all; % start pretty

mainpath='C:\Users\ddblabwu\Documents\MATLAB\MATLAB\code\code\Anthony\Anthony'; %PUT PATH TO THIS FOLDER HERE
addpath(fullfile(mainpath,'bruker_io'));

%% %%
% 
% folders = uigetfile_n_dir;
% 
% k=1;
% for mm = 1:numel(folders)
%     mm
%     %  k=1;
%     r3=@(x)reshape(x,1,1,[]);
%     %         orderedFolders = orderedDirs(ii).cell;
%     %         for j=1:numel(orderedFolders)
%     foldername=folders{mm};
%     acqp{mm}=read_bruker_acqp(fullfile(foldername,'acqp'));
%     try
%         reco=read_bruker_acqp(fullfile(foldername,'pdata/1/reco'));
%     catch %try again in case Box just isn't responsive the first time
%         reco=read_bruker_acqp(fullfile(foldername,'pdata/1/reco'));
%     end
%     im=load_bruker_2dseq(reco.RECO_size(1),fullfile(foldername,'pdata/1/2dseq'),'int32');
%     im=reshape(im,reco.RECO_size(2),reco.RECO_size(1),[]);
%     ims(:,:,1:size(im,3),k)=bsxfun(@rdivide,im,r3(reco.RECO_map_slope)); %store the image, correcting for the storage scaling factor
%     alphas(1,k)=acqp{mm}.ACQ_flip_angle; %also store the flip angle
%     try
%         TRs(1,k)=acqp{mm}.ACQ_repetition_time*1e-3; %also store the repetition time
%     catch
%         acqR = acqp{mm}.ACQ_repetition_time(1)*1e-3;
%         TRs(1,k) = acqR;
%     end
%     k=k+1;
%     %  end
%     %         allIms(mm).im = ims;
% end
% 
% prts = regexp(folders{1},'\','split');
% mainDir=strcat(prts{4},prts{5});
% FA=75;
% imToSeg = squeeze(ims(:,:,:,find(alphas == FA)));
% 
% figure;
% imshow3D(imToSeg,[])
% 
% figure;
% enhancement=max(ims,[],4)-mean(ims(:,:,:,1:4),4); %height of peak
% [~,argmax]=max(ims,[],4); %time of peak
% implay([enhancement/max(enhancement(:)); argmax/size(ims,4)])
% slice=input('Select the Slice of the image that has most placentas\n');
% numPlacentas = input('How many placentas to segment?\n');
% 
% h=figure(1);
% [combinedMasks1,combinedLabels] = deal(zeros(size(squeeze(imToSeg(:,:,slice)))));
% for ii = 1:numPlacentas
%     figure(1);
%     imshow(imToSeg(:,:,slice),[]);
%     placenta =imfreehand;
%     placentaMask = createMask(placenta);
%     masks{ii} = placentaMask.*ii;
%     combinedMasks1 = imadd(combinedMasks1,double(placentaMask));
%     combinedLabels = imadd(combinedLabels, double(masks{ii}));
% end
% 
% figure;
% subplot(2,1,1)
% imshow(combinedMasks1,[]);
% subplot(2,1,2)
% imshow(label2rgb(combinedLabels));


%% Find IntraGate datasets for quantittive perfusion maps

numAnimals=4
%for control = 1:numAnimals
%     clearvars -except numAnimals
%     close all
folders = uigetfile_n_dir;
prts=regexp(folders{1},'\','split');
mainDir = strcat(prts{4},prts{5});
clear ims
clear im
clear acqp

k=1;
r3=@(x)reshape(x,1,1,[]);
for j = 1:numel(folders)
    j
    foldername=folders{j};
    
    % Load IntraGate datasets
    %  k=1;
    acqp{j}=read_bruker_acqp(fullfile(foldername,'acqp'));
    %for j=find(isIntraGate)
    try
        reco=read_bruker_acqp(fullfile(foldername,'pdata/1/reco'));
    catch %try again in case Box just isn't responsive the first time
        reco=read_bruker_acqp(fullfile(foldername,'pdata/1/reco'));
    end
    im=load_bruker_2dseq(reco.RECO_size(1),fullfile(foldername,'pdata/1/2dseq'),'int32');
    im=reshape(im,reco.RECO_size(2),reco.RECO_size(1),[]);
    try
        ims(:,:,1:size(im,3),k)=bsxfun(@rdivide,im,r3(reco.RECO_map_slope)); %store the image, correcting for the storage scaling factor
        ACQ_time{1,k}=acqp{j}.ACQ_time;
        ACQ_abs_time(1,k)=acqp{j}.ACQ_abs_time;
        k=k+1;
    catch
        break;
    end
    
end

%%
figure;
enhancement=max(ims,[],4)-mean(ims(:,:,:,1:4),4); %height of peak
[~,argmax]=max(ims,[],4); %time of peak
implay([enhancement/max(enhancement(:)); argmax/size(ims,4)])
slice=input('Select the Slice of the image that has most placentas\n');
numPlacentas = input('How many placentas to segment?\n');

h=figure(1);
[combinedMasks1,combinedLabels] = deal(zeros(size(squeeze(enhancement(:,:,slice)))));
for ii = 1:numPlacentas
    figure(1);
    imshow(enhancement(:,:,slice),[]);
    placenta =imfreehand;
    placentaMask = createMask(placenta);
    masks{ii} = placentaMask.*ii;
    combinedMasks1 = imadd(combinedMasks1,double(placentaMask));
    combinedLabels = imadd(combinedLabels, double(masks{ii}));
end

figure;
subplot(2,1,1)
imshow(combinedMasks1,[]);
subplot(2,1,2)
imshow(label2rgb(combinedLabels));




% deltat = median(diff(ACQ_abs_time))/60; %in minutes
deltat = diff(ACQ_abs_time);
[~,argmax]=max(deltat);
deltat(argmax)=[]; %remove single long delay
deltat = mean(deltat)/60; %in minutes




% 
% 
% pif=squeeze(sum(sum(ims(:,:,slicePlacenta,:).*placentaMask,1),2)/sum(placentaMask(:)));
% PIFmax=max(pif)-mean(pif(1:4));

%% Steepest slope

% closed-form solution for three-point linear fit is just (last - first)/2:
%syms x1 x2 x3;
%temp=[x1 x2 x3]*pinv([-1 0 1; 1 1 1]);
%slope=temp(1);

% remove signal drops from motion along time dimension ("morphological closing")
processed=max(cat(5,ims(:,:,:,1:end-2),ims(:,:,:,2:end-1),ims(:,:,:,3:end)),[],5);
processed=min(cat(5,processed(:,:,:,1:end-2),processed(:,:,:,2:end-1),processed(:,:,:,3:end)),[],5);
processed=cat(4,ims(:,:,:,1:2),processed,ims(:,:,:,end-1:end));

% spatial median filter (light smoothing, also helpful with motion)
processed=medfilt3(reshape(padarray(processed,[0 0 1 0]),size(ims,1),size(ims,2),[]),[3 3 3]);
processed=reshape(processed,size(ims,1),size(ims,2),[],size(ims,4));
processed=processed(:,:,2:end-1,:);

[steepest_slope,steepest_time]=max((processed(:,:,:,3:end-2)-processed(:,:,:,1:end-4))/2/deltat,[],4);  %ignore last two "unclosed" time points

%% get AIF from kidney blood supply
enhancement=max(ims,[],4)-mean(ims(:,:,:,1:4),4); %height of peak
[~,argmax]=max(ims,[],4); %time of peak
implay([enhancement/max(enhancement(:)); argmax/size(ims,4)])
slice=input('Which slice? ');

h=figure;
imagesc(enhancement(:,:,slice)),colormap(jet(256));
AIF=imfreehand;
AIFmask = createMask(AIF);
close(h)

AIF=squeeze(sum(sum(ims(:,:,slice,:).*AIFmask,1),2)/sum(AIFmask(:)));
AIFmax=max(AIF)-mean(AIF(1:4));

%% get plancental perfusion curve
enhancement=max(ims,[],4)-mean(ims(:,:,:,1:4),4); %height of peak
[~,argmax]=max(ims,[],4); %time of peak
implay([enhancement/max(enhancement(:)); argmax/size(ims,4)])
slicePlacenta=input('Which slice? (Pick a Placenta, peak slice) ');

h=figure;
imagesc(enhancement(:,:,slicePlacenta)),colormap(jet(256));
placenta=imfreehand;
placentaMask = createMask(placenta);
close(h)

pif=squeeze(sum(sum(ims(:,:,slicePlacenta,:).*placentaMask,1),2)/sum(placentaMask(:)));
PIFmax=max(pif)-mean(pif(1:4));


%% Display

combinedMask = permute(combinedMasks1,[2 1 3]);
slopeMap=permute(steepest_slope/AIFmax*100, [2 1 3]);
slopeMap2 = slopeMap(:,:,slice).*combinedMask;
%%
% imageToMP4(slopeMap,fullfile('SlopeMovies_v2',mainDir,strcat('slopeMap','.mp4')));ylabel('Perfusin [100 mL / (min*100mL)]');
figure;
imshow(slopeMap2,[0 prctile(slopeMap2(:),95)]);colorbar;colormap(jet(256));
greydouble =permute(double2rgb(squeeze(enhancement(:,:,slice)),gray(256)),[2 1 3]);
mapDouble = double2rgb(slopeMap2,colormap(jet(256)), [0 prctile(slopeMap2(:),95)]);

outR = greydouble(:,:,1);
outG = greydouble(:,:,2);
outB = greydouble(:,:,3);

mapR = mapDouble(:,:,1);
mapG = mapDouble(:,:,2);
mapB = mapDouble(:,:,3);

outR(slopeMap2> 0) = mapR(slopeMap2 >0);
outG(slopeMap2> 0) = mapG(slopeMap2 >0);
outB(slopeMap2> 0) = mapB(slopeMap2 >0);

finalIm = cat(3,outR,outG,outB);

figure;
imshow(finalIm)

    % aveas(gcf,fullfile('SlopeImages',mainDir,['slopeMap_slice_',num2str(jj),'.fig']),'fig');
    % saveas(gcf,fullfile('SlopeImages',mainDir,['slopeMap_slice_',num2str(jj),'.png']),'png');
    
%%
%mL/min*100mL). Display max is 100 mL/(min*100mL)
%%imshow3D(steepest_time/size(ims,4));colormap(autumn(256));
% pause;
figure,plot((AIF-mean(AIF(1:4)))./max(AIF-mean(AIF(1:4))));title('Perfusion Time Curves');
hold on;
plot((pif-mean(pif(1:4)))./max(pif-mean(pif(1:4))));
xlabel('Time [min]');
legend('Arterial Input Function','Placental Perfusion Time Curve','Location','southwest');

%     saveas(gcf,fullfile('SlopeImages_v2',mainDir,strcat('aif.fig')));
%     saveas(gcf,fullfile('SlopeImages_v2',mainDir,strcat('aif.png')));
%     %  end
%
%     %% Save
%     save(fullfile(foldername,'DCE_test_v2.mat'),'ims','steepest_slope','AIFmax','AIF','AIFmask','slice','enhancement','deltat','ACQ_abs_time','pif','placentaMask','slicePlacenta','PIFmax');
%end %control loop
