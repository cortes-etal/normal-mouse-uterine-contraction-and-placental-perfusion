clear all;
close all hidden;

mainpath='C:\Users\ddblabwu\Documents\MATLAB\MATLAB\code\code\Anthony\Anthony'; %PUT PATH TO THIS FOLDER HERE
addpath(fullfile(mainpath,'bruker_io'));

%% Find IntraGate datasets

numAnimals=4
%for control = 1:numAnimals
%     clearvars -except numAnimals
%     close all
folders = uigetfile_n_dir;
prts=regexp(folders{1},'\','split');
mainDir = strcat(prts{4},prts{5});

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

% deltat = median(diff(ACQ_abs_time))/60; %in minutes
deltat = diff(ACQ_abs_time);
[~,argmax]=max(deltat);
deltat(argmax)=[]; %remove single long delay
deltat = mean(deltat)/60; %in minutes

%% segmenting all placentas
enhancement=max(ims,[],4)-mean(ims(:,:,:,1:4),4); %height of peak
enhancement2 = enhancement;
[~,argmax]=max(ims,[],4); %time of peak
endFlag = 1;
placentaCounter=1;
sliceCounter=1;
h=figure(1);

while endFlag
    implay([enhancement/max(enhancement(:)); argmax/size(ims,4)])
    control = input('Enter "Y" to being segmentation\n','s');
    pMask=zeros(size(enhancement));
    while strcmp(control,'y')
        if sliceCounter ==1 
        fprintf('Beginning segmentation of placenta #%d\n',placentaCounter);
        firstSlice = input('Select First Slice of Placenta\n');
        else 
            firstSlice = firstSlice+1;
        end
        figure(1);
        imagesc(enhancement2(:,:,firstSlice)),colormap(gray(256));
        placenta=drawpolygon;
        placentaMask = createMask(placenta);
        pMask(:,:,firstSlice) = placentaMask;
%          if sliceCounter > 1
%             placentaMask = cat(3,pMask,placentaMask);
%         end
        control = input('Segment another Slice? y/n\n','s');
        if strcmp(control,'y')
            sliceCounter = sliceCounter+1;
        end
    end
%     close(h);
    fMask = pMask .* placentaCounter;
    allMasks{placentaCounter} = fMask;
    enhancement2(pMask == 1) = max(enhancement2(:));
    sliceCounter=1;
    endFlagInput = input('Segment another Placenta?\n','s');
    if ~strcmp(endFlagInput,'y')
        endFlag=0;
    end
end

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
AIF=drawpolygon;
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
placenta=drawpolygon;
placentaMask = createMask(placenta);
close(h)

pif=squeeze(sum(sum(ims(:,:,slicePlacenta,:).*placentaMask,1),2)/sum(placentaMask(:)));
PIFmax=max(pif)-mean(pif(1:4));


%% Display
% implay(enhancement/max(enhancement(:)))
%     implay(steepest_slope/AIFmax*100) %mL/min*100mL). Display max is 100 mL/(min*100mL)
%     implay(steepest_time/size(ims,4))
%     figure,plot(AIF-mean(AIF(1:4)))
%
% 
% mkdir(fullfile('SlopeMovies_v2',mainDir));
% mkdir(fullfile('SlopeImages_v2',mainDir));
% figure;
% %  for ii = 1:numel(slopeCalcs)
%     steepest_slope = slopeCalcs(ii).steepest;
%    AIFmax = aifMax(ii).aifMax;
% folderPrts = regexp(folders{ii},'\','split');
% scanNum = folderPrts{end};
% implay(enhancement/max(enhancement(:)))
combinedMask = zeros(size(enhancement));
for ii = 1:numel(allMasks)
    combinedMask = imadd(combinedMask,allMasks{ii});
end
combinedMask(combinedMask >=1) = 1;
combinedMask = permute(combinedMask,[2 1 3]);
slopeMap=permute(steepest_slope/AIFmax*100, [2 1 3]);
slopeMap2 = slopeMap.*combinedMask;
%%
% imageToMP4(slopeMap,fullfile('SlopeMovies_v2',mainDir,strcat('slopeMap','.mp4')));ylabel('Perfusin [100 mL / (min*100mL)]');
for jj=1:size(slopeMap2,3)
    imshow(slopeMap2(:,:,jj),[0 prctile(slopeMap2(:),95)]);colorbar;colormap(jet(256));
    % saveas(gcf,fullfile('SlopeImages',mainDir,['slopeMap_slice_',num2str(jj),'.fig']),'fig');
    % saveas(gcf,fullfile('SlopeImages',mainDir,['slopeMap_slice_',num2str(jj),'.png']),'png');
    
end
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
