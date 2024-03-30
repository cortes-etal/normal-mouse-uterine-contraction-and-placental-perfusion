

% purpose: create kinetic movies for each animal
%load each scan as a separate image, create movies across the different
%slices to create kinetic movie maps

clear ;close all;clc; % start pretty


%%

files = dir('**/*_v2.mat');
testSig = zeros(50,1);
%%
segs = dir('*.nii.gz');
segNames = extractfield(segs,'name');

outDir = ['opticflow_movies_high_res' ,date] ;
mkdir(outDir)


%%

for ii = 12:15%numel(files)
    ii
    %getting the animal IDs from directory to sort of the segmentaton files
    fprts = regexp(files(ii).folder,'\','split');
    animalID = fprts{4};
    animalE = fprts{5};
    saveDir = [fprts{4} fprts{5}];
    mkdir(fullfile(outDir,saveDir));
    animalIdx = regexp(animalID, '[\d\d\d]');
    animalID = animalID(animalIdx);
    %%
    animalEdx = regexp(animalE,'14'); %str location of animal emrbyonic agge
    
    while isempty(animalEdx) % checking if age is 17 or 14
        animalEdx = regexp(animalE,'17');
        disp('17 not 14');
        
        if isempty(animalEdx)
            animalEdx = regexp(animalE,'15');
            disp('actually 15');
        end
        
        
    end
    
    animalE = animalE(animalEdx:animalEdx+1); % the dataset animal embryonic age
    
    segID = contains(segNames, animalE) & contains(segNames,animalID); % logical vecotr containing idx for seg file
    sFile = fullfile(segs(segID).folder,segs(segID).name); % fname of seg file
    sDat = niftiread(sFile); % load the segmentation file % loading seg fle.\
    
    pMask = sDat > 0;
    
    fname = fullfile(files(ii).folder,files(ii).name);
    fprts2 =regexp(fname,'\','split');
    animal=[fprts2{4} fprts2{5}];
    
    mp4name = [fprts{4}, animalE];
    %     mkdir(fullfile('PerfusionMaps',animal));
    load(fname);
    %%
    opticFlow = opticalFlowHS;
%     opticFlow.NumFrames = 3;
%     opticFlow.NoiseThreshold = 11;
%     opticFlow.GradientFilterSigma = 20;
%     opticFlow.ImageFilterSigma = 25;
     opticFlow.Smoothness=10^(13.23);
% opticFlow.NumIterations = 100;
% opticFlow.FilterSize = 2;
% opticFlow.NumPyramidLevels = 8;
    h = figure;
    movegui(h);
    hViewPanel = uipanel(h,'Position',[0 0 1 1],'Title','Plot of Optical Flow Vectors');
    hPlot = axes(hViewPanel);
    %     label=fullfile(dirname,fileDir,'flowMap.mp4');
    
    slicetracker=1;
    while slicetracker <= size(ims,3)
        
        label =fullfile(outDir,saveDir, [mp4name,'_',num2str(slicetracker),'_.mp4']);
        vidObj = VideoWriter(label,'MPEG-4'); %initializing the vdeio object
        vidObj.FrameRate = 5;
        open(vidObj);
        
        for jj = 1:size(ims,4)
            
            
            flow = estimateFlow(opticFlow,ims(:,:,slicetracker,jj));
            
            imshow(ims(:,:,slicetracker,jj),[]); %plotting image
            hold on
            plot(flow,'DecimationFactor',[5 5],'ScaleFactor',300,'Parent',hPlot); % plot flow vecotrs as quiver arrows
            title(['frame num = ',num2str(jj)])
            q = findobj(gca,'type','Quiver'); %get quiver objects
            % Change color to red
            q.Color = 'y';
            drawnow;
            hold off
            if ii > 1
                currFrame = getframe(gcf);
                img_dat = currFrame.cdata;
                img_dat = imresize(img_dat,4,'Method','bicubic');
                currFrame.cdata=img_dat;
                writeVideo(vidObj,currFrame); % write current frame
            end
            Vxs(:,:,jj) = medfilt2(flow.Vx); %save X-direction of flow
            Vys(:,:,jj) = medfilt2(flow.Vy); %save y-direction of flow
            mag(:,:,jj) = sqrt(Vxs(:,:,jj).^2 + Vys(:,:,jj).^2);
            
        end
        
        close(vidObj);
        save(fullfile(outDir,saveDir,['velocityinfo_slice_',num2str(slicetracker),'_.mat']),'Vxs','Vys','-v7.3');
        
        slicetracker=slicetracker+1;
    end
    
end

%%

%end
