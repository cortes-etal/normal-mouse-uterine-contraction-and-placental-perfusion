% Devin Cortes
% Bioengineering, Department of Developmental Biology
% DCE analysis. this code is to calculate longitudinal statistics for the
% mice in the PCA4 study. These mice were image at E14.5 and E17.5 using
% subQ injection of Gd+ bolus for quantitative perfusion analysis


clear all;close all;clc; % start pretty


%%

files = dir('**/*_v2.mat');
testSig = zeros(50,1);
%%
segs = dir('*.nii.gz');
segNames = extractfield(segs,'name');

outDir = ['longitduinal_' ,date] ;
mkdir(outDir)

arrf=[];
for ii = 1:numel(files)
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
    %     mkdir(fullfile('PerfusionMaps',animal));
    load(fname);
    
    
    sMask = permute(sDat > 0, [2 1 3]);
    
    slopeMap=permute(steepest_slope/AIFmax*100, [2 1 3]);
    slopeMap2 = slopeMap.*sMask;
    
    etohFlag = contains(animal,'EtOH');
    flag17 = contains(animal,'17.5');
    
    if etohFlag && flag17
        %             etoh17_all = [etoh17_all; slopeMap2dat];
        gFlag = 2
    elseif etohFlag && ~flag17
        %             etoh14_all = [etoh14_all; slopeMap2dat];
        gFlag = 4
    elseif ~etohFlag && flag17
        %             cont17_all = [cont17_all; slopeMap2dat];
        gFlag = 1
    elseif ~etohFlag && ~flag17
        %             cont14_all = [cont14_all; slopeMap2dat];
        gFlag = 3
    end
    %
    %     figure;
    %     imshow3D(slopeMap2,[0 prctile(slopeMap2(:), 95)]);colormap(jet(256))
    %     %%
    
    
    labels = unique(sDat);
    labels= labels(2:end); % getting rid of background label
    [wholePerf, highPerf, lowPerf] = deal(zeros(size(labels)));
    allMasks = cell([1,numel(labels)]);
    
    animal
    fname = fullfile(outDir,saveDir,[animal,date]);
    [rgb,slopeLow,slopeHigh] = combineRGB_nD(enhancement,slopeMap2,[fname],animal,animalE);
    
    mean_low(ii) = mean(slopeLow(slopeLow > 0));
    mean_high(ii) = mean(slopeHigh(slopeHigh > 0));
    mean_whole(ii) = mean(slopeMap2(slopeMap2 > 0));
    allIds(ii) = str2double(animalID);
    allages(ii) = str2double(animalE);
    allflags(ii) = gFlag;
    
    
    
    
    % code for placental label splitting, not sure if need to do this for longitudinal analysis
        
        for l = 1:numel(labels)
            %         thisID = animalIDs{ii};
            placentaMasks{l} = sDat == l;
            pMask = sDat == l;
    
            pWhole = permute(pMask, [2 1 3]) .* slopeMap2;
            pHigh = permute(pMask, [2 1 3]) .* slopeHigh;
            pLow = permute(pMask, [2 1 3]) .* slopeLow;
    
            %         pdfDat = histFit_v3(slopeMap2,animal,fullfile(outDir,saveDir,[animal,date,'_',num2str(ii),'_.png']),etohFlag);
    
    
            %%
    
    
    
            highPerfVol = sum(sum(sum(pHigh >0)));
    
    
            lowPerfVol = sum(sum(sum(pLow >0)));
    
            meanHigh = mean(pHigh(pHigh > 0));
            meanLow = mean(pLow(pLow > 0));
            stdHigh =std(pLow(pLow > 0));
            stdLow =std(pLow(pLow > 0));
            meanWhole = mean(pWhole(pWhole > 0));
            stdWhole =std(pWhole(pWhole > 0));
    
    
    
            arr=[l;meanHigh; stdHigh; meanLow; stdLow; meanWhole; stdWhole ;highPerfVol; lowPerfVol; gFlag; str2double(animalID)];
            arrf=[arrf arr];
        end
    
    
end

%%

% labels = ['mean_high'; 'std_high','meanLow','meanWhole','stdWhole','highPerfVol','lowPerfVole','groupFlag','litter']
% array_write = [mean_low; mean_high; mean_whole; allIds; allages; allflags];
xlswrite(fullfile(outDir,['longitudinal_data__whole',date,'_.xlsx']),arrf);




%% calculating actual longitudinal differences


unique(allIds)


