% DEvin Cortes
%this code will load previously created quantitative perfusion maps,
%present the user with the DCE image to segment some placentas, and then
%analyze the different placental compartments: high v low perfusion

clear all;close all;clc; % start pretty


%%
files = dir('**/*_v2.mat');
testSig = zeros(50,1);

segs = dir('*.nii.gz');
segNames = extractfield(segs,'name');

outDir = ['compartment_ROC' ,date] ;
mkdir(outDir);

%%
figure;
[etoh14_all cont14_all etoh17_all cont17_all] = deal([]);
count=1;
for ii = 5%:numel(files)
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
    
    fname = fullfile(files(ii).folder,files(ii).name);
    fprts2 =regexp(fname,'\','split');
    animalIDs{ii}=[fprts2{4} fprts2{5}];
    
    animalE = animalE(animalEdx:animalEdx+1); % the dataset animal embryonic age
    
    segID = contains(segNames, animalE) & contains(segNames,animalID); % logical vecotr containing idx for seg file
    sFile = fullfile(segs(segID).folder,segs(segID).name); % fname of seg file
    sDat = niftiread(sFile); % load the segmentation file % loading seg fle.
    
    fname = fullfile(files(ii).folder,files(ii).name);
    fprts2 =regexp(fname,'\','split');
    animal=[fprts2{4} fprts2{5}];
    %     mkdir(fullfile('PerfusionMaps',animal));
    load(fname);
    
    sMask = permute(sDat > 0, [2 1 3]);
         slopeMap=permute(steepest_slope/AIFmax*100, [2 1 3]);
    slopeMap2 = slopeMap.*sMask;
%     
     [rgb,slopeLow,slopeHigh] = combineRGB_nD(enhancement,slopeMap2,['histo_',num2str(ii),'_.png'],animal);

    
     lowVec = slopeLow(slopeLow >0);
     highVec = slopeHigh(slopeHigh > 0);
     wholeVec = slopeMap2(slopeMap2 > 0);
      ranges = log10([min(lowVec(:)), max(highVec(:))]);
      slopeMap2dat = log10(slopeMap2(slopeMap2 >0));
     cc=1;
     for thresh = ranges(1):1:ranges(2)
         lowROC(cc) = sum(slopeMap2dat < thresh) ./ numel(slopeMap2dat);
         highROC(cc) = sum(slopeMap2dat >= thresh) ./ numel(slopeMap2dat);
         cc=cc+1;
     end
     figure;
     plot(lowROC,highROC);
     pause;
    
   
    
    etohFlag = contains(animal,'EtOH');
    flag17 = contains(animal,'17.5');
    
    if etohFlag && flag17
        etoh17_all = [etoh17_all; slopeMap2dat];
        gFlag = 'etoh17';
    elseif etohFlag && ~flag17
        etoh14_all = [etoh14_all; slopeMap2dat];
        gFlag = 'etoh14';
    elseif ~etohFlag && flag17
        cont17_all = [cont17_all; slopeMap2dat];
        gFlag = 'cont17';
    elseif ~etohFlag && ~flag17
        cont14_all = [cont14_all; slopeMap2dat];
        gFlag = 'cont14'
    end
    
    %% insert placental labeling splitting part
    
    labels = unique(sDat);
    labels= labels(2:end); % getting rid of background label
    [wholePerf, highPerf, lowPerf] = deal(zeros(size(labels)));
    allMasks = cell([1,numel(labels)]);
    
    
    for l = 1:numel(labels)
        thisID = animalIDs{ii};
        placentaMasks{l} = sDat == l;
        pMask = sDat == l;
        
        pWhole = permute(pMask, [2 1 3]) .* slopeMap2;
        pHigh = permute(pMask, [2 1 3]) .* slopeHigh;
        pLow = permute(pMask, [2 1 3]) .* slopeLow;

    end
    
end

%%
variables = {'peakLow','peakHigh','peakWhole','ttpLow','ttpHigh','ttpWhole','etoh','17_5','litter','group'};
cell2_Write = cat(1,variables,finalCell);
writecell(cell2_Write,fullfile(outDir,['allplacentalData_',date,'_.csv']))


