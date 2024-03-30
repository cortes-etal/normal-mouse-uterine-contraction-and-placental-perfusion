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

outDir = ['compartment_CurveDat' ,date] ;
mkdir(outDir);

%%
figure;
[etoh14_all cont14_all etoh17_all cont17_all] = deal([]);
count=1;
for ii = 4%1:numel(files)
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
    
    sMask = (sDat > 0);
    
    slopeMap=(steepest_slope/AIFmax*100);
    slopeMap2 = slopeMap.*sMask;
    
    [rgb,~,slopeLow,slopeHigh] = combineRGB_nD(enhancement,slopeMap2,['histo_',num2str(ii),'_.png'],animal,animalE);
    
    
    slopeMap2dat = slopeMap2(slopeMap2 >0);
    
    
    
    
    etohFlag = contains(animal,'EtOH');
    flag17 = contains(animal,'17.5');
    
    if etohFlag && flag17
        etoh17_all = [etoh17_all; slopeMap2dat];
        gFlag = 2
    elseif etohFlag && ~flag17
        etoh14_all = [etoh14_all; slopeMap2dat];
        gFlag = 4
    elseif ~etohFlag && flag17
        cont17_all = [cont17_all; slopeMap2dat];
        gFlag = 1
    elseif ~etohFlag && ~flag17
        cont14_all = [cont14_all; slopeMap2dat];
        gFlag = 3
    end
    
    %% insert placental labeling splitting part
    
    labels = unique(sDat);
    labels= labels(2:end); % getting rid of background label
    [wholePerf, highPerf, lowPerf] = deal(zeros(size(labels)));
    allMasks = cell([1,numel(labels)]);
    
    highPerfVol(ii,1) = sum(sum(highPerf >0));
    highPerfVol(ii,2) = gFlag;
    lowPerfVol(ii,1) = sum(sum(lowPerf >0));
    lowerPerfVol(ii,2) = gFlag;
    
    
    for l = 1:numel(labels)
        thisID = animalIDs{ii};
        placentaMasks{l} = sDat == l;
        pMask = sDat == l;
        
        pWhole = (pMask) .* slopeMap2;
        pHigh = (pMask) .* slopeHigh;
        pLow = (pMask) .* slopeLow;
        
        wholeVol = numel(pWhole(pWhole > 0));
         highVol = numel(pHigh(pHigh > 0));
          lowVol = numel(pLow(pLow > 0));
        
        
        lowCurve=squeeze(sum(sum(sum(ims.*(pLow > 0),1),2),3)/sum(AIFmask(:)));
        highCurve=squeeze(sum(sum(sum(ims.*(pHigh > 0),1),2),3)/sum(AIFmask(:)));
        wholeCurve=squeeze(sum(sum(sum(ims.*(pWhole > 0),1),2),3)/sum(AIFmask(:)));
        
        lowCurve = lowCurve - mean(lowCurve(1:4));
        highCurve = highCurve - mean(highCurve(1:4));
        wholeCurve = wholeCurve - mean(wholeCurve(1:4));
        
        figure(111);
        subplot(1,2,1)
        plot(ACQ_abs_time./3600,lowCurve./rms(lowCurve));
        hold on;
        plot(ACQ_abs_time./3600,highCurve./rms(highCurve));
        plot(ACQ_abs_time./3600,wholeCurve./rms(wholeCurve));
        hold off;
        xlabel('Minutes [m]');
        title('Intensity Curve');
        legend('Low Perfusion Chamber','High Perfusion Chamber','Whole Placenta','Location','southwest');
        drawnow;
        subplot(1,2,2)
        plot(ACQ_abs_time(1:end-1)./3600,diff(lowCurve./rms(lowCurve)));
        hold on;
        plot(ACQ_abs_time(1:end-1)./3600,diff(highCurve./rms(highCurve)));
        plot(ACQ_abs_time(1:end-1)./3600,diff(wholeCurve./rms(wholeCurve)));
        hold off;
        xlabel('Minutes [m]');
        title(animal);
        legend('Low Perfusion Chamber','High Perfusion Chamber','Whole Placenta','Location','southwest');
        drawnow;
        saveas(gcf,fullfile(outDir,saveDir,['intensityCurves_placenta_',num2str(l),'_.png']),'png')
        slopeCurveLow = diff(lowCurve ./rms(lowCurve));
        slopeCurveHigh = diff(highCurve ./ rms(highCurve));
        slopeCurveWhole = diff(wholeCurve ./ rms(wholeCurve));
        
        if ii == 13
            slopeCurveLow = slopeCurveLow(10:end);
            slopeCurveHigh = slopeCurveHigh(10:end);
            slopeCurveWhole = slopeCurveWhole(10:end);
            
        end
        [peak,ttp] = max(slopeCurveLow);
        [peakH, ttpH] = max(slopeCurveHigh);
        [peakW, ttpW] =max(slopeCurveWhole);
        
        lowVec = pLow(pLow >0);
        highVec = pHigh(pHigh > 0);
        wholeVec = pWhole(pWhole > 0);
    
        
        finalCell{count,1} = peak;
        finalCell{count,2} = peakH;
        finalCell{count,3} =peakW;
        finalCell{count,4}= ttp;
        finalCell{count,5} = ttpH;
        finalCell{count,6} = ttpW;
        finalCell{count,7} = lowVol;
        finalCell{count,8} = highVol;
        finalCell{count,9} = wholeVol;
        finalCell{count,10} = mean(wholeVec);
        finalCell{count,11} = mean(lowVec);
        finalCell{count,12} = mean(highVec);
        finalCell{count,13} = std(wholeVec);
        finalCell{count,14} = std(lowVec);
        finalCell{count,15} = std(highVec);
%         finalCell{count,16} = histcounts(log10(wholeVec),100);
%         finalCell{count,17} = histcounts(log10(lowVec),100);
%         finalCell{count,18} = histcounts(log10(highVec),100);
        finalCell{count,16} = trapz(highCurve);
        finalCell{count,17} = trapz(lowCurve);
        finalCell{count,18} = trapz(wholeCurve);
        finalCell{count,19} = max(diff(AIF));
        finalCell{count,20} = max(AIF) ./ mean(AIF);
        finalCell{count,21} = animal;
        finalCell{count,22} = gFlag;
        count=count+1;
        
        
    end
    
end

%%

variables = {'peakLow','peakHigh','peakWhole','ttpLow','ttpHigh','ttpWhole','volumeLow','volumeHigh','volume Whole','meanWhole','meanLow','meanHigh','stdWhole','stdLow','stdHigh','aoc_high','aoc_low','aoc_whole','AIFMAX','normAIF_max','litter','group'};
%
cell2_Write = cat(1,variables,finalCell);
%
writecell(cell2_Write,fullfile(outDir,['allplacentalData_newThresh_',date,'_.csv']))



