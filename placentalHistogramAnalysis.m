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

outDir = ['compartment3d_histo_' ,date] ;
mkdir(outDir);

%%
figure;
[etoh14_all cont14_all etoh17_all cont17_all] = deal([]);

%%
for ii = 1:numel(files)
    
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
    
    slopeMap2dat = slopeMap2(slopeMap2 >0);
    
    etohFlag = contains(animal,'EtOH');
    flag17 = contains(animal,'17.5');
    
    if etohFlag && flag17
        etoh17_all = [etoh17_all; slopeMap2dat];
        gFlag = 2;
    elseif etohFlag && ~flag17
        etoh14_all = [etoh14_all; slopeMap2dat];
        gFlag = 4;
    elseif ~etohFlag && flag17
        cont17_all = [cont17_all; slopeMap2dat];
        gFlag=1;
    elseif ~etohFlag && ~flag17
        cont14_all = [cont14_all; slopeMap2dat];
        gFlag = 3;
    end
    
    %% insert placental labeling splitting part
    
    labels = unique(sDat);
    labels= labels(2:end); % getting rid of background label
    [wholePerf, highPerf, lowPerf] = deal(zeros(size(labels)));
    allMasks = cell([1,numel(labels)]);
    
    
    %     for l = 1:numel(labels)
    %         thisID = animalIDs{ii};
    %         placentaMasks{l} = sDat == l;
    %         pMask = sDat == l;
    %
    %         pWhole = permute(pMask, [2 1 3]) .* slopeMap2;
    %         pHigh = permute(pMask, [2 1 3]) .* slopeHigh;
    %         pLow = permute(pMask, [2 1 3]) .* slopeLow;
    
    
    [pdfDat, s] = histFit_v3(slopeMap2,animal,fullfile(outDir,saveDir,[animal,date,'_',num2str(ii),'_.png']),etohFlag);
    %     end
    
    %%
    if etohFlag
        sumStats1 = pdfDat(1).objSingle;
        %     sumStats2 = pdfDat(2).objSingle;
        cellData{ii,1} = pdfDat(1).peak;
        cellData{ii,2} = nan;
        cellData{ii,3} = pdfDat(1).halfwidth;
        cellData{ii,4} = nan;
        cellData{ii,5} = pdfDat(1).median;
        cellData{ii,6} = nan;
        cellData{ii,7} = sumStats1.mu;
        cellData{ii,8} = nan;
        cellData{ii,9} = sumStats1.Sigma;
        cellData{ii,10} =nan;
        cellData{ii,11} = pdfDat(1).AOC;
        cellData{ii,12} = nan;
        cellData{ii,13} = s;

        cellData{ii,14} = gFlag;
        
    else
        sumStats1 = pdfDat(1).objSingle;
        sumStats2 = pdfDat(2).objSingle;
        cellData{ii,1} = pdfDat(1).peak;
        cellData{ii,2} = pdfDat(2).peak;
        cellData{ii,3} = pdfDat(1).halfwidth;
        cellData{ii,4} = pdfDat(2).halfwidth;
        cellData{ii,5} = pdfDat(1).median;
        cellData{ii,6} = pdfDat(2).median;
        cellData{ii,7} = sumStats1.mu;
        cellData{ii,8} = sumStats2.mu;
        cellData{ii,9} = sumStats1.Sigma;
        cellData{ii,10} = sumStats2.Sigma;
        cellData{ii,11} = pdfDat(1).AOC;
        cellData{ii,12} = pdfDat(2).AOC;
        cellData{ii,13} = s;
        cellData{ii,14} = gFlag;
        
        
    end
    
    
end
variables = {'peakLowDist','peakHighDist','widthLow','widthHigh','medianLow','medianHigh','meanLow','meanHigh','sigmaLow','sigmaHigh','aocLow','aocHigh','skew','group'};
cell2_Write = cat(1,variables,cellData);
writecell(cell2_Write,fullfile(outDir,['allhistoStats_newThresh',date,'_.csv']))

%%

for ii = 1:size(md,2)
    
    means(ii) = mean(nonzeros(md(:,ii)));
    mean_ssmd(ii) = mean(nonzeros(ssmd(:,ii)));
    mean_overlap(ii) = mean(nonzeros(overlaps(:,ii)));
    
    meanLowMu(ii) = mean(nonzeros(muLow(:,ii)));
    meanHighMu(ii) = mean(nonzeros(muHigh(:,ii)));
    
    meanLowSigma(ii) = mean(nonzeros(varLow(:,ii)));
    meanHighSigma(ii) = mean(nonzeros(varHigh(:,ii)));
    
    meanHighAoc(ii) = mean(nonzeros(aocHigh(:,ii)));
    meanLowAoc(ii) = mean(nonzeros(aocLow(:,ii)));
    
    meanHighCoeff(ii) = mean(nonzeros(ptcoeffLow(:,ii)));
    meanLowCoeff(ii) = mean(nonzeros(ptcoeffHigh(:,ii)));
    
end

%%

figure;
subplot(2,2,1)
histogram(log10(cont14_all))
subplot(2,2,2)
histogram(log10(etoh14_all))
subplot(2,2,3)
histogram(log10(cont17_all))
subplot(2,2,4)
histogram(log10(etoh17_all))

%% get idxs

cont17_idx = contains(animalIDs,'17.5') & contains(animalIDs,'Cont');
cont14_idx = ~contains(animalIDs,'17.5') & contains(animalIDs,'Cont');
etoh17_idx = contains(animalIDs,'17.5') & contains(animalIDs,'EtOH');
etoh14_idx = ~contains(animalIDs,'17.5') & contains(animalIDs,'EtOH');

%%  subplot for embryonic age E14.5
figure;
subplot(1,7,1)

x=categorical({'EtOH E14.5','Cont E14.5'});
x=reordercats(x,{'EtOH E14.5','Cont E14.5'});
bar(x,[mean(means(etoh14_idx)), mean(means(cont14_idx))]);
hold on;
er = errorbar(x,[mean(means(etoh14_idx)), mean(means(cont14_idx))],[std(means(etoh14_idx)) std(means(cont14_idx))]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off;
title('Distribution Differnce of Means');

subplot(1,7,2)

x=categorical({'EtOH E14.5','Cont E14.5'});
x=reordercats(x,{'EtOH E14.5','Cont E14.5'});
bar(x,[mean(mean_ssmd(etoh14_idx)), mean(mean_ssmd(cont14_idx))]);
hold on;
er = errorbar(x,[mean(mean_ssmd(etoh14_idx)), mean(mean_ssmd(cont14_idx))],[std(mean_ssmd(etoh14_idx)) std(mean_ssmd(cont14_idx))]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off;
title('Distribution SSMD');

subplot(1,7,3)

x=categorical({'EtOH E14.5','Cont E14.5'});
x=reordercats(x,{'EtOH E14.5','Cont E14.5'});
bar(x,[mean(mean_overlap(etoh14_idx)), mean(mean_overlap(cont14_idx))]);
hold on;
er = errorbar(x,[mean(mean_overlap(etoh14_idx)), mean(mean_overlap(cont14_idx))],[std(mean_overlap(etoh14_idx)) std(mean_overlap(cont14_idx))]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off;
title('Distribution Overlap');


subplot(1,7,4)

x=categorical({'EtOH E14.5','Cont E14.5'});
x=reordercats(x,{'EtOH E14.5','Cont E14.5'});
bar(x,[mean(meanLowMu(etoh14_idx)), mean(meanLowMu(cont14_idx))]);
hold on;
er = errorbar(x,[mean(meanLowMu(etoh14_idx)), mean(meanLowMu(cont14_idx))],[std(meanLowMu(etoh14_idx)) std(meanLowMu(cont14_idx))]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off;
title('Low Perfusion \mu');


subplot(1,7,5)

x=categorical({'EtOH E14.5','Cont E14.5'});
x=reordercats(x,{'EtOH E14.5','Cont E14.5'});
bar(x,[mean(meanHighMu(etoh14_idx)), mean(meanHighMu(cont14_idx))]);
hold on;
er = errorbar(x,[mean(meanHighMu(etoh14_idx)), mean(meanHighMu(cont14_idx))],[std(meanHighMu(etoh14_idx)) std(meanHighMu(cont14_idx))]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off;
title('High Perfusion \mu');

subplot(1,7,6)

x=categorical({'EtOH E14.5','Cont E14.5'});
x=reordercats(x,{'EtOH E14.5','Cont E14.5'});
bar(x,[mean(meanLowSigma(etoh14_idx)), mean(meanLowSigma(cont14_idx))]);
hold on;
er = errorbar(x,[mean(meanLowSigma(etoh14_idx)), mean(meanLowSigma(cont14_idx))],[std(meanLowSigma(etoh14_idx)) std(meanLowSigma(cont14_idx))]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off;
title('Low Perfusion \sigma');

subplot(1,7,7)

x=categorical({'EtOH E14.5','Cont E14.5'});
x=reordercats(x,{'EtOH E14.5','Cont E14.5'});
bar(x,[mean(meanHighSigma(etoh14_idx)), mean(meanHighSigma(cont14_idx))]);
hold on;
er = errorbar(x,[mean(meanHighSigma(etoh14_idx)), mean(meanHighSigma(cont14_idx))],[std(meanHighSigma(etoh14_idx)) std(meanHighSigma(cont14_idx))]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off;
title('High Pefusion \sigma');
ylabel('Mean \sima Parameter');


%%
figure;
subplot(1,7,1)

x=categorical({'EtOH E17.5','Cont E17.5'});
x=reordercats(x,{'EtOH E17.5','Cont E17.5'});
bar(x,[mean(means(etoh17_idx)), mean(means(cont17_idx))]);
hold on;
er = errorbar(x,[mean(means(etoh17_idx)), mean(means(cont17_idx))],[std(means(etoh17_idx)) std(means(cont17_idx))]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off;
title('Distribution Differnce of Means');

subplot(1,7,2)

x=categorical({'EtOH E17.5','Cont E17.5'});
x=reordercats(x,{'EtOH E17.5','Cont E17.5'});
bar(x,[mean(ptcoeffLow(etoh17_idx)), mean(ptcoeffLow(cont17_idx))]);
hold on;
er = errorbar(x,[mean(ptcoeffLow(etoh17_idx)), mean(ptcoeffLow(cont17_idx))],[std(ptcoeffLow(etoh17_idx)) std(ptcoeffLow(cont17_idx))]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off;
title('pt Coeff');

subplot(1,7,3)

x=categorical({'EtOH E17.5','Cont E17.5'});
x=reordercats(x,{'EtOH E17.5','Cont E17.5'});
bar(x,[mean(ptcoeffHigh(etoh17_idx)), mean(ptcoeffHigh(cont17_idx))]);
hold on;
er = errorbar(x,[mean(ptcoeffHigh(etoh17_idx)), mean(ptcoeffHigh(cont17_idx))],[std(ptcoeffHigh(etoh17_idx)) std(ptcoeffHigh(cont17_idx))]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off;
title('Pt Coeff');


subplot(1,7,4)

x=categorical({'EtOH E17.5','Cont E17.5'});
x=reordercats(x,{'EtOH E17.5','Cont E17.5'});
bar(x,[mean(meanLowMu(etoh17_idx)), mean(meanLowMu(cont17_idx))]);
hold on;
er = errorbar(x,[mean(meanLowMu(etoh17_idx)), mean(meanLowMu(cont17_idx))],[std(meanLowMu(etoh17_idx)) std(meanLowMu(cont17_idx))]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off;
title('Low Perfusion \mu');


subplot(1,7,5)

x=categorical({'EtOH E17.5','Cont E17.5'});
x=reordercats(x,{'EtOH E17.5','Cont E17.5'});
bar(x,[mean(meanHighMu(etoh17_idx)), mean(meanHighMu(cont17_idx))]);
hold on;
er = errorbar(x,[mean(meanHighMu(etoh17_idx)), mean(meanHighMu(cont17_idx))],[std(meanHighMu(etoh17_idx)) std(meanHighMu(cont17_idx))]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off;
title('High Perfusion \mu');

subplot(1,7,6)

x=categorical({'EtOH E17.5','Cont E17.5'});
x=reordercats(x,{'EtOH E17.5','Cont E17.5'});
bar(x,[mean(meanLowAoc(etoh17_idx)), mean(meanLowAoc(cont17_idx))]);
hold on;
er = errorbar(x,[mean(meanLowAoc(etoh17_idx)), mean(meanLowAoc(cont17_idx))],[std(meanLowAoc(etoh17_idx)) std(meanLowAoc(cont17_idx))]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off;
title('Low Perfusion \AOC');

subplot(1,7,7)

x=categorical({'EtOH E17.5','Cont E17.5'});
x=reordercats(x,{'EtOH E17.5','Cont E17.5'});
bar(x,[mean(meanHighAoc(etoh17_idx)), mean(meanHighAoc(cont17_idx))]);
hold on;
er = errorbar(x,[mean(meanHighAoc(etoh17_idx)), mean(meanHighAoc(cont17_idx))],[std(meanHighAoc(etoh17_idx)) std(meanHighAoc(cont17_idx))]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off;
title('High Pefusion \AOC');
xlabel('Mean \Sigma Parameter');
