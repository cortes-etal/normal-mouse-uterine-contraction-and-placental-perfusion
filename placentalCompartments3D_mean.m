% DEvin Cortes
%this code will load previously created quantitative perfusion maps,
%present the user with the DCE image to segment some placentas, and then
%analyze the different placental compartments: high v low perfusion

clear all;close all;clc; % start pretty


%%


files = dir('**/*_v2.mat');
testSig = zeros(50,1);
%%
segs = dir('*.nii.gz');
segNames = extractfield(segs,'name');

outDir = ['compartment3D_meanperf' ,date] ;
mkdir(outDir)

figure;
[counter counter2 counter3 counter4] =deal(1);
[etoh14_all cont14_all etoh17_all cont17_all] = deal([]);
count=1;
%%
for ii = 6%1:numel(files)
    
    %getting the animal IDs from directory to sort of the segmentaton files
    fprts = regexp(files(ii).folder,'\','split');
    animalID = fprts{4};
    animalE = fprts{5};
    saveDir = [fprts{4} fprts{5}];
    mkdir(fullfile(outDir,saveDir));
    animalIdx = regexp(animalID, '[\d][\d][\d]');
    animalID = animalID(animalIdx:animalIdx+2);
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
    if ii == 6
       segID(end-1) = 0;
    end
    
    sFile = fullfile(segs(segID).folder,segs(segID).name); % fname of seg file
    sDat = niftiread(sFile); % load the segmentation file % loading seg fle.
    
    fname = fullfile(files(ii).folder,files(ii).name);
    fprts2 =regexp(fname,'\','split');
    animal=[fprts2{4} fprts2{5}];
    %     mkdir(fullfile('PerfusionMaps',animal));
    load(fname);
    %%
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
    %%
    arrf=[];
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
        
        
        
        arr=[l;meanHigh; stdHigh; meanLow; stdLow; meanWhole; stdWhole ;highPerfVol; lowPerfVol; gFlag];
        arrf=[arrf arr];
    end
    %%
    rowNames = {'placenta_num','mean_perf_high','std_perf_high','mean_perf_low','std_perf_low' ...
        'mean_perf_whole','std_perf_whole','highperf_vol','lowperf_vol','group'};
    writetable(array2table(arrf,'RowNames',rowNames),fullfile(outDir,saveDir,[animal,date,'_.xlsx']));
    saveas(gcf,fullfile(outDir,saveDir,[animal,date,'_.png']),'png');
    save(fullfile(outDir,saveDir,['perfMap_',date,'_.mat']));
    
end

%%
figure;
scatter(ones(size(lowCont_14)),lowCont_14,50,'filled');
hold on;
scatter(3*ones(size(lowEOH_14)),lowEOH_14,50,'filled');
scatter(5*ones(size(highCont_14)),highCont_14,50,'filled');
scatter(7*ones(size(highEOH_14)),highEOH_14,50,'filled');
scatter(9*ones(size(lowCont_17)),lowCont_17,50,'filled');
scatter(11*ones(size(lowEOH_17)),lowEOH_17,50,'filled');
scatter(13*ones(size(highCont_17)),highCont_17,50,'filled');
scatter(15*ones(size(highEOH_17)),highEOH_17,50,'filled');
xlim([0 16])
set(gca,'xtick',[]);
ylabel('Perfusion: ml/min/100ml');

%%
figure;
data = [mean(lowCont_14) mean(lowEOH_14) mean(highCont_14) mean(highEOH_14) ...
    mean(lowCont_17) mean(lowEOH_17) mean(highCont_17) mean(highEOH_17)];
bar(1:8,data)

stdData = [std(lowCont_14) std(lowEOH_14) std(highCont_14) std(highEOH_14) ...
    std(lowCont_17) std(lowEOH_17) std(highCont_17) std(highEOH_17)];
errlow=stdData;
errhigh = stdData;

hold on

er = errorbar(1:8,data,errlow,errhigh);
er.Color = [0 0 0];
er.LineStyle = 'none';

hold off

%%
figure;
data2 = [mean(cont14) mean(etoh14) mean(cont17) mean(etoh17)];
errlow2 = [std(cont14) std(etoh14) std(cont17) std(etoh17)];
errhigh2 =errlow2;

bar(1:4,data2);
hold on

er = errorbar(1:4,data2,errlow2,errhigh2);
er.Color = [0 0 0];
er.LineStyle = 'none';

hold off

%% writing data to excel file

lowarray = zeros(4,32);
lowarray(1,1:27) = lowCont_14;
lowarray(2,1:end) = lowEOH_14;
lowarray(3,1:27) = lowCont_17;
lowarray(4,1:24) = lowEOH_17;


higharray = zeros(4,32);
higharray(1,1:27) = highCont_14;
higharray(2,1:end) = highEOH_14;
higharray(3,1:27) = highCont_17;
higharray(4,1:24) = highEOH_17;


wholearray = zeros(4,32);
wholearray(1,1:27) = cont14;
wholearray(2,1:end) = etoh14;
wholearray(3,1:27) = cont17;
wholearray(4,1:24) = etoh17;

tab = array2table([wholearray; higharray; lowarray]);

writetable(tab,['dce_AnalysisWu_',date,'_.xlsx']);


%%

maps = dir("PerfusionMaps\**\*.mat")

[counter1 counter2 counter3 counter4] = deal(1);
for ii = 1:numel(maps)
    fname=fullfile(maps(ii).folder,maps(ii).name);
    load(fname,'slopeHigh','AIFmax');
    maxHigh(ii) = mode(slopeHigh(slopeHigh>0));
    aifMaxs(ii) = AIFmax;
    
end

figure;
scatter(aifMaxs,maxHigh)
