% DEvin Cortes
%this code will load previously created quantitative perfusion maps,
%present the user with the DCE image to segment some placentas, and then
%analyze the different placental compartments: high v low perfusion

clear all;close all;clc; % start pretty


%%


files = dir('**/*_v2.mat');
testSig = zeros(50,1);
figure;
%%
for ii = 14%1:numel(files)
    fname = fullfile(files(ii).folder,files(ii).name);
    fprts =regexp(fname,'\','split');
    animal=[fprts{4} fprts{5}];
    mkdir(fullfile('PerfusionMaps',animal));
    load(fname);
    
    implay([enhancement/max(enhancement(:))])
    slice=input('Select the Slice of the image that has most placentas\n');
    numPlacentas = input('How many placentas to segment?\n');
    
    h=figure(1);
    [combinedMasks1,combinedLabels] = deal(zeros(size(squeeze(enhancement(:,:,slice)))));
    for jj = 1:numPlacentas
        figure(1);
        imshow(enhancement(:,:,slice),[]);
        placenta =imfreehand;
        placentaMask = createMask(placenta);
        masks{jj} = placentaMask.*jj;
        combinedMasks1 = imadd(combinedMasks1,double(placentaMask));
        combinedLabels = imadd(combinedLabels, double(masks{jj}));
    end
    
    figure;
    subplot(2,1,1)
    imshow(combinedMasks1,[]);
    subplot(2,1,2)
    imshow(label2rgb(combinedLabels));
    
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
    %%
    figure;
    imshow(finalIm)
    
    figure;
    histogram(slopeMap2(slopeMap2 >0))
    
    perfThresh = mean(slopeMap2(slopeMap2 > 0));
    
    lowPerfMask = slopeMap2 < perfThresh;
    highPerfMask = slopeMap2 > perfThresh;
    slopeLow = slopeMap2.*lowPerfMask;
    slopeHigh = slopeMap2.*highPerfMask;
    
    mapLow = double2rgb(slopeLow,cool(256));%,[0 prctile(slopeLow(:),95)]);
    mapHigh = double2rgb(slopeHigh,hot(256));%,[0 prctile(slopeHigh(:),95)]);
    
    lowR =mapLow(:,:,1);
    lowG = mapLow(:,:,2);
    lowB = mapLow(:,:,3);
    
    highR = mapHigh(:,:,1);
    highG = mapHigh(:,:,2);
    highB = mapHigh(:,:,3);
    
    noutR = greydouble(:,:,1);
    noutG = greydouble(:,:,2);
    noutB = greydouble(:,:,3);
    
    noutR(slopeLow > 0) = lowR(slopeLow > 0);
    noutR(slopeHigh > 0) = highR(slopeHigh >0);
    
    noutG(slopeLow >0) = lowG(slopeLow > 0);
    noutG(slopeHigh > 0) = highG(slopeHigh > 0);
    
    noutB(slopeLow > 0) = lowB(slopeLow >0);
    noutB(slopeHigh > 0) = highB(slopeHigh >0);
    
    combinedPerfMaps = cat(3,noutR,noutG,noutB);
    
    %%
    figure
    ax1=axes;
    ax1.Visible = 'off';
    ax2=axes;
    
    linkaxes([ax1,ax2]);
    ax2.XTick = [];
    ax2.YTick = [];
    imagesc(combinedPerfMaps);axis('image');
    
    colormap(ax1,cool(256));
    caxis(ax1,[min(min(slopeLow(slopeLow >0))) max(max(slopeLow(slopeLow>0)))]);
    colormap(ax2,hot(256));
    caxis(ax2,[min(min(slopeHigh(slopeHigh >0))) max(max(slopeHigh(slopeHigh>0)))]);
    
    cb1=colorbar(ax1);
    cb2=colorbar(ax2);
    
    title('Perfusion Map - Placental Chambers');
    
    %%
    meanHigh = mean(slopeHigh(slopeHigh > 0));
    meanLow = mean(slopeLow(slopeLow > 0));
    stdHigh =std(slopeHigh(slopeHigh > 0));
    stdLow =std(slopeLow(slopeLow > 0));
    
    arr=[meanHigh stdHigh; meanLow stdLow];
    
    writetable(array2table(arr),fullfile('PerfusionMaps',animal,['perfDat_',date,'_.xlsx']));
    saveas(gcf,fullfile('PerfusionMaps',animal,['perfMap_',date,'_.png']),'png');
    save(fullfile('PerfusionMaps',animal,['perfMap_',date,'_.mat']));
end


%% evaluating maps

maps = dir("PerfusionMaps\**\*.xlsx")

[counter1 counter2 counter3 counter4] = deal(1);
for ii = 1:numel(maps)
    fname=fullfile(maps(ii).folder,maps(ii).name);
    dat = xlsread(fname);
    
    disp([fname,' : ', num2str(dat(1,1)) ])
    if contains(fname,'EtOH')
        if contains(fname,'E17')
            highEOH_17(counter1) = dat(1,1);
            lowEOH_17(counter1) = dat(2,1);
            counter1=counter1+1;
        else
            highEOH_14(counter2) = dat(1,1);
            lowEOH_14(counter2) = dat(2,1);
            counter2=counter2+1;
        end
        
    else
        if contains(fname,'E14')
            highCont_17(counter3) = dat(1,1);
            lowCont_17(counter3) = dat(2,1);
            counter3=counter3+1;
        else
            highCont_14(counter4) = dat(1,1);
            lowCont_14(counter4) = dat(2,1);
            counter4=counter4+1;
        end
    end
end
%%
figure;
bar([mean(lowCont_14) mean(lowEOH_14) mean(highCont_14) mean(highEOH_14) ...
    mean(lowCont_17) mean(lowEOH_17) mean(highCont_17) mean(highEOH_17)]);
% ['Low 14 Cont', 'Low 14 EtOH','High 14 Cont','High 14 EtOH' ...
%    'Low 17 Cont', 'Low 17 EtOH', 'High 17 Cont', 'High 17 EtOH'], ...

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
