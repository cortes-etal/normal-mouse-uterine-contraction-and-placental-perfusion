% DEvin Cortes
%this code will load previously created quantitative perfusion maps,
%present the user with the DCE image to segment some placentas, and then
%analyze the different placental compartments: high v low perfusion


clear all;close all;clc; % start pretty


%%
files = dir('**/*_v2.mat');
testSig = zeros(50,1);

segs = dir('*kidneys.nii');
segNames = extractfield(segs,'name');

outDir = ['compartment_CurveDat_spect_kidneys' ,date] ;
mkdir(outDir);

%%
figure;
[etoh14_all cont14_all etoh17_all cont17_all] = deal([]);
count=1;
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
    
    fname = fullfile(files(ii).folder,files(ii).name);
    fprts2 =regexp(fname,'\','split');
    animalIDs{ii}=[fprts2{4} fprts2{5}];
    
    animalE = animalE(animalEdx:animalEdx+1); % the dataset animal embryonic age
    
    segID = contains(segNames, animalE) & contains(segNames,animalID); % logical vecotr containing idx for seg file
    
    if sum(segID) == 0
        continue
    end
    
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
    
    [rgb,slopeLow,slopeHigh] = combineRGB_nD(enhancement,slopeMap2,['histo_',num2str(ii),'_.png'],animal,animalE);
    
    
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
        
        
        
        wholeCurve=squeeze(sum(sum(sum(ims.*(pMask > 0),1),2),3)/sum(AIFmask(:)));
        
        wholeCurve = wholeCurve - mean(wholeCurve(1:4));
        
        if ii == 13
            
            slopeCurveWhole = wholeCurve(10:end);
            
        end
        %% Spectrogram Calculation
        
        
        wholeCurve = (wholeCurve - mean(wholeCurve)) ./ rms(wholeCurve);
        Fs = 1/ceil(mean(diff(ACQ_abs_time)));
        
        figure;
        plot(wholeCurve);
        title(animal);
        pause;
        fs=Fs;
        wind= 15;
        Nover = 7;
        NFFT =10*25;
        %
        %         [stft,f,t] = spectrogram(highCurve,wind,Nover,NFFT,Fs);
        %         stftDB_high = 10*log10(abs(stft).^2);
        %
        %         [tmpMom_high,f]= tftmoment(highCurve,fs,2);
        %         [spcMon_high,t] = tfsmoment(highCurve,fs,2);
        %
        %         [stft,f,t] = spectrogram(lowCurve,wind,Nover,NFFT,Fs);
        %         stftDB_low = 10*log10(abs(stft).^2);
        %
        %         [tmpMom_low,f]= tftmoment(lowCurve,fs,2);
        %         [spcMon_low,t] = tfsmoment(lowCurve,fs,2);
        %
        [stft,f,t] = spectrogram(wholeCurve,wind,Nover,NFFT,Fs);
        stftDB_whole = 10*log10(abs(stft).^2);
        figure(111);
        imagesc(t,f,stftDB_whole, [-50 10]);colormap(jet);colorbar;set(gca,'YDir','normal');
        title(animal);
        ylabel('Frequency [Hz]');
        xlabel('Time [s]');
  
        pause;
        %         [tmpMom_w,f]= tftmoment(wholeCurve,fs,2);
        %         [spcMon_w,t] = tfsmoment(wholeCurve,fs,2);
        %         %% PSD Estimation
        %         fs=Fs;
        %                 [SDper_high, Wper_high] = periodogram(highCurve);
        %                 Fper_high = fs*Wper_high./(2*pi);
        %                 SDper(1) = 0;
        
        
        
        %
        %         [SDwelch_high, Wwel_high] = pwelch(highCurve);
        %         Fwel_high = fs*Wwel_high./(2*pi);
        %
        %         [SDthom_high, Wthom_high] = pmtm(highCurve);
        %         Fthom_high = fs*Wthom_high./(2*pi);
        %
        %         n =4; %oreder of AR model
        %         [SDar_high, War_high] = pburg(highCurve,n);
        %         Far_high = fs*War_high./(2*pi);
        %
        %         % now for low curve
        %         [SDper_low, Wper_low] = periodogram(lowCurve);
        %         Fper_low = fs*Wper_low./(2*pi);
        %         SDper(1) = 0;
        %
        %         [SDwelch_low, Wwel_low] = pwelch(lowCurve);
        %         Fwel_low = fs*Wwel_low./(2*pi);
        %
        %         [SDthom_low, Wthom_low] = pmtm(lowCurve);
        %         Fthom_low = fs*Wthom_low./(2*pi);
        %
        %         n =4; %oreder of AR model
        %         [SDar_low, War_low] = pburg(lowCurve,n);
        %         Far_low = fs*War_low./(2*pi);
        %
        % noe for wholeCurve
        [SDper_whole, Wper_whole] = periodogram(wholeCurve);
        Fper_whole = fs*Wper_whole./(2*pi);
        SDper_whole(1) = 0;
        
        figure(1222);
        h=plot(Fper_whole,10*log10(abs(SDper_whole)));
        set(h,'LineWidth',2);
        title(animal);
              ylim([-50 5]);
        pause;
        %
        %         [SDwelch_whole, Wwel_whole] = pwelch(wholeCurve);
        %         Fwel_whole = fs*Wwel_whole./(2*pi);
        %
        %         [SDthom_whole, Wthom_whole] = pmtm(wholeCurve);
        %         Fthom_whole = fs*Wthom_whole./(2*pi);
        %
        %         n =4; %oreder of AR model
        %         [SDar_whole, War_whole] = pburg(wholeCurve,n);
        %         Far_whole = fs*War_whole./(2*pi);
        %
        %
        %         %% doing time and frequency marginals
        %
        %
        %
        %         %% calculating averages
        %
        %         tfmoments_l(:,count) = tmpMom_low;
        %         tfmoments_h(:,count) = tmpMom_high;
        %         tfmoments_w(:,count) = tmpMom_w;
        %
        %         specmoments_l(:,count) = spcMon_low;
        %         specmoments_h(:,count) = spcMon_high;
        %         specmoments_w(:,count) = spcMon_w;
        %
        %         allstft_low(:,:,count) =stftDB_low;
        %         allstft_high(:,:,count) =stftDB_high;
        %         allstft_whole(:,:,count) =stftDB_whole;
        %
        %         allSDper_low(:,count) = SDper_low;
        %         allSDwelech_low(:,count) = SDwelch_low;
        %         allSDthom_low(:,count) = SDthom_low;
        %         allSDar_low(:,count) = SDar_low;
        %
        %         allSDper_h(:,count) = SDper_high;
        %         allSDwelech_h(:,count) = SDwelch_high;
        %         allSDthom_h(:,count) = SDthom_high;
        %         allSDar_h(:,count) = SDar_high;
        %
        %         allSDper_w(:,count) = SDper_whole;
        %         allSDwelech_w(:,count) = SDwelch_whole;
        %         allSDthom_w(:,count) = SDthom_whole;
        %         allSDar_w(:,count) = SDar_whole;
        %
        %         groupidx(count) = gFlag;
        %         count = count +1;
    end
    
end
% displaying and plotting averages
%%
% subplot(1,3,1);
% avg_high17 = mean(allstft_high(:,:,(groupidx == 2 | groupidx == 2)),3);
% imagesc(t,f,avg_high17);colorbar;set(gca,'YDir','normal');
% title('Spectogram of Low Perfusion Curve');
% ylabel('Frequency [Hz]');
% xlabel('Time [s]');
% set(gca,'FontWeight','bold');

%
% %%
% variables = {'peakLow','peakHigh','peakWhole','ttpLow','ttpHigh','ttpWhole','volumeLow','volumeHigh','volume Whole','meanWhole','meanLow','meanHigh','stdWhole','stdLow','stdHigh','wholeHist','lowHist','highHist','litter','group'};
% cell2_Write = cat(1,variables,finalCell);
% writecell(cell2_Write,fullfile(outDir,['allplacentalData_newThresh_2',date,'_.csv']))


