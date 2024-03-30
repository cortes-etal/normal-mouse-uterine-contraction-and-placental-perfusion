% DEvin Cortes
%this code will load previously created quantitative perfusion maps,
%present the user with the DCE image to segment some placentas, and then
%analyze the different placental compartments: high v low perfusion


clear ;close all;clc; % start pretty


%%
files = dir('**/*_v2.mat');
testSig = zeros(50,1);

segs = dir('*.nii.gz');
segNames = extractfield(segs,'name');

outDir = ['compartment_CurveDat_spect_subband_voxel' ,date] ;
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
    
    
    etohFlag = contains(animal,'EtOH');
    flag17 = contains(animal,'17.5');
    
    if etohFlag && flag17
%         etoh17_all = [etoh17_all; slopeMap2dat];
        gFlag = 2
    elseif etohFlag && ~flag17
%         etoh14_all = [etoh14_all; slopeMap2dat];
        gFlag = 4
    elseif ~etohFlag && flag17
%         cont17_all = [cont17_all; slopeMap2dat];
        gFlag = 1
    elseif ~etohFlag && ~flag17
%         cont14_all = [cont14_all; slopeMap2dat];
        gFlag = 3
    end
    
    %% insert placental labeling splitting part
    
    
    singleim = zeros(size(squeeze(ims(:,:,:,1))));
    
    map_rel = zeros(size(singleim));
    map_abs = zeros(size(singleim));

    for kk = 1:numel(singleim)
        if mod(kk,1000) == 0
            disp(kk)
        end
        
        if mod(kk,400000) == 0
            save(['subband_maps_',num2str(kk),'_.mat'],'-v7.3');
        end
           
        tim = singleim;
        tim(kk) = 1;
        wholeCurve=squeeze(sum(sum(sum(ims.*tim,1),2),3)/sum(AIFmask(:)));
        
        wholeCurve = wholeCurve - mean(wholeCurve(1:4));
        
        
        %% Spectrogram Calculation
        
        
        wholeCurve = (wholeCurve - mean(wholeCurve)) ./ rms(wholeCurve);
        Fs = 1/ceil(mean(diff(ACQ_abs_time)));
        fs=Fs;
        wind= 15;
        Nover = 7;
        NFFT =10*25;
        
        %         [stft,f,t] = spectrogram(wholeCurve,wind,Nover,NFFT,Fs);
        %         stftDB_whole = 10*log10(abs(stft).^2);
        
        [SDper_whole, Wper_whole] = periodogram(wholeCurve);
        Fper_whole = fs*Wper_whole./(2*pi);
        SDper_whole(1) = 0;
        
        % doing subbands
        
        totalpower = trapz(Fper_whole,SDper_whole);
        thisband = trapz(Fper_whole(Fper_whole >= 0.007 & Fper_whole <= 0.012),SDper_whole(Fper_whole >= 0.007 & Fper_whole <= 0.012));
        map_rel(kk) = thisband/totalpower;
        map_abs(kk) = thisband;
        
        
        
    end
    
end
% displaying and plotting averages
% %%
% subplot(1,3,1);
% avg_high17 = mean(allstft_high(:,:,(groupidx == 2 | groupidx == 2)),3);
% imagesc(t,f,avg_high17);colormap(jet);colorbar;set(gca,'YDir','normal');
% title('Spectogram of Low Perfusion Curve');
% ylabel('Frequency [Hz]');
% xlabel('Time [s]');
% set(gca,'FontWeight','bold');
%
% %
% %%
% variables = {'peakLow','peakHigh','peakWhole','ttpLow','ttpHigh','ttpWhole','volumeLow','volumeHigh','volume Whole','meanWhole','meanLow','meanHigh','stdWhole','stdLow','stdHigh','wholeHist','lowHist','highHist','litter','group'};
% cell2_Write = cat(1,variables,finalCell);
% writecell(cell2_Write,fullfile(outDir,['allplacentalData_newThresh_2',date,'_.csv']))


