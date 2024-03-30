%Devin Cortes
% Dr. Yijen Wu: The Wu Lab
%^ the purpose of this code is to pally the time serioes placental
%segmenataions across the 13th slice of 599-Cont E17.5



clear all;clc; close all; % start pretty

%%

files = dir('**/*_v2.mat');

dceFile = files(5);
opticFlow_file = 'velocityinfo_slice_13_.mat';
segfile = 'DCESegmentation_Slide13_DWedit.nii.gz';


load(fullfile(dceFile.folder,dceFile.name));
load(fullfile('opticflow_movies17-May-2022','599 Cont-PSFE17.5 (08-26-21)\',opticFlow_file));

sDat = niftiread(segfile);

[X,Y] = meshgrid(0:1:255, 0:1:159);

for jj = 1:size(Vxs,3)
    disp_field(:,:,jj) = divergence(X,Y,Vxs(:,:,jj), Vys(:,:,jj));
end

labels = unique(sDat);
labels(1) = []; % remove background


%% ADD PLACENTAL COMPARTMENT SEGMENTATION HERE....... HOW TO DO THAT......
% FIRST LOOP THROUGH THE SLICE DATA, AND APPLY APRIORI THRESHOLD, BASED
% UPON THE RAW AU INTENSIT VALUES?
ageFlag=1

for ii = 1:size(sDat,3)
    
    perfSlice = ims(:,:,13,ii);
    segSlice = sDat(:,:,ii);
    
    label_one = segSlice == 1;
    label_two = segSlice == 2;
    
    applied_seg_one = perfSlice .* label_one;
    applied_seg_two = perfSlice .* label_two;
    
    if ageFlag
        perfThresh_one = prctile(applied_seg_one(applied_seg_one > 0) ,55); % E17.5 threshold
        perfThresh_two = prctile(applied_seg_two(applied_seg_two > 0) ,55); % E17.5 threshold
        
    else
        perfThresh = prctile(applied_seg(applied_seg > 0,41)); % E14.5 threshold
    end
    
    high_map_one(:,:,ii) = applied_seg_one > perfThresh_one;
    low_map_one(:,:,ii) = applied_seg_one < perfThresh_one;
    
    
    high_map_two(:,:,ii) = applied_seg_two > perfThresh_two;
    low_map_two(:,:,ii) = applied_seg_two < perfThresh_two;
    %
    %     subplot(1,2,1);
    %     imshow(high_map,[]);
    %     subplot(1,2,2);
    %     imshow(low_map,[]);
    %     pause;
    
    
end

high_maps{1} = high_map_one;
high_maps{2} = high_map_two;

low_maps{1} = low_map_one;
low_maps{2} = low_map_two;



%%
slice_13 =squeeze(ims(:,:,13,:));
mag = sqrt(Vxs.^2 + Vys.^2);

for ii = 1:numel(labels)
    
    pmask = sDat == labels(ii);
    
    Vxseg = Vxs .* pmask;
    Vyseg = Vys .* pmask;
    dSeg = disp_field .* pmask;
    perfSeg = slice_13 .* pmask;
    magSeg = mag .* pmask;
    
    Vxseg_h = Vxs .* high_maps{ii};
    Vyseg_h = Vys .* high_maps{ii};
    dSeg_h = disp_field .* high_maps{ii};
    perfSeg_h = slice_13 .* high_maps{ii};
    magSeg_h = mag .*high_maps{ii};
    
    Vxseg_l = Vxs .* low_maps{ii};
    Vyseg_l = Vys .* low_maps{ii};
    dSeg_l = disp_field .* low_maps{ii};
    perfSeg_l = slice_13 .* low_maps{ii};
    magSeg_l = mag .*low_maps{ii};
    
    
    for jj = 1:size(Vxseg,3)
        thresh = prctile(ims(:,:,13,jj),95);
        
        vx_slice = Vxseg(:,:,jj);
        vy_slice = Vyseg(:,:,jj);
        disp_slice = dSeg(:,:,jj);
        perfSlice = perfSeg(:,:,jj);
        magSlice = magSeg(:,:,jj);
        
        vx_slice_h = Vxseg_h(:,:,jj);
        vy_slice_h = Vyseg_h(:,:,jj);
        disp_slice_h = dSeg_h(:,:,jj);
        perfSlice_h = perfSeg_h(:,:,jj);
        magSlice_h = magSeg_h(:,:,jj);
        
        vx_slice_l = Vxseg_l(:,:,jj);
        vy_slice_l = Vyseg_l(:,:,jj);
        disp_slice_l = dSeg_l(:,:,jj);
        perfSlice_l = perfSeg_l(:,:,jj);
        magSlice_l = magSeg_l(:,:,jj);
        
        vel_x(jj) = mean(vx_slice(vx_slice > 0));
        vel_y(jj) = mean(vy_slice(vy_slice > 0));
        disp_vec(jj) =mean(disp_slice(disp_slice > 0));
        perf_curve(jj) = mean(perfSlice(perfSlice > 0));
        magCurve(jj) = mean(magSlice(magSlice > 0));
        
        vel_x_h(jj) = mean(vx_slice_h(vx_slice_h > 0));
        vel_y_h(jj) = mean(vy_slice_h(vy_slice_h > 0));
        disp_vec_h(jj) =mean(disp_slice_h(disp_slice_h > 0));
        perf_curve_h(jj) = mean(perfSlice_h(perfSlice_h > 0));
        magCurve_h(jj) = mean(magSlice_h(magSlice_h > 0));
        
        vel_x_l(jj) = mean(vx_slice_l(vx_slice_l > 0));
        vel_y_l(jj) = mean(vy_slice_l(vy_slice_l > 0));
        disp_vec_l(jj) =mean(disp_slice_l(disp_slice_l > 0));
        perf_curve_l(jj) = mean(perfSlice_l(perfSlice_l > 0));
        magCurve_l(jj) = mean(magSlice_l(magSlice_l > 0));
        
    end
    
    figure(ii);
    subplot(1,5,1);
    plot(vel_x);
    hold on;
    plot(vel_x_h);
    plot(vel_x_l);
    hold off;
    legend('Whole','HP Chamber','LP Chamber');
    title('Velocity in X');
    
    subplot(1,5,2);
    plot(vel_y);
    hold on;
    plot(vel_y_h);
    plot(vel_y_l);
    hold off;
    legend('Whole','HP Chamber','LP Chamber');
    title('Velocity in Y');
    
    subplot(1,5,3);
    plot(disp_vec);
    hold on;
    plot(disp_vec_h);
    plot(disp_vec_l);
    hold off;
    legend('Whole','HP Chamber','LP Chamber');
    title('Displacment');
    
    subplot(1,5,4);
    plot(perf_curve);
    hold on;
    plot(perf_curve_h);
    plot(perf_curve_l);
    hold off;
    legend('Whole','HP Chamber','LP Chamber');
    title('Perfusion');
    
    subplot(1,5,5);
    plot(magCurve);
    hold on;
    plot(magCurve_h);
    plot(magCurve_l);
    hold off;
    legend('Whole','HP Chamber','LP Chamber');
    title('Velocity Magnitude');
    
end

%% 

curves_to_export = [time; perf_curve; perf_curve_h; perf_curve_l; AIF.'; perf_curve./max(perf_curve); ...
    perf_curve_h ./ max(perf_curve_h); perf_curve_l ./ max(perf_curve_l); AIF.'./max(AIF)];

xlswrite('599_e17_placental_cuvres.xlsx',curves_to_export);



%% trying out some psd and spectrogram analysis



figure;
plot(disp_vec);
title('displacemeent spectrogram placent two: animal 599 E175 Control');

fs=1/mean(diff(ACQ_abs_time));
wind= 50;
Nover = 40;
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
[stft,f,t] = spectrogram(disp_vec,wind,Nover,NFFT,fs);
stftDB_whole = 10*log10(abs(stft).^2);
figure(111);
imagesc(t,f,stftDB_whole, [-50 10]);colormap(jet);colorbar;set(gca,'YDir','normal');
title('spectrogram');
ylabel('Frequency [Hz]');
xlabel('Time [s]');
