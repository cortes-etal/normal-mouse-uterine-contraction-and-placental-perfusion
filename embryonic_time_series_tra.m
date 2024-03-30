% DEvin Cortes
% University of Pittsburgh, Dept of Bioengineerring, Dept of Dev Bio
% Dr. Yijen WU
% 28Aug2022, copyright Devin Cortes

% purpose of this code is to load example segmentations of mouse embryos
% and apply them to the mri DCE time series images to see if they incrase
% in intensity, suggesting Gd+ flow into the fetus.

close all; clc; clear all; % start pretty

%%

kinetic_scans_dir = 'F:\DCE Analysis\Wu PCA4\kinetic_movies_nifti_exports26-May-2022';

slices_e14_dir = dir(fullfile(kinetic_scans_dir,'600 Cont-PSFE14.5 (08-22-21)','slice_1*.nii')); % struct containing the chosen mri data
slices_e17_dir = dir(fullfile(kinetic_scans_dir,'600 Cont-PSFE17.5 (08-25-21)','slice_1*.nii')); % struct containing the chosen mri data

slice10_dir = dir(fullfile(kinetic_scans_dir,'600 Cont*'));

e14_segs = dir(fullfile(slices_e14_dir(1).folder,'*.nii.gz')); % finding path of e14 seg
e17_segs = dir(fullfile(slices_e17_dir(1).folder,'*.nii.gz')); % finding path of e17 seg



slices_anal = {'10','11','12','13'};
[combined_slices_14, combined_slices_17] = deal([]);
[comb_seg_14, comb_seg_17] = deal([]);
in_counter = 1;
for ii = 1:numel(slices_e14_dir) % this loop loads, oegnaizes, and 4D concats the seg and mir data
    slice_path_14 = fullfile(slices_e14_dir(ii).folder,slices_e14_dir(ii).name); % get the slice mriu path na,e
    slice_path_17 = fullfile(slices_e17_dir(ii).folder,slices_e17_dir(ii).name);
    
    slice_idx_14 = contains(slice_path_14,slices_anal); % does slice contain the right dir
    slice_idx_17 = contains(slice_path_17,slices_anal); % does slice contain the right dir
    
    if sum(slice_idx_14) && sum(slice_idx_17) % if to load data if above is ok
        e14_seg_p = fullfile(e14_segs(in_counter).folder,e14_segs(in_counter).name);
        e17_seg_p = fullfile(e17_segs(in_counter).folder,e17_segs(in_counter).name);
        
        seg_dat_14= double(niftiread(e14_seg_p));
        seg_dat_17= double(niftiread(e17_seg_p));
        
        comb_seg_14 = cat(4,comb_seg_14,seg_dat_14);
        comb_seg_17 = cat(4,comb_seg_17,seg_dat_17);
        
        
        img_path_14 = fullfile(slices_e14_dir(ii).folder,slices_e14_dir(ii).name); % e14 img path
        img_path_17 = fullfile(slices_e17_dir(ii).folder,slices_e17_dir(ii).name); % e17 img path
        
        mri_dat_14  = niftiread(img_path_14); % load mir data
        mri_dat_17  = niftiread(img_path_17);
        
        combined_slices_14 = cat(4,combined_slices_14,mri_dat_14); % concatenate in 4D
        combined_slices_17 = cat(4,combined_slices_17,mri_dat_17);
        
        in_counter=in_counter+1;
        
    else
        continue;
        
    end
    
    
end

%%
labels = unique(comb_seg_14);
labels(1) = []; % remove background label
for tt = 1:2
    for kk = 1:numel(labels)
        label = labels(kk);
        
        
        
        if tt==1
            img_dat = combined_slices_14;
            seg_dat =comb_seg_14;
            title_str = 'E14 Embryo';
        else
            img_dat = combined_slices_17;
            seg_dat =comb_seg_17;
            title_str = 'E17 Embryo';
        end
        
        for zz = 1:size(img_dat,3)
            d_slice = img_dat(:,:,zz,:);
            seg_slice = seg_dat(:,:,zz,:);
            
            perf_trace(zz,kk) = mean(nonzeros(d_slice(seg_slice == label)));
            
        end
        
        
        
    end
    
    perf_1 = perf_trace(:,1);
    perf_2 = perf_trace(:,2);
    perf_1_prc = (perf_1 - mean(perf_1)) ./ std(perf_1);
    perf_2_prc = (perf_2 - mean(perf_2)) ./ std(perf_2);
    
    proc_perfs = cat(2, perf_1_prc, perf_2_prc);
    
    figure(tt);
    subplot(3,1,1)
    plot((perf_trace));
    title(title_str);
    legend('Embryo 1','Embryo 2');
    subplot(3,1,2)
    plot(movmean(perf_trace,15));
    title('5 Sample Moving Average');
    legend('Embryo 1','Embryo 2');
    subplot(3,1,3)
    
    plot((proc_perfs));
    title('Mean Removed, Std Normalized');
    legend('Embryo 1','Embryo 2');
end



%%
%
% for ii = 1:numel(slice10_dir)
%
%     if ii == 1
%         img_path = fullfile(slice10_dir(1).folder,slice10_dir(1).name);
%         seg_path =e14_seg_path;
%         title_str = 'E14 Embryo';
%     else
%         img_path = fullfile(slice10_dir(2).folder,slice10_dir(2).name);
%         title_str = 'E17 Embryo';
%         seg_path =e17_seg_path;
%
%     end
%
%     img_dat = niftiread(img_path);
%     seg_dat = niftiread(seg_path);
%
%     labels = nonzeros(unique(seg_dat)); % get labels and remove bakcground onev (label = 0)
%
%     for zz=1:numel(labels)
%         label = labels(zz);
%
%         for kk = 1:size(img_dat,3)
%             slice = img_dat(:,:,kk);
%             seg_slice = seg_dat(:,:,kk);
%
%             perf_trace(kk,zz) = mean(nonzeros(slice(seg_slice == label)));
%         end
%
%     end
%
%     figure(ii);
%     subplot(2,1,1)
%     plot((perf_trace));
%     title(title_str);
%     legend('Embryo 1','Embryo 2');
%     subplot(2,1,2)
%     plot(movmean(perf_trace,15));
%     title('5 Sample Moving Average');
%     legend('Embryo 1','Embryo 2');
%
% end
%





