% Devin Cortes
% Dr. Yijen Wu: The Wu Lab
%Purpose if this code is to pally time-series segmentations of the moiuse
%placenta to DCE MRI time series images. the data the segmentations are
%being applied to are optical flow maps as well as dipslacement and
%perfusion maps.

clear all;clc;close all;; % start pretty

%% get and organize directories containing segmentations and data

outDir = ['placental_timeseries_velocity_analysis_',date];
mkdir(outDir);
seg_dir = 'timeSeriesSegmentations';

opticflow_dir = 'opticflow_movies17-May-2022';

dce_files = dir('**/*_v2.mat');
dce_folders = extractfield(dce_files,'folder');

time_series_seg=dir(fullfile(seg_dir,'*.nii.gz'));

opticflow_mat = dir(fullfile(opticflow_dir,'**','*.mat'));

flow_dirs = extractfield(opticflow_mat,'folder');
slice_names= extractfield(opticflow_mat,'name');
% looping over the segmentations (leeast amount of loops) ...
% and trying to automatically pull the correct info from the other folder
% dirs
[X,Y] = meshgrid(0:1:255, 0:1:159); % matrix for displacement field calclation
allplac_trapz_combined =[];
%  begin looping over datasets
for ii = 10%1:numel(time_series_seg)
    seg_path = fullfile(time_series_seg(ii).folder,time_series_seg(ii).name);
    seg_file = time_series_seg(ii).name;
    
    seg_parts = regexp(seg_file,' ','split');
    
    animalId = seg_parts{1};
    age_group_str = seg_parts{2};
    
    if contains(age_group_str,'EtOH')
        continue % skkip the PCA animals
    else
        pattern = {'[1][457][.][5]','slice_[\d][\d]','Slice [\d][\d]'};%,'Slice [\d]','slice_[\d]'}
        identifiers = regexp(seg_file,pattern,'match'); % get segmentation identifiers for age and image slice
        non_empty_cells = ~cellfun(@isempty,identifiers); % get te non empty ones (ie the different pattern optiosn)
        ids = identifiers(non_empty_cells); % LOAD IDS
        age=ids{1}; % separete age
        slice = ids{2}; % separete slice id
        
        ageFlag = contains(age,'17'); % age Flag for threshold selection
        if non_empty_cells(3) == 1 % reformat slice id for proper .mat selection
            slice = ['slice_',slice{1}(end-1:end)];
        else 
            slice=slice{1};
        end
        slice_num = str2double(slice(end-1:end));
        
        flow_dir_idx = contains(flow_dirs,age) & contains(flow_dirs,animalId) & contains(slice_names,slice); % find the proper .mat file based on all the identifiers we pulled
        dce_dir_dix = contains(dce_folders,age) & contains(dce_folders,animalId);
        
        load(fullfile(dce_files(dce_dir_dix).folder,dce_files(dce_dir_dix).name)); % load dce information
        load(fullfile(opticflow_mat(flow_dir_idx).folder,opticflow_mat(flow_dir_idx).name)); % load velocity inforamtion
        segdat = double(niftiread(seg_path)); % read segmentatio file
        
        fs = 1/mean(diff(ACQ_abs_time));
        time = 0:(1/fs):(numel(ACQ_abs_time)*(1/fs) - (1/fs));
        % load dce .mat information
        
        for kk = 1:size(Vxs,3)
            disp_field(:,:,kk) = divergence(X,Y,Vxs(:,:,kk), Vys(:,:,kk));
        end
        
        
        mag = sqrt(Vxs.^2 + Vys.^2);
        
        placenta_labels = nonzeros(unique(segdat)); % GET LABELS And remove bacjground label
        this_slice = squeeze(ims(:,:,slice_num,:)); % get kinetics time series for loaded slice
        [allplac_traces, allplac_trapz] = deal([]);
        for jj = 1:numel(placenta_labels)
            sizes = size(segdat);
            pmask = segdat == placenta_labels(jj); % turn multilabel segmenttion into a single mask of ones for a given placenta and zero everywhere else.
            signal_rep_pmask = reshape(pmask,[sizes(1)*sizes(2),sizes(3)]);
            pix_row_ids = find(sum(signal_rep_pmask,2));
            applied_mask = pmask .* squeeze(ims(:,:,slice_num,:)); % apply the placental mask to the entire time series of a given slice and lose the slice dimension
            
            if ageFlag % different thrshp;ds for different emebryhonic agaes
                pthresh = prctile(applied_mask(applied_mask > 0) ,55); % E17.5 threshold
            else
                pthresh = prctile(applied_mask(applied_mask > 0),41); % E14.5 threshold
            end
            
            high_perf = applied_mask >= pthresh;
            low_perf = applied_mask < pthresh;
            
            high_map = high_perf > 0;
            signal_rep_high = reshape(high_map,[sizes(1)*sizes(2),sizes(3)]);
            high_pix_ids =find(sum(signal_rep_high,2));
            
            low_map = low_perf > 0;
            signal_rep_low = reshape(high_map,[sizes(1)*sizes(2),sizes(3)]);
            low_pix_ids =find(sum(signal_rep_high,2));
            
            % calculating velocity trace in x-direction
            vel_x_high = reshape(high_map .* Vxs, [sizes(1)*sizes(2),sizes(3)]);
            mean_x_high = mean(vel_x_high(high_pix_ids,:),1);
            
            vel_x_low = reshape(low_map .* Vxs, [sizes(1)*sizes(2),sizes(3)]);
            mean_x_low = mean(vel_x_low(low_pix_ids,:),1);
            
            vel_x_whole = reshape(pmask .* Vxs ,[sizes(1)*sizes(2),sizes(3)]);
            mean_x_whole = mean(vel_x_whole(pix_row_ids,:),1);
            
            % velocity cuve in y-direction
            vel_y_high = reshape(high_map .* Vys, [sizes(1)*sizes(2),sizes(3)]);
            mean_y_high = mean(vel_y_high(high_pix_ids,:),1);
            
            vel_y_low = reshape(low_map .* Vys, [sizes(1)*sizes(2),sizes(3)]);
            mean_y_low = mean(vel_y_low(low_pix_ids,:),1);
            
            vel_y_whole = reshape(pmask .* Vys ,[sizes(1)*sizes(2),sizes(3)]);
            mean_y_whole = mean(vel_y_whole(pix_row_ids,:),1);
            
            %veloicty mafnitude trace
            vel_mag_high = reshape(high_map .* mag, [sizes(1)*sizes(2),sizes(3)]);
            mean_mag_high = mean(vel_mag_high(high_pix_ids,:),1);
            
            vel_mag_low = reshape(low_map .* mag, [sizes(1)*sizes(2),sizes(3)]);
            mean_mag_low = mean(vel_mag_low(low_pix_ids,:),1);
            
            vel_mag_whole = reshape(pmask .* mag ,[sizes(1)*sizes(2),sizes(3)]);
            mean_mag_whole = mean(vel_mag_whole(pix_row_ids,:),1);
            
           proc_signal_1 =  process_motion_signals(mean_mag_high,fs);
           proc_signal_2 = process_motion_signals(mean_mag_low,fs);
           proc_signal_3 =  process_motion_signals(mean_mag_whole,fs);
            pause;
            
            %velocity perfusion trace
            perf_high = reshape(high_map .* this_slice, [sizes(1)*sizes(2),sizes(3)]);
            mean_perf_high = mean(perf_high(high_pix_ids,:),1);
            
            perf_low = reshape(low_map .* this_slice, [sizes(1)*sizes(2),sizes(3)]);
            mean_perf_low = mean(perf_low(low_pix_ids,:),1);
            
            perf_whole = reshape(pmask .* this_slice ,[sizes(1)*sizes(2),sizes(3)]);
            mean_perf_whole = mean(perf_whole(pix_row_ids,:),1);
            
            
            %velocity divergence trace
            disp_high = reshape(high_map .* disp_field, [sizes(1)*sizes(2),sizes(3)]);
            mean_disp_high = mean(disp_high(high_pix_ids,:),1);
            
            disp_low = reshape(low_map .* disp_field, [sizes(1)*sizes(2),sizes(3)]);
            mean_disp_low = mean(disp_low(pix_row_ids,:),1);
            
            disp_whole = reshape(pmask .* disp_field ,[sizes(1)*sizes(2),sizes(3)]);
            mean_disp_whole = mean(vel_x_whole(pix_row_ids,:),1);
            
            figure(2);
            title('Divergence Traces');
            plot(time/60,mean_disp_high);hold on;
            plot(time/60,mean_disp_low);
            plot(time/60,mean_disp_whole);
            hold off;
            legend('High','Low','Whole');
            pause;
                        
            %acceleration
            acc_x_high = [diff(mean_x_high)]; %appending zeros for easy epxort
            acc_y_high = [diff(mean_y_high)];
            
            acc_x_low = [diff(mean_x_low)];
            acc_y_low = [diff(mean_y_low)];
            
            acc_x_whole = [diff(mean_x_whole)];
            acc_y_whole = [diff(mean_y_whole)];
            
            acc_mag_high = [diff(mean_mag_high)];
            acc_mag_low = [diff(mean_mag_low)];
            acc_mag_whole = [diff(mean_mag_whole)];
            
            % area under curve calculates: vel -> disp
            total_disp_x_high=trapz(time,mean_x_high);
            total_disp_y_high = trapz(time,mean_y_high);
            
            total_disp_x_low = trapz(time,mean_x_low);
            total_disp_y_low = trapz(time,mean_y_low);
            
            total_disp_x_whole = trapz(time,mean_x_whole);
            total_disp_y_whole = trapz(time,mean_y_whole);
            
            total_disp_mag_high = trapz(time,mean_mag_high);
            total_disp_mag_low = trapz(time,mean_mag_low);
            total_disp_mag_whole = trapz(time,mean_mag_whole);
            
            % area under curve: acc -> vel
            total_vel_x_high=trapz(time(1:end-1),acc_x_high);
            total_vel_y_high = trapz(time(1:end-1),acc_y_high);
            
            total_vel_x_low = trapz(time(1:end-1),acc_x_low);
            total_vel_y_low = trapz(time(1:end-1),acc_y_low);
            
            total_vel_x_whole = trapz(time(1:end-1),acc_x_whole);
            total_vel_y_whole = trapz(time(1:end-1),acc_y_whole);
            
            total_vel_mag_high = trapz(time(1:end-1),acc_mag_high);
            total_vel_mag_low = trapz(time(1:end-1),acc_mag_low);
            total_vel_mag_whole = trapz(time(1:end-1),acc_mag_whole);
            
            label_vec = ones(1,50).*jj;
            
            acc_x_high = [diff(mean_x_high) 0]; %appending zeros for easy epxort
            acc_y_high = [diff(mean_y_high) 0];
            
            acc_x_low = [diff(mean_x_low) 0];
            acc_y_low = [diff(mean_y_low) 0];
            
            acc_x_whole = [diff(mean_x_whole) 0];
            acc_y_whole = [diff(mean_y_whole) 0];
            
            acc_mag_high = [diff(mean_mag_high) 0];
            acc_mag_low = [diff(mean_mag_low) 0];
            acc_mag_whole = [diff(mean_mag_whole) 0];
            
            sngplac_traces = [mean_x_whole; mean_x_low; mean_x_high; mean_y_whole; mean_y_low; mean_y_high; ...
                mean_mag_whole; mean_mag_low; mean_mag_high; mean_perf_whole; mean_perf_low; mean_perf_high; ...
                mean_disp_whole; mean_disp_low; mean_disp_high; acc_x_whole; acc_x_low; acc_x_high; ...
                acc_y_whole; acc_y_low; acc_y_high; acc_mag_whole; acc_mag_low; acc_mag_high; ...
                label_vec];
            
            sngplac_trapz = [total_disp_x_whole; total_disp_x_low; total_disp_x_high; total_disp_y_whole; ...
                total_disp_y_low; total_disp_y_high; total_disp_mag_whole; total_disp_mag_low; total_disp_mag_high; ...
                total_vel_x_whole; total_vel_x_low; total_vel_x_high; total_vel_y_whole; total_vel_y_low; ...
                total_vel_y_high; total_vel_mag_whole; total_vel_mag_low; total_vel_mag_high; jj; str2double(animalId)];
            
            allplac_traces = [allplac_traces; sngplac_traces];
            allplac_trapz = [allplac_trapz; sngplac_trapz];
            
        end %end placenta loop
        
%         xlswrite(fullfile(outDir,['velocity_analysis_',date,'_.xlsx']),allplac_traces,[animalId age{1}]);
%         xlswrite(fullfile(outDir,['trapz_data_',date,'_.xlsx']),allplac_trapz,[animalId age{1}])
%         allplac_trapz_combined=[allplac_trapz_combined allplac_trapz];
    end
    
%     xlswrite(fullfile(outDir,['trapz_curves_',date,'_.xlsx']),allplac_trapz_combined);
    
end


%% %ploting some perfusion cruves to see which one i want for paper

for ii = 10:10:75
    figure(ii)
    plot(allplac_traces(ii,:));
    
end