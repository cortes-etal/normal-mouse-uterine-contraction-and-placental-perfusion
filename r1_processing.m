
% goal of script is to first load R1 maps from the DCE animals and see if
% the images are hopefully the same fucking size as the segmentations;


clear all;clc;close all % strat pretty


%% gettting file data


% post r1_maps

post_r1 = dir(fullfile('R1Maps_post\**\*.mat'));

% pre gd r1 maps

pre_r1 = dir(fullfile('R1Maps_pre\**\*.mat')); % dont think i actually need this anymore?


segs = dir('R1map_seg\**\*.nii.gz');
segNames = extractfield(segs,'name');

outDir = ['r1_processing' ,date] ;
mkdir(outDir);


ratios = [0.024336^3 0.03042^3];


%% doing analysis


for ii= 5:numel(post_r1)
    
    post_gd_path = fullfile(post_r1(ii).folder,post_r1(ii).name)
    
    pre_gd_path = fullfile(pre_r1(ii).folder,pre_r1(ii).name)
    r1_map_post = load(post_gd_path,'R1map');
    r1_map_pre = load(pre_gd_path,'R1map');
    
    r1_map_post = r1_map_post.R1map;
    r1_map_pre = r1_map_pre.R1map;
    
    fname = fullfile(post_r1(ii).folder,post_r1(ii).name);
    %%
    
    fparts=regexp(fname,'\','split');
    
    animalID = fparts{5}
    
    idParts = regexp(animalID,' ','split')
    
    id_num = idParts{1}
    
    ageFlag = contains(animalID,'17')
    mkdir(fullfile(outDir,animalID));
    
    if contains(animalID, '17')
        seg_idx = contains(segNames,id_num) & contains(segNames,'17');
        
    else
        seg_idx = contains(segNames,id_num) & ~contains(segNames,'17');
    end
    
    %%
    if sum(seg_idx) == 0
        continue
    else
        
        segdat = niftiread(fullfile(segs(seg_idx).folder,segs(seg_idx).name));
        
        %%
        %  seg_256 = imresize3(segdat, [256 256 20]);
        seg_256 = segdat;
        R1map = r1_map_post;%-r1_map_pre; % uncomment to do change in R1
        
        r1_diff = r1_map_post-r1_map_pre;
        
        R1map(sign(R1map) ==-1) = 0;
        R1map(R1map < prctile(R1map(:),5)) = 0;
        
        
        r1_diff(sign(r1_diff) ==-1) = 0;
        r1_diff(r1_diff < prctile(r1_diff(:),5)) = 0;
        %         figure(1);
        %         imshow3D(R1map, [ 0 prctile(R1map(:), 95)]);colormap(jet(256));
        %
        %         figure(2);
        %         imshow3D(R1map .* double(seg_256>0), [0 prctile(R1map(:), 95)]);colormap(jet(256));
        %
        %         figure(3);
        %         histogram(nonzeros(R1map(:)));pause;
        %
        %
        [rgb,frgb, r1Low, r1High] = combineRGB_nD(R1map, R1map .* double(seg_256 > 0),animalID,animalID,ageFlag);
        
        [~,~, r1Low_diff, r1High_diff] = combineRGB_nD(r1_diff, r1_diff .* double(seg_256 > 0),animalID,animalID,ageFlag);
        %
        
        p_nums = unique(seg_256);
        p_nums(1) = []; % get rid of backgound label
        
        for pp = 1:numel(p_nums)
            placenta_label = p_nums(pp);
            seg_pp = seg_256 == placenta_label;
            
            r1_low(pp) = mean(nonzeros(r1Low .* seg_pp));
            r1_high(pp) = mean(nonzeros(r1High .* seg_pp));
            r1_whole(pp) = mean(nonzeros(R1map .* seg_pp));
            
            r1_low_diff(pp) = mean(nonzeros(r1Low_diff .* seg_pp));
            r1_high_diff(pp) = mean(nonzeros(r1High_diff .* seg_pp));
            r1_whole_diff(pp) = mean(nonzeros(r1_diff .* seg_pp));
            
            r1_low_v(pp) = numel(nonzeros(r1Low .* seg_pp));
            r1_high_v(pp) = numel(nonzeros(r1High .* seg_pp));
            r1_whole_v(pp) = numel(nonzeros(R1map .* seg_pp));
            
        end
        
        all_r1 = [r1_low; r1_high; r1_whole];
        
        all_r1_diff = [r1_low_diff; r1_high_diff; r1_whole_diff];
        
        all_r1_vol = [r1_low_v; r1_high_v; r1_whole_v];
        
        if ageFlag
            all_gd = [r1_low_diff .* (r1_low_v .* ratios(2)); ...
                r1_high_diff .*(r1_high_v .* ratios(2)); ...
                r1_whole_diff .*(r1_whole_v .* ratios(2))];
            
        else
            all_gd = [r1_low_diff .* (r1_low_v .* ratios(1)); ...
                r1_high_diff .*(r1_high_v .* ratios(1)); ...
                r1_whole_diff .*(r1_whole_v .* ratios(1))];
        end
        
        final_out = [all_r1 all_r1_diff all_gd];
        
        
%         xlswrite(fullfile(outDir,['r1_summary_data_',date,'_.xlsx']),final_out,animalID);
        
        clear r1_low r1_high r1_whole
        clear r1_low_v r1_high_v r1_whole_v
        clear r1_low_diff r1_high_diff r1_whole_diff
    end
    %     R1map(R1map > 15) = 0;
    %     niftiwrite(R1map,fullfile(outDir,animalID,[animalID ,'r1_map_.nii']));
    
end


%% running of cooreigistration between pre gd t1 and post gd t1 (this shit didnt work)
% Options.Registration='affine';
% Options.Penalty=1e5;
% for ii = 4%:numel(post_r1)
%     post_gd_path = fullfile(post_r1(ii).folder,post_r1(ii).name)
%
%
%     pre_gd_path = fullfile(pre_r1(ii).folder,pre_r1(ii).name)
%
%     if contains(post_gd_path,'EtOH')
%         continue
%     else
%         r1_map_post = load(post_gd_path,'R1map');
%         r1_map_pre = load(pre_gd_path,'R1map');
%
%         r1_map_post = r1_map_post.R1map;
%         r1_map_pre = r1_map_pre.R1map;
%
%         [premap_post_reg,~,~,M]= image_registration(rescale(r1_map_post(:,:,10)),rescale(r1_map_pre(:,:,10)),Options);
%     end
%
%
%
% end
%
% imshow(premap_post_reg,[0 prctile(premap_post_reg(:),95)]);