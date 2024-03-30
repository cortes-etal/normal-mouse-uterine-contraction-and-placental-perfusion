
% this code is used to loop through voxels around a suspected aif to see if
% there is a best option to use

clear all;clc;close all;% start pretty

%%

files = dir('**/*_v2.mat');

%%
for ii = 5%1:numel(files)
    fname=fullfile(files(ii).folder,files(ii).name);
    load(fname)
    
    [~,argmax]=max(ims,[],4); %time of peak
    implay([enhancement/max(enhancement(:)); argmax/size(ims,4)])
    slice=input('Which slice? ');
    
    h=figure;
    imagesc(enhancement(:,:,slice)),colormap(jet(256));
    AIF=imfreehand;
    AIFmask = createMask(AIF);
    close(h)
    
    aif_idx = find(AIFmask);
    for jj=1:numel(aif_idx)
        img_null= zeros(size(AIFmask));
        img_null(aif_idx(jj)) = 1;
        
        AIF=squeeze(sum(sum(ims(:,:,slice,:).*img_null,1),2)/sum(img_null(:)));
        AIFmax=max(AIF)-mean(AIF(1:4));
        [~,aifm_idx] = max(AIF);
        
        figure(111);
        plot(AIF);
        hold on;
        plot(aifm_idx,max(AIF),'r+','MarkerSize',8);
        pause;
    end
    hold off;
    
    
    
    
end