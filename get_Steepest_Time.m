clear all;
close all hidden;

mainpath='C:\Users\ddblabwu\Documents\MATLAB\MATLAB\code\code\Anthony\Anthony'; %PUT PATH TO THIS FOLDER HERE
addpath(fullfile(mainpath,'bruker_io'));

%% Find IntraGate datasets

numAnimals=4
%for control = 1:numAnimalse
%     clearvars -except numAnimals
%     close all
    folders = uigetfile_n_dir;
    prts=regexp(folders{1},'\','split');
    mainDir = strcat(prts{4},prts{5});
    
    k=1;
    r3=@(x)reshape(x,1,1,[]);
    for j = 1:numel(folders)
        j
        foldername=folders{j};
        
        % Load IntraGate datasets
        %  k=1;
        acqp{j}=read_bruker_acqp(fullfile(foldername,'acqp'));
        %for j=find(isIntraGate)
        try
            reco=read_bruker_acqp(fullfile(foldername,'pdata/1/reco'));
        catch %try again in case Box just isn't responsive the first time
            reco=read_bruker_acqp(fullfile(foldername,'pdata/1/reco'));
        end
        im=load_bruker_2dseq(reco.RECO_size(1),fullfile(foldername,'pdata/1/2dseq'),'int32');
        im=reshape(im,reco.RECO_size(2),reco.RECO_size(1),[]);
        try
            ims(:,:,1:size(im,3),k)=bsxfun(@rdivide,im,r3(reco.RECO_map_slope)); %store the image, correcting for the storage scaling factor
            ACQ_time{1,k}=acqp{j}.ACQ_time;
            ACQ_abs_time(1,k)=acqp{j}.ACQ_abs_time;
            k=k+1;
        catch
            break;
        end
        
    end
    
    % deltat = median(diff(ACQ_abs_time))/60; %in minutes
    deltat = diff(ACQ_abs_time);
    [~,argmax]=max(deltat);
    deltat(argmax)=[]; %remove single long delay
    deltat = mean(deltat)/60; %in minutes
    
    
    
    %% Steepest slope
    
    % closed-form solution for three-point linear fit is just (last - first)/2:
    %syms x1 x2 x3;
    %temp=[x1 x2 x3]*pinv([-1 0 1; 1 1 1]);
    %slope=temp(1);
    
    % remove signal drops from motion along time dimension ("morphological closing")
    processed=max(cat(5,ims(:,:,:,1:end-2),ims(:,:,:,2:end-1),ims(:,:,:,3:end)),[],5);
    processed=min(cat(5,processed(:,:,:,1:end-2),processed(:,:,:,2:end-1),processed(:,:,:,3:end)),[],5);
    processed=cat(4,ims(:,:,:,1:2),processed,ims(:,:,:,end-1:end));
    
    % spatial median filter (light smoothing, also helpful with motion)
    processed=medfilt3(reshape(padarray(processed,[0 0 1 0]),size(ims,1),size(ims,2),[]),[3 3 3]);
    processed=reshape(processed,size(ims,1),size(ims,2),[],size(ims,4));
    processed=processed(:,:,2:end-1,:);
    
    [steepest_slope,steepest_time]=max((processed(:,:,:,3:end-2)-processed(:,:,:,1:end-4))/2/deltat,[],4);  %ignore last two "unclosed" time points
   