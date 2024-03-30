% purpose of the code is to automate the folder selection,organization and
% saving process of the T1 Mapping/ r1Mapping code of the placental dce

%%CODE EXPECTS A CERTAIN FILE/FOLDER TREE STRUCTURE DEFINED HERE AND IN
%%HELP DOCUMENTATION %%%%%%
clear all;close all;clc;
numAnimals =1% dir(fullfile(pwd,'\*Brain'));

mainpath='C:\Users\ddblabwu\Documents\MATLAB\MATLAB\code\code\Anthony\Anthony'; %PUT PATH TO THIS FOLDER HERE
addpath(fullfile(mainpath,'bruker_io'));

for control = 1:numAnimals
    
    folders = uigetfile_n_dir;
    
    
    %
    %     for ii = 1:numel(folders)
    %         temp = circshift(folders, [0 -ii]);
    %         orderedDirs(ii).cell=temp;
    %     end
    
    
    k=1;
    for mm = 1:numel(folders)
        mm
        %  k=1;
        r3=@(x)reshape(x,1,1,[]);
        %         orderedFolders = orderedDirs(ii).cell;
        %         for j=1:numel(orderedFolders)
        foldername=folders{mm};
        acqp{mm}=read_bruker_acqp(fullfile(foldername,'acqp'));
        try
            reco=read_bruker_acqp(fullfile(foldername,'pdata/1/reco'));
        catch %try again in case Box just isn't responsive the first time
            reco=read_bruker_acqp(fullfile(foldername,'pdata/1/reco'));
        end
        im=load_bruker_2dseq(reco.RECO_size(1),fullfile(foldername,'pdata/1/2dseq'),'int32');
        im=reshape(im,reco.RECO_size(2),reco.RECO_size(1),[]);
        ims(:,:,1:size(im,3),k)=bsxfun(@rdivide,im,r3(reco.RECO_map_slope)); %store the image, correcting for the storage scaling factor
        alphas(1,k)=acqp{mm}.ACQ_flip_angle; %also store the flip angle
        try
            TRs(1,k)=acqp{mm}.ACQ_repetition_time*1e-3; %also store the repetition time
        catch
            acqR = acqp{mm}.ACQ_repetition_time(1)*1e-3;
            TRs(1,k) = acqR;
        end
        k=k+1;
        %  end
        %         allIms(mm).im = ims;
    end
    
    alphas=sort(alphas,'descend');
    Nims = numel(folders);
    %%
    prts = regexp(folders{1},'\','split');
    mainDir=strcat(prts{4},prts{5});
    
    
    %%
    % flips = [75  30 20 10];
    %   for ii =1:numel(allIms)
    %         ims =allIms(end).im;
    im=ims;
    %     if ii < 4
    %         for jj = 1:4 
    %            ims(:,:,:,jj) = im;
    %         end
    %     end
    % Nonlinear fitting
    vec=@(x)x(:);
    e=@(T1)exp(-TRs/T1);
    nonlin=@(e)(1-e)./(1-e.*cos(alphas*pi/180)).*sin(alphas*pi/180) %The combined nonlinear factor
    imvals=reshape(ims,[],Nims);
    mapvals=zeros(size(imvals,1),2);
    Avarpro=@(curve,nonlin)curve*nonlin.'/norm(nonlin)^2; %The linear factor A as a function of T1 (variable projection)
    
    %initial guess
    T1map=(imvals(:,1)*sin(alphas(end)*pi/180)-imvals(:,end)*sin(alphas(1)*pi/180)) ...
        ./(imvals(:,1)*cos(alphas(1)*pi/180)*sin(alphas(end)*pi/180)-imvals(:,end)*cos(alphas(end)*pi/180)*sin(alphas(1)*pi/180));
    T1map=reshape(real(-TRs(1)./log(T1map)),size(ims,1),size(ims,2),[]);
    T1map(T1map>10)=10;
    T1map(T1map<.05)=0.05;
    figure,imagesc(T1map(:,:,ceil(end/2+1)),[0 prctile(T1map(:),95)]),colormap(jet(256)),colorbar
    drawnow
    
    parfor j=1:size(imvals,1)
        cost=@(nonlin)imvals(j,:)-Avarpro(imvals(j,:),nonlin)*nonlin;
        temp=lsqnonlin(@(T1)cost(nonlin(e(T1))),1,.05,10);
        mapvals(j,:)=[temp, Avarpro(imvals(j,:),nonlin(e(temp)))];
    end
    T1map=reshape(mapvals(:,1),size(ims,1),size(ims,2),[]); %place fitted values back in map
    Amap=reshape(mapvals(:,2),size(ims,1),size(ims,2),[]); %place fitted values back in amplitude map
    R1map=1./T1map;
    
    %         allT1(ii).t1=T1map;
    %         allA(ii).A=Amap;
    %         allR1(ii).r1=R1map;
    %
    
    %% displying and saving
    
    mkdir(fullfile('R1Maps_post',mainDir,'R1Maps'));
    imageToMP4(R1map,fullfile('R1Maps_post',mainDir,'R1Maps',strcat('\R1Map','.mp4')));
    for jj = 1:size(R1map,3)
        imagesc(R1map(:,:,jj),[0 prctile(R1map(:),95)]);
        title(['R1 [Hz]:  ', ', Slice ', num2str(jj)],'interpreter','none')
        try
            colormap(turbo(256))
        catch
            colormap(jet(256))
        end
        colorbar
        saveas(gcf,fullfile('R1Maps_post',mainDir,strcat('R1Maps\','R1map_', '_slice_',num2str(jj))),'png');
        saveas(gcf,fullfile('R1Maps_post',mainDir,strcat('R1Maps\','R1map_', '_slice_',num2str(jj))),'fig');
        % end
        
    end
    
    save(fullfile(fullfile('R1Maps_post',mainDir,'R1Maps\'),'T1maps.mat'),'ims','Amap','T1map','R1map');
end

