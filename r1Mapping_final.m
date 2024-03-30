clear all;
close all hidden;

mainpath='C:\Users\ddblabwu\Documents\MATLAB\MATLAB\code\code\Anthony\Anthony'; %PUT PATH TO THIS FOLDER HERE
addpath(fullfile(mainpath,'bruker_io'));

%% Find IntraGate datasets
foldername=uigetdir; %Choose folder which has numbered subfolders 1, 2, 3, ...
not_done=true;
j=10; %STARTING FOLDER
while not_done %go through all numbered subfolders
    if exist(fullfile(foldername,sprintf('%0.d',j)),'dir')
        acqp{j}=read_bruker_acqp(fullfile(foldername,sprintf('%0.d',j),'acqp'));
        isIGF=strfind(acqp{j}.ACQ_method,'IntraGateFLASH');
        if isempty(isIGF)
            isIGF=true;
        end
        isIntraGate(j)=logical(isIGF);
        j=j+1;
    else
        break;
    end
end
%%
Nims=sum(isIntraGate)
% if Nims ~= 5
%     disp('Unexpected number of images found!')
%     keyboard;
% end

%% Load IntraGate datasets
k=1;
r3=@(x)reshape(x,1,1,[]);
counter=1;
for j=find(isIntraGate)
    try
        reco=read_bruker_acqp(fullfile(foldername,sprintf('%0.d',j),'pdata/1/reco'));
    catch %try again in case Box just isn't responsive the first time
        reco=read_bruker_acqp(fullfile(foldername,sprintf('%0.d',j),'pdata/1/reco'));
    end
    im=load_bruker_2dseq(reco.RECO_size(1),fullfile(foldername,sprintf('%0.d',j),'pdata/1/2dseq'),'int32');
    im=reshape(im,reco.RECO_size(2),reco.RECO_size(1),[]);
    ims(:,:,1:size(im,3),k)=bsxfun(@rdivide,im,r3(reco.RECO_map_slope)); %store the image, correcting for the storage scaling factor
    alphas(1,k)=acqp{j}.ACQ_flip_angle; %also store the flip angle
    try
        TRs(1,k)=acqp{j}.ACQ_repetition_time*1e-3; %also store the repetition time
    catch
        acqR = acqp{j}.ACQ_repetition_time(1)*1e-3;
        TRs(1,k) = acqR;
    end
    k=k+1;
    allIms(counter).ims=ims;
    counter=counter+1;
end

%% Fit T1s
%E=exp(-TR/T1)
%y=A*(1-E)/(1-E*cos(alpha))*sin(alpha)
allIms(1:5)=[];
%%
Nims =4;
for ii =4%1:numel(allIms)
    ims =allIms(end).ims;
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
    figure,imagesc(T1map(:,:,ceil(end/2+1)),[0 3]),colormap(jet(256)),colorbar
    drawnow
    
    parfor j=1:size(imvals,1)
        cost=@(nonlin)imvals(j,:)-Avarpro(imvals(j,:),nonlin)*nonlin;
        temp=lsqnonlin(@(T1)cost(nonlin(e(T1))),1,.05,10);
        mapvals(j,:)=[temp, Avarpro(imvals(j,:),nonlin(e(temp)))];
    end
    T1map=reshape(mapvals(:,1),size(ims,1),size(ims,2),[]); %place fitted values back in map
    Amap=reshape(mapvals(:,2),size(ims,1),size(ims,2),[]); %place fitted values back in amplitude map
    R1map=1./T1map;
    
    allT1(ii).t1=T1map;
    allA(ii).A=Amap;
    allR1(ii).r1=R1map;
end
%% Display
mkdir T1Maps
flips = [10 20 30 55];
%figure,imagesc(Amap(:,:,ceil(end/2+1)),[0 prctile(Amap(:),99)]),colormap(gray(256)),colorbar
for ii =1:numel(allT1)
    T1map = allT1(ii).t1;
    imageToMP4(T1map,strcat('T1Maps\fa_',num2str(flips(ii)),'.mp4'));
    for jj = 1:size(T1map,3)
        imagesc(T1map(:,:,jj),[0 5]);
        title(['T1 [s]:  ','FA ' , num2str(flips(ii)), ', Slice ', num2str(jj)],'interpreter','none')
        try
            colormap(turbo(256))
        catch
            colormap(jet(256))
        end
        colorbar
        saveas(gcf,strcat('T1Maps\','T1map_fa_',num2str(flips(ii)), '_slice_',num2str(jj)),'png');
        saveas(gcf,strcat('T1Maps\','T1map_fa_',num2str(flips(ii)), '_slice_',num2str(jj)),'fig');       
    end

end
%% firsdt R1 mapcklc

mkdir R1Maps
flips = [75 55 30 20 10];
figure;
for ii = 1:numel(allR1)
    
    R1map = allR1(ii).r1;
    for jj = 1:size(R1map,3)
        imagesc(R1map(:,:,jj),[0 5]);
        title(['R1 [Hz]:  ','FA ' , num2str(flips(ii)), ', Slice ', num2str(jj)],'interpreter','none')
        try
            colormap(turbo(256))
        catch
            colormap(jet(256))
        end
        colorbar
        saveas(gcf,strcat('R1Maps\','R1map_fa_',num2str(flips(ii)), '_slice_',num2str(jj)),'png');
        saveas(gcf,strcat('R1Maps\','R1map_fa_',num2str(flips(ii)), '_slice_',num2str(jj)),'fig');
        
    end
end

%%
mkdir R1MapComp
figure;
for ii =1:numel(allR1)
    R1map = allR1(ii).r1;
    Amap = allA(ii).A;
    [N,edges]=histcounts(Amap(Amap>median(Amap(:))));
    [~,argmax]=max(N);
    Anorm=(edges(argmax)+edges(argmax+1))/2;
    R1map_comp = R1map .* Amap / Anorm;
    imageToMP4(R1map_comp,strcat('movies\R1Maps\R1map_fa_',num2str(flips(ii)),'.mp4'));
    for jj = 1:size(R1map_comp,3)
        imagesc(R1map_comp(:,:,jj),[0 5])
        title(['R1 [Hz]:  ','FA: ',num2str(flips(ii)),', Slice ', num2str(jj)],'interpreter','none')
        try
            colormap(turbo(256))
        catch
            colormap(jet(256))
        end
        colorbar
        saveas(gcf,strcat('R1MapComp\','R1mapComp_fa_',num2str(flips(ii)), '_slice_',num2str(jj)),'png');
        saveas(gcf,strcat('R1MapComp\','R1mapComp_fa_',num2str(flips(ii)), '_slice_',num2str(jj)),'fig');
        
    end
end
%%
figure,plot(edges(1:end-1),N)

%% Save
save(fullfile(foldername,'T1maps.mat'),'ims','Amap','T1map','R1map','R1map_comp');
