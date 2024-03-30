

% purpose: create kinetic movies for each animal
%load each scan as a separate image, create movies across the different
%slices to create kinetic movie maps

clear all;close all;clc; % start pretty


%%

files = dir('**/*_v2.mat');
testSig = zeros(50,1);
%%
segs = dir('*.nii.gz');
segNames = extractfield(segs,'name');

outDir = ['kinetic_movies_nifti_exports' ,date] ;
mkdir(outDir)


%%

for ii = 5%1:numel(files)
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
    
    animalE = animalE(animalEdx:animalEdx+1); % the dataset animal embryonic age
    
    segID = contains(segNames, animalE) & contains(segNames,animalID); % logical vecotr containing idx for seg file
    sFile = fullfile(segs(segID).folder,segs(segID).name); % fname of seg file
    sDat = niftiread(sFile); % load the segmentation file % loading seg fle.\
    
    pMask = sDat > 0;
    
    fname = fullfile(files(ii).folder,files(ii).name);
    fprts2 =regexp(fname,'\','split');
    animal=[fprts2{4} fprts2{5}];
    %     mkdir(fullfile('PerfusionMaps',animal));
    load(fname);
    
    
    test_slice = 13;
    
     slice_13_og = zeros(size(ims,1),size(ims,2), size(ims,4));
    
    for jj = 1:size(ims,4)
        slice_13_og(:,:,jj) = ims(:,:,test_slice,jj);
        
    end
    
    t1 = slice_13_og(:,:,1);
    t2 = slice_13_og(:,:,12);
    t3 = slice_13_og(:,:,24);
    t4 = slice_13_og(:,:,36);
    t5 = squeeze(slice_13_og(:,:,48));
    
    sz=size(t5);
    xg = 1:sz(1);
    yg = 1:sz(2);
    F = griddedInterpolant({xg,yg},double(t5),'cubic');
    
    xq = (0:5/6:sz(1))';
    yq = (0:5/6:sz(2))';
    vq = uint8(F({xq,yq}));
    imshow(vq)
    title('Higher Resolution')
    %     counter=13;
    %     while counter <= size(ims,3)
    %         for jj = 1:size(ims,4)
    %             slices(:,:,jj) = ims(:,:,counter,jj);%.*pMask(:,:,13);
    %
    %         end
    %
    %         niftiwrite(slices,fullfile(outDir,saveDir,['slice_',num2str(counter),'.nii']))
    %
    %         imageToMP4(slices,fullfile(outDir,saveDir,['slice_',num2str(counter),'_nomask.mp4']),gray(256),[min(slices(:)) max(slices(:))],0);
    %         clear slices;
    %         counter=counter+1;
    %     end
    
end

%%

%end
