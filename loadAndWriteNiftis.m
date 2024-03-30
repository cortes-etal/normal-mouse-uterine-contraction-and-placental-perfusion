

dirs = dir(fullfile('**','*_v2.mat')) 

%%
for jj = 1:numel(dirs)
    files = dir(fullfile(dirs(jj).folder,'*.mat'));
    fname=(fullfile(files(1).folder,files(1).name));
    load(fname);
    fprts=regexp(fname,'\','split');
    niftiwrite(enhancement,[[fprts{4} fprts{5}],'_.nii']);
end