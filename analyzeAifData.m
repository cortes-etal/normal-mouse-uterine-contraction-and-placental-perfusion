
% start pretty
clear all;clc;close all;

%%

files = dir('**/*_v2.mat');
testSig = zeros(50,1);
figure;
for ii = 1:numel(files)
    fname = fullfile(files(ii).folder,files(ii).name);
    fprts =regexp(fname,'\','split');
    load(fname);
    [maxt,idx] = max(AIF);
    testSig(idx) = AIFmax;
    pseudoMax(:,ii) = testSig;
    aifs(:,ii) = AIF;
    labels{ii} = strcat(fprts{4}, ' ', fprts{5});
    
    hold on;
    if contains(labels(ii),'EtOH')
        plot(1:50,AIF,'k');
        plot((idx),maxt,'k+','MarkerSize',8);
    else
        plot(1:50,AIF,'r');
        plot((idx),maxt,'k+','MarkerSize',8);
        
    end
end
title('Arteral Input Functions');
ylabel('Perfusion');

hold off;

%%
figure;
plot(aifs);
legend(labels);
