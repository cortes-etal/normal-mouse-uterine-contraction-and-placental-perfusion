% Deivn Cortes
%creating histogram figure for DCE paper/ 


clear all;clc;close all;

%%
etoh14_dat = load('F:\DCE Analysis\Wu PCA4\PerfusionMaps\578 EtOHE14.5 (08-22-21)\perfMap_03-Dec-2021_.mat');
etoh17_dat = load('F:\DCE Analysis\Wu PCA4\PerfusionMaps\578 EtOHE17.5 (08-25-21)\perfMap_01-Dec-2021_.mat');
cont14_dat = load('F:\DCE Analysis\Wu PCA4\PerfusionMaps\599 Cont-PSFE14.5 (08-23-21)\perfMap_01-Dec-2021_.mat');
cont17_dat = load('F:\DCE Analysis\Wu PCA4\PerfusionMaps\599 Cont-PSFE17.5 (08-26-21)\perfMap_01-Dec-2021_.mat');

%% =loading slope maps

map14e = etoh14_dat.slopeMap2 .* etoh14_dat.*
map17e = etoh17_dat.slopeMap2;
map14c = cont14_dat.slopeMap2;
map17c = cont17_dat.slopeMap2;
imList = {map14e,map14c,map17e,map17c};
figure;
montage(imList,'Size',[1,4])