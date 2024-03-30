function [RGB,fRGB, slopeLow, slopeHigh] = combineRGB_nD(grayIm,regions,fname,animal,ageFlag)
fname
[iG,jG,kG] = size(regions);
RGB = zeros([iG jG 3 kG]); %  3 to represent the 3 color channels
fRGB =RGB;
 v = VideoWriter([fname,'_111_r1_.mp4'],'MPEG-4');
 v.FrameRate=10;
 open(v); figure(1);
for kk = 1:kG
    greydouble =(double2rgb(squeeze(grayIm(:,:,kk)),gray(256)));
    %     greydouble = double2rgb(squeeze(grayIm(:,:,kk)),gray(256));
    s=regions(:,:,kk);
    mapDouble = double2rgb(regions(:,:,kk),colormap(jet(256)), [0 max(s(:))]);
    
    outR = greydouble(:,:,1);
    outG = greydouble(:,:,2);
    outB = greydouble(:,:,3);
    
    mapR = mapDouble(:,:,1);
    mapG = mapDouble(:,:,2);
    mapB = mapDouble(:,:,3);
    
    outR(squeeze(regions(:,:,kk))> 0) = mapR(squeeze(regions(:,:,kk)) >0);
    outG(squeeze(regions(:,:,kk))> 0) = mapG(squeeze(regions(:,:,kk)) >0);
    outB(squeeze(regions(:,:,kk))> 0) = mapB(squeeze(regions(:,:,kk)) >0);
    
    finalIm = cat(3,outR,outG,outB);
    
    figure(1);
    imshow(finalIm)
    pause(0.25);
         currFrame = getframe(gcf);
         writeVideo(v,currFrame);
    
    fRGB(:,:,:,kk) = finalIm;
    
end
 close(v);


%%
%
figure(3)
subplot(1,2,1);
histogram(log10(regions(regions > 0)));
subplot(1,2,2);
histogram((regions(regions > 0)));



%
%
% title('Whole Placenta Perfusion Histogrm');
% ylabel('Frequency');
% xlabel('Perfusion mL/min/100mL');
% subplot(1,2,2);
% histogram(log10(regions(regions >0)))
%
% title('Whole Placenta Perfusion Histogram - Log Transform');
% ylabel('Frequency');
% xlabel('Log10(Perfusion mL/min/100mL)');
% sgtitle(animal);

% saveas(gcf,fname,'png');



% pause(0.24);

% set(gca,'yscale','log');
regionDat = regions(regions> 0);

if strcmp(ageFlag,'17')
    perfThresh = prctile(regionDat,55)
else
    perfThresh = prctile(regionDat,41)
end

% pause;


lowPerfMask = regions < perfThresh;
highPerfMask = regions > perfThresh;



slopeLow = regions.*lowPerfMask;
slopeHigh = regions.*highPerfMask;

% figure(111);
% bar(1:2,[mean(nonzeros(slopeLow)) mean(nonzeros(slopeHigh))]);
% % pause;


 v = VideoWriter([fname,'r1_.mp4'],'MPEG-4');
 v.FrameRate=10;
 open(v);
for kk = 1:kG
    mapLow = double2rgb(slopeLow(:,:,kk),cool(256),[0 1.2]);
    mapHigh = double2rgb(slopeHigh(:,:,kk),hot(256),[1.3 3.5]);
    greydouble =(double2rgb(squeeze(grayIm(:,:,kk)),gray(256)));
    %     greydouble = double2rgb(squeeze(grayIm(:,:,kk)),gray(256));
    
    lowR =mapLow(:,:,1);
    lowG = mapLow(:,:,2);
    lowB = mapLow(:,:,3);
    
    highR = mapHigh(:,:,1);
    highG = mapHigh(:,:,2);
    highB = mapHigh(:,:,3);
    
    noutR = greydouble(:,:,1);
    noutG = greydouble(:,:,2);
    noutB = greydouble(:,:,3);
    
    noutR(squeeze(slopeLow(:,:,kk)) > 0) = lowR(squeeze(slopeLow(:,:,kk)) > 0);
    noutR(squeeze(slopeHigh(:,:,kk)) > 0) = highR(squeeze(slopeHigh(:,:,kk)) >0);
    
    noutG(squeeze(slopeLow(:,:,kk)) >0) = lowG(squeeze(slopeLow(:,:,kk)) > 0);
    noutG(squeeze(slopeHigh(:,:,kk)) > 0) = highG(squeeze(slopeHigh(:,:,kk)) > 0);
    
    noutB(squeeze(slopeLow(:,:,kk)) > 0) = lowB(squeeze(slopeLow(:,:,kk)) >0);
    noutB(squeeze(slopeHigh(:,:,kk)) > 0) = highB(squeeze(slopeHigh(:,:,kk)) >0);
    
    combinedPerfMaps = cat(3,noutR,noutG,noutB);
    
    %
    
         figure(1);
         imshow(combinedPerfMaps);
         pause(0.25);
         currFrame = getframe(gcf);
         writeVideo(v,currFrame);
    %    %
        RGB(:,:,:,kk) = combinedPerfMaps;
    %     pause;
    
    
end
 close(v);
end