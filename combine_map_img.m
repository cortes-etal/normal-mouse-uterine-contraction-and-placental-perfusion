function [RGB, slopeLow, slopeHigh] = combineRGB_nD(grayIm,regions,fname,animal,ageFlag)
fname
[iG,jG,kG] = size(regions);
RGB = zeros([iG jG 3 kG]); %  3 to represent the 3 color channels
fRGB =RGB;
v = VideoWriter([fname,'_111_.mp4'],'MPEG-4');
v.FrameRate=2;
open(v);
figure(1);
for kk = 1:kG
    greydouble =(double2rgb(squeeze(grayIm(:,:,kk)),gray(256)));
    mapDouble = double2rgb(regions(:,:,kk),colormap(jet(256)), [min(regions(:)) max(regions(:))]);
    
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
        currFrame = getframe(gcf);
        writeVideo(v,currFrame);
    
    fRGB(:,:,:,kk) = finalIm;
    
end
close(v);

% close(v);
end