function imageToMP4(DImage,label,map,range,cbar)

vidObj = VideoWriter(label,'MPEG-4'); %initializing the vdeio object
vidObj.FrameRate = 1/.25;
open(vidObj);
figure;
for ii = 1:size(DImage,3)
    ii
    imagesc(squeeze(DImage(:,:,ii,:))); %colormap(jet(256));colorbar;
    axis('image');
    axis('off');
%     colormap(map);
    if cbar
        colorbar;
    end
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
end

close(vidObj);
close all;



end