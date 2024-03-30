function plotting_stftDB(stftDBhigh,stfDBlow,stftDbwhole,t,f,outDir,saveDir,region_num)



figure(111);
set(gcf,'Position',[10 10 1750 400]);
subplot(1,3,1);
imagesc(t,f,stftDBhigh, [-50 10]);colormap(jet);colorbar;set(gca,'YDir','normal');
title('Spectogram of High Perfusion Time Trace');
ylabel('Frequency [Hz]');
xlabel('Time [s]');
set(gca,'FontWeight','bold');

subplot(1,3,2);

imagesc(t,f,stfDBlow, [-50 10]);colormap(jet);colorbar;set(gca,'YDir','normal');
title('Spectogram of Low Perfusion Time Trace');
ylabel('Frequency [Hz]');
xlabel('Time [s]');
set(gca,'FontWeight','bold');

subplot(1,3,3);

imagesc(t,f,stftDbwhole, [-50 10]);colormap(jet);colorbar;set(gca,'YDir','normal');
title('Spectogram of Whole Placenta Perfusion Time Trace');
ylabel('Frequency [Hz]');
xlabel('Time [s]');
set(gca,'FontWeight','bold');


saveas(gcf,fullfile(outDir,saveDir,['spectogram',num2str(region_num),'_.png']),'png')


end