function plotting_combinedPSDs(per,welch,thom,ar,fper,fwel,fthom,far,outDir,saveDir,num_region)


figure(112);
set(gcf,'Position',[10 10 1750 400]);

subplot(1,3,1);
% hold on;
h=plot(fper{1},10*log10(abs(per{1})),'r');
set(h,'LineWidth',2);
% h=plot(fwel{1},10*log10(abs(welch{1})),'k');
% set(h,'LineWidth',2);
% h=plot(fthom{1},10*log10(abs(thom{1})),'b');
% set(h,'LineWidth',2);
% h=plot(far{1},10*log10(abs(ar{1})),'g');
% set(h,'LineWidth',2);
% hold off;
ylim([-50 15]);
xlabel('Frequency [Hz]');
ylabel('Power (dB)');
legend('Periodogram','Welch','Thompson','Burg N =4');
title("PSD Estimates High Perf");
set(gca,'FontWeight','bold','LineWidth',2); %% formatting some plot paramterss

subplot(1,3,2)

% hold on;
h=plot(fper{2},10*log10(abs(per{2})),'r');
set(h,'LineWidth',2);
% h=plot(fwel{2},10*log10(abs(welch{2})),'k');
% set(h,'LineWidth',2);
% h=plot(fthom{2},10*log10(abs(thom{2})),'b');
% set(h,'LineWidth',2);
% h=plot(far{2},10*log10(abs(ar{2})),'g');
% set(h,'LineWidth',2);
% hold off;
ylim([-50 15]);
xlabel('Frequency [Hz]');
ylabel('Power (dB)');
legend('Periodogram','Welch','Thompson','Burg N =4');
title("PSD Estimates High Perf");
set(gca,'FontWeight','bold','LineWidth',2); %% formatting some plot paramterss


subplot(1,3,3)

% hold on;
h=plot(fper{3},10*log10(abs(per{3})),'r');
set(h,'LineWidth',2);
% h=plot(fwel{3},10*log10(abs(welch{3})),'k');
% set(h,'LineWidth',2);
% h=plot(fthom{3},10*log10(abs(thom{3})),'b');
% set(h,'LineWidth',2);
% h=plot(far{3},10*log10(abs(ar{3})),'g');
% set(h,'LineWidth',2);
% hold off;

xlim([0 .03]);
ylim([-50 15]);
xlabel('Frequency [Hz]');
ylabel('Power (dB)');
legend('Periodogram','Welch','Thompson','Burg N =4');
title("PSD Estimates High Perf");
set(gca,'FontWeight','bold','LineWidth',2); %% formatting some plot paramterss


saveas(gcf,fullfile(outDir,saveDir,['psd',num2str(num_region),'_.png']),'png')

end