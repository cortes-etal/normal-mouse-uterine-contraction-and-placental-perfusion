function plotting_perfcurves(high,low,whole)


figure(109);
subplot(1,2,1)
plot(ACQ_abs_time./3600,lowCurve./rms(lowCurve));
hold on;
plot(ACQ_abs_time./3600,highCurve./rms(highCurve));
plot(ACQ_abs_time./3600,wholeCurve./rms(wholeCurve));
hold off;
xlabel('Minutes [m]');
title('Intensity Curve');
legend('Low Perfusion Chamber','High Perfusion Chamber','Whole Placenta','Location','southwest');
drawnow;
subplot(1,2,2)
plot(ACQ_abs_time(1:end-1)./3600,diff(lowCurve./rms(lowCurve)));
hold on;
plot(ACQ_abs_time(1:end-1)./3600,diff(highCurve./rms(highCurve)));
plot(ACQ_abs_time(1:end-1)./3600,diff(wholeCurve./rms(wholeCurve)));
hold off;
xlabel('Minutes [m]');
title(animal);
legend('Low Perfusion Chamber','High Perfusion Chamber','Whole Placenta','Location','southwest');

saveas(gcf,fullfile(outDir,saveDir,['intensityCurves_placenta_',num2str(l),'_.png']),'png')

drawnow;
end