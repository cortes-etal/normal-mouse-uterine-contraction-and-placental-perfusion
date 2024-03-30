

%%
dft = fft(highCurve,10*length(highCurve));
fs = mean(diff(ACQ_abs_time));

dft=dft(1:length(dft)/2+1);

freq = 0:fs/(10*length(highCurve)):fs/2;


dft(1)=0;
figure;
plot(freq,dft);

%%
 fc = 1.5;
 Wn  = fc/(fs/2); %normalzied cutoff frequnecy 
% 
% 
 [b,a] = butter(6,Wn);
 
xFilt = filtfilt(b,a,highCurve);

figure;
plot(highCurve);
hold on;
plot(xFilt);

%%

dft = fft(xFilt,10*length(xFilt));
fs = mean(diff(ACQ_abs_time));

dft=dft(1:length(dft)/2+1);

freq = 0:fs/(10*length(highCurve)):fs/2;


dft(1)=0;
figure;
plot(freq,dft);
