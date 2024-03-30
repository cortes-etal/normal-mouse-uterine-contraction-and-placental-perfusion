function  proc_signal = process_motion_signals(signal,fs)
signal = (signal - mean(signal)) ./ rms(signal);

Wn= [0.01 / (fs/2) 0.02 / (fs/2)];
[b,a] = butter(4,Wn);
proc_signal = filtfilt(b,a,signal);

time = 0:(1/fs):(numel(signal)*(1/fs) - (1/fs));
figure;
plot(time./60,signal);

n_signal =numel(signal);
n_fft = 10* n_signal;


win_cheb_200 = chebwin(n_signal,100).';
win_cheb_50 = chebwin(n_signal,50).';

win_ham = hamming(n_signal).';

win_han = hanning(n_signal).';

fft_cheb200 = fftshift(fft(win_cheb_200 .*signal,n_fft));
fft_cheb50 = fftshift(fft(win_cheb_50 .*signal,n_fft));
fft_win_ham = fftshift(fft(win_ham .*signal,n_fft));
fft_win_han = fftshift(fft(win_han .*signal,n_fft));


fft_cheb200 = fft_cheb200(n_fft/2:end);
fft_cheb50 = fft_cheb50(n_fft/2:end);
fft_win_ham = fft_win_ham(n_fft/2:end);
fft_win_han = fft_win_han(n_fft/2:end);

freq = 0: fs/n_fft : fs/2; % frequency vector 
% 
% figure;
% subplot(1,4,1);
% plot(freq,fft_cheb200);
% subplot(1,4,2);
% plot(freq,fft_cheb50);
% subplot(1,4,3);
% plot(freq,fft_win_ham);
% subplot(1,4,4);
% plot(freq,fft_win_han);










end