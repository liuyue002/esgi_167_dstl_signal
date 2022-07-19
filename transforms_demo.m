omegaf = @(t) 3e2*(t<0.1) + 1e2*(t>=0.1).*(t<0.2) + 6e2*(t>=0.2);
%omegaf = @(t) 1e2;
%omegaf = @(t) 3e2 + 5e2*t;
xf = @(t) cos(2*pi*omegaf(t).*t);

N=3000; % num data pt
T=0.3; %total time
Fs=N/T; % sample frequency
ts=linspace(0,T,N);
xs=xf(ts);
xs=xs+normrnd(0,0.2,size(xs)); %noise
figure;
plot(ts,xs);
title('signal');

%%
yfft=fft(xs);
P2 = abs(yfft/N);
P1 = P2(1:N/2+1)*2;
f = Fs*(0:(N/2))/N;
figure;
plot(f,P1);
xlabel('frequency (Hz)');
title('Fourier transform');

%%
ywht=fwht(xs);
figure;
plot(ywht);
title('Walsh-Hadamard');

%%
figure;
stft(xs,Fs,'Window',hann(200,'periodic'));
title('Short-time FT, Hann 200');

figure;
stft(xs,Fs,'Window',flattopwin(200));
title('Short-time FT, flat top 200');

figure;
stft(xs,Fs,'Window',gausswin(200));
title('Short-time FT, gauss 200');

figure;
stft(xs,Fs,'Window',rectwin(200));
title('Short-time FT, rectangle 200');

%%
figure;
cwt(xs,Fs);
title('Continuous wavelet, morse wavelet');

%%

figure;
cqt(xs,'SamplingFrequency',Fs);
title('Gabor const-Q, default');

%%
figure;
wvd(xs,Fs);
title('Wigner ville, default');