omegaf1 = @(t) 3e2*(t<0.1) + 1e2*(t>=0.1).*(t<0.2) + 6e2*(t>=0.2);
omegaf2 = @(t) 1e2;
omegaf3 = @(t) 1e2 + 2e3*t;
omegaf4 = @(t) 1e2 + 1.5e4*(t.^2);
xf = @(t) 1*cos(2*pi*omegaf4(t).*t) + 2*cos(2*pi*omegaf1(t).*t);

N=3000; % num data pt
T=0.3; %total time
Fs=N/T; % sample frequency
ts=linspace(0,T,N);
xs=xf(ts);
xs=xs+normrnd(0,1,size(xs)); %noise
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
ylim([-2,2])

figure;
stft(xs,Fs,'Window',flattopwin(200));
title('Short-time FT, flat top 200');
ylim([-2,2])

figure;
stft(xs,Fs,'Window',gausswin(200));
title('Short-time FT, gauss 200');
ylim([-2,2])

figure;
stft(xs,Fs,'Window',rectwin(200));
title('Short-time FT, rectangle 200');
ylim([-2,2])

%%
figure;
[cfs,frq] = cwt(xs,Fs);
surface(ts,frq,abs(cfs))
axis tight
shading flat
xlabel("Time (s)")
ylabel("Frequency (Hz)")
%set(gca,"yscale","log")
title('Continuous wavelet, morse wavelet');
ylim([0,2000]);
colorbar;

% figure;
% [a,d] = haart(xs,10,'noninteger');
% plot(d{1});
% title('Haar transform');

%%

figure;
cqt(xs,'SamplingFrequency',Fs);
title('Gabor const-Q, default');
ylim([0,2])

%%
figure;
wvd(xs,Fs);
title('Wigner ville, default');
ylim([0,2]);
caxis([0,500]);
%%
figure;
wsst(xs,Fs);
title('Wavelet synchrosqueezed, default');
ylim([0,2])
colorbar;