%%

% M = 2; % Modulation order
% freqsep = 8; % Frequency separation (Hz)
% nsamp = 40; % Number of samples per symbol
% Fs = 40; % Sample rate (Hz)
% msg1 = randi([0 M-1],100,1);
% xs1 = fskmod(msg1,M,freqsep,nsamp,Fs*3);
% msg2 = randi([0 M-1],100,1);
% xs2 = fskmod(msg2,M,freqsep,nsamp,Fs);
% xs=[xs1;xs2];

msgN=100;
msg = randi([0 1],msgN,1);
Fs = 5*msgN;
ts=0:1/Fs:msgN;
N=length(ts);
msgsig=zeros(1,N);
i=1; 
for j=1:N
    if ts(j)<=i 
        msgsig(j)=msg(i); 
    else 
        i=i+1; 
    end 
end 
x1=sin(2*pi*20*ts);
x2=sin(2*pi*25*ts);
xs1=(msgsig==0).*x1 + (msgsig==1).*x2;
x1=sin(2*pi*40*ts);
x2=sin(2*pi*45*ts);
xs2=(msgsig==0).*x1 + (msgsig==1).*x2;
xs=[xs1,xs2];
xs=xs+normrnd(0,0.05,size(xs)); %noise
%%
%xs=imag(xs');
N=length(xs);
ts=0:1/Fs:(N-1)/Fs;
figure('Position',[100 300 1200 400],'color','w');
plot(ts,xs);
title('signal');
xlabel('t (s)');
xlim([0,ts(end)]);
ylim([-1.1,1.1]);

%%
% for wid=[50,100,200,400]
%     figure;
%     stft(xs,Fs,'Window',hann(wid,'periodic'),'FrequencyRange','onesided');
%     title(['Short-time FT, Hann ',num2str(wid)]);
% end

wid=200;
figure('Position',[100 300 1200 400],'color','w');
[s,f,t]=fsst(xs,Fs,hann(wid,'periodic'),'yaxis');
mesh(t,f,abs(s))
axis tight
view(2)
colorbar
ylabel('Frequency');
xlabel('Time');
title(['FSST, Hann ',num2str(wid)]);
%%
mean_freq=zeros(1,N);
s2=abs(s);
for i=1:N
    mean_freq(i)=mean(s2(:,i).*f);
end
figure;
plot(ts,mean_freq);
title('Mean frequency');
xlabel('t');
ylabel('f');

%% filtered derivative 1
l=30;
%z=[zeros(1,l),xs]-[xs,zeros(1,l)];
z=conv(mean_freq,[-1,zeros(1,l-1),1]);
z=z(1+l:end);
figure;
plot(ts,z);

%% filtered derivative 2
l=200;
%z=[zeros(1,l),xs]-[xs,zeros(1,l)];
z=conv(mean_freq,[-1*ones(1,l),ones(1,l)]);
z=z(1+l*2-1:end)/l^2;
z(end-2*l:end)=0;
figure;
plot(ts,z);

%% Shiryaev

z=cumsum(mean_freq);
figure;
plot(ts,z);