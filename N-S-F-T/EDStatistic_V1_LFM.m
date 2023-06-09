%--------------------------------------------------------------------------
%   Author      : AKRARAI MOHAMED    
%   School      : National Institut Of Posts and Telecommunications
%   Description : In This file we calculate the Probablity of detection for            
%                 the Energy Detection method, for an LFM pulse for Different
%                 Numbers of Samples[100 1000].  
%--------------------------------------------------------------------------

clc;
clear all

Pd = zeros([4 10]);

k = 0;
for SNR = -50:10:-20
    
k = k + 1;
i = 0;
V = zeros([10 1000]);
    
 for N = 100:100:1000
    i = i + 1;
    fs = 100e3;
    TpulseWidth = N * 1e-5; % pulse width variates from 1 ms to 10 ms
    sig = phased.LinearFMWaveform('SampleRate',fs,'SweepBandwidth',4e3,...
    'PulseWidth',TpulseWidth,'PRF',25);

    Xp = sig();
    
    for j =1:1:1000
        
    Y = awgn(Xp,SNR,'measured'); % the signal power -28.3 dBW
    Noise = Y - Xp;
     %decision statistic Tmf can be computed from
     %the squared magnitude of the FFT averaged over N samples
    La  = abs(fft(Noise).^2);
    LambdAmf = sum(La)/N;
    Sta = abs(fft(Y).^2);
    Tmf = sum(Sta)/N;
    
    if  (Tmf > LambdAmf)    
        V( i , j ) = 1;
    end
    end 
    Pd( k , i ) = (sum(V(i,:))/1000);
 end
end

Nsamples = 1:1:10;
figure;
plot(Nsamples*100,Pd(1,:),Nsamples*100,Pd(2,:),Nsamples*100,Pd(3,:),Nsamples*100,Pd(4,:))
title('Pd(N) for Energy Detection Method, LFM signal')
xlabel('N samples')
ylabel('Probability of Detection')
grid
hold off

legend('-40 dB','-30 dB','-20 dB','-10 dB')