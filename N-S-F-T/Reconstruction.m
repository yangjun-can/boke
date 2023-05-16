 %% 压缩感知重构
clear
K=2;      %  稀疏度  
N=256;    %  信号长度
f1=50.05; %  信号频率1
fs=800;   %  采样频率
ts=1/fs;  %  采样间隔
Ts=1:N;   %  采样序列
t= (Ts-1)*ts; %  真实时间序列
f=zeros(N,1);
for i=1:N
    f(i)=Nonuniform_sampling_point(i-1,N);
end
x=0.2*sin(2*pi*f1*Ts*ts); %  完整信号
tfreq=ndftld(x, N, f);    %  原始信号频域
% figure;plot(f,abs(tfreq));

%% 白噪声
SNR=2;
NOISE=randn(size(x));
NOISE=NOISE-mean(NOISE);
signal_power = 1/length(x)*sum(x.*x);
noise_variance = signal_power / ( 10^(SNR/10) );
NOISE=sqrt(noise_variance)/std(NOISE)*NOISE;
x_noise=x+NOISE;
tfreq_noise=ndftld(x_noise, N, f);    %  原始信号频域
%  figure;plot(f,abs(tfreq_noise));
%  figure;plot(abs(fft(NOISE)));
%%  FDM 
[xt_recov_IMFsLowToHigh, xt_recov_IMFsHighToLow, EnergyLekage] = FMD_Low2High_High2LowSacnning(x_noise,fs,t);
figure
subplot(421);plot(xt_recov_IMFsHighToLow(:,1)); xlabel('Signal in the time domain');ylabel('amplitude');title('y1');
subplot(422);plot(xt_recov_IMFsHighToLow(:,2)); xlabel('Signal in the time domain');ylabel('amplitude');title('y2');
subplot(423);plot(xt_recov_IMFsHighToLow(:,3)); xlabel('Signal in the time domain');ylabel('amplitude');title('y3');
subplot(424);plot(xt_recov_IMFsHighToLow(:,4)); xlabel('Signal in the time domain');ylabel('amplitude');title('y4');
subplot(425);plot(xt_recov_IMFsHighToLow(:,5)); xlabel('Signal in the time domain');ylabel('amplitude');title('y5');
subplot(426);plot(xt_recov_IMFsHighToLow(:,6)); xlabel('Signal in the time domain');ylabel('amplitude');title('y6');
subplot(427);plot(xt_recov_IMFsHighToLow(:,7)); xlabel('Signal in the time domain');ylabel('amplitude');title('y7');
subplot(428);plot(xt_recov_IMFsHighToLow(:,8)); xlabel('Signal in the time domain');ylabel('amplitude');title('y8');


% %% direct
%  [hat_x_direct]=csreconstruction(x_noise,N,K); 

% %% nusft
% [out_Large,BB_loc,BB_est,b_loc,b_est]=nusfft(x_noise,K,f); %  估计信号频域
% xnusft=indft1d(out_Large, N, f);  %  估计信号
% [hat_x_nusft]=csreconstruction(xnusft,N,K); 
% 
% %% 18 Low_rank
% out_Large2=chebfun.nufft(x_noise.',f/N,2).';
% x_lowrank=indft1d(out_Large2, N, f);  %  估计信号
% [hat_x_lowrank]=csreconstruction(x_lowrank,N,K); 
% 
% %%  13 Min_Max
% temp=2*pi*f/N;
% st = nufft_init(temp, N, 6, 2*N, 'minmax:kb');
% Min_Max_freq = nufft(x_noise.',st).';
% x_Min_Max=indft1d(Min_Max_freq, N, f);  %  估计信号
% [hat_x_Min_Max]=csreconstruction(x_Min_Max,N,K); 
% 
% %% 14 Fast_Grid
% FGG_NUFFT  = FGG_1d_type2(x_noise,f*2*pi/N);
% FGG_NUFFT = fftshift(FGG_NUFFT).';
% x_FGG=indft1d(FGG_NUFFT, N, f);  %  估计信号
% [hat_x_FGG]=csreconstruction(x_FGG,N,K); 
% 
% % figure;plot(f,abs(out_Large2));
% % figure;
% % plot(hat_x,'.-')                                 %  重建信号
% % hold on;plot(x,'r');                                      %  原始信号
% % legend('Reconstructed signal','Original signal');
% 
% figure; 
% subplot(421);plot(x); xlabel('Sampling point in the time domain');ylabel('amplitude');title('(a)');
% subplot(422);plot(x_noise); xlabel('Sampling point in the time domain');ylabel('amplitude');title('（b）');
% subplot(423);plot(f,real(tfreq_noise)); xlabel('Sampling point in the frequence domain');ylabel('amplitude');title('（c）');
% subplot(424);plot(hat_x_direct,'.-');hold on;plot(x,'r');legend('Reconstructed signal','Original signal'); xlabel('Signal in the time domain');ylabel('amplitude');title('（d）');
% subplot(425);plot(hat_x_FGG,'.-');hold on;plot(x,'r');legend('Reconstructed signal','Original signal'); xlabel('Signal in the time domain');ylabel('amplitude');title('（e）');
% subplot(426);plot(hat_x_lowrank,'.-');hold on;plot(x,'r');legend('Reconstructed signal','Original signal'); xlabel('Signal in the time domain');ylabel('amplitude');title('（f）');
% subplot(427);plot(hat_x_Min_Max,'.-');hold on;plot(x,'r');legend('Reconstructed signal','Original signal'); xlabel('Signal in the time domain');ylabel('amplitude');title('（g）');
% subplot(428);plot(hat_x_nusft,'.-');hold on;plot(x,'r');legend('Reconstructed signal','Original signal'); xlabel('Signal in the time domain');ylabel('amplitude');title('（h）');
% set(gcf,'PaperType','a3');


%% %%%
% norm(hat_x.'-x)/norm(x)                           %  重构误差 
% error=abs(abs(hat_x)-abs(x));
% figure;
% plot(error) 
%% 重构信号
% clc
% clear
% close all
%  
% % Analog signal
% Dt = 0.00005;
% t = - 0.005:Dt:0.005;
% xa = exp(-1000 * abs(t));
%  
% subplot(3,1,1);
% plot(1000*t,xa);
% title('Analog signal');
% xlabel('t in msec');
% ylabel('xa');
%  
% % Discrete-time signal
% Ts = 0.0002;
% Fs = 1/Ts;
% n = -25:25;
% nTs = n*Ts;
% x = exp(-1000*abs(nTs));
%  
% subplot(3,1,2)
% % plot(1000*t,xa);
% % hold on
% stem(n*Ts*1000,x);
% title('Discrete-time signal');
% % hold off
%  
% % Analog signal reconstruction
% xa_r = x * sinc( Fs * ( ones(length(n),1) * t - nTs' * ones(1,length(t) ) ));
% subplot(3,1,3);
% plot(1000*t,xa_r);
% title('Analog signal reconstruction');
% xlabel('t in msec');
% ylabel('xa after reconstruction');
% hold on 
% stem(n*Ts*1000,x)
% hold off
% %check
% error = max(abs(xa_r - xa ))