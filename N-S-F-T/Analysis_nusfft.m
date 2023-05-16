% clc
clear
k=35;   % 信号稀疏度
n=4096;% 信号长度
SNR_dB = 35; % signal SNR in dB

% N=[256,512,1024,2048,4096,8192];
% K=5:5:500; % 信号稀疏度
% SNR_dB_all=-10:0.5:20;
% for qw=1:length(SNR_dB_all)
%     NUsfft_allerr=zeros(length(SNR_dB_all));
%     n=N(qw);      % 信号长度
     disp('********************************');
     disp(['Signal size is ', num2str(n)]);
%     k=K(qw);
     disp(['Signal Sparsity is ', num2str(k)]);
%     SNR_dB=SNR_dB_all(qw);
     disp(['Signal noise is ', num2str(SNR_dB)]);
%% 构造时域上的随机信号
% SNR = 20^(SNR_dB/20);
[Sig,Sig_freq,f, A ,u] = genSigOffGrid(n,k,SNR_dB); % 构造随机的位置和幅值
%% NUDFT
% True_freq = ndftld(Sig, n, f);
%% NUSFFT
[NUSFFT_freq,BB_loc,BB_est,b_loc,b_est]=nusfft(Sig,k,f);
%% 结果图
figure;
% subplot(311);plot(f,abs(True_freq));xlabel('NUDFT frequency');ylabel('amplitude');
subplot(312);plot(f,abs(Sig_freq));xlabel('True frequency');ylabel('amplitude');
subplot(313);plot(f,abs(NUSFFT_freq));xlabel('NUSFFT Estimate frequency');ylabel('amplitude');
% NUsfft_err=0;
% for i=1:n
%     NUsfft_err = NUsfft_err+abs(abs(Sig_freq(i))-abs(NUSFFT_freq(i)));
% end
NUsfft_err=norm(Sig_freq-NUSFFT_freq);
NUsfft_err=NUsfft_err/n;
disp(['NUSFFT的 Error is ',  num2str(NUsfft_err)]);

% end

