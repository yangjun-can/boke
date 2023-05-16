clear
% clc
%% 构造随机信号
% k=1;% 信号稀疏度
% N=512;% 信号长度
% disp(['信号长度为：',num2str(N)]);
% % SNR = 0:5:35; %signal SNR in dB
% Min_Max_l2err=zeros(1,length(SNR) );
% Low_Rank_l2err=zeros(1,length(SNR) );
% FGG_l2err=zeros(1,length(SNR) );
% NUsfft_l2err=zeros(1,length(SNR) );

% % for index=1:length(SNR)    SNR_dB=SNR(index);
% SNR_dB=40;
% disp(['信号信噪比为：',num2str(SNR_dB)]);
% [Sig,Sig_freq,f, A ,u] = genSigOffGrid(N,k,SNR_dB); % 构造随机的位置和幅值
% figure;
% subplot(512); plot(f,abs(Sig_freq)); xlabel('True frequency');ylabel('amplitude');
%% 构造指定频率的信号
% N=512; % 信号长度
% fs=512;     % 采样频率
% ts=1/fs;      % 采样间隔 
% t=(0:N-1)*ts; % 时间向量
% % t=0:N-1; % 采样点
% k=4; % 信号稀疏度
% w1= 51.05*2*pi; w2=128.05*2*pi; % 指定频率
% Am1=200;% 频率幅度1
% Am2 = Am1+10^(32/10);% 频率幅度2   两频率相差的30dB
% x1=Am1*sin(w1*t);x2=Am2*sin(w2*t);
% Sig = x1 + x2; % 信号（时域）
% % SNR_dB = 10; %signal SNR in dB 时域信噪比
% % Sig = awgn(Sig,SNR_dB,'measured');%加信噪比为SNR_dB的高斯白噪声
% % Sig2 = Sig+sqrt(sigma_n^2/2)*(randn(1,N) +1i*randn(1,N));
% f=zeros(N,1); % 非均匀频率点
% for i=1:N
%     f(i)=Nonuniform_sampling_point(i-1,N);
% end
%% 构造Chirp信号
N=512;      % 信号长度
k=3;        % 信号稀疏度
fs=N;          % 采样频率
ts=1/fs;       % 采样间隔 
t1=(0:N-1)*ts; % 时间向量
fd=[100 150 200 200];     % 初始频率
mu=[-1.7 -1.4 -1.2 -2.5]; % 调频率
x1=exp(1j*pi*mu(1)*(t1.^2)).*exp(1j*2*pi*fd(1)*t1);
x2=exp(1j*pi*mu(2)*(t1.^2)).*exp(1j*2*pi*fd(2)*t1);
x3=exp(1j*pi*mu(3)*(t1.^2)).*exp(1j*2*pi*fd(3)*t1);
x4=exp(1j*pi*mu(4)*(t1.^2)).*exp(1j*2*pi*fd(4)*t1);
Sig=x1+x2/2+x3/3+x4/4;      % 线性调频信号（时域）

% SNR_dB=5;Sig = awgn(Sig,SNR_dB,'measured');%加信噪比为SNR_dB的高斯白噪声
f=zeros(N,1); % 非均匀频率点
for i=1:N
    f(i)=Nonuniform_sampling_point(i-1,N);
end
%% NUDFT
% tic; 
% True_freq = ndftld(Sig, N, f);
% disp(['Direct的运行时间为：',num2str(toc),'秒']);
% subplot(511);plot(f,abs(True_freq));xlabel('NUDFT frequency');ylabel('amplitude');
% NUDFT_l2err=norm(Sig_freq.'-True_freq)/N;disp(['NUDFT 的l2 Error is ',  num2str(NUDFT_l2err)]);

%% Min_Max NUFFT
% tic;
%  Min_Max_freq = Min_Max (Sig,N,f);
% temp=2*pi*f/N;
% J=6;
% K=2*N;
% st = nufft_init(temp, N, J, K, 'minmax:kb');
% Min_Max_freq = nufft(Sig.',st);
% disp(['Min_Max_NUFFT的运行时间为：',num2str(toc),'秒']);
% subplot(513);plot(f,abs(Min_Max_freq));xlabel('MinMaxNUFFT Estimate frequency');ylabel('amplitude');
%  Min_Max_l2err=norm(Sig_freq.'-Min_Max_freq)/N;
%  disp(['Min_Max_NUFFT l2 Error is ',  num2str(Min_Max_l2err)]);
% % Min_Max_l2err(index)=norm(Sig_freq.'-Min_Max_freq)/N;

%% Low_Rank_Approximation
% tic;
% Low_Rank_NUFFT = chebfun.nufft(Sig.',f/N,2);
% disp(['Low_Rank_NUFFT的运行时间为：',num2str(toc),'秒']);
% subplot(514);plot(f,abs(Low_Rank_NUFFT));xlabel('LowRankNUFFT Estimate frequency ');ylabel('amplitude');
% Low_Rank_l2err=norm(Sig_freq.'-Low_Rank_NUFFT)/N;
% disp(['Low_Rank_NUFFT l2 Error is ',  num2str(Low_Rank_l2err)]);
% %  Low_Rank_l2err(index)=norm(Sig_freq.'-Low_Rank_NUFFT)/N;

%%  FGG
% tic;
% FGG_NUFFT = FGG_1d_type2(Sig,f*2*pi/N);
% FGG_NUFFT = fftshift(FGG_NUFFT);
% disp(['FGG_NUFFT的运行时间为：',num2str(toc),'秒']);
% subplot(515);plot(f,abs(FGG_NUFFT));xlabel('FGG_NUFFT  Estimate frequency ');ylabel('amplitude');
% FGG_l2err=norm(Sig_freq.'-FGG_NUFFT)/N;
% disp(['FGG_NUFFT l2 Error is ',  num2str(FGG_l2err)]);
% %  FGG_l2err(index)=norm(Sig_freq.'-FGG_NUFFT)/N;

%% nusfft
[NUSFFT_freq,BB_loc,BB_est,b_loc,b_est]=nusfft(Sig,k,f);
% [NUSFFT_freq,BB_loc,b_loc]=nusfft_iter(Sig,k,f);
% subplot(515);plot(f,abs(NUSFFT_freq));xlabel('NUSFFT Estimate frequency');ylabel('amplitude');

%  figure;
% subplot(211);plot(f,abs(True_freq));title('(a) Direct method');xlabel('Non-uniform sampling point in the frequency');ylabel('amplitude');
% subplot(211);plot(abs(Sig));title('(a) Noisy signal');xlabel('sampling point in the time');ylabel('amplitude');
% subplot(212);plot(f,abs(NUSFFT_freq));title('(b) NUSFT method');xlabel('Non-uniform estimation point in the frequency');ylabel('amplitude');
% NUsfft_l2err=norm(Sig_freq-NUSFFT_freq)/N;
% disp(['NUSFFT的 l2 Error is ',  num2str(NUsfft_l2err)]);
% %  NUsfft_l2err(index)=norm(Sig_freq-NUSFFT_freq)/N;

% % end
% % figure;
% % plot(SNR,Min_Max_l2err+0.1e-3,'--b');hold on;
% % plot(SNR,Low_Rank_l2err,'-m');hold on;
% % plot(SNR,FGG_l2err,'-k');hold on;
% % plot(SNR,NUsfft_l2err,':r','LineWidth',2);hold off
% % xlabel('SNR(dB)');ylabel('Average L2 error per entry');
% % legend('Min-Max interpolation','Low rank approximation','Fast gaussian gridding','NUSFT');







 %%
% function [Xk] = dft1d(xn, N)
% % Computes Uniform Discrete Fourier Trans form by multiplication of DFT matrix.
% %  ------------------------------------------------------------
% % xn= N-point 1D input sequence over 0 <=n <= N-1 (row vector)
% % Xk= 1D DFT coefficient array over  0<=k<= N-1
% % N= Length of DFT
% % Usage: Xk= dft1d(xn, N)
% n=0:1 : N-1;  % Index for input data
% k=0:1 : N-1;  % Index for DFT coefficients
% Wn= exp(-1i*2*pi/N) ; % Twiddle factor
% nk=k'*n;% Creates an N X N matrix
% DFTmtx= Wn .^nk; % DFT matrix (NX N)
% Xk= (DFTmtx* xn.'); % DFT coefficients (column vector)
% % xn.' is the non-conjugate transpose of xn and equals to "transpose (xn)"。
% end
% 
% function [xn] = idft1a(Xk, N)
% % Computes Uniform Inverse Discrete Fourier Transform.
% % -------------------------------------------------------
% % xn= N-point 1D input sequence over 0<=n<= N-1 (column vector)
% % Xk= 1D DFT coefficient array over 0 <= k <= N-1
% % N= Length of DFT
% % Usage: xn = idft1d(Xk, N)
% n=0 :1: N-1;% Index for input data
% k=0 :1: N-1; % Index for DFT coefficients
% Wn= exp(1i*2*pi/N); % Twiddle factor
% nk=n'*k; % Creates an NXN matrix
% IDFTmtx= Wn .^ nk; % IDFT matrix (NXN), equivalent to (DFTmtx)
% xn=(IDFTmtx * Xk).'/N; % Reconstructed sequence (row vector)
% end