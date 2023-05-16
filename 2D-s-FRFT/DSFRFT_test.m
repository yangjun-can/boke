clear;
%clc;
% close all;
%% 构造随机稀疏信号
SNR_dB = 80;  %signal SNR in dB 分数域、傅里叶域信噪比
k = 600;  % signal sparsity
N1=256;
N2=256;
N=N1*N2;
p1=0.8;p2=0.8;ts1=0.005;ts2=0.005;
sigma_n = 1; % 噪声的时域标准差 set this as 0 for noiseless signals, otherwise set as 1
if sigma_n == 1
    SNR = 10^(SNR_dB/10);
    sigma_f = sigma_n*sqrt(N);  % 傅里叶域噪声的标准差
    sigma_frft = sigma_n;  % 分数域噪声的标准差
    a_min = sqrt(sigma_frft^2*SNR); % 最小幅度
else
    a_min=10^(SNR_dB/10);
end
[Sig, A, u, v] = genSigOnGrid(N1,N2,k,a_min,sigma_n);
% figure;mesh(0:N1-1,0:N2-1,abs(Sig));zlabel('Amplitude');%title('Sparse signal in the time domain');
% figure;colormap(gray);imagesc(0:N1-1,0:N2-1,abs(Sig));
%% 真实频率
Omega_gt = [u,v];
X_true=zeros(N1,N2);
X_true(sub2ind([N1,N2],u+1,v+1)) = A; % 真实频率
%figure;mesh(0:N1-1,0:N2-1,abs(X_true));
% figure;colormap("colorcube");imagesc(0:N1-1,0:N2-1,abs(X_true));  %二维图，颜色深浅表述幅度大小

%% 2D DFRFT 分解为一维 所得频率
% tic
% X_frft = DFRFT_2D_fft(Sig,p1,p2,ts1,ts2);
% disp([' 按fft对行、列分别DFRFT的时间为 ',  num2str(toc)]);
% disp([' 按fft对行、列分别DFRFT的误差为 ',  num2str(norm(X_true-X_frft))]);
% figure;colormap("hot");mesh(0:N1-1,0:N2-1,abs(X_frft));zlabel('Amplitude');axis([0 N1-1 0 N2-1]);%title("Decomposition method");
% figure;colormap("colorcube");imagesc(0:N1-1,0:N2-1,abs(X_frft));  %二维图，颜色深浅表述幅度大小
% figure;surf(0:N1-1,0:N2-1,abs(X_frft));zlabel('Amplitude');

%% 2D SFRFT 所估计频率
tic
[X_sfrft,us1,us2] = DFRFT_2D_sft(Sig,p1,p2,N1,N2,ts1,ts2,sigma_n,Omega_gt);
disp(['sDFRFT的时间为 ',  num2str(toc)]);
disp(['sDFRFT的误差为 ',  num2str(norm(X_true-X_sfrft))]);
%  figure;colormap("hot");mesh(0:N1-1,0:N2-1,abs(X_sfrft));zlabel('Amplitude');axis([0 N1-1 0 N2-1]);%title("Our methed");
%  figure;colormap("colorcube");imagesc(0:N1-1,0:N2-1,abs(X_sfrft));  %二维图，颜色深浅表述幅度大小
 %% 直接法
%  tic
%  X_frft_direct = frft_direct(Sig,p1,p2,ts1,ts2);
%  disp([' 直接法的时间为 ',  num2str(toc)]);
%  disp([' 直接法的误差为 ',  num2str(norm(X_true-X_frft_direct))]);
% figure;colormap("hot");mesh(0:N1-1,0:N2-1,abs(X_frft_direct));zlabel('Amplitude');axis([0 N1-1 0 N2-1]);%title("Our methed");
% figure;colormap("colorcube");imagesc(0:N1-1,0:N2-1,abs(X_frft_direct));  %二维图，颜色深浅表述幅度大小
%% 2D DFRFT 分解为fft2 所得频率
% tic
% [X_peifrft,u1,u2] = DFRFT_2D_fft2(Sig,p1,p2,N1,N2,ts1,ts2);
% disp([' 按fft2求DFRFT的时间为 ',  num2str(toc)]);
% disp([' 按fft2求DFRFT的误差为 ',  num2str(norm(X_true-X_peifrft))]);
%  figure;mesh(0:N1-1,0:N2-1,abs(X_peifrft));title("DFRFT_pei");