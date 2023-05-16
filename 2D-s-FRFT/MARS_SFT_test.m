%% Test of Robust MARS-SFT
% clear;
close all;
N0=2048;
N1=2048;
N = N0*N1;
L = lcm(N0,N1); % line length

K = 10;  % signal sparsity
gama = 74;  % 定位矫正参数
T = 1;   % 迭代次数

% n_d = 2;  % 内循环次数
% n_s = 3;  % 被投票次数
%% Generate signal
sigma_n = 1; % 噪声的时域标准差 set this as 0 for noiseless signals, otherwise set as 1
sigma_f = sqrt(N)*sigma_n;  % 频域噪声的标准差
if sigma_n==1
%  A = sqrt(sigma_n^2*SNR);     % chirp信号时域的单频峰值
  SNR_dB = 60;  % 最小频率 信噪比
  SNR = 10^(SNR_dB/10);
  a_min = sqrt(sigma_f^2*SNR); % 最小幅度
%     Pd=0.99;  % 大桶检测概率
%     Pfa=Pd^(L*SNR+1); % 大桶虚警概率
     Pfa=1e-2;
     epsilon = -2*L*log(Pfa)*sigma_n^2; %Threshold for frequency detectoin in a slice; set as 1e-10 for noiseless case
     gamma = 1*1e4;   % Threshold for 1-sparse detection; set as 1e-10 for noiseless case
else
    a_min = 10^(SNR_dB/10);
    epsilon = 1e-10;% %Threshold for frequency detectoin in a slice; set as 1e-10 for noiseless case
    gamma = 1e-10; % Threshold for 1-sparse detection; set as 1e-10 for noiseless case
end

% Generate on-grid frequencies
 [Sig, A, u, v] = genSigOnGrid(N0,N1,K,a_min,sigma_n);

% Generate off-grid frequencies
%[Sig, A, u, v] = genSigOffGrid(N0,N1,K,a_min,sigma_n);

% Generate on-grid frequency clusters
%[Sig, A, u, v] = genSigOnGridClusters(N0,N1,K,2,a_min,sigma_n); K = K*25;

%% the groud truth
Omega_gt = [u,v];  % ground truth frequency locations
X_true=zeros(N0,N1);
X_true(sub2ind([N0,N1],u+1,v+1)) = A; % 真实频率
% figure;mesh(0:N0-1,0:N1-1,abs(X_true));
%%  fft
% tic
% X2=fft2(Sig);
% disp([' fft2的时间为 ',  num2str(toc)]);
% figure;mesh(0:N0-1,0:N1-1,abs(X2));title("FFT2结果")
% error2=norm(X2-X_true)
%% Generate window; genearate rect window for on-grid cases
%% Chebyshev window
% att = 70; % PSR in dB
% rho_w = 10^(att/10);
% win0 = chebwin(N0,att);
% win1 = chebwin(N1,att);
%  chebwin(n,r)产生n点的窗函数，其傅里叶变换后的旁瓣波纹比主瓣低rdB，旁瓣是等波纹的。注意：当n为偶数时，窗函数的长度为(n+1)。

%% Rectangular window
win0 = ones(N0,1);
win1 = ones(N1,1);
%% MARS_SFT
Win = win0*win1.';
%A的共轭转置矩阵为A'，转置为A.'
% [Omega, hA, P] = MARS_ISFT(Sig, Win, N0, N1, T, epsilon, gamma, n_d, n_s);
tic
[Omega, hA] = MARS_SFT_corr(Sig, Win, N0, N1, T, epsilon, gamma, gama);
% [Omega, hA] = MARS_SFT_vote(Sig, Win, N0, N1, T, epsilon, gamma, n_d, n_s);
% [Omega, hA] = MARS_SFT_nonoise(Sig, Win, N0, N1, T);
disp([' MARS-SFT的时间为 ',  num2str(toc)]);

%% reconstructed results
[Omega, ind] = sortrows(Omega); % 估计的位置（升序）
%[B,index] = sortrows(A) 默认依据A第一列的数值按升序移动每一行，如果第一列的数值有相同的，依次往右比较。
hA = hA(ind);   % 估计的频率值
X=zeros(N0,N1);
u = Omega(:,1);%（行位置）
v = Omega(:,2);%（列位置）
X(sub2ind([N0,N1],u+1,v+1)) = hA;
% figure;mesh(0:N0-1,0:N1-1,abs(X));

%% Visualize
Omega_gt = round(Omega_gt); % 四舍五入
visual_localization(N0,N1,Omega_gt,Omega);
% figure;spy(X);%稀疏矩阵可视化
error=norm(abs(X-X_true))/N
