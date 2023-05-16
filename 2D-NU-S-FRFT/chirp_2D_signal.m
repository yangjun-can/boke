%clc
%clear
%close
%% 2D chirp 信号
K = 10; %稀疏度
a_max = 18.3;                     
sigma_n = 1; % 噪声方差
delta = 0.1; 
% SNR_dB = 80; % 信噪比
% SNR = 10^(SNR_dB/20);
% a_max = sqrt(sigma_n^2*SNR); % 幅度
fs=256;  %采样率Hz
T=1;   %脉宽：脉冲持续时间
ts=1/fs;  % 采样间隔 
n=round(T*fs); %采样点个数，round是取整
% x=linspace(0,T-1,n); %时间向量
% y=linspace(0,T-1,n); %时间向量
x=(0:n-1)*ts; X=ones(n,1)*x;% 时间向量
y=(0:n-1)*ts; Y=y'*ones(1,n);% 时间向量
A1=a_max;   A2=a_max;  A3=a_max;   % 幅度
M1= 197.9; M2= 28.7;  M3=218.8;   % 行初始频率 【0，fs】之间，且和T的乘积是整数时分数域能量最集中（完全稀疏），否则有泄露 
N1=0.8;  N2=N1;   N3=N1; % 行调频率
R1= 224.2;  R2= 184.1; R3= 54.2;  % 列初始频率
S1=0.24;  S2=S1;   S3=S1;%-7.8 % 列调频率
F1=A1*exp(1j*pi*(2*M1*X+N1*X.^2+2*R1*Y+S1*Y.^2)); %信号1
F2=A2*exp(1j*pi*(2*M2*X+N2*X.^2+2*R2*Y+S2*Y.^2)); %信号2
F3=A3*exp(1j*pi*(2*M3*X+N3*X.^2+2*R3*Y+S3*Y.^2)); %信号2
F=F1+F2+F3;
% F = zeros(n,n);
% A = zeros(K,1);
% M = zeros(K,1);
% N = zeros(K,1);
% R = zeros(K,1);
% S = zeros(K,1);
% for ii = 1:K
%     A(ii) = a_max;% 幅度
%     M(ii) = (randi([0 n-1])+2*delta*rand(1)-delta)*fs/n;% 行初始频率
%     N(ii) = 0.8;% 行调频率
%     R(ii) = (randi([0 n-1])+rand(1)-0.5)*fs/n;% 列初始频率
%     S(ii) = 0.24;% 列调频率
%     F = F + A(ii)*exp(1j*pi*(2*M(ii)*X+N(ii)*X.^2+2*R(ii)*Y+S(ii)*Y.^2)); %信号1
% end
noise = sigma_n*sqrt(2)/2*(randn(n,n)+1i*randn(n,n));  % 时域加复高斯噪声
SNR_t=snr(F,noise)
F = F+noise;%含噪时域信号
% F = awgn(F,15.587,'measured');
%   figure;imagesc(x,y,real(F));colormap(gray(256));colorbar% 添加渐变色条title('2D chirp信号');
%  figure; colormap("hsv");mesh((0:n-1)*T/n,(0:n-1)*T/n,real(F));zlabel('Amplitude');%axis([0 T 0 T]);title("Our methed");
%  colormap("winter");

%% 2D DFRFT
% 最佳旋转角
% alfa=mod(acot(-N(1)),pi); p1=2*alfa/pi;
% beta=mod(acot(-S(1)),pi); p2=2*beta/pi;
alfa=mod(acot(-N1),pi); p1=2*alfa/pi;
beta=mod(acot(-S1),pi); p2=2*beta/pi;

% 2个1D DFRFT
% tic
% X_frft = DFRFT_2D_fft(F,p1,p2,ts,ts);
% disp([' 按fft对行、列分别DFRFT的时间为 ',  num2str(toc)]);
% figure;colormap("hot");mesh((0:n-1)*fs/n,(0:n-1)*fs/n,abs(X_frft)/n);zlabel('Amplitude');axis([0 fs 0 fs]);%title("Decomposition method");
% figure;colormap("colorcube");imagesc((0:n-1)*fs/n,(0:n-1)*fs/n,abs(X_frft)/n);  %二维图，颜色深浅表述幅度大小
  % 2D SFRFT 所估计频率
tic
[X_sfrft,Omega,us1,us2] = DFRFT_2D_sft(F,p1,p2,n,n,ts,ts,sigma_n);
disp(['sDFRFT的时间为 ',  num2str(toc)]);
%  figure;colormap("summer");mesh((0:n-1)*fs/n,(0:n-1)*fs/n,abs(X_sfrft)/n);zlabel('Amplitude');axis([0 fs 0 fs]);%title("Our methed");
 % figure;colormap("summer");imagesc((0:n-1)*fs/n,(0:n-1)*fs/n,abs(X_sfrft)/n);  %二维图，颜色深浅表述幅度大小
% colormap("hot");colormap("summer");colormap("colorcube");
%% 位置可视化 （与真实位置对比）
Omega_gt = [R1 M1
    R2 M2
    R3 M3]; 
Omega_gt = round(Omega_gt); % 四舍五入
% visual_localization(n,n,Omega_gt,Omega*fs/n);
visual_localization(fs,fs,Omega_gt*n/fs,Omega);
