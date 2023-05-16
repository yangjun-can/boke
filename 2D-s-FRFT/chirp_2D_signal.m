clc
clear
%% 构造时域上的1D线性调频信号
% n=256;
% fs=n;          % 采样频率
% ts=1/fs;       % 采样间隔 
% t1=(0:n-1)*ts; % 时间向量
% fd=[100 150 200 200];     % 初始频率
% mu=[-1.7 -1.4 -1.2 -2.5]; % 调频率
% x1=exp(1j*pi*mu(1)*(t1.^2)).*exp(1j*2*pi*fd(1)*t1);
% x2=exp(1j*pi*mu(2)*(t1.^2)).*exp(1j*2*pi*fd(2)*t1);
% x3=exp(1j*pi*mu(3)*(t1.^2)).*exp(1j*2*pi*fd(3)*t1);
% x4=exp(1j*pi*mu(4)*(t1.^2)).*exp(1j*2*pi*fd(4)*t1);
% x=x1+x2/2+x3/3+x4/4;      % 线性调频信号（时域）
%X=fftshift(fft(x)); X1=X/fs;   % 线性调频信号幅值（频域）

%% 2D chirp 信号
fs=2048;  %采样率Hz
T=2^-4;   %脉宽：脉冲持续时间
A1=0.5;  A2=0.5; A3=0.5;   % 幅度
M1= 1088; M2=240; M3=1824; %640 % 行初始频率 【0，fs】之间，且和T的乘积是整数时分数域能量最集中（完全稀疏），否则有泄露 
N1=-0.8; N2=N1;  N3=N1;%-16.4% 行调频率
R1= 416;  R2=1600;R3= 1632;%896 % 列初始频率
S1=0.24; S2=S1;  S3=S1;%-7.8 % 列调频率
ts=1/fs;       % 采样间隔 
n=round(T*fs); %采样点个数，round是取整
% x=linspace(0,T-1,n); %时间向量
% y=linspace(0,T-1,n); %时间向量
x=(0:n-1)*ts; % 时间向量
y=(0:n-1)*ts; % 时间向量
Lx=length(x);
Ly=length(y);
X=ones(Ly,1)*x;
Y=y'*ones(1,Lx);
F1=A1*exp(1j*pi*(2*M1*X+N1*X.^2+2*R1*Y+S1*Y.^2)); %信号1
F2=A2*exp(1j*pi*(2*M2*X+N2*X.^2+2*R2*Y+S2*Y.^2)); %信号2
F3=A3*exp(1j*pi*(2*M3*X+N3*X.^2+2*R3*Y+S3*Y.^2)); %信号2
F=F1+F2+F3;
% F = awgn(F,4.587,'measured');
%  figure;imagesc(x,y,real(F));title('2D chirp信号');colormap(gray(256));colorbar% 添加渐变色条
% figure; colormap("hsv");mesh((0:n-1)*T/n,(0:n-1)*T/n,real(F));zlabel('Amplitude');%axis([0 T 0 T]);title("Our methed");
% colormap("winter");

%% 2D DFRFT
% 最佳旋转角
alfa=mod(acot(-N1),pi); p1=2*alfa/pi;
beta=mod(acot(-S1),pi); p2=2*beta/pi;
  % 2个1D DFRFT
tic
X_frft = DFRFT_2D_fft(F,p1,p2,ts,ts);
disp([' 按fft对行、列分别DFRFT的时间为 ',  num2str(toc)]);
figure;colormap("hot");mesh((0:n-1)*fs/n,(0:n-1)*fs/n,abs(X_frft)/n);zlabel('Amplitude');axis([0 fs 0 fs]);%title("Decomposition method");
figure;colormap("colorcube");imagesc((0:n-1)*fs/n,(0:n-1)*fs/n,abs(X_frft)/n);  %二维图，颜色深浅表述幅度大小
  % 2D SFRFT 所估计频率
% tic
% [X_sfrft,us1,us2] = DFRFT_2D_sft(F,p1,p2,n,n,ts,ts,1);
% disp(['sDFRFT的时间为 ',  num2str(toc)]);
% figure;colormap("summer");mesh((0:n-1)*fs/n,(0:n-1)*fs/n,abs(X_sfrft)/n);zlabel('Amplitude');axis([0 fs 0 fs]);%title("Our methed");
% figure;colormap("colorcube");imagesc((0:n-1)*fs/n,(0:n-1)*fs/n,abs(X_sfrft)/n);  %二维图，颜色深浅表述幅度大小
% % colormap("hot");colormap("summer");
   % MARS_SFT 所估计频率
% tic
% T=3;
% win0 = ones(n,1);
% win1 = ones(n,1);
% Win = win0*win1.';
% epsilon = 2e-6; %Threshold for frequency detectoin in a slice; set as 1e-10 for noiseless case
% gamma = 2e-6;
% [Omega, hA] = MARS_SFT_vote(F, Win, n, n, T, epsilon,gamma, 3, 2);
% [Omega, ind] = sortrows(Omega); % 估计的位置（升序）
% hA = hA(ind)/n;   % 估计的频率值
% X_SFT=zeros(n,n);
% u = Omega(:,1);%（行位置）
% v = Omega(:,2);%（列位置）
% X_SFT(sub2ind([n,n],u+1,v+1)) = hA;
% disp(['MARS_SFT的时间为 ',  num2str(toc)]);
% figure;colormap("summer");mesh((0:n-1)*fs/n,(0:n-1)*fs/n,abs(X_SFT)/n);zlabel('Amplitude');axis([0 fs 0 fs]);%title("Our methed");

