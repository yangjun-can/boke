function [X_frft,Omega,us1,us2] = DFRFT_2D_sft(x,p1,p2,N1,N2,ts1,ts2,sigma_n)
%%
% Input :
%  x 信号
%  p 分数阶数
%  N1*N2 信号大小
%  ts1，ts2  时域采样间隔
%  k 稀疏度
% Output：
%  X_frft：分数阶频率
%  us1,us2 分数域采样间隔
M1=N1;                           % 分数域与时域离散点个数相同
M2=N2;
p1=mod(p1,4);                    % 变换阶数（变换阶以4为周期）
p2=mod(p2,4);
alpha=p1*pi/2;                   % 旋转角度（旋转角以2pi为周期）
beta=p2*pi/2;
S1=sign(sin(alpha));
us1=2*pi*S1*sin(alpha)/(ts1*M1); % 分数域采样间隔
S2=sign(sin(beta));
us2=2*pi*S2*sin(beta)/(ts2*M2);
% sft参数
L = round(lcm(N2,N1)); % 切片长度
T = 1;   % 迭代次数
% if sigma_n==1
%     Pd=0.9;  % 大桶检测概率
%     Pfa=Pd^(L*SNR+1); 
%     Pfa=1e-300;% 大桶虚警概率
%     epsilon = -2*L*log(Pfa)*sigma_n^2 
    epsilon = 1e6;  %Threshold for frequency detectoin in a slice; set as 1e-10 for noiseless case
    gamma = 2e4;   % Threshold for 1-sparse detection; set as 1e-10 for noiseless case
    gama= 50;   % 定位矫正参数
%     n_d = 3;  % 内循环次数
%     n_s = 2;  % 被投票次数
% end
%% Generate window; genearate rect window for on-grid cases
att = 30; % PSR in dB
win1 = chebwin(N1,att);%一维窗
%% 无窗 == 矩形窗
win0 = ones(N1,1);
Win = win0*win0.'; 
%% 乘积法
Win_mulp = win1*win1.';
%% 时域旋转法
Win_rot = getWin2D(win1);
%% 频域旋转
w1_freq = fft(win1);
w2_freq2 = fftshift(getWin2D(fftshift(w1_freq)));
w22 = ifft2(w2_freq2);
Win_rot_fre = mapminmax(w22,0,1);
%% 第一次调制
s=ts1*(0:N1-1);                           % 时域采样点
t=ts2*(0:N2-1);
chirp_s = exp(1i/2*cot(alpha)*s.^2);
chirp_t = exp(1i/2*cot(beta)*t.^2);
x_t=x.*repmat(chirp_s.',1,N2).*repmat(chirp_t,N1,1);
%% 2D DFT
if sin(alpha)>0 && sin(beta)>0
    [Omega, A] = MARS_SFT_corr(x_t, Win_rot_fre, N1, N2, T, epsilon, gamma, gama);  % 2D sft
%     [Omega,A] =MARS_SFT_vote(x_t, Win, N1, N2, T, epsilon, gamma, n_d, n_s);   
    [Omega, ind] = sortrows(Omega); % 估计的位置（升序） 
    A = A(ind);   % 估计的频率值
    X_DFT=zeros(N1,N2);
    u = Omega(:,1);%（行位置）
    v = Omega(:,2);%（列位置）
    X_DFT(sub2ind([N1,N2],u+1,v+1)) = A;
elseif sin(alpha)<0 && sin(beta)<0
    [Omega,A] = MARS_ISFT(x_t,Win,N1,N2,T,epsilon,gamma,n_d,n_s);  % 2D sft
    [Omega, ind] = sortrows(Omega); % 估计的位置（升序） 
    A = A(ind);   % 估计的频率值
    X_DFT=zeros(N1,N2);
    u = Omega(:,1);%（行位置）
    v = Omega(:,2);%（列位置）
    X_DFT(sub2ind([N1,N2],u+1,v+1)) = A;
    X_DFT=X_DFT*N1*N2;                                      % sft逆变换
elseif sin(alpha)>0 && sin(beta)<0
    X_temp=ifft(x_t)*N1;% 各列IDFT
    X_DFT=fft(X_temp,N2,2);
elseif sin(alpha)<0 && sin(beta)>0
    X_temp=ifft(x_t,N2,2)*N2;% 各行IDFT
    X_DFT=fft(X_temp);
end
%% 第二次调制
u=us1*(0:M1-1);
v=us2*(0:M2-1);
k1=sqrt(S1*(sin(alpha)-1i*cos(alpha))/M1);
k2=sqrt(S2*(sin(beta)-1i*cos(beta))/M2);                         % 第二次调制系数
chirp_u = exp(1i/2*cot(alpha)*u.^2);
chirp_v = exp(1i/2*cot(beta)*v.^2);
X_frft=k1*k2*X_DFT.*repmat(chirp_u.',1,M2).*repmat(chirp_v,M1,1);
end