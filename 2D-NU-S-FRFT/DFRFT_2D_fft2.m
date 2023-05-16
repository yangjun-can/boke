function [X_frft,us1,us2] = DFRFT_2D_fft2(x,p1,p2,N1,N2,ts1,ts2)
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

%% 第一次调制
s=ts1*(0:N1-1);                           % 时域采样点
t=ts2*(0:N2-1);
% chirp_s = diag(exp(1i/2*cot(alpha)*s.^2));
% chirp_t = diag(exp(1i/2*cot(beta)*t.^2));
% x_t=chirp_s*x*chirp_t;                    % 第一次调制
chirp_s = exp(1i/2*cot(alpha)*s.^2);
chirp_t = exp(1i/2*cot(beta)*t.^2);
x_t=x.*repmat(chirp_s.',1,N2).*repmat(chirp_t,N1,1);
%% 2D DFT
if sin(alpha)>0 && sin(beta)>0
    X_DFT=fft2(x_t);
elseif sin(alpha)<0 && sin(beta)<0    
    X_DFT=ifft2(x_t)*N1*N2;                                      % sft逆变换
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
% chirp_u = exp(1i*pi*cot(alpha)*u.^2);
k1=sqrt(S1*(sin(alpha)-1i*cos(alpha))/M1);
k2=sqrt(S2*(sin(beta)-1i*cos(beta))/M2);  % 第二次调制系数
% chirp_u = diag(exp(1i/2*cot(alpha)*u.^2));
% chirp_v = diag(exp(1i/2*cot(beta)*v.^2));                       
% X_frft=k1*k2*chirp_u*X_DFT*chirp_v;
chirp_u = exp(1i/2*cot(alpha)*u.^2);
chirp_v = exp(1i/2*cot(beta)*v.^2);
X_frft=k1*k2*X_DFT.*repmat(chirp_u.',1,M2).*repmat(chirp_v,M1,1);
end