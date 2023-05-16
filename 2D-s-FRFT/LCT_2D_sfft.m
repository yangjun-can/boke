function [X_lct,us1,us2] = LCT_2D_sfft(x, a,b,c,d,alfa,beta,gama,delta, N1,N2,ts1,ts2,k)
%%
% Input :
%  x 信号
%  abcd ,alfa,beta,gama,delta 变换参数
%  N1*N2 信号大小
%  ts1,ts2  时域采样间隔
%  k 稀疏度
% Output：
%  X_lct：线性正则域频率
%  us1，us2 线性正则域采样间隔
M1=N1;                           % 线性正则域与时域离散点个数相同
M2=N2;

S1=sign(b);
us1=2*pi*S1*b/(ts1*M1);% 分数域采样间隔
S2=sign(beta);
us2=2*pi*S2*beta/(ts2*M2);
%% 特殊情况
% A=[a,b;c,d];
% r=det(A);
% if r==1   
%    if A==[0,1;-1,0]
%        X_lct=fft(x);
%        return;
%    end
% 
%    if A==[0,-1;1,0]
%        X_temp=fft(x);
%        k = L:-1:1;
%        X_lct = X_temp(k);
%        return;
%    end

%    if A==[cos(alpha),sin(alpha);-sin(alpha),cos(alpha)]
%        X_temp=dfrft(x);
%        k = L:-1:1;
%        X_lct = X_temp(k);
%        return;
%    end
%% 第一次调制
s=ts1*(0:N1-1);
t=ts2*(0:N2-1);
chirp_s = diag(exp(1i*a/(2*b)*s.^2));
chirp_t = diag(exp(1i*alfa/(2*beta)*t.^2));
x_t=chirp_s*x*chirp_t;
%% 2D DFT
if b>0
    L = lcm(N1,N2); % line length
    epsilon = 1e-10; %Threshold for frequency detectoin in a slice; set as 1e-10 for noiseless case % 1e-10表示1*10^-10，科学计数法；
    gamma = 1e-10; % Threshold for 1-sparse detection; set as 1e-10 for noiseless case
    n_s = 1;  % 内循环次数
    n_d = 1;  % 被投票次数
    T = 10;   % 迭代次数
    Win=ones(N1,N2);
    [Omega,A,P] = MARS_SFT(x_t,Win,N1,N2,T,epsilon,gamma,n_d,n_s);  % 2D sft
    X_lct=zeros(N1,N2);
    X_lct(Omega)=A;
else
    X_lct=ifft2(x_t)*N1*N2;
end
%% 第二次调制
u=us1*(0:M1-1);
v=us2*(0:M2-1);
chirp_u = diag(exp(1i*d/(2*b)*u.^2));
chirp_v = diag(exp(1i*delta/(2*beta)*v.^2));
X_lct=sqrt(1/(M1*M2))*chirp_u*X_lct*chirp_v;

%   if b==0
%         k = L:-1:1; 
%         remainder=rem(d,1);
%         if remainder==0
%             X_lct = sqrt(d)*exp(1i/2*c*d*u.^2)*x(d*k);
%         else
%             y=sqrt(1/N)*sfft(x);
%             us=N*abs(a)*ts/M;
%             u=us*(-M/2+1:M/2);
%             z=sqrt(1/M)*y*exp(1i*c/(2*a)*u');
%             if a<0 
%                 X_lct =sfft(z);
%             else
%                 X_lct =isfft(z);
%             return;
%             end
%         return
%         end
%    end
end